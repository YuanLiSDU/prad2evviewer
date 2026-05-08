// replay_add_hv.cpp — append decoded HV trees to a replay ROOT file in place.
//
// For each input file:
//   * Read `scalers` to extract one (event_number, ti_ticks, unix_time)
//     pin per SYNC.  The min/max unix_time defines the run window.
//   * Decode the daily VMDF v2 archive (default /data/prad2/hv_data) over
//     that window via hv::HVDecoder, projecting to optional --channels.
//   * For every HV / booster snapshot, interpolate the bracketing pins
//     to get its (event_number_at_arrival, ti_ticks_at_arrival) — same
//     contract the existing `epics` tree uses.
//   * Open the replay file in UPDATE mode and write five new trees:
//
//       hv                  — one entry per HV snapshot
//                              event_number_at_arrival/I  ti_ticks_at_arrival/L
//                              unix_time/i  (uint32 sec, matches scalers/epics)
//                              t_unix_s/D   (full fractional second)
//                              dv[N_ch]/F   v0set[N_ch]/F
//                              VMon = dv + v0set
//       hv_channels         — one entry per HV channel: channel_id/s, name/C
//       hv_booster          — booster snapshots (only if n_boosters > 0):
//                              event_number_at_arrival/I  ti_ticks_at_arrival/L
//                              unix_time/i  t_unix_s/D
//                              vmon/imon/vset/iset[N_bst]/F
//       hv_booster_channels — one entry per booster: booster_id/s, name/C
//       hv_meta             — one entry: run/HV-archive scalars
//
// Existing replay trees (recon/scalers/epics/...) are preserved.  If the
// file already carries hv/hv_channels/etc., the tool refuses unless -f is
// passed; -f deletes the old keys and writes fresh ones (the freed bytes
// are not reclaimed without a separate `rootcp -O` pass).
//
// CLI:
//   prad2ana_replay_add_hv [-d <hv_dir>] [-c "ch1,ch2,..."] [-p <pad_s>]
//                          [-f] file1.root [file2.root ...]
// Defaults:
//   -d /data/prad2/hv_data   -p 30   (no -c → keep all channels)

#include "HVDecoder.h"

#include <TFile.h>
#include <TKey.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr const char *kDefaultHvDir = "/data/prad2/hv_data";

const std::vector<const char *> kManagedTrees = {
    "hv", "hv_channels", "hv_booster", "hv_booster_channels", "hv_meta",
};


// ─────────────────────────────────────────────────────────────────────────
// CLI
// ─────────────────────────────────────────────────────────────────────────
void usage(const char *prog)
{
    std::fprintf(stderr,
"Usage: %s [options] file1.root [file2.root ...]\n"
"\n"
"  Append decoded HV trees (`hv`, `hv_channels`, `hv_booster`,\n"
"  `hv_booster_channels`, `hv_meta`) to a replay ROOT file in place.\n"
"\n"
"Options:\n"
"  -d <dir>     vmon_*.dat archive directory (default: %s)\n"
"  -c <list>    Comma-separated HV channel names to keep (default: all).\n"
"  -p <sec>     Pad around the scalers window in seconds (default: 30).\n"
"  -f           Overwrite existing HV trees in the file.\n"
"  -h           Show this help.\n",
        prog, kDefaultHvDir);
}

std::vector<std::string> split_comma(const std::string &s)
{
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') {
            if (!cur.empty()) { out.push_back(cur); cur.clear(); }
        } else if (!std::isspace(static_cast<unsigned char>(c))) {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

std::string iso_utc(double unix_s)
{
    std::time_t t = std::time_t(unix_s);
    std::tm utc{};
    gmtime_r(&t, &utc);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &utc);
    return buf;
}


// ─────────────────────────────────────────────────────────────────────────
// SYNC anchor + run window from the `scalers` tree.
// ─────────────────────────────────────────────────────────────────────────
//
// We need two things from scalers:
//
//   * The run window — min/max of `unix_time` (excluding the 0 sentinel
//     for pre-first-SYNC rows).  Used to pick which HV archive window
//     to decode.
//   * One anchor — any (ti_ticks, unix_time) pair from a real SYNC row
//     to convert HV-snapshot wall-clock times back onto the TI ticks
//     axis (1 ms = 250 000 four-ns ticks).
//
// The scalers tree records one row per (source, channel, slot) per
// SYNC, so within one SYNC all rows share (event_number, ti_ticks,
// unix_time).  We deduplicate to one entry per SYNC and pick the
// earliest as the anchor.  Multiple SYNC pins are exposed (just the
// (ti_ticks, unix_time_ms) pairs) so callers that want drift-aware
// fits can still get them.
//
struct SyncRow {
    Int_t    event_number;
    Long64_t ti_ticks;
    UInt_t   unix_time;          // seconds (matches the on-disk u32)
};

bool read_scaler_pins(const std::string &path,
                      std::vector<SyncRow> &pins,
                      std::string &err)
{
    pins.clear();
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        err = "cannot open " + path + " for read";
        return false;
    }
    auto *t = dynamic_cast<TTree *>(f->Get("scalers"));
    if (!t) {
        err = path + ": no `scalers` tree";
        return false;
    }

    Int_t    event_number = 0;
    Long64_t ti_ticks     = 0;
    UInt_t   unix_time    = 0;
    t->SetBranchStatus("*", 0);
    auto enable = [&](const char *name, void *addr) {
        if (!t->GetBranch(name)) return false;
        t->SetBranchStatus(name, 1);
        t->SetBranchAddress(name, addr);
        return true;
    };
    if (!enable("event_number", &event_number) ||
        !enable("ti_ticks",     &ti_ticks)     ||
        !enable("unix_time",    &unix_time))
    {
        err = path + ": scalers tree missing event_number / ti_ticks / unix_time";
        return false;
    }

    std::set<std::tuple<Int_t, Long64_t, UInt_t>> uniq;
    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (unix_time == 0) continue;
        uniq.emplace(event_number, ti_ticks, unix_time);
    }
    if (uniq.empty()) {
        err = path + ": no SYNC scaler rows with unix_time > 0";
        return false;
    }

    pins.reserve(uniq.size());
    for (const auto &tup : uniq) {
        pins.push_back({std::get<0>(tup), std::get<1>(tup), std::get<2>(tup)});
    }
    // Sort by (unix_time, ti_ticks) so multiple SYNCs in the same UTC
    // second are deterministically ordered earliest-first.  This matters
    // when picking the anchor pin: the run-start SYNC (smallest ti_ticks
    // in the smallest unix_time second) is what we want.
    std::sort(pins.begin(), pins.end(),
              [](const SyncRow &a, const SyncRow &b) {
                  if (a.unix_time != b.unix_time) return a.unix_time < b.unix_time;
                  return a.ti_ticks < b.ti_ticks;
              });
    return true;
}


// ─────────────────────────────────────────────────────────────────────────
// Read the events / recon tree's (event_num, timestamp) columns.
// ─────────────────────────────────────────────────────────────────────────
//
// `timestamp` is the per-event 48-bit TI tick (4 ns/tick) — the canonical
// "time within the run" that links physics events, scaler SYNC banks, and
// EPICS arrival rows.  `event_num` is the integer event id used as the
// join key elsewhere in the analysis (recon.event_num,
// scalers.event_number, epics.event_number_at_arrival).
//
// We try `recon` first (the recon-tree convention), then `events` (raw
// replay output).  Both carry the same branch names.
//
bool read_event_table(const std::string &path,
                      std::vector<int32_t> &event_num,
                      std::vector<int64_t> &ti_ticks,
                      std::string &tree_name,
                      std::string &err)
{
    event_num.clear();
    ti_ticks.clear();
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        err = "cannot open " + path + " for read";
        return false;
    }
    TTree *t = nullptr;
    for (const char *cand : {"recon", "events"}) {
        if ((t = dynamic_cast<TTree *>(f->Get(cand)))) {
            tree_name = cand;
            break;
        }
    }
    if (!t) {
        err = path + ": no `recon` or `events` tree";
        return false;
    }

    Int_t    en = 0;
    Long64_t ts = 0;
    t->SetBranchStatus("*", 0);
    auto enable = [&](const char *name, void *addr) {
        if (!t->GetBranch(name)) return false;
        t->SetBranchStatus(name, 1);
        t->SetBranchAddress(name, addr);
        return true;
    };
    if (!enable("event_num", &en) || !enable("timestamp", &ts)) {
        err = path + ": " + tree_name + " missing event_num / timestamp";
        return false;
    }

    Long64_t n = t->GetEntries();
    event_num.reserve(n);
    ti_ticks .reserve(n);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        event_num.push_back(en);
        ti_ticks .push_back(ts);
    }

    // The events should already be in TI-tick order (replay writes them
    // sequentially), but we don't depend on that — sort anyway so
    // upper_bound is well-defined.  Sort by ticks; reorder event_num the
    // same way.
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
                     [&](std::size_t a, std::size_t b) {
                         return ti_ticks[a] < ti_ticks[b];
                     });
    if (!std::is_sorted(idx.begin(), idx.end())) {
        std::vector<int32_t> en_s(n);
        std::vector<int64_t> ts_s(n);
        for (std::size_t k = 0; k < std::size_t(n); ++k) {
            en_s[k] = event_num[idx[k]];
            ts_s[k] = ti_ticks [idx[k]];
        }
        event_num = std::move(en_s);
        ti_ticks  = std::move(ts_s);
    }
    return true;
}


// ─────────────────────────────────────────────────────────────────────────
// V0Set / VSet / ISet matrix reconstruction.
// ─────────────────────────────────────────────────────────────────────────
//
// The segment carries CHTABLE / BOOSTER_TABLE event lists; for analysis
// it's nicer to have the active setpoints aligned with each snapshot row.
// "Most runs only have one CHTABLE event" so the resulting matrix is
// almost a tiled row, and ROOT's compression collapses it.
//
std::vector<float> build_v0set_matrix(const hv::HVSegment &seg)
{
    const int n_snap = seg.n_snapshots();
    const int n_ch   = seg.n_channels();
    std::vector<float> out(std::size_t(n_snap) * n_ch);
    for (int j = 0; j < n_ch; ++j) {
        auto col = seg.v0set_trace(j);
        for (int i = 0; i < n_snap; ++i)
            out[std::size_t(i) * n_ch + j] = col[i];
    }
    return out;
}

void build_booster_setpoint_matrices(const hv::HVSegment &seg,
                                      std::vector<float> &vset,
                                      std::vector<float> &iset)
{
    const int n_bsnap = seg.n_booster_snapshots();
    const int n_bst   = seg.n_boosters();
    vset.assign(std::size_t(n_bsnap) * n_bst,
                std::numeric_limits<float>::quiet_NaN());
    iset.assign(std::size_t(n_bsnap) * n_bst,
                std::numeric_limits<float>::quiet_NaN());
    if (n_bsnap == 0 || n_bst == 0 || seg.booster_events.empty())
        return;

    std::vector<int64_t> ev_ts;
    ev_ts.reserve(seg.booster_events.size());
    for (const auto &ev : seg.booster_events) ev_ts.push_back(ev.abs_ts_ms);

    for (int i = 0; i < n_bsnap; ++i) {
        auto it = std::upper_bound(ev_ts.begin(), ev_ts.end(),
                                   seg.booster_timestamps_ms[i]);
        int k = int(it - ev_ts.begin()) - 1;
        if (k < 0) k = 0;
        if (k >= int(seg.booster_events.size()))
            k = int(seg.booster_events.size()) - 1;
        for (int j = 0; j < n_bst; ++j) {
            vset[std::size_t(i) * n_bst + j] = seg.booster_events[k].vsets[j];
            iset[std::size_t(i) * n_bst + j] = seg.booster_events[k].isets[j];
        }
    }
}


// ─────────────────────────────────────────────────────────────────────────
// per-file driver
// ─────────────────────────────────────────────────────────────────────────
struct Summary {
    int        n_channels       = 0;
    int        n_boosters       = 0;
    int        n_hv_snapshots   = 0;
    int        n_booster_snaps  = 0;
    double     run_t_start_unix = 0.0;
    double     run_t_end_unix   = 0.0;
    int        n_sync_pins      = 0;
    int        n_events         = 0;
    std::string event_tree_name;
    int64_t    anchor_ti_ticks  = 0;
    int64_t    anchor_unix_ms   = 0;
    int32_t    ev_first         = -1;
    int32_t    ev_last          = -1;
    int        interval_ms      = 0;
    std::vector<std::string> source_archive;
    std::vector<std::string> channels;
    std::vector<std::string> booster_names;
    std::vector<std::string> overwrote;
};

bool process_file(const std::string &path,
                  const std::string &hv_dir,
                  const std::vector<std::string> &channels,
                  double pad_s,
                  bool   force,
                  Summary &out_summary,
                  std::string &err)
{
    // ── 1. scaler pins (run window + anchor) ───────────────────────────
    std::vector<SyncRow> pins;
    if (!read_scaler_pins(path, pins, err)) return false;

    double t_start = double(pins.front().unix_time);
    double t_end   = double(pins.back ().unix_time);
    double win_lo  = t_start - pad_s;
    double win_hi  = t_end   + pad_s;

    // Anchor for the unix↔ticks conversion: just take the earliest valid
    // SYNC pin.  TI ticks are exactly linear in time (4 ns/tick) so one
    // anchor is enough for run-scale windows.
    int64_t anchor_ti_ticks = pins.front().ti_ticks;
    int64_t anchor_unix_ms  = static_cast<int64_t>(pins.front().unix_time) * 1000;

    out_summary.run_t_start_unix = t_start;
    out_summary.run_t_end_unix   = t_end;
    out_summary.n_sync_pins      = int(pins.size());
    out_summary.anchor_ti_ticks  = anchor_ti_ticks;
    out_summary.anchor_unix_ms   = anchor_unix_ms;

    std::cout << "  run window  : " << iso_utc(t_start) << " → "
              << iso_utc(t_end) << "  (" << (t_end - t_start)
              << " s, " << pins.size() << " SYNC pins)\n";
    std::cout << "  +pad (" << pad_s << "s) : "
              << iso_utc(win_lo) << " → " << iso_utc(win_hi) << "\n";

    // ── 2. recon/events table (event_num, timestamp) ───────────────────
    std::vector<int32_t> evt_nums;
    std::vector<int64_t> evt_ticks;
    if (!read_event_table(path, evt_nums, evt_ticks,
                          out_summary.event_tree_name, err))
        return false;
    out_summary.n_events = int(evt_nums.size());
    std::cout << "  event tree  : " << out_summary.event_tree_name
              << " (" << evt_nums.size() << " entries)\n";
    std::cout << "  anchor      : ti_ticks=" << anchor_ti_ticks
              << "  unix_ms="   << anchor_unix_ms << "\n";

    // ── 3. existing-tree guard ──────────────────────────────────────────
    {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            err = "cannot open " + path + " for read";
            return false;
        }
        for (const char *n : kManagedTrees)
            if (f->Get(n)) out_summary.overwrote.push_back(n);
    }
    if (!out_summary.overwrote.empty() && !force) {
        std::ostringstream oss;
        oss << path << " already has HV tree(s):";
        for (const auto &s : out_summary.overwrote) oss << " " << s;
        oss << " — pass -f to overwrite";
        err = oss.str();
        return false;
    }

    // ── 4. decode HV ────────────────────────────────────────────────────
    hv::HVDecoder decoder(hv_dir);
    auto seg = decoder.load_window(win_lo, win_hi, channels, /*cache=*/"");
    if (seg.empty()) {
        err = "no HV data overlaps the window — is " + hv_dir
              + " the right archive?";
        return false;
    }

    // ── 5. associate snapshots with the *actual* recon event whose
    //       TI tick falls just before each snapshot's wall-clock time.
    auto hv_assoc = hv::HVSegment::associate_events(
        seg.timestamps_ms, evt_nums, evt_ticks,
        anchor_ti_ticks, anchor_unix_ms);
    auto bst_assoc = hv::HVSegment::associate_events(
        seg.booster_timestamps_ms, evt_nums, evt_ticks,
        anchor_ti_ticks, anchor_unix_ms);

    auto v0_matrix = build_v0set_matrix(seg);
    std::vector<float> bvset_matrix, biset_matrix;
    build_booster_setpoint_matrices(seg, bvset_matrix, biset_matrix);

    // ── 5. write trees in UPDATE mode ───────────────────────────────────
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "UPDATE"));
    if (!f || f->IsZombie()) {
        err = "cannot open " + path + " for UPDATE";
        return false;
    }

    if (force) {
        for (const char *n : kManagedTrees)
            f->Delete((std::string(n) + ";*").c_str());
    }

    const int n_ch    = seg.n_channels();
    const int n_bst   = seg.n_boosters();
    const int n_snap  = seg.n_snapshots();
    const int n_bsnap = seg.n_booster_snapshots();

    out_summary.n_channels      = n_ch;
    out_summary.n_boosters      = n_bst;
    out_summary.n_hv_snapshots  = n_snap;
    out_summary.n_booster_snaps = n_bsnap;
    out_summary.interval_ms     = seg.interval_ms;
    out_summary.source_archive  = seg.source_files;
    out_summary.channels        = seg.channels;
    out_summary.booster_names   = seg.booster_names;

    // ── hv tree ─────────────────────────────────────────────────────────
    {
        Int_t    event_number_at_arrival = 0;
        Long64_t ti_ticks_at_arrival     = 0;
        UInt_t   unix_time               = 0;
        Double_t t_unix_s                = 0.0;
        std::vector<float> dv_row   (std::max(1, n_ch));
        std::vector<float> v0_row   (std::max(1, n_ch));

        TTree *t_hv = new TTree("hv", "HV decoded snapshots (dV + V0Set per channel)");
        t_hv->Branch("event_number_at_arrival", &event_number_at_arrival,
                     "event_number_at_arrival/I");
        t_hv->Branch("ti_ticks_at_arrival", &ti_ticks_at_arrival,
                     "ti_ticks_at_arrival/L");
        t_hv->Branch("unix_time", &unix_time, "unix_time/i");
        t_hv->Branch("t_unix_s",  &t_unix_s,  "t_unix_s/D");
        if (n_ch > 0) {
            t_hv->Branch("dv",    dv_row.data(),
                         Form("dv[%d]/F",    n_ch));
            t_hv->Branch("v0set", v0_row.data(),
                         Form("v0set[%d]/F", n_ch));
        }

        for (int i = 0; i < n_snap; ++i) {
            event_number_at_arrival = hv_assoc.event_number[i];
            ti_ticks_at_arrival     = hv_assoc.ti_ticks[i];
            int64_t ts_ms = seg.timestamps_ms[i];
            unix_time = static_cast<UInt_t>(ts_ms / 1000);
            t_unix_s  = double(ts_ms) / 1000.0;
            for (int j = 0; j < n_ch; ++j) {
                dv_row[j] = seg.dv      [std::size_t(i) * n_ch + j];
                v0_row[j] = v0_matrix   [std::size_t(i) * n_ch + j];
            }
            t_hv->Fill();
        }
        if (!hv_assoc.event_number.empty()) {
            out_summary.ev_first = hv_assoc.event_number.front();
            out_summary.ev_last  = hv_assoc.event_number.back();
        }
        t_hv->Write("hv", TObject::kOverwrite);
    }

    // ── hv_channels tree ────────────────────────────────────────────────
    {
        UShort_t channel_id = 0;
        char     name_buf[64] = {0};
        TTree *t_ch = new TTree("hv_channels",
            "HV channel registry (channel_id ↔ name)");
        t_ch->Branch("channel_id", &channel_id, "channel_id/s");
        t_ch->Branch("name",       name_buf,    "name/C");
        for (int i = 0; i < n_ch; ++i) {
            channel_id = static_cast<UShort_t>(i);
            std::snprintf(name_buf, sizeof(name_buf), "%s",
                          seg.channels[i].c_str());
            t_ch->Fill();
        }
        t_ch->Write("hv_channels", TObject::kOverwrite);
    }

    // ── booster trees (only when present) ───────────────────────────────
    if (n_bst > 0 && n_bsnap > 0) {
        Int_t    event_number_at_arrival = 0;
        Long64_t ti_ticks_at_arrival     = 0;
        UInt_t   unix_time               = 0;
        Double_t t_unix_s                = 0.0;
        std::vector<float> vmon_row(n_bst);
        std::vector<float> imon_row(n_bst);
        std::vector<float> vset_row(n_bst);
        std::vector<float> iset_row(n_bst);

        TTree *t_b = new TTree("hv_booster",
            "Booster (TDK-Lambda) snapshots: VMon/IMon/VSet/ISet per supply");
        t_b->Branch("event_number_at_arrival", &event_number_at_arrival,
                    "event_number_at_arrival/I");
        t_b->Branch("ti_ticks_at_arrival", &ti_ticks_at_arrival,
                    "ti_ticks_at_arrival/L");
        t_b->Branch("unix_time", &unix_time, "unix_time/i");
        t_b->Branch("t_unix_s",  &t_unix_s,  "t_unix_s/D");
        t_b->Branch("vmon", vmon_row.data(), Form("vmon[%d]/F", n_bst));
        t_b->Branch("imon", imon_row.data(), Form("imon[%d]/F", n_bst));
        t_b->Branch("vset", vset_row.data(), Form("vset[%d]/F", n_bst));
        t_b->Branch("iset", iset_row.data(), Form("iset[%d]/F", n_bst));

        for (int i = 0; i < n_bsnap; ++i) {
            event_number_at_arrival = bst_assoc.event_number[i];
            ti_ticks_at_arrival     = bst_assoc.ti_ticks[i];
            int64_t ts_ms = seg.booster_timestamps_ms[i];
            unix_time = static_cast<UInt_t>(ts_ms / 1000);
            t_unix_s  = double(ts_ms) / 1000.0;
            for (int j = 0; j < n_bst; ++j) {
                vmon_row[j] = seg.booster_vmon[std::size_t(i) * n_bst + j];
                imon_row[j] = seg.booster_imon[std::size_t(i) * n_bst + j];
                vset_row[j] = bvset_matrix    [std::size_t(i) * n_bst + j];
                iset_row[j] = biset_matrix    [std::size_t(i) * n_bst + j];
            }
            t_b->Fill();
        }
        t_b->Write("hv_booster", TObject::kOverwrite);
    }

    if (n_bst > 0) {
        UShort_t booster_id = 0;
        char     name_buf[64] = {0};
        TTree *t_bch = new TTree("hv_booster_channels",
            "Booster registry (booster_id ↔ name)");
        t_bch->Branch("booster_id", &booster_id, "booster_id/s");
        t_bch->Branch("name",       name_buf,    "name/C");
        for (int i = 0; i < n_bst; ++i) {
            booster_id = static_cast<UShort_t>(i);
            std::snprintf(name_buf, sizeof(name_buf), "%s",
                          seg.booster_names[i].c_str());
            t_bch->Fill();
        }
        t_bch->Write("hv_booster_channels", TObject::kOverwrite);
    }

    // ── hv_meta tree ────────────────────────────────────────────────────
    {
        Double_t run_t_start_unix = t_start;
        Double_t run_t_end_unix   = t_end;
        Int_t    interval_ms      = seg.interval_ms;
        Int_t    n_channels       = n_ch;
        Int_t    n_boosters       = n_bst;

        TTree *t_meta = new TTree("hv_meta", "HV segment scalars");
        t_meta->Branch("run_t_start_unix", &run_t_start_unix, "run_t_start_unix/D");
        t_meta->Branch("run_t_end_unix",   &run_t_end_unix,   "run_t_end_unix/D");
        t_meta->Branch("interval_ms",      &interval_ms,      "interval_ms/I");
        t_meta->Branch("n_channels",       &n_channels,       "n_channels/I");
        t_meta->Branch("n_boosters",       &n_boosters,       "n_boosters/I");
        t_meta->Fill();
        t_meta->Write("hv_meta", TObject::kOverwrite);
    }

    f->Close();
    return true;
}


// ─────────────────────────────────────────────────────────────────────────
// summary printer
// ─────────────────────────────────────────────────────────────────────────
void print_summary(const std::string &path, const Summary &s)
{
    std::cout << "  → wrote HV trees to " << path << "\n";
    std::cout << "     hv                  : " << s.n_hv_snapshots
              << " entries × dv[" << s.n_channels << "] / v0set["
              << s.n_channels << "] (float32)\n";
    std::cout << "     hv_channels         : " << s.n_channels << " entries\n";
    if (s.n_boosters > 0) {
        std::cout << "     hv_booster          : " << s.n_booster_snaps
                  << " entries × vmon/imon/vset/iset[" << s.n_boosters
                  << "]\n";
        std::cout << "     hv_booster_channels : " << s.n_boosters
                  << " entries (";
        for (std::size_t i = 0; i < s.booster_names.size(); ++i) {
            if (i) std::cout << ", ";
            std::cout << s.booster_names[i];
        }
        std::cout << ")\n";
    }
    std::cout << "     hv_meta             : 1 entry (interval="
              << s.interval_ms << " ms)\n";
    std::cout << "     run window  : " << iso_utc(s.run_t_start_unix)
              << " → " << iso_utc(s.run_t_end_unix)
              << " (" << (s.run_t_end_unix - s.run_t_start_unix) << " s)\n";
    if (s.n_hv_snapshots > 0) {
        std::cout << "     event range : " << s.ev_first << " → "
                  << s.ev_last << " (from " << s.n_events
                  << " events in `" << s.event_tree_name
                  << "`, anchored to " << s.n_sync_pins
                  << " SYNC pin" << (s.n_sync_pins == 1 ? "" : "s") << ")\n";
    }
    std::cout << "     HV archive  : ";
    if (s.source_archive.empty()) std::cout << "(none)";
    else {
        for (std::size_t i = 0; i < s.source_archive.size(); ++i) {
            if (i) std::cout << ", ";
            std::cout << s.source_archive[i];
        }
    }
    std::cout << "\n";
    if (!s.overwrote.empty()) {
        std::cout << "     overwrote   : ";
        for (std::size_t i = 0; i < s.overwrote.size(); ++i) {
            if (i) std::cout << ", ";
            std::cout << s.overwrote[i];
        }
        std::cout << " (old cycles purged)\n";
    }
}

} // anonymous namespace


// ─────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    std::string hv_dir = kDefaultHvDir;
    std::string channels_csv;
    double      pad_s  = 30.0;
    bool        force  = false;

    int opt;
    while ((opt = getopt(argc, argv, "d:c:p:fh")) != -1) {
        switch (opt) {
        case 'd': hv_dir       = optarg;          break;
        case 'c': channels_csv = optarg;          break;
        case 'p': pad_s        = std::atof(optarg); break;
        case 'f': force        = true;            break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    std::vector<std::string> inputs;
    for (int i = optind; i < argc; ++i) inputs.push_back(argv[i]);
    if (inputs.empty()) {
        usage(argv[0]);
        return 1;
    }

    auto channels = split_comma(channels_csv);

    std::cout << "[replay_add_hv] " << inputs.size()
              << " file(s); HV archive: " << hv_dir << "\n";

    int n_ok = 0, n_err = 0;
    for (const auto &fp : inputs) {
        std::cout << "\n[file] " << fp << "\n";
        Summary s;
        std::string err;
        bool ok = false;
        try {
            ok = process_file(fp, hv_dir, channels, pad_s, force, s, err);
        } catch (const std::exception &e) {
            err = e.what();
        }
        if (!ok) {
            std::cerr << "  ✗ " << err << "\n";
            ++n_err;
            continue;
        }
        print_summary(fp, s);
        ++n_ok;
    }
    std::cout << "\n[replay_add_hv] done — " << n_ok << " ok, "
              << n_err << " failed\n";
    return n_err > 0 ? 1 : 0;
}
