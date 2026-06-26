//=============================================================================
// replay_filter.cpp — slow-control filter for replayed ROOT files
//
// Reads one or more replayed ROOT files (raw or recon), applies user-defined
// cuts on the slow streams (DSC2 livetime + EPICS values), and writes a
// single output ROOT file containing:
//   * the events / recon tree with only the physics events bracketed by
//     two adjacent "good" slow-event checkpoints,
//   * the scalers and epics trees concatenated from every input file (no
//     filtering — they are small and useful as run-wide context),
// plus a JSON report with one entry per (cut-channel, slow-event) point so
// downstream tools can plot per-channel value traces with pass/fail status
// and the cut acceptance band.
//
// Cuts JSON schema:
//   {
//     "livetime": {
//       "source":  "ref",            // "ref" | "trg" | "tdc"
//       "channel": 0,                // ignored for "ref"
//       "abs":     { "min": 90, "max": 100 },
//       "rel_rms": 3
//     },
//     "epics": {
//       "<channel_name>": { "abs": {...}, "rel_rms": 3 },
//       ...
//     }
//   }
//
// `rel_rms: N` accepts points within N · σ̂ of the channel's median, where
// σ̂ = 1.4826 · MAD (median absolute deviation).  MAD is robust to heavy
// outliers (one bad reading does not pull the centre or width).
//
// Optional "split" block — classify each event by a slow-control PV level
// (e.g. target cell pressure) and write the "full" and "empty" subsets to
// separate output files.  Every checkpoint is labelled by its PV reading:
//   PV >= full  → full-target        PV <= empty → empty-target
//   empty < PV < full (the ramp), or no reading  → dropped (lands in neither)
// so events taken while the cell is filling/emptying are excluded from both.
//   "split": {
//     "channel": "TGT:PRad:Cell_P",  // EPICS channel to watch
//     "full":  500,                  // PV >= this  → full-target file
//     "empty": 5,                    // PV <= this  → empty-target file
//     "guard_checkpoints": 0,        // optional ± margin dropped at each state edge
//     "labels": ["full", "empty"]    // output suffixes (defaults: full / empty)
//   }
// Classification is per-checkpoint and stateless, so any number of full<->empty
// transitions in one run are handled, and a run that only ever shows one state
// produces only that one file (a pure full run → just <stem>_full.root).  Each
// side still gets every other configured cut and its own report + charge.
//
// Output ROOT file(s):
//   * events / recon — same schema as input, only kept events
//   * scalers / epics — concatenated from every input plus an extra
//     `good` boolean branch per row reflecting that checkpoint's
//     overall verdict (all cuts pass); with split on, `good` also requires
//     the row to belong to that file's side, so live_charge on a side file
//     reproduces that side's post-cut charge directly.
// With split on, one file <stem>_<label>.root + matching report is written per
// target state actually seen (one or two), instead of the single default output.
// JSON report: see the write phase in the source for the full layout (a
// "split" block is added per side when run-splitting is active).
//=============================================================================

#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"     // analysis::get_run_int

#include <TFile.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <getopt.h>
#include <TROOT.h>

using json = nlohmann::json;

namespace {

// ── Cut configuration ────────────────────────────────────────────────────

struct AbsCut {
    bool   has_min = false;
    double min_val = 0;
    bool   has_max = false;
    double max_val = 0;
};

struct ChannelCut {
    AbsCut       abs;
    bool         has_rel_rms = false;
    double       rel_rms_n   = 0;
    // Optional: condition this channel's robust median/MAD on points where
    // every named gating channel's cut passed.  Use when the channel of
    // interest is physically meaningful only in a particular regime —
    // e.g. beam position is only meaningful when current is above some
    // floor.  Accepts either a single string or an array of strings in
    // the JSON; multiple gates are ANDed (a row's value is included in
    // the stats only if ALL listed gating channels passed).  All listed
    // channels must also be configured.  One level of gating is
    // supported (the gating channels themselves use ungated stats).
    // EPICS-on-EPICS only for now.
    std::vector<std::string> gated_by;

    // Robust statistics (filled in phase 2 if has_rel_rms).
    // `center` is the median; `sigma` = 1.4826 · MAD (consistent estimator
    // for a normal distribution, so `rel_rms: N` keeps its intuitive
    // "N standard deviations" meaning).  `n_used` is the input point count
    // (MAD doesn't iterate-and-drop); `n_clipped` is points outside
    // [center − N·sigma, center + N·sigma], reported for traceability.
    bool   stats_valid = false;
    double robust_center = 0;   // median
    double robust_sigma  = 0;   // 1.4826 * MAD
    double mad           = 0;
    int    n_used        = 0;
    int    n_clipped     = 0;
};

struct LivetimeCut {
    bool        enabled = false;
    std::string source  = "ref";
    int         channel = 0;
    ChannelCut  cut;
};

// Live-charge integration over kept slow-event intervals.  Disabled if no
// `charge` block is present in the cut JSON; otherwise sums
//   Σ live_fraction × Δt × ½(I_i + I_{i+1})
// over each adjacent pair of accepted checkpoints, where live_fraction
// is the slice-local DSC2 livetime at the right endpoint and I is the
// configured beam-current EPICS channel.  Output units are
// (beam_current_unit · seconds).
struct ChargeCfg {
    bool        enabled              = false;
    std::string beam_current_channel;   // EPICS channel name to read
};

// Run-splitting by a slow-control PV level (e.g. target cell pressure full vs
// empty).  When enabled the filter classifies every checkpoint by its
// forward-filled `channel` reading — PV >= `full` is full-target, PV <= `empty`
// is empty-target, anything in between (the ramp) or with no reading yet is
// dropped — and writes one output file per state that actually occurs, each
// still subject to every other configured cut and with its own charge integral.
//
// Classification is per-checkpoint and stateless: the in-between dead zone is
// itself the guard band (ramp events land in neither file), any number of
// full<->empty transitions are handled, and a run that only ever shows one
// state yields only that file.  `guard_checkpoints` optionally drops an extra
// ± N checkpoints adjacent to each state edge, for margin beyond the dead zone.
struct SplitCfg {
    bool        enabled  = false;
    std::string channel;                       // EPICS channel to watch

    // Level thresholds (require full > empty).  A checkpoint is full when its
    // PV >= full_thresh, empty when PV <= empty_thresh, dropped otherwise.
    double full_thresh  = 0.0;
    double empty_thresh = 0.0;

    // Optional extra margin: also drop the ± guard_checkpoints points on either
    // side of any state edge (0 = rely on the dead zone alone).
    int    guard_checkpoints = 0;

    // Output suffixes for the full / empty files (and their reports).
    std::string label_full  = "full";
    std::string label_empty = "empty";
};

struct CutConfig {
    LivetimeCut                       livetime;
    std::map<std::string, ChannelCut> epics;
    ChargeCfg                         charge;
    SplitCfg                          split;
    json                              raw;            // echoed in the report
};

void parse_abs(const json &j, AbsCut &a)
{
    if (!j.is_object()) return;
    if (j.contains("min") && j["min"].is_number()) {
        a.has_min = true; a.min_val = j["min"].get<double>();
    }
    if (j.contains("max") && j["max"].is_number()) {
        a.has_max = true; a.max_val = j["max"].get<double>();
    }
}

void parse_channel_cut(const json &j, ChannelCut &c)
{
    if (j.contains("abs")) parse_abs(j["abs"], c.abs);
    if (j.contains("rel_rms") && j["rel_rms"].is_number()) {
        c.has_rel_rms = true;
        c.rel_rms_n   = j["rel_rms"].get<double>();
    }
    if (j.contains("gated_by")) {
        const auto &g = j["gated_by"];
        if (g.is_string()) {
            c.gated_by.push_back(g.get<std::string>());
        } else if (g.is_array()) {
            for (const auto &item : g) {
                if (item.is_string()) c.gated_by.push_back(item.get<std::string>());
            }
        }
    }
}

bool load_cuts(const std::string &path, CutConfig &cfg)
{
    std::ifstream f(path);
    if (!f) {
        std::cerr << "replay_filter: cannot open cuts file: " << path << "\n";
        return false;
    }
    json j;
    try {
        j = json::parse(f, nullptr, true, /*allow_comments=*/true);
    } catch (json::parse_error &e) {
        std::cerr << "replay_filter: cuts JSON parse error: " << e.what() << "\n";
        return false;
    }
    cfg.raw = j;

    if (j.contains("livetime")) {
        const auto &lj = j["livetime"];
        cfg.livetime.enabled = true;
        if (lj.contains("source"))  cfg.livetime.source  = lj["source"].get<std::string>();
        if (lj.contains("channel")) cfg.livetime.channel = lj["channel"].get<int>();
        parse_channel_cut(lj, cfg.livetime.cut);
    }

    if (j.contains("epics") && j["epics"].is_object()) {
        for (auto it = j["epics"].begin(); it != j["epics"].end(); ++it) {
            ChannelCut c;
            parse_channel_cut(it.value(), c);
            cfg.epics[it.key()] = c;
        }
    }

    // Charge integration is opt-in: the cut JSON must name the EPICS
    // channel that carries the beam current.  Anything else (e.g. units)
    // is documented in the report so downstream consumers can convert.
    if (j.contains("charge") && j["charge"].is_object()) {
        const auto &cj = j["charge"];
        if (cj.contains("beam_current") && cj["beam_current"].is_string()) {
            cfg.charge.enabled = true;
            cfg.charge.beam_current_channel = cj["beam_current"].get<std::string>();
        }
    }

    // Run-splitting by PV level (opt-in: needs a `channel` plus `full` and
    // `empty` thresholds, with full > empty).
    if (j.contains("split") && j["split"].is_object()) {
        const auto &sj = j["split"];
        auto &s = cfg.split;
        if (sj.contains("channel") && sj["channel"].is_string())
            s.channel = sj["channel"].get<std::string>();
        const bool has_full  = sj.contains("full")  && sj["full"].is_number();
        const bool has_empty = sj.contains("empty") && sj["empty"].is_number();
        if (has_full)  s.full_thresh  = sj["full"].get<double>();
        if (has_empty) s.empty_thresh = sj["empty"].get<double>();
        if (sj.contains("guard_checkpoints") && sj["guard_checkpoints"].is_number())
            s.guard_checkpoints = std::max(0, sj["guard_checkpoints"].get<int>());
        if (sj.contains("labels") && sj["labels"].is_array()
            && sj["labels"].size() >= 2) {
            s.label_full  = sj["labels"][0].get<std::string>();
            s.label_empty = sj["labels"][1].get<std::string>();
        }
        // Need a channel and a well-ordered (full > empty) threshold pair.
        const bool ok = !s.channel.empty() && has_full && has_empty
                        && s.full_thresh > s.empty_thresh;
        s.enabled = ok;
        if (!s.channel.empty() && !ok) {
            if (!has_full || !has_empty)
                std::cerr << "replay_filter: split needs both `full` and `empty` "
                             "level thresholds — split disabled\n";
            else
                std::cerr << "replay_filter: split `full` (" << s.full_thresh
                          << ") must be > `empty` (" << s.empty_thresh
                          << ") — split disabled\n";
        }
    }
    return true;
}

// ── Robust statistics: median + MAD ──────────────────────────────────────
//
// Uses median absolute deviation (Hampel 1974).  More robust to heavy
// outliers than iterative sigma clipping: a single bad reading shifts the
// median negligibly and inflates MAD only via its own contribution.  The
// 1.4826 factor makes σ̂ a consistent estimator of stddev for a normal
// distribution, so cut thresholds in `rel_rms: N` keep their intuitive
// "N standard deviations" meaning.

struct RobustStats {
    double median    = 0;
    double mad       = 0;        // raw MAD
    double sigma     = 0;        // 1.4826 * MAD
    int    n_used    = 0;
    int    n_clipped = 0;        // points outside the cut band (informational)
};

double median_inplace(std::vector<double> &xs)
{
    if (xs.empty()) return 0;
    size_t n = xs.size();
    auto mid = xs.begin() + n / 2;
    std::nth_element(xs.begin(), mid, xs.end());
    double m = *mid;
    if ((n & 1u) == 0) {
        // even count: average mid with the largest of the lower half
        auto max_lo = std::max_element(xs.begin(), mid);
        m = 0.5 * (m + *max_lo);
    }
    return m;
}

RobustStats robust_mad(const std::vector<double> &xs, double n_sigma_for_clip_count)
{
    RobustStats r;
    if (xs.empty()) return r;
    r.n_used = static_cast<int>(xs.size());

    std::vector<double> tmp = xs;
    r.median = median_inplace(tmp);

    std::vector<double> dev;
    dev.reserve(xs.size());
    for (double x : xs) dev.push_back(std::fabs(x - r.median));
    r.mad   = median_inplace(dev);
    r.sigma = 1.4826 * r.mad;

    if (r.sigma > 0 && n_sigma_for_clip_count > 0) {
        for (double x : xs)
            if (std::fabs(x - r.median) > n_sigma_for_clip_count * r.sigma)
                ++r.n_clipped;
    }
    return r;
}

// ── In-memory slow-event rows ────────────────────────────────────────────

struct ScalerRow {
    int32_t  event_number   = 0;
    int64_t  ti_ticks       = 0;
    int64_t  unix_time      = 0;     // 0 if no SYNC seen yet
    uint32_t sync_counter   = 0;
    uint32_t run_number     = 0;
    uint32_t ref_gated      = 0;
    uint32_t ref_ungated    = 0;
    uint32_t trg_gated[16]   = {};
    uint32_t trg_ungated[16] = {};
    uint32_t tdc_gated[16]   = {};
    uint32_t tdc_ungated[16] = {};
};

struct EpicsArrival {
    int32_t                         event_number = 0;   // event_number_at_arrival
    int64_t                         ti_ticks     = 0;   // ti_ticks_at_arrival,
                                                        // 0 ⇒ not populated
                                                        // (legacy file or no
                                                        // physics seen yet)
    int64_t                         unix_time    = 0;
    uint32_t                        sync_counter = 0;
    uint32_t                        run_number   = 0;
    std::map<std::string, double>   updates;            // sparse
};

// Extract the (gated, ungated) pair the cut targets.  These are cumulative
// counters: the DSC2 increments them since run-start without resetting at
// each readout, so a single row gives the run-average live fraction, not
// the live fraction over the most recent slice.
inline std::pair<uint32_t, uint32_t>
select_scaler_pair(const ScalerRow &r, const LivetimeCut &cfg)
{
    if (cfg.source == "ref") return {r.ref_gated, r.ref_ungated};
    int c = std::clamp(cfg.channel, 0, 15);
    if (cfg.source == "trg") return {r.trg_gated[c], r.trg_ungated[c]};
    if (cfg.source == "tdc") return {r.tdc_gated[c], r.tdc_ungated[c]};
    return {0, 0};
}

// Per-row delta-livetime in percent, indexed by load-order position.
// Walks the rows in event-number order and divides the *change* in gated
// over the change in ungated since the previous reading — i.e. the live
// fraction over the slice between adjacent scaler readouts.  This is what
// quality cuts need: the run-cumulative ratio dilutes a recent dropout
// behind several minutes of good livetime.
//
// The implicit predecessor at run-start is (0, 0), so the first row's
// "delta" equals the cumulative readout over (run_start, first_readout].
// If the counter ever moves backward (DSC2 reset / wrap), the previous
// is rebased to (0, 0) at that row and the delta is taken from there.
// Slots where ungated did not advance return -1 (cannot compute).
std::vector<double>
compute_delta_live_pct(const std::vector<ScalerRow> &scalers,
                       const std::vector<size_t>    &sc_order,
                       const LivetimeCut            &cfg)
{
    std::vector<double> out(scalers.size(), -1.0);
    uint32_t prev_g = 0, prev_u = 0;
    for (size_t k = 0; k < sc_order.size(); ++k) {
        const size_t orig = sc_order[k];
        const auto [g, u] = select_scaler_pair(scalers[orig], cfg);
        // Counter went backward — treat as a reset and rebase the baseline.
        if (g < prev_g || u < prev_u) { prev_g = 0; prev_u = 0; }
        const uint32_t dg = g - prev_g;
        const uint32_t du = u - prev_u;
        if (du > 0 && dg <= du)
            out[orig] = static_cast<double>(dg) / static_cast<double>(du) * 100.0;
        prev_g = g;
        prev_u = u;
    }
    return out;
}

// ── Tree readers ─────────────────────────────────────────────────────────

// load_*() preserves the input-file order; sorting is done downstream via an
// index permutation so phase-5 re-reads (which iterate in input order) can
// look up each row's verdict by its load-order index.
bool load_scalers(const std::vector<std::string> &files, std::vector<ScalerRow> &out)
{
    prad2::RawScalerData sc;
    for (const auto &path : files) {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "replay_filter: cannot open " << path << "\n";
            return false;
        }
        TTree *t = dynamic_cast<TTree *>(f->Get("scalers"));
        if (!t) continue;
        prad2::SetScalerReadBranches(t, sc);
        Long64_t n = t->GetEntries();
        out.reserve(out.size() + n);
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            ScalerRow r;
            r.event_number = sc.event_number;
            r.ti_ticks     = sc.ti_ticks;
            r.unix_time    = sc.unix_time;
            r.sync_counter = sc.sync_counter;
            r.run_number   = sc.run_number;
            r.ref_gated    = sc.ref_gated;
            r.ref_ungated  = sc.ref_ungated;
            std::memcpy(r.trg_gated,   sc.trg_gated,   16 * sizeof(uint32_t));
            std::memcpy(r.trg_ungated, sc.trg_ungated, 16 * sizeof(uint32_t));
            std::memcpy(r.tdc_gated,   sc.tdc_gated,   16 * sizeof(uint32_t));
            std::memcpy(r.tdc_ungated, sc.tdc_ungated, 16 * sizeof(uint32_t));
            out.push_back(r);
        }
    }
    return true;
}

bool load_epics(const std::vector<std::string> &files, std::vector<EpicsArrival> &out)
{
    for (const auto &path : files) {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "replay_filter: cannot open " << path << "\n";
            return false;
        }
        TTree *t = dynamic_cast<TTree *>(f->Get("epics"));
        if (!t) continue;

        prad2::RawEpicsData ep;
        std::vector<std::string> *cp = &ep.channel;
        std::vector<double>      *vp = &ep.value;
        t->SetBranchAddress("event_number_at_arrival", &ep.event_number_at_arrival);
        // ti_ticks_at_arrival was added after the first replays were taken;
        // tolerate its absence so legacy ROOT files still load (we will
        // fall back to the events-tree lookup for those rows).
        const bool has_ticks_at_arrival =
            (t->GetBranch("ti_ticks_at_arrival") != nullptr);
        if (has_ticks_at_arrival)
            t->SetBranchAddress("ti_ticks_at_arrival", &ep.ti_ticks_at_arrival);
        t->SetBranchAddress("unix_time",    &ep.unix_time);
        t->SetBranchAddress("sync_counter", &ep.sync_counter);
        t->SetBranchAddress("run_number",   &ep.run_number);
        t->SetBranchAddress("channel", &cp);
        t->SetBranchAddress("value",   &vp);

        Long64_t n = t->GetEntries();
        out.reserve(out.size() + n);
        for (Long64_t i = 0; i < n; ++i) {
            ep.ti_ticks_at_arrival = 0;
            t->GetEntry(i);
            EpicsArrival a;
            a.event_number = ep.event_number_at_arrival;
            a.ti_ticks     = has_ticks_at_arrival
                              ? static_cast<int64_t>(ep.ti_ticks_at_arrival)
                              : 0;
            a.unix_time    = ep.unix_time;
            a.sync_counter = ep.sync_counter;
            a.run_number   = ep.run_number;
            size_t k_max = std::min(ep.channel.size(), ep.value.size());
            for (size_t k = 0; k < k_max; ++k)
                a.updates[ep.channel[k]] = ep.value[k];
            out.push_back(std::move(a));
        }
    }
    return true;
}

// Index permutation that orders the input vector by event_number.
template <class T>
std::vector<size_t> sort_index_by_event(const std::vector<T> &v)
{
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) { return v[a].event_number < v[b].event_number; });
    return idx;
}

// Pre-scan the events/recon tree across all input files and build a
// lookup event_num → ti_ticks (the 48-bit TI timestamp).  Only `event_num`
// and `timestamp` branches are activated, so this is fast even on
// millions-of-event runs.  Used to pin each report point to the TI tick
// of the physics event it is associated with.
bool build_evn_to_ticks(const std::vector<std::string> &files,
                        const std::string              &tree_name,
                        std::map<int32_t, int64_t>     &out)
{
    int       event_num = 0;
    long long timestamp = 0;
    for (const auto &path : files) {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "replay_filter: cannot open " << path << "\n";
            return false;
        }
        TTree *t = dynamic_cast<TTree *>(f->Get(tree_name.c_str()));
        if (!t) continue;

        // Activate just the two branches we need.
        t->SetBranchStatus("*", 0);
        if (auto *b = t->GetBranch("event_num")) {
            t->SetBranchStatus("event_num", 1);
            t->SetBranchAddress("event_num", &event_num);
        } else {
            std::cerr << "replay_filter: '" << tree_name
                      << "' has no event_num branch in " << path << "\n";
            return false;
        }
        if (auto *b = t->GetBranch("timestamp")) {
            t->SetBranchStatus("timestamp", 1);
            t->SetBranchAddress("timestamp", &timestamp);
        } else {
            std::cerr << "replay_filter: '" << tree_name
                      << "' has no timestamp branch in " << path << "\n";
            return false;
        }

        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            // First-write-wins: in case of duplicate event_num across files
            // (shouldn't happen for a single run), keep the earliest tick.
            out.emplace(event_num, static_cast<int64_t>(timestamp));
        }
    }
    return true;
}

// ── Cut evaluation ───────────────────────────────────────────────────────

bool eval_channel_cut(const ChannelCut &c, double value)
{
    if (c.abs.has_min && !(value >= c.abs.min_val)) return false;
    if (c.abs.has_max && !(value <= c.abs.max_val)) return false;
    if (c.has_rel_rms && c.stats_valid && c.robust_sigma > 0) {
        if (std::fabs(value - c.robust_center) > c.rel_rms_n * c.robust_sigma)
            return false;
    }
    return true;
}

// ── Report point ─────────────────────────────────────────────────────────
//
// Each report point carries:
//   * `associated_evn` — the physics event the slow row is anchored to.
//     Scaler rows: their own event_number (the SYNC physics event whose
//     readout included the scaler bank).  EPICS rows: event_number_at_-
//     arrival (the most recent physics event seen at the time of the
//     EPICS event).  Both are integer keys into the events/recon tree.
//   * `associated_timestamp` (relative seconds) — the TI 48-bit tick of
//     the physics event with event_num == associated_evn, in seconds
//     since the earliest looked-up TI tick.  Looked up once per unique
//     event_number from the events/recon tree, so it is the *exact* time
//     of the physics event the slow row is tied to (no forward-fill /
//     SYNC-interval lag).  null when no physics event matches (e.g.
//     event_number_at_arrival = -1 — EPICS arrived before any physics).
//   * `unix_time` (absolute Unix seconds) — from the 0xE112 HEAD bank.
//     Native for EPICS rows.  Explicitly null for scaler rows: the
//     scaler's cached unix_time can be a SYNC interval old, and emitting
//     it would invite mis-alignment.  Charts that need absolute time
//     should plot associated_timestamp and use any EPICS unix_time as
//     the absolute anchor (one EPICS pin is enough for the whole run).
struct ReportPoint {
    std::string channel;
    bool        pass;
    int32_t     event_number;     // = associated_evn
    bool        has_assoc_t;
    double      assoc_t_rel;      // seconds since the run's earliest event
    bool        has_unix_time;
    int64_t     unix_time;
    double      value;            // NaN ⇒ value not yet seen
    int         side = 0;         // split side this checkpoint belongs to;
                                  // -1 = dropped (ramp/guard). Backfilled in
                                  // phase 5, stays 0 when split is off.
};

// ── Main pipeline ────────────────────────────────────────────────────────

bool detect_event_tree(const std::string &path, std::string &name)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return false;
    if (f->Get("events")) { name = "events"; return true; }
    if (f->Get("recon"))  { name = "recon";  return true; }
    return false;
}

Long64_t tree_entries(const std::string &path, const char *tree_name)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return 0;
    TTree *t = dynamic_cast<TTree *>(f->Get(tree_name));
    return t ? t->GetEntries() : 0;
}

std::string insert_before_root(const std::string &path, const std::string &suffix)
{
    auto p = std::filesystem::path(path);
    std::string name = p.filename().string();
    const std::string ext = ".root";
    if (name.size() >= ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0)
        name.insert(name.size() - ext.size(), suffix);
    else
        name += suffix;
    return name;
}

std::string filtered_output_path(const std::string &out_dir,
                                 const std::string &input_path,
                                 const std::string &label = "")
{
    std::string suffix = "_filter" + (label.empty() ? std::string() : "_" + label);
    return (std::filesystem::path(out_dir) / insert_before_root(input_path, suffix)).string();
}

std::string run_epics_path(const std::string &out_dir, int run_number)
{
    char name[64];
    std::snprintf(name, sizeof(name), "prad_%06d_epics.root", run_number);
    return (std::filesystem::path(out_dir) / name).string();
}

std::string run_report_path(const std::string &out_dir, int run_number)
{
    char name[80];
    std::snprintf(name, sizeof(name), "prad_%06d_filter_report.json", run_number);
    return (std::filesystem::path(out_dir) / name).string();
}

int run(const std::vector<std::string> &input_files,
        const std::string &output_path,
        const std::string &cuts_path,
        const std::string &report_path,
        int run_number_override,
        int num_threads)
{
    CutConfig cuts;
    if (!load_cuts(cuts_path, cuts)) return 1;

    // ---------- Phase 1: load slow streams into memory ----------
    std::vector<ScalerRow>    scalers;
    std::vector<EpicsArrival> epics_rows;
    if (!load_scalers(input_files, scalers))   return 1;
    if (!load_epics  (input_files, epics_rows)) return 1;
    std::cerr << "replay_filter: loaded " << scalers.size() << " scaler + "
              << epics_rows.size() << " epics rows from "
              << input_files.size() << " file(s)\n";

    int run_number = run_number_override;
    if (run_number < 0) run_number = analysis::get_run_int(input_files.front());
    if (run_number < 0 && !scalers.empty())    run_number = (int)scalers.front().run_number;
    if (run_number < 0 && !epics_rows.empty()) run_number = (int)epics_rows.front().run_number;

    struct FileInfo {
        std::string path;
        size_t scaler_offset = 0;
        size_t epics_offset = 0;
        Long64_t scaler_entries = 0;
        Long64_t epics_entries = 0;
    };
    std::vector<FileInfo> file_info;
    file_info.reserve(input_files.size());
    size_t scaler_offset = 0, epics_offset = 0;
    for (const auto &path : input_files) {
        FileInfo info;
        info.path = path;
        info.scaler_offset = scaler_offset;
        info.epics_offset = epics_offset;
        info.scaler_entries = tree_entries(path, "scalers");
        info.epics_entries = tree_entries(path, "epics");
        scaler_offset += static_cast<size_t>(info.scaler_entries);
        epics_offset += static_cast<size_t>(info.epics_entries);
        file_info.push_back(std::move(info));
    }

    // Sort scalers once and precompute delta livetime per row.  Cuts evaluate
    // and report against the slice-local live fraction (Δgated / Δungated),
    // not the run-cumulative ratio cached on each row.
    auto sc_order = sort_index_by_event(scalers);
    std::vector<double> delta_live_pct;
    if (cuts.livetime.enabled)
        delta_live_pct = compute_delta_live_pct(scalers, sc_order, cuts.livetime);

    // ---------- Phase 2: robust stats for rel_rms cuts ----------
    // Ungated channels first (stats from all values), then gated channels
    // (stats restricted to rows where the named gating channel's cut
    // passed).  Gating is one-level: the gating channel itself uses its
    // own ungated stats — gating chains are intentionally not supported
    // to keep the JSON unambiguous.
    auto fill_stats = [&](ChannelCut &c, const std::vector<double> &xs) {
        auto rs = robust_mad(xs, c.rel_rms_n);
        c.stats_valid   = (rs.n_used > 1) && rs.sigma > 0;
        c.robust_center = rs.median;
        c.robust_sigma  = rs.sigma;
        c.mad           = rs.mad;
        c.n_used        = rs.n_used;
        c.n_clipped     = rs.n_clipped;
    };

    // Walk EPICS rows in event-number order so forward-fill of the gating
    // channel reflects the actual time sequence.
    auto ep_order_for_stats = sort_index_by_event(epics_rows);

    // 1. livetime (independent of EPICS).
    if (cuts.livetime.enabled && cuts.livetime.cut.has_rel_rms) {
        std::vector<double> xs;
        xs.reserve(scalers.size());
        for (size_t i = 0; i < scalers.size(); ++i) {
            double v = delta_live_pct[i];
            if (v >= 0) xs.push_back(v);
        }
        fill_stats(cuts.livetime.cut, xs);
    }

    // 2. ungated EPICS channels (stats from all observed values).
    for (auto &kv : cuts.epics) {
        if (!kv.second.has_rel_rms) continue;
        if (!kv.second.gated_by.empty()) continue;
        std::vector<double> xs;
        for (size_t oi : ep_order_for_stats) {
            auto it = epics_rows[oi].updates.find(kv.first);
            if (it != epics_rows[oi].updates.end()) xs.push_back(it->second);
        }
        fill_stats(kv.second, xs);
    }

    // 3. gated EPICS channels (stats from rows where every gate's cut passed).
    for (auto &kv : cuts.epics) {
        if (!kv.second.has_rel_rms) continue;
        if (kv.second.gated_by.empty()) continue;

        // Resolve every gating channel.  If any is missing, fall back to
        // ungated stats (and log) — partial gating would be misleading.
        std::vector<const ChannelCut *> gates;
        gates.reserve(kv.second.gated_by.size());
        bool gates_ok = true;
        for (const auto &gname : kv.second.gated_by) {
            auto it = cuts.epics.find(gname);
            if (it == cuts.epics.end()) {
                std::cerr << "replay_filter: channel '" << kv.first
                          << "' is gated_by '" << gname
                          << "' which is not configured — falling back to ungated stats\n";
                gates_ok = false;
                break;
            }
            if (!it->second.gated_by.empty()) {
                std::cerr << "replay_filter: channel '" << kv.first
                          << "' gated_by '" << gname
                          << "' which is itself gated — chains not supported, "
                             "ignoring the inner gating\n";
            }
            gates.push_back(&it->second);
        }
        if (!gates_ok) {
            std::vector<double> xs;
            for (size_t oi : ep_order_for_stats) {
                auto it = epics_rows[oi].updates.find(kv.first);
                if (it != epics_rows[oi].updates.end()) xs.push_back(it->second);
            }
            fill_stats(kv.second, xs);
            continue;
        }

        std::vector<double> xs;
        std::map<std::string, double> cur_eps;        // forward-fill across rows
        for (size_t oi : ep_order_for_stats) {
            const auto &row = epics_rows[oi];
            for (const auto &up : row.updates) cur_eps[up.first] = up.second;

            auto val_it = row.updates.find(kv.first);
            if (val_it == row.updates.end()) continue;

            bool all_pass = true;
            for (size_t gi = 0; gi < gates.size(); ++gi) {
                auto gv_it = cur_eps.find(kv.second.gated_by[gi]);
                if (gv_it == cur_eps.end()
                    || !eval_channel_cut(*gates[gi], gv_it->second)) {
                    all_pass = false;
                    break;
                }
            }
            if (!all_pass) continue;
            xs.push_back(val_it->second);
        }
        fill_stats(kv.second, xs);
    }

    // ---------- Phase 3: anchor for relative associated_timestamp ----------
    // Slow rows now carry their own TI tick: scalers via `ti_ticks` (from the
    // carrying SYNC event's info.timestamp) and EPICS via `ti_ticks` (the new
    // ti_ticks_at_arrival branch, captured at decode time so it is independent
    // of whether the anchor event was written to the events tree).  We still
    // detect and pre-scan the events tree below — needed for phase 5's
    // physics-filter loop, and as a back-compat fallback for EPICS rows from
    // legacy replays that lack the new branch.
    std::string ev_tree_name;
    if (!detect_event_tree(input_files.front(), ev_tree_name)) {
        std::cerr << "replay_filter: no events/recon tree in "
                  << input_files.front() << "\n";
        return 1;
    }
    std::map<int32_t, int64_t> evn_to_ticks;
    if (!build_evn_to_ticks(input_files, ev_tree_name, evn_to_ticks)) return 1;

    // Anchor = smallest TI tick across every source we have.  Considering
    // the slow rows (not just the events tree) keeps anchor monotonicity
    // when the events tree skips early events (e.g. trigger filter).
    int64_t  ti_anchor    = 0;
    bool     anchor_set   = false;
    auto consider_tick = [&](int64_t t) {
        if (t <= 0) return;
        if (!anchor_set || t < ti_anchor) { ti_anchor = t; anchor_set = true; }
    };
    for (const auto &kv : evn_to_ticks) consider_tick(kv.second);
    for (const auto &s  : scalers)      consider_tick(s.ti_ticks);
    for (const auto &e  : epics_rows)   consider_tick(e.ti_ticks);
    constexpr double TI_TICK_SEC = 4.0e-9;

    // ---------- Phase 4: walk merged timeline, mark good/bad ----------
    // Iterate via index permutations so the parallel verdict vectors stay
    // aligned with the load-order vectors (used in phase 6).  sc_order was
    // built earlier so the delta-livetime precompute could share it.
    auto ep_order = sort_index_by_event(epics_rows);

    std::vector<bool> scaler_verdict(scalers.size(),  false);
    std::vector<bool> epics_verdict (epics_rows.size(), false);

    struct Checkpoint {
        int32_t event_number;
        int64_t unix_time;
        int64_t ti_ticks;        // 0 if unknown — pair contributes no charge
        double  live_fraction;   // [0, 1] from cur_lt/100, NaN if unset
        double  beam_current;    // forward-filled, NaN if not seen yet
        bool    overall_pass;
        double  split_pv;        // forward-filled split channel, NaN if unset
        bool    is_scaler;       // true ⇒ orig indexes `scalers`, else `epics`
        size_t  orig;            // load-order index into the source vector
    };
    std::vector<Checkpoint>  timeline;
    std::vector<ReportPoint> report_points;

    // Forward-fill state for cut evaluation only.
    double                            cur_lt    = -1.0;   // % livetime
    std::map<std::string, double>     cur_eps;
    int64_t                           last_unix = 0;

    size_t i_sc = 0, i_ep = 0;
    while (i_sc < sc_order.size() || i_ep < ep_order.size()) {
        const bool take_sc =
            (i_sc < sc_order.size()) &&
            (i_ep >= ep_order.size() ||
             scalers[sc_order[i_sc]].event_number <=
             epics_rows[ep_order[i_ep]].event_number);

        int32_t cp_evn   = 0;
        int64_t cp_unix  = 0;
        int64_t cp_ticks = 0;        // TI tick captured on the slow row itself
        size_t  orig     = 0;
        bool    is_sc    = take_sc;

        bool emit_unix = false;     // true only for EPICS rows
        if (take_sc) {
            orig = sc_order[i_sc++];
            const auto &s = scalers[orig];
            cp_evn   = s.event_number;
            cp_ticks = s.ti_ticks;   // SYNC event's own info.timestamp
            // Slice-local live fraction (Δgated / Δungated) — see
            // compute_delta_live_pct for why the cumulative row value is
            // not used.  The first row's predecessor is (0, 0).
            cur_lt = cuts.livetime.enabled ? delta_live_pct[orig] : -1.0;
            // Scaler's cached unix_time is intentionally ignored — it lags
            // by up to a SYNC interval and confuses alignment.  Charts that
            // need absolute time should use the EPICS unix_time pins.
            cp_unix = last_unix;
        } else {
            orig = ep_order[i_ep++];
            const auto &e = epics_rows[orig];
            cp_evn   = e.event_number;
            cp_ticks = e.ti_ticks;   // ti_ticks_at_arrival, captured at decode
            for (const auto &kv : e.updates) cur_eps[kv.first] = kv.second;
            if (e.unix_time > 0) last_unix = e.unix_time;
            cp_unix    = last_unix;
            emit_unix  = (e.unix_time > 0);
            // Legacy fallback: if the EPICS row was written before the
            // ti_ticks_at_arrival branch existed, recover the timestamp
            // by joining on event_number_at_arrival.  Misses when the
            // anchor event itself was filtered out at replay time —
            // exactly the failure mode the new branch closes.
            if (cp_ticks <= 0 && cp_evn >= 0) {
                auto it = evn_to_ticks.find(cp_evn);
                if (it != evn_to_ticks.end()) cp_ticks = it->second;
            }
        }

        // associated_timestamp: TI tick on this row, expressed as seconds
        // since the run's earliest seen tick.  null when the row carries
        // no tick (e.g. an EPICS event arriving before the first physics
        // event seen on the channel).
        bool   pt_has_t = (cp_ticks > 0) && anchor_set;
        double pt_t     = pt_has_t
                          ? (cp_ticks - ti_anchor) * TI_TICK_SEC : 0.0;
        bool   pt_has_unix = emit_unix;
        int64_t pt_unix    = emit_unix ? (int64_t)last_unix : 0;

        // Per-channel report points with forward-filled values (dense traces
        // for plotting).
        if (cuts.livetime.enabled) {
            bool has  = (cur_lt >= 0);
            double v  = has ? cur_lt : std::numeric_limits<double>::quiet_NaN();
            bool pass = has && eval_channel_cut(cuts.livetime.cut, cur_lt);
            report_points.push_back({"livetime", pass, cp_evn,
                                     pt_has_t, pt_t,
                                     pt_has_unix, pt_unix, v});
        }
        for (const auto &kv : cuts.epics) {
            auto   it   = cur_eps.find(kv.first);
            bool   has  = (it != cur_eps.end());
            double v    = has ? it->second
                              : std::numeric_limits<double>::quiet_NaN();
            bool   pass = has && eval_channel_cut(kv.second, v);
            report_points.push_back({"epics:" + kv.first, pass, cp_evn,
                                     pt_has_t, pt_t,
                                     pt_has_unix, pt_unix, v});
        }

        // Overall verdict at this checkpoint = AND of every configured cut.
        // Channels that haven't reported yet count as "fail" — the user's
        // spec says events bracketed by an undefined endpoint are dropped.
        bool overall = true;
        if (cuts.livetime.enabled) {
            if (cur_lt < 0 || !eval_channel_cut(cuts.livetime.cut, cur_lt))
                overall = false;
        }
        for (const auto &kv : cuts.epics) {
            auto it = cur_eps.find(kv.first);
            if (it == cur_eps.end() || !eval_channel_cut(kv.second, it->second)) {
                overall = false;
            }
        }
        // Live fraction at this checkpoint (forward-filled %, scaled to
        // [0, 1]); beam current pulled from the configured EPICS channel
        // also via forward-fill.  Missing values stay NaN so the charge
        // integration knows to skip the surrounding pair.
        const double cp_live_fraction = (cur_lt >= 0)
            ? cur_lt * 0.01 : std::numeric_limits<double>::quiet_NaN();
        double cp_current = std::numeric_limits<double>::quiet_NaN();
        if (cuts.charge.enabled) {
            auto it = cur_eps.find(cuts.charge.beam_current_channel);
            if (it != cur_eps.end()) cp_current = it->second;
        }
        // Split PV (forward-filled like the others) so the transition scan
        // below sees a dense trace even on checkpoints carrying no update.
        double cp_split_pv = std::numeric_limits<double>::quiet_NaN();
        if (cuts.split.enabled) {
            auto it = cur_eps.find(cuts.split.channel);
            if (it != cur_eps.end()) cp_split_pv = it->second;
        }

        timeline.push_back({cp_evn, cp_unix, cp_ticks,
                            cp_live_fraction, cp_current, overall,
                            cp_split_pv, is_sc, orig});
        if (is_sc) scaler_verdict[orig] = overall;
        else       epics_verdict [orig] = overall;
    }

    // ---------- Phase 5: classify each checkpoint by PV level ----------
    // When `split` is enabled every checkpoint is labelled by its forward-
    // filled split PV: full (side 0, PV >= full_thresh), empty (side 1,
    // PV <= empty_thresh), or dropped (-1) for the in-between ramp.  The dead
    // zone between the thresholds is itself the guard band; `guard_checkpoints`
    // optionally drops an extra ± margin at each state edge.  Classification is
    // stateless, so any number of full<->empty transitions are handled and a
    // run showing only one state yields only that file.  side 0 ↔ label_full,
    // side 1 ↔ label_empty.  With split off, side[] is 0 everywhere so the
    // keep/charge/output code below is one path that runs for the present
    // side(s).  `split_active` mirrors cuts.split.enabled but degrades to false
    // if the channel never reports / never reaches either level.
    bool split_active = cuts.split.enabled;
    std::vector<int> side(timeline.size(), 0);
    std::vector<int> scaler_side(scalers.size(), 0);
    std::vector<int> epics_side (epics_rows.size(), 0);

    int  n_state_transitions = 0;
    json transitions = json::array();   // {evn, from, to, timestamp} per edge

    if (split_active) {
        const auto &S = cuts.split;
        // After phase 3's forward-fill the only NaNs are the head before the
        // channel's first report; back-fill that head with the first reading
        // (the run's starting state) so a pure run keeps its leading events.
        double first_val = std::numeric_limits<double>::quiet_NaN();
        for (const auto &cp : timeline)
            if (std::isfinite(cp.split_pv)) { first_val = cp.split_pv; break; }

        if (!std::isfinite(first_val)) {
            std::cerr << "replay_filter: split channel '" << S.channel
                      << "' never reported — cannot classify, writing the single "
                         "unsplit output instead\n";
            split_active = false;
        } else {
            // Raw per-checkpoint level classification.
            for (size_t k = 0; k < timeline.size(); ++k) {
                double v = timeline[k].split_pv;
                if (!std::isfinite(v)) v = first_val;          // back-fill head
                side[k] = (v >= S.full_thresh)  ? 0
                        : (v <= S.empty_thresh) ? 1
                        :                         -1;          // ramp / dead zone
            }
            // Optional extra margin: drop ± guard_checkpoints around any edge
            // between two adjacent, different, non-dropped states (a transition
            // faster than the checkpoint spacing leaves no dead-zone point
            // between the two levels).  Measured against the raw labels.
            if (S.guard_checkpoints > 0) {
                const std::vector<int> raw = side;
                const int N = S.guard_checkpoints;
                for (size_t k = 1; k < raw.size(); ++k)
                    if (raw[k] >= 0 && raw[k - 1] >= 0 && raw[k] != raw[k - 1])
                        for (int j = std::max(0, (int)k - N);
                             j < std::min((int)raw.size(), (int)k + N); ++j)
                            side[j] = -1;
            }
            // Count full<->empty transitions (using the last non-dropped state,
            // so a dead-zone gap is one transition not two) and record each
            // edge at the first checkpoint of the new state.
            int prev_state = -1;
            for (size_t k = 0; k < timeline.size(); ++k) {
                const int sd = side[k];
                if (sd < 0) continue;
                if (prev_state >= 0 && sd != prev_state) {
                    ++n_state_transitions;
                    const auto &cp = timeline[k];
                    const bool has_t = (cp.ti_ticks > 0 && anchor_set);
                    transitions.push_back({
                        {"evn",       cp.event_number},
                        {"from",      prev_state == 0 ? "full" : "empty"},
                        {"to",        sd == 0 ? "full" : "empty"},
                        {"timestamp", has_t
                            ? json((cp.ti_ticks - ti_anchor) * TI_TICK_SEC)
                            : json(nullptr)},
                    });
                }
                prev_state = sd;
            }
            // No checkpoint reached either level (PV sat in the dead zone the
            // whole run, or thresholds don't bracket the data) → nothing to
            // split.  Degrade to a single unsplit output rather than emit
            // empty files.
            int nf = 0, ne = 0;
            for (int sd : side) { if (sd == 0) ++nf; else if (sd == 1) ++ne; }
            if (nf == 0 && ne == 0) {
                std::cerr << "replay_filter: split channel '" << S.channel
                          << "' never reached full (>=" << S.full_thresh
                          << ") or empty (<=" << S.empty_thresh << ") — no split, "
                             "writing the single unsplit output instead\n";
                std::fill(side.begin(), side.end(), 0);
                split_active = false;
            }
        }

        // Back-map to load-order rows + tag report points with their side.
        if (split_active) {
            for (size_t k = 0; k < timeline.size(); ++k) {
                if (timeline[k].is_scaler) scaler_side[timeline[k].orig] = side[k];
                else                       epics_side [timeline[k].orig] = side[k];
            }
            const size_t ppc = (cuts.livetime.enabled ? 1u : 0u) + cuts.epics.size();
            if (ppc > 0 && report_points.size() == timeline.size() * ppc)
                for (size_t k = 0; k < timeline.size(); ++k)
                    for (size_t j = 0; j < ppc; ++j)
                        report_points[k * ppc + j].side = side[k];
        }
    }

    // ---------- Phase 5: build per-side keep-intervals (lo, hi] + charge ----------
    // A pair (cp_{i-1}, cp_i) belongs to a side only when both endpoints carry
    // the same non-guard label; pairs that straddle the transition or touch
    // the guard contribute to neither output (those events are dropped).  The
    // charge integration is the same arithmetic as before, bucketed per side:
    //   * gated (live_charge[s]) — both endpoints overall_pass.  Canonical
    //     post-cut number; matches the events written to that side's file.
    //   * ungated (ungated_*[s]) — every valid-data pair on side s regardless
    //     of the cut verdict, so users see how much charge the cuts dropped.
    std::vector<std::pair<int32_t, int32_t>> keep[2];
    std::vector<std::pair<int32_t, int32_t>> span[2];   // ungated per-side ranges
    double live_charge[2]              = {0, 0};
    double live_charge_secs[2]         = {0, 0};
    double real_secs[2]                = {0, 0};
    double ungated_live_charge[2]      = {0, 0};
    double ungated_live_charge_secs[2] = {0, 0};
    double ungated_real_secs[2]        = {0, 0};
    int    n_charge_pairs[2]           = {0, 0};
    int    n_charge_skipped[2]         = {0, 0};
    int    n_ungated_charge_pairs[2]   = {0, 0};
    int    n_ungated_charge_skipped[2] = {0, 0};
    for (size_t i = 1; i < timeline.size(); ++i) {
        const auto &a = timeline[i - 1];
        const auto &b = timeline[i];
        const int sa = side[i - 1], sb = side[i];
        const int ps = (sa >= 0 && sa == sb) ? sa : -1;   // pair's side, or none
        if (ps < 0) continue;                              // straddles / guard
        span[ps].emplace_back(a.event_number, b.event_number);  // ungated range
        const bool good_pair = (a.overall_pass && b.overall_pass);
        if (good_pair)
            keep[ps].emplace_back(a.event_number, b.event_number);
        if (!cuts.charge.enabled) continue;
        const bool data_ok = !(a.ti_ticks <= 0 || b.ti_ticks <= 0 || b.ti_ticks <= a.ti_ticks
            || !std::isfinite(b.live_fraction) || b.live_fraction < 0
            || !std::isfinite(a.beam_current)  || !std::isfinite(b.beam_current));
        if (!data_ok) {
            if (good_pair) ++n_charge_skipped[ps];
            ++n_ungated_charge_skipped[ps];
            continue;
        }
        const double dt = (b.ti_ticks - a.ti_ticks) * TI_TICK_SEC;
        const double I  = 0.5 * (a.beam_current + b.beam_current);
        const double dQ = b.live_fraction * dt * I;
        const double dL = b.live_fraction * dt;
        ungated_live_charge[ps]      += dQ;
        ungated_live_charge_secs[ps] += dL;
        ungated_real_secs[ps]        += dt;
        ++n_ungated_charge_pairs[ps];
        if (good_pair) {
            live_charge[ps]      += dQ;
            live_charge_secs[ps] += dL;
            real_secs[ps]        += dt;
            ++n_charge_pairs[ps];
        }
    }
    auto in_intervals = [](const std::vector<std::pair<int32_t, int32_t>> &iv,
                           int32_t evn) -> bool {
        if (iv.empty()) return false;
        auto it = std::upper_bound(
            iv.begin(), iv.end(), evn,
            [](int32_t e, const std::pair<int32_t, int32_t> &p) { return e < p.first; });
        if (it == iv.begin()) return false;
        --it;
        return evn > it->first && evn <= it->second;
    };
    auto is_kept = [&](int s, int32_t evn) { return in_intervals(keep[s], evn); };
    auto in_span = [&](int s, int32_t evn) { return in_intervals(span[s], evn); };

    // ---------- Phase 6: write the output(s) — one ROOT file + report per side ----------
    // ev_tree_name was detected in phase 3 above.  `write_side` does the whole
    // job for one side; with split off it runs once (side 0 → output_path),
    // with split on it runs once per present target state with the file/report
    // stems suffixed by that state's label.  A row's `good` flag in a side's
    // slow trees is its overall verdict AND-ed with "belongs to this side", so
    // running live_charge on a
    // side file reproduces that side's post-cut charge directly.
    const bool is_recon = (ev_tree_name == "recon");

    // Checkpoints per state (n_cp_side[0]=full, [1]=empty) and dropped
    // (ramp / dead zone / margin) — reported in the split block, and the
    // present-state check (>0) that drives which files get written.
    int n_cp_side[2] = {0, 0}, n_cp_guard = 0;
    for (int sd : side) { if (sd < 0) ++n_cp_guard; else ++n_cp_side[sd]; }

    auto fmt_pct = [](double r) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(2) << (r * 100.0) << "%";
        return o.str();
    };

    if (input_files.size() > 1) {
        std::filesystem::create_directories(output_path);
        std::vector<size_t> all_file_indices(input_files.size());
        for (size_t i = 0; i < all_file_indices.size(); ++i) all_file_indices[i] = i;

        auto write_slow_trees = [&](TFile &out, const std::vector<size_t> &file_indices,
                                    int s, bool restrict_side) -> int {
            {
                prad2::RawScalerData sc;
                bool good = false;
                out.cd();
                TTree *out_sc = new TTree("scalers", "PRad2 DSC2 scaler readouts (concatenated)");
                prad2::SetScalerWriteBranches(out_sc, sc);
                out_sc->Branch("good", &good, "good/O");
                for (size_t fi : file_indices) {
                    std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                    TTree *t = f ? dynamic_cast<TTree *>(f->Get("scalers")) : nullptr;
                    if (!t) continue;
                    prad2::SetScalerReadBranches(t, sc);
                    size_t seq = file_info[fi].scaler_offset;
                    Long64_t n = t->GetEntries();
                    for (Long64_t i = 0; i < n; ++i, ++seq) {
                        t->GetEntry(i);
                        good = (seq < scaler_verdict.size()) ? scaler_verdict[seq] : false;
                        if (good && restrict_side)
                            good = (seq < scaler_side.size() && scaler_side[seq] == s);
                        out.cd();
                        out_sc->Fill();
                    }
                }
                out.cd();
                out_sc->Write();
            }
            {
                prad2::RawEpicsData ep;
                std::vector<std::string> *cp = &ep.channel;
                std::vector<double> *vp = &ep.value;
                bool good = false;
                out.cd();
                TTree *out_ep = new TTree("epics", "PRad2 EPICS slow control (concatenated)");
                prad2::SetEpicsWriteBranches(out_ep, ep);
                out_ep->Branch("good", &good, "good/O");
                for (size_t fi : file_indices) {
                    std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                    TTree *t = f ? dynamic_cast<TTree *>(f->Get("epics")) : nullptr;
                    if (!t) continue;
                    t->SetBranchAddress("event_number_at_arrival", &ep.event_number_at_arrival);
                    const bool has_ticks_in = (t->GetBranch("ti_ticks_at_arrival") != nullptr);
                    if (has_ticks_in)
                        t->SetBranchAddress("ti_ticks_at_arrival", &ep.ti_ticks_at_arrival);
                    t->SetBranchAddress("unix_time",    &ep.unix_time);
                    t->SetBranchAddress("sync_counter", &ep.sync_counter);
                    t->SetBranchAddress("run_number",   &ep.run_number);
                    t->SetBranchAddress("channel", &cp);
                    t->SetBranchAddress("value",   &vp);
                    size_t seq = file_info[fi].epics_offset;
                    Long64_t n = t->GetEntries();
                    for (Long64_t i = 0; i < n; ++i, ++seq) {
                        ep.ti_ticks_at_arrival = 0;
                        t->GetEntry(i);
                        if (ep.ti_ticks_at_arrival <= 0 && ep.event_number_at_arrival >= 0) {
                            auto eit = evn_to_ticks.find(ep.event_number_at_arrival);
                            if (eit != evn_to_ticks.end()) ep.ti_ticks_at_arrival = eit->second;
                        }
                        good = (seq < epics_verdict.size()) ? epics_verdict[seq] : false;
                        if (good && restrict_side)
                            good = (seq < epics_side.size() && epics_side[seq] == s);
                        out.cd();
                        out_ep->Fill();
                    }
                }
                out.cd();
                out_ep->Write();
            }
            {
                prad2::RawRunInfo ri;
                std::string *sp = &ri.daq_config;
                out.cd();
                TTree *out_ri = new TTree("runinfo", "PRad2 control events / DAQ config (concatenated)");
                prad2::SetRunInfoWriteBranches(out_ri, ri);
                for (size_t fi : file_indices) {
                    std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                    TTree *t = f ? dynamic_cast<TTree *>(f->Get("runinfo")) : nullptr;
                    if (!t) continue;
                    prad2::SetRunInfoReadBranches(t, ri);
                    t->SetBranchAddress("daq_config", &sp);
                    Long64_t n = t->GetEntries();
                    for (Long64_t i = 0; i < n; ++i) {
                        ri.daq_config.clear();
                        t->GetEntry(i);
                        out.cd();
                        out_ri->Fill();
                    }
                }
                out.cd();
                out_ri->Write();
            }
            return 0;
        };

        struct WriteStats { int64_t n_in = 0, n_out = 0; };
        auto write_file_side = [&](size_t fi, int s, const std::string &outp) -> WriteStats {
            WriteStats stats;
            std::unique_ptr<TFile> out(TFile::Open(outp.c_str(), "RECREATE"));
            if (!out || out->IsZombie()) {
                std::cerr << "replay_filter: cannot create " << outp << "\n";
                return stats;
            }

            if (!is_recon) {
                prad2::RawEventData ev;
                prad2::RawReadStatus first_status;
                {
                    std::unique_ptr<TFile> f0(TFile::Open(file_info[fi].path.c_str(), "READ"));
                    TTree *t0 = f0 ? dynamic_cast<TTree *>(f0->Get("events")) : nullptr;
                    first_status = prad2::SetRawReadBranches(t0, ev);
                }
                out->cd();
                TTree *out_ev = new TTree("events", "PRad2 filtered replay (raw)");
                prad2::SetRawWriteBranches(out_ev, ev, first_status.has_peaks);
                std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                TTree *t = f ? dynamic_cast<TTree *>(f->Get("events")) : nullptr;
                if (t) {
                    prad2::SetRawReadBranches(t, ev);
                    std::vector<uint32_t> *p_ssp = &ev.ssp_raw;
                    if (t->GetBranch("ssp_raw")) t->SetBranchAddress("ssp_raw", &p_ssp);
                    std::vector<uint32_t> *p_vtp_roc = &ev.vtp_roc_tags, *p_vtp_nw = &ev.vtp_nwords, *p_vtp_w = &ev.vtp_words;
                    if (t->GetBranch("vtp_roc_tags")) t->SetBranchAddress("vtp_roc_tags", &p_vtp_roc);
                    if (t->GetBranch("vtp_nwords"))   t->SetBranchAddress("vtp_nwords",   &p_vtp_nw);
                    if (t->GetBranch("vtp_words"))    t->SetBranchAddress("vtp_words",    &p_vtp_w);
                    std::vector<uint32_t> *p_tdc_roc = &ev.tdc_roc_tags, *p_tdc_nw = &ev.tdc_nwords, *p_tdc_w = &ev.tdc_words;
                    if (t->GetBranch("tdc_roc_tags")) t->SetBranchAddress("tdc_roc_tags", &p_tdc_roc);
                    if (t->GetBranch("tdc_nwords"))   t->SetBranchAddress("tdc_nwords",   &p_tdc_nw);
                    if (t->GetBranch("tdc_words"))    t->SetBranchAddress("tdc_words",    &p_tdc_w);
                    Long64_t n = t->GetEntries();
                    if (!split_active) stats.n_in += n;
                    for (Long64_t i = 0; i < n; ++i) {
                        ev.ssp_raw.clear(); ev.vtp_roc_tags.clear(); ev.vtp_nwords.clear(); ev.vtp_words.clear();
                        ev.tdc_roc_tags.clear(); ev.tdc_nwords.clear(); ev.tdc_words.clear();
                        t->GetEntry(i);
                        if (split_active && in_span(s, ev.event_num)) ++stats.n_in;
                        if (is_kept(s, ev.event_num)) {
                            out->cd();
                            out_ev->Fill();
                            ++stats.n_out;
                        }
                    }
                }
                out->cd();
                out_ev->Write();
            } else {
                prad2::ReconEventData ev;
                out->cd();
                TTree *out_ev = new TTree("recon", "PRad2 filtered replay (recon)");
                prad2::SetReconWriteBranches(out_ev, ev);
                std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                TTree *t = f ? dynamic_cast<TTree *>(f->Get("recon")) : nullptr;
                if (t) {
                    prad2::SetReconReadBranches(t, ev);
                    std::vector<uint32_t> *p_ssp = &ev.ssp_raw;
                    if (t->GetBranch("ssp_raw")) t->SetBranchAddress("ssp_raw", &p_ssp);
                    std::vector<uint32_t> *p_vtp_roc = &ev.vtp_roc_tags, *p_vtp_nw = &ev.vtp_nwords, *p_vtp_w = &ev.vtp_words;
                    if (t->GetBranch("vtp_roc_tags")) t->SetBranchAddress("vtp_roc_tags", &p_vtp_roc);
                    if (t->GetBranch("vtp_nwords"))   t->SetBranchAddress("vtp_nwords",   &p_vtp_nw);
                    if (t->GetBranch("vtp_words"))    t->SetBranchAddress("vtp_words",    &p_vtp_w);
                    Long64_t n = t->GetEntries();
                    if (!split_active) stats.n_in += n;
                    for (Long64_t i = 0; i < n; ++i) {
                        ev.ssp_raw.clear(); ev.vtp_roc_tags.clear(); ev.vtp_nwords.clear(); ev.vtp_words.clear();
                        t->GetEntry(i);
                        if (split_active && in_span(s, ev.event_num)) ++stats.n_in;
                        if (is_kept(s, ev.event_num)) {
                            out->cd();
                            out_ev->Fill();
                            ++stats.n_out;
                        }
                    }
                }
                out->cd();
                out_ev->Write();
            }
            std::vector<size_t> one{fi};
            write_slow_trees(*out, one, s, split_active);
            out->Close();
            return stats;
        };

        std::vector<std::pair<size_t, int>> jobs;
        for (size_t fi = 0; fi < input_files.size(); ++fi) {
            if (!split_active) jobs.emplace_back(fi, 0);
            else {
                if (n_cp_side[0] > 0) jobs.emplace_back(fi, 0);
                if (n_cp_side[1] > 0) jobs.emplace_back(fi, 1);
            }
        }

        std::vector<std::string> output_files(jobs.size());
        std::vector<WriteStats> job_stats(jobs.size());
        std::atomic<size_t> next_job{0};
        std::mutex log_mtx;
        const size_t workers = std::min<size_t>(std::max(1, num_threads), jobs.size());
        std::vector<std::thread> threads;
        threads.reserve(workers);
        for (size_t w = 0; w < workers; ++w) {
            threads.emplace_back([&]() {
                while (true) {
                    size_t ji = next_job.fetch_add(1);
                    if (ji >= jobs.size()) break;
                    const auto [fi, s] = jobs[ji];
                    const std::string label = split_active
                        ? (s == 0 ? cuts.split.label_full : cuts.split.label_empty)
                        : std::string();
                    std::string outp = filtered_output_path(output_path, file_info[fi].path, label);
                    job_stats[ji] = write_file_side(fi, s, outp);
                    output_files[ji] = outp;
                    std::lock_guard<std::mutex> lk(log_mtx);
                    std::cerr << "replay_filter: output ROOT     " << outp << "\n";
                }
            });
        }
        for (auto &t : threads) t.join();

        const std::string slow_out = run_epics_path(output_path, run_number);
        {
            std::unique_ptr<TFile> out(TFile::Open(slow_out.c_str(), "RECREATE"));
            if (!out || out->IsZombie()) {
                std::cerr << "replay_filter: cannot create " << slow_out << "\n";
                return 1;
            }
            write_slow_trees(*out, all_file_indices, 0, false);
            out->Close();
        }
        std::cerr << "replay_filter: run slow ROOT  " << slow_out << "\n";

        int64_t n_in_total = 0, n_pass_phys = 0;
        for (const auto &st : job_stats) {
            n_in_total += st.n_in;
            n_pass_phys += st.n_out;
        }

        auto to_json_or_null = [](bool valid, double v) -> json {
            return valid ? json(v) : json(nullptr);
        };
        auto stats_for = [&](const ChannelCut &c) -> json {
            json s = {
                {"abs_min", c.abs.has_min ? json(c.abs.min_val) : json(nullptr)},
                {"abs_max", c.abs.has_max ? json(c.abs.max_val) : json(nullptr)},
            };
            if (c.has_rel_rms) {
                s["rel_rms"]       = c.rel_rms_n;
                s["robust_center"] = to_json_or_null(c.stats_valid, c.robust_center);
                s["robust_sigma"]  = to_json_or_null(c.stats_valid, c.robust_sigma);
                s["mad"]           = to_json_or_null(c.stats_valid, c.mad);
                s["n_used"]        = c.n_used;
                s["n_clipped"]     = c.n_clipped;
                if (!c.gated_by.empty()) s["gated_by"] = c.gated_by;
            }
            return s;
        };

        json report;
        report["run_number"] = run_number;
        report["input_files"] = input_files;
        report["output_files"] = output_files;
        report["slow_output_file"] = slow_out;
        report["cuts"] = cuts.raw;
        report["robust_method"] = "mad";

        json stats = json::object();
        if (cuts.livetime.enabled) {
            json ls = stats_for(cuts.livetime.cut);
            ls["source"] = cuts.livetime.source;
            ls["channel"] = cuts.livetime.channel;
            stats["livetime"] = std::move(ls);
        }
        for (const auto &kv : cuts.epics)
            stats["epics:" + kv.first] = stats_for(kv.second);
        report["stats"] = std::move(stats);

        int n_pass_cp = 0, n_fail_cp = 0;
        for (const auto &cp : timeline)
            (cp.overall_pass ? n_pass_cp : n_fail_cp)++;
        std::map<std::string, std::pair<int, int>> per_channel;
        for (const auto &p : report_points) {
            auto &c = per_channel[p.channel];
            if (p.pass) ++c.first; else ++c.second;
        }
        json per_channel_json = json::object();
        for (const auto &kv : per_channel) {
            int pass = kv.second.first, fail = kv.second.second, tot = pass + fail;
            per_channel_json[kv.first] = {
                {"n_pass", pass},
                {"n_fail", fail},
                {"pass_rate", tot > 0 ? double(pass) / double(tot) : 0.0},
            };
        }
        json keep_intervals = json::array();
        for (int s = 0; s <= (split_active ? 1 : 0); ++s)
            for (const auto &p : keep[s])
                keep_intervals.push_back({p.first, p.second});

        report["summary"] = {
            {"n_slow_events", (int)timeline.size()},
            {"n_slow_pass", n_pass_cp},
            {"n_slow_reject", n_fail_cp},
            {"slow_pass_rate", timeline.empty() ? 0.0 : double(n_pass_cp) / double(timeline.size())},
            {"n_physics_in", n_in_total},
            {"n_physics_pass", n_pass_phys},
            {"n_physics_reject", n_in_total - n_pass_phys},
            {"physics_pass_rate", n_in_total > 0 ? double(n_pass_phys) / double(n_in_total) : 0.0},
            {"n_keep_intervals", (int)keep_intervals.size()},
            {"per_channel", per_channel_json},
        };
        report["keep_intervals"] = std::move(keep_intervals);
        if (split_active) {
            report["split"] = {
                {"enabled", true},
                {"channel", cuts.split.channel},
                {"full_threshold", cuts.split.full_thresh},
                {"empty_threshold", cuts.split.empty_thresh},
                {"guard_checkpoints", cuts.split.guard_checkpoints},
                {"labels", {cuts.split.label_full, cuts.split.label_empty}},
                {"n_checkpoints_full", n_cp_side[0]},
                {"n_checkpoints_empty", n_cp_side[1]},
                {"n_checkpoints_dropped", n_cp_guard},
                {"n_state_transitions", n_state_transitions},
                {"transitions", transitions},
                {"pure_run", n_state_transitions == 0},
            };
        }
        if (cuts.charge.enabled) {
            report["live_charge"] = {
                {"value_nC", live_charge[0] + live_charge[1]},
                {"unit", "nC"},
                {"beam_current_channel", cuts.charge.beam_current_channel},
                {"beam_current_unit", "nA"},
                {"live_seconds", live_charge_secs[0] + live_charge_secs[1]},
                {"real_seconds", real_secs[0] + real_secs[1]},
                {"ungated_value_nC", ungated_live_charge[0] + ungated_live_charge[1]},
                {"ungated_live_seconds", ungated_live_charge_secs[0] + ungated_live_charge_secs[1]},
                {"ungated_real_seconds", ungated_real_secs[0] + ungated_real_secs[1]},
                {"n_pairs_integrated", n_charge_pairs[0] + n_charge_pairs[1]},
                {"n_pairs_skipped", n_charge_skipped[0] + n_charge_skipped[1]},
                {"n_ungated_pairs_integrated", n_ungated_charge_pairs[0] + n_ungated_charge_pairs[1]},
                {"n_ungated_pairs_skipped", n_ungated_charge_skipped[0] + n_ungated_charge_skipped[1]},
            };
        }
        json pts = json::array();
        pts.get_ptr<json::array_t *>()->reserve(report_points.size());
        for (const auto &p : report_points) {
            pts.push_back({
                {"channel", p.channel},
                {"status", p.pass ? "pass" : "fail"},
                {"associated_evn", p.event_number},
                {"associated_timestamp", p.has_assoc_t ? json(p.assoc_t_rel) : json(nullptr)},
                {"unix_time", p.has_unix_time ? json(p.unix_time) : json(nullptr)},
                {"value", std::isnan(p.value) ? json(nullptr) : json(p.value)},
            });
        }
        report["points"] = std::move(pts);

        std::ofstream of(report_path);
        if (!of) {
            std::cerr << "replay_filter: cannot write " << report_path << "\n";
            return 1;
        }
        of << report.dump(2) << "\n";
        std::cerr << "replay_filter: report written to " << report_path << "\n";
        std::cerr << "  physics      : in=" << n_in_total
                  << "  pass=" << n_pass_phys
                  << "  reject=" << (n_in_total - n_pass_phys)
                  << "  rate=" << fmt_pct(n_in_total > 0 ? double(n_pass_phys) / n_in_total : 0)
                  << "\n";
        return 0;
    }

    auto write_side = [&](int s, const std::string &outp,
                          const std::string &repp) -> int {
    std::unique_ptr<TFile> out(TFile::Open(outp.c_str(), "RECREATE"));
    if (!out || out->IsZombie()) {
        std::cerr << "replay_filter: cannot create " << outp << "\n";
        return 1;
    }

    int64_t n_in = 0, n_out = 0;

    // Filter the events/recon tree.  We use the existing
    // SetRaw{Read,Write}Branches helpers so the output schema matches.
    if (!is_recon) {
        prad2::RawEventData    ev;
        prad2::RawReadStatus   first_status;
        {
            std::unique_ptr<TFile> f0(TFile::Open(input_files.front().c_str(), "READ"));
            TTree *t0 = dynamic_cast<TTree *>(f0->Get("events"));
            first_status = prad2::SetRawReadBranches(t0, ev);
        }
        out->cd();
        TTree *out_ev = new TTree("events", "PRad2 filtered replay (raw)");
        prad2::SetRawWriteBranches(out_ev, ev, first_status.has_peaks);

        for (const auto &path : input_files) {
            std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
            TTree *t = dynamic_cast<TTree *>(f->Get("events"));
            if (!t) continue;
            prad2::SetRawReadBranches(t, ev);
            std::vector<uint32_t> *p_ssp = &ev.ssp_raw;
            if (t->GetBranch("ssp_raw")) t->SetBranchAddress("ssp_raw", &p_ssp);
            // Vector branches need pointer-to-pointer rebinding — same
            // dance as ssp_raw above.  Older replays without these
            // branches just leave the pointers untouched.
            std::vector<uint32_t> *p_vtp_roc = &ev.vtp_roc_tags;
            std::vector<uint32_t> *p_vtp_nw  = &ev.vtp_nwords;
            std::vector<uint32_t> *p_vtp_w   = &ev.vtp_words;
            if (t->GetBranch("vtp_roc_tags")) t->SetBranchAddress("vtp_roc_tags", &p_vtp_roc);
            if (t->GetBranch("vtp_nwords"))   t->SetBranchAddress("vtp_nwords",   &p_vtp_nw);
            if (t->GetBranch("vtp_words"))    t->SetBranchAddress("vtp_words",    &p_vtp_w);
            std::vector<uint32_t> *p_tdc_roc = &ev.tdc_roc_tags;
            std::vector<uint32_t> *p_tdc_nw  = &ev.tdc_nwords;
            std::vector<uint32_t> *p_tdc_w   = &ev.tdc_words;
            if (t->GetBranch("tdc_roc_tags")) t->SetBranchAddress("tdc_roc_tags", &p_tdc_roc);
            if (t->GetBranch("tdc_nwords"))   t->SetBranchAddress("tdc_nwords",   &p_tdc_nw);
            if (t->GetBranch("tdc_words"))    t->SetBranchAddress("tdc_words",    &p_tdc_w);
            Long64_t n = t->GetEntries();
            if (!split_active) n_in += n;       // split off: every event counts
            for (Long64_t i = 0; i < n; ++i) {
                ev.ssp_raw.clear();
                ev.vtp_roc_tags.clear();
                ev.vtp_nwords.clear();
                ev.vtp_words.clear();
                ev.tdc_roc_tags.clear();
                ev.tdc_nwords.clear();
                ev.tdc_words.clear();
                t->GetEntry(i);
                if (split_active && in_span(s, ev.event_num)) ++n_in;
                if (is_kept(s, ev.event_num)) {
                    out->cd();
                    out_ev->Fill();
                    ++n_out;
                }
            }
        }
        out->cd();
        out_ev->Write();
    } else {
        prad2::ReconEventData ev;
        out->cd();
        TTree *out_ev = new TTree("recon", "PRad2 filtered replay (recon)");
        prad2::SetReconWriteBranches(out_ev, ev);

        for (const auto &path : input_files) {
            std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
            TTree *t = dynamic_cast<TTree *>(f->Get("recon"));
            if (!t) continue;
            prad2::SetReconReadBranches(t, ev);
            std::vector<uint32_t> *p_ssp = &ev.ssp_raw;
            if (t->GetBranch("ssp_raw")) t->SetBranchAddress("ssp_raw", &p_ssp);
            // vtp_* vector branches — same pointer-to-pointer dance as
            // ssp_raw; recon files older than 2026-06 don't have them
            // and the cleared vectors just pass through empty.
            std::vector<uint32_t> *p_vtp_roc = &ev.vtp_roc_tags;
            std::vector<uint32_t> *p_vtp_nw  = &ev.vtp_nwords;
            std::vector<uint32_t> *p_vtp_w   = &ev.vtp_words;
            if (t->GetBranch("vtp_roc_tags")) t->SetBranchAddress("vtp_roc_tags", &p_vtp_roc);
            if (t->GetBranch("vtp_nwords"))   t->SetBranchAddress("vtp_nwords",   &p_vtp_nw);
            if (t->GetBranch("vtp_words"))    t->SetBranchAddress("vtp_words",    &p_vtp_w);
            Long64_t n = t->GetEntries();
            if (!split_active) n_in += n;       // split off: every event counts
            for (Long64_t i = 0; i < n; ++i) {
                ev.ssp_raw.clear();
                ev.vtp_roc_tags.clear();
                ev.vtp_nwords.clear();
                ev.vtp_words.clear();
                t->GetEntry(i);
                if (split_active && in_span(s, ev.event_num)) ++n_in;
                if (is_kept(s, ev.event_num)) {
                    out->cd();
                    out_ev->Fill();
                    ++n_out;
                }
            }
        }
        out->cd();
        out_ev->Write();
    }

    // Concatenate scalers tree from every input.  Adds a `good` boolean
    // (per-checkpoint overall verdict from phase 3) so downstream tools can
    // colour the run's livetime trace by pass/fail without recomputing.
    {
        prad2::RawScalerData sc;
        bool                 good = false;
        out->cd();
        TTree *out_sc = new TTree("scalers", "PRad2 DSC2 scaler readouts (concatenated)");
        prad2::SetScalerWriteBranches(out_sc, sc);
        out_sc->Branch("good", &good, "good/O");

        size_t seq = 0;
        for (const auto &path : input_files) {
            std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
            TTree *t = dynamic_cast<TTree *>(f->Get("scalers"));
            if (!t) continue;
            prad2::SetScalerReadBranches(t, sc);
            Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                good = (seq < scaler_verdict.size()) ? scaler_verdict[seq] : false;
                if (good && split_active)
                    good = (seq < scaler_side.size() && scaler_side[seq] == s);
                ++seq;
                out->cd();
                out_sc->Fill();
            }
        }
        out->cd();
        out_sc->Write();
    }

    // Concatenate epics tree from every input, tagged the same way.
    // We also resolve ti_ticks_at_arrival here: bind it from the input
    // when present, otherwise fill it from the events-tree lookup so the
    // output is always self-contained — downstream consumers (live-charge
    // recomputation, etc.) read the row's tick directly without needing
    // to know whether the upstream replay carried the new branch.
    {
        prad2::RawEpicsData ep;
        std::vector<std::string> *cp = &ep.channel;
        std::vector<double>      *vp = &ep.value;
        bool good = false;
        out->cd();
        TTree *out_ep = new TTree("epics", "PRad2 EPICS slow control (concatenated)");
        prad2::SetEpicsWriteBranches(out_ep, ep);
        out_ep->Branch("good", &good, "good/O");

        size_t seq = 0;
        for (const auto &path : input_files) {
            std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
            TTree *t = dynamic_cast<TTree *>(f->Get("epics"));
            if (!t) continue;
            t->SetBranchAddress("event_number_at_arrival", &ep.event_number_at_arrival);
            const bool has_ticks_in =
                (t->GetBranch("ti_ticks_at_arrival") != nullptr);
            if (has_ticks_in)
                t->SetBranchAddress("ti_ticks_at_arrival", &ep.ti_ticks_at_arrival);
            t->SetBranchAddress("unix_time",    &ep.unix_time);
            t->SetBranchAddress("sync_counter", &ep.sync_counter);
            t->SetBranchAddress("run_number",   &ep.run_number);
            t->SetBranchAddress("channel", &cp);
            t->SetBranchAddress("value",   &vp);
            Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                ep.ti_ticks_at_arrival = 0;
                t->GetEntry(i);
                if (ep.ti_ticks_at_arrival <= 0
                    && ep.event_number_at_arrival >= 0) {
                    auto eit = evn_to_ticks.find(ep.event_number_at_arrival);
                    if (eit != evn_to_ticks.end())
                        ep.ti_ticks_at_arrival = eit->second;
                }
                good = (seq < epics_verdict.size()) ? epics_verdict[seq] : false;
                if (good && split_active)
                    good = (seq < epics_side.size() && epics_side[seq] == s);
                ++seq;
                out->cd();
                out_ep->Fill();
            }
        }
        out->cd();
        out_ep->Write();
    }

    // Pass through the runinfo tree (1 row per CODA control event,
    // including the long DAQ-config text on PRESTART).  No filtering
    // applied — the whole point of this tree is run-scoped metadata.
    {
        prad2::RawRunInfo ri;
        std::string      *sp = &ri.daq_config;
        out->cd();
        TTree *out_ri = new TTree("runinfo",
                                  "PRad2 control events / DAQ config (concatenated)");
        prad2::SetRunInfoWriteBranches(out_ri, ri);
        for (const auto &path : input_files) {
            std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
            TTree *t = dynamic_cast<TTree *>(f->Get("runinfo"));
            if (!t) continue;
            prad2::SetRunInfoReadBranches(t, ri);
            t->SetBranchAddress("daq_config", &sp);
            Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                ri.daq_config.clear();
                t->GetEntry(i);
                out->cd();
                out_ri->Fill();
            }
        }
        out->cd();
        out_ri->Write();
    }

    out->Close();

    // ---------- Phase 6: write the JSON report ----------
    auto to_json_or_null = [](bool valid, double v) -> json {
        return valid ? json(v) : json(nullptr);
    };
    auto stats_for = [&](const ChannelCut &c) -> json {
        json s = {
            {"abs_min", c.abs.has_min ? json(c.abs.min_val) : json(nullptr)},
            {"abs_max", c.abs.has_max ? json(c.abs.max_val) : json(nullptr)},
        };
        if (c.has_rel_rms) {
            s["rel_rms"]       = c.rel_rms_n;
            s["robust_center"] = to_json_or_null(c.stats_valid, c.robust_center);
            s["robust_sigma"]  = to_json_or_null(c.stats_valid, c.robust_sigma);
            s["mad"]           = to_json_or_null(c.stats_valid, c.mad);
            // n_used is the count *after* gating (if any) — useful for
            // sanity-checking that the gating restriction left enough data
            // to compute meaningful stats.
            s["n_used"]        = c.n_used;
            s["n_clipped"]     = c.n_clipped;
            if (!c.gated_by.empty()) {
                // Always emit as array — even single-gate cases — so
                // downstream tools can iterate without checking type.
                s["gated_by"] = c.gated_by;
            }
        }
        return s;
    };

    json report;
    report["run_number"]   = run_number;
    report["input_files"]  = input_files;
    report["output_file"]  = outp;
    report["cuts"]         = cuts.raw;
    report["robust_method"] = "mad";   // 1.4826 * MAD as σ̂

    json stats = json::object();
    if (cuts.livetime.enabled) {
        json ls = stats_for(cuts.livetime.cut);
        ls["source"]  = cuts.livetime.source;
        ls["channel"] = cuts.livetime.channel;
        stats["livetime"] = std::move(ls);
    }
    for (const auto &kv : cuts.epics)
        stats["epics:" + kv.first] = stats_for(kv.second);
    report["stats"] = std::move(stats);

    // Slow-event counts — restricted to this side's checkpoints when splitting,
    // so each side report describes only its own slice of the run.
    int n_pass_cp = 0, n_fail_cp = 0, n_slow_side = 0;
    for (size_t k = 0; k < timeline.size(); ++k) {
        if (split_active && side[k] != s) continue;
        ++n_slow_side;
        (timeline[k].overall_pass ? n_pass_cp : n_fail_cp)++;
    }
    json keep_intervals = json::array();
    for (const auto &p : keep[s]) keep_intervals.push_back({p.first, p.second});

    // Per-channel breakdown — number of slow-event checkpoints where this
    // channel's cut accepted vs rejected the value.  Helps the user see
    // immediately which cut is doing the rejecting (e.g. "beam current
    // killed 80% of points, livetime barely matters").
    std::map<std::string, std::pair<int, int>> per_channel;   // ch → {pass, fail}
    for (const auto &p : report_points) {
        if (split_active && p.side != s) continue;    // this side's points only
        auto &c = per_channel[p.channel];
        if (p.pass) ++c.first; else ++c.second;
    }
    json per_channel_json = json::object();
    for (const auto &kv : per_channel) {
        int pass = kv.second.first, fail = kv.second.second;
        int tot  = pass + fail;
        per_channel_json[kv.first] = {
            {"n_pass",    pass},
            {"n_fail",    fail},
            {"pass_rate", tot > 0 ? double(pass) / double(tot) : 0.0},
        };
    }

    const int64_t n_in_total  = n_in;
    const int64_t n_pass_phys = n_out;
    const int64_t n_rej_phys  = n_in_total - n_pass_phys;
    const int     n_slow      = n_slow_side;

    report["summary"] = {
        // Slow-event checkpoint counts.
        {"n_slow_events",      n_slow},
        {"n_slow_pass",        n_pass_cp},
        {"n_slow_reject",      n_fail_cp},
        {"slow_pass_rate",     n_slow > 0 ? double(n_pass_cp) / double(n_slow) : 0.0},
        // Physics-event counts.
        {"n_physics_in",       n_in_total},
        {"n_physics_pass",     n_pass_phys},
        {"n_physics_reject",   n_rej_phys},
        {"physics_pass_rate",  n_in_total > 0
                                ? double(n_pass_phys) / double(n_in_total) : 0.0},
        // Keep-interval count (each is a (lo, hi] range of accepted events).
        {"n_keep_intervals",   (int)keep[s].size()},
        // Per-cut breakdown — which channel rejected how often.
        {"per_channel",        per_channel_json},
    };
    report["keep_intervals"] = std::move(keep_intervals);

    // Split metadata — present only when run-splitting is active.  Each side
    // file is self-describing: it records the level thresholds, how many
    // checkpoints fell in each state vs the dropped ramp, and every
    // full<->empty transition seen, so neither side needs the other's report.
    if (split_active) {
        report["split"] = {
            {"enabled",               true},
            {"channel",               cuts.split.channel},
            {"side",                  s},
            {"state",                 s == 0 ? "full" : "empty"},
            {"label",                 s == 0 ? cuts.split.label_full
                                             : cuts.split.label_empty},
            {"full_threshold",        cuts.split.full_thresh},
            {"empty_threshold",       cuts.split.empty_thresh},
            {"guard_checkpoints",     cuts.split.guard_checkpoints},
            {"n_checkpoints_this",    n_cp_side[s]},
            {"n_checkpoints_full",    n_cp_side[0]},
            {"n_checkpoints_empty",   n_cp_side[1]},
            {"n_checkpoints_dropped", n_cp_guard},
            {"n_state_transitions",   n_state_transitions},
            {"transitions",           transitions},
            {"pure_run",              n_state_transitions == 0},
        };
    }

    // Live-charge integration over kept intervals.  Units: assume the
    // configured EPICS beam-current channel publishes in nA (true for
    // hallb_IPM2C21A_CUR and the other Hall B IPM scalers), so
    // value = Σ live_fraction · Δt · I  ⇒  nA · s = nC.  Also emit the
    // accumulated live time so the average current is recoverable.
    //
    // value_nC / live_seconds / real_seconds are the gated (post-cut)
    // numbers — only adjacent pairs where every cut passed on both
    // endpoints contribute, matching the keep-interval physics events
    // written to the output ROOT file.  ungated_* are the same
    // quantities accumulated over every valid-data pair (no cut
    // applied), reported so users can see how much charge the cuts
    // threw away.  When the cuts reject nothing, ungated_* equals the
    // gated values.
    if (cuts.charge.enabled) {
        report["live_charge"] = {
            {"value_nC",                   live_charge[s]},
            {"unit",                       "nC"},
            {"beam_current_channel",       cuts.charge.beam_current_channel},
            {"beam_current_unit",          "nA"},
            {"live_seconds",               live_charge_secs[s]},
            {"real_seconds",               real_secs[s]},
            {"ungated_value_nC",           ungated_live_charge[s]},
            {"ungated_live_seconds",       ungated_live_charge_secs[s]},
            {"ungated_real_seconds",       ungated_real_secs[s]},
            {"n_pairs_integrated",         n_charge_pairs[s]},
            {"n_pairs_skipped",            n_charge_skipped[s]},
            {"n_ungated_pairs_integrated", n_ungated_charge_pairs[s]},
            {"n_ungated_pairs_skipped",    n_ungated_charge_skipped[s]},
        };
    }

    json pts = json::array();
    pts.get_ptr<json::array_t *>()->reserve(report_points.size());
    for (const auto &p : report_points) {
        json e = {
            {"channel",              p.channel},
            {"status",               p.pass ? "pass" : "fail"},
            {"associated_evn",       p.event_number},
            // TI ticks of the associated physics event, in seconds since the
            // run's earliest event.  Null when no physics event matches
            // (associated_evn=-1, or input has no events tree).
            {"associated_timestamp", p.has_assoc_t ? json(p.assoc_t_rel) : json(nullptr)},
            // Native EPICS unix_time on EPICS rows; null on scaler rows.
            {"unix_time",            p.has_unix_time ? json(p.unix_time) : json(nullptr)},
            {"value",                std::isnan(p.value) ? json(nullptr) : json(p.value)},
        };
        pts.push_back(std::move(e));
    }
    report["points"] = std::move(pts);

    std::ofstream of(repp);
    if (!of) {
        std::cerr << "replay_filter: cannot write " << repp << "\n";
        return 1;
    }
    of << report.dump(2) << "\n";

    const std::string tag = split_active
        ? "[" + (s == 0 ? cuts.split.label_full : cuts.split.label_empty) + "] " : "";
    std::cerr << "replay_filter: " << tag << "report written to " << repp << "\n";
    std::cerr << "replay_filter: " << tag << "output ROOT     " << outp << "\n";
    std::cerr << "  slow events  : " << n_slow
              << "  pass=" << n_pass_cp << "  reject=" << n_fail_cp
              << "  rate=" << fmt_pct(n_slow > 0 ? double(n_pass_cp) / n_slow : 0)
              << "\n";
    std::cerr << "  keep intervals: " << keep[s].size() << "\n";
    std::cerr << "  physics      : in=" << n_in_total
              << "  pass="  << n_pass_phys
              << "  reject=" << n_rej_phys
              << "  rate=" << fmt_pct(n_in_total > 0
                                       ? double(n_pass_phys) / n_in_total : 0)
              << "\n";
    if (cuts.charge.enabled)
        std::cerr << "  live charge  : " << std::fixed << std::setprecision(3)
                  << live_charge[s] << " nC over " << live_charge_secs[s]
                  << " s live" << std::defaultfloat << "\n";
    if (!per_channel.empty()) {
        std::cerr << "  per-channel reject:\n";
        for (const auto &kv : per_channel) {
            std::cerr << "    " << kv.first << ": "
                      << kv.second.second << " / "
                      << (kv.second.first + kv.second.second)
                      << " ("
                      << fmt_pct(double(kv.second.second) /
                                 std::max(1, kv.second.first + kv.second.second))
                      << ")\n";
        }
    }
    return 0;
    };   // end write_side

    // ---------- Dispatch: single output, or one file per present target state ----------
    // with_label suffixes the file/report stems with a side label.  Handles the
    // compound ".report.json" suffix so a report becomes
    // "<stem>_<label>.report.json" rather than "<stem>.report_<label>.json".
    auto with_label = [](const std::string &path, const std::string &label) {
        static const std::vector<std::string> compound = {".report.json"};
        for (const auto &suf : compound) {
            if (path.size() > suf.size() &&
                path.compare(path.size() - suf.size(), suf.size(), suf) == 0)
                return path.substr(0, path.size() - suf.size())
                       + "_" + label + suf;
        }
        auto dot = path.rfind('.');
        return (dot == std::string::npos)
            ? path + "_" + label
            : path.substr(0, dot) + "_" + label + path.substr(dot);
    };

    if (!split_active)        // split off or degraded → the single default output
        return write_side(0, output_path, report_path);

    // Level split: write one labelled file per target state that actually
    // occurred (n_cp_side[s] > 0).  A pure run thus yields a single labelled
    // file (e.g. <stem>_full.root); a run with a transition yields both.
    const char *names[2]       = {"full", "empty"};
    const std::string labels[2] = {cuts.split.label_full, cuts.split.label_empty};
    std::cerr << "replay_filter: split on '" << cuts.split.channel
              << "' by level (full >= " << cuts.split.full_thresh
              << ", empty <= " << cuts.split.empty_thresh << "): "
              << n_cp_side[0] << " full + " << n_cp_side[1] << " empty checkpoints, "
              << n_cp_guard << " dropped, " << n_state_transitions
              << " transition(s)\n";
    int rc = 0;
    for (int s = 0; s <= 1 && rc == 0; ++s) {
        if (n_cp_side[s] == 0) {
            std::cerr << "replay_filter: no " << names[s]
                      << "-target checkpoints — skipping " << names[s] << " file\n";
            continue;
        }
        rc = write_side(s, with_label(output_path, labels[s]),
                           with_label(report_path, labels[s]));
    }
    return rc;
}

void usage(const char *prog)
{
    std::cerr <<
        "Usage: " << prog << " <input.root> [more.root ...]\n"
        "       -o <output.root|output_dir>  -c <cuts.json> [-j <report.json>]\n"
        "       [-r <run_num>] [-t threads]\n"
        "\n"
        "Filters replayed ROOT files by slow-control cuts (DSC2 livetime\n"
        "+ EPICS).  Writes a single ROOT file with the kept physics events\n"
        "and the full scaler/epics streams concatenated, plus a JSON report\n"
        "with per-(cut, slow-event) pass/fail status for chart plotting.\n"
        "With multiple inputs, -o is an output directory; one filtered ROOT\n"
        "is written per input by inserting _filter before the final .root,\n"
        "plus prad_<run>_epics.root and one run-level JSON report.\n"
        "\n"
        "With a \"split\" block events are instead classified by a PV level\n"
        "(e.g. target cell pressure) and the full/empty subsets written to\n"
        "separate files <stem>_<full>.root / <stem>_<empty>.root.  Single-input\n"
        "mode writes side reports; multi-input mode writes one run-level report.\n"
        "PV>=full and PV<=empty\n"
        "select the two states; the in-between ramp lands in neither.  Only\n"
        "the states that occur are written, so a pure run yields one file.\n"
        "\n"
        "Cut JSON example:\n"
        "  {\n"
        "    \"livetime\": { \"source\": \"ref\", \"abs\": { \"min\": 90 } },\n"
        "    \"epics\": {\n"
        "      \"hallb_IPM2C21A_CUR\":  { \"abs\": { \"min\": 3 } },\n"
        "      \"hallb_IPM2C21A_XPOS\": { \"rel_rms\": 3 },\n"
        "      \"hallb_IPM2C21A_YPOS\": { \"rel_rms\": 3 }\n"
        "    },\n"
        "    \"split\": {\n"
        "      \"channel\": \"TGT:PRad:Cell_P\", \"full\": 500, \"empty\": 5,\n"
        "      \"guard_checkpoints\": 0, \"labels\": [\"full\", \"empty\"]\n"
        "    }\n"
        "  }\n"
        "  (split needs full > empty; guard_checkpoints drops an extra +/- N\n"
        "   points at each state edge, beyond the dead zone.)\n";
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    std::vector<std::string> inputs;
    std::string output, cuts_path, report_path;
    int run_override = -1;
    int num_threads = 1;

    ROOT::EnableThreadSafety();

    int opt;
    while ((opt = getopt(argc, argv, "o:c:j:r:t:h")) != -1) {
        switch (opt) {
        case 'o': output       = optarg;           break;
        case 'c': cuts_path    = optarg;           break;
        case 'j': report_path  = optarg;           break;
        case 'r': run_override = std::atoi(optarg); break;
        case 't': num_threads  = std::atoi(optarg); break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    for (int i = optind; i < argc; ++i) inputs.push_back(argv[i]);

    if (inputs.empty() || output.empty() || cuts_path.empty()) {
        usage(argv[0]);
        return 1;
    }
    num_threads = std::max(1, num_threads);
    if (report_path.empty()) {
        if (inputs.size() > 1) {
            int rn = run_override >= 0 ? run_override : analysis::get_run_int(inputs.front());
            report_path = run_report_path(output, rn);
        } else {
            auto dot = output.rfind('.');
            report_path = (dot == std::string::npos)
                ? output + ".report.json"
                : output.substr(0, dot) + ".report.json";
        }
    }
    return run(inputs, output, cuts_path, report_path, run_override, num_threads);
}
