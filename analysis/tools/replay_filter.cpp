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
// Optional "split" block — partition the run at a slow-control PV transition
// (target cell pressure stepping empty<->full mid-run) into TWO outputs, the
// events before the crossing (side A) and after it (side B), each still
// subject to every other cut and each with its own report + charge.  A guard
// band brackets the crossing so events taken while the PV ramps land in
// neither side:
//   "split": {
//     "channel":  "TGT:PRad:Cell_P",   // EPICS channel to watch
//     "threshold": 1.0,                 // single crossing, OR …
//     "level_low": 0.2, "level_high": 2.0,   // … hysteresis band
//     "guard": { "mode": "checkpoints", "checkpoints": 2 },  // OR
//     "guard": { "mode": "stability", "epsilon": 0.1, "consecutive": 3 },
//     "labels": ["full", "empty"]       // suffixes for the two files
//   }
//
// Output ROOT file(s):
//   * events / recon — same schema as input, only kept events
//   * scalers / epics — concatenated from every input plus an extra
//     `good` boolean branch per row reflecting that checkpoint's
//     overall verdict (all cuts pass); with split on, `good` also requires
//     the row to belong to that file's side, so live_charge on a side file
//     reproduces that side's post-cut charge directly.
// With split on, two files <stem>_<labelA>.root / <stem>_<labelB>.root and
// two matching reports are written instead of one.
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
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <getopt.h>

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

// Run-splitting on a slow-control PV transition (e.g. target cell pressure
// stepping between full and empty within one DAQ run).  When enabled the
// filter detects the first sustained crossing of `channel` past `threshold`
// (or, with `level_low`/`level_high` set, a hysteresis band) and writes TWO
// output files instead of one — the events before the transition (side A)
// and after it (side B), each still subject to every other configured cut
// and each with its own charge integral.
//
// A guard band brackets the transition so events taken while the PV is still
// ramping land in neither output.  Two guard modes:
//   * "checkpoints" — drop `guard_checkpoints` slow-control points on each
//     side of the crossing (cheap, predictable; checkpoints are ~regularly
//     spaced so this is effectively a time guard).  This is the example.
//   * "stability"   — keep extending the guard outward from the crossing
//     until the PV has settled within `guard_epsilon` of each side's level
//     for `guard_consecutive` consecutive readings.  Most faithful to the
//     physical ramp when the settle time is not constant.
struct SplitCfg {
    bool        enabled  = false;
    std::string channel;                       // EPICS channel to watch

    // Detection.  If level_low/level_high are both set, use them as a
    // hysteresis band (robust to noise); otherwise a single `threshold`
    // crossing (with the guard supplying the deadband).
    bool   has_threshold = false;
    double threshold     = 0.0;
    bool   has_levels    = false;
    double level_low     = 0.0;
    double level_high    = 0.0;

    // Guard band around the crossing.
    enum class Guard { Checkpoints, Stability };
    Guard  guard            = Guard::Checkpoints;
    int    guard_checkpoints = 2;              // ± N points (the example)
    double guard_epsilon     = 0.0;            // "settled within" band
    int    guard_consecutive = 3;              // K consecutive settled points

    // Output labels appended to the stems of the two files / reports.
    std::string label_a = "A";
    std::string label_b = "B";
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

    // Run-splitting on a PV transition (opt-in: needs a `channel` plus either
    // a `threshold` or a `level_low`/`level_high` hysteresis band).
    if (j.contains("split") && j["split"].is_object()) {
        const auto &sj = j["split"];
        auto &s = cfg.split;
        if (sj.contains("channel") && sj["channel"].is_string())
            s.channel = sj["channel"].get<std::string>();
        if (sj.contains("threshold") && sj["threshold"].is_number()) {
            s.has_threshold = true;
            s.threshold     = sj["threshold"].get<double>();
        }
        if (sj.contains("level_low")  && sj["level_low"].is_number() &&
            sj.contains("level_high") && sj["level_high"].is_number()) {
            s.has_levels  = true;
            s.level_low   = sj["level_low"].get<double>();
            s.level_high  = sj["level_high"].get<double>();
            if (s.level_low > s.level_high) std::swap(s.level_low, s.level_high);
            // A single threshold midway between the levels is the fallback
            // crossing test when hysteresis can't latch (e.g. PV starts
            // inside the band).
            if (!s.has_threshold) {
                s.has_threshold = true;
                s.threshold     = 0.5 * (s.level_low + s.level_high);
            }
        }
        if (sj.contains("guard") && sj["guard"].is_object()) {
            const auto &gj = sj["guard"];
            std::string mode = gj.value("mode", "checkpoints");
            s.guard = (mode == "stability") ? SplitCfg::Guard::Stability
                                            : SplitCfg::Guard::Checkpoints;
            if (gj.contains("checkpoints") && gj["checkpoints"].is_number())
                s.guard_checkpoints = gj["checkpoints"].get<int>();
            if (gj.contains("epsilon") && gj["epsilon"].is_number())
                s.guard_epsilon = gj["epsilon"].get<double>();
            if (gj.contains("consecutive") && gj["consecutive"].is_number())
                s.guard_consecutive = gj["consecutive"].get<int>();
        }
        if (sj.contains("labels") && sj["labels"].is_array()
            && sj["labels"].size() >= 2) {
            s.label_a = sj["labels"][0].get<std::string>();
            s.label_b = sj["labels"][1].get<std::string>();
        }
        // Need both a channel and a way to detect the crossing.
        s.enabled = !s.channel.empty() && s.has_threshold;
        if (!s.channel.empty() && !s.has_threshold)
            std::cerr << "replay_filter: split.channel set but no threshold / "
                         "level_low+level_high — split disabled\n";
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

int run(const std::vector<std::string> &input_files,
        const std::string &output_path,
        const std::string &cuts_path,
        const std::string &report_path,
        int run_number_override)
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

    // ---------- Phase 5: detect the split transition, label each side ----------
    // When `split` is enabled the run is partitioned into side A (before the
    // PV transition) and side B (after it), with a guard band around the
    // crossing whose checkpoints land in neither side.  `side[k]` ∈ {0, 1}
    // for A/B, or −1 for guard/excluded; with split off it is 0 everywhere,
    // so the keep/charge/output code below is one code path that simply runs
    // once (n_sides == 1) or twice (n_sides == 2).
    const int n_sides = cuts.split.enabled ? 2 : 1;
    std::vector<int> side(timeline.size(), 0);
    std::vector<int> scaler_side(scalers.size(), 0);
    std::vector<int> epics_side (epics_rows.size(), 0);

    bool   split_found     = false;
    int    split_cross_idx = -1;            // first B-side checkpoint
    int    guard_lo = (int)timeline.size(); // [guard_lo, guard_hi) = excluded
    int    guard_hi = (int)timeline.size();
    int32_t split_cross_evn = -1;
    double  split_cross_t   = std::numeric_limits<double>::quiet_NaN();
    double  split_pre_level = std::numeric_limits<double>::quiet_NaN();
    double  split_post_level= std::numeric_limits<double>::quiet_NaN();

    if (cuts.split.enabled) {
        const auto &S = cuts.split;
        // First checkpoint that actually carries a PV reading.
        size_t first_valid = timeline.size();
        for (size_t i = 0; i < timeline.size(); ++i)
            if (std::isfinite(timeline[i].split_pv)) { first_valid = i; break; }

        if (first_valid < timeline.size()) {
            const double v0 = timeline[first_valid].split_pv;
            split_pre_level = v0;
            const bool start_high = (v0 >= S.threshold);
            // Raw crossing: first checkpoint definitively on the far side.
            // With a hysteresis band the PV must actually reach the opposite
            // level; with a bare threshold a single crossing suffices (the
            // guard supplies the deadband against bounce).
            for (size_t i = first_valid; i < timeline.size(); ++i) {
                const double v = timeline[i].split_pv;
                if (!std::isfinite(v)) continue;
                bool now_far;
                if (S.has_levels)
                    now_far = start_high ? (v <= S.level_low) : (v >= S.level_high);
                else
                    now_far = start_high ? (v <  S.threshold) : (v >= S.threshold);
                if (now_far) { split_cross_idx = (int)i; break; }
            }
        }

        if (split_cross_idx >= 0) {
            split_found     = true;
            split_cross_evn = timeline[split_cross_idx].event_number;
            split_post_level= timeline[split_cross_idx].split_pv;
            if (timeline[split_cross_idx].ti_ticks > 0 && anchor_set)
                split_cross_t = (timeline[split_cross_idx].ti_ticks - ti_anchor)
                                * TI_TICK_SEC;

            if (S.guard == SplitCfg::Guard::Checkpoints) {
                const int N = std::max(0, S.guard_checkpoints);
                guard_lo = std::max(0, split_cross_idx - N);
                guard_hi = std::min((int)timeline.size(), split_cross_idx + N);
            } else {
                // Stability: extend the guard outward from the crossing until
                // the PV has settled within ε of each side's level.  pre_level
                // is the run's starting value; post_level is the value at the
                // crossing.  Symmetric ε band on both sides.
                const double eps = S.guard_epsilon;
                const double pre = split_pre_level, post = split_post_level;
                int lo = split_cross_idx, hi = split_cross_idx;
                // walk left: last A-side checkpoint still settled at `pre`
                for (int i = split_cross_idx - 1; i >= 0; --i) {
                    const double v = timeline[i].split_pv;
                    if (std::isfinite(v) && std::fabs(v - pre) <= eps) { lo = i + 1; break; }
                    if (i == 0) lo = 0;
                }
                // walk right: first B-side checkpoint settled at `post`
                for (int i = split_cross_idx; i < (int)timeline.size(); ++i) {
                    const double v = timeline[i].split_pv;
                    if (std::isfinite(v) && std::fabs(v - post) <= eps) { hi = i; break; }
                    if (i == (int)timeline.size() - 1) hi = (int)timeline.size();
                }
                guard_lo = std::max(0, lo);
                guard_hi = std::min((int)timeline.size(), std::max(hi, lo));
            }
        } else {
            std::cerr << "replay_filter: split channel '" << cuts.split.channel
                      << "' never crossed the configured threshold — side B "
                         "will be empty (no transition seen in this input)\n";
        }

        // Label every checkpoint, then back-map to load-order rows so the
        // scaler/epics writers can tag each row with its side.
        for (size_t k = 0; k < timeline.size(); ++k) {
            int sd;
            if      ((int)k <  guard_lo) sd = 0;   // before the ramp → A
            else if ((int)k >= guard_hi) sd = 1;   // after the ramp  → B
            else                         sd = -1;  // inside the guard → drop
            side[k] = sd;
            if (timeline[k].is_scaler) scaler_side[timeline[k].orig] = sd;
            else                       epics_side [timeline[k].orig] = sd;
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
    auto is_kept = [&](int s, int32_t evn) -> bool {
        const auto &kk = keep[s];
        if (kk.empty()) return false;
        auto it = std::upper_bound(
            kk.begin(), kk.end(), evn,
            [](int32_t e, const std::pair<int32_t, int32_t> &p) { return e < p.first; });
        if (it == kk.begin()) return false;
        --it;
        return evn > it->first && evn <= it->second;
    };

    // ---------- Phase 6: write the output(s) — one ROOT file + report per side ----------
    // ev_tree_name was detected in phase 3 above.  `write_side` does the whole
    // job for one side; with split off it runs once (side 0 → output_path),
    // with split on it runs twice with the file/report stems suffixed by the
    // side labels.  A row's `good` flag in a side's slow trees is its overall
    // verdict AND-ed with "belongs to this side", so running live_charge on a
    // side file reproduces that side's post-cut charge directly.
    const bool is_recon = (ev_tree_name == "recon");

    // Checkpoints per side / in the guard band — reported in the split block.
    int n_cp_side[2] = {0, 0}, n_cp_guard = 0;
    for (int sd : side) { if (sd < 0) ++n_cp_guard; else ++n_cp_side[sd]; }

    auto fmt_pct = [](double r) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(2) << (r * 100.0) << "%";
        return o.str();
    };

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
            n_in += n;
            for (Long64_t i = 0; i < n; ++i) {
                ev.ssp_raw.clear();
                ev.vtp_roc_tags.clear();
                ev.vtp_nwords.clear();
                ev.vtp_words.clear();
                ev.tdc_roc_tags.clear();
                ev.tdc_nwords.clear();
                ev.tdc_words.clear();
                t->GetEntry(i);
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
            Long64_t n = t->GetEntries();
            n_in += n;
            for (Long64_t i = 0; i < n; ++i) {
                ev.ssp_raw.clear();
                t->GetEntry(i);
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
                if (good && n_sides > 1)
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
                if (good && n_sides > 1)
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

    int n_pass_cp = 0, n_fail_cp = 0;
    for (const auto &cp : timeline) (cp.overall_pass ? n_pass_cp : n_fail_cp)++;
    json keep_intervals = json::array();
    for (const auto &p : keep[s]) keep_intervals.push_back({p.first, p.second});

    // Per-channel breakdown — number of slow-event checkpoints where this
    // channel's cut accepted vs rejected the value.  Helps the user see
    // immediately which cut is doing the rejecting (e.g. "beam current
    // killed 80% of points, livetime barely matters").
    std::map<std::string, std::pair<int, int>> per_channel;   // ch → {pass, fail}
    for (const auto &p : report_points) {
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
    const int     n_slow      = static_cast<int>(timeline.size());

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

    // Split metadata — present only when run-splitting is active.  Records
    // the transition this side was cut at, plus how many checkpoints landed
    // in each side and in the guard band, so every side file is self-
    // describing without needing the other side's report.
    if (cuts.split.enabled) {
        json sp = {
            {"enabled",             true},
            {"channel",             cuts.split.channel},
            {"side",                s},
            {"label",               s == 0 ? cuts.split.label_a : cuts.split.label_b},
            {"role",                s == 0 ? "before_transition" : "after_transition"},
            {"transition_found",    split_found},
            {"guard_mode",          cuts.split.guard == SplitCfg::Guard::Stability
                                     ? "stability" : "checkpoints"},
            {"n_checkpoints_side",  n_cp_side[s]},
            {"n_checkpoints_guard", n_cp_guard},
        };
        if (cuts.split.has_levels) {
            sp["level_low"]  = cuts.split.level_low;
            sp["level_high"] = cuts.split.level_high;
        }
        if (cuts.split.has_threshold) sp["threshold"] = cuts.split.threshold;
        if (cuts.split.guard == SplitCfg::Guard::Checkpoints) {
            sp["guard_checkpoints"] = cuts.split.guard_checkpoints;
        } else {
            sp["guard_epsilon"]     = cuts.split.guard_epsilon;
            sp["guard_consecutive"] = cuts.split.guard_consecutive;
        }
        if (split_found) {
            sp["transition_evn"]       = split_cross_evn;
            sp["transition_timestamp"] = std::isnan(split_cross_t)
                                          ? json(nullptr) : json(split_cross_t);
            sp["pre_level"]            = std::isnan(split_pre_level)
                                          ? json(nullptr) : json(split_pre_level);
            sp["post_level"]           = std::isnan(split_post_level)
                                          ? json(nullptr) : json(split_post_level);
        }
        report["split"] = std::move(sp);
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

    const std::string tag = (n_sides > 1)
        ? "[" + (s == 0 ? cuts.split.label_a : cuts.split.label_b) + "] " : "";
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

    // ---------- Dispatch: one output, or two split at the transition ----------
    if (n_sides == 1)
        return write_side(0, output_path, report_path);

    // Suffix the file/report stems with the side labels.  Handles the
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
    std::cerr << "replay_filter: split on '" << cuts.split.channel
              << "' → side " << cuts.split.label_a << " (before) + side "
              << cuts.split.label_b << " (after)";
    if (split_found) std::cerr << ", transition at event " << split_cross_evn;
    else             std::cerr << " — NO transition detected, side "
                               << cuts.split.label_b << " will be empty";
    std::cerr << "\n";
    int rc = write_side(0, with_label(output_path, cuts.split.label_a),
                           with_label(report_path, cuts.split.label_a));
    if (rc) return rc;
    return write_side(1, with_label(output_path, cuts.split.label_b),
                          with_label(report_path, cuts.split.label_b));
}

void usage(const char *prog)
{
    std::cerr <<
        "Usage: " << prog << " <input.root> [more.root ...]\n"
        "       -o <output.root>  -c <cuts.json> [-j <report.json>] [-r <run_num>]\n"
        "\n"
        "Filters replayed ROOT files by slow-control cuts (DSC2 livetime\n"
        "+ EPICS).  Writes a single ROOT file with the kept physics events\n"
        "and the full scaler/epics streams concatenated, plus a JSON report\n"
        "with per-(cut, slow-event) pass/fail status for chart plotting.\n"
        "\n"
        "With a \"split\" block the run is instead partitioned at a PV\n"
        "transition (e.g. target cell pressure stepping empty<->full) into\n"
        "TWO files <stem>_<labelA>.root and <stem>_<labelB>.root, each with\n"
        "its own report and live-charge integral; a guard band brackets the\n"
        "crossing so ramp events land in neither.\n"
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
        "      \"channel\": \"TGT:PRad:Cell_P\", \"threshold\": 1.0,\n"
        "      \"guard\": { \"mode\": \"checkpoints\", \"checkpoints\": 2 },\n"
        "      \"labels\": [\"full\", \"empty\"]\n"
        "    }\n"
        "  }\n"
        "  (split detection: a single \"threshold\", or a \"level_low\"+\n"
        "   \"level_high\" hysteresis band.  guard mode: \"checkpoints\" (+/- N\n"
        "   points) or \"stability\" (\"epsilon\"+\"consecutive\").)\n";
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    std::vector<std::string> inputs;
    std::string output, cuts_path, report_path;
    int run_override = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:c:j:r:h")) != -1) {
        switch (opt) {
        case 'o': output       = optarg;           break;
        case 'c': cuts_path    = optarg;           break;
        case 'j': report_path  = optarg;           break;
        case 'r': run_override = std::atoi(optarg); break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    for (int i = optind; i < argc; ++i) inputs.push_back(argv[i]);

    if (inputs.empty() || output.empty() || cuts_path.empty()) {
        usage(argv[0]);
        return 1;
    }
    if (report_path.empty()) {
        auto dot = output.rfind('.');
        report_path = (dot == std::string::npos)
            ? output + ".report.json"
            : output.substr(0, dot) + ".report.json";
    }
    return run(inputs, output, cuts_path, report_path, run_override);
}
