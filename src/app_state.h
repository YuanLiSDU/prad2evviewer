#pragma once
//=============================================================================
// app_state.h — shared application state for the event viewer/monitor
//
// Owns all configuration, accumulated data (histograms, LMS), and HyCal system.
// The ViewerServer maintains two AppState instances (file and online) with
// separate accumulators but identical configuration.
//=============================================================================

// forward declaration (full definition in data_source.h)
struct ReconEventData;

#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "EpicsStore.h"
#include "DaqConfig.h"
#include "DscData.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "GemSystem.h"
#include "GemCluster.h"
#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"
#include "viewer_utils.h"

#include <nlohmann/json.hpp>

#include <array>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <mutex>
#include <atomic>
#include <iostream>

// Resolve a JSON array of bit numbers / names to a uint32_t mask.
// Names are looked up in `bits_defs` ([{bit,name,label},…]).  Used by both
// TriggerFilter (trigger_bits.json) and PeakFilter (peak_quality_bits).
inline uint32_t bitsMaskFromArray(const nlohmann::json &arr,
                                   const nlohmann::json &bits_defs,
                                   const char *who = "BitMask") {
    if (!arr.is_array() || arr.empty()) return 0;
    uint32_t m = 0;
    for (auto &item : arr) {
        if (item.is_number()) {
            m |= (1u << item.get<int>());
        } else if (item.is_string()) {
            auto s = item.get<std::string>();
            bool found = false;
            for (auto &def : bits_defs) {
                if (def.value("name", "") == s || def.value("label", "") == s) {
                    m |= (1u << def.value("bit", 0));
                    found = true;
                    break;
                }
            }
            if (!found)
                std::cerr << who << ": unknown bit name '" << s << "'\n";
        }
    }
    return m;
}

inline uint32_t bitsMaskFromKey(const nlohmann::json &section, const char *key,
                                 const nlohmann::json &bits_defs,
                                 const char *who = "BitMask") {
    if (!section.contains(key)) return 0;
    return bitsMaskFromArray(section[key], bits_defs, who);
}

// Reverse-lookup: encode a bitmask as a JSON array of bit names against
// bits_defs.  Bits not found in bits_defs are emitted as numeric indices so
// nothing is silently dropped.
inline nlohmann::json bitsMaskToNames(uint32_t mask,
                                       const nlohmann::json &bits_defs) {
    nlohmann::json out = nlohmann::json::array();
    if (!mask) return out;
    for (int b = 0; b < 32; ++b) {
        if (!(mask & (1u << b))) continue;
        bool found = false;
        for (auto &def : bits_defs) {
            if (def.value("bit", -1) == b) {
                out.push_back(def.value("name", ""));
                found = true;
                break;
            }
        }
        if (!found) out.push_back(b);
    }
    return out;
}

// Trigger filter: reject overrides accept. accept==0 means accept all.
struct TriggerFilter {
    uint32_t accept = 0;  // 0 = accept all
    uint32_t reject = 0;  // 0 = reject none

    bool operator()(uint32_t bits) const {
        if (reject && (bits & reject)) return false;
        return accept == 0 || (bits & accept);
    }

    // parse from JSON section containing accept_trigger_bits / reject_trigger_bits
    // values can be: bit numbers (8, 24), or names ("LMS", "Pulser") resolved
    // against the trigger_bits_def lookup table from trigger_bits.json
    void parse(const nlohmann::json &section,
               const nlohmann::json &bits_defs = nlohmann::json::array()) {
        accept = bitsMaskFromKey(section, "accept_trigger_bits", bits_defs, "TriggerFilter");
        reject = bitsMaskFromKey(section, "reject_trigger_bits", bits_defs, "TriggerFilter");
    }

    nlohmann::json toJson() const {
        return {{"trigger_accept", accept}, {"trigger_reject", reject}};
    }

    friend std::ostream &operator<<(std::ostream &os, const TriggerFilter &f) {
        return os << "trigger accept=0x" << std::hex << f.accept
                  << " reject=0x" << f.reject << std::dec;
    }
};

// Per-peak filter for the Waveform Tab histograms (and any consumer that opts
// into it).  Each axis is an optional [min, max] range — missing bound means
// no constraint on that side.  Quality bits use accept/reject masks resolved
// against AppState::peak_quality_bits_def.  `enable=false` makes the filter a
// no-op (predicate returns true unconditionally) — the GUI's "apply" checkbox
// drives this flag.
struct PeakFilter {
    // Default `enable=true` so the JSON-configured filter is active on
    // startup — preserves today's behavior (time-cut always applied) without
    // needing a runtime flag in monitor_config.json.  The GUI "apply"
    // checkbox toggles this at runtime.
    bool enable = true;
    std::optional<float> time_min, time_max;
    std::optional<float> integral_min, integral_max;
    std::optional<float> height_min, height_max;
    uint32_t q_accept = 0;     // 0 = accept any
    uint32_t q_reject = 0;     // 0 = reject none

    bool operator()(const fdec::Peak &pk) const {
        if (time_min     && pk.time     < *time_min)     return false;
        if (time_max     && pk.time     > *time_max)     return false;
        if (integral_min && pk.integral < *integral_min) return false;
        if (integral_max && pk.integral > *integral_max) return false;
        if (height_min   && pk.height   < *height_min)   return false;
        if (height_max   && pk.height   > *height_max)   return false;
        if (q_reject && (pk.quality & q_reject))         return false;
        if (q_accept && !(pk.quality & q_accept))        return false;
        return true;
    }

    void parse(const nlohmann::json &filter,
               const nlohmann::json &quality_bits_def);
    nlohmann::json toJson(const nlohmann::json &quality_bits_def) const;
};

struct AppState {
    // ---- Configuration (set once at startup, then read-only) ---------------
    HistConfig hist_cfg;
    TriggerFilter waveform_trigger;
    int hist_nbins   = 0;
    int pos_nbins    = 0;
    int height_nbins = 0;

    // Waveform-Tab peak filter (height/integral/time hist + show overlays).
    // `enable` is the GUI "apply" checkbox; runtime-only, not in monitor_config.json.
    // `peak_filter_default` snapshots the JSON-configured filter at startup
    // so the Cut-Settings "Reset" button can restore the file values without
    // a server round-trip.
    // peak_quality_bits_def is the [{bit,name,label},…] palette exposed via /api/config.
    PeakFilter      peak_filter;
    PeakFilter      peak_filter_default;
    nlohmann::json  peak_quality_bits_def = nlohmann::json::array();

    // Gain-monitoring (LMS/Alpha) time window — distinct from peak_filter so
    // the calibration tracks aren't perturbed by Waveform-Tab edits.
    float lms_time_min = 160.0f;
    float lms_time_max = 220.0f;

    // reference lines for plots (pass-through to frontend)
    // key = plot name, value = array of {axis, pos, lw, ls, color}
    nlohmann::json ref_lines = nlohmann::json::object();

    // trigger definitions (pass-through to frontend)
    nlohmann::json trigger_bits_def = nlohmann::json::array();  // FP bits: [{bit, name, label}, ...]
    nlohmann::json trigger_type_def = nlohmann::json::array();  // main trigger: [{type, tag, name, label, primary_bit}, ...]

    evc::DaqConfig daq_cfg;
    fdec::HyCalSystem hycal;
    fdec::ClusterConfig cluster_cfg;

    // Per-channel pulse-template store for the NNLS pile-up deconvolver.
    // Loaded by init() from `daq_cfg.wave_cfg.nnls_deconv.template_file`
    // (resolved against db_dir).  When invalid (file missing / parse
    // failure) the deconv path silently falls back to the local-maxima
    // peak heights — every WaveAnalyzer in the app picks this up via
    // SetTemplateStore() during the same init.
    fdec::PulseTemplateStore template_store;

    // GEM system
    gem::GemSystem gem_sys;
    gem::GemCluster gem_clusterer;
    bool gem_enabled = false;       // true if gem_map.json loaded successfully

    // GEM per-detector lab-frame transform (same type as HyCal)
    std::vector<DetectorTransform> gem_transforms;  // indexed by detector id

    // GEM occupancy (accumulated per-detector 2D histograms)
    static constexpr int GEM_OCC_NX = 50;
    static constexpr int GEM_OCC_NY = 30;
    std::vector<Histogram2D> gem_occupancy;  // one per detector

    std::unordered_map<int, int> roc_to_crate;  // ROC tag → crate index
    nlohmann::json crate_roc_json;              // crate→ROC tag JSON
    nlohmann::json base_config;                 // modules, daq, crate_roc for /api/config

    // LMS config
    TriggerFilter lms_trigger;
    TriggerFilter alpha_trigger;          // Am-241 alpha-source trigger (LMS ref channels)
    float    lms_warn_thresh  = 0.1f;
    float    lms_warn_min_mean = 100.f;  // warn if mean below this
    int      lms_max_history  = 5000;

    // Gain-drift detection: compare current LMS mean against a baseline LMS dat
    // file (produced by prad2ana_gain_monitor for a known-good calibration run).
    // Drift is reference-normalized so a global lamp-output / FADC-integration
    // change between baseline and current run does not register as drift.
    //
    // For each module M with a baseline:
    //   lamp_scale = current_LMS_mean[ref] / baseline_LMS_peak[ref]   (using
    //                lms_drift_ref_channel; falls back to first reference)
    //   expected   = baseline_LMS_peak[M] * lamp_scale
    //   drift      = current_LMS_mean[M] / expected
    //
    // A module is flagged "drift" when drift < lms_drift_low or > lms_drift_high.
    // Drift is treated as a top-priority error in the LMS table / report / tab dot
    // (ranks above the existing stability/floor warn).
    std::string lms_drift_baseline_file;        // path resolved at init (empty = disabled)
    std::string lms_drift_ref_channel;          // ref-channel name (e.g. "LMS2"); empty = first ref
    // Per-module-type drift thresholds.  PbWO4 (W) and PbGlass (G) have very
    // different LMS-coupling and PMT aging behaviours; in the 024352→024459
    // baseline span, healthy W modules sat within ±5% of 1.0 while the entire
    // PbGlass population shifted up by ~21%.  A single threshold either misses
    // genuine W drift or floods the alarm list with the G systematic, so each
    // type gets its own [low, high] pair.  Defaults are tight on W (where
    // healthy width is ~5%, so ±15% is ~3σ from the bulk) and looser on G to
    // accommodate the systematic — re-tune in monitor_config.json once a fresh
    // baseline lands and the G shift is folded into the reference.
    float       lms_drift_low_w  = 0.85f;       // PbWO4: flag drift below this
    float       lms_drift_high_w = 1.15f;       // PbWO4: flag drift above this
    float       lms_drift_low_g  = 0.7f;        // PbGlass: flag drift below this
    float       lms_drift_high_g = 1.5f;        // PbGlass: flag drift above this
    // Suppress drift warnings for these module types entirely.  Drift values
    // are still computed and shown in the table (so operators can inspect them
    // on demand), but the module's state never escalates to "drift" — useful
    // when a known systematic between baseline calibrations would otherwise
    // flood the alarm list.  Type names follow hycal_map.json "t" field
    // (PbGlass / PbWO4 / LMS / Veto); for raw typed JSON we parse via
    // HyCalSystem::parse_type at init.  Stored as the ModuleType enum cast
    // to int for O(1) lookup.
    std::set<int> lms_drift_suppress_types;
    std::vector<std::string> lms_drift_suppress_type_names;  // echoed in /api/config
    std::unordered_map<std::string, float> lms_baseline_peak;       // module name → baseline LMS peak
    std::unordered_map<std::string, float> lms_baseline_ref_peak;   // ref-channel name → baseline LMS peak

    // LMS reference channels (for normalization)
    struct LmsRefChannel {
        std::string name;
        int module_index = -1;
    };
    std::vector<LmsRefChannel> lms_ref_channels;

    // Latest per-channel integrals for ref correction (LMS_signal / Alpha_signal).
    // Updated unconditionally regardless of lms_max_history saturation, so the
    // correction factor always reflects the most recent readings.  Only the ref
    // module entries are read by the correction code; others are bookkeeping.
    std::map<int, float> latest_lms_integral;     // module_index → latest LMS-trigger integral
    std::map<int, float> latest_alpha_integral;   // module_index → latest Alpha-trigger integral

    // Monitor-mode status panel: shell-poll metrics shown in the header.
    // All three (livetime, beam energy, beam current) follow the same pattern
    // — run `command`, parse first float from stdout, display with `unit`.
    // ViewerServer reads these after init() and runs one shared poll thread
    // that ticks each metric on its own poll_sec; an empty command skips
    // that metric.  Avoids a build-time EPICS dependency by shelling out to
    // whatever tool the host provides (typically caget).
    //
    // livetime: healthy/warning are percent thresholds for color
    // (≥ healthy → green, ≥ warning → orange, otherwise red).
    std::string livetime_cmd;
    std::string livetime_unit       = "%";
    int         livetime_poll_sec   = 30;
    float       livetime_healthy    = 90.f;
    float       livetime_warning    = 80.f;

    // Beam status (energy, current).  trip_warn for current colors red when
    // the reading drops below the threshold (beam-trip indicator).
    struct ShellMetric {
        std::string command;
        std::string unit;
        int         poll_sec        = 5;
        bool        has_trip_warn   = false;
        float       trip_warn_below = 0.f;
    };
    ShellMetric beam_energy_status;
    ShellMetric beam_current_status;

    // Measured DAQ livetime from DSC2 scalers in the EVIO stream
    // (gated / ungated, expressed as percent).  Bank tag, slot, source,
    // channel live in daq_cfg.dsc_scaler so the decoder/format details stay
    // in prad2dec.  dsc_prev_* hold the gated/ungated values from the
    // previous SYNC so the displayed live time reflects the most recent
    // ~2 s window (DSC2 counts accumulate; ratio of cumulative values
    // averages over the whole run).
    std::atomic<double>   measured_livetime{-1.0};
    std::atomic<uint32_t> dsc_prev_gated{0};
    std::atomic<uint32_t> dsc_prev_ungated{0};

    // online refresh rates (ms), served to frontend
    int refresh_event_ms = 200;
    int refresh_ring_ms  = 500;
    int refresh_hist_ms  = 2000;
    int refresh_lms_ms   = 2000;

    // Physics / coordinate config
    float target_x=0, target_y=0, target_z=0;  // target position in lab frame (mm)
    DetectorTransform hycal_transform;           // HyCal position + tilting
    float ea_angle_min=0.f, ea_angle_max=8.f, ea_angle_step=0.2f;   // degrees
    float ea_energy_min=0.f, ea_energy_max=3000.f, ea_energy_step=100.f; // MeV
    // Single-source beam energy: MBSY2C_energy from EPICS overrides; runinfo is fallback.
    // Read in physics paths (cluster filling, plots); written by init() (runinfo) and
    // processEpics() (EPICS). Atomic so processEpics can update without holding data_mtx.
    std::atomic<float> beam_energy{2200.f};                  // MeV
    float       beam_energy_runinfo = 0.f;                   // fallback (loaded from runinfo)
    std::string beam_energy_epics_channel = "MBSY2C_energy"; // EPICS channel name
    float       beam_energy_min_valid = 100.f;               // ignore stale/zero EPICS reads
    TriggerFilter physics_trigger;

    // Møller selection config
    float moller_energy_tol = 0.1f;     // energy sum within this fraction of beam_energy
    float moller_angle_min  = 1.0f;     // deg — require one cluster in this range
    float moller_angle_max  = 1.2f;     // deg
    // Møller XY histogram
    float moller_xy_x_min=-600.f, moller_xy_x_max=600.f, moller_xy_x_step=5.f;  // mm
    float moller_xy_y_min=-600.f, moller_xy_y_max=600.f, moller_xy_y_step=5.f;  // mm

    // HyCal cluster-hit XY (single-cluster ep-elastic candidates) — cuts + hist
    int   hxy_n_clusters      = 1;        // require Ncl == this
    float hxy_energy_frac_min = 0.9f;     // require E_cl >= frac * beam_energy
    int   hxy_nblocks_min     = 5;
    int   hxy_nblocks_max     = 20;
    float hxy_x_min=-600.f, hxy_x_max=600.f, hxy_x_step=5.f;  // mm
    float hxy_y_min=-600.f, hxy_y_max=600.f, hxy_y_step=5.f;  // mm

    // GEM↔HyCal matching: per-detector residuals filled when ep candidate fires.
    // The cut is parametric: cut = match_nsigma * sqrt(sigma_HC² + sigma_GEM²),
    // where sigma_HC = hycal.PositionResolution(E) and sigma_GEM = gem_pos_res[d]
    // (both projected to the residual plane).  See reconstruction_config.json:matching.
    bool  gem_match_require_ep = true;    // gate on hxy_* selection (clean track)
    float gem_match_nsigma     = 3.f;     // residual cut in σ_total
    float gem_resid_min = -25.f, gem_resid_max = 25.f, gem_resid_step = 0.5f;  // mm

    // Per-detector GEM position resolution (mm), parsed from
    // reconstruction_config.json:matching:gem_pos_res.  HyCal's energy-
    // dependent resolution lives on HyCalSystem (PositionResolution(E)).
    std::vector<float> gem_pos_res;

    // GEM tracking-efficiency monitor — leave-one-out per detector.
    //
    // For each test detector D, the OTHER 3 GEMs + HyCal anchor a straight
    // line; the line is then projected to D and a hit is searched within the
    // matching window.  D contributes nothing to the anchor, so the
    // prediction at D is unbiased.  Three modes select how the anchor is
    // built:
    //   loo-target-seed  single seed = (target → HyCal); cheapest.  Default.
    //   loo              for each OTHER GEM hit, draw seed (HyCal → that
    //                    hit), pick lowest-χ² anchor.  No vertex assumption.
    //   loo-target-in    same GEM-seeding as `loo` but the fit also includes
    //                    (target_x, target_y, target_z) as a soft constraint
    //                    weighted by σ_target_{x,y} ⊕ slope·σ_target_z.
    // Per-detector counters: `gem_eff_den[D]` is the number of anchors that
    // succeeded for test detector D, `gem_eff_num[D]` is the subset of those
    // where D actually had a hit at the prediction.
    enum class GemEffLooMode { TargetSeed = 0, Loo = 1, LooTargetIn = 2 };
    GemEffLooMode gem_eff_loo_mode  = GemEffLooMode::TargetSeed;
    float gem_eff_min_cluster_energy = 500.f;  // MeV — HyCal cluster gate for the seed
    float gem_eff_match_nsigma       = 3.f;
    float gem_eff_max_chi2           = 3.5f;
    int   gem_eff_max_hits_per_det   = 50;
    int   gem_eff_min_denom          = 20;
    float gem_eff_healthy            = 90.f;
    float gem_eff_warning            = 70.f;
    // Beam-spot transverse + target-length sigmas, used only by loo-target-in
    // (and only matter when the fit actually pulls toward the target).
    float gem_eff_target_sigma_x    = 1.0f;   // mm
    float gem_eff_target_sigma_y    = 1.0f;   // mm
    float gem_eff_target_sigma_z    = 20.0f;  // mm
    // Per-detector efficiency-vs-position grid: predicted (local_x, local_y)
    // at each LOO test detector, binned over the detector's active area.
    // Numerator counts "matched" tests; denominator counts every successful
    // anchor.  GUI divides per-bin to show local efficiency (where the
    // overlap region is high vs. low).
    int gem_eff_grid_nx = 30;
    int gem_eff_grid_ny = 30;

    // EPICS config
    int   epics_max_history = 5000;
    float epics_warn_thresh  = 0.1f;
    float epics_alert_thresh = 0.2f;
    int   epics_min_avg_pts  = 10;
    int   epics_mean_window  = 20;   // compute mean from most recent N snapshots
    std::vector<std::vector<std::string>> epics_default_slots;  // per-slot channel lists

    // Elog config
    std::string elog_url;
    std::string elog_logbook;
    std::string elog_author;
    std::vector<std::string> elog_tags;
    std::string elog_cert;         // SSL client certificate path
    std::string elog_key;          // SSL client key path

    // Master enable for the on-demand auto-report pipeline. Reflected
    // in the header status pill on every connected client; gates
    // dispatchCapture on END / run-change in the ET reader thread.
    bool        auto_report_enabled = true;
    // Dry-run guard. When false, /api/elog/post still saves the report
    // XML + image attachments under auto_report_local_save_dir but
    // skips the curl upload — used to validate the pipeline end-to-end
    // without polluting the (non-deletable) logbook.
    bool        auto_report_post_to_elog = false;
    // Required.  Each capture is mirrored under
    // <local_save_dir>/run_NNNNNN/ so the body can be inspected and
    // replayed manually with curl.
    std::string auto_report_local_save_dir;
    // Server-side dedup window (ms) between auto posts for the same
    // run; absorbs the END + run-change double-fire and multi-browser
    // races. Default 15 min.
    int         auto_report_min_interval_ms = 900000;

    // color range defaults: key "tab:metric" → [min, max]
    std::map<std::string, std::pair<float, float>> color_range_defaults;

    // cluster config
    TriggerFilter cluster_trigger;
    float    adc_to_mev        = 1.0f;
    float    cl_hist_min       = 0.f;
    float    cl_hist_max       = 3000.f;
    float    cl_hist_step      = 10.f;
    // nclusters_hist range is float so the user can shift bin edges by half
    // a step (default 0.5 .. 10.5 / 1 → bin centers land on 1, 2, …, 10).
    // The bucket index that a given Ncl event falls in is reused as the
    // index into cluster_energy_hist_by_ncl / nblocks_hist_by_ncl, so the
    // dependent histograms can be filtered to "events with this many
    // clusters" by clicking a bar in the GUI.
    float    nclusters_hist_min  = 0.5f;
    float    nclusters_hist_max  = 10.5f;
    float    nclusters_hist_step = 1.0f;
    int      nblocks_hist_min    = 0;
    int      nblocks_hist_max    = 40;
    int      nblocks_hist_step   = 1;
    // Raw (per-event) HyCal energy sum: total energy deposited across all
    // modules before clustering.  Wider range than energy_hist because a
    // single Moller event can deposit ~Eb across two clusters.
    float    raw_energy_hist_min  = 0.f;
    float    raw_energy_hist_max  = 6000.f;
    float    raw_energy_hist_step = 20.f;

    // ---- Event filters (loaded from external JSON via loadFilter) -----------
    // trigger_type filter: if enabled, only events with trigger_type in accept pass
    struct TriggerTypeFilter {
        bool enable = false;
        std::vector<uint8_t> accept;  // empty when disabled
        bool operator()(uint8_t tt) const {
            if (!enable || accept.empty()) return true;
            for (auto t : accept) if (t == tt) return true;
            return false;
        }
    } trigger_type_filter;

    WaveformFilter waveform_filter;
    ClusterFilter  cluster_filter;
    // resolved indices for fast filter evaluation
    std::unordered_set<std::string> filter_wf_keys;   // DAQ keys ("roc_slot_ch")
    std::unordered_set<int> filter_cl_includes;        // module indices for includes check
    std::unordered_set<int> filter_cl_centers;          // module IDs (PrimEx) for center check

    // ---- Accumulated data (guarded by data_mtx) ----------------------------
    mutable std::mutex data_mtx;
    std::map<std::string, Histogram> histograms;
    std::map<std::string, Histogram> pos_histograms;
    std::map<std::string, Histogram> height_histograms;
    std::map<std::string, int>       occupancy;
    std::map<std::string, int>       occupancy_tcut;
    std::atomic<int>                 events_processed{0};

    Histogram cluster_energy_hist;
    Histogram nclusters_hist;
    Histogram nblocks_hist;
    Histogram raw_energy_hist;     // sum of all module energies per event
    // Per-Ncl-bucket dependent histograms.  Index = the bucket of the
    // event's Ncl in nclusters_hist (i.e. floor((Ncl - min) / step)),
    // so size matches nclusters_hist.bins.size() once init() runs.
    // Filled alongside the unfiltered hists; the GUI selects which one
    // to plot when the user clicks an Ncl bar.
    std::vector<Histogram> cluster_energy_hist_by_ncl;
    std::vector<Histogram> nblocks_hist_by_ncl;
    Histogram2D energy_angle_hist;
    Histogram2D moller_xy_hist;
    Histogram2D hycal_xy_hist;       // single-cluster ep-elastic candidates (cuts in hxy_*)
    int         hycal_xy_events = 0; // events passing hycal_cluster_hit cuts

    // GEM↔HyCal residuals (one Histogram per GEM detector for each axis).
    // gem_match_events = events that contributed; gem_match_hits[d] = matched hits per det.
    std::vector<Histogram> gem_dx_hist;   // size = nDetectors
    std::vector<Histogram> gem_dy_hist;
    int                    gem_match_events = 0;
    std::vector<int>       gem_match_hits;  // per-det count of in-window hits

    // GEM efficiency counters: per-detector numerator, single shared
    // denominator (incremented once per good track).  See class-level comment
    // above for the definition of a good track.
    std::vector<int> gem_eff_num;   // per-detector numerator
    std::vector<int> gem_eff_den;   // per-detector denominator (LOO test count)
    // Per-stage breakdown for the loo-target-seed mode: lets us diff the
    // server's anchor pipeline against the offline gem_eff_audit.py at
    // each gate.  Only the target-seeded path is instrumented because
    // the GEM-seeded modes try multiple seeds per (HyCal, test_d) and
    // the per-stage counts would not be 1:1 comparable.
    int gem_eff_diag_call[4]       = {0,0,0,0};
    int gem_eff_diag_3matched[4]   = {0,0,0,0};
    int gem_eff_diag_pass_chi2[4]  = {0,0,0,0};
    int gem_eff_diag_pass_resid[4] = {0,0,0,0};
    static constexpr int GEM_EFF_MAX_DETS = 4;
    // Snapshot of the last good track for the "last good event" panel.
    // Stores the single fit + per-detector status (used in fit, prediction,
    // residual) so the frontend can draw the track and per-detector markers.
    struct GemEffSnapshot {
        bool  valid    = false;
        int   event_id = -1;
        // HyCal cluster lab-frame xyz — anchor of the fit.
        float hycal_x = 0.f, hycal_y = 0.f, hycal_z = 0.f;
        // Lab-frame fit line: x(z) = ax + bx·z, y(z) = ay + by·z
        float chi2_per_dof = -1.f;
        float ax = 0.f, bx = 0.f, ay = 0.f, by = 0.f;
        struct Det {
            bool  used_in_fit = false;       // matched within match_window of seed line
            bool  hit_present = false;       // hit_lab_* is valid
            float hit_lab_x = 0.f, hit_lab_y = 0.f, hit_lab_z = 0.f;
            bool  inside = false;            // predicted point inside active area
            float predicted_lab_x = 0.f, predicted_lab_y = 0.f, predicted_lab_z = 0.f;
            float predicted_local_x = 0.f, predicted_local_y = 0.f;
            float resid_dx = 0.f, resid_dy = 0.f; // hit_local - predicted_local (only if used_in_fit)
        };
        Det dets[GEM_EFF_MAX_DETS];
        // Closest approach of the fit line to the lab z-axis: minimizes
        // r²(z)=x(z)²+y(z)², so z_lab = -(ax·bx + ay·by) / (bx² + by²).
        // z_offset = z_lab - target_z is the inferred vertex displacement
        // from nominal target z; this is what the histogram accumulates.
        float z_target_lab    = 0.f;
        float z_target_offset = 0.f;
        bool  z_target_valid  = false;
    };
    GemEffSnapshot gem_eff_snapshot;
    // Per-detector efficiency-vs-position grid (predicted local coords of the
    // LOO test point).  Each detector has the same (nx,ny) bin count; bin
    // step = det.size / N so the grid spans the full active area.
    std::vector<Histogram2D> gem_eff_grid_num;
    std::vector<Histogram2D> gem_eff_grid_den;
    // Inferred vertex-z spread (DOCA to z-axis − target_z), one entry per
    // matched event.  Range/step are configurable via
    // monitor_config.json: gem.efficiency.z_target_hist.{min,max,step}.
    Histogram gem_eff_z_target_hist;
    float gem_eff_z_target_min  = -50.f;   // mm
    float gem_eff_z_target_max  =  50.f;   // mm
    float gem_eff_z_target_step =   1.f;   // mm
    int         moller_events = 0;
    int       cluster_events_processed = 0;

    // ---- LMS data (guarded by lms_mtx) -------------------------------------
    mutable std::mutex lms_mtx;
    std::map<int, std::vector<LmsEntry>> lms_history;
    std::atomic<int> lms_events{0};
    uint64_t lms_first_ts = 0;

    // Sync reference point for absolute time display
    // sync_unix = absolute time, sync_rel_sec = relative time on LMS axis
    uint32_t sync_unix    = 0;
    double   sync_rel_sec = 0.;
    uint32_t pending_sync_unix = 0;  // stashed until first LMS event arrives
    uint64_t pending_sync_ti   = 0;

    // ---- EPICS data (guarded by epics_mtx) ----------------------------------
    mutable std::mutex epics_mtx;
    epics::EpicsStore epics;
    std::atomic<int> epics_events{0};

    // ---- GEM calibration revision -------------------------------------------
    // Bumped on any change to per-APV pedestal noise (e.g. LoadPedestals).
    // The frontend caches /api/gem/calib (which carries this rev) and
    // refetches when it sees a different value embedded in /api/gem/apv/<n>.
    // Threshold (zs_sigma) changes do NOT bump this — the frontend reads
    // zs_sigma per event so band tracks the encoded hits.
    std::atomic<int> gem_calib_rev{0};

    // ---- Initialization (call once at startup) -----------------------------

    // Load all configs from db_dir.  Empty filename ⇒ auto-find in db_dir:
    //   daq_config_file   → daq_config.json   (DAQ + raw decoding)
    //   monitor_config_file → monitor_config.json (GUI / online server)
    //   recon_config_file → reconstruction_config.json (runinfo + clustering)
    void init(const std::string &db_dir,
              const std::string &daq_config_file,
              const std::string &monitor_config_file = "",
              const std::string &recon_config_file = "");

    // ---- Per-event processing ----------------------------------------------

    // Run the GEM clear + ProcessEvent + Reconstruct pipeline without
    // touching accumulators.  Used by the on-demand APV endpoint so a
    // re-request for an event already accumulated still leaves gem_sys
    // populated with that event's per-APV working buffers.
    void prepareGemForView(const ssp::SspEventData &ssp_evt);

    // Process GEM SSP data for one event. Call after DecodeEvent with ssp_evt.
    // Calls prepareGemForView, then accumulates occupancy + histograms.
    void processGemEvent(const ssp::SspEventData &ssp_evt);

    // Process one fully-decoded event: histograms + clustering + LMS.
    // Single pass over all channels (analyzes each channel once).
    // Thread-safe (acquires data_mtx + lms_mtx internally).
    void processEvent(fdec::EventData &event,
                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Process a pre-computed recon event (from ROOT recon files).
    // Fills cluster/physics histograms from pre-computed clusters.
    void processReconEvent(const struct ReconEventData &recon);

    // Project a lab-frame point through the target onto the HyCal local
    // plane (z_local = 0). Returns HyCal-local (px, py).
    // Used for GEM hit overlays and matching residuals.
    void projectToHyCalLocal(float Gx, float Gy, float Gz,
                             float &px, float &py) const;

    // Encode one decoded event as JSON.
    // include_samples=false: summaries only (peaks + pedestal, ~20KB).
    // include_samples=true:  full waveforms included (~800KB, used for ring buffer).
    nlohmann::json encodeEventJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres,
                                   bool include_samples = false);

    // Encode raw waveform for a single channel (key = "roc_slot_ch").
    nlohmann::json encodeWaveformJson(fdec::EventData &event, const std::string &chan_key,
                                      fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Compute clusters for one decoded event, return JSON response.
    nlohmann::json computeClustersJson(fdec::EventData &event, int ev_id,
                                       fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

    // Encode pre-computed recon clusters (from ROOT recon files) as JSON.
    nlohmann::json encodeReconClustersJson(const struct ReconEventData &recon, int ev_id);

    // Record a sync event's absolute time. Call when a Sync event is scanned.
    // last_ti_ts is the TI timestamp of the most recent physics event.
    void recordSyncTime(uint32_t unix_time, uint64_t last_ti_ts);

    // ---- EPICS processing ---------------------------------------------------
    void processEpics(const std::string &text, int32_t event_number, uint64_t timestamp);
    void clearEpics();        // locks epics_mtx

    // ---- GEM tracking efficiency monitor (called per HyCal cluster) --------
    // Pass A: GEM0 seed → tests {1,2,3}.  Pass B: GEM1 seed → tests {0}.
    // hits_by_det[d] = lab-frame (x,y,z) of every reconstructed GEM-d hit
    // available for this event (capped internally to gem_eff_max_hits_per_det).
    // Updates gem_eff_num/den, residual histograms, and gem_eff_snapshot.
    using LabHit = std::array<float, 3>;
    void runGemEfficiency(int event_id,
                          float hcx, float hcy, float hcz, float hc_energy,
                          const std::vector<std::vector<LabHit>> &hits_by_det);
    void clearGemEfficiency();   // counters + snapshot (data_mtx already held)
    void initGemEfficiency();    // size num/den/residuals (called from init())
    nlohmann::json gemEffSnapshotJson() const;  // assumes data_mtx held

    // ---- DSC2 scaler processing --------------------------------------------
    // Update measured_livetime from an already-decoded DSC2 record.  Caller
    // is expected to obtain it via EvChannel::Dsc() (which routes through
    // dsc::Dsc2Decoder for both bank-format detection and the configured
    // source/channel selection).  No-op when `dsc.present` is false or no
    // ungated counts have been seen yet.
    void processDsc(const dsc::DscEventData &dsc);

    // ---- Clearing ----------------------------------------------------------
    void clearHistograms();   // locks data_mtx
    void clearLms();          // locks lms_mtx

    // ---- API response builders (thread-safe) -------------------------------
    // type: 0=integral, 1=position, 2=height
    nlohmann::json apiHist(int type, const std::string &key) const;
    nlohmann::json apiClusterHist() const;
    nlohmann::json apiOccupancy() const;
    nlohmann::json apiColorRanges() const;
    nlohmann::json apiLmsSummary(int ref_index = -1) const;
    nlohmann::json apiLmsModule(int module_index, int ref_index = -1) const;
    nlohmann::json apiLmsRefChannels() const;
    // Drift-detection is gated on having both per-module and ref-channel
    // baselines: lamp_scale needs the ref to cancel global LMS-pulser drift.
    bool driftEnabled() const {
        return !lms_baseline_peak.empty() && !lms_baseline_ref_peak.empty();
    }
    nlohmann::json apiEnergyAngle() const;
    nlohmann::json apiMoller() const;
    nlohmann::json apiHycalXY() const;
    nlohmann::json apiGemResiduals() const;
    nlohmann::json apiGemEfficiency() const;
    nlohmann::json apiEpicsChannels() const;
    nlohmann::json apiEpicsChannel(const std::string &name) const;
    nlohmann::json apiEpicsBatch(const std::vector<std::string> &names) const;
    nlohmann::json apiEpicsLatest() const;
    nlohmann::json apiGemHits() const;
    nlohmann::json apiGemConfig() const;
    nlohmann::json apiGemOccupancy() const;
    // Per-event APV waveform dump (for the GEM APV monitor tab).
    // Caller must have just populated gem_sys with this event (e.g. via
    // processGemEvent or accumulate); this method only reads.  Raw ADC
    // samples come from ssp_evt; processed values + ZS hit mask come
    // from gem_sys's per-APV working buffer.
    nlohmann::json apiGemApv(const ssp::SspEventData &ssp_evt, int evnum) const;

    // One-shot calibration payload for the GEM APV tab: returns
    //   {rev, zs_sigma, apvs:[{id, noise:[128]}, ...]}
    // The frontend caches this and refetches only when the calib_rev
    // embedded in /api/gem/apv/<n> diverges from the cached value.
    nlohmann::json apiGemCalib() const;

    // Update the software N-sigma cut on this AppState's gem_sys.  Does
    // NOT bump gem_calib_rev — noise is unchanged, and frontend uses the
    // per-event zs_sigma so band always tracks encoded hits.
    void setGemZsSigma(float v);

    // ---- Filters ---------------------------------------------------------------
    std::string loadFilter(const nlohmann::json &j);
    void unloadFilter();
    nlohmann::json filterToJson() const;
    bool filterActive() const { return trigger_type_filter.enable || waveform_filter.enable || cluster_filter.enable; }
    bool evaluateFilter(fdec::EventData &event, ssp::SspEventData *ssp) const;

    // Fill common config fields into a JSON object (used by both viewer and monitor).
    void fillConfigJson(nlohmann::json &cfg) const;

    // Handle a read-only API route. Returns {handled, response_json}.
    // Does NOT handle /api/config, clear endpoints, or mode-specific routes.
    struct ApiResult { bool handled; std::string body; };
    ApiResult handleReadApi(const std::string &uri) const;

private:
    void resolveFilterKeys();
};
