#include "app_state.h"

#include <cmath>

using json = nlohmann::json;

// Serialize a Histogram to JSON.
static json histToJson(const Histogram &h, float mn, float mx, float st)
{
    if (h.bins.empty())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"min", mn}, {"max", mx}, {"step", st}};
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"min", mn}, {"max", mx}, {"step", st}};
}

//=============================================================================
// API response builders
//=============================================================================

json AppState::apiColorRanges() const
{
    json obj = json::object();
    for (auto &[k, v] : color_range_defaults)
        obj[k] = {v.first, v.second};
    return obj;
}

json AppState::apiHist(int type, const std::string &key) const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    auto &hmap = (type == 0) ? histograms : (type == 1) ? pos_histograms : height_histograms;
    int nbins  = (type == 0) ? hist_nbins : (type == 1) ? pos_nbins : height_nbins;
    auto it = hmap.find(key);
    if (it == hmap.end())
        return {{"bins", std::vector<int>(nbins, 0)}, {"underflow", 0}, {"overflow", 0},
                {"events", events_processed.load()}};
    auto &h = it->second;
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"events", events_processed.load()}};
}

json AppState::apiClusterHist() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    json r = histToJson(cluster_energy_hist, cl_hist_min, cl_hist_max, cl_hist_step);
    r["events"] = cluster_events_processed.load();
    r["nclusters"] = histToJson(nclusters_hist,
        nclusters_hist_min, nclusters_hist_max, nclusters_hist_step);
    r["nblocks"] = histToJson(nblocks_hist,
        (float)nblocks_hist_min, (float)nblocks_hist_max, (float)nblocks_hist_step);
    r["raw_energy"] = histToJson(raw_energy_hist,
        raw_energy_hist_min, raw_energy_hist_max, raw_energy_hist_step);
    // Per-Ncl bucket dependent histograms — bins_by_ncl[i] is the
    // bins array of the i-th bucket (same indexing as nclusters_hist).
    // Frontend uses these to redraw the energy / blocks histos when the
    // user clicks a particular Ncl bar.
    json energy_by_ncl = json::array();
    json blocks_by_ncl = json::array();
    for (auto &h : cluster_energy_hist_by_ncl) energy_by_ncl.push_back(h.bins);
    for (auto &h : nblocks_hist_by_ncl)        blocks_by_ncl.push_back(h.bins);
    r["bins_by_ncl"]            = energy_by_ncl;
    r["nblocks"]["bins_by_ncl"] = blocks_by_ncl;
    return r;
}

json AppState::apiEnergyAngle() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {{"bins", energy_angle_hist.bins},
            {"nx", energy_angle_hist.nx}, {"ny", energy_angle_hist.ny},
            {"angle_min", ea_angle_min}, {"angle_max", ea_angle_max}, {"angle_step", ea_angle_step},
            {"energy_min", ea_energy_min}, {"energy_max", ea_energy_max}, {"energy_step", ea_energy_step},
            {"target", {target_x, target_y, target_z}},
            {"hycal_z", hycal_transform.z},
            {"beam_energy", beam_energy.load()},
            {"events", cluster_events_processed.load()}};
}

json AppState::apiMoller() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {{"xy_bins", moller_xy_hist.bins},
            {"xy_nx", moller_xy_hist.nx}, {"xy_ny", moller_xy_hist.ny},
            {"xy_x_min", moller_xy_x_min}, {"xy_x_max", moller_xy_x_max}, {"xy_x_step", moller_xy_x_step},
            {"xy_y_min", moller_xy_y_min}, {"xy_y_max", moller_xy_y_max}, {"xy_y_step", moller_xy_y_step},
            {"moller_events", moller_events},
            {"total_events", cluster_events_processed.load()},
            {"target",  {target_x, target_y, target_z}},
            {"hycal_z", hycal_transform.z},
            {"cuts", {{"energy_tolerance", moller_energy_tol},
                      {"angle_min", moller_angle_min}, {"angle_max", moller_angle_max}}}};
}

json AppState::apiHycalXY() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {{"xy_bins", hycal_xy_hist.bins},
            {"xy_nx", hycal_xy_hist.nx}, {"xy_ny", hycal_xy_hist.ny},
            {"xy_x_min", hxy_x_min}, {"xy_x_max", hxy_x_max}, {"xy_x_step", hxy_x_step},
            {"xy_y_min", hxy_y_min}, {"xy_y_max", hxy_y_max}, {"xy_y_step", hxy_y_step},
            {"events", hycal_xy_events},
            {"total_events", cluster_events_processed.load()},
            {"beam_energy", beam_energy.load()},
            {"cuts", {{"n_clusters", hxy_n_clusters},
                      {"energy_frac_min", hxy_energy_frac_min},
                      {"nblocks_min", hxy_nblocks_min},
                      {"nblocks_max", hxy_nblocks_max}}}};
}

json AppState::apiGemResiduals() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    json dets = json::array();
    int n = (int)gem_dx_hist.size();
    int n_dets_runtime = std::min(n, gem_sys.GetNDetectors());
    for (int d = 0; d < n; ++d) {
        std::string name = (d < n_dets_runtime)
            ? gem_sys.GetDetectors()[d].name
            : ("GEM" + std::to_string(d));
        dets.push_back({
            {"id", d},
            {"name", name},
            {"dx_hist", histToJson(gem_dx_hist[d], gem_resid_min, gem_resid_max, gem_resid_step)},
            {"dy_hist", histToJson(gem_dy_hist[d], gem_resid_min, gem_resid_max, gem_resid_step)},
            {"matched_hits", gem_match_hits[d]},
        });
    }
    return {{"enabled", gem_enabled},
            {"detectors", dets},
            {"events", gem_match_events},
            {"cuts", {{"match_nsigma", gem_match_nsigma},
                      {"require_ep_candidate", gem_match_require_ep}}}};
}

json AppState::apiGemEfficiency() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    int n = (int)gem_eff_num.size();
    int n_dets_runtime = std::min(n, gem_sys.GetNDetectors());

    auto loo_mode_name = [](GemEffLooMode m) -> const char* {
        switch (m) {
            case GemEffLooMode::Loo:           return "loo";
            case GemEffLooMode::LooTargetIn:   return "loo-target-in";
            case GemEffLooMode::TargetSeed:    return "loo-target-seed";
        }
        return "loo-target-seed";
    };

    json counters = json::array();
    json detectors = json::array();
    for (int d = 0; d < n; ++d) {
        std::string name = (d < n_dets_runtime)
            ? gem_sys.GetDetectors()[d].name
            : ("GEM" + std::to_string(d));
        int num = gem_eff_num[d];
        int den = (d < (int)gem_eff_den.size()) ? gem_eff_den[d] : 0;
        float eff_pct = (den > 0) ? (100.f * num / den) : 0.f;
        // Per-detector LOO: each card has its own denominator (anchors that
        // succeeded with detector d as the test detector).
        counters.push_back({
            {"id", d}, {"name", name},
            {"num", num}, {"den", den}, {"eff_pct", eff_pct},
        });
        json info = {{"id", d}, {"name", name}};
        if (d < n_dets_runtime) {
            const auto &det = gem_sys.GetDetectors()[d];
            info["x_size"] = det.planes[0].size;
            info["y_size"] = det.planes[1].size;
            // Active strip extent in detector-local coords (mm).  Tighter
            // than the bbox on the inner-edge side because pos=11 shares
            // strip numbers with pos=10 via shared_pos and doesn't extend
            // the readout — frontend uses these to draw the dashed frame
            // so the heatmap fits flush against it.
            auto xr = gem_sys.GetActiveExtent(d, 0);
            auto yr = gem_sys.GetActiveExtent(d, 1);
            info["x_active"] = json::array({xr.first, xr.second});
            info["y_active"] = json::array({yr.first, yr.second});
        }
        if (d < (int)gem_transforms.size()) {
            const auto &t = gem_transforms[d];
            info["position"] = json::array({t.x, t.y, t.z});
            info["tilting"]  = json::array({t.rx, t.ry, t.rz});
        }
        // Per-detector efficiency-vs-position grid: bins of predicted local
        // (x,y) at this detector when it served as the LOO test.  Frontend
        // computes per-bin eff = num/den (masking 0-den bins).
        if (d < n_dets_runtime
            && d < (int)gem_eff_grid_num.size()
            && d < (int)gem_eff_grid_den.size()) {
            const auto &gn  = gem_eff_grid_num[d];
            const auto &gd  = gem_eff_grid_den[d];
            auto xr = gem_sys.GetActiveExtent(d, 0);
            auto yr = gem_sys.GetActiveExtent(d, 1);
            info["eff_grid"] = {
                {"nx",     gn.nx},
                {"ny",     gn.ny},
                {"x_min",  xr.first},
                {"x_max",  xr.second},
                {"y_min",  yr.first},
                {"y_max",  yr.second},
                {"num",    gn.bins},
                {"den",    gd.bins},
            };
        }
        detectors.push_back(info);
    }
    json diag = json::array();
    for (int d = 0; d < 4; ++d) {
        diag.push_back({
            {"test_d",       d},
            {"n_call",       gem_eff_diag_call[d]},
            {"n_3matched",   gem_eff_diag_3matched[d]},
            {"n_pass_chi2",  gem_eff_diag_pass_chi2[d]},
            {"n_pass_resid", gem_eff_diag_pass_resid[d]},
        });
    }
    return {
        {"enabled",   gem_enabled},
        {"counters",  counters},
        {"detectors", detectors},
        {"diag",      diag},
        {"snapshot",  gemEffSnapshotJson()},
        {"hycal_z",   hycal_transform.z},
        {"target_z",  target_z},
        {"z_target_hist", histToJson(gem_eff_z_target_hist,
                                     gem_eff_z_target_min,
                                     gem_eff_z_target_max,
                                     gem_eff_z_target_step)},
        {"config", {
            {"loo_mode",              loo_mode_name(gem_eff_loo_mode)},
            {"min_cluster_energy",    gem_eff_min_cluster_energy},
            {"match_nsigma",          gem_eff_match_nsigma},
            {"max_chi2_per_dof",      gem_eff_max_chi2},
            {"max_hits_per_detector", gem_eff_max_hits_per_det},
            {"min_denom_for_eff",     gem_eff_min_denom},
            {"healthy",               gem_eff_healthy},
            {"warning",               gem_eff_warning},
            {"target_sigma",          {gem_eff_target_sigma_x,
                                       gem_eff_target_sigma_y,
                                       gem_eff_target_sigma_z}},
        }},
    };
}

json AppState::apiOccupancy() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    json jocc = json::object(), jtcut = json::object();
    for (auto &[k,v] : occupancy) jocc[k] = v;
    for (auto &[k,v] : occupancy_tcut) jtcut[k] = v;
    return {{"occ", jocc}, {"occ_tcut", jtcut}, {"total", events_processed.load()}};
}

// Reference correction: scalar factor = Alpha_latest / LMS_latest for the chosen
// ref channel.  Both readings are the most recent of their respective triggers
// (Am-241 alpha source vs. LMS pulser).  Applied uniformly to all history
// entries: corrected = signal * factor.  The Alpha source is stable, so this
// removes the LMS pulser's own drift while leaving real channel-gain changes.
struct RefCorrection {
    float factor = 1.f;
    float lms = 0.f;          // latest LMS integral on the ref module (for telemetry)
    float alpha = 0.f;        // latest Alpha integral on the ref module
    bool  active = false;
};

static RefCorrection buildRefCorrection(
    const std::map<int, float> &latest_lms,
    const std::map<int, float> &latest_alpha,
    const std::vector<AppState::LmsRefChannel> &refs, int ref_index)
{
    RefCorrection rc;
    if (ref_index < 0 || ref_index >= static_cast<int>(refs.size())) return rc;
    int ri = refs[ref_index].module_index;
    if (ri < 0) return rc;
    auto lit = latest_lms.find(ri);
    auto ait = latest_alpha.find(ri);
    if (lit == latest_lms.end() || ait == latest_alpha.end()) return rc;
    if (lit->second <= 0 || ait->second <= 0) return rc;
    rc.lms    = lit->second;
    rc.alpha  = ait->second;
    rc.factor = rc.alpha / rc.lms;
    rc.active = true;
    return rc;
}

// Apply correction: returns signal * factor (uniform across history).
static float applyRefCorrection(float val, const RefCorrection &rc)
{
    if (!rc.active) return val;
    return val * rc.factor;
}

json AppState::apiLmsSummary(int ref_index) const
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    auto rc = buildRefCorrection(latest_lms_integral, latest_alpha_integral,
                                  lms_ref_channels, ref_index);

    // ---- Gain-drift lamp scale --------------------------------------------
    // lamp_scale = current_LMS_mean[ref] / baseline_LMS_peak[ref], using the
    // ref channel named in lms_drift_ref_channel (falls back to the first
    // ref channel that has both a current history and a baseline entry).
    // This cancels the LMS pulser / FADC scale change between baseline and
    // current run so the drift ratio reflects PMT gain only.  Computed once
    // per request and shared across modules.
    bool drift_enabled = driftEnabled();
    double lamp_scale = 0.;
    std::string lamp_ref_used;
    auto compute_curr_ref_mean = [&](int mod_idx) -> double {
        auto it = lms_history.find(mod_idx);
        if (it == lms_history.end() || it->second.empty()) return 0.;
        double s = 0; int n = 0;
        for (auto &e : it->second) { s += e.integral; ++n; }
        return n > 0 ? s / n : 0.;
    };
    if (drift_enabled) {
        // Try the configured ref channel first (e.g. "LMS2"), then fall
        // back to whichever lms_ref_channels entry has a usable pair.
        auto try_ref = [&](const std::string &name) {
            if (lamp_scale > 0) return;
            auto bit = lms_baseline_ref_peak.find(name);
            if (bit == lms_baseline_ref_peak.end() || bit->second <= 0) return;
            for (auto &rc : lms_ref_channels) {
                if (rc.name != name || rc.module_index < 0) continue;
                double curr = compute_curr_ref_mean(rc.module_index);
                if (curr > 0) {
                    lamp_scale = curr / bit->second;
                    lamp_ref_used = name;
                }
                break;
            }
        };
        if (!lms_drift_ref_channel.empty()) try_ref(lms_drift_ref_channel);
        if (lamp_scale <= 0)
            for (auto &rc : lms_ref_channels) try_ref(rc.name);
    }

    json mods = json::object();
    for (auto &[idx, hist] : lms_history) {
        if (hist.empty()) continue;
        double sum = 0, sum2 = 0;
        int count = 0;
        for (auto &e : hist) {
            float v = applyRefCorrection(e.integral, rc);
            sum += v; sum2 += v * v;
            count++;
        }
        if (count == 0) continue;
        double mean = sum / count;
        double var = sum2 / count - mean * mean;
        double rms = var > 0 ? std::sqrt(var) : 0;

        bool warn = (mean > 0 && rms / mean > lms_warn_thresh) ||
                    (mean < lms_warn_min_mean);

        // ---- drift-from-baseline (gain monitor) ----
        // drift = current mean / (baseline_lms_peak * lamp_scale)
        // Use the UNCORRECTED current mean (independent of the ref-correction
        // toggle) so the drift state never changes when the user flips the
        // LMS-tab Ref dropdown for visual normalization.
        double raw_sum = 0;
        for (auto &e : hist) raw_sum += e.integral;
        double raw_mean = raw_sum / count;
        double drift_val = std::nan("");
        bool   drift_flag = false;
        bool   drift_suppressed = false;
        if (drift_enabled && lamp_scale > 0
            && idx >= 0 && idx < hycal.module_count())
        {
            auto &mod = hycal.module(idx);
            auto bit = lms_baseline_peak.find(mod.name);
            // Skip channels whose baseline LMS peak is implausibly low —
            // gain_fitter occasionally returns a noise-floor fit (~few ADC)
            // for dead/saturated modules, which would inflate raw_mean/expected
            // into a phantom huge drift.  30 ADC is well below any healthy
            // LMS peak (typically 200+) but above noise.
            if (bit != lms_baseline_peak.end() && bit->second > 30.f) {
                double expected = bit->second * lamp_scale;
                // Same alive-channel gate on the post-scale expected value:
                // if lamp_scale collapsed the baseline below the warn floor,
                // the channel is effectively dead this run and would warn
                // anyway — don't double-flag it.
                if (expected > lms_warn_min_mean) {
                    drift_val = raw_mean / expected;
                    // Pick threshold pair by module type.  W/G have different
                    // healthy widths so they get separate bounds; everything
                    // else (V, LMS refs themselves, etc.) uses the W bounds
                    // as a sensible default.
                    bool is_glass = mod.is_glass();
                    float lo = is_glass ? lms_drift_low_g  : lms_drift_low_w;
                    float hi = is_glass ? lms_drift_high_g : lms_drift_high_w;
                    if (drift_val < lo || drift_val > hi)
                        drift_flag = true;
                    // Suppress the WARN (not the value) for whole module
                    // types listed in lms_drift_suppress_types.  Operators
                    // can still see the drift number in the table column;
                    // it just doesn't escalate state to "drift" or move
                    // the row to the top.
                    if (drift_flag &&
                        lms_drift_suppress_types.count(static_cast<int>(mod.type)))
                    {
                        drift_flag = false;
                        drift_suppressed = true;
                    }
                }
            }
        }

        // Tri-state, top-priority first: drift > warn > ok.  The single
        // 'state' string is what the GUI / report sort and color on; the
        // legacy 'warn' bool stays true for warn OR drift so any older UI
        // stops on a problem.
        const char *state = drift_flag ? "drift"
                          : (warn ? "warn" : "ok");

        if (idx >= 0 && idx < hycal.module_count()) {
            auto &mod = hycal.module(idx);
            json entry = {
                {"name", mod.name}, {"mean", std::round(mean * 10) / 10},
                {"rms",  std::round(rms  * 100) / 100},
                {"count", count}, {"warn", warn || drift_flag},
                {"state", state}};
            if (!std::isnan(drift_val))
                entry["drift"] = std::round(drift_val * 1000) / 1000;
            else
                entry["drift"] = nullptr;
            if (drift_suppressed) entry["drift_suppressed"] = true;
            mods[std::to_string(idx)] = std::move(entry);
        }
    }
    return {{"modules", mods}, {"events", lms_events.load()},
            {"trigger", lms_trigger.toJson()},
            {"ref_index", ref_index},
            {"ref_factor", rc.factor},
            {"ref_lms", rc.lms},
            {"ref_alpha", rc.alpha},
            {"drift_enabled", drift_enabled},
            {"drift_lamp_scale", lamp_scale > 0 ? json(lamp_scale) : json(nullptr)},
            {"drift_lamp_ref",   lamp_ref_used},
            {"drift_low_w",      lms_drift_low_w},
            {"drift_high_w",     lms_drift_high_w},
            {"drift_low_g",      lms_drift_low_g},
            {"drift_high_g",     lms_drift_high_g},
            {"drift_suppress_types", lms_drift_suppress_type_names},
            {"sync_unix", sync_unix}, {"sync_rel_sec", sync_rel_sec}};
}

json AppState::apiLmsModule(int mod_idx, int ref_index) const
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    auto it = lms_history.find(mod_idx);
    if (it == lms_history.end() || it->second.empty())
        return {{"time", json::array()}, {"integral", json::array()}, {"events", 0}};

    auto rc = buildRefCorrection(latest_lms_integral, latest_alpha_integral,
                                  lms_ref_channels, ref_index);

    auto &hist = it->second;
    json t_arr = json::array(), v_arr = json::array();
    for (auto &e : hist) {
        float v = applyRefCorrection(e.integral, rc);
        t_arr.push_back(std::round(e.time_sec * 100) / 100);
        v_arr.push_back(std::round(v * 10) / 10);
    }
    std::string name = (mod_idx >= 0 && mod_idx < hycal.module_count())
        ? hycal.module(mod_idx).name : "";
    return {{"name", name}, {"time", t_arr}, {"integral", v_arr},
            {"events", (int)t_arr.size()},
            {"ref_index", ref_index},
            {"ref_factor", rc.factor},
            {"ref_lms", rc.lms},
            {"ref_alpha", rc.alpha},
            {"sync_unix", sync_unix}, {"sync_rel_sec", sync_rel_sec}};
}

json AppState::apiLmsRefChannels() const
{
    json arr = json::array();
    for (size_t i = 0; i < lms_ref_channels.size(); ++i) {
        arr.push_back({
            {"index", (int)i},
            {"name", lms_ref_channels[i].name},
            {"module_index", lms_ref_channels[i].module_index},
        });
    }
    return arr;
}

//=============================================================================
// EPICS
//=============================================================================

void AppState::processEpics(const std::string &text, int32_t event_number, uint64_t timestamp)
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    epics.Feed(event_number, timestamp, text);
    epics.Trim(epics_max_history);
    epics_events++;

    // Single-source beam energy: latest valid MBSY2C_energy reading overrides the
    // runinfo fallback. Skip values below min_valid (zero/garbage during beam trips).
    if (!beam_energy_epics_channel.empty()) {
        int id = epics.GetChannelId(beam_energy_epics_channel);
        if (id >= 0) {
            int n = epics.GetSnapshotCount();
            if (n > 0) {
                const auto &snap = epics.GetSnapshot(n - 1);
                if (id < (int)snap.values.size()) {
                    float v = snap.values[id];
                    if (v > beam_energy_min_valid) beam_energy.store(v);
                }
            }
        }
    }
}

void AppState::clearEpics()
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    epics.Clear();
    epics_events = 0;
}

// ---------- DSC2 scaler bank → measured livetime --------------------------
// Bank parsing + (source, channel) selection live in dsc::Dsc2Decoder
// (see prad2dec/include/Dsc2Decoder.h).  This method only adds the
// AppState-side concern: convert two consecutive cumulative readings into a
// realtime live-time fraction and broadcast it.  Convention: gated counts
// LIVE time in this DAQ (Group A enabled while NOT busy), so
// live = gated / ungated.  Until a second read arrives we display the
// cumulative ratio.
void AppState::processDsc(const dsc::DscEventData &dsc)
{
    if (!dsc.present || dsc.ungated == 0) return;

    const uint32_t prev_g = dsc_prev_gated.load();
    const uint32_t prev_u = dsc_prev_ungated.load();
    const double lt = (dsc.gated >= prev_g && dsc.ungated > prev_u)
        ? (double)(dsc.gated   - prev_g) /
          (double)(dsc.ungated - prev_u) * 100.0
        : (double)dsc.gated / (double)dsc.ungated * 100.0;
    measured_livetime.store(lt);
    dsc_prev_gated.store(dsc.gated);
    dsc_prev_ungated.store(dsc.ungated);
}

json AppState::apiEpicsChannels() const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    json names = json::array();
    for (auto &n : epics.GetChannelNames()) names.push_back(n);
    json slots = json::array();
    for (auto &s : epics_default_slots) slots.push_back(s);
    return {{"channels", names}, {"slots", slots},
            {"events", epics_events.load()}};
}

// Time anchor for EPICS snapshots.  We pick the EARLIEST snapshot that
// actually carries a TI tick (timestamp != 0) — snapshots with timestamp == 0
// are EPICS events that arrived before any physics event (typical at the
// start of a run, or right after an ET reconnect when last_ti_ts is reset
// to 0 and EPICS comes through before any new physics event).  Using the
// first non-zero snapshot guarantees ti_delta_sec never has to special-case
// "anchor came later than data"; combined with ti_delta_sec's own guards
// against now==0 and now<base, this is the fix for the famous
// 73,786,976,288 s underflow displayed in the EPICS monitor.
static uint64_t epics_anchor_ts(const epics::EpicsStore &store)
{
    int n = store.GetSnapshotCount();
    for (int i = 0; i < n; ++i) {
        uint64_t t = store.GetSnapshot(i).timestamp;
        if (t != 0) return t;
    }
    return 0;
}

json AppState::apiEpicsChannel(const std::string &name) const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    int id = epics.GetChannelId(name);
    if (id < 0)
        return {{"name", name}, {"time", json::array()}, {"value", json::array()}, {"count", 0}};

    int nsnap = epics.GetSnapshotCount();
    json t_arr = json::array(), v_arr = json::array();

    uint64_t t0 = epics_anchor_ts(epics);
    for (int i = 0; i < nsnap; ++i) {
        auto &snap = epics.GetSnapshot(i);
        double t_sec = ti_delta_sec(snap.timestamp, t0);
        float val = (id < (int)snap.values.size()) ? snap.values[id] : 0.f;
        t_arr.push_back(std::round(t_sec * 100) / 100);
        v_arr.push_back(val);
    }
    return {{"name", name}, {"time", t_arr}, {"value", v_arr}, {"count", nsnap}};
}

json AppState::apiEpicsBatch(const std::vector<std::string> &names) const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    int nsnap = epics.GetSnapshotCount();
    uint64_t t0 = epics_anchor_ts(epics);

    // build shared time array once
    json t_arr = json::array();
    for (int i = 0; i < nsnap; ++i) {
        double t_sec = ti_delta_sec(epics.GetSnapshot(i).timestamp, t0);
        t_arr.push_back(std::round(t_sec * 100) / 100);
    }

    json channels = json::array();
    for (auto &name : names) {
        int id = epics.GetChannelId(name);
        if (id < 0) {
            channels.push_back({{"name", name}, {"value", json::array()}, {"count", 0}});
            continue;
        }
        json v_arr = json::array();
        for (int i = 0; i < nsnap; ++i) {
            auto &snap = epics.GetSnapshot(i);
            v_arr.push_back((id < (int)snap.values.size()) ? snap.values[id] : 0.f);
        }
        channels.push_back({{"name", name}, {"value", v_arr}, {"count", nsnap}});
    }
    return {{"time", t_arr}, {"channels", channels}};
}

json AppState::apiEpicsLatest() const
{
    std::lock_guard<std::mutex> lk(epics_mtx);
    json channels = json::array();
    int nsnap = epics.GetSnapshotCount();
    int nch = epics.GetChannelCount();
    if (nsnap == 0 || nch == 0)
        return {{"channels", channels}, {"events", epics_events.load()}};

    auto &latest = epics.GetSnapshot(nsnap - 1);

    // compute per-channel mean from most recent mean_window snapshots
    int win_start = std::max(0, nsnap - epics_mean_window);
    std::vector<double> sums(nch, 0.0);
    std::vector<int> counts(nch, 0);
    for (int i = win_start; i < nsnap; ++i) {
        auto &snap = epics.GetSnapshot(i);
        for (int ch = 0; ch < std::min(nch, (int)snap.values.size()); ++ch) {
            sums[ch] += snap.values[ch];
            counts[ch]++;
        }
    }

    for (int ch = 0; ch < nch; ++ch) {
        float val = (ch < (int)latest.values.size()) ? latest.values[ch] : 0.f;
        float mean = (counts[ch] > 0) ? static_cast<float>(sums[ch] / counts[ch]) : val;
        channels.push_back({
            {"name", epics.GetChannelName(ch)},
            {"value", std::round(val * 1000) / 1000},
            {"mean", std::round(mean * 1000) / 1000},
            {"count", counts[ch]},
        });
    }
    return {{"channels", channels}, {"events", epics_events.load()}};
}

//=============================================================================
// Shared config + API routing (used by both viewer and monitor)
//=============================================================================

void AppState::fillConfigJson(json &cfg) const
{
    cfg["hist"] = {
        {"bin_min", hist_cfg.bin_min}, {"bin_max", hist_cfg.bin_max},
        {"bin_step", hist_cfg.bin_step},
        {"pos_min", hist_cfg.pos_min}, {"pos_max", hist_cfg.pos_max},
        {"pos_step", hist_cfg.pos_step},
        {"height_min", hist_cfg.height_min}, {"height_max", hist_cfg.height_max},
        {"height_step", hist_cfg.height_step},
    };
    cfg["waveform_filter"]         = peak_filter.toJson(peak_quality_bits_def);
    cfg["waveform_filter_active"]  = peak_filter.enable;
    cfg["waveform_filter_default"] = peak_filter_default.toJson(peak_quality_bits_def);
    cfg["quality_bits"]            = peak_quality_bits_def;
    cfg["ref_lines"] = ref_lines;
    cfg["trigger_bits"] = trigger_bits_def;
    cfg["trigger_type"] = trigger_type_def;
    cfg["trigger_filter"] = {
        {"dq",      waveform_trigger.toJson()},
        {"cluster", cluster_trigger.toJson()},
        {"lms",     lms_trigger.toJson()},
        {"physics", physics_trigger.toJson()},
    };
    cfg["cluster_hist"] = {{"min", cl_hist_min}, {"max", cl_hist_max}, {"step", cl_hist_step}};
    cfg["nclusters_hist"] = {{"min", nclusters_hist_min}, {"max", nclusters_hist_max}, {"step", nclusters_hist_step}};
    cfg["nblocks_hist"] = {{"min", nblocks_hist_min}, {"max", nblocks_hist_max}, {"step", nblocks_hist_step}};
    cfg["raw_energy_hist"] = {{"min", raw_energy_hist_min}, {"max", raw_energy_hist_max}, {"step", raw_energy_hist_step}};
    cfg["color_ranges"] = apiColorRanges();
    cfg["refresh_ms"] = {{"event", refresh_event_ms}, {"ring", refresh_ring_ms},
                         {"histogram", refresh_hist_ms}, {"lms", refresh_lms_ms}};
    cfg["lms"] = {
        {"trigger", lms_trigger.toJson()},
        {"warn_threshold", lms_warn_thresh},
        {"events", lms_events.load()}, {"ref_channels", apiLmsRefChannels()},
        {"drift_enabled",  driftEnabled()},
        {"drift_baseline", lms_drift_baseline_file},
        {"drift_ref",      lms_drift_ref_channel},
        {"drift_low_w",    lms_drift_low_w},
        {"drift_high_w",   lms_drift_high_w},
        {"drift_low_g",    lms_drift_low_g},
        {"drift_high_g",   lms_drift_high_g},
        {"drift_suppress_types", lms_drift_suppress_type_names},
    };
    auto metric_cfg = [](const ShellMetric &m) {
        nlohmann::json j = {
            {"enabled",  !m.command.empty()},
            {"unit",     m.unit},
            {"poll_sec", m.poll_sec},
        };
        if (m.has_trip_warn) j["trip_warn_below"] = m.trip_warn_below;
        return j;
    };
    cfg["monitor_status"] = {
        {"livetime", {
            {"enabled",          !livetime_cmd.empty()},
            {"measured_enabled", daq_cfg.dsc_scaler.enabled()},
            {"unit",             livetime_unit},
            {"poll_sec",         livetime_poll_sec},
            {"healthy",          livetime_healthy},
            {"warning",          livetime_warning},
        }},
        {"beam", {
            {"energy",  metric_cfg(beam_energy_status)},
            {"current", metric_cfg(beam_current_status)},
        }},
    };
    cfg["runinfo"] = {
        {"beam_energy", beam_energy.load()},
        {"beam_energy_runinfo", beam_energy_runinfo},
        {"calibration", {{"default_adc2mev", adc_to_mev}}},
        {"target", {target_x, target_y, target_z}},
        {"hycal", {
            {"position", {hycal_transform.x, hycal_transform.y, hycal_transform.z}},
            {"tilting", {hycal_transform.rx, hycal_transform.ry, hycal_transform.rz}},
        }},
    };
    cfg["physics"] = {
        {"trigger", physics_trigger.toJson()},
        {"beam_energy", {
            {"epics_channel", beam_energy_epics_channel},
            {"min_valid", beam_energy_min_valid},
        }},
        {"energy_angle_hist", {
            {"angle_min", ea_angle_min}, {"angle_max", ea_angle_max}, {"angle_step", ea_angle_step},
            {"energy_min", ea_energy_min}, {"energy_max", ea_energy_max}, {"energy_step", ea_energy_step},
        }},
        {"moller", {
            {"energy_tolerance", moller_energy_tol},
            {"angle_min", moller_angle_min}, {"angle_max", moller_angle_max},
        }},
        {"hycal_cluster_hit", {
            {"n_clusters", hxy_n_clusters},
            {"energy_frac_min", hxy_energy_frac_min},
            {"nblocks_min", hxy_nblocks_min},
            {"nblocks_max", hxy_nblocks_max},
        }},
    };
    cfg["auto_report"] = {
        {"enabled",                  auto_report_enabled},
        {"post_to_elog",             auto_report_post_to_elog},
        {"local_save_dir",           auto_report_local_save_dir},
        {"min_interval_ms",          auto_report_min_interval_ms},
        {"schedule_minutes",         auto_report_schedule_minutes},
        // Surfaced to the client so report.js can flag low/partial-data
        // reports inline; report.js falls back to its own defaults if
        // missing, so the client and server stay backward-compatible.
        {"min_events_for_schedule",  auto_report_min_events_for_schedule},
        {"schedule_max_wait_min",    auto_report_schedule_max_wait_min},
        {"partial_threshold_events", auto_report_partial_threshold_events},
        {"elog", {
            {"url", elog_url}, {"logbook", elog_logbook},
            {"author", elog_author}, {"tags", elog_tags},
        }},
    };
    cfg["epics"] = {
        {"max_history", epics_max_history},
        {"warn_threshold", epics_warn_thresh}, {"alert_threshold", epics_alert_thresh},
        {"min_avg_points", epics_min_avg_pts}, {"mean_window", epics_mean_window},
        {"slots", epics_default_slots},
    };
    // GEM tab owns its configuration: detector geometry (from apiGemConfig)
    // plus the diagnostic configs that used to live under physics.
    cfg["gem"] = apiGemConfig();
    cfg["gem"]["hycal_match"] = {
        {"require_ep_candidate", gem_match_require_ep},
        {"match_nsigma",         gem_match_nsigma},
        {"residual_hist", {
            {"min", gem_resid_min}, {"max", gem_resid_max}, {"step", gem_resid_step},
        }},
    };
    {
        const char *loo_name = "loo-target-seed";
        switch (gem_eff_loo_mode) {
            case GemEffLooMode::Loo:           loo_name = "loo"; break;
            case GemEffLooMode::LooTargetIn:   loo_name = "loo-target-in"; break;
            case GemEffLooMode::TargetSeed:    loo_name = "loo-target-seed"; break;
        }
        cfg["gem"]["efficiency"] = {
            {"loo_mode",              loo_name},
            {"min_cluster_energy",    gem_eff_min_cluster_energy},
            {"match_nsigma",          gem_eff_match_nsigma},
            {"max_chi2_per_dof",      gem_eff_max_chi2},
            {"max_hits_per_detector", gem_eff_max_hits_per_det},
            {"min_denom_for_eff",     gem_eff_min_denom},
            {"healthy",               gem_eff_healthy},
            {"warning",               gem_eff_warning},
            {"target_sigma",          json::array({gem_eff_target_sigma_x,
                                                   gem_eff_target_sigma_y,
                                                   gem_eff_target_sigma_z})},
        };
    }
    cfg["gem"]["pos_res"] = gem_pos_res;
    cfg["gem"]["hycal_pos_res"] = json::array({
        hycal.GetPositionResolutionA(),
        hycal.GetPositionResolutionB(),
        hycal.GetPositionResolutionC()
    });
}

AppState::ApiResult AppState::handleReadApi(const std::string &uri) const
{
    if (uri == "/api/occupancy")
        return {true, apiOccupancy().dump()};
    if (uri == "/api/physics/energy_angle")
        return {true, apiEnergyAngle().dump()};
    if (uri == "/api/physics/moller")
        return {true, apiMoller().dump()};
    if (uri == "/api/physics/hycal_xy")
        return {true, apiHycalXY().dump()};
    if (uri == "/api/gem/residuals")
        return {true, apiGemResiduals().dump()};
    if (uri == "/api/gem/efficiency")
        return {true, apiGemEfficiency().dump()};
    if (uri == "/api/cluster_hist")
        return {true, apiClusterHist().dump()};
    if (uri.rfind("/api/hist/", 0) == 0)
        return {true, apiHist(0, uri.substr(10)).dump()};
    if (uri.rfind("/api/poshist/", 0) == 0)
        return {true, apiHist(1, uri.substr(13)).dump()};
    if (uri.rfind("/api/heighthist/", 0) == 0)
        return {true, apiHist(2, uri.substr(16)).dump()};
    if (uri == "/api/lms/refs")
        return {true, apiLmsRefChannels().dump()};
    if (uri.rfind("/api/lms/", 0) == 0) {
        int ref = -1;
        auto qpos = uri.find('?');
        std::string path = (qpos != std::string::npos) ? uri.substr(9, qpos - 9) : uri.substr(9);
        if (qpos != std::string::npos) {
            std::string q = uri.substr(qpos + 1);
            if (q.rfind("ref=", 0) == 0) ref = std::atoi(q.c_str() + 4);
        }
        if (path == "summary") return {true, apiLmsSummary(ref).dump()};
        if (path == "clear")   return {false, ""};  // clear handled by caller
        return {true, apiLmsModule(std::atoi(path.c_str()), ref).dump()};
    }
    if (uri.rfind("/api/epics/", 0) == 0) {
        std::string path = uri.substr(11);
        if (path == "channels") return {true, apiEpicsChannels().dump()};
        if (path == "latest")   return {true, apiEpicsLatest().dump()};
        if (path == "clear")    return {false, ""};  // clear handled by caller
        if (path.rfind("batch?", 0) == 0) {
            // /api/epics/batch?ch=name1&ch=name2&...
            std::string query = path.substr(6);
            std::vector<std::string> names;
            for (size_t pos = 0; pos < query.size();) {
                size_t amp = query.find('&', pos);
                if (amp == std::string::npos) amp = query.size();
                std::string kv = query.substr(pos, amp - pos);
                if (kv.rfind("ch=", 0) == 0) {
                    // URL-decode
                    std::string raw = kv.substr(3), name;
                    for (size_t i = 0; i < raw.size(); ++i) {
                        if (raw[i] == '%' && i + 2 < raw.size()) {
                            int hi = 0, lo = 0;
                            if (std::sscanf(raw.c_str() + i + 1, "%1x%1x", &hi, &lo) == 2) {
                                name += static_cast<char>((hi << 4) | lo);
                                i += 2; continue;
                            }
                        }
                        if (raw[i] == '+') name += ' ';
                        else name += raw[i];
                    }
                    names.push_back(name);
                }
                pos = amp + 1;
            }
            return {true, apiEpicsBatch(names).dump()};
        }
        if (path.rfind("channel/", 0) == 0) {
            // URL-decode the channel name (e.g. %3A → :)
            std::string raw = path.substr(8), name;
            for (size_t i = 0; i < raw.size(); ++i) {
                if (raw[i] == '%' && i + 2 < raw.size()) {
                    int hi = 0, lo = 0;
                    if (std::sscanf(raw.c_str() + i + 1, "%1x%1x", &hi, &lo) == 2) {
                        name += static_cast<char>((hi << 4) | lo);
                        i += 2;
                        continue;
                    }
                }
                name += raw[i];
            }
            return {true, apiEpicsChannel(name).dump()};
        }
    }
    if (uri == "/api/gem/hits")
        return {true, apiGemHits().dump()};
    if (uri == "/api/gem/config")
        return {true, apiGemConfig().dump()};
    if (uri == "/api/gem/occupancy")
        return {true, apiGemOccupancy().dump()};
    return {false, ""};
}

