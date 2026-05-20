#include "app_state.h"
#include "data_source.h"
#include "load_daq_config.h"
#include "Fadc250FwAnalyzer.h"

#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdlib>

using json = nlohmann::json;

//=============================================================================
// PeakFilter — JSON parse / serialize
//=============================================================================
namespace {

// Set out to filter[axis][bound] iff that field exists and is numeric.
void setOpt(const json &range, const char *bound, std::optional<float> &out)
{
    if (range.is_object() && range.contains(bound) && range[bound].is_number())
        out = range[bound].get<float>();
}

void parseAxis(const json &filter, const char *axis,
               std::optional<float> &mn, std::optional<float> &mx)
{
    mn.reset(); mx.reset();
    if (!filter.is_object() || !filter.contains(axis)) return;
    const auto &r = filter[axis];
    setOpt(r, "min", mn);
    setOpt(r, "max", mx);
}

json axisToJson(const std::optional<float> &mn, const std::optional<float> &mx)
{
    json o = json::object();
    if (mn) o["min"] = *mn;
    if (mx) o["max"] = *mx;
    return o;
}

} // namespace

void PeakFilter::parse(const json &filter, const json &quality_bits_def)
{
    parseAxis(filter, "time",     time_min,     time_max);
    parseAxis(filter, "integral", integral_min, integral_max);
    parseAxis(filter, "height",   height_min,   height_max);
    q_accept = q_reject = 0;
    if (filter.is_object() && filter.contains("quality_bits")
        && filter["quality_bits"].is_object()) {
        const auto &qb = filter["quality_bits"];
        if (qb.contains("accept"))
            q_accept = bitsMaskFromArray(qb["accept"], quality_bits_def, "PeakFilter");
        if (qb.contains("reject"))
            q_reject = bitsMaskFromArray(qb["reject"], quality_bits_def, "PeakFilter");
    }
}

json PeakFilter::toJson(const json &quality_bits_def) const
{
    json out = json::object();
    json t = axisToJson(time_min,     time_max);     if (!t.empty()) out["time"]     = t;
    json i = axisToJson(integral_min, integral_max); if (!i.empty()) out["integral"] = i;
    json h = axisToJson(height_min,   height_max);   if (!h.empty()) out["height"]   = h;
    if (q_accept || q_reject) {
        out["quality_bits"] = {
            {"accept", bitsMaskToNames(q_accept, quality_bits_def)},
            {"reject", bitsMaskToNames(q_reject, quality_bits_def)}
        };
    }
    return out;
}

static json histToJson(const Histogram &h, float mn, float mx, float st)
{
    if (h.bins.empty())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"min", mn}, {"max", mx}, {"step", st}};
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"min", mn}, {"max", mx}, {"step", st}};
}

//=============================================================================
// GEM efficiency monitor — internal helpers
//=============================================================================
namespace {

struct Line3D {
    float ax = 0.f, bx = 0.f;        // x(z) = ax + bx·z
    float ay = 0.f, by = 0.f;        // y(z) = ay + by·z
    float chi2_per_dof = 0.f;
};

// Two-point seed line in lab frame.  Caller guarantees z1 != z2 (HyCal vs GEM).
Line3D seedLine(float x1, float y1, float z1, float x2, float y2, float z2)
{
    Line3D L{};
    float dz = z2 - z1;
    if (std::abs(dz) < 1e-6f) { L.ax = x1; L.ay = y1; return L; }
    L.bx = (x2 - x1) / dz;  L.ax = x1 - L.bx * z1;
    L.by = (y2 - y1) / dz;  L.ay = y1 - L.by * z1;
    return L;
}

// Independent weighted LSQ fits in (z, x) and (z, y) — 4-parameter line.
// `wy = nullptr` reuses `wx` for both axes (the common case where σ_x = σ_y);
// pass distinct arrays for anisotropic per-point uncertainties (e.g. the
// target point in loo-target-in, where σ_target_z couples differently into
// σ_x_eff and σ_y_eff via the slope).  dof = 2N - 4.
bool fitWeightedLine(int N,
                     const float *z, const float *x, const float *y,
                     const float *wx, const float *wy,
                     Line3D &out)
{
    if (N < 2) return false;
    if (wy == nullptr) wy = wx;
    // x-fit
    double Swx=0, Szx=0, Szzx=0, Sx=0, Sxz=0;
    for (int i = 0; i < N; ++i) {
        double wi = wx[i];
        Swx  += wi;
        Szx  += wi * z[i];
        Szzx += wi * z[i] * z[i];
        Sx   += wi * x[i];
        Sxz  += wi * x[i] * z[i];
    }
    double Dx = Swx * Szzx - Szx * Szx;
    if (std::abs(Dx) < 1e-9) return false;
    double bx = (Swx * Sxz - Szx * Sx) / Dx;
    double ax = (Sx - bx * Szx) / Swx;
    // y-fit
    double Swy=0, Szy=0, Szzy=0, Sy=0, Syz=0;
    for (int i = 0; i < N; ++i) {
        double wi = wy[i];
        Swy  += wi;
        Szy  += wi * z[i];
        Szzy += wi * z[i] * z[i];
        Sy   += wi * y[i];
        Syz  += wi * y[i] * z[i];
    }
    double Dy = Swy * Szzy - Szy * Szy;
    if (std::abs(Dy) < 1e-9) return false;
    double by = (Swy * Syz - Szy * Sy) / Dy;
    double ay = (Sy - by * Szy) / Swy;
    out.ax = (float)ax; out.bx = (float)bx;
    out.ay = (float)ay; out.by = (float)by;
    int dof = 2 * N - 4;
    if (dof > 0) {
        double chi2 = 0;
        for (int i = 0; i < N; ++i) {
            double dxp = (ax + bx * z[i]) - x[i];
            double dyp = (ay + by * z[i]) - y[i];
            chi2 += wx[i] * dxp * dxp + wy[i] * dyp * dyp;
        }
        out.chi2_per_dof = (float)(chi2 / dof);
    } else {
        out.chi2_per_dof = 0.f;
    }
    return true;
}

// Convenience overload: same weight for x and y at every point.
inline bool fitWeightedLine(int N,
                            const float *z, const float *x, const float *y,
                            const float *w, Line3D &out)
{
    return fitWeightedLine(N, z, x, y, w, nullptr, out);
}

// Project a lab-frame line onto a detector's local plane (z_local = 0) using
// the labToLocal-then-1D-interpolate trick.
void projectLineToLocal(const DetectorTransform &xform, const Line3D &L,
                        float &px, float &py)
{
    float ax1 = L.ax,                 ay1 = L.ay,                 z1 = 0.f;
    float ax2 = L.ax + L.bx * 1000.f, ay2 = L.ay + L.by * 1000.f, z2 = 1000.f;
    float l1x, l1y, l1z, l2x, l2y, l2z;
    xform.labToLocal(ax1, ay1, z1, l1x, l1y, l1z);
    xform.labToLocal(ax2, ay2, z2, l2x, l2y, l2z);
    float dz = l2z - l1z;
    if (std::abs(dz) < 1e-6f) { px = l1x; py = l1y; return; }
    float s = -l1z / dz;
    px = l1x + s * (l2x - l1x);
    py = l1y + s * (l2y - l1y);
}

}  // anonymous namespace

//=============================================================================
// Per-event processing
//=============================================================================

// Encode peak array for one channel.
static json encodePeaks(const fdec::WaveResult &wres)
{
    json parr = json::array();
    for (int p = 0; p < wres.npeaks; ++p) {
        auto &pk = wres.peaks[p];
        parr.push_back({
            {"p", pk.pos}, {"t", std::round(pk.time * 10) / 10},
            {"h", std::round(pk.height * 10) / 10},
            {"i", std::round(pk.integral * 10) / 10},
            {"l", pk.left}, {"r", pk.right},
            {"o", pk.overflow ? 1 : 0},
            {"q", pk.quality},
        });
    }
    return parr;
}

// Run the firmware-faithful FADC250 emulator on one channel and serialize the
// result.  Used by both the on-demand /api/waveform path (file mode) and the
// ring-buffer encoder (online mode) so the DAQ overlay works in both modes.
static json encodeChannelDaq(const fdec::ChannelData &cd, float ped_mean,
                             const evc::DaqConfig::Fadc250FwConfig &fw_cfg)
{
    fdec::Fadc250FwAnalyzer fw_ana(fw_cfg);
    fdec::DaqWaveResult daq_res;
    fw_ana.Analyze(cd.samples, cd.nsamples, ped_mean, daq_res);

    json daq_pulses = json::array();
    for (int p = 0; p < daq_res.npeaks; ++p) {
        const auto &pk = daq_res.peaks[p];
        daq_pulses.push_back({
            {"n",       pk.pulse_id},
            {"vmin",    std::round(pk.vmin  * 10) / 10},
            {"vp",      std::round(pk.vpeak * 10) / 10},
            {"va",      std::round(pk.va    * 10) / 10},
            {"coarse",  pk.coarse},
            {"fine",    pk.fine},
            {"t",       std::round(pk.time_ns * 100) / 100},
            {"cross",   pk.cross_sample},
            {"vp_pos",  pk.peak_sample},
            {"i",       std::round(pk.integral * 10) / 10},
            {"wlo",     pk.window_lo},
            {"whi",     pk.window_hi},
            {"q",       pk.quality},
        });
    }
    return {
        {"vnoise",     std::round(daq_res.vnoise * 10) / 10},
        {"ped_used",   std::round(ped_mean       * 10) / 10},
        {"tet",        fw_cfg.TET},
        {"nsb",        fw_cfg.NSB},
        {"nsa",        fw_cfg.NSA},
        {"max_pulses", fw_cfg.MAX_PULSES},
        {"nsat",       fw_cfg.NSAT},
        {"nped",       fw_cfg.NPED},
        {"maxped",     fw_cfg.MAXPED},
        {"clk_ns",     fw_cfg.CLK_NS},
        {"pk",         daq_pulses},
    };
}

json AppState::encodeEventJson(fdec::EventData &event, int ev_id,
                               fdec::WaveAnalyzer &ana, fdec::WaveResult &wres,
                               bool include_samples)
{
    json channels = json::object();
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.SetChannelKey(roc.tag, s, c);
                ana.Analyze(cd.samples, cd.nsamples, wres);
                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                json ch_j = {
                    {"pm", std::round(wres.ped.mean * 10) / 10},
                    {"pr", std::round(wres.ped.rms * 10) / 10},
                    {"pk", encodePeaks(wres)},
                };
                if (include_samples) {
                    json sarr = json::array();
                    for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);
                    ch_j["s"] = std::move(sarr);
                    // Embed firmware-emulator output alongside the samples so
                    // the waveform tab's DAQ overlay works for ring-buffer
                    // events (online mode) — /api/waveform is file-mode only.
                    ch_j["daq"] = encodeChannelDaq(cd, wres.ped.mean,
                                                   daq_cfg.fadc250_fw);
                }
                channels[key] = std::move(ch_j);
            }
        }
    }
    return {{"event", ev_id}, {"channels", channels},
            {"event_number", event.info.event_number},
            {"trigger_type", event.info.trigger_type},
            {"trigger_bits", event.info.trigger_bits},
            {"run_number", event.info.run_number}};
}

json AppState::encodeWaveformJson(fdec::EventData &event, const std::string &chan_key,
                                  fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // parse "roc_slot_ch" key
    int roc_tag = 0, sl = 0, ch = 0;
    if (std::sscanf(chan_key.c_str(), "%d_%d_%d", &roc_tag, &sl, &ch) != 3)
        return {{"error", "invalid channel key"}};

    // find the channel in the event
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present || roc.tag != roc_tag) continue;
        if (!roc.slots[sl].present) break;
        if (!(roc.slots[sl].channel_mask & (1ull << ch))) break;
        auto &cd = roc.slots[sl].channels[ch];
        if (cd.nsamples <= 0) break;

        ana.SetChannelKey(roc_tag, sl, ch);
        ana.Analyze(cd.samples, cd.nsamples, wres);

        json sarr = json::array();
        for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);

        return {{"key", chan_key}, {"s", sarr},
                {"pm", std::round(wres.ped.mean * 10) / 10},
                {"pr", std::round(wres.ped.rms * 10) / 10},
                {"pk", encodePeaks(wres)},
                {"daq", encodeChannelDaq(cd, wres.ped.mean,
                                         daq_cfg.fadc250_fw)}};
    }
    return {{"error", "channel not found"}};
}

void AppState::projectToHyCalLocal(float Gx, float Gy, float Gz,
                                   float &px, float &py) const
{
    // Transform target and source points into HyCal-local frame, then linearly
    // interpolate along the line to the local z=0 plane.  Equivalent to
    // intersecting the lab-frame line with the tilted HyCal plane, but cleaner
    // because the math reduces to a 1D parametric solve.
    float Tx, Ty, Tz, gx, gy, gz;
    hycal_transform.labToLocal(target_x, target_y, target_z, Tx, Ty, Tz);
    hycal_transform.labToLocal(Gx, Gy, Gz, gx, gy, gz);
    float dz = gz - Tz;
    if (std::abs(dz) < 1e-6f) { px = gx; py = gy; return; }
    float s = -Tz / dz;
    px = Tx + s * (gx - Tx);
    py = Ty + s * (gy - Ty);
}

json AppState::computeClustersJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    if (!cluster_trigger(event.info.trigger_bits))
        return {{"event", ev_id}, {"hits", json::object()}, {"clusters", json::array()},
                {"info", "trigger filtered"}};

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);

    int nmod = hycal.module_count();
    std::vector<float> mod_energy(nmod, 0.f);

    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        auto cit = roc_to_crate.find(roc.tag);
        if (cit == roc_to_crate.end()) continue;
        int crate = cit->second;

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                const auto *mod = hycal.module_by_daq(crate, s, c);
                if (!mod || !mod->is_hycal()) continue;

                if (is_adc1881m) {
                    float adc_val = cd.samples[0];
                    if (adc_val <= 0) continue;
                    float energy = (mod->cal_factor > 0.)
                        ? static_cast<float>(mod->energize(adc_val))
                        : adc_val * adc_to_mev;
                    mod_energy[mod->index] = energy;
                    clusterer.AddHit(mod->index, energy, 0.f);
                    continue;
                }

                ana.SetChannelKey(roc.tag, s, c);
                ana.Analyze(cd.samples, cd.nsamples, wres);
                if (wres.npeaks <= 0) continue;

                if (cluster_cfg.seed_time_window > 0.f) {
                    // Multi-pulse mode — push every detected peak; the
                    // clusterer applies the seed-anchored time gate.
                    float total = 0.f;
                    for (int p = 0; p < wres.npeaks; ++p) {
                        const auto &pk = wres.peaks[p];
                        float adc_val = pk.integral;
                        if (adc_val <= 0) continue;
                        float energy = (mod->cal_factor > 0.)
                            ? static_cast<float>(mod->energize(adc_val))
                            : adc_val * adc_to_mev;
                        clusterer.AddHit(mod->index, energy, pk.time);
                        total += energy;
                    }
                    if (total > 0.f) mod_energy[mod->index] = total;
                } else {
                    // Legacy: largest-integral peak across all detected
                    // peaks (no time gate) — Waveform-Tab peak_filter is
                    // decoupled.
                    float adc_val = bestPeak(wres);
                    if (adc_val <= 0) continue;
                    float energy = (mod->cal_factor > 0.)
                        ? static_cast<float>(mod->energize(adc_val))
                        : adc_val * adc_to_mev;
                    mod_energy[mod->index] = energy;
                    clusterer.AddHit(mod->index, energy, 0.f);
                }
            }
        }
    }

    clusterer.FormClusters();

    json hits_j = json::object();
    for (int i = 0; i < nmod; ++i)
        if (mod_energy[i] > 0.f)
            hits_j[std::to_string(i)] = std::round(mod_energy[i] * 100) / 100;

    std::vector<fdec::HyCalCluster::RecoResult> reco;
    clusterer.ReconstructMatched(reco);

    json cl_arr = json::array();
    for (auto &r : reco) {
        auto &cmod = hycal.module(r.cluster->center.index);
        json indices = json::array();
        for (auto &h : r.cluster->hits) indices.push_back(h.index);
        cl_arr.push_back({
            {"id", static_cast<int>(cl_arr.size())},
            {"center", cmod.name}, {"center_id", cmod.id},
            {"x", std::round(r.hit.x * 10) / 10},
            {"y", std::round(r.hit.y * 10) / 10},
            {"energy", std::round(r.hit.energy * 10) / 10},
            {"nblocks", r.hit.nblocks}, {"npos", r.hit.npos},
            {"modules", indices},
        });
    }

    return {{"event", ev_id}, {"hits", hits_j}, {"clusters", cl_arr}};
}

void AppState::recordSyncTime(uint32_t unix_time, uint64_t last_ti_ts)
{
    if (unix_time == 0) return;
    std::lock_guard<std::mutex> lk(lms_mtx);
    if (sync_unix != 0) return;   // already have a sync reference

    if (lms_first_ts == 0) {
        // No LMS events yet — stash for later.
        // Will be applied when the first LMS event sets lms_first_ts.
        pending_sync_unix = unix_time;
        pending_sync_ti = last_ti_ts;
        return;
    }

    sync_unix = unix_time;
    // ti_delta_sec guards against underflow when last_ti_ts < lms_first_ts
    // (sync arrived earlier than the first LMS event but after some other
    // physics event that already advanced last_ti_ts).
    sync_rel_sec = ti_delta_sec(last_ti_ts, lms_first_ts);
}

void AppState::processEvent(fdec::EventData &event,
                            fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // --- check which consumers need this event ---
    uint32_t tb = event.info.trigger_bits;
    bool do_hist    = waveform_trigger(tb);
    bool do_cluster = cluster_trigger(tb);
    bool do_lms     = lms_trigger.accept != 0 && lms_trigger(tb);
    bool do_alpha   = alpha_trigger.accept != 0 && alpha_trigger(tb);

    if (!do_hist && !do_cluster && !do_lms && !do_alpha) {
        std::lock_guard<std::mutex> lk(data_mtx);
        events_processed++;
        return;
    }

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");

    // clustering setup (stack-allocated, per-event)
    fdec::HyCalCluster clusterer(hycal);
    if (do_cluster) clusterer.SetConfig(cluster_cfg);
    // Sum of every module's reconstructed energy this event — fed into
    // the "Raw Energy Sum" histogram below (sum is taken before clustering
    // so isolated hits below the cluster threshold still count).
    float total_module_energy = 0.f;

    // LMS timing
    double lms_time = 0;

    // acquire both locks for the merged pass
    std::unique_lock<std::mutex> lk1(data_mtx, std::defer_lock);
    std::unique_lock<std::mutex> lk2(lms_mtx, std::defer_lock);
    std::lock(lk1, lk2);

    if (do_lms) {
        if (lms_first_ts == 0) {
            lms_first_ts = event.info.timestamp;
            // apply stashed sync time from a control event that arrived before LMS data
            if (pending_sync_unix != 0 && sync_unix == 0) {
                sync_unix = pending_sync_unix;
                // PRESTART/GO arrives before physics events, so pending_sync_ti is
                // typically 0 → ti_delta_sec returns 0 (run start = LMS start).
                // Also guards the rare reorder case (pending_sync_ti < lms_first_ts).
                sync_rel_sec = ti_delta_sec(pending_sync_ti, lms_first_ts);
            }
        }
        // ti_delta_sec keeps lms_time honest if event.info.timestamp ever
        // comes through as 0 (decoder fall-through) or earlier than the
        // first LMS event captured (out-of-order ET delivery).
        lms_time = ti_delta_sec(event.info.timestamp, lms_first_ts);
    }

    // --- single pass: analyze once per channel, feed all consumers ---
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;

        // crate lookup (needed by cluster + LMS + Alpha consumers)
        int crate = -1;
        if (do_cluster || do_lms || do_alpha) {
            auto cit = roc_to_crate.find(roc.tag);
            if (cit != roc_to_crate.end()) crate = cit->second;
        }

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                // ── analyze ONCE ──
                // peak_for_lms_alpha keeps the historical time-window semantics
                // (now driven by lms_time_min/max from `lms_monitor.time_cut`).
                // peak_for_cluster has no time cut — Waveform-Tab filter is
                // intentionally decoupled from clustering input; per-cluster
                // cuts will be added back later in their own knob.
                float peak_for_lms_alpha = -1;
                float peak_for_cluster   = -1;
                if (!is_adc1881m) {
                    ana.SetChannelKey(roc.tag, s, c);
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    peak_for_lms_alpha = bestPeakInWindow(wres, lms_time_min, lms_time_max);
                    peak_for_cluster   = bestPeak(wres);
                } else {
                    wres.npeaks = 0;
                    peak_for_lms_alpha = peak_for_cluster = cd.samples[0];
                }

                // ── histogram consumer ──
                // Peaks are already gated by the analyzer's height threshold.
                // Apply peak_filter (when enabled) on top.  Time hist gets
                // every passing peak; height/integral hists get the
                // best-integral passing peak (preserves per-event semantics).
                if (do_hist && !is_adc1881m) {
                    std::string key = std::to_string(roc.tag) + "_"
                                   + std::to_string(s) + "_" + std::to_string(c);
                    bool any_peak = false, any_passing = false;
                    float bestI = -1, bestH = -1;
                    for (int p = 0; p < wres.npeaks; ++p) {
                        auto &pk = wres.peaks[p];
                        any_peak = true;
                        if (peak_filter.enable && !peak_filter(pk)) continue;
                        any_passing = true;
                        auto &ph = pos_histograms[key];
                        if (ph.bins.empty()) ph.init(pos_nbins);
                        ph.fill(pk.time, hist_cfg.pos_min, hist_cfg.pos_step);
                        if (pk.integral > bestI) { bestI = pk.integral; bestH = pk.height; }
                    }
                    if (bestI >= 0) {
                        auto &h = histograms[key];
                        if (h.bins.empty()) h.init(hist_nbins);
                        h.fill(bestI, hist_cfg.bin_min, hist_cfg.bin_step);
                        auto &hh = height_histograms[key];
                        if (hh.bins.empty()) hh.init(height_nbins);
                        hh.fill(bestH, hist_cfg.height_min, hist_cfg.height_step);
                    }
                    if (any_peak)    occupancy[key]++;
                    if (any_passing) occupancy_tcut[key]++;
                }

                // ── cluster consumer ──
                if (do_cluster && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod && mod->is_hycal()) {
                        if (is_adc1881m) {
                            float adc_val = (float)cd.samples[0];
                            if (adc_val > 0) {
                                float energy = (mod->cal_factor > 0.)
                                    ? static_cast<float>(mod->energize(adc_val))
                                    : adc_val * adc_to_mev;
                                clusterer.AddHit(mod->index, energy, 0.f);
                                total_module_energy += energy;
                            }
                        } else if (cluster_cfg.seed_time_window > 0.f) {
                            // Multi-pulse: hand every analyzer-detected peak
                            // to the clusterer (which applies the seed-anchored
                            // time gate).  No pre-window — the gate is the cut.
                            for (int p = 0; p < wres.npeaks; ++p) {
                                const auto &pk = wres.peaks[p];
                                if (pk.integral <= 0) continue;
                                float energy = (mod->cal_factor > 0.)
                                    ? static_cast<float>(mod->energize(pk.integral))
                                    : pk.integral * adc_to_mev;
                                clusterer.AddHit(mod->index, energy, pk.time);
                                total_module_energy += energy;
                            }
                        } else if (peak_for_cluster > 0) {
                            float energy = (mod->cal_factor > 0.)
                                ? static_cast<float>(mod->energize(peak_for_cluster))
                                : peak_for_cluster * adc_to_mev;
                            clusterer.AddHit(mod->index, energy, 0.f);
                            total_module_energy += energy;
                        }
                    }
                }

                // ── LMS consumer ──
                if (do_lms && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod) {
                        float val = is_adc1881m ? (float)cd.samples[0] : peak_for_lms_alpha;
                        if (val > 0) {
                            auto &hist = lms_history[mod->index];
                            if (static_cast<int>(hist.size()) < lms_max_history)
                                hist.push_back({lms_time, val});
                            // Always track the latest reading, even after history saturates,
                            // so the LMS/Alpha ref correction stays current.
                            latest_lms_integral[mod->index] = val;
                        }
                    }
                }

                // ── Alpha consumer (Am-241 reference) ──
                if (do_alpha && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod) {
                        float val = is_adc1881m ? (float)cd.samples[0] : peak_for_lms_alpha;
                        if (val > 0) latest_alpha_integral[mod->index] = val;
                    }
                }
            }
        }
    }

    // --- post-loop: clustering + physics histograms ---
    if (do_cluster) {
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> reco_hits;
        clusterer.ReconstructHits(reco_hits);

        struct ClusterInfo { float lx, ly, lz, theta; };
        std::vector<ClusterInfo> cinfo(reco_hits.size());
        for (size_t i = 0; i < reco_hits.size(); ++i) {
            auto &rh = reco_hits[i];
            auto &ci = cinfo[i];
            // Use shower depth as z_local so the lab cluster sits at the
            // shower-max plane, not the front face.  GEM-projection lever
            // arms (hcz / z_gem) depend on this — at PRad scale the
            // ~10 cm depth shifts predicted GEM positions by ~1-2 mm.
            const float z_local = fdec::shower_depth(rh.center_id, rh.energy);
            hycal_transform.toLab(rh.x, rh.y, z_local, ci.lx, ci.ly, ci.lz);
            float dx = ci.lx - target_x, dy = ci.ly - target_y, dz = ci.lz - target_z;
            float rv = std::sqrt(dx*dx + dy*dy);
            ci.theta = std::atan2(rv, dz) * (180.f / 3.14159265f);
        }

        // Per-Ncl bucket index (-1 if Ncl falls outside the nclusters_hist
        // range, in which case the bucketed hists get no fill — same
        // semantics as Histogram::fill underflow/overflow).
        int ncl_bucket = -1;
        {
            float fb = ((float)reco_hits.size() - nclusters_hist_min) / nclusters_hist_step;
            if (fb >= 0.f) {
                int b = (int)fb;
                if (b < (int)cluster_energy_hist_by_ncl.size())
                    ncl_bucket = b;
            }
        }
        for (size_t i = 0; i < reco_hits.size(); ++i) {
            cluster_energy_hist.fill(reco_hits[i].energy, cl_hist_min, cl_hist_step);
            nblocks_hist.fill(reco_hits[i].nblocks, nblocks_hist_min, nblocks_hist_step);
            if (ncl_bucket >= 0) {
                cluster_energy_hist_by_ncl[ncl_bucket].fill(
                    reco_hits[i].energy, cl_hist_min, cl_hist_step);
                nblocks_hist_by_ncl[ncl_bucket].fill(
                    reco_hits[i].nblocks, nblocks_hist_min, nblocks_hist_step);
            }
        }
        nclusters_hist.fill(reco_hits.size(), nclusters_hist_min, nclusters_hist_step);
        raw_energy_hist.fill(total_module_energy, raw_energy_hist_min, raw_energy_hist_step);
        cluster_events_processed++;

        bool physics_accept = physics_trigger(tb);
        if (physics_accept) {
            float Eb = beam_energy.load();
            for (size_t i = 0; i < reco_hits.size(); ++i) {
                energy_angle_hist.fill(cinfo[i].theta, reco_hits[i].energy,
                    ea_angle_min, ea_angle_step, ea_energy_min, ea_energy_step);
            }
            if (reco_hits.size() == 2 && Eb > 0) {
                float esum = reco_hits[0].energy + reco_hits[1].energy;
                bool energy_ok = std::abs(esum - Eb) < moller_energy_tol * Eb;
                bool angle_ok = false;
                for (int j = 0; j < 2; ++j)
                    if (cinfo[j].theta >= moller_angle_min && cinfo[j].theta <= moller_angle_max)
                        angle_ok = true;
                if (energy_ok && angle_ok) {
                    moller_events++;
                    for (int j = 0; j < 2; ++j) {
                        moller_xy_hist.fill(cinfo[j].lx, cinfo[j].ly,
                            moller_xy_x_min, moller_xy_x_step, moller_xy_y_min, moller_xy_y_step);
                    }
                }
            }
            // HyCal cluster-hit XY: single-cluster ep-elastic candidates.
            // Same gate is reused (when configured) for GEM↔HyCal residuals.
            bool ep_cand = false;
            if ((int)reco_hits.size() == hxy_n_clusters && Eb > 0) {
                const auto &cl = reco_hits[0];
                bool nb_ok = cl.nblocks >= hxy_nblocks_min && cl.nblocks <= hxy_nblocks_max;
                bool e_ok  = cl.energy  >= hxy_energy_frac_min * Eb;
                if (nb_ok && e_ok) {
                    ep_cand = true;
                    hycal_xy_hist.fill(cinfo[0].lx, cinfo[0].ly,
                        hxy_x_min, hxy_x_step, hxy_y_min, hxy_y_step);
                    hycal_xy_events++;
                }
            }
            // GEM↔HyCal matching residuals.  Reference is the FIRST cluster's
            // HyCal-local xy — for ep candidates that's the only cluster, for
            // multi-cluster events it's the leading reconstructed hit.  The
            // residual lives at the HyCal plane, so σ_GEM is projected through
            // the target onto that plane: σ_total² = σ_HC(E)² + (σ_GEM·z_hc/z_gem)².
            if (gem_enabled && (ep_cand || !gem_match_require_ep) && !reco_hits.empty()) {
                const float ref_x = reco_hits[0].x, ref_y = reco_hits[0].y;
                const float sigma_hc = hycal.PositionResolution(reco_hits[0].energy);
                const float z_hc = cinfo[0].lz;
                const int n_dets = std::min<int>(gem_sys.GetNDetectors(),
                                                 (int)gem_dx_hist.size());
                for (int d = 0; d < n_dets; ++d) {
                    auto &xform = gem_transforms[d];
                    const float z_gem  = (xform.z != 0.f) ? xform.z : 1.f;
                    const float s_gem  = (d < (int)gem_pos_res.size())
                                            ? gem_pos_res[d] : 0.1f;
                    const float s_gem_at_hc = s_gem * std::abs(z_hc / z_gem);
                    const float s_total = std::sqrt(sigma_hc*sigma_hc
                                                  + s_gem_at_hc*s_gem_at_hc);
                    const float cut = gem_match_nsigma * s_total;
                    for (auto &h : gem_sys.GetHits(d)) {
                        float lx, ly, lz;
                        xform.toLab(h.x, h.y, lx, ly, lz);
                        float px, py;
                        projectToHyCalLocal(lx, ly, lz, px, py);
                        float dxr = px - ref_x, dyr = py - ref_y;
                        if (std::sqrt(dxr*dxr + dyr*dyr) < cut) {
                            gem_dx_hist[d].fill(dxr, gem_resid_min, gem_resid_step);
                            gem_dy_hist[d].fill(dyr, gem_resid_min, gem_resid_step);
                            gem_match_hits[d]++;
                        }
                    }
                }
                gem_match_events++;
            }
            // GEM tracking efficiency — per HyCal cluster, no target assumption.
            // Builds per-detector lab-frame hit lists once and runs Pass A / B
            // for each cluster passing min_cluster_energy.
            if (gem_enabled && !reco_hits.empty()) {
                const int n_gem = std::min<int>(gem_sys.GetNDetectors(),
                                                (int)gem_transforms.size());
                std::vector<std::vector<LabHit>> hits_by_det(n_gem);
                for (int d = 0; d < n_gem; ++d) {
                    auto &xform = gem_transforms[d];
                    for (auto &h : gem_sys.GetHits(d)) {
                        float lx, ly, lz;
                        xform.toLab(h.x, h.y, lx, ly, lz);
                        hits_by_det[d].push_back({lx, ly, lz});
                    }
                }
                for (size_t i = 0; i < reco_hits.size(); ++i) {
                    if (reco_hits[i].energy < gem_eff_min_cluster_energy) continue;
                    runGemEfficiency((int)event.info.event_number,
                                     cinfo[i].lx, cinfo[i].ly, cinfo[i].lz,
                                     reco_hits[i].energy,
                                     hits_by_det);
                }
            }
        }
    }

    events_processed++;
    if (do_lms) lms_events++;
}

void AppState::processReconEvent(const ReconEventData &recon)
{
    uint32_t tb = recon.trigger_bits;
    bool do_cluster = cluster_trigger(tb);
    bool do_physics = physics_trigger(tb);

    std::lock_guard<std::mutex> lk(data_mtx);
    events_processed++;

    if (do_cluster && !recon.clusters.empty()) {
        int ncl_bucket = -1;
        {
            float fb = ((float)recon.clusters.size() - nclusters_hist_min)
                       / nclusters_hist_step;
            if (fb >= 0.f) {
                int b = (int)fb;
                if (b < (int)cluster_energy_hist_by_ncl.size())
                    ncl_bucket = b;
            }
        }
        // Recon path has no per-module energies, so the closest analog
        // to the raw energy sum is the sum of cluster energies (the same
        // value the GUI shows in the cluster-tab summary as "ECl Sum").
        float total_cluster_energy = 0.f;
        for (auto &cl : recon.clusters) {
            cluster_energy_hist.fill(cl.energy, cl_hist_min, cl_hist_step);
            nblocks_hist.fill(cl.nblocks, nblocks_hist_min, nblocks_hist_step);
            total_cluster_energy += cl.energy;
            if (ncl_bucket >= 0) {
                cluster_energy_hist_by_ncl[ncl_bucket].fill(
                    cl.energy, cl_hist_min, cl_hist_step);
                nblocks_hist_by_ncl[ncl_bucket].fill(
                    cl.nblocks, nblocks_hist_min, nblocks_hist_step);
            }
        }
        nclusters_hist.fill(recon.clusters.size(), nclusters_hist_min, nclusters_hist_step);
        raw_energy_hist.fill(total_cluster_energy, raw_energy_hist_min, raw_energy_hist_step);
        cluster_events_processed++;
    }

    if (do_physics && !recon.clusters.empty()) {
        struct CI { float lx, ly, lz, theta; };
        std::vector<CI> cinfo(recon.clusters.size());
        for (size_t i = 0; i < recon.clusters.size(); ++i) {
            auto &cl = recon.clusters[i];
            auto &ci = cinfo[i];
            // Same shower-depth correction as the EVIO path above.
            const float z_local = fdec::shower_depth(cl.center_id, cl.energy);
            hycal_transform.toLab(cl.x, cl.y, z_local, ci.lx, ci.ly, ci.lz);
            float dx = ci.lx - target_x, dy = ci.ly - target_y, dz = ci.lz - target_z;
            float r = std::sqrt(dx*dx + dy*dy);
            ci.theta = std::atan2(r, dz) * (180.f / 3.14159265f);
        }
        for (size_t i = 0; i < recon.clusters.size(); ++i)
            energy_angle_hist.fill(cinfo[i].theta, recon.clusters[i].energy,
                ea_angle_min, ea_angle_step, ea_energy_min, ea_energy_step);

        float Eb = beam_energy.load();
        if (recon.clusters.size() == 2 && Eb > 0) {
            float esum = recon.clusters[0].energy + recon.clusters[1].energy;
            bool energy_ok = std::abs(esum - Eb) < moller_energy_tol * Eb;
            bool angle_ok = false;
            for (int j = 0; j < 2; ++j)
                if (cinfo[j].theta >= moller_angle_min && cinfo[j].theta <= moller_angle_max)
                    angle_ok = true;
            if (energy_ok && angle_ok) {
                moller_events++;
                for (int j = 0; j < 2; ++j) {
                    moller_xy_hist.fill(cinfo[j].lx, cinfo[j].ly,
                        moller_xy_x_min, moller_xy_x_step, moller_xy_y_min, moller_xy_y_step);
                }
            }
        }
        // HyCal cluster-hit XY: single-cluster ep-elastic candidates.
        bool ep_cand = false;
        if ((int)recon.clusters.size() == hxy_n_clusters && Eb > 0) {
            const auto &cl = recon.clusters[0];
            bool nb_ok = cl.nblocks >= hxy_nblocks_min && cl.nblocks <= hxy_nblocks_max;
            bool e_ok  = cl.energy  >= hxy_energy_frac_min * Eb;
            if (nb_ok && e_ok) {
                ep_cand = true;
                hycal_xy_hist.fill(cinfo[0].lx, cinfo[0].ly,
                    hxy_x_min, hxy_x_step, hxy_y_min, hxy_y_step);
                hycal_xy_events++;
            }
        }
        // GEM↔HyCal matching residuals (ROOT recon path uses recon.gem_hits,
        // which carry detector-local x,y just like the live gem_sys hits).
        // Same parametric cut as the live path:
        //   σ_total² = σ_HC(E)² + (σ_GEM·z_hc/z_gem)²,  cut = nsigma · σ_total.
        if (gem_enabled && (ep_cand || !gem_match_require_ep) && !recon.clusters.empty()) {
            const float ref_x = recon.clusters[0].x, ref_y = recon.clusters[0].y;
            const float sigma_hc = hycal.PositionResolution(recon.clusters[0].energy);
            const float z_hc = cinfo[0].lz;
            const int n_dets = (int)gem_dx_hist.size();
            for (auto &gh : recon.gem_hits) {
                if (gh.det_id < 0 || gh.det_id >= n_dets) continue;
                if (gh.det_id >= (int)gem_transforms.size()) continue;
                auto &xform = gem_transforms[gh.det_id];
                float lx, ly, lz;
                xform.toLab(gh.x, gh.y, lx, ly, lz);
                float px, py;
                projectToHyCalLocal(lx, ly, lz, px, py);
                float dxr = px - ref_x, dyr = py - ref_y;
                const float z_gem = (xform.z != 0.f) ? xform.z : 1.f;
                const float s_gem = (gh.det_id < (int)gem_pos_res.size())
                                        ? gem_pos_res[gh.det_id] : 0.1f;
                const float s_gem_at_hc = s_gem * std::abs(z_hc / z_gem);
                const float s_total = std::sqrt(sigma_hc*sigma_hc
                                              + s_gem_at_hc*s_gem_at_hc);
                const float cut = gem_match_nsigma * s_total;
                if (std::sqrt(dxr*dxr + dyr*dyr) < cut) {
                    gem_dx_hist[gh.det_id].fill(dxr, gem_resid_min, gem_resid_step);
                    gem_dy_hist[gh.det_id].fill(dyr, gem_resid_min, gem_resid_step);
                    gem_match_hits[gh.det_id]++;
                }
            }
            gem_match_events++;
        }
        // GEM tracking efficiency (recon path mirrors the live-data path).
        if (gem_enabled && !recon.clusters.empty()) {
            const int n_gem = std::min<int>(gem_sys.GetNDetectors(),
                                            (int)gem_transforms.size());
            std::vector<std::vector<LabHit>> hits_by_det(n_gem);
            for (auto &gh : recon.gem_hits) {
                if (gh.det_id < 0 || gh.det_id >= n_gem) continue;
                auto &xform = gem_transforms[gh.det_id];
                float lx, ly, lz;
                xform.toLab(gh.x, gh.y, lx, ly, lz);
                hits_by_det[gh.det_id].push_back({lx, ly, lz});
            }
            for (size_t i = 0; i < recon.clusters.size(); ++i) {
                if (recon.clusters[i].energy < gem_eff_min_cluster_energy) continue;
                runGemEfficiency(recon.event_num,
                                 cinfo[i].lx, cinfo[i].ly, cinfo[i].lz,
                                 recon.clusters[i].energy,
                                 hits_by_det);
            }
        }
    }
}

json AppState::encodeReconClustersJson(const ReconEventData &recon, int ev_id)
{
    json hits_j = json::object();
    json cl_arr = json::array();

    for (size_t i = 0; i < recon.clusters.size(); ++i) {
        auto &cl = recon.clusters[i];
        std::string center_name;
        if (cl.center_id >= 0 && cl.center_id < hycal.module_count())
            center_name = hycal.module(cl.center_id).name;
        hits_j[std::to_string(cl.center_id)] =
            std::round(cl.energy * 100) / 100;
        cl_arr.push_back({
            {"id", (int)i}, {"center", center_name},
            {"center_id", cl.center_id},
            {"x", std::round(cl.x * 10) / 10},
            {"y", std::round(cl.y * 10) / 10},
            {"energy", std::round(cl.energy * 10) / 10},
            {"nblocks", cl.nblocks}, {"npos", 0},
            {"modules", json::array({cl.center_id})},
        });
    }
    return {{"event", ev_id}, {"hits", hits_j}, {"clusters", cl_arr}};
}

void AppState::prepareGemForView(const ssp::SspEventData &ssp_evt)
{
    if (!gem_enabled || ssp_evt.nmpds == 0) return;
    gem_sys.Clear();
    gem_sys.ProcessEvent(ssp_evt);
    gem_sys.Reconstruct(gem_clusterer);
}

void AppState::processGemEvent(const ssp::SspEventData &ssp_evt)
{
    if (!gem_enabled || ssp_evt.nmpds == 0) return;
    prepareGemForView(ssp_evt);

    // Strip-level diagnostic: fill per-detector occupancy over the active
    // strip extent (smaller than PlaneConfig.size on the beam-hole side
    // because pos=11 reuses pos=10 via shared_pos) so the heatmap matches
    // the dashed detector frame drawn in the GUI.  No target assumption is
    // made here — lab-frame plots live in the matching/efficiency views.
    std::lock_guard<std::mutex> lk(data_mtx);
    const int n_dets = std::min<int>(gem_sys.GetNDetectors(),
                                     (int)gem_transforms.size());
    for (int d = 0; d < n_dets; ++d) {
        auto xr = gem_sys.GetActiveExtent(d, 0);
        auto yr = gem_sys.GetActiveExtent(d, 1);
        const float xStep = (xr.second - xr.first) / GEM_OCC_NX;
        const float yStep = (yr.second - yr.first) / GEM_OCC_NY;
        for (auto &h : gem_sys.GetHits(d))
            gem_occupancy[d].fill(h.x, h.y, xr.first, xStep, yr.first, yStep);
    }
}

//=============================================================================
// GEM API builders
//=============================================================================

nlohmann::json AppState::apiGemHits() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    result["n_detectors"] = gem_sys.GetNDetectors();
    json detectors = json::array();
    json all_hits = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json dj;
        dj["name"] = det.name;
        dj["id"] = det.id;

        // 1D clusters per plane
        for (int p = 0; p < 2; ++p) {
            std::string pname = (p == 0) ? "x_clusters" : "y_clusters";
            json clusters = json::array();
            for (auto &cl : gem_sys.GetPlaneClusters(d, p)) {
                clusters.push_back({
                    {"position", cl.position},
                    {"peak_charge", cl.peak_charge},
                    {"total_charge", cl.total_charge},
                    {"size", (int)cl.hits.size()},
                    {"max_timebin", cl.max_timebin}
                });
            }
            dj[pname] = clusters;
        }

        // 2D hits (transformed to lab frame) — build per-det and all_hits in one pass.
        // proj_x/proj_y are the lab→target line projected onto the HyCal local
        // plane; the cluster tab overlays these on the geo view.
        auto &xform = gem_transforms[d];
        json hits = json::array();
        for (auto &h : gem_sys.GetHits(d)) {
            float lx, ly, lz;
            xform.toLab(h.x, h.y, lx, ly, lz);
            float px, py;
            projectToHyCalLocal(lx, ly, lz, px, py);
            hits.push_back({
                {"x", lx}, {"y", ly},
                {"proj_x", px}, {"proj_y", py},
                {"x_charge", h.x_charge}, {"y_charge", h.y_charge},
                {"x_size", h.x_size}, {"y_size", h.y_size}
            });
            all_hits.push_back({
                {"x", lx}, {"y", ly},
                {"proj_x", px}, {"proj_y", py}, {"det", d},
                {"x_charge", h.x_charge}, {"y_charge", h.y_charge}
            });
        }
        dj["hits_2d"] = hits;
        detectors.push_back(dj);
    }
    result["detectors"] = detectors;
    result["all_hits"] = all_hits;
    return result;
}

nlohmann::json AppState::apiGemConfig() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    result["n_detectors"] = gem_sys.GetNDetectors();
    json layers = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json lj = {
            {"id", det.id},
            {"name", det.name},
            {"type", det.type},
            {"x_pitch", det.planes[0].pitch},
            {"y_pitch", det.planes[1].pitch},
            {"x_apvs", det.planes[0].n_apvs},
            {"y_apvs", det.planes[1].n_apvs},
            {"x_size", det.planes[0].size},
            {"y_size", det.planes[1].size}
        };
        auto &t = gem_transforms[d];
        lj["position"] = json::array({t.x, t.y, t.z});
        lj["tilting"]  = json::array({t.rx, t.ry, t.rz});
        layers.push_back(lj);
    }
    result["layers"] = layers;
    result["occ_nx"] = GEM_OCC_NX;
    result["occ_ny"] = GEM_OCC_NY;
    return result;
}

nlohmann::json AppState::apiGemOccupancy() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    std::lock_guard<std::mutex> lk(data_mtx);
    json dets = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        // Active strip extent in detector-local coords (mm).  See
        // GemSystem::GetActiveExtent — tighter than PlaneConfig.size on the
        // inner-edge side; matches the bin range used in processGemEvent.
        auto xr = gem_sys.GetActiveExtent(d, 0);
        auto yr = gem_sys.GetActiveExtent(d, 1);
        json dj;
        dj["id"] = det.id;
        dj["name"] = det.name;
        dj["x_size"] = det.planes[0].size;
        dj["y_size"] = det.planes[1].size;
        dj["x_active"] = json::array({xr.first, xr.second});
        dj["y_active"] = json::array({yr.first, yr.second});
        dj["nx"] = GEM_OCC_NX;
        dj["ny"] = GEM_OCC_NY;
        dj["bins"] = gem_occupancy[d].bins;
        dets.push_back(dj);
    }
    result["detectors"] = dets;
    result["total"] = events_processed.load();
    return result;
}

nlohmann::json AppState::apiGemApv(const ssp::SspEventData &ssp_evt, int evnum,
                                   bool skip_sw_zs,
                                   bool *any_full_readout) const
{
    if (any_full_readout) *any_full_readout = false;
    json result = json::object();
    result["enabled"] = gem_enabled;
    result["event"]   = evnum;
    if (!gem_enabled) {
        result["detectors"] = json::array();
        result["apvs"]      = json::array();
        return result;
    }
    // Global software N-sigma multiplier and pedestal calibration revision.
    // The frontend pairs this zs_sigma with the per-APV noise[] from
    // /api/gem/calib to draw the threshold band; gem_calib_rev lets it
    // detect when the cached calib payload is stale and needs re-fetching.
    result["zs_sigma"]  = gem_sys.GetZeroSupThreshold();
    result["calib_rev"] = gem_calib_rev.load();

    // Detector summary — det_id → name + APV count, used by the frontend
    // to lay out one section per GEM with consistent ordering.
    auto &dets = gem_sys.GetDetectors();
    json det_arr = json::array();
    std::vector<int> apv_counts(dets.size(), 0);
    for (int i = 0; i < gem_sys.GetNApvs(); ++i) {
        auto &cfg = gem_sys.GetApvConfig(i);
        if (cfg.det_id >= 0 && cfg.det_id < (int)apv_counts.size())
            apv_counts[cfg.det_id]++;
    }
    for (size_t d = 0; d < dets.size(); ++d) {
        det_arr.push_back({
            {"id",      dets[d].id},
            {"name",    dets[d].name},
            {"n_apvs",  apv_counts[d]},
        });
    }
    result["detectors"] = det_arr;

    // Per-APV dump.  Each APV carries:
    //   raw[128][6]        — int16 firmware samples (0 if APV not in event)
    //   processed[128][6]  — pedestal + CM corrected float (0 if not present)
    //   hits[128]          — software ZS survivor (post-cut), 0/1
    //   fw_hits[128]       — firmware survivor (pre-software-cut), 0/1
    //   cm[6] | null       — firmware online common mode per time sample
    //   no_hit_fr          — firmware full-readout (nstrips==128) but no survivors
    //   full_readout       — firmware sent all 128 channels (nstrips==128)
    //   present            — APV showed up in this event's SSP data
    // Per-APV pedestal noise lives on /api/gem/calib (one-shot, cached
    // by the frontend until calib_rev changes).
    constexpr int N_STRIPS = 128;
    constexpr int N_TS     = 6;
    json apvs = json::array();
    for (int i = 0; i < gem_sys.GetNApvs(); ++i) {
        auto &cfg = gem_sys.GetApvConfig(i);
        if (cfg.crate_id < 0 || cfg.mpd_id < 0 || cfg.adc_ch < 0)
            continue;   // unmapped slot

        const ssp::ApvData *raw = ssp_evt.findApv(cfg.crate_id, cfg.mpd_id, cfg.adc_ch);
        bool present = (raw != nullptr) && raw->present;
        // Firmware full-readout: every channel present in the SSP stream.
        // Used by the "Latest full-readout" mode to single out the
        // prescaled monitoring events that bypass online ZS.
        bool full_readout = present && raw->nstrips >= N_STRIPS;
        if (full_readout && any_full_readout) *any_full_readout = true;

        json raw_arr = json::array();
        json proc_arr = json::array();
        json hit_arr = json::array();
        json fw_hit_arr = json::array();
        bool any_hit = false;

        for (int s = 0; s < N_STRIPS; ++s) {
            json raw_row  = json::array();
            json proc_row = json::array();
            for (int t = 0; t < N_TS; ++t) {
                if (present)
                    raw_row.push_back(static_cast<int>(raw->strips[s][t]));
                else
                    raw_row.push_back(0);
                // 1-decimal rounding keeps payload tight without losing
                // anything visible at the panel's sub-pixel resolution.
                float v = present ? gem_sys.GetProcessedAdc(i, s, t) : 0.f;
                proc_row.push_back(std::round(v * 10.f) / 10.f);
            }
            raw_arr.push_back(std::move(raw_row));
            proc_arr.push_back(std::move(proc_row));
            // In snapshot (skip_sw_zs) mode, every channel of a full-readout
            // APV is marked as a hit so the client's signal-only filter
            // doesn't hide it.  Firmware-ZS'd APVs keep their normal mask
            // (mixed events are rare in practice but handled correctly).
            bool sw_hit = present && gem_sys.IsChannelHit(i, s);
            bool hit = (skip_sw_zs && full_readout) ? true : sw_hit;
            if (hit) any_hit = true;
            hit_arr.push_back(hit ? 1 : 0);
            fw_hit_arr.push_back(present && raw->hasStrip(s) ? 1 : 0);
        }

        // Firmware online CM (6 samples) — only present when the MPD emitted
        // type-0xD debug-header words; otherwise null so the frontend can
        // skip the overlay rather than draw zeros.
        json cm_val = nullptr;
        if (present && raw->has_online_cm) {
            json cm_arr = json::array();
            for (int t = 0; t < N_TS; ++t)
                cm_arr.push_back(static_cast<int>(raw->online_cm[t]));
            cm_val = std::move(cm_arr);
        }

        std::string det_name;
        if (cfg.det_id >= 0 && cfg.det_id < (int)dets.size())
            det_name = dets[cfg.det_id].name;
        const char *plane = (cfg.plane_type == 0) ? "X"
                          : (cfg.plane_type == 1) ? "Y" : "?";

        // Firmware full-readout warning: APV reported every strip but none
        // survived ZS.  Mirrors the gem_event_viewer.py "no hits" badge.
        // any_hit here reflects the (possibly overridden) hits[] mask, so
        // in snapshot mode a full-readout APV with no real survivors will
        // still report no_hit_fr=false — which is exactly what the operator
        // wants (the badge is for diagnosing "where did the hits go?",
        // which is moot when we explicitly disabled the cut).
        bool no_hit_fr = full_readout && !any_hit;

        apvs.push_back({
            {"id",            i},
            {"det_id",        cfg.det_id},
            {"det_name",      det_name},
            {"plane",         plane},
            {"det_pos",       cfg.det_pos},
            {"crate",         cfg.crate_id},
            {"mpd",           cfg.mpd_id},
            {"adc",           cfg.adc_ch},
            {"present",       present},
            {"no_hit_fr",     no_hit_fr},
            {"full_readout",  full_readout},
            {"raw",           std::move(raw_arr)},
            {"processed",     std::move(proc_arr)},
            {"hits",          std::move(hit_arr)},
            {"fw_hits",       std::move(fw_hit_arr)},
            {"cm",            std::move(cm_val)},
        });
    }
    result["apvs"] = apvs;
    return result;
}

nlohmann::json AppState::apiGemCalib() const
{
    json result = json::object();
    result["enabled"]   = gem_enabled;
    result["rev"]       = gem_calib_rev.load();
    result["zs_sigma"]  = gem_enabled ? gem_sys.GetZeroSupThreshold() : 0.f;
    json apvs = json::array();
    if (gem_enabled) {
        constexpr int N_STRIPS = 128;
        for (int i = 0; i < gem_sys.GetNApvs(); ++i) {
            auto &cfg = gem_sys.GetApvConfig(i);
            if (cfg.crate_id < 0 || cfg.mpd_id < 0 || cfg.adc_ch < 0)
                continue;
            json noise_arr = json::array();
            for (int s = 0; s < N_STRIPS; ++s)
                noise_arr.push_back(std::round(cfg.pedestal[s].noise * 10.f) / 10.f);
            apvs.push_back({{"id", i}, {"noise", std::move(noise_arr)}});
        }
    }
    result["apvs"] = std::move(apvs);
    return result;
}

void AppState::setGemZsSigma(float v)
{
    if (v < 0.f) v = 0.f;
    gem_sys.SetZeroSupThreshold(v);
}

//=============================================================================
// GEM efficiency monitor — main entry, init, clear, snapshot serialization
//=============================================================================

void AppState::initGemEfficiency()
{
    int n_gem = gem_enabled ? (int)gem_transforms.size() : 0;
    gem_eff_num.assign(n_gem, 0);
    gem_eff_den.assign(n_gem, 0);
    gem_eff_grid_num.assign(n_gem, Histogram2D{});
    gem_eff_grid_den.assign(n_gem, Histogram2D{});
    for (int d = 0; d < n_gem; ++d) {
        gem_eff_grid_num[d].init(gem_eff_grid_nx, gem_eff_grid_ny);
        gem_eff_grid_den[d].init(gem_eff_grid_nx, gem_eff_grid_ny);
    }
    gem_eff_snapshot = GemEffSnapshot{};
    int nbins = (int)std::lround((gem_eff_z_target_max - gem_eff_z_target_min)
                                  / gem_eff_z_target_step);
    gem_eff_z_target_hist.init(nbins);
}

void AppState::clearGemEfficiency()
{
    for (auto &n : gem_eff_num) n = 0;
    for (auto &n : gem_eff_den) n = 0;
    for (auto &h : gem_eff_grid_num) h.clear();
    for (auto &h : gem_eff_grid_den) h.clear();
    for (int d = 0; d < 4; ++d) {
        gem_eff_diag_call[d] = 0;
        gem_eff_diag_3matched[d] = 0;
        gem_eff_diag_pass_chi2[d] = 0;
        gem_eff_diag_pass_resid[d] = 0;
    }
    gem_eff_snapshot = GemEffSnapshot{};
    gem_eff_z_target_hist.clear();
}

void AppState::runGemEfficiency(int event_id,
                                float hcx, float hcy, float hcz, float hc_energy,
                                const std::vector<std::vector<LabHit>> &hits_by_det)
{
    if (!gem_enabled) return;
    int n_dets = std::min((int)hits_by_det.size(),
                  std::min((int)gem_transforms.size(), gem_sys.GetNDetectors()));
    n_dets = std::min(n_dets, GEM_EFF_MAX_DETS);
    if (n_dets < 3) return;

    // Resolutions (mm) — σ_HC(E) at HyCal face, σ_GEM[d] at GEM plane.
    // Target-seeded line goes (target → HyCal), so the projection error at
    // depth z scales linearly from 0 at the target to σ_HC at HyCal:
    //   σ_HC@gem = σ_HC · |(z_gem - z_target) / (z_hc - z_target)|
    //   σ_total  = sqrt(σ_HC@gem² + σ_GEM²),  cut = match_nsigma · σ_total
    const float sigma_hc = hycal.PositionResolution(hc_energy);
    const float w_h      = 1.f / (sigma_hc * sigma_hc);
    const float lever_hc = (hcz != target_z) ? (hcz - target_z) : 1.f;
    auto sigmaGem = [&](int d) -> float {
        return (d >= 0 && d < (int)gem_pos_res.size()) ? gem_pos_res[d] : 0.1f;
    };

    // Find closest GEM-d hit (in detector-local coords) within
    // match_nsigma · σ_total of a predicted local point.  -1 if none in window.
    auto findClosest = [&](int d, float pred_lx, float pred_ly,
                           int &out_idx,
                           float &out_lab_x, float &out_lab_y, float &out_lab_z) {
        out_idx = -1;
        if (d < 0 || d >= n_dets) return;
        const auto &hits = hits_by_det[d];
        int max_n = std::min((int)hits.size(), gem_eff_max_hits_per_det);
        if (max_n == 0) return;
        const auto &xform = gem_transforms[d];
        const float z_gem        = xform.z;
        const float s_hc_at_gem  = sigma_hc * std::abs((z_gem - target_z) / lever_hc);
        const float s_gem        = sigmaGem(d);
        const float s_total      = std::sqrt(s_hc_at_gem*s_hc_at_gem
                                            + s_gem*s_gem);
        const float cut          = gem_eff_match_nsigma * s_total;
        float best_d2 = cut * cut;
        for (int i = 0; i < max_n; ++i) {
            const auto &h = hits[i];
            float lx, ly, lz;
            xform.labToLocal(h[0], h[1], h[2], lx, ly, lz);
            float dx = lx - pred_lx, dy = ly - pred_ly;
            float d2 = dx*dx + dy*dy;
            if (d2 < best_d2) {
                best_d2 = d2;
                out_idx = i;
                out_lab_x = h[0]; out_lab_y = h[1]; out_lab_z = h[2];
            }
        }
    };

    auto inside = [&](int d, float lx, float ly) -> bool {
        if (d < 0 || d >= gem_sys.GetNDetectors()) return false;
        const auto &det = gem_sys.GetDetectors()[d];
        float xmax = det.planes[0].size * 0.5f;
        float ymax = det.planes[1].size * 0.5f;
        return std::abs(lx) <= xmax && std::abs(ly) <= ymax;
    };

    // ---- LOO anchor finder ------------------------------------------------
    // Try a single seed (HyCal, S_hit) and return the resulting fit if it
    // passes χ² + per-detector residual gates.  S = -1 means target-seeded
    // (line goes target → HyCal, no GEM in the seed).  candidate_dets is a
    // bitmask: bit d set ⇒ detector d is allowed as a match candidate.  The
    // anchor never accesses the bit-cleared detectors (i.e. the test
    // detector is invisible to the anchor).
    constexpr int CAP = GEM_EFF_MAX_DETS + 1 + 1;  // +HyCal +target
    struct Anchor {
        bool   valid = false;
        Line3D fit;
        bool   matched[GEM_EFF_MAX_DETS] = {false,false,false,false};
        float  cand_lx[GEM_EFF_MAX_DETS] = {0,0,0,0};
        float  cand_ly[GEM_EFF_MAX_DETS] = {0,0,0,0};
        float  cand_lz[GEM_EFF_MAX_DETS] = {0,0,0,0};
    };

    // Target-seeded LOO is the production-default path; we use a passed-in
    // test_d to bump the per-stage diag counters for that mode only.  GEM-
    // seeded modes try many seeds per (HyCal, test_d), so per-seed gates
    // don't map 1:1 onto the audit script's per-anchor breakdown.
    auto trySeed = [&](int S, int seed_idx, unsigned cand_mask,
                       bool target_in_fit, int diag_test_d) -> Anchor {
        Anchor a{};
        Line3D seed;
        if (S < 0) {
            seed = seedLine(target_x, target_y, target_z, hcx, hcy, hcz);
        } else {
            const auto &hits_s = hits_by_det[S];
            if (seed_idx < 0 || seed_idx >= (int)hits_s.size()) return a;
            const auto &g0 = hits_s[seed_idx];
            seed = seedLine(hcx, hcy, hcz, g0[0], g0[1], g0[2]);
            a.matched[S] = true;
            a.cand_lx[S] = g0[0]; a.cand_ly[S] = g0[1]; a.cand_lz[S] = g0[2];
        }
        for (int d = 0; d < n_dets; ++d) {
            if (!(cand_mask & (1u << d))) continue;
            if (a.matched[d]) continue;          // S already matched
            float pred_lx, pred_ly;
            projectLineToLocal(gem_transforms[d], seed, pred_lx, pred_ly);
            int idx; float lab_x, lab_y, lab_z;
            findClosest(d, pred_lx, pred_ly, idx, lab_x, lab_y, lab_z);
            if (idx >= 0) {
                a.matched[d] = true;
                a.cand_lx[d] = lab_x; a.cand_ly[d] = lab_y; a.cand_lz[d] = lab_z;
            }
        }
        // Need all 3 candidate detectors to have matched (LOO requires the
        // OTHER 3 GEMs all delivered a hit; this is the "clean basis" cut).
        int nmatch = 0;
        for (int d = 0; d < n_dets; ++d) if (a.matched[d]) ++nmatch;
        if (nmatch < 3) return a;
        if (S < 0 && diag_test_d >= 0 && diag_test_d < 4)
            gem_eff_diag_3matched[diag_test_d]++;
        // Build fit arrays: HyCal (+ optionally target) + matched GEMs.
        float zarr[CAP], xarr[CAP], yarr[CAP], wxarr[CAP], wyarr[CAP];
        int   N = 0;
        zarr[N] = hcz; xarr[N] = hcx; yarr[N] = hcy;
        wxarr[N] = wyarr[N] = w_h; ++N;
        if (target_in_fit) {
            // σ_target_z couples to the transverse measurement at z=target_z
            // through the slope: σ_x_eff² = σ_target_x² + (bx_est·σ_z)².
            // Slope estimate from the target → HyCal lever arm.
            const float lever = (hcz != target_z) ? (hcz - target_z) : 1.f;
            const float bx_est = (hcx - target_x) / lever;
            const float by_est = (hcy - target_y) / lever;
            const float sx2 = gem_eff_target_sigma_x * gem_eff_target_sigma_x
                              + (bx_est * gem_eff_target_sigma_z)
                              * (bx_est * gem_eff_target_sigma_z);
            const float sy2 = gem_eff_target_sigma_y * gem_eff_target_sigma_y
                              + (by_est * gem_eff_target_sigma_z)
                              * (by_est * gem_eff_target_sigma_z);
            zarr[N] = target_z; xarr[N] = target_x; yarr[N] = target_y;
            wxarr[N] = 1.f / sx2; wyarr[N] = 1.f / sy2;
            ++N;
        }
        for (int d = 0; d < n_dets; ++d) {
            if (!a.matched[d]) continue;
            zarr[N] = a.cand_lz[d]; xarr[N] = a.cand_lx[d]; yarr[N] = a.cand_ly[d];
            const float s = sigmaGem(d);
            wxarr[N] = wyarr[N] = 1.f / (s * s);
            ++N;
        }
        if (!fitWeightedLine(N, zarr, xarr, yarr, wxarr, wyarr, a.fit)) return a;
        if (a.fit.chi2_per_dof > gem_eff_max_chi2) {
            a.matched[0] = a.matched[1] = a.matched[2] = a.matched[3] = false;
            return a;
        }
        if (S < 0 && diag_test_d >= 0 && diag_test_d < 4)
            gem_eff_diag_pass_chi2[diag_test_d]++;
        // Per-detector fit-residual gate on the 3 anchors.
        for (int d = 0; d < n_dets; ++d) {
            if (!a.matched[d]) continue;
            const float px = a.fit.ax + a.fit.bx * a.cand_lz[d];
            const float py = a.fit.ay + a.fit.by * a.cand_lz[d];
            const float dx = a.cand_lx[d] - px;
            const float dy = a.cand_ly[d] - py;
            const float s  = sigmaGem(d);
            const float c  = gem_eff_match_nsigma * s;
            if (dx*dx + dy*dy > c*c) {
                a.matched[0] = a.matched[1] = a.matched[2] = a.matched[3] = false;
                return a;
            }
        }
        if (S < 0 && diag_test_d >= 0 && diag_test_d < 4)
            gem_eff_diag_pass_resid[diag_test_d]++;
        a.valid = true;
        return a;
    };

    // Build the anchor for one test_d under the configured LOO mode.
    auto buildAnchor = [&](int test_d) -> Anchor {
        unsigned cand_mask = 0;
        for (int d = 0; d < n_dets; ++d) {
            if (d != test_d) cand_mask |= (1u << d);
        }
        const bool target_in_fit =
            (gem_eff_loo_mode == GemEffLooMode::LooTargetIn);
        if (gem_eff_loo_mode == GemEffLooMode::TargetSeed) {
            return trySeed(-1, -1, cand_mask, /*target_in_fit=*/false,
                           /*diag_test_d=*/test_d);
        }
        // GEM-seeded modes (loo, loo-target-in): try every (S, hit) in the
        // candidate detectors, keep the lowest-χ² anchor.
        Anchor best{};
        best.fit.chi2_per_dof = std::numeric_limits<float>::infinity();
        for (int S = 0; S < n_dets; ++S) {
            if (S == test_d) continue;
            const auto &hits_s = hits_by_det[S];
            int max_seeds = std::min((int)hits_s.size(),
                                     gem_eff_max_hits_per_det);
            for (int si = 0; si < max_seeds; ++si) {
                Anchor a = trySeed(S, si, cand_mask, target_in_fit,
                                   /*diag_test_d=*/-1);
                if (!a.valid) continue;
                if (a.fit.chi2_per_dof < best.fit.chi2_per_dof) best = a;
            }
        }
        if (!best.valid) best = Anchor{};
        return best;
    };

    // Per-test-detector LOO — run the anchor finder excluding D, then probe D.
    // The snapshot records the LAST valid anchor (whichever test_d succeeds
    // last); per-detector counters accumulate independently.
    GemEffSnapshot &snap = gem_eff_snapshot;
    bool any_valid = false;
    for (int test_d = 0; test_d < n_dets; ++test_d) {
        if (gem_eff_loo_mode == GemEffLooMode::TargetSeed
            && test_d >= 0 && test_d < 4)
            gem_eff_diag_call[test_d]++;
        Anchor a = buildAnchor(test_d);
        if (!a.valid) continue;
        gem_eff_den[test_d]++;
        // Project the anchor fit to test_d, search for a hit within
        // match_nsigma · σ_GEM[test_d].
        float pred_lx, pred_ly;
        projectLineToLocal(gem_transforms[test_d], a.fit, pred_lx, pred_ly);
        int idx; float lab_x, lab_y, lab_z;
        const float s_test = sigmaGem(test_d);
        // Reuse findClosest's gate by temporarily using a tighter window:
        // it normally uses σ_total but here we want σ_GEM only.  Inline
        // the search to keep findClosest unchanged.
        idx = -1;
        const auto &hits_t = hits_by_det[test_d];
        const int max_n = std::min((int)hits_t.size(), gem_eff_max_hits_per_det);
        const float cut_t = gem_eff_match_nsigma * s_test;
        float best_d2 = cut_t * cut_t;
        for (int i = 0; i < max_n; ++i) {
            const auto &h = hits_t[i];
            float lx, ly, lz;
            gem_transforms[test_d].labToLocal(h[0], h[1], h[2], lx, ly, lz);
            float dxr = lx - pred_lx, dyr = ly - pred_ly;
            float d2 = dxr*dxr + dyr*dyr;
            if (d2 < best_d2) {
                best_d2 = d2;
                idx = i;
                lab_x = h[0]; lab_y = h[1]; lab_z = h[2];
            }
        }
        if (idx >= 0) gem_eff_num[test_d]++;
        // Local-coord eff grid: bin the predicted point on test_d's plane
        // over the *active* strip extent (smaller than PlaneConfig.size on
        // the beam-hole side because the inner-edge APV reuses strip numbers
        // via shared_pos — see GemSystem::GetActiveExtent).  Predictions
        // outside the active extent fall in the Histogram2D out-of-range
        // branch (no fill), so the heatmap matches the dashed detector
        // frame drawn in the GUI.
        if (test_d < (int)gem_eff_grid_den.size()
            && test_d < gem_sys.GetNDetectors()) {
            auto xr = gem_sys.GetActiveExtent(test_d, 0);
            auto yr = gem_sys.GetActiveExtent(test_d, 1);
            const float xStep = (xr.second - xr.first) / gem_eff_grid_nx;
            const float yStep = (yr.second - yr.first) / gem_eff_grid_ny;
            gem_eff_grid_den[test_d].fill(pred_lx, pred_ly,
                                           xr.first, xStep, yr.first, yStep);
            if (idx >= 0)
                gem_eff_grid_num[test_d].fill(pred_lx, pred_ly,
                                               xr.first, xStep, yr.first, yStep);
        }

        // Snapshot — record the latest successful LOO test for the GUI.
        snap = GemEffSnapshot{};
        snap.valid = true;
        snap.event_id = event_id;
        snap.hycal_x = hcx; snap.hycal_y = hcy; snap.hycal_z = hcz;
        snap.chi2_per_dof = a.fit.chi2_per_dof;
        snap.ax = a.fit.ax; snap.bx = a.fit.bx;
        snap.ay = a.fit.ay; snap.by = a.fit.by;
        for (int R = 0; R < GEM_EFF_MAX_DETS; ++R) {
            auto &dx = snap.dets[R];
            dx = {};
            if (R >= n_dets) continue;
            // R is "in the anchor fit" iff R != test_d and R was matched.
            dx.used_in_fit = (R != test_d) && a.matched[R];
            dx.hit_present = (R == test_d) ? (idx >= 0) : a.matched[R];
            if (R == test_d && idx >= 0) {
                dx.hit_lab_x = lab_x; dx.hit_lab_y = lab_y; dx.hit_lab_z = lab_z;
            } else if (R != test_d && a.matched[R]) {
                dx.hit_lab_x = a.cand_lx[R];
                dx.hit_lab_y = a.cand_ly[R];
                dx.hit_lab_z = a.cand_lz[R];
            }
            float plx, ply;
            projectLineToLocal(gem_transforms[R], a.fit, plx, ply);
            dx.predicted_local_x = plx;
            dx.predicted_local_y = ply;
            float pX, pY, pZ;
            gem_transforms[R].toLab(plx, ply, pX, pY, pZ);
            dx.predicted_lab_x = pX;
            dx.predicted_lab_y = pY;
            dx.predicted_lab_z = pZ;
            dx.inside = inside(R, plx, ply);
            if (dx.hit_present) {
                float hlx, hly, hlz;
                gem_transforms[R].labToLocal(dx.hit_lab_x, dx.hit_lab_y, dx.hit_lab_z,
                                              hlx, hly, hlz);
                dx.resid_dx = hlx - plx;
                dx.resid_dy = hly - ply;
            }
        }
        any_valid = true;
    }
    if (any_valid) {
        // Closest approach of the fit line to the lab z-axis:
        //   minimize r²(z)=(ax+bx·z)²+(ay+by·z)² → z = -(ax·bx + ay·by)/(bx²+by²)
        // Using the LAST valid LOO anchor's fit as the per-event representative.
        const float bx = snap.bx, by = snap.by;
        const float den = bx*bx + by*by;
        if (den > 1e-12f) {
            const float z_lab    = -(snap.ax*bx + snap.ay*by) / den;
            const float z_offset = z_lab - target_z;
            snap.z_target_lab    = z_lab;
            snap.z_target_offset = z_offset;
            snap.z_target_valid  = true;
            gem_eff_z_target_hist.fill(z_offset,
                                       gem_eff_z_target_min,
                                       gem_eff_z_target_step);
        }
    }
}

nlohmann::json AppState::gemEffSnapshotJson() const
{
    using nlohmann::json;
    const auto &s = gem_eff_snapshot;
    if (!s.valid) return json(nullptr);
    json dets = json::array();
    for (int R = 0; R < GEM_EFF_MAX_DETS; ++R) {
        const auto &d = s.dets[R];
        json e = {
            {"id", R},
            {"used_in_fit", d.used_in_fit},
            {"hit_present", d.hit_present},
            {"inside",      d.inside},
            {"predicted_lab",   {d.predicted_lab_x, d.predicted_lab_y, d.predicted_lab_z}},
            {"predicted_local", {d.predicted_local_x, d.predicted_local_y}},
            {"resid",           {d.resid_dx, d.resid_dy}},
        };
        if (d.hit_present) e["hit_lab"] = json::array({d.hit_lab_x, d.hit_lab_y, d.hit_lab_z});
        dets.push_back(e);
    }
    json out = json{
        {"event_id",     s.event_id},
        {"hycal_lab",    json::array({s.hycal_x, s.hycal_y, s.hycal_z})},
        {"chi2_per_dof", s.chi2_per_dof},
        {"fit",          {{"ax", s.ax}, {"bx", s.bx}, {"ay", s.ay}, {"by", s.by}}},
        {"dets",         dets},
    };
    if (s.z_target_valid) {
        out["z_target_lab"]    = s.z_target_lab;
        out["z_target_offset"] = s.z_target_offset;
    }
    return out;
}

//=============================================================================
// Clearing
//=============================================================================

void AppState::clearHistograms()
{
    std::lock_guard<std::mutex> lk(data_mtx);
    for (auto &[k, h] : histograms)        h.clear();
    for (auto &[k, h] : pos_histograms)   h.clear();
    for (auto &[k, h] : height_histograms) h.clear();
    occupancy.clear();
    occupancy_tcut.clear();
    events_processed = 0;
    cluster_energy_hist.clear();
    nclusters_hist.clear();
    nblocks_hist.clear();
    raw_energy_hist.clear();
    for (auto &h : cluster_energy_hist_by_ncl) h.clear();
    for (auto &h : nblocks_hist_by_ncl)        h.clear();
    energy_angle_hist.clear();
    moller_xy_hist.clear();
    moller_events = 0;
    hycal_xy_hist.clear();
    hycal_xy_events = 0;
    for (auto &h : gem_dx_hist) h.clear();
    for (auto &h : gem_dy_hist) h.clear();
    for (auto &n : gem_match_hits) n = 0;
    gem_match_events = 0;
    clearGemEfficiency();
    cluster_events_processed = 0;
    for (auto &h : gem_occupancy) h.clear();
}

void AppState::clearLms()
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    lms_history.clear();
    latest_lms_integral.clear();
    latest_alpha_integral.clear();
    lms_events = 0;
    lms_first_ts = 0;
    sync_unix = 0;
    sync_rel_sec = 0.;
    pending_sync_unix = 0;
    pending_sync_ti = 0;
}

