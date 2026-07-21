#include "GemSystem.h"
#include "GemCluster.h"
#include "SspData.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstring>
#include <sstream>

using namespace gem;

//=============================================================================
// Construction / destruction
//=============================================================================

GemSystem::GemSystem() = default;
GemSystem::~GemSystem() = default;

//=============================================================================
// Init — load GEM map from JSON
//=============================================================================

void GemSystem::Init(const std::string &map_file)
{
    std::ifstream f(map_file);
    if (!f.is_open()) {
        std::cerr << "GemSystem::Init: cannot open " << map_file << std::endl;
        return;
    }

    nlohmann::json j;
    try { j = nlohmann::json::parse(f, nullptr, true, true); }
    catch (const nlohmann::json::parse_error &e) {
        std::cerr << "GemSystem::Init: parse error: " << e.what() << std::endl;
        return;
    }

    // --- parse layers/detectors ---
    detectors_.clear();
    if (j.contains("layers")) {
        for (auto &layer : j["layers"]) {
            DetectorConfig det;
            det.id   = layer.value("id", static_cast<int>(detectors_.size()));
            det.name = layer.value("name", "GEM" + std::to_string(det.id));
            det.type = layer.value("type", "PRADGEM");

            det.planes[0].type   = 0; // X
            det.planes[0].n_apvs = layer.value("x_apvs", 12);
            det.planes[0].pitch  = layer.value("x_pitch", 0.4f);
            det.planes[0].size   = det.planes[0].n_apvs * APV_STRIP_SIZE * det.planes[0].pitch;

            det.planes[1].type   = 1; // Y
            det.planes[1].n_apvs = layer.value("y_apvs", 24);
            det.planes[1].pitch  = layer.value("y_pitch", 0.4f);
            det.planes[1].size   = det.planes[1].n_apvs * APV_STRIP_SIZE * det.planes[1].pitch;

            detectors_.push_back(det);
        }
    }

    // --- parse global parameters ---
    apv_channels_    = j.value("apv_channels", 128);
    readout_center_  = j.value("readout_center", 32);
    common_thres_    = j.value("common_mode_threshold", 20.f);
    zerosup_thres_   = j.value("zero_suppression_threshold", 5.f);
    crosstalk_thres_ = j.value("cross_talk_threshold", 8.f);

    // Reconstruction config (clustering + XY matching) is supplied by the
    // application layer via PipelineBuilder/reconstruction_config.json:
    //  - one ClusterConfig per detector via SetReconConfigs()
    //  - strip-level cuts (reject_first/last_timebin, min_peak/sum_adc)
    // Until then each detector starts with library defaults.

    // --- parse APV entries ---
    apvs_.clear();
    apv_map_.clear();
    if (j.contains("apvs")) {
        for (auto &entry : j["apvs"]) {
            ApvConfig apv;
            apv.crate_id    = entry.value("crate", -1);
            apv.mpd_id      = entry.value("mpd", -1);
            apv.adc_ch      = entry.value("adc", -1);
            apv.det_id      = entry.value("det", 0);

            std::string plane_str = entry.value("plane", "X");
            apv.plane_type = (plane_str == "Y" || plane_str == "1") ? 1 : 0;

            apv.orient      = entry.value("orient", 0);
            apv.plane_index = entry.value("pos", 0);
            apv.det_pos     = entry.value("det_pos", 0);

            apv.pin_rotate   = entry.value("pin_rotate", 0);
            apv.shared_pos   = entry.value("shared_pos", -1);
            apv.hybrid_board = entry.value("hybrid_board", true);
            apv.match        = entry.value("match", "");

            int idx = static_cast<int>(apvs_.size());
            apv_map_[packApvKey(apv.crate_id, apv.mpd_id, apv.adc_ch)] = idx;
            apvs_.push_back(apv);
        }
    }

    // --- allocate working data ---
    apv_work_.resize(apvs_.size());
    for (size_t i = 0; i < apvs_.size(); ++i)
        buildStripMap(static_cast<int>(i));

    // --- allocate per-plane and per-detector storage ---
    plane_data_.resize(detectors_.size());
    det_hits_.resize(detectors_.size());
    per_det_cfgs_.assign(detectors_.size(), ClusterConfig{});
}

//=============================================================================
// Helper — translate a hardware crate ID via the optional remap.
//=============================================================================
static inline int remapCrate(int crate, const std::map<int, int> &m)
{
    if (m.empty()) return crate;
    auto it = m.find(crate);
    return (it != m.end()) ? it->second : crate;
}

//=============================================================================
// LoadPedestals — per-strip pedestal from upstream APV-block text format
//
// Each APV block:
//   APV <crate> <slot> <fiber> <adc>            (header line)
//   <strip> <offset> <noise>                     (128 strip lines)
//
// "<slot>" is hardware metadata (physical MPD slot in the VME crate); it
// is ignored for matching — APVs are looked up by (crate, fiber, adc).
//
// crate_remap: hardware crate ID -> logical crate ID expected by
// gem_map.json (e.g. 146 -> 1, 147 -> 2). Empty = identity.
//=============================================================================

void GemSystem::LoadPedestals(const std::string &ped_file,
                              const std::map<int, int> &crate_remap)
{
    std::ifstream f(ped_file);
    if (!f.is_open()) {
        std::cerr << "GemSystem::LoadPedestals: cannot open " << ped_file << std::endl;
        return;
    }

    int n_apvs_loaded = 0;
    int n_apvs_unmapped = 0;
    int cur_idx = -1;       // index into apvs_ for the current APV block
    int strips_in_block = 0;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        // strip leading whitespace and check first non-whitespace char
        size_t firstc = line.find_first_not_of(" \t");
        if (firstc == std::string::npos) continue;
        if (line[firstc] == '#') continue;

        std::istringstream iss(line);
        std::string tok;
        iss >> tok;
        if (tok == "APV") {
            int crate, slot, fiber, adc;
            if (!(iss >> crate >> slot >> fiber >> adc)) {
                std::cerr << "GemSystem::LoadPedestals: malformed APV header: "
                          << line << "\n";
                cur_idx = -1;
                continue;
            }
            cur_idx = FindApvIndex(remapCrate(crate, crate_remap), fiber, adc);
            strips_in_block = 0;
            if (cur_idx < 0) ++n_apvs_unmapped;
            else             ++n_apvs_loaded;
        } else {
            // strip line: "<strip> <offset> <noise>"
            if (cur_idx < 0) continue;   // current APV not in our map
            int strip = std::stoi(tok);
            float offset, noise;
            if (!(iss >> offset >> noise)) continue;
            if (strip < 0 || strip >= APV_STRIP_SIZE) continue;
            apvs_[cur_idx].pedestal[strip].offset = offset;
            apvs_[cur_idx].pedestal[strip].noise  = noise;
            ++strips_in_block;
        }
    }

    std::cerr << "GemSystem::LoadPedestals: " << n_apvs_loaded
              << " APVs loaded";
    if (n_apvs_unmapped)
        std::cerr << " (" << n_apvs_unmapped << " unmapped — wrong crate IDs?)";
    std::cerr << " from " << ped_file << "\n";
}

//=============================================================================
// LoadCommonModeRange — per-APV common mode range from upstream text
//
// Format: <crate> <slot> <fiber> <adc> <cm_min> <cm_max>     (one per line)
//
// Same crate_remap and (crate, fiber, adc) keying as LoadPedestals.
//=============================================================================

void GemSystem::LoadCommonModeRange(const std::string &cm_file,
                                    const std::map<int, int> &crate_remap)
{
    std::ifstream f(cm_file);
    if (!f.is_open()) {
        std::cerr << "GemSystem::LoadCommonModeRange: cannot open " << cm_file << std::endl;
        return;
    }

    int n_loaded = 0;
    int n_unmapped = 0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        size_t firstc = line.find_first_not_of(" \t");
        if (firstc == std::string::npos || line[firstc] == '#') continue;

        std::istringstream iss(line);
        int crate, slot, fiber, adc;
        float cm_min, cm_max;
        if (!(iss >> crate >> slot >> fiber >> adc >> cm_min >> cm_max)) continue;

        int idx = FindApvIndex(remapCrate(crate, crate_remap), fiber, adc);
        if (idx < 0) { ++n_unmapped; continue; }

        apvs_[idx].cm_range_min = cm_min;
        apvs_[idx].cm_range_max = cm_max;
        ++n_loaded;
    }
    std::cerr << "GemSystem::LoadCommonModeRange: " << n_loaded
              << " APVs loaded";
    if (n_unmapped)
        std::cerr << " (" << n_unmapped << " unmapped — wrong crate IDs?)";
    std::cerr << " from " << cm_file << "\n";
}

//=============================================================================
// GetHoleXOffset — beam hole X offset from detector center
//=============================================================================

float GemSystem::GetHoleXOffset() const
{
    if (detectors_.empty()) return 0.f;
    // scan match APVs on detector 0, X plane — collect their mapped strip numbers
    int ref_det = detectors_[0].id;
    float pitch = detectors_[0].planes[0].pitch;
    float size  = detectors_[0].planes[0].size;
    int smin = INT_MAX, smax = INT_MIN;
    for (size_t i = 0; i < apvs_.size(); ++i) {
        auto &a = apvs_[i];
        if (a.det_id != ref_det || a.plane_type != 0 || a.match.empty()) continue;
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
            int s = apv_work_[i].strip_map[ch];
            if (s >= 0) { smin = std::min(smin, s); smax = std::max(smax, s); }
        }
    }
    if (smin > smax) return 0.f;
    // hole center in detector-local coords (centered on detector midpoint)
    float hole_center = (smin + smax + 1) * 0.5f * pitch;
    return hole_center - size * 0.5f;
}

//=============================================================================
// GetActiveExtent — true bounding box of mapped strips in local coords
//=============================================================================
//
// PlaneConfig.size = n_apvs * APV_STRIP_SIZE * pitch is the bounding box
// assumption.  When the inner-edge APV uses `shared_pos` (e.g. pos=11 with
// shared_pos=10 on PRAD GEMs), it doesn't extend the strip range — one APV's
// worth of bbox stays empty.  Visualization code wants the tight extent so
// the dashed detector frame and the eff/occupancy bins line up with the
// real readout instead of leaving a margin at the beam-hole side.
//
// Strip-to-local conversion follows the same convention as
// GemSystem::ProcessApv (line 599):
//     pos(strip) = strip * pitch - size/2 + pitch/2
// so the returned extent is the LEFT/BOTTOM edge of the lowest strip and
// the RIGHT/TOP edge of the highest strip (i.e. (smin)*pitch - size/2 and
// (smax+1)*pitch - size/2).

std::pair<float, float> GemSystem::GetActiveExtent(int det_id, int plane) const
{
    if (plane != 0 && plane != 1) return {0.f, 0.f};
    if (det_id < 0 || det_id >= (int)detectors_.size())
        return {0.f, 0.f};
    const float pitch = detectors_[det_id].planes[plane].pitch;
    const float size  = detectors_[det_id].planes[plane].size;
    const float full_min = -size * 0.5f;
    const float full_max =  size * 0.5f;
    int smin = INT_MAX, smax = INT_MIN;
    for (size_t i = 0; i < apvs_.size(); ++i) {
        const auto &a = apvs_[i];
        if (a.det_id != detectors_[det_id].id) continue;
        if (a.plane_type != plane) continue;
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
            int s = apv_work_[i].strip_map[ch];
            if (s >= 0) { smin = std::min(smin, s); smax = std::max(smax, s); }
        }
    }
    if (smin > smax) return {full_min, full_max};
    const float lo = smin * pitch + full_min;
    const float hi = (smax + 1) * pitch + full_min;
    return {lo, hi};
}

//=============================================================================
// Clear — reset per-event data
//=============================================================================

void GemSystem::Clear()
{
    for (auto &w : apv_work_) {
        std::memset(w.hit_pos, 0, sizeof(w.hit_pos));
        // Also zero the raw waveform buffer.  Strict reading of the ZS
        // path says raw is invisible to consumers when hit_pos is false,
        // but in practice not zeroing it makes the per-event result
        // depend on which prior events were processed (see commit msg
        // for the gem_eff_audit/server divergence we chased down).
        std::memset(w.raw, 0, sizeof(w.raw));
    }
    for (auto &pd : plane_data_) {
        pd[0].hits.clear();
        pd[0].clusters.clear();
        pd[1].hits.clear();
        pd[1].clusters.clear();
    }
    for (auto &dh : det_hits_)
        dh.clear();
    all_hits_.clear();
}

//=============================================================================
// ProcessEvent — decode SSP data → strip hits
//=============================================================================

void GemSystem::ProcessEvent(const ssp::SspEventData &evt)
{
    for (int mi = 0; mi < evt.nmpds; ++mi) {
        auto &mpd = evt.mpds[mi];
        if (!mpd.present) continue;

        for (int ai = 0; ai < ssp::MAX_APVS_PER_MPD; ++ai) {
            auto &apv = mpd.apvs[ai];
            if (!apv.present) continue;

            int idx = FindApvIndex(apv.addr.crate_id, apv.addr.mpd_id, apv.addr.adc_ch);
            if (idx < 0) continue;

            processApv(idx, apv);
        }
    }
}

//=============================================================================
// SetReconConfigs — install per-detector clustering / XY-matching params.
// Clamps to detectors_.size(); pads short input with library defaults so
// callers can supply any number of entries (typical: one per detector).
//=============================================================================

void GemSystem::SetReconConfigs(std::vector<ClusterConfig> cfgs)
{
    cfgs.resize(detectors_.size(), ClusterConfig{});
    per_det_cfgs_ = std::move(cfgs);
}

//=============================================================================
// Reconstruct — run clustering on all planes, then 2D matching
//=============================================================================

void GemSystem::Reconstruct(GemCluster &clusterer)
{
    for (int d = 0; d < static_cast<int>(detectors_.size()); ++d) {
        // apply this detector's clustering / matching config
        clusterer.SetConfig(per_det_cfgs_[d]);

        // cluster X and Y planes
        for (int p = 0; p < 2; ++p) {
            auto &pd = plane_data_[d][p];
            clusterer.FormClusters(pd.hits, pd.clusters);
        }

        // Cartesian reconstruction: match X and Y clusters
        auto &xc = plane_data_[d][0].clusters;
        auto &yc = plane_data_[d][1].clusters;
        clusterer.CartesianReconstruct(xc, yc, det_hits_[d], d);

        // accumulate all hits
        all_hits_.insert(all_hits_.end(), det_hits_[d].begin(), det_hits_[d].end());
    }
}

//=============================================================================
// FindApvIndex — O(1) lookup
//=============================================================================

int GemSystem::FindApvIndex(int crate, int mpd, int adc) const
{
    auto it = apv_map_.find(packApvKey(crate, mpd, adc));
    return (it != apv_map_.end()) ? it->second : -1;
}

//=============================================================================
// Accessors
//=============================================================================

const std::vector<StripHit>& GemSystem::GetPlaneHits(int det, int plane) const
{
    return plane_data_[det][plane].hits;
}

const std::vector<StripCluster>& GemSystem::GetPlaneClusters(int det, int plane) const
{
    return plane_data_[det][plane].clusters;
}

const std::vector<GEMHit>& GemSystem::GetHits(int det) const
{
    return det_hits_[det];
}

//=============================================================================
// processApv — per-APV: pedestal subtraction, common mode, zero suppression
//
// Two paths, picked per-APV from the data itself:
//
//   data.nstrips < APV_STRIP_SIZE → firmware zero-suppressed.  Only the
//     surviving strips are in the bank; values are pedestal + CM subtracted
//     by firmware.  Re-apply pedestal.noise × zerosup_thres_ on the
//     surviving strips so the viewer's threshold tracks the offline
//     pedestals even if the firmware threshold drifts.
//
//   data.nstrips == APV_STRIP_SIZE → full readout.  Every channel is in the
//     bank (calibration / debug mode, possibly with firmware emitting CM
//     debug headers but without actually applying ZS).  Run the full
//     offline pipeline: pedestal subtract → sorting common-mode → per-strip
//     ZS with pedestal.noise × zerosup_thres_.
//
// NOTE: we do NOT use `data.has_online_cm` as the discriminator.  The MPD
// firmware can emit type-0xD debug-header words (which set has_online_cm)
// while still sending all 128 strips raw — i.e. CM computation was done
// but ZS was not applied.  `nstrips` is the only signal that actually
// tells us whether the data coming in is zero-suppressed.
//=============================================================================

// Flat-buffer index: raw[ts * APV_STRIP_SIZE + ch].
#define RAW_IDX(ch, ts) ((ts) * APV_STRIP_SIZE + (ch))

void GemSystem::processApv(int apv_idx, const ssp::ApvData &data)
{
    auto &cfg = apvs_[apv_idx];
    auto &work = apv_work_[apv_idx];

    if (data.nstrips < APV_STRIP_SIZE) {
        // Online-ZS path: firmware already pedestal + CM subtracted.  Apply
        // a software N-sigma cut on the surviving strips so absent firmware
        // pedestals can't leak sub-threshold strips into reconstruction.
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
            if (!data.hasStrip(ch)) {
                work.hit_pos[ch] = false;
                continue;
            }
            float avg = 0.f;
            for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
                float v = static_cast<float>(data.strips[ch][ts]);
                work.raw[RAW_IDX(ch, ts)] = v;
                avg += v;
            }
            avg /= SSP_TIME_SAMPLES;
            work.hit_pos[ch] = (avg > cfg.pedestal[ch].noise * zerosup_thres_);
        }
        collectHits(apv_idx);
        return;
    }

    // --- offline pipeline for full-readout (pedestal-calibration) runs ---

    // --- copy raw data into working buffer ---
    for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
        for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
            work.raw[RAW_IDX(ch, ts)] = static_cast<float>(data.strips[ch][ts]);
        }
    }

    // --- common mode correction for each time sample ---
    for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
        float *buf = &work.raw[ts * APV_STRIP_SIZE];

        // subtract pedestal offset
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch)
            buf[ch] -= cfg.pedestal[ch].offset;

        // compute and subtract common mode (sorting algorithm)
        float cm = commonModeSorting(buf, APV_STRIP_SIZE, apv_idx);
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch)
            buf[ch] -= cm;
    }

    // --- zero suppression ---
    for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
        float avg = 0.f;
        for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts)
            avg += work.raw[RAW_IDX(ch, ts)];
        avg /= SSP_TIME_SAMPLES;

        work.hit_pos[ch] = (avg > cfg.pedestal[ch].noise * zerosup_thres_);
    }

    // --- collect hits to plane ---
    collectHits(apv_idx);
}

//=============================================================================
// commonModeSorting — MPD version: remove top N high-ADC strips from average
//=============================================================================

float GemSystem::commonModeSorting(float *buf, int size, [[maybe_unused]] int apv_idx)
{
    float sum = 0.f;
    int count = 0;

    // Track the top NUM_HIGH_STRIPS highest values to exclude
    // Use a stack array to avoid heap allocation in this hot-path function.
    float high_vals[NUM_HIGH_STRIPS];
    std::fill(high_vals, high_vals + NUM_HIGH_STRIPS, -9999.f);

    for (int i = 0; i < size; ++i) {
        sum += buf[i];
        count++;

        // Maintain sorted list of top N highest values
        if (buf[i] > high_vals[0]) {
            // Find insertion point and shift
            int pos = 0;
            while (pos < NUM_HIGH_STRIPS - 1 && high_vals[pos + 1] < buf[i])
                pos++;
            for (int j = 0; j < pos; ++j)
                high_vals[j] = high_vals[j + 1];
            high_vals[pos] = buf[i];
        }
    }

    // Subtract the top N values from sum
    for (int i = 0; i < NUM_HIGH_STRIPS && count > 1; ++i) {
        sum -= high_vals[i];
        count--;
    }

    return (count > 0) ? sum / static_cast<float>(count) : 0.f;
}

//=============================================================================
// commonModeDanning — Danning algorithm with common mode range
//=============================================================================

float GemSystem::commonModeDanning(float *buf, int size, int apv_idx)
{
    auto &cfg = apvs_[apv_idx];

    // Step 1: average A — only values within common mode range
    float avgA = 0.f;
    int countA = 0;
    for (int i = 0; i < size; ++i) {
        if (buf[i] >= cfg.cm_range_min && buf[i] <= cfg.cm_range_max) {
            avgA += buf[i];
            countA++;
        }
    }
    if (countA == 0) return 0.f;
    avgA /= static_cast<float>(countA);

    // Step 2: average B — values below avgA + RMS_THRESHOLD * noise
    static constexpr float RMS_THRESHOLD = 3.f;
    float avgB = 0.f;
    int countB = 0;
    for (int i = 0; i < size; ++i) {
        if (buf[i] < avgA + RMS_THRESHOLD * cfg.pedestal[i].noise) {
            avgB += buf[i];
            countB++;
        }
    }

    return (countB > 0) ? avgB / static_cast<float>(countB) : 0.f;
}

//=============================================================================
// collectHits — gather zero-suppressed hits into plane data
//=============================================================================

void GemSystem::collectHits(int apv_idx)
{
    auto &cfg = apvs_[apv_idx];
    auto &work = apv_work_[apv_idx];

    if (cfg.det_id < 0 || cfg.det_id >= static_cast<int>(detectors_.size()))
        return;
    if (cfg.plane_type < 0 || cfg.plane_type > 1)
        return;

    auto &det = detectors_[cfg.det_id];
    auto &plane = det.planes[cfg.plane_type];
    auto &hits = plane_data_[cfg.det_id][cfg.plane_type].hits;

    for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
        if (!work.hit_pos[ch]) continue;

        int plane_strip = work.strip_map[ch];
        if (plane_strip < 0) continue;

        // Find max charge and max timebin
        float max_charge = -1e9f;
        float sum_adc = 0.f;
        short max_tb = 0;
        // Use a stack array to avoid heap allocation in this hot-path loop.
        float ts_adc[SSP_TIME_SAMPLES];
        for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
            float val = work.raw[RAW_IDX(ch, ts)];
            ts_adc[ts] = val;
            sum_adc += val;
            if (val > max_charge) {
                max_charge = val;
                max_tb = static_cast<short>(ts);
            }
        }

        // Strip-level cuts
        if (reject_first_timebin_ && max_tb == 0) continue;
        if (reject_last_timebin_  && max_tb == SSP_TIME_SAMPLES - 1) continue;
        if (min_peak_adc_ > 0.f   && max_charge < min_peak_adc_) continue;
        if (min_sum_adc_  > 0.f   && sum_adc < min_sum_adc_) continue;

        // Calculate physical position
        float pos = static_cast<float>(plane_strip) * plane.pitch
                    - plane.size * 0.5f + plane.pitch * 0.5f;

        // Check cross-talk
        bool xtalk = (max_charge < cfg.pedestal[ch].noise * crosstalk_thres_)
                     && (max_charge > cfg.pedestal[ch].noise * zerosup_thres_);

        StripHit hit;
        hit.strip       = plane_strip;
        hit.charge      = max_charge;
        hit.max_timebin = max_tb;
        hit.position    = pos;
        hit.cross_talk  = xtalk;
        hit.ts_adc.assign(ts_adc, ts_adc + SSP_TIME_SAMPLES);
        hits.push_back(std::move(hit));
    }
}

#undef RAW_IDX

//=============================================================================
// buildStripMap — compute APV channel → plane strip mapping
//
// Implements the full MapStrip pipeline (from PRadAnalyzer/mpd_gem_view_ssp):
//   1. APV25 internal channel mapping (chip wiring, universal)
//   2. Hybrid board pin conversion (MPD electronics only)
//   3. Readout strip scaling (configurable offset: 32 normal, 48 for special APVs)
//   4. 7-bit mask
//   5. Orient flip
//   6. Plane-wide strip number with configurable offset
//=============================================================================

// Public, stateless — declared in GemSystem.h.  buildStripMap() delegates
// here so on-line reconstruction and off-line analyses share one impl.
namespace gem {

int MapStrip(int ch, int plane_index, int orient,
             int pin_rotate, int shared_pos, bool hybrid_board,
             int apv_channels, int readout_center)
{
    const int N = apv_channels;
    const int readout_off = readout_center + pin_rotate;
    const int eff_pos = (shared_pos >= 0) ? shared_pos : plane_index;
    const int plane_shift = (eff_pos - plane_index) * N - pin_rotate;

    // Step 1: APV25 internal channel mapping (chip wiring, universal).
    int strip = 32 * (ch % 4) + 8 * (ch / 4) - 31 * (ch / 16);

    // Step 2: hybrid board pin conversion (MPD electronics).
    if (hybrid_board)
        strip = strip + 1 + strip % 4 - 5 * ((strip / 4) % 2);

    // Step 3: readout strip mapping (odd/even fan-out around center).
    // readout_off <= 0 → skip (steps 1+2 already give the final strip).
    if (readout_off > 0) {
        if (strip & 1)
            strip = readout_off - (strip + 1) / 2;
        else
            strip = readout_off + strip / 2;
    }

    // Step 4: channel mask.
    strip &= (N - 1);

    // Step 5: orient flip.
    if (orient == 1)
        strip = (N - 1) - strip;

    // Step 6: plane-wide strip number.
    strip += plane_shift + plane_index * N;

    return strip;
}

std::vector<int> MapApvStrips(int plane_index, int orient,
                              int pin_rotate, int shared_pos,
                              bool hybrid_board,
                              int apv_channels, int readout_center)
{
    std::vector<int> out(apv_channels);
    for (int ch = 0; ch < apv_channels; ++ch) {
        out[ch] = MapStrip(ch, plane_index, orient,
                           pin_rotate, shared_pos, hybrid_board,
                           apv_channels, readout_center);
    }
    return out;
}

} // namespace gem

void GemSystem::buildStripMap(int apv_idx)
{
    auto &cfg  = apvs_[apv_idx];
    auto &work = apv_work_[apv_idx];
    const int N = apv_channels_;

    for (int ch = 0; ch < N; ++ch) {
        work.strip_map[ch] = MapStrip(ch, cfg.plane_index, cfg.orient,
                                      cfg.pin_rotate, cfg.shared_pos,
                                      cfg.hybrid_board, N, readout_center_);
    }
}
