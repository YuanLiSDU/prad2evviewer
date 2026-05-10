#pragma once
//=============================================================================
// viewer_utils.h — shared utilities for the event viewer/monitor
//=============================================================================

#include "DetectorTransform.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

#include <fstream>
#include <string>
#include <vector>
#include <cmath>

// --- TI timestamp conversion ------------------------------------------------
// TI clock runs at 250 MHz → 4 ns per tick
static constexpr double TI_TICK_SEC = 4e-9;

// Safe (now − base) seconds.  Both args are uint64 TI ticks; this guards
// against the classic sign-of-life bug:
//   uint64_t a - uint64_t b  with a < b  → wraps to ~2^64
//   ~2^64 × 4 ns = ~73,786,976,288 s, which is what shows up in the
//   EPICS monitor when a snapshot has no TI anchor (timestamp = 0) and
//   the anchor t0 was captured later.
// Treats now == 0 as "no anchor yet" → returns 0 (snapshot stays at the
// origin instead of producing a fake huge offset).  Negative deltas
// (now < base, can happen across an ET reconnect or out-of-order events)
// return as honest negatives instead of wrapping.
inline double ti_delta_sec(uint64_t now, uint64_t base) {
    if (now == 0) return 0.0;
    if (base == 0) return static_cast<double>(now) * TI_TICK_SEC;
    return (now >= base)
        ?  static_cast<double>(now  - base) * TI_TICK_SEC
        : -static_cast<double>(base - now ) * TI_TICK_SEC;
}

// --- file I/O helpers -------------------------------------------------------
inline std::string readFile(const std::string &path) {
    std::ifstream f(path);
    if (!f) return "";
    return {std::istreambuf_iterator<char>(f), {}};
}

inline std::string findFile(const std::string &name, const std::string &base) {
    { std::ifstream f(name); if (f.good()) return name; }
    std::string p = base + "/" + name;
    { std::ifstream f(p); if (f.good()) return p; }
    return "";
}

inline std::string contentType(const std::string &path) {
    if (path.size() >= 5 && path.substr(path.size()-5) == ".html") return "text/html; charset=utf-8";
    if (path.size() >= 4 && path.substr(path.size()-4) == ".css")  return "text/css; charset=utf-8";
    if (path.size() >= 3 && path.substr(path.size()-3) == ".js")   return "application/javascript; charset=utf-8";
    return "application/octet-stream";
}

// --- Histogram (used by both viewer and monitor) ----------------------------
struct Histogram {
    int underflow = 0, overflow = 0;
    std::vector<int> bins;
    void init(int n) { bins.assign(n, 0); underflow = overflow = 0; }
    void fill(float v, float bmin, float bstep) {
        if (v < bmin) { ++underflow; return; }
        int b = (int)((v - bmin) / bstep);
        if (b >= (int)bins.size()) { ++overflow; return; }
        ++bins[b];
    }
    void clear() { std::fill(bins.begin(), bins.end(), 0); underflow = overflow = 0; }
};

struct Histogram2D {
    int nx = 0, ny = 0;
    std::vector<int> bins;  // row-major: bins[iy*nx + ix]
    void init(int nx_, int ny_) { nx = nx_; ny = ny_; bins.assign(nx * ny, 0); }
    void fill(float vx, float vy, float xmin, float xstep, float ymin, float ystep) {
        int ix = (int)((vx - xmin) / xstep);
        int iy = (int)((vy - ymin) / ystep);
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) return;
        bins[iy * nx + ix]++;
    }
    void clear() { std::fill(bins.begin(), bins.end(), 0); }
};

// --- Histogram config -------------------------------------------------------
// Pure binning info — peak-detection thresholds live in WaveConfig
// (daq_config.json fadc250_waveform.analyzer); per-tab Waveform-Tab cuts
// live in PeakFilter (monitor_config.json waveform.filter).
struct HistConfig {
    float bin_min     = 0;
    float bin_max     = 20000;
    float bin_step    = 100;
    float pos_min     = 0;
    float pos_max     = 400;
    float pos_step    = 4;
    float height_min  = 0;
    float height_max  = 4000;
    float height_step = 10;
};

// --- Event-level filters (loaded from external JSON, applied per-event) ------
// Each filter has enable=false by default; disabled filters are skipped.

struct WaveformFilter {
    bool  enable       = false;
    std::vector<std::string> modules;   // HyCal module names; empty = no module restriction
    int   n_peaks_min  = 1;             // qualifying-peak count range
    int   n_peaks_max  = 999999;
    float time_min     = -1e30f;        // peak time range (omit = no cut)
    float time_max     =  1e30f;
    float integral_min = -1e30f;        // peak integral range
    float integral_max =  1e30f;
    float height_min   = -1e30f;        // peak height range
    float height_max   =  1e30f;
};

struct ClusterFilter {
    bool  enable       = false;
    int   n_min        = 0;             // qualifying-cluster count range
    int   n_max        = 999999;
    float energy_min   = 0;             // per-cluster energy range
    float energy_max   = 1e30f;
    int   size_min     = 1;             // per-cluster nblocks range
    int   size_max     = 999999;
    std::vector<std::string> includes_modules;  // cluster must contain >= includes_min of these
    int   includes_min = 1;
    std::vector<std::string> center_modules;    // cluster center must be in this list
};

// --- LMS entry (shared between viewer FileData and monitor globals) ---------
struct LmsEntry {
    double time_sec;    // seconds since first LMS event (from TI timestamp)
    float  integral;    // peak integral within timing cut (or raw ADC for ADC1881M)
};

// --- Peak extraction helpers ------------------------------------------------
// Peaks in `wres` are already gated by the analyzer's
// max(peak_nsigma·ped.rms, min_peak_height) detection cut, so no extra
// height filter is needed in these picks.

// Best peak integral within a time window. Returns -1 if no peak qualifies.
inline float bestPeakInWindow(const fdec::WaveResult &wres,
                               float time_min, float time_max)
{
    float best = -1;
    for (int p = 0; p < wres.npeaks; ++p) {
        auto &pk = wres.peaks[p];
        if (pk.time >= time_min && pk.time <= time_max)
            if (pk.integral > best) best = pk.integral;
    }
    return best;
}

// Best peak integral across all detected peaks (no time cut) — clustering
// input after the per-tab cut decoupling.
inline float bestPeak(const fdec::WaveResult &wres)
{
    float best = -1;
    for (int p = 0; p < wres.npeaks; ++p)
        if (wres.peaks[p].integral > best) best = wres.peaks[p].integral;
    return best;
}
