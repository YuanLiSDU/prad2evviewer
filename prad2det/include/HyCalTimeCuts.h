#pragma once
//=============================================================================
// HyCalTimeCuts.h — per-module HyCal peak-time window with a default
// fallback.
//
// Built once per run by PipelineBuilder (or by analysis tools directly via
// LoadHyCalTimeCuts).  Sized to hycal.module_count() so the per-event
// inner loop can index it by `mod->index` without bounds work in the hot
// path — `at(idx)` for direct access, `in_window(idx, t)` for the common
// "is this peak inside the cut" call.
//
// File format (under database/hycal_time_cut/, e.g. cut_2p1.json):
//
//     {
//       "default": [lo_ns, hi_ns],
//       "modules": [
//         {"name": "G123", "window": [lo, hi]},
//         {"name": "W735", "window": [lo, hi]}
//       ]
//     }
//
// Resolution per module:
//   * Module listed in `modules`  -> per-module window.
//   * Else if `default` in file   -> file's default window.
//   * Else                        -> [def_lo, def_hi] passed in by the caller
//                                    (typically hc_time_window from runinfo).
//
// Keyed by module *name* (not PrimEx ID) to match the established
// calibration_factor_*.json convention.
//=============================================================================

#include "HyCalSystem.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace prad2 {

struct HyCalTimeCuts {
    float default_lo = 100.f;          // ns
    float default_hi = 200.f;          // ns
    std::vector<float> lo;             // per-module-index; size = module_count()
    std::vector<float> hi;
    int   n_overrides = 0;             // count of per-module rows applied

    struct Window { float lo, hi; };

    Window at(int module_index) const {
        if (module_index < 0 || module_index >= static_cast<int>(lo.size()))
            return {default_lo, default_hi};
        return {lo[module_index], hi[module_index]};
    }

    bool in_window(int module_index, float t) const {
        const Window w = at(module_index);
        return t > w.lo && t < w.hi;
    }
};

// Build a HyCalTimeCuts table sized to hycal.module_count().  When `path`
// is empty or unreadable, returns a uniform table at [def_lo, def_hi] so
// every caller gets a usable lookup either way (single code path at the
// call site).
//
// `def_lo`/`def_hi` come from RunConfig::hc_time_win_lo/_hi (i.e. the
// runinfo `time_cuts.hc_time_window`).  A "default" inside the file
// overrides them.
inline HyCalTimeCuts LoadHyCalTimeCuts(const std::string &path,
                                       const fdec::HyCalSystem &hycal,
                                       float def_lo, float def_hi)
{
    HyCalTimeCuts cuts;
    cuts.default_lo = def_lo;
    cuts.default_hi = def_hi;
    const int n = hycal.module_count();
    cuts.lo.assign(n, def_lo);
    cuts.hi.assign(n, def_hi);

    if (path.empty()) return cuts;

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Warning: cannot open HyCal time-cut file " << path
                  << ", using uniform window.\n";
        return cuts;
    }
    auto j = nlohmann::json::parse(f, nullptr, false, true);
    if (j.is_discarded()) {
        std::cerr << "Warning: failed to parse " << path
                  << ", using uniform window.\n";
        return cuts;
    }

    if (j.contains("default") && j["default"].is_array() && j["default"].size() >= 2) {
        float dlo = j["default"][0].get<float>();
        float dhi = j["default"][1].get<float>();
        cuts.default_lo = dlo;
        cuts.default_hi = dhi;
        std::fill(cuts.lo.begin(), cuts.lo.end(), dlo);
        std::fill(cuts.hi.begin(), cuts.hi.end(), dhi);
    }

    int unknown = 0;
    if (j.contains("modules") && j["modules"].is_array()) {
        for (const auto &m : j["modules"]) {
            if (!m.contains("name")) continue;
            if (!m.contains("window") || !m["window"].is_array() || m["window"].size() < 2) continue;
            const auto *mod = hycal.module_by_name(m["name"].get<std::string>());
            if (!mod) { ++unknown; continue; }
            cuts.lo[mod->index] = m["window"][0].get<float>();
            cuts.hi[mod->index] = m["window"][1].get<float>();
            ++cuts.n_overrides;
        }
    }

    std::cerr << "HyCal time: " << cuts.n_overrides << " module overrides";
    if (unknown) std::cerr << " (" << unknown << " unknown names skipped)";
    std::cerr << ", default=[" << cuts.default_lo << ", " << cuts.default_hi
              << "] ns from " << path << "\n";
    return cuts;
}

} // namespace prad2
