#pragma once
//=============================================================================
// HyCalRfOffsets.h — per-module HyCal→RF time offset table.
//
// Built once per run by PipelineBuilder (or directly by analysis tools via
// LoadHyCalRfOffsets).  Sized to hycal.module_count() so the per-event
// inner loop can index by `mod->index` without bounds work — `at(idx)`
// for direct access, `apply(idx, dt)` for the common "subtract offset
// and re-fold" call.
//
// File format (under database/hycal_rf_offsets/, e.g. 24340.json):
//
//     {
//       "default": 0.0,
//       "modules": [
//         {"name": "G123", "offset_ns":  1.21},
//         {"name": "W735", "offset_ns": -0.74}
//       ]
//     }
//
// Resolution per module:
//   * Module listed in `modules`  -> per-module offset.
//   * Else if `default` in file   -> file's default offset.
//   * Else                        -> `def_off` passed in by the caller
//                                    (typically 0.0).
//
// Keyed by module *name* (not PrimEx ID) to match the calibration_factor_*
// and hycal_time_cut conventions.  Offsets live on (−T_RF/2, T_RF/2] —
// anything outside that range can't be recovered from RF alignment alone
// (see docs/analysis_notes/rf_time_reconstruction_plan.md).
//=============================================================================

#include "HyCalSystem.h"
#include "RfTime.h"          // RF_PERIOD_NS, FoldRfDelta

#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace prad2 {

struct HyCalRfOffsets {
    float default_off = 0.f;           // ns
    std::vector<float> off;            // per-module-index; size = module_count()
    int   n_overrides = 0;             // count of per-module rows applied

    // Offset for a given module index (default if out of range).
    float at(int module_index) const {
        if (module_index < 0 || module_index >= static_cast<int>(off.size()))
            return default_off;
        return off[module_index];
    }

    // Apply per-module offset to an already-folded Δt and re-fold so the
    // result stays on (−T_RF/2, T_RF/2].  NaN passes through unchanged.
    float apply(int module_index, float folded_dt) const {
        return FoldRfDelta(folded_dt - at(module_index));
    }
};

// Build a HyCalRfOffsets table sized to hycal.module_count().  When `path`
// is empty or unreadable, returns a uniform table at `def_off` so callers
// can use a single code path either way.
inline HyCalRfOffsets LoadHyCalRfOffsets(const std::string &path,
                                         const fdec::HyCalSystem &hycal,
                                         float def_off = 0.f)
{
    HyCalRfOffsets cuts;
    cuts.default_off = def_off;
    const int n = hycal.module_count();
    cuts.off.assign(n, def_off);

    if (path.empty()) return cuts;

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Warning: cannot open HyCal RF-offset file " << path
                  << ", using uniform " << def_off << " ns.\n";
        return cuts;
    }
    auto j = nlohmann::json::parse(f, nullptr, false, true);
    if (j.is_discarded()) {
        std::cerr << "Warning: failed to parse " << path
                  << ", using uniform " << def_off << " ns.\n";
        return cuts;
    }

    if (j.contains("default") && j["default"].is_number()) {
        const float dof = j["default"].get<float>();
        cuts.default_off = dof;
        std::fill(cuts.off.begin(), cuts.off.end(), dof);
    }

    int unknown = 0;
    if (j.contains("modules") && j["modules"].is_array()) {
        for (const auto &m : j["modules"]) {
            if (!m.contains("name")) continue;
            if (!m.contains("offset_ns") || !m["offset_ns"].is_number()) continue;
            const auto *mod = hycal.module_by_name(m["name"].get<std::string>());
            if (!mod) { ++unknown; continue; }
            cuts.off[mod->index] = m["offset_ns"].get<float>();
            ++cuts.n_overrides;
        }
    }

    std::cerr << "HyCal RF: " << cuts.n_overrides << " module overrides";
    if (unknown) std::cerr << " (" << unknown << " unknown names skipped)";
    std::cerr << ", default=" << cuts.default_off << " ns from " << path
              << "\n";
    return cuts;
}

} // namespace prad2
