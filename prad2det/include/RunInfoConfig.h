#pragma once
//=============================================================================
// RunInfoConfig.h — run-period detector geometry / calibration metadata
//
// Shared across the viewer/server, analysis tools, Python bindings, and
// ROOT scripts.  Loads the three-tier runinfo file:
//
//     {
//       "defaults": {                         // invariants — merged under
//         "target": [...], "matching": {...}, // every configurations entry.
//         "gem": {"detectors": [{"id":0,"tilting":[...]}, ...]}
//       },
//       "configurations": [                   // sparse period entries.
//         {"from_run":     0, "beam_energy": 2100.0, ...},
//         {"from_run": 24185, ...}
//       ],
//       "gem_pedestals": [                    // independent per-run table.
//         {"from_run":     0, "pedestal_file": "...", "common_mode_file": "..."},
//         {"from_run": 24023, ...}
//       ]
//     }
//
// Lookup rules:
//   * configurations: overlay every entry with from_run <= run_num in
//     ascending from_run order on top of `defaults` (chained inheritance).
//     Sparse later entries only need to list the fields that differ from
//     the previous period — e.g. a 0.7 GeV entry that reuses the 3.5 GeV
//     survey only specifies beam_energy.  When run_num < 0, all entries
//     are applied (i.e. you get the latest chained state).
//   * gem_pedestals: pick one entry (largest from_run <= run_num),
//     override gem.pedestal_file / gem.common_mode_file on the merged
//     result.  No chaining — each pedestal table row is self-contained.
//
// `run_number` is accepted as an alias for `from_run`, and either of the
// new blocks may be omitted — older runinfo files (no defaults, no
// gem_pedestals; gem.pedestal_file inline in each entry) still load.
//
// Header-only (nlohmann::json + std).  Lives in prad2det/include/ so all
// libraries can pull it in without dragging analysis/ROOT dependencies.
//=============================================================================

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace prad2 {

// Holds all run-specific detector geometry and beam parameters.
// Defaults match the historical 2.2 GeV setup so single-run tools that
// fail to load a runinfo file still produce sensible numbers.
struct RunConfig {
    std::string energy_calib_file;
    float default_adc2mev = 0.12f;
    float Ebeam     = 0.f;
    float target_x  = 0.f;
    float target_y  = 0.f;
    float target_z  = 0.f;
    float hycal_x      = 0.f;
    float hycal_y      = 0.f;
    float hycal_z      = 6225.0f;
    float hycal_tilt_x = 0.f;
    float hycal_tilt_y = 0.f;
    float hycal_tilt_z = 0.f;
    float gem_x[4]      = {-252.8f,  252.8f, -252.8f,  252.8f};
    float gem_y[4]      = {0.f, 0.f, 0.f, 0.f};
    float gem_z[4]      = {5423.0f, 5384.0f, 5823.0f, 5784.0f};
    float gem_tilt_x[4] = {0.f, 0.f, 0.f, 0.f};
    float gem_tilt_y[4] = {0.f, 0.f, 0.f, 0.f};
    float gem_tilt_z[4] = {0.f, 180.f, 0.f, 180.f};
    // Per-run GEM calibration files (paths resolved by the caller against
    // the database dir). Empty -> skip loading (use defaults baked into
    // GemSystem: zero pedestals, permissive CM range).
    std::string gem_pedestal_file;
    std::string gem_common_mode_file;
    // Analysis-specific extras (viewer ignores).
    float hc_time_win_lo = 100.f;  // ns
    float hc_time_win_hi = 200.f;  // ns
    float matching_radius     = 15.f;
    bool  matching_use_square = true;
    bool  matching_energy_dependent = true;
    float matching_sigma = 5.f;
    // For gain correction: which run to use as reference for computing the correction factors.  If negative, use the latest run with gain factors available.
    std::string gain_data_dir = "";
    int gain_ref_run = 23915;
    // Optional per-module HyCal peak-time window file (relative to database
    // dir).  Empty -> uniform [hc_time_win_lo, hc_time_win_hi] from runinfo
    // is used for every module.  See prad2det/include/HyCalTimeCuts.h.
    std::string hycal_time_cut_file;
    // Optional per-module HyCal→RF offset file (relative to database dir).
    // Empty -> uniform 0 ns offset for every module.  See
    // prad2det/include/HyCalRfOffsets.h and
    // docs/analysis_notes/rf_time_reconstruction_plan.md.
    std::string hycal_rf_offset_file;
};

// Returns a RunConfig populated by chaining all matching `configurations`
// entries (and one `gem_pedestals` entry) from `path`.
//
// Selection rule:
//   * configurations: every entry with from_run <= run_num (or every
//     entry when run_num < 0), applied in ascending from_run order.
//   * gem_pedestals: the single entry with the largest from_run <= run_num
//     (latest overall when run_num < 0).
// `run_number` is accepted as an alias for `from_run`.
//
// On any failure (file missing, parse error, no "configurations" array,
// no matching entry) emits a warning and returns the default-constructed
// RunConfig — callers that need to know success/failure should check
// e.g. `result.Ebeam > 0`.
inline RunConfig LoadRunConfig(const std::string &path, int run_num)
{
    RunConfig result;

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Warning: cannot open runinfo file " << path
                  << ", using defaults.\n";
        return result;
    }
    auto cfg = nlohmann::json::parse(f, nullptr, false, true);
    if (cfg.is_discarded()) {
        std::cerr << "Warning: failed to parse " << path
                  << ", using defaults.\n";
        return result;
    }
    if (!cfg.contains("configurations") || !cfg["configurations"].is_array()) {
        std::cerr << "Warning: " << path
                  << " has no \"configurations\" array, using defaults.\n";
        return result;
    }

    // Pull the trigger-run integer off an entry; supports `from_run` (new)
    // and `run_number` (legacy alias).  Returns -1 if neither is present.
    auto entry_from_run = [](const nlohmann::json &e) -> int {
        if (e.contains("from_run"))   return e["from_run"].get<int>();
        if (e.contains("run_number")) return e["run_number"].get<int>();
        return -1;
    };

    // Pick the largest entry with from_run <= run_num (or the largest
    // overall when run_num < 0).  Used by the gem_pedestals lookup.
    auto pick_best = [&](const nlohmann::json &arr) -> std::pair<const nlohmann::json *, int> {
        const nlohmann::json *best = nullptr;
        int best_run = -1;
        for (const auto &e : arr) {
            int rn = entry_from_run(e);
            if (rn < 0) continue;
            if (run_num < 0) {
                if (rn > best_run) { best = &e; best_run = rn; }
            } else if (rn <= run_num && rn > best_run) {
                best = &e; best_run = rn;
            }
        }
        return {best, best_run};
    };

    // Collect all entries with from_run <= run_num (or all entries when
    // run_num < 0), sorted ascending so chained application overlays in
    // the right order.  Used by the configurations chain.
    auto collect_chain = [&](const nlohmann::json &arr)
        -> std::vector<std::pair<int, const nlohmann::json *>>
    {
        std::vector<std::pair<int, const nlohmann::json *>> chain;
        for (const auto &e : arr) {
            int rn = entry_from_run(e);
            if (rn < 0) continue;
            if (run_num >= 0 && rn > run_num) continue;
            chain.emplace_back(rn, &e);
        }
        std::sort(chain.begin(), chain.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
        return chain;
    };

    // Field-by-field overlay.  Called with `defaults` first (if present),
    // then the picked configurations entry — each pass only writes fields
    // it actually contains, so the period entry overlays defaults.
    auto apply_entry = [&](const nlohmann::json &c) {
        if (c.contains("beam_energy")) result.Ebeam = c["beam_energy"].get<float>();
        if (c.contains("calibration")) {
            const auto &cal = c["calibration"];
            if (cal.contains("file"))            result.energy_calib_file = cal["file"].get<std::string>();
            if (cal.contains("default_adc2mev")) result.default_adc2mev   = cal["default_adc2mev"].get<float>();
        }
        if (c.contains("target") && c["target"].is_array() && c["target"].size() >= 3) {
            result.target_x = c["target"][0].get<float>();
            result.target_y = c["target"][1].get<float>();
            result.target_z = c["target"][2].get<float>();
        }
        if (c.contains("hycal")) {
            const auto &h = c["hycal"];
            if (h.contains("position") && h["position"].is_array() && h["position"].size() >= 3) {
                result.hycal_x = h["position"][0].get<float>();
                result.hycal_y = h["position"][1].get<float>();
                result.hycal_z = h["position"][2].get<float>();
            }
            if (h.contains("tilting") && h["tilting"].is_array() && h["tilting"].size() >= 3) {
                result.hycal_tilt_x = h["tilting"][0].get<float>();
                result.hycal_tilt_y = h["tilting"][1].get<float>();
                result.hycal_tilt_z = h["tilting"][2].get<float>();
            }
        }
        if (c.contains("gem") && c["gem"].is_object()) {
            const auto &g = c["gem"];
            if (g.contains("pedestal_file"))    result.gem_pedestal_file    = g["pedestal_file"].get<std::string>();
            if (g.contains("common_mode_file")) result.gem_common_mode_file = g["common_mode_file"].get<std::string>();
            if (g.contains("detectors") && g["detectors"].is_array()) {
                for (const auto &d : g["detectors"]) {
                    if (!d.contains("id")) continue;
                    int id = d["id"].get<int>();
                    if (id < 0 || id >= 4) continue;
                    if (d.contains("position") && d["position"].is_array() && d["position"].size() >= 3) {
                        result.gem_x[id] = d["position"][0].get<float>();
                        result.gem_y[id] = d["position"][1].get<float>();
                        result.gem_z[id] = d["position"][2].get<float>();
                    }
                    if (d.contains("tilting") && d["tilting"].is_array() && d["tilting"].size() >= 3) {
                        result.gem_tilt_x[id] = d["tilting"][0].get<float>();
                        result.gem_tilt_y[id] = d["tilting"][1].get<float>();
                        result.gem_tilt_z[id] = d["tilting"][2].get<float>();
                    }
                }
            }
        }
        if (c.contains("time_cuts")) {
            const auto &tc = c["time_cuts"];
            if (tc.contains("hc_time_window") && tc["hc_time_window"].is_array()
                    && tc["hc_time_window"].size() >= 2) {
                result.hc_time_win_lo = tc["hc_time_window"][0].get<float>();
                result.hc_time_win_hi = tc["hc_time_window"][1].get<float>();
            }
            if (tc.contains("hycal_module_file"))
                result.hycal_time_cut_file = tc["hycal_module_file"].get<std::string>();
            if (tc.contains("hycal_rf_offsets"))
                result.hycal_rf_offset_file = tc["hycal_rf_offsets"].get<std::string>();
        }
        if (c.contains("matching")) {
            const auto &m = c["matching"];
            if (m.contains("radius"))           result.matching_radius           = m["radius"].get<float>();
            if (m.contains("use_square_cut"))   result.matching_use_square       = m["use_square_cut"].get<bool>();
            if (m.contains("energy_dependent")) result.matching_energy_dependent = m["energy_dependent"].get<bool>();
            if (m.contains("sigma"))            result.matching_sigma            = m["sigma"].get<float>();
        }
        if (c.contains("gain_factor") && c["gain_factor"].is_object()) {
            const auto &gf = c["gain_factor"];
            if (gf.contains("data_dir")) result.gain_data_dir = gf["data_dir"].get<std::string>();
            if (gf.contains("ref_run"))  result.gain_ref_run  = gf["ref_run"].get<int>();
        }
    };

    if (run_num < 0) {
        std::cerr << "Warning: unknown run number, picking the entry with "
                     "the largest from_run from " << path << ".\n";
    }

    // 1) defaults first (no-op if absent).
    if (cfg.contains("defaults") && cfg["defaults"].is_object())
        apply_entry(cfg["defaults"]);

    // 2) configurations chain: every entry with from_run <= run_num,
    //    applied in ascending order so the latest overrides earlier ones.
    auto chain = collect_chain(cfg["configurations"]);
    if (chain.empty()) {
        std::cerr << "Warning: no matching configuration in " << path
                  << " for run " << run_num << ", using defaults.\n";
        return result;
    }
    for (const auto &kv : chain) apply_entry(*kv.second);
    int best_run = chain.back().first;

    // 3) gem_pedestals: independent lookup; overrides the gem.*_file fields
    //    so a new pedestal can be added without touching configurations.
    int best_ped_run = -1;
    if (cfg.contains("gem_pedestals") && cfg["gem_pedestals"].is_array()) {
        auto [ped, ped_run] = pick_best(cfg["gem_pedestals"]);
        if (ped != nullptr) {
            if (ped->contains("pedestal_file"))
                result.gem_pedestal_file    = (*ped)["pedestal_file"].get<std::string>();
            if (ped->contains("common_mode_file"))
                result.gem_common_mode_file = (*ped)["common_mode_file"].get<std::string>();
            best_ped_run = ped_run;
        }
    }

    // 4） beam and target positions: independent lookup; overrides the target_* fields
    int best_beam_run = -1;
    if (cfg.contains("beam_target_positions") && cfg["beam_target_positions"].is_array()) {
        auto [btp, btp_run] = pick_best(cfg["beam_target_positions"]);
        if (btp != nullptr && btp->contains("target") && (*btp)["target"].is_array() && (*btp)["target"].size() >= 3) {
            result.target_x = (*btp)["target"][0].get<float>();
            result.target_y = (*btp)["target"][1].get<float>();
            result.target_z = (*btp)["target"][2].get<float>();
            best_beam_run = btp_run;
        }
    }

    std::cerr << "RunInfo   : chained " << chain.size()
              << " config(s), last from_run=" << best_run;
    if (best_ped_run >= 0) std::cerr << "  ped_from_run=" << best_ped_run;
    if (best_beam_run >= 0) std::cerr << "  beam_from_run=" << best_beam_run;
    std::cerr << " from " << path << "\n";

    // If gain_data_dir is empty (not set in JSON, or explicitly set to ""),
    // fall back to <db>/gain_factor derived from the runinfo file location.
    if (result.gain_data_dir.empty()) {
        result.gain_data_dir =
            std::filesystem::path(path).parent_path().parent_path().string()
            + "/gain_factor";
    }

    // Apply the beam and target offset to the detector positions so that the caller gets
    // coordinates in the beam center frame (analysis frame).
    // Detector positions have only one entry for all the runs, in the HyCal center 
    // coordinate frame.  When a run-specific target position is given, we shift the 
    // HyCal and GEM positions to the beam center coordinate frame. 
    result.hycal_x -= result.target_x;
    result.hycal_y -= result.target_y;
    result.hycal_z -= result.target_z;
    for (int i = 0; i < 4; ++i) {
        result.gem_x[i] -= result.target_x;
        result.gem_y[i] -= result.target_y;
        result.gem_z[i] -= result.target_z;
    }

    return result;
}

// Append (or overwrite by run_number) an entry into the "configurations"
// array of `path`. Creates the file if missing. Sorts the array by
// run_number for readability. Atomic-rename via tmp file.
//
// NOTE: this writes a flat, fully-specified entry — it does NOT split into
// the defaults / gem_pedestals tiers.  Round-trips through LoadRunConfig
// correctly (run_number is accepted as a from_run alias, and inline gem
// pedestal paths are still merged), but the resulting file mixes tier
// shapes if it already had a `defaults` / `gem_pedestals` block.  For
// programmatic pedestal table updates, edit `gem_pedestals` directly.
inline bool WriteRunConfig(const std::string &path, int run_num,
                           const RunConfig &geo)
{
    nlohmann::json cfg;
    {
        std::ifstream in(path);
        if (in) {
            cfg = nlohmann::json::parse(in, nullptr, false, true);
            if (cfg.is_discarded()) {
                std::cerr << "Warning: failed to parse " << path
                          << ", overwriting with new data.\n";
                cfg = nlohmann::json::object();
            }
        } else {
            cfg = nlohmann::json::object();
        }
    }

    if (!cfg.contains("configurations") || !cfg["configurations"].is_array())
        cfg["configurations"] = nlohmann::json::array();

    nlohmann::json entry;
    entry["run_number"]                      = run_num;
    entry["beam_energy"]                     = geo.Ebeam;
    entry["calibration"]["file"]             = geo.energy_calib_file;
    entry["calibration"]["default_adc2mev"]  = geo.default_adc2mev;
    entry["target"] = nlohmann::json::array({geo.target_x, geo.target_y, geo.target_z});
    entry["hycal"]["position"] = nlohmann::json::array({geo.hycal_x, geo.hycal_y, geo.hycal_z});
    entry["hycal"]["tilting"]  = nlohmann::json::array({geo.hycal_tilt_x, geo.hycal_tilt_y, geo.hycal_tilt_z});
    {
        nlohmann::json gem = nlohmann::json::object();
        if (!geo.gem_pedestal_file.empty())    gem["pedestal_file"]    = geo.gem_pedestal_file;
        if (!geo.gem_common_mode_file.empty()) gem["common_mode_file"] = geo.gem_common_mode_file;
        nlohmann::json detectors = nlohmann::json::array();
        for (int i = 0; i < 4; ++i) {
            nlohmann::json g;
            g["id"]       = i;
            g["position"] = nlohmann::json::array({geo.gem_x[i], geo.gem_y[i], geo.gem_z[i]});
            g["tilting"]  = nlohmann::json::array({geo.gem_tilt_x[i], geo.gem_tilt_y[i], geo.gem_tilt_z[i]});
            detectors.push_back(g);
        }
        gem["detectors"] = detectors;
        entry["gem"] = gem;
    }
    entry["time_cuts"]["hc_time_window"] = nlohmann::json::array({geo.hc_time_win_lo, geo.hc_time_win_hi});
    if (!geo.hycal_time_cut_file.empty())
        entry["time_cuts"]["hycal_module_file"] = geo.hycal_time_cut_file;
    if (!geo.hycal_rf_offset_file.empty())
        entry["time_cuts"]["hycal_rf_offsets"] = geo.hycal_rf_offset_file;
    entry["matching"]["radius"]          = geo.matching_radius;
    entry["matching"]["use_square_cut"]  = geo.matching_use_square;
    if (!geo.gain_data_dir.empty() || geo.gain_ref_run >= 0) {
        if (!geo.gain_data_dir.empty()) entry["gain_factor"]["data_dir"] = geo.gain_data_dir;
        if (geo.gain_ref_run >= 0)      entry["gain_factor"]["ref_run"]  = geo.gain_ref_run;
    }

    auto &arr = cfg["configurations"];
    bool replaced = false;
    for (auto &e : arr) {
        if (e.contains("run_number") && e["run_number"].get<int>() == run_num) {
            // Merge entry into e field-by-field: existing keys are updated
            // in-place (preserving their original position); new keys are
            // appended at the end.
            for (auto it = entry.begin(); it != entry.end(); ++it)
                e[it.key()] = it.value();
            replaced = true;
            break;
        }
    }
    if (!replaced) arr.push_back(entry);

    std::sort(arr.begin(), arr.end(), [](const nlohmann::json &a, const nlohmann::json &b) {
        int ra = a.contains("run_number") ? a["run_number"].get<int>() : -1;
        int rb = b.contains("run_number") ? b["run_number"].get<int>() : -1;
        return ra < rb;
    });

    std::string tmp = path + ".tmp";
    {
        std::ofstream out(tmp);
        if (!out) {
            std::cerr << "Error: cannot write " << tmp << "\n";
            return false;
        }
        out << cfg.dump(4) << "\n";
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::cerr << "Error: failed to rename " << tmp << " -> " << path << "\n";
        return false;
    }
    std::cerr << (replaced ? "Updated" : "Appended") << " run_number=" << run_num
              << " in " << path << "\n";
    return true;
}

} // namespace prad2
