#pragma once
//=============================================================================
// PipelineBuilder.h — one-stop wiring for HyCal + GEM detector setup.
//
// Three callers used to hand-wire the same Init / LoadCalibration /
// LoadPedestals / LoadCommonModeRange / SetReconConfigs sequence — analysis
// scripts, the C++ server, and the Python bindings.  Forgetting any step
// (most painfully the GEM crate remap from daq_cfg.roc_tags) silently
// dropped data.  PipelineBuilder consolidates all of it so the wiring is
// impossible to forget.
//
// Boundary: the builder owns "detectors ready" — initialized + calibrated
// HyCalSystem + GemSystem, prepared DetectorTransforms, and the resolved
// configs.  It does NOT own per-event scratch (HyCalCluster, GemCluster,
// WaveAnalyzer) or server-only concerns (monitor_config, trigger filters,
// EPICS, livetime, histograms) — those stay in the caller.
//
// Usage (typical analysis script):
//
//   auto p = prad2::PipelineBuilder()
//       .set_run_number_from_evio(evio_path)
//       .build();
//
//   fdec::HyCalCluster hc(p.hycal);
//   hc.SetConfig(p.hycal_cluster_cfg);
//   gem::GemCluster   gem_cl;
//   fdec::WaveAnalyzer wa(p.daq_cfg.wave_cfg);
//   // ... per-event loop using p.hycal, p.gem, p.hycal_transform, ...
//=============================================================================

#include "DetectorTransform.h"
#include "GemCluster.h"      // gem::ClusterConfig (per-detector knobs)
#include "GemSystem.h"
#include "HyCalCluster.h"    // fdec::ClusterConfig (HyCal clusterer knobs)
#include "HyCalSystem.h"
#include "HyCalTimeCuts.h"   // prad2::HyCalTimeCuts (per-module time window)
#include "HyCalRfOffsets.h"  // prad2::HyCalRfOffsets (per-module RF offset)
#include "RunInfoConfig.h"

#include "DaqConfig.h"       // evc::DaqConfig — pulled in from prad2dec

#include <array>
#include <functional>
#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace prad2 {

// --- result bundle --------------------------------------------------------
//
// Returned by PipelineBuilder::build().  Move-only in practice (HyCalSystem
// and GemSystem hold large internal buffers).  Callers either move this
// into long-lived storage (server: into AppState members) or hold it on
// the stack across the per-event loop (scripts).
struct Pipeline {
    // Inputs as loaded — kept for diagnostics, server status panels, and
    // downstream code that wants to see the resolved values.
    evc::DaqConfig                     daq_cfg;
    prad2::RunConfig                   run_cfg;
    int                                run_number = -1;

    // Detectors — initialized + calibrated + ready per-event.
    fdec::HyCalSystem                  hycal;
    gem::GemSystem                     gem;

    // HyCal cluster config — caller constructs HyCalCluster(p.hycal) and
    // calls SetConfig(p.hycal_cluster_cfg).  GEM per-detector configs are
    // already installed inside `gem` via SetReconConfigs.
    fdec::ClusterConfig                hycal_cluster_cfg;

    // Per-module HyCal peak-time windows.  Always sized to hycal.module_count()
    // when build() succeeds (uniform default when no per-module file).  Use
    // `at(mod->index)` or `in_window(mod->index, t)` in the per-event loop.
    HyCalTimeCuts                      hycal_time_cuts;

    // Per-module HyCal→RF time offsets (ns, folded to (−T_RF/2, T_RF/2]).
    // Always sized to hycal.module_count() (uniform 0.0 when no file).
    // Use `apply(mod->index, folded_dt)` to subtract + re-fold.
    HyCalRfOffsets                     hycal_rf_offsets;

    // Geometry — set() already called, rotation matrix cached.
    DetectorTransform                  hycal_transform;
    std::array<DetectorTransform, 4>   gem_transforms;

    // From reconstruction_config.json:matching.  Defaults reproduce the
    // historical 2.6/sqrt(E_GeV) face sigma and 0.1 mm GEM resolution.
    std::array<float, 3>               hycal_pos_res  = {2.6f, 0.f, 0.f};
    std::vector<float>                 gem_pos_res    = {0.1f, 0.1f, 0.1f, 0.1f};
    std::array<float, 3>               target_pos_res = {1.f, 1.f, 20.f};

    // GEM file-side hardware crate ID -> logical crate ID, derived from
    // daq_cfg.roc_tags entries with type=="gem".  Exposed so callers
    // loading extra GEM-keyed files don't have to re-derive it.
    std::map<int, int>                 gem_crate_remap;

    // Resolved absolute paths (logging + status displays).  Empty when
    // the corresponding step was skipped (e.g. no calibration file).
    std::string                        daq_config_path;
    std::string                        recon_config_path;
    std::string                        runinfo_path;
    std::string                        hycal_map_path;
    std::string                        gem_map_path;
    std::string                        hycal_calib_path;
    std::string                        hycal_time_cut_path;     // optional per-module window file
    std::string                        hycal_rf_offset_path;    // optional per-module RF offset file
    std::string                        gem_pedestal_path;
    std::string                        gem_common_mode_path;

    Pipeline() = default;
    Pipeline(Pipeline &&) = default;
    Pipeline &operator=(Pipeline &&) = default;
    Pipeline(const Pipeline &) = delete;
    Pipeline &operator=(const Pipeline &) = delete;
};

// --- builder --------------------------------------------------------------
//
// Fluent — every setter returns *this so calls chain.  Empty path strings
// fall back to defaults; non-empty paths go through the path resolver
// (default: "{database_dir}/{relative}", or absolute pass-through).
class PipelineBuilder {
public:
    PipelineBuilder() = default;

    // --- path overrides ---------------------------------------------------
    PipelineBuilder &set_database_dir(std::string p);
    PipelineBuilder &set_daq_config(std::string p);

    // Skip the load_daq_config step and use a caller-supplied DaqConfig
    // instead.  Server uses this because it loads daq_cfg early (for
    // legacy ADC pedestals + NNLS template store) before reaching the
    // builder.  When set, set_daq_config(path) is ignored and
    // pipeline.daq_config_path will be empty unless the caller also
    // calls set_daq_config(known_path) for logging purposes.
    PipelineBuilder &set_loaded_daq_config(evc::DaqConfig cfg);

    PipelineBuilder &set_recon_config(std::string p);
    PipelineBuilder &set_runinfo(std::string p);
    PipelineBuilder &set_hycal_map(std::string p);
    PipelineBuilder &set_gem_map(std::string p);
    PipelineBuilder &set_hycal_calib(std::string p);
    PipelineBuilder &set_hycal_time_cut(std::string p);
    PipelineBuilder &set_hycal_rf_offset(std::string p);
    PipelineBuilder &set_gem_pedestal(std::string p);
    PipelineBuilder &set_gem_common_mode(std::string p);

    // --- run number selection --------------------------------------------
    // Default: -1 (LoadRunConfig picks the largest entry).  set_*_from_evio
    // sniffs an integer from "prad_NNNNNN.evio.*" / "run_NNNNNN.*" patterns;
    // -1 if the basename doesn't match.
    PipelineBuilder &set_run_number(int n);
    PipelineBuilder &set_run_number_from_evio(const std::string &evio_path);

    // --- logging ---------------------------------------------------------
    // Default: std::cerr.  Pass nullptr to silence.  Format mirrors the
    // existing [setup] / [GEMSYS] / [GEMCFG] / [PEDSUM] lines so the
    // GEM-efficiency parity audit (Python vs C++ server) keeps diffing clean.
    PipelineBuilder &set_log_stream(std::ostream *s);
    PipelineBuilder &set_log_pedestal_checksum(bool b);

    // --- path resolver ---------------------------------------------------
    // Optional override for resolving relative paths.  Default joins with
    // database_dir (matches analysis scripts' resolve_db_path).  Server
    // plugs in findFile(p, db_dir) for multi-dir search.  Receives any
    // non-absolute path; expected to return an absolute path or "" if not
    // found.
    PipelineBuilder &set_path_resolver(
        std::function<std::string(const std::string &)> resolver);

    // --- build -----------------------------------------------------------
    // Throws std::runtime_error on hard failures (missing daq_config,
    // recon_config, runinfo, hycal_map, gem_map).  Soft failures (missing
    // calib / pedestal / CM) emit warnings via log_stream and leave the
    // corresponding fields empty.
    Pipeline build();

private:
    // Path overrides — empty means "use default" (see header docstring).
    std::string database_dir_;
    std::string daq_config_path_;
    evc::DaqConfig daq_config_loaded_;
    bool        have_loaded_daq_config_ = false;
    std::string recon_config_path_;
    std::string runinfo_path_;
    std::string hycal_map_path_;
    std::string gem_map_path_;
    std::string hycal_calib_path_;
    std::string hycal_time_cut_path_;
    std::string hycal_rf_offset_path_;
    std::string gem_pedestal_path_;
    std::string gem_common_mode_path_;

    int run_number_ = -1;

    std::ostream *log_stream_ = nullptr;          // resolved in build()
    bool          log_stream_default_ = true;     // true → use std::cerr
    bool          log_pedestal_checksum_ = true;

    std::function<std::string(const std::string &)> path_resolver_;
};

} // namespace prad2
