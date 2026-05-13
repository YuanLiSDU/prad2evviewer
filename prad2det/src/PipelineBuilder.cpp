#include "PipelineBuilder.h"

#include "load_daq_config.h"

#include <nlohmann/json.hpp>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <utility>

using nlohmann::json;

namespace prad2 {

// =========================================================================
// internal helpers
// =========================================================================

namespace {

bool is_absolute_path(const std::string &p)
{
    if (p.empty()) return false;
    if (p[0] == '/' || p[0] == '\\') return true;
    if (p.size() >= 2 && p[1] == ':') return true;  // Windows drive letter
    return false;
}

int sniff_run_from_basename(const std::string &path)
{
    if (path.empty()) return -1;
    // Strip directory.
    auto slash = path.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? path
                                                    : path.substr(slash + 1);
    static const std::regex pat(R"((?:prad|run)_0*(\d+))",
                                std::regex_constants::icase);
    std::smatch m;
    if (std::regex_search(base, m, pat)) {
        try { return std::stoi(m[1].str()); } catch (...) {}
    }
    return -1;
}

json parse_json_file(const std::string &path)
{
    std::ifstream f(path);
    if (!f) return json::object();
    auto j = json::parse(f, nullptr, /*allow_exceptions=*/false,
                         /*ignore_comments=*/true);
    return j.is_discarded() ? json::object() : j;
}

void apply_gem_cluster_overrides(const json &j, gem::ClusterConfig &cfg)
{
    if (j.contains("min_cluster_hits"))    cfg.min_cluster_hits    = j["min_cluster_hits"];
    if (j.contains("max_cluster_hits"))    cfg.max_cluster_hits    = j["max_cluster_hits"];
    if (j.contains("consecutive_thres"))   cfg.consecutive_thres   = j["consecutive_thres"];
    if (j.contains("split_thres"))         cfg.split_thres         = j["split_thres"];
    if (j.contains("cross_talk_width"))    cfg.cross_talk_width    = j["cross_talk_width"];
    if (j.contains("charac_dists") && j["charac_dists"].is_array()) {
        cfg.charac_dists.clear();
        for (auto &v : j["charac_dists"])
            cfg.charac_dists.push_back(v.get<float>());
    }
    if (j.contains("match_mode"))          cfg.match_mode          = j["match_mode"];
    if (j.contains("match_adc_asymmetry")) cfg.match_adc_asymmetry = j["match_adc_asymmetry"];
    if (j.contains("match_time_diff"))     cfg.match_time_diff     = j["match_time_diff"];
    if (j.contains("match_ts_period"))     cfg.ts_period           = j["match_ts_period"];
}

void apply_hycal_cluster_overrides(const json &j, fdec::ClusterConfig &cfg)
{
    if (j.contains("min_module_energy"))  cfg.min_module_energy  = j["min_module_energy"];
    if (j.contains("min_center_energy"))  cfg.min_center_energy  = j["min_center_energy"];
    if (j.contains("min_cluster_energy")) cfg.min_cluster_energy = j["min_cluster_energy"];
    if (j.contains("min_cluster_size"))   cfg.min_cluster_size   = j["min_cluster_size"];
    if (j.contains("corner_conn"))        cfg.corner_conn        = j["corner_conn"];
    if (j.contains("split_iter"))         cfg.split_iter         = j["split_iter"];
    if (j.contains("least_split"))        cfg.least_split        = j["least_split"];
    if (j.contains("log_weight_thres"))   cfg.log_weight_thres   = j["log_weight_thres"];
    if (j.contains("seed_time_window"))   cfg.seed_time_window   = j["seed_time_window"];
}

} // namespace

// =========================================================================
// fluent setters
// =========================================================================

PipelineBuilder &PipelineBuilder::set_database_dir(std::string p)      { database_dir_         = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_daq_config(std::string p)        { daq_config_path_      = std::move(p); return *this; }

PipelineBuilder &PipelineBuilder::set_loaded_daq_config(evc::DaqConfig cfg)
{
    daq_config_loaded_      = std::move(cfg);
    have_loaded_daq_config_ = true;
    return *this;
}
PipelineBuilder &PipelineBuilder::set_recon_config(std::string p)      { recon_config_path_    = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_runinfo(std::string p)           { runinfo_path_         = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_map(std::string p)         { hycal_map_path_       = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_gem_map(std::string p)           { gem_map_path_         = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_calib(std::string p)       { hycal_calib_path_     = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_time_cut(std::string p)    { hycal_time_cut_path_  = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_gem_pedestal(std::string p)      { gem_pedestal_path_    = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_gem_common_mode(std::string p)   { gem_common_mode_path_ = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_run_number(int n)                { run_number_ = n; return *this; }

PipelineBuilder &PipelineBuilder::set_run_number_from_evio(const std::string &p)
{
    int n = sniff_run_from_basename(p);
    if (n > 0) run_number_ = n;
    return *this;
}

PipelineBuilder &PipelineBuilder::set_log_stream(std::ostream *s)
{
    log_stream_         = s;
    log_stream_default_ = false;
    return *this;
}

PipelineBuilder &PipelineBuilder::set_log_pedestal_checksum(bool b)
{
    log_pedestal_checksum_ = b;
    return *this;
}

PipelineBuilder &PipelineBuilder::set_path_resolver(
    std::function<std::string(const std::string &)> resolver)
{
    path_resolver_ = std::move(resolver);
    return *this;
}

// =========================================================================
// build()
// =========================================================================

Pipeline PipelineBuilder::build()
{
    // --- resolve effective log target -------------------------------------
    std::ostream *log = log_stream_default_ ? &std::cerr : log_stream_;
    auto LOG = [log](const std::string &line) {
        if (log) (*log) << line << '\n';
    };

    // --- resolve database dir + path resolver -----------------------------
    std::string db_dir = database_dir_;
    if (db_dir.empty()) {
        if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))
            db_dir = env;
        else
            db_dir = "database";
    }
    auto resolve = [&](const std::string &p) -> std::string {
        if (p.empty()) return p;
        if (is_absolute_path(p)) return p;
        if (path_resolver_) return path_resolver_(p);
        return db_dir + "/" + p;
    };

    Pipeline out;

    // --- 1. DAQ config ----------------------------------------------------
    if (have_loaded_daq_config_) {
        out.daq_cfg = std::move(daq_config_loaded_);
        // If the caller knows the source path, surface it; otherwise leave
        // empty so logs make it obvious nothing was loaded by the builder.
        out.daq_config_path = daq_config_path_;
        std::ostringstream oss;
        oss << "[setup] DAQ config : (caller-supplied)"
            << (daq_config_path_.empty() ? "" : " " + daq_config_path_);
        LOG(oss.str());
    } else {
        std::string daq_path = daq_config_path_.empty()
            ? resolve("daq_config.json")
            : resolve(daq_config_path_);
        if (daq_path.empty() || !evc::load_daq_config(daq_path, out.daq_cfg)) {
            throw std::runtime_error(
                "PipelineBuilder: cannot load DAQ config '" + daq_path + "'");
        }
        out.daq_config_path = daq_path;
        std::ostringstream oss;
        oss << "[setup] DAQ config : " << daq_path;
        LOG(oss.str());
    }

    // --- 2. GEM crate remap from daq_cfg.roc_tags ------------------------
    for (const auto &re : out.daq_cfg.roc_tags)
        if (re.type == "gem")
            out.gem_crate_remap[(int)re.tag] = re.crate;

    // --- 3. recon config (soft — falls back to library defaults) ----------
    std::string recon_path = recon_config_path_.empty()
        ? resolve("reconstruction_config.json")
        : resolve(recon_config_path_);
    json recon = parse_json_file(recon_path);
    if (recon.empty()) {
        std::ostringstream oss;
        oss << "[WARN] PipelineBuilder: cannot load reconstruction_config '"
            << recon_path << "' — proceeding with library defaults.";
        LOG(oss.str());
    } else {
        out.recon_config_path = recon_path;
    }

    // --- 4. runinfo (soft — defaults to RunConfig{} if absent) ------------
    std::string ri_path = runinfo_path_;
    if (ri_path.empty() && recon.contains("runinfo")
            && recon["runinfo"].is_string())
        ri_path = recon["runinfo"].get<std::string>();
    ri_path = resolve(ri_path);

    out.run_number = run_number_;
    if (!ri_path.empty()) {
        if (run_number_ > 0) {
            std::ostringstream oss;
            oss << "[setup] Run number : " << run_number_;
            LOG(oss.str());
        } else {
            LOG("[setup] Run number : (latest entry — no run-number override)");
        }
        out.run_cfg      = prad2::LoadRunConfig(ri_path, run_number_);
        out.runinfo_path = ri_path;
        std::ostringstream oss;
        oss << "[setup] RunInfo    : " << ri_path
            << "  beam=" << static_cast<int>(out.run_cfg.Ebeam)
            << " MeV  hycal_z=" << std::fixed << std::setprecision(1)
            << out.run_cfg.hycal_z << " mm";
        LOG(oss.str());
    } else {
        LOG("[WARN] PipelineBuilder: no runinfo path resolved — using default "
            "RunConfig (geometry/calibration paths empty).");
    }

    // --- 5. detector transforms ------------------------------------------
    out.hycal_transform.set(
        out.run_cfg.hycal_x, out.run_cfg.hycal_y, out.run_cfg.hycal_z,
        out.run_cfg.hycal_tilt_x, out.run_cfg.hycal_tilt_y, out.run_cfg.hycal_tilt_z);
    for (int d = 0; d < 4; ++d) {
        out.gem_transforms[d].set(
            out.run_cfg.gem_x[d], out.run_cfg.gem_y[d], out.run_cfg.gem_z[d],
            out.run_cfg.gem_tilt_x[d], out.run_cfg.gem_tilt_y[d], out.run_cfg.gem_tilt_z[d]);
    }

    // --- 6. HyCal map (soft — leaves hycal default-constructed on failure) -
    std::string hc_map = hycal_map_path_.empty()
        ? resolve("hycal_map.json")
        : resolve(hycal_map_path_);
    if (hc_map.empty()) {
        LOG("[WARN] PipelineBuilder: no HyCal map resolved.");
    } else if (!out.hycal.Init(hc_map)) {
        std::ostringstream oss;
        oss << "[WARN] PipelineBuilder: HyCalSystem.Init failed for '"
            << hc_map << "'.";
        LOG(oss.str());
    } else {
        out.hycal_map_path = hc_map;
    }

    // --- 7. HyCal calibration --------------------------------------------
    std::string hc_calib = hycal_calib_path_.empty()
        ? resolve(out.run_cfg.energy_calib_file)
        : resolve(hycal_calib_path_);
    if (!hc_calib.empty()) {
        int n = out.hycal.LoadCalibration(hc_calib);
        out.hycal_calib_path = hc_calib;
        std::ostringstream oss;
        oss << "[setup] HC calib   : " << hc_calib << " (" << n << " modules)";
        LOG(oss.str());
    } else {
        LOG("[WARN] no HyCal calibration file — energies will be wrong.");
    }

    // --- 7b. HyCal per-module time-cut table -----------------------------
    // Always populate `out.hycal_time_cuts` (uniform default when no file)
    // so per-event callers can use a single code path.  The path comes
    // from runinfo's `time_cuts.hycal_module_file` unless overridden.
    {
        std::string hc_time_path = hycal_time_cut_path_.empty()
            ? resolve(out.run_cfg.hycal_time_cut_file)
            : resolve(hycal_time_cut_path_);
        out.hycal_time_cuts = prad2::LoadHyCalTimeCuts(
            hc_time_path, out.hycal,
            out.run_cfg.hc_time_win_lo, out.run_cfg.hc_time_win_hi);
        out.hycal_time_cut_path = hc_time_path;

        std::ostringstream oss;
        oss << "[setup] HC time   : default=["
            << out.hycal_time_cuts.default_lo << ", "
            << out.hycal_time_cuts.default_hi << "] ns";
        if (out.hycal_time_cuts.n_overrides > 0) {
            oss << "  per-module=" << out.hycal_time_cuts.n_overrides
                << " (" << hc_time_path << ")";
        }
        LOG(oss.str());
    }

    // --- 8. matching (HyCal sigma + GEM sigma + target sigma) ------------
    if (recon.contains("matching")) {
        const auto &m = recon["matching"];
        if (m.contains("hycal_pos_res") && m["hycal_pos_res"].is_array()
                && m["hycal_pos_res"].size() >= 3) {
            out.hycal_pos_res[0] = m["hycal_pos_res"][0].get<float>();
            out.hycal_pos_res[1] = m["hycal_pos_res"][1].get<float>();
            out.hycal_pos_res[2] = m["hycal_pos_res"][2].get<float>();
        }
        if (m.contains("gem_pos_res") && m["gem_pos_res"].is_array()) {
            out.gem_pos_res.clear();
            for (auto &v : m["gem_pos_res"])
                out.gem_pos_res.push_back(v.get<float>());
        }
        if (m.contains("target_pos_res") && m["target_pos_res"].is_array()
                && m["target_pos_res"].size() >= 3) {
            out.target_pos_res[0] = m["target_pos_res"][0].get<float>();
            out.target_pos_res[1] = m["target_pos_res"][1].get<float>();
            out.target_pos_res[2] = m["target_pos_res"][2].get<float>();
        }
    }
    out.hycal.SetPositionResolutionParams(
        out.hycal_pos_res[0], out.hycal_pos_res[1], out.hycal_pos_res[2]);
    {
        std::ostringstream oss;
        oss << "[setup] HC sigma(E)= sqrt(("
            << std::fixed << std::setprecision(3) << out.hycal_pos_res[0]
            << "/sqrt(E_GeV))^2+(" << out.hycal_pos_res[1]
            << "/E_GeV)^2+" << out.hycal_pos_res[2] << "^2) mm";
        LOG(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "[setup] GEM sigma  : [";
        for (size_t i = 0; i < out.gem_pos_res.size(); ++i) {
            if (i) oss << ",";
            oss << std::fixed << std::setprecision(3) << out.gem_pos_res[i];
        }
        oss << "] mm";
        LOG(oss.str());
    }

    // --- 9. HyCal cluster config -----------------------------------------
    if (recon.contains("hycal") && recon["hycal"].is_object())
        apply_hycal_cluster_overrides(recon["hycal"], out.hycal_cluster_cfg);
    {
        std::ostringstream oss;
        oss << "[setup] HC cluster : min_mod_E=" << out.hycal_cluster_cfg.min_module_energy
            << "  min_ctr_E=" << out.hycal_cluster_cfg.min_center_energy
            << "  min_cl_E=" << out.hycal_cluster_cfg.min_cluster_energy
            << "  split_iter=" << out.hycal_cluster_cfg.split_iter
            << "  seed_t_win=" << out.hycal_cluster_cfg.seed_time_window << "ns"
            << (out.hycal_cluster_cfg.seed_time_window > 0.f ? " (gated)" : " (off)");
        LOG(oss.str());
    }

    // --- 10. GEM map (soft — gem stays default-constructed on failure) ---
    std::string gem_map = gem_map_path_.empty()
        ? resolve("gem_map.json")
        : resolve(gem_map_path_);
    if (gem_map.empty()) {
        LOG("[WARN] PipelineBuilder: no GEM map resolved — GEM disabled.");
    } else {
        out.gem.Init(gem_map);
        out.gem_map_path = gem_map;
        std::ostringstream oss;
        oss << "[setup] GEM map    : " << gem_map
            << "  (" << out.gem.GetNDetectors() << " detectors)";
        LOG(oss.str());
    }

    if (!out.gem_crate_remap.empty()) {
        std::ostringstream oss;
        oss << "[setup] GEM crate remap: {";
        bool first = true;
        for (const auto &kv : out.gem_crate_remap) {
            if (!first) oss << ", ";
            oss << kv.first << ": " << kv.second;
            first = false;
        }
        oss << "}";
        LOG(oss.str());
    }

    // --- 11. GEM pedestals -----------------------------------------------
    std::string ped_path = gem_pedestal_path_.empty()
        ? resolve(out.run_cfg.gem_pedestal_file)
        : resolve(gem_pedestal_path_);
    if (!ped_path.empty()) {
        out.gem.LoadPedestals(ped_path, out.gem_crate_remap);
        out.gem_pedestal_path = ped_path;
        std::ostringstream oss;
        oss << "[setup] GEM peds   : " << ped_path;
        LOG(oss.str());
    } else {
        LOG("[WARN] no GEM pedestal file — full-readout data reconstructs empty.");
    }

    // --- 12. GEM common-mode range ---------------------------------------
    std::string cm_path = gem_common_mode_path_.empty()
        ? resolve(out.run_cfg.gem_common_mode_file)
        : resolve(gem_common_mode_path_);
    if (!cm_path.empty()) {
        out.gem.LoadCommonModeRange(cm_path, out.gem_crate_remap);
        out.gem_common_mode_path = cm_path;
        std::ostringstream oss;
        oss << "[setup] GEM CM     : " << cm_path;
        LOG(oss.str());
    }

    // --- pedestal checksum (matches Python audit's [PEDSUM] line) --------
    if (log_pedestal_checksum_) {
        int n_apvs = out.gem.GetNApvs();
        double sum_noise = 0.0, sum_off = 0.0;
        long n_strips = 0;
        for (int ai = 0; ai < n_apvs; ++ai) {
            const auto &apv = out.gem.GetApvConfig(ai);
            for (int ch = 0; ch < 128; ++ch) {
                sum_noise += apv.pedestal[ch].noise;
                sum_off   += apv.pedestal[ch].offset;
                ++n_strips;
            }
        }
        std::ostringstream oss;
        oss << "[PEDSUM] n_apvs=" << n_apvs
            << " n_strips=" << n_strips
            << " sum_noise=" << std::fixed << std::setprecision(6) << sum_noise
            << " sum_offset=" << sum_off;
        LOG(oss.str());
    }

    // --- 13. GEMSYS dump (post-Init globals) -----------------------------
    {
        std::ostringstream oss;
        oss << "[GEMSYS] common_mode_thr=" << out.gem.GetCommonModeThreshold()
            << " zero_sup_thr="    << out.gem.GetZeroSupThreshold()
            << " cross_talk_thr="  << out.gem.GetCrossTalkThreshold()
            << " min_peak="        << out.gem.GetMinPeakAdc()
            << " min_sum="         << out.gem.GetMinSumAdc()
            << " rej_first="       << (int)out.gem.GetRejectFirstTimebin()
            << " rej_last="        << (int)out.gem.GetRejectLastTimebin();
        LOG(oss.str());
    }

    // --- 14. GEM per-detector cluster configs ----------------------------
    {
        gem::ClusterConfig def;
        if (recon.contains("gem") && recon["gem"].is_object()) {
            const auto &gemr = recon["gem"];
            if (gemr.contains("default") && gemr["default"].is_object())
                apply_gem_cluster_overrides(gemr["default"], def);
            std::vector<gem::ClusterConfig> per(out.gem.GetNDetectors(), def);
            for (int d = 0; d < out.gem.GetNDetectors(); ++d) {
                std::string key = std::to_string(d);
                if (gemr.contains(key) && gemr[key].is_object())
                    apply_gem_cluster_overrides(gemr[key], per[d]);
            }
            // [GEMCFG] dump (matches Python audit byte-for-byte).
            for (int d = 0; d < (int)per.size(); ++d) {
                const auto &c = per[d];
                std::ostringstream oss;
                oss << "[GEMCFG] d" << d
                    << " min_hits=" << c.min_cluster_hits
                    << " max_hits=" << c.max_cluster_hits
                    << " consec="   << c.consecutive_thres
                    << " split="    << c.split_thres
                    << " xtalk="    << c.cross_talk_width
                    << " match_mode=" << c.match_mode
                    << " asym="     << c.match_adc_asymmetry
                    << " tdiff="    << c.match_time_diff
                    << " tperiod="  << c.ts_period;
                LOG(oss.str());
            }
            out.gem.SetReconConfigs(std::move(per));
        }
    }

    return out;
}

} // namespace prad2
