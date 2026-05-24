#include "app_state.h"
#include "data_source.h"
#include "load_daq_config.h"
#include "PipelineBuilder.h"
#include "RunInfoConfig.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

using json = nlohmann::json;

namespace {

// Parse a prad_NNNNNN_LMS.dat file (produced by prad2ana_gain_monitor /
// gain_fitter).  Format (7 whitespace-separated cols, no header):
//   LMS reference rows (LMS1/2/3):
//       name alpha_peak alpha_sigma alpha_chi2 lms_peak lms_sigma lms_chi2
//       → store vals[3] (lms_peak) under the channel name.
//   HyCal module rows (W*/G*):
//       name lms_peak lms_sigma lms_chi2 gain_factor1 gain_factor2 gain_factor3
//       → store vals[0] (lms_peak) under the module name.
// Rows are classified by name prefix (not line number) so blank lines or
// reordering can't misclassify a row.  Returns true if any row parsed.
static bool parse_lms_dat(const std::string &path,
                          std::unordered_map<std::string, float> &mod_peak,
                          std::unordered_map<std::string, float> &ref_peak)
{
    std::ifstream f(path);
    if (!f.is_open()) return false;

    mod_peak.clear();
    ref_peak.clear();
    std::string line;
    while (std::getline(f, line)) {
        // strip commas so we accept both space- and comma-separated rows
        for (auto &c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        std::string name;
        if (!(iss >> name)) continue;
        std::vector<float> vals;
        float v;
        while (iss >> v) vals.push_back(v);
        if (vals.size() < 6) continue;

        if (name.rfind("LMS", 0) == 0)
            ref_peak[name] = vals[3];
        else if (fdec::HyCalSystem::name_to_id(name) >= 0)
            mod_peak[name] = vals[0];
        // unknown prefix: silently skip (keeps the parser tolerant of
        // future header lines / extra channels)
    }
    return !mod_peak.empty() || !ref_peak.empty();
}

} // namespace

//=============================================================================
// Initialization
//
// Three top-level configs are involved:
//   daq_config.json            DAQ + raw decoding (event tags, bank tags,
//                              ROC layout, sync format, file pointers).
//                              Loaded once via evc::load_daq_config(); the
//                              file pointers (hycal_map_file, gem_map_file,
//                              pedestal_file) tell us where the merged
//                              detector maps + pedestal JSON live.
//   monitor_config.json        GUI + online server (waveform/hycal_hist
//                              binning, lms_monitor, monitor_status, epics,
//                              physics display cuts, gem diagnostics, elog).
//   reconstruction_config.json runinfo pointer + cluster/hit reco knobs
//                              (hycal clustering, gem per-detector
//                              ClusterConfig with default + per-id overrides).
//=============================================================================

void AppState::init(const std::string &db_dir,
                    const std::string &daq_config_file,
                    const std::string &monitor_config_file,
                    const std::string &recon_config_file)
{
    // --- DAQ config (required, single source of truth for file pointers) ---
    std::string daq_cfg_path = daq_config_file;
    if (daq_cfg_path.empty())
        daq_cfg_path = findFile("daq_config.json", db_dir);

    if (daq_cfg_path.empty() || !evc::load_daq_config(daq_cfg_path, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_cfg_path.empty() ? " (not found)" : ": " + daq_cfg_path)
                  << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::cerr << "DAQ config: " << daq_cfg_path
              << " (adc_format=" << daq_cfg.adc_format << ")\n";

    // optional ADC1881M pedestals (PRad legacy)
    if (!daq_cfg.pedestal_file.empty()) {
        std::string ped_file = findFile(daq_cfg.pedestal_file, db_dir);
        if (!ped_file.empty() && evc::load_pedestals(ped_file, daq_cfg))
            std::cerr << "Pedestals : " << ped_file
                      << " (" << daq_cfg.pedestals.size() << " channels)\n";
    }

    // optional NNLS pile-up deconv template store.  Loaded only when the
    // analyzer config has nnls_deconv.enabled and a template_file path —
    // otherwise the store stays invalid() and every WaveAnalyzer that
    // borrows it falls back to local-maxima peak heights silently.
    if (daq_cfg.wave_cfg.nnls_deconv.enabled
        && !daq_cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        std::string tmpl_path = findFile(
            daq_cfg.wave_cfg.nnls_deconv.template_file, db_dir);
        if (tmpl_path.empty())
            tmpl_path = daq_cfg.wave_cfg.nnls_deconv.template_file;
        // LoadFromFile prints its own warning on failure.
        template_store.LoadFromFile(tmpl_path, daq_cfg.wave_cfg);
    }

    // --- resolve monitor + reconstruction config paths ---------------------
    std::string monitor_path = monitor_config_file;
    if (monitor_path.empty())
        monitor_path = findFile("monitor_config.json", db_dir);

    std::string recon_path = recon_config_file;
    if (recon_path.empty())
        recon_path = findFile("reconstruction_config.json", db_dir);

    // --- trigger definitions (needed for trigger filter parsing) -----------
    {
        std::string tbpath = findFile("trigger_bits.json", db_dir);
        std::string tbs = readFile(tbpath);
        if (!tbs.empty()) {
            auto tb = json::parse(tbs, nullptr, false);
            if (tb.is_array()) {
                trigger_bits_def = tb;
            } else if (tb.is_object()) {
                if (tb.contains("trigger_bits"))
                    trigger_bits_def = tb["trigger_bits"];
                if (tb.contains("trigger_type"))
                    trigger_type_def = tb["trigger_type"];
            }
        }
    }

    // --- load monitor config -----------------------------------------------
    json mcfg = json::object();
    if (!monitor_path.empty()) {
        std::string s = readFile(monitor_path);
        if (!s.empty()) {
            auto j = json::parse(s, nullptr, false);
            if (!j.is_discarded()) mcfg = std::move(j);
        }
    }

    // Per-peak quality bit palette (mirrors Q_PEAK_* in Fadc250Data.h).
    // Exposed via /api/config so the GUI populates the Cut-Settings dialog
    // dropdowns from a single source of truth.  Each push_back arg is a
    // 3-pair brace-init: nlohmann's auto-detector sees 3 string-keyed pairs
    // and builds an object, then push_back appends it to the array.
    peak_quality_bits_def = json::array();
    peak_quality_bits_def.push_back({{"bit", 0}, {"name", "PILED"},       {"label", "Pile-up"}});
    peak_quality_bits_def.push_back({{"bit", 1}, {"name", "DECONVOLVED"}, {"label", "Deconvolved"}});

    // waveform binning + trigger filter (monitor side)
    if (mcfg.contains("waveform")) {
        auto &w = mcfg["waveform"];
        waveform_trigger.parse(w, trigger_bits_def);
        if (w.contains("filter") && w["filter"].is_object())
            peak_filter.parse(w["filter"], peak_quality_bits_def);
        // Snapshot for the Cut-Settings "Reset" button.  `enable` defaults
        // to true so Reset both restores the JSON ranges and re-arms the
        // apply toggle (matches startup state).
        peak_filter_default = peak_filter;
        if (w.contains("integral_hist")) {
            auto &ih = w["integral_hist"];
            if (ih.contains("min"))  hist_cfg.bin_min  = ih["min"];
            if (ih.contains("max"))  hist_cfg.bin_max  = ih["max"];
            if (ih.contains("step")) hist_cfg.bin_step = ih["step"];
        }
        if (w.contains("time_hist")) {
            auto &th = w["time_hist"];
            if (th.contains("min"))  hist_cfg.pos_min  = th["min"];
            if (th.contains("max"))  hist_cfg.pos_max  = th["max"];
            if (th.contains("step")) hist_cfg.pos_step = th["step"];
        }
        if (w.contains("height_hist")) {
            auto &hh = w["height_hist"];
            if (hh.contains("min"))  hist_cfg.height_min  = hh["min"];
            if (hh.contains("max"))  hist_cfg.height_max  = hh["max"];
            if (hh.contains("step")) hist_cfg.height_step = hh["step"];
        }
        // Note: peak detection thresholds (peak_nsigma, min_peak_height,
        // min_peak_ratio) live in daq_config.json `fadc250_waveform.analyzer`
        // and are loaded into wave_cfg by the prad2dec config loader.
    }
    // ref_lines: assemble the flat key→[lines] map the frontend expects
    // (`refShapes(key)` consumes it).  Two sources, in order:
    //   1. Nested form — any object like `<section>.<subkey>.ref_lines` is
    //      stored under `<subkey>`.  This lets each `*_hist` block keep its
    //      reference lines next to its own min/max/step.
    //   2. Top-level `ref_lines` — for plots that don't have a natural
    //      sub-block (raw waveform, lms, physics scatter plots).  Top-level
    //      entries override nested ones if both are defined.
    ref_lines = json::object();
    for (auto it = mcfg.begin(); it != mcfg.end(); ++it) {
        if (it.key() == "ref_lines" || !it.value().is_object()) continue;
        for (auto sit = it.value().begin(); sit != it.value().end(); ++sit) {
            if (!sit.value().is_object()) continue;
            if (sit.value().contains("ref_lines")
                && sit.value()["ref_lines"].is_array())
                ref_lines[sit.key()] = sit.value()["ref_lines"];
        }
    }
    if (mcfg.contains("ref_lines") && mcfg["ref_lines"].is_object())
        for (auto it = mcfg["ref_lines"].begin(); it != mcfg["ref_lines"].end(); ++it)
            ref_lines[it.key()] = it.value();

    hist_nbins = std::max(1, (int)std::ceil(
        (hist_cfg.bin_max - hist_cfg.bin_min) / hist_cfg.bin_step));
    pos_nbins = std::max(1, (int)std::ceil(
        (hist_cfg.pos_max - hist_cfg.pos_min) / hist_cfg.pos_step));
    height_nbins = std::max(1, (int)std::ceil(
        (hist_cfg.height_max - hist_cfg.height_min) / hist_cfg.height_step));
    {
        std::cerr << "Waveform  : peak_nsigma=" << daq_cfg.wave_cfg.peak_nsigma
                  << " min_peak_height=" << daq_cfg.wave_cfg.min_peak_height
                  << " filter[" << (peak_filter.enable ? "on" : "off") << "]="
                  << peak_filter.toJson(peak_quality_bits_def).dump()
                  << " " << waveform_trigger << "\n";
    }

    // --- detector pipeline (HyCal + GEM, runinfo, recon config) -----------
    // PipelineBuilder consolidates HyCal/GEM Init + LoadCalibration +
    // LoadPedestals + LoadCommonModeRange + per-detector ClusterConfig +
    // DetectorTransform construction + matching parameters into one call.
    // We hand it our already-loaded daq_cfg (the legacy ADC pedestal +
    // NNLS template store steps above need it earlier than the builder).
    // findFile is wired in as the path resolver so the multi-dir search
    // semantics are preserved for non-default file locations.
    {
        // Capture daq_cfg-sourced overrides before the std::move below.
        std::string hycal_map_override = daq_cfg.hycal_map_file;
        std::string gem_map_override   = daq_cfg.gem_map_file;

        prad2::Pipeline pipeline = prad2::PipelineBuilder()
            .set_database_dir(db_dir)
            .set_loaded_daq_config(std::move(daq_cfg))
            .set_daq_config(daq_cfg_path)            // logging only
            .set_recon_config(recon_path)
            .set_hycal_map(std::move(hycal_map_override))
            .set_gem_map(std::move(gem_map_override))
            .set_path_resolver(
                [&](const std::string &p) { return findFile(p, db_dir); })
            .set_log_stream(&std::cerr)
            .build();

        // Move pipeline contents into AppState members.
        daq_cfg              = std::move(pipeline.daq_cfg);
        hycal                = std::move(pipeline.hycal);
        gem_sys              = std::move(pipeline.gem);
        cluster_cfg          = pipeline.hycal_cluster_cfg;
        hycal_transform      = pipeline.hycal_transform;
        gem_transforms.assign(pipeline.gem_transforms.begin(),
                              pipeline.gem_transforms.end());
        gem_pos_res          = std::move(pipeline.gem_pos_res);

        beam_energy_runinfo  = pipeline.run_cfg.Ebeam;
        beam_energy.store(pipeline.run_cfg.Ebeam);
        target_x             = pipeline.run_cfg.target_x;
        target_y             = pipeline.run_cfg.target_y;
        target_z             = pipeline.run_cfg.target_z;
        adc_to_mev           = pipeline.run_cfg.default_adc2mev;

        gem_eff_target_sigma_x = pipeline.target_pos_res[0];
        gem_eff_target_sigma_y = pipeline.target_pos_res[1];
        gem_eff_target_sigma_z = pipeline.target_pos_res[2];

        gem_enabled = (gem_sys.GetNDetectors() > 0);
        if (gem_enabled) {
            gem_occupancy.assign(gem_sys.GetNDetectors(), Histogram2D{});
            for (auto &h : gem_occupancy) h.init(GEM_OCC_NX, GEM_OCC_NY);
        }
        std::cerr << "HyCal     : " << hycal.module_count() << " modules\n";
        if (gem_enabled)
            std::cerr << "GEM       : " << gem_sys.GetNDetectors() << " detectors\n";
    }

    // --- crate_roc map (directly from daq_cfg.roc_tags) -------------------
    crate_roc_json = json::object();
    for (const auto &re : daq_cfg.roc_tags) {
        // only data ROCs (type "roc"/"gem"); ti_slaves share crate numbers
        // but have different tags and must not overwrite the data ROC entry.
        if (!re.type.empty() && re.type != "roc" && re.type != "gem") continue;
        if (re.crate < 0) continue;
        crate_roc_json[std::to_string(re.crate)] = re.tag;
    }
    if (crate_roc_json.empty())
        crate_roc_json = {{"0",0x80},{"1",0x82},{"2",0x84},{"3",0x86},{"4",0x88},{"5",0x8a},{"6",0x8c}};

    roc_to_crate.clear();
    for (auto &[k, v] : crate_roc_json.items())
        roc_to_crate[v.get<int>()] = std::stoi(k);

    // --- monitor config: remaining sections -------------------------------
    if (mcfg.contains("online")) {
        auto &on = mcfg["online"];
        if (on.contains("refresh_ms")) {
            auto &r = on["refresh_ms"];
            if (r.contains("event"))     refresh_event_ms = r["event"];
            if (r.contains("ring"))      refresh_ring_ms  = r["ring"];
            if (r.contains("histogram")) refresh_hist_ms  = r["histogram"];
            if (r.contains("lms"))       refresh_lms_ms   = r["lms"];
        }
    }

    // hycal_hist: trigger filter + display-histogram binning for the cluster
    // monitor.  Cluster-reco knobs (min_*_energy, split_iter, ...) come from
    // reconstruction_config.json:hycal further below.
    if (mcfg.contains("hycal_hist")) {
        auto &hh = mcfg["hycal_hist"];
        cluster_trigger.parse(hh, trigger_bits_def);
        if (hh.contains("energy_hist")) {
            auto &eh = hh["energy_hist"];
            if (eh.contains("min"))  cl_hist_min  = eh["min"];
            if (eh.contains("max"))  cl_hist_max  = eh["max"];
            if (eh.contains("step")) cl_hist_step = eh["step"];
        }
        if (hh.contains("nclusters_hist")) {
            auto &nh = hh["nclusters_hist"];
            if (nh.contains("min"))  nclusters_hist_min  = nh["min"];
            if (nh.contains("max"))  nclusters_hist_max  = nh["max"];
            if (nh.contains("step")) nclusters_hist_step = nh["step"];
        }
        if (hh.contains("nblocks_hist")) {
            auto &bh = hh["nblocks_hist"];
            if (bh.contains("min"))  nblocks_hist_min  = bh["min"];
            if (bh.contains("max"))  nblocks_hist_max  = bh["max"];
            if (bh.contains("step")) nblocks_hist_step = bh["step"];
        }
        if (hh.contains("raw_energy_hist")) {
            auto &rh = hh["raw_energy_hist"];
            if (rh.contains("min"))  raw_energy_hist_min  = rh["min"];
            if (rh.contains("max"))  raw_energy_hist_max  = rh["max"];
            if (rh.contains("step")) raw_energy_hist_step = rh["step"];
        }
    }

    if (mcfg.contains("lms_monitor")) {
        auto &lm = mcfg["lms_monitor"];
        lms_trigger.parse(lm, trigger_bits_def);
        if (lm.contains("time_cut")) {
            auto &tc = lm["time_cut"];
            if (tc.contains("min")) lms_time_min = tc["min"];
            if (tc.contains("max")) lms_time_max = tc["max"];
        }
        if (lm.contains("warn_threshold")) lms_warn_thresh    = lm["warn_threshold"];
        if (lm.contains("warn_min_mean"))  lms_warn_min_mean  = lm["warn_min_mean"];
        if (lm.contains("max_history"))    lms_max_history    = lm["max_history"];
        if (lm.contains("reference_channels")) {
            for (auto &name : lm["reference_channels"]) {
                std::string n = name.get<std::string>();
                const auto *mod = hycal.module_by_name(n);
                lms_ref_channels.push_back({n, mod ? mod->index : -1});
            }
        }
        if (lm.contains("alpha"))
            alpha_trigger.parse(lm["alpha"], trigger_bits_def);

        // Gain-drift baseline (optional). When configured, apiLmsSummary
        // computes drift = current_LMS / (baseline_LMS * lamp_scale) per
        // module and flags drift outside [drift_low, drift_high] as a
        // top-priority error (above the existing rms/floor warn).
        if (lm.contains("drift_baseline_file")) {
            std::string p = lm["drift_baseline_file"].get<std::string>();
            if (!p.empty() && p[0] != '/') p = db_dir + "/" + p;
            lms_drift_baseline_file = p;
        }
        if (lm.contains("drift_ref_channel"))
            lms_drift_ref_channel = lm["drift_ref_channel"].get<std::string>();
        // Common-knob fallback: 'drift_warn_low/high' applies to both W and G
        // unless overridden by the per-type keys below.
        if (lm.contains("drift_warn_low")) {
            float v = lm["drift_warn_low"];
            lms_drift_low_w = v; lms_drift_low_g = v;
        }
        if (lm.contains("drift_warn_high")) {
            float v = lm["drift_warn_high"];
            lms_drift_high_w = v; lms_drift_high_g = v;
        }
        // Per-module-type overrides (recommended).  Use these to handle the
        // different healthy widths of PbWO4 vs PbGlass.
        if (lm.contains("drift_warn_low_w"))  lms_drift_low_w  = lm["drift_warn_low_w"];
        if (lm.contains("drift_warn_high_w")) lms_drift_high_w = lm["drift_warn_high_w"];
        if (lm.contains("drift_warn_low_g"))  lms_drift_low_g  = lm["drift_warn_low_g"];
        if (lm.contains("drift_warn_high_g")) lms_drift_high_g = lm["drift_warn_high_g"];

        // Suppress drift WARNINGS (not the drift value itself) for these
        // module types — names follow hycal_map.json's "t" field.  Useful when
        // a known systematic between baseline calibrations would otherwise
        // flood the alarm list.  Unknown type names are ignored with a warning.
        if (lm.contains("drift_warn_suppress_types")
            && lm["drift_warn_suppress_types"].is_array())
        {
            for (auto &el : lm["drift_warn_suppress_types"]) {
                if (!el.is_string()) continue;
                std::string tname = el.get<std::string>();
                fdec::ModuleType t = fdec::HyCalSystem::parse_type(tname);
                if (t == fdec::ModuleType::Unknown) {
                    std::cerr << "[WARN] LMS drift_warn_suppress_types: unknown type '"
                              << tname << "' — accepted: PbGlass, PbWO4, LMS, Veto\n";
                    continue;
                }
                lms_drift_suppress_types.insert(static_cast<int>(t));
                lms_drift_suppress_type_names.push_back(tname);
            }
        }

        bool drift_loaded = false;
        if (!lms_drift_baseline_file.empty()) {
            drift_loaded = parse_lms_dat(lms_drift_baseline_file,
                                         lms_baseline_peak,
                                         lms_baseline_ref_peak);
            if (!drift_loaded)
                std::cerr << "[WARN] LMS drift baseline could not be loaded: "
                          << lms_drift_baseline_file << "\n";
            else if (lms_baseline_peak.empty() || lms_baseline_ref_peak.empty())
                std::cerr << "[WARN] LMS drift baseline missing "
                          << (lms_baseline_peak.empty() ? "module" : "reference")
                          << " rows — drift detection will be disabled. File: "
                          << lms_drift_baseline_file << "\n";
        }

        std::cerr << "LMS       : " << lms_trigger
                  << " time_cut=[" << lms_time_min << "," << lms_time_max << "]"
                  << " warn=" << lms_warn_thresh
                  << " refs=" << lms_ref_channels.size()
                  << " alpha=" << alpha_trigger;
        if (drift_loaded) {
            std::cerr << " drift_baseline=" << lms_drift_baseline_file
                      << " (" << lms_baseline_peak.size() << " mods, "
                      << lms_baseline_ref_peak.size() << " refs)"
                      << " drift_W=[" << lms_drift_low_w << "," << lms_drift_high_w << "]"
                      << " drift_G=[" << lms_drift_low_g << "," << lms_drift_high_g << "]"
                      << " ref_ch='" << lms_drift_ref_channel << "'";
            if (!lms_drift_suppress_type_names.empty()) {
                std::cerr << " suppress=[";
                for (size_t i = 0; i < lms_drift_suppress_type_names.size(); ++i) {
                    if (i) std::cerr << ",";
                    std::cerr << lms_drift_suppress_type_names[i];
                }
                std::cerr << "]";
            }
        }
        std::cerr << "\n";
    }

    // monitor_status: nested {livetime, beam: {energy, current}} — header
    // status panel for online mode (livetime + beam shell-poll metrics).
    if (mcfg.contains("monitor_status")) {
        auto parse_shell_metric = [](const nlohmann::json &j, ShellMetric &m) {
            if (j.contains("command"))  m.command  = j["command"].get<std::string>();
            if (j.contains("unit"))     m.unit     = j["unit"].get<std::string>();
            if (j.contains("poll_sec")) m.poll_sec = std::max(1, (int)j["poll_sec"]);
            if (j.contains("trip_warn_below")) {
                m.has_trip_warn   = true;
                m.trip_warn_below = j["trip_warn_below"];
            }
        };
        auto &ms = mcfg["monitor_status"];
        if (ms.contains("livetime")) {
            auto &lt = ms["livetime"];
            if (lt.contains("command"))  livetime_cmd      = lt["command"].get<std::string>();
            if (lt.contains("unit"))     livetime_unit     = lt["unit"].get<std::string>();
            if (lt.contains("poll_sec")) livetime_poll_sec = std::max(1, (int)lt["poll_sec"]);
            if (lt.contains("healthy"))  livetime_healthy  = lt["healthy"];
            if (lt.contains("warning"))  livetime_warning  = lt["warning"];
        }
        if (ms.contains("beam")) {
            auto &beam = ms["beam"];
            if (beam.contains("energy"))  parse_shell_metric(beam["energy"],  beam_energy_status);
            if (beam.contains("current")) parse_shell_metric(beam["current"], beam_current_status);
        }
    }

    const auto &ds = daq_cfg.dsc_scaler;
    using DSrc = evc::DaqConfig::DscScaler::Source;
    const char *src_name = (ds.source == DSrc::Ref) ? "ref"
                         : (ds.source == DSrc::Trg) ? "trg" : "tdc";
    std::cerr << "Livetime  : "
              << (livetime_cmd.empty() ? "disabled"
                                       : ("'" + livetime_cmd + "' every "
                                          + std::to_string(livetime_poll_sec) + "s"))
              << " healthy>=" << livetime_healthy
              << " warn>=" << livetime_warning;
    if (ds.enabled()) {
        std::cerr << " | DSC2 bank=0x" << std::hex << ds.bank_tag << std::dec
                  << " slot=" << ds.slot << " src=" << src_name;
        if (ds.source != DSrc::Ref)
            std::cerr << " ch=" << ds.channel;
    } else {
        std::cerr << " | DSC2 disabled";
    }
    std::cerr << "\n";

    auto log_metric = [](const char *label, const ShellMetric &m) {
        std::cerr << label
                  << (m.command.empty() ? "disabled"
                                        : ("'" + m.command + "' every "
                                           + std::to_string(m.poll_sec) + "s"
                                           + " unit=" + m.unit));
        if (m.has_trip_warn)
            std::cerr << " trip_warn<" << m.trip_warn_below;
        std::cerr << "\n";
    };
    log_metric("BeamE     : ", beam_energy_status);
    log_metric("BeamI     : ", beam_current_status);

    if (mcfg.contains("color_ranges")) {
        for (auto &[key, val] : mcfg["color_ranges"].items()) {
            if (val.is_array() && val.size() == 2)
                color_range_defaults[key] = {val[0].get<float>(), val[1].get<float>()};
        }
        std::cerr << "Color ranges: " << color_range_defaults.size() << " entries\n";
    }

    // Auto-report config (the elog block lives nested under auto_report
    // since elog writes are now driven exclusively by the auto-report
    // pipeline — no manual Post-to-Elog dialog).
    if (mcfg.contains("auto_report")) {
        auto &ar = mcfg["auto_report"];
        if (ar.contains("enabled"))
            auto_report_enabled = ar["enabled"].get<bool>();
        if (ar.contains("post_to_elog"))
            auto_report_post_to_elog = ar["post_to_elog"].get<bool>();
        if (ar.contains("local_save_dir"))
            auto_report_local_save_dir = ar["local_save_dir"].get<std::string>();
        if (ar.contains("min_interval_ms"))
            auto_report_min_interval_ms = ar["min_interval_ms"].get<int>();
        if (ar.contains("schedule_minutes"))
            auto_report_schedule_minutes = ar["schedule_minutes"].get<int>();
        if (ar.contains("min_events_for_schedule"))
            auto_report_min_events_for_schedule = ar["min_events_for_schedule"].get<int>();
        if (ar.contains("schedule_max_wait_min"))
            auto_report_schedule_max_wait_min = ar["schedule_max_wait_min"].get<int>();
        if (ar.contains("partial_threshold_events"))
            auto_report_partial_threshold_events = ar["partial_threshold_events"].get<int>();
        if (ar.contains("elog")) {
            auto &el = ar["elog"];
            if (el.contains("url"))     elog_url     = el["url"];
            if (el.contains("logbook")) elog_logbook = el["logbook"];
            if (el.contains("author"))  elog_author  = el["author"];
            if (el.contains("tags"))
                for (auto &t : el["tags"]) elog_tags.push_back(t);
            if (el.contains("cert")) elog_cert = el["cert"];
            if (el.contains("key"))  elog_key  = el["key"];
        }
    }
    std::cerr << "AutoReport: enabled=" << (auto_report_enabled ? "ON" : "OFF")
              << " post_to_elog=" << (auto_report_post_to_elog ? "yes" : "no")
              << " min_interval_ms=" << auto_report_min_interval_ms
              << " schedule_minutes=" << auto_report_schedule_minutes
              << " min_events_for_schedule=" << auto_report_min_events_for_schedule
              << " schedule_max_wait_min=" << auto_report_schedule_max_wait_min
              << " partial_threshold_events=" << auto_report_partial_threshold_events
              << (auto_report_local_save_dir.empty() ? std::string()
                  : " local_save=" + auto_report_local_save_dir)
              << "\n";
    std::cerr << "Elog      : " << elog_url
              << " logbook=" << elog_logbook
              << (elog_cert.empty() ? "" : " cert=" + elog_cert)
              << "\n";

    if (mcfg.contains("physics")) {
        auto &ph = mcfg["physics"];
        physics_trigger.parse(ph, trigger_bits_def);
        if (ph.contains("beam_energy")) {
            auto &be = ph["beam_energy"];
            if (be.contains("epics_channel")) beam_energy_epics_channel = be["epics_channel"];
            if (be.contains("min_valid"))     beam_energy_min_valid     = be["min_valid"];
        }
        if (ph.contains("energy_angle_hist")) {
            auto &ea = ph["energy_angle_hist"];
            if (ea.contains("angle_min"))   ea_angle_min   = ea["angle_min"];
            if (ea.contains("angle_max"))   ea_angle_max   = ea["angle_max"];
            if (ea.contains("angle_step"))  ea_angle_step  = ea["angle_step"];
            if (ea.contains("energy_min"))  ea_energy_min  = ea["energy_min"];
            if (ea.contains("energy_max"))  ea_energy_max  = ea["energy_max"];
            if (ea.contains("energy_step")) ea_energy_step = ea["energy_step"];
        }
        if (ph.contains("moller")) {
            auto &ml = ph["moller"];
            if (ml.contains("energy_tolerance")) moller_energy_tol = ml["energy_tolerance"];
            if (ml.contains("angle_min"))        moller_angle_min  = ml["angle_min"];
            if (ml.contains("angle_max"))        moller_angle_max  = ml["angle_max"];
            if (ml.contains("xy_hist")) {
                auto &xy = ml["xy_hist"];
                if (xy.contains("x_min"))  moller_xy_x_min  = xy["x_min"];
                if (xy.contains("x_max"))  moller_xy_x_max  = xy["x_max"];
                if (xy.contains("x_step")) moller_xy_x_step = xy["x_step"];
                if (xy.contains("y_min"))  moller_xy_y_min  = xy["y_min"];
                if (xy.contains("y_max"))  moller_xy_y_max  = xy["y_max"];
                if (xy.contains("y_step")) moller_xy_y_step = xy["y_step"];
            }
        }
        if (ph.contains("hycal_cluster_hit")) {
            auto &hc = ph["hycal_cluster_hit"];
            if (hc.contains("n_clusters"))      hxy_n_clusters      = hc["n_clusters"];
            if (hc.contains("energy_frac_min")) hxy_energy_frac_min = hc["energy_frac_min"];
            if (hc.contains("nblocks_min"))     hxy_nblocks_min     = hc["nblocks_min"];
            if (hc.contains("nblocks_max"))     hxy_nblocks_max     = hc["nblocks_max"];
            if (hc.contains("xy_hist")) {
                auto &xy = hc["xy_hist"];
                if (xy.contains("x_min"))  hxy_x_min  = xy["x_min"];
                if (xy.contains("x_max"))  hxy_x_max  = xy["x_max"];
                if (xy.contains("x_step")) hxy_x_step = xy["x_step"];
                if (xy.contains("y_min"))  hxy_y_min  = xy["y_min"];
                if (xy.contains("y_max"))  hxy_y_max  = xy["y_max"];
                if (xy.contains("y_step")) hxy_y_step = xy["y_step"];
            }
        }
        std::cerr << "Physics   : " << physics_trigger
                  << " Moller: tol=" << moller_energy_tol
                  << " angle=[" << moller_angle_min << "," << moller_angle_max << "]"
                  << " HyCalXY: Ncl=" << hxy_n_clusters
                  << " E>=" << hxy_energy_frac_min << "*Eb"
                  << " blocks=[" << hxy_nblocks_min << "," << hxy_nblocks_max << "]"
                  << " beam_src='" << beam_energy_epics_channel << "'\n";
    }

    if (mcfg.contains("epics")) {
        auto &ep = mcfg["epics"];
        if (ep.contains("max_history"))     epics_max_history  = ep["max_history"];
        if (ep.contains("warn_threshold"))  epics_warn_thresh  = ep["warn_threshold"];
        if (ep.contains("alert_threshold")) epics_alert_thresh = ep["alert_threshold"];
        if (ep.contains("min_avg_points"))  epics_min_avg_pts  = ep["min_avg_points"];
        if (ep.contains("mean_window"))     epics_mean_window  = ep["mean_window"];
        if (ep.contains("slots")) {
            for (auto &slot : ep["slots"]) {
                std::vector<std::string> names;
                for (auto &ch : slot) names.push_back(ch);
                epics_default_slots.push_back(std::move(names));
            }
        }
        std::cerr << "EPICS     : max_history=" << epics_max_history
                  << " slots=" << epics_default_slots.size() << "\n";
    }

    // GEM diagnostics (HyCal-anchored matching residuals + tracking efficiency).
    // Not to be confused with reconstruction_config.json:gem (per-detector
    // ClusterConfig); this section is monitor-side only.
    if (mcfg.contains("gem")) {
        auto &gemcfg = mcfg["gem"];
        if (gemcfg.contains("hycal_match")) {
            auto &gm = gemcfg["hycal_match"];
            if (gm.contains("require_ep_candidate")) gem_match_require_ep = gm["require_ep_candidate"];
            if (gm.contains("match_nsigma"))         gem_match_nsigma     = gm["match_nsigma"];
            if (gm.contains("residual_hist")) {
                auto &rh = gm["residual_hist"];
                if (rh.contains("min"))  gem_resid_min  = rh["min"];
                if (rh.contains("max"))  gem_resid_max  = rh["max"];
                if (rh.contains("step")) gem_resid_step = rh["step"];
            }
        }
        if (gemcfg.contains("efficiency")) {
            auto &ge = gemcfg["efficiency"];
            if (ge.contains("min_cluster_energy"))    gem_eff_min_cluster_energy = ge["min_cluster_energy"];
            if (ge.contains("match_nsigma"))          gem_eff_match_nsigma       = ge["match_nsigma"];
            if (ge.contains("max_chi2_per_dof"))      gem_eff_max_chi2           = ge["max_chi2_per_dof"];
            if (ge.contains("max_hits_per_detector")) gem_eff_max_hits_per_det   = ge["max_hits_per_detector"];
            if (ge.contains("min_denom_for_eff"))     gem_eff_min_denom          = ge["min_denom_for_eff"];
            if (ge.contains("healthy"))               gem_eff_healthy            = ge["healthy"];
            if (ge.contains("warning"))               gem_eff_warning            = ge["warning"];
            if (ge.contains("loo_mode")) {
                std::string m = ge["loo_mode"].get<std::string>();
                if      (m == "loo")              gem_eff_loo_mode = GemEffLooMode::Loo;
                else if (m == "loo-target-in")    gem_eff_loo_mode = GemEffLooMode::LooTargetIn;
                else if (m == "loo-target-seed")  gem_eff_loo_mode = GemEffLooMode::TargetSeed;
                else std::cerr << "[WARN] gem.efficiency.loo_mode='" << m
                               << "' unknown; falling back to loo-target-seed\n";
            }
            if (ge.contains("z_target_hist")) {
                auto &zh = ge["z_target_hist"];
                if (zh.contains("min"))  gem_eff_z_target_min  = zh["min"];
                if (zh.contains("max"))  gem_eff_z_target_max  = zh["max"];
                if (zh.contains("step")) gem_eff_z_target_step = zh["step"];
            }
            if (ge.contains("local_grid")) {
                auto &lg = ge["local_grid"];
                if (lg.contains("nx")) gem_eff_grid_nx = lg["nx"];
                if (lg.contains("ny")) gem_eff_grid_ny = lg["ny"];
                if (gem_eff_grid_nx < 1) gem_eff_grid_nx = 1;
                if (gem_eff_grid_ny < 1) gem_eff_grid_ny = 1;
            }
        }
        std::cerr << "GEM cfg   : match=" << gem_match_nsigma << "σ"
                  << "  efficiency=" << gem_eff_match_nsigma << "σ"
                  << "  chi2/dof<=" << gem_eff_max_chi2 << "\n";
    }

    std::cerr << "Monitor   : " << (monitor_path.empty() ? "(none)" : monitor_path) << "\n";
    std::cerr << "Clustering: min_mod=" << cluster_cfg.min_module_energy
              << " min_center=" << cluster_cfg.min_center_energy
              << " min_cluster=" << cluster_cfg.min_cluster_energy
              << " " << cluster_trigger
              << " hist=[" << cl_hist_min << "," << cl_hist_max
              << "]/" << cl_hist_step << "\n";
    std::cerr << "Reco      : " << (recon_path.empty() ? "(none)" : recon_path)
              << " (adc_to_mev=" << adc_to_mev << ")\n";

    // --- init derived histograms ------------------------------------------
    int cl_nbins = std::max(1, (int)std::ceil((cl_hist_max - cl_hist_min) / cl_hist_step));
    cluster_energy_hist.init(cl_nbins);
    int nb_nclusters = std::max(1, (int)std::ceil(
        (nclusters_hist_max - nclusters_hist_min) / nclusters_hist_step));
    nclusters_hist.init(nb_nclusters);
    int nb_blocks = std::max(1, (nblocks_hist_max - nblocks_hist_min) / nblocks_hist_step);
    nblocks_hist.init(nb_blocks);
    int re_nbins = std::max(1, (int)std::ceil(
        (raw_energy_hist_max - raw_energy_hist_min) / raw_energy_hist_step));
    raw_energy_hist.init(re_nbins);
    cluster_energy_hist_by_ncl.assign(nb_nclusters, Histogram{});
    nblocks_hist_by_ncl.assign(nb_nclusters, Histogram{});
    for (auto &h : cluster_energy_hist_by_ncl) h.init(cl_nbins);
    for (auto &h : nblocks_hist_by_ncl)        h.init(nb_blocks);
    int ea_nx = std::max(1, (int)std::ceil((ea_angle_max - ea_angle_min) / ea_angle_step));
    int ea_ny = std::max(1, (int)std::ceil((ea_energy_max - ea_energy_min) / ea_energy_step));
    energy_angle_hist.init(ea_nx, ea_ny);
    int ml_nx = std::max(1, (int)std::ceil((moller_xy_x_max - moller_xy_x_min) / moller_xy_x_step));
    int ml_ny = std::max(1, (int)std::ceil((moller_xy_y_max - moller_xy_y_min) / moller_xy_y_step));
    moller_xy_hist.init(ml_nx, ml_ny);
    int hxy_nx = std::max(1, (int)std::ceil((hxy_x_max - hxy_x_min) / hxy_x_step));
    int hxy_ny = std::max(1, (int)std::ceil((hxy_y_max - hxy_y_min) / hxy_y_step));
    hycal_xy_hist.init(hxy_nx, hxy_ny);
    {
        int n_gem = gem_enabled ? (int)gem_transforms.size() : 0;
        int resid_nbins = std::max(1, (int)std::ceil((gem_resid_max - gem_resid_min) / gem_resid_step));
        gem_dx_hist.assign(n_gem, Histogram{});
        gem_dy_hist.assign(n_gem, Histogram{});
        gem_match_hits.assign(n_gem, 0);
        for (auto &h : gem_dx_hist) h.init(resid_nbins);
        for (auto &h : gem_dy_hist) h.init(resid_nbins);
    }
    initGemEfficiency();
    // hycal_transform is already prepared by setTransform() above (or, if
    // runinfo wasn't loaded, will lazy-prepare on first toLab/rotate/matrix
    // call) — no eager prepare needed here.
}
