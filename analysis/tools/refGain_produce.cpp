// refGain_produce.cpp — EVIO -> LMS/alpha gain factors (single-pass, global fit)
//
// Reads EVIO files directly, classifies events by channel count,
// accumulates histograms over ALL events (one global fit), then writes:
//   * <db>/gain_factor/ref_gain/prad_XXXXXX_LMS.dat  (fit results)
//   * <db>/gain_factor/ref_gain/prad_XXXXXX_LMS_hists.root  (histograms)
//
// Usage:
//   refGain_produce <evio_file_or_dir> [more files/dirs...]
//              [-o output.dat]    default: <db>/gain_factor/ref_gain/prad_XXXXXX_LMS.dat
//              [-r hists.root]    default: <db>/gain_factor/ref_gain/prad_XXXXXX_LMS_hists.root
//              [-c daq_config.json]
//              [-d hycal_map.json]
//              [-f max_files]
//              [-n max_events]

#include "Replay.h"
#include "EventData.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "load_daq_config.h"
#include "DaqConfig.h"
#include "EvChannel.h"
#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"
#include "Fadc250Data.h"

#include <TClass.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
// gain_factor.h uses TH1F -- include TH1F.h first
#include "gain_factor.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;

// ── constants (mirror replay_gainCorr.cpp) ───────────────────────────────────
static constexpr int   N_W         = 1156;
static constexpr int   N_LMS       = 3;
static constexpr int   LMS_ID_BASE = 3101;   // LMS1=3101, LMS2=3102, LMS3=3103
static constexpr int   W_ID_BASE   = 1000;   // PbWO4: module_id = W-id + 1000
static constexpr int   HIST_BINS   = 600;
static constexpr float HIST_MIN    = 0.f;
static constexpr float HIST_MAX    = 15000.f;

// Create a histogram detached from any ROOT directory (heap-only).
static TH1F *makeH(const char *name)
{
    auto *h = new TH1F(name, name, HIST_BINS, HIST_MIN, HIST_MAX);
    h->SetDirectory(nullptr);  // detach from gROOT/any TFile
    return h;
}

static std::vector<std::string> collectEvioFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (std::filesystem::is_directory(path)) {
        for (auto &e : std::filesystem::directory_iterator(path)) {
            if (e.is_regular_file() &&
                e.path().filename().string().find(".evio") != std::string::npos)
                files.push_back(e.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

int main(int argc, char *argv[])
{
    TClass::GetClass("TFile");

    std::string output_dat, output_root, daq_config, daq_map;
    int max_events = -1;
    int max_files  = -1;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    daq_config = db_dir + "/daq_config.json";

    int opt;
    while ((opt = getopt(argc, argv, "o:r:c:d:f:n:")) != -1) {
        switch (opt) {
            case 'o': output_dat  = optarg; break;
            case 'r': output_root = optarg; break;
            case 'c': daq_config  = optarg; break;
            case 'd': daq_map     = optarg; break;
            case 'f': max_files   = std::atoi(optarg); break;
            case 'n': max_events  = std::atoi(optarg); break;
        }
    }

    std::vector<std::string> evio_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectEvioFiles(argv[i]);
        evio_files.insert(evio_files.end(), f.begin(), f.end());
    }

    if (evio_files.empty()) {
        std::cerr <<
            "Usage: refGain_produce <evio_file_or_dir> [...]\n"
            "       [-o output.dat]    default: <db>/gain_factor/ref_gain/prad_XXXXXX_LMS.dat\n"
            "       [-r hists.root]    default: <db>/gain_factor/ref_gain/prad_XXXXXX_LMS_hists.root\n"
            "       [-c daq_config.json] [-d hycal_map.json] [-f max_files] [-n max_events]\n";
        return 1;
    }

    if (daq_map.empty()) daq_map = db_dir + "/hycal_map.json";
    if (max_files > 0 && (int)evio_files.size() > max_files)
        evio_files.resize(max_files);

    int run = get_run_int(evio_files[0]);
    std::string ref_gain_dir = db_dir + "/gain_factor/ref_gain";
    {
        std::error_code ec;
        std::filesystem::create_directories(ref_gain_dir, ec);
        if (ec)
            std::cerr << "Warning: cannot create " << ref_gain_dir << ": " << ec.message() << "\n";
    }
    if (output_dat.empty())
        output_dat  = ref_gain_dir + "/" + Form("prad_%06d_LMS.dat",       run);
    if (output_root.empty())
        output_root = ref_gain_dir + "/" + Form("prad_%06d_LMS_hists.root", run);

    std::cout << "  Files   : " << evio_files.size() << "\n"
              << "  DAQ cfg : " << daq_config << "\n"
              << "  HyCal   : " << daq_map << "\n"
              << "  Out dat : " << output_dat << "\n"
              << "  Out root: " << output_root << "\n";

    // ── Load DAQ config ───────────────────────────────────────────────────────
    evc::DaqConfig daq_cfg;
    evc::load_daq_config(daq_config, daq_cfg);

    // ── Load module map ───────────────────────────────────────────────────────
    analysis::Replay replay;
    replay.LoadDaqConfig(daq_config);
    replay.LoadHyCalMap(daq_map);

    // ── ROC tag -> crate mapping ──────────────────────────────────────────────
    std::unordered_map<int, int> roc_to_crate;
    {
        std::ifstream dcf(daq_config);
        if (dcf.is_open()) {
            auto dcj = nlohmann::json::parse(dcf, nullptr, false, true);
            if (dcj.contains("roc_tags") && dcj["roc_tags"].is_array()) {
                for (auto &entry : dcj["roc_tags"]) {
                    int tag   = std::stoi(entry.at("tag").get<std::string>(), nullptr, 16);
                    int crate = entry.at("crate").get<int>();
                    roc_to_crate[tag] = crate;
                }
            }
        }
    }

    // ── Histograms (heap-only, detached from ROOT directory system) ───────────
    TH1F *mod_lms  [N_W];
    TH1F *ref_lms  [N_LMS];
    TH1F *ref_alpha[N_LMS];
    for (int i = 0; i < N_W;   ++i) mod_lms[i]   = makeH(Form("mod_lms_%d",   i + 1));
    for (int i = 0; i < N_LMS; ++i) ref_lms[i]   = makeH(Form("ref_lms_%d",   i + 1));
    for (int i = 0; i < N_LMS; ++i) ref_alpha[i]  = makeH(Form("ref_alpha_%d", i + 1));

    // ── EVIO reading loop ─────────────────────────────────────────────────────
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);

    fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
    fdec::PulseTemplateStore tpl_store;
    if (daq_cfg.wave_cfg.nnls_deconv.enabled &&
        !daq_cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        tpl_store.LoadFromFile(
            db_dir + "/" + daq_cfg.wave_cfg.nnls_deconv.template_file,
            daq_cfg.wave_cfg);
    }
    ana.SetTemplateStore(&tpl_store);

    // Heap-allocate the large EventData struct (following codebase convention).
    auto event = std::make_unique<fdec::EventData>();
    fdec::WaveResult wres;

    long total = 0, n_lms_ev = 0, n_alpha_ev = 0;
    bool done_early = false;

    for (const auto &evio : evio_files) {
        if (done_early) break;
        if (ch.OpenAuto(evio) != evc::status::success) {
            std::cerr << "Cannot open " << evio << "\n";
            continue;
        }
        std::cout << "  Processing " << evio << "\n";

        while (ch.Read() == evc::status::success) {
            if (done_early) break;
            if (!ch.Scan()) continue;
            if (ch.GetEventType() != evc::EventType::Physics) continue;

            for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
                if (max_events > 0 && total >= max_events) { done_early = true; break; }

                event->clear();
                if (!ch.DecodeEvent(ie, *event, nullptr)) continue;

                // ── Pass 1: count channels with data (for classification) ──────
                int nch = 0;
                for (int r = 0; r < event->nrocs; ++r) {
                    auto &roc = event->rocs[r];
                    if (!roc.present) continue;
                    for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                        if (!roc.slots[s].present) continue;
                        for (int c = 0; c < 16; ++c) {
                            if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                            if (roc.slots[s].channels[c].nsamples > 0) ++nch;
                        }
                    }
                }

                bool trig_lms      = (event->info.trigger_bits & prad2::TBIT_lms)   != 0;
                bool trig_alpha    = (event->info.trigger_bits & prad2::TBIT_alpha) != 0;

                // Classify by channel count (trigger_bits unreliable):
                //   LMS flash illuminates all ~1735 channels -> nch > 1000
                //   Alpha source fires only a few crystals  -> 1 <= nch < 50
                const bool is_lms   = (nch > 1000) && trig_lms;
                const bool is_alpha = (nch >= 1 && nch < 50) && trig_alpha;
                if (!is_lms && !is_alpha) continue;

                if (is_lms)   ++n_lms_ev;
                if (is_alpha) ++n_alpha_ev;

                // ── Pass 2: decode waveforms, fill histograms ─────────────────
                for (int r = 0; r < event->nrocs; ++r) {
                    auto &roc = event->rocs[r];
                    if (!roc.present) continue;
                    auto cit   = roc_to_crate.find(roc.tag);
                    int  crate = (cit == roc_to_crate.end()) ? (int)roc.tag : cit->second;

                    for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                        if (!roc.slots[s].present) continue;
                        for (int c = 0; c < 16; ++c) {
                            if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                            auto &cd = roc.slots[s].channels[c];
                            if (cd.nsamples <= 0) continue;

                            int mid = replay.moduleID(crate, s, c);
                            if (mid < 0) continue;
                            auto mtyp = replay.moduleType(crate, s, c);
                            if (mtyp != prad2::MOD_PbWO4 && mtyp != prad2::MOD_LMS) continue;

                            ana.SetChannelKey(roc.tag, s, c);
                            ana.Analyze(cd.samples, cd.nsamples, wres);
                            if (wres.npeaks != 1) continue;
                            float integral = wres.peaks[0].integral;
                            if (integral <= 0.f) continue;

                            if (mtyp == prad2::MOD_PbWO4) {
                                int wid = mid - W_ID_BASE;   // 1..N_W
                                if (wid >= 1 && wid <= N_W && is_lms)
                                    mod_lms[wid - 1]->Fill(integral);
                            } else {
                                int lid = mid - LMS_ID_BASE; // 0..N_LMS-1
                                if (lid >= 0 && lid < N_LMS) {
                                    if (is_lms)   ref_lms[lid]->Fill(integral);
                                    if (is_alpha) ref_alpha[lid]->Fill(integral);
                                }
                            }
                        }
                    }
                }
                ++total;
                if (total % 1000 == 0)
                    std::cerr << "\r  Events: " << total
                              << "  LMS: " << n_lms_ev
                              << "  alpha: " << n_alpha_ev << std::flush;
            }
        }
        ch.Close();
    }
    std::cerr << "\n  Total: " << total
              << "  LMS events: " << n_lms_ev
              << "  alpha events: " << n_alpha_ev << "\n";

    // ── Save raw histograms to ROOT file (before fitting) ─────────────────────
    // Write first so the TFile never sees post-fit TF1 functions,
    // and so closing the file cannot affect in-memory histogram ownership.
    {
        TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
        if (!outfile || !outfile->IsOpen()) {
            std::cerr << "Warning: cannot create " << output_root << "\n";
        } else {
            outfile->mkdir("ref");
            outfile->cd("ref");
            for (int i = 0; i < N_LMS; ++i) {
                ref_lms[i]->Write();
                ref_alpha[i]->Write();
            }
            outfile->mkdir("modules");
            outfile->cd("modules");
            for (int i = 0; i < N_W; ++i)
                if (mod_lms[i]->GetEntries() > 0)
                    mod_lms[i]->Write();
            delete outfile;   // Close() + free; histograms stay in memory (SetDirectory(nullptr))
            std::cout << "Histograms saved to " << output_root << "\n";
        }
    }
    gROOT->cd();  // restore global directory after TFile scope

    // ── Fit reference PMTs ────────────────────────────────────────────────────
    prad2::FitResult fit_lms[N_LMS], fit_alpha[N_LMS];
    for (int i = 0; i < N_LMS; ++i) {
        fit_lms[i]   = prad2::gain_hist_fitter(ref_lms[i],   0.1f);
        fit_alpha[i] = prad2::gain_hist_fitter(ref_alpha[i], 0.1f);
    }

    // ── Write .dat file ───────────────────────────────────────────────────────
    std::ofstream out(output_dat);
    if (!out.is_open()) {
        std::cerr << "Cannot open output file: " << output_dat << "\n";
        return 1;
    }

    out << std::left
        << std::setw(8)  << "Name"
        << std::setw(12) << "lms_peak"
        << std::setw(12) << "lms_sigma"
        << std::setw(14) << "lms_chi2/ndf"
        << std::setw(16) << "alpha_peak (g1)"
        << std::setw(16) << "alpha_sigma (g2)"
        << std::setw(16) << "alpha_chi2/ndf (g3)"
        << "\n";

    for (int i = 0; i < N_LMS; ++i) {
        out << std::setw(8)  << ("LMS" + std::to_string(i + 1))
            << std::setw(12) << std::fixed << std::setprecision(3) << fit_lms[i].mean
            << std::setw(12) << std::fixed << std::setprecision(3) << fit_lms[i].sigma
            << std::setw(14) << std::fixed << std::setprecision(3) << fit_lms[i].chi2pndf
            << std::setw(16) << std::fixed << std::setprecision(3) << fit_alpha[i].mean
            << std::setw(16) << std::fixed << std::setprecision(3) << fit_alpha[i].sigma
            << std::setw(16) << std::fixed << std::setprecision(3) << fit_alpha[i].chi2pndf
            << "\n";
    }

    int n_mod = 0;
    for (int i = 0; i < N_W; ++i) {
        if (mod_lms[i]->GetEntries() < 10) continue;
        auto fr = prad2::gain_hist_fitter(mod_lms[i], 0.1f);

        // g[j] = W_peak * alpha_ref[j] / lms_ref[j]
        float g[N_LMS] = {};
        for (int j = 0; j < N_LMS; ++j) {
            if (fit_lms[j].mean > 0.f && fr.mean > 0.f && fit_alpha[j].mean > 0.f)
                g[j] = fr.mean * fit_alpha[j].mean / fit_lms[j].mean;
        }

        out << std::setw(8)  << ("W" + std::to_string(i + 1))
            << std::setw(12) << std::fixed << std::setprecision(3) << fr.mean
            << std::setw(12) << std::fixed << std::setprecision(3) << fr.sigma
            << std::setw(14) << std::fixed << std::setprecision(3) << fr.chi2pndf
            << std::setw(16) << std::fixed << std::setprecision(3) << g[0]
            << std::setw(16) << std::fixed << std::setprecision(3) << g[1]
            << std::setw(16) << std::fixed << std::setprecision(3) << g[2]
            << "\n";
        ++n_mod;
    }
    out.close();
    std::cout << "Fit results written to " << output_dat
              << " (" << n_mod << " W modules)\n";

    // Cleanup
    for (int i = 0; i < N_W;   ++i) delete mod_lms[i];
    for (int i = 0; i < N_LMS; ++i) { delete ref_lms[i]; delete ref_alpha[i]; }

    return 0;
}
