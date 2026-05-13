// epCalib.cpp: multi-threaded version of epCalib.cpp
// Files are processed in rounds; each round assigns one file per thread.
// Each thread opens its own TFile, runs HyCal clustering, and fills per-module
// energy histograms independently (reused across rounds).
// After all rounds finish, histograms are merged and the calibration fitting
// is performed single-threaded.
//=============================================================================
//
// Usage: epCalib <input_raw.root|dir> [more files/dirs...]
//                  [-i iteration] [-o output_root_file]
//                  [-E Ebeam] [-D daq_config.json] [-n max_events]
//                  [-j num_threads] [-f (use firmware peaks)]
//   - input_raw.root|dir: input ROOT file(s) or directory with *_raw.root files
//   - iteration: calibration iteration (default: 1)
//   - output_root_file: output ROOT file (default: auto from db_dir)
//   - Ebeam: beam energy in MeV (default: 2100)
//   - daq_config.json: DAQ config file (default: db_dir/daq_config.json)
//   - max_events: max total events to process (default: all)
//   - num_threads: number of parallel threads (default: 4)
//   - -f: use firmware peak analysis instead of DAQ-mode peaks (default: false)
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "HyCalTimeCuts.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "load_daq_config.h"
#include "RunInfoConfig.h"
#include "gain_factor.h"

#include <TFile.h>
#include <TTree.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TClass.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <thread>
#include <mutex>
#include <memory>
#include <algorithm>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars = prad2::RawEventData;
using namespace analysis;

// ── Per-thread accumulated results ──────────────────────────────────────────
struct ThreadResult {
    fdec::HyCalSystem                        hycal;
    std::vector<std::unique_ptr<TH1F>>       module_hists;   // indexed by module index
    std::unique_ptr<TH2F>                    h2_energy_theta;
    std::unique_ptr<TH2F>                    hit_pos;
    std::unique_ptr<TH1F>                    h_E_1cl;
    long long                                events_processed = 0;
};

// ── File collection helper ───────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_raw.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    // ROOT multi-thread safety (must be called before any ROOT object creation)
    ROOT::EnableThreadSafety();
    // Force dictionary loading in main thread
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");
    TClass::GetClass("TH1F");
    TClass::GetClass("TH2F");

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_root_file, daq_config_file;
    int  iteration   = 1;
    int  max_events  = -1;
    int  num_threads = 4;
    float Ebeam      = 2100.f;
    float hycal_z    = 6269.f;
    bool firmware_peaks = false;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:E:D:n:j:f")) != -1) {
        switch (opt) {
            case 'i': iteration        = std::atoi(optarg); break;
            case 'o': output_root_file = optarg; break;
            case 'E': Ebeam            = std::atof(optarg); break;
            case 'D': daq_config_file  = optarg; break;
            case 'f': firmware_peaks   = true; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'j': num_threads      = std::atoi(optarg); break;
        }
    }

    // Collect all input files
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: epCalib <input_raw.root|dir> [more...] "
                     "[-i iter] [-o out.root] [-E Ebeam] [-D daq.json] "
                     "[-n max_events] [-j threads] [-f (use firmware peaks)]\n";
        return 1;
    }

    // ── Run number / output paths ─────────────────────────────────────────────
    std::string run_str = "unknown";
    {
        std::string fname = fs::path(root_files[0]).filename().string();
        auto ppos = fname.find("prad_");
        if (ppos != std::string::npos) {
            size_t s = ppos + 5, e = s;
            while (e < fname.size() && std::isdigit((unsigned char)fname[e])) e++;
            if (e > s) run_str = std::to_string(std::stoul(fname.substr(s, e - s)));
        }
    }
    std::string run_out_dir = "Physics_calib/" + run_str;
    fs::create_directories(run_out_dir);
    std::cerr << "Output directory: " << run_out_dir << "\n";

    std::string input_calib_file, output_calib_file;
    if (iteration == 1)
        input_calib_file = db_dir + "/calibration/calibration_factor_2_0.json";
    else if (iteration > 1)
        input_calib_file = run_out_dir + Form("/calib_iter%d.json", iteration - 1);
    else {
        std::cerr << "Invalid iteration number: " << iteration << ". Must be >= 1.\n";
        return 1;
    }
    output_calib_file = run_out_dir + Form("/calib_iter%d.json", iteration);

    if (output_root_file.empty())
        output_root_file = run_out_dir + Form("/CalibResult_iter%d.root", iteration);

    if (daq_config_file.empty())
        daq_config_file = db_dir + "/daq_config.json";

    // ── Thread count ─────────────────────────────────────────────────────────
    int n_files    = static_cast<int>(root_files.size());
    num_threads    = std::max(1, std::min(num_threads, n_files));
    int num_rounds = (n_files + num_threads - 1) / num_threads;
    std::cout << "Processing " << n_files << " file(s) with "
              << num_threads << " thread(s), " << num_rounds << " round(s)\n";

    // ── Initialize per-thread results (once, reused across rounds) ───────────
    std::vector<std::unique_ptr<ThreadResult>> results(num_threads);
    std::mutex io_mtx;

    for (int tid = 0; tid < num_threads; ++tid) {
        auto res = std::make_unique<ThreadResult>();
        // Each thread initializes its own HyCalSystem
        res->hycal.Init(db_dir + "/hycal_map.json");
        int nmatched = res->hycal.LoadCalibration(input_calib_file);
        std::cerr << "[thread " << tid << "] calibration: "
                  << input_calib_file << " (" << nmatched << " modules)\n";

        int nmod_t = res->hycal.module_count();
        res->module_hists.resize(nmod_t);
        for (int i = 0; i < nmod_t; ++i) {
            const auto &mod = res->hycal.module(i);
            res->module_hists[i] = std::make_unique<TH1F>(
                Form("h_%s_%d", mod.name.c_str(), tid),
                Form("%s cluster energy;Energy (MeV);Counts", mod.name.c_str()),
                250, 0, 5000);
            res->module_hists[i]->SetDirectory(nullptr);
        }
        res->h2_energy_theta = std::make_unique<TH2F>(
            Form("h2_energy_theta_%d", tid),
            "Energy vs Theta;Theta (deg);Energy (MeV)",
            80, 0, 8, 2000, 0, 4000);
        res->h2_energy_theta->SetDirectory(nullptr);
        res->hit_pos = std::make_unique<TH2F>(
            Form("hit_pos_%d", tid),
            "One Cluster Event Hit positions;hycal X (mm);hycal Y (mm)",
            250, -500, 500, 250, -500, 500);
        res->hit_pos->SetDirectory(nullptr);
        res->h_E_1cl = std::make_unique<TH1F>(
            Form("one_cluster_energy_%d", tid),
            "Single-cluster Event energy;E (MeV);Counts",
            1000, 0, 5000);
        res->h_E_1cl->SetDirectory(nullptr);
        results[tid] = std::move(res);
    }

    // ── Process files in rounds: num_threads files per round, 1 file/thread ──
    static constexpr uint32_t TBIT_sum = (1u << 8);
    static constexpr uint32_t TBIT_lms = (1u << 24);

    for (int round = 0; round < num_rounds; ++round) {
        int round_start        = round * num_threads;
        int round_end          = std::min(round_start + num_threads, n_files);
        int threads_this_round = round_end - round_start;

        std::cout << "\nRound " << (round + 1) << "/" << num_rounds
                  << ": files [" << round_start << ", " << round_end - 1 << "]\n";

        std::vector<std::thread> threads;
        threads.reserve(threads_this_round);

        for (int t = 0; t < threads_this_round; ++t) {
            threads.emplace_back([&, t, round]() {
                int fi = round * num_threads + t;
                ThreadResult *res = results[t].get();

                // Open a single file (no TChain)
                TFile *rfile = TFile::Open(root_files[fi].c_str(), "READ");
                if (!rfile || rfile->IsZombie()) {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "[thread " << t << "] cannot open "
                              << root_files[fi] << "\n";
                    if (rfile) { rfile->Close(); delete rfile; }
                    return;
                }
                TTree *tree = dynamic_cast<TTree *>(rfile->Get("events"));
                if (!tree) {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "[thread " << t << "] no 'events' tree in "
                              << root_files[fi] << "\n";
                    rfile->Close(); delete rfile;
                    return;
                }

                long long nentries = tree->GetEntries();
                if (max_events > 0) {
                    // Per-file cap: distribute max_events evenly across all files.
                    // Use ceiling division so small max_events still processes >= 1 event/file.
                    long long cap = ((long long)max_events + n_files - 1) / n_files;
                    nentries = std::min(nentries, cap);
                }

                EventVars ev;
                prad2::SetRawReadBranches(tree, ev);

                int run_num = get_run_int(root_files[fi]);
                auto localConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
                //auto gain_correction = prad2::ComputeGainCorrection(localConfig.gain_data_dir, run_num, localConfig.gain_ref_run);

                // Per-module HyCal time-cut table.  When the runinfo entry
                // doesn't reference a per-module file the table is uniform.
                std::string hc_time_path = localConfig.hycal_time_cut_file.empty()
                    ? std::string()
                    : db_dir + "/" + localConfig.hycal_time_cut_file;
                auto hc_time_cuts = prad2::LoadHyCalTimeCuts(
                    hc_time_path, res->hycal,
                    localConfig.hc_time_win_lo, localConfig.hc_time_win_hi);

                fdec::HyCalCluster clusterer(res->hycal);
                fdec::ClusterConfig cl_cfg;
                clusterer.SetConfig(cl_cfg);

                for (long long i = 0; i < nentries; ++i) {
                    tree->GetEntry(i);

                    if (i % 5000 == 0) {
                        std::lock_guard<std::mutex> lk(io_mtx);
                        std::cout << "[thread " << t << "] file[" << fi << "] event "
                                  << i + 1 << "/" << nentries << "\r" << std::flush;
                    }

                    if (!(ev.trigger_bits & TBIT_sum)) continue;
                    if (  ev.trigger_bits & TBIT_lms ) continue;

                    if (ev.nch > 500) continue;

                    clusterer.Clear();
                    for (int j = 0; j < ev.nch; ++j) {
                        const auto *mod = res->hycal.module_by_id(ev.module_id[j]);
                        if (!mod || !mod->is_hycal()) continue;

                        float adc = 0.f;

                        if(!firmware_peaks){
                            const auto hc_win = hc_time_cuts.at(mod->index);
                            int bestIdx = -1;
                            float bestHeight = -1.f;
                            for(int p = 0; p < ev.npeaks[j]; ++p){
                                if(ev.peak_time[j][p] > hc_win.lo &&
                                    ev.peak_time[j][p] < hc_win.hi) {
                                    if(ev.peak_integral[j][p] > bestHeight) {
                                        bestHeight = ev.peak_integral[j][p];
                                        bestIdx = p;
                                    }
                                }
                            }
                            if (bestIdx < 0) continue;
                            adc    = ev.peak_integral[j][bestIdx];
                            int mod_id = ev.module_id[j];
                            // gain correction for HyCal modules
                            //if(mod_id>1000) adc *= gain_correction.w[mod_id-1000].corr[1]; // Use g2-based correction for PbWO4 (matches LMS2)
                            //else adc *= gain_correction.g[mod_id].corr[1]; // Use g2-based correction for PbGlass (matches LMS2)
                            if(mod_id > 1000) adc *= ev.gain_factor[j];
                            else adc *= ev.gain_factor[j];
                        }
                        else{
                            if(ev.daq_npeaks[j] <= 0) continue;
                            adc = ev.daq_peak_integral[j][0]; // take the first DAQ peak (should be only one)
                        }

                        float energy = (mod->cal_factor > 0)
                            ? static_cast<float>(mod->energize(adc))
                            : adc * 0.f;
                        clusterer.AddHit(mod->index, energy, 0.f); // time info not available in raw branches, set to 0
                    }

                    clusterer.FormClusters();
                    std::vector<fdec::ClusterHit> hits;
                    clusterer.ReconstructHits(hits);

                    if (hits.size() == 1 && hits[0].nblocks > 5) {
                        int midx = res->hycal.id_to_index(hits[0].center_id);
                        if (midx >= 0 && midx < (int)res->module_hists.size())
                            res->module_hists[midx]->Fill(hits[0].energy);
                        res->hit_pos->Fill(hits[0].x, hits[0].y);
                        res->h_E_1cl->Fill(hits[0].energy);
                        float theta = std::atan(std::sqrt(hits[0].x * hits[0].x +
                                                          hits[0].y * hits[0].y) / hycal_z)
                                      * 180.f / 3.14159265f;
                        res->h2_energy_theta->Fill(theta, hits[0].energy);
                    }
                }

                res->events_processed += nentries;
                rfile->Close(); delete rfile;
                {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cout << "\n[thread " << t << "] done: " << root_files[fi]
                              << " (" << nentries << " events)\n";
                }
            });
        }

        for (auto &th : threads)
            th.join();
        std::cout << "Round " << (round + 1) << " complete.\n";
    }

    // ── Merge histograms (single-threaded) ────────────────────────────────────
    std::cout << "\nAll rounds finished. Merging histograms...\n";

    // Initialize main HyCalSystem for the fitting stage
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    int nmatched = hycal.LoadCalibration(input_calib_file);
    std::cerr << "Main: calibration loaded (" << nmatched << " modules)\n";

    analysis::PhysicsTools physics(hycal);

    // Aggregate: merge per-module energy histograms from each thread
    int nmod = hycal.module_count();
    for (int m = 0; m < nmod; ++m) {
        int mod_id = hycal.module(m).id;
        TH1F *main_h = physics.GetModuleEnergyHist(mod_id);
        if (!main_h) continue;
        for (int tid = 0; tid < num_threads; ++tid) {
            if (!results[tid]) continue;
            if (m < (int)results[tid]->module_hists.size() && results[tid]->module_hists[m])
                main_h->Add(results[tid]->module_hists[m].get());
        }
    }

    // Merge E-vs-theta 2D histogram
    TH2F *main_etheta = physics.GetEnergyVsThetaHist();
    for (int tid = 0; tid < num_threads; ++tid) {
        if (!results[tid]) continue;
        if (main_etheta && results[tid]->h2_energy_theta)
            main_etheta->Add(results[tid]->h2_energy_theta.get());
    }

    // Merge global histograms
    TH2F *hit_pos = new TH2F("hit_pos",
        "One Cluster Event Hit positions;hycal X (mm);hycal Y (mm)",
        250, -500, 500, 250, -500, 500);
    hit_pos->SetDirectory(nullptr);
    TH1F *h_E_1cl = new TH1F("one_cluster_energy",
        "Single-cluster Event energy;E (MeV);Counts",
        1000, 0, 5000);
    h_E_1cl->SetDirectory(nullptr);

    long long total_events = 0;
    for (int tid = 0; tid < num_threads; ++tid) {
        if (!results[tid]) continue;
        if (results[tid]->hit_pos) hit_pos->Add(results[tid]->hit_pos.get());
        if (results[tid]->h_E_1cl) h_E_1cl->Add(results[tid]->h_E_1cl.get());
        total_events += results[tid]->events_processed;
    }
    std::cerr << "Total events processed: " << total_events << "\n";

    // ── Single-threaded analysis: fit peaks, update calibration ──────────────
    std::cerr << "Fitting peaks and calculating calibration constants...\n";

    TH1F *h_mearured_peak = new TH1F("measured_peak",
        "Measured Peak Position;Energy (MeV);Counts", 1000, 0, 5000);
    TH1F *h_recon_sigma = new TH1F("recon_sigma",
        "Reconstructed Cluster Energy Resolution;Sigma (MeV);Counts", 100, 0, 200);
    TH1F *h_recon_chi2 = new TH1F("recon_chi2/ndf",
        "Reconstructed E hist Fit Chi2/ndf;Chi2/ndf;Counts", 100, 0, 50);
    TH1F *ratio_module_all = new TH1F("ratio_all",
        "Ratio of Expected/Measured Peak Position for All Modules;Ratio;Modules", 200, 0, 4);
    TH2F *module_ratio = new TH2F("#cbar#bar{E_{recon}} - E_{expect}#cbar #/ E_{expect}",
        "#cbar#bar{E_{recon}} - E_{expect}#cbar #/ E_{expect}",
        34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);

    std::string dat_out_path = run_out_dir + Form("/fitting_parameters_iter%d.dat", iteration);
    std::ofstream dat_out(dat_out_path);
    if (!dat_out.is_open()) {
        std::cerr << "Cannot open output file " << dat_out_path << "\n";
        return 1;
    }

    dat_out << std::left;
    dat_out << std::setw(8)  << "Module"
            << std::setw(16) << "ExpectedPeak"
            << std::setw(16) << "MeasuredPeak"
            << std::setw(16) << "oldFactor"
            << std::setw(16) << "Ratio"
            << std::setw(16) << "Sigma"
            << std::setw(16) << "Chi2/ndf" << "\n";

    int n_calibrated = 0;
    for (int m = 0; m < nmod; ++m) {
        const std::string &name = hycal.module(m).name;
        if (name[0] != 'W') continue;
        int mod_id = hycal.module(m).id;
        auto [peak, sigma, chi2] = physics.FitPeakResolution(mod_id);
        if (peak <= 0 || sigma <= 0 || sigma > 5 * 0.026*peak || chi2 >= 2.f) {
            std::cout << "Check!!! Module " << hycal.module(m).name
                 << ": fit failed (peak=" << peak
                 << ", sigma=" << sigma
                 << ", chi2/ndf=" << chi2 << ")\n";
        }

        if(physics.GetModuleEnergyHist(mod_id)->GetEntries() < 1.) continue; // skip modules with no entries
        if(peak <= 0) peak = physics.GetModuleEnergyHist(mod_id)->GetMean(); // fallback to mean if fit failed
        if(peak <= 0) continue; // still no valid peak after fallback — skip to avoid Inf factor
        
        float theta_deg = std::atan(std::sqrt(hycal.module(m).x * hycal.module(m).x +
                                                hycal.module(m).y * hycal.module(m).y)
                                    / hycal_z) * 180.f / 3.14159265f;
        float expected_peak = physics.ExpectedEnergy(theta_deg, Ebeam, "ep");
        float ratio         = expected_peak / peak;
        if(ratio < 0.5) ratio = 0.5; 
        if(ratio > 2.0) ratio = 2.0;
        ratio_module_all->Fill(ratio);

        double current_factor = hycal.GetCalibConstant(mod_id);
        double new_factor = current_factor * ratio;
        hycal.SetCalibConstant(mod_id, new_factor);
        hycal.SetCalibBaseEnergy(mod_id, expected_peak);

        module_ratio->Fill(hycal.module(m).x, hycal.module(m).y,
                            std::abs(1.f - 1.f / ratio));

        h_mearured_peak->Fill(peak);
        h_recon_sigma->Fill(sigma);
        h_recon_chi2->Fill(chi2);

        dat_out << std::setw(8)  << name
                << std::setw(16) << expected_peak
                << std::setw(16) << peak
                << std::setw(16) << current_factor
                << std::setw(16) << ratio
                << std::setw(16) << sigma
                << std::setw(16) << chi2 << "\n";
        n_calibrated++;
    }
    std::cerr << "Calibrated " << n_calibrated << " modules. "
              << "Results written to " << dat_out_path << "\n";

    // Write new calibration constants
    hycal.PrintCalibConstants(output_calib_file);

    // ── Write output ROOT file ────────────────────────────────────────────────
    TFile outfile(output_root_file.c_str(), "RECREATE");

    TCanvas *c = new TCanvas("c", "Calibration", 1200, 1200);
    module_ratio->SetStats(0);
    module_ratio->Draw("COLZ");
    module_ratio->Write();

    TLatex t;
    t.SetTextSize(0.0122);
    t.SetTextColor(kBlack);
    for (int m = 0; m < nmod; ++m) {
        const std::string &name = hycal.module(m).name;
        if (name[0] != 'W') continue;
        t.DrawLatex(hycal.module(m).x - 6., hycal.module(m).y - 2., name.c_str());
    }
    c->Update();
    c->Write();

    TCanvas *c_hit = new TCanvas("c_hit", "Hit Position", 1200, 1200);
    hit_pos->Draw("COLZ");
    hit_pos->Write();
    outfile.cd();
    c_hit->Write();

    physics.GetEnergyVsThetaHist()->Write();
    h_E_1cl->Write();
    h_mearured_peak->Write();
    h_recon_sigma->Write();
    h_recon_chi2->Write();
    ratio_module_all->Write();

    physics.FillNeventsModuleMap();
    TH2F *h_map = physics.GetNeventsModuleMapHist();
    TCanvas *c_map = new TCanvas("c_map", "Number of Events per Module", 1200, 1200);
    h_map->Draw("COLZ");
    c_map->Write();
    h_map->Write();

    outfile.mkdir("module_energy");
    outfile.cd("module_energy");
    for (int i = 0; i < nmod; ++i) {
        int mod_id = hycal.module(i).id;
        TH1F *h = physics.GetModuleEnergyHist(mod_id);
        if (h && h->GetEntries() > 0) h->Write();
    }

    outfile.cd();
    outfile.Close();
    std::cout << "Output written to " << output_root_file << "\n";

    return 0;
}
