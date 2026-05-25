
//=============================================================================
// det_calib — detector position calibration via Møller scattering
//
// Phase 1: replay EVIO files (multi-threaded) → per-file *_recon.root
// Phase 2: chain all *_recon.root, select Møller events, fit vertex-z and
//          beam-center distributions for HyCal and each GEM plane.
//
// Usage:
//   det_calib <evio_file_or_dir> [more files/dirs...]
//             -o output_dir
//             [-f max_files] [-j num_threads] [-n max_events]
//             [-c daq_config.json] [-d hycal_map.json]
//             [-g gem_pedestal_file] [-z zerosup_threshold]
//             [--hx val] [--hy val]
//             [--g0x val] [--g0y val] ... [--g3x val] [--g3y val]
//             [-s]  save intermediate *_recon.root files (default: delete)
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TClass.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TF1.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

// ── forward declarations ──────────────────────────────────────────────────
static std::vector<std::string> collectEvioFiles(const std::string &path);
static std::string makeOutputFile(const std::string &evio_path);
static void analyzeMollers(const std::vector<std::string> &recon_files,
                            const std::string &output,
                            const std::string &run_str,
                            RunConfig &geo,
                            int max_events,
                            float HyCal_shift_x, float HyCal_shift_y,
                            const float GEM_shift_x[4], const float GEM_shift_y[4]);
static float fitAndDraw(TH1F *hist, const std::string &out_path,
                        float survey_position);

// ── file helpers ──────────────────────────────────────────────────────────
static std::vector<std::string> collectEvioFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &e : fs::directory_iterator(path)) {
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

static std::string makeOutputFile(const std::string &evio_path)
{
    std::string out = fs::path(evio_path).filename().string();
    auto pos = out.find(".evio");
    if (pos != std::string::npos)
        out = out.substr(0, pos) + out.substr(pos + 5);
    out += "_recon.root";
    return out;
}

// ── Main ──────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string daq_config, daq_map, gem_ped_file, output_dir, run_config;
    float zerosup_override = 0.f;
    int  max_files   = -1;
    int  num_threads = 4;
    int  max_events  = -1;
    bool save_recon  = false;

    float HyCal_shift_x    = 0.f;
    float HyCal_shift_y    = 0.f;
    float GEM_shift_x[4]   = {0.f, 0.f, 0.f, 0.f};
    float GEM_shift_y[4]   = {0.f, 0.f, 0.f, 0.f};

    enum {
        OPT_HX=256, OPT_HY,
        OPT_G0X, OPT_G0Y, OPT_G1X, OPT_G1Y,
        OPT_G2X, OPT_G2Y, OPT_G3X, OPT_G3Y
    };
    static const struct option long_opts[] = {
        {"hx",  required_argument, nullptr, OPT_HX},
        {"hy",  required_argument, nullptr, OPT_HY},
        {"g0x", required_argument, nullptr, OPT_G0X},
        {"g0y", required_argument, nullptr, OPT_G0Y},
        {"g1x", required_argument, nullptr, OPT_G1X},
        {"g1y", required_argument, nullptr, OPT_G1Y},
        {"g2x", required_argument, nullptr, OPT_G2X},
        {"g2y", required_argument, nullptr, OPT_G2Y},
        {"g3x", required_argument, nullptr, OPT_G3X},
        {"g3y", required_argument, nullptr, OPT_G3Y},
        {nullptr, 0, nullptr, 0}
    };

    int longidx = 0, opt;
    while ((opt = getopt_long(argc, argv, "o:f:n:j:c:C:d:g:z:s", long_opts, &longidx)) != -1) {
        switch (opt) {
            case 'o': output_dir       = optarg; break;
            case 'f': max_files        = std::atoi(optarg); break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'j': num_threads      = std::atoi(optarg); break;
            case 'c': run_config       = optarg; break;
            case 'C': daq_config       = optarg; break;
            case 'd': daq_map          = optarg; break;
            case 'g': gem_ped_file     = optarg; break;
            case 'z': zerosup_override = std::atof(optarg); break;
            case 's': save_recon       = true; break;
            case OPT_HX:  HyCal_shift_x  = std::atof(optarg); break;
            case OPT_HY:  HyCal_shift_y  = std::atof(optarg); break;
            case OPT_G0X: GEM_shift_x[0] = std::atof(optarg); break;
            case OPT_G0Y: GEM_shift_y[0] = std::atof(optarg); break;
            case OPT_G1X: GEM_shift_x[1] = std::atof(optarg); break;
            case OPT_G1Y: GEM_shift_y[1] = std::atof(optarg); break;
            case OPT_G2X: GEM_shift_x[2] = std::atof(optarg); break;
            case OPT_G2Y: GEM_shift_y[2] = std::atof(optarg); break;
            case OPT_G3X: GEM_shift_x[3] = std::atof(optarg); break;
            case OPT_G3Y: GEM_shift_y[3] = std::atof(optarg); break;
        }
    }

    // collect input EVIO files (files, directories, or mixed)
    std::vector<std::string> evio_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectEvioFiles(argv[i]);
        evio_files.insert(evio_files.end(), f.begin(), f.end());
    }

    if (evio_files.empty() || output_dir.empty()) {
        std::cerr <<
            "Usage: det_calib <evio_file_or_dir> [more files/dirs...] -o output_dir\n"
            "       [-f max_files] [-j threads] [-n max_events]\n"
            "       [-C daq_config.json] [-d hycal_map.json]\n"
            "       [-c run_config.json] [-g gem_ped_file] [-z zerosup]\n"
            "       [-s] (save intermediate *_recon.root; default: delete)\n"
            "       [--hx val] [--hy val]                  HyCal shift (mm)\n"
            "       [--g0x val] [--g0y val] .. [--g3x val] [--g3y val]  GEM shift (mm)\n";
        return 1;
    }

    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (daq_config.empty()) daq_config = db_dir + "/daq_config.json";
    if (daq_map.empty())    daq_map    = db_dir + "/hycal_map.json";
    if (run_config.empty()) run_config = db_dir + "/runinfo/general.json";

    int run_num = get_run_int(evio_files[0]);
    std::string run_str = get_run_str(evio_files[0]);
    gRunConfig = LoadRunConfig(run_config, run_num);

    std::cout << "=== Phase 1: EVIO replay ===\n"
              << "  Files   : " << num_files << "\n"
              << "  Threads : " << num_threads << "\n"
              << "  Output  : " << output_dir << "\n"
              << "  DAQ cfg : " << daq_config << "\n"
              << "  HyCal   : " << daq_map << "\n";

    // ── Phase 1: multi-threaded EVIO → _recon.root ───────────────────────
    std::atomic<int>         next_file{0};
    std::mutex               io_mtx;
    std::atomic<int>         replay_errors{0};
    std::vector<std::string> recon_out_files;
    std::mutex               out_files_mtx;

    auto worker = [&]() {
        analysis::Replay replay;
        if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
        replay.LoadHyCalMap(daq_map);

        while (true) {
            int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            std::string out = output_dir + "/" + makeOutputFile(evio_files[idx]);
            bool ok = replay.ProcessWithRecon(evio_files[idx], out, gRunConfig, db_dir,
                                              daq_config, gem_ped_file, zerosup_override);

            std::lock_guard<std::mutex> lk(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                          << evio_files[idx] << " -> " << out << "\n";
                std::lock_guard<std::mutex> lk2(out_files_mtx);
                recon_out_files.push_back(out);
            } else {
                std::cerr << "  [" << (idx + 1) << "/" << num_files << "] FAILED: "
                          << evio_files[idx] << "\n";
                replay_errors++;
            }
        }
    };

    {
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (int i = 0; i < num_threads; ++i)
            threads.emplace_back(worker);
        for (auto &t : threads) t.join();
    }

    std::cout << "Phase 1 done: " << num_files << " files, "
              << replay_errors.load() << " errors\n";

    if (recon_out_files.empty()) {
        std::cerr << "No recon files produced; aborting.\n";
        return replay_errors.load() > 0 ? 1 : 0;
    }
    std::sort(recon_out_files.begin(), recon_out_files.end());

    // ── Phase 2: Møller selection and detector calibration ───────────────
    std::cout << "\n=== Phase 2: Møller analysis ===\n";
    std::string calib_out = output_dir + "/prad_" + run_str + "_posCalib.root";
    analyzeMollers(recon_out_files, calib_out, run_str, gRunConfig,
                   max_events, HyCal_shift_x, HyCal_shift_y, GEM_shift_x, GEM_shift_y);

    // ── Clean up intermediate recon files ────────────────────────────────
    if (!save_recon) {
        std::cout << "\n=== Cleaning up intermediate recon files ===\n";
        int n_del = 0;
        for (const auto &f : recon_out_files) {
            std::error_code ec;
            fs::remove(f, ec);
            if (!ec) { ++n_del; }
            else std::cerr << "  Warning: cannot remove " << f << ": " << ec.message() << "\n";
        }
        std::cout << "  Deleted " << n_del << " file(s).\n";
    }

    std::cout << "All done.\n";
    return replay_errors.load() > 0 ? 1 : 0;
}

// ── Phase 2: Møller analysis ──────────────────────────────────────────────
static void analyzeMollers(const std::vector<std::string> &recon_files,
                            const std::string &output,
                            const std::string &run_str,
                            RunConfig &geo,
                            int max_events,
                            float HyCal_shift_x, float HyCal_shift_y,
                            const float GEM_shift_x[4], const float GEM_shift_y[4])
{
    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR) + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // --- setup TChain ---
    TChain chain("recon");
    for (const auto &f : recon_files) {
        chain.Add(f.c_str());
        std::cout << "  + " << f << "\n";
    }
    std::cout << "  Total entries: " << chain.GetEntries() << "\n";

    EventVars_Recon ev;
    prad2::SetReconReadBranches(&chain, ev);

    // --- output histograms ---
    TH1F *vertex_hycal   = new TH1F("vertex_hycal",   "Moller vertex z HyCal;Z (mm);Counts",         600,  5600, 6800);
    TH2F *center_hycal   = new TH2F("center_hycal",   "Moller center HyCal;X (mm);Y (mm)",            200,  -100,  100, 200, -100, 100);
    TH1F *center_hycal_x = new TH1F("center_hycal_x", "Moller center X HyCal;X (mm);Counts",         320,   -20,   20);
    TH1F *center_hycal_y = new TH1F("center_hycal_y", "Moller center Y HyCal;Y (mm);Counts",         320,   -20,   20);
    TH2F *hits_hycal     = new TH2F("hits_hycal",     "Moller hits HyCal;X (mm);Y (mm)",             400,  -400,  400, 400, -400, 400);

    TH1F *vertex_gem[4], *center_gem_x[4], *center_gem_y[4];
    TH2F *center_gem[4], *hits_gem[4];
    for (int d = 0; d < 4; d++) {
        vertex_gem[d]   = new TH1F(Form("vertex_gem%d",   d), Form("Moller vertex z GEM%d;Z (mm);Counts",    d), 500, 5200, 6200);
        center_gem[d]   = new TH2F(Form("center_gem%d",   d), Form("Moller center GEM%d;X (mm);Y (mm)",      d),  200, -100,  100, 200, -100, 100);
        center_gem_x[d] = new TH1F(Form("center_gem_x%d", d), Form("Moller center X GEM%d;X (mm);Counts",    d),  400,  -20,   20);
        center_gem_y[d] = new TH1F(Form("center_gem_y%d", d), Form("Moller center Y GEM%d;Y (mm);Counts",    d),  100,  -20,   20);
        hits_gem[d]     = new TH2F(Form("hits_gem%d",     d), Form("Moller hits GEM%d;X (mm);Y (mm)",        d),  400, -400,  400, 400, -400, 400);
    }

    TH2F *E_theta = new TH2F("E_theta", "Moller E vs #theta;#theta (deg);E (MeV)", 120, 0, 6, int(geo.Ebeam)/2, 0, geo.Ebeam);

    MollerData hycal_mollers;
    MollerData gem_mollers[4];

    // --- event loop: select Møller events ---
    Long64_t N = chain.GetEntries();
    if (max_events > 0 && max_events < (int)N) N = max_events;

    for (Long64_t i = 0; i < N; i++) {
        chain.GetEntry(i);
        if (i % 10000 == 0)
            std::cerr << "\r  Event: " << i << " / " << N << std::flush;
        if (ev.n_clusters != 2) continue;
        if (ev.matchNum != 2) continue;
        if (ev.cl_nblocks[0] < 3 || ev.cl_nblocks[1] < 3) continue;
        if (fabs(ev.cl_x[0]) < 20.75f * 2.5f && fabs(ev.cl_y[0]) < 20.75f * 2.5f) continue;
        if (fabs(ev.cl_x[1]) < 20.75f * 2.5f && fabs(ev.cl_y[1]) < 20.75f * 2.5f) continue;
        if (fabs(ev.cl_x[0]) > 20.75f * 16.f || fabs(ev.cl_y[0]) > 20.75f * 16.f) continue;
        if (fabs(ev.cl_x[1]) > 20.75f * 16.f || fabs(ev.cl_y[1]) > 20.75f * 16.f) continue;

        float Epair = ev.mHit_E[0] + ev.mHit_E[1];
        if (geo.Ebeam <= 0.f) {
            std::cerr << "\nError: Ebeam not set, cannot apply Moller energy cut.\n";
            break;
        }

        MollerEvent h_m(
            {ev.mHit_x[0], ev.mHit_y[0], ev.mHit_z[0], ev.mHit_E[0]},
            {ev.mHit_x[1], ev.mHit_y[1], ev.mHit_z[1], ev.mHit_E[1]});

        float theta1 = atan2(std::sqrt(h_m.first.x*h_m.first.x + h_m.first.y*h_m.first.y), h_m.first.z) * 180.f / M_PI;
        float theta2 = atan2(std::sqrt(h_m.second.x*h_m.second.x + h_m.second.y*h_m.second.y), h_m.second.z) * 180.f / M_PI;
        if (!physics.isMoller_kinematic(theta1, h_m.first.E, theta2, h_m.second.E, geo.Ebeam, 0.035f))
            continue;
        if(fabs(physics.GetMollerPhiDiff(h_m)) > 10.f) continue;

        //add some scattering angle cuts for 0.7GeV
        if(geo.Ebeam > 0.f && geo.Ebeam < 1000.f) {
            if(theta1 < 1.5f || theta2 < 1.5f) continue;
        }

        E_theta->Fill(theta1, h_m.first.E);
        E_theta->Fill(theta2, h_m.second.E);

        hycal_mollers.push_back(h_m);
        hits_hycal->Fill(h_m.first.x,  h_m.first.y);
        hits_hycal->Fill(h_m.second.x, h_m.second.y);

        // upstream GEM: both hits on the same chamber
        if (ev.mHit_gid[0][0] == ev.mHit_gid[1][0]) {
            int det_id = ev.mHit_gid[0][0];
            if (det_id >= 0 && det_id < 4) {
                gem_mollers[det_id].push_back(MollerEvent(
                    {ev.mHit_gx[0][0], ev.mHit_gy[0][0], ev.mHit_gz[0][0], ev.mHit_E[0]},
                    {ev.mHit_gx[1][0], ev.mHit_gy[1][0], ev.mHit_gz[1][0], ev.mHit_E[1]}));
            }
        }

        // downstream GEM: both hits on the same chamber
        if (ev.mHit_gid[0][1] == ev.mHit_gid[1][1]) {
            int det_id = ev.mHit_gid[0][1];
            if (det_id >= 0 && det_id < 4) {
                gem_mollers[det_id].push_back(MollerEvent(
                    {ev.mHit_gx[0][1], ev.mHit_gy[0][1], ev.mHit_gz[0][1], ev.mHit_E[0]},
                    {ev.mHit_gx[1][1], ev.mHit_gy[1][1], ev.mHit_gz[1][1], ev.mHit_E[1]}));
            }
        }
    }
    std::cerr << "\n";

    std::cerr << "Summary of selected Moller events:\n";
    std::cerr << "  HyCal: " << hycal_mollers.size() << "\n";
    for (int d = 0; d < 4; d++)
        std::cerr << "  GEM " << d << ": " << gem_mollers[d].size() << "\n";

    // --- fill histograms ---
    TransformDetData(hycal_mollers, HyCal_shift_x, HyCal_shift_y, 0.f);
    for (size_t i = 0; i < hycal_mollers.size(); i++) {
        vertex_hycal->Fill(physics.GetMollerZdistance(hycal_mollers[i], geo.Ebeam));
        if (i >= 1) {
            auto c = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            if(c[0] == c[1] && c[0] == 0.f) continue; // skip zero-center events (likely bad recon)
            center_hycal->Fill(c[0], c[1]);
            center_hycal_x->Fill(c[0]);
            center_hycal_y->Fill(c[1]);
        }
    }

    for (int d = 0; d < 4; d++) {
        TransformDetData(gem_mollers[d], GEM_shift_x[d], GEM_shift_y[d], 0.f);
        for (size_t i = 0; i < gem_mollers[d].size(); i++) {
            vertex_gem[d]->Fill(physics.GetMollerZdistance(gem_mollers[d][i], geo.Ebeam));
            if (i >= 1) {
                auto c = physics.GetMollerCenter(gem_mollers[d][i-1], gem_mollers[d][i]);
                if(c[0] == c[1] && c[0] == 0.f) continue; // skip zero-center events (likely bad recon)
                center_gem[d]->Fill(c[0], c[1]);
                center_gem_x[d]->Fill(c[0]);
                center_gem_y[d]->Fill(c[1]);
            }
            hits_gem[d]->Fill(gem_mollers[d][i].first.x,  gem_mollers[d][i].first.y);
            hits_gem[d]->Fill(gem_mollers[d][i].second.x, gem_mollers[d][i].second.y);
        }
    }

    // --- fit and print results ---
    std::string plot_dir = "Poscalib_result/" + run_str;
    float hycal_vertex_z = fitAndDraw(vertex_hycal,   plot_dir + "/hycal_vertex_z",  geo.hycal_z);
    float hycal_center_x = fitAndDraw(center_hycal_x, plot_dir + "/hycal_center_x",  geo.hycal_x + HyCal_shift_x);
    float hycal_center_y = fitAndDraw(center_hycal_y, plot_dir + "/hycal_center_y",  geo.hycal_y + HyCal_shift_y);

    float gem_vertex_z[4], gem_center_x[4], gem_center_y[4];
    for (int d = 0; d < 4; d++) {
        gem_vertex_z[d] = fitAndDraw(vertex_gem[d],   plot_dir + "/gem" + std::to_string(d) + "_vertex_z",  geo.gem_z[d]);
        gem_center_x[d] = fitAndDraw(center_gem_x[d], plot_dir + "/gem" + std::to_string(d) + "_center_x",  geo.gem_x[d] + GEM_shift_x[d]);
        gem_center_y[d] = fitAndDraw(center_gem_y[d], plot_dir + "/gem" + std::to_string(d) + "_center_y",  geo.gem_y[d] + GEM_shift_y[d]);
    }

    std::cerr << "HyCal vertex z: " << hycal_vertex_z << " mm  (survey: " << geo.hycal_z << " mm)\n";
    std::cerr << "HyCal center x: " << hycal_center_x << " mm  (survey: " << geo.hycal_x << " mm)\n";
    std::cerr << "HyCal center y: " << hycal_center_y << " mm  (survey: " << geo.hycal_y << " mm)\n";
    for (int d = 0; d < 4; d++) {
        std::cerr << "GEM " << d << " vertex z: " << gem_vertex_z[d] << " mm  (survey: " << geo.gem_z[d] << " mm)\n";
        std::cerr << "GEM " << d << " center x: " << gem_center_x[d] << " mm  (survey: " << geo.gem_x[d] << " mm)\n";
        std::cerr << "GEM " << d << " center y: " << gem_center_y[d] << " mm  (survey: " << geo.gem_y[d] << " mm)\n";
    }

    float target_delta_z[4];
    for (int d = 0; d < 4; d++) {
        target_delta_z[d] = geo.gem_z[d] - gem_vertex_z[d];
        std::cerr << "target Z position from GEM\n";
        std::cerr << "  GEM " << d << ": " << target_delta_z[d] << " mm\n";
    }

    // --- summary in HyCal-center frame ---
    // geo.* positions are in beam-center frame (LoadRunConfig subtracted target).
    // HyCal is the reference: its x/y are 0 by definition in this frame.
    // Beam center in HyCal frame = target + Moller-center correction.
    // GEM positions in HyCal frame = beam-frame position + target offset.
    const float tx = geo.target_x, ty = geo.target_y, tz = geo.target_z;
    std::cerr << "\n";
    std::cerr << "========================================\n";
    std::cerr << " Summary in HyCal-center frame (mm)\n";
    std::cerr << "========================================\n";
    std::cerr << " HyCal z (calibrated):  " << (hycal_vertex_z + tz) << "\n";
    std::cerr << "\n";
    std::cerr << " Beam center (calibrated):\n";
    std::cerr << "   \"target\": [" << (tx + hycal_center_x)
              << ", " << (ty + hycal_center_y)
              << ", " << (target_delta_z[0]+target_delta_z[1]+target_delta_z[2]+target_delta_z[3]) / 4.0
              << "]\n";
    std::cerr << "\n";
    std::cerr << " GEM positions (calibrated):\n";
    for (int d = 0; d < 4; d++) {
        std::cerr << "   {\"id\": " << d
                  << ", \"position\": [" << (geo.gem_x[d] - gem_center_x[d] + tx)
                  << ", " << (geo.gem_y[d] - gem_center_y[d] + ty)
                  << ", " << (gem_vertex_z[d] + tz) << "]"
                  << ", \"tilting\": [" << geo.gem_tilt_x[d]
                  << ", " << geo.gem_tilt_y[d]
                  << ", " << geo.gem_tilt_z[d] << "]}\n";
    }
    std::cerr << "========================================\n";

    // --- save histograms ---
    TFile outfile(output.c_str(), "RECREATE");
    E_theta->Write();
    vertex_hycal->Write();   center_hycal->Write();
    center_hycal_x->Write(); center_hycal_y->Write();
    hits_hycal->Write();
    for (int d = 0; d < 4; d++) {
        vertex_gem[d]->Write();   center_gem[d]->Write();
        center_gem_x[d]->Write(); center_gem_y[d]->Write();
        hits_gem[d]->Write();
    }
    outfile.Close();
    std::cerr << "Calibration histograms saved to " << output << "\n";
}

// ── fitAndDraw ────────────────────────────────────────────────────────────
// Fit range is determined automatically from the half-maximum (50%) points:
// scan outward from the peak bin to find the FWHM, then fit ± 1.5 * HWHM.
static float fitAndDraw(TH1F *hist, const std::string &out_path,
                        float survey_position)
{
    TCanvas c("", "", 800, 600);

    int    peak_bin = hist->GetMaximumBin();
    double peak_x   = hist->GetBinCenter(peak_bin);
    double peak_val = hist->GetMaximum();
    double half_max = peak_val * 0.5;
    int    nb       = hist->GetNbinsX();

    // scan left of peak for half-maximum crossing
    double lo_x = hist->GetBinCenter(1);
    for (int b = peak_bin - 1; b >= 1; --b) {
        if (hist->GetBinContent(b) <= half_max) {
            lo_x = hist->GetBinCenter(b);
            break;
        }
    }
    // scan right of peak for half-maximum crossing
    double hi_x = hist->GetBinCenter(nb);
    for (int b = peak_bin + 1; b <= nb; ++b) {
        if (hist->GetBinContent(b) <= half_max) {
            hi_x = hist->GetBinCenter(b);
            break;
        }
    }

    // half-width at half-maximum; guard against degenerate (flat/empty) histograms
    double hwhm = 0.5 * (hi_x - lo_x);
    if (hwhm < hist->GetBinWidth(1)) hwhm = hist->GetBinWidth(1);

    double fit_lo = peak_x - 1. * hwhm;
    double fit_hi = peak_x + 1. * hwhm;

    hist->Fit("gaus", "rq", "", fit_lo, fit_hi);
    hist->Draw();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04f);
    TF1 *gf = hist->GetFunction("gaus");
    if (gf) {
        latex.DrawLatex(0.15, 0.85, Form("%.2f mm +- %.2f mm",
                                          gf->GetParameter(1), gf->GetParError(1)));
        latex.DrawLatex(0.15, 0.80, Form("Survey position: %.2f mm", survey_position));
        latex.DrawLatex(0.15, 0.75, Form("FWHM: %.2f mm", 2.0 * hwhm));
    }
    fs::create_directories(fs::path(out_path).parent_path());
    c.SaveAs((out_path + ".png").c_str());
    return gf ? gf->GetParameter(1) : static_cast<float>(peak_x);
}