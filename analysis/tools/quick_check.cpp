// quick_check.C — ROOT script version 
//
// Reads reconstructed ROOT tree (output of replay_recon), runs physics
// analysis using PhysicsTools and MatchingTools from prad2det, and saves
// histograms to an output ROOT file.
// Usage:
//   quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]
//   -o  output ROOT file (default: input filename with _quick_check.root suffix)
//   -n  max events to process (default: all)
// Example:
//   quick_check recon.root -o recon_check.root -n 10000
//   quick_check recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TCanvas.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Forward declaration
float fitAndDraw(TH1F* hist, const std::string& out_path, float fit_range);

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 2108.f;
    int run_id = 12345;

    // --- geometry constants (can be made configurable) ---
    const float gem_z[4] = {5407.f + 39.71f/2, 5407.f - 39.71f/2,
                            5807.f + 39.71f/2, 5807.f - 39.71f/2};
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }
    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- load run config: assign run_id and Ebeam from gRunConfig ---
    run_id = analysis::get_run_int(root_files[0]);
    gRunConfig = analysis::LoadRunConfig(dbDir + "/runinfo/general.json", run_id);
    Ebeam = gRunConfig.Ebeam > 0.f ? gRunConfig.Ebeam : Ebeam;

    std::cout << "Processing run " << run_id << " with Ebeam = " << Ebeam << " MeV\n";

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);
    MatchingTools matching;
    
    // --- setup TChain and branches ---
    TChain *chain = new TChain("recon");
    for (const auto &f : root_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);

    // --- output file ---
    TString outName = output;
    if (outName.IsNull()) {
        outName = root_files[0];
        outName.ReplaceAll("_recon.root", "_quick_check.root");
    }
    TFile outfile(outName, "RECREATE");

    // --- histograms ---
    TH2F *hit_pos = new TH2F("hit_pos",
        "Hit positions;X (mm);Y (mm)", 250, -500, 500, 250, -500, 500);
    TH1F *h_1cl = new TH1F("one_cluster_energy",
        "Single-cluster energy;E (MeV);Counts", 1000, 0, 4000);
    TH1F *h_2cl = new TH1F("two_cluster_energy",
        "Two-cluster energy;E (MeV);Counts", 1000, 0, 4000);
    TH1F *h_all = new TH1F("clusters_energy",
        "All clusters;E (MeV);Counts", 1000, 0, 4000);
    TH1F *h_tot = new TH1F("total_energy",
        "Total energy per event;E (MeV);Counts", 1000, 0, 4000);
    TH2F *h2_energy_theta_ep_ee = new TH2F("energy_vs_theta",
        "Energy vs Theta(1 cluster);Theta (deg);Energy (MeV)", 160, 0, 8, 4000, 0, 4000);

    MollerData hycal_mollers;

    // --- event loop : basic distributions + Moller candidates on HyCal ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;

    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 100000 == 0)
            std::cerr << "\rPass 1: " << i << " / " << N << std::flush;

        for (int j = 0; j < ev.n_clusters; j++) {
            float r = std::sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]);
            float theta = std::atan(r / ev.cl_z[j]) * 180.f / M_PI;

            physics.FillEnergyVsModule(ev.cl_center[j], ev.cl_energy[j]);
            physics.FillEnergyVsTheta(theta, ev.cl_energy[j]);
            hit_pos->Fill(ev.cl_x[j], ev.cl_y[j]);
            h_all->Fill(ev.cl_energy[j]);

            if (ev.cl_nblocks[j] > 2) {
                int mod_id = ev.cl_center[j];
                auto mod = hycal.module_by_id(mod_id);
                if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
                // require hit to be in central 3x3 of a 5x5 grid (|xd|,|yd| < 0.3)
                float xd = (ev.cl_x[j] - (float)mod->x) / (float)mod->size_x;
                float yd = (ev.cl_y[j] - (float)mod->y) / (float)mod->size_y;
                if (std::abs(xd) >= 0.3f || std::abs(yd) >= 0.3f) continue;
                physics.FillEnergyVsTheta(theta, ev.cl_energy[j]);
            }
        }
        h_tot->Fill(ev.total_energy);

        if (ev.n_clusters == 1) {
            physics.FillModuleEnergy(ev.cl_center[0], ev.cl_energy[0]);
            h_1cl->Fill(ev.cl_energy[0]);
            h2_energy_theta_ep_ee->Fill(std::atan(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180.f / M_PI,
                                      ev.cl_energy[0]);
        }

        if (ev.n_clusters == 2) {
            h_2cl->Fill(ev.cl_energy[0]);
            h_2cl->Fill(ev.cl_energy[1]);

            //h2_energy_theta_ep_ee->Fill(std::atan(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180.f / M_PI,
             //                         ev.cl_energy[0]);
            //h2_energy_theta_ep_ee->Fill(std::atan(std::sqrt(ev.cl_x[1]*ev.cl_x[1] + ev.cl_y[1]*ev.cl_y[1]) / ev.cl_z[1]) * 180.f / M_PI,
              //                        ev.cl_energy[1]);

            float Epair = ev.cl_energy[0] + ev.cl_energy[1];
            float sigma = Ebeam * 0.025f / std::sqrt(Ebeam / 1000.f);
            if (std::abs(Epair - Ebeam) < 4. * sigma) {
                MollerEvent mp(
                    {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                    {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
                physics.FillMollerPhiDiff(physics.GetMollerPhiDiff(mp));
                if(fabs(physics.GetMollerPhiDiff(mp)) > 5.f ) continue;
                hycal_mollers.push_back(mp);
                physics.Fill2armMollerPosHist(mp.first.x, mp.first.y);
                physics.Fill2armMollerPosHist(mp.second.x, mp.second.y);
            }
        }
    }

    // --- Moller vertex analysis ---
    for (size_t i = 0; i < hycal_mollers.size(); i++) {
        physics.FillMollerZ(physics.GetMollerZdistance(hycal_mollers[i], Ebeam));
        if (i >= 1) {
            auto c = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            physics.FillMollerXY(c[0], c[1]);
        }
    }

    float hycal_vertex_z = fitAndDraw(physics.GetMollerZHist(), "Poscalib_result/hycal_vertex_z", 100.);
    float hycal_center_x = fitAndDraw(physics.GetMollerXHist(), "Poscalib_result/hycal_center_x", 2.);
    float hycal_center_y = fitAndDraw(physics.GetMollerYHist(), "Poscalib_result/hycal_center_y", 2.);

    // --- write output ---
    outfile.cd();
    hit_pos->Write();

    outfile.mkdir("energy_plots"); outfile.cd("energy_plots");
    if (physics.GetEnergyVsModuleHist()) physics.GetEnergyVsModuleHist()->Write();
    if (physics.GetEnergyVsThetaHist())  physics.GetEnergyVsThetaHist()->Write();
    h_1cl->Write(); h_2cl->Write(); h_all->Write(); h_tot->Write();
    h2_energy_theta_ep_ee->Write();

    outfile.cd();
    outfile.mkdir("physics_yields"); outfile.cd("physics_yields");
    auto h_ep = physics.GetEpYieldHist(physics.GetEnergyVsThetaHist(), Ebeam);
    auto h_ee = physics.GetEeYieldHist(physics.GetEnergyVsThetaHist(), Ebeam);
    auto h_ratio = physics.GetYieldRatioHist(h_ep.get(), h_ee.get());
    if (h_ep) h_ep->Write();
    if (h_ee) h_ee->Write();
    if (h_ratio) h_ratio->Write();

    outfile.cd();
    outfile.mkdir("moller_analysis"); outfile.cd("moller_analysis");
    if (physics.Get2armMollerPosHist()) physics.Get2armMollerPosHist()->Write();
    if (physics.GetMollerPhiDiffHist()) physics.GetMollerPhiDiffHist()->Write();
    if (physics.GetMollerXHist()) physics.GetMollerXHist()->Write();
    if (physics.GetMollerYHist()) physics.GetMollerYHist()->Write();
    if (physics.GetMollerZHist()) physics.GetMollerZHist()->Write();

    physics.FillNeventsModuleMap();
    TH2F *h_map = physics.GetNeventsModuleMapHist();
    TCanvas *c_map = new TCanvas("c_map", "Number of Events per Module", 1200, 1200);
    h_map->Draw("COLZ");
    c_map->Write();
    h_map->Write();

    outfile.mkdir("module_energy"); outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); i++) {
        TH1F *h = physics.GetModuleEnergyHist(i);
        if (h && h->GetEntries() > 0) h->Write();
    }

    outfile.Close();
    
    // physics.Resolution2Database(run_id); // example run ID

    std::cerr << "Result saved -> " << outName.Data() << "\n";
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_recon.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

float fitAndDraw(TH1F* hist, const std::string& out_path, const float fit_range){
    TCanvas *c = new TCanvas("", "", 800, 600);
    float mean = hist->GetBinCenter(hist->GetMaximumBin());
    hist->Fit("gaus", "rq", "", mean-fit_range, mean+fit_range);
    hist->Draw();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.15, 0.85, Form("%.2f mm +- %.2f mm", hist->GetFunction("gaus")->GetParameter(1), hist->GetFunction("gaus")->GetParError(1)));
    fs::create_directories(fs::path(out_path).parent_path());
    c->SaveAs((out_path + ".png").c_str());
    delete c;

    return hist->GetFunction("gaus")->GetParameter(1);
}