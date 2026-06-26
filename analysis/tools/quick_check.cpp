// quick_check.C — ROOT script version 
//
// Reads reconstructed ROOT tree (output of replay_recon), runs physics
// analysis using PhysicsTools and MatchingTools from prad2det, and saves
// histograms to an output ROOT file.
// Usage:
//   quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]
//   -o  output ROOT file (default: input filename with _quick_check.root suffix)
//   -n  max events to process (default: all)
//   -j  number of input-file worker threads (default: 4)
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
#include <TROOT.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
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
static std::string makeDefaultOutput(const std::string &input_path);

bool inHyCal(double xmm, double ymm) {
    const double module = 20.75; // mm
    return (fabs(xmm) > module * 2.5 || fabs(ymm) > module * 2.5)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
}

struct QuickResult {
    std::unique_ptr<PhysicsTools> physics;
    std::unique_ptr<TH2F> hit_pos;
    std::unique_ptr<TH1F> h_1cl;
    std::unique_ptr<TH1F> h_2cl;
    std::unique_ptr<TH1F> h_all;
    std::unique_ptr<TH1F> h_tot;
    std::unique_ptr<TH2F> h2_energy_theta_ep_ee;
    MollerData mollers;
    Long64_t processed = 0;
};

static void detach(TH1 *h)
{
    if (h) h->SetDirectory(nullptr);
}

static std::unique_ptr<QuickResult> makeResult(fdec::HyCalSystem &hycal)
{
    auto r = std::make_unique<QuickResult>();
    r->physics = std::make_unique<PhysicsTools>(hycal);
    r->hit_pos = std::make_unique<TH2F>("hit_pos",
        "Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h_1cl = std::make_unique<TH1F>("one_cluster_energy",
        "Single-cluster energy;E (MeV);Counts", 4000, 0, 4000);
    r->h_2cl = std::make_unique<TH1F>("two_cluster_energy",
        "Two-cluster energy;E (MeV);Counts", 4000, 0, 4000);
    r->h_all = std::make_unique<TH1F>("clusters_energy",
        "All clusters;E (MeV);Counts", 4000, 0, 4000);
    r->h_tot = std::make_unique<TH1F>("total_energy",
        "Total energy per event;E (MeV);Counts", 4000, 0, 4000);
    r->h2_energy_theta_ep_ee = std::make_unique<TH2F>("energy_vs_theta",
        "Energy vs Theta(1 cluster);Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    detach(r->hit_pos.get());
    detach(r->h_1cl.get());
    detach(r->h_2cl.get());
    detach(r->h_all.get());
    detach(r->h_tot.get());
    detach(r->h2_energy_theta_ep_ee.get());
    return r;
}

static Long64_t reconEntries(const std::string &path)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return 0;
    TTree *t = dynamic_cast<TTree *>(f->Get("recon"));
    return t ? t->GetEntries() : 0;
}

static bool processFile(const std::string &path,
                        Long64_t max_entries,
                        fdec::HyCalSystem &hycal,
                        float Ebeam,
                        QuickResult &out)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << path << "\n";
        return false;
    }
    TTree *tree = dynamic_cast<TTree *>(f->Get("recon"));
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in " << path << "\n";
        return false;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);
    Long64_t n = tree->GetEntries();
    if (max_entries >= 0 && max_entries < n) n = max_entries;

    auto &physics = *out.physics;
    for (Long64_t i = 0; i < n; i++) {
        tree->GetEntry(i);
        for (int j = 0; j < ev.n_clusters; j++) {
            float r = std::sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]);
            float theta = std::atan(r / ev.cl_z[j]) * 180.f / M_PI;

            physics.FillEnergyVsModule(ev.cl_center[j], ev.cl_energy[j]);
            out.hit_pos->Fill(ev.cl_x[j], ev.cl_y[j]);
            out.h_all->Fill(ev.cl_energy[j]);

            if (ev.cl_nblocks[j] > 3 && inHyCal(ev.cl_x[j], ev.cl_y[j])) {
                physics.FillEnergyVsTheta(theta, ev.cl_energy[j]);
            }
        }
        out.h_tot->Fill(ev.total_energy);

        if (ev.n_clusters == 1 && inHyCal(ev.cl_x[0], ev.cl_y[0])) {
            physics.FillModuleEnergy(ev.cl_center[0], ev.cl_energy[0]);
            out.h_1cl->Fill(ev.cl_energy[0]);
            out.h2_energy_theta_ep_ee->Fill(
                std::atan(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180.f / M_PI,
                ev.cl_energy[0]);
        }

        if (ev.n_clusters == 2 && inHyCal(ev.cl_x[0], ev.cl_y[0]) && inHyCal(ev.cl_x[1], ev.cl_y[1])) {
            out.h_2cl->Fill(ev.cl_energy[0]);
            out.h_2cl->Fill(ev.cl_energy[1]);

            float Epair = ev.cl_energy[0] + ev.cl_energy[1];
            float sigma = Ebeam * 0.025f / std::sqrt(Ebeam / 1000.f);
            if (std::abs(Epair - Ebeam) < 4. * sigma) {
                MollerEvent mp(
                    {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                    {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
                physics.FillMollerPhiDiff(physics.GetMollerPhiDiff(mp));
                if (fabs(physics.GetMollerPhiDiff(mp)) > 5.f) continue;
                out.mollers.push_back(mp);
                physics.Fill2armMollerPosHist(mp.first.x, mp.first.y);
                physics.Fill2armMollerPosHist(mp.second.x, mp.second.y);
            }
        }
    }
    out.processed += n;
    return true;
}

static void mergeResult(QuickResult &dst, const QuickResult &src, fdec::HyCalSystem &hycal)
{
    dst.hit_pos->Add(src.hit_pos.get());
    dst.h_1cl->Add(src.h_1cl.get());
    dst.h_2cl->Add(src.h_2cl.get());
    dst.h_all->Add(src.h_all.get());
    dst.h_tot->Add(src.h_tot.get());
    dst.h2_energy_theta_ep_ee->Add(src.h2_energy_theta_ep_ee.get());

    dst.physics->GetEnergyVsModuleHist()->Add(src.physics->GetEnergyVsModuleHist());
    dst.physics->GetEnergyVsThetaHist()->Add(src.physics->GetEnergyVsThetaHist());
    dst.physics->Get2armMollerPosHist()->Add(src.physics->Get2armMollerPosHist());
    dst.physics->GetMollerPhiDiffHist()->Add(src.physics->GetMollerPhiDiffHist());
    for (int i = 0; i < hycal.module_count(); ++i) {
        int module_id = hycal.module(i).id;
        TH1F *d = dst.physics->GetModuleEnergyHist(module_id);
        TH1F *s = src.physics->GetModuleEnergyHist(module_id);
        if (d && s) d->Add(s);
    }
    dst.processed += src.processed;
}

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 2108.f;
    int run_id = 12345;
    
    int max_events = -1;
    int num_threads = 4;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:j:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'j': num_threads = std::atoi(optarg); break;
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
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]\n";
        return 1;
    }
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(root_files.size())));
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);

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
    std::cout << "Processing " << root_files.size() << " file(s) with "
              << num_threads << " thread(s)\n";

    std::vector<Long64_t> file_limits(root_files.size(), -1);
    if (max_events > 0) {
        Long64_t remaining = max_events;
        for (size_t i = 0; i < root_files.size(); ++i) {
            Long64_t n = reconEntries(root_files[i]);
            file_limits[i] = std::min(n, remaining);
            remaining -= file_limits[i];
            if (remaining <= 0) {
                for (size_t j = i + 1; j < root_files.size(); ++j)
                    file_limits[j] = 0;
                break;
            }
        }
    }

    std::vector<std::unique_ptr<QuickResult>> results(root_files.size());
    std::atomic<size_t> next_file{0};
    std::atomic<int> errors{0};
    std::mutex io_mtx;
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&]() {
            while (true) {
                size_t idx = next_file.fetch_add(1);
                if (idx >= root_files.size()) break;
                auto res = makeResult(hycal);
                {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "Processing file [" << (idx + 1) << "/"
                              << root_files.size() << "]: " << root_files[idx] << "\n";
                }
                if (!processFile(root_files[idx], file_limits[idx], hycal, Ebeam, *res))
                    ++errors;
                results[idx] = std::move(res);
            }
        });
    }
    for (auto &t : threads) t.join();
    if (errors > 0) return 1;

    auto merged = makeResult(hycal);
    MollerData hycal_mollers;
    for (auto &res : results) {
        if (!res) continue;
        mergeResult(*merged, *res, hycal);
        hycal_mollers.insert(hycal_mollers.end(), res->mollers.begin(), res->mollers.end());
    }
    PhysicsTools &physics = *merged->physics;
    TH2F *hit_pos = merged->hit_pos.get();
    TH1F *h_1cl = merged->h_1cl.get();
    TH1F *h_2cl = merged->h_2cl.get();
    TH1F *h_all = merged->h_all.get();
    TH1F *h_tot = merged->h_tot.get();
    TH2F *h2_energy_theta_ep_ee = merged->h2_energy_theta_ep_ee.get();

    TString outName = output;
    if (outName.IsNull())
        outName = makeDefaultOutput(root_files[0]);
    TFile outfile(outName, "RECREATE");

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
            std::string name = entry.path().filename().string();
            if (entry.is_regular_file() &&
                name.find("_recon") != std::string::npos &&
                name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

static std::string makeDefaultOutput(const std::string &input_path)
{
    fs::path p(input_path);
    std::string name = p.filename().string();
    const std::string ext = ".root";
    if (name.size() >= ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0) {
        name.insert(name.size() - ext.size(), "_quick_check");
    } else {
        name += "_quick_check.root";
    }
    return (p.parent_path() / name).string();
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
