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
#include <limits>
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

const int Nbins = 33;
const float binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

struct QuickResult {
    std::unique_ptr<PhysicsTools> physics;
    std::unique_ptr<TH2F> hit_pos;
    std::unique_ptr<TH1F> h_1cl;
    std::unique_ptr<TH1F> h_2cl;
    std::unique_ptr<TH1F> h_all;
    std::unique_ptr<TH1F> h_tot;
    std::unique_ptr<TH2F> h2_energy_theta_ep_ee;
    MollerData mollers;
    MollerData mollers_hc;
    Long64_t processed = 0;

    std::unique_ptr<TH2F> h2_ep_hits;
    std::unique_ptr<TH2F> h2_ee_hits;
    std::unique_ptr<TH2F> h2_ep_E_angle;
    std::unique_ptr<TH2F> h2_ee_E_angle;

    std::unique_ptr<TH1F> h_ep_yield;
    std::unique_ptr<TH1F> h_ee_yield;
    std::unique_ptr<TH1F> h_ep_ee_ratio;

    std::unique_ptr<TH1F> h_ee_center_x;
    std::unique_ptr<TH1F> h_ee_center_y;
    std::unique_ptr<TH1F> h_ee_vertex_z;

    std::unique_ptr<TH2F> h2_ep_hits_hc;
    std::unique_ptr<TH2F> h2_ee_hits_hc;
    std::unique_ptr<TH2F> h2_ep_E_angle_hc;
    std::unique_ptr<TH2F> h2_ee_E_angle_hc;

    std::unique_ptr<TH1F> h_ee_center_x_hc;
    std::unique_ptr<TH1F> h_ee_center_y_hc;
    std::unique_ptr<TH1F> h_ee_vertex_z_hc;
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
    
    r->h2_ep_hits = std::make_unique<TH2F>("ep_hits",
        "EP Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ee_hits = std::make_unique<TH2F>("ee_hits",
        "EE Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ep_E_angle = std::make_unique<TH2F>("ep_E_angle",
        "EP Energy vs Angle;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    r->h2_ee_E_angle = std::make_unique<TH2F>("ee_E_angle",
        "EE Energy vs Angle;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h_ep_yield = std::make_unique<TH1F>("ep_yield",
        "EP Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ee_yield = std::make_unique<TH1F>("ee_yield",
        "EE Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ep_ee_ratio = std::make_unique<TH1F>("ep_ee_ratio",
        "EP/EE Yield Ratio;Scattering Angle (deg);Counts", Nbins, binEdge);

    r->h_ee_center_x = std::make_unique<TH1F>("ee_center_x",
        "EE Center X;X (mm);Counts", 800, -20, 20);
    r->h_ee_center_y = std::make_unique<TH1F>("ee_center_y",
        "EE Center Y;Y (mm);Counts", 800, -20, 20);
    r->h_ee_vertex_z = std::make_unique<TH1F>("ee_vertex_z",
        "EE Vertex Z;Z (mm);Counts", 8000, 5000, 9000);

    r->h2_ep_hits_hc = std::make_unique<TH2F>("ep_hits_hc",
        "EP Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ee_hits_hc = std::make_unique<TH2F>("ee_hits_hc",
        "EE Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ep_E_angle_hc = std::make_unique<TH2F>("ep_E_angle_hc",
        "EP Energy vs Angle hycal;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    r->h2_ee_E_angle_hc = std::make_unique<TH2F>("ee_E_angle_hc",
        "EE Energy vs Angle hycal;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h_ee_center_x_hc = std::make_unique<TH1F>("ee_center_x_hc",
        "EE Center X hycal;X (mm);Counts", 800, -20, 20);
    r->h_ee_center_y_hc = std::make_unique<TH1F>("ee_center_y_hc",
        "EE Center Y hycal;Y (mm);Counts", 800, -20, 20);
    r->h_ee_vertex_z_hc = std::make_unique<TH1F>("ee_vertex_z_hc",
        "EE Vertex Z hycal;Z (mm);Counts", 8000, 5000, 9000);


    detach(r->hit_pos.get());
    detach(r->h_1cl.get());
    detach(r->h_2cl.get());
    detach(r->h_all.get());
    detach(r->h_tot.get());
    detach(r->h2_energy_theta_ep_ee.get());
    detach(r->h2_ep_hits.get());
    detach(r->h2_ee_hits.get());
    detach(r->h2_ep_E_angle.get());
    detach(r->h2_ee_E_angle.get());
    detach(r->h_ep_yield.get());
    detach(r->h_ee_yield.get());
    detach(r->h_ep_ee_ratio.get());
    detach(r->h_ee_center_x.get());
    detach(r->h_ee_center_y.get());
    detach(r->h_ee_vertex_z.get());
    detach(r->h2_ep_hits_hc.get());
    detach(r->h2_ee_hits_hc.get());
    detach(r->h2_ep_E_angle_hc.get());
    detach(r->h2_ee_E_angle_hc.get());
    detach(r->h_ee_center_x_hc.get());
    detach(r->h_ee_center_y_hc.get());
    detach(r->h_ee_vertex_z_hc.get());
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
            out.h2_ep_hits_hc->Fill(ev.cl_x[0], ev.cl_y[0]);
            out.h2_ep_E_angle_hc->Fill(
                std::atan(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180.f / M_PI,
                ev.cl_energy[0]);
        }

        if (ev.n_clusters == 2 && inHyCal(ev.cl_x[0], ev.cl_y[0]) && inHyCal(ev.cl_x[1], ev.cl_y[1])) {
            out.h_2cl->Fill(ev.cl_energy[0]);
            out.h_2cl->Fill(ev.cl_energy[1]);

            float Epair = ev.cl_energy[0] + ev.cl_energy[1];
            float sigma = Ebeam * 0.035f / std::sqrt(Ebeam / 1000.f);
            if (std::abs(Epair - Ebeam) < 4. * sigma) {
                MollerEvent mp(
                    {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                    {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
                physics.FillMollerPhiDiff(physics.GetMollerPhiDiff(mp));
                if(physics.GetMollerPhiDiff(mp) < 10.f) {
                    out.h2_ee_hits_hc->Fill(ev.cl_x[0], ev.cl_y[0]);
                    out.h2_ee_hits_hc->Fill(ev.cl_x[1], ev.cl_y[1]);
                    float t1 = std::atan2(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]), ev.cl_z[0]) * 180.f / M_PI;
                    float t2 = std::atan2(std::sqrt(ev.cl_x[1]*ev.cl_x[1] + ev.cl_y[1]*ev.cl_y[1]), ev.cl_z[1]) * 180.f / M_PI;
                    out.h2_ee_E_angle_hc->Fill(t1, ev.cl_energy[0]);
                    out.h2_ee_E_angle_hc->Fill(t2, ev.cl_energy[1]);
                    out.mollers_hc.push_back(mp);

                    if (out.mollers_hc.size() > 1) {
                        auto center = physics.GetMollerCenter(out.mollers_hc[out.mollers_hc.size() - 2], mp);
                        out.h_ee_center_x_hc->Fill(center[0]);
                        out.h_ee_center_y_hc->Fill(center[1]);
                        if (out.mollers_hc.size() > 2) {
                            auto center2 = physics.GetMollerCenter(out.mollers_hc[out.mollers_hc.size() - 3], mp);
                            out.h_ee_center_x_hc->Fill(center2[0]);
                            out.h_ee_center_y_hc->Fill(center2[1]);
                        }
                    }
                    float vertex = physics.GetMollerZdistance(mp, Ebeam);
                    out.h_ee_vertex_z_hc->Fill(vertex);

                }
            }
        }

        //loop over GEM matched results
        std::vector<analysis::HCHit> moller_Hits_candidate;
        for (int j = 0; j < ev.matchNum; j++) {
            float x = ev.mHit_gx[j][1];
            float y = ev.mHit_gy[j][1];
            float z = ev.mHit_gz[j][1];
            float E = ev.mHit_E[j];
            //project to HyCal plane
            float scale = ev.mHit_z[j] / z;
            x *= scale;
            y *= scale;
            z *= scale;
            float theta = std::atan(std::sqrt(x*x + y*y) / z) * 180.f / M_PI;

            if(!inHyCal(x, y)) continue;

            float expectE = physics.ExpectedEnergy(theta, Ebeam, "ep");
            if (fabs(E - expectE) < 3.f * expectE * 0.035f / std::sqrt(E/1000.f)) {
                out.h2_ep_hits->Fill(x, y);
                out.h2_ep_E_angle->Fill(theta, E);
                out.h_ep_yield->Fill(theta);
            }

            if (E > 60. && E < Ebeam - 0.035 * Ebeam / sqrt(Ebeam/1000.))
                moller_Hits_candidate.push_back({x, y, z, E});
        }

        std::sort(moller_Hits_candidate.begin(), moller_Hits_candidate.end(),
                  [](const analysis::HCHit &a, const analysis::HCHit &b){ return a.energy > b.energy; });

        analysis::MollerData mollerData_event;
        int nCand = moller_Hits_candidate.size();

        for (int j = 0; j < nCand; j++) {
            const auto &hit1 = moller_Hits_candidate[j];
            float theta1 = std::atan2(std::sqrt(hit1.x*hit1.x + hit1.y*hit1.y), hit1.z)
                         * 180.f / M_PI;
            for (int k = j + 1; k < nCand; k++) {
                const auto &hit2 = moller_Hits_candidate[k];
                float theta2 = std::atan2(std::sqrt(hit2.x*hit2.x + hit2.y*hit2.y), hit2.z)
                             * 180.f / M_PI;
                MollerEvent mev(
                    {hit1.x, hit1.y, hit1.z, hit1.energy},
                    {hit2.x, hit2.y, hit2.z, hit2.energy});
                if (!physics.isMoller_kinematic(theta1, hit1.energy,
                                                theta2, hit2.energy,
                                                Ebeam, 0.035f))
                    continue;
                if (fabs(physics.GetMollerPhiDiff(mev)) > 10.f) continue;
                mollerData_event.push_back(mev);
            }
        }

        if (mollerData_event.empty()) continue;

        if (mollerData_event.size() > 1) {
            auto getPt = [](const MollerEvent &mev) -> float {
                float sin_t1 = std::sqrt(mev.first.x*mev.first.x + mev.first.y*mev.first.y)
                             / std::sqrt(mev.first.z*mev.first.z + mev.first.x*mev.first.x + mev.first.y*mev.first.y);
                float sin_t2 = std::sqrt(mev.second.x*mev.second.x + mev.second.y*mev.second.y)
                             / std::sqrt(mev.second.z*mev.second.z + mev.second.x*mev.second.x + mev.second.y*mev.second.y);
                return std::fabs(mev.first.E * sin_t1 - mev.second.E * sin_t2);
            };
            auto best = std::min_element(
                mollerData_event.begin(), mollerData_event.end(),
                [&](const MollerEvent &a, const MollerEvent &b) {
                    return getPt(a) < getPt(b);
                });
            MollerEvent bestPair = *best;
            mollerData_event.clear();
            mollerData_event.push_back(bestPair);
        }

        MollerEvent &mev = mollerData_event.front();
        out.mollers.push_back(mev);

        float t1 = std::atan2(std::sqrt(mev.first.x*mev.first.x + mev.first.y*mev.first.y),
                              mev.first.z) * 180.f / M_PI;
        float t2 = std::atan2(std::sqrt(mev.second.x*mev.second.x + mev.second.y*mev.second.y),
                              mev.second.z) * 180.f / M_PI;
        out.h2_ee_hits->Fill(mev.first.x, mev.first.y);
        out.h2_ee_hits->Fill(mev.second.x, mev.second.y);
        out.h2_ee_E_angle->Fill(t1, mev.first.E);
        out.h2_ee_E_angle->Fill(t2, mev.second.E);
        out.h_ee_yield->Fill(t1);
        out.h_ee_yield->Fill(t2);

        float vertex = physics.GetMollerZdistance(mev, Ebeam);
        out.h_ee_vertex_z->Fill(vertex);

        if (out.mollers.size() > 1) {
            auto center = physics.GetMollerCenter(out.mollers[out.mollers.size() - 2], mev);
            out.h_ee_center_x->Fill(center[0]);
            out.h_ee_center_y->Fill(center[1]);
            if (out.mollers.size() > 2) {
                auto center2 = physics.GetMollerCenter(out.mollers[out.mollers.size() - 3], mev);
                out.h_ee_center_x->Fill(center2[0]);
                out.h_ee_center_y->Fill(center2[1]);
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
    dst.h2_ep_hits->Add(src.h2_ep_hits.get());
    dst.h2_ee_hits->Add(src.h2_ee_hits.get());
    dst.h2_ep_E_angle->Add(src.h2_ep_E_angle.get());
    dst.h2_ee_E_angle->Add(src.h2_ee_E_angle.get());
    dst.h_ep_yield->Add(src.h_ep_yield.get());
    dst.h_ee_yield->Add(src.h_ee_yield.get());
    dst.h_ee_center_x->Add(src.h_ee_center_x.get());
    dst.h_ee_center_y->Add(src.h_ee_center_y.get());
    dst.h_ee_vertex_z->Add(src.h_ee_vertex_z.get());
    dst.h2_ep_hits_hc->Add(src.h2_ep_hits_hc.get());
    dst.h2_ee_hits_hc->Add(src.h2_ee_hits_hc.get());
    dst.h2_ep_E_angle_hc->Add(src.h2_ep_E_angle_hc.get());
    dst.h2_ee_E_angle_hc->Add(src.h2_ee_E_angle_hc.get());
    dst.h_ee_center_x_hc->Add(src.h_ee_center_x_hc.get());
    dst.h_ee_center_y_hc->Add(src.h_ee_center_y_hc.get());
    dst.h_ee_vertex_z_hc->Add(src.h_ee_vertex_z_hc.get());

    dst.physics->GetEnergyVsModuleHist()->Add(src.physics->GetEnergyVsModuleHist());
    dst.physics->GetEnergyVsThetaHist()->Add(src.physics->GetEnergyVsThetaHist());
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
    for (int i = 1; i <= merged->h_ep_ee_ratio->GetNbinsX(); i++) {
        double ep = merged->h_ep_yield->GetBinContent(i);
        double ee = merged->h_ee_yield->GetBinContent(i);
        merged->h_ep_ee_ratio->SetBinContent(i, ee > 0. ? ep / ee : 0.);
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

    float hycal_vertex_z = fitAndDraw(merged->h_ee_vertex_z.get(), "Poscalib_result/hycal_vertex_z", 100.);
    float hycal_center_x = fitAndDraw(merged->h_ee_center_x.get(), "Poscalib_result/hycal_center_x", 2.);
    float hycal_center_y = fitAndDraw(merged->h_ee_center_y.get(), "Poscalib_result/hycal_center_y", 2.);

    // --- write output ---
    outfile.cd();
    hit_pos->Write();

    merged->h2_ep_hits->Write();
    merged->h2_ee_hits->Write();
    merged->h_ee_center_x->Write();
    merged->h_ee_center_y->Write();
    merged->h_ee_vertex_z->Write();
    merged->h2_ep_hits_hc->Write();
    merged->h2_ee_hits_hc->Write();
    merged->h_ee_center_x_hc->Write();
    merged->h_ee_center_y_hc->Write();
    merged->h_ee_vertex_z_hc->Write();

    outfile.mkdir("energy_plots"); outfile.cd("energy_plots");
    if (physics.GetEnergyVsModuleHist()) physics.GetEnergyVsModuleHist()->Write();
    if (physics.GetEnergyVsThetaHist())  physics.GetEnergyVsThetaHist()->Write();
    h_1cl->Write(); h_2cl->Write(); h_all->Write(); h_tot->Write();
    h2_energy_theta_ep_ee->Write();
    merged->h2_ep_E_angle->Write();
    merged->h2_ee_E_angle->Write();
    merged->h2_ep_E_angle_hc->Write();
    merged->h2_ee_E_angle_hc->Write();

    outfile.cd();
    outfile.mkdir("physics_yields"); outfile.cd("physics_yields");
    merged->h_ep_yield->Write();
    merged->h_ee_yield->Write();
    merged->h_ep_ee_ratio->Write();

    outfile.cd();
    outfile.mkdir("moller_analysis"); outfile.cd("moller_analysis");
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
    if (!hist) {
        std::cerr << "Cannot fit a null histogram for " << out_path << "\n";
        return std::numeric_limits<float>::quiet_NaN();
    }

    TCanvas c("", "", 800, 600);
    const float peak = hist->GetBinCenter(hist->GetMaximumBin());
    int fit_status = -1;
    if (hist->GetEntries() > 0)
        fit_status = hist->Fit("gaus", "rq", "", peak-fit_range, peak+fit_range);
    TF1 *gaus = hist->GetFunction("gaus");
    const bool fit_succeeded = fit_status == 0 && gaus;

    hist->Draw();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    if (fit_succeeded) {
        latex.DrawLatex(0.15, 0.85,
            Form("%.2f mm +- %.2f mm", gaus->GetParameter(1), gaus->GetParError(1)));
    } else {
        latex.DrawLatex(0.15, 0.85, "Gaussian fit failed");
        std::cerr << "Gaussian fit failed for " << hist->GetName()
                  << " (status " << fit_status << "); using peak position "
                  << peak << " mm\n";
    }
    fs::create_directories(fs::path(out_path).parent_path());
    c.SaveAs((out_path + ".png").c_str());

    return fit_succeeded ? gaus->GetParameter(1) : peak;
}
