// A tool to transform the Geant4 ouput into the replay_recon format,
// include matching between HyCal clusters and GEM hits
// can be used to compare with the replay reconstruction result and test your analysis
// scrirt, 
// Usage:
//   sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_luminosity(nb^-1)> <ee_luminosity(nb^-1)> [options]
//   -o  output ROOT file (default: sim_recon.root)
// Example:
//   sim2replay ep.root ee.root 1e6 1e6 -o sim_recon.root -n 100000 

#include "MatchingTools.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TChain.h>

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

// --- geometry constants (can be made configurable) ---
const float gem_z[4] = {5407.f + 39.71f/2, 5407.f - 39.71f/2,
                        5807.f + 39.71f/2, 5807.f - 39.71f/2};
const float hycal_z = 6225.f;

double beamE = 3500.0; // MeV, can be made configurable

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;
//event data structure for the input Geant4 simulation output
struct SimEventData {
    int VD_n = 0; // number of hits in the virtual detector
    double VD_x[500] = {};
    double VD_y[500] = {};
    double VD_z[500] = {};
    double VD_E[500] = {};
    int GEM_n = 0; // number of hits in the GEMs
    int GEM_id[500] = {};
    double GEM_x_in[500] = {};
    double GEM_y_in[500] = {};
    double GEM_z_in[500] = {};
    double GEM_x_out[500] = {};
    double GEM_y_out[500] = {};
    double GEM_z_out[500] = {};
    double GEM_edep[500] = {};
};
void setupSimBranches(TTree *tree, SimEventData &ev)
{   
    tree->SetBranchAddress("VD.N", &ev.VD_n);
    tree->SetBranchAddress("VD.X", ev.VD_x);
    tree->SetBranchAddress("VD.Y", ev.VD_y);
    tree->SetBranchAddress("VD.Z", ev.VD_z);
    tree->SetBranchAddress("VD.P", ev.VD_E);
    tree->SetBranchAddress("GEM.N", &ev.GEM_n);
    tree->SetBranchAddress("GEM.DID", ev.GEM_id);
    tree->SetBranchAddress("GEM.X", ev.GEM_x_in);
    tree->SetBranchAddress("GEM.Y", ev.GEM_y_in);
    tree->SetBranchAddress("GEM.Z", ev.GEM_z_in);
    tree->SetBranchAddress("GEM.Xout", ev.GEM_x_out);
    tree->SetBranchAddress("GEM.Yout", ev.GEM_y_out);
    tree->SetBranchAddress("GEM.Zout", ev.GEM_z_out);
    tree->SetBranchAddress("GEM.Edep", ev.GEM_edep);
}

static std::vector<std::string> collectRootFiles(const std::string &path);
float EResolution(float E);

// Find the HyCal module ID (PrimEx ID) whose area contains (x, y).
// Returns the module's PrimEx ID, or -1 if not found.
static int findModuleID(const fdec::HyCalSystem &hycal, double x, double y)
{
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        if (std::abs(x - m.x) <= m.size_x / 2. &&
            std::abs(y - m.y) <= m.size_y / 2.)
            return m.id;
    }
    return -1;
}

int main (int argc, char *argv[])
{
    std::string ep_file_dir, ee_file_dir;
    double ep_lumi = 0., ee_lumi = 0.;
    std::string outName = "sim_recon.root";
    int opt;
    while ((opt = getopt(argc, argv, "o:")) != -1) {
        switch (opt) {
            case 'o': outName = optarg; break;
        }
    }
    if (optind + 4 > argc) {
        std::cerr << "Usage: sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_lumi(nb^-1)> <ee_lumi(nb^-1)> [options]\n";
        std::cerr << "Options:\n  -o  output ROOT file (default: sim_recon.root)\n";
        return 1;
    }
    ep_file_dir = argv[optind];
    ee_file_dir = argv[optind + 1];
    ep_lumi = std::stod(argv[optind + 2]);
    ee_lumi = std::stod(argv[optind + 3]);

    // collect input files (can be files, directories)
    std::vector<std::string> ep_files, ee_files;
    auto f = collectRootFiles(ep_file_dir);
    ep_files.insert(ep_files.end(), f.begin(), f.end());
    f = collectRootFiles(ee_file_dir);
    ee_files.insert(ee_files.end(), f.begin(), f.end());

    if (ep_files.empty() || ee_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_lumi(nb^-1)> <ee_lumi(nb^-1)> [options]\n";
        return 1;
    }
    if(ep_lumi <= 0. || ee_lumi <= 0.) {
        std::cerr << "Invalid luminosity values.\n";
        std::cerr << "Usage: sim2replay <dir of input_ep.root> <dir of input_ee.root> <ep_lumi(nb^-1)> <ee_lumi(nb^-1)> [options]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);
    MatchingTools matching;

    // --- setup TChain and branches ---
    auto sim = std::make_unique<SimEventData>();
    TChain *chain_ee = new TChain("T");
    for (const auto &f : ee_files) {
        chain_ee->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree_ee = chain_ee;
    if (!tree_ee) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }
    setupSimBranches(tree_ee, *sim);
    TChain *chain_ep = new TChain("T");
    for (const auto &f : ep_files) {
        chain_ep->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree_ep = chain_ep;
    if (!tree_ep) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }
    setupSimBranches(tree_ep, *sim);

    TFile *outfile = TFile::Open(outName.c_str(), "RECREATE");
    // create TTree and branches for reconstructed data
    TTree *tree_out = new TTree("recon", "PRad2 replay reconstruction from G4");
    auto ev = std::make_unique<EventVars_Recon>();
    prad2::SetReconWriteBranches(tree_out, *ev, false); // false indicates not x17_mode

    // caculate luminosity and number of events to process for ep and ee
    double lumi = std::min(ep_lumi, ee_lumi);
    int N_ep = lumi / ep_lumi * tree_ep->GetEntries();
    int N_ee = lumi / ee_lumi * tree_ee->GetEntries();

    std::cerr << "Processing " << N_ep << " ep events and " << N_ee << " ee events (luminosity = " << lumi << " nb^-1)\n";

    // --- loop over events and fill output tree ---
    int ep_count = 0, ee_count = 0;
    for(int i = 0; i < N_ep + N_ee; i++){
        if(i + 1 % 4 == 0 && ep_count < N_ep) {
            tree_ep->GetEntry(ep_count);
            ep_count++;
        } else if(ee_count < N_ee) {
            tree_ee->GetEntry(ee_count);
            ee_count++;
        }
        if(ep_count == N_ep && ee_count < N_ee) {
            tree_ee->GetEntry(ee_count);
            ee_count++;
        }
        if(ee_count == N_ee && ep_count < N_ep) {
            tree_ep->GetEntry(ep_count);
            ep_count++;
        }
        if(ee_count == N_ee && ep_count == N_ep) break;
        if((i + 1) % 1000 == 0) {
            std::cerr << "Processing event " << i + 1 << "/" << N_ep + N_ee << "\r" << std::flush;
        }

        // Here you can fill the ev structure with the HyCal virtual plane hit and GEM hit information
        // For example, you can use the matching tools to match HyCal clusters and GEM hits
        // and fill the ev.matchG_x, etc. arrays accordingly
        // This part will depend on how you want to do the reconstruction and matching based on the input simulation data structure
        // You can also use the physics tools to calculate expected energies, angles, etc. for the clusters and hits, and fill the ev structure with those values as well

        //first total energy cut
        *ev = EventVars_Recon{}; // reset all fields for this event
        float total_energy = 0.f;
        for(int j = 0; j < sim->VD_n; j++){
            if(sim->VD_E[j] < 1. / 300.f * beamE) continue; // add some energy threshold to reduce noise
            sim->VD_E[j] += gRandom->Gaus(0, EResolution(sim->VD_E[j])); // add some energy smearing
            total_energy += sim->VD_E[j];
        }
        if(total_energy < 0.5f * beamE) continue; // only keep events with total energy above 0.5*beamE

        ev->event_num = i;
        ev->trigger_type = 8; // primary_bit 8
        ev->trigger_bits = prad2::TBIT_sum; // SSP raw sum
        ev->timestamp = 0;
        ev->total_energy = total_energy;

        for(int j = 0; j < sim->VD_n; j++){
            if(sim->VD_E[j] < 1. / 300.f * beamE) continue; // add some energy threshold to reduce noise
            ev->cl_x[ev->n_clusters] = float(sim->VD_x[j] + gRandom->Gaus(0, 2.6/sqrt(sim->VD_E[j] / 1000.f)));
            ev->cl_y[ev->n_clusters] = float(sim->VD_y[j] + gRandom->Gaus(0, 2.6/sqrt(sim->VD_E[j] / 1000.f)));
            ev->cl_z[ev->n_clusters] = float(hycal_z);
            ev->cl_energy[ev->n_clusters] = float(sim->VD_E[j]);
            ev->cl_nblocks[ev->n_clusters] = 1;
            ev->cl_center[ev->n_clusters] = findModuleID(hycal, sim->VD_x[j], sim->VD_y[j]);
            ev->cl_flag[ev->n_clusters] = ( 1 << 1); // set flag to PbWO4
            ev->n_clusters++;
        }

        for(int j = 0; j < sim->GEM_n; j++){
            if(sim->GEM_edep[j] < 26.e-6*2.) continue; // add some energy threshold to reduce noise
            ev->det_id[ev->n_gem_hits] = sim->GEM_id[j];
            ev->gem_x[ev->n_gem_hits] = float(0.5*(sim->GEM_x_in[j] + sim->GEM_x_out[j]) + gRandom->Gaus(0, 0.07)); // add some position smearing
            ev->gem_y[ev->n_gem_hits] = float(0.5*(sim->GEM_y_in[j] + sim->GEM_y_out[j]) + gRandom->Gaus(0, 0.07));
            ev->gem_x_charge[ev->n_gem_hits] = 0.f;
            ev->gem_y_charge[ev->n_gem_hits] = 0.f;
            ev->gem_x_peak[ev->n_gem_hits] = 0.f;
            ev->gem_y_peak[ev->n_gem_hits] = 0.f;
            ev->gem_x_size[ev->n_gem_hits] = 3; // set some default cluster size
            ev->gem_y_size[ev->n_gem_hits] = 3;
            ev->n_gem_hits++;
        }

        //do matching between HyCal clusters and GEM hits
        //store all the hits on HyCal and GEMs in this event
        std::vector<HCHit> hc_hits;
        std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
        for (int i = 0; i < ev->n_clusters; i++)
            hc_hits.push_back({ev->cl_x[i], ev->cl_y[i], ev->cl_z[i], ev->cl_energy[i], ev->cl_center[i], ev->cl_flag[i]});
        for (int i = 0; i < ev->n_gem_hits; ++i)
            gem_hits[ev->det_id[i]].push_back(GEMHit{ev->gem_x[i], ev->gem_y[i], gem_z[ev->det_id[i]], ev->det_id[i]});

        //GetProjection(hc_hits, 6225.f);

        matching.SetMatchRange(10.f); // matching radius in mm, 15mm default
        //matching.SetSquareSelection(true); // use square cut instead of circular cut
        std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);

        // save transformed HyCal positions for all clusters (regardless of match)
        for (int i = 0; i < ev->n_clusters; ++i) {
            for (int j = 0; j < 2; j++) {
                ev->matchGEMx[i][j] = -999.f;
                ev->matchGEMy[i][j] = -999.f;
                ev->matchGEMz[i][j] = -999.f;
            }
            ev->matchFlag[i] = 0;
        }
        ev->matchNum = std::min((int)matched_hits.size(), prad2::kMaxClusters);
        for (int i = 0; i < ev->matchNum; i++) {
            int cl_idx = matched_hits[i].hycal_idx;
            ev->matchFlag[cl_idx] = matched_hits[i].mflag;
            for (int j = 0; j < 2; j++) {
                ev->matchGEMx[cl_idx][j] = matched_hits[i].gem[j].x;
                ev->matchGEMy[cl_idx][j] = matched_hits[i].gem[j].y;
                ev->matchGEMz[cl_idx][j] = matched_hits[i].gem[j].z;
            }
            // quick access arrays
            ev->mHit_E[i]  = matched_hits[i].hycal_hit.energy;
            ev->mHit_x[i]  = matched_hits[i].hycal_hit.x;
            ev->mHit_y[i]  = matched_hits[i].hycal_hit.y;
            ev->mHit_z[i]  = matched_hits[i].hycal_hit.z;
            for (int j = 0; j < 2; j++) {
                ev->mHit_gx[i][j]  = matched_hits[i].gem[j].x;
                ev->mHit_gy[i][j]  = matched_hits[i].gem[j].y;
                ev->mHit_gz[i][j]  = matched_hits[i].gem[j].z;
                ev->mHit_gid[i][j] = matched_hits[i].gem[j].det_id;
            }
        }

        tree_out->Fill();
    }
    outfile->Write();
    delete outfile;
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find(".root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

float EResolution(float E)
{
    // Example energy resolution function, can be modified based on actual detector performance
    return E * 0.026f / sqrt(E / 1000.f); // 2.6% at 1 GeV, scaling with sqrt(E)
}

