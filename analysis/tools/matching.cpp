// A quick check tool to test the matching result between HyCal clusters and GEM hits in the replay output
// Usage:
//   matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]
//   -o  output ROOT file (default: input filename with _matching.root suffix)
//   -n  max events to process (default: all)
// Example:
//   matching recon.root -o recon_matching.root -n 10000
//   matching recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
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

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    
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

    // --- load detector geometry config from JSON ---
    std::string run_str = get_run_str(root_files[0]);
    int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(dbDir + "/runinfo/general.json", run_num);

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
        outName.ReplaceAll("_recon.root", "_matching.root");
    }
    TFile outfile(outName, "RECREATE");
    int hitN = 0;
    float HC_X[100], HC_Y[100], HC_Z[100], HC_Energy[100];
    float Gup_X[100], Gup_Y[100], Gup_Z[100];
    float Gdown_X[100], Gdown_Y[100], Gdown_Z[100];
    uint16_t center_id[100];
    uint32_t flag[100];

    TTree *outTree = new TTree("matching", "HyCal-GEM matching results");
    outTree->Branch("hitN", &hitN, "hitN/I");
    outTree->Branch("HC_X", HC_X, "HC_X[hitN]/F");
    outTree->Branch("HC_Y", HC_Y, "HC_Y[hitN]/F");
    outTree->Branch("HC_Z", HC_Z, "HC_Z[hitN]/F");
    outTree->Branch("HC_Energy", HC_Energy, "HC_Energy[hitN]/F");
    outTree->Branch("Gup_X", Gup_X, "Gup_X[hitN]/F");
    outTree->Branch("Gup_Y", Gup_Y, "Gup_Y[hitN]/F");
    outTree->Branch("Gup_Z", Gup_Z, "Gup_Z[hitN]/F");
    outTree->Branch("Gdown_X", Gdown_X, "Gdown_X[hitN]/F");
    outTree->Branch("Gdown_Y", Gdown_Y, "Gdown_Y[hitN]/F");
    outTree->Branch("Gdown_Z", Gdown_Z, "Gdown_Z[hitN]/F");
    outTree->Branch("center_id", center_id, "center_id[hitN]/s");
    outTree->Branch("flag", flag, "flag[hitN]/i");

    // --- event loop  ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;
    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "Reading " << i << " events / " << N << " total events\r" << std::flush;

        //store all the hits on HyCal and GEMs in this event
        std::vector<HCHit> hc_hits;
        std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
        for( int j = 0; j < ev.n_clusters; j++) {
            hc_hits.push_back(HCHit{ev.cl_x[j], ev.cl_y[j], ev.cl_z[j],
                               ev.cl_energy[j], ev.cl_center[j], ev.cl_flag[j]});
        }
        for (int j = 0; j < ev.n_gem_hits; j++) {
            gem_hits[ev.det_id[j]].push_back(GEMHit{ev.gem_x[j], ev.gem_y[j], 0.f, ev.det_id[j]});
        }

        // ev.cl_x / ev.gem_x are already lab-frame (Replay applied the
        // DetectorTransform before writing the tree), so no further frame
        // change here.

        //then matching between GEM hits and HyCal clusters
            //optional settings
        //matching.SetMatchRange(5.0f); // matching radius in mm
        //matching.SetSquareSelection(true); // use square cut instead of circular cut
        std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
        //show how to access the matching result
        for (auto &m : matched_hits) {
            HCHit hycal_hit = m.hycal_hit;  //the HyCal cluster be matched
            GEMHit gem_up_hit = m.gem[0];    //best-matched GEM hit from upstream pair (GEM1/GEM2)
            GEMHit gem_down_hit = m.gem[1]; //best-matched GEM hit from downstream pair (GEM3/GEM4)
            std::vector<GEMHit> gem1_matches = m.gem1_hits;
            std::vector<GEMHit> gem2_matches = m.gem2_hits;
            std::vector<GEMHit> gem3_matches = m.gem3_hits;
            std::vector<GEMHit> gem4_matches = m.gem4_hits;

            int hycal_idx = m.hycal_idx;  //index of the cluster in the original vector

            //projection to same z plane (e.g. Hycal crystal surface)
            GetProjection(hycal_hit, gRunConfig.hycal_z);
            GetProjection(gem_up_hit, gRunConfig.hycal_z);
            GetProjection(gem_down_hit, gRunConfig.hycal_z);

            //save the matching result into output tree
            if (hitN < 100) { // safety check for array size
                HC_X[hitN] = hycal_hit.x;
                HC_Y[hitN] = hycal_hit.y;
                HC_Z[hitN] = hycal_hit.z;
                HC_Energy[hitN] = hycal_hit.energy;
                Gup_X[hitN] = gem_up_hit.x;
                Gup_Y[hitN] = gem_up_hit.y;
                Gup_Z[hitN] = gem_up_hit.z;
                Gdown_X[hitN] = gem_down_hit.x;
                Gdown_Y[hitN] = gem_down_hit.y;
                Gdown_Z[hitN] = gem_down_hit.z;
                center_id[hitN] = hycal_hit.center_id;
                flag[hitN] = hycal_hit.flag;
                hitN++;
            }
        }
        outTree->Fill();
        hitN = 0; // reset for next event
    }
    std::cerr << "Finished processing " << N << " events.\n";
    outTree->Write();
    outfile.Close();
    std::cerr << "Output saved to " << outName << "\n";
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