// epCalib.cpp: tool to get calibration constants from photon beam in calibration runs
// on each module, read the rawdata.root(peak mode) file, reconstruction and fit the peak,
// get the ratio of expected/measured peak position, and write to a database file. 
//=============================================================================
//
// Usage: gamma_calib <input_raw.root|dir> input_calib_file output_calib_file
//                         --option       [-o output_root_file] [-n max_events]
//
// Reads rawdata(adc level).root (peak mode), runs HyCal clustering, fills per-module energy histograms
//=============================================================================

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "EventData.h"
#include "EventData_io.h"
#include "load_daq_config.h"

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <algorithm>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars       = prad2::RawEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);

int main(int argc, char *argv[])
{
    std::string in_calib_file, out_calib_file, output_root_file;
    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    int max_events = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:")) != -1) {
        switch (opt) {
            case 'o': output_root_file = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
    }

    if(optind + 2 != argc) {
        std::cerr << "Usage: " << argv[0] << " <input_raw.root|dir> <input_calib_file> <output_calib_file>\n"
                  << "Options:\n"
                  << "  -o output_root_file   Output ROOT file for reconstructed variables (optional)\n"
                  << "  -n max_events         Maximum number of events to process (default: all)\n";
        return 1;
    }
    std::string input = argv[optind];
    in_calib_file = argv[optind + 1];
    out_calib_file = argv[optind + 2];

    // collect input files (can be files, directories)
    std::vector<std::string> input_files = collectRootFiles(input);
    if (input_files.empty()) {
        std::cerr << "No input files found in " << input << "\n";
        return 1;
    }

    // --- setup TChain and branches ---
    TChain *chain = new TChain("events");
    for (const auto &f : input_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'events' in input files\n";
        return 1;
    }
    auto ev = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(tree, *ev);

    // --- setup output ROOT file
    if(output_root_file.empty()) output_root_file = "gamma_calib_info.root";
    TFile *outfile = TFile::Open(output_root_file.c_str(), "RECREATE");

    //setup for reconstruction
    fdec::HyCalSystem hycal;
    evc::DaqConfig daq_cfg;
    std::string daq_config_file = db_dir + "/daq_config.json"; // default DAQ config for PRad2
    evc::load_daq_config(daq_config_file, daq_cfg);
    hycal.Init(db_dir + "/hycal_map.json");

    int nmatched = hycal.LoadCalibration(in_calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << in_calib_file << " (" << nmatched << " modules)\n";

    analysis::PhysicsTools physics(hycal);
    fdec::ClusterConfig cl_cfg;

    TH1F *ratio_module_all = new TH1F("ratio_all", "Ratio of Expected/Measured Peak Position for All Modules;Ratio;Modules", 100, 0, 3);


    //loop over events
    int N_events = tree->GetEntries();
    N_events = (max_events > 0 && max_events < N_events) ? max_events : N_events;
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cl_cfg);
    for (int i = 0; i < N_events; i++) {
        tree->GetEntry(i);

        // TODO: use config-driven trigger filter (monitor_config.json "physics" section
        // accept_trigger_bits/reject_trigger_bits) instead of hardcoded bit check.
        // Currently drops all non-SSP_RawSum events, including LMS.
        if (!(ev->trigger_bits & prad2::TBIT_sum)) continue;
        if (ev->trigger_bits & prad2::TBIT_lms) continue;

        // --- HyCal clustering ---
        clusterer.Clear();
        for(int j = 0; j < ev->nch; j++){
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev->npeaks[j] <= 0) continue;
            float adc = ev->peak_integral[j][0];
            float energy = (mod->cal_factor > 0.) ?
                static_cast<float>(mod->energize(adc)) : adc * 0.2f;
            clusterer.AddHit(mod->index, energy, 0.f); // time info not available in raw branches, set to 0
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        //TO DO: event selection
        bool tagger_t1 = true;
        if(hits.size() == 1 && tagger_t1) {
            physics.FillModuleEnergy(hits[0].center_id, hits[0].energy);
        }
    }
    //calculate calibration constants for each module, output calibration file to database
    for (int m = 0; m < hycal.module_count(); m++) {
        auto mod = hycal.module(m);
        auto [peak, sigma, chi2] = physics.FitPeakResolution(mod->id);
        if (peak > 0 && sigma > 0) {
            std::string name = mod->name;
            if(name[0] != 'W') continue; // only calibrate PbWO4
            float expected_peak = physics.ExpectedPeakPosition(name);
            if (expected_peak > 0) {
                float ratio = expected_peak / peak;
                float calib_const = ratio * mod->cal_factor;
                hycal.SetCalibConstant(mod->id, calib_const);
                ratio_module_all->Fill(ratio);
                std::cerr << "Module " << name << ": peak = " << peak << ", expected = " << expected_peak << ", ratio = " << ratio << ", new calib const = " << calib_const << "\n";
            }
        }
    }

    //write the new calibration constants to database file
    hycal.PrintCalibConstants(out_calib_file);
    outfile->Close();
}


// ── Helpers ──────────────────────────────────────────────────────────────
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
