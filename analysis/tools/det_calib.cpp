
#include "PhysicsTools.h"
#include "HyCalSystem.h"
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
#include <TLatex.h>
#include <TCanvas.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);
void projectToHyCalSurface(MollerData &m_data, float hycal_z);
float fitAndDraw(TH1F* hist, const std::string& out_path, const float survey_position, const float fit_range = 4.);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    std::string run_config;

    float HyCal_shift_x = 0.f;
    float HyCal_shift_y = 0.f;
    float GEM_shift_x[4] = {0.f, 0.f, 0.f, 0.f};
    float GEM_shift_y[4] = {0.f, 0.f, 0.f, 0.f};

    int max_events = -1;

    enum { OPT_HX=256, OPT_HY,
           OPT_G0X, OPT_G0Y, OPT_G1X, OPT_G1Y,
           OPT_G2X, OPT_G2Y, OPT_G3X, OPT_G3Y };
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
    while ((opt = getopt_long(argc, argv, "o:n:c:", long_opts, &longidx)) != -1) {
        switch (opt) {
            case 'o': output     = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'c': run_config = optarg; break;
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
    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: det_calib <input_recon.root|dir> [more files...]\n"
                     "    [-o out.root] [-n max_events] [-c run_config.json]\n"
                     "    [--hx val] [--hy val]                  HyCal shift (mm)\n"
                     "    [--g0x val] [--g0y val] .. [--g3x val] [--g3y val]  GEM shift (mm)\n";
        return 1;
    }
    // extract run number from first input file name (e.g. prad_023626.00000_recon.root -> 23626)
    std::string run_str = get_run_str(root_files[0]);
    int run_num = get_run_int(root_files[0]);

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- load detector geometry config from JSON ---
    if (run_config.empty()) {
        run_config = dbDir + "/runinfo/general.json";
    }
    RunConfig geo = LoadRunConfig(run_config, run_num);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

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

    //output histograms for calibration results
    TH1F *vertex_hycal = new TH1F("vertex_hycal", "Moller vertex z distance HyCal;Z (mm);Counts", 600, 5600, 6800);
    TH2F *center_hycal = new TH2F("center_hycal", "Moller center distribution HyCal;X (mm);Y (mm)", 200, -100, 100, 200, -100, 100);
    TH1F *center_hycal_x = new TH1F("center_hycal_x", "Moller center X distribution HyCal;X (mm);Counts", 80*4, -20, 20);
    TH1F *center_hycal_y = new TH1F("center_hycal_y", "Moller center Y distribution HyCal;Y (mm);Counts", 80*4, -20, 20);

    TH1F *vertex_gem[4];
    TH2F *center_gem[4];
    TH1F *center_gem_x[4];
    TH1F *center_gem_y[4];
    for (int d = 0; d < 4; d++) {
        vertex_gem[d] = new TH1F(Form("vertex_gem%d", d), Form("Moller vertex z distance GEM%d;Z (mm);Counts", d), 1000, 5200, 6200);
        center_gem[d] = new TH2F(Form("center_gem%d", d), Form("Moller center distribution GEM%d;X (mm);Y (mm)", d), 200, -100, 100, 200, -100, 100);
        center_gem_x[d] = new TH1F(Form("center_gem_x%d", d), Form("Moller center X distribution GEM%d;X (mm);Counts", d), 80*5, -20, 20);
        center_gem_y[d] = new TH1F(Form("center_gem_y%d", d), Form("Moller center Y distribution GEM%d;Y (mm);Counts", d), 80*4, -20, 20);
    }

    TH2F *hits_hycal = new TH2F("hits_hycal", "Moller hit HyCal positions;X (mm);Y (mm)", 400, -400, 400, 400, -400, 400);
    TH2F *hits_gem[4];
    for (int d = 0; d < 4; d++) {
        hits_gem[d] = new TH2F(Form("hits_gem%d", d), Form("Moller hit GEM%d positions;X (mm);Y (mm)", d), 400, -400, 400, 400, -400, 400);
    }

    // --- output file ---
    TString outName = output;
    if (outName.IsNull()) {
        outName = root_files[0];
        outName.ReplaceAll("_recon.root", "_posCalib.root");
    }
    TFile outfile(outName, "RECREATE");

    MollerData hycal_mollers;
    MollerData gem_mollers[4];

    // --- event loop : select Moller events on HyCal and each GEM plane ---
    int N = tree->GetEntries();
    if (max_events > 0 && max_events < N) N = max_events;

    for (int i = 0; i < N; i++) {
        tree->GetEntry(i);
        if (i % 1000 == 0)
            std::cerr << "\rPass 1: " << i << " / " << N << std::flush;

        bool good_moller = false;
        if(ev.matchNum == 2){
            float Epair = ev.mHit_E[0] + ev.mHit_E[1];
            if(geo.Ebeam <= 0.f){
                std::cerr << "Error: Ebeam not set, cannot apply Moller energy cut.\n";
                break;
            }
            if (std::abs(Epair - geo.Ebeam) < 4.f * geo.Ebeam * 0.025f / std::sqrt(geo.Ebeam / 1000.f)) {
                good_moller = true;
            }
        }
        if(!good_moller) continue;

        //have selected good Moller events for further analysis
        MollerEvent h_m;
        MollerEvent g_m;
        
        h_m = MollerEvent(
                {ev.mHit_x[0], ev.mHit_y[0], ev.mHit_z[0], ev.mHit_E[0]},
                {ev.mHit_x[1], ev.mHit_y[1], ev.mHit_z[1], ev.mHit_E[1]});
        hycal_mollers.push_back(h_m);
        
        // select two moller on one chamber for upstream GEMs 
        if(ev.mHit_gid[0][0] == ev.mHit_gid[1][0]){
            g_m = MollerEvent(
                {ev.mHit_gx[0][0], ev.mHit_gy[0][0], ev.mHit_gz[0][0], ev.mHit_E[0]},
                {ev.mHit_gx[1][0], ev.mHit_gy[1][0], ev.mHit_gz[1][0], ev.mHit_E[1]});
            int det_id = ev.mHit_gid[0][0];
            if(det_id >= 0 && det_id < 4) gem_mollers[det_id].push_back(g_m);
            else std::cerr << "Warning: Invalid GEM det_id " << det_id << " in event " << ev.event_num << "\n";
        }

        // select two moller on one chamber for downstream GEMs
        if(ev.mHit_gid[0][1] == ev.mHit_gid[1][1]){
            g_m = MollerEvent(
                {ev.mHit_gx[0][1], ev.mHit_gy[0][1], ev.mHit_gz[0][1], ev.mHit_E[0]},
                {ev.mHit_gx[1][1], ev.mHit_gy[1][1], ev.mHit_gz[1][1], ev.mHit_E[1]});
            int det_id = ev.mHit_gid[0][1];
            if(det_id >= 0 && det_id < 4) gem_mollers[det_id].push_back(g_m);
            else std::cerr << "Warning: Invalid GEM det_id " << det_id << " in event " << ev.event_num << "\n";
        }

        hits_hycal->Fill(h_m.first.x, h_m.first.y);
        hits_hycal->Fill(h_m.second.x, h_m.second.y);
    }

    // After collecting Moller events, analyze them for detector calibration
    //summary of Moller events on each detector plane
    std::cerr << "\nSummary of selected Moller events:\n";
    std::cerr << "HyCal: " << hycal_mollers.size() << " events\n";
    for (int d = 0; d < 4; d++) {
        std::cerr << "GEM " << d << ": " << gem_mollers[d].size() << " events\n";
    }

    //hycal Moller events
    //projectToHyCalSurface(hycal_mollers, hycal_z); //project to HyCal surface
    //move to beam center coordinates
    TransformDetData(hycal_mollers, HyCal_shift_x, HyCal_shift_y, 0.f);
    for (int i = 0; i < hycal_mollers.size(); i++) {
        vertex_hycal->Fill(physics.GetMollerZdistance(hycal_mollers[i], geo.Ebeam));
        if (i >= 1) {
            auto c = physics.GetMollerCenter(hycal_mollers[i-1], hycal_mollers[i]);
            center_hycal->Fill(c[0], c[1]);
            center_hycal_x->Fill(c[0]);
            center_hycal_y->Fill(c[1]);
        }
    }

    //gem Moller events
    for (int d = 0; d < 4; d++) {
        TransformDetData(gem_mollers[d], GEM_shift_x[d], GEM_shift_y[d], 0.f);
        for (int i = 0; i < gem_mollers[d].size(); i++) {
            vertex_gem[d]->Fill(physics.GetMollerZdistance(gem_mollers[d][i], geo.Ebeam));
            if (i >= 1) {
                auto c = physics.GetMollerCenter(gem_mollers[d][i-1], gem_mollers[d][i]);
                center_gem[d]->Fill(c[0], c[1]);
                center_gem_x[d]->Fill(c[0]);
                center_gem_y[d]->Fill(c[1]);
            }
            const auto &m = gem_mollers[d][i];
            hits_gem[d]->Fill(m.first.x, m.first.y);
            hits_gem[d]->Fill(m.second.x, m.second.y);
        }
    }

    //fit histograms, and get the beam position and vertex distance for each detector plane
    float hycal_vertex_z = fitAndDraw(vertex_hycal, "Poscalib_result/" + run_str +"/hycal_vertex_z", geo.hycal_z,  100.);
    float hycal_center_x = fitAndDraw(center_hycal_x, "Poscalib_result/" + run_str +"/hycal_center_x", geo.hycal_x+HyCal_shift_x, 2.);
    float hycal_center_y = fitAndDraw(center_hycal_y, "Poscalib_result/" + run_str +"/hycal_center_y", geo.hycal_y+HyCal_shift_y, 2.);
    float gem_vertex_z[4];
    float gem_center_x[4];
    float gem_center_y[4];
    for (int d = 0; d < 4; d++) {
        gem_vertex_z[d] = fitAndDraw(vertex_gem[d], "Poscalib_result/" + run_str + "/gem" + std::to_string(d) + "_vertex_z", geo.gem_z[d], 25.);
        gem_center_x[d] = fitAndDraw(center_gem_x[d], "Poscalib_result/" + run_str + "/gem" + std::to_string(d) + "_center_x", geo.gem_x[d]+GEM_shift_x[d], 0.3);
        gem_center_y[d] = fitAndDraw(center_gem_y[d], "Poscalib_result/" + run_str + "/gem" + std::to_string(d) + "_center_y", geo.gem_y[d]+GEM_shift_y[d], 1.);
    }
    //print summary of calibration results
    std::cerr << "HyCal vertex z distance: " << hycal_vertex_z << " mm (survey position " << geo.hycal_z << " mm)" << "\n";
    std::cerr << "HyCal center x: " << hycal_center_x << " mm (survey position " << geo.hycal_x << " mm)\n";
    std::cerr << "HyCal center y: " << hycal_center_y << " mm (survey position " << geo.hycal_y << " mm)\n";
    
    for (int d = 0; d < 4; d++) {
        std::cerr << "GEM " << d << " vertex z distance: " << gem_vertex_z[d] << " mm (survey position " << geo.gem_z[d] << " mm)\n";
        std::cerr << "GEM " << d << " center x: " << gem_center_x[d] << " mm (survey position " << geo.gem_x[d] << " mm)\n";
        std::cerr << "GEM " << d << " center y: " << gem_center_y[d] << " mm (survey position " << geo.gem_y[d] << " mm)\n";
    }

    geo.hycal_z = hycal_vertex_z;
    geo.hycal_x += hycal_center_x;
    geo.hycal_y += hycal_center_y;
    for (int d = 0; d < 4; d++) {
        geo.gem_z[d] = gem_vertex_z[d];
        geo.gem_x[d] += gem_center_x[d];
        geo.gem_y[d] += gem_center_y[d];
    }
    //write back the updated geometry config to JSON file
    //WriteRunConfig(run_config, run_num, geo);
    

    //save histograms
    outfile.cd();
    vertex_hycal->Write();
    center_hycal->Write();
    center_hycal_x->Write();
    center_hycal_y->Write();
    for (int d = 0; d < 4; d++) {
        vertex_gem[d]->Write();
        center_gem[d]->Write();
        center_gem_x[d]->Write();
        center_gem_y[d]->Write();
    }
    hits_hycal->Write();
    for (int d = 0; d < 4; d++) {
        hits_gem[d]->Write();
    }
    outfile.Close();
    std::cerr << "Calibration histograms saved to " << outName << "\n";

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

void projectToHyCalSurface(MollerData &m_data, float hycal_z)
{
    //project the Moller event from target center(z = 0) to the HyCal surface (z = hycal_z)
    for (auto &evt : m_data) {
        for (auto *dp : {&evt.first, &evt.second}) {
            float scale = hycal_z / dp->z;
            dp->x = dp->x * scale;
            dp->y = dp->y * scale;
            dp->z = hycal_z;
        }
    }
}

float fitAndDraw(TH1F* hist, const std::string& out_path, const float survey_position, const float fit_range){
    TCanvas *c = new TCanvas("", "", 800, 600);
    float mean = hist->GetBinCenter(hist->GetMaximumBin());
    hist->Fit("gaus", "rq", "", mean-fit_range, mean+fit_range);
    hist->Draw();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.15, 0.85, Form("%.2f mm +- %.2f mm", hist->GetFunction("gaus")->GetParameter(1), hist->GetFunction("gaus")->GetParError(1)));
    latex->DrawLatex(0.15, 0.80, Form("Calibrated position(except z): %.2f mm", survey_position));
    fs::create_directories(fs::path(out_path).parent_path());
    c->SaveAs((out_path + ".png").c_str());
    delete c;

    return hist->GetFunction("gaus")->GetParameter(1);
}