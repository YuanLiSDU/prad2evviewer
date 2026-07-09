
//=============================================================================
// det_calib — detector position calibration via Møller scattering
//
// Process: Read evio files, decode and reconstruct, judge if it's 
// Moller events, analyze and fill histgroams, find out the detector alignment
//
// Usage:
//   det_calib <evio_file_or_dir> [more files/dirs...]
//             -o output_dir
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"
#include "MatchingTools.h"
#include "PipelineBuilder.h"
#include "PulseTemplateStore.h"
#include "gain_factor.h"

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
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

struct TaggedMollerEvent {
    DataPoint first;
    DataPoint second;
    int event_num = -1;

    TaggedMollerEvent() = default;
    TaggedMollerEvent(const DataPoint &first_, const DataPoint &second_,
                      int event_num_)
        : first(first_), second(second_), event_num(event_num_) {}
};

using TaggedMollerData = std::vector<TaggedMollerEvent>;

static MollerEvent ToMollerEvent(const TaggedMollerEvent &event)
{
    return MollerEvent(event.first, event.second);
}

// ── forward declarations ──────────────────────────────────────────────────
static std::vector<std::string> collectEvioFiles(const std::string &path);
static bool ProcessEVIO(const std::string &input_evio,
                        RunConfig &run_config,
                        const std::string &db_dir,
                        const std::string &run_config_file,
                        const std::string &daq_config_file,
                        const std::string &hycal_map_file,
                        const std::string &gem_ped_file,
                        float zerosup_override,
                        int max_events,
                        TaggedMollerData &hycal_mollers,
                        std::array<TaggedMollerData, 4> &gem_mollers);

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

// ── Main ──────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string daq_config, daq_map, gem_ped_file, output_dir, run_config;
    float zerosup_override = 5.f;
    int  max_files   = -1;
    int  num_threads = 4;
    int  max_events  = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:f:n:j:c:C:d:g:z:")) != -1) {
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
            default: return 1;
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
            "       [-c run_config.json] [-g gem_pedestal.json]\n";
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
    gRunConfig = LoadRunConfig(run_config, run_num);

    // Each worker writes only to the result slot belonging to its EVIO file.
    // This avoids locking the (potentially large) Moller vectors while events
    // are being reconstructed.  The main thread merges the slots after all
    // workers have finished.
    std::vector<TaggedMollerData> hycal_results(num_files);
    std::vector<std::array<TaggedMollerData, 4>> gem_results(num_files);
    std::vector<char> processed_ok(num_files, 0);
    std::atomic<int> next_file{0};
    std::atomic<int> errors{0};
    std::mutex io_mtx;

    auto worker = [&]() {
        while (true) {
            const int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            RunConfig file_config = gRunConfig;
            const bool ok = ProcessEVIO(
                evio_files[idx], file_config, db_dir, run_config, daq_config,
                daq_map, gem_ped_file, zerosup_override, max_events,
                hycal_results[idx], gem_results[idx]);
            processed_ok[idx] = ok ? 1 : 0;

            std::lock_guard<std::mutex> lock(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                          << evio_files[idx] << ": "
                          << hycal_results[idx].size() << " HyCal Mollers\n";
            } else {
                ++errors;
                std::cerr << "  [" << (idx + 1) << "/" << num_files
                          << "] FAILED: " << evio_files[idx] << "\n";
            }
        }
    };

    std::cout << "Processing " << num_files << " EVIO files with "
              << num_threads << " threads\n";
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i)
        threads.emplace_back(worker);
    for (auto &thread : threads)
        thread.join();

    // Aggregate in input-file order, so results are deterministic regardless
    // of the order in which worker threads completed.
    TaggedMollerData hycal_mollers;
    std::array<TaggedMollerData, 4> gem_mollers;
    size_t total_hycal = 0;
    std::array<size_t, 4> total_gem{};
    for (int i = 0; i < num_files; ++i) {
        if (!processed_ok[i]) continue;
        total_hycal += hycal_results[i].size();
        for (size_t det = 0; det < gem_mollers.size(); ++det)
            total_gem[det] += gem_results[i][det].size();
    }
    hycal_mollers.reserve(total_hycal);
    for (size_t det = 0; det < gem_mollers.size(); ++det)
        gem_mollers[det].reserve(total_gem[det]);

    for (int i = 0; i < num_files; ++i) {
        if (!processed_ok[i]) continue;
        hycal_mollers.insert(hycal_mollers.end(),
                             hycal_results[i].begin(), hycal_results[i].end());
        for (size_t det = 0; det < gem_mollers.size(); ++det) {
            gem_mollers[det].insert(gem_mollers[det].end(),
                                    gem_results[i][det].begin(),
                                    gem_results[i][det].end());
        }
    }

    std::cout << "Collected " << hycal_mollers.size()
              << " HyCal Moller events\n";
    for (size_t det = 0; det < gem_mollers.size(); ++det)
        std::cout << "  GEM " << det << ": " << gem_mollers[det].size()
                  << " Moller events\n";

    //Creat the histgrams
    TH2F *h2_hycal_hits = new TH2F("h2_hycal_hits", "HyCal Hits", 720, -360, 360, 720, -360, 360);
    TH2F *h2_gem0_hits  = new TH2F("h2_gem0_hits",  "GEM0 Hits",  720, -360, 360, 720, -360, 360);
    TH2F *h2_gem1_hits  = new TH2F("h2_gem1_hits",  "GEM1 Hits",  720, -360, 360, 720, -360, 360);
    TH2F *h2_gem2_hits  = new TH2F("h2_gem2_hits",  "GEM2 Hits",  720, -360, 360, 720, -360, 360);
    TH2F *h2_gem3_hits  = new TH2F("h2_gem3_hits",  "GEM3 Hits",  720, -360, 360, 720, -360, 360);

    TH2F *h2_hycal_mollerCenter = new TH2F("h2_hycal_mollerCenterX", "HyCal Moller Center X", 720, -360, 360, 720, -360, 360);
    TH2F *h2_gem0_mollerCenter  = new TH2F("h2_gem0_mollerCenter",   "GEM0 Moller Center",    720, -360, 360, 720, -360, 360);
    TH2F *h2_gem1_mollerCenter  = new TH2F("h2_gem1_mollerCenter",   "GEM1 Moller Center",    720, -360, 360, 720, -360, 360);
    TH2F *h2_gem2_mollerCenter  = new TH2F("h2_gem2_mollerCenter",   "GEM2 Moller Center",    720, -360, 360, 720, -360, 360);
    TH2F *h2_gem3_mollerCenter  = new TH2F("h2_gem3_mollerCenter",   "GEM3 Moller Center",    720, -360, 360, 720, -360, 360);

    TH1F *h1_hycal_mollerCenterX   = new TH1F("h1_hycal_mollerCenterX",   "HyCal Moller Center X",   400, -50, 50);
    TH1F *h1_hycal_mollerCenterY   = new TH1F("h1_hycal_mollerCenterY",   "HyCal Moller Center Y",   400, -50, 50);
    TH1F *h1_hycal_mollerZdistance = new TH1F("h1_hycal_mollerZdistance", "HyCal Moller Z distance", 4000, 0, 10000);

    TH1F *h1_gem0_mollerCenterX   = new TH1F("h1_gem0_mollerCenterX",   "GEM0 Moller Center X",   800, -50, 50);
    TH1F *h1_gem0_mollerCenterY   = new TH1F("h1_gem0_mollerCenterY",   "GEM0 Moller Center Y",   200, -50, 50);
    TH1F *h1_gem0_mollerZdistance = new TH1F("h1_gem0_mollerZdistance", "GEM0 Moller Z distance", 4000, 0, 10000);

    TH1F *h1_gem1_mollerCenterX   = new TH1F("h1_gem1_mollerCenterX",   "GEM1 Moller Center X",   800, -50, 50);
    TH1F *h1_gem1_mollerCenterY   = new TH1F("h1_gem1_mollerCenterY",   "GEM1 Moller Center Y",   200, -50, 50);
    TH1F *h1_gem1_mollerZdistance = new TH1F("h1_gem1_mollerZdistance", "GEM1 Moller Z distance", 4000, 0, 10000);

    TH1F *h1_gem2_mollerCenterX   = new TH1F("h1_gem2_mollerCenterX",   "GEM2 Moller Center X",   800, -50, 50);
    TH1F *h1_gem2_mollerCenterY   = new TH1F("h1_gem2_mollerCenterY",   "GEM2 Moller Center Y",   200, -50, 50);
    TH1F *h1_gem2_mollerZdistance = new TH1F("h1_gem2_mollerZdistance", "GEM2 Moller Z distance", 4000, 0, 10000);

    TH1F *h1_gem3_mollerCenterX   = new TH1F("h1_gem3_mollerCenterX",   "GEM3 Moller Center X",   800, -50, 50);
    TH1F *h1_gem3_mollerCenterY   = new TH1F("h1_gem3_mollerCenterY",   "GEM3 Moller Center Y",   200, -50, 50);
    TH1F *h1_gem3_mollerZdistance = new TH1F("h1_gem3_mollerZdistance", "GEM3 Moller Z distance", 4000, 0, 10000);

    TH1F *h1_phi_diff_hycal_gem0 = new TH1F("h1_phi_diff_hycal_gem0", "Phi Difference HyCal-GEM0", 80, -20, 20);
    TH1F *h1_phi_diff_hycal_gem1 = new TH1F("h1_phi_diff_hycal_gem1", "Phi Difference HyCal-GEM1", 80, -20, 20);
    TH1F *h1_phi_diff_hycal_gem2 = new TH1F("h1_phi_diff_hycal_gem2", "Phi Difference HyCal-GEM2", 80, -20, 20);
    TH1F *h1_phi_diff_hycal_gem3 = new TH1F("h1_phi_diff_hycal_gem3", "Phi Difference HyCal-GEM3", 80, -20, 20);

    auto fill_moller_data = [&](const TaggedMollerData &mollers, TH2F *hits,
                                TH2F *centers, TH1F *center_x,
                                TH1F *center_y, TH1F *z_distance) {
        const size_t n = mollers.size();
        for (size_t i = 0; i < n; ++i) {
            const TaggedMollerEvent &event = mollers[i];
            const MollerEvent current = ToMollerEvent(event);
            hits->Fill(event.first.x, event.first.y);
            hits->Fill(event.second.x, event.second.y);
            z_distance->Fill(PhysicsTools::GetMollerZdistance(current, gRunConfig.Ebeam));

            if (i == 0) continue;
            MollerEvent previous = ToMollerEvent(mollers[i - 1]);
            auto center = PhysicsTools::GetMollerCenter(previous, current);
            if(center[0] != 0 || center[1] != 0) {
                centers->Fill(center[0], center[1]);
                center_x->Fill(center[0]);
                center_y->Fill(center[1]);
            }

            if (i > 1) {
                previous = ToMollerEvent(mollers[i - 2]);
                center = PhysicsTools::GetMollerCenter(previous, current);
                if(center[0] != 0 || center[1] != 0) {
                    centers->Fill(center[0], center[1]);
                    center_x->Fill(center[0]);
                    center_y->Fill(center[1]);
                }
            }
        }
    };

    fill_moller_data(hycal_mollers, h2_hycal_hits, h2_hycal_mollerCenter,
                     h1_hycal_mollerCenterX, h1_hycal_mollerCenterY,
                     h1_hycal_mollerZdistance);

    std::array<TH2F *, 4> gem_hits = {
        h2_gem0_hits, h2_gem1_hits, h2_gem2_hits, h2_gem3_hits};
    std::array<TH2F *, 4> gem_centers = {
        h2_gem0_mollerCenter, h2_gem1_mollerCenter,
        h2_gem2_mollerCenter, h2_gem3_mollerCenter};
    std::array<TH1F *, 4> gem_center_x = {
        h1_gem0_mollerCenterX, h1_gem1_mollerCenterX,
        h1_gem2_mollerCenterX, h1_gem3_mollerCenterX};
    std::array<TH1F *, 4> gem_center_y = {
        h1_gem0_mollerCenterY, h1_gem1_mollerCenterY,
        h1_gem2_mollerCenterY, h1_gem3_mollerCenterY};
    std::array<TH1F *, 4> gem_z_distance = {
        h1_gem0_mollerZdistance, h1_gem1_mollerZdistance,
        h1_gem2_mollerZdistance, h1_gem3_mollerZdistance};
    std::array<TH1F *, 4> hycal_gem_phi = {
        h1_phi_diff_hycal_gem0, h1_phi_diff_hycal_gem1,
        h1_phi_diff_hycal_gem2, h1_phi_diff_hycal_gem3};

    auto track_phi_diff = [](const DataPoint &hc, const DataPoint &gem) {
        float diff = PhysicsTools::GetPhiAngle(gem.x, gem.y) - PhysicsTools::GetPhiAngle(hc.x, hc.y);
        return diff;
    };

    std::unordered_map<int, const TaggedMollerEvent *> hycal_by_event;
    hycal_by_event.reserve(hycal_mollers.size());
    for (const TaggedMollerEvent &hc : hycal_mollers)
        hycal_by_event.emplace(hc.event_num, &hc);

    for (size_t det = 0; det < gem_mollers.size(); ++det) {
        const TaggedMollerData &gm = gem_mollers[det];
        fill_moller_data(gm, gem_hits[det], gem_centers[det],
                         gem_center_x[det], gem_center_y[det],
                         gem_z_distance[det]);

        TH1F *phi_hist = hycal_gem_phi[det];
        for (const TaggedMollerEvent &gem : gm) {
            const auto hc_it = hycal_by_event.find(gem.event_num);
            if (hc_it == hycal_by_event.end()) continue;
            const TaggedMollerEvent &hc = *hc_it->second;
            phi_hist->Fill(track_phi_diff(hc.first, gem.first));
            phi_hist->Fill(track_phi_diff(hc.second, gem.second));
        }
    }

    fs::create_directories(output_dir);
    const std::string fit_run_tag = run_num >= 0 ? std::to_string(run_num) : "unknown";
    const fs::path fit_plot_dir =
        fs::path(output_dir) / ("det_calib_run" + fit_run_tag + "_fits");
    fs::create_directories(fit_plot_dir);

    struct PeakFitResult {
        std::string detector;
        std::string parameter;
        TH1F *hist = nullptr;
        double entries = 0.;
        double center = 0.;
        double fwhm = 0.;
        double sigma = 0.;
        double left = 0.;
        double right = 0.;
        double peak = 0.;
        std::string status = "empty";
        bool ok = false;
    };

    auto extract_peak = [](TH1F *hist, const std::string &detector,
                           const std::string &parameter) {
        PeakFitResult result;
        result.detector = detector;
        result.parameter = parameter;
        result.hist = hist;
        if (hist) result.entries = hist->GetEntries();
        if (!hist || hist->GetEntries() <= 0.) return result;

        const int nbins = hist->GetNbinsX();
        const int max_bin = hist->GetMaximumBin();
        const double peak = hist->GetBinContent(max_bin);
        result.center = hist->GetXaxis()->GetBinCenter(max_bin);
        result.peak = peak;
        if (nbins <= 0 || peak <= 0.) return result;

        const double half_peak = 0.5 * peak;
        int left_bin = max_bin;
        while (left_bin > 1 && hist->GetBinContent(left_bin) >= half_peak)
            --left_bin;

        int right_bin = max_bin;
        while (right_bin < nbins && hist->GetBinContent(right_bin) >= half_peak)
            ++right_bin;

        auto crossing = [&](int bin0, int bin1) {
            const double x0 = hist->GetXaxis()->GetBinCenter(bin0);
            const double x1 = hist->GetXaxis()->GetBinCenter(bin1);
            const double y0 = hist->GetBinContent(bin0);
            const double y1 = hist->GetBinContent(bin1);
            if (std::abs(y1 - y0) < 1e-12) return 0.5 * (x0 + x1);
            double frac = (half_peak - y0) / (y1 - y0);
            frac = std::max(0.0, std::min(1.0, frac));
            return x0 + frac * (x1 - x0);
        };

        const bool left_at_edge = left_bin == 1 &&
                                  hist->GetBinContent(left_bin) >= half_peak;
        const bool right_at_edge = right_bin == nbins &&
                                   hist->GetBinContent(right_bin) >= half_peak;
        result.left = left_at_edge
            ? hist->GetXaxis()->GetBinLowEdge(1)
            : crossing(left_bin, left_bin + 1);
        result.right = right_at_edge
            ? hist->GetXaxis()->GetBinUpEdge(nbins)
            : crossing(right_bin - 1, right_bin);
        result.fwhm = std::max(0.0, result.right - result.left);
        if (result.fwhm <= 0.)
            result.fwhm = hist->GetXaxis()->GetBinWidth(max_bin);
        result.sigma = result.fwhm / 2.354820045;

        double fit_lo = result.left;
        double fit_hi = result.right;
        const double bin_width = hist->GetXaxis()->GetBinWidth(max_bin);
        if (hist->GetXaxis()->FindBin(fit_hi) -
                hist->GetXaxis()->FindBin(fit_lo) + 1 < 3) {
            fit_lo = result.center - 2.0 * bin_width;
            fit_hi = result.center + 2.0 * bin_width;
        }
        fit_lo = std::max(fit_lo, hist->GetXaxis()->GetXmin());
        fit_hi = std::min(fit_hi, hist->GetXaxis()->GetXmax());

        TF1 fit_func((std::string(hist->GetName()) + "_gaus_fit").c_str(),
                     "gaus", fit_lo, fit_hi);
        fit_func.SetParameters(peak, result.center, std::max(result.sigma, bin_width));
        fit_func.SetParLimits(0, 0., std::max(peak * 10., 1.));
        fit_func.SetParLimits(1, fit_lo, fit_hi);
        fit_func.SetParLimits(2, 0.1 * bin_width, std::max(fit_hi - fit_lo, bin_width));

        const int fit_status = hist->Fit(&fit_func, "QRN");
        if (fit_status == 0 &&
            std::isfinite(fit_func.GetParameter(1)) &&
            std::isfinite(fit_func.GetParameter(2)) &&
            fit_func.GetParameter(2) > 0.) {
            result.peak = fit_func.GetParameter(0);
            result.center = fit_func.GetParameter(1);
            result.sigma = std::abs(fit_func.GetParameter(2));
            result.fwhm = 2.354820045 * result.sigma;
            result.left = result.center - 0.5 * result.fwhm;
            result.right = result.center + 0.5 * result.fwhm;
            result.status = "ok";
            result.ok = true;
        } else {
            result.status = "fit-failed";
        }
        return result;
    };

    std::vector<PeakFitResult> fit_results;
    fit_results.reserve(19);
    gROOT->SetBatch(kTRUE);
    TCanvas fit_canvas("c_det_calib_fit", "Detector calibration fit", 900, 650);

    auto fit_draw_save = [&](TH1F *hist, const std::string &detector,
                             const std::string &parameter) {
        PeakFitResult result = extract_peak(hist, detector, parameter);
        fit_results.push_back(result);

        fit_canvas.Clear();
        if (hist) hist->Draw();
        TLatex latex;
        latex.SetNDC(true);
        latex.SetTextSize(0.035);
        if (result.ok) {
            TF1 *peak_shape = new TF1(
                (std::string(hist->GetName()) + "_fwhm_peak").c_str(),
                "gaus", result.left, result.right);
            peak_shape->SetParameters(result.peak, result.center, result.sigma);
            peak_shape->SetLineColor(kRed + 1);
            peak_shape->SetLineWidth(2);
            hist->GetListOfFunctions()->Add(peak_shape);
            peak_shape->Draw("same");
            latex.DrawLatex(0.15, 0.84, Form("center = %.4g", result.center));
            latex.DrawLatex(0.15, 0.79, Form("FWHM = %.4g", result.fwhm));
        } else {
            latex.DrawLatex(0.15, 0.84, "No peak found");
        }
        if (hist) {
            const std::string png =
                (fit_plot_dir / (std::string(hist->GetName()) + ".png")).string();
            fit_canvas.SaveAs(png.c_str());
        }
        return result;
    };

    PeakFitResult hycal_fit_x =
        fit_draw_save(h1_hycal_mollerCenterX, "HyCal", "center_x_mm");
    PeakFitResult hycal_fit_y =
        fit_draw_save(h1_hycal_mollerCenterY, "HyCal", "center_y_mm");
    PeakFitResult hycal_fit_z =
        fit_draw_save(h1_hycal_mollerZdistance, "HyCal", "z_distance_mm");
    std::array<PeakFitResult, 4> gem_fit_x;
    std::array<PeakFitResult, 4> gem_fit_y;
    std::array<PeakFitResult, 4> gem_fit_z;
    std::array<PeakFitResult, 4> gem_fit_phi;
    for (size_t det = 0; det < gem_center_x.size(); ++det) {
        const std::string det_name = "GEM" + std::to_string(det);
        gem_fit_x[det] = fit_draw_save(gem_center_x[det], det_name, "center_x_mm");
        gem_fit_y[det] = fit_draw_save(gem_center_y[det], det_name, "center_y_mm");
        gem_fit_z[det] = fit_draw_save(gem_z_distance[det], det_name, "z_distance_mm");
        gem_fit_phi[det] =
            fit_draw_save(hycal_gem_phi[det], det_name, "phi_hycal_gem_deg");
    }

    const fs::path fit_log_path =
        fs::path(output_dir) / ("det_calib_run" + fit_run_tag + "_fit_summary.txt");
    std::ofstream fit_log(fit_log_path.string());
    auto write_fit_table = [&](std::ostream &os) {
        os << std::left
           << std::setw(10) << "Detector"
           << std::setw(20) << "Parameter"
           << std::right
           << std::setw(12) << "Entries"
           << std::setw(14) << "Center"
           << std::setw(14) << "FWHM"
           << std::setw(14) << "Sigma"
           << std::setw(12) << "Status" << '\n';
        os << std::fixed << std::setprecision(4);
        for (const PeakFitResult &result : fit_results) {
            os << std::left
               << std::setw(10) << result.detector
               << std::setw(20) << result.parameter
               << std::right
               << std::setw(12) << result.entries
               << std::setw(14) << result.center
               << std::setw(14) << result.fwhm
               << std::setw(14) << result.sigma
               << std::setw(12) << result.status << '\n';
        }
    };
    auto write_position_summary = [&](std::ostream &os) {
        os << "\nPeak-center detector calibration parameters\n";
        os << std::left
           << std::setw(10) << "Detector"
           << std::right
           << std::setw(14) << "x_mm"
           << std::setw(14) << "y_mm"
           << std::setw(14) << "z_mm"
           << std::setw(14) << "phi_deg" << '\n';
        os << std::fixed << std::setprecision(4);
        os << std::left << std::setw(10) << "HyCal" << std::right
           << std::setw(14) << hycal_fit_x.center
           << std::setw(14) << hycal_fit_y.center
           << std::setw(14) << hycal_fit_z.center
           << std::setw(14) << 0.0 << '\n';
        for (size_t det = 0; det < gem_fit_x.size(); ++det) {
            os << std::left << std::setw(10) << ("GEM" + std::to_string(det))
               << std::right
               << std::setw(14) << gem_fit_x[det].center
               << std::setw(14) << gem_fit_y[det].center
               << std::setw(14) << gem_fit_z[det].center
               << std::setw(14) << gem_fit_phi[det].center << '\n';
        }

        const double target_x = gRunConfig.target_x;
        const double target_y = gRunConfig.target_y;
        const double target_z = gRunConfig.target_z;
        const double hycal_json_x = gRunConfig.hycal_x + target_x;
        const double hycal_json_y = gRunConfig.hycal_y + target_y;
        const double hycal_dx = hycal_fit_x.ok ? hycal_fit_x.center : 0.;
        const double hycal_dy = hycal_fit_y.ok ? hycal_fit_y.center : 0.;
        const double corrected_target_x = target_x + hycal_dx;
        const double corrected_target_y = target_y + hycal_dy;
        const double corrected_target_z = 0.;
        const double corrected_hycal_z =
            hycal_fit_z.ok ? hycal_fit_z.center : gRunConfig.hycal_z + target_z;

        os << "\nCorrected runinfo geometry snippet\n";
        os << std::fixed << std::setprecision(6);
        os << "{\n";
        if (run_num >= 0)
            os << "    \"from_run\": " << run_num << ",\n";
        os << "    \"target\": ["
           << corrected_target_x << ", "
           << corrected_target_y << ", "
           << corrected_target_z << "],\n";
        os << "    \"hycal\": {\n";
        os << "        \"position\": ["
           << hycal_json_x << ", "
           << hycal_json_y << ", "
           << corrected_hycal_z << "], "
           << "\"tilting\": ["
           << gRunConfig.hycal_tilt_x << ", "
           << gRunConfig.hycal_tilt_y << ", "
           << gRunConfig.hycal_tilt_z << "]\n";
        os << "    },\n";
        os << "    \"gem\": {\n";
        os << "        \"detectors\": [\n";
        for (size_t det = 0; det < gem_fit_x.size(); ++det) {
            const double old_gem_x = gRunConfig.gem_x[det] + target_x;
            const double old_gem_y = gRunConfig.gem_y[det] + target_y;
            const double old_gem_z = gRunConfig.gem_z[det] + target_z;
            const double corrected_gem_x = gem_fit_x[det].ok
                ? old_gem_x - gem_fit_x[det].center
                : old_gem_x;
            const double corrected_gem_y = gem_fit_y[det].ok
                ? old_gem_y - gem_fit_y[det].center
                : old_gem_y;
            const double corrected_gem_z = gem_fit_z[det].ok
                ? gem_fit_z[det].center
                : old_gem_z;
            os << "            {\"id\": " << det
               << ", \"position\": ["
               << corrected_gem_x << ", "
               << corrected_gem_y << ", "
               << corrected_gem_z << "], \"tilting\": ["
               << gRunConfig.gem_tilt_x[det] << ", "
               << gRunConfig.gem_tilt_y[det] << ", "
               << gRunConfig.gem_tilt_z[det] << "]}";
            os << (det + 1 == gem_fit_x.size() ? "\n" : ",\n");
        }
        os << "        ]\n";
        os << "    }\n";
        os << "}\n";
    };

    std::cout << "Detector position calibration fit summary:\n";
    write_fit_table(std::cout);
    write_position_summary(std::cout);
    if (fit_log.is_open()) {
        fit_log << "Detector position calibration fit summary\n";
        write_fit_table(fit_log);
        write_position_summary(fit_log);
        std::cout << "Fit summary written to " << fit_log_path << "\n";
    } else {
        std::cerr << "Warning: cannot write fit summary to "
                  << fit_log_path << "\n";
    }

    fs::create_directories(output_dir);
    const std::string run_tag = run_num >= 0 ? std::to_string(run_num) : "unknown";
    const std::string output_root =
        (fs::path(output_dir) / ("det_calib_run" + run_tag + ".root")).string();
    TFile outfile(output_root.c_str(), "RECREATE");
    h2_hycal_hits->Write();
    h2_gem0_hits->Write();
    h2_gem1_hits->Write();
    h2_gem2_hits->Write();
    h2_gem3_hits->Write();
    h2_hycal_mollerCenter->Write();
    h2_gem0_mollerCenter->Write();
    h2_gem1_mollerCenter->Write();
    h2_gem2_mollerCenter->Write();
    h2_gem3_mollerCenter->Write();
    h1_hycal_mollerCenterX->Write();
    h1_hycal_mollerCenterY->Write();
    h1_hycal_mollerZdistance->Write();
    h1_gem0_mollerCenterX->Write();
    h1_gem0_mollerCenterY->Write();
    h1_gem0_mollerZdistance->Write();
    h1_gem1_mollerCenterX->Write();
    h1_gem1_mollerCenterY->Write();
    h1_gem1_mollerZdistance->Write();
    h1_gem2_mollerCenterX->Write();
    h1_gem2_mollerCenterY->Write();
    h1_gem2_mollerZdistance->Write();
    h1_gem3_mollerCenterX->Write();
    h1_gem3_mollerCenterY->Write();
    h1_gem3_mollerZdistance->Write();
    h1_phi_diff_hycal_gem0->Write();
    h1_phi_diff_hycal_gem1->Write();
    h1_phi_diff_hycal_gem2->Write();
    h1_phi_diff_hycal_gem3->Write();
    outfile.Close();
    std::cout << "Histograms written to " << output_root << "\n";

    return errors.load() == 0 ? 0 : 1;
}

static bool ProcessEVIO (const std::string &input_evio, RunConfig &gRunConfig,
                                const std::string &db_dir,
                                const std::string &run_config_file,
                                const std::string &daq_config_file,
                                const std::string &hycal_map_file,
                                const std::string &gem_ped_file,
                                const float zerosup_override,
                                const int max_events,
                                TaggedMollerData &hycal_mollers,
                                std::array<TaggedMollerData, 4> &gem_mollers)
{
    fdec::HyCalSystem                 hycal;
    gem::GemSystem                    gem_sys;
    fdec::ClusterConfig               cluster_cfg;
    prad2::HyCalTimeCuts              hc_time_cuts;
    DetectorTransform                 hycal_transform;
    std::array<DetectorTransform, 4>  gem_transforms;
    std::unordered_map<int, int>      roc_to_crate;

    prad2::Pipeline pipeline;
    try {
        pipeline = prad2::PipelineBuilder()
            .set_database_dir(db_dir)
            .set_daq_config(daq_config_file)
            .set_runinfo(run_config_file)
            .set_hycal_map(hycal_map_file)
            .set_gem_pedestal(gem_ped_file)
            .set_run_number_from_evio(input_evio)
            .set_log_stream(&std::cerr)
            .build();
    } catch (const std::exception &e) {
        std::cerr << "Replay: setup failed for " << input_evio
                  << ": " << e.what() << "\n";
        return false;
    }

    evc::DaqConfig daq_cfg = std::move(pipeline.daq_cfg);
    gRunConfig      = pipeline.run_cfg;
    hycal            = std::move(pipeline.hycal);
    gem_sys          = std::move(pipeline.gem);
    cluster_cfg      = pipeline.hycal_cluster_cfg;
    hc_time_cuts     = std::move(pipeline.hycal_time_cuts);
    hycal_transform  = pipeline.hycal_transform;
    gem_transforms   = pipeline.gem_transforms;

    // ROC→crate map from the same DAQ config the builder consumed.
    for (const auto &re : daq_cfg.roc_tags) {
        if (re.crate < 0) continue;
        if (!re.type.empty() && re.type != "roc" && re.type != "gem") continue;
        roc_to_crate[re.tag] = re.crate;
    }

    if (zerosup_override > 0.f)
        gem_sys.SetZeroSupThreshold(zerosup_override);

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching;
    matching.SetMatchRange(gRunConfig.matching_radius);
    matching.SetSquareSelection(gRunConfig.matching_use_square);
    matching.SetEnergyDependent(gRunConfig.matching_energy_dependent);
    matching.SetMatchSigma(gRunConfig.matching_sigma);
    PhysicsTools physics(hycal);
    //open EVIO file and output ROOT file
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    auto ev = std::make_unique<EventVars_Recon>();

    //initialize tools for event decoder and cluster reconstruction
    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
    fdec::PulseTemplateStore template_store;
    if (daq_cfg.wave_cfg.nnls_deconv.enabled
        && !daq_cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        template_store.LoadFromFile(
            db_dir + "/" + daq_cfg.wave_cfg.nnls_deconv.template_file,
            daq_cfg.wave_cfg);
    }
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;

    int total = 0;
    int processed_events = 0;

    int run_num = get_run_int(input_evio);
    std::string gain_data_dir = gRunConfig.gain_data_dir;
    if (gain_data_dir.empty())
        gain_data_dir = db_dir + "/gain_factor";
    else if (!fs::path(gain_data_dir).is_absolute())
        gain_data_dir = db_dir + "/" + gain_data_dir;
    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(
        gain_data_dir + "/gain_correction", run_num);

    // Per-detector lab transforms — set up by either branch of the detector
    // wiring above (PipelineBuilder for PRad-II, BuildLabTransforms for PRad-1).
    const auto &hc_xform = hycal_transform;
    const auto &g_xform  = gem_transforms;

    while ((max_events <= 0 || processed_events < max_events)
           && ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            if (max_events > 0 && processed_events >= max_events) break;
            event->clear();
            ssp_evt->clear();
            clusterer.Clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;
            ++processed_events;

            *ev = EventVars_Recon{};
            ev->event_num    = event->info.event_number;
            ev->trigger_bits = event->info.trigger_bits;

            bool is_sum      = (ev->trigger_bits & prad2::TBIT_sum)   != 0;
            if (!is_sum) continue;

            // Per-event gain correction (time-series lookup by event number).
            const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

            // decode FADC250 and reconstruct HyCal data
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 64; ++c) { //should be 16, a bigger number to adapt PRad1 data
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0) continue;

                        const auto *mod = hycal.module_by_daq(crate, s, c);
                        if (!mod || !mod->is_hycal()) continue;
                        // Per-ID gain correction: average of LMS 2/3 channels.
                        const float gain = (mod->id > 1000)
                            ? (gain_corr.w[mod->id - 1000].corr[1] + gain_corr.w[mod->id - 1000].corr[2]) / 2.0f
                            : gain_corr.g[mod->id].avg;

                        ana.SetChannelKey(roc.tag, s, c);
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        if (wres.npeaks <= 0) continue;

                        const auto hc_win = hc_time_cuts.at(mod->index);
                        if (cluster_cfg.seed_time_window > 0.f) {
                            // Multi-pulse mode: push every peak inside the trigger
                            // window into the clusterer; the seed-anchored timing
                            // coincidence cut is applied inside HyCalCluster.
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                const auto &pk = wres.peaks[p];
                                if (pk.time <= hc_win.lo) continue;
                                if (pk.time >= hc_win.hi) continue;
                                float adc = pk.integral * gain;
                                float energy = static_cast<float>(mod->energize(adc));
                                clusterer.AddHit(mod->index, energy, pk.time);
                                ev->total_energy += energy;
                                nch++;
                            }
                        } else {
                            // Legacy: pick the largest in-window peak as the single
                            // module hit, time field unused downstream.
                            int bestIdx = -1;
                            float bestHeight = -1.f;
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                const auto &pk = wres.peaks[p];
                                if (pk.time > hc_win.lo &&
                                    pk.time < hc_win.hi &&
                                    pk.height > bestHeight) {
                                    bestHeight = pk.height;
                                    bestIdx = p;
                                }
                            }
                            if (bestIdx < 0) continue;
                            float adc = wres.peaks[bestIdx].integral * gain;
                            float energy = static_cast<float>(mod->energize(adc));
                            clusterer.AddHit(mod->index, energy, wres.peaks[bestIdx].time);
                            ev->total_energy += energy;
                            nch++;
                        }
                    }
                }
            }
            if(nch > 200) continue; // too many hits, likely noise, skip the event

            clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hits;
            clusterer.ReconstructHits(hits);
            //HyCal event reconstrued
            ev->n_clusters = std::min((int)hits.size(), prad2::kMaxClusters);
            for (int i = 0; i < ev->n_clusters; ++i) {
                ev->cl_nblocks[i] = hits[i].nblocks;
                ev->cl_time[i]    = hits[i].time;
                //transform the cluster positions to the lab coordinate
                HCHit local_hit = {hits[i].x, hits[i].y, fdec::shower_depth(hits[i].center_id, hits[i].energy),
                    hits[i].energy, static_cast<uint16_t>(hits[i].center_id), hits[i].flag};
                analysis::ApplyToLab(hc_xform, local_hit);
                GetProjection(local_hit, gRunConfig.hycal_z);
                ev->cl_x[i] = local_hit.x;
                ev->cl_y[i] = local_hit.y;
                ev->cl_z[i] = local_hit.z;
                ev->cl_energy[i] = local_hit.energy;
                ev->cl_linear_corr[i] = hits[i].linear_corr;
                ev->cl_center[i] = local_hit.center_id;
                ev->cl_flag[i] = local_hit.flag;
            }

            //decode GEM data and reconstruct GEM hits
            if(gem_sys.GetNDetectors() > 0){
                gem_sys.Clear();
                gem_sys.ProcessEvent(*ssp_evt);
                gem_sys.Reconstruct(gem_clusterer);
                auto &all_hits = gem_sys.GetAllHits();
                ev->n_gem_hits = 0;
                for (const auto &h : all_hits) {
                    if (ev->n_gem_hits >= prad2::kMaxGemHits) break;
                    if (h.det_id < 0 || h.det_id >= 4) continue;
                    const int i = ev->n_gem_hits++;
                    ev->det_id[i] = h.det_id;
                    ev->gem_x_charge[i] = h.x_charge;
                    ev->gem_y_charge[i] = h.y_charge;
                    ev->gem_x_peak[i] = h.x_peak;
                    ev->gem_y_peak[i] = h.y_peak;
                    ev->gem_x_size[i] = h.x_size;
                    ev->gem_y_size[i] = h.y_size;
                    ev->gem_x_mTbin[i] = h.x_max_timebin;
                    ev->gem_y_mTbin[i] = h.y_max_timebin;
                    //transform the GEM hit positions to the lab coordinate
                    GEMHit local_hit = {h.x, h.y, 0.f, static_cast<uint8_t>(h.det_id)};
                    int d = local_hit.det_id;
                    if (d >= 0 && d < 4) {
                        analysis::ApplyToLab(g_xform[d], local_hit);
                    }
                    ev->gem_x[i] = local_hit.x;
                    ev->gem_y[i] = local_hit.y;
                    ev->gem_z[i] = local_hit.z;
                }

                // Perform matching between HyCal clusters and GEM hits
                //store all the hits on HyCal and GEMs in this event
                std::vector<HCHit> hc_hits;
                std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
                for (int i = 0; i < ev->n_clusters; ++i)
                    hc_hits.push_back({ev->cl_x[i], ev->cl_y[i], ev->cl_z[i], ev->cl_energy[i], ev->cl_center[i], ev->cl_flag[i]});
                for (int i = 0; i < ev->n_gem_hits; ++i) {
                    const int det_id = ev->det_id[i];
                    if (det_id < 0 || det_id >= 4) continue;
                    gem_hits[det_id].push_back(GEMHit{
                        ev->gem_x[i], ev->gem_y[i], ev->gem_z[i],
                        static_cast<uint8_t>(det_id)});
                }
                
                // already transform to the coordinates

                std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
                std::vector<MatchHit_perChamber> matched_hits_chamber = matching.MatchPerChamber(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]); 
                
                for (const auto &m : matched_hits_chamber) {
                    const int cl_idx = m.hycal_idx;
                    if (cl_idx < 0 || cl_idx >= prad2::kMaxClusters) continue;
                    for(int j = 0; j < 4; j++){
                        ev->matchGEMx[cl_idx][j] = m.gem_hits[j][0];
                        ev->matchGEMy[cl_idx][j] = m.gem_hits[j][1];
                        ev->matchGEMz[cl_idx][j] = m.gem_hits[j][2];
                    }
                    ev->matchFlag[cl_idx] = m.mflag;
                }

                ev->matchNum = std::min((int)matched_hits.size(), prad2::kMaxClusters);
                for (int i = 0; i < ev->matchNum; i++){
                    // save the matched GEM hit (must 2 matchings) info in mHit_ arrays for quick check
                    ev->mHit_E[i] = matched_hits[i].hycal_hit.energy;
                    ev->mHit_x[i] = matched_hits[i].hycal_hit.x;
                    ev->mHit_y[i] = matched_hits[i].hycal_hit.y;
                    ev->mHit_z[i] = matched_hits[i].hycal_hit.z;
                    for(int j = 0; j < 2; j++) {
                        ev->mHit_gx[i][j] =  matched_hits[i].gem[j].x;
                        ev->mHit_gy[i][j] =  matched_hits[i].gem[j].y;
                        ev->mHit_gz[i][j] =  matched_hits[i].gem[j].z;
                        ev->mHit_gid[i][j] = matched_hits[i].gem[j].det_id; // placeholder for GEM hit ID if needed
                    }
                }
            }

            // select events and analyze
            if (ev->n_clusters != 2) continue;
            if (ev->matchNum != 2) continue;
            if (ev->cl_nblocks[0] < 3 || ev->cl_nblocks[1] < 3) continue;
            if (std::fabs(ev->cl_x[0]) < 20.75f * 2.5f && std::fabs(ev->cl_y[0]) < 20.75f * 2.5f) continue;
            if (std::fabs(ev->cl_x[1]) < 20.75f * 2.5f && std::fabs(ev->cl_y[1]) < 20.75f * 2.5f) continue;
            if (std::fabs(ev->cl_x[0]) > 20.75f * 16.f || std::fabs(ev->cl_y[0]) > 20.75f * 16.f) continue;
            if (std::fabs(ev->cl_x[1]) > 20.75f * 16.f || std::fabs(ev->cl_y[1]) > 20.75f * 16.f) continue;

            MollerEvent m_hc(
                {ev->cl_x[0], ev->cl_y[0], ev->cl_z[0], ev->cl_energy[0]},
                {ev->cl_x[1], ev->cl_y[1], ev->cl_z[1], ev->cl_energy[1]});

            constexpr float kRadToDeg = 57.29577951308232f;
            float theta1 = std::atan2(std::hypot(m_hc.first.x, m_hc.first.y), m_hc.first.z) * kRadToDeg;
            float theta2 = std::atan2(std::hypot(m_hc.second.x, m_hc.second.y), m_hc.second.z) * kRadToDeg;
            if (!physics.isMoller_kinematic(theta1, m_hc.first.E,
                                            theta2, m_hc.second.E,
                                            gRunConfig.Ebeam, 0.035f))
                continue;
            if (std::fabs(physics.GetMollerPhiDiff(m_hc)) > 6.f) continue;

            //add some scattering angle cuts for 0.7GeV
            if(gRunConfig.Ebeam > 0.f && gRunConfig.Ebeam < 1000.f) {
                if(theta1 < 1.5f || theta2 < 1.5f) continue;
            }

            hycal_mollers.emplace_back(m_hc.first, m_hc.second, ev->event_num);

            for(int did = 0; did < 4; did ++){
                if (((ev->matchFlag[0] & (1u << did)) != 0)
                    && ((ev->matchFlag[1] & (1u << did)) != 0)) {
                    gem_mollers[did].emplace_back(
                        DataPoint(ev->matchGEMx[0][did], ev->matchGEMy[0][did],
                                  ev->matchGEMz[0][did], ev->cl_energy[0]),
                        DataPoint(ev->matchGEMx[1][did], ev->matchGEMy[1][did],
                                  ev->matchGEMz[1][did], ev->cl_energy[1]),
                        ev->event_num);
                }
            }

            total++;
            if (total % 1000 == 0)
                std::cerr << "\rFind: " << total << " moller events on HyCal" << std::flush;
        }
    }
    std::cerr << "\rReplay: " << total << " moller events reconstructed on HyCal\n";

    return true;
}
