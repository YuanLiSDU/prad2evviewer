//=============================================================================
// replay_gainCorr — extract LMS/alpha events from EVIO, then compute gain corrections
//
// Phase 1: replay EVIO files (multi-threaded) → per-file *_lms.root
// Phase 2: chain all *_lms.root, accumulate histograms in windows of N LMS
//          events, fit them, write gain-correction rows to gain_corr.root.
//          Each row includes the batch-midpoint unix_time derived from the
//          EPICS tree (with the scaler tree as fallback).
//
// Usage:
//   replay_gainCorr <evio_file_or_dir> [more files/dirs...]
//              -o output_dir
//              [-f max_files] [-j num_threads]
//              [-c daq_config.json] [-d hycal_map.json]
//              [-b batch_size] [-r ref_run]
//              [-s]  save intermediate *_lms.root files (hadd into one)
//   -o  output directory (REQUIRED)
//   -f  max files to process (default: all)
//   -j  number of threads (default: 4)
//   -c  DAQ configuration file
//   -d  HyCal map file (default: <db>/hycal_map.json)
//   -b  LMS events per gain-correction batch (default: 4000)
//   -r  reference run number for gain table (default: from general.json)
//   -s  save intermediate *_lms.root (merged via hadd); default: delete them
//=============================================================================

#include "Replay.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "GainCorrCompute.h"

#include <TClass.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;

// ── file discovery ───────────────────────────────────────────────────────────
static std::vector<std::string> collectEvioFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (std::filesystem::is_directory(path)) {
        for (auto &e : std::filesystem::directory_iterator(path)) {
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

// Pick up to max evenly-spaced elements from a vector.
static std::vector<TH1F*> subsample(const std::vector<TH1F*> &all, int max)
{
    if ((int)all.size() <= max) return all;
    std::vector<TH1F*> out;
    out.reserve(max);
    for (int i = 0; i < max; ++i) {
        int idx = (int)((i + 0.5) * (double)all.size() / max);
        out.push_back(all[idx]);
    }
    return out;
}

static const int kPlotColors[10] = {
    kBlack, kRed, kBlue, kGreen+2, kMagenta+1,
    kOrange+7, kCyan+1, kViolet+1, kTeal+1, kPink+1
};

// Draw overlaid histograms onto the current pad.
static void drawOverlay(const std::vector<TH1F*> &all, int max_hists, const char *title)
{
    if (all.empty()) return;
    auto sel = subsample(all, max_hists);

    // Determine x range from peak position and FWHM across all selected histograms.
    // Display window: [peak - 4*FWHM, peak + 4*FWHM], unioned across histograms.
    double ymax    = 0.;
    double xlo_acc =  1.e15;
    double xhi_acc = -1.e15;
    for (auto *h : sel) {
        ymax = std::max(ymax, h->GetMaximum());
        int    nb       = h->GetNbinsX();
        int    peak_bin = h->GetMaximumBin();
        double peak_val = h->GetBinContent(peak_bin);
        if (peak_val <= 0.) continue;

        // Scan outward from peak to estimate FWHM.
        double half_max = peak_val * 0.5;
        int lo_hm = 1;
        for (int b = peak_bin; b >= 1; --b)
            if (h->GetBinContent(b) < half_max) { lo_hm = b; break; }
        int hi_hm = nb;
        for (int b = peak_bin; b <= nb; ++b)
            if (h->GetBinContent(b) < half_max) { hi_hm = b; break; }
        double fwhm = h->GetBinCenter(hi_hm) - h->GetBinCenter(lo_hm);
        if (fwhm < h->GetBinWidth(1)) fwhm = h->GetBinWidth(1) * 10.;

        double center = h->GetBinCenter(peak_bin);
        xlo_acc = std::min(xlo_acc, center - 4. * fwhm);
        xhi_acc = std::max(xhi_acc, center + 4. * fwhm);
    }
    double xlo = 0., xhi = sel[0]->GetXaxis()->GetXmax();
    if (xlo_acc < xhi_acc) {
        xlo = std::max(0., xlo_acc);
        xhi = std::min(xhi, xhi_acc);
    }

    for (int k = 0; k < (int)sel.size(); ++k) {
        sel[k]->SetLineColor(kPlotColors[k % 10]);
        sel[k]->SetLineWidth(2);
        if (k == 0) {
            sel[k]->SetTitle(title);
            sel[k]->GetXaxis()->SetTitle("ADC counts");
            sel[k]->GetXaxis()->SetRangeUser(xlo, xhi);
            sel[k]->GetYaxis()->SetTitle("Entries");
            sel[k]->GetYaxis()->SetRangeUser(0., ymax * 1.15);
            sel[k]->Draw("HIST");
        } else {
            sel[k]->Draw("HIST SAME");
        }
    }
}

static void savePlots(const std::string &pdf_path,
                      GainPlotStore     &ps,
                      const GainPlotConfig &cfg)
{
    gROOT->SetBatch(kTRUE);
    TCanvas c("c_lmsplots", "LMS plots", 900, 700);
    c.Print((pdf_path + "[").c_str());   // open PDF

    // ── LMS reference PMT pages (one page per LMS channel) ──────────────
    for (int i = 0; i < kGainNLMS; ++i) {
        c.Clear();
        c.Divide(1, 2);

        c.cd(1);
        gPad->SetLeftMargin(0.12);
        drawOverlay(ps.ref_lms[i],   cfg.max_hists,
                    Form("LMS%d  LMS signal", i + 1));

        c.cd(2);
        gPad->SetLeftMargin(0.12);
        drawOverlay(ps.ref_alpha[i], cfg.max_hists,
                    Form("LMS%d  alpha signal", i + 1));

        c.Print(pdf_path.c_str());
    }

    // ── W module pages ────────────────────────────────────────────────────
    for (int wid : cfg.w_ids) {
        auto it = ps.mod_w.find(wid);
        if (it == ps.mod_w.end() || it->second.empty()) continue;
        c.Clear();
        gPad->SetLeftMargin(0.12);
        drawOverlay(it->second, cfg.max_hists,
                    Form("W%d  LMS signal", wid));
        c.Print(pdf_path.c_str());
    }

    c.Print((pdf_path + "]").c_str());   // close PDF
    std::cout << "  [plot] saved to " << pdf_path << "\n";
}

// ── main ─────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string daq_config, daq_map, output_dir;
    int  max_files   = -1;
    int  num_threads = 4;
    int  batch_size  = 4000;
    int  ref_run     = -1;
    bool save_lms    = false;
    GainPlotConfig plot_cfg;
    GainPlotStore  ps;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    daq_config = db_dir + "/daq_config.json";

    int opt;
    while ((opt = getopt(argc, argv, "o:f:c:d:j:b:r:pw:s")) != -1) {
        switch (opt) {
            case 'o': output_dir       = optarg; break;
            case 'f': max_files        = std::atoi(optarg); break;
            case 'c': daq_config       = optarg; break;
            case 'd': daq_map          = optarg; break;
            case 'j': num_threads      = std::atoi(optarg); break;
            case 'b': batch_size       = std::atoi(optarg); break;
            case 'r': ref_run          = std::atoi(optarg); break;
            case 'p': plot_cfg.enabled = true; break;
            case 's': save_lms         = true;  break;
            case 'w': {
                // comma-separated W module IDs, e.g. -w 565,892
                std::istringstream ss(optarg);
                std::string tok;
                while (std::getline(ss, tok, ',')) {
                    int id = std::atoi(tok.c_str());
                    if (id >= 1 && id <= kGainNW) plot_cfg.w_ids.insert(id);
                }
                break;
            }
        }
    }

    std::vector<std::string> evio_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectEvioFiles(argv[i]);
        evio_files.insert(evio_files.end(), f.begin(), f.end());
    }

    if (evio_files.empty() || output_dir.empty()) {
        std::cerr <<
            "Usage: replay_gainCorr <evio_file_or_dir> [...] -o output_dir\n"
            "       [-f max_files] [-j threads] [-c daq_config.json] [-d hycal_map.json]\n"
            "       [-b batch_size (4000)] [-r ref_run]\n"
            "       [-s] (save intermediate *_lms.root merged via hadd; default: delete)\n"
            "       [-p] [-w id1,id2,...] (enable plots; -w selects extra W modules)\n";
        return 1;
    }

    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    if (daq_map.empty()) daq_map = db_dir + "/hycal_map.json";

    int run_num = get_run_int(evio_files[0]);
    gRunConfig  = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    if (ref_run < 0) ref_run = gRunConfig.gain_ref_run;

    std::cout << "=== Phase 1: EVIO replay ===\n"
              << "  Files   : " << num_files << "\n"
              << "  Threads : " << num_threads << "\n"
              << "  Output  : " << output_dir << "\n"
              << "  DAQ cfg : " << daq_config << "\n"
              << "  HyCal   : " << daq_map << "\n";

    // ── Phase 1: multi-threaded EVIO → _lms.root ────────────────────────────
    std::atomic<int>         next_file{0};
    std::mutex               io_mtx;
    std::atomic<int>         replay_errors{0};
    std::vector<std::string> lms_out_files;
    std::mutex               out_files_mtx;

    auto worker = [&]() {
        analysis::Replay replay;
        if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
        replay.LoadHyCalMap(daq_map);

        while (true) {
            int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            std::string out = output_dir + "/" + MakeLMSOutputFile(evio_files[idx]);
            bool ok = replay.Process_LMSgainFactor(evio_files[idx], out,
                                                   db_dir, daq_config);
            {
                std::lock_guard<std::mutex> lk(io_mtx);
                if (ok) {
                    std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                              << evio_files[idx] << " -> " << out << "\n";
                } else {
                    std::cerr << "  [" << (idx + 1) << "/" << num_files << "] FAILED: "
                              << evio_files[idx] << "\n";
                    ++replay_errors;
                }
            }
            if (ok) {
                std::lock_guard<std::mutex> lk(out_files_mtx);
                lms_out_files.push_back(out);
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

    if (lms_out_files.empty()) {
        std::cerr << "No LMS files produced; skipping gain calculation.\n";
        return replay_errors.load() > 0 ? 1 : 0;
    }

    // ── Phase 2: gain corrections ────────────────────────────────────────────
    std::sort(lms_out_files.begin(), lms_out_files.end());

    auto ref_tbl = prad2::LoadRefGain(gRunConfig.gain_data_dir + "/ref_gain", ref_run);
    if (!ref_tbl.loaded) {
        std::cerr << "Warning: reference gain table not loaded"
                  << " (ref_run=" << ref_run
                  << ", dir=" << gRunConfig.gain_data_dir << ")\n";
    }

    // Save gain-correction output to the project database directory.
    std::string gain_corr_dir = db_dir + "/gain_factor/gain_correction";
    {
        std::error_code ec;
        std::filesystem::create_directories(gain_corr_dir, ec);
        if (ec)
            std::cerr << "Warning: cannot create " << gain_corr_dir
                      << ": " << ec.message() << "\n";
    }
    std::string gain_out = gain_corr_dir + "/" +
        std::string(Form("prad_%06d_gain_corr.root", run_num));

    std::cout << "\n=== Phase 2: gain corrections ===\n"
              << "  Batch size : " << batch_size << " LMS events\n"
              << "  Ref run    : " << ref_run << "\n"
              << "  Output     : " << gain_out << "\n";

    ComputeGainCorrections(lms_out_files, gain_out, batch_size, ref_run, ref_tbl, &plot_cfg, &ps);

    if (plot_cfg.enabled) {
        std::string pdf_out = gain_out.substr(0, gain_out.rfind('.')) + "_plots.pdf";
        std::cout << "\n=== Phase 3: saving plots ===\n"
                  << "  Output : " << pdf_out << "\n";
        if (!plot_cfg.w_ids.empty()) {
            std::cout << "  W mods :";
            for (int id : plot_cfg.w_ids) std::cout << " W" << id;
            std::cout << "\n";
        }
        savePlots(pdf_out, ps, plot_cfg);
    }

    // ── Clean up / merge intermediate LMS root files ────────────────────────
    if (!save_lms) {
        std::cout << "\n=== Cleaning up intermediate LMS root files ===\n";
        int n_del = 0;
        for (const auto &f : lms_out_files) {
            std::error_code ec;
            if (std::filesystem::remove(f, ec))
                ++n_del;
            else
                std::cerr << "  Warning: cannot delete " << f << ": " << ec.message() << "\n";
        }
        std::cout << "  Deleted " << n_del << " file(s).\n";
    } else {
        std::string merged = output_dir + "/" +
            std::string(Form("prad_%06d_lms.root", run_num));
        std::cout << "\n=== Merging intermediate LMS root files (hadd) ===\n"
                  << "  Output : " << merged << "\n";

        // Build hadd command: hadd -f <merged> <file1> <file2> ...
        std::string cmd = "hadd -f " + merged;
        for (const auto &f : lms_out_files)
            cmd += " " + f;
        int rc = std::system(cmd.c_str());
        if (rc != 0) {
            std::cerr << "  Warning: hadd returned " << rc
                      << "; individual files kept in " << output_dir << "\n";
        } else {
            // Remove individual files after successful merge.
            int n_del = 0;
            for (const auto &f : lms_out_files) {
                std::error_code ec;
                if (std::filesystem::remove(f, ec)) ++n_del;
            }
            std::cout << "  Merged and removed " << n_del << " individual file(s).\n";
        }
    }

    std::cout << "All done.\n";
    return replay_errors.load() > 0 ? 1 : 0;
}
