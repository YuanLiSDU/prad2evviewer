//=============================================================================
// hycal_rf_offset_calib — fit per-module HyCal→RF time offsets.
//
// Reads recon ROOT files produced by prad2ana_replay_recon (which now carry
// `cl_dt_rf`, `cl_center`, `cl_time`, and the RF arrays) and emits a
// database/hycal_rf_offsets/<run>.json that the next recon pass picks up.
//
// Folded `cl_dt_rf` already lives on (−T_RF/2, T_RF/2] = (−2.004, 2.004) ns.
// On a first pass with no offsets loaded, this distribution carries both
// (a) the per-crystal timing resolution σ_t and (b) the per-crystal
// constant offset.  The peak fit recovers (b); the residual width is (a).
//
// Cuts to keep only events with a trustworthy `cl_time` reference:
//   * n_clusters == 1                     (no pile-up ambiguity)
//   * cl_energy  > min_cluster_energy_MeV (waveform SNR → reliable timing)
//   * rf_n_a > 0                          (need an RF tick to fold against)
//   * cl_dt_rf == cl_dt_rf                (drop NaN)
//
// Per-module fit: Gaussian + flat background in a ±fit_window_ns range
// around the bin maximum.  A module passes when entries ≥ min_entries AND
// the fit converges; otherwise the JSON gets `offset_ns: 0.0` and the
// sidecar CSV flags it as "uncalibrated".
//
// Usage:
//   prad2ana_hycal_rf_offset_calib                                   \
//       -i prad_024840.*_recon.root                                  \
//       -o database/hycal_rf_offsets/24840.json                      \
//       [--min-energy 200]   [--min-entries 100]                     \
//       [--fit-window 1.5]   [--canvas qc_24840.pdf]                 \
//       [--db database/hycal_map.json]
//=============================================================================

#include "EventData.h"
#include "EventData_io.h"
#include "HyCalSystem.h"
#include "RfTime.h"
#include "InstallPaths.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using nlohmann::json;

namespace {

struct Config {
    std::vector<std::string> inputs;
    std::string              output_json;
    std::string              hycal_map = "hycal_map.json";
    std::string              canvas_pdf;
    std::string              csv_out;
    float                    min_cluster_energy = 200.f;   // MeV
    int                      min_entries        = 100;
    float                    fit_window_ns      = 1.5f;
    bool                     require_single_cluster = true;
};

void usage(const char *prog)
{
    std::cerr << "Usage: " << prog << " -i <recon.root>... -o <out.json>\n"
              << "  -i  one or more recon ROOT files (repeatable, also accepts "
                 "comma-separated)\n"
              << "  -o  output JSON path (REQUIRED)\n"
              << "  --hycal-map        HyCal map JSON (default: search database)\n"
              << "  --min-energy       min cluster energy MeV (default 200)\n"
              << "  --min-entries      min hist entries to fit (default 100)\n"
              << "  --fit-window       ± fit range around peak in ns "
                 "(default 1.5)\n"
              << "  --allow-multi      allow events with n_clusters > 1\n"
              << "  --canvas           write per-module QC plots to a PDF\n"
              << "  --csv              write per-module summary CSV\n";
}

bool parse_args(int argc, char **argv, Config &cfg)
{
    static option long_opts[] = {
        {"hycal-map",   required_argument, 0, 1},
        {"min-energy",  required_argument, 0, 2},
        {"min-entries", required_argument, 0, 3},
        {"fit-window",  required_argument, 0, 4},
        {"allow-multi", no_argument,       0, 5},
        {"canvas",      required_argument, 0, 6},
        {"csv",         required_argument, 0, 7},
        {0,0,0,0}};
    int c;
    while ((c = getopt_long(argc, argv, "i:o:h", long_opts, nullptr)) != -1) {
        switch (c) {
        case 'i': {
            std::string s = optarg;
            size_t pos = 0, comma;
            while ((comma = s.find(',', pos)) != std::string::npos) {
                cfg.inputs.emplace_back(s.substr(pos, comma - pos));
                pos = comma + 1;
            }
            cfg.inputs.emplace_back(s.substr(pos));
            break;
        }
        case 'o': cfg.output_json = optarg; break;
        case 1:   cfg.hycal_map   = optarg; break;
        case 2:   cfg.min_cluster_energy = std::stof(optarg); break;
        case 3:   cfg.min_entries  = std::stoi(optarg); break;
        case 4:   cfg.fit_window_ns = std::stof(optarg); break;
        case 5:   cfg.require_single_cluster = false; break;
        case 6:   cfg.canvas_pdf   = optarg; break;
        case 7:   cfg.csv_out      = optarg; break;
        case 'h':
        default:  usage(argv[0]); return false;
        }
    }
    if (cfg.inputs.empty() || cfg.output_json.empty()) {
        usage(argv[0]);
        return false;
    }
    return true;
}

// Fit one Gaussian + constant around the histogram max.  Returns true
// when the fit converges and mean lies inside the (−T_RF/2, T_RF/2] range.
// Populates `mean`, `sigma`, `chi2_ndf`.
bool fit_peak(TH1 *h, float fit_window, float &mean, float &sigma,
              float &chi2_ndf, std::string &status)
{
    if (h->GetEntries() < 1) { status = "no-entries"; return false; }
    const int max_bin = h->GetMaximumBin();
    const float center = h->GetXaxis()->GetBinCenter(max_bin);

    const float lo = center - fit_window;
    const float hi = center + fit_window;

    // Gaussian + constant background.  Starting values from the histogram
    // moments inside the fit window so unconverged fits are rare.
    TF1 fit("rf_peak", "gaus(0)+[3]", lo, hi);
    fit.SetParameter(0, h->GetBinContent(max_bin));
    fit.SetParameter(1, center);
    fit.SetParameter(2, 0.5f);            // expect σ ~ T_RF/8 .. T_RF/4
    fit.SetParameter(3, h->GetMinimum());
    fit.SetParLimits(1, lo, hi);
    fit.SetParLimits(2, 0.05f, prad2::RF_PERIOD_NS / 2.f);
    fit.SetParLimits(3, 0.f, std::max(1.0, h->GetMaximum()));

    const int rc = h->Fit(&fit, "QRNL");  // Q=quiet R=range N=no-draw L=likelihood
    if (rc != 0) { status = "fit-failed"; return false; }

    mean     = fit.GetParameter(1);
    sigma    = std::abs(fit.GetParameter(2));
    chi2_ndf = (fit.GetNDF() > 0) ? fit.GetChisquare() / fit.GetNDF() : -1.f;
    if (!std::isfinite(mean) || !std::isfinite(sigma)) {
        status = "non-finite";
        return false;
    }
    // Keep the offset within one RF period; anything else is a fit
    // wandering into the noise floor.
    if (std::abs(mean) > prad2::RF_PERIOD_NS / 2.f) {
        status = "out-of-range";
        return false;
    }
    status = "ok";
    return true;
}

} // namespace

int main(int argc, char **argv)
{
    Config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;

    // Load HyCal map so we have module names + types for the per-module loop.
    fdec::HyCalSystem hycal;
    std::string map_path = cfg.hycal_map;
    if (map_path.find('/') == std::string::npos) {
        std::string db_dir = prad2::resolve_data_dir(
            "PRAD2_DATABASE_DIR",
            {"../share/prad2evviewer/database"},
            DATABASE_DIR);
        if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;
        map_path = db_dir + "/" + map_path;
    }
    if (!hycal.Init(map_path)) {
        std::cerr << "Cannot load HyCal map from " << map_path << "\n";
        return 1;
    }
    std::cerr << "[setup] HyCal map : " << map_path
              << " (" << hycal.module_count() << " modules)\n";

    // Chain the inputs.
    TChain chain("recon");
    for (const auto &f : cfg.inputs) {
        const int n = chain.Add(f.c_str());
        std::cerr << "[setup] added " << f << "  files=" << n << "\n";
    }
    if (chain.GetEntries() == 0) {
        std::cerr << "No entries in the input chain.\n";
        return 1;
    }
    std::cerr << "[setup] total entries: " << chain.GetEntries() << "\n";

    // Bind only the branches we need (cheap when reading 150k+ events).
    prad2::ReconEventData ev;
    chain.SetBranchStatus("*", 0);
    for (auto *b : {"n_clusters", "cl_energy", "cl_center", "cl_time",
                    "cl_dt_rf", "rf_n_a"}) {
        chain.SetBranchStatus(b, 1);
    }
    chain.SetBranchAddress("n_clusters", &ev.n_clusters);
    chain.SetBranchAddress("cl_energy",  ev.cl_energy);
    chain.SetBranchAddress("cl_center",  ev.cl_center);
    chain.SetBranchAddress("cl_time",    ev.cl_time);
    chain.SetBranchAddress("cl_dt_rf",   ev.cl_dt_rf);
    chain.SetBranchAddress("rf_n_a",     &ev.rf_n_a);

    // Per-module histograms.
    const int   n_mods = hycal.module_count();
    const float lo     = -prad2::RF_PERIOD_NS / 2.f;
    const float hi     =  prad2::RF_PERIOD_NS / 2.f;
    const int   nbins  = 100;
    std::vector<TH1F *> hist(n_mods, nullptr);

    // One aggregate "all clusters" histogram and one "high-E single cluster"
    // for the QC summary.
    TH1F h_all  ("h_dt_all",   "All clusters;#Deltat (ns);entries", nbins, lo, hi);
    TH1F h_hiE  ("h_dt_hiE",   "High-E single-cluster;#Deltat (ns);entries",
                 nbins, lo, hi);
    TH1F h_unf  ("h_unfolded", "(cl_time-rf_ns_a[0]) high-E;ns;entries",
                 200, -150, 150);

    // Counters.
    long long n_seen = 0, n_kept = 0;

    const Long64_t N = chain.GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        chain.GetEntry(i);
        ++n_seen;
        if (ev.rf_n_a == 0) continue;
        if (cfg.require_single_cluster && ev.n_clusters != 1) continue;
        for (int k = 0; k < ev.n_clusters && k < prad2::kMaxClusters; ++k) {
            const float dt = ev.cl_dt_rf[k];
            if (!std::isfinite(dt)) continue;
            if (ev.cl_energy[k] < cfg.min_cluster_energy) continue;
            const int center_id = ev.cl_center[k];
            const auto *mod = hycal.module_by_id(center_id);
            if (!mod) continue;
            const int idx = mod->index;
            if (idx < 0 || idx >= n_mods) continue;

            if (!hist[idx]) {
                std::string name = "h_dt_" + mod->name;
                std::string ttl  = mod->name + ";#Deltat (ns);entries";
                hist[idx] = new TH1F(name.c_str(), ttl.c_str(), nbins, lo, hi);
                hist[idx]->SetDirectory(nullptr);  // not associated with any file
            }
            hist[idx]->Fill(dt);
            h_all.Fill(dt);
            if (ev.n_clusters == 1 && ev.cl_energy[k] > 500.f) {
                h_hiE.Fill(dt);
            }
            ++n_kept;
        }
        if ((i & 0xFFFF) == 0 && i > 0) {
            std::cerr << "  scanned " << i << " / " << N << "\n";
        }
    }
    std::cerr << "[fill] seen=" << n_seen << " kept=" << n_kept << "\n";
    std::cerr << "[fill] aggregate folded:  mean=" << h_all.GetMean()
              << "  RMS=" << h_all.GetRMS()
              << "  N=" << (Long64_t)h_all.GetEntries() << "\n";
    std::cerr << "[fill] high-E folded:     mean=" << h_hiE.GetMean()
              << "  RMS=" << h_hiE.GetRMS()
              << "  N=" << (Long64_t)h_hiE.GetEntries() << "\n";

    // Per-module fit.
    struct FitRow { std::string name; int id; long long entries;
                    float mean, sigma, chi2_ndf; std::string status; };
    std::vector<FitRow> rows;
    rows.reserve(n_mods);
    int n_fit_ok = 0, n_fit_skip = 0;
    for (int idx = 0; idx < n_mods; ++idx) {
        const auto &mod = hycal.module(idx);
        if (!mod.is_hycal()) continue;                  // skip Veto/LMS
        FitRow row{mod.name, mod.id, 0, 0.f, 0.f, -1.f, "no-hist"};
        if (hist[idx]) {
            row.entries = static_cast<long long>(hist[idx]->GetEntries());
            if (row.entries >= cfg.min_entries) {
                if (fit_peak(hist[idx], cfg.fit_window_ns,
                             row.mean, row.sigma, row.chi2_ndf, row.status)) {
                    ++n_fit_ok;
                } else {
                    ++n_fit_skip;
                }
            } else {
                row.status = "low-stats";
                ++n_fit_skip;
            }
        } else {
            ++n_fit_skip;
        }
        rows.push_back(std::move(row));
    }
    std::cerr << "[fit ] ok=" << n_fit_ok << "  skipped=" << n_fit_skip << "\n";

    // Emit JSON.
    json j_out;
    j_out["default"] = 0.0;
    json modules = json::array();
    for (const auto &r : rows) {
        if (r.status != "ok") continue;
        json m;
        m["name"]      = r.name;
        m["offset_ns"] = r.mean;
        modules.push_back(m);
    }
    j_out["modules"] = modules;

    // Pretty-print with stable ordering.
    std::ofstream of(cfg.output_json);
    if (!of) {
        std::cerr << "Cannot write " << cfg.output_json << "\n";
        return 1;
    }
    of << j_out.dump(2) << "\n";
    of.close();
    std::cerr << "[out ] wrote " << modules.size() << " offsets to "
              << cfg.output_json << "\n";

    // Optional CSV — full per-module status, including unfitted modules.
    if (!cfg.csv_out.empty()) {
        std::ofstream csv(cfg.csv_out);
        csv << "name,id,entries,mean_ns,sigma_ns,chi2_ndf,status\n";
        csv << std::fixed << std::setprecision(4);
        for (const auto &r : rows) {
            csv << r.name << ',' << r.id << ',' << r.entries << ','
                << r.mean << ',' << r.sigma << ',' << r.chi2_ndf << ','
                << r.status << '\n';
        }
        std::cerr << "[out ] wrote " << cfg.csv_out << "\n";
    }

    // Optional QC PDF: aggregate plot, then one page per fitted module.
    if (!cfg.canvas_pdf.empty()) {
        gStyle->SetOptStat(2200);
        gStyle->SetOptFit(111);
        TCanvas c("c_rf", "RF Δt QC", 1000, 700);
        c.Print((cfg.canvas_pdf + "[").c_str());

        c.cd();
        h_all.SetLineColor(kBlack);
        h_all.Draw();
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.035);
        lat.DrawLatex(0.55, 0.85, Form("N = %lld", (Long64_t)h_all.GetEntries()));
        lat.DrawLatex(0.55, 0.80, Form("mean = %.3f ns", h_all.GetMean()));
        lat.DrawLatex(0.55, 0.75, Form("RMS  = %.3f ns", h_all.GetRMS()));
        c.Print(cfg.canvas_pdf.c_str());

        c.cd();
        h_hiE.SetLineColor(kBlue);
        h_hiE.Draw();
        c.Print(cfg.canvas_pdf.c_str());

        // Per-module pages (up to 100 modules — keep PDF small).
        int n_drawn = 0;
        for (int idx = 0; idx < n_mods && n_drawn < 100; ++idx) {
            if (!hist[idx]) continue;
            if (hist[idx]->GetEntries() < cfg.min_entries) continue;
            c.cd();
            hist[idx]->Draw();
            c.Print(cfg.canvas_pdf.c_str());
            ++n_drawn;
        }
        c.Print((cfg.canvas_pdf + "]").c_str());
        std::cerr << "[out ] wrote " << n_drawn << " module pages to "
                  << cfg.canvas_pdf << "\n";
    }

    // Free histograms (we don't own them through a TFile).
    for (auto *h : hist) delete h;
    return 0;
}
