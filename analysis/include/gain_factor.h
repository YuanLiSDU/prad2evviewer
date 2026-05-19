#pragma once
//=============================================================================
// gain_factor.h — per-module LMS gain factors and time-dependent corrections
//
// ── Reference gain (.dat files, produced by refGain_produce) ─────────────────
//   Directory: <db>/gain_factor/
//   File:      prad_XXXXXX_LMS.dat  (7 columns: Name lms_peak lms_sigma
//              lms_chi2/ndf g1 g2 g3)
//   W/G module lines only; LMS header lines ignored.
//
//   auto tbl = prad2::LoadRefGain(db + "/gain_factor", run_num);
//   float g1_W1 = tbl.w[1].g[0];
//
// ── Time-dependent correction (.root files, produced by replay_gainCorr) ──────
//   Directory: <db>/gain_factor/gain_correction/
//   File:      prad_XXXXXX_gain_corr.root  (TTree "gain_corr", one entry/batch)
//   Only W (PbWO4) modules; G (PbGlass) corrections not stored in these files.
//
//   // One-time setup (single-threaded):
//   auto ts = prad2::LoadGainCorrTimeSeries(
//                 db + "/gain_factor/gain_correction", run_num);
//   // Per-event lookup (read-only → safe from multiple threads after init):
//   const auto& corr = ts.GetCorr(event_num);
//   new_adc2mev = old_adc2mev * corr.w[module_id].avg;
//
// ── Selection rule (both variants) ───────────────────────────────────────────
//   run_num >= 0 -> largest file run number that is <= run_num
//   run_num <  0 -> file with the largest run number ("latest")
//=============================================================================

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

namespace prad2 {

// Single module's three gain factors (g1, g2, g3 from the LMS/alpha fit).
// Zero-initialised by default so missing entries are safe to use.
struct RefGainFactor {
    float g[3] = {0.f, 0.f, 0.f};  // g[0]=g1, g[1]=g2, g[2]=g3
};

// Full table for one run.  Arrays indexed by module numeric ID; index 0 is
// unused (modules are 1-based in the dat file).
struct RefGainTable {
    static constexpr int MAX_W = 1157;  // W1 .. W1156
    static constexpr int MAX_G = 901;   // G1 .. G900

    RefGainFactor w[MAX_W];  // PWO crystal modules
    RefGainFactor g[MAX_G];  // PbGlass modules
    int  run_number = -1;
    bool loaded     = false;
};

// Returns the path to the gain factor file whose run number exactly matches run_num.
// Returns an empty string if no matching file is found or the directory is inaccessible.
inline std::string FindRefGainFile(const std::string &dir, int run_num)
{
    const std::regex pat(R"(prad_(\d{6})_LMS\.dat)");

    std::error_code ec;
    for (auto &entry : std::filesystem::directory_iterator(dir, ec)) {
        if (!entry.is_regular_file()) continue;
        std::smatch m;
        std::string fname = entry.path().filename().string();
        if (!std::regex_match(fname, m, pat)) continue;

        if (std::stoi(m[1].str()) == run_num)
            return entry.path().string();
    }
    if (ec)
        std::cerr << "Warning: cannot iterate gain_factor dir " << dir
                  << ": " << ec.message() << "\n";

    return {};
}

// Load reference gain factors for the given run number.
// On any failure returns a default-constructed (unloaded) table.
inline RefGainTable LoadRefGain(const std::string &dir, int run_num)
{
    RefGainTable tbl;

    std::string path = FindRefGainFile(dir, run_num);
    if (path.empty()) {
        std::cerr << "Warning: no gain factor file found in " << dir
                  << " for run " << run_num << "\n";
        return tbl;
    }

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Warning: cannot open gain factor file " << path << "\n";
        return tbl;
    }

    // Parse the embedded run number from the file name.
    {
        std::regex pat(R"(prad_(\d{6})_LMS\.dat)");
        std::smatch m;
        std::string fname = std::filesystem::path(path).filename().string();
        if (std::regex_search(fname, m, pat))
            tbl.run_number = std::stoi(m[1].str());
    }

    // Skip the header line (contains string column names, not numeric data).
    { std::string header; std::getline(f, header); }

    std::string name;
    float col2, col3, col4, g1, g2, g3;
    while (f >> name >> col2 >> col3 >> col4 >> g1 >> g2 >> g3) {
        if (name.empty()) continue;

        if (name[0] == 'W') {
            int id = std::stoi(name.substr(1));
            if (id >= 1 && id < RefGainTable::MAX_W) {
                tbl.w[id].g[0] = g1;
                tbl.w[id].g[1] = g2;
                tbl.w[id].g[2] = g3;
            }
        } else if (name[0] == 'G') {
            int id = std::stoi(name.substr(1));
            if (id >= 1 && id < RefGainTable::MAX_G) {
                tbl.g[id].g[0] = g1;
                tbl.g[id].g[1] = g2;
                tbl.g[id].g[2] = g3;
            }
        }
        // LMS lines and anything else are silently skipped.
    }

    tbl.loaded = true;
    std::cerr << "RefGain: loaded run " << tbl.run_number
              << " from " << path << "\n";
    return tbl;
}

// Per-module gain correction factor: correction[id] = g_ref / g_current.
// Applying this to the current ADC->MeV scale compensates for gain drift.
// A value of 1.0 means no correction needed; > 1.0 means gain dropped.
struct GainCorrTable {
    static constexpr int MAX_W = RefGainTable::MAX_W;
    static constexpr int MAX_G = RefGainTable::MAX_G;

    // correction[id][j] = ref.g[j] / cur.g[j]  (j = 0,1,2 for g1,g2,g3)
    // avg[id]           = mean of the valid (non-zero) per-LMS corrections
    struct Entry {
        float corr[3] = {1.f, 1.f, 1.f};  // per-LMS correction
        float avg      = 1.f;              // average over valid LMS channels
    };

    Entry w[MAX_W];
    Entry g[MAX_G];

    int ref_run = -1;
    int cur_run = -1;
};

// ── Time-dependent gain correction from replay_gainCorr ROOT output ───────────

// Locate prad_XXXXXX_gain_corr.root for the given run.
// Same selection rule as FindRefGainFile.
inline std::string FindGainCorrRootFile(const std::string &dir, int run_num)
{
    const std::regex pat(R"(prad_(\d{6})_gain_corr\.root)");

    std::error_code ec;
    for (auto &entry : std::filesystem::directory_iterator(dir, ec)) {
        if (!entry.is_regular_file()) continue;
        std::smatch m;
        std::string fname = entry.path().filename().string();
        if (!std::regex_match(fname, m, pat)) continue;

        if (std::stoi(m[1].str()) == run_num)
            return entry.path().string();
    }
    if (ec)
        std::cerr << "Warning: cannot iterate gain_correction dir " << dir
                  << ": " << ec.message() << "\n";

    return {};
}

// A time series of gain correction tables loaded from a replay_gainCorr ROOT
// file (one entry per LMS batch).  After construction the object is read-only
// and therefore safe to access concurrently from multiple threads.
struct GainCorrTimeSeries {
    struct Batch {
        int           event_num_start = 0;
        int           event_num_end   = 0;
        GainCorrTable corr;
    };

    std::vector<Batch> batches;   // sorted ascending by event_num_start
    int  run_num = -1;
    bool loaded  = false;

    // Return the correction table for the given event number.
    // Finds the last batch whose event_num_start <= event_num.
    // Falls back to the first batch if event_num precedes all batches.
    // Thread-safe after LoadGainCorrTimeSeries() returns: purely read-only.
    const GainCorrTable &GetCorr(int event_num) const noexcept
    {
        // C++11 guarantees thread-safe initialisation of function-scope statics.
        static const GainCorrTable kIdentity{};
        if (batches.empty()) return kIdentity;
        for (auto it = batches.rbegin(); it != batches.rend(); ++it)
            if (event_num >= it->event_num_start) return it->corr;
        return batches.front().corr;
    }
};

// Load the gain correction time series from the replay_gainCorr ROOT output.
//
// corr_dir : directory containing prad_XXXXXX_gain_corr.root files
//            (typically <db>/gain_factor/gain_correction)
// run_num  : select file by run number (see selection rule in header)
//
// Thread safety
//   Call once from a single thread during setup.
//   The returned GainCorrTimeSeries may then be shared across threads for
//   read-only access via GetCorr().  Requires ROOT::EnableThreadSafety() to
//   have been called before spawning worker threads.
inline GainCorrTimeSeries LoadGainCorrTimeSeries(const std::string &corr_dir,
                                                  int               run_num)
{
    GainCorrTimeSeries ts;
    ts.run_num = run_num;

    std::string path = FindGainCorrRootFile(corr_dir, run_num);
    if (path.empty()) {
        std::cerr << "Warning: no gain_corr root file in " << corr_dir
                  << " for run " << run_num << "\n";
        return ts;
    }

    // Open privately — no gDirectory/gFile side-effects on the caller's context.
    TFile *f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Warning: cannot open " << path << "\n";
        delete f;
        return ts;
    }

    TTree *tree = nullptr;
    f->GetObject("gain_corr", tree);
    if (!tree) {
        std::cerr << "Warning: no 'gain_corr' tree in " << path << "\n";
        delete f;
        return ts;
    }

    // Written with gain_corr_W[N_W][N_LMS], N_W = 1156, N_LMS = 3.
    static constexpr int kNW   = GainCorrTable::MAX_W - 1;  // 1156
    static constexpr int kNLMS = 3;

    int   ev_start = 0, ev_end = 0;
    float corr_W[kNW][kNLMS];

    tree->SetBranchAddress("event_num_start", &ev_start);
    tree->SetBranchAddress("event_num_end",   &ev_end);
    tree->SetBranchAddress("gain_corr_W",      corr_W);

    const Long64_t nentries = tree->GetEntries();
    ts.batches.reserve(static_cast<size_t>(nentries));

    for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);

        GainCorrTimeSeries::Batch b;
        b.event_num_start = ev_start;
        b.event_num_end   = ev_end;
        b.corr.cur_run    = run_num;

        for (int wi = 0; wi < kNW; ++wi) {
            float sum = 0.f;
            for (int j = 0; j < kNLMS; ++j) {
                // 0 in the ROOT file signals a failed fit — treat as identity.
                float v = (corr_W[wi][j] > 0.f) ? corr_W[wi][j] : 1.f;
                b.corr.w[wi + 1].corr[j] = v;
                sum += v;
            }
            b.corr.w[wi + 1].avg = sum / static_cast<float>(kNLMS);
        }
        ts.batches.push_back(std::move(b));
    }

    delete f;   // TFile owns TTree; both are freed here

    // Ensure ascending order (defensive; writer should already sort).
    std::sort(ts.batches.begin(), ts.batches.end(),
              [](const GainCorrTimeSeries::Batch &a,
                 const GainCorrTimeSeries::Batch &b) {
                  return a.event_num_start < b.event_num_start;
              });

    ts.loaded = true;
    std::cerr << "GainCorrTS: run " << run_num << ": "
              << ts.batches.size() << " batches from " << path << "\n";
    return ts;
}

//reuse the code in gain_fitter.cpp
struct FitResult{
    float mean     = 0.;
    float sigma    = 0.;
    float chi2pndf = 0.;
};

inline FitResult gain_hist_fitter(TH1F* h, const float & fac)
{
    FitResult result;

    if (!h) {
        return result;
    }

    const int nBins = h->GetNbinsX();
    if (nBins <= 0) {
        return result;
    }

    const int maxBin = h->GetMaximumBin();
    const double maxContent = h->GetBinContent(maxBin);

    if (maxContent <= 0.) {
        return result;
    }

    const double threshold = fac * maxContent;

    int leftBin  = maxBin;
    int rightBin = maxBin;

    // Find first bin to the left below 10% of max
    for (int ibin = maxBin; ibin >= 1; --ibin) {
        if (h->GetBinContent(ibin) < threshold) {
            leftBin = ibin;
            break;
        }
        if (ibin == 1) {
            leftBin = 1;
        }
    }

    // Find first bin to the right below 10% of max
    for (int ibin = maxBin; ibin <= nBins; ++ibin) {
        if (h->GetBinContent(ibin) < threshold) {
            rightBin = ibin;
            break;
        }
        if (ibin == nBins) {
            rightBin = nBins;
        }
    }

    // Optional: if you want the fit range to stay inside the above-threshold region,
    // shift inward by one bin when the threshold-crossing bin itself is below threshold.
    if (leftBin < maxBin && h->GetBinContent(leftBin) < threshold) {
        leftBin++;
    }
    if (rightBin > maxBin && h->GetBinContent(rightBin) < threshold) {
        rightBin--;
    }

    // Ensure at least 5 bins around the maximum are included so that a
    // 3-parameter Gaussian fit is always well-constrained, even when the
    // peak is very narrow (only a few filled bins).
    if (leftBin  > maxBin - 2) leftBin  = maxBin - 2;
    if (rightBin < maxBin + 2) rightBin = maxBin + 2;

    // Safety check
    if (leftBin < 1) leftBin = 1;
    if (rightBin > nBins) rightBin = nBins;
    if (leftBin >= rightBin) {
        return result;
    }
    
    const double xLow  = h->GetXaxis()->GetBinLowEdge(leftBin);
    const double xHigh = h->GetXaxis()->GetBinUpEdge(rightBin);

    const double peakX = h->GetXaxis()->GetBinCenter(maxBin);

    // A reasonable initial sigma guess from fit window width
    double sigmaGuess = 0.5 * (xHigh - xLow) / 2.0;
    if (sigmaGuess <= 0.) {
        sigmaGuess = h->GetRMS();
    }
    if (sigmaGuess <= 0.) {
        sigmaGuess = h->GetBinWidth(maxBin);
    }

    std::string fitName = std::string(h->GetName()) + "_gaus_fit";
    TF1 * gausFit = new TF1(fitName.c_str(), "gaus", xLow, xHigh);
    gausFit->SetParameters(maxContent, peakX, sigmaGuess);

    // R = fit in function range, Q = quiet, N = do not store function in histogram
    // "N" is required to avoid double-ownership: new TF1 is registered in gROOT's
    // global list; without "N", Fit() also stores it in the histogram's list,
    // causing a double-free during ROOT cleanup at program exit.
    int fitStatus = h->Fit(gausFit, "RQN");

    if (fitStatus != 0) {
        delete gausFit;
        return result;
    }

    const float mean  = static_cast<float>(gausFit->GetParameter(1));
    const float sigma = static_cast<float>(gausFit->GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.f) {
        delete gausFit;
        return result;
    }
    result.mean  = mean;
    result.sigma = sigma;

    const double ndf = gausFit->GetNDF();
    if (ndf > 0) {
        result.chi2pndf = static_cast<float>(gausFit->GetChisquare() / ndf);
    }

    // Transfer ownership from gROOT's global list to the histogram so the fit
    // curve is saved in the output file and there is only one owner (no double-free).
    gROOT->GetListOfFunctions()->Remove(gausFit);
    h->GetListOfFunctions()->Add(gausFit);

    return result;
}

} // namespace prad2

