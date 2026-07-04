//=============================================================================
// GainCorrCompute.cpp — shared LMS/alpha gain-correction batch calculation
//=============================================================================

#include "GainCorrCompute.h"
#include "EventData.h"
#include "EventData_io.h"

#include <TChain.h>
#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>

namespace analysis {
namespace {

TH1F *makeH(const char *name)
{
    return new TH1F(name, name, kGainHistBins, kGainHistMin, kGainHistMax);
}

void resetH(TH1F *h)
{
    h->Reset("ICESM");
}

struct UnixTimeAnchor {
    int       event_number = -1;
    long long ti_ticks     = 0;
    uint32_t  unix_time    = 0;
};

void loadUnixTimeAnchors(const std::vector<std::string> &files,
                         std::vector<UnixTimeAnchor>     &anchors,
                         bool                            use_epics)
{
    for (const auto &path : files) {
        TFile file(path.c_str(), "READ");
        if (file.IsZombie()) continue;

        if (use_epics) {
            auto *tree = file.Get<TTree>("epics");
            if (!tree) continue;

            prad2::RawEpicsData row;
            prad2::SetEpicsReadBranches(tree, row);
            for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
                row = prad2::RawEpicsData{};
                tree->GetEntry(i);
                if (row.unix_time == 0) continue;
                anchors.push_back({row.event_number_at_arrival,
                                   row.ti_ticks_at_arrival,
                                   row.unix_time});
            }
        } else {
            auto *tree = file.Get<TTree>("scalers");
            if (!tree) continue;

            prad2::RawScalerData row;
            prad2::SetScalerReadBranches(tree, row);
            for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
                row = prad2::RawScalerData{};
                tree->GetEntry(i);
                if (row.unix_time == 0) continue;
                anchors.push_back({row.event_number, row.ti_ticks,
                                   row.unix_time});
            }
        }
    }
}

uint32_t batchUnixTime(const std::vector<UnixTimeAnchor> &anchors,
                       long long batch_ti_ticks,
                       int       batch_event_number)
{
    const UnixTimeAnchor *best = nullptr;

    // TI timestamps give the precise relative position on the run timeline.
    if (batch_ti_ticks > 0) {
        uint64_t best_distance = std::numeric_limits<uint64_t>::max();
        for (const auto &anchor : anchors) {
            if (anchor.ti_ticks <= 0) continue;
            const uint64_t distance = static_cast<uint64_t>(
                std::llabs(batch_ti_ticks - anchor.ti_ticks));
            if (distance < best_distance) {
                best = &anchor;
                best_distance = distance;
            }
        }
    }

    // Legacy files may have unix_time but no ti_ticks_at_arrival branch.
    if (!best && batch_event_number >= 0) {
        long long best_distance = std::numeric_limits<long long>::max();
        for (const auto &anchor : anchors) {
            if (anchor.event_number < 0) continue;
            const long long distance = std::llabs(
                static_cast<long long>(batch_event_number) - anchor.event_number);
            if (distance < best_distance) {
                best = &anchor;
                best_distance = distance;
            }
        }
    }

    if (!best) return 0;

    constexpr double kTiTickSeconds = 4.0e-9;  // 250 MHz TI clock
    double unix_seconds = static_cast<double>(best->unix_time);
    if (batch_ti_ticks > 0 && best->ti_ticks > 0)
        unix_seconds += (batch_ti_ticks - best->ti_ticks) * kTiTickSeconds;

    if (!std::isfinite(unix_seconds) || unix_seconds <= 0.0 ||
        unix_seconds > static_cast<double>(std::numeric_limits<uint32_t>::max()))
        return 0;
    return static_cast<uint32_t>(std::llround(unix_seconds));
}

} // namespace

GainPlotStore::~GainPlotStore()
{
    for (int i = 0; i < kGainNLMS; ++i) {
        for (auto *h : ref_lms[i])   delete h;
        for (auto *h : ref_alpha[i]) delete h;
    }
    for (auto &kv : mod_w)
        for (auto *h : kv.second) delete h;
}

std::string MakeLMSOutputFile(const std::string &evio_path)
{
    std::string out = std::filesystem::path(evio_path).filename().string();
    auto pos = out.find(".evio");
    if (pos != std::string::npos)
        out = out.substr(0, pos) + out.substr(pos + 5);
    out += "_lms.root";
    return out;
}

void SetupGainBranches(TTree *tree, GainBatch &b)
{
    tree->Branch("batch_id",        &b.batch_id,        "batch_id/I");
    tree->Branch("event_num_start", &b.event_num_start, "event_num_start/I");
    tree->Branch("event_num_end",   &b.event_num_end,   "event_num_end/I");
    tree->Branch("n_lms_events",    &b.n_lms_events,    "n_lms_events/I");
    tree->Branch("n_alpha_events",  &b.n_alpha_events,  "n_alpha_events/I");
    tree->Branch("ref_run",         &b.ref_run,         "ref_run/I");
    tree->Branch("unix_time",       &b.unix_time,       "unix_time/i");

    tree->Branch("refPMT_ratio",       b.refPMT_ratio,
                 Form("refPMT_ratio[%d]/F",       kGainNLMS));
    tree->Branch("gain_W",             b.gain_W,
                 Form("gain_W[%d][%d]/F",         kGainNW, kGainNLMS));
    tree->Branch("gain_W_ref",         b.gain_W_ref,
                 Form("gain_W_ref[%d][%d]/F",     kGainNW, kGainNLMS));
    tree->Branch("gain_corr_W",        b.gain_corr_W,
                 Form("gain_corr_W[%d][%d]/F",    kGainNW, kGainNLMS));
    tree->Branch("fit_mean_ref_lms",   b.fit_mean_ref_lms,
                 Form("fit_mean_ref_lms[%d]/F",   kGainNLMS));
    tree->Branch("fit_mean_ref_alpha", b.fit_mean_ref_alpha,
                 Form("fit_mean_ref_alpha[%d]/F", kGainNLMS));
    tree->Branch("fit_mean_W_lms",     b.fit_mean_W_lms,
                 Form("fit_mean_W_lms[%d]/F",     kGainNW));
}

void FlushGainBatch(GainBatch &b, TTree *tree,
                    TH1F *mod_lms[kGainNW],
                    TH1F *ref_lms[kGainNLMS],
                    TH1F *ref_alpha[kGainNLMS],
                    const prad2::RefGainTable &ref_tbl)
{
    prad2::FitResult fit_ref_lms[kGainNLMS], fit_ref_alpha[kGainNLMS];
    for (int i = 0; i < kGainNLMS; ++i) {
        fit_ref_lms[i]   = prad2::gain_hist_fitter(ref_lms[i],   0.1f);
        fit_ref_alpha[i] = prad2::gain_hist_fitter(ref_alpha[i], 0.1f);
        b.fit_mean_ref_lms[i]   = fit_ref_lms[i].mean;
        b.fit_mean_ref_alpha[i] = fit_ref_alpha[i].mean;
        b.refPMT_ratio[i] = (fit_ref_alpha[i].mean > 0.f)
            ? fit_ref_lms[i].mean / fit_ref_alpha[i].mean : 0.f;
    }

    for (int i = 0; i < kGainNW; ++i) {
        prad2::FitResult fr = prad2::gain_hist_fitter(mod_lms[i], 0.1f);
        b.fit_mean_W_lms[i] = fr.mean;
        for (int j = 0; j < kGainNLMS; ++j) {
            b.gain_W[i][j] = (b.refPMT_ratio[j] > 0.f)
                ? fr.mean / b.refPMT_ratio[j] : 0.f;
            b.gain_W_ref[i][j]  = ref_tbl.w[i + 1].g[j];
            b.gain_corr_W[i][j] = (b.gain_W[i][j] > 0.f)
                ? b.gain_W_ref[i][j] / b.gain_W[i][j] : 1.f;
        }
    }

    tree->Fill();
}

bool ComputeGainCorrections(const std::vector<std::string> &lms_files,
                            const std::string             &gain_out,
                            int                            batch_size,
                            int                            ref_run_num,
                            const prad2::RefGainTable     &ref_tbl,
                            const GainPlotConfig          *plot_cfg,
                            GainPlotStore                 *plot_store)
{
    if (lms_files.empty()) {
        std::cerr << "No LMS files provided for gain calculation.\n";
        return false;
    }
    if (batch_size <= 0) {
        std::cerr << "Invalid batch size: " << batch_size << "\n";
        return false;
    }

    TChain chain("lms_gain");
    for (auto &f : lms_files) {
        chain.Add(f.c_str());
        std::cout << "  [gain] + " << f << "\n";
    }
    std::cout << "  [gain] total entries: " << chain.GetEntries() << "\n";

    prad2::LMSEventData ev;
    prad2::SetLMSReadBranches(&chain, ev);

    // Process_LMSgainFactor writes EPICS and scaler side trees alongside
    // lms_gain.  EPICS provides the preferred absolute-time pins; scaler
    // pins are a fallback for runs/files without usable EPICS records.
    std::vector<UnixTimeAnchor> time_anchors;
    loadUnixTimeAnchors(lms_files, time_anchors, true);
    const char *time_anchor_source = "EPICS";
    if (time_anchors.empty()) {
        loadUnixTimeAnchors(lms_files, time_anchors, false);
        time_anchor_source = "scaler";
    }
    std::cout << "  [gain] time anchors: " << time_anchors.size();
    if (!time_anchors.empty()) std::cout << " (" << time_anchor_source << ")";
    std::cout << "\n";

    TH1F *mod_lms  [kGainNW];
    TH1F *ref_lms  [kGainNLMS];
    TH1F *ref_alpha[kGainNLMS];
    for (int i = 0; i < kGainNW;   ++i) mod_lms[i]   = makeH(Form("mod_lms_%d",   i + 1));
    for (int i = 0; i < kGainNLMS; ++i) ref_lms[i]   = makeH(Form("ref_lms_%d",   i + 1));
    for (int i = 0; i < kGainNLMS; ++i) ref_alpha[i] = makeH(Form("ref_alpha_%d", i + 1));

    TFile *outfile = TFile::Open(gain_out.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Cannot create gain output file: " << gain_out << "\n";
        delete outfile;
        for (int i = 0; i < kGainNW;   ++i) delete mod_lms[i];
        for (int i = 0; i < kGainNLMS; ++i) { delete ref_lms[i]; delete ref_alpha[i]; }
        return false;
    }
    TTree *out_tree = new TTree("gain_corr", "LMS gain correction time series");
    GainBatch batch;
    SetupGainBranches(out_tree, batch);
    batch.ref_run = ref_run_num;

    Long64_t nentries    = chain.GetEntries();
    int      lms_count   = 0;
    int      alpha_count = 0;
    int      batch_id    = 0;
    int      ev_start    = 0;
    long long ti_start   = 0;

    auto captureForPlot = [&]() {
        if (!plot_cfg || !plot_cfg->enabled || !plot_store) return;
        for (int i = 0; i < kGainNLMS; ++i) {
            auto *hl = (TH1F*)ref_lms[i]->Clone(
                Form("ref_lms%d_b%d", i + 1, batch.batch_id));
            hl->SetDirectory(nullptr);
            plot_store->ref_lms[i].push_back(hl);
            auto *ha = (TH1F*)ref_alpha[i]->Clone(
                Form("ref_alpha%d_b%d", i + 1, batch.batch_id));
            ha->SetDirectory(nullptr);
            plot_store->ref_alpha[i].push_back(ha);
        }
        for (int wid : plot_cfg->w_ids) {
            if (wid >= 1 && wid <= kGainNW) {
                auto *hw = (TH1F*)mod_lms[wid - 1]->Clone(
                    Form("mod_lms%d_b%d", wid, batch.batch_id));
                hw->SetDirectory(nullptr);
                plot_store->mod_w[wid].push_back(hw);
            }
        }
    };

    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        chain.GetEntry(ientry);

        const bool is_lms   = (ev.event_type == 0);
        const bool is_alpha = (ev.event_type == 1);
        if (!is_lms && !is_alpha) continue;

        for (int ich = 0; ich < ev.nch; ++ich) {
            if (ev.npeaks[ich] != 1) continue;
            float h = ev.peak_integral[ich][0];
            if (h <= 0.f) continue;

            const int mid  = ev.module_id[ich];
            const int mtyp = ev.module_type[ich];

            if (mtyp == prad2::MOD_PbWO4) {
                int wid = mid - kGainWIDBase;
                if (wid < 1 || wid > kGainNW) continue;
                if (is_lms) mod_lms[wid - 1]->Fill(h);
            } else if (mtyp == prad2::MOD_LMS) {
                int lid = mid - kGainLMSIDBase;
                if (lid < 0 || lid >= kGainNLMS) continue;
                if (is_lms)   ref_lms[lid]->Fill(h);
                if (is_alpha) ref_alpha[lid]->Fill(h);
            }
        }

        if (is_lms) {
            if (lms_count == 0) {
                ev_start = ev.event_num;
                ti_start = ev.timestamp;
            }
            ++lms_count;
        }
        if (is_alpha) ++alpha_count;

        if (lms_count >= batch_size) {
            batch.batch_id        = batch_id++;
            batch.event_num_start = ev_start;
            batch.event_num_end   = ev.event_num;
            batch.n_lms_events    = lms_count;
            batch.n_alpha_events  = alpha_count;
            const long long ti_mid = ti_start + (ev.timestamp - ti_start) / 2;
            const int event_mid = ev_start + (ev.event_num - ev_start) / 2;
            batch.unix_time = batchUnixTime(time_anchors, ti_mid, event_mid);

            FlushGainBatch(batch, out_tree, mod_lms, ref_lms, ref_alpha, ref_tbl);
            captureForPlot();

            for (int i = 0; i < kGainNW;   ++i) resetH(mod_lms[i]);
            for (int i = 0; i < kGainNLMS; ++i) { resetH(ref_lms[i]); resetH(ref_alpha[i]); }
            lms_count   = 0;
            alpha_count = 0;

            std::cerr << "\r  [gain] batch " << batch_id << " written" << std::endl;
        }
    }

    std::cerr << "\n";
    outfile->cd();
    out_tree->Write();
    delete outfile;
    std::cout << "  [gain] " << batch_id << " batches written to " << gain_out << "\n";

    for (int i = 0; i < kGainNW;   ++i) delete mod_lms[i];
    for (int i = 0; i < kGainNLMS; ++i) { delete ref_lms[i]; delete ref_alpha[i]; }
    return true;
}

} // namespace analysis
