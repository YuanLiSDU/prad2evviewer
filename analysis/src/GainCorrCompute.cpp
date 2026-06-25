//=============================================================================
// GainCorrCompute.cpp — shared LMS/alpha gain-correction batch calculation
//=============================================================================

#include "GainCorrCompute.h"
#include "EventData.h"
#include "EventData_io.h"

#include <TChain.h>
#include <TFile.h>

#include <algorithm>
#include <filesystem>
#include <iostream>

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
            if (lms_count == 0) ev_start = ev.event_num;
            ++lms_count;
        }
        if (is_alpha) ++alpha_count;

        if (lms_count >= batch_size) {
            batch.batch_id        = batch_id++;
            batch.event_num_start = ev_start;
            batch.event_num_end   = ev.event_num;
            batch.n_lms_events    = lms_count;
            batch.n_alpha_events  = alpha_count;

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
