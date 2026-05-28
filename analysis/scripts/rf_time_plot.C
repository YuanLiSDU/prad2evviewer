//============================================================================
// rf_time_plot.C — folded vs. unfolded HyCal→RF time on a recon file.
//
// Reads `cl_dt_rf` (folded, with per-module offsets applied) and
// re-derives the unfolded version from `cl_time` + `rf_ns_a[]` so you can
// see the bunch-resolved peak collapse onto the 4.008 ns lattice.
//
// Usage:
//   cd <build_dir>
//   root -l ../analysis/scripts/rootlogon.C
//   .x ../analysis/scripts/rf_time_plot.C+( \
//       "/path/to/prad_NNNNNN.NNNNN_recon.root", \
//       "rf_plot.pdf",  // output PDF (optional)
//       200.f,           // min cluster energy MeV
//       -1)              // module id to single out, or -1 for all
//============================================================================

#include "EventData.h"
#include "EventData_io.h"
#include "RfTime.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <string>

void rf_time_plot(const char *infile,
                  const char *outfile = "rf_plot.pdf",
                  float min_cluster_energy = 200.f,
                  int   single_center_id = -1)
{
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[rf_time_plot] cannot open " << infile << "\n";
        return;
    }
    TTree *t = dynamic_cast<TTree *>(fin->Get("recon"));
    if (!t) {
        std::cerr << "[rf_time_plot] no 'recon' tree in " << infile << "\n";
        return;
    }
    if (!t->GetBranch("cl_dt_rf")) {
        std::cerr << "[rf_time_plot] no cl_dt_rf branch — re-replay with "
                     "the RF reconstruction changes (see "
                     "docs/analysis_notes/rf_time_reconstruction_plan.md).\n";
        return;
    }

    // Bind branches.  We only need a handful — turn the rest off for speed.
    prad2::ReconEventData ev;
    t->SetBranchStatus("*", 0);
    for (auto *b : {"n_clusters", "cl_energy", "cl_center", "cl_time",
                    "cl_dt_rf", "rf_n_a", "rf_ns_a"}) {
        t->SetBranchStatus(b, 1);
    }
    t->SetBranchAddress("n_clusters", &ev.n_clusters);
    t->SetBranchAddress("cl_energy",  ev.cl_energy);
    t->SetBranchAddress("cl_center",  ev.cl_center);
    t->SetBranchAddress("cl_time",    ev.cl_time);
    t->SetBranchAddress("cl_dt_rf",   ev.cl_dt_rf);
    t->SetBranchAddress("rf_n_a",     &ev.rf_n_a);
    t->SetBranchAddress("rf_ns_a",    ev.rf_ns_a);

    // Histograms.
    const float T_RF  = prad2::RF_PERIOD_NS;
    TH1F h_fold("h_fold",
                "Folded #Deltat (cl_time #minus nearest RF, with per-module offsets);"
                "#Deltat (ns);clusters",
                100, -T_RF / 2.f, T_RF / 2.f);
    TH1F h_unfold("h_unfold",
                  "Unfolded #Deltat (cl_time #minus first RF tick);"
                  "#Deltat (ns);clusters",
                  300, -prad2::RF_DIV_NS, prad2::RF_DIV_NS);
    TH1F h_lattice("h_lattice",
                   "Unfolded with 4.008 ns lattice highlighted;#Deltat (ns);clusters",
                   400, -2.f * prad2::RF_DIV_NS, 2.f * prad2::RF_DIV_NS);

    const Long64_t N = t->GetEntries();
    long long n_kept = 0;
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        if (ev.rf_n_a == 0) continue;
        if (ev.n_clusters != 1) continue;
        for (int k = 0; k < ev.n_clusters; ++k) {
            if (ev.cl_energy[k] < min_cluster_energy) continue;
            if (single_center_id >= 0 && ev.cl_center[k] != single_center_id)
                continue;
            const float dt_folded = ev.cl_dt_rf[k];
            if (!std::isfinite(dt_folded)) continue;
            h_fold.Fill(dt_folded);
            // Unfolded: cl_time - earliest RF tick.  Range will be wide
            // (~ -T_div, +T_div) since cluster time can sit on any of the
            // 32 hidden bunches.
            const float dt_unf = ev.cl_time[k] - ev.rf_ns_a[0];
            h_unfold.Fill(dt_unf);
            h_lattice.Fill(dt_unf);
            ++n_kept;
        }
    }

    std::cout << "[rf_time_plot] scanned=" << N << " kept=" << n_kept << "\n";
    std::cout << "[rf_time_plot] folded:    mean=" << h_fold.GetMean()
              << " RMS=" << h_fold.GetRMS()
              << " N=" << (Long64_t)h_fold.GetEntries() << "\n";
    std::cout << "[rf_time_plot] unfolded:  mean=" << h_unfold.GetMean()
              << " RMS=" << h_unfold.GetRMS()
              << " N=" << (Long64_t)h_unfold.GetEntries() << "\n";

    // Draw.
    gStyle->SetOptStat(2200);
    TCanvas c("c_rf", "RF folding", 1200, 800);
    c.Divide(2, 2);

    c.cd(1);
    h_unfold.SetLineColor(kBlue + 1);
    h_unfold.Draw();
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.15, 0.84, "Unfolded (cl_time #minus rf_ns_a[0])");
    lat.DrawLatex(0.15, 0.78, Form("RMS = %.2f ns", h_unfold.GetRMS()));

    c.cd(2);
    h_fold.SetLineColor(kRed + 1);
    h_fold.Draw();
    lat.DrawLatex(0.15, 0.84, "Folded mod 4.008 ns (with offsets)");
    lat.DrawLatex(0.15, 0.78, Form("RMS = %.2f ns", h_fold.GetRMS()));
    if (single_center_id >= 0)
        lat.DrawLatex(0.15, 0.72, Form("center id = %d", single_center_id));

    c.cd(3);
    h_lattice.SetLineColor(kBlue + 1);
    h_lattice.Draw();
    lat.DrawLatex(0.15, 0.84, "Unfolded wide range — 4.008 ns lattice visible");
    // Draw vertical lines on the 4.008 ns lattice (every 8th line for
    // readability).
    {
        for (int n = -15; n <= 15; ++n) {
            const float x = n * 8.f * T_RF;
            TLine *ln = new TLine(x, 0, x, h_lattice.GetMaximum() * 0.6);
            ln->SetLineColor(kRed + 1);
            ln->SetLineStyle(2);
            ln->Draw();
        }
    }

    c.cd(4);
    // Log-y of the folded distribution to expose the background floor.
    gPad->SetLogy();
    h_fold.Draw();
    lat.DrawLatex(0.15, 0.84, "Folded, log-y");

    if (outfile && outfile[0]) {
        c.SaveAs(outfile);
        std::cout << "[rf_time_plot] wrote " << outfile << "\n";
    }

    fin->Close();
}
