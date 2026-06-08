
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMarker.h>
#include <TLegend.h>

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

using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);
void process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal, 
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events = -1);

float resolution = 0.035; // pre-defined energy resolution

bool Vetoed(float cl_time, float sci_time, float sci_int){
    // Simple veto logic: if the cluster time is within a certain window of the scintillator time, and the scintillator signal is above a threshold, we consider it a vetoed event.
    const float time_shift = 35.f; // ns
    const float time_window = 7.f; // ns
    const float int_threshold = 2000.f; // arbitrary units
    return (fabs(cl_time - sci_time - time_shift) < time_window) && (sci_int > int_threshold);
}

int main(int argc, char *argv[]){

    std::string output, input_3p5, input_0p7;
    std::string pngDir = "module_hists";
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:n:p:")) != -1) {
        switch (opt) {
            case 'a': input_3p5 = optarg; break;
            case 'b': input_0p7 = optarg; break;
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'p': pngDir = optarg; break;
        }
    }
    

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // Energy histogram for each crystal
    std::map<int, TH1F*> energy_hists_3p5;
    std::map<int, TH1F*> energy_hists_0p7;
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        if (!m.is_pwo4()) continue;
        std::string hname = "h_energy_" + m.name;
        std::string htitle = "Energy " + m.name + ";E (MeV);Counts";
        energy_hists_3p5[m.id] = new TH1F((hname + "_3p5").c_str(), (htitle + " (3.5)").c_str(), 400, 0, 4000);
        energy_hists_0p7[m.id] = new TH1F((hname + "_0p7").c_str(), (htitle + " (0.7)").c_str(), 400, 0, 4000);
    }

    // --- setup TChain and branches ---
    TChain *chain = new TChain("recon");
        chain->Add(input_3p5.c_str());
        std::cerr << "Added file: " << input_3p5 << "\n";
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);

    process_event(tree, ev, hycal, energy_hists_3p5, physics, 3485.f, max_events);

    // --- repeat for 0.7 GeV data ---
    TChain *chain2 = new TChain("recon");
        chain2->Add(input_0p7.c_str());
        std::cerr << "Added file: " << input_0p7 << "\n";
    TTree *tree2 = chain2;
    if (!tree2) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev2;
    prad2::SetReconReadBranches(tree2, ev2);
    process_event(tree2, ev2, hycal, energy_hists_0p7, physics, 729.f, max_events);

    std::string calib_file = dbDir + "/" + "calibration/calibration_factor_3p5_June1.json";
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    TH1F *h_energy_peak_3p5 = new TH1F("h_energy_peak_3p5", "Energy Peak Distribution;Energy (MeV);Counts", 4000, 0, 4000);

    // calculate non-linearity module by module and save to output file
    gSystem->mkdir(pngDir.c_str(), true);
    TFile outFile(output.empty() ? "nonlinearity_results.root" : output.c_str(), "RECREATE");
    for (int i = 0; i < hycal.module_count(); i++) {
        const auto &mod = hycal.module(i);
        int mod_id = mod.id;
        if (!mod.is_pwo4()) continue; // only look at PbWO4 crystals

        auto it_3p5 = energy_hists_3p5.find(mod_id);
        if (it_3p5 == energy_hists_3p5.end() || !it_3p5->second) continue;
        auto hist_3p5 = it_3p5->second;
        auto it_0p7 = energy_hists_0p7.find(mod_id);
        if (it_0p7 == energy_hists_0p7.end() || !it_0p7->second) continue;
        auto hist_0p7 = it_0p7->second;

        float x = mod.x, y = mod.y, z = 6270.f;
        float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
        float e_p_exp_3p5 = physics.ExpectedEnergy(theta, 3485.f, "ep");
        float e_e_exp_3p5 = physics.ExpectedEnergy(theta, 3485.f, "ee");
        float e_p_exp_0p7 = physics.ExpectedEnergy(theta, 729.f, "ep");
        float e_e_exp_0p7 = physics.ExpectedEnergy(theta, 729.f, "ee");

        float sigma_ep_3p5 = resolution * e_p_exp_3p5 / sqrt(e_p_exp_3p5/1000.f);
        float sigma_ee_3p5 = resolution * e_e_exp_3p5 / sqrt(e_e_exp_3p5/1000.f);
        float sigma_ep_0p7 = resolution * e_p_exp_0p7 / sqrt(e_p_exp_0p7/1000.f);
        float sigma_ee_0p7 = resolution * e_e_exp_0p7 / sqrt(e_e_exp_0p7/1000.f);

        int _fit_uid = mod_id * 10;

        // Centroid-based Gaussian fit in [Eexp±3σ]:
        // weighted-average centroid → initial mean; walk from peak bin to 50% height for fit range.
        auto fitPeakAndDraw = [&](TH1F *h, double Eexp, double sigma, int color) -> float {
            ++_fit_uid;
            int b0 = std::max(1, h->FindBin(Eexp - 3.*sigma));
            int b1 = std::min(h->GetNbinsX(), h->FindBin(Eexp + 3.*sigma));
            double ypad_min = 0;
            double ypad = h->GetMaximum() * 1.2;
            // search-window: dashed, full y-axis height
            TLine *wl = new TLine(Eexp - 3.*sigma, ypad_min, Eexp - 3.*sigma, ypad);
            wl->SetLineColor(6); wl->SetLineStyle(7); wl->SetLineWidth(2); wl->Draw();
            TLine *wr = new TLine(Eexp + 3.*sigma, ypad_min, Eexp + 3.*sigma, ypad);
            wr->SetLineColor(6); wr->SetLineStyle(7); wr->SetLineWidth(2); wr->Draw();
            // weighted centroid as initial mean
            double wsum = 0., wpos = 0.;
            for (int ib = b0; ib <= b1; ++ib) {
                double c = h->GetBinContent(ib);
                if (c > 0.) { wsum += c; wpos += c * h->GetBinCenter(ib); }
            }
            if (wsum <= 0.) return 0.f;
            double center = wpos / wsum;
            // find peak bin (max content) in window
            int pb = b0;
            double ph = 0.;
            for (int ib = b0; ib <= b1; ++ib)
                if (h->GetBinContent(ib) > ph) { ph = h->GetBinContent(ib); pb = ib; }
            if (ph <= 0.) return 0.f;
            // walk from peak bin to 50% threshold for fit range
            double thr = 0.2 * ph;
            int lo = pb, hi = pb;
            while (lo > b0 && h->GetBinContent(lo - 1) > thr) --lo;
            while (hi < b1 && h->GetBinContent(hi + 1) > thr) ++hi;
            double fl = h->GetBinLowEdge(lo);
            double fh = h->GetBinLowEdge(hi + 1);
            if (fh - fl < 2. * h->GetBinWidth(pb)) return 0.f;
            // fit-range bracket: solid verticals + dotted cap at 14% pad height
            double bkh = ypad;
            TLine *bkl = new TLine(fl, 0., fl, bkh);
            bkl->SetLineColor(color); bkl->SetLineWidth(2); bkl->Draw();
            TLine *bkr = new TLine(fh, 0., fh, bkh);
            bkr->SetLineColor(color); bkr->SetLineWidth(2); bkr->Draw();
            TLine *bkc = new TLine(fl, bkh, fh, bkh);
            bkc->SetLineColor(color); bkc->SetLineWidth(1); bkc->SetLineStyle(3); bkc->Draw();
            // Gaussian fit
            TF1 *g = new TF1(Form("_gfit_%d", _fit_uid), "gaus", fl, fh);
            g->SetParameters(ph, center, (fh - fl) / 4.);
            h->Fit(g, "RQ0", "", fl, fh);
            g->SetLineColor(color); g->SetLineWidth(3); g->Draw("same");
            // draw a marker at the fitted center, at Gaussian peak height
            double fit_mean = center;
            double fit_amp  = g->GetParameter(0);
            TMarker *mk = new TMarker(fit_mean, fit_amp, 29); // star symbol
            mk->SetMarkerColor(color); mk->SetMarkerSize(2.5); mk->Draw();
            return static_cast<float>(fit_mean);
        };

        // --- Draw both beam-energy histograms on a two-pad canvas, save PNG ---
        TCanvas *ch = new TCanvas(Form("ch_mod_W%d", mod_id-1000),
            Form("Module W%d Histograms", mod_id-1000), 800, 800);
        ch->Divide(1, 2, 0, 0);

        // --- top pad: 3.5 GeV ---
        ch->cd(1);
        gPad->SetBottomMargin(0.005);
        gPad->SetTopMargin(0.10);
        gPad->SetLeftMargin(0.12);
        hist_3p5->GetXaxis()->SetLabelSize(0);
        hist_3p5->GetXaxis()->SetTitleSize(0);
        hist_3p5->SetTitle(Form("Module W%d;  ;Counts", mod_id-1000));
        hist_3p5->SetLineColor(kBlack);
        hist_3p5->SetLineWidth(2);
        hist_3p5->SetStats(0);
        hist_3p5->Draw("HIST");
        float peak_ep_3p5 = fitPeakAndDraw(hist_3p5, e_p_exp_3p5, sigma_ep_3p5, kRed);
        float peak_ee_3p5 = fitPeakAndDraw(hist_3p5, e_e_exp_3p5, sigma_ee_3p5, kBlue);
        {
            TLatex lat;
            lat.SetNDC(); lat.SetTextSize(0.050);
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.50, 0.86, Form("e-p: exp=%.0f  meas=%s",
                (double)e_p_exp_3p5, peak_ep_3p5 > 0.f ? Form("%.0f MeV", (double)peak_ep_3p5) : "N/A"));
            lat.SetTextColor(kBlue);
            lat.DrawLatex(0.50, 0.78, Form("e-e: exp=%.0f  meas=%s",
                (double)e_e_exp_3p5, peak_ee_3p5 > 0.f ? Form("%.0f MeV", (double)peak_ee_3p5) : "N/A"));
            lat.SetTextColor(kBlack);
            lat.DrawLatex(0.15, 0.86, "E_{beam} = 3.5 GeV");
        }

        // --- bottom pad: 0.7 GeV ---
        ch->cd(2);
        gPad->SetTopMargin(0.005);
        gPad->SetBottomMargin(0.14);
        gPad->SetLeftMargin(0.12);
        hist_0p7->SetTitle(";Energy (MeV);Counts");
        hist_0p7->SetLineColor(kBlack);
        hist_0p7->SetLineWidth(2);
        hist_0p7->SetStats(0);
        hist_0p7->Draw("HIST");
        float peak_ep_0p7 = 0.f, peak_ee_0p7 = 0.f;
        peak_ep_0p7 = fitPeakAndDraw(hist_0p7, e_p_exp_0p7, sigma_ep_0p7, kRed);
        if (std::abs(e_p_exp_0p7 - e_e_exp_0p7) > 170.f) {
            peak_ee_0p7 = fitPeakAndDraw(hist_0p7, e_e_exp_0p7, sigma_ee_0p7, kBlue);
        }
        {
            TLatex lat;
            lat.SetNDC(); lat.SetTextSize(0.050);
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.50, 0.86, Form("e-p: exp=%.0f  meas=%s",
                (double)e_p_exp_0p7, peak_ep_0p7 > 0.f ? Form("%.0f MeV", (double)peak_ep_0p7) : "N/A"));
            lat.SetTextColor(kBlue);
            lat.DrawLatex(0.50, 0.78, Form("e-e: exp=%.0f  meas=%s",
                (double)e_e_exp_0p7, peak_ee_0p7 > 0.f ? Form("%.0f MeV", (double)peak_ee_0p7) : "N/A"));
            lat.SetTextColor(kBlack);
            lat.DrawLatex(0.15, 0.86, "E_{beam} = 0.7 GeV");
        }

        ch->SaveAs(Form("%s/mod_W%d.png", pngDir.c_str(), mod_id-1000));
        delete ch;

        // if the anchor point (3.5 GeV e-p) has no clean peak, skip this module
        if (peak_ep_3p5 == 0.f) continue;

        h_energy_peak_3p5->Fill(peak_ep_3p5);

        // make a canvas, measured E vs expected E; only add points with valid peaks
        TCanvas *c = new TCanvas(Form("c_mod_W%d", mod_id-1000), Form("Module W%d Non-linearity", mod_id-1000), 1400, 800);
        c->SetGrid();
        TGraph *g = new TGraph();
        int np = 0;
        auto addPoint = [&](double Eexp, float peak) {
            if (peak != 0.f) g->SetPoint(np++, Eexp, peak);
        };
        float scale = e_p_exp_3p5 / peak_ep_3p5;
        peak_ep_3p5 *= scale; // convert to GeV for fitting
        peak_ee_3p5 *= scale;
        peak_ep_0p7 *= scale;
        peak_ee_0p7 *= scale;
        addPoint(e_p_exp_3p5, peak_ep_3p5);
        addPoint(e_e_exp_3p5, peak_ee_3p5);
        addPoint(e_p_exp_0p7, peak_ep_0p7);
        //if(729. - e_e_exp_0p7 > 6. * 729.*resolution/sqrt(729./1000.f)) 
            addPoint(e_e_exp_0p7, peak_ee_0p7);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetTitle(Form("Module W%d Non-linearity;E_{incident} (MeV);E_{measure} (MeV)", mod_id-1000));
        g->Draw("AP");

        // perfect linearity reference line (y = x)
        double xmin = g->GetXaxis()->GetXmin();
        double xmax = g->GetXaxis()->GetXmax();
        TLine *ref = new TLine(xmin, xmin, xmax, xmax);
        ref->SetLineColor(kRed);
        ref->SetLineStyle(2);
        ref->SetLineWidth(2);
        ref->Draw();

        // fitting: E_incident = E_meas + nl * (E_meas - E_base) /1000.
        // => E_meas = (E_incident + nl/1000. * E_base) / (1 + nl/1000.)
        TF1 *fitLine = new TF1("fitLine",
            [](double *x, double *p){ return (x[0] + p[0] / 1000. * p[1]) / (1.0 + p[0] / 1000.); },
            xmin, xmax, 2);
        fitLine->SetParameter(0, 0.01); // nl ~ small non-linearity
        fitLine->FixParameter(1, e_p_exp_3p5); // E_base fixed at anchor
        fitLine->SetLineColor(kBlue);
        fitLine->SetLineWidth(2);
        g->Fit(fitLine, "Q");
        fitLine->Draw("same");

        double nl        = fitLine->GetParameter(0);
        double nl_err    = fitLine->GetParError(0);
        double chi2      = fitLine->GetChisquare();
        int    ndf       = fitLine->GetNDF();

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.035);
        tex->DrawLatex(0.15, 0.82, "E_{incident} = E_{meas} + nl #times (E_{meas} - E_{base}) / 1000.");
        tex->DrawLatex(0.15, 0.76, Form("nl = %.4f #pm %.4f", nl, nl_err));
        tex->DrawLatex(0.15, 0.70, Form("#chi^{2}/ndf = %.2f / %d = %.3f", chi2, ndf, ndf > 0 ? chi2/ndf : 0.));
        tex->DrawLatex(0.15, 0.64, Form("E_{base} = %.4f GeV", (double)e_p_exp_3p5));
        // corrected points: E_corr = E_meas + nl * (E_meas - E_base) / 1000.
        TGraph *gCorr = new TGraph();
        gCorr->SetMarkerStyle(24); // open circle
        gCorr->SetMarkerSize(1.0);
        gCorr->SetMarkerColor(kGreen+2);
        {
            int nc = 0;
            auto addCorrPoint = [&](double Eexp, float peak) {
                if (peak != 0.f) {
                    double E_corr = peak + nl / 1000. * (peak - e_p_exp_3p5);
                    gCorr->SetPoint(nc++, Eexp, E_corr);
                }
            };
            addCorrPoint(e_p_exp_3p5, peak_ep_3p5);
            addCorrPoint(e_e_exp_3p5, peak_ee_3p5);
            addCorrPoint(e_p_exp_0p7, peak_ep_0p7);
            addCorrPoint(e_e_exp_0p7, peak_ee_0p7);
        }
        gCorr->Draw("P same");

        // legend
        TLegend *leg = new TLegend(0.60, 0.15, 0.88, 0.42);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.032);
        leg->AddEntry(g,       "Measured points",           "p");
        leg->AddEntry(fitLine, "Non-linear fit",            "l");
        leg->AddEntry(ref,     "Perfect linearity (y = x)", "l");
        leg->AddEntry(gCorr,   "Corrected points",          "p");
        leg->Draw();

        c->Write();
        delete leg;
        delete gCorr;
        delete tex;
        delete fitLine;
        delete ref;
        delete c;
        delete g;

        // outermost ring or absorber/beam-hole region: no non-linearity correction
        bool is_outer    = (mod.row == 1-1 || mod.row == 34-1 || mod.column == 1-1 || mod.column == 34-1);
        bool is_absorber = (mod.row >= 16-1 && mod.row <= 19-1 && mod.column >= 16-1 && mod.column <= 19-1);
        if (is_outer || is_absorber) nl = 0.0;
        hycal.SetCalibNonLinearity(mod_id, nl);
    }
    h_energy_peak_3p5->Write();
    std::cout << "Results saved to " << outFile.GetName() << "\n";
    outFile.Close();

    hycal.PrintCalibConstants("new_calibration_NonLinearity.json");

}

void process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal, 
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events)
{   
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if( i % 1000 == 0) {
            std::cerr << "Processing event " << i << "/" << tree->GetEntries() << "\r" << std::flush;
        }
        if ( ev.trigger_bits & ( 1 << 8) == 0) continue;
        if (max_events > 0 && i >= max_events) {
            std::cerr << "Reached max events limit: " << max_events << "\n";
            break;
        }
        
        for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] < 4) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals

            // require hit to be in central 3x3 of a 5x5 grid (|xd|,|yd| < 0.3)
            float xd = (ev.cl_x[j] - (float)mod->x) / (float)mod->size_x;
            float yd = (ev.cl_y[j] - (float)mod->y) / (float)mod->size_y;
            if (std::abs(xd) >= 0.3f || std::abs(yd) >= 0.3f) continue;

            float x = ev.cl_x[j], y = ev.cl_y[j], z = ev.cl_z[j];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            float energy = ev.cl_energy[j];

            bool veto = false;
                float sci_time, sci_int;
                for(int k = 0; k < ev.veto_nch; k++){
                    for(int p = 0; p < ev.veto_npeaks[k]; p++){
                        sci_time = ev.veto_peak_time[k][p];
                        sci_int = ev.veto_peak_integral[k][p];
                        veto = Vetoed(ev.cl_time[j], sci_time, sci_int);
                        if(veto) break;
                    }
                    if(veto) break;
                }
                if(theta > 1.8) veto = false;

            if(veto && energy > 550. && Ebeam < 1000.f) continue;

            energy_hists[mod_id]->Fill(energy);
        }
    }
}