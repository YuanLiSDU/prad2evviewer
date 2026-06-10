
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
// returns the number of events passing the sum-trigger selection
long long process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
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

    std::string output;
    std::vector<std::string> input_3p5, input_0p7;
    std::string pngDir = "module_hists";

    int max_events = -1;
    {
        std::vector<std::string> *cur = nullptr;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-a")       { cur = &input_3p5; }
            else if (arg == "-b")  { cur = &input_0p7; }
            else if (arg == "-o")  { cur = nullptr; if (++i < argc) output = argv[i]; }
            else if (arg == "-n")  { cur = nullptr; if (++i < argc) max_events = std::atoi(argv[i]); }
            else if (arg == "-p")  { cur = nullptr; if (++i < argc) pngDir = argv[i]; }
            else if (!arg.empty() && arg[0] != '-' && cur) { cur->push_back(arg); }
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
    for (const auto &f : input_3p5) {
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

    long long n_sum_3p5 = process_event(tree, ev, hycal, energy_hists_3p5, physics, 3485.f, max_events);

    // --- repeat for 0.7 GeV data ---
    TChain *chain2 = new TChain("recon");
    for (const auto &f : input_0p7) {
        chain2->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
    }
    TTree *tree2 = chain2;
    if (!tree2) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev2;
    prad2::SetReconReadBranches(tree2, ev2);
    long long n_sum_0p7 = process_event(tree2, ev2, hycal, energy_hists_0p7, physics, 729.f, max_events);

    // a file with no sum-trigger events would silently produce an all-default
    // calibration (every module skipped); refuse to write output in that case
    if (n_sum_3p5 <= 0 || n_sum_0p7 <= 0) {
        std::cerr << "No sum-trigger events selected (3.5 GeV: " << n_sum_3p5
                  << ", 0.7 GeV: " << n_sum_0p7
                  << "); check trigger_bits in the inputs. Aborting before writing calibration.\n";
        return 1;
    }

    std::string calib_file = dbDir + "/" + "calibration/calibration_factor_3p5_June1.json";
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    TH1F *h_energy_peak_3p5 = new TH1F("h_energy_peak_3p5", "Energy Peak Distribution;Energy (MeV);Counts", 4000, 0, 4000);

    // calculate non-linearity module by module and save to output file
    gSystem->mkdir(pngDir.c_str(), true);
    TFile outFile(output.empty() ? "nonlinearity_results.root" : output.c_str(), "RECREATE");
    TDirectory *dir_3p5 = outFile.mkdir("energy_3p5GeV");
    TDirectory *dir_0p7 = outFile.mkdir("energy_0p7GeV");
    TDirectory *dir_lin = outFile.mkdir("linearity");
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

        // Peak fit in [Eexp±6σ]: the weighted centroid seeds a two-stage
        // Gaussian fit (mean±2σ, then mean±1σ); returns the final fitted mean.
        auto fitPeakAndDraw = [&](TH1F *h, double Eexp, double sigma, int color) -> float {
            ++_fit_uid;
            int b0 = std::max(1, h->FindBin(Eexp - 6.*sigma));
            int b1 = std::min(h->GetNbinsX(), h->FindBin(Eexp + 6.*sigma));
            double ypad_min = 0;
            double ypad = h->GetMaximum() * 1.2;
            // search-window: dashed, full y-axis height
            TLine *wl = new TLine(Eexp - 6.*sigma, ypad_min, Eexp - 6.*sigma, ypad);
            wl->SetLineColor(6); wl->SetLineStyle(7); wl->SetLineWidth(2); wl->Draw();
            TLine *wr = new TLine(Eexp + 6.*sigma, ypad_min, Eexp + 6.*sigma, ypad);
            wr->SetLineColor(6); wr->SetLineStyle(7); wr->SetLineWidth(2); wr->Draw();
            // Step 2: weighted mean within search window
            double wsum = 0., wpos = 0.;
            for (int ib = b0; ib <= b1; ++ib) {
                double c = h->GetBinContent(ib);
                if (c > 0.) { wsum += c; wpos += c * h->GetBinCenter(ib); }
            }
            if (wsum <= 0.) return 0.f;
            double mean = wpos / wsum;
            // Step 3: first Gaussian fit within [mean ± 2σ]
            auto estimateSig = [](double E) -> double {
                return (E > 0.) ? E * 0.035 / std::sqrt(E / 1000.) : 1.;
            };
            double sig = estimateSig(mean);
            {
                TF1 g1(Form("_fpk_g1_%d", _fit_uid), "gaus", mean - 2.*sig, mean + 2.*sig);
                g1.SetParameters(h->GetMaximum(), mean, sig);
                h->Fit(&g1, "RQ0");
                mean = g1.GetParameter(1);
            }
            // Step 4: final Gaussian fit within [mean ± 1σ]
            sig = estimateSig(mean);
            TF1 *g = new TF1(Form("_gfit_%d", _fit_uid), "gaus", mean - sig, mean + sig);
            g->SetParameters(h->GetMaximum(), mean, sig);
            h->Fit(g, "RQ0");
            mean = g->GetParameter(1);
            double fit_amp = g->GetParameter(0);
            g->SetLineColor(color); g->SetLineWidth(3); g->Draw("same");
            // draw a marker at the fitted mean, at Gaussian peak height
            TMarker *mk = new TMarker(mean, fit_amp, 29); // star symbol
            mk->SetMarkerColor(color); mk->SetMarkerSize(2.5); mk->Draw();
            return static_cast<float>(mean);
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

        // make a canvas, E_rec/E_exp vs E_rec; only add points with valid peaks
        TCanvas *c = new TCanvas(Form("c_mod_W%d", mod_id-1000), Form("Module W%d Non-linearity", mod_id-1000), 1400, 800);
        c->SetGrid();
        TGraph *g = new TGraph();
        int np = 0;
        auto addPoint = [&](double Eexp, float peak) {
            if (peak != 0.f) g->SetPoint(np++, peak, peak / Eexp);
        };
        float scale = e_p_exp_3p5 / peak_ep_3p5;
        //peak_ep_3p5 *= scale;
        //peak_ee_3p5 *= scale;
        //peak_ep_0p7 *= scale;
        //peak_ee_0p7 *= scale;
        addPoint(e_p_exp_3p5, peak_ep_3p5);
        addPoint(e_e_exp_3p5, peak_ee_3p5);
        addPoint(e_p_exp_0p7, peak_ep_0p7);
        addPoint(e_e_exp_0p7, peak_ee_0p7);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetTitle(Form("Module W%d Non-linearity;E_{rec} (MeV);E_{rec}/E_{exp}", mod_id-1000));
        g->Draw("AP");
        g->GetYaxis()->SetRangeUser(0.85, 1.15);
        g->Draw("AP");

        // perfect linearity reference line (y = 1)
        double xmin = g->GetXaxis()->GetXmin();
        double xmax = g->GetXaxis()->GetXmax();
        TLine *ref = new TLine(xmin, 1.0, xmax, 1.0);
        ref->SetLineColor(kRed);
        ref->SetLineStyle(2);
        ref->SetLineWidth(2);
        ref->Draw();

        // 1st order fit: E_rec/E_exp = 1 + nl1 * (E_rec - E_base)/1000
        TF1 *fitLine = new TF1("fitLine",
            [](double *x, double *p){ return 1.0 + p[0] * (x[0] - p[1])/1000.0; },
            xmin, xmax, 2);
        fitLine->SetParameter(0, 0.01);
        fitLine->FixParameter(1, e_p_exp_3p5);
        fitLine->SetLineColor(kBlue);
        fitLine->SetLineWidth(2);
        g->Fit(fitLine, "RQ0");
        fitLine->Draw("same");

        double nl        = fitLine->GetParameter(0);
        double nl_err    = fitLine->GetParError(0);
        double chi2      = fitLine->GetChisquare();
        int    ndf       = fitLine->GetNDF();

        // 2nd order fit: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000 + nl2*((E_rec-E_base)/1000)^2
        TF1 *fitLine2 = new TF1("fitLine2",
            [](double *x, double *p){
                double t = (x[0] - p[2]) / 1000.0;
                return 1.0 + p[0] * t + p[1] * t * t;
            },
            xmin, xmax, 3);
        fitLine2->SetParameter(0, nl);
        fitLine2->SetParameter(1, 0.0);
        fitLine2->FixParameter(2, e_p_exp_3p5);
        fitLine2->SetLineColor(kMagenta+1);
        fitLine2->SetLineWidth(2);
        fitLine2->SetLineStyle(7);
        g->Fit(fitLine2, "RQ0+");
        fitLine2->Draw("same");

        double nl2_1     = fitLine2->GetParameter(0);
        double nl2_1_err = fitLine2->GetParError(0);
        double nl2_2     = fitLine2->GetParameter(1);
        double nl2_2_err = fitLine2->GetParError(1);
        double chi2_2    = fitLine2->GetChisquare();
        int    ndf_2     = fitLine2->GetNDF();

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.030);
        tex->SetTextColor(kBlue);
        tex->DrawLatex(0.15, 0.88, "1st: E_{rec}/E_{exp} = 1 + nl_{1} #times (E_{rec}-E_{base})/1000");
        tex->DrawLatex(0.15, 0.83, Form("nl_{1} = %.4f #pm %.4f, #chi^{2}/ndf = %.2f/%d", nl, nl_err, chi2, ndf));
        tex->SetTextColor(kMagenta+1);
        tex->DrawLatex(0.15, 0.78, "2nd: + nl_{2} #times ((E_{rec}-E_{base})/1000)^{2}");
        tex->DrawLatex(0.15, 0.73, Form("nl_{1} = %.4f #pm %.4f, nl_{2} = %.4f #pm %.4f", nl2_1, nl2_1_err, nl2_2, nl2_2_err));
        tex->DrawLatex(0.15, 0.68, Form("#chi^{2}/ndf = %.2f/%d,  E_{base} = %.1f MeV", chi2_2, ndf_2, (double)e_p_exp_3p5));
        tex->SetTextColor(kBlack);

        // corrected points using 1st order: E_corr = E_rec / (1 + nl*(E_rec-E_base)/1000)
        TGraph *gCorr = new TGraph();
        gCorr->SetMarkerStyle(24);
        gCorr->SetMarkerSize(1.0);
        gCorr->SetMarkerColor(kGreen+2);
        {
            int nc = 0;
            auto addCorrPoint = [&](double Eexp, float peak) {
                if (peak != 0.f) {
                    double denom = 1.0 + nl * (peak - e_p_exp_3p5)/1000.0;
                    double E_corr = (denom != 0.) ? peak / denom : peak;
                    gCorr->SetPoint(nc++, E_corr, E_corr / Eexp);
                }
            };
            addCorrPoint(e_p_exp_3p5, peak_ep_3p5);
            addCorrPoint(e_e_exp_3p5, peak_ee_3p5);
            addCorrPoint(e_p_exp_0p7, peak_ep_0p7);
            addCorrPoint(e_e_exp_0p7, peak_ee_0p7);
        }
        gCorr->Draw("P same");

        // corrected points using 2nd order fit
        TGraph *gCorr2 = new TGraph();
        gCorr2->SetMarkerStyle(25);
        gCorr2->SetMarkerSize(1.0);
        gCorr2->SetMarkerColor(kOrange+2);
        {
            int nc = 0;
            auto addCorrPoint2 = [&](double Eexp, float peak) {
                if (peak != 0.f) {
                    double t = (peak - e_p_exp_3p5) / 1000.0;
                    double denom = 1.0 + nl2_1 * t + nl2_2 * t * t;
                    double E_corr = (denom != 0.) ? peak / denom : peak;
                    gCorr2->SetPoint(nc++, E_corr, E_corr / Eexp);
                }
            };
            addCorrPoint2(e_p_exp_3p5, peak_ep_3p5);
            addCorrPoint2(e_e_exp_3p5, peak_ee_3p5);
            addCorrPoint2(e_p_exp_0p7, peak_ep_0p7);
            addCorrPoint2(e_e_exp_0p7, peak_ee_0p7);
        }
        gCorr2->Draw("P same");

        // legend
        TLegend *leg = new TLegend(0.60, 0.12, 0.92, 0.52);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.026);
        leg->AddEntry(g,        "Measured points",         "p");
        leg->AddEntry(fitLine,  "1st order fit",           "l");
        leg->AddEntry(fitLine2, "2nd order fit",           "l");
        leg->AddEntry(ref,      "Perfect linearity (y=1)", "l");
        leg->AddEntry(gCorr,    "Corrected points (1st)",  "p");
        leg->AddEntry(gCorr2,   "Corrected points (2nd)",  "p");
        leg->Draw();

        dir_lin->cd();
        c->Write();
        outFile.cd();
        delete leg;
        delete gCorr2;
        delete gCorr;
        delete tex;
        delete fitLine2;
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
    dir_3p5->cd();
    for (auto &[id, h] : energy_hists_3p5) if (h) h->Write();
    dir_0p7->cd();
    for (auto &[id, h] : energy_hists_0p7) if (h) h->Write();
    outFile.cd();
    h_energy_peak_3p5->Write();
    std::cout << "Results saved to " << outFile.GetName() << "\n";
    outFile.Close();

    hycal.PrintCalibConstants("new_calibration_NonLinearity.json");

}

long long process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events)
{
    long long n_accepted = 0;
    for (int i = 0; i < tree->GetEntries(); i++) {
        // entry-count cap; checked before the trigger cut so -n N stops at entry N
        if (max_events > 0 && i >= max_events) {
            std::cerr << "Reached max events limit: " << max_events << "\n";
            break;
        }
        tree->GetEntry(i);
        if( i % 1000 == 0) {
            std::cerr << "Processing event " << i << "/" << tree->GetEntries() << "\r" << std::flush;
        }
        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        n_accepted++;

        for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] < 4) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals

            // require hit to be in central 3x3 of a 5x5 grid (|xd|,|yd| < 0.3)
            float xd = (ev.cl_x[j] - (float)mod->x) / (float)mod->size_x;
            float yd = (ev.cl_y[j] - (float)mod->y) / (float)mod->size_y;
            //if (std::abs(xd) >= 0.25f || std::abs(yd) >= 0.25f) continue;

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
                if(theta > 1.3) veto = false;

            if(veto && energy > 600. && Ebeam < 1000.f) continue;

            energy_hists[mod_id]->Fill(energy);
        }
    }
    return n_accepted;
}
