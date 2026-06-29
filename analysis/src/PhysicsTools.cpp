//=============================================================================
// PhysicsTools.cpp — physics analysis tools
//=============================================================================

#include "PhysicsTools.h"
#include "InstallPaths.h"
#include <TF1.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace analysis {

// Physical constants
static constexpr float M_PROTON  = 938.272f;   // MeV
static constexpr float M_ELECTRON = 0.511f;    // MeV
static constexpr float DEG2RAD = 3.14159265f / 180.f;

PhysicsTools::PhysicsTools(fdec::HyCalSystem &hycal)
    : hycal_(hycal)
{
    int nmod = hycal_.module_count();
    module_hists_.resize(nmod);
    for (int i = 0; i < nmod; ++i) {
        auto &mod = hycal_.module(i);
        std::string name = "h_" + mod.name;
        std::string title = mod.name + " cluster energy;Energy (MeV);Counts";
        module_hists_[i] = std::make_unique<TH1F>(name.c_str(), title.c_str(), 500, 0, 5000);
    }
    h2_energy_module_ = std::make_unique<TH2F>(
        "h2_energy_module", "Energy vs Module;Module Index;Energy (MeV)",
        nmod, 0, nmod, 2000, 0, 4000);
    h2_energy_theta_ = std::make_unique<TH2F>(
        "h2_energy_theta", "Energy vs Theta;Theta (deg);Energy (MeV)",
        160, 0, 8, 4000, 0, 4000);
    h2_Nevents_moduleMap_ = std::make_unique<TH2F>(
        "h2_Nevents_moduleMap", "Number of Events per Module;Column;Row",
        34, 0.5, 34.5, 34, -34.5, -0.5);

    h2_moller_pos_ = std::make_unique<TH2F>(
        "h2_moller_pos", "Moller 2-arm Hit Position;X (mm);Y (mm)",
        200, -500, 500, 200, -500, 500);

    moller_phi_diff_ = std::make_unique<TH1F>(
        "h_moller_phi_diff", "Moller Phi Difference;Phi_{e1} - Phi_{e2} (deg);Counts",
        40, -20, 20);
    moller_x_ = std::make_unique<TH1F>(
        "h_moller_x", "Moller Center X Position (HyCal);X (mm);Counts",
        100, -10, 10);
    moller_y_ = std::make_unique<TH1F>(
        "h_moller_y", "Moller Center Y Position (HyCal);Y (mm);Counts",
        100, -10, 10);
    moller_z_ = std::make_unique<TH1F>(
        "h_moller_z", "Moller Z Position (HyCal);Z (mm);Counts",
        1000, 5000, 8000);

    // histograms for gain monitoring replay
    for (int ch = 0; ch < 4; ++ch) {
        h_lmsCH_lmsHeight_[ch] = std::make_unique<TH1F>(
            Form("h_lmsCH%d_lmsHeight", ch), Form("LMS%d Peak Height;Height (ADC);Counts", ch), 1000, 0, 4000);
        h_lmsCH_lmsIntegral_[ch] = std::make_unique<TH1F>(
            Form("h_lmsCH%d_lmsIntegral", ch), Form("LMS%d Peak Integral;Integral (ADC*ns);Counts", ch), 1000, 0, 40000);
        h_lmsCH_alphaHeight_[ch] = std::make_unique<TH1F>(
            Form("h_lmsCH%d_alphaHeight", ch), Form("LMS%d Alpha Peak Height;Height (ADC);Counts", ch), 1000, 0, 4000);
        h_lmsCH_alphaIntegral_[ch] = std::make_unique<TH1F>(
            Form("h_lmsCH%d_alphaIntegral", ch), Form("LMS%d Alpha Peak Integral;Integral (ADC*ns);Counts", ch), 1000, 0, 400000);
    }
    module_gains_.resize(nmod);
    for (int i = 0; i < nmod; ++i)
        module_gains_[i].name = hycal_.module(i).name;

    h_modCH_lmsHeight_.resize(nmod);
    h_modCH_lmsIntegral_.resize(nmod);
    for (int i = 0; i < nmod; ++i) {
        auto &mod = hycal_.module(i);
        std::string name_height = "h_mod" + mod.name + "_lmsHeight";
        std::string title_height = mod.name + " LMS Peak Height;Height (ADC);Counts";
        h_modCH_lmsHeight_[i] = std::make_unique<TH1F>(name_height.c_str(), title_height.c_str(), 1000, 0, 4000);
        std::string name_integral = "h_mod" + mod.name + "_lmsIntegral";
        std::string title_integral = mod.name + " LMS Peak Integral;Integral (ADC*ns);Counts";
        h_modCH_lmsIntegral_[i] = std::make_unique<TH1F>(name_integral.c_str(), title_integral.c_str(), 1000, 0, 40000);
    }
}

PhysicsTools::~PhysicsTools() = default;

// Crystal Ball: p[0]=amp, p[1]=mean, p[2]=sigma, p[3]=alpha, p[4]=n
static double crystalBallFunc(double *x, double *p)
{
    double amp   = p[0];
    double mu    = p[1];
    double sigma = p[2];
    double alpha = p[3];
    double n     = p[4];
    double t = (x[0] - mu) / sigma;
    if (t > -std::abs(alpha)) {
        return amp * std::exp(-0.5 * t * t);
    } else {
        double a = std::pow(n / std::abs(alpha), n) * std::exp(-0.5 * alpha * alpha);
        double b = n / std::abs(alpha) - std::abs(alpha);
        return amp * a * std::pow(b - t, -n);
    }
}

void PhysicsTools::FillModuleEnergy(int module_id, float energy)
{   
    if (module_id >= 0){
        int module_index = hycal_.id_to_index(module_id);
        if (module_index >= 0 && module_index < (int)module_hists_.size())
            module_hists_[module_index]->Fill(energy);
    }
}

TH1F *PhysicsTools::GetModuleEnergyHist(int module_id) const
{
    int module_index = hycal_.id_to_index(module_id);
    if (module_index >= 0 && module_index < (int)module_hists_.size())
        return module_hists_[module_index].get();
    return nullptr;
}

void PhysicsTools::FillEnergyVsModule(int module_id, float energy)
{
    int module_index = hycal_.id_to_index(module_id);
    if (module_index >= 0 && module_index < (int)module_hists_.size())
        h2_energy_module_->Fill(module_index, energy);
}

void PhysicsTools::FillEnergyVsTheta(float theta_deg, float energy)
{
    if (h2_energy_theta_)
        h2_energy_theta_->Fill(theta_deg, energy);
}

std::unique_ptr<TH1F> PhysicsTools::GetEpYieldHist(TH2F *energy_theta, float Ebeam)
{
    if (!energy_theta) return nullptr;

    auto h_ep = std::make_unique<TH1F>("ep_yield", "Elastic e-p Yield;Scattering Angle (deg);Counts", 80, 0, 8);
    h_ep->SetDirectory(nullptr);
    for (int i = 1; i <= energy_theta->GetNbinsX(); i++) {
        for (int j = 1; j <= energy_theta->GetNbinsY(); j++) {
            float theta = energy_theta->GetXaxis()->GetBinCenter(i);
            float E = energy_theta->GetYaxis()->GetBinCenter(j);
            float E_expected = ExpectedEnergy(theta, Ebeam, "ep");
            if (std::abs(E - E_expected) < E_expected*0.026f/std::sqrt(E_expected/1000.f)) {
                float count = energy_theta->GetBinContent(i, j);
                h_ep->Fill(theta, count);
            }
        }
    }
    return h_ep;
}

std::unique_ptr<TH1F> PhysicsTools::GetEeYieldHist(TH2F *energy_theta, float Ebeam)
{
    if (!energy_theta) return nullptr;

    auto h_ee = std::make_unique<TH1F>("ee_yield", "Elastic e-e Yield;Scattering Angle (deg);Counts", 80, 0, 8);
    h_ee->SetDirectory(nullptr);
    for (int i = 1; i <= energy_theta->GetNbinsX(); i++) {
        for (int j = 1; j <= energy_theta->GetNbinsY(); j++) {
            float theta = energy_theta->GetXaxis()->GetBinCenter(i);
            float E = energy_theta->GetYaxis()->GetBinCenter(j);
            float E_expected = ExpectedEnergy(theta, Ebeam, "ee");
            if (std::abs(E - E_expected) < E_expected*0.026f/std::sqrt(E_expected/1000.f)) {
                float count = energy_theta->GetBinContent(i, j);
                h_ee->Fill(theta, count);
            }
        }
    }
    return h_ee;
}

std::unique_ptr<TH1F> PhysicsTools::GetYieldRatioHist(TH1F *ep_hist, TH1F *ee_hist)
{
    if (!ep_hist || !ee_hist) return nullptr;

    auto h_ratio = std::make_unique<TH1F>("yield_ratio", "Yield Ratio (e-p / e-e);Scattering Angle (deg);Ratio", 80, 0, 8);
    h_ratio->SetDirectory(nullptr);
    for (int i = 1; i <= ep_hist->GetNbinsX(); i++) {
        float theta = ep_hist->GetXaxis()->GetBinCenter(i);
        float ep_count = ep_hist->GetBinContent(i);
        float ee_count = ee_hist->GetBinContent(i);
        if (ee_count > 0) {
            h_ratio->Fill(theta, ep_count / ee_count);
        }
    }
    return h_ratio;
}

void PhysicsTools::Fill2armMollerPosHist(float x, float y)
{
    if (h2_moller_pos_)
        h2_moller_pos_->Fill(x, y);
}

std::array<float, 3> PhysicsTools::FitPeakResolution(int module_id) const
{
    int module_index = hycal_.id_to_index(module_id);
    if (module_index < 0 || module_index >= (int)module_hists_.size())
        return {0.f, 0.f, 0.f};

    TH1F *h = module_hists_[module_index].get();
    if (!h || h->GetEntries() < 1) return {0.f, 0.f, 100.f};

    const double resolution = 0.035; // 3.5% / sqrt(E/1000) energy resolution

    // estimate sigma from energy resolution: sigma = E * resolution / sqrt(E/1000)
    auto estimateSigma = [&](double E) -> double {
        return (E > 0.) ? E * resolution / std::sqrt(E / 1000.) : 1.;
    };

    // Step 1: find histogram maximum as initial center
    double center = h->GetBinCenter(h->GetMaximumBin());

    // Step 2: weighted mean within [center ± 3σ]
    double sigma = estimateSigma(center);
    double sumW = 0., sumWx = 0.;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double x = h->GetBinCenter(b);
        if (x < center - 3.*sigma || x > center + 3.*sigma) continue;
        double w = h->GetBinContent(b);
        sumW  += w;
        sumWx += w * x;
    }
    double mean = (sumW > 0.) ? sumWx / sumW : center;

    // Step 3: first Gaussian fit within [mean ± 2σ], σ re-estimated from new mean
    sigma = estimateSigma(mean);
    {
        TF1 g1("_fpk_g1_", "gaus", mean - 2.*sigma, mean + 2.*sigma);
        g1.SetParameters(h->GetMaximum(), mean, sigma);
        h->Fit(&g1, "RQ0");
        mean  = g1.GetParameter(1);
        // keep sigma from resolution for next iteration boundary
    }

    // Step 4: final Gaussian fit within [mean ± 1σ], σ re-estimated from new mean
    sigma = estimateSigma(mean);
    TF1 g2("_fpk_g2_", "gaus", mean - sigma, mean + sigma);
    g2.SetParameters(h->GetMaximum(), mean, sigma);
    h->Fit(&g2, "RQ0");
    mean  = g2.GetParameter(1);
    sigma = std::abs(g2.GetParameter(2));
    double chi2 = (g2.GetNDF() > 0) ? g2.GetChisquare() / g2.GetNDF() : 0.;

    return {static_cast<float>(mean), static_cast<float>(sigma), static_cast<float>(chi2)};
}

void PhysicsTools::Resolution2Database(int run_id)
{
    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    std::string filename = db_dir + Form("/recon/run_%d.dat", run_id);

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    int module_count = hycal_.module_count();
    for (int m = 0; m < module_count; m++) {
        int module_id = hycal_.module(m).id;
        auto [peak, sigma, chi2] = FitPeakResolution(module_id);
        if (peak > 0 && sigma > 0) {
            std::string name = hycal_.module(m).name;
            out << name << " " << peak << " " << sigma << " " << chi2 << "\n";
        }
    }
}

float PhysicsTools::ExpectedEnergy(float theta_deg, float Ebeam, const std::string &type)
{
    float theta = theta_deg * DEG2RAD;
    float cos_t = std::cos(theta);
    float sin_t = std::sin(theta);

    if (type == "ep") {
        // elastic e-p: E' = E * M / (M + E*(1 - cos_t))
        // where M = proton mass
        float expectE = Ebeam * M_PROTON / (M_PROTON + Ebeam * (1.f - cos_t));
        float eloss = EnergyLoss(theta_deg, expectE);
        return expectE - eloss;
    }
    if (type == "ee") {
        // Moller scattering: exact lab-frame formula from 4-momentum conservation
        // E' = m * [(gamma+1) + (gamma-1)*cos^2(theta)] / [(gamma+1) - (gamma-1)*cos^2(theta)]
        float gamma = Ebeam / M_ELECTRON;
        float num = (gamma + 1.f) + (gamma - 1.f) * cos_t * cos_t;
        float den = (gamma + 1.f) - (gamma - 1.f) * cos_t * cos_t;
        if (den <= 0) return 0.f;
        float expectE = M_ELECTRON * num / den;
        float eloss = EnergyLoss(theta_deg, expectE);
        return expectE - eloss;
    }
    return 0.f;
}

float PhysicsTools::EnergyLoss(float theta_deg, float E)
{
    // simplified energy loss through target materials
    // path lengths scale as 1/cos(theta) for small angles
    float theta = theta_deg * DEG2RAD;
    float sec = (std::cos(theta) > 0.01f) ? (1.f / std::cos(theta)) : 100.f;

    // material thicknesses (mm) and dE/dx (MeV/mm) — approximate values
    // aluminum window: 0.5 mm, dE/dx ~ 1.6 MeV/mm, vacuum window
    // GEM foils: ~0.05 mm effective per GEM, dE/dx ~ 2.0 MeV/mm
    // GEM win Al foils: ~0.06 mm effective per GEM, dE/dx ~ 1.6 MeV/mm
    // kapton window: ~0.24 mm per GEM, dE/dx ~ 1.8 MeV/mm,
    float eloss = 0.f;
    eloss += 0.500f * 1.6f * sec;  // Al window
    eloss += 0.120f * 1.6f * sec;  // GEM win Al foils (2 GEMs)
    eloss += 0.100f * 2.0f * sec;  // GEM foils (2 GEMs)
    eloss += 0.480f * 1.8f * sec;  // kapton cover

    return eloss;  // total energy loss in MeV
}

bool PhysicsTools::isMoller_kinematic(float theta_deg1, float energy1, float theta_deg2, float energy2, float EBeam, float resolution)
{
    float expectE1 = ExpectedEnergy(theta_deg1, EBeam, "ee");
    float expectE2 = ExpectedEnergy(theta_deg2, EBeam, "ee");

    bool E_sum = false, E1_ok = false, E2_ok = false, phi_ok = false;

    if(fabs(energy1 + energy2 - EBeam) < 5.f * resolution * EBeam / sqrt(EBeam/1000.f)) 
        E_sum = true;
    if(fabs(energy1 - expectE1) < 3.5f * expectE1 * resolution / sqrt(expectE1/1000.f)) 
        E1_ok = true;
    if(fabs(energy2 - expectE2) < 3.5f * expectE2 * resolution / sqrt(expectE2/1000.f)) 
        E2_ok = true;

    return E_sum && E1_ok && E2_ok;
}

std::array<float, 2> PhysicsTools::GetMollerCenter(MollerEvent &event1, MollerEvent &event2)
{
    float x1[2], y1[2];
    float x2[2], y2[2];

    x1[0] = event1.first.x; y1[0] = event1.first.y;
    x1[1] = event1.second.x; y1[1] = event1.second.y;
    x2[0] = event2.first.x; y2[0] = event2.first.y;
    x2[1] = event2.second.x; y2[1] = event2.second.y;

    //two lines: y = ax + b, y = cx + d
    float dx1 = x1[0] - x1[1];
    float dx2 = x2[0] - x2[1];
    if (std::abs(dx1) < 1e-6f || std::abs(dx2) < 1e-6f)
        return {0.f, 0.f};  // vertical line — degenerate

    float a = (y1[0] - y1[1]) / dx1;
    float b = y1[0] - a * x1[0];
    float c = (y2[0] - y2[1]) / dx2;
    float d = y2[0] - c * x2[0];

    if (std::abs(a - c) < 1e-6f)
        return {0.f, 0.f};  // parallel lines — no intersection

    float x_cross = (d - b) / (a - c);
    float y_cross = a * x_cross + b;

    return {x_cross, y_cross};

}

float PhysicsTools::GetMollerZdistance(MollerEvent &event, float Ebeam)
{
    float R1 = sqrt(event.first.x*event.first.x + event.first.y*event.first.y);
    float R2 = sqrt(event.second.x*event.second.x + event.second.y*event.second.y);
    float z = sqrt( (Ebeam + M_ELECTRON) * R1 * R2 / (2.*M_ELECTRON) );
    return z;
}

float PhysicsTools::GetMollerPhiDiff(MollerEvent &event1)
{
    // Calculate the azimuthal angle difference (phi) for a Moller event
    float x1 = event1.first.x, y1 = event1.first.y;
    float x2 = event1.second.x, y2 = event1.second.y;
    float phi1 = GetPhiAngle(x1, y1);
    float phi2 = GetPhiAngle(x2, y2);
    float phi_diff = fabs(phi1 - phi2) - 180.f; // Expecting back-to-back, so difference should be around 180 degrees
    return phi_diff;
}

float PhysicsTools::GetPhiAngle(float x, float y)
{
    // atan2 handles all quadrants and x==0 correctly
    float phi = std::atan2(y, x) * 180.f / static_cast<float>(TMath::Pi());
    if (phi < 0) phi += 360.f;
    return phi;
}

//for gain factor monitoring

// Helper: fit a TH1F with a Gaussian and return {mean, sigma, chi2/ndf}.
// Returns {0,0,0} if histogram is null or has too few entries.
static std::array<double, 3> fitGaus(TH1F *h)
{
    if (!h || h->GetEntries() < 10) return {0., 0., 0.};
    double peak0 = h->GetBinCenter(h->GetMaximumBin());
    double rms0  = h->GetRMS();
    if (rms0 <= 0.) rms0 = peak0 * 0.1;
    double lo = peak0 - 2. * rms0, hi = peak0 + 2. * rms0;
    TF1 gaus("_fg_", "gaus", lo, hi);
    gaus.SetParameters(h->GetMaximum(), peak0, rms0);
    h->Fit(&gaus, "RQ0");
    double chi2 = (gaus.GetNDF() > 0) ? gaus.GetChisquare() / gaus.GetNDF() : 0.;
    return {gaus.GetParameter(1), std::abs(gaus.GetParameter(2)), chi2};
}

void PhysicsTools::ComputeModuleGains()
{
    // --- fit LMS and alpha reference channels (index 1..3) ---
    double lms_ref[4]   = {}, alpha_ref[4] = {};
    for (int i = 1; i <= 3; ++i) {
        auto r = fitGaus(h_lmsCH_lmsIntegral_[i].get());
        lms_ref[i] = r[0];
        auto a = fitGaus(h_lmsCH_alphaIntegral_[i].get());
        alpha_ref[i] = a[0];
    }

    // --- update module_gains_ in-place ---
    int nmod = hycal_.module_count();
    for (int i = 0; i < nmod; ++i) {
        auto &mod = hycal_.module(i);
        auto &res = module_gains_[i];
        // reset numeric fields
        res.lms_peak = res.lms_sigma = res.lms_chi2 = 0.f;
        res.g[0] = res.g[1] = res.g[2] = res.g[3] = 0.f;

        if (!mod.is_hycal()) continue;
        if (mod.name.empty() || mod.name[0] != 'W') continue;

        TH1F *h = (i < (int)h_modCH_lmsIntegral_.size())
                  ? h_modCH_lmsIntegral_[i].get() : nullptr;
        auto r = fitGaus(h);
        if (r[0] <= 0.) continue;

        res.lms_peak  = static_cast<float>(r[0]);
        res.lms_sigma = static_cast<float>(r[1]);
        res.lms_chi2  = static_cast<float>(r[2]);
        for (int j = 1; j <= 3; ++j) {
            res.g[j] = (lms_ref[j] > 0. && alpha_ref[j] > 0.)
                       ? static_cast<float>(r[0] * alpha_ref[j] / lms_ref[j])
                       : 1.f;
        }
    }

    // --- reset histograms for next fill cycle ---
    for (int i = 1; i <= 3; ++i) {
        if (h_lmsCH_lmsIntegral_[i])   h_lmsCH_lmsIntegral_[i]->Reset();
        if (h_lmsCH_alphaIntegral_[i]) h_lmsCH_alphaIntegral_[i]->Reset();
    }
    for (auto &h : h_modCH_lmsIntegral_)
        if (h) h->Reset();
}

}
