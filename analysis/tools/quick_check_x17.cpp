// quick_check.C — ROOT script version 
//
// Reads reconstructed ROOT tree (output of replay_recon), runs physics
// analysis using PhysicsTools and MatchingTools from prad2det, and saves
// histograms to an output ROOT file.
// Usage:
//   quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]
//   -o  output ROOT file (default: input filename with _quick_check.root suffix)
//   -n  max events to process (default: all)
//   -j  number of input-file worker threads (default: 4)
// Example:
//   quick_check recon.root -o recon_check.root -n 10000
//   quick_check recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <limits>
#include <thread>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);
static std::string makeDefaultOutput(const std::string &input_path);

bool inHyCal(float xmm, float ymm) {
    const float module = 20.75; // mm
    return (fabs(xmm) > module * 2. || fabs(ymm) > module * 2.)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
}

static float electronPairInvariantMass(
    float x1, float y1, float z1, float E1,
    float x2, float y2, float z2, float E2)
{
    constexpr float electron_mass = 0.51099895f; // MeV
    const float position[2][3] = {{x1, y1, z1}, {x2, y2, z2}};
    const float energy[2] = {E1, E2};
    float momentum[2][3] = {};

    for (int i = 0; i < 2; ++i) {
        const float norm = std::sqrt(
            position[i][0] * position[i][0]
            + position[i][1] * position[i][1]
            + position[i][2] * position[i][2]);
        if (energy[i] < electron_mass || norm <= 0.f)
            return std::numeric_limits<float>::quiet_NaN();
        const float p = std::sqrt(std::max(
            0.f, energy[i] * energy[i] - electron_mass * electron_mass));
        momentum[i][0] = p * position[i][0] / norm;
        momentum[i][1] = p * position[i][1] / norm;
        momentum[i][2] = p * position[i][2] / norm;
    }

    const float total_energy = energy[0] + energy[1];
    const float px = momentum[0][0] + momentum[1][0];
    const float py = momentum[0][1] + momentum[1][1];
    const float pz = momentum[0][2] + momentum[1][2];
    const float mass2 = total_energy * total_energy - px * px - py * py - pz * pz;
    return std::sqrt(std::max(0.f, mass2));
}

struct GEMTrackPair {
    float xu, yu, zu;
    float xd, yd, zd;
};

static bool solveGEMVertex(const GEMTrackPair tracks[3], std::array<float, 3> &vertex)
{
    float A[3][3] = {};
    float b[3] = {};

    for (int i = 0; i < 3; ++i) {
        const float dx = tracks[i].xd - tracks[i].xu;
        const float dy = tracks[i].yd - tracks[i].yu;
        const float dz = tracks[i].zd - tracks[i].zu;
        const float norm = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (norm <= 0.f)
            return false;

        const float ux = dx / norm;
        const float uy = dy / norm;
        const float uz = dz / norm;
        const float p0[3] = {tracks[i].xu, tracks[i].yu, tracks[i].zu};

        const float P[3][3] = {
            {1.f - ux * ux,     -ux * uy,        -ux * uz},
            {    -uy * ux, 1.f - uy * uy,        -uy * uz},
            {    -uz * ux,     -uz * uy,    1.f - uz * uz}
        };

        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c)
                A[r][c] += P[r][c];
            b[r] += P[r][0] * p0[0] + P[r][1] * p0[1] + P[r][2] * p0[2];
        }
    }

    float M[3][4] = {
        {A[0][0], A[0][1], A[0][2], b[0]},
        {A[1][0], A[1][1], A[1][2], b[1]},
        {A[2][0], A[2][1], A[2][2], b[2]}
    };

    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 3; ++row) {
            if (std::fabs(M[row][col]) > std::fabs(M[pivot][col]))
                pivot = row;
        }
        if (std::fabs(M[pivot][col]) < 1e-6f)
            return false;
        if (pivot != col) {
            for (int k = col; k < 4; ++k)
                std::swap(M[col][k], M[pivot][k]);
        }

        const float inv = 1.f / M[col][col];
        for (int k = col; k < 4; ++k)
            M[col][k] *= inv;

        for (int row = 0; row < 3; ++row) {
            if (row == col)
                continue;
            const float factor = M[row][col];
            for (int k = col; k < 4; ++k)
                M[row][k] -= factor * M[col][k];
        }
    }

    vertex = {M[0][3], M[1][3], M[2][3]};
    return std::isfinite(vertex[0]) && std::isfinite(vertex[1]) && std::isfinite(vertex[2]);
}

static bool isMollerAngleEnergyMatch(const PhysicsTools &physics,
                                     float theta,
                                     float energy,
                                     float ebeam,
                                     float nsigma = 3.5f)
{
    const float expectE = physics.ExpectedEnergy(theta, ebeam, "ee");
    if (expectE <= 0.f)
        return false;
    return std::fabs(energy - expectE)
           < nsigma * 0.031f * expectE / std::sqrt(expectE / 1000.f);
}

static bool isMollerPairCandidate(const PhysicsTools &physics,
                                  float xa, float ya, float za, float Ea, float theta_a,
                                  float xb, float yb, float zb, float Eb, float theta_b,
                                  float ebeam)
{
    if (!isMollerAngleEnergyMatch(physics, theta_a, Ea, ebeam) ||
        !isMollerAngleEnergyMatch(physics, theta_b, Eb, ebeam)) {
        return false;
    }

    const float sigma_pair = ebeam * 0.035f / std::sqrt(ebeam / 1000.f);
    if (std::fabs((Ea + Eb) - ebeam) >= 4.f * sigma_pair)
        return false;

    MollerEvent mev({xa, ya, za, Ea}, {xb, yb, zb, Eb});
    return std::fabs(physics.GetMollerPhiDiff(mev)) < 10.f;
}

static bool hasMollerPairInTriplet(const PhysicsTools &physics,
                                   float x1, float y1, float z1, float E1, float theta1,
                                   float x2, float y2, float z2, float E2, float theta2,
                                   float x3, float y3, float z3, float E3, float theta3,
                                   float ebeam)
{
    return isMollerPairCandidate(physics,
                                 x1, y1, z1, E1, theta1,
                                 x2, y2, z2, E2, theta2,
                                 ebeam)
           || isMollerPairCandidate(physics,
                                    x1, y1, z1, E1, theta1,
                                    x3, y3, z3, E3, theta3,
                                    ebeam)
           || isMollerPairCandidate(physics,
                                    x2, y2, z2, E2, theta2,
                                    x3, y3, z3, E3, theta3,
                                    ebeam);
}

const int Nbins = 33;
const float binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

struct QuickResult {
    std::unique_ptr<PhysicsTools> physics;
    Long64_t processed = 0;
    // For X17, only HyCal now
    std::unique_ptr<TH1F> h_3cl_totalE;
    std::unique_ptr<TH2F> h2_3cl_hits;
    std::unique_ptr<TH2F> h2_3cl_E_angle;
    std::unique_ptr<TH1F> h_3cl_E;
    std::unique_ptr<TH1F> h_3cl_yield;
    std::unique_ptr<TH1F> h_3cl_mass;
    std::unique_ptr<TH1F> h_3cl_cluster_num;
    std::unique_ptr<TH1F> h_3cl_ptx;
    std::unique_ptr<TH1F> h_3cl_pty;
    std::unique_ptr<TH2F> h2_3cl_Pt;
    std::unique_ptr<TH1F> h_3cl_tDiff;

    std::unique_ptr<TH1F> h_3cl_totalE_cut;
    std::unique_ptr<TH2F> h2_3cl_hits_cut;
    std::unique_ptr<TH2F> h2_3cl_E_angle_cut;
    std::unique_ptr<TH1F> h_3cl_E_cut;
    std::unique_ptr<TH1F> h_3cl_yield_cut;
    std::unique_ptr<TH1F> h_3cl_mass_cut;
    std::unique_ptr<TH1F> h_3cl_ptx_cut;
    std::unique_ptr<TH1F> h_3cl_pty_cut;
    std::unique_ptr<TH2F> h2_3cl_Pt_cut;
    std::unique_ptr<TH1F> h_3cl_tDiff_cut;

    std::unique_ptr<TH1F> h_3cl_cluster_gem_num;
    std::unique_ptr<TH1F> h_3cl_totalE_gem;
    std::unique_ptr<TH2F> h2_3cl_hits_gem;
    std::unique_ptr<TH2F> h2_3cl_E_angle_gem;
    std::unique_ptr<TH1F> h_3cl_E_gem;
    std::unique_ptr<TH1F> h_3cl_yield_gem;
    std::unique_ptr<TH1F> h_3cl_mass_gem;
    std::unique_ptr<TH1F> h_3cl_ptx_gem;
    std::unique_ptr<TH1F> h_3cl_pty_gem;
    std::unique_ptr<TH2F> h2_3cl_Pt_gem;
    std::unique_ptr<TH1F> h_3cl_tDiff_gem;

    std::unique_ptr<TH1F> h_3cl_totalE_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_hits_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_E_angle_gem_cut;
    std::unique_ptr<TH1F> h_3cl_E_gem_cut;
    std::unique_ptr<TH1F> h_3cl_yield_gem_cut;
    std::unique_ptr<TH1F> h_3cl_mass_gem_cut;
    std::unique_ptr<TH1F> h_3cl_ptx_gem_cut;
    std::unique_ptr<TH1F> h_3cl_pty_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_Pt_gem_cut;
    std::unique_ptr<TH1F> h_3cl_tDiff_gem_cut;
    std::unique_ptr<TH1F> h_3cl_Vz_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_Vxy_gem_cut;

    //cut step by step
    std::unique_ptr<TH2F> h_3cl_Pt_totalE_cut;
    std::unique_ptr<TH2F> h_3cl_Pt_totalE_time_cut;
    std::unique_ptr<TH2F> h_3cl_Pt_totalE_time_clusterE_cut;
    std::unique_ptr<TH2F> h_3cl_Pt_totalE_time_clusterE_moller_cut;
    std::unique_ptr<TH1F> h_3cl_mass_totalE_cut;
    std::unique_ptr<TH1F> h_3cl_mass_totalE_time_cut;
    std::unique_ptr<TH1F> h_3cl_mass_totalE_time_clusterE_cut;
    std::unique_ptr<TH1F> h_3cl_mass_totalE_time_clusterE_Pt_cut;
    std::unique_ptr<TH1F> h_3cl_mass_totalE_time_clusterE_Pt_moller_cut;
};

static void detach(TH1 *h)
{
    if (h) h->SetDirectory(nullptr);
}

static std::unique_ptr<QuickResult> makeResult(fdec::HyCalSystem &hycal)
{
    auto r = std::make_unique<QuickResult>();
    r->physics = std::make_unique<PhysicsTools>(hycal);

    //For X17, only HyCal now
    r->h2_3cl_hits = std::make_unique<TH2F>("3cl_hits",
        "3-Cluster Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_3cl_E_angle = std::make_unique<TH2F>("3cl_E_angle",
        "3-Cluster Energy vs Angle hycal;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_3cl_E = std::make_unique<TH1F>("3cl_E",
        "3-Cluster Energy;Energy (MeV);Counts", 3750, 0, 2500);
    r->h_3cl_totalE = std::make_unique<TH1F>("3cl_totalE",
        "3-Cluster Total Energy;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_yield = std::make_unique<TH1F>("3cl_yield",
        "3-Cluster Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_3cl_mass = std::make_unique<TH1F>("3cl_mass",
        "3-Cluster Inv. Mass;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_cluster_num = std::make_unique<TH1F>("3cl_cluster_num",
        "3-Cluster Number;Number of Clusters;Counts", 20, 0, 20);
    r->h_3cl_ptx = std::make_unique<TH1F>("3cl_ptx",
        "3-Cluster Ptx;Ptx (MeV);Counts", 200, -50, 50);
    r->h_3cl_pty = std::make_unique<TH1F>("3cl_pty",
        "3-Cluster Pty;Pty (MeV);Counts", 200, -50, 50);
    r->h2_3cl_Pt = std::make_unique<TH2F>("3cl_Pt",
        "3-Cluster Pt hycal;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_tDiff = std::make_unique<TH1F>("3cl_tDiff",
        "3-Cluster Time Difference;Time Difference (ns);Counts", 200, -50, 50);

    r->h2_3cl_hits_cut = std::make_unique<TH2F>("3cl_hits_cut",
        "3-Cluster Hit positions hycal - Cut;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_3cl_E_angle_cut = std::make_unique<TH2F>("3cl_E_angle_cut",
        "3-Cluster Energy vs Angle hycal - Cut;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_3cl_E_cut = std::make_unique<TH1F>("3cl_E_cut",
        "3-Cluster Energy - Cut;Energy (MeV);Counts", 3750, 0, 2500);
    r->h_3cl_totalE_cut = std::make_unique<TH1F>("3cl_totalE_cut",
        "3-Cluster Total Energy - Cut;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_yield_cut = std::make_unique<TH1F>("3cl_yield_cut",
        "3-Cluster Yield - Cut;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_3cl_mass_cut = std::make_unique<TH1F>("3cl_mass_cut",
        "3-Cluster Inv. Mass - Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_ptx_cut = std::make_unique<TH1F>("3cl_ptx_cut",
        "3-Cluster Ptx - Cut;Ptx (MeV);Counts", 200, -50, 50);
    r->h_3cl_pty_cut = std::make_unique<TH1F>("3cl_pty_cut",
        "3-Cluster Pty - Cut;Pty (MeV);Counts", 200, -50, 50);
    r->h2_3cl_Pt_cut = std::make_unique<TH2F>("3cl_Pt_cut",
        "3-Cluster Pt hycal - Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_tDiff_cut = std::make_unique<TH1F>("3cl_tDiff_cut",
        "3-Cluster Time Difference - Cut;Time Difference (ns);Counts", 200, -50, 50);

    // use gem matching to cut the 3-cluster events
    r->h_3cl_cluster_gem_num = std::make_unique<TH1F>( "3cl_cluster_gem_num",
        "GEM-matched Cluster Number;Number of Clusters;Counts", 20, 0, 20);
    r->h2_3cl_hits_gem = std::make_unique<TH2F>("3cl_hits_gem",
        "3-Cluster Hit positions hycal with GEM matching;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_3cl_E_angle_gem = std::make_unique<TH2F>("3cl_E_angle_gem",
        "3-Cluster Energy vs Angle hycal with GEM matching;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_3cl_E_gem = std::make_unique<TH1F>("3cl_E_gem",
        "3-Cluster Energy with GEM matching;Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_totalE_gem = std::make_unique<TH1F>("3cl_totalE_gem",
        "3-Cluster Total Energy with GEM matching;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_yield_gem = std::make_unique<TH1F>("3cl_yield_gem",
        "3-Cluster Yield with GEM matching;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_3cl_mass_gem = std::make_unique<TH1F>("3cl_mass_gem",
        "3-Cluster Inv. Mass with GEM matching;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_ptx_gem = std::make_unique<TH1F>("3cl_ptx_gem",
        "3-Cluster Ptx with GEM matching;Ptx (MeV);Counts", 200, -50, 50);
    r->h_3cl_pty_gem = std::make_unique<TH1F>("3cl_pty_gem",
        "3-Cluster Pty with GEM matching;Pty (MeV);Counts", 200, -50, 50);
    r->h_3cl_tDiff_gem = std::make_unique<TH1F>("3cl_tDiff_gem",
        "3-Cluster Time Difference with GEM matching;Time Difference (ns);Counts", 200, -50, 50);
    r->h2_3cl_Pt_gem = std::make_unique<TH2F>("3cl_Pt_gem",
        "3-Cluster Pt hycal with GEM matching;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);

    r->h2_3cl_hits_gem_cut = std::make_unique<TH2F>("3cl_hits_gem_cut",
        "3-Cluster Hit positions with GEM matching - Cut;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_3cl_E_angle_gem_cut = std::make_unique<TH2F>("3cl_E_angle_gem_cut",
        "3-Cluster Energy vs Angle with GEM matching - Cut;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_3cl_E_gem_cut = std::make_unique<TH1F>("3cl_E_gem_cut",
        "3-Cluster Energy with GEM matching - Cut;Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_totalE_gem_cut = std::make_unique<TH1F>("3cl_totalE_gem_cut",
        "3-Cluster Total Energy with GEM matching - Cut;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_yield_gem_cut = std::make_unique<TH1F>("3cl_yield_gem_cut",
        "3-Cluster Yield with GEM matching - Cut;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_3cl_mass_gem_cut = std::make_unique<TH1F>("3cl_mass_gem_cut",
        "3-Cluster Inv. Mass with GEM matching - Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_ptx_gem_cut = std::make_unique<TH1F>("3cl_ptx_gem_cut",
        "3-Cluster Ptx with GEM matching - Cut;Ptx (MeV);Counts", 200, -50, 50);
    r->h_3cl_pty_gem_cut = std::make_unique<TH1F>("3cl_pty_gem_cut",
        "3-Cluster Pty with GEM matching - Cut;Pty (MeV);Counts", 200, -50, 50);
    r->h2_3cl_Pt_gem_cut = std::make_unique<TH2F>("3cl_Pt_gem_cut",
        "3-Cluster Pt hycal with GEM matching - Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_tDiff_gem_cut = std::make_unique<TH1F>("3cl_tDiff_gem_cut",
        "3-Cluster Time Difference with GEM matching - Cut;Time Difference (ns);Counts", 200, -50, 50);
    r->h_3cl_Vz_gem_cut = std::make_unique<TH1F>("3cl_Vz_gem_cut",
        "3-Cluster GEM Vertex Z with Cut;Z (mm);Counts", 14000, -5000, 9000);
    r->h2_3cl_Vxy_gem_cut = std::make_unique<TH2F>("3cl_Vxy_gem_cut",
        "3-Cluster GEM Vertex XY with Cut;X (mm);Y (mm)", 800, -40, 40, 800, -40, 40);
    // step by step cuts
    r->h_3cl_Pt_totalE_cut = std::make_unique<TH2F>("3cl_Pt_totalE_cut",
        "3-Cluster Pt with Total Energy Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_Pt_totalE_time_cut = std::make_unique<TH2F>("3cl_Pt_totalE_time_cut",
        "3-Cluster Pt with Total Energy and Time Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_Pt_totalE_time_clusterE_cut = std::make_unique<TH2F>("3cl_Pt_totalE_time_clusterE_cut",
        "3-Cluster Pt with Total Energy, Time, and Cluster Energy Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_Pt_totalE_time_clusterE_moller_cut = std::make_unique<TH2F>("3cl_Pt_totalE_time_clusterE_moller_cut",
        "3-Cluster Pt with Total Energy, Time, Cluster Energy, and Moller Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_mass_totalE_cut = std::make_unique<TH1F>("3cl_mass_totalE_cut",
        "3-Cluster Inv. Mass with Total Energy Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_mass_totalE_time_cut = std::make_unique<TH1F>("3cl_mass_totalE_time_cut",
        "3-Cluster Inv. Mass with Total Energy and Time Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_mass_totalE_time_clusterE_cut = std::make_unique<TH1F>("3cl_mass_totalE_time_clusterE_cut",
        "3-Cluster Inv. Mass with Total Energy, Time, and Cluster Energy Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_mass_totalE_time_clusterE_Pt_cut = std::make_unique<TH1F>("3cl_mass_totalE_time_clusterE_Pt_cut",
        "3-Cluster Inv. Mass with Total Energy, Time, Cluster Energy, and Pt Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_mass_totalE_time_clusterE_Pt_moller_cut = std::make_unique<TH1F>("3cl_mass_totalE_time_clusterE_Pt_moller_cut",
        "3-Cluster Inv. Mass with Total Energy, Time, Cluster Energy, Pt, and Moller Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);

    detach(r->h2_3cl_hits.get());
    detach(r->h2_3cl_E_angle.get());
    detach(r->h_3cl_E.get());
    detach(r->h_3cl_totalE.get());
    detach(r->h_3cl_yield.get());
    detach(r->h_3cl_mass.get());
    detach(r->h_3cl_cluster_num.get());
    detach(r->h_3cl_ptx.get());
    detach(r->h_3cl_pty.get());
    detach(r->h2_3cl_Pt.get());
    detach(r->h_3cl_tDiff.get());

    detach(r->h2_3cl_hits_cut.get());
    detach(r->h2_3cl_E_angle_cut.get());
    detach(r->h_3cl_E_cut.get());
    detach(r->h_3cl_totalE_cut.get());
    detach(r->h_3cl_yield_cut.get());
    detach(r->h_3cl_mass_cut.get());
    detach(r->h_3cl_ptx_cut.get());
    detach(r->h_3cl_pty_cut.get());
    detach(r->h2_3cl_Pt_cut.get());
    detach(r->h_3cl_tDiff_cut.get());
    
    detach(r->h_3cl_cluster_gem_num.get());
    detach(r->h2_3cl_hits_gem.get());
    detach(r->h2_3cl_E_angle_gem.get());
    detach(r->h_3cl_E_gem.get());
    detach(r->h_3cl_totalE_gem.get());
    detach(r->h_3cl_yield_gem.get());
    detach(r->h_3cl_mass_gem.get());
    detach(r->h_3cl_ptx_gem.get());
    detach(r->h_3cl_pty_gem.get());
    detach(r->h2_3cl_Pt_gem.get());
    detach(r->h_3cl_tDiff_gem.get());
    detach(r->h2_3cl_hits_gem_cut.get());
    detach(r->h2_3cl_E_angle_gem_cut.get());
    detach(r->h_3cl_E_gem_cut.get());
    detach(r->h_3cl_totalE_gem_cut.get());
    detach(r->h_3cl_yield_gem_cut.get());
    detach(r->h_3cl_mass_gem_cut.get());
    detach(r->h_3cl_ptx_gem_cut.get());
    detach(r->h_3cl_pty_gem_cut.get());
    detach(r->h2_3cl_Pt_gem_cut.get());
    detach(r->h_3cl_tDiff_gem_cut.get());
    detach(r->h_3cl_Vz_gem_cut.get());
    detach(r->h2_3cl_Vxy_gem_cut.get());
    detach(r->h_3cl_Pt_totalE_cut.get());
    detach(r->h_3cl_Pt_totalE_time_cut.get());
    detach(r->h_3cl_Pt_totalE_time_clusterE_cut.get());
    detach(r->h_3cl_Pt_totalE_time_clusterE_moller_cut.get());
    detach(r->h_3cl_mass_totalE_cut.get());
    detach(r->h_3cl_mass_totalE_time_cut.get());
    detach(r->h_3cl_mass_totalE_time_clusterE_cut.get());
    detach(r->h_3cl_mass_totalE_time_clusterE_Pt_cut.get());
    detach(r->h_3cl_mass_totalE_time_clusterE_Pt_moller_cut.get());
    return r;
}

static Long64_t reconEntries(const std::string &path)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return 0;
    TTree *t = dynamic_cast<TTree *>(f->Get("recon"));
    return t ? t->GetEntries() : 0;
}

static bool processFile(const std::string &path,
                        Long64_t max_entries,
                        fdec::HyCalSystem &hycal,
                        float Ebeam,
                        QuickResult &out)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << path << "\n";
        return false;
    }
    TTree *tree = dynamic_cast<TTree *>(f->Get("recon"));
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in " << path << "\n";
        return false;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);
    Long64_t n = tree->GetEntries();
    if (max_entries >= 0 && max_entries < n) n = max_entries;

    auto &physics = *out.physics;
    for (Long64_t i = 0; i < n; i++) {
        tree->GetEntry(i);

        // trigger selection
        bool is_3cluster = (ev.trigger_bits & prad2::TBIT_3cl) != 0;
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;

        if (!is_3cluster && !is_sum) continue;

        // x17 trigger selection
        if (is_3cluster) {
            out.h_3cl_cluster_num->Fill(ev.n_clusters);
            if(ev.n_clusters == 3) {
                float x1 = ev.cl_x[0], y1 = ev.cl_y[0], z1 = ev.cl_z[0];
                float x2 = ev.cl_x[1], y2 = ev.cl_y[1], z2 = ev.cl_z[1];
                float x3 = ev.cl_x[2], y3 = ev.cl_y[2], z3 = ev.cl_z[2];
                float E1 = ev.cl_energy[0], E2 = ev.cl_energy[1], E3 = ev.cl_energy[2];

                float t1 = ev.cl_time[0], t2 = ev.cl_time[1], t3 = ev.cl_time[2];
                float tDiff = std::max({std::fabs(t1 - t2), std::fabs(t1 - t3), std::fabs(t2 - t3)});

                bool nblocks_ok = true;
                if(ev.cl_nblocks[0] < 1 || ev.cl_nblocks[1] < 1 || ev.cl_nblocks[2] < 1)
                    nblocks_ok = false;

                float theta1 = std::atan2(std::sqrt(x1*x1 + y1*y1), z1) * 180.f / M_PI;
                float theta2 = std::atan2(std::sqrt(x2*x2 + y2*y2), z2) * 180.f / M_PI;
                float theta3 = std::atan2(std::sqrt(x3*x3 + y3*y3), z3) * 180.f / M_PI;

                float mass1 = electronPairInvariantMass(x1, y1, z1, E1, x2, y2, z2, E2);
                float mass2 = electronPairInvariantMass(x1, y1, z1, E1, x3, y3, z3, E3);
                float mass3 = electronPairInvariantMass(x2, y2, z2, E2, x3, y3, z3, E3);

                // Pt x and Pt y calculation
                auto get_pt = [](float x, float y, float z, float energy) {
                    constexpr float electron_mass = 0.51099895f; // MeV
                    const float norm = std::sqrt(x*x + y*y + z*z);
                    if (norm <= 0.f || energy < electron_mass)
                        return std::pair<float, float>{0.f, 0.f};
                    const float p = std::sqrt(std::max(
                        0.f, energy*energy - electron_mass*electron_mass));
                    return std::pair<float, float>{p*x/norm, p*y/norm};
                };
                const auto [px1, py1] = get_pt(x1, y1, z1, E1);
                const auto [px2, py2] = get_pt(x2, y2, z2, E2);
                const auto [px3, py3] = get_pt(x3, y3, z3, E3);
                const float ptx = px1 + px2 + px3;
                const float pty = py1 + py2 + py3;

                bool time_cut = false, totalE_cut = false, Pt_cut = false, clusterE_cut = false, inHyCal_cut = false;
                bool moller_cut = false;

                if(tDiff < 16.f)
                    time_cut = true;
                if(E1 + E2 + E3 < Ebeam + 250.f && E1 + E2 + E3 > 0.8 * Ebeam)
                    totalE_cut = true;
                if(std::sqrt(ptx*ptx + pty*pty) < 5.f)
                    Pt_cut = true;
                if(inHyCal(x1, y1) && inHyCal(x2, y2) && inHyCal(x3, y3))
                    inHyCal_cut = true;
                if(E1 > 70.f && E2 > 70.f && E3 > 70.f && E1 < 0.75 * Ebeam && E2 < 0.75 * Ebeam && E3 < 0.75 * Ebeam)
                    clusterE_cut = true;

                const bool has_moller_pair = hasMollerPairInTriplet(
                    physics,
                    x1, y1, z1, E1, theta1,
                    x2, y2, z2, E2, theta2,
                    x3, y3, z3, E3, theta3,
                    Ebeam);
                //moller_cut = !has_moller_pair;
                moller_cut = !(isMollerAngleEnergyMatch(physics, theta1, E1, Ebeam) ||
                             isMollerAngleEnergyMatch(physics, theta2, E2, Ebeam) ||
                             isMollerAngleEnergyMatch(physics, theta3, E3, Ebeam));

                out.h_3cl_totalE->Fill(E1 + E2 + E3);
                out.h_3cl_ptx->Fill(ptx);
                out.h_3cl_pty->Fill(pty);
                out.h2_3cl_Pt->Fill(ptx, pty);
                out.h_3cl_tDiff->Fill(tDiff);
                out.h_3cl_E->Fill(E1);
                out.h_3cl_E->Fill(E2);
                out.h_3cl_E->Fill(E3);
                out.h2_3cl_hits->Fill(x1, y1);
                out.h2_3cl_hits->Fill(x2, y2);
                out.h2_3cl_hits->Fill(x3, y3);
                out.h2_3cl_E_angle->Fill(theta1, E1);
                out.h2_3cl_E_angle->Fill(theta2, E2);
                out.h2_3cl_E_angle->Fill(theta3, E3);
                out.h_3cl_yield->Fill(theta1);
                out.h_3cl_yield->Fill(theta2);
                out.h_3cl_yield->Fill(theta3);
                if (std::isfinite(mass1)) out.h_3cl_mass->Fill(mass1);
                if (std::isfinite(mass2)) out.h_3cl_mass->Fill(mass2);
                if (std::isfinite(mass3)) out.h_3cl_mass->Fill(mass3);

                if(nblocks_ok && totalE_cut && time_cut && Pt_cut && inHyCal_cut && clusterE_cut && moller_cut) {
                    out.h_3cl_tDiff_cut->Fill(tDiff);
                    out.h_3cl_totalE_cut->Fill(E1 + E2 + E3);
                    out.h_3cl_ptx_cut->Fill(ptx);
                    out.h_3cl_pty_cut->Fill(pty);
                    out.h2_3cl_Pt_cut->Fill(ptx, pty);
                    out.h_3cl_E_cut->Fill(E1);
                    out.h_3cl_E_cut->Fill(E2);
                    out.h_3cl_E_cut->Fill(E3);
                    out.h2_3cl_hits_cut->Fill(x1, y1);
                    out.h2_3cl_hits_cut->Fill(x2, y2);
                    out.h2_3cl_hits_cut->Fill(x3, y3);
                    out.h2_3cl_E_angle_cut->Fill(theta1, E1);
                    out.h2_3cl_E_angle_cut->Fill(theta2, E2);
                    out.h2_3cl_E_angle_cut->Fill(theta3, E3);
                    out.h_3cl_yield_cut->Fill(theta1);
                    out.h_3cl_yield_cut->Fill(theta2);
                    out.h_3cl_yield_cut->Fill(theta3);
                    if (std::isfinite(mass1)) out.h_3cl_mass_cut->Fill(mass1);
                    if (std::isfinite(mass2)) out.h_3cl_mass_cut->Fill(mass2);
                    if (std::isfinite(mass3)) out.h_3cl_mass_cut->Fill(mass3);
                }

            }

            //loop over all clusters for GEM matching
            struct Hits{
                float xu, yu, zu; // upstream
                float xd, yd, zd; // downstream
                float x, y, z; // projected to HyCal plane
                float E, t; // cluster energy and time
                int index; // cluster index
                bool gem_matched[4] = {false, false, false, false}; // 0,1 downstream, 2,3 upstream
            };

            std::vector<Hits> hits_candidate;

            for (int j = 0; j < ev.n_clusters; j++) {
                if(ev.cl_nblocks[j] < 1) continue;
                if(fdec::test_bit(ev.cl_flag[j], fdec::kInnerBound)) continue;
                if(fdec::test_bit(ev.cl_flag[j], fdec::kOuterBound)) continue;
                if(ev.cl_energy[j] < 70.f || ev.cl_energy[j] > 0.75 * Ebeam) continue;

                bool gem[4] = {false, false, false, false}; // 0,1 downstream, 2,3 upstream
                if ((ev.matchFlag[j] & 1u << 0) != 0) gem[0] = true;
                if ((ev.matchFlag[j] & 1u << 1) != 0) gem[1] = true;
                if ((ev.matchFlag[j] & 1u << 2) != 0) gem[2] = true;
                if ((ev.matchFlag[j] & 1u << 3) != 0) gem[3] = true;

                const bool has_downstream = gem[0] || gem[1];
                const bool has_upstream = gem[2] || gem[3];
                if (!has_downstream || !has_upstream) continue;

                const int u = gem[2] ? 2 : 3;
                const int d = gem[0] ? 0 : 1;

                Hits hit{};
                hit.xu = ev.matchGEMx[j][u];
                hit.yu = ev.matchGEMy[j][u];
                hit.zu = ev.matchGEMz[j][u];
                hit.xd = ev.matchGEMx[j][d];
                hit.yd = ev.matchGEMy[j][d];
                hit.zd = ev.matchGEMz[j][d];
                float scale = ev.cl_z[j] / hit.zu;
                hit.x = hit.xu * scale;
                hit.y = hit.yu * scale;
                hit.z = hit.zu * scale;
                hit.E = ev.cl_energy[j];
                hit.t = ev.cl_time[j];
                hit.index = j;
                hit.gem_matched[0] = gem[0];
                hit.gem_matched[1] = gem[1];
                hit.gem_matched[2] = gem[2];
                hit.gem_matched[3] = gem[3];
                hits_candidate.push_back(hit);
            }

            //HyCal cluster multi-pulse mode may have 2 events in one event
            //need to check the time to split the clusters into two groups by time
            //if the time difference is larger than 16 ns, we consider them as two different events
            //start from the first cluster(already sorted by energy), loop over the rest clusters, if the time
            //difference is smaller than 16 ns, we consider them as the same event, otherwise we consider them as different events

            int gem_matched_count = 0;
            float x1 = 0.f, y1 = 0.f, z1 = 0.f, x2 = 0.f, y2 = 0.f, z2 = 0.f, x3 = 0.f, y3 = 0.f, z3 = 0.f;
            float E1 = 0.f, E2 = 0.f, E3 = 0.f, t1 = 0.f, t2 = 0.f, t3 = 0.f;
            float x1u = 0.f, y1u = 0.f, z1u = 0.f, x2u = 0.f, y2u = 0.f, z2u = 0.f, x3u = 0.f, y3u = 0.f, z3u = 0.f;
            float x1d = 0.f, y1d = 0.f, z1d = 0.f, x2d = 0.f, y2d = 0.f, z2d = 0.f, x3d = 0.f, y3d = 0.f, z3d = 0.f;
            int index1 = -1, index2 = -1, index3 = -1;

            std::sort(hits_candidate.begin(), hits_candidate.end(),
                      [](const Hits &a, const Hits &b) { return a.E > b.E; });

            std::vector<Hits> time_group_1;
            std::vector<Hits> time_group_2;
            if (!hits_candidate.empty()) {
                const float t0 = hits_candidate.front().t;
                std::vector<Hits> remaining_hits;
                for (const auto &h : hits_candidate) {
                    if (std::fabs(h.t - t0) < 16.f)
                        time_group_1.push_back(h);
                    else
                        remaining_hits.push_back(h);
                }

                if (!remaining_hits.empty()) {
                    const float t1 = remaining_hits.front().t;
                    for (const auto &h : remaining_hits) {
                        if (std::fabs(h.t - t1) < 16.f)
                            time_group_2.push_back(h);
                    }
                }
            }
            out.h_3cl_cluster_gem_num->Fill(time_group_1.size());
            if (!time_group_2.empty()) out.h_3cl_cluster_gem_num->Fill(time_group_2.size());

            if (time_group_1.size() == 3) {
                const Hits &h1 = time_group_1[0];
                const Hits &h2 = time_group_1[1];
                const Hits &h3 = time_group_1[2];

                x1 = h1.x; y1 = h1.y; z1 = h1.z; E1 = h1.E; t1 = h1.t; index1 = h1.index;
                x2 = h2.x; y2 = h2.y; z2 = h2.z; E2 = h2.E; t2 = h2.t; index2 = h2.index;
                x3 = h3.x; y3 = h3.y; z3 = h3.z; E3 = h3.E; t3 = h3.t; index3 = h3.index;

                x1u = h1.xu; y1u = h1.yu; z1u = h1.zu; x1d = h1.xd; y1d = h1.yd; z1d = h1.zd;
                x2u = h2.xu; y2u = h2.yu; z2u = h2.zu; x2d = h2.xd; y2d = h2.yd; z2d = h2.zd;
                x3u = h3.xu; y3u = h3.yu; z3u = h3.zu; x3d = h3.xd; y3d = h3.yd; z3d = h3.zd;

                float tDiff = std::max({std::fabs(t1 - t2), std::fabs(t1 - t3), std::fabs(t2 - t3)});

                float theta1 = std::atan2(std::sqrt(x1*x1 + y1*y1), z1) * 180.f / M_PI;
                float theta2 = std::atan2(std::sqrt(x2*x2 + y2*y2), z2) * 180.f / M_PI;
                float theta3 = std::atan2(std::sqrt(x3*x3 + y3*y3), z3) * 180.f / M_PI;

                const GEMTrackPair gem_tracks[3] = {
                    {x1u, y1u, z1u, x1d, y1d, z1d},
                    {x2u, y2u, z2u, x2d, y2d, z2d},
                    {x3u, y3u, z3u, x3d, y3d, z3d}
                };
                std::array<float, 3> gem_vertex = {0.f, 0.f, 0.f};
                const bool gem_vertex_ok = solveGEMVertex(gem_tracks, gem_vertex);
                const float vx = gem_vertex_ok ? gem_vertex[0] : -999.f;
                const float vy = gem_vertex_ok ? gem_vertex[1] : -999.f;
                const float vz = gem_vertex_ok ? gem_vertex[2] : -9999.f;

                auto rel_coord = [&](float x, float y, float z) {
                    return std::array<float, 3>{x - vx, y - vy, z - vz};
                };
                auto theta_from_vertex = [&](float x, float y, float z) {
                    const auto rel = rel_coord(x, y, z);
                    return std::atan2(std::sqrt(rel[0] * rel[0] + rel[1] * rel[1]), rel[2]) * 180.f / M_PI;
                };

                float mass1 = electronPairInvariantMass(x1, y1, z1, E1, x2, y2, z2, E2);
                float mass2 = electronPairInvariantMass(x1, y1, z1, E1, x3, y3, z3, E3);
                float mass3 = electronPairInvariantMass(x2, y2, z2, E2, x3, y3, z3, E3);

                // Pt x and Pt y calculation
                auto get_pt = [](float x, float y, float z, float energy) {
                    constexpr float electron_mass = 0.51099895f; // MeV
                    const float norm = std::sqrt(x*x + y*y + z*z);
                    if (norm <= 0.f || energy < electron_mass)
                        return std::pair<float, float>{0.f, 0.f};
                    const float p = std::sqrt(std::max(
                        0.f, energy*energy - electron_mass*electron_mass));
                    return std::pair<float, float>{p*x/norm, p*y/norm};
                };
                const auto [px1, py1] = get_pt(x1, y1, z1, E1);
                const auto [px2, py2] = get_pt(x2, y2, z2, E2);
                const auto [px3, py3] = get_pt(x3, y3, z3, E3);
                const float ptx = px1 + px2 + px3;
                const float pty = py1 + py2 + py3;

                bool time_cut = false, totalE_cut = false, Pt_cut = false, clusterE_cut = false, inHyCal_cut = false;

                bool moller_cut = false;
                const bool has_moller_pair = hasMollerPairInTriplet(
                    physics,
                    x1, y1, z1, E1, theta1,
                    x2, y2, z2, E2, theta2,
                    x3, y3, z3, E3, theta3,
                    Ebeam);
                //moller_cut = !has_moller_pair;
                moller_cut = !(isMollerAngleEnergyMatch(physics, theta1, E1, Ebeam) ||
                             isMollerAngleEnergyMatch(physics, theta2, E2, Ebeam) ||
                             isMollerAngleEnergyMatch(physics, theta3, E3, Ebeam));

                if(tDiff < 16.f)
                    time_cut = true;
                if(E1 + E2 + E3 < Ebeam + 250.f && E1 + E2 + E3 > 0.8 * Ebeam)
                    totalE_cut = true;
                if(std::sqrt(ptx*ptx + pty*pty) < 5.f)
                    Pt_cut = true;
                if(inHyCal(x1, y1) && inHyCal(x2, y2) && inHyCal(x3, y3))
                    inHyCal_cut = true;
                if(E1 > 70.f && E2 > 70.f && E3 > 70.f && E1 < 0.75 * Ebeam && E2 < 0.75 * Ebeam && E3 < 0.75 * Ebeam)
                    clusterE_cut = true;

                out.h_3cl_totalE_gem->Fill(E1 + E2 + E3);
                out.h_3cl_tDiff_gem->Fill(tDiff);
                out.h_3cl_ptx_gem->Fill(ptx);
                out.h_3cl_pty_gem->Fill(pty);
                out.h2_3cl_Pt_gem->Fill(ptx, pty);
                out.h_3cl_E_gem->Fill(E1);
                out.h_3cl_E_gem->Fill(E2);
                out.h_3cl_E_gem->Fill(E3);
                out.h2_3cl_hits_gem->Fill(x1, y1);
                out.h2_3cl_hits_gem->Fill(x2, y2);
                out.h2_3cl_hits_gem->Fill(x3, y3);
                out.h2_3cl_E_angle_gem->Fill(theta1, E1);
                out.h2_3cl_E_angle_gem->Fill(theta2, E2);
                out.h2_3cl_E_angle_gem->Fill(theta3, E3);
                out.h_3cl_yield_gem->Fill(theta1);
                out.h_3cl_yield_gem->Fill(theta2);
                out.h_3cl_yield_gem->Fill(theta3);
                //if (std::isfinite(mass1)) out.h_3cl_mass_gem->Fill(mass1);
                if (std::isfinite(mass2)) out.h_3cl_mass_gem->Fill(mass2);
                if (std::isfinite(mass3)) out.h_3cl_mass_gem->Fill(mass3);

                auto fill_mass_step = [&](TH1F *hist) {
                    //if (std::isfinite(mass1)) hist->Fill(mass1);
                    if (std::isfinite(mass2)) hist->Fill(mass2);
                    if (std::isfinite(mass3)) hist->Fill(mass3);
                };

                if (totalE_cut) {
                    out.h_3cl_Pt_totalE_cut->Fill(ptx, pty);
                    fill_mass_step(out.h_3cl_mass_totalE_cut.get());
                }
                if (totalE_cut && time_cut) {
                    out.h_3cl_Pt_totalE_time_cut->Fill(ptx, pty);
                    fill_mass_step(out.h_3cl_mass_totalE_time_cut.get());
                }
                if (totalE_cut && time_cut && clusterE_cut) {
                    out.h_3cl_Pt_totalE_time_clusterE_cut->Fill(ptx, pty);
                    fill_mass_step(out.h_3cl_mass_totalE_time_clusterE_cut.get());
                }
                if (totalE_cut && time_cut && clusterE_cut && Pt_cut) {
                    fill_mass_step(out.h_3cl_mass_totalE_time_clusterE_Pt_cut.get());
                }
                if (totalE_cut && time_cut && clusterE_cut && moller_cut) {
                    out.h_3cl_Pt_totalE_time_clusterE_moller_cut->Fill(ptx, pty);
                }
                if (totalE_cut && time_cut && clusterE_cut && Pt_cut && moller_cut) {
                    fill_mass_step(out.h_3cl_mass_totalE_time_clusterE_Pt_moller_cut.get());
                }

                if(time_cut && totalE_cut && Pt_cut && clusterE_cut && inHyCal_cut && moller_cut) {
                    out.h_3cl_totalE_gem_cut->Fill(E1 + E2 + E3);
                    out.h_3cl_ptx_gem_cut->Fill(ptx);
                    out.h_3cl_pty_gem_cut->Fill(pty);
                    out.h2_3cl_Pt_gem_cut->Fill(ptx, pty);
                    out.h_3cl_Vz_gem_cut->Fill(vz);
                    out.h2_3cl_Vxy_gem_cut->Fill(vx, vy);
                    out.h_3cl_E_gem_cut->Fill(E1);
                    out.h_3cl_E_gem_cut->Fill(E2);
                    out.h_3cl_E_gem_cut->Fill(E3);
                    out.h_3cl_tDiff_gem_cut->Fill(tDiff);
                    out.h2_3cl_hits_gem_cut->Fill(x1, y1);
                    out.h2_3cl_hits_gem_cut->Fill(x2, y2);
                    out.h2_3cl_hits_gem_cut->Fill(x3, y3);
                    out.h2_3cl_E_angle_gem_cut->Fill(theta1, E1);
                    out.h2_3cl_E_angle_gem_cut->Fill(theta2, E2);
                    out.h2_3cl_E_angle_gem_cut->Fill(theta3, E3);
                    out.h_3cl_yield_gem_cut->Fill(theta1);
                    out.h_3cl_yield_gem_cut->Fill(theta2);
                    out.h_3cl_yield_gem_cut->Fill(theta3);
                    //if (std::isfinite(mass1)) out.h_3cl_mass_gem_cut->Fill(mass1);
                    if (std::isfinite(mass2)) out.h_3cl_mass_gem_cut->Fill(mass2);
                    if (std::isfinite(mass3)) out.h_3cl_mass_gem_cut->Fill(mass3);
                }
            }
            if (time_group_2.size() == 3) {
                const Hits &h1 = time_group_2[0];
                const Hits &h2 = time_group_2[1];
                const Hits &h3 = time_group_2[2];

                x1 = h1.x; y1 = h1.y; z1 = h1.z; E1 = h1.E; t1 = h1.t; index1 = h1.index;
                x2 = h2.x; y2 = h2.y; z2 = h2.z; E2 = h2.E; t2 = h2.t; index2 = h2.index;
                x3 = h3.x; y3 = h3.y; z3 = h3.z; E3 = h3.E; t3 = h3.t; index3 = h3.index;

                x1u = h1.xu; y1u = h1.yu; z1u = h1.zu; x1d = h1.xd; y1d = h1.yd; z1d = h1.zd;
                x2u = h2.xu; y2u = h2.yu; z2u = h2.zu; x2d = h2.xd; y2d = h2.yd; z2d = h2.zd;
                x3u = h3.xu; y3u = h3.yu; z3u = h3.zu; x3d = h3.xd; y3d = h3.yd; z3d = h3.zd;

                float tDiff = std::max({std::fabs(t1 - t2), std::fabs(t1 - t3), std::fabs(t2 - t3)});

                float theta1 = std::atan2(std::sqrt(x1*x1 + y1*y1), z1) * 180.f / M_PI;
                float theta2 = std::atan2(std::sqrt(x2*x2 + y2*y2), z2) * 180.f / M_PI;
                float theta3 = std::atan2(std::sqrt(x3*x3 + y3*y3), z3) * 180.f / M_PI;

                const GEMTrackPair gem_tracks[3] = {
                    {x1u, y1u, z1u, x1d, y1d, z1d},
                    {x2u, y2u, z2u, x2d, y2d, z2d},
                    {x3u, y3u, z3u, x3d, y3d, z3d}
                };
                std::array<float, 3> gem_vertex = {0.f, 0.f, 0.f};
                const bool gem_vertex_ok = solveGEMVertex(gem_tracks, gem_vertex);
                const float vx = gem_vertex_ok ? gem_vertex[0] : -999.f;
                const float vy = gem_vertex_ok ? gem_vertex[1] : -999.f;
                const float vz = gem_vertex_ok ? gem_vertex[2] : -9999.f;

                auto rel_coord = [&](float x, float y, float z) {
                    return std::array<float, 3>{x - vx, y - vy, z - vz};
                };
                auto theta_from_vertex = [&](float x, float y, float z) {
                    const auto rel = rel_coord(x, y, z);
                    return std::atan2(std::sqrt(rel[0] * rel[0] + rel[1] * rel[1]), rel[2]) * 180.f / M_PI;
                };

                float mass1 = electronPairInvariantMass(x1, y1, z1, E1, x2, y2, z2, E2);
                float mass2 = electronPairInvariantMass(x1, y1, z1, E1, x3, y3, z3, E3);
                float mass3 = electronPairInvariantMass(x2, y2, z2, E2, x3, y3, z3, E3);

                // Pt x and Pt y calculation
                auto get_pt = [](float x, float y, float z, float energy) {
                    constexpr float electron_mass = 0.51099895f; // MeV
                    const float norm = std::sqrt(x*x + y*y + z*z);
                    if (norm <= 0.f || energy < electron_mass)
                        return std::pair<float, float>{0.f, 0.f};
                    const float p = std::sqrt(std::max(
                        0.f, energy*energy - electron_mass*electron_mass));
                    return std::pair<float, float>{p*x/norm, p*y/norm};
                };
                const auto [px1, py1] = get_pt(x1, y1, z1, E1);
                const auto [px2, py2] = get_pt(x2, y2, z2, E2);
                const auto [px3, py3] = get_pt(x3, y3, z3, E3);
                const float ptx = px1 + px2 + px3;
                const float pty = py1 + py2 + py3;

                bool time_cut = false, totalE_cut = false, Pt_cut = false, clusterE_cut = false, inHyCal_cut = false;

                bool moller_cut = false;
                const bool has_moller_pair = hasMollerPairInTriplet(
                    physics,
                    x1, y1, z1, E1, theta1,
                    x2, y2, z2, E2, theta2,
                    x3, y3, z3, E3, theta3,
                    Ebeam);
                //moller_cut = !has_moller_pair;
                moller_cut = !(isMollerAngleEnergyMatch(physics, theta1, E1, Ebeam) ||
                             isMollerAngleEnergyMatch(physics, theta2, E2, Ebeam) ||
                             isMollerAngleEnergyMatch(physics, theta3, E3, Ebeam));

                if(tDiff < 16.f)
                    time_cut = true;
                if(E1 + E2 + E3 < Ebeam + 250.f && E1 + E2 + E3 > 0.8 * Ebeam)
                    totalE_cut = true;
                if(std::sqrt(ptx*ptx + pty*pty) < 5.f)
                    Pt_cut = true;
                if(inHyCal(x1, y1) && inHyCal(x2, y2) && inHyCal(x3, y3))
                    inHyCal_cut = true;
                if(E1 > 70.f && E2 > 70.f && E3 > 70.f && E1 < 0.75 * Ebeam && E2 < 0.75 * Ebeam && E3 < 0.75 * Ebeam)
                    clusterE_cut = true;

                out.h_3cl_totalE_gem->Fill(E1 + E2 + E3);
                out.h_3cl_tDiff_gem->Fill(tDiff);
                out.h_3cl_ptx_gem->Fill(ptx);
                out.h_3cl_pty_gem->Fill(pty);
                out.h2_3cl_Pt_gem->Fill(ptx, pty);
                out.h_3cl_E_gem->Fill(E1);
                out.h_3cl_E_gem->Fill(E2);
                out.h_3cl_E_gem->Fill(E3);
                out.h2_3cl_hits_gem->Fill(x1, y1);
                out.h2_3cl_hits_gem->Fill(x2, y2);
                out.h2_3cl_hits_gem->Fill(x3, y3);
                out.h2_3cl_E_angle_gem->Fill(theta1, E1);
                out.h2_3cl_E_angle_gem->Fill(theta2, E2);
                out.h2_3cl_E_angle_gem->Fill(theta3, E3);
                out.h_3cl_yield_gem->Fill(theta1);
                out.h_3cl_yield_gem->Fill(theta2);
                out.h_3cl_yield_gem->Fill(theta3);
                //if (std::isfinite(mass1)) out.h_3cl_mass_gem->Fill(mass1);
                if (std::isfinite(mass2)) out.h_3cl_mass_gem->Fill(mass2);
                if (std::isfinite(mass3)) out.h_3cl_mass_gem->Fill(mass3);

                auto fill_mass_step = [&](TH1F *hist) {
                    //if (std::isfinite(mass1)) hist->Fill(mass1);
                    if (std::isfinite(mass2)) hist->Fill(mass2);
                    if (std::isfinite(mass3)) hist->Fill(mass3);
                };

                if (totalE_cut) {
                    out.h_3cl_Pt_totalE_cut->Fill(ptx, pty);
                    fill_mass_step(out.h_3cl_mass_totalE_cut.get());
                }
                if (totalE_cut && time_cut) {
                    out.h_3cl_Pt_totalE_time_cut->Fill(ptx, pty);
                    fill_mass_step(out.h_3cl_mass_totalE_time_cut.get());
                }
                if (totalE_cut && time_cut && clusterE_cut) {
                    out.h_3cl_Pt_totalE_time_clusterE_cut->Fill(ptx, pty);
                    fill_mass_step(out.h_3cl_mass_totalE_time_clusterE_cut.get());
                }
                if (totalE_cut && time_cut && clusterE_cut && Pt_cut) {
                    fill_mass_step(out.h_3cl_mass_totalE_time_clusterE_Pt_cut.get());
                }
                if (totalE_cut && time_cut && clusterE_cut && moller_cut) {
                    out.h_3cl_Pt_totalE_time_clusterE_moller_cut->Fill(ptx, pty);
                }
                if (totalE_cut && time_cut && clusterE_cut && Pt_cut && moller_cut) {
                    fill_mass_step(out.h_3cl_mass_totalE_time_clusterE_Pt_moller_cut.get());
                }

                if(time_cut && totalE_cut && Pt_cut && clusterE_cut && inHyCal_cut && moller_cut) {
                    out.h_3cl_totalE_gem_cut->Fill(E1 + E2 + E3);
                    out.h_3cl_ptx_gem_cut->Fill(ptx);
                    out.h_3cl_pty_gem_cut->Fill(pty);
                    out.h2_3cl_Pt_gem_cut->Fill(ptx, pty);
                    out.h_3cl_Vz_gem_cut->Fill(vz);
                    out.h2_3cl_Vxy_gem_cut->Fill(vx, vy);
                    out.h_3cl_E_gem_cut->Fill(E1);
                    out.h_3cl_E_gem_cut->Fill(E2);
                    out.h_3cl_E_gem_cut->Fill(E3);
                    out.h_3cl_tDiff_gem_cut->Fill(tDiff);
                    out.h2_3cl_hits_gem_cut->Fill(x1, y1);
                    out.h2_3cl_hits_gem_cut->Fill(x2, y2);
                    out.h2_3cl_hits_gem_cut->Fill(x3, y3);
                    out.h2_3cl_E_angle_gem_cut->Fill(theta1, E1);
                    out.h2_3cl_E_angle_gem_cut->Fill(theta2, E2);
                    out.h2_3cl_E_angle_gem_cut->Fill(theta3, E3);
                    out.h_3cl_yield_gem_cut->Fill(theta1);
                    out.h_3cl_yield_gem_cut->Fill(theta2);
                    out.h_3cl_yield_gem_cut->Fill(theta3);
                    //if (std::isfinite(mass1)) out.h_3cl_mass_gem_cut->Fill(mass1);
                    if (std::isfinite(mass2)) out.h_3cl_mass_gem_cut->Fill(mass2);
                    if (std::isfinite(mass3)) out.h_3cl_mass_gem_cut->Fill(mass3);
                }
            }
        }
    }
    out.processed += n;
    return true;
}

static void mergeResult(QuickResult &dst, const QuickResult &src, fdec::HyCalSystem &hycal)
{
    // X17 three-cluster histograms.
    dst.h_3cl_totalE->Add(src.h_3cl_totalE.get());
    dst.h2_3cl_hits->Add(src.h2_3cl_hits.get());
    dst.h2_3cl_E_angle->Add(src.h2_3cl_E_angle.get());
    dst.h_3cl_E->Add(src.h_3cl_E.get());
    dst.h_3cl_yield->Add(src.h_3cl_yield.get());
    dst.h_3cl_mass->Add(src.h_3cl_mass.get());
    dst.h_3cl_cluster_num->Add(src.h_3cl_cluster_num.get());
    dst.h_3cl_ptx->Add(src.h_3cl_ptx.get());
    dst.h_3cl_pty->Add(src.h_3cl_pty.get());
    dst.h2_3cl_Pt->Add(src.h2_3cl_Pt.get());
    dst.h_3cl_tDiff->Add(src.h_3cl_tDiff.get());

    dst.h_3cl_totalE_cut->Add(src.h_3cl_totalE_cut.get());
    dst.h2_3cl_hits_cut->Add(src.h2_3cl_hits_cut.get());
    dst.h2_3cl_E_angle_cut->Add(src.h2_3cl_E_angle_cut.get());
    dst.h_3cl_E_cut->Add(src.h_3cl_E_cut.get());
    dst.h_3cl_yield_cut->Add(src.h_3cl_yield_cut.get());
    dst.h_3cl_mass_cut->Add(src.h_3cl_mass_cut.get());
    dst.h_3cl_ptx_cut->Add(src.h_3cl_ptx_cut.get());
    dst.h_3cl_pty_cut->Add(src.h_3cl_pty_cut.get());
    dst.h2_3cl_Pt_cut->Add(src.h2_3cl_Pt_cut.get());
    dst.h_3cl_tDiff_cut->Add(src.h_3cl_tDiff_cut.get());

    // X17 three-cluster histograms with GEM matching.
    dst.h_3cl_cluster_gem_num->Add(src.h_3cl_cluster_gem_num.get());
    dst.h_3cl_totalE_gem->Add(src.h_3cl_totalE_gem.get());
    dst.h2_3cl_hits_gem->Add(src.h2_3cl_hits_gem.get());
    dst.h2_3cl_E_angle_gem->Add(src.h2_3cl_E_angle_gem.get());
    dst.h_3cl_E_gem->Add(src.h_3cl_E_gem.get());
    dst.h_3cl_yield_gem->Add(src.h_3cl_yield_gem.get());
    dst.h_3cl_mass_gem->Add(src.h_3cl_mass_gem.get());
    dst.h_3cl_ptx_gem->Add(src.h_3cl_ptx_gem.get());
    dst.h_3cl_pty_gem->Add(src.h_3cl_pty_gem.get());
    dst.h2_3cl_Pt_gem->Add(src.h2_3cl_Pt_gem.get());
    dst.h_3cl_tDiff_gem->Add(src.h_3cl_tDiff_gem.get());

    dst.h_3cl_totalE_gem_cut->Add(src.h_3cl_totalE_gem_cut.get());
    dst.h2_3cl_hits_gem_cut->Add(src.h2_3cl_hits_gem_cut.get());
    dst.h2_3cl_E_angle_gem_cut->Add(src.h2_3cl_E_angle_gem_cut.get());
    dst.h_3cl_E_gem_cut->Add(src.h_3cl_E_gem_cut.get());
    dst.h_3cl_yield_gem_cut->Add(src.h_3cl_yield_gem_cut.get());
    dst.h_3cl_mass_gem_cut->Add(src.h_3cl_mass_gem_cut.get());
    dst.h_3cl_ptx_gem_cut->Add(src.h_3cl_ptx_gem_cut.get());
    dst.h_3cl_pty_gem_cut->Add(src.h_3cl_pty_gem_cut.get());
    dst.h2_3cl_Pt_gem_cut->Add(src.h2_3cl_Pt_gem_cut.get());
    dst.h_3cl_tDiff_gem_cut->Add(src.h_3cl_tDiff_gem_cut.get());
    dst.h_3cl_Vz_gem_cut->Add(src.h_3cl_Vz_gem_cut.get());
    dst.h2_3cl_Vxy_gem_cut->Add(src.h2_3cl_Vxy_gem_cut.get());
    dst.h_3cl_Pt_totalE_cut->Add(src.h_3cl_Pt_totalE_cut.get());
    dst.h_3cl_Pt_totalE_time_cut->Add(src.h_3cl_Pt_totalE_time_cut.get());
    dst.h_3cl_Pt_totalE_time_clusterE_cut->Add(src.h_3cl_Pt_totalE_time_clusterE_cut.get());
    dst.h_3cl_Pt_totalE_time_clusterE_moller_cut->Add(src.h_3cl_Pt_totalE_time_clusterE_moller_cut.get());
    dst.h_3cl_mass_totalE_cut->Add(src.h_3cl_mass_totalE_cut.get());
    dst.h_3cl_mass_totalE_time_cut->Add(src.h_3cl_mass_totalE_time_cut.get());
    dst.h_3cl_mass_totalE_time_clusterE_cut->Add(src.h_3cl_mass_totalE_time_clusterE_cut.get());
    dst.h_3cl_mass_totalE_time_clusterE_Pt_cut->Add(src.h_3cl_mass_totalE_time_clusterE_Pt_cut.get());
    dst.h_3cl_mass_totalE_time_clusterE_Pt_moller_cut->Add(src.h_3cl_mass_totalE_time_clusterE_Pt_moller_cut.get());

    dst.processed += src.processed;
}

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 2108.f;
    int run_id = 12345;
    
    int max_events = -1;
    int num_threads = 4;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:j:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'j': num_threads = std::atoi(optarg); break;
        }
    }
    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]\n";
        return 1;
    }
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(root_files.size())));
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- load run config: assign run_id and Ebeam from gRunConfig ---
    run_id = analysis::get_run_int(root_files[0]);
    gRunConfig = analysis::LoadRunConfig(dbDir + "/runinfo/general.json", run_id);
    Ebeam = gRunConfig.Ebeam > 0.f ? gRunConfig.Ebeam : Ebeam;

    std::cout << "Processing run " << run_id << " with Ebeam = " << Ebeam << " MeV\n";

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    std::cout << "Processing " << root_files.size() << " file(s) with "
              << num_threads << " thread(s)\n";

    std::vector<Long64_t> file_limits(root_files.size(), -1);
    if (max_events > 0) {
        Long64_t remaining = max_events;
        for (size_t i = 0; i < root_files.size(); ++i) {
            Long64_t n = reconEntries(root_files[i]);
            file_limits[i] = std::min(n, remaining);
            remaining -= file_limits[i];
            if (remaining <= 0) {
                for (size_t j = i + 1; j < root_files.size(); ++j)
                    file_limits[j] = 0;
                break;
            }
        }
    }

    std::vector<std::unique_ptr<QuickResult>> results(root_files.size());
    std::atomic<size_t> next_file{0};
    std::atomic<int> errors{0};
    std::mutex io_mtx;
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&]() {
            while (true) {
                size_t idx = next_file.fetch_add(1);
                if (idx >= root_files.size()) break;
                auto res = makeResult(hycal);
                {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "Processing file [" << (idx + 1) << "/"
                              << root_files.size() << "]: " << root_files[idx] << "\n";
                }
                if (!processFile(root_files[idx], file_limits[idx], hycal, Ebeam, *res))
                    ++errors;
                results[idx] = std::move(res);
            }
        });
    }
    for (auto &t : threads) t.join();
    if (errors > 0) return 1;

    auto merged = makeResult(hycal);

    TString outName = output;
    if (outName.IsNull())
        outName = makeDefaultOutput(root_files[0]);
    TFile outfile(outName, "RECREATE");

    // --- write output ---
    outfile.cd();
    outfile.mkdir("x17"); outfile.cd("x17");
    merged->h_3cl_totalE->Write();
    merged->h2_3cl_hits->Write();
    merged->h2_3cl_E_angle->Write();
    merged->h_3cl_E->Write();
    merged->h_3cl_yield->Write();
    merged->h_3cl_mass->Write();
    merged->h_3cl_cluster_num->Write();
    merged->h_3cl_ptx->Write();
    merged->h_3cl_pty->Write();
    merged->h2_3cl_Pt->Write();
    merged->h_3cl_tDiff->Write();

    outfile.cd();
    outfile.mkdir("x17_cut"); outfile.cd("x17_cut");
    merged->h_3cl_totalE_cut->Write();
    merged->h2_3cl_hits_cut->Write();
    merged->h2_3cl_E_angle_cut->Write();
    merged->h_3cl_E_cut->Write();
    merged->h_3cl_yield_cut->Write();
    merged->h_3cl_mass_cut->Write();
    merged->h_3cl_ptx_cut->Write();
    merged->h_3cl_pty_cut->Write();
    merged->h2_3cl_Pt_cut->Write();
    merged->h_3cl_tDiff_cut->Write();

    outfile.cd();
    outfile.mkdir("x17_gem"); outfile.cd("x17_gem");
    merged->h_3cl_cluster_gem_num->Write();
    merged->h_3cl_totalE_gem->Write();
    merged->h2_3cl_hits_gem->Write();
    merged->h2_3cl_E_angle_gem->Write();
    merged->h_3cl_E_gem->Write();
    merged->h_3cl_yield_gem->Write();
    merged->h_3cl_mass_gem->Write();
    merged->h_3cl_ptx_gem->Write();
    merged->h_3cl_pty_gem->Write();
    merged->h2_3cl_Pt_gem->Write();
    merged->h_3cl_tDiff_gem->Write();

    outfile.cd();
    outfile.mkdir("x17_gem_cut"); outfile.cd("x17_gem_cut");
    merged->h_3cl_totalE_gem_cut->Write();
    merged->h2_3cl_hits_gem_cut->Write();
    merged->h2_3cl_E_angle_gem_cut->Write();
    merged->h_3cl_E_gem_cut->Write();
    merged->h_3cl_yield_gem_cut->Write();
    merged->h_3cl_mass_gem_cut->Write();
    merged->h_3cl_ptx_gem_cut->Write();
    merged->h_3cl_pty_gem_cut->Write();
    merged->h2_3cl_Pt_gem_cut->Write();
    merged->h_3cl_tDiff_gem_cut->Write();
    merged->h_3cl_Vz_gem_cut->Write();
    merged->h2_3cl_Vxy_gem_cut->Write();

    outfile.cd();
    outfile.mkdir("x17_gem_cut_steps"); outfile.cd("x17_gem_cut_steps");
    merged->h_3cl_Pt_totalE_cut->Write();
    merged->h_3cl_Pt_totalE_time_cut->Write();
    merged->h_3cl_Pt_totalE_time_clusterE_cut->Write();
    merged->h_3cl_Pt_totalE_time_clusterE_moller_cut->Write();
    merged->h_3cl_mass_totalE_cut->Write();
    merged->h_3cl_mass_totalE_time_cut->Write();
    merged->h_3cl_mass_totalE_time_clusterE_cut->Write();
    merged->h_3cl_mass_totalE_time_clusterE_Pt_cut->Write();
    merged->h_3cl_mass_totalE_time_clusterE_Pt_moller_cut->Write();

    outfile.Close();

    std::cerr << "Result saved -> " << outName.Data() << "\n";
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            std::string name = entry.path().filename().string();
            if (entry.is_regular_file() &&
                name.find("_recon") != std::string::npos &&
                name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

static std::string makeDefaultOutput(const std::string &input_path)
{
    fs::path p(input_path);
    std::string name = p.filename().string();
    const std::string ext = ".root";
    if (name.size() >= ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0) {
        name.insert(name.size() - ext.size(), "_quick_check_x17");
    } else {
        name += "_quick_check_x17.root";
    }
    return (p.parent_path() / name).string();
}
