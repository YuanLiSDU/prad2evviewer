//=============================================================================
// sim2replay_hc - convert G4 HC.ModuleEdep + GEM point hits to X17 replay recon
//
// Usage:
//   prad2ana_sim2replay_hc <input.root|dir> [more inputs...]
//       [-o output.root] [-n max_events] [-r run_number]
//       [-m hycal_module.txt] [-s hycal_module_shuffled.dat]
//=============================================================================

#include "EventData.h"
#include "EventData_io.h"
#include "HyCalCluster.h"
#include "InstallPaths.h"
#include "MatchingTools.h"
#include "PipelineBuilder.h"

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

namespace {

constexpr double kGemEdepMinMeV = 1.e-4;  // 100 eV
constexpr double kTargetCenterMm = -4550.0;
constexpr double kCrystalSurfMm = kTargetCenterMm + 7506.464;
constexpr double kHyCalZFromTargetMm = kCrystalSurfMm - kTargetCenterMm;
constexpr double kInvalidModuleEdep = 1.e30;

// Local-max 3x3 VTP-like cluster cuts.  Edit these values by hand as needed.
constexpr float kVtpSeedEnergyMinMeV = 40.f;
constexpr float kVtpSeedEnergyMaxMeV = 10000.f;
constexpr float kVtpTotalEnergyMinMeV = 70.f;
constexpr float kVtpTotalEnergyMaxMeV = 1600.f;
constexpr float kVtpModuleFiredMinEnergyMeV = 8.f;
constexpr int kVtpBlocksMin = 2;
constexpr int kVtpBlocksMax = 9;
constexpr bool kVtpRejectCenter4x4Seeds = true;
constexpr double kVtpCenter4x4AbsXMaxMm = 41.6;
constexpr double kVtpCenter4x4AbsYMaxMm = 41.6;

// Simulated trigger-bit emulation.  Bit numbers match database/trigger_bits.json.
constexpr float kTriggerRawSumTotalEnergyMinMeV = 900.f;
constexpr float kTrigger3ClusterTotalEnergyMinMeV = 1100.f;
constexpr float kTrigger3ClusterTotalEnergyMaxMeV = 10000.f;
constexpr int kTrigger3ClusterMinClusters = 3;

// Smear reconstructed HyCal cluster energy.  Resolution = a / sqrt(E[GeV]).
constexpr bool kHyCalSmearClusterEnergy = true;
constexpr float kHyCalClusterEnergyResolutionA = 0.03f;
constexpr unsigned kHyCalSmearSeed = 12345u;

// G4 DID is built from volume copy numbers:
// layer0 L/R -> 0/1, layer1 L/R -> 2/3.  Replay/general.json uses
// downstream R/L -> 0/1 and upstream R/L -> 2/3.
constexpr int kG4GemDidToReplayDetId[4] = {3, 2, 1, 0};

struct ModuleTableRow {
    int input_index = -1;
    std::string name;
    std::string type;
    double sx = 0.;
    double sy = 0.;
    double x = 0.;
    double y = 0.;
};

struct ModuleMapEntry {
    int input_index = -1;
    int module_index = -1;
    std::string name;
};

struct ModuleMap {
    int table_rows = 0;
    std::vector<ModuleMapEntry> w_modules;
};

long long coordKey(double value)
{
    return static_cast<long long>(std::llround(value * 1000.));
}

std::string moduleKey(double sx, double sy, double x, double y)
{
    std::ostringstream os;
    os << coordKey(sx) << ':' << coordKey(sy) << ':'
       << coordKey(x) << ':' << coordKey(y);
    return os.str();
}

void printUsage(const char *prog)
{
    std::cerr
        << "Usage: " << prog << " <input.root|dir> [more inputs...] "
        << "[-o output.root] [-n max_events] [-r run_number] "
        << "[-m hycal_module.txt] [-s hycal_module_shuffled.dat]\n"
        << "  -o  output ROOT file (default: sim_hc_recon.root)\n"
        << "  -n  maximum events to process (default: all)\n"
        << "  -r  run number for general.json selection (default: -1/latest)\n"
        << "  -m  replay HyCal module table (default: <db>/hycal_module.txt)\n"
        << "  -s  G4 HyCal geometry order file (default: auto-detect "
        << "hycal_module_shuffled.dat)\n";
}

std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (const auto &entry : fs::directory_iterator(path)) {
            if (!entry.is_regular_file()) continue;
            const auto name = entry.path().filename().string();
            if (name.find(".root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

bool hasBranch(TTree &tree, const char *name)
{
    if (tree.GetBranch(name))
        return true;
    auto *branches = tree.GetListOfBranches();
    return branches && branches->FindObject(name);
}

bool requireBranches(TTree &tree, const std::vector<const char *> &names)
{
    bool ok = true;
    for (const char *name : names) {
        if (!hasBranch(tree, name)) {
            std::cerr << "sim2replay_hc: missing input branch '" << name << "'\n";
            ok = false;
        }
    }
    return ok;
}

std::vector<ModuleTableRow> loadReplayModuleRows(const std::string &path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open HyCal module table: " + path);
    }

    std::vector<ModuleTableRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#')
            continue;

        std::istringstream iss(line);
        std::string name, type;
        double sx = 0., sy = 0., sz = 0., x = 0., y = 0., z = 0.;
        if (!(iss >> name >> type >> sx >> sy >> sz >> x >> y >> z))
            continue;

        rows.push_back({
            static_cast<int>(rows.size()), name, type, sx, sy, x, y});
    }

    if (rows.empty())
        throw std::runtime_error("no modules found in " + path);

    return rows;
}

std::string findG4ModuleOrder(const std::vector<std::string> &input_files,
                              const std::string &db_dir,
                              const std::string &module_table)
{
    std::vector<fs::path> candidates;
    if (!module_table.empty())
        candidates.push_back(fs::path(module_table).parent_path() /
                             "hycal_module_shuffled.dat");
    candidates.push_back(fs::path(db_dir) / "hycal_module_shuffled.dat");

    for (const auto &file : input_files) {
        fs::path parent = fs::path(file).parent_path();
        candidates.push_back(parent / "hycal_module_shuffled.dat");
        candidates.push_back(parent / "PRadSim_X17_shiftShielding" /
                             "database" / "hycal_module_shuffled.dat");
    }

    for (const auto &path : candidates) {
        if (!path.empty() && fs::exists(path))
            return path.string();
    }

    for (const auto &file : input_files) {
        fs::path parent = fs::path(file).parent_path();
        if (parent.empty() || !fs::exists(parent))
            continue;
        try {
            for (const auto &entry : fs::recursive_directory_iterator(parent)) {
                if (!entry.is_regular_file())
                    continue;
                if (entry.path().filename() == "hycal_module_shuffled.dat")
                    return entry.path().string();
            }
        } catch (const std::exception &) {
        }
    }

    return {};
}

ModuleMap loadModuleMap(const std::string &module_table,
                        const std::string &g4_module_order,
                        const fdec::HyCalSystem &hycal)
{
    const auto replay_rows = loadReplayModuleRows(module_table);
    ModuleMap out;
    out.table_rows = static_cast<int>(replay_rows.size());

    std::unordered_map<std::string, const ModuleTableRow *> w_by_geometry;
    for (const auto &row : replay_rows) {
        if (row.type != "PbWO4" || row.name.empty() || row.name[0] != 'W')
            continue;
        w_by_geometry.emplace(moduleKey(row.sx, row.sy, row.x, row.y), &row);
    }

    if (!g4_module_order.empty()) {
        std::ifstream in(g4_module_order);
        if (!in) {
            throw std::runtime_error("cannot open G4 module order file: " +
                                     g4_module_order);
        }

        std::string line;
        int row_index = 0;
        while (std::getline(in, line)) {
            const auto first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || line[first] == '#')
                continue;

            std::istringstream iss(line);
            double sx = 0., sy = 0., x = 0., y = 0.;
            if (!(iss >> sx >> sy >> x >> y))
                continue;

            const int input_index = row_index++;
            if (input_index < 576)
                continue;

            auto it = w_by_geometry.find(moduleKey(sx, sy, x, y));
            if (it == w_by_geometry.end()) {
                throw std::runtime_error(
                    "G4 PWO geometry row has no matching replay module in " +
                    module_table + ": input index " +
                    std::to_string(input_index));
            }

            const auto *mod = hycal.module_by_name(it->second->name);
            if (!mod || !mod->is_pwo4()) {
                throw std::runtime_error("W module '" + it->second->name +
                                         "' from " + module_table +
                                         " is not present as PbWO4 in HyCalSystem");
            }
            out.w_modules.push_back({
                input_index, mod->index, it->second->name});
        }

        if (row_index > 0)
            out.table_rows = row_index;
    } else {
        std::cerr
            << "sim2replay_hc: warning: G4 hycal_module_shuffled.dat was not "
            << "found; falling back to replay table row order\n";
        for (const auto &row : replay_rows) {
            if (row.type != "PbWO4" || row.name.empty() || row.name[0] != 'W')
                continue;

            const auto *mod = hycal.module_by_name(row.name);
            if (!mod || !mod->is_pwo4()) {
                throw std::runtime_error("W module '" + row.name +
                                         "' from " + module_table +
                                         " is not present as PbWO4 in HyCalSystem");
            }
            out.w_modules.push_back({row.input_index, mod->index, row.name});
        }
    }

    if (out.w_modules.empty())
        throw std::runtime_error("no PbWO4/W modules found in " + module_table);

    return out;
}

void clearReconEvent(prad2::ReconEventData &ev)
{
    ev = prad2::ReconEventData{};
    std::fill(std::begin(ev.cl_linear_corr), std::end(ev.cl_linear_corr), 1.f);
    std::fill(std::begin(ev.cl_dt_rf), std::end(ev.cl_dt_rf),
              std::numeric_limits<float>::quiet_NaN());
}

uint16_t clampToU16(float value)
{
    if (!std::isfinite(value) || value <= 0.f)
        return 0;
    if (value >= static_cast<float>(std::numeric_limits<uint16_t>::max()))
        return std::numeric_limits<uint16_t>::max();
    return static_cast<uint16_t>(std::lround(value));
}

bool rejectVtpSeedModule(const fdec::Module &module)
{
    if (!kVtpRejectCenter4x4Seeds)
        return false;
    return std::abs(module.x) < kVtpCenter4x4AbsXMaxMm &&
           std::abs(module.y) < kVtpCenter4x4AbsYMaxMm;
}

void fillLocalMaxVtpClusters(const fdec::HyCalSystem &hycal,
                             const std::vector<float> &module_energy,
                             prad2::ReconEventData &ev)
{
    std::vector<int> seeds;
    seeds.reserve(module_energy.size());

    for (int i = 0; i < hycal.module_count(); ++i) {
        const float seed_e = module_energy[i];
        if (seed_e <= 0.f ||
            seed_e < kVtpModuleFiredMinEnergyMeV ||
            seed_e < kVtpSeedEnergyMinMeV ||
            seed_e > kVtpSeedEnergyMaxMeV)
            continue;

        const auto &seed = hycal.module(i);
        if (!seed.is_pwo4())
            continue;
        if (rejectVtpSeedModule(seed))
            continue;

        bool is_local_max = true;
        hycal.for_each_neighbor(i, true, [&](int ni) {
            if (!is_local_max)
                return;
            const float neighbor_e = module_energy[ni];
            if (neighbor_e > seed_e ||
                (neighbor_e == seed_e && neighbor_e > 0.f && ni < i)) {
                is_local_max = false;
            }
        });
        if (is_local_max)
            seeds.push_back(i);
    }

    std::sort(seeds.begin(), seeds.end(), [&](int a, int b) {
        const float ea = module_energy[a];
        const float eb = module_energy[b];
        if (ea != eb)
            return ea > eb;
        return a < b;
    });

    for (int seed_idx : seeds) {
        if (ev.vtp_cl_n >= vtp::MAX_PRAD_CLUSTERS)
            break;

        float e3x3 = 0.f;
        int blocks = 0;

        auto add_if_fired = [&](int module_idx) {
            const float e = module_energy[module_idx];
            if (e <= 0.f || e < kVtpModuleFiredMinEnergyMeV)
                return;
            e3x3 += e;
            ++blocks;
        };

        add_if_fired(seed_idx);
        hycal.for_each_neighbor(seed_idx, true, add_if_fired);

        if (e3x3 < kVtpTotalEnergyMinMeV ||
            e3x3 > kVtpTotalEnergyMaxMeV ||
            blocks < kVtpBlocksMin ||
            blocks > kVtpBlocksMax) {
            continue;
        }

        const auto &seed = hycal.module(seed_idx);
        ev.vtp_cl_time[ev.vtp_cl_n] = 0;
        ev.vtp_cl_energy[ev.vtp_cl_n] = clampToU16(e3x3);
        ev.vtp_cl_center[ev.vtp_cl_n] = static_cast<uint16_t>(seed.id);
        ev.vtp_cl_blocks[ev.vtp_cl_n] = static_cast<uint8_t>(
            std::min(blocks, 255));
        ++ev.vtp_cl_n;
    }
}

uint32_t simulatedTriggerBits(float total_energy_mev, int vtp_cluster_count)
{
    uint32_t bits = 0;
    if (total_energy_mev > kTriggerRawSumTotalEnergyMinMeV)
        bits |= prad2::TBIT_sum;
    if (vtp_cluster_count >= kTrigger3ClusterMinClusters &&
        total_energy_mev > kTrigger3ClusterTotalEnergyMinMeV &&
        total_energy_mev < kTrigger3ClusterTotalEnergyMaxMeV) {
        bits |= prad2::TBIT_3cl;
    }
    return bits;
}

float smearHyCalClusterEnergy(float energy_mev, std::mt19937 &rng)
{
    if (!kHyCalSmearClusterEnergy || energy_mev <= 0.f)
        return energy_mev;

    const float energy_gev = energy_mev / 1000.f;
    const float sigma = energy_mev * kHyCalClusterEnergyResolutionA /
                        std::sqrt(energy_gev);
    std::normal_distribution<float> gaussian(energy_mev, sigma);
    return std::max(0.f, gaussian(rng));
}

int replayGemDetId(int g4_did)
{
    if (g4_did < 0 || g4_did >= 4)
        return -1;
    return kG4GemDidToReplayDetId[g4_did];
}

} // namespace

int main(int argc, char *argv[])
{
    std::string output = "sim_hc_recon.root";
    std::string module_table;
    std::string g4_module_order;
    Long64_t max_events = -1;
    int run_number = -1;

    int opt = 0;
    while ((opt = getopt(argc, argv, "o:n:r:m:s:h")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoll(optarg); break;
            case 'r': run_number = std::atoi(optarg); break;
            case 'm': module_table = optarg; break;
            case 's': g4_module_order = optarg; break;
            case 'h':
            default:
                printUsage(argv[0]);
                return (opt == 'h') ? 0 : 1;
        }
    }

    std::vector<std::string> input_files;
    for (int i = optind; i < argc; ++i) {
        auto files = collectRootFiles(argv[i]);
        input_files.insert(input_files.end(), files.begin(), files.end());
    }

    if (input_files.empty()) {
        printUsage(argv[0]);
        return 1;
    }

    const std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (module_table.empty())
        module_table = db_dir + "/hycal_module.txt";
    if (g4_module_order.empty())
        g4_module_order = findG4ModuleOrder(input_files, db_dir, module_table);

    prad2::Pipeline pipeline;
    try {
        pipeline = prad2::PipelineBuilder()
            .set_database_dir(db_dir)
            .set_run_number(run_number)
            .set_log_stream(&std::cerr)
            .build();
    } catch (const std::exception &e) {
        std::cerr << "sim2replay_hc: detector setup failed: " << e.what() << "\n";
        return 1;
    }

    if (pipeline.hycal.module_count() <= 0) {
        std::cerr << "sim2replay_hc: HyCalSystem has no modules\n";
        return 1;
    }

    ModuleMap module_map;
    try {
        module_map = loadModuleMap(module_table, g4_module_order,
                                   pipeline.hycal);
    } catch (const std::exception &e) {
        std::cerr << "sim2replay_hc: " << e.what() << "\n";
        return 1;
    }

    std::cerr << "[setup] DB dir       : " << db_dir << "\n"
              << "[setup] module table : " << module_table
              << " (" << module_map.table_rows << " rows, "
              << module_map.w_modules.size() << " W modules)\n"
              << "[setup] G4 module map: "
              << (g4_module_order.empty() ? "<replay row order fallback>"
                                          : g4_module_order)
              << "\n"
              << "[setup] run number   : " << run_number << "\n"
              << "[setup] target z     : " << kTargetCenterMm << " mm\n"
              << "[setup] crystal surf : " << kCrystalSurfMm << " mm (G4)\n"
              << "[setup] hycal_z      : " << kHyCalZFromTargetMm
              << " mm (from target)\n"
              << "[setup] VTP 3x3 cuts : seed=[" << kVtpSeedEnergyMinMeV
              << ", " << kVtpSeedEnergyMaxMeV << "] MeV"
              << " total=[" << kVtpTotalEnergyMinMeV
              << ", " << kVtpTotalEnergyMaxMeV << "] MeV"
              << " fired>=" << kVtpModuleFiredMinEnergyMeV << " MeV"
              << " blocks=[" << kVtpBlocksMin << ", " << kVtpBlocksMax
              << "]"
              << " reject_center4x4="
              << (kVtpRejectCenter4x4Seeds ? "on" : "off")
              << " |x|<" << kVtpCenter4x4AbsXMaxMm
              << " |y|<" << kVtpCenter4x4AbsYMaxMm << " mm\n"
              << "[setup] trigger bits  : RawSum if totalE>"
              << kTriggerRawSumTotalEnergyMinMeV
              << " MeV; 3Cluster if vtp_cl_n>="
              << kTrigger3ClusterMinClusters
              << " and totalE=(" << kTrigger3ClusterTotalEnergyMinMeV
              << ", " << kTrigger3ClusterTotalEnergyMaxMeV << ") MeV\n"
              << "[setup] HC E smear    : "
              << (kHyCalSmearClusterEnergy ? "on" : "off")
              << " sigma/E=" << kHyCalClusterEnergyResolutionA
              << "/sqrt(E_GeV), seed=" << kHyCalSmearSeed << "\n"
              << "[setup] GEM DID map   : G4 0,1,2,3 -> replay "
              << kG4GemDidToReplayDetId[0] << ","
              << kG4GemDidToReplayDetId[1] << ","
              << kG4GemDidToReplayDetId[2] << ","
              << kG4GemDidToReplayDetId[3] << "\n";

    TChain chain("T");
    for (const auto &file : input_files) {
        chain.Add(file.c_str());
        std::cerr << "[input] " << file << "\n";
    }

    if (chain.GetEntries() <= 0) {
        std::cerr << "sim2replay_hc: no entries found in input tree 'T'\n";
        return 1;
    }

    if (!requireBranches(chain, {
            "HC.ModuleEdep",
            "GEM.N", "GEM.X", "GEM.Y", "GEM.Z",
            "GEM.Xout", "GEM.Yout", "GEM.Zout",
            "GEM.Edep", "GEM.DID"})) {
        return 1;
    }

    TTreeReader reader(&chain);
    std::unique_ptr<TTreeReaderValue<int>> event_id;
    if (hasBranch(chain, "EventID"))
        event_id = std::make_unique<TTreeReaderValue<int>>(reader, "EventID");

    TTreeReaderArray<double> hc_module_edep(reader, "HC.ModuleEdep");
    TTreeReaderValue<int> gem_n(reader, "GEM.N");
    TTreeReaderArray<double> gem_x(reader, "GEM.X");
    TTreeReaderArray<double> gem_y(reader, "GEM.Y");
    TTreeReaderArray<double> gem_z(reader, "GEM.Z");
    TTreeReaderArray<double> gem_xout(reader, "GEM.Xout");
    TTreeReaderArray<double> gem_yout(reader, "GEM.Yout");
    TTreeReaderArray<double> gem_zout(reader, "GEM.Zout");
    TTreeReaderArray<double> gem_edep(reader, "GEM.Edep");
    TTreeReaderArray<int> gem_did(reader, "GEM.DID");

    std::unique_ptr<TFile> outfile(TFile::Open(output.c_str(), "RECREATE"));
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "sim2replay_hc: cannot create " << output << "\n";
        return 1;
    }

    TTree tree_out("recon", "X17 replay reconstruction from G4 HC.ModuleEdep");
    auto ev = std::make_unique<prad2::ReconEventData>();
    prad2::SetReconWriteBranches(&tree_out, *ev, true);

    fdec::HyCalCluster clusterer(pipeline.hycal);
    clusterer.SetConfig(pipeline.hycal_cluster_cfg);
    analysis::MatchingTools matching;
    matching.SetMatchRange(pipeline.run_cfg.matching_radius);
    matching.SetSquareSelection(pipeline.run_cfg.matching_use_square);
    matching.SetEnergyDependent(pipeline.run_cfg.matching_energy_dependent);
    matching.SetMatchSigma(pipeline.run_cfg.matching_sigma);

    std::vector<float> module_energy(pipeline.hycal.module_count(), 0.f);
    std::mt19937 smear_rng(kHyCalSmearSeed);

    Long64_t processed = 0;
    bool checked_hc_size = false;
    while (reader.Next()) {
        if (max_events >= 0 && processed >= max_events)
            break;

        clearReconEvent(*ev);
        ev->event_num = event_id ? **event_id : static_cast<int>(processed);
        std::fill(module_energy.begin(), module_energy.end(), 0.f);

        clusterer.Clear();
        const auto hc_size = hc_module_edep.GetSize();
        if (!checked_hc_size) {
            checked_hc_size = true;
            if (hc_size < module_map.table_rows) {
                std::cerr << "sim2replay_hc: HC.ModuleEdep has " << hc_size
                          << " entries but module table has "
                          << module_map.table_rows << " rows\n";
                return 1;
            }
        }

        for (const auto &entry : module_map.w_modules) {
            if (entry.input_index < 0 || entry.input_index >= hc_size)
                continue;
            const double edep = hc_module_edep[entry.input_index];
            if (!std::isfinite(edep) || edep <= 0. ||
                edep >= kInvalidModuleEdep)
                continue;
            const float energy = static_cast<float>(edep);
            module_energy[entry.module_index] += energy;
            clusterer.AddHit(entry.module_index, energy, 0.f);
            ev->total_energy += energy;
        }

        fillLocalMaxVtpClusters(pipeline.hycal, module_energy, *ev);
        ev->trigger_bits = simulatedTriggerBits(ev->total_energy,
                                                ev->vtp_cl_n);

        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hc_hits;
        clusterer.ReconstructHits(hc_hits);
        ev->n_clusters = std::min(static_cast<int>(hc_hits.size()),
                                  prad2::kMaxClusters);
        for (int i = 0; i < ev->n_clusters; ++i) {
            const auto &hit = hc_hits[i];
            ev->cl_x[i] = hit.x;
            ev->cl_y[i] = hit.y;
            ev->cl_z[i] = kHyCalZFromTargetMm;
            ev->cl_energy[i] = smearHyCalClusterEnergy(hit.energy,
                                                       smear_rng);
            ev->cl_linear_corr[i] = hit.linear_corr;
            ev->cl_time[i] = hit.time;
            ev->cl_nblocks[i] = static_cast<uint8_t>(
                std::min(hit.nblocks, 255));
            ev->cl_center[i] = static_cast<uint16_t>(hit.center_id);
            ev->cl_flag[i] = hit.flag;
        }

        Long64_t n_gem_arrays = static_cast<Long64_t>(*gem_n);
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_x.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_y.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_z.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_xout.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_yout.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_zout.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_edep.GetSize()));
        n_gem_arrays = std::min(
            n_gem_arrays, static_cast<Long64_t>(gem_did.GetSize()));

        for (Long64_t i = 0; i < n_gem_arrays; ++i) {
            if (gem_edep[i] <= kGemEdepMinMeV)
                continue;
            const int did = replayGemDetId(gem_did[i]);
            if (did < 0)
                continue;
            if (ev->n_gem_hits >= prad2::kMaxGemHits)
                break;

            const int out_idx = ev->n_gem_hits++;
            ev->det_id[out_idx] = static_cast<uint8_t>(did);
            ev->gem_x[out_idx] = static_cast<float>(
                0.5 * (gem_x[i] + gem_xout[i]));
            ev->gem_y[out_idx] = static_cast<float>(
                0.5 * (gem_y[i] + gem_yout[i]));
            ev->gem_z[out_idx] = static_cast<float>(
                gem_z[i] - kTargetCenterMm);
        }

        std::vector<analysis::HCHit> hc_match_hits;
        std::vector<analysis::GEMHit> gem_match_hits[4];
        hc_match_hits.reserve(ev->n_clusters);
        for (int i = 0; i < ev->n_clusters; ++i) {
            hc_match_hits.push_back({
                ev->cl_x[i], ev->cl_y[i], ev->cl_z[i],
                ev->cl_energy[i], ev->cl_center[i], ev->cl_flag[i]});
        }
        for (int i = 0; i < ev->n_gem_hits; ++i) {
            const int d = ev->det_id[i];
            if (d < 0 || d >= 4)
                continue;
            gem_match_hits[d].push_back({
                ev->gem_x[i], ev->gem_y[i], ev->gem_z[i], ev->det_id[i]});
        }

        const auto matched_hits = matching.Match(
            hc_match_hits,
            gem_match_hits[0], gem_match_hits[1],
            gem_match_hits[2], gem_match_hits[3]);
        const auto matched_hits_chamber = matching.MatchPerChamber(
            hc_match_hits,
            gem_match_hits[0], gem_match_hits[1],
            gem_match_hits[2], gem_match_hits[3]);

        for (const auto &m : matched_hits_chamber) {
            const int cl_idx = m.hycal_idx;
            if (cl_idx < 0 || cl_idx >= ev->n_clusters)
                continue;
            for (int d = 0; d < 4; ++d) {
                ev->matchGEMx[cl_idx][d] = m.gem_hits[d][0];
                ev->matchGEMy[cl_idx][d] = m.gem_hits[d][1];
                ev->matchGEMz[cl_idx][d] = m.gem_hits[d][2];
            }
            ev->matchFlag[cl_idx] = m.mflag;
        }

        ev->matchNum = std::min(static_cast<int>(matched_hits.size()),
                                prad2::kMaxClusters);
        for (int i = 0; i < ev->matchNum; ++i) {
            const auto &m = matched_hits[i];
            ev->mHit_E[i] = m.hycal_hit.energy;
            ev->mHit_x[i] = m.hycal_hit.x;
            ev->mHit_y[i] = m.hycal_hit.y;
            ev->mHit_z[i] = m.hycal_hit.z;
            for (int d = 0; d < 2; ++d) {
                ev->mHit_gx[i][d] = m.gem[d].x;
                ev->mHit_gy[i][d] = m.gem[d].y;
                ev->mHit_gz[i][d] = m.gem[d].z;
                ev->mHit_gid[i][d] = m.gem[d].det_id;
            }
        }

        tree_out.Fill();
        ++processed;
        if (processed % 10000 == 0) {
            std::cerr << "\rsim2replay_hc: " << processed
                      << " events processed" << std::flush;
        }
    }

    std::cerr << "\rsim2replay_hc: " << processed
              << " events written to " << output << "\n";
    outfile->cd();
    tree_out.Write();
    outfile->Close();
    return 0;
}
