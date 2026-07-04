//=============================================================================
// online_gain_monitor_base — incremental EVIO→LMS replay plus gain-corr update
//
// This tool is intended for online monitoring.  It preserves intermediate
// *_lms.root files and skips Phase-1 replay for EVIO files whose LMS output
// already exists.  Phase 2 is recomputed from all available LMS files so the
// output gain_corr.root stays deterministic and uses the same calculation as
// replay_gainCorr.
//=============================================================================

#include "Replay.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "GainCorrCompute.h"

#include <TClass.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;

static std::vector<std::string> collectEvioFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (std::filesystem::is_directory(path)) {
        for (auto &e : std::filesystem::directory_iterator(path)) {
            if (e.is_regular_file() &&
                e.path().filename().string().find(".evio") != std::string::npos)
                files.push_back(e.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

static std::vector<std::string> collectLMSFiles(const std::string &path)
{
    std::vector<std::string> files;
    std::error_code ec;
    for (auto &e : std::filesystem::directory_iterator(path, ec)) {
        if (e.is_regular_file() &&
            e.path().filename().string().find("_lms.root") != std::string::npos)
            files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    return files;
}

static bool reportBatchTimeCoverage(const std::string &path)
{
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) {
        std::cerr << "Cannot reopen gain output for time-axis validation: "
                  << path << "\n";
        return false;
    }
    auto *tree = file.Get<TTree>("gain_corr");
    if (!tree || !tree->GetBranch("unix_time")) {
        std::cerr << "Gain output is missing gain_corr/unix_time: " << path << "\n";
        return false;
    }

    uint32_t unix_time = 0;
    tree->SetBranchAddress("unix_time", &unix_time);
    Long64_t valid = 0;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (unix_time > 0) ++valid;
    }
    std::cout << "  Batch time  : " << valid << "/" << tree->GetEntries()
              << " batch(es) have Unix time\n";
    if (tree->GetEntries() > 0 && valid == 0)
        std::cerr << "Warning: no valid EPICS/scaler time anchors; date/time axes will be empty.\n";
    return true;
}

static void usage()
{
    std::cerr <<
        "Usage: prad2ana_online_gain_monitor_base [evio_file_or_dir ...] -w lms_work_dir\n"
        "       [-o gain_corr.root] [-f max_files] [-j threads]\n"
        "       [-c daq_config.json] [-d hycal_map.json]\n"
        "       [-b batch_size (1000)] [-r ref_run] [-R ref_gain_file.dat] [-a]\n"
        "\n"
        "  -w  directory for persistent per-file *_lms.root files (REQUIRED)\n"
        "  -o  output gain_corr.root (default: <db>/gain_factor/gain_correction/prad_RUN_gain_corr.root)\n"
        "  -a  re-analyze existing *_lms.root files only; no EVIO input required\n"
        "  -R  explicit reference gain file; overrides -r/default ref run\n";
}

int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    std::string daq_config = db_dir + "/daq_config.json";
    std::string daq_map;
    std::string work_dir;
    std::string gain_out;
    int max_files   = -1;
    int num_threads = 4;
    int batch_size  = 1000;
    int ref_run     = -1;
    std::string ref_gain_file;
    bool reanalyze_only = false;

    int opt;
    while ((opt = getopt(argc, argv, "aw:o:f:c:d:j:b:r:R:")) != -1) {
        switch (opt) {
            case 'a': reanalyze_only = true; break;
            case 'w': work_dir    = optarg; break;
            case 'o': gain_out    = optarg; break;
            case 'f': max_files   = std::atoi(optarg); break;
            case 'c': daq_config  = optarg; break;
            case 'd': daq_map     = optarg; break;
            case 'j': num_threads = std::atoi(optarg); break;
            case 'b': batch_size  = std::atoi(optarg); break;
            case 'r': ref_run     = std::atoi(optarg); break;
            case 'R': ref_gain_file = optarg; break;
            default:
                usage();
                return 1;
        }
    }

    std::vector<std::string> evio_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectEvioFiles(argv[i]);
        evio_files.insert(evio_files.end(), f.begin(), f.end());
    }
    std::sort(evio_files.begin(), evio_files.end());

    if (work_dir.empty() || (!reanalyze_only && evio_files.empty())) {
        usage();
        return 1;
    }

    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    evio_files.resize(num_files);
    if (num_files > 0)
        num_threads = std::max(1, std::min(num_threads, num_files));

    if (daq_map.empty()) daq_map = db_dir + "/hycal_map.json";

    std::error_code ec;
    std::filesystem::create_directories(work_dir, ec);
    if (ec) {
        std::cerr << "Cannot create work directory " << work_dir
                  << ": " << ec.message() << "\n";
        return 1;
    }

    auto lms_files = collectLMSFiles(work_dir);
    if (reanalyze_only && lms_files.empty()) {
        std::cerr << "No LMS files found in " << work_dir
                  << "; cannot re-analyze.\n";
        return 1;
    }

    int run_num = !evio_files.empty()
        ? get_run_int(evio_files[0])
        : get_run_int(lms_files[0]);
    if (run_num < 0) {
        std::cerr << "Cannot determine run number from input files.\n";
        return 1;
    }
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    if (ref_run < 0) ref_run = gRunConfig.gain_ref_run;

    if (gain_out.empty()) {
        std::string gain_corr_dir = db_dir + "/gain_factor/gain_correction";
        std::filesystem::create_directories(gain_corr_dir, ec);
        if (ec)
            std::cerr << "Warning: cannot create " << gain_corr_dir
                      << ": " << ec.message() << "\n";
        gain_out = gain_corr_dir + "/" +
            std::string(Form("prad_%06d_gain_corr.root", run_num));
    } else {
        auto parent = std::filesystem::path(gain_out).parent_path();
        if (!parent.empty())
            std::filesystem::create_directories(parent, ec);
    }

    std::vector<std::string> todo;
    if (!reanalyze_only) {
        todo.reserve(evio_files.size());
        for (const auto &evio : evio_files) {
            std::string lms = work_dir + "/" + MakeLMSOutputFile(evio);
            if (!std::filesystem::exists(lms))
                todo.push_back(evio);
        }
    }

    std::cout << (reanalyze_only
                  ? "=== Re-analyze existing LMS files ===\n"
                  : "=== Incremental LMS replay ===\n")
              << "  EVIO files : " << evio_files.size() << "\n"
              << "  To replay  : " << todo.size() << "\n"
              << "  Threads    : " << num_threads << "\n"
              << "  Work dir   : " << work_dir << "\n"
              << "  DAQ cfg    : " << daq_config << "\n"
              << "  HyCal      : " << daq_map << "\n";

    std::atomic<int> next_file{0};
    std::atomic<int> replay_errors{0};
    std::mutex io_mtx;

    auto worker = [&]() {
        analysis::Replay replay;
        if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
        replay.LoadHyCalMap(daq_map);

        while (true) {
            int idx = next_file.fetch_add(1);
            if (idx >= static_cast<int>(todo.size())) break;

            std::string out = work_dir + "/" + MakeLMSOutputFile(todo[idx]);
            bool ok = replay.Process_LMSgainFactor(todo[idx], out, db_dir, daq_config);
            std::lock_guard<std::mutex> lk(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << todo.size() << "] "
                          << todo[idx] << " -> " << out << "\n";
            } else {
                std::cerr << "  [" << (idx + 1) << "/" << todo.size() << "] FAILED: "
                          << todo[idx] << "\n";
                ++replay_errors;
            }
        }
    };

    if (!todo.empty()) {
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (int i = 0; i < num_threads; ++i)
            threads.emplace_back(worker);
        for (auto &t : threads) t.join();
    }

    lms_files = collectLMSFiles(work_dir);
    if (lms_files.empty()) {
        std::cerr << "No LMS files found in " << work_dir
                  << "; skipping gain calculation.\n";
        return replay_errors.load() > 0 ? 1 : 0;
    }

    auto ref_tbl = ref_gain_file.empty()
        ? prad2::LoadRefGain(gRunConfig.gain_data_dir + "/ref_gain", ref_run)
        : prad2::LoadRefGainFile(ref_gain_file);
    if (!ref_tbl.loaded) {
        std::cerr << "Warning: reference gain table not loaded"
                  << " (ref_run=" << ref_run
                  << ", dir=" << gRunConfig.gain_data_dir
                  << ", file=" << ref_gain_file << ")\n";
    }
    if (!ref_gain_file.empty() && ref_tbl.run_number >= 0)
        ref_run = ref_tbl.run_number;

    std::cout << "\n=== Gain correction update ===\n"
              << "  LMS files  : " << lms_files.size() << "\n"
              << "  Batch size : " << batch_size << " LMS events\n"
              << "  Ref run    : " << ref_run << "\n"
              << "  Ref file   : " << (ref_gain_file.empty() ? "(auto)" : ref_gain_file) << "\n"
              << "  Output     : " << gain_out << "\n";

    bool ok = ComputeGainCorrections(lms_files, gain_out, batch_size, ref_run, ref_tbl);
    if (ok) ok = reportBatchTimeCoverage(gain_out);
    std::cout << "All done.\n";
    return (!ok || replay_errors.load() > 0) ? 1 : 0;
}
