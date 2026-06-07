//=============================================================================
// replay_recon — convert multiple EVIO files to reconstructed ROOT trees (multi-threaded)
//
// Usage: replay_recon <evio_file_or_dir> [more files/dirs...]
//                     -o output_dir [-f max_files] [-n max_events] [-p] [-j num_threads]
//                     [-c daq_config.json] [-d daq_map.json]
//                     [-g gem_pedestal.json] [-z zerosup_threshold]
//   -o  output directory (REQUIRED)
//   -f  max files to process (default: all)
//   -n  max events per file (default: all)
//   -p  read prad1 data and do not include GEM
//   -j  number of threads (default: 4)
//   -c  DAQ configuration file
//   -d  HyCal map file (default: <db>/hycal_map.json)
//   -g  GEM pedestal file
//   -z  zero-suppression threshold override
//=============================================================================

#include "Replay.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <algorithm>
#include <vector>
#include <map>
#include <thread>
#include <atomic>
#include <mutex>

#include <TFileMerger.h>
#include <TClass.h>
#include <TROOT.h>
#include <TH1F.h>
#include "gain_factor.h"

#include <cerrno>
#include <cstring>
#include <sys/wait.h>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;

// Auto-run replay_gainCorr if no gain-correction file exists for this run.
// Returns true when the file is ready; prints a warning and returns false on
// failure (the replay continues with identity gain correction).
static bool ensureGainCorr(int run_num,
                            const std::string &db_dir,
                            const std::vector<std::string> &evio_files,
                            int num_files,
                            const std::string &daq_config,
                            const std::string &daq_map,
                            int num_threads)
{
    std::string gain_corr_dir = db_dir + "/gain_factor/gain_correction";
    if (!prad2::FindGainCorrRootFile(gain_corr_dir, run_num).empty())
        return true;   // already exists

    std::cout << "[gain_corr] No gain-correction file for run " << run_num
              << "; launching prad2ana_replay_gainCorr...\n";

    // Locate the binary beside this executable.
    std::string gainCorr_exe = prad2::module_dir() + "/prad2ana_replay_gainCorr";

    // Temporary directory for intermediate *_lms.root files.
    // replay_gainCorr deletes them (no -s flag); we clean up the dir itself.
    char tmpl[] = "./prad2_lms_XXXXXX";
    char *tmp = mkdtemp(tmpl);
    if (!tmp) {
        std::cerr << "[gain_corr] mkdtemp failed: " << std::strerror(errno) << "\n";
        return false;
    }
    std::string tmp_dir(tmp);

    // Build argv for execvp — no shell, so no command-injection risk.
    std::vector<std::string> arg_strs;
    arg_strs.push_back(gainCorr_exe);
    for (int i = 0; i < num_files; ++i)
        arg_strs.push_back(evio_files[i]);
    arg_strs.push_back("-o"); arg_strs.push_back(tmp_dir);
    arg_strs.push_back("-j"); arg_strs.push_back(std::to_string(num_threads));
    if (!daq_config.empty()) { arg_strs.push_back("-c"); arg_strs.push_back(daq_config); }
    if (!daq_map.empty())    { arg_strs.push_back("-d"); arg_strs.push_back(daq_map); }

    std::vector<char *> argv_vec;
    for (auto &s : arg_strs) argv_vec.push_back(const_cast<char *>(s.c_str()));
    argv_vec.push_back(nullptr);

    pid_t pid = fork();
    if (pid == 0) {
        execvp(gainCorr_exe.c_str(), argv_vec.data());
        std::cerr << "[gain_corr] execvp failed: " << std::strerror(errno) << "\n";
        _exit(1);
    }
    if (pid < 0) {
        std::cerr << "[gain_corr] fork failed: " << std::strerror(errno) << "\n";
        std::filesystem::remove_all(tmp_dir);
        return false;
    }

    int wstatus = 0;
    waitpid(pid, &wstatus, 0);
    std::filesystem::remove_all(tmp_dir);

    if (!WIFEXITED(wstatus) || WEXITSTATUS(wstatus) != 0) {
        std::cerr << "[gain_corr] replay_gainCorr exited with code "
                  << WEXITSTATUS(wstatus) << "\n";
        return false;
    }
    if (prad2::FindGainCorrRootFile(gain_corr_dir, run_num).empty()) {
        std::cerr << "[gain_corr] Output not found in " << gain_corr_dir << "\n";
        return false;
    }
    std::cout << "[gain_corr] Gain-correction file ready.\n";
    return true;
}

static std::vector<std::string> collectEvioFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (std::filesystem::is_directory(path)) {
        for (auto &entry : std::filesystem::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find(".evio") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

static std::string makeOutputFile(const std::string &evio_path)
{
    std::string out = std::filesystem::path(evio_path).filename().string();
    auto pos = out.find(".evio");
    if (pos != std::string::npos)
        out = out.substr(0, pos) + out.substr(pos + 5);
    out += "_recon.root";
    return out;
}

int main(int argc, char *argv[])
{
    // Initialize ROOT for multi-threading
    ROOT::EnableThreadSafety();

    // Force ROOT dictionary initialization in main thread
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string daq_config, daq_map, gem_ped_file, output_dir;
    float zerosup_override = 0.f;
    int max_files = -1;
    int num_threads = 4;
    bool prad1 = false;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:f:c:d:j:g:z:p")) != -1) {
        switch (opt) {
            case 'o': output_dir = optarg; break;
            case 'f': max_files = std::atoi(optarg); break;
            case 'c': daq_config = optarg; break;
            case 'd': daq_map = optarg; break;
            case 'j': num_threads = std::atoi(optarg); break;
            case 'g': gem_ped_file = optarg; break;
            case 'z': zerosup_override = std::atof(optarg); break;
            case 'p': prad1 = true; break;
        }
    }

    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> evio_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectEvioFiles(argv[i]);
        evio_files.insert(evio_files.end(), f.begin(), f.end());
    }

    if (evio_files.empty() || output_dir.empty()) {
        std::cerr << "Usage: replay_recon <evio_file_or_dir> [more files/dirs...] -o output_dir\n"
                  << "       [-f max_files] [-j threads] [-D daq_config.json]\n"
                  << "       [-g gem_ped.json] [-z threshold] [-p]\n";
        std::cerr << "  -o  output directory (REQUIRED)\n";
        std::cerr << "  -f  max files to process (default: all)\n";
        std::cerr << "  -j  number of threads (default: 4)\n";
        std::cerr << "  -D  DAQ config JSON (default: <db>/daq_config.json)\n";
        std::cerr << "  -g  GEM pedestal JSON\n";
        std::cerr << "  -z  zero-suppression threshold override\n";
        std::cerr << "  -p  PRad1 mode (no GEM)\n";
        return 1;
    }
    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::cout << "Processing " << num_files << " files with "
              << num_threads << " threads\n";
    
    if(daq_map.empty()) daq_map = db_dir + "/hycal_map.json";

    // Group files by run number; ensure gain correction for every distinct run.
    {
        std::map<int, std::vector<std::string>> run_files_map;
        for (int i = 0; i < num_files; ++i)
            run_files_map[get_run_int(evio_files[i])].push_back(evio_files[i]);

        std::cout << "Detected " << run_files_map.size() << " run(s):";
        for (auto &[rn, rf] : run_files_map)
            std::cout << " run" << rn << " (" << rf.size() << " file(s))";
        std::cout << "\n";

        for (auto &[rn, rf] : run_files_map)
            ensureGainCorr(rn, db_dir, rf, static_cast<int>(rf.size()),
                           daq_config, daq_map, num_threads);
    }

    int run_num = get_run_int(evio_files[0]);
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    // shared work queue: atomic index into file list
    std::atomic<int> next_file{0};
    std::mutex io_mtx;
    std::atomic<int> errors{0};

    auto worker = [&]() {
        // each thread gets its own Replay instance (own EvChannel, own buffers)
        analysis::Replay replay;
        if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
        replay.LoadHyCalMap(daq_map);
        std::cerr << "Using HyCal map: " << daq_map << "\n";

        while (true) {
            int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            std::string out = output_dir + "/" + makeOutputFile(evio_files[idx]);
            bool ok = replay.ProcessWithRecon(evio_files[idx], out, gRunConfig, db_dir, daq_config,
                                              gem_ped_file, zerosup_override, prad1);

            std::lock_guard<std::mutex> lk(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                          << evio_files[idx] << " -> " << out << "\n";
            } else {
                std::cerr << "  [" << (idx + 1) << "/" << num_files << "] FAILED: "
                          << evio_files[idx] << "\n";
                errors++;
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i)
        threads.emplace_back(worker);
    for (auto &t : threads)
        t.join();

    std::cout << "Done: " << num_files << " files"
              << (errors > 0 ? ", " + std::to_string(errors.load()) + " errors" : "")
              << "\n";

    return errors > 0 ? 1 : 0;
}
