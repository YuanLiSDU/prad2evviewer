#pragma once
//=============================================================================
// Replay.h — convert raw DAQ data (EVIO) to ROOT trees
//
// Decodes EVIO events and writes per-channel waveform/peak data to a TTree.
// Depends on prad2dec (decoder) and ROOT (TFile/TTree).
//=============================================================================

#include "EvChannel.h"
#include "EventData.h"
#include "WaveAnalyzer.h"
#include "Fadc250FwAnalyzer.h"
#include "DaqConfig.h"
#include "ConfigSetup.h"
#include "load_daq_config.h"

#include <TFile.h>
#include <TTree.h>

#include <string>
#include <unordered_map>

namespace analysis {

// Aliases for the shared replay data structures
using EventVars       = prad2::RawEventData;
using EventVars_Recon = prad2::ReconEventData;
using LMSEventVars    = prad2::LMSEventData;

class Replay
{
public:
    Replay() = default;

    // Load DAQ configuration (event tags, ADC format, etc.).
    void LoadDaqConfig(const std::string &json_path) { evc::load_daq_config(json_path, daq_cfg_); }

    // Load the merged HyCal map.  Populates both the (crate,slot,ch)→name
    // DAQ lookup used by moduleName() and the name→ModuleType lookup used
    // by moduleType().  The "t" field in each record ("PbGlass" / "PbWO4" /
    // "Veto" / "LMS") is the single source of truth for category dispatch;
    // entries without a "daq" block contribute to module_types_ but not
    // daq_map_.  Calling this is strongly recommended — without it every
    // channel returns MOD_UNKNOWN and module_id encoding falls back to
    // HyCal-only conventions.
    void LoadHyCalMap(const std::string &json_path);

    std::string moduleName(int roc, int slot, int ch) const;
    // Returns the prad2::ModuleType enum for the channel, or MOD_UNKNOWN
    // if (a) no DAQ-map entry exists, or (b) hycal_map.json wasn't loaded
    // / doesn't contain this module.
    prad2::ModuleType moduleType(int roc, int slot, int ch) const;
    // Returns the globally-unique module_id (see RawEventData docs):
    //   PbGlass : 1..1156      (G-module ID)
    //   PbWO4   : 1001..2152   (W-module ID + 1000)
    //   VETO    : 3001..3004   (V1..V4)
    //   LMS     : 3100..3103   (LMSPin=3100, LMS1..3=3101..3103)
    // Returns -1 if the module name is unknown.
    int moduleID(int roc, int slot, int ch) const;

    // Convert an EVIO file to a ROOT file with a TTree.
    // max_events <= 0 means process all. peaks=true adds peak branches.
    bool Process(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                 const std::string &db_dir,
                 int max_events = -1, bool write_peaks = false, const std::string &daq_config_file = "");

    bool ProcessWithRecon(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                            const std::string &db_dir,
                            const std::string &daq_config_file = "",
                            const std::string &gem_ped_file = "", float zerosup_override = 0.f,
                            bool prad1 = false);
    
    bool ProcessWithReconX17(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                            const std::string &db_dir,
                            const std::string &daq_config_file = "",
                            const std::string &gem_ped_file = "", float zerosup_override = 0.f);
    
    bool Process_LMSgainFactor(const std::string &input_evio, const std::string &output_root,
                                const std::string &db_dir, const std::string &daq_config_file);

private:
    void setupBranches(TTree *tree, EventVars &ev, bool write_peaks);
    void clearEvent(EventVars &ev);

    void setupReconBranches(TTree *tree, EventVars_Recon &ev, bool x17_mode);
    void clearReconEvent(EventVars_Recon &ev);

    void setupLMSBranches(TTree *tree, LMSEventVars &ev);
    void clearLMSEvent(LMSEventVars &ev);


    using DaqMap = std::unordered_map<std::string, std::string>;  // "roc_slot_ch" -> name
    DaqMap daq_map_;
    // name → ModuleType, populated by LoadHyCalMap().  Empty if the
    // modules JSON wasn't loaded; moduleType() then returns MOD_UNKNOWN.
    std::unordered_map<std::string, prad2::ModuleType> module_types_;
    evc::DaqConfig daq_cfg_;
};

} // namespace analysis
