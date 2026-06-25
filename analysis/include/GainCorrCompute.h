#pragma once
//=============================================================================
// GainCorrCompute.h — shared LMS/alpha gain-correction batch calculation
//=============================================================================

#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>

#include "gain_factor.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace analysis {

static constexpr int   kGainNW        = 1156;
static constexpr int   kGainNLMS      = 3;
static constexpr int   kGainLMSIDBase = 3101;  // 3101=LMS1, 3102=LMS2, 3103=LMS3
static constexpr int   kGainWIDBase   = 1000;  // PbWO4 module_id = W-id + 1000
static constexpr int   kGainHistBins  = 600;
static constexpr float kGainHistMin   = 0.f;
static constexpr float kGainHistMax   = 15000.f;

struct GainBatch {
    int   batch_id         = 0;
    int   event_num_start  = 0;
    int   event_num_end    = 0;
    int   n_lms_events     = 0;
    int   n_alpha_events   = 0;
    int   ref_run          = 0;

    float refPMT_ratio       [kGainNLMS]            = {};
    float gain_W             [kGainNW][kGainNLMS]   = {};
    float gain_W_ref         [kGainNW][kGainNLMS]   = {};
    float gain_corr_W        [kGainNW][kGainNLMS]   = {};
    float fit_mean_ref_lms   [kGainNLMS]            = {};
    float fit_mean_ref_alpha [kGainNLMS]            = {};
    float fit_mean_W_lms     [kGainNW]              = {};
};

struct GainPlotConfig {
    bool          enabled   = false;
    int           max_hists = 10;
    std::set<int> w_ids;
};

struct GainPlotStore {
    std::vector<TH1F*> ref_lms  [kGainNLMS];
    std::vector<TH1F*> ref_alpha[kGainNLMS];
    std::map<int, std::vector<TH1F*>> mod_w;  // key = W-id (1-based)

    ~GainPlotStore();
};

std::string MakeLMSOutputFile(const std::string &evio_path);

void SetupGainBranches(TTree *tree, GainBatch &b);

void FlushGainBatch(GainBatch &b, TTree *tree,
                    TH1F *mod_lms[kGainNW],
                    TH1F *ref_lms[kGainNLMS],
                    TH1F *ref_alpha[kGainNLMS],
                    const prad2::RefGainTable &ref_tbl);

bool ComputeGainCorrections(const std::vector<std::string> &lms_files,
                            const std::string             &gain_out,
                            int                            batch_size,
                            int                            ref_run_num,
                            const prad2::RefGainTable     &ref_tbl,
                            const GainPlotConfig          *plot_cfg = nullptr,
                            GainPlotStore                 *plot_store = nullptr);

} // namespace analysis
