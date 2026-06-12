#pragma once
//=============================================================================
// EventData.h — Shared data structures for ROOT replay trees
//
// Used by:
//   - analysis/Replay (writer: EVIO → ROOT)
//   - root_data_source (reader: ROOT → viewer)
//   - analysis tools (reader: ROOT → physics analysis)
//
// These structs define the branch layout of ROOT TTrees produced by
// replay_rawdata ("events" tree) and replay_recon ("recon" tree).
// Changing a struct here automatically updates all readers and writers.
//
// NOTE: No ROOT headers needed — uses standard C++ types only.
//       TTree branch setup uses these as plain arrays.
//=============================================================================

#include "Fadc250Data.h"   // MAX_SAMPLES, MAX_PEAKS, MAX_ROCS, MAX_SLOTS
#include "SspData.h"       // MAX_MPDS, MAX_APVS_PER_MPD, APV_STRIP_SIZE, SSP_TIME_SAMPLES
#include "TdcData.h"       // RfTimeData::MAX_HITS_PER_CH

#include <cstdint>
#include <string>
#include <vector>

namespace prad2 {

// ── Slow-control / scaler tree capacity constants ───────────────────────
//
// DSC2 has 16 TRG + 16 TDC channels per slot; we mirror that fixed layout
// in the side tree so analysis can index by channel directly without
// de-serializing a vector.
static constexpr int kDscChannels = 16;

// ── Capacity constants ───────────────────────────────────────────────────

static constexpr int kMaxChannels  = fdec::MAX_ROCS * fdec::MAX_SLOTS * 16;
static constexpr int kMaxGemStrips = ssp::MAX_MPDS * ssp::MAX_APVS_PER_MPD * ssp::APV_STRIP_SIZE;
static constexpr int kMaxClusters  = 100;
static constexpr int kMaxGemHits   = 400;

// ── Front-panel trigger bits ───────────────────────────────────────────────
//
// Single source of truth for the FP trigger-bit masks carried in the
// `trigger_bits` field (multi-bit, from TI master d[5]).  Mind operator
// precedence when testing — the mask must be parenthesised:
//   if ((ev.trigger_bits & TBIT_sum) == 0) continue;   // require sum trigger
static constexpr uint32_t TBIT_sum   = (1u << 8);    // total-energy (sum) trigger
static constexpr uint32_t TBIT_lms   = (1u << 24);   // LMS light-monitoring
static constexpr uint32_t TBIT_alpha = (1u << 25);   // alpha / pulser

// ── Module type categorisation ────────────────────────────────────────────
//
// Single source of truth at the data-tree level.  Values come from the "t"
// field of hycal_map.json (PbGlass / PbWO4 / Veto / LMS), parsed at
// load time and stored per-channel in RawEventData::module_type.  Numeric
// values are arbitrary but stable — kept as a uint8_t so the TTree branch
// stays compact (1 byte per channel).
//
// Consumers categorise channels by this field; they should not parse the
// module name.  An unrecognised type falls through to MOD_UNKNOWN and the
// channel is still written to the tree (so nothing is silently dropped).
enum ModuleType : uint8_t {
    MOD_UNKNOWN = 0,
    MOD_PbGlass = 1,
    MOD_PbWO4   = 2,
    MOD_VETO    = 3,   // Veto scintillators (V1..V4)
    MOD_LMS     = 4,   // LMS reference PMTs (LMSPin, LMS1..3)
};

// ── Raw replay ("events" tree) ───────────────────────────────────────────
//
// One flat FADC250 channel array.  HyCal, Veto, and LMS channels all live
// in the same arrays — distinguish via `module_type[i]`.  Branch prefix
// stays "hycal.*" for backwards compatibility with existing consumers; they
// do `hycal.module_by_id(...)` then `if (!mod || !mod->is_hycal()) continue;`,
// which transparently skips Veto/LMS entries because their module_id values
// (3001..3004 / 3100..3103) are not registered in HyCalSystem.
//
// module_id encoding (globally unique across types):
//   MOD_PbGlass : 1..900       (HyCalSystem G-module IDs; 576 modules, sparse)
//   MOD_PbWO4   : 1001..2156   (HyCal W-module IDs + 1000; 1152 modules, sparse)
//   MOD_VETO    : 3001..3004   (V1..V4)
//   MOD_LMS     : 3100..3103   (LMSPin=3100, LMS1..3 = 3101..3103)
struct RawEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits      = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    // FADC250 per-channel data (HyCal + Veto + LMS, distinguished by module_type).
    // nsamples is uint8 because PTW max ≤ MAX_SAMPLES = 200 (firmware register
    // 0x0007); npeaks (below) is uint8 because MAX_PEAKS = 8.  Both fit
    // comfortably in 8 bits and shave ~9 KB/event vs int.
    int          nch = 0;
    uint16_t     module_id[kMaxChannels]   = {};
    uint8_t      module_type[kMaxChannels] = {};   // ModuleType enum value
    uint8_t      nsamples[kMaxChannels]    = {};
    uint16_t     samples[kMaxChannels][fdec::MAX_SAMPLES] = {};
    float        gain_factor[kMaxChannels] = {};   // 1.0 for non-HyCal types

    // Optional soft-analyzer peak data (gated on -p flag in replay_rawdata).
    //
    // Pedestal-quality fields (ped_nused / ped_quality / ped_slope) are
    // produced by WaveAnalyzer alongside ped_mean / ped_rms — see the
    // Q_PED_* bitmask in Fadc250Data.h for ped_quality semantics.
    float   ped_mean[kMaxChannels]                       = {};
    float   ped_rms[kMaxChannels]                        = {};
    uint8_t ped_nused[kMaxChannels]                      = {};
    uint8_t ped_quality[kMaxChannels]                    = {};
    float   ped_slope[kMaxChannels]                      = {};
    uint8_t npeaks[kMaxChannels]                         = {};
    float   peak_height[kMaxChannels][fdec::MAX_PEAKS]   = {};
    float   peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    float   peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};
    uint8_t peak_quality[kMaxChannels][fdec::MAX_PEAKS]  = {};   // Q_PEAK_* bitmask (currently just Q_PEAK_PILED)

    // Optional firmware-mode (FADC250 Modes 1/2/3) peak data — also gated on -p.
    // Produced by Fadc250FwAnalyzer using the soft pedestal mean as PED.
    //   daq_npeaks       — number of pulses kept (≤ Fadc250FwConfig.MAX_PULSES)
    //   daq_peak_vp      — Vpeak (pedestal-subtracted ADC counts)
    //   daq_peak_integral— Σ over [cross−NSB, cross+NSA] (Mode 2 integral)
    //   daq_peak_time    — interpolated mid-amplitude time (ns)
    //   daq_peak_cross   — Tcross sample index (Mode 1 "first sample number")
    //   daq_peak_pos     — sample index of Vp itself (different from Tcross
    //                      whenever the leading edge spans multiple samples)
    //   daq_peak_coarse  — 4-ns clock index of Vba (10-bit firmware field)
    //   daq_peak_fine    — sub-sample fine bits, 0..63 (62.5 ps LSB)
    //   daq_peak_quality — bitmask: Q_DAQ_PEAK_AT_BOUNDARY|Q_DAQ_NSB_TRUNCATED|
    //                      Q_DAQ_NSA_TRUNCATED|Q_DAQ_VA_OUT_OF_RANGE (see Fadc250Data.h)
    uint8_t daq_npeaks[kMaxChannels] = {};
    float   daq_peak_vp[kMaxChannels][fdec::MAX_PEAKS]       = {};
    float   daq_peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};
    float   daq_peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    int     daq_peak_cross[kMaxChannels][fdec::MAX_PEAKS]    = {};
    int     daq_peak_pos[kMaxChannels][fdec::MAX_PEAKS]      = {};
    int     daq_peak_coarse[kMaxChannels][fdec::MAX_PEAKS]   = {};
    int     daq_peak_fine[kMaxChannels][fdec::MAX_PEAKS]     = {};
    uint8_t daq_peak_quality[kMaxChannels][fdec::MAX_PEAKS]  = {};

    // GEM per-strip data
    int        gem_nch = 0;
    uint8_t mpd_crate[kMaxGemStrips]  = {};
    uint8_t mpd_fiber[kMaxGemStrips]  = {};
    uint8_t apv[kMaxGemStrips]        = {};
    uint8_t strip[kMaxGemStrips]      = {};
    int16_t ssp_samples[kMaxGemStrips][ssp::SSP_TIME_SAMPLES] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;

    // Raw 0xE122 VTP bank words across every VTP ROC in this event.
    // Stored as a flat triple of vectors (parallel arrays) so ROOT can
    // serialize without a custom dictionary for nested vectors:
    //   vtp_roc_tags[i] — parent ROC bank tag of bank i (e.g. 0x96 for
    //                     hycal1vtp, 0x90 for gem1vtp).
    //   vtp_nwords[i]   — number of 32-bit words in bank i.
    //   vtp_words       — concatenated payload, bank i occupies
    //                     vtp_words[off..off+vtp_nwords[i]) where off =
    //                     Σ vtp_nwords[0..i-1].
    // Stored raw so future record-type additions (PRAD_CLUSTER — TAG_EXP
    // 0x1CC, decoded by prad2dec/src/VtpDecoder.cpp — or the still-
    // unspecified PRad TRIGGER 0x1D bit fields) can be re-decoded offline
    // without rerunning the replay.
    std::vector<uint32_t> vtp_roc_tags;
    std::vector<uint32_t> vtp_nwords;
    std::vector<uint32_t> vtp_words;

    // Raw 0xE107 V1190/V1290 TDC bank words across every TDC ROC in this
    // event.  Same flat-triple-of-vectors layout as vtp_*; offset of
    // bank i is Σ tdc_nwords[0..i-1].  In current PRad-II runs there
    // is one TDC ROC per event, ROC 0x40 ("rf"), carrying ~10–12 hits
    // (slot 16, channels 0 and 8 — the divided CEBAF RF reference).
    // Each word packs one hit:
    //   bits 31:27 slot  bit 26 edge (0=lead, 1=trail)
    //   bits 25:19 chan  bits 18:00 TDC value (LSB ≈ 24 ps after rol2
    //                                          V1190→V1290 normalization)
    // Tagger TDC (ROC 0x8E) — when re-enabled — lands in the same
    // branches with its own roc_tag; offline tools split by tag.
    std::vector<uint32_t> tdc_roc_tags;
    std::vector<uint32_t> tdc_nwords;
    std::vector<uint32_t> tdc_words;
};

// ── Reconstructed replay ("recon" tree) ──────────────────────────────────

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    // HyCal clusters
    float total_energy = 0.f;
    int     n_clusters = 0;
    float cl_x[kMaxClusters]       = {};
    float cl_y[kMaxClusters]       = {};
    float cl_z[kMaxClusters]       = {};
    float cl_energy[kMaxClusters]  = {};
    float cl_time[kMaxClusters]    = {};
    uint8_t cl_nblocks[kMaxClusters] = {};
    uint16_t cl_center[kMaxClusters]  = {};
    uint32_t cl_flag[kMaxClusters]    = {};
    // Matching results
    uint32_t matchFlag[kMaxClusters] = {};
    float    matchGEMx[kMaxClusters][4] = {};
    float    matchGEMy[kMaxClusters][4] = {};
    float    matchGEMz[kMaxClusters][4] = {};
    int      matchNum = 0; // number of clusters with matches (for quick access, can be derived from matchFlag)
    //for quick simple access to each matched hit on HC and GEM planes
    // HC_Energy, HC_x/y/z, GEM_x/y/z (in mm, beam center and target center coordinate)
    float    mHit_E[kMaxClusters] = {};
    float    mHit_x[kMaxClusters] = {};
    float    mHit_y[kMaxClusters] = {};
    float    mHit_z[kMaxClusters] = {};
    float    mHit_gx[kMaxClusters][2] = {};
    float    mHit_gy[kMaxClusters][2] = {};
    float    mHit_gz[kMaxClusters][2] = {};
    float    mHit_gid[kMaxClusters][2] = {}; //det_id for matched GEM hits

    // GEM reconstructed hits
    int        n_gem_hits = 0;
    uint8_t det_id[kMaxGemHits]       = {};
    float   gem_x[kMaxGemHits]        = {};
    float   gem_y[kMaxGemHits]        = {};
    float   gem_z[kMaxGemHits]        = {};
    float   gem_x_charge[kMaxGemHits] = {};
    float   gem_y_charge[kMaxGemHits] = {};
    float   gem_x_peak[kMaxGemHits]   = {};
    float   gem_y_peak[kMaxGemHits]   = {};
    uint8_t gem_x_size[kMaxGemHits]   = {};
    uint8_t gem_y_size[kMaxGemHits]   = {};
    uint8_t gem_x_mTbin[kMaxGemHits]   = {};
    uint8_t gem_y_mTbin[kMaxGemHits]   = {};

    //veto information
    int      veto_nch = 0;
    uint8_t veto_id[4]   = {}; // 1,2,3,4 for veto1-4
    int veto_npeaks[4] = {};
    float veto_peak_time[4][fdec::MAX_PEAKS]     = {};
    float veto_peak_height[4][fdec::MAX_PEAKS]   = {};
    float veto_peak_integral[4][fdec::MAX_PEAKS] = {};

    //LMS reference PMT information
    int      lms_nch = 0;
    uint8_t lms_id[4]   = {}; // 1,2,3 for lms 1-3, 0 for lms Pin
    int lms_npeaks[4] = {};
    float lms_peak_time[4][fdec::MAX_PEAKS]     = {};
    float lms_peak_height[4][fdec::MAX_PEAKS]   = {};
    float lms_peak_integral[4][fdec::MAX_PEAKS] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;

    // Raw 0xE122 VTP bank words — same flat triple-of-vectors layout as
    // RawEventData above (bank i occupies vtp_words[off..off+vtp_nwords[i])
    // where off = Σ vtp_nwords[0..i-1]).  Carried on the recon tree so the
    // PRAD_CLUSTER (TAG_EXP 0x1CC — trigger-level cluster, see
    // prad2dec/include/VtpData.h) and still-unspecified TRIGGER (0x1D)
    // payloads can be studied against reconstructed quantities without a
    // co-replayed raw file.  Cheap: PRad-II VTP banks are 3–7 words per ROC.
    std::vector<uint32_t> vtp_roc_tags;
    std::vector<uint32_t> vtp_nwords;
    std::vector<uint32_t> vtp_words;

    // RF reference (decoded from 0xE107 ROC 0x40 slot 16 ch 0 / ch 8).
    // Same content as tdc::RfTimeData but flattened for ROOT storage:
    // rf_ns_a[0..rf_n_a) and rf_ns_b[0..rf_n_b) are the leading-edge ns
    // arrays for the two RF channels (TDC_LSB_NS already applied,
    // trailing edges dropped).  Cheap: ≤16 floats × 2 channels per event
    // ≈ 128 bytes worst case, and most events only fill ~6 entries each.
    // See prad2det/include/RfTime.h for the folding rule and the
    // analysis-side helper that consumes these arrays.
    uint8_t rf_n_a = 0;
    uint8_t rf_n_b = 0;
    float   rf_ns_a[tdc::RfTimeData::MAX_HITS_PER_CH] = {};
    float   rf_ns_b[tdc::RfTimeData::MAX_HITS_PER_CH] = {};

    // Per-cluster RF Δt, folded onto (−T_RF/2, T_RF/2] using channel A as
    // the reference (`prad2::ClusterDeltaRf`).  Per-module offsets from
    // database/hycal_rf_offsets/*.json have already been applied and the
    // result re-folded.  NaN when rf_n_a == 0 for this event.
    float cl_dt_rf[kMaxClusters] = {};
};

// ── Scaler ("scalers" tree) ──────────────────────────────────────────────
//
// One row per DSC2 readout, written when a physics event carrying a SYNC
// flag arrives.  Counts ACCUMULATE from the GO transition; offline code
// gets a windowed live time by diffing two consecutive rows.
//
// Join key for the events tree: `event_number` is the physics event_number
// of the wrapping physics event — analysis can look up scaler rows for any
// physics event N by finding the row with `event_number ≤ N`.
struct RawScalerData {
    int       event_number = 0;       // physics event_number this DSC2 lives inside
    long long ti_ticks     = 0;       // 48-bit TI timestamp of the carrying physics event
    uint32_t  unix_time    = 0;       // most recent SYNC/EPICS unix_time (Sync().unix_time)
    uint32_t  sync_counter = 0;       // most recent SYNC counter
    uint32_t  run_number   = 0;
    uint8_t   trigger_type = 0;       // trigger_type of the carrying physics event

    int       slot         = -1;      // DSC2 physical slot
    uint32_t  gated        = 0;       // selected source: gated   (live)
    uint32_t  ungated      = 0;       // selected source: ungated (total)
    float     live_ratio   = -1.f;    // cached gated/ungated; -1 if ungated == 0
    uint8_t   source       = 0;       // 0=ref, 1=trg, 2=tdc — matches DscScaler::Source
    uint8_t   channel      = 0;       // 0..15 (ignored for ref)

    // Full counter set, kept for diagnostics and so analysis can compute
    // per-channel live time without rerunning the decoder.
    uint32_t ref_gated     = 0;
    uint32_t ref_ungated   = 0;
    uint32_t trg_gated[kDscChannels]   = {};
    uint32_t trg_ungated[kDscChannels] = {};
    uint32_t tdc_gated[kDscChannels]   = {};
    uint32_t tdc_ungated[kDscChannels] = {};
};

// ── EPICS ("epics" tree) ─────────────────────────────────────────────────
//
// One row per EPICS event (top-level tag 0x001F).  Channel readings are
// parallel name/value vectors so analysis can iterate without per-row
// std::pair overhead.  `event_number_at_arrival` is the most-recent
// physics event_number the channel had seen before this EPICS event —
// used as the join key to the events tree.  `ti_ticks_at_arrival` is the
// TI 48-bit tick of that same physics event, captured at decode time so
// the row carries its own anchor on the run timeline (analysis does not
// have to back-look it up via the events tree, which may have dropped
// the anchor event during replay-time filtering).
struct RawEpicsData {
    int      event_number_at_arrival = -1;
    long long ti_ticks_at_arrival    = 0;
    uint32_t unix_time               = 0;
    uint32_t sync_counter            = 0;
    uint32_t run_number              = 0;

    std::vector<std::string> channel;
    std::vector<double>      value;
};

// ── Run-start metadata ("runinfo" tree) ──────────────────────────────────
//
// One row per CODA control event encountered (PRESTART / GO / END), but
// in normal practice only the PRESTART row carries a non-empty
// `daq_config` payload — that's the 0xE10E STRING bank holding the full
// concatenated DAQ configuration text the run was started with.  GO and
// END are emitted with `daq_config` empty so analysis can reconstruct
// the run's start/stop times even when no PRESTART is in the input
// (e.g. processing a non-zero split file in isolation).
//
// `event_tag` distinguishes which control event this row came from:
//   0x11 = PRESTART, 0x12 = GO, 0x14 = END (CODA2 tags; CODA3 0xFFD1/2/4
// are normalized down to these by SyncInfo).
struct RawRunInfo {
    uint32_t    run_number  = 0;
    uint32_t    unix_time   = 0;     // absolute Unix seconds
    uint8_t     run_type    = 0;     // 0 for GO/END; non-zero on PRESTART
    uint8_t     event_tag   = 0;     // 0x11 / 0x12 / 0x14
    std::string daq_config;          // full 0xE10E text (empty on GO/END)
};

// only save the LMS and gain-correction-related fields to the gain tree
struct LMSEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits      = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    int event_type = 0; // 0 for LMS events, 1 for alpha events

    // FADC250 per-channel data — HyCal + LMS in one flat array,
    // distinguished by module_type (MOD_LMS = 4).
    int          nch = 0;
    uint16_t     module_id[kMaxChannels]   = {};
    uint8_t      module_type[kMaxChannels] = {};   // ModuleType enum value

    uint8_t npeaks[kMaxChannels]                         = {};
    float   peak_height[kMaxChannels][fdec::MAX_PEAKS]   = {};
    float   peak_time[kMaxChannels][fdec::MAX_PEAKS]     = {};
    float   peak_integral[kMaxChannels][fdec::MAX_PEAKS] = {};
};

} // namespace prad2
