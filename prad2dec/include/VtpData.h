#pragma once
//=============================================================================
// VtpData.h — pre-allocated event data for VTP (0xE122) readout
//
// The VTP trigger processor ships self-describing 32-bit records. PRad-II
// firmware reports trigger-level HyCal clusters as PRAD_CLUSTER records
// (a TAG_EXP expansion, 9-bit tag 0x1CC); the CLAS12-style EC_PEAK /
// EC_CLUSTER records are also decoded in case future firmware emits them.
// Other CLAS12 records are parsed past but not stored.
//=============================================================================

#include <cstdint>
#include <cstddef>

namespace vtp
{

// --- capacity limits --------------------------------------------------------
static constexpr int MAX_EC_PEAKS      = 512;  // per event, across all ROCs
static constexpr int MAX_EC_CLUSTERS   = 64;
static constexpr int MAX_PRAD_CLUSTERS = 256;  // LMS flashes report >100 (139 seen)
static constexpr int MAX_BLOCKS        = 16;   // one per VTP ROC per event

// --- record types -----------------------------------------------------------

// EC_PEAK (XML dict: 0x14)
//   w0: [26] inst, [25:24] view, [23:16] time
//   w1: [25:16] coord, [15:0] energy
struct EcPeak {
    uint32_t roc_tag;    // parent ROC bank tag
    uint8_t  inst;       // 0 or 1 (PCal vs HyCal instance)
    uint8_t  view;       // 0..3 (U/V/W strip view)
    uint8_t  time;       // 8-bit
    uint16_t coord;      // 10-bit strip/xtal coordinate
    uint16_t energy;     // 16-bit energy sum
};

// EC_CLUSTER (XML dict: 0x15)
//   w0: [26] inst, [23:16] time, [15:0] energy
//   w1: [29:20] coordW, [19:10] coordV, [9:0] coordU
struct EcCluster {
    uint32_t roc_tag;
    uint8_t  inst;
    uint8_t  time;
    uint16_t energy;
    uint16_t coordU;     // 10-bit
    uint16_t coordV;
    uint16_t coordW;
};

// PRAD_CLUSTER (TAG_EXP expansion — 9-bit tag bits[31:23] = 0x1CC, written
// "0x10+0x0CC" in $CLON_PARMS/clonbanks.xml).  Trigger-level HyCal cluster
// found by the crate's VTP: seed module, 14-bit energy, hit count, time.
//   w0: [13:0] energy
//   w1: [26:15] id, [14:11] nhits, [10:0] time
// id encodes the seed module: bit 11 set → PbWO4 W(id & 0x7FF) (clonbanks
// displays these as 10001-11156 → W1-1156), else PbGlass G(id) for 1-900.
// Verified against offline HyCal clustering on run 24340: 88% of records
// share the reconstructed cluster's center module, energy correlation 0.99
// (VTP E ≈ 0.75-0.8 × calibrated MeV — trigger-level gains).
struct PradCluster {
    uint32_t roc_tag;    // parent ROC bank tag (which HyCal crate VTP)
    uint16_t energy;     // 14-bit trigger-level energy sum
    uint16_t id;         // 12-bit raw module field
    uint8_t  nhits;      // 4-bit number of modules in the cluster
    uint16_t time;       // 11-bit time within the readout window

    bool     is_pbwo4() const { return (id & 0x800) != 0; }
    uint16_t module()   const { return is_pbwo4() ? (id & 0x7FF) : id; }
};

// VTP block (0x10 BLKHDR + 0x11 BLKTLR pair; 0x12 / 0x13 header data)
struct VtpBlock {
    uint32_t roc_tag;
    uint8_t  slot;            // BLKHDR [26:22]
    uint8_t  module_id;       // BLKHDR may carry module id in upper bits
    uint16_t block_number;    // BLKHDR [17:08]
    uint8_t  block_level;     // BLKHDR [7:0]
    uint32_t nwords;          // BLKTLR [21:0]
    uint32_t event_number;    // EVTHDR [26:0]
    uint64_t trigger_time;    // TRGTIME (48-bit)
    bool     has_trailer;     // set when BLKTLR was seen
    bool     trailer_mismatch; // trailer slot != header slot
};

// --- full event data --------------------------------------------------------
struct VtpEventData {
    int n_peaks         = 0;
    int n_clusters      = 0;
    int n_prad_clusters = 0;
    int n_blocks        = 0;

    EcPeak      peaks[MAX_EC_PEAKS];
    EcCluster   clusters[MAX_EC_CLUSTERS];
    PradCluster prad_clusters[MAX_PRAD_CLUSTERS];
    VtpBlock    blocks[MAX_BLOCKS];

    void clear()
    {
        n_peaks = 0;
        n_clusters = 0;
        n_prad_clusters = 0;
        n_blocks = 0;
    }
};

} // namespace vtp
