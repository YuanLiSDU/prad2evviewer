#include "VtpDecoder.h"

using namespace vtp;

// Record type codes (value of bits[31:27] on a defining word).
// All have bit 31 = 1, so valid codes are 0x10..0x1F.
namespace {

enum RecordType : uint8_t {
    T_BLKHDR      = 0x10,
    T_BLKTLR      = 0x11,
    T_EVTHDR      = 0x12,
    T_TRGTIME     = 0x13,
    T_EC_PEAK     = 0x14,
    T_EC_CLUSTER  = 0x15,
    T_HTCC        = 0x16,
    T_FT          = 0x17,
    T_FTOF        = 0x18,
    T_CTOF        = 0x19,
    T_CND         = 0x1A,
    T_PCU         = 0x1B,
    T_TAG_EXP     = 0x1C,
    T_TRIGGER     = 0x1D,
    T_DNV         = 0x1E,
    T_FILLER      = 0x1F,
    T_NONE        = 0x00,
};

inline bool is_defining(uint32_t w)       { return (w >> 31) & 0x1; }
inline uint8_t defining_type(uint32_t w)  { return (w >> 27) & 0x1F; }

// TAG_EXP records extend the type with bits [26:23]; the full 9-bit tag is
// bits [31:23].  PRAD_CLUSTER is 0x1CC ("0x10+0x0CC" in clonbanks.xml).
inline uint16_t defining_tag9(uint32_t w) { return (w >> 23) & 0x1FF; }
constexpr uint16_t TAG9_PRAD_CLUSTER = 0x1CC;

} // namespace

int VtpDecoder::DecodeRoc(const uint32_t *data, size_t nwords,
                          uint32_t roc_tag, VtpEventData &evt)
{
    int records_decoded = 0;
    int cur_block = -1;          // index into evt.blocks
    uint8_t cur_type = T_NONE;   // type of the most recent defining word
    int cur_cont = 0;            // 1-based continuation-word counter within a record

    // Saved fields from the defining word of a 2-word record, used when the
    // continuation word arrives.
    uint8_t  peak_inst = 0, peak_view = 0, peak_time = 0;
    uint8_t  cluster_inst = 0, cluster_time = 0;
    uint16_t cluster_energy = 0;
    uint16_t pradcl_energy = 0;
    bool     pradcl_active = false;   // TAG_EXP defining word was PRAD_CLUSTER
    uint32_t trgtime_lo = 0;

    for (size_t i = 0; i < nwords; ++i) {
        uint32_t w = data[i];

        if (is_defining(w)) {
            cur_type = defining_type(w);
            cur_cont = 0;

            switch (cur_type) {

            case T_BLKHDR: {
                if (evt.n_blocks < MAX_BLOCKS) {
                    cur_block = evt.n_blocks++;
                    VtpBlock &b = evt.blocks[cur_block];
                    b.roc_tag        = roc_tag;
                    b.slot           = (w >> 22) & 0x1F;
                    b.module_id      = (w >> 18) & 0xF;
                    b.block_number   = (w >> 8)  & 0x3FF;
                    b.block_level    = w & 0xFF;
                    b.nwords         = 0;
                    b.event_number   = 0;
                    b.trigger_time   = 0;
                    b.has_trailer    = false;
                    b.trailer_mismatch = false;
                } else {
                    cur_block = -1;
                }
                break;
            }

            case T_BLKTLR: {
                if (cur_block >= 0) {
                    VtpBlock &b = evt.blocks[cur_block];
                    uint32_t trailer_nwords = w & 0x3FFFFF;
                    b.nwords       = trailer_nwords;
                    b.has_trailer  = true;
                    uint8_t trailer_slot = (w >> 22) & 0x1F;
                    if (trailer_slot != b.slot) b.trailer_mismatch = true;
                }
                break;
            }

            case T_EVTHDR: {
                if (cur_block >= 0)
                    evt.blocks[cur_block].event_number = w & 0x7FFFFFF;
                break;
            }

            case T_TRGTIME: {
                trgtime_lo = w & 0xFFFFFF;
                if (cur_block >= 0)
                    evt.blocks[cur_block].trigger_time = trgtime_lo;
                break;
            }

            case T_EC_PEAK: {
                peak_inst = (w >> 26) & 0x1;
                peak_view = (w >> 24) & 0x3;
                peak_time = (w >> 16) & 0xFF;
                break;
            }

            case T_EC_CLUSTER: {
                cluster_inst   = (w >> 26) & 0x1;
                cluster_time   = (w >> 16) & 0xFF;
                cluster_energy = w & 0xFFFF;
                break;
            }

            case T_TAG_EXP: {
                // PRAD_CLUSTER w0: [13:0] energy.  Other expansion tags
                // are stepped over.
                pradcl_active = (defining_tag9(w) == TAG9_PRAD_CLUSTER);
                if (pradcl_active)
                    pradcl_energy = w & 0x3FFF;
                break;
            }

            // Types we parse past but don't store (CLAS12/HPS records,
            // PRad TRIGGER summary, DNV, FILLER).
            default:
                break;
            }
        }
        else {
            // Continuation word — dispatch by the last defining type.
            ++cur_cont;

            switch (cur_type) {

            case T_TRGTIME: {
                if (cur_cont == 1 && cur_block >= 0) {
                    uint64_t hi = w & 0xFFFFFF;
                    evt.blocks[cur_block].trigger_time =
                        (hi << 24) | trgtime_lo;
                }
                break;
            }

            case T_EC_PEAK: {
                if (cur_cont == 1 && evt.n_peaks < MAX_EC_PEAKS) {
                    EcPeak &p = evt.peaks[evt.n_peaks++];
                    p.roc_tag = roc_tag;
                    p.inst    = peak_inst;
                    p.view    = peak_view;
                    p.time    = peak_time;
                    p.coord   = (w >> 16) & 0x3FF;
                    p.energy  = w & 0xFFFF;
                    ++records_decoded;
                }
                break;
            }

            case T_EC_CLUSTER: {
                if (cur_cont == 1 && evt.n_clusters < MAX_EC_CLUSTERS) {
                    EcCluster &c = evt.clusters[evt.n_clusters++];
                    c.roc_tag = roc_tag;
                    c.inst    = cluster_inst;
                    c.time    = cluster_time;
                    c.energy  = cluster_energy;
                    c.coordW  = (w >> 20) & 0x3FF;
                    c.coordV  = (w >> 10) & 0x3FF;
                    c.coordU  = w & 0x3FF;
                    ++records_decoded;
                }
                break;
            }

            case T_TAG_EXP: {
                // PRAD_CLUSTER w1: [26:15] id, [14:11] nhits, [10:0] time.
                if (pradcl_active && cur_cont == 1
                    && evt.n_prad_clusters < MAX_PRAD_CLUSTERS) {
                    PradCluster &c = evt.prad_clusters[evt.n_prad_clusters++];
                    c.roc_tag = roc_tag;
                    c.energy  = pradcl_energy;
                    c.id      = (w >> 15) & 0xFFF;
                    c.nhits   = (w >> 11) & 0xF;
                    c.time    = w & 0x7FF;
                    ++records_decoded;
                }
                break;
            }

            // For other types, continuation words are stepped over without
            // interpretation — we stay in the current type until the next
            // defining word flips it.
            default:
                break;
            }
        }
    }

    return records_decoded;
}
