#pragma once
//=============================================================================
// TdcData.h — pre-allocated event data for V1190 TDC readout (0xE107)
//
// The 0xE107 bank is emitted by rol2.c after reformatting the raw v1190/v1290
// hardware stream. Each 32-bit word is a single TDC hit packed as:
//
//   bits 31:27 — slot       (5 bits, 0-31)
//   bit  26    — edge       (0 = leading, 1 = trailing)
//   bits 25:19 — channel    (7 bits, 0-127)
//   bits 18:00 — TDC value  (19 bits)
//
// For V1190 boards, rol2.c left-shifts the native 19-bit value by 2 to match
// the v1290's 25 ps LSB, truncating the result back to 19 bits. V1290 data
// enters the 0xE107 stream already in this layout with a different native
// channel width (5 bits) — the XML dictionary documents only the V1190 form
// used in PRad-II.
//
// Tagger readout: 1 VME crate with 3 V1190 boards × 128 channels = 384 ch.
//=============================================================================

#include <cstdint>
#include <cstddef>

namespace tdc
{

// --- TDC tick calibration --------------------------------------------------
// Single source of truth for converting the 19-bit TDC value to ns. Raffaella's
// dedicated calibration gives 23.436 ps per TDC channel; earlier estimates from
// Sergey Boyarinov (2026-05-05) were ~24 ps for the rol2-normalized stream, with
// the underlying V1190 LSB after the V1190→V1290 left-shift at 25 ps. Use
// TDC_LSB_NS rather than embedding a magic constant in analysis code.
static constexpr double TDC_LSB_NS = 23.436e-3; // 23.436 ps per TDC ch, from Raffaella

// --- PRad-II RF reference cabling ------------------------------------------
// Channels 0 and 8 of slot 16 in ROC 0x40 carry the divided CEBAF RF
// reference (run 24386 onward, period ≈ 131.3 ns ≈ 7.61 MHz). Cabling
// constants live here so analysis doesn't have to remember the layout.
static constexpr uint32_t RF_ROC_TAG = 0x40;
static constexpr uint8_t  RF_SLOT    = 16;
static constexpr uint8_t  RF_CH_A    = 0;
static constexpr uint8_t  RF_CH_B    = 8;

// --- capacity limits --------------------------------------------------------
// 3 boards × 128 channels × worst-case multi-hit burst × 2 edges.
// 4096 comfortably covers a saturated tagger event without heap growth.
static constexpr int MAX_TDC_HITS = 4096;

// Optional per-slot/channel index for fast lookup.
static constexpr int MAX_TDC_SLOTS    = 32;   // V1190 slot field is 5 bits
static constexpr int MAX_TDC_CHANNELS = 128;  // V1190 has 128 channels/board

// --- one TDC hit ------------------------------------------------------------
struct TdcHit
{
    uint32_t roc_tag;   // parent ROC bank tag (e.g. 0x008E for the tagger crate)
    uint8_t  slot;      // 5-bit slot, V1190 board position in the VME crate
    uint8_t  channel;   // 7-bit channel, 0-127
    uint8_t  edge;      // 0 = leading, 1 = trailing
    uint32_t value;     // 19-bit TDC value (LSB = 25 ps after rol2 shift)
};

// --- full event data --------------------------------------------------------
struct TdcEventData
{
    int    n_hits = 0;
    TdcHit hits[MAX_TDC_HITS];

    void clear() { n_hits = 0; }

    // Count hits for a given slot (linear scan — fine for tagger-sized events).
    int countSlot(uint8_t slot) const
    {
        int n = 0;
        for (int i = 0; i < n_hits; ++i)
            if (hits[i].slot == slot) ++n;
        return n;
    }

    // Count hits for a given (slot, channel).
    int countChannel(uint8_t slot, uint8_t channel) const
    {
        int n = 0;
        for (int i = 0; i < n_hits; ++i)
            if (hits[i].slot == slot && hits[i].channel == channel) ++n;
        return n;
    }
};

// --- RF time (per event, both reference channels) --------------------------
// Compact, analysis-friendly view of the RF reference: just the leading-edge
// tick lists for the two PRad-II RF channels (RF_CH_A=0 and RF_CH_B=8 of
// slot RF_SLOT in ROC RF_ROC_TAG).  Each channel carries roughly 5–6 hits
// per trigger window (period ≈ 131.3 ns); MAX_HITS_PER_CH=16 covers the
// observed range with margin. Units are ns (TDC_LSB_NS already applied).
struct RfTimeData
{
    static constexpr int MAX_HITS_PER_CH = 16;

    int   n_a = 0;
    int   n_b = 0;
    float ns_a[MAX_HITS_PER_CH] = {};
    float ns_b[MAX_HITS_PER_CH] = {};

    void clear() { n_a = 0; n_b = 0; }

    // Pick the RF tick nearest a reference time (e.g. trigger latency).
    // Returns NaN if the requested channel has no hits in this event.
    float nearest_a(float t_ref_ns) const { return nearest(ns_a, n_a, t_ref_ns); }
    float nearest_b(float t_ref_ns) const { return nearest(ns_b, n_b, t_ref_ns); }

private:
    static float nearest(const float *v, int n, float t_ref);
};

} // namespace tdc
