#pragma once
//=============================================================================
// RfTime.h — CEBAF RF folding helpers for HyCal timing alignment.
//
// PRad-II reads a downsampled CEBAF RF reference on slot 16 / ch 0 + ch 8
// of ROC 0x40 (see prad2dec/include/TdcData.h).  The recorded leading
// edges are spaced ~131 ns (the downsampler is ×32), but the underlying
// CEBAF RF lattice has period 4.008 ns.  Each recorded edge sits on the
// same 4.008 ns lattice as all the hidden edges, so any per-event time
// difference  (t_ref − t_RF)  can be folded modulo the underlying 4.008 ns
// period to collapse all 32 hidden bunches onto one peak — *provided* the
// per-crystal timing resolution is better than 4.008 ns.
//
// Folded onto (−T_RF/2, T_RF/2] so the natural cut "|Δt| < N·σ" works
// straight off the branch.  See docs/analysis_notes/rf_time_reconstruction_plan.md
// for the full physics summary (R. Demichelis, 2026-05).
//=============================================================================

#include "TdcData.h"

#include <cmath>

namespace prad2 {

// Underlying CEBAF RF period in Hall B (249.5 MHz), as measured for PRad-II
// (R. Demichelis, 2026-05).  The recorded "divided" reference is this
// multiplied by RF_DIVIDER.
static constexpr float RF_PERIOD_NS = 4.008f;
static constexpr int   RF_DIVIDER   = 32;
static constexpr float RF_DIV_NS    = RF_PERIOD_NS * RF_DIVIDER;  // ~128.256

// Fold a raw delta time onto (−T_RF/2, T_RF/2].  NaN passes through
// unchanged so callers can detect "no RF tick this event" without having
// to special-case the branch value.
inline float FoldRfDelta(float dt_ns)
{
    if (!std::isfinite(dt_ns)) return dt_ns;
    return dt_ns - RF_PERIOD_NS * std::nearbyint(dt_ns / RF_PERIOD_NS);
}

// Compute the folded Δt between a reference time `t_ref_ns` (typically
// ClusterHit.time) and the nearest RF tick on channel A (or B if
// `use_b`).  Returns NaN when the requested channel has no hits.
inline float ClusterDeltaRf(float t_ref_ns, const tdc::RfTimeData &rf,
                            bool use_b = false)
{
    const float t_rf = use_b ? rf.nearest_b(t_ref_ns)
                              : rf.nearest_a(t_ref_ns);
    return FoldRfDelta(t_ref_ns - t_rf);
}

} // namespace prad2
