# RF time reconstruction — implementation plan

**Status.** Implemented on branch `RF_recon` (2026-05-28). PR1+PR2+PR3
landed together; vertex-Z / tic-jitter (PR4) deferred. Synthesizes
Raffaella Demichelis's 2026-05 guidance with the existing decoder /
cluster infrastructure. Companion to
[`extract_rf_time.md`](extract_rf_time.md), which covers the end-user
read-side after the changes here.

## Result on run 24840 (first segment, 151 500 events)

| Distribution | mean (ns) | RMS (ns) | N |
|---|---|---|---|
| Unfolded `cl_time − rf_ns_a[0]` | 46.6 | 38.0 | 66 788 |
| Folded mod 4.008 ns, **no** offsets | −0.016 | 1.164 | 66 788 |
| Folded mod 4.008 ns, **84-module** offsets | +0.003 | 1.034 | 66 788 |
| Single PWO W493, no offsets | +0.206 | 1.226 | 3 155 |
| Single PWO W493, with offsets | +0.301 | 0.967 | 3 155 |

Folding alone collapses 38 ns → 1.2 ns (a factor ~32, as expected from
the divider). Per-module offsets pull single-crystal RMS to ≤ 1 ns,
demonstrating that per-PWO timing resolution σ_t < T_RF / 4.
First-pass per-module offsets span (−2, +2) ns with stddev 1.06 ns —
i.e. a substantial chunk of the spread *is* a per-crystal constant,
not the resolution.

---

## 1. What we have, what's missing

| Piece | Status | Location |
|---|---|---|
| 0xE107 TDC bank decoded into per-event hit list | ✅ done | `prad2dec/src/TdcDecoder.cpp` |
| RF-only compact view (`RfTimeData`, leading-edge ns arrays on ch 0 / ch 8) | ✅ done | `prad2dec/include/TdcData.h:96-114`, `TdcDecoder.cpp:61-114` |
| `tdc_roc_tags` / `tdc_nwords` / `tdc_words` on the **raw** events tree | ✅ done | `prad2det/include/EventData.h:168-170` |
| `RfTimeData` + `decode_rf_replay` exposed in Python | ✅ done | `python/bind_dec.cpp:925-967` |
| Per-module HyCal pulse time → `ClusterHit.time` (center module only, ns from sample 0) | ✅ done | `prad2dec/src/WaveAnalyzer.cpp:315`, `prad2det/src/HyCalCluster.cpp:434`, `analysis/src/Replay.cpp:869,890` |
| Per-cluster time on the **recon** tree (`cl_time`) | ✅ done | `EventData.h:188`, `Replay.cpp:910` |
| **TDC / RF bank on the recon tree** | ❌ missing | comment at `EventData.h:243-247` |
| **Per-module HyCal → RF offsets (calibration)** | ❌ missing | no file under `database/` |
| **Per-cluster Δt(HyCal, RF) folded mod 4.008 ns** | ❌ missing | nothing computes it |
| **Calibration tool** that fills the offset table from real data | ❌ missing | — |

So the gap is: a tick-folding rule, a per-module offset table, a calibration
tool that fills it, and a small set of recon-tree branches so downstream
analysis (timing cuts, accidentals subtraction) can use the result.

---

## 2. Physics summary (from Raffaella, 2026-05)

1. The CEBAF RF signal in Hall B has period `T_RF = 4.008 ns` (f = 249.5 MHz).
   What we read on `slot 16 / ch 0` and `ch 8` of ROC `0x40` is downsampled
   by **32** in the front-end, so the recorded leading edges are spaced
   `32 · 4.008 = 128.256 ns` (we see ~131 ns in the data, which is consistent
   given V1190 LSB calibration uncertainty — to be cross-checked against the
   nominal once offsets are in).
2. **The underlying 4.008 ns lattice is what matters.** Each recorded RF edge
   is one of the (hidden) 32× more frequent edges, so for any per-event
   reference time `t` we should fold `(t − t_RF) mod T_RF` — not mod 131 ns.
   That folding collapses all 32 hidden bunches onto a single peak, provided
   the per-crystal timing resolution σ_t < T_RF.
3. RF channel timing resolution is **~30-40 ps** (CLAS12 measurement on the
   same hardware). The HyCal crystal resolution is what we have to *measure*
   — if σ_t < T_RF/2 we'll see a peak after folding; if σ_t > T_RF we'll
   see the flat distribution we currently see on the un-folded ch 0/8 data.
4. **Trigger latency cancels in `(t_cluster − t_RF)`** because both times come
   from the same FADC window / TDC stop reference. So we never need to know
   absolute latency — only the per-crystal offset that brings the folded
   peak to zero.
5. Per-crystal offsets recoverable from RF alignment are **modulo 4.008 ns**.
   Anything larger (cable lengths, FADC-board phase offsets) must be removed
   first by other means (waveform alignment, cosmic timing, etc.) — RF
   calibration is the *last*, fine step.
6. Two future corrections, to be added in a second pass:
   * Vertex-Z correction (path length to each module depends on z_vtx).
   * Jitter correction on the FADC clock tic.
   Both can plug into the same per-cluster Δt formula and don't change the
   structure of this plan.

---

## 3. Design goals

* **Per-cluster Δt_RF folded mod 4.008 ns is a recon-tree branch.** Analysis
  code can apply a timing cut without re-doing the decoding work.
* **Per-module offsets live in a JSON file under `database/`**, mirroring
  the existing `database/hycal_time_cut/example.json` convention so the
  same runinfo wiring (period entries → file path) handles it.
* **The raw RF ticks survive on the recon tree** (small: ≤16 floats per
  channel per event ≈ 128 bytes). This is what lets us re-derive Δt after
  the fact when we change offsets, add ToF corrections, or audit a peak.
* **Calibration is offline.** A new tool under `analysis/tools/` (or a
  ROOT macro under `analysis/scripts/`) eats a recon ROOT file, fills
  per-module folded-Δt histograms, fits the peak, writes JSON. Re-running
  recon then picks up the new offsets.
* **Zero impact on hot path.** Decoding RF and folding Δt is O(hits in
  event) ≤ 12 — already trivial compared to clustering.

---

## 4. Folding rule (the one piece of math)

Given a per-event reference time `t_ref` (typically `ClusterHit.time` in
ns from FADC sample 0) and the nearest RF tick `t_RF` from channel A:

```
raw_dt   = t_ref − t_RF                         # ns, ranges ~[−T_div/2, T_div/2]
folded   = raw_dt − T_RF · round(raw_dt / T_RF) # ns, always in (−T_RF/2, T_RF/2]
```

where `T_RF = 4.008 ns` and `T_div = 32 · T_RF ≈ 128.256 ns` is the recorded
period. The first step pulls `raw_dt` into the recorded bunch (cheap because
`nearest_a()` already does the closest-tick selection); the second step
collapses to a single bunch.

After per-module offsets are applied,

```
delta_t  = folded − offset[center_id]
delta_t  = delta_t − T_RF · round(delta_t / T_RF)   # re-fold to (−T_RF/2, T_RF/2]
```

`delta_t` is what we cut on. The expected shape on a clean elastic sample is
a Gaussian centered at zero with σ ≈ σ_HyCal_time ⊕ σ_RF (RF ~35 ps
negligible).

A small detail: `nearest_a()` returns `NaN` for empty events. For those we
write `delta_t = NaN` and the analysis cut treats `isnan` as fail-cut. No
silent zero.

---

## 5. New / changed code

### 5.1 New header: `prad2det/include/RfTime.h`

Small, self-contained — folding + the constants. Lives in `prad2det` (not
`prad2dec`) so the recon side has the period constant in one place and
`prad2dec` doesn't have to know about clustering.

```cpp
namespace prad2 {

// Underlying CEBAF RF period (Hall B, 2026 measurement: 249.5 MHz).
// The recorded "divided" reference is this multiplied by RF_DIVIDER.
static constexpr float RF_PERIOD_NS  = 4.008f;
static constexpr int   RF_DIVIDER    = 32;
static constexpr float RF_DIV_NS     = RF_PERIOD_NS * RF_DIVIDER;  // ~128.256

// Fold a raw (t_ref − t_RF) onto (−T_RF/2, T_RF/2].  Returns NaN unchanged.
inline float FoldRfDelta(float dt_ns)
{
    if (!std::isfinite(dt_ns)) return dt_ns;
    return dt_ns - RF_PERIOD_NS * std::nearbyint(dt_ns / RF_PERIOD_NS);
}

// Convenience: compute folded Δt from cluster time and the RfTimeData
// channel-A array. Returns NaN if rf is empty.
inline float ClusterDeltaRf(float t_cluster, const tdc::RfTimeData &rf,
                            bool use_b = false)
{
    const float t_rf = use_b ? rf.nearest_b(t_cluster)
                              : rf.nearest_a(t_cluster);
    return FoldRfDelta(t_cluster - t_rf);
}

} // namespace prad2
```

### 5.2 Per-module offset table: `prad2det/include/HyCalRfOffsets.h`

Copy-paste of `HyCalTimeCuts.h` re-stamped for a per-module scalar offset
(ns, in `(−T_RF/2, T_RF/2]`). Same JSON shape so the runinfo loader treats
it identically:

```
{
  "default": 0.0,
  "modules": [
    {"name": "G123", "offset_ns":  1.21},
    {"name": "W735", "offset_ns": -0.74}
  ]
}
```

Public API:

```cpp
struct HyCalRfOffsets {
    float default_off = 0.f;
    std::vector<float> off;        // size = module_count()
    int   n_overrides = 0;

    float at(int module_index) const;     // returns default_off if OOB
    // Already-folded Δt with offset applied + re-fold.
    float apply(int module_index, float folded_dt) const;
};
HyCalRfOffsets LoadHyCalRfOffsets(const std::string &path,
                                  const fdec::HyCalSystem &hycal,
                                  float def_off = 0.f);
```

### 5.3 Wire into `PipelineBuilder` and `RunInfoConfig`

`RunInfoConfig.h` gains a `hycal_rf_offset_file` string (parallel to
`hycal_time_cut_file`); `PipelineBuilder.h` gains a `HyCalRfOffsets
hycal_rf_offsets` member, loaded the same way `hycal_time_cuts` is loaded
(`PipelineBuilder.cpp` — find the existing `LoadHyCalTimeCuts(...)` call
and add a sibling). New runinfo key:

```
"time_cuts": {
    "hc_time_window":    [130.0, 200.0],
    "hycal_module_file": "hycal_time_cut/example.json",
    "hycal_rf_offsets":  "hycal_rf_offsets/24340.json"
}
```

### 5.4 Recon-tree branches: `ReconEventData`

Extend `EventData.h` after the cluster arrays:

```cpp
// RF reference (decoded from 0xE107 ROC 0x40 slot 16 ch 0 / ch 8).
// Same content as tdc::RfTimeData but flattened for ROOT storage.
uint8_t  rf_n_a = 0;                          // ≤ MAX_HITS_PER_CH (16)
uint8_t  rf_n_b = 0;
float    rf_ns_a[tdc::RfTimeData::MAX_HITS_PER_CH] = {};   // leading-edge ns
float    rf_ns_b[tdc::RfTimeData::MAX_HITS_PER_CH] = {};

// Per-cluster RF delta-t, folded to (−T_RF/2, T_RF/2] AFTER per-module
// offset is applied.  NaN when the event has no RF hits in ch A.
float    cl_dt_rf[kMaxClusters] = {};
```

`EventData_io.h` gets the corresponding `tree->Branch(...)` /
`SetBranchAddress` lines. `rf_ns_a/b` are fixed-size float arrays keyed
to `rf_n_a/b`, same pattern as `veto_peak_time`. The total cost is
≤ 8 + 2·64 + 4·100 = 536 bytes per event for an unmatched-cluster
event (~10× cheaper than today's `ssp_raw`).

### 5.5 Per-event fill in `analysis/src/Replay.cpp`

Two new blocks inside the `ProcessWithRecon` event loop. Both run after
the existing `FormClusters()` / cluster-fill block (`Replay.cpp:903-922`).

**Block A — decode RF once per event:**

```cpp
tdc::RfTimeData rf;
tdc::RfTimeDecoder::DecodeReplay(raw.tdc_roc_tags, raw.tdc_nwords,
                                 raw.tdc_words, rf);
ev->rf_n_a = static_cast<uint8_t>(rf.n_a);
ev->rf_n_b = static_cast<uint8_t>(rf.n_b);
std::copy(rf.ns_a, rf.ns_a + rf.n_a, ev->rf_ns_a);
std::copy(rf.ns_b, rf.ns_b + rf.n_b, ev->rf_ns_b);
```

**Block B — per-cluster folded Δt with offset:**

```cpp
const auto &rf_off = pipeline.hycal_rf_offsets;
for (int i = 0; i < ev->n_clusters; ++i) {
    const int center_idx = hycal_sys.module_index_by_id(ev->cl_center[i]);
    const float dt0 = prad2::ClusterDeltaRf(ev->cl_time[i], rf);
    ev->cl_dt_rf[i] = rf_off.apply(center_idx, dt0);
}
```

(`module_index_by_id` may need adding to `HyCalSystem`; it's a single-line
lookup since the system already has `module_by_id`.)

### 5.6 Python bindings

`python/bind_det.cpp` exports `HyCalRfOffsets`, `prad2::FoldRfDelta`,
`prad2::ClusterDeltaRf`. The recon-tree branches show up automatically
once `EventData_io.h` is updated.

### 5.7 EVIO viewer

Optional and out of scope for the first PR — the live event viewer can
keep showing the existing RF arrays. A small per-event "Δt(cluster, RF)"
overlay can land later.

---

## 6. Calibration tool

### 6.1 What it does

Reads a recon ROOT file (which already has `cl_time`, `cl_center`,
`rf_ns_a/b` after the changes in §5), histograms folded
`dt0 = ClusterDeltaRf(cl_time, rf)` per `cl_center`, fits a Gaussian
peak per module, and writes `database/hycal_rf_offsets/<runinfo>.json`.

```
analysis/tools/hycal_rf_offset_calib.cpp
```

mirroring the structure of `det_calib.cpp` / `epCalib.cpp`. Invocation:

```
hycal_rf_offset_calib                                    \
    -i prad_024340.*_recon.root                          \
    -o database/hycal_rf_offsets/24340.json              \
    --min-cluster-energy 200                             \
    --require-single-cluster                             \
    [--fit-window 1.5]
```

Cuts that filter out bad timing references:

| Cut | Rationale |
|---|---|
| `n_clusters == 1` (configurable) | Avoid multi-cluster pile-up bias on `cl_time`. |
| `cl_energy > 200 MeV` (configurable) | High-E clusters have the best per-module timing (waveform SNR). |
| `rf_n_a > 0` (mandatory) | Need an RF tick to subtract. |
| Trigger-bit mask (configurable) | Keep only physics triggers we trust the timing of. |

For each module:

1. Fill 1D histogram of `FoldRfDelta(cl_time − nearest_a(cl_time))` with
   range `(−T_RF/2, T_RF/2]` and ~100 bins (40 ps each).
2. Fit Gaussian + constant background in a `±fit_window` ns range around
   the bin maximum.
3. If the fit χ²/ndf is reasonable AND the entries exceed a minimum,
   record `offset = fit_mean`. Otherwise leave at `default=0` and add a
   diagnostic note.
4. Emit the JSON. Also emit a sidecar CSV with `(module_name, entries,
   mean, sigma, chi2_ndf, status)` for quick QC.

### 6.2 Iteration story

* First pass: re-replay the run with the new code (offsets all zero).
  Δt histograms will be wide and possibly bimodal — that tells us
  whether σ_HyCal < T_RF.
* Run the calibration tool, get a JSON.
* Re-replay (or just re-derive in analysis from `rf_ns_a` + `cl_time`)
  with the offsets loaded. Histograms should narrow to single Gaussians
  centered at 0 per module.
* Iterate once if the Gaussian fits look biased by the first pass's wide
  windows.

### 6.3 Sanity: aggregated-detector check

Before fitting per-module, the tool prints two integrated histograms:

* **All clusters, unfolded** (`cl_time − nearest_a`, range `[−T_div/2,
  T_div/2]`) — should look like the current flat distribution.
* **All clusters, folded** (`FoldRfDelta(...)`) — if σ_HyCal < T_RF this
  shows a *single* peak even before per-module offsets, since most
  modules sit in a narrow offset range and the bulk averages out the
  spread.

If the folded plot is still flat the timing resolution isn't there yet
— stop, fix HyCal time reconstruction (CFD vs quadratic, gain-vs-timing,
walk correction), come back.

---

## 7. Validation plan

1. **Folded vs un-folded comparison on a single high-E PWO crystal**
   (Raffaella's suggestion: "select a single crystal" first). The
   calibration tool produces this as one of the QC plots.
2. **Per-sector aggregation**: average folded Δt by HyCal sector
   (Top/Right/Bottom/Left/Center). A coherent offset between sectors
   would point at cable-length differences > T_RF that the modulo
   calibration can't fix — flagged in the QC report.
3. **PWO vs PbGlass split.** Different module types likely have
   different σ_t. The calibration JSON has per-module granularity so
   this is automatic, but the QC summary prints the two distributions.
4. **Stability across runs.** Re-run the tool on a second physics run
   from the same period; cross-correlate the two offset JSONs. RMS of
   the differences should be ≪ σ_HyCal_time.
5. **End-to-end timing cut on elastic e-p**: tighten `|delta_t| < 3σ`
   on the elastic peak in `epCalib`; the survival rate of true elastic
   events vs accidentals goes in the validation note.

---

## 8. Open questions for Rafo

* The recorded period is `~131 ns` in data, vs the nominal `128.256 ns`.
  Is this expected (TDC LSB calibration uncertainty), or should we adjust
  `TDC_LSB_NS` from `0.024` toward something closer to `0.0245`? A few
  events of `(t[i+1] − t[i])` averaged should pin this down. If LSB is
  off, the folding period also shifts proportionally.
* Should we always use channel A for the reference, or take
  `mean(nearest_a, nearest_b)` to reduce single-channel jitter? Two
  channels should be redundant copies; averaging gains ~√2 if they're
  uncorrelated.
* Do we want a *per-cluster* Δt in the recon tree, or *per-module*
  (one Δt for each contributing module, not just the seed)? The latter
  gives access to per-hit timing in the cluster — useful for split-pulse
  studies — at the cost of `kMaxClusters · kMaxHitsPerCluster` floats.
* The vertex-Z and tic-jitter corrections — do they belong in the
  recon-time `cl_dt_rf`, or in an analysis-level helper that uses the
  raw `rf_ns_a/b` arrays plus a reconstructed vertex Z? Probably the
  latter (vertex isn't known at HyCal recon time), but worth confirming.

---

## 9. Implemented PRs (branch `RF_recon`)

| PR | Scope | Status |
|---|---|---|
| 1 | `RfTime.h`, recon-tree `rf_n_a/b`, `rf_ns_a/b[16]`, `Replay.cpp` snapshot + decode, `EventData_io.h` branches, read-status flag | ✅ |
| 2 | `HyCalRfOffsets.h` (mirrors `HyCalTimeCuts.h`), runinfo key `time_cuts.hycal_rf_offsets`, `PipelineBuilder` step 7c, `cl_dt_rf` branch filled post-clustering, Python binding | ✅ |
| 3 | `analysis/tools/hycal_rf_offset_calib.cpp` (Gauss+const fit), `analysis/scripts/rf_time_plot.C` QC macro, first JSON checked in at `database/hycal_rf_offsets/24840.json` (84 modules) | ✅ |
| 4 | Vertex-Z + tic-jitter corrections | deferred |

### Files touched

* New: `prad2det/include/RfTime.h`, `prad2det/include/HyCalRfOffsets.h`,
  `analysis/tools/hycal_rf_offset_calib.cpp`,
  `analysis/scripts/rf_time_plot.C`,
  `database/hycal_rf_offsets/24840.json`.
* Modified: `prad2det/include/EventData.h`,
  `prad2det/include/EventData_io.h`,
  `prad2det/include/RunInfoConfig.h`,
  `prad2det/include/PipelineBuilder.h`,
  `prad2det/src/PipelineBuilder.cpp`,
  `analysis/src/Replay.cpp`,
  `analysis/CMakeLists.txt`,
  `python/bind_det.cpp`,
  `database/runinfo/general.json` (period entry for run 24840),
  `docs/analysis_notes/extract_rf_time.md`.

---

## 10. References

* Email from R. Demichelis, 2026-05.
* [`docs/analysis_notes/extract_rf_time.md`](extract_rf_time.md) — current
  raw-tree access pattern, mirrors §5.1 here on the read side.
* [`prad2dec/include/TdcData.h`](../../prad2dec/include/TdcData.h),
  [`prad2dec/src/TdcDecoder.cpp`](../../prad2dec/src/TdcDecoder.cpp) —
  decoder.
* [`prad2det/include/HyCalTimeCuts.h`](../../prad2det/include/HyCalTimeCuts.h)
  — template the new offset table copies.
* [`analysis/src/Replay.cpp:781-998`](../../analysis/src/Replay.cpp) —
  `ProcessWithRecon` loop where the new fills land.
