# Accessing RF time from PRad-II replayed raw files

**Status.** RF time (V1190 TDC, ROC `0x40`) is captured by the
`prad2dec` decoder and written to **both** trees:

* The `events` tree of every `*_raw.root` carries the raw 0xE107 words
  (`tdc_roc_tags` / `tdc_nwords` / `tdc_words`), suitable for re-decoding.
* The `recon` tree of every `*_recon.root` (since 2026-05) carries the
  decoded leading-edge ns arrays plus a per-cluster folded Δt:
  * `rf_n_a` / `rf_n_b` — count of leading edges on CH_A / CH_B
  * `rf_ns_a[16]` / `rf_ns_b[16]` — ns arrays (TDC_LSB_NS already applied)
  * `cl_dt_rf[n_clusters]` — `(cl_time − nearest RF tick)` folded onto
    `(−T_RF/2, T_RF/2]` with per-module HyCal→RF offsets applied
    (NaN when there's no RF tick that event).

The folding rule (Raffaella Demichelis, 2026-05): the recorded RF is a
÷32 of the underlying 249.5 MHz CEBAF lattice, so any `Δt` between a
cluster and an RF tick must be folded **mod 4.008 ns**, not mod 131 ns.
See [`rf_time_reconstruction_plan.md`](rf_time_reconstruction_plan.md)
for the full physics summary and implementation plan.

Files replayed before midnight 2026-05-06 do not have the raw tdc_*
branches; files replayed before 2026-05-28 do not have the recon-tree
`rf_*` / `cl_dt_rf` branches. Re-run `prad2ana_replay_rawdata` (raw) or
`prad2ana_replay_recon` (recon) against the original EVIO to pick them
up.

---

## 1. What's in the replayed raw file

Three branches on the `events` tree, all `vector<uint32_t>`, one entry
per physics event:

| Branch | Meaning |
|---|---|
| `tdc_roc_tags` | Parent ROC tag of each 0xE107 bank. Currently only `[0x40]` (the "rf" ROC) is populated. |
| `tdc_nwords`   | Word count per bank, parallel to `tdc_roc_tags`. |
| `tdc_words`    | Concatenated TDC hit words (one hit per word, V1190 bit layout). |

Each hit word packs:

```
bits 31:27  slot     (5 bits)
bit  26     edge     (0 = leading, 1 = trailing)
bits 25:19  channel  (7 bits)
bits 18:00  TDC      (19 bits, LSB ≈ 24 ps — Sergey, 2026-05-05)
```

The LSB calibration constant is `tdc::TDC_LSB_NS = 0.024` (defined
once in `prad2dec/include/TdcData.h`, mirrored to Python as
`prad2py.dec.TDC_LSB_NS`). Use it rather than embedding the number in
analysis code — if the calibration ever changes you only have to fix
it in one place.

In current PRad-II runs the populated channels are **slot 16 ch 0**
and **slot 16 ch 8** (a divided CEBAF RF reference), with ~6
leading-edge hits per channel per event and a period of ~131.3 ns
(~7.6 MHz). Cabling constants — `RF_ROC_TAG`, `RF_SLOT`, `RF_CH_A`,
`RF_CH_B` — also live in `TdcData.h` and are exposed to Python.

---

## 2. Quick inspection from raw EVIO (no replay needed)

```
evio_dump <prad_NNNNNN.evio.NNNNN> -m rf -n 5
```

Prints the first 5 events with decoded ns + per-channel deltas, then a
per-channel hit-rate / mean-period summary across the whole file.
Useful for sanity-checking new runs before committing to a replay.

---

## 3. Worked example macro (ROOT)

A runnable example lives at `analysis/scripts/rf_time_example.C`. Run
it from the build directory:

```
cd <build_dir>
root -l ../analysis/scripts/rootlogon.C
.x ../analysis/scripts/rf_time_example.C+("PATH_TO_RAW.root", 5)
```

For each of the first N events it prints the per-channel ns arrays
(channels A=0 and B=8) plus three sample `nearest_a()` / `nearest_b()`
lookups so you can see the helper API in action.

---

## 4. C++ analysis snippet

Bind the three vector branches with the usual pointer-to-pointer
dance, then decode per event. The decoder takes care of the bit
layout, the LSB calibration, and the ROC/slot filtering — analysis
code never touches the raw bits.

```cpp
#include "TdcDecoder.h"
#include "EventData.h"
#include "EventData_io.h"

prad2::RawEventData ev;
prad2::SetRawReadBranches(t, ev);

std::vector<uint32_t> *p_roc = &ev.tdc_roc_tags;
std::vector<uint32_t> *p_nw  = &ev.tdc_nwords;
std::vector<uint32_t> *p_w   = &ev.tdc_words;
t->SetBranchAddress("tdc_roc_tags", &p_roc);
t->SetBranchAddress("tdc_nwords",   &p_nw);
t->SetBranchAddress("tdc_words",    &p_w);

tdc::RfTimeData rf;
for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    tdc::RfTimeDecoder::DecodeReplay(*p_roc, *p_nw, *p_w, rf);

    // rf.n_a, rf.ns_a[0..n_a)  — leading-edge ns, channel 0
    // rf.n_b, rf.ns_b[0..n_b)  — leading-edge ns, channel 8
    // rf.nearest_a(t_ref_ns)   — pick the tick closest to a reference
}
```

Need raw ticks or trailing edges as well? Use the lower-level
`TdcDecoder` instead:

```cpp
tdc::TdcEventData hits;
tdc::TdcDecoder::DecodeReplay(*p_roc, *p_nw, *p_w, hits);

for (int i = 0; i < hits.n_hits; ++i) {
    const auto &h = hits.hits[i];
    // h.value (raw ticks), h.slot, h.channel, h.edge, h.roc_tag
}
```

`tdc::TDC_LSB_NS` and the cabling constants
(`RF_ROC_TAG` / `RF_SLOT` / `RF_CH_A` / `RF_CH_B`) come from the same
header; if the cabling moves they only need to change in
`prad2dec/include/TdcData.h`.

---

## 5. Python analysis snippet

The same helpers are exposed via `prad2py`:

```python
import uproot
from prad2py import dec as D

t = uproot.open("PATH_TO_RAW.root")["events"]
arrs = t.arrays(["tdc_roc_tags", "tdc_nwords", "tdc_words"])

rf = D.RfTimeData()
for i in range(len(arrs)):
    D.decode_rf_replay(
        list(map(int, arrs.tdc_roc_tags[i])),
        list(map(int, arrs.tdc_nwords[i])),
        list(map(int, arrs.tdc_words[i])),
        rf,
    )
    a = rf.ns_a                 # numpy float32, length rf.n_a
    b = rf.ns_b                 # numpy float32, length rf.n_b
    # rf.nearest_a(400.0)       # nearest tick to 400 ns
```

For the full TDC view (ticks, edge, slot, channel, ROC tag):

```python
hits = D.TdcEventData()
D.decode_tdc_replay(roc_tags_list, nwords_list, words_list, hits)
arr = hits.hits_numpy
# structured array: roc_tag, slot, channel, edge, _pad, value
```

The same constants are exposed at module level:

```python
D.TDC_LSB_NS, D.RF_ROC_TAG, D.RF_SLOT, D.RF_CH_A, D.RF_CH_B
```

---

## 6. Folded view on the recon tree (since 2026-05-28)

For end-user analysis the canonical entry point is the per-cluster
`cl_dt_rf` branch:

```python
import uproot
arr = uproot.open("PATH_TO_RECON.root")["recon"].arrays(
    ["cl_dt_rf", "cl_energy", "cl_center", "n_clusters", "rf_n_a"])
# Cut: single-cluster events with an RF tick, drop NaN
import numpy as np
mask = (arr["n_clusters"] == 1) & (arr["rf_n_a"] > 0)
dt = arr["cl_dt_rf"][mask, 0]
dt = dt[np.isfinite(dt)]
# Folded onto (-2.004, 2.004) — a Gaussian peak means timing resolution
# is better than the 4.008 ns CEBAF RF period.
```

The C++ macro `analysis/scripts/rf_time_plot.C` produces a four-panel
QC PDF showing the folded distribution vs the unfolded one (mod 131 ns)
and the wide-range unfolded plot with the 4.008 ns lattice highlighted.

The per-module offsets that `cl_dt_rf` already includes are calibrated
by `prad2ana_hycal_rf_offset_calib`:

```
prad2ana_hycal_rf_offset_calib                              \
    -i prad_024840.*_recon.root                             \
    -o database/hycal_rf_offsets/24840.json                 \
    --min-energy 200 --min-entries 100 --fit-window 1.5     \
    --csv rf_offsets_24840.csv
```

Wire the result into runinfo via a period entry:

```
{
    "from_run":  24840,
    "time_cuts": {"hycal_rf_offsets": "hycal_rf_offsets/24840.json"}
}
```

## 7. Known limitations / open questions

* **Bunch-resolved RF requires folding mod 4.008 ns.** The recorded
  signal is the ÷32 divided reference (period ~131 ns); a naive
  `cl_time − nearest_RF` distribution looks ~uniform over 131 ns.
  After mod-4.008 folding the bunch structure becomes a Gaussian peak
  (RMS ~ 1 ns on a single PWO crystal in run 24840), confirming
  per-crystal timing resolution σ_t < T_RF.
* **Only two channels** (slot 16 ch 0 and ch 8), leading edges only.
  No trailing edges, no other slots. Worth checking whether a faster
  RF (or a bunch-tagging signal) is available in the DAQ that we
  should pull into the readout.
* **Offsets recoverable from RF are modulo 4.008 ns.** Larger constant
  offsets (cable lengths, FADC-board phase differences) must come from
  another calibration before the RF fold is applied.

Questions or issues — reply / comment on the corresponding logbook
post. More technical info is in
[`docs/REPLAYED_DATA.md`](../REPLAYED_DATA.md).
