# Accessing RF time from PRad-II replayed raw files

**Status.** RF time (V1190 TDC, ROC `0x40`) is captured by the
`prad2dec` decoder and written to the `events` tree of every
`*_raw.root` produced by `prad2ana_replay_rawdata`.

The `recon` tree does **not** carry it yet — a reconstructed per-event
RF scalar will land there later (need to work with Rafo on the
reconstruction algorithm). For now, read the RF info from the matching
`*_raw.root` of the same run.

Files replayed before midnight 2026-05-06 do not have these branches —
re-run `prad2ana_replay_rawdata` against the original EVIO to pick them
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

## 6. Known limitations / open questions

* **No bunch-resolved RF in run 24386.** The populated 7.6 MHz signal
  is the coarse divided reference (period ~131 ns). The folded
  `pulse - nearest_RF` distribution is essentially uniform across
  that period, so the data does **not** currently support sub-ns
  time alignment via this signal alone.
* **Only two channels** (slot 16 ch 0 and ch 8), leading edges only.
  No trailing edges, no other slots. Worth checking whether a faster
  RF (or a bunch-tagging signal) is available in the DAQ that we
  should pull into the readout.
* **Recon tree doesn't carry `tdc_*` yet.** Read them from the
  matching `*_raw.root` for now.

Questions or issues — reply / comment on the corresponding logbook
post. More technical info is in
[`docs/REPLAYED_DATA.md`](../REPLAYED_DATA.md).
