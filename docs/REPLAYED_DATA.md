# Replayed Data Trees

Each replay tool writes one main per-event tree, two slow-control side
trees, and a small per-run metadata tree in the same ROOT file:

| Tool | Main tree | Side trees | Contents |
|---|---|---|---|
| `prad2ana_replay_rawdata` | `events` | `scalers`, `epics`, `runinfo` | Raw FADC250 waveforms + GEM strip data + raw VTP / TDC banks + optional per-channel peak analysis |
| `prad2ana_replay_recon`   | `recon`  | `scalers`, `epics`, `runinfo` | HyCal clusters + GEM hits (lab frame) + HyCal↔GEM matches |
| `prad2ana_replay_filter`  | `events` *or* `recon` | `scalers`, `epics` (with `good`), `runinfo` | Subset of the input main tree retained by slow-control cuts; full slow streams concatenated and tagged; `runinfo` passed through verbatim |

`scalers` / `epics` fire on a different cadence than the main tree
(typically once every 1–2 s versus per-trigger), so they share no row
indexing with `events`/`recon`.  Join them by `event_number`
(`event_number_at_arrival` for `epics`) — see the per-tree sections
below.  `runinfo` is run-scoped (one row per CODA control event:
PRESTART / GO / END), with no row-level join key.

# `events` tree (raw)

Written by `prad2ana_replay_rawdata`.  Per-event scalars and per-channel
arrays sized by `hycal.nch` / `gem.nch`.  The `hycal.*` arrays cover **all**
FADC250 channels (HyCal + Veto + LMS) — distinguish via `hycal.module_type`:

| Value | Type | Source |
|---|---|---|
| 0 | Unknown | — |
| 1 | PbGlass | `hycal_map.json` `t="PbGlass"` |
| 2 | PbWO4   | `hycal_map.json` `t="PbWO4"`   |
| 3 | VETO    | `hycal_map.json` `t="Veto"` (V1..V4) |
| 4 | LMS     | `hycal_map.json` `t="LMS"` (LMSPin, LMS1..3) |

`hycal.module_id` encoding (globally unique):
PbGlass = 1..1156, PbWO4 = 1001..2152, VETO = 3001..3004,
LMS = 3100 (Pin) / 3101..3103.

## Branches

### Event header — always written

| Branch | Type | Meaning |
|---|---|---|
| `event_num`    | `int`   | Event number (TI / trigger bank) |
| `trigger_type` | `uint8` | Main trigger type |
| `trigger_bits` | `uint32` | FP trigger inputs (32-bit bitmask) |
| `timestamp`    | `int64` | 48-bit TI timestamp (250 MHz ticks) |
| `ssp_raw`      | `vector<uint32>` | Raw 0xE10C SSP trigger bank words |
| `vtp_roc_tags` | `vector<uint32>` | Parent ROC tag for each 0xE122 VTP bank (see "VTP raw banks" below) |
| `vtp_nwords`   | `vector<uint32>` | Word count per VTP bank, parallel to `vtp_roc_tags` |
| `vtp_words`    | `vector<uint32>` | Concatenated VTP bank payload — bank `i` occupies `vtp_words[Σ vtp_nwords[0..i-1] .. +vtp_nwords[i])` |
| `tdc_roc_tags` | `vector<uint32>` | Parent ROC tag for each 0xE107 TDC bank (see "TDC raw banks" below) |
| `tdc_nwords`   | `vector<uint32>` | Word count per TDC bank, parallel to `tdc_roc_tags` |
| `tdc_words`    | `vector<uint32>` | Concatenated TDC bank payload (one V1190/V1290 hit per word) |

### HyCal / Veto / LMS FADC250 — always written

| Branch | Type | Meaning |
|---|---|---|
| `hycal.nch`         | `int`            | Number of FADC250 channels in event |
| `hycal.module_id`   | `uint16[nch]`    | See ranges above |
| `hycal.module_type` | `uint8[nch]`     | Category enum (legend above) |
| `hycal.nsamples`    | `uint8[nch]`     | Samples per channel (≤ 200) |
| `hycal.samples`     | `uint16[nch][200]` | Raw 12-bit ADC samples |
| `hycal.gain_factor` | `float[nch]`     | Gain correction (1.0 for Veto/LMS) |

### Soft-analyzer outputs — only with `-p`

Local-maxima search with median/MAD-bootstrapped iterative-outlier-rejection
pedestal (`fdec::WaveAnalyzer`).  Up to `MAX_PEAKS = 8` peaks per channel.
Without `-p`, the soft analyzer is skipped entirely (the pedestal estimate
that the firmware analyzer uses is also gated on `-p`).

| Branch | Type | Meaning |
|---|---|---|
| `hycal.ped_mean`       | `float[nch]`      | Pedestal mean (post-rejection) |
| `hycal.ped_rms`        | `float[nch]`      | Pedestal RMS  (post-rejection) |
| `hycal.ped_nused`      | `uint8[nch]`      | # samples surviving outlier rejection |
| `hycal.ped_quality`    | `uint8[nch]`      | `Q_PED_*` bitmask (legend below) |
| `hycal.ped_slope`      | `float[nch]`      | Linear drift across surviving samples (ADC/sample) |
| `hycal.npeaks`         | `uint8[nch]`      | Soft peaks found |
| `hycal.peak_height`    | `float[nch][8]`   | Peak height above pedestal |
| `hycal.peak_time`      | `float[nch][8]`   | Peak time (ns) — quadratic-vertex sub-sample interpolation around the raw peak |
| `hycal.peak_integral`  | `float[nch][8]`   | Peak integral over `[peak.left, peak.right]` (INCLUSIVE) |
| `hycal.peak_quality`   | `uint8[nch][8]`   | `Q_PEAK_*` bitmask: `1` = `Q_PEAK_PILED` (this peak's integration window touches/overlaps an adjacent peak's, within `cfg.peak_pileup_gap` samples) |

`hycal.ped_quality` bits (defined in `prad2dec/include/Fadc250Data.h`):

| Bit | Flag | Meaning |
|---|---|---|
| `0`     | `Q_PED_GOOD`             | clean pedestal — no flags set |
| `1<<0`  | `Q_PED_NOT_CONVERGED`    | `ped_max_iter` exhausted, kept-mask still moving |
| `1<<1`  | `Q_PED_FLOOR_ACTIVE`     | `rms < ped_flatness` — `ped_flatness` was the active band (typical for very quiet channels; informational) |
| `1<<2`  | `Q_PED_TOO_FEW_SAMPLES`  | < 5 samples survived rejection — estimate is unreliable |
| `1<<3`  | `Q_PED_PULSE_IN_WINDOW`  | a peak landed inside the pedestal window we used |
| `1<<4`  | `Q_PED_OVERFLOW`         | a raw sample in the window hit the 12-bit overflow (4095) |
| `1<<5`  | `Q_PED_TRAILING_WINDOW`  | adaptive logic preferred trailing samples over the leading window (informational, not a problem flag) |

A clean event filter is just `hycal.ped_quality == 0`; the
`PULSE_IN_WINDOW` and `TRAILING_WINDOW` flags are useful diagnostics for
events that fail it.

### Firmware (FADC250 Mode 1/2/3) peaks — only with `-p`

Bit-faithful emulation of the JLab FADC250 firmware (Hall-D V3
extensions: NSAT/NPED/MAXPED).  Configured via the
[`fadc250_waveform.firmware`](../database/daq_config.json) block in
`daq_config.json`.  See
[`docs/clas_fadc/FADC250_algorithms.md`](../docs/clas_fadc/FADC250_algorithms.md)
for the algorithm spec.

| Branch | Type | Meaning |
|---|---|---|
| `hycal.daq_npeaks`        | `uint8[nch]`    | Firmware pulses (≤ NPEAK) |
| `hycal.daq_peak_vp`       | `float[nch][8]` | Vpeak (pedestal-subtracted) |
| `hycal.daq_peak_integral` | `float[nch][8]` | Σ over [cross−NSB, cross+NSA] (Mode 2) |
| `hycal.daq_peak_time`     | `float[nch][8]` | Mid-amplitude time, ns (62.5 ps LSB) |
| `hycal.daq_peak_cross`    | `int[nch][8]`   | Tcross sample index |
| `hycal.daq_peak_pos`      | `int[nch][8]`   | Sample of Vpeak |
| `hycal.daq_peak_coarse`   | `int[nch][8]`   | 4-ns clock index of Vba (10-bit) |
| `hycal.daq_peak_fine`     | `int[nch][8]`   | Fine bits 0..63 (6-bit) |
| `hycal.daq_peak_quality`  | `uint8[nch][8]` | `Q_DAQ_*` bitmask: `1` = peak@boundary, `2` = NSB-trunc, `4` = NSA-trunc, `8` = Va out-of-range |

### GEM strips — always written

| Branch | Type | Meaning |
|---|---|---|
| `gem.nch`         | `int`             | Number of GEM strips fired |
| `gem.mpd_crate`   | `uint8[nch]`      | MPD crate ID |
| `gem.mpd_fiber`   | `uint8[nch]`      | MPD fiber ID |
| `gem.apv`         | `uint8[nch]`      | APV ADC channel |
| `gem.strip`       | `uint8[nch]`      | Strip number on the APV |
| `gem.ssp_samples` | `int16[nch][6]`   | 6 SSP time samples per strip |

### VTP raw banks — always written

PRad-II reads up to 9 VTP banks per physics event (`0xE122`): one per
HyCal-side ti_slave (ROCs `0x83 0x85 0x87 0x89 0x8B 0x8D 0x96`) and, on
GEM-bearing runs, one per GEM-side VTP slave (ROCs `0x90 0x91`).  The
banks are stored as raw 32-bit words so any new VTP record type added
by future firmware (PRad TRIGGER `0x1D`, TAG_EXP `0x1C`, EC_PEAK /
EC_CLUSTER for the CLAS12-style decoder, …) can be re-decoded offline
without rerunning the replay.

A flat triple-of-vectors encoding sidesteps the need for a custom ROOT
dictionary on `vector<vector<uint32>>`:

| Branch | Type | Meaning |
|---|---|---|
| `vtp_roc_tags` | `vector<uint32>` | Parent ROC bank tag of each VTP bank (e.g. `0x96` = `hycal1vtp`, `0x90` = `gem1vtp`) |
| `vtp_nwords`   | `vector<uint32>` | Number of 32-bit words in each VTP bank, parallel to `vtp_roc_tags` |
| `vtp_words`    | `vector<uint32>` | Concatenated payload across all VTP banks in this event |

To iterate VTP banks for one event:

```cpp
size_t off = 0;
for (size_t i = 0; i < ev.vtp_roc_tags.size(); ++i) {
    uint32_t roc = ev.vtp_roc_tags[i];
    const uint32_t *p = ev.vtp_words.data() + off;
    size_t n = ev.vtp_nwords[i];
    // ... decode p[0..n) here ...
    off += n;
}
```

Decoder reference: each defining word has bit 31 = 1 with the record
type in bits `[30:27]` (`0x10` = BLKHDR, `0x11` = BLKTLR, `0x12` =
EVTHDR, `0x13` = TRGTIME, `0x14` = EC_PEAK, `0x15` = EC_CLUSTER,
`0x1C` = TAG_EXP, `0x1D` = PRad TRIGGER summary, `0x1F` = FILLER).
PRad-II runs ship EVTHDR + TRGTIME on every HyCal VTP ROC; crates whose
VTP found a trigger-level cluster additionally ship **PRAD_CLUSTER**
(a TAG_EXP expansion with 9-bit tag `[31:23]` = `0x1CC`: w0 `[13:0]`
energy; w1 `[26:15]` seed module id, `[14:11]` nhits, `[10:0]` time —
see `docs/rols/banktags.md` for the module-id encoding) plus a TRIGGER
summary.  `VtpDecoder.cpp` stores PRAD_CLUSTER in
`VtpEventData.prad_clusters`; validated against offline clustering at
0.99 energy correlation.  `evio_dump -m vtp <file>` is the quick
inspector and prints the decoded fields inline.

### TDC raw banks — always written

V1190/V1290 TDC hits arrive in `0xE107` banks (output of `rol2.c` — the
raw `0xE10B` hardware stream is already stripped of TDC headers / EOB
markers / global trailer, so the payload is a flat array of hits).  In
PRad-II run 24386 onward there is one active TDC ROC per physics event,
ROC `0x40` ("rf"), carrying ~10–12 hits/event on slot 16, channels 0
and 8 — the divided CEBAF RF reference (period ≈ 131.3 ns ≈ 7.61 MHz).
The tagger TDC (ROC `0x8E`) lands in the same branches if/when it
comes back online; offline tools split by `tdc_roc_tags[i]`.

Same flat triple-of-vectors encoding as VTP raw banks:

| Branch | Type | Meaning |
|---|---|---|
| `tdc_roc_tags` | `vector<uint32>` | Parent ROC bank tag of each TDC bank (e.g. `0x40` = `rf`, `0x8E` = `tagger`) |
| `tdc_nwords`   | `vector<uint32>` | Number of 32-bit words in each TDC bank, parallel to `tdc_roc_tags` |
| `tdc_words`    | `vector<uint32>` | Concatenated payload — one hit per word |

Per-hit bit layout (all words in `tdc_words`):

| Bits | Field | Meaning |
|---|---|---|
| `[31:27]` | `slot`    | V1190/V1290 board slot (0–31) |
| `[26]`    | `edge`    | 0 = leading, 1 = trailing |
| `[25:19]` | `channel` | Channel within the board (0–127) |
| `[18:00]` | `value`   | TDC tick count (LSB = 23.436 ps, calibrated, after `rol2`'s V1190 → V1290 normalization) |

#### Decoding — use the prad2dec helpers, not raw bit shifts

The bit fields above are an implementation detail; analysis code should
go through `prad2dec` so the calibration constant (`tdc::TDC_LSB_NS =
23.436e-3`) and the PRad-II RF cabling (`tdc::RF_ROC_TAG / RF_SLOT /
RF_CH_A / RF_CH_B`) live in one place.

**All TDC hits** (RF + future tagger), via `tdc::TdcDecoder::DecodeReplay`:

```cpp
#include "TdcDecoder.h"
tdc::TdcEventData hits;
tdc::TdcDecoder::DecodeReplay(ev.tdc_roc_tags, ev.tdc_nwords,
                              ev.tdc_words, hits);
for (int i = 0; i < hits.n_hits; ++i) {
    const tdc::TdcHit &h = hits.hits[i];
    double t_ns = h.value * tdc::TDC_LSB_NS;
    if (h.roc_tag == tdc::RF_ROC_TAG && h.channel == tdc::RF_CH_A)
        ; // ... RF channel A hit at t_ns ...
}
```

**RF only** — compact per-channel ns arrays via
`tdc::RfTimeDecoder::DecodeReplay`:

```cpp
#include "TdcDecoder.h"
tdc::RfTimeData rf;
tdc::RfTimeDecoder::DecodeReplay(ev.tdc_roc_tags, ev.tdc_nwords,
                                 ev.tdc_words, rf);
// rf.n_a / rf.ns_a[k]  — leading-edge times on RF_CH_A (already in ns)
// rf.n_b / rf.ns_b[k]  — leading-edge times on RF_CH_B
float t_rf = rf.nearest_a(/*trigger latency*/ 400.f);   // closest tick
```

Same helpers in Python (`prad2py.dec`):

```python
from prad2py import dec as D
rf = D.RfTimeData()
D.decode_rf_replay(roc_tags, nwords, words, rf)   # vectors from uproot
t_rf = rf.nearest_a(400.0)                        # ns
```

Quick inspector for raw EVIO (no replay needed): `evio_dump -m rf <file>`.

# `recon` tree (reconstructed)

Written by `prad2ana_replay_recon`.  HyCal clustering, GEM hit
reconstruction, and per-cluster HyCal↔GEM matching.  All positions are in
the lab frame (target-centered, beam-aligned, mm).  Trigger filter applied
upstream — only physics events reach the tree.

### Event header

| Branch | Type | Meaning |
|---|---|---|
| `event_num`    | `int`            | Event number |
| `trigger_type` | `uint8`          | Main trigger type |
| `trigger_bits` | `uint32`         | FP trigger bits (32-bit bitmask) |
| `timestamp`    | `int64`          | 48-bit TI timestamp (250 MHz ticks) |
| `total_energy` | `float`          | Σ HyCal cluster energy (MeV) |
| `ssp_raw`      | `vector<uint32>` | Raw 0xE10C SSP trigger bank words |
| `vtp_roc_tags` | `vector<uint32>` | Parent ROC tag for each 0xE122 VTP bank (same encoding as the events tree — see "VTP raw banks" above) |
| `vtp_nwords`   | `vector<uint32>` | Word count per VTP bank, parallel to `vtp_roc_tags` |
| `vtp_words`    | `vector<uint32>` | Concatenated VTP bank payload — bank `i` occupies `vtp_words[Σ vtp_nwords[0..i-1] .. +vtp_nwords[i])` |

The `vtp_*` triple was added to the recon tree in 2026-06 (PRad-II VTP
banks are only 3–7 words per ROC, so the cost is negligible) so the
still-unspecified TRIGGER `0x1D` / TAG_EXP `0x1C` payloads can be studied
against reconstructed quantities; recon files replayed earlier don't have
it.  Raw `tdc_*` words are still NOT carried here — the RF reference they
hold arrives decoded in the `rf_n_a/rf_ns_a/...` branches below; for the
raw hits use a co-replayed `*_raw.root` events tree.

### HyCal clusters

`n_clusters` ≤ 100.  Lab frame (target-centered, beam-aligned, mm).

| Branch | Type | Meaning |
|---|---|---|
| `n_clusters` | `int`                | Number of reconstructed clusters |
| `cl_x`       | `float[n_clusters]`  | Cluster x at HyCal face + shower depth |
| `cl_y`       | `float[n_clusters]`  | Cluster y |
| `cl_z`       | `float[n_clusters]`  | `hycal_z` + shower-depth correction |
| `cl_energy`  | `float[n_clusters]`  | Cluster energy (MeV) |
| `cl_nblocks` | `uint8[n_clusters]`  | Modules summed into the cluster |
| `cl_center`  | `uint16[n_clusters]` | Center module ID (same numbering as `hycal.module_id`) |
| `cl_flag`    | `uint32[n_clusters]` | HyCal cluster flags (bit field) |

### Per-cluster HyCal↔GEM match (all 4 GEMs)

For each HyCal cluster, the closest GEM hit on each of the 4 detectors
within the matching window is recorded (or `0`/`NaN` if none).  Use this
when you want a fixed-shape `[n_clusters][4]` view.

| Branch | Type | Meaning |
|---|---|---|
| `matchFlag` | `uint32[n_clusters]`    | Per-cluster match flags (which GEMs matched) |
| `matchGEMx` | `float[n_clusters][4]`  | Matched GEM x (det 0..3) |
| `matchGEMy` | `float[n_clusters][4]`  | Matched GEM y |
| `matchGEMz` | `float[n_clusters][4]`  | Matched GEM z |

### Quick-access matched pairs (clusters with ≥2 GEMs matched)

`match_num` ≤ 100. Convenient `[match_num][2]` view for analyses that only
care about clusters confirmed on at least two GEM planes.

| Branch | Type | Meaning |
|---|---|---|
| `match_num` | `int`                  | Number of clusters with ≥2 GEMs matched |
| `mHit_E`    | `float[match_num]`     | HyCal cluster energy (MeV) |
| `mHit_x`    | `float[match_num]`     | HyCal cluster x |
| `mHit_y`    | `float[match_num]`     | HyCal cluster y |
| `mHit_z`    | `float[match_num]`     | HyCal cluster z |
| `mHit_gx`   | `float[match_num][2]`  | First 2 matched GEM x |
| `mHit_gy`   | `float[match_num][2]`  | First 2 matched GEM y |
| `mHit_gz`   | `float[match_num][2]`  | First 2 matched GEM z |
| `mHit_gid`  | `float[match_num][2]`  | det_id (0..3) of those 2 GEM hits |

### GEM reconstructed hits

`n_gem_hits` ≤ 400.  All hits across all 4 detectors, lab frame.

| Branch | Type | Meaning |
|---|---|---|
| `n_gem_hits`   | `int`                  | Total GEM hits across all detectors |
| `det_id`       | `uint8[n_gem_hits]`    | GEM detector ID (0..3) |
| `gem_x`        | `float[n_gem_hits]`    | Hit x (lab) |
| `gem_y`        | `float[n_gem_hits]`    | Hit y |
| `gem_z`        | `float[n_gem_hits]`    | Hit z (per-detector plane) |
| `gem_x_charge` | `float[n_gem_hits]`    | Total ADC of the X cluster |
| `gem_y_charge` | `float[n_gem_hits]`    | Total ADC of the Y cluster |
| `gem_x_peak`   | `float[n_gem_hits]`    | Max-strip ADC, X plane |
| `gem_y_peak`   | `float[n_gem_hits]`    | Max-strip ADC, Y plane |
| `gem_x_size`   | `uint8[n_gem_hits]`    | Strips in X cluster |
| `gem_y_size`   | `uint8[n_gem_hits]`    | Strips in Y cluster |
| `gem_x_mTbin`  | `uint8[n_gem_hits]`    | Time-sample bin of max-ADC strip, X |
| `gem_y_mTbin`  | `uint8[n_gem_hits]`    | Time-sample bin of max-ADC strip, Y |

### Veto + LMS (peak summaries)

Lightweight tag of the best soft peak per Veto / LMS channel — full
waveforms live in the `events` tree, not here.

| Branch | Type | Meaning |
|---|---|---|
| `veto_nch`          | `int`              | Number of Veto channels with data |
| `veto_id`           | `uint8[veto_nch]`  | 0..3 (V1..V4) |
| `veto_npeaks`       | `int[veto_nch]`    | Soft peaks found |
| `veto_peak_time`    | `float[veto_nch][8]` | Peak time (ns) |
| `veto_peak_integral`| `float[veto_nch][8]` | Peak integral |
| `lms_nch`           | `int`              | Number of LMS channels with data |
| `lms_id`            | `uint8[lms_nch]`   | 0=Pin, 1..3 = LMS1..3 |
| `lms_npeaks`        | `int[lms_nch]`     | Soft peaks found |
| `lms_peak_time`     | `float[lms_nch][8]` | Peak time (ns) |
| `lms_peak_integral` | `float[lms_nch][8]` | Peak integral |

# `scalers` tree (DSC2 livetime)

Written by every replay tool.  One row per DSC2 SYNC readout (~once per
SYNC interval, typically 1–2 s).  Counts accumulate from the GO
transition and are not reset between rows; instantaneous quantities
require the difference between two consecutive entries.  See
[`docs/rols/banktags.md`](rols/banktags.md) §0xE115 and
`prad2dec/include/Dsc2Decoder.h` for the bank format and the in-DAQ
gating convention (Group A counts during live; Group B is free-running).

Join key: a scaler row is emitted inside a particular physics event;
its `event_number` matches that physics event's `event_num` in the
main tree.

| Branch | Type | Meaning |
|---|---|---|
| `event_number`  | `int`     | Physics event number this DSC2 read lives inside |
| `ti_ticks`      | `int64`   | TI 48-bit timestamp at this read (250 MHz ticks) |
| `unix_time`     | `uint32`  | Most recent EPICS unix_time observed at this read (0 before any EPICS arrived) |
| `sync_counter`  | `uint32`  | Most recent EPICS HEAD-bank counter |
| `run_number`    | `uint32`  | CODA run number |
| `trigger_type`  | `uint8`   | Trigger type of the carrying physics event |
| `slot`          | `int`     | DSC2 slot in its VME crate |
| `gated`         | `uint32`  | Selected source counter, group A (live) |
| `ungated`       | `uint32`  | Selected source counter, group B (total) |
| `live_ratio`    | `float`   | `gated / ungated` — cumulative live fraction at this read; -1 if `ungated == 0` |
| `source`        | `uint8`   | Selection: `0` = ref, `1` = trg, `2` = tdc (matches `daq_config.json:dsc_scaler.source`) |
| `channel`       | `uint8`   | Selected channel index 0–15; ignored when `source == ref` |
| `ref_gated`     | `uint32`  | DSC2 reference counter, group A |
| `ref_ungated`   | `uint32`  | DSC2 reference counter, group B |
| `trg_gated`     | `uint32[16]` | Per-channel TRG counter, group A |
| `trg_ungated`   | `uint32[16]` | Per-channel TRG counter, group B |
| `tdc_gated`     | `uint32[16]` | Per-channel TDC counter, group A |
| `tdc_ungated`   | `uint32[16]` | Per-channel TDC counter, group B |
| `good`          | `bool`    | **Only in `replay_filter` output** — overall slow-control verdict at this checkpoint (all configured cuts passed) |

# `epics` tree (slow control)

Written by every replay tool.  One row per EPICS event (top-level EVIO
tag `0x001F`); each row carries the channel/value pairs from a single
0xE114 string bank parsed via `epics::ParseEpicsText`.  Channel names
are heterogeneous between rows: only those that updated in this EPICS
event are listed.  Persistent run-wide channel registry and snapshot
indexing live in `epics::EpicsStore` (monitor-server side, not in the
replay tree).

Join key: `event_number_at_arrival` is the most recent physics
`event_num` observed by the decoder at the time this EPICS event
arrived (`-1` for EPICS that arrived before any physics event).
`ti_ticks_at_arrival` is the TI 48-bit tick of that same physics
event, captured at decode time so analysis does not have to look it
back up via the events tree (which may have dropped the anchor event
during replay-time filtering).

| Branch | Type | Meaning |
|---|---|---|
| `event_number_at_arrival` | `int`               | Most recent physics `event_num` at EPICS arrival; `-1` if none |
| `ti_ticks_at_arrival`     | `long long`         | TI tick of that physics event (4 ns / tick); `0` if none |
| `unix_time`               | `uint32`            | Absolute Unix seconds (from the 0xE112 HEAD bank) |
| `sync_counter`            | `uint32`            | Monotonic HEAD-bank counter |
| `run_number`              | `uint32`            | CODA run number |
| `channel`                 | `vector<string>`    | Channel names that updated in this EPICS event |
| `value`                   | `vector<double>`    | Parallel values; `value[i]` is `channel[i]`'s reading |
| `good`                    | `bool`              | **Only in `replay_filter` output** — overall slow-control verdict at this checkpoint |

ROOT vector branches require a stable pointer-to-pointer for
`SetBranchAddress`; see `prad2det/include/EventData_io.h` ::
`SetEpicsReadBranches` for the canonical reader skeleton.

# `runinfo` tree (CODA control events / DAQ config)

Written by every replay tool.  One row per CODA control event encountered
in the input EVIO stream — typically PRESTART, GO, and END once each per
run.  `replay_filter` passes the input runinfo rows through verbatim
(concatenated across input files; no filter applied — the whole point
of this tree is run-scoped metadata).

The PRESTART row is the one that carries the long DAQ-config text from
the `0xE10E` STRING bank: TI / DSC2 / FADC250 / FAV3 / TDC1190 / SSP /
VTP / TS settings, per-channel pedestals + gains, trigger masks — i.e.
the full set of `.cnf` files CODA concatenated to start the run.  GO
and END rows have `daq_config` empty (the bank is PRESTART-only) but
preserve the absolute Unix time for run start/end accounting.

Distinguish rows by `event_tag`: `0x11` = PRESTART, `0x12` = GO,
`0x14` = END.  GO/END events do not always carry a populated
`run_number` slot in the underlying CODA control bank, so `run_number`
on those rows may be `0`; PRESTART is authoritative.

| Branch | Type | Meaning |
|---|---|---|
| `run_number` | `uint32`  | CODA run number (PRESTART authoritative; may be 0 on GO/END) |
| `unix_time`  | `uint32`  | Absolute Unix seconds of the control event |
| `run_type`   | `uint8`   | CODA run type (non-zero on PRESTART; 0 on GO/END) |
| `event_tag`  | `uint8`   | `0x11` PRESTART / `0x12` GO / `0x14` END |
| `daq_config` | `string`  | Full DAQ configuration text from the 0xE10E bank — empty on GO/END |

The `daq_config` payload is large (~70 KB raw text on PRad-II runs)
but compresses to ~25 KB; storage cost is one-time per run.  Recipe
for extracting the text from a replayed file:

```cpp
TTree *t = (TTree*)f->Get("runinfo");
prad2::RawRunInfo ri;
prad2::SetRunInfoReadBranches(t, ri);
std::string *sp = &ri.daq_config;
t->SetBranchAddress("daq_config", &sp);
for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    ri.daq_config.clear();
    t->GetEntry(i);
    if (ri.event_tag == 0x11) {
        // ri.daq_config is the PRESTART config text
    }
}
```

The DAQ config text is the same `.cnf`-style format CODA writes to
disk; grep it for `VTP_CRATE`, `SSP_PRAD_TRGBIT`, `FAV3_TET`, etc. to
read out the run's trigger / readout settings.

# Run example

```bash
ssh clasrun@clonfarm11
source ~/prad2_daq/prad2_env.csh
cd /data/replay_raw/
prad2ana_replay_rawdata /data/evio/data/prad_024154/prad_024154.evio.00000 -o ./ -p
```

Output: `/data/replay_raw/prad_024154.00000_raw.root`.

DAQ-emulation knobs (`TET` / `NSB` / `NSA` / `NPEAK` / `NSAT` / `NPED` /
`MAXPED`) are read from the `fadc250_waveform.firmware` block in
`daq_config.json`; the offline soft analyzer (`WaveAnalyzer`) is configured
from the sibling `fadc250_waveform.analyzer` block — override either to
match the actual run.
