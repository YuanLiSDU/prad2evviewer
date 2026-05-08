# Adding HV data to a replay file & filtering events by HV stability

This is a hands-on walkthrough of the two-step workflow that joins the
PRad-II HV archive (the daily `vmon_YYYYMMDD.dat` files written by
`prad2hvd`) onto a replayed run, then keeps only the events whose HV was
stable on a chosen set of channels.

```
   vmon_*.dat archive ──┐
                        ├──► prad2ana_replay_add_hv ──► hv* trees inside
   prad_NNNNNN_*.root ──┘                              the same replay file
                                                          │
                                                          ▼
                                                hv_event_filter.py ──► stable-event list,
                                                                       JSON report
```

Step 1 (`prad2ana_replay_add_hv`) is a compiled C++ analysis tool that
appends five trees to an existing replay ROOT file in place — no second
ROOT file is produced.  Step 2 (`hv_event_filter.py`) is a Python script
that reads those trees, computes per-channel rolling-std on dV, and tags
each physics event with a stability flag.

The workflow is idempotent on the replay file (use `-f` to re-run with
different channel filters) and the reports are reproducible.

---

## Prerequisites

* The replay file must already contain a `scalers` tree with the standard
  branches (`event_number`, `ti_ticks`, `unix_time`).  All
  `prad2ana_replay_*` outputs satisfy this — `_recon.root`,
  `_filtered.root`, raw `_raw.root` all work.
* `vmon_YYYYMMDD.dat` archive files for the run's date(s).  The default
  search location is `/data/prad2/hv_data` (override with `-d`).
* `prad2evviewer` build with the analysis target enabled
  (`-DBUILD_ANALYSIS=ON`, the default).  After `make install` the binary
  lands at `<prefix>/bin/prad2ana_replay_add_hv`.
* For step 2: Python 3 with `numpy` and `uproot` installed
  (`pip install uproot`).

---

## Step 1 — append HV trees to the replay file

### CLI

```
prad2ana_replay_add_hv [-d <hv_dir>] [-c "ch1,ch2,..."] [-p <pad_s>]
                       [-f] file1.root [file2.root ...]
```

| Option | Default | Meaning |
|---|---|---|
| `-d <dir>`  | `/data/prad2/hv_data` | Directory holding `vmon_*.dat`. |
| `-c <list>` | _all_                 | Comma-separated HV channel names to keep (e.g. `"W1,W100,G1,G2"`).  Filtering at write time keeps the file small when you only care about a few channels. |
| `-p <sec>`  | `30`                  | Pad around the scalers `unix_time` window in seconds. |
| `-f`        | _off_                 | Overwrite existing `hv*` trees in the file (deletes the old keys first; the file is then re-populated with fresh cycles). |

### What it writes

Five trees, all in the same ROOT file as the input:

| Tree | Entries | Purpose |
|---|---|---|
| `hv` | one per HV snapshot (200 ms cadence) | per-channel `dv[N_ch]/F` and `v0set[N_ch]/F`, plus event-key columns |
| `hv_channels` | one per HV channel | `channel_id/s`, `name/C` — name lookup table for the columns of `hv.dv` / `hv.v0set` |
| `hv_booster` | one per booster snapshot (only when `n_boosters > 0`) | TDK-Lambda booster monitoring: `vmon/imon/vset/iset[N_bst]/F` |
| `hv_booster_channels` | one per booster | name lookup for the booster columns |
| `hv_meta` | exactly one entry | `run_t_start_unix/D`, `run_t_end_unix/D`, `interval_ms/I`, `n_channels/I`, `n_boosters/I` |

### `hv` branch reference

| Branch | Type | Meaning |
|---|---|---|
| `event_number_at_arrival` | `int32`/I | Event whose TI tick is the largest `≤` this snapshot's converted ticks.  `-1` for snapshots before the first physics event. |
| `ti_ticks_at_arrival`     | `int64`/L | That event's `recon.timestamp` (4-ns ticks). |
| `unix_time`               | `uint32`/i | Snapshot wall-clock in epoch seconds (matches `scalers.unix_time` / `epics.unix_time`). |
| `t_unix_s`                | `double`/D | Same instant in fractional seconds — keep this if you need ms precision (the recorder's poll cadence is 200 ms). |
| `dv[N_ch]`                | `float32`/F | `VMon - V0Set` for each kept channel. |
| `v0set[N_ch]`             | `float32`/F | Active V0Set for each channel.  Reconstructed `VMon = dv + v0set`. |

The `event_number_at_arrival` / `ti_ticks_at_arrival` follow the same
contract as the existing `epics` tree, so the join idiom is the same:
look up the most recent event where `recon.event_num ≤
event_number_at_arrival`.

### Example: tag run 24340

```bash
# Input: prad_024340_filtered.root (already contains recon/scalers/epics)
prad2ana_replay_add_hv \
    -d /data/prad2/hv_data \
    -c "W1,W100,G1,G2,PRIMARY1_0" \
    prad_024340_filtered.root
```

Expected stdout:

```
[replay_add_hv] 1 file(s); HV archive: /data/prad2/hv_data

[file] prad_024340_filtered.root
  run window  : 2026-05-05T01:33:21 → 2026-05-05T01:41:00  (459 s, 72 SYNC pins)
  +pad (30s) : 2026-05-05T01:32:51 → 2026-05-05T01:41:30
  event tree  : recon (836741 entries)
  anchor      : ti_ticks=5509953108  unix_ms=1777944801000
  → wrote HV trees to prad_024340_filtered.root
     hv                  : 2150 entries × dv[5] / v0set[5] (float32)
     hv_channels         : 5 entries
     hv_booster          : 8 entries × vmon/imon/vset/iset[3]
     hv_booster_channels : 3 entries (pradbst4, pradbst5, pradbst6)
     hv_meta             : 1 entry (interval=200 ms)
     run window  : 2026-05-05T01:33:21 → 2026-05-05T01:41:00 (459 s)
     event range : -1 → 5299427 (from 836741 events in `recon`, anchored to 72 SYNC pins)
     HV archive  : vmon_20260504.dat

[replay_add_hv] done — 1 ok, 0 failed
```

The `event_range -1 → 5299427` line confirms the anchor-based lookup
worked: the `-1` rows are HV snapshots that fall in the pre-run pad
(before the first physics event), and `5299427` is the event_num of the
last physics event tagged.

### Example: full-archive snapshot (every channel)

Drop `-c` to keep every channel in the archive (1792 channels for a
PRad-II run).  Output file grows by ~1 MB per minute of run; for an
8-minute run with 1792 channels you get a `hv` tree of about 30 MB.

```bash
prad2ana_replay_add_hv -d /data/prad2/hv_data prad_024340_filtered.root
```

Use this when you don't yet know which channels you'll cut on — the
filter step (below) can then scope the analysis to any subset, type, or
single channel without re-running step 1.

### Example: re-tagging with `-f`

```bash
prad2ana_replay_add_hv -f -c "W1,W100" prad_024340_filtered.root
```

This deletes the existing `hv*` keys and writes fresh cycles with the
narrower channel set.  The freed bytes stay in the file until you repack
it (e.g. `rootcp -O <input> <output>`), so an `-f` cycle does not shrink
the file on disk.

### Sanity-check the result

```bash
root -l -q -e '
TFile *f = TFile::Open("prad_024340_filtered.root");
auto *t = (TTree*)f->Get("hv");
t->Print();
t->Draw("dv[1]:t_unix_s-1777944801", "", "P");
'
```

You should see a sub-volt fluctuation band around 0 V if HV was stable.

---

## Step 2 — filter events by HV stability

The Python script lives at
`<repo>/analysis/pyscripts/hv_event_filter.py` (also installed under
`<prefix>/share/prad2evviewer/analysis/pyscripts/`).  It reads the `hv`
trees written in step 1 and emits:

1. a stdout summary,
2. a list of stable physics-event numbers (optional),
3. a JSON report with per-channel rejection counts (optional).

### CLI

```
python3 hv_event_filter.py <replay.root> \
    [--scope all|type|channel] [--type-name PbWO4] [--channel-name W100] \
    [--abs-threshold V] [--rel-n-mad N] [--require-v0set] \
    [--window-s SECS] [--top-n N] \
    [--out-events FILE] [--out-json FILE]
```

| Option | Default | Meaning |
|---|---|---|
| `--scope`         | `all`     | `all` → every channel found in `hycal_map.json`. `type` → restrict to one of `PbWO4` / `PbGlass` / `Veto` / `LMS`.  `channel` → a single channel by name. |
| `--type-name`     | _required for `--scope type`_     | Module type. |
| `--channel-name`  | _required for `--scope channel`_  | Channel name (e.g. `W100`). |
| `--abs-threshold` | _disabled_ | Absolute fluctuation gate: rolling-std(dV) above this many volts marks the channel unstable. |
| `--rel-n-mad`     | _disabled_ | Relative fluctuation gate: rolling-std(dV) above `median + N · 1.4826 · MAD` (over the run, computed from the rolling-std distribution) marks it unstable.  Robust against a few transient spikes. |
| `--require-v0set` | _off_     | Also require `|dv| ≤ --abs-threshold` (channel is on setpoint).  Off by default — channels held intentionally low don't fail. |
| `--window-s`      | `5`       | Rolling-std window in seconds. |
| `--top-n`         | `20`      | How many worst-rejecting channels to print. |
| `--out-events`    | _none_    | Stable `event_num` list, one per line. |
| `--out-json`      | _none_    | Full report (parameters, totals, per-channel histogram). |
| `--hycal-map`     | `database/hycal_map.json` | The map used to resolve channel name → module type. |
| `--event-tree`    | auto      | Tree name carrying physics events.  Tries `recon` then `events`. |

You must enable at least one of `--abs-threshold` / `--rel-n-mad`.

### How the gates work

For each HV snapshot and each monitored channel:

1. Compute a rolling population-std of dV over a `--window-s` window
   (NaN-aware so an offline channel doesn't poison its neighbours).
2. The absolute gate trips if `rolling_std > --abs-threshold`.
3. The relative gate trips if
   `rolling_std > median(rolling_std) + --rel-n-mad · 1.4826 · MAD(rolling_std)`.
   The `1.4826` factor turns MAD into a Gaussian-σ-equivalent, so
   `--rel-n-mad 5` means "more than 5σ above the typical fluctuation
   level for this channel".
4. With `--require-v0set`, also trip if `|dV| > --abs-threshold` (i.e.
   channel is off setpoint).

A snapshot is **stable** when no monitored channel trips.

The flag is propagated to physics events through
`event_number_at_arrival`: each event picks up the stability of the most
recent snapshot whose `event_number_at_arrival ≤ event.event_num`.
Events that precede the first HV snapshot are counted as unstable and
reported separately in `n_events_pre_hv`.

### Example A — single channel (`W100`), absolute + relative

```bash
python3 hv_event_filter.py prad_024340_filtered.root \
    --scope channel --channel-name W100 \
    --abs-threshold 0.5 --rel-n-mad 5 \
    --window-s 5 \
    --out-events stable_events_W100.txt
```

Expected output (run 24340):

```
[hv_event_filter] prad_024340_filtered.root
  hv tree         : 2,150 snapshots × 5 channels (period=200 ms)
  event tree      : recon, 836,741 entries
  hycal map       : .../database/hycal_map.json (1735 entries)
  scope           : channel (channel=W100)
  monitored ch    : 1 channel(s)
  rolling window  : 5.0 s (25 samples)
  absolute gate   : rolling-std ≤ 0.5 V
  relative gate   : rolling-std ≤ median + 5 · 1.4826 · MAD

  HV snapshots    : 2,150 total, 2,076 stable (96.56%), 74 unstable
  events          : 836,741 total, 823,373 stable (98.40%)

  Top 20 rejecting channel(s):
    Channel       type        n_unstable     frac
    ----------------------------------------------
    W100          PbWO4               74    3.44%

  → stable event list → stable_events_W100.txt (823,373 events)
```

### Example B — by module type (every PbWO4 channel)

This requires step 1 to have been run **without** `-c` (so all 1792
channels are present).

```bash
python3 hv_event_filter.py prad_024340_filtered.root \
    --scope type --type-name PbWO4 \
    --abs-threshold 0.5 --rel-n-mad 5 \
    --top-n 8
```

The "top rejecting channels" report tells you which modules are
chronically noisy:

```
  Top 8 rejecting channel(s):
    Channel       type        n_unstable     frac
    ----------------------------------------------
    W25           PbWO4             2150  100.00%
    W887          PbWO4             2111   98.19%
    W888          PbWO4             2073   96.42%
    W890          PbWO4             2066   96.09%
    W922          PbWO4             2046   95.16%
    W820          PbWO4             2041   94.93%
    W683          PbWO4             2028   94.33%
    W821          PbWO4             2021   94.00%
```

`W25` triggering on every snapshot is a smoking gun for a chronically
high-noise channel — masking it before tightening the cuts is a typical
next step.

### Example C — every HyCal channel, robust (relative-only)

When you want the cut to learn what's normal for each channel rather
than imposing a global volt threshold:

```bash
python3 hv_event_filter.py prad_024340_filtered.root \
    --rel-n-mad 5 \
    --out-json hv_filter.json
```

The relative gate uses each channel's own median + MAD as the baseline,
so a quiet channel is held to a tight threshold and a noisier channel
isn't penalised for its baseline noise level.  Combine with
`--abs-threshold` if you also want a hard ceiling.

### Example D — feed the stable-event list into a downstream cut

```bash
python3 hv_event_filter.py prad_024340_filtered.root \
    --scope type --type-name PbWO4 \
    --abs-threshold 1.0 \
    --out-events stable_24340.txt

# In analysis ROOT macro:
#   std::set<int> good = read_event_list("stable_24340.txt");
#   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
#       t->GetEntry(i);
#       if (!good.count(event_num)) continue;
#       ...
#   }
```

`stable_24340.txt` is a plain newline-separated list of `event_num`
values that passed the gate; load it into a `std::set<int>` (or
`numpy.fromfile(..., dtype=int, sep="\n")` in Python) and use it as a
membership test.

### Inspecting the JSON report

`--out-json` writes a full report:

```json
{
  "input": "prad_024340_filtered.root",
  "scope": "type",
  "type_name": "PbWO4",
  "abs_threshold": 0.5,
  "rel_n_mad": 5.0,
  "window_s": 5.0,
  "interval_ms": 200,
  "n_snapshots": 2150,
  "n_snapshots_stable": 0,
  "n_snapshots_unstable": 2150,
  "n_events": 836741,
  "n_events_stable": 0,
  "rejection_per_channel": { "W25": 2150, "W887": 2111, ... },
  "per_channel_stats": {
    "W25": { "type": "PbWO4", "n_unstable_snapshots": 2150,
              "rs_median": 0.71, "rs_mad": 0.06, "rs_thresh_rel": 1.155 },
    ...
  }
}
```

`per_channel_stats[*].rs_median` and `rs_mad` summarise each channel's
rolling-std distribution; `rs_thresh_rel` is the actual relative cut
applied (`median + N·1.4826·MAD`).  Use these to pick a sensible
`--abs-threshold` after seeing the per-channel scatter.

---

## End-to-end test

A repeatable test against the in-tree `prad_024340_filtered.root`:

```bash
# (1) Make a working copy so we don't mutate the canonical file.
cp /home/cpeng/Projects/prad/prad2evviewer/prad_024340_filtered.root \
   /tmp/test_run24340.root

# (2) Tag with HV — keep all channels.
prad2ana_replay_add_hv \
    -d /mnt/hgfs/Data/PRad2/hv_data \
    /tmp/test_run24340.root

# Expect: 836741 events tagged, anchor=ti_ticks=5509953108 unix_ms=1777944801000

# (3) Quick stability filter on the W100 channel.
python3 analysis/pyscripts/hv_event_filter.py \
    /tmp/test_run24340.root \
    --scope channel --channel-name W100 \
    --abs-threshold 0.5 --rel-n-mad 5 \
    --out-events /tmp/test_24340_stable.txt

# Expect: 823,373 / 836,741 events stable (98.40%); W100 contributes 74 unstable snapshots.

# (4) Cross-check via PyROOT — pick a stable event and confirm.
python3 -c "
import uproot, numpy as np
f = uproot.open('/tmp/test_run24340.root')
hv = f['hv']
ts = hv['t_unix_s'].array(library='np')
ev = hv['event_number_at_arrival'].array(library='np')
dv = hv['dv'].array(library='np')
print('hv n_snap =', len(ts))
ci = list(map(str, f['hv_channels']['name'].array(library='np'))).index('W100')
print('W100 dv first 5:', dv[:5, ci])
"
```

Replace step (1)'s source path with your own replay file.

---

## What the trees look like in `recon` queries

Once `hv*` is in the file, it composes naturally with other replay
trees.  A few common idioms:

```c++
// Avg HyCal cluster energy on stable HV events only
TTree *recon = (TTree*)f->Get("recon");
TTree *hv    = (TTree*)f->Get("hv");
recon->AddFriend(...);   // hv is keyed by event_number_at_arrival, not event_num,
                         // so a strict TTree::AddFriend won't work directly.
                         // Instead, build a (event_num → stable) map up front:
//   1) Stream `hv` once, mark each entry's stability flag.
//   2) Run upper_bound on `event_number_at_arrival` for each recon row.
// (See the Python implementation in hv_event_filter.py for the canonical version.)

// Plot dV trace for one channel
hv->Draw("dv[ci]:t_unix_s", Form("event_number_at_arrival>=%d", run_start_evn));
```

For ROOT-side analyses, you can also run `hv_event_filter.py
--out-events` and then filter `recon` rows by membership in that list —
it's the simpler workflow when the cut definition lives in Python.
