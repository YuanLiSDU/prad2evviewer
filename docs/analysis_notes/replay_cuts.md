# PRad-II replay + slow-control filter — quick start

End-to-end walkthrough for one run on `clasrun@clonfarm11`: stage the
EVIO from the DAQ disk, replay it (raw or recon), apply slow-control
cuts (livetime + EPICS) and inspect the per-checkpoint pass/fail trace.
The example uses run **024327**; substitute the run number you actually
want to process.

```
   clondaq2:/data/stage2 ──► clonfarm11:/data/evio ──► prad2ana_replay_rawdata  ─┐
                                                  └─► prad2ana_replay_recon ────┤
                                                                                ▼
                                                                prad2ana_replay_filter
                                                                                │
                                                                                ▼
                                      *_filter.root + prad_<run>_epics.root + *.report.json
                                                                                │
                                                                                ▼
                                                          prad2ana_replay_report_viewer
```

The four `prad2ana_replay_*` binaries are installed by the standard
`make install` (build with `-DBUILD_ANALYSIS=ON`, the default) and land
under `<prefix>/bin/`. The `prad2_env.csh` shell snippet adds them to
`PATH` and sets `PRAD2_*` data-dir env vars — sourcing it is enough,
no per-tool flags are needed for the canonical layout on
`clonfarm11`.

---

## 1. Stage the raw EVIO from clondaq2

```
ssh clasrun@clonfarm11
scp -r clondaq2:/data/stage2/prad_024327 /data/evio/
```

A run directory contains all `prad_NNNNNN.evio.NNNNN` split files. The
replay tools accept either a directory or a list of files, so copying
the whole directory is the simplest workflow.

---

## 2. Set up the environment

```
source /home/clasrun/prad2_daq/prad2_env.csh
```

This puts the build's `bin/` first on `PATH` and points the analysis
tools at the installed `database/` (DAQ config, HyCal map, GEM map).

---

## 3. Replay (raw or recon)

Each split file is replayed independently and produces one ROOT file
with three trees in common — `events` (raw) or `recon`
(reconstructed), `scalers`, and `epics`. The `scalers` and `epics`
trees are what the next step's cuts read; the events/recon tree is
what they filter.

### Raw (no clustering / matching)

```
mkdir -p /data/replay_raw/prad_024327
prad2ana_replay_rawdata /data/evio/prad_024327 \
    -o /data/replay_raw/prad_024327 -p
```

`-p` adds the FADC peak-analysis branches (soft + firmware DAQ-mode
peaks) — handy for low-level checks; drop it for the leanest output.
Default thread count is 4; pass `-j N` to change it.

### Reconstructed (HyCal clusters + GEM hits + HyCal↔GEM match)

```
mkdir -p /data/replay_recon/prad_024327
prad2ana_replay_recon /data/evio/prad_024327 \
    -o /data/replay_recon/prad_024327 \
    -j 8
```

The recon output is what physics analyses normally consume; the raw
output is mostly for diagnostics or if you need access to the FADC
samples that the recon step has already collapsed into clusters. Pick
one — there is no need to run both. By default replay recon also merges
successful split outputs in groups of 80 with `hadd`; add `-m 0` if a
cut workflow should keep only the per-split `_recon.root` files.

---

## 4. Apply slow-control cuts (livetime / EPICS)

Pass every replayed ROOT file for the run in **one** invocation so
livetime deltas and EPICS forward-fill are computed across the full
run, not per split file. The cut JSON template ships at
`analysis/cuts/prad2_default.json`; copy it and tune the thresholds
for your analysis.

```
cd /data/replay_recon/prad_024327
prad2ana_replay_filter prad_024327_recon_*.root \
    -o . \
    -c /home/clasrun/prad2_daq/prad2evviewer/analysis/cuts/prad2_default.json \
    -t 8
```

Output:

* `prad_024327_recon_000_filter.root`, `prad_024327_recon_001_filter.root`,
  ... — one filtered ROOT per input ROOT, named by inserting `_filter`
  before the final `.root`; each file keeps its own `scalers`, `epics`,
  and `runinfo` trees with an extra `good` bool on the slow trees.
* `prad_024327_epics.root` — run-level `scalers`, `epics`, and `runinfo`
  only, concatenated from every input file.
* `prad_024327_filter_report.json` — per-(channel, checkpoint) pass/fail
  trace, robust median + MAD per channel, `keep_intervals`, and (when
  `charge` is configured) live-charge integration. Suitable for plotting
  in a quality dashboard.

Then compute post-cut live charge from all filtered ROOTs together:

```
prad2ana_live_charge prad_024327*_filter*.root \
    -j prad_024327_live_charge.json
```

The convenience scripts in `scripts/shell/` run this sequence as
`replay_recon -> replay_filter -> live_charge -> quick_check`, write all
products directly under `<output_base>/prad_<run>/`, and reuse the same
CPU count for replay, filter, and quick check.

### Cut JSON in 30 seconds

```jsonc
{
  "livetime": {
    "source":  "trg",        // "ref" | "trg" | "tdc"
    "channel": 2,            // ignored for "ref"
    "abs":     { "min": 90, "max": 100 }
  },
  "epics": {
    "hallb_IPM2C21A_CUR":  { "abs": { "min": 0.3 } },
    "hallb_IPM2C21A_XPOS": { "rel_rms": 3, "gated_by": "hallb_IPM2C21A_CUR" },
    "hallb_IPM2C21A_YPOS": { "rel_rms": 3, "gated_by": "hallb_IPM2C21A_CUR" }
  },
  "charge": { "beam_current": "hallb_IPM2C21A_CUR" }
}
```

* `abs` is a hard min/max gate (drop the side you don't care about).
* `rel_rms: N` keeps points within N · σ̂ of the channel's median,
  where σ̂ = 1.4826 · MAD — robust to a few transient spikes.
* `gated_by` restricts the rel_rms median/MAD calculation to points
  where the named channel's cut already passed. Useful when a channel
  is only meaningful in a regime — e.g. beam position is meaningless
  while the current is below the floor. The cut itself is still
  applied to every checkpoint; only the stats are conditioned.
* `charge` is opt-in. When present, the report adds a `live_charge`
  block with both gated and ungated `value_nC` (assuming the named
  channel reports nA, the Hall B IPM convention) plus the live and
  real seconds used.

> **DSC2 cabling note.** PRad-II runs sampled so far (024246, 024340)
> have the Faraday signal on `trg`/`tdc` channel 2 (channels 0 and 1
> are unwired). Switch the channel here if a future run moves the
> cable.

---

## 5. Inspect the cut effects

### GUI

```
prad2ana_replay_report_viewer
```

Click **Open Report** and pick the report JSON. Three linked rows share
the x-axis (associated_timestamp by default, switchable to
associated_evn): cut status across channels, livetime + data rate, and
EPICS values with the accept band drawn in. Pan/zoom via the matplotlib
toolbar; **Save** dumps the current view as PNG/PDF.

### Headless

```
prad2ana_replay_report_viewer --cli prad_024327_filter_report.json
```

Renders one static figure with a row per channel — status on top, then
livetime + rate, then one row per EPICS channel. Figure height scales
with the channel count. Add `--evn` to plot against `associated_evn`
instead of seconds.

---

## Caveat — old replay outputs

ROOT files replayed **before** the `scalers` and `epics` trees were
added will not work with `prad2ana_replay_filter`. The tool builds its
keep-interval timeline from rows of those two trees; with the trees
absent the timeline is empty, every physics event falls outside any
keep-interval, and the output `events` tree comes out empty. **No
error is raised** — the report's `summary.n_physics_pass` will simply
be `0`.

Re-run step 3 (`prad2ana_replay_rawdata` / `prad2ana_replay_recon`)
against the original EVIO to regenerate the side trees before
filtering. A quick sanity check on a candidate file:

```bash
root -l -q -e '
TFile *f = TFile::Open("prad_024327.00000_recon.root");
std::cout << "scalers: " << (f->Get("scalers") ? "yes" : "no") << "\n";
std::cout << "epics  : " << (f->Get("epics")   ? "yes" : "no") << "\n";
'
```

If either prints `no`, the file predates the side-tree change and
must be re-replayed.
