# Online LMS Gain Monitor

This document describes the online HyCal PbWO4 gain-monitoring workflow
implemented by:

- `scripts/online_gain_monitor.py` — PyQt6 monitor, downloader, replay
  queue manager, and visualisation front end;
- `analysis/tools/online_gain_monitor_base.cpp` — incremental EVIO → LMS
  replay and gain-correction back end;
- `analysis/src/GainCorrCompute.cpp` — shared batch calculation used to
  produce the `gain_corr` ROOT tree.

The monitor is intended to run while DAQ data are still arriving on
`clondaq2`. It copies only completed EVIO splits, reuses persistent
intermediate LMS files, and refreshes module gain maps and time-series
plots after each successful update.

```
 clondaq2:/data/stage2/prad_RUN
                 │
                 │ ssh directory scan + scp
                 ▼
 <storage>/evio/prad_RUN
                 │
                 │ move completed download set into a durable queue
                 ▼
 <storage>/replay_queue/prad_RUN/<snapshot>
                 │
                 │ prad2ana_online_gain_monitor_base
                 ▼
 <storage>/lms/prad_RUN/*_lms.root
                 │
                 │ recompute from every available LMS file
                 ▼
 <storage>/gain/prad_RUN/prad_RUN_gain_corr.root
                 │
                 └──► online_gain_monitor.py
                       HyCal map + selected-module chart + ref-PMT chart
```

---

## 1. Prerequisites

Build the project with analysis and Python tools enabled:

```bash
cmake -S . -B build -DBUILD_ANALYSIS=ON -DBUILD_PYTHON=ON
cmake --build build -j$(nproc)
```

Required runtime components:

| Component | Purpose |
|---|---|
| ROOT | C++ replay and gain fitting |
| Python 3 + PyQt6 | GUI |
| NumPy + uproot | Read `gain_corr.root` without blocking the GUI |
| `ssh` and `scp` | Scan and copy EVIO files from the DAQ host |
| `hycal_map.json` | HyCal geometry and module names |
| reference gain `.dat` file | Reference values used to form gain corrections |

Install the Python dependencies with:

```bash
python3 -m pip install -r scripts/requirements.txt uproot
```

The GUI honours `PRAD2_DATABASE_DIR`. If it is unset, it uses the
repository's `database/` directory.

Linux scheduling tools are optional:

- `taskset` separates the GUI/download work from replay CPUs;
- `nice` lowers replay CPU priority;
- `ionice` lowers replay disk priority.

The monitor continues without CPU affinity when these facilities are
unavailable.

The default online locations are:

```text
remote host : clondaq2
remote base : /data/stage2
storage base: /data/gain_monitor
```

The account running the monitor must be able to use non-interactive
`ssh`/`scp` to the selected host.

---

## 2. Starting the monitor

From an installed build:

```bash
online_gain_monitor 24929
```

From the source tree:

```bash
python3 scripts/online_gain_monitor.py 24929
```

Multiple initial runs may be supplied:

```bash
online_gain_monitor 24929 24930 24931
```

Select a colour theme with:

```bash
online_gain_monitor --theme dark 24929
```

The available values for `--theme` come from
`scripts/hycal_geoview.py`.

Supplying a run number loads an existing gain ROOT file when one is
already present, but it does not start polling automatically. Press
**Start** to begin online scanning.

---

## 3. Online workflow

### 3.1 Remote scan

For each scan, the GUI:

1. lists `prad_NNNNNN` directories below the remote base;
2. starts with the requested run and includes later runs discovered on
   the remote host;
3. lists the run's `.evio.` split files in version order;
4. deliberately excludes the newest remote split;
5. skips files that already have a persistent `_lms.root`, are already
   queued, or are already downloaded;
6. copies up to the configured limit with `scp`.

The newest split is held back because CODA may still be writing it. It
becomes eligible after a newer split appears.

When a later run is discovered and queued, the scan floor advances to
that run. Earlier runs are no longer polled automatically.

### 3.2 Durable replay queue

Downloaded files first land in:

```text
<storage>/evio/prad_RUN/
```

At the end of a scan they are moved into a timestamped snapshot:

```text
<storage>/replay_queue/prad_RUN/YYYYMMDD_HHMMSS_PID/
```

The queue is intentionally disk-backed. If the GUI or machine stops,
pressing **Start** later recovers existing snapshot directories and
continues processing them.

Snapshots are deleted only after the C++ update tool exits
successfully. A failed update leaves the original EVIO files in the
queue for inspection or retry.

All pending snapshots for the same run are combined into one back-end
invocation where possible.

### 3.3 Incremental LMS replay

For every queued EVIO split, the C++ back end calculates its expected
LMS filename with `MakeLMSOutputFile()`:

```text
prad_024929.evio.00012
    └──► prad_024929.00012_lms.root
```

If that file already exists in the work directory, phase 1 is skipped.
Otherwise `Replay::Process_LMSgainFactor()` creates it.

This makes EVIO replay incremental: old splits are not decoded again.
The intermediate files remain under:

```text
<storage>/lms/prad_RUN/
```

### 3.4 Gain-correction update

After phase 1, the back end sorts and chains every `_lms.root` in the
run's work directory. `ComputeGainCorrections()` then recreates the
complete output:

```text
<storage>/gain/prad_RUN/prad_RUN_gain_corr.root
```

Recomputing phase 2 from all persistent LMS files keeps the result
deterministic and prevents a partially updated ROOT tree.

A batch closes after the configured number of LMS events. Alpha events
encountered between those LMS events are included in the same batch.
An incomplete final LMS batch is not written until enough later events
arrive to reach the batch size.

### 3.5 GUI refresh

After a successful back-end update, the GUI loads changed ROOT files in
a worker thread using uproot. Unchanged files are retained from the
in-memory cache. The newest configured runs are then redrawn in:

- the HyCal W-module map;
- the selected module's batch chart;
- the three reference-PMT LMS/Alpha ratio panels;
- the run-average status strip.

---

## 4. Top controls

### Run and process controls

| Control | Meaning |
|---|---|
| **Start run** | One or more numeric run numbers separated by spaces or commas. The first/minimum run seeds online discovery. |
| **Start** | Recover queued snapshots, start replay if needed, and begin periodic remote scans. |
| **Stop** | Stop the timer, kill active download/replay processes, and clear the in-memory queue. Disk snapshots are preserved if replay did not finish. |
| **Scan Now** | Run one remote scan immediately. It may also be used before **Start**, but periodic polling is enabled only while the monitor is running. |
| **Open Folder** | Open a directory containing `prad_NNNNNN_gain_corr.root` files and use it as an offline data source. |
| **every** | Poll interval in seconds when no immediately downloadable file is found or a retry is needed. Range: 5–3600 s; default: 30 s. |
| **batch** | LMS events per gain-fit batch. Range: 100–100000; default: 1000. At the nominal 10 Hz LMS rate, 1000 events are roughly 100 s. |
| **evio skip** | Sampling factor for remote EVIO splits. `0` keeps every eligible split; `N` keeps one split and skips the next `N`. |
| **max evio** | Maximum number of files copied in one scan for each run. Default: 10. |
| **queue cap** | Maximum pending EVIO files across the download and replay-queue trees. Downloads pause at this limit. Default: 100. |
| **runs shown** | Number of latest runs loaded and displayed together. Default: 5. |
| **threads** | Requested C++ replay threads. The GUI may cap this to the CPUs available after reserving resources for the viewer/download work. |

Changing **batch** affects the next back-end recomputation. Existing
LMS intermediates are reused, so changing batch size does not require
downloading or replaying EVIO again.

### Reference and display controls

| Control | Meaning |
|---|---|
| **ref run** | Reference-gain run passed as `-r`. `auto` (`-1`) uses `gain_ref_run` from `runinfo/general.json`. |
| **Ref File** | Select an explicit reference gain `.dat` file. This takes precedence over **ref run** and is passed with `-R`. |
| **quantity** | Select the value shown by the map and selected-module chart; definitions are below. |
| **map range** | Symmetric display range around 1 for gain ratios, or around 0 for change quantities. This changes colours/axes only, not calculations. |
| **drop warn** | Show a warning when the latest short-term gain change of a W module is below the negative threshold. `0%` disables warnings. |
| **ref** | Display Ref1, Ref2, Ref3, or all three plus their valid average. |
| **first 5 → 1** | For non-change quantities, divide each module/reference series by its average over the first five visible batches. |
| **ratio tol** | Reject a reference PMT for a batch when its LMS/Alpha ratio differs from that reference's visible-run mean by more than this percentage. |

Reference-ratio rejection is applied independently to Ref1, Ref2, and
Ref3. A rejected point appears as a red cross in the lower chart and is
excluded from module averages and change calculations.

---

## 5. Quantity definitions

Let

```text
G(run,batch,module,ref) = gain_W / gain_W_ref
```

for a valid module/reference value.

| Quantity | Map meaning | Batch-chart meaning |
|---|---|---|
| **Gain / Ref Gain** | Latest `G` value | `G` for every visible batch |
| **Change** | Latest batch relative to the mean of up to five preceding batches | Each batch relative to its preceding five-batch mean |
| **Long Change** | Mean of the latest five long-change samples | Each batch relative to the first five visible batches |
| **Run-to-Run Change** | Mean change of the latest run relative to the immediately previous run mean | Every batch in a run relative to the previous run mean |
| **Current vs All Runs** | Mean change of the latest run relative to the mean of all earlier visible runs | Every batch in a run relative to the accumulated earlier-run mean |

Change values use:

```text
(current - baseline) / baseline
```

and are displayed as percentages.

The **drop warn** calculation always uses short-term **Change**,
regardless of the quantity currently selected. It compares the latest
batch with up to five preceding batches and averages the three valid
reference channels for each W module.

---

## 6. Main display

### HyCal map

Only the 1156 PbWO4 modules (`W1`–`W1156`) are coloured. Click a module
to select it for the chart on the right.

- ratio quantities use a scale centred on 1;
- change quantities use a blue/green/orange/red scale centred on 0;
- hovering a module shows its current numeric value;
- the selected module receives a highlighted border.

### Selected-module chart

The upper-right chart displays batch values for the selected W module.
Run boundaries are dashed and labelled with six-digit run numbers.

With **Avg** selected, Ref1/2/3 points are shown separately and a line
connects their valid average. The x-axis is a continuous batch index
across visible runs; its label estimates batch duration from
`batch_size / 10 Hz`.

### Reference-PMT ratio chart

The lower-right display has one panel per LMS reference PMT. It plots:

```text
fit_mean_ref_lms / fit_mean_ref_alpha
```

The dotted line is the visible-data mean used by **ratio tol**.
Accepted points are coloured dots; rejected points are red crosses.

### Run-average strip

The strip above the charts reports:

- previous visible run average;
- current run average;
- average across all visible runs.

It follows the currently selected quantity and reference choice.

### Log and status

The lower-left log records scan commands, copied/queued files, replay
output, warnings, and failures. It is capped to avoid unlimited GUI
memory growth.

The status text reports states such as `Idle`, `Scanning...`,
`Replaying...`, or the number of loaded runs and batches.

---

## 7. Path controls and File menu

| Field/action | Meaning |
|---|---|
| **host** | SSH host containing online EVIO files. |
| **remote** | Parent directory containing `prad_NNNNNN` run directories. |
| **storage base** | Local root for downloads, queue snapshots, LMS files, and gain ROOT files. |
| **...** | Browse for a new storage base. |
| **preview work** | Read-only preview of the selected run's persistent LMS directory. |
| **preview out** | Read-only preview of the selected run's gain-correction ROOT file. |
| **File → Open Output Folder...** | Load gain-correction ROOT files recursively from an arbitrary directory without starting online scanning. |
| **File → Reload ROOT** | Re-read changed gain ROOT files from the open folder or monitor storage. |

Opening an output folder is useful for offline comparison of previous
runs. Pressing **Start** returns the application to online monitor
storage.

---

## 8. Back-end command-line tool

The C++ tool can be used without the GUI:

```bash
prad2ana_online_gain_monitor_base \
    /data/gain_monitor/replay_queue/prad_024929/snapshot1 \
    /data/gain_monitor/replay_queue/prad_024929/snapshot2 \
    -w /data/gain_monitor/lms/prad_024929 \
    -o /data/gain_monitor/gain/prad_024929/prad_024929_gain_corr.root \
    -b 1000 -j 16
```

Options:

| Option | Default | Meaning |
|---|---|---|
| `-w <dir>` | required | Persistent directory for per-EVIO `_lms.root` files. |
| `-o <file>` | database gain-correction directory | Output `gain_corr.root`. |
| `-f <N>` | all | Use at most `N` collected EVIO files. |
| `-j <N>` | 4 | Parallel EVIO replay threads. |
| `-c <file>` | `<db>/daq_config.json` | DAQ configuration. |
| `-d <file>` | `<db>/hycal_map.json` | HyCal DAQ map. |
| `-b <N>` | 1000 | LMS events per output batch. |
| `-r <run>` | runinfo reference run | Reference-gain run. |
| `-R <file>` | unset | Explicit reference gain table; overrides `-r`. |

Inputs may be individual EVIO files, directories, or a mixture. The
run number is inferred from the first sorted EVIO path, so one
invocation should contain only one run.

The back end does not delete EVIO inputs. Queue deletion is performed
by the Python monitor only after a successful exit and only for paths
under its expected `replay_queue/prad_RUN/` tree.

---

## 9. `gain_corr` ROOT tree

The output contains one `gain_corr` entry per completed LMS batch:

| Branch | Shape | Meaning |
|---|---|---|
| `batch_id` | scalar | Zero-based batch number |
| `event_num_start` | scalar | First LMS event number in the batch |
| `event_num_end` | scalar | Event number at batch closure |
| `n_lms_events` | scalar | LMS events accumulated |
| `n_alpha_events` | scalar | Alpha events accumulated in the same interval |
| `ref_run` | scalar | Reference-gain run |
| `refPMT_ratio` | `[3]` | Fitted LMS/Alpha ratio for Ref1–3 |
| `gain_W` | `[1156][3]` | Current W-module gain from each reference |
| `gain_W_ref` | `[1156][3]` | Reference gain table copied into the batch |
| `gain_corr_W` | `[1156][3]` | `gain_W_ref / gain_W` correction |
| `fit_mean_ref_lms` | `[3]` | Fitted LMS peak of each reference PMT |
| `fit_mean_ref_alpha` | `[3]` | Fitted Alpha peak of each reference PMT |
| `fit_mean_W_lms` | `[1156]` | Fitted LMS peak of each W module |

See [gain_correction.md](gain_correction.md) for the underlying gain
definitions and fitting details.

---

## 10. Recovery and troubleshooting

### The GUI shows no data

Check that:

1. `numpy` and `uproot` are installed;
2. the expected `prad_NNNNNN_gain_corr.root` exists;
3. it contains a non-empty `gain_corr` tree;
4. at least one complete LMS batch has accumulated.

### Downloads do not start

Test the connection outside the GUI:

```bash
ssh clondaq2 'ls /data/stage2'
```

Also check **queue cap**. Downloads pause while the number of EVIO
files in `evio/` plus `replay_queue/` is at or above that limit.

### The newest EVIO is always missing

This is intentional. The monitor never downloads the newest remote
split because it may still be open. It is picked up after the next
split appears.

### Replay failed

The queued snapshot is retained. Inspect the GUI log, correct the
configuration or reference-gain problem, then press **Start** to
recover and retry the on-disk queue.

### Changing batch size appears to lose the final data

Only complete LMS batches are written. The final partial batch remains
absent until later LMS events bring it to the selected batch size.

### Disk usage grows

Successful queued EVIO snapshots are deleted, but `_lms.root` files are
deliberately persistent because they make future recomputation cheap.
Archive or remove old run directories under `<storage>/lms/` only when
they are no longer needed.

