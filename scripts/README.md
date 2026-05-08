# Scripts

Python utilities for HyCal / tagger / GEM-HyCal monitoring and operator
workflows.

All HyCal GUIs share `hycal_geoview.py` — a reusable HyCal module-map
widget (rectangles in physical coordinates, optional colour bar,
zoom/pan, hover tooltips) and a shared theme palette.  New GUIs
subclass `HyCalMapWidget` and override the paint hooks; the extension
points are documented in `hycal_geoview.py`.

The module exposes `set_theme("dark" | "light")` and an
`apply_theme_palette(window)` helper.  Scripts typically wire a
`--theme` argument in `main()` before constructing any widget.

Installed wrappers are generated when `-DBUILD_PYTHON=ON` is set at
build time — every GUI listed below can then be launched from `$PATH`
by its bare name.  The `python scripts/...` examples shown work
equally well from a source checkout.

## HyCal GUIs

### hycal_event_viewer.py

Event-by-event EVIO browser with two tabs:

- **Waveform** — FADC250 peak finding and per-module histograms.  HyCal geo picker plus raw-waveform plot on the left, four stacked histograms on the right (peak height, integral, time, peaks/event).  Clicking a module in the geo picker switches the selection; **Process next 10 k** fills histograms in a background pass.
- **Cluster** — HyCal energy heatmap with live clustering (`prad2py.det.HyCalCluster`), cluster overlays, and a cluster table.

Opens in **random-access** mode (native event-pointer table via
`evc::EvChannel::OpenRandomAccess`); a single Scan-only pass indexes
physics sub-events so **Prev** / **Next** / **Jump** are fast in
either direction.  **File → Save** writes the current per-module
histograms to JSON for later inspection.

```bash
hycal_event_viewer                        # File → Open…
hycal_event_viewer run.evio.00000
```

### hycal_scaler_map.py

Live colour-coded HyCal FADC scaler map.  Polls
`B_DET_HYCAL_FADC_<name>` EPICS channels every few seconds.
Requires `pyepics` for real EPICS; `--sim` runs without it.

```bash
hycal_scaler_map          # real EPICS (default)
hycal_scaler_map --sim    # simulation (random values)
```

### hycal_pedestal_monitor.py

Measures and monitors FADC250 pedestals on all seven HyCal crates.
Reads pedestal means from `.cnf` files, parses per-channel RMS from
`faV3peds` stdout, and flags irregular channels.

```bash
hycal_pedestal_monitor          # view existing data
hycal_pedestal_monitor --sim    # test with simulated data
```

**Shift pedestal check (operator procedure).** Pedestals must be
measured before the first DAQ run of each shift, while DAQ is idle.

1. Make sure the DAQ is **stopped**.
2. Launch `hycal_pedestal_monitor`.
3. Click **Measure Pedestals** and confirm.  The tool SSHs to `adchycal1`..`adchycal7` and runs `faV3peds` (a few minutes).
4. Inspect the two maps (left: current mean; right: difference from configured) and the report panel for flagged channels:
   - `DEAD` — avg < 1, rms < 0.1.
   - `OUT OF RANGE` — mean outside 50..300.
   - `HIGH RMS` — sigma > 1.5.
   - `DRIFT` — shifted > 3 counts from configured.
5. Click **Save Report** for the shift log.
6. If new issues appear, notify the run coordinator before starting data taking.  Thresholds are defined at the top of the script.

### hycal_gain_monitor.py

Per-module LMS-based gain factors.  Loads the text
`prad_{run:06d}_LMS.dat` files produced by the offline gain analysis
and plots the HyCal geo map, the LMS reference-channel stability,
a run-to-run drift view, and a table of irregular `(module, run)`
outliers.

```bash
hycal_gain_monitor
```

### hycal_map_builder.py

Generic HyCal geometry viewer that colour-maps user data loaded from
JSON or plain text.  Supports day/night theme toggle, PbGlass alpha
slider, zoom/pan, and palette cycling (click the colour bar).

- **JSON**: `{"<module_name>": {"<field>": <value>, ...}, ...}` — the last entry of a history list is used; nested dicts are flattened with dot notation; non-numeric fields (timestamps) are ignored.
- **Text**: whitespace/comma/tab-delimited rows `<module> <v1> <v2> ...` with optional header.

```bash
hycal_map_builder                     # empty map
hycal_map_builder mydata.json         # auto-load
hycal_map_builder mydata.txt --field rms
```

## GEM ↔ HyCal

### gem_hycal_match_viewer.py

Per-event GEM↔HyCal matching browser.  Shares `analysis/pyscripts/_common.py`
with the offline analysis, so the reconstruction (HyCal waveform →
energy → island clusters; GEM pedestal → CM → ZS → 2-D hits) and the
parametric matching cut match the offline `gem_hycal_matching.py` /
`.C` outputs bit for bit.

Two views side by side:

- **Front view** — HyCal geometry (modules + cluster centroids) with GEM hits projected through the target onto the HyCal plane, colour-coded per detector.  Dashed lines connect each best-matched HC × GEM pair.
- **Side view** (Z-Y) — target / four GEM planes / HyCal face with hit markers and matched-pair lines.

The toolbar adds **Find next ▶▶** driven by two thresholds: minimum
matched hits per detector (*N*) and minimum detectors satisfied (*K*).
The search runs as a foreground scan with a cancellable progress
dialog.  An `nσ` spinbox tweaks the matching window without
re-decoding the current event.

```bash
gem_hycal_match_viewer                       # File → Open…
gem_hycal_match_viewer run.evio.00000        # auto-load
gem_hycal_match_viewer run.evio.00000 -r 23867
```

### coincidence_monitor.py

Connects to a running `prad2_server` (HTTP REST API, port 5051 by
default), iterates every event in the loaded file, and accumulates
per-module coincidence statistics between the four upstream
scintillators (V1–V4) and every HyCal module.  Two event-selection
modes are provided:

- **AND** — `rate(V_i, M) = N(V_i fired AND M is best) / N(M is best)`; only events with at least one scintillator firing are counted.  Reports *P(V_i fired | M was the highest-ADC module)*.
- **OR** — `rate(V_i, M) = N(V_i fired AND M is best) / N(V_i fired)`; all HyCal cluster events counted regardless of scintillators.  Reports *P(M is best | V_i fired)*.

A channel is "fired" when at least one FADC peak height (above
pedestal) exceeds the user threshold.  The bottom half of the window
shows the currently selected scintillator waveform and the HyCal
module last clicked on the map, both fetched from the server for the
event entered in the Event Browser.

```bash
python scripts/coincidence_monitor.py [--url http://HOST:PORT] [--theme dark|light]
```

## Tagger / TDC

### tagger_viewer.py

Viewer for V1190 TDC hits from the tagger crate (ROC `0x008E`,
bank `0xE107`).  Shows a per-slot bar chart of hits/channel, a
single-channel TDC-value histogram, and event-wise correlation tabs
(Δt = A − B and a 2-D `tdc(A)` vs `tdc(B)` heatmap).  The bar chart
auto-sizes its x-axis to 16 / 32 / 64 / 128 channels based on the
highest channel actually hit.  Human-readable counter names come
from `database/tagger_map.json`.

**Offline — EVIO file** (decoded in process by `prad2py`):

```bash
tagger_viewer /data/stage6/prad_023667/prad_023667.evio.00000 \
              -n 200000          # cap physics events (optional)
              --roc 0x8E         # restrict to the tagger ROC (optional)
```

**Online — live ET stream** from a running `prad2_server`.  The server
decodes tagger TDC hits only when at least one client is subscribed,
so regular monitoring is unaffected:

```bash
# DAQ machine (one-time)
prad2_server --online --port 5051

# Viewer (anywhere with PyQt6 + QtWebSockets installed)
tagger_viewer --live ws://clondaq6:5051
```

On startup with `--live`, a fast subscribe/ack round-trip is run
*before* the main window opens — the viewer exits with a clear error
rather than showing an empty window if the server is unreachable or
the protocol does not match.  Pass `--no-smoke-test` to skip it.
The **File** menu adds **Connect to prad2_server…** (Ctrl+L) and
**Disconnect**.  **Pause** / **Clear** sit next to the **Bins**
spinner.  Memory is capped at 10 M hits (rolling — the oldest half
is dropped).

Wire protocol — WebSocket JSON messages `tagger_subscribe` /
`tagger_subscribed` / `tagger_unsubscribe`.  The binary frame format
(useful for writing a different client):

```
magic        "TGR1"   (4 bytes)
flags        u32      (bit 0 = some frames have been dropped)
n_hits       u32
first_seq    u32
last_seq     u32
dropped      u32      (cumulative since server start)
records      n_hits × 16-byte packed BinHit
               u32 event_num, u32 trigger_bits, u16 roc_tag,
               u8 slot, u8 channel_edge (bit 7 = edge, bits 6:0 = channel),
               u32 tdc
```

## Calibration / analysis viewers

### ep_calib_viewer.py

Visualises `epCalib` output.  Reads
`Physics_calib/{run}/calib_iter{N}.json` together with the
companion `CalibResult_iter{N}.root` and renders the per-module
energy distribution, the run-wide ratio, the energy-θ correlation,
the (col, row) event map, and the measured-peak / σ histograms in a
single window.

```bash
python scripts/ep_calib_viewer.py [Physics_calib_dir]
python scripts/ep_calib_viewer.py build/Physics_calib
```

### json_flattener.py

Generic helper that flattens an arbitrary JSON document into a flat
`(key, value)` table.  Nested objects become dotted paths; arrays
use numeric indices.  `--max-depth N` caps recursion (containers at
the cap are kept as compact JSON strings).

```bash
python scripts/json_flattener.py input.json
python scripts/json_flattener.py input.json --max-depth 2
python scripts/json_flattener.py input.json --csv > out.csv
cat input.json | python scripts/json_flattener.py
```

## DAQ dev tools (`scripts/daq_tool/`)

Dev-only GUIs that do not ship with the installation (run them from
the source checkout).  The CMake install step excludes this directory.

### trigger_mask_editor.py

Visual editor for FAV3 trigger masks.  Displays a HyCal geo view
(with LMS1-3, V1-V4 below) and lets you toggle channels off/on by
clicking or dragging modules.  Generates trigger-mask `.cnf` files —
only slots with disabled channels are written; unmapped DAQ channels
(slot positions with no module) are always masked off.

```bash
python scripts/daq_tool/trigger_mask_editor.py                     # start fresh
python scripts/daq_tool/trigger_mask_editor.py -i existing.cnf     # load existing mask
python scripts/daq_tool/trigger_mask_editor.py -o output.cnf       # set default save path
```

### fadc_gain_config.py

Generates a text-based `adchycal_gain.cnf` for the FADC250 DAQ
(`FAV3_ALLCH_GAIN` entries, one per 16-channel slot, grouped by
crate / slot).  Gains come from a calibration JSON
(`-c path.json`, `{"name", "factor"}` entries) or from uniform
values per module type (`--pbwo4-gain` / `--pbglass-gain`).

```bash
python scripts/daq_tool/fadc_gain_config.py
python scripts/daq_tool/fadc_gain_config.py -c database/calibration/adc_to_mev_factors_cosmic.json
python scripts/daq_tool/fadc_gain_config.py --pbwo4-gain 0.15 --pbglass-gain 0.12
```

## Dev one-shot helpers (`scripts/dev_tool/`)

Also not installed — one-off generators kept for reproducibility.

| Script | Purpose |
|---|---|
| `convert_gem_map.py` | Converts the upstream APV mapping text format (`docs/gem/gem_map_prad2*.txt`) into the `gem_map.json` schema we consume. |
| `merge_hycal_map.py` | Joins the legacy `hycal_modules.json` + `hycal_daq_map.json` into the unified `hycal_map.json` schema. |
| `extract_tagger_map.py` | Parses `docs/Tagger_translation_0.xlsx` into `database/tagger_map.json` (the counter-name lookup used by `tagger_viewer`).  Re-run whenever the tagger cabling changes. |

## Shell scripts (`scripts/shell/`)

Not installed as a directory — `prad2_setup.sh` and
`prad2_setup.csh.in` are installed individually to `<prefix>/bin/`
by CMake; everything else is meant to be used from the source
checkout on DAQ / operator machines.

- `prad2_setup.sh` — bash/zsh environment setup, sourced at runtime from `<prefix>/bin/`.
- `prad2_setup.csh.in` — CMake template; `configure_file` bakes the install prefix in and writes `prad2_setup.csh` to `<prefix>/bin/`.
- `run_gain_monitor.sh` — wrapper that parallelises `prad2ana_gain_monitor` across the sub-files of a run and merges the outputs via `hadd`.  Requires the installed analysis binaries on `$PATH` (source `prad2_setup.sh` first).

  ```bash
  scripts/shell/run_gain_monitor.sh <run_number> <num_cpus> [subfile_min] [subfile_max]
  ```

- `start_prad2mon` — tmux-session template for running `prad2_server` under tmux with a log tee.  Copy it, edit the site-specific config block at the top, and `chmod +x`.

## Using prad2py directly

`prad2py` exposes the decoder through a `dec` submodule and the
reconstruction through `det` — useful for custom offline analysis
that goes beyond what the GUIs above do.  Build it once with
`-DBUILD_PYTHON=ON`, then:

```python
import prad2py
from prad2py import dec                 # evio reader + event types

cfg = dec.load_daq_config()             # installed daq_config.json
ch  = dec.EvChannel(); ch.set_config(cfg)
ch.open("/data/.../prad_023671.evio.00000")

while ch.read() == dec.Status.success:
    if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
        continue
    for i in range(ch.get_n_events()):
        ch.select_event(i)              # picks sub-event + clears cache
        info     = ch.info()            # cheapest: TI/trigger metadata
        fadc_evt = ch.fadc()            # FADC250 waveforms, cached
        tdc_evt  = ch.tdc()             # V1190 tagger hits, cached
        # gem_evt = ch.gem();  vtp_evt = ch.vtp();  dsc_evt = ch.dsc()
```

Random-access mode (also used internally by `hycal_event_viewer` and
`prad2_server`'s file mode):

```python
ch = dec.EvChannel(); ch.set_config(cfg)
ch.open_random_access("run.evio.00000")
n = ch.get_random_access_event_count()
# jump to any evio event in O(1)
for i in (0, n // 2, n - 1):
    ch.read_event_by_index(i)
    ch.scan(); ch.select_event(0)
    print(int(ch.info().event_number))
```

Full reconstruction helpers live in `det`:

```python
from prad2py import dec, det

cfg = dec.load_daq_config()
ch  = dec.EvChannel(); ch.set_config(cfg); ch.open("run.evio.00000")

# --- GEM reconstruction -----------------------------------------------
gsys = det.GemSystem()
gsys.init("database/gem_map.json")
gsys.load_pedestals("database/gem_ped.json")    # optional
gcl  = det.GemCluster()

# --- HyCal reconstruction ---------------------------------------------
hsys = det.HyCalSystem()
hsys.init("database/hycal_map.json")
hsys.load_calibration("database/hycal_calib.json")
hcl  = det.HyCalCluster(hsys)

while ch.read() == dec.Status.success:
    if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
        continue
    for i in range(ch.get_n_events()):
        ch.select_event(i)

        # GEM 2-D hits
        gsys.clear()
        gsys.process_event(ch.gem())
        gsys.reconstruct(gcl)
        for h in gsys.get_all_hits():
            print("GEM", h.det_id, h.x, h.y, h.x_charge, h.y_charge)

        # HyCal — feed per-module energies (e.g. from ch.fadc() + your
        # calibration), then cluster:
        # hcl.clear()
        # for module_idx, energy_mev, time_ns in my_hycal_hits(ch.fadc()):
        #     hcl.add_hit(module_idx, energy_mev, time_ns)
        # hcl.form_clusters()
        # for c in hcl.reconstruct_hits():
        #     print("HyCal", c.center_id, c.x, c.y, c.energy, c.time)
```

`prad2py.load_tdc_hits(path, ...)` remains the convenience entry point
for the common "one-shot flat table of hits" workflow.

## GEM tools

GEM-specific scripts and the `gem_dump` C++ binary live in the
top-level [`gem/`](../gem/README.md) directory — that README also
contains the GEM detector reference notes.

## Tagger ↔ HyCal coincidence (ROOT macro)

See `analysis/scripts/tagger_hycal_correlation.C` — a self-contained
ROOT/ACLiC macro that builds ΔT histograms for (T10R, E49…E58) pairs,
Gauss-fits each coincidence peak, applies a ±*N*σ timing cut, and plots
the W1156 peak height/integral for the selected events.

```bash
cd build
root -l ../analysis/scripts/rootlogon.C
.x ../analysis/scripts/tagger_hycal_correlation.C+( \
     "/data/stage6/prad_023671/prad_023671.evio.00000", \
     "tagger_w1156_corr.root", 500000)
```
