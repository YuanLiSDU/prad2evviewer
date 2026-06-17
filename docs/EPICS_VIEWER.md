# EPICS Viewer

`analysis/tools/epics_viewer.py` is a pure PyQt6 desktop viewer for the
slow-control side trees written by the PRad-II replay tools.  It reads
`epics`, `scalers`, and `runinfo` from one replayed ROOT file and shows
four interactive scatter plots on a shared run-time axis.

The viewer is meant for quick run-quality inspection: beam current and
position, target pressure, DSC2 livetime, scaler counters, and CODA
control-event metadata can be checked without producing a
`replay_filter` JSON report first.

# Quick Start

```bash
python3 analysis/tools/epics_viewer.py /path/to/prad_025215_filter.root
```

From the GUI, use `Open ROOT file...` to load another file and `Reload`
to reread the current one.

Optional fallback:

```bash
python3 analysis/tools/epics_viewer.py --uproot /path/to/file.root
```

`--uproot` enables an uproot/awkward fallback if the default ROOT-based
readers are not available.  The normal reader keeps ROOT outside the Qt
GUI process: it first tries the external `root` executable, then a
short-lived PyROOT subprocess.

Last-resort fallback:

```bash
python3 analysis/tools/epics_viewer.py --pyroot /path/to/file.root
```

`--pyroot` enables the older in-process PyROOT reader.  Use it only when
the external `root` command, the PyROOT subprocess reader, and the
uproot fallback are unavailable.

# Dependencies

The viewer intentionally avoids Qt WebEngine.  It uses the same style as
the other scripts viewers: standard Qt widgets and QPainter drawing.

Required:

| Dependency | Purpose |
|---|---|
| `python3` | Runtime |
| `PyQt6` | GUI widgets and drawing |
| `root` command or PyROOT | ROOT file reader |

Optional:

| Dependency | Purpose |
|---|---|
| `uproot` + `awkward` | Fallback ROOT reader when `--uproot` is enabled |

On machines with mixed Qt installations, the script prefers PyQt6's
bundled Qt libraries before importing PyQt6.  This avoids common
`Qt_6_PRIVATE_API` symbol errors from loading mismatched system Qt
libraries.

# Input ROOT Trees

The viewer reads the replay side trees documented in
[`REPLAYED_DATA.md`](REPLAYED_DATA.md):

| Tree | Use in the viewer |
|---|---|
| `epics` | Slow-control PV name/value updates |
| `scalers` | DSC2 counters and derived DAQ livetime |
| `runinfo` | CODA PRESTART/GO/END metadata and DAQ-config presence |

All numeric series are aligned to a common x-axis:

```text
x = (unix_time - earliest_valid_unix_time) / 60
```

The top-right page header shows the absolute Unix-time range as local
date/time, while every plot uses minutes from the first slow-control
timestamp.

# Interface

The main page contains four plot panels arranged side by side.  Each
panel has its own channel selector and zoom state.

| Control | Behavior |
|---|---|
| Channel combo box | Selects the quantity plotted in that panel |
| Type in combo box | Case-insensitive substring search over all channels |
| Left-drag inside plot | Zooms to the selected x/y rectangle |
| Right-click plot | Resets that plot to the full range |
| Double-click plot | Also resets that plot to the full range |
| `Open ROOT file...` | Opens another replayed ROOT file |
| `Reload` | Rereads the currently selected file |
| `allow uproot fallback` | Allows the slower uproot path if ROOT is unavailable |

Each plot draws all available points for the selected channel as a
scatter plot.  Points whose source row has `good == false` are drawn in
red; all other points are drawn in green.  The `good` flag is present in
`prad2ana_replay_filter` outputs and absent in unfiltered replay files.

## Plot Contents

Each panel shows one numeric time series from the selected ROOT file.
The selected channel name is shown in the combo box above the plot.

| Plot element | Meaning |
|---|---|
| x-axis | Minutes from the first valid slow-control Unix timestamp in the file |
| y-axis | Numeric value of the selected channel |
| green point | A normal point, or a point from a file/tree without a `good` flag |
| red point | A point whose source slow-control row has `good == false` |
| top-right page date range | Absolute local date/time span corresponding to the plotted run |
| `points | full` | Number of points in this channel; the plot is at the full range |
| `points | zoomed` | Number of points in this channel; the plot is currently zoomed |

In `prad2ana_replay_filter` outputs, red points mark slow-control
checkpoints rejected by the filter's configured cuts.  In unfiltered
`prad2ana_replay_rawdata` or `prad2ana_replay_recon` files, no `good`
branch exists, so all points are treated as normal and drawn green.

## Zooming

Zooming is per panel: changing or resetting one plot does not affect the
other three plots.

| Action | Result |
|---|---|
| Left-click and drag inside the plot area | Draws a zoom rectangle |
| Release the left mouse button | Zooms to the selected x/y rectangle |
| Right-click inside the plot | Resets that plot to the full range |
| Double-click inside the plot | Also resets that plot to the full range |
| Select another channel in the combo box | Resets that panel's zoom and draws the new channel |

The zoom rectangle changes both x and y ranges.  Use a narrow horizontal
drag to inspect a short time interval, or a rectangle around a cluster of
points to inspect both time and value structure.

# Channel Names

The combo boxes contain a flat list of plottable numeric series.  Prefixes
identify their source tree.

## EPICS Channels

EPICS PVs keep their original names.  Examples:

| Example | Meaning |
|---|---|
| `hallb_IPM2C21A_CUR` | Beam current PV |
| `hallb_IPM2C21A_XPOS` | Beam x-position PV |
| `hallb_IPM2C21A_YPOS` | Beam y-position PV |
| `TGT:PRad:Cell_P` | Target-cell pressure PV |

EPICS rows are sparse: a PV only appears when that EPICS event updated
it.  The viewer therefore plots the actual update points and does not
forward-fill EPICS values.

## Scaler Channels

`scalers:` channels come from the DSC2 scaler tree.

| Channel | Meaning |
|---|---|
| `scalers:gated` | Selected DSC2 gated counter |
| `scalers:ungated` | Selected DSC2 ungated counter |
| `scalers:live_ratio` | Cumulative `gated / ungated` stored in the tree |
| `scalers:ref_gated` | Reference gated counter |
| `scalers:ref_ungated` | Reference ungated counter |
| `scalers:trg_gated[N]` | TRG channel `N` gated counter |
| `scalers:trg_ungated[N]` | TRG channel `N` ungated counter |
| `scalers:tdc_gated[N]` | TDC channel `N` gated counter |
| `scalers:tdc_ungated[N]` | TDC channel `N` ungated counter |
| `scalers:event_number` | Physics event number carrying the scaler readout |
| `scalers:slot` | DSC2 slot |
| `scalers:source` | Livetime source selector: ref/trg/tdc encoded as 0/1/2 |
| `scalers:channel` | Selected DSC2 channel for trg/tdc sources |

The full 16 channels of each DSC2 counter array are exposed as separate
series, for example `scalers:trg_gated[0]` through
`scalers:trg_gated[15]`.

## Derived DAQ Livetime

Two DAQ livetime channels are provided:

| Channel | Definition |
|---|---|
| `DAQ live time [%]` | Slice-local `delta_gated / delta_ungated * 100` |
| `DAQ cumulative live time [%]` | `scalers:live_ratio * 100` |

`DAQ live time [%]` matches the livetime shown in `replay_filter`
reports and in `replay_report_viewer.py`.  It is computed from adjacent
DSC2 scaler rows sorted by `scalers:event_number`.  If the counter moves
backward, the baseline is reset to `(0, 0)`, matching
`replay_filter::compute_delta_live_pct()`.  Rows with no ungated
advance are skipped.

This is the quantity to use for normal DAQ live-time quality checks.  It
should typically sit near the expected run value, for example around
91 percent in stable production running.

## Runinfo Channels

`runinfo:` channels come from CODA control events.  A run normally has
PRESTART, GO, and END rows.

| Channel | Meaning |
|---|---|
| `runinfo:event_tag` | CODA control tag: `0x11` PRESTART, `0x12` GO, `0x14` END |
| `runinfo:run_type` | CODA run type, usually nonzero on PRESTART |
| `runinfo:run_number` | CODA run number |
| `runinfo:has_daq_config` | 1 when the row carries DAQ config text, else 0 |
| `runinfo:daq_config_bytes` | Size of the DAQ config text in bytes |

The raw `daq_config` branch is large text and is not plotted directly.
The viewer exposes whether it is present and how large it is.

# Default Plots

When a file is loaded, the viewer tries to choose useful defaults:

1. `DAQ live time [%]`
2. `scalers:live_ratio`
3. `hallb_IPM2C21A_CUR`
4. `hallb_IPM2C21A_XPOS`

If a preferred channel is absent, the next available channel is used.

# Reading Strategy

The loader tries readers in this order:

1. The command-line `root` executable, using a temporary C++ macro that
   dumps the plottable values to TSV.
2. A short-lived PyROOT subprocess that dumps the same TSV format and
   exits with `os._exit(0)` to avoid ROOT shutdown destructors in the Qt
   GUI process.
3. uproot/awkward, only when `--uproot` or `allow uproot fallback` is
   enabled.
4. PyROOT in the current Python process, only when `--pyroot` is enabled.

Loading runs in a `QThread`, so the GUI remains responsive while a file
is being read.

# Troubleshooting

## PyQt6 import fails with a Qt private-API symbol error

This usually means PyQt6 and system Qt libraries are mixed.  The viewer
attempts to prepend PyQt6's bundled Qt library directory to
`LD_LIBRARY_PATH` and re-exec itself before importing PyQt6.  This is the
same workaround used by other PyQt6 viewers in `scripts/`.

## Closing the viewer prints a ROOT segmentation violation

Do not use `--pyroot` unless it is necessary.  Some ROOT/PyROOT builds
can crash while the Python process exits and ROOT global objects are
destroyed after Qt has already started shutting down.  The normal command

```bash
python3 analysis/tools/epics_viewer.py /path/to/file.root
```

uses either the external `root` command or a PyROOT subprocess to read
the file, so ROOT is kept out of the GUI process and the window can close
cleanly.

## The file opens but a selected channel has no points

Check the page summary for the total point count, then try searching the
combo box for another channel.  Some EPICS PVs only update in certain
runs, and `runinfo` has only a few points by design.

## DAQ live time does not match `scalers:live_ratio`

That is expected.  `scalers:live_ratio` is cumulative.  `DAQ live time
[%]` is slice-local and is the report-viewer quantity:

```text
DAQ live time [%] = (gated_i - gated_{i-1}) /
                    (ungated_i - ungated_{i-1}) * 100
```

Use `DAQ live time [%]` for the standard livetime trace.
