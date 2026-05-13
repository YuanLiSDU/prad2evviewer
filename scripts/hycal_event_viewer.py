#!/usr/bin/env python3
"""
HyCal Event Viewer
==================
Browses an evio file event-by-event.  Opens the file in evio's
random-access mode (no event processing at open time), indexes physics
sub-events, and lets the user step through them.  Two tabs:

* **Waveform** — per-module FADC display, stacked waveform, four
  accumulating histograms (peak height, integral, time, n-peaks).
  "Process next 10k" fills histograms in a background pass.
* **Cluster** — HyCal heatmap of per-module energy with cluster
  overlays (crosshair + energy label), cluster table, selector.
  Clustering uses ``prad2py.det.HyCalCluster`` on live ADC data;
  calibration comes from ``HyCalSystem::Init``.

Usage
-----
    python scripts/hycal_event_viewer.py RUN.evio.00000
    python scripts/hycal_event_viewer.py             # File → Open…
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtWidgets import (
    QAbstractItemView, QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout, QGridLayout,
    QLabel, QComboBox, QCheckBox, QCompleter, QFileDialog, QMessageBox,
    QProgressDialog, QSizePolicy, QStatusBar, QToolTip, QPushButton,
    QSpinBox, QDoubleSpinBox, QSplitter, QTabWidget, QTableWidget,
    QTableWidgetItem, QHeaderView, QDockWidget, QGroupBox,
    QDialog, QDialogButtonBox, QFrame, QLineEdit,
)
from PyQt6.QtCore import (
    Qt, QObject, QPointF, QRectF, QThread, pyqtSignal, QTimer,
)
from PyQt6.QtGui import (
    QAction, QBrush, QKeySequence, QPainter, QColor, QPen, QFont, QPolygonF,
    QShortcut, QDoubleValidator,
)

from hycal_geoview import (
    load_modules as load_geo_modules,
    HyCalMapWidget, ColorRangeController, cmap_qcolor,
    apply_theme_palette, set_theme,
    available_themes, THEME, themed,
)


# ===========================================================================
#  prad2py discovery (mirrors tagger_viewer.py)
# ===========================================================================

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR   = _SCRIPT_DIR.parent

for _cand in (
    _REPO_DIR / "build" / "python",
    _REPO_DIR / "build-release" / "python",
    _REPO_DIR / "build" / "Release" / "python",
):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

try:
    import prad2py                                # type: ignore
    _HAVE_PRAD2PY = True
    _PRAD2PY_ERR  = ""
except Exception as _exc:
    prad2py = None                                # type: ignore
    _HAVE_PRAD2PY = False
    _PRAD2PY_ERR  = f"{type(_exc).__name__}: {_exc}"


def _check_evchannel_support() -> Optional[str]:
    """Return None if the loaded prad2py has the expected EvChannel API,
    else an error string suitable for showing to the user."""
    if not _HAVE_PRAD2PY:
        return _PRAD2PY_ERR
    try:
        ch = prad2py.dec.EvChannel()
        if not hasattr(ch, "open_auto"):
            return ("prad2py is missing open_auto — rebuild prad2py after "
                    "the EvChannel.cpp changes:\n"
                    "  cmake -DBUILD_PYTHON=ON -S . -B build && "
                    "cmake --build build --target prad2py")
    except Exception as e:
        return f"{type(e).__name__}: {e}"
    return None


# ===========================================================================
#  WaveAnalyzer — direct exports of the C++ implementation in prad2py.dec.
#  The server uses the same code; ~50× faster than the previous Python port
#  which mattered a lot when "Accumulate all modules" analyses 1700+
#  channels per event.
# ===========================================================================

WaveConfig = prad2py.dec.WaveConfig if _HAVE_PRAD2PY else None
Peak       = prad2py.dec.Peak       if _HAVE_PRAD2PY else None


def analyze(samples, cfg):
    """Run the C++ WaveAnalyzer on one channel's samples.
    Returns ``(ped_mean, ped_rms, peaks_list)``."""
    return prad2py.dec.WaveAnalyzer(cfg).analyze(samples)


# ===========================================================================
#  Histogram accumulator
# ===========================================================================

@dataclass
class Hist1D:
    nbins:   int
    bmin:    float = 0.0
    bstep:   float = 1.0
    bins:    np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    under:   int = 0
    over:    int = 0

    def __post_init__(self):
        if self.bins.size == 0:
            self.bins = np.zeros(self.nbins, dtype=np.int64)

    def fill(self, v: float):
        if v < self.bmin:
            self.under += 1
            return
        b = int((v - self.bmin) / self.bstep)
        if b >= self.nbins:
            self.over += 1
            return
        self.bins[b] += 1

    def reset(self):
        self.bins[:] = 0
        self.under = 0
        self.over = 0

    def to_json(self) -> Dict:
        return {"bins": self.bins.tolist(),
                "underflow": int(self.under), "overflow": int(self.over)}


@dataclass
class ChannelHists:
    roc:         int
    slot:        int
    channel:     int
    module:      Optional[str]
    events:      int = 0
    peak_events: int = 0
    height:      Optional[Hist1D] = None
    integral:    Optional[Hist1D] = None
    position:    Optional[Hist1D] = None
    npeaks:      Optional[Hist1D] = None


# ===========================================================================
#  Waveform peak filter — mirrors C++ PeakFilter (src/app_state.h).
# ===========================================================================
#
# Each axis is an optional [min, max] range; missing bound = no constraint.
# Quality bits use accept / reject masks resolved against PEAK_QUALITY_BITS.
# `enable=False` makes the filter a no-op — driven by the GUI "apply" toggle.
# Same JSON shape as the web monitor's `waveform_filter`, so the same
# monitor_config.json `waveform.filter` block configures both viewers.

# Mirrors AppState::peak_quality_bits_def in src/app_state_init.cpp.  Bit
# values must agree with prad2py.dec.Q_PEAK_* (we read those at runtime so
# the names stay authoritative if the C++ side adds another flag).
def _resolve_quality_bits() -> List[Dict[str, object]]:
    if not _HAVE_PRAD2PY:
        return []
    out: List[Dict[str, object]] = []
    for name, label in (("PILED", "Pile-up"), ("DECONVOLVED", "Deconvolved")):
        bit = getattr(prad2py.dec, f"Q_PEAK_{name}", None)
        if bit is None:
            continue
        # Q_PEAK_* are masks (1 << bit_index).  Store the mask directly —
        # the filter compares peak.quality & mask, no shift needed.
        out.append({"mask": int(bit), "name": name, "label": label})
    return out


PEAK_QUALITY_BITS: List[Dict[str, object]] = _resolve_quality_bits()


@dataclass
class WaveformFilter:
    """Per-peak filter for the Waveform tab — replicates ``PeakFilter`` from
    ``src/app_state.h`` so monitor_config.json's ``waveform.filter`` JSON
    drives both the Python viewer and the live web monitor identically.

    A ``None`` bound means "no constraint" (matches the web monitor's empty
    input fields).  ``enable=False`` short-circuits ``passes()`` so the user
    can toggle filtering off without losing their tuned values.
    """
    enable:       bool = True
    time_min:     Optional[float] = None
    time_max:     Optional[float] = None
    integral_min: Optional[float] = None
    integral_max: Optional[float] = None
    height_min:   Optional[float] = None
    height_max:   Optional[float] = None
    q_accept:     int = 0   # bitmask, 0 = accept any
    q_reject:     int = 0   # bitmask, 0 = reject none

    def passes(self, pk) -> bool:
        if not self.enable:
            return True
        if self.time_min     is not None and pk.time     < self.time_min:     return False
        if self.time_max     is not None and pk.time     > self.time_max:     return False
        if self.integral_min is not None and pk.integral < self.integral_min: return False
        if self.integral_max is not None and pk.integral > self.integral_max: return False
        if self.height_min   is not None and pk.height   < self.height_min:   return False
        if self.height_max   is not None and pk.height   > self.height_max:   return False
        q = int(getattr(pk, "quality", 0))
        if self.q_reject and (q & self.q_reject):    return False
        if self.q_accept and not (q & self.q_accept): return False
        return True

    def copy(self) -> "WaveformFilter":
        return WaveformFilter(**self.__dict__)

    @classmethod
    def from_json(cls, j: Optional[Dict]) -> "WaveformFilter":
        """Parse a ``waveform.filter`` JSON object into a filter.  Mirrors
        ``PeakFilter::parse``; unknown keys are ignored so future schema
        additions don't break loading."""
        f = cls()
        j = j or {}

        def _axis(key: str):
            ax = j.get(key) or {}
            lo = ax.get("min")
            hi = ax.get("max")
            return (float(lo) if lo is not None else None,
                    float(hi) if hi is not None else None)

        f.time_min,     f.time_max     = _axis("time")
        f.integral_min, f.integral_max = _axis("integral")
        f.height_min,   f.height_max   = _axis("height")

        qb = j.get("quality_bits") or {}
        f.q_accept = _names_to_mask(qb.get("accept") or [])
        f.q_reject = _names_to_mask(qb.get("reject") or [])
        return f

    def to_json(self) -> Dict:
        out: Dict = {}

        def _axis(lo, hi):
            ax: Dict = {}
            if lo is not None: ax["min"] = lo
            if hi is not None: ax["max"] = hi
            return ax

        if (a := _axis(self.time_min,     self.time_max))     : out["time"]     = a
        if (a := _axis(self.integral_min, self.integral_max)) : out["integral"] = a
        if (a := _axis(self.height_min,   self.height_max))   : out["height"]   = a
        if self.q_accept or self.q_reject:
            out["quality_bits"] = {
                "accept": _mask_to_names(self.q_accept),
                "reject": _mask_to_names(self.q_reject),
            }
        return out


def _names_to_mask(names) -> int:
    m = 0
    for n in names or ():
        for d in PEAK_QUALITY_BITS:
            if d["name"] == n:
                m |= int(d["mask"])
                break
    return m


def _mask_to_names(mask: int) -> List[str]:
    if not mask:
        return []
    return [d["name"] for d in PEAK_QUALITY_BITS if mask & int(d["mask"])]


# ===========================================================================
#  Config / map loaders
# ===========================================================================

def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    """(crate, slot, channel) -> module_name from hycal_map.json.

    Records without a "daq" block (boosters, PRad-1 V1-V4) are skipped.
    """
    with open(path, encoding="utf-8") as f:
        entries = json.load(f)
    out: Dict[Tuple[int, int, int], str] = {}
    for e in entries:
        d = e.get("daq")
        if not d:
            continue
        out[(int(d["crate"]), int(d["slot"]), int(d["channel"]))] = e["n"]
    return out


def load_roc_tag_map(path: Path) -> Dict[int, int]:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    out = {}
    for r in cfg.get("roc_tags", []):
        tag = int(r["tag"], 16) if isinstance(r["tag"], str) else int(r["tag"])
        if r.get("type") == "roc":
            out[tag] = int(r["crate"])
    return out


def load_hist_config(path: Path) -> Dict:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return cfg.get("waveform", {})


def load_trigger_bit_map(path: Path) -> Dict[str, int]:
    if not path.is_file():
        return {}
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return {e["name"]: int(e["bit"]) for e in cfg.get("trigger_bits", [])}


def _mask_from_names(names: List[str], bitmap: Dict[str, int]) -> int:
    m = 0
    for n in names:
        if n in bitmap:
            m |= (1 << bitmap[n])
        else:
            print(f"  warning: trigger bit {n!r} not in trigger_bits.json",
                  file=sys.stderr)
    return m


# ===========================================================================
#  Hist bin builder helpers
# ===========================================================================

def _nbins(c: Dict) -> int:
    span = c["max"] - c["min"]
    return max(1, int(np.ceil(span / c["step"])))


def _make_hists(h_cfg: Dict, i_cfg: Dict, p_cfg: Dict, n_cfg: Dict,
                roc: int, slot: int, channel: int,
                module: Optional[str]) -> ChannelHists:
    return ChannelHists(
        roc=roc, slot=slot, channel=channel, module=module,
        height  =Hist1D(_nbins(h_cfg), h_cfg["min"], h_cfg["step"]),
        integral=Hist1D(_nbins(i_cfg), i_cfg["min"], i_cfg["step"]),
        position=Hist1D(_nbins(p_cfg), p_cfg["min"], p_cfg["step"]),
        npeaks  =Hist1D(_nbins(n_cfg), n_cfg["min"], n_cfg["step"]),
    )


# ===========================================================================
#  Indexer — background pass to locate all physics sub-events
# ===========================================================================

class IndexerWorker(QObject):
    """Scans the file once in RA mode to record (evio_idx, sub_idx) per
    physics sub-event.  No waveform decoding — Scan() only."""

    progressed = pyqtSignal(int, int)   # (evio_events_scanned, total_evio_events)
    finished   = pyqtSignal(object)     # {"path": str, "index": list, "total_evio": int}
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str, daq_config_path: str):
        super().__init__()
        self._path = evio_path
        self._daq_cfg_path = daq_config_path
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            self.finished.emit(self._run())
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> Dict:
        dec = prad2py.dec
        cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
               else dec.load_daq_config())
        ch  = dec.EvChannel()
        ch.set_config(cfg)
        st = ch.open_auto(self._path)
        if st != dec.Status.success:
            raise RuntimeError(f"cannot open {self._path}: {st}")
        is_ra = ch.is_random_access()

        # In RA mode we know the total upfront; sequential mode walks to EOF
        # so we just report a rolling count.
        total_evio = ch.get_random_access_event_count() if is_ra else 0
        index: List[Tuple[int, int]] = []

        progress_every = max(1, total_evio // 200) if total_evio else 500
        ei = 0
        while ch.read() == dec.Status.success:
            if self._cancel:
                break
            if ch.scan() and ch.get_event_type() == dec.EventType.Physics:
                for si in range(ch.get_n_events()):
                    index.append((ei, si))
            ei += 1
            if (ei % progress_every) == 0:
                self.progressed.emit(ei, total_evio or ei)

        ch.close()
        self.progressed.emit(ei, total_evio or ei)
        return {"path": self._path, "index": index,
                "total_evio": ei, "cancelled": self._cancel,
                "random_access": is_ra}


# ===========================================================================
#  Batch processor — fills all-module hists for the next N events
# ===========================================================================

class BatchWorker(QObject):
    """Reads events start_idx .. start_idx + n - 1 (no display updates) and
    fills histograms for every channel present.  Runs in its own thread with
    its own EvChannel handle (separate from the UI's)."""

    progressed = pyqtSignal(int, int, int)  # (done, target, peaks_found)
    finished   = pyqtSignal(int)            # events_processed
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str, daq_config_path: str,
                 index: List[Tuple[int, int]],
                 start_idx: int, count: int,
                 channels: Dict[Tuple[int, int, int], ChannelHists],
                 wcfg: WaveConfig, peak_filter: "WaveformFilter",
                 accept_mask: int, reject_mask: int,
                 accumulated: Optional[np.ndarray] = None):
        super().__init__()
        self._path = evio_path
        self._daq_cfg_path = daq_config_path
        self._index = index
        self._start = start_idx
        self._count = count
        self._channels = channels
        self._wcfg = wcfg
        # Snapshot the filter so user edits during the batch don't change
        # the cuts mid-run (each batch should fill against a stable cut).
        self._filter = peak_filter.copy()
        self._accept = accept_mask
        self._reject = reject_mask
        self._accumulated = accumulated     # shared bool array, mutated in place
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            self.finished.emit(self._run())
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> int:
        dec = prad2py.dec
        cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
               else dec.load_daq_config())
        ch  = dec.EvChannel()
        ch.set_config(cfg)
        st = ch.open_auto(self._path)
        if st != dec.Status.success:
            raise RuntimeError(f"cannot open {self._path}: {st}")
        is_ra = ch.is_random_access()

        n_done = 0
        peaks_found = 0
        progress_every = max(1, self._count // 100)
        wcfg = self._wcfg
        flt = self._filter
        channels = self._channels

        def _fold_event(phys_idx: int, sub_idx: int) -> int:
            """Scan + select current event, accumulate hists for every
            channel, return peaks found this event (or -1 if the event
            is rejected by trigger mask or dedup)."""
            nonlocal n_done
            if (self._accumulated is not None
                    and 0 <= phys_idx < self._accumulated.size
                    and self._accumulated[phys_idx]):
                n_done += 1
                return -1
            if not ch.scan():
                return -1
            ch.select_event(sub_idx)
            info = ch.info()
            tb = int(info.trigger_bits)
            if self._accept and (tb & self._accept) == 0:
                n_done += 1
                return -1
            if self._reject and (tb & self._reject):
                n_done += 1
                return -1
            fadc_evt = ch.fadc()
            pfound = 0
            for r in range(fadc_evt.nrocs):
                roc = fadc_evt.roc(r)
                roc_tag = int(roc.tag)
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        key = (roc_tag, s, c)
                        hits = channels.get(key)
                        if hits is None:
                            continue
                        samples = slot.channel(c).samples
                        if samples.size < 10:
                            continue
                        _, _, peaks = analyze(samples, wcfg)
                        np_kept = 0
                        for p in peaks:
                            if not flt.passes(p):
                                continue
                            hits.height.fill(p.height)
                            hits.integral.fill(p.integral)
                            hits.position.fill(p.time)
                            np_kept += 1
                            pfound += 1
                        hits.npeaks.fill(np_kept)
                        hits.events += 1
                        if np_kept > 0:
                            hits.peak_events += 1
            if self._accumulated is not None and 0 <= phys_idx < self._accumulated.size:
                self._accumulated[phys_idx] = True
            n_done += 1
            return pfound

        try:
            if is_ra:
                # RA: jump directly to each phys event's evio block.
                for i in range(self._count):
                    if self._cancel: break
                    phys_idx = self._start + i
                    if phys_idx >= len(self._index): break
                    ev_idx, sub_idx = self._index[phys_idx]
                    if ch.read_event_by_index(ev_idx) != dec.Status.success:
                        continue
                    pf = _fold_event(phys_idx, sub_idx)
                    if pf > 0: peaks_found += pf
                    if (i % progress_every) == 0:
                        self.progressed.emit(n_done, self._count, peaks_found)
            else:
                # Sequential: walk forward through the file, processing the
                # index entries in order.  self._index is already in
                # evio-order so consecutive phys entries only ever require
                # more Read()s, never a rewind.
                cur_evio = -1
                for i in range(self._count):
                    if self._cancel: break
                    phys_idx = self._start + i
                    if phys_idx >= len(self._index): break
                    need_evio, sub_idx = self._index[phys_idx]
                    while cur_evio < need_evio:
                        if ch.read() != dec.Status.success:
                            raise RuntimeError(
                                f"EOF before reaching evio event {need_evio}")
                        cur_evio += 1
                    pf = _fold_event(phys_idx, sub_idx)
                    if pf > 0: peaks_found += pf
                    if (i % progress_every) == 0:
                        self.progressed.emit(n_done, self._count, peaks_found)
        finally:
            ch.close()

        self.progressed.emit(n_done, self._count, peaks_found)
        return n_done


def _find_channel_samples(fadc_evt, roc_tag: int, slot: int, channel: int):
    """Return samples array for (roc, slot, ch), or None if not present."""
    for r in range(fadc_evt.nrocs):
        roc = fadc_evt.roc(r)
        if int(roc.tag) != roc_tag:
            continue
        if slot not in roc.present_slots():
            continue
        slot_data = roc.slot(slot)
        if channel not in slot_data.present_channels():
            continue
        return slot_data.channel(channel).samples
    return None


# ===========================================================================
#  Small themed overlay controls for plot widgets
# ===========================================================================

def _overlay_checkbox_qss() -> str:
    """QSS for a compact checkbox drawn on top of a plot canvas."""
    return (
        f"QCheckBox{{color:{THEME.TEXT_DIM};background:{THEME.PANEL};"
        f"padding:2px 6px;border:1px solid {THEME.BORDER};border-radius:6px;}}"
        f"QCheckBox:hover{{color:{THEME.TEXT};"
        f"border:1px solid {THEME.ACCENT};}}"
        f"QCheckBox:checked{{color:{THEME.TEXT};}}"
        f"QCheckBox::indicator{{width:12px;height:12px;"
        f"border:1px solid {THEME.BORDER};border-radius:3px;"
        f"background:{THEME.BG};}}"
        f"QCheckBox::indicator:hover{{border:1px solid {THEME.ACCENT};}}"
        f"QCheckBox::indicator:checked{{background:{THEME.ACCENT};"
        f"border:1px solid {THEME.ACCENT};}}"
    )


def _overlay_button_qss() -> str:
    """QSS for a compact pushbutton drawn on top of a plot canvas."""
    return (
        f"QPushButton{{color:{THEME.TEXT_DIM};background:{THEME.PANEL};"
        f"padding:2px 8px;border:1px solid {THEME.BORDER};border-radius:6px;"
        f"font:bold 9pt Monospace;}}"
        f"QPushButton:hover{{color:{THEME.TEXT};"
        f"border:1px solid {THEME.ACCENT};}}"
        f"QPushButton:disabled{{color:{THEME.TEXT_MUTED};}}"
    )


# ===========================================================================
#  Hist1DWidget — QPainter bar chart with optional log Y
# ===========================================================================

class Hist1DWidget(QWidget):
    PAD_L, PAD_R, PAD_T, PAD_B = 58, 14, 20, 20

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(140)
        self._bins: np.ndarray = np.zeros(0, dtype=np.int64)
        self._bmin: float = 0.0
        self._bstep: float = 1.0
        self._under: int = 0
        self._over: int = 0
        self._title: str = ""
        self._xlabel: str = ""
        self._color: QColor = QColor(THEME.ACCENT)
        self._log_y: bool = False
        self._hover_idx: int = -1

        self._logy_cb = QCheckBox("log Y", self)
        self._logy_cb.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._logy_cb.setStyleSheet(_overlay_checkbox_qss())
        self._logy_cb.toggled.connect(self.set_log_y)
        self._logy_cb.adjustSize()
        self._logy_cb.raise_()

    def set_data(self, bins, bmin: float, bstep: float,
                 under: int = 0, over: int = 0,
                 title: str = "", xlabel: str = "",
                 color: Optional[str] = None):
        self._bins = np.asarray(bins, dtype=np.int64)
        self._bmin = float(bmin)
        self._bstep = float(bstep)
        self._under = int(under)
        self._over  = int(over)
        self._title = title
        self._xlabel = xlabel
        if color:
            self._color = QColor(color)
        self._hover_idx = -1
        self.update()

    def set_log_y(self, on: bool):
        if on != self._log_y:
            self._log_y = on
            self.update()

    def clear(self, title: str = ""):
        self._bins = np.zeros(0, dtype=np.int64)
        self._title = title
        self._under = self._over = 0
        self._hover_idx = -1
        self.update()

    def _plot_rect(self) -> QRectF:
        w, h = self.width(), self.height()
        return QRectF(self.PAD_L, self.PAD_T,
                      max(1.0, w - self.PAD_L - self.PAD_R),
                      max(1.0, h - self.PAD_T - self.PAD_B))

    def resizeEvent(self, ev):
        cb = self._logy_cb
        cb.adjustSize()
        cb.move(self.width() - cb.width() - 6, 4)
        super().resizeEvent(ev)

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = self._plot_rect()
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        if self._title:
            f = QFont("Monospace", 10); f.setBold(True)
            p.setFont(f)
            p.setPen(QColor(THEME.TEXT))
            p.drawText(int(r.left()), int(r.top() - 8), self._title)

        n = self._bins.size
        if n == 0 or self._bins.sum() == 0:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "(no data)")
            return

        if self._log_y:
            vals = np.where(self._bins > 0,
                            np.log10(self._bins.astype(np.float64)), 0.0)
            ymin, ymax = 0.0, float(vals.max())
        else:
            vals = self._bins.astype(np.float64)
            ymin, ymax = 0.0, float(vals.max())
        if ymax <= 0:
            ymax = 1.0

        bar_w = r.width() / n
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(self._color)
        for i in range(n):
            v = vals[i]
            if v <= 0:
                continue
            h = (v - ymin) / (ymax - ymin) * r.height()
            if h < 1.0:
                continue
            x0 = r.left() + i * bar_w
            y0 = r.bottom() - h
            p.fillRect(QRectF(x0, y0, max(bar_w, 1.0), h), self._color)

        if 0 <= self._hover_idx < n:
            x0 = r.left() + self._hover_idx * bar_w
            p.setPen(QPen(QColor(THEME.SELECT_BORDER), 1.2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(QRectF(x0, r.top(), max(bar_w, 1.0), r.height()))

        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            if self._log_y:
                val = 10 ** (ymin + frac * (ymax - ymin))
            else:
                val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4),
                       _fmt_count(val))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            x = r.left() + frac * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            val = self._bmin + frac * n * self._bstep
            p.drawText(int(x - 24), int(r.bottom() + 14), f"{val:g}")

        entries = int(self._bins.sum())
        info = f"N={entries:,}"
        if self._under:
            info += f"  under={self._under:,}"
        if self._over:
            info += f"  over={self._over:,}"
        p.setFont(QFont("Monospace", 9))
        p.setPen(QColor(THEME.TEXT_DIM))
        info_rect = QRectF(r.left(), r.top() - 20,
                           max(1.0, self.width() - self._logy_cb.width()
                               - 20 - r.left()),
                           14)
        p.drawText(info_rect,
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   info)

    def mouseMoveEvent(self, ev):
        r = self._plot_rect()
        n = self._bins.size
        if n == 0 or not r.contains(ev.position()):
            if self._hover_idx != -1:
                self._hover_idx = -1
                QToolTip.hideText()
                self.update()
            return
        idx = int((ev.position().x() - r.left()) / r.width() * n)
        idx = max(0, min(n - 1, idx))
        if idx != self._hover_idx:
            self._hover_idx = idx
            self.update()
        lo = self._bmin + idx * self._bstep
        hi = lo + self._bstep
        QToolTip.showText(ev.globalPosition().toPoint(),
                          f"[{lo:g}, {hi:g})  count={int(self._bins[idx]):,}",
                          self)

    def leaveEvent(self, _ev):
        if self._hover_idx != -1:
            self._hover_idx = -1
            QToolTip.hideText()
            self.update()


def _fmt_count(v: float) -> str:
    av = abs(v)
    if av == 0:
        return "0"
    if av >= 1e6:
        return f"{v/1e6:.2f}M"
    if av >= 1e3:
        return f"{v/1e3:.1f}k"
    if av >= 10:
        return f"{v:.0f}"
    return f"{v:.2g}"


# ===========================================================================
#  WaveformPlotWidget — draws the current event's raw FADC samples
# ===========================================================================

class WaveformPlotWidget(QWidget):
    PAD_L, PAD_R, PAD_T, PAD_B = 52, 14, 22, 30
    MAX_STACK = 200

    # Same peak colour palette the web frontend uses (resources/viewer.js PC).
    _PEAK_PALETTE = (
        "#00b4d8", "#ff6b6b", "#51cf66", "#ffd43b",
        "#cc5de8", "#ff922b", "#20c997", "#f06595",
    )

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(150)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self._samples: np.ndarray = np.zeros(0, dtype=np.float32)
        self._peaks: List[Peak] = []
        self._ped_mean: float = 0.0
        self._ped_rms: float = 0.0
        self._title: str = ""
        self._clk_mhz: float = 250.0
        # Active waveform filter (mirrors web monitor's PeakFilter).  Peaks
        # rejected by it are drawn dimmer and the cut regions get a faint
        # overlay rectangle on the plot — same visual language as
        # resources/waveform.js.  ``_filter_show`` mirrors the "show"
        # toggle: when False the overlays are hidden but peaks are still
        # dimmed (so the user can still see which ones are filtered).
        self._filter: Optional[WaveformFilter] = None
        self._filter_show: bool = True

        # --- stack mode state ---
        self._stack_enabled: bool = False
        self._stack_traces: List[np.ndarray] = []   # bounded by MAX_STACK
        self._stack_key: str = ""                   # (module,roc,slot,ch) key

        # --- overlay controls (top-right) ---
        self._stack_cb = QCheckBox("Stack", self)
        self._stack_cb.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._stack_cb.setStyleSheet(_overlay_checkbox_qss())
        self._stack_cb.setToolTip(
            f"Overlay waveforms across events (up to {self.MAX_STACK}). "
            "Peaks and integral shading hidden in stack mode.")
        self._stack_cb.toggled.connect(self._on_stack_toggled)
        self._stack_cb.adjustSize()

        self._stack_clear_btn = QPushButton("Clear", self)
        self._stack_clear_btn.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._stack_clear_btn.setStyleSheet(_overlay_button_qss())
        self._stack_clear_btn.setToolTip("Drop all stacked waveforms")
        self._stack_clear_btn.clicked.connect(self.clear_stack)
        self._stack_clear_btn.setVisible(False)
        self._stack_clear_btn.adjustSize()

        self._stack_count_lbl = QLabel("", self)
        self._stack_count_lbl.setFont(QFont("Monospace", 9))
        self._stack_count_lbl.setStyleSheet(
            f"color:{THEME.TEXT_DIM};background:transparent;")
        self._stack_count_lbl.setVisible(False)
        self._stack_count_lbl.adjustSize()

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def set_data(self, samples: np.ndarray, peaks: List[Peak],
                 ped_mean: float, ped_rms: float,
                 title: str, clk_mhz: float = 250.0,
                 stack_key: Optional[str] = None,
                 filter: Optional["WaveformFilter"] = None,
                 filter_show: bool = True):
        samples = np.asarray(samples, dtype=np.float32)
        self._samples = samples
        self._peaks = peaks
        self._ped_mean = ped_mean
        self._ped_rms  = ped_rms
        self._title = title
        self._clk_mhz = clk_mhz
        self._filter = filter
        self._filter_show = bool(filter_show)

        if self._stack_enabled:
            key = stack_key if stack_key is not None else title
            if key != self._stack_key:
                self._stack_traces = []
                self._stack_key = key
            if samples.size >= 2:
                self._stack_traces.append(samples.copy())
                if len(self._stack_traces) > self.MAX_STACK:
                    self._stack_traces = self._stack_traces[-self.MAX_STACK:]
            self._update_stack_counter()
        self.update()

    def clear(self, title: str = ""):
        self._samples = np.zeros(0, dtype=np.float32)
        self._peaks = []
        self._title = title
        self.update()

    def clear_stack(self):
        """Drop every accumulated trace but keep the current waveform."""
        self._stack_traces = []
        self._stack_key = ""
        self._update_stack_counter()
        self.update()

    def reset_stack_if_new_key(self, key: str):
        """Reset traces when the caller switches to a different channel.

        Lets _display_waveform report a module change even when the new
        module has no samples in the current event — otherwise the empty
        early-return path would leave the previous module's stacks behind.
        """
        if self._stack_enabled and key != self._stack_key:
            self._stack_traces = []
            self._stack_key = key
            self._update_stack_counter()
            self.update()

    def is_stacking(self) -> bool:
        return self._stack_enabled

    # ------------------------------------------------------------------
    #  Internals
    # ------------------------------------------------------------------

    def _on_stack_toggled(self, on: bool):
        self._stack_enabled = on
        self._stack_clear_btn.setVisible(on)
        self._stack_count_lbl.setVisible(on)
        if not on:
            self._stack_traces = []
            self._stack_key = ""
        self._update_stack_counter()
        self._layout_overlays()
        self.update()

    def _update_stack_counter(self):
        self._stack_count_lbl.setText(
            f"{len(self._stack_traces)}/{self.MAX_STACK}")
        self._stack_count_lbl.adjustSize()
        self._layout_overlays()

    def _layout_overlays(self):
        # top-right: [count]  [Clear]  [Stack]
        margin = 6
        x = self.width() - margin
        y = 4
        x -= self._stack_cb.width()
        self._stack_cb.move(x, y)
        if self._stack_clear_btn.isVisible():
            x -= self._stack_clear_btn.width() + 4
            self._stack_clear_btn.move(x, y)
        if self._stack_count_lbl.isVisible():
            x -= self._stack_count_lbl.width() + 6
            self._stack_count_lbl.move(x, y + 2)

    def _plot_rect(self) -> QRectF:
        w, h = self.width(), self.height()
        return QRectF(self.PAD_L, self.PAD_T,
                      max(1.0, w - self.PAD_L - self.PAD_R),
                      max(1.0, h - self.PAD_T - self.PAD_B))

    def resizeEvent(self, ev):
        self._stack_cb.adjustSize()
        self._stack_clear_btn.adjustSize()
        self._stack_count_lbl.adjustSize()
        self._layout_overlays()
        super().resizeEvent(ev)

    # ------------------------------------------------------------------
    #  Painting
    # ------------------------------------------------------------------

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = self._plot_rect()
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        if self._title:
            f = QFont("Monospace", 10); f.setBold(True)
            p.setFont(f)
            p.setPen(QColor(THEME.TEXT))
            suffix = f" — Stacked ({len(self._stack_traces)})" if self._stack_enabled else ""
            p.drawText(int(r.left()), int(r.top() - 6), self._title + suffix)

        if self._stack_enabled:
            self._paint_stacked(p, r)
        else:
            self._paint_single(p, r)

    # --- single-event (default) view ---------------------------------

    def _peak_passes(self, pk) -> bool:
        """Delegate to the active WaveformFilter.passes() so the dimming
        logic stays in lockstep with the histogram / geo filtering."""
        if self._filter is None:
            return True
        return self._filter.passes(pk)

    def _paint_single(self, p: QPainter, r: QRectF):
        n = self._samples.size
        if n < 2:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(no waveform — click Next to load an event)")
            return

        ymin = float(self._samples.min())
        ymax = float(self._samples.max())
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad_y = (ymax - ymin) * 0.05
        ymin -= pad_y; ymax += pad_y

        def to_sx(i: float) -> float:
            return r.left() + (i / (n - 1)) * r.width()

        def to_sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        # ---- cut-region overlays (drawn first, low layer) ---------------
        # Mirrors xRangeShapes / yRangeShapes in resources/waveform.js.  Only
        # drawn when the user has enabled "show" — independent of the
        # "apply" toggle so users can hide the chrome while keeping the
        # filter active (or vice versa).
        ns_total = (n - 1) * 1000.0 / self._clk_mhz if self._clk_mhz > 0 else 0.0
        cut_fill = QColor(THEME.TEXT_MUTED); cut_fill.setAlphaF(0.18)
        cut_edge = QColor(THEME.HIGHLIGHT)
        f = self._filter
        if f is not None and self._filter_show:
            # Time cut (left/right shaded bands outside [time_min, time_max]).
            if ns_total > 0 and f.time_min is not None and f.time_min > 0:
                x1 = to_sx((f.time_min * self._clk_mhz / 1000.0))
                x1 = min(max(x1, r.left()), r.right())
                if x1 > r.left():
                    p.setPen(Qt.PenStyle.NoPen)
                    p.fillRect(QRectF(r.left(), r.top(),
                                      x1 - r.left(), r.height()), cut_fill)
                    p.setPen(QPen(cut_edge, 1, Qt.PenStyle.DashLine))
                    p.drawLine(int(x1), int(r.top()), int(x1), int(r.bottom()))
            if ns_total > 0 and f.time_max is not None and f.time_max < ns_total:
                x2 = to_sx((f.time_max * self._clk_mhz / 1000.0))
                x2 = min(max(x2, r.left()), r.right())
                if x2 < r.right():
                    p.setPen(Qt.PenStyle.NoPen)
                    p.fillRect(QRectF(x2, r.top(),
                                      r.right() - x2, r.height()), cut_fill)
                    p.setPen(QPen(cut_edge, 1, Qt.PenStyle.DashLine))
                    p.drawLine(int(x2), int(r.top()), int(x2), int(r.bottom()))

        # pedestal baseline
        y_ped = to_sy(self._ped_mean) if self._ped_mean != 0 else None
        if y_ped is not None:
            p.setPen(QPen(QColor(THEME.TEXT_DIM), 1, Qt.PenStyle.DashLine))
            p.drawLine(int(r.left()), int(y_ped), int(r.right()), int(y_ped))
            # threshold line (same formula as waveform.js: pm + max(5*pr, 3))
            thr_v = self._ped_mean + max(5.0 * self._ped_rms, 3.0)
            y_thr = to_sy(thr_v)
            p.setPen(QPen(QColor(THEME.TEXT_MUTED), 1, Qt.PenStyle.DotLine))
            p.drawLine(int(r.left()), int(y_thr), int(r.right()), int(y_thr))
            # Height-cut shading — top/bottom bands for samples above
            # height_max + ped or below height_min + ped (filter is on
            # sample-pedestal units; we paint in raw ADC).  Drawn only when
            # "show" is on.
            if f is not None and self._filter_show:
                if f.height_min is not None:
                    hcut_v = self._ped_mean + f.height_min
                    if hcut_v > thr_v:
                        y_hcut = to_sy(hcut_v)
                        y_hcut = min(max(y_hcut, r.top()), r.bottom())
                        if y_hcut > r.top():
                            p.setPen(Qt.PenStyle.NoPen)
                            p.fillRect(QRectF(r.left(), y_hcut,
                                              r.width(), r.bottom() - y_hcut),
                                       cut_fill)
                            p.setPen(QPen(cut_edge, 1, Qt.PenStyle.DashLine))
                            p.drawLine(int(r.left()), int(y_hcut),
                                       int(r.right()), int(y_hcut))
                if f.height_max is not None:
                    hcut_v = self._ped_mean + f.height_max
                    y_hcut = to_sy(hcut_v)
                    y_hcut = min(max(y_hcut, r.top()), r.bottom())
                    if y_hcut < r.bottom():
                        p.setPen(Qt.PenStyle.NoPen)
                        p.fillRect(QRectF(r.left(), r.top(),
                                          r.width(), y_hcut - r.top()),
                                   cut_fill)
                        p.setPen(QPen(cut_edge, 1, Qt.PenStyle.DashLine))
                        p.drawLine(int(r.left()), int(y_hcut),
                                   int(r.right()), int(y_hcut))

        # Fill the integral area (between pedestal and waveform) per peak,
        # colour-coded from _PEAK_PALETTE. Mirrors resources/waveform.js.
        # Peaks rejected by the active filter are drawn with reduced alpha
        # so the user can tell at a glance which peaks the geo / hists use.
        if self._peaks and y_ped is not None:
            for i, pk in enumerate(self._peaks):
                passes = self._peak_passes(pk)
                base = QColor(self._PEAK_PALETTE[i % len(self._PEAK_PALETTE)])
                fill = QColor(base)
                fill.setAlphaF(0.18 if passes else 0.06)
                poly = QPolygonF()
                j = max(0, int(pk.left))
                j_end = min(n - 1, int(pk.right))
                for k in range(j, j_end + 1):
                    poly.append(QPointF(to_sx(k),
                                        to_sy(float(self._samples[k]))))
                # close along the pedestal baseline
                poly.append(QPointF(to_sx(j_end), y_ped))
                poly.append(QPointF(to_sx(j), y_ped))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(fill)
                p.drawPolygon(poly)
                # outline the peak section with the solid palette colour
                outline = QColor(base)
                if not passes:
                    outline.setAlphaF(0.45)
                p.setPen(QPen(outline, 2))
                p.setBrush(Qt.BrushStyle.NoBrush)
                for k in range(j, j_end):
                    p.drawLine(QPointF(to_sx(k),
                                       to_sy(float(self._samples[k]))),
                               QPointF(to_sx(k + 1),
                                       to_sy(float(self._samples[k + 1]))))

        # waveform line (default accent, drawn under peak outlines)
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        for i in range(n - 1):
            p.drawLine(int(to_sx(i)),     int(to_sy(float(self._samples[i]))),
                       int(to_sx(i + 1)), int(to_sy(float(self._samples[i + 1]))))

        # peak markers (diamonds, coloured per peak; rejected peaks get a
        # hollow diamond to keep them readable but visually distinct).
        if self._peaks:
            for i, pk in enumerate(self._peaks):
                if pk.pos < 0 or pk.pos >= n:
                    continue
                passes = self._peak_passes(pk)
                col = QColor(self._PEAK_PALETTE[i % len(self._PEAK_PALETTE)])
                if not passes:
                    col.setAlphaF(0.55)
                p.setPen(QPen(col, 1.2))
                p.setBrush(col if passes else Qt.BrushStyle.NoBrush)
                cx = to_sx(pk.pos)
                cy = to_sy(float(self._samples[pk.pos]))
                diamond = QPolygonF([
                    QPointF(cx,     cy - 4),
                    QPointF(cx + 4, cy),
                    QPointF(cx,     cy + 4),
                    QPointF(cx - 4, cy),
                ])
                p.drawPolygon(diamond)

        self._paint_axes(p, r, ymin, ymax)

        # ped/rms/peak-count readout — drawn inside the plot at top-right to
        # stay clear of the Stack checkbox / Clear button in the widget's
        # top-right margin.
        info = (f"ped={self._ped_mean:.1f}  rms={self._ped_rms:.2f}  "
                f"peaks={len(self._peaks)}")
        p.setFont(QFont("Monospace", 9))
        fm = p.fontMetrics()
        tw = fm.horizontalAdvance(info)
        th = fm.height()
        pad = 4
        box = QRectF(r.right() - tw - 2 * pad - 2, r.top() + 4,
                     tw + 2 * pad, th + 2)
        bg = QColor(THEME.BG); bg.setAlphaF(0.70)
        p.fillRect(box, bg)
        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(box,
                   Qt.AlignmentFlag.AlignCenter, info)

    # --- stacked overlay view -----------------------------------------

    def _paint_stacked(self, p: QPainter, r: QRectF):
        traces = self._stack_traces
        if not traces:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(stack is empty — step through events to accumulate)")
            return

        # Compute common y-range across all traces.
        ymin = min(float(w.min()) for w in traces)
        ymax = max(float(w.max()) for w in traces)
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad_y = (ymax - ymin) * 0.05
        ymin -= pad_y; ymax += pad_y

        # Width uses the max length so shorter traces still fit left-aligned.
        n_max = max(w.size for w in traces)

        def to_sx(i: float, n: int) -> float:
            return r.left() + (i / max(1, n - 1)) * r.width()

        def to_sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        # Dimmed stacked traces.
        dim = QColor(THEME.ACCENT); dim.setAlphaF(0.18)
        p.setPen(QPen(dim, 1))
        for w in traces[:-1]:
            n = w.size
            for i in range(n - 1):
                p.drawLine(int(to_sx(i, n)),     int(to_sy(float(w[i]))),
                           int(to_sx(i + 1, n)), int(to_sy(float(w[i + 1]))))

        # Latest trace drawn on top at full colour.
        latest = traces[-1]
        n = latest.size
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        for i in range(n - 1):
            p.drawLine(int(to_sx(i, n)),     int(to_sy(float(latest[i]))),
                       int(to_sx(i + 1, n)), int(to_sy(float(latest[i + 1]))))

        self._paint_axes(p, r, ymin, ymax, n=n_max)

        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(QRectF(r.left(), r.top() - 20,
                          max(1.0, r.width() - 8), 14),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   f"stack={len(traces)}/{self.MAX_STACK}")

    # --- shared axis/tick drawing -------------------------------------

    def _paint_axes(self, p: QPainter, r: QRectF,
                    ymin: float, ymax: float, n: Optional[int] = None):
        if n is None:
            n = self._samples.size
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4), f"{val:.0f}")
        tick_every = max(1, n // 8)
        for i in range(0, n, tick_every):
            x = r.left() + (i / max(1, n - 1)) * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            ns = i * 1e3 / self._clk_mhz
            p.drawText(int(x - 18), int(r.bottom() + 14), f"{ns:g}")

        p.setFont(QFont("Monospace", 9))
        p.drawText(int(r.left() + r.width() / 2 - 10),
                   int(r.bottom() + 26), "ns")


# ===========================================================================
#  WaveformGeoView — small HyCal overview for module selection
# ===========================================================================

# Module types that should get a small name label painted on top of the cell
# (so the tiny LMS / V blocks off to the left of HyCal are identifiable).
_LABEL_TYPES = {"LMS", "Veto"}

class WaveformGeoView(HyCalMapWidget):
    """Compact HyCal geo view with two colour-coding modes.

    * ``current``  — module colour = max peak integral in the current event.
    * ``overall``  — module colour = occupancy (events-with-peak / accumulated
      events) across all events the user has browsed / batched.

    Modules that have never been seen in an event are drawn in a flat grey
    so "no-data" stays visually distinct from "data, low value".  Clicking
    any module emits moduleClicked with its name.
    """

    MODE_CURRENT = "current"
    MODE_OVERALL = "overall"

    # Resolved at paint time so the active theme wins; see :class:`THEME`.
    @property
    def UNAVAIL_COLOR(self) -> QColor:
        return QColor(THEME.BORDER)

    @property
    def SELECT_COLOR(self) -> QColor:
        return QColor(THEME.SELECT_BORDER)

    def __init__(self, parent=None):
        # margin_bottom must exceed the base's colour-bar anchor (cb_y =
        # h - 40) so the module rects clear the bar — leave ~16 px gap.
        super().__init__(parent, show_colorbar=True, include_lms=True,
                         margin_top=4, margin_bottom=56,
                         min_size=(220, 280), shrink=0.90)
        self._available: set = set()
        self._selected_name: Optional[str] = None
        self._label_names: set = set()          # filled in set_modules()
        self._mode = self.MODE_CURRENT
        self._current_vals: Dict[str, float] = {}
        self._overall_vals: Dict[str, float] = {}
        # Modules where the analyser found peaks in the current event but
        # all of them were rejected by the threshold / time-window cut —
        # i.e. they appear as shaded peaks on the waveform plot but don't
        # contribute to the geo's max-integral.  Used by _tooltip_text.
        self._rejected_current: set = set()
        # Headless range controller: handles auto-fit logic.  Both vmin and
        # vmax are inline-editable on the colorbar so the user can pin the
        # palette to a fixed range when comparing events.
        self._range_ctrl = ColorRangeController(
            self, auto_fit="minmax", parent=self)

        # Top-left mode toggle.  Default label matches MODE_CURRENT.
        self._mode_btn = QPushButton("Current", self)
        self._mode_btn.setFixedSize(74, 22)
        _f = QFont("Consolas", 9); _f.setBold(True)
        self._mode_btn.setFont(_f)
        self._mode_btn.setToolTip(
            "Colour coding:\n"
            "  Current — max peak integral in the currently viewed event\n"
            "  Overall — occupancy (events-with-peak / accumulated events)")
        self._mode_btn.setStyleSheet(themed(
            "QPushButton{background:rgba(29,29,31,220);color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:4px;}"
            "QPushButton:hover{background:#28282a;color:#e6edf3;}"))
        self._mode_btn.clicked.connect(self._toggle_mode)

    def set_modules(self, modules):
        super().set_modules(modules)
        self._label_names = {m.name for m in self._modules
                             if m.mod_type in _LABEL_TYPES}

    def set_available(self, names):
        self._available = set(names)
        self.update()

    def set_selected_module(self, name: Optional[str]):
        if name != self._selected_name:
            self._selected_name = name
            self.update()

    def set_current_values(self, vals: Dict[str, float]):
        self._current_vals = vals
        if self._mode == self.MODE_CURRENT:
            self._apply_mode_values()

    def set_rejected_current(self, names: set):
        """Modules where the current event had peaks rejected by the cut.
        Pass an empty set to clear.  Drives the tooltip wording so users can
        tell "no peak" apart from "peak rejected by filter"."""
        self._rejected_current = set(names) if names else set()

    def set_overall_values(self, vals: Dict[str, float]):
        self._overall_vals = vals
        if self._mode == self.MODE_OVERALL:
            self._apply_mode_values()

    def _toggle_mode(self):
        self._mode = (self.MODE_OVERALL if self._mode == self.MODE_CURRENT
                      else self.MODE_CURRENT)
        self._mode_btn.setText(
            "Overall" if self._mode == self.MODE_OVERALL else "Current")
        self._apply_mode_values()

    def _apply_mode_values(self):
        if self._mode == self.MODE_OVERALL:
            self.set_values(self._overall_vals)
            self.set_range(0.0, 1.0)          # occupancy fraction
        else:
            self.set_values(self._current_vals)
            # Re-fit per event; user can override via inline colorbar edit.
            self._range_ctrl.auto_fit(self._current_vals)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._mode_btn.move(6, 6)

    def _colorbar_center_text(self) -> str:
        return ("occupancy" if self._mode == self.MODE_OVERALL
                else "max peak integral")

    def _paint_modules(self, p):
        # Keep the "not yet seen" grey distinct from the colormap's low end
        # so unused modules don't masquerade as "low value".
        avail = self._available
        u_col = self.UNAVAIL_COLOR
        no_data = self.NO_DATA_COLOR
        stops = self.palette_stops()
        vmin, vmax = self._vmin, self._vmax
        vals = self._values
        for name, rect in self._rects.items():
            if name not in avail:
                p.fillRect(rect, u_col)
                continue
            v = vals.get(name)
            if v is None:
                p.fillRect(rect, no_data)
                continue
            t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
            p.fillRect(rect, cmap_qcolor(t, stops))

    def _paint_overlays(self, p, w, h):
        p.setPen(QColor(THEME.TEXT))
        p.setFont(QFont("Monospace", 7, QFont.Weight.Bold))
        for name in self._label_names:
            r = self._rects.get(name)
            if r is not None:
                p.drawText(r, Qt.AlignmentFlag.AlignCenter, name)

        if self._selected_name and self._selected_name in self._rects:
            p.setPen(QPen(self.SELECT_COLOR, 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._selected_name])
        super()._paint_overlays(p, w, h)   # hover border

    def _tooltip_text(self, name: str) -> str:
        if name not in self._available:
            return f"{name}  (not seen yet)"
        v = self._values.get(name)
        unit = ("occupancy" if self._mode == self.MODE_OVERALL
                else "max integral")
        if v is None:
            # In Current mode, distinguish "no peak found" from "peak found
            # but rejected by the threshold/time cut" — the latter is what
            # confused the user when the waveform plot shows shaded peaks
            # but the geo / tooltip shows nothing.
            if (self._mode == self.MODE_CURRENT
                    and name in self._rejected_current):
                return f"{name}  ({unit}: — peak outside cut)"
            return f"{name}  ({unit}: —)"
        return f"{name}  {unit}={v:.3g}"


# ===========================================================================
#  Cluster display wrappers
# ===========================================================================


class _DisplayCluster:
    """Display-friendly wrapper around a bound ClusterHit.  We can't set
    arbitrary attributes on pybind11 objects, so carry centre-name +
    member-name lists alongside the raw hit."""
    __slots__ = ("energy", "x", "y", "nblocks", "npos",
                 "center_id", "center_name", "members")

    def __init__(self, hit, center_name: str, members: List[str]):
        self.energy      = float(hit.energy)
        self.x           = float(hit.x)
        self.y           = float(hit.y)
        self.nblocks     = int(hit.nblocks)
        self.npos        = int(hit.npos)
        self.center_id   = int(hit.center_id)
        self.center_name = center_name or ""
        self.members     = list(members)


# ===========================================================================
#  Cluster map widget — HyCal heatmap + cluster overlays
# ===========================================================================


# Per-cluster frame/chip colours (mirrors ``PC`` in resources/cluster.js so
# the Python viewer and the web monitor agree on "cluster #3 is green").
CLUSTER_PALETTE: Tuple[str, ...] = (
    "#00b4d8", "#ff6b6b", "#51cf66", "#ffd43b",
    "#cc5de8", "#ff922b", "#20c997", "#f06595",
)


def cluster_color(idx: int) -> QColor:
    return QColor(CLUSTER_PALETTE[idx % len(CLUSTER_PALETTE)])


class HyCalClusterMap(HyCalMapWidget):
    """HyCal map coloured by per-module energy (MeV) with cluster overlays
    (crosshair + small circle + energy label) mirroring the web monitor."""

    def __init__(self, parent=None):
        super().__init__(parent, include_lms=True, enable_zoom_pan=True,
                         show_colorbar=True)
        self._clusters: List = []               # List[ClusterHit]
        self._selected_cluster: Optional[int] = None
        self._member_modules: set = set()       # module names of selected cluster

    # -- public API ------------------------------------------------------

    def set_clusters(self, clusters):
        self._clusters = list(clusters) if clusters else []
        self._recompute_membership()
        self.update()

    def set_selected_cluster(self, idx: Optional[int]):
        """Highlight one cluster; pass None for 'show all'."""
        if idx is not None and not (0 <= idx < len(self._clusters)):
            idx = None
        self._selected_cluster = idx
        self._recompute_membership()
        self.update()

    # -- internals -------------------------------------------------------

    def _recompute_membership(self):
        """Build set of module names that belong to the selected cluster.
        ClusterHit doesn't carry a module list; we use ``center_id`` plus
        the ``members`` list if the caller attached one.  If not, only the
        centre module lights up for a selected cluster — still useful."""
        self._member_modules.clear()
        if self._selected_cluster is None:
            return
        cl = self._clusters[self._selected_cluster]
        members = getattr(cl, "members", None)
        if members:
            self._member_modules.update(members)
        elif getattr(cl, "center_name", ""):
            self._member_modules.add(cl.center_name)

    def _paint_modules(self, p):
        """Default colormap paint, but dim non-members when a cluster is
        selected.  Calls the base implementation via a per-module alpha."""
        dim = self._selected_cluster is not None and self._member_modules
        if not dim:
            super()._paint_modules(p)
            return
        stops = self.palette_stops()
        vmin, vmax = self._vmin, self._vmax
        no_data = self.NO_DATA_COLOR
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is None:
                col = QColor(no_data)
            else:
                t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                col = cmap_qcolor(t, stops)
            if name not in self._member_modules:
                col = QColor(col.red(), col.green(), col.blue(), 60)
            p.fillRect(rect, col)

    def _paint_cluster_frames(self, p):
        """Draw a coloured border around every member module of every
        cluster (or just the selected cluster).  Mirrors the per-module
        border colouring in resources/cluster.js."""
        if not self._clusters:
            return
        sel = self._selected_cluster
        p.save()
        p.setBrush(Qt.BrushStyle.NoBrush)
        for i, cl in enumerate(self._clusters):
            if sel is not None and i != sel:
                continue
            members = getattr(cl, "members", None) or ()
            if not members:
                continue
            width = 2.5 if (sel is not None and i == sel) else 1.5
            p.setPen(QPen(cluster_color(i), width))
            for name in members:
                rect = self._rects.get(name)
                if rect is not None:
                    p.drawRect(rect)
        p.restore()

    def _paint_overlays(self, p, w, h):
        # Per-cluster coloured frames first, so the hover border (drawn by
        # super) and the cluster crosshairs/labels stay visible on top.
        self._paint_cluster_frames(p)
        super()._paint_overlays(p, w, h)   # hover border
        if not self._clusters:
            return

        p.save()
        cross_pen = QPen(QColor("#ffd166"), 1.6)
        circle_pen = QPen(QColor("#ffd166"), 1.4)
        label_pen = QPen(QColor("#ffffff" if THEME.BG.startswith("#0")
                                or THEME.BG == "#000000" else "#1d1d1f"))
        font = QFont("Monospace", 9, QFont.Weight.Bold)
        p.setFont(font)
        fm = p.fontMetrics()

        for i, cl in enumerate(self._clusters):
            if self._selected_cluster is not None and i != self._selected_cluster:
                continue
            pt = self.geo_to_canvas(cl.x, cl.y)
            # crosshair
            p.setPen(cross_pen)
            L = 9
            p.drawLine(QPointF(pt.x() - L, pt.y()),
                       QPointF(pt.x() + L, pt.y()))
            p.drawLine(QPointF(pt.x(), pt.y() - L),
                       QPointF(pt.x(), pt.y() + L))
            # circle scaled by log-ish energy
            r = max(6.0, min(20.0, 4.0 + 2.0 * (cl.energy ** 0.33)))
            p.setPen(circle_pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawEllipse(pt, r, r)
            # label: "C{id}  {E:.0f} MeV"
            text = f"C{i}  {cl.energy:.0f}"
            tw = fm.horizontalAdvance(text)
            tx = pt.x() + r + 3
            ty = pt.y() - 3
            # background halo for readability
            halo = QColor(0, 0, 0, 140)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(halo)
            p.drawRect(QRectF(tx - 2, ty - fm.ascent(),
                              tw + 4, fm.height()))
            p.setPen(QColor("#ffd166"))
            p.drawText(QPointF(tx, ty), text)
        p.restore()


# ===========================================================================
#  Cluster panel — selector + table + footer stats
# ===========================================================================


class ClusterPanel(QWidget):
    """Right-side pane for the Cluster tab: combo + table + summary line."""

    clusterSelected = pyqtSignal(object)    # int | None  (None = show all)

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(4)

        self._selector = QComboBox()
        self._selector.addItem("All clusters", None)
        self._selector.currentIndexChanged.connect(self._on_selector_changed)
        lay.addWidget(self._selector)

        self._table = QTableWidget(0, 6)
        self._table.setHorizontalHeaderLabels(
            ["#", "center", "E [MeV]", "x [mm]", "y [mm]", "nblocks"])
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        hh.setStretchLastSection(True)
        self._table.itemSelectionChanged.connect(self._on_table_selection)
        lay.addWidget(self._table, stretch=1)

        self._summary = QLabel("no event")
        self._summary.setFont(QFont("Monospace", 9))
        lay.addWidget(self._summary)

        self._clusters: List = []

    def set_clusters(self, clusters):
        self._clusters = list(clusters) if clusters else []
        # Repopulate selector
        self._selector.blockSignals(True)
        self._selector.clear()
        self._selector.addItem("All clusters", None)
        for i, cl in enumerate(self._clusters):
            label = f"C{i}  E={cl.energy:.0f} MeV  {getattr(cl, 'center_name', '?')}"
            self._selector.addItem(label, i)
        self._selector.setCurrentIndex(0)
        self._selector.blockSignals(False)

        # Populate table
        self._table.blockSignals(True)
        self._table.setRowCount(len(self._clusters))
        chip_font = QFont("Monospace", 9, QFont.Weight.Bold)
        text_black = QBrush(QColor("#000000"))
        for i, cl in enumerate(self._clusters):
            row = [
                f"{i}",
                getattr(cl, "center_name", "?"),
                f"{cl.energy:.1f}",
                f"{cl.x:.1f}",
                f"{cl.y:.1f}",
                f"{cl.nblocks}",
            ]
            for c, v in enumerate(row):
                item = QTableWidgetItem(v)
                if c == 0:
                    # Colour chip linking the row to its cluster colour.
                    item.setBackground(QBrush(cluster_color(i)))
                    item.setForeground(text_black)
                    item.setFont(chip_font)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                elif c in (2, 3, 4, 5):
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                self._table.setItem(i, c, item)
        self._table.blockSignals(False)

        n = len(self._clusters)
        tot = sum(c.energy for c in self._clusters)
        self._summary.setText(f"{n} clusters   ΣE = {tot:.0f} MeV")

    def _on_selector_changed(self, _idx: int):
        data = self._selector.currentData()
        # Sync table selection without echoing
        self._table.blockSignals(True)
        if data is None:
            self._table.clearSelection()
        else:
            self._table.selectRow(int(data))
        self._table.blockSignals(False)
        self.clusterSelected.emit(data)

    def _on_table_selection(self):
        rows = self._table.selectionModel().selectedRows()
        if not rows:
            return
        i = rows[0].row()
        # Sync combo
        self._selector.blockSignals(True)
        self._selector.setCurrentIndex(i + 1)   # +1 for leading "All"
        self._selector.blockSignals(False)
        self.clusterSelected.emit(i)


# ===========================================================================
#  Cut Settings dialog — modal editor for the active WaveformFilter.
#  Mirrors resources/cut_dialog.js + viewer.html's #cut-dialog markup so
#  users moving between the web monitor and this viewer see the same
#  fields with the same semantics.
# ===========================================================================


class _RangeRow(QWidget):
    """Two QLineEdits (min/max) with a validator.  Empty value = no
    constraint, matching the web monitor's <input type='number'> + null
    parsing in cut_dialog.js."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)
        self.min_edit = QLineEdit()
        self.max_edit = QLineEdit()
        for w, ph in ((self.min_edit, "min"), (self.max_edit, "max")):
            w.setPlaceholderText(ph)
            w.setValidator(QDoubleValidator(self))
            w.setMinimumWidth(80)
        lay.addWidget(QLabel("min"))
        lay.addWidget(self.min_edit, 1)
        lay.addSpacing(6)
        lay.addWidget(QLabel("max"))
        lay.addWidget(self.max_edit, 1)

    def set_values(self, lo: Optional[float], hi: Optional[float]):
        self.min_edit.setText("" if lo is None else _fmt_filter_num(lo))
        self.max_edit.setText("" if hi is None else _fmt_filter_num(hi))

    def values(self) -> Tuple[Optional[float], Optional[float]]:
        def _parse(s: str) -> Optional[float]:
            s = s.strip()
            if not s:
                return None
            try:
                return float(s)
            except ValueError:
                return None
        return _parse(self.min_edit.text()), _parse(self.max_edit.text())


def _fmt_filter_num(v: float) -> str:
    """Compact representation that round-trips through float() — mirrors
    JSON.stringify(num) on the web side."""
    if v == int(v):
        return f"{int(v)}"
    return f"{v:g}"


class CutSettingsDialog(QDialog):
    """Modal "Cut Settings" editor — replicates resources/cut_dialog.js.

    The user edits the active filter's per-axis ranges and quality bits,
    then clicks Save to commit.  Cancel discards changes.  Reset reverts
    the form (without committing) to the defaults snapshotted at startup
    from monitor_config.json's ``waveform.filter`` block.
    """

    def __init__(self, parent: QWidget,
                 current: WaveformFilter,
                 default: WaveformFilter):
        super().__init__(parent)
        self.setWindowTitle("Cut Settings")
        self.setModal(True)
        self._default = default

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(8)

        # Time / Integral / Height range groups (one QGroupBox each, mirroring
        # the web's <fieldset class="cut-fieldset"> blocks).
        self._rows: Dict[str, _RangeRow] = {}
        for axis, label, suffix in (("time",     "Time",     " (ns)"),
                                    ("integral", "Integral", ""),
                                    ("height",   "Height",   "")):
            gb = QGroupBox(label + suffix)
            gl = QVBoxLayout(gb)
            gl.setContentsMargins(8, 4, 8, 6)
            row = _RangeRow()
            self._rows[axis] = row
            gl.addWidget(row)
            root.addWidget(gb)

        # Quality bits — two columns of checkboxes (accept / reject).
        # Mirrors buildBitList() in cut_dialog.js.
        qg = QGroupBox("Quality bits")
        qg_lay = QGridLayout(qg)
        qg_lay.setContentsMargins(8, 4, 8, 6)
        qg_lay.setHorizontalSpacing(20)
        accept_lbl = QLabel("Accept")
        reject_lbl = QLabel("Reject")
        for lbl in (accept_lbl, reject_lbl):
            f = QFont(); f.setBold(True)
            lbl.setFont(f)
        qg_lay.addWidget(accept_lbl, 0, 0)
        qg_lay.addWidget(reject_lbl, 0, 1)

        self._accept_checks: Dict[str, QCheckBox] = {}
        self._reject_checks: Dict[str, QCheckBox] = {}
        if not PEAK_QUALITY_BITS:
            note = QLabel("(no quality bits exposed by prad2py)")
            note.setStyleSheet(f"color:{THEME.TEXT_DIM}; font-style:italic;")
            qg_lay.addWidget(note, 1, 0, 1, 2)
        else:
            for r, b in enumerate(PEAK_QUALITY_BITS, start=1):
                acc = QCheckBox(str(b["label"]))
                rej = QCheckBox(str(b["label"]))
                self._accept_checks[str(b["name"])] = acc
                self._reject_checks[str(b["name"])] = rej
                acc.toggled.connect(self._sync_bit_mutex)
                rej.toggled.connect(self._sync_bit_mutex)
                qg_lay.addWidget(acc, r, 0)
                qg_lay.addWidget(rej, r, 1)
        hint = QLabel(
            "Empty = no constraint.  Accept: peak's set bits must overlap "
            "the accepted flags.  Reject: peak fails if any rejected flag is set.")
        hint.setStyleSheet(f"color:{THEME.TEXT_DIM}; font-size:10px;")
        hint.setWordWrap(True)
        qg_lay.addWidget(hint, qg_lay.rowCount(), 0, 1, 2)
        root.addWidget(qg)

        # Standard Cancel / Reset / Save buttonbox.
        bb = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Reset)
        bb.button(QDialogButtonBox.StandardButton.Save).setDefault(True)
        bb.button(QDialogButtonBox.StandardButton.Reset).setToolTip(
            "Restore the form to monitor_config.json's filter values "
            "(does not commit until you click Save).")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(
            self._on_reset)
        root.addWidget(bb)

        self._populate(current)

    # ------------------------------------------------------------------

    def _populate(self, f: WaveformFilter):
        self._rows["time"]    .set_values(f.time_min,     f.time_max)
        self._rows["integral"].set_values(f.integral_min, f.integral_max)
        self._rows["height"]  .set_values(f.height_min,   f.height_max)
        for d in PEAK_QUALITY_BITS:
            name = str(d["name"]); mask = int(d["mask"])
            acc = self._accept_checks.get(name)
            rej = self._reject_checks.get(name)
            if acc is not None: acc.setChecked(bool(f.q_accept & mask))
            if rej is not None: rej.setChecked(bool(f.q_reject & mask))
        self._sync_bit_mutex()

    def _on_reset(self):
        self._populate(self._default)

    def _sync_bit_mutex(self):
        """A bit can be in Accept OR Reject but not both — disable the
        twin checkbox when its sibling is checked.  Mirrors
        syncBitMutualExclusion() in cut_dialog.js."""
        for name, acc in self._accept_checks.items():
            rej = self._reject_checks.get(name)
            if rej is None:
                continue
            acc.setEnabled(not rej.isChecked())
            rej.setEnabled(not acc.isChecked())

    def result_filter(self) -> WaveformFilter:
        """Build a WaveformFilter from the current form contents.  ``enable``
        is preserved from the input filter (controlled by the toolbar's
        "apply" toggle, not this dialog)."""
        f = WaveformFilter()
        f.time_min,     f.time_max     = self._rows["time"].values()
        f.integral_min, f.integral_max = self._rows["integral"].values()
        f.height_min,   f.height_max   = self._rows["height"].values()
        f.q_accept = _names_to_mask([n for n, cb in self._accept_checks.items()
                                     if cb.isChecked()])
        f.q_reject = _names_to_mask([n for n, cb in self._reject_checks.items()
                                     if cb.isChecked()])
        return f


# ===========================================================================
#  Main window
# ===========================================================================

_NATKEY_RE = re.compile(r"(\d+)")
def _natural_sort_key(s: str):
    return [int(p) if p.isdigit() else p.lower()
            for p in _NATKEY_RE.split(s or "")]


class HyCalEventViewer(QMainWindow):

    def __init__(self,
                 *,
                 hist_config: Dict,
                 daq_map: Dict,
                 roc_to_crate: Dict,
                 accept_mask: int,
                 reject_mask: int,
                 daq_config_path: str,
                 hycal_modules: Optional[List] = None,
                 hycal_map_path: Optional[str] = None,
                 recon_config_path: Optional[str] = None):
        super().__init__()
        self._hist_config   = hist_config
        self._daq_map       = daq_map
        self._roc_to_crate  = roc_to_crate
        self._accept_mask   = accept_mask
        self._reject_mask   = reject_mask
        self._daq_cfg_path  = daq_config_path
        self._hycal_modules = hycal_modules or []
        self._hycal_map_path = hycal_map_path
        self._recon_cfg_path = recon_config_path

        # Bin configs — merge user config with defaults, add n-peaks hist
        self._h_cfg = hist_config.get("height_hist",
                                      {"min": 0, "max": 4000,  "step": 10})
        self._i_cfg = hist_config.get("integral_hist",
                                      {"min": 0, "max": 20000, "step": 100})
        self._p_cfg = hist_config.get("time_hist",
                                      {"min": 0, "max": 400,   "step": 4})
        # Left edge at -0.5 so integer n_peaks values (0, 1, 2 …) sit on bin
        # centres rather than at the left edge of each bar.
        self._n_cfg = {"min": -0.5, "max": 10.5, "step": 1}

        # Seed the analyzer config from daq_config.json's
        # `fadc250_waveform.analyzer` block.  This is the single source of
        # truth for peak-detection knobs (peak_nsigma, min_peak_height,
        # min_peak_ratio); monitor_config.json no longer overrides them.
        # Falls back to plain defaults if the daq_config can't be loaded —
        # opening files later will fail loudly anyway in that case.
        try:
            _dc = (prad2py.dec.load_daq_config(self._daq_cfg_path)
                   if self._daq_cfg_path
                   else prad2py.dec.load_daq_config())
            self._wcfg = WaveConfig(_dc.wave_cfg)
        except Exception:
            self._wcfg = WaveConfig()

        # Waveform peak filter — same JSON shape and semantics as the web
        # monitor's `waveform.filter` (see resources/cut_dialog.js + C++
        # PeakFilter).  Defaults parsed from monitor_config.json; the
        # snapshot in `_filter_default` lets the Cut-Settings "Reset" button
        # restore the file values without a server round-trip.  ``enable``
        # is the GUI "apply" toggle (filter is a no-op when False); ``show``
        # is the client-side overlay toggle on the waveform plot.
        _flt_json = hist_config.get("filter") or {}
        self._filter         = WaveformFilter.from_json(_flt_json)
        self._filter_default = WaveformFilter.from_json(_flt_json)
        self._filter_show    = True

        # Debounce timer for Advanced-dock re-runs: coalesce rapid
        # slider drags into a single re-read + re-analyse.  Must exist
        # before the dock widgets are built (they connect to it via
        # _on_advanced_changed).
        self._adv_debounce_ms = 150
        self._adv_redraw_timer = QTimer(self)
        self._adv_redraw_timer.setSingleShot(True)
        self._adv_redraw_timer.timeout.connect(self._rerun_current_event)

        # File state
        self._evio_path: Optional[Path] = None
        self._index: List[Tuple[int, int]] = []
        self._current_idx: int = -1
        # One bool per physics sub-event — True once its peaks have been
        # folded into self._channels hists.  Lets Prev/Next re-display an
        # event without double-counting.  np.bool = 1 byte/event, so 1 M
        # events ≈ 1 MB, 10 M ≈ 10 MB — negligible.
        self._accumulated: Optional[np.ndarray] = None

        # Per-channel accumulated hists, keyed by (roc, slot, ch)
        self._channels: Dict[Tuple[int, int, int], ChannelHists] = {}
        self._selected_key: Optional[Tuple[int, int, int]] = None

        # Reader state — open EvChannel in RA mode, kept alive across browse
        self._ch: Optional["prad2py.dec.EvChannel"] = None
        self._reader_path: Optional[str] = None
        self._reader_is_ra: bool = False
        self._reader_pos: int = -1

        # Worker threads
        self._idx_worker: Optional[IndexerWorker] = None
        self._idx_thread: Optional[QThread] = None
        self._batch_worker: Optional[BatchWorker] = None
        self._batch_thread: Optional[QThread] = None

        # HyCal clustering — wired up via prad2py.det.PipelineBuilder so the
        # daq / recon / runinfo configs and per-run HyCal calibration are all
        # loaded automatically (matches the live server / analysis scripts).
        # `_pipeline` keeps the C++ Pipeline object alive; the borrowed
        # `hycal` reference is stored separately as `_hcsys` for the existing
        # call sites (`module_by_daq`, `load_calibration`, …).  Per-channel
        # DAQ → module lookups are cached by (crate, slot, ch) so the
        # per-event hot loop doesn't go through pybind11 every hit.
        self._pipeline = None
        self._hcsys = None
        self._hccl  = None
        self._hc_cache: Dict[Tuple[int, int, int], object] = {}
        # True once a calibration file has actually been loaded (either via
        # PipelineBuilder's runinfo lookup or the manual menu override).
        # Gates the "no calibration" warning banner on the cluster tab.
        self._hycal_calib_loaded = False
        if not self._build_hycal_pipeline():
            self._hycal_init_fallback()

        apply_theme_palette(self)
        self._build_ui()
        self._make_menu()

    # -- UI --

    def _build_ui(self):
        self.setWindowTitle("HyCal Event Viewer")
        self.resize(1500, 1000)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(3)

        self._file_lbl = QLabel("(no file loaded)")
        self._file_lbl.setFont(QFont("Monospace", 10))
        self._file_lbl.setStyleSheet(themed("color:#8b949e;"))
        root.addWidget(self._file_lbl)

        # -- top control bar: navigation on the left, module picker on the right --
        top = QHBoxLayout()

        self._prev_btn = self._small_btn("◀ Prev", self._on_prev)
        self._next_btn = self._small_btn("Next ▶", self._on_next)
        self._prev_btn.setEnabled(False)
        self._next_btn.setEnabled(False)
        top.addWidget(self._prev_btn)
        top.addWidget(self._next_btn)

        top.addSpacing(12)
        top.addWidget(self._mk_label("Event:"))
        self._event_spin = QSpinBox()
        self._event_spin.setFont(QFont("Monospace", 10))
        self._event_spin.setMinimum(0)
        self._event_spin.setMaximum(0)
        self._event_spin.setStyleSheet(themed(
            "QSpinBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:6px;padding:2px 6px;}"))
        self._event_spin.editingFinished.connect(self._on_spin_jump)
        self._event_spin.setEnabled(False)
        top.addWidget(self._event_spin)
        self._total_lbl = QLabel(" / 0")
        self._total_lbl.setFont(QFont("Monospace", 10))
        self._total_lbl.setStyleSheet(themed("color:#8b949e;"))
        top.addWidget(self._total_lbl)

        top.addSpacing(18)
        self._batch_btn = self._small_btn("Process next 10k",
                                          self._on_batch_10k, primary=True)
        self._batch_btn.setEnabled(False)
        top.addWidget(self._batch_btn)

        self._batch_status = QLabel("")
        self._batch_status.setFont(QFont("Monospace", 10))
        self._batch_status.setStyleSheet(themed("color:#8b949e;"))
        top.addSpacing(8)
        top.addWidget(self._batch_status)

        top.addStretch(1)

        # Right cluster: module dropdown + reset hist.
        mod_lbl = QLabel("Module:")
        mod_lbl.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
        mod_lbl.setStyleSheet(themed("color:#c9d1d9;"))
        top.addWidget(mod_lbl)
        self._combo = QComboBox()
        self._combo.setEditable(True)
        self._combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._combo.setFont(QFont("Monospace", 11))
        self._combo.setMinimumContentsLength(32)
        self._combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:6px;padding:2px 6px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "selection-background-color:#1f6feb;}"))
        comp = self._combo.completer()
        if comp is not None:
            comp.setFilterMode(Qt.MatchFlag.MatchContains)
            comp.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            comp.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        self._combo.currentIndexChanged.connect(self._on_combo_changed)
        top.addWidget(self._combo)
        self._reset_btn = self._small_btn("Reset hist", self._reset_current_hists)
        self._reset_btn.setEnabled(False)
        top.addWidget(self._reset_btn)

        root.addLayout(top)

        self._info = QLabel("")
        self._info.setFont(QFont("Monospace", 10))
        self._info.setStyleSheet(themed("color:#8b949e;"))
        root.addWidget(self._info)

        # -- tabbed central area: Waveform + Cluster -----------------------
        self._tabs = QTabWidget()
        self._tabs.addTab(self._build_waveform_tab(), "Waveform")
        self._tabs.addTab(self._build_cluster_tab(),  "Cluster")
        root.addWidget(self._tabs, stretch=1)

        # -- advanced dock (hidden by default) -----------------------------
        self._adv_dock = self._build_advanced_dock()

        self.setStatusBar(QStatusBar())
        self._clear_plots()

        # Keyboard: ← / → to navigate prev / next.
        QShortcut(QKeySequence(Qt.Key.Key_Left),  self, activated=self._on_prev)
        QShortcut(QKeySequence(Qt.Key.Key_Right), self, activated=self._on_next)

    # ---- Tabs -----------------------------------------------------------

    def _build_waveform_tab(self) -> QWidget:
        """Original waveform viewer body: geo + waveform on the left, four
        histograms stacked on the right."""
        tab = QWidget()
        lay = QVBoxLayout(tab)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(3)

        # Cut-Settings toolbar — mirrors the web monitor's <div class="tcut-bar">
        # with the "Cut Settings…" / apply / show controls.  apply toggles
        # WaveformFilter.enable; show toggles overlay drawing (no filter
        # change).  Both fire a re-run so changes are visible immediately.
        cut_bar = QHBoxLayout()
        cut_bar.setContentsMargins(2, 0, 2, 0)
        cut_bar.setSpacing(8)
        self._cut_settings_btn = QPushButton("Cut Settings…")
        self._cut_settings_btn.setToolTip("Edit waveform peak filter ranges")
        self._cut_settings_btn.setStyleSheet(themed(
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;"
            "padding:3px 10px;font:bold 9pt Monospace;}"
            "QPushButton:hover{background:#30363d;color:#e6edf3;}"))
        self._cut_settings_btn.clicked.connect(self._open_cut_settings_dialog)
        cut_bar.addWidget(self._cut_settings_btn)
        self._cut_apply_cb = QCheckBox("apply")
        self._cut_apply_cb.setChecked(self._filter.enable)
        self._cut_apply_cb.setToolTip(
            "Apply the peak filter to histograms / geo coloring / clustering. "
            "When off, every analyzed peak counts.")
        self._cut_apply_cb.toggled.connect(self._on_cut_apply_toggled)
        cut_bar.addWidget(self._cut_apply_cb)
        self._cut_show_cb = QCheckBox("show")
        self._cut_show_cb.setChecked(self._filter_show)
        self._cut_show_cb.setToolTip(
            "Show cut-range overlays on the waveform plot.  Independent of "
            "apply — overlays can be hidden while the filter is active.")
        self._cut_show_cb.toggled.connect(self._on_cut_show_toggled)
        cut_bar.addWidget(self._cut_show_cb)
        cut_bar.addStretch(1)
        lay.addLayout(cut_bar)

        split = QSplitter(Qt.Orientation.Horizontal)

        # Left: geo view (square, top) + waveform plot (bottom)
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        left_lay.setSpacing(4)
        self._geo = WaveformGeoView()
        if self._hycal_modules:
            self._geo.set_modules(self._hycal_modules)
        self._geo.moduleClicked.connect(self._on_geo_clicked)
        self._geo.setSizePolicy(QSizePolicy.Policy.Expanding,
                                QSizePolicy.Policy.Expanding)
        left_lay.addWidget(self._geo, stretch=3)
        self._wave = WaveformPlotWidget()
        left_lay.addWidget(self._wave, stretch=1)
        split.addWidget(left)

        # Right: four histograms stacked vertically (each wide, long in x).
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(0)
        self._h_height   = Hist1DWidget()
        self._h_integral = Hist1DWidget()
        self._h_position = Hist1DWidget()
        self._h_npeaks   = Hist1DWidget()
        for hist in (self._h_height, self._h_integral,
                     self._h_position, self._h_npeaks):
            right_lay.addWidget(hist, stretch=1)
        split.addWidget(right)

        # Even 50/50 split between the geo+waveform column and the hist stack.
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setSizes([750, 750])
        lay.addWidget(split, stretch=1)
        return tab

    def _build_cluster_tab(self) -> QWidget:
        """HyCal heatmap + cluster panel.  Populated each event by
        ``_display_clusters``."""
        tab = QWidget()
        root = QVBoxLayout(tab)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(2)

        # Calibration warning banner — visible when no gain file has been
        # loaded (cal_factor==0 → every energize() returns 0 → no clusters).
        self._calib_warn_lbl = QLabel(
            "⚠ No HyCal calibration loaded — cluster energies will be 0 "
            "and no clusters will form.  Use File → Load HyCal calibration…")
        self._calib_warn_lbl.setStyleSheet(
            f"background:{THEME.DANGER}; color:#ffffff; "
            f"padding:4px 8px; font-weight: bold;")
        self._calib_warn_lbl.setWordWrap(True)
        self._calib_warn_lbl.setVisible(not self._hycal_calib_loaded)
        root.addWidget(self._calib_warn_lbl)

        body = QWidget()
        lay = QHBoxLayout(body)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)

        self._cluster_map = HyCalClusterMap()
        if self._hycal_modules:
            self._cluster_map.set_modules(self._hycal_modules)
        self._cluster_map.setSizePolicy(QSizePolicy.Policy.Expanding,
                                        QSizePolicy.Policy.Expanding)
        lay.addWidget(self._cluster_map, stretch=3)

        self._cluster_panel = ClusterPanel()
        self._cluster_panel.clusterSelected.connect(
            self._cluster_map.set_selected_cluster)
        lay.addWidget(self._cluster_panel, stretch=1)

        root.addWidget(body, stretch=1)
        return tab

    # ---- Advanced tuning dock ------------------------------------------

    def _build_advanced_dock(self):
        """Right-side collapsible dock exposing WaveConfig + HyCalClusterConfig.

        Both configs live in prad2py bindings; changes trigger a re-run of
        the current event (``_rerun_current``) so the effect is immediate.
        """
        dock = QDockWidget("Advanced tuning", self)
        dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea
                             | Qt.DockWidgetArea.LeftDockWidgetArea)
        dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable
                         | QDockWidget.DockWidgetFeature.DockWidgetMovable
                         | QDockWidget.DockWidgetFeature.DockWidgetFloatable)

        root = QWidget()
        rlay = QVBoxLayout(root)
        rlay.setContentsMargins(6, 6, 6, 6)

        # ---- WaveConfig ------------------------------------------------
        wg = QGroupBox("Waveform analyser")
        wf = QFormLayout(wg)
        self._adv_wave_spins: Dict[str, QDoubleSpinBox] = {}
        wave_fields = [
            ("peak_nsigma",      0.0,    50.0, 0.5,  "peak detection threshold (× pedestal RMS)"),
            ("min_peak_height",  0.0,  4096.0, 1.0,  "absolute floor on detected peak height (ADC)"),
            ("min_peak_ratio",   0.0,     1.0, 0.01, "secondary/primary peak ratio"),
            ("int_tail_ratio",   0.0,     1.0, 0.01, "tail integration cut"),
            ("ped_flatness",     0.0,  1000.0, 0.5,  "pedestal RMS ceiling"),
            ("clk_mhz",          1.0,  1000.0, 1.0,  "FADC clock (MHz)"),
        ]
        for name, lo, hi, step, tip in wave_fields:
            sp = QDoubleSpinBox()
            sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(3)
            val = getattr(self._wcfg, name, 0.0)
            sp.setValue(float(val))
            sp.setToolTip(tip)
            sp.valueChanged.connect(self._on_advanced_changed)
            wf.addRow(name, sp)
            self._adv_wave_spins[name] = sp

        self._adv_wave_int_spins: Dict[str, QSpinBox] = {}
        for name, lo, hi, tip in [
                ("smooth_order", 1, 16,  "kernel order (1 = identity, N gives 2N-1 taps)"),
                ("ped_nsamples", 1, 64,  "samples to use for pedestal"),
                ("ped_max_iter", 1, 100, "pedestal iteration cap"),
                ("overflow",     0, 65535, "overflow cutoff (ADC)")]:
            sp = QSpinBox()
            sp.setRange(lo, hi); sp.setValue(int(getattr(self._wcfg, name, 0)))
            sp.setToolTip(tip)
            sp.valueChanged.connect(self._on_advanced_changed)
            wf.addRow(name, sp)
            self._adv_wave_int_spins[name] = sp

        # Note: peak filter (time / integral / height / quality_bits) lives
        # in the Cut Settings dialog above the waveform plot — same UX as
        # the web monitor.  Keeping the dock focused on analyser knobs
        # avoids two places editing the same numbers.

        rlay.addWidget(wg)

        # ---- HyCalClusterConfig ---------------------------------------
        cg = QGroupBox("HyCal clustering")
        cf = QFormLayout(cg)
        self._adv_cluster_spins: Dict[str, object] = {}
        if self._hccl is not None:
            ccfg = self._hccl.get_config()
            cluster_fields = [
                ("min_module_energy",  0.0, 1e4, 0.1, "single-module threshold (MeV)"),
                ("min_center_energy",  0.0, 1e4, 0.1, "seed threshold (MeV)"),
                ("min_cluster_energy", 0.0, 1e5, 0.1, "total cluster threshold (MeV)"),
                ("log_weight_thres",   0.0, 20.0, 0.1, "log-weight offset"),
                ("least_split",        0.0, 1.0, 0.01, "min fraction to keep a split hit"),
                ("seed_time_window",   -1.0, 200.0, 0.5,
                    "Multi-pulse seed-time gate (ns).  ≤0 disables timing "
                    "gating (legacy single-pulse-per-module mode).  >0 lets "
                    "AddHit() be called once per pulse; FormClusters then "
                    "groups neighbours within ±this window of the seed pulse.\n"
                    "Persistent default: 'seed_time_window' under the 'hycal' "
                    "block in database/reconstruction_config.json."),
            ]
            for name, lo, hi, step, tip in cluster_fields:
                sp = QDoubleSpinBox()
                sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(3)
                sp.setValue(float(getattr(ccfg, name, 0.0)))
                sp.setToolTip(tip)
                sp.valueChanged.connect(self._on_advanced_changed)
                cf.addRow(name, sp)
                self._adv_cluster_spins[name] = sp

            for name, lo, hi, tip in [
                    ("min_cluster_size", 1, 100, "min modules in cluster"),
                    ("split_iter",       0, 100, "island-split iteration cap")]:
                sp = QSpinBox()
                sp.setRange(lo, hi); sp.setValue(int(getattr(ccfg, name, 0)))
                sp.setToolTip(tip)
                sp.valueChanged.connect(self._on_advanced_changed)
                cf.addRow(name, sp)
                self._adv_cluster_spins[name] = sp

            cbx = QCheckBox("corner_conn (include diagonal neighbors)")
            cbx.setChecked(bool(getattr(ccfg, "corner_conn", False)))
            cbx.toggled.connect(self._on_advanced_changed)
            cf.addRow(cbx)
            self._adv_cluster_spins["corner_conn"] = cbx
        else:
            cf.addRow(QLabel("(HyCalSystem not initialized)"))
        rlay.addWidget(cg)

        rlay.addStretch(1)

        # Snapshot initial widget values (from monitor_config.json + WaveConfig
        # defaults + HyCalClusterConfig defaults) so "Reset to defaults"
        # can restore them after arbitrary tuning.
        self._adv_defaults: Dict[str, object] = {}
        for name, sp in self._adv_wave_spins.items():
            self._adv_defaults[f"wave_{name}"] = sp.value()
        for name, sp in self._adv_wave_int_spins.items():
            self._adv_defaults[f"wave_{name}"] = sp.value()
        for name, w in self._adv_cluster_spins.items():
            if isinstance(w, QCheckBox):
                self._adv_defaults[f"cluster_{name}"] = w.isChecked()
            else:
                self._adv_defaults[f"cluster_{name}"] = w.value()

        reset_btn = QPushButton("Reset to defaults")
        reset_btn.clicked.connect(self._reset_advanced_defaults)
        rlay.addWidget(reset_btn)

        dock.setWidget(root)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        dock.hide()
        return dock

    def _reset_advanced_defaults(self):
        """Restore every Advanced-dock widget to its initial (post-load)
        value and re-run the current event once.  The peak filter lives in
        the Cut Settings dialog and has its own Reset button."""
        if not hasattr(self, "_adv_defaults"):
            return
        for name, sp in self._adv_wave_spins.items():
            sp.blockSignals(True)
            sp.setValue(self._adv_defaults[f"wave_{name}"])
            sp.blockSignals(False)
        for name, sp in self._adv_wave_int_spins.items():
            sp.blockSignals(True)
            sp.setValue(self._adv_defaults[f"wave_{name}"])
            sp.blockSignals(False)
        for name, w in self._adv_cluster_spins.items():
            w.blockSignals(True)
            if isinstance(w, QCheckBox):
                w.setChecked(self._adv_defaults[f"cluster_{name}"])
            else:
                w.setValue(self._adv_defaults[f"cluster_{name}"])
            w.blockSignals(False)
        # Single re-run after all widgets are restored.
        self._on_advanced_changed()

    def _on_advanced_changed(self, *_):
        """Push every dock value back into WaveConfig + HyCalClusterConfig,
        then schedule a debounced re-run of the current event.  Setattrs
        are wrapped in try/except so a type mismatch on any single field
        (e.g. pybind's strict int↔float) doesn't take down the whole dock."""
        def _safe(obj, name, value):
            try:
                setattr(obj, name, value)
            except Exception as exc:     # noqa: BLE001
                print(f"[advanced] {type(obj).__name__}.{name} "
                      f"setattr failed: {exc}", file=sys.stderr)

        for name, sp in self._adv_wave_spins.items():
            _safe(self._wcfg, name, float(sp.value()))
        for name, sp in self._adv_wave_int_spins.items():
            _safe(self._wcfg, name, int(sp.value()))

        if self._hccl is not None and self._adv_cluster_spins:
            ccfg = self._hccl.get_config()
            for name, widget in self._adv_cluster_spins.items():
                if isinstance(widget, QCheckBox):
                    _safe(ccfg, name, bool(widget.isChecked()))
                elif isinstance(widget, QSpinBox):
                    _safe(ccfg, name, int(widget.value()))
                else:
                    _safe(ccfg, name, float(widget.value()))
            self._hccl.set_config(ccfg)

        # Debounce: coalesce slider drags into one re-run.
        self._adv_redraw_timer.start(self._adv_debounce_ms)

    def _rerun_current_event(self):
        """Re-read + re-analyse the current event with the latest config."""
        if self._current_idx >= 0:
            self._goto(self._current_idx)

    # ---- Cut Settings dialog --------------------------------------------

    def _open_cut_settings_dialog(self):
        """Modal editor for the active peak filter — replicates the web
        monitor's Cut Settings dialog (resources/cut_dialog.js).  On Save,
        the new filter values replace ``self._filter`` (preserving
        ``enable``, which stays under the toolbar's "apply" toggle) and
        the current event is re-analysed."""
        dlg = CutSettingsDialog(self, self._filter, self._filter_default)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        new_filter = dlg.result_filter()
        new_filter.enable = self._filter.enable
        self._filter = new_filter
        self._rerun_current_event()

    def _on_cut_apply_toggled(self, on: bool):
        """User flipped the toolbar's "apply" checkbox — toggles
        WaveformFilter.enable so every peak counts when off, every peak
        is filtered when on.  Re-runs the current event to refresh hists
        / geo / clusters."""
        self._filter.enable = bool(on)
        self._rerun_current_event()

    def _on_cut_show_toggled(self, on: bool):
        """User flipped the "show" checkbox — toggles overlay drawing on
        the waveform plot only (filter values + apply state unchanged).
        Cheap repaint, no full re-analysis."""
        self._filter_show = bool(on)
        if self._current_idx >= 0 and self._ch is not None:
            self._display_waveform(self._ch.fadc())
        else:
            self._wave.update()

    def _small_btn(self, text: str, slot, primary: bool = False) -> QPushButton:
        btn = QPushButton(text)
        bg = "#1f6feb" if primary else "#21262d"
        fg = "#ffffff" if primary else "#c9d1d9"
        btn.setStyleSheet(themed(
            f"QPushButton{{background:{bg};color:{fg};"
            f"border:1px solid #30363d;padding:5px 14px;"
            f"font:bold 10pt Monospace;border-radius:3px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{background:#161b22;color:#484f58;}}"))
        btn.clicked.connect(slot)
        return btn

    def _mk_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 10))
        lbl.setStyleSheet(themed("color:#c9d1d9;"))
        return lbl

    def _make_menu(self):
        mb = self.menuBar()
        mf = mb.addMenu("&File")

        a_open = QAction("Open &evio…", self)
        a_open.setShortcut("Ctrl+O")
        a_open.triggered.connect(self._open_evio_dialog)
        mf.addAction(a_open)

        a_calib = QAction("Load HyCal &calibration…", self)
        a_calib.triggered.connect(self._load_hycal_calib_dialog)
        mf.addAction(a_calib)

        self._a_save = QAction("&Save histograms as JSON…", self)
        self._a_save.setShortcut("Ctrl+S")
        self._a_save.triggered.connect(self._save_json_dialog)
        self._a_save.setEnabled(False)
        mf.addAction(self._a_save)

        mf.addSeparator()
        a_quit = QAction("&Quit", self)
        a_quit.setShortcut("Ctrl+Q")
        a_quit.triggered.connect(self.close)
        mf.addAction(a_quit)

        mv = mb.addMenu("&View")
        a_adv = self._adv_dock.toggleViewAction()
        a_adv.setText("&Advanced tuning")
        a_adv.setShortcut("Ctrl+T")
        mv.addAction(a_adv)

    # -- HyCal pipeline (PipelineBuilder) ------------------------------

    def _guess_database_dir(self) -> Optional[str]:
        """Pick a database dir for resolving recon-config-internal paths
        (e.g., 'runinfo': 'runinfo/general.json').  We prefer the
        parent of whichever full-path config the user actually pointed at
        — recon, daq, hycal-map, in that order — so the resolved paths
        stay self-consistent with the rest of their tree."""
        for p in (self._recon_cfg_path, self._daq_cfg_path,
                  self._hycal_map_path):
            if p:
                parent = Path(p).expanduser().resolve().parent
                if parent.is_dir():
                    return str(parent)
        return None

    def _build_hycal_pipeline(self, evio_path: Optional[str] = None) -> bool:
        """Wire up HyCal via prad2py.det.PipelineBuilder.

        Pulls daq / recon / runinfo / map / per-run calibration from the
        configured paths so cluster energies are populated without a
        manual ``Load HyCal calibration…`` step.  ``evio_path`` lets the
        builder pick the right run number from the file name (otherwise
        it uses the latest entry in runinfo).

        Returns True on success, False if the build raised — callers
        should then drop back to ``_hycal_init_fallback`` so the rest
        of the viewer still functions for waveform-only browsing.
        """
        if not _HAVE_PRAD2PY:
            return False
        try:
            b = prad2py.det.PipelineBuilder()
            # Anchor relative paths inside recon_config (e.g. "runinfo": …)
            # to whichever directory the configs we were given live in,
            # otherwise the builder resolves them against CWD and silently
            # drops calibration when run from outside the repo.
            db_dir = self._guess_database_dir()
            if db_dir:
                b.set_database_dir(db_dir)
            if self._daq_cfg_path:
                b.set_daq_config(self._daq_cfg_path)
            if self._recon_cfg_path:
                b.set_recon_config(self._recon_cfg_path)
            if self._hycal_map_path:
                b.set_hycal_map(self._hycal_map_path)
            if evio_path:
                b.set_run_number_from_evio(str(evio_path))
            pipeline = b.build()
        except Exception as exc:                      # noqa: BLE001
            print(f"[hycal] PipelineBuilder.build failed: "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
            return False

        self._pipeline = pipeline
        self._hcsys = pipeline.hycal
        self._hccl = prad2py.det.HyCalCluster(pipeline.hycal)
        try:
            self._hccl.set_config(pipeline.hycal_cluster_cfg)
        except Exception as exc:                      # noqa: BLE001
            print(f"[hycal] HyCalCluster.set_config failed: "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
        self._hc_cache.clear()
        # PipelineBuilder leaves hycal_calib_path empty when no calibration
        # file resolved (no runinfo / no energy_calib_file).  Truthy =
        # calibration actually loaded → suppress the cluster-tab warning.
        self._hycal_calib_loaded = bool(getattr(pipeline, "hycal_calib_path", "")
                                        or "")
        # Sync the Advanced dock's cluster spinboxes to the new config so
        # the user sees the values that actually drive reconstruction.
        # `hasattr` guard: this also runs from __init__ before the dock exists.
        if hasattr(self, "_adv_cluster_spins"):
            self._sync_advanced_dock_cluster_cfg()
        if hasattr(self, "_calib_warn_lbl"):
            self._calib_warn_lbl.setVisible(not self._hycal_calib_loaded)
        return True

    def _hycal_init_fallback(self) -> bool:
        """Fallback to a plain HyCalSystem.init when PipelineBuilder fails
        (e.g., daq_config missing).  Provides geometry-only operation: no
        calibration, so cluster energies stay at 0 until the user picks a
        file via File → Load HyCal calibration…"""
        self._pipeline = None
        if not (_HAVE_PRAD2PY and self._hycal_map_path):
            self._hcsys = None
            self._hccl = None
            return False
        try:
            sys_obj = prad2py.det.HyCalSystem()
            if not sys_obj.init(str(self._hycal_map_path)):
                print("[hycal] HyCalSystem.init returned False — "
                      "cluster tab will be empty", file=sys.stderr)
                self._hcsys = None
                self._hccl = None
                return False
            self._hcsys = sys_obj
            self._hccl = prad2py.det.HyCalCluster(sys_obj)
        except Exception as exc:                      # noqa: BLE001
            print(f"[hycal] init failed: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            self._hcsys = None
            self._hccl = None
            return False
        self._hc_cache.clear()
        self._hycal_calib_loaded = False
        return True

    def _sync_advanced_dock_cluster_cfg(self):
        """Refresh the Advanced-dock cluster spinboxes from _hccl's config.
        Called after a pipeline rebuild so the dock matches the new run's
        recon-config defaults instead of stale values."""
        if self._hccl is None or not getattr(self, "_adv_cluster_spins", None):
            return
        try:
            ccfg = self._hccl.get_config()
        except Exception:
            return
        for name, widget in self._adv_cluster_spins.items():
            try:
                v = getattr(ccfg, name)
            except AttributeError:
                continue
            widget.blockSignals(True)
            try:
                if isinstance(widget, QCheckBox):
                    widget.setChecked(bool(v))
                elif isinstance(widget, QSpinBox):
                    widget.setValue(int(v))
                else:
                    widget.setValue(float(v))
            finally:
                widget.blockSignals(False)

    # -- file open --

    def _open_evio_dialog(self):
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open evio file", str(Path.cwd()),
            "evio files (*.evio *.evio.*);;All files (*)")
        if path_str:
            self.open_path(Path(path_str))

    def _load_hycal_calib_dialog(self):
        """Load HyCal per-module calibration constants (cal_factor,
        cal_base_energy, cal_non_linear) from a JSON file via
        HyCalSystem.LoadCalibration.  Success hides the "no calibration"
        warning banner on the cluster tab and re-runs the current event
        so cluster energies reflect the new calibration."""
        if self._hcsys is None:
            QMessageBox.warning(self, "HyCal system not ready",
                "Cannot load calibration — HyCalSystem.init() did not "
                "complete.  Check that hycal_map.json is present.")
            return
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Load HyCal calibration", str(Path.cwd()),
            "JSON files (*.json);;All files (*)")
        if not path_str:
            return
        try:
            nmatched = self._hcsys.load_calibration(path_str)
        except Exception as exc:             # noqa: BLE001
            QMessageBox.critical(self, "Calibration load failed",
                                 f"{type(exc).__name__}: {exc}")
            return
        if nmatched is None or int(nmatched) <= 0:
            QMessageBox.warning(self, "Calibration load failed",
                f"load_calibration returned {nmatched!r} — 0 modules "
                f"matched.  Check that the file format matches "
                f"HyCalSystem::LoadCalibration's expected schema.")
            return
        self._hycal_calib_loaded = True
        if hasattr(self, "_calib_warn_lbl"):
            self._calib_warn_lbl.setVisible(False)
        self.statusBar().showMessage(
            f"Calibration loaded: {Path(path_str).name} "
            f"({int(nmatched)} modules matched)", 5000)
        # Re-run the current event so cluster energies reflect the new
        # cal_factors.  No-op if no event has been loaded yet.
        if self._current_idx >= 0:
            self._goto(self._current_idx)

    def open_path(self, path: Path):
        err = _check_evchannel_support()
        if err:
            QMessageBox.critical(self, "prad2py issue", err)
            return
        if self._idx_thread is not None:
            QMessageBox.information(self, "Busy", "Already indexing.")
            return

        # Tear down any previous reader / hists
        self._close_reader()
        self._channels.clear()
        self._selected_key = None
        self._combo.blockSignals(True); self._combo.clear(); self._combo.blockSignals(False)
        self._index = []
        self._current_idx = -1
        self._accumulated = None
        self._clear_plots()

        # Rebuild the HyCal pipeline against this evio's run number so the
        # right per-run calibration is loaded.  Failures fall back to
        # whatever pipeline / map-only state was already in place — the
        # waveform side keeps working regardless.
        self._build_hycal_pipeline(evio_path=str(path))

        # Start indexer
        dlg = QProgressDialog(f"Indexing {path.name} …", "Cancel", 0, 100, self)
        dlg.setWindowTitle("Indexing")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()

        worker = IndexerWorker(str(path), self._daq_cfg_path)
        thread = QThread(self)
        worker.moveToThread(thread)

        def _on_progress(done: int, total: int):
            if total > 0:
                dlg.setMaximum(total)
                dlg.setValue(done)
            dlg.setLabelText(f"Indexing {path.name}\n"
                             f"evio events: {done:,} / {total:,}")

        thread.started.connect(worker.run)
        worker.progressed.connect(_on_progress)
        worker.finished.connect(lambda res: self._on_index_done(path, res))
        worker.failed.connect(lambda msg: self._on_index_failed(path, msg))
        dlg.canceled.connect(worker.request_cancel)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(dlg.close)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(lambda: setattr(self, "_idx_thread", None))

        self._idx_worker = worker
        self._idx_thread = thread
        thread.start()

    def _on_index_done(self, path: Path, res: Dict):
        self._evio_path = path
        self._index = res["index"]
        n_phys = len(self._index)
        cancelled = bool(res.get("cancelled"))
        self._accumulated = np.zeros(n_phys, dtype=bool)

        # Open reader handle for browse use
        ok, err = self._open_reader(str(path))
        if not ok:
            QMessageBox.critical(self, "Open failed",
                                 f"Indexing finished but reader open failed:\n{err}")
            return

        mode_note = ("" if self._reader_is_ra
                     else "   [sequential mode — Prev is slow]")
        self._file_lbl.setText(
            f"{path.name}   physics events: {n_phys:,}   "
            f"(evio blocks: {res['total_evio']:,})"
            + mode_note
            + ("   [indexing cancelled]" if cancelled else ""))
        self._info.setText("Select a module, then click Next to start browsing.")

        self._event_spin.setMaximum(max(0, n_phys - 1))
        self._event_spin.setValue(0)
        if n_phys > 0:
            self._prev_btn.setEnabled(True)
            self._next_btn.setEnabled(True)
            self._event_spin.setEnabled(True)
            self._batch_btn.setEnabled(True)
            self._a_save.setEnabled(True)
            self._reset_btn.setEnabled(True)
        self._total_lbl.setText(f" / {max(0, n_phys - 1):,}")

        self.statusBar().showMessage(
            f"Indexed {n_phys:,} physics events from {path.name}")

    def _on_index_failed(self, path: Path, msg: str):
        QMessageBox.critical(self, "Indexing failed", f"{path}\n\n{msg}")
        self.statusBar().showMessage(f"Failed to index {path.name}")

    # -- reader (browse handle) --

    def _open_reader(self, path: str) -> Tuple[bool, str]:
        try:
            dec = prad2py.dec
            cfg = (dec.load_daq_config(self._daq_cfg_path) if self._daq_cfg_path
                   else dec.load_daq_config())
            ch  = dec.EvChannel()
            ch.set_config(cfg)
            st = ch.open_auto(path)
            if st != dec.Status.success:
                return False, f"status = {st}"
            self._ch = ch
            self._reader_path = path
            self._reader_is_ra = bool(ch.is_random_access())
            self._reader_pos = -1     # sequential-mode cursor
            return True, ""
        except Exception as e:
            return False, f"{type(e).__name__}: {e}"

    def _close_reader(self):
        if self._ch is not None:
            try:
                self._ch.close()
            except Exception:
                pass
        self._ch = None
        self._reader_path = None
        self._reader_is_ra = False
        self._reader_pos = -1

    # -- navigation --

    def _on_prev(self):
        if self._current_idx > 0:
            self._goto(self._current_idx - 1)
        elif self._current_idx == -1 and self._index:
            self._goto(0)

    def _on_next(self):
        if self._current_idx < len(self._index) - 1:
            self._goto(self._current_idx + 1)

    def _on_spin_jump(self):
        v = self._event_spin.value()
        if 0 <= v < len(self._index) and v != self._current_idx:
            self._goto(v)

    def _goto(self, phys_idx: int):
        if self._ch is None or not (0 <= phys_idx < len(self._index)):
            return
        ev_idx, sub_idx = self._index[phys_idx]
        dec = prad2py.dec
        # RA: jump in O(1).  Sequential: close/reopen on backward jumps,
        # walk forward to the target.
        if self._reader_is_ra:
            st = self._ch.read_event_by_index(ev_idx)
            if st != dec.Status.success:
                self.statusBar().showMessage(
                    f"read_event_by_index({ev_idx}) → {st}")
                return
        else:
            if self._reader_pos > ev_idx:
                # Backward seek — reopen and walk forward from start.
                self._close_reader()
                ok, err = self._open_reader(str(self._evio_path))
                if not ok:
                    self.statusBar().showMessage(
                        f"reopen for backward seek failed: {err}")
                    return
            while self._reader_pos < ev_idx:
                if self._ch.read() != dec.Status.success:
                    self.statusBar().showMessage(
                        f"EOF before evio event {ev_idx}")
                    return
                self._reader_pos += 1
        if not self._ch.scan():
            self.statusBar().showMessage(f"scan() failed at physics #{phys_idx}")
            return
        self._ch.select_event(sub_idx)
        info = self._ch.info()
        tb = int(info.trigger_bits)

        self._current_idx = phys_idx
        self._event_spin.blockSignals(True)
        self._event_spin.setValue(phys_idx)
        self._event_spin.blockSignals(False)

        # Apply trigger filter: skip updating hists but still show waveform
        trig_ok = True
        if self._accept_mask and (tb & self._accept_mask) == 0: trig_ok = False
        if self._reject_mask and (tb & self._reject_mask):       trig_ok = False

        fadc_evt = self._ch.fadc()
        self._update_channel_list_from_event(fadc_evt)
        self._accumulate_and_display(fadc_evt, info, trig_ok)

    def _update_channel_list_from_event(self, fadc_evt):
        """Add any new (roc, slot, ch) seen in this event to the combo."""
        added = False
        for r in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            if roc_tag not in self._roc_to_crate:
                continue
            crate = self._roc_to_crate[roc_tag]
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = (roc_tag, s, c)
                    if key not in self._channels:
                        module = self._daq_map.get((crate, s, c))
                        self._channels[key] = _make_hists(
                            self._h_cfg, self._i_cfg, self._p_cfg, self._n_cfg,
                            roc_tag, s, c, module)
                        added = True
        if added:
            self._refresh_combo()
            self._geo.set_available({c.module for c in self._channels.values()
                                     if c.module})

    def _refresh_combo(self):
        """Re-populate combo from discovered channels, preserving selection."""
        prev_key = self._selected_key
        items: List[Tuple[str, Tuple[int, int, int], str]] = []  # (sort_key, key, label)
        for key, ch in self._channels.items():
            mod = ch.module or "(unmapped)"
            label = (f"{mod:<8}  roc=0x{ch.roc:02X}  s={ch.slot:>2}  "
                     f"ch={ch.channel:>2}")
            sort_key = (0 if ch.module else 1, _natural_sort_key(mod))
            items.append((sort_key, key, label))
        items.sort(key=lambda x: x[0])

        self._combo.blockSignals(True)
        self._combo.clear()
        self._combo_keys: List[Tuple[int, int, int]] = []
        sel_idx = 0
        for i, (_, key, label) in enumerate(items):
            self._combo.addItem(label)
            self._combo_keys.append(key)
            if key == prev_key:
                sel_idx = i
        if self._combo_keys:
            self._combo.setCurrentIndex(sel_idx)
            self._selected_key = self._combo_keys[sel_idx]
        self._combo.blockSignals(False)

    def _on_combo_changed(self, idx: int):
        if not hasattr(self, "_combo_keys"):
            return
        if 0 <= idx < len(self._combo_keys):
            self._selected_key = self._combo_keys[idx]
            hits = self._channels.get(self._selected_key)
            self._geo.set_selected_module(hits.module if hits else None)
            # Re-render: hists from cache, waveform from current event
            self._display_hists_for_selected()
            if self._current_idx >= 0 and self._ch is not None:
                # Re-pull the waveform for the new module from current event
                fadc_evt = self._ch.fadc()
                self._display_waveform(fadc_evt)

    def _on_geo_clicked(self, name: str):
        """Geo-view click: switch combo to the (first) channel for this module."""
        if not name:
            return
        for i, key in enumerate(getattr(self, "_combo_keys", [])):
            hits = self._channels.get(key)
            if hits and hits.module == name:
                self._combo.setCurrentIndex(i)   # triggers _on_combo_changed
                return
        # Module exists in geometry but hasn't been seen yet — ignore the click
        self.statusBar().showMessage(
            f"{name}: no events seen yet for this module", 2000)

    # -- accumulate + display --

    def _accumulate_and_display(self, fadc_evt, info, trig_ok: bool):
        # Analyse every channel present in the event.  The C++ WaveAnalyzer
        # makes this cheap (~50× faster than the old Python port), so we
        # no longer need a "selected-channel-only" fast path.  Histogram
        # fills are dedup'd via self._accumulated: an event already folded
        # in still gets re-analysed for display, but isn't counted again.
        sel_peaks: List[Peak] = []
        wcfg = self._wcfg
        flt = self._filter
        sel_key = self._selected_key
        idx = self._current_idx
        already = (self._accumulated is not None and 0 <= idx < self._accumulated.size
                   and bool(self._accumulated[idx]))
        do_fill = trig_ok and not already

        current_vals: Dict[str, float] = {}   # module_name -> max peak integral
        module_energies: Dict[str, float] = {}  # name -> MeV (cluster tab)
        # Modules where the analyser found at least one peak rejected by
        # the active filter — drives the geo tooltip so users can tell
        # "no signal" apart from "rejected by filter".
        rejected_current: set = set()

        # Reset clustering state for this event.
        if self._hccl is not None:
            self._hccl.clear()

        for r in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            crate = self._roc_to_crate.get(roc_tag)
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = (roc_tag, s, c)
                    hits = self._channels.get(key)
                    if hits is None:
                        continue
                    samples = slot.channel(c).samples
                    if samples.size < 10:
                        continue
                    _, _, peaks = analyze(samples, wcfg)
                    if key == sel_key:
                        sel_peaks = peaks
                    # Pick the best (highest-integral) peak that passes the
                    # active filter; matches viewer_utils.h::bestPeakInWindow
                    # plus the integral/height/quality cuts in PeakFilter.
                    # ``best_time`` is needed for HyCalCluster.add_hit's
                    # multi-pulse seed-time gating (commit 0dafbae).
                    best_int = 0.0
                    best_time = 0.0
                    any_peak = len(peaks) > 0
                    for p in peaks:
                        if not flt.passes(p):
                            continue
                        if p.integral > best_int:
                            best_int = p.integral
                            best_time = float(p.time)
                    max_int = best_int
                    if max_int > 0 and hits.module:
                        current_vals[hits.module] = max_int
                    elif any_peak and hits.module:
                        rejected_current.add(hits.module)

                    # Feed the cluster tab: resolve channel → HyCal module
                    # and push (module_idx, energy_MeV, time_ns).  Without
                    # calibration loaded, cal_factor==0 → energize returns
                    # 0 → no clusters form; the cluster tab surfaces a
                    # warning banner (File → Load HyCal calibration…).
                    if self._hccl is not None and crate is not None and max_int > 0:
                        mod = self._resolve_hycal_module(crate, s, c, key)
                        if mod is not None:
                            energy = mod.energize(max_int)
                            if energy > 0:
                                self._hccl.add_hit(mod.index, energy, best_time)
                                module_energies[mod.name] = energy

                    if not do_fill:
                        continue
                    kept = 0
                    for p in peaks:
                        if not flt.passes(p):
                            continue
                        hits.height.fill(p.height)
                        hits.integral.fill(p.integral)
                        hits.position.fill(p.time)
                        kept += 1
                    hits.npeaks.fill(kept)
                    hits.events += 1
                    if kept > 0:
                        hits.peak_events += 1

        if do_fill and self._accumulated is not None and 0 <= idx < self._accumulated.size:
            self._accumulated[idx] = True

        self._geo.set_rejected_current(rejected_current)
        self._geo.set_current_values(current_vals)
        # Overall occupancy only needs refreshing when hists actually changed.
        if do_fill:
            self._geo.set_overall_values(self._compute_overall_occupancy())

        if sel_key is None:
            self._set_info_line(info, peaks=None)
            self._wave.clear("(select a module to view its waveform)")
        else:
            self._set_info_line(info, peaks=sel_peaks)
        self._display_hists_for_selected()
        self._display_waveform(fadc_evt)
        self._display_clusters(module_energies)

    def _resolve_hycal_module(self, crate: int, slot: int, ch: int,
                              cache_key: Tuple[int, int, int]):
        """Look up a HyCal Module by DAQ address, caching on the ROC-tag key.
        Returns the Module or None if this channel isn't a HyCal module
        (LMS / Veto / scaler channels return None and are silently ignored)."""
        if self._hcsys is None:
            return None
        if cache_key in self._hc_cache:
            return self._hc_cache[cache_key]
        m = self._hcsys.module_by_daq(int(crate), int(slot), int(ch))
        # Filter: only HyCal PWO4/PbGlass modules contribute to clustering.
        if m is not None and not m.is_hycal():
            m = None
        self._hc_cache[cache_key] = m
        return m

    def _display_clusters(self, module_energies: Dict[str, float]):
        """Push per-event energies + reconstructed clusters to the Cluster tab."""
        if self._hccl is None:
            return
        # Push energies to the heatmap
        if module_energies:
            vmax = max(module_energies.values())
        else:
            vmax = 1.0
        self._cluster_map.set_values(module_energies)
        self._cluster_map.set_range(0.0, vmax if vmax > 0 else 1.0)

        # Reconstruct + wrap for display
        self._hccl.form_clusters()
        matched = self._hccl.reconstruct_matched()
        display: List["_DisplayCluster"] = []
        for rr in matched:
            mc = rr.cluster
            hit = rr.hit
            centre = mc.center
            cname = ""
            members = []
            try:
                cmod = self._hcsys.module(centre.index)
                cname = cmod.name
            except Exception:
                pass
            for h in mc.hits:
                try:
                    members.append(self._hcsys.module(h.index).name)
                except Exception:
                    pass
            display.append(_DisplayCluster(hit, cname, members))

        self._cluster_map.set_clusters(display)
        self._cluster_panel.set_clusters(display)

    def _compute_overall_occupancy(self) -> Dict[str, float]:
        """module_name -> events_with_peak / events_accumulated (skip empty)."""
        out: Dict[str, float] = {}
        for hits in self._channels.values():
            if hits.module and hits.events > 0:
                out[hits.module] = hits.peak_events / hits.events
        return out

    def _set_info_line(self, info, peaks: Optional[List[Peak]]):
        pieces = [
            f"event #{self._current_idx:,}",
            f"tb=0x{int(info.trigger_bits):X}",
            f"evnum={int(info.event_number)}",
        ]
        if peaks is not None:
            pieces.append(f"peaks(this view)={len(peaks)}")
        key = self._selected_key
        if key:
            hits = self._channels.get(key)
            if hits:
                pieces.append(f"accum={hits.events:,}")
                pieces.append(f"w/peak={hits.peak_events:,}")
        self._info.setText("   ".join(pieces))

    def _display_hists_for_selected(self):
        key = self._selected_key
        hits = self._channels.get(key) if key else None
        if hits is None:
            self._clear_hists()
            return
        mod = hits.module or "(unmapped)"
        self._h_height.set_data(hits.height.bins, hits.height.bmin,
                                hits.height.bstep,
                                under=hits.height.under, over=hits.height.over,
                                title=f"{mod}  —  Peak Height [ADC]",
                                color="#e599f7")
        self._h_integral.set_data(hits.integral.bins, hits.integral.bmin,
                                  hits.integral.bstep,
                                  under=hits.integral.under, over=hits.integral.over,
                                  title=f"{mod}  —  Peak Integral [ADC·sample]",
                                  color="#00b4d8")
        self._h_position.set_data(hits.position.bins, hits.position.bmin,
                                  hits.position.bstep,
                                  under=hits.position.under, over=hits.position.over,
                                  title=f"{mod}  —  Peak Time [ns]",
                                  color="#51cf66")
        self._h_npeaks.set_data(hits.npeaks.bins, hits.npeaks.bmin,
                                hits.npeaks.bstep,
                                under=hits.npeaks.under, over=hits.npeaks.over,
                                title=f"{mod}  —  Peaks / Event",
                                color="#ffa657")

    def _display_waveform(self, fadc_evt):
        key = self._selected_key
        if not key:
            self._wave.clear("(select a module)")
            return
        roc_tag, slot, ch = key
        stack_key = f"{roc_tag:02X}_{slot}_{ch}"
        samples = _find_channel_samples(fadc_evt, roc_tag, slot, ch)
        if samples is None or samples.size == 0:
            # Keep the stacker intact on empty events — matches
            # resources/waveform.js:105 — but still clear it when the
            # user has switched to a different module.
            if self._wave.is_stacking():
                self._wave.reset_stack_if_new_key(stack_key)
                return
            hits = self._channels.get(key)
            mod = hits.module if hits and hits.module else "(unmapped)"
            self._wave.clear(f"{mod} not present in event #{self._current_idx}")
            return
        ped_mean, ped_rms, peaks = analyze(samples, self._wcfg)
        hits = self._channels.get(key)
        mod = hits.module if hits and hits.module else "(unmapped)"
        self._wave.set_data(samples, peaks, ped_mean, ped_rms,
                            title=(f"{mod}   roc=0x{roc_tag:02X}  "
                                   f"slot={slot}  ch={ch}"),
                            clk_mhz=self._wcfg.clk_mhz,
                            stack_key=stack_key,
                            filter=self._filter,
                            filter_show=self._filter_show)

    def _clear_hists(self):
        self._h_height.clear("Peak Height")
        self._h_integral.clear("Peak Integral")
        self._h_position.clear("Peak Time")
        self._h_npeaks.clear("Peaks / Event")

    def _clear_plots(self):
        self._clear_hists()
        self._wave.clear("(open an evio file and click Next)")

    def _reset_current_hists(self):
        key = self._selected_key
        if not key: return
        hits = self._channels.get(key)
        if not hits: return
        hits.height.reset()
        hits.integral.reset()
        hits.position.reset()
        hits.npeaks.reset()
        hits.events = 0
        hits.peak_events = 0
        self._display_hists_for_selected()
        self._geo.set_overall_values(self._compute_overall_occupancy())
        self.statusBar().showMessage(
            f"Reset histograms for {hits.module or '(unmapped)'}")

    # -- batch 10k --

    def _on_batch_10k(self):
        if self._batch_thread is not None:
            QMessageBox.information(self, "Busy", "Batch already running.")
            return
        start_idx = max(0, self._current_idx + 1)
        if start_idx >= len(self._index):
            QMessageBox.information(self, "End of file",
                                    "Already at the last physics event.")
            return
        remaining = len(self._index) - start_idx
        count = min(10_000, remaining)

        worker = BatchWorker(
            evio_path=str(self._evio_path),
            daq_config_path=self._daq_cfg_path,
            index=self._index, start_idx=start_idx, count=count,
            channels=self._channels, wcfg=self._wcfg,
            peak_filter=self._filter,
            accept_mask=self._accept_mask, reject_mask=self._reject_mask,
            accumulated=self._accumulated,
        )
        thread = QThread(self)
        worker.moveToThread(thread)

        # Modal progress dialog — blocks input to the main window until the
        # batch finishes (or the user cancels), so they can't switch modules
        # / reload / Prev / Next while hists are being filled underneath.
        dlg = QProgressDialog(
            f"Processing {count:,} events…", "Cancel", 0, count, self)
        dlg.setWindowTitle("Accumulating")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setAutoReset(False)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()

        def _on_progress(done: int, target: int, peaks: int):
            dlg.setValue(done)
            dlg.setLabelText(
                f"Processing {target:,} events\n"
                f"done: {done:,} / {target:,}   peaks found: {peaks:,}")
            self._batch_status.setText(
                f"batch: {done:,}/{target:,}  peaks={peaks:,}")
            # refresh hists + geo overall map incrementally
            self._display_hists_for_selected()
            self._geo.set_overall_values(self._compute_overall_occupancy())

        def _on_finished(n: int):
            self._current_idx = start_idx + n - 1
            self._event_spin.blockSignals(True)
            self._event_spin.setValue(max(0, self._current_idx))
            self._event_spin.blockSignals(False)
            self._batch_status.setText(
                f"batch done: {n:,} events processed")
            self._display_hists_for_selected()
            self._geo.set_overall_values(self._compute_overall_occupancy())
            # advance one more to show the next waveform
            if self._current_idx + 1 < len(self._index):
                self._goto(self._current_idx + 1)

        def _on_failed(msg: str):
            QMessageBox.critical(self, "Batch failed", msg)
            self._batch_status.setText("batch failed")

        def _cleanup():
            dlg.close()
            self._batch_btn.setEnabled(True)
            self._batch_worker = None
            self._batch_thread = None

        thread.started.connect(worker.run)
        worker.progressed.connect(_on_progress)
        worker.finished.connect(_on_finished)
        worker.failed.connect(_on_failed)
        dlg.canceled.connect(worker.request_cancel)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(_cleanup)

        self._batch_worker = worker
        self._batch_thread = thread
        self._batch_btn.setEnabled(False)
        self._batch_status.setText(f"batch: 0/{count:,}")
        thread.start()

    # -- JSON save --

    def _save_json_dialog(self):
        if not self._channels:
            return
        default = (self._evio_path.name + ".waveform.json"
                   if self._evio_path else "waveform_hist.json")
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save histograms as JSON", str(Path.cwd() / default),
            "JSON files (*.json)")
        if not path_str:
            return
        try:
            out = {
                "source_file": str(self._evio_path) if self._evio_path else "",
                "height_hist":   {"min": self._h_cfg["min"], "max": self._h_cfg["max"],
                                  "step": self._h_cfg["step"]},
                "integral_hist": {"min": self._i_cfg["min"], "max": self._i_cfg["max"],
                                  "step": self._i_cfg["step"]},
                "position_hist": {"min": self._p_cfg["min"], "max": self._p_cfg["max"],
                                  "step": self._p_cfg["step"]},
                "npeaks_hist":   {"min": self._n_cfg["min"], "max": self._n_cfg["max"],
                                  "step": self._n_cfg["step"]},
                "filter":        self._filter.to_json(),
                "filter_active": self._filter.enable,
                "wave_config":   self._wcfg.__dict__.copy(),
                "channels":      {},
            }
            for (roc, slot, ch), hits in sorted(self._channels.items()):
                out["channels"][f"{roc}_{slot}_{ch}"] = {
                    "module":      hits.module,
                    "roc":         hits.roc,
                    "slot":        hits.slot,
                    "channel":     hits.channel,
                    "events":      hits.events,
                    "peak_events": hits.peak_events,
                    "height_hist":   hits.height.to_json(),
                    "integral_hist": hits.integral.to_json(),
                    "position_hist": hits.position.to_json(),
                    "npeaks_hist":   hits.npeaks.to_json(),
                }
            Path(path_str).parent.mkdir(parents=True, exist_ok=True)
            with open(path_str, "w", encoding="utf-8") as f:
                json.dump(out, f)
            self.statusBar().showMessage(f"Saved {path_str}")
        except Exception as ex:
            QMessageBox.warning(self, "Save failed", f"{path_str}\n\n{ex}")

    # -- close --

    def closeEvent(self, ev):
        if self._idx_worker is not None:
            self._idx_worker.request_cancel()
        if self._batch_worker is not None:
            self._batch_worker.request_cancel()
        for thr in (self._idx_thread, self._batch_thread):
            if thr is not None and thr.isRunning():
                thr.quit()
                thr.wait(3000)
        self._close_reader()
        super().closeEvent(ev)


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="HyCal Event Viewer — browse an evio file event-by-event "
                    "with Waveform and Cluster tabs.")
    ap.add_argument("path", nargs="?", type=Path,
                    help="evio file to open (otherwise use File → Open…).")
    ap.add_argument("--config", type=Path,
                    default=_REPO_DIR / "database" / "monitor_config.json",
                    help="monitor_config.json (waveform binning).")
    ap.add_argument("--daq-config", type=Path,
                    default=_REPO_DIR / "database" / "daq_config.json",
                    help="daq_config.json (ROC-tag → crate mapping).")
    ap.add_argument("--hycal-map", type=Path,
                    default=_REPO_DIR / "database" / "hycal_map.json",
                    help="hycal_map.json (module geometry + DAQ map).")
    ap.add_argument("--recon-config", type=Path,
                    default=_REPO_DIR / "database" / "reconstruction_config.json",
                    help="reconstruction_config.json (runinfo pointer + "
                         "HyCal cluster config).  Resolved by PipelineBuilder "
                         "to load the per-run HyCal calibration automatically.")
    ap.add_argument("--trigger-bits", type=Path,
                    default=_REPO_DIR / "database" / "trigger_bits.json",
                    help="trigger_bits.json (for --accept/--reject-trigger).")
    ap.add_argument("--accept-trigger", action="append", default=[],
                    metavar="NAME",
                    help="Require at least one of these trigger bits (repeatable).")
    ap.add_argument("--reject-trigger", action="append", default=None,
                    metavar="NAME",
                    help="Drop events with any of these trigger bits (repeatable). "
                         "Default: uses monitor_config.json setting.")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

    hist_cfg      = load_hist_config(args.config)      if args.config.is_file()      else {}
    roc_to_crate  = load_roc_tag_map(args.daq_config)  if args.daq_config.is_file()  else {}
    daq_map       = load_daq_map(args.hycal_map)       if args.hycal_map.is_file()   else {}
    bit_map       = load_trigger_bit_map(args.trigger_bits)
    hycal_modules = (load_geo_modules(args.hycal_map)
                     if args.hycal_map.is_file() else [])

    accept_names = args.accept_trigger or hist_cfg.get("accept_trigger_bits", []) or []
    reject_names = (args.reject_trigger if args.reject_trigger is not None
                    else hist_cfg.get("reject_trigger_bits", []) or [])
    accept_mask = _mask_from_names(accept_names, bit_map) if accept_names else 0
    reject_mask = _mask_from_names(reject_names, bit_map) if reject_names else 0

    app = QApplication(sys.argv)
    win = HyCalEventViewer(
        hist_config       = hist_cfg,
        daq_map           = daq_map,
        roc_to_crate      = roc_to_crate,
        accept_mask       = accept_mask,
        reject_mask       = reject_mask,
        daq_config_path   = str(args.daq_config) if args.daq_config.is_file() else "",
        hycal_modules     = hycal_modules,
        hycal_map_path    = (str(args.hycal_map)
                             if args.hycal_map.is_file() else None),
        recon_config_path = (str(args.recon_config)
                             if args.recon_config.is_file() else None),
    )
    win.show()
    if args.path is not None:
        QTimer.singleShot(0, lambda: win.open_path(args.path))
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
