#!/usr/bin/env python3
"""
Online LMS Gain Monitor (PyQt6)
===============================

Polls clondaq2 for new EVIO files for a user-selected run, downloads missing
files, updates LMS gain-correction ROOT output through
``prad2ana_online_gain_monitor_base``, and displays per-batch W-module gain
time series.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
import re
import shlex
import shutil
import sys
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from zoneinfo import ZoneInfo

try:
    import numpy as np
    import uproot
    from uproot.source.file import MemmapSource
    from uproot.source.futures import TrivialExecutor
    HAS_UPROOT = True
except ImportError:
    np = None
    uproot = None
    MemmapSource = None
    TrivialExecutor = None
    HAS_UPROOT = False

from PyQt6.QtCore import QObject, QProcess, QThread, QTimer, Qt, pyqtSignal, pyqtSlot
from PyQt6.QtGui import QAction, QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import (
    QApplication, QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QHBoxLayout,
    QLabel, QLineEdit, QMessageBox, QMainWindow, QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter,
    QTextEdit, QVBoxLayout, QWidget,
)

from hycal_geoview import (
    HyCalMapWidget,
    THEME,
    apply_theme_palette,
    available_themes,
    load_modules,
    set_theme,
    themed,
)


SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR", SCRIPT_DIR / ".." / "database")).resolve()
MODULES_JSON = DB_DIR / "hycal_map.json"

REMOTE_HOST = "clondaq2"
REMOTE_BASE = "/data/stage2"
STORAGE_BASE = "/data/gain_monitor"
DEFAULT_INTERVAL_S = 30
DEFAULT_BATCH_SIZE = 1000
DEFAULT_REPLAY_THREADS = 50
DEFAULT_THEME = "light"
DEFAULT_RUNS_SHOWN = 5
DEFAULT_GAIN_DROP_WARN_PCT = 3.0
DEFAULT_MAX_DOWNLOAD_EVIO = 10
DEFAULT_QUEUE_CAP_EVIO = 100
NEW_YORK_TZ = ZoneInfo("America/New_York")
LOG_MAX_BLOCKS = 200
LOG_BUFFER_LIMIT = 250
LOG_FLUSH_MS = 200
MAINTENANCE_INTERVAL_MS = 60 * 1000
QUIET_PROCESS_PREFIXES = (
    "Already replayed:",
    "Already queued:",
    "Safety hold:",
)

QUANTITIES = [
    ("gain_norm_W", "Gain / Ref Gain"),
    ("gain_change", "Change"),
    ("gain_long_change", "Long Change"),
    ("gain_run_change", "Run-to-Run Change"),
    ("gain_all_run_change", "Current vs All Runs"),
]
REF_CHOICES = ["Ref1", "Ref2", "Ref3", "Avg"]
REF_COLORS = [
    QColor("#2997ff"),
    QColor("#30d158"),
    QColor("#ff9f0a"),
    QColor("#bf5af2"),
]
DEFAULT_REF_RATIO_BAD_PCT = 20
CHART_BG = QColor("#ffffff")
CHART_TEXT = QColor("#1f2328")
CHART_TEXT_DIM = QColor("#57606a")
CHART_GRID = QColor("#d8dee4")
CHART_BORDER = QColor("#8c959f")
CHART_LEGEND_BG = QColor(255, 255, 255, 230)
CHANGE_STOPS = [
    (0.00, (5, 10, 24)),
    (0.25, (29, 78, 216)),
    (0.50, (34, 197, 94)),
    (0.75, (245, 158, 11)),
    (1.00, (220, 38, 38)),
]
CHANGE_DEFAULT_RANGE = 0.10

GAIN_CORR_FILE_RE = re.compile(r"^prad_(\d{6})_gain_corr(?:\.|。)root$")

set_theme(DEFAULT_THEME)

_CPU_LAYOUT: Optional[Tuple[Optional[int], List[int]]] = None


def _format_cpu_list(cpus: List[int]) -> str:
    if not cpus:
        return ""
    ranges = []
    start = prev = cpus[0]
    for cpu in cpus[1:]:
        if cpu == prev + 1:
            prev = cpu
            continue
        ranges.append(f"{start}-{prev}" if start != prev else str(start))
        start = prev = cpu
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ",".join(ranges)


def configure_cpu_layout() -> Tuple[Optional[int], List[int]]:
    global _CPU_LAYOUT
    if _CPU_LAYOUT is not None:
        return _CPU_LAYOUT

    try:
        available = sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        available = list(range(os.cpu_count() or 1))

    if len(available) <= 1 or not hasattr(os, "sched_setaffinity"):
        _CPU_LAYOUT = (available[0] if available else None, available)
        return _CPU_LAYOUT

    viewer_cpu = available[0]
    replay_cpus = available[1:]
    try:
        os.sched_setaffinity(0, {viewer_cpu})
    except OSError:
        _CPU_LAYOUT = (None, available)
        return _CPU_LAYOUT

    _CPU_LAYOUT = (viewer_cpu, replay_cpus)
    return _CPU_LAYOUT


def start_bash_process(proc: QProcess, bash: str, cpu_list: str = ""):
    if cpu_list:
        proc.start("taskset", ["-c", cpu_list, "bash", "-lc", bash])
    else:
        proc.start("bash", ["-lc", bash])


@dataclass
class GainData:
    run_number: int
    path: Path
    event_start: object
    event_end: object
    unix_time: object
    gain_w: object
    gain_w_ref: object
    gain_corr_w: object
    fit_mean_w_lms: object
    ref_ratio: object

    @property
    def nbatches(self) -> int:
        return int(len(self.event_start))


def find_update_tool() -> str:
    # Prefer the executable built from the same source tree as this GUI.
    # A PATH installation may be older and lack newer options such as -a.
    candidates = [
        str((SCRIPT_DIR / ".." / "build" / "bin" / "prad2ana_online_gain_monitor_base").resolve()),
        str((SCRIPT_DIR / ".." / "build-clang" / "bin" / "prad2ana_online_gain_monitor_base").resolve()),
        shutil.which("prad2ana_online_gain_monitor_base"),
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return "prad2ana_online_gain_monitor_base"


def run_tag(run_number: int) -> str:
    return f"prad_{run_number:06d}"


def run_storage_paths(run_number: int, storage_base: str) -> Tuple[Path, Path, Path]:
    base = Path(storage_base).expanduser()
    tag = run_tag(run_number)
    evio_dir = base / "evio" / tag
    work_dir = base / "lms" / tag
    out_root = base / "gain" / tag / f"{tag}_gain_corr.root"
    return evio_dir, work_dir, out_root


def load_gain_root(path: Path, run_number: int) -> GainData:
    if not HAS_UPROOT:
        raise RuntimeError("uproot/numpy not installed. Install with: pip install uproot numpy")
    executor = TrivialExecutor()
    with uproot.open(
        str(path),
        handler=MemmapSource,
        decompression_executor=executor,
        interpretation_executor=executor,
        object_cache=None,
        array_cache=None,
    ) as f:
        if "gain_corr" not in f:
            raise RuntimeError(f"No gain_corr tree in {path}")
        t = f["gain_corr"]
        event_start = t["event_num_start"].array(library="np")
        event_end = t["event_num_end"].array(library="np")
        unix_time = (
            t["unix_time"].array(library="np")
            if "unix_time" in t.keys()
            else np.zeros(len(event_start), dtype=np.uint32)
        )
        if len(event_start) == 0:
            raise RuntimeError(f"No batch entries in gain_corr tree: {path}")
        gain_w = t["gain_W"].array(library="np")
        gain_corr_w = t["gain_corr_W"].array(library="np")
        gain_w_ref = (
            t["gain_W_ref"].array(library="np")
            if "gain_W_ref" in t.keys()
            else np.where(gain_corr_w > 0.0, gain_w * gain_corr_w, np.nan)
        )
        return GainData(
            run_number=run_number,
            path=path,
            event_start=event_start,
            event_end=event_end,
            unix_time=unix_time,
            gain_w=gain_w,
            gain_w_ref=gain_w_ref,
            gain_corr_w=gain_corr_w,
            fit_mean_w_lms=t["fit_mean_W_lms"].array(library="np"),
            ref_ratio=t["refPMT_ratio"].array(library="np"),
        )


def _parse_lms_dat(path: Path) -> Optional[np.ndarray]:
    """Parse prad_XXXXXX_LMS.dat; return ref gain array [1156, 3].

    The .dat file columns are:
      0:Name  1:lms_peak  2:lms_sigma  3:lms_chi2/ndf  4:g1  5:g2  6:g3

    For W module rows, columns 4-6 (g1/g2/g3) store the reference gain values
    corresponding to Ref PMT 1/2/3, in the same units as gain_W in the ROOT file.
    Dividing gain_W by these values gives a ratio ~1 for an unchanged gain.
    """
    entries: Dict[str, List[float]] = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 7 or parts[0] == "Name":
                    continue
                try:
                    entries[parts[0]] = [float(parts[4]), float(parts[5]), float(parts[6])]
                except (ValueError, IndexError):
                    continue
    except OSError:
        return None
    result = np.full((1156, 3), np.nan, dtype=float)
    for i in range(1156):
        vals = entries.get(f"W{i + 1}")
        if vals is None:
            continue
        for j in range(3):
            if vals[j] > 0.0:
                result[i, j] = vals[j]
    return result


def build_gain_summary_groups(path: Path) -> Dict[str, List[int]]:
    """Return W-array indices for the absorber region and two open rings."""
    with open(path) as stream:
        entries = json.load(stream)
    grouped = {"under absorber": [], "open layer 1": [], "open layer 2": []}
    positioned = []
    for entry in entries:
        name = entry.get("n", "")
        geo = entry.get("geo", {})
        if not name.startswith("W") or "row" not in geo or "col" not in geo:
            continue
        try:
            index = int(name[1:]) - 1
        except ValueError:
            continue
        row, col = int(geo["row"]), int(geo["col"])
        # The absorber covers rows/columns 15..18. Each open layer is the
        # next Chebyshev ring around that square.
        distance = max(max(15 - row, 0, row - 18), max(15 - col, 0, col - 18))
        if distance == 0:
            group = "under absorber"
        elif distance == 1:
            group = "open layer 1"
        elif distance == 2:
            group = "open layer 2"
        else:
            continue
        positioned.append((group, row, col, index))
    for group, _, _, index in sorted(positioned, key=lambda item: (item[0], item[1], item[2])):
        grouped[group].append(index)
    return grouped


class GainRootLoader(QObject):
    finished = pyqtSignal(object, object)

    def __init__(self, files: List[Tuple[int, Path]], worker_cpus: List[int]):
        super().__init__()
        self._files = files
        self._worker_cpus = worker_cpus

    @pyqtSlot()
    def run(self):
        if self._worker_cpus and hasattr(os, "sched_setaffinity"):
            try:
                os.sched_setaffinity(0, set(self._worker_cpus))
            except OSError:
                pass
        loaded: Dict[int, GainData] = {}
        errors: List[Tuple[int, Path, str]] = []
        for run, path in self._files:
            if QThread.currentThread().isInterruptionRequested():
                break
            try:
                loaded[run] = load_gain_root(path, run)
            except Exception as exc:
                errors.append((run, path, str(exc)))
        self.finished.emit(loaded, errors)


def positive_avg(vals) -> float:
    good = [
        float(v) for v in vals
        if math.isfinite(float(v)) and float(v) > 0.0 and float(v) != 1.0
    ]
    return sum(good) / len(good) if good else math.nan


def finite_avg(vals) -> float:
    good = [float(v) for v in vals if math.isfinite(float(v))]
    return sum(good) / len(good) if good else math.nan


def safe_divide(numerator, denominator, valid=None):
    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    if valid is None:
        valid = np.isfinite(num) & np.isfinite(den) & (den != 0.0)
    return np.divide(
        num,
        den,
        out=np.full(np.broadcast_shapes(num.shape, den.shape), np.nan, dtype=float),
        where=valid,
    )


def masked_mean(samples: List[object], masks: List[List[bool]]):
    if not samples:
        return None
    arr = np.asarray(samples, dtype=float)
    mask = np.asarray(masks, dtype=bool)
    valid = mask[:, None, :] & np.isfinite(arr) & (arr > 0.0)
    count = np.sum(valid, axis=0)
    total = np.sum(np.where(valid, arr, 0.0), axis=0)
    return safe_divide(total, count, count > 0)


def ref_ratio_ok(value: float, center: float, bad_rel: float) -> bool:
    return (
        math.isfinite(value) and value > 0.0
        and math.isfinite(center) and center > 0.0
        and abs(value - center) / center <= bad_rel
    )


def nice_ticks(y_min: float, y_max: float, target: int = 5) -> Tuple[List[float], float, float, int]:
    if not (math.isfinite(y_min) and math.isfinite(y_max)):
        return [0.0, 1.0], 0.0, 1.0, 0
    if y_min == y_max:
        d = abs(y_min) * 0.05 if y_min else 1.0
        y_min -= d
        y_max += d
    raw_step = abs(y_max - y_min) / max(1, target - 1)
    if raw_step <= 0.0 or not math.isfinite(raw_step):
        raw_step = 1.0
    mag = 10 ** math.floor(math.log10(raw_step))
    norm = raw_step / mag
    if norm <= 1.0:
        step = 1.0 * mag
    elif norm <= 2.0:
        step = 2.0 * mag
    elif norm <= 2.5:
        step = 2.5 * mag
    elif norm <= 5.0:
        step = 5.0 * mag
    else:
        step = 10.0 * mag
    tick_min = math.floor(y_min / step) * step
    tick_max = math.ceil(y_max / step) * step
    ticks: List[float] = []
    cur = tick_min
    limit = 0
    while cur <= tick_max + step * 0.5 and limit < 40:
        ticks.append(0.0 if abs(cur) < step * 1e-6 else cur)
        cur += step
        limit += 1
    step_text = f"{step:.10f}".rstrip("0").rstrip(".")
    decimals = len(step_text.split(".", 1)[1]) if "." in step_text else 0
    return ticks, tick_min, tick_max, decimals


def tick_label(value: float, decimals: int) -> str:
    if decimals <= 0:
        return f"{int(round(value))}"
    if abs(value) < 0.5 * 10 ** (-decimals):
        value = 0.0
    return f"{value:.{decimals}f}"


def decimate_points(points, max_points: int):
    """Bound paint work while retaining the first and last samples."""
    if len(points) <= max_points:
        return points
    step = max(1, math.ceil((len(points) - 1) / max(1, max_points - 1)))
    reduced = points[::step]
    if reduced[-1] != points[-1]:
        reduced.append(points[-1])
    return reduced


class GainMapWidget(HyCalMapWidget):
    def __init__(self, parent=None):
        super().__init__(parent, enable_zoom_pan=True, include_lms=False)
        self._selected: Optional[str] = None
        self._change_mode = False

    def set_selected(self, name: Optional[str]):
        self._selected = name
        self.update()

    def set_change_mode(self, enabled: bool):
        self._change_mode = enabled
        self.update()

    def palette_stops(self):
        return CHANGE_STOPS if self._change_mode else super().palette_stops()

    def _colorbar_center_text(self) -> str:
        return "down  stable  up" if self._change_mode else super()._colorbar_center_text()

    def _fmt_value(self, v: float) -> str:
        if not math.isfinite(v):
            return "nan"
        if self._change_mode:
            return f"{v * 100:+.2f}%"
        return f"{v:.4f}"

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        if v is None or not math.isfinite(v):
            return name
        return f"{name}: {v:.6f}"

    def _paint_overlays(self, p: QPainter, w: int, h: int):
        super()._paint_overlays(p, w, h)
        if self._selected and self._selected in self._rects:
            p.setPen(QPen(QColor(THEME.SELECT_BORDER), 2.5))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._selected])


class BatchChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._title = "No module selected"
        self._x_label = "batch id"
        self._x_mode = "batch"
        self._x: List[float] = []
        self._series: List[Tuple[str, List[float], QColor, bool]] = []
        self._run_spans: List[Tuple[int, float, float]] = []
        self._default_y_range: Optional[Tuple[float, float]] = None
        self.setMinimumSize(220, 180)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def set_data(self, title: str, x: List[float], series: List[Tuple[str, List[float], QColor, bool]],
                 run_spans: Optional[List[Tuple[int, float, float]]] = None,
                 default_y_range: Optional[Tuple[float, float]] = None,
                 x_label: str = "batch id", x_mode: str = "batch"):
        self._title = title
        self._x_label = x_label
        self._x_mode = x_mode
        self._x = x
        self._series = series
        self._run_spans = run_spans or []
        self._default_y_range = default_y_range
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, CHART_BG)

        pad_l = 64 if w < 420 else 74
        pad_r, pad_t = 12, 46
        pad_b = 72 if self._x_mode == "unix" else 50
        plot_w = max(1, w - pad_l - pad_r)
        plot_h = max(1, h - pad_t - pad_b)

        p.setPen(CHART_TEXT)
        p.setFont(QFont("Consolas", 11, QFont.Weight.Bold))
        title = p.fontMetrics().elidedText(self._title, Qt.TextElideMode.ElideRight, max(40, w - 20))
        p.drawText(10, 22, title)

        finite_vals: List[float] = []
        for _, values, _, _ in self._series:
            finite_vals.extend(
                v for x, v in zip(self._x, values)
                if math.isfinite(x) and math.isfinite(v)
            )
        finite_x = [x for x in self._x if math.isfinite(x)]
        if not finite_x or not finite_vals:
            p.setPen(CHART_TEXT_DIM)
            p.setFont(QFont("Consolas", 12))
            message = "No Unix time data" if self._x_mode in {"unix", "minutes"} else "No batch data"
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, message)
            return

        x_min, x_max = min(finite_x), max(finite_x)
        if self._x_mode != "unix":
            x_min = min(0.0, x_min)
        if x_min == x_max:
            x_min -= 1
            x_max += 1
        x_tick_target = max(10, min(20, plot_w // 45)) if self._x_mode == "unix" else 24
        x_ticks, x_min, x_max, x_decimals = nice_ticks(x_min, x_max, x_tick_target)
        x_ticks = [t for t in x_ticks if x_min <= t <= x_max]

        data_y_min, data_y_max = min(finite_vals), max(finite_vals)
        if self._default_y_range is not None:
            y_min = min(self._default_y_range[0], data_y_min)
            y_max = max(self._default_y_range[1], data_y_max)
            if data_y_min < self._default_y_range[0]:
                y_min -= (self._default_y_range[0] - data_y_min) * 0.08
            if data_y_max > self._default_y_range[1]:
                y_max += (data_y_max - self._default_y_range[1]) * 0.08
            y_ticks_all, _, _, y_decimals = nice_ticks(y_min, y_max, 24)
            y_ticks = [t for t in y_ticks_all if y_min <= t <= y_max]
        else:
            y_min, y_max = data_y_min, data_y_max
            if y_min == y_max:
                d = abs(y_min) * 0.05 if y_min != 0 else 1.0
                y_min -= d
                y_max += d
            else:
                d = (y_max - y_min) * 0.12
                y_min -= d
                y_max += d
            y_ticks, y_min, y_max, y_decimals = nice_ticks(y_min, y_max, 24)

        def sx(x):
            return pad_l + (x - x_min) / (x_max - x_min) * plot_w

        def sy(y):
            return pad_t + (y_max - y) / (y_max - y_min) * plot_h

        p.setPen(QPen(CHART_GRID, 1))
        for tick in y_ticks:
            yy = sy(tick)
            p.drawLine(pad_l, int(yy), pad_l + plot_w, int(yy))
        for tick in x_ticks:
            xx = sx(tick)
            p.drawLine(int(xx), pad_t, int(xx), pad_t + plot_h)
        p.setPen(QPen(CHART_BORDER, 1))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(pad_l, pad_t, plot_w, plot_h)

        if self._run_spans:
            dash_pen = QPen(CHART_TEXT_DIM, 1)
            dash_pen.setStyle(Qt.PenStyle.DashLine)
            p.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
            for run, start, end in self._run_spans:
                for bx in (start, end):
                    if x_min <= bx <= x_max:
                        xx = sx(bx)
                        p.setPen(dash_pen)
                        p.drawLine(int(xx), pad_t, int(xx), pad_t + plot_h)
                mid = (start + end) * 0.5
                if x_min <= mid <= x_max:
                    p.setPen(CHART_TEXT_DIM)
                    label = f"{run:06d}"
                    p.drawText(int(sx(mid) - 28), pad_t + 18, label)

        p.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
        p.setPen(CHART_TEXT)
        for tick in y_ticks:
            yy = sy(tick)
            p.drawText(4, int(yy + 4), tick_label(tick, y_decimals))
        for tick in x_ticks:
            xx = sx(tick)
            if self._x_mode == "unix":
                p.setFont(QFont("Consolas", 7, QFont.Weight.Bold))
                label = datetime.fromtimestamp(tick, NEW_YORK_TZ).strftime("%m-%d\n%H:%M")
                label_x = max(pad_l, min(int(xx - 25), pad_l + plot_w - 50))
                p.drawText(
                    label_x, pad_t + plot_h + 4, 50, 30,
                    Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop,
                    label,
                )
            else:
                label = tick_label(tick, x_decimals)
                tw = p.fontMetrics().horizontalAdvance(label)
                p.drawText(int(xx - tw / 2), pad_t + plot_h + 16, label)
        p.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
        xlabel = p.fontMetrics().elidedText(self._x_label, Qt.TextElideMode.ElideRight, max(40, plot_w))
        tw = p.fontMetrics().horizontalAdvance(xlabel)
        p.drawText(int(pad_l + plot_w / 2 - tw / 2), h - 10, xlabel)

        legend_items = []
        for label, values, color, draw_line in self._series:
            pts = [
                (sx(x), sy(v))
                for x, v in zip(self._x, values)
                if math.isfinite(x) and math.isfinite(v) and y_min <= v <= y_max
            ]
            pts = decimate_points(pts, max(100, plot_w * 2))
            if len(pts) < 1:
                continue
            legend_items.append((label, color, draw_line, pts))
            if draw_line:
                p.setPen(QPen(color, 2.0))
                for a, b in zip(pts, pts[1:]):
                    p.drawLine(int(a[0]), int(a[1]), int(b[0]), int(b[1]))
            p.setBrush(color)
            for x, y in pts:
                p.drawEllipse(int(x - 2), int(y - 2), 4, 4)

        if legend_items:
            p.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
            fm = p.fontMetrics()
            item_w = [fm.horizontalAdvance(label) + 22 for label, _, _, _ in legend_items]
            total_w = min(sum(item_w) + 12, max(80, plot_w - 12))
            lx = pad_l + plot_w - total_w - 6
            # Keep the Ref1/Ref2/Ref3 legend outside the plotting rectangle,
            # in the chart header's upper-right corner.
            ly = 5
            lh = 20
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(CHART_LEGEND_BG)
            p.drawRoundedRect(lx, ly, total_w, lh, 4, 4)
            p.setPen(QPen(CHART_BORDER, 1))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRoundedRect(lx, ly, total_w, lh, 4, 4)
            xcur = lx + 8
            for (label, color, draw_line, _), width in zip(legend_items, item_w):
                if xcur + width > lx + total_w:
                    break
                p.setPen(QPen(color, 2))
                cy = ly + lh // 2
                if draw_line:
                    p.drawLine(xcur, cy, xcur + 12, cy)
                p.setBrush(color)
                p.drawEllipse(int(xcur + 5), int(cy - 3), 6, 6)
                p.setPen(color)
                p.drawText(int(xcur + 18), int(ly + 14), label)
                xcur += width


class RefRatioChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._x: List[float] = []
        self._x_mode = "batch"
        self._ratios: List[List[float]] = [[], [], []]
        self._ok: List[List[bool]] = [[], [], []]
        self._run_spans: List[Tuple[int, float, float]] = []
        self._centers: List[float] = [math.nan, math.nan, math.nan]
        self.setMinimumSize(220, 240)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def set_data(self, x: List[float], ratios: List[List[float]],
                 ok: List[List[bool]], run_spans: List[Tuple[int, float, float]],
                 centers: Optional[List[float]] = None, x_mode: str = "batch"):
        self._x = x
        self._x_mode = x_mode
        self._ratios = ratios
        self._ok = ok
        self._run_spans = run_spans
        self._centers = centers or [math.nan, math.nan, math.nan]
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, CHART_BG)

        p.setPen(CHART_TEXT)
        p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
        p.drawText(8, 18, "Ref PMT LMS/Alpha ratio")

        finite = [
            v for series in self._ratios for x, v in zip(self._x, series)
            if math.isfinite(x) and math.isfinite(v) and v > 0.0
        ]
        finite_x = [x for x in self._x if math.isfinite(x)]
        if not finite_x or not finite:
            p.setPen(CHART_TEXT_DIM)
            p.setFont(QFont("Consolas", 10))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No ref ratio data")
            return

        x_min, x_max = min(finite_x), max(finite_x)
        if self._x_mode != "unix":
            x_min = min(0.0, x_min)
        if x_min == x_max:
            x_min -= 1.0
            x_max += 1.0
        pad_l, pad_r, pad_t, pad_b = 68, 10, 28, 24
        gap = 28
        plot_w = max(1, w - pad_l - pad_r)
        panel_w = max(1, (plot_w - 2 * gap) / 3)
        plot_h = max(1, h - pad_t - pad_b)
        y_ranges: List[Tuple[float, float]] = []
        y_ticks_by_panel: List[List[float]] = []
        y_decimals: List[int] = []
        for panel in range(3):
            vals = [v for v in self._ratios[panel] if math.isfinite(v) and v > 0.0]
            if not vals:
                y_ranges.append((0.0, 1.0))
                y_ticks_by_panel.append([0.0, 0.5, 1.0])
                y_decimals.append(1)
                continue
            y_min, y_max = min(vals), max(vals)
            center = self._centers[panel] if panel < len(self._centers) else math.nan
            if math.isfinite(center) and center > 0.0:
                y_min = min(y_min, center)
                y_max = max(y_max, center)
            if y_min == y_max:
                d = abs(y_min) * 0.02 if y_min else 1.0
                y_min -= d
                y_max += d
            else:
                d = (y_max - y_min) * 0.15
                y_min -= d
                y_max += d
            ticks, y_min, y_max, decimals = nice_ticks(y_min, y_max, 20)
            y_ranges.append((y_min, y_max))
            y_ticks_by_panel.append(ticks)
            y_decimals.append(decimals)

        def sx(panel, x):
            base = pad_l + panel * (panel_w + gap)
            return base + (x - x_min) / (x_max - x_min) * panel_w

        def sy(panel, y):
            y_min, y_max = y_ranges[panel]
            return pad_t + (y_max - y) / (y_max - y_min) * plot_h

        dash_pen = QPen(CHART_TEXT_DIM, 1)
        dash_pen.setStyle(Qt.PenStyle.DashLine)
        bad_pen = QPen(QColor(THEME.DANGER), 1.6)
        center_pen = QPen(CHART_TEXT_DIM, 1)
        center_pen.setStyle(Qt.PenStyle.DotLine)
        for panel in range(3):
            x0 = pad_l + panel * (panel_w + gap)
            p.setPen(QPen(CHART_GRID, 1))
            for i in range(21):
                xx = x0 + panel_w * i / 20
                p.drawLine(int(xx), pad_t, int(xx), pad_t + plot_h)
            p.setPen(QPen(CHART_BORDER, 1))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(int(x0), pad_t, int(panel_w), plot_h)
            if panel > 0:
                p.drawLine(int(x0 - gap / 2), pad_t, int(x0 - gap / 2), pad_t + plot_h)
            p.setPen(REF_COLORS[panel])
            p.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
            p.drawText(int(x0 + 6), pad_t + 14, REF_CHOICES[panel])
            p.setFont(QFont("Consolas", 8, QFont.Weight.Bold))
            for tick in y_ticks_by_panel[panel]:
                yy = sy(panel, tick)
                p.setPen(QPen(CHART_GRID, 1))
                p.drawLine(int(x0), int(yy), int(x0 + panel_w), int(yy))
                p.setPen(CHART_TEXT)
                p.drawText(int(x0 + 4), int(yy - 3), tick_label(tick, y_decimals[panel]))
            center = self._centers[panel] if panel < len(self._centers) else math.nan
            if math.isfinite(center) and center > 0.0:
                yy = sy(panel, center)
                p.setPen(center_pen)
                p.drawLine(int(x0), int(yy), int(x0 + panel_w), int(yy))
                p.setPen(CHART_TEXT_DIM)
                p.drawText(int(x0 + panel_w - 58), pad_t + plot_h - 4, f"avg {center:.3g}")

            pts = []
            for x, y, ok in zip(self._x, self._ratios[panel], self._ok[panel]):
                if not (math.isfinite(x) and math.isfinite(y) and y > 0.0):
                    continue
                px, py = sx(panel, x), sy(panel, y)
                pts.append((px, py, ok))
            pts = decimate_points(pts, max(80, int(panel_w) * 2))
            p.setPen(QPen(REF_COLORS[panel], 1.3))
            last_good = None
            for px, py, ok in pts:
                if ok and last_good is not None:
                    p.drawLine(int(last_good[0]), int(last_good[1]), int(px), int(py))
                last_good = (px, py) if ok else None
            for px, py, ok in pts:
                if ok:
                    p.setPen(QPen(REF_COLORS[panel], 1))
                    p.setBrush(REF_COLORS[panel])
                    p.drawEllipse(int(px - 2), int(py - 2), 4, 4)
                else:
                    p.setPen(bad_pen)
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawLine(int(px - 4), int(py - 4), int(px + 4), int(py + 4))
                    p.drawLine(int(px - 4), int(py + 4), int(px + 4), int(py - 4))

            p.setPen(dash_pen)
            for _, start, end in self._run_spans:
                for bx in (start, end):
                    if x_min <= bx <= x_max:
                        xx = sx(panel, bx)
                        p.drawLine(int(xx), pad_t, int(xx), pad_t + plot_h)


class OnlineGainMonitor(QMainWindow):
    def __init__(self):
        super().__init__()
        self._modules = load_modules(MODULES_JSON)
        self._gain_summary_groups = build_gain_summary_groups(MODULES_JSON)
        self._selected_module = "W1"
        self._x_axis_mode = "batch"
        self._data: Optional[GainData] = None
        self._runs_data: Dict[int, GainData] = {}
        self._known_runs: set[int] = set()
        self._opened_output_folder: Optional[Path] = None
        self._folder_gain_files: Dict[int, Path] = {}
        self._updating_run_range = False
        self._ref_gain_file: Optional[Path] = None
        self._ref_file_gain_array: Optional[np.ndarray] = None
        self._ref_gain_cache: Dict[int, Optional[np.ndarray]] = {}
        self._last_gain_drop_warning_key: Optional[Tuple[int, int, float]] = None
        self._replay_queue: List[Tuple[int, Path]] = []
        self._queued_snapshots: set[Path] = set()
        self._active_replay_snapshots: List[Path] = []
        self._active_replay_file_count = 0
        self._reanalyze_pending = False
        self._scan_copied_count = 0
        self._scan_floor_run: Optional[int] = None
        self._last_pending_evio_count = 0
        self._root_signatures: Dict[int, Tuple[str, int, int]] = {}
        self._root_load_thread: Optional[QThread] = None
        self._root_load_worker: Optional[GainRootLoader] = None
        self._pending_root_load: Optional[Tuple[List[Tuple[int, Path]], str, int]] = None
        self._root_load_context = None
        self._view_context = None
        self._process_output_remainders: Dict[Tuple[int, bool], str] = {}
        self._download_process = QProcess(self)
        self._download_process.readyReadStandardOutput.connect(self._on_stdout)
        self._download_process.readyReadStandardError.connect(self._on_stderr)
        self._download_process.finished.connect(self._on_download_finished)
        self._replay_process = QProcess(self)
        self._replay_process.readyReadStandardOutput.connect(self._on_stdout)
        self._replay_process.readyReadStandardError.connect(self._on_stderr)
        self._replay_process.finished.connect(self._on_replay_finished)
        self._reanalyze_process = QProcess(self)
        self._reanalyze_process.readyReadStandardOutput.connect(self._on_stdout)
        self._reanalyze_process.readyReadStandardError.connect(self._on_stderr)
        self._reanalyze_process.finished.connect(self._on_reanalyze_finished)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._scan_once)
        self._batch_reanalyze_timer = QTimer(self)
        self._batch_reanalyze_timer.setSingleShot(True)
        self._batch_reanalyze_timer.setInterval(600)
        self._batch_reanalyze_timer.timeout.connect(self._reanalyze_lms)

        self._build_ui()
        self._log_buffer: List[Tuple[str, str]] = []
        self._log_flush_timer = QTimer(self)
        self._log_flush_timer.setSingleShot(True)
        self._log_flush_timer.timeout.connect(self._flush_log)
        self._maintenance_timer = QTimer(self)
        self._maintenance_timer.timeout.connect(self._maintenance)
        self._maintenance_timer.start(MAINTENANCE_INTERVAL_MS)
        self._viewer_cpu, self._replay_cpus = configure_cpu_layout()
        self._download_cpu_list = _format_cpu_list(self._replay_cpus[:1]) if self._replay_cpus else ""
        self._replay_worker_cpus = self._replay_cpus[1:] if len(self._replay_cpus) > 1 else self._replay_cpus
        if self._replay_cpus and self._viewer_cpu is not None:
            self._append(
                f"CPU isolation: viewer on CPU {self._viewer_cpu}; download on CPU(s) {self._download_cpu_list or 'all'}; replay can use {len(self._replay_worker_cpus)} CPU(s).",
                "ok",
            )
        else:
            self._append("CPU isolation unavailable; using thread cap only.", "warn")
        self._map.set_modules(self._modules)
        self._map.set_selected(self._selected_module)
        self._map.moduleClicked.connect(self._on_module_clicked)
        self._update_paths_from_run()

    def _build_ui(self):
        self.setWindowTitle("Online LMS Gain Monitor")
        self.resize(1400, 850)
        self.setMinimumSize(760, 520)
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        top_widget = QWidget()
        top_rows = QVBoxLayout(top_widget)
        top_rows.setContentsMargins(0, 0, 0, 0)
        top_rows.setSpacing(4)
        top_row1 = QWidget()
        top_row2 = QWidget()
        top = QHBoxLayout(top_row1)
        top.setContentsMargins(0, 0, 0, 0)
        top.setSpacing(5)
        top2 = QHBoxLayout(top_row2)
        top2.setContentsMargins(0, 0, 0, 0)
        top2.setSpacing(5)
        title = QLabel("ONLINE LMS GAIN MONITOR")
        title.setFont(QFont("Consolas", 14, QFont.Weight.Bold))
        title.setStyleSheet(themed(f"color:{THEME.ACCENT};"))
        top.addWidget(title)

        top.addWidget(self._label("Start run"))
        self._run_edit = QLineEdit("")
        self._run_edit.setPlaceholderText("24929")
        self._run_edit.setMinimumWidth(150)
        self._run_edit.editingFinished.connect(self._update_paths_from_run)
        top.addWidget(self._run_edit)

        self._start_btn = self._button("Start", THEME.SUCCESS, self._start)
        self._stop_btn = self._button("Stop", THEME.DANGER, self._stop)
        self._stop_btn.setEnabled(False)
        self._scan_btn = self._button("Scan Now", THEME.ACCENT, self._scan_once)
        self._open_folder_btn = self._button("Open Folder", THEME.WARN, self._open_output_folder)
        top.addWidget(self._start_btn)
        top.addWidget(self._stop_btn)
        top.addWidget(self._scan_btn)
        top.addWidget(self._open_folder_btn)

        top.addWidget(self._label("every"))
        self._interval = QSpinBox()
        self._interval.setRange(5, 3600)
        self._interval.setValue(DEFAULT_INTERVAL_S)
        self._interval.setSuffix(" s")
        self._interval.setFixedWidth(80)
        top.addWidget(self._interval)

        top.addWidget(self._label("batch"))
        self._batch = QSpinBox()
        self._batch.setRange(100, 100000)
        self._batch.setValue(DEFAULT_BATCH_SIZE)
        self._batch.setSingleStep(100)
        self._batch.setFixedWidth(95)
        self._batch.valueChanged.connect(self._schedule_batch_reanalysis)
        top.addWidget(self._batch)
        self._reanalyze_btn = self._button("Re-analyze", THEME.ACCENT, self._reanalyze_lms)
        top.addWidget(self._reanalyze_btn)

        top.addWidget(self._label("evio skip"))
        self._evio_skip = QSpinBox()
        self._evio_skip.setRange(0, 999)
        self._evio_skip.setValue(0)
        self._evio_skip.setFixedWidth(70)
        top.addWidget(self._evio_skip)

        top.addWidget(self._label("max evio"))
        self._max_download_evio = QSpinBox()
        self._max_download_evio.setRange(1, 100000)
        self._max_download_evio.setValue(DEFAULT_MAX_DOWNLOAD_EVIO)
        self._max_download_evio.setFixedWidth(75)
        top.addWidget(self._max_download_evio)

        top.addWidget(self._label("queue cap"))
        self._queue_cap_evio = QSpinBox()
        self._queue_cap_evio.setRange(1, 1000000)
        self._queue_cap_evio.setValue(DEFAULT_QUEUE_CAP_EVIO)
        self._queue_cap_evio.setFixedWidth(85)
        top.addWidget(self._queue_cap_evio)

        top.addWidget(self._label("runs shown"))
        self._runs_shown = QSpinBox()
        self._runs_shown.setRange(1, 999)
        self._runs_shown.setValue(DEFAULT_RUNS_SHOWN)
        self._runs_shown.setFixedWidth(70)
        self._runs_shown.valueChanged.connect(self._runs_shown_changed)
        top.addWidget(self._runs_shown)
        self._valid_runs_count = QLabel("0 loaded")
        self._valid_runs_count.setFont(QFont("Consolas", 9))
        self._valid_runs_count.setStyleSheet(themed(f"color:{THEME.TEXT_DIM};"))
        top.addWidget(self._valid_runs_count)

        top.addWidget(self._label("threads"))
        self._threads = QSpinBox()
        self._threads.setRange(1, 1024)
        self._threads.setValue(DEFAULT_REPLAY_THREADS)
        self._threads.setFixedWidth(60)
        top.addWidget(self._threads)

        top2.addWidget(self._label("ref run"))
        self._ref_run = QSpinBox()
        self._ref_run.setRange(-1, 999999)
        self._ref_run.setValue(-1)
        self._ref_run.setSpecialValueText("auto")
        self._ref_run.setFixedWidth(90)
        self._ref_run.valueChanged.connect(self._on_ref_run_changed)
        top2.addWidget(self._ref_run)

        self._ref_file_btn = self._button("Ref File", THEME.ACCENT, self._browse_ref_gain_file)
        top2.addWidget(self._ref_file_btn)

        top2.addWidget(self._label("start run"))
        self._display_run_start = QComboBox()
        self._display_run_start.setMinimumWidth(92)
        self._display_run_start.setEnabled(False)
        self._display_run_start.currentIndexChanged.connect(self._display_run_range_changed)
        top2.addWidget(self._display_run_start)

        top2.addWidget(self._label("end run"))
        self._display_run_end = QComboBox()
        self._display_run_end.setMinimumWidth(92)
        self._display_run_end.setEnabled(False)
        self._display_run_end.currentIndexChanged.connect(self._display_run_range_changed)
        top2.addWidget(self._display_run_end)

        top2.addWidget(self._label("quantity"))
        self._quantity = QComboBox()
        self._quantity.addItems([label for _, label in QUANTITIES])
        self._quantity.currentIndexChanged.connect(self._invalidate_and_refresh_views)
        top2.addWidget(self._quantity)

        top2.addWidget(self._label("map range"))
        self._change_range = QDoubleSpinBox()
        self._change_range.setRange(0.1, 100.0)
        self._change_range.setDecimals(1)
        self._change_range.setSingleStep(0.5)
        self._change_range.setValue(CHANGE_DEFAULT_RANGE * 100.0)
        self._change_range.setSuffix(" %")
        self._change_range.setFixedWidth(85)
        self._change_range.valueChanged.connect(self._refresh_views)
        top2.addWidget(self._change_range)

        top2.addWidget(self._label("drop warn"))
        self._drop_warn = QDoubleSpinBox()
        self._drop_warn.setRange(0.0, 100.0)
        self._drop_warn.setDecimals(1)
        self._drop_warn.setSingleStep(0.5)
        self._drop_warn.setValue(DEFAULT_GAIN_DROP_WARN_PCT)
        self._drop_warn.setSuffix(" %")
        self._drop_warn.setFixedWidth(85)
        top2.addWidget(self._drop_warn)

        top2.addWidget(self._label("ref"))
        self._ref = QComboBox()
        self._ref.addItems(REF_CHOICES)
        self._ref.setCurrentIndex(3)
        self._ref.currentIndexChanged.connect(self._refresh_views)
        top2.addWidget(self._ref)

        self._norm_first5 = QCheckBox("first 5 -> 1")
        self._norm_first5.setFont(QFont("Consolas", 10))
        self._norm_first5.setStyleSheet(themed(
            f"QCheckBox{{color:{THEME.TEXT_DIM};spacing:4px;}}"
            f"QCheckBox::indicator{{width:14px;height:14px;}}"
        ))
        self._norm_first5.stateChanged.connect(self._invalidate_and_refresh_views)
        top2.addWidget(self._norm_first5)

        top2.addWidget(self._label("ratio tol"))
        self._ref_ratio_tol = QSpinBox()
        self._ref_ratio_tol.setRange(0, 100)
        self._ref_ratio_tol.setValue(DEFAULT_REF_RATIO_BAD_PCT)
        self._ref_ratio_tol.setSuffix(" %")
        self._ref_ratio_tol.setFixedWidth(80)
        self._ref_ratio_tol.valueChanged.connect(self._invalidate_and_refresh_views)
        top2.addWidget(self._ref_ratio_tol)

        top2.addStretch()
        self._status = QLabel("Idle")
        self._status.setFont(QFont("Consolas", 10))
        self._status.setStyleSheet(themed(f"color:{THEME.TEXT_DIM};"))
        top2.addWidget(self._status)
        top_rows.addWidget(top_row1)
        top_rows.addWidget(top_row2)
        root.addWidget(self._hscroll(top_widget))

        dirs_widget = QWidget()
        dirs = QHBoxLayout(dirs_widget)
        dirs.setContentsMargins(0, 0, 0, 0)
        dirs.setSpacing(5)
        self._host = self._path_edit(REMOTE_HOST)
        self._remote_base = self._path_edit(REMOTE_BASE)
        self._storage_base = self._path_edit(STORAGE_BASE)
        self._work_dir = self._path_edit("")
        self._out_root = self._path_edit("")
        self._work_dir.setReadOnly(True)
        self._out_root.setReadOnly(True)
        dirs.addWidget(self._label("host"))
        dirs.addWidget(self._host, 1)
        dirs.addWidget(self._label("remote"))
        dirs.addWidget(self._remote_base, 2)
        dirs.addWidget(self._label("storage base"))
        dirs.addWidget(self._storage_base, 3)
        browse_base = self._button("...", THEME.TEXT, self._browse_base)
        browse_base.setFixedWidth(30)
        dirs.addWidget(browse_base)
        dirs.addWidget(self._label("preview work"))
        dirs.addWidget(self._work_dir, 3)
        dirs.addWidget(self._label("preview out"))
        dirs.addWidget(self._out_root, 3)
        root.addWidget(self._hscroll(dirs_widget))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        self._map = GainMapWidget()
        self._map.setMinimumWidth(220)
        self._chart = BatchChart()
        self._ratio_chart = RefRatioChart()
        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.document().setMaximumBlockCount(LOG_MAX_BLOCKS)
        self._log.setMaximumHeight(125)
        self._log.setStyleSheet(themed(
            f"QTextEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};font-family:Consolas;font-size:9pt;}}"))

        left = QWidget()
        left.setMinimumWidth(260)
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        left_lay.setSpacing(6)
        left_lay.addWidget(self._map, stretch=4)
        left_lay.addWidget(self._log, stretch=2)

        right = QWidget()
        right.setMinimumWidth(420)
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(6)
        self._run_avg_label = QLabel("Run avg: n/a")
        self._run_avg_label.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
        self._run_avg_label.setStyleSheet(themed(
            f"QLabel{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};padding:4px 6px;}}"
        ))
        self._run_avg_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        right_lay.addWidget(self._run_avg_label)

        self._gain_summary = QTextEdit()
        self._gain_summary.setReadOnly(True)
        self._gain_summary.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self._gain_summary.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._gain_summary.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._gain_summary.setFixedHeight(82)
        self._gain_summary.setFont(QFont("Consolas", 9))
        self._gain_summary.setPlainText(
            "under absorber: n/a\nopen layer 1: n/a\nopen layer 2: n/a"
        )
        self._gain_summary.setStyleSheet(themed(
            f"QTextEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};padding:3px 5px;}}"
        ))
        right_lay.addWidget(self._gain_summary)

        axis_bar = QWidget()
        axis_lay = QHBoxLayout(axis_bar)
        axis_lay.setContentsMargins(0, 0, 0, 0)
        axis_lay.setSpacing(4)
        axis_lay.addStretch()
        axis_lay.addWidget(self._label("X axis"))
        self._x_axis_group = QButtonGroup(self)
        self._x_axis_group.setExclusive(True)
        self._x_axis_buttons: Dict[str, QPushButton] = {}
        for mode, text in (("batch", "Batch ID"), ("unix", "Date/Time"), ("minutes", "Minutes")):
            btn = QPushButton(text)
            btn.setCheckable(True)
            btn.setChecked(mode == self._x_axis_mode)
            btn.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
            btn.setStyleSheet(themed(
                f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT_DIM};"
                f"border:1px solid {THEME.BORDER};border-radius:3px;padding:3px 8px;}}"
                f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}"
                f"QPushButton:checked{{background:{THEME.ACCENT};color:white;"
                f"border-color:{THEME.ACCENT};}}"
            ))
            btn.clicked.connect(lambda checked, m=mode: checked and self._set_x_axis_mode(m))
            self._x_axis_group.addButton(btn)
            self._x_axis_buttons[mode] = btn
            axis_lay.addWidget(btn)
        right_lay.addWidget(axis_bar)
        right_lay.addWidget(self._chart, stretch=3)
        right_lay.addWidget(self._ratio_chart, stretch=2)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setChildrenCollapsible(False)
        splitter.setStretchFactor(0, 4)
        splitter.setStretchFactor(1, 6)
        splitter.setSizes([560, 840])
        root.addWidget(splitter, stretch=1)

        self._style_inputs()
        file_menu = self.menuBar().addMenu("File")
        open_action = QAction("Open Output Folder...", self)
        open_action.triggered.connect(self._open_output_folder)
        file_menu.addAction(open_action)
        reload_action = QAction("Reload ROOT", self)
        reload_action.triggered.connect(self._load_root)
        file_menu.addAction(reload_action)

    def _label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Consolas", 10))
        lbl.setStyleSheet(themed(f"color:{THEME.TEXT_DIM};"))
        lbl.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        return lbl

    def _path_edit(self, text: str) -> QLineEdit:
        edit = QLineEdit(text)
        edit.setMinimumWidth(70)
        edit.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        return edit

    def _button(self, text: str, color: str, slot):
        btn = QPushButton(text)
        btn.clicked.connect(slot)
        btn.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
        btn.setStyleSheet(
            f"QPushButton{{background:{THEME.BUTTON};color:{color};"
            f"border:1px solid {THEME.BORDER};border-radius:3px;padding:4px 10px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}"
            f"QPushButton:disabled{{color:{THEME.TEXT_MUTED};}}")
        return btn

    def _hscroll(self, widget: QWidget) -> QScrollArea:
        scroll = QScrollArea()
        widget.setStyleSheet(themed(f"QWidget{{background:{THEME.BG};}}"))
        scroll.setWidget(widget)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        scroll.setMinimumHeight(widget.sizeHint().height() + 4)
        scroll.viewport().setStyleSheet(themed(f"background:{THEME.BG};"))
        scroll.setStyleSheet(themed(f"""
            QScrollArea {{
                background: {THEME.BG};
                border: none;
            }}
            QScrollBar:horizontal {{
                background: {THEME.BG};
                height: 12px;
                margin: 0px;
            }}
            QScrollBar::handle:horizontal {{
                background: {THEME.BORDER};
                min-width: 32px;
            }}
            QScrollBar::add-line:horizontal,
            QScrollBar::sub-line:horizontal {{
                width: 0px;
                height: 0px;
            }}
            QScrollBar::add-page:horizontal,
            QScrollBar::sub-page:horizontal {{
                background: transparent;
            }}
        """))
        return scroll

    def _style_inputs(self):
        ss = themed(
            f"QLineEdit,QSpinBox,QDoubleSpinBox,QComboBox{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};border-radius:3px;padding:2px 5px;"
            "font-family:Consolas;font-size:10pt;}"
            f"QComboBox QAbstractItemView{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"selection-background-color:{THEME.ACCENT_STRONG};}}")
        for cls in (QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox):
            for w in self.findChildren(cls):
                w.setStyleSheet(ss)

    def _browse_base(self):
        d = QFileDialog.getExistingDirectory(self, "Select storage base directory")
        if d:
            self._storage_base.setText(d)
            self._opened_output_folder = None
            self._folder_gain_files = {}
            self._set_folder_run_choices([], preserve=False)
            self._work_dir.clear()
            self._out_root.clear()
            self._update_paths_from_run()

    def _open_output_folder(self):
        start = str(self._opened_output_folder or Path(self._storage_base.text().strip() or STORAGE_BASE))
        d = QFileDialog.getExistingDirectory(self, "Select gain_corr output folder", start)
        if not d:
            return
        self._load_output_folder(Path(d))

    def _browse_ref_gain_file(self):
        start = str(self._ref_gain_file.parent if self._ref_gain_file else DB_DIR / "gain_factor" / "ref_gain")
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select reference gain file",
            start,
            "Reference gain (*.dat);;All files (*)",
        )
        if not path:
            return
        self._ref_gain_file = Path(path)
        self._ref_file_btn.setText(f"Ref: {self._ref_gain_file.name}")
        self._append(f"Using reference gain file: {self._ref_gain_file}", "ok")
        self._ref_file_gain_array = _parse_lms_dat(self._ref_gain_file)
        if self._ref_file_gain_array is None:
            self._append(f"Warning: failed to parse ref gain from {self._ref_gain_file.name}", "warn")
        else:
            self._append(f"Ref gain loaded: {int(np.sum(np.isfinite(self._ref_file_gain_array[:, 0])))} valid W modules", "ok")
        self._ref_gain_cache.clear()
        self._view_context = None
        if hasattr(self, "_runs_data") and self._runs_data:
            self._refresh_views()

    def _on_ref_run_changed(self):
        self._ref_gain_cache.clear()
        self._view_context = None
        if hasattr(self, "_runs_data") and self._runs_data:
            self._refresh_views()

    def _get_ref_run_gain(self) -> Optional[np.ndarray]:
        """Return ref gain array [1156, 3], or None (use embedded gain_W_ref).

        Priority:
          1. Explicitly selected .dat file via Ref File button
          2. prad_XXXXXX_LMS.dat in DB_DIR/gain_factor/ref_gain/ for ref run spinbox
          3. Average gain_W from gain_corr.root in storage
        """
        # Priority 1: directly selected ref gain file
        if self._ref_gain_file is not None:
            return self._ref_file_gain_array

        ref_run_num = self._ref_run.value()
        if ref_run_num < 0:
            return None
        if ref_run_num in self._ref_gain_cache:
            return self._ref_gain_cache[ref_run_num]

        # 1. Try .dat file
        dat_path = DB_DIR / "gain_factor" / "ref_gain" / f"prad_{ref_run_num:06d}_LMS.dat"
        if dat_path.exists():
            result = _parse_lms_dat(dat_path)
            if result is not None:
                self._ref_gain_cache[ref_run_num] = result
                self._append(f"Loaded ref gain from {dat_path.name}", "ok")
                return result
            self._append(f"Failed to parse ref gain file: {dat_path.name}", "warn")

        # 2. Fall back to gain_corr.root
        storage = self._storage_base.text().strip() or STORAGE_BASE
        _, _, path = run_storage_paths(ref_run_num, storage)
        if not path.exists():
            self._append(f"Ref run {ref_run_num:06d}: no .dat or gain_corr.root found", "warn")
            self._ref_gain_cache[ref_run_num] = None
            return None
        try:
            ref_data = load_gain_root(path, ref_run_num)
        except Exception as exc:
            self._append(f"Failed to load ref run {ref_run_num:06d}: {exc}", "error")
            self._ref_gain_cache[ref_run_num] = None
            return None
        arr = np.asarray(ref_data.gain_w, dtype=float)
        valid = np.isfinite(arr) & (arr > 0.0)
        count = np.sum(valid, axis=0)
        total = np.sum(np.where(valid, arr, 0.0), axis=0)
        result = safe_divide(total, count, count > 0)
        self._ref_gain_cache[ref_run_num] = result
        self._append(f"Loaded ref gain from run {ref_run_num:06d}: {path.name}", "ok")
        return result

    def _runs_shown_changed(self):
        if self._opened_output_folder is not None:
            self._load_output_folder(self._opened_output_folder)
        else:
            self._load_root()

    def _set_folder_run_choices(self, runs: List[int], preserve: bool = True):
        old_start = self._display_run_start.currentData() if preserve else None
        old_end = self._display_run_end.currentData() if preserve else None
        self._updating_run_range = True
        try:
            for combo in (self._display_run_start, self._display_run_end):
                combo.clear()
                for run in runs:
                    combo.addItem(f"{run:06d}", run)
                combo.setEnabled(bool(runs))
            if not runs:
                return
            default_start = runs[max(0, len(runs) - self._runs_shown.value())]
            start = old_start if old_start in self._folder_gain_files else default_start
            end = old_end if old_end in self._folder_gain_files else runs[-1]
            if start > end:
                start, end = end, start
            self._display_run_start.setCurrentIndex(runs.index(start))
            self._display_run_end.setCurrentIndex(runs.index(end))
        finally:
            self._updating_run_range = False

    def _display_run_range_changed(self):
        if self._updating_run_range or self._opened_output_folder is None:
            return
        start = self._display_run_start.currentData()
        end = self._display_run_end.currentData()
        if start is None or end is None:
            return
        self._updating_run_range = True
        try:
            if start > end:
                if self.sender() is self._display_run_start:
                    self._display_run_end.setCurrentIndex(self._display_run_start.currentIndex())
                else:
                    self._display_run_start.setCurrentIndex(self._display_run_end.currentIndex())
        finally:
            self._updating_run_range = False
        self._load_selected_folder_runs()

    def _load_selected_folder_runs(self):
        if not self._folder_gain_files:
            return
        start = self._display_run_start.currentData()
        end = self._display_run_end.currentData()
        if start is None or end is None:
            return
        files = [
            (run, path)
            for run, path in sorted(self._folder_gain_files.items())
            if start <= run <= end
        ]
        self._run_edit.setText(" ".join(f"{run:06d}" for run, _ in files))
        self._request_root_load(
            files,
            str(self._opened_output_folder),
            limit=max(1, len(files)),
        )

    def _parse_runs(self) -> List[int]:
        text = self._run_edit.text().strip()
        if not text:
            return []
        runs: List[int] = []
        for tok in re.split(r"[,\s]+", text):
            if not tok:
                continue
            if not tok.isdigit():
                self._append(f"Run number must be numeric: {tok}", "error")
                return []
            runs.append(int(tok))
        return runs

    def _parse_run(self) -> Optional[int]:
        runs = self._parse_runs()
        return runs[0] if runs else None

    def _runs_to_load(self) -> List[int]:
        runs = set(self._parse_runs())
        runs.update(self._known_runs)
        return sorted(runs)

    def _update_paths_from_run(self):
        run = self._parse_run()
        if run is None:
            return
        storage = self._storage_base.text().strip() or STORAGE_BASE
        _, work_dir, out_root = run_storage_paths(run, storage)
        self._work_dir.setText(str(work_dir))
        self._out_root.setText(str(out_root))

    def _start(self):
        seed_runs = self._parse_runs()
        if not seed_runs:
            return
        self._scan_floor_run = min(seed_runs)
        self._opened_output_folder = None
        self._folder_gain_files = {}
        self._set_folder_run_choices([], preserve=False)
        self._update_paths_from_run()
        self._start_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._enqueue_existing_replay_snapshots()
        self._start_next_replay()
        self._scan_once()

    def _stop(self):
        self._timer.stop()
        self._batch_reanalyze_timer.stop()
        self._reanalyze_pending = False
        if self._download_process.state() != QProcess.ProcessState.NotRunning:
            self._download_process.kill()
        if self._replay_process.state() != QProcess.ProcessState.NotRunning:
            self._replay_process.kill()
        if self._reanalyze_process.state() != QProcess.ProcessState.NotRunning:
            self._reanalyze_process.kill()
        self._replay_queue.clear()
        self._queued_snapshots.clear()
        self._active_replay_snapshots.clear()
        self._start_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)
        self._status.setText("Stopped")

    def _scan_once(self):
        seed_runs = self._parse_runs()
        if not seed_runs:
            return
        self._opened_output_folder = None
        if self._download_process.state() != QProcess.ProcessState.NotRunning:
            self._append("Previous download scan still running; skipping this tick.", "warn")
            return
        self._timer.stop()
        self._scan_copied_count = 0
        self._update_paths_from_run()

        host = self._host.text().strip() or REMOTE_HOST
        remote_base = self._remote_base.text().strip() or REMOTE_BASE
        storage_base = self._storage_base.text().strip() or STORAGE_BASE

        seed_run = max(min(seed_runs), self._scan_floor_run or min(seed_runs))
        explicit_runs = " ".join(str(r) for r in seed_runs)
        bash = f"""
set -u
SEED_RUN={seed_run}
EXPLICIT_RUNS={shlex.quote(explicit_runs)}
HOST={shlex.quote(host)}
REMOTE_BASE={shlex.quote(remote_base.rstrip('/'))}
STORAGE_BASE={shlex.quote(storage_base)}
EVIO_SKIP={self._evio_skip.value()}
MAX_DOWNLOAD_EVIO={self._max_download_evio.value()}
QUEUE_CAP_EVIO={self._queue_cap_evio.value()}
echo "Download settings: storage=$STORAGE_BASE"

REMOTE_RUNS=$(ssh "$HOST" "bash -c 'ls \"$REMOTE_BASE\" 2>/dev/null | grep -E \"^prad_[0-9]{{6}}$\" | sort'" || true)
RUNS=""
for r in $EXPLICIT_RUNS; do
    [ "$r" -lt "$SEED_RUN" ] && continue
    case " $RUNS " in *" $r "*) ;; *) RUNS="$RUNS $r" ;; esac
done
while IFS= read -r d; do
    [ -z "$d" ] && continue
    R=${{d#prad_}}
    N=$((10#$R))
    [ "$N" -lt "$SEED_RUN" ] && continue
    case " $RUNS " in *" $N "*) ;; *) RUNS="$RUNS $N" ;; esac
done <<< "$REMOTE_RUNS"

for RUN in $RUNS; do
RUN_TAG=$(printf 'prad_%06d' "$RUN")
REMOTE_DIR="$REMOTE_BASE/$RUN_TAG"
LOCAL_DIR="$STORAGE_BASE/evio/$RUN_TAG"
WORK_DIR="$STORAGE_BASE/lms/$RUN_TAG"
OUT_ROOT="$STORAGE_BASE/gain/$RUN_TAG/${{RUN_TAG}}_gain_corr.root"
QUEUE_ROOT="$STORAGE_BASE/replay_queue/$RUN_TAG"
mkdir -p "$LOCAL_DIR" "$WORK_DIR" "$(dirname "$OUT_ROOT")" "$QUEUE_ROOT"
echo "Scanning $HOST:$REMOTE_DIR"
ALL_FILES=$(ssh "$HOST" "bash -c 'ls \"$REMOTE_DIR\" 2>/dev/null | sort'" || true)
EVIO_FILES=$(printf '%s\\n' "$ALL_FILES" | grep '\\.evio\\.' | sort -V || true)
REMOTE_COUNT=$(printf '%s\\n' "$EVIO_FILES" | grep -c '\\.evio\\.' || true)
LATEST_EVIO=$(printf '%s\\n' "$EVIO_FILES" | tail -n 1)
SAFE_FILES=$(printf '%s\\n' "$EVIO_FILES" | sed '$d')
if [ -n "$LATEST_EVIO" ]; then
    echo "Safety hold: newest remote EVIO will not be downloaded: $LATEST_EVIO"
fi
COPIED=0
ALREADY=0
QUEUED=0
SKIPPED_SAMPLE=0
EVIO_INDEX=0
PENDING_TOTAL=$(find "$STORAGE_BASE/replay_queue" "$STORAGE_BASE/evio" -type f -name '*.evio.*' 2>/dev/null | wc -l | tr -d ' ')
if [ "$PENDING_TOTAL" -ge "$QUEUE_CAP_EVIO" ]; then
    echo "Pending EVIO queue cap reached: $PENDING_TOTAL/$QUEUE_CAP_EVIO; pausing downloads until replay clears files."
fi
while IFS= read -r f; do
    [ -z "$f" ] && continue
    case "$f" in
        *.evio.*) ;;
        *) continue ;;
    esac
    TAKE=1
    if [ "$EVIO_SKIP" -gt 0 ]; then
        MOD=$((EVIO_INDEX % (EVIO_SKIP + 1)))
        if [ "$MOD" -ne 0 ]; then TAKE=0; fi
    fi
    EVIO_INDEX=$((EVIO_INDEX+1))
    if [ "$TAKE" -eq 0 ]; then
        SKIPPED_SAMPLE=$((SKIPPED_SAMPLE+1))
        continue
    fi
    LMS_NAME="${{f/.evio/}}_lms.root"
    if [ -f "$WORK_DIR/$LMS_NAME" ]; then
        echo "Already replayed: $f (skipping download)"
        ALREADY=$((ALREADY+1))
    elif find "$QUEUE_ROOT" -type f -name "$f" -print -quit | grep -q .; then
        echo "Already queued: $f (skipping download)"
        QUEUED=$((QUEUED+1))
    elif [ -f "$LOCAL_DIR/$f" ]; then
        ALREADY=$((ALREADY+1))
    else
        if [ "$PENDING_TOTAL" -ge "$QUEUE_CAP_EVIO" ]; then
            echo "Pending EVIO queue cap reached: $PENDING_TOTAL/$QUEUE_CAP_EVIO; stopping download for now."
            break
        fi
        if [ "$COPIED" -ge "$MAX_DOWNLOAD_EVIO" ]; then
            echo "Reached $MAX_DOWNLOAD_EVIO files limit, stopping download for now."
            break
        fi

        echo "Copying $f"
        scp "$HOST:$REMOTE_DIR/$f" "$LOCAL_DIR/"
        COPIED=$((COPIED+1))
        PENDING_TOTAL=$((PENDING_TOTAL+1))
    fi
done <<< "$SAFE_FILES"

LOCAL_COUNT=$(find "$LOCAL_DIR" -maxdepth 1 -type f -name '*.evio.*' | wc -l | tr -d ' ')
LMS_COUNT=$(find "$WORK_DIR" -maxdepth 1 -type f -name '*_lms.root' | wc -l | tr -d ' ')
PENDING_TOTAL=$(find "$STORAGE_BASE/replay_queue" "$STORAGE_BASE/evio" -type f -name '*.evio.*' 2>/dev/null | wc -l | tr -d ' ')
echo "remote=$REMOTE_COUNT held_latest=$LATEST_EVIO local=$LOCAL_COUNT copied=$COPIED already=$ALREADY queued=$QUEUED pending=$PENDING_TOTAL/$QUEUE_CAP_EVIO sampled_skip=$SKIPPED_SAMPLE lms=$LMS_COUNT"
if [ "$LOCAL_COUNT" -gt 0 ]; then
    SNAP_DIR="$QUEUE_ROOT/$(date +%Y%m%d_%H%M%S)_$$"
    mkdir -p "$SNAP_DIR"
    find "$LOCAL_DIR" -maxdepth 1 -type f -name '*.evio.*' -exec mv -t "$SNAP_DIR" {{}} +
    SNAP_COUNT=$(find "$SNAP_DIR" -maxdepth 1 -type f -name '*.evio.*' | wc -l | tr -d ' ')
    echo "Queued $SNAP_COUNT EVIO file(s) for replay: $SNAP_DIR"
    echo "__ONLINE_GAIN_QUEUE__ run=$RUN snapshot=$SNAP_DIR"
else
    echo "No new local EVIO files queued."
fi
LOCAL_COUNT=$(find "$LOCAL_DIR" -maxdepth 1 -type f -name '*.evio.*' | wc -l | tr -d ' ')
echo "__ONLINE_GAIN_STATUS__ run=$RUN remote=$REMOTE_COUNT local=$LOCAL_COUNT copied=$COPIED lms=$LMS_COUNT out=$OUT_ROOT"
done
exit 0
"""
        self._status.setText("Scanning...")
        self._append(f"$ auto-scan from run {seed_run:06d}", "cmd")
        start_bash_process(self._download_process, bash, getattr(self, "_download_cpu_list", ""))

    def _on_stdout(self):
        proc = self.sender()
        if proc is None:
            return
        self._consume_process_output(proc, is_stderr=False)

    def _on_stderr(self):
        proc = self.sender()
        if proc is None:
            return
        self._consume_process_output(proc, is_stderr=True)

    def _consume_process_output(self, proc: QProcess, is_stderr: bool,
                                final: bool = False):
        key = (id(proc), is_stderr)
        chunk = (
            bytes(proc.readAllStandardError())
            if is_stderr
            else bytes(proc.readAllStandardOutput())
        ).decode(errors="replace")
        text = self._process_output_remainders.get(key, "") + chunk
        if final:
            lines = text.splitlines()
            self._process_output_remainders.pop(key, None)
        else:
            parts = text.split("\n")
            lines = parts[:-1]
            self._process_output_remainders[key] = parts[-1]
        for line in lines:
            line = line.rstrip("\r")
            if not line:
                continue
            if not is_stderr:
                self._record_status_line(line)
                if line.startswith(QUIET_PROCESS_PREFIXES):
                    continue
            self._append(line, "error" if is_stderr else "normal")

    def _record_status_line(self, line: str):
        if "__ONLINE_GAIN_QUEUE__" in line:
            m_run = re.search(r"\brun=(\d+)\b", line)
            m_snap = re.search(r"\bsnapshot=(.+)$", line)
            if m_run and m_snap:
                self._enqueue_replay_snapshot(int(m_run.group(1)), Path(m_snap.group(1).strip()))
                self._start_next_replay()
            return
        if "__ONLINE_GAIN_STATUS__" in line:
            m_run = re.search(r"\brun=(\d+)\b", line)
            if m_run:
                self._known_runs.add(int(m_run.group(1)))
            m_copied = re.search(r"\bcopied=(\d+)\b", line)
            if m_copied:
                self._scan_copied_count += int(m_copied.group(1))
            m_pending = re.search(r"\bpending=(\d+)/\d+\b", line)
            if m_pending:
                self._last_pending_evio_count = int(m_pending.group(1))

    def _on_download_finished(self, code, status):
        self._consume_process_output(self._download_process, is_stderr=False, final=True)
        self._consume_process_output(self._download_process, is_stderr=True, final=True)
        color = "ok" if code == 0 else "error"
        self._append(f"[download scan finished: exit {code}]", color)
        self._start_next_replay()
        if self._stop_btn.isEnabled():
            interval_ms = self._interval.value() * 1000
            if self._pending_evio_count() >= self._queue_cap_evio.value():
                self._append(
                    f"Pending EVIO cap reached; next scan delayed by {self._interval.value()} s.",
                    "warn",
                )
                self._timer.start(interval_ms)
            elif code != 0:
                self._append(
                    f"Download scan failed; retrying in {self._interval.value()} s.",
                    "warn",
                )
                self._timer.start(interval_ms)
            elif self._scan_copied_count > 0:
                self._append(
                    f"Downloaded {self._scan_copied_count} new EVIO file(s); scanning again now.",
                    "ok",
                )
                QTimer.singleShot(0, self._scan_once)
            else:
                self._append(
                    f"No new EVIO files; next scan in {self._interval.value()} s.",
                    "normal",
                )
                self._timer.start(interval_ms)

    def _on_replay_finished(self, code, status):
        self._consume_process_output(self._replay_process, is_stderr=False, final=True)
        self._consume_process_output(self._replay_process, is_stderr=True, final=True)
        color = "ok" if code == 0 else "error"
        self._append(f"[replay finished: exit {code}]", color)
        for snapshot in self._active_replay_snapshots:
            self._queued_snapshots.discard(snapshot)
        self._active_replay_snapshots.clear()
        if code == 0:
            self._last_pending_evio_count = max(
                0, self._last_pending_evio_count - self._active_replay_file_count,
            )
        self._active_replay_file_count = 0
        self._load_root()
        if self._reanalyze_pending:
            self._reanalyze_lms()
        else:
            self._start_next_replay()

    def _enqueue_replay_snapshot(self, run: int, snapshot: Path) -> bool:
        snapshot = snapshot.resolve()
        if snapshot in self._queued_snapshots:
            return False
        if not snapshot.exists() or not snapshot.is_dir():
            return False
        self._advance_scan_floor(run)
        self._queued_snapshots.add(snapshot)
        self._replay_queue.append((run, snapshot))
        return True

    def _advance_scan_floor(self, run: int):
        if self._scan_floor_run is not None and run <= self._scan_floor_run:
            return
        previous = self._scan_floor_run
        self._scan_floor_run = run
        if previous is not None:
            self._append(
                f"Scan advanced to run {run:06d}; runs before it will no longer be scanned.",
                "ok",
            )

    def _enqueue_existing_replay_snapshots(self):
        storage_base = Path(self._storage_base.text().strip() or STORAGE_BASE).expanduser()
        queue_root = storage_base / "replay_queue"
        if not queue_root.exists():
            return
        added = 0
        recovered_files = 0
        for run_dir in sorted(queue_root.glob("prad_[0-9][0-9][0-9][0-9][0-9][0-9]")):
            if not run_dir.is_dir():
                continue
            try:
                run = int(run_dir.name.split("_", 1)[1])
            except (IndexError, ValueError):
                continue
            for snapshot in sorted(p for p in run_dir.iterdir() if p.is_dir()):
                if self._enqueue_replay_snapshot(run, snapshot):
                    added += 1
                    recovered_files += sum(
                        1 for path in snapshot.glob("*.evio.*") if path.is_file()
                    )
        self._last_pending_evio_count = max(
            self._last_pending_evio_count, recovered_files,
        )
        if added:
            self._append(f"Recovered {added} replay snapshot(s) from disk queue.", "warn")

    def _pending_evio_count(self) -> int:
        return self._last_pending_evio_count

    def _schedule_batch_reanalysis(self, _value=None):
        self._batch_reanalyze_timer.start()

    def _reanalyze_lms(self):
        self._batch_reanalyze_timer.stop()
        if (
            self._replay_process.state() != QProcess.ProcessState.NotRunning
            or self._reanalyze_process.state() != QProcess.ProcessState.NotRunning
        ):
            self._reanalyze_pending = True
            self._status.setText("Re-analysis queued")
            return

        self._reanalyze_pending = False
        runs = sorted(self._visible_runs_data())
        if not runs:
            runs = sorted(set(self._parse_runs()))
        storage_base = self._storage_base.text().strip() or STORAGE_BASE
        jobs: List[Tuple[int, Path, Path]] = []
        for run in runs:
            _, work_dir, out_root = run_storage_paths(run, storage_base)
            if work_dir.is_dir() and any(work_dir.glob("*_lms.root")):
                jobs.append((run, work_dir, out_root))
        if not jobs:
            self._status.setText("No LMS ROOT files to re-analyze")
            self._append("No existing *_lms.root files found for the selected run(s).", "warn")
            self._start_next_replay()
            return

        tool = find_update_tool()
        ref_arg = ""
        if self._ref_gain_file is not None:
            ref_arg = f" -R {shlex.quote(str(self._ref_gain_file))}"
        elif self._ref_run.value() >= 0:
            ref_arg = f" -r {self._ref_run.value()}"

        commands = []
        for run, work_dir, out_root in jobs:
            commands.append(
                f'echo "Re-analyzing run {run:06d}: batch={self._batch.value()}"\n'
                f'mkdir -p {shlex.quote(str(out_root.parent))}\n'
                f'nice -n 10 {shlex.quote(tool)} -a '
                f'-w {shlex.quote(str(work_dir))} '
                f'-o {shlex.quote(str(out_root))} '
                f'-b {self._batch.value()}{ref_arg}'
            )
        bash = "set -euo pipefail\n" + "\n".join(commands)

        self._reanalyze_btn.setEnabled(False)
        self._status.setText("Re-analyzing LMS ROOT files...")
        self._append(
            f"$ {tool} -a  # re-analyze {len(jobs)} run(s), batch size {self._batch.value()}",
            "cmd",
        )
        replay_cpus = getattr(self, "_replay_worker_cpus", [])
        cpu_list = _format_cpu_list(replay_cpus[:1]) if replay_cpus else ""
        start_bash_process(self._reanalyze_process, bash, cpu_list)

    def _on_reanalyze_finished(self, code, status):
        self._consume_process_output(self._reanalyze_process, is_stderr=False, final=True)
        self._consume_process_output(self._reanalyze_process, is_stderr=True, final=True)
        self._reanalyze_btn.setEnabled(True)
        color = "ok" if code == 0 else "error"
        self._append(f"[re-analysis finished: exit {code}]", color)
        self._status.setText("Re-analysis finished" if code == 0 else "Re-analysis failed")
        self._load_root()
        if self._reanalyze_pending:
            self._reanalyze_lms()
        else:
            self._start_next_replay()

    def _start_next_replay(self):
        if (
            self._replay_process.state() != QProcess.ProcessState.NotRunning
            or self._reanalyze_process.state() != QProcess.ProcessState.NotRunning
        ):
            return
        while self._replay_queue:
            run, snapshot = self._replay_queue.pop(0)
            if snapshot.exists():
                break
            self._queued_snapshots.discard(snapshot)
            self._append(f"Replay snapshot missing, skipped: {snapshot}", "warn")
        else:
            return
        self._advance_scan_floor(run)
        snapshots = [snapshot]
        kept_queue: List[Tuple[int, Path]] = []
        for qrun, qsnap in self._replay_queue:
            if qrun == run and qsnap.exists():
                snapshots.append(qsnap)
            else:
                kept_queue.append((qrun, qsnap))
        self._replay_queue = kept_queue
        self._active_replay_snapshots = snapshots
        self._active_replay_file_count = sum(
            1
            for snapshot_dir in snapshots
            for path in snapshot_dir.glob("*.evio.*")
            if path.is_file()
        )

        storage_base = self._storage_base.text().strip() or STORAGE_BASE
        _, work_dir, out_root = run_storage_paths(run, storage_base)
        tool = find_update_tool()
        requested_threads = self._threads.value()
        replay_cpus = getattr(self, "_replay_worker_cpus", [])
        replay_cpu_capacity = len(replay_cpus) if replay_cpus else max(1, (os.cpu_count() or 2) - 2)
        replay_threads = min(
            requested_threads,
            max(1, replay_cpu_capacity),
            max(1, self._active_replay_file_count),
        )
        replay_cpu_list = ""
        if replay_cpus and self._viewer_cpu is not None:
            replay_cpu_list = _format_cpu_list(replay_cpus[:replay_threads])
        if replay_threads < requested_threads:
            self._append(
                f"Replay uses {replay_threads}/{requested_threads} configured threads "
                f"for {self._active_replay_file_count} queued file(s); "
                "unused workers are not created.",
                "normal",
            )

        ref_arg = ""
        if self._ref_gain_file is not None:
            ref_arg = f" -R {shlex.quote(str(self._ref_gain_file))}"
        elif self._ref_run.value() >= 0:
            ref_arg = f" -r {self._ref_run.value()}"

        snapshot_args = " ".join(shlex.quote(str(p)) for p in snapshots)
        bash = f"""
set -u
RUN={run}
RUN_TAG=$(printf 'prad_%06d' "$RUN")
SNAP_DIRS=({snapshot_args})
STORAGE_BASE={shlex.quote(storage_base)}
WORK_DIR={shlex.quote(str(work_dir))}
OUT_ROOT={shlex.quote(str(out_root))}
TOOL={shlex.quote(tool)}
BATCH={self._batch.value()}
THREADS={replay_threads}
REPLAY_CPU_LIST={shlex.quote(replay_cpu_list)}
mkdir -p "$WORK_DIR" "$(dirname "$OUT_ROOT")"
SNAP_COUNT=0
for d in "${{SNAP_DIRS[@]}}"; do
    C=$(find "$d" -maxdepth 1 -type f -name '*.evio.*' 2>/dev/null | wc -l | tr -d ' ')
    SNAP_COUNT=$((SNAP_COUNT + C))
done
echo "Replay settings: run=$RUN_TAG batch=$BATCH replay_threads=$THREADS replay_cpus=${{REPLAY_CPU_LIST:-all}} snapshots=${{#SNAP_DIRS[@]}} snapshot_files=$SNAP_COUNT"
if [ "$SNAP_COUNT" -eq 0 ]; then
    echo "Empty replay snapshots, removing empty directories."
    for d in "${{SNAP_DIRS[@]}}"; do rmdir "$d" 2>/dev/null || true; done
    exit 0
fi
RUN_CMD=()
if [ -n "$REPLAY_CPU_LIST" ] && command -v taskset >/dev/null 2>&1; then
    RUN_CMD+=(taskset -c "$REPLAY_CPU_LIST")
fi
RUN_CMD+=(nice -n 10)
if command -v ionice >/dev/null 2>&1; then
    RUN_CMD+=(ionice -c2 -n7)
fi
RUN_CMD+=("$TOOL")
if "${{RUN_CMD[@]}}" -w "$WORK_DIR" -o "$OUT_ROOT" -b "$BATCH" -j "$THREADS"{ref_arg} "${{SNAP_DIRS[@]}}"; then
    echo "Replay/update succeeded for $RUN_TAG; deleting queued EVIO snapshots."
    for d in "${{SNAP_DIRS[@]}}"; do
        case "$d" in
            "$STORAGE_BASE"/replay_queue/prad_[0-9][0-9][0-9][0-9][0-9][0-9]/*)
                rm -rf "$d"
                ;;
            *)
                echo "Refusing to delete unexpected snapshot path: $d" >&2
                ;;
        esac
    done
else
    echo "ERROR: gain update failed for $RUN_TAG; keeping snapshots" >&2
    exit 1
fi
LMS_COUNT=$(find "$WORK_DIR" -maxdepth 1 -type f -name '*_lms.root' | wc -l | tr -d ' ')
echo "__ONLINE_GAIN_STATUS__ run=$RUN remote=0 local=0 copied=0 lms=$LMS_COUNT out=$OUT_ROOT"
"""
        self._status.setText("Replaying...")
        self._append(f"$ replay {len(snapshots)} queued snapshot(s) for run {run:06d}", "cmd")
        start_bash_process(self._replay_process, bash, replay_cpu_list)

    def _append(self, text: str, kind: str = "normal"):
        colors = {
            "normal": THEME.TEXT,
            "cmd": THEME.TEXT_DIM,
            "warn": THEME.WARN,
            "error": THEME.DANGER,
            "ok": THEME.SUCCESS,
        }
        safe = (
            text.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
        )
        self._log_buffer.append((safe, colors.get(kind, THEME.TEXT)))
        if len(self._log_buffer) > LOG_BUFFER_LIMIT:
            del self._log_buffer[:-LOG_BUFFER_LIMIT]
        if not self._log_flush_timer.isActive():
            self._log_flush_timer.start(LOG_FLUSH_MS)

    def _flush_log(self):
        if not self._log_buffer:
            return
        lines = self._log_buffer
        self._log_buffer = []
        html = "".join(
            f"<span style='color:{color}'>{text}</span><br>"
            for text, color in lines
        )
        self._log.moveCursor(self._log.textCursor().MoveOperation.End)
        self._log.insertHtml(html)
        self._log.moveCursor(self._log.textCursor().MoveOperation.End)

    def _maintenance(self):
        self._flush_log()
        keep = (
            len(self._runs_data)
            if self._opened_output_folder is not None
            else max(DEFAULT_RUNS_SHOWN, self._runs_shown.value())
        )
        retained = sorted(self._runs_data)[-keep:]
        retained_set = set(retained)
        if len(retained_set) < len(self._runs_data):
            self._runs_data = {
                run: self._runs_data[run]
                for run in retained
            }
            self._root_signatures = {
                run: sig
                for run, sig in self._root_signatures.items()
                if run in retained_set
            }
            self._data = self._runs_data[max(self._runs_data)] if self._runs_data else None
        gc.collect(0)

    def _request_root_load(self, files: List[Tuple[int, Path]], label: str,
                           missing: int = 0, limit: Optional[int] = None):
        by_run = {run: path for run, path in files}
        limit = max(1, self._runs_shown.value() if limit is None else limit)
        selected = sorted(by_run.items())[-limit:]
        if not selected:
            self._status.setText("No gain_corr.root yet")
            return
        if self._root_load_thread is not None:
            self._pending_root_load = (selected, label, missing)
            return
        self._begin_root_load(selected, label, missing)

    def _begin_root_load(self, files: List[Tuple[int, Path]], label: str,
                         missing: int):
        signatures: Dict[int, Tuple[str, int, int]] = {}
        cached: Dict[int, GainData] = {}
        fallback: Dict[int, GainData] = {}
        changed: List[Tuple[int, Path]] = []
        for run, path in files:
            try:
                stat = path.stat()
            except OSError:
                missing += 1
                continue
            signature = (str(path.resolve()), stat.st_mtime_ns, stat.st_size)
            signatures[run] = signature
            old = self._runs_data.get(run)
            if old is not None and old.path == path:
                fallback[run] = old
            if old is not None and self._root_signatures.get(run) == signature:
                cached[run] = old
            else:
                changed.append((run, path))

        self._root_load_context = (cached, fallback, signatures, label, missing)
        if not changed:
            self._apply_root_load({}, [])
            return

        self._status.setText(f"Loading {len(changed)} updated ROOT file(s)...")
        thread = QThread(self)
        background_cpus = getattr(self, "_replay_cpus", [])
        worker = GainRootLoader(
            changed,
            background_cpus[-1:] if background_cpus else [],
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._apply_root_load)
        worker.finished.connect(thread.quit, Qt.ConnectionType.DirectConnection)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(self._root_thread_finished)
        thread.finished.connect(thread.deleteLater)
        self._root_load_thread = thread
        self._root_load_worker = worker
        thread.start()

    def _apply_root_load(self, loaded, errors):
        if self._root_load_context is None:
            return
        cached, fallback, signatures, label, missing = self._root_load_context
        merged = dict(cached)
        merged.update(loaded)
        failed_runs = {run for run, _, _ in errors}
        for run in failed_runs:
            if run in fallback:
                merged[run] = fallback[run]
        for run, path, message in errors:
            self._append(f"{path}: {message}", "error")
        if not merged:
            self._status.setText("No readable gain_corr ROOT files")
            self._root_load_context = None
            return

        self._runs_data = dict(sorted(merged.items()))
        self._root_signatures = {
            run: signatures[run]
            for run in self._runs_data
            if run in signatures and run not in failed_runs
        }
        for run in failed_runs:
            if run in self._runs_data and run in self._root_signatures:
                self._root_signatures.pop(run, None)
        self._data = self._runs_data[max(self._runs_data)]
        nbatches = sum(data.nbatches for data in self._runs_data.values())
        self._valid_runs_count.setText(f"{len(self._runs_data)} loaded")
        miss = f" | missing {missing} run(s)" if missing else ""
        self._status.setText(
            f"{len(self._runs_data)} run(s), {nbatches} batches loaded from {label}{miss}"
        )
        self._view_context = None
        self._refresh_views()
        self._check_gain_drop_warning()
        self._root_load_context = None

    def _root_thread_finished(self):
        self._root_load_thread = None
        self._root_load_worker = None
        pending = self._pending_root_load
        self._pending_root_load = None
        if pending is not None:
            QTimer.singleShot(0, lambda request=pending: self._begin_root_load(*request))

    def _load_root(self):
        if self._opened_output_folder is not None:
            self._load_output_folder(self._opened_output_folder)
            return
        storage = self._storage_base.text().strip() or STORAGE_BASE
        files_by_run: Dict[int, Path] = {}

        # Discover actual output files first. Run-number gaps must not consume
        # the "runs shown" limit.
        gain_dir = Path(storage).expanduser() / "gain"
        for path in gain_dir.glob("prad_*/prad_*_gain_corr.root"):
            if not path.is_file():
                continue
            match = GAIN_CORR_FILE_RE.match(path.name)
            if match:
                files_by_run[int(match.group(1))] = path

        # Also cover explicit/just-discovered runs when a nonstandard layout
        # or a newly-created file has not appeared in the directory scan yet.
        missing: List[int] = []
        for run in self._runs_to_load():
            if run in files_by_run:
                continue
            _, _, path = run_storage_paths(run, storage)
            if path.is_file():
                files_by_run[run] = path
            else:
                missing.append(run)

        limit = max(1, self._runs_shown.value())
        files = sorted(files_by_run.items())[-limit:]
        if not files:
            self._valid_runs_count.setText("0 loaded")
        self._request_root_load(
            files, "monitor storage", len(missing),
        )

    def _load_output_folder(self, folder: Path):
        if not folder.exists() or not folder.is_dir():
            self._status.setText("Output folder not found")
            self._append(f"Output folder not found: {folder}", "error")
            return
        files_by_run: Dict[int, Path] = {}
        for path in sorted(folder.rglob("*gain_corr*root")):
            if not path.is_file():
                continue
            m = GAIN_CORR_FILE_RE.match(path.name)
            if not m:
                continue
            files_by_run[int(m.group(1))] = path
        if not files_by_run:
            self._folder_gain_files = {}
            self._set_folder_run_choices([], preserve=False)
            self._valid_runs_count.setText("0 loaded")
            self._status.setText("No gain_corr ROOT files")
            self._append(f"No prad_XXXXXX_gain_corr.root files under {folder}", "warn")
            return

        previous_folder = self._opened_output_folder
        self._opened_output_folder = folder
        self._folder_gain_files = files_by_run
        runs = sorted(files_by_run)
        self._known_runs.update(runs)
        self._set_folder_run_choices(runs, preserve=previous_folder == folder)
        self._append(f"Loading gain_corr files from {folder}", "ok")
        self._load_selected_folder_runs()

    def _on_module_clicked(self, name: str):
        if name:
            self._selected_module = name
            self._map.set_selected(name)
            # The selected module only changes the time-series chart. Map
            # values, gain arrays and run averages are unchanged, so avoid a
            # full view refresh on this latency-sensitive click path.
            self._refresh_chart()

    def _module_index(self, name: str) -> Optional[int]:
        m = re.match(r"^W(\d+)$", name)
        if not m:
            return None
        idx = int(m.group(1)) - 1
        return idx if 0 <= idx < 1156 else None

    def _ref_ratio_bad_rel(self) -> float:
        return float(self._ref_ratio_tol.value()) / 100.0

    def _quantity_array(self, data: GainData):
        if self._quantity.currentIndex() == 0:
            ref_gain = self._get_ref_run_gain()
            if ref_gain is not None:
                gain_w = np.asarray(data.gain_w, dtype=float)
                valid = np.isfinite(gain_w) & np.isfinite(ref_gain) & (ref_gain > 0.0)
                return safe_divide(gain_w, ref_gain, valid)
            valid = (
                np.isfinite(data.gain_w)
                & np.isfinite(data.gain_w_ref)
                & (data.gain_w_ref > 0.0)
            )
            return safe_divide(data.gain_w, data.gain_w_ref, valid)
        return data.gain_w

    def _quantity_key(self) -> str:
        idx = self._quantity.currentIndex()
        if idx < 0 or idx >= len(QUANTITIES):
            return QUANTITIES[0][0]
        return QUANTITIES[idx][0]

    def _base_gain_array(self, data: GainData):
        ref_gain = self._get_ref_run_gain()
        if ref_gain is not None:
            gain_w = np.asarray(data.gain_w, dtype=float)
            valid = np.isfinite(gain_w) & np.isfinite(ref_gain) & (ref_gain > 0.0)
            return safe_divide(gain_w, ref_gain, valid)
        valid = (
            np.isfinite(data.gain_w)
            & np.isfinite(data.gain_w_ref)
            & (data.gain_w_ref > 0.0)
        )
        return safe_divide(data.gain_w, data.gain_w_ref, valid)

    def _is_change_quantity(self) -> bool:
        return self._quantity_key() in {
            "gain_change",
            "gain_long_change",
            "gain_run_change",
            "gain_all_run_change",
        }

    def _change_map_range(self) -> float:
        if not hasattr(self, "_change_range"):
            return CHANGE_DEFAULT_RANGE
        return max(0.001, float(self._change_range.value()) / 100.0)

    def _visible_runs_data(self) -> Dict[int, GainData]:
        if not self._runs_data:
            return {}
        limit = (
            len(self._runs_data)
            if self._opened_output_folder is not None
            else self._runs_shown.value() if hasattr(self, "_runs_shown") else DEFAULT_RUNS_SHOWN
        )
        runs = [run for run in sorted(self._runs_data) if self._runs_data[run].nbatches > 0][-limit:]
        return {run: self._runs_data[run] for run in runs}

    def _visible_quantity_arrays(self, runs_data: Dict[int, GainData]):
        arrays = {
            run: np.asarray(self._quantity_array(data), dtype=float)
            for run, data in runs_data.items()
        }
        if not arrays or not self._norm_first5.isChecked():
            return arrays

        first_batches = []
        remaining = 5
        for run in sorted(runs_data):
            arr = arrays[run]
            if arr.size == 0:
                continue
            take = min(remaining, arr.shape[0])
            if take > 0:
                first_batches.append(arr[:take])
                remaining -= take
            if remaining <= 0:
                break
        if not first_batches:
            return arrays

        base_samples = np.concatenate(first_batches, axis=0)
        valid = np.isfinite(base_samples) & (base_samples > 0.0)
        count = np.sum(valid, axis=0)
        total = np.sum(np.where(valid, base_samples, 0.0), axis=0)
        baseline = safe_divide(total, count, count > 0)
        baseline = np.where(np.isfinite(baseline) & (baseline > 0.0), baseline, np.nan)
        for run in arrays:
            arrays[run] = safe_divide(arrays[run], baseline)
        return arrays

    def _visible_base_gain_arrays(self, runs_data: Dict[int, GainData]):
        return {
            run: np.asarray(self._base_gain_array(data), dtype=float)
            for run, data in runs_data.items()
        }

    def _flatten_batches(self, runs_data: Dict[int, GainData], arrays_by_run, masks_by_run):
        batches = []
        for run, data in sorted(runs_data.items()):
            arr = arrays_by_run[run]
            masks = masks_by_run.get(run, [])
            for b in range(data.nbatches):
                mask = masks[b] if b < len(masks) else [False, False, False]
                batches.append((run, b, arr[b], mask))
        return batches

    def _change_array(self, runs_data: Dict[int, GainData], masks_by_run):
        arrays_by_run = self._visible_base_gain_arrays(runs_data)
        batches = self._flatten_batches(runs_data, arrays_by_run, masks_by_run)
        if len(batches) < 2:
            return None

        quantity = self._quantity_key()
        if quantity == "gain_change":
            return self._short_change_array_from_batches(batches)
        if quantity in {"gain_run_change", "gain_all_run_change"}:
            return self._latest_run_average_change(
                runs_data,
                arrays_by_run,
                masks_by_run,
                use_all_previous=quantity == "gain_all_run_change",
            )

        n = min(5, len(batches))
        early = batches[:n]
        late = batches[-n:]
        early_mean = masked_mean([b[2] for b in early], [b[3] for b in early])
        late_mean = masked_mean([b[2] for b in late], [b[3] for b in late])
        if early_mean is None or late_mean is None:
            return None
        valid = np.isfinite(early_mean) & np.isfinite(late_mean) & (early_mean > 0.0)
        return safe_divide(late_mean - early_mean, early_mean, valid)

    def _run_mean_arrays(self, runs_data: Dict[int, GainData], arrays_by_run, masks_by_run):
        means = {}
        for run in sorted(runs_data):
            arr = arrays_by_run[run]
            masks = masks_by_run.get(run, [])
            means[run] = masked_mean(
                [arr[b] for b in range(arr.shape[0])],
                [masks[b] if b < len(masks) else [False, False, False]
                 for b in range(arr.shape[0])],
            )
        return means

    def _previous_run_baseline(self, runs: List[int], run_means, run_index: int,
                               use_all_previous: bool):
        previous = [
            run_means[run]
            for run in runs[:run_index]
            if run_means.get(run) is not None
        ]
        if not previous:
            return None
        if not use_all_previous:
            return previous[-1]
        arr = np.asarray(previous, dtype=float)
        valid = np.isfinite(arr)
        count = np.sum(valid, axis=0)
        total = np.sum(np.where(valid, arr, 0.0), axis=0)
        return safe_divide(total, count, count > 0)

    def _latest_run_average_change(self, runs_data: Dict[int, GainData], arrays_by_run,
                                   masks_by_run, use_all_previous: bool):
        runs = sorted(runs_data)
        if len(runs) < 2:
            return None
        run_means = self._run_mean_arrays(runs_data, arrays_by_run, masks_by_run)
        current = run_means.get(runs[-1])
        baseline = self._previous_run_baseline(
            runs, run_means, len(runs) - 1, use_all_previous,
        )
        if current is None or baseline is None:
            return None
        valid = np.isfinite(current) & np.isfinite(baseline) & (baseline > 0.0)
        return safe_divide(current - baseline, baseline, valid)

    def _short_change_array_from_batches(self, batches):
        if len(batches) < 2:
            return None
        current = batches[-1][2]
        current_mask = batches[-1][3]
        prev = batches[max(0, len(batches) - 6):-1]
        baseline = masked_mean([b[2] for b in prev], [b[3] for b in prev])
        if baseline is None:
            return None
        valid = (
            np.asarray(current_mask, dtype=bool)[None, :]
            & np.isfinite(current)
            & np.isfinite(baseline)
            & (baseline > 0.0)
        )
        return safe_divide(current - baseline, baseline, valid)

    def _short_change_array(self, runs_data: Dict[int, GainData], masks_by_run):
        arrays_by_run = self._visible_base_gain_arrays(runs_data)
        batches = self._flatten_batches(runs_data, arrays_by_run, masks_by_run)
        return self._short_change_array_from_batches(batches)

    def _chart_arrays(self, runs_data: Dict[int, GainData], masks_by_run):
        if not self._is_change_quantity():
            return self._visible_quantity_arrays(runs_data)

        base_by_run = self._visible_base_gain_arrays(runs_data)
        batches = self._flatten_batches(runs_data, base_by_run, masks_by_run)
        results = {
            run: np.full_like(base_by_run[run], np.nan, dtype=float)
            for run in runs_data
        }
        if len(batches) < 2:
            return results

        quantity = self._quantity_key()
        if quantity in {"gain_run_change", "gain_all_run_change"}:
            runs = sorted(runs_data)
            run_means = self._run_mean_arrays(runs_data, base_by_run, masks_by_run)
            use_all_previous = quantity == "gain_all_run_change"
            for run_index, run in enumerate(runs):
                baseline = self._previous_run_baseline(
                    runs, run_means, run_index, use_all_previous,
                )
                if baseline is None:
                    continue
                arr = base_by_run[run]
                masks = np.asarray(masks_by_run.get(run, []), dtype=bool)
                if masks.shape != (arr.shape[0], 3):
                    continue
                valid = (
                    masks[:, None, :] & np.isfinite(arr)
                    & np.isfinite(baseline)[None, :, :]
                    & (baseline[None, :, :] > 0.0)
                )
                results[run] = safe_divide(
                    arr - baseline[None, :, :], baseline[None, :, :], valid,
                )
            return results

        early_n = min(5, len(batches))
        early_baseline = masked_mean(
            [b[2] for b in batches[:early_n]],
            [b[3] for b in batches[:early_n]],
        )
        ordered_runs = sorted(runs_data)
        values = np.concatenate([base_by_run[run] for run in ordered_runs], axis=0)
        ref_masks = np.concatenate(
            [np.asarray(masks_by_run.get(run, []), dtype=bool) for run in ordered_runs],
            axis=0,
        )
        if quantity == "gain_change":
            sample_valid = ref_masks[:, None, :] & np.isfinite(values) & (values > 0.0)
            sample_values = np.where(sample_valid, values, 0.0)
            sums = np.concatenate(
                [np.zeros((1,) + values.shape[1:], dtype=float), np.cumsum(sample_values, axis=0)],
                axis=0,
            )
            counts = np.concatenate(
                [np.zeros((1,) + values.shape[1:], dtype=np.int64), np.cumsum(sample_valid, axis=0)],
                axis=0,
            )
            ends = np.arange(values.shape[0])
            starts = np.maximum(0, ends - 5)
            window_count = counts[ends] - counts[starts]
            baseline = safe_divide(
                sums[ends] - sums[starts], window_count, window_count > 0,
            )
        else:
            baseline = np.broadcast_to(early_baseline, values.shape) if early_baseline is not None else None
        if baseline is not None:
            valid = (
                ref_masks[:, None, :] & np.isfinite(values)
                & np.isfinite(baseline) & (baseline > 0.0)
            )
            changes = safe_divide(values - baseline, baseline, valid)
            offset = 0
            for run in ordered_runs:
                size = base_by_run[run].shape[0]
                results[run] = changes[offset:offset + size]
                offset += size
        return results

    def _ref_ratio_centers_and_masks(self):
        runs_data = self._visible_runs_data()
        ratios = [np.asarray(data.ref_ratio, dtype=float) for data in runs_data.values()]
        stacked = np.concatenate(ratios, axis=0) if ratios else np.empty((0, 3))
        valid = np.isfinite(stacked) & (stacked > 0.0)
        counts = np.sum(valid, axis=0)
        centers_array = safe_divide(
            np.sum(np.where(valid, stacked, 0.0), axis=0), counts, counts > 0,
        )
        centers = centers_array.tolist()
        bad_rel = self._ref_ratio_bad_rel()
        masks = {}
        for run, data in runs_data.items():
            values = np.asarray(data.ref_ratio, dtype=float)
            center_valid = np.isfinite(centers_array) & (centers_array > 0.0)
            masks[run] = (
                np.isfinite(values) & (values > 0.0) & center_valid[None, :]
                & (np.abs(values - centers_array[None, :])
                   <= bad_rel * centers_array[None, :])
            )
        return centers, masks

    def _get_view_context(self):
        if self._view_context is not None:
            return self._view_context
        runs_data = self._visible_runs_data()
        if not runs_data:
            self._view_context = ({}, [], {}, {})
            return self._view_context
        centers, masks = self._ref_ratio_centers_and_masks()
        arrays = self._chart_arrays(runs_data, masks)
        self._view_context = (runs_data, centers, masks, arrays)
        return self._view_context

    def _latest_map_array(self, runs_data: Dict[int, GainData], arrays_by_run):
        latest_run = max(runs_data)
        quantity = self._quantity_key()
        if not self._is_change_quantity() or quantity == "gain_change":
            arr = arrays_by_run[latest_run]
            return arr[-1] if arr.shape[0] else None

        samples = []
        if quantity == "gain_long_change":
            remaining = 5
            for run in reversed(sorted(runs_data)):
                arr = arrays_by_run[run]
                take = min(remaining, arr.shape[0])
                if take:
                    samples.append(arr[-take:])
                    remaining -= take
                if remaining <= 0:
                    break
        else:
            samples.append(arrays_by_run[latest_run])
        if not samples:
            return None
        stacked = np.concatenate(list(reversed(samples)), axis=0)
        valid = np.isfinite(stacked)
        count = np.sum(valid, axis=0)
        total = np.sum(np.where(valid, stacked, 0.0), axis=0)
        return np.divide(
            total,
            count,
            out=np.full_like(total, np.nan, dtype=float),
            where=count > 0,
        )

    def _selected_values_for_batch(self) -> Dict[str, float]:
        runs_data, _, _, arrays_by_run = self._get_view_context()
        if not runs_data:
            return {}
        latest = self._latest_map_array(runs_data, arrays_by_run)
        if latest is None:
            return {}
        ref_idx = self._ref.currentIndex()
        if ref_idx < 3:
            selected = np.asarray(latest[:, ref_idx], dtype=float)
        else:
            latest = np.asarray(latest, dtype=float)
            valid = np.isfinite(latest)
            if not self._is_change_quantity():
                valid &= (latest > 0.0) & (latest != 1.0)
            count = np.sum(valid, axis=1)
            selected = safe_divide(
                np.sum(np.where(valid, latest, 0.0), axis=1), count, count > 0,
            )
        return {f"W{i + 1}": float(value) for i, value in enumerate(selected)}

    def _run_average_value(self, arr, masks) -> float:
        ref_idx = self._ref.currentIndex()
        values = np.asarray(arr, dtype=float)
        batch_masks = np.asarray(masks, dtype=bool)
        if values.size == 0 or batch_masks.shape != (values.shape[0], 3):
            return math.nan
        if ref_idx < 3:
            sample = values[:, :, ref_idx]
            valid = np.isfinite(sample) & batch_masks[:, ref_idx, None]
            if not self._is_change_quantity():
                valid &= sample > 0.0
            return float(np.mean(sample[valid])) if np.any(valid) else math.nan

        valid = np.isfinite(values) & batch_masks[:, None, :]
        if not self._is_change_quantity():
            valid &= (values > 0.0) & (values != 1.0)
        count = np.sum(valid, axis=2)
        per_module = safe_divide(
            np.sum(np.where(valid, values, 0.0), axis=2), count, count > 0,
        )
        finite = np.isfinite(per_module)
        return float(np.mean(per_module[finite])) if np.any(finite) else math.nan

    def _fmt_run_average(self, value: float) -> str:
        if not math.isfinite(value):
            return "n/a"
        if self._is_change_quantity():
            return f"{value * 100:+.2f}%"
        return f"{value:.4f}"

    def _refresh_run_average_monitor(self):
        if not hasattr(self, "_run_avg_label"):
            return
        runs_data, _, masks_by_run, arrays_by_run = self._get_view_context()
        if not runs_data:
            self._run_avg_label.setText("Run avg: n/a")
            return
        runs = sorted(runs_data)
        run_avgs = {
            run: self._run_average_value(arrays_by_run[run], masks_by_run.get(run, []))
            for run in runs
        }
        all_vals = [v for v in run_avgs.values() if math.isfinite(v)]
        all_avg = finite_avg(all_vals)
        current = runs[-1]
        previous = runs[-2] if len(runs) >= 2 else None
        prev_text = f"{previous:06d}: {self._fmt_run_average(run_avgs[previous])}" if previous else "n/a"
        cur_text = f"{current:06d}: {self._fmt_run_average(run_avgs[current])}"
        qlabel = QUANTITIES[self._quantity.currentIndex()][1]
        ref_label = REF_CHOICES[self._ref.currentIndex()]
        self._run_avg_label.setText(
            f"Run avg {qlabel} [{ref_label}]   Prev {prev_text}   Current {cur_text}   All {self._fmt_run_average(all_avg)}"
        )

    def _refresh_gain_summary(self):
        if not hasattr(self, "_gain_summary"):
            return
        runs_data, _, masks_by_run, _ = self._get_view_context()
        labels = ("under absorber", "open layer 1", "open layer 2")
        if not runs_data:
            self._gain_summary.setPlainText("\n".join(f"{label}: n/a" for label in labels))
            return

        remaining = 5
        samples = []
        sample_masks = []
        ref_gain = self._get_ref_run_gain()
        for run in reversed(sorted(runs_data)):
            data = runs_data[run]
            take = min(remaining, data.nbatches)
            if take <= 0:
                continue
            gain = np.asarray(data.gain_w[-take:], dtype=float)
            if ref_gain is not None:
                valid = np.isfinite(gain) & np.isfinite(ref_gain)[None, :, :] & (ref_gain[None, :, :] > 0.0)
                ratio = safe_divide(gain, ref_gain[None, :, :], valid)
            else:
                embedded_ref = np.asarray(data.gain_w_ref[-take:], dtype=float)
                valid = np.isfinite(gain) & np.isfinite(embedded_ref) & (embedded_ref > 0.0)
                ratio = safe_divide(gain, embedded_ref, valid)
            masks = np.asarray(masks_by_run.get(run, []), dtype=bool)
            if masks.shape != (data.nbatches, 3):
                masks = np.zeros((data.nbatches, 3), dtype=bool)
            samples.append(ratio)
            sample_masks.append(masks[-take:])
            remaining -= take
            if remaining == 0:
                break

        if not samples:
            self._gain_summary.setPlainText("\n".join(f"{label}: n/a" for label in labels))
            return
        values = np.concatenate(list(reversed(samples)), axis=0)
        masks = np.concatenate(list(reversed(sample_masks)), axis=0)
        ref_index = self._ref.currentIndex()
        if ref_index < 3:
            selected = values[:, :, ref_index]
            valid = np.isfinite(selected) & (selected > 0.0) & masks[:, ref_index, None]
            count = np.sum(valid, axis=0)
            module_means = safe_divide(
                np.sum(np.where(valid, selected, 0.0), axis=0), count, count > 0,
            )
        else:
            valid = (
                np.isfinite(values) & (values > 0.0) & (values != 1.0)
                & masks[:, None, :]
            )
            count = np.sum(valid, axis=(0, 2))
            module_means = safe_divide(
                np.sum(np.where(valid, values, 0.0), axis=(0, 2)), count, count > 0,
            )

        lines = []
        for label in labels:
            group_values = [
                float(module_means[index])
                for index in self._gain_summary_groups[label]
                if 0 <= index < len(module_means)
            ]
            finite_values = sorted(value for value in group_values if math.isfinite(value))
            numbers = [f"{value:.4f}" for value in finite_values]
            numbers.extend("nan" for value in group_values if not math.isfinite(value))
            lines.append(f"{label}: {' '.join(numbers)}")
        self._gain_summary.setPlainText("\n".join(lines))

    def _invalidate_and_refresh_views(self):
        self._view_context = None
        self._refresh_views()

    def _refresh_views(self):
        values = self._selected_values_for_batch()
        finite = [v for v in values.values() if math.isfinite(v)]
        is_change = self._is_change_quantity()
        self._map.set_change_mode(is_change)
        if finite:
            span = self._change_map_range()
            if is_change:
                lo, hi = -span, span
            else:
                lo, hi = 1.0 - span, 1.0 + span
            self._map.set_values(values)
            self._map.set_range(lo, hi)
        else:
            self._map.set_values({})

        self._refresh_chart()
        self._refresh_run_average_monitor()
        self._refresh_gain_summary()

    def _set_x_axis_mode(self, mode: str):
        if mode not in {"batch", "unix", "minutes"}:
            return
        self._x_axis_mode = mode
        self._refresh_chart()

    def _check_gain_drop_warning(self):
        runs_data = self._visible_runs_data()
        if not runs_data:
            return
        latest_run = max(runs_data)
        latest_data = runs_data[latest_run]
        if latest_data.nbatches < 2:
            return
        threshold = self._drop_warn.value() / 100.0
        if threshold <= 0.0:
            return

        try:
            latest_event = int(latest_data.event_end[-1])
        except Exception:
            latest_event = latest_data.nbatches
        key = (latest_run, latest_event, round(threshold, 5))
        if key == self._last_gain_drop_warning_key:
            return

        _, _, masks, _ = self._get_view_context()
        change = self._short_change_array(runs_data, masks)
        if change is None:
            return

        valid = np.isfinite(change)
        count = np.sum(valid, axis=1)
        module_change = safe_divide(
            np.sum(np.where(valid, change, 0.0), axis=1), count, count > 0,
        )
        dropped = [
            (f"W{i + 1}", float(value))
            for i, value in enumerate(module_change)
            if math.isfinite(float(value)) and value <= -threshold
        ]
        if not dropped:
            return

        self._last_gain_drop_warning_key = key
        dropped.sort(key=lambda item: item[1])
        names = ", ".join(name for name, _ in dropped[:30])
        more = "" if len(dropped) <= 30 else f"\n... and {len(dropped) - 30} more"
        message = (
            f"{len(dropped)} W module(s) have short-term gain drop > {threshold * 100:.1f}%.\n\n"
            f"{names}{more}"
        )
        self._append(message.replace("\n", " "), "warn")
        QMessageBox.warning(self, "Gain Drop Warning", message)

    def _refresh_chart(self):
        runs_data, centers, masks_by_run, arrays_by_run = self._get_view_context()
        if not runs_data:
            self._chart.set_data("No data", [], [])
            self._ratio_chart.set_data([], [[], [], []], [[], [], []], [])
            return
        idx = self._module_index(self._selected_module)
        if idx is None:
            self._chart.set_data("Select a W module", [], [])
            return
        x: List[float] = []
        run_spans: List[Tuple[int, float, float]] = []
        per_ref_values = [[], [], []]
        ref_ratio_values = [[], [], []]
        ref_ratio_ok = [[], [], []]
        avg_values: List[float] = []
        x_mode = self._x_axis_mode
        valid_unix_times = [
            float(value)
            for data in runs_data.values()
            for value in np.asarray(data.unix_time).reshape(-1)
            if math.isfinite(float(value)) and float(value) > 0.0
        ]
        unix_origin = min(valid_unix_times) if valid_unix_times else math.nan
        offset = 0.0
        for run, data in sorted(runs_data.items()):
            arr = arrays_by_run[run]
            start = offset
            run_x: List[float] = []
            for b in range(data.nbatches):
                if x_mode == "batch":
                    x_value = offset
                else:
                    unix_value = float(data.unix_time[b]) if b < len(data.unix_time) else 0.0
                    if not math.isfinite(unix_value) or unix_value <= 0.0:
                        x_value = math.nan
                    elif x_mode == "minutes":
                        x_value = (unix_value - unix_origin) / 60.0
                    else:
                        x_value = unix_value
                x.append(x_value)
                if math.isfinite(x_value):
                    run_x.append(x_value)
                mask = masks_by_run.get(run, [])[b] if b < len(masks_by_run.get(run, [])) else [False, False, False]
                for j in range(3):
                    ref_ratio_values[j].append(float(data.ref_ratio[b, j]))
                    ref_ratio_ok[j].append(mask[j])
                    per_ref_values[j].append(float(arr[b, idx, j]) if mask[j] else math.nan)
                vals = [arr[b, idx, j] for j in range(3) if mask[j]]
                avg_values.append(finite_avg(vals) if self._is_change_quantity() else positive_avg(vals))
                offset += 1.0
            end = offset
            if x_mode == "batch" and end > start:
                run_spans.append((run, start - 0.5, end - 0.5))
            elif run_x:
                run_spans.append((run, min(run_x), max(run_x)))
        ref_idx = self._ref.currentIndex()
        series = []
        if ref_idx < 3:
            series.append((REF_CHOICES[ref_idx], per_ref_values[ref_idx], REF_COLORS[ref_idx], False))
        else:
            for j in range(3):
                series.append((REF_CHOICES[j], per_ref_values[j], REF_COLORS[j], False))
            series.append(("Avg", avg_values, REF_COLORS[3], True))
        qlabel = QUANTITIES[self._quantity.currentIndex()][1]
        if self._norm_first5.isChecked() and not self._is_change_quantity():
            qlabel += " (first 5 avg = 1)"
        span = self._change_map_range()
        default_y = (-span, span) if self._is_change_quantity() else (1.0 - span, 1.0 + span)
        if x_mode == "unix":
            xlabel = "New York date / time (ET)"
        elif x_mode == "minutes":
            xlabel = "minutes from first valid batch"
        else:
            batch_seconds = self._batch.value() / 10.0
            batch_minutes = batch_seconds / 60.0
            if batch_minutes >= 1.0:
                xlabel = f"batch id (~{batch_minutes:.1f}min)"
            else:
                xlabel = f"batch id (~{batch_seconds:.0f}s)"
        self._chart.set_data(
            f"{self._selected_module} {qlabel}",
            x, series, run_spans, default_y, xlabel, x_mode,
        )
        self._ratio_chart.set_data(
            x, ref_ratio_values, ref_ratio_ok, run_spans, centers, x_mode,
        )

    def closeEvent(self, event):
        self._timer.stop()
        self._batch_reanalyze_timer.stop()
        self._maintenance_timer.stop()
        self._log_flush_timer.stop()
        self._flush_log()
        if self._download_process.state() != QProcess.ProcessState.NotRunning:
            self._download_process.kill()
            self._download_process.waitForFinished(2000)
        if self._replay_process.state() != QProcess.ProcessState.NotRunning:
            self._replay_process.kill()
            self._replay_process.waitForFinished(2000)
        if self._reanalyze_process.state() != QProcess.ProcessState.NotRunning:
            self._reanalyze_process.kill()
            self._reanalyze_process.waitForFinished(2000)
        if self._root_load_thread is not None:
            self._root_load_thread.requestInterruption()
            self._root_load_thread.quit()
            self._root_load_thread.wait()
        super().closeEvent(event)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Online PyQt6 LMS gain monitor")
    parser.add_argument("--theme", choices=available_themes(), default=DEFAULT_THEME)
    parser.add_argument("run", nargs="*", help="initial run number(s)")
    args = parser.parse_args(argv)

    set_theme(args.theme)
    app = QApplication(sys.argv[:1])
    win = OnlineGainMonitor()
    if args.run:
        win._run_edit.setText(" ".join(args.run))
        win._update_paths_from_run()
        win._load_root()
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
