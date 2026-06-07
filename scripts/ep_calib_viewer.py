#!/usr/bin/env python3
"""
ep_calib_viewer.py  —  epCalib output visualization GUI
========================================================
Reads Physics_calib/{run}/calib_iter{N}.json  +  CalibResult_iter{N}.root

ROOT file layout assumed:
  module_energy/h_{name}   — per-module energy TH1F
  ratio_all                — ratio distribution TH1F
  h2_energy_theta          — 2D energy-theta TH2F
  h2_Nevents_moduleMap     — 2D (col, row) event map TH2F
  measured_peak            — distribution of measured peak values
  recon_sigma              — distribution of sigma values

Usage:
  python scripts/ep_calib_viewer.py [Physics_calib_dir]
  python scripts/ep_calib_viewer.py build/Physics_calib
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from PyQt6.QtCore import Qt, QRectF, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QFont
from PyQt6.QtWidgets import (
    QApplication, QButtonGroup, QCheckBox, QComboBox, QFileDialog,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton,
    QSizePolicy, QSplitter, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

try:
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ── find hycal_geoview  ───────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from hycal_geoview import (          # noqa: E402  (after sys.path tweak)
    HyCalMapWidget, Module, load_modules, cmap_qcolor,
    PALETTE_NAMES, apply_theme_palette, set_theme, available_themes, THEME,
    ColorRangeControl,
)

# ── database path  ────────────────────────────────────────────────────────────
_DB_CANDIDATES = [
    _SCRIPT_DIR.parent / "database" / "hycal_map.json",
    _SCRIPT_DIR.parent / "build"    / "database" / "hycal_map.json",
]
MODULES_JSON = next((p for p in _DB_CANDIDATES if p.is_file()), None)


# =============================================================================
#  Data structures
# =============================================================================

class ModuleMetrics:
    """Per-module computed metrics for one (run, iteration) pair."""
    __slots__ = ("name", "stats", "peak", "sigma", "chi2", "amp",
                 "factor", "base_energy", "old_factor")

    def __init__(self, name: str):
        self.name        = name
        self.stats       = 0.0    # event count (histogram integral)
        self.peak        = 0.0    # measured peak (MeV)
        self.sigma       = 0.0    # Gaussian sigma (MeV)
        self.chi2        = 0.0    # chi2/ndf of Gaussian fit
        self.amp         = 0.0    # fitted Gaussian amplitude (counts)
        self.factor      = 1.0    # calibration factor from JSON (after this iteration)
        self.base_energy = 0.0    # expected peak from JSON (MeV)
        self.old_factor  = 0.0    # factor before this iteration (from .dat oldFactor column)


class IterData:
    """All data for one (run, iteration) pair."""
    def __init__(self):
        self.metrics:        Dict[str, ModuleMetrics] = {}
        self.root_path:      Optional[Path]           = None
        self.factors:        Dict[str, dict]          = {}   # raw JSON entries
        self.expected_peaks: Dict[str, float]         = {}   # from .dat ExpectedPeak
        self.old_factors:    Dict[str, float]         = {}   # from .dat oldFactor
        self.modified_modules: set = set()            # modules manually edited this session
        self._global_hists:  dict = {}   # preloaded by worker: ratio_all etc.

    @property
    def json_path(self) -> Optional[Path]:
        """Derive the calib_iter{N}.json path from root_path."""
        if self.root_path is None:
            return None
        return self.root_path.with_name(
            self.root_path.name
                .replace("CalibResult_iter", "calib_iter")
                .replace(".root", ".json")
        )

    @property
    def dat_path(self) -> Optional[Path]:
        """Derive the fitting_parameters_iter{N}.dat path from root_path."""
        if self.root_path is None:
            return None
        return self.root_path.with_name(
            self.root_path.name
                .replace("CalibResult_iter", "fitting_parameters_iter")
                .replace(".root", ".dat")
        )


# =============================================================================
#  I/O helpers
# =============================================================================

def scan_calib_dir(base: Path) -> Dict[str, Dict[int, IterData]]:
    """Return {run_str: {iter_num: IterData}} discovered under *base*."""
    result: Dict[str, Dict[int, IterData]] = {}
    if not base.is_dir():
        return result
    for run_dir in sorted(base.iterdir()):
        if not run_dir.is_dir():
            continue
        iters: Dict[int, IterData] = {}
        for jf in sorted(run_dir.glob("calib_iter*.json")):
            m = re.search(r"calib_iter(\d+)\.json", jf.name)
            if not m:
                continue
            it   = int(m.group(1))
            data = IterData()
            with open(jf) as f:
                entries = json.load(f)
            data.factors = {e["name"]: e for e in entries}
            rf = run_dir / f"CalibResult_iter{it}.root"
            if rf.is_file():
                data.root_path = rf
            # parse fitting_parameters_iter{N}.dat for ExpectedPeak
            df = run_dir / f"fitting_parameters_iter{it}.dat"
            if df.is_file():
                data.expected_peaks = _parse_dat_expected_peaks(df)
                data.old_factors    = _parse_dat_old_factors(df)
            iters[it] = data
        if iters:
            result[run_dir.name] = iters
    return result


def _parse_dat_expected_peaks(dat_path: Path) -> Dict[str, float]:
    """Parse fitting_parameters_iter{N}.dat and return {module: ExpectedPeak}."""
    result: Dict[str, float] = {}
    try:
        for line in dat_path.read_text().splitlines():
            tokens = line.split()
            # skip header and blank lines
            if len(tokens) < 2 or tokens[0] == "Module":
                continue
            try:
                result[tokens[0]] = float(tokens[1])
            except ValueError:
                pass
    except Exception:
        pass
    return result


def _parse_dat_old_factors(dat_path: Path) -> Dict[str, float]:
    """Parse fitting_parameters_iter{N}.dat and return {module: oldFactor} (column 3)."""
    result: Dict[str, float] = {}
    try:
        for line in dat_path.read_text().splitlines():
            tokens = line.split()
            # columns: Module ExpectedPeak MeasuredPeak oldFactor Ratio Sigma Chi2/ndf
            if len(tokens) < 4 or tokens[0] == "Module":
                continue
            try:
                result[tokens[0]] = float(tokens[3])
            except ValueError:
                pass
    except Exception:
        pass
    return result


def _gaussian(x, amp, mu, sigma):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _fit_histogram(counts: np.ndarray, edges: np.ndarray
                   ) -> Tuple[float, float, float, float]:
    """Return (peak_MeV, sigma_MeV, chi2_ndf, amp).

    Mirrors the iterative algorithm in PhysicsTools::FitPeakResolution():
      1. Find histogram maximum as initial center.
      2. Weighted mean within [center ± 3σ]  (σ from 3.5%/√(E/1000) resolution).
      3. First Gaussian fit in [mean ± 2σ]  (σ re-estimated from new mean).
      4. Second (final) Gaussian fit in [mean ± 1σ]  (σ re-estimated again).
    Falls back gracefully when scipy is unavailable or fits diverge.
    """
    centers = 0.5 * (edges[:-1] + edges[1:])
    total   = counts.sum()
    if total < 5:
        return 0.0, 0.0, 0.0, 0.0

    RESOLUTION = 0.035  # 3.5% / sqrt(E/1000)

    def estimate_sigma(E: float) -> float:
        return E * RESOLUTION / math.sqrt(max(E, 1.0) / 1000.0)

    # Step 1: histogram maximum as initial center
    center = float(centers[int(np.argmax(counts))])

    # Step 2: weighted mean within [center ± 3σ]
    sigma = estimate_sigma(center)
    mask  = (centers >= center - 3.0 * sigma) & (centers <= center + 3.0 * sigma)
    w = counts[mask]
    x = centers[mask]
    mean = float((x * w).sum() / w.sum()) if w.sum() > 0 else center

    if not HAS_SCIPY:
        return mean, estimate_sigma(mean), 0.0, 0.0

    def _gauss_fit(mu0: float, half_width: float):
        """Fit Gaussian in [mu0 ± half_width]; return (mu, sigma, chi2_ndf, amp) or None."""
        lo, hi = mu0 - half_width, mu0 + half_width
        m = (centers >= lo) & (centers <= hi) & (counts > 0)
        if m.sum() < 4:
            return None
        xf = centers[m]
        yf = counts[m].astype(float)
        sig0 = max(half_width / 3.0, float(edges[1] - edges[0]))
        try:
            popt, _ = curve_fit(
                _gaussian, xf, yf,
                p0=[float(yf.max()), mu0, sig0],
                bounds=([0, lo, float(edges[1] - edges[0])],
                        [float(yf.max()) * 5, hi, half_width * 2]),
                maxfev=5000,
            )
            amp_f, mu_f, sig_f = float(popt[0]), float(popt[1]), abs(float(popt[2]))
            resid = yf - _gaussian(xf, *popt)
            chi2  = float((resid ** 2 / np.maximum(yf, 1.0)).sum())
            ndf   = max(len(xf) - 3, 1)
            return mu_f, sig_f, chi2 / ndf, amp_f
        except Exception:
            return None

    # Step 3: first Gaussian fit within [mean ± 2σ]
    sigma = estimate_sigma(mean)
    r1 = _gauss_fit(mean, 2.0 * sigma)
    if r1 is not None:
        mean = r1[0]

    # Step 4: final Gaussian fit within [mean ± 1σ]
    sigma = estimate_sigma(mean)
    r2 = _gauss_fit(mean, sigma)
    if r2 is not None:
        return r2
    if r1 is not None:
        return r1
    return mean, estimate_sigma(mean), 0.0, 0.0


def compute_metrics(data: IterData) -> None:
    """Fill *data.metrics* from ROOT histograms + JSON factors."""
    factors = data.factors

    # --- histogram-based metrics ---
    if HAS_UPROOT and data.root_path is not None:
        with uproot.open(str(data.root_path)) as f:
            # collect available module_energy histogram names
            me_keys = {k.split(";")[0].replace("module_energy/", ""): k
                       for k in f.keys() if k.startswith("module_energy/")}

            for hname, raw_key in me_keys.items():
                # hname like "h_W1" -> module name "W1"
                mod_name = hname[len("h_"):]
                counts, edges = f[raw_key].to_numpy()
                mm = data.metrics.setdefault(mod_name, ModuleMetrics(mod_name))
                mm.stats = float(counts.sum())
                mm.peak, mm.sigma, mm.chi2, mm.amp = _fit_histogram(counts, edges)

    # --- JSON factor data ---
    for name, entry in factors.items():
        mm = data.metrics.setdefault(name, ModuleMetrics(name))
        mm.factor = float(entry.get("factor") or 1.0)
        # Prefer ExpectedPeak from .dat (per-run, per-module value);
        # fall back to base_energy in JSON only if .dat was not loaded.
        if name in data.expected_peaks:
            mm.base_energy = data.expected_peaks[name]
        else:
            mm.base_energy = float(entry.get("base_energy") or 0.0)
        # oldFactor from .dat (the factor before this iteration)
        if name in data.old_factors:
            mm.old_factor = data.old_factors[name]
        else:
            mm.old_factor = 0.0  # not available (e.g. no .dat file)

    # --- compute delta_E = (peak - base_energy) / base_energy ---
    # handled on-demand in _refresh_map()


# =============================================================================
#  Map widget
# =============================================================================

class CalibMapWidget(HyCalMapWidget):
    """HyCalMapWidget specialised for calibration quality display."""

    MODE_DELTA_E = "ΔE/E"
    MODE_SIGMA   = "σ (MeV)"
    MODE_CHI2    = "χ²/ndf"
    MODE_STATS   = "Statistics"
    MODE_FACTOR  = "Calib Factor"

    multiSelectionChanged = pyqtSignal(set)  # emitted when selected set changes

    def __init__(self, parent=None):
        super().__init__(parent, enable_zoom_pan=True, min_size=(460, 460))
        self._label = ""
        self._marked_modules: set = set()    # red ✕ (manually edited)
        self._selected_modules: set = set()  # blue ○ (multi-select pending)
        self._multi_select_mode: bool = False

    # ── marked (red ✕) ──────────────────────────────────────────────────────

    def set_marked_modules(self, names: set) -> None:
        self._marked_modules = set(names)
        self.update()

    # ── multi-select (blue ○) ────────────────────────────────────────────────

    def set_multi_select_mode(self, enabled: bool) -> None:
        self._multi_select_mode = enabled
        if not enabled:
            self._selected_modules.clear()
            self.multiSelectionChanged.emit(set())
        self.update()

    def clear_selection(self) -> None:
        self._selected_modules.clear()
        self.multiSelectionChanged.emit(set())
        self.update()

    def get_selected_modules(self) -> set:
        return set(self._selected_modules)

    def _handle_click(self, pos):
        if self._check_inline_range_edit_click(pos):
            return
        if self._cb_rect and self._cb_rect.contains(pos):
            self.paletteClicked.emit()
            return
        name = self._hit(pos)
        if self._multi_select_mode and name:
            if name in self._selected_modules:
                self._selected_modules.discard(name)
            else:
                self._selected_modules.add(name)
            self.multiSelectionChanged.emit(set(self._selected_modules))
            self.update()
        else:
            self.moduleClicked.emit(name or "")

    # ── tooltip ─────────────────────────────────────────────────────────────

    def set_map_label(self, label: str):
        self._label = label

    def _colorbar_center_text(self) -> str:
        return self._label or PALETTE_NAMES[self._palette_idx]

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        base = name if v is None else f"{name}: {self._fmt_value(v)}"
        tags = []
        if name in self._marked_modules:
            tags.append("✕ edited")
        if name in self._selected_modules:
            tags.append("◯ selected")
        if tags:
            base += f"  [{', '.join(tags)}]"
        return base

    # ── overlays ─────────────────────────────────────────────────────────────

    def _paint_overlays(self, p, w, h):
        super()._paint_overlays(p, w, h)
        from PyQt6.QtGui import QPen, QColor
        from PyQt6.QtCore import Qt
        # red ✕ for manually edited modules
        if self._marked_modules:
            pen = QPen(QColor("#ff3b30"), 1.5)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            p.setPen(pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            for name in self._marked_modules:
                rect = self._rects.get(name)
                if rect is None:
                    continue
                m = max(rect.width(), rect.height()) * 0.18
                x0, y0 = rect.x() + m, rect.y() + m
                x1, y1 = rect.right() - m, rect.bottom() - m
                p.drawLine(int(x0), int(y0), int(x1), int(y1))
                p.drawLine(int(x1), int(y0), int(x0), int(y1))
        # blue ○ for multi-selected modules
        if self._selected_modules:
            pen = QPen(QColor("#3b9eff"), 2.0)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            p.setPen(pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            for name in self._selected_modules:
                rect = self._rects.get(name)
                if rect is None:
                    continue
                m = max(rect.width(), rect.height()) * 0.08
                p.drawEllipse(
                    int(rect.x() + m), int(rect.y() + m),
                    int(rect.width() - 2 * m), int(rect.height() - 2 * m)
                )


# =============================================================================
#  Matplotlib canvas
# =============================================================================

class MplCanvas(FigureCanvas):
    def __init__(self, width=4, height=3, dpi=90):
        fc = THEME.CANVAS
        self.fig = Figure(figsize=(width, height), dpi=dpi,
                          facecolor=fc, tight_layout=True)
        self.ax  = self.fig.add_subplot(111)
        self.ax.set_facecolor(fc)
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

    def _style_ax(self, title="", xlabel="", ylabel=""):
        ax = self.ax
        ax.set_title(title, color=THEME.TEXT, fontsize=12)
        ax.set_xlabel(xlabel, color=THEME.TEXT_DIM, fontsize=10)
        ax.set_ylabel(ylabel, color=THEME.TEXT_DIM, fontsize=10)
        ax.tick_params(colors=THEME.TEXT_DIM, labelsize=9)
        for sp in ax.spines.values():
            sp.set_edgecolor(THEME.BORDER)

    def clear_ax(self):
        self.ax.cla()
        self.ax.set_facecolor(THEME.CANVAS)

    def redraw(self):
        self.fig.canvas.draw_idle()


# =============================================================================
#  Module detail panel
# =============================================================================

class ModuleDetailPanel(QWidget):
    refitApplied = pyqtSignal(str)   # emitted with module name after Apply+Save

    def __init__(self, parent=None):
        super().__init__(parent)
        self._root_path:       Optional[Path]          = None
        self._mm:              Optional[ModuleMetrics] = None
        self._iter_data:       Optional[IterData]      = None
        self._cur_module_name: str   = ""
        self._hist_cache: "OrderedDict[str, Tuple[np.ndarray, np.ndarray]]" = OrderedDict()
        self._refit_peak:      float = 0.0
        self._refit_sigma:     float = 0.0
        self._refit_chi2:      float = 0.0
        self._refit_new_factor:float = 0.0
        self._span_selector:   object = None   # matplotlib SpanSelector instance
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        self._name_lbl = QLabel("← click a module on the map")
        self._name_lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        self._name_lbl.setStyleSheet(f"color: {THEME.ACCENT};")
        lay.addWidget(self._name_lbl)

        self._info_lbl = QLabel("")
        self._info_lbl.setFont(QFont("Monospace", 11))
        self._info_lbl.setStyleSheet(f"color: {THEME.TEXT};")
        self._info_lbl.setWordWrap(True)
        lay.addWidget(self._info_lbl)

        self._canvas = MplCanvas(width=4, height=3)
        lay.addWidget(self._canvas)

        # ── Manual re-fit controls ──────────────────────────────────────
        rf_box = QGroupBox("Manual Gauss Re-fit")
        rf_lay = QVBoxLayout(rf_box)
        rf_lay.setContentsMargins(6, 12, 6, 4)
        rf_lay.setSpacing(4)

        range_row = QHBoxLayout()
        range_row.setSpacing(4)
        range_row.addWidget(QLabel("Fit range:"))
        self._rf_xmin = QLineEdit()
        self._rf_xmin.setFixedWidth(72)
        self._rf_xmin.setPlaceholderText("Xmin")
        range_row.addWidget(self._rf_xmin)
        range_row.addWidget(QLabel("–"))
        self._rf_xmax = QLineEdit()
        self._rf_xmax.setFixedWidth(72)
        self._rf_xmax.setPlaceholderText("Xmax")
        range_row.addWidget(self._rf_xmax)
        range_row.addSpacing(16)
        range_row.addWidget(QLabel("Init μ₀:"))
        self._rf_mu0 = QLineEdit()
        self._rf_mu0.setFixedWidth(72)
        self._rf_mu0.setPlaceholderText("μ₀ (MeV)")
        range_row.addWidget(self._rf_mu0)
        range_row.addSpacing(16)
        range_row.addWidget(QLabel("On drag:"))
        self._span_action = QComboBox()
        self._span_action.addItems(["Fill range only", "Auto Refit", "Auto Mean"])
        self._span_action.setToolTip(
            "Action to perform after dragging a range on the histogram:\n"
            "  Fill range only — populate Xmin/Xmax/\u03bc\u2080 without triggering a calculation\n"
            "  Auto Refit     — immediately run Gaussian refit\n"
            "  Auto Mean      — immediately compute weighted mean")
        range_row.addWidget(self._span_action)
        range_row.addStretch()
        rf_lay.addLayout(range_row)

        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        self._refit_btn = QPushButton("Run Refit")
        self._refit_btn.setEnabled(False)
        self._refit_btn.clicked.connect(self._do_refit)
        btn_row.addWidget(self._refit_btn)
        self._mean_btn = QPushButton("Use Mean (χ²=1)")
        self._mean_btn.setEnabled(False)
        self._mean_btn.setToolTip(
            "Compute weighted mean/σ from histogram (in fit range if set)\n"
            "and use as peak value; χ²/ndf is fixed at 1.0")
        self._mean_btn.clicked.connect(self._do_use_mean)
        btn_row.addWidget(self._mean_btn)
        self._apply_btn = QPushButton("Apply && Save JSON")
        self._apply_btn.setEnabled(False)
        self._apply_btn.clicked.connect(self._apply_and_save)
        btn_row.addWidget(self._apply_btn)
        btn_row.addStretch()
        rf_lay.addLayout(btn_row)

        self._refit_status = QLabel("")
        self._refit_status.setFont(QFont("Monospace", 10))
        self._refit_status.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        self._refit_status.setWordWrap(True)
        rf_lay.addWidget(self._refit_status)

        lay.addWidget(rf_box)

        # ── Set Factor directly ─────────────────────────────────────────
        sf_box = QGroupBox("Set Factor Directly")
        sf_lay = QHBoxLayout(sf_box)
        sf_lay.setContentsMargins(6, 12, 6, 4)
        sf_lay.setSpacing(6)
        sf_lay.addWidget(QLabel("Preset:"))
        for preset in (0.122, 0.15):
            btn = QPushButton(f"{preset}")
            btn.setFixedWidth(60)
            btn.setToolTip(f"Set factor to {preset} and save JSON")
            btn.clicked.connect(lambda _, v=preset: self._do_set_factor(v))
            sf_lay.addWidget(btn)
        sf_lay.addSpacing(16)
        sf_lay.addWidget(QLabel("Custom:"))
        self._sf_custom = QLineEdit()
        self._sf_custom.setFixedWidth(80)
        self._sf_custom.setPlaceholderText("factor")
        sf_lay.addWidget(self._sf_custom)
        self._sf_apply_btn = QPushButton("Set && Save")
        self._sf_apply_btn.setToolTip("Apply the custom factor value and save JSON")
        self._sf_apply_btn.clicked.connect(self._do_set_factor_custom)
        sf_lay.addWidget(self._sf_apply_btn)
        sf_lay.addSpacing(16)
        self._restore_btn = QPushButton("↩ Restore Old Factor")
        self._restore_btn.setToolTip("Write old_factor (from .dat) back to JSON and keep the mark")
        self._restore_btn.clicked.connect(self._do_restore_old_factor)
        sf_lay.addWidget(self._restore_btn)
        self._clear_mark_btn = QPushButton("✕ Clear Mark")
        self._clear_mark_btn.setToolTip("Remove the red-cross mark from this module (does not change factor)")
        self._clear_mark_btn.clicked.connect(self._do_clear_mark)
        sf_lay.addWidget(self._clear_mark_btn)
        sf_lay.addStretch()
        lay.addWidget(sf_box)

    def set_root_path(self, path: Optional[Path]):
        self._root_path = path

    def show_module(self, name: str, mm: Optional[ModuleMetrics]):
        self._mm = mm
        self._cur_module_name  = name
        self._refit_peak       = 0.0
        self._refit_new_factor = 0.0
        self._apply_btn.setEnabled(False)
        self._refit_status.setText("")
        self._refit_btn.setEnabled(bool(name) and HAS_SCIPY and HAS_UPROOT)
        self._mean_btn.setEnabled(bool(name) and HAS_UPROOT)

        # Pre-populate fit range from current fit result
        if mm and mm.peak > 0 and mm.sigma > 0:
            self._rf_xmin.setText(f"{mm.peak - 3 * mm.sigma:.1f}")
            self._rf_xmax.setText(f"{mm.peak + 3 * mm.sigma:.1f}")
            self._rf_mu0.setText(f"{mm.peak:.1f}")
        elif mm and mm.base_energy > 0:
            self._rf_mu0.setText(f"{mm.base_energy:.1f}")

        self._name_lbl.setText(f"Module: {name}")
        if mm:
            delta_e = ""
            if mm.base_energy > 0 and mm.peak > 0:
                de = (mm.peak - mm.base_energy) / mm.base_energy * 100
                delta_e = f"  ΔE/E={de:+.2f}%"
            self._info_lbl.setText(
                f"Factor: {mm.factor:.4f}   Base energy: {mm.base_energy:.1f} MeV\n"
                f"Events: {mm.stats:.0f}   "
                f"Peak: {mm.peak:.1f} MeV   σ: {mm.sigma:.1f} MeV   "
                f"χ²/ndf: {mm.chi2:.3f}{delta_e}"
            )
        else:
            self._info_lbl.setText("")

        self._canvas.clear_ax()
        ax = self._canvas.ax

        hist_data = self._read_module_hist(name)
        if hist_data is not None:
            counts, edges = hist_data
            centers = 0.5 * (edges[:-1] + edges[1:])
            # crop to non-zero range + some margin
            nz = np.where(counts > 0)[0]
            if nz.size:
                lo_i = max(0, nz[0] - 5)
                hi_i = min(len(counts), nz[-1] + 6)
                counts  = counts[lo_i:hi_i]
                centers = centers[lo_i:hi_i]
                edges_c = edges[lo_i:hi_i + 1]
            else:
                edges_c = edges

            ax.bar(centers, counts, width=np.diff(edges_c),
                   color=THEME.ACCENT, alpha=0.75, linewidth=0)

            # overlay Gaussian fit
            if mm and mm.peak > 0 and mm.sigma > 0 and HAS_SCIPY:
                x   = np.linspace(edges_c[0], edges_c[-1], 500)
                bw  = edges_c[1] - edges_c[0] if len(edges_c) > 1 else 1.0
                amp = mm.amp if mm.amp > 0 else (
                    counts.sum() * bw / (mm.sigma * math.sqrt(2 * math.pi)))
                gaus = amp * np.exp(-0.5 * ((x - mm.peak) / mm.sigma) ** 2)
                ax.plot(x, gaus, color=THEME.WARN, linewidth=1.5,
                        label=f"Gauss fit  μ={mm.peak:.0f}")

            if mm and mm.base_energy > 0:
                ax.axvline(mm.base_energy, color=THEME.SUCCESS,
                           linewidth=1.2, linestyle="--",
                           label=f"Expected {mm.base_energy:.0f} MeV")

            ax.legend(fontsize=9, labelcolor=THEME.TEXT,
                      facecolor=THEME.PANEL, edgecolor=THEME.BORDER)
        else:
            msg = ("uproot not installed\npip install uproot"
                   if not HAS_UPROOT else "No histogram found in ROOT file")
            ax.text(0.5, 0.5, msg, transform=ax.transAxes,
                    ha="center", va="center",
                    color=THEME.TEXT_DIM, fontsize=9)

        self._canvas._style_ax(f"{name}  energy histogram",
                                "Energy (MeV)", "Counts")
        self._canvas.redraw()
        self._setup_span_selector()

    # ── Interactive span selector ──────────────────────────────────────────

    def _setup_span_selector(self) -> None:
        """Attach (or re-attach) a SpanSelector to the current axes.

        The selector is re-created every time the histogram is redrawn so that
        it stays bound to the live Axes object.
        """
        # Disconnect the old selector to avoid ghost callbacks
        if self._span_selector is not None:
            try:
                self._span_selector.disconnect_events()
            except Exception:
                pass
            self._span_selector = None

        ax = self._canvas.ax
        self._span_selector = SpanSelector(
            ax,
            self._on_span_select,
            direction="horizontal",
            useblit=True,
            props=dict(alpha=0.25, facecolor=THEME.WARN),
            interactive=True,
            drag_from_anywhere=True,
        )

    def _on_span_select(self, xmin: float, xmax: float) -> None:
        """Called by SpanSelector after the user finishes dragging."""
        if xmax - xmin < 1.0:   # ignore accidental single-clicks
            return
        mu0 = 0.5 * (xmin + xmax)
        self._rf_xmin.setText(f"{xmin:.1f}")
        self._rf_xmax.setText(f"{xmax:.1f}")
        self._rf_mu0.setText(f"{mu0:.1f}")

        action = self._span_action.currentIndex()
        if action == 1:    # Auto Refit
            self._do_refit()
        elif action == 2:  # Auto Mean
            self._do_use_mean()

    _HIST_CACHE_MAX = 30

    def _read_module_hist(self, name: str
                          ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        # Fast path: bounded LRU cache (cleared when iteration changes)
        if name in self._hist_cache:
            self._hist_cache.move_to_end(name)
            return self._hist_cache[name]
        if not HAS_UPROOT or self._root_path is None:
            return None
        key = f"module_energy/h_{name}"
        try:
            with uproot.open(str(self._root_path)) as f:
                avail = {k.split(";")[0] for k in f.keys()}
                if key not in avail:
                    return None
                result = f[key].to_numpy()
                self._hist_cache[name] = result
                if len(self._hist_cache) > self._HIST_CACHE_MAX:
                    self._hist_cache.popitem(last=False)
                return result
        except Exception:
            return None

    # ──────────────────────────────────────────────────────────────────────

    def set_iter_data(self, data: Optional["IterData"]) -> None:
        """Store reference to current IterData (needed for JSON write-back)."""
        self._iter_data = data
        self._hist_cache.clear()

    def _do_set_factor(self, value: float) -> None:
        """Directly set factor to *value* and write to JSON."""
        name = self._cur_module_name
        if not name or self._iter_data is None:
            self._refit_status.setText("No module selected")
            return
        mm = self._iter_data.metrics.get(name)
        if mm is None:
            self._refit_status.setText(f"Module {name!r} not in metrics")
            return

        old_factor = mm.factor
        mm.factor  = value

        entry = self._iter_data.factors.get(name)
        if entry is not None:
            entry["factor"] = value

        jpath = self._iter_data.json_path
        if jpath is not None:
            try:
                entries = list(self._iter_data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                note = f"  →  saved {jpath.name}"
            except Exception as exc:
                self._refit_status.setText(f"Save failed: {exc}")
                return
        else:
            note = "  (JSON path not found — not saved)"

        self._iter_data.modified_modules.add(name)
        self._mm = mm
        self._info_lbl.setText(
            f"Factor: {mm.factor:.5f}  (was {old_factor:.5f})   "
            f"Base energy: {mm.base_energy:.1f} MeV\n"
            f"Events: {mm.stats:.0f}   "
            f"Peak: {mm.peak:.1f} MeV   σ: {mm.sigma:.1f} MeV   "
            f"χ²/ndf: {mm.chi2:.3f}"
        )
        self._refit_status.setText(f"Factor set to {value:.5f}  (was {old_factor:.5f}){note}")
        self.refitApplied.emit(name)

    def _do_set_factor_custom(self) -> None:
        """Read the custom factor QLineEdit and call _do_set_factor."""
        txt = self._sf_custom.text().strip()
        try:
            value = float(txt)
        except ValueError:
            self._refit_status.setText("Invalid factor value — enter a number")
            return
        if value <= 0:
            self._refit_status.setText("Factor must be positive")
            return
        self._do_set_factor(value)

    def _do_restore_old_factor(self) -> None:
        """Restore old_factor from .dat to JSON and remove the mark."""
        name = self._cur_module_name
        if not name or self._iter_data is None:
            self._refit_status.setText("No module selected")
            return
        mm = self._iter_data.metrics.get(name)
        if mm is None:
            self._refit_status.setText(f"Module {name!r} not in metrics")
            return
        if mm.old_factor <= 0.0:
            self._refit_status.setText("old_factor not available (.dat not loaded)")
            return

        old_cur = mm.factor
        mm.factor = mm.old_factor

        entry = self._iter_data.factors.get(name)
        if entry is not None:
            entry["factor"] = mm.old_factor

        jpath = self._iter_data.json_path
        if jpath is not None:
            try:
                entries = list(self._iter_data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                note = f"  →  saved {jpath.name}"
            except Exception as exc:
                self._refit_status.setText(f"Save failed: {exc}")
                return
        else:
            note = "  (JSON path not found — not saved)"

        self._iter_data.modified_modules.add(name)
        self._mm = mm
        self._info_lbl.setText(
            f"Factor: {mm.factor:.5f}  (restored from {old_cur:.5f})   "
            f"Base energy: {mm.base_energy:.1f} MeV\n"
            f"Events: {mm.stats:.0f}   "
            f"Peak: {mm.peak:.1f} MeV   σ: {mm.sigma:.1f} MeV   "
            f"χ²/ndf: {mm.chi2:.3f}"
        )
        self._refit_status.setText(f"Restored old_factor={mm.old_factor:.5f}{note}")
        self.refitApplied.emit(name)

    def _do_clear_mark(self) -> None:
        """Remove red-cross mark without changing factor."""
        name = self._cur_module_name
        if not name or self._iter_data is None:
            return
        self._iter_data.modified_modules.discard(name)
        self._refit_status.setText(f"Mark cleared for {name}")
        self.refitApplied.emit(name)

    def _do_use_mean(self) -> None:
        """Use histogram weighted mean as peak value; chi2/ndf is fixed at 1.0."""
        name = self._cur_module_name
        if not name or self._root_path is None:
            self._refit_status.setText("No module selected")
            return

        hist_data = self._read_module_hist(name)
        if hist_data is None:
            self._refit_status.setText("Histogram not available in ROOT file")
            return

        counts, edges = hist_data
        centers = 0.5 * (edges[:-1] + edges[1:])

        # Apply fit range if provided
        xmin_txt = self._rf_xmin.text().strip()
        xmax_txt = self._rf_xmax.text().strip()
        try:
            xmin = float(xmin_txt) if xmin_txt else None
        except ValueError:
            xmin = None
        try:
            xmax = float(xmax_txt) if xmax_txt else None
        except ValueError:
            xmax = None

        mask = counts > 0
        if xmin is not None:
            mask &= centers >= xmin
        if xmax is not None:
            mask &= centers <= xmax

        c_sel = counts[mask]
        x_sel = centers[mask]
        total = c_sel.sum()
        if total < 1:
            self._refit_status.setText("No data in selected range")
            return

        mu0  = float((x_sel * c_sel).sum() / total)
        var  = float(((x_sel - mu0) ** 2 * c_sel).sum() / total)
        bw   = float(edges[1] - edges[0]) if len(edges) > 1 else 1.0
        sig0 = math.sqrt(max(var, bw ** 2 / 12.0))

        old_peak   = self._mm.peak   if (self._mm and self._mm.peak > 0.0) else 0.0
        # Use oldFactor from .dat (factor before this iteration);
        # fall back to current JSON factor if .dat was not available.
        old_factor = (self._mm.old_factor if (self._mm and self._mm.old_factor > 0.0)
                      else (self._mm.factor if self._mm else 1.0))
        base_energy = self._mm.base_energy if (self._mm and self._mm.base_energy > 0.0) else 0.0
        # Formula: new_factor = old_factor * expected_peak / refit_peak
        new_factor = old_factor * base_energy / mu0 if (mu0 > 0.0 and base_energy > 0.0) else 0.0

        self._refit_peak       = mu0
        self._refit_sigma      = sig0
        self._refit_chi2       = 1.0
        self._refit_new_factor = new_factor

        range_note = ""
        if xmin is not None or xmax is not None:
            lo_s = f"{xmin:.1f}" if xmin is not None else "-inf"
            hi_s = f"{xmax:.1f}" if xmax is not None else "+inf"
            range_note = f"  [range: {lo_s} – {hi_s}]"

        msg = (
            f"Mean:  μ = {mu0:.2f} MeV     "
            f"σ = {sig0:.2f} MeV     χ²/ndf = 1.000  (fixed){range_note}"
        )
        if new_factor > 0.0:
            msg += f"\nNew factor: {new_factor:.5f}   (oldFactor={old_factor:.5f},  expected={base_energy:.1f} MeV)"
        self._refit_status.setText(msg)
        self._apply_btn.setEnabled(new_factor > 0.0)
        bw_mean = float(edges[1] - edges[0]) if len(edges) > 1 else 1.0
        amp_mean = counts.sum() * bw_mean / (sig0 * math.sqrt(2 * math.pi))
        self._redraw_refit_overlay(name, counts, edges, mu0, sig0, amp_mean)

    def _do_refit(self) -> None:
        """Gaussian refit using user-specified range and initial peak center."""
        if not HAS_SCIPY or not HAS_UPROOT:
            self._refit_status.setText("scipy and uproot required")
            return
        name = self._cur_module_name
        if not name or self._root_path is None:
            self._refit_status.setText("No module selected")
            return
        try:
            xmin = float(self._rf_xmin.text())
            xmax = float(self._rf_xmax.text())
            mu0  = float(self._rf_mu0.text())
        except ValueError:
            self._refit_status.setText("Invalid input — enter numeric values in all three fields")
            return
        if xmin >= xmax:
            self._refit_status.setText("Fit range error: Xmin must be < Xmax")
            return

        hist_data = self._read_module_hist(name)
        if hist_data is None:
            self._refit_status.setText("Histogram not available in ROOT file")
            return

        counts, edges = hist_data
        centers = 0.5 * (edges[:-1] + edges[1:])
        bw = float(edges[1] - edges[0]) if len(edges) > 1 else 1.0

        mask = (centers >= xmin) & (centers <= xmax) & (counts > 0)
        if mask.sum() < 4:
            self._refit_status.setText(
                f"Only {mask.sum()} non-empty bins in fit range — widen the range")
            return

        xf   = centers[mask]
        yf   = counts[mask].astype(float)
        amp0 = float(yf.max())
        sig0 = max((xmax - xmin) / 5.0, bw)

        try:
            popt, pcov = curve_fit(
                _gaussian, xf, yf,
                p0=[amp0, mu0, sig0],
                bounds=([0.0, xmin, bw], [amp0 * 3.0, xmax, xmax - xmin]),
                maxfev=6000,
            )
            amp_f, mu_f, sig_f = popt
            resid    = yf - _gaussian(xf, *popt)
            chi2     = float((resid ** 2 / np.maximum(yf, 1.0)).sum())
            ndf      = max(len(xf) - 3, 1)
            chi2_ndf = chi2 / ndf

            self._refit_peak  = float(mu_f)
            self._refit_sigma = float(abs(sig_f))
            self._refit_chi2  = chi2_ndf

            # Use oldFactor from .dat (factor before this iteration);
            # fall back to current JSON factor if .dat was not available.
            old_factor  = (self._mm.old_factor if (self._mm and self._mm.old_factor > 0.0)
                           else (self._mm.factor if self._mm else 1.0))
            base_energy = self._mm.base_energy if (self._mm and self._mm.base_energy > 0.0) else 0.0
            # Formula: new_factor = old_factor * expected_peak / refit_peak
            self._refit_new_factor = (
                old_factor * base_energy / mu_f if (mu_f > 0.0 and base_energy > 0.0) else 0.0
            )

            msg = (
                f"Refit:  μ = {mu_f:.2f} MeV     "
                f"σ = {abs(sig_f):.2f} MeV     χ²/ndf = {chi2_ndf:.3f}"
            )
            if self._refit_new_factor > 0.0:
                msg += (
                    f"\nNew factor: {self._refit_new_factor:.5f}"
                    f"   (oldFactor={old_factor:.5f},  expected={base_energy:.1f} MeV)"
                )
            self._refit_status.setText(msg)
            self._apply_btn.setEnabled(self._refit_new_factor > 0.0)
            self._redraw_refit_overlay(
                name, counts, edges, float(mu_f), float(abs(sig_f)), float(amp_f))

        except Exception as exc:
            self._refit_status.setText(f"Fit failed: {exc}")
            self._refit_peak       = 0.0
            self._refit_new_factor = 0.0

    def _redraw_refit_overlay(
            self,
            name: str,
            counts: np.ndarray,
            edges: np.ndarray,
            mu_f: float,
            sig_f: float,
            amp_f: float = 0.0,
    ) -> None:
        """Redraw histogram with old (dashed) + new (solid) Gaussian overlaid."""
        self._canvas.clear_ax()
        ax = self._canvas.ax

        centers = 0.5 * (edges[:-1] + edges[1:])
        nz = np.where(counts > 0)[0]
        if nz.size:
            lo_i = max(0, nz[0] - 5)
            hi_i = min(len(counts), nz[-1] + 6)
            c, ce, ee = counts[lo_i:hi_i], centers[lo_i:hi_i], edges[lo_i:hi_i + 1]
        else:
            c, ce, ee = counts, centers, edges

        bw = float(ee[1] - ee[0]) if len(ee) > 1 else 1.0
        ax.bar(ce, c, width=bw, color=THEME.ACCENT, alpha=0.6, linewidth=0)

        x = np.linspace(float(ee[0]), float(ee[-1]), 600)
        # original fit — dim dashed line
        mm = self._mm
        if mm and mm.peak > 0.0 and mm.sigma > 0.0:
            orig_amp = mm.amp if mm.amp > 0.0 else (
                float(c.sum()) * bw / (mm.sigma * math.sqrt(2 * math.pi)))
            ax.plot(x, orig_amp * np.exp(-0.5 * ((x - mm.peak) / mm.sigma) ** 2),
                    color=THEME.TEXT_DIM, linewidth=1.0, linestyle="--",
                    label=f"Original  μ={mm.peak:.0f}")
        # new refit — highlighted solid line
        if amp_f <= 0.0:
            amp_f = float(c.sum()) * bw / (sig_f * math.sqrt(2 * math.pi))
        ax.plot(x, amp_f * np.exp(-0.5 * ((x - mu_f) / sig_f) ** 2),
                color=THEME.WARN, linewidth=2.0,
                label=f"Refit  μ={mu_f:.0f}")

        if mm and mm.base_energy > 0.0:
            ax.axvline(mm.base_energy, color=THEME.SUCCESS,
                       linewidth=1.2, linestyle="--",
                       label=f"Expected {mm.base_energy:.0f} MeV")

        ax.legend(fontsize=9, labelcolor=THEME.TEXT,
                  facecolor=THEME.PANEL, edgecolor=THEME.BORDER)
        self._canvas._style_ax(f"{name}  energy histogram (refit)",
                                "Energy (MeV)", "Counts")
        self._canvas.redraw()
        self._setup_span_selector()

    def _apply_and_save(self) -> None:
        """Update in-memory IterData metrics and write calibration JSON to disk."""
        name = self._cur_module_name
        if (
            not name
            or self._iter_data is None
            or self._refit_peak <= 0.0
            or self._refit_new_factor <= 0.0
        ):
            return
        mm = self._iter_data.metrics.get(name)
        if mm is None:
            return

        old_factor = mm.factor
        mm.peak   = self._refit_peak
        mm.sigma  = self._refit_sigma
        mm.chi2   = self._refit_chi2
        mm.factor = self._refit_new_factor

        # keep factors dict in sync so JSON reflects the update
        entry = self._iter_data.factors.get(name)
        if entry is not None:
            entry["factor"] = self._refit_new_factor

        # write JSON to disk
        jpath = self._iter_data.json_path
        if jpath is not None:
            try:
                entries = list(self._iter_data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                save_note = f"  →  saved {jpath.name}"
            except Exception as exc:
                self._refit_status.setText(f"Save failed: {exc}")
                return
        else:
            save_note = "  (JSON path not found — not saved to disk)"

        # write .dat to disk
        dpath = self._iter_data.dat_path
        dat_note = ""
        if dpath is not None and dpath.is_file():
            try:
                lines = dpath.read_text().splitlines(keepends=True)
                new_lines = []
                updated = False
                for line in lines:
                    tokens = line.split()
                    if tokens and tokens[0] == name:
                        # reconstruct line preserving column widths
                        ratio = (mm.base_energy / mm.peak
                                 if mm.peak > 0.0 else 0.0)
                        new_line = (
                            f"{name:<16}{mm.base_energy:<16.6g}"
                            f"{mm.peak:<16.6g}{ratio:<16.6g}"
                            f"{mm.sigma:<16.6g}{mm.chi2:<16.6g}\n"
                        )
                        new_lines.append(new_line)
                        updated = True
                    else:
                        new_lines.append(line)
                if updated:
                    dpath.write_text("".join(new_lines))
                    dat_note = f", {dpath.name}"
                else:
                    dat_note = f" (module {name!r} not found in .dat)"
            except Exception as exc:
                dat_note = f" (.dat write failed: {exc})"
        elif dpath is not None:
            dat_note = " (.dat file not found)"

        self._mm = mm
        delta_e = ""
        if mm.base_energy > 0 and mm.peak > 0:
            de = (mm.peak - mm.base_energy) / mm.base_energy * 100
            delta_e = f"  ΔE/E={de:+.2f}%"
        self._info_lbl.setText(
            f"Factor: {mm.factor:.5f}  (was {old_factor:.5f})   "
            f"Base energy: {mm.base_energy:.1f} MeV\n"
            f"Events: {mm.stats:.0f}   "
            f"Peak: {mm.peak:.1f} MeV   σ: {mm.sigma:.1f} MeV   "
            f"χ²/ndf: {mm.chi2:.3f}{delta_e}"
        )
        first_line = self._refit_status.text().split("\n")[0]
        self._refit_status.setText(
            first_line + f"\nApplied.{save_note}{dat_note}")
        self._apply_btn.setEnabled(False)
        if self._iter_data is not None:
            self._iter_data.modified_modules.add(name)
        self.refitApplied.emit(name)


# =============================================================================
#  Global stats panel
# =============================================================================

class GlobalStatsPanel(QWidget):

    _RIGHT_VIEWS = ["Energy vs θ", "One Cluster Energy", "Resolution (σ)"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self._root_path: Optional[Path] = None
        self._right_view_idx = 0
        self._cached_hists: dict = {}
        self._build_ui()

    def _build_ui(self):
        lay = QHBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.setSpacing(4)
        self._cv_ratio  = MplCanvas(width=3, height=2.4)
        self._cv_etheta = MplCanvas(width=3, height=2.4)
        self._etheta_cbar = None
        lay.addWidget(self._cv_ratio)

        # right panel: toggle button + canvas
        right_box = QWidget()
        right_lay = QVBoxLayout(right_box)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(2)
        self._right_toggle = QPushButton(self._RIGHT_VIEWS[self._right_view_idx])
        self._right_toggle.setToolTip("Click to switch between views")
        self._right_toggle.clicked.connect(self._cycle_right_view)
        right_lay.addWidget(self._right_toggle)
        right_lay.addWidget(self._cv_etheta)
        lay.addWidget(right_box)

    def _cycle_right_view(self):
        self._right_view_idx = (self._right_view_idx + 1) % len(self._RIGHT_VIEWS)
        self._right_toggle.setText(self._RIGHT_VIEWS[self._right_view_idx])
        self._draw_right()

    def set_root_path(self, path: Optional[Path]):
        """Called when no pre-loaded cache is available (legacy / refresh path)."""
        self._root_path = path
        self._cached_hists.clear()
        if HAS_UPROOT and path is not None:
            try:
                with uproot.open(str(path)) as f:
                    avail = {k.split(";")[0] for k in f.keys()}
                    for key in ("ratio_all", "h2_energy_theta",
                                "one_cluster_energy", "recon_sigma"):
                        if key in avail:
                            try:
                                self._cached_hists[key] = f[key].to_numpy()
                            except Exception:
                                pass
            except Exception:
                pass
        self._draw_ratio()
        self._draw_right()

    def set_preloaded(self, path: Optional[Path], hists: dict) -> None:
        """Use already-loaded histogram cache from the background worker."""
        self._root_path = path
        self._cached_hists = hists
        self._draw_ratio()
        self._draw_right()

    def _draw_right(self):
        idx = self._right_view_idx
        if idx == 0:
            self._draw_etheta()
        elif idx == 1:
            self._draw_one_cluster_energy()
        else:
            self._draw_resolution()

    def _open_root(self):
        if not HAS_UPROOT or self._root_path is None:
            return None
        try:
            return uproot.open(str(self._root_path))
        except Exception:
            return None

    def _draw_ratio(self):
        ax = self._cv_ratio.ax
        self._cv_ratio.clear_ax()
        data = self._cached_hists.get("ratio_all")
        if data is not None:
            try:
                counts, edges = data
                centers = 0.5 * (edges[:-1] + edges[1:])
                nz = counts > 0
                if nz.any():
                    ax.bar(centers[nz], counts[nz],
                           width=(edges[1] - edges[0]),
                           color=THEME.SUCCESS, alpha=0.8, linewidth=0)
                ax.axvline(1.0, color=THEME.WARN, linewidth=1.2,
                           linestyle="--", label="ratio = 1")
                lo = max(0, centers[nz][0] - 0.3) if nz.any() else 0
                hi = min(4, centers[nz][-1] + 0.3) if nz.any() else 4
                ax.set_xlim(lo, hi)
                ax.legend(fontsize=9, labelcolor=THEME.TEXT,
                          facecolor=THEME.PANEL, edgecolor=THEME.BORDER)
            except Exception:
                pass
        else:
            ax.text(0.5, 0.5, "ratio_all\nnot found",
                    transform=ax.transAxes, ha="center", va="center",
                    color=THEME.TEXT_DIM, fontsize=9)
        self._cv_ratio._style_ax("Peak ratio distribution",
                                  "Ratio (exp / meas)", "Modules")
        self._cv_ratio.redraw()

    def _draw_etheta(self):
        ax = self._cv_etheta.ax
        # remove previous colorbar before clearing axes
        if self._etheta_cbar is not None:
            self._etheta_cbar.remove()
            self._etheta_cbar = None
        self._cv_etheta.clear_ax()
        data = self._cached_hists.get("h2_energy_theta")
        if data is not None:
            try:
                counts, xedges, yedges = data
                pcm = ax.pcolormesh(xedges, yedges, counts.T,
                                    cmap="viridis", shading="flat")
                self._etheta_cbar = self._cv_etheta.fig.colorbar(pcm, ax=ax, pad=0.01)
            except Exception:
                ax.text(0.5, 0.5, "h2_energy_theta\nnot found",
                        transform=ax.transAxes, ha="center", va="center",
                        color=THEME.TEXT_DIM, fontsize=9)
        self._cv_etheta._style_ax("Energy vs θ", "θ (deg)", "Energy (MeV)")
        self._cv_etheta.redraw()

    def _draw_1d_hist(self, hist_key: str, title: str, xlabel: str,
                      color: str, axvline: Optional[float] = None):
        """Generic helper: draw a 1D histogram from the ROOT file."""
        if self._etheta_cbar is not None:
            self._etheta_cbar.remove()
            self._etheta_cbar = None
        ax = self._cv_etheta.ax
        self._cv_etheta.clear_ax()
        data = self._cached_hists.get(hist_key)
        if data is not None:
            try:
                counts, edges = data
                centers = 0.5 * (edges[:-1] + edges[1:])
                nz = counts > 0
                if nz.any():
                    ax.bar(centers, counts, width=(edges[1] - edges[0]),
                           color=color, alpha=0.85, linewidth=0)
                    lo = centers[nz][0] - (edges[1] - edges[0])
                    hi = centers[nz][-1] + (edges[1] - edges[0])
                    ax.set_xlim(lo, hi)
                if axvline is not None:
                    ax.axvline(axvline, color=THEME.WARN,
                               linewidth=1.2, linestyle="--")
            except Exception:
                ax.text(0.5, 0.5, f"{hist_key}\nnot found",
                        transform=ax.transAxes, ha="center", va="center",
                        color=THEME.TEXT_DIM, fontsize=9)
        self._cv_etheta._style_ax(title, xlabel, "Counts")
        self._cv_etheta.redraw()

    def _draw_one_cluster_energy(self):
        self._draw_1d_hist("one_cluster_energy",
                           "One Cluster Energy", "Energy (MeV)",
                           color=THEME.ACCENT)

    def _draw_resolution(self):
        self._draw_1d_hist("recon_sigma",
                           "Resolution (σ)", "σ (MeV)",
                           color=THEME.SUCCESS)


# =============================================================================
#  Background worker
# =============================================================================

class _MetricsWorker(QThread):
    """Run compute_metrics() and preload global hists in a background thread."""
    finished = pyqtSignal(object)   # emits the IterData when done

    _GLOBAL_KEYS = ("ratio_all", "h2_energy_theta",
                    "one_cluster_energy", "recon_sigma")

    def __init__(self, data: "IterData", parent=None):
        super().__init__(parent)
        self._data = data

    def run(self) -> None:
        compute_metrics(self._data)          # fills _hist_cache + metrics
        # also preload global histograms so the main thread has zero ROOT I/O
        if (HAS_UPROOT
                and self._data.root_path is not None
                and not self._data._global_hists):
            try:
                with uproot.open(str(self._data.root_path)) as f:
                    avail = {k.split(";")[0] for k in f.keys()}
                    for key in self._GLOBAL_KEYS:
                        if key in avail:
                            try:
                                self._data._global_hists[key] = f[key].to_numpy()
                            except Exception:
                                pass
            except Exception:
                pass
        self.finished.emit(self._data)


# =============================================================================
#  W-module layer helper
# =============================================================================

def _build_w_layer_map(modules_json_path: Optional[Path]) -> Dict[str, int]:
    """Return {module_name: layer_number} for all PbWO4 (W) modules.

    Layer 1 = outermost ring (adjacent to the PbGlass border), increasing
    inward.  Computed as min(row-1, max_row-row, col-1, max_col-col) + 1
    using the 1-based row/col grid indices stored in hycal_map.json.
    """
    if modules_json_path is None or not modules_json_path.is_file():
        return {}
    try:
        with open(modules_json_path) as f:
            data = json.load(f)
    except Exception:
        return {}
    entries = [
        (e["n"], e["geo"]["row"], e["geo"]["col"])
        for e in data
        if e.get("t") == "PbWO4"
        and "geo" in e
        and "row" in e["geo"]
        and "col" in e["geo"]
    ]
    if not entries:
        return {}
    max_row = max(r for _, r, _ in entries)
    max_col = max(c for _, _, c in entries)
    return {
        name: min(row - 1, max_row - row, col - 1, max_col - col) + 1
        for name, row, col in entries
    }


# =============================================================================
#  Main window
# =============================================================================

class EpCalibViewerWindow(QMainWindow):

    _MODES = [
        CalibMapWidget.MODE_DELTA_E,
        CalibMapWidget.MODE_STATS,
        CalibMapWidget.MODE_SIGMA,
        CalibMapWidget.MODE_CHI2,
        CalibMapWidget.MODE_FACTOR,
    ]

    def __init__(self, initial_dir: Optional[Path] = None):
        super().__init__()
        self._calib_dir: Optional[Path] = None
        self._scan_data: Dict[str, Dict[int, IterData]] = {}
        self._modules:   List[Module] = []
        self._cur_data:  Optional[IterData] = None
        self._map_mode = CalibMapWidget.MODE_DELTA_E
        self._worker:    Optional[_MetricsWorker] = None
        self._abandoned_workers: List[_MetricsWorker] = []
        self._w_layer_map: Dict[str, int] = _build_w_layer_map(MODULES_JSON)

        self._build_ui()
        self.setWindowTitle("epCalib Viewer")
        self.resize(1380, 900)
        self._apply_stylesheet()

        if MODULES_JSON:
            self._modules = load_modules(MODULES_JSON)
            self._map.set_modules(self._modules)

        if initial_dir:
            self._load_calib_dir(initial_dir)

    # ── stylesheet ────────────────────────────────────────────────────────────

    def _apply_stylesheet(self):
        t = THEME
        self.setStyleSheet(f"""
            QMainWindow, QWidget {{ background: {t.BG}; color: {t.TEXT}; }}
            QComboBox {{ background: {t.PANEL}; border: 1px solid {t.BORDER};
                         color: {t.TEXT}; padding: 2px 6px; }}
            QComboBox QAbstractItemView {{ background: {t.PANEL}; color: {t.TEXT}; }}
            QPushButton {{ background: {t.BUTTON}; border: 1px solid {t.BORDER};
                           color: {t.TEXT}; padding: 4px 10px; border-radius: 4px; }}
            QPushButton:hover   {{ background: {t.BUTTON_HOVER}; }}
            QPushButton:checked {{ background: {t.ACCENT_STRONG};
                                   border: 1px solid {t.ACCENT_BORDER}; }}
            QGroupBox {{ color: {t.TEXT_DIM}; font-size: 12px;
                         border: 1px solid {t.BORDER}; margin-top: 6px; }}
            QGroupBox::title {{ subcontrol-origin: margin; left: 8px; }}
            QLineEdit {{ background: {t.PANEL}; border: 1px solid {t.BORDER};
                         color: {t.TEXT}; padding: 2px 4px; border-radius: 3px; }}
            QLineEdit:focus {{ border: 1px solid {t.ACCENT_BORDER}; }}
            QCheckBox {{ color: {t.TEXT}; spacing: 5px; }}
            QCheckBox::indicator {{ width: 14px; height: 14px;
                                    background: {t.PANEL}; border: 1px solid {t.BORDER};
                                    border-radius: 3px; }}
            QCheckBox::indicator:checked {{ background: {t.ACCENT_STRONG};
                                            border: 1px solid {t.ACCENT_BORDER}; }}
            QSplitter::handle {{ background: {t.PANEL}; }}
            QStatusBar {{ color: {t.TEXT_DIM}; font-size: 11px; }}
        """)

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(5)

        # ── toolbar ───────────────────────────────────────────────────────────
        top = QHBoxLayout()
        top.setSpacing(8)
        title = QLabel("epCalib Viewer")
        title.setFont(QFont("Monospace", 15, QFont.Weight.Bold))
        title.setStyleSheet(f"color: {THEME.ACCENT};")
        top.addWidget(title)
        self._dir_btn = QPushButton("Open Physics_calib dir…")
        self._dir_btn.clicked.connect(self._browse_dir)
        top.addWidget(self._dir_btn)
        self._refresh_btn = QPushButton("⟳ Refresh")
        self._refresh_btn.setToolTip("Rescan the current directory for new/updated files")
        self._refresh_btn.clicked.connect(self._refresh_dir)
        self._refresh_btn.setEnabled(False)
        top.addWidget(self._refresh_btn)
        self._dir_lbl = QLabel("(no directory)")
        self._dir_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 12px;")
        top.addWidget(self._dir_lbl)
        top.addStretch()
        top.addWidget(QLabel("Run:"))
        self._run_cb = QComboBox()
        self._run_cb.setMinimumWidth(100)
        self._run_cb.currentTextChanged.connect(self._on_run_changed)
        top.addWidget(self._run_cb)
        top.addWidget(QLabel("Iteration:"))
        self._iter_cb = QComboBox()
        self._iter_cb.setMinimumWidth(70)
        self._iter_cb.currentTextChanged.connect(self._on_iter_changed)
        top.addWidget(self._iter_cb)
        root.addLayout(top)

        # ── map mode buttons ──────────────────────────────────────────────────
        mode_bar = QHBoxLayout()
        mode_bar.addWidget(QLabel("Map:"))
        self._mode_grp = QButtonGroup(self)
        self._mode_grp.setExclusive(True)
        for mode in self._MODES:
            btn = QPushButton(mode)
            btn.setCheckable(True)
            if mode == self._map_mode:
                btn.setChecked(True)
            btn.clicked.connect(lambda _, m=mode: self._set_map_mode(m))
            self._mode_grp.addButton(btn)
            mode_bar.addWidget(btn)
        mode_bar.addStretch()
        self._bad_lbl = QLabel("")
        self._bad_lbl.setStyleSheet(f"color: {THEME.DANGER}; font-size: 11px;")
        mode_bar.addWidget(self._bad_lbl)
        root.addLayout(mode_bar)

        # ── Z-axis range / log controls ───────────────────────────────────────
        # Reusable widget from hycal_geoview.  Click "Auto" for a one-shot
        # range fit; double-click to keep auto-fitting on every refresh.
        # Mode switches still trigger an explicit one-shot fit via auto_fit().
        self._map = CalibMapWidget()
        self._map.moduleClicked.connect(self._on_module_clicked)
        self._map.moduleHovered.connect(self._on_module_hovered)
        self._map.multiSelectionChanged.connect(self._on_multi_selection_changed)

        z_bar = QHBoxLayout()
        z_bar.setSpacing(6)
        self._range_ctrl = ColorRangeControl(
            self._map,
            auto_fit="minmax",
            include_log=True,
        )
        z_bar.addWidget(self._range_ctrl)
        z_bar.addStretch()
        root.addLayout(z_bar)

        # ── main splitter ─────────────────────────────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(4)
        root.addWidget(splitter, stretch=1)

        # left: HyCal map
        map_box = QGroupBox("HyCal Map  (scroll to zoom · drag to pan · click module)")
        ml = QVBoxLayout(map_box)
        ml.setContentsMargins(2, 14, 2, 2)
        ml.addWidget(self._map)

        # ── batch multi-select bar ────────────────────────────────────────
        batch_bar = QHBoxLayout()
        batch_bar.setSpacing(6)
        self._multiselect_btn = QPushButton("◯ Multi-Select")
        self._multiselect_btn.setCheckable(True)
        self._multiselect_btn.setToolTip(
            "Toggle multi-select mode\n"
            "Click modules to add/remove from selection (blue circle)\n"
            "Then set a factor for all selected modules at once,\n"
            "or restore old factors (from .dat) for all selected modules")
        self._multiselect_btn.toggled.connect(self._map.set_multi_select_mode)
        batch_bar.addWidget(self._multiselect_btn)
        self._sel_count_lbl = QLabel("0 selected")
        self._sel_count_lbl.setStyleSheet(
            f"color: {THEME.TEXT_DIM}; font-size: 11px; min-width: 72px;")
        batch_bar.addWidget(self._sel_count_lbl)
        batch_bar.addSpacing(8)
        batch_bar.addWidget(QLabel("Factor:"))
        self._batch_factor_edit = QLineEdit()
        self._batch_factor_edit.setFixedWidth(80)
        self._batch_factor_edit.setPlaceholderText("value")
        batch_bar.addWidget(self._batch_factor_edit)
        for _preset in (0.122, 0.15):
            _btn = QPushButton(f"{_preset}")
            _btn.setFixedWidth(54)
            _btn.setToolTip(f"Fill factor field with {_preset}")
            _btn.clicked.connect(
                lambda _, v=_preset: self._batch_factor_edit.setText(str(v)))
            batch_bar.addWidget(_btn)
        self._batch_apply_btn = QPushButton("Apply to All Selected")
        self._batch_apply_btn.setEnabled(False)
        self._batch_apply_btn.setToolTip(
            "Write the given factor to all selected modules and save JSON")
        self._batch_apply_btn.clicked.connect(self._do_batch_set_factor)
        batch_bar.addWidget(self._batch_apply_btn)
        self._batch_restore_btn = QPushButton("\u21a9 Restore Old Factor")
        self._batch_restore_btn.setEnabled(False)
        self._batch_restore_btn.setToolTip(
            "Restore old_factor (pre-iteration value from .dat) for all selected modules\n"
            "and save JSON.  Modules without a .dat entry are skipped.")
        self._batch_restore_btn.clicked.connect(self._do_batch_restore_old_factor)
        batch_bar.addWidget(self._batch_restore_btn)
        self._batch_clear_btn = QPushButton("Clear Selection")
        self._batch_clear_btn.setEnabled(False)
        self._batch_clear_btn.clicked.connect(self._do_batch_clear_selection)
        batch_bar.addWidget(self._batch_clear_btn)
        # ── outer W-layer restore ────────────────────────────────────────
        from PyQt6.QtWidgets import QSpinBox
        sep = QLabel("|")  # visual separator
        sep.setStyleSheet(f"color: {THEME.BORDER}; margin: 0 4px;")
        batch_bar.addWidget(sep)
        batch_bar.addWidget(QLabel("Outer W layers:"))
        self._outer_layer_spin = QSpinBox()
        self._outer_layer_spin.setRange(1, 3)
        self._outer_layer_spin.setValue(1)
        self._outer_layer_spin.setFixedWidth(44)
        self._outer_layer_spin.setToolTip(
            "Number of outermost PbWO4 rings to include\n"
            "(1 = single border ring, 2 = two outer rings, 3 = three)")
        batch_bar.addWidget(self._outer_layer_spin)
        self._outer_restore_btn = QPushButton("\u21a9 Restore Outer W Old Factor")
        self._outer_restore_btn.setToolTip(
            "Restore old_factor (from .dat) for the outermost N rings of W (PbWO4) modules\n"
            "and save JSON.  Affected modules are marked with a red cross.")
        self._outer_restore_btn.clicked.connect(self._do_restore_outer_w_layers)
        batch_bar.addWidget(self._outer_restore_btn)
        batch_bar.addStretch()
        ml.addLayout(batch_bar)

        splitter.addWidget(map_box)

        # right panel
        right = QSplitter(Qt.Orientation.Vertical)
        right.setHandleWidth(4)
        splitter.addWidget(right)

        detail_box = QGroupBox("Module Detail")
        dl = QVBoxLayout(detail_box)
        dl.setContentsMargins(2, 14, 2, 2)
        self._detail = ModuleDetailPanel()
        self._detail.refitApplied.connect(self._on_refit_applied)
        dl.addWidget(self._detail)
        right.addWidget(detail_box)

        stats_box = QGroupBox("Global Statistics")
        sl = QVBoxLayout(stats_box)
        sl.setContentsMargins(2, 14, 2, 2)
        self._stats = GlobalStatsPanel()
        sl.addWidget(self._stats)
        right.addWidget(stats_box)

        splitter.setSizes([580, 670])
        right.setSizes([490, 290])

        # status bar
        self._hover_lbl = QLabel("")
        self.statusBar().addWidget(self._hover_lbl)
        if not HAS_UPROOT:
            self.statusBar().showMessage(
                "uproot not installed — histograms unavailable.  "
                "Run: pip install uproot", 0)
        elif not HAS_SCIPY:
            self.statusBar().showMessage(
                "scipy not installed — Gaussian fitting disabled.  "
                "Run: pip install scipy", 0)

    # ── directory loading ─────────────────────────────────────────────────────

    def _browse_dir(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select Physics_calib directory",
            str(self._calib_dir or Path.cwd()))
        if path:
            self._load_calib_dir(Path(path))

    def _refresh_dir(self):
        if self._calib_dir:
            cur_run = self._run_cb.currentText()
            cur_iter = self._iter_cb.currentText()
            self._load_calib_dir(self._calib_dir)
            # try to restore previous run/iter selection
            if cur_run in self._scan_data:
                self._run_cb.setCurrentText(cur_run)
                self._on_run_changed(cur_run)
                if cur_iter:
                    self._iter_cb.setCurrentText(cur_iter)
                    self._on_iter_changed(cur_iter)

    def _load_calib_dir(self, path: Path):
        self._calib_dir = path
        self._dir_lbl.setText(str(path))
        self._scan_data = scan_calib_dir(path)
        self._refresh_btn.setEnabled(True)

        self._run_cb.blockSignals(True)
        self._run_cb.clear()
        for run in sorted(self._scan_data.keys()):
            self._run_cb.addItem(run)
        self._run_cb.blockSignals(False)

        if self._scan_data:
            self._run_cb.setCurrentIndex(0)
            self._on_run_changed(self._run_cb.currentText())

    # ── combo callbacks ───────────────────────────────────────────────────────

    def _on_run_changed(self, run: str):
        if run not in self._scan_data:
            return
        self._iter_cb.blockSignals(True)
        self._iter_cb.clear()
        for it in sorted(self._scan_data[run].keys()):
            self._iter_cb.addItem(str(it))
        self._iter_cb.blockSignals(False)
        if self._iter_cb.count():
            self._iter_cb.setCurrentIndex(self._iter_cb.count() - 1)
            self._on_iter_changed(self._iter_cb.currentText())

    def _on_iter_changed(self, it_str: str):
        if not it_str:
            return
        run = self._run_cb.currentText()
        it  = int(it_str)
        data = self._scan_data.get(run, {}).get(it)
        if data is None:
            return
        self._cur_data = data
        self._detail.set_root_path(data.root_path)
        self._detail.set_iter_data(data)
        # if metrics + global hists are already cached, update UI immediately
        if data.metrics and data._global_hists:
            self._stats.set_preloaded(data.root_path, data._global_hists)
            self._refresh_map()
            self._map.set_marked_modules(data.modified_modules)
            return
        if data.metrics and not data._global_hists:
            # metrics done but global hists not yet - rare; use blocking path
            self._stats.set_root_path(data.root_path)
            self._refresh_map()
            self._map.set_marked_modules(data.modified_modules)
            return
        # Disconnect previous worker (if still running) without blocking the UI.
        # Keep a reference in _abandoned_workers so the QThread object stays alive
        # until it finishes naturally.
        if self._worker is not None and self._worker.isRunning():
            try:
                self._worker.finished.disconnect(self._on_metrics_ready)
            except TypeError:
                pass
            self._abandoned_workers.append(self._worker)
        # Prune any already-finished abandoned workers
        self._abandoned_workers = [w for w in self._abandoned_workers if w.isRunning()]
        self.statusBar().showMessage(
            f"Loading run {run} iter {it} …")
        self._worker = _MetricsWorker(data, self)
        self._worker.finished.connect(self._on_metrics_ready)
        self._worker.start()

    def _on_metrics_ready(self, data: IterData) -> None:
        """Called from background worker when compute_metrics + global hist load finish."""
        self.statusBar().clearMessage()
        if data is self._cur_data:
            self._stats.set_preloaded(data.root_path, data._global_hists)
            self._refresh_map()
            self._map.set_marked_modules(data.modified_modules)
        self._worker = None

    # ── map mode ──────────────────────────────────────────────────────────────

    def _set_map_mode(self, mode: str):
        self._map_mode = mode
        # Mode switch always re-fits — different modes have wildly different
        # value scales, so keeping the previous range would be useless.
        if self._cur_data is not None:
            self._refresh_map()
            self._range_ctrl.auto_fit()

    def _refresh_map(self):
        data = self._cur_data
        if data is None:
            return

        values: Dict[str, float] = {}
        mode = self._map_mode

        if mode == CalibMapWidget.MODE_DELTA_E:
            for name, mm in data.metrics.items():
                if mm.base_energy > 0 and mm.peak > 0:
                    values[name] = abs(mm.peak - mm.base_energy) / mm.base_energy
            label = "ΔE/E = |peak − expected| / expected"

        elif mode == CalibMapWidget.MODE_STATS:
            for name, mm in data.metrics.items():
                values[name] = mm.stats
            label = "Event count per module"

        elif mode == CalibMapWidget.MODE_SIGMA:
            for name, mm in data.metrics.items():
                if mm.sigma > 0:
                    values[name] = mm.sigma
            label = "Gaussian σ (MeV)"

        elif mode == CalibMapWidget.MODE_CHI2:
            for name, mm in data.metrics.items():
                if mm.chi2 > 0:
                    values[name] = mm.chi2
            label = "χ²/ndf  (Gaussian fit)"

        elif mode == CalibMapWidget.MODE_FACTOR:
            for name, mm in data.metrics.items():
                values[name] = mm.factor
            label = "Calibration factor"

        self._map.set_map_label(label)
        self._map.set_values(values)
        # Range control re-fits when pinned; otherwise the user's manual
        # values stay put.
        self._range_ctrl.notify_values_changed(values)

        # count "bad" modules: chi2 > 2 or |ΔE/E| > 2%
        n_bad = sum(
            1 for mm in data.metrics.values()
            if mm.chi2 > 2.0 or (
                mm.base_energy > 0 and mm.peak > 0
                and abs(mm.peak - mm.base_energy) / mm.base_energy > 0.02
            )
        )
        total = len(data.metrics)
        self._bad_lbl.setText(
            f"  Flagged (χ²>2 or |ΔE/E|>2%): {n_bad}/{total}" if total else "")

    # ── module interaction ────────────────────────────────────────────────────

    def _on_module_clicked(self, name: str):
        if not name or self._cur_data is None:
            return
        mm = self._cur_data.metrics.get(name)
        self._detail.show_module(name, mm)

    def _on_module_hovered(self, name: str):
        if not name or self._cur_data is None:
            self._hover_lbl.setText("")
            return
        mm = self._cur_data.metrics.get(name)
        if mm:
            de = ""
            if mm.base_energy > 0 and mm.peak > 0:
                de = f"  ΔE/E={(mm.peak - mm.base_energy)/mm.base_energy*100:+.2f}%"
            self._hover_lbl.setText(
                f"{name}  |  factor={mm.factor:.4f}  "
                f"events={mm.stats:.0f}  "
                f"peak={mm.peak:.1f} MeV  σ={mm.sigma:.1f} MeV  "
                f"χ²/ndf={mm.chi2:.3f}{de}")
        else:
            self._hover_lbl.setText(name)

    def _on_refit_applied(self, _name: str) -> None:
        """Called after a manual refit is saved; refresh map colors and marks."""
        self._refresh_map()
        if self._cur_data is not None:
            self._map.set_marked_modules(self._cur_data.modified_modules)

    def _on_multi_selection_changed(self, selected: set) -> None:
        n = len(selected)
        self._sel_count_lbl.setText(f"{n} selected")
        self._batch_apply_btn.setEnabled(n > 0)
        self._batch_restore_btn.setEnabled(n > 0)
        self._batch_clear_btn.setEnabled(n > 0)

    def _do_batch_clear_selection(self) -> None:
        self._map.clear_selection()

    def _do_batch_set_factor(self) -> None:
        """Apply the batch factor to all multi-selected modules and save JSON."""
        txt = self._batch_factor_edit.text().strip()
        try:
            value = float(txt)
        except ValueError:
            self.statusBar().showMessage(
                "Batch apply: invalid factor value — enter a number", 4000)
            return
        if value <= 0:
            self.statusBar().showMessage(
                "Batch apply: factor must be positive", 4000)
            return
        data = self._cur_data
        if data is None:
            return
        selected = self._map.get_selected_modules()
        if not selected:
            return

        updated: List[str] = []
        for name in selected:
            mm = data.metrics.get(name)
            if mm is None:
                continue
            mm.factor = value
            entry = data.factors.get(name)
            if entry is not None:
                entry["factor"] = value
            data.modified_modules.add(name)
            updated.append(name)

        jpath = data.json_path
        if jpath is not None:
            try:
                entries = list(data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                note = f"saved {jpath.name}"
            except Exception as exc:
                self.statusBar().showMessage(
                    f"Batch apply: save failed: {exc}", 6000)
                return
        else:
            note = "JSON path not found — not saved"

        self.statusBar().showMessage(
            f"Set factor={value:.5f} for {len(updated)} module(s).  {note}", 6000)

        # turn off multi-select mode, clear selection, refresh display
        self._multiselect_btn.setChecked(False)
        self._map.clear_selection()
        self._refresh_map()
        self._map.set_marked_modules(data.modified_modules)

    def _do_batch_restore_old_factor(self) -> None:
        """Restore old_factor (from .dat) for all multi-selected modules and save JSON."""
        data = self._cur_data
        if data is None:
            return
        selected = self._map.get_selected_modules()
        if not selected:
            return

        updated: List[str] = []
        skipped: List[str] = []
        for name in selected:
            mm = data.metrics.get(name)
            if mm is None or mm.old_factor <= 0.0:
                skipped.append(name)
                continue
            mm.factor = mm.old_factor
            entry = data.factors.get(name)
            if entry is not None:
                entry["factor"] = mm.old_factor
            data.modified_modules.add(name)
            updated.append(name)

        if not updated:
            self.statusBar().showMessage(
                "Batch restore: no modules had old_factor available (.dat not loaded?)", 6000)
            return

        jpath = data.json_path
        if jpath is not None:
            try:
                entries = list(data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                note = f"saved {jpath.name}"
            except Exception as exc:
                self.statusBar().showMessage(
                    f"Batch restore: save failed: {exc}", 6000)
                return
        else:
            note = "JSON path not found — not saved"

        msg = f"Restored old_factor for {len(updated)} module(s).  {note}"
        if skipped:
            msg += f"  ({len(skipped)} skipped — no old_factor)"
        self.statusBar().showMessage(msg, 6000)

        # turn off multi-select mode, clear selection, refresh display
        self._multiselect_btn.setChecked(False)
        self._map.clear_selection()
        self._refresh_map()
        self._map.set_marked_modules(data.modified_modules)

    def _do_restore_outer_w_layers(self) -> None:
        """Restore old_factor for the outermost N rings of PbWO4 (W) modules."""
        data = self._cur_data
        if data is None:
            self.statusBar().showMessage("No calibration data loaded", 4000)
            return
        if not self._w_layer_map:
            self.statusBar().showMessage(
                "W layer map not available — check hycal_map.json path", 4000)
            return

        n_layers = self._outer_layer_spin.value()
        target_names = {
            name for name, layer in self._w_layer_map.items() if layer <= n_layers
        }

        updated: List[str] = []
        skipped: List[str] = []
        for name in target_names:
            mm = data.metrics.get(name)
            if mm is None or mm.old_factor <= 0.0:
                skipped.append(name)
                continue
            mm.factor = mm.old_factor
            entry = data.factors.get(name)
            if entry is not None:
                entry["factor"] = mm.old_factor
            data.modified_modules.add(name)
            updated.append(name)

        if not updated:
            self.statusBar().showMessage(
                f"Outer W restore: no modules had old_factor available "
                f"(.dat not loaded?), checked {len(target_names)} modules", 6000)
            return

        jpath = data.json_path
        if jpath is not None:
            try:
                entries = list(data.factors.values())
                with open(jpath, "w") as fj:
                    json.dump(entries, fj, indent=2)
                note = f"saved {jpath.name}"
            except Exception as exc:
                self.statusBar().showMessage(
                    f"Outer W restore: save failed: {exc}", 6000)
                return
        else:
            note = "JSON path not found — not saved"

        msg = (
            f"Restored old_factor for {len(updated)} outer W module(s) "
            f"(≤{n_layers} layer(s)).  {note}"
        )
        if skipped:
            msg += f"  ({len(skipped)} skipped — no old_factor)"
        self.statusBar().showMessage(msg, 8000)
        self._refresh_map()
        self._map.set_marked_modules(data.modified_modules)


# =============================================================================
#  Entry point
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description="epCalib output GUI viewer")
    ap.add_argument("calib_dir", nargs="?", type=Path,
                    help="Path to Physics_calib directory")
    ap.add_argument("--theme", choices=available_themes(), default="dark")
    args = ap.parse_args()

    set_theme(args.theme)

    app = QApplication(sys.argv)
    app.setApplicationName("epCalib Viewer")
    win = EpCalibViewerWindow(initial_dir=args.calib_dir)
    apply_theme_palette(win)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
