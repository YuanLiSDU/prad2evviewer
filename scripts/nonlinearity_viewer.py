#!/usr/bin/env python3
"""
nonlinearity_viewer.py — HyCal PWO4 non-linearity analysis viewer
==================================================================

Reads the ROOT output of prad2ana_nonlinearity:
  energy_3p5GeV/h_energy_{name}_3p5   per-module TH1F at 3.5 GeV
  energy_0p7GeV/h_energy_{name}_0p7   per-module TH1F at 0.7 GeV

Left:  HyCal map (W modules only), coloured by non-linearity parameter nl.
Right: three stacked plots for the selected module:
         • 3.5 GeV energy histogram
         • 0.7 GeV energy histogram
         • Linearity plot (E_measured vs E_expected)
       Plus controls to re-fit individual peaks and re-run the nl fit.

All fit algorithms match nonlinearity.cpp exactly.

Usage:
  python scripts/nonlinearity_viewer.py [root_file] [--theme dark|light]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QApplication, QButtonGroup, QComboBox, QFileDialog, QGridLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton,
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

# ── path setup for shared scripts ─────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from hycal_geoview import (
    HyCalMapWidget, Module, load_modules,
    apply_theme_palette, set_theme, available_themes, THEME,
    ColorRangeControl,
)

# ── database paths ─────────────────────────────────────────────────────────
_DB_CANDIDATES = [
    _SCRIPT_DIR.parent / "database" / "hycal_map.json",
    _SCRIPT_DIR.parent / "build" / "database" / "hycal_map.json",
]
MODULES_JSON = next((p for p in _DB_CANDIDATES if p.is_file()), None)


# =============================================================================
# Physics constants & helpers  (match PhysicsTools.cpp)
# =============================================================================

M_PROTON   = 938.272   # MeV
M_ELECTRON = 0.511     # MeV
Z_HYCAL    = 6270.0    # mm
RESOLUTION = 0.035     # 3.5% / sqrt(E/1000)


def _energy_loss(theta_deg: float) -> float:
    """Material energy loss in MeV — identical to PhysicsTools::EnergyLoss."""
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sec = 1.0 / cos_t if cos_t > 0.01 else 100.0
    eloss  = 0.500 * 1.6 * sec   # Al window
    eloss += 0.120 * 1.6 * sec   # GEM Al foils (2 GEMs)
    eloss += 0.100 * 2.0 * sec   # GEM foils (2 GEMs)
    eloss += 0.480 * 1.8 * sec   # kapton cover
    return eloss


def expected_energy(theta_deg: float, Ebeam: float, kind: str) -> float:
    """Expected scattered energy in MeV — identical to PhysicsTools::ExpectedEnergy."""
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    if kind == "ep":
        E = Ebeam * M_PROTON / (M_PROTON + Ebeam * (1.0 - cos_t))
    elif kind == "ee":
        gamma = Ebeam / M_ELECTRON
        num = (gamma + 1.0) + (gamma - 1.0) * cos_t * cos_t
        den = (gamma + 1.0) - (gamma - 1.0) * cos_t * cos_t
        if den <= 0.0:
            return 0.0
        E = M_ELECTRON * num / den
    else:
        return 0.0
    return max(0.0, E - _energy_loss(theta_deg))


def estimate_sigma(E: float) -> float:
    """Energy resolution estimate: σ = E × 3.5% / √(E/1000)."""
    return E * RESOLUTION / math.sqrt(max(E, 1.0) / 1000.0)


# =============================================================================
# Gaussian fit helpers  (match nonlinearity.cpp fitPeakAndDraw)
# =============================================================================

def _gauss(x: np.ndarray, amp: float, mu: float, sig: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def _gauss_fit_window(
    centers: np.ndarray,
    counts: np.ndarray,
    mu0: float,
    half_width: float,
) -> Optional[Tuple[float, float, float]]:
    """Fit Gaussian in [mu0 ± half_width].  Returns (mean, sigma, amp) or None."""
    lo, hi = mu0 - half_width, mu0 + half_width
    mask = (centers >= lo) & (centers <= hi) & (counts > 0)
    if mask.sum() < 4:
        return None
    xf, yf = centers[mask], counts[mask].astype(float)
    amp0 = float(yf.max())
    bw   = float(centers[1] - centers[0]) if len(centers) > 1 else 1.0
    sig0 = max(half_width / 3.0, bw)
    try:
        popt, _ = curve_fit(
            _gauss, xf, yf,
            p0=[amp0, mu0, sig0],
            bounds=([0.0, lo, 1e-6], [np.inf, hi, half_width * 2.0]),
            maxfev=2000,
        )
        return float(popt[1]), abs(float(popt[2])), float(popt[0])
    except Exception:
        return None


def fit_peak_full(
    counts: np.ndarray,
    edges: np.ndarray,
    Eexp: float,
    sigma_exp: float,
) -> Tuple[float, float, float]:
    """
    3-step Gaussian fit matching nonlinearity.cpp fitPeakAndDraw.

    Step 2: weighted mean in search window [Eexp ± 6σ_exp]
    Step 3: first Gaussian fit in [mean ± 2σ(mean)]
    Step 4: final Gaussian fit in [mean ± 1σ(mean)]

    Returns (mean_MeV, sigma_MeV, amplitude).  mean=0 on failure.
    """
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Step 2: weighted mean in search window
    mask = (centers >= Eexp - 6.0 * sigma_exp) & (centers <= Eexp + 6.0 * sigma_exp)
    w = counts[mask]
    x = centers[mask]
    if w.sum() <= 0.0:
        return 0.0, 0.0, 0.0

    mean = float((x * w).sum() / w.sum())

    if not HAS_SCIPY:
        return mean, estimate_sigma(mean), float(w.max())

    best: Optional[Tuple[float, float, float]] = None

    # Step 3: first fit in [mean ± 2σ]
    r = _gauss_fit_window(centers, counts, mean, 2.0 * estimate_sigma(mean))
    if r is not None:
        mean = r[0]
        best = r

    # Step 4: final fit in [mean ± 1σ]
    r = _gauss_fit_window(centers, counts, mean, estimate_sigma(mean))
    if r is not None:
        return r
    if best is not None:
        return best
    return mean, estimate_sigma(mean), float(w.max())


# =============================================================================
# Non-linearity model & fit  (match nonlinearity.cpp)
# =============================================================================
# Model: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000  [1st order]
#        E_rec/E_exp = 1 + nl1*t + nl2*t^2           [2nd order], t=(E_rec-E_base)/1000

def nl_ratio_1st(E_rec: np.ndarray, nl1: float, E_base: float) -> np.ndarray:
    """1st order: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000"""
    return 1.0 + nl1 * (E_rec - E_base) / 1000.0


def nl_ratio_2nd(E_rec: np.ndarray, nl1: float, nl2: float, E_base: float) -> np.ndarray:
    """2nd order: E_rec/E_exp = 1 + nl1*t + nl2*t^2, t=(E_rec-E_base)/1000"""
    t = (E_rec - E_base) / 1000.0
    return 1.0 + nl1 * t + nl2 * t * t


def fit_nonlinearity(
    E_exp: np.ndarray,
    E_rec: np.ndarray,
    E_base: float,
) -> Tuple[float, float, float, int]:
    """Fit 1st-order nl model to E_rec/E_exp vs E_rec.
    Returns (nl1, nl1_err, chi2, ndf)."""
    if not HAS_SCIPY or len(E_exp) < 2:
        return 0.0, 0.0, 0.0, 0
    ratio = E_rec / np.where(E_exp > 0, E_exp, 1.0)
    try:
        popt, pcov = curve_fit(
            lambda x, nl1: nl_ratio_1st(x, nl1, E_base),
            E_rec, ratio,
            p0=[0.01],
            maxfev=5000,
        )
        nl1     = float(popt[0])
        nl1_err = float(np.sqrt(max(pcov[0, 0], 0.0)))
        resid   = ratio - nl_ratio_1st(E_rec, nl1, E_base)
        chi2    = float(np.sum(resid ** 2))
        ndf     = len(E_rec) - 1
        return nl1, nl1_err, chi2, ndf
    except Exception:
        return 0.0, 0.0, 0.0, 0


def fit_nonlinearity_2nd(
    E_exp: np.ndarray,
    E_rec: np.ndarray,
    E_base: float,
    nl1_init: float = 0.01,
) -> Tuple[float, float, float, float, float, int]:
    """Fit 2nd-order nl model to E_rec/E_exp vs E_rec.
    Returns (nl1, nl1_err, nl2, nl2_err, chi2, ndf)."""
    if not HAS_SCIPY or len(E_exp) < 3:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0
    ratio = E_rec / np.where(E_exp > 0, E_exp, 1.0)
    try:
        popt, pcov = curve_fit(
            lambda x, nl1, nl2: nl_ratio_2nd(x, nl1, nl2, E_base),
            E_rec, ratio,
            p0=[nl1_init, 0.0],
            maxfev=5000,
        )
        nl1     = float(popt[0])
        nl2     = float(popt[1])
        perr    = np.sqrt(np.maximum(np.diag(pcov), 0.0))
        nl1_err = float(perr[0])
        nl2_err = float(perr[1])
        resid   = ratio - nl_ratio_2nd(E_rec, nl1, nl2, E_base)
        chi2    = float(np.sum(resid ** 2))
        ndf     = len(E_rec) - 2
        return nl1, nl1_err, nl2, nl2_err, chi2, ndf
    except Exception:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0


# =============================================================================
# Per-module data  (peak keys: ep_3p5, ee_3p5, ep_0p7, ee_0p7)
# =============================================================================

_PEAK_KEYS   = ["ep_3p5", "ee_3p5", "ep_0p7", "ee_0p7"]
_PEAK_LABELS = {
    "ep_3p5": "e-p  3.5 GeV",
    "ee_3p5": "e-e  3.5 GeV",
    "ep_0p7": "e-p  0.7 GeV",
    "ee_0p7": "e-e  0.7 GeV",
}
_PEAK_EBEAM = {"ep_3p5": 3485.0, "ee_3p5": 3485.0, "ep_0p7": 729.0, "ee_0p7": 729.0}
_PEAK_KIND  = {"ep_3p5": "ep",   "ee_3p5": "ee",   "ep_0p7": "ep",  "ee_0p7": "ee"}
_PEAK_COLOR = {"ep_3p5": "tomato",  "ee_3p5": "cornflowerblue",
               "ep_0p7": "tomato",  "ee_0p7": "cornflowerblue"}
_PEAK_TAG   = {"ep_3p5": "3p5", "ee_3p5": "3p5", "ep_0p7": "0p7", "ee_0p7": "0p7"}


class ModuleData:
    """All non-linearity data for one PbWO4 module."""

    def __init__(self, name: str, x: float, y: float, row: int = -1, col: int = -1):
        self.name  = name
        self.x     = float(x)
        self.y     = float(y)
        self.theta = math.degrees(math.atan2(math.hypot(x, y), Z_HYCAL))

        # expected energies and sigmas
        self.e_exp:     Dict[str, float] = {}
        self.sigma_exp: Dict[str, float] = {}
        for key in _PEAK_KEYS:
            E = expected_energy(self.theta, _PEAK_EBEAM[key], _PEAK_KIND[key])
            self.e_exp[key]     = E
            self.sigma_exp[key] = estimate_sigma(E) if E > 0.0 else 1.0

        # fitted peak results (0 = not found / invalid)
        self.peaks:      Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.peak_sigma: Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.peak_amp:   Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}

        # calibration scale: e_exp["ep_3p5"] / peaks["ep_3p5"]
        self.scale: float = 1.0

        # non-linearity fit result (1st order)
        self.nl:      float = 0.0
        self.nl_err:  float = 0.0
        self.chi2_nl: float = 0.0
        self.ndf_nl:  int   = 0
        # non-linearity fit result (2nd order)
        self.nl2_1:      float = 0.0
        self.nl2_1_err:  float = 0.0
        self.nl2_2:      float = 0.0
        self.nl2_2_err:  float = 0.0
        self.chi2_nl2:   float = 0.0
        self.ndf_nl2:    int   = 0

        # geometry flags (0-based row/col from JSON)
        self.is_outer    = (row in (0, 33) or col in (0, 33)) if row >= 0 else False
        self.is_absorber = (15 <= row <= 18 and 15 <= col <= 18) if row >= 0 else False

        # raw histograms loaded from ROOT (counts, edges)
        self.hist_3p5: Optional[Tuple[np.ndarray, np.ndarray]] = None
        self.hist_0p7: Optional[Tuple[np.ndarray, np.ndarray]] = None

    def scaled(self, key: str) -> float:
        """Return peak value scaled by calibration factor."""
        return self.peaks[key] * self.scale

    def recompute_scale(self) -> None:
        """Recompute calibration scale from ep_3p5 anchor peak."""
        p = self.peaks.get("ep_3p5", 0.0)
        E = self.e_exp.get("ep_3p5", 0.0)
        self.scale = E / p if (p > 0.0 and E > 0.0) else 1.0

    def linearity_xy(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return (E_expected, E_rec_scaled) arrays for all valid peaks."""
        xs, ys = [], []
        for key in _PEAK_KEYS:
            p = self.scaled(key)
            if p > 0.0:
                xs.append(self.e_exp[key])
                ys.append(p)
        return np.array(xs), np.array(ys)

    def refit_nl(self) -> None:
        """Re-run the 1st and 2nd order non-linearity fits from current scaled peaks."""
        E_exp, E_rec = self.linearity_xy()
        E_base = self.e_exp["ep_3p5"]
        if len(E_exp) >= 2 and E_base > 0.0:
            nl, err, chi2, ndf = fit_nonlinearity(E_exp, E_rec, E_base)
            self.nl      = nl if not (self.is_outer or self.is_absorber) else 0.0
            self.nl_err  = err
            self.chi2_nl = chi2
            self.ndf_nl  = ndf
        else:
            self.nl = self.nl_err = self.chi2_nl = 0.0
            self.ndf_nl = 0
        if len(E_exp) >= 3 and E_base > 0.0:
            nl1, nl1_err, nl2, nl2_err, chi2_2, ndf_2 = fit_nonlinearity_2nd(
                E_exp, E_rec, E_base, nl1_init=self.nl)
            self.nl2_1     = nl1 if not (self.is_outer or self.is_absorber) else 0.0
            self.nl2_1_err = nl1_err
            self.nl2_2     = nl2 if not (self.is_outer or self.is_absorber) else 0.0
            self.nl2_2_err = nl2_err
            self.chi2_nl2  = chi2_2
            self.ndf_nl2   = ndf_2
        else:
            self.nl2_1 = self.nl2_1_err = self.nl2_2 = self.nl2_2_err = self.chi2_nl2 = 0.0
            self.ndf_nl2 = 0

        # Mean |corrected ratio - 1| for map colouring
        E_base = self.e_exp["ep_3p5"]
        self.mean_corr_resid_1st: float = 0.0
        self.mean_corr_resid_2nd: float = 0.0
        if len(E_exp) > 0 and E_base > 0.0:
            resids1, resids2 = [], []
            for Ee, Er in zip(E_exp, E_rec):
                if Ee <= 0.0:
                    continue
                d1 = 1.0 + self.nl * (Er - E_base) / 1000.0
                if d1 != 0.0:
                    Ec1 = Er / d1
                    resids1.append(abs(Ec1 / Ee - 1.0))
                t = (Er - E_base) / 1000.0
                d2 = 1.0 + self.nl2_1 * t + self.nl2_2 * t * t
                if d2 != 0.0:
                    Ec2 = Er / d2
                    resids2.append(abs(Ec2 / Ee - 1.0))
            if resids1:
                self.mean_corr_resid_1st = float(np.mean(resids1))
            if resids2:
                self.mean_corr_resid_2nd = float(np.mean(resids2))


# =============================================================================
# Data loading
# =============================================================================

def _load_json_geometry(path: Path) -> Dict[str, Tuple[float, float, int, int]]:
    """Read hycal_map.json and return {name: (x, y, row_0based, col_0based)}
    for PbWO4 modules only."""
    result: Dict[str, Tuple[float, float, int, int]] = {}
    try:
        with open(path) as f:
            data = json.load(f)
        for e in data:
            if e.get("t") != "PbWO4":
                continue
            g = e.get("geo", {})
            name = e.get("n", "")
            x    = float(g.get("x", 0.0))
            y    = float(g.get("y", 0.0))
            row  = int(g.get("row", 0)) - 1   # JSON is 1-based → convert to 0-based
            col  = int(g.get("col", 0)) - 1
            result[name] = (x, y, row, col)
    except Exception as e:
        print(f"Warning: cannot read geometry from {path}: {e}", file=sys.stderr)
    return result


def load_all(root_path: Optional[Path]) -> Dict[str, ModuleData]:
    """Load module geometry, histograms from ROOT, and fit all peaks."""
    if MODULES_JSON is None:
        print("Warning: hycal_map.json not found", file=sys.stderr)
        return {}

    geo = _load_json_geometry(MODULES_JSON)
    data: Dict[str, ModuleData] = {}
    for name, (x, y, row, col) in geo.items():
        data[name] = ModuleData(name, x, y, row, col)

    if root_path is not None and root_path.is_file() and HAS_UPROOT:
        try:
            with uproot.open(str(root_path)) as f:
                dir3 = f.get("energy_3p5GeV")
                dir0 = f.get("energy_0p7GeV")
                for name, d in data.items():
                    if dir3 is not None:
                        try:
                            h = dir3[f"h_energy_{name}_3p5"]
                            d.hist_3p5 = h.to_numpy()
                        except Exception:
                            pass
                    if dir0 is not None:
                        try:
                            h = dir0[f"h_energy_{name}_0p7"]
                            d.hist_0p7 = h.to_numpy()
                        except Exception:
                            pass
        except Exception as e:
            print(f"ROOT read error: {e}", file=sys.stderr)
    elif root_path is not None and not HAS_UPROOT:
        print("Warning: uproot not installed — cannot read ROOT file.", file=sys.stderr)

    # Fit all peaks
    for d in data.values():
        if d.hist_3p5 is not None:
            counts, edges = d.hist_3p5
            for key in ("ep_3p5", "ee_3p5"):
                E, sig = d.e_exp[key], d.sigma_exp[key]
                if E > 0.0:
                    mu, sigma_f, amp = fit_peak_full(counts, edges, E, sig)
                    d.peaks[key]      = mu
                    d.peak_sigma[key] = sigma_f
                    d.peak_amp[key]   = amp

        if d.hist_0p7 is not None:
            counts, edges = d.hist_0p7
            for key in ("ep_0p7", "ee_0p7"):
                E, sig = d.e_exp[key], d.sigma_exp[key]
                # Skip ee_0p7 when peaks overlap (same cut as nonlinearity.cpp)
                if key == "ee_0p7" and abs(d.e_exp["ep_0p7"] - E) <= 170.0:
                    continue
                if E > 0.0:
                    mu, sigma_f, amp = fit_peak_full(counts, edges, E, sig)
                    d.peaks[key]      = mu
                    d.peak_sigma[key] = sigma_f
                    d.peak_amp[key]   = amp

        d.recompute_scale()
        d.refit_nl()

    return data


# =============================================================================
# Matplotlib canvas  (4 peak histograms on top row + 1 large linearity plot)
# =============================================================================

class TripleCanvas(FigureCanvas):
    def __init__(self):
        fc = THEME.CANVAS
        self.fig = Figure(facecolor=fc)
        # Top row: 4 individual peak histograms; bottom: linearity (2× height)
        gs = self.fig.add_gridspec(2, 4, height_ratios=[1, 2],
                                   hspace=0.65, wspace=0.50,
                                   left=0.07, right=0.97, top=0.96, bottom=0.07)
        self.axes: List = [self.fig.add_subplot(gs[0, i]) for i in range(4)]
        self.axes.append(self.fig.add_subplot(gs[1, :]))
        for ax in self.axes:
            ax.set_facecolor(fc)
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def _style(self, ax, title="", xlabel="", ylabel=""):
        tc = THEME.TEXT
        dc = THEME.TEXT_DIM
        ax.set_title(title, color=tc, fontsize=9, pad=3)
        ax.set_xlabel(xlabel, color=dc, fontsize=8)
        ax.set_ylabel(ylabel, color=dc, fontsize=8)
        ax.tick_params(colors=dc, labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(dc)

    def redraw(self):
        self.fig.canvas.draw_idle()


# =============================================================================
# Main window
# =============================================================================

class NLViewerWindow(QMainWindow):

    def __init__(self, root_path: Optional[Path] = None):
        super().__init__()
        self.setWindowTitle("HyCal Non-linearity Viewer")
        self.resize(1480, 940)

        self._data: Dict[str, ModuleData] = {}
        self._cur:  Optional[ModuleData]  = None
        self._w_modules: List[Module] = []   # W-only Module list for the map

        # SpanSelector state: per peak key → (xmin, xmax) or None
        self._span_range: Dict[str, Optional[Tuple[float, float]]] = {
            k: None for k in _PEAK_KEYS}
        self._span_selectors: List = []   # keep references to avoid GC

        # scan state
        self._scan_timer  = QTimer(self)
        self._scan_timer.setInterval(500)   # ms per module
        self._scan_timer.timeout.connect(self._scan_step)
        self._scan_index  = 0
        self._scan_names: List[str] = []

        # track active map mode to detect changes and preserve user range
        self._current_map_mode: str = ""

        self._build_ui()
        if root_path is not None:
            self._load(root_path)

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root_lay = QVBoxLayout(central)
        root_lay.setContentsMargins(4, 4, 4, 4)
        root_lay.setSpacing(4)

        # top bar
        top = QHBoxLayout()
        btn = QPushButton("Open ROOT file…")
        btn.clicked.connect(self._open_file)
        top.addWidget(btn)
        self._file_lbl = QLabel("No file loaded")
        self._file_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        top.addWidget(self._file_lbl, stretch=1)
        top.addWidget(QLabel("Color map:"))
        self._mode_cb = QComboBox()
        self._mode_cb.addItems([
            "Non-linearity (nl)",
            "Mean |corr/exp-1| (1st order)",
            "Mean |corr/exp-1| (2nd order)",
            "Flat",
        ])
        self._mode_cb.currentIndexChanged.connect(self._on_mode_changed)
        top.addWidget(self._mode_cb)

        top.addSpacing(12)
        self._scan_btn = QPushButton("▶ Scan")
        self._scan_btn.setToolTip(
            "Auto-scan through all W modules.\n"
            "Press Enter or click the button again to stop.")
        self._scan_btn.setEnabled(False)
        self._scan_btn.clicked.connect(self._toggle_scan)
        top.addWidget(self._scan_btn)

        top.addWidget(QLabel("Speed:"))
        self._scan_speed_cb = QComboBox()
        self._scan_speed_cb.setFixedWidth(90)
        for label, ms in [("0.2 s", 200), ("0.5 s", 500), ("1 s", 1000),
                          ("2 s", 2000), ("5 s", 5000)]:
            self._scan_speed_cb.addItem(label, ms)
        self._scan_speed_cb.setCurrentIndex(1)   # default 0.5 s
        self._scan_speed_cb.currentIndexChanged.connect(self._on_scan_speed_changed)
        top.addWidget(self._scan_speed_cb)

        self._scan_lbl = QLabel("")
        self._scan_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 11px;")
        top.addWidget(self._scan_lbl)

        top.addSpacing(12)
        self._export_btn = QPushButton("↩ Export nl → JSON")
        self._export_btn.setToolTip(
            "Read a base calibration JSON, fill in non_linear values from the current\n"
            "analysis results, and save a new JSON in the current directory.")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._export_nl)
        top.addWidget(self._export_btn)

        root_lay.addLayout(top)

        # Enter key stops scan
        stop_sc = QShortcut(QKeySequence(Qt.Key.Key_Return), self)
        stop_sc.activated.connect(self._stop_scan)
        stop_sc2 = QShortcut(QKeySequence(Qt.Key.Key_Enter), self)
        stop_sc2.activated.connect(self._stop_scan)

        # horizontal splitter: map | plots+controls
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_lay.addWidget(splitter, stretch=1)

        # ── left: map ─────────────────────────────────────────────────────
        map_wrap = QWidget()
        map_lay  = QVBoxLayout(map_wrap)
        map_lay.setContentsMargins(0, 0, 0, 0)
        self._map = HyCalMapWidget(map_wrap, enable_zoom_pan=True, min_size=(440, 440))
        self._map.moduleClicked.connect(self._on_click)
        map_lay.addWidget(self._map)
        self._range_ctrl = ColorRangeControl(self._map, orientation="horizontal",
                                             parent=map_wrap)
        map_lay.addWidget(self._range_ctrl)
        splitter.addWidget(map_wrap)

        # ── right: plots + controls ────────────────────────────────────────
        right     = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(4, 0, 4, 4)
        right_lay.setSpacing(4)

        self._mod_lbl = QLabel("← click a W module on the map")
        self._mod_lbl.setFont(QFont("Monospace", 12, QFont.Weight.Bold))
        self._mod_lbl.setStyleSheet(f"color: {THEME.ACCENT};")
        right_lay.addWidget(self._mod_lbl)

        self._canvas = TripleCanvas()
        right_lay.addWidget(self._canvas, stretch=1)

        right_lay.addWidget(self._build_refit_panel())

        splitter.addWidget(right)
        splitter.setSizes([460, 900])

    def _build_refit_panel(self) -> QGroupBox:
        box = QGroupBox("Peak re-fit  (drag on histogram to auto-fit)")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(8, 10, 8, 6)
        lay.setSpacing(4)

        row1 = QHBoxLayout()

        hint_lbl = QLabel("Drag on any peak histogram to auto-fit that peak.")
        hint_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 11px;")
        row1.addWidget(hint_lbl)
        row1.addSpacing(12)
        row1.addWidget(QLabel("Range:"))
        self._range_lbl = QLabel("(drag to set → auto-fit)")
        self._range_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 11px;")
        row1.addWidget(self._range_lbl)

        self._clear_range_btn = QPushButton("✕ Clear")
        self._clear_range_btn.setFixedWidth(60)
        self._clear_range_btn.setEnabled(False)
        self._clear_range_btn.setToolTip("Clear the drag-selected fit range and revert to auto (Eexp ± 6σ)")
        self._clear_range_btn.clicked.connect(self._clear_span_range)
        row1.addWidget(self._clear_range_btn)

        row1.addSpacing(16)

        self._refit_peak_btn = QPushButton("Refit All Peaks")
        self._refit_peak_btn.setEnabled(False)
        self._refit_peak_btn.setToolTip(
            "Re-run 3-step auto Gaussian fit on all 4 peaks (clears manual drag ranges).\n"
            "Updates the linearity plot.")
        self._refit_peak_btn.clicked.connect(self._do_refit_all_peaks)
        row1.addWidget(self._refit_peak_btn)

        row1.addSpacing(16)

        self._refit_nl_btn = QPushButton("Refit Linearity")
        self._refit_nl_btn.setEnabled(False)
        self._refit_nl_btn.setToolTip(
            "Re-run the non-linearity fit using current measurement points.\n"
            "1st: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000\n"
            "2nd: + nl2*((E_rec-E_base)/1000)^2")
        self._refit_nl_btn.clicked.connect(self._do_refit_nl)
        row1.addWidget(self._refit_nl_btn)

        row1.addStretch()
        lay.addLayout(row1)

        self._status_lbl = QLabel("(select a module first)")
        self._status_lbl.setFont(QFont("Monospace", 9))
        self._status_lbl.setWordWrap(True)
        self._status_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        lay.addWidget(self._status_lbl)

        return box

    # ── file loading ───────────────────────────────────────────────────────

    def _open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open nonlinearity ROOT file", "",
            "ROOT files (*.root);;All files (*)")
        if path:
            self._load(Path(path))

    def _after_load(self):
        """Enable buttons once data is loaded."""
        self._scan_btn.setEnabled(bool(self._data))
        self._export_btn.setEnabled(bool(self._data))

    def _load(self, path: Path):
        self._file_lbl.setText(f"Loading {path.name} …")
        QApplication.processEvents()
        self._data   = load_all(path)
        self._cur    = None
        self._file_lbl.setText(str(path))
        self._mod_lbl.setText("← click a W module on the map")
        self._refit_peak_btn.setEnabled(False)
        self._refit_nl_btn.setEnabled(False)
        self._status_lbl.setText("")
        self._span_range = {k: None for k in _PEAK_KEYS}

        # Build W-only module list for the map
        if MODULES_JSON is not None:
            all_mods = load_modules(MODULES_JSON)
            self._w_modules = [m for m in all_mods if m.mod_type == "PbWO4"]
        else:
            # Fallback: build synthetic Module objects from ModuleData
            self._w_modules = [
                Module(d.name, "PbWO4", d.x, d.y, 20.8, 20.8)
                for d in self._data.values()
            ]

        self._map.set_modules(self._w_modules)
        self._refresh_map()
        self._clear_plots()
        self._after_load()

    # ── map ───────────────────────────────────────────────────────────────

    def _on_mode_changed(self):
        """Called only when the user switches the colour-map combo.
        Forces an auto range-reset."""
        self._current_map_mode = ""   # flag: treat as mode change
        self._refresh_map()

    def _refresh_map(self):
        """Update map values.  Auto-resets the colour range only when the
        mode has actually changed; otherwise preserves the user-set range
        so manual edits survive data updates (peak refit, etc.)."""
        if not self._data:
            return
        mode = self._mode_cb.currentText()
        mode_changed = (mode != self._current_map_mode)
        self._current_map_mode = mode

        if "Non-linearity" in mode:
            values = {name: d.nl for name, d in self._data.items()}
            if mode_changed:
                vals = [v for v in values.values() if v != 0.0]
                if vals:
                    arr  = np.array(vals)
                    mean = float(np.mean(arr))
                    std  = float(np.std(arr))
                    self._range_ctrl.set_range(mean - 3.0 * std,
                                               mean + 3.0 * std)
        elif "1st order" in mode:
            values = {name: d.mean_corr_resid_1st for name, d in self._data.items()}
            if mode_changed:
                vals = [v for v in values.values() if v > 0.0]
                vmax = float(np.percentile(np.array(vals), 95)) if vals else 0.01
                self._range_ctrl.set_range(0.0, vmax)
        elif "2nd order" in mode:
            values = {name: d.mean_corr_resid_2nd for name, d in self._data.items()}
            if mode_changed:
                vals = [v for v in values.values() if v > 0.0]
                vmax = float(np.percentile(np.array(vals), 95)) if vals else 0.01
                self._range_ctrl.set_range(0.0, vmax)
        else:
            values = {name: 0.0 for name in self._data}
            if mode_changed:
                self._range_ctrl.set_range(0.0, 1.0)
        self._map.set_values(values)

    def _on_click(self, name: str):
        if not name:
            return
        d = self._data.get(name)
        if d is None:
            return
        self._cur = d
        # Clear span ranges when switching modules
        self._span_range = {k: None for k in _PEAK_KEYS}
        self._range_lbl.setText("(auto — drag to set)")
        self._clear_range_btn.setEnabled(False)
        info = (f"Module: {name}   θ={d.theta:.3f}°   "
                f"nl = {d.nl:.4f} ± {d.nl_err:.4f}   "
                f"χ²/ndf = {d.chi2_nl:.2f}/{d.ndf_nl}")
        self._mod_lbl.setText(info)
        self._refit_peak_btn.setEnabled(True)
        self._refit_nl_btn.setEnabled(True)
        self._status_lbl.setText("")
        self._draw_module(d)

    def _clear_span_range(self):
        """Clear all drag-selected fit ranges."""
        self._span_range = {k: None for k in _PEAK_KEYS}
        self._range_lbl.setText("(auto — drag to set)")
        self._clear_range_btn.setEnabled(False)

    # ── drawing ────────────────────────────────────────────────────────────

    def _clear_plots(self):
        for ax in self._canvas.axes:
            ax.cla()
            ax.set_facecolor(THEME.CANVAS)
        self._canvas.redraw()

    def _draw_module(self, d: ModuleData):
        axes = self._canvas.axes
        for ax in axes:
            ax.cla()
            ax.set_facecolor(THEME.CANVAS)
            ax.tick_params(colors=THEME.TEXT_DIM, labelsize=7)
            for sp in ax.spines.values():
                sp.set_color(THEME.BORDER)

        # 4 individual peak histograms on top row
        ebeam_lbl = {"3p5": "3.5 GeV", "0p7": "0.7 GeV"}
        for i, key in enumerate(_PEAK_KEYS):
            self._draw_energy_hist_single(axes[i], d, key)
            self._canvas._style(axes[i],
                title=f"{_PEAK_LABELS[key]}  ({ebeam_lbl[_PEAK_TAG[key]]})",
                xlabel="E (MeV)", ylabel="Counts")

        # Large linearity plot on bottom row
        self._draw_linearity(axes[4], d)
        self._canvas._style(axes[4], title="Non-linearity",
                            xlabel="E_rec (MeV)", ylabel="E_rec / E_exp")
        self._canvas.redraw()
        self._setup_span_selectors(d)

    def _draw_energy_hist_single(self, ax, d: ModuleData, key: str):
        """Draw histogram for one peak key with x-range [Eexp ± 8σ]."""
        tag  = _PEAK_TAG[key]
        hist = d.hist_3p5 if tag == "3p5" else d.hist_0p7
        if hist is None:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, color=THEME.TEXT_DIM, fontsize=8)
            return

        counts, edges = hist
        centers = 0.5 * (edges[:-1] + edges[1:])
        color = _PEAK_COLOR[key]
        Eexp  = d.e_exp[key]
        sig   = d.sigma_exp[key]

        # X range: [Eexp ± 8σ], clipped to histogram bounds
        xlo = max(float(edges[0]),  Eexp - 8.0 * sig)
        xhi = min(float(edges[-1]), Eexp + 8.0 * sig)
        ax.set_xlim(xlo, xhi)

        ax.step(centers, counts, where="mid",
                color=THEME.TEXT, linewidth=0.9, alpha=0.85)

        # Drag-selected span range band
        span = self._span_range.get(key)
        if span is not None:
            ax.axvspan(span[0], span[1], alpha=0.18,
                       color=THEME.ACCENT, linewidth=0, zorder=2)
            ax.axvline(span[0], color=THEME.ACCENT,
                       linewidth=1.2, linestyle="-", alpha=0.8, zorder=3)
            ax.axvline(span[1], color=THEME.ACCENT,
                       linewidth=1.2, linestyle="-", alpha=0.8, zorder=3)

        if Eexp > 0.0:
            # Search window shading [Eexp ± 6σ]
            ax.axvspan(Eexp - 6.0 * sig, Eexp + 6.0 * sig,
                       alpha=0.07, color=color, linewidth=0)
            # Expected energy marker
            ax.axvline(Eexp, color=color, linewidth=0.8, linestyle=":", alpha=0.5)
            ax.text(0.97, 0.97, f"exp={Eexp:.0f}\nσ={sig:.0f}",
                    transform=ax.transAxes, fontsize=6, ha="right", va="top",
                    color=THEME.TEXT_DIM)

        # Fitted Gaussian overlay
        mu   = d.peaks[key]
        amp  = d.peak_amp[key]
        sigf = d.peak_sigma[key]
        if mu > 0.0 and amp > 0.0 and sigf > 0.0:
            xg = np.linspace(max(xlo, mu - 4.0 * sigf),
                             min(xhi, mu + 4.0 * sigf), 300)
            yg = _gauss(xg, amp, mu, sigf)
            ax.plot(xg, yg, color=color, linewidth=2.0)
            ax.axvline(mu, color=color, linewidth=1.2, linestyle="--", alpha=0.7)
            ax.annotate(f"{mu:.1f} MeV",
                        xy=(mu, amp), xytext=(0, 6), textcoords="offset points",
                        color=color, fontsize=7, ha="center", va="bottom")
        elif Eexp > 0.0:
            ax.axvline(Eexp, color=color, linewidth=1.0, linestyle=":", alpha=0.5)

    def _draw_linearity(self, ax, d: ModuleData):
        E_base = d.e_exp["ep_3p5"]

        # Collect valid (key, E_exp, E_rec_scaled) pairs
        pts: List[Tuple[str, float, float]] = []
        for key in _PEAK_KEYS:
            p = d.scaled(key)
            if p > 0.0:
                pts.append((key, d.e_exp[key], p))

        if not pts:
            ax.text(0.5, 0.5, "No valid peaks found", ha="center", va="center",
                    transform=ax.transAxes, color=THEME.TEXT_DIM, fontsize=9)
            return

        E_rec_arr = np.array([p[2] for p in pts])
        E_exp_arr = np.array([p[1] for p in pts])
        ratio_arr = E_rec_arr / np.where(E_exp_arr > 0, E_exp_arr, 1.0)

        # Fixed y range [0.85, 1.15]
        ax.set_ylim(0.85, 1.15)

        # x plot range
        margin = (E_rec_arr.max() - E_rec_arr.min()) * 0.15 if len(E_rec_arr) > 1 else E_rec_arr.max() * 0.15
        xplot = np.linspace(max(0, E_rec_arr.min() - margin), E_rec_arr.max() + margin, 300)

        # Perfect linearity reference y = 1
        ax.axhline(1.0, color="tomato", linestyle="--",
                   linewidth=1.5, alpha=0.75, label="y = 1 (ideal)")

        # 1st order fit curve
        if E_base > 0.0 and len(pts) >= 2:
            y1fit = nl_ratio_1st(xplot, d.nl, E_base)
            ax.plot(xplot, y1fit, color="cornflowerblue", linewidth=2.0,
                    label="1st order fit")

        # 2nd order fit curve
        if E_base > 0.0 and len(pts) >= 3 and d.ndf_nl2 >= 0:
            y2fit = nl_ratio_2nd(xplot, d.nl2_1, d.nl2_2, E_base)
            ax.plot(xplot, y2fit, color="violet", linewidth=2.0,
                    linestyle="--", label="2nd order fit")

        # Scatter: measured points
        for key, xe, ye_rec in pts:
            ratio = ye_rec / xe if xe > 0 else 0.0
            ax.scatter(ye_rec, ratio, color=_PEAK_COLOR[key], s=70, zorder=5,
                       edgecolors="white", linewidths=0.5)
            ax.annotate(_PEAK_LABELS[key], xy=(ye_rec, ratio),
                        xytext=(6, 4), textcoords="offset points",
                        color=_PEAK_COLOR[key], fontsize=6.5)

        # Corrected points — 1st order: E_corr = E_rec / (1 + nl*(E_rec-E_base)/1000)
        if E_base > 0.0 and len(pts) >= 2:
            for key, xe, ye_rec in pts:
                denom = 1.0 + d.nl * (ye_rec - E_base) / 1000.0
                E_corr = ye_rec / denom if denom != 0.0 else ye_rec
                ratio_corr = E_corr / xe if xe > 0 else 0.0
                ax.scatter(E_corr, ratio_corr, color="limegreen", s=55, zorder=6,
                           marker="o", edgecolors="white", linewidths=0.5)

        # Corrected points — 2nd order
        if E_base > 0.0 and len(pts) >= 3 and d.ndf_nl2 >= 0:
            for key, xe, ye_rec in pts:
                t = (ye_rec - E_base) / 1000.0
                denom = 1.0 + d.nl2_1 * t + d.nl2_2 * t * t
                E_corr2 = ye_rec / denom if denom != 0.0 else ye_rec
                ratio_corr2 = E_corr2 / xe if xe > 0 else 0.0
                ax.scatter(E_corr2, ratio_corr2, color="orange", s=55, zorder=6,
                           marker="s", edgecolors="white", linewidths=0.5)

        # Dummy legend entries for point types
        ax.scatter([], [], color="tomato",         s=50, label="e-p (meas.)")
        ax.scatter([], [], color="cornflowerblue", s=50, label="e-e (meas.)")
        ax.scatter([], [], color="limegreen",      s=50, marker="o", label="corrected (1st)")
        ax.scatter([], [], color="orange",         s=50, marker="s", label="corrected (2nd)")

        # Fit result annotation
        if E_base > 0.0:
            chi1_str = (f"χ²/ndf={d.chi2_nl:.2f}/{d.ndf_nl}"
                        if d.ndf_nl > 0 else "")
            chi2_str = (f"χ²/ndf={d.chi2_nl2:.2f}/{d.ndf_nl2}"
                        if d.ndf_nl2 > 0 else "")
            note = (
                f"1st: nl₁={d.nl:.4f}±{d.nl_err:.4f}  {chi1_str}\n"
                f"2nd: nl₁={d.nl2_1:.4f}±{d.nl2_1_err:.4f}"
                f"  nl₂={d.nl2_2:.4f}±{d.nl2_2_err:.4f}\n"
                f"     {chi2_str}\n"
                f"E_base={E_base:.1f} MeV"
            )
            ax.text(0.03, 0.97, note,
                    transform=ax.transAxes, fontsize=7.0,
                    va="top", ha="left", color="lightcyan",
                    bbox=dict(boxstyle="round,pad=0.35",
                              fc="#1a2a3a", ec=THEME.BORDER, alpha=0.85))

        # Two side-by-side legends in the upper-right corner
        _leg_kw = dict(facecolor=THEME.PANEL, edgecolor=THEME.BORDER,
                       labelcolor=THEME.TEXT, fontsize=13)
        all_h, all_l = ax.get_legend_handles_labels()
        _fit_labels  = {"y = 1 (ideal)", "1st order fit", "2nd order fit"}
        _pt_labels   = {"e-p (meas.)", "e-e (meas.)", "corrected (1st)", "corrected (2nd)"}
        h_fit = [h for h, l in zip(all_h, all_l) if l in _fit_labels]
        l_fit = [l for l in all_l if l in _fit_labels]
        h_pt  = [h for h, l in zip(all_h, all_l) if l in _pt_labels]
        l_pt  = [l for l in all_l if l in _pt_labels]
        # Left legend: fits + reference  (anchored at 77% from left)
        leg_left = ax.legend(h_fit, l_fit,
                             loc="upper right", bbox_to_anchor=(0.77, 1.0),
                             **_leg_kw)
        ax.add_artist(leg_left)
        # Right legend: point types  (anchored at right edge)
        ax.legend(h_pt, l_pt,
                  loc="upper right", bbox_to_anchor=(1.0, 1.0),
                  **_leg_kw)

    # ── span selector ──────────────────────────────────────────────────────

    def _setup_span_selectors(self, d: ModuleData):
        """Attach SpanSelectors to all 4 peak histogram axes.
        Must be called after every canvas redraw (axes are recreated)."""
        self._span_selectors.clear()
        for ax, key in zip(self._canvas.axes[:4], _PEAK_KEYS):
            tag  = _PEAK_TAG[key]
            hist = d.hist_3p5 if tag == "3p5" else d.hist_0p7
            if hist is None:
                continue
            sel = SpanSelector(
                ax,
                lambda xmin, xmax, _key=key: self._on_span_select(_key, xmin, xmax),
                direction="horizontal",
                useblit=True,
                props=dict(alpha=0.25,
                           facecolor=THEME.ACCENT if hasattr(THEME, 'ACCENT') else "#007aff"),
                interactive=True,
                drag_from_anywhere=True,
            )
            self._span_selectors.append(sel)

    def _on_span_select(self, key: str, xmin: float, xmax: float):
        if xmax - xmin < 1.0:
            return
        self._span_range[key] = (xmin, xmax)
        self._range_lbl.setText(
            f"{_PEAK_LABELS.get(key, key)}:  [{xmin:.1f},  {xmax:.1f}] MeV")
        self._clear_range_btn.setEnabled(True)
        ok = self._do_refit_peak(key)
        d = self._cur
        if d is None:
            return
        if ok:
            d.refit_nl()
            self._refresh_map()
        self._draw_module(d)
        self._mod_lbl.setText(
            f"Module: {d.name}   θ={d.theta:.3f}°   "
            f"nl={d.nl:.4f}±{d.nl_err:.4f}   "
            f"χ²/ndf={d.chi2_nl:.2f}/{d.ndf_nl}")

    # ── re-fit actions ─────────────────────────────────────────────────────

    def _do_refit_peak(self, key: str) -> bool:
        """Fit peak for *key*. Returns True on success. Does NOT redraw."""
        d = self._cur
        if d is None:
            return False
        tag  = _PEAK_TAG[key]
        hist = d.hist_3p5 if tag == "3p5" else d.hist_0p7
        if hist is None:
            self._status_lbl.setText(f"No histogram for {_PEAK_LABELS[key]}.")
            return False
        Eexp = d.e_exp[key]
        sig  = d.sigma_exp[key]
        if Eexp <= 0.0:
            return False
        counts, edges = hist
        span = self._span_range.get(key)
        if span is not None and HAS_SCIPY:
            xmin, xmax = span
            centers = 0.5 * (edges[:-1] + edges[1:])
            mask = (centers >= xmin) & (centers <= xmax) & (counts > 0)
            if mask.sum() < 4:
                self._status_lbl.setText(
                    f"Not enough bins in [{xmin:.0f}, {xmax:.0f}] MeV — try wider.")
                return False
            xf, yf = centers[mask], counts[mask].astype(float)
            amp0 = float(yf.max())
            mu0  = float((xf * yf).sum() / yf.sum())
            sig0 = max((xmax - xmin) / 4.0, float(edges[1] - edges[0]))
            try:
                from scipy.optimize import curve_fit as _cf
                popt, _ = _cf(
                    _gauss, xf, yf,
                    p0=[amp0, mu0, sig0],
                    bounds=([0., xmin, 1e-6], [np.inf, xmax, (xmax - xmin)]),
                    maxfev=3000,
                )
                mu, sigma_f, amp = float(popt[1]), abs(float(popt[2])), float(popt[0])
            except Exception as exc:
                self._status_lbl.setText(f"Fit failed: {exc}")
                return False
            range_note = f"  (drag [{xmin:.0f}, {xmax:.0f}] MeV)"
        else:
            mu, sigma_f, amp = fit_peak_full(counts, edges, Eexp, sig)
            range_note = f"  (auto [{Eexp - 6*sig:.0f}, {Eexp + 6*sig:.0f}] MeV)"
        if mu <= 0.0:
            self._status_lbl.setText(
                f"Fit FAILED for {_PEAK_LABELS[key]}{range_note}.")
            return False
        d.peaks[key]      = mu
        d.peak_sigma[key] = sigma_f
        d.peak_amp[key]   = amp
        d.recompute_scale()
        self._status_lbl.setText(
            f"{_PEAK_LABELS[key]}: μ={mu:.2f} MeV, σ={sigma_f:.2f} MeV{range_note}")
        return True

    def _do_refit_all_peaks(self):
        """Refit all 4 peaks with auto windows (clears manual ranges)."""
        d = self._cur
        if d is None:
            return
        self._span_range = {k: None for k in _PEAK_KEYS}
        any_ok = False
        for key in _PEAK_KEYS:
            any_ok |= self._do_refit_peak(key)
        if any_ok:
            d.refit_nl()
            self._refresh_map()
            self._draw_module(d)
            self._mod_lbl.setText(
                f"Module: {d.name}   θ={d.theta:.3f}°   "
                f"nl={d.nl:.4f}±{d.nl_err:.4f}   "
                f"χ²/ndf={d.chi2_nl:.2f}/{d.ndf_nl}")
            self._status_lbl.setText("All 4 peaks refit with auto windows.")

    # ── scan ──────────────────────────────────────────────────────────────

    def _export_nl(self):
        """Export nl coefficients to a new calibration JSON.

        Workflow:
          1. Ask user to select a *base* calibration JSON (for factor/base_energy values).
          2. Patch non_linear field for every W module that has a valid nl result.
          3. Save as a new file (same dir, same name with _nonLinear suffix, or user-chosen path).
        """
        if not self._data:
            return

        # ── Step 1: choose base file ──────────────────────────────────────────
        calib_dir = str(
            _SCRIPT_DIR.parent / "database" / "calibration"
            if (_SCRIPT_DIR.parent / "database" / "calibration").is_dir()
            else Path(".")
        )
        src_path, _ = QFileDialog.getOpenFileName(
            self, "Select base calibration JSON", calib_dir,
            "JSON files (*.json);;All files (*)")
        if not src_path:
            return
        src_path = Path(src_path)

        try:
            with open(src_path) as f:
                entries = json.load(f)
        except Exception as exc:
            self._status_lbl.setText(f"Failed to read {src_path.name}: {exc}")
            return

        # ── Step 2: patch non_linear / nl0 / nl1 / nl2 ───────────────────────
        n_patched = 0
        for entry in entries:
            name = entry.get("name", "")
            d = self._data.get(name)
            if d is not None and d.nl != 0.0:
                entry["non_linear"] = round(float(d.nl), 6)
                entry["nl0"] = round(float(d.nl), 6)
                entry["nl1"] = round(float(d.nl2_1), 6)
                entry["nl2"] = round(float(d.nl2_2), 6)
                n_patched += 1

        # ── Step 3: choose output path ────────────────────────────────────────
        stem = src_path.stem
        if not stem.endswith("_nonLinear"):
            stem += "_nonLinear"
        default_out = str(src_path.parent / (stem + ".json"))
        out_path, _ = QFileDialog.getSaveFileName(
            self, "Save new calibration JSON", default_out,
            "JSON files (*.json);;All files (*)")
        if not out_path:
            return

        try:
            with open(out_path, "w") as f:
                json.dump(entries, f, indent=2)
            self._status_lbl.setText(
                f"Exported nl for {n_patched} W modules (nl0/nl1/nl2)  →  {Path(out_path).name}")
        except Exception as exc:
            self._status_lbl.setText(f"Save failed: {exc}")

    def _on_scan_speed_changed(self, _idx: int):
        ms = self._scan_speed_cb.currentData()
        self._scan_timer.setInterval(ms)

    def _toggle_scan(self):
        if self._scan_timer.isActive():
            self._stop_scan()
        else:
            self._start_scan()

    def _start_scan(self):
        # Build sorted name list from all loaded W modules
        all_names = sorted(self._data.keys(),
                           key=lambda n: int(n[1:]) if n[1:].isdigit() else 0)
        if not all_names:
            return
        # Start from current module if possible
        cur_name = self._cur.name if self._cur else ""
        try:
            self._scan_index = all_names.index(cur_name)
        except ValueError:
            self._scan_index = 0
        self._scan_names = all_names
        self._scan_btn.setText("■ Stop")
        self._scan_timer.setInterval(self._scan_speed_cb.currentData())
        self._scan_timer.start()
        self._scan_step()   # show first module immediately

    def _stop_scan(self):
        if not self._scan_timer.isActive():
            return
        self._scan_timer.stop()
        self._scan_btn.setText("▶ Scan")
        self._scan_lbl.setText("")

    def _scan_step(self):
        names = self._scan_names
        if not names or self._scan_index >= len(names):
            self._stop_scan()
            return
        name = names[self._scan_index]
        # Highlight on map and draw module plots via the normal click handler
        self._map._hovered = name
        self._map.update()
        self._on_click(name)
        total = len(names)
        self._scan_lbl.setText(f"{name}  ({self._scan_index + 1}/{total})")
        self._scan_index += 1
        if self._scan_index >= total:
            self._scan_timer.stop()
            self._scan_btn.setText("▶ Scan")
            self._scan_lbl.setText(f"Done ({total} modules)")

    def _do_refit_nl(self):
        d = self._cur
        if d is None:
            return

        d.refit_nl()
        self._refresh_map()
        self._draw_module(d)
        self._status_lbl.setText(
            f"Linearity refit:  nl = {d.nl:.4f} ± {d.nl_err:.4f},  "
            f"χ²/ndf = {d.chi2_nl:.2f}/{d.ndf_nl}")
        self._mod_lbl.setText(
            f"Module: {d.name}   θ={d.theta:.3f}°   "
            f"nl = {d.nl:.4f} ± {d.nl_err:.4f}   "
            f"χ²/ndf = {d.chi2_nl:.2f}/{d.ndf_nl}")


# =============================================================================
# Entry point
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description="HyCal PWO4 non-linearity analysis viewer")
    ap.add_argument("root_file", nargs="?", type=Path,
                    help="Path to nonlinearity_results.root (output of prad2ana_nonlinearity)")
    ap.add_argument("--theme", choices=available_themes(), default="dark")
    args = ap.parse_args()

    set_theme(args.theme)

    app = QApplication(sys.argv)
    app.setApplicationName("NL Viewer")
    win = NLViewerWindow(root_path=args.root_file)
    apply_theme_palette(win)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
