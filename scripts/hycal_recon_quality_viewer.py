#!/usr/bin/env python3
"""
HyCal Reconstruction Quality Viewer

Features:
- Scan a folder of quick_check ROOT files (prad_*_quick_check.root).
- Fit per-module energy spectra (module_energy/h_<module>) with Gaussian.
- Cache fit results per run to JSON files in a dedicated cache folder.
- Load from cache on later startup without refitting.
- Full refit button for all runs and all W modules.
- Fit settings dialog and threshold dialog.
- Left panel: W-module map with selectable metric.
- Right panel: selected-module run-by-run center/resolution and current-run spectrum.
- Manual drag-refit on spectrum updates the selected run+module cache entry.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import site
import sys
import tempfile
import traceback
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
try:
    import uproot
    HAS_UPROOT = True
except Exception:
    HAS_UPROOT = False

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


def _fix_qt_lib_path() -> None:
    """Match replay_viewer startup workaround for Qt ABI conflicts.

    Re-exec once with PyQt6 bundled Qt6/lib prepended to LD_LIBRARY_PATH so
    the dynamic linker consistently resolves Qt symbols from one runtime.
    """
    sp_list: list[str] = []
    try:
        sp_list += site.getsitepackages()
    except AttributeError:
        pass
    try:
        sp_list.append(site.getusersitepackages())
    except Exception:
        pass

    for sp in sp_list:
        qt6_lib = os.path.join(sp, "PyQt6", "Qt6", "lib")
        if os.path.isdir(qt6_lib):
            current = os.environ.get("LD_LIBRARY_PATH", "")
            entries = [e for e in current.split(":") if e]
            if qt6_lib not in entries:
                new_path = qt6_lib + ((":" + current) if current else "")
                env = os.environ.copy()
                env["LD_LIBRARY_PATH"] = new_path
                os.execve(sys.executable, [sys.executable] + sys.argv, env)
            return


_fix_qt_lib_path()

try:
    from PyQt6.QtCore import QThread, Qt, pyqtSignal
    from PyQt6.QtGui import QFont
    from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
        QSpinBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QProgressBar,
    QSplitter,
    QVBoxLayout,
    QWidget,
    )
except Exception as exc:
    print("Failed to import PyQt6 runtime.", file=sys.stderr)
    print(f"Reason: {exc}", file=sys.stderr)
    print("Try:", file=sys.stderr)
    print("  python3 -m pip install --user --upgrade PyQt6 PyQt6-Qt6 PyQt6-sip", file=sys.stderr)
    print("Or run with:", file=sys.stderr)
    print("  LD_LIBRARY_PATH=$HOME/.local/lib/python3.12/site-packages/PyQt6/Qt6/lib python3 scripts/hycal_recon_quality_viewer.py", file=sys.stderr)
    raise

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hycal_geoview import (  # noqa: E402
    HyCalMapWidget,
    THEME,
    apply_theme_palette,
    available_themes,
    load_modules,
    set_theme,
)


RUN_FILE_RE = re.compile(r"prad_(\d+)_quick_check\.root$")
_HIST_KEY_CACHE: Dict[str, Dict[str, str]] = {}
FIT_CACHE_SCHEMA_VERSION = 5


@dataclass
class FitSettings:
    window_fraction: float = 0.35
    min_half_window: float = 40.0
    max_half_window: float = 250.0
    min_entries: int = 30


@dataclass
class ThresholdSettings:
    offset_pct: float = 3.0
    resolution_change_pct: float = 5.0
    jump_pct: float = 5.0
    slope_abs_threshold: float = 0.4


@dataclass
class RunMeta:
    run: int
    root_path: str
    root_mtime: float
    ebeam: float
    hycal_z: float


@dataclass
class RunGroup:
    key_run: int
    label: str
    members: List[RunMeta]
    ebeam: float
    hycal_z: float
    root_mtime: float


class QualityMapWidget(HyCalMapWidget):
    def __init__(self, parent=None):
        super().__init__(parent, enable_zoom_pan=True, min_size=(520, 520))
        self._metric_label = ""

    def set_metric_label(self, label: str):
        self._metric_label = label

    def _colorbar_center_text(self) -> str:
        return self._metric_label if self._metric_label else super()._colorbar_center_text()


class MplCanvas(FigureCanvas):
    def __init__(self, title: str, xlabel: str, ylabel: str):
        self.fig = Figure(figsize=(5, 3), dpi=100, tight_layout=True)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self._title = title
        self._xlabel = xlabel
        self._ylabel = ylabel
        self.style_ax()

    def style_ax(self):
        self.ax.set_title(self._title, color=THEME.TEXT, fontsize=11)
        self.ax.set_xlabel(self._xlabel, color=THEME.TEXT_DIM, fontsize=9)
        self.ax.set_ylabel(self._ylabel, color=THEME.TEXT_DIM, fontsize=9)
        self.ax.tick_params(colors=THEME.TEXT_DIM, labelsize=8)
        for sp in self.ax.spines.values():
            sp.set_edgecolor(THEME.BORDER)
        self.ax.set_facecolor(THEME.CANVAS)
        self.fig.patch.set_facecolor(THEME.CANVAS)

    def clear(self):
        self.ax.cla()
        self.style_ax()


class FitSettingsDialog(QDialog):
    def __init__(self, settings: FitSettings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fit Settings")
        self._settings = settings
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        form = QFormLayout()

        self.window_fraction = QDoubleSpinBox()
        self.window_fraction.setRange(0.05, 1.0)
        self.window_fraction.setSingleStep(0.05)
        self.window_fraction.setValue(self._settings.window_fraction)
        form.addRow("Window fraction", self.window_fraction)

        self.min_half_window = QDoubleSpinBox()
        self.min_half_window.setRange(5.0, 500.0)
        self.min_half_window.setSingleStep(5.0)
        self.min_half_window.setValue(self._settings.min_half_window)
        form.addRow("Min half-window (MeV)", self.min_half_window)

        self.max_half_window = QDoubleSpinBox()
        self.max_half_window.setRange(10.0, 1000.0)
        self.max_half_window.setSingleStep(10.0)
        self.max_half_window.setValue(self._settings.max_half_window)
        form.addRow("Max half-window (MeV)", self.max_half_window)

        self.min_entries = QDoubleSpinBox()
        self.min_entries.setDecimals(0)
        self.min_entries.setRange(1, 1_000_000)
        self.min_entries.setSingleStep(10)
        self.min_entries.setValue(self._settings.min_entries)
        form.addRow("Min entries", self.min_entries)

        lay.addLayout(form)
        row = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        row.addStretch()
        row.addWidget(ok_btn)
        row.addWidget(cancel_btn)
        lay.addLayout(row)

    def values(self) -> FitSettings:
        return FitSettings(
            window_fraction=float(self.window_fraction.value()),
            min_half_window=float(self.min_half_window.value()),
            max_half_window=float(self.max_half_window.value()),
            min_entries=int(self.min_entries.value()),
        )


class ThresholdDialog(QDialog):
    def __init__(self, settings: ThresholdSettings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Threshold Settings")
        self._settings = settings
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        form = QFormLayout()

        self.offset_pct = QDoubleSpinBox()
        self.offset_pct.setRange(0.1, 100.0)
        self.offset_pct.setSingleStep(0.5)
        self.offset_pct.setValue(self._settings.offset_pct)
        form.addRow("Offset threshold (%)", self.offset_pct)

        self.res_change_pct = QDoubleSpinBox()
        self.res_change_pct.setRange(0.1, 100.0)
        self.res_change_pct.setSingleStep(0.5)
        self.res_change_pct.setValue(self._settings.resolution_change_pct)
        form.addRow("Resolution change threshold (%)", self.res_change_pct)

        self.jump_pct = QDoubleSpinBox()
        self.jump_pct.setRange(0.1, 100.0)
        self.jump_pct.setSingleStep(0.5)
        self.jump_pct.setValue(self._settings.jump_pct)
        form.addRow("Jump threshold (%)", self.jump_pct)

        self.slope_abs_threshold = QDoubleSpinBox()
        self.slope_abs_threshold.setRange(0.001, 100.0)
        self.slope_abs_threshold.setSingleStep(0.1)
        self.slope_abs_threshold.setValue(self._settings.slope_abs_threshold)
        form.addRow("Slope threshold (MeV/index)", self.slope_abs_threshold)

        lay.addLayout(form)
        row = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        row.addStretch()
        row.addWidget(ok_btn)
        row.addWidget(cancel_btn)
        lay.addLayout(row)

    def values(self) -> ThresholdSettings:
        return ThresholdSettings(
            offset_pct=float(self.offset_pct.value()),
            resolution_change_pct=float(self.res_change_pct.value()),
            jump_pct=float(self.jump_pct.value()),
            slope_abs_threshold=float(self.slope_abs_threshold.value()),
        )


class MapRangeDialog(QDialog):
    def __init__(self, vmin: float, vmax: float, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Map Range")
        self._auto = False
        lay = QVBoxLayout(self)
        form = QFormLayout()

        self.vmin = QDoubleSpinBox()
        self.vmin.setRange(-1e9, 1e9)
        self.vmin.setDecimals(6)
        self.vmin.setValue(float(vmin))
        form.addRow("Min", self.vmin)

        self.vmax = QDoubleSpinBox()
        self.vmax.setRange(-1e9, 1e9)
        self.vmax.setDecimals(6)
        self.vmax.setValue(float(vmax))
        form.addRow("Max", self.vmax)
        lay.addLayout(form)

        row = QHBoxLayout()
        auto_btn = QPushButton("Auto")
        apply_btn = QPushButton("Apply")
        cancel_btn = QPushButton("Cancel")
        auto_btn.clicked.connect(self._on_auto)
        apply_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        row.addWidget(auto_btn)
        row.addStretch(1)
        row.addWidget(apply_btn)
        row.addWidget(cancel_btn)
        lay.addLayout(row)

    def _on_auto(self):
        self._auto = True
        self.accept()

    def values(self) -> Tuple[bool, float, float]:
        return self._auto, float(self.vmin.value()), float(self.vmax.value())


def _find_database_file(rel_path: str) -> Optional[Path]:
    candidates = [
        SCRIPT_DIR.parent / "database" / rel_path,
        Path.cwd() / "database" / rel_path,
    ]
    for c in candidates:
        if c.is_file():
            return c.resolve()
    return None


def _settings_hash(settings: FitSettings) -> str:
    payload = json.dumps(asdict(settings), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _load_run_config(runinfo_path: Path, run: int) -> Tuple[float, float]:
    """Return (ebeam, hycal_z) for run by chained config semantics."""
    with open(runinfo_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    def deep_merge(dst: dict, src: dict):
        for k, v in src.items():
            if isinstance(v, dict) and isinstance(dst.get(k), dict):
                deep_merge(dst[k], v)
            else:
                dst[k] = v

    merged = {}
    defaults = data.get("defaults", {})
    if isinstance(defaults, dict):
        deep_merge(merged, defaults)

    entries = data.get("configurations", [])
    chosen = []
    for e in entries:
        if not isinstance(e, dict):
            continue
        fr = e.get("from_run", e.get("run_number", -1))
        try:
            fr_i = int(fr)
        except Exception:
            continue
        if fr_i <= run:
            chosen.append((fr_i, e))
    chosen.sort(key=lambda x: x[0])
    for _, e in chosen:
        deep_merge(merged, e)

    ebeam = float(merged.get("beam_energy", 0.0))
    hycal = merged.get("hycal", {})
    pos = hycal.get("position", [0.0, 0.0, 6225.0])
    hycal_z = float(pos[2]) if isinstance(pos, list) and len(pos) >= 3 else 6225.0
    return ebeam, hycal_z


def _scan_quick_check(folder: Path, runinfo_path: Path) -> List[RunMeta]:
    runs: List[RunMeta] = []
    for p in sorted(folder.glob("prad_*_quick_check.root")):
        m = RUN_FILE_RE.search(p.name)
        if not m:
            continue
        run = int(m.group(1))
        ebeam, hycal_z = _load_run_config(runinfo_path, run)
        runs.append(RunMeta(run=run, root_path=str(p.resolve()), root_mtime=p.stat().st_mtime,
                            ebeam=ebeam, hycal_z=hycal_z))
    runs.sort(key=lambda r: r.run)
    return runs


def _build_run_groups(runs: List[RunMeta], merge_n: int) -> List[RunGroup]:
    if merge_n <= 1:
        out = []
        for r in runs:
            out.append(RunGroup(
                key_run=r.run,
                label=str(r.run),
                members=[r],
                ebeam=r.ebeam,
                hycal_z=r.hycal_z,
                root_mtime=r.root_mtime,
            ))
        return out

    out: List[RunGroup] = []
    for i in range(0, len(runs), merge_n):
        grp = runs[i:i + merge_n]
        if not grp:
            continue
        key = grp[0].run
        label = f"{grp[0].run}-{grp[-1].run}" if len(grp) > 1 else str(grp[0].run)
        out.append(RunGroup(
            key_run=key,
            label=label,
            members=grp,
            ebeam=float(np.mean(np.array([g.ebeam for g in grp], dtype=float))),
            hycal_z=float(np.mean(np.array([g.hycal_z for g in grp], dtype=float))),
            root_mtime=max(g.root_mtime for g in grp),
        ))
    return out


def _energy_loss(theta_deg: float) -> float:
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sec = 1.0 / cos_t if cos_t > 0.01 else 100.0
    eloss = 0.500 * 1.6 * sec
    eloss += 0.120 * 1.6 * sec
    eloss += 0.100 * 2.0 * sec
    eloss += 0.480 * 1.8 * sec
    return eloss


def _expected_energy_ep(theta_deg: float, ebeam: float) -> float:
    m_proton = 938.2720813
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    e = ebeam * m_proton / (m_proton + ebeam * (1.0 - cos_t))
    return max(0.0, e - _energy_loss(theta_deg))


def _theta_from_module(x: float, y: float, hycal_z: float) -> float:
    r = math.sqrt(x * x + y * y)
    return math.degrees(math.atan2(r, max(1e-6, hycal_z)))


def _resolution_percent(peak_center_mev: float, sigma_mev: float) -> float:
    """Return stochastic-term resolution in percent.

    a(%) = (sigma / E) * sqrt(E_GeV) * 100
    where E and sigma are provided in MeV.
    """
    # Always use fitted peak center as E.
    e_mev = max(1e-6, float(peak_center_mev))
    e_gev = e_mev / 1000.0
    return float(sigma_mev) / e_mev * math.sqrt(max(1e-12, e_gev)) * 100.0


def _compute_center_4x4_modules(modules: List[object]) -> Set[str]:
    """Return module names inside the central signed 4x4 rows/columns.

    This selects the 4 x-columns and 4 y-rows closest to zero (signed
    coordinates), which corresponds to the center block containing the beam
    hole region.
    """
    if not modules:
        return set()

    ux = sorted({float(m.x) for m in modules}, key=lambda v: (abs(v), v))
    uy = sorted({float(m.y) for m in modules}, key=lambda v: (abs(v), v))
    if len(ux) < 4 or len(uy) < 4:
        return set()

    x_sel = sorted(ux[:4])
    y_sel = sorted(uy[:4])
    # Geometric coordinates are exact map values, but keep small tolerance.
    x_tol = 1e-6
    y_tol = 1e-6

    excluded: Set[str] = set()
    for m in modules:
        x = float(m.x)
        y = float(m.y)
        in_x = any(abs(x - xv) <= x_tol for xv in x_sel)
        in_y = any(abs(y - yv) <= y_tol for yv in y_sel)
        if in_x and in_y:
            excluded.add(str(m.name))
    return excluded


def _compute_outer_8layer_modules(modules: List[object]) -> Set[str]:
    """Return module names in outermost 8 geometry layers."""
    if not modules:
        return set()
    ux = sorted({abs(float(m.x)) for m in modules}, reverse=True)
    uy = sorted({abs(float(m.y)) for m in modules}, reverse=True)
    idx_x = {v: i for i, v in enumerate(ux)}
    idx_y = {v: i for i, v in enumerate(uy)}

    out: Set[str] = set()
    for m in modules:
        dx = idx_x.get(abs(float(m.x)), 9999)
        dy = idx_y.get(abs(float(m.y)), 9999)
        depth = min(dx, dy)
        if depth < 8:
            out.add(str(m.name))
    return out


def _rebin2_hist(counts: np.ndarray, edges: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = int(counts.size)
    if n < 2:
        return counts, edges
    m = (n // 2) * 2
    c2 = counts[:m].reshape(-1, 2).sum(axis=1)
    e2 = edges[:m + 1:2]
    if m < n:
        c2 = np.concatenate([c2, np.array([counts[-1]], dtype=float)])
        e2 = np.concatenate([e2, np.array([edges[-1]], dtype=float)])
    return np.asarray(c2, dtype=float), np.asarray(e2, dtype=float)


def _run_cache_file(cache_root: Path, run: int) -> Path:
    return cache_root / f"run_{run}.json"


def _normalize_module_name(name: str) -> str:
    s = str(name).strip().upper()
    if s.startswith("W") and s[1:].isdigit():
        return f"W{int(s[1:])}"
    return s


def _hist_key_candidates(module_name: str) -> List[str]:
    raw = str(module_name).strip()
    up = raw.upper()
    low = raw.lower()
    out = [
        f"module_energy/h_{raw}",
        f"module_energy/h_{up}",
        f"module_energy/h_{low}",
        f"module_energy/{raw}",
        f"module_energy/{up}",
        f"module_energy/{low}",
        f"h_{raw}",
        f"h_{up}",
        f"h_{low}",
    ]
    norm = _normalize_module_name(raw)
    if norm.startswith("W") and norm[1:].isdigit():
        n = int(norm[1:])
        for w in (3, 4, 5):
            zp = f"W{n:0{w}d}"
            out.extend([
                f"module_energy/h_{zp}",
                f"module_energy/{zp}",
                f"h_{zp}",
            ])
    return list(dict.fromkeys(out))


def _get_hist_key_map(root_path: str) -> Dict[str, str]:
    cached = _HIST_KEY_CACHE.get(root_path)
    if cached is not None:
        return cached
    out: Dict[str, str] = {}
    if not HAS_UPROOT:
        _HIST_KEY_CACHE[root_path] = out
        return out
    try:
        with uproot.open(root_path) as f:
            keys = f.keys(recursive=True, cycle=False)
        for key in keys:
            k = str(key)
            if "module_energy" not in k:
                continue
            leaf = k.rsplit("/", 1)[-1]
            if leaf.startswith("h_"):
                leaf = leaf[2:]
            norm = _normalize_module_name(leaf)
            out.setdefault(leaf, k)
            out.setdefault(leaf.upper(), k)
            out.setdefault(norm, k)
    except Exception:
        out = {}
    _HIST_KEY_CACHE[root_path] = out
    return out


def _atomic_write_json(path: Path, payload: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp_fit_", suffix=".json", dir=str(path.parent))
    os.close(fd)
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _gaussian(x: np.ndarray, amp: float, mu: float, sigma: float) -> np.ndarray:
    sigma = max(1e-9, float(sigma))
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _fit_gaussian_array(
    centers: np.ndarray,
    counts: np.ndarray,
    fit_xmin: float,
    fit_xmax: float,
) -> Tuple[bool, Optional[float], Optional[float], Optional[float]]:
    """Fit gaussian on histogram arrays; returns (ok, mean, sigma, chi2_ndf)."""
    mask = (centers >= fit_xmin) & (centers <= fit_xmax) & (counts >= 0)
    x = centers[mask]
    y = counts[mask]
    if x.size < 4 or float(np.sum(y)) <= 0:
        return False, None, None, None

    amp0 = float(np.max(y))
    mu0 = float(x[int(np.argmax(y))])
    span = max(1e-3, fit_xmax - fit_xmin)
    sigma0 = max(1.0, span / 6.0)

    # Preferred: scipy non-linear least squares.
    if HAS_SCIPY:
        try:
            popt, _ = curve_fit(
                _gaussian,
                x,
                y,
                p0=[amp0, mu0, sigma0],
                bounds=([
                    0.0,
                    fit_xmin,
                    1e-6,
                ], [
                    np.inf,
                    fit_xmax,
                    span,
                ]),
                maxfev=10000,
            )
            amp, mu, sigma = float(popt[0]), float(popt[1]), abs(float(popt[2]))
            y_fit = _gaussian(x, amp, mu, sigma)
            denom = np.clip(y_fit, 1.0, None)
            ndf = int(x.size - 3)
            chi2_ndf = float(np.sum((y - y_fit) ** 2 / denom) / ndf) if ndf > 0 else 0.0
            return True, mu, sigma, chi2_ndf
        except Exception:
            pass

    # Fallback: weighted moments (works without scipy).
    w = np.clip(y, 0.0, None)
    sw = float(np.sum(w))
    if sw <= 0:
        return False, None, None, None
    mu = float(np.sum(w * x) / sw)
    var = float(np.sum(w * (x - mu) ** 2) / sw)
    sigma = math.sqrt(max(1e-9, var))
    y_fit = _gaussian(x, amp0, mu, sigma)
    denom = np.clip(y_fit, 1.0, None)
    ndf = int(x.size - 3)
    chi2_ndf = float(np.sum((y - y_fit) ** 2 / denom) / ndf) if ndf > 0 else 0.0
    return True, mu, sigma, chi2_ndf


def _fit_hist_for_module(root_path: str, module_name: str,
                         settings: FitSettings,
                         manual_window: Optional[Tuple[float, float]] = None,
                         rebin2: bool = False) -> dict:
    out = {
        "status": "unknown",
        "entries": 0,
        "mean": None,
        "sigma": None,
        "chi2_ndf": None,
        "fit_xmin": None,
        "fit_xmax": None,
    }
    arrays = _load_hist_arrays(root_path, module_name, rebin2=rebin2)
    if arrays is None:
        out["status"] = "missing_hist"
        return out

    counts, edges = arrays
    if counts.size == 0 or edges.size < 2:
        out["status"] = "missing_hist"
        return out

    centers = 0.5 * (edges[:-1] + edges[1:])
    entries = int(np.sum(np.clip(counts, 0.0, None)))
    out["entries"] = entries
    if entries < settings.min_entries:
        out["status"] = "low_stats"
        return out

    peak = float(centers[int(np.argmax(counts))])
    xmin_axis = float(edges[0])
    xmax_axis = float(edges[-1])

    if manual_window is not None:
        fit_xmin = max(xmin_axis, min(manual_window))
        fit_xmax = min(xmax_axis, max(manual_window))
        if fit_xmax <= fit_xmin:
            out["status"] = "bad_manual_window"
            return out
    else:
        half = abs(peak) * settings.window_fraction
        half = max(settings.min_half_window, min(settings.max_half_window, half))
        fit_xmin = max(xmin_axis, peak - half)
        fit_xmax = min(xmax_axis, peak + half)

    ok, mean, sigma, chi2 = _fit_gaussian_array(centers, counts, fit_xmin, fit_xmax)
    if not ok or mean is None or sigma is None:
        out["status"] = "fit_failed"
        out["fit_xmin"] = fit_xmin
        out["fit_xmax"] = fit_xmax
        return out

    out.update({
        "status": "ok",
        "mean": float(mean),
        "sigma": float(sigma),
        "chi2_ndf": float(chi2) if chi2 is not None else None,
        "fit_xmin": fit_xmin,
        "fit_xmax": fit_xmax,
    })
    return out


def _fit_hist_for_module_multi(root_paths: List[str], module_name: str,
                               settings: FitSettings,
                               manual_window: Optional[Tuple[float, float]] = None,
                               rebin2: bool = False) -> dict:
    out = {
        "status": "unknown",
        "entries": 0,
        "mean": None,
        "sigma": None,
        "chi2_ndf": None,
        "fit_xmin": None,
        "fit_xmax": None,
    }
    arrays = _load_hist_arrays_multi(root_paths, module_name, rebin2=rebin2)
    if arrays is None:
        out["status"] = "missing_hist"
        return out

    counts, edges = arrays
    if counts.size == 0 or edges.size < 2:
        out["status"] = "missing_hist"
        return out

    centers = 0.5 * (edges[:-1] + edges[1:])
    entries = int(np.sum(np.clip(counts, 0.0, None)))
    out["entries"] = entries
    if entries < settings.min_entries:
        out["status"] = "low_stats"
        return out

    peak = float(centers[int(np.argmax(counts))])
    xmin_axis = float(edges[0])
    xmax_axis = float(edges[-1])

    if manual_window is not None:
        fit_xmin = max(xmin_axis, min(manual_window))
        fit_xmax = min(xmax_axis, max(manual_window))
        if fit_xmax <= fit_xmin:
            out["status"] = "bad_manual_window"
            return out
    else:
        half = abs(peak) * settings.window_fraction
        half = max(settings.min_half_window, min(settings.max_half_window, half))
        fit_xmin = max(xmin_axis, peak - half)
        fit_xmax = min(xmax_axis, peak + half)

    ok, mean, sigma, chi2 = _fit_gaussian_array(centers, counts, fit_xmin, fit_xmax)
    if not ok or mean is None or sigma is None:
        out["status"] = "fit_failed"
        out["fit_xmin"] = fit_xmin
        out["fit_xmax"] = fit_xmax
        return out

    out.update({
        "status": "ok",
        "mean": float(mean),
        "sigma": float(sigma),
        "chi2_ndf": float(chi2) if chi2 is not None else None,
        "fit_xmin": fit_xmin,
        "fit_xmax": fit_xmax,
    })
    return out


def _load_hist_arrays(root_path: str, module_name: str, rebin2: bool = False) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    if not HAS_UPROOT:
        return None
    try:
        with uproot.open(root_path) as f:
            obj = None
            for key in _hist_key_candidates(module_name):
                obj = f.get(key)
                if obj is not None:
                    break
            if obj is None:
                key_map = _get_hist_key_map(root_path)
                resolved = key_map.get(_normalize_module_name(module_name))
                if resolved:
                    obj = f.get(resolved)
            if obj is None:
                return None
            counts, edges = obj.to_numpy(flow=False)
            counts = np.asarray(counts, dtype=float)
            edges = np.asarray(edges, dtype=float)
            if rebin2:
                counts, edges = _rebin2_hist(counts, edges)
            return counts, edges
    except Exception:
        return None


def _load_hist_arrays_multi(root_paths: List[str], module_name: str, rebin2: bool = False) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    ref_edges = None
    sum_counts = None
    for rp in root_paths:
        arr = _load_hist_arrays(rp, module_name, rebin2=rebin2)
        if arr is None:
            continue
        counts, edges = arr
        if ref_edges is None:
            ref_edges = edges
            sum_counts = counts.astype(float)
        else:
            if edges.shape != ref_edges.shape or not np.allclose(edges, ref_edges, atol=1e-9, rtol=1e-9):
                continue
            sum_counts = sum_counts + counts
    if ref_edges is None or sum_counts is None:
        return None
    return np.asarray(sum_counts, dtype=float), np.asarray(ref_edges, dtype=float)


class FitWorker(QThread):
    progress = pyqtSignal(int, int, str)
    runReady = pyqtSignal(int, dict)
    allDone = pyqtSignal(bool, str)

    def __init__(
        self,
        run_groups: List[RunGroup],
        module_xy: Dict[str, Tuple[float, float]],
        rebin_modules: Set[str],
        cache_root: Path,
        fit_settings: FitSettings,
        force_refit: bool,
        parent=None,
    ):
        super().__init__(parent)
        self.run_groups = run_groups
        self.module_xy = module_xy
        self.rebin_modules = rebin_modules
        self.cache_root = cache_root
        self.fit_settings = fit_settings
        self.force_refit = force_refit

    def _is_cache_valid(self, payload: dict, run_group: RunGroup) -> bool:
        if not isinstance(payload, dict):
            return False
        if int(payload.get("fit_cache_schema", 0)) != FIT_CACHE_SCHEMA_VERSION:
            return False
        if payload.get("settings_hash") != _settings_hash(self.fit_settings):
            return False
        root_mtime = float(payload.get("root_mtime", -1))
        if abs(root_mtime - run_group.root_mtime) > 1e-6:
            return False
        if payload.get("source_runs") != [m.run for m in run_group.members]:
            return False
        return isinstance(payload.get("modules"), dict)

    def run(self):
        if not HAS_UPROOT:
            self.allDone.emit(False, "uproot is not available")
            return
        try:
            total = len(self.run_groups)
            for idx, run_group in enumerate(self.run_groups, start=1):
                self.progress.emit(idx, total, f"Run {run_group.label}")
                cpath = _run_cache_file(self.cache_root, run_group.key_run)
                run_payload = None
                if (not self.force_refit) and cpath.is_file():
                    try:
                        with open(cpath, "r", encoding="utf-8") as f:
                            loaded = json.load(f)
                        if self._is_cache_valid(loaded, run_group):
                            run_payload = loaded
                    except Exception:
                        run_payload = None

                if run_payload is None:
                    modules = {}
                    root_paths = [m.root_path for m in run_group.members]
                    for mname, (mx, my) in self.module_xy.items():
                        one = _fit_hist_for_module_multi(
                            root_paths,
                            mname,
                            self.fit_settings,
                            rebin2=(mname in self.rebin_modules),
                        )
                        if one.get("status") == "ok":
                            mean = float(one["mean"])
                            sigma = float(one["sigma"])
                            theta = _theta_from_module(mx, my, run_group.hycal_z)
                            expected = _expected_energy_ep(theta, run_group.ebeam)
                            delta = mean - expected
                            resolution = _resolution_percent(mean, sigma)
                            one.update({
                                "theta_deg": theta,
                                "expected_energy": expected,
                                "delta_energy": delta,
                                "resolution": resolution,
                            })
                        modules[mname] = one
                    run_payload = {
                        "run": run_group.key_run,
                        "label": run_group.label,
                        "source_runs": [m.run for m in run_group.members],
                        "root_mtime": run_group.root_mtime,
                        "ebeam": run_group.ebeam,
                        "hycal_z": run_group.hycal_z,
                        "fit_cache_schema": FIT_CACHE_SCHEMA_VERSION,
                        "settings_hash": _settings_hash(self.fit_settings),
                        "updated_at": datetime.now().isoformat(timespec="seconds"),
                        "modules": modules,
                    }
                    _atomic_write_json(cpath, run_payload)

                self.runReady.emit(run_group.key_run, run_payload)
            self.allDone.emit(True, "Done")
        except Exception as exc:
            self.allDone.emit(False, f"Fit worker failed: {exc}\n{traceback.format_exc()}")


class MainWindow(QMainWindow):
    METRICS = [
        ("delta_energy", "DeltaE (MeV)"),
        ("chi2_ndf", "Chi2/NDF"),
        ("resolution", "Resolution (%)"),
        ("resolution_var", "Resolution variance"),
        ("entries", "Statistics"),
        ("offset_count", "Offset count"),
        ("res_change_count", "Resolution change count"),
        ("center_slope", "Center slope"),
        ("jump_count", "Jump count"),
    ]

    def __init__(self, initial_folder: Optional[Path] = None):
        super().__init__()
        self.setWindowTitle("HyCal Recon Quality Viewer")
        self.resize(1850, 1050)
        apply_theme_palette(self)

        self.modules_json = _find_database_file("hycal_map.json")
        self.runinfo_json = _find_database_file("runinfo/general.json")

        self.fit_settings = FitSettings()
        self.thresholds = ThresholdSettings()

        self.quick_folder: Optional[Path] = initial_folder
        self.cache_root: Optional[Path] = None

        self.w_modules = []
        self.module_xy: Dict[str, Tuple[float, float]] = {}
        self.fit_module_xy: Dict[str, Tuple[float, float]] = {}
        self.excluded_modules: Set[str] = set()
        self.rebin_outer_modules: Set[str] = set()
        self.run_list: List[RunMeta] = []
        self.run_groups: List[RunGroup] = []

        self.run_payloads: Dict[int, dict] = {}
        self.trends: Dict[str, dict] = {}
        self.current_module: Optional[str] = None
        self.map_range_override_by_metric: Dict[str, Tuple[float, float]] = {}
        self.last_map_auto_range_by_metric: Dict[str, Tuple[float, float]] = {}

        self.worker: Optional[FitWorker] = None
        self.span_selector = None

        self._build_ui()
        self._load_geometry()
        self._bind_map_events()

        if self.quick_folder is not None and self.quick_folder.is_dir():
            self.folder_edit.setText(str(self.quick_folder))
            self.scan_and_fit(force_refit=False)

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        top = QGridLayout()
        row = 0

        title = QLabel("HyCal Recon Quality Viewer")
        title.setFont(QFont("Consolas", 14, QFont.Weight.Bold))
        title.setStyleSheet(f"color:{THEME.ACCENT};")
        top.addWidget(title, row, 0, 1, 4)

        row += 1
        top.addWidget(QLabel("Quick-check folder:"), row, 0)
        self.folder_edit = QLineEdit()
        top.addWidget(self.folder_edit, row, 1, 1, 2)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.on_browse_folder)
        top.addWidget(browse_btn, row, 3)

        row += 1
        top.addWidget(QLabel("Cache folder:"), row, 0)
        self.cache_edit = QLineEdit()
        self.cache_edit.setReadOnly(True)
        top.addWidget(self.cache_edit, row, 1, 1, 2)

        self.scan_btn = QPushButton("Load / Fit")
        self.scan_btn.clicked.connect(lambda: self.scan_and_fit(force_refit=False))
        top.addWidget(self.scan_btn, row, 3)

        row += 1
        top.addWidget(QLabel("Run:"), row, 0)
        self.run_combo = QComboBox()
        self.run_combo.currentIndexChanged.connect(self.refresh_all_views)
        top.addWidget(self.run_combo, row, 1)

        top.addWidget(QLabel("Map metric:"), row, 2)
        self.metric_combo = QComboBox()
        for key, label in self.METRICS:
            self.metric_combo.addItem(label, key)
        self.metric_combo.currentIndexChanged.connect(self.refresh_map)
        top.addWidget(self.metric_combo, row, 3)
        self.map_range_btn = QPushButton("Map Range")
        self.map_range_btn.clicked.connect(self.on_map_range)
        top.addWidget(self.map_range_btn, row, 4)

        row += 1
        self.refit_btn = QPushButton("Refit All Runs")
        self.refit_btn.clicked.connect(lambda: self.scan_and_fit(force_refit=True))
        top.addWidget(self.refit_btn, row, 0)

        settings_btn = QPushButton("Fit Settings")
        settings_btn.clicked.connect(self.on_fit_settings)
        top.addWidget(settings_btn, row, 1)

        threshold_btn = QPushButton("Thresholds")
        threshold_btn.clicked.connect(self.on_threshold_settings)
        top.addWidget(threshold_btn, row, 2)

        self.status_label = QLabel("Ready")
        top.addWidget(self.status_label, row, 3)

        top.addWidget(QLabel("Merge adjacent runs:"), row, 4)
        self.merge_spin = QSpinBox()
        self.merge_spin.setRange(1, 20)
        self.merge_spin.setValue(1)
        top.addWidget(self.merge_spin, row, 5)

        row += 1
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        top.addWidget(self.progress, row, 0, 1, 6)

        root.addLayout(top)

        split = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(split, stretch=1)

        # Left panel
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        self.map_widget = QualityMapWidget()
        left_lay.addWidget(self.map_widget)
        split.addWidget(left)

        # Right panel
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)

        self.module_label = QLabel("Click a W module on map")
        self.module_label.setFont(QFont("Consolas", 12, QFont.Weight.Bold))
        self.module_label.setStyleSheet(f"color:{THEME.ACCENT};")
        right_lay.addWidget(self.module_label)

        self.center_canvas = MplCanvas("Peak center trend", "Run index delta", "Center (MeV)")
        right_lay.addWidget(self.center_canvas, stretch=1)

        self.res_canvas = MplCanvas("Resolution trend", "Run index delta", "Resolution (%)")
        right_lay.addWidget(self.res_canvas, stretch=1)

        spec_head = QHBoxLayout()
        spec_head.addWidget(QLabel("Spectrum run:"))
        self.prev_run_btn = QPushButton("Prev")
        self.prev_run_btn.setFixedWidth(64)
        self.prev_run_btn.clicked.connect(self.on_prev_run)
        spec_head.addWidget(self.prev_run_btn)
        self.next_run_btn = QPushButton("Next")
        self.next_run_btn.setFixedWidth(64)
        self.next_run_btn.clicked.connect(self.on_next_run)
        spec_head.addWidget(self.next_run_btn)
        spec_head.addStretch(1)
        right_lay.addLayout(spec_head)

        self.spec_canvas = MplCanvas("Energy spectrum", "Energy (MeV)", "Counts")
        right_lay.addWidget(self.spec_canvas, stretch=1)

        self.detail_label = QLabel("")
        self.detail_label.setWordWrap(True)
        self.detail_label.setStyleSheet(f"color:{THEME.TEXT_DIM};")
        right_lay.addWidget(self.detail_label)

        split.addWidget(right)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)

    def _load_geometry(self):
        if self.modules_json is None:
            QMessageBox.critical(self, "Error", "Cannot find database/hycal_map.json")
            return
        mods = load_modules(self.modules_json)
        self.w_modules = [m for m in mods if m.mod_type == "PbWO4"]
        self.module_xy = {m.name: (float(m.x), float(m.y)) for m in self.w_modules}
        self.excluded_modules = _compute_center_4x4_modules(self.w_modules)
        self.rebin_outer_modules = _compute_outer_8layer_modules(self.w_modules)
        self.fit_module_xy = {
            name: xy for name, xy in self.module_xy.items()
            if name not in self.excluded_modules
        }
        self.map_widget.set_modules(self.w_modules)
        self.status_label.setText(
            f"Excluded center 4x4: {len(self.excluded_modules)}; rebin2 outer8: {len(self.rebin_outer_modules)}"
        )

    def _bind_map_events(self):
        self.map_widget.moduleClicked.connect(self.on_module_clicked)

    def on_browse_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose quick-check folder", str(Path.cwd()))
        if not folder:
            return
        self.folder_edit.setText(folder)
        self.scan_and_fit(force_refit=False)

    def on_fit_settings(self):
        dlg = FitSettingsDialog(self.fit_settings, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        self.fit_settings = dlg.values()
        self.status_label.setText("Fit settings updated. Click Load/Fit or Refit.")

    def on_threshold_settings(self):
        dlg = ThresholdDialog(self.thresholds, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        self.thresholds = dlg.values()
        self.compute_trends()
        self.refresh_all_views()

    def on_map_range(self):
        metric = str(self.metric_combo.currentData())
        override = self.map_range_override_by_metric.get(metric)
        auto = self.last_map_auto_range_by_metric.get(metric, (0.0, 1.0))
        vmin, vmax = override if override is not None else auto
        dlg = MapRangeDialog(vmin, vmax, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        is_auto, lo, hi = dlg.values()
        if is_auto:
            self.map_range_override_by_metric.pop(metric, None)
        else:
            if hi < lo:
                lo, hi = hi, lo
            if math.isclose(lo, hi):
                hi = lo + 1.0
            self.map_range_override_by_metric[metric] = (lo, hi)
        self.refresh_map()

    def scan_and_fit(self, force_refit: bool):
        if not HAS_UPROOT:
            QMessageBox.critical(
                self,
                "uproot missing",
                "Need uproot to read quick_check ROOT files.\n"
                "Install with: python3 -m pip install --user uproot",
            )
            return
        if self.runinfo_json is None:
            QMessageBox.critical(self, "Error", "Cannot find database/runinfo/general.json")
            return

        folder_text = self.folder_edit.text().strip()
        if not folder_text:
            QMessageBox.warning(self, "Input needed", "Please select quick-check folder")
            return
        qdir = Path(folder_text)
        if not qdir.is_dir():
            QMessageBox.warning(self, "Invalid path", f"Not a folder: {qdir}")
            return

        self.quick_folder = qdir
        self.cache_root = qdir / "hycal_recon_fit_cache"
        self.cache_edit.setText(str(self.cache_root))

        try:
            self.run_list = _scan_quick_check(qdir, self.runinfo_json)
        except Exception as exc:
            QMessageBox.critical(self, "Scan failed", str(exc))
            return

        if not self.run_list:
            QMessageBox.warning(self, "No files", "No prad_*_quick_check.root found")
            return

        self.run_payloads.clear()
        merge_n = int(self.merge_spin.value())
        self.run_groups = _build_run_groups(self.run_list, merge_n)
        self.current_module = None
        self.module_label.setText("Click a W module on map")

        self._set_busy(True)
        self.progress.setValue(0)
        self.status_label.setText("Fitting...")

        self.worker = FitWorker(
            run_groups=self.run_groups,
            module_xy=self.fit_module_xy,
            rebin_modules=self.rebin_outer_modules,
            cache_root=self.cache_root,
            fit_settings=self.fit_settings,
            force_refit=force_refit,
            parent=self,
        )
        self.worker.progress.connect(self.on_fit_progress)
        self.worker.runReady.connect(self.on_run_ready)
        self.worker.allDone.connect(self.on_fit_done)
        self.worker.start()

    def _set_busy(self, busy: bool):
        self.scan_btn.setEnabled(not busy)
        self.refit_btn.setEnabled(not busy)
        self.metric_combo.setEnabled(not busy)
        self.run_combo.setEnabled(not busy)

    def on_fit_progress(self, idx: int, total: int, msg: str):
        self.status_label.setText(msg)
        if total > 0:
            self.progress.setValue(int(100 * idx / total))

    def on_run_ready(self, run: int, payload: dict):
        self.run_payloads[run] = payload

    def on_fit_done(self, ok: bool, message: str):
        self._set_busy(False)
        if not ok:
            QMessageBox.critical(self, "Fit failed", message)
            self.status_label.setText("Fit failed")
            return
        self.progress.setValue(100)
        self.status_label.setText("Loaded")

        self.run_combo.blockSignals(True)
        self.run_combo.clear()
        for rg in self.run_groups:
            self.run_combo.addItem(rg.label, rg.key_run)
        self.run_combo.blockSignals(False)

        self.compute_trends()
        self.refresh_all_views()

    def on_prev_run(self):
        idx = self.run_combo.currentIndex()
        if idx > 0:
            self.run_combo.setCurrentIndex(idx - 1)

    def on_next_run(self):
        idx = self.run_combo.currentIndex()
        if idx >= 0 and idx + 1 < self.run_combo.count():
            self.run_combo.setCurrentIndex(idx + 1)

    def compute_trends(self):
        self.trends = {}
        if not self.run_list:
            return

        runs = sorted(self.run_payloads.keys())
        for mname in self.fit_module_xy:
            centers = []
            resolutions = []
            valid_runs = []
            for r in runs:
                mod = self.run_payloads.get(r, {}).get("modules", {}).get(mname, {})
                if mod.get("status") != "ok":
                    continue
                centers.append(float(mod["mean"]))
                resolutions.append(float(mod["resolution"]))
                valid_runs.append(r)

            if not valid_runs:
                self.trends[mname] = {
                    "offset_count": 0,
                    "res_change_count": 0,
                    "center_slope": 0.0,
                    "jump_count": 0,
                }
                continue

            offset_count = 0
            for r in valid_runs:
                mod = self.run_payloads[r]["modules"][mname]
                exp = float(mod.get("expected_energy", 0.0))
                delta = float(mod.get("delta_energy", 0.0))
                if exp > 0 and abs(delta) / exp * 100.0 > self.thresholds.offset_pct:
                    offset_count += 1

            res_change_count = 0
            jump_count = 0
            for i in range(1, len(centers)):
                prev_res = resolutions[i - 1]
                cur_res = resolutions[i]
                if prev_res > 0 and abs(cur_res - prev_res) / prev_res * 100.0 > self.thresholds.resolution_change_pct:
                    res_change_count += 1

                prev_c = centers[i - 1]
                cur_c = centers[i]
                if prev_c != 0 and abs(cur_c - prev_c) / abs(prev_c) * 100.0 > self.thresholds.jump_pct:
                    jump_count += 1

            if len(centers) >= 2:
                x = np.arange(len(centers), dtype=float)
                slope = float(np.polyfit(x, np.array(centers, dtype=float), 1)[0])
            else:
                slope = 0.0

            self.trends[mname] = {
                "offset_count": offset_count,
                "res_change_count": res_change_count,
                "center_slope": slope,
                "jump_count": jump_count,
            }

    def _selected_run(self) -> Optional[int]:
        if self.run_combo.currentIndex() < 0:
            return None
        return int(self.run_combo.currentData())

    def refresh_all_views(self):
        self.refresh_map()
        self.refresh_module_views()

    def refresh_map(self):
        run = self._selected_run()
        if run is None:
            self.map_widget.set_values({})
            return

        metric = self.metric_combo.currentData()
        label = self.metric_combo.currentText()
        self.map_widget.set_metric_label(label)

        vals: Dict[str, float] = {}
        if metric in {"delta_energy", "chi2_ndf", "resolution", "entries"}:
            payload = self.run_payloads.get(run, {})
            mod_map = payload.get("modules", {})
            for mname in self.fit_module_xy:
                mod = mod_map.get(mname, {})
                if mod.get("status") != "ok":
                    continue
                v = mod.get(metric)
                if v is not None:
                    vals[mname] = float(v)
        elif metric == "resolution_var":
            payload = self.run_payloads.get(run, {})
            mod_map = payload.get("modules", {})
            res_vals: Dict[str, float] = {}
            for mname in self.fit_module_xy:
                mod = mod_map.get(mname, {})
                if mod.get("status") != "ok":
                    continue
                rv = mod.get("resolution")
                if rv is not None:
                    res_vals[mname] = float(rv)
            if res_vals:
                avg_res = float(np.mean(np.array(list(res_vals.values()), dtype=float)))
                vals = {k: (v - avg_res) ** 2 for k, v in res_vals.items()}
        else:
            for mname in self.fit_module_xy:
                vals[mname] = float(self.trends.get(mname, {}).get(metric, 0.0))

        self.map_widget.set_values(vals)
        override = self.map_range_override_by_metric.get(str(metric))
        if vals:
            vmin = min(vals.values())
            vmax = max(vals.values())
            if math.isclose(vmin, vmax):
                vmax = vmin + 1.0
            self.last_map_auto_range_by_metric[str(metric)] = (vmin, vmax)
            if override is None:
                self.map_widget.set_range(vmin, vmax)
            else:
                self.map_widget.set_range(override[0], override[1])
        else:
            self.last_map_auto_range_by_metric[str(metric)] = (0.0, 1.0)
            if override is None:
                self.map_widget.set_range(0.0, 1.0)
            else:
                self.map_widget.set_range(override[0], override[1])

    def _set_adaptive_ylim(self, ax, arrays: List[np.ndarray], min_span: float):
        vals = []
        for arr in arrays:
            if arr is None:
                continue
            a = np.asarray(arr, dtype=float)
            if a.size == 0:
                continue
            a = a[np.isfinite(a)]
            if a.size:
                vals.append(a)
        if not vals:
            return
        allv = np.concatenate(vals)
        lo = float(np.min(allv))
        hi = float(np.max(allv))
        if not np.isfinite(lo) or not np.isfinite(hi):
            return
        span = hi - lo
        if span < min_span:
            c = 0.5 * (lo + hi)
            lo = c - 0.5 * min_span
            hi = c + 0.5 * min_span
        else:
            pad = 0.12 * span
            lo -= pad
            hi += pad
        ax.set_ylim(lo, hi)

    def on_module_clicked(self, name: str):
        self.current_module = name if name else None
        self.refresh_module_views()

    def refresh_module_views(self):
        run = self._selected_run()
        mname = self.current_module
        if run is None or not mname:
            self.module_label.setText("Click a W module on map")
            self.center_canvas.clear()
            self.center_canvas.draw_idle()
            self.res_canvas.clear()
            self.res_canvas.draw_idle()
            self.spec_canvas.clear()
            self.spec_canvas.draw_idle()
            self.detail_label.setText("")
            self._reset_span_selector()
            return

        if mname in self.excluded_modules:
            self.module_label.setText(f"Module: {mname} (excluded center 4x4)")
            self.center_canvas.clear()
            self.center_canvas.draw_idle()
            self.res_canvas.clear()
            self.res_canvas.draw_idle()
            self.spec_canvas.clear()
            self.spec_canvas.draw_idle()
            self.detail_label.setText("This module is excluded from fit and trend calculations.")
            self._reset_span_selector()
            return

        self.module_label.setText(f"Module: {mname}   Run: {run}")

        valid = []
        for r in sorted(self.run_payloads.keys()):
            mod = self.run_payloads.get(r, {}).get("modules", {}).get(mname, {})
            if mod.get("status") != "ok":
                continue
            valid.append((r, mod))

        self._draw_center_trend(valid)
        self._draw_resolution_trend(valid)
        self._draw_spectrum(run, mname)

        current_mod = self.run_payloads.get(run, {}).get("modules", {}).get(mname, {})
        status = current_mod.get("status", "unknown")
        mean = current_mod.get("mean")
        sigma = current_mod.get("sigma")
        chi2 = current_mod.get("chi2_ndf")
        exp = current_mod.get("expected_energy")
        de = current_mod.get("delta_energy")
        res = current_mod.get("resolution")
        self.detail_label.setText(
            f"status={status}  mean={_fmt(mean)}  sigma={_fmt(sigma)}  "
            f"chi2/ndf={_fmt(chi2)}  expected={_fmt(exp)}  "
            f"delta={_fmt(de)}  resolution={_fmt(res)}%"
        )

    def _draw_center_trend(self, valid: List[Tuple[int, dict]]):
        self.center_canvas.clear()
        if not valid:
            self.center_canvas.draw_idle()
            return
        first_run = valid[0][0]
        x = np.array([r - first_run for r, _ in valid], dtype=float)
        y = np.array([float(m["mean"]) for _, m in valid], dtype=float)
        yerr = np.array([float(m.get("sigma", 0.0)) for _, m in valid], dtype=float)
        exp = np.array([float(m.get("expected_energy", np.nan)) for _, m in valid], dtype=float)
        self.center_canvas.ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o-",
            lw=1.2,
            ms=4.0,
            elinewidth=0.9,
            capsize=2.0,
            color="#58a6ff",
        )
        if np.isfinite(exp).any():
            off = self.thresholds.offset_pct / 100.0
            self.center_canvas.ax.plot(x, exp, lw=1.0, ls="--", color="#9aa4af")
            self.center_canvas.ax.plot(x, exp * (1.0 + off), lw=1.0, ls="--", color="#f97316")
            self.center_canvas.ax.plot(x, exp * (1.0 - off), lw=1.0, ls="--", color="#f97316")
            self._set_adaptive_ylim(
                self.center_canvas.ax,
                [y + yerr, y - yerr, exp, exp * (1.0 + off), exp * (1.0 - off)],
                min_span=40.0,
            )
        else:
            self._set_adaptive_ylim(self.center_canvas.ax, [y + yerr, y - yerr], min_span=40.0)
        self.center_canvas.ax.set_title(f"Peak center trend (first run={first_run})", color=THEME.TEXT)
        self.center_canvas.ax.grid(alpha=0.25)
        self.center_canvas.draw_idle()

    def _draw_resolution_trend(self, valid: List[Tuple[int, dict]]):
        self.res_canvas.clear()
        if not valid:
            self.res_canvas.draw_idle()
            return
        first_run = valid[0][0]
        x = np.array([r - first_run for r, _ in valid], dtype=float)
        y = np.array([float(m["resolution"]) for _, m in valid], dtype=float)
        self.res_canvas.ax.plot(x, y, marker="o", lw=1.4, color="#ff9f0a")
        nominal = float(np.mean(y)) if y.size else 0.0
        band = self.thresholds.resolution_change_pct / 100.0
        self.res_canvas.ax.axhline(nominal, color="#9aa4af", ls="--", lw=1.0)
        self.res_canvas.ax.axhline(nominal * (1.0 + band), color="#f97316", ls="--", lw=1.0)
        self.res_canvas.ax.axhline(nominal * (1.0 - band), color="#f97316", ls="--", lw=1.0)
        self._set_adaptive_ylim(
            self.res_canvas.ax,
            [
                y,
                np.array([nominal]),
                np.array([nominal * (1.0 + band)]),
                np.array([nominal * (1.0 - band)]),
            ],
            min_span=0.3,
        )
        self.res_canvas.ax.set_title(f"Resolution trend (first run={first_run})", color=THEME.TEXT)
        self.res_canvas.ax.grid(alpha=0.25)
        self.res_canvas.draw_idle()

    def _draw_spectrum(self, run: int, mname: str):
        self.spec_canvas.clear()
        group = next((g for g in self.run_groups if g.key_run == run), None)
        if group is None:
            self.spec_canvas.draw_idle()
            self._reset_span_selector()
            return

        arrays = _load_hist_arrays_multi([m.root_path for m in group.members], mname, rebin2=(mname in self.rebin_outer_modules))
        if arrays is None:
            self.spec_canvas.ax.text(0.5, 0.5, "No histogram", ha="center", va="center",
                                     transform=self.spec_canvas.ax.transAxes, color=THEME.TEXT_DIM)
            self.spec_canvas.draw_idle()
            self._reset_span_selector()
            return

        counts, edges = arrays
        centers = 0.5 * (edges[:-1] + edges[1:])
        self.spec_canvas.ax.step(centers, counts, where="mid", color="#cccccc", lw=1.0)

        mod = self.run_payloads.get(run, {}).get("modules", {}).get(mname, {})
        if mod.get("status") == "ok":
            mean = float(mod["mean"])
            sigma = float(mod["sigma"])
            expected = mod.get("expected_energy")
            amp = float(np.max(counts)) if counts.size else 0.0
            if sigma > 0 and amp > 0:
                x = np.linspace(np.min(centers), np.max(centers), 500)
                y = amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
                self.spec_canvas.ax.plot(x, y, color="#f97316", lw=1.6)
            if expected is not None:
                self.spec_canvas.ax.axvline(float(expected), color="#9aa4af", ls=":", lw=1.0)
            fx0 = mod.get("fit_xmin")
            fx1 = mod.get("fit_xmax")
            if fx0 is not None and fx1 is not None:
                self.spec_canvas.ax.axvline(float(fx0), color="#58a6ff", ls="--", lw=1.0)
                self.spec_canvas.ax.axvline(float(fx1), color="#58a6ff", ls="--", lw=1.0)

        self.spec_canvas.ax.grid(alpha=0.25)
        self.spec_canvas.ax.set_title(f"Spectrum: {mname} @ run {run} (drag to refit)", color=THEME.TEXT)
        self.spec_canvas.ax.set_xlim(1800.0, 2600.0)
        self.spec_canvas.draw_idle()

        self._reset_span_selector()
        self.span_selector = SpanSelector(
            self.spec_canvas.ax,
            self.on_manual_span,
            "horizontal",
            useblit=True,
            props=dict(alpha=0.2, facecolor="#58a6ff"),
            interactive=False,
            drag_from_anywhere=True,
        )

    def _reset_span_selector(self):
        if self.span_selector is not None:
            try:
                self.span_selector.set_active(False)
            except Exception:
                pass
        self.span_selector = None

    def on_manual_span(self, xmin: float, xmax: float):
        run = self._selected_run()
        mname = self.current_module
        if run is None or not mname:
            return
        if mname in self.excluded_modules:
            return
        if abs(xmax - xmin) < 1e-3:
            return

        group = next((g for g in self.run_groups if g.key_run == run), None)
        if group is None:
            return

        one = _fit_hist_for_module_multi(
            [m.root_path for m in group.members],
            mname,
            self.fit_settings,
            manual_window=(xmin, xmax),
            rebin2=(mname in self.rebin_outer_modules),
        )
        if one.get("status") == "ok":
            mx, my = self.module_xy[mname]
            theta = _theta_from_module(mx, my, group.hycal_z)
            expected = _expected_energy_ep(theta, group.ebeam)
            mean = float(one["mean"])
            sigma = float(one["sigma"])
            one.update({
                "theta_deg": theta,
                "expected_energy": expected,
                "delta_energy": mean - expected,
                "resolution": _resolution_percent(mean, sigma),
            })

        payload = self.run_payloads.get(run)
        if not payload:
            return
        payload.setdefault("modules", {})[mname] = one
        payload["fit_cache_schema"] = FIT_CACHE_SCHEMA_VERSION
        payload["updated_at"] = datetime.now().isoformat(timespec="seconds")

        if self.cache_root is not None:
            cpath = _run_cache_file(self.cache_root, run)
            _atomic_write_json(cpath, payload)

        self.compute_trends()
        self.refresh_all_views()
        self.status_label.setText(f"Manual refit updated: run {run} {mname}")


def _fmt(v) -> str:
    if v is None:
        return "NA"
    try:
        x = float(v)
    except Exception:
        return "NA"
    if math.isnan(x) or math.isinf(x):
        return "NA"
    if x == 0:
        return "0"
    return f"{x:.5g}"


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="HyCal recon quality viewer")
    ap.add_argument("folder", nargs="?", default="", help="quick_check folder")
    ap.add_argument("--theme", choices=available_themes(), default="dark")
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    set_theme(args.theme)

    app = QApplication(sys.argv)
    initial_folder = Path(args.folder).resolve() if args.folder else None
    w = MainWindow(initial_folder=initial_folder)
    w.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
