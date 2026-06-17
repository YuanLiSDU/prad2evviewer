#!/usr/bin/env python3
"""
EPICS Tree Viewer
=================

Pure PyQt6 viewer for the slow-control side trees in PRad-II ROOT files:
``epics``, ``scalers`` and ``runinfo``.  The implementation intentionally
follows the style of the scripts/ viewers: standard Qt widgets, QPainter
drawing, and a background loader thread.

Usage
-----
    python3 analysis/tools/epics_viewer.py input.root
    python3 analysis/tools/epics_viewer.py
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Workaround for mixed Qt installations.
#
# Some systems load PyQt6 from ~/.local but resolve parts of Qt (for example
# libQt6DBus.so.6) from /usr/lib.  If those Qt builds do not match, imports can
# fail with Qt_6_PRIVATE_API undefined-symbol errors.  scripts/replay_viewer.py
# uses the same idea: prefer PyQt6's bundled Qt6 libraries, then re-exec before
# any PyQt shared library is loaded.
# ---------------------------------------------------------------------------
def _fix_qt_lib_path() -> None:
    import site

    site_dirs: List[str] = []
    try:
        site_dirs += site.getsitepackages()
    except AttributeError:
        pass
    try:
        site_dirs.append(site.getusersitepackages())
    except AttributeError:
        pass

    for sp in site_dirs:
        qt6_lib = Path(sp) / "PyQt6" / "Qt6" / "lib"
        if not qt6_lib.is_dir():
            continue
        qt6_lib_s = str(qt6_lib)
        current = os.environ.get("LD_LIBRARY_PATH", "")
        if qt6_lib_s in current.split(":"):
            return
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = qt6_lib_s + (":" + current if current else "")
        os.execve(sys.executable, [sys.executable] + sys.argv, env)


if __name__ == "__main__":
    _fix_qt_lib_path()

from PyQt6.QtCore import QRectF, QThread, Qt, pyqtSignal
from PyQt6.QtGui import QAction, QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QCompleter,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


BG = "#0f1419"
PANEL = "#151b22"
PANEL_2 = "#111820"
BORDER = "#2f3a46"
GRID = "#293440"
TEXT = "#d8dee9"
TEXT_DIM = "#8b98a8"
ACCENT = "#53b6ff"
POINT = "#51d18f"
POINT_BAD = "#ff6f6f"

SCALER_SCALARS = [
    "event_number",
    "ti_ticks",
    "sync_counter",
    "run_number",
    "trigger_type",
    "slot",
    "gated",
    "ungated",
    "live_ratio",
    "source",
    "channel",
    "ref_gated",
    "ref_ungated",
]
SCALER_ARRAYS = ["trg_gated", "trg_ungated", "tdc_gated", "tdc_ungated"]
RUNINFO_SCALARS = ["run_number", "run_type", "event_tag"]


# ===========================================================================
#  Data loading
# ===========================================================================

class EpicsLoadError(RuntimeError):
    pass


class EpicsLoader(QThread):
    finished = pyqtSignal(object, str)

    def __init__(
        self,
        path: str,
        allow_uproot: bool = False,
        allow_pyroot: bool = False,
        parent=None,
    ):
        super().__init__(parent)
        self._path = path
        self._allow_uproot = allow_uproot
        self._allow_pyroot = allow_pyroot

    def run(self):
        try:
            self.finished.emit(
                load_epics_file(self._path, self._allow_uproot, self._allow_pyroot),
                "",
            )
        except Exception as exc:
            self.finished.emit({}, str(exc))


def load_epics_file(
    path: str,
    allow_uproot: bool = False,
    allow_pyroot: bool = False,
) -> Dict[str, Any]:
    root_path = Path(path).expanduser()
    if not root_path.is_file():
        raise EpicsLoadError(f"File does not exist: {root_path}")

    errors: List[str] = []
    try:
        return _load_with_root_cli(root_path)
    except ImportError as exc:
        errors.append(str(exc))
    except EpicsLoadError:
        raise
    except Exception as exc:
        errors.append(f"_load_with_root_cli: {type(exc).__name__}: {exc}")

    try:
        return _load_with_pyroot_subprocess(root_path)
    except ImportError as exc:
        errors.append(str(exc))
    except EpicsLoadError:
        raise
    except Exception as exc:
        errors.append(f"_load_with_pyroot_subprocess: {type(exc).__name__}: {exc}")

    if allow_uproot:
        try:
            return _load_with_uproot(root_path)
        except ImportError as exc:
            errors.append(str(exc))
        except Exception as exc:
            errors.append(f"_load_with_uproot: {type(exc).__name__}: {exc}")

    if allow_pyroot:
        try:
            return _load_with_pyroot(root_path)
        except ImportError as exc:
            errors.append(str(exc))
        except EpicsLoadError:
            raise
        except Exception as exc:
            errors.append(f"_load_with_pyroot: {type(exc).__name__}: {exc}")

    hint = (
        "Could not read the epics/scalers/runinfo trees. Please source a ROOT "
        "environment so either the 'root' command or PyROOT is available."
    )
    if not allow_uproot:
        hint += " You can also try --uproot if uproot and awkward are installed."
    if not allow_pyroot:
        hint += (
            " As a last resort, try --pyroot, but PyROOT is disabled by default "
            "because it can crash during Qt shutdown on some machines."
        )
    if errors:
        hint += "\n\nTried:\n" + "\n".join(f"  - {e}" for e in errors)
    raise EpicsLoadError(hint)


def _append_point(
    series: Dict[str, Dict[str, List[Any]]],
    channel: str,
    unix_time: float,
    value: float,
    good: Optional[bool],
) -> None:
    if not math.isfinite(unix_time) or unix_time <= 0 or not math.isfinite(value):
        return
    rec = series.setdefault(channel, {"x": [], "y": [], "good": []})
    rec["x"].append(float(unix_time))
    rec["y"].append(float(value))
    rec["good"].append(True if good is None else bool(good))


def _finalise_payload(
    path: Path,
    series: Dict[str, Dict[str, List[Any]]],
    n_entries: int,
    t0: Optional[float],
    t1: Optional[float],
    loader: str,
) -> Dict[str, Any]:
    if not series:
        raise EpicsLoadError(f"{path}: epics tree has no channel/value points")

    all_x = [
        float(x)
        for rec in series.values()
        for x in rec.get("x", [])
        if math.isfinite(float(x)) and float(x) > 0
    ]
    if not all_x:
        raise EpicsLoadError(f"{path}: no time-stamped points in epics/scalers/runinfo")

    gated = series.get("scalers:gated")
    ungated = series.get("scalers:ungated")
    evn = series.get("scalers:event_number")
    if gated and ungated and "DAQ live time [%]" not in series:
        rows = []
        n = min(len(gated.get("x", [])), len(gated.get("y", [])),
                len(ungated.get("y", [])))
        evns = evn.get("y", []) if evn else []
        good = gated.get("good", [])
        for i in range(n):
            try:
                order_key = float(evns[i]) if i < len(evns) else float(gated["x"][i])
                rows.append((
                    order_key,
                    float(gated["x"][i]),
                    int(float(gated["y"][i])),
                    int(float(ungated["y"][i])),
                    True if i >= len(good) else bool(good[i]),
                ))
            except (TypeError, ValueError):
                continue
        rows.sort(key=lambda r: r[0])

        xs: List[float] = []
        ys: List[float] = []
        good_vals: List[bool] = []
        prev_g = 0
        prev_u = 0
        for _order_key, x, g, u, is_good in rows:
            if g < prev_g or u < prev_u:
                prev_g = 0
                prev_u = 0
            dg = g - prev_g
            du = u - prev_u
            if du > 0 and 0 <= dg <= du:
                xs.append(float(x))
                ys.append(float(dg) / float(du) * 100.0)
                good_vals.append(is_good)
            prev_g = g
            prev_u = u
        if xs:
            series["DAQ live time [%]"] = {"x": xs, "y": ys, "good": good_vals}

    live_ratio = series.get("scalers:live_ratio")
    if live_ratio and "DAQ cumulative live time [%]" not in series:
        xs = []
        ys = []
        good_vals = []
        good = live_ratio.get("good", [])
        for i, (x, y) in enumerate(zip(live_ratio.get("x", []),
                                       live_ratio.get("y", []))):
            try:
                yf = float(y)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(yf) or yf < 0:
                continue
            xs.append(float(x))
            ys.append(100.0 * yf)
            good_vals.append(True if i >= len(good) else bool(good[i]))
        if xs:
            series["DAQ cumulative live time [%]"] = {
                "x": xs, "y": ys, "good": good_vals
            }

    t0_calc = min(all_x) if t0 is None else t0
    t1_calc = max(all_x) if t1 is None else t1
    for rec in series.values():
        rec["x"] = [round((float(x) - t0_calc) / 60.0, 6) for x in rec.get("x", [])]

    channels = sorted(series.keys(), key=str.casefold)
    n_points = sum(len(s["x"]) for s in series.values())
    duration_min = max(0.0, (t1_calc - t0_calc) / 60.0)

    return {
        "path": str(path),
        "title": path.name,
        "loader": loader,
        "n_entries": int(n_entries),
        "n_points": int(n_points),
        "duration_min": float(duration_min),
        "t0_unix": int(t0_calc),
        "t1_unix": int(t1_calc),
        "channels": channels,
        "series": series,
    }


def _load_with_pyroot(path: Path) -> Dict[str, Any]:
    try:
        import ROOT  # type: ignore
    except Exception as exc:
        raise ImportError("PyROOT is not available in this Python.") from exc

    ROOT.gROOT.SetBatch(True)
    fin = ROOT.TFile.Open(str(path), "READ")
    if not fin or fin.IsZombie():
        raise EpicsLoadError(f"Cannot open ROOT file: {path}")

    try:
        series: Dict[str, Dict[str, List[Any]]] = {}
        n_entries = 0
        t0: Optional[float] = None
        t1: Optional[float] = None

        def note_time(unix_time: float) -> None:
            nonlocal t0, t1
            if not math.isfinite(unix_time) or unix_time <= 0:
                return
            t0 = unix_time if t0 is None else min(t0, unix_time)
            t1 = unix_time if t1 is None else max(t1, unix_time)

        tree = fin.Get("epics")
        if tree:
            branches = {br.GetName() for br in tree.GetListOfBranches()}
            if {"unix_time", "channel", "value"} <= branches:
                has_good = "good" in branches
                n_epics = int(tree.GetEntries())
                n_entries += n_epics
                for i in range(n_epics):
                    tree.GetEntry(i)
                    unix_time = float(getattr(tree, "unix_time"))
                    note_time(unix_time)
                    good = bool(getattr(tree, "good")) if has_good else None
                    channels = getattr(tree, "channel")
                    values = getattr(tree, "value")
                    n = min(int(channels.size()), int(values.size()))
                    for j in range(n):
                        _append_point(
                            series, str(channels[j]), unix_time, float(values[j]), good
                        )

        tree = fin.Get("scalers")
        if tree:
            branches = {br.GetName() for br in tree.GetListOfBranches()}
            if "unix_time" in branches:
                has_good = "good" in branches
                n_scalers = int(tree.GetEntries())
                n_entries += n_scalers
                for i in range(n_scalers):
                    tree.GetEntry(i)
                    unix_time = float(getattr(tree, "unix_time"))
                    note_time(unix_time)
                    good = bool(getattr(tree, "good")) if has_good else None
                    for name in SCALER_SCALARS:
                        if name in branches:
                            _append_point(
                                series,
                                f"scalers:{name}",
                                unix_time,
                                float(getattr(tree, name)),
                                good,
                            )
                    for name in SCALER_ARRAYS:
                        if name not in branches:
                            continue
                        arr = getattr(tree, name)
                        for ch in range(16):
                            _append_point(
                                series,
                                f"scalers:{name}[{ch}]",
                                unix_time,
                                float(arr[ch]),
                                good,
                            )

        tree = fin.Get("runinfo")
        if tree:
            branches = {br.GetName() for br in tree.GetListOfBranches()}
            if "unix_time" in branches:
                n_runinfo = int(tree.GetEntries())
                n_entries += n_runinfo
                for i in range(n_runinfo):
                    tree.GetEntry(i)
                    unix_time = float(getattr(tree, "unix_time"))
                    note_time(unix_time)
                    for name in RUNINFO_SCALARS:
                        if name in branches:
                            _append_point(
                                series,
                                f"runinfo:{name}",
                                unix_time,
                                float(getattr(tree, name)),
                                None,
                            )
                    cfg_size = 0
                    if "daq_config" in branches:
                        cfg = getattr(tree, "daq_config")
                        cfg_size = len(str(cfg))
                    _append_point(
                        series, "runinfo:has_daq_config", unix_time,
                        1.0 if cfg_size > 0 else 0.0, None
                    )
                    _append_point(
                        series, "runinfo:daq_config_bytes", unix_time,
                        float(cfg_size), None
                    )

        return _finalise_payload(path, series, n_entries, t0, t1, "PyROOT")
    finally:
        fin.Close()


def _root_quote(path: Path) -> str:
    text = str(path)
    return text.replace("\\", "\\\\").replace('"', '\\"')


def _load_with_root_cli(path: Path) -> Dict[str, Any]:
    root_exe = shutil.which("root")
    if not root_exe:
        raise ImportError("'root' command is not in PATH.")

    macro = r'''
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void dump_epics_for_viewer(const char *input, const char *output)
{
    TFile *fin = TFile::Open(input, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "cannot open " << input << std::endl;
        return;
    }
    std::ofstream fout(output);
    if (!fout.is_open()) {
        std::cerr << "cannot open output " << output << std::endl;
        fin->Close();
        return;
    }
    Long64_t nentries = 0;
    if (fin->Get("epics"))   nentries += ((TTree*)fin->Get("epics"))->GetEntries();
    if (fin->Get("scalers")) nentries += ((TTree*)fin->Get("scalers"))->GetEntries();
    if (fin->Get("runinfo")) nentries += ((TTree*)fin->Get("runinfo"))->GetEntries();
    fout << "#entries\t" << nentries << "\n";
    fout << std::setprecision(17);

    auto emit = [&](UInt_t unix_time, bool good,
                    const std::string &name, double value) {
        if (unix_time == 0) return;
        fout << unix_time << "\t" << (good ? 1 : 0) << "\t"
             << name << "\t" << value << "\n";
    };

    TTree *tree = (TTree*)fin->Get("epics");
    if (tree && tree->GetBranch("unix_time") &&
        tree->GetBranch("channel") && tree->GetBranch("value")) {
        UInt_t unix_time = 0;
        Bool_t good = true;
        std::vector<std::string> *channel = nullptr;
        std::vector<double> *value = nullptr;
        tree->SetBranchAddress("unix_time", &unix_time);
        tree->SetBranchAddress("channel", &channel);
        tree->SetBranchAddress("value", &value);
        bool has_good = tree->GetBranch("good") != nullptr;
        if (has_good) tree->SetBranchAddress("good", &good);
        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            if (!channel || !value) continue;
            size_t n = std::min(channel->size(), value->size());
            for (size_t j = 0; j < n; ++j)
                emit(unix_time, has_good ? good : true, (*channel)[j], (*value)[j]);
        }
    }

    tree = (TTree*)fin->Get("scalers");
    if (tree && tree->GetBranch("unix_time")) {
        Int_t event_number = 0, slot = 0;
        Long64_t ti_ticks = 0;
        UInt_t unix_time = 0, sync_counter = 0, run_number = 0;
        UChar_t trigger_type = 0, source = 0, channel = 0;
        UInt_t gated = 0, ungated = 0, ref_gated = 0, ref_ungated = 0;
        Float_t live_ratio = 0.f;
        UInt_t trg_gated[16] = {}, trg_ungated[16] = {};
        UInt_t tdc_gated[16] = {}, tdc_ungated[16] = {};
        Bool_t good = true;
        auto bind = [&](const char *name, void *addr) {
            if (tree->GetBranch(name)) tree->SetBranchAddress(name, addr);
        };
        bind("event_number", &event_number);
        bind("ti_ticks", &ti_ticks);
        bind("unix_time", &unix_time);
        bind("sync_counter", &sync_counter);
        bind("run_number", &run_number);
        bind("trigger_type", &trigger_type);
        bind("slot", &slot);
        bind("gated", &gated);
        bind("ungated", &ungated);
        bind("live_ratio", &live_ratio);
        bind("source", &source);
        bind("channel", &channel);
        bind("ref_gated", &ref_gated);
        bind("ref_ungated", &ref_ungated);
        bind("trg_gated", trg_gated);
        bind("trg_ungated", trg_ungated);
        bind("tdc_gated", tdc_gated);
        bind("tdc_ungated", tdc_ungated);
        bind("good", &good);

        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            emit(unix_time, good, "scalers:event_number", event_number);
            emit(unix_time, good, "scalers:ti_ticks", (double)ti_ticks);
            emit(unix_time, good, "scalers:sync_counter", sync_counter);
            emit(unix_time, good, "scalers:run_number", run_number);
            emit(unix_time, good, "scalers:trigger_type", trigger_type);
            emit(unix_time, good, "scalers:slot", slot);
            emit(unix_time, good, "scalers:gated", gated);
            emit(unix_time, good, "scalers:ungated", ungated);
            emit(unix_time, good, "scalers:live_ratio", live_ratio);
            emit(unix_time, good, "scalers:source", source);
            emit(unix_time, good, "scalers:channel", channel);
            emit(unix_time, good, "scalers:ref_gated", ref_gated);
            emit(unix_time, good, "scalers:ref_ungated", ref_ungated);
            for (int ch = 0; ch < 16; ++ch) {
                emit(unix_time, good, "scalers:trg_gated[" + std::to_string(ch) + "]", trg_gated[ch]);
                emit(unix_time, good, "scalers:trg_ungated[" + std::to_string(ch) + "]", trg_ungated[ch]);
                emit(unix_time, good, "scalers:tdc_gated[" + std::to_string(ch) + "]", tdc_gated[ch]);
                emit(unix_time, good, "scalers:tdc_ungated[" + std::to_string(ch) + "]", tdc_ungated[ch]);
            }
        }
    }

    tree = (TTree*)fin->Get("runinfo");
    if (tree && tree->GetBranch("unix_time")) {
        UInt_t run_number = 0, unix_time = 0;
        UChar_t run_type = 0, event_tag = 0;
        std::string daq_config;
        std::string *daq_config_p = &daq_config;
        auto bind = [&](const char *name, void *addr) {
            if (tree->GetBranch(name)) tree->SetBranchAddress(name, addr);
        };
        bind("run_number", &run_number);
        bind("unix_time", &unix_time);
        bind("run_type", &run_type);
        bind("event_tag", &event_tag);
        if (tree->GetBranch("daq_config"))
            tree->SetBranchAddress("daq_config", &daq_config_p);
        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            daq_config.clear();
            tree->GetEntry(i);
            emit(unix_time, true, "runinfo:run_number", run_number);
            emit(unix_time, true, "runinfo:run_type", run_type);
            emit(unix_time, true, "runinfo:event_tag", event_tag);
            emit(unix_time, true, "runinfo:has_daq_config", daq_config.empty() ? 0.0 : 1.0);
            emit(unix_time, true, "runinfo:daq_config_bytes", (double)daq_config.size());
        }
    }
    fout.close();
    fin->Close();
}
'''

    with tempfile.TemporaryDirectory(prefix="epics_viewer_") as tmp:
        tmp_dir = Path(tmp)
        macro_path = tmp_dir / "dump_epics_for_viewer.C"
        out_path = tmp_dir / "epics.tsv"
        macro_path.write_text(macro)

        call = (
            f'{macro_path}("'
            f'{_root_quote(path)}","{_root_quote(out_path)}")'
        )
        proc = subprocess.run(
            [root_exe, "-l", "-b", "-q", call],
            text=True,
            capture_output=True,
            timeout=180,
        )
        if proc.returncode != 0 or not out_path.exists():
            msg = (proc.stderr or proc.stdout or "").strip()
            raise EpicsLoadError(f"ROOT command failed while reading EPICS tree: {msg}")

        return _parse_root_cli_tsv(path, out_path)


def _load_with_pyroot_subprocess(path: Path) -> Dict[str, Any]:
    helper = r'''
import math
import os
import sys
import traceback

SCALER_SCALARS = __SCALER_SCALARS__
SCALER_ARRAYS = __SCALER_ARRAYS__
RUNINFO_SCALARS = __RUNINFO_SCALARS__


def emit(fout, unix_time, good, name, value):
    try:
        t = float(unix_time)
        v = float(value)
    except (TypeError, ValueError):
        return
    if not math.isfinite(t) or t <= 0 or not math.isfinite(v):
        return
    fout.write("%d\t%d\t%s\t%.17g\n" % (int(t), 1 if good else 0, name, v))


def branch_names(tree):
    return set(br.GetName() for br in tree.GetListOfBranches())


def main(input_path, output_path):
    import ROOT  # noqa: F401

    ROOT.gROOT.SetBatch(True)
    fin = ROOT.TFile.Open(input_path, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError("cannot open %s" % input_path)

    with open(output_path, "w") as fout:
        n_entries = 0
        for tree_name in ("epics", "scalers", "runinfo"):
            tree = fin.Get(tree_name)
            if tree:
                n_entries += int(tree.GetEntries())
        fout.write("#entries\t%d\n" % n_entries)

        tree = fin.Get("epics")
        if tree:
            branches = branch_names(tree)
            if {"unix_time", "channel", "value"} <= branches:
                has_good = "good" in branches
                for i in range(int(tree.GetEntries())):
                    tree.GetEntry(i)
                    unix_time = getattr(tree, "unix_time")
                    good = bool(getattr(tree, "good")) if has_good else True
                    channels = getattr(tree, "channel")
                    values = getattr(tree, "value")
                    n = min(int(channels.size()), int(values.size()))
                    for j in range(n):
                        emit(fout, unix_time, good, str(channels[j]), values[j])

        tree = fin.Get("scalers")
        if tree:
            branches = branch_names(tree)
            if "unix_time" in branches:
                has_good = "good" in branches
                for i in range(int(tree.GetEntries())):
                    tree.GetEntry(i)
                    unix_time = getattr(tree, "unix_time")
                    good = bool(getattr(tree, "good")) if has_good else True
                    for name in SCALER_SCALARS:
                        if name in branches:
                            emit(fout, unix_time, good,
                                 "scalers:" + name, getattr(tree, name))
                    for name in SCALER_ARRAYS:
                        if name not in branches:
                            continue
                        arr = getattr(tree, name)
                        for ch in range(16):
                            emit(fout, unix_time, good,
                                 "scalers:%s[%d]" % (name, ch), arr[ch])

        tree = fin.Get("runinfo")
        if tree:
            branches = branch_names(tree)
            if "unix_time" in branches:
                for i in range(int(tree.GetEntries())):
                    tree.GetEntry(i)
                    unix_time = getattr(tree, "unix_time")
                    for name in RUNINFO_SCALARS:
                        if name in branches:
                            emit(fout, unix_time, True,
                                 "runinfo:" + name, getattr(tree, name))
                    cfg_size = 0
                    if "daq_config" in branches:
                        cfg_size = len(str(getattr(tree, "daq_config")))
                    emit(fout, unix_time, True,
                         "runinfo:has_daq_config", 1.0 if cfg_size > 0 else 0.0)
                    emit(fout, unix_time, True,
                         "runinfo:daq_config_bytes", float(cfg_size))

    fin.Close()


if __name__ == "__main__":
    try:
        main(sys.argv[1], sys.argv[2])
    except Exception:
        traceback.print_exc()
        sys.stderr.flush()
        os._exit(2)
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)
'''
    helper = helper.replace("__SCALER_SCALARS__", repr(SCALER_SCALARS))
    helper = helper.replace("__SCALER_ARRAYS__", repr(SCALER_ARRAYS))
    helper = helper.replace("__RUNINFO_SCALARS__", repr(RUNINFO_SCALARS))

    with tempfile.TemporaryDirectory(prefix="epics_viewer_pyroot_") as tmp:
        tmp_dir = Path(tmp)
        script_path = tmp_dir / "dump_epics_for_viewer.py"
        out_path = tmp_dir / "epics.tsv"
        script_path.write_text(helper)

        proc = subprocess.run(
            [sys.executable, str(script_path), str(path), str(out_path)],
            text=True,
            capture_output=True,
            timeout=180,
        )
        if proc.returncode != 0 or not out_path.exists():
            msg = (proc.stderr or proc.stdout or "").strip()
            if "No module named" in msg or "ImportError" in msg:
                raise ImportError("PyROOT is not available in this Python.")
            raise EpicsLoadError(
                f"PyROOT subprocess failed while reading EPICS tree: {msg}"
            )

        return _parse_root_cli_tsv(path, out_path, "PyROOT subprocess")


def _parse_root_cli_tsv(
    path: Path,
    out_path: Path,
    loader: str = "root command",
) -> Dict[str, Any]:
    series: Dict[str, Dict[str, List[Any]]] = {}
    n_entries = 0
    t0: Optional[float] = None
    t1: Optional[float] = None

    with out_path.open() as fin:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#entries\t"):
                try:
                    n_entries = int(line.split("\t", 1)[1])
                except ValueError:
                    n_entries = 0
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            unix_time = float(parts[0])
            good = parts[1] != "0"
            channel = parts[2]
            value = float(parts[3])
            if math.isfinite(unix_time) and unix_time > 0:
                t0 = unix_time if t0 is None else min(t0, unix_time)
                t1 = unix_time if t1 is None else max(t1, unix_time)
            _append_point(series, channel, unix_time, value, good)

    return _finalise_payload(path, series, n_entries, t0, t1, loader)


def _load_with_uproot(path: Path) -> Dict[str, Any]:
    try:
        import awkward as ak  # type: ignore
        import uproot  # type: ignore
    except Exception as exc:
        raise ImportError("uproot/awkward is not available.") from exc

    series: Dict[str, Dict[str, List[Any]]] = {}
    t0: Optional[float] = None
    t1: Optional[float] = None
    n_entries = 0

    def note_time(unix_time: float) -> None:
        nonlocal t0, t1
        if not math.isfinite(unix_time) or unix_time <= 0:
            return
        t0 = unix_time if t0 is None else min(t0, unix_time)
        t1 = unix_time if t1 is None else max(t1, unix_time)

    with uproot.open(str(path)) as fin:
        if "epics" in fin:
            tree = fin["epics"]
            branches = set(tree.keys())
            if {"unix_time", "channel", "value"} <= branches:
                wanted = ["unix_time", "channel", "value"]
                if "good" in branches:
                    wanted.append("good")
                arrays = tree.arrays(wanted, library="ak")
                times = ak.to_list(arrays["unix_time"])
                channel_rows = ak.to_list(arrays["channel"])
                value_rows = ak.to_list(arrays["value"])
                good_rows = (
                    ak.to_list(arrays["good"])
                    if "good" in arrays.fields else [None] * len(times)
                )
                n_entries += len(times)
                for unix_time_raw, channels, values, good in zip(
                    times, channel_rows, value_rows, good_rows
                ):
                    unix_time = float(unix_time_raw)
                    note_time(unix_time)
                    for channel, value in zip(channels, values):
                        _append_point(series, str(channel), unix_time, float(value), good)

        if "scalers" in fin:
            tree = fin["scalers"]
            branches = set(tree.keys())
            if "unix_time" in branches:
                wanted = ["unix_time"]
                wanted += [name for name in SCALER_SCALARS if name in branches]
                wanted += [name for name in SCALER_ARRAYS if name in branches]
                if "good" in branches:
                    wanted.append("good")
                arrays = tree.arrays(wanted, library="ak")
                times = ak.to_list(arrays["unix_time"])
                good_rows = (
                    ak.to_list(arrays["good"])
                    if "good" in arrays.fields else [None] * len(times)
                )
                n_entries += len(times)
                scalar_data = {
                    name: ak.to_list(arrays[name])
                    for name in SCALER_SCALARS if name in arrays.fields
                }
                array_data = {
                    name: ak.to_list(arrays[name])
                    for name in SCALER_ARRAYS if name in arrays.fields
                }
                for i, unix_time_raw in enumerate(times):
                    unix_time = float(unix_time_raw)
                    note_time(unix_time)
                    good = good_rows[i]
                    for name, values in scalar_data.items():
                        _append_point(series, f"scalers:{name}", unix_time,
                                      float(values[i]), good)
                    for name, rows in array_data.items():
                        row = rows[i]
                        for ch, value in enumerate(row[:16]):
                            _append_point(series, f"scalers:{name}[{ch}]",
                                          unix_time, float(value), good)

        if "runinfo" in fin:
            tree = fin["runinfo"]
            branches = set(tree.keys())
            if "unix_time" in branches:
                wanted = ["unix_time"]
                wanted += [name for name in RUNINFO_SCALARS if name in branches]
                if "daq_config" in branches:
                    wanted.append("daq_config")
                arrays = tree.arrays(wanted, library="ak")
                times = ak.to_list(arrays["unix_time"])
                n_entries += len(times)
                scalar_data = {
                    name: ak.to_list(arrays[name])
                    for name in RUNINFO_SCALARS if name in arrays.fields
                }
                cfg_rows = (
                    ak.to_list(arrays["daq_config"])
                    if "daq_config" in arrays.fields else [""] * len(times)
                )
                for i, unix_time_raw in enumerate(times):
                    unix_time = float(unix_time_raw)
                    note_time(unix_time)
                    for name, values in scalar_data.items():
                        _append_point(series, f"runinfo:{name}", unix_time,
                                      float(values[i]), None)
                    cfg = "" if cfg_rows[i] is None else str(cfg_rows[i])
                    _append_point(series, "runinfo:has_daq_config", unix_time,
                                  1.0 if cfg else 0.0, None)
                    _append_point(series, "runinfo:daq_config_bytes", unix_time,
                                  float(len(cfg)), None)

    return _finalise_payload(path, series, n_entries, t0, t1, "uproot")


# ===========================================================================
#  Plot helpers
# ===========================================================================

def _nice_ticks(lo: float, hi: float, max_ticks: int = 6) -> List[float]:
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return [lo if math.isfinite(lo) else 0.0]
    raw = (hi - lo) / max(max_ticks - 1, 1)
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    step = mag
    for c in (1, 2, 2.5, 5, 10):
        if c * mag >= raw:
            step = c * mag
            break
    start = math.ceil(lo / step) * step
    out: List[float] = []
    v = start
    guard = 0
    while v <= hi + step * 0.01 and guard < 100:
        out.append(v)
        v += step
        guard += 1
    return out


def _fmt_num(v: float) -> str:
    if not math.isfinite(v):
        return ""
    a = abs(v)
    if a >= 1e4 or (0 < a < 1e-3):
        return f"{v:.2e}"
    if a >= 100:
        return f"{v:.1f}"
    if a >= 10:
        return f"{v:.2f}"
    return f"{v:.4g}"


def _fmt_date_range(data: Dict[str, Any]) -> str:
    try:
        t0 = int(data.get("t0_unix"))
        t1 = int(data.get("t1_unix"))
    except (TypeError, ValueError):
        return ""
    if t0 <= 0 or t1 <= 0:
        return ""
    d0 = datetime.fromtimestamp(t0)
    d1 = datetime.fromtimestamp(t1)
    if d0.date() == d1.date():
        return f"{d0:%Y-%m-%d %H:%M} - {d1:%H:%M}"
    return f"{d0:%Y-%m-%d %H:%M} - {d1:%Y-%m-%d %H:%M}"


def _default_channels(channels: List[str]) -> List[str]:
    preferred = [
        "DAQ live time [%]",
        "scalers:live_ratio",
        "hallb_IPM2C21A_CUR",
        "hallb_IPM2C21A_XPOS",
        "hallb_IPM2C21A_YPOS",
        "TGT:PRad:Cell_P",
        "runinfo:event_tag",
    ]
    out: List[str] = []
    for ch in preferred:
        if ch in channels:
            out.append(ch)
    for ch in channels:
        if len(out) >= 4:
            break
        if ch not in out:
            out.append(ch)
    while len(out) < 4:
        out.append(channels[0] if channels else "")
    return out


class ScatterPlotWidget(QWidget):
    PAD_L, PAD_R, PAD_T, PAD_B = 104, 18, 20, 58

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(230, 170)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)
        self._data: Dict[str, Any] = {}
        self._channel = ""
        self._x_range: Optional[Tuple[float, float]] = None
        self._y_range: Optional[Tuple[float, float]] = None
        self._drag_start: Optional[Tuple[float, float]] = None
        self._drag_now: Optional[Tuple[float, float]] = None

    def set_data(self, data: Dict[str, Any], channel: str) -> None:
        self._data = data
        self._channel = channel
        self.reset_zoom()

    def set_channel(self, channel: str) -> None:
        self._channel = channel
        self.reset_zoom()

    def reset_zoom(self) -> None:
        self._x_range = None
        self._y_range = None
        self._drag_start = None
        self._drag_now = None
        self.update()

    def _plot_rect(self) -> QRectF:
        return QRectF(
            self.PAD_L,
            self.PAD_T,
            max(10.0, self.width() - self.PAD_L - self.PAD_R),
            max(10.0, self.height() - self.PAD_T - self.PAD_B),
        )

    def _series(self) -> Dict[str, Any]:
        return (self._data.get("series") or {}).get(self._channel) or {}

    def _full_range(self) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        xs = self._series().get("x") or []
        ys = self._series().get("y") or []
        xmax = float(self._data.get("duration_min") or 0.0)
        ymin = float("inf")
        ymax = float("-inf")
        for x, y in zip(xs, ys):
            try:
                xf = float(x)
                yf = float(y)
            except (TypeError, ValueError):
                continue
            if math.isfinite(xf):
                xmax = max(xmax, xf)
            if math.isfinite(yf):
                ymin = min(ymin, yf)
                ymax = max(ymax, yf)
        if xmax <= 0:
            xmax = 1.0
        if not math.isfinite(ymin) or not math.isfinite(ymax):
            ymin, ymax = 0.0, 1.0
        elif ymax <= ymin:
            pad = max(1.0, abs(ymin) * 0.05)
            ymin -= pad
            ymax += pad
        else:
            pad = (ymax - ymin) * 0.08
            ymin -= pad
            ymax += pad
        return (0.0, xmax), (ymin, ymax)

    def _active_range(self) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        xr, yr = self._full_range()
        return self._x_range or xr, self._y_range or yr

    @staticmethod
    def _to_sx(x: float, rect: QRectF, xr: Tuple[float, float]) -> float:
        return rect.left() + (x - xr[0]) / (xr[1] - xr[0]) * rect.width()

    @staticmethod
    def _to_sy(y: float, rect: QRectF, yr: Tuple[float, float]) -> float:
        return rect.bottom() - (y - yr[0]) / (yr[1] - yr[0]) * rect.height()

    @staticmethod
    def _from_sx(sx: float, rect: QRectF, xr: Tuple[float, float]) -> float:
        return xr[0] + (sx - rect.left()) / rect.width() * (xr[1] - xr[0])

    @staticmethod
    def _from_sy(sy: float, rect: QRectF, yr: Tuple[float, float]) -> float:
        return yr[0] + (rect.bottom() - sy) / rect.height() * (yr[1] - yr[0])

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.RightButton:
            self.reset_zoom()
            event.accept()
            return
        if event.button() != Qt.MouseButton.LeftButton:
            return
        pos = event.position()
        if self._plot_rect().contains(pos):
            self._drag_start = (pos.x(), pos.y())
            self._drag_now = self._drag_start
            self.update()

    def mouseMoveEvent(self, event) -> None:
        if self._drag_start is not None:
            pos = event.position()
            self._drag_now = (pos.x(), pos.y())
            self.update()

    def mouseReleaseEvent(self, event) -> None:
        if event.button() != Qt.MouseButton.LeftButton or self._drag_start is None:
            return
        pos = event.position()
        end = (pos.x(), pos.y())
        start = self._drag_start
        self._drag_start = None
        self._drag_now = None

        if abs(end[0] - start[0]) > 8 and abs(end[1] - start[1]) > 8:
            rect = self._plot_rect()
            xr, yr = self._active_range()
            x0 = max(rect.left(), min(rect.right(), start[0]))
            x1 = max(rect.left(), min(rect.right(), end[0]))
            y0 = max(rect.top(), min(rect.bottom(), start[1]))
            y1 = max(rect.top(), min(rect.bottom(), end[1]))
            xa = self._from_sx(min(x0, x1), rect, xr)
            xb = self._from_sx(max(x0, x1), rect, xr)
            ya = self._from_sy(max(y0, y1), rect, yr)
            yb = self._from_sy(min(y0, y1), rect, yr)
            if xb > xa and yb > ya:
                self._x_range = (xa, xb)
                self._y_range = (ya, yb)
        self.update()

    def mouseDoubleClickEvent(self, event) -> None:
        self.reset_zoom()
        event.accept()

    def contextMenuEvent(self, event) -> None:
        self.reset_zoom()
        event.accept()

    def paintEvent(self, _event) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(PANEL))

        rect = self._plot_rect()
        xr, yr = self._active_range()
        p.setPen(QPen(QColor(BORDER), 1))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(rect)

        p.setFont(QFont("Monospace", 12))
        for tick in _nice_ticks(xr[0], xr[1], 7):
            sx = self._to_sx(tick, rect, xr)
            p.setPen(QPen(QColor(GRID), 1))
            p.drawLine(int(sx), int(rect.top()), int(sx), int(rect.bottom()))
            p.setPen(QColor(TEXT_DIM))
            p.drawText(
                QRectF(sx - 48, rect.bottom() + 7, 96, 22),
                Qt.AlignmentFlag.AlignHCenter,
                _fmt_num(tick),
            )

        for tick in _nice_ticks(yr[0], yr[1], 10):
            sy = self._to_sy(tick, rect, yr)
            p.setPen(QPen(QColor(GRID), 1))
            p.drawLine(int(rect.left()), int(sy), int(rect.right()), int(sy))
            p.setPen(QColor(TEXT_DIM))
            p.drawText(
                QRectF(24, sy - 11, self.PAD_L - 32, 22),
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                _fmt_num(tick),
            )

        series = self._series()
        xs = series.get("x") or []
        ys = series.get("y") or []
        good = series.get("good") or []

        p.setClipRect(rect)
        for i, (x, y) in enumerate(zip(xs, ys)):
            try:
                xf = float(x)
                yf = float(y)
            except (TypeError, ValueError):
                continue
            if xf < xr[0] or xf > xr[1] or yf < yr[0] or yf > yr[1]:
                continue
            sx = self._to_sx(xf, rect, xr)
            sy = self._to_sy(yf, rect, yr)
            is_good = i >= len(good) or bool(good[i])
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(POINT if is_good else POINT_BAD))
            p.drawEllipse(QRectF(sx - 2.0, sy - 2.0, 4.0, 4.0))
        p.setClipping(False)

        if self._drag_start and self._drag_now:
            x0, y0 = self._drag_start
            x1, y1 = self._drag_now
            x0 = max(rect.left(), min(rect.right(), x0))
            x1 = max(rect.left(), min(rect.right(), x1))
            y0 = max(rect.top(), min(rect.bottom(), y0))
            y1 = max(rect.top(), min(rect.bottom(), y1))
            zr = QRectF(min(x0, x1), min(y0, y1), abs(x1 - x0), abs(y1 - y0))
            p.fillRect(zr, QColor(83, 182, 255, 45))
            p.setPen(QPen(QColor(ACCENT), 1))
            p.drawRect(zr)

        p.setFont(QFont("Monospace", 12))
        p.setPen(QColor(TEXT_DIM))
        tag = "zoomed" if self._x_range else "full"
        p.drawText(
            QRectF(rect.left(), 0, rect.width(), 18),
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
            f"{len(xs)} points | {tag}",
        )
        p.setFont(QFont("Monospace", 10, QFont.Weight.Bold))
        p.drawText(
            QRectF(rect.left(), rect.bottom() + 31, rect.width(), 18),
            Qt.AlignmentFlag.AlignHCenter,
            "time from first EPICS point [min]",
        )


class PlotPanel(QFrame):
    def __init__(self, index: int, parent=None):
        super().__init__(parent)
        self._data: Dict[str, Any] = {}
        self.setStyleSheet(
            f"QFrame{{background:{PANEL};border:1px solid {BORDER};border-radius:6px;}}"
            f"QComboBox{{background:#0d1117;color:{TEXT};border:1px solid #33404d;"
            "border-radius:4px;padding:2px 6px;}"
            f"QLabel{{color:{TEXT_DIM};border:0;background:transparent;}}"
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        top = QHBoxLayout()
        idx = QLabel(str(index + 1))
        idx.setAlignment(Qt.AlignmentFlag.AlignCenter)
        idx.setFixedWidth(22)
        idx.setStyleSheet(f"color:{ACCENT};font-weight:bold;")
        top.addWidget(idx)

        self.combo = QComboBox()
        self.combo.setEditable(True)
        self.combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.combo.setToolTip("Type any keyword to search channels, case-insensitive.")
        if self.combo.lineEdit() is not None:
            self.combo.lineEdit().setPlaceholderText("search channel...")
        self.combo.currentTextChanged.connect(self._on_channel_changed)
        top.addWidget(self.combo, 1)
        layout.addLayout(top)

        self.plot = ScatterPlotWidget()
        layout.addWidget(self.plot, 1)

        hint = QLabel("left-drag zoom | right-click/double-click reset")
        hint.setFont(QFont("Monospace", 8))
        layout.addWidget(hint)

    def set_data(self, data: Dict[str, Any], channels: List[str], channel: str) -> None:
        self._data = data
        self.combo.blockSignals(True)
        self.combo.clear()
        self.combo.addItems(channels)
        completer = QCompleter(channels, self.combo)
        completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        completer.setFilterMode(Qt.MatchFlag.MatchContains)
        completer.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        self.combo.setCompleter(completer)
        if channel in channels:
            self.combo.setCurrentText(channel)
        self.combo.blockSignals(False)
        self.plot.set_data(data, self.combo.currentText())

    def _on_channel_changed(self, channel: str) -> None:
        if channel in (self._data.get("series") or {}):
            self.plot.set_channel(channel)


class EpicsPlotGrid(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        header = QHBoxLayout()
        layout.addLayout(header)

        self.summary = QLabel("No file loaded")
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.summary.setStyleSheet(f"color:{TEXT_DIM};")
        header.addWidget(self.summary, 1)

        self.date_range = QLabel("")
        self.date_range.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.date_range.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.date_range.setFont(QFont("Monospace", 13, QFont.Weight.Bold))
        self.date_range.setStyleSheet(f"color:{ACCENT};")
        header.addWidget(self.date_range, 0)

        grid = QGridLayout()
        grid.setSpacing(10)
        layout.addLayout(grid, 1)

        self.panels = [PlotPanel(i) for i in range(4)]
        for i, panel in enumerate(self.panels):
            grid.addWidget(panel, 0, i)
            grid.setColumnStretch(i, 1)

    def set_data(self, data: Dict[str, Any]) -> None:
        channels = list(data.get("channels") or [])
        defaults = _default_channels(channels)
        duration = float(data.get("duration_min") or 0.0)
        self.summary.setText(
            f"{data.get('title', 'EPICS ROOT file')} | "
            f"{len(channels)} channels | {data.get('n_entries', 0)} entries | "
            f"{data.get('n_points', 0)} points | {duration:.2f} min | "
            f"loader: {data.get('loader', 'unknown')}"
        )
        self.date_range.setText(_fmt_date_range(data))
        for panel, channel in zip(self.panels, defaults):
            panel.set_data(data, channels, channel)

    def set_message(self, message: str) -> None:
        self.summary.setText(message)
        self.date_range.setText("")


# ===========================================================================
#  Main window
# ===========================================================================

class EpicsViewer(QMainWindow):
    def __init__(
        self,
        initial_path: Optional[str] = None,
        allow_uproot: bool = False,
        allow_pyroot: bool = False,
    ):
        super().__init__()
        self.setWindowTitle("PRad-II EPICS viewer")
        self.resize(1550, 760)

        self._current_path: Optional[str] = None
        self._loader: Optional[EpicsLoader] = None
        self._allow_uproot = allow_uproot
        self._allow_pyroot = allow_pyroot

        self._build_ui()
        self._build_menu()
        if initial_path:
            self.open_root_file(initial_path)

    def _build_ui(self) -> None:
        central = QWidget()
        central.setStyleSheet(f"background:{BG};color:{TEXT};")
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        top = QHBoxLayout()
        layout.addLayout(top)

        self.open_btn = QPushButton("Open ROOT file...")
        self.open_btn.clicked.connect(self.on_open)
        top.addWidget(self.open_btn)

        self.reload_btn = QPushButton("Reload")
        self.reload_btn.setEnabled(False)
        self.reload_btn.clicked.connect(self.on_reload)
        top.addWidget(self.reload_btn)

        self.uproot_cb = QCheckBox("allow uproot fallback")
        self.uproot_cb.setChecked(self._allow_uproot)
        self.uproot_cb.toggled.connect(self._set_uproot)
        top.addWidget(self.uproot_cb)

        self.path_label = QLabel("No file loaded")
        self.path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.path_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        top.addWidget(self.path_label, 1)

        self.grid = EpicsPlotGrid()
        layout.addWidget(self.grid, 1)

        self.status = QLabel("Ready")
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.status.setStyleSheet(f"color:{TEXT_DIM};")
        layout.addWidget(self.status)

    def _build_menu(self) -> None:
        menu = self.menuBar().addMenu("&File")
        open_act = QAction("Open ROOT file...", self)
        open_act.triggered.connect(self.on_open)
        menu.addAction(open_act)
        reload_act = QAction("Reload", self)
        reload_act.triggered.connect(self.on_reload)
        menu.addAction(reload_act)
        menu.addSeparator()
        quit_act = QAction("Quit", self)
        quit_act.triggered.connect(self.close)
        menu.addAction(quit_act)

    def _set_uproot(self, checked: bool) -> None:
        self._allow_uproot = checked

    def on_open(self) -> None:
        start = str(Path(self._current_path).parent) if self._current_path else ""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open ROOT file with epics tree",
            start,
            "ROOT files (*.root);;All files (*)",
        )
        if path:
            self.open_root_file(path)

    def on_reload(self) -> None:
        if self._current_path:
            self.open_root_file(self._current_path)

    def open_root_file(self, path: str) -> None:
        if self._loader is not None and self._loader.isRunning():
            QMessageBox.information(self, "Still loading", "Please wait for loading to finish.")
            return
        self._current_path = path
        self.path_label.setText(path)
        self.status.setText(f"Loading {path} ...")
        self.open_btn.setEnabled(False)
        self.reload_btn.setEnabled(False)

        self._loader = EpicsLoader(
            path,
            allow_uproot=self._allow_uproot,
            allow_pyroot=self._allow_pyroot,
            parent=self,
        )
        self._loader.finished.connect(self._on_loaded)
        self._loader.start()

    def _on_loaded(self, data: Dict[str, Any], error: str) -> None:
        self.open_btn.setEnabled(True)
        self.reload_btn.setEnabled(self._current_path is not None)

        if error:
            self.status.setText(f"Error: {error}")
            self.grid.set_message(error)
            QMessageBox.critical(self, "Failed to load EPICS tree", error)
            return

        title = data.get("title", "EPICS")
        channels = data.get("channels") or []
        n_points = data.get("n_points", 0)
        duration = float(data.get("duration_min") or 0.0)
        self.status.setText(
            f"Loaded {title}: {len(channels)} channels, {n_points} points, "
            f"{duration:.2f} min"
        )
        self.grid.set_data(data)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Pure PyQt6 viewer for PRad-II epics/scalers/runinfo ROOT trees."
    )
    parser.add_argument("root_file", nargs="?", help="ROOT file containing replay side trees.")
    parser.add_argument(
        "--uproot",
        action="store_true",
        help="Allow uproot fallback if ROOT-based readers are unavailable.",
    )
    parser.add_argument(
        "--pyroot",
        action="store_true",
        help=(
            "Allow PyROOT fallback inside the GUI process. This is disabled by "
            "default because some ROOT/Qt combinations crash during shutdown."
        ),
    )
    args = parser.parse_args(argv)

    app = QApplication.instance() or QApplication(sys.argv)
    win = EpicsViewer(
        args.root_file,
        allow_uproot=args.uproot,
        allow_pyroot=args.pyroot,
    )
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
