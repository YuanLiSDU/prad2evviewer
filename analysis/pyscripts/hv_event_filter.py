#!/usr/bin/env python3
"""
hv_event_filter.py — flag physics events by HV stability.

Reads a replay ROOT file that already has the `hv` / `hv_channels` trees
appended (see `prad2ana_replay_add_hv`).  For every HV snapshot we compute
a per-channel rolling std of dV; a snapshot is "stable" iff every monitored
channel passes BOTH gates the user enabled:

  * absolute  : rolling-std(dV) ≤ --abs-threshold   (volts)
  * relative  : rolling-std(dV) ≤ median(rolling-std)
                                + --rel-n-mad · 1.4826 · MAD(rolling-std)
                MAD is the median absolute deviation of the per-snapshot
                rolling-std distribution; 1.4826 makes it a Gaussian-σ
                equivalent.  The MAD method is robust to a few transient
                spikes — it lets us learn "what's normal" for each channel
                without a few HV trips poisoning the threshold.

Snapshots that fail are propagated to physics events through the
`event_number_at_arrival` column the C++ tool writes: each event picks up
the stability flag of the most recent HV snapshot whose
`event_number_at_arrival ≤ event.event_num`.

Outputs:
  * stdout summary — total events, stable fraction, top-N rejecting channels.
  * --out-events FILE      : one stable event_num per line.
  * --out-json    FILE     : full report as JSON (parameters, counts,
                             per-channel rejection histogram, monitored
                             channel list).

Channel scope (one of):
  --scope all                          : every channel that lives in the
                                         hycal map (default).
  --scope type   --type-name PbWO4     : module type filter.
  --scope channel --channel-name W100  : single channel.

Usage:
  python3 hv_event_filter.py prad_024340_filtered.root \\
      --abs-threshold 0.5 --rel-n-mad 5 --window-s 5 \\
      --out-json hv_filter.json --out-events stable_events.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

try:
    import uproot
except ImportError as e:
    raise SystemExit(
        f"[hv_event_filter] uproot required (pip install uproot): {e}")


DEFAULT_HYCAL_MAP = (Path(__file__).resolve().parents[2]
                     / "database" / "hycal_map.json")


# ───────────────────────────────────────────────────────────────────────
# helpers
# ───────────────────────────────────────────────────────────────────────

def _load_module_types(path: Path) -> dict:
    """name → type ('PbWO4' | 'PbGlass' | 'Veto' | 'LMS')."""
    with open(path) as f:
        records = json.load(f)
    out = {}
    for r in records:
        name = r.get("n")
        mtype = r.get("t")
        if name and mtype:
            out[name] = mtype
    return out


def _rolling_std(arr: np.ndarray, win: int) -> np.ndarray:
    """Centered rolling population-std with window `win`.  NaN-aware:
    samples where arr is NaN are dropped from the running mean/var, so
    an offline channel does not poison its neighbours' rolling stats."""
    n = arr.size
    if n == 0 or win <= 1:
        return np.zeros(n, dtype=np.float32)
    valid = ~np.isnan(arr)
    a = np.where(valid, arr.astype(np.float64, copy=False), 0.0)
    cum  = np.concatenate(([0.0], np.cumsum(a, dtype=np.float64)))
    cum2 = np.concatenate(([0.0], np.cumsum(a * a, dtype=np.float64)))
    cumv = np.concatenate(([0],   np.cumsum(valid.astype(np.int64))))
    half = win // 2
    idx = np.arange(n)
    lo = np.maximum(idx - half, 0)
    hi = np.minimum(idx - half + win, n)
    cnt  = np.maximum(cumv[hi] - cumv[lo], 1)
    s    = cum [hi] - cum [lo]
    s2   = cum2[hi] - cum2[lo]
    mean = s / cnt
    var  = np.maximum(s2 / cnt - mean * mean, 0.0)
    return np.sqrt(var).astype(np.float32)


def _resolve_channels(hv_names: list,
                      module_types: dict,
                      scope: str,
                      type_name: str | None,
                      channel_name: str | None) -> list:
    """Return a list of (idx_in_hv_tree, name, mtype-or-None) tuples."""
    if scope == "all":
        return [(i, n, module_types.get(n))
                for i, n in enumerate(hv_names)
                if n in module_types]
    if scope == "type":
        if not type_name:
            raise SystemExit("[hv_event_filter] --scope type requires --type-name")
        wanted = type_name
        return [(i, n, module_types[n])
                for i, n in enumerate(hv_names)
                if n in module_types and module_types[n] == wanted]
    if scope == "channel":
        if not channel_name:
            raise SystemExit("[hv_event_filter] --scope channel requires --channel-name")
        return [(i, n, module_types.get(n))
                for i, n in enumerate(hv_names)
                if n == channel_name]
    raise SystemExit(f"[hv_event_filter] unknown --scope {scope!r}")


# ───────────────────────────────────────────────────────────────────────
# main
# ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input",
                    help="Replay ROOT file with `hv` tree (run "
                         "prad2ana_replay_add_hv on it first).")
    ap.add_argument("--hycal-map", default=str(DEFAULT_HYCAL_MAP),
                    help=f"hycal_map.json (default: {DEFAULT_HYCAL_MAP})")
    ap.add_argument("--scope", choices=["all", "type", "channel"],
                    default="all",
                    help='Channel scope (default: all = every HV channel '
                         'present in hycal_map.json).')
    ap.add_argument("--type-name", default=None,
                    help="Module type to filter on (PbWO4, PbGlass, Veto, "
                         "LMS) — used with --scope type.")
    ap.add_argument("--channel-name", default=None,
                    help="Single channel to monitor (e.g. W100) — used "
                         "with --scope channel.")
    ap.add_argument("--abs-threshold", type=float, default=None,
                    help="Absolute fluctuation gate: rolling-std(dV) above "
                         "this many volts marks the channel unstable. "
                         "Disabled if not given.")
    ap.add_argument("--rel-n-mad", type=float, default=None,
                    help="Relative fluctuation gate (MAD): rolling-std(dV) "
                         "above median + N · 1.4826 · MAD marks the "
                         "channel unstable.  Disabled if not given.")
    ap.add_argument("--window-s", type=float, default=5.0,
                    help="Rolling-std window in seconds (default 5).")
    ap.add_argument("--require-v0set", action="store_true",
                    help="Also require channel to be at its setpoint "
                         "(|dv| ≤ --abs-threshold).  Off by default — "
                         "channels held intentionally low don't fail.")
    ap.add_argument("--top-n", type=int, default=20,
                    help="How many worst-rejecting channels to print "
                         "(default 20).")
    ap.add_argument("--out-events", default=None,
                    help="Write stable event_num list (one per line) here.")
    ap.add_argument("--out-json", default=None,
                    help="Write full report as JSON here.")
    ap.add_argument("--event-tree", default=None,
                    help="Tree name carrying the physics events "
                         "(default: try `recon`, then `events`).")
    args = ap.parse_args()

    if args.abs_threshold is None and args.rel_n_mad is None:
        sys.exit("[hv_event_filter] need at least one of --abs-threshold "
                 "or --rel-n-mad — neither is set")

    # ── 1. read HV trees ───────────────────────────────────────────────
    with uproot.open(args.input) as f:
        if "hv" not in f:
            sys.exit(f"[hv_event_filter] {args.input}: no `hv` tree "
                     "(run prad2ana_replay_add_hv first)")
        hv_tree = f["hv"]
        ch_tree = f["hv_channels"]

        ev_at_arrival = hv_tree["event_number_at_arrival"].array(library="np").astype(np.int64)
        t_unix_s      = hv_tree["t_unix_s"].array(library="np")
        dv_2d         = hv_tree["dv"].array(library="np")
        names         = ch_tree["name"].array(library="np").tolist()
        names         = [str(n) for n in names]

        # Pick the events-side tree
        evt_tree_name = args.event_tree
        if evt_tree_name is None:
            for cand in ("recon", "events"):
                if cand in f:
                    evt_tree_name = cand
                    break
        if evt_tree_name is None or evt_tree_name not in f:
            sys.exit(f"[hv_event_filter] {args.input}: no `recon` or `events` "
                     f"tree (override with --event-tree).")
        evt_tree    = f[evt_tree_name]
        event_nums  = evt_tree["event_num"].array(library="np").astype(np.int64)

        # Sample period from hv_meta (falls back to median diff)
        if "hv_meta" in f:
            interval_ms = int(f["hv_meta"]["interval_ms"].array(library="np")[0])
        else:
            ts = hv_tree["unix_time"].array(library="np")
            interval_ms = (200 if ts.size < 2
                           else int(np.median(np.diff(ts)) * 1000) or 200)

    win = max(1, int(round(args.window_s * 1000.0 / max(interval_ms, 1))))
    print(f"[hv_event_filter] {args.input}")
    print(f"  hv tree         : {dv_2d.shape[0]:,} snapshots × "
          f"{dv_2d.shape[1]} channels (period={interval_ms} ms)")
    print(f"  event tree      : {evt_tree_name}, "
          f"{event_nums.size:,} entries")

    # ── 2. resolve channel scope ────────────────────────────────────────
    module_types = _load_module_types(Path(args.hycal_map))
    selected = _resolve_channels(names, module_types,
                                  args.scope, args.type_name,
                                  args.channel_name)
    if not selected:
        sys.exit(f"[hv_event_filter] no channels match scope={args.scope}; "
                 "nothing to monitor.")
    print(f"  hycal map       : {args.hycal_map} ({len(module_types)} entries)")
    print(f"  scope           : {args.scope}"
          + (f' (type={args.type_name})' if args.scope == 'type' else '')
          + (f' (channel={args.channel_name})' if args.scope == 'channel' else ''))
    print(f"  monitored ch    : {len(selected)} channel(s)")
    print(f"  rolling window  : {args.window_s} s ({win} samples)")
    if args.abs_threshold is not None:
        print(f"  absolute gate   : rolling-std ≤ {args.abs_threshold:g} V")
        if args.require_v0set:
            print(f"                   AND |dv| ≤ {args.abs_threshold:g} V")
    if args.rel_n_mad is not None:
        print(f"  relative gate   : rolling-std ≤ median + "
              f"{args.rel_n_mad:g} · 1.4826 · MAD")

    # ── 3. per-snapshot stability + per-channel rejection counts ───────
    n_snap = dv_2d.shape[0]
    unstable_overall = np.zeros(n_snap, dtype=bool)
    rejection = OrderedDict()
    per_channel_stats = OrderedDict()

    for ch_idx, ch_name, mtype in selected:
        col = dv_2d[:, ch_idx].astype(np.float64)
        rs  = _rolling_std(col, win)
        bad = np.zeros(n_snap, dtype=bool)
        valid = ~np.isnan(rs)

        if args.abs_threshold is not None:
            bad |= (rs > args.abs_threshold)
            if args.require_v0set:
                bad |= (np.abs(col) > args.abs_threshold)

        med = mad = thr = None
        if args.rel_n_mad is not None and valid.any():
            med = float(np.median(rs[valid]))
            mad = float(np.median(np.abs(rs[valid] - med)))
            thr = med + args.rel_n_mad * 1.4826 * mad
            bad |= (rs > thr)

        bad |= np.isnan(col)        # offline samples are unstable
        n_bad = int(bad.sum())
        rejection[ch_name] = n_bad
        per_channel_stats[ch_name] = {
            "type": mtype,
            "n_unstable_snapshots": n_bad,
            "rs_median": med,
            "rs_mad":    mad,
            "rs_thresh_rel": thr,
        }
        unstable_overall |= bad

    stable_per_snap = ~unstable_overall

    n_stable_snap   = int(stable_per_snap.sum())
    n_unstable_snap = int(unstable_overall.sum())
    print()
    print(f"  HV snapshots    : {n_snap:,} total, "
          f"{n_stable_snap:,} stable ({100.0*n_stable_snap/n_snap:.2f}%), "
          f"{n_unstable_snap:,} unstable")

    # ── 4. propagate to physics events ─────────────────────────────────
    # `ev_at_arrival[i]` is the event_num closest to HV snapshot i; sort
    # by it so we can searchsorted for each physics event.
    order = np.argsort(ev_at_arrival, kind="mergesort")
    sorted_ev_at_arrival = ev_at_arrival[order]
    sorted_stable        = stable_per_snap[order]

    # For each event: find largest HV-snapshot index whose
    # event_number_at_arrival ≤ event_num.  That snapshot's stability
    # is the event's HV-quality flag.
    j = np.searchsorted(sorted_ev_at_arrival, event_nums, side="right") - 1
    out_of_range = j < 0
    j_safe = np.maximum(j, 0)
    event_stable = np.where(out_of_range, False, sorted_stable[j_safe])

    n_events     = int(event_nums.size)
    n_evt_stable = int(event_stable.sum())
    n_evt_pre_hv = int(out_of_range.sum())
    print(f"  events          : {n_events:,} total, "
          f"{n_evt_stable:,} stable ({100.0*n_evt_stable/n_events:.2f}%)")
    if n_evt_pre_hv:
        print(f"                    "
              f"({n_evt_pre_hv:,} before first HV snapshot — counted unstable)")

    # ── 5. report on most-rejecting channels ───────────────────────────
    print()
    print(f"  Top {args.top_n} rejecting channel(s):")
    print(f"    {'Channel':<12s}  {'type':<8s}  "
          f"{'n_unstable':>12s}  {'frac':>7s}")
    print("    " + "-" * 46)
    sorted_rej = sorted(rejection.items(), key=lambda kv: -kv[1])
    shown = 0
    for ch, n_bad in sorted_rej:
        if shown >= args.top_n:
            break
        if n_bad == 0:
            break       # everything else is also zero
        mtype = per_channel_stats[ch]["type"] or "-"
        frac  = 100.0 * n_bad / n_snap
        print(f"    {ch:<12s}  {mtype:<8s}  "
              f"{n_bad:>12d}  {frac:>6.2f}%")
        shown += 1
    if shown == 0:
        print("    (no channels triggered any rejection)")

    # ── 6. optional outputs ────────────────────────────────────────────
    if args.out_events:
        good_events = event_nums[event_stable]
        with open(args.out_events, "w") as f:
            f.write("\n".join(str(int(e)) for e in good_events))
            f.write("\n")
        print(f"\n  → stable event list → {args.out_events} "
              f"({good_events.size:,} events)")

    if args.out_json:
        report = {
            "input": str(args.input),
            "event_tree": evt_tree_name,
            "hycal_map":  str(args.hycal_map),
            "scope":      args.scope,
            "type_name":  args.type_name,
            "channel_name": args.channel_name,
            "abs_threshold": args.abs_threshold,
            "rel_n_mad":     args.rel_n_mad,
            "require_v0set": args.require_v0set,
            "window_s":      args.window_s,
            "interval_ms":   interval_ms,
            "monitored_channels": [c[1] for c in selected],
            "n_snapshots":          n_snap,
            "n_snapshots_stable":   n_stable_snap,
            "n_snapshots_unstable": n_unstable_snap,
            "n_events":          n_events,
            "n_events_stable":   n_evt_stable,
            "n_events_pre_hv":   n_evt_pre_hv,
            "rejection_per_channel": dict(rejection),
            "per_channel_stats":     dict(per_channel_stats),
        }
        with open(args.out_json, "w") as f:
            json.dump(report, f, indent=2, default=lambda o: None)
        print(f"  → JSON report   → {args.out_json}")


if __name__ == "__main__":
    main()
