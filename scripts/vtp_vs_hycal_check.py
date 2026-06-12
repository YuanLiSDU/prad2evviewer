#!/usr/bin/env python3
"""Audit VTP 0xE122 bank record content across a replayed run and compare any
EC_PEAK / EC_CLUSTER energy fields to a HyCal FADC waveform-integral total.

This is a diagnostic for the question: "does the VTP bank carry raw-sum /
cluster energy info that we can sanity-check against HyCal total energy?"
The answer is YES (June 2026): the TAG_EXP records are PRAD_CLUSTER —
trigger-level HyCal clusters, 9-bit tag 0x1CC per $CLON_PARMS/clonbanks.xml:
w0 [13:0] E; w1 [26:15] seed module (bit 11 → W(id&0x7FF), else G(id)),
[14:11] N hits, [10:0] T.  Decoded by prad2dec VtpDecoder into
VtpEventData.prad_clusters (prad2py.dec.PradCluster in Python); validated
vs offline clustering at 0.99 energy correlation on run 24340.  The
CLAS12-style EC_PEAK / EC_CLUSTER records this script also looks for are
never emitted by the PRad-II firmware.

Usage:
    python3 scripts/vtp_vs_hycal_check.py test_data/prad_024555*_raw.root
"""

import argparse
import sys
from collections import Counter
from glob import glob

import ROOT


# ---------- VTP record-type decoder (matches prad2dec/src/VtpDecoder.cpp) ----

T_BLKHDR, T_BLKTLR, T_EVTHDR, T_TRGTIME = 0x10, 0x11, 0x12, 0x13
T_EC_PEAK, T_EC_CLUSTER = 0x14, 0x15
T_TAG_EXP, T_TRIGGER, T_DNV, T_FILLER = 0x1C, 0x1D, 0x1E, 0x1F

TYPE_NAMES = {
    T_BLKHDR: "BLKHDR", T_BLKTLR: "BLKTLR",
    T_EVTHDR: "EVTHDR", T_TRGTIME: "TRGTIME",
    T_EC_PEAK: "EC_PEAK", T_EC_CLUSTER: "EC_CLUSTER",
    0x16: "HTCC", 0x17: "FT", 0x18: "FTOF", 0x19: "CTOF",
    0x1A: "CND", 0x1B: "PCU",
    T_TAG_EXP: "TAG_EXP", T_TRIGGER: "TRIGGER",
    T_DNV: "DNV", T_FILLER: "FILLER",
}


def decode_vtp_bank(words):
    """Single-pass replica of VtpDecoder::DecodeRoc.

    Returns:
        record_type_counts: Counter of defining-word type codes
        peaks: list of dicts {inst, view, time, coord, energy} from EC_PEAK
        clusters: list of dicts {inst, time, energy, coordU, coordV, coordW}
    """
    counts = Counter()
    peaks = []
    clusters = []

    cur_type = 0
    cur_cont = 0
    peak_inst = peak_view = peak_time = 0
    cl_inst = cl_time = 0
    cl_energy = 0

    for w in words:
        defining = (w >> 31) & 1
        if defining:
            cur_type = (w >> 27) & 0x1F
            cur_cont = 0
            counts[cur_type] += 1
            if cur_type == T_EC_PEAK:
                peak_inst = (w >> 26) & 1
                peak_view = (w >> 24) & 3
                peak_time = (w >> 16) & 0xFF
            elif cur_type == T_EC_CLUSTER:
                cl_inst = (w >> 26) & 1
                cl_time = (w >> 16) & 0xFF
                cl_energy = w & 0xFFFF
        else:
            cur_cont += 1
            if cur_type == T_EC_PEAK and cur_cont == 1:
                peaks.append(dict(
                    inst=peak_inst, view=peak_view, time=peak_time,
                    coord=(w >> 16) & 0x3FF, energy=w & 0xFFFF))
            elif cur_type == T_EC_CLUSTER and cur_cont == 1:
                clusters.append(dict(
                    inst=cl_inst, time=cl_time, energy=cl_energy,
                    coordW=(w >> 20) & 0x3FF,
                    coordV=(w >> 10) & 0x3FF,
                    coordU=w & 0x3FF))
    return counts, peaks, clusters


# ---------- HyCal proxy energy from raw FADC samples ------------------------

# Run 24555 uses FADC W_OFFSET = 2950 ns (one of the late TET overrides in the
# config) with a soft window we'll fix here.  Without -p we have no firmware
# pedestal — we estimate a per-channel pedestal from the first samples of each
# waveform.
PED_NSAMP = 4           # use first 4 samples for pedestal
PULSE_LO, PULSE_HI = 10, 60   # sample range to integrate for a rough peak

# Default ADC->MeV factor when per-module calibration is not loaded.
# (database/runinfo/general.json defaults.calibration.default_adc2mev = 0.12)
DEFAULT_ADC2MEV = 0.12


def hycal_proxy_energy(samples_array, nsamples_array, module_type_array,
                       gain_array):
    """Σ over HyCal modules of (integral - pedestal_baseline) * gain * adc2mev.

    samples_array: shape (nch, 200) uint16
    nsamples_array: shape (nch,) uint8
    module_type_array: shape (nch,) — only types 1 (PbGlass) and 2 (PbWO4)
    gain_array: shape (nch,) — LMS-relative gain correction
    """
    import numpy as np
    nch = len(module_type_array)
    if nch == 0:
        return 0.0
    samples = np.asarray(samples_array, dtype=np.float32).reshape(nch, -1)
    nsamp = np.asarray(nsamples_array, dtype=np.int32)
    mtype = np.asarray(module_type_array, dtype=np.int32)
    gain = np.asarray(gain_array, dtype=np.float32)

    mask = (mtype == 1) | (mtype == 2)   # HyCal modules only
    if not mask.any():
        return 0.0

    samples = samples[mask]
    nsamp = nsamp[mask]
    gain = gain[mask]

    ped = samples[:, :PED_NSAMP].mean(axis=1)
    # Integrate samples in [PULSE_LO, min(PULSE_HI, nsamp)] minus pedestal.
    hi = np.minimum(nsamp, PULSE_HI)
    integral = np.zeros(samples.shape[0], dtype=np.float32)
    for i in range(samples.shape[0]):
        a, b = PULSE_LO, int(hi[i])
        if b <= a:
            continue
        integral[i] = samples[i, a:b].sum() - ped[i] * (b - a)
    integral = np.clip(integral, 0, None)
    adc_after_gain = integral * gain
    energy_mev = adc_after_gain * DEFAULT_ADC2MEV
    return float(energy_mev.sum())


# ---------- top-level driver -------------------------------------------------

def audit_files(file_patterns, max_events_per_file=None, hycal_for_events=20):
    files = []
    for p in file_patterns:
        files += sorted(glob(p))
    if not files:
        sys.exit(f"no files match: {file_patterns}")

    # Per-ROC record-type tally, plus collections of energy values + flag
    # showing whether any non-stub (EC_PEAK / EC_CLUSTER) records exist.
    roc_type_counts = {}          # roc -> Counter(type -> count)
    roc_event_counts = Counter()  # roc -> # events with that ROC present
    n_events_total = 0
    n_events_with_ec = 0
    fat_events = []               # (file, idx, ev_num, vtp_nwords_total)

    # Cache for the EC-bearing events: also collect HyCal proxy energy.
    sample_records = []   # list of dicts with VTP+HyCal sums

    for fp in files:
        f = ROOT.TFile.Open(fp)
        if not f or f.IsZombie():
            print(f"  [skip] cannot open {fp}", file=sys.stderr)
            continue
        t = f.Get("events")
        nentries = t.GetEntries()
        if max_events_per_file is not None:
            nentries = min(nentries, max_events_per_file)

        for i in range(nentries):
            t.GetEntry(i)
            n_events_total += 1
            roc_tags = list(t.vtp_roc_tags)
            nwords = list(t.vtp_nwords)
            words = list(t.vtp_words)

            off = 0
            ev_has_ec = False
            ev_ec_peak_sum = 0.0
            ev_ec_cluster_sum = 0.0
            ev_n_peaks = 0
            ev_n_clusters = 0

            for bi, roc in enumerate(roc_tags):
                nw = nwords[bi]
                block = words[off:off + nw]
                off += nw

                roc_event_counts[roc] += 1
                counts, peaks, clusters = decode_vtp_bank(block)
                if roc not in roc_type_counts:
                    roc_type_counts[roc] = Counter()
                roc_type_counts[roc].update(counts)

                if peaks:
                    ev_has_ec = True
                    ev_n_peaks += len(peaks)
                    ev_ec_peak_sum += sum(p["energy"] for p in peaks)
                if clusters:
                    ev_has_ec = True
                    ev_n_clusters += len(clusters)
                    ev_ec_cluster_sum += sum(c["energy"] for c in clusters)

            if ev_has_ec:
                n_events_with_ec += 1
                if len(sample_records) < hycal_for_events:
                    # Pull HyCal proxy for this event.
                    nch = getattr(t, "hycal.nch")
                    try:
                        import numpy as np
                        samples = np.frombuffer(
                            getattr(t, "hycal.samples"), dtype=np.uint16)
                        # samples is flat (kMaxChannels*200) — but
                        # we only use first nch*200
                    except Exception:
                        samples = None
                    # Easier: just use ROOT array conversion via list.
                    sm = list(getattr(t, "hycal.samples"))
                    ns = list(getattr(t, "hycal.nsamples"))
                    mt = list(getattr(t, "hycal.module_type"))
                    gn = list(getattr(t, "hycal.gain_factor"))
                    e_hycal = hycal_proxy_energy(sm, ns, mt, gn)
                    sample_records.append(dict(
                        file=fp, idx=i,
                        event_num=int(t.event_num),
                        trigger_bits=int(t.trigger_bits),
                        n_peaks=ev_n_peaks,
                        n_clusters=ev_n_clusters,
                        vtp_peak_sum=ev_ec_peak_sum,
                        vtp_cluster_sum=ev_ec_cluster_sum,
                        hycal_proxy_mev=e_hycal,
                        hycal_nch=nch,
                    ))

            if sum(nwords) > 30:
                fat_events.append((fp, i, int(t.event_num), sum(nwords)))

        f.Close()

    # ---------- Report ---------------------------------------------------
    print("=" * 78)
    print(f"Scanned {len(files)} files, {n_events_total:,} physics events")
    print("=" * 78)
    print()
    print(f"Events with any EC_PEAK or EC_CLUSTER record: {n_events_with_ec}")
    print(f"Events with >30 total VTP words (fat VTP):    {len(fat_events)}")
    print()

    print("Per-ROC record-type tally")
    print("-" * 78)
    header = f"{'ROC':>6}  {'events':>10}  " + "  ".join(
        f"{n:>9}" for n in (
            "BLKHDR", "BLKTLR", "EVTHDR", "TRGTIME",
            "EC_PEAK", "EC_CLUSTER", "TAG_EXP", "TRIGGER", "DNV", "FILLER"))
    print(header)
    for roc in sorted(roc_type_counts.keys()):
        ct = roc_type_counts[roc]
        row = [f"0x{roc:04x}", f"{roc_event_counts[roc]:>10,}"]
        for tc in (T_BLKHDR, T_BLKTLR, T_EVTHDR, T_TRGTIME,
                   T_EC_PEAK, T_EC_CLUSTER, T_TAG_EXP, T_TRIGGER,
                   T_DNV, T_FILLER):
            row.append(f"{ct.get(tc, 0):>9,}")
        print("  ".join(row))
    print()

    if n_events_with_ec == 0:
        print("=> No EC_PEAK or EC_CLUSTER records were shipped by any VTP ROC.")
        print("   VTP banks carry trigger metadata only (EVTHDR / TRGTIME /")
        print("   TAG_EXP / TRIGGER); the trigger-time raw-sum and cluster-sum")
        print("   that the firmware computes are NOT written to the readout —")
        print("   only the SSP front-panel trigger bits derived from them.")
        print("   Comparison with HyCal total energy is therefore not possible")
        print("   from this run's VTP banks.")
    else:
        print(f"=> Sample of {len(sample_records)} EC-bearing events:")
        for r in sample_records:
            print(f"   ev#{r['event_num']:>8}  "
                  f"trig=0x{r['trigger_bits']:08x}  "
                  f"VTP n_peaks={r['n_peaks']:>3}  n_cl={r['n_clusters']:>3}  "
                  f"Σ_peak={r['vtp_peak_sum']:>8.1f}  "
                  f"Σ_cluster={r['vtp_cluster_sum']:>8.1f}  "
                  f"HyCal_proxy={r['hycal_proxy_mev']:>9.1f} MeV  "
                  f"hycal.nch={r['hycal_nch']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="paths or globs to *_raw.root")
    ap.add_argument("-n", "--max-events-per-file", type=int, default=None)
    ap.add_argument("--hycal-samples", type=int, default=20,
                    help="how many EC-bearing events to compute HyCal proxy for")
    args = ap.parse_args()
    audit_files(args.files, args.max_events_per_file, args.hycal_samples)


if __name__ == "__main__":
    main()
