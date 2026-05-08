#pragma once
//=============================================================================
// HVDecoder.h — Reader for VMDF v2 high-voltage archive files.
//
// VMDF v2 is the binary format produced by prad2hvd
// (prad2hvmon/scripts/vmon_reader.py owns the canonical Python reference).
// The on-disk encoding is delta-based:
//
//   * One CHTABLE record names every channel and stores its V0Set.
//   * Each DV record holds (VMon - V0Set) for every channel at one snapshot.
//   * Booster lanes (TDK-Lambda) live in parallel BOOSTER_TABLE / BOOSTER
//     records that carry already-absolute VMon/IMon.
//   * Reconstructed VMon = dV + V0Set(latest CHTABLE at-or-before snapshot).
//
// What this class adds:
//   * mmap-based parser that allocates only the requested rows × columns,
//     so a few-channel / few-minute query against the multi-GB daily archive
//     stays bounded in memory.
//   * (channel name, unix_time) → nearest / nearest_next lookup.
//   * Stable-interval finder that excludes HV ramps / glitches.
//   * Compact binary cache so re-running an analysis avoids re-scanning
//     the daily files.
//
// Use through HVDecoder (multi-file) or — once you've loaded a window —
// through the resulting HVSegment.  Both are wrapped by the prad2py.dec
// Python bindings (see python/bind_dec.cpp::bind_hv).
//=============================================================================

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace hv
{

// -------------------------------------------------------------------------
// CHTABLE / BOOSTER_TABLE history events.  Stored projected to the kept
// channel set so V0Set / VSet can be reconstructed at any timestamp via
// the same column ordering as `dv` / `booster_vmon`.
// -------------------------------------------------------------------------
struct ChEvent {
    int64_t            abs_ts_ms = 0;
    std::vector<float> v0sets;       // size == HVSegment::n_channels()
};

struct BstEvent {
    int64_t            abs_ts_ms = 0;
    std::vector<float> vsets;        // size == HVSegment::n_boosters()
    std::vector<float> isets;
};

// -------------------------------------------------------------------------
// Lookup descriptors (used by HVSegment::value_at).
// -------------------------------------------------------------------------
enum class Kind {
    VMon,        // HV: dV + V0Set
    DV,          // HV: VMon - V0Set
    V0Set,       // HV: setpoint at this instant (CHTABLE at-or-before)
    BoosterVMon, // booster
    BoosterIMon,
    BoosterVSet, // booster setpoint (BOOSTER_TABLE at-or-before)
    BoosterISet,
};

enum class Side {
    Nearest,     // closest snapshot in absolute time
    Next,        // first snapshot at-or-after the query
};

struct LookupResult {
    bool   ok;          // false → query was out of range / NaN
    double t_unix_s;    // the picked snapshot's wall-clock time
    double value;
};

struct Interval {
    double t_start_unix;
    double t_end_unix;
};

// -------------------------------------------------------------------------
// HVSegment — one materialized time × channel window.
//
// Layout mirrors VMDF v2: dv is row-major (n_snapshots × n_channels) and
// is held in delta form.  Reconstructed VMon at snapshot i for channel ci
// is `dv[i, ci] + v0set_trace(ci)[i]`.
//
// All numeric arrays are owned by the segment (vectors).  The Python
// bindings expose them as numpy copies — never views onto the C++ vector,
// since the HVSegment may go out of scope while a numpy reference lingers.
// -------------------------------------------------------------------------
class HVSegment
{
public:
    HVSegment() = default;

    // ── HV ──────────────────────────────────────────────────────────────
    std::vector<int64_t>      timestamps_ms;     // n_snap
    std::vector<std::string>  channels;          // n_ch (canonical order)
    std::vector<float>        dv;                // n_snap × n_ch (row-major)
    std::vector<ChEvent>      ch_events;
    int                       interval_ms = 0;

    // ── Booster ─────────────────────────────────────────────────────────
    std::vector<int64_t>      booster_timestamps_ms;
    std::vector<std::string>  booster_names;
    std::vector<float>        booster_vmon;      // n_bsnap × n_bst
    std::vector<float>        booster_imon;      // n_bsnap × n_bst
    std::vector<BstEvent>     booster_events;

    // ── meta (carried in cache; not load-bearing for analysis) ─────────
    std::vector<std::string>  source_files;
    double                    window_start_unix = 0.0;
    double                    window_end_unix   = 0.0;

    // ── shape ───────────────────────────────────────────────────────────
    int  n_channels()           const { return static_cast<int>(channels.size()); }
    int  n_snapshots()          const { return static_cast<int>(timestamps_ms.size()); }
    int  n_boosters()           const { return static_cast<int>(booster_names.size()); }
    int  n_booster_snapshots()  const { return static_cast<int>(booster_timestamps_ms.size()); }
    bool empty()                const { return n_snapshots() == 0 && n_booster_snapshots() == 0; }

    int  channel_index(const std::string &name) const;
    int  booster_index(const std::string &name) const;

    // ── reconstructed traces (returned as fresh vectors) ───────────────
    // V0Set per snapshot for a kept channel (size == n_snapshots()).
    // For HV channels not present in any CHTABLE, returns NaN.
    std::vector<float> v0set_trace(int ch_idx) const;

    // ── (channel, unix_time) lookup ─────────────────────────────────────
    LookupResult value_at(const std::string &name,
                          double unix_time,
                          Kind   kind,
                          Side   side) const;

    LookupResult nearest      (const std::string &name, double unix_time,
                               Kind kind = Kind::VMon) const
    { return value_at(name, unix_time, kind, Side::Nearest); }

    LookupResult nearest_next (const std::string &name, double unix_time,
                               Kind kind = Kind::VMon) const
    { return value_at(name, unix_time, kind, Side::Next); }

    // ── stable-interval finder ──────────────────────────────────────────
    //
    // Returns (t_start, t_end) ranges of the segment in which every named
    // channel is "stable":
    //   * rolling-std(dV) over `window_s` <= std_threshold
    //   * dV is not NaN
    //   * |dV| <= dv_threshold (if dv_threshold has a value)
    //
    // Stable runs shorter than `min_duration_s` are dropped; `guard_s` is
    // trimmed from each end so post-transition settle time is excluded.
    std::vector<Interval> find_stable_intervals(
        const std::vector<std::string> &channels_to_check,
        double window_s        = 5.0,
        double std_threshold   = 0.5,
        std::optional<double> dv_threshold = std::nullopt,
        double min_duration_s  = 5.0,
        double guard_s         = 1.0) const;

    // ── persistence ─────────────────────────────────────────────────────
    //
    // save() writes a real VMDF v2 file (same format the recorder produces),
    // sized to this segment: header.n_channels = projected channel count,
    // header.t0_epoch_ms = min over kept timestamps so every record's
    // uint32 dt is in range.  The result is a drop-in for vmon_reader.py
    // and any other VMDF tool, just smaller.
    //
    // load() pulls the entire file in via the same mmap/projection pipeline
    // load_window() uses (no time filter, no channel filter).  So a saved
    // segment round-trips verbatim, but you can also point load() at any
    // VMDF v2 file (a daily archive, or another segment).
    //
    // Window-meta fields (window_start_unix / window_end_unix /
    // source_files) are NOT persisted — they're derivable from the data
    // (timestamps_ms.front()/back()) or the original analysis context.
    void              save(const std::string &path) const;
    static HVSegment  load(const std::string &path);
};


// -------------------------------------------------------------------------
// HVDecoder — open one or more vmon_*.dat files; load_window() returns a
// segment.  Construction does not read data; load_window() does.
//
// Files named `vmon_YYYYMMDD.dat` are matched against the requested
// window's UTC date ± 1 day (the recorder rotates around local midnight,
// so a file labelled day D may carry data from D-1 through D+0).  Files
// that don't match the naming convention are tried unconditionally.
// -------------------------------------------------------------------------
class HVDecoder
{
public:
    // `source` may be a directory (auto-discover vmon_*.dat) or a single
    // .dat path.  Use the vector overload for a curated subset.
    explicit HVDecoder(const std::string &source);
    explicit HVDecoder(const std::vector<std::string> &files);

    const std::vector<std::string> &files() const { return files_; }

    // Load a windowed segment.  Empty channels list = all channels.
    // Empty cache_path = no caching.  If cache_path exists, the cache is
    // returned verbatim — channel/window args are not re-validated.  Delete
    // the cache file to force a refetch.
    HVSegment load_window(double t_start_unix,
                          double t_end_unix,
                          const std::vector<std::string> &channels = {},
                          const std::string &cache_path = "") const;

private:
    std::vector<std::string> files_;

    // Conservative date-based pruning (UTC date ± 1 day).
    std::vector<std::string>
        candidate_files_for(int64_t t_start_ms, int64_t t_end_ms) const;

    static std::vector<std::string>
        discover_files(const std::string &source);
};

} // namespace hv
