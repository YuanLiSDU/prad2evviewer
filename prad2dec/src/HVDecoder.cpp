// HVDecoder.cpp — VMDF v2 archive parser, windowed segment, lookup helpers.
//
// The on-disk layout this file walks is owned by prad2hvmon (see
// prad2hvmon/scripts/vmon_reader.py for the canonical Python implementation).
// Treat it as the authoritative reference for record sizes / tags; this
// file mirrors them under HV_FORMAT_*.
//
// Parsing strategy: mmap the file (POSIX) then walk records in two
// conceptual passes.  Pass 1 collects (timestamp, byte-offset) for every
// DV / BOOSTER record whose timestamp falls in the requested window plus
// every CHTABLE / BOOSTER_TABLE in the file (small; needed for V0Set
// reconstruction across boundaries).  Pass 2 copies floats out of the
// mapped pages into freshly-sized result vectors (n_kept_rows × n_kept_cols).

#include "HVDecoder.h"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace hv {
namespace {

// ── on-disk constants ────────────────────────────────────────────────────
// Mirror of prad2hvmon/scripts/vmon_reader.py — keep them lock-step.
constexpr char     HV_MAGIC[4]          = {'V','M','D','2'};
constexpr uint8_t  HV_TAG_CHTABLE       = 0x01;
constexpr uint8_t  HV_TAG_DV            = 0x02;
constexpr uint8_t  HV_TAG_BOOSTER_TABLE = 0x03;
constexpr uint8_t  HV_TAG_BOOSTER       = 0x04;
constexpr int      HV_NAME_LEN          = 12;
constexpr int      HV_CHREC_BYTES       = HV_NAME_LEN + 4;          // name + V0Set
constexpr int      HV_BST_TABLE_BYTES   = HV_NAME_LEN + 4 + 4;      // name + VSet + ISet
constexpr int      HV_BST_SNAP_BYTES    = 4 + 4;                    // VMon + IMon

// HVSegment::save() / load() write VMDF v2 itself (same format the recorder
// produces) — so a saved segment is a valid drop-in for vmon_reader.py and
// any other VMDF-aware tool, just smaller.

// ── little-endian readers (the recorder writes LE) ───────────────────────
inline uint16_t rd_u16(const uint8_t *p) {
    uint16_t v; std::memcpy(&v, p, 2); return v;
}
inline uint32_t rd_u32(const uint8_t *p) {
    uint32_t v; std::memcpy(&v, p, 4); return v;
}
inline int64_t rd_i64(const uint8_t *p) {
    int64_t v; std::memcpy(&v, p, 8); return v;
}
inline float   rd_f32(const uint8_t *p) {
    float v;   std::memcpy(&v, p, 4); return v;
}

// Decode a 12-byte channel name field (NUL-terminated, ASCII).
std::string decode_name(const uint8_t *p)
{
    int n = 0;
    while (n < HV_NAME_LEN && p[n] != 0) ++n;
    return std::string(reinterpret_cast<const char *>(p), n);
}

// RAII wrapper around a POSIX mmap of a file.
class MmapFile {
public:
    explicit MmapFile(const std::string &path)
    {
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0)
            throw std::runtime_error("HVDecoder: open(" + path + "): "
                                     + std::strerror(errno));
        struct stat st;
        if (::fstat(fd_, &st) != 0) {
            ::close(fd_);
            throw std::runtime_error("HVDecoder: fstat(" + path + "): "
                                     + std::strerror(errno));
        }
        size_ = static_cast<std::size_t>(st.st_size);
        if (size_ == 0) {
            ::close(fd_);
            throw std::runtime_error("HVDecoder: empty file " + path);
        }
        void *m = ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (m == MAP_FAILED) {
            ::close(fd_);
            throw std::runtime_error("HVDecoder: mmap(" + path + "): "
                                     + std::strerror(errno));
        }
        // Sequential read pattern: hint the kernel to do read-ahead and
        // drop pages as we move past them, keeping RSS low across one
        // pass over a multi-GB file.
        ::madvise(m, size_, MADV_SEQUENTIAL);
        data_ = static_cast<const uint8_t *>(m);
    }
    ~MmapFile() {
        if (data_) ::munmap(const_cast<uint8_t *>(data_), size_);
        if (fd_ >= 0) ::close(fd_);
    }
    MmapFile(const MmapFile &)            = delete;
    MmapFile &operator=(const MmapFile &) = delete;

    const uint8_t *data() const { return data_; }
    std::size_t    size() const { return size_; }

private:
    int            fd_   = -1;
    const uint8_t *data_ = nullptr;
    std::size_t    size_ = 0;
};

// ── per-file scan result (intermediate; never exposed) ───────────────────
struct FileChunk {
    int                       interval_ms = 0;
    int                       n_ch_in_file = 0;
    int                       n_bst_in_file = 0;

    std::vector<std::string>  all_channels;       // canonical order in this file
    std::vector<std::string>  booster_names;
    std::vector<int>          col_idx_in_file;    // kept_channels[i] → column in file
    std::vector<std::string>  kept_channels;

    std::vector<int64_t>      timestamps_ms;
    std::vector<float>        dv;                 // n_dv × kept_channels.size()
    std::vector<ChEvent>      ch_events;

    std::vector<int64_t>      booster_timestamps_ms;
    std::vector<float>        booster_vmon;       // n_bsnap × n_bst_in_file
    std::vector<float>        booster_imon;
    std::vector<BstEvent>     booster_events;
};

// Read the canonical channel-name list from a CHTABLE record's data area.
std::vector<std::string>
decode_name_table(const uint8_t *p, int n, int stride)
{
    std::vector<std::string> out;
    out.reserve(n);
    for (int i = 0; i < n; ++i)
        out.push_back(decode_name(p + i * stride));
    return out;
}

// Project `src_row` (size == src_cols.size()) into `dst_row` (size ==
// dst_cols.size()) via the [src_idx[c]] mapping; missing → NaN.
inline void project_row(const float *src,
                        const std::vector<int> &src_idx,
                        float *dst, int n_dst)
{
    for (int j = 0; j < n_dst; ++j) {
        int s = src_idx[j];
        dst[j] = (s >= 0) ? src[s] : std::numeric_limits<float>::quiet_NaN();
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Parse one VMDF file, keeping only DV/booster records in [t0_ms, t1_ms),
// projected onto the requested channel filter (empty = keep all).
// ─────────────────────────────────────────────────────────────────────────
FileChunk read_window(const std::string &path,
                      int64_t t0_ms, int64_t t1_ms,
                      const std::vector<std::string> &channel_filter)
{
    MmapFile mm(path);
    const uint8_t *base = mm.data();
    const std::size_t size = mm.size();

    if (size < 20)
        throw std::runtime_error(path + ": file too small");
    if (std::memcmp(base, HV_MAGIC, 4) != 0) {
        if (std::memcmp(base, "VMD1", 4) == 0)
            throw std::runtime_error(path + ": legacy VMD1 not supported");
        throw std::runtime_error(path + ": bad magic");
    }
    uint16_t ver       = rd_u16(base + 4);
    int      n_ch      = rd_u16(base + 6);
    int      interval  = rd_u16(base + 8);
    int      n_bst     = rd_u16(base + 10);
    int64_t  t0_epoch  = rd_i64(base + 12);
    if (ver != 2)
        throw std::runtime_error(path + ": unsupported VMDF version "
                                 + std::to_string(ver));

    const std::size_t chtable_size   = 1 + 4 + std::size_t(n_ch)  * HV_CHREC_BYTES;
    const std::size_t dv_size        = 1 + 4 + std::size_t(n_ch)  * 4;
    const std::size_t bst_table_size = 1 + 4 + std::size_t(n_bst) * HV_BST_TABLE_BYTES;
    const std::size_t bst_snap_size  = 1 + 4 + std::size_t(n_bst) * HV_BST_SNAP_BYTES;

    // Pass 1: walk every record, capturing offsets for things we'll keep.
    std::vector<std::pair<int64_t, std::size_t>> dv_recs;     // (ts_ms, data_off)
    std::vector<std::pair<int64_t, std::size_t>> bst_recs;
    std::vector<std::pair<int64_t, std::size_t>> ch_tbls;     // every CHTABLE
    std::vector<std::pair<int64_t, std::size_t>> bst_tbls;    // every BOOSTER_TABLE

    std::size_t pos = 20;
    while (pos < size) {
        uint8_t tag = base[pos];
        if (tag == HV_TAG_CHTABLE) {
            if (pos + chtable_size > size) break;
            uint32_t dt = rd_u32(base + pos + 1);
            ch_tbls.emplace_back(t0_epoch + int64_t(dt), pos + 5);
            pos += chtable_size;
        }
        else if (tag == HV_TAG_DV) {
            if (pos + dv_size > size) break;
            uint32_t dt = rd_u32(base + pos + 1);
            int64_t  ts = t0_epoch + int64_t(dt);
            if (ts >= t0_ms && ts < t1_ms)
                dv_recs.emplace_back(ts, pos + 5);
            pos += dv_size;
        }
        else if (tag == HV_TAG_BOOSTER_TABLE) {
            if (pos + bst_table_size > size) break;
            uint32_t dt = rd_u32(base + pos + 1);
            bst_tbls.emplace_back(t0_epoch + int64_t(dt), pos + 5);
            pos += bst_table_size;
        }
        else if (tag == HV_TAG_BOOSTER) {
            if (pos + bst_snap_size > size) break;
            uint32_t dt = rd_u32(base + pos + 1);
            int64_t  ts = t0_epoch + int64_t(dt);
            if (ts >= t0_ms && ts < t1_ms)
                bst_recs.emplace_back(ts, pos + 5);
            pos += bst_snap_size;
        }
        else {
            // Unknown tag = corruption / EOF; bail without a hard error.
            break;
        }
    }

    // Resolve channel + booster names from the latest tables.
    std::vector<std::string> all_channels;
    if (!ch_tbls.empty())
        all_channels = decode_name_table(base + ch_tbls.back().second,
                                         n_ch, HV_CHREC_BYTES);
    else {
        all_channels.reserve(n_ch);
        for (int i = 0; i < n_ch; ++i)
            all_channels.push_back("ch" + std::to_string(i));
    }
    std::vector<std::string> booster_names;
    if (!bst_tbls.empty())
        booster_names = decode_name_table(base + bst_tbls.back().second,
                                          n_bst, HV_BST_TABLE_BYTES);
    else {
        booster_names.reserve(n_bst);
        for (int i = 0; i < n_bst; ++i)
            booster_names.push_back("bst" + std::to_string(i));
    }

    // Resolve channel filter to file-local column indices.
    std::vector<int>         col_idx;
    std::vector<std::string> kept;
    if (channel_filter.empty()) {
        col_idx.reserve(n_ch);
        kept = all_channels;
        for (int i = 0; i < n_ch; ++i) col_idx.push_back(i);
    } else {
        for (const auto &cn : channel_filter) {
            auto it = std::find(all_channels.begin(), all_channels.end(), cn);
            if (it == all_channels.end()) continue;     // not in this file
            col_idx.push_back(int(it - all_channels.begin()));
            kept.push_back(cn);
        }
    }

    FileChunk chunk;
    chunk.interval_ms     = interval;
    chunk.n_ch_in_file    = n_ch;
    chunk.n_bst_in_file   = n_bst;
    chunk.all_channels    = std::move(all_channels);
    chunk.booster_names   = std::move(booster_names);
    chunk.col_idx_in_file = col_idx;
    chunk.kept_channels   = kept;

    // ── Pass 2: copy / project ──────────────────────────────────────────
    const int n_dv = int(dv_recs.size());
    const int n_keep_ch = int(col_idx.size());
    chunk.timestamps_ms.resize(n_dv);
    chunk.dv.resize(std::size_t(n_dv) * n_keep_ch);

    for (int i = 0; i < n_dv; ++i) {
        chunk.timestamps_ms[i] = dv_recs[i].first;
        const float *row_src = reinterpret_cast<const float *>(
                                   base + dv_recs[i].second);
        float *row_dst = chunk.dv.data() + std::size_t(i) * n_keep_ch;
        if (n_keep_ch == n_ch) {
            std::memcpy(row_dst, row_src, std::size_t(n_ch) * 4);
        } else {
            for (int j = 0; j < n_keep_ch; ++j)
                row_dst[j] = row_src[col_idx[j]];
        }
    }

    // CHTABLE history (always small; project columns).
    chunk.ch_events.reserve(ch_tbls.size());
    for (auto &kv : ch_tbls) {
        ChEvent ev;
        ev.abs_ts_ms = kv.first;
        ev.v0sets.resize(n_keep_ch);
        const uint8_t *tbl = base + kv.second;
        for (int j = 0; j < n_keep_ch; ++j) {
            int s = col_idx[j];
            ev.v0sets[j] = rd_f32(tbl + s * HV_CHREC_BYTES + HV_NAME_LEN);
        }
        chunk.ch_events.push_back(std::move(ev));
    }

    // Booster snapshots (n_bst is always small, < 10 in practice; keep all
    // booster columns and let the merge sort it out).
    const int n_bsnap = int(bst_recs.size());
    chunk.booster_timestamps_ms.resize(n_bsnap);
    chunk.booster_vmon.resize(std::size_t(n_bsnap) * n_bst);
    chunk.booster_imon.resize(std::size_t(n_bsnap) * n_bst);
    for (int i = 0; i < n_bsnap; ++i) {
        chunk.booster_timestamps_ms[i] = bst_recs[i].first;
        const uint8_t *p = base + bst_recs[i].second;
        for (int j = 0; j < n_bst; ++j) {
            chunk.booster_vmon[std::size_t(i) * n_bst + j] =
                rd_f32(p + j * HV_BST_SNAP_BYTES);
            chunk.booster_imon[std::size_t(i) * n_bst + j] =
                rd_f32(p + j * HV_BST_SNAP_BYTES + 4);
        }
    }

    // BOOSTER_TABLE history.
    chunk.booster_events.reserve(bst_tbls.size());
    for (auto &kv : bst_tbls) {
        BstEvent ev;
        ev.abs_ts_ms = kv.first;
        ev.vsets.resize(n_bst);
        ev.isets.resize(n_bst);
        const uint8_t *tbl = base + kv.second;
        for (int j = 0; j < n_bst; ++j) {
            ev.vsets[j] = rd_f32(tbl + j * HV_BST_TABLE_BYTES + HV_NAME_LEN);
            ev.isets[j] = rd_f32(tbl + j * HV_BST_TABLE_BYTES + HV_NAME_LEN + 4);
        }
        chunk.booster_events.push_back(std::move(ev));
    }
    return chunk;
}

// ─────────────────────────────────────────────────────────────────────────
// Merge per-file chunks into one segment.
// ─────────────────────────────────────────────────────────────────────────
HVSegment merge_chunks(const std::vector<FileChunk> &chunks_in,
                       const std::vector<std::string> &channel_filter)
{
    HVSegment seg;
    if (chunks_in.empty()) return seg;

    // Sort chunks chronologically (use first DV ts; fall back to first booster).
    std::vector<const FileChunk *> chunks;
    chunks.reserve(chunks_in.size());
    for (const auto &c : chunks_in) chunks.push_back(&c);
    std::sort(chunks.begin(), chunks.end(),
              [](const FileChunk *a, const FileChunk *b) {
                  int64_t ta = !a->timestamps_ms.empty() ? a->timestamps_ms.front()
                             : !a->booster_timestamps_ms.empty()
                                    ? a->booster_timestamps_ms.front() : 0;
                  int64_t tb = !b->timestamps_ms.empty() ? b->timestamps_ms.front()
                             : !b->booster_timestamps_ms.empty()
                                    ? b->booster_timestamps_ms.front() : 0;
                  return ta < tb;
              });

    // Canonical channel order.
    std::vector<std::string> canon;
    if (!channel_filter.empty()) {
        canon = channel_filter;
    } else {
        std::set<std::string> seen;
        for (const auto *c : chunks)
            for (const auto &n : c->kept_channels)
                if (seen.insert(n).second) canon.push_back(n);
    }
    seg.channels = canon;
    const int n_canon = int(canon.size());

    // Booster canonical order = union (preserving first-seen).
    std::vector<std::string> bst_canon;
    {
        std::set<std::string> seen;
        for (const auto *c : chunks)
            for (const auto &n : c->booster_names)
                if (seen.insert(n).second) bst_canon.push_back(n);
    }
    seg.booster_names = bst_canon;
    const int n_bcanon = int(bst_canon.size());

    seg.interval_ms = chunks.front()->interval_ms;

    // For each chunk, build its column-projection map onto canonical order.
    auto build_proj = [](const std::vector<std::string> &src,
                         const std::vector<std::string> &dst) {
        std::vector<int> map(dst.size(), -1);
        for (std::size_t j = 0; j < dst.size(); ++j) {
            auto it = std::find(src.begin(), src.end(), dst[j]);
            if (it != src.end()) map[j] = int(it - src.begin());
        }
        return map;
    };

    std::vector<std::vector<int>> ch_proj(chunks.size());
    std::vector<std::vector<int>> bst_proj(chunks.size());
    for (std::size_t i = 0; i < chunks.size(); ++i) {
        ch_proj [i] = build_proj(chunks[i]->kept_channels, canon);
        bst_proj[i] = build_proj(chunks[i]->booster_names, bst_canon);
    }

    // Pre-size totals.
    std::size_t tot_dv = 0, tot_ch_ev = 0, tot_bsnap = 0, tot_bst_ev = 0;
    for (const auto *c : chunks) {
        tot_dv     += c->timestamps_ms.size();
        tot_ch_ev  += c->ch_events.size();
        tot_bsnap  += c->booster_timestamps_ms.size();
        tot_bst_ev += c->booster_events.size();
    }
    seg.timestamps_ms.reserve(tot_dv);
    seg.dv.resize(tot_dv * n_canon);
    seg.ch_events.reserve(tot_ch_ev);
    seg.booster_timestamps_ms.reserve(tot_bsnap);
    seg.booster_vmon.resize(tot_bsnap * n_bcanon);
    seg.booster_imon.resize(tot_bsnap * n_bcanon);
    seg.booster_events.reserve(tot_bst_ev);

    std::size_t row_out = 0;
    std::size_t bsnap_out = 0;
    for (std::size_t ci = 0; ci < chunks.size(); ++ci) {
        const FileChunk &c = *chunks[ci];
        const int n_src_ch = int(c.kept_channels.size());

        // dv rows
        for (std::size_t i = 0; i < c.timestamps_ms.size(); ++i, ++row_out) {
            seg.timestamps_ms.push_back(c.timestamps_ms[i]);
            const float *src = c.dv.data() + i * n_src_ch;
            float       *dst = seg.dv.data() + row_out * n_canon;
            project_row(src, ch_proj[ci], dst, n_canon);
        }
        // ch_events
        for (const auto &ev : c.ch_events) {
            ChEvent out;
            out.abs_ts_ms = ev.abs_ts_ms;
            out.v0sets.resize(n_canon);
            project_row(ev.v0sets.data(), ch_proj[ci],
                        out.v0sets.data(), n_canon);
            seg.ch_events.push_back(std::move(out));
        }
        // booster snapshots
        const int n_src_bst = int(c.booster_names.size());
        for (std::size_t i = 0; i < c.booster_timestamps_ms.size();
             ++i, ++bsnap_out)
        {
            seg.booster_timestamps_ms.push_back(c.booster_timestamps_ms[i]);
            const float *vsrc = c.booster_vmon.data() + i * n_src_bst;
            const float *isrc = c.booster_imon.data() + i * n_src_bst;
            float *vdst = seg.booster_vmon.data() + bsnap_out * n_bcanon;
            float *idst = seg.booster_imon.data() + bsnap_out * n_bcanon;
            project_row(vsrc, bst_proj[ci], vdst, n_bcanon);
            project_row(isrc, bst_proj[ci], idst, n_bcanon);
        }
        // booster events
        for (const auto &ev : c.booster_events) {
            BstEvent out;
            out.abs_ts_ms = ev.abs_ts_ms;
            out.vsets.resize(n_bcanon);
            out.isets.resize(n_bcanon);
            project_row(ev.vsets.data(), bst_proj[ci],
                        out.vsets.data(), n_bcanon);
            project_row(ev.isets.data(), bst_proj[ci],
                        out.isets.data(), n_bcanon);
            seg.booster_events.push_back(std::move(out));
        }
    }

    // Sort each block by timestamp (handles overlapping daily files).
    auto reorder_2d = [](std::vector<float> &flat,
                         std::size_t n_rows, std::size_t row_stride,
                         const std::vector<std::size_t> &order)
    {
        std::vector<float> tmp(flat.size());
        for (std::size_t i = 0; i < n_rows; ++i) {
            const float *src = flat.data() + order[i] * row_stride;
            float       *dst = tmp.data()  + i        * row_stride;
            std::memcpy(dst, src, row_stride * sizeof(float));
        }
        flat = std::move(tmp);
    };

    if (!seg.timestamps_ms.empty()) {
        std::vector<std::size_t> order(seg.timestamps_ms.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) {
                             return seg.timestamps_ms[a] < seg.timestamps_ms[b];
                         });
        std::vector<int64_t> ts_sorted(order.size());
        for (std::size_t i = 0; i < order.size(); ++i)
            ts_sorted[i] = seg.timestamps_ms[order[i]];
        seg.timestamps_ms = std::move(ts_sorted);
        reorder_2d(seg.dv, order.size(), n_canon, order);
    }
    if (!seg.ch_events.empty()) {
        std::stable_sort(seg.ch_events.begin(), seg.ch_events.end(),
                         [](const ChEvent &a, const ChEvent &b) {
                             return a.abs_ts_ms < b.abs_ts_ms;
                         });
    }
    if (!seg.booster_timestamps_ms.empty()) {
        std::vector<std::size_t> order(seg.booster_timestamps_ms.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) {
                             return seg.booster_timestamps_ms[a]
                                  < seg.booster_timestamps_ms[b];
                         });
        std::vector<int64_t> ts_sorted(order.size());
        for (std::size_t i = 0; i < order.size(); ++i)
            ts_sorted[i] = seg.booster_timestamps_ms[order[i]];
        seg.booster_timestamps_ms = std::move(ts_sorted);
        reorder_2d(seg.booster_vmon, order.size(), n_bcanon, order);
        reorder_2d(seg.booster_imon, order.size(), n_bcanon, order);
    }
    if (!seg.booster_events.empty()) {
        std::stable_sort(seg.booster_events.begin(), seg.booster_events.end(),
                         [](const BstEvent &a, const BstEvent &b) {
                             return a.abs_ts_ms < b.abs_ts_ms;
                         });
    }
    return seg;
}

// ─────────────────────────────────────────────────────────────────────────
// Centered rolling population-std with NaN-aware count.
// ─────────────────────────────────────────────────────────────────────────
std::vector<float> rolling_std(const std::vector<float> &x, int win)
{
    const int n = int(x.size());
    std::vector<float> out(n, 0.0f);
    if (n == 0 || win <= 1) return out;

    // Use cumulative sums on (x, x*x, valid_mask).
    std::vector<double>   cx (n + 1, 0.0);
    std::vector<double>   cxx(n + 1, 0.0);
    std::vector<int64_t>  cv (n + 1, 0);
    for (int i = 0; i < n; ++i) {
        bool   v   = !std::isnan(x[i]);
        double xi  = v ? double(x[i]) : 0.0;
        cx [i+1] = cx [i] + xi;
        cxx[i+1] = cxx[i] + xi * xi;
        cv [i+1] = cv [i] + (v ? 1 : 0);
    }
    int half = win / 2;
    for (int i = 0; i < n; ++i) {
        int lo = std::max(i - half, 0);
        int hi = std::min(i - half + win, n);
        int64_t cnt = std::max<int64_t>(cv[hi] - cv[lo], 1);
        double  s   = cx [hi] - cx [lo];
        double  s2  = cxx[hi] - cxx[lo];
        double  m   = s / double(cnt);
        double  v   = std::max(s2 / double(cnt) - m * m, 0.0);
        out[i] = float(std::sqrt(v));
    }
    return out;
}

// Pick nearest / next index into `ts` for query `ut_ms`.
int pick_index(const std::vector<int64_t> &ts, int64_t ut_ms, Side side)
{
    if (ts.empty()) return -1;
    auto it = std::lower_bound(ts.begin(), ts.end(), ut_ms);
    int j = int(it - ts.begin());
    if (side == Side::Next) {
        return (j < int(ts.size())) ? j : -1;
    }
    // Nearest
    if (j == 0)              return 0;
    if (j == int(ts.size())) return int(ts.size()) - 1;
    int64_t d_prev = ut_ms - ts[j - 1];
    int64_t d_next = ts[j] - ut_ms;
    return (d_prev <= d_next) ? (j - 1) : j;
}

// ─────────────────────────────────────────────────────────────────────────
// VMDF v2 writer helpers (used by HVSegment::save).
// ─────────────────────────────────────────────────────────────────────────
template <class T>
inline void write_pod(std::ostream &os, const T &v)
{
    os.write(reinterpret_cast<const char *>(&v), sizeof(T));
}

// Write the fixed 12-byte channel-name field used by VMDF.  Match the
// recorder (prad2hvmon/include/vmon_recorder.h::buildCurrentTable):
// at most NAME_LEN-1 = 11 bytes copied, the rest NUL-padded — guarantees
// at least one trailing NUL so split-on-NUL parsers stay safe.
inline void write_name(std::ostream &os, const std::string &name)
{
    char buf[HV_NAME_LEN] = {0};
    std::memcpy(buf, name.data(),
                std::min<std::size_t>(name.size(), HV_NAME_LEN - 1));
    os.write(buf, HV_NAME_LEN);
}

// ─────────────────────────────────────────────────────────────────────────
// Filename helpers (vmon_YYYYMMDD.dat date pruning).
// ─────────────────────────────────────────────────────────────────────────
const std::regex &filename_re()
{
    static const std::regex re(R"(^vmon_(\d{8})\.dat$)",
                               std::regex::icase);
    return re;
}

bool parse_filename_date(const std::string &name, int &y, int &m, int &d)
{
    std::smatch sm;
    if (!std::regex_match(name, sm, filename_re())) return false;
    const std::string &s = sm[1].str();
    y = std::stoi(s.substr(0, 4));
    m = std::stoi(s.substr(4, 2));
    d = std::stoi(s.substr(6, 2));
    return true;
}

// Compare a filename-derived YMD with the UTC date of unix_ms.
// Returns -1 if file_date < unix_date, 0 if equal, +1 if greater.
int compare_filename_date(int fy, int fm, int fd, int64_t unix_ms)
{
    std::time_t t = std::time_t(unix_ms / 1000);
    std::tm utc{};
    gmtime_r(&t, &utc);
    int ty = utc.tm_year + 1900, tm_ = utc.tm_mon + 1, td = utc.tm_mday;
    if (fy != ty) return (fy < ty) ? -1 : +1;
    if (fm != tm_) return (fm < tm_) ? -1 : +1;
    if (fd != td) return (fd < td) ? -1 : +1;
    return 0;
}

} // namespace


// =========================================================================
// HVSegment members
// =========================================================================

int HVSegment::channel_index(const std::string &name) const
{
    auto it = std::find(channels.begin(), channels.end(), name);
    return (it == channels.end()) ? -1 : int(it - channels.begin());
}

int HVSegment::booster_index(const std::string &name) const
{
    auto it = std::find(booster_names.begin(), booster_names.end(), name);
    return (it == booster_names.end()) ? -1 : int(it - booster_names.begin());
}

std::vector<float> HVSegment::v0set_trace(int ch_idx) const
{
    const int n = n_snapshots();
    std::vector<float> out(n, std::numeric_limits<float>::quiet_NaN());
    if (n == 0 || ch_events.empty() ||
        ch_idx < 0 || ch_idx >= n_channels())
        return out;

    std::vector<int64_t> ev_ts;
    ev_ts.reserve(ch_events.size());
    for (const auto &ev : ch_events) ev_ts.push_back(ev.abs_ts_ms);

    for (int i = 0; i < n; ++i) {
        // searchsorted(ts, side=right) - 1
        auto it = std::upper_bound(ev_ts.begin(), ev_ts.end(),
                                   timestamps_ms[i]);
        int idx = int(it - ev_ts.begin()) - 1;
        if (idx < 0) idx = 0;
        if (idx >= int(ev_ts.size())) idx = int(ev_ts.size()) - 1;
        out[i] = ch_events[idx].v0sets[ch_idx];
    }
    return out;
}

LookupResult HVSegment::value_at(const std::string &name,
                                  double unix_time,
                                  Kind kind, Side side) const
{
    LookupResult miss{false, 0.0, 0.0};
    int64_t ut_ms = int64_t(std::llround(unix_time * 1000.0));

    int hv_i = channel_index(name);
    if (hv_i >= 0) {
        if (timestamps_ms.empty()) return miss;
        int j = pick_index(timestamps_ms, ut_ms, side);
        if (j < 0) return miss;
        const int n_ch = n_channels();
        float dv_v = dv[std::size_t(j) * n_ch + hv_i];
        if (std::isnan(dv_v)) return miss;

        float val = 0.f;
        if (kind == Kind::DV) {
            val = dv_v;
        } else if (kind == Kind::V0Set || kind == Kind::VMon) {
            // Look up V0Set at this snapshot's time via ch_events.
            float v0 = std::numeric_limits<float>::quiet_NaN();
            if (!ch_events.empty()) {
                int64_t ts = timestamps_ms[j];
                int lo = 0, hi = int(ch_events.size()) - 1, k = 0;
                // upper_bound + step back
                while (lo <= hi) {
                    int m = (lo + hi) / 2;
                    if (ch_events[m].abs_ts_ms <= ts) { k = m; lo = m + 1; }
                    else { hi = m - 1; }
                }
                v0 = ch_events[k].v0sets[hv_i];
            }
            if (kind == Kind::V0Set) val = v0;
            else                     val = dv_v + v0;
        } else {
            return miss;   // booster kind on HV channel
        }
        if (std::isnan(val)) return miss;
        return LookupResult{true, double(timestamps_ms[j]) / 1000.0, double(val)};
    }

    int b_i = booster_index(name);
    if (b_i < 0) return miss;
    if (booster_timestamps_ms.empty()) return miss;
    int j = pick_index(booster_timestamps_ms, ut_ms, side);
    if (j < 0) return miss;
    const int n_b = n_boosters();
    float val = std::numeric_limits<float>::quiet_NaN();
    if      (kind == Kind::BoosterVMon) val = booster_vmon[std::size_t(j) * n_b + b_i];
    else if (kind == Kind::BoosterIMon) val = booster_imon[std::size_t(j) * n_b + b_i];
    else if (kind == Kind::BoosterVSet || kind == Kind::BoosterISet) {
        if (!booster_events.empty()) {
            int64_t ts = booster_timestamps_ms[j];
            int lo = 0, hi = int(booster_events.size()) - 1, k = 0;
            while (lo <= hi) {
                int m = (lo + hi) / 2;
                if (booster_events[m].abs_ts_ms <= ts) { k = m; lo = m + 1; }
                else { hi = m - 1; }
            }
            val = (kind == Kind::BoosterVSet) ? booster_events[k].vsets[b_i]
                                              : booster_events[k].isets[b_i];
        }
    } else {
        return miss; // HV kind on booster
    }
    if (std::isnan(val)) return miss;
    return LookupResult{true, double(booster_timestamps_ms[j]) / 1000.0,
                        double(val)};
}

std::vector<Interval>
HVSegment::find_stable_intervals(const std::vector<std::string> &cn_check,
                                  double window_s,
                                  double std_threshold,
                                  std::optional<double> dv_threshold,
                                  double min_duration_s,
                                  double guard_s) const
{
    std::vector<Interval> out;
    const int n = n_snapshots();
    if (n == 0) return out;

    // Sample period (ms): use header value when available, else the
    // median diff (handles caches that may not carry interval_ms).
    double dt_ms = 0.0;
    if (interval_ms > 0) {
        dt_ms = interval_ms;
    } else if (n > 1) {
        std::vector<int64_t> diffs;
        diffs.reserve(n - 1);
        for (int i = 1; i < n; ++i)
            diffs.push_back(timestamps_ms[i] - timestamps_ms[i - 1]);
        std::nth_element(diffs.begin(),
                         diffs.begin() + diffs.size() / 2, diffs.end());
        dt_ms = double(diffs[diffs.size() / 2]);
    } else {
        dt_ms = 200.0;
    }
    int win = std::max(1, int(std::round(window_s * 1000.0 / dt_ms)));

    std::vector<bool> unstable(n, false);
    const int n_ch = n_channels();
    for (const auto &cn : cn_check) {
        int ci = channel_index(cn);
        if (ci < 0)
            throw std::runtime_error("HVSegment::find_stable_intervals: "
                                     "channel '" + cn + "' not in segment");
        // Project this channel's column.
        std::vector<float> col(n);
        for (int i = 0; i < n; ++i)
            col[i] = dv[std::size_t(i) * n_ch + ci];
        std::vector<float> rs = rolling_std(col, win);

        for (int i = 0; i < n; ++i) {
            if (std::isnan(col[i])) { unstable[i] = true; continue; }
            if (rs[i] > std_threshold) { unstable[i] = true; continue; }
            if (dv_threshold && std::fabs(col[i]) > *dv_threshold)
                unstable[i] = true;
        }
    }

    int i = 0;
    while (i < n) {
        if (unstable[i]) { ++i; continue; }
        int j = i;
        while (j < n && !unstable[j]) ++j;
        double t0 = double(timestamps_ms[i])     / 1000.0;
        double t1 = double(timestamps_ms[j - 1]) / 1000.0;
        if ((t1 - t0) >= min_duration_s) {
            double a = t0 + guard_s;
            double b = t1 - guard_s;
            if (b > a) out.push_back(Interval{a, b});
        }
        i = j;
    }
    return out;
}

// =========================================================================
// VMDF v2 export — write the segment as a real VMDF file.
//
// What goes where:
//   * 20-byte header: magic "VMD2", version=2, n_channels (= projected
//     count), interval_ms, n_boosters, t0_epoch_ms (chosen as the
//     minimum across all kept records so every record's dt fits in
//     uint32 ms — VMDF's record timestamps are int64-equivalent only
//     via t0+dt).  The format itself does NOT carry an explicit
//     snapshot count — readers walk records to count them.
//   * One CHTABLE per ChEvent.  Each CHTABLE carries the projected
//     channel names + V0Set vector.  When the segment is empty of
//     CHTABLEs (no source CHTABLE was kept) we synthesize one at t0
//     so the reader still gets channel names; V0Set is then NaN.
//   * One DV record per snapshot, with float32 dV per kept channel.
//   * One BOOSTER_TABLE per BstEvent (synth at t0 if empty), one
//     BOOSTER record per booster snapshot.
//
// Records are written in chronological order with stable sort, so a
// CHTABLE update at the same instant as a DV gets emitted first (matches
// the recorder's own convention and lets the reader's `names` reflect
// the new V0Set before the DV is consumed).
// =========================================================================
void HVSegment::save(const std::string &path) const
{
    if (n_channels() == 0 && n_boosters() == 0)
        throw std::runtime_error("HVSegment::save: empty segment "
                                 "(no channels and no boosters)");

    // Pick t0 = min of every kept timestamp (DV / CHTABLE / BST / BST_TABLE).
    // Guarantees uint32 dt is always non-negative and within range as long
    // as the segment span is < 49 days (always true for run-scale windows).
    int64_t t0 = std::numeric_limits<int64_t>::max();
    if (!timestamps_ms.empty())
        t0 = std::min(t0, timestamps_ms.front());
    for (const auto &ev : ch_events)
        t0 = std::min(t0, ev.abs_ts_ms);
    if (!booster_timestamps_ms.empty())
        t0 = std::min(t0, booster_timestamps_ms.front());
    for (const auto &ev : booster_events)
        t0 = std::min(t0, ev.abs_ts_ms);
    if (t0 == std::numeric_limits<int64_t>::max()) t0 = 0;

    int64_t span = 0;
    auto bump_span = [&](int64_t ts) { span = std::max(span, ts - t0); };
    if (!timestamps_ms.empty())         bump_span(timestamps_ms.back());
    for (const auto &ev : ch_events)    bump_span(ev.abs_ts_ms);
    if (!booster_timestamps_ms.empty()) bump_span(booster_timestamps_ms.back());
    for (const auto &ev : booster_events) bump_span(ev.abs_ts_ms);
    if (span > int64_t(0xFFFFFFFFLL))
        throw std::runtime_error("HVSegment::save: span > 49 days; "
                                 "uint32 dt would overflow");

    // If we have no CHTABLE (or BOOSTER_TABLE), synthesize one so the
    // reader can resolve channel names.  V0Set = NaN in that case.
    std::vector<ChEvent> chs = ch_events;
    if (chs.empty() && n_channels() > 0) {
        ChEvent ev;
        ev.abs_ts_ms = t0;
        ev.v0sets.assign(n_channels(),
                         std::numeric_limits<float>::quiet_NaN());
        chs.push_back(std::move(ev));
    }
    std::vector<BstEvent> bsts = booster_events;
    if (bsts.empty() && n_boosters() > 0) {
        BstEvent ev;
        ev.abs_ts_ms = t0;
        ev.vsets.assign(n_boosters(),
                        std::numeric_limits<float>::quiet_NaN());
        ev.isets.assign(n_boosters(),
                        std::numeric_limits<float>::quiet_NaN());
        bsts.push_back(std::move(ev));
    }

    // Build a chronologically-sorted record list.  Stable sort on
    // (ts, type) puts CHTABLE/BOOSTER_TABLE before DV/BOOSTER at the same
    // instant — ensures the reader's `names` is current when the DV
    // arrives.
    enum class Rec : uint8_t {
        ChTable    = 0,
        BstTable   = 1,
        DvSnap     = 2,
        BstSnap    = 3,
    };
    struct Item { int64_t ts; Rec type; std::size_t idx; };
    std::vector<Item> items;
    items.reserve(chs.size() + timestamps_ms.size()
                  + bsts.size() + booster_timestamps_ms.size());
    for (std::size_t i = 0; i < chs.size(); ++i)
        items.push_back({chs[i].abs_ts_ms, Rec::ChTable, i});
    for (std::size_t i = 0; i < bsts.size(); ++i)
        items.push_back({bsts[i].abs_ts_ms, Rec::BstTable, i});
    for (std::size_t i = 0; i < timestamps_ms.size(); ++i)
        items.push_back({timestamps_ms[i], Rec::DvSnap, i});
    for (std::size_t i = 0; i < booster_timestamps_ms.size(); ++i)
        items.push_back({booster_timestamps_ms[i], Rec::BstSnap, i});
    std::stable_sort(items.begin(), items.end(),
                     [](const Item &a, const Item &b) {
                         if (a.ts != b.ts) return a.ts < b.ts;
                         return uint8_t(a.type) < uint8_t(b.type);
                     });

    // ── write file ───────────────────────────────────────────────────────
    std::ofstream os(path, std::ios::binary);
    if (!os) throw std::runtime_error("HVSegment::save: open " + path);

    os.write(HV_MAGIC, 4);
    write_pod<uint16_t>(os, 2);                              // version
    write_pod<uint16_t>(os, uint16_t(n_channels()));         // n_channels
    write_pod<uint16_t>(os, uint16_t(std::max(0, interval_ms)));
    write_pod<uint16_t>(os, uint16_t(n_boosters()));         // n_boosters
    write_pod<int64_t> (os, t0);                             // t0_epoch_ms

    const int n_ch = n_channels();
    const int n_b  = n_boosters();

    for (const auto &it : items) {
        uint32_t dt = uint32_t(it.ts - t0);
        switch (it.type) {
        case Rec::ChTable: {
            os.put(char(HV_TAG_CHTABLE));
            write_pod(os, dt);
            const auto &ev = chs[it.idx];
            for (int j = 0; j < n_ch; ++j) {
                write_name(os, channels[j]);
                write_pod<float>(os, ev.v0sets[j]);
            }
            break;
        }
        case Rec::BstTable: {
            os.put(char(HV_TAG_BOOSTER_TABLE));
            write_pod(os, dt);
            const auto &ev = bsts[it.idx];
            for (int j = 0; j < n_b; ++j) {
                write_name(os, booster_names[j]);
                write_pod<float>(os, ev.vsets[j]);
                write_pod<float>(os, ev.isets[j]);
            }
            break;
        }
        case Rec::DvSnap: {
            os.put(char(HV_TAG_DV));
            write_pod(os, dt);
            const float *row = dv.data() + it.idx * n_ch;
            os.write(reinterpret_cast<const char *>(row),
                     std::size_t(n_ch) * sizeof(float));
            break;
        }
        case Rec::BstSnap: {
            os.put(char(HV_TAG_BOOSTER));
            write_pod(os, dt);
            const std::size_t row_off = it.idx * n_b;
            for (int j = 0; j < n_b; ++j) {
                write_pod<float>(os, booster_vmon[row_off + j]);
                write_pod<float>(os, booster_imon[row_off + j]);
            }
            break;
        }
        }
    }
    if (!os) throw std::runtime_error("HVSegment::save: write failed " + path);
}

// =========================================================================
// VMDF v2 load — return a segment populated from the entire file.
//
// load_window already does the heavy lifting (mmap, walk records, project
// columns, merge across files).  Here we just route a single-file path
// through it with an unbounded time window and no channel filter.
// =========================================================================
HVSegment HVSegment::load(const std::string &path)
{
    if (!std::filesystem::exists(path))
        throw std::runtime_error("HVSegment::load: file not found " + path);

    HVDecoder dec(std::vector<std::string>{path});
    // Wide bounds so every record passes the time filter.  4e10 unix sec ≈
    // year 3237 — plenty.  Internally these get rounded to int64 ms.
    HVSegment seg = dec.load_window(0.0, 4.0e10, /*channels=*/{},
                                    /*cache_path=*/"");
    return seg;
}


// =========================================================================
// HVDecoder members
// =========================================================================

HVDecoder::HVDecoder(const std::string &source)
    : files_(discover_files(source))
{
    if (files_.empty())
        throw std::runtime_error("HVDecoder: no vmon_*.dat files found at "
                                 + source);
}

HVDecoder::HVDecoder(const std::vector<std::string> &files)
    : files_(files)
{
    std::sort(files_.begin(), files_.end());
    if (files_.empty())
        throw std::runtime_error("HVDecoder: empty file list");
}

std::vector<std::string>
HVDecoder::discover_files(const std::string &source)
{
    std::vector<std::string> out;
    namespace fs = std::filesystem;
    fs::path p(source);
    if (fs::is_directory(p)) {
        for (const auto &e : fs::directory_iterator(p)) {
            if (!e.is_regular_file()) continue;
            int y, mo, d;
            if (parse_filename_date(e.path().filename().string(), y, mo, d))
                out.push_back(e.path().string());
        }
        std::sort(out.begin(), out.end());
    } else if (fs::is_regular_file(p)) {
        out.push_back(p.string());
    } else {
        throw std::runtime_error("HVDecoder: not a file or directory: " + source);
    }
    return out;
}

std::vector<std::string>
HVDecoder::candidate_files_for(int64_t t_start_ms, int64_t t_end_ms) const
{
    // Pad the query date range by 1 day on each side; the recorder rotates
    // around local midnight, so a file labelled day D may legitimately
    // carry data from D-1 (late evening) through D+0.
    int64_t pad_ms = 86400LL * 1000;
    int64_t lo = t_start_ms - pad_ms;
    int64_t hi = t_end_ms   + pad_ms;

    std::vector<std::string> out;
    for (const auto &f : files_) {
        std::string name = std::filesystem::path(f).filename().string();
        int y, mo, d;
        if (!parse_filename_date(name, y, mo, d)) {
            out.push_back(f);   // unrecognized → try it
            continue;
        }
        if (compare_filename_date(y, mo, d, lo) >= 0 &&
            compare_filename_date(y, mo, d, hi) <= 0)
            out.push_back(f);
    }
    return out;
}

HVSegment HVDecoder::load_window(double t_start_unix, double t_end_unix,
                                  const std::vector<std::string> &channels,
                                  const std::string &cache_path) const
{
    if (!cache_path.empty() && std::filesystem::exists(cache_path)) {
        return HVSegment::load(cache_path);
    }
    if (t_end_unix <= t_start_unix)
        throw std::runtime_error("HVDecoder::load_window: end <= start");

    int64_t t0 = int64_t(std::llround(t_start_unix * 1000.0));
    int64_t t1 = int64_t(std::llround(t_end_unix   * 1000.0));

    std::vector<FileChunk> chunks;
    std::vector<std::string> used_files;
    for (const auto &f : candidate_files_for(t0, t1)) {
        try {
            FileChunk c = read_window(f, t0, t1, channels);
            if (!c.timestamps_ms.empty() ||
                !c.booster_timestamps_ms.empty()) {
                chunks.push_back(std::move(c));
                used_files.push_back(std::filesystem::path(f).filename().string());
            }
        } catch (const std::exception &e) {
            std::cerr << "[HVDecoder] skipping " << f << ": " << e.what()
                      << std::endl;
        }
    }
    HVSegment seg = merge_chunks(chunks, channels);
    seg.window_start_unix = t_start_unix;
    seg.window_end_unix   = t_end_unix;
    seg.source_files      = std::move(used_files);

    if (!cache_path.empty()) seg.save(cache_path);
    return seg;
}

} // namespace hv
