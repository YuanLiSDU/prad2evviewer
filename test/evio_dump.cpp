// test/evio_dump.cpp
// Diagnostic tool to inspect EVIO file structure.
//
// Modes:
//   evio_dump <file>                          -- summary: count events by type/tag
//   evio_dump <file> --tree [--num N]         -- print bank tree for first N events
//   evio_dump <file> --tags                   -- list all unique bank tags with counts
//   evio_dump <file> --epics                  -- dump EPICS text from all EPICS events
//   evio_dump <file> --event N                -- detailed dump of event N (1-based)

#include "EvChannel.h"
#include "EvStruct.h"
#include "Fadc250Data.h"
#include "TdcData.h"
#include "VtpData.h"
#include "load_daq_config.h"
#include "InstallPaths.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <tuple>
#include <getopt.h>
#include <bitset>
#include <cstdlib>
#include <ctime>

using namespace evc;

// --- helpers ----------------------------------------------------------------
static std::string hex(uint32_t v)
{
    char buf[16];
    snprintf(buf, sizeof(buf), "0x%04X", v);
    return buf;
}

static std::string tag_label(uint32_t tag)
{
    // CODA3 physics events (page 27/29)
    if (tag == 0xFF50) return "PHYSICS(PEB)";
    if (tag == 0xFF58) return "PHYSICS(PEB+sync)";
    if (tag == 0xFF70) return "PHYSICS(SEB)";
    if (tag == 0xFF78) return "PHYSICS(SEB+sync)";
    if (tag >= 0xFF50 && tag <= 0xFF8F) return "PHYSICS";

    // CODA3 trigger banks (page 26)
    if (tag >= 0xFF10 && tag <= 0xFF1F) return "TRIGGER(raw)";
    if (tag >= 0xFF20 && tag <= 0xFF2F) return "TRIGGER(built)";
    if (tag == 0xFF4F) return "TRIGGER(bad)";

    // CODA3 control events (page 20/29)
    if (tag == 0xFFD0) return "SYNC";
    if (tag == 0xFFD1) return "PRESTART";
    if (tag == 0xFFD2) return "GO";
    if (tag == 0xFFD3) return "PAUSE";
    if (tag == 0xFFD4) return "END";
    if (tag >= 0xFFD0 && tag <= 0xFFDF) return "CONTROL";

    // Legacy CODA2 tags (may appear in older data)
    if (tag == 0x11) return "PRESTART(legacy)";
    if (tag == 0x12) return "GO(legacy)";
    if (tag == 0x14) return "END(legacy)";
    if (tag == 0xC1) return "SYNC(legacy)";

    // JLab single-event physics
    if ((tag >= 0x00A0 && tag <= 0x00BF) || (tag >= 0xFF50 && tag <= 0xFF8F)) return "PHYSICS(built)";

    // JLab-specific banks
    if (tag == 0xC000) return "TRIGGER_BANK";

    // EPICS
    if (tag == 0x1F) return "EPICS";

    return "";
}

// --- mode: summary ----------------------------------------------------------
static int doSummary(EvChannel &ch)
{
    std::map<uint32_t, int> tag_counts;
    int total = 0;

    while (ch.Read() == status::success) {
        auto hdr = ch.GetEvHeader();
        tag_counts[hdr.tag]++;
        total++;
    }

    std::cout << "=== Event Summary ===\n";
    std::cout << "Total EVIO records: " << total << "\n\n";
    std::cout << std::setw(12) << "Tag" << std::setw(12) << "Count"
              << "  " << "Label" << "\n";
    std::cout << std::string(40, '-') << "\n";

    for (auto &[tag, cnt] : tag_counts) {
        std::cout << std::setw(12) << hex(tag) << std::setw(12) << cnt
                  << "  " << tag_label(tag) << "\n";
    }
    return 0;
}

// --- mode: tree -------------------------------------------------------------
static int doTree(EvChannel &ch, int num)
{
    int count = 0;
    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (++count > num) break;

        auto hdr = ch.GetEvHeader();
        std::cout << "========== Record " << count
                  << "  tag=" << hex(hdr.tag)
                  << " (" << tag_label(hdr.tag) << ")"
                  << "  num=" << hdr.num
                  << "  len=" << hdr.length << "w"
                  << " ==========\n";
        ch.PrintTree(std::cout);
        std::cout << "\n";
    }
    std::cout << "Printed " << std::min(count, num) << " record(s).\n";
    return 0;
}

// --- mode: tags (tree-structured bank-tag hierarchy) ------------------------
// Builds a trie of (root → parent → ... → child) paths across every record
// in the file and prints it depth-first. Because each subtree is keyed by its
// full path from root (not just by parent tag), the same data bank under
// different ROC parents shows independent counts.
//
// Line format per node:
//     <indent><tag>  <TYPE>   <count>/<parent_count> (<ratio>)  <desc>  [<size>]
//
// ratio = child count / parent count
//   1.000 : exactly one child per parent occurrence
//  ~1.086 : most events have one, some have multiple
//   0.000 : rare optional bank

struct TrieNode {
    uint32_t tag  = 0;
    uint32_t type = 0;
    int      count = 0;
    size_t   min_words = SIZE_MAX;
    size_t   max_words = 0;
    std::map<uint32_t, TrieNode> children;
};

// Informative names for bank tags known by JLab ROLs but not held in
// DaqConfig (e.g. no decoder yet).  Reference: docs/rols/clonbanks_20260406.xml.
static const char *known_bank_name(uint32_t tag)
{
    switch (tag) {
    case 0xE122: return "VTP Hardware Data";
    case 0xE10B: return "V1190/V1290 Hardware Data";
    case 0xE141: return "FAV3 Hardware Data";
    case 0xE104: return "VSCM Hardware Data";
    case 0xE105: return "DCRB Hardware Data";
    case 0xE115: return "DSC2 Scalers";
    case 0xE112: return "HEAD bank (raw)";
    case 0xE123: return "SSP-RICH";
    case 0xE125: return "SIS3801 Scalers";
    case 0xE131: return "VFTDC";
    case 0xE133: return "Helicity Decoder";
    case 0xE140: return "MPD raw (reserved)";
    case 0xE114: return "EPICS string data";
    default:     return nullptr;
    }
}

// Depth-aware description of a bank tag (event / ROC / data-bank / composite-inner).
static std::string tag_description(uint32_t tag, int depth, const DaqConfig &cfg)
{
    // depth 0: event bank (trigger/class identifier)
    if (depth == 0) {
        std::string lbl = tag_label(tag);
        if (!lbl.empty()) return lbl;
        // Recognize any tag in [physics_base, physics_base+0x7F] as a physics
        // trigger even if not enumerated in physics_tags — e.g. 0x80 = trigger 0.
        if (cfg.is_physics(tag) ||
            (tag >= cfg.physics_base && tag <= cfg.physics_base + 0x7F))
        {
            char buf[48];
            snprintf(buf, sizeof(buf), "PHYSICS(trg=0x%02X)", tag - cfg.physics_base);
            return buf;
        }
        if (cfg.is_monitoring(tag)) return "MONITORING";
        if (cfg.is_epics(tag))      return "EPICS";
        if (cfg.is_control(tag))    return "CONTROL";
    }

    // depth 1: ROC crate or depth-1 data bank (e.g. 0xC000 trigger bank)
    if (depth == 1) {
        if (tag == cfg.trigger_bank_tag) return "CODA trigger bank";
        if (tag == cfg.ti_master_tag)    return "ti_master";
        for (auto &re : cfg.roc_tags) {
            if (re.tag == tag) {
                std::string s = re.name;
                if (!re.type.empty()) s += " [" + re.type + "]";
                return s;
            }
        }
    }

    // depth 2: data banks inside ROCs
    if (depth == 2) {
        if (tag == cfg.ti_bank_tag)        return "TI Hardware Data";
        if (tag == cfg.trigger_bank_tag)   return "CODA trigger bank";
        if (tag == cfg.fadc_composite_tag) return "FADC250 composite";
        if (tag == cfg.fadc_raw_tag)       return "FADC250 raw";
        if (tag == cfg.adc1881m_bank_tag)  return "ADC1881M raw";
        if (cfg.is_ssp_bank(tag))          return "SSP/MPD data";
        if (tag == cfg.run_info_tag)       return "Run info";
        if (tag == cfg.daq_config_tag)     return "DAQ config string";
        if (tag == cfg.tdc_bank_tag)       return "V1190 TDC data";
        if (tag == cfg.epics_bank_tag)     return "EPICS data";
        if (auto *k = known_bank_name(tag)) return k;
    }

    // depth 3: composite internals
    if (depth == 3) {
        if (tag == 0x000D) return "composite metadata";
        if (tag == 0x0000) return "composite data words";
    }

    return "";
}

static void printTrie(const TrieNode &node, int indent_level, int parent_count,
                      int depth, const DaqConfig &cfg)
{
    // Sort children by tag (std::map already sorts, but sort a vector view for
    // flexibility if we later want to sort by count).
    std::vector<const TrieNode*> ordered;
    ordered.reserve(node.children.size());
    for (auto &kv : node.children) ordered.push_back(&kv.second);

    std::string indent(indent_level * 2, ' ');
    for (auto *c : ordered) {
        std::cout << indent << hex(c->tag)
                  << "  " << std::setw(9) << std::left << TypeName(c->type)
                  << std::right;

        char cbuf[96];
        double ratio = parent_count > 0 ? static_cast<double>(c->count) / parent_count : 0.0;
        snprintf(cbuf, sizeof(cbuf), "  %9d/%-9d (%.3f)", c->count, parent_count, ratio);
        std::cout << cbuf;

        std::string desc = tag_description(c->tag, depth, cfg);
        if (!desc.empty()) std::cout << "  " << desc;

        std::cout << "  [";
        if (c->min_words == c->max_words) std::cout << c->min_words << "w";
        else std::cout << c->min_words << "-" << c->max_words << "w";
        std::cout << "]\n";

        printTrie(*c, indent_level + 1, c->count, depth + 1, cfg);
    }
}

static int doTags(EvChannel &ch)
{
    TrieNode root;
    int nrecords = 0;
    std::vector<TrieNode*> trie_pos;   // reused across events

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        nrecords++;

        const auto &nodes = ch.GetNodes();
        trie_pos.assign(nodes.size(), nullptr);

        // Each node in the scanned bank tree contributes +1 to exactly one
        // trie position, determined by its full path from root. Because
        // nodes[] is DFS-ordered, every node's parent has already been
        // visited and its trie position cached in trie_pos[].
        for (size_t ni = 0; ni < nodes.size(); ++ni) {
            const auto &n = nodes[ni];
            TrieNode *parent_trie = (n.parent >= 0) ? trie_pos[n.parent] : &root;

            auto &child = parent_trie->children[n.tag];
            child.tag  = n.tag;
            child.type = n.type;
            child.count++;
            child.min_words = std::min(child.min_words, n.data_words);
            child.max_words = std::max(child.max_words, n.data_words);

            trie_pos[ni] = &child;
        }
    }

    std::cout << "=== Bank Tag Hierarchy (across " << nrecords << " records) ===\n\n";
    std::cout << "Legend: <tag>  <TYPE>  <count>/<parent_count> (ratio)  <description>  [<size>]\n";
    std::cout << "        ratio = child count / parent count"
                 " (1.000 = exactly-one-per-parent)\n\n";

    printTrie(root, 0, nrecords, 0, ch.GetConfig());
    return 0;
}

// --- mode: epics ------------------------------------------------------------
// In PRad-II these are SYNC events wrapping EPICS + an absolute-time / live-
// time HEAD bank (0xE112 in this DAQ config).  Dumps:
//   - bank tree of each record
//   - raw words + CODA-style decode of every non-string leaf bank
//   - first N bytes of each string bank (EPICS payload)
//
// `max_events=0` means unlimited; otherwise stop after that many records.
static int doEpics(EvChannel &ch, int max_events)
{
    int count = 0, record = 0;

    while (ch.Read() == status::success) {
        record++;
        if (max_events > 0 && count >= max_events) break;
        if (!ch.Scan()) continue;

        auto hdr = ch.GetEvHeader();
        // check for EPICS: common tags are 0x1F, but also scan for string banks
        bool is_epics = (hdr.tag == 0x1F || hdr.tag == 0x1f);

        if (!is_epics) {
            // also check if any child bank has string data
            for (auto &n : ch.GetNodes()) {
                if (n.depth == 1 &&
                    (n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) &&
                    n.data_words > 4)
                {
                    is_epics = true;
                    break;
                }
            }
        }

        if (!is_epics) continue;

        count++;
        std::cout << "--- SYNC/EPICS record " << count
                  << " (file record " << record
                  << ", tag=" << hex(hdr.tag)
                  << ", num=" << hdr.num << ") ---\n";

        // print tree structure
        ch.PrintTree(std::cout);

        // Dump every non-string leaf bank inside this SYNC event.  Layout is
        // DAQ-specific (the PRad-II 0xE112 "HEAD" bank has 5 words with an
        // inner header in d[0..1], an absolute unix time around d[3], and
        // live-time counters in the remaining words) — rather than commit to
        // a guess, print each word as decimal + hex and flag any slot whose
        // value falls in the 2001-2100 unix-time range.
        std::cout << "\n  Non-string data banks (raw + per-word decode):\n";
        for (auto &n : ch.GetNodes()) {
            if (n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) continue;
            if (IsContainer(n.type)) continue;
            if (n.data_words == 0) continue;

            const uint32_t *d = ch.GetData(n);
            std::cout << "    tag=" << hex(n.tag)
                      << " type=" << TypeName(n.type)
                      << " words=" << n.data_words << "\n";

            // raw hex dump (up to 12 words — enough for 5-word HEAD bank).
            std::cout << "      raw:";
            size_t nshow = std::min<size_t>(n.data_words, 12);
            for (size_t i = 0; i < nshow; ++i)
                std::cout << " " << std::hex << std::setw(8) << std::setfill('0')
                          << d[i] << std::setfill(' ') << std::dec;
            if (n.data_words > nshow) std::cout << " ...";
            std::cout << "\n";

            // Plausible unix-time window: 2001-01-01 .. 2100-01-01.
            constexpr uint32_t UNIX_MIN = 978307200u;    // 2001-01-01 UTC
            constexpr uint32_t UNIX_MAX = 4102444800u;   // 2100-01-01 UTC

            for (size_t i = 0; i < n.data_words; ++i) {
                std::cout << "      d[" << i << "]="
                          << std::setw(12) << d[i]
                          << "  (0x" << std::hex << std::setw(8)
                          << std::setfill('0') << d[i]
                          << std::setfill(' ') << std::dec << ")";
                if (d[i] >= UNIX_MIN && d[i] <= UNIX_MAX) {
                    time_t t = static_cast<time_t>(d[i]);
                    char tbuf[64] = "";
                    std::tm *gm = std::gmtime(&t);
                    if (gm) std::strftime(tbuf, sizeof(tbuf),
                                          "%Y-%m-%d %H:%M:%S UTC", gm);
                    std::cout << "  [unix time? " << tbuf << "]";
                }
                std::cout << "\n";
            }
        }

        // extract and print all string data
        for (auto &n : ch.GetNodes()) {
            if ((n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) && n.data_words > 0) {
                const char *raw = reinterpret_cast<const char*>(ch.GetData(n));
                size_t max_len = n.data_words * 4;
                size_t len = 0;
                while (len < max_len && raw[len] != '\0') ++len;

                std::cout << "\n  [String data from tag=" << hex(n.tag)
                          << ", " << len << " bytes]:\n";

                // first ~800 chars — enough to see several channel lines
                // without flooding the terminal.  Full payload is at offset 0
                // in the bank if callers want to re-dump it.
                size_t show = std::min(len, size_t(800));
                std::string text(raw, show);
                std::cout << text;
                if (show < len) std::cout << "\n  ... (" << len - show << " more bytes)";
                std::cout << "\n";
            }
        }
        std::cout << "\n";
    }

    std::cout << "Printed " << count << " SYNC/EPICS record(s)"
              << " (scanned " << record << " total records)";
    if (max_events > 0 && count >= max_events)
        std::cout << "; stopped at -n limit";
    std::cout << ".\n";
    return 0;
}

// --- mode: vtp --------------------------------------------------------------
// Walks physics events, dumps the raw 0xE122 VTP banks word-by-word with
// per-word record-type decoding (BLKHDR/BLKTLR/EVTHDR/TRGTIME/EC_PEAK/
// EC_CLUSTER/PRAD_CLUSTER/PRAD_TRIGGER/...), and runs ch.Vtp() to show what
// the current decoder turns each bank into.  The raw view is the ground
// truth — the PRad-II VTP firmware is a CLAS12 derivative; trigger-level
// HyCal clusters ship as PRAD_CLUSTER (a TAG_EXP expansion, 9-bit tag
// bits[31:23] = 0x1CC per $CLON_PARMS/clonbanks.xml), while PRAD_TRIGGER
// 0x1D still ships through the decoder unstored.
//
// `max_events=0` falls back to 5 events.

namespace {
const char *vtp_record_name(uint8_t t)
{
    switch (t) {
    case 0x10: return "BLKHDR";
    case 0x11: return "BLKTLR";
    case 0x12: return "EVTHDR";
    case 0x13: return "TRGTIME";
    case 0x14: return "EC_PEAK";
    case 0x15: return "EC_CLUSTER";
    case 0x16: return "HTCC";
    case 0x17: return "FT";
    case 0x18: return "FTOF";
    case 0x19: return "CTOF";
    case 0x1A: return "CND";
    case 0x1B: return "PCU";
    case 0x1C: return "TAG_EXP";
    case 0x1D: return "TRIGGER";    // PRad summary
    case 0x1E: return "DNV";
    case 0x1F: return "FILLER";
    default:   return "?";
    }
}
} // anon

static int doVtp(EvChannel &ch, int max_events)
{
    if (max_events <= 0) max_events = 5;

    int physics_seen = 0, physics_with_vtp = 0, dumped = 0;
    int total_blocks = 0, total_peaks = 0, total_clusters = 0;
    int total_prad_clusters = 0;
    std::map<uint8_t, int> rectype_counts;       // across all dumped events

    while (ch.Read() == status::success) {
        if (dumped >= max_events) break;
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;
        ++physics_seen;

        const int n_sub = ch.GetNEvents();
        for (int ie = 0; ie < n_sub && dumped < max_events; ++ie) {
            ch.SelectEvent(ie);

            // Collect the 0xE122 banks for this sub-event (parented by a
            // ROC bank — there are typically 7 VTP ROCs in PRad-II runs).
            std::vector<const EvNode *> vtp_banks;
            for (auto &n : ch.GetNodes()) {
                if (n.tag != 0xE122 || n.data_words == 0) continue;
                if (n.parent >= 0
                    && ch.GetNodes()[n.parent].type == DATA_COMPOSITE) continue;
                vtp_banks.push_back(&n);
            }
            if (vtp_banks.empty()) continue;
            ++physics_with_vtp;
            ++dumped;

            const auto &info = ch.Info();
            std::cout << "=== Physics event #" << info.event_number
                      << " (sub " << ie << "/" << n_sub
                      << ", trigger_type=" << int(info.trigger_type)
                      << ", trigger_bits=0x" << std::hex << info.trigger_bits
                      << std::dec
                      << ", " << vtp_banks.size() << " VTP bank(s)) ===\n";

            for (const EvNode *np : vtp_banks) {
                uint32_t roc = (np->parent >= 0)
                               ? ch.GetNodes()[np->parent].tag : 0;
                const uint32_t *d = ch.GetData(*np);
                std::cout << "  roc=" << hex(roc)
                          << "  words=" << np->data_words << "\n";

                // Walk every word.  Defining words have bit31=1; we print
                // the type name and a short field summary.  Continuation
                // words are printed indented under the most recent
                // defining word so the structure reads top-down.
                uint8_t cur_type = 0;
                bool cur_pradcl = false;  // last TAG_EXP was PRAD_CLUSTER (tag9 0x1CC)
                for (size_t i = 0; i < np->data_words; ++i) {
                    uint32_t w = d[i];
                    bool defining = (w >> 31) & 0x1;
                    std::cout << "    [" << std::setw(2) << i << "] "
                              << std::hex << std::setw(8) << std::setfill('0')
                              << w << std::setfill(' ') << std::dec << "  ";
                    if (defining) {
                        cur_type = (w >> 27) & 0x1F;
                        rectype_counts[cur_type]++;
                        if (cur_type == 0x1C)
                            cur_pradcl = (((w >> 23) & 0x1FF) == 0x1CC);
                        std::cout << vtp_record_name(cur_type)
                                  << " (0x" << std::hex << int(cur_type)
                                  << std::dec << ")";
                        switch (cur_type) {
                        case 0x10:  // BLKHDR
                            std::cout << "  slot=" << ((w >> 22) & 0x1F)
                                      << " mod=" << ((w >> 18) & 0xF)
                                      << " blk#=" << ((w >> 8) & 0x3FF)
                                      << " level=" << (w & 0xFF);
                            break;
                        case 0x11:  // BLKTLR
                            std::cout << "  slot=" << ((w >> 22) & 0x1F)
                                      << " nwords=" << (w & 0x3FFFFF);
                            break;
                        case 0x12:  // EVTHDR
                            std::cout << "  evn=" << (w & 0x7FFFFFF);
                            break;
                        case 0x13:  // TRGTIME (lo 24 bits; hi in next word)
                            std::cout << "  time_lo=" << (w & 0xFFFFFF);
                            break;
                        case 0x14:  // EC_PEAK (defining)
                            std::cout << "  inst=" << ((w >> 26) & 0x1)
                                      << " view=" << ((w >> 24) & 0x3)
                                      << " time=" << ((w >> 16) & 0xFF);
                            break;
                        case 0x15:  // EC_CLUSTER (defining)
                            std::cout << "  inst=" << ((w >> 26) & 0x1)
                                      << " time=" << ((w >> 16) & 0xFF)
                                      << " energy=" << (w & 0xFFFF);
                            break;
                        case 0x1C:  // TAG_EXP — PRAD_CLUSTER if tag9 == 0x1CC
                            if (cur_pradcl)
                                std::cout << "  PRAD_CLUSTER  E="
                                          << (w & 0x3FFF);
                            else
                                std::cout << "  tag9=0x" << std::hex
                                          << ((w >> 23) & 0x1FF) << std::dec;
                            break;
                        case 0x1D:  // PRad TRIGGER summary
                            std::cout << "  payload=0x" << std::hex
                                      << (w & 0x07FFFFFF) << std::dec;
                            break;
                        default:
                            break;
                        }
                    } else {
                        std::cout << "  cont(" << vtp_record_name(cur_type)
                                  << ")";
                        switch (cur_type) {
                        case 0x13:  // TRGTIME hi
                            std::cout << "  time_hi=" << (w & 0xFFFFFF);
                            break;
                        case 0x14:
                            std::cout << "  coord=" << ((w >> 16) & 0x3FF)
                                      << " energy=" << (w & 0xFFFF);
                            break;
                        case 0x15:
                            std::cout << "  U=" << (w & 0x3FF)
                                      << " V=" << ((w >> 10) & 0x3FF)
                                      << " W=" << ((w >> 20) & 0x3FF);
                            break;
                        case 0x1C:  // PRAD_CLUSTER w1: id / nhits / time
                            if (cur_pradcl) {
                                uint32_t id = (w >> 15) & 0xFFF;
                                std::cout << "  module="
                                          << ((id & 0x800) ? "W" : "G")
                                          << ((id & 0x800) ? (id & 0x7FF) : id)
                                          << " N=" << ((w >> 11) & 0xF)
                                          << " T=" << (w & 0x7FF);
                            }
                            break;
                        case 0x1D:
                            std::cout << "  payload=0x" << std::hex
                                      << w << std::dec;
                            break;
                        default:
                            break;
                        }
                    }
                    std::cout << "\n";
                }
            }

            // Show what the current decoder kept.
            const auto &v = ch.Vtp();
            std::cout << "  decoder kept: blocks=" << v.n_blocks
                      << " peaks=" << v.n_peaks
                      << " clusters=" << v.n_clusters
                      << " prad_clusters=" << v.n_prad_clusters << "\n";
            for (int b = 0; b < v.n_blocks; ++b) {
                const auto &bk = v.blocks[b];
                std::cout << "    block[" << b << "] roc=" << hex(bk.roc_tag)
                          << " slot=" << int(bk.slot)
                          << " evn="  << bk.event_number
                          << " trgT=" << bk.trigger_time
                          << " nwords=" << bk.nwords
                          << " trailer=" << (bk.has_trailer ? "y" : "n")
                          << "\n";
                total_blocks++;
            }
            for (int p = 0; p < v.n_peaks; ++p) {
                const auto &pk = v.peaks[p];
                std::cout << "    peak[" << p << "] roc=" << hex(pk.roc_tag)
                          << " inst=" << int(pk.inst)
                          << " view=" << int(pk.view)
                          << " coord=" << pk.coord
                          << " energy=" << pk.energy << "\n";
                total_peaks++;
            }
            for (int c = 0; c < v.n_clusters; ++c) {
                const auto &cl = v.clusters[c];
                std::cout << "    cluster[" << c << "] roc=" << hex(cl.roc_tag)
                          << " UVW=(" << cl.coordU << "," << cl.coordV
                          << "," << cl.coordW << ")"
                          << " energy=" << cl.energy << "\n";
                total_clusters++;
            }
            for (int c = 0; c < v.n_prad_clusters; ++c) {
                const auto &cl = v.prad_clusters[c];
                std::cout << "    prad_cluster[" << c << "] roc="
                          << hex(cl.roc_tag)
                          << " module=" << (cl.is_pbwo4() ? "W" : "G")
                          << cl.module()
                          << " E=" << cl.energy
                          << " N=" << int(cl.nhits)
                          << " T=" << cl.time << "\n";
                total_prad_clusters++;
            }
            std::cout << "\n";
        }
    }

    std::cout << "Scanned " << physics_seen << " physics events ("
              << physics_with_vtp << " with VTP banks); dumped first "
              << dumped << ".\n"
              << "Decoder totals across dumped events: blocks=" << total_blocks
              << " peaks=" << total_peaks
              << " clusters=" << total_clusters
              << " prad_clusters=" << total_prad_clusters << "\n"
              << "Record-type tally (defining words seen):\n";
    for (const auto &kv : rectype_counts) {
        std::cout << "  0x" << std::hex << int(kv.first) << std::dec
                  << " " << vtp_record_name(kv.first)
                  << ": " << kv.second << "\n";
    }
    return 0;
}

// --- mode: rf ---------------------------------------------------------------
// Inspect the RF timing TDC bank delivered by ROC 0x40, bank 0xE107 (V1190
// hit format -- bits 31:27 slot, 26 edge, 25:19 channel, 18:0 TDC value with
// ~24 ps LSB after rol2.c's V1190->V1290 normalization).  In run 24386 there
// are 2 active channels (slot=16, ch=0 and ch=8), each firing roughly five
// hits per event spaced ~131 ns apart -- a divided CEBAF RF reference.
//
// Per-event lines list raw TDC ticks followed by ns.  The summary tracks
// per-channel hit counts and the mean inter-hit interval (1/period gives the
// generator frequency seen by the TDC, useful for sanity-checking cabling).
//
// `max_events=0` falls back to 5; passing -n N caps the printed events while
// the summary still scans the whole file.

static constexpr double TDC_LSB_NS  = tdc::TDC_LSB_NS; // 23.436 ps/tick, single source: TdcData.h
static constexpr uint32_t RF_ROC_TAG = 0x0040;

static int doRf(EvChannel &ch, int max_events)
{
    if (max_events <= 0) max_events = 5;

    int physics_seen = 0, with_rf = 0, dumped = 0;
    std::map<std::pair<uint8_t,uint8_t>, uint64_t> ch_hits;   // (slot,ch) -> count
    std::map<uint8_t, uint64_t> diff_sum;                     // ch -> sum of consecutive diffs
    std::map<uint8_t, uint64_t> diff_n;                       // ch -> diff count
    std::map<uint8_t, uint32_t> diff_min, diff_max;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;
        ++physics_seen;
        ch.SelectEvent(0);

        const auto &t = ch.Tdc();

        // Group hits in this event by (slot, channel) preserving arrival order;
        // V1190 emits hits in time order within a channel, so we don't need to
        // resort before computing inter-hit deltas.
        std::map<std::pair<uint8_t,uint8_t>, std::vector<uint32_t>> per_chan;
        for (int i = 0; i < t.n_hits; ++i) {
            const auto &h = t.hits[i];
            if (h.roc_tag != RF_ROC_TAG) continue;
            per_chan[{h.slot, h.channel}].push_back(h.value);
        }
        if (per_chan.empty()) continue;
        ++with_rf;

        for (auto &kv : per_chan) {
            ch_hits[kv.first] += kv.second.size();
            for (size_t i = 1; i < kv.second.size(); ++i) {
                uint32_t d = kv.second[i] - kv.second[i-1];
                uint8_t c = kv.first.second;
                diff_sum[c] += d;
                diff_n[c]++;
                if (!diff_min.count(c) || d < diff_min[c]) diff_min[c] = d;
                if (!diff_max.count(c) || d > diff_max[c]) diff_max[c] = d;
            }
        }

        if (dumped >= max_events) continue;
        ++dumped;

        const auto &info = ch.Info();
        std::cout << "=== Physics event #" << info.event_number
                  << "  ts=" << info.timestamp
                  << "  trigger_bits=0x" << std::hex << info.trigger_bits
                  << std::dec << " ===\n";
        for (auto &kv : per_chan) {
            std::cout << "  slot=" << int(kv.first.first)
                      << " ch="    << std::setw(3) << int(kv.first.second)
                      << "  nhits=" << kv.second.size() << "  values:";
            for (uint32_t v : kv.second)
                std::cout << "  " << v << "(" << std::fixed
                          << std::setprecision(2) << v * TDC_LSB_NS << "ns)";
            std::cout << "\n";
            if (kv.second.size() >= 2) {
                std::cout << "    deltas:";
                for (size_t i = 1; i < kv.second.size(); ++i) {
                    uint32_t d = kv.second[i] - kv.second[i-1];
                    std::cout << "  " << d << "(" << std::fixed
                              << std::setprecision(2) << d * TDC_LSB_NS << "ns)";
                }
                std::cout << "\n";
            }
        }
        std::cout << std::resetiosflags(std::ios::fixed) << "\n";
    }

    std::cout << "=== RF Time Summary (ROC 0x40, bank 0xE107) ===\n";
    std::cout << "  physics events scanned: " << physics_seen
              << "  (" << with_rf << " carry RF hits, "
              << std::fixed << std::setprecision(2)
              << (physics_seen ? 100.0 * with_rf / physics_seen : 0.0)
              << "%)" << std::resetiosflags(std::ios::fixed) << "\n\n";

    std::cout << "  Per-channel hit rate:\n";
    std::cout << "    " << std::setw(5) << "slot" << std::setw(5) << "ch"
              << std::setw(14) << "total_hits" << std::setw(14) << "hits/event"
              << "\n";
    for (auto &kv : ch_hits)
        std::cout << "    " << std::setw(5) << int(kv.first.first)
                  << std::setw(5) << int(kv.first.second)
                  << std::setw(14) << kv.second
                  << std::setw(14) << std::fixed << std::setprecision(3)
                  << (with_rf ? double(kv.second) / with_rf : 0.0) << "\n";

    std::cout << std::resetiosflags(std::ios::fixed)
              << "\n  Inter-hit interval per channel (LSB="
              << TDC_LSB_NS << " ns):\n";
    std::cout << "    " << std::setw(5) << "ch"
              << std::setw(12) << "mean(ticks)" << std::setw(10) << "mean(ns)"
              << std::setw(12) << "min(ticks)"  << std::setw(12) << "max(ticks)"
              << std::setw(12) << "freq(MHz)"   << std::setw(10) << "N\n";
    for (auto &kv : diff_n) {
        if (kv.second == 0) continue;
        double mean_ticks = double(diff_sum[kv.first]) / kv.second;
        double mean_ns    = mean_ticks * TDC_LSB_NS;
        double freq_mhz   = mean_ns > 0 ? 1000.0 / mean_ns : 0.0;
        std::cout << "    " << std::setw(5) << int(kv.first)
                  << std::setw(12) << std::fixed << std::setprecision(2) << mean_ticks
                  << std::setw(10) << std::setprecision(3) << mean_ns
                  << std::setw(12) << diff_min[kv.first]
                  << std::setw(12) << diff_max[kv.first]
                  << std::setw(12) << std::setprecision(4) << freq_mhz
                  << std::setw(10) << kv.second << "\n";
    }
    return 0;
}

// --- mode: single event detail ----------------------------------------------
static int doEvent(EvChannel &ch, int target)
{
    int record = 0;

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;

        if (record != target) continue;

        auto hdr = ch.GetEvHeader();
        std::cout << "=== Record " << record
                  << "  tag=" << hex(hdr.tag)
                  << " (" << tag_label(hdr.tag) << ")"
                  << "  type=0x" << std::hex << hdr.type << std::dec
                  << "  num=" << hdr.num
                  << "  length=" << hdr.length << "w"
                  << " ===\n\n";

        // full tree
        std::cout << "--- Bank Tree ---\n";
        ch.PrintTree(std::cout);

        // if physics, try decoding
        if (ch.GetNEvents() > 0) {
            std::cout << "\n--- Physics Decode ---\n";
            std::cout << "Sub-events in block: " << ch.GetNEvents() << "\n";

            fdec::EventData evt;
            for (int i = 0; i < ch.GetNEvents(); ++i) {
                if (!ch.DecodeEvent(i, evt)) {
                    std::cout << "  sub-event " << i << ": decode failed\n";
                    continue;
                }

                std::cout << "  sub-event " << i
                          << ": event#=" << evt.info.event_number
                          << " trigger#=" << evt.info.trigger_number
                          << " type=0x" << std::hex << (int)evt.info.trigger_type
                          << " trigger_bits=0x"
                          << evt.info.trigger_bits << std::dec
                          << " timestamp=" << evt.info.timestamp
                          << " run=" << evt.info.run_number
                          << " unix_time=" << evt.info.unix_time
                          << " rocs=" << evt.nrocs << "\n";

                for (int r = 0; r < evt.nrocs; ++r) {
                    auto &roc = evt.rocs[r];
                    if (!roc.present) continue;
                    std::cout << "    ROC tag=" << hex(roc.tag)
                              << " slots=" << roc.nslots << "\n";

                    for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                        auto &slot = roc.slots[s];
                        if (!slot.present) continue;

                        int nch = 0;
                        for (int c = 0; c < fdec::MAX_CHANNELS; ++c)
                            if (slot.channel_mask & (1ull << c)) nch++;

                        std::cout << "      slot=" << std::setw(2) << s
                                  << " trigger=" << slot.trigger
                                  << " timestamp=" << slot.timestamp
                                  << " channels=" << nch << " [";

                        for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                            if (!(slot.channel_mask & (1ull << c))) continue;
                            std::cout << c << "(" << slot.channels[c].nsamples << ")";
                            if (c < fdec::MAX_CHANNELS - 1) std::cout << " ";
                        }
                        std::cout << "]\n";
                    }
                }
            }
        }

        // dump raw words of any non-container leaf banks
        std::cout << "\n--- Leaf Bank Data (first 16 words) ---\n";
        for (auto &n : ch.GetNodes()) {
            if (n.child_count > 0 || n.data_words == 0) continue;
            if (IsContainer(n.type) || n.type == DATA_COMPOSITE) continue;

            std::cout << "  tag=" << hex(n.tag)
                      << " type=" << TypeName(n.type)
                      << " depth=" << n.depth
                      << " words=" << n.data_words << " |";

            const uint32_t *d = ch.GetData(n);
            int show = std::min<int>(n.data_words, 16);
            for (int i = 0; i < show; ++i)
                std::cout << " " << std::hex << std::setw(8) << std::setfill('0')
                          << d[i] << std::setfill(' ') << std::dec;
            if (n.data_words > 16) std::cout << " ...";
            std::cout << "\n";
        }

        return 0;
    }

    std::cerr << "Record " << target << " not found (file has " << record << " records).\n";
    return 1;
}

// --- mode: trig-debug -------------------------------------------------------
// Cross-correlates three trigger layers:
//   1. Event tag (top-level) = 0x80 + TI_event_type
//   2. TI event_type (d[0] bits 31:24) = TS trigger decision
//   3. FP trigger bits (TI master d[5]) = raw 32-bit front panel input snapshot
//
// FP bit assignments (from prad_v0.trg):
//   Bits 8-15:  SSP PRAD TRGBIT 0-7  (P2 outputs)
//   Bit 23:     v1495 OR from SD/FADC
//   Bit 24:     v1495 LMS
//   Bit 25:     v1495 alpha
//   Bit 26:     v1495 Faraday
//   Bit 27:     v1495 Master OR
//
// TS monitors mask: 0x0F00FF00 (bits 8-15 and 24-27)

struct TrigDebugEntry {
    uint32_t event_tag;         // top-level bank tag
    uint32_t ti_event_type;     // d[0] >> 24 (from any TI bank)
    uint32_t ti_nwords;         // nwords from TI event header
    uint32_t fp_trigger_bits;   // d[5] from TI master (0 if unavailable)
    int32_t  event_number;      // from 0xC000 or TI d[1]
    int      nrocs;             // number of FADC composite banks found
    bool     has_ti_master;     // TI master bank found
    bool     tag_matches;       // event_tag == 0x80 + ti_event_type
};

static const uint32_t TS_FP_MASK = 0x0F00FF00;

struct FpBitInfo {
    int bit;
    uint32_t mask;
    const char *name;
};

static const FpBitInfo fp_bits[] = {
    {  8, 0x00000100, "SSP TRGBIT0 (RawSum>1000)" },
    {  9, 0x00000200, "SSP TRGBIT1 (1clus>1GeV)"  },
    { 10, 0x00000400, "SSP TRGBIT2 (2clus>1GeV)"  },
    { 11, 0x00000800, "SSP TRGBIT3 (3clus>1GeV)"  },
    { 12, 0x00001000, "SSP TRGBIT4 (disabled)"     },
    { 13, 0x00002000, "SSP TRGBIT5 (disabled)"     },
    { 14, 0x00004000, "SSP TRGBIT6 (disabled)"     },
    { 15, 0x00008000, "SSP TRGBIT7 (100Hz pulser)" },
    { 23, 0x00800000, "v1495 OR from SD/FADC"      },
    { 24, 0x01000000, "v1495 LMS"                   },
    { 25, 0x02000000, "v1495 alpha"                 },
    { 26, 0x04000000, "v1495 Faraday"               },
    { 27, 0x08000000, "v1495 Master OR"             },
};
static const int N_FP_BITS = sizeof(fp_bits) / sizeof(fp_bits[0]);

static int doTrigDebug(EvChannel &ch, bool verbose)
{
    const auto &cfg = ch.GetConfig();
    std::vector<TrigDebugEntry> entries;

    // per-tag accumulators
    struct TagStats {
        int count = 0;
        int has_fadc = 0;
        uint32_t fp_or = 0;         // OR of all FP bits seen with this tag
        uint32_t fp_and = 0xFFFFFFFF; // AND of all FP bits
        int fp_bit_counts[32] = {};
    };
    std::map<uint32_t, TagStats> tag_stats;

    int record = 0, physics = 0;

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;
        if (ch.GetNEvents() == 0) continue;

        auto hdr = ch.GetEvHeader();
        auto &nodes = ch.GetNodes();

        for (int iev = 0; iev < ch.GetNEvents(); ++iev) {
            TrigDebugEntry e = {};
            e.event_tag = hdr.tag;

            // --- extract from 0xC000 trigger bank ---
            for (auto &n : nodes) {
                if (n.tag == cfg.trigger_bank_tag && n.type == DATA_UINT32 && n.data_words >= 1) {
                    e.event_number = static_cast<int32_t>(ch.GetData(n)[0]);
                    break;
                }
            }

            // --- extract TI event_type from first 0xE10A bank ---
            for (auto &n : nodes) {
                if (n.tag == cfg.ti_bank_tag && n.type == DATA_UINT32 && n.data_words >= 1) {
                    const uint32_t *d = ch.GetData(n);
                    e.ti_event_type = (d[0] >> 24) & 0xFF;
                    e.ti_nwords = d[0] & 0xFFFF;
                    if (n.data_words >= 2 && e.event_number == 0)
                        e.event_number = static_cast<int32_t>(d[1]);
                    break;
                }
            }

            // --- extract FP trigger bits from TI master's 0xE10A ---
            for (auto &n : nodes) {
                if (n.depth == 1 && n.tag == cfg.ti_master_tag) {
                    e.has_ti_master = true;
                    for (size_t ci = 0; ci < n.child_count; ++ci) {
                        auto &child = nodes[n.child_first + ci];
                        if (child.tag == cfg.ti_bank_tag && child.type == DATA_UINT32) {
                            const uint32_t *d = ch.GetData(child);
                            size_t nw = child.data_words;
                            // d[4] = 8-bit trigger type byte
                            // d[5] = 32-bit FP trigger inputs (if FP readout enabled)
                            if (nw > 5)
                                e.fp_trigger_bits = d[5];
                            break;
                        }
                    }
                    break;
                }
            }

            // --- count FADC composite banks ---
            for (auto &n : nodes) {
                if (n.tag == cfg.fadc_composite_tag && n.type == DATA_COMPOSITE)
                    e.nrocs++;
            }

            // --- verify tag = 0x80 + event_type ---
            e.tag_matches = (e.event_tag == (0x80u + e.ti_event_type));

            // --- accumulate ---
            physics++;
            entries.push_back(e);

            auto &ts = tag_stats[e.event_tag];
            ts.count++;
            if (e.nrocs > 0) ts.has_fadc++;
            ts.fp_or |= e.fp_trigger_bits;
            ts.fp_and &= e.fp_trigger_bits;
            for (int b = 0; b < 32; ++b)
                if (e.fp_trigger_bits & (1u << b)) ts.fp_bit_counts[b]++;
        }
    }

    // === Per-event detail (-v) ===
    if (verbose) {
        std::cout << std::setw(8) << "event#"
                  << std::setw(8) << "tag"
                  << std::setw(8) << "TItype"
                  << std::setw(6) << "chk"
                  << std::setw(12) << "FP_bits"
                  << std::setw(12) << "FP_masked"
                  << std::setw(6) << "FADC"
                  << "  active FP signals"
                  << "\n";
        std::cout << std::string(90, '-') << "\n";

        for (auto &e : entries) {
            uint32_t masked = e.fp_trigger_bits & TS_FP_MASK;
            std::cout << std::setw(8) << e.event_number
                      << "  0x" << std::hex << std::setw(4) << std::setfill('0') << e.event_tag
                      << "    0x" << std::setw(2) << e.ti_event_type
                      << std::setfill(' ') << std::dec
                      << std::setw(6) << (e.tag_matches ? "OK" : "FAIL")
                      << "  0x" << std::hex << std::setw(8) << std::setfill('0') << e.fp_trigger_bits
                      << "  0x" << std::setw(8) << masked
                      << std::setfill(' ') << std::dec
                      << std::setw(6) << e.nrocs;

            // list active FP bit names
            std::cout << "  ";
            bool first = true;
            for (int i = 0; i < N_FP_BITS; ++i) {
                if (e.fp_trigger_bits & fp_bits[i].mask) {
                    if (!first) std::cout << ", ";
                    std::cout << "b" << fp_bits[i].bit;
                    first = false;
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // === Per-tag summary ===
    std::cout << "=== Trigger Debug Summary (" << physics << " physics events, "
              << record << " records) ===\n\n";

    std::cout << "--- Event Tag Summary ---\n";
    std::cout << std::setw(8) << "tag"
              << std::setw(8) << "TItype"
              << std::setw(8) << "count"
              << std::setw(8) << "w/FADC"
              << std::setw(12) << "FP_OR"
              << std::setw(12) << "FP_AND"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    int mismatch_total = 0;
    for (auto &[tag, ts] : tag_stats) {
        uint32_t expected_type = tag - 0x80;
        std::cout << "  0x" << std::hex << std::setw(4) << std::setfill('0') << tag
                  << "    0x" << std::setw(2) << expected_type
                  << std::setfill(' ') << std::dec
                  << std::setw(8) << ts.count
                  << std::setw(8) << ts.has_fadc
                  << "  0x" << std::hex << std::setw(8) << std::setfill('0') << ts.fp_or
                  << "  0x" << std::setw(8) << ts.fp_and
                  << std::setfill(' ') << std::dec
                  << "\n";
    }

    // === Tag verification ===
    for (auto &e : entries)
        if (!e.tag_matches) mismatch_total++;

    std::cout << "\n--- Tag Verification ---\n";
    std::cout << "  tag == 0x80 + TI_event_type: "
              << (physics - mismatch_total) << " OK, "
              << mismatch_total << " MISMATCH\n";
    if (mismatch_total > 0) {
        std::cout << "  First mismatches:\n";
        int shown = 0;
        for (auto &e : entries) {
            if (!e.tag_matches && shown < 10) {
                std::cout << "    event#=" << e.event_number
                          << " tag=0x" << std::hex << e.event_tag
                          << " TI_type=0x" << e.ti_event_type
                          << " expected_tag=0x" << (0x80 + e.ti_event_type)
                          << std::dec << "\n";
                shown++;
            }
        }
    }

    // === FP bit activity per tag ===
    std::cout << "\n--- FP Bit Activity (TS mask = 0x"
              << std::hex << std::setw(8) << std::setfill('0') << TS_FP_MASK
              << std::setfill(' ') << std::dec << ") ---\n";

    // header
    std::cout << std::setw(5) << "bit" << std::setw(12) << "hex"
              << "  " << std::setw(30) << std::left << "name" << std::right
              << std::setw(8) << "total";
    for (auto &[tag, ts] : tag_stats)
        std::cout << std::setw(8) << ("0x" + hex(tag).substr(2));
    std::cout << "\n";
    std::cout << std::string(55 + 8 * tag_stats.size(), '-') << "\n";

    for (int i = 0; i < N_FP_BITS; ++i) {
        int bit = fp_bits[i].bit;
        bool in_mask = (TS_FP_MASK & fp_bits[i].mask) != 0;

        // count total across all tags
        int total = 0;
        for (auto &[tag, ts] : tag_stats)
            total += ts.fp_bit_counts[bit];

        if (total == 0 && !in_mask) continue; // skip unused bits not in mask

        std::cout << std::setw(5) << bit
                  << "  0x" << std::hex << std::setw(8) << std::setfill('0')
                  << fp_bits[i].mask << std::setfill(' ') << std::dec
                  << (in_mask ? "* " : "  ")
                  << std::setw(30) << std::left << fp_bits[i].name << std::right
                  << std::setw(8) << total;

        for (auto &[tag, ts] : tag_stats)
            std::cout << std::setw(8) << ts.fp_bit_counts[bit];
        std::cout << "\n";
    }
    std::cout << "\n  * = in TS_FP_INPUT_MASK\n";

    // === Check d[4] vs d[0] consistency ===
    std::cout << "\n--- TI Event Header vs d[4] ---\n";
    int d4_match = 0, d4_mismatch = 0, d4_unavail = 0;
    for (auto &e : entries) {
        if (e.ti_nwords <= 3) { d4_unavail++; continue; }
        // We don't have d[4] stored; but we can verify TI type from d[0]
        // matches event_tag. Already done above.
    }
    std::cout << "  (d[4] check requires direct bank access — use -m event -n N for individual inspection)\n";
    std::cout << "  TI event header nwords distribution:\n";
    std::map<uint32_t, int> nwords_dist;
    for (auto &e : entries) nwords_dist[e.ti_nwords]++;
    for (auto &[nw, cnt] : nwords_dist)
        std::cout << "    nwords=" << nw << ": " << cnt << " events\n";

    return 0;
}

// --- mode: triggers ---------------------------------------------------------
// Uses the lazy Info() accessor so we skip Fadc250/SSP/VTP/TDC decoding — we
// only need the TI/trigger-bank metadata.  Typically 5-10× faster than a full
// DecodeEvent() on a 1.9M-event run.
//
// Per-trigger-bits counts are tracked per event type so EPICS / Sync / other
// non-physics events (which have no 0xE10A TI bank and therefore default to
// trigger_bits=0) don't get lumped in with real physics triggers.
//
// Rate estimates use the TI timestamp span across physics events (250 MHz /
// 4 ns ticks) as the wall-clock denominator — the same span is applied to
// non-physics rows since they share the run.
static int doTriggers(EvChannel &ch, bool verbose)
{
    static constexpr double TI_TICK_SEC = 4e-9;

    int record = 0, decoded = 0;
    // key: (event_type, trigger_bits) → count
    std::map<std::pair<EventType, uint32_t>, int> trig_counts;

    // Physics-event TI timestamps drive the run span; non-physics events have
    // info.timestamp == 0 so they don't contribute.
    uint64_t first_ts = 0, last_ts = 0;

    if (verbose) {
        std::cout << std::setw(10) << "evtype"
                  << std::setw(8) << "event#"
                  << std::setw(10) << "trigger#"
                  << std::setw(14) << "trigger_bits"
                  << std::setw(18) << "timestamp"
                  << "\n";
        std::cout << std::string(60, '-') << "\n";
    }

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;
        if (ch.GetNEvents() == 0) continue;

        EventType et = ch.GetEventType();
        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ch.SelectEvent(i);
            const auto &info = ch.Info();
            decoded++;
            trig_counts[{et, info.trigger_bits}]++;

            if (et == EventType::Physics && info.timestamp != 0) {
                if (first_ts == 0) first_ts = info.timestamp;
                last_ts = info.timestamp;
            }

            if (verbose) {
                const char *et_name =
                    et == EventType::Physics ? "Physics" :
                    et == EventType::Epics   ? "Epics"   :
                    et == EventType::Sync    ? "Sync"    :
                    et == EventType::Unknown ? "Unknown" : "Other";
                std::cout << std::setw(10) << et_name
                          << std::setw(8) << info.event_number
                          << std::setw(10) << info.trigger_number
                          << "    0x" << std::hex << std::setw(8)
                          << std::setfill('0') << info.trigger_bits
                          << std::dec << std::setfill(' ')
                          << std::setw(18) << info.timestamp
                          << "\n";
            }
        }
    }

    double duration_sec = (last_ts > first_ts)
        ? static_cast<double>(last_ts - first_ts) * TI_TICK_SEC
        : 0.0;

    std::cout << "=== Trigger Bits Summary (" << decoded << " events";
    if (duration_sec > 0)
        std::cout << ", span " << std::fixed << std::setprecision(3)
                  << duration_sec << " s";
    std::cout << ") ===\n";

    std::cout << "  " << std::setw(10) << std::left << "evtype"
              << std::setw(14) << "trigger_bits"
              << std::setw(10) << "count"
              << std::setw(14) << "rate(Hz)" << std::right << "\n";
    std::cout << "  " << std::string(48, '-') << "\n";
    for (auto &[key, cnt] : trig_counts) {
        const char *et_name =
            key.first == EventType::Physics ? "Physics" :
            key.first == EventType::Epics   ? "Epics"   :
            key.first == EventType::Sync    ? "Sync"    :
            key.first == EventType::Unknown ? "Unknown" : "Other";
        std::cout << "  " << std::setw(10) << std::left << et_name << std::right
                  << "  0x" << std::hex << std::setw(8) << std::setfill('0')
                  << key.second << std::dec << std::setfill(' ') << "  "
                  << std::setw(10) << cnt;
        if (duration_sec > 0) {
            char rbuf[32];
            snprintf(rbuf, sizeof(rbuf), "%.3f", cnt / duration_sec);
            std::cout << std::setw(14) << rbuf;
        } else {
            std::cout << std::setw(14) << "--";
        }
        std::cout << "\n";
    }
    std::cout << "\nNotes:\n"
                 "  - Rates use the TI-timestamp span (250 MHz / 4 ns ticks)\n"
                 "    of physics events as the wall-clock denominator; non-\n"
                 "    physics rows share the same denominator.  '--' means no\n"
                 "    decoded physics events with a valid timestamp were found.\n"
                 "  - Non-physics rows (Epics/Sync/Unknown) have no 0xE10A TI\n"
                 "    bank, so their trigger_bits default to 0.  Inspect their\n"
                 "    bank structure with `-m tree -n <N>` or `-m epics`.\n";

    return 0;
}

// --- mode: bank-debug -------------------------------------------------------
// Walks the event tree the same way DecodeEvent does, verifying the dispatch
// logic: which banks appear at which depth, container vs flat layout, and
// whether every bank is recognized by the decoder.

struct BankStat {
    int count = 0;
    size_t min_words = SIZE_MAX, max_words = 0;
    std::set<uint32_t> parent_tags;
    std::set<int> depths;
};

// Known data-bank tags (mirrors EvChannel.cpp table + infrastructure banks)
struct BankTagInfo {
    uint32_t    tag;
    const char *name;
    const char *status;   // "decoded", "skipped", "info"
};
static const BankTagInfo bank_tag_info[] = {
    { 0xC000, "CODA trigger bank",    "decoded" },
    { 0xE10A, "TI/TS",                "decoded" },
    { 0xE101, "FADC250 composite",    "decoded" },
    { 0xE109, "FADC250 raw",          "decoded" },
    { 0xE120, "ADC1881M (Fastbus)",   "decoded" },
    { 0xE10C, "SSP",                  "decoded" },
    { 0x0DE9, "VTP/MPD (GEM)",        "decoded" },
    { 0xE10F, "Run info",             "decoded" },
    { 0xE10E, "DAQ config string",    "info"    },
    { 0xE10B, "V1190/V1290 TDC",      "no decoder" },
    { 0xE141, "FAV3 (FADC v3)",       "no decoder" },
    { 0xE104, "VSCM",                 "no decoder" },
    { 0xE105, "DCRB/DC/Vetroc",       "no decoder" },
    { 0xE115, "DSC2 scaler",          "no decoder" },
    { 0xE112, "HEAD bank",            "no decoder" },
    { 0xE123, "SSP-RICH",             "no decoder" },
    { 0xE125, "Per-slot data",        "no decoder" },
    { 0xE131, "VFTDC",                "no decoder" },
    { 0xE133, "Helicity Decoder",     "no decoder" },
    { 0xE140, "Special (pid=0)",      "no decoder" },
};

static const BankTagInfo *findBankInfo(uint32_t tag)
{
    for (auto &b : bank_tag_info)
        if (b.tag == tag) return &b;
    return nullptr;
}

static int doBankDebug(EvChannel &ch, bool verbose)
{
    const auto &cfg = ch.GetConfig();

    // --- per-depth, per-tag stats ---
    // key: (depth << 16 | is_container << 15 | 0) for containers, (depth << 16 | tag) for leaf
    std::map<uint32_t, BankStat> depth0_tags;   // event-level tags
    std::map<uint32_t, BankStat> depth1_tags;   // ROC or flat data banks
    std::map<uint64_t, BankStat> depth2_tags;   // data banks keyed by (parent_tag << 32 | tag)

    int nrecords = 0, nphysics = 0;
    int n_built = 0, n_flat = 0, n_mixed = 0;  // event structure types

    while (ch.Read() == status::success) {
        nrecords++;
        if (!ch.Scan()) continue;
        if (ch.GetNEvents() == 0) continue;
        nphysics++;

        auto &nodes = ch.GetNodes();
        if (nodes.empty()) continue;
        auto &ev = nodes[0];

        // track event tag
        auto &es = depth0_tags[ev.tag];
        es.count++;
        es.min_words = std::min(es.min_words, ev.data_words);
        es.max_words = std::max(es.max_words, ev.data_words);
        es.depths.insert(0);

        // classify event structure: built (all depth-1 children are containers)
        // vs flat (some depth-1 children are data banks)
        bool has_container = false, has_leaf = false;
        for (size_t ci = 0; ci < ev.child_count; ++ci) {
            auto &child = nodes[ev.child_first + ci];
            bool is_cont = IsContainer(child.type) || child.type == DATA_BANK2;
            if (is_cont) has_container = true;
            else         has_leaf = true;
        }
        if (has_container && !has_leaf)      n_built++;
        else if (!has_container && has_leaf) n_flat++;
        else if (has_container && has_leaf)  n_mixed++;

        // walk depth-1 children
        for (size_t ci = 0; ci < ev.child_count; ++ci) {
            auto &d1 = nodes[ev.child_first + ci];
            bool is_cont = IsContainer(d1.type) || d1.type == DATA_BANK2;

            auto &d1s = depth1_tags[d1.tag];
            d1s.count++;
            d1s.min_words = std::min(d1s.min_words, d1.data_words);
            d1s.max_words = std::max(d1s.max_words, d1.data_words);
            d1s.depths.insert(is_cont ? 1 : -1); // -1 = flat leaf at depth 1
            d1s.parent_tags.insert(ev.tag);

            // if container, walk depth-2 children
            if (is_cont) {
                for (size_t di = 0; di < d1.child_count; ++di) {
                    auto &d2 = nodes[d1.child_first + di];
                    uint64_t key = (uint64_t(d1.tag) << 32) | d2.tag;
                    auto &d2s = depth2_tags[key];
                    d2s.count++;
                    d2s.min_words = std::min(d2s.min_words, d2.data_words);
                    d2s.max_words = std::max(d2s.max_words, d2.data_words);
                    d2s.depths.insert(2);
                    d2s.parent_tags.insert(d1.tag);
                }
            }
        }
    }

    // === Output ===
    std::cout << "=== Bank Structure Debug (" << nphysics << " physics events, "
              << nrecords << " records) ===\n\n";

    // event structure classification
    std::cout << "--- Event Structure ---\n"
              << "  Built (ROC containers at depth 1): " << n_built << "\n"
              << "  Flat  (data banks at depth 1):     " << n_flat << "\n"
              << "  Mixed (both):                      " << n_mixed << "\n\n";

    // depth 0: event tags
    std::cout << "--- Depth 0: Event Tags ---\n"
              << std::setw(10) << "tag" << std::setw(8) << "count"
              << std::setw(16) << "size (words)" << "  label\n"
              << std::string(52, '-') << "\n";
    for (auto &[tag, s] : depth0_tags) {
        std::cout << "  " << hex(tag) << std::setw(8) << s.count << "  ";
        if (s.min_words == s.max_words)
            std::cout << std::setw(14) << s.min_words;
        else
            std::cout << std::setw(6) << s.min_words << " - " << std::setw(6) << s.max_words;
        std::cout << "  " << tag_label(tag) << "\n";
    }

    // depth 1: ROC crates and flat data banks
    std::cout << "\n--- Depth 1: ROC Crates / Flat Data Banks ---\n"
              << std::setw(10) << "tag" << std::setw(8) << "count"
              << std::setw(16) << "size (words)" << std::setw(10) << "layout"
              << "  identity\n"
              << std::string(75, '-') << "\n";
    for (auto &[tag, s] : depth1_tags) {
        bool as_container = s.depths.count(1) > 0;
        bool as_flat = s.depths.count(-1) > 0;

        std::cout << "  " << hex(tag) << std::setw(8) << s.count << "  ";
        if (s.min_words == s.max_words)
            std::cout << std::setw(14) << s.min_words;
        else
            std::cout << std::setw(6) << s.min_words << " - " << std::setw(6) << s.max_words;

        // layout column
        if (as_container && as_flat)      std::cout << std::setw(10) << "BOTH";
        else if (as_container)            std::cout << std::setw(10) << "container";
        else                              std::cout << std::setw(10) << "flat";

        // identity
        std::cout << "  ";
        if (tag == cfg.trigger_bank_tag)  std::cout << "TRIGGER BANK";
        else if (tag == cfg.ti_master_tag) std::cout << "TI MASTER";
        else {
            // check roc_tags
            bool found = false;
            for (auto &re : cfg.roc_tags) {
                if (re.tag == tag) {
                    std::cout << re.name;
                    if (!re.type.empty()) std::cout << " [" << re.type << "]";
                    found = true;
                    break;
                }
            }
            if (!found) {
                // check known data bank tags (for flat events)
                auto *bi = findBankInfo(tag);
                if (bi) std::cout << bi->name << " (" << bi->status << ")";
                else    std::cout << "*** UNKNOWN ***";
            }
        }
        std::cout << "\n";
    }

    // depth 2: data banks inside ROC containers
    std::cout << "\n--- Depth 2: Data Banks Inside ROC Crates ---\n"
              << std::setw(10) << "parent" << std::setw(10) << "tag"
              << std::setw(8) << "count" << std::setw(16) << "size (words)"
              << "  identity / dispatch\n"
              << std::string(80, '-') << "\n";

    for (auto &[key, s] : depth2_tags) {
        uint32_t parent = static_cast<uint32_t>(key >> 32);
        uint32_t tag = static_cast<uint32_t>(key & 0xFFFFFFFF);

        std::cout << "  " << hex(parent) << "  " << hex(tag)
                  << std::setw(8) << s.count << "  ";
        if (s.min_words == s.max_words)
            std::cout << std::setw(14) << s.min_words;
        else
            std::cout << std::setw(6) << s.min_words << " - " << std::setw(6) << s.max_words;

        // dispatch status
        std::cout << "  ";
        if (tag == cfg.ti_bank_tag)           std::cout << "-> decodeTIBank()";
        else if (tag == cfg.run_info_tag)     std::cout << "-> decodeRunInfo()";
        else if (tag == cfg.fadc_composite_tag) std::cout << "-> Fadc250Decoder";
        else if (tag == cfg.fadc_raw_tag)     std::cout << "-> Fadc250RawDecoder";
        else if (cfg.is_ssp_bank(tag))        std::cout << "-> SspDecoder";
        else if (tag == 0xE122)               std::cout << "-> VtpDecoder";
        else if (tag == cfg.tdc_bank_tag)     std::cout << "-> TdcDecoder";
        else if (tag == cfg.adc1881m_bank_tag) std::cout << "-> Adc1881mDecoder";
        else if (tag == cfg.daq_config_tag)   std::cout << "-> skip (config string)";
        else {
            auto *bi = findBankInfo(tag);
            if (bi)  std::cout << bi->name << " (" << bi->status << ")";
            else     std::cout << "*** UNKNOWN — not dispatched ***";
        }
        // flag if parent is not a known ROC
        bool parent_known = (parent == cfg.ti_master_tag);
        for (auto &re : cfg.roc_tags)
            if (re.tag == parent) { parent_known = true; break; }
        if (!parent_known)
            std::cout << " [parent ROC unknown]";

        std::cout << "\n";
    }

    // summary of dispatch coverage
    std::cout << "\n--- Dispatch Coverage ---\n";
    int total_d2 = 0, dispatched_d2 = 0, skipped_d2 = 0, unknown_d2 = 0;
    for (auto &[key, s] : depth2_tags) {
        uint32_t tag = static_cast<uint32_t>(key & 0xFFFFFFFF);
        total_d2 += s.count;
        if (tag == cfg.ti_bank_tag || tag == cfg.run_info_tag ||
            tag == cfg.fadc_composite_tag || tag == cfg.fadc_raw_tag ||
            cfg.is_ssp_bank(tag) || tag == cfg.adc1881m_bank_tag ||
            tag == cfg.tdc_bank_tag || tag == 0xE122)
            dispatched_d2 += s.count;
        else if (tag == cfg.daq_config_tag ||
                 (findBankInfo(tag) && std::string(findBankInfo(tag)->status) == "info"))
            skipped_d2 += s.count;
        else {
            auto *bi = findBankInfo(tag);
            if (bi)  skipped_d2 += s.count;   // known but no decoder
            else     unknown_d2 += s.count;
        }
    }
    std::cout << "  Total depth-2 banks:  " << total_d2 << "\n"
              << "  Dispatched (decoded): " << dispatched_d2 << "\n"
              << "  Known (no decoder):   " << skipped_d2 << "\n"
              << "  Unknown:              " << unknown_d2 << "\n";

    if (verbose) {
        // per-event detail: print the first few events showing their structure
        std::cout << "\n--- First 5 Event Structures ---\n";
        EvChannel ch2;
        ch2.SetConfig(cfg);
        // reopen not practical; this would need a separate pass.
        // Instead, point user to -m tree mode.
        std::cout << "  (use -m tree -n 5 for per-event bank tree)\n";
    }

    return 0;
}

// --- main -------------------------------------------------------------------
static void usage(const char *prog)
{
    std::cerr
        << "EVIO file structure diagnostic tool\n\n"
        << "Usage:\n"
        << "  " << prog << " <file> [options]\n\n"
        << "Modes (default: summary):\n"
        << "  -m tree       Print bank tree\n"
        << "  -m tags       List all unique bank tags with stats\n"
        << "  -m epics      Dump all EPICS event text\n"
        << "  -m vtp        Decode 0xE122 VTP banks (blocks / EC peaks / EC clusters)\n"
        << "  -m rf         Decode RF time TDC hits (ROC 0x40, bank 0xE107)\n"
        << "  -m event      Detailed dump of a single record\n"
        << "  -m triggers   List trigger bit counts (add -v for per-event detail)\n"
        << "  -m trig-debug Cross-correlate event tags, TI event type, and FP trigger bits\n"
        << "  -m bank-debug Verify bank structure, depth layout, and dispatch coverage\n\n"
        << "Options:\n"
        << "  -D <file>     DAQ configuration (auto-searches daq_config.json if omitted)\n"
        << "  -n <N>        Number of events (tree mode, default 5) or event number (event mode)\n"
        << "  -v            Verbose output (triggers mode: print every event)\n";
}

int main(int argc, char *argv[])
{
    std::string daq_config_file;
    std::string mode;
    int num = 5;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "D:m:n:vh")) != -1) {
        switch (opt) {
        case 'D': daq_config_file = optarg; break;
        case 'm': mode = optarg; break;
        case 'n': num = std::atoi(optarg); break;
        case 'v': verbose = true; break;
        default:  usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string path = argv[optind];

    // auto-search for daq_config.json if not specified.  Prefer the
    // install-aware resolver (env var → exe-relative → compile default),
    // then fall back to CWD-relative paths for dev-in-tree runs.
    if (daq_config_file.empty()) {
        std::string db_dir = prad2::resolve_data_dir(
            "PRAD2_DATABASE_DIR",
            {"../share/prad2evviewer/database"},
            DATABASE_DIR);
        if (!db_dir.empty()) {
            std::string cand = db_dir + "/daq_config.json";
            std::ifstream f(cand);
            if (f.good()) daq_config_file = std::move(cand);
        }
    }
    if (daq_config_file.empty()) {
        for (auto p : {"daq_config.json", "database/daq_config.json", "../database/daq_config.json"}) {
            std::ifstream f(p);
            if (f.good()) { daq_config_file = p; break; }
        }
    }

    evc::DaqConfig daq_cfg;
    if (daq_config_file.empty() || !evc::load_daq_config(daq_config_file, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_config_file.empty() ? " (not found)" : ": " + daq_config_file)
                  << "\n";
        return 1;
    }
    std::cerr << "DAQ config: " << daq_config_file
              << " (adc_format=" << daq_cfg.adc_format << ")\n";

    EvChannel ch;
    ch.SetConfig(daq_cfg);
    if (ch.OpenAuto(path) != status::success) {
        std::cerr << "Failed to open: " << path << "\n";
        return 1;
    }

    int rc;
    if      (mode.empty())         rc = doSummary(ch);
    else if (mode == "tree")       rc = doTree(ch, num);
    else if (mode == "tags")       rc = doTags(ch);
    else if (mode == "epics")      rc = doEpics(ch, num);
    else if (mode == "vtp")        rc = doVtp(ch, num);
    else if (mode == "rf")         rc = doRf(ch, num);
    else if (mode == "triggers")   rc = doTriggers(ch, verbose);
    else if (mode == "trig-debug") rc = doTrigDebug(ch, verbose);
    else if (mode == "bank-debug") rc = doBankDebug(ch, verbose);
    else if (mode == "event")      rc = doEvent(ch, num);
    else if (mode == "summary")    rc = doSummary(ch);
    else {
        std::cerr << "Error: unknown mode '" << mode << "'.\n"
                  << "Valid modes: summary, tree, tags, epics, vtp, rf, event, "
                     "triggers, trig-debug, bank-debug.\n";
        ch.Close();
        return 1;
    }

    ch.Close();
    return rc;
}
