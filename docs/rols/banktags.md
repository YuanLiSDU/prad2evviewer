# PRad-II EVIO Bank Tag Reference

Reference for EVIO bank tags in PRad-II data.

**Authoritative source:** [`clonbanks_20260406.xml`](clonbanks_20260406.xml) — the
official CLAS/Hall-B EVIO bank dictionary. Always cross-check against this file
when adding new bank support; this Markdown is a curated PRad-II view of it,
augmented with what we have actually observed in data and ROL sources
(rol1.c, rol2.c, vtp1mpd.c) plus the trigger configuration (prad_v0.trg).

The "Banks present in run 023453" table reflects an early-stage tag dump and
may not show banks introduced in later runs (e.g. `0xE122` VTP).

## Trigger System Overview

```
FADC250 → VTP (per-crate clustering/sums) → SSP (cross-crate trigger logic)
                                                  ↓
                                            TRGBIT 0-7 via P2
                                                  ↓
                              TI/TS Front Panel inputs (bits 8-15)
                              v1495 → FP inputs (bits 24-27)
                                                  ↓
                              TS: mask + prescale + trigger table → event_type
                                                  ↓
                              Event Builder: event_tag = 0x80 + event_type
```

### FP Trigger Input Bit Map (from prad_v0.trg)

The 32-bit FP trigger inputs (TI master d[5]) have these assignments:

**SSP PRAD trigger bits (bits 8-15), routed via SSP P2 outputs:**

| Bit | Hex        | SSP P2 | Source     | Condition                              | Prescale |
|-----|------------|--------|------------|----------------------------------------|----------|
| 8   | 0x00000100 | P2OUT0 | TRGBIT 0   | Raw energy sum > 1000                  | ÷2       |
| 9   | 0x00000200 | P2OUT1 | TRGBIT 1   | >= 1 cluster, cluster E > 1000 MeV     | ÷1       |
| 10  | 0x00000400 | P2OUT2 | TRGBIT 2   | >= 2 clusters, cluster E > 1000 MeV    | ÷1       |
| 11  | 0x00000800 | P2OUT3 | TRGBIT 3   | >= 3 clusters, cluster E > 1000 MeV    | ÷1       |
| 12  | 0x00001000 | P2OUT4 | TRGBIT 4   | (disabled, mult >= 255)                | —        |
| 13  | 0x00002000 | P2OUT5 | TRGBIT 5   | (disabled, mult >= 255)                | —        |
| 14  | 0x00004000 | P2OUT6 | TRGBIT 6   | (disabled, mult >= 255)                | —        |
| 15  | 0x00008000 | P2OUT7 | TRGBIT 7   | **100 Hz pulser** (monitoring)         | ÷1       |

**v1495 external trigger bits (bits 16-27):**

| Bit | Hex        | Source   | TS prescale           |
|-----|------------|----------|-----------------------|
| 23  | 0x00800000 | v1495    | OR from SD/FADC       |
| 24  | 0x01000000 | v1495    | **LMS** — ÷1 (no prescale) |
| 25  | 0x02000000 | v1495    | alpha — ÷16385        |
| 26  | 0x04000000 | v1495    | Faraday — ÷16385      |
| 27  | 0x08000000 | v1495    | Master OR — ÷16385    |

**TS FP input mask:** `0x0F00FF00` (monitors bits 8-15 and 24-27).

### Event Tag Convention

PRad-II uses CODA single-event mode:

```
event_tag = 0x80 + TI_event_type
```

`TI_event_type` is set by `tiLoadTriggerTable(3)` based on the TS trigger
decision from the FP inputs above.

**Important:** The event_type (8-bit, in TI event header d[0] bits 31:24) and
the FP trigger inputs (32-bit, in TI master d[5]) are **different signals**:
- `event_type` = the TS trigger DECISION (which trigger type was formed)
- `d[5]` = raw FP input SNAPSHOT at trigger time (all 32 bits, including signals
  that didn't contribute to THIS trigger due to prescaling or masking)

They are related but not identical — an FP bit can be active in d[5] without
contributing to the trigger (e.g., prescaled away), and the event_type encoding
depends on the trigger table mapping.

### CODA3 Built-Trigger Tags (not used in PRad-II)

| Tag    | Description                    |
|--------|--------------------------------|
| 0xFF50 | PEB (Primary Event Builder)    |
| 0xFF58 | PEB with sync flag             |
| 0xFF70 | SEB (Secondary Event Builder)  |
| 0xFF78 | SEB with sync flag             |
| 0xFF60 | Streaming ROC Raw              |
| 0xFF62 | Streaming Data Concentrator    |
| 0xFF64 | Streaming SEB                  |
| 0xFF66 | Streaming PEB                  |

---

## Depth 0 — Top-Level Events

### Physics Events (tag = 0x80 + TI_event_type)

| Tag    | TI type | FP source | Count | w/FADC | Size (words)  | Description                      |
|--------|---------|-----------|-------|--------|---------------|----------------------------------|
| 0x00A9 | 0x29    | bit 8: SSP TRGBIT0 | 796 | 795 | 240 - 71046 | **SSP raw sum > 1000** (÷2 prescale) |
| 0x00B0 | 0x30    | bit 15: SSP TRGBIT7 | 16144 | 180 | 240 - 2105 | **100 Hz pulser** (monitoring) |
| 0x00B9 | 0x39    | bit 24: v1495 LMS | 773 | 773 | 24660 - 54586 | **LMS trigger** (no prescale) |
| 0x00BA | 0x3A    | bit 25: v1495 alpha | 2 | 0 | 240 | **alpha** (÷16385 prescale) |
| 0x00BC | 0x3C    | bit 27: v1495 Master OR | 443 | 5 | 240 - 305 | **Master OR** (÷16385 prescale) |
| 0x00FA | 0xFA    | bits 9+10: TRGBIT1+2 | 1 | 1 | 1063 | **cluster trigger** (special, tag≠0x80+type) |

### Control Events (CODA2 legacy)

| Tag    | Count | Description               |
|--------|-------|---------------------------|
| 0x0011 | 2     | Prestart                  |
| 0x0012 | 2     | Go                        |
| 0x0014 | 2     | End                       |

CODA3 equivalents: 0xFFD1 (Prestart), 0xFFD2 (Go), 0xFFD4 (End).

### Other Event Types

| Tag    | Description                              |
|--------|------------------------------------------|
| 0x001F | EPICS slow control                       |
| 0x00C1 | Sync event                               |

---

## Depth 1 — ROC Banks (inside physics events)

Each physics event contains 16 children: 1 trigger bank + 1 TI master + 7 FADC crates + 7 TI slaves.

### Trigger Bank

| Tag    | Type   | Words | Description                                          |
|--------|--------|-------|------------------------------------------------------|
| 0xC000 | UINT32 | 3     | CODA trigger bank: `[event_number, event_tag, reserved]` |

### TI Master Crate

| Tag    | Type | ROC ID | Words     | Description                                |
|--------|------|--------|-----------|--------------------------------------------|
| 0x0027 | BANK | 39     | 121 - 636 | TI master — trigger info, SSP, run info    |

### HyCal FADC Crates (even tags — contain FADC waveform data)

| Tag    | Type | ROC ID | Name      | Crate | Words        |
|--------|------|--------|-----------|-------|--------------|
| 0x0080 | BANK | 128    | adchycal1 | 0     | 6 - 11768    |
| 0x0082 | BANK | 130    | adchycal2 | 1     | 6 - 12947    |
| 0x0084 | BANK | 132    | adchycal3 | 2     | 6 - 12896    |
| 0x0086 | BANK | 134    | adchycal4 | 3     | 6 - 12230    |
| 0x0088 | BANK | 136    | adchycal5 | 4     | 6 - 12845    |
| 0x008A | BANK | 138    | adchycal6 | 5     | 6 - 12896    |
| 0x008C | BANK | 140    | adchycal7 | 6     | 6 - 11871    |

Min size (6 words) = monitoring events with TI data only.
Max size (~12K words) = physics events with full FADC waveforms.

### TI Slave Crates (odd tags — TI timing data only)

| Tag    | Type | ROC ID | Name      | Crate | Words   |
|--------|------|--------|-----------|-------|---------|
| 0x0081 | BANK | 129    | tihycal1  | 0     | 6 - 89  |
| 0x0083 | BANK | 131    | tihycal2  | 1     | 6 - 90  |
| 0x0085 | BANK | 133    | tihycal3  | 2     | 6 - 90  |
| 0x0087 | BANK | 135    | tihycal4  | 3     | 6 - 90  |
| 0x0089 | BANK | 137    | tihycal5  | 4     | 6 - 90  |
| 0x008B | BANK | 139    | tihycal6  | 5     | 6 - 90  |
| 0x008D | BANK | 141    | tihycal7  | 6     | 6 - 90  |

### Tagger V1190 TDC Crate

| Tag    | Type | ROC ID | Name    | Crate | Words     | Description                                        |
|--------|------|--------|---------|-------|-----------|----------------------------------------------------|
| 0x008E | BANK | 142    | tagger  | 7     | 8 - 300   | One VME crate with three V1190 TDC boards. Carries `0xE10A` TI data and `0xE107` TDC hits. |

### GEM VTP/MPD Crates (not in run 023453)

| Tag    | Type | ROC ID | Name    | Crate | ROL source  |
|--------|------|--------|---------|-------|-------------|
| 0x0031 | BANK | 49     | gemroc1 | 49    | vtp1mpd.c   |
| 0x0034 | BANK | 52     | gemroc2 | 52    | vtp1mpd.c   |

---

## Depth 2 — Data Banks (inside ROC crates)

> **Note:** Tables below were originally compiled from the early-stage run 023453.
> Newer runs may include additional banks (e.g. `0xE122`); see "Banks observed in
> later runs" below.

### Banks present in run 023453

| Tag    | Type      | Count  | Words       | Parents               | ROL source                          | Description                        |
|--------|-----------|--------|-------------|-----------------------|-------------------------------------|------------------------------------|
| 0xE10A | UINT32    | 272385 | 4 - 7       | all 15 ROC crates     | `BANKOPEN(0xe10a,1,rol->pid)` → `tiReadBlock()`, stripped by rol2.c | TI hardware data         |
| 0xE101 | COMPOSITE | 11163  | 63 - 12939  | even FADC crates only | rol2.c reformats 0xE109 → 0xE101    | FADC250 composite waveforms `c,i,l,N(c,Ns)` |
| 0xE10C | UINT32    | 18159  | 102         | 0x0027 only           | `BANKOPEN(0xe10C,1,rol->pid)`       | SSP trigger processor data         |
| 0xE10F | UINT32    | 18159  | 6           | 0x0027 only           | CODA/EB                             | Run info                           |
| 0xE10E | STRING    | 16     | 81 - 6565   | 0x0011 + all ROCs     | `BANKOPEN(0xe10E,3,rol->pid)`       | DAQ config readback (first events) |

### Banks observed in later runs (post-023453)

| Tag    | Type   | Words | Parents                          | Description                                                |
|--------|--------|-------|----------------------------------|------------------------------------------------------------|
| 0xE122 | UINT32 | 3-7   | TI slave crates (0x81/83/85/87/89/8B/8D) + 0x96 | "VTP Hardware Data" (per `clonbanks_20260406.xml`). Each HyCal crate's VTP reports through its TI-slave ROC: always EVTHDR (`0x12`) + TRGTIME (`0x13`, 2 words); crates whose VTP found a trigger-level cluster additionally ship **PRAD_CLUSTER** (TAG_EXP `0x1C` with 9-bit tag `0x1CC`, 2 words — seed module / energy / nhits / time) and a TRIGGER (`0x1D`, 2 words) summary. The CLAS12-style EC_PEAK / EC_CLUSTER records are never emitted by the PRad-II firmware. Full format below. |

#### `0xE122` VTP Hardware Data — full format (from official dictionary)

Self-describing 32-bit words; type code in bits[31:27] (block-level types) or
bits[31:23] (tag-expansion subtypes). PRad-II only cares about the ECAL records.

| Type | Name | Layout |
|------|------|--------|
| `0x10` | BLKHDR | `[26:22]` slot, `[17:08]` block#, `[7:0]` block_level |
| `0x11` | BLKTLR | `[26:22]` slot, `[21:0]` nwords |
| `0x12` | EVTHDR | `[26:0]` event_number |
| `0x13` | TRGTIME | word0 `[23:0]` time_lo, word1 `[23:0]` time_hi |
| `0x14` | **EC_PEAK** (2 words) | w0: `[26]` inst, `[25:24]` view, `[23:16]` time; w1: `[25:16]` coord, `[15:0]` energy |
| `0x15` | **EC_CLUSTER** (2 words) | w0: `[26]` inst, `[23:16]` time, `[15:0]` energy; w1: `[29:20]` coordW, `[19:10]` coordV, `[9:0]` coordU |
| `0x16` | HTCC_CLUSTER | (CLAS12, not used in PRad) |
| `0x17` | FT_CLUSTER | (CLAS12, not used) |
| `0x18` | FTOF_CLUSTER | (CLAS12, not used) |
| `0x19` | CTOF_CLUSTER | (CLAS12, not used) |
| `0x1A` | CND_CLUSTER | (CLAS12, not used) |
| `0x1B` | PCU_CLUSTER | (CLAS12, not used) |
| `0x1C` | Tag-expansion — subtype in `[26:23]`, full 9-bit tag in `[31:23]` | **PRAD_CLUSTER** = tag `0x1CC` ("0x10+0x0CC" in clonbanks.xml), 2 words: w0 `[13:0]` E; w1 `[26:15]` ID (1-900 → G1-900; bit 11 set → W(`ID & 0x7FF`), shown as 10001-11156 → W1-1156), `[14:11]` N (hits), `[10:0]` T (time). Other subtypes (DC_ROAD, HPS_*) are CLAS12/HPS, not used. |
| `0x1D` | TRIGGER | `[26:16]` trig_time, `[15:0]` trig_pattern_lo, w1 `[15:0]` trig_pattern_hi |
| `0x1E` | DNV (data not valid) | skip |
| `0x1F` | FILLER | skip |

Note: the `EC_PEAK` and `EC_CLUSTER` formats here come from the XML dictionary
and supersede the older single-word layouts in `rol2.c` (which appears to be an
earlier version of the parser).

**Decoder:** `prad2dec/src/VtpDecoder.cpp` parses block headers/trailers,
EVTHDR, TRGTIME, and stores EC_PEAK / EC_CLUSTER / PRAD_CLUSTER records in
`vtp::VtpEventData`. Other CLAS12 records (HTCC, FT, FTOF, CTOF, CND, PCU,
non-0x1CC tag-expansions, TRIGGER) are parsed past but not stored. The
decoder is invoked by `EvChannel::DecodeEvent` whenever a caller passes a
`vtp::VtpEventData*` (3rd optional argument).

PRAD_CLUSTER was validated against offline HyCal clustering (run 24340,
19k events): 88% of records share the reconstructed cluster's center
module and the energy correlation is 0.99 with VTP E ≈ 0.75-0.8 × the
calibrated cluster energy (trigger-level gains); T peaks sharply at 43
(fixed trigger latency).  `evio_dump -m vtp <file>` prints the decoded
fields inline.

### Banks defined in the official dictionary

Names and field layouts taken from `clonbanks_20260406.xml`. Banks marked
"observed in PRad-II" are the ones we have actually decoded from PRad-II data
or referenced in the ROLs; the rest are listed for completeness so unexpected
tags can be looked up quickly.

| Tag    | Type      | Official name                                 | Notes / PRad-II usage                                              |
|--------|-----------|-----------------------------------------------|--------------------------------------------------------------------|
| 0xE101 | COMPOSITE | FADC250 Window Raw Data (mode 1)              | **observed** — primary HyCal waveform bank, format `c,i,l,N(c,Ns)` |
| 0xE102 | COMPOSITE | FADC250 Pulse Integral Data with Timing (mode 7) | (mode 7, not used)                                              |
| 0xE103 | COMPOSITE | FADC250 Pulse Integral Data (mode 3)          | (mode 3, not used)                                                 |
| 0xE104 | UINT32    | VSCM Hardware Data / VSCM Raw Data            | reserved in rol1.c, not used in PRad-II                            |
| 0xE105 | UINT32    | DCRB Hardware Data                            | reserved in rol1.c, not used in PRad-II (was misnamed "DCRB/DC/Vetroc") |
| 0xE107 | UINT32    | V1190 TDC Data                                | **observed** — tagger crate `0x008E` emits 0-22 words/event. Each word is one hit: `slot[31:27] \| edge[26] \| ch[25:19] \| tdc[18:0]`. Decoded by `prad2dec/src/TdcDecoder.cpp`. |
| 0xE108 | UINT32    | DCRBGTP Hardware Data                         | not used                                                           |
| 0xE109 | UINT32    | FADC250 Hardware Data (raw)                   | **observed in some configs** — pre-rol2.c reformatting; rol2.c converts → 0xE101 |
| 0xE10A | UINT32    | TI/TS Hardware Data                           | **observed** — TI master/slave timing, FP trigger pattern in d[5]  |
| 0xE10B | UINT32    | V1190/V1290 Hardware Data                     | reserved in rol1.c, not used                                       |
| 0xE10C | UINT32    | **SSP Hardware Data**                         | **observed** — official format is HPS-style (HPS_CLUSTER, HPS_TRIGGER, TRIGGER pattern). For PRad-II this is the **trigger processor data** in the TI master crate (0x0027). ⚠️ NOT the MPD/APV strip format we use for GEMs — that's a separate `0x0DE9`-style bank. |
| 0xE10D | UINT32    | GTP Hardware Data                             | not used                                                           |
| 0xE10E | CHAR      | Run Config File                               | **observed** — DAQ config readback string (first events only)      |
| 0xE10F | UINT32    | **HEAD bank**                                 | **observed** — version, run#, event#, unix_time, event_type, roc_pattern, classification, trigger_bits. (We previously called this "Run info" in code; the data fields we extract — run number, event count, unix time — are correct, just under the wrong name.) |
| 0xE111 | COMPOSITE | SVT Composite Data                            | not used                                                           |
| 0xE112 | UINT32    | HEAD bank raw format                          | block-level form of 0xE10F, not used in PRad-II                    |
| 0xE113 | UINT32    | SLAC SVT before disent                        | not used                                                           |
| 0xE114 | CHAR      | EPICS data                                    | **observed** — inner data bank inside top-level `0x001F` EPICS events |
| 0xE115 | UINT32    | DSC2 Scalers raw format                       | **observed** — used by `livetime` tool for gated/ungated counts (see DSC2 docs) |
| 0xE116 | COMPOSITE | DCRB Composite Data                           | not used                                                           |
| 0xE117 | UINT32    | ADC v792 raw format                           | not used                                                           |
| 0xE118 | UINT32    | MVT raw format                                | not used                                                           |
| 0xE119 | UINT32    | FTT raw format                                | not used                                                           |
| 0xE11A | UINT32    | TDC v775 raw format                           | not used                                                           |
| 0xE11B | COMPOSITE | MVT composite data                            | not used                                                           |
| 0xE11E | UINT32    | DAQ performance timers                        | not used                                                           |
| 0xE11F | UINT32    | SRS Raw Data                                  | not used                                                           |
| 0xE120 | UINT32    | FASTBUS Raw Data                              | **observed** — used by PRad-1 ADC1881M legacy decoder              |
| 0xE121 | UINT32    | PGEM crate Raw Data                           | not used (separate GEM crate readout)                              |
| 0xE122 | UINT32    | VTP Hardware Data                             | **observed in TI slave crates (3-7 words): EVTHDR + TRGTIME + optional PRAD_CLUSTER / TRIGGER; decoder in `prad2dec/src/VtpDecoder.cpp`**. Full format documented above. |
| 0xE123 | UINT32    | SSP-RICH Hardware Data                        | reserved in rol1.c, not used in PRad-II                            |
| 0xE124 | COMPOSITE | SSP-RICH Composite Data                       | not used                                                           |
| 0xE125 | UINT32    | **SIS3801 Scalers raw format**                | reserved in rol1.c. (Previously misnamed "Per-slot data" — corrected per official dictionary.) |
| 0xE126 | COMPOSITE | FADC250 Window Raw Data (mode 1 packed/compressed) | not used                                                       |
| 0xE127 | COMPOSITE | SVT Composite Data with tag 9 (TDCEVT)         | not used                                                           |
| 0xE128 | COMPOSITE | MVT composite data packed                     | not used                                                           |
| 0xE129 | COMPOSITE | TPC composite data packed                     | not used                                                           |
| 0xE130 | COMPOSITE | DCRB Composite Data with width                | not used                                                           |
| 0xE131 | UINT32    | VFTDC Hardware Data                           | not used                                                           |
| 0xE132 | COMPOSITE | VFTDC Composite Data                          | not used                                                           |
| 0xE133 | UINT32    | Helicity Decoder Hardware Data                | not used                                                           |
| 0xE134 | UINT32    | SAMPA/FELIX raw format                        | not used                                                           |
| 0xE135 | UINT32    | VMM raw format                                | not used                                                           |
| 0xE136 | UINT32    | MAROC raw format                              | not used                                                           |
| 0xE137 | COMPOSITE | MAROC Composite Data                          | not used                                                           |
| 0xE138 | UINT32    | PETIROC raw format                            | not used                                                           |
| 0xE139 | COMPOSITE | PETIROC Composite Data                        | not used                                                           |
| 0xE140 | UINT32    | (reserved for MPD raw format)                 | XML notes "use Hall A tag 0x0DE9, but reserve 0xE140". (Previously misnamed "Special data bank" — corrected.) |
| 0xE141 | UINT32    | FAV3 Hardware Data                            | reserved, not used                                                 |
| 0x0DE9 | UINT32    | **MPD raw format (Hall A tag)**               | **observed** — confirmed in run 001137 (gemroc1/gemroc2, 28066 words/event). PRad-II GEM readout uses this standard Hall-A tag. |

---

## Depth 3 — Inside FADC Composite Banks (0xE101)

| Tag    | Type      | Count | Words       | Description                                  |
|--------|-----------|-------|-------------|----------------------------------------------|
| 0x000D | CHAR8     | 11163 | 4           | Composite format descriptor: `"c,i,l,N(c,Ns)"` |
| 0x0000 | UNKNOWN32 | 11163 | 56 - 12932  | Composite data payload (slot headers + channel samples) |

---

## Data Word Formats

### TI Data (0xE10A) — after rol2.c block header stripping

The raw TI block (from `tiReadBlock()`) includes block header/trailer.
rol2.c strips these, outputting per-event data starting with the event header.

**TI Slave (4 words, nwords=3):**

| Word | Example      | Content                                    | Decoder config         |
|------|--------------|-------------------------------------------|------------------------|
| d[0] | `0x30010003` | Event header: `event_type(8) \| 0x01(8) \| nwords(16)` | baseline trigger_bits: `(d[0]>>24) & 0xFF` |
| d[1] | `0x000046EC` | Event number (32-bit)                      | `trigger_word = 1`     |
| d[2] | `0x100901FE` | Timestamp low (32-bit)                     | `time_low_word = 2`    |
| d[3] | `0x0000000B` | `evnum_high[19:16] \| ts_high[15:0]`       | `time_high_word = 3`, mask = 0xFFFF |

**TI Master (7 words, nwords=6):**

| Word | Example      | Content                                    | Decoder config         |
|------|--------------|-------------------------------------------|------------------------|
| d[0]-d[3] | | Same as slave                              |                        |
| d[4] | `0x00000000` | Trigger type byte (**always zero** in current config) | — (unused)   |
| d[5] | `0x00008000` | **32-bit FP trigger inputs** (see bit map above) | `trigger_type_word = 5`, mask = 0xFFFFFFFF |
| d[6] | `0x80000000` | Additional TI flags                        | —                      |

The `event_type` (d[0] bits 31:24) and FP trigger inputs (d[5]) are **related
through the trigger table** (`tiLoadTriggerTable(3)`), not by direct bit encoding:
- `event_type` = trigger table OUTPUT → determines event tag via `tag = 0x80 + event_type`
- `d[5]` = raw FP input SNAPSHOT → which detector signals were active at trigger time

**Confirmed mapping** (from trig-debug on run 023453):

| FP bit | FP signal              | → event_type | → event_tag | Count |
|--------|------------------------|-------------|-------------|-------|
| 8      | SSP TRGBIT0 (RawSum)   | 0x29        | 0x00A9      | 796   |
| 15     | SSP TRGBIT7 (Pulser)   | 0x30        | 0x00B0      | 16144 |
| 24     | v1495 LMS              | 0x39        | 0x00B9      | 773   |
| 25     | v1495 alpha            | 0x3A        | 0x00BA      | 2     |
| 27     | v1495 Master OR        | 0x3C        | 0x00BC      | 443   |
| 9+10   | TRGBIT1+2 (clusters)   | 0xFA        | 0x00FA      | 1     |

### Trigger Bank (0xC000)

| Word | Content                | Decoder config              |
|------|------------------------|-----------------------------|
| d[0] | Event number           | `trig_event_number_word = 0` |
| d[1] | Event tag (= top-level tag, e.g. 0xB0) | `trig_event_type_word = 1` |
| d[2] | Reserved (0)           |                             |

### Run Info Bank (0xE10F) — in TI master only

| Word | Content                | Decoder config              |
|------|------------------------|-----------------------------|
| d[0] | Header                 |                             |
| d[1] | Run number             | `ri_run_number_word = 1`    |
| d[2] | Event count            | `ri_event_count_word = 2`   |
| d[3] | Unix timestamp         | `ri_unix_time_word = 3`     |
| d[4] - d[5] | Additional info |                             |

### FADC250 Composite Format (0xE101)

Format string: `"c,i,l,N(c,Ns)"` — packed byte stream (little-endian).

Per slot:
- `c` (1 byte): slot ID
- `i` (4 bytes): trigger / event number
- `l` (8 bytes): 48-bit timestamp
- `N` (4 bytes): channel count

Per channel (repeated N times):
- `c` (1 byte): channel ID
- `N` (4 bytes): sample count
- `s` (2 bytes × N): ADC samples (uint16, native endian)

### FADC250 Raw Hardware Format (0xE109)

Self-describing 32-bit words with type code in bits[31:27]:

| Type (bits 31:27) | Name             | Key fields                            |
|-------------------|------------------|---------------------------------------|
| 0x00              | Block Header     | slot(5), module_id(4), block#(10), nevents(8) |
| 0x01              | Block Trailer    | slot(5), nwords(22)                   |
| 0x02              | Event Header     | slot(5), trigger#(22)                 |
| 0x03              | Trigger Time     | time(24), continuation: time_high(24) |
| 0x04              | Window Raw Data  | channel(4), width(12), then 2 ADC samples/word |
| 0x1F              | Filler           | skip                                  |

Sample data words (continuation after type 0x04 header):
- Bits 28:16 = ADC sample i (13 bits)
- Bits 12:0 = ADC sample i+1 (13 bits)

### V1190 TDC Format (0xE107)

rol2.c reformats the raw `0xE10B` v1190/v1290 hardware stream (block/global
headers, TDC headers, EOB markers, global trailer) into a flat array of hits.
Each 32-bit word is one TDC hit:

| Bits   | Width | Field | Notes                                                |
|--------|-------|-------|------------------------------------------------------|
| 31:27  | 5     | SLOT  | V1190 board slot in the VME crate (0-31)             |
| 26     | 1     | EDGE  | 0 = leading, 1 = trailing                            |
| 25:19  | 7     | CH    | Channel within the board (0-127)                     |
| 18:00  | 19    | TDC   | Timing value. For v1190, rol2.c left-shifts the native 100 ps count by 2 to match the v1290 25 ps LSB before truncating to 19 bits. |

No framing words are present in the payload — `nwords == 0` simply means no
hits for this event. PRad-II tagger observed range: 0-22 words per event.

**Decoder:** `prad2dec/src/TdcDecoder.cpp` — appends hits to `tdc::TdcEventData`
when the caller passes a `TdcEventData*` as the 4th optional argument to
`EvChannel::DecodeEvent`.

### SSP/MPD Bitfield Format (0xE10C, 0x0DE9)

Self-describing 32-bit words. Bit 31 = type-defining flag, bits 30:27 = type tag.

| Type | Name              | Key fields                              |
|------|-------------------|-----------------------------------------|
| 0    | Block Header      | slot(5), module_id(4), block#(10), nevents(8) |
| 1    | Block Trailer     | words_in_block(22), slot(5)             |
| 2    | Event Header      | trigger_number(27)                      |
| 3    | Trigger Time      | low(24) + continuation high(24)         |
| 5    | MPD Frame / APV   | Defining: mpd_id(5), fiber(6), flags(5). Continuation: 3 words per strip with 6 time samples (13-bit signed each) |
| 4, 6-12 | MPD Timestamp  | timestamp_fine(8), coarse(16+24), event_count(20) |
| 0xD  | MPD Debug         | Online common mode: 6 × 13-bit CM values across 3 words |
| 0xE  | Data Not Valid     | skip                                    |
| 0xF  | Filler            | skip                                    |

---

## Event Structure Diagram

```
[Physics Event: tag = 0x80 + TI_event_type] (depth 0)
 ├── [0xC000 UINT32, 3w] — trigger bank
 ├── [0x0027 BANK] — TI master crate
 │    ├── [0xE10A UINT32, 7w] — TI data (trigger#, timestamp, FP trigger bits)
 │    ├── [0xE10C UINT32, 102w] — SSP trigger processor data
 │    ├── [0xE10E STRING] — DAQ config string (first events only)
 │    └── [0xE10F UINT32, 6w] — run info (run#, unix time)
 ├── [0x0080 BANK] — HyCal FADC crate 1
 │    ├── [0xE10A UINT32, 4w] — TI data
 │    └── [0xE101 COMPOSITE] — FADC250 waveforms (physics triggers only)
 ├── [0x0081 BANK] — TI slave for crate 1
 │    └── [0xE10A UINT32, 4w] — TI data only
 ├── ... (7 FADC crates 0x80-0x8C + 7 TI slaves 0x81-0x8D)
 └── [0x008D BANK] — TI slave for crate 7
      └── [0xE10A UINT32, 4w]

[Control Event: tag = 0x0011/0x0012/0x0014] (depth 0)
 └── [0xE10E STRING] — DAQ config / run parameters
```
