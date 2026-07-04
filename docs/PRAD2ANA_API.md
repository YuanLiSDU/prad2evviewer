# prad2ana — C++ API Reference

`libprad2ana.a` — the offline replay and physics-analysis library.
Depends on [`prad2dec`](PRAD2DEC_API.md), [`prad2det`](PRAD2DET_API.md),
and ROOT 6.0+. All public symbols live in namespace `analysis`, except
the gain-factor helpers under `prad2`.

For a high-level walkthrough of the executables built on top of this
library (`prad2ana_replay_rawdata`, `prad2ana_replay_recon`,
`prad2ana_epCalib`, …) see [`analysis/README.md`](../analysis/README.md).
This document is the symbol reference for callers that link directly
against `libprad2ana.a` (analysis tools, ACLiC scripts in
`analysis/scripts/`, downstream user code).

| Header | Public symbols |
|---|---|
| [`Replay.h`](#replayh) | `analysis::Replay`, type aliases `EventVars`, `EventVars_Recon` |
| [`ConfigSetup.h`](#configsetuph) | `analysis::RunConfig` (alias to `prad2::RunConfig`), `gRunConfig`, `LabTransforms`, `BuildLabTransforms`, `ApplyToLab`, `TransformDetData`, `get_run_str`, `get_run_int` |
| [`PhysicsTools.h`](#physicstoolsh) | `analysis::PhysicsTools`, `GEMHit`, `HCHit`, `DataPoint`, `MollerEvent`, `MollerData` |
| [`MatchingTools.h`](#matchingtoolsh) | `analysis::MatchingTools`, `MatchHit`, `MatchHit_perChamber`, `MatchFlag`, `ProjectHit`, `GetProjection*`, `GetProjectionHits` |
| [`gain_factor.h`](#gain_factorh) | `prad2::GainFactor`, `GainFactorTable`, `GainCorrTable`, `FindGainFactorFile`, `LoadGainFactors`, `ComputeGainCorrection` |

---

## `Replay.h`

`analysis::Replay` — converts raw DAQ data (EVIO) into ROOT trees.
Three replay entry points:

1. **Raw replay** (`Process`) — per-channel waveform/peak data,
   `events` tree (`prad2::RawEventData`).
2. **Recon replay** (`ProcessWithRecon`) — full reconstruction
   (HyCal clusters, GEM hits, HyCal↔GEM matches), `recon` tree
   (`prad2::ReconEventData`).
3. **X17 recon replay** (`ProcessWithReconX17`) — cluster-trigger
   reconstruction with the X17 blind-sample selection, also written to
   a `recon` tree.

All replay entry points also write the side trees `scalers` (DSC2) and `epics`
(0x001F text banks) — see [`prad2det`](PRAD2DET_API.md#eventdata_ioh).

### Type aliases

```cpp
using EventVars       = prad2::RawEventData;
using EventVars_Recon = prad2::ReconEventData;
```

### Construction and configuration

```cpp
Replay r;
r.LoadDaqConfig(json_path);          // delegates to evc::load_daq_config
r.LoadHyCalMap(hycal_map_json_path); // populates daq_map_ + module_types_
```

Calling `LoadHyCalMap` is strongly recommended — without it every
channel reports `MOD_UNKNOWN` and `module_id` encoding falls back to
HyCal-only conventions.

### Channel introspection

```cpp
std::string moduleName(int roc, int slot, int ch) const;
prad2::ModuleType moduleType(int roc, int slot, int ch) const;
int moduleID(int roc, int slot, int ch) const;
```

`moduleID` returns the globally-unique encoding documented in
[`EventData.h`](PRAD2DET_API.md#eventdatah) (G: 1..1156, W: 1001..2152,
Veto: 3001..3004, LMS: 3100..3103). Returns `-1` for unknown channels.

### Replay drivers

```cpp
bool Process(const std::string &input_evio,
             const std::string &output_root,
             RunConfig &gRunConfig,
             const std::string &db_dir,
             int max_events = -1,
             bool write_peaks = false,
             const std::string &daq_config_file = "");

bool ProcessWithRecon(const std::string &input_evio,
                      const std::string &output_root,
                      RunConfig &gRunConfig,
                      const std::string &db_dir,
                      const std::string &daq_config_file = "",
                      const std::string &gem_ped_file    = "",
                      float zerosup_override            = 0.f,
                      bool prad1                        = false);

bool ProcessWithReconX17(const std::string &input_evio,
                         const std::string &output_root,
                         RunConfig &gRunConfig,
                         const std::string &db_dir,
                         const std::string &daq_config_file = "",
                         const std::string &gem_ped_file    = "",
                         float zerosup_override            = 0.f);
```

`max_events <= 0` ⇒ process every event. `write_peaks=true` adds the
optional soft-analyzer + firmware-mode peak branches to the `events`
tree (`hycal.peak_*`, `hycal.daq_peak_*`).

`ProcessWithRecon` runs the full pipeline via
[`prad2::PipelineBuilder`](PRAD2DET_API.md#pipelinebuilderh):
HyCal clustering, GEM strip clustering, HyCal↔GEM matching with the
runinfo σ parameters; `gem_ped_file` overrides the per-run pedestal,
`zerosup_override > 0` overrides the GEM ZS threshold, `prad1=true`
selects the legacy ADC1881M readout path.

`ProcessWithReconX17` uses the same detector pipeline and output data
type, but accepts events carrying `TBIT_1cl`, `TBIT_2cl`, `TBIT_3cl`,
`TBIT_lms`, or `TBIT_alpha`. It keeps every 1-/2-cluster event and only
the deterministic 10% of 3-cluster events for which
`event_num % 10 == 8`.

---

## `ConfigSetup.h`

Analysis-side helpers around `prad2::RunConfig`. Header-only.

### Re-exports

```cpp
using RunConfig = ::prad2::RunConfig;
using ::prad2::LoadRunConfig;
using ::prad2::WriteRunConfig;

inline RunConfig gRunConfig;   // single-run global, multi-run code should use locals
```

See the [`prad2det`](PRAD2DET_API.md#runinfoconfigh) reference for
`RunConfig` fields and `LoadRunConfig` selection rules.

### Lab-frame transforms

```cpp
struct LabTransforms {
    DetectorTransform                hycal;
    std::array<DetectorTransform, 4> gem;
};

LabTransforms BuildLabTransforms(const RunConfig &geo = gRunConfig);

template <typename Hit>
void ApplyToLab(const DetectorTransform &xform, Hit &h);   // in-place lab = R*[h.x,h.y,h.z] + t
```

`BuildLabTransforms` populates each `DetectorTransform` via `set(...)`,
so the rotation matrices are precomputed before `toLab` is called per
hit. Build once per run and reuse.

### Translation-only calibration shift

```cpp
void TransformDetData(MollerData &mollers,
                      float detX, float detY, float ZfromTarget);
```

Used by `det_calib` to apply per-detector alignment offsets to fitted
Moller pairs. Not a detector-frame transform — plain arithmetic.

### Run-number filename parsers

```cpp
std::string get_run_str(const std::string &file_name);   // "unknown" on failure
int         get_run_int(const std::string &file_name);   // -1 on failure
```

Pulls digits from `prad_<digits>...` style filenames.

---

## `PhysicsTools.h`

`analysis::PhysicsTools` — owns the per-module energy histograms,
Moller geometry, gain-monitoring histograms, and the kinematic
calculations used by epCalib / det_calib / gain_replay / matching.

### Per-event types

`GEMHit { x, y, z, det_id }` — `det_id ∈ 0..3` for GEM1..GEM4 (5 = unset).

`HCHit { x, y, z, energy, center_id, flag }`.

`DataPoint { x, y, z, E }` — used inside `MollerEvent`.

```cpp
typedef std::pair<DataPoint, DataPoint> MollerEvent;
typedef std::vector<MollerEvent>        MollerData;
```

### Construction

```cpp
explicit PhysicsTools(fdec::HyCalSystem &hycal);
```

The reference must outlive the `PhysicsTools` instance. Histograms are
constructed in the body of the constructor (one `TH1F` per HyCal module
plus the 2-D / Moller / gain-monitor histograms).

### Per-module energy histograms

```cpp
void  FillModuleEnergy(int module_id, float energy);
TH1F *GetModuleEnergyHist(int module_id) const;

void  FillEnergyVsModule(int module_id, float energy);
TH2F *GetEnergyVsModuleHist() const;

void  FillEnergyVsTheta(float theta_deg, float energy);
TH2F *GetEnergyVsThetaHist() const;

void  FillNeventsModuleMap(int module_id);
void  FillNeventsModuleMap();    // populate from filled per-module hists
TH2F *GetNeventsModuleMapHist() const;
```

### Yield analysis (caller owns the returned histograms)

```cpp
std::unique_ptr<TH1F> GetEpYieldHist (TH2F *energy_theta, float Ebeam);
std::unique_ptr<TH1F> GetEeYieldHist (TH2F *energy_theta, float Ebeam);
std::unique_ptr<TH1F> GetYieldRatioHist(TH1F *ep_hist, TH1F *ee_hist);
```

### Moller-event histograms

```cpp
void  Fill2armMollerPosHist(float x, float y);
TH2F *Get2armMollerPosHist() const;

void  FillMollerPhiDiff(float phi_diff);
void  FillMollerXY    (float x, float y);
void  FillMollerZ     (float z);
TH1F *GetMollerPhiDiffHist() const;
TH1F *GetMollerXHist()       const;
TH1F *GetMollerYHist()       const;
TH1F *GetMollerZHist()       const;
```

### Resolution / calibration

```cpp
std::array<float, 3> FitPeakResolution(int module_id) const;   // {peak, sigma, chi2}
void                 Resolution2Database(int run_id);
TF1                  nonLinearity_func_;
```

### Gain analysis

```cpp
struct GainResult {
    std::string name;
    float lms_peak  = 0;
    float lms_sigma = 0;
    float lms_chi2  = 0;
    float g[4]      = {};   // g[1..3] = mod_lms * alpha_ref[j] / lms_ref[j]
};

void                    ComputeModuleGains();
std::vector<GainResult> module_gains_;   // indexed by module index
float                   GetModuleGainFactor(int module_id) const;   // mean of g[1..3]
```

`ComputeModuleGains` fits LMS / α reference channels and every W-module,
populates `module_gains_`, and resets the source histograms.

### Gain-monitoring fillers (replay-time histograms)

```cpp
void Fill_lmsCH_lmsHeight   (int lms_id, float height);
void Fill_lmsCH_lmsIntegral (int lms_id, float integral);
void Fill_lmsCH_alphaHeight (int lms_id, float height);
void Fill_lmsCH_alphaIntegral(int lms_id, float integral);
void Fill_modCH_lmsHeight   (int module_id, float height);
void Fill_modCH_lmsIntegral (int module_id, float integral);

TH1F *Get_lmsCH_lmsHeightHist  (int lms_id) const;
TH1F *Get_lmsCH_lmsIntegralHist(int lms_id) const;
TH1F *Get_lmsCH_alphaHeightHist(int lms_id) const;
TH1F *Get_lmsCH_alphaIntegralHist(int lms_id) const;
TH1F *Get_modCH_lmsHeightHist  (int module_id) const;
TH1F *Get_modCH_lmsIntegralHist(int module_id) const;
```

### Kinematics (static)

```cpp
static float ExpectedEnergy(float theta_deg, float Ebeam,
                            const std::string &type);   // "ep" or "ee"
static float EnergyLoss   (float theta_deg, float E);   // target + windows
```

### Moller geometry

```cpp
std::array<float, 2> GetMollerCenter(MollerEvent &e1, MollerEvent &e2);
float                GetMollerZdistance(MollerEvent &e, float Ebeam);
float                GetMollerPhiDiff (MollerEvent &e1);   // ≈ 180° for elastic ee
float                GetPhiAngle      (float x, float y);
```

---

## `MatchingTools.h`

`analysis::MatchingTools` — HyCal cluster ↔ GEM hit matching with the
projection-plus-cut algorithm ported from
`PRadAnalyzer/PRadDetMatch.cpp`.

### Matching flag enum

```cpp
enum MatchFlag : uint32_t {
    kGEM1Match = 0,  // bit 0
    kGEM2Match = 1,
    kGEM3Match = 2,
    kGEM4Match = 3,
};
```

### Projection helpers

```cpp
struct ProjectHit { float x_proj, y_proj, z_proj; };

ProjectHit GetProjectionHits(float x, float y, float z, float projection_z);
void GetProjection(HCHit &hc,                       float projection_z);
void GetProjection(std::vector<HCHit> &hc,          float projection_z);
void GetProjection(GEMHit &gem,                     float projection_z);
void GetProjection(std::vector<GEMHit> &gem,        float projection_z);
```

`GetProjection` updates `(x, y, z)` of each hit to project a straight
line from the target through the hit to the requested `projection_z`.

### Match outputs

`MatchHit` — one HyCal cluster paired with all candidate GEM hits per
detector, plus the chosen "best" pair:

```cpp
class MatchHit {
public:
    HCHit                hycal_hit;
    std::vector<GEMHit>  gem1_hits, gem2_hits, gem3_hits, gem4_hits;

    GEMHit               gem[2];      // best-matched upstream and downstream
    uint32_t             mflag = 0;   // OR of MatchFlag bits
    uint16_t             hycal_idx = 0;

    MatchHit(const HCHit &, std::vector<GEMHit> &g1, std::vector<GEMHit> &g2,
             const std::vector<GEMHit> &g3, const std::vector<GEMHit> &g4);
};
```

`MatchHit_perChamber` — the per-chamber variant; stores the best match
per detector as a flat `[det_id][x/y/z]` array for analyses that don't
collapse to one upstream/downstream pair:

```cpp
class MatchHit_perChamber {
public:
    HCHit       hycal_hit;
    float       gem_hits[4][3] = {};
    uint32_t    mflag = 0;
    uint16_t    hycal_idx = 0;

    explicit MatchHit_perChamber(const HCHit &);
};
```

### `MatchingTools` methods

```cpp
std::vector<MatchHit> Match(
    std::vector<HCHit> &hycalHits,
    const std::vector<GEMHit> &gem1, const std::vector<GEMHit> &gem2,
    const std::vector<GEMHit> &gem3, const std::vector<GEMHit> &gem4) const;

std::vector<MatchHit_perChamber> MatchPerChamber(
    std::vector<HCHit> &hycalHits,
    const std::vector<GEMHit> &gem1, const std::vector<GEMHit> &gem2,
    const std::vector<GEMHit> &gem3, const std::vector<GEMHit> &gem4) const;

void SetMatchRange     (float range);   // mm; default 15
void SetSquareSelection(bool sq);       // true = square cut, false = circular
```

---

## `gain_factor.h`

Header-only loaders for the per-module LMS gain-factor database. Lives
under namespace `prad2` (not `analysis`) so other libraries can share it
without dragging ROOT.

### File format

`<dir>/prad_XXXXXX_LMS.dat` — whitespace-delimited:

```
Name  lms_peak  lms_sigma  lms_chi2/ndf  g1  g2  g3
```

Only `W*` and `G*` lines are read; LMS header rows and other prefixes
are silently skipped.

### Selection rule

Same as `LoadRunConfig`: `run_num >= 0` ⇒ largest run ≤ requested;
`run_num < 0` ⇒ latest. If no file satisfies `run ≤ run_num`,
`FindGainFactorFile` falls back to the nearest available file with a
warning.

### Types

```cpp
struct GainFactor { float g[3] = {0, 0, 0}; };   // g1, g2, g3

struct GainFactorTable {
    static constexpr int MAX_W = 1157;   // W1..W1156
    static constexpr int MAX_G = 901;    // G1..G900
    GainFactor w[MAX_W];
    GainFactor g[MAX_G];
    int  run_number = -1;
    bool loaded     = false;
};

struct GainCorrTable {
    struct Entry {
        float corr[3] = {1.f, 1.f, 1.f};   // ref.g[j] / cur.g[j]
        float avg     = 1.f;               // mean over non-zero corr[]
    };
    Entry w[GainFactorTable::MAX_W];
    Entry g[GainFactorTable::MAX_G];
    int   ref_run = -1;
    int   cur_run = -1;
};
```

### Functions

```cpp
std::string     prad2::FindGainFactorFile(const std::string &dir, int run_num);
GainFactorTable prad2::LoadGainFactors   (const std::string &dir, int run_num);

GainCorrTable prad2::ComputeGainCorrection(const GainFactorTable &ref_tbl,
                                           const GainFactorTable &cur_tbl);
GainCorrTable prad2::ComputeGainCorrection(const std::string &dir,
                                           int cur_run, int ref_run);
```

`new_adc2mev = old_adc2mev * corr.w[id].avg` is the typical applier
(see top of header for the worked example).

---

## Build / link

```cmake
target_link_libraries(your_tool PRIVATE prad2ana)
# transitively pulls in: prad2dec, prad2det, nlohmann_json,
#                        ROOT::{Core, Tree, RIO, Hist, Graf, Gpad, Spectrum}
```

`libprad2ana.a` is built with PIC so it can also be linked into ACLiC-
built shared objects (the `analysis/scripts/*.C` macros in installed
mode).

## Dependencies

- [`prad2dec`](PRAD2DEC_API.md) — EVIO reader, decoders, soft + firmware
  waveform analyzers.
- [`prad2det`](PRAD2DET_API.md) — `HyCalSystem`, `GemSystem`,
  `PipelineBuilder`, `RunConfig`, replay tree schema.
- ROOT 6.0+ — `TFile`, `TTree`, `TH1F`, `TH2F`, `TF1`, `TSpectrum`.
- [nlohmann/json](https://github.com/nlohmann/json) — fetched
  automatically by the top-level CMake.
