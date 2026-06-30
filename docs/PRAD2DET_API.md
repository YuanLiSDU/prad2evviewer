# prad2det — C++ API Reference

Static library that consumes decoded front-end data from
[`prad2dec`](PRAD2DEC_API.md) and produces physics-level outputs:
HyCal per-module hits, GEM strip clusters, 2-D GEM hits, and lab-frame
geometry transforms.

For a high-level usage walkthrough see [`prad2det/README.md`](../prad2det/README.md).
This document is the symbol reference: every public class, struct, free
function, and namespace constant exported by the library.

| Header | Namespace | Public symbols |
|---|---|---|
| [`HyCalSystem.h`](#hycalsystemh) | `fdec` | `HyCalSystem`, `Module`, `ModuleType`, `Sector`, `SectorInfo`, `NeighborInfo`, `DaqAddr`, `LayoutFlag` |
| [`HyCalCluster.h`](#hycalclusterh) | `fdec` | `HyCalCluster`, `ClusterConfig`, `ModuleHit`, `ModuleCluster`, `ClusterHit`, `IClusterProfile`, `SimpleProfile`, `shower_depth()` |
| [`GemSystem.h`](#gemsystemh) | `gem` | `GemSystem`, `DetectorConfig`, `PlaneConfig`, `ApvConfig`, `ApvPedestal`, `ClusterConfig`, `StripHit`, `StripCluster`, `GEMHit`, `MapStrip()`, `MapApvStrips()` |
| [`GemCluster.h`](#gemclusterh) | `gem` | `GemCluster` |
| [`GemPedestal.h`](#gempedestalh) | `gem` | `GemPedestal` |
| [`DetectorTransform.h`](#detectortransformh) | (global) | `DetectorTransform` |
| [`RunInfoConfig.h`](#runinfoconfigh) | `prad2` | `RunConfig`, `LoadRunConfig()`, `WriteRunConfig()` |
| [`PipelineBuilder.h`](#pipelinebuilderh) | `prad2` | `PipelineBuilder`, `Pipeline` |
| [`EventData.h`](#eventdatah) | `prad2` | `RawEventData`, `ReconEventData`, `RawScalerData`, `RawEpicsData`, `ModuleType`, capacity constants |
| [`EventData_io.h`](#eventdata_ioh) | `prad2` | `SetRawWriteBranches`, `SetRawReadBranches`, `SetReconWriteBranches`, `SetReconReadBranches`, `SetScalerWriteBranches`, `SetScalerReadBranches`, `SetEpicsWriteBranches`, `SetEpicsReadBranches`, `FillScalerRow`, `FillEpicsRow`, `RawReadStatus`, `ReconReadStatus` |

---

## `HyCalSystem.h`

`fdec::HyCalSystem` — geometry, calibration, and neighbour lookup for the
HyCal calorimeter. Initialized once from JSON; immutable per event.

### Constants

| Constant | Value | Meaning |
|---|---|---|
| `HYCAL_MAX_MODULES` | 1800 | hard cap on module array (1728 real + margin) |
| `MAX_NEIGHBORS`    | 12   | max cross-sector neighbours per module |
| `PWO_ID0`          | 1000 | PRimEx ID offset for crystal modules (W) |
| `SECTOR_GRID_MAX_CELLS` | 34×34 | sector-grid capacity |

### Enums

| Enum | Values |
|---|---|
| `ModuleType` | `PbGlass=0`, `PbWO4=1`, `LMS=2`, `Veto=3`, `Unknown=-1` |
| `Sector`     | `Center=0`, `Top=1`, `Right=2`, `Bottom=3`, `Left=4`, `Max=5` |
| `LayoutFlag` (bit positions) | `kPbGlass`, `kPbWO4`, `kTransition`, `kInnerBound`, `kOuterBound`, `kDeadModule`, `kDeadNeighbor`, `kSplit`, `kLeakCorr` |

Helpers: `set_bit`, `clear_bit`, `test_bit` operate on the packed `Module::flag`.

### Structs

`Module` — one calorimeter cell.

| Field | Description |
|---|---|
| `name`, `id`, `index` | identity (`id` is PRimEx 1..1156 G / 1001..2152 W; `index` is contiguous) |
| `type` | `ModuleType` |
| `x`, `y`, `size_x`, `size_y` | geometry (mm) |
| `flag`, `sector`, `row`, `column` | layout |
| `daq` | `DaqAddr{crate, slot, channel}` |
| `cal_factor`, `cal_base_energy`, `cal_non_linear` | per-module calibration |
| `neighbors[]`, `neighbor_count` | pre-computed cross-sector neighbours |

Member helpers: `energize(adc)` → MeV (incl. non-linear correction);
`is_pwo4()`, `is_glass()`, `is_hycal()`; `is_neighbor(other, include_corners=true)`.

`SectorInfo` — per-sector geometry. `set_boundary(x1,y1,x2,y2)` /
`get_boundary(...)`.

`NeighborInfo { index, dx, dy, dist }` — quantized neighbour entry.

`DaqAddr { crate, slot, channel }`.

`SectorGrid` — flat 2-D grid for O(1) same-sector neighbour lookup.

### `HyCalSystem` methods

| Method | Description |
|---|---|
| `bool Init(const std::string &map_path)` | Load merged HyCal map JSON (one record per module: `n,t,geo,daq`). Returns `false` on error. |
| `int LoadCalibration(const std::string &calib_path)` | Load `[{name,factor,base_energy,non_linear}, …]`. Returns matched count or `-1`. |
| `int module_count() const` | |
| `const Module &module(int index) const` / `Module &module(int)` | |
| `const Module *module_by_name(const std::string&) const` | |
| `const Module *module_by_id(int primex_id) const` | |
| `const Module *module_by_daq(int crate, int slot, int ch) const` | |
| `double GetCalibConstant(int primex_id) const` | also `GetCalibBaseEnergy`, `GetCalibNonLinearity` |
| `void SetCalibConstant(int primex_id, double)` | also `SetCalibBaseEnergy`, `SetCalibNonLinearity` |
| `void PrintCalibConstants(const std::string &output_file) const` | |
| `const SectorInfo &sector_info(int s) const` | |
| `int get_sector_id(double x, double y) const` | |
| `void qdist(x1,y1,s1, x2,y2,s2, &dx,&dy) const` / `qdist(m1, m2, &dx,&dy)` | quantized inter-module distance |
| `template<typename Fn> void for_each_neighbor(int idx, bool include_corners, Fn) const` | walks ±1 grid + cross-sector list |
| `void SetPositionResolutionParams(float A, float B, float C)` | |
| `float GetPositionResolutionA/B/C() const` | |
| `float PositionResolution(float energy_mev) const` | σ(E) = √((A/√Eg)² + (B/Eg)² + C²), in mm |
| `static int name_to_id(const std::string&)` | |
| `int id_to_index(int id) const` | |
| `static std::string id_to_name(int id)` | |
| `static ModuleType parse_type(const std::string&)` | |

---

## `HyCalCluster.h`

`fdec::HyCalCluster` — Island clustering on top of a `HyCalSystem`.
Algorithm details: `docs/technical_notes/hycal_clustering/`.

### `ClusterConfig`

| Field | Default | Meaning |
|---|---|---|
| `min_module_energy` | 0 MeV | hit threshold |
| `min_center_energy` | 10 MeV | seed threshold |
| `min_cluster_energy`| 50 MeV | total-energy cut |
| `min_cluster_size`  | 1 | minimum hit count |
| `corner_conn`       | false | include diagonal neighbours in grouping |
| `split_iter`        | 6 | shower-split fixed-point iterations |
| `least_split`       | 0.01 | minimum split fraction |
| `log_weight_thres`  | 3.6 | W = max(0, thres + ln(E_i / E_tot)) |
| `non_linear_corr`   | true | apply per-module energy non-linearity correction |

### Per-event types

`ModuleHit { index, energy, time }`.

`ModuleCluster { center, hits, energy, flag }` with `add_hit(h)`.

`ClusterHit { center_id, x, y, energy, time, nblocks, npos, flag }` —
reconstructed cluster output.

`shower_depth(int center_id, float energy_mev)` — shower-max depth (mm)
discriminating PWO4 vs PbGlass via `PWO_ID0`.

### Cluster profile

`IClusterProfile` — abstract interface; one method
`GetFraction(ModuleType, dist, energy)`.

`SimpleProfile` — analytical fall-off (default).

### `HyCalCluster` methods

| Method | Description |
|---|---|
| `explicit HyCalCluster(const HyCalSystem &sys)` | |
| `void SetConfig(const ClusterConfig&)` / `const ClusterConfig &GetConfig() const` | |
| `void SetProfile(IClusterProfile *prof)` | takes ownership |
| `void Clear()` | |
| `void AddHit(int module_index, float energy, float time)` | |
| `void FormClusters()` | |
| `void ReconstructHits(std::vector<ClusterHit> &out) const` | |
| `void ReconstructMatched(std::vector<RecoResult> &out) const` | pairs `{const ModuleCluster*, ClusterHit}` for clusters passing thresholds |
| `const std::vector<ModuleCluster> &GetClusters() const` | |

Non-copyable (owns the profile pointer).

`SPLIT_MAX_HITS = 100`, `SPLIT_MAX_MAXIMA = 10`, `SplitContainer` —
internal split buffers, exposed for embedding.

---

## `GemSystem.h`

`gem::GemSystem` — GEM detector hierarchy and the per-event chain
(pedestal subtraction → common-mode → ZS → strip mapping). Auto-detects
full-readout vs online-ZS data per APV.

### Per-event types

`StripHit { strip, charge, max_timebin, position, cross_talk, ts_adc }`.

`StripCluster { position, peak_charge, total_charge, max_timebin,
cross_talk, hits }` — see comment in header about value-init importance.

`GEMHit { x, y, z, det_id, x_charge, y_charge, x_peak, y_peak,
x_max_timebin, y_max_timebin, x_size, y_size }`.

### Configuration types

`ApvPedestal { offset = 0, noise = 5000 }` (large default ⇒ no hits until calibrated).

`ApvConfig` — DAQ address (`crate_id, mpd_id, adc_ch`), detector mapping
(`det_id, plane_type ∈ {0=X,1=Y}, orient, plane_index, det_pos`), strip
mapping (`pin_rotate, shared_pos, hybrid_board, match`), per-strip
pedestal table, and per-APV CM range.

`PlaneConfig { type, size, n_apvs, pitch=0.4 mm }`.

`DetectorConfig { name, id, type="PRADGEM", planes[2] }`.

`ClusterConfig` — per-detector knobs consumed by `GemCluster`:

| Field | Default | |
|---|---|---|
| `min_cluster_hits`  | 1   | |
| `max_cluster_hits`  | 20  | |
| `consecutive_thres` | 1   | max gap between strips |
| `split_thres`       | 14  | charge-valley depth for splitting |
| `cross_talk_width`  | 2 mm | |
| `charac_dists`      | {} | cross-talk characteristic distances |
| `match_mode`        | 1 | 0 = ADC-sorted 1:1, 1 = Cartesian-with-cuts |
| `match_adc_asymmetry` | 0.8 | max \|Qx-Qy\|/(Qx+Qy); <0 disables |
| `match_time_diff`     | 50 ns | <0 disables |
| `ts_period`           | 25 ns | |

### Strip-mapping helpers (pure)

```cpp
int MapStrip(int ch, int plane_index, int orient,
             int pin_rotate=0, int shared_pos=-1,
             bool hybrid_board=true, int apv_channels=128,
             int readout_center=32);

std::vector<int> MapApvStrips(int plane_index, int orient,
                              int pin_rotate=0, int shared_pos=-1,
                              bool hybrid_board=true, int apv_channels=128,
                              int readout_center=32);
```

Both are stateless — Python bindings reproduce the layout without
instantiating a `GemSystem`. `GemSystem::buildStripMap()` delegates to
`MapStrip` (single source of truth).

### `GemSystem` methods

| Method | Description |
|---|---|
| `void Init(const std::string &map_file)` | Load `gem_map.json` (detectors + APVs). |
| `void LoadPedestals(const std::string &ped_file, const std::map<int,int> &crate_remap = {})` | Upstream APV-block text format. `crate_remap` maps file-side hardware crate IDs to logical IDs (e.g. 146→1, 147→2). Slot field in headers is ignored. |
| `void LoadCommonModeRange(const std::string &cm_file, const std::map<int,int> &crate_remap = {})` | Per-APV CM bounds for the Danning algorithm. |
| `void SetReconConfigs(std::vector<ClusterConfig>)` | Per-detector cluster/match parameters; vector clamped/padded to `GetNDetectors()`. |
| `const std::vector<ClusterConfig>& GetReconConfigs() const` | |
| `void Clear()` | |
| `void ProcessEvent(const ssp::SspEventData&)` | pedestal/CM/ZS + APV→strip mapping |
| `void Reconstruct(GemCluster &)` | applies per-detector config, fills `det_hits_` and `all_hits_` |
| `int GetNDetectors() const` | |
| `const std::vector<DetectorConfig>& GetDetectors() const` | |
| `const std::vector<StripHit>& GetPlaneHits(int det, int plane) const` | |
| `const std::vector<StripCluster>& GetPlaneClusters(int det, int plane) const` | |
| `const std::vector<GEMHit>& GetHits(int det) const` | |
| `const std::vector<GEMHit>& GetAllHits() const` | |
| `int FindApvIndex(int crate, int mpd, int adc) const` | O(1) |
| `int GetNApvs() const`, `const ApvConfig& GetApvConfig(int idx) const` | |
| `float GetHoleXOffset() const` | beam-hole X offset from match-APV strip positions, `0` if absent |
| `std::pair<float,float> GetActiveExtent(int det_id, int plane) const` | tight bbox of mapped strips |
| `bool IsChannelHit(int apv_idx, int ch) const`, `bool HasApvZsHits(int apv_idx) const`, `float GetProcessedAdc(int apv_idx, int ch, int ts) const` | per-APV ZS results, valid after `ProcessEvent` |
| `Get/Set CommonModeThreshold`, `ZeroSupThreshold`, `Get CrossTalkThreshold`, `RejectFirstTimebin`, `RejectLastTimebin`, `MinPeakAdc`, `MinSumAdc` | |

Internal capacity: `APV_STRIP_SIZE = 128`, `SSP_TIME_SAMPLES = 6`,
`NUM_HIGH_STRIPS = 20` (sorting CM).

---

## `GemCluster.h`

`gem::GemCluster` — per-plane strip clustering and X/Y matching.
Algorithm details: `docs/technical_notes/gem_clustering/`.

| Method | Description |
|---|---|
| `void SetConfig(const ClusterConfig&)` / `const ClusterConfig &GetConfig() const` | |
| `void FormClusters(std::vector<StripHit>&, std::vector<StripCluster>&) const` | sorts hits by strip, groups, splits, charge-weights, and filters |
| `void CartesianReconstruct(const std::vector<StripCluster>& xc, const std::vector<StripCluster>& yc, std::vector<GEMHit>& hits, int det_id) const` | match X/Y to 2-D hits |

`ClusterConfig` is shared with `GemSystem.h`.

---

## `GemPedestal.h`

`gem::GemPedestal` — per-strip pedestal accumulator with the same CM
correction as the live GEM pipeline.

| Method | Description |
|---|---|
| `void Clear()` | drop accumulators |
| `void Accumulate(const ssp::SspEventData &evt)` | fold one event; APVs in online-ZS (nstrips ≠ 128) are skipped |
| `int NumApvs() const`, `int NumStrips() const` | |
| `int Write(const std::string &output_path) const` | JSON in the format `GemSystem::LoadPedestals` reads. Returns APV count or `<0` on failure. |

Non-copyable (`std::unique_ptr<Impl>`).

---

## `DetectorTransform.h`

Header-only POD. 3×3 rotation (Rx · Ry · Rz, intrinsic) + translation.
Cached internally; mutators auto-invalidate.

```cpp
DetectorTransform t;
t.set(x, y, z, rx_deg, ry_deg, rz_deg);   // preferred — set + cache rebuild
t.toLab(dx, dy, lx, ly, lz);              // 2-D point on detector plane
t.toLab(dx, dy, dz, lx, ly, lz);          // 3-D (e.g. shower depth)
t.labToLocal(lx, ly, lz, dx, dy, dz);     // inverse
t.rotate(dx, dy, ox, oy);                 // rotation only (drawing)
t.normal(nx, ny, nz);                     // surface normal in lab
const auto &m = t.matrix();               // cached 3×3 + translation
t.invalidate();                           // call after raw field writes
```

Direct field access (`x, y, z, rx, ry, rz`) is allowed but the caller
must `invalidate()` afterwards. `set()` does both in one shot.

---

## `RunInfoConfig.h`

Header-only loader for the run-period JSON
(`{ "configurations": [...] }`). `LoadRunConfig` picks the entry with
the largest `run_number` ≤ the requested run (or the latest if `< 0`).

### `prad2::RunConfig` fields

| Field | Default | |
|---|---|---|
| `energy_calib_file` | "" | HyCal calibration file path |
| `default_adc2mev` | 0.078 | fallback gain |
| `Ebeam` | 0 | MeV |
| `target_x/y/z` | 0,0,0 | mm |
| `hycal_x/y/z` | 0, 0, 6225 | mm |
| `hycal_tilt_x/y/z` | 0,0,0 | deg |
| `gem_x[4], gem_y[4], gem_z[4]` | (-252.8, +252.8, -252.8, +252.8) / 0 / (5423, 5384, 5823, 5784) | mm |
| `gem_tilt_x/y/z[4]` | 0 / 0 / (0, 180, 0, 180) | deg |
| `gem_pedestal_file`, `gem_common_mode_file` | "" | per-run paths |
| `hc_time_win_lo`, `hc_time_win_hi` | 100, 200 | ns |
| `matching_radius` | 15 mm | |
| `matching_use_square` | true | square vs circular cut |
| `gain_data_dir` | "gain_factor" | |
| `gain_ref_run` | 23915 | reference run for gain factors |

### Functions

```cpp
RunConfig LoadRunConfig(const std::string &path, int run_num);
bool      WriteRunConfig(const std::string &path, int run_num, const RunConfig &);
```

`WriteRunConfig` updates an existing run entry or appends a new one,
keeping the array sorted, and writes atomically via `<path>.tmp` rename.

---

## `PipelineBuilder.h`

`prad2::PipelineBuilder` — fluent helper that consolidates the
`Init` → `LoadCalibration` → `LoadPedestals` → `LoadCommonModeRange` →
`SetReconConfigs` sequence. Three callers (analysis scripts, the live
server, Python bindings) share the same wiring.

### `prad2::Pipeline` (build result)

| Field | Description |
|---|---|
| `daq_cfg` | `evc::DaqConfig` from prad2dec |
| `run_cfg` | `prad2::RunConfig` |
| `run_number` | resolved run (from `set_run_number*`) |
| `hycal` | initialized + calibrated `fdec::HyCalSystem` |
| `gem`   | initialized GEM with pedestals / CM / per-detector configs installed |
| `hycal_cluster_cfg` | `fdec::ClusterConfig` (caller wires into `HyCalCluster::SetConfig`) |
| `hycal_transform`, `gem_transforms[4]` | poses already `set()`, matrices cached |
| `hycal_pos_res[3]` | (A,B,C) face σ; default {2.6, 0, 0} |
| `gem_pos_res` | per-detector σ (mm) |
| `target_pos_res[3]` | (σx, σy, σz) at target |
| `gem_crate_remap` | derived from `daq_cfg.roc_tags[type=="gem"]` |
| `daq_config_path`, `recon_config_path`, `runinfo_path`, `hycal_map_path`, `gem_map_path`, `hycal_calib_path`, `gem_pedestal_path`, `gem_common_mode_path` | resolved absolute paths (empty if step skipped) |

Move-only — `HyCalSystem` and `GemSystem` hold large internal buffers.

### `PipelineBuilder` setters

All return `*this` for chaining. Empty paths fall back to defaults; non-empty
paths go through the resolver (default joins with `database_dir`).

```cpp
PipelineBuilder &set_database_dir(std::string);
PipelineBuilder &set_daq_config(std::string);
PipelineBuilder &set_loaded_daq_config(evc::DaqConfig);   // skip parse step
PipelineBuilder &set_recon_config(std::string);
PipelineBuilder &set_runinfo(std::string);
PipelineBuilder &set_hycal_map(std::string);
PipelineBuilder &set_gem_map(std::string);
PipelineBuilder &set_hycal_calib(std::string);
PipelineBuilder &set_gem_pedestal(std::string);
PipelineBuilder &set_gem_common_mode(std::string);

PipelineBuilder &set_run_number(int);
PipelineBuilder &set_run_number_from_evio(const std::string &evio_path);

PipelineBuilder &set_log_stream(std::ostream *);   // nullptr to silence
PipelineBuilder &set_log_pedestal_checksum(bool);

PipelineBuilder &set_path_resolver(std::function<std::string(const std::string&)>);

Pipeline build();   // throws std::runtime_error on hard failures
```

Hard failures (throws): missing daq_config / recon_config / runinfo /
hycal_map / gem_map. Soft failures (warns): missing calibration /
pedestal / common-mode files leave the corresponding `Pipeline` paths
empty.

---

## `EventData.h`

Plain-data structs that define the branch layout of the ROOT trees
produced by `analysis/Replay`. No ROOT headers — pure C++ types.

### Capacity constants

| Constant | Value |
|---|---|
| `kDscChannels` | 16 (DSC2 TRG + TDC channels per slot) |
| `kMaxChannels` | `fdec::MAX_ROCS · MAX_SLOTS · 16` |
| `kMaxGemStrips` | `ssp::MAX_MPDS · MAX_APVS_PER_MPD · APV_STRIP_SIZE` |
| `kMaxClusters` | 100 |
| `kMaxGemHits` | 400 |

### `prad2::ModuleType` (data-tree level)

`MOD_UNKNOWN = 0`, `MOD_PbGlass = 1`, `MOD_PbWO4 = 2`, `MOD_VETO = 3`,
`MOD_LMS = 4`. `module_id` ranges:

| Type | id range |
|---|---|
| MOD_PbGlass | 1..1156 |
| MOD_PbWO4   | 1001..2152 |
| MOD_VETO    | 3001..3004 (V1..V4) |
| MOD_LMS     | 3100..3103 (LMSPin, LMS1..3) |

### `RawEventData` ("events" tree)

Per-event header (`event_num`, `trigger_type`, `trigger_bits`,
`timestamp`), unified FADC250 channel array
(`nch`, `module_id`, `module_type`, `nsamples`, `samples`,
`gain_factor`), optional soft-analyzer peaks (`ped_*`, `npeaks`,
`peak_*`), optional firmware-mode peaks (`daq_npeaks`, `daq_peak_*`),
and GEM strip data (`gem_nch`, `mpd_crate`, `mpd_fiber`, `apv`,
`strip`, `ssp_samples`). `ssp_raw` is a `std::vector<uint32_t>` of
the raw 0xE10C SSP trigger bank words.

### `ReconEventData` ("recon" tree)

HyCal clusters (`n_clusters`, `cl_x/y/z/energy/time/nblocks/center/flag`),
HyCal↔GEM matches (`matchFlag`, `matchGEM[xyz]`), quick-access matched
pairs (`match_num`, `mHit_*`), GEM hits (`n_gem_hits`, `det_id`, `gem_*`),
Veto and LMS soft-peak summaries, and `ssp_raw`.

### `RawScalerData` ("scalers" tree)

DSC2 scaler row written when a physics event with a SYNC flag arrives.
Counts accumulate from GO; differencing two rows gives a windowed
live-time. Carries the wrapping physics event's `event_number` and
`ti_ticks` as join keys.

### `RawEpicsData` ("epics" tree)

One row per EPICS event (top-level tag 0x001F). Channel readings as
parallel `channel`/`value` vectors. `event_number_at_arrival` and
`ti_ticks_at_arrival` are the join keys to the events tree.

---

## `EventData_io.h`

Header-only TTree branch I/O — included only by ROOT-aware consumers
(replay tools, viewer's root data source, sim2replay). The library has
no link-time ROOT dependency.

### Writers (always set up the same branch list)

```cpp
SetRawWriteBranches(TTree*, RawEventData&, bool with_peaks);
SetReconWriteBranches(TTree*, ReconEventData&);
SetScalerWriteBranches(TTree*, RawScalerData&);
SetEpicsWriteBranches(TTree*, RawEpicsData&);
```

### Readers (skip missing branches; return optional-group flags)

```cpp
RawReadStatus   SetRawReadBranches(TTree*, RawEventData&);
ReconReadStatus SetReconReadBranches(TTree*, ReconEventData&);
SetScalerReadBranches(TTree*, RawScalerData&);   // void
SetEpicsReadBranches(TTree*, RawEpicsData&);     // void
```

`RawReadStatus { has_peaks, has_daq_peaks, has_gem, has_ssp_raw }`.
`ReconReadStatus { has_match_num, has_per_cl_match, has_veto, has_lms, has_ssp_raw }`.

`ssp_raw` is a `std::vector<uint32_t>` branch — readers must bind their
own held `vector<uint32_t>**` (see comment in header).

### Fillers — copy decoded `prad2dec` records into the side-tree PODs

```cpp
FillScalerRow(const dsc::DscEventData&, const psync::SyncInfo&,
              const fdec::EventInfo&, const evc::DaqConfig::DscScaler&,
              RawScalerData&);

FillEpicsRow(const epics::EpicsRecord&, RawEpicsData&);
```

Schema reference: `analysis/REPLAYED_DATA.md` is the human-readable
counterpart; this header is the executable schema.

---

## Dependencies

- [`prad2dec`](PRAD2DEC_API.md) — provides SSP / FADC / DSC / EPICS
  event-data types, `evc::DaqConfig`, and the `evc::EvChannel` reader.
- [nlohmann/json](https://github.com/nlohmann/json) — fetched
  automatically; used by every JSON-backed loader.
- ROOT — only needed by callers of `EventData_io.h`; the library itself
  has no ROOT link-time dependency.
