// bind_det.cpp — pybind11 bindings for prad2det (prad2py.det submodule).
//
// Phase 2a: GEM reconstruction.
//   prad2py.det.GemSystem            (detector hierarchy + processing)
//   prad2py.det.GemCluster           (clustering algorithm)
//   prad2py.det.ClusterConfig        (tuning knobs for GemCluster)
//   prad2py.det.StripHit / StripCluster / GEMHit   (per-event outputs)
//   prad2py.det.ApvConfig / PlaneConfig / DetectorConfig / ApvPedestal
//
// Phase 2b (HyCal) and 2c (DetectorTransform) plug in alongside these —
// each gets its own `py::class_<...>` block below.  EpicsStore moved to
// the dec submodule (prad2py.dec.EpicsStore) when it migrated to prad2dec.
//
// Usage sketch (matches the C++ driver in test/gem_dump.cpp):
//
//     from prad2py import dec, det
//     cfg = dec.load_daq_config()
//     ch  = dec.EvChannel(); ch.set_config(cfg); ch.open(path)
//
//     gsys = det.GemSystem()
//     gsys.init("database/gem_map.json")
//     gsys.load_pedestals("database/gem_ped.json")    # optional
//     gcl = det.GemCluster()
//
//     while ch.read() == dec.Status.success:
//         if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
//             continue
//         for i in range(ch.get_n_events()):
//             ch.select_event(i)
//             ssp = ch.gem()                  # SspEventData
//             gsys.clear()
//             gsys.process_event(ssp)
//             gsys.reconstruct(gcl)
//             for h in gsys.get_all_hits():   # list[GEMHit]
//                 print(h.det_id, h.x, h.y, h.x_charge, h.y_charge)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "GemSystem.h"
#include "GemCluster.h"
#include "GemPedestal.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "HyCalTimeCuts.h"
#include "HyCalRfOffsets.h"
#include "RfTime.h"
#include "DetectorTransform.h"
#include "PipelineBuilder.h"
#include "RunInfoConfig.h"
#include "SspData.h"

#include <cstring>
#include <cstdio>
#include <string>

namespace py = pybind11;

// -------------------------------------------------------------------------
// GEM bindings
// -------------------------------------------------------------------------
static void bind_gem(py::module_ &m)
{
    // --- stateless strip-mapping (shared with GemSystem::buildStripMap) -----
    // Exposed at module level so Python layout / diagnostic scripts
    // (gem/gem_strip_map.py, gem/gem_layout.py) hit the same C++
    // implementation the reconstruction uses — no more drifting duplicates.

    m.def("map_strip", &gem::MapStrip,
          py::arg("ch"), py::arg("plane_index"), py::arg("orient"),
          py::arg("pin_rotate")     = 0,
          py::arg("shared_pos")     = -1,
          py::arg("hybrid_board")   = true,
          py::arg("apv_channels")   = 128,
          py::arg("readout_center") = 32,
          "Map one APV25 channel index to the plane-wide strip number.  "
          "Pure function — same 6-step pipeline used by "
          "GemSystem::buildStripMap.");

    m.def("map_apv_strips", &gem::MapApvStrips,
          py::arg("plane_index"), py::arg("orient"),
          py::arg("pin_rotate")     = 0,
          py::arg("shared_pos")     = -1,
          py::arg("hybrid_board")   = true,
          py::arg("apv_channels")   = 128,
          py::arg("readout_center") = 32,
          "Compute plane-wide strip numbers for every channel of one APV.  "
          "Returns a Python list of length `apv_channels`.");

    // --- configuration leaves ------------------------------------------------

    py::class_<gem::ApvPedestal>(m, "ApvPedestal",
        "Per-strip pedestal offset + noise RMS.")
        .def(py::init<>())
        .def_readwrite("offset", &gem::ApvPedestal::offset)
        .def_readwrite("noise",  &gem::ApvPedestal::noise);

    py::class_<gem::PlaneConfig>(m, "PlaneConfig",
        "One plane (X or Y) of a GEM detector.")
        .def_readwrite("type",   &gem::PlaneConfig::type)
        .def_readwrite("size",   &gem::PlaneConfig::size)
        .def_readwrite("n_apvs", &gem::PlaneConfig::n_apvs)
        .def_readwrite("pitch",  &gem::PlaneConfig::pitch);

    py::class_<gem::DetectorConfig>(m, "DetectorConfig",
        "One GEM detector — name, id, type, and its two planes.")
        .def_readwrite("name",   &gem::DetectorConfig::name)
        .def_readwrite("id",     &gem::DetectorConfig::id)
        .def_readwrite("type",   &gem::DetectorConfig::type)
        .def_property_readonly("plane_x",
            [](const gem::DetectorConfig &d) -> const gem::PlaneConfig& { return d.planes[0]; },
            py::return_value_policy::reference_internal)
        .def_property_readonly("plane_y",
            [](const gem::DetectorConfig &d) -> const gem::PlaneConfig& { return d.planes[1]; },
            py::return_value_policy::reference_internal);

    py::class_<gem::ApvConfig>(m, "ApvConfig",
        "Per-APV configuration: DAQ address + detector mapping + strip-mapping "
        "params + per-strip pedestals.  One entry per APV in the gem_map.json.")
        .def_readwrite("crate_id",     &gem::ApvConfig::crate_id)
        .def_readwrite("mpd_id",       &gem::ApvConfig::mpd_id)
        .def_readwrite("adc_ch",       &gem::ApvConfig::adc_ch)
        .def_readwrite("det_id",       &gem::ApvConfig::det_id)
        .def_readwrite("plane_type",   &gem::ApvConfig::plane_type)
        .def_readwrite("orient",       &gem::ApvConfig::orient)
        .def_readwrite("plane_index",  &gem::ApvConfig::plane_index)
        .def_readwrite("det_pos",      &gem::ApvConfig::det_pos)
        .def_readwrite("pin_rotate",   &gem::ApvConfig::pin_rotate)
        .def_readwrite("shared_pos",   &gem::ApvConfig::shared_pos)
        .def_readwrite("hybrid_board", &gem::ApvConfig::hybrid_board)
        .def_readwrite("match",        &gem::ApvConfig::match)
        .def_readwrite("cm_range_min", &gem::ApvConfig::cm_range_min)
        .def_readwrite("cm_range_max", &gem::ApvConfig::cm_range_max)
        .def("pedestal",
            [](const gem::ApvConfig &c, int ch) -> const gem::ApvPedestal& {
                if (ch < 0 || ch >= 128)
                    throw py::index_error("APV channel out of range [0,128)");
                return c.pedestal[ch];
            },
            py::arg("ch"),
            py::return_value_policy::reference_internal,
            "Pedestal (offset + noise) for one of 128 APV channels.");

    // --- per-event output types ---------------------------------------------

    py::class_<gem::StripHit>(m, "StripHit",
        "One strip's zero-suppressed pulse — one entry per channel with charge "
        "above the ZS threshold.  `ts_adc` holds the 6 pedestal/CM-corrected "
        "time samples.")
        .def_readonly("strip",       &gem::StripHit::strip)
        .def_readonly("charge",      &gem::StripHit::charge)
        .def_readonly("max_timebin", &gem::StripHit::max_timebin)
        .def_readonly("position",    &gem::StripHit::position)
        .def_readonly("cross_talk",  &gem::StripHit::cross_talk)
        .def_property_readonly("ts_adc",
            [](const gem::StripHit &h) {
                // Copy into a fresh numpy array — safer than viewing into
                // the C++ vector, which could dangle if the StripHit is
                // moved / the owning GemSystem is cleared.
                auto arr = py::array_t<float>(
                    static_cast<py::ssize_t>(h.ts_adc.size()));
                if (!h.ts_adc.empty())
                    std::memcpy(arr.mutable_data(),
                                h.ts_adc.data(),
                                h.ts_adc.size() * sizeof(float));
                return arr;
            },
            "Time-sample ADC values after pedestal + common-mode correction "
            "(numpy float32, one entry per SSP time sample).");

    py::class_<gem::StripCluster>(m, "StripCluster",
        "One 1-D plane cluster: a group of consecutive StripHits with "
        "charge-weighted position.")
        .def_readonly("position",     &gem::StripCluster::position)
        .def_readonly("peak_charge",  &gem::StripCluster::peak_charge)
        .def_readonly("total_charge", &gem::StripCluster::total_charge)
        .def_readonly("max_timebin",  &gem::StripCluster::max_timebin)
        .def_readonly("cross_talk",   &gem::StripCluster::cross_talk)
        .def_readonly("hits",         &gem::StripCluster::hits);

    py::class_<gem::GEMHit>(m, "GEMHit",
        "2-D reconstructed GEM hit (X cluster × Y cluster match).")
        .def_readonly("x",             &gem::GEMHit::x)
        .def_readonly("y",             &gem::GEMHit::y)
        .def_readonly("z",             &gem::GEMHit::z)
        .def_readonly("det_id",        &gem::GEMHit::det_id)
        .def_readonly("x_charge",      &gem::GEMHit::x_charge)
        .def_readonly("y_charge",      &gem::GEMHit::y_charge)
        .def_readonly("x_peak",        &gem::GEMHit::x_peak)
        .def_readonly("y_peak",        &gem::GEMHit::y_peak)
        .def_readonly("x_max_timebin", &gem::GEMHit::x_max_timebin)
        .def_readonly("y_max_timebin", &gem::GEMHit::y_max_timebin)
        .def_readonly("x_size",        &gem::GEMHit::x_size)
        .def_readonly("y_size",        &gem::GEMHit::y_size)
        .def("__repr__", [](const gem::GEMHit &h) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<GEMHit det=%d x=%.3f y=%.3f z=%.3f qx=%.1f qy=%.1f>",
                h.det_id, h.x, h.y, h.z, h.x_charge, h.y_charge);
            return std::string(buf);
        });

    // --- GemCluster (configurable clustering algorithm) ---------------------

    py::class_<gem::ClusterConfig>(m, "ClusterConfig",
        "Tuning knobs for GemCluster.  Defaults reproduce the mpd_gem_view_ssp "
        "reconstruction chain.")
        .def(py::init<>())
        .def_readwrite("min_cluster_hits",   &gem::ClusterConfig::min_cluster_hits)
        .def_readwrite("max_cluster_hits",   &gem::ClusterConfig::max_cluster_hits)
        .def_readwrite("consecutive_thres",  &gem::ClusterConfig::consecutive_thres)
        .def_readwrite("split_thres",        &gem::ClusterConfig::split_thres)
        .def_readwrite("cross_talk_width",   &gem::ClusterConfig::cross_talk_width)
        .def_readwrite("charac_dists",       &gem::ClusterConfig::charac_dists)
        .def_readwrite("match_mode",         &gem::ClusterConfig::match_mode)
        .def_readwrite("match_adc_asymmetry",&gem::ClusterConfig::match_adc_asymmetry)
        .def_readwrite("match_time_diff",    &gem::ClusterConfig::match_time_diff)
        .def_readwrite("ts_period",          &gem::ClusterConfig::ts_period);

    py::class_<gem::GemCluster>(m, "GemCluster",
        "Strip clustering + X/Y cluster matching.  Used by GemSystem during "
        "Reconstruct(); can also be driven directly for custom analysis.")
        .def(py::init<>())
        .def("set_config", &gem::GemCluster::SetConfig, py::arg("cfg"))
        .def("get_config", &gem::GemCluster::GetConfig,
             py::return_value_policy::reference_internal)
        .def("form_clusters",
            [](const gem::GemCluster &self, std::vector<gem::StripHit> &hits) {
                std::vector<gem::StripCluster> out;
                self.FormClusters(hits, out);
                return out;
            },
            py::arg("hits"),
            "Group the supplied StripHit list into StripClusters.  Modifies "
            "`hits` in place (sorted by strip).")
        .def("cartesian_reconstruct",
            [](const gem::GemCluster &self,
               const std::vector<gem::StripCluster> &x_clusters,
               const std::vector<gem::StripCluster> &y_clusters,
               int det_id) {
                std::vector<gem::GEMHit> out;
                self.CartesianReconstruct(x_clusters, y_clusters, out, det_id);
                return out;
            },
            py::arg("x_clusters"), py::arg("y_clusters"),
            py::arg("det_id"),
            "Match X clusters with Y clusters on one detector to produce "
            "2-D GEMHits.");

    // --- GemSystem (the main entry point) -----------------------------------

    py::class_<gem::GemSystem>(m, "GemSystem",
        "PRad-II GEM detector system: loads gem_map.json / gem_ped.json, "
        "processes SspEventData (pedestal subtraction, common-mode correction, "
        "zero suppression, strip mapping), and hands off to GemCluster for "
        "2-D reconstruction.")
        .def(py::init<>())

        // initialization
        .def("init",
            [](gem::GemSystem &self, const std::string &path) {
                py::gil_scoped_release rel;
                self.Init(path);
            },
            py::arg("map_file"),
            "Load the detector hierarchy and APV mapping from a JSON file "
            "(typically database/gem_map.json).")
        .def("load_pedestals",
            [](gem::GemSystem &self, const std::string &path,
               const std::map<int, int> &crate_remap) {
                py::gil_scoped_release rel;
                self.LoadPedestals(path, crate_remap);
            },
            py::arg("ped_file"),
            py::arg("crate_remap") = std::map<int, int>{},
            "Load per-strip pedestal mean/RMS.  Required before ProcessEvent "
            "for real zero suppression — defaults keep all strips silent.  "
            "`crate_remap` (file-side hardware crate ID → logical crate ID "
            "in gem_map.json) defaults to identity; pass {tag: crate, ...} "
            "from daq_cfg.roc_tags when the pedestal file uses raw EVIO "
            "bank crate IDs (e.g. 146 → 1, 147 → 2 for PRad-II).")
        .def("load_common_mode_range",
            [](gem::GemSystem &self, const std::string &path,
               const std::map<int, int> &crate_remap) {
                py::gil_scoped_release rel;
                self.LoadCommonModeRange(path, crate_remap);
            },
            py::arg("cm_file"),
            py::arg("crate_remap") = std::map<int, int>{},
            "Optional per-APV common-mode suppression window file.  "
            "`crate_remap` semantics match load_pedestals.")

        // per-event processing
        .def("clear", &gem::GemSystem::Clear,
            "Reset per-event working buffers.  Call before every "
            "process_event().")
        .def("process_event",
            [](gem::GemSystem &self, const ssp::SspEventData &evt) {
                py::gil_scoped_release rel;
                self.ProcessEvent(evt);
            },
            py::arg("ssp_evt"),
            "Run pedestal + common-mode + zero-suppression over every APV in "
            "the given SspEventData.  Results feed reconstruct().")
        .def("reconstruct",
            [](gem::GemSystem &self, gem::GemCluster &cl) {
                py::gil_scoped_release rel;
                self.Reconstruct(cl);
            },
            py::arg("clusterer"),
            "Cluster the per-plane strip hits and match X/Y clusters into "
            "2-D GEMHits via the supplied GemCluster instance.")
        .def("set_recon_configs", &gem::GemSystem::SetReconConfigs,
             py::arg("cfgs"),
             "Install per-detector ClusterConfig (clustering + XY matching). "
             "Vector size is clamped/padded to GetNDetectors() — entry [d] "
             "is applied before clustering detector d in reconstruct().")
        .def("get_recon_configs", &gem::GemSystem::GetReconConfigs,
             py::return_value_policy::reference_internal,
             "Per-detector ClusterConfig list installed via set_recon_configs.")

        // accessors
        .def("get_n_detectors", &gem::GemSystem::GetNDetectors)
        .def("get_detectors",   &gem::GemSystem::GetDetectors,
             py::return_value_policy::reference_internal)
        .def("get_plane_hits", &gem::GemSystem::GetPlaneHits,
             py::arg("det"), py::arg("plane"),
             py::return_value_policy::reference_internal,
             "Per-plane StripHit list (plane=0 for X, 1 for Y).")
        .def("get_plane_clusters", &gem::GemSystem::GetPlaneClusters,
             py::arg("det"), py::arg("plane"),
             py::return_value_policy::reference_internal,
             "Per-plane StripCluster list (plane=0 for X, 1 for Y).  Populated "
             "after reconstruct().")
        .def("get_hits", &gem::GemSystem::GetHits, py::arg("det"),
             py::return_value_policy::reference_internal,
             "Per-detector reconstructed 2-D GEMHits.")
        .def("get_all_hits", &gem::GemSystem::GetAllHits,
             py::return_value_policy::reference_internal,
             "All GEMHits across every detector (flattened view).")

        // APV diagnostics
        .def("get_n_apvs", &gem::GemSystem::GetNApvs)
        .def("get_apv_config", &gem::GemSystem::GetApvConfig, py::arg("index"),
             py::return_value_policy::reference_internal)
        .def("find_apv_index", &gem::GemSystem::FindApvIndex,
             py::arg("crate"), py::arg("mpd"), py::arg("adc"),
             "O(1) DAQ→APV index lookup; -1 if no matching APV in the map.")
        .def("get_hole_x_offset", &gem::GemSystem::GetHoleXOffset,
             "Beam-hole X offset (mm) inferred from the `match` APVs in the map.")
        .def("get_active_extent", &gem::GemSystem::GetActiveExtent,
             py::arg("det_id"), py::arg("plane"),
             "Active strip extent for (det_id, plane) in detector-local "
             "coords (mm).  Returns (lo, hi) — tighter than PlaneConfig.size "
             "on the inner-edge side when split APVs (shared_pos) reuse "
             "strip numbers.  plane: 0=X, 1=Y.")

        // post-process-event diagnostics
        .def("is_channel_hit", &gem::GemSystem::IsChannelHit,
             py::arg("apv_index"), py::arg("ch"),
             "True if this strip survived zero suppression in the last "
             "process_event().")
        .def("has_apv_zs_hits", &gem::GemSystem::HasApvZsHits,
             py::arg("apv_index"),
             "True if any channel in this APV survived zero suppression.")
        .def("get_processed_adc", &gem::GemSystem::GetProcessedAdc,
             py::arg("apv_index"), py::arg("ch"), py::arg("ts"),
             "Pedestal + common-mode-corrected ADC for (APV, channel, time "
             "sample); valid after process_event().")
        .def("get_apv_frame",
            [](const gem::GemSystem &self, int apv_idx) {
                // Copy into a fresh (128, 6) float32 array in strip-major
                // order.  Underlying storage is time-major (raw[ts*128+ch])
                // — single-APV transposed copy is ~3 kB, negligible.
                constexpr int NS = ssp::APV_STRIP_SIZE;
                constexpr int NT = ssp::SSP_TIME_SAMPLES;
                py::array_t<float> out({NS, NT});
                auto buf = out.mutable_unchecked<2>();
                for (int ch = 0; ch < NS; ++ch)
                    for (int ts = 0; ts < NT; ++ts)
                        buf(ch, ts) = self.GetProcessedAdc(apv_idx, ch, ts);
                return out;
            },
            py::arg("apv_index"),
            "(128, 6) float32 array of pedestal+CM-subtracted ADC values for "
            "every strip of this APV.  Valid after process_event().")
        .def("get_apv_hit_mask",
            [](const gem::GemSystem &self, int apv_idx) {
                constexpr int NS = ssp::APV_STRIP_SIZE;
                py::array_t<bool> out(NS);
                auto buf = out.mutable_unchecked<1>();
                for (int ch = 0; ch < NS; ++ch)
                    buf(ch) = self.IsChannelHit(apv_idx, ch);
                return out;
            },
            py::arg("apv_index"),
            "(128,) bool array: True where the channel survived ZS in the "
            "last process_event().")
        .def("get_apv_ped_noise",
            [](const gem::GemSystem &self, int apv_idx) {
                constexpr int NS = ssp::APV_STRIP_SIZE;
                py::array_t<float> out(NS);
                auto buf = out.mutable_unchecked<1>();
                const auto &cfg = self.GetApvConfig(apv_idx);
                for (int ch = 0; ch < NS; ++ch)
                    buf(ch) = cfg.pedestal[ch].noise;
                return out;
            },
            py::arg("apv_index"),
            "(128,) float32 per-strip pedestal RMS (noise) from the "
            "loaded pedestal file.  Multiply by zero_sup_threshold to get "
            "the per-channel ZS cutoff.")

        // threshold knobs
        .def_property("common_mode_threshold",
            &gem::GemSystem::GetCommonModeThreshold,
            &gem::GemSystem::SetCommonModeThreshold)
        .def_property("zero_sup_threshold",
            &gem::GemSystem::GetZeroSupThreshold,
            &gem::GemSystem::SetZeroSupThreshold)
        .def_property_readonly("cross_talk_threshold",
            &gem::GemSystem::GetCrossTalkThreshold)
        // Per-strip filters loaded from gem_map.json — exposed read-only so
        // analysis scripts can verify they match the live monitor.
        .def_property_readonly("reject_first_timebin",
            &gem::GemSystem::GetRejectFirstTimebin)
        .def_property_readonly("reject_last_timebin",
            &gem::GemSystem::GetRejectLastTimebin)
        .def_property_readonly("min_peak_adc",
            &gem::GemSystem::GetMinPeakAdc)
        .def_property_readonly("min_sum_adc",
            &gem::GemSystem::GetMinSumAdc);

    // --- GemPedestal -------------------------------------------------------
    py::class_<gem::GemPedestal>(m, "GemPedestal",
        "Accumulate GEM per-strip pedestals from SSP raw data, then write "
        "a JSON file that GemSystem.load_pedestals can consume.  Same "
        "algorithm as `gem_dump -m ped` — both call this class.")
        .def(py::init<>())
        .def("clear", &gem::GemPedestal::Clear,
             "Drop all accumulated stats.")
        .def("accumulate",
            [](gem::GemPedestal &self, const ssp::SspEventData &evt) {
                py::gil_scoped_release rel;
                self.Accumulate(evt);
            },
            py::arg("ssp_event"),
            "Fold one event's SSP data into the running pedestal "
            "accumulators.  APVs with nstrips != 128 (online-ZS) are "
            "silently skipped.")
        .def_property_readonly("num_apvs", &gem::GemPedestal::NumApvs,
             "Number of APVs with at least one contribution.")
        .def_property_readonly("num_strips", &gem::GemPedestal::NumStrips,
             "Number of strips (across all APVs) with at least one "
             "contribution.")
        .def("write",
            [](const gem::GemPedestal &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.Write(path);
            },
            py::arg("output_path"),
            "Serialize the accumulated mean/RMS to JSON.  Returns the "
            "number of APVs written, or a negative value on I/O failure.");
}

// -------------------------------------------------------------------------
// HyCal bindings
// -------------------------------------------------------------------------
static void bind_hycal(py::module_ &m)
{
    // --- free helpers (module-level under prad2py.det) ----------------------
    m.def("shower_depth", &fdec::shower_depth,
          py::arg("center_id"), py::arg("energy_mev"),
          "Maximum-shower-development depth (mm) for an EM shower of energy "
          "`energy_mev` in the HyCal module with PrimEx id `center_id`.  "
          "Discriminates W/G modules by the PWO_ID0 boundary.  Returns 0 for "
          "energy ≤ 0.  Use to set the z coordinate of a HyCal cluster "
          "centroid before lab-frame projection.");

    // --- enums --------------------------------------------------------------

    py::enum_<fdec::ModuleType>(m, "ModuleType")
        .value("PbGlass", fdec::ModuleType::PbGlass)
        .value("PbWO4",   fdec::ModuleType::PbWO4)
        .value("LMS",     fdec::ModuleType::LMS)
        .value("Unknown", fdec::ModuleType::Unknown);

    py::enum_<fdec::Sector>(m, "Sector")
        .value("Center", fdec::Sector::Center)
        .value("Top",    fdec::Sector::Top)
        .value("Right",  fdec::Sector::Right)
        .value("Bottom", fdec::Sector::Bottom)
        .value("Left",   fdec::Sector::Left);

    // --- small leaf types ---------------------------------------------------

    py::class_<fdec::DaqAddr>(m, "DaqAddr",
        "HyCal module DAQ mapping: (crate, slot, channel).")
        .def_readwrite("crate",   &fdec::DaqAddr::crate)
        .def_readwrite("slot",    &fdec::DaqAddr::slot)
        .def_readwrite("channel", &fdec::DaqAddr::channel);

    py::class_<fdec::NeighborInfo>(m, "NeighborInfo",
        "One cross-sector neighbor of a HyCal module.")
        .def_readonly("index", &fdec::NeighborInfo::index)
        .def_readonly("dx",    &fdec::NeighborInfo::dx)
        .def_readonly("dy",    &fdec::NeighborInfo::dy)
        .def_readonly("dist",  &fdec::NeighborInfo::dist);

    py::class_<fdec::SectorInfo>(m, "SectorInfo",
        "HyCal sector (Center / Top / Right / Bottom / Left) geometry.")
        .def_readonly("id",      &fdec::SectorInfo::id)
        .def_readonly("mtype",   &fdec::SectorInfo::mtype)
        .def_readonly("msize_x", &fdec::SectorInfo::msize_x)
        .def_readonly("msize_y", &fdec::SectorInfo::msize_y)
        .def("boundary",
            [](const fdec::SectorInfo &s) {
                double x1, y1, x2, y2;
                s.get_boundary(x1, y1, x2, y2);
                return py::make_tuple(x1, y1, x2, y2);
            },
            "Return (x_min, y_min, x_max, y_max) of this sector's rectangle.");

    // --- Module -------------------------------------------------------------

    py::class_<fdec::Module>(m, "Module",
        "One HyCal module: identity, geometry, DAQ mapping, calibration, and "
        "pre-computed neighbor list.")
        .def_readonly("name",   &fdec::Module::name)
        .def_readonly("id",     &fdec::Module::id)
        .def_readonly("index",  &fdec::Module::index)
        .def_readonly("type",   &fdec::Module::type)
        .def_readonly("x",      &fdec::Module::x)
        .def_readonly("y",      &fdec::Module::y)
        .def_readonly("size_x", &fdec::Module::size_x)
        .def_readonly("size_y", &fdec::Module::size_y)
        .def_readonly("flag",   &fdec::Module::flag)
        .def_readonly("sector", &fdec::Module::sector)
        .def_readonly("row",    &fdec::Module::row)
        .def_readonly("column", &fdec::Module::column)
        .def_readonly("daq",    &fdec::Module::daq)
        .def_readonly("cal_factor",      &fdec::Module::cal_factor)
        .def_readonly("cal_base_energy", &fdec::Module::cal_base_energy)
        .def_readonly("cal_non_linear_1",  &fdec::Module::cal_non_linear_1)
        .def_readonly("cal_non_linear_2",  &fdec::Module::cal_non_linear_2)
        .def("energize", &fdec::Module::energize, py::arg("adc"),
             "Convert a pedestal-subtracted ADC value to MeV, including the "
             "non-linear correction term.")
        .def("is_pwo4",  &fdec::Module::is_pwo4)
        .def("is_glass", &fdec::Module::is_glass)
        .def("is_hycal", &fdec::Module::is_hycal)
        .def_property_readonly("neighbors",
            [](const fdec::Module &m) {
                py::list out;
                for (int i = 0; i < m.neighbor_count; ++i)
                    out.append(m.neighbors[i]);
                return out;
            },
            "Pre-computed cross-sector neighbor list (NeighborInfo[]).  "
            "Same-sector neighbors are resolved via the sector grid instead.")
        .def("__repr__", [](const fdec::Module &m) {
            char buf[96];
            std::snprintf(buf, sizeof(buf),
                "<Module %s id=%d idx=%d x=%.1f y=%.1f>",
                m.name.c_str(), m.id, m.index, m.x, m.y);
            return std::string(buf);
        });

    // --- HyCalSystem --------------------------------------------------------

    py::class_<fdec::HyCalSystem>(m, "HyCalSystem",
        "HyCal detector geometry + DAQ map + calibration.  Initialized once "
        "per job and then immutable — no per-event state lives here.")
        .def(py::init<>())
        .def("init",
            [](fdec::HyCalSystem &self, const std::string &map_path) {
                py::gil_scoped_release rel;
                return self.Init(map_path);
            },
            py::arg("map_path"),
            "Load HyCal module geometry + DAQ map from hycal_map.json.  "
            "Returns True on success.")
        .def("load_calibration",
            [](fdec::HyCalSystem &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.LoadCalibration(path);
            },
            py::arg("calib_path"),
            "Load per-module calibration constants from a JSON file.  "
            "Returns the number of modules matched, or -1 on error.")

        .def("module_count", &fdec::HyCalSystem::module_count)
        .def("module",
            static_cast<const fdec::Module &(fdec::HyCalSystem::*)(int) const>(
                &fdec::HyCalSystem::module),
            py::arg("index"),
            py::return_value_policy::reference_internal,
            "Module by array index (0 .. module_count()-1).")
        .def("module_by_name",
            [](const fdec::HyCalSystem &self, const std::string &name)
                -> const fdec::Module* {
                return self.module_by_name(name);
            },
            py::arg("name"),
            py::return_value_policy::reference_internal,
            "Module by name (e.g. 'W735'); None if not found.")
        .def("module_by_id",
            [](const fdec::HyCalSystem &self, int id) -> const fdec::Module* {
                return self.module_by_id(id);
            },
            py::arg("primex_id"),
            py::return_value_policy::reference_internal,
            "Module by PrimEx id (G: 1-576, W: 1001-2152); None if not found.")
        .def("module_by_daq",
            [](const fdec::HyCalSystem &self, int crate, int slot, int ch)
                -> const fdec::Module* {
                return self.module_by_daq(crate, slot, ch);
            },
            py::arg("crate"), py::arg("slot"), py::arg("channel"),
            py::return_value_policy::reference_internal,
            "Module by DAQ address; None if not found.")

        .def("get_calib_constant", &fdec::HyCalSystem::GetCalibConstant,
             py::arg("primex_id"))
        .def("set_calib_constant", &fdec::HyCalSystem::SetCalibConstant,
             py::arg("primex_id"), py::arg("factor"))
        .def("print_calib_constants", &fdec::HyCalSystem::PrintCalibConstants,
             py::arg("output_file"))

        .def("set_position_resolution_params",
             &fdec::HyCalSystem::SetPositionResolutionParams,
             py::arg("A"), py::arg("B"), py::arg("C"),
             "Set the [A, B, C] coefficients of the HyCal-face position "
             "resolution formula. See PositionResolution(E).")
        .def("position_resolution",
             &fdec::HyCalSystem::PositionResolution,
             py::arg("energy_mev"),
             "sigma(E) at HyCal face (mm): "
             "sqrt((A/sqrt(E_GeV))^2 + (B/E_GeV)^2 + C^2).")

        .def("sector_info",
            &fdec::HyCalSystem::sector_info, py::arg("sector"),
            py::return_value_policy::reference_internal)
        .def("get_sector_id", &fdec::HyCalSystem::get_sector_id,
             py::arg("x"), py::arg("y"),
             "Return the sector id containing the (x, y) point, or -1 if none.")

        .def("qdist",
            [](const fdec::HyCalSystem &self,
               double x1, double y1, int s1,
               double x2, double y2, int s2) {
                double dx, dy;
                self.qdist(x1, y1, s1, x2, y2, s2, dx, dy);
                return py::make_tuple(dx, dy);
            },
            py::arg("x1"), py::arg("y1"), py::arg("s1"),
            py::arg("x2"), py::arg("y2"), py::arg("s2"),
            "Quantized (dx, dy) distance between two points across "
            "(potentially different) sectors.")
        .def("qdist_modules",
            [](const fdec::HyCalSystem &self,
               const fdec::Module &m1, const fdec::Module &m2) {
                double dx, dy;
                self.qdist(m1, m2, dx, dy);
                return py::make_tuple(dx, dy);
            },
            py::arg("m1"), py::arg("m2"),
            "Quantized (dx, dy) distance between two modules.")

        .def("for_each_neighbor",
            [](const fdec::HyCalSystem &self, int module_index,
               bool include_corners, const py::function &fn) {
                self.for_each_neighbor(module_index, include_corners,
                    [&](int ni) { fn(ni); });
            },
            py::arg("module_index"), py::arg("include_corners"),
            py::arg("fn"),
            "Walk every neighbor of `module_index`, same-sector (grid) + "
            "cross-sector (pre-computed list), calling fn(neighbor_index).")

        .def_static("name_to_id",  &fdec::HyCalSystem::name_to_id,
                    py::arg("name"))
        .def_static("id_to_name",  &fdec::HyCalSystem::id_to_name,
                    py::arg("id"))
        .def_static("parse_type",  &fdec::HyCalSystem::parse_type,
                    py::arg("type_string"));

    // --- clustering types ---------------------------------------------------

    py::class_<fdec::ClusterConfig>(m, "HyCalClusterConfig",
        "HyCal island-clustering tuning parameters.  Note: named "
        "`HyCalClusterConfig` in Python so it doesn't collide with the GEM "
        "`ClusterConfig`.")
        .def(py::init<>())
        .def_readwrite("min_module_energy",  &fdec::ClusterConfig::min_module_energy)
        .def_readwrite("min_center_energy",  &fdec::ClusterConfig::min_center_energy)
        .def_readwrite("min_cluster_energy", &fdec::ClusterConfig::min_cluster_energy)
        .def_readwrite("min_cluster_size",   &fdec::ClusterConfig::min_cluster_size)
        .def_readwrite("corner_conn",        &fdec::ClusterConfig::corner_conn)
        .def_readwrite("non_linear_corr",    &fdec::ClusterConfig::non_linear_corr)
        .def_readwrite("split_iter",         &fdec::ClusterConfig::split_iter)
        .def_readwrite("least_split",        &fdec::ClusterConfig::least_split)
        .def_readwrite("log_weight_thres",   &fdec::ClusterConfig::log_weight_thres)
        .def_readwrite("seed_time_window",   &fdec::ClusterConfig::seed_time_window);

    py::class_<fdec::ModuleHit>(m, "ModuleHit",
        "One HyCal module hit fed into the clusterer.")
        .def(py::init<>())
        .def_readwrite("index",  &fdec::ModuleHit::index)
        .def_readwrite("energy", &fdec::ModuleHit::energy)
        .def_readwrite("time",   &fdec::ModuleHit::time);

    py::class_<fdec::ModuleCluster>(m, "ModuleCluster",
        "Module-level cluster: seed + constituent hits + total energy + flags.")
        .def_readonly("center", &fdec::ModuleCluster::center)
        .def_readonly("hits",   &fdec::ModuleCluster::hits)
        .def_readonly("energy", &fdec::ModuleCluster::energy)
        .def_readonly("flag",   &fdec::ModuleCluster::flag);

    py::class_<fdec::ClusterHit>(m, "ClusterHit",
        "Reconstructed HyCal cluster with (x, y) position in mm.")
        .def_readonly("center_id", &fdec::ClusterHit::center_id)
        .def_readonly("x",         &fdec::ClusterHit::x)
        .def_readonly("y",         &fdec::ClusterHit::y)
        .def_readonly("energy",    &fdec::ClusterHit::energy)
        .def_readonly("time",      &fdec::ClusterHit::time)
        .def_readonly("nblocks",   &fdec::ClusterHit::nblocks)
        .def_readonly("npos",      &fdec::ClusterHit::npos)
        .def_readonly("flag",      &fdec::ClusterHit::flag)
        .def("__repr__", [](const fdec::ClusterHit &c) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<ClusterHit center_id=%d x=%.2f y=%.2f E=%.1f nblocks=%d>",
                c.center_id, c.x, c.y, c.energy, c.nblocks);
            return std::string(buf);
        });

    py::class_<fdec::HyCalCluster::RecoResult>(m, "HyCalRecoResult",
        "Paired ModuleCluster + reconstructed ClusterHit from "
        "HyCalCluster.reconstruct_matched().  `cluster` is a reference back "
        "into HyCalCluster's internal list — invalidated on clear().")
        .def_property_readonly("cluster",
            [](const fdec::HyCalCluster::RecoResult &r) -> const fdec::ModuleCluster& {
                return *r.cluster;
            },
            py::return_value_policy::reference_internal)
        .def_readonly("hit", &fdec::HyCalCluster::RecoResult::hit);

    py::class_<fdec::HyCalCluster::SeedNeighborTiming>(m, "SeedNeighborTiming",
        "One (seed, neighbour-pulse) row from HyCalCluster.collect_neighbor_timing(). "
        "Used to histogram dt vs. spatial distance / energy on real data and "
        "pick a value for HyCalClusterConfig.seed_time_window.")
        .def_readonly("seed_module",     &fdec::HyCalCluster::SeedNeighborTiming::seed_module)
        .def_readonly("neighbor_module", &fdec::HyCalCluster::SeedNeighborTiming::neighbor_module)
        .def_readonly("seed_time",       &fdec::HyCalCluster::SeedNeighborTiming::seed_time)
        .def_readonly("neighbor_time",   &fdec::HyCalCluster::SeedNeighborTiming::neighbor_time)
        .def_readonly("dt",              &fdec::HyCalCluster::SeedNeighborTiming::dt)
        .def_readonly("seed_energy",     &fdec::HyCalCluster::SeedNeighborTiming::seed_energy)
        .def_readonly("neighbor_energy", &fdec::HyCalCluster::SeedNeighborTiming::neighbor_energy)
        .def_readonly("dx_q",            &fdec::HyCalCluster::SeedNeighborTiming::dx_q)
        .def_readonly("dy_q",            &fdec::HyCalCluster::SeedNeighborTiming::dy_q);

    // --- HyCalCluster -------------------------------------------------------

    py::class_<fdec::HyCalCluster>(m, "HyCalCluster",
        "Island clustering on top of a HyCalSystem: accepts (module_index, "
        "energy) hits per event, groups connected modules, splits multi-"
        "maximum groups, and produces ClusterHit positions via log-weighted "
        "center-of-gravity.")
        .def(py::init<const fdec::HyCalSystem &>(), py::arg("system"),
             py::keep_alive<1, 2>(),
             "Construct with a reference to an initialized HyCalSystem — "
             "the system must outlive this clusterer.")
        .def("set_config", &fdec::HyCalCluster::SetConfig, py::arg("cfg"))
        .def("get_config", &fdec::HyCalCluster::GetConfig,
             py::return_value_policy::reference_internal)
        .def("clear", &fdec::HyCalCluster::Clear,
             "Reset per-event state.  Call before add_hit()s for each event.")
        .def("add_hit", &fdec::HyCalCluster::AddHit,
             py::arg("module_index"), py::arg("energy"), py::arg("time"),
             "Add a hit for the module at `module_index` with the given "
             "calibrated energy (MeV) and time (ns).")
        .def("form_clusters",
            [](fdec::HyCalCluster &self) {
                py::gil_scoped_release rel;
                self.FormClusters();
            },
            "Run the island grouping + splitting algorithm over the hits "
            "accumulated via add_hit().")
        .def("reconstruct_hits",
            [](const fdec::HyCalCluster &self) {
                std::vector<fdec::ClusterHit> out;
                { py::gil_scoped_release rel; self.ReconstructHits(out); }
                return out;
            },
            "Return the list of ClusterHits (x, y, energy, nblocks, …) for "
            "every cluster that passed min_cluster_energy / min_cluster_size.")
        .def("reconstruct_matched",
            [](const fdec::HyCalCluster &self) {
                std::vector<fdec::HyCalCluster::RecoResult> out;
                { py::gil_scoped_release rel; self.ReconstructMatched(out); }
                return out;
            },
            "Return [(ModuleCluster, ClusterHit)] pairs for every cluster "
            "that passed thresholds — avoids the fragile parallel-iteration "
            "between get_clusters() and reconstruct_hits().")
        .def("get_clusters", &fdec::HyCalCluster::GetClusters,
             py::return_value_policy::reference_internal,
             "Low-level module-level clusters (ModuleCluster[]).  Only valid "
             "between form_clusters() and the next clear().")
        .def("collect_neighbor_timing",
            [](const fdec::HyCalCluster &self, double max_quantized_dist) {
                std::vector<fdec::HyCalCluster::SeedNeighborTiming> out;
                { py::gil_scoped_release rel;
                  self.CollectNeighborTiming(out, max_quantized_dist); }
                return out;
            },
            py::arg("max_quantized_dist") = 5.0,
            "Identify seed candidates (largest pulse satisfying "
            "min_center_energy that hasn't already seeded another cluster in "
            "this scan) and return one SeedNeighborTiming row per neighbouring "
            "pulse within `max_quantized_dist` module units of the seed — "
            "WITHOUT applying any timing cut and WITHOUT consuming neighbour "
            "pulses.  Use to inform HyCalClusterConfig.seed_time_window.");
}

// -------------------------------------------------------------------------
// Helper bindings (Phase 2c)
// -------------------------------------------------------------------------
static void bind_transform(py::module_ &m)
{
    py::class_<DetectorTransform::Matrix>(m, "TransformMatrix",
        "Pre-computed 3×3 rotation + translation matrix cached inside a "
        "DetectorTransform.  Row-major: rIJ.")
        .def_readonly("r00", &DetectorTransform::Matrix::r00)
        .def_readonly("r01", &DetectorTransform::Matrix::r01)
        .def_readonly("r02", &DetectorTransform::Matrix::r02)
        .def_readonly("r10", &DetectorTransform::Matrix::r10)
        .def_readonly("r11", &DetectorTransform::Matrix::r11)
        .def_readonly("r12", &DetectorTransform::Matrix::r12)
        .def_readonly("r20", &DetectorTransform::Matrix::r20)
        .def_readonly("r21", &DetectorTransform::Matrix::r21)
        .def_readonly("r22", &DetectorTransform::Matrix::r22)
        .def_readonly("tx",  &DetectorTransform::Matrix::tx)
        .def_readonly("ty",  &DetectorTransform::Matrix::ty)
        .def_readonly("tz",  &DetectorTransform::Matrix::tz);

    // Property setter that invalidates the cached rotation matrix on
    // every write, so `t.x = 5` does the obvious thing — the next to_lab
    // call rebuilds with the new field.  Macro keeps the six bindings
    // identical without the boilerplate.
#define PRAD2_BIND_TRANSFORM_AXIS(NAME, MEMBER)                            \
    .def_property(NAME,                                                    \
        [](const DetectorTransform &t) { return t.MEMBER; },               \
        [](DetectorTransform &t, float v) {                                \
            t.MEMBER = v; t.invalidate();                                  \
        })

    py::class_<DetectorTransform>(m, "DetectorTransform",
        "Planar detector pose (origin + tilts in degrees).  to_lab(x, y[, z]) "
        "returns a lab-frame 3-vector; rotate(x, y) skips the translation. "
        "Field setters auto-invalidate the cached rotation matrix; prefer "
        "set(x, y, z, rx, ry, rz) when writing all six at once.")
        .def(py::init<>())
        PRAD2_BIND_TRANSFORM_AXIS("x",  x)
        PRAD2_BIND_TRANSFORM_AXIS("y",  y)
        PRAD2_BIND_TRANSFORM_AXIS("z",  z)
        PRAD2_BIND_TRANSFORM_AXIS("rx", rx)
        PRAD2_BIND_TRANSFORM_AXIS("ry", ry)
        PRAD2_BIND_TRANSFORM_AXIS("rz", rz)
        .def("set", &DetectorTransform::set,
             py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("rx"), py::arg("ry"), py::arg("rz"),
             "Set translation (mm) and tilts (degrees), then rebuild the "
             "cached rotation matrix.  Preferred over per-field writes.")
        .def("invalidate", &DetectorTransform::invalidate,
             "Mark the cached rotation matrix stale; the next to_lab / "
             "rotate / matrix call will rebuild it.")
        .def("prepare", &DetectorTransform::prepare,
             "Force matrix precomputation (idempotent — first call after "
             "construction or invalidate()).  Called implicitly by "
             "to_lab / rotate / matrix().")
        .def("to_lab",
            [](const DetectorTransform &self,
               float dx, float dy, float dz) {
                float lx, ly, lz;
                self.toLab(dx, dy, dz, lx, ly, lz);
                return py::make_tuple(lx, ly, lz);
            },
            py::arg("dx"), py::arg("dy"), py::arg("dz") = 0.0f,
            "Map a local point to lab-frame (x, y, z).  dz defaults to 0 "
            "for planar GEM hits; pass shower depth for HyCal clusters.")
        .def("lab_to_local",
            [](const DetectorTransform &self,
               float lx, float ly, float lz) {
                float dx, dy, dz;
                self.labToLocal(lx, ly, lz, dx, dy, dz);
                return py::make_tuple(dx, dy, dz);
            },
            py::arg("lx"), py::arg("ly"), py::arg("lz"),
            "Inverse of to_lab: map a lab-frame point to detector-local "
            "(dx, dy, dz).  Useful for histogramming predicted hit "
            "positions on the detector plane.")
        .def("rotate",
            [](const DetectorTransform &self, float dx, float dy) {
                float ox, oy;
                self.rotate(dx, dy, ox, oy);
                return py::make_tuple(ox, oy);
            },
            py::arg("dx"), py::arg("dy"),
            "Rotate a (dx, dy) vector by the detector tilts — no "
            "translation.  Handy for drawing in detector-local coords.")
        .def("matrix", &DetectorTransform::matrix,
             py::return_value_policy::reference_internal);

#undef PRAD2_BIND_TRANSFORM_AXIS
}

// -------------------------------------------------------------------------
// PipelineBuilder bindings — one-stop wiring of HyCal + GEM detectors.
// -------------------------------------------------------------------------
//
// Mirrors the C++ side at prad2det/include/PipelineBuilder.h: the builder
// loads daq_config + reconstruction_config + runinfo, initializes both
// detectors with calibration / pedestals / per-detector cluster configs,
// and constructs the lab-frame DetectorTransforms.  Replaces ~100 LOC of
// JSON parsing + Init/Load* orchestration in analysis/pyscripts/_common.py.
//
// Usage (Python):
//
//   from prad2py import det
//   p = (det.PipelineBuilder()
//          .set_run_number_from_evio(evio_path)
//          .build())
//   hc = det.HyCalCluster(p.hycal); hc.set_config(p.hycal_cluster_cfg)
//   gem_cl = det.GemCluster()
//   wa = dec.WaveAnalyzer(p.daq_cfg.wave_cfg)
//   # ... per-event loop using p.hycal, p.gem, p.hycal_transform, ...
static void bind_pipeline(py::module_ &m)
{
    // RunConfig — geometry + per-run calibration metadata produced by
    // LoadRunConfig().  Bound as a small read-only view; Python builds
    // its own lab transforms via Pipeline.hycal_transform / gem_transforms
    // rather than touching these fields directly.
    py::class_<prad2::RunConfig>(m, "RunConfig",
        "Run-period geometry + calibration paths.  Read off Pipeline.run_cfg.")
        .def(py::init<>())
        .def_readonly("energy_calib_file",    &prad2::RunConfig::energy_calib_file)
        .def_readonly("default_adc2mev",      &prad2::RunConfig::default_adc2mev)
        .def_readonly("Ebeam",                &prad2::RunConfig::Ebeam)
        .def_readonly("target_x",             &prad2::RunConfig::target_x)
        .def_readonly("target_y",             &prad2::RunConfig::target_y)
        .def_readonly("target_z",             &prad2::RunConfig::target_z)
        .def_readonly("hycal_x",              &prad2::RunConfig::hycal_x)
        .def_readonly("hycal_y",              &prad2::RunConfig::hycal_y)
        .def_readonly("hycal_z",              &prad2::RunConfig::hycal_z)
        .def_readonly("hycal_tilt_x",         &prad2::RunConfig::hycal_tilt_x)
        .def_readonly("hycal_tilt_y",         &prad2::RunConfig::hycal_tilt_y)
        .def_readonly("hycal_tilt_z",         &prad2::RunConfig::hycal_tilt_z)
        .def_property_readonly("gem_x",
            [](const prad2::RunConfig &r) {
                return std::vector<float>(r.gem_x, r.gem_x + 4);
            })
        .def_property_readonly("gem_y",
            [](const prad2::RunConfig &r) {
                return std::vector<float>(r.gem_y, r.gem_y + 4);
            })
        .def_property_readonly("gem_z",
            [](const prad2::RunConfig &r) {
                return std::vector<float>(r.gem_z, r.gem_z + 4);
            })
        .def_readonly("gem_pedestal_file",    &prad2::RunConfig::gem_pedestal_file)
        .def_readonly("gem_common_mode_file", &prad2::RunConfig::gem_common_mode_file)
        .def_readonly("hc_time_win_lo",       &prad2::RunConfig::hc_time_win_lo)
        .def_readonly("hc_time_win_hi",       &prad2::RunConfig::hc_time_win_hi)
        .def_readonly("hycal_time_cut_file",  &prad2::RunConfig::hycal_time_cut_file);

    // HyCalTimeCuts — per-module HyCal peak-time window with a default
    // fallback (uniform when no per-module file is present).  Sized to
    // hycal.module_count(); use `at(mod.index)` in the per-event loop.
    py::class_<prad2::HyCalTimeCuts::Window>(m, "HyCalTimeWindow")
        .def_readonly("lo", &prad2::HyCalTimeCuts::Window::lo)
        .def_readonly("hi", &prad2::HyCalTimeCuts::Window::hi)
        .def("__repr__", [](const prad2::HyCalTimeCuts::Window &w) {
            return "<HyCalTimeWindow lo=" + std::to_string(w.lo)
                 + " hi=" + std::to_string(w.hi) + ">";
        });

    py::class_<prad2::HyCalTimeCuts>(m, "HyCalTimeCuts",
        "Per-module HyCal peak-time window.  Read off Pipeline.hycal_time_cuts.")
        .def(py::init<>())
        .def_readonly("default_lo",  &prad2::HyCalTimeCuts::default_lo)
        .def_readonly("default_hi",  &prad2::HyCalTimeCuts::default_hi)
        .def_readonly("n_overrides", &prad2::HyCalTimeCuts::n_overrides)
        .def("at",        &prad2::HyCalTimeCuts::at,        py::arg("module_index"))
        .def("in_window", &prad2::HyCalTimeCuts::in_window,
             py::arg("module_index"), py::arg("t"));

    // HyCalRfOffsets — per-module HyCal→RF time offset, with a default
    // fallback.  Used by Replay to fill ReconEventData.cl_dt_rf.
    py::class_<prad2::HyCalRfOffsets>(m, "HyCalRfOffsets",
        "Per-module HyCal→RF time offset (ns).  Read off "
        "Pipeline.hycal_rf_offsets.")
        .def(py::init<>())
        .def_readonly("default_off", &prad2::HyCalRfOffsets::default_off)
        .def_readonly("n_overrides", &prad2::HyCalRfOffsets::n_overrides)
        .def("at",    &prad2::HyCalRfOffsets::at,    py::arg("module_index"))
        .def("apply", &prad2::HyCalRfOffsets::apply,
             py::arg("module_index"), py::arg("folded_dt"),
             "Subtract per-module offset from folded Δt and re-fold onto "
             "(-T_RF/2, T_RF/2].");

    // RF folding helpers — module-level so analysis scripts can use them
    // without instantiating a Pipeline.
    m.attr("RF_PERIOD_NS") = prad2::RF_PERIOD_NS;
    m.attr("RF_DIVIDER")   = prad2::RF_DIVIDER;
    m.attr("RF_DIV_NS")    = prad2::RF_DIV_NS;
    m.def("fold_rf_delta", &prad2::FoldRfDelta, py::arg("dt_ns"),
          "Fold a raw Δt onto (-T_RF/2, T_RF/2].  NaN passes through "
          "unchanged.");
    m.def("cluster_delta_rf", &prad2::ClusterDeltaRf,
          py::arg("t_ref_ns"), py::arg("rf"), py::arg("use_b") = false,
          "Compute the folded Δt between a reference time (typically "
          "ClusterHit.time) and the nearest RF tick on channel A "
          "(B if use_b=True).  Returns NaN when the channel has no hits.");

    // Pipeline — result bundle.  Owned by Python.  HyCalSystem and GemSystem
    // are exposed by reference (their lifetime is tied to Pipeline) so the
    // user can pass them directly into HyCalCluster / GEM reconstruction.
    py::class_<prad2::Pipeline>(m, "Pipeline",
        "Wired-up HyCal + GEM detectors plus the resolved configs.  Returned "
        "by PipelineBuilder.build().  Move-only on the C++ side; on the "
        "Python side a Pipeline owns its detectors and you typically hold one "
        "per analysis script.")
        .def_readonly("daq_cfg",            &prad2::Pipeline::daq_cfg,
                      py::return_value_policy::reference_internal)
        .def_readonly("run_cfg",            &prad2::Pipeline::run_cfg,
                      py::return_value_policy::reference_internal)
        .def_readonly("run_number",         &prad2::Pipeline::run_number)
        .def_readonly("hycal",              &prad2::Pipeline::hycal,
                      py::return_value_policy::reference_internal)
        .def_readonly("gem",                &prad2::Pipeline::gem,
                      py::return_value_policy::reference_internal)
        .def_readonly("hycal_cluster_cfg",  &prad2::Pipeline::hycal_cluster_cfg)
        .def_readonly("hycal_time_cuts",    &prad2::Pipeline::hycal_time_cuts,
                      py::return_value_policy::reference_internal,
                      "Per-module HyCal peak-time window (sized to hycal.module_count()).")
        .def_readonly("hycal_rf_offsets",   &prad2::Pipeline::hycal_rf_offsets,
                      py::return_value_policy::reference_internal,
                      "Per-module HyCal→RF time offset, ns.")
        .def_readonly("hycal_transform",    &prad2::Pipeline::hycal_transform,
                      py::return_value_policy::reference_internal)
        .def_property_readonly("gem_transforms",
            [](const prad2::Pipeline &p) {
                return std::vector<DetectorTransform>(p.gem_transforms.begin(),
                                                      p.gem_transforms.end());
            },
            "Per-detector lab transforms (list of 4 DetectorTransform).")
        .def_property_readonly("hycal_pos_res",
            [](const prad2::Pipeline &p) {
                return std::vector<float>(p.hycal_pos_res.begin(),
                                          p.hycal_pos_res.end());
            },
            "[A, B, C] coefficients of HyCal-face position resolution.")
        .def_readonly("gem_pos_res",        &prad2::Pipeline::gem_pos_res)
        .def_property_readonly("target_pos_res",
            [](const prad2::Pipeline &p) {
                return std::vector<float>(p.target_pos_res.begin(),
                                          p.target_pos_res.end());
            },
            "[sigma_x, sigma_y, sigma_z] of the target gas distribution.")
        .def_readonly("gem_crate_remap",    &prad2::Pipeline::gem_crate_remap,
                      "Hardware crate ID -> logical crate ID for GEM (from "
                      "daq_cfg.roc_tags entries with type=='gem').")
        .def_readonly("daq_config_path",    &prad2::Pipeline::daq_config_path)
        .def_readonly("recon_config_path",  &prad2::Pipeline::recon_config_path)
        .def_readonly("runinfo_path",       &prad2::Pipeline::runinfo_path)
        .def_readonly("hycal_map_path",     &prad2::Pipeline::hycal_map_path)
        .def_readonly("gem_map_path",       &prad2::Pipeline::gem_map_path)
        .def_readonly("hycal_calib_path",     &prad2::Pipeline::hycal_calib_path)
        .def_readonly("hycal_time_cut_path",  &prad2::Pipeline::hycal_time_cut_path)
        .def_readonly("gem_pedestal_path",    &prad2::Pipeline::gem_pedestal_path)
        .def_readonly("gem_common_mode_path", &prad2::Pipeline::gem_common_mode_path);

    // PipelineBuilder — fluent setters return *this so calls chain.
    // pybind11 needs reference_internal to keep the chain working from
    // Python (otherwise each setter would copy a fresh builder).
    py::class_<prad2::PipelineBuilder>(m, "PipelineBuilder",
        "Fluent builder.  Empty path strings fall back to defaults; relative "
        "paths resolve via PRAD2_DATABASE_DIR (or set_database_dir override). "
        "build() throws on missing daq_config; missing runinfo / hycal_map / "
        "gem_map / calibration / pedestals warn and proceed.")
        .def(py::init<>())
        .def("set_database_dir",     &prad2::PipelineBuilder::set_database_dir,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_daq_config",       &prad2::PipelineBuilder::set_daq_config,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_recon_config",     &prad2::PipelineBuilder::set_recon_config,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_runinfo",          &prad2::PipelineBuilder::set_runinfo,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_hycal_map",        &prad2::PipelineBuilder::set_hycal_map,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_gem_map",          &prad2::PipelineBuilder::set_gem_map,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_hycal_calib",      &prad2::PipelineBuilder::set_hycal_calib,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_hycal_time_cut",   &prad2::PipelineBuilder::set_hycal_time_cut,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_gem_pedestal",     &prad2::PipelineBuilder::set_gem_pedestal,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_gem_common_mode",  &prad2::PipelineBuilder::set_gem_common_mode,
             py::arg("path"), py::return_value_policy::reference_internal)
        .def("set_run_number",       &prad2::PipelineBuilder::set_run_number,
             py::arg("n"), py::return_value_policy::reference_internal)
        .def("set_run_number_from_evio",
             &prad2::PipelineBuilder::set_run_number_from_evio,
             py::arg("evio_path"), py::return_value_policy::reference_internal)
        .def("set_log_pedestal_checksum",
             &prad2::PipelineBuilder::set_log_pedestal_checksum,
             py::arg("enabled"), py::return_value_policy::reference_internal)
        .def("build",
            [](prad2::PipelineBuilder &self) {
                py::gil_scoped_release rel;
                return self.build();
            },
            "Run all the wiring; returns a Pipeline.  Releases the GIL "
            "during the (potentially slow) Init/LoadCalibration/LoadPedestals "
            "calls so Python threads can keep running.");
}

// -------------------------------------------------------------------------
// Submodule entry point — called from prad2py.cpp
// -------------------------------------------------------------------------
void register_det(py::module_ &m)
{
    auto det = m.def_submodule("det",
        "prad2det bindings — GEM + HyCal reconstruction and slow-control "
        "helpers.");

    bind_gem(det);       // 2a
    bind_hycal(det);     // 2b
    bind_transform(det); // 2c
    bind_pipeline(det);  // 2d — one-stop wiring (PipelineBuilder)
}
