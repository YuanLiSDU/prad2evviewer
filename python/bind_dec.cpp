// bind_dec.cpp — pybind11 bindings for prad2dec (prad2py.dec submodule)
//
// Exposes:
//   prad2py.dec.EventType            (enum)
//   prad2py.dec.Status               (enum)
//   prad2py.dec.DaqConfig            (config struct)
//   prad2py.dec.EventInfo            (per-event metadata)
//   prad2py.dec.ChannelData / SlotData / RocData / EventData  (fdec)
//   prad2py.dec.ApvAddress / ApvData / MpdData / SspEventData (ssp)
//   prad2py.dec.EcPeak / EcCluster / PradCluster / VtpBlock / VtpEventData  (vtp)
//   prad2py.dec.TdcHit / TdcEventData                         (tdc)
//   prad2py.dec.SyncInfo                                      (sync)
//   prad2py.dec.EpicsStore / EpicsSnapshot                    (slow control)
//   prad2py.dec.EvChannel            (evio reader)
//   prad2py.dec.load_daq_config(path) -> DaqConfig
//
// The bulk per-channel arrays (ChannelData.samples, ApvData.strips, etc.)
// are returned as numpy arrays — copies by default so the buffer stays
// valid after the next DecodeEvent call.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"
#include "Fadc250FwAnalyzer.h"
#include "SspData.h"
#include "VtpData.h"
#include "TdcData.h"
#include "TdcDecoder.h"
#include "SyncData.h"
#include "EpicsStore.h"
#include "HVDecoder.h"

#include <cstdlib>
#include <memory>
#include <string>

#include "InstallPaths.h"

namespace py = pybind11;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace {

std::string default_daq_config_path()
{
    // Same lookup policy as prad2py.cpp's copy (these two TUs each have
    // a `load_daq_config("")` overload that needs the default path).
    std::string dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../../share/prad2evviewer/database",
         "../share/prad2evviewer/database"},
        DATABASE_DIR);
    return dir + "/daq_config.json";
}

// -------------------------------------------------------------------------
// FADC250 data types (fdec::)
// -------------------------------------------------------------------------
void bind_fadc(py::module_ &m)
{
    // ChannelData — bulk samples are exposed as a numpy view over the
    // first `nsamples` elements.  The buffer belongs to the owning event,
    // so callers must copy (default) to keep data across DecodeEvent calls.
    py::class_<fdec::ChannelData>(m, "ChannelData",
        "Per-channel FADC250 samples (nsamples uint16 values).")
        .def_readonly("nsamples", &fdec::ChannelData::nsamples)
        .def_property_readonly("samples",
            [](const fdec::ChannelData &c) {
                // Always copy — the owning EventData may be reused/overwritten.
                // Explicit copy: allocate a fresh numpy array with no base
                // pointer, then std::copy_n into it.  The (shape,strides,ptr)
                // constructor without a base is interpreted by pybind11 as
                // "wrap this pointer", which silently corrupts memory when the
                // C++ object is freed (we hit this exact bug — accessing
                // .samples poisoned unrelated SspEventData buffers stored
                // alongside).  Allocating-then-copying is unambiguously safe.
                py::array_t<uint16_t> arr(c.nsamples);
                std::copy_n(c.samples, c.nsamples, arr.mutable_data());
                return arr;
            },
            "16-bit ADC samples as a fresh numpy array (copy of nsamples values).");

    py::class_<fdec::SlotData>(m, "SlotData",
        "One FADC250 slot: trigger info plus per-channel samples.")
        .def_readonly("present",      &fdec::SlotData::present)
        .def_readonly("trigger",      &fdec::SlotData::trigger)
        .def_readonly("timestamp",    &fdec::SlotData::timestamp)
        .def_readonly("nchannels",    &fdec::SlotData::nchannels)
        .def_readonly("channel_mask", &fdec::SlotData::channel_mask)
        .def("channel",
            [](const fdec::SlotData &s, int ch) -> const fdec::ChannelData& {
                if (ch < 0 || ch >= fdec::MAX_CHANNELS)
                    throw py::index_error("channel out of range");
                return s.channels[ch];
            },
            py::arg("channel"),
            py::return_value_policy::reference_internal,
            "Access ChannelData by channel index (no presence check — use channel_mask).")
        .def("present_channels",
            [](const fdec::SlotData &s) {
                std::vector<int> out;
                out.reserve(s.nchannels);
                for (int c = 0; c < fdec::MAX_CHANNELS; ++c)
                    if (s.channel_mask & (1ULL << c)) out.push_back(c);
                return out;
            },
            "List of channel indices with nsamples > 0 this event.");

    py::class_<fdec::RocData>(m, "RocData",
        "One ROC crate worth of FADC data.")
        .def_readonly("present", &fdec::RocData::present)
        .def_readonly("tag",     &fdec::RocData::tag)
        .def_readonly("nslots",  &fdec::RocData::nslots)
        .def("slot",
            [](const fdec::RocData &r, int sl) -> const fdec::SlotData& {
                if (sl < 0 || sl >= fdec::MAX_SLOTS)
                    throw py::index_error("slot out of range");
                return r.slots[sl];
            },
            py::arg("slot"),
            py::return_value_policy::reference_internal)
        .def("present_slots",
            [](const fdec::RocData &r) {
                std::vector<int> out;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s)
                    if (r.slots[s].present) out.push_back(s);
                return out;
            });

    py::class_<fdec::EventInfo>(m, "EventInfo",
        "Event-level metadata (TI + trigger bank).")
        .def_readwrite("type",            &fdec::EventInfo::type)
        .def_readwrite("trigger_type",    &fdec::EventInfo::trigger_type)
        .def_readwrite("trigger_bits",    &fdec::EventInfo::trigger_bits)
        .def_readwrite("event_tag",       &fdec::EventInfo::event_tag)
        .def_readwrite("event_number",    &fdec::EventInfo::event_number)
        .def_readwrite("trigger_number",  &fdec::EventInfo::trigger_number)
        .def_readwrite("timestamp",       &fdec::EventInfo::timestamp)
        .def_readwrite("run_number",      &fdec::EventInfo::run_number)
        .def_readwrite("unix_time",       &fdec::EventInfo::unix_time)
        .def("__repr__", [](const fdec::EventInfo &i) {
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                "<EventInfo evt=%d trig=%d tag=0x%x bits=0x%x ts=%llu>",
                i.event_number, i.trigger_number, i.event_tag, i.trigger_bits,
                (unsigned long long)i.timestamp);
            return std::string(buf);
        });

    py::class_<psync::SyncInfo>(m, "SyncInfo",
        "Absolute-time / run-state snapshot populated from SYNC/EPICS 0xE112 "
        "HEAD banks and from PRESTART/GO/END control-event payloads.  Persists "
        "across events in EvChannel — see EvChannel.sync().")
        .def_readwrite("run_number",   &psync::SyncInfo::run_number)
        .def_readwrite("sync_counter", &psync::SyncInfo::sync_counter)
        .def_readwrite("unix_time",    &psync::SyncInfo::unix_time)
        .def_readwrite("event_tag",    &psync::SyncInfo::event_tag)
        .def_readwrite("run_type",     &psync::SyncInfo::run_type)
        .def("valid", &psync::SyncInfo::valid,
             "True if unix_time has been populated (any SYNC/control event seen).")
        .def("__repr__", [](const psync::SyncInfo &s) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<SyncInfo run=%u counter=%u unix=%u tag=0x%x type=%u>",
                s.run_number, s.sync_counter, s.unix_time, s.event_tag,
                (unsigned)s.run_type);
            return std::string(buf);
        });

    py::class_<fdec::EventData, std::shared_ptr<fdec::EventData>>(m, "EventData",
        "Full decoded event: EventInfo + FADC ROCs.")
        .def(py::init<>())
        .def_readonly("info",   &fdec::EventData::info)
        .def_readonly("nrocs",  &fdec::EventData::nrocs)
        .def("roc",
            [](const fdec::EventData &e, int i) -> const fdec::RocData& {
                if (i < 0 || i >= e.nrocs)
                    throw py::index_error("ROC index out of range");
                return e.rocs[e.roc_index[i]];
            },
            py::arg("index"),
            py::return_value_policy::reference_internal,
            "Access the i-th active ROC by slot in roc_index[].")
        .def("find_roc",
            [](const fdec::EventData &e, uint32_t tag) -> const fdec::RocData* {
                return e.findRoc(tag);   // nullptr → None via pybind11
            },
            py::arg("tag"),
            py::return_value_policy::reference_internal,
            "Locate a ROC by its bank tag (e.g. 0x80). Returns None if absent.")
        .def("clear", &fdec::EventData::clear);

    // ----- WaveAnalyzer types --------------------------------------------
    py::class_<fdec::Peak>(m, "Peak",
        "One FADC pulse found by WaveAnalyzer.\n\n"
        "  left / right — INCLUSIVE integration bounds (both samples are in `integral`)\n"
        "  pos          — raw-sample maximum near the smoothed peak\n"
        "  quality      — Q_PEAK_* bitmask (currently just Q_PEAK_PILED)")
        .def(py::init<>())
        .def_readwrite("height",   &fdec::Peak::height)
        .def_readwrite("integral", &fdec::Peak::integral)
        .def_readwrite("time",     &fdec::Peak::time)
        .def_readwrite("pos",      &fdec::Peak::pos)
        .def_readwrite("left",     &fdec::Peak::left)
        .def_readwrite("right",    &fdec::Peak::right)
        .def_readwrite("overflow", &fdec::Peak::overflow)
        .def_readwrite("quality",  &fdec::Peak::quality);

    py::class_<fdec::Pedestal>(m, "Pedestal",
        "Per-channel pedestal estimate from WaveAnalyzer.\n\n"
        "Fields:\n"
        "  mean / rms — final values after iterative outlier rejection\n"
        "  nused      — samples that survived rejection (≤ ped_nsamples)\n"
        "  quality    — Q_PED_* bitmask (see module-level Q_PED_* attrs)\n"
        "  slope      — least-squares slope across the survivors (ADC/sample)")
        .def(py::init<>())
        .def_readwrite("mean",    &fdec::Pedestal::mean)
        .def_readwrite("rms",     &fdec::Pedestal::rms)
        .def_readwrite("nused",   &fdec::Pedestal::nused)
        .def_readwrite("quality", &fdec::Pedestal::quality)
        .def_readwrite("slope",   &fdec::Pedestal::slope);

    py::class_<fdec::WaveConfig::NnlsDeconvConfig>(m, "NnlsDeconvConfig",
        "Nested config (fdec::WaveConfig::NnlsDeconvConfig) — knobs for the "
        "NNLS pile-up deconvolution.  Mirrors the JSON layout in "
        "database/daq_config.json under fadc250_waveform.analyzer.nnls_deconv.")
        .def(py::init<>())
        .def_readwrite("enabled",            &fdec::WaveConfig::NnlsDeconvConfig::enabled)
        .def_readwrite("template_file",      &fdec::WaveConfig::NnlsDeconvConfig::template_file)
        .def_readwrite("apply_to_all_peaks", &fdec::WaveConfig::NnlsDeconvConfig::apply_to_all_peaks)
        .def_readwrite("tau_r_min_ns",    &fdec::WaveConfig::NnlsDeconvConfig::tau_r_min_ns)
        .def_readwrite("tau_r_max_ns",    &fdec::WaveConfig::NnlsDeconvConfig::tau_r_max_ns)
        .def_readwrite("tau_f_min_ns",    &fdec::WaveConfig::NnlsDeconvConfig::tau_f_min_ns)
        .def_readwrite("tau_f_max_ns",    &fdec::WaveConfig::NnlsDeconvConfig::tau_f_max_ns)
        .def_readwrite("shape_window_factor", &fdec::WaveConfig::NnlsDeconvConfig::shape_window_factor)
        .def_readwrite("t0_window_ns",        &fdec::WaveConfig::NnlsDeconvConfig::t0_window_ns)
        .def_readwrite("amp_max_factor",      &fdec::WaveConfig::NnlsDeconvConfig::amp_max_factor)
        .def_readwrite("pre_samples",         &fdec::WaveConfig::NnlsDeconvConfig::pre_samples)
        .def_readwrite("post_samples",        &fdec::WaveConfig::NnlsDeconvConfig::post_samples);

    py::class_<fdec::PulseTemplateStore>(m, "PulseTemplateStore",
        "Per-type pulse-template store loaded from the JSON written by "
        "fit_pulse_template.py.  Bind via WaveAnalyzer.set_template_store() "
        "and the analyzer's Analyze() method will auto-deconvolve piled "
        "events using the per-type template (PbGlass / PbWO4 / LMS / Veto) "
        "for each channel's category.")
        .def(py::init<>())
        .def("load_from_file", &fdec::PulseTemplateStore::LoadFromFile,
             py::arg("path"), py::arg("wave_cfg"),
             "Load per-type templates from the JSON's `_by_type` block "
             "and the (roc_tag, slot, channel) → module_type lookup from "
             "the per-channel records.  wave_cfg's `nnls_deconv` τ-range "
             "gates are applied to each per-type entry; `clk_mhz` sets "
             "the precomputed grid period.  Returns False on file-not-"
             "found / parse error / empty contents (caller falls back to "
             "non-deconv mode).")
        .def("lookup",
             [](const fdec::PulseTemplateStore &self,
                int roc_tag, int slot, int channel) -> py::object {
                 const fdec::PulseTemplate *t = self.Lookup(
                     roc_tag, slot, channel);
                 if (!t) return py::none();
                 return py::cast(*t);
             },
             py::arg("roc_tag"), py::arg("slot"), py::arg("channel"),
             "Resolve the channel's module type and return the matching "
             "per-type template. Returns None when the channel is not in "
             "the loaded JSON or its type lacks a usable per-type entry.")
        .def("clear", &fdec::PulseTemplateStore::Clear)
        .def_property_readonly("valid", &fdec::PulseTemplateStore::valid)
        .def_property_readonly("n_channels_known",
            &fdec::PulseTemplateStore::n_channels_known)
        .def_property_readonly("n_types_loaded",
            &fdec::PulseTemplateStore::n_types_loaded)
        .def("type_template",
            [](const fdec::PulseTemplateStore &self,
               const std::string &type_name) -> py::object {
                const fdec::PulseTemplate *t = self.type_template(type_name);
                if (!t) return py::none();
                return py::cast(*t);
            },
            py::arg("type_name"),
            "Per-type template (`type_name` ∈ "
            "{'PbGlass', 'PbWO4', 'LMS', 'Veto'}) parsed from the JSON's "
            "`_by_type` block.  Returns None if the type isn't present.");

    py::class_<fdec::WaveConfig>(m, "WaveConfig",
        "Knobs for WaveAnalyzer (smoothing, thresholds, pedestal window, ...).")
        .def(py::init<>())
        .def(py::init<const fdec::WaveConfig &>(), py::arg("other"),
             "Copy-construct from another WaveConfig (use this when "
             "snapshotting cfg.wave_cfg from a DaqConfig — direct "
             "assignment would alias the underlying field).")
        .def_readwrite("smooth_order",     &fdec::WaveConfig::smooth_order)
        .def_readwrite("peak_nsigma",      &fdec::WaveConfig::peak_nsigma)
        .def_readwrite("min_peak_height",  &fdec::WaveConfig::min_peak_height)
        .def_readwrite("min_peak_ratio",   &fdec::WaveConfig::min_peak_ratio)
        .def_readwrite("int_tail_ratio",   &fdec::WaveConfig::int_tail_ratio)
        .def_readwrite("tail_break_n",     &fdec::WaveConfig::tail_break_n)
        .def_readwrite("peak_pileup_gap",  &fdec::WaveConfig::peak_pileup_gap)
        .def_readwrite("ped_nsamples",     &fdec::WaveConfig::ped_nsamples)
        .def_readwrite("ped_flatness",   &fdec::WaveConfig::ped_flatness)
        .def_readwrite("ped_max_iter",   &fdec::WaveConfig::ped_max_iter)
        .def_readwrite("overflow",       &fdec::WaveConfig::overflow)
        .def_readwrite("clk_mhz",        &fdec::WaveConfig::clk_mhz)
        .def_readwrite("nnls_deconv",    &fdec::WaveConfig::nnls_deconv);

    py::class_<fdec::PulseTemplate>(m, "PulseTemplate",
        "Per-channel two-tau pulse template used by WaveAnalyzer.deconvolve(). "
        "tau_r_ns / tau_f_ns are the median values from "
        "fit_pulse_template.py.  is_global=True marks a per-type aggregate "
        "(the only kind currently produced by PulseTemplateStore).")
        .def(py::init<>())
        .def(py::init([](float tr, float tf, bool is_global) {
                fdec::PulseTemplate t{tr, tf, is_global};
                return t;
             }),
             py::arg("tau_r_ns"), py::arg("tau_f_ns"),
             py::arg("is_global") = false)
        .def_readwrite("tau_r_ns",  &fdec::PulseTemplate::tau_r_ns)
        .def_readwrite("tau_f_ns",  &fdec::PulseTemplate::tau_f_ns)
        .def_readwrite("is_global", &fdec::PulseTemplate::is_global);

    py::class_<fdec::WaveAnalyzer::PulseFitResult>(m, "PulseFitResult",
        "Output of WaveAnalyzer.fit_pulse_shape() — three-parameter "
        "Levenberg-Marquardt fit of the unit-amplitude two-tau model to "
        "a single normalised pulse.  Inspect `.ok` before reading the "
        "shape parameters.")
        .def_readonly("ok",            &fdec::WaveAnalyzer::PulseFitResult::ok)
        .def_readonly("t0_ns",         &fdec::WaveAnalyzer::PulseFitResult::t0_ns)
        .def_readonly("tau_r_ns",      &fdec::WaveAnalyzer::PulseFitResult::tau_r_ns)
        .def_readonly("tau_f_ns",      &fdec::WaveAnalyzer::PulseFitResult::tau_f_ns)
        .def_readonly("peak_amp",      &fdec::WaveAnalyzer::PulseFitResult::peak_amp)
        .def_readonly("chi2_per_dof",  &fdec::WaveAnalyzer::PulseFitResult::chi2_per_dof)
        .def_readonly("n_iter",        &fdec::WaveAnalyzer::PulseFitResult::n_iter);

    py::class_<fdec::WaveAnalyzer::PulseFitTwoTauPResult>(m, "PulseFitTwoTauPResult",
        "Output of WaveAnalyzer.fit_pulse_shape_two_tau_p() — same as "
        "PulseFitResult plus the rise-edge exponent `p` "
        "(T = [1-exp(-(t-t0)/τ_r)]^p · exp(-(t-t0)/τ_f)).")
        .def_readonly("ok",            &fdec::WaveAnalyzer::PulseFitTwoTauPResult::ok)
        .def_readonly("t0_ns",         &fdec::WaveAnalyzer::PulseFitTwoTauPResult::t0_ns)
        .def_readonly("tau_r_ns",      &fdec::WaveAnalyzer::PulseFitTwoTauPResult::tau_r_ns)
        .def_readonly("tau_f_ns",      &fdec::WaveAnalyzer::PulseFitTwoTauPResult::tau_f_ns)
        .def_readonly("p",             &fdec::WaveAnalyzer::PulseFitTwoTauPResult::p)
        .def_readonly("peak_amp",      &fdec::WaveAnalyzer::PulseFitTwoTauPResult::peak_amp)
        .def_readonly("chi2_per_dof",  &fdec::WaveAnalyzer::PulseFitTwoTauPResult::chi2_per_dof)
        .def_readonly("n_iter",        &fdec::WaveAnalyzer::PulseFitTwoTauPResult::n_iter);

    py::class_<fdec::WaveResult>(m, "WaveResult",
        "Output of WaveAnalyzer.analyze_result(): {ped, npeaks, peaks}. "
        "Use this as the second argument to WaveAnalyzer.deconvolve() to "
        "avoid re-running the peak finder.")
        .def(py::init<>())
        .def_readwrite("ped",    &fdec::WaveResult::ped)
        .def_readonly("npeaks",  &fdec::WaveResult::npeaks)
        .def_property_readonly("peaks", [](const fdec::WaveResult &self) {
            py::list out;
            for (int i = 0; i < self.npeaks; ++i) out.append(self.peaks[i]);
            return out;
        });

    py::class_<fdec::DeconvOutput>(m, "DeconvOutput",
        "NNLS pile-up deconvolution output.  amplitude/height/integral are "
        "Python lists of length n (= WaveResult.npeaks).  state is a "
        "single-valued enum (Q_DECONV_*); on any non-converged outcome the "
        "lists are empty.")
        .def(py::init<>())
        .def_readonly("state",        &fdec::DeconvOutput::state)
        .def_readonly("n",            &fdec::DeconvOutput::n)
        .def_readonly("chi2_per_dof", &fdec::DeconvOutput::chi2_per_dof)
        .def_property_readonly("amplitude", [](const fdec::DeconvOutput &self) {
            return std::vector<float>(self.amplitude, self.amplitude + self.n);
        })
        .def_property_readonly("height", [](const fdec::DeconvOutput &self) {
            return std::vector<float>(self.height, self.height + self.n);
        })
        .def_property_readonly("integral", [](const fdec::DeconvOutput &self) {
            return std::vector<float>(self.integral, self.integral + self.n);
        })
        .def_property_readonly("t0_ns", [](const fdec::DeconvOutput &self) {
            return std::vector<float>(self.t0_ns, self.t0_ns + self.n);
        })
        .def_property_readonly("tau_r_ns", [](const fdec::DeconvOutput &self) {
            return std::vector<float>(self.tau_r_ns, self.tau_r_ns + self.n);
        })
        .def_property_readonly("tau_f_ns", [](const fdec::DeconvOutput &self) {
            return std::vector<float>(self.tau_f_ns, self.tau_f_ns + self.n);
        });

    py::class_<fdec::WaveAnalyzer>(m, "WaveAnalyzer",
        "Fast C++ pedestal + peak finder used by the server.  Call "
        "analyze(samples_numpy) to get (ped_mean, ped_rms, [Peak, ...]).")
        .def(py::init<const fdec::WaveConfig &>(),
             py::arg("cfg") = fdec::WaveConfig{})
        .def_readwrite("cfg", &fdec::WaveAnalyzer::cfg)
        .def("analyze",
            [](const fdec::WaveAnalyzer &self, py::array_t<uint16_t> samples) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::WaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), res);
                }
                py::list peaks;
                for (int i = 0; i < res.npeaks; ++i)
                    peaks.append(res.peaks[i]);
                return py::make_tuple(res.ped.mean, res.ped.rms, peaks);
            },
            py::arg("samples"),
            "Analyze one channel's FADC waveform (uint16 numpy array). "
            "Returns (pedestal_mean, pedestal_rms, peaks_list).  Peaks "
            "above ``max(cfg.peak_nsigma × ped.rms, cfg.min_peak_height)`` "
            "are kept; up to ``MAX_PEAKS`` per call.")
        .def("analyze_full",
            [](const fdec::WaveAnalyzer &self, py::array_t<uint16_t> samples) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::WaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), res);
                }
                py::list peaks;
                for (int i = 0; i < res.npeaks; ++i)
                    peaks.append(res.peaks[i]);
                return py::make_tuple(res.ped, peaks);
            },
            py::arg("samples"),
            "Like analyze() but returns the full (Pedestal, [Peak, ...]) "
            "tuple — Pedestal exposes nused / quality / slope, each Peak "
            "exposes its quality bitmask.  Use this when documenting the "
            "analyzer or when downstream code needs the quality flags.")
        .def("smooth",
            [](const fdec::WaveAnalyzer &self, py::array_t<uint16_t> samples) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                const int n = static_cast<int>(buf.shape[0]);
                if (n > fdec::MAX_SAMPLES)
                    throw py::value_error("samples length exceeds MAX_SAMPLES");
                py::array_t<float> out(n);
                {
                    py::gil_scoped_release rel;
                    self.smooth(static_cast<const uint16_t*>(buf.ptr), n,
                                static_cast<float*>(out.request().ptr));
                }
                return out;
            },
            py::arg("samples"),
            "Run only the triangular-kernel smoothing pass (no pedestal / "
            "peak finding).  Returns the smoothed buffer as a float32 "
            "numpy array — useful for plotting the curve the peak finder "
            "actually sees.")
        .def_static("fit_pulse_shape",
            [](py::array_t<uint16_t> slice, int peak_idx,
               float ped, float ped_rms, float clk_ns,
               float model_err_floor) {
                py::buffer_info buf = slice.request();
                if (buf.ndim != 1)
                    throw py::value_error("slice must be a 1-D uint16 array");
                fdec::WaveAnalyzer::PulseFitResult res;
                {
                    py::gil_scoped_release rel;
                    res = fdec::WaveAnalyzer::FitPulseShape(
                        static_cast<const uint16_t*>(buf.ptr),
                        static_cast<int>(buf.shape[0]),
                        peak_idx, ped, ped_rms, clk_ns,
                        model_err_floor);
                }
                return res;
            },
            py::arg("slice"), py::arg("peak_idx"),
            py::arg("ped"), py::arg("ped_rms"), py::arg("clk_ns"),
            py::arg("model_err_floor") = 0.01f,
            "Static C++ Levenberg-Marquardt fit of the unit-amplitude "
            "two-tau pulse-shape model to one normalised pulse — orders "
            "of magnitude faster than scipy.optimize.curve_fit on the "
            "equivalent model.  `slice` is a uint16 waveform window "
            "around the peak; `peak_idx` is the position of the maximum "
            "within the slice; `ped` / `ped_rms` come from the same "
            "WaveAnalyzer pass that produced the peak.  Returns a "
            "PulseFitResult; check `.ok` before reading the params.")
        .def_static("fit_pulse_shape_two_tau_p",
            [](py::array_t<uint16_t> slice, int peak_idx,
               float ped, float ped_rms, float clk_ns,
               float model_err_floor) {
                py::buffer_info buf = slice.request();
                if (buf.ndim != 1)
                    throw py::value_error("slice must be a 1-D uint16 array");
                fdec::WaveAnalyzer::PulseFitTwoTauPResult res;
                {
                    py::gil_scoped_release rel;
                    res = fdec::WaveAnalyzer::FitPulseShapeTwoTauP(
                        static_cast<const uint16_t*>(buf.ptr),
                        static_cast<int>(buf.shape[0]),
                        peak_idx, ped, ped_rms, clk_ns,
                        model_err_floor);
                }
                return res;
            },
            py::arg("slice"), py::arg("peak_idx"),
            py::arg("ped"), py::arg("ped_rms"), py::arg("clk_ns"),
            py::arg("model_err_floor") = 0.01f,
            "Same machinery as fit_pulse_shape() but with the four-"
            "parameter two-tau-with-rise-exponent model: "
            "T = [1-exp(-(t-t0)/τ_r)]^p · exp(-(t-t0)/τ_f).  Allows the "
            "rise *shape* to vary independently of the rise *timescale* — "
            "addresses the systematic onset mismatch the standard two-tau "
            "form has against multi-stage-shaped PMT pulses.")
        .def("set_template_store",
            [](fdec::WaveAnalyzer &self,
               const fdec::PulseTemplateStore *store) {
                self.SetTemplateStore(store);
            },
            py::arg("store").none(true),
            py::keep_alive<1, 2>(),
            "Bind a PulseTemplateStore so subsequent analyze() calls can "
            "auto-deconvolve.  Pass None to detach.  The analyzer keeps "
            "a non-owning reference; the store must outlive the analyzer.")
        .def("set_channel_key", &fdec::WaveAnalyzer::SetChannelKey,
            py::arg("roc_tag"), py::arg("slot"), py::arg("channel"),
            "Tell the analyzer which channel the next analyze() call is "
            "for, so the bound PulseTemplateStore can look up the right "
            "template.  Pass any negative value to disable deconv.")
        .def("clear_channel_key", &fdec::WaveAnalyzer::ClearChannelKey)
        .def("analyze_result",
            [](const fdec::WaveAnalyzer &self, py::array_t<uint16_t> samples) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::WaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), res);
                }
                return res;
            },
            py::arg("samples"),
            "Like analyze() but returns the full WaveResult object — pass "
            "this back to deconvolve() so the deconvolver sees the same "
            "pedestal and peak set the analyzer found.")
        .def("deconvolve",
            [](const fdec::WaveAnalyzer &self,
               py::array_t<uint16_t> samples,
               const fdec::WaveResult &wres,
               const fdec::PulseTemplate &tmpl) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::DeconvOutput out;
                {
                    py::gil_scoped_release rel;
                    self.Deconvolve(static_cast<const uint16_t*>(buf.ptr),
                                    static_cast<int>(buf.shape[0]),
                                    wres, tmpl, out);
                }
                return out;
            },
            py::arg("samples"), py::arg("wres"), py::arg("template_"),
            "Run pile-up deconvolution against the supplied template "
            "`template_` (typically a per-type entry from "
            "PulseTemplateStore.lookup() / type_template()).  Caller must "
            "have already obtained `wres` from analyze_result() with the "
            "same samples.  Returns a DeconvOutput; check `.state` "
            "(Q_DECONV_*) before reading the height/integral arrays.");

    // ----- Firmware-faithful (Mode 1/2/3) analyzer -----------------------
    // Quality bitmask constants — exposed at module scope so callers can
    // test pulse.quality & dec.Q_DAQ_PEAK_AT_BOUNDARY etc.
    m.attr("Q_DAQ_GOOD")             = py::int_(fdec::Q_DAQ_GOOD);
    m.attr("Q_DAQ_PEAK_AT_BOUNDARY") = py::int_(fdec::Q_DAQ_PEAK_AT_BOUNDARY);
    m.attr("Q_DAQ_NSB_TRUNCATED")    = py::int_(fdec::Q_DAQ_NSB_TRUNCATED);
    m.attr("Q_DAQ_NSA_TRUNCATED")    = py::int_(fdec::Q_DAQ_NSA_TRUNCATED);
    m.attr("Q_DAQ_VA_OUT_OF_RANGE")  = py::int_(fdec::Q_DAQ_VA_OUT_OF_RANGE);

    // Soft-analyzer peak quality bitmask — currently just the pile-up
    // flag.  Both peaks in a piled-up pair get the bit set.
    m.attr("Q_PEAK_GOOD")            = py::int_(fdec::Q_PEAK_GOOD);
    m.attr("Q_PEAK_PILED")           = py::int_(fdec::Q_PEAK_PILED);
    m.attr("Q_PEAK_DECONVOLVED")     = py::int_(fdec::Q_PEAK_DECONVOLVED);

    // Pedestal-fit quality bitmask — exposed so callers can filter on
    // Pedestal.quality (e.g. dec.Q_PED_PULSE_IN_WINDOW).
    m.attr("Q_PED_GOOD")             = py::int_(fdec::Q_PED_GOOD);
    m.attr("Q_PED_NOT_CONVERGED")    = py::int_(fdec::Q_PED_NOT_CONVERGED);
    m.attr("Q_PED_FLOOR_ACTIVE")     = py::int_(fdec::Q_PED_FLOOR_ACTIVE);
    m.attr("Q_PED_TOO_FEW_SAMPLES")  = py::int_(fdec::Q_PED_TOO_FEW_SAMPLES);
    m.attr("Q_PED_PULSE_IN_WINDOW")  = py::int_(fdec::Q_PED_PULSE_IN_WINDOW);
    m.attr("Q_PED_OVERFLOW")         = py::int_(fdec::Q_PED_OVERFLOW);
    m.attr("Q_PED_TRAILING_WINDOW")  = py::int_(fdec::Q_PED_TRAILING_WINDOW);

    // NNLS deconv state — single-valued enum (not a bitmask).  Test
    // dec_out.state == dec.Q_DECONV_APPLIED (or _FALLBACK_GLOBAL) to know
    // whether to read dec_out.height / .integral.
    m.attr("Q_DECONV_NOT_RUN")          = py::int_(fdec::Q_DECONV_NOT_RUN);
    m.attr("Q_DECONV_NO_TEMPLATE")      = py::int_(fdec::Q_DECONV_NO_TEMPLATE);
    m.attr("Q_DECONV_BAD_TEMPLATE")     = py::int_(fdec::Q_DECONV_BAD_TEMPLATE);
    m.attr("Q_DECONV_LM_NOT_CONVERGED") = py::int_(fdec::Q_DECONV_LM_NOT_CONVERGED);
    m.attr("Q_DECONV_APPLIED")          = py::int_(fdec::Q_DECONV_APPLIED);
    m.attr("Q_DECONV_FALLBACK_GLOBAL")  = py::int_(fdec::Q_DECONV_FALLBACK_GLOBAL);

    py::class_<fdec::DaqPeak>(m, "DaqPeak",
        "One firmware-mode pulse (Mode 1 + Mode 2 + Mode 3 result).")
        .def(py::init<>())
        .def_readwrite("pulse_id",     &fdec::DaqPeak::pulse_id)
        .def_readwrite("vmin",         &fdec::DaqPeak::vmin)
        .def_readwrite("vpeak",        &fdec::DaqPeak::vpeak)
        .def_readwrite("va",           &fdec::DaqPeak::va)
        .def_readwrite("coarse",       &fdec::DaqPeak::coarse)
        .def_readwrite("fine",         &fdec::DaqPeak::fine)
        .def_readwrite("time_units",   &fdec::DaqPeak::time_units)
        .def_readwrite("time_ns",      &fdec::DaqPeak::time_ns)
        .def_readwrite("cross_sample", &fdec::DaqPeak::cross_sample)
        .def_readwrite("peak_sample",  &fdec::DaqPeak::peak_sample)
        .def_readwrite("integral",     &fdec::DaqPeak::integral)
        .def_readwrite("window_lo",    &fdec::DaqPeak::window_lo)
        .def_readwrite("window_hi",    &fdec::DaqPeak::window_hi)
        .def_readwrite("quality",      &fdec::DaqPeak::quality);

    py::class_<evc::DaqConfig::Fadc250FwConfig>(m, "Fadc250FwConfig",
        "Firmware Mode 1/2/3 emulation knobs.  Names follow the Hall-D V3 "
        "firmware API (faV3HallDSetProcMode).  NPEAK is an alias for "
        "MAX_PULSES.")
        .def(py::init<>())
        .def_readwrite("TET",        &evc::DaqConfig::Fadc250FwConfig::TET)
        .def_readwrite("NSB",        &evc::DaqConfig::Fadc250FwConfig::NSB)
        .def_readwrite("NSA",        &evc::DaqConfig::Fadc250FwConfig::NSA)
        .def_readwrite("MAX_PULSES", &evc::DaqConfig::Fadc250FwConfig::MAX_PULSES)
        .def_readwrite("NSAT",       &evc::DaqConfig::Fadc250FwConfig::NSAT)
        .def_readwrite("NPED",       &evc::DaqConfig::Fadc250FwConfig::NPED)
        .def_readwrite("MAXPED",     &evc::DaqConfig::Fadc250FwConfig::MAXPED)
        .def_readwrite("CLK_NS",     &evc::DaqConfig::Fadc250FwConfig::CLK_NS);

    py::class_<fdec::Fadc250FwAnalyzer>(m, "Fadc250FwAnalyzer",
        "Firmware-faithful FADC250 Mode 1/2/3 emulator.  Mirrors the on-board "
        "TDC + pulse-windowing exactly (per FADC250 User's Manual).  Use this "
        "to compare against firmware-reported pulse data, or for DAQ signal "
        "studies where the soft (local-maxima) analyzer would diverge.")
        .def(py::init<const evc::DaqConfig::Fadc250FwConfig &>(),
             py::arg("cfg") = evc::DaqConfig::Fadc250FwConfig{})
        .def_readwrite("cfg", &fdec::Fadc250FwAnalyzer::cfg)
        .def("analyze",
            [](const fdec::Fadc250FwAnalyzer &self,
               py::array_t<uint16_t> samples, float ped) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::DaqWaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), ped, res);
                }
                py::list peaks;
                for (int i = 0; i < res.npeaks; ++i)
                    peaks.append(res.peaks[i]);
                return py::make_tuple(res.vnoise, peaks);
            },
            py::arg("samples"), py::arg("ped"),
            "Run firmware Mode 3 (TDC) → Mode 1 + Mode 2 windowing on a "
            "uint16 sample array.  Returns (vnoise, [DaqPeak, ...]).  ``ped`` "
            "is the per-channel pedestal (firmware register or soft-analyzer "
            "estimate); samples have not been pedestal-subtracted.")
        .def("analyze_full",
            [](const fdec::Fadc250FwAnalyzer &self,
               py::array_t<uint16_t> samples, float ped) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::DaqWaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), ped, res);
                }
                py::list peaks;
                for (int i = 0; i < res.npeaks; ++i)
                    peaks.append(res.peaks[i]);
                return py::make_tuple(res.vnoise, peaks);
            },
            py::arg("samples"), py::arg("ped"),
            "Same return shape as analyze() — provided for parity with "
            "WaveAnalyzer.analyze_full().  DaqPeak already exposes every "
            "field including quality, so there's no extra info to surface "
            "here today.");
}

// -------------------------------------------------------------------------
// SSP / MPD / APV (ssp::)
// -------------------------------------------------------------------------
void bind_ssp(py::module_ &m)
{
    py::class_<ssp::ApvAddress>(m, "ApvAddress",
        "APV hardware identifier (crate, MPD/fiber, ADC channel).")
        .def(py::init<>())
        .def_readwrite("crate_id", &ssp::ApvAddress::crate_id)
        .def_readwrite("mpd_id",   &ssp::ApvAddress::mpd_id)
        .def_readwrite("adc_ch",   &ssp::ApvAddress::adc_ch)
        .def("__repr__", [](const ssp::ApvAddress &a) {
            char buf[96];
            std::snprintf(buf, sizeof(buf),
                "<ApvAddress crate=%d mpd=%d adc=%d>",
                a.crate_id, a.mpd_id, a.adc_ch);
            return std::string(buf);
        });

    py::class_<ssp::ApvData>(m, "ApvData",
        "One APV chip's strip/time-sample matrix.")
        .def_readonly("addr",          &ssp::ApvData::addr)
        .def_readonly("present",       &ssp::ApvData::present)
        .def_readonly("nstrips",       &ssp::ApvData::nstrips)
        .def_readonly("flags",         &ssp::ApvData::flags)
        .def_readonly("has_online_cm", &ssp::ApvData::has_online_cm)
        .def_property_readonly("strips",
            [](const ssp::ApvData &a) {
                // [APV_STRIP_SIZE][SSP_TIME_SAMPLES] int16 — explicit
                // allocate-then-copy.  The 3-arg array_t(shape,strides,ptr)
                // constructor has ambiguous ownership and was observed to
                // poison unrelated buffers (see ChannelData.samples note).
                py::array_t<int16_t> arr({
                    static_cast<py::ssize_t>(ssp::APV_STRIP_SIZE),
                    static_cast<py::ssize_t>(ssp::SSP_TIME_SAMPLES)});
                std::copy_n(&a.strips[0][0],
                            ssp::APV_STRIP_SIZE * ssp::SSP_TIME_SAMPLES,
                            arr.mutable_data());
                return arr;
            },
            "(128, 6) int16 numpy array of raw ADC samples (fresh copy).")
        .def_property_readonly("online_cm",
            [](const ssp::ApvData &a) {
                py::array_t<int16_t> arr(static_cast<py::ssize_t>(ssp::SSP_TIME_SAMPLES));
                std::copy_n(a.online_cm, ssp::SSP_TIME_SAMPLES, arr.mutable_data());
                return arr;
            },
            "6-element online common-mode values (fresh copy).")
        .def("has_strip", &ssp::ApvData::hasStrip);

    py::class_<ssp::MpdData>(m, "MpdData",
        "One MPD card: a vector of APVs.")
        .def_readonly("crate_id", &ssp::MpdData::crate_id)
        .def_readonly("mpd_id",   &ssp::MpdData::mpd_id)
        .def_readonly("present",  &ssp::MpdData::present)
        .def_readonly("napvs",    &ssp::MpdData::napvs)
        .def("apv",
            [](const ssp::MpdData &m, int adc) -> const ssp::ApvData& {
                if (adc < 0 || adc >= ssp::MAX_APVS_PER_MPD)
                    throw py::index_error("APV index out of range");
                return m.apvs[adc];
            },
            py::return_value_policy::reference_internal);

    py::class_<ssp::SspEventData, std::shared_ptr<ssp::SspEventData>>(m, "SspEventData",
        "Per-event GEM strip data grouped by (crate, MPD).")
        .def(py::init<>())
        .def_readonly("nmpds", &ssp::SspEventData::nmpds)
        .def("mpd",
            [](const ssp::SspEventData &e, int i) -> const ssp::MpdData& {
                if (i < 0 || i >= e.nmpds)
                    throw py::index_error("MPD index out of range");
                return e.mpds[i];
            },
            py::return_value_policy::reference_internal)
        .def("find_apv",
            [](const ssp::SspEventData &e, int crate, int mpd, int adc)
                -> const ssp::ApvData* {
                // Return a raw pointer and let the .def handler do the
                // cast — manual py::cast(..., reference_internal) loses
                // the parent linkage and keep_alive<0,1> fails with
                // "Could not activate keep_alive".  nullptr → None is
                // handled automatically by pybind11.
                return e.findApv(crate, mpd, adc);
            },
            py::arg("crate"), py::arg("mpd"), py::arg("adc"),
            py::return_value_policy::reference_internal)
        .def("clear", &ssp::SspEventData::clear);
}

// -------------------------------------------------------------------------
// VTP (vtp::)
// -------------------------------------------------------------------------
void bind_vtp(py::module_ &m)
{
    py::class_<vtp::EcPeak>(m, "EcPeak")
        .def_readonly("roc_tag", &vtp::EcPeak::roc_tag)
        .def_readonly("inst",    &vtp::EcPeak::inst)
        .def_readonly("view",    &vtp::EcPeak::view)
        .def_readonly("time",    &vtp::EcPeak::time)
        .def_readonly("coord",   &vtp::EcPeak::coord)
        .def_readonly("energy",  &vtp::EcPeak::energy);

    py::class_<vtp::EcCluster>(m, "EcCluster")
        .def_readonly("roc_tag", &vtp::EcCluster::roc_tag)
        .def_readonly("inst",    &vtp::EcCluster::inst)
        .def_readonly("time",    &vtp::EcCluster::time)
        .def_readonly("energy",  &vtp::EcCluster::energy)
        .def_readonly("coordU",  &vtp::EcCluster::coordU)
        .def_readonly("coordV",  &vtp::EcCluster::coordV)
        .def_readonly("coordW",  &vtp::EcCluster::coordW);

    py::class_<vtp::PradCluster>(m, "PradCluster")
        .def_readonly("roc_tag",  &vtp::PradCluster::roc_tag)
        .def_readonly("energy",   &vtp::PradCluster::energy)
        .def_readonly("id",       &vtp::PradCluster::id)
        .def_readonly("nhits",    &vtp::PradCluster::nhits)
        .def_readonly("time",     &vtp::PradCluster::time)
        .def_property_readonly("is_pbwo4", &vtp::PradCluster::is_pbwo4)
        .def_property_readonly("module",   &vtp::PradCluster::module)
        .def("__repr__", [](const vtp::PradCluster &c) {
            char buf[96];
            snprintf(buf, sizeof(buf),
                     "PradCluster(roc=0x%X, module=%s%u, E=%u, N=%u, T=%u)",
                     c.roc_tag, c.is_pbwo4() ? "W" : "G", c.module(),
                     c.energy, c.nhits, c.time);
            return std::string(buf);
        });

    py::class_<vtp::VtpBlock>(m, "VtpBlock")
        .def_readonly("roc_tag",          &vtp::VtpBlock::roc_tag)
        .def_readonly("slot",             &vtp::VtpBlock::slot)
        .def_readonly("module_id",        &vtp::VtpBlock::module_id)
        .def_readonly("block_number",     &vtp::VtpBlock::block_number)
        .def_readonly("block_level",      &vtp::VtpBlock::block_level)
        .def_readonly("nwords",           &vtp::VtpBlock::nwords)
        .def_readonly("event_number",     &vtp::VtpBlock::event_number)
        .def_readonly("trigger_time",     &vtp::VtpBlock::trigger_time)
        .def_readonly("has_trailer",      &vtp::VtpBlock::has_trailer)
        .def_readonly("trailer_mismatch", &vtp::VtpBlock::trailer_mismatch);

    py::class_<vtp::VtpEventData, std::shared_ptr<vtp::VtpEventData>>(m, "VtpEventData")
        .def(py::init<>())
        .def_readonly("n_peaks",         &vtp::VtpEventData::n_peaks)
        .def_readonly("n_clusters",      &vtp::VtpEventData::n_clusters)
        .def_readonly("n_prad_clusters", &vtp::VtpEventData::n_prad_clusters)
        .def_readonly("n_blocks",        &vtp::VtpEventData::n_blocks)
        .def("peak",
            [](const vtp::VtpEventData &e, int i) -> const vtp::EcPeak& {
                if (i < 0 || i >= e.n_peaks)
                    throw py::index_error("peak index out of range");
                return e.peaks[i];
            },
            py::return_value_policy::reference_internal)
        .def("cluster",
            [](const vtp::VtpEventData &e, int i) -> const vtp::EcCluster& {
                if (i < 0 || i >= e.n_clusters)
                    throw py::index_error("cluster index out of range");
                return e.clusters[i];
            },
            py::return_value_policy::reference_internal)
        .def("prad_cluster",
            [](const vtp::VtpEventData &e, int i) -> const vtp::PradCluster& {
                if (i < 0 || i >= e.n_prad_clusters)
                    throw py::index_error("prad_cluster index out of range");
                return e.prad_clusters[i];
            },
            py::return_value_policy::reference_internal)
        .def("block",
            [](const vtp::VtpEventData &e, int i) -> const vtp::VtpBlock& {
                if (i < 0 || i >= e.n_blocks)
                    throw py::index_error("block index out of range");
                return e.blocks[i];
            },
            py::return_value_policy::reference_internal)
        .def("clear", &vtp::VtpEventData::clear);
}

// -------------------------------------------------------------------------
// TDC (tdc::)
// -------------------------------------------------------------------------
void bind_tdc(py::module_ &m)
{
    py::class_<tdc::TdcHit>(m, "TdcHit")
        .def_readonly("roc_tag", &tdc::TdcHit::roc_tag)
        .def_readonly("slot",    &tdc::TdcHit::slot)
        .def_readonly("channel", &tdc::TdcHit::channel)
        .def_readonly("edge",    &tdc::TdcHit::edge)
        .def_readonly("value",   &tdc::TdcHit::value)
        .def("__repr__", [](const tdc::TdcHit &h) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "<TdcHit roc=0x%x slot=%u ch=%u edge=%u value=%u>",
                h.roc_tag, h.slot, h.channel, h.edge, h.value);
            return std::string(buf);
        });

    py::class_<tdc::TdcEventData, std::shared_ptr<tdc::TdcEventData>>(m, "TdcEventData")
        .def(py::init<>())
        .def_readonly("n_hits", &tdc::TdcEventData::n_hits)
        .def("hit",
            [](const tdc::TdcEventData &e, int i) -> const tdc::TdcHit& {
                if (i < 0 || i >= e.n_hits)
                    throw py::index_error("hit index out of range");
                return e.hits[i];
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly("hits_numpy",
            [](const tdc::TdcEventData &e) {
                // Bulk accessor used by tight Python loops: returns the
                // first n_hits entries of the hits[] buffer as a numpy
                // structured array (copy).  Layout matches the in-memory
                // tdc::TdcHit struct exactly — 12 bytes per row, with one
                // byte of padding between ``edge`` and ``value`` so the
                // uint32 stays 4-byte aligned.
                py::list fields;
                fields.append(py::make_tuple("roc_tag", "<u4"));
                fields.append(py::make_tuple("slot",    "u1"));
                fields.append(py::make_tuple("channel", "u1"));
                fields.append(py::make_tuple("edge",    "u1"));
                fields.append(py::make_tuple("_pad",    "u1"));
                fields.append(py::make_tuple("value",   "<u4"));
                py::dtype dt = py::dtype::from_args(fields);
                // Pass the buffer address; py::array copies when no base
                // handle is given, so the returned array is independent of
                // the TdcEventData (safe across the next decode call).
                const void *src = (e.n_hits > 0) ? (const void*)e.hits
                                                  : nullptr;
                return py::array(dt, (py::ssize_t)e.n_hits, src);
            },
            "Hits as a numpy structured array (fresh copy, length n_hits).\n"
            "Fields: roc_tag <u4, slot u1, channel u1, edge u1, _pad u1, "
            "value <u4.  Use this for bulk per-event loops to avoid "
            "per-attribute Python-to-C++ call overhead.")
        .def("clear", &tdc::TdcEventData::clear);

    // Module-level constants — single source of truth for the TDC LSB +
    // PRad-II RF cabling. Mirror the C++ constexpr names from TdcData.h.
    m.attr("TDC_LSB_NS") = tdc::TDC_LSB_NS;
    m.attr("RF_ROC_TAG") = tdc::RF_ROC_TAG;
    m.attr("RF_SLOT")    = tdc::RF_SLOT;
    m.attr("RF_CH_A")    = tdc::RF_CH_A;
    m.attr("RF_CH_B")    = tdc::RF_CH_B;

    // Compact RF-time view per event.
    py::class_<tdc::RfTimeData, std::shared_ptr<tdc::RfTimeData>>(m, "RfTimeData")
        .def(py::init<>())
        .def_readonly("n_a", &tdc::RfTimeData::n_a)
        .def_readonly("n_b", &tdc::RfTimeData::n_b)
        .def_property_readonly("ns_a", [](const tdc::RfTimeData &r) {
            return py::array(py::dtype("<f4"), (py::ssize_t)r.n_a, r.ns_a);
        }, "Leading-edge times (ns) on RF_CH_A — fresh copy, length n_a.")
        .def_property_readonly("ns_b", [](const tdc::RfTimeData &r) {
            return py::array(py::dtype("<f4"), (py::ssize_t)r.n_b, r.ns_b);
        }, "Leading-edge times (ns) on RF_CH_B — fresh copy, length n_b.")
        .def("nearest_a", &tdc::RfTimeData::nearest_a, py::arg("t_ref_ns"),
             "RF_CH_A tick nearest the reference time (ns); NaN if no hits.")
        .def("nearest_b", &tdc::RfTimeData::nearest_b, py::arg("t_ref_ns"),
             "RF_CH_B tick nearest the reference time (ns); NaN if no hits.")
        .def("clear", &tdc::RfTimeData::clear);

    // Free helpers — use as `prad2py.dec.decode_tdc_replay(roc_tags, nwords,
    // words, evt)` and `prad2py.dec.decode_rf_replay(roc_tags, nwords,
    // words, rf)`. Inputs are array-likes of uint32 (numpy is fine — we
    // copy into vectors at the C++ boundary).
    m.def("decode_tdc_replay",
        [](std::vector<uint32_t> roc_tags,
           std::vector<uint32_t> nwords,
           std::vector<uint32_t> words,
           tdc::TdcEventData &out) {
            return tdc::TdcDecoder::DecodeReplay(roc_tags, nwords, words, out);
        },
        py::arg("roc_tags"), py::arg("nwords"),
        py::arg("words"),    py::arg("out"),
        "Decode the replay tree's flat-triple TDC representation into a "
        "TdcEventData. `out` is cleared first; returns the hit count.");

    m.def("decode_rf_replay",
        [](std::vector<uint32_t> roc_tags,
           std::vector<uint32_t> nwords,
           std::vector<uint32_t> words,
           tdc::RfTimeData &out) {
            tdc::RfTimeDecoder::DecodeReplay(roc_tags, nwords, words, out);
        },
        py::arg("roc_tags"), py::arg("nwords"),
        py::arg("words"),    py::arg("out"),
        "Decode the replay tree's TDC vectors directly into RfTimeData "
        "(filtered to RF_ROC_TAG / RF_SLOT, leading edges only).");
}

// -------------------------------------------------------------------------
// DaqConfig + helpers
// -------------------------------------------------------------------------
void bind_config(py::module_ &m)
{
    py::class_<evc::DaqConfig::RocEntry>(m, "RocEntry")
        .def_readwrite("tag",   &evc::DaqConfig::RocEntry::tag)
        .def_readwrite("name",  &evc::DaqConfig::RocEntry::name)
        .def_readwrite("crate", &evc::DaqConfig::RocEntry::crate)
        .def_readwrite("type",  &evc::DaqConfig::RocEntry::type);

    py::class_<evc::DaqConfig>(m, "DaqConfig",
        "DAQ configuration (bank tags, ROC map, TI layout, ...)")
        .def(py::init<>())
        .def_readwrite("physics_tags",      &evc::DaqConfig::physics_tags)
        .def_readwrite("physics_base",      &evc::DaqConfig::physics_base)
        .def_readwrite("monitoring_tags",   &evc::DaqConfig::monitoring_tags)
        .def_readwrite("prestart_tag",      &evc::DaqConfig::prestart_tag)
        .def_readwrite("go_tag",            &evc::DaqConfig::go_tag)
        .def_readwrite("end_tag",           &evc::DaqConfig::end_tag)
        .def_readwrite("sync_tag",          &evc::DaqConfig::sync_tag)
        .def_readwrite("epics_tag",         &evc::DaqConfig::epics_tag)
        .def_readwrite("adc_format",        &evc::DaqConfig::adc_format)
        .def_readwrite("sparsify_sigma",    &evc::DaqConfig::sparsify_sigma)
        .def_readwrite("fadc_composite_tag",&evc::DaqConfig::fadc_composite_tag)
        .def_readwrite("adc1881m_bank_tag", &evc::DaqConfig::adc1881m_bank_tag)
        .def_readwrite("ti_bank_tag",       &evc::DaqConfig::ti_bank_tag)
        .def_readwrite("trigger_bank_tag",  &evc::DaqConfig::trigger_bank_tag)
        .def_readwrite("run_info_tag",      &evc::DaqConfig::run_info_tag)
        .def_readwrite("daq_config_tag",    &evc::DaqConfig::daq_config_tag)
        .def_readwrite("epics_bank_tag",    &evc::DaqConfig::epics_bank_tag)
        .def_readwrite("ssp_bank_tags",     &evc::DaqConfig::ssp_bank_tags)
        .def_readwrite("fadc_raw_tag",      &evc::DaqConfig::fadc_raw_tag)
        .def_readwrite("tdc_bank_tag",      &evc::DaqConfig::tdc_bank_tag)
        .def_readwrite("roc_tags",          &evc::DaqConfig::roc_tags)
        .def_readwrite("ti_master_tag",     &evc::DaqConfig::ti_master_tag)
        .def_readwrite("verbose_decode",    &evc::DaqConfig::verbose_decode)
        .def_readwrite("fadc250_fw",        &evc::DaqConfig::fadc250_fw)
        .def_readwrite("wave_cfg",          &evc::DaqConfig::wave_cfg)
        .def("is_physics",    &evc::DaqConfig::is_physics)
        .def("is_monitoring", &evc::DaqConfig::is_monitoring)
        .def("is_control",    &evc::DaqConfig::is_control)
        .def("is_sync",       &evc::DaqConfig::is_sync)
        .def("is_epics",      &evc::DaqConfig::is_epics)
        .def("is_ssp_bank",   &evc::DaqConfig::is_ssp_bank);

    m.def("load_daq_config",
        [](const std::string &path) {
            evc::DaqConfig cfg;
            std::string p = path.empty() ? default_daq_config_path() : path;
            if (!evc::load_daq_config(p, cfg))
                throw std::runtime_error("Failed to load DAQ config: " + p);
            return cfg;
        },
        py::arg("path") = std::string(""),
        "Load a DaqConfig from JSON. Empty path uses the installed default.");
}

// -------------------------------------------------------------------------
// EventType / Status enums
// -------------------------------------------------------------------------
void bind_enums(py::module_ &m)
{
    py::enum_<evc::EventType>(m, "EventType")
        .value("Unknown",  evc::EventType::Unknown)
        .value("Physics",  evc::EventType::Physics)
        .value("Sync",     evc::EventType::Sync)
        .value("Epics",    evc::EventType::Epics)
        .value("Prestart", evc::EventType::Prestart)
        .value("Go",       evc::EventType::Go)
        .value("End",      evc::EventType::End)
        .value("Control",  evc::EventType::Control);

    py::enum_<evc::status>(m, "Status")
        .value("failure",    evc::status::failure)
        .value("success",    evc::status::success)
        .value("incomplete", evc::status::incomplete)
        .value("empty",      evc::status::empty)
        .value("eof",        evc::status::eof);
}

// -------------------------------------------------------------------------
// EvChannel
// -------------------------------------------------------------------------
void bind_channel(py::module_ &m)
{
    py::class_<evc::EvChannel>(m, "EvChannel",
        "Evio event reader and scanner.")
        .def(py::init<size_t>(),
            py::arg("buflen") = 1024u * 2000u,
            "Construct with an internal buffer of `buflen` uint32 words.")
        .def("set_config", &evc::EvChannel::SetConfig, py::arg("cfg"))
        .def("get_config", &evc::EvChannel::GetConfig,
             py::return_value_policy::reference_internal)
        .def("open_sequential",
            [](evc::EvChannel &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.OpenSequential(path);
            },
            py::arg("path"),
            "Open an evio file in sequential mode.  Pairs with "
            "open_random_access() and open_auto() — most callers should "
            "prefer open_auto() which picks the best mode automatically.")
        .def("close", &evc::EvChannel::Close)
        .def("read",
            [](evc::EvChannel &self) {
                py::gil_scoped_release rel;
                return self.Read();
            },
            "Read the next record into the internal buffer. Returns Status.")
        // ---- Random-access mode (evio "ra") -------------------------------
        .def("open_random_access",
            [](evc::EvChannel &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.OpenRandomAccess(path);
            },
            py::arg("path"),
            "Open an evio file in random-access mode.  evio mmaps the file "
            "and builds an event pointer table during open — after this you "
            "can jump to any event via read_event_by_index().  Use "
            "get_random_access_event_count() to get the total.")
        .def("get_random_access_event_count",
            &evc::EvChannel::GetRandomAccessEventCount,
            "Total number of events in the random-access table.  Returns 0 "
            "if not opened in random-access mode.")
        .def("read_event_by_index",
            [](evc::EvChannel &self, int i) {
                py::gil_scoped_release rel;
                return self.ReadEventByIndex(i);
            },
            py::arg("index"),
            "Read the event at 0-based evio index `i` into the internal "
            "buffer.  scan() / select_event() / info() / fadc() / ... then "
            "work identically to the sequential path.  Returns Status.")
        .def("open_auto",
            [](evc::EvChannel &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.OpenAuto(path);
            },
            py::arg("path"),
            "Open ``path`` with random-access mode if the file supports it, "
            "otherwise fall back to the sequential mode.  After success, "
            "``is_random_access()`` reports which mode was selected; callers "
            "dispatch between ``read_event_by_index()`` (RA) and ``read()`` "
            "(sequential).  Recommended default when you don't care which "
            "mode you get.")
        .def("is_random_access", &evc::EvChannel::IsRandomAccess,
            "True iff the current handle was opened in random-access mode.")
        .def("scan", &evc::EvChannel::Scan,
            "Scan the currently-held record. Call after a successful Read().")
        .def("get_event_type", &evc::EvChannel::GetEventType)
        .def("get_n_events",   &evc::EvChannel::GetNEvents)
        // ---- Lazy per-product accessors (new API) -------------------------
        // select_event() picks the sub-event; info/fadc/gem/tdc/vtp each
        // decode on first call and return the cached result on repeat calls.
        .def("select_event", &evc::EvChannel::SelectEvent, py::arg("i") = 0,
            "Select the sub-event index for subsequent info/fadc/gem/tdc/vtp "
            "calls.  Clears the product cache if the index changed.  Must be "
            "called after scan(); for PRad-II single-event data use i=0.")
        .def("info",
            [](const evc::EvChannel &self) {
                fdec::EventInfo out;
                {
                    py::gil_scoped_release rel;
                    out = self.Info();
                }
                return out;
            },
            "Decode (or return cached) event info for the currently-selected "
            "sub-event.  Cheapest accessor — skips FADC/SSP/VTP/TDC work.")
        .def("fadc",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<fdec::EventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Fadc();
                }
                return evt;
            },
            "Decode (or return cached) FADC250/ADC1881M waveforms as "
            "EventData (contains event info + per-ROC/slot/channel samples).")
        .def("gem",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<ssp::SspEventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Gem();
                }
                return evt;
            },
            "Decode (or return cached) SSP/MPD GEM strip data as SspEventData.")
        .def("tdc",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<tdc::TdcEventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Tdc();
                }
                return evt;
            },
            "Decode (or return cached) V1190 TDC tagger hits as TdcEventData.")
        .def("vtp",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<vtp::VtpEventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Vtp();
                }
                return evt;
            },
            "Decode (or return cached) VTP ECAL peaks/clusters as VtpEventData.")
        .def("decode_event",
            [](const evc::EvChannel &self, int i,
               bool with_ssp, bool with_vtp, bool with_tdc) -> py::dict {
                auto evt = std::make_shared<fdec::EventData>();
                std::shared_ptr<ssp::SspEventData> ssp_ptr;
                std::shared_ptr<vtp::VtpEventData> vtp_ptr;
                std::shared_ptr<tdc::TdcEventData> tdc_ptr;
                if (with_ssp) ssp_ptr = std::make_shared<ssp::SspEventData>();
                if (with_vtp) vtp_ptr = std::make_shared<vtp::VtpEventData>();
                if (with_tdc) tdc_ptr = std::make_shared<tdc::TdcEventData>();
                bool ok;
                {
                    py::gil_scoped_release rel;
                    ok = self.DecodeEvent(i, *evt,
                                          ssp_ptr.get(), vtp_ptr.get(), tdc_ptr.get());
                }
                py::dict out;
                out["ok"]    = ok;
                out["event"] = py::cast(evt);
                out["ssp"]   = ssp_ptr ? py::cast(ssp_ptr) : py::none();
                out["vtp"]   = vtp_ptr ? py::cast(vtp_ptr) : py::none();
                out["tdc"]   = tdc_ptr ? py::cast(tdc_ptr) : py::none();
                return out;
            },
            py::arg("i") = 0,
            py::kw_only(),
            py::arg("with_ssp") = false,
            py::arg("with_vtp") = false,
            py::arg("with_tdc") = false,
            "Full decode. Returns {'ok': bool, 'event': EventData, "
            "'ssp': SspEventData|None, 'vtp': ..., 'tdc': ...}.")
        .def("sync",
            [](const evc::EvChannel &self) {
                psync::SyncInfo out;
                {
                    py::gil_scoped_release rel;
                    out = self.Sync();
                }
                return out;
            },
            "Absolute-time / run-state snapshot.  Persists across events — "
            "refreshed only when the channel lands on a SYNC/EPICS or "
            "control event (PRESTART/GO/END), otherwise returns the prior "
            "snapshot.  Diff `sync_counter` against your last-seen value to "
            "detect new SYNC ticks; for control events use `event_tag`.")
        .def("extract_epics_text", &evc::EvChannel::ExtractEpicsText,
            "Raw EPICS payload for the current event (empty if not EPICS).");
}

// -------------------------------------------------------------------------
// EpicsStore — slow-control snapshot accumulator.  Fed via Feed() with the
// raw text from EvChannel::ExtractEpicsText() (or the higher-level
// EvChannel::Epics() result), then queried by event_number for the most
// recent value of a channel.  See prad2dec/include/EpicsStore.h.
// -------------------------------------------------------------------------
void bind_epics(py::module_ &m)
{
    py::class_<epics::EpicsStore::Snapshot>(m, "EpicsSnapshot",
        "One EPICS snapshot (event number + timestamp + channel values).")
        .def_readonly("event_number", &epics::EpicsStore::Snapshot::event_number)
        .def_readonly("timestamp",    &epics::EpicsStore::Snapshot::timestamp)
        .def_readonly("values",       &epics::EpicsStore::Snapshot::values);

    py::class_<epics::EpicsStore>(m, "EpicsStore",
        "Accumulator for EPICS slow-control snapshots.  Feed it raw EPICS "
        "text via feed() whenever an EPICS event arrives (channel discovery "
        "is automatic), then query get_value() for any subsequent physics "
        "event — returns the most recent snapshot at or before that event.")
        .def(py::init<>())
        .def("feed", &epics::EpicsStore::Feed,
             py::arg("event_number"), py::arg("timestamp"), py::arg("text"),
             "Parse a raw EPICS payload and store a snapshot.  Text format: "
             "one `value  channel_name` line per channel.")
        .def("get_value",
            [](const epics::EpicsStore &self, int32_t event_number,
               const std::string &channel) -> py::object {
                float v = 0.f;
                if (self.GetValue(event_number, channel, v))
                    return py::float_(v);
                return py::none();
            },
            py::arg("event_number"), py::arg("channel"),
            "Return the most recent value of `channel` at or before "
            "`event_number`, or None if unknown / not yet seen.")
        .def("find_snapshot",
            [](const epics::EpicsStore &self, int32_t event_number)
                -> const epics::EpicsStore::Snapshot* {
                return self.FindSnapshot(event_number);
            },
            py::arg("event_number"),
            py::return_value_policy::reference_internal,
            "Return the whole snapshot at or before `event_number`, or None.")
        .def("get_channel_count", &epics::EpicsStore::GetChannelCount)
        .def("get_channel_id",    &epics::EpicsStore::GetChannelId,
             py::arg("name"),
             "Return the integer id assigned to `name`, or -1 if unknown.")
        .def("get_channel_name",  &epics::EpicsStore::GetChannelName,
             py::arg("id"), py::return_value_policy::reference_internal)
        .def("get_channel_names", &epics::EpicsStore::GetChannelNames,
             py::return_value_policy::reference_internal)
        .def("get_snapshot_count", &epics::EpicsStore::GetSnapshotCount)
        .def("get_snapshot",       &epics::EpicsStore::GetSnapshot,
             py::arg("index"), py::return_value_policy::reference_internal)
        .def("trim", &epics::EpicsStore::Trim, py::arg("max_count"),
             "Drop oldest snapshots so at most `max_count` remain.")
        .def("clear", &epics::EpicsStore::Clear);
}

// -------------------------------------------------------------------------
// HV archive (hv::HVDecoder, HVSegment)
// -------------------------------------------------------------------------
void bind_hv(py::module_ &m)
{
    // Lookup descriptor enums.  Bound as nested-style attributes so the
    // call sites read `dec.HVKind.VMon` / `dec.HVSide.Nearest`.
    py::enum_<hv::Kind>(m, "HVKind",
        "Trace selector for HVSegment.value_at / nearest / nearest_next. "
        "VMon/DV/V0Set apply to HV channels; the four Booster* values to "
        "booster lanes.")
        .value("VMon",         hv::Kind::VMon)
        .value("DV",           hv::Kind::DV)
        .value("V0Set",        hv::Kind::V0Set)
        .value("BoosterVMon",  hv::Kind::BoosterVMon)
        .value("BoosterIMon",  hv::Kind::BoosterIMon)
        .value("BoosterVSet",  hv::Kind::BoosterVSet)
        .value("BoosterISet",  hv::Kind::BoosterISet);

    py::enum_<hv::Side>(m, "HVSide",
        "Side selector for nearest/next lookup.  Nearest = closest in time; "
        "Next = first snapshot at-or-after the query.")
        .value("Nearest", hv::Side::Nearest)
        .value("Next",    hv::Side::Next);

    py::class_<hv::LookupResult>(m, "HVLookupResult",
        "(channel, unix_time) lookup result. ok=False when out-of-range or NaN.")
        .def_readonly("ok",       &hv::LookupResult::ok)
        .def_readonly("t_unix_s", &hv::LookupResult::t_unix_s)
        .def_readonly("value",    &hv::LookupResult::value)
        .def("__repr__", [](const hv::LookupResult &r) {
            char buf[128];
            if (!r.ok) std::snprintf(buf, sizeof(buf),
                                     "<HVLookupResult ok=False>");
            else       std::snprintf(buf, sizeof(buf),
                                     "<HVLookupResult t=%.3f v=%.3f>",
                                     r.t_unix_s, r.value);
            return std::string(buf);
        });

    py::class_<hv::Interval>(m, "HVInterval",
        "Closed time interval (unix seconds) returned by find_stable_intervals.")
        .def_readonly("t_start_unix", &hv::Interval::t_start_unix)
        .def_readonly("t_end_unix",   &hv::Interval::t_end_unix)
        .def_property_readonly("duration_s",
            [](const hv::Interval &iv) { return iv.t_end_unix - iv.t_start_unix; })
        .def("__repr__", [](const hv::Interval &iv) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<HVInterval %.3f → %.3f (%.2fs)>",
                iv.t_start_unix, iv.t_end_unix,
                iv.t_end_unix - iv.t_start_unix);
            return std::string(buf);
        });

    py::class_<hv::ChEvent>(m, "HVChEvent",
        "One CHTABLE event: an absolute timestamp plus the V0Set vector "
        "(projected to the segment's kept-channel order).")
        .def_readonly("abs_ts_ms", &hv::ChEvent::abs_ts_ms)
        .def_property_readonly("v0sets",
            [](const hv::ChEvent &e) {
                py::array_t<float> arr(e.v0sets.size());
                std::copy(e.v0sets.begin(), e.v0sets.end(), arr.mutable_data());
                return arr;
            },
            "V0Set per channel as a fresh numpy array.");

    py::class_<hv::BstEvent>(m, "HVBstEvent",
        "One BOOSTER_TABLE event: timestamp + setpoints (VSet / ISet).")
        .def_readonly("abs_ts_ms", &hv::BstEvent::abs_ts_ms)
        .def_property_readonly("vsets",
            [](const hv::BstEvent &e) {
                py::array_t<float> arr(e.vsets.size());
                std::copy(e.vsets.begin(), e.vsets.end(), arr.mutable_data());
                return arr;
            })
        .def_property_readonly("isets",
            [](const hv::BstEvent &e) {
                py::array_t<float> arr(e.isets.size());
                std::copy(e.isets.begin(), e.isets.end(), arr.mutable_data());
                return arr;
            });

    // ── HVSegment ────────────────────────────────────────────────────────
    // Bulk arrays go through allocate-then-copy (fresh numpy buffers, no
    // base pointer).  This is the same pattern ChannelData.samples uses;
    // exposing a pybind11 view onto the underlying C++ vector silently
    // corrupts memory once the segment is freed (see the codebase
    // "py::array_t binding ownership trap" memory).
    py::class_<hv::HVSegment>(m, "HVSegment",
        "Time-windowed HV + booster snapshot block.\n\n"
        "Numeric arrays returned as fresh numpy copies — the underlying "
        "C++ vectors are not aliased into the numpy buffer.")
        .def(py::init<>())
        .def_readonly("interval_ms",       &hv::HVSegment::interval_ms)
        .def_readonly("window_start_unix", &hv::HVSegment::window_start_unix)
        .def_readonly("window_end_unix",   &hv::HVSegment::window_end_unix)
        .def_readonly("source_files",      &hv::HVSegment::source_files)
        .def_property_readonly("channels",
            [](const hv::HVSegment &s) { return s.channels; })
        .def_property_readonly("booster_names",
            [](const hv::HVSegment &s) { return s.booster_names; })
        .def_property_readonly("n_channels",          &hv::HVSegment::n_channels)
        .def_property_readonly("n_snapshots",         &hv::HVSegment::n_snapshots)
        .def_property_readonly("n_boosters",          &hv::HVSegment::n_boosters)
        .def_property_readonly("n_booster_snapshots", &hv::HVSegment::n_booster_snapshots)
        .def_property_readonly("empty",               &hv::HVSegment::empty)
        .def_property_readonly("t_start_unix",
            [](const hv::HVSegment &s) -> py::object {
                if (s.n_snapshots() == 0) return py::none();
                return py::float_(double(s.timestamps_ms.front()) / 1000.0);
            })
        .def_property_readonly("t_end_unix",
            [](const hv::HVSegment &s) -> py::object {
                if (s.n_snapshots() == 0) return py::none();
                return py::float_(double(s.timestamps_ms.back()) / 1000.0);
            })

        // 1-D numpy arrays
        .def_property_readonly("timestamps_ms",
            [](const hv::HVSegment &s) {
                py::array_t<int64_t> arr(s.timestamps_ms.size());
                std::copy(s.timestamps_ms.begin(), s.timestamps_ms.end(),
                          arr.mutable_data());
                return arr;
            },
            "Per-snapshot epoch-ms timestamps (int64).")
        .def_property_readonly("booster_timestamps_ms",
            [](const hv::HVSegment &s) {
                py::array_t<int64_t> arr(s.booster_timestamps_ms.size());
                std::copy(s.booster_timestamps_ms.begin(),
                          s.booster_timestamps_ms.end(), arr.mutable_data());
                return arr;
            })

        // 2-D numpy arrays — row-major (n_rows, n_cols)
        .def_property_readonly("dv",
            [](const hv::HVSegment &s) {
                py::array_t<float> arr({(py::ssize_t)s.n_snapshots(),
                                        (py::ssize_t)s.n_channels()});
                if (!s.dv.empty())
                    std::copy(s.dv.begin(), s.dv.end(), arr.mutable_data());
                return arr;
            },
            "VMon - V0Set per (snapshot, channel) as a fresh (N_snap × N_ch) "
            "float32 array.")
        .def_property_readonly("booster_vmon",
            [](const hv::HVSegment &s) {
                py::array_t<float> arr({(py::ssize_t)s.n_booster_snapshots(),
                                        (py::ssize_t)s.n_boosters()});
                if (!s.booster_vmon.empty())
                    std::copy(s.booster_vmon.begin(), s.booster_vmon.end(),
                              arr.mutable_data());
                return arr;
            })
        .def_property_readonly("booster_imon",
            [](const hv::HVSegment &s) {
                py::array_t<float> arr({(py::ssize_t)s.n_booster_snapshots(),
                                        (py::ssize_t)s.n_boosters()});
                if (!s.booster_imon.empty())
                    std::copy(s.booster_imon.begin(), s.booster_imon.end(),
                              arr.mutable_data());
                return arr;
            })

        // Tables (lists of HVChEvent / HVBstEvent — each event already
        // exposes its arrays as numpy copies)
        .def_readonly("ch_events",      &hv::HVSegment::ch_events)
        .def_readonly("booster_events", &hv::HVSegment::booster_events)

        // Index resolution
        .def("channel_index", &hv::HVSegment::channel_index,
             py::arg("name"),
             "Return position of `name` in channels[], or -1 if absent.")
        .def("booster_index", &hv::HVSegment::booster_index,
             py::arg("name"),
             "Return position of `name` in booster_names[], or -1 if absent.")

        // Reconstructed traces
        .def("v0set_trace",
            [](const hv::HVSegment &s, int ch_idx) {
                std::vector<float> v = s.v0set_trace(ch_idx);
                py::array_t<float> arr(v.size());
                std::copy(v.begin(), v.end(), arr.mutable_data());
                return arr;
            },
            py::arg("ch_idx"),
            "V0Set per snapshot for kept channel `ch_idx` (NaN if no "
            "CHTABLE seen for this channel).")
        .def("vmon_trace",
            [](const hv::HVSegment &s, const std::string &name) {
                int i = s.channel_index(name);
                if (i < 0)
                    throw py::key_error("HV channel not found: " + name);
                std::vector<float> v0 = s.v0set_trace(i);
                py::array_t<int64_t> ts(s.timestamps_ms.size());
                py::array_t<float>   vmon(v0.size());
                std::copy(s.timestamps_ms.begin(), s.timestamps_ms.end(),
                          ts.mutable_data());
                const int n_ch = s.n_channels();
                for (std::size_t k = 0; k < v0.size(); ++k)
                    vmon.mutable_data()[k] =
                        s.dv[k * n_ch + i] + v0[k];
                return py::make_tuple(ts, vmon);
            },
            py::arg("name"),
            "Return (timestamps_ms, vmon) numpy arrays for HV channel `name`. "
            "VMon is reconstructed via the most-recent CHTABLE V0Set.")
        .def("dv_trace",
            [](const hv::HVSegment &s, const std::string &name) {
                int i = s.channel_index(name);
                if (i < 0)
                    throw py::key_error("HV channel not found: " + name);
                py::array_t<int64_t> ts(s.timestamps_ms.size());
                py::array_t<float>   dv(s.timestamps_ms.size());
                std::copy(s.timestamps_ms.begin(), s.timestamps_ms.end(),
                          ts.mutable_data());
                const int n_ch = s.n_channels();
                for (std::size_t k = 0; k < s.timestamps_ms.size(); ++k)
                    dv.mutable_data()[k] = s.dv[k * n_ch + i];
                return py::make_tuple(ts, dv);
            },
            py::arg("name"),
            "Return (timestamps_ms, dV) numpy arrays for HV channel `name`.")

        // (channel, unix_time) lookup
        .def("value_at", &hv::HVSegment::value_at,
             py::arg("name"), py::arg("unix_time"),
             py::arg("kind") = hv::Kind::VMon,
             py::arg("side") = hv::Side::Nearest,
             "Look up a channel's value at (or near) a wall-clock time.\n"
             "Returns HVLookupResult — check .ok before reading "
             ".t_unix_s / .value.")
        .def("nearest", &hv::HVSegment::nearest,
             py::arg("name"), py::arg("unix_time"),
             py::arg("kind") = hv::Kind::VMon,
             "Snapshot closest in absolute time. HVLookupResult; .ok=False "
             "when out of range or NaN.")
        .def("nearest_next", &hv::HVSegment::nearest_next,
             py::arg("name"), py::arg("unix_time"),
             py::arg("kind") = hv::Kind::VMon,
             "First snapshot at-or-after the query time. HVLookupResult; "
             ".ok=False if the query is past the last snapshot.")

        // (event_num, event_ti_ticks, anchor) → per-snapshot event tag.
        .def_static("associate_events",
            [](py::array_t<int64_t> snapshot_ts_ms,
               py::array_t<int32_t> event_num,
               py::array_t<int64_t> event_ti_ticks,
               int64_t anchor_ti_ticks,
               int64_t anchor_unix_time_ms) {
                auto sn = snapshot_ts_ms .unchecked<1>();
                auto en = event_num      .unchecked<1>();
                auto tt = event_ti_ticks .unchecked<1>();
                if (en.shape(0) != tt.shape(0))
                    throw std::invalid_argument(
                        "associate_events: event_num and event_ti_ticks "
                        "must be the same length");
                std::vector<int64_t> snaps(sn.shape(0));
                for (py::ssize_t i = 0; i < sn.shape(0); ++i)
                    snaps[i] = sn(i);
                std::vector<int32_t> ens(en.shape(0));
                std::vector<int64_t> tts(tt.shape(0));
                for (py::ssize_t i = 0; i < en.shape(0); ++i) {
                    ens[i] = en(i);
                    tts[i] = tt(i);
                }
                auto assoc = hv::HVSegment::associate_events(
                    snaps, ens, tts,
                    anchor_ti_ticks, anchor_unix_time_ms);

                py::array_t<int32_t> ev_arr(assoc.event_number.size());
                py::array_t<int64_t> tk_arr(assoc.ti_ticks.size());
                std::copy(assoc.event_number.begin(),
                          assoc.event_number.end(),  ev_arr.mutable_data());
                std::copy(assoc.ti_ticks.begin(),
                          assoc.ti_ticks.end(),      tk_arr.mutable_data());
                return py::make_tuple(ev_arr, tk_arr);
            },
            py::arg("snapshot_ts_ms"),
            py::arg("event_num"),
            py::arg("event_ti_ticks"),
            py::arg("anchor_ti_ticks"),
            py::arg("anchor_unix_time_ms"),
            "Associate each entry of `snapshot_ts_ms` (epoch ms — typically "
            "seg.timestamps_ms or seg.booster_timestamps_ms) with the most "
            "recent physics event.\n\n"
            "Inputs:\n"
            "  snapshot_ts_ms      : int64 1-D array of times to tag.\n"
            "  event_num           : int32 1-D array, recon `event_num` branch.\n"
            "  event_ti_ticks      : int64 1-D array, recon `timestamp` branch\n"
            "                        (TI 4-ns ticks, sorted ascending).\n"
            "  anchor_ti_ticks,\n"
            "  anchor_unix_time_ms : any SYNC scaler row's (ti_ticks, unix_time*1000)\n"
            "                        — pins the linear ti_ticks ↔ unix mapping.\n\n"
            "Returns (event_number_per_snap int32, ti_ticks_per_snap int64).\n"
            "Snapshots before the first event get -1 sentinels; snapshots "
            "after the last event clamp to the last event.")

        // Stable-interval finder
        .def("find_stable_intervals",
            [](const hv::HVSegment &self,
               const std::vector<std::string> &channels,
               double window_s,
               double std_threshold,
               py::object dv_threshold,
               double min_duration_s,
               double guard_s) {
                std::optional<double> dvt;
                if (!dv_threshold.is_none())
                    dvt = dv_threshold.cast<double>();
                return self.find_stable_intervals(channels, window_s,
                                                  std_threshold, dvt,
                                                  min_duration_s, guard_s);
            },
            py::arg("channels"),
            py::kw_only(),
            py::arg("window_s") = 5.0,
            py::arg("std_threshold") = 0.5,
            py::arg("dv_threshold") = py::none(),
            py::arg("min_duration_s") = 5.0,
            py::arg("guard_s") = 1.0,
            "Find time intervals where ALL named channels are stable.\n\n"
            "A snapshot is unstable if for ANY of `channels`:\n"
            "  - rolling-std(dV) over `window_s` > std_threshold, OR\n"
            "  - dV is NaN, OR\n"
            "  - |dV| > dv_threshold (only if dv_threshold is not None).\n"
            "Stable runs shorter than `min_duration_s` are dropped; "
            "`guard_s` is trimmed from each end.")

        // Persistence
        .def("save", &hv::HVSegment::save, py::arg("path"),
             "Write the segment to a compact binary cache (magic VMHV0001).")
        .def_static("load", &hv::HVSegment::load, py::arg("path"),
             "Load a cache previously written by save().")
        .def("__repr__", [](const hv::HVSegment &s) {
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                "<HVSegment n_ch=%d n_snap=%d n_bst=%d n_bsnap=%d "
                "interval=%dms>",
                s.n_channels(), s.n_snapshots(),
                s.n_boosters(), s.n_booster_snapshots(),
                s.interval_ms);
            return std::string(buf);
        });

    // ── HVDecoder ────────────────────────────────────────────────────────
    py::class_<hv::HVDecoder>(m, "HVDecoder",
        "Open a directory (or list) of vmon_*.dat files and load windowed "
        "HVSegment instances.\n\n"
        "Construction is cheap (filename enumeration only); "
        "load_window() does the actual I/O.  The mmap-based parser only "
        "materializes the requested rows × columns, so a few-channel "
        "few-minute query against a multi-GB daily file stays bounded "
        "in memory.")
        .def(py::init<const std::string &>(), py::arg("source"),
             "Open a directory of vmon_*.dat files (auto-discover) or a "
             "single .dat path.")
        .def(py::init<const std::vector<std::string> &>(), py::arg("files"),
             "Open a curated list of vmon_*.dat paths.")
        .def_property_readonly("files", &hv::HVDecoder::files,
             "Discovered file list (sorted).")
        .def("load_window", &hv::HVDecoder::load_window,
             py::arg("t_start_unix"), py::arg("t_end_unix"),
             py::arg("channels") = std::vector<std::string>{},
             py::arg("cache_path") = std::string(),
             "Load HV data for [t_start_unix, t_end_unix) seconds.  "
             "Empty `channels` keeps every channel in the archive.  When "
             "`cache_path` is given and exists, the cache is returned "
             "verbatim (channels/window args not re-validated).  Delete "
             "the cache file to force a refetch.");
}

} // anonymous namespace

// -------------------------------------------------------------------------
// Entry point for the main module (prad2py.cpp calls this).
// -------------------------------------------------------------------------
void register_dec(py::module_ &m)
{
    auto dec = m.def_submodule("dec",
        "prad2dec bindings: evio reader, event data types, TDC/SSP/VTP accessors.");

    bind_enums(dec);
    bind_config(dec);
    bind_fadc(dec);
    bind_ssp(dec);
    bind_vtp(dec);
    bind_tdc(dec);
    bind_epics(dec);
    bind_channel(dec);
    bind_hv(dec);
}
