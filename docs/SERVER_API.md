# PRad-II Event Viewer — Server API Reference

Base URL: `http://localhost:<port>` (default port 5051).

All responses are JSON unless noted otherwise. The server also pushes real-time updates over a WebSocket connection on the same port.

Large payloads are gzipped on demand when the client advertises `Accept-Encoding: gzip` (browsers always do; bare `urllib`/`curl` clients receive plain bytes).

---

## Server

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/config` | Server configuration and capabilities (mode, ET status, source caps, file metadata, filter state, detector geometry) |
| GET | `/api/progress` | File-loading progress (`loading`, `phase`, `current`, `total`, `file`) |
| GET | `/api/monitor_status` | Header-panel status: `livetime.{ts,measured}` and `beam.{energy,current}`. Each value is `<0` when unavailable |
| POST | `/api/hist_config` | Update waveform peak filter at runtime. JSON body: `{"waveform_filter": {...}, "waveform_filter_active": <bool>}`. Clears accumulated histograms and broadcasts `hist_config_updated` + `hist_cleared` |

## Mode Switching

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/mode/online` | Switch to ET/online mode. Optional JSON body: `{"host","port","et_file","station"}` |
| POST | `/api/mode/file` | Switch to file/offline mode |

## Events

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/event/<n>` | Decoded event `n` (1-based in file mode; seq number in online mode) |
| GET | `/api/event/latest` | Latest event from the ring buffer (online mode) |
| GET | `/api/waveform/<n>/<roc_slot_ch>` | On-demand waveform samples for a single channel (file mode) |
| GET | `/api/clusters/<n>` | Cluster reconstruction for event `n` |
| GET | `/api/ring` | Ring buffer summary: `{"ring": [seq…], "latest": <seq>}` |

## File Browser

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/files[?dir=<subdir>]` | List entries under `data_dir` (or under `<subdir>` relative to it). Each entry is `{type:"dir",name,count}` or `{type:"file",name,size,…}` |
| GET | `/api/load?file=<path>&hist=0\|1` | Load an evio file (relative to `data_dir`). `hist=1` enables histogram preprocessing |

## Histograms & Occupancy

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/occupancy` | Per-module occupancy (hit counts and integrals) |
| GET | `/api/cluster_hist` | Cluster-level histograms |
| GET | `/api/hist/<module_key>` | Integral (amplitude) histogram for a module |
| GET | `/api/poshist/<module_key>` | Position histogram for a module |
| GET | `/api/heighthist/<module_key>` | Peak-height histogram for a module |
| POST | `/api/hist/clear` | Clear accumulated histograms (broadcasts `hist_cleared`) |

## LMS (Laser Monitoring)

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/lms/summary[?ref=<ch>]` | LMS summary across all modules. Optional `ref` channel for normalization |
| GET | `/api/lms/<module>[?ref=<ch>]` | LMS time series for a single module |
| GET | `/api/lms/refs` | List available LMS reference channels |
| POST | `/api/lms/clear` | Clear LMS accumulator (broadcasts `lms_cleared`) |

## EPICS

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/epics/channels` | List of known EPICS channel names |
| GET | `/api/epics/latest` | Latest values for all EPICS channels |
| GET | `/api/epics/channel/<name>` | Time series for a single EPICS channel (URL-decoded) |
| GET | `/api/epics/batch?ch=<n1>&ch=<n2>` | Batch fetch multiple channels |
| POST | `/api/epics/clear` | Clear EPICS history (broadcasts `epics_cleared`) |

## GEM

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/gem/config` | GEM detector geometry and strip mapping |
| GET | `/api/gem/hits` | GEM hits for the current event |
| GET | `/api/gem/occupancy` | GEM strip occupancy histograms |
| GET | `/api/gem/residuals` | GEM↔HyCal matching residuals |
| GET | `/api/gem/efficiency` | Per-detector tracking-efficiency counters + last-good-event snapshot |
| GET | `/api/gem/calib` | Per-APV pedestal noise + global `zs_sigma` (one-shot; cached on the client and refetched only when `calib_rev` changes) |
| GET | `/api/gem/apv/<n>` | Per-event GEM APV waveforms for event `n` |
| POST | `/api/gem/threshold` | Update zero-suppression threshold. JSON body: `{"zs_sigma": N}`. Broadcasts `gem_threshold_updated` |

## Physics

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/physics/energy_angle` | Energy vs. angle distribution |
| GET | `/api/physics/moller` | Moller scattering analysis. Gated by its own trigger filter (`physics.moller.accept/reject_trigger_bits` in monitor_config.json — X17: 2-cluster trigger; inherits the `physics` filter when unset). Response carries the active masks in `trigger` + `trigger_accept_names` |
| GET | `/api/physics/hycal_xy` | Single-cluster HyCal hit map |

## Elog

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/elog/post` | Post to the electronic logbook. JSON body: `{"xml": "<elog XML>", "run_number": N, "auto": <bool>, "request_id": "<id>"}`. Saves locally first, then uploads. With `auto:true` the `request_id` must match an in-flight `capture_request`. Server-side dedup absorbs duplicates within `auto_report.min_interval_ms` |

## Event Filters

Filters control which events are navigable and accumulated in file mode.
Load from a file (`-f filter.json`) or at runtime via the API / UI panel.

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/filter` | Current filter state |
| POST | `/api/filter/load` | Load filter (JSON body). Rebuilds indices and histograms |
| POST | `/api/filter/unload` | Remove all filters. Rebuilds |
| GET | `/api/filter/indices` | List of 1-based event indices passing the filter |

### Filter JSON format

```json
{
  "waveform": {
    "enable": true,
    "modules": ["W100", "W101"],
    "n_peaks_min": 1, "n_peaks_max": 999999,
    "time_min": 160, "time_max": 220,
    "integral_min": 100, "integral_max": 15000,
    "height_min": 20, "height_max": 4000
  },
  "clustering": {
    "enable": true,
    "n_min": 1, "n_max": 999999,
    "energy_min": 50, "energy_max": 2500,
    "size_min": 1, "size_max": 999999,
    "includes_modules": ["W100", "G25"], "includes_min": 1,
    "center_modules": ["W100", "W101"]
  }
}
```

Each section has `enable: false` by default. All other fields optional.
Events pass if ALL enabled filters return true. Loading/unloading clears
accumulated data and rebuilds histograms if the file was preprocessed.

---

## WebSocket Messages

Connect to `ws://localhost:<port>` for real-time push notifications.

### Client → Server

| `type` | Fields | Description |
|--------|--------|-------------|
| `client_hello` | `capabilities: ["auto_report", …]` | Advertise client capabilities. Clients that include `"auto_report"` become eligible for `capture_request` dispatch |
| `tagger_subscribe` | — | Subscribe to live tagger binary stream |
| `tagger_unsubscribe` | — | Unsubscribe from live tagger stream |

### Server → Client

| `type` | Fields | Description |
|--------|--------|-------------|
| `server_hello` | `capabilities` | Reply to `client_hello` |
| `status` | `connected`, `waiting?`, `retries?` | ET connection status change |
| `new_event` | `seq` | New event available in ring buffer |
| `mode_changed` | `mode` | Mode changed (`file`, `online`, `idle`) |
| `file_loaded` | (full config) | File finished loading |
| `load_progress` | `phase`, `current`, `total` | File load progress update |
| `hist_cleared` | — | Histograms were cleared |
| `hist_config_updated` | `waveform_filter`, `waveform_filter_active` | Waveform peak filter changed |
| `lms_cleared` | — | LMS data was cleared |
| `lms_event` | `count` | New LMS trigger event |
| `epics_cleared` | — | EPICS data was cleared |
| `epics_event` | `count` | New EPICS event received |
| `gem_threshold_updated` | `zs_sigma` | GEM zero-suppression threshold changed |
| `control_event` | `kind` (`prestart`/`end`), `run_number`, `unix_time` | EVIO control event seen on the stream |
| `capture_request` | `request_id`, `run`, `reason` | Server asks an `auto_report`-capable client to take and post an auto-report. Client must reply via `POST /api/elog/post` with `auto:true` and the same `request_id` |
| `auto_capture_done` | `run`, `posted`, `saved_xml` | Auto-report capture handshake finished |
| `tagger_subscribed` | `subscribers` | Confirms tagger subscription |
| `tagger_unsubscribed` | `subscribers` | Confirms tagger unsubscription |

### Tagger binary frames

While subscribed, the server pushes binary WebSocket frames directly to
each tagger subscriber. Each frame begins with a 24-byte header:

```
offset  bytes  field
  0      4    magic "TGR1"
  4      4    flags  (bit 0 = drops occurred)
  8      4    n_hits
 12      4    first_seq
 16      4    last_seq
 20      4    dropped_frames
 24+     n_hits × TaggerBinHit (event_num, trigger_bits, roc_tag, slot, channel|edge, tdc)
```

Hits are batched up to `TAGGER_BATCH_MAX` or the time-based flush interval, whichever comes first.

---

## CLI Interactive Mode

Start the server with `-i` (or `--interactive`) to enable the stdin command interface:

```
prad2_server data.evio -H -i
```

| Command | Description |
|---------|-------------|
| `status` | Show current mode, file info, ET connection |
| `load <path> [1]` | Load an evio file (append `1` for histograms) |
| `online` | Switch to ET/online mode |
| `offline` | Switch to file/offline mode |
| `clear hist\|lms\|epics` | Clear accumulators |
| `filter` | Show current filter state |
| `filter load <f>` | Load event filter from JSON file |
| `filter unload` | Remove all filters |
| `quit` / `exit` | Stop the server |
| `help` / `?` | Show command list |
