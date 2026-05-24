#pragma once
// =========================================================================
// viewer_server.h — Unified HTTP/WebSocket server for PRad-II event viewer
//
// Combines file-based viewing and online ET monitoring into a single server.
// Mode switching between "idle", "file", and "online" via API or user actions.
// =========================================================================

#include "data_source.h"
#include "app_state.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "TdcData.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <atomic>
#include <chrono>
#include <ctime>
#include <deque>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

// ── Shared types ─────────────────────────────────────────────────────────

struct FileData {
    std::string filepath;
    int event_count = 0;
    DataSourceCaps caps;
};

struct Progress {
    std::atomic<bool> loading{false};
    std::atomic<int>  phase{0};      // 0=idle, 1=indexing, 2=histograms
    std::atomic<int>  current{0};
    std::atomic<int>  total{0};
    std::string       target_file;
    mutable std::mutex mtx;

    nlohmann::json toJson() const;
    void setFile(const std::string &f);
};

struct RingEntry {
    int seq;
    std::string json_str;       // pre-encoded event JSON
    std::string cluster_str;    // pre-encoded cluster JSON
    std::string gem_apv_str;    // pre-encoded GEM per-APV waveform JSON
    // Pre-compressed gzip bytes for the big payload — deflated once per
    // event in the ET reader thread so each viewer's HTTP request just
    // copies the cached blob instead of re-running zlib.  Empty when
    // gem_apv_str is empty (e.g. GEM disabled).
    std::string gem_apv_gz;

    // Raw event copies kept so /api/hist_config can recompute cluster_str
    // (and re-encode json_str) under a new time/threshold window without
    // waiting for new events to arrive.  ~8MB per entry — bounded by ring_size_.
    std::shared_ptr<fdec::EventData>    event_data;
    std::shared_ptr<ssp::SspEventData>  ssp_data;
};

// ── ViewerServer ─────────────────────────────────────────────────────────

class ViewerServer {
public:
    using WsServer = websocketpp::server<websocketpp::config::asio>;

    struct Config {
        std::string database_dir;
        std::string resource_dir;
        std::string data_dir;           // file browsing root (empty = disabled)
        std::string daq_config_file;          // empty = auto-find in database_dir
        std::string monitor_config_file;      // empty = auto-find (monitor_config.json)
        std::string reconstruction_config_file; // empty = auto-find (reconstruction_config.json)
        std::string initial_file;       // .evio file to open on startup
        int    port         = 5051;
        bool   hist_enabled = false;
        bool   start_online = false;    // connect ET on startup
        bool   interactive  = false;    // enable stdin command loop
        std::string filter_file;        // external filter JSON (-f)
    };

    ViewerServer();
    ~ViewerServer();

    // Initialize application state. Must be called before run/startAsync.
    void init(const Config &cfg);

    // Run the server (blocking). Loads initial file, then serves.
    void run();

    // Start server in a background thread. Returns the actual port.
    int startAsync(int port = 0);

    // Stop the server and all background threads.
    void stop();

    // Load a file by absolute path. Switches to file mode.
    // Non-blocking: spawns a background load thread.
    void loadFile(const std::string &path, bool hist);

    int  port() const { return port_; }
    std::string mode() const;
    bool isLoading() const { return progress_.loading.load(); }
    nlohmann::json getProgress() const { return progress_.toJson(); }

    // Active AppState for the current mode.
    AppState &activeApp();

private:
    // ── Mode ─────────────────────────────────────────────────────────────
    enum class Mode { Idle, File, Online };
    std::atomic<Mode> mode_{Mode::Idle};
    std::mutex mode_mtx_;       // serialises mode transitions

    void setMode(Mode m);       // set + broadcast

    // ── Dual AppState (file vs online, never mixed) ──────────────────────
    AppState    app_file_;
    AppState    app_online_;
    Config      cfg_;
    std::string res_dir_;
    int         port_ = 0;
    std::atomic<bool> running_{true};

    // ── WebSocket ────────────────────────────────────────────────────────
    std::unique_ptr<WsServer> server_;
    std::thread server_thread_;
    std::set<websocketpp::connection_hdl,
             std::owner_less<websocketpp::connection_hdl>> ws_clients_;
    // Clients that have advertised the on-demand auto-report protocol via
    // a client_hello message. Pre-update tabs (no hello sent) stay out of
    // this set, so dispatchCapture skips them — the watchdog never burns
    // 30 s on a client that wouldn't know what to do with capture_request.
    std::set<websocketpp::connection_hdl,
             std::owner_less<websocketpp::connection_hdl>> reporter_capable_;
    std::mutex ws_mtx_;

    void wsBroadcast(const std::string &msg);
    void handleWsMessage(websocketpp::connection_hdl hdl,
                         const std::string &payload);

    // ── Tagger live stream (binary broadcast to subscribed WebSocket clients) ─
    // Zero-cost when no one is subscribed: the ET reader only runs the TDC
    // decoder when tagger_subs_count_ > 0.  Frames are a 24-byte header
    // followed by N × 16B packed BinHit records — see viewer_server_et.cpp
    // for the exact layout (the Python client in scripts/tagger_viewer.py
    // mirrors it).
    std::set<websocketpp::connection_hdl,
             std::owner_less<websocketpp::connection_hdl>> tagger_subs_;
    std::mutex                 tagger_subs_mtx_;
    std::atomic<int>           tagger_subs_count_{0};
    std::atomic<uint64_t>      tagger_dropped_frames_{0};  // incremented on per-subscriber send failure

    void taggerSubscribe(websocketpp::connection_hdl hdl);
    void taggerUnsubscribe(websocketpp::connection_hdl hdl);
    void taggerBroadcastBinary(const void *data, size_t nbytes);

    // ── File mode ────────────────────────────────────────────────────────
    std::shared_ptr<FileData> file_data_;
    std::mutex file_data_mtx_;

    std::unique_ptr<DataSource> data_source_;
    std::mutex data_source_mtx_;
    std::unordered_map<int, uint32_t> crate_to_roc_;  // for ROOT data sources

    mutable Progress progress_;
    std::atomic<bool> hist_enabled_{false};
    std::thread load_thread_;
    std::mutex load_mtx_;

    // On-demand accumulation: mirrors the online-mode logic for file browsing.
    // - Preprocessed (hist_enabled_): all events already processed by
    //   buildHistograms(); further calls are no-ops.
    // - Not preprocessed: processEvent/processGemEvent are called once per
    //   event as the user browses (deduped by ondemand_processed_).
    // Any new accumulation added to processEvent() automatically follows
    // this pattern — no per-endpoint code is needed.
    std::unordered_set<int> ondemand_processed_;
    std::mutex ondemand_mtx_;
    void accumulate(int ev1, fdec::EventData &event, ssp::SspEventData *ssp);

    // ── Filters ──────────────────────────────────────────────────────────
    std::vector<int> filtered_indices_;   // 1-based event indices passing filter
    void buildFilteredIndex();
    std::string applyFilter(const nlohmann::json &fj);
    void clearFilter();

    void buildHistograms();
    void loadFileInternal(const std::string &filepath);

    std::string decodeRawEvent(int ev1, fdec::EventData &event,
                               ssp::SspEventData *ssp_evt = nullptr);
    nlohmann::json decodeEvent(int ev1);
    nlohmann::json computeClusters(int ev1);

    // ── Online mode (ET) ─────────────────────────────────────────────────
#ifdef WITH_ET
    struct EtCfg {
        std::string host    = "localhost";
        int         port    = 11111;
        std::string et_file = "/tmp/et_sys_prad2";
        std::string station = "prad2_monitor";
    } et_cfg_;

    int ring_size_ = 20;
    std::deque<RingEntry> ring_;
    std::mutex ring_mtx_;

    // Latest full-readout snapshot — most recent monitoring event where the
    // firmware bypassed online ZS for every APV (or at least one).  Encoded
    // by the ET reader with skip_sw_zs=true so the client's signal-only
    // filter doesn't hide the entire pedestal/noise spectrum.  Served by
    // /api/gem/apv/latest_full and announced over WS as gem_apv_full_event
    // so the gem_apv tab can refresh selectively (it ignores the regular
    // per-event new_event in 'Latest full-readout' source mode).
    std::mutex  latest_full_apv_mtx_;
    std::string latest_full_apv_json_;
    std::string latest_full_apv_gz_;

    std::atomic<bool> et_active_{false};
    std::atomic<bool> et_connected_{false};
    std::atomic<int>  et_generation_{0};    // bumped to trigger reconnect
    std::thread et_thread_;

    void etReaderThread();
    void sleepMs(int ms);

    // Monitor-mode header status (livetime + beam energy + current).  All
    // three are <0 when unavailable so the frontend can hide each cell
    // independently.  A single background thread polls every configured
    // metric on its own cadence (per-metric poll_sec); a metric whose
    // command is empty stays at -1.0.  Avoids a build-time EPICS dependency
    // by shelling out to whatever tool the host provides (typically caget).
    //   - livetime_       ← AppState::livetime_cmd        (TS, percent)
    //   - beam_energy_    ← AppState::beam_energy_status  (MeV)
    //   - beam_current_   ← AppState::beam_current_status (nA)
    // The DSC2-derived "measured" livetime companion lives on
    // AppState::measured_livetime, populated from the EVIO stream.
    std::atomic<double> livetime_{-1.0};
    std::atomic<double> beam_energy_{-1.0};
    std::atomic<double> beam_current_{-1.0};
    std::thread         monitor_status_thread_;
    void                monitorStatusPollThread();
#endif

    void joinAll();      // join all background threads (safe to call multiple times)
    void commandLoop();  // interactive stdin command loop

    // ── HTTP / resource handling ─────────────────────────────────────────
    void setupServer(int port);
    bool serveResource(const std::string &uri, WsServer::connection_ptr con);
    nlohmann::json listFiles(const std::string &subdir = "");
    std::string resolveDataFile(const std::string &relpath);
    nlohmann::json buildConfig();
    nlohmann::json handleElogPost(const std::string &body);
    void onHttp(WsServer *srv, websocketpp::connection_hdl hdl);

    // ---- Auto-report dispatch (on-demand) ----------------------------------
    // The auto-post pipeline is server-driven: on a trigger event (END or
    // run-number change), the server picks ONE alive WS client, asks it
    // to capture screenshots + POST, then assembles the report.  Only one
    // capture is in flight at a time.  If the chosen client doesn't
    // respond inside ELOG_CAPTURE_TIMEOUT_S, the watchdog forwards the
    // request to the next alive client.  When all candidates are
    // exhausted the failure is recorded in summary.json.
    std::mutex elog_post_mtx_;
    struct PendingCapture {
        std::string request_id;
        uint32_t    run     = 0;
        std::string reason;
        std::time_t started = 0;     // unix seconds, last dispatch attempt
        std::set<websocketpp::connection_hdl,
                 std::owner_less<websocketpp::connection_hdl>> tried;
    };
    std::mutex                   pending_capture_mtx_;
    std::optional<PendingCapture> pending_capture_;
    bool save_dir_writable_ = false;

    // Pick an alive WS client (excluding `tried`), dispatch a
    // capture_request, and update pending_capture_.  Returns false if no
    // candidate available — caller should record the failure.
    bool dispatchCapture(uint32_t run, const std::string &reason,
                         const std::string &request_id_in = "");
    // Watchdog poll — call periodically (piggy-backed on the monitor
    // status thread) to retry stale requests.
    void autoReportWatchdog();
    // Write/refresh <local_save_dir>/summary.json after a save attempt.
    // low_data marks a "suspicious" report (all accumulators were empty
    // at capture time); recorded in the per-run record so operators
    // grep'ing summary.json can spot bad reports without opening the XML.
    void appendAutoReportSummary(uint32_t run, const std::string &saved_xml,
                                 bool posted, const std::string &lognumber,
                                 const std::string &reason,
                                 const std::string &error,
                                 bool low_data = false);
    // Verify we can write to local_save_dir at startup.  Sets
    // save_dir_writable_ + logs result.
    void checkSaveDirWritable();

    // ---- Deferred autoclear (run-boundary path) ----------------------------
    // scheduleAutoClear() — called from etReaderThread on END / PRESTART /
    // run-change — starts a server-side wipe of hist+lms+epics after a
    // delay (5 s at every current call site).  The countdown PAUSES while
    // pending_capture_ is set, so a run-change-fallback capture firing on
    // the first physics of the new run still snapshots the prior run's
    // data before the wipe.  Last-call-wins: repeated calls re-anchor the
    // countdown to the most recent boundary.  Manual /api/hist/clear etc.
    // bypass this entirely so operator presses are immediate.
    struct AutoClearState {
        bool                                  pending      = false;
        int                                   remaining_ms = 0;
        std::chrono::steady_clock::time_point last_tick;
    };
    std::mutex      autoclear_mtx_;
    AutoClearState  autoclear_state_;

    void scheduleAutoClear(int delay_ms);
    void tickAutoClear();
    void runAutoClearNow();
    // Wipe the most-recent monitoring-event GEM APV snapshot.  Called from
    // runAutoClearNow so the GEM APV tab's "Source: Monitoring event" view
    // doesn't show an evt number left over from the previous run.
    void clearLatestFullApv();

    // ---- Auto-report schedule (per-run 45-min timer) ----------------------
    // armScheduleForRun() is idempotent — re-arming for the same run is a
    // no-op; the timer only restarts when the observed run number changes.
    // tickAutoReportSchedule() runs in the same TICK_MS poll as
    // tickAutoClear / autoReportWatchdog; when elapsed >= schedule_minutes
    // it fires a dispatchCapture(run, "schedule").  Per-run on-disk dedup
    // (inside dispatchCapture) drops the dispatch if a report for that run
    // already exists, so the 45-min trigger races run-change / END
    // harmlessly — first-saved wins.
    struct AutoReportSchedule {
        uint32_t                              run = 0;       // current run being timed
        std::chrono::steady_clock::time_point started;       // when first observed
        std::chrono::steady_clock::time_point last_attempt;  // throttles retry-on-failure
        bool                                  armed = false; // ↦ true after armScheduleForRun
        bool                                  fired = false; // ↦ true once dispatch landed (sent or dedup hit)
    };
    std::mutex          ar_sched_mtx_;
    AutoReportSchedule  ar_sched_;

    void armScheduleForRun(uint32_t run);
    void tickAutoReportSchedule();
    // True iff <local_save_dir>/run_NNNNNN/ contains any .xml file.  Used
    // by dispatchCapture for the per-run "one report per run" guarantee.
    bool hasSavedReportForRun(uint32_t run);
};
