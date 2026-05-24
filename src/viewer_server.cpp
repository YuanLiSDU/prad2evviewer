// =========================================================================
// viewer_server.cpp — Unified server implementation
//
// Merges file-based event viewing and online ET monitoring into one server.
// =========================================================================

#include "viewer_server.h"

#ifdef WITH_ET
#include "EtChannel.h"
#endif

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <chrono>

using json = nlohmann::json;
namespace fs = std::filesystem;
using namespace evc;

// =========================================================================
// Progress
// =========================================================================

json Progress::toJson() const
{
    std::lock_guard<std::mutex> lk(mtx);
    return {{"loading", loading.load()},
            {"phase", phase == 1 ? "indexing" : phase == 2 ? "histograms" : "idle"},
            {"current", current.load()}, {"total", total.load()}, {"file", target_file}};
}

void Progress::setFile(const std::string &f)
{
    std::lock_guard<std::mutex> lk(mtx);
    target_file = f;
}


// =========================================================================
// ViewerServer — lifecycle
// =========================================================================

ViewerServer::ViewerServer() = default;

ViewerServer::~ViewerServer()
{
    stop();
    joinAll();
}

void ViewerServer::joinAll()
{
#ifdef WITH_ET
    if (et_thread_.joinable()) et_thread_.join();
    if (monitor_status_thread_.joinable()) monitor_status_thread_.join();
#endif
    { std::lock_guard<std::mutex> lk(load_mtx_);
      if (load_thread_.joinable()) load_thread_.join(); }
    if (server_thread_.joinable()) server_thread_.join();
}

void ViewerServer::init(const Config &cfg)
{
    cfg_ = cfg;
    res_dir_ = cfg.resource_dir;
    hist_enabled_ = cfg.hist_enabled;

    const auto &db_dir = cfg.database_dir;

    // --- resolve DAQ config path ---
    std::string daq_cfg_file = cfg.daq_config_file;
    if (daq_cfg_file.empty())
        daq_cfg_file = findFile("daq_config.json", db_dir);

    // --- load ET config from monitor_config.json "online" section ---
#ifdef WITH_ET
    {
        std::string cf = cfg.monitor_config_file;
        if (cf.empty()) cf = findFile("monitor_config.json", db_dir);
        std::string s = readFile(cf);
        if (!s.empty()) {
            auto j = json::parse(s, nullptr, false);
            if (j.contains("online")) {
                auto &e = j["online"];
                if (e.contains("et_host"))    et_cfg_.host    = e["et_host"];
                if (e.contains("et_port"))    et_cfg_.port    = e["et_port"];
                if (e.contains("et_file"))    et_cfg_.et_file = e["et_file"];
                if (e.contains("et_station")) et_cfg_.station = e["et_station"];
                if (e.contains("ring_buffer_size")) ring_size_ = e["ring_buffer_size"];
            }
        }
    }
#endif

    // --- initialize both AppState instances (same config, separate accumulators) ---
    app_file_.init(db_dir, daq_cfg_file,
                   cfg.monitor_config_file, cfg.reconstruction_config_file);
    app_online_.init(db_dir, daq_cfg_file,
                     cfg.monitor_config_file, cfg.reconstruction_config_file);

    // --- load external filter file if specified ---
    if (!cfg.filter_file.empty()) {
        std::string s = readFile(cfg.filter_file);
        if (!s.empty()) {
            auto fj = json::parse(s, nullptr, false);
            if (!fj.is_discarded()) {
                app_file_.loadFilter(fj);
                app_online_.loadFilter(fj);
            } else {
                std::cerr << "Warning: failed to parse filter file: " << cfg.filter_file << "\n";
            }
        } else {
            std::cerr << "Warning: cannot read filter file: " << cfg.filter_file << "\n";
        }
    }

    // --- build base_config JSON for /api/config ---
    // File pointers come from the parsed daq_cfg (already loaded by AppState::init).
    json hycal_map_j = json::array();
    {
        const std::string map_name = app_file_.daq_cfg.hycal_map_file.empty()
            ? std::string("hycal_map.json") : app_file_.daq_cfg.hycal_map_file;
        std::string s = readFile(findFile(map_name, db_dir));
        if (!s.empty()) hycal_map_j = json::parse(s, nullptr, false);
    }
    json base_cfg = {
        {"hycal_map", hycal_map_j}, {"crate_roc", app_file_.crate_roc_json},
    };
    app_file_.base_config = std::move(base_cfg);
    app_online_.base_config = app_file_.base_config;

    // --- build crate→ROC tag map for ROOT data sources ---
    for (auto &[k, v] : app_file_.crate_roc_json.items())
        crate_to_roc_[std::stoi(k)] = v.get<uint32_t>();

    // --- validate resources ---
    if (readFile(res_dir_ + "/viewer.html").empty())
        std::cerr << "Warning: viewer.html not found in " << res_dir_ << "\n";

    std::cerr << "Database  : " << db_dir << "\n"
              << "Resources : " << res_dir_ << "\n";
    if (!cfg.data_dir.empty())
        std::cerr << "Data dir  : " << cfg.data_dir << "\n";
#ifdef WITH_ET
    std::cerr << "ET target : " << et_cfg_.host << ":" << et_cfg_.port << "\n"
              << "Ring buf  : " << ring_size_ << " events\n";
#endif
}

// =========================================================================
// Server setup & run
// =========================================================================

void ViewerServer::setupServer(int port)
{
    server_ = std::make_unique<WsServer>();
    server_->set_access_channels(websocketpp::log::alevel::none);
    server_->set_error_channels(websocketpp::log::elevel::warn |
                                websocketpp::log::elevel::rerror);
    server_->init_asio();
    server_->set_reuse_addr(true);

    server_->set_http_handler([this](websocketpp::connection_hdl hdl) {
        onHttp(server_.get(), hdl);
    });
    server_->set_open_handler([this](websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(ws_mtx_);
        ws_clients_.insert(hdl);
    });
    server_->set_close_handler([this](websocketpp::connection_hdl hdl) {
        {
            std::lock_guard<std::mutex> lk(ws_mtx_);
            ws_clients_.erase(hdl);
            reporter_capable_.erase(hdl);
        }
        // Always remove from tagger subscribers on disconnect so the gate
        // flips cleanly when the last subscriber drops.
        taggerUnsubscribe(hdl);
    });
    server_->set_message_handler(
        [this](websocketpp::connection_hdl hdl,
               WsServer::message_ptr msg) {
            handleWsMessage(hdl, msg->get_payload());
        });

    // bind — retry a few ports if port==0
    if (port == 0) {
        for (int p = 15050; p < 15150; ++p) {
            try { server_->listen(p); port_ = p; break; }
            catch (...) { continue; }
        }
        if (port_ == 0) { server_->listen(15050); port_ = 15050; }
    } else {
        server_->listen(port);
        port_ = port;
    }

    server_->start_accept();
}

void ViewerServer::run()
{
    setupServer(cfg_.port);
    checkSaveDirWritable();

    // start ET reader thread (sleeps until et_active_)
#ifdef WITH_ET
    et_thread_ = std::thread([this]() { etReaderThread(); });
    // Always spawn the poll thread — it now also drives the auto-report
    // watchdog, which must run even if no shell metrics are configured.
    monitor_status_thread_ = std::thread([this]() { monitorStatusPollThread(); });
#endif

    // load initial file (blocking)
    if (!cfg_.initial_file.empty()) {
        mode_ = Mode::File;
        hist_enabled_ = cfg_.hist_enabled;
        loadFileInternal(cfg_.initial_file);
    }

    // if --et, switch to online mode
    if (cfg_.start_online) {
        mode_ = Mode::Online;
#ifdef WITH_ET
        et_active_ = true;
#endif
    }

    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    std::cout << "Server at http://localhost:" << port_ << "\n"
              << "  Mode: " << mode() << "\n"
              << "  " << (data ? data->event_count : 0) << " events"
              << (hist_enabled_ ? ", histograms enabled" : "")
              << (!cfg_.data_dir.empty() ? ", file browser enabled" : "")
              << (cfg_.interactive ? "\n  Type 'help' for commands, " : "\n  ")
              << "Ctrl+C to stop\n";

    // interactive stdin command loop (detached — will be killed on exit)
    if (cfg_.interactive)
        std::thread([this]() { commandLoop(); }).detach();

    server_->run();

    // cleanup — stop() already set flags; joinAll() waits for threads
    stop();
    joinAll();
}

int ViewerServer::startAsync(int port)
{
    if (port != 0) cfg_.port = port;
    setupServer(cfg_.port);
    checkSaveDirWritable();

    // start ET reader thread
#ifdef WITH_ET
    et_thread_ = std::thread([this]() { etReaderThread(); });
    // Always spawn the poll thread — it now also drives the auto-report
    // watchdog, which must run even if no shell metrics are configured.
    monitor_status_thread_ = std::thread([this]() { monitorStatusPollThread(); });
#endif

    // load initial file (async)
    if (!cfg_.initial_file.empty())
        loadFile(cfg_.initial_file, cfg_.hist_enabled);

    // start online if requested
    if (cfg_.start_online) {
        mode_ = Mode::Online;
#ifdef WITH_ET
        et_active_ = true;
#endif
    }

    // run server in background
    server_thread_ = std::thread([this]() { server_->run(); });
    return port_;
}

void ViewerServer::stop()
{
    running_ = false;
#ifdef WITH_ET
    et_active_ = false;
#endif

    if (server_) {
        try { server_->stop_listening(); server_->stop(); } catch (...) {}
    }

    // Thread joins are handled by run()/startAsync() cleanup — not here.
    // stop() may be called from a signal handler where blocking on join()
    // would hang if a thread is stuck in a library call (e.g. et_open).
}

// =========================================================================
// Mode switching
// =========================================================================

std::string ViewerServer::mode() const
{
    switch (mode_.load()) {
    case Mode::File:   return "file";
    case Mode::Online: return "online";
    default:           return "idle";
    }
}

AppState &ViewerServer::activeApp()
{
    return mode_.load() == Mode::Online ? app_online_ : app_file_;
}

void ViewerServer::setMode(Mode m)
{
    mode_ = m;
    wsBroadcast(json({{"type", "mode_changed"}, {"mode", mode()}}).dump());
}

// =========================================================================
// WebSocket broadcast
// =========================================================================

void ViewerServer::wsBroadcast(const std::string &msg)
{
    std::lock_guard<std::mutex> lk(ws_mtx_);
    for (auto &hdl : ws_clients_) {
        try { server_->send(hdl, msg, websocketpp::frame::opcode::text); }
        catch (...) {}
    }
}

// ── Incoming WebSocket messages (JSON text) ────────────────────────────────

void ViewerServer::handleWsMessage(websocketpp::connection_hdl hdl,
                                   const std::string &payload)
{
    // Ignore empty / binary frames.  We only speak a tiny control vocabulary.
    if (payload.empty() || payload[0] != '{') return;

    // Guard against any json type_error / bad_alloc leaking into asio — an
    // uncaught exception here would tear down the io_context and crash the
    // server.
    try {

    auto j = json::parse(payload, nullptr, false);
    if (j.is_discarded() || !j.is_object() ||
        !j.contains("type") || !j["type"].is_string()) return;
    const std::string t = j["type"].get<std::string>();

    if (t == "client_hello") {
        // The on-demand auto-report rollout: only clients that
        // advertise the "auto_report" capability are eligible to be
        // picked by dispatchCapture.  Old (pre-update) tabs never send
        // this message and stay out of the candidate pool — they keep
        // receiving every other broadcast unchanged.
        bool can_capture = false;
        if (j.contains("capabilities") && j["capabilities"].is_array()) {
            for (auto &c : j["capabilities"]) {
                if (c.is_string() && c.get<std::string>() == "auto_report") {
                    can_capture = true; break;
                }
            }
        }
        if (can_capture) {
            std::lock_guard<std::mutex> lk(ws_mtx_);
            reporter_capable_.insert(hdl);
        }
        try {
            server_->send(hdl,
                json({{"type", "server_hello"},
                      {"capabilities", json::array({"auto_report"})}}).dump(),
                websocketpp::frame::opcode::text);
        } catch (...) {}
    }
    else if (t == "tagger_subscribe") {
        taggerSubscribe(hdl);
        // Acknowledge so the client can confirm its subscription took effect.
        try {
            server_->send(hdl,
                json({{"type", "tagger_subscribed"},
                      {"subscribers", tagger_subs_count_.load()}}).dump(),
                websocketpp::frame::opcode::text);
        } catch (...) {}
    }
    else if (t == "tagger_unsubscribe") {
        taggerUnsubscribe(hdl);
        try {
            server_->send(hdl,
                json({{"type", "tagger_unsubscribed"},
                      {"subscribers", tagger_subs_count_.load()}}).dump(),
                websocketpp::frame::opcode::text);
        } catch (...) {}
    }
    // Unknown types are silently ignored — old clients stay happy.

    } catch (const std::exception &e) {
        std::cerr << "[ws] malformed message dropped: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "[ws] malformed message dropped (unknown)\n";
    }
}

// ── Tagger subscription registry ───────────────────────────────────────────

void ViewerServer::taggerSubscribe(websocketpp::connection_hdl hdl)
{
    std::lock_guard<std::mutex> lk(tagger_subs_mtx_);
    if (tagger_subs_.insert(hdl).second)
        tagger_subs_count_.store(static_cast<int>(tagger_subs_.size()));
}

void ViewerServer::taggerUnsubscribe(websocketpp::connection_hdl hdl)
{
    std::lock_guard<std::mutex> lk(tagger_subs_mtx_);
    if (tagger_subs_.erase(hdl))
        tagger_subs_count_.store(static_cast<int>(tagger_subs_.size()));
}

void ViewerServer::taggerBroadcastBinary(const void *data, size_t nbytes)
{
    std::lock_guard<std::mutex> lk(tagger_subs_mtx_);
    for (auto &hdl : tagger_subs_) {
        try {
            server_->send(hdl, data, nbytes,
                          websocketpp::frame::opcode::binary);
        }
        catch (...) {
            tagger_dropped_frames_.fetch_add(1, std::memory_order_relaxed);
        }
    }
}

// =========================================================================
// File mode — loading
// =========================================================================

void ViewerServer::loadFile(const std::string &path, bool hist)
{
    {
        std::lock_guard<std::mutex> lk(mode_mtx_);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online)
            et_active_ = false;
#endif
        hist_enabled_ = hist;
        setMode(Mode::File);
    }

    std::lock_guard<std::mutex> lk(load_mtx_);
    if (load_thread_.joinable()) load_thread_.join();
    load_thread_ = std::thread([this, path]() { loadFileInternal(path); });
}

void ViewerServer::loadFileInternal(const std::string &filepath)
{
    progress_.loading = true;
    progress_.setFile(filepath);
    progress_.phase = 0; progress_.current = 0; progress_.total = 0;

    auto data = std::make_shared<FileData>();
    data->filepath = filepath;

    // reset on-demand accumulation tracking
    { std::lock_guard<std::mutex> lk(ondemand_mtx_); ondemand_processed_.clear(); }

    std::cerr << "Loading: " << filepath << "\n";

    // create and open the appropriate data source
    progress_.phase = 1;
    auto source = createDataSource(filepath, app_file_.daq_cfg, crate_to_roc_, &app_file_.hycal);
    if (!source) {
        std::cerr << "  Error: unsupported file type\n";
        progress_.loading = false; return;
    }
    std::string err = source->open(filepath);
    if (!err.empty()) {
        std::cerr << "  Error: " << err << "\n";
        progress_.loading = false; return;
    }

    data->event_count = source->eventCount();
    data->caps = source->capabilities();
    std::cerr << "  Indexed " << data->event_count << " events"
              << " (source: " << data->caps.source_type << ")\n";

    // install the new data source before building histograms
    // (buildHistograms reads data_source_; loadFileInternal runs on the load
    // thread, and HTTP threads acquire data_source_mtx_ for random access)
    { std::lock_guard<std::mutex> lk(data_source_mtx_); data_source_ = std::move(source); }

    // always reset accumulators when opening a new file
    app_file_.clearHistograms();
    app_file_.clearLms();
    app_file_.clearEpics();

    // build filtered index if filter is active
    if (app_file_.filterActive())
        buildFilteredIndex();

    if (hist_enabled_) {
        progress_.total = data->event_count;
        buildHistograms();
    }

    { std::lock_guard<std::mutex> lk(file_data_mtx_); file_data_ = data; }

    progress_.loading = false;
    progress_.phase = 0;
    std::cerr << "  Ready\n";
}

void ViewerServer::buildHistograms()
{
    app_file_.clearHistograms();
    app_file_.clearLms();
    app_file_.clearEpics();

    if (!data_source_) return;

    fdec::WaveAnalyzer ana(app_file_.daq_cfg.wave_cfg);
    ana.SetTemplateStore(&app_file_.template_store);
    fdec::WaveResult wres;

    progress_.phase = 2; progress_.current = 0;

    data_source_->iterateAll(
        // physics events (EVIO / ROOT raw)
        [&](int idx, fdec::EventData &event, ssp::SspEventData *ssp) {
            progress_.current = app_file_.events_processed.load() + 1;
            // skip events that don't pass the filter
            if (app_file_.filterActive() && !app_file_.evaluateFilter(event, ssp))
                return;
            // GEM reco BEFORE processEvent — runGemEfficiency (inside
            // processEvent) reads gem_sys.GetHits(d), which only holds the
            // current event's hits after processGemEvent runs.  Online and
            // single-event accumulate paths already use this order; the
            // preprocess loop must too, otherwise efficiency reads stale
            // GEM hits from the previous event.
            if (ssp) app_file_.processGemEvent(*ssp);
            app_file_.processEvent(event, ana, wres);
        },
        // recon events (ROOT recon)
        [&](int idx, const ReconEventData &recon) {
            progress_.current = app_file_.events_processed.load() + 1;
            app_file_.processReconEvent(recon);
        },
        // control events (sync/prestart/go)
        [&](uint32_t unix_time, uint64_t last_ti_ts) {
            if (app_file_.sync_unix == 0)
                app_file_.recordSyncTime(unix_time, last_ti_ts);
        },
        // EPICS events
        [&](const std::string &text, int32_t ev_num, uint64_t ts) {
            app_file_.processEpics(text, app_file_.events_processed.load(), ts);
        },
        // DSC2 scaler bank → measured livetime
        [&](const dsc::DscEventData &dsc) {
            app_file_.processDsc(dsc);
        },
        app_file_.daq_cfg.dsc_scaler.bank_tag
    );

    std::cerr << "  Histograms: " << app_file_.events_processed.load() << " events"
              << ", clusters: " << app_file_.cluster_events_processed.load()
              << ", LMS: " << app_file_.lms_events.load() << "\n";
}

// =========================================================================
// File mode — event decoding
// =========================================================================

std::string ViewerServer::decodeRawEvent(int ev1, fdec::EventData &event,
                                         ssp::SspEventData *ssp_evt)
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    if (!data) return "no file loaded";
    int idx = ev1 - 1;
    if (idx < 0 || idx >= data->event_count) return "event out of range";

    std::lock_guard<std::mutex> lk(data_source_mtx_);
    if (!data_source_) return "no data source";
    return data_source_->decodeEvent(idx, event, ssp_evt);
}

// Central accumulation gate for file mode.
// Mirrors online-mode logic: processEvent + processGemEvent are called once
// per event.  Preprocessed files skip this (buildHistograms already did it).
// Any new accumulation added to processEvent() works automatically.
void ViewerServer::accumulate(int ev1, fdec::EventData &event,
                              ssp::SspEventData *ssp)
{
    if (hist_enabled_) return;  // preprocessed — already done

    std::lock_guard<std::mutex> lk(ondemand_mtx_);
    if (!ondemand_processed_.insert(ev1).second) return;  // already seen

    // skip events that don't pass the filter
    if (app_file_.filterActive() && !app_file_.evaluateFilter(event, ssp))
        return;

    if (ssp) app_file_.processGemEvent(*ssp);

    fdec::WaveAnalyzer ana(app_file_.daq_cfg.wave_cfg);
    ana.SetTemplateStore(&app_file_.template_store);
    fdec::WaveResult wres;
    app_file_.processEvent(event, ana, wres);
}

// =========================================================================
// Filters
// =========================================================================

std::string ViewerServer::applyFilter(const json &fj)
{
    std::string err = app_file_.loadFilter(fj);
    if (!err.empty()) return err;
    app_online_.loadFilter(fj);

    // clear + rebuild
    app_file_.clearHistograms();
    app_file_.clearLms();
    app_online_.clearHistograms();
    app_online_.clearLms();
    { std::lock_guard<std::mutex> lk(ondemand_mtx_); ondemand_processed_.clear(); }

    buildFilteredIndex();

    if (hist_enabled_) {
        std::shared_ptr<FileData> data;
        { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
        if (data) buildHistograms();
    }

    wsBroadcast("{\"type\":\"hist_cleared\"}");
    return "";
}

void ViewerServer::clearFilter()
{
    app_file_.unloadFilter();
    app_online_.unloadFilter();

    app_file_.clearHistograms();
    app_file_.clearLms();
    app_online_.clearHistograms();
    app_online_.clearLms();
    { std::lock_guard<std::mutex> lk(ondemand_mtx_); ondemand_processed_.clear(); }

    filtered_indices_.clear();

    if (hist_enabled_) {
        std::shared_ptr<FileData> data;
        { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
        if (data) buildHistograms();
    }

    wsBroadcast("{\"type\":\"hist_cleared\"}");
}

void ViewerServer::buildFilteredIndex()
{
    filtered_indices_.clear();
    if (!app_file_.filterActive()) return;

    std::lock_guard<std::mutex> lk(data_source_mtx_);
    if (!data_source_) return;

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    auto &ssp_evt = *ssp_ptr;
    int idx = 0;

    progress_.loading = true;
    progress_.setFile("Filtering events");
    progress_.phase = 1;
    progress_.current = 0;
    progress_.total = 0;

    data_source_->iterateAll(
        [&](int /*i*/, fdec::EventData &ev, ssp::SspEventData *ssp) {
            idx++;
            progress_.current = idx;
            if (app_file_.evaluateFilter(ev, ssp))
                filtered_indices_.push_back(idx);  // 1-based
        },
        nullptr, nullptr, nullptr
    );

    progress_.loading = false;
    progress_.phase = 0;
    std::cerr << "Filter: " << filtered_indices_.size() << " / " << idx << " events pass\n";
}

json ViewerServer::decodeEvent(int ev1)
{
    // check if this is a recon source (no per-channel data)
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    if (data && data->caps.source_type == "root_recon") {
        // return minimal event info + empty channels
        ReconEventData recon;
        { std::lock_guard<std::mutex> lk(data_source_mtx_);
          if (!data_source_ || !data_source_->decodeReconEvent(ev1 - 1, recon))
              return {{"error", "decode error"}}; }
        return {{"event", ev1}, {"channels", json::object()},
                {"event_number", recon.event_num},
                {"trigger_bits", recon.trigger_bits},
                {"run_number", recon.run_number}};
    }

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    std::string err = decodeRawEvent(ev1, event, ssp_ptr.get());
    if (!err.empty()) return {{"error", err}};

    accumulate(ev1, event, ssp_ptr.get());

    fdec::WaveAnalyzer ana(app_file_.daq_cfg.wave_cfg);
    ana.SetTemplateStore(&app_file_.template_store);
    fdec::WaveResult wres;
    json result = app_file_.encodeEventJson(event, ev1, ana, wres);

    // Tag the event type so the frontend can label non-Physics samples
    // (Sync / EPICS / control) in the status bar instead of showing them as
    // "0 channels, no trigger".  Default is "physics".
    {
        std::lock_guard<std::mutex> lk(data_source_mtx_);
        if (data_source_) {
            using ET = evc::EventType;
            ET et = data_source_->eventTypeAt(ev1 - 1);
            const char *kind = "physics";
            switch (et) {
                case ET::Physics:  kind = "physics";  break;
                case ET::Sync:     kind = "sync";     break;
                case ET::Epics:    kind = "epics";    break;
                case ET::Prestart: kind = "prestart"; break;
                case ET::Go:       kind = "go";       break;
                case ET::End:      kind = "end";      break;
                case ET::Control:  kind = "control";  break;
                case ET::Unknown:  kind = "unknown";  break;
            }
            result["event_kind"] = kind;
        }
    }
    return result;
}

json ViewerServer::computeClusters(int ev1)
{
    // recon source: return pre-computed clusters
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
    if (data && data->caps.source_type == "root_recon") {
        ReconEventData recon;
        { std::lock_guard<std::mutex> lk(data_source_mtx_);
          if (!data_source_ || !data_source_->decodeReconEvent(ev1 - 1, recon))
              return {{"error", "decode error"}}; }
        return app_file_.encodeReconClustersJson(recon, ev1);
    }

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    std::string err = decodeRawEvent(ev1, event, ssp_ptr.get());
    if (!err.empty()) return {{"error", err}};

    accumulate(ev1, event, ssp_ptr.get());

    fdec::WaveAnalyzer ana(app_file_.daq_cfg.wave_cfg);
    ana.SetTemplateStore(&app_file_.template_store);
    fdec::WaveResult wres;
    return app_file_.computeClustersJson(event, ev1, ana, wres);
}

// =========================================================================
// Interactive command loop (stdin)
// =========================================================================

void ViewerServer::commandLoop()
{
    std::string line;
    while (running_ && std::getline(std::cin, line)) {
        // trim whitespace
        auto s = line.find_first_not_of(" \t");
        if (s == std::string::npos) continue;
        line = line.substr(s, line.find_last_not_of(" \t") - s + 1);
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string cmd;
        iss >> cmd;

        if (cmd == "help" || cmd == "?") {
            std::cout << "Commands:\n"
                      << "  status            — show server status\n"
                      << "  load <file> [H]   — load an evio file (H=1 for histograms)\n"
#ifdef WITH_ET
                      << "  online            — switch to ET/online mode\n"
#endif
                      << "  offline           — switch to file mode\n"
                      << "  clear <what>      — clear hist, lms, or epics\n"
                      << "  filter            — show current filter\n"
                      << "  filter load <f>   — load filter from JSON file\n"
                      << "  filter unload     — remove all filters\n"
                      << "  quit / exit       — stop the server\n";
        }
        else if (cmd == "status") {
            std::shared_ptr<FileData> data;
            { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }
            std::cout << "Mode: " << mode() << "\n";
            if (data)
                std::cout << "File: " << data->filepath << " (" << data->event_count << " events)\n";
#ifdef WITH_ET
            std::cout << "ET: " << (et_connected_.load() ? "connected" : "not connected")
                      << " (" << et_cfg_.host << ":" << et_cfg_.port << ")\n";
            std::cout << "Events processed: " << app_online_.events_processed.load() << "\n";
#endif
        }
        else if (cmd == "load") {
            std::string path;
            int hist = 0;
            iss >> path >> hist;
            if (path.empty()) {
                std::cout << "Usage: load <filepath> [hist=0|1]\n";
            } else {
                loadFile(path, hist != 0);
                std::cout << "Loading " << path << (hist ? " (with histograms)" : "") << "\n";
            }
        }
#ifdef WITH_ET
        else if (cmd == "online") {
            std::lock_guard<std::mutex> lk(mode_mtx_);
            if (mode_.load() == Mode::Online)
                et_generation_++;
            et_active_ = true;
            setMode(Mode::Online);
            std::cout << "Switched to online mode\n";
        }
#endif
        else if (cmd == "offline") {
            std::lock_guard<std::mutex> lk(mode_mtx_);
            if (mode_.load() == Mode::Online) {
#ifdef WITH_ET
                et_active_ = false;
#endif
                std::shared_ptr<FileData> data;
                { std::lock_guard<std::mutex> lk2(file_data_mtx_); data = file_data_; }
                setMode(data ? Mode::File : Mode::Idle);
            }
            std::cout << "Switched to " << mode() << " mode\n";
        }
        else if (cmd == "clear") {
            std::string what;
            iss >> what;
            if (what == "hist") {
                activeApp().clearHistograms();
                wsBroadcast("{\"type\":\"hist_cleared\"}");
                std::cout << "Histograms cleared\n";
            } else if (what == "lms") {
                activeApp().clearLms();
                wsBroadcast("{\"type\":\"lms_cleared\"}");
                std::cout << "LMS cleared\n";
            } else if (what == "epics") {
                activeApp().clearEpics();
                wsBroadcast("{\"type\":\"epics_cleared\"}");
                std::cout << "EPICS cleared\n";
            } else {
                std::cout << "Usage: clear hist|lms|epics\n";
            }
        }
        else if (cmd == "filter") {
            std::string sub;
            iss >> sub;
            if (sub == "load") {
                std::string path;
                iss >> path;
                if (path.empty()) {
                    std::cout << "Usage: filter load <path.json>\n";
                } else {
                    std::string s = readFile(path);
                    if (s.empty()) { std::cout << "Cannot read: " << path << "\n"; }
                    else {
                        auto fj = json::parse(s, nullptr, false);
                        if (fj.is_discarded()) { std::cout << "Invalid JSON\n"; }
                        else {
                            std::string err = applyFilter(fj);
                            if (err.empty())
                                std::cout << "Filter loaded: " << filtered_indices_.size()
                                          << " events pass\n";
                            else std::cout << "Error: " << err << "\n";
                        }
                    }
                }
            } else if (sub == "unload") {
                clearFilter();
                std::cout << "Filters unloaded\n";
            } else {
                std::cout << activeApp().filterToJson().dump(2) << "\n";
            }
        }
        else if (cmd == "quit" || cmd == "exit") {
            std::cout << "Stopping server...\n";
            stop();
            break;
        }
        else {
            std::cout << "Unknown command: " << cmd << " (type 'help')\n";
        }
    }
}

// =========================================================================
// Online mode — ET reader thread
// =========================================================================


