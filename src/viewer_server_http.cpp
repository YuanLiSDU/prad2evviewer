#include "viewer_server.h"
#include "http_compress.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>

namespace fs = std::filesystem;
using json = nlohmann::json;

// =========================================================================
// Resource serving
// =========================================================================

bool ViewerServer::serveResource(const std::string &uri,
                                 WsServer::connection_ptr con)
{
    if (res_dir_.empty()) return false;
    std::string relpath = (uri == "/") ? "viewer.html" : uri.substr(1);
    if (relpath.find("..") != std::string::npos || relpath[0] == '/')
        return false;
    std::string fullpath = res_dir_ + "/" + relpath;
    std::string content = readFile(fullpath);
    if (content.empty()) return false;
    con->set_status(websocketpp::http::status_code::ok);
    con->set_body(content);
    con->append_header("Content-Type", contentType(fullpath));
    return true;
}

json ViewerServer::listFiles(const std::string &subdir)
{
    json entries = json::array();
    if (cfg_.data_dir.empty()) return entries;
    try {
        fs::path root(cfg_.data_dir);
        fs::path dir = subdir.empty() ? root : root / subdir;
        // security: ensure dir is under root
        auto canon_root = fs::canonical(root);
        auto canon_dir  = fs::canonical(dir);
        if (canon_dir.string().rfind(canon_root.string(), 0) != 0)
            return entries;

        for (auto &entry : fs::directory_iterator(
                 canon_dir, fs::directory_options::skip_permission_denied)) {
            auto rel = fs::relative(entry.path(), root).string();
            if (entry.is_directory()) {
                // count data files inside (non-recursive quick scan)
                int count = 0;
                try {
                    for (auto &child : fs::recursive_directory_iterator(
                             entry.path(), fs::directory_options::skip_permission_denied)) {
                        if (!child.is_regular_file()) continue;
                        auto fn = child.path().filename().string();
                        if (fn.find(".evio") != std::string::npos ||
                            fn.find(".root") != std::string::npos)
                            count++;
                    }
                } catch (...) {}
                if (count > 0)
                    entries.push_back(json{{"type", "dir"}, {"name", rel}, {"count", count}});
            } else if (entry.is_regular_file()) {
                auto fn = entry.path().filename().string();
                if (fn.find(".evio") == std::string::npos &&
                    fn.find(".root") == std::string::npos)
                    continue;
                auto sz = entry.file_size();
                entries.push_back(json{{"type", "file"}, {"name", rel},
                                       {"size", sz},
                                       {"size_mb", std::round(sz / 1048576.0 * 10) / 10}});
            }
        }
    } catch (...) {}
    std::sort(entries.begin(), entries.end(),
              [](const json &a, const json &b) {
                  // dirs first, then by name
                  bool da = a["type"] == "dir", db = b["type"] == "dir";
                  if (da != db) return da > db;
                  return a["name"] < b["name"];
              });
    return entries;
}

std::string ViewerServer::resolveDataFile(const std::string &relpath)
{
    if (cfg_.data_dir.empty()) return "";
    try {
        fs::path full = fs::canonical(fs::path(cfg_.data_dir) / relpath);
        fs::path root = fs::canonical(fs::path(cfg_.data_dir));
        if (full.string().rfind(root.string(), 0) != 0) return "";
        if (!fs::is_regular_file(full)) return "";
        return full.string();
    } catch (...) { return ""; }
}

// =========================================================================
// Config JSON
// =========================================================================

json ViewerServer::buildConfig()
{
    std::shared_ptr<FileData> data;
    { std::lock_guard<std::mutex> lk(file_data_mtx_); data = file_data_; }

    auto &app = activeApp();
    json cfg = app.base_config;
    cfg["mode"] = mode();
#ifdef WITH_ET
    cfg["et_available"] = true;
    cfg["et_connected"] = et_connected_.load();
    cfg["ring_buffer_size"] = ring_size_;
    cfg["et_config"] = {
        {"host", et_cfg_.host}, {"port", et_cfg_.port},
        {"et_file", et_cfg_.et_file}, {"station", et_cfg_.station},
    };
#else
    cfg["et_available"] = false;
    cfg["et_connected"] = false;
    cfg["ring_buffer_size"] = 0;
    cfg["et_config"] = json::object();
#endif
    cfg["file_available"] = !cfg_.data_dir.empty() || (data != nullptr);
    cfg["total_events"] = data ? data->event_count : 0;
    cfg["current_file"] = data ? data->filepath : "";
    cfg["data_dir_enabled"] = !cfg_.data_dir.empty();
    cfg["data_dir"] = cfg_.data_dir;
    cfg["hist_enabled"] = (mode_.load() == Mode::Online) ? true : hist_enabled_.load();
    cfg["filter_active"] = app_file_.filterActive();
    cfg["filtered_count"] = filtered_indices_.empty() ? (data ? data->event_count : 0)
                                                       : (int)filtered_indices_.size();

    // data source capabilities
    // In online mode without a file, report EVIO-native capabilities
    DataSourceCaps caps;
    if (data) {
        caps = data->caps;
    } else if (mode_.load() == Mode::Online) {
        caps.source_type   = "evio";
        caps.has_waveforms = true;
        caps.has_peaks     = true;
        caps.has_pedestals = true;
        caps.has_epics     = true;
        caps.has_sync      = true;
    }
    cfg["source"] = {
        {"type", caps.source_type},
        {"has_waveforms", caps.has_waveforms},
        {"has_peaks", caps.has_peaks},
        {"has_pedestals", caps.has_pedestals},
        {"has_clusters", caps.has_clusters},
        {"has_gem_raw", caps.has_gem_raw},
        {"has_gem_hits", caps.has_gem_hits},
        {"has_epics", caps.has_epics},
        {"has_sync", caps.has_sync},
    };

    app.fillConfigJson(cfg);
    return cfg;
}

// =========================================================================
// Elog post
// =========================================================================

namespace {
// Minimal base64 decoder (RFC4648). Skips whitespace and '='.
std::string _b64Decode(const std::string &s)
{
    static int8_t TBL[256];
    static bool init = false;
    if (!init) {
        for (int i = 0; i < 256; ++i) TBL[i] = -1;
        const char *a = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        for (int i = 0; i < 64; ++i) TBL[(unsigned char)a[i]] = static_cast<int8_t>(i);
        init = true;
    }
    std::string out;
    out.reserve(s.size() * 3 / 4);
    int val = 0, valb = -8;
    for (unsigned char c : s) {
        if (TBL[c] < 0) continue;
        val = (val << 6) | TBL[c];
        valb += 6;
        if (valb >= 0) {
            out.push_back(static_cast<char>((val >> valb) & 0xFF));
            valb -= 8;
        }
    }
    return out;
}

// Pull <Attachment> blocks out of the elog XML body, base64-decode each,
// and write the bytes to <dir>/<filename>. The XML is the canonical
// payload — there is no separate JSON attachment field anymore.
void _extractXmlAttachments(const std::string &xml, const fs::path &dir)
{
    static const std::string A_OPEN  = "<Attachment>";
    static const std::string A_CLOSE = "</Attachment>";
    static const std::string F_OPEN  = "<filename>";
    static const std::string F_CLOSE = "</filename>";
    static const std::string D_OPEN  = "<data encoding=\"base64\">";
    static const std::string D_CLOSE = "</data>";
    size_t pos = 0;
    while (true) {
        auto a_beg = xml.find(A_OPEN, pos);
        if (a_beg == std::string::npos) break;
        auto a_end = xml.find(A_CLOSE, a_beg);
        if (a_end == std::string::npos) break;
        const std::string block = xml.substr(a_beg, a_end - a_beg);
        pos = a_end + A_CLOSE.size();

        auto fb = block.find(F_OPEN);
        auto fe = block.find(F_CLOSE);
        auto db = block.find(D_OPEN);
        auto de = block.find(D_CLOSE);
        if (fb == std::string::npos || fe == std::string::npos ||
            db == std::string::npos || de == std::string::npos) continue;
        std::string filename = block.substr(fb + F_OPEN.size(),
                                            fe - fb - F_OPEN.size());
        std::string b64      = block.substr(db + D_OPEN.size(),
                                            de - db - D_OPEN.size());
        // strip any path component (don't trust the client filename)
        auto slash = filename.find_last_of("/\\");
        if (slash != std::string::npos) filename = filename.substr(slash + 1);
        if (filename.empty()) continue;

        std::ofstream f(dir / filename, std::ios::binary);
        if (!f) {
            std::cerr << "Elog local-save: open "
                      << (dir / filename) << " failed\n";
            continue;
        }
        std::string raw = _b64Decode(b64);
        f.write(raw.data(), static_cast<std::streamsize>(raw.size()));
    }
}

// Result of the local archive write — both the directory and the XML
// path so the post step can read back from the same file.
struct LocalSaveResult {
    fs::path    dir;
    fs::path    xml;
    bool        ok = false;
    std::string error;
};

LocalSaveResult _saveReportLocally(const std::string &local_save_dir,
                                   uint32_t run, const std::string &xml_body)
{
    LocalSaveResult r;
    if (local_save_dir.empty()) {
        r.error = "local_save_dir not configured";
        return r;
    }
    char run_name[32];
    std::snprintf(run_name, sizeof(run_name), "run_%06u", run);
    r.dir = fs::path(local_save_dir) / run_name;
    std::error_code ec;
    fs::create_directories(r.dir, ec);
    if (ec) {
        r.error = "mkdir " + r.dir.string() + ": " + ec.message();
        return r;
    }
    std::time_t t = std::time(nullptr);
    char ts[32];
    std::strftime(ts, sizeof(ts), "%Y%m%dT%H%M%SZ", std::gmtime(&t));
    r.xml = r.dir / (std::string("report_") + ts + ".xml");
    {
        std::ofstream f(r.xml);
        if (!f) { r.error = "open " + r.xml.string() + " failed"; return r; }
        f << xml_body;
        if (!f.good()) {
            r.error = "write " + r.xml.string() + " failed";
            return r;
        }
    }
    _extractXmlAttachments(xml_body, r.dir);
    r.ok = true;
    return r;
}
}  // namespace

json ViewerServer::handleElogPost(const std::string &body)
{
    if (body.empty())
        return {{"ok", false}, {"error", "Empty body"}};
    auto &app = activeApp();
    if (app.elog_url.empty())
        return {{"ok", false}, {"error", "No elog URL configured"}};

    auto req = json::parse(body, nullptr, false);
    if (req.is_discarded() || !req.contains("xml"))
        return {{"ok", false}, {"error", "Invalid request"}};

    std::string xml_body  = req["xml"].get<std::string>();
    bool        is_auto   = req.value("auto", false);
    uint32_t    run       = req.value("run_number", 0u);
    std::string request_id = req.value("request_id", std::string());

    // On-demand auto-post validation: when auto:true, request_id MUST
    // be present and match the in-flight pending capture.  Empty or
    // stale ids are rejected so a buggy / malicious client can't
    // bypass the dispatch handshake.
    if (is_auto) {
        std::lock_guard<std::mutex> lk(pending_capture_mtx_);
        if (request_id.empty()) {
            return {{"ok", false}, {"error", "auto post requires request_id"}};
        }
        if (!pending_capture_ ||
            pending_capture_->request_id != request_id) {
            return {{"ok", false}, {"error", "stale request_id"}};
        }
    }

    // Hard precondition: local archive must be configured.  Every elog
    // post is mirrored on disk first; the actual upload reads back from
    // that file so the local copy is the canonical artifact.
    if (app.auto_report_local_save_dir.empty()) {
        return {{"ok", false},
                {"error", "auto_report.local_save_dir is required"}};
    }
    if (!save_dir_writable_) {
        return {{"ok", false},
                {"error", "local_save_dir not writable"}};
    }

    // Hold elog_post_mtx_ for the entire dedup-check + save sequence so
    // two concurrent POSTs (multi-browser race) can't both pass the
    // on-disk check before either has written. Releases as soon as the
    // file is on disk; the curl upload is allowed to run unlocked.
    std::unique_lock<std::mutex> post_lock(elog_post_mtx_);

    // Server-side dedup: if run_NNNNNN/ already holds an XML modified
    // within the rate window, skip this duplicate. Catches multi-browser
    // races and quick double-fires from the same client (END + run-change
    // fallback) without burning disk on near-identical snapshots.
    if (run > 0 && app.auto_report_min_interval_ms > 0) {
        char run_name[32];
        std::snprintf(run_name, sizeof(run_name), "run_%06u", run);
        fs::path run_dir = fs::path(app.auto_report_local_save_dir) / run_name;
        if (fs::exists(run_dir) && fs::is_directory(run_dir)) {
            std::time_t now_t = std::time(nullptr);
            int window_sec = app.auto_report_min_interval_ms / 1000;
            fs::path most_recent;
            std::time_t most_recent_t = 0;
            std::error_code ec;
            for (auto &entry : fs::directory_iterator(run_dir, ec)) {
                if (!entry.is_regular_file()) continue;
                if (entry.path().extension() != ".xml") continue;
                struct stat st;
                if (::stat(entry.path().c_str(), &st) == 0 &&
                    st.st_mtime > most_recent_t)
                {
                    most_recent_t = st.st_mtime;
                    most_recent   = entry.path();
                }
            }
            if (most_recent_t && (now_t - most_recent_t) < window_sec) {
                std::cerr << "Elog post: server-side dup skip for run "
                          << run << " (recent: " << most_recent
                          << ", " << (now_t - most_recent_t) << "s old)\n";
                return {{"ok", true}, {"skipped", true},
                        {"detail", "server-side dedup: recent save within "
                                  + std::to_string(window_sec / 60) + " min"},
                        {"saved_dir", run_dir.string()},
                        {"saved_xml", most_recent.string()}};
            }
        }
    }

    auto sr = _saveReportLocally(app.auto_report_local_save_dir, run, xml_body);
    if (!sr.ok) {
        std::cerr << "Elog local-save: " << sr.error << "\n";
        return {{"ok", false},
                {"error", "local save failed: " + sr.error}};
    }
    const std::string saved_dir = sr.dir.string();
    const std::string saved_xml = sr.xml.string();

    // Dry-run gate: keep the local copy, skip the upload.
    if (is_auto && !app.auto_report_post_to_elog) {
        std::cerr << "Elog post: auto-mode dry run (post_to_elog=false)"
                  << " saved=" << saved_xml << "\n";
        return {{"ok", true}, {"posted", false}, {"dry_run", true},
                {"saved_dir", saved_dir}, {"saved_xml", saved_xml},
                {"status", "dry_run"}};
    }

    // Save is committed; release the lock so concurrent runs aren't
    // blocked on the curl upload below.
    post_lock.unlock();

    // Upload the saved XML — single source of truth.  The file on disk is
    // byte-identical to what hits elog, so a manual replay can use the
    // same `curl --upload-file <saved_xml>` command.
    std::string cert_flag;
    if (!app.elog_cert.empty())
        cert_flag = " --cert '" + app.elog_cert + "' --key '" + app.elog_key + "'";
    // JLab elog uses /incoming/<name> as a one-shot key — each filename
    // can be processed at most once.  Embedding the run number makes
    // the filename itself a per-run dedup guard: whichever monitor
    // instance / replay tool posts first for a run wins, every other
    // attempt for the same run hits "already processed" at the elog
    // and is rejected without creating a duplicate entry.  Run 0
    // (manual posts, missing run info) falls back to the saved file's
    // timestamped basename.
    std::string upload_name;
    if (run > 0) {
        char buf[40];
        std::snprintf(buf, sizeof(buf), "prad2_run_%06u.xml", run);
        upload_name = buf;
    } else {
        upload_name = "prad2_" + sr.xml.filename().string();
    }
    std::string cmd = "curl -s -o /dev/null -w '%{http_code}'" + cert_flag
                    + " --upload-file '" + saved_xml + "' '"
                    + app.elog_url + "/incoming/" + upload_name
                    + "' 2>/dev/null";
    std::string http_code;
    FILE *p = popen(cmd.c_str(), "r");
    if (p) {
        char buf[256] = {};
        if (fgets(buf, sizeof(buf), p)) http_code = buf;
        while (!http_code.empty() &&
               (http_code.back() == '\n' || http_code.back() == '\r'))
            http_code.pop_back();
        pclose(p);
    }
    bool ok = (http_code.find("200") != std::string::npos ||
               http_code.find("201") != std::string::npos);
    std::cerr << "Elog post: " << app.elog_url << " <- " << saved_xml
              << " -> HTTP " << http_code << (ok ? " OK" : " FAIL") << "\n";

    // On-demand auto-post handshake: clear the pending slot and tell
    // every connected client to drop their "Auto-reporting…" badge.
    if (is_auto && !request_id.empty()) {
        std::lock_guard<std::mutex> lk(pending_capture_mtx_);
        if (pending_capture_ &&
            pending_capture_->request_id == request_id) {
            pending_capture_.reset();
        }
        json done = {{"type", "auto_capture_done"},
                     {"run", run}, {"posted", ok},
                     {"saved_xml", saved_xml}};
        wsBroadcast(done.dump());
        appendAutoReportSummary(run, saved_xml, ok, std::string(),
                                std::string("upload"), std::string());
    }
    return {{"ok", ok}, {"posted", ok}, {"status", http_code},
            {"saved_dir", saved_dir}, {"saved_xml", saved_xml}};
}

// =========================================================================
// On-demand auto-report dispatch
// =========================================================================

namespace {
// Compose summary.json next to local_save_dir.  Whole-file rewrite per
// update — the file stays small (one JSON record per run), and atomic
// O_TRUNC is fine for our cadence.
fs::path _summaryPath(const std::string &local_save_dir)
{
    return fs::path(local_save_dir) / "summary.json";
}

std::string _isoNowZ()
{
    std::time_t t = std::time(nullptr);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
    return std::string(buf);
}
}  // namespace

void ViewerServer::checkSaveDirWritable()
{
    auto &app = activeApp();
    save_dir_writable_ = false;
    if (app.auto_report_local_save_dir.empty()) {
        std::cerr << "AutoReport: local_save_dir not configured — disabled\n";
        return;
    }
    fs::path dir = app.auto_report_local_save_dir;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        std::cerr << "AutoReport: cannot create " << dir
                  << " (" << ec.message() << ") — disabled\n";
        return;
    }
    fs::path probe = dir / ".prad2_write_check";
    {
        std::ofstream f(probe);
        if (!f) {
            std::cerr << "AutoReport: " << dir
                      << " is not writable — disabled\n";
            return;
        }
        f << "ok\n";
    }
    fs::remove(probe, ec);
    save_dir_writable_ = true;
    std::cerr << "AutoReport: save_dir " << dir << " is writable\n";
}

void ViewerServer::appendAutoReportSummary(uint32_t run,
                                           const std::string &saved_xml,
                                           bool posted,
                                           const std::string &lognumber,
                                           const std::string &reason,
                                           const std::string &error)
{
    auto &app = activeApp();
    if (app.auto_report_local_save_dir.empty()) return;
    fs::path path = _summaryPath(app.auto_report_local_save_dir);
    json doc = json::object();
    {
        std::ifstream f(path);
        if (f) {
            std::stringstream ss; ss << f.rdbuf();
            auto parsed = json::parse(ss.str(), nullptr, false);
            if (!parsed.is_discarded() && parsed.is_object()) doc = parsed;
        }
    }
    doc["updated"]            = _isoNowZ();
    doc["auto_post_enabled"]  = app.auto_report_enabled;
    doc["post_to_elog"]       = app.auto_report_post_to_elog;
    doc["save_dir_writable"]  = save_dir_writable_;
    doc["min_interval_ms"]    = app.auto_report_min_interval_ms;
    if (!doc.contains("runs") || !doc["runs"].is_array())
        doc["runs"] = json::array();

    json rec = {
        {"run",       run},
        {"saved_xml", saved_xml},
        {"saved_at",  _isoNowZ()},
        {"posted",    posted},
        {"lognumber", lognumber},
        {"reason",    reason},
    };
    if (!error.empty()) rec["error"] = error;
    doc["runs"].push_back(rec);

    std::ofstream f(path, std::ios::trunc);
    if (!f) {
        std::cerr << "AutoReport summary: cannot write " << path << "\n";
        return;
    }
    f << doc.dump(2) << "\n";
}

bool ViewerServer::dispatchCapture(uint32_t run, const std::string &reason,
                                   const std::string &request_id_in)
{
    auto &app = activeApp();
    if (!app.auto_report_enabled) return false;
    if (!save_dir_writable_) {
        std::cerr << "AutoReport: dispatch skipped (save_dir not writable)\n";
        appendAutoReportSummary(run, std::string(), false, std::string(),
                                reason, "save_dir not writable");
        return false;
    }
    if (!run) return false;

    // Per-run dedup: one report per run.  If any XML for this run is
    // already saved under local_save_dir, drop the dispatch silently —
    // the first trigger to land on disk (schedule / run-change / END)
    // wins.  Watchdog retries (request_id_in non-empty) bypass this so
    // an in-flight retry can still complete.
    if (request_id_in.empty() && hasSavedReportForRun(run)) {
        std::cerr << "AutoReport: dispatch skipped for run " << run
                  << " (" << reason << ") — report already saved\n";
        return false;
    }

    // Only clients that advertised the on-demand auto-report
    // capability via client_hello are eligible.  Old (pre-update) tabs
    // never appear here so the watchdog won't waste 30 s timing them
    // out — operators see a clean "no responsive client" record
    // instead of a long stutter.
    std::set<websocketpp::connection_hdl,
             std::owner_less<websocketpp::connection_hdl>> alive;
    size_t total_clients = 0;
    {
        std::lock_guard<std::mutex> lk(ws_mtx_);
        total_clients = ws_clients_.size();
        alive = reporter_capable_;
    }
    if (alive.empty()) {
        std::string detail = total_clients
            ? "no on-demand-capable client (" +
              std::to_string(total_clients) +
              " legacy tab(s) connected — refresh required)"
            : "no connected client";
        std::cerr << "AutoReport: " << detail << " for run "
                  << run << " (" << reason << ")\n";
        appendAutoReportSummary(run, std::string(), false, std::string(),
                                reason, detail);
        return false;
    }

    std::lock_guard<std::mutex> plk(pending_capture_mtx_);

    constexpr int CAPTURE_TIMEOUT_S = 30;
    std::time_t now_t = std::time(nullptr);

    PendingCapture pc;
    if (!request_id_in.empty() && pending_capture_ &&
        pending_capture_->request_id == request_id_in) {
        // watchdog retry — same request, advance to next candidate
        pc = *pending_capture_;
    } else if (pending_capture_ &&
               (now_t - pending_capture_->started) < CAPTURE_TIMEOUT_S) {
        // Another capture is already in flight (different request_id);
        // don't clobber it. Per-run on-disk dedup will absorb any
        // genuine duplicate when both finish.
        std::cerr << "AutoReport: dispatch skipped for run " << run
                  << " — capture for run " << pending_capture_->run
                  << " (request " << pending_capture_->request_id
                  << ") still in flight\n";
        return false;
    } else {
        char buf[40];
        std::snprintf(buf, sizeof(buf), "%lx-%u", (long)now_t, run);
        pc.request_id = buf;
        pc.run        = run;
        pc.reason     = reason;
        pc.tried.clear();
    }
    pc.started = now_t;

    // First alive client not in pc.tried.
    websocketpp::connection_hdl chosen;
    bool found = false;
    for (auto &h : alive) {
        if (pc.tried.find(h) == pc.tried.end()) {
            chosen = h;
            found  = true;
            break;
        }
    }
    if (!found) {
        std::cerr << "AutoReport: all clients exhausted for run " << run
                  << " (request " << pc.request_id << ")\n";
        appendAutoReportSummary(run, std::string(), false, std::string(),
                                reason, "all clients timed out");
        pending_capture_.reset();
        return false;
    }
    pc.tried.insert(chosen);
    pending_capture_ = pc;

    json msg = {
        {"type",       "capture_request"},
        {"request_id", pc.request_id},
        {"run",        pc.run},
        {"reason",     pc.reason},
    };
    try {
        server_->send(chosen, msg.dump(),
                      websocketpp::frame::opcode::text);
    } catch (const std::exception &e) {
        std::cerr << "AutoReport: send failed (" << e.what()
                  << "), retrying via watchdog\n";
        return false;
    }
    std::cerr << "AutoReport: dispatched capture_request for run " << run
              << " (reason " << pc.reason
              << ", request " << pc.request_id
              << ", attempt " << pc.tried.size() << ")\n";
    return true;
}

void ViewerServer::autoReportWatchdog()
{
    constexpr int CAPTURE_TIMEOUT_S = 30;
    std::optional<PendingCapture> stale;
    {
        std::lock_guard<std::mutex> lk(pending_capture_mtx_);
        if (pending_capture_ &&
            std::time(nullptr) - pending_capture_->started > CAPTURE_TIMEOUT_S)
        {
            stale = pending_capture_;
        }
    }
    if (!stale) return;
    std::cerr << "AutoReport: watchdog timeout for run " << stale->run
              << " (request " << stale->request_id << ")\n";
    dispatchCapture(stale->run, stale->reason, stale->request_id);
}

// =========================================================================
// Deferred autoclear — gates the PRESTART data wipe on capture state
// =========================================================================
// Why: with multiple browser tabs connected, every non-chosen tab calling
// /api/hist/clear on PRESTART wipes the server's histograms while the
// chosen reporter is still snapping screenshots, and the run-change
// fallback (when END is missed) dispatches its capture AFTER PRESTART has
// already broadcast.  Both races land empty data + a low Samples count
// in the auto-report.  Funnelling PRESTART clears through this delayed,
// capture-aware scheduler keeps the server's state intact long enough
// for either the END-driven or run-change-driven capture to complete.
// Manual /api/hist/clear is intentionally NOT routed here — operator
// presses must take effect immediately.
// =========================================================================

void ViewerServer::scheduleAutoClear(int delay_ms)
{
    if (delay_ms < 0) delay_ms = 5000;
    std::lock_guard<std::mutex> lk(autoclear_mtx_);
    // Last-call-wins.  Run-boundary signals fire in sequence (e.g. END
    // schedules 5 s, then PRESTART arrives and re-schedules 5 s); each
    // refresh anchors the countdown to the most recent boundary.  Tick
    // pauses while pending_capture_ is set, so the actual fire time is
    // (last schedule) + delay + (paused capture duration).
    autoclear_state_.pending      = true;
    autoclear_state_.remaining_ms = delay_ms;
    autoclear_state_.last_tick    = std::chrono::steady_clock::now();
    std::cerr << "AutoClear: scheduled in " << delay_ms << " ms\n";
}

void ViewerServer::tickAutoClear()
{
    bool should_run = false;
    {
        std::lock_guard<std::mutex> lk(autoclear_mtx_);
        if (!autoclear_state_.pending) return;

        auto now = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now - autoclear_state_.last_tick).count();
        autoclear_state_.last_tick = now;

        // Pause the countdown while a capture is in flight.  The chosen
        // reporter's screenshots fetch /api/occupancy etc. mid-capture; if
        // the timer expired and we cleared, the run-change-driven capture
        // would land empty data.
        bool capture_in_flight;
        {
            std::lock_guard<std::mutex> lkc(pending_capture_mtx_);
            capture_in_flight = pending_capture_.has_value();
        }
        if (capture_in_flight) return;

        autoclear_state_.remaining_ms -= static_cast<int>(elapsed_ms);
        if (autoclear_state_.remaining_ms <= 0) {
            autoclear_state_.pending      = false;
            autoclear_state_.remaining_ms = 0;
            should_run = true;
        }
    }
    if (should_run) runAutoClearNow();
}

void ViewerServer::runAutoClearNow()
{
    std::cerr << "AutoClear: firing — clearing histograms / lms / epics\n";
    auto &app = activeApp();
    app.clearHistograms();
    app.clearLms();
    app.clearEpics();
    // Per-domain broadcasts keep existing handlers (which do partial UI
    // resets) in step.  The autoclear_done broadcast piggy-backs on top
    // so clients can run a single full clearFrontend in one place
    // instead of inferring it from the three partial resets.
    wsBroadcast("{\"type\":\"hist_cleared\"}");
    wsBroadcast("{\"type\":\"lms_cleared\"}");
    wsBroadcast("{\"type\":\"epics_cleared\"}");
    wsBroadcast("{\"type\":\"autoclear_done\"}");
}

// =========================================================================
// Auto-report schedule — per-run 45-min checkpoint
// =========================================================================
// Two-trigger model (plus END as a last-resort): whichever fires first
// dispatches the capture for that run; subsequent triggers are dropped by
// the on-disk per-run dedup inside dispatchCapture.
//   1) armScheduleForRun(run) — call whenever a control event or physics
//      event surfaces a run_number.  Idempotent on the same run; restarts
//      the timer when the run changes.
//   2) tickAutoReportSchedule() — polled from the same TICK_MS loop as
//      tickAutoClear / autoReportWatchdog.  Fires dispatchCapture once
//      elapsed >= schedule_minutes; flips ar_sched_.fired so the same run
//      isn't re-dispatched on subsequent ticks.
// =========================================================================

void ViewerServer::armScheduleForRun(uint32_t run)
{
    if (run == 0) return;
    std::lock_guard<std::mutex> lk(ar_sched_mtx_);
    if (ar_sched_.armed && ar_sched_.run == run) return;
    auto now = std::chrono::steady_clock::now();
    ar_sched_.run          = run;
    ar_sched_.started      = now;
    // Park last_attempt at epoch so the first post-elapsed tick fires
    // immediately rather than waiting another retry-throttle window.
    ar_sched_.last_attempt = std::chrono::steady_clock::time_point{};
    ar_sched_.armed        = true;
    ar_sched_.fired        = false;
    std::cerr << "AutoReport: schedule armed for run " << run << "\n";
}

void ViewerServer::tickAutoReportSchedule()
{
    constexpr auto RETRY_THROTTLE = std::chrono::seconds(5);

    uint32_t run_to_fire = 0;
    {
        std::lock_guard<std::mutex> lk(ar_sched_mtx_);
        if (!ar_sched_.armed || ar_sched_.fired || ar_sched_.run == 0)
            return;
        int sched_min = activeApp().auto_report_schedule_minutes;
        if (sched_min <= 0) return;
        auto now = std::chrono::steady_clock::now();
        if (now - ar_sched_.started < std::chrono::minutes(sched_min))
            return;
        // Retry-on-failure throttle: don't hammer dispatchCapture every
        // TICK_MS when there's no capable client.  5 s is plenty fast
        // to catch a tab that just connected, and quiet enough to keep
        // the log readable.
        if (now - ar_sched_.last_attempt < RETRY_THROTTLE) return;
        run_to_fire           = ar_sched_.run;
        ar_sched_.last_attempt = now;
    }
    std::cerr << "AutoReport: schedule timer fired for run "
              << run_to_fire << "\n";
    bool dispatched = dispatchCapture(run_to_fire, "schedule");
    // fired ⇒ "we're done with this run".  dispatched=true means a
    // capture_request is on the wire (handleElogPost will reset
    // pending_capture_ when the POST lands).  Even if dispatched
    // returned false, hasSavedReportForRun catches the "dedup'd
    // because END or run-change beat us" case so we don't loop on
    // a run that's already been reported.  Otherwise (no client +
    // no save), leave fired=false so the next throttled tick retries.
    bool done = dispatched || hasSavedReportForRun(run_to_fire);
    if (done) {
        std::lock_guard<std::mutex> lk(ar_sched_mtx_);
        // Guard against a concurrent armScheduleForRun(next_run)
        // having flipped run/started — only mark fired if it's still
        // the same run we just acted on.
        if (ar_sched_.run == run_to_fire) ar_sched_.fired = true;
    }
}

bool ViewerServer::hasSavedReportForRun(uint32_t run)
{
    auto &app = activeApp();
    if (app.auto_report_local_save_dir.empty() || run == 0) return false;
    char run_name[32];
    std::snprintf(run_name, sizeof(run_name), "run_%06u", run);
    fs::path run_dir = fs::path(app.auto_report_local_save_dir) / run_name;
    std::error_code ec;
    if (!fs::exists(run_dir, ec) || !fs::is_directory(run_dir, ec))
        return false;
    // Manual iterator drive — the range-for syntax uses the throwing
    // overload of operator++, which would propagate (e.g. on a
    // concurrent rmdir or a permission flap) all the way out through
    // dispatchCapture and tickAutoReportSchedule and kill the monitor
    // thread.  Drive via the (ec)-overload increment instead, and treat
    // any unexpected error as "no saved report" — the worst case is
    // one spurious dispatch, vs. losing all auto-clears / watchdogs.
    try {
        fs::directory_iterator it(run_dir, ec), end;
        if (ec) return false;
        for (; it != end; it.increment(ec)) {
            if (ec) return false;
            if (it->is_regular_file(ec) && !ec &&
                it->path().extension() == ".xml")
                return true;
        }
    } catch (const std::exception &e) {
        std::cerr << "AutoReport: hasSavedReportForRun(" << run
                  << ") threw — treating as no save: " << e.what() << "\n";
        return false;
    }
    return false;
}

// =========================================================================
// HTTP handler
// =========================================================================

void ViewerServer::onHttp(WsServer *srv, websocketpp::connection_hdl hdl)
{
    auto con = srv->get_con_from_hdl(hdl);
    std::string uri = con->get_resource();

    // gzip if the client advertised it (browsers always do; the urllib-
    // based Python tools don't and keep getting plain bytes).  Read once
    // per request — cheap header lookup, but no point doing it twice.
    bool wants_gzip = prad2::client_accepts_gzip(
        con->get_request_header("Accept-Encoding"));

    // reply(body, [content_type], [pre_gz]) — when pre_gz is non-null and
    // the client accepts gzip, the cached compressed bytes are served
    // verbatim (saves re-deflating the same payload for each viewer).
    // Otherwise, large bodies are compressed on demand if the client
    // accepts gzip; small bodies and gzip-disabled clients get plain.
    auto reply = [&](const std::string &body,
                     const std::string &ct = "application/json",
                     const std::string *pre_gz = nullptr) {
        con->set_status(websocketpp::http::status_code::ok);
        con->append_header("Content-Type", ct);
        if (wants_gzip && pre_gz && !pre_gz->empty()) {
            con->set_body(*pre_gz);
            con->append_header("Content-Encoding", "gzip");
            return;
        }
        if (wants_gzip && body.size() >= prad2::kGzipMinBytes) {
            try {
                con->set_body(prad2::gzip_compress(body));
                con->append_header("Content-Encoding", "gzip");
                return;
            } catch (...) {
                // Fall through to plain body on any zlib failure.
            }
        }
        con->set_body(body);
    };

    // Any exception thrown below (e.g. std::stoi on a malformed %xx in the
    // URL, json type_error, std::bad_alloc) would otherwise unwind into
    // websocketpp / asio and terminate the server.  Convert to a 400/500
    // response and keep the io_context alive.
    try {

    // --- static resources ---
    if (serveResource(uri, con)) return;

    // --- config ---
    if (uri == "/api/config") { reply(buildConfig().dump()); return; }

    // --- runtime peak_filter updates (Waveform-Tab cut settings) ---
    // POST {waveform_filter?, waveform_filter_active?}.  Updates both file &
    // online AppStates, clears accumulated histograms (so the new filter
    // takes effect on the next refresh in online mode), and broadcasts
    // hist_config_updated + hist_cleared.  Peak detection thresholds live in
    // daq_config.json (fadc250_waveform.analyzer) and are not runtime-mutable.
    if (uri == "/api/hist_config") {
        std::string body = con->get_request_body();
        auto j = json::parse(body, nullptr, false);
        if (j.is_discarded() || !j.is_object()) {
            reply("{\"error\":\"invalid JSON\"}"); return;
        }
        bool changed = false;
        auto applyTo = [&](AppState &app) {
            if (j.contains("waveform_filter") && j["waveform_filter"].is_object()) {
                app.peak_filter.parse(j["waveform_filter"], app.peak_quality_bits_def);
                changed = true;
            }
            if (j.contains("waveform_filter_active")
                && j["waveform_filter_active"].is_boolean()) {
                app.peak_filter.enable = j["waveform_filter_active"].get<bool>();
                changed = true;
            }
        };
        applyTo(app_file_);
        applyTo(app_online_);
        if (changed) {
            // Drop accumulated bins so future fills reflect the new filter
            // (online mode would otherwise show old + new entries mixed).
            app_file_.clearHistograms();
            app_online_.clearHistograms();
            wsBroadcast("{\"type\":\"hist_cleared\"}");
        }
        json payload = {
            {"type",                   "hist_config_updated"},
            {"waveform_filter",        app_file_.peak_filter.toJson(app_file_.peak_quality_bits_def)},
            {"waveform_filter_active", app_file_.peak_filter.enable}
        };
        wsBroadcast(payload.dump());
        json resp = payload;
        resp["ok"] = true;
        resp.erase("type");
        reply(resp.dump());
        return;
    }

    // --- mode switching ---
    if (uri == "/api/mode/online") {
#ifdef WITH_ET
        {
            std::lock_guard<std::mutex> lk(mode_mtx_);
            // apply optional ET config overrides (serialised by mode_mtx_)
            std::string body = con->get_request_body();
            if (!body.empty()) {
                auto j = json::parse(body, nullptr, false);
                if (!j.is_discarded()) {
                    if (j.contains("host"))    et_cfg_.host    = j["host"];
                    if (j.contains("port"))    et_cfg_.port    = j["port"];
                    if (j.contains("et_file")) et_cfg_.et_file = j["et_file"];
                    if (j.contains("station")) et_cfg_.station = j["station"];
                }
            }
            // bump generation so ET reader reconnects with new config
            if (mode_.load() == Mode::Online)
                et_generation_++;
            et_active_ = true;
            setMode(Mode::Online);
        }
        reply(json({{"mode", "online"}}).dump());
#else
        reply(json({{"error", "ET support not compiled"}}).dump());
#endif
        return;
    }
    if (uri == "/api/mode/file") {
        std::lock_guard<std::mutex> lk(mode_mtx_);
        if (mode_.load() == Mode::Online) {
#ifdef WITH_ET
            et_active_ = false;
#endif
            std::shared_ptr<FileData> data;
            { std::lock_guard<std::mutex> lk2(file_data_mtx_); data = file_data_; }
            setMode(data ? Mode::File : Mode::Idle);
        }
        reply(json({{"mode", mode()}}).dump());
        return;
    }

    // --- event/latest (online only) ---
    if (uri == "/api/event/latest") {
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            if (!ring_.empty()) { reply(ring_.back().json_str); return; }
            reply("{\"error\":\"no events yet\"}"); return;
        }
#endif
        reply("{\"error\":\"not in online mode\"}"); return;
    }

    // --- event/<n> (mode-dependent) ---
    if (uri.rfind("/api/event/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 11);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &e : ring_) {
                if (e.seq == evnum) { reply(e.json_str); return; }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        reply(decodeEvent(evnum).dump()); return;
    }

    // --- waveform/<n>/<key> (file mode only — on-demand single-channel samples) ---
    if (uri.rfind("/api/waveform/", 0) == 0) {
        // parse /api/waveform/<evnum>/<roc_slot_ch>
        std::string rest = uri.substr(14);
        auto slash = rest.find('/');
        if (slash == std::string::npos) {
            reply("{\"error\":\"usage: /api/waveform/<event>/<roc_slot_ch>\"}"); return;
        }
        int evnum = std::atoi(rest.substr(0, slash).c_str());
        std::string chan_key = rest.substr(slash + 1);

        auto event_ptr = std::make_unique<fdec::EventData>();
        auto ssp_ptr = std::make_unique<ssp::SspEventData>();
        std::string err = decodeRawEvent(evnum, *event_ptr, ssp_ptr.get());
        if (!err.empty()) { reply(json({{"error", err}}).dump()); return; }

        accumulate(evnum, *event_ptr, ssp_ptr.get());

        fdec::WaveAnalyzer ana(activeApp().daq_cfg.wave_cfg);
        ana.SetTemplateStore(&activeApp().template_store);
        fdec::WaveResult wres;
        reply(activeApp().encodeWaveformJson(*event_ptr, chan_key, ana, wres).dump());
        return;
    }

    // --- clusters/<n> (mode-dependent) ---
    if (uri.rfind("/api/clusters/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 14);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &e : ring_) {
                if (e.seq == evnum) {
                    if (e.cluster_str.empty()) {
                        // Cache was invalidated by /api/hist_config; we can't
                        // recompute here without the raw EventData, so report
                        // pending and let the next live event refill.
                        reply("{\"error\":\"clusters pending — config changed, "
                              "wait for next event\"}");
                    } else {
                        reply(e.cluster_str);
                    }
                    return;
                }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        reply(computeClusters(evnum).dump()); return;
    }

    // --- gem/calib (one-shot per-APV pedestal noise + global zs_sigma) ---
    // The frontend caches this and refetches only when the calib_rev
    // embedded in /api/gem/apv/<n> diverges from the cached value.
    if (uri == "/api/gem/calib") {
        reply(activeApp().apiGemCalib().dump());
        return;
    }

    // --- gem/threshold (POST {zs_sigma:N} — applies to all consumers) ---
    // Updates both file & online AppStates so the new threshold takes
    // effect for the active mode immediately and for the other mode on
    // the next event it processes.  Old ring entries keep their encoded
    // hits[] (cosmetic lag of ~ring_size events); each new event reflects
    // the new threshold via its zs_sigma field.  Broadcasts a WS notice
    // so other open viewers can refresh their toolbar input.
    if (uri == "/api/gem/threshold") {
        std::string body = con->get_request_body();
        auto j = json::parse(body, nullptr, false);
        if (j.is_discarded() || !j.is_object() ||
            !j.contains("zs_sigma") || !j["zs_sigma"].is_number()) {
            reply("{\"error\":\"expected {\\\"zs_sigma\\\":N}\"}"); return;
        }
        float new_sigma = j["zs_sigma"].get<float>();
        app_file_.setGemZsSigma(new_sigma);
        app_online_.setGemZsSigma(new_sigma);
        wsBroadcast(json({{"type", "gem_threshold_updated"},
                          {"zs_sigma", new_sigma}}).dump());
        reply(json({{"ok", true},
                    {"zs_sigma", new_sigma}}).dump());
        return;
    }

    // --- gem/apv/latest_full (most recent monitoring event snapshot) ---
    // Holds the JSON of the last event where any APV came in firmware
    // full-readout (nstrips==128).  Encoded by the ET reader with
    // skip_sw_zs=true so the entire pedestal/noise spectrum is visible
    // across all 128 channels regardless of the software σ cut.  Updated
    // on the prescaled DAQ monitoring events; clients learn about updates
    // through the gem_apv_full_event WS broadcast.  Returns an error blob
    // when no full-readout event has been captured yet (typical right
    // after startup or in file mode, where the ET reader doesn't run).
    if (uri == "/api/gem/apv/latest_full") {
#ifdef WITH_ET
        std::string body, gz;
        {
            std::lock_guard<std::mutex> lk(latest_full_apv_mtx_);
            body = latest_full_apv_json_;
            gz   = latest_full_apv_gz_;
        }
        if (body.empty()) {
            reply("{\"error\":\"no full-readout event yet\"}");
        } else {
            reply(body, "application/json", gz.empty() ? nullptr : &gz);
        }
#else
        reply("{\"error\":\"online mode disabled\"}");
#endif
        return;
    }

    // --- gem/apv/<n> (per-event GEM APV waveforms, mode-dependent) ---
    // Online: served from a per-ring-entry pre-encoded string so older
    // events don't disturb the live gem_sys state.
    // File: decode + accumulate (which fills gem_sys), then build JSON.
    if (uri.rfind("/api/gem/apv/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 13);
#ifdef WITH_ET
        if (mode_.load() == Mode::Online) {
            std::lock_guard<std::mutex> lk(ring_mtx_);
            for (auto &e : ring_) {
                if (e.seq == evnum) {
                    if (e.gem_apv_str.empty())
                        reply("{\"error\":\"gem apv pending\"}");
                    else
                        reply(e.gem_apv_str, "application/json", &e.gem_apv_gz);
                    return;
                }
            }
            reply("{\"error\":\"event not in ring buffer\"}"); return;
        }
#endif
        auto event_ptr = std::make_unique<fdec::EventData>();
        auto ssp_ptr   = std::make_unique<ssp::SspEventData>();
        std::string err = decodeRawEvent(evnum, *event_ptr, ssp_ptr.get());
        if (!err.empty()) { reply(json({{"error", err}}).dump()); return; }
        accumulate(evnum, *event_ptr, ssp_ptr.get());
        // accumulate() dedupes by event id, so on a re-request gem_sys
        // may still hold a different event's working buffers.  Force a
        // re-process (no histogram side effects) so apiGemApv reads the
        // requested event regardless of cache state.
        activeApp().prepareGemForView(*ssp_ptr);
        reply(activeApp().apiGemApv(*ssp_ptr, evnum).dump());
        return;
    }

    // --- ring buffer ---
    if (uri == "/api/ring") {
#ifdef WITH_ET
        std::lock_guard<std::mutex> lk(ring_mtx_);
        json arr = json::array();
        for (auto &e : ring_) arr.push_back(e.seq);
        reply(json({{"ring", arr},
                     {"latest", ring_.empty() ? 0 : ring_.back().seq}}).dump());
#else
        reply(json({{"ring", json::array()}, {"latest", 0}}).dump());
#endif
        return;
    }

    // --- progress ---
    if (uri == "/api/progress") { reply(progress_.toJson().dump()); return; }

    // MONITOR STATUS — header panel for online mode.  Each value may be <0
    // to mean "not available"; the frontend hides cells independently.
    //   livetime.ts       ← caget poll (AppState::livetime_cmd)
    //   livetime.measured ← DSC2 scaler in EVIO stream (activeApp())
    //   beam.energy       ← caget poll (AppState::beam_energy_status)
    //   beam.current      ← caget poll (AppState::beam_current_status)
    if (uri == "/api/monitor_status") {
        double ts = -1.0, be = -1.0, bc = -1.0;
        double meas = activeApp().measured_livetime.load();
#ifdef WITH_ET
        ts = livetime_.load();
        be = beam_energy_.load();
        bc = beam_current_.load();
#endif
        reply(json({
            {"livetime", {{"ts", ts}, {"measured", meas}}},
            {"beam",     {{"energy", be}, {"current", bc}}},
        }).dump());
        return;
    }

    // --- clear endpoints (always available, clears active mode's data) ---
    // --- filter endpoints ---
    if (uri == "/api/filter") {
        reply(activeApp().filterToJson().dump()); return;
    }
    if (uri == "/api/filter/load") {
        std::string body = con->get_request_body();
        auto fj = json::parse(body, nullptr, false);
        if (fj.is_discarded()) { reply("{\"error\":\"invalid JSON\"}"); return; }
        std::string err = applyFilter(fj);
        if (!err.empty()) { reply(json({{"error", err}}).dump()); return; }
        reply(json({{"status", "ok"}, {"filter", activeApp().filterToJson()},
                     {"filtered_count", (int)filtered_indices_.size()}}).dump());
        return;
    }
    if (uri == "/api/filter/unload") {
        clearFilter();
        reply(json({{"status", "ok"}, {"filter_active", false}}).dump());
        return;
    }
    if (uri == "/api/filter/indices") {
        reply(json(filtered_indices_).dump()); return;
    }

    // --- clear endpoints (always available, clears active mode's data) ---
    if (uri == "/api/hist/clear") {
        activeApp().clearHistograms();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"hist_cleared\"}");
        return;
    }
    if (uri == "/api/lms/clear") {
        activeApp().clearLms();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"lms_cleared\"}");
        return;
    }
    if (uri == "/api/epics/clear") {
        activeApp().clearEpics();
        reply("{\"cleared\":true}");
        wsBroadcast("{\"type\":\"epics_cleared\"}");
        return;
    }
    // --- shared read-only API routes ---
    auto result = activeApp().handleReadApi(uri);
    if (result.handled) { reply(result.body); return; }

    // --- file browser ---
    if (uri == "/api/files" || uri.rfind("/api/files?", 0) == 0) {
        std::string subdir;
        auto qpos = uri.find('?');
        if (qpos != std::string::npos) {
            std::string q = uri.substr(qpos + 1);
            if (q.rfind("dir=", 0) == 0) {
                subdir = q.substr(4);
                // URL-decode
                std::string dec;
                for (size_t i = 0; i < subdir.size(); ++i) {
                    if (subdir[i] == '%' && i + 2 < subdir.size()) {
                        dec += (char)std::stoi(subdir.substr(i + 1, 2), nullptr, 16);
                        i += 2;
                    } else if (subdir[i] == '+') dec += ' ';
                    else dec += subdir[i];
                }
                subdir = dec;
            }
        }
        reply(json({{"entries", listFiles(subdir)}}).dump()); return;
    }

    // --- load file (relative path from data_dir) ---
    if (uri.rfind("/api/load?", 0) == 0) {
        auto qpos = uri.find('?');
        std::string query = uri.substr(qpos + 1);
        std::string relpath;
        bool do_hist = false;
        for (size_t pos = 0; pos < query.size();) {
            size_t amp = query.find('&', pos);
            if (amp == std::string::npos) amp = query.size();
            std::string kv = query.substr(pos, amp - pos);
            auto eq = kv.find('=');
            if (eq != std::string::npos) {
                std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
                if (k == "file") relpath = v;
                if (k == "hist") do_hist = (v == "1");
            }
            pos = amp + 1;
        }
        // URL-decode %xx
        std::string decoded;
        for (size_t i = 0; i < relpath.size(); ++i) {
            if (relpath[i] == '%' && i + 2 < relpath.size()) {
                decoded += (char)std::stoi(relpath.substr(i + 1, 2), nullptr, 16);
                i += 2;
            } else if (relpath[i] == '+') decoded += ' ';
            else decoded += relpath[i];
        }
        relpath = decoded;

        std::string fullpath = resolveDataFile(relpath);
        if (fullpath.empty()) { reply("{\"error\":\"invalid path\"}"); return; }

        loadFile(fullpath, do_hist);
        reply(json({{"status", "loading"}, {"file", relpath},
                     {"hist_enabled", do_hist}}).dump());
        return;
    }

    // --- elog post ---
    if (uri == "/api/elog/post") {
        std::string body = con->get_request_body();
        json result = handleElogPost(body);
        if (!result.value("ok", true)) {
            con->set_status(result.contains("error") &&
                            result["error"] == "Empty body"
                            ? websocketpp::http::status_code::bad_request
                            : websocketpp::http::status_code::ok);
        } else {
            con->set_status(websocketpp::http::status_code::ok);
        }
        con->set_body(result.dump());
        con->append_header("Content-Type", "application/json");
        return;
    }

    // --- 404 ---
    con->set_status(websocketpp::http::status_code::not_found);
    con->set_body("404 Not Found");

    } catch (const std::exception &e) {
        con->set_status(websocketpp::http::status_code::bad_request);
        con->set_body(json({{"error", e.what()}, {"uri", uri}}).dump());
        con->append_header("Content-Type", "application/json");
        std::cerr << "[http] " << uri << " → " << e.what() << "\n";
    } catch (...) {
        con->set_status(websocketpp::http::status_code::internal_server_error);
        con->set_body("{\"error\":\"unknown exception\"}");
        con->append_header("Content-Type", "application/json");
        std::cerr << "[http] " << uri << " → unknown exception\n";
    }
}
