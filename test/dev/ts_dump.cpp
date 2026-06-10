// ts_dump — dump timestamp information from EVIO events for debugging
//
// Usage: ts_dump <input.evio> [-n max_events] [-D daq_config.json]
//
// Prints for each physics event:
//   event#, trigger_bits (hex), TI timestamp (raw ticks), time (sec from first),
//   unix_time (if available from sync/control events)

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <getopt.h>

using namespace evc;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

int main(int argc, char *argv[])
{
    std::string input, daq_config_file;
    int max_events = 50;

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;
    daq_config_file = db_dir + "/daq_config.json";

    int opt;
    while ((opt = getopt(argc, argv, "n:D:h")) != -1) {
        switch (opt) {
        case 'n': max_events = std::atoi(optarg); break;
        case 'D': daq_config_file = optarg; break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " <input.evio> [-n max_events] [-D daq_config.json]\n";
            return opt == 'h' ? 0 : 1;
        }
    }
    if (optind < argc) input = argv[optind];
    if (input.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.evio> [-n max_events] [-D daq_config.json]\n";
        return 1;
    }

    DaqConfig cfg;
    if (!load_daq_config(daq_config_file, cfg)) {
        std::cerr << "Failed to load DAQ config: " << daq_config_file << "\n";
        return 1;
    }

    EvChannel ch;
    ch.SetConfig(cfg);
    if (ch.OpenAuto(input) != status::success) {
        std::cerr << "Cannot open: " << input << "\n";
        return 1;
    }

    static constexpr double TI_TICK_SEC = 4e-9;

    auto event = std::make_unique<fdec::EventData>();
    uint64_t first_ts = 0;
    uint32_t last_sync_unix = 0;
    int total = 0;
    int buf = 0;
    int lms_count = 0;

    std::cout << std::left
              << std::setw(8)  << "Event#"
              << std::setw(12) << "Trigger"
              << std::setw(18) << "TI_Timestamp"
              << std::setw(14) << "dt(sec)"
              << std::setw(12) << "SyncUnix"
              << std::setw(8)  << "Type"
              << "\n"
              << std::string(72, '-') << "\n";

    while (ch.Read() == status::success) {
        ++buf;
        if (!ch.Scan()) continue;

        auto evtype = ch.GetEventType();

        // print control events (sync/prestart/go/end)
        if (evtype == EventType::Sync || evtype == EventType::Prestart ||
            evtype == EventType::Go || evtype == EventType::End) {
            uint32_t ct = ch.Sync().unix_time;
            if (ct != 0) last_sync_unix = ct;
            std::cout << std::setw(8)  << "--"
                      << std::setw(12) << "--"
                      << std::setw(18) << "--"
                      << std::setw(14) << "--"
                      << std::setw(12) << ct
                      << std::setw(8);
            if (evtype == EventType::Sync) std::cout << "SYNC";
            else if (evtype == EventType::Prestart) std::cout << "PRSTART";
            else if (evtype == EventType::Go) std::cout << "GO";
            else std::cout << "END";
            if (ct != 0)
                std::cout << "  unix=" << ct;
            std::cout << "\n";
            continue;
        }

        if (evtype == EventType::Epics) {
            std::cout << std::setw(8) << "--"
                      << std::setw(12) << "--"
                      << std::setw(18) << "--"
                      << std::setw(14) << "--"
                      << std::setw(12) << "--"
                      << std::setw(8) << "EPICS"
                      << "\n";
            continue;
        }

        if (evtype != EventType::Physics) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            if (!ch.DecodeEvent(ie, *event)) continue;

            if (first_ts == 0 && event->info.timestamp != 0)
                first_ts = event->info.timestamp;

            double dt = (first_ts != 0 && event->info.timestamp != 0)
                ? static_cast<double>(event->info.timestamp - first_ts) * TI_TICK_SEC
                : 0.0;

            bool is_lms = (event->info.trigger_bits & (1u << 24)) != 0; // LMS = bit 24 (database/trigger_bits.json)
            if (is_lms) lms_count++;

            std::cout << std::setw(8)  << event->info.event_number
                      << std::setw(12) << ("0x" + [](uint32_t v){
                             char buf[16]; snprintf(buf, sizeof(buf), "%08X", v);
                             return std::string(buf);
                         }(event->info.trigger_bits))
                      << std::setw(18) << event->info.timestamp
                      << std::setw(14) << std::fixed << std::setprecision(4) << dt
                      << std::setw(12) << (last_sync_unix ? std::to_string(last_sync_unix) : "--")
                      << std::setw(8) << (is_lms ? "LMS" : "")
                      << "\n";

            total++;
            if (max_events > 0 && total >= max_events) goto done;
        }
    }
done:
    ch.Close();

    std::cout << "\n--- Summary ---\n"
              << "Total physics events: " << total << "\n"
              << "LMS events: " << lms_count << "\n"
              << "First TI timestamp: " << first_ts << "\n"
              << "TI tick: " << TI_TICK_SEC << " sec (250 MHz)\n";
    if (first_ts != 0 && total > 0) {
        double total_time = static_cast<double>(event->info.timestamp - first_ts) * TI_TICK_SEC;
        std::cout << "Time span: " << std::fixed << std::setprecision(2)
                  << total_time << " sec\n";
    }
    if (last_sync_unix != 0)
        std::cout << "Last sync unix time: " << last_sync_unix << "\n";

    // also dump per-slot timestamps from the first event for comparison
    std::cout << "\n--- First event per-slot timestamps ---\n";
    ch.OpenAuto(input);
    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;
        event->clear();
        if (!ch.DecodeEvent(0, *event)) continue;
        std::cout << "EventInfo.timestamp = " << event->info.timestamp << "\n";
        for (int r = 0; r < event->nrocs; ++r) {
            auto &roc = event->rocs[r];
            if (!roc.present) continue;
            for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                if (!roc.slots[s].present) continue;
                std::cout << "  ROC 0x" << std::hex << roc.tag << std::dec
                          << " slot " << s
                          << " timestamp = " << roc.slots[s].timestamp << "\n";
            }
        }
        break;
    }
    ch.Close();

    return 0;
}
