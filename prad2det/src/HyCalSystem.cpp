//=============================================================================
// HyCalSystem.cpp — HyCal geometry, neighbor computation, and module lookup
//
// Ported from PRadHyCalDetector (PRadAnalyzer) with modernized data structures.
// The quantized-distance algorithm and neighbor criteria are preserved exactly.
//
// Chao Peng (original PRadAnalyzer), adapted for prad2decoder.
//=============================================================================

#include "HyCalSystem.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <iomanip>

using json = nlohmann::json;

namespace fdec
{

//=============================================================================
// Static helpers
//=============================================================================

ModuleType HyCalSystem::parse_type(const std::string &t)
{
    if (t == "PbGlass") return ModuleType::PbGlass;
    if (t == "PbWO4")   return ModuleType::PbWO4;
    if (t == "LMS")     return ModuleType::LMS;
    if (t == "Veto")    return ModuleType::Veto;
    return ModuleType::Unknown;
}

int HyCalSystem::name_to_id(const std::string &name)
{
    if (name.empty()) return -1;
    char prefix = name[0];
    try {
        if (prefix == 'W' || prefix == 'w')
            return std::stoi(name.substr(1)) + PWO_ID0;
        if (prefix == 'G' || prefix == 'g')
            return std::stoi(name.substr(1));
    } catch (...) {}
    return -1;
}

int HyCalSystem::id_to_index(int id) const
{
    auto it = id_map_.find(id);
    return (it != id_map_.end()) ? it->second : -1;
}

std::string HyCalSystem::id_to_name(int id)
{
    if (id < 0 || id > 2156) return "UNKNOWN";
    if (id >= PWO_ID0)
        return "W" + std::to_string(id - PWO_ID0);
    else
        return "G" + std::to_string(id);
}

// Line-segment intersection (from cana::intersection)
// Returns: -1 parallel, 0 within both segments, 1/2/3 outside segment(s)
int HyCalSystem::line_intersect(double x1, double y1, double x2, double y2,
                                double x3, double y3, double x4, double y4,
                                double &xc, double &yc)
{
    constexpr double inf = 0.001;
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (std::abs(denom) < inf) return -1;

    xc = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x3 * y4 - y3 * x4) * (x1 - x2)) / denom;
    yc = ((x1 * y2 - y1 * x2) * (y3 - y4) - (x3 * y4 - y3 * x4) * (y1 - y2)) / denom;

    bool out1 = ((x1 - xc) * (xc - x2) < -inf) || ((y1 - yc) * (yc - y2) < -inf);
    bool out2 = ((x3 - xc) * (xc - x4) < -inf) || ((y3 - yc) * (yc - y4) < -inf);

    return (out1 ? 1 : 0) + (out2 ? 2 : 0);
}

//=============================================================================
// Init — load modules + daq from one merged JSON, compute sectors, build neighbors
//=============================================================================

bool HyCalSystem::Init(const std::string &map_path)
{
    std::ifstream f(map_path);
    if (!f.is_open()) {
        std::cerr << "HyCalSystem: cannot open " << map_path << std::endl;
        return false;
    }

    json jdata;
    try { jdata = json::parse(f); }
    catch (const json::parse_error &e) {
        std::cerr << "HyCalSystem: JSON parse error in " << map_path
                  << ": " << e.what() << std::endl;
        return false;
    }

    if (!jdata.is_array()) {
        std::cerr << "HyCalSystem: " << map_path
                  << " is not a JSON array" << std::endl;
        return false;
    }

    modules_.clear();
    modules_.reserve(jdata.size());
    name_map_.clear();
    id_map_.clear();
    daq_map_.clear();

    for (auto &jm : jdata) {
        Module m;
        m.name   = jm.at("n").get<std::string>();
        m.type   = parse_type(jm.at("t").get<std::string>());
        const auto &geo = jm.at("geo");
        m.size_x = geo.at("sx").get<double>();
        m.size_y = geo.at("sy").get<double>();
        m.x      = geo.at("x").get<double>();
        m.y      = geo.at("y").get<double>();
        m.id     = name_to_id(m.name);
        m.index  = static_cast<int>(modules_.size());

        // DAQ block is optional — entries without it (e.g. boosters in
        // prad2hvmon, V1-V4 in PRad-1) leave m.daq at the default {-1,-1,-1}.
        if (jm.contains("daq")) {
            const auto &d = jm["daq"];
            m.daq = {d.at("crate").get<int>(),
                     d.at("slot").get<int>(),
                     d.at("channel").get<int>()};
        }

        modules_.push_back(std::move(m));
    }

    n_modules_ = static_cast<int>(modules_.size());

    // build lookup maps now that all modules are in their final positions
    for (int i = 0; i < n_modules_; ++i) {
        name_map_[modules_[i].name] = i;
        if (modules_[i].id >= 0)
            id_map_[modules_[i].id] = i;
        const auto &d = modules_[i].daq;
        if (d.crate >= 0)
            daq_map_[pack_daq(d.crate, d.slot, d.channel)] = i;
    }

    // --- compute sectors, layout, grids, and neighbors ----------------------
    compute_sectors();

    for (auto &m : modules_)
        assign_layout(m);

    build_sector_grids();
    build_neighbors();

    return true;
}

//=============================================================================
// Sector ID from position relative to center boundary
//=============================================================================

static int sector_from_center_bounds(double x, double y,
                                     double cx1, double cy1, double cx2, double cy2)
{
    if (x > cx2) return (y < cy1) ? static_cast<int>(Sector::Bottom) : static_cast<int>(Sector::Right);
    if (x < cx1) return (y > cy2) ? static_cast<int>(Sector::Top) : static_cast<int>(Sector::Left);
    if (y > cy2) return static_cast<int>(Sector::Top);
    if (y < cy1) return static_cast<int>(Sector::Bottom);
    return static_cast<int>(Sector::Center);
}

int HyCalSystem::get_sector_id(double x, double y) const
{
    double x1, y1, x2, y2;
    sectors_[static_cast<int>(Sector::Center)].get_boundary(x1, y1, x2, y2);
    return sector_from_center_bounds(x, y, x1, y1, x2, y2);
}

//=============================================================================
// Sector computation — find bounding boxes for each sector
//=============================================================================

void HyCalSystem::compute_sectors()
{
    constexpr int NS = static_cast<int>(Sector::Max);
    std::array<bool, NS> inited{};
    std::array<double, NS> xmin{}, ymin{}, xmax{}, ymax{};

    // First pass: find PbWO4 (center) bounding box
    double cx1 = 1e9, cy1 = 1e9, cx2 = -1e9, cy2 = -1e9;
    for (auto &m : modules_) {
        if (m.type != ModuleType::PbWO4) continue;
        double hx = m.size_x / 2., hy = m.size_y / 2.;
        cx1 = std::min(cx1, m.x - hx);
        cy1 = std::min(cy1, m.y - hy);
        cx2 = std::max(cx2, m.x + hx);
        cy2 = std::max(cy2, m.y + hy);
    }

    // Second pass: assign sector to each module and accumulate boundaries
    for (auto &m : modules_) {
        if (!m.is_hycal()) {
            m.sector = -1;
            continue;
        }

        if (m.type == ModuleType::PbWO4) {
            m.sector = static_cast<int>(Sector::Center);
        } else {
            m.sector = sector_from_center_bounds(m.x, m.y, cx1, cy1, cx2, cy2);
        }

        int s = m.sector;
        if (s < 0 || s >= NS) continue;

        double hx = m.size_x / 2., hy = m.size_y / 2.;
        if (inited[s]) {
            xmin[s] = std::min(xmin[s], m.x - hx);
            ymin[s] = std::min(ymin[s], m.y - hy);
            xmax[s] = std::max(xmax[s], m.x + hx);
            ymax[s] = std::max(ymax[s], m.y + hy);
        } else {
            inited[s] = true;
            auto &si   = sectors_[s];
            si.id      = s;
            si.mtype   = m.type;
            si.msize_x = m.size_x;
            si.msize_y = m.size_y;
            xmin[s] = m.x - hx;
            ymin[s] = m.y - hy;
            xmax[s] = m.x + hx;
            ymax[s] = m.y + hy;
        }
    }

    for (int i = 0; i < NS; ++i) {
        if (inited[i])
            sectors_[i].set_boundary(xmin[i], ymin[i], xmax[i], ymax[i]);
    }
}

//=============================================================================
// Layout assignment — set flags, row, column (PrimEx-specific)
//=============================================================================

void HyCalSystem::assign_layout(Module &m) const
{
    if (!m.is_hycal()) return;

    m.flag = 0;

    if (m.is_pwo4()) {
        set_bit(m.flag, kPbWO4);
        int pid = m.id - PWO_ID0 - 1;
        m.row    = pid / 34;
        m.column = pid % 34;

        // inner boundary (beam hole): rows/cols 15-18 (0-based)
        if (m.row >= 15 && m.row <= 18 && m.column >= 15 && m.column <= 18)
            set_bit(m.flag, kInnerBound);
        // transition: outer ring of crystal array
        if (m.row == 0 || m.row == 33 || m.column == 0 || m.column == 33)
            set_bit(m.flag, kTransition);

    } else {
        set_bit(m.flag, kPbGlass);
        int pid = m.id - 1;
        int g_row = pid / 30;       // 0-based
        int g_col = pid % 30;

        // outer boundary
        if (g_row == 0 || g_row == 29 || g_col == 0 || g_col == 29)
            set_bit(m.flag, kOuterBound);

        // determine which glass sector and set transition flags
        // Top sector: g_col < 24, g_row < 6
        if (g_col < 24 && g_row < 6) {
            m.row = g_row;
            m.column = g_col;
            if (g_row == 5 && g_col >= 5)
                set_bit(m.flag, kTransition);
        }
        // Right sector: g_col >= 24, g_row < 24
        else if (g_col >= 24 && g_row < 24) {
            m.row = g_row;
            m.column = g_col - 24;
            if (g_col == 24 && g_row >= 5)
                set_bit(m.flag, kTransition);
        }
        // Bottom sector: g_col >= 6, g_row >= 24
        else if (g_col >= 6 && g_row >= 24) {
            m.row = g_row - 24;
            m.column = g_col - 6;
            if (g_row == 24 && g_col - 6 < 19)
                set_bit(m.flag, kTransition);
        }
        // Left sector: g_col < 6, g_row >= 6
        else if (g_col < 6 && g_row >= 6) {
            m.row = g_row - 6;
            m.column = g_col;
            if (g_col == 5 && g_row - 6 < 19)
                set_bit(m.flag, kTransition);
        }
    }
}

//=============================================================================
// Quantized distance — the core algorithm from PRadHyCalDetector
//=============================================================================

void HyCalSystem::qdist(double x1, double y1, int s1,
                        double x2, double y2, int s2,
                        double &dx, double &dy) const
{
    const auto &sec1 = sectors_[s1];
    const auto &sec2 = sectors_[s2];

    if (s1 != s2) {
        const auto &center = sectors_[static_cast<int>(Sector::Center)];
        const auto &boundary = center.boundpts;

        // different sectors with different module types
        if (sec1.mtype != sec2.mtype) {
            for (size_t i = 0; i < boundary.size(); ++i) {
                size_t ip = (i == 0) ? boundary.size() - 1 : i - 1;
                auto &p1 = boundary[ip];
                auto &p2 = boundary[i];

                double xc = 0., yc = 0.;
                int inter = line_intersect(p1.x, p1.y, p2.x, p2.y,
                                           x1, y1, x2, y2, xc, yc);
                if (inter == 0) {
                    dx = (x2 - xc) / sec2.msize_x + (xc - x1) / sec1.msize_x;
                    dy = (y2 - yc) / sec2.msize_y + (yc - y1) / sec1.msize_y;
                    return;
                }
            }
        }
        // different sectors, same module type — may cross through center
        else {
            double xc[2], yc[2];
            size_t ic = 0;
            for (size_t i = 0; i < boundary.size(); ++i) {
                size_t ip = (i == 0) ? boundary.size() - 1 : i - 1;
                auto &p1 = boundary[ip];
                auto &p2 = boundary[i];

                int inter = line_intersect(p1.x, p1.y, p2.x, p2.y,
                                           x1, y1, x2, y2,
                                           xc[ic], yc[ic]);
                if (inter == 0 && ic++ > 0) break;
            }

            if (ic > 1) {
                double dxt = x2 - x1, dyt = y2 - y1;
                double dxc = xc[0] - xc[1], dyc = yc[0] - yc[1];
                double sign = (dxt * dxc > 0.) ? 1. : -1.;
                dxc *= sign;
                dyc *= sign;

                dx = (dxt - dxc) / sec1.msize_x + dxc / center.msize_x;
                dy = (dyt - dyc) / sec1.msize_y + dyc / center.msize_y;
                return;
            }
        }
    }

    // same sector or fallback
    dx = (x2 - x1) / sec1.msize_x;
    dy = (y2 - y1) / sec1.msize_y;
}

//=============================================================================
// Build sector grids — map (row, col) → module index per sector
//=============================================================================

void HyCalSystem::build_sector_grids()
{
    // Grid dimensions per sector (must match assign_layout row/col ranges)
    static constexpr int dims[][2] = {
        {34, 34},   // Center (PbWO4)
        { 6, 24},   // Top
        {24,  6},   // Right
        { 6, 24},   // Bottom
        {24,  6},   // Left
    };

    for (int s = 0; s < static_cast<int>(Sector::Max); ++s)
        sector_grids_[s].init(dims[s][0], dims[s][1]);

    for (int i = 0; i < n_modules_; ++i) {
        auto &m = modules_[i];
        if (m.sector < 0) continue;
        auto &grid = sector_grids_[m.sector];
        if (grid.valid(m.row, m.column))
            grid.at(m.row, m.column) = i;
    }
}

//=============================================================================
// Build cross-sector neighbors — only edge modules need qdist
//=============================================================================

void HyCalSystem::build_neighbors()
{
    for (auto &m : modules_) m.neighbor_count = 0;

    // Identify edge modules: any HyCal module with a missing/empty grid neighbor
    std::vector<int> edge_mods;
    edge_mods.reserve(300);

    for (int i = 0; i < n_modules_; ++i) {
        auto &m = modules_[i];
        if (!m.is_hycal() || m.sector < 0) continue;
        const auto &grid = sector_grids_[m.sector];
        bool edge = false;
        for (int dr = -1; dr <= 1 && !edge; ++dr)
            for (int dc = -1; dc <= 1 && !edge; ++dc) {
                if (dr == 0 && dc == 0) continue;
                int nr = m.row + dr, nc = m.column + dc;
                if (!grid.valid(nr, nc) || grid.at(nr, nc) < 0)
                    edge = true;
            }
        if (edge) edge_mods.push_back(i);
    }

    // Compute cross-sector neighbors among edge module pairs
    for (size_t a = 0; a < edge_mods.size(); ++a) {
        auto &m1 = modules_[edge_mods[a]];
        for (size_t b = a + 1; b < edge_mods.size(); ++b) {
            auto &m2 = modules_[edge_mods[b]];
            if (m1.sector == m2.sector) continue;

            double dx, dy;
            qdist(m1.x, m1.y, m1.sector, m2.x, m2.y, m2.sector, dx, dy);

            if (std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
                float fdx = static_cast<float>(dx);
                float fdy = static_cast<float>(dy);
                float dist = std::sqrt(fdx * fdx + fdy * fdy);

                int i = edge_mods[a], j = edge_mods[b];
                if (m1.neighbor_count < MAX_NEIGHBORS)
                    m1.neighbors[m1.neighbor_count++] = {j, fdx, fdy, dist};
                if (m2.neighbor_count < MAX_NEIGHBORS)
                    m2.neighbors[m2.neighbor_count++] = {i, -fdx, -fdy, dist};
            }
        }
    }
}

//=============================================================================
// Module lookup
//=============================================================================

const Module *HyCalSystem::module_by_name(const std::string &name) const
{
    auto it = name_map_.find(name);
    return (it != name_map_.end()) ? &modules_[it->second] : nullptr;
}

const Module *HyCalSystem::module_by_id(int primex_id) const
{
    auto it = id_map_.find(primex_id);
    return (it != id_map_.end()) ? &modules_[it->second] : nullptr;
}

const Module *HyCalSystem::module_by_daq(int crate, int slot, int ch) const
{
    auto it = daq_map_.find(pack_daq(crate, slot, ch));
    return (it != daq_map_.end()) ? &modules_[it->second] : nullptr;
}

double HyCalSystem::GetCalibConstant(int primex_id) const
{
    const Module *m = module_by_id(primex_id);
    return m ? m->cal_factor : 0.;
}

double HyCalSystem::GetCalibBaseEnergy(int primex_id) const
{
    const Module *m = module_by_id(primex_id);
    return m ? m->cal_base_energy : 0.;
}

double HyCalSystem::GetCalibNonLinearity1(int primex_id) const
{
    const Module *m = module_by_id(primex_id);
    return m ? m->cal_non_linear_1 : 0.;
}

double HyCalSystem::GetCalibNonLinearity2(int primex_id) const
{
    const Module *m = module_by_id(primex_id);
    return m ? m->cal_non_linear_2 : 0.;
}

void HyCalSystem::SetCalibConstant(int primex_id, double factor)
{
    auto it = id_map_.find(primex_id);
    if (it != id_map_.end())
        modules_[it->second].cal_factor = factor;
}

void HyCalSystem::SetCalibBaseEnergy(int primex_id, double energy)
{
    auto it = id_map_.find(primex_id);
    if (it != id_map_.end())
        modules_[it->second].cal_base_energy = energy;
}

void HyCalSystem::SetCalibNonLinearity(int primex_id, double nl)
{
    auto it = id_map_.find(primex_id);
    if (it != id_map_.end()) {
        modules_[it->second].cal_non_linear_1 = nl;
        modules_[it->second].cal_non_linear_2 = 0.;
    }
}

void HyCalSystem::SetCalibNonLinearity(int primex_id, double nl1, double nl2)
{
    auto it = id_map_.find(primex_id);
    if (it != id_map_.end()) {
        modules_[it->second].cal_non_linear_1 = nl1;
        modules_[it->second].cal_non_linear_2 = nl2;
    }
}

void HyCalSystem::PrintCalibConstants(const std::string &output_file) const
{
    std::ofstream f(output_file);
    if (!f.is_open()) {
        std::cerr << "HyCalSystem::PrintCalibConstants: cannot open " << output_file << "\n";
        return;
    }

    json j = json::array();
    for (const auto &m : modules_) {
        if (m.id < 0) continue;
        j.push_back({
            {"name",        m.name},
            {"factor",      m.cal_factor},
            {"base_energy", m.cal_base_energy},
            {"nl1",         m.cal_non_linear_1},
            {"nl2",         m.cal_non_linear_2}
        });
    }
    f << j.dump(2) << "\n";
}

int HyCalSystem::LoadCalibration(const std::string &calib_path)
{
    std::ifstream f(calib_path);
    if (!f.is_open()) {
        std::cerr << "HyCalSystem::LoadCalibration: cannot open " << calib_path << "\n";
        return -1;
    }

    json j;
    try { j = json::parse(f, nullptr, true, true); }
    catch (const json::parse_error &e) {
        std::cerr << "HyCalSystem::LoadCalibration: parse error: " << e.what() << "\n";
        return -1;
    }

    int matched = 0;
    for (auto &entry : j) {
        std::string name = entry.value("name", "");
        if (name.empty()) continue;

        const Module *m = module_by_name(name);
        if (!m) continue;

        Module &mod = modules_[m->index];
        mod.cal_factor      = entry.value("factor", 0.0);
        mod.cal_base_energy = entry.value("base_energy", 0.0);
        mod.cal_non_linear_1 = entry.value("nl1", 0.0);
        mod.cal_non_linear_2 = entry.value("nl2", 0.0);
        ++matched;
    }

    return matched;
}

} // namespace fdec
