#pragma once
//=============================================================================
// HyCalSystem.h — HyCal detector geometry, module lookup, and neighbor system
//
// Initialized once from JSON config files. Immutable after Init().
// Provides O(1) module lookup by name, ID, or DAQ address, plus pre-computed
// neighbor lists with quantized distances for clustering.
//
// No per-event state lives here — callers provide flat energy arrays.
//=============================================================================

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <cmath>

namespace fdec
{

// --- capacity limits --------------------------------------------------------
static constexpr int HYCAL_MAX_MODULES  = 1800;  // 1728 real + margin
static constexpr int MAX_NEIGHBORS      = 12;    // max neighbors per module

// --- module types -----------------------------------------------------------
enum class ModuleType : int8_t { PbGlass = 0, PbWO4 = 1, LMS = 2, Veto = 3, Unknown = -1 };

// crystal module IDs start from 1000 (PrimEx convention)
static constexpr int PWO_ID0 = 1000;

// --- layout flags (bit positions) -------------------------------------------
enum LayoutFlag : uint32_t {
    kPbGlass     = 0,
    kPbWO4       = 1,
    kTransition  = 2,   // boundary between PbWO4 and PbGlass
    kInnerBound  = 3,   // inner boundary of HyCal (beam hole)
    kOuterBound  = 4,   // outer boundary of HyCal
    kDeadModule  = 5,
    kDeadNeighbor = 6,
    // cluster flags (used on ModuleCluster::flag, not module layout)
    kSplit       = 7,   // cluster was split from a multi-maximum group
    kLeakCorr    = 8,   // leakage correction applied
};

inline void     set_bit(uint32_t &f, uint32_t b) { f |= (1u << b); }
inline void     clear_bit(uint32_t &f, uint32_t b) { f &= ~(1u << b); }
inline bool     test_bit(uint32_t f, uint32_t b) { return (f >> b) & 1u; }

// --- sector types -----------------------------------------------------------
enum class Sector : int { Center = 0, Top = 1, Right = 2, Bottom = 3, Left = 4, Max = 5 };

struct SectorInfo {
    int         id      = -1;
    ModuleType  mtype   = ModuleType::Unknown;
    double      msize_x = 0.;
    double      msize_y = 0.;
    // boundary rectangle: [x_min, y_min] to [x_max, y_max]
    // stored as 4 corner points (BL, TL, TR, BR) for intersection tests
    struct Point { double x, y; };
    std::array<Point, 4> boundpts{};

    void set_boundary(double x1, double y1, double x2, double y2)
    {
        // counter-clockwise: BL, TL, TR, BR  (matches old code: [3]=TL,[0]=BL,[1]=TR,[2]=BR)
        // old code order: (x1,y2), (x1,y1), (x2,y1), (x2,y2)
        boundpts[0] = {x1, y2};  // top-left
        boundpts[1] = {x1, y1};  // bottom-left
        boundpts[2] = {x2, y1};  // bottom-right
        boundpts[3] = {x2, y2};  // top-right
    }

    void get_boundary(double &x1, double &y1, double &x2, double &y2) const
    {
        x1 = boundpts[1].x;  y1 = boundpts[1].y;
        x2 = boundpts[3].x;  y2 = boundpts[3].y;
    }
};

// --- neighbor info ----------------------------------------------------------
struct NeighborInfo {
    int   index;        // index into Module array
    float dx, dy;       // quantized distances
    float dist;         // sqrt(dx*dx + dy*dy)
};

// --- sector grid for O(1) neighbor lookup -----------------------------------
static constexpr int SECTOR_GRID_MAX_CELLS = 34 * 34;  // Center is largest

struct SectorGrid {
    int nrows = 0, ncols = 0;
    int cells[SECTOR_GRID_MAX_CELLS];

    void init(int r, int c) {
        nrows = r; ncols = c;
        for (int i = 0; i < r * c; ++i) cells[i] = -1;
    }
    bool valid(int r, int c) const { return r >= 0 && r < nrows && c >= 0 && c < ncols; }
    int  at(int r, int c) const    { return cells[r * ncols + c]; }
    int& at(int r, int c)          { return cells[r * ncols + c]; }
};

// --- DAQ address ------------------------------------------------------------
struct DaqAddr {
    int crate   = -1;
    int slot    = -1;
    int channel = -1;
};

// --- module -----------------------------------------------------------------
struct Module {
    // identity
    std::string name;
    int         id          = -1;       // PrimEx ID (G: 1-576, W: 1001-2152)
    int         index       = -1;       // contiguous array index

    // geometry
    ModuleType  type        = ModuleType::Unknown;
    double      x           = 0.;
    double      y           = 0.;
    double      size_x      = 0.;
    double      size_y      = 0.;

    // layout
    uint32_t    flag        = 0;
    int         sector      = -1;       // Sector enum as int
    int         row         = 0;
    int         column      = 0;

    // DAQ mapping
    DaqAddr     daq;

    // calibration (per-module)
    double      cal_factor      = 0.;   // MeV per ADC count (or per integral unit)
    double      cal_base_energy = 0.;   // calibration beam energy (MeV)
    double      cal_non_linear_1  = 0.;   // 1st order non-linear correction coefficient
    double      cal_non_linear_2  = 0.;   // 2nd order non-linear correction coefficient

    // convert ADC value (pedestal-subtracted) to energy in MeV
    double energize(double adc) const {
        if (adc < 0.) return 0.;
        double E = cal_factor * adc;
        return E;
    }

    // pre-computed cross-sector neighbors (same-sector uses grid lookup)
    int          neighbor_count = 0;
    NeighborInfo neighbors[MAX_NEIGHBORS];

    // helpers
    bool is_pwo4()   const { return type == ModuleType::PbWO4; }
    bool is_glass()  const { return type == ModuleType::PbGlass; }
    bool is_hycal()  const { return is_pwo4() || is_glass(); }

    bool is_neighbor(const Module &other, bool include_corners = true) const
    {
        // O(1) same-sector: row/col adjacency on uniform grid
        if (sector >= 0 && sector == other.sector) {
            int dr = std::abs(row - other.row);
            int dc = std::abs(column - other.column);
            if (dr > 1 || dc > 1 || (dr == 0 && dc == 0)) return false;
            return include_corners || (dr + dc <= 1);
        }
        // cross-sector: scan pre-computed list
        for (int i = 0; i < neighbor_count; ++i) {
            if (neighbors[i].index == other.index)
                return include_corners || neighbors[i].dist < 1.2f;
        }
        return false;
    }
};

// --- HyCalSystem ------------------------------------------------------------
class HyCalSystem
{
public:
    HyCalSystem() = default;

    // Initialize from the merged HyCal map JSON.
    // map_path: hycal_map.json — array of per-module records of the form
    //   {"n": ..., "t": ..., "geo": {sx,sy,x,y,sec,row,col},
    //                          "daq": {crate,slot,channel}}.
    // The "daq" block is optional (e.g. PRad-1 has no daq for V1-V4).
    // Returns false on error.
    bool Init(const std::string &map_path);

    // Load per-module calibration from JSON file.
    // Format: [{"name":"W735","factor":0.37,"base_energy":2138.67,"nl1":0.006,"nl2":0.0}, ...]
    // Returns number of modules matched, or -1 on error.
    int LoadCalibration(const std::string &calib_path);

    // --- module access (all O(1)) -------------------------------------------
    int             module_count()                     const { return n_modules_; }
    const Module   &module(int index)                  const { return modules_[index]; }
    Module         &module(int index)                        { return modules_[index]; }
    const Module   *module_by_name(const std::string &name) const;
    const Module   *module_by_id(int primex_id)             const;
    const Module   *module_by_daq(int crate, int slot, int ch) const;

    // --- calibration accessors ----------------------------------------------
    double GetCalibConstant(int primex_id) const;
    double GetCalibBaseEnergy(int primex_id) const;
    double GetCalibNonLinearity1(int primex_id) const;
    double GetCalibNonLinearity2(int primex_id) const;
    void   SetCalibConstant(int primex_id, double factor);
    void   SetCalibBaseEnergy(int primex_id, double energy);
    void   SetCalibNonLinearity(int primex_id, double nl);
    void   SetCalibNonLinearity(int primex_id, double nl1, double nl2);
    void   PrintCalibConstants(const std::string &output_file) const;

    // --- sector info --------------------------------------------------------
    const SectorInfo &sector_info(int s)               const { return sectors_[s]; }
    int  get_sector_id(double x, double y)             const;

    // --- quantized distance -------------------------------------------------
    void qdist(double x1, double y1, int s1,
               double x2, double y2, int s2,
               double &dx, double &dy)                 const;

    void qdist(const Module &m1, const Module &m2,
               double &dx, double &dy)                 const
    {
        qdist(m1.x, m1.y, m1.sector, m2.x, m2.y, m2.sector, dx, dy);
    }

    // --- neighbor enumeration (grid + cross-sector) --------------------------
    // Calls fn(neighbor_module_index) for each neighbor of the given module.
    // Same-sector: walks ±1 in the sector grid (O(8)).
    // Cross-sector: iterates pre-computed list (typically 0-3 entries).
    template<typename Fn>
    void for_each_neighbor(int module_index, bool include_corners, Fn &&fn) const
    {
        const auto &m = modules_[module_index];
        if (m.sector < 0) return;
        const auto &grid = sector_grids_[m.sector];
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                if (dr == 0 && dc == 0) continue;
                if (!include_corners && dr != 0 && dc != 0) continue;
                int nr = m.row + dr, nc = m.column + dc;
                if (grid.valid(nr, nc)) {
                    int ni = grid.at(nr, nc);
                    if (ni >= 0) fn(ni);
                }
            }
        }
        for (int i = 0; i < m.neighbor_count; ++i) {
            if (!include_corners && m.neighbors[i].dist >= 1.2f) continue;
            fn(m.neighbors[i].index);
        }
    }

    // --- position resolution at HyCal face ----------------------------------
    // sigma(E) = sqrt( (A / sqrt(E_GeV))^2 + (B / E_GeV)^2 + C^2 )   [mm]
    //   A — stochastic term [mm·sqrt(GeV)]
    //   B — noise term       [mm·GeV]
    //   C — constant term    [mm]
    // Defaults are zero except C=5 mm so that an uninitialized HyCalSystem
    // returns a sensible (conservative) value.  Loaded from
    // reconstruction_config.json:matching:hycal_pos_res = [A, B, C].
    void  SetPositionResolutionParams(float A, float B, float C)
    {
        pos_res_A_ = A;
        pos_res_B_ = B;
        pos_res_C_ = C;
    }
    float GetPositionResolutionA() const { return pos_res_A_; }
    float GetPositionResolutionB() const { return pos_res_B_; }
    float GetPositionResolutionC() const { return pos_res_C_; }
    float PositionResolution(float energy_mev) const
    {
        const float E_GeV = (energy_mev > 0.f) ? energy_mev / 1000.f : 1e-6f;
        const float a = pos_res_A_ / std::sqrt(E_GeV);
        const float b = pos_res_B_ / E_GeV;
        return std::sqrt(a * a + b * b + pos_res_C_ * pos_res_C_);
    }

    // --- static helpers -----------------------------------------------------
    static int          name_to_id(const std::string &name);
    int                 id_to_index(int id) const;
    static std::string  id_to_name(int id);
    static ModuleType   parse_type(const std::string &t);

private:
    void  compute_sectors();
    void  assign_layout(Module &m) const;
    void  build_sector_grids();
    void  build_neighbors();

    // line-segment intersection (ported from cana::intersection)
    static int line_intersect(double x1, double y1, double x2, double y2,       
                              double x3, double y3, double x4, double y4,
                              double &xc, double &yc);

    int                  n_modules_ = 0;
    std::vector<Module>  modules_;

    // position resolution coefficients (see PositionResolution above)
    float                pos_res_A_ = 0.f;   // mm·sqrt(GeV)
    float                pos_res_B_ = 0.f;   // mm·GeV
    float                pos_res_C_ = 5.f;   // mm

    // lookup maps
    std::unordered_map<std::string, int> name_map_;     // name → index
    std::unordered_map<int, int>         id_map_;       // primex_id → index
    std::unordered_map<uint64_t, int>    daq_map_;      // packed daq addr → index

    // sector info
    std::array<SectorInfo, static_cast<int>(Sector::Max)> sectors_;

    // sector grids for O(1) same-sector neighbor lookup
    std::array<SectorGrid, static_cast<int>(Sector::Max)> sector_grids_;

    // pack DAQ address into a single key
    static uint64_t pack_daq(int crate, int slot, int ch)
    {
        return (static_cast<uint64_t>(crate) << 32) |
               (static_cast<uint64_t>(slot)  << 16) |
               static_cast<uint64_t>(ch);
    }
};

} // namespace fdec
