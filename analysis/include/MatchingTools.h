#pragma once
//=============================================================================
// MatchingTools.h — detector matching tools for PRad2
//
// adapted from PRadAnalyzer/PRadDetMatch.cpp
//=============================================================================

#include "PhysicsTools.h"
#include <vector>
#include <cstdint>
#include <cmath>

namespace analysis {

// --- matching flag bit positions --------------------------------------------
enum MatchFlag : uint32_t {
    kGEM1Match = 0,
    kGEM2Match = 1,
    kGEM3Match = 2,
    kGEM4Match = 3,
};

struct ProjectHit
{
    float x_proj;
    float y_proj;
    float z_proj;

    ProjectHit(float x, float y, float z) : x_proj(x), y_proj(y), z_proj(z) {};
};

ProjectHit GetProjectionHits(float x, float y, float z, float projection_z);
void GetProjection(HCHit &hc, float projection_z);
void GetProjection(std::vector<HCHit> &hc, float projection_z);
void GetProjection(GEMHit &gem, float projection_z);
void GetProjection(std::vector<GEMHit> &gem, float projection_z);

class MatchHit
{
    public:
        analysis::HCHit hycal_hit;
        std::vector<analysis::GEMHit> gem1_hits;
        std::vector<analysis::GEMHit> gem2_hits;
        std::vector<analysis::GEMHit> gem3_hits;
        std::vector<analysis::GEMHit> gem4_hits;

        MatchHit(const analysis::HCHit &hycal_hit, const std::vector<analysis::GEMHit> &g1, const std::vector<analysis::GEMHit> &g2,
                 const std::vector<analysis::GEMHit> &g3, const std::vector<analysis::GEMHit> &g4)
            : hycal_hit(hycal_hit), gem1_hits(g1), gem2_hits(g2), gem3_hits(g3), gem4_hits(g4) {}

        // --- added for matching logic ----------------------------------------
        analysis::GEMHit gem[2];       // best-matched upstream and downstream GEM hits
        uint32_t   mflag = 0;     //E.g. 0101 matching flags (see MatchFlag enum)
        uint16_t    hycal_idx = 0; // index into original hycal vector
};

class MatchHit_perChamber
{
    public:
        analysis::HCHit hycal_hit;

        MatchHit_perChamber(const analysis::HCHit &hycal_hit)
            : hycal_hit(hycal_hit){}

        // --- added for matching logic ----------------------------------------
        float gem_hits[4][3] = {}; // store the best matches for each GEM chamber, format: [det_id][x/y/z]
        uint32_t   mflag = 0;     // matching flag
        uint16_t    hycal_idx = 0; // index into original hycal vector
};

class MatchingTools
{
public:
    MatchingTools() = default;

    std::vector<MatchHit> Match(const std::vector<analysis::HCHit> &hycalHits,
                            const std::vector<analysis::GEMHit> &gem1_hits,
                            const std::vector<analysis::GEMHit> &gem2_hits,
                            const std::vector<analysis::GEMHit> &gem3_hits,
                            const std::vector<analysis::GEMHit> &gem4_hits) const;

    std::vector<MatchHit_perChamber> MatchPerChamber(const std::vector<analysis::HCHit> &hycalHits,
                            const std::vector<analysis::GEMHit> &gem1_hits,
                            const std::vector<analysis::GEMHit> &gem2_hits,
                            const std::vector<analysis::GEMHit> &gem3_hits,
                            const std::vector<analysis::GEMHit> &gem4_hits) const;

    // configuration setters
    void SetMatchRange(float range)    { matchRange_ = range; }
    void SetSquareSelection(bool sq)   { squareSel_ = sq; }

private:
    float matchRange_ = 15.f;   // mm, spatial matching window
    bool  squareSel_  = true;   // true = square window, false = circular

    float ProjectionDistance(const analysis::HCHit &h, const analysis::GEMHit &g) const;
    float ProjectionDistance(const analysis::GEMHit &g1, const analysis::GEMHit &g2, float ref_z) const;
    bool  PreMatch(const analysis::HCHit &hycal, const analysis::GEMHit &gem) const;
    void  PostMatch(MatchHit &h) const;
};

} // namespace analysis
