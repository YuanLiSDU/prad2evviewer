//=============================================================================
// HyCalCluster.cpp — Island clustering for HyCal
//
// Ported from PRadIslandCluster / PRadHyCalReconstructor (PRadAnalyzer).
// Uses pre-computed neighbor lists from HyCalSystem for fast adjacency checks.
//
// Chao Peng (original PRadAnalyzer), adapted for prad2decoder.
//=============================================================================

#include "HyCalCluster.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace fdec
{

static constexpr int ISLAND_GROUP_RESERVE = 50;
static constexpr int POS_RECON_HITS       = 15;

//=============================================================================
// shower_depth (free function — declared in HyCalCluster.h)
//
// Maximum-shower-development depth into the calorimeter face for an EM
// shower of energy E, t = X0 · (ln(E/Ec) − Cf), with Cf = 0.5 for photons.
// Constants verbatim from the legacy analysis::PhysicsTools::GetShowerDepth
// (PRadAnalyzer lineage) and not expected to change — they're physical
// properties of the W and G modules.
//=============================================================================
float shower_depth(int center_id, float energy_mev)
{
    if (energy_mev <= 0.f) return 0.f;
    if (center_id >= PWO_ID0)            // PbWO4 (W-modules)
        return 8.6f  * (std::log(energy_mev / 1.1f)  - 0.5f);
    return            26.7f * (std::log(energy_mev / 2.84f) - 0.5f);  // PbGlass
}

//=============================================================================
// Construction / setup
//=============================================================================

HyCalCluster::HyCalCluster(const HyCalSystem &sys)
    : sys_(sys)
    , profile_(new SimpleProfile())
    , owns_profile_(true)
{
}

HyCalCluster::~HyCalCluster()
{
    if (owns_profile_) delete profile_;
}

void HyCalCluster::SetProfile(IClusterProfile *prof)
{
    if (owns_profile_) delete profile_;
    profile_ = prof;
    owns_profile_ = false;
}

//=============================================================================
// Per-event interface
//=============================================================================

void HyCalCluster::Clear()
{
    hits_.clear();
    groups_.clear();
    clusters_.clear();
}

void HyCalCluster::AddHit(int module_index, float energy, float time)
{
    if (module_index < 0 || module_index >= sys_.module_count()) return;
    if (energy > config_.min_module_energy)
        hits_.push_back({module_index, energy, time});
}

void HyCalCluster::FormClusters()
{
    clusters_.clear();
    groups_.clear();

    // step 1: group adjacent hits using DFS
    group_hits();

    // step 2: find maxima and split each group
    for (auto &group : groups_)
        split_cluster(group);
}

void HyCalCluster::ReconstructHits(std::vector<ClusterHit> &out) const
{
    out.clear();
    out.reserve(clusters_.size());

    for (auto &cl : clusters_) {
        if (cl.energy < config_.min_cluster_energy) continue;
        if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) continue;
        out.push_back(reconstruct_pos(cl));
    }
    std::sort(out.begin(), out.end(), [](const ClusterHit &a, const ClusterHit &b) {
        return a.energy > b.energy;
    });
}

void HyCalCluster::ReconstructMatched(std::vector<RecoResult> &out) const
{
    out.clear();
    out.reserve(clusters_.size());

    for (auto &cl : clusters_) {
        if (cl.energy < config_.min_cluster_energy) continue;
        if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) continue;
        out.push_back({&cl, reconstruct_pos(cl)});
    }
}

//=============================================================================
// Seed-driven BFS grouping (multi-pulse aware)
//
// Pulses are sorted by energy descending; the largest unconsumed pulse that
// passes min_center_energy becomes a cluster seed and grows an island via
// BFS through module neighbours.  For each neighbour module, the LARGEST
// unconsumed pulse whose time lies within ±seed_time_window of the seed is
// added to the group and marked consumed.  Pulses that fall outside the
// seed's window stay alive in the pool and can seed a later cluster at a
// different timing within the same event.
//
// When seed_time_window <= 0 the dt gate is bypassed and the algorithm
// degenerates to the legacy "all connected pulses" behaviour — single-pulse
// callers see no change.
//=============================================================================

void HyCalCluster::group_hits()
{
    // Build per-module pulse lists.  Multiple pulses on the same module
    // are common with FADC waveform data; mod_to_hits_[m] holds every
    // hit index whose ModuleHit::index == m.
    mod_to_hits_.assign(sys_.module_count(), {});
    for (int i = 0; i < static_cast<int>(hits_.size()); ++i)
        mod_to_hits_[hits_[i].index].push_back(i);

    hit_group_id_.assign(hits_.size(), -1);
    consumed_.assign(hits_.size(), false);

    // Energy-descending seed order — the global largest pulse seeds first,
    // ensuring that if multiple pulses on a module are time-coincident with
    // the seed, the most energetic shower claims the cluster.
    std::vector<int> order(hits_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return hits_[a].energy > hits_[b].energy; });

    for (int seed : order) {
        if (consumed_[seed]) continue;
        // sorted desc — once we drop below the seed threshold no later
        // pulse can satisfy it.
        if (hits_[seed].energy < config_.min_center_energy) break;

        int gid = static_cast<int>(groups_.size());
        groups_.emplace_back();
        groups_.back().reserve(ISLAND_GROUP_RESERVE);

        consumed_[seed] = true;
        hit_group_id_[seed] = gid;
        grow_island(seed, gid, groups_.back());
    }
}

void HyCalCluster::grow_island(int seed_idx, int group_id, std::vector<int> &group)
{
    const float seed_time = hits_[seed_idx].time;
    const bool  use_time  = config_.seed_time_window > 0.f;
    const float dt_max    = config_.seed_time_window;

    // BFS using the group vector itself as the queue: indices [qi..size())
    // are still to be expanded.  Seed pushed by caller; we expand the
    // frontier in insertion order.
    group.push_back(seed_idx);
    for (size_t qi = 0; qi < group.size(); ++qi) {
        int hi = group[qi];
        sys_.for_each_neighbor(hits_[hi].index, config_.corner_conn, [&](int ni) {
            // Per neighbour MODULE, pick at most one pulse to add to this
            // group: the LARGEST-energy unconsumed pulse whose time lies
            // within ±seed_time_window of the seed (or any unconsumed
            // pulse when gating is off).  Largest-amplitude is more
            // reliable than closest-in-time — small pulses have noisy
            // peak times.  Other pulses on the same module remain in
            // the pool for a different seed at a different timing.
            int best_k = -1;
            for (int k : mod_to_hits_[ni]) {
                if (consumed_[k]) continue;
                if (use_time && std::fabs(hits_[k].time - seed_time) > dt_max)
                    continue;
                if (best_k < 0 || hits_[k].energy > hits_[best_k].energy)
                    best_k = k;
            }
            if (best_k >= 0) {
                consumed_[best_k] = true;
                hit_group_id_[best_k] = group_id;
                group.push_back(best_k);
            }
        });
    }
}

//=============================================================================
// Split cluster — find local maxima, distribute hits
//=============================================================================

void HyCalCluster::split_cluster(const std::vector<int> &group)
{
    auto maxima = find_maxima(group);
    if (maxima.empty()) return;

    if (maxima.size() == 1 ||
        group.size() >= SPLIT_MAX_HITS ||
        maxima.size() >= SPLIT_MAX_MAXIMA)
    {
        // single cluster from this group
        auto &seed = hits_[maxima[0]];
        clusters_.emplace_back();
        auto &cl = clusters_.back();
        cl.center = seed;
        cl.flag   = sys_.module(seed.index).flag;

        for (int hi : group)
            cl.add_hit(hits_[hi]);
    }
    else {
        split_hits(maxima, group);
    }
}

std::vector<int> HyCalCluster::find_maxima(const std::vector<int> &group) const
{
    std::vector<int> local_max;
    local_max.reserve(20);
    if (group.empty()) return local_max;

    const int gid = hit_group_id_[group[0]];

    for (int hi : group) {
        auto &hit = hits_[hi];
        if (hit.energy < config_.min_center_energy)
            continue;

        bool is_max = true;
        // include corners when checking for maxima (same as old code).  With
        // multi-pulse modules, only the pulse that joined this group counts
        // — others on the same module belong to a different (later) seed.
        sys_.for_each_neighbor(hit.index, true, [&](int ni) {
            if (!is_max) return;
            for (int hj : mod_to_hits_[ni]) {
                if (hit_group_id_[hj] == gid && hits_[hj].energy > hit.energy) {
                    is_max = false;
                    return;
                }
            }
        });

        if (is_max)
            local_max.push_back(hi);
    }

    return local_max;
}

//=============================================================================
// Hit splitting — distribute shared hits among multiple maxima
//=============================================================================

void HyCalCluster::split_hits(const std::vector<int> &maxima,
                               const std::vector<int> &group)
{
    SplitContainer split;  // ~4KB on stack, safe per-call

    int nmax  = static_cast<int>(maxima.size());
    int nhits = static_cast<int>(group.size());

    // initialize fractions from profile
    for (int i = 0; i < nmax; ++i) {
        auto &center = hits_[maxima[i]];
        for (int j = 0; j < nhits; ++j) {
            auto &hit = hits_[group[j]];
            split.frac[j][i] = get_profile_frac(center, hit) * center.energy;
        }
    }

    // iterative refinement
    eval_fraction(maxima, group, split);

    // create clusters from final fractions
    for (int i = 0; i < nmax; ++i) {
        clusters_.emplace_back();
        auto &cl = clusters_.back();
        cl.center = hits_[maxima[i]];
        cl.flag   = sys_.module(cl.center.index).flag;

        for (int j = 0; j < nhits; ++j) {
            if (split.frac[j][i] == 0.f) continue;

            float nf = split.norm_frac(i, j);
            if (nf < config_.least_split) {
                split.total[j] -= split.frac[j][i];
                continue;
            }

            ModuleHit new_hit = hits_[group[j]];
            new_hit.energy *= nf;
            cl.add_hit(new_hit);

            if (new_hit.index == cl.center.index)
                cl.center.energy = new_hit.energy;

            set_bit(cl.flag, kSplit);
        }
    }
}

void HyCalCluster::eval_fraction(const std::vector<int> &maxima,
                                  const std::vector<int> &group,
                                  SplitContainer &split) const
{
    int nmax  = static_cast<int>(maxima.size());
    int nhits = static_cast<int>(group.size());

    struct BaseHit { float x, y, E; };
    BaseHit temp[POS_RECON_HITS];

    int iters = config_.split_iter;
    while (iters-- > 0) {
        split.sum_frac(nhits, nmax);

        for (int i = 0; i < nmax; ++i) {
            auto &center_hit = hits_[maxima[i]];
            const auto &center_mod = sys_.module(center_hit.index);

            // gather 3x3 neighbors for position reconstruction
            float tot_E = center_hit.energy;
            int count = 0;

            for (int j = 0; j < nhits; ++j) {
                auto &hit = hits_[group[j]];
                if (hit.index == center_hit.index || split.frac[j][i] == 0.f)
                    continue;

                // check if within 3x3 using pre-computed neighbors
                const auto &hit_mod = sys_.module(hit.index);
                double dx, dy;
                sys_.qdist(center_mod, hit_mod, dx, dy);

                if (std::abs(dx) < 1.01 && std::abs(dy) < 1.01 && count < POS_RECON_HITS) {
                    float frac_E = hit.energy * split.norm_frac(i, j);
                    temp[count] = {static_cast<float>(dx), static_cast<float>(dy), frac_E};
                    tot_E += frac_E;
                    count++;
                }
            }

            // reconstruct position (log-weighted)
            float wx = 0.f, wy = 0.f;
            float wtot = get_weight(center_hit.energy, tot_E);

            for (int k = 0; k < count; ++k) {
                float w = get_weight(temp[k].E, tot_E);
                if (w > 0.f) {
                    wx += temp[k].x * w;
                    wy += temp[k].y * w;
                    wtot += w;
                }
            }

            float cx = center_mod.x, cy = center_mod.y;
            if (wtot > 0.f) {
                cx += (wx / wtot) * center_mod.size_x;
                cy += (wy / wtot) * center_mod.size_y;
            }

            // update fractions with new center position
            for (int j = 0; j < nhits; ++j) {
                auto &hit = hits_[group[j]];
                split.frac[j][i] = get_profile_frac_at(cx, cy, tot_E, hit) * tot_E;
            }
        }
    }
    split.sum_frac(nhits, nmax);
}

//=============================================================================
// Position reconstruction — log-weighted centroid
//=============================================================================

ClusterHit HyCalCluster::reconstruct_pos(const ModuleCluster &cl) const
{
    const auto &center_mod = sys_.module(cl.center.index);

    struct BaseHit { float x, y, E; };
    BaseHit temp[POS_RECON_HITS];
    int count = 0;

    // gather 3x3 neighbors
    for (auto &hit : cl.hits) {
        if (hit.index == cl.center.index) continue;
        if (count >= POS_RECON_HITS) break;

        double dx, dy;
        sys_.qdist(center_mod, sys_.module(hit.index), dx, dy);
        if (std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
            temp[count++] = {static_cast<float>(dx), static_cast<float>(dy), hit.energy};
        }
    }

    // total energy
    float tot_E = cl.energy;

    // weighted position
    float wx = 0.f, wy = 0.f;
    float wtot = get_weight(cl.center.energy, tot_E);
    int npos = (wtot > 0.f) ? 1 : 0;

    for (int i = 0; i < count; ++i) {
        float w = get_weight(temp[i].E, tot_E);
        if (w > 0.f) {
            wx += temp[i].x * w;
            wy += temp[i].y * w;
            wtot += w;
            npos++;
        }
    }

    ClusterHit result;
    result.center_id = center_mod.id;
    result.energy    = cl.energy;
    result.time      = cl.center.time;
    result.nblocks   = static_cast<int>(cl.hits.size());
    result.flag      = cl.flag;
    result.linear_corr = 1.f;

    if (config_.non_linear_corr) {
        // 1/linear_corr = E_rec/E_exp
        // = 1 + nl1*(E_rec-E_base)/1000 + nl2*((E_rec-E_base)/1000)^2
        const float nl1 = sys_.GetCalibNonLinearity1(center_mod.id);
        const float nl2 = sys_.GetCalibNonLinearity2(center_mod.id);
        const float base_energy = sys_.GetCalibBaseEnergy(center_mod.id);
        const float delta_gev = (cl.energy - base_energy) / 1000.f;
        result.linear_corr = 1.f / (1.f + nl1 * delta_gev
                                          + nl2 * delta_gev * delta_gev);
        if (cl.energy > 3800.f || result.linear_corr < 0.7f || result.linear_corr > 1.3f)
            result.linear_corr = 1.f;
    }

    result.energy = cl.energy * result.linear_corr;

    if (wtot > 0.f) {
        result.x = center_mod.x + (wx / wtot) * center_mod.size_x;
        result.y = center_mod.y + (wy / wtot) * center_mod.size_y;
    } else {
        result.x = center_mod.x;
        result.y = center_mod.y;
    }
    result.npos = npos;

    return result;
}

float HyCalCluster::get_weight(float E, float E_total) const
{
    if (E_total <= 0.f) return 0.f;
    float w = config_.log_weight_thres + std::log(E / E_total);
    return (w > 0.f) ? w : 0.f;
}

//=============================================================================
// Profile helpers
//=============================================================================

float HyCalCluster::get_profile_frac(const ModuleHit &center, const ModuleHit &hit) const
{
    const auto &m1 = sys_.module(center.index);
    const auto &m2 = sys_.module(hit.index);
    double dx, dy;
    sys_.qdist(m1, m2, dx, dy);
    float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    // center module contains ~78% of energy, scale up to estimate total
    return profile_->GetFraction(m1.type, dist, center.energy / 0.78f);
}

float HyCalCluster::get_profile_frac_at(float cx, float cy, float cE,
                                          const ModuleHit &hit) const
{
    const auto &m = sys_.module(hit.index);
    int sid = sys_.get_sector_id(cx, cy);
    double dx, dy;
    sys_.qdist(cx, cy, sid, m.x, m.y, m.sector, dx, dy);
    float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    ModuleType type = sys_.sector_info(sid).mtype;
    return profile_->GetFraction(type, dist, cE);
}

//=============================================================================
// CollectNeighborTiming — emit (seed, neighbour, dt) rows without applying any
// timing cut, for picking the production value of seed_time_window from real
// data.
//
// Seed selection mirrors group_hits(): largest-energy pulse passing
// min_center_energy that hasn't already been claimed AS A SEED in this
// scan.  After selecting a seed, every other pulse within ±max_qdist on
// each axis is emitted — neighbour pulses themselves are NOT consumed,
// so a pulse may appear as a neighbour of more than one seed (intentional:
// the study is about the dt landscape, not about cluster assignment).
//
// Const-friendly: works on local mod_to_hits / consumed scratch so it
// doesn't disturb the per-event state used by FormClusters().
//=============================================================================
void HyCalCluster::CollectNeighborTiming(std::vector<SeedNeighborTiming> &out,
                                          double max_qdist) const
{
    out.clear();
    if (hits_.empty()) return;

    std::vector<std::vector<int>> mod_hits(sys_.module_count());
    for (int i = 0; i < static_cast<int>(hits_.size()); ++i)
        mod_hits[hits_[i].index].push_back(i);

    std::vector<int> order(hits_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return hits_[a].energy > hits_[b].energy; });

    // A seed is "claimed" only against its own (module, time) so a single
    // physics shower with several pulses doesn't double-seed.  Use a small
    // claim window matching min_center_energy regions (5 ns is shorter than
    // any reasonable pulse separation) to keep the seed list compact.
    constexpr float SEED_CLAIM_WINDOW_NS = 5.f;

    std::vector<bool> seed_claimed(hits_.size(), false);

    for (int seed : order) {
        if (seed_claimed[seed]) continue;
        if (hits_[seed].energy < config_.min_center_energy) break;

        // Claim this seed and any pulses on the SAME module within the
        // tight seed-claim window so we don't emit multiple seeds for the
        // same physics event in the same module.
        const float st = hits_[seed].time;
        const int   sm = hits_[seed].index;
        for (int k : mod_hits[sm])
            if (std::fabs(hits_[k].time - st) <= SEED_CLAIM_WINDOW_NS)
                seed_claimed[k] = true;

        const auto &smod = sys_.module(sm);
        for (int k = 0; k < static_cast<int>(hits_.size()); ++k) {
            if (k == seed) continue;
            const auto &nmod = sys_.module(hits_[k].index);
            double dx, dy;
            sys_.qdist(smod, nmod, dx, dy);
            if (std::fabs(dx) > max_qdist || std::fabs(dy) > max_qdist) continue;

            SeedNeighborTiming s;
            s.seed_module     = sm;
            s.neighbor_module = hits_[k].index;
            s.seed_time       = st;
            s.neighbor_time   = hits_[k].time;
            s.dt              = hits_[k].time - st;
            s.seed_energy     = hits_[seed].energy;
            s.neighbor_energy = hits_[k].energy;
            s.dx_q            = dx;
            s.dy_q            = dy;
            out.push_back(s);
        }
    }
}

} // namespace fdec
