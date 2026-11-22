#ifndef XGALOIS_ALGEBRA_SCHREIER_SIMS_HPP_
#define XGALOIS_ALGEBRA_SCHREIER_SIMS_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/algebra/permutation.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace algebra {

namespace detail {

inline uint32_t FirstMovedPoint(const std::vector<uint32_t>& img) noexcept {
  for (uint32_t i = 0; i < static_cast<uint32_t>(img.size()); ++i) {
    if (img[i] != i) {
      return i;
    }
  }
  return static_cast<uint32_t>(img.size());
}

inline bool IsIdentityMapping(const std::vector<uint32_t>& img) noexcept {
  for (uint32_t i = 0; i < static_cast<uint32_t>(img.size()); ++i) {
    if (img[i] != i) {
      return false;
    }
  }
  return true;
}

inline std::vector<uint32_t> InvertMapping(const std::vector<uint32_t>& img) {
  std::vector<uint32_t> inv(img.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(img.size()); ++i) {
    inv[img[i]] = i;
  }
  return inv;
}

// Jerrum's filter: maintains an equivalent generating set of a permutation
// subgroup on {0, ..., degree - 1} whose size is strictly < degree by keeping
// a forest whose edges are {first_moved(g), g(first_moved(g))}.
class JerrumFilter {
 public:
  explicit JerrumFilter(uint32_t degree = 0) : degree_(degree), adj_(degree) {}

  // Inserts permutation `g_img` (mapping of length degree_) into the filter.
  // Returns true if the stored generating set changed.
  bool Insert(std::vector<uint32_t> g_img) {
    bool modified = false;
    std::vector<int32_t> parent_node(degree_, -1);
    std::vector<int32_t> parent_edge(degree_, -1);
    std::vector<uint32_t> bfs_queue;
    bfs_queue.reserve(degree_);

    while (!IsIdentityMapping(g_img)) {
      const uint32_t a = FirstMovedPoint(g_img);
      const uint32_t b = g_img[a];  // Always b > a.

      // Search for a path from `a` to `b` in the current forest.
      std::fill(parent_node.begin(), parent_node.end(), -1);
      std::fill(parent_edge.begin(), parent_edge.end(), -1);
      bfs_queue.clear();
      parent_node[a] = static_cast<int32_t>(a);
      bfs_queue.push_back(a);

      bool found = false;
      for (std::size_t q_idx = 0; q_idx < bfs_queue.size() && !found; ++q_idx) {
        const uint32_t u = bfs_queue[q_idx];
        for (const auto& [v, edge_id] : adj_[u]) {
          if (!edges_[edge_id].active) {
            continue;
          }
          if (parent_node[v] == -1) {
            parent_node[v] = static_cast<int32_t>(u);
            parent_edge[v] = static_cast<int32_t>(edge_id);
            if (v == b) {
              found = true;
              break;
            }
            bfs_queue.push_back(v);
          }
        }
      }

      if (!found) {
        AddEdge(a, b, std::move(g_img));
        return true;
      }

      // Reconstruct the path p[0] = a, p[1], ..., p[t] = b and forest edge ids
      // along the path.
      std::vector<uint32_t> path_nodes;
      std::vector<int32_t> path_edges;
      {
        uint32_t cur = b;
        while (cur != a) {
          path_nodes.push_back(cur);
          path_edges.push_back(parent_edge[cur]);
          cur = static_cast<uint32_t>(parent_node[cur]);
        }
        path_nodes.push_back(a);
        std::reverse(path_nodes.begin(), path_nodes.end());
        std::reverse(path_edges.begin(), path_edges.end());
      }
      // Close the cycle by adding the new edge from b back to a (id = -1).
      path_edges.push_back(-1);

      const std::size_t cycle_len = path_nodes.size();
      std::size_t min_pos = 0;
      for (std::size_t idx = 1; idx < cycle_len; ++idx) {
        if (path_nodes[idx] < path_nodes[min_pos]) {
          min_pos = idx;
        }
      }
      const uint32_t v_min = path_nodes[min_pos];

      // Compute the cycle product w starting at v_min so w fixes all points
      // <= v_min: w = factor_{cycle_len-1} * ... * factor_0 (left-action).
      const std::vector<uint32_t> g_inv = InvertMapping(g_img);
      std::vector<uint32_t> w(degree_);
      std::iota(w.begin(), w.end(), uint32_t{0});

      for (std::size_t step = 0; step < cycle_len; ++step) {
        const std::size_t cur_idx = (min_pos + step) % cycle_len;
        const std::size_t next_idx = (cur_idx + 1) % cycle_len;
        const uint32_t u = path_nodes[cur_idx];
        const uint32_t v = path_nodes[next_idx];
        const int32_t edge_id = path_edges[cur_idx];

        const std::vector<uint32_t>* factor = nullptr;
        if (edge_id == -1) {
          // New edge connects a < b with g_img(a) = b.
          factor = (u < v) ? &g_img : &g_inv;
        } else {
          const Edge& e = edges_[static_cast<std::size_t>(edge_id)];
          factor = (u < v) ? &e.img : &e.inv;
        }
        for (uint32_t x = 0; x < degree_; ++x) {
          w[x] = (*factor)[w[x]];
        }
      }

      if (v_min < a) {
        // Neither endpoint of the new edge (a, b) is v_min, so
        // path_edges[min_pos] is guaranteed to be a forest edge incident to
        // v_min. Replace it with g.
        const std::size_t drop_edge =
            static_cast<std::size_t>(path_edges[min_pos]);
        RemoveEdge(drop_edge);
        AddEdge(a, b, std::move(g_img));
        modified = true;
      }

      g_img = std::move(w);
    }

    return modified;
  }

  std::vector<std::vector<uint32_t>> ActiveMappings() const {
    std::vector<std::vector<uint32_t>> out;
    for (const Edge& e : edges_) {
      if (e.active) {
        out.push_back(e.img);
      }
    }
    return out;
  }

  std::vector<std::vector<uint32_t>> ActiveInverses() const {
    std::vector<std::vector<uint32_t>> out;
    for (const Edge& e : edges_) {
      if (e.active) {
        out.push_back(e.inv);
      }
    }
    return out;
  }

  std::size_t Size() const noexcept { return active_count_; }

 private:
  struct Edge {
    uint32_t u = 0;
    uint32_t v = 0;
    std::vector<uint32_t> img;
    std::vector<uint32_t> inv;
    bool active = true;
  };

  void AddEdge(uint32_t u, uint32_t v, std::vector<uint32_t> img) {
    const std::size_t id = edges_.size();
    std::vector<uint32_t> inv = InvertMapping(img);
    edges_.push_back(Edge{u, v, std::move(img), std::move(inv), true});
    adj_[u].emplace_back(v, id);
    adj_[v].emplace_back(u, id);
    ++active_count_;
  }

  void RemoveEdge(std::size_t id) {
    Edge& e = edges_[id];
    if (!e.active) {
      return;
    }
    e.active = false;
    --active_count_;
    auto erase_from = [&](uint32_t node) {
      auto& vec = adj_[node];
      vec.erase(
          std::remove_if(vec.begin(), vec.end(),
                         [id](const auto& item) { return item.second == id; }),
          vec.end());
    };
    erase_from(e.u);
    erase_from(e.v);
  }

  uint32_t degree_ = 0;
  std::vector<Edge> edges_;
  std::vector<std::vector<std::pair<uint32_t, std::size_t>>> adj_;
  std::size_t active_count_ = 0;
};

}  // namespace detail

// Base and Strong Generating Set (BSGS) stabilizer chain for a finite
// permutation group G <= S_n using left-action composition (a * b)(x) =
// a(b(x)).
class BsgsChain {
 public:
  BsgsChain() = default;

  explicit BsgsChain(uint32_t degree) : degree_(degree) {}

  uint32_t Degree() const noexcept { return degree_; }

  const std::vector<uint32_t>& Base() const noexcept { return base_; }

  std::vector<Permutation> StrongGenerators() const {
    std::unordered_set<Permutation> seen;
    std::vector<Permutation> out;
    for (const Level& lvl : levels_) {
      for (const auto& img : lvl.gens) {
        Permutation p = Permutation::FromMapping(img);
        if (seen.insert(p).second) {
          out.push_back(std::move(p));
        }
      }
    }
    return out;
  }

  std::vector<uint32_t> OrbitSizes() const {
    std::vector<uint32_t> sizes;
    sizes.reserve(levels_.size());
    for (const Level& lvl : levels_) {
      sizes.push_back(static_cast<uint32_t>(lvl.orbit.size()));
    }
    return sizes;
  }

  BigInt Order() const {
    BigInt order(1);
    for (const Level& lvl : levels_) {
      order *= BigInt::FromUint64(static_cast<uint64_t>(lvl.orbit.size()));
    }
    return order;
  }

  // Sifts permutation `g` through the stabilizer chain.
  // Returns the identity permutation if and only if g is in the group.
  Permutation Sift(const Permutation& g) const {
    if (g.Degree() != degree_) {
      throw std::invalid_argument(
          "BsgsChain::Sift: permutation degree does not match BSGS degree");
    }
    std::vector<uint32_t> cur = g.Mapping();
    SiftInPlace(cur, 0);
    return Permutation::FromMapping(std::move(cur));
  }

  bool Contains(const Permutation& g) const {
    if (g.Degree() != degree_) {
      return false;
    }
    std::vector<uint32_t> cur = g.Mapping();
    const std::size_t drop = SiftInPlace(cur, 0);
    return drop == levels_.size() && detail::IsIdentityMapping(cur);
  }

  // Generates a uniformly random element of G as u_0 * u_1 * ... * u_{k-1}.
  Permutation RandomElement(std::mt19937_64& rng) const {
    std::vector<uint32_t> acc(degree_);
    std::iota(acc.begin(), acc.end(), uint32_t{0});
    std::vector<uint32_t> u_buf(degree_);
    std::vector<uint32_t> next_buf(degree_);

    for (std::size_t i = 0; i < levels_.size(); ++i) {
      const Level& lvl = levels_[i];
      std::uniform_int_distribution<std::size_t> dist(0, lvl.orbit.size() - 1);
      const uint32_t gamma = lvl.orbit[dist(rng)];
      if (gamma == lvl.beta) {
        continue;
      }
      ComputeCosetRep(lvl, gamma, u_buf);
      // acc = acc * u_buf (left-action).
      for (uint32_t x = 0; x < degree_; ++x) {
        next_buf[x] = acc[u_buf[x]];
      }
      acc.swap(next_buf);
    }
    return Permutation::FromMapping(std::move(acc));
  }

  // Enumerates all elements of G via the Cartesian product of transversals.
  // Throws std::length_error if |G| exceeds `max_elements`.
  std::vector<Permutation> EnumerateElements(
      std::size_t max_elements = 100000) const {
    const BigInt ord = Order();
    if (ord > BigInt::FromUint64(static_cast<uint64_t>(max_elements))) {
      throw std::length_error(
          "BsgsChain::EnumerateElements: group order exceeds max_elements");
    }
    std::vector<Permutation> result;
    result.reserve(static_cast<std::size_t>(ord.ToUint64()));

    std::vector<uint32_t> id(degree_);
    std::iota(id.begin(), id.end(), uint32_t{0});
    EnumerateRecursive(0, id, result);
    return result;
  }

 private:
  struct Level {
    uint32_t beta = 0;
    detail::JerrumFilter filter;
    std::vector<std::vector<uint32_t>> gens;
    std::vector<std::vector<uint32_t>> gen_invs;
    std::vector<uint32_t> orbit;
    std::vector<int32_t> parent;
    std::vector<int32_t> gen_idx;
    std::size_t scan_pos = 0;

    Level(uint32_t degree, uint32_t base_pt)
        : beta(base_pt),
          filter(degree),
          parent(degree, -1),
          gen_idx(degree, -1) {
      orbit.push_back(base_pt);
      parent[base_pt] = static_cast<int32_t>(base_pt);
    }

    void SyncAndRebuildOrbit() {
      gens = filter.ActiveMappings();
      gen_invs = filter.ActiveInverses();
      orbit.clear();
      std::fill(parent.begin(), parent.end(), -1);
      std::fill(gen_idx.begin(), gen_idx.end(), -1);
      orbit.push_back(beta);
      parent[beta] = static_cast<int32_t>(beta);
      for (std::size_t idx = 0; idx < orbit.size(); ++idx) {
        const uint32_t pt = orbit[idx];
        for (std::size_t g_i = 0; g_i < gens.size(); ++g_i) {
          const uint32_t nxt = gens[g_i][pt];
          if (parent[nxt] == -1) {
            parent[nxt] = static_cast<int32_t>(pt);
            gen_idx[nxt] = static_cast<int32_t>(g_i);
            orbit.push_back(nxt);
          }
        }
      }
      scan_pos = 0;
    }
  };

  // Reconstructs coset representative u_gamma mapping lvl.beta -> gamma.
  void ComputeCosetRep(const Level& lvl, uint32_t gamma,
                       std::vector<uint32_t>& out) const {
    out.resize(degree_);
    std::iota(out.begin(), out.end(), uint32_t{0});
    if (gamma == lvl.beta) {
      return;
    }
    // Walk backwards from gamma to lvl.beta; if path is
    // beta = x_0 --(s_1)--> x_1 --(s_2)--> ... --(s_m)--> gamma,
    // then u_gamma = s_m * ... * s_2 * s_1 (left-action).
    std::vector<int32_t> steps;
    uint32_t cur = gamma;
    while (cur != lvl.beta) {
      steps.push_back(lvl.gen_idx[cur]);
      cur = static_cast<uint32_t>(lvl.parent[cur]);
    }
    // Apply from s_1 to s_m on the left: out[x] = s_j[out[x]].
    for (std::size_t k = steps.size(); k-- > 0;) {
      const auto& s = lvl.gens[static_cast<std::size_t>(steps[k])];
      for (uint32_t x = 0; x < degree_; ++x) {
        out[x] = s[out[x]];
      }
    }
  }

  // Sifts `cur` in-place starting from `start_level`.
  // Returns the first level where `cur(beta_i)` is not in the orbit, or
  // `levels_.size()` if `cur` passed all levels.
  std::size_t SiftInPlace(std::vector<uint32_t>& cur,
                          std::size_t start_level) const {
    for (std::size_t i = start_level; i < levels_.size(); ++i) {
      const Level& lvl = levels_[i];
      uint32_t img_beta = cur[lvl.beta];
      if (lvl.parent[img_beta] == -1) {
        return i;
      }
      // Left-multiply cur by u_{img_beta}^{-1} step-by-step along the
      // Schreier vector back to lvl.beta.
      while (img_beta != lvl.beta) {
        const auto& s_inv =
            lvl.gen_invs[static_cast<std::size_t>(lvl.gen_idx[img_beta])];
        for (uint32_t x = 0; x < degree_; ++x) {
          cur[x] = s_inv[cur[x]];
        }
        img_beta = static_cast<uint32_t>(lvl.parent[img_beta]);
      }
    }
    return levels_.size();
  }

  void EnumerateRecursive(std::size_t level_idx,
                          const std::vector<uint32_t>& prefix,
                          std::vector<Permutation>& out) const {
    if (level_idx == levels_.size()) {
      out.push_back(Permutation::FromMapping(prefix));
      return;
    }
    const Level& lvl = levels_[level_idx];
    std::vector<uint32_t> u_buf(degree_);
    std::vector<uint32_t> next_prefix(degree_);
    for (uint32_t gamma : lvl.orbit) {
      if (gamma == lvl.beta) {
        EnumerateRecursive(level_idx + 1, prefix, out);
      } else {
        ComputeCosetRep(lvl, gamma, u_buf);
        for (uint32_t x = 0; x < degree_; ++x) {
          next_prefix[x] = prefix[u_buf[x]];
        }
        EnumerateRecursive(level_idx + 1, next_prefix, out);
      }
    }
  }

  void AppendLevel(uint32_t base_pt) {
    base_.push_back(base_pt);
    levels_.emplace_back(degree_, base_pt);
  }

  friend BsgsChain BuildBsgsDeterministic(
      uint32_t degree, const std::vector<Permutation>& generators);
  friend BsgsChain BuildBsgsRandomized(
      uint32_t degree, const std::vector<Permutation>& generators,
      uint64_t seed);

  uint32_t degree_ = 0;
  std::vector<uint32_t> base_;
  std::vector<Level> levels_;
};

// Deterministic Schreier-Sims algorithm using Schreier vectors and Jerrum's
// filter to keep strong generator counts O(n) per level.
inline BsgsChain BuildBsgsDeterministic(
    uint32_t degree, const std::vector<Permutation>& generators) {
  BsgsChain chain(degree);

  // Validate and collect non-identity generators.
  std::vector<std::vector<uint32_t>> non_id_gens;
  non_id_gens.reserve(generators.size());
  for (const Permutation& g : generators) {
    if (g.Degree() != degree) {
      throw std::invalid_argument(
          "BuildBsgsDeterministic: generator degree mismatch");
    }
    if (!g.IsIdentity()) {
      non_id_gens.push_back(g.Mapping());
    }
  }

  if (non_id_gens.empty()) {
    return chain;
  }

  // Initialize base and initial stabilizer levels from input generators.
  for (const auto& g_img : non_id_gens) {
    std::size_t deepest = 0;
    while (deepest < chain.levels_.size() &&
           g_img[chain.levels_[deepest].beta] == chain.levels_[deepest].beta) {
      ++deepest;
    }
    if (deepest == chain.levels_.size()) {
      chain.AppendLevel(detail::FirstMovedPoint(g_img));
    }
    for (std::size_t l = 0; l <= deepest; ++l) {
      chain.levels_[l].filter.Insert(g_img);
    }
  }
  for (auto& lvl : chain.levels_) {
    lvl.SyncAndRebuildOrbit();
  }

  // Main Schreier-Sims closure loop.
  // Invariant: upon entering level `i`, levels i+1 .. k-1 form a complete BSGS
  // for <S^{(i+1)}>.
  std::vector<uint32_t> u_gamma(degree);
  std::vector<uint32_t> prod(degree);

  int64_t i = static_cast<int64_t>(chain.levels_.size()) - 1;
  while (i >= 0) {
    const std::size_t ui = static_cast<std::size_t>(i);
    bool added_generator = false;
    std::size_t restart_level = ui;

    const std::size_t num_orbit = chain.levels_[ui].orbit.size();
    const std::size_t num_gens = chain.levels_[ui].gens.size();
    const std::size_t total_pairs = num_orbit * num_gens;

    std::size_t cached_orbit_idx = num_orbit;  // invalid sentinel

    while (chain.levels_[ui].scan_pos < total_pairs) {
      const std::size_t pos = chain.levels_[ui].scan_pos++;
      const std::size_t o_idx = pos / num_gens;
      const std::size_t g_idx = pos % num_gens;

      const auto& lvl = chain.levels_[ui];
      const uint32_t gamma = lvl.orbit[o_idx];
      const uint32_t delta = lvl.gens[g_idx][gamma];

      // Skip Schreier tree edges: if s(gamma) == delta is the defining tree
      // edge for delta, then u_delta = s * u_gamma, so the Schreier generator
      // u_delta^{-1} * s * u_gamma is trivially the identity.
      if (lvl.parent[delta] == static_cast<int32_t>(gamma) &&
          lvl.gen_idx[delta] == static_cast<int32_t>(g_idx)) {
        continue;
      }

      if (cached_orbit_idx != o_idx) {
        chain.ComputeCosetRep(lvl, gamma, u_gamma);
        cached_orbit_idx = o_idx;
      }

      // Form prod = s * u_gamma (left-action).
      const auto& s = lvl.gens[g_idx];
      for (uint32_t x = 0; x < degree; ++x) {
        prod[x] = s[u_gamma[x]];
      }

      // Sifting from level `ui` first applies u_delta^{-1} on the left to
      // produce the Schreier generator g_{gamma, s} = u_delta^{-1} * s *
      // u_gamma, then sifts it through levels ui+1 .. k-1.
      const std::size_t drop = chain.SiftInPlace(prod, ui);
      if (drop < chain.levels_.size() || !detail::IsIdentityMapping(prod)) {
        std::size_t target_level = drop;
        if (target_level == chain.levels_.size()) {
          chain.AppendLevel(detail::FirstMovedPoint(prod));
        }
        // Insert the non-trivial residue into levels ui+1 .. target_level.
        for (std::size_t l = ui + 1; l <= target_level; ++l) {
          if (chain.levels_[l].filter.Insert(prod)) {
            chain.levels_[l].SyncAndRebuildOrbit();
          }
        }
        added_generator = true;
        restart_level = target_level;
        break;
      }
    }

    if (added_generator) {
      i = static_cast<int64_t>(restart_level);
    } else {
      --i;
    }
  }

  return chain;
}

// Fast randomized Schreier-Sims algorithm using Product Replacement random
// elements and Jerrum's filter.
inline BsgsChain BuildBsgsRandomized(uint32_t degree,
                                     const std::vector<Permutation>& generators,
                                     uint64_t seed = 42) {
  BsgsChain chain(degree);

  std::vector<std::vector<uint32_t>> non_id_gens;
  non_id_gens.reserve(generators.size());
  for (const Permutation& g : generators) {
    if (g.Degree() != degree) {
      throw std::invalid_argument(
          "BuildBsgsRandomized: generator degree mismatch");
    }
    if (!g.IsIdentity()) {
      non_id_gens.push_back(g.Mapping());
    }
  }

  if (non_id_gens.empty()) {
    return chain;
  }

  auto insert_into_chain = [&](std::vector<uint32_t> elem,
                               bool include_level_zero) -> bool {
    const std::size_t drop = chain.SiftInPlace(elem, 0);
    if (drop == chain.levels_.size() && detail::IsIdentityMapping(elem)) {
      return false;
    }
    std::size_t target = drop;
    if (target == chain.levels_.size()) {
      chain.AppendLevel(detail::FirstMovedPoint(elem));
    }
    const std::size_t start_l = (include_level_zero || target == 0) ? 0 : 1;
    for (std::size_t l = start_l; l <= target; ++l) {
      if (chain.levels_[l].filter.Insert(elem)) {
        chain.levels_[l].SyncAndRebuildOrbit();
      }
    }
    return true;
  };

  // Seed level 0 with all input generators so <S^{(0)}> == G.
  for (const auto& g_img : non_id_gens) {
    std::size_t deepest = 0;
    while (deepest < chain.levels_.size() &&
           g_img[chain.levels_[deepest].beta] == chain.levels_[deepest].beta) {
      ++deepest;
    }
    if (deepest == chain.levels_.size()) {
      chain.AppendLevel(detail::FirstMovedPoint(g_img));
    }
    for (std::size_t l = 0; l <= deepest; ++l) {
      chain.levels_[l].filter.Insert(g_img);
    }
  }
  for (auto& lvl : chain.levels_) {
    lvl.SyncAndRebuildOrbit();
  }

  // Initialize Product Replacement pool for near-uniform random elements of G.
  const std::size_t pool_size =
      std::max<std::size_t>(12, non_id_gens.size() + 4);
  std::vector<std::vector<uint32_t>> pool;
  pool.reserve(pool_size);
  for (std::size_t idx = 0; idx < pool_size; ++idx) {
    pool.push_back(non_id_gens[idx % non_id_gens.size()]);
  }

  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<std::size_t> pool_dist(0, pool_size - 1);
  std::vector<uint32_t> scratch(degree);

  auto step_pool = [&]() -> const std::vector<uint32_t>& {
    const std::size_t a = pool_dist(rng);
    std::size_t b = pool_dist(rng);
    if (b == a) {
      b = (b + 1) % pool_size;
    }
    const uint64_t mode = rng() & 3ULL;
    const auto& pa = pool[a];
    const auto& pb = pool[b];
    if (mode == 0) {
      for (uint32_t x = 0; x < degree; ++x) {
        scratch[x] = pa[pb[x]];
      }
    } else if (mode == 1) {
      for (uint32_t x = 0; x < degree; ++x) {
        scratch[x] = pb[pa[x]];
      }
    } else if (mode == 2) {
      const auto pb_inv = detail::InvertMapping(pb);
      for (uint32_t x = 0; x < degree; ++x) {
        scratch[x] = pa[pb_inv[x]];
      }
    } else {
      const auto pb_inv = detail::InvertMapping(pb);
      for (uint32_t x = 0; x < degree; ++x) {
        scratch[x] = pb_inv[pa[x]];
      }
    }
    pool[a] = scratch;
    return pool[a];
  };

  // Burn-in mixing steps.
  const std::size_t burn_in = std::max<std::size_t>(60, 4 * pool_size);
  for (std::size_t step = 0; step < burn_in; ++step) {
    step_pool();
  }

  // Sift random elements until kRequiredConsecutiveSifts consecutive random
  // elements sift cleanly (failure probability <= 2^{-60}).
  constexpr std::size_t kRequiredConsecutiveSifts = 60;
  std::size_t consecutive = 0;
  std::vector<uint32_t> candidate(degree);

  while (consecutive < kRequiredConsecutiveSifts) {
    const auto& r = step_pool();
    candidate = r;
    if (insert_into_chain(std::move(candidate), /*include_level_zero=*/false)) {
      consecutive = 0;
    } else {
      ++consecutive;
    }
  }

  return chain;
}

}  // namespace algebra

using algebra::BsgsChain;
using algebra::BuildBsgsDeterministic;
using algebra::BuildBsgsRandomized;

}  // namespace xg

#endif  // XGALOIS_ALGEBRA_SCHREIER_SIMS_HPP_
