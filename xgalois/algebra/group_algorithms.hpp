#ifndef XGALOIS_ALGEBRA_GROUP_ALGORITHMS_HPP_
#define XGALOIS_ALGEBRA_GROUP_ALGORITHMS_HPP_

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/group_base.hpp"
#include "xgalois/algebra/snf.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace algebra {

namespace detail {

inline bool IsPrimeU64(uint64_t n) noexcept {
  if (n < 2) {
    return false;
  }
  if (n == 2 || n == 3) {
    return true;
  }
  if (n % 2 == 0 || n % 3 == 0) {
    return false;
  }
  for (uint64_t d = 5; d <= n / d; d += 6) {
    if (n % d == 0 || n % (d + 2) == 0) {
      return false;
    }
  }
  return true;
}

template <typename T>
struct GenericValueHash {
  static constexpr bool kSupported =
      requires(const T& v) { std::hash<T>{}(v); };

  std::size_t operator()(const T& v) const noexcept {
    if constexpr (kSupported) {
      return std::hash<T>{}(v);
    } else {
      return 0;
    }
  }
};

template <typename A, typename B>
struct GenericValueHash<std::pair<A, B>> {
  static constexpr bool kSupported =
      GenericValueHash<A>::kSupported && GenericValueHash<B>::kSupported;

  std::size_t operator()(const std::pair<A, B>& p) const noexcept {
    const std::size_t h1 = GenericValueHash<A>{}(p.first);
    const std::size_t h2 = GenericValueHash<B>{}(p.second);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

template <Group G>
inline typename G::element_type CanonicalizeElement(
    const G& group, const typename G::element_type& elem) {
  return group.Op(elem, group.Identity());
}

template <Group G>
struct GroupElemHash {
  const G* group = nullptr;

  std::size_t operator()(const typename G::element_type& elem) const {
    using E = typename G::element_type;
    const E canon = CanonicalizeElement(*group, elem);
    if constexpr (GenericValueHash<E>::kSupported) {
      return GenericValueHash<E>{}(canon);
    } else {
      return std::hash<std::string>{}(GroupElementToString(*group, canon));
    }
  }
};

template <Group G>
struct GroupElemEqual {
  const G* group = nullptr;

  bool operator()(const typename G::element_type& a,
                  const typename G::element_type& b) const {
    return static_cast<bool>(group->Equals(a, b));
  }
};

template <Group G>
using GroupElemIndexMap =
    std::unordered_map<typename G::element_type, std::size_t, GroupElemHash<G>,
                       GroupElemEqual<G>>;

template <Group G>
inline GroupElemIndexMap<G> MakeGroupIndexMap(const G& group) {
  return GroupElemIndexMap<G>(16, GroupElemHash<G>{&group},
                              GroupElemEqual<G>{&group});
}

template <Group G>
struct ClosureResult {
  std::vector<typename G::element_type> elements;
  GroupElemIndexMap<G> index_map;
};

template <Group G>
inline ClosureResult<G> ComputeClosure(
    const G& group, const std::vector<typename G::element_type>& generators,
    std::size_t max_elements = 100000) {
  using E = typename G::element_type;
  ClosureResult<G> res{std::vector<E>{}, MakeGroupIndexMap(group)};

  const E id = group.Identity();
  res.elements.push_back(id);
  res.index_map.emplace(id, 0);

  // Collect deduplicated non-identity generators and their inverses.
  std::vector<E> step_gens;
  auto step_seen = MakeGroupIndexMap(group);
  step_seen.emplace(id, 0);
  for (const E& raw_g : generators) {
    const E g = CanonicalizeElement(group, raw_g);
    if (step_seen.emplace(g, step_gens.size()).second) {
      step_gens.push_back(g);
    }
    const E g_inv = CanonicalizeElement(group, group.Inv(g));
    if (step_seen.emplace(g_inv, step_gens.size()).second) {
      step_gens.push_back(g_inv);
    }
  }

  for (std::size_t q_idx = 0; q_idx < res.elements.size(); ++q_idx) {
    const E cur = res.elements[q_idx];
    for (const E& s : step_gens) {
      const E nxt = CanonicalizeElement(group, group.Op(cur, s));
      if (res.index_map.find(nxt) == res.index_map.end()) {
        if (res.elements.size() >= max_elements) {
          throw std::length_error(
              "ComputeClosure: group/subgroup order exceeds limit");
        }
        const std::size_t idx = res.elements.size();
        res.elements.push_back(nxt);
        res.index_map.emplace(res.elements.back(), idx);
      }
    }
  }
  return res;
}

template <Group G>
inline ClosureResult<G> EnumerateGroup(const G& group,
                                       std::size_t max_elements = 100000) {
  static_assert(
      requires(const G& g) { g.Generators(); },
      "Group algorithm requires G::Generators()");
  return ComputeClosure(group, group.Generators(), max_elements);
}

// Lightweight integer-indexed group representation for fast exact Sylow
// recursion (quotients and centralizers).
struct IndexedGroup {
  std::size_t n = 0;  // Element 0 is always the identity.
  std::vector<std::vector<uint32_t>> mul;
  std::vector<uint32_t> inv;

  uint32_t Pow(uint32_t a, uint64_t exp) const noexcept {
    uint32_t acc = 0;
    uint32_t base = a;
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = mul[acc][base];
      }
      exp >>= 1;
      if (exp > 0) {
        base = mul[base][base];
      }
    }
    return acc;
  }

  uint64_t ElementOrder(uint32_t a) const noexcept {
    uint32_t cur = a;
    uint64_t ord = 1;
    while (cur != 0) {
      cur = mul[cur][a];
      ++ord;
    }
    return ord;
  }
};

// Exact recursive Sylow p-subgroup finder following the constructive proof of
// Sylow's First Theorem (center quotient + class equation centralizer descent).
// Precondition: p is prime and p_pow = p^v divides ig.n with gcd(p, ig.n/p_pow)
// == 1 and p_pow > 1.
inline std::vector<uint32_t> FindSylowInIndexedGroup(const IndexedGroup& ig,
                                                     uint64_t p,
                                                     uint64_t p_pow) {
  if (ig.n == p_pow) {
    std::vector<uint32_t> all(ig.n);
    for (std::size_t i = 0; i < ig.n; ++i) {
      all[i] = static_cast<uint32_t>(i);
    }
    return all;
  }

  // Fast base case: p_pow == p (cyclic subgroup of prime order p).
  if (p_pow == p) {
    for (uint32_t x = 1; x < static_cast<uint32_t>(ig.n); ++x) {
      const uint64_t ord = ig.ElementOrder(x);
      if (ord % p == 0) {
        const uint32_t z = ig.Pow(x, ord / p);
        std::vector<uint32_t> sub;
        sub.reserve(static_cast<std::size_t>(p));
        uint32_t cur = 0;
        for (uint64_t k = 0; k < p; ++k) {
          sub.push_back(cur);
          cur = ig.mul[cur][z];
        }
        return sub;
      }
    }
  }

  // Compute the center Z(H).
  std::vector<uint32_t> center;
  std::vector<uint8_t> in_center(ig.n, 0);
  for (uint32_t x = 0; x < static_cast<uint32_t>(ig.n); ++x) {
    bool commutes = true;
    for (uint32_t y = 1; y < static_cast<uint32_t>(ig.n); ++y) {
      if (ig.mul[x][y] != ig.mul[y][x]) {
        commutes = false;
        break;
      }
    }
    if (commutes) {
      center.push_back(x);
      in_center[x] = 1;
    }
  }

  // Case 1: p divides |Z(H)|. Quotient by a central subgroup N of order p and
  // pull back a Sylow p-subgroup of H / N.
  if (center.size() % p == 0) {
    uint32_t z = 0;
    for (uint32_t c : center) {
      if (c == 0) {
        continue;
      }
      const uint64_t ord = ig.ElementOrder(c);
      if (ord % p == 0) {
        z = ig.Pow(c, ord / p);
        break;
      }
    }
    // Build N = <z>.
    std::vector<uint32_t> normal_elems;
    normal_elems.reserve(static_cast<std::size_t>(p));
    uint32_t cur = 0;
    for (uint64_t k = 0; k < p; ++k) {
      normal_elems.push_back(cur);
      cur = ig.mul[cur][z];
    }

    // Partition H into cosets of N.
    const std::size_t q_n = ig.n / static_cast<std::size_t>(p);
    std::vector<int32_t> coset_of(ig.n, -1);
    std::vector<uint32_t> coset_rep;
    coset_rep.reserve(q_n);

    for (uint32_t u = 0; u < static_cast<uint32_t>(ig.n); ++u) {
      if (coset_of[u] != -1) {
        continue;
      }
      const int32_t cid = static_cast<int32_t>(coset_rep.size());
      coset_rep.push_back(u);
      for (uint32_t n_el : normal_elems) {
        coset_of[ig.mul[u][n_el]] = cid;
      }
    }

    IndexedGroup q_ig;
    q_ig.n = q_n;
    q_ig.mul.assign(q_n, std::vector<uint32_t>(q_n, 0));
    q_ig.inv.resize(q_n, 0);
    for (std::size_t i = 0; i < q_n; ++i) {
      q_ig.inv[i] = static_cast<uint32_t>(coset_of[ig.inv[coset_rep[i]]]);
      for (std::size_t j = 0; j < q_n; ++j) {
        q_ig.mul[i][j] =
            static_cast<uint32_t>(coset_of[ig.mul[coset_rep[i]][coset_rep[j]]]);
      }
    }

    const std::vector<uint32_t> bar_p =
        FindSylowInIndexedGroup(q_ig, p, p_pow / p);
    std::vector<uint8_t> in_bar_p(q_n, 0);
    for (uint32_t c : bar_p) {
      in_bar_p[c] = 1;
    }

    std::vector<uint32_t> pullback;
    pullback.reserve(static_cast<std::size_t>(p_pow));
    for (uint32_t u = 0; u < static_cast<uint32_t>(ig.n); ++u) {
      if (in_bar_p[static_cast<std::size_t>(coset_of[u])] != 0) {
        pullback.push_back(u);
      }
    }
    return pullback;
  }

  // Case 2: p does not divide |Z(H)|. By the class equation there exists a
  // non-central element x whose centralizer C_H(x) is a proper subgroup of H
  // with p^v dividing |C_H(x)|.
  std::vector<uint8_t> visited(ig.n, 0);
  for (uint32_t x = 1; x < static_cast<uint32_t>(ig.n); ++x) {
    if (in_center[x] != 0 || visited[x] != 0) {
      continue;
    }
    std::vector<uint32_t> centralizer;
    for (uint32_t y = 0; y < static_cast<uint32_t>(ig.n); ++y) {
      if (ig.mul[x][y] == ig.mul[y][x]) {
        centralizer.push_back(y);
      }
    }
    if (centralizer.size() < ig.n && centralizer.size() % p_pow == 0) {
      // Build IndexedGroup for C_H(x) (0 is always in centralizer at index 0).
      const std::size_t m = centralizer.size();
      std::vector<int32_t> orig_to_sub(ig.n, -1);
      for (std::size_t idx = 0; idx < m; ++idx) {
        orig_to_sub[centralizer[idx]] = static_cast<int32_t>(idx);
      }
      IndexedGroup sub_ig;
      sub_ig.n = m;
      sub_ig.mul.assign(m, std::vector<uint32_t>(m, 0));
      sub_ig.inv.resize(m, 0);
      for (std::size_t i = 0; i < m; ++i) {
        sub_ig.inv[i] =
            static_cast<uint32_t>(orig_to_sub[ig.inv[centralizer[i]]]);
        for (std::size_t j = 0; j < m; ++j) {
          sub_ig.mul[i][j] = static_cast<uint32_t>(
              orig_to_sub[ig.mul[centralizer[i]][centralizer[j]]]);
        }
      }
      const std::vector<uint32_t> sub_p =
          FindSylowInIndexedGroup(sub_ig, p, p_pow);
      std::vector<uint32_t> lifted;
      lifted.reserve(sub_p.size());
      for (uint32_t idx : sub_p) {
        lifted.push_back(centralizer[idx]);
      }
      return lifted;
    }
    // Mark entire conjugacy class of x as visited.
    for (uint32_t g_el = 0; g_el < static_cast<uint32_t>(ig.n); ++g_el) {
      const uint32_t conj = ig.mul[ig.mul[g_el][x]][ig.inv[g_el]];
      visited[conj] = 1;
    }
  }

  throw std::logic_error(
      "FindSylowInIndexedGroup: internal invariant violated");
}

}  // namespace detail

// ============================================================================
// ElementOrder
// ============================================================================

template <Group G>
inline uint64_t ElementOrder(const G& group,
                             const typename G::element_type& x) {
  if constexpr (requires(const G& g, const typename G::element_type& e) {
                  { g.ElementOrder(e) } -> std::convertible_to<uint64_t>;
                }) {
    return static_cast<uint64_t>(group.ElementOrder(x));
  } else {
    if (group.IsIdentity(x)) {
      return 1;
    }
    typename G::element_type cur = x;
    uint64_t ord = 1;
    while (!group.IsIdentity(cur)) {
      cur = group.Op(cur, x);
      ++ord;
    }
    return ord;
  }
}

// ============================================================================
// Subgroup<G>
// ============================================================================

template <Group G>
class Subgroup {
 public:
  using group_type = G;
  using element_type = typename G::element_type;

  Subgroup(std::shared_ptr<const G> group, std::vector<element_type> generators)
      : group_(std::move(group)),
        generators_(std::move(generators)),
        index_map_(detail::MakeGroupIndexMap(RequireGroupRef(group_))) {
    auto closure = detail::ComputeClosure(*group_, generators_);
    elements_ = std::move(closure.elements);
    index_map_ = std::move(closure.index_map);
    // Canonicalize stored generators and keep at least identity if empty.
    for (auto& g : generators_) {
      g = detail::CanonicalizeElement(*group_, g);
    }
    if (generators_.empty()) {
      generators_.push_back(group_->Identity());
    }
  }

  static Subgroup FromElements(std::shared_ptr<const G> group,
                               std::vector<element_type> elements,
                               std::vector<element_type> generators = {}) {
    if (!group) {
      throw std::invalid_argument("Subgroup: group pointer cannot be null");
    }
    Subgroup sub(group);
    sub.elements_.reserve(elements.size());
    for (const auto& raw : elements) {
      const element_type canon = detail::CanonicalizeElement(*group, raw);
      if (sub.index_map_.find(canon) == sub.index_map_.end()) {
        const std::size_t idx = sub.elements_.size();
        sub.elements_.push_back(canon);
        sub.index_map_.emplace(sub.elements_.back(), idx);
      }
    }
    if (sub.elements_.empty()) {
      const element_type id = group->Identity();
      sub.elements_.push_back(id);
      sub.index_map_.emplace(id, 0);
    }
    if (generators.empty()) {
      sub.generators_ = sub.ExtractCompactGenerators();
    } else {
      sub.generators_.reserve(generators.size());
      for (const auto& g : generators) {
        sub.generators_.push_back(detail::CanonicalizeElement(*group, g));
      }
    }
    return sub;
  }

  std::shared_ptr<const G> AmbientGroup() const noexcept { return group_; }

  std::shared_ptr<const G> Group() const noexcept { return group_; }

  const std::vector<element_type>& Elements() const noexcept {
    return elements_;
  }

  uint64_t Order() const noexcept {
    return static_cast<uint64_t>(elements_.size());
  }

  bool Contains(const element_type& e) const {
    const element_type canon = detail::CanonicalizeElement(*group_, e);
    return index_map_.find(canon) != index_map_.end();
  }

  const std::vector<element_type>& Generators() const noexcept {
    return generators_;
  }

  // Index [G : H] = |G| / |H| by Lagrange's theorem.
  BigInt Index() const {
    if constexpr (FiniteGroup<G>) {
      const BigInt ambient_order = group_->Order();
      const BigInt sub_order = BigInt::FromUint64(Order());
      return ambient_order / sub_order;
    } else {
      const auto all = detail::EnumerateGroup(*group_);
      return BigInt::FromUint64(
          static_cast<uint64_t>(all.elements.size() / elements_.size()));
    }
  }

  bool IsTrivial() const noexcept { return elements_.size() == 1; }

  bool operator==(const Subgroup<G>& other) const {
    if (Order() != other.Order()) {
      return false;
    }
    for (const auto& e : elements_) {
      if (!other.Contains(e)) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const Subgroup<G>& other) const { return !(*this == other); }

 private:
  static const G& RequireGroupRef(const std::shared_ptr<const G>& group) {
    if (!group) {
      throw std::invalid_argument("Subgroup: group pointer cannot be null");
    }
    return *group;
  }

  explicit Subgroup(std::shared_ptr<const G> group)
      : group_(std::move(group)),
        index_map_(detail::MakeGroupIndexMap(RequireGroupRef(group_))) {}

  std::vector<element_type> ExtractCompactGenerators() const {
    std::vector<element_type> gens;
    auto covered = detail::ComputeClosure(*group_, gens);
    for (const auto& e : elements_) {
      if (covered.index_map.find(e) == covered.index_map.end()) {
        gens.push_back(e);
        covered = detail::ComputeClosure(*group_, gens);
        if (covered.elements.size() == elements_.size()) {
          break;
        }
      }
    }
    if (gens.empty()) {
      gens.push_back(group_->Identity());
    }
    return gens;
  }

  std::shared_ptr<const G> group_;
  std::vector<element_type> generators_;
  std::vector<element_type> elements_;
  detail::GroupElemIndexMap<G> index_map_;
};

// ============================================================================
// Cosets & Normality
// ============================================================================

template <Group G>
inline std::vector<std::vector<typename G::element_type>> LeftCosets(
    const G& group, const Subgroup<G>& h) {
  using E = typename G::element_type;
  const auto all = detail::EnumerateGroup(group);
  auto visited = detail::MakeGroupIndexMap(group);
  std::vector<std::vector<E>> cosets;
  cosets.reserve(all.elements.size() / h.Elements().size());

  for (const E& g_el : all.elements) {
    if (visited.find(g_el) != visited.end()) {
      continue;
    }
    std::vector<E> coset;
    coset.reserve(h.Elements().size());
    for (const E& h_el : h.Elements()) {
      const E prod = detail::CanonicalizeElement(group, group.Op(g_el, h_el));
      visited.emplace(prod, 0);
      coset.push_back(prod);
    }
    cosets.push_back(std::move(coset));
  }
  return cosets;
}

template <Group G>
inline std::vector<std::vector<typename G::element_type>> RightCosets(
    const G& group, const Subgroup<G>& h) {
  using E = typename G::element_type;
  const auto all = detail::EnumerateGroup(group);
  auto visited = detail::MakeGroupIndexMap(group);
  std::vector<std::vector<E>> cosets;
  cosets.reserve(all.elements.size() / h.Elements().size());

  for (const E& g_el : all.elements) {
    if (visited.find(g_el) != visited.end()) {
      continue;
    }
    std::vector<E> coset;
    coset.reserve(h.Elements().size());
    for (const E& h_el : h.Elements()) {
      const E prod = detail::CanonicalizeElement(group, group.Op(h_el, g_el));
      visited.emplace(prod, 0);
      coset.push_back(prod);
    }
    cosets.push_back(std::move(coset));
  }
  return cosets;
}

// Tests whether H is normal in G: s * y * s^{-1} in H and s^{-1} * y * s in H
// for all generators s of G and all generators y of H.
template <Group G>
inline bool IsNormalSubgroup(const G& group, const Subgroup<G>& h) {
  using E = typename G::element_type;
  const std::vector<E> g_gens = group.Generators();
  const std::vector<E>& h_gens = h.Generators();
  for (const E& s : g_gens) {
    const E s_inv = group.Inv(s);
    for (const E& y : h_gens) {
      const E conj1 = group.Op(group.Op(s, y), s_inv);
      if (!h.Contains(conj1)) {
        return false;
      }
      const E conj2 = group.Op(group.Op(s_inv, y), s);
      if (!h.Contains(conj2)) {
        return false;
      }
    }
  }
  return true;
}

// ============================================================================
// QuotientGroup<G>: G / N with canonical coset representatives
// ============================================================================

template <Group G>
class QuotientGroup {
 public:
  using ambient_group_type = G;
  using element_type = typename G::element_type;

  QuotientGroup(std::shared_ptr<const G> group, Subgroup<G> normal_subgroup)
      : group_(std::move(group)),
        normal_(std::move(normal_subgroup)),
        elem_to_rep_(16, detail::GroupElemHash<G>{&RequireGroupRef(group_)},
                     detail::GroupElemEqual<G>{&RequireGroupRef(group_)}) {
    if (!IsNormalSubgroup(*group_, normal_)) {
      throw std::invalid_argument(
          "QuotientGroup: subgroup must be normal in ambient group");
    }

    const auto cosets = LeftCosets(*group_, normal_);
    const element_type id = group_->Identity();
    coset_reps_.reserve(cosets.size());

    for (const auto& coset : cosets) {
      element_type rep = coset.front();
      if (normal_.Contains(rep)) {
        rep = id;
      } else if constexpr (requires(const element_type& a,
                                    const element_type& b) {
                             { a < b } -> std::convertible_to<bool>;
                           }) {
        for (const auto& x : coset) {
          if (x < rep) {
            rep = x;
          }
        }
      }
      coset_reps_.push_back(rep);
      for (const auto& x : coset) {
        elem_to_rep_.emplace(x, rep);
      }
    }

    // Build deduplicated canonical generators for G / N.
    auto gen_seen = detail::MakeGroupIndexMap(*group_);
    gen_seen.emplace(id, 0);
    for (const auto& s : group_->Generators()) {
      const element_type proj = Canonicalize(s);
      if (gen_seen.emplace(proj, 0).second) {
        generators_.push_back(proj);
      }
    }
    if (generators_.empty()) {
      generators_.push_back(id);
    }
  }

  std::shared_ptr<const G> AmbientGroup() const noexcept { return group_; }

  const Subgroup<G>& NormalSubgroup() const noexcept { return normal_; }

  const std::vector<element_type>& CosetRepresentatives() const noexcept {
    return coset_reps_;
  }

  element_type Canonicalize(const element_type& x) const {
    const element_type canon = detail::CanonicalizeElement(*group_, x);
    const auto it = elem_to_rep_.find(canon);
    if (it == elem_to_rep_.end()) {
      throw std::invalid_argument(
          "QuotientGroup::Canonicalize: element is not in ambient group");
    }
    return it->second;
  }

  element_type Project(const element_type& x) const { return Canonicalize(x); }

  element_type Op(const element_type& a, const element_type& b) const {
    return Canonicalize(group_->Op(a, b));
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return normal_.Contains(group_->Div(a, b));
  }

  element_type Identity() const { return group_->Identity(); }

  bool IsIdentity(const element_type& a) const { return normal_.Contains(a); }

  element_type Inv(const element_type& a) const {
    return Canonicalize(group_->Inv(a));
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return Canonicalize(group_->Div(a, b));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    return Canonicalize(group_->Pow(a, exp));
  }

  BigInt Order() const {
    return BigInt::FromUint64(static_cast<uint64_t>(coset_reps_.size()));
  }

  std::vector<element_type> Generators() const { return generators_; }

  std::string ToString(const element_type& a) const {
    return detail::GroupElementToString(*group_, Canonicalize(a)) + "N";
  }

 private:
  static const G& RequireGroupRef(const std::shared_ptr<const G>& group) {
    if (!group) {
      throw std::invalid_argument(
          "QuotientGroup: group pointer cannot be null");
    }
    return *group;
  }

  std::shared_ptr<const G> group_;
  Subgroup<G> normal_;
  std::vector<element_type> coset_reps_;
  std::vector<element_type> generators_;
  std::unordered_map<element_type, element_type, detail::GroupElemHash<G>,
                     detail::GroupElemEqual<G>>
      elem_to_rep_;
};

// ============================================================================
// Center, CommutatorSubgroup, DerivedSeries, IsSolvable
// ============================================================================

template <Group G>
inline Subgroup<G> Center(std::shared_ptr<const G> group) {
  if (!group) {
    throw std::invalid_argument("Center: group pointer cannot be null");
  }
  using E = typename G::element_type;
  const auto all = detail::EnumerateGroup(*group);
  const std::vector<E> gens = group->Generators();
  std::vector<E> center_elems;
  for (const E& z : all.elements) {
    bool commutes = true;
    for (const E& s : gens) {
      if (!group->Equals(group->Op(z, s), group->Op(s, z))) {
        commutes = false;
        break;
      }
    }
    if (commutes) {
      center_elems.push_back(z);
    }
  }
  return Subgroup<G>::FromElements(std::move(group), std::move(center_elems));
}

template <Group G>
inline Subgroup<G> Center(std::shared_ptr<G> group) {
  return Center(std::shared_ptr<const G>(std::move(group)));
}

template <Group G>
inline Subgroup<G> Center(const G& group) {
  return Center(std::make_shared<const G>(group));
}

namespace detail {

// Computes the normal closure of `seed` inside the subgroup generated by
// `ambient_gens` within `*group`.
template <Group G>
inline Subgroup<G> NormalClosureInSubgroup(
    std::shared_ptr<const G> group,
    const std::vector<typename G::element_type>& ambient_gens,
    std::vector<typename G::element_type> seed) {
  using E = typename G::element_type;
  std::vector<E> cur_gens;
  for (const E& x : seed) {
    if (!group->IsIdentity(x)) {
      cur_gens.push_back(CanonicalizeElement(*group, x));
    }
  }
  Subgroup<G> k(group, cur_gens);
  if (k.IsTrivial()) {
    return k;
  }

  std::vector<E> conj_ops;
  conj_ops.reserve(ambient_gens.size() * 2);
  for (const E& s : ambient_gens) {
    if (!group->IsIdentity(s)) {
      conj_ops.push_back(s);
      conj_ops.push_back(group->Inv(s));
    }
  }

  while (true) {
    bool added = false;
    const std::size_t gen_count = cur_gens.size();
    for (std::size_t gi = 0; gi < gen_count; ++gi) {
      const E y = cur_gens[gi];
      for (const E& s : conj_ops) {
        const E c = CanonicalizeElement(
            *group, group->Op(group->Op(s, y), group->Inv(s)));
        if (!k.Contains(c)) {
          cur_gens.push_back(c);
          k = Subgroup<G>(group, cur_gens);
          added = true;
          break;
        }
      }
      if (added) {
        break;
      }
    }
    if (!added) {
      break;
    }
  }
  return k;
}

template <Group G>
inline Subgroup<G> CommutatorOfGenerators(
    std::shared_ptr<const G> group,
    const std::vector<typename G::element_type>& gens) {
  using E = typename G::element_type;
  std::vector<E> comm_seeds;
  for (std::size_t i = 0; i < gens.size(); ++i) {
    const E a_inv = group->Inv(gens[i]);
    for (std::size_t j = i + 1; j < gens.size(); ++j) {
      const E b_inv = group->Inv(gens[j]);
      // [a, b] = a * b * a^{-1} * b^{-1}
      const E comm =
          group->Op(group->Op(gens[i], gens[j]), group->Op(a_inv, b_inv));
      if (!group->IsIdentity(comm)) {
        comm_seeds.push_back(comm);
      }
    }
  }
  return NormalClosureInSubgroup(std::move(group), gens, std::move(comm_seeds));
}

}  // namespace detail

template <Group G>
inline Subgroup<G> CommutatorSubgroup(std::shared_ptr<const G> group) {
  if (!group) {
    throw std::invalid_argument(
        "CommutatorSubgroup: group pointer cannot be null");
  }
  const auto gens = group->Generators();
  return detail::CommutatorOfGenerators(std::move(group), gens);
}

template <Group G>
inline Subgroup<G> CommutatorSubgroup(std::shared_ptr<G> group) {
  return CommutatorSubgroup(std::shared_ptr<const G>(std::move(group)));
}

template <Group G>
inline Subgroup<G> CommutatorSubgroup(const G& group) {
  return CommutatorSubgroup(std::make_shared<const G>(group));
}

template <Group G>
inline Subgroup<G> CommutatorSubgroup(const Subgroup<G>& subgroup) {
  return detail::CommutatorOfGenerators(subgroup.AmbientGroup(),
                                        subgroup.Generators());
}

template <Group G>
inline std::vector<Subgroup<G>> DerivedSeries(std::shared_ptr<const G> group) {
  if (!group) {
    throw std::invalid_argument("DerivedSeries: group pointer cannot be null");
  }
  std::vector<Subgroup<G>> series;
  series.emplace_back(group, group->Generators());
  while (!series.back().IsTrivial()) {
    Subgroup<G> next = CommutatorSubgroup(series.back());
    if (next.Order() == series.back().Order()) {
      break;
    }
    series.push_back(std::move(next));
  }
  return series;
}

template <Group G>
inline std::vector<Subgroup<G>> DerivedSeries(std::shared_ptr<G> group) {
  return DerivedSeries(std::shared_ptr<const G>(std::move(group)));
}

template <Group G>
inline std::vector<Subgroup<G>> DerivedSeries(const G& group) {
  return DerivedSeries(std::make_shared<const G>(group));
}

template <Group G>
inline bool IsSolvable(const G& group) {
  const auto series = DerivedSeries(group);
  return !series.empty() && series.back().IsTrivial();
}

// ============================================================================
// SylowSubgroup: exact p-Sylow subgroup of order p^v || |G|
// ============================================================================

template <Group G>
inline Subgroup<G> SylowSubgroup(std::shared_ptr<const G> group, uint64_t p) {
  if (!group) {
    throw std::invalid_argument("SylowSubgroup: group pointer cannot be null");
  }
  if (!detail::IsPrimeU64(p)) {
    throw std::invalid_argument("SylowSubgroup: p must be a prime number");
  }
  using E = typename G::element_type;

  const auto closure = detail::EnumerateGroup(*group);
  const uint64_t n = static_cast<uint64_t>(closure.elements.size());

  uint64_t rem = n;
  uint64_t p_pow = 1;
  while (rem % p == 0) {
    rem /= p;
    p_pow *= p;
  }

  if (p_pow == 1) {
    return Subgroup<G>(std::move(group), {});
  }
  if (p_pow == n) {
    return Subgroup<G>::FromElements(group, closure.elements,
                                     group->Generators());
  }

  // Fast path when p^1 || |G|: locate an element of order divisible by p.
  if (p_pow == p) {
    for (std::size_t i = 1; i < closure.elements.size(); ++i) {
      const E& x = closure.elements[i];
      const uint64_t ord = ElementOrder(*group, x);
      if (ord % p == 0) {
        const E z = group->Pow(x, ord / p);
        return Subgroup<G>(std::move(group), {z});
      }
    }
  }

  // Build integer-indexed Cayley table and invoke exact constructive Sylow.
  detail::IndexedGroup ig;
  ig.n = closure.elements.size();
  ig.mul.assign(ig.n, std::vector<uint32_t>(ig.n, 0));
  ig.inv.resize(ig.n, 0);
  for (std::size_t i = 0; i < ig.n; ++i) {
    const E inv_i =
        detail::CanonicalizeElement(*group, group->Inv(closure.elements[i]));
    ig.inv[i] = static_cast<uint32_t>(closure.index_map.at(inv_i));
    for (std::size_t j = 0; j < ig.n; ++j) {
      const E prod = detail::CanonicalizeElement(
          *group, group->Op(closure.elements[i], closure.elements[j]));
      ig.mul[i][j] = static_cast<uint32_t>(closure.index_map.at(prod));
    }
  }

  const std::vector<uint32_t> sylow_indices =
      detail::FindSylowInIndexedGroup(ig, p, p_pow);

  std::vector<E> sylow_elems;
  sylow_elems.reserve(sylow_indices.size());
  for (uint32_t idx : sylow_indices) {
    sylow_elems.push_back(closure.elements[idx]);
  }
  return Subgroup<G>::FromElements(std::move(group), std::move(sylow_elems));
}

template <Group G>
inline Subgroup<G> SylowSubgroup(std::shared_ptr<G> group, uint64_t p) {
  return SylowSubgroup(std::shared_ptr<const G>(std::move(group)), p);
}

template <Group G>
inline Subgroup<G> SylowSubgroup(const G& group, uint64_t p) {
  return SylowSubgroup(std::make_shared<const G>(group), p);
}

// ============================================================================
// GroupHomomorphism<G, H> & VerifyFirstIsomorphismTheorem
// ============================================================================

template <Group G, Group H>
class GroupHomomorphism {
 public:
  using domain_type = G;
  using codomain_type = H;
  using domain_element_type = typename G::element_type;
  using codomain_element_type = typename H::element_type;
  using MapFn =
      std::function<codomain_element_type(const domain_element_type&)>;

  GroupHomomorphism(std::shared_ptr<const G> domain,
                    std::shared_ptr<const H> codomain, MapFn map)
      : domain_(std::move(domain)),
        codomain_(std::move(codomain)),
        map_(std::move(map)) {
    if (!domain_ || !codomain_) {
      throw std::invalid_argument(
          "GroupHomomorphism: domain and codomain pointers cannot be null");
    }
    if (!map_) {
      throw std::invalid_argument("GroupHomomorphism: map cannot be empty");
    }
  }

  // Constructs a homomorphism uniquely determined by specifying images of
  // generators of G. Throws std::invalid_argument if the assignment does not
  // extend to a well-defined homomorphism on G.
  static GroupHomomorphism FromGeneratorImages(
      std::shared_ptr<const G> domain, std::shared_ptr<const H> codomain,
      const std::vector<std::pair<domain_element_type, codomain_element_type>>&
          generator_images) {
    if (!domain || !codomain) {
      throw std::invalid_argument(
          "GroupHomomorphism::FromGeneratorImages: pointers cannot be null");
    }

    struct TableState {
      std::shared_ptr<const G> dom;
      std::shared_ptr<const H> cod;
      std::unordered_map<domain_element_type, codomain_element_type,
                         detail::GroupElemHash<G>, detail::GroupElemEqual<G>>
          table;

      TableState(std::shared_ptr<const G> d, std::shared_ptr<const H> c)
          : dom(std::move(d)),
            cod(std::move(c)),
            table(16, detail::GroupElemHash<G>{dom.get()},
                  detail::GroupElemEqual<G>{dom.get()}) {}
    };

    auto state = std::make_shared<TableState>(domain, codomain);

    std::vector<std::pair<domain_element_type, codomain_element_type>> steps;
    steps.reserve(generator_images.size() * 2);
    for (const auto& [g_raw, h_raw] : generator_images) {
      const auto g = detail::CanonicalizeElement(*domain, g_raw);
      const auto h = detail::CanonicalizeElement(*codomain, h_raw);
      steps.emplace_back(g, h);
      steps.emplace_back(
          detail::CanonicalizeElement(*domain, domain->Inv(g)),
          detail::CanonicalizeElement(*codomain, codomain->Inv(h)));
    }

    std::vector<domain_element_type> queue;
    const auto id_g = domain->Identity();
    const auto id_h = codomain->Identity();
    queue.push_back(id_g);
    state->table.emplace(id_g, id_h);

    for (std::size_t q_idx = 0; q_idx < queue.size(); ++q_idx) {
      const domain_element_type cur_g = queue[q_idx];
      const codomain_element_type cur_h = state->table.at(cur_g);
      for (const auto& [sg, sh] : steps) {
        const domain_element_type nxt_g =
            detail::CanonicalizeElement(*domain, domain->Op(cur_g, sg));
        const codomain_element_type exp_h =
            detail::CanonicalizeElement(*codomain, codomain->Op(cur_h, sh));
        const auto it = state->table.find(nxt_g);
        if (it == state->table.end()) {
          state->table.emplace(nxt_g, exp_h);
          queue.push_back(nxt_g);
        } else if (!codomain->Equals(it->second, exp_h)) {
          throw std::invalid_argument(
              "GroupHomomorphism::FromGeneratorImages: generator assignment "
              "does not define a valid homomorphism");
        }
      }
    }

    // Verify the provided generators actually generate all of domain.
    const auto full = detail::EnumerateGroup(*domain);
    if (state->table.size() != full.elements.size()) {
      throw std::invalid_argument(
          "GroupHomomorphism::FromGeneratorImages: provided elements do not "
          "generate the full domain group");
    }

    MapFn lookup = [state](const domain_element_type& x) {
      const auto canon = detail::CanonicalizeElement(*state->dom, x);
      const auto it = state->table.find(canon);
      if (it == state->table.end()) {
        throw std::invalid_argument(
            "GroupHomomorphism: element not in domain group");
      }
      return it->second;
    };
    return GroupHomomorphism(std::move(domain), std::move(codomain),
                             std::move(lookup));
  }

  static GroupHomomorphism FromGeneratorImages(
      std::shared_ptr<const G> domain, std::shared_ptr<const H> codomain,
      const std::vector<codomain_element_type>& images) {
    if (!domain) {
      throw std::invalid_argument(
          "GroupHomomorphism::FromGeneratorImages: domain cannot be null");
    }
    const auto dom_gens = domain->Generators();
    if (dom_gens.size() != images.size()) {
      throw std::invalid_argument(
          "GroupHomomorphism::FromGeneratorImages: image count must match "
          "domain generator count");
    }
    std::vector<std::pair<domain_element_type, codomain_element_type>> pairs;
    pairs.reserve(dom_gens.size());
    for (std::size_t i = 0; i < dom_gens.size(); ++i) {
      pairs.emplace_back(dom_gens[i], images[i]);
    }
    return FromGeneratorImages(std::move(domain), std::move(codomain), pairs);
  }

  std::shared_ptr<const G> Domain() const noexcept { return domain_; }

  std::shared_ptr<const H> Codomain() const noexcept { return codomain_; }

  codomain_element_type Apply(const domain_element_type& x) const {
    return detail::CanonicalizeElement(*codomain_, map_(x));
  }

  codomain_element_type operator()(const domain_element_type& x) const {
    return Apply(x);
  }

  // Verifies phi(e_G) == e_H and phi(a * s) == phi(a) * phi(s) for all a in G
  // and all generators s of G (equivalent to phi(ab) == phi(a)phi(b) on all of
  // G by induction on generator word length).
  bool VerifyHomomorphism() const {
    if (!codomain_->IsIdentity(Apply(domain_->Identity()))) {
      return false;
    }
    const auto all = detail::EnumerateGroup(*domain_);
    const auto gens = domain_->Generators();
    std::vector<codomain_element_type> gen_imgs;
    gen_imgs.reserve(gens.size());
    for (const auto& s : gens) {
      gen_imgs.push_back(Apply(s));
    }
    for (const auto& a : all.elements) {
      const auto phi_a = Apply(a);
      for (std::size_t gi = 0; gi < gens.size(); ++gi) {
        const auto lhs = Apply(domain_->Op(a, gens[gi]));
        const auto rhs = codomain_->Op(phi_a, gen_imgs[gi]);
        if (!codomain_->Equals(lhs, rhs)) {
          return false;
        }
      }
    }
    return true;
  }

  Subgroup<G> Kernel() const {
    const auto all = detail::EnumerateGroup(*domain_);
    std::vector<domain_element_type> ker_elems;
    for (const auto& x : all.elements) {
      if (codomain_->IsIdentity(Apply(x))) {
        ker_elems.push_back(x);
      }
    }
    return Subgroup<G>::FromElements(domain_, std::move(ker_elems));
  }

  Subgroup<H> Image() const {
    const auto gens = domain_->Generators();
    std::vector<codomain_element_type> img_gens;
    img_gens.reserve(gens.size());
    for (const auto& s : gens) {
      img_gens.push_back(Apply(s));
    }
    return Subgroup<H>(codomain_, std::move(img_gens));
  }

  bool IsInjective() const {
    return VerifyHomomorphism() && Kernel().Order() == 1;
  }

  bool IsSurjective() const {
    if (!VerifyHomomorphism()) {
      return false;
    }
    const uint64_t img_ord = Image().Order();
    if constexpr (FiniteGroup<H>) {
      return BigInt::FromUint64(img_ord) == codomain_->Order();
    } else {
      const auto cod_all = detail::EnumerateGroup(*codomain_);
      return img_ord == cod_all.elements.size();
    }
  }

  bool IsIsomorphism() const { return IsInjective() && IsSurjective(); }

 private:
  std::shared_ptr<const G> domain_;
  std::shared_ptr<const H> codomain_;
  MapFn map_;
};

// Verifies the First Isomorphism Theorem for phi: G -> H:
//   1. phi is a valid homomorphism,
//   2. ker(phi) is normal in G,
//   3. |G / ker(phi)| == |im(phi)|,
//   4. The induced map phi_bar: G / ker(phi) -> H is well-defined on cosets,
//      is an injective homomorphism, and has image equal to im(phi).
template <Group G, Group H>
inline bool VerifyFirstIsomorphismTheorem(const GroupHomomorphism<G, H>& phi) {
  if (!phi.VerifyHomomorphism()) {
    return false;
  }
  const Subgroup<G> ker = phi.Kernel();
  const Subgroup<H> img = phi.Image();
  if (!IsNormalSubgroup(*phi.Domain(), ker)) {
    return false;
  }

  auto quot = std::make_shared<QuotientGroup<G>>(phi.Domain(), ker);
  if (quot->Order() != BigInt::FromUint64(img.Order())) {
    return false;
  }

  // Verify well-definedness on every coset of ker(phi).
  const auto cosets = LeftCosets(*phi.Domain(), ker);
  for (const auto& coset : cosets) {
    const auto rep_img = phi.Apply(quot->Canonicalize(coset.front()));
    for (const auto& x : coset) {
      if (!phi.Codomain()->Equals(phi.Apply(x), rep_img)) {
        return false;
      }
    }
  }

  // Induced homomorphism phi_bar: G / ker(phi) -> H.
  GroupHomomorphism<QuotientGroup<G>, H> induced(
      quot, phi.Codomain(),
      [&phi, quot](const typename G::element_type& coset_rep) {
        return phi.Apply(quot->Canonicalize(coset_rep));
      });

  if (!induced.VerifyHomomorphism()) {
    return false;
  }
  if (!induced.IsInjective()) {
    return false;
  }
  return induced.Image() == img;
}

// ============================================================================
// CayleyTable
// ============================================================================

template <Group G>
struct CayleyTableResult {
  std::vector<typename G::element_type> elements;
  std::vector<std::vector<std::size_t>> table;
};

template <Group G>
inline CayleyTableResult<G> CayleyTable(const G& group) {
  const auto closure = detail::EnumerateGroup(group);
  const std::size_t n = closure.elements.size();
  CayleyTableResult<G> res;
  res.elements = closure.elements;
  res.table.assign(n, std::vector<std::size_t>(n, 0));

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const auto prod = detail::CanonicalizeElement(
          group, group.Op(res.elements[i], res.elements[j]));
      res.table[i][j] = closure.index_map.at(prod);
    }
  }
  return res;
}

// ============================================================================
// BigInt Smith Normal Form Abelian Group Decomposition
// ============================================================================

inline AbelianGroupDecompositionResult<BigInt> DecomposeAbelianGroupBigInt(
    const std::vector<std::vector<BigInt>>& relations) {
  return DecomposeFinitelyGeneratedAbelianGroup<BigInt>(relations);
}

}  // namespace algebra

template <Group G>
struct AlgebraTraits<algebra::QuotientGroup<G>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = detail::TraitIsCommutativeOp<G>();
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

using algebra::CayleyTable;
using algebra::CayleyTableResult;
using algebra::Center;
using algebra::CommutatorSubgroup;
using algebra::DecomposeAbelianGroupBigInt;
using algebra::DerivedSeries;
using algebra::ElementOrder;
using algebra::GroupHomomorphism;
using algebra::IsNormalSubgroup;
using algebra::IsSolvable;
using algebra::LeftCosets;
using algebra::QuotientGroup;
using algebra::RightCosets;
using algebra::Subgroup;
using algebra::SylowSubgroup;
using algebra::VerifyFirstIsomorphismTheorem;

}  // namespace xg

#endif  // XGALOIS_ALGEBRA_GROUP_ALGORITHMS_HPP_
