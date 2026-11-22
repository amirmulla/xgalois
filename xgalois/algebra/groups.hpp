#ifndef XGALOIS_ALGEBRA_GROUPS_HPP_
#define XGALOIS_ALGEBRA_GROUPS_HPP_

#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/group_base.hpp"
#include "xgalois/algebra/permutation.hpp"
#include "xgalois/algebra/schreier_sims.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace algebra {

namespace detail {

__extension__ using group_u128 = unsigned __int128;

inline uint64_t AddModU64(uint64_t a, uint64_t b, uint64_t mod) noexcept {
  return static_cast<uint64_t>(
      (static_cast<group_u128>(a % mod) + static_cast<group_u128>(b % mod)) %
      mod);
}

inline uint64_t SubModU64(uint64_t a, uint64_t b, uint64_t mod) noexcept {
  const uint64_t ar = a % mod;
  const uint64_t br = b % mod;
  return (ar >= br) ? (ar - br) : (mod - (br - ar));
}

inline uint64_t MulModU64(uint64_t a, uint64_t b, uint64_t mod) noexcept {
  return static_cast<uint64_t>(
      (static_cast<group_u128>(a % mod) * static_cast<group_u128>(b % mod)) %
      mod);
}

}  // namespace detail

// ============================================================================
// CyclicGroup: (Z / nZ, +)
// ============================================================================

class CyclicGroup {
 public:
  using element_type = uint64_t;

  explicit CyclicGroup(uint64_t order) : order_(order) {
    if (order_ == 0) {
      throw std::invalid_argument("CyclicGroup: order must be >= 1");
    }
  }

  uint64_t Modulus() const noexcept { return order_; }

  element_type Op(const element_type& a, const element_type& b) const noexcept {
    return detail::AddModU64(a, b, order_);
  }

  bool Equals(const element_type& a, const element_type& b) const noexcept {
    return (a % order_) == (b % order_);
  }

  element_type Identity() const noexcept { return 0; }

  bool IsIdentity(const element_type& a) const noexcept {
    return (a % order_) == 0;
  }

  element_type Inv(const element_type& a) const noexcept {
    const uint64_t r = a % order_;
    return (r == 0) ? 0 : (order_ - r);
  }

  element_type Div(const element_type& a,
                   const element_type& b) const noexcept {
    return detail::SubModU64(a, b, order_);
  }

  element_type Pow(const element_type& a, uint64_t exp) const noexcept {
    return detail::MulModU64(a, exp, order_);
  }

  uint64_t ElementOrder(const element_type& a) const noexcept {
    const uint64_t r = a % order_;
    if (r == 0) {
      return 1;
    }
    return order_ / std::gcd(r, order_);
  }

  BigInt Order() const { return BigInt::FromUint64(order_); }

  std::vector<element_type> Generators() const {
    if (order_ == 1) {
      return {0};
    }
    return {1};
  }

  std::string ToString(const element_type& a) const {
    return std::to_string(a % order_);
  }

 private:
  uint64_t order_ = 1;
};

// ============================================================================
// SymmetricGroup: S_n
// ============================================================================

class SymmetricGroup {
 public:
  using element_type = Permutation;

  explicit SymmetricGroup(uint32_t degree) : degree_(degree) {
    if (degree_ == 0) {
      throw std::invalid_argument("SymmetricGroup: degree must be >= 1");
    }
    if (degree_ == 1) {
      generators_ = {Permutation::Identity(1)};
    } else if (degree_ == 2) {
      generators_ = {Permutation::FromCycles(2, {{0, 1}})};
    } else {
      std::vector<uint32_t> n_cycle(degree_);
      std::iota(n_cycle.begin(), n_cycle.end(), uint32_t{0});
      generators_ = {
          Permutation::FromCycles(degree_, {{0, 1}}),
          Permutation::FromCycles(degree_, {n_cycle}),
      };
    }
    bsgs_ = BuildBsgsDeterministic(degree_, generators_);
  }

  uint32_t Degree() const noexcept { return degree_; }

  const BsgsChain& Bsgs() const noexcept { return bsgs_; }

  bool Contains(const element_type& p) const noexcept {
    return p.Degree() == degree_;
  }

  element_type Op(const element_type& a, const element_type& b) const {
    RequireMember(a);
    RequireMember(b);
    return a.Compose(b);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    RequireMember(a);
    RequireMember(b);
    return a == b;
  }

  element_type Identity() const { return Permutation::Identity(degree_); }

  bool IsIdentity(const element_type& a) const {
    RequireMember(a);
    return a.IsIdentity();
  }

  element_type Inv(const element_type& a) const {
    RequireMember(a);
    return a.Inverse();
  }

  element_type Div(const element_type& a, const element_type& b) const {
    RequireMember(a);
    RequireMember(b);
    return a.Compose(b.Inverse());
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    RequireMember(a);
    return a.PowUnsigned(exp);
  }

  uint64_t ElementOrder(const element_type& a) const {
    RequireMember(a);
    return a.Order().ToUint64();
  }

  BigInt Order() const { return BigInt::Factorial(degree_); }

  std::vector<element_type> Generators() const { return generators_; }

  std::string ToString(const element_type& a) const {
    RequireMember(a);
    return a.ToString();
  }

 private:
  void RequireMember(const element_type& p) const {
    if (p.Degree() != degree_) {
      throw std::invalid_argument(
          "SymmetricGroup: permutation degree mismatch");
    }
  }

  uint32_t degree_ = 1;
  std::vector<element_type> generators_;
  BsgsChain bsgs_;
};

// ============================================================================
// AlternatingGroup: A_n
// ============================================================================

class AlternatingGroup {
 public:
  using element_type = Permutation;

  explicit AlternatingGroup(uint32_t degree) : degree_(degree) {
    if (degree_ == 0) {
      throw std::invalid_argument("AlternatingGroup: degree must be >= 1");
    }
    if (degree_ <= 2) {
      generators_ = {Permutation::Identity(degree_)};
    } else {
      generators_.reserve(degree_ - 2);
      for (uint32_t k = 2; k < degree_; ++k) {
        generators_.push_back(Permutation::FromCycles(degree_, {{0, 1, k}}));
      }
    }
    bsgs_ = BuildBsgsDeterministic(degree_, generators_);
  }

  uint32_t Degree() const noexcept { return degree_; }

  const BsgsChain& Bsgs() const noexcept { return bsgs_; }

  bool Contains(const element_type& p) const {
    return p.Degree() == degree_ && p.IsEven();
  }

  element_type Op(const element_type& a, const element_type& b) const {
    RequireMember(a);
    RequireMember(b);
    return a.Compose(b);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    RequireMember(a);
    RequireMember(b);
    return a == b;
  }

  element_type Identity() const { return Permutation::Identity(degree_); }

  bool IsIdentity(const element_type& a) const {
    RequireMember(a);
    return a.IsIdentity();
  }

  element_type Inv(const element_type& a) const {
    RequireMember(a);
    return a.Inverse();
  }

  element_type Div(const element_type& a, const element_type& b) const {
    RequireMember(a);
    RequireMember(b);
    return a.Compose(b.Inverse());
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    RequireMember(a);
    return a.PowUnsigned(exp);
  }

  uint64_t ElementOrder(const element_type& a) const {
    RequireMember(a);
    return a.Order().ToUint64();
  }

  BigInt Order() const {
    if (degree_ <= 2) {
      return BigInt(1);
    }
    return BigInt::Factorial(degree_) / BigInt(2);
  }

  std::vector<element_type> Generators() const { return generators_; }

  std::string ToString(const element_type& a) const {
    RequireMember(a);
    return a.ToString();
  }

 private:
  void RequireMember(const element_type& p) const {
    if (p.Degree() != degree_) {
      throw std::invalid_argument(
          "AlternatingGroup: permutation degree mismatch");
    }
    if (!p.IsEven()) {
      throw std::domain_error(
          "AlternatingGroup: odd permutation is not in A_n");
    }
  }

  uint32_t degree_ = 1;
  std::vector<element_type> generators_;
  BsgsChain bsgs_;
};

// ============================================================================
// DihedralGroup: D_n of order 2n
// ============================================================================

struct DihedralElement {
  uint64_t rot = 0;
  bool flip = false;

  bool operator==(const DihedralElement&) const noexcept = default;
  std::strong_ordering operator<=>(const DihedralElement&) const noexcept =
      default;

  std::string ToString() const {
    if (!flip) {
      if (rot == 0) {
        return "e";
      }
      return (rot == 1) ? "r" : ("r^" + std::to_string(rot));
    }
    if (rot == 0) {
      return "s";
    }
    return (rot == 1) ? "r*s" : ("r^" + std::to_string(rot) + "*s");
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const DihedralElement& elem) {
    return os << elem.ToString();
  }
};

// Dihedral group D_n of order 2n with elements r^k * s^b satisfying
// r^n = e, s^2 = e, and s * r^k = r^{-k} * s.
class DihedralGroup {
 public:
  using element_type = DihedralElement;

  explicit DihedralGroup(uint64_t n) : n_(n) {
    if (n_ == 0) {
      throw std::invalid_argument("DihedralGroup: n must be >= 1");
    }
  }

  uint64_t N() const noexcept { return n_; }

  element_type Rotation(uint64_t k = 1) const noexcept {
    return DihedralElement{k % n_, false};
  }

  element_type Reflection(uint64_t k = 0) const noexcept {
    return DihedralElement{k % n_, true};
  }

  element_type Op(const element_type& a, const element_type& b) const noexcept {
    if (!a.flip) {
      return DihedralElement{detail::AddModU64(a.rot, b.rot, n_), b.flip};
    }
    return DihedralElement{detail::SubModU64(a.rot, b.rot, n_), !b.flip};
  }

  bool Equals(const element_type& a, const element_type& b) const noexcept {
    return (a.rot % n_) == (b.rot % n_) && a.flip == b.flip;
  }

  element_type Identity() const noexcept { return DihedralElement{0, false}; }

  bool IsIdentity(const element_type& a) const noexcept {
    return (a.rot % n_) == 0 && !a.flip;
  }

  element_type Inv(const element_type& a) const noexcept {
    const uint64_t r = a.rot % n_;
    if (a.flip) {
      return DihedralElement{r, true};
    }
    return DihedralElement{(r == 0) ? 0 : (n_ - r), false};
  }

  element_type Div(const element_type& a,
                   const element_type& b) const noexcept {
    return Op(a, Inv(b));
  }

  element_type Pow(const element_type& a, uint64_t exp) const noexcept {
    if (a.flip) {
      return ((exp & uint64_t{1}) != 0) ? DihedralElement{a.rot % n_, true}
                                        : Identity();
    }
    return DihedralElement{detail::MulModU64(a.rot, exp, n_), false};
  }

  uint64_t ElementOrder(const element_type& a) const noexcept {
    if (a.flip) {
      return 2;
    }
    const uint64_t r = a.rot % n_;
    if (r == 0) {
      return 1;
    }
    return n_ / std::gcd(r, n_);
  }

  BigInt Order() const { return BigInt::FromUint64(n_) * BigInt(2); }

  std::vector<element_type> Generators() const {
    if (n_ == 1) {
      return {DihedralElement{0, true}};
    }
    return {DihedralElement{1, false}, DihedralElement{0, true}};
  }

  std::string ToString(const element_type& a) const {
    return DihedralElement{a.rot % n_, a.flip}.ToString();
  }

 private:
  uint64_t n_ = 1;
};

// ============================================================================
// PermutationGroup: arbitrary subgroup of S_n backed by BsgsChain
// ============================================================================

class PermutationGroup {
 public:
  using element_type = Permutation;

  PermutationGroup(uint32_t degree, std::vector<Permutation> generators)
      : degree_(degree),
        generators_(std::move(generators)),
        bsgs_(BuildBsgsDeterministic(degree_, generators_)) {
    if (degree_ == 0) {
      throw std::invalid_argument("PermutationGroup: degree must be >= 1");
    }
    if (generators_.empty()) {
      generators_.push_back(Permutation::Identity(degree_));
    }
  }

  static PermutationGroup FromBsgs(BsgsChain bsgs) {
    if (bsgs.Degree() == 0) {
      throw std::invalid_argument(
          "PermutationGroup::FromBsgs: degree must be >= 1");
    }
    PermutationGroup group;
    group.degree_ = bsgs.Degree();
    group.generators_ = bsgs.StrongGenerators();
    if (group.generators_.empty()) {
      group.generators_.push_back(Permutation::Identity(group.degree_));
    }
    group.bsgs_ = std::move(bsgs);
    return group;
  }

  uint32_t Degree() const noexcept { return degree_; }

  const BsgsChain& Bsgs() const noexcept { return bsgs_; }

  bool Contains(const element_type& p) const { return bsgs_.Contains(p); }

  element_type Op(const element_type& a, const element_type& b) const {
    RequireDegree(a);
    RequireDegree(b);
    return a.Compose(b);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    RequireDegree(a);
    RequireDegree(b);
    return a == b;
  }

  element_type Identity() const { return Permutation::Identity(degree_); }

  bool IsIdentity(const element_type& a) const {
    RequireDegree(a);
    return a.IsIdentity();
  }

  element_type Inv(const element_type& a) const {
    RequireDegree(a);
    return a.Inverse();
  }

  element_type Div(const element_type& a, const element_type& b) const {
    RequireDegree(a);
    RequireDegree(b);
    return a.Compose(b.Inverse());
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    RequireDegree(a);
    return a.PowUnsigned(exp);
  }

  uint64_t ElementOrder(const element_type& a) const {
    RequireDegree(a);
    return a.Order().ToUint64();
  }

  BigInt Order() const { return bsgs_.Order(); }

  std::vector<element_type> Generators() const { return generators_; }

  std::string ToString(const element_type& a) const {
    RequireDegree(a);
    return a.ToString();
  }

 private:
  PermutationGroup() = default;

  void RequireDegree(const element_type& p) const {
    if (p.Degree() != degree_) {
      throw std::invalid_argument(
          "PermutationGroup: permutation degree mismatch");
    }
  }

  uint32_t degree_ = 1;
  std::vector<element_type> generators_;
  BsgsChain bsgs_;
};

// ============================================================================
// DirectProductGroup<G1, G2>: G1 x G2
// ============================================================================

template <Group G1, Group G2>
class DirectProductGroup {
 public:
  using first_group_type = G1;
  using second_group_type = G2;
  using element_type =
      std::pair<typename G1::element_type, typename G2::element_type>;

  DirectProductGroup(G1 g1, G2 g2) : g1_(std::move(g1)), g2_(std::move(g2)) {}

  const G1& First() const noexcept { return g1_; }
  const G2& Second() const noexcept { return g2_; }

  element_type Op(const element_type& a, const element_type& b) const {
    return element_type{g1_.Op(a.first, b.first), g2_.Op(a.second, b.second)};
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return static_cast<bool>(g1_.Equals(a.first, b.first)) &&
           static_cast<bool>(g2_.Equals(a.second, b.second));
  }

  element_type Identity() const {
    return element_type{g1_.Identity(), g2_.Identity()};
  }

  bool IsIdentity(const element_type& a) const {
    return static_cast<bool>(g1_.IsIdentity(a.first)) &&
           static_cast<bool>(g2_.IsIdentity(a.second));
  }

  element_type Inv(const element_type& a) const {
    return element_type{g1_.Inv(a.first), g2_.Inv(a.second)};
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return element_type{g1_.Div(a.first, b.first), g2_.Div(a.second, b.second)};
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    return element_type{g1_.Pow(a.first, exp), g2_.Pow(a.second, exp)};
  }

  BigInt Order() const
    requires(FiniteGroup<G1> && FiniteGroup<G2>)
  {
    return g1_.Order() * g2_.Order();
  }

  std::vector<element_type> Generators() const
    requires(FiniteGroup<G1> && FiniteGroup<G2>)
  {
    std::vector<element_type> gens;
    const auto id1 = g1_.Identity();
    const auto id2 = g2_.Identity();
    for (const auto& x : g1_.Generators()) {
      if (!g1_.IsIdentity(x)) {
        gens.emplace_back(x, id2);
      }
    }
    for (const auto& y : g2_.Generators()) {
      if (!g2_.IsIdentity(y)) {
        gens.emplace_back(id1, y);
      }
    }
    if (gens.empty()) {
      gens.emplace_back(id1, id2);
    }
    return gens;
  }

  std::string ToString(const element_type& a) const {
    return "(" + detail::GroupElementToString(g1_, a.first) + ", " +
           detail::GroupElementToString(g2_, a.second) + ")";
  }

 private:
  G1 g1_;
  G2 g2_;
};

// ============================================================================
// SemidirectProductGroup<N, H, ActionFn>: N |x|_phi H
//   (n1, h1) * (n2, h2) = (n1 * phi(h1)(n2), h1 * h2)
// ============================================================================

template <Group N, Group H, typename ActionFn>
class SemidirectProductGroup {
 public:
  using normal_group_type = N;
  using acting_group_type = H;
  using action_type = ActionFn;
  using element_type =
      std::pair<typename N::element_type, typename H::element_type>;

  SemidirectProductGroup(N normal, H acting, ActionFn action)
      : normal_(std::move(normal)),
        acting_(std::move(acting)),
        action_(std::move(action)) {}

  const N& NormalGroup() const noexcept { return normal_; }
  const H& ActingGroup() const noexcept { return acting_; }

  typename N::element_type ApplyAction(
      const typename H::element_type& h,
      const typename N::element_type& n) const {
    if constexpr (std::invocable<const ActionFn&,
                                 const typename H::element_type&,
                                 const typename N::element_type&>) {
      return action_(h, n);
    } else {
      return action_(h)(n);
    }
  }

  element_type Op(const element_type& a, const element_type& b) const {
    return element_type{normal_.Op(a.first, ApplyAction(a.second, b.first)),
                        acting_.Op(a.second, b.second)};
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return static_cast<bool>(normal_.Equals(a.first, b.first)) &&
           static_cast<bool>(acting_.Equals(a.second, b.second));
  }

  element_type Identity() const {
    return element_type{normal_.Identity(), acting_.Identity()};
  }

  bool IsIdentity(const element_type& a) const {
    return static_cast<bool>(normal_.IsIdentity(a.first)) &&
           static_cast<bool>(acting_.IsIdentity(a.second));
  }

  element_type Inv(const element_type& a) const {
    const auto h_inv = acting_.Inv(a.second);
    return element_type{ApplyAction(h_inv, normal_.Inv(a.first)), h_inv};
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return Op(a, Inv(b));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type acc = Identity();
    if (exp == 0) {
      return acc;
    }
    element_type base = a;
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Op(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Op(base, base);
      }
    }
    return acc;
  }

  BigInt Order() const
    requires(FiniteGroup<N> && FiniteGroup<H>)
  {
    return normal_.Order() * acting_.Order();
  }

  std::vector<element_type> Generators() const
    requires(FiniteGroup<N> && FiniteGroup<H>)
  {
    std::vector<element_type> gens;
    const auto id_n = normal_.Identity();
    const auto id_h = acting_.Identity();
    for (const auto& n : normal_.Generators()) {
      if (!normal_.IsIdentity(n)) {
        gens.emplace_back(n, id_h);
      }
    }
    for (const auto& h : acting_.Generators()) {
      if (!acting_.IsIdentity(h)) {
        gens.emplace_back(id_n, h);
      }
    }
    if (gens.empty()) {
      gens.emplace_back(id_n, id_h);
    }
    return gens;
  }

  std::string ToString(const element_type& a) const {
    return "(" + detail::GroupElementToString(normal_, a.first) + ", " +
           detail::GroupElementToString(acting_, a.second) + ")";
  }

 private:
  N normal_;
  H acting_;
  ActionFn action_;
};

}  // namespace algebra

// ============================================================================
// AlgebraTraits specializations in namespace xg
// ============================================================================

template <>
struct AlgebraTraits<algebra::CyclicGroup> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = true;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <>
struct AlgebraTraits<algebra::SymmetricGroup> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <>
struct AlgebraTraits<algebra::AlternatingGroup> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <>
struct AlgebraTraits<algebra::DihedralGroup> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <>
struct AlgebraTraits<algebra::PermutationGroup> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <Group G1, Group G2>
struct AlgebraTraits<algebra::DirectProductGroup<G1, G2>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp =
      detail::TraitIsCommutativeOp<G1>() && detail::TraitIsCommutativeOp<G2>();
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <Group N, Group H, typename ActionFn>
struct AlgebraTraits<algebra::SemidirectProductGroup<N, H, ActionFn>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

using algebra::AlternatingGroup;
using algebra::CyclicGroup;
using algebra::DihedralElement;
using algebra::DihedralGroup;
using algebra::DirectProductGroup;
using algebra::PermutationGroup;
using algebra::SemidirectProductGroup;
using algebra::SymmetricGroup;

}  // namespace xg

namespace std {

template <>
struct hash<xg::algebra::DihedralElement> {
  size_t operator()(const xg::algebra::DihedralElement& elem) const noexcept {
    const size_t h1 = std::hash<uint64_t>{}(elem.rot);
    const size_t h2 = std::hash<bool>{}(elem.flip);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

}  // namespace std

#endif  // XGALOIS_ALGEBRA_GROUPS_HPP_
