#ifndef XGALOIS_RING_QUADRATIC_INTEGERS_HPP_
#define XGALOIS_RING_QUADRATIC_INTEGERS_HPP_

#include <array>
#include <cstdint>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace ring {

namespace detail {

// Exact nearest-integer quotient q minimizing |num - q * den| for den > 0.
// Guarantees |num - q * den| <= den / 2.
inline BigInt RoundDiv(const BigInt& num, const BigInt& den) {
  auto [q, r] = BigInt::EuclideanDivRem(num, den);
  if (r + r > den) {
    q += BigInt(1);
  }
  return q;
}

}  // namespace detail

// ============================================================================
// Gaussian Integers: Z[i] = {re + im * i : re, im in Z}
// ============================================================================

struct GaussianInt {
  BigInt re{0};
  BigInt im{0};

  bool operator==(const GaussianInt& other) const noexcept = default;
};

inline std::ostream& operator<<(std::ostream& os, const GaussianInt& z) {
  os << z.re.ToString();
  if (z.im.IsNegative()) {
    os << " - " << (-z.im).ToString() << "i";
  } else {
    os << " + " << z.im.ToString() << "i";
  }
  return os;
}

class GaussianIntegerRing {
 public:
  using element_type = GaussianInt;

  GaussianIntegerRing() = default;

  element_type Add(const element_type& x, const element_type& y) const {
    return {x.re + y.re, x.im + y.im};
  }

  element_type Sub(const element_type& x, const element_type& y) const {
    return {x.re - y.re, x.im - y.im};
  }

  element_type Neg(const element_type& x) const { return {-x.re, -x.im}; }

  element_type Mul(const element_type& x, const element_type& y) const {
    return {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re};
  }

  element_type Pow(const element_type& x, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = x;
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  element_type Conjugate(const element_type& x) const { return {x.re, -x.im}; }

  BigInt Norm(const element_type& x) const { return x.re * x.re + x.im * x.im; }

  BigInt EuclideanNorm(const element_type& x) const { return Norm(x); }

  element_type AdditiveIdentity() const { return {BigInt(0), BigInt(0)}; }

  element_type MultiplicativeIdentity() const { return {BigInt(1), BigInt(0)}; }

  bool IsZero(const element_type& x) const noexcept {
    return x.re.IsZero() && x.im.IsZero();
  }

  bool IsOne(const element_type& x) const noexcept {
    return x.re.IsOne() && x.im.IsZero();
  }

  bool Equals(const element_type& x, const element_type& y) const noexcept {
    return x == y;
  }

  std::string ToString(const element_type& x) const {
    if (x.im.IsZero()) {
      return x.re.ToString();
    }
    if (x.re.IsZero()) {
      return x.im.ToString() + "i";
    }
    if (x.im.IsNegative()) {
      return x.re.ToString() + " - " + (-x.im).ToString() + "i";
    }
    return x.re.ToString() + " + " + x.im.ToString() + "i";
  }

  bool IsUnit(const element_type& x) const { return Norm(x).IsOne(); }

  std::optional<element_type> TryInv(const element_type& x) const {
    if (!IsUnit(x)) {
      return std::nullopt;
    }
    // For N(x) == 1, x * conj(x) == 1.
    return Conjugate(x);
  }

  // Canonical associate in the first quadrant (re > 0, im >= 0).
  // Returns {canonical, unit} such that canonical == unit * x.
  std::pair<element_type, element_type> AssociateCanonical(
      const element_type& x) const {
    if (IsZero(x)) {
      return {AdditiveIdentity(), MultiplicativeIdentity()};
    }
    if (x.re.IsPositive() && !x.im.IsNegative()) {
      return {x, {BigInt(1), BigInt(0)}};
    }
    if (!x.re.IsPositive() && x.im.IsPositive()) {
      // Multiply by -i: (a + bi)(-i) = b - ai.
      const element_type unit{BigInt(0), BigInt(-1)};
      return {{x.im, -x.re}, unit};
    }
    if (x.re.IsNegative() && !x.im.IsPositive()) {
      // Multiply by -1: -(a + bi) = -a - bi.
      const element_type unit{BigInt(-1), BigInt(0)};
      return {{-x.re, -x.im}, unit};
    }
    // x.re >= 0 && x.im < 0: multiply by i: (a + bi)i = -b + ai.
    const element_type unit{BigInt(0), BigInt(1)};
    return {{-x.im, x.re}, unit};
  }

  // Euclidean division in Z[i] via nearest-integer complex rounding.
  // Guarantees x == q * y + r and N(r) <= N(y) / 2 < N(y).
  std::pair<element_type, element_type> DivRem(const element_type& x,
                                               const element_type& y) const {
    if (IsZero(y)) {
      throw std::domain_error("GaussianIntegerRing::DivRem: division by zero");
    }
    const BigInt norm_y = Norm(y);
    const BigInt num_re = x.re * y.re + x.im * y.im;
    const BigInt num_im = x.im * y.re - x.re * y.im;
    element_type q{detail::RoundDiv(num_re, norm_y),
                   detail::RoundDiv(num_im, norm_y)};
    element_type r = Sub(x, Mul(q, y));
    return {std::move(q), std::move(r)};
  }

  // Extended Euclidean Algorithm: returns {g, s, t} with g canonical in the
  // first quadrant (or 0 if a == b == 0) and s * a + t * b == g.
  std::tuple<element_type, element_type, element_type> ExtendedGcd(
      const element_type& a, const element_type& b) const {
    element_type old_r = a;
    element_type r = b;
    element_type old_s = MultiplicativeIdentity();
    element_type s = AdditiveIdentity();
    element_type old_t = AdditiveIdentity();
    element_type t = MultiplicativeIdentity();

    while (!IsZero(r)) {
      auto [q, rem] = DivRem(old_r, r);
      old_r = std::move(r);
      r = std::move(rem);

      element_type next_s = Sub(old_s, Mul(q, s));
      old_s = std::move(s);
      s = std::move(next_s);

      element_type next_t = Sub(old_t, Mul(q, t));
      old_t = std::move(t);
      t = std::move(next_t);
    }

    auto [g_canon, unit] = AssociateCanonical(old_r);
    return {std::move(g_canon), Mul(unit, old_s), Mul(unit, old_t)};
  }
};

// ============================================================================
// Eisenstein Integers: Z[omega] = {a + b * omega : a, b in Z},
// omega^2+omega+1=0
// ============================================================================

struct EisensteinInt {
  BigInt a{0};
  BigInt b{0};

  bool operator==(const EisensteinInt& other) const noexcept = default;
};

inline std::ostream& operator<<(std::ostream& os, const EisensteinInt& z) {
  os << z.a.ToString();
  if (z.b.IsNegative()) {
    os << " - " << (-z.b).ToString() << "w";
  } else {
    os << " + " << z.b.ToString() << "w";
  }
  return os;
}

class EisensteinIntegerRing {
 public:
  using element_type = EisensteinInt;

  EisensteinIntegerRing() = default;

  element_type Add(const element_type& x, const element_type& y) const {
    return {x.a + y.a, x.b + y.b};
  }

  element_type Sub(const element_type& x, const element_type& y) const {
    return {x.a - y.a, x.b - y.b};
  }

  element_type Neg(const element_type& x) const { return {-x.a, -x.b}; }

  // (a + b*w)(c + d*w) = (a*c - b*d) + (a*d + b*c - b*d)*w
  element_type Mul(const element_type& x, const element_type& y) const {
    const BigInt bd = x.b * y.b;
    return {x.a * y.a - bd, x.a * y.b + x.b * y.a - bd};
  }

  element_type Pow(const element_type& x, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = x;
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  // conj(a + b*w) = a + b*w^2 = (a - b) - b*w.
  element_type Conjugate(const element_type& x) const {
    return {x.a - x.b, -x.b};
  }

  // N(a + b*w) = a^2 - a*b + b^2.
  BigInt Norm(const element_type& x) const {
    return x.a * x.a - x.a * x.b + x.b * x.b;
  }

  BigInt EuclideanNorm(const element_type& x) const { return Norm(x); }

  element_type AdditiveIdentity() const { return {BigInt(0), BigInt(0)}; }

  element_type MultiplicativeIdentity() const { return {BigInt(1), BigInt(0)}; }

  bool IsZero(const element_type& x) const noexcept {
    return x.a.IsZero() && x.b.IsZero();
  }

  bool IsOne(const element_type& x) const noexcept {
    return x.a.IsOne() && x.b.IsZero();
  }

  bool Equals(const element_type& x, const element_type& y) const noexcept {
    return x == y;
  }

  std::string ToString(const element_type& x) const {
    if (x.b.IsZero()) {
      return x.a.ToString();
    }
    if (x.a.IsZero()) {
      return x.b.ToString() + "w";
    }
    if (x.b.IsNegative()) {
      return x.a.ToString() + " - " + (-x.b).ToString() + "w";
    }
    return x.a.ToString() + " + " + x.b.ToString() + "w";
  }

  // Units in Z[w] are the 6 elements with N(x) == 1: +/-1, +/-w, +/-(1+w).
  bool IsUnit(const element_type& x) const { return Norm(x).IsOne(); }

  std::optional<element_type> TryInv(const element_type& x) const {
    if (!IsUnit(x)) {
      return std::nullopt;
    }
    // For N(x) == 1, x * conj(x) == 1.
    return Conjugate(x);
  }

  // Canonical associate in the fundamental sector [0, pi/3): a > b >= 0.
  // Returns {canonical, unit} such that canonical == unit * x.
  std::pair<element_type, element_type> AssociateCanonical(
      const element_type& x) const {
    if (IsZero(x)) {
      return {AdditiveIdentity(), MultiplicativeIdentity()};
    }
    const std::array<element_type, 6> kUnits = {
        element_type{BigInt(1), BigInt(0)},    // 1
        element_type{BigInt(1), BigInt(1)},    // -w^2 = 1 + w
        element_type{BigInt(0), BigInt(1)},    // w
        element_type{BigInt(-1), BigInt(0)},   // -1
        element_type{BigInt(-1), BigInt(-1)},  // w^2 = -1 - w
        element_type{BigInt(0), BigInt(-1)},   // -w
    };
    for (const auto& u : kUnits) {
      element_type cand = Mul(u, x);
      if (!cand.b.IsNegative() && cand.a > cand.b) {
        return {std::move(cand), u};
      }
    }
    // Unreachable for nonzero x; fallback preserves invariant.
    return {x, MultiplicativeIdentity()};
  }

  // Euclidean division in Z[w] using exact 6-neighbor Voronoi cell check
  // around the rounded quotient, guaranteeing N(r) <= floor(N(y)/3) < N(y).
  std::pair<element_type, element_type> DivRem(const element_type& x,
                                               const element_type& y) const {
    if (IsZero(y)) {
      throw std::domain_error(
          "EisensteinIntegerRing::DivRem: division by zero");
    }
    const BigInt norm_y = Norm(y);
    const element_type prod = Mul(x, Conjugate(y));
    const element_type q0{detail::RoundDiv(prod.a, norm_y),
                          detail::RoundDiv(prod.b, norm_y)};

    element_type best_q = q0;
    element_type best_r = Sub(x, Mul(best_q, y));
    BigInt best_norm = Norm(best_r);

    // Check the 6 unit Voronoi neighbors around q0 in the hexagonal lattice.
    static constexpr std::array<std::pair<int64_t, int64_t>, 6> kNeighbors = {
        {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {-1, -1}}};
    for (const auto& [da, db] : kNeighbors) {
      element_type cand_q{q0.a + BigInt(da), q0.b + BigInt(db)};
      element_type cand_r = Sub(x, Mul(cand_q, y));
      BigInt cand_norm = Norm(cand_r);
      if (cand_norm < best_norm) {
        best_q = std::move(cand_q);
        best_r = std::move(cand_r);
        best_norm = std::move(cand_norm);
      }
    }
    return {std::move(best_q), std::move(best_r)};
  }

  // Extended Euclidean Algorithm: returns {g, s, t} with g canonical in the
  // fundamental sector (or 0 if a == b == 0) and s * a + t * b == g.
  std::tuple<element_type, element_type, element_type> ExtendedGcd(
      const element_type& a, const element_type& b) const {
    element_type old_r = a;
    element_type r = b;
    element_type old_s = MultiplicativeIdentity();
    element_type s = AdditiveIdentity();
    element_type old_t = AdditiveIdentity();
    element_type t = MultiplicativeIdentity();

    while (!IsZero(r)) {
      auto [q, rem] = DivRem(old_r, r);
      old_r = std::move(r);
      r = std::move(rem);

      element_type next_s = Sub(old_s, Mul(q, s));
      old_s = std::move(s);
      s = std::move(next_s);

      element_type next_t = Sub(old_t, Mul(q, t));
      old_t = std::move(t);
      t = std::move(next_t);
    }

    auto [g_canon, unit] = AssociateCanonical(old_r);
    return {std::move(g_canon), Mul(unit, old_s), Mul(unit, old_t)};
  }
};

}  // namespace ring

template <>
struct AlgebraTraits<ring::GaussianIntegerRing> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = false;
};

template <>
struct AlgebraTraits<ring::EisensteinIntegerRing> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_QUADRATIC_INTEGERS_HPP_
