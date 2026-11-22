#ifndef XGALOIS_EC_ELLIPTIC_CURVE_HPP_
#define XGALOIS_EC_ELLIPTIC_CURVE_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace ec {

enum class CurveForm : uint8_t {
  kShortWeierstrass,
  kBinaryWeierstrass,
  kMontgomery,
  kTwistedEdwards,
};

// ============================================================================
// Point Representations
// ============================================================================

template <Field F>
struct AffinePoint {
  using field_type = F;
  using element_type = typename F::element_type;

  element_type x{};
  element_type y{};
  bool is_infinity = true;

  constexpr AffinePoint() = default;

  constexpr AffinePoint(element_type x_val, element_type y_val,
                        bool infinity = false)
      : x(std::move(x_val)), y(std::move(y_val)), is_infinity(infinity) {}

  static constexpr AffinePoint Identity() { return AffinePoint{}; }

  bool operator==(const AffinePoint& other) const {
    if (is_infinity && other.is_infinity) {
      return true;
    }
    if (is_infinity != other.is_infinity) {
      return false;
    }
    return x == other.x && y == other.y;
  }

  bool operator!=(const AffinePoint& other) const { return !(*this == other); }
};

// Short Weierstrass Jacobian (X/Z^2, Y/Z^3) or Binary Lopez-Dahab (X/Z, Y/Z^2).
template <Field F>
struct JacobianPoint {
  using field_type = F;
  using element_type = typename F::element_type;

  element_type x{};
  element_type y{};
  element_type z{};
  bool is_infinity = true;

  constexpr JacobianPoint() = default;

  constexpr JacobianPoint(element_type x_val, element_type y_val,
                          element_type z_val, bool infinity = false)
      : x(std::move(x_val)),
        y(std::move(y_val)),
        z(std::move(z_val)),
        is_infinity(infinity) {}

  static constexpr JacobianPoint Identity() { return JacobianPoint{}; }
};

// Twisted Edwards extended coordinates (X : Y : Z : T) with x = X/Z, y = Y/Z,
// and T = X * Y / Z.
template <Field F>
struct ExtendedEdwardsPoint {
  using field_type = F;
  using element_type = typename F::element_type;

  element_type x{};
  element_type y{};
  element_type z{};
  element_type t{};

  constexpr ExtendedEdwardsPoint() = default;

  constexpr ExtendedEdwardsPoint(element_type x_val, element_type y_val,
                                 element_type z_val, element_type t_val)
      : x(std::move(x_val)),
        y(std::move(y_val)),
        z(std::move(z_val)),
        t(std::move(t_val)) {}
};

// ============================================================================
// Internal Field & Number-Theory Helpers
// ============================================================================

namespace detail {

template <typename F>
inline BigInt FieldBigCharacteristic(const F& field) {
  if constexpr (requires { field.BigCharacteristic(); }) {
    return field.BigCharacteristic();
  } else if constexpr (requires { field.Characteristic(); }) {
    return BigInt::FromUint64(static_cast<uint64_t>(field.Characteristic()));
  } else {
    return BigInt(0);
  }
}

template <typename F>
inline std::size_t FieldAbsDegree(const F& field) {
  if constexpr (requires { field.AbsoluteDegree(); }) {
    return static_cast<std::size_t>(field.AbsoluteDegree());
  } else if constexpr (requires { field.Degree(); }) {
    return static_cast<std::size_t>(field.Degree());
  } else {
    return 1;
  }
}

template <typename F>
inline BigInt FieldBigOrder(const F& field) {
  if constexpr (requires { field.BigOrder(); }) {
    return field.BigOrder();
  } else if constexpr (requires { field.Order(); }) {
    const uint64_t ord = static_cast<uint64_t>(field.Order());
    if (ord > 0) {
      return BigInt::FromUint64(ord);
    }
    const BigInt ch = FieldBigCharacteristic(field);
    if (ch > BigInt(0)) {
      return ch.Pow(static_cast<uint64_t>(FieldAbsDegree(field)));
    }
    return BigInt(0);
  } else {
    return BigInt(0);
  }
}

template <Field F>
inline typename F::element_type FromSmallInt(const F& field, uint64_t n) {
  using E = typename F::element_type;
  if (n == 0) {
    return field.AdditiveIdentity();
  }
  E acc = field.AdditiveIdentity();
  E base = field.MultiplicativeIdentity();
  while (n > 0) {
    if ((n & 1ULL) != 0) {
      acc = field.Add(acc, base);
    }
    n >>= 1;
    if (n > 0) {
      base = field.Add(base, base);
    }
  }
  return acc;
}

template <Field F>
inline typename F::element_type FromBigInt(const F& field, const BigInt& n) {
  using E = typename F::element_type;
  if (n.IsZero()) {
    return field.AdditiveIdentity();
  }
  if (n.IsNegative()) {
    return field.Neg(FromBigInt(field, -n));
  }
  E acc = field.AdditiveIdentity();
  E base = field.MultiplicativeIdentity();
  const std::size_t bits = n.BitLength();
  for (std::size_t i = 0; i < bits; ++i) {
    if (n.GetBit(i)) {
      acc = field.Add(acc, base);
    }
    if (i + 1 < bits) {
      base = field.Add(base, base);
    }
  }
  return acc;
}

template <Field F>
inline typename F::element_type FieldSqr(const F& field,
                                         const typename F::element_type& a) {
  if constexpr (requires { field.Sqr(a); }) {
    return field.Sqr(a);
  } else {
    return field.Mul(a, a);
  }
}

template <Field F>
inline typename F::element_type FieldPowBig(const F& field,
                                            const typename F::element_type& a,
                                            const BigInt& exp) {
  using E = typename F::element_type;
  if constexpr (requires {
                  { field.PowBig(a, exp) } -> std::same_as<E>;
                }) {
    return field.PowBig(a, exp);
  } else {
    if (exp.IsZero()) {
      return field.MultiplicativeIdentity();
    }
    if (exp.IsNegative()) {
      return FieldPowBig(field, field.Inv(a), -exp);
    }
    if (field.IsZero(a)) {
      return field.AdditiveIdentity();
    }
    if (field.IsOne(a)) {
      return field.MultiplicativeIdentity();
    }
    E result = field.MultiplicativeIdentity();
    E base = a;
    const std::size_t bits = exp.BitLength();
    for (std::size_t i = 0; i < bits; ++i) {
      if (exp.GetBit(i)) {
        result = field.Mul(result, base);
      }
      if (i + 1 < bits) {
        base = FieldSqr(field, base);
      }
    }
    return result;
  }
}

// Deterministic pseudo-random field element generator from a uint64_t seed.
// Works for prime fields, BigInt prime fields, binary extensions, and recursive
// FieldTower extensions.
template <Field F>
inline typename F::element_type SampleFieldElement(const F& field,
                                                   uint64_t& state) {
  using E = typename F::element_type;
  auto next_u64 = [&state]() noexcept -> uint64_t {
    state += 0x9e3779b97f4a7c15ULL;
    uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  };

  if constexpr (requires {
                  typename F::base_field_type;
                  field.BaseField();
                  field.RelativeDegree();
                  field.FromCoeffs(
                      std::vector<typename F::base_element_type>{});
                }) {
    const std::size_t deg = field.RelativeDegree();
    std::vector<typename F::base_element_type> coeffs;
    coeffs.reserve(deg);
    for (std::size_t i = 0; i < deg; ++i) {
      coeffs.push_back(SampleFieldElement(field.BaseField(), state));
    }
    return field.FromCoeffs(coeffs);
  } else if constexpr (requires {
                         field.Modulus();
                         field.Reduce(std::declval<E>());
                       } && std::same_as<E, BigInt>) {
    const BigInt& p = field.Modulus();
    BigInt val(0);
    for (std::size_t i = 0; i < p.LimbCount(); ++i) {
      val = (val << 64) + BigInt::FromUint64(next_u64());
    }
    return field.Reduce(val);
  } else if constexpr (requires {
                         F::kBits;
                         field.Modulus();
                         field.Reduce(std::declval<E>());
                       }) {
    BigInt val(0);
    const std::size_t limbs = F::kBits / 64;
    for (std::size_t i = 0; i < limbs; ++i) {
      val = (val << 64) + BigInt::FromUint64(next_u64());
    }
    const BigInt p = BigInt::FromBigUInt(field.Modulus());
    return (val % p).template ToBigUInt<F::kBits>();
  } else if constexpr (std::unsigned_integral<E>) {
    const BigInt q = FieldBigOrder(field);
    if (q > BigInt(0) && q.LimbCount() == 1) {
      const uint64_t q64 = q.ToUint64();
      return static_cast<E>(next_u64() % q64);
    }
    return static_cast<E>(next_u64());
  } else {
    return FromSmallInt(field, next_u64());
  }
}

template <Field F>
inline bool FieldHasSqrt(const F& field, const typename F::element_type& a) {
  if (field.IsZero(a)) {
    return true;
  }
  if constexpr (requires {
                  { field.HasSqrt(a) } -> std::convertible_to<bool>;
                }) {
    return static_cast<bool>(field.HasSqrt(a));
  } else {
    const BigInt ch = FieldBigCharacteristic(field);
    if (ch == BigInt(2)) {
      return true;
    }
    const BigInt q = FieldBigOrder(field);
    if (q <= BigInt(0)) {
      throw std::domain_error(
          "FieldHasSqrt requires a finite field with known order");
    }
    const BigInt exp = (q - BigInt(1)) >> 1;
    return field.IsOne(FieldPowBig(field, a, exp));
  }
}

template <Field F>
inline std::optional<typename F::element_type> FieldSqrt(
    const F& field, const typename F::element_type& a) {
  using E = typename F::element_type;
  if (field.IsZero(a)) {
    return field.AdditiveIdentity();
  }
  if constexpr (requires {
                  { field.Sqrt(a) } -> std::same_as<std::optional<E>>;
                }) {
    return field.Sqrt(a);
  } else {
    const BigInt q = FieldBigOrder(field);
    if (q <= BigInt(0)) {
      throw std::domain_error(
          "FieldSqrt requires a finite field with known order");
    }
    const BigInt ch = FieldBigCharacteristic(field);
    if (ch == BigInt(2)) {
      return FieldPowBig(field, a, q >> 1);
    }

    const BigInt q_minus_1 = q - BigInt(1);
    const BigInt half = q_minus_1 >> 1;
    if (!field.IsOne(FieldPowBig(field, a, half))) {
      return std::nullopt;
    }

    // Fast path for q == 3 (mod 4).
    if ((q.Limb(0) & 3ULL) == 3ULL) {
      return FieldPowBig(field, a, (q + BigInt(1)) >> 2);
    }

    // General Tonelli-Shanks over any odd-order finite field F_q.
    const std::size_t s = q_minus_1.CountTrailingZeroBits();
    const BigInt odd_q = q_minus_1 >> s;

    uint64_t seed = 0x123456789abcdef0ULL;
    E z = field.AdditiveIdentity();
    while (true) {
      E cand = SampleFieldElement(field, seed);
      if (field.IsZero(cand)) {
        continue;
      }
      if (!field.IsOne(FieldPowBig(field, cand, half))) {
        z = cand;
        break;
      }
    }

    std::size_t m = s;
    E c = FieldPowBig(field, z, odd_q);
    E t = FieldPowBig(field, a, odd_q);
    E r = FieldPowBig(field, a, (odd_q + BigInt(1)) >> 1);

    while (!field.IsOne(t)) {
      std::size_t i = 0;
      E tmp = t;
      while (!field.IsOne(tmp)) {
        tmp = FieldSqr(field, tmp);
        ++i;
        if (i >= m) {
          return std::nullopt;
        }
      }
      E b = c;
      for (std::size_t j = 0; j < m - i - 1; ++j) {
        b = FieldSqr(field, b);
      }
      m = i;
      c = FieldSqr(field, b);
      t = field.Mul(t, c);
      r = field.Mul(r, b);
    }
    return r;
  }
}

// Computes signed w-NAF representation of a non-negative BigInt scalar.
// Digits are in {0, +/-1, +/-3, ..., +/-(2^(w-1) - 1)} in LSB-first order.
inline std::vector<int8_t> ComputeWnaf(const BigInt& scalar, int w = 4) {
  std::vector<int8_t> naf;
  if (scalar.IsZero()) {
    return naf;
  }
  const int64_t mod = int64_t{1} << w;
  const int64_t half = mod >> 1;
  BigInt k = scalar;
  naf.reserve(k.BitLength() + 1);
  while (!k.IsZero()) {
    if (k.IsOdd()) {
      int64_t d =
          static_cast<int64_t>(k.Limb(0) & static_cast<uint64_t>(mod - 1));
      if (d >= half) {
        d -= mod;
      }
      naf.push_back(static_cast<int8_t>(d));
      if (d > 0) {
        k -= BigInt::FromInt64(d);
      } else {
        k += BigInt::FromInt64(-d);
      }
    } else {
      naf.push_back(0);
    }
    k >>= 1;
  }
  return naf;
}

inline uint32_t ExtractWindowBits(const BigInt& scalar, std::size_t bit_offset,
                                  std::size_t width) noexcept {
  uint32_t val = 0;
  for (std::size_t b = 0; b < width; ++b) {
    if (scalar.GetBit(bit_offset + b)) {
      val |= (uint32_t{1} << b);
    }
  }
  return val;
}

template <Field F>
inline std::string CanonicalPointKey(const F& field, const AffinePoint<F>& pt) {
  if (pt.is_infinity) {
    return "O";
  }
  return field.ToString(pt.x) + "," + field.ToString(pt.y);
}

inline BigInt BigLcm(const BigInt& a, const BigInt& b) {
  if (a.IsZero() || b.IsZero()) {
    return BigInt(0);
  }
  return (a / Gcd(a, b)) * b;
}

// Finds the unique multiple of `step` in [lo, hi], or std::nullopt if 0 or >=2
// multiples lie in the interval.
inline std::optional<BigInt> UniqueMultipleInInterval(const BigInt& step,
                                                      const BigInt& lo,
                                                      const BigInt& hi) {
  if (step <= BigInt(0) || lo > hi) {
    return std::nullopt;
  }
  const BigInt first = ((lo + step - BigInt(1)) / step) * step;
  if (first > hi) {
    return std::nullopt;
  }
  if (first + step <= hi) {
    return std::nullopt;
  }
  return first;
}

// Solves x == r1 (mod m1) and x == r2 (mod m2), returning {x mod lcm, lcm} or
// std::nullopt if inconsistent.
inline std::optional<std::pair<BigInt, BigInt>> CombineCongruences(
    const BigInt& r1, const BigInt& m1, const BigInt& r2, const BigInt& m2) {
  const auto [g, u, v] = ExtendedGcd(m1, m2);
  (void)v;
  const BigInt diff = r2 - r1;
  const auto [q, rem] = BigInt::EuclideanDivRem(diff, g);
  if (!rem.IsZero()) {
    return std::nullopt;
  }
  const BigInt m2_div_g = m2 / g;
  const BigInt k = BigInt::EuclideanDivRem(q * u, m2_div_g).second;
  const BigInt lcm = m1 * m2_div_g;
  const BigInt x = BigInt::EuclideanDivRem(r1 + m1 * k, lcm).second;
  return std::make_pair(x, lcm);
}

inline std::optional<BigInt> UniqueCongruentInInterval(const BigInt& rem,
                                                       const BigInt& mod,
                                                       const BigInt& lo,
                                                       const BigInt& hi) {
  if (mod <= BigInt(0) || lo > hi) {
    return std::nullopt;
  }
  const BigInt r = BigInt::EuclideanDivRem(rem, mod).second;
  const BigInt lo_rem = BigInt::EuclideanDivRem(lo, mod).second;
  BigInt first = lo + (r - lo_rem);
  if (first < lo) {
    first += mod;
  }
  if (first > hi) {
    return std::nullopt;
  }
  if (first + mod <= hi) {
    return std::nullopt;
  }
  return first;
}

}  // namespace detail

// ============================================================================
// EllipticCurve<F, Form>
// ============================================================================

template <Field F, CurveForm Form = CurveForm::kShortWeierstrass>
class EllipticCurve {
 public:
  static constexpr CurveForm kForm = Form;
  using field_type = F;
  using field_element_type = typename F::element_type;
  using element_type = AffinePoint<F>;

  // Constructs a curve over `field`:
  // - kShortWeierstrass: y^2 = x^3 + param1*x + param2  (a, b)
  // - kBinaryWeierstrass: y^2 + xy = x^3 + param1*x^2 + param2  (a, b)
  // - kMontgomery: param2*y^2 = x^3 + param1*x^2 + x  (A, B)
  // - kTwistedEdwards: param1*x^2 + y^2 = 1 + param2*x^2*y^2  (a, d)
  EllipticCurve(F field, field_element_type param1, field_element_type param2)
      : field_(std::move(field)),
        p1_(std::move(param1)),
        p2_(std::move(param2)) {
    ValidateAndPrecompute();
  }

  EllipticCurve(const std::shared_ptr<F>& field, field_element_type param1,
                field_element_type param2)
      : EllipticCurve(RequireField(field), std::move(param1),
                      std::move(param2)) {}

  const F& Field() const noexcept { return field_; }
  const field_element_type& A() const noexcept { return p1_; }
  const field_element_type& B() const noexcept { return p2_; }
  const field_element_type& D() const noexcept
    requires(Form == CurveForm::kTwistedEdwards)
  {
    return p2_;
  }

  // --------------------------------------------------------------------------
  // Group & Curve Primitives
  // --------------------------------------------------------------------------

  AffinePoint<F> Identity() const {
    if constexpr (Form == CurveForm::kTwistedEdwards) {
      return AffinePoint<F>(field_.AdditiveIdentity(),
                            field_.MultiplicativeIdentity(), false);
    } else {
      return AffinePoint<F>(field_.AdditiveIdentity(),
                            field_.AdditiveIdentity(), true);
    }
  }

  bool IsIdentity(const AffinePoint<F>& p) const {
    if constexpr (Form == CurveForm::kTwistedEdwards) {
      return p.is_infinity || (field_.IsZero(p.x) && field_.IsOne(p.y));
    } else {
      return p.is_infinity;
    }
  }

  bool Equals(const AffinePoint<F>& p, const AffinePoint<F>& q) const {
    const bool p_id = IsIdentity(p);
    const bool q_id = IsIdentity(q);
    if (p_id && q_id) {
      return true;
    }
    if (p_id != q_id) {
      return false;
    }
    return alg::Equals(field_, p.x, q.x) && alg::Equals(field_, p.y, q.y);
  }

  bool IsOnCurve(const AffinePoint<F>& p) const {
    if (p.is_infinity) {
      return true;
    }
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      const auto y2 = detail::FieldSqr(field_, p.y);
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto x3 = field_.Mul(x2, p.x);
      const auto rhs = field_.Add(field_.Add(x3, field_.Mul(p1_, p.x)), p2_);
      return alg::Equals(field_, y2, rhs);
    } else if constexpr (Form == CurveForm::kBinaryWeierstrass) {
      const auto y2 = detail::FieldSqr(field_, p.y);
      const auto xy = field_.Mul(p.x, p.y);
      const auto lhs = field_.Add(y2, xy);
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto x3 = field_.Mul(x2, p.x);
      const auto rhs = field_.Add(field_.Add(x3, field_.Mul(p1_, x2)), p2_);
      return alg::Equals(field_, lhs, rhs);
    } else if constexpr (Form == CurveForm::kMontgomery) {
      const auto y2 = detail::FieldSqr(field_, p.y);
      const auto lhs = field_.Mul(p2_, y2);
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto x3 = field_.Mul(x2, p.x);
      const auto rhs = field_.Add(field_.Add(x3, field_.Mul(p1_, x2)), p.x);
      return alg::Equals(field_, lhs, rhs);
    } else {
      // Twisted Edwards: a*x^2 + y^2 = 1 + d*x^2*y^2
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto y2 = detail::FieldSqr(field_, p.y);
      const auto lhs = field_.Add(field_.Mul(p1_, x2), y2);
      const auto rhs = field_.Add(field_.MultiplicativeIdentity(),
                                  field_.Mul(p2_, field_.Mul(x2, y2)));
      return alg::Equals(field_, lhs, rhs);
    }
  }

  AffinePoint<F> Neg(const AffinePoint<F>& p) const {
    if (IsIdentity(p)) {
      return Identity();
    }
    if constexpr (Form == CurveForm::kShortWeierstrass ||
                  Form == CurveForm::kMontgomery) {
      return AffinePoint<F>(p.x, field_.Neg(p.y), false);
    } else if constexpr (Form == CurveForm::kBinaryWeierstrass) {
      return AffinePoint<F>(p.x, field_.Add(p.x, p.y), false);
    } else {
      return AffinePoint<F>(field_.Neg(p.x), p.y, false);
    }
  }

  AffinePoint<F> Double(const AffinePoint<F>& p) const {
    if (IsIdentity(p)) {
      return Identity();
    }
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      if (field_.IsZero(p.y)) {
        return Identity();
      }
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto three_x2 = field_.Add(field_.Add(x2, x2), x2);
      const auto num = field_.Add(three_x2, p1_);
      const auto den = field_.Add(p.y, p.y);
      const auto lam = field_.Div(num, den);
      const auto x3 =
          field_.Sub(detail::FieldSqr(field_, lam), field_.Add(p.x, p.x));
      const auto y3 = field_.Sub(field_.Mul(lam, field_.Sub(p.x, x3)), p.y);
      return AffinePoint<F>(x3, y3, false);
    } else if constexpr (Form == CurveForm::kBinaryWeierstrass) {
      if (field_.IsZero(p.x)) {
        return Identity();
      }
      const auto lam = field_.Add(p.x, field_.Div(p.y, p.x));
      const auto x3 =
          field_.Add(field_.Add(detail::FieldSqr(field_, lam), lam), p1_);
      const auto lam_plus_1 = field_.Add(lam, field_.MultiplicativeIdentity());
      const auto y3 =
          field_.Add(detail::FieldSqr(field_, p.x), field_.Mul(lam_plus_1, x3));
      return AffinePoint<F>(x3, y3, false);
    } else if constexpr (Form == CurveForm::kMontgomery) {
      if (field_.IsZero(p.y)) {
        return Identity();
      }
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto three_x2 = field_.Add(field_.Add(x2, x2), x2);
      const auto two_ax = field_.Mul(field_.Add(p1_, p1_), p.x);
      const auto num = field_.Add(field_.Add(three_x2, two_ax),
                                  field_.MultiplicativeIdentity());
      const auto den = field_.Mul(field_.Add(p2_, p2_), p.y);
      const auto lam = field_.Div(num, den);
      const auto b_lam2 = field_.Mul(p2_, detail::FieldSqr(field_, lam));
      const auto x3 = field_.Sub(field_.Sub(b_lam2, p1_), field_.Add(p.x, p.x));
      const auto y3 = field_.Sub(field_.Mul(lam, field_.Sub(p.x, x3)), p.y);
      return AffinePoint<F>(x3, y3, false);
    } else {
      return ToAffine(DoubleExtended(ToExtended(p)));
    }
  }

  AffinePoint<F> Add(const AffinePoint<F>& p, const AffinePoint<F>& q) const {
    if (IsIdentity(p)) {
      return IsIdentity(q) ? Identity() : q;
    }
    if (IsIdentity(q)) {
      return p;
    }
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      if (alg::Equals(field_, p.x, q.x)) {
        if (alg::Equals(field_, p.y, q.y)) {
          return Double(p);
        }
        return Identity();
      }
      const auto lam = field_.Div(field_.Sub(q.y, p.y), field_.Sub(q.x, p.x));
      const auto x3 =
          field_.Sub(field_.Sub(detail::FieldSqr(field_, lam), p.x), q.x);
      const auto y3 = field_.Sub(field_.Mul(lam, field_.Sub(p.x, x3)), p.y);
      return AffinePoint<F>(x3, y3, false);
    } else if constexpr (Form == CurveForm::kBinaryWeierstrass) {
      if (alg::Equals(field_, p.x, q.x)) {
        if (alg::Equals(field_, p.y, q.y)) {
          return Double(p);
        }
        return Identity();
      }
      const auto lam = field_.Div(field_.Add(q.y, p.y), field_.Add(q.x, p.x));
      const auto x3 =
          field_.Add(field_.Add(field_.Add(detail::FieldSqr(field_, lam), lam),
                                field_.Add(p.x, q.x)),
                     p1_);
      const auto y3 =
          field_.Add(field_.Add(field_.Mul(lam, field_.Add(p.x, x3)), x3), p.y);
      return AffinePoint<F>(x3, y3, false);
    } else if constexpr (Form == CurveForm::kMontgomery) {
      if (alg::Equals(field_, p.x, q.x)) {
        if (alg::Equals(field_, p.y, q.y)) {
          return Double(p);
        }
        return Identity();
      }
      const auto lam = field_.Div(field_.Sub(q.y, p.y), field_.Sub(q.x, p.x));
      const auto b_lam2 = field_.Mul(p2_, detail::FieldSqr(field_, lam));
      const auto x3 = field_.Sub(field_.Sub(field_.Sub(b_lam2, p1_), p.x), q.x);
      const auto y3 = field_.Sub(field_.Mul(lam, field_.Sub(p.x, x3)), p.y);
      return AffinePoint<F>(x3, y3, false);
    } else {
      return ToAffine(AddExtended(ToExtended(p), ToExtended(q)));
    }
  }

  AffinePoint<F> Sub(const AffinePoint<F>& p, const AffinePoint<F>& q) const {
    return Add(p, Neg(q));
  }

  // xg::AbelianGroup interface aliases.
  AffinePoint<F> Op(const AffinePoint<F>& p, const AffinePoint<F>& q) const {
    return Add(p, q);
  }

  AffinePoint<F> Inv(const AffinePoint<F>& p) const { return Neg(p); }

  AffinePoint<F> Div(const AffinePoint<F>& p, const AffinePoint<F>& q) const {
    return Sub(p, q);
  }

  AffinePoint<F> Pow(const AffinePoint<F>& p, uint64_t exp) const {
    return ScalarMul(p, BigInt::FromUint64(exp));
  }

  std::string ToString(const AffinePoint<F>& p) const {
    if (IsIdentity(p)) {
      return "O";
    }
    return "(" + field_.ToString(p.x) + ", " + field_.ToString(p.y) + ")";
  }

  // --------------------------------------------------------------------------
  // Projective / Jacobian / Extended Coordinate Operations
  // --------------------------------------------------------------------------

  JacobianPoint<F> ToJacobian(const AffinePoint<F>& p) const
    requires(Form == CurveForm::kShortWeierstrass ||
             Form == CurveForm::kBinaryWeierstrass)
  {
    if (p.is_infinity) {
      return JacobianPoint<F>(field_.MultiplicativeIdentity(),
                              field_.MultiplicativeIdentity(),
                              field_.AdditiveIdentity(), true);
    }
    return JacobianPoint<F>(p.x, p.y, field_.MultiplicativeIdentity(), false);
  }

  AffinePoint<F> ToAffine(const JacobianPoint<F>& j) const
    requires(Form == CurveForm::kShortWeierstrass ||
             Form == CurveForm::kBinaryWeierstrass)
  {
    if (j.is_infinity || field_.IsZero(j.z)) {
      return Identity();
    }
    const auto z_inv = field_.Inv(j.z);
    const auto z_inv2 = detail::FieldSqr(field_, z_inv);
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      const auto z_inv3 = field_.Mul(z_inv2, z_inv);
      return AffinePoint<F>(field_.Mul(j.x, z_inv2), field_.Mul(j.y, z_inv3),
                            false);
    } else {
      return AffinePoint<F>(field_.Mul(j.x, z_inv), field_.Mul(j.y, z_inv2),
                            false);
    }
  }

  // Montgomery's trick: converts N Jacobian/Lopez-Dahab points to affine using
  // a single field inversion.
  std::vector<AffinePoint<F>> BatchToAffine(
      std::span<const JacobianPoint<F>> pts) const
    requires(Form == CurveForm::kShortWeierstrass ||
             Form == CurveForm::kBinaryWeierstrass)
  {
    const std::size_t n = pts.size();
    std::vector<AffinePoint<F>> result(n, Identity());
    std::vector<std::size_t> active;
    active.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      if (!pts[i].is_infinity && !field_.IsZero(pts[i].z)) {
        active.push_back(i);
      }
    }
    const std::size_t m = active.size();
    if (m == 0) {
      return result;
    }

    std::vector<field_element_type> prefix(m);
    prefix[0] = pts[active[0]].z;
    for (std::size_t k = 1; k < m; ++k) {
      prefix[k] = field_.Mul(prefix[k - 1], pts[active[k]].z);
    }

    auto inv_acc = field_.Inv(prefix[m - 1]);
    for (std::size_t k = m; k-- > 0;) {
      const std::size_t idx = active[k];
      const auto z_inv = (k > 0) ? field_.Mul(inv_acc, prefix[k - 1]) : inv_acc;
      if (k > 0) {
        inv_acc = field_.Mul(inv_acc, pts[idx].z);
      }
      const auto z_inv2 = detail::FieldSqr(field_, z_inv);
      if constexpr (Form == CurveForm::kShortWeierstrass) {
        const auto z_inv3 = field_.Mul(z_inv2, z_inv);
        result[idx] = AffinePoint<F>(field_.Mul(pts[idx].x, z_inv2),
                                     field_.Mul(pts[idx].y, z_inv3), false);
      } else {
        result[idx] = AffinePoint<F>(field_.Mul(pts[idx].x, z_inv),
                                     field_.Mul(pts[idx].y, z_inv2), false);
      }
    }
    return result;
  }

  JacobianPoint<F> DoubleJacobian(const JacobianPoint<F>& p) const
    requires(Form == CurveForm::kShortWeierstrass ||
             Form == CurveForm::kBinaryWeierstrass)
  {
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      if (p.is_infinity || field_.IsZero(p.z) || field_.IsZero(p.y)) {
        return JacobianPoint<F>::Identity();
      }
      const auto xx = detail::FieldSqr(field_, p.x);
      const auto yy = detail::FieldSqr(field_, p.y);
      const auto yyyy = detail::FieldSqr(field_, yy);
      const auto zz = detail::FieldSqr(field_, p.z);

      const auto x_yy = field_.Mul(p.x, yy);
      const auto two_x_yy = field_.Add(x_yy, x_yy);
      const auto s = field_.Add(two_x_yy, two_x_yy);  // 4 * X * Y^2

      const auto three_xx = field_.Add(field_.Add(xx, xx), xx);
      const auto m =
          field_.IsZero(p1_)
              ? three_xx
              : field_.Add(three_xx,
                           field_.Mul(p1_, detail::FieldSqr(field_, zz)));

      const auto x3 = field_.Sub(detail::FieldSqr(field_, m), field_.Add(s, s));
      const auto two_yyyy = field_.Add(yyyy, yyyy);
      const auto four_yyyy = field_.Add(two_yyyy, two_yyyy);
      const auto eight_yyyy = field_.Add(four_yyyy, four_yyyy);
      const auto y3 = field_.Sub(field_.Mul(m, field_.Sub(s, x3)), eight_yyyy);
      const auto yz = field_.Mul(p.y, p.z);
      const auto z3 = field_.Add(yz, yz);
      return JacobianPoint<F>(x3, y3, z3, false);
    } else {
      // Binary Lopez-Dahab doubling:
      // Z3 = X1^2 * Z1^2, X3 = X1^4 + b*Z1^4,
      // Y3 = b*Z1^4*Z3 + X3*(a*Z3 + Y1^2 + b*Z1^4).
      if (p.is_infinity || field_.IsZero(p.z) || field_.IsZero(p.x)) {
        return JacobianPoint<F>::Identity();
      }
      const auto x2 = detail::FieldSqr(field_, p.x);
      const auto z2 = detail::FieldSqr(field_, p.z);
      const auto z3 = field_.Mul(x2, z2);
      const auto x4 = detail::FieldSqr(field_, x2);
      const auto z4 = detail::FieldSqr(field_, z2);
      const auto bz4 = field_.Mul(p2_, z4);
      const auto x3 = field_.Add(x4, bz4);
      const auto y2 = detail::FieldSqr(field_, p.y);
      const auto term = field_.Add(field_.Add(field_.Mul(p1_, z3), y2), bz4);
      const auto y3 = field_.Add(field_.Mul(bz4, z3), field_.Mul(x3, term));
      return JacobianPoint<F>(x3, y3, z3, false);
    }
  }

  JacobianPoint<F> AddJacobian(const JacobianPoint<F>& p,
                               const JacobianPoint<F>& q) const
    requires(Form == CurveForm::kShortWeierstrass ||
             Form == CurveForm::kBinaryWeierstrass)
  {
    if (p.is_infinity || field_.IsZero(p.z)) {
      return q;
    }
    if (q.is_infinity || field_.IsZero(q.z)) {
      return p;
    }
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      const auto z1z1 = detail::FieldSqr(field_, p.z);
      const auto z2z2 = detail::FieldSqr(field_, q.z);
      const auto u1 = field_.Mul(p.x, z2z2);
      const auto u2 = field_.Mul(q.x, z1z1);
      const auto s1 = field_.Mul(p.y, field_.Mul(q.z, z2z2));
      const auto s2 = field_.Mul(q.y, field_.Mul(p.z, z1z1));

      if (alg::Equals(field_, u1, u2)) {
        if (alg::Equals(field_, s1, s2)) {
          return DoubleJacobian(p);
        }
        return JacobianPoint<F>::Identity();
      }

      const auto h = field_.Sub(u2, u1);
      const auto r = field_.Sub(s2, s1);
      const auto hh = detail::FieldSqr(field_, h);
      const auto hhh = field_.Mul(h, hh);
      const auto v = field_.Mul(u1, hh);

      const auto x3 = field_.Sub(field_.Sub(detail::FieldSqr(field_, r), hhh),
                                 field_.Add(v, v));
      const auto y3 =
          field_.Sub(field_.Mul(r, field_.Sub(v, x3)), field_.Mul(s1, hhh));
      const auto z3 = field_.Mul(h, field_.Mul(p.z, q.z));
      return JacobianPoint<F>(x3, y3, z3, false);
    } else {
      // Binary Lopez-Dahab general addition via affine fallback for simplicity.
      return ToJacobian(Add(ToAffine(p), ToAffine(q)));
    }
  }

  JacobianPoint<F> AddMixedJacobian(const JacobianPoint<F>& p,
                                    const AffinePoint<F>& q) const
    requires(Form == CurveForm::kShortWeierstrass ||
             Form == CurveForm::kBinaryWeierstrass)
  {
    if (q.is_infinity) {
      return p;
    }
    if (p.is_infinity || field_.IsZero(p.z)) {
      return ToJacobian(q);
    }
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      const auto z1z1 = detail::FieldSqr(field_, p.z);
      const auto u2 = field_.Mul(q.x, z1z1);
      const auto s2 = field_.Mul(q.y, field_.Mul(p.z, z1z1));
      if (alg::Equals(field_, p.x, u2)) {
        if (alg::Equals(field_, p.y, s2)) {
          return DoubleJacobian(p);
        }
        return JacobianPoint<F>::Identity();
      }
      const auto h = field_.Sub(u2, p.x);
      const auto r = field_.Sub(s2, p.y);
      const auto hh = detail::FieldSqr(field_, h);
      const auto hhh = field_.Mul(h, hh);
      const auto v = field_.Mul(p.x, hh);
      const auto x3 = field_.Sub(field_.Sub(detail::FieldSqr(field_, r), hhh),
                                 field_.Add(v, v));
      const auto y3 =
          field_.Sub(field_.Mul(r, field_.Sub(v, x3)), field_.Mul(p.y, hhh));
      const auto z3 = field_.Mul(p.z, h);
      return JacobianPoint<F>(x3, y3, z3, false);
    } else {
      // Lopez-Dahab mixed addition (Z2 = 1).
      const auto z1z1 = detail::FieldSqr(field_, p.z);
      const auto a_val = field_.Add(field_.Mul(q.y, z1z1), p.y);
      const auto b_val = field_.Add(field_.Mul(q.x, p.z), p.x);
      if (field_.IsZero(b_val)) {
        if (field_.IsZero(a_val)) {
          return DoubleJacobian(p);
        }
        return JacobianPoint<F>::Identity();
      }
      const auto c_val = field_.Mul(p.z, b_val);
      const auto d_val = field_.Mul(detail::FieldSqr(field_, b_val),
                                    field_.Add(c_val, field_.Mul(p1_, z1z1)));
      const auto z3 = detail::FieldSqr(field_, c_val);
      const auto e_val = field_.Mul(a_val, c_val);
      const auto x3 =
          field_.Add(field_.Add(detail::FieldSqr(field_, a_val), d_val), e_val);
      const auto f_val = field_.Add(x3, field_.Mul(q.x, z3));
      const auto g_val = field_.Add(x3, field_.Mul(q.y, z3));
      const auto y3 =
          field_.Add(field_.Mul(e_val, f_val), field_.Mul(z3, g_val));
      return JacobianPoint<F>(x3, y3, z3, false);
    }
  }

  ExtendedEdwardsPoint<F> ToExtended(const AffinePoint<F>& p) const
    requires(Form == CurveForm::kTwistedEdwards)
  {
    if (IsIdentity(p)) {
      return ExtendedEdwardsPoint<F>(
          field_.AdditiveIdentity(), field_.MultiplicativeIdentity(),
          field_.MultiplicativeIdentity(), field_.AdditiveIdentity());
    }
    return ExtendedEdwardsPoint<F>(p.x, p.y, field_.MultiplicativeIdentity(),
                                   field_.Mul(p.x, p.y));
  }

  AffinePoint<F> ToAffine(const ExtendedEdwardsPoint<F>& e) const
    requires(Form == CurveForm::kTwistedEdwards)
  {
    if (field_.IsZero(e.z)) {
      return Identity();
    }
    const auto z_inv = field_.Inv(e.z);
    return AffinePoint<F>(field_.Mul(e.x, z_inv), field_.Mul(e.y, z_inv),
                          false);
  }

  ExtendedEdwardsPoint<F> AddExtended(const ExtendedEdwardsPoint<F>& p1,
                                      const ExtendedEdwardsPoint<F>& p2) const
    requires(Form == CurveForm::kTwistedEdwards)
  {
    const auto a = field_.Mul(p1.x, p2.x);
    const auto b = field_.Mul(p1.y, p2.y);
    const auto c = field_.Mul(p2_, field_.Mul(p1.t, p2.t));
    const auto d = field_.Mul(p1.z, p2.z);
    const auto e = field_.Sub(
        field_.Sub(field_.Mul(field_.Add(p1.x, p1.y), field_.Add(p2.x, p2.y)),
                   a),
        b);
    const auto f = field_.Sub(d, c);
    const auto g = field_.Add(d, c);
    const auto h = field_.Sub(b, field_.Mul(p1_, a));
    return ExtendedEdwardsPoint<F>(field_.Mul(e, f), field_.Mul(g, h),
                                   field_.Mul(f, g), field_.Mul(e, h));
  }

  ExtendedEdwardsPoint<F> DoubleExtended(
      const ExtendedEdwardsPoint<F>& p1) const
    requires(Form == CurveForm::kTwistedEdwards)
  {
    const auto a = detail::FieldSqr(field_, p1.x);
    const auto b = detail::FieldSqr(field_, p1.y);
    const auto z2 = detail::FieldSqr(field_, p1.z);
    const auto c = field_.Add(z2, z2);
    const auto d = field_.Mul(p1_, a);
    const auto e = field_.Sub(
        field_.Sub(detail::FieldSqr(field_, field_.Add(p1.x, p1.y)), a), b);
    const auto g = field_.Add(d, b);
    const auto f = field_.Sub(g, c);
    const auto h = field_.Sub(d, b);
    return ExtendedEdwardsPoint<F>(field_.Mul(e, f), field_.Mul(g, h),
                                   field_.Mul(f, g), field_.Mul(e, h));
  }

  // --------------------------------------------------------------------------
  // Montgomery XZ Ladder
  // --------------------------------------------------------------------------

  field_element_type LadderMulX(const field_element_type& x_u,
                                const BigInt& k) const
    requires(Form == CurveForm::kMontgomery)
  {
    const BigInt k_abs = k.Abs();
    if (k_abs.IsZero()) {
      return field_.AdditiveIdentity();
    }
    const auto [x2, z2, x3, z3] = MontgomeryLadderXZ(x_u, k_abs);
    (void)x3;
    (void)z3;
    if (field_.IsZero(z2)) {
      return field_.AdditiveIdentity();
    }
    return field_.Mul(x2, field_.Inv(z2));
  }

  // --------------------------------------------------------------------------
  // Scalar Multiplication & Multi-Scalar Multiplication
  // --------------------------------------------------------------------------

  AffinePoint<F> ScalarMul(const AffinePoint<F>& p, const BigInt& k) const {
    if (k.IsZero() || IsIdentity(p)) {
      return Identity();
    }
    if (k.IsNegative()) {
      return ScalarMul(Neg(p), -k);
    }
    if (k.IsOne()) {
      return p;
    }

    if constexpr (Form == CurveForm::kMontgomery) {
      return ScalarMulMontgomery(p, k);
    } else if constexpr (Form == CurveForm::kTwistedEdwards) {
      // Signed w-NAF (w = 4) in Extended Edwards coordinates.
      const auto naf = detail::ComputeWnaf(k, 4);
      const auto p_ext = ToExtended(p);
      const auto p2 = DoubleExtended(p_ext);
      // Precompute odd multiples: table[j] = (2*j + 1) * P for j = 0..3.
      ExtendedEdwardsPoint<F> pos[4];
      ExtendedEdwardsPoint<F> neg[4];
      pos[0] = p_ext;
      for (int j = 1; j < 4; ++j) {
        pos[j] = AddExtended(pos[j - 1], p2);
      }
      for (int j = 0; j < 4; ++j) {
        neg[j] = ExtendedEdwardsPoint<F>(field_.Neg(pos[j].x), pos[j].y,
                                         pos[j].z, field_.Neg(pos[j].t));
      }
      auto acc = ToExtended(Identity());
      for (std::size_t i = naf.size(); i-- > 0;) {
        acc = DoubleExtended(acc);
        const int8_t d = naf[i];
        if (d > 0) {
          acc = AddExtended(acc, pos[(d - 1) >> 1]);
        } else if (d < 0) {
          acc = AddExtended(acc, neg[((-d) - 1) >> 1]);
        }
      }
      return ToAffine(acc);
    } else {
      // Signed w-NAF (w = 4) in Jacobian / Lopez-Dahab coordinates.
      const auto naf = detail::ComputeWnaf(k, 4);
      const auto p_jac = ToJacobian(p);
      const auto p2 = DoubleJacobian(p_jac);
      JacobianPoint<F> odd_jac[4];
      odd_jac[0] = p_jac;
      for (int j = 1; j < 4; ++j) {
        odd_jac[j] = AddJacobian(odd_jac[j - 1], p2);
      }
      // Convert odd multiples to affine with a single inversion for fast mixed
      // additions!
      const auto pos =
          BatchToAffine(std::span<const JacobianPoint<F>>(odd_jac, 4));
      AffinePoint<F> neg[4];
      for (int j = 0; j < 4; ++j) {
        neg[j] = Neg(pos[j]);
      }

      auto acc = JacobianPoint<F>::Identity();
      for (std::size_t i = naf.size(); i-- > 0;) {
        acc = DoubleJacobian(acc);
        const int8_t d = naf[i];
        if (d > 0) {
          acc = AddMixedJacobian(acc, pos[(d - 1) >> 1]);
        } else if (d < 0) {
          acc = AddMixedJacobian(acc, neg[((-d) - 1) >> 1]);
        }
      }
      return ToAffine(acc);
    }
  }

  // Constant-time scalar multiplication (Montgomery ladder with uniform
  // add+double per bit).
  AffinePoint<F> ScalarMulConstantTime(const AffinePoint<F>& p,
                                       const BigInt& k) const {
    if (k.IsZero() || IsIdentity(p)) {
      return Identity();
    }
    if (k.IsNegative()) {
      return ScalarMulConstantTime(Neg(p), -k);
    }
    if constexpr (Form == CurveForm::kMontgomery) {
      return ScalarMulMontgomery(p, k);
    } else if constexpr (Form == CurveForm::kTwistedEdwards) {
      auto r0 = ToExtended(Identity());
      auto r1 = ToExtended(p);
      const std::size_t bits = k.BitLength();
      uint8_t prev_bit = 0;
      for (std::size_t i = bits; i-- > 0;) {
        const uint8_t bit = k.GetBit(i) ? 1U : 0U;
        if ((prev_bit ^ bit) != 0U) {
          std::swap(r0, r1);
        }
        r1 = AddExtended(r0, r1);
        r0 = DoubleExtended(r0);
        prev_bit = bit;
      }
      if (prev_bit != 0U) {
        std::swap(r0, r1);
      }
      return ToAffine(r0);
    } else {
      auto r0 = JacobianPoint<F>::Identity();
      auto r1 = ToJacobian(p);
      const std::size_t bits = k.BitLength();
      uint8_t prev_bit = 0;
      for (std::size_t i = bits; i-- > 0;) {
        const uint8_t bit = k.GetBit(i) ? 1U : 0U;
        if ((prev_bit ^ bit) != 0U) {
          std::swap(r0, r1);
        }
        r1 = AddJacobian(r0, r1);
        r0 = DoubleJacobian(r0);
        prev_bit = bit;
      }
      if (prev_bit != 0U) {
        std::swap(r0, r1);
      }
      return ToAffine(r0);
    }
  }

  // Multi-scalar multiplication: sum_i scalars[i] * pts[i].
  // Uses Pippenger's bucket algorithm for N >= 16 and Straus/Shamir for N < 16.
  AffinePoint<F> MultiScalarMul(std::span<const AffinePoint<F>> pts,
                                std::span<const BigInt> scalars) const {
    if (pts.size() != scalars.size()) {
      throw std::invalid_argument(
          "MultiScalarMul: points and scalars spans must have equal length");
    }
    std::vector<AffinePoint<F>> clean_pts;
    std::vector<BigInt> clean_scalars;
    clean_pts.reserve(pts.size());
    clean_scalars.reserve(scalars.size());
    std::size_t max_bits = 0;

    for (std::size_t i = 0; i < pts.size(); ++i) {
      if (scalars[i].IsZero() || IsIdentity(pts[i])) {
        continue;
      }
      if (scalars[i].IsNegative()) {
        clean_pts.push_back(Neg(pts[i]));
        clean_scalars.push_back(-scalars[i]);
      } else {
        clean_pts.push_back(pts[i]);
        clean_scalars.push_back(scalars[i]);
      }
      max_bits = std::max(max_bits, clean_scalars.back().BitLength());
    }

    const std::size_t n = clean_pts.size();
    if (n == 0) {
      return Identity();
    }
    if (n == 1) {
      return ScalarMul(clean_pts[0], clean_scalars[0]);
    }

    if constexpr (Form == CurveForm::kMontgomery) {
      // Montgomery affine accumulation via Straus / individual ladders.
      AffinePoint<F> acc = Identity();
      for (std::size_t i = 0; i < n; ++i) {
        acc = Add(acc, ScalarMul(clean_pts[i], clean_scalars[i]));
      }
      return acc;
    } else if constexpr (Form == CurveForm::kTwistedEdwards) {
      if (n < 16) {
        // Straus / Shamir simultaneous binary double-and-add.
        std::vector<ExtendedEdwardsPoint<F>> ext_pts;
        ext_pts.reserve(n);
        for (const auto& pt : clean_pts) {
          ext_pts.push_back(ToExtended(pt));
        }
        auto acc = ToExtended(Identity());
        for (std::size_t bit = max_bits; bit-- > 0;) {
          acc = DoubleExtended(acc);
          for (std::size_t i = 0; i < n; ++i) {
            if (clean_scalars[i].GetBit(bit)) {
              acc = AddExtended(acc, ext_pts[i]);
            }
          }
        }
        return ToAffine(acc);
      }

      // Pippenger bucket algorithm in Extended Edwards coordinates.
      std::size_t c = 2;
      while ((std::size_t{1} << (c + 1)) <= n && c < 12) {
        ++c;
      }
      const std::size_t num_buckets = (std::size_t{1} << c) - 1;
      const std::size_t num_windows = (max_bits + c - 1) / c;
      const auto ext_id = ToExtended(Identity());

      std::vector<ExtendedEdwardsPoint<F>> ext_pts;
      ext_pts.reserve(n);
      for (const auto& pt : clean_pts) {
        ext_pts.push_back(ToExtended(pt));
      }

      auto total = ext_id;
      for (std::size_t w = num_windows; w-- > 0;) {
        if (w + 1 < num_windows) {
          for (std::size_t d = 0; d < c; ++d) {
            total = DoubleExtended(total);
          }
        }
        std::vector<ExtendedEdwardsPoint<F>> buckets(num_buckets, ext_id);
        const std::size_t bit_offset = w * c;
        for (std::size_t i = 0; i < n; ++i) {
          const uint32_t digit =
              detail::ExtractWindowBits(clean_scalars[i], bit_offset, c);
          if (digit > 0) {
            buckets[digit - 1] = AddExtended(buckets[digit - 1], ext_pts[i]);
          }
        }
        auto running = ext_id;
        auto window_sum = ext_id;
        for (std::size_t b = num_buckets; b-- > 0;) {
          running = AddExtended(running, buckets[b]);
          window_sum = AddExtended(window_sum, running);
        }
        total = AddExtended(total, window_sum);
      }
      return ToAffine(total);
    } else {
      // Short Weierstrass or Binary Weierstrass in Jacobian/Lopez-Dahab.
      if (n < 16) {
        auto acc = JacobianPoint<F>::Identity();
        for (std::size_t bit = max_bits; bit-- > 0;) {
          acc = DoubleJacobian(acc);
          for (std::size_t i = 0; i < n; ++i) {
            if (clean_scalars[i].GetBit(bit)) {
              acc = AddMixedJacobian(acc, clean_pts[i]);
            }
          }
        }
        return ToAffine(acc);
      }

      // Pippenger bucket algorithm.
      std::size_t c = 2;
      while ((std::size_t{1} << (c + 1)) <= n && c < 12) {
        ++c;
      }
      const std::size_t num_buckets = (std::size_t{1} << c) - 1;
      const std::size_t num_windows = (max_bits + c - 1) / c;
      const auto jac_id = JacobianPoint<F>::Identity();

      auto total = jac_id;
      for (std::size_t w = num_windows; w-- > 0;) {
        if (w + 1 < num_windows) {
          for (std::size_t d = 0; d < c; ++d) {
            total = DoubleJacobian(total);
          }
        }
        std::vector<JacobianPoint<F>> buckets(num_buckets, jac_id);
        const std::size_t bit_offset = w * c;
        for (std::size_t i = 0; i < n; ++i) {
          const uint32_t digit =
              detail::ExtractWindowBits(clean_scalars[i], bit_offset, c);
          if (digit > 0) {
            buckets[digit - 1] =
                AddMixedJacobian(buckets[digit - 1], clean_pts[i]);
          }
        }
        auto running = jac_id;
        auto window_sum = jac_id;
        for (std::size_t b = num_buckets; b-- > 0;) {
          if (!buckets[b].is_infinity) {
            running = AddJacobian(running, buckets[b]);
          }
          if (!running.is_infinity) {
            window_sum = AddJacobian(window_sum, running);
          }
        }
        if (!window_sum.is_infinity) {
          total = AddJacobian(total, window_sum);
        }
      }
      return ToAffine(total);
    }
  }

  // --------------------------------------------------------------------------
  // Short Weierstrass Invariants & Number-Theoretic Methods
  // --------------------------------------------------------------------------

  field_element_type Discriminant() const
    requires(Form == CurveForm::kShortWeierstrass)
  {
    const auto a2 = detail::FieldSqr(field_, p1_);
    const auto a3 = field_.Mul(a2, p1_);
    const auto b2 = detail::FieldSqr(field_, p2_);
    const auto four_a3 = field_.Mul(detail::FromSmallInt(field_, 4), a3);
    const auto twenty_seven_b2 =
        field_.Mul(detail::FromSmallInt(field_, 27), b2);
    const auto inner = field_.Add(four_a3, twenty_seven_b2);
    return field_.Neg(field_.Mul(detail::FromSmallInt(field_, 16), inner));
  }

  field_element_type JInvariant() const
    requires(Form == CurveForm::kShortWeierstrass)
  {
    const auto a2 = detail::FieldSqr(field_, p1_);
    const auto a3 = field_.Mul(a2, p1_);
    const auto b2 = detail::FieldSqr(field_, p2_);
    const auto four_a3 = field_.Mul(detail::FromSmallInt(field_, 4), a3);
    const auto twenty_seven_b2 =
        field_.Mul(detail::FromSmallInt(field_, 27), b2);
    const auto den = field_.Add(four_a3, twenty_seven_b2);
    const auto num = field_.Mul(detail::FromSmallInt(field_, 1728), four_a3);
    return field_.Div(num, den);
  }

  // Recursive division polynomial psi_n(x, y) on y^2 = x^3 + a*x + b.
  field_element_type DivisionPolynomial(uint32_t n, const field_element_type& x,
                                        const field_element_type& y) const
    requires(Form == CurveForm::kShortWeierstrass)
  {
    if (n == 0) {
      return field_.AdditiveIdentity();
    }
    if (n == 1) {
      return field_.MultiplicativeIdentity();
    }
    const auto two_y = field_.Add(y, y);
    if (n == 2) {
      return two_y;
    }

    const auto x2 = detail::FieldSqr(field_, x);
    const auto x3 = field_.Mul(x2, x);
    const auto x4 = detail::FieldSqr(field_, x2);
    const auto x6 = field_.Mul(x4, x2);
    const auto a2 = detail::FieldSqr(field_, p1_);
    const auto a3 = field_.Mul(a2, p1_);
    const auto b2 = detail::FieldSqr(field_, p2_);
    const auto ab = field_.Mul(p1_, p2_);

    // f2 = (2*y)^2 = 4*y^2, f2_sq = (2*y)^4.
    const auto f2 = detail::FieldSqr(field_, two_y);
    const auto f2_sq = detail::FieldSqr(field_, f2);

    // Reduced division polynomial psi_tilde_n where psi_n = psi_tilde_n (odd n)
    // and psi_n = 2*y * psi_tilde_n (even n).
    std::unordered_map<uint32_t, field_element_type> memo;
    memo.emplace(0U, field_.AdditiveIdentity());
    memo.emplace(1U, field_.MultiplicativeIdentity());
    memo.emplace(2U, field_.MultiplicativeIdentity());

    // psi_3 = 3*x^4 + 6*a*x^2 + 12*b*x - a^2
    const auto psi3 = field_.Sub(
        field_.Add(
            field_.Add(field_.Mul(detail::FromSmallInt(field_, 3), x4),
                       field_.Mul(detail::FromSmallInt(field_, 6),
                                  field_.Mul(p1_, x2))),
            field_.Mul(detail::FromSmallInt(field_, 12), field_.Mul(p2_, x))),
        a2);
    memo.emplace(3U, psi3);

    // psi_tilde_4 = 2*(x^6 + 5*a*x^4 + 20*b*x^3 - 5*a^2*x^2 - 4*a*b*x - 8*b^2 -
    // a^3)
    const auto inner4 = field_.Sub(
        field_.Sub(
            field_.Sub(
                field_.Sub(
                    field_.Add(field_.Add(x6, field_.Mul(detail::FromSmallInt(
                                                             field_, 5),
                                                         field_.Mul(p1_, x4))),
                               field_.Mul(detail::FromSmallInt(field_, 20),
                                          field_.Mul(p2_, x3))),
                    field_.Mul(detail::FromSmallInt(field_, 5),
                               field_.Mul(a2, x2))),
                field_.Mul(detail::FromSmallInt(field_, 4), field_.Mul(ab, x))),
            field_.Mul(detail::FromSmallInt(field_, 8), b2)),
        a3);
    memo.emplace(4U, field_.Add(inner4, inner4));

    std::function<field_element_type(uint32_t)> eval_tilde =
        [&](uint32_t k) -> field_element_type {
      const auto it = memo.find(k);
      if (it != memo.end()) {
        return it->second;
      }
      field_element_type val;
      if ((k & 1U) == 0U) {
        const uint32_t m = k >> 1;
        const auto pm = eval_tilde(m);
        const auto pm_p2 = eval_tilde(m + 2);
        const auto pm_m1 = eval_tilde(m - 1);
        const auto pm_m2 = eval_tilde(m - 2);
        const auto pm_p1 = eval_tilde(m + 1);
        const auto t1 = field_.Mul(pm_p2, detail::FieldSqr(field_, pm_m1));
        const auto t2 = field_.Mul(pm_m2, detail::FieldSqr(field_, pm_p1));
        val = field_.Mul(pm, field_.Sub(t1, t2));
      } else {
        const uint32_t m = k >> 1;
        const auto pm = eval_tilde(m);
        const auto pm_p2 = eval_tilde(m + 2);
        const auto pm_m1 = eval_tilde(m - 1);
        const auto pm_p1 = eval_tilde(m + 1);
        const auto pm3 = field_.Mul(detail::FieldSqr(field_, pm), pm);
        const auto pmp1_3 = field_.Mul(detail::FieldSqr(field_, pm_p1), pm_p1);
        auto t1 = field_.Mul(pm_p2, pm3);
        auto t2 = field_.Mul(pm_m1, pmp1_3);
        if ((m & 1U) == 1U) {
          t2 = field_.Mul(f2_sq, t2);
        } else {
          t1 = field_.Mul(f2_sq, t1);
        }
        val = field_.Sub(t1, t2);
      }
      memo.emplace(k, val);
      return val;
    };

    const auto res = eval_tilde(n);
    return ((n & 1U) == 0U) ? field_.Mul(two_y, res) : res;
  }

  // Exact order of a point P on E(F_q). If group_order_hint > 0, factors the
  // hint and reduces; otherwise runs Hasse-interval BSGS.
  BigInt PointOrder(const AffinePoint<F>& p,
                    const BigInt& group_order_hint = BigInt(0)) const
    requires(Form == CurveForm::kShortWeierstrass)
  {
    if (IsIdentity(p)) {
      return BigInt(1);
    }
    if (group_order_hint > BigInt(0)) {
      if (!IsIdentity(ScalarMul(p, group_order_hint))) {
        throw std::invalid_argument(
            "PointOrder: group_order_hint does not annihilate P");
      }
      return ReduceMultipleToExactOrder(p, group_order_hint);
    }

    const BigInt q = detail::FieldBigOrder(field_);
    if (q <= BigInt(0)) {
      throw std::domain_error(
          "PointOrder without hint requires a finite field with known order");
    }

    const BigInt h_bound = ISqrt(q * BigInt(4)) + BigInt(1);
    const BigInt center = q + BigInt(1);
    // Baby-step table size s = ceil(sqrt(2 * h_bound + 1)).
    const BigInt s_big = ISqrt(h_bound * BigInt(2) + BigInt(1)) + BigInt(1);
    const uint64_t s = s_big.ToUint64();

    std::unordered_map<std::string, uint64_t> baby;
    baby.reserve(static_cast<std::size_t>(s) + 1);
    AffinePoint<F> cur = Identity();
    baby.emplace(detail::CanonicalPointKey(field_, cur), 0ULL);
    for (uint64_t j = 1; j < s; ++j) {
      cur = Add(cur, p);
      if (IsIdentity(cur)) {
        return BigInt::FromUint64(j);
      }
      baby.emplace(detail::CanonicalPointKey(field_, cur), j);
    }

    // Giant steps from Q = (q + 1) * P with step S = s * P.
    const AffinePoint<F> q_pt = ScalarMul(p, center);
    const AffinePoint<F> step_pt = Add(cur, p);  // s * P
    const AffinePoint<F> neg_step_pt = Neg(step_pt);

    const uint64_t max_i = (h_bound / s_big).ToUint64() + 2;
    AffinePoint<F> g_pos = q_pt;
    AffinePoint<F> g_neg = q_pt;

    for (uint64_t i = 0; i <= max_i; ++i) {
      // Check g_pos = (q + 1 + i*s) * P == j * P  =>  M = q + 1 + i*s - j.
      auto it = baby.find(detail::CanonicalPointKey(field_, g_pos));
      if (it != baby.end()) {
        const BigInt m = center + BigInt::FromUint64(i) * s_big -
                         BigInt::FromUint64(it->second);
        if (m > BigInt(0) && IsIdentity(ScalarMul(p, m))) {
          return ReduceMultipleToExactOrder(p, m);
        }
      }
      // Check g_pos == -j * P  =>  M = q + 1 + i*s + j.
      auto it_neg = baby.find(detail::CanonicalPointKey(field_, Neg(g_pos)));
      if (it_neg != baby.end()) {
        const BigInt m = center + BigInt::FromUint64(i) * s_big +
                         BigInt::FromUint64(it_neg->second);
        if (m > BigInt(0) && IsIdentity(ScalarMul(p, m))) {
          return ReduceMultipleToExactOrder(p, m);
        }
      }

      if (i > 0) {
        // Check g_neg = (q + 1 - i*s) * P == +/- j * P.
        auto it2 = baby.find(detail::CanonicalPointKey(field_, g_neg));
        if (it2 != baby.end()) {
          const BigInt m = center - BigInt::FromUint64(i) * s_big -
                           BigInt::FromUint64(it2->second);
          if (m > BigInt(0) && IsIdentity(ScalarMul(p, m))) {
            return ReduceMultipleToExactOrder(p, m);
          }
        }
        auto it2_neg = baby.find(detail::CanonicalPointKey(field_, Neg(g_neg)));
        if (it2_neg != baby.end()) {
          const BigInt m = center - BigInt::FromUint64(i) * s_big +
                           BigInt::FromUint64(it2_neg->second);
          if (m > BigInt(0) && IsIdentity(ScalarMul(p, m))) {
            return ReduceMultipleToExactOrder(p, m);
          }
        }
      }

      g_pos = Add(g_pos, step_pt);
      g_neg = Add(g_neg, neg_step_pt);
    }

    throw std::runtime_error("PointOrder: BSGS failed to locate point order");
  }

  // Exact point counting #E(F_q) for Short Weierstrass curves over finite
  // fields: direct Legendre-symbol sum for small prime fields (q <= 10000), and
  // deterministic Shanks-Mestre BSGS with quadratic twist fallback otherwise.
  BigInt CurveOrder() const
    requires(Form == CurveForm::kShortWeierstrass)
  {
    const BigInt q = detail::FieldBigOrder(field_);
    if (q <= BigInt(0)) {
      throw std::domain_error(
          "CurveOrder requires a finite field with known order");
    }

    // Fast exact Legendre-symbol sum for small prime fields (deg == 1, q <=
    // 10000).
    if (detail::FieldAbsDegree(field_) == 1 && q <= BigInt(10000)) {
      const uint64_t p = q.ToUint64();
      uint64_t count = 1;  // Point at infinity.
      field_element_type x = field_.AdditiveIdentity();
      const field_element_type one = field_.MultiplicativeIdentity();
      for (uint64_t i = 0; i < p; ++i) {
        const auto x2 = detail::FieldSqr(field_, x);
        const auto x3 = field_.Mul(x2, x);
        const auto rhs = field_.Add(field_.Add(x3, field_.Mul(p1_, x)), p2_);
        if (field_.IsZero(rhs)) {
          count += 1;
        } else if (detail::FieldHasSqrt(field_, rhs)) {
          count += 2;
        }
        x = field_.Add(x, one);
      }
      return BigInt::FromUint64(count);
    }

    // Shanks-Mestre BSGS across Hasse interval [q + 1 - H, q + 1 + H] using
    // both E and its quadratic twist E^d.
    const BigInt h_bound = ISqrt(q * BigInt(4)) + BigInt(1);
    const BigInt lo = q + BigInt(1) - h_bound;
    const BigInt hi = q + BigInt(1) + h_bound;
    const BigInt two_q_plus_2 = q * BigInt(2) + BigInt(2);

    // Find a quadratic non-residue d in F_q for the quadratic twist E^d.
    uint64_t seed = 0x517cc1b727220a95ULL;
    field_element_type d = field_.AdditiveIdentity();
    while (true) {
      const auto cand = detail::SampleFieldElement(field_, seed);
      if (!field_.IsZero(cand) && !detail::FieldHasSqrt(field_, cand)) {
        d = cand;
        break;
      }
    }
    const auto d2 = detail::FieldSqr(field_, d);
    const auto d3 = field_.Mul(d2, d);
    const EllipticCurve<F, CurveForm::kShortWeierstrass> twist(
        field_, field_.Mul(p1_, d2), field_.Mul(p2_, d3));

    BigInt lcm_e(1);
    BigInt lcm_t(1);

    for (uint64_t attempt = 0; attempt < 256; ++attempt) {
      const auto x = detail::SampleFieldElement(field_, seed);
      const auto x2 = detail::FieldSqr(field_, x);
      const auto x3 = field_.Mul(x2, x);
      const auto rhs = field_.Add(field_.Add(x3, field_.Mul(p1_, x)), p2_);

      const auto y_opt = detail::FieldSqrt(field_, rhs);
      if (y_opt.has_value()) {
        const AffinePoint<F> pt(x, *y_opt, false);
        const BigInt ord = PointOrder(pt);
        lcm_e = detail::BigLcm(lcm_e, ord);
        if (auto uniq = detail::UniqueMultipleInInterval(lcm_e, lo, hi)) {
          return *uniq;
        }
      } else {
        // d * rhs is a square, giving point (d*x, sqrt(d^3 * rhs)) on E^d.
        const auto rhs_t = field_.Mul(d3, rhs);
        const auto yt_opt = detail::FieldSqrt(field_, rhs_t);
        if (yt_opt.has_value()) {
          const AffinePoint<F> pt_t(field_.Mul(d, x), *yt_opt, false);
          const BigInt ord_t = twist.PointOrder(pt_t);
          lcm_t = detail::BigLcm(lcm_t, ord_t);
          if (auto uniq_t = detail::UniqueMultipleInInterval(lcm_t, lo, hi)) {
            return two_q_plus_2 - *uniq_t;
          }
        }
      }

      // Combine congruences N == 0 (mod lcm_e) and N == 2q+2 (mod lcm_t).
      if (lcm_e > BigInt(1) && lcm_t > BigInt(1)) {
        if (auto comb = detail::CombineCongruences(BigInt(0), lcm_e,
                                                   two_q_plus_2, lcm_t)) {
          if (auto uniq = detail::UniqueCongruentInInterval(
                  comb->first, comb->second, lo, hi)) {
            return *uniq;
          }
        }
      }
    }

    throw std::runtime_error(
        "CurveOrder: failed to uniquely determine curve order");
  }

  // Smallest k >= 1 such that r | (q^k - 1), bounded by max_k (returns 0 if
  // none found within max_k).
  uint64_t EmbeddingDegree(const BigInt& r, uint64_t max_k = 64) const {
    if (r <= BigInt(1)) {
      throw std::invalid_argument("EmbeddingDegree: r must be > 1");
    }
    const BigInt q = detail::FieldBigOrder(field_);
    if (q <= BigInt(0)) {
      throw std::domain_error(
          "EmbeddingDegree requires a finite field with known order");
    }
    if (!Gcd(q, r).IsOne()) {
      throw std::invalid_argument("EmbeddingDegree: r must be coprime to q");
    }
    const BigInt q_mod = q % r;
    BigInt cur = q_mod;
    for (uint64_t k = 1; k <= max_k; ++k) {
      if (cur.IsOne()) {
        return k;
      }
      cur = (cur * q_mod) % r;
    }
    return 0;
  }

  // Finds two independent r-torsion points (P, Q) forming a basis of E[r] over
  // the current field when E[r] is contained in E(F_q).
  std::pair<AffinePoint<F>, AffinePoint<F>> TorsionBasis(
      const BigInt& r, const BigInt& curve_order_hint = BigInt(0)) const
    requires(Form == CurveForm::kShortWeierstrass)
  {
    if (r <= BigInt(1)) {
      throw std::invalid_argument("TorsionBasis: r must be > 1");
    }
    const BigInt n =
        (curve_order_hint > BigInt(0)) ? curve_order_hint : CurveOrder();
    const BigInt r2 = r * r;
    if (!(n % r2).IsZero()) {
      throw std::domain_error(
          "TorsionBasis: r^2 does not divide #E(F_q); E[r] not rational over "
          "field");
    }

    // Factor out prime factors of r from n to obtain cofactor m coprime to r.
    const auto r_factors = FactorInteger(r);
    BigInt m = n;
    for (const auto& [prime, exp] : r_factors) {
      (void)exp;
      while ((m % prime).IsZero()) {
        m /= prime;
      }
    }

    auto has_exact_order_r = [&](const AffinePoint<F>& pt) -> bool {
      if (IsIdentity(pt) || !IsIdentity(ScalarMul(pt, r))) {
        return false;
      }
      for (const auto& [prime, exp] : r_factors) {
        (void)exp;
        if (IsIdentity(ScalarMul(pt, r / prime))) {
          return false;
        }
      }
      return true;
    };

    // Given P and Q of exact order r, checks <P> \cap <Q> == {O} by verifying
    // that for every prime l | r, (r/l)*Q is not in <(r/l)*P>.
    auto are_independent = [&](const AffinePoint<F>& p,
                               const AffinePoint<F>& q) -> bool {
      for (const auto& [prime, exp] : r_factors) {
        (void)exp;
        const BigInt cofactor = r / prime;
        const AffinePoint<F> lp = ScalarMul(p, cofactor);
        const AffinePoint<F> lq = ScalarMul(q, cofactor);
        // Check if lq in <lp> via BSGS of size sqrt(prime).
        const BigInt s_big = ISqrt(prime) + BigInt(1);
        const uint64_t s = s_big.ToUint64();
        std::unordered_map<std::string, uint64_t> table;
        table.reserve(static_cast<std::size_t>(s));
        AffinePoint<F> cur = Identity();
        for (uint64_t j = 0; j < s; ++j) {
          table.emplace(detail::CanonicalPointKey(field_, cur), j);
          cur = Add(cur, lp);
        }
        const AffinePoint<F> neg_step = Neg(cur);  // -s * lp
        AffinePoint<F> gamma = lq;
        for (uint64_t i = 0; i <= s; ++i) {
          if (table.find(detail::CanonicalPointKey(field_, gamma)) !=
              table.end()) {
            return false;
          }
          gamma = Add(gamma, neg_step);
        }
      }
      return true;
    };

    auto project_to_r_torsion =
        [&](const AffinePoint<F>& raw) -> std::optional<AffinePoint<F>> {
      AffinePoint<F> cur = ScalarMul(raw, m);
      if (IsIdentity(cur)) {
        return std::nullopt;
      }
      // Push down into E[r] by multiplying by prime factors of r while order
      // still exceeds r.
      while (!IsIdentity(ScalarMul(cur, r))) {
        bool progressed = false;
        for (const auto& [prime, exp] : r_factors) {
          (void)exp;
          const AffinePoint<F> nxt = ScalarMul(cur, prime);
          if (!IsIdentity(nxt)) {
            cur = nxt;
            progressed = true;
            break;
          }
        }
        if (!progressed) {
          return std::nullopt;
        }
      }
      if (has_exact_order_r(cur)) {
        return cur;
      }
      return std::nullopt;
    };

    uint64_t seed = 0xa0761d6478bd642fULL;
    std::optional<AffinePoint<F>> p_opt;
    for (uint64_t attempt = 0; attempt < 4096; ++attempt) {
      const auto x = detail::SampleFieldElement(field_, seed);
      const auto x2 = detail::FieldSqr(field_, x);
      const auto rhs =
          field_.Add(field_.Add(field_.Mul(x2, x), field_.Mul(p1_, x)), p2_);
      const auto y_opt = detail::FieldSqrt(field_, rhs);
      if (!y_opt.has_value()) {
        continue;
      }
      const auto cand_opt =
          project_to_r_torsion(AffinePoint<F>(x, *y_opt, false));
      if (!cand_opt.has_value()) {
        continue;
      }
      if (!p_opt.has_value()) {
        p_opt = *cand_opt;
      } else if (are_independent(*p_opt, *cand_opt)) {
        return {*p_opt, *cand_opt};
      }
    }

    throw std::runtime_error(
        "TorsionBasis: failed to find two independent r-torsion points");
  }

 private:
  static F RequireField(const std::shared_ptr<F>& ptr) {
    if (!ptr) {
      throw std::invalid_argument(
          "EllipticCurve: field pointer cannot be null");
    }
    return *ptr;
  }

  void ValidateAndPrecompute() {
    const BigInt ch = detail::FieldBigCharacteristic(field_);
    if constexpr (Form == CurveForm::kShortWeierstrass) {
      if (ch == BigInt(2) || ch == BigInt(3)) {
        throw std::invalid_argument(
            "EllipticCurve<kShortWeierstrass>: field characteristic must not "
            "be 2 or 3");
      }
      if (field_.IsZero(Discriminant())) {
        throw std::invalid_argument(
            "EllipticCurve<kShortWeierstrass>: curve is singular (Delta == 0)");
      }
    } else if constexpr (Form == CurveForm::kBinaryWeierstrass) {
      if (ch != BigInt(2)) {
        throw std::invalid_argument(
            "EllipticCurve<kBinaryWeierstrass>: field must have characteristic "
            "2");
      }
      if (field_.IsZero(p2_)) {
        throw std::invalid_argument(
            "EllipticCurve<kBinaryWeierstrass>: parameter b must be non-zero");
      }
    } else if constexpr (Form == CurveForm::kMontgomery) {
      if (ch == BigInt(2)) {
        throw std::invalid_argument(
            "EllipticCurve<kMontgomery>: field characteristic must not be 2");
      }
      const auto four = detail::FromSmallInt(field_, 4);
      const auto disc_factor = field_.Sub(detail::FieldSqr(field_, p1_), four);
      if (field_.IsZero(p2_) || field_.IsZero(disc_factor)) {
        throw std::invalid_argument(
            "EllipticCurve<kMontgomery>: curve is singular (B*(A^2 - 4) == 0)");
      }
      // Precompute a24 = (A - 2) / 4 for RFC 7748 Montgomery XZ ladder.
      const auto two = detail::FromSmallInt(field_, 2);
      mont_a24_ = field_.Div(field_.Sub(p1_, two), four);
    } else {
      if (ch == BigInt(2)) {
        throw std::invalid_argument(
            "EllipticCurve<kTwistedEdwards>: field characteristic must not be "
            "2");
      }
      if (field_.IsZero(p1_) || field_.IsZero(p2_) ||
          alg::Equals(field_, p1_, p2_)) {
        throw std::invalid_argument(
            "EllipticCurve<kTwistedEdwards>: require a != 0, d != 0, a != d");
      }
    }
  }

  struct MontgomeryXZResult {
    field_element_type x2;
    field_element_type z2;
    field_element_type x3;
    field_element_type z3;
  };

  MontgomeryXZResult MontgomeryLadderXZ(const field_element_type& x_u,
                                        const BigInt& k_pos) const
    requires(Form == CurveForm::kMontgomery)
  {
    field_element_type x2 = field_.MultiplicativeIdentity();
    field_element_type z2 = field_.AdditiveIdentity();
    field_element_type x3 = x_u;
    field_element_type z3 = field_.MultiplicativeIdentity();

    uint8_t swap = 0;
    const std::size_t bits = k_pos.BitLength();
    for (std::size_t t = bits; t-- > 0;) {
      const uint8_t k_t = k_pos.GetBit(t) ? 1U : 0U;
      swap ^= k_t;
      if (swap != 0U) {
        std::swap(x2, x3);
        std::swap(z2, z3);
      }
      swap = k_t;

      const auto a0 = field_.Add(x2, z2);
      const auto aa = detail::FieldSqr(field_, a0);
      const auto b0 = field_.Sub(x2, z2);
      const auto bb = detail::FieldSqr(field_, b0);
      const auto e0 = field_.Sub(aa, bb);
      const auto c0 = field_.Add(x3, z3);
      const auto d0 = field_.Sub(x3, z3);
      const auto da = field_.Mul(d0, a0);
      const auto cb = field_.Mul(c0, b0);

      x3 = detail::FieldSqr(field_, field_.Add(da, cb));
      z3 = field_.Mul(x_u, detail::FieldSqr(field_, field_.Sub(da, cb)));
      x2 = field_.Mul(aa, bb);
      z2 = field_.Mul(e0, field_.Add(aa, field_.Mul(mont_a24_, e0)));
    }
    if (swap != 0U) {
      std::swap(x2, x3);
      std::swap(z2, z3);
    }
    return {x2, z2, x3, z3};
  }

  // Full affine scalar multiplication on Montgomery curves via XZ ladder +
  // Okeya-Sakurai y-coordinate recovery.
  AffinePoint<F> ScalarMulMontgomery(const AffinePoint<F>& p,
                                     const BigInt& k_pos) const
    requires(Form == CurveForm::kMontgomery)
  {
    if (field_.IsZero(p.y)) {
      // 2-torsion point: 2*P == O.
      return k_pos.IsOdd() ? p : Identity();
    }
    const auto [x2_proj, z2_proj, x3_proj, z3_proj] =
        MontgomeryLadderXZ(p.x, k_pos);
    if (field_.IsZero(z2_proj)) {
      return Identity();
    }
    if (field_.IsZero(z3_proj)) {
      // (k + 1) * P == O  =>  k * P == -P.
      return Neg(p);
    }

    const auto x2 = field_.Div(x2_proj, z2_proj);
    const auto x3 = field_.Div(x3_proj, z3_proj);
    const auto one = field_.MultiplicativeIdentity();
    const auto two_a = field_.Add(p1_, p1_);
    const auto term1 = field_.Mul(field_.Add(field_.Mul(p.x, x2), one),
                                  field_.Add(field_.Add(p.x, x2), two_a));
    const auto dx = field_.Sub(p.x, x2);
    const auto term2 = field_.Mul(detail::FieldSqr(field_, dx), x3);
    const auto num = field_.Sub(field_.Sub(term1, two_a), term2);
    const auto den = field_.Mul(field_.Add(p2_, p2_), p.y);
    const auto y2 = field_.Div(num, den);
    return AffinePoint<F>(x2, y2, false);
  }

  BigInt ReduceMultipleToExactOrder(const AffinePoint<F>& p,
                                    BigInt mult) const {
    const auto factors = FactorInteger(mult);
    for (const auto& [prime, exp] : factors) {
      for (uint64_t e = 0; e < exp; ++e) {
        const BigInt cand = mult / prime;
        if (IsIdentity(ScalarMul(p, cand))) {
          mult = cand;
        } else {
          break;
        }
      }
    }
    return mult;
  }

  F field_;
  field_element_type p1_{};
  field_element_type p2_{};
  field_element_type mont_a24_{};
};

}  // namespace ec

// ============================================================================
// AlgebraTraits specialization in namespace xg
// ============================================================================

template <Field F, ec::CurveForm Form>
struct AlgebraTraits<ec::EllipticCurve<F, Form>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = true;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_EC_ELLIPTIC_CURVE_HPP_
