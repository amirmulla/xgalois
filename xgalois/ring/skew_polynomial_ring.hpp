#ifndef XGALOIS_RING_SKEW_POLYNOMIAL_RING_HPP_
#define XGALOIS_RING_SKEW_POLYNOMIAL_RING_HPP_

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"

namespace xg {
namespace ring {

namespace detail {

struct NoInvAuto {};

template <typename Fn, typename F, typename E>
inline E InvokeSkewMap(const Fn& fn, const F& field, const E& a) {
  if constexpr (std::invocable<const Fn&, const F&, const E&>) {
    return fn(field, a);
  } else {
    return fn(a);
  }
}

template <Field F>
struct ZeroDerivation {
  using E = typename F::element_type;
  E operator()(const F& field, const E&) const {
    return field.AdditiveIdentity();
  }
};

template <Field F>
struct FrobeniusStepMap {
  using E = typename F::element_type;

  uint32_t steps = 0;
  uint64_t fallback_exp = 1;

  E operator()(const F& field, const E& a) const {
    if (steps == 0) {
      return a;
    }
    if constexpr (requires(const F& f, const E& x) {
                    { f.Frobenius(x) } -> std::same_as<E>;
                  }) {
      E cur = a;
      for (uint32_t i = 0; i < steps; ++i) {
        cur = field.Frobenius(cur);
      }
      return cur;
    } else if constexpr (requires(const F& f, const E& x, uint64_t e) {
                           { f.Pow(x, e) } -> std::same_as<E>;
                         }) {
      return field.Pow(a, fallback_exp);
    } else {
      E acc = field.MultiplicativeIdentity();
      E base = a;
      uint64_t exp = fallback_exp;
      while (exp > 0) {
        if ((exp & uint64_t{1}) != 0) {
          acc = field.Mul(acc, base);
        }
        exp >>= 1;
        if (exp > 0) {
          base = field.Mul(base, base);
        }
      }
      return acc;
    }
  }
};

inline uint64_t SafePowU64(uint64_t base, uint32_t exp) noexcept {
  uint64_t acc = 1;
  uint64_t b = base;
  while (exp > 0) {
    if ((exp & 1U) != 0) {
      acc *= b;
    }
    exp >>= 1;
    if (exp > 0) {
      b *= b;
    }
  }
  return acc;
}

}  // namespace detail

// Ore / Skew polynomial ring F[x; sigma, delta] over a field F with
// commutation law: x * a = sigma(a) * x + delta(a).
// Elements are ascending-degree coefficient vectors with trailing zeros
// trimmed. Satisfies xg::Ring (non-commutative).
template <Field F, typename AutoFn, typename DerivFn,
          typename InvAutoFn = detail::NoInvAuto>
class SkewPolynomialRing {
 public:
  using field_type = F;
  using coeff_type = typename F::element_type;
  using element_type = std::vector<coeff_type>;

  static constexpr bool kHasInverseAuto =
      !std::same_as<std::remove_cvref_t<InvAutoFn>, detail::NoInvAuto>;

  SkewPolynomialRing(F field, AutoFn sigma, DerivFn delta)
    requires(!kHasInverseAuto)
      : field_(std::move(field)),
        sigma_(std::move(sigma)),
        delta_(std::move(delta)) {}

  SkewPolynomialRing(F field, AutoFn sigma, DerivFn delta, InvAutoFn sigma_inv)
    requires kHasInverseAuto
      : field_(std::move(field)),
        sigma_(std::move(sigma)),
        delta_(std::move(delta)),
        sigma_inv_(std::move(sigma_inv)) {}

  const F& BaseField() const noexcept { return field_; }

  coeff_type ApplyAuto(const coeff_type& a) const {
    return detail::InvokeSkewMap(sigma_, field_, a);
  }

  coeff_type ApplyDeriv(const coeff_type& a) const {
    return detail::InvokeSkewMap(delta_, field_, a);
  }

  coeff_type ApplyInvAuto(const coeff_type& a) const
    requires kHasInverseAuto
  {
    return detail::InvokeSkewMap(sigma_inv_, field_, a);
  }

  coeff_type ApplyAutoPower(const coeff_type& a, std::size_t k) const {
    coeff_type cur = a;
    for (std::size_t i = 0; i < k; ++i) {
      cur = ApplyAuto(cur);
    }
    return cur;
  }

  coeff_type ApplyInvAutoPower(const coeff_type& a, std::size_t k) const
    requires kHasInverseAuto
  {
    coeff_type cur = a;
    for (std::size_t i = 0; i < k; ++i) {
      cur = ApplyInvAuto(cur);
    }
    return cur;
  }

  element_type Normalize(element_type f) const {
    Trim(f);
    return f;
  }

  element_type Constant(const coeff_type& c) const {
    if (field_.IsZero(c)) {
      return {};
    }
    return element_type{c};
  }

  element_type Monomial(const coeff_type& c, std::size_t deg) const {
    if (field_.IsZero(c)) {
      return {};
    }
    element_type out(deg + 1, field_.AdditiveIdentity());
    out[deg] = c;
    return out;
  }

  element_type X() const {
    return Monomial(field_.MultiplicativeIdentity(), 1);
  }

  element_type AdditiveIdentity() const { return {}; }

  element_type MultiplicativeIdentity() const {
    return Constant(field_.MultiplicativeIdentity());
  }

  int64_t Degree(const element_type& f) const {
    for (std::size_t i = f.size(); i-- > 0;) {
      if (!field_.IsZero(f[i])) {
        return static_cast<int64_t>(i);
      }
    }
    return -1;
  }

  coeff_type LeadingCoeff(const element_type& f) const {
    const int64_t deg = Degree(f);
    if (deg < 0) {
      return field_.AdditiveIdentity();
    }
    return f[static_cast<std::size_t>(deg)];
  }

  bool IsZero(const element_type& f) const { return Degree(f) < 0; }

  bool IsOne(const element_type& f) const {
    return Degree(f) == 0 && field_.IsOne(f[0]);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    const int64_t da = Degree(a);
    const int64_t db = Degree(b);
    if (da != db) {
      return false;
    }
    for (int64_t i = 0; i <= da; ++i) {
      const auto idx = static_cast<std::size_t>(i);
      if (!alg::Equals(field_, a[idx], b[idx])) {
        return false;
      }
    }
    return true;
  }

  element_type Add(const element_type& a, const element_type& b) const {
    const std::size_t max_len = std::max(a.size(), b.size());
    const coeff_type zero = field_.AdditiveIdentity();
    element_type out;
    out.reserve(max_len);
    for (std::size_t i = 0; i < max_len; ++i) {
      const coeff_type& ca = (i < a.size()) ? a[i] : zero;
      const coeff_type& cb = (i < b.size()) ? b[i] : zero;
      out.push_back(field_.Add(ca, cb));
    }
    Trim(out);
    return out;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    const std::size_t max_len = std::max(a.size(), b.size());
    const coeff_type zero = field_.AdditiveIdentity();
    element_type out;
    out.reserve(max_len);
    for (std::size_t i = 0; i < max_len; ++i) {
      const coeff_type& ca = (i < a.size()) ? a[i] : zero;
      const coeff_type& cb = (i < b.size()) ? b[i] : zero;
      out.push_back(field_.Sub(ca, cb));
    }
    Trim(out);
    return out;
  }

  element_type Neg(const element_type& a) const {
    element_type out;
    out.reserve(a.size());
    for (const auto& c : a) {
      out.push_back(field_.Neg(c));
    }
    Trim(out);
    return out;
  }

  // Left-multiplication by a constant scalar: (c * f)_i = c * f_i.
  element_type LeftScalarMul(const coeff_type& c, const element_type& f) const {
    if (field_.IsZero(c) || IsZero(f)) {
      return {};
    }
    element_type out;
    out.reserve(f.size());
    for (const auto& coeff : f) {
      out.push_back(field_.Mul(c, coeff));
    }
    Trim(out);
    return out;
  }

  // Non-commutative skew polynomial multiplication f * g using
  // x * a = sigma(a) * x + delta(a).
  element_type Mul(const element_type& f, const element_type& g) const {
    const int64_t df = Degree(f);
    const int64_t dg = Degree(g);
    if (df < 0 || dg < 0) {
      return {};
    }
    const auto nf = static_cast<std::size_t>(df + 1);
    const auto ng = static_cast<std::size_t>(dg + 1);
    element_type out(nf + ng - 1, field_.AdditiveIdentity());

    // Maintain h_k = x^k * g iteratively via h_{k+1} = x * h_k.
    element_type h(g.begin(), g.begin() + static_cast<std::ptrdiff_t>(ng));
    for (std::size_t k = 0; k < nf; ++k) {
      if (!field_.IsZero(f[k])) {
        const coeff_type& fk = f[k];
        if (out.size() < h.size()) {
          out.resize(h.size(), field_.AdditiveIdentity());
        }
        for (std::size_t j = 0; j < h.size(); ++j) {
          if (!field_.IsZero(h[j])) {
            out[j] = field_.Add(out[j], field_.Mul(fk, h[j]));
          }
        }
      }
      if (k + 1 < nf) {
        h = MulByXOnLeft(h);
      }
    }
    Trim(out);
    return out;
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = Normalize(a);
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

  // Right Euclidean division: returns {q_R, r_R} such that
  // a == q_R * b + r_R and deg(r_R) < deg(b).
  std::pair<element_type, element_type> RightDivRem(
      const element_type& a, const element_type& b) const {
    element_type bn = Normalize(b);
    const int64_t db = Degree(bn);
    if (db < 0) {
      throw std::domain_error(
          "SkewPolynomialRing::RightDivRem: division by zero");
    }

    element_type rem = Normalize(a);
    const int64_t da = Degree(rem);
    if (da < db) {
      return {AdditiveIdentity(), std::move(rem)};
    }

    const auto udb = static_cast<std::size_t>(db);
    const auto q_len = static_cast<std::size_t>(da - db + 1);

    // Precompute x^k * b for k = 0 .. da - db.
    std::vector<element_type> x_pow_b(q_len);
    x_pow_b[0] = bn;
    for (std::size_t k = 1; k < q_len; ++k) {
      x_pow_b[k] = MulByXOnLeft(x_pow_b[k - 1]);
    }

    element_type q(q_len, field_.AdditiveIdentity());
    for (int64_t shift = da - db; shift >= 0; --shift) {
      const auto k = static_cast<std::size_t>(shift);
      if (k + udb >= rem.size() || field_.IsZero(rem[k + udb])) {
        continue;
      }
      const element_type& shifted = x_pow_b[k];
      const coeff_type lc_shifted = shifted[k + udb];
      const coeff_type c = field_.Mul(rem[k + udb], field_.Inv(lc_shifted));
      q[k] = c;
      for (std::size_t j = 0; j < shifted.size(); ++j) {
        if (!field_.IsZero(shifted[j])) {
          rem[j] = field_.Sub(rem[j], field_.Mul(c, shifted[j]));
        }
      }
    }

    Trim(q);
    Trim(rem);
    return {std::move(q), std::move(rem)};
  }

  // Left Euclidean division: returns {q_L, r_L} such that
  // a == b * q_L + r_L and deg(r_L) < deg(b). Requires sigma^{-1}.
  std::pair<element_type, element_type> LeftDivRem(const element_type& a,
                                                   const element_type& b) const
    requires kHasInverseAuto
  {
    element_type bn = Normalize(b);
    const int64_t db = Degree(bn);
    if (db < 0) {
      throw std::domain_error(
          "SkewPolynomialRing::LeftDivRem: division by zero");
    }

    element_type rem = Normalize(a);
    const int64_t da = Degree(rem);
    if (da < db) {
      return {AdditiveIdentity(), std::move(rem)};
    }

    const auto udb = static_cast<std::size_t>(db);
    const coeff_type inv_lc_b = field_.Inv(bn[udb]);
    const auto q_len = static_cast<std::size_t>(da - db + 1);
    element_type q(q_len, field_.AdditiveIdentity());

    for (int64_t shift = da - db; shift >= 0; --shift) {
      const auto k = static_cast<std::size_t>(shift);
      if (k + udb >= rem.size() || field_.IsZero(rem[k + udb])) {
        continue;
      }
      // Leading term of b * (c * x^k) is lc(b) * sigma^{db}(c) * x^{k + db}.
      const coeff_type c =
          ApplyInvAutoPower(field_.Mul(inv_lc_b, rem[k + udb]), udb);
      q[k] = c;
      // b * (c * x^k) == (b * c) * x^k by associativity.
      const element_type bc = Mul(bn, Constant(c));
      for (std::size_t j = 0; j < bc.size(); ++j) {
        if (!field_.IsZero(bc[j])) {
          rem[k + j] = field_.Sub(rem[k + j], bc[j]);
        }
      }
    }

    Trim(q);
    Trim(rem);
    return {std::move(q), std::move(rem)};
  }

  // Right Extended Euclidean Algorithm (GCRD):
  // Returns {g, u, v} such that u * a + v * b == g with g monic (when
  // non-zero). When stop_deg >= 0, terminates early as soon as deg(r) <=
  // stop_deg (partial Skew-EEA).
  std::tuple<element_type, element_type, element_type> RightExtendedGcd(
      const element_type& a, const element_type& b,
      int64_t stop_deg = -1) const {
    element_type r0 = Normalize(a);
    element_type r1 = Normalize(b);
    element_type u0 = MultiplicativeIdentity();
    element_type u1 = AdditiveIdentity();
    element_type v0 = AdditiveIdentity();
    element_type v1 = MultiplicativeIdentity();

    while (!IsZero(r1) && (stop_deg < 0 || Degree(r0) > stop_deg)) {
      auto [q, rem] = RightDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(rem);

      element_type next_u = Sub(u0, Mul(q, u1));
      u0 = std::move(u1);
      u1 = std::move(next_u);

      element_type next_v = Sub(v0, Mul(q, v1));
      v0 = std::move(v1);
      v1 = std::move(next_v);
    }

    if (!IsZero(r0)) {
      const coeff_type inv_lc = field_.Inv(LeadingCoeff(r0));
      r0 = LeftScalarMul(inv_lc, r0);
      u0 = LeftScalarMul(inv_lc, u0);
      v0 = LeftScalarMul(inv_lc, v0);
    }
    return {std::move(r0), std::move(u0), std::move(v0)};
  }

  // Left Extended Euclidean Algorithm (GCLD):
  // Returns {g, u, v} such that a * u + b * v == g with g monic (when
  // non-zero). When stop_deg >= 0, terminates early as soon as deg(r) <=
  // stop_deg.
  std::tuple<element_type, element_type, element_type> LeftExtendedGcd(
      const element_type& a, const element_type& b, int64_t stop_deg = -1) const
    requires kHasInverseAuto
  {
    element_type r0 = Normalize(a);
    element_type r1 = Normalize(b);
    element_type u0 = MultiplicativeIdentity();
    element_type u1 = AdditiveIdentity();
    element_type v0 = AdditiveIdentity();
    element_type v1 = MultiplicativeIdentity();

    while (!IsZero(r1) && (stop_deg < 0 || Degree(r0) > stop_deg)) {
      auto [q, rem] = LeftDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(rem);

      element_type next_u = Sub(u0, Mul(u1, q));
      u0 = std::move(u1);
      u1 = std::move(next_u);

      element_type next_v = Sub(v0, Mul(v1, q));
      v0 = std::move(v1);
      v1 = std::move(next_v);
    }

    if (!IsZero(r0)) {
      // Right-scale by scalar c = sigma^{-deg(r0)}(lc(r0)^{-1}) so r0 * c is
      // monic and a * (u0 * c) + b * (v0 * c) == r0 * c.
      const auto deg_r0 = static_cast<std::size_t>(Degree(r0));
      const coeff_type c =
          ApplyInvAutoPower(field_.Inv(LeadingCoeff(r0)), deg_r0);
      const element_type c_poly = Constant(c);
      r0 = Mul(r0, c_poly);
      u0 = Mul(u0, c_poly);
      v0 = Mul(v0, c_poly);
    }
    return {std::move(r0), std::move(u0), std::move(v0)};
  }

  // Least Common Left Multiple (monic): unique monic skew polynomial m of
  // minimal degree right-divisible by both a and b.
  element_type Lclm(const element_type& a, const element_type& b) const {
    element_type r0 = Normalize(a);
    element_type r1 = Normalize(b);
    if (IsZero(r0) || IsZero(r1)) {
      return AdditiveIdentity();
    }
    element_type u0 = MultiplicativeIdentity();
    element_type u1 = AdditiveIdentity();

    while (!IsZero(r1)) {
      auto [q, rem] = RightDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(rem);

      element_type next_u = Sub(u0, Mul(q, u1));
      u0 = std::move(u1);
      u1 = std::move(next_u);
    }

    // At termination r1 == 0 => u1 * a + v1 * b == 0, so m = u1 * a is a
    // minimal-degree common left multiple of a and b.
    element_type m = Mul(u1, a);
    if (!IsZero(m)) {
      m = LeftScalarMul(field_.Inv(LeadingCoeff(m)), m);
    }
    return m;
  }

  // Least Common Right Multiple (monic): unique monic skew polynomial m of
  // minimal degree left-divisible by both a and b.
  element_type Lcrm(const element_type& a, const element_type& b) const
    requires kHasInverseAuto
  {
    element_type r0 = Normalize(a);
    element_type r1 = Normalize(b);
    if (IsZero(r0) || IsZero(r1)) {
      return AdditiveIdentity();
    }
    element_type u0 = MultiplicativeIdentity();
    element_type u1 = AdditiveIdentity();

    while (!IsZero(r1)) {
      auto [q, rem] = LeftDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(rem);

      element_type next_u = Sub(u0, Mul(u1, q));
      u0 = std::move(u1);
      u1 = std::move(next_u);
    }

    // At termination r1 == 0 => a * u1 + b * v1 == 0, so m = a * u1 is a
    // minimal-degree common right multiple of a and b.
    element_type m = Mul(a, u1);
    if (!IsZero(m)) {
      const auto deg_m = static_cast<std::size_t>(Degree(m));
      const coeff_type c =
          ApplyInvAutoPower(field_.Inv(LeadingCoeff(m)), deg_m);
      m = Mul(m, Constant(c));
    }
    return m;
  }

  // Operator evaluation (for delta == 0): f(x) = sum_i f_i * sigma^i(x).
  coeff_type EvalOperator(const element_type& f, const coeff_type& x) const {
    coeff_type acc = field_.AdditiveIdentity();
    if (IsZero(f) || field_.IsZero(x)) {
      return acc;
    }
    coeff_type cur = x;
    for (std::size_t i = 0; i < f.size(); ++i) {
      if (!field_.IsZero(f[i])) {
        acc = field_.Add(acc, field_.Mul(f[i], cur));
      }
      if (i + 1 < f.size()) {
        cur = ApplyAuto(cur);
      }
    }
    return acc;
  }

  // O(k^2) minimal monic skew polynomial vanishing on span_{F^sigma}(points)
  // (valid when delta == 0).
  element_type MinimalPolynomial(const std::vector<coeff_type>& points) const {
    element_type m = MultiplicativeIdentity();
    for (const auto& pt : points) {
      const coeff_type beta = EvalOperator(m, pt);
      if (field_.IsZero(beta)) {
        continue;
      }
      // Left-multiply by (x - sigma(beta) / beta).
      const coeff_type c = field_.Mul(ApplyAuto(beta), field_.Inv(beta));
      element_type lin{field_.Neg(c), field_.MultiplicativeIdentity()};
      m = Mul(lin, m);
    }
    return m;
  }

  // O(k^2) Newton-style operator interpolation (valid when delta == 0):
  // returns skew polynomial f with deg(f) < k and
  // EvalOperator(f, points[i]) == values[i] for all i.
  element_type InterpolateOperator(
      const std::vector<coeff_type>& points,
      const std::vector<coeff_type>& values) const {
    if (points.size() != values.size()) {
      throw std::invalid_argument(
          "SkewPolynomialRing::InterpolateOperator: points and values size "
          "mismatch");
    }
    element_type f = AdditiveIdentity();
    element_type m = MultiplicativeIdentity();
    const std::size_t k = points.size();
    for (std::size_t i = 0; i < k; ++i) {
      const coeff_type beta = EvalOperator(m, points[i]);
      const coeff_type cur_val = EvalOperator(f, points[i]);
      const coeff_type diff = field_.Sub(values[i], cur_val);
      if (field_.IsZero(beta)) {
        if (!field_.IsZero(diff)) {
          throw std::invalid_argument(
              "SkewPolynomialRing::InterpolateOperator: inconsistent values on "
              "linearly dependent points");
        }
        continue;
      }
      const coeff_type lambda = field_.Mul(diff, field_.Inv(beta));
      if (!field_.IsZero(lambda)) {
        f = Add(f, LeftScalarMul(lambda, m));
      }
      if (i + 1 < k) {
        const coeff_type c = field_.Mul(ApplyAuto(beta), field_.Inv(beta));
        element_type lin{field_.Neg(c), field_.MultiplicativeIdentity()};
        m = Mul(lin, m);
      }
    }
    return f;
  }

  std::string ToString(const element_type& f) const {
    const int64_t deg = Degree(f);
    if (deg < 0) {
      return "0";
    }
    std::string out;
    bool first = true;
    for (int64_t i = deg; i >= 0; --i) {
      const auto idx = static_cast<std::size_t>(i);
      if (field_.IsZero(f[idx])) {
        continue;
      }
      if (!first) {
        out += " + ";
      }
      first = false;
      const std::string c_str = field_.ToString(f[idx]);
      if (i == 0) {
        out += c_str;
      } else if (i == 1) {
        out += c_str + "*x";
      } else {
        out += c_str + "*x^" + std::to_string(i);
      }
    }
    return out;
  }

 private:
  // Computes x * h using x * c = sigma(c) * x + delta(c).
  element_type MulByXOnLeft(const element_type& h) const {
    if (h.empty()) {
      return {};
    }
    element_type out(h.size() + 1, field_.AdditiveIdentity());
    for (std::size_t i = 0; i < h.size(); ++i) {
      if (field_.IsZero(h[i])) {
        continue;
      }
      out[i + 1] = field_.Add(out[i + 1], ApplyAuto(h[i]));
      if constexpr (!std::same_as<std::remove_cvref_t<DerivFn>,
                                  detail::ZeroDerivation<F>>) {
        out[i] = field_.Add(out[i], ApplyDeriv(h[i]));
      }
    }
    Trim(out);
    return out;
  }

  void Trim(element_type& f) const {
    while (!f.empty() && field_.IsZero(f.back())) {
      f.pop_back();
    }
  }

  F field_;
  AutoFn sigma_;
  DerivFn delta_;
  [[no_unique_address]] InvAutoFn sigma_inv_{};
};

// Specialization of SkewPolynomialRing with delta = 0 and Frobenius
// automorphism sigma(a) = a^{p^s} (with precomputed inverse sigma^{-1} =
// a^{p^{m - s}}).
template <Field F>
class FrobeniusSkewRing
    : public SkewPolynomialRing<F, detail::FrobeniusStepMap<F>,
                                detail::ZeroDerivation<F>,
                                detail::FrobeniusStepMap<F>> {
  using Base = SkewPolynomialRing<F, detail::FrobeniusStepMap<F>,
                                  detail::ZeroDerivation<F>,
                                  detail::FrobeniusStepMap<F>>;

 public:
  using field_type = typename Base::field_type;
  using coeff_type = typename Base::coeff_type;
  using element_type = typename Base::element_type;

  // Constructs F[x; sigma] where sigma(a) = a^{p^{frobenius_step}} (default
  // frobenius_step = 1 gives the base Frobenius a -> a^p).
  explicit FrobeniusSkewRing(F field, uint32_t frobenius_step = 1)
      : Base(field, MakeForwardMap(field, frobenius_step),
             detail::ZeroDerivation<F>{},
             MakeInverseMap(field, frobenius_step)),
        frobenius_step_(ComputeReducedStep(field, frobenius_step)) {}

  // Constructs F[x; sigma] where sigma(a) = a^q for subfield order q = p^s.
  static FrobeniusSkewRing FromSubfieldOrder(F field, uint64_t q) {
    uint64_t p = 0;
    if constexpr (requires(const F& f) { f.Characteristic(); }) {
      p = static_cast<uint64_t>(field.Characteristic());
    }
    if (p < 2 || q < 2) {
      throw std::invalid_argument(
          "FrobeniusSkewRing::FromSubfieldOrder: invalid characteristic or q");
    }
    uint32_t s = 0;
    uint64_t cur = 1;
    while (cur < q) {
      if (cur > q / p) {
        throw std::invalid_argument(
            "FrobeniusSkewRing::FromSubfieldOrder: q must be a power of "
            "characteristic p");
      }
      cur *= p;
      ++s;
    }
    if (cur != q) {
      throw std::invalid_argument(
          "FrobeniusSkewRing::FromSubfieldOrder: q must be a power of "
          "characteristic p");
    }
    return FrobeniusSkewRing(std::move(field), s);
  }

  uint32_t FrobeniusStep() const noexcept { return frobenius_step_; }

 private:
  static uint32_t FieldDegree(const F& field) {
    if constexpr (requires(const F& f) { f.Degree(); }) {
      return static_cast<uint32_t>(field.Degree());
    } else {
      return 1;
    }
  }

  static uint64_t FieldCharacteristic(const F& field) {
    if constexpr (requires(const F& f) { f.Characteristic(); }) {
      return static_cast<uint64_t>(field.Characteristic());
    } else {
      return 1;
    }
  }

  static uint32_t ComputeReducedStep(const F& field, uint32_t s) {
    const uint32_t m = FieldDegree(field);
    return (m == 0) ? 0 : (s % m);
  }

  static detail::FrobeniusStepMap<F> MakeForwardMap(const F& field,
                                                    uint32_t s) {
    const uint32_t step = ComputeReducedStep(field, s);
    const uint64_t p = FieldCharacteristic(field);
    return detail::FrobeniusStepMap<F>{step, detail::SafePowU64(p, step)};
  }

  static detail::FrobeniusStepMap<F> MakeInverseMap(const F& field,
                                                    uint32_t s) {
    const uint32_t m = FieldDegree(field);
    const uint32_t step = ComputeReducedStep(field, s);
    const uint32_t inv_step = (m == 0 || step == 0) ? 0 : (m - step);
    const uint64_t p = FieldCharacteristic(field);
    return detail::FrobeniusStepMap<F>{inv_step,
                                       detail::SafePowU64(p, inv_step)};
  }

  uint32_t frobenius_step_;
};

}  // namespace ring

template <Field F, typename AutoFn, typename DerivFn, typename InvAutoFn>
struct AlgebraTraits<ring::SkewPolynomialRing<F, AutoFn, DerivFn, InvAutoFn>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <Field F>
struct AlgebraTraits<ring::FrobeniusSkewRing<F>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_SKEW_POLYNOMIAL_RING_HPP_
