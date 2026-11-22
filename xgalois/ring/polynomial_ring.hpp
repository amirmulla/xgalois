#ifndef XGALOIS_RING_POLYNOMIAL_RING_HPP_
#define XGALOIS_RING_POLYNOMIAL_RING_HPP_

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"

namespace xg {
namespace ring {

// Polynomial element over an arbitrary ring R, stored in ascending degree order
// (coeffs[i] is the coefficient of x^i) with canonical trailing zeros trimmed.
template <Ring R>
struct RingPoly {
  using coeff_type = typename R::element_type;

  std::vector<coeff_type> coeffs;

  RingPoly() = default;
  explicit RingPoly(std::vector<coeff_type> c) : coeffs(std::move(c)) {}
  RingPoly(std::initializer_list<coeff_type> init) : coeffs(init) {}

  bool Empty() const noexcept { return coeffs.empty(); }
  std::size_t Size() const noexcept { return coeffs.size(); }

  const coeff_type& operator[](std::size_t i) const { return coeffs[i]; }
  coeff_type& operator[](std::size_t i) { return coeffs[i]; }

  bool operator==(const RingPoly& other) const
    requires requires(const coeff_type& a, const coeff_type& b) {
      { a == b } -> std::convertible_to<bool>;
    }
  = default;
};

// Univariate polynomial ring R[x] over any xg::Ring R.
// Satisfies CommutativeRing when CommutativeRing<R>, IntegralDomain when
// IntegralDomain<R>, and EuclideanDomain when Field<R>.
template <Ring R>
class PolynomialRing {
 public:
  using base_ring_type = R;
  using coeff_type = typename R::element_type;
  using scalar_type = coeff_type;
  using element_type = RingPoly<R>;

  PolynomialRing()
    requires std::default_initializable<R>
  = default;

  explicit PolynomialRing(R base_ring) : ring_(std::move(base_ring)) {}

  const R& BaseRing() const noexcept { return ring_; }

  element_type Normalize(element_type p) const {
    Trim(p);
    return p;
  }

  element_type FromCoeffs(std::vector<coeff_type> coeffs) const {
    element_type p(std::move(coeffs));
    Trim(p);
    return p;
  }

  element_type Constant(const coeff_type& c) const {
    if (ring_.IsZero(c)) {
      return {};
    }
    return element_type(std::vector<coeff_type>{c});
  }

  element_type Monomial(const coeff_type& c, std::size_t degree) const {
    if (ring_.IsZero(c)) {
      return {};
    }
    std::vector<coeff_type> coeffs(degree + 1, ring_.AdditiveIdentity());
    coeffs[degree] = c;
    return element_type(std::move(coeffs));
  }

  element_type X() const { return Monomial(ring_.MultiplicativeIdentity(), 1); }

  element_type AdditiveIdentity() const { return {}; }

  element_type MultiplicativeIdentity() const {
    return Constant(ring_.MultiplicativeIdentity());
  }

  // Returns -1 for the zero polynomial, otherwise the highest non-zero degree.
  int64_t Degree(const element_type& p) const {
    for (std::size_t i = p.coeffs.size(); i-- > 0;) {
      if (!ring_.IsZero(p.coeffs[i])) {
        return static_cast<int64_t>(i);
      }
    }
    return -1;
  }

  coeff_type LeadingCoeff(const element_type& p) const {
    const int64_t deg = Degree(p);
    if (deg < 0) {
      return ring_.AdditiveIdentity();
    }
    return p.coeffs[static_cast<std::size_t>(deg)];
  }

  bool IsZero(const element_type& p) const { return Degree(p) < 0; }

  bool IsOne(const element_type& p) const {
    if (ring_.IsZero(ring_.MultiplicativeIdentity())) {
      return IsZero(p);
    }
    return Degree(p) == 0 && ring_.IsOne(p.coeffs[0]);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    const int64_t da = Degree(a);
    const int64_t db = Degree(b);
    if (da != db) {
      return false;
    }
    for (int64_t i = 0; i <= da; ++i) {
      const auto idx = static_cast<std::size_t>(i);
      if (!alg::Equals(ring_, a.coeffs[idx], b.coeffs[idx])) {
        return false;
      }
    }
    return true;
  }

  element_type Add(const element_type& a, const element_type& b) const {
    const std::size_t max_len = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<coeff_type> out;
    out.reserve(max_len);
    const coeff_type zero = ring_.AdditiveIdentity();
    for (std::size_t i = 0; i < max_len; ++i) {
      const coeff_type& ca = (i < a.coeffs.size()) ? a.coeffs[i] : zero;
      const coeff_type& cb = (i < b.coeffs.size()) ? b.coeffs[i] : zero;
      out.push_back(ring_.Add(ca, cb));
    }
    element_type res(std::move(out));
    Trim(res);
    return res;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    const std::size_t max_len = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<coeff_type> out;
    out.reserve(max_len);
    const coeff_type zero = ring_.AdditiveIdentity();
    for (std::size_t i = 0; i < max_len; ++i) {
      const coeff_type& ca = (i < a.coeffs.size()) ? a.coeffs[i] : zero;
      const coeff_type& cb = (i < b.coeffs.size()) ? b.coeffs[i] : zero;
      out.push_back(ring_.Sub(ca, cb));
    }
    element_type res(std::move(out));
    Trim(res);
    return res;
  }

  element_type Neg(const element_type& a) const {
    std::vector<coeff_type> out;
    out.reserve(a.coeffs.size());
    for (const auto& c : a.coeffs) {
      out.push_back(ring_.Neg(c));
    }
    element_type res(std::move(out));
    Trim(res);
    return res;
  }

  element_type ScalarMul(const coeff_type& scalar,
                         const element_type& p) const {
    if (ring_.IsZero(scalar) || IsZero(p)) {
      return {};
    }
    std::vector<coeff_type> out;
    out.reserve(p.coeffs.size());
    for (const auto& c : p.coeffs) {
      out.push_back(ring_.Mul(scalar, c));
    }
    element_type res(std::move(out));
    Trim(res);
    return res;
  }

  element_type Scale(const coeff_type& scalar, const element_type& p) const {
    return ScalarMul(scalar, p);
  }

  // Multiplies a and b preserving coefficient multiplication order a[i] * b[j]
  // and trimming trailing zeros (including zero-divisor cancellations).
  element_type Mul(const element_type& a, const element_type& b) const {
    const int64_t da = Degree(a);
    const int64_t db = Degree(b);
    if (da < 0 || db < 0) {
      return {};
    }
    const auto na = static_cast<std::size_t>(da + 1);
    const auto nb = static_cast<std::size_t>(db + 1);
    std::vector<coeff_type> out(na + nb - 1, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < na; ++i) {
      if (ring_.IsZero(a.coeffs[i])) {
        continue;
      }
      for (std::size_t j = 0; j < nb; ++j) {
        if (ring_.IsZero(b.coeffs[j])) {
          continue;
        }
        out[i + j] = ring_.Add(out[i + j], ring_.Mul(a.coeffs[i], b.coeffs[j]));
      }
    }
    element_type res(std::move(out));
    Trim(res);
    return res;
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

  // Evaluates p(x) using Horner's scheme.
  coeff_type Eval(const element_type& p, const coeff_type& x) const {
    const int64_t deg = Degree(p);
    if (deg < 0) {
      return ring_.AdditiveIdentity();
    }
    coeff_type acc = p.coeffs[static_cast<std::size_t>(deg)];
    for (int64_t i = deg - 1; i >= 0; --i) {
      acc = ring_.Add(ring_.Mul(acc, x), p.coeffs[static_cast<std::size_t>(i)]);
    }
    return acc;
  }

  // Formal derivative p'(x) = sum_{k=1}^{deg} (k * p_k) x^{k-1}.
  element_type Derivative(const element_type& p) const {
    const int64_t deg = Degree(p);
    if (deg <= 0) {
      return {};
    }
    const AdditiveGroupView<R> add_group(ring_);
    std::vector<coeff_type> out;
    out.reserve(static_cast<std::size_t>(deg));
    for (int64_t k = 1; k <= deg; ++k) {
      out.push_back(add_group.Pow(p.coeffs[static_cast<std::size_t>(k)],
                                  static_cast<uint64_t>(k)));
    }
    element_type res(std::move(out));
    Trim(res);
    return res;
  }

  // Exact Euclidean division by a monic polynomial over ANY CommutativeRing R.
  // Returns {q, r} with a == q * monic_b + r and deg(r) < deg(monic_b).
  std::pair<element_type, element_type> MonicDivRem(
      const element_type& a, const element_type& monic_b) const
    requires CommutativeRing<R>
  {
    element_type b = Normalize(monic_b);
    const int64_t db = Degree(b);
    if (db < 0) {
      throw std::domain_error("PolynomialRing::MonicDivRem: division by zero");
    }
    if (!ring_.IsOne(b.coeffs[static_cast<std::size_t>(db)])) {
      throw std::invalid_argument(
          "PolynomialRing::MonicDivRem: divisor must be monic");
    }

    element_type rem = Normalize(a);
    const int64_t da = Degree(rem);
    if (da < db) {
      return {AdditiveIdentity(), std::move(rem)};
    }

    const auto q_len = static_cast<std::size_t>(da - db + 1);
    const auto udb = static_cast<std::size_t>(db);
    element_type q(std::vector<coeff_type>(q_len, ring_.AdditiveIdentity()));

    for (int64_t shift = da - db; shift >= 0; --shift) {
      const auto ushift = static_cast<std::size_t>(shift);
      const coeff_type c = rem.coeffs[ushift + udb];
      if (ring_.IsZero(c)) {
        continue;
      }
      q.coeffs[ushift] = c;
      for (std::size_t j = 0; j <= udb; ++j) {
        rem.coeffs[ushift + j] =
            ring_.Sub(rem.coeffs[ushift + j], ring_.Mul(c, b.coeffs[j]));
      }
    }

    Trim(q);
    Trim(rem);
    return {std::move(q), std::move(rem)};
  }

  // Division-free pseudo-division over any CommutativeRing / IntegralDomain R.
  // Returns {q, r} such that lc(b)^(max(0, deg(a) - deg(b) + 1)) * a
  // == q * b + r with deg(r) < deg(b).
  std::pair<element_type, element_type> PseudoDivRem(
      const element_type& a, const element_type& b) const
    requires CommutativeRing<R>
  {
    element_type bn = Normalize(b);
    const int64_t db = Degree(bn);
    if (db < 0) {
      throw std::domain_error("PolynomialRing::PseudoDivRem: division by zero");
    }

    element_type rem = Normalize(a);
    const int64_t da = Degree(rem);
    if (da < db) {
      return {AdditiveIdentity(), std::move(rem)};
    }

    const auto udb = static_cast<std::size_t>(db);
    const coeff_type lc = bn.coeffs[udb];
    const auto q_len = static_cast<std::size_t>(da - db + 1);
    element_type q(std::vector<coeff_type>(q_len, ring_.AdditiveIdentity()));
    rem.coeffs.resize(static_cast<std::size_t>(da + 1),
                      ring_.AdditiveIdentity());

    for (int64_t shift = da - db; shift >= 0; --shift) {
      const auto ushift = static_cast<std::size_t>(shift);
      const coeff_type c = rem.coeffs[ushift + udb];

      for (std::size_t k = ushift + 1; k < q_len; ++k) {
        q.coeffs[k] = ring_.Mul(lc, q.coeffs[k]);
      }
      q.coeffs[ushift] = c;

      for (std::size_t k = 0; k <= ushift + udb; ++k) {
        rem.coeffs[k] = ring_.Mul(lc, rem.coeffs[k]);
      }
      if (!ring_.IsZero(c)) {
        for (std::size_t j = 0; j <= udb; ++j) {
          rem.coeffs[ushift + j] =
              ring_.Sub(rem.coeffs[ushift + j], ring_.Mul(c, bn.coeffs[j]));
        }
      }
    }

    Trim(q);
    Trim(rem);
    return {std::move(q), std::move(rem)};
  }

  // Euclidean division over a Field R.
  std::pair<element_type, element_type> DivRem(const element_type& a,
                                               const element_type& b) const
    requires Field<R>
  {
    element_type bn = Normalize(b);
    const int64_t db = Degree(bn);
    if (db < 0) {
      throw std::domain_error("PolynomialRing::DivRem: division by zero");
    }

    element_type rem = Normalize(a);
    const int64_t da = Degree(rem);
    if (da < db) {
      return {AdditiveIdentity(), std::move(rem)};
    }

    const auto udb = static_cast<std::size_t>(db);
    const coeff_type inv_lc = ring_.Inv(bn.coeffs[udb]);
    const auto q_len = static_cast<std::size_t>(da - db + 1);
    element_type q(std::vector<coeff_type>(q_len, ring_.AdditiveIdentity()));

    for (int64_t shift = da - db; shift >= 0; --shift) {
      const auto ushift = static_cast<std::size_t>(shift);
      const coeff_type top = rem.coeffs[ushift + udb];
      if (ring_.IsZero(top)) {
        continue;
      }
      const coeff_type factor = ring_.Mul(top, inv_lc);
      q.coeffs[ushift] = factor;
      for (std::size_t j = 0; j <= udb; ++j) {
        rem.coeffs[ushift + j] =
            ring_.Sub(rem.coeffs[ushift + j], ring_.Mul(factor, bn.coeffs[j]));
      }
    }

    Trim(q);
    Trim(rem);
    return {std::move(q), std::move(rem)};
  }

  // Euclidean norm: 0 for zero polynomial, deg(p) + 1 for non-zero p.
  uint64_t EuclideanNorm(const element_type& p) const
    requires Field<R>
  {
    const int64_t deg = Degree(p);
    return (deg < 0) ? uint64_t{0} : static_cast<uint64_t>(deg + 1);
  }

  bool IsUnit(const element_type& p) const
    requires Field<R>
  {
    return Degree(p) == 0;
  }

  std::optional<element_type> TryInv(const element_type& p) const
    requires Field<R>
  {
    if (Degree(p) != 0) {
      return std::nullopt;
    }
    return Constant(ring_.Inv(p.coeffs[0]));
  }

  // Returns {monic_p, unit_const} with monic_p == Mul(unit_const, p).
  std::pair<element_type, element_type> AssociateCanonical(
      const element_type& p) const
    requires Field<R>
  {
    element_type pn = Normalize(p);
    if (IsZero(pn)) {
      return {AdditiveIdentity(), MultiplicativeIdentity()};
    }
    const coeff_type inv_lc = ring_.Inv(LeadingCoeff(pn));
    return {ScalarMul(inv_lc, pn), Constant(inv_lc)};
  }

  // Extended Euclidean Algorithm over Field R: returns {g, s, t} with g monic
  // (or zero if a == b == 0) and s * a + t * b == g.
  std::tuple<element_type, element_type, element_type> ExtendedGcd(
      const element_type& a, const element_type& b) const
    requires Field<R>
  {
    element_type old_r = Normalize(a);
    element_type r = Normalize(b);
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

    auto [g_monic, unit] = AssociateCanonical(old_r);
    return {std::move(g_monic), Mul(unit, old_s), Mul(unit, old_t)};
  }

  std::string ToString(const element_type& p) const {
    const int64_t deg = Degree(p);
    if (deg < 0) {
      return "0";
    }
    std::string out;
    bool first = true;
    for (int64_t i = deg; i >= 0; --i) {
      const auto idx = static_cast<std::size_t>(i);
      if (ring_.IsZero(p.coeffs[idx])) {
        continue;
      }
      if (!first) {
        out += " + ";
      }
      first = false;
      const std::string c_str = ring_.ToString(p.coeffs[idx]);
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
  void Trim(element_type& p) const {
    while (!p.coeffs.empty() && ring_.IsZero(p.coeffs.back())) {
      p.coeffs.pop_back();
    }
  }

  R ring_;
};

}  // namespace ring

template <Ring R>
struct AlgebraTraits<ring::PolynomialRing<R>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = detail::TraitIsCommutativeMul<R>();
  static constexpr bool kHasNoZeroDivisors =
      detail::TraitHasNoZeroDivisors<R>();
  static constexpr bool kIsPid = detail::TraitIsField<R>();
  static constexpr bool kIsEuclideanDomain = detail::TraitIsField<R>();
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_POLYNOMIAL_RING_HPP_
