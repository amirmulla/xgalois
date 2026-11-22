#ifndef XGALOIS_POLY_HALF_GCD_HPP_
#define XGALOIS_POLY_HALF_GCD_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace poly {

inline constexpr int64_t kHalfGcdCutoff = 32;
inline constexpr int64_t kNewtonDivCutoff = 64;

template <FiniteField F>
struct HalfGcdMatrix {
  using Poly = PolynomialDense<F>;

  Poly m00;
  Poly m01;
  Poly m10;
  Poly m11;

  static HalfGcdMatrix Identity(const std::shared_ptr<F>& field,
                                const std::string& var = "x") {
    return HalfGcdMatrix{Poly::One(field, var), Poly::Zero(field, var),
                         Poly::Zero(field, var), Poly::One(field, var)};
  }

  std::pair<Poly, Poly> Apply(const Poly& a, const Poly& b) const {
    return {m00 * a + m01 * b, m10 * a + m11 * b};
  }

  // Computes lhs * rhs (2x2 matrix multiplication over F[x]).
  static HalfGcdMatrix Multiply(const HalfGcdMatrix& lhs,
                                const HalfGcdMatrix& rhs) {
    return HalfGcdMatrix{lhs.m00 * rhs.m00 + lhs.m01 * rhs.m10,
                         lhs.m00 * rhs.m01 + lhs.m01 * rhs.m11,
                         lhs.m10 * rhs.m00 + lhs.m11 * rhs.m10,
                         lhs.m10 * rhs.m01 + lhs.m11 * rhs.m11};
  }

  // Left-multiplies *this in-place by the Euclidean step matrix
  //   Q = [[0, 1], [1, -q]]
  // corresponding to (r_prev, r_curr) -> (r_curr, r_prev - q * r_curr).
  void LeftMulStep(const Poly& q) {
    Poly next_m10 = m00 - q * m10;
    Poly next_m11 = m01 - q * m11;
    m00 = std::move(m10);
    m01 = std::move(m11);
    m10 = std::move(next_m10);
    m11 = std::move(next_m11);
  }
};

namespace detail {

template <FiniteField F>
inline void ValidateSamePolyField(const PolynomialDense<F>& a,
                                  const PolynomialDense<F>& b) {
  if (!a.HasField() || !b.HasField()) {
    throw std::invalid_argument(
        "HalfGcd: polynomials must have an associated field");
  }
  if (a.Field() != b.Field()) {
    throw std::invalid_argument(
        "HalfGcd: polynomials must belong to the same field");
  }
}

// Returns floor(p / x^shift) (drops terms of degree < shift).
template <FiniteField F>
inline PolynomialDense<F> TopShift(const PolynomialDense<F>& p, int64_t shift) {
  if (shift <= 0) {
    return p;
  }
  const std::size_t s = static_cast<std::size_t>(shift);
  if (p.IsZero() || s >= p.Size()) {
    return PolynomialDense<F>::Zero(p.Field(), p.GetVariable());
  }
  const auto& raw = p.RawCoeffs();
  std::vector<typename F::element_type> out(
      raw.begin() + static_cast<std::ptrdiff_t>(s), raw.end());
  return PolynomialDense<F>::FromRaw(std::move(out), p.Field(),
                                     p.GetVariable());
}

// Returns p mod x^len (keeps terms of degree < len).
template <FiniteField F>
inline PolynomialDense<F> TruncatePoly(const PolynomialDense<F>& p,
                                       std::size_t len) {
  if (len == 0 || p.IsZero()) {
    return PolynomialDense<F>::Zero(p.Field(), p.GetVariable());
  }
  const auto& raw = p.RawCoeffs();
  if (raw.size() <= len) {
    return p;
  }
  std::vector<typename F::element_type> out(
      raw.begin(), raw.begin() + static_cast<std::ptrdiff_t>(len));
  return PolynomialDense<F>::FromRaw(std::move(out), p.Field(),
                                     p.GetVariable());
}

// Computes g(x)^(-1) mod x^precision via Newton iteration assuming g(0) == 1.
// Valid across all field characteristics (including characteristic 2).
template <FiniteField F>
inline PolynomialDense<F> PowerSeriesInvModX(const PolynomialDense<F>& g,
                                             std::size_t precision) {
  using Poly = PolynomialDense<F>;
  using Raw = typename F::element_type;
  const auto field = g.Field();
  const std::string& var = g.GetVariable();
  if (precision == 0) {
    return Poly::Zero(field, var);
  }
  const Raw one = field->MultiplicativeIdentity();
  const Raw two = field->Add(one, one);

  Poly h = Poly::One(field, var);
  std::size_t cur = 1;
  while (cur < precision) {
    const std::size_t next = std::min(cur * 2, precision);
    const Poly g_trunc = TruncatePoly(g, next);
    const Poly gh = TruncatePoly(g_trunc * h, next);

    // Compute (2 - gh) mod x^next.
    const auto& gh_raw = gh.RawCoeffs();
    std::vector<Raw> two_minus_gh(std::max<std::size_t>(1, gh_raw.size()),
                                  field->AdditiveIdentity());
    const Raw gh0 = gh_raw.empty() ? field->AdditiveIdentity() : gh_raw[0];
    two_minus_gh[0] = field->Sub(two, gh0);
    for (std::size_t i = 1; i < gh_raw.size(); ++i) {
      two_minus_gh[i] = field->Neg(gh_raw[i]);
    }
    const Poly factor = Poly::FromRaw(std::move(two_minus_gh), field, var);
    h = TruncatePoly(h * factor, next);
    cur = next;
  }
  return h;
}

// Fast polynomial division with remainder using Newton power-series inversion
// above kNewtonDivCutoff, falling back to PolynomialDense::DivRem below.
template <FiniteField F>
inline std::pair<PolynomialDense<F>, PolynomialDense<F>> FastDivRem(
    const PolynomialDense<F>& a, const PolynomialDense<F>& b) {
  using Poly = PolynomialDense<F>;
  using Raw = typename F::element_type;
  if (b.IsZero()) {
    throw std::invalid_argument("FastDivRem: division by zero polynomial");
  }
  const int64_t deg_a = static_cast<int64_t>(a.Degree());
  const int64_t deg_b = static_cast<int64_t>(b.Degree());
  if (deg_a < deg_b) {
    return {Poly::Zero(a.Field(), a.GetVariable()), a};
  }
  const int64_t deg_q = deg_a - deg_b;
  if (deg_b < kNewtonDivCutoff || deg_q < kNewtonDivCutoff) {
    return a.DivRem(b);
  }

  const auto field = a.Field();
  const Raw lc_b = b.LeadingRawCoeff();
  const bool is_monic = field->IsOne(lc_b);
  const auto inv_lc =
      GaloisFieldElement<F>::FromInternal(field->Inv(lc_b), field);
  const Poly b_monic = is_monic ? b : (b * inv_lc);

  const std::size_t k = static_cast<std::size_t>(deg_q) + 1;
  const Poly rev_b = b_monic.Reverse(static_cast<std::size_t>(deg_b));
  const Poly inv_rev_b = PowerSeriesInvModX(rev_b, k);
  const Poly rev_a =
      TruncatePoly(a.Reverse(static_cast<std::size_t>(deg_a)), k);
  const Poly rev_q = TruncatePoly(rev_a * inv_rev_b, k);
  Poly q = rev_q.Reverse(k - 1);
  if (!is_monic) {
    q *= inv_lc;
  }
  Poly r = a - q * b;
  return {std::move(q), std::move(r)};
}

// Core recursive Half-GCD assuming a != 0 and deg(a) > deg(b).
// Returns unimodular matrix M such that (r_j, r_{j+1})^T = M * (a, b)^T is
// the unique consecutive EEA remainder pair satisfying:
//   deg(r_j) >= ceil(deg(a) / 2) > deg(r_{j+1}).
template <FiniteField F>
inline HalfGcdMatrix<F> HalfGcdRec(const PolynomialDense<F>& a,
                                   const PolynomialDense<F>& b) {
  using Poly = PolynomialDense<F>;
  const auto field = a.Field();
  const std::string& var = a.GetVariable();
  const int64_t d = static_cast<int64_t>(a.Degree());
  const int64_t m = (d + 1) / 2;

  if (b.IsZero() || static_cast<int64_t>(b.Degree()) < m) {
    return HalfGcdMatrix<F>::Identity(field, var);
  }

  if (d <= kHalfGcdCutoff) {
    Poly r0 = a;
    Poly r1 = b;
    HalfGcdMatrix<F> mat = HalfGcdMatrix<F>::Identity(field, var);
    while (!r1.IsZero() && static_cast<int64_t>(r1.Degree()) >= m) {
      auto [q, r_next] = r0.DivRem(r1);
      r0 = std::move(r1);
      r1 = std::move(r_next);
      mat.LeftMulStep(q);
    }
    return mat;
  }

  // First half-size recursive call on top (d - m) degrees.
  const Poly a1 = TopShift(a, m);
  const Poly b1 = TopShift(b, m);
  HalfGcdMatrix<F> r_mat = HalfGcdRec(a1, b1);

  auto [c, d_poly] = r_mat.Apply(a, b);
  if (d_poly.IsZero() || static_cast<int64_t>(d_poly.Degree()) < m) {
    return r_mat;
  }

  // Single bridging Euclidean division step.
  auto [q, e] = FastDivRem(c, d_poly);
  r_mat.LeftMulStep(q);
  if (e.IsZero() || static_cast<int64_t>(e.Degree()) < m) {
    return r_mat;
  }

  // Second half-size recursive call with exact shift s2 = 2*m - deg(d_poly)
  // so that ceil((deg(d_poly) - s2) / 2) + s2 == m.
  const int64_t k = static_cast<int64_t>(d_poly.Degree());
  const int64_t s2 = 2 * m - k;
  const Poly d2 = TopShift(d_poly, s2);
  const Poly e2 = TopShift(e, s2);
  HalfGcdMatrix<F> s_mat = HalfGcdRec(d2, e2);

  return HalfGcdMatrix<F>::Multiply(s_mat, r_mat);
}

}  // namespace detail

// Moller / von zur Gathen-Gerhard recursive 2x2 matrix Half-GCD over F[x].
// Returns unimodular matrix M mapping (a, b)^T to a consecutive EEA remainder
// pair (r_j, r_{j+1})^T whose degrees straddle ceil(deg(a) / 2) (when
// deg(a) > deg(b)).
template <FiniteField F>
inline HalfGcdMatrix<F> HalfGcd(const PolynomialDense<F>& a,
                                const PolynomialDense<F>& b) {
  using Poly = PolynomialDense<F>;
  detail::ValidateSamePolyField(a, b);
  const auto field = a.Field();
  const std::string& var = a.GetVariable();

  if (b.IsZero()) {
    return HalfGcdMatrix<F>::Identity(field, var);
  }
  if (a.IsZero()) {
    return HalfGcdMatrix<F>{Poly::Zero(field, var), Poly::One(field, var),
                            Poly::One(field, var), Poly::Zero(field, var)};
  }
  if (a.Degree() <= b.Degree()) {
    auto [q, r_next] = detail::FastDivRem(a, b);
    HalfGcdMatrix<F> step = HalfGcdMatrix<F>::Identity(field, var);
    step.LeftMulStep(q);
    if (r_next.IsZero()) {
      return step;
    }
    HalfGcdMatrix<F> h = detail::HalfGcdRec(b, r_next);
    return HalfGcdMatrix<F>::Multiply(h, step);
  }
  return detail::HalfGcdRec(a, b);
}

// Fast O(M(d) log d) Extended Euclidean Algorithm via Half-GCD.
// Returns {g, u, v} with monic g (or {0, 0, 0} if a == b == 0) such that
//   u * a + v * b == g.
template <FiniteField F>
inline std::tuple<PolynomialDense<F>, PolynomialDense<F>, PolynomialDense<F>>
FastExtendedGcd(const PolynomialDense<F>& a, const PolynomialDense<F>& b) {
  using Poly = PolynomialDense<F>;
  using Elem = typename Poly::ElementType;
  detail::ValidateSamePolyField(a, b);
  const auto field = a.Field();
  const std::string& var = a.GetVariable();
  const Poly zero = Poly::Zero(field, var);
  const Poly one = Poly::One(field, var);

  if (a.IsZero() && b.IsZero()) {
    return {zero, zero, zero};
  }
  if (b.IsZero()) {
    const Elem lc = a.LeadingCoeff();
    if (!lc.IsOne()) {
      const Elem inv_lc = lc.Inv();
      return {a * inv_lc, Poly::FromRaw({inv_lc.RawValue()}, field, var), zero};
    }
    return {a, one, zero};
  }
  if (a.IsZero()) {
    const Elem lc = b.LeadingCoeff();
    if (!lc.IsOne()) {
      const Elem inv_lc = lc.Inv();
      return {b * inv_lc, zero, Poly::FromRaw({inv_lc.RawValue()}, field, var)};
    }
    return {b, zero, one};
  }

  Poly r0 = a;
  Poly r1 = b;
  HalfGcdMatrix<F> mat = HalfGcdMatrix<F>::Identity(field, var);

  while (!r1.IsZero()) {
    if (r0.Degree() <= r1.Degree() ||
        static_cast<int64_t>(r0.Degree()) <= kHalfGcdCutoff) {
      auto [q, r_next] = detail::FastDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(r_next);
      mat.LeftMulStep(q);
    } else {
      HalfGcdMatrix<F> h = detail::HalfGcdRec(r0, r1);
      auto [next_r0, next_r1] = h.Apply(r0, r1);
      r0 = std::move(next_r0);
      r1 = std::move(next_r1);
      mat = HalfGcdMatrix<F>::Multiply(h, mat);
      if (!r1.IsZero()) {
        auto [q, r_next] = detail::FastDivRem(r0, r1);
        r0 = std::move(r1);
        r1 = std::move(r_next);
        mat.LeftMulStep(q);
      }
    }
  }

  Poly g = std::move(r0);
  Poly u = std::move(mat.m00);
  Poly v = std::move(mat.m01);
  if (!g.IsZero() && !g.IsMonic()) {
    const Elem inv_lc = g.LeadingCoeff().Inv();
    g *= inv_lc;
    u *= inv_lc;
    v *= inv_lc;
  }
  return {std::move(g), std::move(u), std::move(v)};
}

// Fast O(M(d) log d) Partial Extended Euclidean Algorithm via Half-GCD.
// Advances the Euclidean sequence starting at (r_0, r_1) = (a, b) and stops
// at the FIRST remainder r (starting with r_1 = b) satisfying
//   r.IsZero() || r.Degree() < stop_deg,
// returning {r, u, v} such that u * a + v * b == r.
// Matches Sugiyama's exact stopping condition for Reed-Solomon key equations.
template <FiniteField F>
inline std::tuple<PolynomialDense<F>, PolynomialDense<F>, PolynomialDense<F>>
PartialExtendedGcd(const PolynomialDense<F>& a, const PolynomialDense<F>& b,
                   int64_t stop_deg) {
  using Poly = PolynomialDense<F>;
  detail::ValidateSamePolyField(a, b);
  const auto field = a.Field();
  const std::string& var = a.GetVariable();
  const Poly zero = Poly::Zero(field, var);
  const Poly one = Poly::One(field, var);

  if (b.IsZero() || static_cast<int64_t>(b.Degree()) < stop_deg) {
    return {b, zero, one};
  }
  if (a.IsZero()) {
    return {zero, one, zero};
  }

  Poly r0 = a;
  Poly r1 = b;
  HalfGcdMatrix<F> mat = HalfGcdMatrix<F>::Identity(field, var);

  if (r0.Degree() <= r1.Degree()) {
    auto [q, r_next] = detail::FastDivRem(r0, r1);
    r0 = std::move(r1);
    r1 = std::move(r_next);
    mat.LeftMulStep(q);
  }

  while (!r1.IsZero() && static_cast<int64_t>(r1.Degree()) >= stop_deg) {
    const int64_t d = static_cast<int64_t>(r0.Degree());
    if (d <= kHalfGcdCutoff) {
      auto [q, r_next] = r0.DivRem(r1);
      r0 = std::move(r1);
      r1 = std::move(r_next);
      mat.LeftMulStep(q);
    } else if (stop_deg > 0 && 2 * stop_deg >= d) {
      // Exact single-shot Half-GCD jump directly to threshold stop_deg:
      // s = 2 * stop_deg - d satisfies 0 <= s < stop_deg <= deg(r1) < d and
      // ceil((d - s) / 2) + s == stop_deg.
      const int64_t s = 2 * stop_deg - d;
      HalfGcdMatrix<F> h =
          detail::HalfGcdRec(detail::TopShift(r0, s), detail::TopShift(r1, s));
      auto [next_r0, next_r1] = h.Apply(r0, r1);
      r0 = std::move(next_r0);
      r1 = std::move(next_r1);
      mat = HalfGcdMatrix<F>::Multiply(h, mat);
    } else {
      HalfGcdMatrix<F> h = detail::HalfGcdRec(r0, r1);
      auto [next_r0, next_r1] = h.Apply(r0, r1);
      r0 = std::move(next_r0);
      r1 = std::move(next_r1);
      mat = HalfGcdMatrix<F>::Multiply(h, mat);
      if (!r1.IsZero() && static_cast<int64_t>(r1.Degree()) >= stop_deg) {
        auto [q, r_next] = detail::FastDivRem(r0, r1);
        r0 = std::move(r1);
        r1 = std::move(r_next);
        mat.LeftMulStep(q);
      }
    }
  }

  return {std::move(r1), std::move(mat.m10), std::move(mat.m11)};
}

}  // namespace poly

using poly::FastExtendedGcd;
using poly::HalfGcd;
using poly::HalfGcdMatrix;
using poly::kHalfGcdCutoff;
using poly::PartialExtendedGcd;

}  // namespace xg

#endif  // XGALOIS_POLY_HALF_GCD_HPP_
