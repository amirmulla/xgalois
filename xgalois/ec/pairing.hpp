#ifndef XGALOIS_EC_PAIRING_HPP_
#define XGALOIS_EC_PAIRING_HPP_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/ec/elliptic_curve.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace ec {

namespace detail {

// Evaluates Miller's function f_{r, P}(Q) on a Short Weierstrass curve E.
// Returns std::nullopt if Q falls into the support of any intermediate line
// function (zero numerator or zero denominator).
template <Field F>
inline std::optional<typename F::element_type> MillerLoopTry(
    const EllipticCurve<F, CurveForm::kShortWeierstrass>& curve,
    const AffinePoint<F>& p, const AffinePoint<F>& q, const BigInt& r) {
  using E = typename F::element_type;
  if (r <= BigInt(0)) {
    throw std::invalid_argument("MillerLoop: order r must be positive");
  }
  const F& field = curve.Field();
  if (curve.IsIdentity(p) || curve.IsIdentity(q) || r.IsOne()) {
    return field.MultiplicativeIdentity();
  }

  AffinePoint<F> t = p;
  E num = field.MultiplicativeIdentity();
  E den = field.MultiplicativeIdentity();

  auto tangent_step = [&]() -> bool {
    if (curve.IsIdentity(t)) {
      return true;
    }
    if (field.IsZero(t.y)) {
      // Vertical tangent line at 2-torsion point T: x - x_T.
      const E v = field.Sub(q.x, t.x);
      if (field.IsZero(v)) {
        return false;
      }
      num = field.Mul(num, v);
      t = curve.Identity();
      return true;
    }
    const E x2 = FieldSqr(field, t.x);
    const E three_x2 = field.Add(field.Add(x2, x2), x2);
    const E lam =
        field.Div(field.Add(three_x2, curve.A()), field.Add(t.y, t.y));
    const E x_next = field.Sub(FieldSqr(field, lam), field.Add(t.x, t.x));
    const E y_next = field.Sub(field.Mul(lam, field.Sub(t.x, x_next)), t.y);

    const E line_val =
        field.Sub(field.Sub(q.y, t.y), field.Mul(lam, field.Sub(q.x, t.x)));
    const E vert_val = field.Sub(q.x, x_next);
    if (field.IsZero(line_val) || field.IsZero(vert_val)) {
      return false;
    }
    num = field.Mul(num, line_val);
    den = field.Mul(den, vert_val);
    t = AffinePoint<F>(x_next, y_next, false);
    return true;
  };

  auto chord_step = [&]() -> bool {
    if (curve.IsIdentity(t)) {
      t = p;
      return true;
    }
    if (alg::Equals(field, t.x, p.x)) {
      if (alg::Equals(field, t.y, p.y)) {
        return tangent_step();
      }
      // T == -P: vertical chord x - x_T, result T + P = O.
      const E v = field.Sub(q.x, t.x);
      if (field.IsZero(v)) {
        return false;
      }
      num = field.Mul(num, v);
      t = curve.Identity();
      return true;
    }
    const E lam = field.Div(field.Sub(p.y, t.y), field.Sub(p.x, t.x));
    const E x_next = field.Sub(field.Sub(FieldSqr(field, lam), t.x), p.x);
    const E y_next = field.Sub(field.Mul(lam, field.Sub(t.x, x_next)), t.y);

    const E line_val =
        field.Sub(field.Sub(q.y, t.y), field.Mul(lam, field.Sub(q.x, t.x)));
    const E vert_val = field.Sub(q.x, x_next);
    if (field.IsZero(line_val) || field.IsZero(vert_val)) {
      return false;
    }
    num = field.Mul(num, line_val);
    den = field.Mul(den, vert_val);
    t = AffinePoint<F>(x_next, y_next, false);
    return true;
  };

  const std::size_t bits = r.BitLength();
  for (std::size_t i = bits - 1; i-- > 0;) {
    num = FieldSqr(field, num);
    den = FieldSqr(field, den);
    if (!tangent_step()) {
      return std::nullopt;
    }
    if (r.GetBit(i)) {
      if (!chord_step()) {
        return std::nullopt;
      }
    }
  }

  if (field.IsZero(num) || field.IsZero(den)) {
    return std::nullopt;
  }
  return field.Div(num, den);
}

template <Field F>
inline std::optional<AffinePoint<F>> SampleCurvePoint(
    const EllipticCurve<F, CurveForm::kShortWeierstrass>& curve,
    uint64_t& seed) {
  const F& field = curve.Field();
  for (uint32_t tries = 0; tries < 256; ++tries) {
    const auto x = SampleFieldElement(field, seed);
    const auto x2 = FieldSqr(field, x);
    const auto rhs = field.Add(
        field.Add(field.Mul(x2, x), field.Mul(curve.A(), x)), curve.B());
    const auto y_opt = FieldSqrt(field, rhs);
    if (y_opt.has_value()) {
      return AffinePoint<F>(x, *y_opt, false);
    }
  }
  return std::nullopt;
}

}  // namespace detail

// ============================================================================
// Miller Loop: f_{r, P}(Q)
// ============================================================================

template <Field F>
inline typename F::element_type MillerLoop(
    const EllipticCurve<F, CurveForm::kShortWeierstrass>& curve,
    const AffinePoint<F>& p, const AffinePoint<F>& q, const BigInt& r) {
  const auto res = detail::MillerLoopTry(curve, p, q, r);
  if (!res.has_value()) {
    throw std::domain_error(
        "MillerLoop: evaluation point Q lies in the support of f_{r, P}");
  }
  return *res;
}

// ============================================================================
// Weil Pairing: e_r(P, Q) for P, Q in E[r]
// ============================================================================

template <Field F>
inline typename F::element_type WeilPairing(
    const EllipticCurve<F, CurveForm::kShortWeierstrass>& curve,
    const AffinePoint<F>& p, const AffinePoint<F>& q, const BigInt& r) {
  if (r <= BigInt(0)) {
    throw std::invalid_argument("WeilPairing: order r must be positive");
  }
  const F& field = curve.Field();
  if (curve.IsIdentity(p) || curve.IsIdentity(q) || curve.Equals(p, q)) {
    return field.MultiplicativeIdentity();
  }

  // Fast path: unshifted formula e_r(P, Q) = (-1)^r * f_{r, P}(Q) / f_{r, Q}(P)
  // whenever Q is disjoint from supp(f_{r, P}) and P from supp(f_{r, Q}).
  const auto fp_opt = detail::MillerLoopTry(curve, p, q, r);
  const auto fq_opt = detail::MillerLoopTry(curve, q, p, r);
  if (fp_opt.has_value() && fq_opt.has_value()) {
    auto ratio = field.Div(*fp_opt, *fq_opt);
    if (r.IsOdd()) {
      ratio = field.Neg(ratio);
    }
    return ratio;
  }

  // Fallback: shifted divisor formula
  //   e_r(P, Q) = (f_{r,P}(Q + S) / f_{r,P}(S)) / (f_{r,Q}(P - S) /
  //   f_{r,Q}(-S))
  uint64_t seed = 0x94d049bb133111ebULL;
  for (uint32_t attempt = 0; attempt < 256; ++attempt) {
    const auto s_opt = detail::SampleCurvePoint(curve, seed);
    if (!s_opt.has_value()) {
      continue;
    }
    const AffinePoint<F>& s = *s_opt;
    const AffinePoint<F> neg_s = curve.Neg(s);
    const AffinePoint<F> q_plus_s = curve.Add(q, s);
    const AffinePoint<F> p_minus_s = curve.Add(p, neg_s);
    if (curve.IsIdentity(s) || curve.IsIdentity(q_plus_s) ||
        curve.IsIdentity(p_minus_s)) {
      continue;
    }
    const auto f_p_qs = detail::MillerLoopTry(curve, p, q_plus_s, r);
    const auto f_p_s = detail::MillerLoopTry(curve, p, s, r);
    const auto f_q_ps = detail::MillerLoopTry(curve, q, p_minus_s, r);
    const auto f_q_ns = detail::MillerLoopTry(curve, q, neg_s, r);
    if (f_p_qs.has_value() && f_p_s.has_value() && f_q_ps.has_value() &&
        f_q_ns.has_value()) {
      const auto num = field.Mul(*f_p_qs, *f_q_ns);
      const auto den = field.Mul(*f_p_s, *f_q_ps);
      return field.Div(num, den);
    }
  }

  throw std::runtime_error(
      "WeilPairing: failed to find a valid auxiliary shift point S");
}

// ============================================================================
// Reduced Tate Pairing: t_r(P, Q) = f_{r, P}(Q)^((q^k - 1) / r)
// ============================================================================

template <Field F>
inline typename F::element_type TatePairing(
    const EllipticCurve<F, CurveForm::kShortWeierstrass>& curve,
    const AffinePoint<F>& p, const AffinePoint<F>& q, const BigInt& r,
    const BigInt& field_order = BigInt(0)) {
  if (r <= BigInt(0)) {
    throw std::invalid_argument("TatePairing: order r must be positive");
  }
  const F& field = curve.Field();
  const BigInt q_ord =
      (field_order > BigInt(0)) ? field_order : detail::FieldBigOrder(field);
  if (q_ord <= BigInt(1)) {
    throw std::domain_error(
        "TatePairing: finite field order must be > 1 (pass field_order "
        "explicitly if unavailable)");
  }
  const auto [exp, rem] = BigInt::DivRem(q_ord - BigInt(1), r);
  if (!rem.IsZero()) {
    throw std::invalid_argument("TatePairing: r must divide (field_order - 1)");
  }
  if (curve.IsIdentity(p) || curve.IsIdentity(q)) {
    return field.MultiplicativeIdentity();
  }

  // Fast path: direct f_{r, P}(Q)^((q^k - 1)/r).
  if (const auto direct = detail::MillerLoopTry(curve, p, q, r);
      direct.has_value()) {
    return detail::FieldPowBig(field, *direct, exp);
  }

  // Fallback: shifted divisor f_{r, P}(Q + S) / f_{r, P}(S).
  uint64_t seed = 0xbf58476d1ce4e5b9ULL;
  for (uint32_t attempt = 0; attempt < 256; ++attempt) {
    const auto s_opt = detail::SampleCurvePoint(curve, seed);
    if (!s_opt.has_value()) {
      continue;
    }
    const AffinePoint<F>& s = *s_opt;
    const AffinePoint<F> q_plus_s = curve.Add(q, s);
    if (curve.IsIdentity(s) || curve.IsIdentity(q_plus_s)) {
      continue;
    }
    const auto f_qs = detail::MillerLoopTry(curve, p, q_plus_s, r);
    const auto f_s = detail::MillerLoopTry(curve, p, s, r);
    if (f_qs.has_value() && f_s.has_value()) {
      return detail::FieldPowBig(field, field.Div(*f_qs, *f_s), exp);
    }
  }

  throw std::runtime_error(
      "TatePairing: failed to find a valid auxiliary shift point S");
}

// ============================================================================
// Elliptic Riemann-Roch Evaluation Matrix for AG Codes
// ============================================================================

// Constructs the evaluation matrix of the Riemann-Roch space L(m * O) on a
// Short Weierstrass elliptic curve E(F) at affine points `eval_points`.
// Basis functions are ordered by pole order 2*i + 3*j <= m (j in {0, 1}):
//   {1, x, y, x^2, x*y, x^3, x^2*y, ...}
// Dimension is 1 for m <= 1 and m for m >= 2.
template <Field F>
inline GFMatrix<F> EllipticRiemannRochMatrix(
    const EllipticCurve<F, CurveForm::kShortWeierstrass>& curve,
    const std::vector<AffinePoint<F>>& eval_points, uint32_t m) {
  using RawArray = typename GFMatrix<F>::raw_array_type;
  using ShapeType = typename GFMatrix<F>::shape_type;

  const F& field = curve.Field();
  for (const auto& pt : eval_points) {
    if (curve.IsIdentity(pt) || !curve.IsOnCurve(pt)) {
      throw std::invalid_argument(
          "EllipticRiemannRochMatrix: all evaluation points must be affine "
          "points on the curve");
    }
  }

  // Enumerate (i, j) pairs by pole order 2*i + 3*j <= m with j in {0, 1}.
  std::vector<std::pair<uint32_t, uint32_t>> basis;
  basis.emplace_back(0U, 0U);  // Pole order 0: 1.
  for (uint32_t pole = 2; pole <= m; ++pole) {
    if ((pole & 1U) == 0U) {
      basis.emplace_back(pole / 2U, 0U);
    } else {
      basis.emplace_back((pole - 3U) / 2U, 1U);
    }
  }

  const std::size_t rows = basis.size();
  const std::size_t cols = eval_points.size();
  RawArray raw = RawArray::from_shape(ShapeType{rows, cols});
  auto* data = raw.data();

  const uint32_t max_x_pow = (m >= 2) ? (m / 2U) : 0U;
  std::vector<typename F::element_type> x_pows(max_x_pow + 1);
  for (std::size_t c = 0; c < cols; ++c) {
    const auto& pt = eval_points[c];
    x_pows[0] = field.MultiplicativeIdentity();
    for (uint32_t p = 1; p <= max_x_pow; ++p) {
      x_pows[p] = field.Mul(x_pows[p - 1], pt.x);
    }
    for (std::size_t r = 0; r < rows; ++r) {
      const auto [i_pow, j_pow] = basis[r];
      data[r * cols + c] =
          (j_pow == 0U) ? x_pows[i_pow] : field.Mul(x_pows[i_pow], pt.y);
    }
  }

  auto field_ptr = std::make_shared<F>(field);
  return GFMatrix<F>::FromRaw(std::move(field_ptr), std::move(raw));
}

}  // namespace ec
}  // namespace xg

#endif  // XGALOIS_EC_PAIRING_HPP_
