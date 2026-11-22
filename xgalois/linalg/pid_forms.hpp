#ifndef XGALOIS_LINALG_PID_FORMS_HPP_
#define XGALOIS_LINALG_PID_FORMS_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace linalg {

// ============================================================================
// Hermite Normal Form & Smith Normal Form over Arbitrary Euclidean Domains
// ============================================================================

template <EuclideanDomain R>
struct HnfResult {
  using element_type = typename R::element_type;

  std::vector<std::vector<element_type>> H;
  std::vector<std::vector<element_type>> U;
  std::vector<std::size_t> pivot_cols;
  std::size_t rank = 0;
};

template <EuclideanDomain R>
struct PidSnfResult {
  using element_type = typename R::element_type;

  std::vector<std::vector<element_type>> S;
  std::vector<std::vector<element_type>> U;
  std::vector<std::vector<element_type>> V;
  std::vector<element_type> invariant_factors;
  std::size_t rank = 0;
};

namespace detail {

template <EuclideanDomain R>
inline std::vector<std::vector<typename R::element_type>> MakeRingIdentity(
    const R& ring, std::size_t n) {
  using E = typename R::element_type;
  const E zero = ring.AdditiveIdentity();
  const E one = ring.MultiplicativeIdentity();
  std::vector<std::vector<E>> eye(n, std::vector<E>(n, zero));
  for (std::size_t i = 0; i < n; ++i) {
    eye[i][i] = one;
  }
  return eye;
}

template <EuclideanDomain R>
inline void ValidateRectangularRingMatrix(
    const std::vector<std::vector<typename R::element_type>>& matrix,
    const char* fn_name) {
  if (matrix.empty()) {
    return;
  }
  const std::size_t cols = matrix[0].size();
  for (std::size_t i = 1; i < matrix.size(); ++i) {
    if (matrix[i].size() != cols) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": input matrix must be rectangular");
    }
  }
}

// Lightweight raw polynomial operations over a FiniteField F for Frobenius
// Normal Form, MinimalPoly, and CharacteristicPoly without depending on
// xgalois/ring/polynomial_ring.hpp.
template <FiniteField F>
struct RawPolyOps {
  using E = typename F::element_type;
  using Poly = std::vector<E>;

  static void Trim(const F& field, Poly& p) {
    while (!p.empty() && field.IsZero(p.back())) {
      p.pop_back();
    }
  }

  static bool IsZero(const Poly& p) noexcept { return p.empty(); }

  static int64_t Degree(const Poly& p) noexcept {
    return static_cast<int64_t>(p.size()) - 1;
  }

  static Poly Constant(const F& field, const E& c) {
    if (field.IsZero(c)) {
      return {};
    }
    return {c};
  }

  static Poly One(const F& field) { return {field.MultiplicativeIdentity()}; }

  static Poly X(const F& field) {
    return {field.AdditiveIdentity(), field.MultiplicativeIdentity()};
  }

  static Poly Add(const F& field, const Poly& a, const Poly& b) {
    const std::size_t len = std::max(a.size(), b.size());
    const E zero = field.AdditiveIdentity();
    Poly out(len, zero);
    for (std::size_t i = 0; i < len; ++i) {
      const E& ca = (i < a.size()) ? a[i] : zero;
      const E& cb = (i < b.size()) ? b[i] : zero;
      out[i] = field.Add(ca, cb);
    }
    Trim(field, out);
    return out;
  }

  static Poly Sub(const F& field, const Poly& a, const Poly& b) {
    const std::size_t len = std::max(a.size(), b.size());
    const E zero = field.AdditiveIdentity();
    Poly out(len, zero);
    for (std::size_t i = 0; i < len; ++i) {
      const E& ca = (i < a.size()) ? a[i] : zero;
      const E& cb = (i < b.size()) ? b[i] : zero;
      out[i] = field.Sub(ca, cb);
    }
    Trim(field, out);
    return out;
  }

  static Poly Neg(const F& field, const Poly& a) {
    Poly out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      out[i] = field.Neg(a[i]);
    }
    Trim(field, out);
    return out;
  }

  static Poly Scale(const F& field, const Poly& a, const E& scalar) {
    if (field.IsZero(scalar) || a.empty()) {
      return {};
    }
    Poly out(a.size());
    field.ScalarMulArray(std::span<E>(out), std::span<const E>(a), scalar);
    Trim(field, out);
    return out;
  }

  static Poly Mul(const F& field, const Poly& a, const Poly& b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    const E zero = field.AdditiveIdentity();
    Poly out(a.size() + b.size() - 1, zero);
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (field.IsZero(a[i])) {
        continue;
      }
      field.FmaArray(std::span<E>(out.data() + i, b.size()),
                     std::span<const E>(b), a[i]);
    }
    Trim(field, out);
    return out;
  }

  static std::pair<Poly, Poly> DivRem(const F& field, const Poly& a,
                                      const Poly& b) {
    if (b.empty()) {
      throw std::domain_error("RawPolyOps::DivRem: division by zero");
    }
    Poly rem = a;
    Trim(field, rem);
    if (rem.size() < b.size()) {
      return {{}, std::move(rem)};
    }
    const std::size_t db = b.size() - 1;
    const E inv_lc = field.Inv(b[db]);
    const std::size_t q_len = rem.size() - b.size() + 1;
    Poly q(q_len, field.AdditiveIdentity());

    for (std::size_t shift = q_len; shift-- > 0;) {
      const E top = rem[shift + db];
      if (field.IsZero(top)) {
        continue;
      }
      const E factor = field.Mul(top, inv_lc);
      q[shift] = factor;
      const E neg_factor = field.Neg(factor);
      field.FmaArray(std::span<E>(rem.data() + shift, db + 1),
                     std::span<const E>(b), neg_factor);
    }
    Trim(field, q);
    Trim(field, rem);
    return {std::move(q), std::move(rem)};
  }

  static Poly MakeMonic(const F& field, Poly p) {
    Trim(field, p);
    if (p.empty()) {
      return p;
    }
    if (field.IsOne(p.back())) {
      return p;
    }
    const E inv_lc = field.Inv(p.back());
    field.ScalarMulArray(std::span<E>(p), std::span<const E>(p), inv_lc);
    p.back() = field.MultiplicativeIdentity();
    return p;
  }

  static Poly Gcd(const F& field, Poly a, Poly b) {
    Trim(field, a);
    Trim(field, b);
    while (!b.empty()) {
      Poly rem = DivRem(field, a, b).second;
      a = std::move(b);
      b = std::move(rem);
    }
    return MakeMonic(field, std::move(a));
  }

  static Poly Lcm(const F& field, const Poly& a, const Poly& b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    const Poly g = Gcd(field, a, b);
    const Poly q = DivRem(field, a, g).first;
    return MakeMonic(field, Mul(field, q, b));
  }
};

// Computes y = A * x where A is n x n row-major raw buffer.
template <FiniteField F>
inline void RawMatVec(const F& field, const typename F::element_type* a_data,
                      std::size_t n,
                      std::span<const typename F::element_type> x,
                      std::span<typename F::element_type> y) {
  for (std::size_t i = 0; i < n; ++i) {
    y[i] = field.DotProduct(
        std::span<const typename F::element_type>(a_data + i * n, n), x);
  }
}

// Evaluates p(A) * v using Horner's scheme: p(x) = sum c_k x^k.
template <FiniteField F>
inline std::vector<typename F::element_type> EvalPolyMatVec(
    const F& field, const typename F::element_type* a_data, std::size_t n,
    const std::vector<typename F::element_type>& poly,
    std::span<const typename F::element_type> v) {
  using E = typename F::element_type;
  const E zero = field.AdditiveIdentity();
  std::vector<E> acc(n, zero);
  if (poly.empty()) {
    return acc;
  }
  std::vector<E> tmp(n, zero);
  field.ScalarMulArray(std::span<E>(acc), v, poly.back());
  for (std::size_t k = poly.size() - 1; k-- > 0;) {
    RawMatVec(field, a_data, n, std::span<const E>(acc), std::span<E>(tmp));
    acc.swap(tmp);
    if (!field.IsZero(poly[k])) {
      field.FmaArray(std::span<E>(acc), v, poly[k]);
    }
  }
  return acc;
}

// Incremental echelon subspace tracker over F^n.
template <FiniteField F>
class IncrementalSubspace {
 public:
  using E = typename F::element_type;

  IncrementalSubspace(const F& field, std::size_t dim)
      : field_(&field),
        dim_(dim),
        pivot_row_for_col_(dim, kNoPivot),
        basis_rows_() {}

  // Reduces v against the current basis. Returns true if v is already in the
  // span (reduces to 0), false if v is linearly independent (and if insert is
  // true, adds the reduced vector normalized to unit pivot).
  bool ReduceAndMaybeInsert(std::vector<E> v, bool insert) {
    for (std::size_t col = 0; col < dim_; ++col) {
      if (field_->IsZero(v[col])) {
        continue;
      }
      const std::size_t p_idx = pivot_row_for_col_[col];
      if (p_idx != kNoPivot) {
        const E neg_factor = field_->Neg(v[col]);
        v[col] = field_->AdditiveIdentity();
        if (col + 1 < dim_) {
          field_->FmaArray(
              std::span<E>(v.data() + col + 1, dim_ - col - 1),
              std::span<const E>(basis_rows_[p_idx].data() + col + 1,
                                 dim_ - col - 1),
              neg_factor);
        }
      } else {
        if (insert) {
          const E inv_piv = field_->Inv(v[col]);
          v[col] = field_->MultiplicativeIdentity();
          if (col + 1 < dim_) {
            std::span<E> tail(v.data() + col + 1, dim_ - col - 1);
            field_->ScalarMulArray(tail, tail, inv_piv);
          }
          pivot_row_for_col_[col] = basis_rows_.size();
          basis_rows_.push_back(std::move(v));
        }
        return false;
      }
    }
    return true;
  }

  std::size_t Dimension() const noexcept { return basis_rows_.size(); }

 private:
  static constexpr std::size_t kNoPivot = static_cast<std::size_t>(-1);

  const F* field_;
  std::size_t dim_;
  std::vector<std::size_t> pivot_row_for_col_;
  std::vector<std::vector<E>> basis_rows_;
};

// Computes the exact minimal annihilating polynomial mu_{A, v}(x) of vector v
// under A, and optionally appends all Krylov vectors {v, Av, ..., A^{d-1}v}
// into the provided subspace tracker.
template <FiniteField F>
inline std::vector<typename F::element_type> VectorMinimalPoly(
    const F& field, const typename F::element_type* a_data, std::size_t n,
    const std::vector<typename F::element_type>& v0,
    IncrementalSubspace<F>* global_subspace = nullptr) {
  using E = typename F::element_type;
  const E zero = field.AdditiveIdentity();
  const E one = field.MultiplicativeIdentity();

  constexpr std::size_t kNoPivot = static_cast<std::size_t>(-1);
  std::vector<std::size_t> pivot_row_for_col(n, kNoPivot);
  std::vector<std::vector<E>> echelon_rows;
  std::vector<std::vector<E>>
      comb_rows;  // Combination coefficients in Krylov basis.

  std::vector<E> cur = v0;
  std::vector<E> next(n, zero);

  for (std::size_t d = 0; d <= n; ++d) {
    std::vector<E> reduced = cur;
    std::vector<E> coeffs(d + 1, zero);
    coeffs[d] = one;

    std::size_t new_pivot_col = kNoPivot;
    for (std::size_t col = 0; col < n; ++col) {
      if (field.IsZero(reduced[col])) {
        continue;
      }
      const std::size_t p_idx = pivot_row_for_col[col];
      if (p_idx != kNoPivot) {
        const E neg_factor = field.Neg(reduced[col]);
        reduced[col] = zero;
        if (col + 1 < n) {
          field.FmaArray(std::span<E>(reduced.data() + col + 1, n - col - 1),
                         std::span<const E>(
                             echelon_rows[p_idx].data() + col + 1, n - col - 1),
                         neg_factor);
        }
        const auto& c_row = comb_rows[p_idx];
        field.FmaArray(std::span<E>(coeffs.data(), c_row.size()),
                       std::span<const E>(c_row), neg_factor);
      } else {
        new_pivot_col = col;
        break;
      }
    }

    if (new_pivot_col == kNoPivot) {
      // cur is a linear combination of previous Krylov vectors:
      // sum_{k=0}^{d} coeffs[k] * A^k v0 == 0 with coeffs[d] == 1.
      RawPolyOps<F>::Trim(field, coeffs);
      return coeffs;
    }

    // Normalize reduced row and its combination coefficients.
    const E inv_piv = field.Inv(reduced[new_pivot_col]);
    reduced[new_pivot_col] = one;
    if (new_pivot_col + 1 < n) {
      std::span<E> tail(reduced.data() + new_pivot_col + 1,
                        n - new_pivot_col - 1);
      field.ScalarMulArray(tail, tail, inv_piv);
    }
    field.ScalarMulArray(std::span<E>(coeffs), std::span<const E>(coeffs),
                         inv_piv);

    pivot_row_for_col[new_pivot_col] = echelon_rows.size();
    echelon_rows.push_back(std::move(reduced));
    comb_rows.push_back(std::move(coeffs));

    if (global_subspace != nullptr) {
      global_subspace->ReduceAndMaybeInsert(cur, /*insert=*/true);
    }

    RawMatVec(field, a_data, n, std::span<const E>(cur), std::span<E>(next));
    cur.swap(next);
  }

  throw std::runtime_error("VectorMinimalPoly: unexpected non-termination");
}

}  // namespace detail

// Exact row Hermite Normal Form over ANY EuclideanDomain R.
// Computes unimodular U (m x m) and row-HNF H (m x n) such that U * A == H,
// every pivot H[r][pivot_cols[r]] is canonicalized via alg::AssociateCanonical,
// and every entry above a pivot is reduced modulo that pivot via alg::DivRem.
template <EuclideanDomain R>
HnfResult<R> HermiteNormalForm(
    const R& ring,
    const std::vector<std::vector<typename R::element_type>>& matrix) {
  using E = typename R::element_type;

  detail::ValidateRectangularRingMatrix<R>(matrix, "HermiteNormalForm");

  HnfResult<R> result;
  const std::size_t m = matrix.size();
  if (m == 0) {
    return result;
  }
  const std::size_t n = matrix[0].size();
  result.H = matrix;
  result.U = detail::MakeRingIdentity(ring, m);
  if (n == 0) {
    return result;
  }

  auto& H = result.H;
  auto& U = result.U;

  auto swap_rows = [&](std::size_t r1, std::size_t r2) {
    if (r1 != r2) {
      std::swap(H[r1], H[r2]);
      std::swap(U[r1], U[r2]);
    }
  };

  auto scale_row = [&](std::size_t r, const E& unit) {
    if (ring.IsOne(unit)) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!ring.IsZero(H[r][j])) {
        H[r][j] = ring.Mul(unit, H[r][j]);
      }
    }
    for (std::size_t j = 0; j < m; ++j) {
      if (!ring.IsZero(U[r][j])) {
        U[r][j] = ring.Mul(unit, U[r][j]);
      }
    }
  };

  auto sub_row_mul = [&](std::size_t dst, std::size_t src, const E& factor) {
    if (ring.IsZero(factor)) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!ring.IsZero(H[src][j])) {
        H[dst][j] = ring.Sub(H[dst][j], ring.Mul(factor, H[src][j]));
      }
    }
    for (std::size_t j = 0; j < m; ++j) {
      if (!ring.IsZero(U[src][j])) {
        U[dst][j] = ring.Sub(U[dst][j], ring.Mul(factor, U[src][j]));
      }
    }
  };

  std::size_t pivot_row = 0;
  for (std::size_t col = 0; col < n && pivot_row < m; ++col) {
    // Reduce column `col` in rows [pivot_row, m) down to a single GCD pivot.
    while (true) {
      bool found = false;
      std::size_t best_r = pivot_row;
      decltype(alg::EuclideanNorm(ring, H[0][0])) best_norm{};

      for (std::size_t i = pivot_row; i < m; ++i) {
        if (!ring.IsZero(H[i][col])) {
          auto norm = alg::EuclideanNorm(ring, H[i][col]);
          if (!found || norm < best_norm) {
            found = true;
            best_r = i;
            best_norm = std::move(norm);
          }
        }
      }

      if (!found) {
        break;
      }

      swap_rows(pivot_row, best_r);

      bool any_remaining = false;
      for (std::size_t i = pivot_row + 1; i < m; ++i) {
        if (!ring.IsZero(H[i][col])) {
          auto [q, rem] = alg::DivRem(ring, H[i][col], H[pivot_row][col]);
          (void)rem;
          sub_row_mul(i, pivot_row, q);
          if (!ring.IsZero(H[i][col])) {
            any_remaining = true;
          }
        }
      }

      if (!any_remaining) {
        break;
      }
    }

    if (ring.IsZero(H[pivot_row][col])) {
      continue;
    }

    // Canonicalize pivot entry via unit multiplication.
    auto [canon_piv, unit] = alg::AssociateCanonical(ring, H[pivot_row][col]);
    (void)canon_piv;
    scale_row(pivot_row, unit);

    // Reduce entries above the pivot modulo the canonical pivot.
    for (std::size_t i = 0; i < pivot_row; ++i) {
      if (!ring.IsZero(H[i][col])) {
        auto [q, rem] = alg::DivRem(ring, H[i][col], H[pivot_row][col]);
        (void)rem;
        if (!ring.IsZero(q)) {
          sub_row_mul(i, pivot_row, q);
        }
      }
    }

    result.pivot_cols.push_back(col);
    ++pivot_row;
  }

  result.rank = pivot_row;
  return result;
}

// Exact Smith Normal Form over ANY EuclideanDomain R.
// Computes unimodular U (m x m), unimodular V (n x n), and diagonal S (m x n)
// such that U * A * V == S with canonicalized invariant factors
// d_1 | d_2 | ... | d_r on the diagonal.
template <EuclideanDomain R>
PidSnfResult<R> SmithNormalForm(
    const R& ring,
    const std::vector<std::vector<typename R::element_type>>& matrix) {
  using E = typename R::element_type;

  detail::ValidateRectangularRingMatrix<R>(matrix, "SmithNormalForm");

  PidSnfResult<R> result;
  const std::size_t m = matrix.size();
  if (m == 0) {
    return result;
  }
  const std::size_t n = matrix[0].size();
  result.S = matrix;
  result.U = detail::MakeRingIdentity(ring, m);
  result.V = detail::MakeRingIdentity(ring, n);
  if (n == 0) {
    return result;
  }

  auto& S = result.S;
  auto& U = result.U;
  auto& V = result.V;

  auto swap_rows = [&](std::size_t r1, std::size_t r2) {
    if (r1 != r2) {
      std::swap(S[r1], S[r2]);
      std::swap(U[r1], U[r2]);
    }
  };

  auto swap_cols = [&](std::size_t c1, std::size_t c2) {
    if (c1 != c2) {
      for (std::size_t i = 0; i < m; ++i) {
        std::swap(S[i][c1], S[i][c2]);
      }
      for (std::size_t i = 0; i < n; ++i) {
        std::swap(V[i][c1], V[i][c2]);
      }
    }
  };

  auto scale_row = [&](std::size_t r, const E& unit) {
    if (ring.IsOne(unit)) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!ring.IsZero(S[r][j])) {
        S[r][j] = ring.Mul(unit, S[r][j]);
      }
    }
    for (std::size_t j = 0; j < m; ++j) {
      if (!ring.IsZero(U[r][j])) {
        U[r][j] = ring.Mul(unit, U[r][j]);
      }
    }
  };

  auto add_row_mul = [&](std::size_t dst, std::size_t src, const E& factor) {
    if (ring.IsZero(factor)) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!ring.IsZero(S[src][j])) {
        S[dst][j] = ring.Add(S[dst][j], ring.Mul(factor, S[src][j]));
      }
    }
    for (std::size_t j = 0; j < m; ++j) {
      if (!ring.IsZero(U[src][j])) {
        U[dst][j] = ring.Add(U[dst][j], ring.Mul(factor, U[src][j]));
      }
    }
  };

  auto add_col_mul = [&](std::size_t dst, std::size_t src, const E& factor) {
    if (ring.IsZero(factor)) {
      return;
    }
    for (std::size_t i = 0; i < m; ++i) {
      if (!ring.IsZero(S[i][src])) {
        S[i][dst] = ring.Add(S[i][dst], ring.Mul(factor, S[i][src]));
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (!ring.IsZero(V[i][src])) {
        V[i][dst] = ring.Add(V[i][dst], ring.Mul(factor, V[i][src]));
      }
    }
  };

  const std::size_t limit = std::min(m, n);
  std::size_t r = 0;

  while (r < limit) {
    // 1. Select the nonzero entry in S[r..m-1][r..n-1] with minimal norm.
    bool found_pivot = false;
    std::size_t best_r = r;
    std::size_t best_c = r;
    decltype(alg::EuclideanNorm(ring, S[0][0])) best_norm{};

    for (std::size_t i = r; i < m; ++i) {
      for (std::size_t j = r; j < n; ++j) {
        if (!ring.IsZero(S[i][j])) {
          auto norm = alg::EuclideanNorm(ring, S[i][j]);
          if (!found_pivot || norm < best_norm) {
            found_pivot = true;
            best_r = i;
            best_c = j;
            best_norm = std::move(norm);
          }
        }
      }
    }

    if (!found_pivot) {
      break;
    }

    swap_rows(r, best_r);
    swap_cols(r, best_c);

    // 2. Eliminate column r and row r and enforce divisibility of trailing
    // block.
    while (true) {
      while (true) {
        // Clear column r below row r.
        bool col_has_nonzero = false;
        std::size_t min_row = r;
        decltype(alg::EuclideanNorm(ring, S[0][0])) min_col_norm{};
        for (std::size_t i = r + 1; i < m; ++i) {
          if (!ring.IsZero(S[i][r])) {
            auto [q, rem] = alg::DivRem(ring, S[i][r], S[r][r]);
            (void)rem;
            add_row_mul(i, r, ring.Neg(q));
            if (!ring.IsZero(S[i][r])) {
              auto norm = alg::EuclideanNorm(ring, S[i][r]);
              if (!col_has_nonzero || norm < min_col_norm) {
                col_has_nonzero = true;
                min_row = i;
                min_col_norm = std::move(norm);
              }
            }
          }
        }
        if (col_has_nonzero) {
          swap_rows(r, min_row);
          continue;
        }

        // Clear row r right of column r.
        bool row_has_nonzero = false;
        std::size_t min_col = r;
        decltype(alg::EuclideanNorm(ring, S[0][0])) min_row_norm{};
        for (std::size_t j = r + 1; j < n; ++j) {
          if (!ring.IsZero(S[r][j])) {
            auto [q, rem] = alg::DivRem(ring, S[r][j], S[r][r]);
            (void)rem;
            add_col_mul(j, r, ring.Neg(q));
            if (!ring.IsZero(S[r][j])) {
              auto norm = alg::EuclideanNorm(ring, S[r][j]);
              if (!row_has_nonzero || norm < min_row_norm) {
                row_has_nonzero = true;
                min_col = j;
                min_row_norm = std::move(norm);
              }
            }
          }
        }
        if (row_has_nonzero) {
          swap_cols(r, min_col);
          continue;
        }

        break;
      }

      // Check divisibility of trailing submatrix S[r+1..m-1][r+1..n-1].
      bool found_non_multiple = false;
      std::size_t bad_row = r;
      for (std::size_t i = r + 1; i < m && !found_non_multiple; ++i) {
        for (std::size_t j = r + 1; j < n; ++j) {
          if (!ring.IsZero(S[i][j])) {
            auto [q, rem] = alg::DivRem(ring, S[i][j], S[r][r]);
            (void)q;
            if (!ring.IsZero(rem)) {
              found_non_multiple = true;
              bad_row = i;
              break;
            }
          }
        }
      }

      if (!found_non_multiple) {
        break;
      }

      add_row_mul(r, bad_row, ring.MultiplicativeIdentity());
    }

    // Canonicalize pivot S[r][r].
    auto [canon, unit] = alg::AssociateCanonical(ring, S[r][r]);
    (void)canon;
    scale_row(r, unit);

    result.invariant_factors.push_back(S[r][r]);
    ++r;
  }

  result.rank = r;
  return result;
}

// ============================================================================
// Canonical Forms over Finite Fields: MinimalPoly, CharacteristicPoly,
// FrobeniusNormalForm (Rational Canonical Form)
// ============================================================================

template <FiniteField F>
struct FrobeniusNormalFormResult {
  std::vector<PolynomialDense<F>> invariant_factors;
  GFMatrix<F> rational_canonical_form;
  GFMatrix<F> similarity_matrix;
};

// O(n^3) exact deterministic Krylov minimal polynomial over any FiniteField F.
template <FiniteField F>
PolynomialDense<F> MinimalPoly(const GFMatrix<F>& a,
                               const std::string& var = "x") {
  using E = typename F::element_type;
  detail::RequireValidArray(a, "MinimalPoly");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("MinimalPoly requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  if (n == 0) {
    return PolynomialDense<F>::One(std::move(field), var);
  }

  const E* a_data = a.RawData().data();
  const E zero = field->AdditiveIdentity();
  const E one = field->MultiplicativeIdentity();

  detail::IncrementalSubspace<F> spanned(*field, n);
  std::vector<E> mu = detail::RawPolyOps<F>::One(*field);

  // Start from a deterministic pseudo-random vector so generic matrices finish
  // in a single Krylov pass.
  std::mt19937_64 rng(0x9E3779B97F4A7C15ULL ^ static_cast<uint64_t>(n));
  std::vector<E> seed_v(n);
  bool seed_nonzero = false;
  for (std::size_t i = 0; i < n; ++i) {
    seed_v[i] = field->Random(rng);
    if (!field->IsZero(seed_v[i])) {
      seed_nonzero = true;
    }
  }
  if (seed_nonzero) {
    mu = detail::VectorMinimalPoly(*field, a_data, n, seed_v, &spanned);
  }

  // Sweep standard basis vectors e_i not yet in the A-invariant span W.
  for (std::size_t i = 0; i < n && spanned.Dimension() < n; ++i) {
    if (mu.size() == n + 1) {
      break;
    }
    std::vector<E> ei(n, zero);
    ei[i] = one;
    if (spanned.ReduceAndMaybeInsert(ei, /*insert=*/false)) {
      continue;
    }
    auto mu_i = detail::VectorMinimalPoly(*field, a_data, n, ei, &spanned);
    mu = detail::RawPolyOps<F>::Lcm(*field, mu, mu_i);
  }

  return PolynomialDense<F>::FromRaw(std::move(mu), std::move(field), var);
}

// O(n^3) exact Hessenberg-reduction characteristic polynomial det(xI - A) over
// any FiniteField F. Strictly more constrained than the unconstrained Berkowitz
// overload in linalg.hpp so overload resolution unambiguously selects this
// O(n^3) implementation for all FiniteField types.
template <FiniteField F>
PolynomialDense<F> CharacteristicPoly(const GFMatrix<F>& a,
                                      const std::string& var = "x") {
  using E = typename F::element_type;
  detail::RequireValidArray(a, "CharacteristicPoly");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("CharacteristicPoly requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  if (n == 0) {
    return PolynomialDense<F>::One(std::move(field), var);
  }

  const E zero = field->AdditiveIdentity();
  const E one = field->MultiplicativeIdentity();

  // Copy A into a contiguous row-major buffer and reduce to upper Hessenberg
  // form via exact similarity transformations in O(n^3).
  std::vector<E> h(a.RawData().data(), a.RawData().data() + n * n);

  for (std::size_t k = 0; k + 2 < n; ++k) {
    std::size_t pivot = k + 1;
    while (pivot < n && field->IsZero(h[pivot * n + k])) {
      ++pivot;
    }
    if (pivot == n) {
      continue;
    }
    if (pivot != k + 1) {
      // Similarity swap: swap rows (k+1, pivot) and columns (k+1, pivot).
      for (std::size_t j = k; j < n; ++j) {
        std::swap(h[(k + 1) * n + j], h[pivot * n + j]);
      }
      for (std::size_t i = 0; i < n; ++i) {
        std::swap(h[i * n + (k + 1)], h[i * n + pivot]);
      }
    }

    const E inv_piv = field->Inv(h[(k + 1) * n + k]);
    const std::size_t tail_len = n - k;
    std::span<const E> piv_row(h.data() + (k + 1) * n + k, tail_len);

    for (std::size_t i = k + 2; i < n; ++i) {
      const E entry = h[i * n + k];
      if (field->IsZero(entry)) {
        continue;
      }
      const E c = field->Mul(entry, inv_piv);
      const E neg_c = field->Neg(c);
      // Left operation: row[i] -= c * row[k + 1]
      field->FmaArray(std::span<E>(h.data() + i * n + k, tail_len), piv_row,
                      neg_c);
      h[i * n + k] = zero;
      // Right inverse operation: col[k + 1] += c * col[i]
      for (std::size_t row = 0; row < n; ++row) {
        const E val = h[row * n + i];
        if (!field->IsZero(val)) {
          h[row * n + (k + 1)] =
              field->Add(h[row * n + (k + 1)], field->Mul(c, val));
        }
      }
    }
  }

  // Compute p_k(x) = det(xI_k - H_{0..k-1, 0..k-1}) in O(n^3):
  //   p_0(x) = 1
  //   p_k(x) = (x - H[k-1][k-1]) * p_{k-1}(x)
  //            - sum_{r=0}^{k-2} H[r][k-1] * (prod_{j=r+1}^{k-1} H[j][j-1]) *
  //            p_r(x)
  std::vector<std::vector<E>> p(n + 1);
  p[0] = {one};

  for (std::size_t k = 1; k <= n; ++k) {
    std::vector<E> pk(k + 1, zero);
    // x * p_{k-1}(x)
    for (std::size_t d = 0; d < k; ++d) {
      pk[d + 1] = p[k - 1][d];
    }
    // - H[k-1][k-1] * p_{k-1}(x)
    const E neg_diag = field->Neg(h[(k - 1) * n + (k - 1)]);
    if (!field->IsZero(neg_diag)) {
      field->FmaArray(std::span<E>(pk.data(), k), std::span<const E>(p[k - 1]),
                      neg_diag);
    }

    E subdiag_prod = one;
    for (std::size_t r = k - 1; r-- > 0;) {
      subdiag_prod = field->Mul(subdiag_prod, h[(r + 1) * n + r]);
      if (field->IsZero(subdiag_prod)) {
        break;
      }
      const E h_rk = h[r * n + (k - 1)];
      if (field->IsZero(h_rk)) {
        continue;
      }
      const E neg_coeff = field->Neg(field->Mul(h_rk, subdiag_prod));
      field->FmaArray(std::span<E>(pk.data(), r + 1), std::span<const E>(p[r]),
                      neg_coeff);
    }
    p[k] = std::move(pk);
  }

  return PolynomialDense<F>::FromRaw(std::move(p[n]), std::move(field), var);
}

// Exact Frobenius Normal Form (Rational Canonical Form) over any FiniteField F.
// Computes monic invariant factors f_1 | f_2 | ... | f_s, block-diagonal
// companion matrix RCF, and invertible similarity matrix P with
// P^{-1} * A * P == RCF.
template <FiniteField F>
FrobeniusNormalFormResult<F> FrobeniusNormalForm(const GFMatrix<F>& a,
                                                 const std::string& var = "x") {
  using E = typename F::element_type;
  using Poly = std::vector<E>;
  using Ops = detail::RawPolyOps<F>;

  detail::RequireValidArray(a, "FrobeniusNormalForm");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("FrobeniusNormalForm requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  FrobeniusNormalFormResult<F> result;
  if (n == 0) {
    result.rational_canonical_form = GFMatrix<F>::Zeros(field, {0, 0});
    result.similarity_matrix = GFMatrix<F>::Zeros(field, {0, 0});
    return result;
  }

  const E zero = field->AdditiveIdentity();
  const E one = field->MultiplicativeIdentity();
  const E* a_data = a.RawData().data();

  // Build M(x) = xI - A over F[x].
  std::vector<std::vector<Poly>> s(n, std::vector<Poly>(n));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const E neg_aij = field->Neg(a_data[i * n + j]);
      if (i == j) {
        s[i][j] = {neg_aij, one};
      } else if (!field->IsZero(neg_aij)) {
        s[i][j] = {neg_aij};
      }
    }
  }

  // Track W_eval = Phi(U(x)^{-1}) in F^{n x n} directly as columns w_cols[k]
  // in F^n, where Phi(sum c_t x^t) = sum A^t c_t. Initially U(x) = I_n so
  // w_cols[k] = e_k.
  std::vector<std::vector<E>> w_cols(n, std::vector<E>(n, zero));
  for (std::size_t k = 0; k < n; ++k) {
    w_cols[k][k] = one;
  }

  auto swap_rows = [&](std::size_t r1, std::size_t r2) {
    if (r1 != r2) {
      std::swap(s[r1], s[r2]);
      std::swap(w_cols[r1], w_cols[r2]);
    }
  };

  auto swap_cols = [&](std::size_t c1, std::size_t c2) {
    if (c1 != c2) {
      for (std::size_t i = 0; i < n; ++i) {
        std::swap(s[i][c1], s[i][c2]);
      }
    }
  };

  auto scale_row_by_scalar = [&](std::size_t r, const E& unit) {
    if (field->IsOne(unit)) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!s[r][j].empty()) {
        s[r][j] = Ops::Scale(*field, s[r][j], unit);
      }
    }
    const E inv_unit = field->Inv(unit);
    field->ScalarMulArray(std::span<E>(w_cols[r]),
                          std::span<const E>(w_cols[r]), inv_unit);
  };

  // Row operation: row[dst] += factor(x) * row[src]
  // Inverse column operation on W: w_col[src] -= factor(A) * w_col[dst]
  auto add_row_mul = [&](std::size_t dst, std::size_t src, const Poly& factor) {
    if (factor.empty()) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!s[src][j].empty()) {
        s[dst][j] =
            Ops::Add(*field, s[dst][j], Ops::Mul(*field, factor, s[src][j]));
      }
    }
    const std::vector<E> delta = detail::EvalPolyMatVec(
        *field, a_data, n, factor, std::span<const E>(w_cols[dst]));
    field->SubArray(std::span<E>(w_cols[src]), std::span<const E>(w_cols[src]),
                    std::span<const E>(delta));
  };

  auto add_col_mul = [&](std::size_t dst, std::size_t src, const Poly& factor) {
    if (factor.empty()) {
      return;
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (!s[i][src].empty()) {
        s[i][dst] =
            Ops::Add(*field, s[i][dst], Ops::Mul(*field, factor, s[i][src]));
      }
    }
  };

  // Smith Normal Form of xI - A over F[x].
  for (std::size_t r = 0; r < n; ++r) {
    std::size_t best_r = r;
    std::size_t best_c = r;
    std::size_t best_size = static_cast<std::size_t>(-1);
    for (std::size_t i = r; i < n; ++i) {
      for (std::size_t j = r; j < n; ++j) {
        if (!s[i][j].empty() && s[i][j].size() < best_size) {
          best_size = s[i][j].size();
          best_r = i;
          best_c = j;
        }
      }
    }
    swap_rows(r, best_r);
    swap_cols(r, best_c);

    while (true) {
      while (true) {
        // Clear column r below row r.
        bool col_has_nonzero = false;
        std::size_t min_row = r;
        std::size_t min_col_size = static_cast<std::size_t>(-1);
        for (std::size_t i = r + 1; i < n; ++i) {
          if (!s[i][r].empty()) {
            auto [q, rem] = Ops::DivRem(*field, s[i][r], s[r][r]);
            (void)rem;
            add_row_mul(i, r, Ops::Neg(*field, q));
            if (!s[i][r].empty() && s[i][r].size() < min_col_size) {
              col_has_nonzero = true;
              min_row = i;
              min_col_size = s[i][r].size();
            }
          }
        }
        if (col_has_nonzero) {
          swap_rows(r, min_row);
          continue;
        }

        // Clear row r right of column r.
        bool row_has_nonzero = false;
        std::size_t min_col = r;
        std::size_t min_row_size = static_cast<std::size_t>(-1);
        for (std::size_t j = r + 1; j < n; ++j) {
          if (!s[r][j].empty()) {
            auto [q, rem] = Ops::DivRem(*field, s[r][j], s[r][r]);
            (void)rem;
            add_col_mul(j, r, Ops::Neg(*field, q));
            if (!s[r][j].empty() && s[r][j].size() < min_row_size) {
              row_has_nonzero = true;
              min_col = j;
              min_row_size = s[r][j].size();
            }
          }
        }
        if (row_has_nonzero) {
          swap_cols(r, min_col);
          continue;
        }

        break;
      }

      // Check divisibility of trailing block S[r+1..n-1][r+1..n-1].
      bool found_non_multiple = false;
      std::size_t bad_row = r;
      for (std::size_t i = r + 1; i < n && !found_non_multiple; ++i) {
        for (std::size_t j = r + 1; j < n; ++j) {
          if (!s[i][j].empty()) {
            if (!Ops::DivRem(*field, s[i][j], s[r][r]).second.empty()) {
              found_non_multiple = true;
              bad_row = i;
              break;
            }
          }
        }
      }

      if (!found_non_multiple) {
        break;
      }

      add_row_mul(r, bad_row, Ops::One(*field));
    }

    // Normalize S[r][r] to monic.
    if (!field->IsOne(s[r][r].back())) {
      scale_row_by_scalar(r, field->Inv(s[r][r].back()));
    }
  }

  // Extract non-trivial invariant factors (degree >= 1) and construct the
  // block-diagonal companion matrix RCF and similarity matrix P.
  using RawArray = typename GFMatrix<F>::raw_array_type;
  using ShapeType = typename GFMatrix<F>::shape_type;
  RawArray rcf_raw = RawArray::from_shape(ShapeType{n, n});
  RawArray p_raw = RawArray::from_shape(ShapeType{n, n});
  rcf_raw.fill(zero);
  p_raw.fill(zero);
  E* rcf_data = rcf_raw.data();
  E* p_data = p_raw.data();

  std::size_t offset = 0;
  std::vector<E> cur(n, zero);
  std::vector<E> next(n, zero);

  for (std::size_t k = 0; k < n; ++k) {
    const Poly& f = s[k][k];
    if (f.size() <= 1) {
      continue;  // Unit invariant factor (1).
    }
    const std::size_t deg = f.size() - 1;
    result.invariant_factors.push_back(
        PolynomialDense<F>::FromRaw(f, field, var));

    // Populate companion block C(f) at [offset .. offset + deg - 1].
    for (std::size_t i = 0; i + 1 < deg; ++i) {
      rcf_data[(offset + i + 1) * n + (offset + i)] = one;
    }
    for (std::size_t i = 0; i < deg; ++i) {
      rcf_data[(offset + i) * n + (offset + deg - 1)] = field->Neg(f[i]);
    }

    // Populate columns [offset .. offset + deg - 1] of P with
    // {v, Av, A^2 v, ..., A^{deg-1} v} where v = w_cols[k].
    cur = w_cols[k];
    for (std::size_t j = 0; j < deg; ++j) {
      for (std::size_t row = 0; row < n; ++row) {
        p_data[row * n + (offset + j)] = cur[row];
      }
      if (j + 1 < deg) {
        detail::RawMatVec(*field, a_data, n, std::span<const E>(cur),
                          std::span<E>(next));
        cur.swap(next);
      }
    }

    offset += deg;
  }

  result.rational_canonical_form =
      GFMatrix<F>::FromRaw(field, std::move(rcf_raw));
  result.similarity_matrix = GFMatrix<F>::FromRaw(field, std::move(p_raw));
  return result;
}

}  // namespace linalg
}  // namespace xg

#endif  // XGALOIS_LINALG_PID_FORMS_HPP_
