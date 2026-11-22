#ifndef XGALOIS_ALGEBRA_SNF_HPP_
#define XGALOIS_ALGEBRA_SNF_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace xg {
namespace algebra {

namespace detail {

template <typename IntT>
inline IntT SnfAbs(const IntT& x) {
  return x < IntT{0} ? -x : x;
}

// Nearest-integer quotient q such that |a - q * b| <= |b| / 2 (for b > 0),
// minimizing intermediate entry growth without overflowing signed integers.
template <typename IntT>
inline IntT SnfNearestQuotient(const IntT& a, const IntT& b) {
  IntT q = a / b;
  IntT r = a - q * b;
  if (r > IntT{0}) {
    if (r > b - r) {
      q = q + IntT{1};
    }
  } else if (r < IntT{0}) {
    const IntT nr = -r;
    if (nr > b - nr) {
      q = q - IntT{1};
    }
  }
  return q;
}

template <typename IntT>
inline std::vector<std::vector<IntT>> MakeIdentityMatrix(std::size_t n) {
  std::vector<std::vector<IntT>> eye(n, std::vector<IntT>(n, IntT{0}));
  for (std::size_t i = 0; i < n; ++i) {
    eye[i][i] = IntT{1};
  }
  return eye;
}

template <typename IntT>
inline void FactorPositiveInteger(IntT n,
                                  std::vector<std::pair<IntT, uint64_t>>& out) {
  if (n <= IntT{1}) {
    return;
  }

  auto extract_factor = [&](const IntT& p) {
    if (n % p == IntT{0}) {
      uint64_t exp = 0;
      while (n % p == IntT{0}) {
        n = n / p;
        ++exp;
      }
      out.emplace_back(p, exp);
    }
  };

  extract_factor(IntT{2});
  extract_factor(IntT{3});

  // 6k +/- 1 wheel trial division; uses d <= n / d to avoid squaring overflow
  // across int64_t, __int128, and BigInt.
  IntT d{5};
  while (d <= n / d) {
    extract_factor(d);
    const IntT d2 = d + IntT{2};
    if (d2 <= n / d2) {
      extract_factor(d2);
    }
    d = d + IntT{6};
  }

  if (n > IntT{1}) {
    out.emplace_back(std::move(n), uint64_t{1});
  }
}

}  // namespace detail

// Result of Smith Normal Form decomposition over Z:
//   U * A * V == S
// where U (m x m) and V (n x n) are unimodular (|det(U)| == |det(V)| == 1),
// S (m x n) is diagonal with strictly positive invariant_factors
// d_1 | d_2 | ... | d_r on S[0][0], ..., S[r-1][r-1] and zeros elsewhere.
template <typename IntT>
struct SmithNormalFormResult {
  std::vector<std::vector<IntT>> S;
  std::vector<std::vector<IntT>> U;
  std::vector<std::vector<IntT>> V;
  std::vector<IntT> invariant_factors;
  std::size_t rank = 0;
};

template <typename IntT>
SmithNormalFormResult<IntT> ComputeSmithNormalForm(
    const std::vector<std::vector<IntT>>& matrix) {
  SmithNormalFormResult<IntT> result;
  const std::size_t m = matrix.size();
  if (m == 0) {
    return result;
  }
  const std::size_t n = matrix[0].size();
  for (std::size_t i = 1; i < m; ++i) {
    if (matrix[i].size() != n) {
      throw std::invalid_argument(
          "ComputeSmithNormalForm: input matrix must be rectangular");
    }
  }

  result.S = matrix;
  result.U = detail::MakeIdentityMatrix<IntT>(m);
  result.V = detail::MakeIdentityMatrix<IntT>(n);

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

  auto negate_row = [&](std::size_t r) {
    for (std::size_t j = 0; j < n; ++j) {
      S[r][j] = -S[r][j];
    }
    for (std::size_t j = 0; j < m; ++j) {
      U[r][j] = -U[r][j];
    }
  };

  auto add_row_mul = [&](std::size_t dst, std::size_t src, const IntT& factor) {
    if (factor == IntT{0}) {
      return;
    }
    for (std::size_t j = 0; j < n; ++j) {
      S[dst][j] = S[dst][j] + factor * S[src][j];
    }
    for (std::size_t j = 0; j < m; ++j) {
      U[dst][j] = U[dst][j] + factor * U[src][j];
    }
  };

  auto add_col_mul = [&](std::size_t dst, std::size_t src, const IntT& factor) {
    if (factor == IntT{0}) {
      return;
    }
    for (std::size_t i = 0; i < m; ++i) {
      S[i][dst] = S[i][dst] + factor * S[i][src];
    }
    for (std::size_t i = 0; i < n; ++i) {
      V[i][dst] = V[i][dst] + factor * V[i][src];
    }
  };

  const std::size_t limit = std::min(m, n);
  std::size_t r = 0;

  while (r < limit) {
    // 1. Select the nonzero pivot in S[r..m-1][r..n-1] with minimal |value|.
    bool found_pivot = false;
    std::size_t best_r = r;
    std::size_t best_c = r;
    IntT best_abs{0};
    for (std::size_t i = r; i < m; ++i) {
      for (std::size_t j = r; j < n; ++j) {
        if (S[i][j] != IntT{0}) {
          const IntT mag = detail::SnfAbs(S[i][j]);
          if (!found_pivot || mag < best_abs) {
            found_pivot = true;
            best_r = i;
            best_c = j;
            best_abs = mag;
          }
        }
      }
    }

    if (!found_pivot) {
      break;
    }

    swap_rows(r, best_r);
    swap_cols(r, best_c);

    // 2. Reduce row r and column r to zero off-diagonal and ensure S[r][r]
    //    divides all entries in S[r+1..m-1][r+1..n-1].
    while (true) {
      if (S[r][r] < IntT{0}) {
        negate_row(r);
      }

      // Iteratively clear column r below r and row r right of r.
      while (true) {
        // Clear column r below r.
        bool col_has_nonzero = false;
        std::size_t min_row = r;
        IntT min_col_abs{0};
        for (std::size_t i = r + 1; i < m; ++i) {
          if (S[i][r] != IntT{0}) {
            const IntT q = detail::SnfNearestQuotient(S[i][r], S[r][r]);
            add_row_mul(i, r, -q);
            if (S[i][r] != IntT{0}) {
              const IntT mag = detail::SnfAbs(S[i][r]);
              if (!col_has_nonzero || mag < min_col_abs) {
                col_has_nonzero = true;
                min_row = i;
                min_col_abs = mag;
              }
            }
          }
        }
        if (col_has_nonzero) {
          swap_rows(r, min_row);
          if (S[r][r] < IntT{0}) {
            negate_row(r);
          }
          continue;
        }

        // Clear row r right of r.
        bool row_has_nonzero = false;
        std::size_t min_col = r;
        IntT min_row_abs{0};
        for (std::size_t j = r + 1; j < n; ++j) {
          if (S[r][j] != IntT{0}) {
            const IntT q = detail::SnfNearestQuotient(S[r][j], S[r][r]);
            add_col_mul(j, r, -q);
            if (S[r][j] != IntT{0}) {
              const IntT mag = detail::SnfAbs(S[r][j]);
              if (!row_has_nonzero || mag < min_row_abs) {
                row_has_nonzero = true;
                min_col = j;
                min_row_abs = mag;
              }
            }
          }
        }
        if (row_has_nonzero) {
          swap_cols(r, min_col);
          if (S[r][r] < IntT{0}) {
            negate_row(r);
          }
          continue;
        }

        break;
      }

      if (S[r][r] < IntT{0}) {
        negate_row(r);
      }

      // Verify divisibility of the trailing submatrix S[r+1..m-1][r+1..n-1].
      bool found_non_multiple = false;
      std::size_t bad_row = r;
      for (std::size_t i = r + 1; i < m && !found_non_multiple; ++i) {
        for (std::size_t j = r + 1; j < n; ++j) {
          if (S[i][j] % S[r][r] != IntT{0}) {
            found_non_multiple = true;
            bad_row = i;
            break;
          }
        }
      }

      if (!found_non_multiple) {
        break;
      }

      // Add the row containing the non-multiple to row r; the next Euclidean
      // pass strictly decreases |S[r][r]|.
      add_row_mul(r, bad_row, IntT{1});
    }

    result.invariant_factors.push_back(S[r][r]);
    ++r;
  }

  result.rank = r;
  return result;
}

// Decomposition of a finitely generated abelian group
//   G = Z^n / rowspace(R)
// into:
//   G ~= Z^(free_rank) (+) (Z/d_1 Z) (+) ... (+) (Z/d_k Z)
// where each invariant factor satisfies d_i > 1 and d_1 | d_2 | ... | d_k,
// and elementary_divisors lists the corresponding prime powers (p, e) for p^e.
template <typename IntT>
struct AbelianGroupDecompositionResult {
  std::size_t free_rank = 0;
  std::vector<IntT> invariant_factors;
  std::vector<std::pair<IntT, uint64_t>> elementary_divisors;
};

// Decomposes G = Z^n / rowspace(relation_matrix), where each row of
// relation_matrix represents a relation on the n generators (n = column count).
template <typename IntT>
AbelianGroupDecompositionResult<IntT> DecomposeFinitelyGeneratedAbelianGroup(
    const std::vector<std::vector<IntT>>& relation_matrix) {
  AbelianGroupDecompositionResult<IntT> result;
  if (relation_matrix.empty()) {
    return result;
  }

  const std::size_t num_generators = relation_matrix[0].size();
  const auto snf = ComputeSmithNormalForm(relation_matrix);

  result.free_rank = num_generators - snf.rank;
  for (const IntT& d : snf.invariant_factors) {
    if (d > IntT{1}) {
      result.invariant_factors.push_back(d);
      detail::FactorPositiveInteger(d, result.elementary_divisors);
    }
  }

  std::sort(result.elementary_divisors.begin(),
            result.elementary_divisors.end(),
            [](const auto& lhs, const auto& rhs) {
              if (lhs.first != rhs.first) {
                return lhs.first < rhs.first;
              }
              return lhs.second < rhs.second;
            });

  return result;
}

}  // namespace algebra

using algebra::AbelianGroupDecompositionResult;
using algebra::ComputeSmithNormalForm;
using algebra::DecomposeFinitelyGeneratedAbelianGroup;
using algebra::SmithNormalFormResult;

}  // namespace xg

#endif  // XGALOIS_ALGEBRA_SNF_HPP_
