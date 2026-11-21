#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

namespace {

// Independent naive byte-by-byte GF(2) reference helpers for differential
// testing against the bit-packed uint64_t M4RI engine.
std::vector<uint8_t> NaiveMatMulGF2(const uint8_t *a, std::size_t m,
                                    std::size_t k, const uint8_t *b,
                                    std::size_t n) {
  std::vector<uint8_t> c(m * n, 0);
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t kk = 0; kk < k; ++kk) {
      if ((a[i * k + kk] & 1U) == 0) {
        continue;
      }
      for (std::size_t j = 0; j < n; ++j) {
        c[i * n + j] ^= (b[kk * n + j] & 1U);
      }
    }
  }
  return c;
}

std::pair<std::vector<uint8_t>, std::size_t> NaiveRrefGF2(const uint8_t *src,
                                                          std::size_t rows,
                                                          std::size_t cols) {
  std::vector<uint8_t> r(src, src + rows * cols);
  std::size_t pivot_row = 0;
  for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
    std::size_t found = pivot_row;
    while (found < rows && (r[found * cols + col] & 1U) == 0) {
      ++found;
    }
    if (found == rows) {
      continue;
    }
    if (found != pivot_row) {
      for (std::size_t j = 0; j < cols; ++j) {
        std::swap(r[pivot_row * cols + j], r[found * cols + j]);
      }
    }
    for (std::size_t k = 0; k < rows; ++k) {
      if (k != pivot_row && (r[k * cols + col] & 1U) != 0) {
        for (std::size_t j = col; j < cols; ++j) {
          r[k * cols + j] ^= r[pivot_row * cols + j];
        }
      }
    }
    ++pivot_row;
  }
  return {std::move(r), pivot_row};
}

// Construct a guaranteed-invertible n x n GF(2) matrix by multiplying random
// unit lower-triangular and unit upper-triangular matrices (with random row
// permutations).
GFArray<GaloisFieldBinary> MakeInvertibleGF2Matrix(
    const std::shared_ptr<GaloisFieldBinary> &field, std::size_t n,
    std::mt19937_64 &rng) {
  auto L = linalg::eye(field, n);
  auto U = linalg::eye(field, n);
  std::uniform_int_distribution<unsigned int> bit_dist(0U, 1U);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      L.Raw(i, j) = static_cast<uint8_t>(bit_dist(rng));
    }
    for (std::size_t j = i + 1; j < n; ++j) {
      U.Raw(i, j) = static_cast<uint8_t>(bit_dist(rng));
    }
  }
  auto A = linalg::dot(L, U);
  // Apply a few row swaps to exercise pivoting across word boundaries.
  for (std::size_t i = 0; i < n; ++i) {
    std::uniform_int_distribution<std::size_t> row_dist(0, n - 1);
    const std::size_t j = row_dist(rng);
    if (i != j) {
      for (std::size_t c = 0; c < n; ++c) {
        std::swap(A.Raw(i, c), A.Raw(j, c));
      }
    }
  }
  return A;
}

GFArray<GaloisFieldBinary> EvalPolyAtMatrixGF2(
    const PolynomialDense<GaloisFieldBinary> &p,
    const GFArray<GaloisFieldBinary> &A) {
  auto field = A.Field();
  const std::size_t n = A.Shape(0);
  const int deg = p.Degree();
  auto I = linalg::eye(field, n);
  if (deg < 0) {
    return linalg::zeros(field, {n, n});
  }
  auto acc = I * p[static_cast<std::size_t>(deg)];
  for (std::size_t k = static_cast<std::size_t>(deg); k-- > 0;) {
    acc = linalg::dot(acc, A) + (I * p[k]);
  }
  return acc;
}

}  // namespace

class LinalgGF2Test : public ::testing::Test {
 protected:
  void SetUp() override { field = std::make_shared<GaloisFieldBinary>(); }
  std::shared_ptr<GaloisFieldBinary> field;
};

TEST_F(LinalgGF2Test, DotProductDifferentialAcrossWordBoundaries) {
  std::mt19937_64 rng(2026);
  const std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> shapes =
      {
          {1, 1, 1},
          {7, 7, 7},
          {64, 64, 64},
          {65, 65, 65},
          {130, 97, 73},
          {3, 65, 130},  // Exercises m < 4 direct bit-scan fast
                         // path
          {65, 130, 3},
      };

  for (const auto &[m, k, n] : shapes) {
    auto A = linalg::random_matrix(field, {m, k}, rng);
    auto B = linalg::random_matrix(field, {k, n}, rng);

    auto C = linalg::dot(A, B);
    ASSERT_EQ(C.Shape(0), m);
    ASSERT_EQ(C.Shape(1), n);

    auto ref = NaiveMatMulGF2(A.RawData().data(), m, k, B.RawData().data(), n);
    for (std::size_t i = 0; i < m * n; ++i) {
      ASSERT_EQ(C.RawData().data()[i], ref[i])
          << "Mismatch for shape (" << m << "," << k << "," << n << ") at "
          << i;
    }

    // 1D * 2D
    auto u = linalg::random_matrix(field, {m}, rng);
    auto uA = linalg::dot(u, A);
    ASSERT_EQ(uA.Dimension(), 1);
    ASSERT_EQ(uA.Shape(0), k);
    auto ref_uA =
        NaiveMatMulGF2(u.RawData().data(), 1, m, A.RawData().data(), k);
    for (std::size_t j = 0; j < k; ++j) {
      ASSERT_EQ(uA.RawData().data()[j], ref_uA[j]);
    }

    // 2D * 1D
    auto v = linalg::random_matrix(field, {k}, rng);
    auto Av = linalg::dot(A, v);
    ASSERT_EQ(Av.Dimension(), 1);
    ASSERT_EQ(Av.Shape(0), m);
    auto ref_Av =
        NaiveMatMulGF2(A.RawData().data(), m, k, v.RawData().data(), 1);
    for (std::size_t i = 0; i < m; ++i) {
      ASSERT_EQ(Av.RawData().data()[i], ref_Av[i]);
    }

    // 1D * 1D
    auto v2 = linalg::random_matrix(field, {k}, rng);
    auto d = linalg::dot(v, v2);
    auto ref_d =
        NaiveMatMulGF2(v.RawData().data(), 1, k, v2.RawData().data(), 1);
    EXPECT_EQ(d.RawData()(), ref_d[0]);
    EXPECT_EQ(linalg::vdot(v, v2).Value(), ref_d[0]);
  }
}

TEST_F(LinalgGF2Test, RrefRankAndDetDifferentialAcrossWordBoundaries) {
  std::mt19937_64 rng(777);
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {1, 1}, {7, 7}, {64, 64}, {65, 65}, {130, 97}, {97, 130},
  };

  for (const auto &[rows, cols] : shapes) {
    auto A = linalg::random_matrix(field, {rows, cols}, rng);
    auto R = linalg::rref(A);
    const std::size_t rank = linalg::matrix_rank(A);

    auto [ref_rref, ref_rank] = NaiveRrefGF2(A.RawData().data(), rows, cols);
    EXPECT_EQ(rank, ref_rank) << "Rank mismatch for " << rows << "x" << cols;
    for (std::size_t i = 0; i < rows * cols; ++i) {
      ASSERT_EQ(R.RawData().data()[i], ref_rref[i])
          << "RREF mismatch for " << rows << "x" << cols << " at " << i;
    }

    // Row echelon form must preserve rank and be upper trapezoidal.
    auto E = linalg::row_echelon(A);
    EXPECT_EQ(linalg::matrix_rank(E), ref_rank);

    if (rows == cols) {
      const uint8_t expected_det = (ref_rank == rows) ? 1 : 0;
      EXPECT_EQ(linalg::det(A).Value(), expected_det);
    }
  }
}

TEST_F(LinalgGF2Test, InversionAndSolveAcrossWordBoundaries) {
  std::mt19937_64 rng(999);
  for (std::size_t n : {std::size_t{1}, std::size_t{7}, std::size_t{64},
                        std::size_t{65}, std::size_t{100}}) {
    auto A = MakeInvertibleGF2Matrix(field, n, rng);
    EXPECT_EQ(linalg::det(A).Value(), 1);
    EXPECT_EQ(linalg::matrix_rank(A), n);

    auto A_inv = linalg::inv(A);
    auto I = linalg::eye(field, n);
    EXPECT_EQ(linalg::dot(A, A_inv), I);
    EXPECT_EQ(linalg::dot(A_inv, A), I);

    // 1D solve
    auto b1 = linalg::random_matrix(field, {n}, rng);
    auto x1 = linalg::solve(A, b1);
    EXPECT_EQ(linalg::dot(A, x1), b1);

    // Multi-RHS 2D solve (width crossing 64-bit word boundary: 70 columns)
    const std::size_t rhs_cols = (n >= 64) ? 70 : 5;
    auto B2 = linalg::random_matrix(field, {n, rhs_cols}, rng);
    auto X2 = linalg::solve(A, B2);
    EXPECT_EQ(linalg::dot(A, X2), B2);
  }

  // Singular matrix must throw for inv and solve.
  auto S = linalg::zeros(field, {65, 65});
  S.Raw(0, 0) = 1;
  EXPECT_EQ(linalg::det(S).Value(), 0);
  EXPECT_THROW(linalg::inv(S), std::runtime_error);
  EXPECT_THROW(linalg::solve(S, linalg::ones(field, {65})), std::runtime_error);
}

TEST_F(LinalgGF2Test, NullspaceSystematicFormAndSubspaces) {
  std::mt19937_64 rng(31415);
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {7, 12},
      {64, 64},
      {65, 80},
      {130, 97},
  };

  for (const auto &[m, n] : shapes) {
    auto A = linalg::random_matrix(field, {m, n}, rng);
    const std::size_t rank = linalg::matrix_rank(A);

    // Right nullspace: A * N^T == 0 and rank(A) + nullity == n
    auto N = linalg::nullspace(A);
    ASSERT_EQ(N.Shape(0), n - rank);
    ASSERT_EQ(N.Shape(1), n);
    EXPECT_EQ(linalg::matrix_rank(N), n - rank);
    EXPECT_EQ(linalg::dot(A, linalg::transpose(N)),
              linalg::zeros(field, {m, n - rank}));

    // Left nullspace: L * A == 0
    auto L = linalg::left_nullspace(A);
    ASSERT_EQ(L.Shape(0), m - rank);
    ASSERT_EQ(L.Shape(1), m);
    EXPECT_EQ(linalg::dot(L, A), linalg::zeros(field, {m - rank, n}));

    // Row space & column space
    auto RS = linalg::row_space(A);
    EXPECT_EQ(RS.Shape(0), rank);
    EXPECT_EQ(linalg::matrix_rank(RS), rank);

    auto CS = linalg::column_space(A);
    EXPECT_EQ(CS.Shape(1), rank);
    EXPECT_EQ(linalg::matrix_rank(CS), rank);

    // Systematic form: [I_k | P] and orthogonal dual H with A * H^T == 0
    auto sys = linalg::systematic_form(A);
    const std::size_t k = sys.pivot_cols.size();
    EXPECT_EQ(k, rank);
    EXPECT_EQ(sys.free_cols.size(), n - k);
    EXPECT_EQ(sys.col_perm.size(), n);
    ASSERT_EQ(sys.G_sys.Shape(0), k);
    ASSERT_EQ(sys.G_sys.Shape(1), n);
    for (std::size_t i = 0; i < k; ++i) {
      for (std::size_t j = 0; j < k; ++j) {
        ASSERT_EQ(sys.G_sys.Raw(i, j), i == j ? 1U : 0U);
      }
    }

    auto H = linalg::zeros(field, {n - k, n});
    for (std::size_t r = 0; r < n - k; ++r) {
      H.Raw(r, sys.free_cols[r]) = 1;
      for (std::size_t c = 0; c < k; ++c) {
        H.Raw(r, sys.pivot_cols[c]) = sys.G_sys.Raw(c, k + r);
      }
    }
    EXPECT_EQ(linalg::dot(A, linalg::transpose(H)),
              linalg::zeros(field, {m, n - k}));
  }

  // Full-column-rank matrix has {0, n} nullspace.
  auto I65 = linalg::eye(field, 65);
  auto N0 = linalg::nullspace(I65);
  EXPECT_EQ(N0.Shape(0), 0);
  EXPECT_EQ(N0.Shape(1), 65);
}

TEST_F(LinalgGF2Test, PLUDecompositionAcrossWordBoundaries) {
  std::mt19937_64 rng(27182);
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {1, 1}, {7, 7}, {64, 64}, {65, 65}, {130, 97},
  };

  for (const auto &[m, n] : shapes) {
    auto A = linalg::random_matrix(field, {m, n}, rng);
    auto res = linalg::plu(A);

    EXPECT_EQ(linalg::dot(res.P, A), linalg::dot(res.L, res.U))
        << "PLU failed for shape " << m << "x" << n;

    // L must be unit lower triangular, U must be upper trapezoidal.
    for (std::size_t i = 0; i < m; ++i) {
      ASSERT_EQ(res.L.Raw(i, i), 1);
      for (std::size_t j = i + 1; j < m; ++j) {
        ASSERT_EQ(res.L.Raw(i, j), 0);
      }
      const std::size_t u_zero_lim = std::min(i, n);
      for (std::size_t j = 0; j < u_zero_lim; ++j) {
        ASSERT_EQ(res.U.Raw(i, j), 0);
      }
    }
  }
}

TEST_F(LinalgGF2Test, CharacteristicPolyCayleyHamiltonOverChar2) {
  std::mt19937_64 rng(8080);
  for (std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{4},
                        std::size_t{8}, std::size_t{12}}) {
    auto A = linalg::random_matrix(field, {n, n}, rng);
    auto p = linalg::characteristic_poly(A);

    EXPECT_EQ(p.Degree(), static_cast<int>(n));
    EXPECT_EQ(p[n].Value(), 1);
    EXPECT_EQ(p[n - 1].Value(), linalg::trace(A).Value());
    EXPECT_EQ(p[0].Value(), linalg::det(A).Value());

    // Cayley-Hamilton: p_A(A) == 0 over GF(2)!
    EXPECT_EQ(EvalPolyAtMatrixGF2(p, A), linalg::zeros(field, {n, n}));
  }
}

TEST_F(LinalgGF2Test, MatrixPowerAndGArrayCompatibility) {
  std::mt19937_64 rng(5150);
  auto A = MakeInvertibleGF2Matrix(field, 10, rng);

  EXPECT_EQ(linalg::matrix_power(A, 0), linalg::eye(field, 10));
  EXPECT_EQ(linalg::matrix_power(A, 1), A);
  EXPECT_EQ(linalg::matrix_power(A, -1), linalg::inv(A));
  EXPECT_EQ(
      linalg::dot(linalg::matrix_power(A, 7), linalg::matrix_power(A, -7)),
      linalg::eye(field, 10));

  // Verify legacy garray API routes cleanly through bit-packed GF(2) engine.
  auto A_ga = A.ToElementArray();
  auto A_inv_ga = linalg::inv(A_ga);
  auto I_ga = linalg::dot(A_ga, A_inv_ga);
  for (std::size_t i = 0; i < 10; ++i) {
    for (std::size_t j = 0; j < 10; ++j) {
      EXPECT_EQ(I_ga(i, j).Value(), i == j ? 1U : 0U);
    }
  }
}
