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
std::vector<uint8_t> NaiveMatMulGf2(const uint8_t *a, std::size_t m,
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

std::pair<std::vector<uint8_t>, std::size_t> NaiveRrefGf2(const uint8_t *src,
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
GFArray<GaloisFieldBinary> MakeInvertibleGf2Matrix(
    const std::shared_ptr<GaloisFieldBinary> &field, std::size_t n,
    std::mt19937_64 &rng) {
  auto l_mat = linalg::Eye(field, n);
  auto u_mat = linalg::Eye(field, n);
  std::uniform_int_distribution<unsigned int> bit_dist(0U, 1U);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      l_mat.Raw(i, j) = static_cast<uint8_t>(bit_dist(rng));
    }
    for (std::size_t j = i + 1; j < n; ++j) {
      u_mat.Raw(i, j) = static_cast<uint8_t>(bit_dist(rng));
    }
  }
  auto a = linalg::Dot(l_mat, u_mat);
  // Apply a few row swaps to exercise pivoting across word boundaries.
  for (std::size_t i = 0; i < n; ++i) {
    std::uniform_int_distribution<std::size_t> row_dist(0, n - 1);
    const std::size_t j = row_dist(rng);
    if (i != j) {
      for (std::size_t c = 0; c < n; ++c) {
        std::swap(a.Raw(i, c), a.Raw(j, c));
      }
    }
  }
  return a;
}

GFArray<GaloisFieldBinary> EvalPolyAtMatrixGf2(
    const PolynomialDense<GaloisFieldBinary> &p,
    const GFArray<GaloisFieldBinary> &a) {
  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  const int deg = p.Degree();
  auto eye_mat = linalg::Eye(field, n);
  if (deg < 0) {
    return linalg::Zeros(field, {n, n});
  }
  auto acc = eye_mat * p[static_cast<std::size_t>(deg)];
  for (std::size_t k = static_cast<std::size_t>(deg); k-- > 0;) {
    acc = linalg::Dot(acc, a) + (eye_mat * p[k]);
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
    auto a = linalg::RandomMatrix(field, {m, k}, rng);
    auto b = linalg::RandomMatrix(field, {k, n}, rng);

    auto c = linalg::Dot(a, b);
    ASSERT_EQ(c.Shape(0), m);
    ASSERT_EQ(c.Shape(1), n);

    auto ref = NaiveMatMulGf2(a.RawData().data(), m, k, b.RawData().data(), n);
    for (std::size_t i = 0; i < m * n; ++i) {
      ASSERT_EQ(c.RawData().data()[i], ref[i])
          << "Mismatch for shape (" << m << "," << k << "," << n << ") at "
          << i;
    }

    // 1D * 2D
    auto u = linalg::RandomMatrix(field, {m}, rng);
    auto u_a = linalg::Dot(u, a);
    ASSERT_EQ(u_a.Dimension(), std::size_t{1});
    ASSERT_EQ(u_a.Shape(0), k);
    auto ref_ua =
        NaiveMatMulGf2(u.RawData().data(), 1, m, a.RawData().data(), k);
    for (std::size_t j = 0; j < k; ++j) {
      ASSERT_EQ(u_a.RawData().data()[j], ref_ua[j]);
    }

    // 2D * 1D
    auto v = linalg::RandomMatrix(field, {k}, rng);
    auto a_v = linalg::Dot(a, v);
    ASSERT_EQ(a_v.Dimension(), std::size_t{1});
    ASSERT_EQ(a_v.Shape(0), m);
    auto ref_av =
        NaiveMatMulGf2(a.RawData().data(), m, k, v.RawData().data(), 1);
    for (std::size_t i = 0; i < m; ++i) {
      ASSERT_EQ(a_v.RawData().data()[i], ref_av[i]);
    }

    // 1D * 1D
    auto v2 = linalg::RandomMatrix(field, {k}, rng);
    auto d = linalg::Dot(v, v2);
    auto ref_d =
        NaiveMatMulGf2(v.RawData().data(), 1, k, v2.RawData().data(), 1);
    EXPECT_EQ(d.RawData()(), ref_d[0]);
    EXPECT_EQ(linalg::Vdot(v, v2).Value(), ref_d[0]);
  }
}

TEST_F(LinalgGF2Test, RrefRankAndDetDifferentialAcrossWordBoundaries) {
  std::mt19937_64 rng(777);
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {1, 1}, {7, 7}, {64, 64}, {65, 65}, {130, 97}, {97, 130},
  };

  for (const auto &[rows, cols] : shapes) {
    auto a = linalg::RandomMatrix(field, {rows, cols}, rng);
    auto r = linalg::Rref(a);
    const std::size_t rank = linalg::MatrixRank(a);

    auto [ref_rref, ref_rank] = NaiveRrefGf2(a.RawData().data(), rows, cols);
    EXPECT_EQ(rank, ref_rank) << "Rank mismatch for " << rows << "x" << cols;
    for (std::size_t i = 0; i < rows * cols; ++i) {
      ASSERT_EQ(r.RawData().data()[i], ref_rref[i])
          << "RREF mismatch for " << rows << "x" << cols << " at " << i;
    }

    // Row echelon form must preserve rank and be upper trapezoidal.
    auto e = linalg::RowEchelon(a);
    EXPECT_EQ(linalg::MatrixRank(e), ref_rank);

    if (rows == cols) {
      const uint8_t expected_det = (ref_rank == rows) ? 1 : 0;
      EXPECT_EQ(linalg::Det(a).Value(), expected_det);
    }
  }
}

TEST_F(LinalgGF2Test, InversionAndSolveAcrossWordBoundaries) {
  std::mt19937_64 rng(999);
  for (std::size_t n : {std::size_t{1}, std::size_t{7}, std::size_t{64},
                        std::size_t{65}, std::size_t{100}}) {
    auto a = MakeInvertibleGf2Matrix(field, n, rng);
    EXPECT_EQ(linalg::Det(a).Value(), uint8_t{1});
    EXPECT_EQ(linalg::MatrixRank(a), n);

    auto a_inv = linalg::Inv(a);
    auto eye_mat = linalg::Eye(field, n);
    EXPECT_EQ(linalg::Dot(a, a_inv), eye_mat);
    EXPECT_EQ(linalg::Dot(a_inv, a), eye_mat);

    // 1D solve
    auto b1 = linalg::RandomMatrix(field, {n}, rng);
    auto x1 = linalg::Solve(a, b1);
    EXPECT_EQ(linalg::Dot(a, x1), b1);

    // Multi-RHS 2D solve (width crossing 64-bit word boundary: 70 columns)
    const std::size_t rhs_cols = (n >= 64) ? 70 : 5;
    auto b2 = linalg::RandomMatrix(field, {n, rhs_cols}, rng);
    auto x2 = linalg::Solve(a, b2);
    EXPECT_EQ(linalg::Dot(a, x2), b2);
  }

  // Singular matrix must throw for Inv and Solve.
  auto s = linalg::Zeros(field, {65, 65});
  s.Raw(0, 0) = 1;
  EXPECT_EQ(linalg::Det(s).Value(), uint8_t{0});
  EXPECT_THROW(linalg::Inv(s), std::runtime_error);
  EXPECT_THROW(linalg::Solve(s, linalg::Ones(field, {65})), std::runtime_error);
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
    auto a = linalg::RandomMatrix(field, {m, n}, rng);
    const std::size_t rank = linalg::MatrixRank(a);

    // Right nullspace: A * N^T == 0 and rank(A) + nullity == n
    auto ns = linalg::Nullspace(a);
    ASSERT_EQ(ns.Shape(0), n - rank);
    ASSERT_EQ(ns.Shape(1), n);
    EXPECT_EQ(linalg::MatrixRank(ns), n - rank);
    EXPECT_EQ(linalg::Dot(a, linalg::Transpose(ns)),
              linalg::Zeros(field, {m, n - rank}));

    // Left nullspace: L * A == 0
    auto lns = linalg::LeftNullspace(a);
    ASSERT_EQ(lns.Shape(0), m - rank);
    ASSERT_EQ(lns.Shape(1), m);
    EXPECT_EQ(linalg::Dot(lns, a), linalg::Zeros(field, {m - rank, n}));

    // Row space & column space
    auto rs = linalg::RowSpace(a);
    EXPECT_EQ(rs.Shape(0), rank);
    EXPECT_EQ(linalg::MatrixRank(rs), rank);

    auto cs = linalg::ColumnSpace(a);
    EXPECT_EQ(cs.Shape(1), rank);
    EXPECT_EQ(linalg::MatrixRank(cs), rank);

    // Systematic form: [I_k | P] and orthogonal dual H with A * H^T == 0
    auto sys = linalg::SystematicForm(a);
    const std::size_t k = sys.pivot_cols.size();
    EXPECT_EQ(k, rank);
    EXPECT_EQ(sys.free_cols.size(), n - k);
    EXPECT_EQ(sys.col_perm.size(), n);
    ASSERT_EQ(sys.g_sys.Shape(0), k);
    ASSERT_EQ(sys.g_sys.Shape(1), n);
    for (std::size_t i = 0; i < k; ++i) {
      for (std::size_t j = 0; j < k; ++j) {
        ASSERT_EQ(sys.g_sys.Raw(i, j), i == j ? uint8_t{1} : uint8_t{0});
      }
    }

    auto h = linalg::Zeros(field, {n - k, n});
    for (std::size_t r = 0; r < n - k; ++r) {
      h.Raw(r, sys.free_cols[r]) = 1;
      for (std::size_t c = 0; c < k; ++c) {
        h.Raw(r, sys.pivot_cols[c]) = sys.g_sys.Raw(c, k + r);
      }
    }
    EXPECT_EQ(linalg::Dot(a, linalg::Transpose(h)),
              linalg::Zeros(field, {m, n - k}));
  }

  // Full-column-rank matrix has {0, n} nullspace.
  auto eye_65 = linalg::Eye(field, 65);
  auto ns_empty = linalg::Nullspace(eye_65);
  EXPECT_EQ(ns_empty.Shape(0), std::size_t{0});
  EXPECT_EQ(ns_empty.Shape(1), std::size_t{65});
}

TEST_F(LinalgGF2Test, PluDecompositionAcrossWordBoundaries) {
  std::mt19937_64 rng(27182);
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {1, 1}, {7, 7}, {64, 64}, {65, 65}, {130, 97},
  };

  for (const auto &[m, n] : shapes) {
    auto a = linalg::RandomMatrix(field, {m, n}, rng);
    auto res = linalg::Plu(a);

    EXPECT_EQ(linalg::Dot(res.p, a), linalg::Dot(res.l, res.u))
        << "PLU failed for shape " << m << "x" << n;

    // L must be unit lower triangular, U must be upper trapezoidal.
    for (std::size_t i = 0; i < m; ++i) {
      ASSERT_EQ(res.l.Raw(i, i), uint8_t{1});
      for (std::size_t j = i + 1; j < m; ++j) {
        ASSERT_EQ(res.l.Raw(i, j), uint8_t{0});
      }
      const std::size_t u_zero_lim = std::min(i, n);
      for (std::size_t j = 0; j < u_zero_lim; ++j) {
        ASSERT_EQ(res.u.Raw(i, j), uint8_t{0});
      }
    }
  }
}

TEST_F(LinalgGF2Test, CharacteristicPolyCayleyHamiltonOverChar2) {
  std::mt19937_64 rng(8080);
  for (std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{4},
                        std::size_t{8}, std::size_t{12}}) {
    auto a = linalg::RandomMatrix(field, {n, n}, rng);
    auto p = linalg::CharacteristicPoly(a);

    EXPECT_EQ(p.Degree(), static_cast<int>(n));
    EXPECT_EQ(p[n].Value(), uint8_t{1});
    EXPECT_EQ(p[n - 1].Value(), linalg::Trace(a).Value());
    EXPECT_EQ(p[0].Value(), linalg::Det(a).Value());

    // Cayley-Hamilton: p_A(A) == 0 over GF(2)!
    EXPECT_EQ(EvalPolyAtMatrixGf2(p, a), linalg::Zeros(field, {n, n}));
  }
}

TEST_F(LinalgGF2Test, MatrixPower) {
  std::mt19937_64 rng(5150);
  auto a = MakeInvertibleGf2Matrix(field, 10, rng);

  EXPECT_EQ(linalg::MatrixPower(a, 0), linalg::Eye(field, 10));
  EXPECT_EQ(linalg::MatrixPower(a, 1), a);
  EXPECT_EQ(linalg::MatrixPower(a, -1), linalg::Inv(a));
  EXPECT_EQ(linalg::Dot(linalg::MatrixPower(a, 7), linalg::MatrixPower(a, -7)),
            linalg::Eye(field, 10));
}
