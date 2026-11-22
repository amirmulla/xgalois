#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <tuple>
#include <vector>

#include "xgalois/coding/ldpc.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/fast_matmul.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/linalg/pid_forms.hpp"
#include "xgalois/linalg/sparse_linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/ring/integer_ring.hpp"
#include "xgalois/ring/polynomial_ring.hpp"
#include "xgalois/utils/bigint.hpp"

using namespace xg;

namespace {

template <EuclideanDomain R>
std::vector<std::vector<typename R::element_type>> RingMatMul(
    const R& ring, const std::vector<std::vector<typename R::element_type>>& a,
    const std::vector<std::vector<typename R::element_type>>& b) {
  using E = typename R::element_type;
  const std::size_t m = a.size();
  const std::size_t k = b.size();
  const std::size_t n = (k == 0) ? 0 : b[0].size();
  std::vector<std::vector<E>> c(m, std::vector<E>(n, ring.AdditiveIdentity()));
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t t = 0; t < k; ++t) {
      if (ring.IsZero(a[i][t])) {
        continue;
      }
      for (std::size_t j = 0; j < n; ++j) {
        if (!ring.IsZero(b[t][j])) {
          c[i][j] = ring.Add(c[i][j], ring.Mul(a[i][t], b[t][j]));
        }
      }
    }
  }
  return c;
}

template <EuclideanDomain R>
void ExpectRingMatEq(
    const R& ring, const std::vector<std::vector<typename R::element_type>>& a,
    const std::vector<std::vector<typename R::element_type>>& b) {
  ASSERT_EQ(a.size(), b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    ASSERT_EQ(a[i].size(), b[i].size());
    for (std::size_t j = 0; j < a[i].size(); ++j) {
      EXPECT_TRUE(alg::Equals(ring, a[i][j], b[i][j]))
          << "Mismatch at (" << i << ", " << j
          << "): " << ring.ToString(a[i][j]) << " vs "
          << ring.ToString(b[i][j]);
    }
  }
}

template <EuclideanDomain R>
void ExpectUnimodular(
    const R& ring,
    const std::vector<std::vector<typename R::element_type>>& u) {
  const auto snf = linalg::SmithNormalForm(ring, u);
  ASSERT_EQ(snf.rank, u.size());
  for (const auto& d : snf.invariant_factors) {
    EXPECT_TRUE(alg::IsUnit(ring, d))
        << "Non-unit invariant factor in unimodular matrix: "
        << ring.ToString(d);
  }
}

template <FiniteField F>
GFMatrix<F> MakeGFMatrix(const std::shared_ptr<F>& field, std::size_t rows,
                         std::size_t cols,
                         const std::vector<typename F::element_type>& vals) {
  typename GFMatrix<F>::raw_array_type raw =
      GFMatrix<F>::raw_array_type::from_shape({rows, cols});
  for (std::size_t i = 0; i < rows * cols; ++i) {
    raw.data()[i] = vals[i];
  }
  return GFMatrix<F>(field, raw);
}

template <FiniteField F>
GFMatrix<F> EvalPolyAtMatrix(const PolynomialDense<F>& p,
                             const GFMatrix<F>& a) {
  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  const int deg = p.Degree();
  auto eye = linalg::Eye(field, n);
  if (deg < 0) {
    return linalg::Zeros(field, {n, n});
  }
  auto acc = eye * p[static_cast<std::size_t>(deg)];
  for (std::size_t k = static_cast<std::size_t>(deg); k-- > 0;) {
    acc = linalg::Dot(acc, a) + (eye * p[k]);
  }
  return acc;
}

template <FiniteField F>
void VerifyFrobeniusAndPolys(const GFMatrix<F>& a) {
  auto field = a.Field();
  const std::size_t n = a.Shape(0);

  const auto chi = linalg::CharacteristicPoly(a);
  const auto mu = linalg::MinimalPoly(a);
  const auto frob = linalg::FrobeniusNormalForm(a);

  // 1. Cayley-Hamilton: chi(A) == 0 and mu(A) == 0.
  const auto zero_mat = linalg::Zeros(field, {n, n});
  EXPECT_EQ(EvalPolyAtMatrix(chi, a), zero_mat);
  EXPECT_EQ(EvalPolyAtMatrix(mu, a), zero_mat);
  EXPECT_TRUE(chi.IsMonic());
  EXPECT_TRUE(mu.IsMonic());
  EXPECT_EQ(chi.Degree(), static_cast<int>(n));

  // 2. mu | chi.
  EXPECT_TRUE((chi % mu).IsZero());

  // 3. Invariant factors: monic, f_i | f_{i+1}, prod(f_i) == chi, f_s == mu.
  auto prod = PolynomialDense<F>::One(field);
  for (std::size_t i = 0; i < frob.invariant_factors.size(); ++i) {
    const auto& fi = frob.invariant_factors[i];
    EXPECT_TRUE(fi.IsMonic());
    EXPECT_GE(fi.Degree(), 1);
    prod = prod * fi;
    if (i + 1 < frob.invariant_factors.size()) {
      EXPECT_TRUE((frob.invariant_factors[i + 1] % fi).IsZero());
    }
  }
  EXPECT_EQ(prod, chi);
  if (n > 0) {
    ASSERT_FALSE(frob.invariant_factors.empty());
    EXPECT_EQ(frob.invariant_factors.back(), mu);
  }

  // 4. Similarity: P invertible and P^{-1} * A * P == RCF.
  const auto det_p = linalg::Det(frob.similarity_matrix);
  EXPECT_FALSE(field->IsZero(det_p.RawValue()));
  const auto p_inv = linalg::Inv(frob.similarity_matrix);
  const auto transformed =
      linalg::Dot(linalg::Dot(p_inv, a), frob.similarity_matrix);
  EXPECT_EQ(transformed, frob.rational_canonical_form);
}

}  // namespace

// ============================================================================
// 1. HermiteNormalForm & SmithNormalForm over Z (IntegerRing<BigInt>)
// ============================================================================

TEST(PidFormsTest, HermiteAndSmithOverBigInt) {
  const ring::IntegerRing<BigInt> z_ring;
  std::vector<std::vector<BigInt>> a = {
      {BigInt(2), BigInt(4), BigInt(4)},
      {BigInt(-6), BigInt(6), BigInt(12)},
      {BigInt(10), BigInt(-4), BigInt(-16)},
  };

  const auto hnf = linalg::HermiteNormalForm(z_ring, a);
  EXPECT_EQ(hnf.rank, std::size_t{3});
  EXPECT_EQ(hnf.pivot_cols, (std::vector<std::size_t>{0, 1, 2}));
  ExpectRingMatEq(z_ring, RingMatMul(z_ring, hnf.U, a), hnf.H);
  ExpectUnimodular(z_ring, hnf.U);

  // Verify HNF canonical pivots and reduced upper-triangular entries.
  for (std::size_t r = 0; r < hnf.rank; ++r) {
    const std::size_t pc = hnf.pivot_cols[r];
    EXPECT_GT(hnf.H[r][pc], BigInt(0));
    for (std::size_t i = 0; i < r; ++i) {
      EXPECT_GE(hnf.H[i][pc], BigInt(0));
      EXPECT_LT(hnf.H[i][pc], hnf.H[r][pc]);
    }
    for (std::size_t i = r + 1; i < a.size(); ++i) {
      EXPECT_EQ(hnf.H[i][pc], BigInt(0));
    }
  }

  const auto snf = linalg::SmithNormalForm(z_ring, a);
  EXPECT_EQ(snf.rank, std::size_t{3});
  ExpectRingMatEq(
      z_ring, RingMatMul(z_ring, RingMatMul(z_ring, snf.U, a), snf.V), snf.S);
  ExpectUnimodular(z_ring, snf.U);
  ExpectUnimodular(z_ring, snf.V);

  // Invariant factors for this classic matrix are [2, 6, 12].
  ASSERT_EQ(snf.invariant_factors.size(), std::size_t{3});
  EXPECT_EQ(snf.invariant_factors[0], BigInt(2));
  EXPECT_EQ(snf.invariant_factors[1], BigInt(6));
  EXPECT_EQ(snf.invariant_factors[2], BigInt(12));
  for (std::size_t i = 0; i + 1 < snf.invariant_factors.size(); ++i) {
    EXPECT_TRUE(z_ring.IsZero(alg::DivRem(z_ring, snf.invariant_factors[i + 1],
                                          snf.invariant_factors[i])
                                  .second));
  }
}

TEST(PidFormsTest, RectangularRankDeficientOverBigInt) {
  const ring::IntegerRing<BigInt> z_ring;
  std::vector<std::vector<BigInt>> a = {
      {BigInt(6), BigInt(10), BigInt(14), BigInt(2)},
      {BigInt(9), BigInt(15), BigInt(21), BigInt(3)},
      {BigInt(0), BigInt(4), BigInt(8), BigInt(12)},
  };

  const auto hnf = linalg::HermiteNormalForm(z_ring, a);
  EXPECT_EQ(hnf.rank, std::size_t{2});
  ExpectRingMatEq(z_ring, RingMatMul(z_ring, hnf.U, a), hnf.H);
  ExpectUnimodular(z_ring, hnf.U);

  const auto snf = linalg::SmithNormalForm(z_ring, a);
  EXPECT_EQ(snf.rank, std::size_t{2});
  ExpectRingMatEq(
      z_ring, RingMatMul(z_ring, RingMatMul(z_ring, snf.U, a), snf.V), snf.S);
  ExpectUnimodular(z_ring, snf.U);
  ExpectUnimodular(z_ring, snf.V);
  ASSERT_EQ(snf.invariant_factors.size(), std::size_t{2});
  EXPECT_TRUE(z_ring.IsZero(
      alg::DivRem(z_ring, snf.invariant_factors[1], snf.invariant_factors[0])
          .second));
}

// ============================================================================
// 2. HermiteNormalForm & SmithNormalForm over F_7[x]
// ============================================================================

TEST(PidFormsTest, HermiteAndSmithOverGF7PolyRing) {
  using GF7 = GaloisFieldPrime<uint32_t>;
  using PolyRing7 = ring::PolynomialRing<GF7>;
  using P = typename PolyRing7::element_type;

  const GF7 gf7(7);
  const PolyRing7 poly_ring(gf7);

  // Build 3x3 polynomial matrix M(x) = xI - A over F_7[x] for a matrix with
  // repeated eigenvalue 2 and non-trivial Jordan structure.
  const P x = poly_ring.X();
  const P c0 = poly_ring.AdditiveIdentity();
  const P c1 = poly_ring.MultiplicativeIdentity();
  const P c2 = poly_ring.Constant(2);
  const P c5 = poly_ring.Constant(5);  // -2 mod 7
  const P c6 = poly_ring.Constant(6);  // -1 mod 7

  // A = [[2, 1, 0], [0, 2, 0], [0, 0, 2]] -> xI - A has invariant factors
  // [1, (x - 2), (x - 2)^2].
  std::vector<std::vector<P>> mat = {
      {poly_ring.Add(x, c5), c6, c0},
      {c0, poly_ring.Add(x, c5), c0},
      {c0, c0, poly_ring.Add(x, c5)},
  };
  (void)c1;
  (void)c2;

  const auto hnf = linalg::HermiteNormalForm(poly_ring, mat);
  EXPECT_EQ(hnf.rank, std::size_t{3});
  ExpectRingMatEq(poly_ring, RingMatMul(poly_ring, hnf.U, mat), hnf.H);
  ExpectUnimodular(poly_ring, hnf.U);

  for (std::size_t r = 0; r < hnf.rank; ++r) {
    const std::size_t pc = hnf.pivot_cols[r];
    const auto [canon, unit] = alg::AssociateCanonical(poly_ring, hnf.H[r][pc]);
    (void)unit;
    EXPECT_TRUE(poly_ring.Equals(hnf.H[r][pc], canon));
    for (std::size_t i = 0; i < r; ++i) {
      EXPECT_LT(poly_ring.Degree(hnf.H[i][pc]), poly_ring.Degree(hnf.H[r][pc]));
    }
  }

  const auto snf = linalg::SmithNormalForm(poly_ring, mat);
  EXPECT_EQ(snf.rank, std::size_t{3});
  ExpectRingMatEq(
      poly_ring,
      RingMatMul(poly_ring, RingMatMul(poly_ring, snf.U, mat), snf.V), snf.S);
  ExpectUnimodular(poly_ring, snf.U);
  ExpectUnimodular(poly_ring, snf.V);

  const P xm2 = poly_ring.Sub(x, poly_ring.Constant(2));
  const P xm2_sq = poly_ring.Mul(xm2, xm2);
  ASSERT_EQ(snf.invariant_factors.size(), std::size_t{3});
  EXPECT_TRUE(poly_ring.Equals(snf.invariant_factors[0],
                               poly_ring.MultiplicativeIdentity()));
  EXPECT_TRUE(poly_ring.Equals(snf.invariant_factors[1], xm2));
  EXPECT_TRUE(poly_ring.Equals(snf.invariant_factors[2], xm2_sq));
}

// ============================================================================
// 3. FrobeniusNormalForm, MinimalPoly, CharacteristicPoly over GF(5) & GF(2^4)
// ============================================================================

TEST(PidFormsTest, FrobeniusJordanBlocksOverGF5) {
  using GF5 = GaloisFieldPrime<uint32_t>;
  auto gf5 = std::make_shared<GF5>(5);

  // Non-diagonalizable 5x5 matrix with repeated eigenvalue 3:
  // J_2(3) (+) J_2(3) (+) J_1(1), conjugated by an invertible matrix Q.
  auto j = MakeGFMatrix(gf5, 5, 5,
                        {
                            3, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3,
                            1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1,
                        });
  auto q = MakeGFMatrix(gf5, 5, 5,
                        {
                            1, 2, 0, 1, 3, 0, 1, 4, 2, 1, 2, 0, 1,
                            3, 0, 1, 1, 2, 1, 4, 3, 2, 1, 0, 2,
                        });
  ASSERT_NE(linalg::Det(q).RawValue(), uint32_t{0});
  auto a = linalg::Dot(linalg::Dot(q, j), linalg::Inv(q));

  VerifyFrobeniusAndPolys(a);

  const auto frob = linalg::FrobeniusNormalForm(a);
  // Invariant factors must be: f_1 = (x - 3)^2, f_2 = (x - 3)^2 * (x - 1).
  ASSERT_EQ(frob.invariant_factors.size(), std::size_t{2});
  EXPECT_EQ(frob.invariant_factors[0].Degree(), 2);
  EXPECT_EQ(frob.invariant_factors[1].Degree(), 3);
}

TEST(PidFormsTest, FrobeniusAndPolysOverGF16) {
  using GF16 = GaloisFieldBinaryExtension<uint16_t>;
  auto gf16 = std::make_shared<GF16>(4);

  // Non-diagonalizable 4x4 matrix over GF(2^4) with repeated eigenvalue alpha
  // (raw value 2): J_2(2) (+) J_2(2), conjugated by an invertible matrix.
  auto j = MakeGFMatrix(gf16, 4, 4,
                        {
                            2,
                            1,
                            0,
                            0,
                            0,
                            2,
                            0,
                            0,
                            0,
                            0,
                            2,
                            1,
                            0,
                            0,
                            0,
                            2,
                        });
  auto q = MakeGFMatrix(gf16, 4, 4,
                        {
                            1,
                            3,
                            5,
                            7,
                            0,
                            1,
                            2,
                            4,
                            0,
                            0,
                            1,
                            9,
                            6,
                            0,
                            0,
                            1,
                        });
  ASSERT_NE(linalg::Det(q).RawValue(), uint16_t{0});
  auto a = linalg::Dot(linalg::Dot(q, j), linalg::Inv(q));

  VerifyFrobeniusAndPolys(a);

  const auto frob = linalg::FrobeniusNormalForm(a);
  ASSERT_EQ(frob.invariant_factors.size(), std::size_t{2});
  EXPECT_EQ(frob.invariant_factors[0], frob.invariant_factors[1]);
  EXPECT_EQ(frob.invariant_factors[0].Degree(), 2);

  // Also test a random 6x6 matrix over GF(2^4).
  std::mt19937_64 rng(2026);
  auto rand_a = linalg::RandomMatrix(gf16, {6, 6}, rng);
  VerifyFrobeniusAndPolys(rand_a);
}

// ============================================================================
// 4. SparseMatrix, BlockWiedemannNullspace, BlockLanczosNullspaceBinary
// ============================================================================

TEST(SparseLinalgTest, SparseMatrixOperationsAndLdpcInterop) {
  using GF17 = GaloisFieldPrime<uint32_t>;
  auto gf17 = std::make_shared<GF17>(17);

  std::mt19937_64 rng(77);
  auto dense = linalg::RandomMatrix(gf17, {12, 15}, rng);
  auto sp = linalg::SparseMatrix<GF17>::FromDense(dense);
  EXPECT_EQ(sp.NumRows(), std::size_t{12});
  EXPECT_EQ(sp.NumCols(), std::size_t{15});
  EXPECT_EQ(sp.ToDense(), dense);

  auto sp_t = sp.Transpose();
  EXPECT_EQ(sp_t.ToDense(), linalg::Transpose(dense));

  auto x = linalg::RandomMatrix(gf17, {15}, rng);
  EXPECT_EQ(sp.SpMv(x), linalg::Dot(dense, x));

  auto y = linalg::RandomMatrix(gf17, {12}, rng);
  EXPECT_EQ(sp.SpMvTranspose(y), linalg::Dot(linalg::Transpose(dense), y));

  // Interop with coding::SparseBinaryMatrix via duck-typed factory.
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  std::vector<std::pair<std::size_t, std::size_t>> edges = {
      {0, 0}, {0, 2}, {1, 1}, {1, 3}, {2, 0}, {2, 3}};
  coding::SparseBinaryMatrix sbm(3, 4, edges);
  auto sp2 =
      linalg::SparseMatrix<GaloisFieldBinary>::FromSparseBinaryMatrix(gf2, sbm);
  EXPECT_EQ(sp2.ToDense(), sbm.ToDense(gf2));
}

TEST(SparseLinalgTest, BlockWiedemannNullspaceOverGF17AndGF256) {
  // 60 x 64 sparse matrix over GF(17).
  {
    using GF17 = GaloisFieldPrime<uint32_t>;
    auto gf17 = std::make_shared<GF17>(17);
    std::mt19937_64 rng(101);
    std::vector<std::tuple<std::size_t, std::size_t, uint32_t>> triplets;
    for (std::size_t r = 0; r < 60; ++r) {
      for (std::size_t k = 0; k < 5; ++k) {
        const std::size_t c = static_cast<std::size_t>(rng() % 64);
        uint32_t val = gf17->Random(rng);
        if (val == 0) {
          val = 1;
        }
        triplets.emplace_back(r, c, val);
      }
    }
    auto sp = linalg::SparseMatrix<GF17>::FromTriplets(gf17, 60, 64, triplets);
    auto ns = linalg::BlockWiedemannNullspace(sp, 42);
    ASSERT_GE(ns.Shape(0), std::size_t{1});
    ASSERT_EQ(ns.Shape(1), std::size_t{64});

    auto dense = sp.ToDense();
    auto prod = linalg::Dot(dense, linalg::Transpose(ns));
    EXPECT_EQ(prod, linalg::Zeros(gf17, {60, ns.Shape(0)}));

    // Every returned nullspace vector must be non-zero.
    for (std::size_t i = 0; i < ns.Shape(0); ++i) {
      bool any_nz = false;
      for (std::size_t j = 0; j < 64; ++j) {
        if (ns.Raw(i, j) != 0) {
          any_nz = true;
          break;
        }
      }
      EXPECT_TRUE(any_nz);
    }
  }

  // 60 x 64 sparse matrix over GF(2^8).
  {
    using GF256 = GaloisFieldBinaryExtension<uint8_t>;
    auto gf256 = std::make_shared<GF256>(8);
    std::mt19937_64 rng(202);
    std::vector<std::tuple<std::size_t, std::size_t, uint8_t>> triplets;
    for (std::size_t r = 0; r < 60; ++r) {
      for (std::size_t k = 0; k < 5; ++k) {
        const std::size_t c = static_cast<std::size_t>(rng() % 64);
        uint8_t val = gf256->Random(rng);
        if (val == 0) {
          val = 1;
        }
        triplets.emplace_back(r, c, val);
      }
    }
    auto sp =
        linalg::SparseMatrix<GF256>::FromTriplets(gf256, 60, 64, triplets);
    auto ns = linalg::BlockWiedemannNullspace(sp, 99);
    ASSERT_GE(ns.Shape(0), std::size_t{1});
    ASSERT_EQ(ns.Shape(1), std::size_t{64});

    auto dense = sp.ToDense();
    auto prod = linalg::Dot(dense, linalg::Transpose(ns));
    EXPECT_EQ(prod, linalg::Zeros(gf256, {60, ns.Shape(0)}));
  }
}

TEST(SparseLinalgTest, BlockLanczosNullspaceBinary120x128) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  std::mt19937_64 rng(303);

  std::vector<std::vector<std::size_t>> row_cols(120);
  for (std::size_t r = 0; r < 120; ++r) {
    // 6 distinct column entries per parity check.
    while (row_cols[r].size() < 6) {
      const std::size_t c = static_cast<std::size_t>(rng() % 128);
      if (std::find(row_cols[r].begin(), row_cols[r].end(), c) ==
          row_cols[r].end()) {
        row_cols[r].push_back(c);
      }
    }
  }

  auto sp = linalg::SparseMatrix<GaloisFieldBinary>::FromRowIndices(
      gf2, 120, 128, row_cols);
  const auto null_vecs = linalg::BlockLanczosNullspaceBinary(sp, 42);
  ASSERT_GE(null_vecs.size(), std::size_t{1});

  // Pack returned vectors into a GF(2) matrix to check full linear independence
  // and exact A * v == 0.
  GFMatrix<GaloisFieldBinary>::raw_array_type ns_raw =
      GFMatrix<GaloisFieldBinary>::raw_array_type::from_shape(
          {null_vecs.size(), std::size_t{128}});
  for (std::size_t i = 0; i < null_vecs.size(); ++i) {
    ASSERT_EQ(null_vecs[i].size(), std::size_t{128});
    bool any_one = false;
    for (std::size_t j = 0; j < 128; ++j) {
      const uint8_t bit = null_vecs[i][j];
      EXPECT_LE(bit, uint8_t{1});
      if (bit != 0) {
        any_one = true;
      }
      ns_raw.data()[i * 128 + j] = bit;
    }
    EXPECT_TRUE(any_one);
  }
  auto ns_mat = GFMatrix<GaloisFieldBinary>::FromRaw(gf2, std::move(ns_raw));
  EXPECT_EQ(linalg::MatrixRank(ns_mat), null_vecs.size());

  auto prod = linalg::Dot(sp.ToDense(), linalg::Transpose(ns_mat));
  EXPECT_EQ(prod, linalg::Zeros(gf2, {120, null_vecs.size()}));
}

// ============================================================================
// 5. M4riMatMul, M4riRref, StrassenMatMul, FastMatMul
// ============================================================================

TEST(FastMatMulTest, M4riMatMulAndM4riRrefMatchReferenceBitForBit) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  std::mt19937_64 rng(404);

  // Odd and multi-word dimensions across 64-bit word boundaries.
  auto a = linalg::RandomMatrix(gf2, {73, 131}, rng);
  auto b = linalg::RandomMatrix(gf2, {131, 97}, rng);
  EXPECT_EQ(linalg::M4riMatMul(a, b), linalg::Dot(a, b));
  EXPECT_EQ(linalg::FastMatMul(a, b), linalg::Dot(a, b));

  // Full-rank and rank-deficient matrices for M4riRref.
  auto r_full = linalg::RandomMatrix(gf2, {85, 140}, rng);
  auto [rref_m4ri, pivots_m4ri] = linalg::M4riRref(r_full);
  auto rref_ref = linalg::Rref(r_full);
  EXPECT_EQ(rref_m4ri, rref_ref);
  EXPECT_EQ(pivots_m4ri.size(), linalg::MatrixRank(r_full));

  // Explicitly rank-deficient matrix (product of 90x37 and 37x135).
  auto low_a = linalg::RandomMatrix(gf2, {90, 37}, rng);
  auto low_b = linalg::RandomMatrix(gf2, {37, 135}, rng);
  auto low_rank = linalg::Dot(low_a, low_b);
  auto [low_rref, low_pivs] = linalg::M4riRref(low_rank);
  EXPECT_EQ(low_rref, linalg::Rref(low_rank));
  EXPECT_EQ(low_pivs.size(), linalg::MatrixRank(low_rank));
}

TEST(FastMatMulTest, StrassenAndFastMatMulOddAndPowerOfTwoDims) {
  std::mt19937_64 rng(505);

  // GF(2): odd dimensions 129 x 131 x 127.
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto a = linalg::RandomMatrix(gf2, {129, 131}, rng);
    auto b = linalg::RandomMatrix(gf2, {131, 127}, rng);
    auto ref = linalg::Dot(a, b);
    EXPECT_EQ(linalg::StrassenMatMul(a, b, 32), ref);
    EXPECT_EQ(linalg::FastMatMul(a, b), ref);
  }

  // GF(2^8): odd dimensions 129 x 131 x 127 and power-of-2 128 x 128 x 128.
  {
    using GF256 = GaloisFieldBinaryExtension<uint8_t>;
    auto gf256 = std::make_shared<GF256>(8);

    auto a_odd = linalg::RandomMatrix(gf256, {129, 131}, rng);
    auto b_odd = linalg::RandomMatrix(gf256, {131, 127}, rng);
    auto ref_odd = linalg::Dot(a_odd, b_odd);
    EXPECT_EQ(linalg::StrassenMatMul(a_odd, b_odd, 64), ref_odd);
    EXPECT_EQ(linalg::FastMatMul(a_odd, b_odd), ref_odd);

    auto a_pow2 = linalg::RandomMatrix(gf256, {128, 128}, rng);
    auto b_pow2 = linalg::RandomMatrix(gf256, {128, 128}, rng);
    auto ref_pow2 = linalg::Dot(a_pow2, b_pow2);
    EXPECT_EQ(linalg::StrassenMatMul(a_pow2, b_pow2, 64), ref_pow2);
    EXPECT_EQ(linalg::FastMatMul(a_pow2, b_pow2), ref_pow2);
  }

  // GF(257): odd dimensions 129 x 131 x 127 and medium 45 x 53 x 49.
  {
    using GF257 = GaloisFieldPrime<uint32_t>;
    auto gf257 = std::make_shared<GF257>(257);

    auto a_odd = linalg::RandomMatrix(gf257, {129, 131}, rng);
    auto b_odd = linalg::RandomMatrix(gf257, {131, 127}, rng);
    auto ref_odd = linalg::Dot(a_odd, b_odd);
    EXPECT_EQ(linalg::StrassenMatMul(a_odd, b_odd, 64), ref_odd);
    EXPECT_EQ(linalg::FastMatMul(a_odd, b_odd), ref_odd);

    auto a_med = linalg::RandomMatrix(gf257, {45, 53}, rng);
    auto b_med = linalg::RandomMatrix(gf257, {53, 49}, rng);
    EXPECT_EQ(linalg::FastMatMul(a_med, b_med), linalg::Dot(a_med, b_med));
  }
}
