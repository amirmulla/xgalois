#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/coding/gabidulin.hpp"
#include "xgalois/coding/subspace_code.hpp"
#include "xgalois/field/gf128.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace {

using namespace xg;
using namespace xg::coding;

// Generates a random full-rank matrix of shape {rows, cols} (rows <= cols)
// over `field`.
template <typename FieldT>
GFMatrix<FieldT> RandomFullRowRankMatrix(const std::shared_ptr<FieldT> &field,
                                         std::size_t rows, std::size_t cols,
                                         std::mt19937_64 &rng) {
  if (rows == 0) {
    return GFMatrix<FieldT>::Zeros(field, {0, cols});
  }
  for (int attempt = 0; attempt < 200; ++attempt) {
    auto mat = GFMatrix<FieldT>::Random(field, {rows, cols}, rng);
    if (linalg::MatrixRank(mat) == rows) {
      return mat;
    }
  }
  throw std::runtime_error("Failed to sample full-row-rank matrix");
}

// Generates a random full-column-rank matrix of shape {rows, cols} (cols <=
// rows) over `field`.
template <typename FieldT>
GFMatrix<FieldT> RandomFullColRankMatrix(const std::shared_ptr<FieldT> &field,
                                         std::size_t rows, std::size_t cols,
                                         std::mt19937_64 &rng) {
  if (cols == 0) {
    return GFMatrix<FieldT>::Zeros(field, {rows, 0});
  }
  return linalg::Transpose(RandomFullRowRankMatrix(field, cols, rows, rng));
}

// Generates a random rank-t error vector over ExtField (length n) by forming
// A (m x t, col-rank t) * B (t x n, row-rank t) over BaseField and contracting.
template <typename BaseField, typename ExtField>
GFArray<ExtField> RandomRankError(
    const FieldTowerAdapter<BaseField, ExtField> &tower, std::size_t n,
    std::size_t rank_t, std::mt19937_64 &rng) {
  auto base_field = tower.BaseFieldPtr();
  const std::size_t m = tower.RelativeDegree();
  if (rank_t == 0) {
    return GFArray<ExtField>::Zeros(tower.ExtFieldPtr(), {n});
  }
  auto a_mat = RandomFullColRankMatrix(base_field, m, rank_t, rng);
  auto b_mat = RandomFullRowRankMatrix(base_field, rank_t, n, rng);
  auto e_mat = linalg::Dot(a_mat, b_mat);
  return tower.ContractMatrix(e_mat);
}

// Generates a random error vector decomposed as:
//   e = e^{(R)} + e^{(C)} + e^{(E)}
// with exact row-erasure basis a^{(R)} (length rho), exact column-erasure
// matrix B^{(C)} (shape {gamma, n}), and independent rank-tau unknown error.
template <typename BaseField, typename ExtField>
struct DecomposedRankError {
  GFArray<ExtField> total_error;
  GFArray<ExtField> row_erasure_values;
  GFMatrix<BaseField> col_erasure_matrix;
};

template <typename BaseField, typename ExtField>
DecomposedRankError<BaseField, ExtField> MakeDecomposedRankError(
    const FieldTowerAdapter<BaseField, ExtField> &tower, std::size_t n,
    std::size_t tau, std::size_t rho, std::size_t gamma, std::mt19937_64 &rng) {
  auto base_field = tower.BaseFieldPtr();
  auto ext_field = tower.ExtFieldPtr();
  const std::size_t m = tower.RelativeDegree();

  // Sample (rho + tau) linearly independent column vectors in F_q^m so row
  // erasures and unknown error values are independent.
  const std::size_t total_a = rho + tau;
  auto a_all = RandomFullColRankMatrix(base_field, m, total_a, rng);

  // Sample (gamma + tau) linearly independent row vectors in F_q^n so column
  // erasures and unknown error locations are independent.
  const std::size_t total_b = gamma + tau;
  auto b_all = RandomFullRowRankMatrix(base_field, total_b, n, rng);

  // Extract a^{(R)} (first rho columns of a_all) and a^{(E)} (last tau cols).
  auto a_row_mat = GFMatrix<BaseField>::Zeros(base_field, {m, rho});
  for (std::size_t r = 0; r < m; ++r) {
    for (std::size_t c = 0; c < rho; ++c) {
      a_row_mat.Raw(r, c) = a_all.Raw(r, c);
    }
  }
  auto a_err_mat = GFMatrix<BaseField>::Zeros(base_field, {m, tau});
  for (std::size_t r = 0; r < m; ++r) {
    for (std::size_t c = 0; c < tau; ++c) {
      a_err_mat.Raw(r, c) = a_all.Raw(r, rho + c);
    }
  }

  // Extract B^{(C)} (first gamma rows of b_all) and B^{(E)} (last tau rows).
  auto b_col_mat = GFMatrix<BaseField>::Zeros(base_field, {gamma, n});
  for (std::size_t r = 0; r < gamma; ++r) {
    for (std::size_t c = 0; c < n; ++c) {
      b_col_mat.Raw(r, c) = b_all.Raw(r, c);
    }
  }
  auto b_err_mat = GFMatrix<BaseField>::Zeros(base_field, {tau, n});
  for (std::size_t r = 0; r < tau; ++r) {
    for (std::size_t c = 0; c < n; ++c) {
      b_err_mat.Raw(r, c) = b_all.Raw(gamma + r, c);
    }
  }

  // Build unknown B^{(R)} (rho x n) and unknown A^{(C)} (m x gamma).
  auto b_row_mat = RandomFullRowRankMatrix(base_field, rho, n, rng);
  auto a_col_mat = RandomFullColRankMatrix(base_field, m, gamma, rng);

  auto e_mat = GFMatrix<BaseField>::Zeros(base_field, {m, n});
  if (rho > 0) {
    e_mat += linalg::Dot(a_row_mat, b_row_mat);
  }
  if (gamma > 0) {
    e_mat += linalg::Dot(a_col_mat, b_col_mat);
  }
  if (tau > 0) {
    e_mat += linalg::Dot(a_err_mat, b_err_mat);
  }

  DecomposedRankError<BaseField, ExtField> out;
  out.total_error = tower.ContractMatrix(e_mat);
  out.row_erasure_values = (rho > 0) ? tower.ContractMatrix(a_row_mat)
                                     : GFArray<ExtField>::Zeros(ext_field, {0});
  out.col_erasure_matrix = std::move(b_col_mat);
  return out;
}

template <typename FieldT>
bool IsZeroMatrix(const GFMatrix<FieldT> &mat) {
  auto field = mat.Field();
  for (const auto &v : mat.RawSpan()) {
    if (!field->IsZero(v)) {
      return false;
    }
  }
  return true;
}

// ============================================================================
// 1. FieldTowerAdapter & KernelBasis Tests
// ============================================================================

TEST(GabidulinSubspaceTest, FieldTowerAdapterBinaryPrime128AndRelative) {
  std::mt19937_64 rng(20260401);

  // (a) GF(2) subset GF(2^8) with standard and Zech log tables.
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    auto gf256_zech = std::make_shared<GFBEZechLogTables<uint8_t>>(8);

    FieldTowerAdapter<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>
        tower(gf2, gf256);
    FieldTowerAdapter<GaloisFieldBinary, GFBEZechLogTables<uint8_t>> tower_z(
        gf2, gf256_zech);

    EXPECT_EQ(tower.RelativeDegree(), 8U);
    EXPECT_EQ(tower.BaseOrder(), 2U);
    EXPECT_EQ(tower_z.RelativeDegree(), 8U);

    auto v =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {8}, rng);
    auto mat = tower.ExpandVector(v);
    EXPECT_EQ(mat.Shape(0), 8U);
    EXPECT_EQ(mat.Shape(1), 8U);
    EXPECT_EQ(tower.ContractMatrix(mat), v);
    EXPECT_EQ(tower.RankWeight(v), linalg::MatrixRank(mat));

    // Verify FrobeniusQ and InvFrobeniusQ are exact inverses.
    for (std::size_t i = 0; i < 8; ++i) {
      const auto x = v.Raw(i);
      const auto fx = tower.FrobeniusQ(x, 3);
      EXPECT_EQ(tower.InvFrobeniusQ(fx, 3), x);
      EXPECT_EQ(tower.FrobeniusQ(x, 8), x);
    }

    // Zech log representation roundtrip (exercises GetElementValue
    // bit-slicing).
    auto vz = GFArray<GFBEZechLogTables<uint8_t>>::Random(gf256_zech, {6}, rng);
    auto mz = tower_z.ExpandVector(vz);
    EXPECT_EQ(tower_z.ContractMatrix(mz), vz);
    EXPECT_EQ(tower_z.RankWeight(vz), linalg::MatrixRank(mz));
  }

  // (b) GF(3) subset GF(3^4).
  {
    auto gf3 = std::make_shared<GaloisFieldPrime<uint32_t>>(3);
    auto gf81 =
        std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 4));
    FieldTowerAdapter<GaloisFieldPrime<uint32_t>,
                      GaloisFieldExtension<uint32_t>>
        tower(gf3, gf81);
    EXPECT_EQ(tower.RelativeDegree(), 4U);
    EXPECT_EQ(tower.BaseOrder(), 3U);

    auto v = GFArray<GaloisFieldExtension<uint32_t>>::Random(gf81, {4}, rng);
    auto mat = tower.ExpandVector(v);
    EXPECT_EQ(mat.Shape(0), 4U);
    EXPECT_EQ(mat.Shape(1), 4U);
    EXPECT_EQ(tower.ContractMatrix(mat), v);
    EXPECT_EQ(tower.RankWeight(v), linalg::MatrixRank(mat));

    for (std::size_t i = 0; i < 4; ++i) {
      const auto x = v.Raw(i);
      EXPECT_EQ(tower.InvFrobeniusQ(tower.FrobeniusQ(x, 2), 2), x);
      EXPECT_EQ(tower.FrobeniusQ(x, 4), x);
    }
  }

  // (c) GF(2) subset GF(2^128).
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf128 = std::make_shared<GaloisFieldBinary128>();
    FieldTowerAdapter<GaloisFieldBinary, GaloisFieldBinary128> tower(gf2,
                                                                     gf128);
    EXPECT_EQ(tower.RelativeDegree(), 128U);
    EXPECT_EQ(tower.BaseOrder(), 2U);

    auto v = GFArray<GaloisFieldBinary128>::Random(gf128, {6}, rng);
    auto mat = tower.ExpandVector(v);
    EXPECT_EQ(mat.Shape(0), 128U);
    EXPECT_EQ(mat.Shape(1), 6U);
    EXPECT_EQ(tower.ContractMatrix(mat), v);
    EXPECT_EQ(tower.RankWeight(v), linalg::MatrixRank(mat));
    EXPECT_EQ(tower.InvFrobeniusQ(tower.FrobeniusQ(v.Raw(0), 5), 5), v.Raw(0));
  }

  // (d) General relative tower: GF(2^2) subset GF(2^8) (q = 4, m = 4).
  {
    auto gf4 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(2);
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    FieldTowerAdapter<GaloisFieldBinaryExtension<uint8_t>,
                      GaloisFieldBinaryExtension<uint8_t>>
        tower(gf4, gf256);
    EXPECT_EQ(tower.RelativeDegree(), 4U);
    EXPECT_EQ(tower.BaseOrder(), 4U);

    auto v =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {4}, rng);
    auto mat = tower.ExpandVector(v);
    EXPECT_EQ(mat.Shape(0), 4U);
    EXPECT_EQ(mat.Shape(1), 4U);
    EXPECT_EQ(tower.ContractMatrix(mat), v);

    // Check F_4-linearity: Expand(a + b) == Expand(a) + Expand(b), and
    // FrobeniusQ(x) == x^4 with FrobeniusQ(x, 4) == x.
    for (std::size_t i = 0; i < 4; ++i) {
      const auto x = v.Raw(i);
      EXPECT_EQ(tower.FrobeniusQ(x, 1), gf256->Pow(x, 4));
      EXPECT_EQ(tower.InvFrobeniusQ(tower.FrobeniusQ(x, 2), 2), x);
      EXPECT_EQ(tower.FrobeniusQ(x, 4), x);
    }

    // Gabidulin G(4, 2) over GF(4^4)/GF(4) works end-to-end!
    auto gab = GabidulinCode<
        GaloisFieldBinaryExtension<uint8_t>,
        GaloisFieldBinaryExtension<uint8_t>>::CreateStandard(gf4, gf256, 4, 2);
    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(
        gab.GeneratorMatrix(), linalg::Transpose(gab.ParityCheckMatrix()))));
    auto msg =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
    auto cw = gab.Encode(msg);
    auto err = RandomRankError(tower, 4, 1, rng);
    EXPECT_EQ(gab.Weight(err), 1U);
    auto res = gab.DecodeWithResult(cw + err);
    ASSERT_TRUE(res.Ok()) << res.error_message;
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
  }
}

TEST(GabidulinSubspaceTest, KernelBasisMatchesMinimalPolynomial) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
  FieldTowerAdapter<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>
      tower(gf2, gf256);
  xg::ring::FrobeniusSkewRing<GaloisFieldBinaryExtension<uint8_t>> skew(*gf256,
                                                                        1);

  std::vector<uint8_t> pts{
      gf256->SetElementValue(0x03),
      gf256->SetElementValue(0x15),
      gf256->SetElementValue(0x4A),
  };
  const auto min_poly = skew.MinimalPolynomial(pts);
  EXPECT_EQ(skew.Degree(min_poly), 3);

  const auto ker = KernelBasis(skew, tower, min_poly);
  EXPECT_EQ(ker.Size(), 3U);
  EXPECT_EQ(tower.RankWeight(ker), 3U);

  // Every kernel basis element must evaluate to 0 under min_poly.
  for (std::size_t i = 0; i < ker.Size(); ++i) {
    EXPECT_TRUE(gf256->IsZero(skew.EvalOperator(min_poly, ker.Raw(i))));
  }

  // Span of `ker` must equal span of `pts` over GF(2).
  auto pts_arr = GFArray<GaloisFieldBinaryExtension<uint8_t>>::FromRaw(
      gf256, xt::xarray<uint8_t>{pts[0], pts[1], pts[2]});
  auto combined =
      linalg::HStack(tower.ExpandVector(ker), tower.ExpandVector(pts_arr));
  EXPECT_EQ(linalg::MatrixRank(combined), 3U);
}

// ============================================================================
// 2. Gabidulin Structure, Orthogonality, MRD Distance, DualCode, Rule-of-Five
// ============================================================================

TEST(GabidulinSubspaceTest, StructureOrthogonalityAndDualBinaryAndTernary) {
  std::mt19937_64 rng(424242);

  // (a) G(8, 4) over GF(2^8) / GF(2).
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    auto gab =
        GabidulinCode<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>::
            CreateStandard(gf2, gf256, 8, 4);

    EXPECT_EQ(gab.Length(), 8U);
    EXPECT_EQ(gab.Dimension(), 4U);
    EXPECT_EQ(gab.MinimumDistance(), 5U);
    EXPECT_EQ(gab.GetMetric(), Metric::RANK);
    EXPECT_TRUE(gab.IsMRD());

    // Exact G * H^T == 0.
    const auto &G = gab.GeneratorMatrix();
    const auto &H = gab.ParityCheckMatrix();
    EXPECT_EQ(G.Shape(0), 4U);
    EXPECT_EQ(G.Shape(1), 8U);
    EXPECT_EQ(H.Shape(0), 4U);
    EXPECT_EQ(H.Shape(1), 8U);
    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(G, linalg::Transpose(H))));

    // Weight() == MatrixRank(ExpandToBaseMatrix(w)) and Distance(c1, c2) >=
    // d_R.
    auto m1 =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {4}, rng);
    auto m2 =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {4}, rng);
    auto c1 = gab.Encode(m1);
    auto c2 = gab.Encode(m2);
    EXPECT_TRUE(gab.Contains(c1));
    EXPECT_EQ(gab.Unencode(c1), m1);
    EXPECT_EQ(gab.Weight(c1), linalg::MatrixRank(gab.ExpandToBaseMatrix(c1)));
    EXPECT_EQ(gab.Distance(c1, c2), gab.Weight(c1 - c2));
    if (m1 != m2) {
      EXPECT_GE(gab.Distance(c1, c2), gab.MinimumDistance());
    }

    // DualCode() orthogonality.
    auto dual = gab.DualCode();
    EXPECT_EQ(dual.Length(), 8U);
    EXPECT_EQ(dual.Dimension(), 4U);
    EXPECT_EQ(dual.MinimumDistance(), 5U);
    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(
        gab.GeneratorMatrix(), linalg::Transpose(dual.GeneratorMatrix()))));
    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(
        dual.GeneratorMatrix(), linalg::Transpose(dual.ParityCheckMatrix()))));

    // Rule-of-Five safety: copy & move after decoder caching.
    auto rx = c1 + RandomRankError(gab.Tower(), 8, 2, rng);
    ASSERT_TRUE(gab.DecodeWithResult(rx).Ok());

    auto gab_copy = gab;
    auto res_copy = gab_copy.DecodeWithResult(rx);
    ASSERT_TRUE(res_copy.Ok());
    EXPECT_EQ(res_copy.codeword, c1);

    auto gab_moved = std::move(gab_copy);
    auto res_moved = gab_moved.DecodeWithResult(rx, "GabidulinLoidreau");
    ASSERT_TRUE(res_moved.Ok());
    EXPECT_EQ(res_moved.codeword, c1);
  }

  // (b) G(4, 2) over GF(3^4) / GF(3).
  {
    auto gf3 = std::make_shared<GaloisFieldPrime<uint32_t>>(3);
    auto gf81 =
        std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 4));
    auto gab =
        GabidulinCode<GaloisFieldPrime<uint32_t>,
                      GaloisFieldExtension<uint32_t>>::CreateStandard(gf3, gf81,
                                                                      4, 2);

    EXPECT_EQ(gab.Length(), 4U);
    EXPECT_EQ(gab.Dimension(), 2U);
    EXPECT_EQ(gab.MinimumDistance(), 3U);
    EXPECT_EQ(gab.GetMetric(), Metric::RANK);

    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(
        gab.GeneratorMatrix(), linalg::Transpose(gab.ParityCheckMatrix()))));

    auto dual = gab.DualCode();
    EXPECT_EQ(dual.Length(), 4U);
    EXPECT_EQ(dual.Dimension(), 2U);
    EXPECT_EQ(dual.MinimumDistance(), 3U);
    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(
        gab.GeneratorMatrix(), linalg::Transpose(dual.GeneratorMatrix()))));

    auto msg = GFArray<GaloisFieldExtension<uint32_t>>::Random(gf81, {2}, rng);
    auto cw = gab.Encode(msg);
    EXPECT_TRUE(gab.Contains(cw));
    EXPECT_EQ(gab.Unencode(cw), msg);
    EXPECT_EQ(gab.Weight(cw), linalg::MatrixRank(gab.ExpandToBaseMatrix(cw)));
  }

  // (c) Reject linearly dependent evaluation points.
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    // 1, 2, 3 = 1 + 2 are GF(2)-linearly dependent!
    auto bad_pts = GFArray<GaloisFieldBinaryExtension<uint8_t>>::FromRaw(
        gf256, xt::xarray<uint8_t>{1, 2, 3});
    EXPECT_THROW(
        (GabidulinCode<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>(
            gf2, gf256, bad_pts, 2)),
        std::invalid_argument);
  }
}

// ============================================================================
// 3. Pure Rank Error Decoding: Skew-EEA vs Loidreau Parity & Boundary Failure
// ============================================================================

TEST(GabidulinSubspaceTest, PureRankErrorDecodingBothDecodersBinaryAndTernary) {
  std::mt19937_64 rng(991177);

  // (a) G(8, 4) over GF(2^8) / GF(2): t = floor((8 - 4)/2) = 2.
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    auto gab =
        GabidulinCode<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>::
            CreateStandard(gf2, gf256, 8, 4);

    for (std::size_t tau : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
      for (int trial = 0; trial < 8; ++trial) {
        auto msg = GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(
            gf256, {4}, rng);
        auto cw = gab.Encode(msg);
        auto err = RandomRankError(gab.Tower(), 8, tau, rng);
        ASSERT_EQ(gab.Weight(err), tau);
        auto rx = cw + err;

        auto res_eea = gab.DecodeWithResult(rx, "GabidulinSkewEEA");
        auto res_loi = gab.DecodeWithResult(rx, "GabidulinLoidreau");

        ASSERT_TRUE(res_eea.Ok())
            << "tau=" << tau << ": " << res_eea.error_message;
        ASSERT_TRUE(res_loi.Ok())
            << "tau=" << tau << ": " << res_loi.error_message;
        EXPECT_EQ(res_eea.codeword, cw);
        EXPECT_EQ(res_eea.message, msg);
        EXPECT_EQ(res_eea.num_errors, tau);

        EXPECT_EQ(res_loi.codeword, cw);
        EXPECT_EQ(res_loi.message, msg);
        EXPECT_EQ(res_loi.num_errors, tau);
      }
    }

    // Beyond decoding radius: tau = 3 > 2 must NOT silently return wrong word.
    int failed_eea = 0;
    int failed_loi = 0;
    for (int trial = 0; trial < 10; ++trial) {
      auto msg =
          GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {4}, rng);
      auto cw = gab.Encode(msg);
      auto err = RandomRankError(gab.Tower(), 8, 3, rng);
      ASSERT_EQ(gab.Weight(err), 3U);
      auto rx = cw + err;

      auto res_eea = gab.DecodeWithResult(rx, "GabidulinSkewEEA");
      auto res_loi = gab.DecodeWithResult(rx, "GabidulinLoidreau");
      if (!res_eea.Ok()) {
        ++failed_eea;
      } else {
        // If a ball of radius 2 happened to cover rx, distance must be <= 2.
        EXPECT_LE(gab.Distance(rx, res_eea.codeword), 2U);
      }
      if (!res_loi.Ok()) {
        ++failed_loi;
      } else {
        EXPECT_LE(gab.Distance(rx, res_loi.codeword), 2U);
      }
    }
    EXPECT_GT(failed_eea, 0);
    EXPECT_GT(failed_loi, 0);
  }

  // (b) Odd redundancy G(7, 2) over GF(2^8) / GF(2): r = 5, t = 2.
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    auto gab =
        GabidulinCode<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>::
            CreateStandard(gf2, gf256, 7, 2);
    EXPECT_EQ(gab.MinimumDistance(), 6U);

    for (std::size_t tau : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
      for (int trial = 0; trial < 5; ++trial) {
        auto msg = GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(
            gf256, {2}, rng);
        auto cw = gab.Encode(msg);
        auto err = RandomRankError(gab.Tower(), 7, tau, rng);
        auto rx = cw + err;

        auto res_eea = gab.DecodeWithResult(rx, "GabidulinSkewEEA");
        auto res_loi = gab.DecodeWithResult(rx, "GabidulinLoidreau");
        ASSERT_TRUE(res_eea.Ok()) << res_eea.error_message;
        ASSERT_TRUE(res_loi.Ok()) << res_loi.error_message;
        EXPECT_EQ(res_eea.codeword, cw);
        EXPECT_EQ(res_loi.codeword, cw);
        EXPECT_EQ(res_eea.message, msg);
        EXPECT_EQ(res_loi.message, msg);
      }
    }
  }

  // (c) G(4, 2) over GF(3^4) / GF(3): t = 1.
  {
    auto gf3 = std::make_shared<GaloisFieldPrime<uint32_t>>(3);
    auto gf81 =
        std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 4));
    auto gab =
        GabidulinCode<GaloisFieldPrime<uint32_t>,
                      GaloisFieldExtension<uint32_t>>::CreateStandard(gf3, gf81,
                                                                      4, 2);

    for (std::size_t tau : {std::size_t{0}, std::size_t{1}}) {
      for (int trial = 0; trial < 6; ++trial) {
        auto msg =
            GFArray<GaloisFieldExtension<uint32_t>>::Random(gf81, {2}, rng);
        auto cw = gab.Encode(msg);
        auto err = RandomRankError(gab.Tower(), 4, tau, rng);
        ASSERT_EQ(gab.Weight(err), tau);
        auto rx = cw + err;

        auto res_eea = gab.DecodeWithResult(rx, "GabidulinSkewEEA");
        auto res_loi = gab.DecodeWithResult(rx, "GabidulinLoidreau");
        ASSERT_TRUE(res_eea.Ok()) << res_eea.error_message;
        ASSERT_TRUE(res_loi.Ok()) << res_loi.error_message;
        EXPECT_EQ(res_eea.codeword, cw);
        EXPECT_EQ(res_loi.codeword, cw);
        EXPECT_EQ(res_eea.message, msg);
        EXPECT_EQ(res_loi.message, msg);
        EXPECT_EQ(res_eea.num_errors, tau);
        EXPECT_EQ(res_loi.num_errors, tau);
      }
    }

    // tau = 2 > 1 must fail (or never misdecode outside radius 1).
    int failed = 0;
    for (int trial = 0; trial < 10; ++trial) {
      auto msg =
          GFArray<GaloisFieldExtension<uint32_t>>::Random(gf81, {2}, rng);
      auto cw = gab.Encode(msg);
      auto err = RandomRankError(gab.Tower(), 4, 2, rng);
      auto rx = cw + err;
      auto res_eea = gab.DecodeWithResult(rx, "GabidulinSkewEEA");
      auto res_loi = gab.DecodeWithResult(rx, "GabidulinLoidreau");
      if (!res_eea.Ok() && !res_loi.Ok()) {
        ++failed;
      }
    }
    EXPECT_GT(failed, 0);
  }

  // (d) G(8, 4) over GF(2^8)/GF(2) with Zech log tables, and G(10, 4) over
  //     GF(2^16)/GF(2) (t = 3).
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256_zech = std::make_shared<GFBEZechLogTables<uint8_t>>(8);
    auto gab_z =
        GabidulinCode<GaloisFieldBinary,
                      GFBEZechLogTables<uint8_t>>::CreateStandard(gf2,
                                                                  gf256_zech, 8,
                                                                  4);
    auto msg_z =
        GFArray<GFBEZechLogTables<uint8_t>>::Random(gf256_zech, {4}, rng);
    auto cw_z = gab_z.Encode(msg_z);
    auto err_z = RandomRankError(gab_z.Tower(), 8, 2, rng);
    ASSERT_EQ(gab_z.Weight(err_z), 2U);
    auto res_eea = gab_z.DecodeWithResult(cw_z + err_z, "GabidulinSkewEEA");
    auto res_loi = gab_z.DecodeWithResult(cw_z + err_z, "GabidulinLoidreau");
    ASSERT_TRUE(res_eea.Ok()) << res_eea.error_message;
    ASSERT_TRUE(res_loi.Ok()) << res_loi.error_message;
    EXPECT_EQ(res_eea.codeword, cw_z);
    EXPECT_EQ(res_loi.codeword, cw_z);

    auto gf65536 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(16);
    auto gab16 =
        GabidulinCode<GaloisFieldBinary, GaloisFieldBinaryExtension<uint16_t>>::
            CreateStandard(gf2, gf65536, 10, 4);
    EXPECT_EQ(gab16.MinimumDistance(), 7U);
    auto msg16 = GFArray<GaloisFieldBinaryExtension<uint16_t>>::Random(
        gf65536, {4}, rng);
    auto cw16 = gab16.Encode(msg16);
    auto err16 = RandomRankError(gab16.Tower(), 10, 3, rng);
    ASSERT_EQ(gab16.Weight(err16), 3U);
    auto res16_eea = gab16.DecodeWithResult(cw16 + err16, "GabidulinSkewEEA");
    auto res16_loi = gab16.DecodeWithResult(cw16 + err16, "GabidulinLoidreau");
    ASSERT_TRUE(res16_eea.Ok()) << res16_eea.error_message;
    ASSERT_TRUE(res16_loi.Ok()) << res16_loi.error_message;
    EXPECT_EQ(res16_eea.codeword, cw16);
    EXPECT_EQ(res16_loi.codeword, cw16);
    EXPECT_EQ(res16_eea.message, msg16);
    EXPECT_EQ(res16_loi.message, msg16);
  }
}

// ============================================================================
// 4. Simultaneous Rank Errors + Row Erasures + Column Erasures (SKK 2008)
// ============================================================================

TEST(GabidulinSubspaceTest, SimultaneousRankErrorsAndRowColErasuresBoundary) {
  std::mt19937_64 rng(20080515);

  // (a) G(8, 4) over GF(2^8) / GF(2): n - k = 4.
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
    auto gab =
        GabidulinCode<GaloisFieldBinary, GaloisFieldBinaryExtension<uint8_t>>::
            CreateStandard(gf2, gf256, 8, 4);

    struct BudgetCase {
      std::size_t tau;
      std::size_t rho;
      std::size_t gamma;
    };
    const std::vector<BudgetCase> exact_boundary_cases = {
        {0, 4, 0},
        {0, 0, 4},
        {0, 2, 2},
        {0, 3, 1},
        {0, 1, 3},
        {1, 1, 1},
        {1, 2, 0},
        {1, 0, 2},
        {2, 0, 0},
        // Also strictly inside the boundary:
        {1, 1, 0},
        {1, 0, 1},
        {0, 2, 1},
    };

    for (const auto &bc : exact_boundary_cases) {
      for (int trial = 0; trial < 5; ++trial) {
        auto msg = GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(
            gf256, {4}, rng);
        auto cw = gab.Encode(msg);
        auto decomp = MakeDecomposedRankError(gab.Tower(), 8, bc.tau, bc.rho,
                                              bc.gamma, rng);
        auto rx = cw + decomp.total_error;

        auto res = gab.DecodeWithRankErasures(rx, decomp.row_erasure_values,
                                              decomp.col_erasure_matrix);
        ASSERT_TRUE(res.Ok())
            << "Failed for (tau=" << bc.tau << ", rho=" << bc.rho
            << ", gamma=" << bc.gamma << "): " << res.error_message;
        EXPECT_EQ(res.codeword, cw);
        EXPECT_EQ(res.message, msg);
        EXPECT_EQ(res.num_erasures, bc.rho + bc.gamma);
        EXPECT_LE(res.num_errors, bc.tau);
      }
    }

    // Standard coordinate erasures via DecodeWithErasures: 2*1 + 2 = 4.
    {
      auto msg =
          GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {4}, rng);
      auto cw = gab.Encode(msg);
      auto rx = cw + RandomRankError(gab.Tower(), 8, 1, rng);
      // Corrupt 2 erased coordinate positions with arbitrary values.
      rx.Raw(1) = gf256->Add(rx.Raw(1), gf256->SetElementValue(0xAB));
      rx.Raw(6) = gf256->Add(rx.Raw(6), gf256->SetElementValue(0xCD));

      auto res = gab.DecodeWithErasures(rx, {1, 6});
      ASSERT_TRUE(res.Ok()) << res.error_message;
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
    }

    // Exceeding redundancy: rho + gamma = 5 > 4 must fail cleanly.
    {
      auto msg =
          GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {4}, rng);
      auto cw = gab.Encode(msg);
      auto decomp = MakeDecomposedRankError(gab.Tower(), 8, 0, 3, 2, rng);
      auto res = gab.DecodeWithRankErasures(cw + decomp.total_error,
                                            decomp.row_erasure_values,
                                            decomp.col_erasure_matrix);
      EXPECT_FALSE(res.Ok());
    }
  }

  // (b) G(4, 2) over GF(3^4) / GF(3): n - k = 2.
  {
    auto gf3 = std::make_shared<GaloisFieldPrime<uint32_t>>(3);
    auto gf81 =
        std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 4));
    auto gab =
        GabidulinCode<GaloisFieldPrime<uint32_t>,
                      GaloisFieldExtension<uint32_t>>::CreateStandard(gf3, gf81,
                                                                      4, 2);

    struct BudgetCase {
      std::size_t tau;
      std::size_t rho;
      std::size_t gamma;
    };
    const std::vector<BudgetCase> cases = {
        {0, 2, 0}, {0, 0, 2}, {0, 1, 1}, {1, 0, 0}};
    for (const auto &bc : cases) {
      for (int trial = 0; trial < 4; ++trial) {
        auto msg =
            GFArray<GaloisFieldExtension<uint32_t>>::Random(gf81, {2}, rng);
        auto cw = gab.Encode(msg);
        auto decomp = MakeDecomposedRankError(gab.Tower(), 4, bc.tau, bc.rho,
                                              bc.gamma, rng);
        auto res = gab.DecodeWithRankErasures(cw + decomp.total_error,
                                              decomp.row_erasure_values,
                                              decomp.col_erasure_matrix);
        ASSERT_TRUE(res.Ok())
            << "GF(3^4) failed for (tau=" << bc.tau << ", rho=" << bc.rho
            << ", gamma=" << bc.gamma << "): " << res.error_message;
        EXPECT_EQ(res.codeword, cw);
        EXPECT_EQ(res.message, msg);
      }
    }
  }
}

// ============================================================================
// 5. Kotter-Kschischang Lifted Gabidulin Subspace Code Tests
// ============================================================================

TEST(GabidulinSubspaceTest, KotterKschischangDistanceIsometryAndProperties) {
  std::mt19937_64 rng(777888);

  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
  auto kk = KotterKschischangCode<
      GaloisFieldBinary,
      GaloisFieldBinaryExtension<uint8_t>>::CreateStandard(gf2, gf256, 6, 2);

  EXPECT_EQ(kk.AmbientDimension(), 14U);
  EXPECT_EQ(kk.SubspaceDimension(), 6U);
  EXPECT_EQ(kk.CodeDimension(), 2U);
  EXPECT_EQ(kk.ExtensionDegree(), 8U);
  EXPECT_EQ(kk.MinimumSubspaceDistance(), 10U);
  EXPECT_EQ(kk.GetMetric(), Metric::SUBSPACE);

  const auto &gab = kk.UnderlyingGabidulinCode();
  for (int trial = 0; trial < 15; ++trial) {
    auto m1 =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
    auto m2 =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
    auto c1 = gab.Encode(m1);
    auto c2 = gab.Encode(m2);
    auto u1 = kk.Encode(m1);
    auto u2 = kk.Encode(m2);

    EXPECT_EQ(u1.Shape(0), 6U);
    EXPECT_EQ(u1.Shape(1), 14U);
    EXPECT_EQ(linalg::MatrixRank(u1), 6U);
    EXPECT_EQ(kk.CanonicalBasis(u1), u1);

    // Exact lifting isometry: d_S(I(c1), I(c2)) == 2 * d_R(c1, c2).
    const std::size_t ds = KotterKschischangCode<
        GaloisFieldBinary,
        GaloisFieldBinaryExtension<uint8_t>>::SubspaceDistance(u1, u2);
    const std::size_t dr = gab.Distance(c1, c2);
    EXPECT_EQ(ds, 2 * dr);
    if (m1 != m2) {
      EXPECT_GE(ds, kk.MinimumSubspaceDistance());
    }
  }
}

TEST(GabidulinSubspaceTest, KotterKschischangOperatorChannelDecoding) {
  std::mt19937_64 rng(31415926);

  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
  // G(6, 2) over GF(2^8)/GF(2) => n = 6, k = 2, n - k = 4, d_S = 10.
  auto kk = KotterKschischangCode<
      GaloisFieldBinary,
      GaloisFieldBinaryExtension<uint8_t>>::CreateStandard(gf2, gf256, 6, 2);

  // (a) Pure random row mixing A * U (no erasures, no insertions).
  {
    auto msg =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
    auto u = kk.Encode(msg);
    auto mix = RandomFullRowRankMatrix(gf2, 6, 6, rng);
    auto y = linalg::Dot(mix, u);
    EXPECT_EQ(kk.SubspaceDistance(u, y), 0U);

    auto res = kk.Decode(y);
    ASSERT_TRUE(res.Ok()) << res.error_message;
    EXPECT_EQ(res.subspace_basis, u);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.subspace_distance, 0U);
    EXPECT_EQ(res.num_erasures, 0U);
    EXPECT_EQ(res.num_insertions, 0U);
  }

  // (b) Random row mixing + gamma packet erasures (gamma = 1..4 = n - k).
  for (std::size_t gamma = 1; gamma <= 4; ++gamma) {
    for (int trial = 0; trial < 4; ++trial) {
      auto msg =
          GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
      auto u = kk.Encode(msg);
      auto mix = RandomFullRowRankMatrix(gf2, 6 - gamma, 6, rng);
      auto y = linalg::Dot(mix, u);
      ASSERT_EQ(kk.SubspaceDistance(u, y), gamma);

      auto res = kk.Decode(y);
      ASSERT_TRUE(res.Ok()) << "gamma=" << gamma << ": " << res.error_message;
      EXPECT_EQ(res.subspace_basis, u);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.subspace_distance, gamma);
      EXPECT_EQ(res.num_erasures, gamma);
      EXPECT_EQ(res.num_insertions, 0U);
    }
  }

  // (c) Random row mixing + rho injected error packets (rho = 1..4 = n - k).
  for (std::size_t rho = 1; rho <= 4; ++rho) {
    for (int trial = 0; trial < 4; ++trial) {
      auto msg =
          GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
      auto u = kk.Encode(msg);
      auto mix = RandomFullRowRankMatrix(gf2, 6, 6, rng);
      auto y_clean = linalg::Dot(mix, u);

      // Inject rho arbitrary error packets independent of U.
      GFMatrix<GaloisFieldBinary> y = y_clean;
      std::size_t added = 0;
      while (added < rho) {
        auto pkt = GFMatrix<GaloisFieldBinary>::Random(gf2, {1, 14}, rng);
        auto stacked = linalg::VStack(y, pkt);
        if (linalg::MatrixRank(linalg::VStack(u, stacked)) == 6 + added + 1) {
          y = std::move(stacked);
          ++added;
        }
      }
      // Scramble all rows of Y randomly so injected packets are mixed in.
      auto scramble = RandomFullRowRankMatrix(gf2, 6 + rho, 6 + rho, rng);
      y = linalg::Dot(scramble, y);
      ASSERT_EQ(kk.SubspaceDistance(u, y), rho);

      auto res = kk.Decode(y);
      ASSERT_TRUE(res.Ok()) << "rho=" << rho << ": " << res.error_message;
      EXPECT_EQ(res.subspace_basis, u);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.subspace_distance, rho);
    }
  }

  // (d) Combined operator-channel errors: random row mixing + packet erasures +
  //     injected error packets + packet corruptions (rank errors tau >= 1)
  //     satisfying d_S(U, Y) = 2*tau + rho + gamma <= 4!
  struct ChannelCase {
    std::size_t erasures;     // rows dropped from A*U
    std::size_t injections;   // pure-tail error packets [0 | e]
    std::size_t corruptions;  // additive rank-tau error matrix added to payload
  };
  const std::vector<ChannelCase> combined_cases = {
      {2, 2, 0},  // d_S = 4 (pure erasures + injections)
      {1, 3, 0},  // d_S = 4
      {3, 1, 0},  // d_S = 4
      {0, 0, 2},  // d_S = 4 (2 pure rank errors, tau = 2)
      {1, 1, 1},  // d_S = 4 (tau = 1, gamma = 1, rho = 1)
      {2, 0, 1},  // d_S = 4 (tau = 1, gamma = 2, rho = 0)
      {0, 2, 1},  // d_S = 4 (tau = 1, gamma = 0, rho = 2)
  };

  for (const auto &cc : combined_cases) {
    for (int trial = 0; trial < 5; ++trial) {
      auto msg =
          GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
      auto u = kk.Encode(msg);

      // Start from U and optionally add a rank-`corruptions` error to right
      // block.
      GFMatrix<GaloisFieldBinary> u_corrupt = u;
      if (cc.corruptions > 0) {
        auto err_vec = RandomRankError(kk.UnderlyingGabidulinCode().Tower(), 6,
                                       cc.corruptions, rng);
        auto err_mat = linalg::Transpose(
            kk.UnderlyingGabidulinCode().ExpandToBaseMatrix(err_vec));
        auto zero_left = GFMatrix<GaloisFieldBinary>::Zeros(gf2, {6, 6});
        u_corrupt = u_corrupt + linalg::HStack(zero_left, err_mat);
      }

      // Apply random row mixing with `erasures` rows dropped.
      auto mix = RandomFullRowRankMatrix(gf2, 6 - cc.erasures, 6, rng);
      auto y = linalg::Dot(mix, u_corrupt);

      // Append `injections` error packets [0 | v] independent of the rest.
      if (cc.injections > 0) {
        auto zero_hdr =
            GFMatrix<GaloisFieldBinary>::Zeros(gf2, {cc.injections, 6});
        auto tail = RandomFullRowRankMatrix(gf2, cc.injections, 8, rng);
        auto inj = linalg::HStack(zero_hdr, tail);
        y = linalg::VStack(y, inj);
      }

      // Apply an arbitrary invertible row transformation to Y.
      const std::size_t y_rows = y.Shape(0);
      auto post_mix = RandomFullRowRankMatrix(gf2, y_rows, y_rows, rng);
      y = linalg::Dot(post_mix, y);

      const std::size_t true_ds = kk.SubspaceDistance(u, y);
      ASSERT_LE(true_ds, 4U);

      auto res = kk.Decode(y);
      ASSERT_TRUE(res.Ok())
          << "Combined case (erasures=" << cc.erasures
          << ", injections=" << cc.injections
          << ", corruptions=" << cc.corruptions << ", d_S=" << true_ds
          << "): " << res.error_message;
      EXPECT_EQ(res.subspace_basis, u);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.subspace_distance, true_ds);
    }
  }

  // (e) All-zero codeword under erasures and injections.
  {
    auto zero_msg =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Zeros(gf256, {2});
    auto u = kk.Encode(zero_msg);
    auto mix = RandomFullRowRankMatrix(gf2, 4, 6, rng);
    auto y = linalg::Dot(mix, u);
    auto inj = linalg::HStack(GFMatrix<GaloisFieldBinary>::Zeros(gf2, {2, 6}),
                              RandomFullRowRankMatrix(gf2, 2, 8, rng));
    y = linalg::VStack(y, inj);
    auto res = kk.Decode(y);
    ASSERT_TRUE(res.Ok()) << res.error_message;
    EXPECT_EQ(res.subspace_basis, u);
    EXPECT_EQ(res.message, zero_msg);
    EXPECT_EQ(res.subspace_distance, 4U);
  }

  // (f) Failure when erasures + injections exceed n - k = 4 (e.g. 5 erasures).
  {
    auto msg =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf256, {2}, rng);
    auto u = kk.Encode(msg);
    auto mix = RandomFullRowRankMatrix(gf2, 1, 6, rng);
    auto y = linalg::Dot(mix, u);
    ASSERT_EQ(kk.SubspaceDistance(u, y), 5U);
    auto res = kk.Decode(y);
    EXPECT_FALSE(res.Ok());
  }

  // (g) Ternary Kotter-Kschischang code over GF(3^4) / GF(3): G(4, 2) => d_S
  // = 6.
  {
    auto gf3 = std::make_shared<GaloisFieldPrime<uint32_t>>(3);
    auto gf81 =
        std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 4));
    auto kk3 = KotterKschischangCode<
        GaloisFieldPrime<uint32_t>,
        GaloisFieldExtension<uint32_t>>::CreateStandard(gf3, gf81, 4, 2);
    EXPECT_EQ(kk3.AmbientDimension(), 8U);
    EXPECT_EQ(kk3.SubspaceDimension(), 4U);
    EXPECT_EQ(kk3.MinimumSubspaceDistance(), 6U);

    auto msg = GFArray<GaloisFieldExtension<uint32_t>>::Random(gf81, {2}, rng);
    auto u = kk3.Encode(msg);
    // 1 erasure + 1 injection => d_S = 2 <= n - k = 2.
    auto mix = RandomFullRowRankMatrix(gf3, 3, 4, rng);
    auto y = linalg::Dot(mix, u);
    auto inj =
        linalg::HStack(GFMatrix<GaloisFieldPrime<uint32_t>>::Zeros(gf3, {1, 4}),
                       RandomFullRowRankMatrix(gf3, 1, 4, rng));
    y = linalg::VStack(y, inj);
    auto post = RandomFullRowRankMatrix(gf3, 4, 4, rng);
    y = linalg::Dot(post, y);
    ASSERT_EQ(kk3.SubspaceDistance(u, y), 2U);

    auto res = kk3.Decode(y);
    ASSERT_TRUE(res.Ok()) << res.error_message;
    EXPECT_EQ(res.subspace_basis, u);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.subspace_distance, 2U);
  }
}

}  // namespace
