#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

namespace {

template <typename Field>
GFArray<Field> EvalPolyAtMatrix(const PolynomialDense<Field> &p,
                                const GFArray<Field> &a) {
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

template <typename Field>
GFArray<Field> MakeGFMatrix(
    const std::shared_ptr<Field> &field, std::size_t rows, std::size_t cols,
    const std::vector<typename Field::element_type> &v) {
  typename GFArray<Field>::raw_array_type raw =
      GFArray<Field>::raw_array_type::from_shape({rows, cols});
  for (std::size_t i = 0; i < rows * cols; ++i) {
    raw.data()[i] = v[i];
  }
  return GFArray<Field>(field, raw);
}

template <typename Field>
GFArray<Field> MakeGFVector(
    const std::shared_ptr<Field> &field,
    const std::vector<typename Field::element_type> &v) {
  typename GFArray<Field>::raw_array_type raw =
      GFArray<Field>::raw_array_type::from_shape({v.size()});
  for (std::size_t i = 0; i < v.size(); ++i) {
    raw.data()[i] = v[i];
  }
  return GFArray<Field>(field, raw);
}

}  // namespace

class LinalgPrimeTest : public ::testing::Test {
 protected:
  using FieldType = GaloisFieldPrime<uint32_t>;

  void SetUp() override { field = std::make_shared<FieldType>(7); }

  std::shared_ptr<FieldType> field;
};

TEST_F(LinalgPrimeTest, VectorDotProduct) {
  auto v1 = MakeGFVector(field, {2, 3, 4});
  auto v2 = MakeGFVector(field, {1, 5, 6});

  auto result = linalg::Dot(v1, v2);

  ASSERT_EQ(result.Dimension(), std::size_t{0});
  EXPECT_EQ(result.RawData()(), uint32_t{6});
}

TEST_F(LinalgPrimeTest, VectorDotProductMismatchedSizes) {
  auto v1 = linalg::Zeros(field, {3});
  auto v2 = linalg::Zeros(field, {4});

  EXPECT_THROW(linalg::Dot(v1, v2), std::invalid_argument);
}

TEST_F(LinalgPrimeTest, VectorMatrixMultiplication1D2D) {
  // [2, 5] * [[1, 2, 3], [4, 5, 6]] = [2+20, 4+25, 6+30] mod 7 = [1, 1, 1]
  auto v = MakeGFVector(field, {2, 5});
  auto a = MakeGFMatrix(field, 2, 3, {1, 2, 3, 4, 5, 6});
  auto result = linalg::Dot(v, a);

  ASSERT_EQ(result.Dimension(), std::size_t{1});
  ASSERT_EQ(result.Shape(0), std::size_t{3});
  EXPECT_EQ(result.At(0).Value(), uint32_t{1});
  EXPECT_EQ(result.At(1).Value(), uint32_t{1});
  EXPECT_EQ(result.At(2).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, MatrixVectorMultiplication) {
  auto a = MakeGFMatrix(field, 2, 3, {1, 2, 3, 4, 5, 6});
  auto v = MakeGFVector(field, {1, 2, 3});

  auto result = linalg::Dot(a, v);

  ASSERT_EQ(result.Dimension(), std::size_t{1});
  ASSERT_EQ(result.Shape(0), std::size_t{2});
  EXPECT_EQ(result.At(0).Value(), uint32_t{0});
  EXPECT_EQ(result.At(1).Value(), uint32_t{4});
}

TEST_F(LinalgPrimeTest, MatrixMultiplication) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 4});
  auto b = MakeGFMatrix(field, 2, 2, {5, 6, 1, 2});

  auto result = linalg::Dot(a, b);

  EXPECT_EQ(result.At(0, 0).Value(), uint32_t{0});
  EXPECT_EQ(result.At(0, 1).Value(), uint32_t{3});
  EXPECT_EQ(result.At(1, 0).Value(), uint32_t{5});
  EXPECT_EQ(result.At(1, 1).Value(), uint32_t{5});
}

TEST_F(LinalgPrimeTest, OuterProduct) {
  auto v1 = MakeGFVector(field, {2, 3});
  auto v2 = MakeGFVector(field, {4, 5, 6});

  auto result = linalg::Outer(v1, v2);

  EXPECT_EQ(result.At(0, 0).Value(), uint32_t{1});
  EXPECT_EQ(result.At(0, 1).Value(), uint32_t{3});
  EXPECT_EQ(result.At(0, 2).Value(), uint32_t{5});
  EXPECT_EQ(result.At(1, 0).Value(), uint32_t{5});
  EXPECT_EQ(result.At(1, 1).Value(), uint32_t{1});
  EXPECT_EQ(result.At(1, 2).Value(), uint32_t{4});
}

TEST_F(LinalgPrimeTest, Trace) {
  auto a = MakeGFMatrix(field, 3, 3, {1, 2, 3, 4, 5, 6, 0, 1, 2});

  auto result = linalg::Trace(a);

  EXPECT_EQ(result.Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, TraceNonSquareMatrix) {
  auto a = linalg::Zeros(field, {2, 3});

  EXPECT_THROW(linalg::Trace(a), std::invalid_argument);
}

TEST_F(LinalgPrimeTest, Determinant2x2) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 4});

  auto result = linalg::Det(a);

  EXPECT_EQ(result.Value(), uint32_t{5});
}

TEST_F(LinalgPrimeTest, Determinant3x3) {
  auto a = MakeGFMatrix(field, 3, 3, {1, 2, 3, 0, 1, 4, 5, 6, 0});

  auto result = linalg::Det(a);

  EXPECT_EQ(result.Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, DeterminantZeroPivotRowSwapSignFlip) {
  // [[0, 1], [1, 0]] has det = -1 == 6 mod 7
  auto p2 = MakeGFMatrix(field, 2, 2, {0, 1, 1, 0});
  EXPECT_EQ(linalg::Det(p2).Value(), uint32_t{6});

  // [[0, 1, 2], [1, 3, 4], [2, 5, 1]]:
  // det = -1*(1 - 8) + 2*(5 - 6) = 7 - 2 = 5 mod 7
  auto a3 = MakeGFMatrix(field, 3, 3, {0, 1, 2, 1, 3, 4, 2, 5, 1});
  EXPECT_EQ(linalg::Det(a3).Value(), uint32_t{5});

  // [[0, 0, 1], [0, 1, 0], [1, 0, 0]] has 1 swap => det = -1 == 6 mod 7
  auto p3 = MakeGFMatrix(field, 3, 3, {0, 0, 1, 0, 1, 0, 1, 0, 0});
  EXPECT_EQ(linalg::Det(p3).Value(), uint32_t{6});
}

TEST_F(LinalgPrimeTest, IdentityMatrix) {
  auto eye_mat = linalg::Eye(field, 3);

  EXPECT_EQ(eye_mat.At(0, 0).Value(), uint32_t{1});
  EXPECT_EQ(eye_mat.At(0, 1).Value(), uint32_t{0});
  EXPECT_EQ(eye_mat.At(0, 2).Value(), uint32_t{0});
  EXPECT_EQ(eye_mat.At(1, 0).Value(), uint32_t{0});
  EXPECT_EQ(eye_mat.At(1, 1).Value(), uint32_t{1});
  EXPECT_EQ(eye_mat.At(1, 2).Value(), uint32_t{0});
  EXPECT_EQ(eye_mat.At(2, 0).Value(), uint32_t{0});
  EXPECT_EQ(eye_mat.At(2, 1).Value(), uint32_t{0});
  EXPECT_EQ(eye_mat.At(2, 2).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, MatrixPower) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 4});

  auto a2 = linalg::MatrixPower(a, 2);

  EXPECT_EQ(a2.At(0, 0).Value(), uint32_t{0});
  EXPECT_EQ(a2.At(0, 1).Value(), uint32_t{3});
  EXPECT_EQ(a2.At(1, 0).Value(), uint32_t{1});
  EXPECT_EQ(a2.At(1, 1).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, MatrixPowerZero) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 4});

  auto a0 = linalg::MatrixPower(a, 0);

  EXPECT_EQ(a0.At(0, 0).Value(), uint32_t{1});
  EXPECT_EQ(a0.At(0, 1).Value(), uint32_t{0});
  EXPECT_EQ(a0.At(1, 0).Value(), uint32_t{0});
  EXPECT_EQ(a0.At(1, 1).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, MatrixPowerNegative) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 5});
  auto a_inv = linalg::Inv(a);
  auto a_neg1 = linalg::MatrixPower(a, -1);
  EXPECT_EQ(a_neg1, a_inv);

  auto a_pos5 = linalg::MatrixPower(a, 5);
  auto a_neg5 = linalg::MatrixPower(a, -5);
  auto prod = linalg::Dot(a_pos5, a_neg5);
  EXPECT_EQ(prod, linalg::Eye(field, 2));

  // Singular matrix negative power must throw.
  auto s = MakeGFMatrix(field, 2, 2, {1, 2, 2, 4});
  EXPECT_THROW(linalg::MatrixPower(s, -2), std::runtime_error);
}

TEST_F(LinalgPrimeTest, LinearSystemSolve) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 4});
  auto b = MakeGFVector(field, {5, 6});

  auto x = linalg::Solve(a, b);

  auto verification = linalg::Dot(a, x);
  EXPECT_EQ(verification.At(0).Value(), b.At(0).Value());
  EXPECT_EQ(verification.At(1).Value(), b.At(1).Value());
}

TEST_F(LinalgPrimeTest, LinearSystemSolveMultiRHSAndSingular) {
  auto a = MakeGFMatrix(field, 3, 3, {0, 1, 2, 1, 3, 4, 2, 5, 1});
  auto b = MakeGFMatrix(field, 3, 2, {1, 6, 2, 5, 3, 4});

  auto x = linalg::Solve(a, b);
  ASSERT_EQ(x.Dimension(), std::size_t{2});
  ASSERT_EQ(x.Shape(0), std::size_t{3});
  ASSERT_EQ(x.Shape(1), std::size_t{2});
  EXPECT_EQ(linalg::Dot(a, x), b);

  auto s = MakeGFMatrix(field, 2, 2, {1, 2, 2, 4});
  auto rhs = MakeGFVector(field, {3, 6});
  EXPECT_THROW(linalg::Solve(s, rhs), std::runtime_error);
}

TEST_F(LinalgPrimeTest, MatrixInversion) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 5});

  auto a_inv = linalg::Inv(a);
  auto verification = linalg::Dot(a, a_inv);

  EXPECT_EQ(verification.At(0, 0).Value(), uint32_t{1});
  EXPECT_EQ(verification.At(0, 1).Value(), uint32_t{0});
  EXPECT_EQ(verification.At(1, 0).Value(), uint32_t{0});
  EXPECT_EQ(verification.At(1, 1).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, SingularMatrixInversion) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 2, 4});

  EXPECT_THROW(linalg::Inv(a), std::runtime_error);
}

TEST_F(LinalgPrimeTest, VdotProduct) {
  auto a = MakeGFMatrix(field, 2, 2, {1, 2, 3, 4});
  auto b = MakeGFMatrix(field, 2, 2, {5, 6, 0, 1});

  auto result = linalg::Vdot(a, b);

  EXPECT_EQ(result.Value(), uint32_t{0});
}

TEST_F(LinalgPrimeTest, KroneckerProduct) {
  auto a = MakeGFMatrix(field, 1, 2, {1, 2});
  auto b = MakeGFMatrix(field, 2, 2, {3, 4, 5, 6});

  auto result = linalg::Kron(a, b);

  EXPECT_EQ(result.Shape(0), std::size_t{2});
  EXPECT_EQ(result.Shape(1), std::size_t{4});

  EXPECT_EQ(result.At(0, 0).Value(), uint32_t{3});
  EXPECT_EQ(result.At(0, 1).Value(), uint32_t{4});
  EXPECT_EQ(result.At(0, 2).Value(), uint32_t{6});
  EXPECT_EQ(result.At(0, 3).Value(), uint32_t{1});

  EXPECT_EQ(result.At(1, 0).Value(), uint32_t{5});
  EXPECT_EQ(result.At(1, 1).Value(), uint32_t{6});
  EXPECT_EQ(result.At(1, 2).Value(), uint32_t{3});
  EXPECT_EQ(result.At(1, 3).Value(), uint32_t{5});
}

TEST_F(LinalgPrimeTest, RowEchelonForm) {
  auto a = MakeGFMatrix(field, 3, 4, {1, 2, 1, 0, 2, 4, 3, 1, 3, 6, 0, 5});

  auto r = linalg::RowEchelon(a);

  EXPECT_EQ(r.At(0, 0).Value(), uint32_t{1});
  EXPECT_EQ(r.At(0, 1).Value(), uint32_t{2});
  EXPECT_EQ(r.At(0, 2).Value(), uint32_t{1});
  EXPECT_EQ(r.At(0, 3).Value(), uint32_t{0});

  EXPECT_EQ(r.At(1, 0).Value(), uint32_t{0});
  EXPECT_EQ(r.At(1, 1).Value(), uint32_t{0});
  EXPECT_EQ(r.At(1, 2).Value(), uint32_t{1});
  EXPECT_EQ(r.At(1, 3).Value(), uint32_t{1});

  EXPECT_EQ(r.At(2, 0).Value(), uint32_t{0});
  EXPECT_EQ(r.At(2, 1).Value(), uint32_t{0});
  EXPECT_EQ(r.At(2, 2).Value(), uint32_t{0});
  EXPECT_EQ(r.At(2, 3).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, MatrixRank) {
  auto a = MakeGFMatrix(field, 3, 4, {1, 2, 1, 0, 2, 4, 3, 1, 3, 6, 0, 5});
  EXPECT_EQ(linalg::MatrixRank(a), std::size_t{3});

  auto b = MakeGFMatrix(field, 2, 3, {1, 2, 3, 2, 4, 6});
  EXPECT_EQ(linalg::MatrixRank(b), std::size_t{1});
}

TEST_F(LinalgPrimeTest, CrossProduct3D) {
  auto v1 = MakeGFVector(field, {1, 2, 3});
  auto v2 = MakeGFVector(field, {4, 5, 6});

  auto result = linalg::Cross(v1, v2);

  EXPECT_EQ(result.Shape(0), std::size_t{3});
  EXPECT_EQ(result.At(0).Value(), uint32_t{4});
  EXPECT_EQ(result.At(1).Value(), uint32_t{6});
  EXPECT_EQ(result.At(2).Value(), uint32_t{4});
}

TEST_F(LinalgPrimeTest, CrossProduct2D) {
  auto v1 = MakeGFVector(field, {1, 2});
  auto v2 = MakeGFVector(field, {3, 4});

  auto result = linalg::Cross(v1, v2);

  EXPECT_EQ(result.Shape(0), std::size_t{3});
  EXPECT_EQ(result.At(0).Value(), uint32_t{0});
  EXPECT_EQ(result.At(1).Value(), uint32_t{0});
  EXPECT_EQ(result.At(2).Value(), uint32_t{5});
}

TEST_F(LinalgPrimeTest, ZerosMatrix) {
  auto z = linalg::Zeros(field, {2, 3});

  EXPECT_EQ(z.Shape(0), std::size_t{2});
  EXPECT_EQ(z.Shape(1), std::size_t{3});

  EXPECT_EQ(z.At(0, 0).Value(), uint32_t{0});
  EXPECT_EQ(z.At(0, 1).Value(), uint32_t{0});
  EXPECT_EQ(z.At(0, 2).Value(), uint32_t{0});
  EXPECT_EQ(z.At(1, 0).Value(), uint32_t{0});
  EXPECT_EQ(z.At(1, 1).Value(), uint32_t{0});
  EXPECT_EQ(z.At(1, 2).Value(), uint32_t{0});
}

TEST_F(LinalgPrimeTest, ZerosVector) {
  auto z = linalg::Zeros(field, {4});

  EXPECT_EQ(z.Shape(0), std::size_t{4});
  EXPECT_EQ(z.Dimension(), std::size_t{1});

  EXPECT_EQ(z.At(0).Value(), uint32_t{0});
  EXPECT_EQ(z.At(1).Value(), uint32_t{0});
  EXPECT_EQ(z.At(2).Value(), uint32_t{0});
  EXPECT_EQ(z.At(3).Value(), uint32_t{0});
}

TEST_F(LinalgPrimeTest, ReducedRowEchelonForm) {
  auto a = MakeGFMatrix(field, 3, 4, {1, 2, 1, 0, 2, 4, 3, 1, 3, 6, 0, 5});

  auto r = linalg::Rref(a);

  // Verify ALL 4 columns of RREF:
  // [[1, 2, 0, 0],
  //  [0, 0, 1, 0],
  //  [0, 0, 0, 1]]
  EXPECT_EQ(r.At(0, 0).Value(), uint32_t{1});
  EXPECT_EQ(r.At(0, 1).Value(), uint32_t{2});
  EXPECT_EQ(r.At(0, 2).Value(), uint32_t{0});
  EXPECT_EQ(r.At(0, 3).Value(), uint32_t{0});

  EXPECT_EQ(r.At(1, 0).Value(), uint32_t{0});
  EXPECT_EQ(r.At(1, 1).Value(), uint32_t{0});
  EXPECT_EQ(r.At(1, 2).Value(), uint32_t{1});
  EXPECT_EQ(r.At(1, 3).Value(), uint32_t{0});

  EXPECT_EQ(r.At(2, 0).Value(), uint32_t{0});
  EXPECT_EQ(r.At(2, 1).Value(), uint32_t{0});
  EXPECT_EQ(r.At(2, 2).Value(), uint32_t{0});
  EXPECT_EQ(r.At(2, 3).Value(), uint32_t{1});
}

TEST_F(LinalgPrimeTest, NullspaceAndFundamentalSubspaces) {
  // Rank-2 matrix of shape 3x5 over GF(7)
  auto a =
      MakeGFMatrix(field, 3, 5, {1, 2, 3, 4, 5, 2, 4, 6, 1, 3, 0, 1, 2, 3, 4});
  const std::size_t rank = linalg::MatrixRank(a);
  EXPECT_EQ(rank, std::size_t{2});

  auto ns = linalg::Nullspace(a);
  EXPECT_EQ(ns.Shape(0), 5 - rank);
  EXPECT_EQ(ns.Shape(1), std::size_t{5});
  EXPECT_EQ(linalg::MatrixRank(ns), 5 - rank);
  // A * N^T == 0
  auto prod = linalg::Dot(a, linalg::Transpose(ns));
  EXPECT_EQ(prod, linalg::Zeros(field, {3, 5 - rank}));

  // Left nullspace: L * A == 0
  auto lns = linalg::LeftNullspace(a);
  EXPECT_EQ(lns.Shape(0), 3 - rank);
  EXPECT_EQ(lns.Shape(1), std::size_t{3});
  EXPECT_EQ(linalg::Dot(lns, a), linalg::Zeros(field, {3 - rank, 5}));

  // Row space & column space dimensions
  auto rs = linalg::RowSpace(a);
  EXPECT_EQ(rs.Shape(0), rank);
  EXPECT_EQ(rs.Shape(1), std::size_t{5});
  EXPECT_EQ(linalg::MatrixRank(rs), rank);

  auto cs = linalg::ColumnSpace(a);
  EXPECT_EQ(cs.Shape(0), std::size_t{3});
  EXPECT_EQ(cs.Shape(1), rank);
  EXPECT_EQ(linalg::MatrixRank(cs), rank);

  // Full column rank matrix -> empty nullspace of shape {0, n}
  auto full_col = MakeGFMatrix(field, 3, 2, {1, 0, 0, 1, 2, 3});
  auto ns_empty = linalg::Nullspace(full_col);
  EXPECT_EQ(ns_empty.Shape(0), std::size_t{0});
  EXPECT_EQ(ns_empty.Shape(1), std::size_t{2});
}

TEST_F(LinalgPrimeTest, SystematicForm) {
  // Generator matrix with non-leading pivot columns {1, 3}
  auto g = MakeGFMatrix(field, 2, 5, {0, 1, 2, 0, 3, 0, 0, 0, 2, 4});
  auto res = linalg::SystematicForm(g);

  EXPECT_EQ(res.pivot_cols, (std::vector<std::size_t>{1, 3}));
  EXPECT_EQ(res.free_cols, (std::vector<std::size_t>{0, 2, 4}));
  EXPECT_EQ(res.col_perm, (std::vector<std::size_t>{1, 3, 0, 2, 4}));
  ASSERT_EQ(res.g_sys.Shape(0), std::size_t{2});
  ASSERT_EQ(res.g_sys.Shape(1), std::size_t{5});

  // Left block must be I_2
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      EXPECT_EQ(res.g_sys.At(i, j).Value(), i == j ? uint32_t{1} : uint32_t{0});
    }
  }

  // Construct parity check matrix H in original coordinates and check G * H^T =
  // 0
  const std::size_t k = res.pivot_cols.size();
  const std::size_t n = g.Shape(1);
  auto h = linalg::Zeros(field, {n - k, n});
  for (std::size_t r = 0; r < n - k; ++r) {
    h.Raw(r, res.free_cols[r]) = field->MultiplicativeIdentity();
    for (std::size_t c = 0; c < k; ++c) {
      h.Raw(r, res.pivot_cols[c]) = field->Neg(res.g_sys.Raw(c, k + r));
    }
  }
  EXPECT_EQ(linalg::Dot(g, linalg::Transpose(h)),
            linalg::Zeros(field, {2, n - k}));
}

TEST_F(LinalgPrimeTest, PluDecomposition) {
  // Matrix requiring row swaps (A(0,0) == 0)
  auto a = MakeGFMatrix(field, 4, 4,
                        {0, 2, 1, 3, 1, 4, 2, 5, 2, 1, 0, 6, 3, 5, 4, 2});
  auto res = linalg::Plu(a);

  // Verify P * A == L * U
  EXPECT_EQ(linalg::Dot(res.p, a), linalg::Dot(res.l, res.u));

  // Verify L is unit lower triangular and U is upper triangular
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(res.l.At(i, i).Value(), uint32_t{1});
    for (std::size_t j = i + 1; j < 4; ++j) {
      EXPECT_EQ(res.l.At(i, j).Value(), uint32_t{0});
    }
    for (std::size_t j = 0; j < i; ++j) {
      EXPECT_EQ(res.u.At(i, j).Value(), uint32_t{0});
    }
  }

  // Verify det(A) == (-1)^num_swaps * prod(diag(U))
  uint32_t diag_prod = 1;
  for (std::size_t i = 0; i < 4; ++i) {
    diag_prod = field->Mul(diag_prod, res.u.Raw(i, i));
  }
  if (res.num_swaps % 2 == 1) {
    diag_prod = field->Neg(diag_prod);
  }
  EXPECT_EQ(diag_prod, linalg::Det(a).Value());

  // Rectangular matrix PLU
  auto rect = MakeGFMatrix(field, 3, 4, {0, 1, 2, 3, 2, 3, 4, 5, 4, 6, 1, 2});
  auto rect_plu = linalg::Plu(rect);
  EXPECT_EQ(linalg::Dot(rect_plu.p, rect), linalg::Dot(rect_plu.l, rect_plu.u));
}

TEST_F(LinalgPrimeTest, CharacteristicPolyCayleyHamilton) {
  std::mt19937_64 rng(42);
  for (std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{3},
                        std::size_t{4}, std::size_t{5}}) {
    auto a = linalg::RandomMatrix(field, {n, n}, rng);
    auto p = linalg::CharacteristicPoly(a);

    EXPECT_EQ(p.Degree(), static_cast<int>(n));
    EXPECT_EQ(p[n].Value(), uint32_t{1});

    // Degree n-1 coefficient == -trace(A)
    auto neg_tr = -linalg::Trace(a);
    EXPECT_EQ(p[n - 1], neg_tr);

    // Constant term == (-1)^n * det(A)
    auto det_a = linalg::Det(a);
    auto expected_c0 = (n % 2 == 1) ? -det_a : det_a;
    EXPECT_EQ(p[0], expected_c0);

    // Cayley-Hamilton: p_A(A) == 0
    auto p_a = EvalPolyAtMatrix(p, a);
    EXPECT_EQ(p_a, linalg::Zeros(field, {n, n}));
  }
}

TEST_F(LinalgPrimeTest, GFArrayDirectAPIUtilities) {
  auto a = MakeGFMatrix(field, 2, 3, {1, 2, 3, 4, 5, 6});
  auto b = MakeGFMatrix(field, 2, 2, {6, 5, 4, 3});
  auto c = MakeGFMatrix(field, 1, 3, {2, 0, 1});

  auto a_t = linalg::Transpose(a);
  EXPECT_EQ(a_t.Shape(0), std::size_t{3});
  EXPECT_EQ(a_t.Shape(1), std::size_t{2});
  EXPECT_EQ(linalg::Transpose(a_t), a);

  auto h = linalg::HStack(a, b);
  EXPECT_EQ(h.Shape(0), std::size_t{2});
  EXPECT_EQ(h.Shape(1), std::size_t{5});
  EXPECT_EQ(h.At(0, 3).Value(), uint32_t{6});
  EXPECT_EQ(h.At(1, 4).Value(), uint32_t{3});

  auto v = linalg::VStack(a, c);
  EXPECT_EQ(v.Shape(0), std::size_t{3});
  EXPECT_EQ(v.Shape(1), std::size_t{3});
  EXPECT_EQ(v.At(2, 0).Value(), uint32_t{2});
  EXPECT_EQ(v.At(2, 2).Value(), uint32_t{1});

  auto ones_mat = linalg::Ones(field, {2, 2});
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      EXPECT_EQ(ones_mat.At(i, j).Value(), uint32_t{1});
    }
  }
}

TEST_F(LinalgPrimeTest, EmptyArrayAndNullFieldErrorHandling) {
  std::shared_ptr<FieldType> null_field;
  EXPECT_THROW(linalg::Eye(null_field, 2), std::invalid_argument);
  EXPECT_THROW(linalg::Zeros(null_field, {2, 2}), std::invalid_argument);
  EXPECT_THROW(linalg::Ones(null_field, {2, 2}), std::invalid_argument);
  EXPECT_THROW(linalg::RandomMatrix(null_field, {2, 2}), std::invalid_argument);

  // Default-constructed GFArray has null field pointer.
  GFArray<FieldType> empty_gf;
  EXPECT_THROW(linalg::Det(empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::Inv(empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::Rref(empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::Dot(empty_gf, empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::CharacteristicPoly(empty_gf), std::invalid_argument);

  // Mismatched fields.
  auto other_field = std::make_shared<FieldType>(11);
  auto a7 = linalg::Eye(field, 2);
  auto a11 = linalg::Eye(other_field, 2);
  EXPECT_THROW(linalg::Dot(a7, a11), std::invalid_argument);
  EXPECT_THROW(linalg::Solve(a7, linalg::Ones(other_field, {2})),
               std::invalid_argument);
}
