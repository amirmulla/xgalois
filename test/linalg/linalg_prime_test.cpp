#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

namespace {

template <typename Field>
GFArray<Field> EvalPolyAtMatrix(const PolynomialDense<Field> &p,
                                const GFArray<Field> &A) {
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
  using FieldElement = GaloisFieldElement<FieldType>;

  void SetUp() override { field = std::make_shared<FieldType>(7); }

  std::shared_ptr<FieldType> field;
};

TEST_F(LinalgPrimeTest, VectorDotProduct) {
  xt::xarray<FieldElement> v1({3});
  xt::xarray<FieldElement> v2({3});

  v1(0) = FieldElement(2, field);
  v1(1) = FieldElement(3, field);
  v1(2) = FieldElement(4, field);

  v2(0) = FieldElement(1, field);
  v2(1) = FieldElement(5, field);
  v2(2) = FieldElement(6, field);

  auto result = linalg::dot(v1, v2);

  EXPECT_EQ(result().Value(), 6);
}

TEST_F(LinalgPrimeTest, VectorDotProductMismatchedSizes) {
  auto v1 = linalg::zeros({3}, field);
  auto v2 = linalg::zeros({4}, field);

  EXPECT_THROW(linalg::dot(v1, v2), std::invalid_argument);
}

TEST_F(LinalgPrimeTest, VectorMatrixMultiplication1D2D) {
  xt::xarray<FieldElement> v({2});
  xt::xarray<FieldElement> A({2, 3});

  v(0) = FieldElement(2, field);
  v(1) = FieldElement(5, field);

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(3, field);
  A(1, 0) = FieldElement(4, field);
  A(1, 1) = FieldElement(5, field);
  A(1, 2) = FieldElement(6, field);

  // [2, 5] * [[1, 2, 3], [4, 5, 6]] = [2+20, 4+25, 6+30] mod 7 = [1, 1, 1]
  auto result = linalg::dot(v, A);
  ASSERT_EQ(result.dimension(), 1);
  ASSERT_EQ(result.shape(0), 3);
  EXPECT_EQ(result(0).Value(), 1);
  EXPECT_EQ(result(1).Value(), 1);
  EXPECT_EQ(result(2).Value(), 1);

  auto gv = MakeGFVector(field, {2, 5});
  auto gA = MakeGFMatrix(field, 2, 3, {1, 2, 3, 4, 5, 6});
  auto gres = linalg::dot(gv, gA);
  ASSERT_EQ(gres.Dimension(), 1);
  ASSERT_EQ(gres.Shape(0), 3);
  EXPECT_EQ(gres.At(0).Value(), 1);
  EXPECT_EQ(gres.At(1).Value(), 1);
  EXPECT_EQ(gres.At(2).Value(), 1);
}

TEST_F(LinalgPrimeTest, MatrixVectorMultiplication) {
  xt::xarray<FieldElement> A({2, 3});
  xt::xarray<FieldElement> v({3});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(3, field);
  A(1, 0) = FieldElement(4, field);
  A(1, 1) = FieldElement(5, field);
  A(1, 2) = FieldElement(6, field);

  v(0) = FieldElement(1, field);
  v(1) = FieldElement(2, field);
  v(2) = FieldElement(3, field);

  auto result = linalg::dot(A, v);

  EXPECT_EQ(result(0).Value(), 0);
  EXPECT_EQ(result(1).Value(), 4);
}

TEST_F(LinalgPrimeTest, MatrixMultiplication) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  B(0, 0) = FieldElement(5, field);
  B(0, 1) = FieldElement(6, field);
  B(1, 0) = FieldElement(1, field);
  B(1, 1) = FieldElement(2, field);

  auto result = linalg::dot(A, B);

  EXPECT_EQ(result(0, 0).Value(), 0);
  EXPECT_EQ(result(0, 1).Value(), 3);
  EXPECT_EQ(result(1, 0).Value(), 5);
  EXPECT_EQ(result(1, 1).Value(), 5);
}

TEST_F(LinalgPrimeTest, OuterProduct) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({3});

  v1(0) = FieldElement(2, field);
  v1(1) = FieldElement(3, field);

  v2(0) = FieldElement(4, field);
  v2(1) = FieldElement(5, field);
  v2(2) = FieldElement(6, field);

  auto result = linalg::outer(v1, v2);

  EXPECT_EQ(result(0, 0).Value(), 1);
  EXPECT_EQ(result(0, 1).Value(), 3);
  EXPECT_EQ(result(0, 2).Value(), 5);
  EXPECT_EQ(result(1, 0).Value(), 5);
  EXPECT_EQ(result(1, 1).Value(), 1);
  EXPECT_EQ(result(1, 2).Value(), 4);
}

TEST_F(LinalgPrimeTest, Trace) {
  xt::xarray<FieldElement> A({3, 3});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(3, field);
  A(1, 0) = FieldElement(4, field);
  A(1, 1) = FieldElement(5, field);
  A(1, 2) = FieldElement(6, field);
  A(2, 0) = FieldElement(0, field);
  A(2, 1) = FieldElement(1, field);
  A(2, 2) = FieldElement(2, field);

  auto result = linalg::trace(A);

  EXPECT_EQ(result.Value(), 1);
}

TEST_F(LinalgPrimeTest, TraceNonSquareMatrix) {
  auto A = linalg::zeros({2, 3}, field);

  EXPECT_THROW(linalg::trace(A), std::invalid_argument);
}

TEST_F(LinalgPrimeTest, Determinant2x2) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  auto result = linalg::det(A);

  EXPECT_EQ(result.Value(), 5);
}

TEST_F(LinalgPrimeTest, Determinant3x3) {
  xt::xarray<FieldElement> A({3, 3});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(3, field);
  A(1, 0) = FieldElement(0, field);
  A(1, 1) = FieldElement(1, field);
  A(1, 2) = FieldElement(4, field);
  A(2, 0) = FieldElement(5, field);
  A(2, 1) = FieldElement(6, field);
  A(2, 2) = FieldElement(0, field);

  auto result = linalg::det(A);

  EXPECT_EQ(result.Value(), 1);
}

TEST_F(LinalgPrimeTest, DeterminantZeroPivotRowSwapSignFlip) {
  // [[0, 1], [1, 0]] has det = -1 == 6 mod 7
  auto P2 = MakeGFMatrix(field, 2, 2, {0, 1, 1, 0});
  EXPECT_EQ(linalg::det(P2).Value(), 6);

  // [[0, 1, 2], [1, 3, 4], [2, 5, 1]]:
  // det = -1*(1 - 8) + 2*(5 - 6) = 7 - 2 = 5 mod 7
  auto A3 = MakeGFMatrix(field, 3, 3, {0, 1, 2, 1, 3, 4, 2, 5, 1});
  EXPECT_EQ(linalg::det(A3).Value(), 5);

  // [[0, 0, 1], [0, 1, 0], [1, 0, 0]] has 1 swap => det = -1 == 6 mod 7
  auto P3 = MakeGFMatrix(field, 3, 3, {0, 0, 1, 0, 1, 0, 1, 0, 0});
  EXPECT_EQ(linalg::det(P3).Value(), 6);
}

TEST_F(LinalgPrimeTest, IdentityMatrix) {
  auto I = linalg::eye(3, field);

  EXPECT_EQ(I(0, 0).Value(), 1);
  EXPECT_EQ(I(0, 1).Value(), 0);
  EXPECT_EQ(I(0, 2).Value(), 0);
  EXPECT_EQ(I(1, 0).Value(), 0);
  EXPECT_EQ(I(1, 1).Value(), 1);
  EXPECT_EQ(I(1, 2).Value(), 0);
  EXPECT_EQ(I(2, 0).Value(), 0);
  EXPECT_EQ(I(2, 1).Value(), 0);
  EXPECT_EQ(I(2, 2).Value(), 1);
}

TEST_F(LinalgPrimeTest, MatrixPower) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  auto A2 = linalg::matrix_power(A, 2);

  EXPECT_EQ(A2(0, 0).Value(), 0);
  EXPECT_EQ(A2(0, 1).Value(), 3);
  EXPECT_EQ(A2(1, 0).Value(), 1);
  EXPECT_EQ(A2(1, 1).Value(), 1);
}

TEST_F(LinalgPrimeTest, MatrixPowerZero) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  auto A0 = linalg::matrix_power(A, 0);

  EXPECT_EQ(A0(0, 0).Value(), 1);
  EXPECT_EQ(A0(0, 1).Value(), 0);
  EXPECT_EQ(A0(1, 0).Value(), 0);
  EXPECT_EQ(A0(1, 1).Value(), 1);
}

TEST_F(LinalgPrimeTest, MatrixPowerNegative) {
  auto A = MakeGFMatrix(field, 2, 2, {1, 2, 3, 5});
  auto A_inv = linalg::inv(A);
  auto A_neg1 = linalg::matrix_power(A, -1);
  EXPECT_EQ(A_neg1, A_inv);

  auto A_pos5 = linalg::matrix_power(A, 5);
  auto A_neg5 = linalg::matrix_power(A, -5);
  auto prod = linalg::dot(A_pos5, A_neg5);
  EXPECT_EQ(prod, linalg::eye(field, 2));

  // Singular matrix negative power must throw.
  auto S = MakeGFMatrix(field, 2, 2, {1, 2, 2, 4});
  EXPECT_THROW(linalg::matrix_power(S, -2), std::runtime_error);
}

TEST_F(LinalgPrimeTest, LinearSystemSolve) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> b({2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  b(0) = FieldElement(5, field);
  b(1) = FieldElement(6, field);

  auto x = linalg::solve(A, b);

  auto verification = linalg::dot(A, x);
  EXPECT_EQ(verification(0).Value(), b(0).Value());
  EXPECT_EQ(verification(1).Value(), b(1).Value());
}

TEST_F(LinalgPrimeTest, LinearSystemSolveMultiRHSAndSingular) {
  auto A = MakeGFMatrix(field, 3, 3, {0, 1, 2, 1, 3, 4, 2, 5, 1});
  auto B = MakeGFMatrix(field, 3, 2, {1, 6, 2, 5, 3, 4});

  auto X = linalg::solve(A, B);
  ASSERT_EQ(X.Dimension(), 2);
  ASSERT_EQ(X.Shape(0), 3);
  ASSERT_EQ(X.Shape(1), 2);
  EXPECT_EQ(linalg::dot(A, X), B);

  auto S = MakeGFMatrix(field, 2, 2, {1, 2, 2, 4});
  auto b = MakeGFVector(field, {3, 6});
  EXPECT_THROW(linalg::solve(S, b), std::runtime_error);
}

TEST_F(LinalgPrimeTest, MatrixInversion) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(5, field);

  auto A_inv = linalg::inv(A);

  auto verification = linalg::dot(A, A_inv);

  EXPECT_EQ(verification(0, 0).Value(), 1);
  EXPECT_EQ(verification(0, 1).Value(), 0);
  EXPECT_EQ(verification(1, 0).Value(), 0);
  EXPECT_EQ(verification(1, 1).Value(), 1);
}

TEST_F(LinalgPrimeTest, SingularMatrixInversion) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(2, field);
  A(1, 1) = FieldElement(4, field);

  EXPECT_THROW(linalg::inv(A), std::runtime_error);
}

TEST_F(LinalgPrimeTest, VdotProduct) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  B(0, 0) = FieldElement(5, field);
  B(0, 1) = FieldElement(6, field);
  B(1, 0) = FieldElement(0, field);
  B(1, 1) = FieldElement(1, field);

  auto result = linalg::vdot(A, B);

  EXPECT_EQ(result.Value(), 0);
}

TEST_F(LinalgPrimeTest, KroneckerProduct) {
  xt::xarray<FieldElement> A({1, 2});
  xt::xarray<FieldElement> B({2, 2});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);

  B(0, 0) = FieldElement(3, field);
  B(0, 1) = FieldElement(4, field);
  B(1, 0) = FieldElement(5, field);
  B(1, 1) = FieldElement(6, field);

  auto result = linalg::kron(A, B);

  EXPECT_EQ(result.shape(0), 2);
  EXPECT_EQ(result.shape(1), 4);

  EXPECT_EQ(result(0, 0).Value(), 3);
  EXPECT_EQ(result(0, 1).Value(), 4);
  EXPECT_EQ(result(0, 2).Value(), 6);
  EXPECT_EQ(result(0, 3).Value(), 1);

  EXPECT_EQ(result(1, 0).Value(), 5);
  EXPECT_EQ(result(1, 1).Value(), 6);
  EXPECT_EQ(result(1, 2).Value(), 3);
  EXPECT_EQ(result(1, 3).Value(), 5);
}

TEST_F(LinalgPrimeTest, RowEchelonForm) {
  xt::xarray<FieldElement> A({3, 4});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(1, field);
  A(0, 3) = FieldElement(0, field);
  A(1, 0) = FieldElement(2, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(3, field);
  A(1, 3) = FieldElement(1, field);
  A(2, 0) = FieldElement(3, field);
  A(2, 1) = FieldElement(6, field);
  A(2, 2) = FieldElement(0, field);
  A(2, 3) = FieldElement(5, field);

  auto R = linalg::row_echelon(A);

  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(0, 1).Value(), 2);
  EXPECT_EQ(R(0, 2).Value(), 1);
  EXPECT_EQ(R(0, 3).Value(), 0);

  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(1, 1).Value(), 0);
  EXPECT_EQ(R(1, 2).Value(), 1);
  EXPECT_EQ(R(1, 3).Value(), 1);

  EXPECT_EQ(R(2, 0).Value(), 0);
  EXPECT_EQ(R(2, 1).Value(), 0);
  EXPECT_EQ(R(2, 2).Value(), 0);
  EXPECT_EQ(R(2, 3).Value(), 1);
}

TEST_F(LinalgPrimeTest, MatrixRank) {
  xt::xarray<FieldElement> A({3, 4});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(1, field);
  A(0, 3) = FieldElement(0, field);
  A(1, 0) = FieldElement(2, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(3, field);
  A(1, 3) = FieldElement(1, field);
  A(2, 0) = FieldElement(3, field);
  A(2, 1) = FieldElement(6, field);
  A(2, 2) = FieldElement(0, field);
  A(2, 3) = FieldElement(5, field);
  EXPECT_EQ(linalg::matrix_rank(A), 3);

  xt::xarray<FieldElement> B({2, 3});

  B(0, 0) = FieldElement(1, field);
  B(0, 1) = FieldElement(2, field);
  B(0, 2) = FieldElement(3, field);
  B(1, 0) = FieldElement(2, field);
  B(1, 1) = FieldElement(4, field);
  B(1, 2) = FieldElement(6, field);
  EXPECT_EQ(linalg::matrix_rank(B), 1);
}

TEST_F(LinalgPrimeTest, CrossProduct3D) {
  xt::xarray<FieldElement> v1({3});
  xt::xarray<FieldElement> v2({3});

  v1(0) = FieldElement(1, field);
  v1(1) = FieldElement(2, field);
  v1(2) = FieldElement(3, field);

  v2(0) = FieldElement(4, field);
  v2(1) = FieldElement(5, field);
  v2(2) = FieldElement(6, field);

  auto result = linalg::cross(v1, v2);

  EXPECT_EQ(result.shape(0), 3);
  EXPECT_EQ(result(0).Value(), 4);
  EXPECT_EQ(result(1).Value(), 6);
  EXPECT_EQ(result(2).Value(), 4);
}

TEST_F(LinalgPrimeTest, CrossProduct2D) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({2});

  v1(0) = FieldElement(1, field);
  v1(1) = FieldElement(2, field);

  v2(0) = FieldElement(3, field);
  v2(1) = FieldElement(4, field);

  auto result = linalg::cross(v1, v2);

  EXPECT_EQ(result.shape(0), 3);
  EXPECT_EQ(result(0).Value(), 0);
  EXPECT_EQ(result(1).Value(), 0);
  EXPECT_EQ(result(2).Value(), 5);
}

TEST_F(LinalgPrimeTest, ZerosMatrix) {
  auto Z = linalg::zeros({2, 3}, field);

  EXPECT_EQ(Z.shape(0), 2);
  EXPECT_EQ(Z.shape(1), 3);

  EXPECT_EQ(Z(0, 0).Value(), 0);
  EXPECT_EQ(Z(0, 1).Value(), 0);
  EXPECT_EQ(Z(0, 2).Value(), 0);
  EXPECT_EQ(Z(1, 0).Value(), 0);
  EXPECT_EQ(Z(1, 1).Value(), 0);
  EXPECT_EQ(Z(1, 2).Value(), 0);
}

TEST_F(LinalgPrimeTest, ZerosVector) {
  auto z = linalg::zeros({4}, field);

  EXPECT_EQ(z.shape(0), 4);
  EXPECT_EQ(z.dimension(), 1);

  EXPECT_EQ(z(0).Value(), 0);
  EXPECT_EQ(z(1).Value(), 0);
  EXPECT_EQ(z(2).Value(), 0);
  EXPECT_EQ(z(3).Value(), 0);
}

TEST_F(LinalgPrimeTest, ReducedRowEchelonForm) {
  xt::xarray<FieldElement> A({3, 4});

  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(1, field);
  A(0, 3) = FieldElement(0, field);
  A(1, 0) = FieldElement(2, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(3, field);
  A(1, 3) = FieldElement(1, field);
  A(2, 0) = FieldElement(3, field);
  A(2, 1) = FieldElement(6, field);
  A(2, 2) = FieldElement(0, field);
  A(2, 3) = FieldElement(5, field);

  auto R = linalg::rref(A);

  // Verify ALL 4 columns of RREF:
  // [[1, 2, 0, 0],
  //  [0, 0, 1, 0],
  //  [0, 0, 0, 1]]
  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(0, 1).Value(), 2);
  EXPECT_EQ(R(0, 2).Value(), 0);
  EXPECT_EQ(R(0, 3).Value(), 0);

  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(1, 1).Value(), 0);
  EXPECT_EQ(R(1, 2).Value(), 1);
  EXPECT_EQ(R(1, 3).Value(), 0);

  EXPECT_EQ(R(2, 0).Value(), 0);
  EXPECT_EQ(R(2, 1).Value(), 0);
  EXPECT_EQ(R(2, 2).Value(), 0);
  EXPECT_EQ(R(2, 3).Value(), 1);
}

TEST_F(LinalgPrimeTest, NullspaceAndFundamentalSubspaces) {
  // Rank-2 matrix of shape 3x5 over GF(7)
  auto A =
      MakeGFMatrix(field, 3, 5, {1, 2, 3, 4, 5, 2, 4, 6, 1, 3, 0, 1, 2, 3, 4});
  const std::size_t rank = linalg::matrix_rank(A);
  EXPECT_EQ(rank, 2);

  auto N = linalg::nullspace(A);
  EXPECT_EQ(N.Shape(0), 5 - rank);
  EXPECT_EQ(N.Shape(1), 5);
  EXPECT_EQ(linalg::matrix_rank(N), 5 - rank);
  // A * N^T == 0
  auto prod = linalg::dot(A, linalg::transpose(N));
  EXPECT_EQ(prod, linalg::zeros(field, {3, 5 - rank}));

  // Left nullspace: L * A == 0
  auto L = linalg::left_nullspace(A);
  EXPECT_EQ(L.Shape(0), 3 - rank);
  EXPECT_EQ(L.Shape(1), 3);
  EXPECT_EQ(linalg::dot(L, A), linalg::zeros(field, {3 - rank, 5}));

  // Row space & column space dimensions
  auto RS = linalg::row_space(A);
  EXPECT_EQ(RS.Shape(0), rank);
  EXPECT_EQ(RS.Shape(1), 5);
  EXPECT_EQ(linalg::matrix_rank(RS), rank);

  auto CS = linalg::column_space(A);
  EXPECT_EQ(CS.Shape(0), 3);
  EXPECT_EQ(CS.Shape(1), rank);
  EXPECT_EQ(linalg::matrix_rank(CS), rank);

  // Full column rank matrix -> empty nullspace of shape {0, n}
  auto FullCol = MakeGFMatrix(field, 3, 2, {1, 0, 0, 1, 2, 3});
  auto N_empty = linalg::nullspace(FullCol);
  EXPECT_EQ(N_empty.Shape(0), 0);
  EXPECT_EQ(N_empty.Shape(1), 2);
}

TEST_F(LinalgPrimeTest, SystematicForm) {
  // Generator matrix with non-leading pivot columns {1, 3}
  auto G = MakeGFMatrix(field, 2, 5, {0, 1, 2, 0, 3, 0, 0, 0, 2, 4});
  auto res = linalg::systematic_form(G);

  EXPECT_EQ(res.pivot_cols, (std::vector<std::size_t>{1, 3}));
  EXPECT_EQ(res.free_cols, (std::vector<std::size_t>{0, 2, 4}));
  EXPECT_EQ(res.col_perm, (std::vector<std::size_t>{1, 3, 0, 2, 4}));
  ASSERT_EQ(res.G_sys.Shape(0), 2);
  ASSERT_EQ(res.G_sys.Shape(1), 5);

  // Left block must be I_2
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      EXPECT_EQ(res.G_sys.At(i, j).Value(), i == j ? 1U : 0U);
    }
  }

  // Construct parity check matrix H in original coordinates and check G * H^T =
  // 0
  const std::size_t k = res.pivot_cols.size();
  const std::size_t n = G.Shape(1);
  auto H = linalg::zeros(field, {n - k, n});
  for (std::size_t r = 0; r < n - k; ++r) {
    H.Raw(r, res.free_cols[r]) = field->MultiplicativeIdentity();
    for (std::size_t c = 0; c < k; ++c) {
      H.Raw(r, res.pivot_cols[c]) = field->Neg(res.G_sys.Raw(c, k + r));
    }
  }
  EXPECT_EQ(linalg::dot(G, linalg::transpose(H)),
            linalg::zeros(field, {2, n - k}));
}

TEST_F(LinalgPrimeTest, PLUDecomposition) {
  // Matrix requiring row swaps (A(0,0) == 0)
  auto A = MakeGFMatrix(field, 4, 4,
                        {0, 2, 1, 3, 1, 4, 2, 5, 2, 1, 0, 6, 3, 5, 4, 2});
  auto res = linalg::plu(A);

  // Verify P * A == L * U
  EXPECT_EQ(linalg::dot(res.P, A), linalg::dot(res.L, res.U));

  // Verify L is unit lower triangular and U is upper triangular
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(res.L.At(i, i).Value(), 1);
    for (std::size_t j = i + 1; j < 4; ++j) {
      EXPECT_EQ(res.L.At(i, j).Value(), 0);
    }
    for (std::size_t j = 0; j < i; ++j) {
      EXPECT_EQ(res.U.At(i, j).Value(), 0);
    }
  }

  // Verify det(A) == (-1)^num_swaps * prod(diag(U))
  uint32_t diag_prod = 1;
  for (std::size_t i = 0; i < 4; ++i) {
    diag_prod = field->Mul(diag_prod, res.U.Raw(i, i));
  }
  if (res.num_swaps % 2 == 1) {
    diag_prod = field->Neg(diag_prod);
  }
  EXPECT_EQ(diag_prod, linalg::det(A).Value());

  // Rectangular matrix PLU
  auto Rect = MakeGFMatrix(field, 3, 4, {0, 1, 2, 3, 2, 3, 4, 5, 4, 6, 1, 2});
  auto rect_plu = linalg::plu(Rect);
  EXPECT_EQ(linalg::dot(rect_plu.P, Rect), linalg::dot(rect_plu.L, rect_plu.U));
}

TEST_F(LinalgPrimeTest, CharacteristicPolyCayleyHamilton) {
  std::mt19937_64 rng(42);
  for (std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{3},
                        std::size_t{4}, std::size_t{5}}) {
    auto A = linalg::random_matrix(field, {n, n}, rng);
    auto p = linalg::characteristic_poly(A);

    EXPECT_EQ(p.Degree(), static_cast<int>(n));
    EXPECT_EQ(p[n].Value(), 1);

    // Degree n-1 coefficient == -trace(A)
    auto neg_tr = -linalg::trace(A);
    EXPECT_EQ(p[n - 1], neg_tr);

    // Constant term == (-1)^n * det(A)
    auto det_a = linalg::det(A);
    auto expected_c0 = (n % 2 == 1) ? -det_a : det_a;
    EXPECT_EQ(p[0], expected_c0);

    // Cayley-Hamilton: p_A(A) == 0
    auto pA = EvalPolyAtMatrix(p, A);
    EXPECT_EQ(pA, linalg::zeros(field, {n, n}));
  }
}

TEST_F(LinalgPrimeTest, GFArrayDirectAPIUtilities) {
  auto A = MakeGFMatrix(field, 2, 3, {1, 2, 3, 4, 5, 6});
  auto B = MakeGFMatrix(field, 2, 2, {6, 5, 4, 3});
  auto C = MakeGFMatrix(field, 1, 3, {2, 0, 1});

  auto AT = linalg::transpose(A);
  EXPECT_EQ(AT.Shape(0), 3);
  EXPECT_EQ(AT.Shape(1), 2);
  EXPECT_EQ(linalg::transpose(AT), A);

  auto H = linalg::hstack(A, B);
  EXPECT_EQ(H.Shape(0), 2);
  EXPECT_EQ(H.Shape(1), 5);
  EXPECT_EQ(H.At(0, 3).Value(), 6);
  EXPECT_EQ(H.At(1, 4).Value(), 3);

  auto V = linalg::vstack(A, C);
  EXPECT_EQ(V.Shape(0), 3);
  EXPECT_EQ(V.Shape(1), 3);
  EXPECT_EQ(V.At(2, 0).Value(), 2);
  EXPECT_EQ(V.At(2, 2).Value(), 1);

  auto O = linalg::ones(field, {2, 2});
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      EXPECT_EQ(O.At(i, j).Value(), 1);
    }
  }
}

TEST_F(LinalgPrimeTest, EmptyArrayAndNullFieldErrorHandling) {
  std::shared_ptr<FieldType> null_field;
  EXPECT_THROW(linalg::eye(null_field, 2), std::invalid_argument);
  EXPECT_THROW(linalg::eye(2, null_field), std::invalid_argument);
  EXPECT_THROW(linalg::zeros(null_field, {2, 2}), std::invalid_argument);
  EXPECT_THROW(linalg::zeros({2, 2}, null_field), std::invalid_argument);
  EXPECT_THROW(linalg::ones(null_field, {2, 2}), std::invalid_argument);
  EXPECT_THROW(linalg::random_matrix(null_field, {2, 2}),
               std::invalid_argument);

  // Default-constructed GFArray has null field pointer.
  GFArray<FieldType> empty_gf;
  EXPECT_THROW(linalg::det(empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::inv(empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::rref(empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::dot(empty_gf, empty_gf), std::invalid_argument);
  EXPECT_THROW(linalg::characteristic_poly(empty_gf), std::invalid_argument);

  // Empty garray (0 elements) cannot infer field.
  xt::xarray<FieldElement> empty_garray =
      xt::xarray<FieldElement>::from_shape({0, 0});
  EXPECT_THROW(linalg::det(empty_garray), std::invalid_argument);
  EXPECT_THROW(linalg::dot(empty_garray, empty_garray), std::invalid_argument);
  EXPECT_THROW(linalg::inv(empty_garray), std::invalid_argument);

  // garray with uninitialized (null-field) elements.
  xt::xarray<FieldElement> uninit_garray({2, 2});
  EXPECT_THROW(linalg::det(uninit_garray), std::invalid_argument);

  // Mismatched fields.
  auto other_field = std::make_shared<FieldType>(11);
  auto A7 = linalg::eye(field, 2);
  auto A11 = linalg::eye(other_field, 2);
  EXPECT_THROW(linalg::dot(A7, A11), std::invalid_argument);
  EXPECT_THROW(linalg::solve(A7, linalg::ones(other_field, {2})),
               std::invalid_argument);
}
