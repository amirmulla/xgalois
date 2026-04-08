#include <gtest/gtest.h>

#include <memory>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

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
  xt::xarray<FieldElement> v1({3});
  xt::xarray<FieldElement> v2({4});

  EXPECT_THROW(linalg::dot(v1, v2), std::invalid_argument);
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
  xt::xarray<FieldElement> A({2, 3});

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
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(2, 0).Value(), 0);
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

  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(2, 0).Value(), 0);

  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(2, 0).Value(), 0);
}
