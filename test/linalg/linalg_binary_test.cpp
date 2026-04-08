#include <gtest/gtest.h>

#include <memory>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

class LinalgBinaryExtensionTest : public ::testing::Test {
 protected:
  using FieldType = GaloisFieldBinaryExtension<uint8_t>;
  using FieldElement = GaloisFieldElement<FieldType>;

  void SetUp() override {
    field = std::make_shared<FieldType>(4, "int", "x^4 + x + 1");
  }

  std::shared_ptr<FieldType> field;
};

TEST_F(LinalgBinaryExtensionTest, BasicTest) { EXPECT_TRUE(field != nullptr); }

TEST_F(LinalgBinaryExtensionTest, VectorDotProduct) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({2});

  v1(0) = FieldElement(2, field);
  v1(1) = FieldElement(3, field);

  v2(0) = FieldElement(5, field);
  v2(1) = FieldElement(2, field);

  auto result = linalg::dot(v1, v2);

  EXPECT_EQ(result().Value(), 12);
}

TEST_F(LinalgBinaryExtensionTest, MatrixMultiplication) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(3, field);
  A(1, 0) = FieldElement(1, field);
  A(1, 1) = FieldElement(4, field);

  B(0, 0) = FieldElement(5, field);
  B(0, 1) = FieldElement(2, field);
  B(1, 0) = FieldElement(8, field);
  B(1, 1) = FieldElement(3, field);

  auto result = linalg::dot(A, B);

  EXPECT_EQ(result(0, 0).Value(), 1);
  EXPECT_EQ(result(0, 1).Value(), 1);
  EXPECT_EQ(result(1, 0).Value(), 3);
  EXPECT_EQ(result(1, 1).Value(), 14);
}

TEST_F(LinalgBinaryExtensionTest, OuterProduct) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({2});

  v1(0) = FieldElement(2, field);
  v1(1) = FieldElement(3, field);

  v2(0) = FieldElement(4, field);
  v2(1) = FieldElement(1, field);

  auto result = linalg::outer(v1, v2);

  EXPECT_EQ(result.shape(0), 2);
  EXPECT_EQ(result.shape(1), 2);

  EXPECT_EQ(result(0, 0).Value(), 8);
  EXPECT_EQ(result(0, 1).Value(), 2);
  EXPECT_EQ(result(1, 0).Value(), 12);
  EXPECT_EQ(result(1, 1).Value(), 3);
}

TEST_F(LinalgBinaryExtensionTest, Trace) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(3, field);
  A(1, 0) = FieldElement(1, field);
  A(1, 1) = FieldElement(4, field);

  auto result = linalg::trace(A);

  EXPECT_EQ(result.Value(), 6);
}

TEST_F(LinalgBinaryExtensionTest, Determinant2x2) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(3, field);
  A(1, 0) = FieldElement(1, field);
  A(1, 1) = FieldElement(4, field);

  auto result = linalg::det(A);

  EXPECT_EQ(result.Value(), 11);
}

TEST_F(LinalgBinaryExtensionTest, IdentityMatrix) {
  auto I = linalg::eye(2, field);

  EXPECT_EQ(I.shape(0), 2);
  EXPECT_EQ(I.shape(1), 2);

  EXPECT_EQ(I(0, 0).Value(), 1);
  EXPECT_EQ(I(0, 1).Value(), 0);
  EXPECT_EQ(I(1, 0).Value(), 0);
  EXPECT_EQ(I(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, MatrixPower) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(1, 0) = FieldElement(0, field);
  A(1, 1) = FieldElement(3, field);

  auto A0 = linalg::matrix_power(A, 0);
  EXPECT_EQ(A0(0, 0).Value(), 1);
  EXPECT_EQ(A0(0, 1).Value(), 0);
  EXPECT_EQ(A0(1, 0).Value(), 0);
  EXPECT_EQ(A0(1, 1).Value(), 1);

  auto A1 = linalg::matrix_power(A, 1);
  EXPECT_EQ(A1(0, 0).Value(), 2);
  EXPECT_EQ(A1(0, 1).Value(), 1);
  EXPECT_EQ(A1(1, 0).Value(), 0);
  EXPECT_EQ(A1(1, 1).Value(), 3);

  auto A2 = linalg::matrix_power(A, 2);
  EXPECT_EQ(A2(0, 0).Value(), 4);
  EXPECT_EQ(A2(0, 1).Value(), 1);
  EXPECT_EQ(A2(1, 0).Value(), 0);
  EXPECT_EQ(A2(1, 1).Value(), 5);
}

TEST_F(LinalgBinaryExtensionTest, VdotProduct) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(1, 0) = FieldElement(0, field);
  A(1, 1) = FieldElement(3, field);

  B(0, 0) = FieldElement(4, field);
  B(0, 1) = FieldElement(2, field);
  B(1, 0) = FieldElement(1, field);
  B(1, 1) = FieldElement(0, field);

  auto result = linalg::vdot(A, B);

  EXPECT_EQ(result.Value(), 10);
}

TEST_F(LinalgBinaryExtensionTest, ZerosMatrix) {
  auto Z = linalg::zeros({2, 2}, field);

  EXPECT_EQ(Z.shape(0), 2);
  EXPECT_EQ(Z.shape(1), 2);

  EXPECT_EQ(Z(0, 0).Value(), 0);
  EXPECT_EQ(Z(0, 1).Value(), 0);
  EXPECT_EQ(Z(1, 0).Value(), 0);
  EXPECT_EQ(Z(1, 1).Value(), 0);
}

TEST_F(LinalgBinaryExtensionTest, LinearSystemSolve) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> b({2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  b(0) = FieldElement(6, field);
  b(1) = FieldElement(9, field);

  auto x_sol = linalg::solve(A, b);

  EXPECT_EQ(x_sol(0).Value(), 10);
  EXPECT_EQ(x_sol(1).Value(), 1);

  auto verification = linalg::dot(A, x_sol);
  EXPECT_EQ(verification(0).Value(), b(0).Value());
  EXPECT_EQ(verification(1).Value(), b(1).Value());
}

TEST_F(LinalgBinaryExtensionTest, MatrixInversion) {
  xt::xarray<FieldElement> A({2, 2});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  auto A_inv = linalg::inv(A);

  auto I = linalg::dot(A, A_inv);
  EXPECT_EQ(I(0, 0).Value(), 1);
  EXPECT_EQ(I(0, 1).Value(), 0);
  EXPECT_EQ(I(1, 0).Value(), 0);
  EXPECT_EQ(I(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, RowEchelonForm) {
  xt::xarray<FieldElement> A({2, 3});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(0, 2) = FieldElement(4, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(8, field);

  auto REF = linalg::row_echelon(A);

  EXPECT_EQ(REF(0, 0).Value(), 1);
  EXPECT_EQ(REF(1, 0).Value(), 0);
  EXPECT_EQ(REF(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, RREF) {
  xt::xarray<FieldElement> A({2, 3});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(0, 2) = FieldElement(4, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(8, field);

  auto R = linalg::rref(A);

  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(0, 1).Value(), 0);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, MatrixRank) {
  xt::xarray<FieldElement> A({2, 3});

  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(0, 2) = FieldElement(4, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(8, field);
  EXPECT_EQ(linalg::matrix_rank(A), 2);

  xt::xarray<FieldElement> B({2, 2});

  B(0, 0) = FieldElement(2, field);
  B(0, 1) = FieldElement(4, field);
  B(1, 0) = FieldElement(3, field);
  B(1, 1) = FieldElement(10, field);
  EXPECT_EQ(linalg::matrix_rank(B), 2);

  xt::xarray<FieldElement> C({2, 2});

  C(0, 0) = FieldElement(2, field);
  C(0, 1) = FieldElement(4, field);
  C(1, 0) = FieldElement(3, field);
  C(1, 1) = FieldElement(6, field);
  EXPECT_EQ(linalg::matrix_rank(C), 1);
}
