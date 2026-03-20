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

  // 2*1 + 3*5 + 4*6 = 2 + 15 + 24 = 41 ≡ 6 (mod 7)
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

  // A = [[1, 2, 3], [4, 5, 6]]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(3, field);
  A(1, 0) = FieldElement(4, field);
  A(1, 1) = FieldElement(5, field);
  A(1, 2) = FieldElement(6, field);

  // v = [1, 2, 3]
  v(0) = FieldElement(1, field);
  v(1) = FieldElement(2, field);
  v(2) = FieldElement(3, field);

  auto result = linalg::dot(A, v);

  // [1*1 + 2*2 + 3*3, 4*1 + 5*2 + 6*3] = [14, 32] ≡ [0, 4] (mod 7)
  EXPECT_EQ(result(0).Value(), 0);
  EXPECT_EQ(result(1).Value(), 4);
}

TEST_F(LinalgPrimeTest, MatrixMultiplication) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  // A = [[1, 2], [3, 4]]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  // B = [[5, 6], [1, 2]]
  B(0, 0) = FieldElement(5, field);
  B(0, 1) = FieldElement(6, field);
  B(1, 0) = FieldElement(1, field);
  B(1, 1) = FieldElement(2, field);

  auto result = linalg::dot(A, B);

  // [[1*5 + 2*1, 1*6 + 2*2], [3*5 + 4*1, 3*6 + 4*2]]
  // = [[7, 10], [19, 26]] ≡ [[0, 3], [5, 5]] (mod 7)
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

  // [[2*4, 2*5, 2*6], [3*4, 3*5, 3*6]]
  // = [[8, 10, 12], [12, 15, 18]] ≡ [[1, 3, 5], [5, 1, 4]] (mod 7)
  EXPECT_EQ(result(0, 0).Value(), 1);
  EXPECT_EQ(result(0, 1).Value(), 3);
  EXPECT_EQ(result(0, 2).Value(), 5);
  EXPECT_EQ(result(1, 0).Value(), 5);
  EXPECT_EQ(result(1, 1).Value(), 1);
  EXPECT_EQ(result(1, 2).Value(), 4);
}

TEST_F(LinalgPrimeTest, Trace) {
  xt::xarray<FieldElement> A({3, 3});

  // A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(0, 2) = FieldElement(3, field);
  A(1, 0) = FieldElement(4, field);
  A(1, 1) = FieldElement(5, field);
  A(1, 2) = FieldElement(6, field);
  A(2, 0) = FieldElement(0, field);  // 7 ≡ 0 (mod 7)
  A(2, 1) = FieldElement(1, field);  // 8 ≡ 1 (mod 7)
  A(2, 2) = FieldElement(2, field);  // 9 ≡ 2 (mod 7)

  auto result = linalg::trace(A);

  // trace = 1 + 5 + 2 = 8 ≡ 1 (mod 7)
  EXPECT_EQ(result.Value(), 1);
}

TEST_F(LinalgPrimeTest, TraceNonSquareMatrix) {
  xt::xarray<FieldElement> A({2, 3});

  EXPECT_THROW(linalg::trace(A), std::invalid_argument);
}

TEST_F(LinalgPrimeTest, Determinant2x2) {
  xt::xarray<FieldElement> A({2, 2});

  // A = [[1, 2], [3, 4]]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  auto result = linalg::det(A);

  // det = 1*4 - 2*3 = 4 - 6 = -2 ≡ 5 (mod 7)
  EXPECT_EQ(result.Value(), 5);
}

TEST_F(LinalgPrimeTest, Determinant3x3) {
  xt::xarray<FieldElement> A({3, 3});

  // A = [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
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

  // Using cofactor expansion along first row:
  // det = 1*(1*0 - 4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5)
  //     = 1*(-24) - 2*(-20) + 3*(-5)
  //     = -24 + 40 - 15 = 1 ≡ 1 (mod 7)
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

  // A = [[1, 2], [3, 4]]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  auto A2 = linalg::matrix_power(A, 2);

  // A^2 = [[1*1 + 2*3, 1*2 + 2*4], [3*1 + 4*3, 3*2 + 4*4]]
  //     = [[7, 10], [15, 22]] ≡ [[0, 3], [1, 1]] (mod 7)
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

  // A^0 should be identity matrix
  EXPECT_EQ(A0(0, 0).Value(), 1);
  EXPECT_EQ(A0(0, 1).Value(), 0);
  EXPECT_EQ(A0(1, 0).Value(), 0);
  EXPECT_EQ(A0(1, 1).Value(), 1);
}

TEST_F(LinalgPrimeTest, LinearSystemSolve) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> b({2});

  // A = [[1, 2], [3, 4]], b = [5, 6]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);

  b(0) = FieldElement(5, field);
  b(1) = FieldElement(6, field);

  auto x = linalg::solve(A, b);

  // Verify Ax = b
  auto verification = linalg::dot(A, x);
  EXPECT_EQ(verification(0).Value(), b(0).Value());
  EXPECT_EQ(verification(1).Value(), b(1).Value());
}

TEST_F(LinalgPrimeTest, MatrixInversion) {
  xt::xarray<FieldElement> A({2, 2});

  // A = [[1, 2], [3, 5]] (det = 1*5 - 2*3 = -1 ≡ 6 (mod 7), invertible)
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(5, field);

  auto A_inv = linalg::inv(A);

  // Verify A * A^(-1) = I
  auto verification = linalg::dot(A, A_inv);

  EXPECT_EQ(verification(0, 0).Value(), 1);
  EXPECT_EQ(verification(0, 1).Value(), 0);
  EXPECT_EQ(verification(1, 0).Value(), 0);
  EXPECT_EQ(verification(1, 1).Value(), 1);
}

TEST_F(LinalgPrimeTest, SingularMatrixInversion) {
  xt::xarray<FieldElement> A({2, 2});

  // Singular matrix: A = [[1, 2], [2, 4]] (det = 0)
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
  B(1, 0) = FieldElement(0, field);  // 7 ≡ 0 (mod 7)
  B(1, 1) = FieldElement(1, field);  // 8 ≡ 1 (mod 7)

  auto result = linalg::vdot(A, B);

  // vdot flattens both arrays and computes dot product
  // A_flat = [1, 2, 3, 4], B_flat = [5, 6, 0, 1]
  // result = 1*5 + 2*6 + 3*0 + 4*1 = 5 + 12 + 0 + 4 = 21 ≡ 0 (mod 7)
  EXPECT_EQ(result.Value(), 0);
}

TEST_F(LinalgPrimeTest, KroneckerProduct) {
  xt::xarray<FieldElement> A({1, 2});
  xt::xarray<FieldElement> B({2, 2});

  // A = [[1, 2]]
  A(0, 0) = FieldElement(1, field);
  A(0, 1) = FieldElement(2, field);

  // B = [[3, 4], [5, 6]]
  B(0, 0) = FieldElement(3, field);
  B(0, 1) = FieldElement(4, field);
  B(1, 0) = FieldElement(5, field);
  B(1, 1) = FieldElement(6, field);

  auto result = linalg::kron(A, B);

  // Expected: [[1*B, 2*B]]
  // [[1*[3,4], 2*[3,4]],
  //  [1*[5,6], 2*[5,6]]]
  // = [[[3,4], [6,1]],
  //    [[5,6], [3,5]]]
  // Reshaped to 2x4:
  // [[3, 4, 6, 1],
  //  [5, 6, 3, 5]]

  EXPECT_EQ(result.shape(0), 2);
  EXPECT_EQ(result.shape(1), 4);

  EXPECT_EQ(result(0, 0).Value(), 3);  // 1*3
  EXPECT_EQ(result(0, 1).Value(), 4);  // 1*4
  EXPECT_EQ(result(0, 2).Value(), 6);  // 2*3
  EXPECT_EQ(result(0, 3).Value(), 1);  // 2*4 = 8 % 7 = 1

  EXPECT_EQ(result(1, 0).Value(), 5);  // 1*5
  EXPECT_EQ(result(1, 1).Value(), 6);  // 1*6
  EXPECT_EQ(result(1, 2).Value(), 3);  // 2*5 = 10 % 7 = 3
  EXPECT_EQ(result(1, 3).Value(), 5);  // 2*6 = 12 % 7 = 5
}

TEST_F(LinalgPrimeTest, RowEchelonForm) {
  xt::xarray<FieldElement> A({3, 4});

  // A = [[1, 2, 1, 0],
  //      [2, 4, 3, 1],
  //      [3, 6, 0, 5]]
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

  // Check that it's in row echelon form (leading elements)
  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(2, 0).Value(), 0);
}

TEST_F(LinalgPrimeTest, MatrixRank) {
  xt::xarray<FieldElement> A({3, 4});
  // A = [[1, 2, 1, 0],
  //      [2, 4, 3, 1],
  //      [3, 6, 0, 5]]
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
  // B = [[1, 2, 3],
  //      [2, 4, 6]]
  // Rank is 1
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

  // v1 = [1, 2, 3]
  v1(0) = FieldElement(1, field);
  v1(1) = FieldElement(2, field);
  v1(2) = FieldElement(3, field);

  // v2 = [4, 5, 6]
  v2(0) = FieldElement(4, field);
  v2(1) = FieldElement(5, field);
  v2(2) = FieldElement(6, field);

  auto result = linalg::cross(v1, v2);

  // c0 = a1*b2 - a2*b1 = 2*6 - 3*5 = 12 - 15 = -3 ≡ 4 (mod 7)
  // c1 = a2*b0 - a0*b2 = 3*4 - 1*6 = 12 - 6 = 6 ≡ 6 (mod 7)
  // c2 = a0*b1 - a1*b0 = 1*5 - 2*4 = 5 - 8 = -3 ≡ 4 (mod 7)
  // result = [4, 6, 4]
  EXPECT_EQ(result.shape(0), 3);
  EXPECT_EQ(result(0).Value(), 4);
  EXPECT_EQ(result(1).Value(), 6);
  EXPECT_EQ(result(2).Value(), 4);
}

TEST_F(LinalgPrimeTest, CrossProduct2D) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({2});

  // v1 = [1, 2] (effectively [1, 2, 0])
  v1(0) = FieldElement(1, field);
  v1(1) = FieldElement(2, field);

  // v2 = [3, 4] (effectively [3, 4, 0])
  v2(0) = FieldElement(3, field);
  v2(1) = FieldElement(4, field);

  auto result = linalg::cross(v1, v2);
  // a0=1, a1=2, a2=0
  // b0=3, b1=4, b2=0
  // c0 = a1*b2 - a2*b1 = 2*0 - 0*4 = 0
  // c1 = a2*b0 - a0*b2 = 0*3 - 1*0 = 0
  // c2 = a0*b1 - a1*b0 = 1*4 - 2*3 = 4 - 6 = -2 ≡ 5 (mod 7)
  // result = [0, 0, 5]
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

  // A = [[1, 2, 1, 0],
  //      [2, 4, 3, 1],
  //      [3, 6, 0, 5]]
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

  // Check that it's in reduced row echelon form
  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(2, 0).Value(), 0);

  // First pivot column should have 1 in first row, 0s below
  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(2, 0).Value(), 0);
}
