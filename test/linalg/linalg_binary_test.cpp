#include <gtest/gtest.h>

#include <memory>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

class LinalgBinaryExtensionTest : public ::testing::Test {
 protected:
  // Using GF(2^4) with irreducible polynomial x^4 + x + 1
  // Elements are represented by uint8_t
  using FieldType = GaloisFieldBinaryExtension<uint8_t>;
  using FieldElement = GaloisFieldElement<FieldType>;

  void SetUp() override {
    // Degree 4, irreducible polynomial x^4 + x + 1 (binary string)
    field = std::make_shared<FieldType>(4, "int", "x^4 + x + 1");
  }

  std::shared_ptr<FieldType> field;
};

// Add a simple test to verify the framework is working
TEST_F(LinalgBinaryExtensionTest, BasicTest) { EXPECT_TRUE(field != nullptr); }

TEST_F(LinalgBinaryExtensionTest, VectorDotProduct) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({2});

  // v1 = [x, x+1] (values [2, 3])
  v1(0) = FieldElement(2, field);  // x
  v1(1) = FieldElement(3, field);  // x+1

  // v2 = [x^2+1, x] (values [5, 2])
  v2(0) = FieldElement(5, field);  // x^2+1
  v2(1) = FieldElement(2, field);  // x

  auto result = linalg::dot(v1, v2);

  // Calculation:
  // (x * (x^2+1)) + ((x+1) * x)
  // x * (x^2+1) = x^3 + x (binary 1010, decimal 10)
  // (x+1) * x = x^2 + x (binary 0110, decimal 6)
  // Sum: (x^3 + x) + (x^2 + x) = x^3 + x^2 (binary 1100, decimal 12)
  EXPECT_EQ(result().Value(), 12);
}

TEST_F(LinalgBinaryExtensionTest, MatrixMultiplication) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  // A = [[x, x+1], [1, x^2]] (values [[2,3],[1,4]])
  A(0, 0) = FieldElement(2, field);  // x
  A(0, 1) = FieldElement(3, field);  // x+1
  A(1, 0) = FieldElement(1, field);  // 1
  A(1, 1) = FieldElement(4, field);  // x^2

  // B = [[x^2+1, x], [x^3, x+1]] (values [[5,2],[8,3]])
  B(0, 0) = FieldElement(5, field);  // x^2+1
  B(0, 1) = FieldElement(2, field);  // x
  B(1, 0) = FieldElement(8, field);  // x^3
  B(1, 1) = FieldElement(3, field);  // x+1

  auto result = linalg::dot(A, B);

  // Expected result matrix:
  // [[3, 1],
  //  [3, 14]]
  // Calculations:
  // A*B[0,0]: (x*(x^2+1)) + ((x+1)*x^3) = (x^3+x) + (x^4+x^3)
  //          x^4 = x+1 (mod x^4+x+1)
  //          = (x^3+x) + (x+1+x^3) = 2x^3+2x+1 = 1 (in GF(2^m), 2a=0)
  // A*B[0,1]: (x*x) + ((x+1)*(x+1)) = x^2 + (x^2+1) = 1
  // A*B[1,0]: (1*(x^2+1)) + (x^2*x^3) = (x^2+1) + x^5
  //          x^5 = x*x^4 = x(x+1) = x^2+x
  //          = (x^2+1) + (x^2+x) = x+1 = 3
  // A*B[1,1]: (1*x) + (x^2*(x+1)) = x + (x^3+x^2) = x^3+x^2+x (binary 1110,
  // decimal 14)

  EXPECT_EQ(result(0, 0).Value(), 1);
  EXPECT_EQ(result(0, 1).Value(), 1);
  EXPECT_EQ(result(1, 0).Value(), 3);
  EXPECT_EQ(result(1, 1).Value(), 14);
}

TEST_F(LinalgBinaryExtensionTest, OuterProduct) {
  xt::xarray<FieldElement> v1({2});
  xt::xarray<FieldElement> v2({2});

  // v1 = [x, x+1] (values [2, 3])
  v1(0) = FieldElement(2, field);  // x
  v1(1) = FieldElement(3, field);  // x+1

  // v2 = [x^2, 1] (values [4, 1])
  v2(0) = FieldElement(4, field);  // x^2
  v2(1) = FieldElement(1, field);  // 1

  auto result = linalg::outer(v1, v2);

  // Expected: [[x*x^2, x*1], [(x+1)*x^2, (x+1)*1]]
  //         = [[x^3,     x], [x^3+x^2, x+1]]
  // Values:  [[8,       2], [12,      3]]
  EXPECT_EQ(result.shape(0), 2);
  EXPECT_EQ(result.shape(1), 2);

  EXPECT_EQ(result(0, 0).Value(), 8);   // x^3
  EXPECT_EQ(result(0, 1).Value(), 2);   // x
  EXPECT_EQ(result(1, 0).Value(), 12);  // x^3+x^2
  EXPECT_EQ(result(1, 1).Value(), 3);   // x+1
}

TEST_F(LinalgBinaryExtensionTest, Trace) {
  xt::xarray<FieldElement> A({2, 2});

  // A = [[x, x+1], [1, x^2]] (values [[2,3],[1,4]])
  A(0, 0) = FieldElement(2, field);  // x
  A(0, 1) = FieldElement(3, field);  // x+1
  A(1, 0) = FieldElement(1, field);  // 1
  A(1, 1) = FieldElement(4, field);  // x^2

  auto result = linalg::trace(A);

  // trace = x + x^2 = x^2+x (binary 0110, decimal 6)
  EXPECT_EQ(result.Value(), 6);
}

TEST_F(LinalgBinaryExtensionTest, Determinant2x2) {
  xt::xarray<FieldElement> A({2, 2});

  // A = [[x, x+1], [1, x^2]] (values [[2,3],[1,4]])
  A(0, 0) = FieldElement(2, field);  // x
  A(0, 1) = FieldElement(3, field);  // x+1
  A(1, 0) = FieldElement(1, field);  // 1
  A(1, 1) = FieldElement(4, field);  // x^2

  auto result = linalg::det(A);

  // det = (x * x^2) + ((x+1) * 1)  (since subtraction is addition in GF(2^m))
  //     = x^3 + (x+1) = x^3+x+1 (binary 1011, decimal 11)
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

  // A = [[x, 1], [0, x+1]] (values [[2,1],[0,3]])
  A(0, 0) = FieldElement(2, field);  // x
  A(0, 1) = FieldElement(1, field);  // 1
  A(1, 0) = FieldElement(0, field);  // 0
  A(1, 1) = FieldElement(3, field);  // x+1

  // Test A^0 (Identity)
  auto A0 = linalg::matrix_power(A, 0);
  EXPECT_EQ(A0(0, 0).Value(), 1);
  EXPECT_EQ(A0(0, 1).Value(), 0);
  EXPECT_EQ(A0(1, 0).Value(), 0);
  EXPECT_EQ(A0(1, 1).Value(), 1);

  // Test A^1 (A itself)
  auto A1 = linalg::matrix_power(A, 1);
  EXPECT_EQ(A1(0, 0).Value(), 2);
  EXPECT_EQ(A1(0, 1).Value(), 1);
  EXPECT_EQ(A1(1, 0).Value(), 0);
  EXPECT_EQ(A1(1, 1).Value(), 3);

  // Test A^2
  // A^2 = [[x, 1], [0, x+1]] * [[x, 1], [0, x+1]]
  //     = [[x*x + 1*0,   x*1 + 1*(x+1)],
  //        [0*x + (x+1)*0, 0*1 + (x+1)*(x+1)]]
  //     = [[x^2,     x + x+1],
  //        [0,       (x+1)^2]]
  //     = [[x^2,     1],
  //        [0,       x^2+1]]
  // Values: [[4, 1], [0, 5]]
  auto A2 = linalg::matrix_power(A, 2);
  EXPECT_EQ(A2(0, 0).Value(), 4);  // x^2
  EXPECT_EQ(A2(0, 1).Value(), 1);  // 1
  EXPECT_EQ(A2(1, 0).Value(), 0);  // 0
  EXPECT_EQ(A2(1, 1).Value(), 5);  // x^2+1
}

TEST_F(LinalgBinaryExtensionTest, VdotProduct) {
  xt::xarray<FieldElement> A({2, 2});
  xt::xarray<FieldElement> B({2, 2});

  // A = [[x, 1], [0, x+1]] (values [[2,1],[0,3]])
  // A_flat = [x, 1, 0, x+1]
  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(1, 0) = FieldElement(0, field);
  A(1, 1) = FieldElement(3, field);

  // B = [[x^2, x], [1, 0]] (values [[4,2],[1,0]])
  // B_flat = [x^2, x, 1, 0]
  B(0, 0) = FieldElement(4, field);
  B(0, 1) = FieldElement(2, field);
  B(1, 0) = FieldElement(1, field);
  B(1, 1) = FieldElement(0, field);

  auto result = linalg::vdot(A, B);

  // vdot = (x*x^2) + (1*x) + (0*1) + ((x+1)*0)
  //      = x^3 + x + 0 + 0 = x^3+x (binary 1010, decimal 10)
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

  // A = [[x, 1], [x+1, x^2]] (values [[2,1],[3,4]])
  // det(A) = x*x^2 + 1*(x+1) = x^3 + x+1 = 1011 (11) != 0, so invertible
  A(0, 0) = FieldElement(2, field);  // x
  A(0, 1) = FieldElement(1, field);  // 1
  A(1, 0) = FieldElement(3, field);  // x+1
  A(1, 1) = FieldElement(4, field);  // x^2

  // b = [x^2+x, x^3+1] (values [6, 9])
  b(0) = FieldElement(6, field);  // x^2+x
  b(1) = FieldElement(9, field);  // x^3+1

  auto x_sol = linalg::solve(A, b);

  // Expected x_sol = [x^3+x, 1] (values [10, 1])
  EXPECT_EQ(x_sol(0).Value(), 10);
  EXPECT_EQ(x_sol(1).Value(), 1);

  // Double check with A*x_sol = b
  auto verification = linalg::dot(A, x_sol);
  EXPECT_EQ(verification(0).Value(), b(0).Value());
  EXPECT_EQ(verification(1).Value(), b(1).Value());
}

TEST_F(LinalgBinaryExtensionTest, MatrixInversion) {
  xt::xarray<FieldElement> A({2, 2});
  // A = [[x, 1], [x+1, x^2]] (values [[2,1],[3,4]])
  A(0, 0) = FieldElement(2, field);  // x
  A(0, 1) = FieldElement(1, field);  // 1
  A(1, 0) = FieldElement(3, field);  // x+1
  A(1, 1) = FieldElement(4, field);  // x^2

  auto A_inv = linalg::inv(A);

  // Verify A * A_inv = I
  auto I = linalg::dot(A, A_inv);
  EXPECT_EQ(I(0, 0).Value(), 1);
  EXPECT_EQ(I(0, 1).Value(), 0);
  EXPECT_EQ(I(1, 0).Value(), 0);
  EXPECT_EQ(I(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, RowEchelonForm) {
  xt::xarray<FieldElement> A({2, 3});
  // A = [[x, 1, x^2], [x+1, x^2, x^3]]
  // Values: [[2, 1, 4], [3, 4, 8]]
  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(0, 2) = FieldElement(4, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(8, field);

  auto REF = linalg::row_echelon(A);

  // Check that it's in row echelon form (leading elements)
  EXPECT_EQ(REF(0, 0).Value(), 1);
  EXPECT_EQ(REF(1, 0).Value(), 0);
  EXPECT_EQ(REF(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, RREF) {
  xt::xarray<FieldElement> A({2, 3});
  // A = [[x, 1, x^2], [x+1, x^2, x^3]]
  // Values: [[2, 1, 4], [3, 4, 8]]
  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(0, 2) = FieldElement(4, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(8, field);

  auto R = linalg::rref(A);

  // Check that it's in reduced row echelon form
  EXPECT_EQ(R(0, 0).Value(), 1);
  EXPECT_EQ(R(0, 1).Value(), 0);
  EXPECT_EQ(R(1, 0).Value(), 0);
  EXPECT_EQ(R(1, 1).Value(), 1);
}

TEST_F(LinalgBinaryExtensionTest, MatrixRank) {
  xt::xarray<FieldElement> A({2, 3});
  // A = [[x, 1, x^2], [x+1, x^2, x^3]] -> Rank 2 from RREF
  A(0, 0) = FieldElement(2, field);
  A(0, 1) = FieldElement(1, field);
  A(0, 2) = FieldElement(4, field);
  A(1, 0) = FieldElement(3, field);
  A(1, 1) = FieldElement(4, field);
  A(1, 2) = FieldElement(8, field);
  EXPECT_EQ(linalg::matrix_rank(A), 2);

  xt::xarray<FieldElement> B({2, 2});
  // B = [[x, x^2], [x+1, x^3+x]] -> det != 0. Rank 2
  B(0, 0) = FieldElement(2, field);
  B(0, 1) = FieldElement(4, field);
  B(1, 0) = FieldElement(3, field);
  B(1, 1) = FieldElement(10, field);
  EXPECT_EQ(linalg::matrix_rank(B), 2);

  xt::xarray<FieldElement> C({2, 2});
  // C = [[x, x^2], [x+1, x^2+x]] -> det = 0. Rank 1
  C(0, 0) = FieldElement(2, field);
  C(0, 1) = FieldElement(4, field);
  C(1, 0) = FieldElement(3, field);
  C(1, 1) = FieldElement(6, field);
  EXPECT_EQ(linalg::matrix_rank(C), 1);
}
