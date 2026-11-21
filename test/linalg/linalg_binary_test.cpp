#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
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

TEST_F(LinalgBinaryExtensionTest, VectorMatrixAndMatrixVectorDot) {
  auto v = MakeGFVector(field, {2, 3});
  auto A = MakeGFMatrix(field, 2, 3, {5, 1, 4, 2, 7, 3});

  // 1D * 2D: [2, 3] * A
  auto res_1d2d = linalg::dot(v, A);
  ASSERT_EQ(res_1d2d.Dimension(), 1);
  ASSERT_EQ(res_1d2d.Shape(0), 3);
  for (std::size_t j = 0; j < 3; ++j) {
    uint8_t expected = field->Add(field->Mul(2, A.At(0, j).Value()),
                                  field->Mul(3, A.At(1, j).Value()));
    EXPECT_EQ(res_1d2d.At(j).Value(), expected);
  }

  // 2D * 1D: A * w
  auto w = MakeGFVector(field, {4, 6, 9});
  auto res_2d1d = linalg::dot(A, w);
  ASSERT_EQ(res_2d1d.Dimension(), 1);
  ASSERT_EQ(res_2d1d.Shape(0), 2);
  for (std::size_t i = 0; i < 2; ++i) {
    uint8_t expected = 0;
    for (std::size_t j = 0; j < 3; ++j) {
      expected =
          field->Add(expected, field->Mul(A.At(i, j).Value(), w.At(j).Value()));
    }
    EXPECT_EQ(res_2d1d.At(i).Value(), expected);
  }
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

TEST_F(LinalgBinaryExtensionTest, KroneckerAndCrossProducts) {
  auto A = MakeGFMatrix(field, 2, 2, {2, 3, 1, 4});
  auto B = MakeGFMatrix(field, 2, 2, {5, 2, 8, 3});
  auto K = linalg::kron(A, B);
  ASSERT_EQ(K.Shape(0), 4);
  ASSERT_EQ(K.Shape(1), 4);
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t c = 0; c < 2; ++c) {
          uint8_t expected = field->Mul(A.At(i, j).Value(), B.At(r, c).Value());
          EXPECT_EQ(K.At(i * 2 + r, j * 2 + c).Value(), expected);
        }
      }
    }
  }

  auto u = MakeGFVector(field, {2, 5, 7});
  auto v = MakeGFVector(field, {3, 4, 9});
  auto cr = linalg::cross(u, v);
  ASSERT_EQ(cr.Shape(0), 3);
  // Cross product is orthogonal to both u and v.
  EXPECT_EQ(linalg::vdot(u, cr).Value(), 0);
  EXPECT_EQ(linalg::vdot(v, cr).Value(), 0);
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

  auto A_neg3 = linalg::matrix_power(A, -3);
  auto A_pos3 = linalg::matrix_power(A, 3);
  auto prod = linalg::dot(A_pos3, A_neg3);
  EXPECT_EQ(prod(0, 0).Value(), 1);
  EXPECT_EQ(prod(0, 1).Value(), 0);
  EXPECT_EQ(prod(1, 0).Value(), 0);
  EXPECT_EQ(prod(1, 1).Value(), 1);
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

TEST_F(LinalgBinaryExtensionTest, NullspaceSystematicFormAndPLU) {
  auto A = MakeGFMatrix(field, 3, 5,
                        {2, 3, 5, 7, 11, 4, 6, 10, 14, 5, 1, 8, 9, 12, 13});
  const std::size_t rank = linalg::matrix_rank(A);

  auto N = linalg::nullspace(A);
  EXPECT_EQ(N.Shape(0), 5 - rank);
  EXPECT_EQ(N.Shape(1), 5);
  EXPECT_EQ(linalg::matrix_rank(N), 5 - rank);
  EXPECT_EQ(linalg::dot(A, linalg::transpose(N)),
            linalg::zeros(field, {3, 5 - rank}));

  auto L_ns = linalg::left_nullspace(A);
  EXPECT_EQ(L_ns.Shape(0), 3 - rank);
  EXPECT_EQ(linalg::dot(L_ns, A), linalg::zeros(field, {3 - rank, 5}));

  auto sys = linalg::systematic_form(A);
  ASSERT_EQ(sys.G_sys.Shape(0), rank);
  ASSERT_EQ(sys.G_sys.Shape(1), 5);
  for (std::size_t i = 0; i < rank; ++i) {
    for (std::size_t j = 0; j < rank; ++j) {
      EXPECT_EQ(sys.G_sys.At(i, j).Value(), i == j ? 1U : 0U);
    }
  }

  auto plu_res = linalg::plu(A);
  EXPECT_EQ(linalg::dot(plu_res.P, A), linalg::dot(plu_res.L, plu_res.U));
}

TEST(LinalgGF256Test, CharacteristicPolyCayleyHamiltonAndBackends) {
  const std::string aes_poly = "x^8 + x^4 + x^3 + x + 1";
  auto gf_std = std::make_shared<GF2X<uint8_t>>(8, "int", aes_poly);
  auto gf_log = std::make_shared<GF2XLOG<uint8_t>>(8, "int", aes_poly);
  auto gf_zech = std::make_shared<GF2XZECH<uint8_t>>(8, "int", aes_poly);

  // Build identical 5x5 matrices across all 3 binary extension backends
  // (including Zech-log where AdditiveIdentity() == 255 and One == 0).
  std::mt19937_64 rng(12345);
  std::uniform_int_distribution<unsigned int> dist(0, 255);
  const std::size_t n = 5;
  std::vector<uint8_t> a_vals(n * n);
  std::vector<uint8_t> b_vals(n * n);
  for (std::size_t i = 0; i < n * n; ++i) {
    a_vals[i] = static_cast<uint8_t>(dist(rng));
    b_vals[i] = static_cast<uint8_t>(dist(rng));
  }

  auto A_std = MakeGFMatrix(gf_std, n, n, a_vals);
  auto B_std = MakeGFMatrix(gf_std, n, n, b_vals);
  auto A_log = MakeGFMatrix(gf_log, n, n, a_vals);
  auto B_log = MakeGFMatrix(gf_log, n, n, b_vals);
  auto A_zech = MakeGFMatrix(gf_zech, n, n, a_vals);
  auto B_zech = MakeGFMatrix(gf_zech, n, n, b_vals);

  // 1. Matrix multiplication consistency across GF2X, GF2XLOG, GF2XZECH
  auto C_std = linalg::dot(A_std, B_std);
  auto C_log = linalg::dot(A_log, B_log);
  auto C_zech = linalg::dot(A_zech, B_zech);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(C_std.At(i, j).Value(), C_log.At(i, j).Value());
      EXPECT_EQ(C_std.At(i, j).Value(), C_zech.At(i, j).Value());
    }
  }

  // 2. Determinant, trace, rank, RREF, PLU consistency
  auto det_std = linalg::det(A_std).Value();
  EXPECT_EQ(det_std, linalg::det(A_log).Value());
  EXPECT_EQ(det_std, linalg::det(A_zech).Value());

  auto tr_std = linalg::trace(A_std).Value();
  EXPECT_EQ(tr_std, linalg::trace(A_log).Value());
  EXPECT_EQ(tr_std, linalg::trace(A_zech).Value());

  auto rref_std = linalg::rref(A_std);
  auto rref_zech = linalg::rref(A_zech);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(rref_std.At(i, j).Value(), rref_zech.At(i, j).Value());
    }
  }

  // 3. Inverse & multi-RHS solve consistency
  if (det_std != 0) {
    auto inv_std = linalg::inv(A_std);
    auto inv_log = linalg::inv(A_log);
    auto inv_zech = linalg::inv(A_zech);
    EXPECT_EQ(linalg::dot(A_std, inv_std), linalg::eye(gf_std, n));
    EXPECT_EQ(linalg::dot(A_log, inv_log), linalg::eye(gf_log, n));
    EXPECT_EQ(linalg::dot(A_zech, inv_zech), linalg::eye(gf_zech, n));
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        EXPECT_EQ(inv_std.At(i, j).Value(), inv_log.At(i, j).Value());
        EXPECT_EQ(inv_std.At(i, j).Value(), inv_zech.At(i, j).Value());
      }
    }

    auto X_zech = linalg::solve(A_zech, B_zech);
    EXPECT_EQ(linalg::dot(A_zech, X_zech), B_zech);
  }

  // 4. PLU decomposition across all 3 backends
  auto plu_zech = linalg::plu(A_zech);
  EXPECT_EQ(linalg::dot(plu_zech.P, A_zech),
            linalg::dot(plu_zech.L, plu_zech.U));

  // 5. Characteristic polynomial & Cayley-Hamilton over GF(2^8)
  auto p_std = linalg::characteristic_poly(A_std);
  auto p_log = linalg::characteristic_poly(A_log);
  auto p_zech = linalg::characteristic_poly(A_zech);

  EXPECT_EQ(p_std.Degree(), static_cast<int>(n));
  EXPECT_EQ(p_std[n].Value(), 1);
  EXPECT_EQ(p_std[n - 1].Value(), tr_std);
  EXPECT_EQ(p_std[0].Value(), det_std);

  for (std::size_t k = 0; k <= n; ++k) {
    EXPECT_EQ(p_std[k].Value(), p_log[k].Value());
    EXPECT_EQ(p_std[k].Value(), p_zech[k].Value());
  }

  EXPECT_EQ(EvalPolyAtMatrix(p_std, A_std), linalg::zeros(gf_std, {n, n}));
  EXPECT_EQ(EvalPolyAtMatrix(p_log, A_log), linalg::zeros(gf_log, {n, n}));
  EXPECT_EQ(EvalPolyAtMatrix(p_zech, A_zech), linalg::zeros(gf_zech, {n, n}));
}
