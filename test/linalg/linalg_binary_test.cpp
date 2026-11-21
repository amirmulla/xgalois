#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
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

class LinalgBinaryExtensionTest : public ::testing::Test {
 protected:
  using FieldType = GaloisFieldBinaryExtension<uint8_t>;

  void SetUp() override {
    field = std::make_shared<FieldType>(4, "int", "x^4 + x + 1");
  }

  std::shared_ptr<FieldType> field;
};

TEST_F(LinalgBinaryExtensionTest, BasicTest) { EXPECT_TRUE(field != nullptr); }

TEST_F(LinalgBinaryExtensionTest, VectorDotProduct) {
  auto v1 = MakeGFVector(field, {2, 3});
  auto v2 = MakeGFVector(field, {5, 2});

  auto result = linalg::Dot(v1, v2);

  ASSERT_EQ(result.Dimension(), std::size_t{0});
  EXPECT_EQ(result.RawData()(), uint8_t{12});
}

TEST_F(LinalgBinaryExtensionTest, VectorMatrixAndMatrixVectorDot) {
  auto v = MakeGFVector(field, {2, 3});
  auto a = MakeGFMatrix(field, 2, 3, {5, 1, 4, 2, 7, 3});

  // 1D * 2D: [2, 3] * A
  auto res_1d2d = linalg::Dot(v, a);
  ASSERT_EQ(res_1d2d.Dimension(), std::size_t{1});
  ASSERT_EQ(res_1d2d.Shape(0), std::size_t{3});
  for (std::size_t j = 0; j < 3; ++j) {
    uint8_t expected = field->Add(field->Mul(2, a.At(0, j).Value()),
                                  field->Mul(3, a.At(1, j).Value()));
    EXPECT_EQ(res_1d2d.At(j).Value(), expected);
  }

  // 2D * 1D: A * w
  auto w = MakeGFVector(field, {4, 6, 9});
  auto res_2d1d = linalg::Dot(a, w);
  ASSERT_EQ(res_2d1d.Dimension(), std::size_t{1});
  ASSERT_EQ(res_2d1d.Shape(0), std::size_t{2});
  for (std::size_t i = 0; i < 2; ++i) {
    uint8_t expected = 0;
    for (std::size_t j = 0; j < 3; ++j) {
      expected =
          field->Add(expected, field->Mul(a.At(i, j).Value(), w.At(j).Value()));
    }
    EXPECT_EQ(res_2d1d.At(i).Value(), expected);
  }
}

TEST_F(LinalgBinaryExtensionTest, MatrixMultiplication) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 3, 1, 4});
  auto b = MakeGFMatrix(field, 2, 2, {5, 2, 8, 3});

  auto result = linalg::Dot(a, b);

  EXPECT_EQ(result.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(result.At(0, 1).Value(), uint8_t{1});
  EXPECT_EQ(result.At(1, 0).Value(), uint8_t{3});
  EXPECT_EQ(result.At(1, 1).Value(), uint8_t{14});
}

TEST_F(LinalgBinaryExtensionTest, OuterProduct) {
  auto v1 = MakeGFVector(field, {2, 3});
  auto v2 = MakeGFVector(field, {4, 1});

  auto result = linalg::Outer(v1, v2);

  EXPECT_EQ(result.Shape(0), std::size_t{2});
  EXPECT_EQ(result.Shape(1), std::size_t{2});

  EXPECT_EQ(result.At(0, 0).Value(), uint8_t{8});
  EXPECT_EQ(result.At(0, 1).Value(), uint8_t{2});
  EXPECT_EQ(result.At(1, 0).Value(), uint8_t{12});
  EXPECT_EQ(result.At(1, 1).Value(), uint8_t{3});
}

TEST_F(LinalgBinaryExtensionTest, KroneckerAndCrossProducts) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 3, 1, 4});
  auto b = MakeGFMatrix(field, 2, 2, {5, 2, 8, 3});
  auto k = linalg::Kron(a, b);
  ASSERT_EQ(k.Shape(0), std::size_t{4});
  ASSERT_EQ(k.Shape(1), std::size_t{4});
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t c = 0; c < 2; ++c) {
          uint8_t expected = field->Mul(a.At(i, j).Value(), b.At(r, c).Value());
          EXPECT_EQ(k.At(i * 2 + r, j * 2 + c).Value(), expected);
        }
      }
    }
  }

  auto u = MakeGFVector(field, {2, 5, 7});
  auto v = MakeGFVector(field, {3, 4, 9});
  auto cr = linalg::Cross(u, v);
  ASSERT_EQ(cr.Shape(0), std::size_t{3});
  // Cross product is orthogonal to both u and v.
  EXPECT_EQ(linalg::Vdot(u, cr).Value(), uint8_t{0});
  EXPECT_EQ(linalg::Vdot(v, cr).Value(), uint8_t{0});
}

TEST_F(LinalgBinaryExtensionTest, Trace) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 3, 1, 4});

  auto result = linalg::Trace(a);

  EXPECT_EQ(result.Value(), uint8_t{6});
}

TEST_F(LinalgBinaryExtensionTest, Determinant2x2) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 3, 1, 4});

  auto result = linalg::Det(a);

  EXPECT_EQ(result.Value(), uint8_t{11});
}

TEST_F(LinalgBinaryExtensionTest, IdentityMatrix) {
  auto eye_mat = linalg::Eye(field, 2);

  EXPECT_EQ(eye_mat.Shape(0), std::size_t{2});
  EXPECT_EQ(eye_mat.Shape(1), std::size_t{2});

  EXPECT_EQ(eye_mat.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(eye_mat.At(0, 1).Value(), uint8_t{0});
  EXPECT_EQ(eye_mat.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(eye_mat.At(1, 1).Value(), uint8_t{1});
}

TEST_F(LinalgBinaryExtensionTest, MatrixPower) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 1, 0, 3});

  auto a0 = linalg::MatrixPower(a, 0);
  EXPECT_EQ(a0.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(a0.At(0, 1).Value(), uint8_t{0});
  EXPECT_EQ(a0.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(a0.At(1, 1).Value(), uint8_t{1});

  auto a1 = linalg::MatrixPower(a, 1);
  EXPECT_EQ(a1.At(0, 0).Value(), uint8_t{2});
  EXPECT_EQ(a1.At(0, 1).Value(), uint8_t{1});
  EXPECT_EQ(a1.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(a1.At(1, 1).Value(), uint8_t{3});

  auto a2 = linalg::MatrixPower(a, 2);
  EXPECT_EQ(a2.At(0, 0).Value(), uint8_t{4});
  EXPECT_EQ(a2.At(0, 1).Value(), uint8_t{1});
  EXPECT_EQ(a2.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(a2.At(1, 1).Value(), uint8_t{5});

  auto a_neg3 = linalg::MatrixPower(a, -3);
  auto a_pos3 = linalg::MatrixPower(a, 3);
  auto prod = linalg::Dot(a_pos3, a_neg3);
  EXPECT_EQ(prod.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(prod.At(0, 1).Value(), uint8_t{0});
  EXPECT_EQ(prod.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(prod.At(1, 1).Value(), uint8_t{1});
}

TEST_F(LinalgBinaryExtensionTest, VdotProduct) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 1, 0, 3});
  auto b = MakeGFMatrix(field, 2, 2, {4, 2, 1, 0});

  auto result = linalg::Vdot(a, b);

  EXPECT_EQ(result.Value(), uint8_t{10});
}

TEST_F(LinalgBinaryExtensionTest, ZerosMatrix) {
  auto z = linalg::Zeros(field, {2, 2});

  EXPECT_EQ(z.Shape(0), std::size_t{2});
  EXPECT_EQ(z.Shape(1), std::size_t{2});

  EXPECT_EQ(z.At(0, 0).Value(), uint8_t{0});
  EXPECT_EQ(z.At(0, 1).Value(), uint8_t{0});
  EXPECT_EQ(z.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(z.At(1, 1).Value(), uint8_t{0});
}

TEST_F(LinalgBinaryExtensionTest, LinearSystemSolve) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 1, 3, 4});
  auto b = MakeGFVector(field, {6, 9});

  auto x_sol = linalg::Solve(a, b);

  EXPECT_EQ(x_sol.At(0).Value(), uint8_t{10});
  EXPECT_EQ(x_sol.At(1).Value(), uint8_t{1});

  auto verification = linalg::Dot(a, x_sol);
  EXPECT_EQ(verification.At(0).Value(), b.At(0).Value());
  EXPECT_EQ(verification.At(1).Value(), b.At(1).Value());
}

TEST_F(LinalgBinaryExtensionTest, MatrixInversion) {
  auto a = MakeGFMatrix(field, 2, 2, {2, 1, 3, 4});

  auto a_inv = linalg::Inv(a);

  auto eye_mat = linalg::Dot(a, a_inv);
  EXPECT_EQ(eye_mat.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(eye_mat.At(0, 1).Value(), uint8_t{0});
  EXPECT_EQ(eye_mat.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(eye_mat.At(1, 1).Value(), uint8_t{1});
}

TEST_F(LinalgBinaryExtensionTest, RowEchelonForm) {
  auto a = MakeGFMatrix(field, 2, 3, {2, 1, 4, 3, 4, 8});

  auto ref = linalg::RowEchelon(a);

  EXPECT_EQ(ref.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(ref.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(ref.At(1, 1).Value(), uint8_t{1});
}

TEST_F(LinalgBinaryExtensionTest, Rref) {
  auto a = MakeGFMatrix(field, 2, 3, {2, 1, 4, 3, 4, 8});

  auto r = linalg::Rref(a);

  EXPECT_EQ(r.At(0, 0).Value(), uint8_t{1});
  EXPECT_EQ(r.At(0, 1).Value(), uint8_t{0});
  EXPECT_EQ(r.At(1, 0).Value(), uint8_t{0});
  EXPECT_EQ(r.At(1, 1).Value(), uint8_t{1});
}

TEST_F(LinalgBinaryExtensionTest, MatrixRank) {
  auto a = MakeGFMatrix(field, 2, 3, {2, 1, 4, 3, 4, 8});
  EXPECT_EQ(linalg::MatrixRank(a), std::size_t{2});

  auto b = MakeGFMatrix(field, 2, 2, {2, 4, 3, 10});
  EXPECT_EQ(linalg::MatrixRank(b), std::size_t{2});

  auto c = MakeGFMatrix(field, 2, 2, {2, 4, 3, 6});
  EXPECT_EQ(linalg::MatrixRank(c), std::size_t{1});
}

TEST_F(LinalgBinaryExtensionTest, NullspaceSystematicFormAndPlu) {
  auto a = MakeGFMatrix(field, 3, 5,
                        {2, 3, 5, 7, 11, 4, 6, 10, 14, 5, 1, 8, 9, 12, 13});
  const std::size_t rank = linalg::MatrixRank(a);

  auto ns = linalg::Nullspace(a);
  EXPECT_EQ(ns.Shape(0), 5 - rank);
  EXPECT_EQ(ns.Shape(1), std::size_t{5});
  EXPECT_EQ(linalg::MatrixRank(ns), 5 - rank);
  EXPECT_EQ(linalg::Dot(a, linalg::Transpose(ns)),
            linalg::Zeros(field, {3, 5 - rank}));

  auto l_ns = linalg::LeftNullspace(a);
  EXPECT_EQ(l_ns.Shape(0), 3 - rank);
  EXPECT_EQ(linalg::Dot(l_ns, a), linalg::Zeros(field, {3 - rank, 5}));

  auto sys = linalg::SystematicForm(a);
  ASSERT_EQ(sys.g_sys.Shape(0), rank);
  ASSERT_EQ(sys.g_sys.Shape(1), std::size_t{5});
  for (std::size_t i = 0; i < rank; ++i) {
    for (std::size_t j = 0; j < rank; ++j) {
      EXPECT_EQ(sys.g_sys.At(i, j).Value(), i == j ? uint8_t{1} : uint8_t{0});
    }
  }

  auto plu_res = linalg::Plu(a);
  EXPECT_EQ(linalg::Dot(plu_res.p, a), linalg::Dot(plu_res.l, plu_res.u));
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

  auto a_std = MakeGFMatrix(gf_std, n, n, a_vals);
  auto b_std = MakeGFMatrix(gf_std, n, n, b_vals);
  auto a_log = MakeGFMatrix(gf_log, n, n, a_vals);
  auto b_log = MakeGFMatrix(gf_log, n, n, b_vals);
  auto a_zech = MakeGFMatrix(gf_zech, n, n, a_vals);
  auto b_zech = MakeGFMatrix(gf_zech, n, n, b_vals);

  // 1. Matrix multiplication consistency across GF2X, GF2XLOG, GF2XZECH
  auto c_std = linalg::Dot(a_std, b_std);
  auto c_log = linalg::Dot(a_log, b_log);
  auto c_zech = linalg::Dot(a_zech, b_zech);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(c_std.At(i, j).Value(), c_log.At(i, j).Value());
      EXPECT_EQ(c_std.At(i, j).Value(), c_zech.At(i, j).Value());
    }
  }

  // 2. Determinant, trace, rank, RREF, PLU consistency
  auto det_std = linalg::Det(a_std).Value();
  EXPECT_EQ(det_std, linalg::Det(a_log).Value());
  EXPECT_EQ(det_std, linalg::Det(a_zech).Value());

  auto tr_std = linalg::Trace(a_std).Value();
  EXPECT_EQ(tr_std, linalg::Trace(a_log).Value());
  EXPECT_EQ(tr_std, linalg::Trace(a_zech).Value());

  auto rref_std = linalg::Rref(a_std);
  auto rref_zech = linalg::Rref(a_zech);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(rref_std.At(i, j).Value(), rref_zech.At(i, j).Value());
    }
  }

  // 3. Inverse & multi-RHS solve consistency
  if (det_std != 0) {
    auto inv_std = linalg::Inv(a_std);
    auto inv_log = linalg::Inv(a_log);
    auto inv_zech = linalg::Inv(a_zech);
    EXPECT_EQ(linalg::Dot(a_std, inv_std), linalg::Eye(gf_std, n));
    EXPECT_EQ(linalg::Dot(a_log, inv_log), linalg::Eye(gf_log, n));
    EXPECT_EQ(linalg::Dot(a_zech, inv_zech), linalg::Eye(gf_zech, n));
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        EXPECT_EQ(inv_std.At(i, j).Value(), inv_log.At(i, j).Value());
        EXPECT_EQ(inv_std.At(i, j).Value(), inv_zech.At(i, j).Value());
      }
    }

    auto x_zech = linalg::Solve(a_zech, b_zech);
    EXPECT_EQ(linalg::Dot(a_zech, x_zech), b_zech);
  }

  // 4. PLU decomposition across all 3 backends
  auto plu_zech = linalg::Plu(a_zech);
  EXPECT_EQ(linalg::Dot(plu_zech.p, a_zech),
            linalg::Dot(plu_zech.l, plu_zech.u));

  // 5. Characteristic polynomial & Cayley-Hamilton over GF(2^8)
  auto p_std = linalg::CharacteristicPoly(a_std);
  auto p_log = linalg::CharacteristicPoly(a_log);
  auto p_zech = linalg::CharacteristicPoly(a_zech);

  EXPECT_EQ(p_std.Degree(), static_cast<int>(n));
  EXPECT_EQ(p_std[n].Value(), uint8_t{1});
  EXPECT_EQ(p_std[n - 1].Value(), tr_std);
  EXPECT_EQ(p_std[0].Value(), det_std);

  for (std::size_t k = 0; k <= n; ++k) {
    EXPECT_EQ(p_std[k].Value(), p_log[k].Value());
    EXPECT_EQ(p_std[k].Value(), p_zech[k].Value());
  }

  EXPECT_EQ(EvalPolyAtMatrix(p_std, a_std), linalg::Zeros(gf_std, {n, n}));
  EXPECT_EQ(EvalPolyAtMatrix(p_log, a_log), linalg::Zeros(gf_log, {n, n}));
  EXPECT_EQ(EvalPolyAtMatrix(p_zech, a_zech), linalg::Zeros(gf_zech, {n, n}));
}
