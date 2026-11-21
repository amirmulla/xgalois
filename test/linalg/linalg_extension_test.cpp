#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

namespace {

using ExtField = GaloisFieldExtension<uint8_t>;

GFArray<ExtField> MakeExtMatrixFromInts(const std::shared_ptr<ExtField> &field,
                                        std::size_t rows, std::size_t cols,
                                        const std::vector<uint64_t> &vals) {
  GFArray<ExtField>::raw_array_type raw =
      GFArray<ExtField>::raw_array_type::from_shape({rows, cols});
  for (std::size_t i = 0; i < rows * cols; ++i) {
    raw.data()[i] = field->FromInt(vals[i]);
  }
  return GFArray<ExtField>::FromRaw(field, std::move(raw));
}

GFArray<ExtField> MakeExtVectorFromInts(const std::shared_ptr<ExtField> &field,
                                        const std::vector<uint64_t> &vals) {
  GFArray<ExtField>::raw_array_type raw =
      GFArray<ExtField>::raw_array_type::from_shape({vals.size()});
  for (std::size_t i = 0; i < vals.size(); ++i) {
    raw.data()[i] = field->FromInt(vals[i]);
  }
  return GFArray<ExtField>::FromRaw(field, std::move(raw));
}

GFArray<ExtField> EvalPolyAtMatrixExt(const PolynomialDense<ExtField> &p,
                                      const GFArray<ExtField> &a) {
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

void RunFullExtensionLinalgSuite(const std::shared_ptr<ExtField> &field,
                                 uint64_t seed) {
  std::mt19937_64 rng(seed);
  const auto zero = field->AdditiveIdentity();
  const auto one = field->MultiplicativeIdentity();

  // 1. Structural utilities: Eye, Zeros, Ones, Transpose, HStack, VStack
  auto eye_4 = linalg::Eye(field, 4);
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      EXPECT_EQ(eye_4.Raw(i, j), i == j ? one : zero);
    }
  }

  auto a = linalg::RandomMatrix(field, {3, 4}, rng);
  auto a_t = linalg::Transpose(a);
  EXPECT_EQ(a_t.Shape(0), std::size_t{4});
  EXPECT_EQ(a_t.Shape(1), std::size_t{3});
  EXPECT_EQ(linalg::Transpose(a_t), a);

  auto b = linalg::RandomMatrix(field, {3, 2}, rng);
  auto h = linalg::HStack(a, b);
  EXPECT_EQ(h.Shape(0), std::size_t{3});
  EXPECT_EQ(h.Shape(1), std::size_t{6});
  EXPECT_EQ(h.Raw(2, 5), b.Raw(2, 1));

  auto c = linalg::RandomMatrix(field, {2, 4}, rng);
  auto v = linalg::VStack(a, c);
  EXPECT_EQ(v.Shape(0), std::size_t{5});
  EXPECT_EQ(v.Shape(1), std::size_t{4});
  EXPECT_EQ(v.Raw(4, 3), c.Raw(1, 3));

  // 2. Dot (1D*1D, 1D*2D, 2D*1D, 2D*2D), Vdot, Outer, Cross, Kron
  auto u = linalg::RandomMatrix(field, {3}, rng);
  auto vec = linalg::RandomMatrix(field, {3}, rng);
  auto d11 = linalg::Dot(u, vec);
  EXPECT_EQ(d11.RawData()(), linalg::Vdot(u, vec).RawValue());

  auto u_a = linalg::Dot(u, a);
  ASSERT_EQ(u_a.Dimension(), std::size_t{1});
  ASSERT_EQ(u_a.Shape(0), std::size_t{4});
  for (std::size_t j = 0; j < 4; ++j) {
    auto acc = zero;
    for (std::size_t i = 0; i < 3; ++i) {
      acc = field->Add(acc, field->Mul(u.Raw(i), a.Raw(i, j)));
    }
    EXPECT_EQ(u_a.Raw(j), acc);
  }

  auto w = linalg::RandomMatrix(field, {4}, rng);
  auto a_w = linalg::Dot(a, w);
  ASSERT_EQ(a_w.Dimension(), std::size_t{1});
  ASSERT_EQ(a_w.Shape(0), std::size_t{3});
  for (std::size_t i = 0; i < 3; ++i) {
    auto acc = zero;
    for (std::size_t j = 0; j < 4; ++j) {
      acc = field->Add(acc, field->Mul(a.Raw(i, j), w.Raw(j)));
    }
    EXPECT_EQ(a_w.Raw(i), acc);
  }

  auto out_uv = linalg::Outer(u, vec);
  ASSERT_EQ(out_uv.Shape(0), std::size_t{3});
  ASSERT_EQ(out_uv.Shape(1), std::size_t{3});
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(out_uv.Raw(i, j), field->Mul(u.Raw(i), vec.Raw(j)));
    }
  }

  auto cr = linalg::Cross(u, vec);
  ASSERT_EQ(cr.Shape(0), std::size_t{3});
  EXPECT_TRUE(linalg::Vdot(u, cr).IsZero());
  EXPECT_TRUE(linalg::Vdot(vec, cr).IsZero());

  auto k2 = MakeExtMatrixFromInts(field, 2, 2, {1, 2, 3, 4});
  auto k3 = MakeExtMatrixFromInts(field, 2, 2, {2, 0, 5, 1});
  auto kron_res = linalg::Kron(k2, k3);
  ASSERT_EQ(kron_res.Shape(0), std::size_t{4});
  ASSERT_EQ(kron_res.Shape(1), std::size_t{4});
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t col = 0; col < 2; ++col) {
          EXPECT_EQ(kron_res.Raw(i * 2 + r, j * 2 + col),
                    field->Mul(k2.Raw(i, j), k3.Raw(r, col)));
        }
      }
    }
  }

  // 3. Invertible 4x4 matrix: Det, Trace, Inv, Solve (1D & 2D), MatrixPower
  GFArray<ExtField> m;
  for (int attempt = 0; attempt < 100; ++attempt) {
    auto candidate = linalg::RandomMatrix(field, {4, 4}, rng);
    if (!linalg::Det(candidate).IsZero()) {
      m = std::move(candidate);
      break;
    }
  }
  ASSERT_TRUE(m.HasField());
  ASSERT_FALSE(linalg::Det(m).IsZero());

  // Multiplicativity of determinant: det(M * N) == det(M) * det(N)
  auto n_sq = linalg::RandomMatrix(field, {4, 4}, rng);
  EXPECT_EQ(linalg::Det(linalg::Dot(m, n_sq)),
            linalg::Det(m) * linalg::Det(n_sq));

  // Zero-pivot row swap sign flip: [[0, a], [b, c]] has det = -a*b
  auto alpha = field->MultiplicativeGenerator();
  auto alpha2 = field->Mul(alpha, alpha);
  auto zp = GFArray<ExtField>::Zeros(field, {2, 2});
  zp.Raw(0, 1) = alpha;
  zp.Raw(1, 0) = alpha2;
  zp.Raw(1, 1) = one;
  EXPECT_EQ(linalg::Det(zp).RawValue(), field->Neg(field->Mul(alpha, alpha2)));

  // Inversion
  auto m_inv = linalg::Inv(m);
  EXPECT_EQ(linalg::Dot(m, m_inv), eye_4);
  EXPECT_EQ(linalg::Dot(m_inv, m), eye_4);

  // Matrix powers (0, positive, negative)
  EXPECT_EQ(linalg::MatrixPower(m, 0), eye_4);
  EXPECT_EQ(linalg::MatrixPower(m, 1), m);
  EXPECT_EQ(linalg::MatrixPower(m, -1), m_inv);
  EXPECT_EQ(linalg::Dot(linalg::MatrixPower(m, 5), linalg::MatrixPower(m, -5)),
            eye_4);

  // Solve 1D & multi-RHS 2D
  auto rhs1 = linalg::RandomMatrix(field, {4}, rng);
  auto sol1 = linalg::Solve(m, rhs1);
  EXPECT_EQ(linalg::Dot(m, sol1), rhs1);

  auto rhs2 = linalg::RandomMatrix(field, {4, 3}, rng);
  auto sol2 = linalg::Solve(m, rhs2);
  EXPECT_EQ(linalg::Dot(m, sol2), rhs2);

  // Singular matrix throws
  auto sing = GFArray<ExtField>::Zeros(field, {3, 3});
  sing.Raw(0, 0) = alpha;
  sing.Raw(1, 0) = alpha2;
  EXPECT_TRUE(linalg::Det(sing).IsZero());
  EXPECT_THROW(linalg::Inv(sing), std::runtime_error);
  EXPECT_THROW(linalg::Solve(sing, linalg::Ones(field, {3})),
               std::runtime_error);

  // 4. Rref, MatrixRank, Nullspace, LeftNullspace, RowSpace, ColumnSpace,
  //    SystematicForm
  auto rect = linalg::RandomMatrix(field, {3, 5}, rng);
  const std::size_t rank = linalg::MatrixRank(rect);
  auto e = linalg::RowEchelon(rect);
  auto r = linalg::Rref(rect);
  EXPECT_EQ(linalg::MatrixRank(e), rank);
  EXPECT_EQ(linalg::MatrixRank(r), rank);

  auto ns = linalg::Nullspace(rect);
  ASSERT_EQ(ns.Shape(0), 5 - rank);
  ASSERT_EQ(ns.Shape(1), std::size_t{5});
  EXPECT_EQ(linalg::MatrixRank(ns), 5 - rank);
  EXPECT_EQ(linalg::Dot(rect, linalg::Transpose(ns)),
            linalg::Zeros(field, {3, 5 - rank}));

  auto lns = linalg::LeftNullspace(rect);
  ASSERT_EQ(lns.Shape(0), 3 - rank);
  EXPECT_EQ(linalg::Dot(lns, rect), linalg::Zeros(field, {3 - rank, 5}));

  EXPECT_EQ(linalg::MatrixRank(linalg::RowSpace(rect)), rank);
  EXPECT_EQ(linalg::MatrixRank(linalg::ColumnSpace(rect)), rank);

  auto sys = linalg::SystematicForm(rect);
  ASSERT_EQ(sys.pivot_cols.size(), rank);
  ASSERT_EQ(sys.g_sys.Shape(0), rank);
  ASSERT_EQ(sys.g_sys.Shape(1), std::size_t{5});
  for (std::size_t i = 0; i < rank; ++i) {
    for (std::size_t j = 0; j < rank; ++j) {
      EXPECT_EQ(sys.g_sys.Raw(i, j), i == j ? one : zero);
    }
  }

  // 5. PLU decomposition: P * A == L * U
  auto plu_res = linalg::Plu(m);
  EXPECT_EQ(linalg::Dot(plu_res.p, m), linalg::Dot(plu_res.l, plu_res.u));
  auto diag_prod = one;
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(plu_res.l.Raw(i, i), one);
    diag_prod = field->Mul(diag_prod, plu_res.u.Raw(i, i));
  }
  if (plu_res.num_swaps % 2 == 1) {
    diag_prod = field->Neg(diag_prod);
  }
  EXPECT_EQ(diag_prod, linalg::Det(m).RawValue());

  auto plu_rect = linalg::Plu(rect);
  EXPECT_EQ(linalg::Dot(plu_rect.p, rect), linalg::Dot(plu_rect.l, plu_rect.u));

  // 6. Characteristic polynomial & Cayley-Hamilton
  auto p = linalg::CharacteristicPoly(m);
  EXPECT_EQ(p.Degree(), 4);
  EXPECT_EQ(p[4].RawValue(), one);
  EXPECT_EQ(p[3], -linalg::Trace(m));
  EXPECT_EQ(p[0], linalg::Det(m));  // (-1)^4 == +1
  EXPECT_EQ(EvalPolyAtMatrixExt(p, m), linalg::Zeros(field, {4, 4}));

  // 3x3 odd dimension constant term sign check: (-1)^3 * det == -det
  auto m3 = linalg::RandomMatrix(field, {3, 3}, rng);
  auto p3 = linalg::CharacteristicPoly(m3);
  EXPECT_EQ(p3.Degree(), 3);
  EXPECT_EQ(p3[3].RawValue(), one);
  EXPECT_EQ(p3[2], -linalg::Trace(m3));
  EXPECT_EQ(p3[0], -linalg::Det(m3));
  EXPECT_EQ(EvalPolyAtMatrixExt(p3, m3), linalg::Zeros(field, {3, 3}));
}

}  // namespace

TEST(LinalgExtensionTest, GF3SquaredFullSuite) {
  auto field = std::make_shared<ExtField>(std::make_pair(3ULL, 2ULL));
  ASSERT_EQ(field->Characteristic(), uint64_t{3});
  ASSERT_EQ(field->Order(), uint64_t{9});
  RunFullExtensionLinalgSuite(field, 1001);
}

TEST(LinalgExtensionTest, GF5SquaredFullSuite) {
  auto field = std::make_shared<ExtField>(std::make_pair(5ULL, 2ULL));
  ASSERT_EQ(field->Characteristic(), uint64_t{5});
  ASSERT_EQ(field->Order(), uint64_t{25});
  RunFullExtensionLinalgSuite(field, 2002);
}

TEST(LinalgExtensionTest, ExplicitSmallGF9ArithmeticVerification) {
  auto field = std::make_shared<ExtField>(std::make_pair(3ULL, 2ULL));
  // Test explicit 2x2 matrix over GF(3^2) with non-base-field polynomial
  // elements (e.g., alpha, alpha + 1, 2*alpha + 1).
  auto a = MakeExtMatrixFromInts(field, 2, 2, {3, 4, 7, 5});
  auto b = MakeExtVectorFromInts(field, {5, 8});

  // det(A) = A(0,0)*A(1,1) - A(0,1)*A(1,0)
  auto expected_det = field->Sub(field->Mul(a.Raw(0, 0), a.Raw(1, 1)),
                                 field->Mul(a.Raw(0, 1), a.Raw(1, 0)));
  EXPECT_EQ(linalg::Det(a).RawValue(), expected_det);
  ASSERT_FALSE(field->IsZero(expected_det));

  auto x = linalg::Solve(a, b);
  EXPECT_EQ(linalg::Dot(a, x), b);

  auto a_inv = linalg::Inv(a);
  EXPECT_EQ(linalg::Dot(a, a_inv), linalg::Eye(field, 2));
}
