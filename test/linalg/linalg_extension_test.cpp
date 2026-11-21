#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

namespace {

using ExtField = GaloisFieldExtension<uint8_t>;
using ExtElem = GaloisFieldElement<ExtField>;

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
                                      const GFArray<ExtField> &A) {
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

void RunFullExtensionLinalgSuite(const std::shared_ptr<ExtField> &field,
                                 uint64_t seed) {
  std::mt19937_64 rng(seed);
  const auto zero = field->AdditiveIdentity();
  const auto one = field->MultiplicativeIdentity();

  // 1. Structural utilities: eye, zeros, ones, transpose, hstack, vstack
  auto I4 = linalg::eye(field, 4);
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      EXPECT_EQ(I4.Raw(i, j), i == j ? one : zero);
    }
  }

  auto A = linalg::random_matrix(field, {3, 4}, rng);
  auto AT = linalg::transpose(A);
  EXPECT_EQ(AT.Shape(0), 4);
  EXPECT_EQ(AT.Shape(1), 3);
  EXPECT_EQ(linalg::transpose(AT), A);

  auto B = linalg::random_matrix(field, {3, 2}, rng);
  auto H = linalg::hstack(A, B);
  EXPECT_EQ(H.Shape(0), 3);
  EXPECT_EQ(H.Shape(1), 6);
  EXPECT_EQ(H.Raw(2, 5), B.Raw(2, 1));

  auto C = linalg::random_matrix(field, {2, 4}, rng);
  auto V = linalg::vstack(A, C);
  EXPECT_EQ(V.Shape(0), 5);
  EXPECT_EQ(V.Shape(1), 4);
  EXPECT_EQ(V.Raw(4, 3), C.Raw(1, 3));

  // 2. dot (1D*1D, 1D*2D, 2D*1D, 2D*2D), vdot, outer, cross, kron
  auto u = linalg::random_matrix(field, {3}, rng);
  auto v = linalg::random_matrix(field, {3}, rng);
  auto d11 = linalg::dot(u, v);
  EXPECT_EQ(d11.RawData()(), linalg::vdot(u, v).RawValue());

  auto uA = linalg::dot(u, A);
  ASSERT_EQ(uA.Dimension(), 1);
  ASSERT_EQ(uA.Shape(0), 4);
  for (std::size_t j = 0; j < 4; ++j) {
    auto acc = zero;
    for (std::size_t i = 0; i < 3; ++i) {
      acc = field->Add(acc, field->Mul(u.Raw(i), A.Raw(i, j)));
    }
    EXPECT_EQ(uA.Raw(j), acc);
  }

  auto w = linalg::random_matrix(field, {4}, rng);
  auto Aw = linalg::dot(A, w);
  ASSERT_EQ(Aw.Dimension(), 1);
  ASSERT_EQ(Aw.Shape(0), 3);
  for (std::size_t i = 0; i < 3; ++i) {
    auto acc = zero;
    for (std::size_t j = 0; j < 4; ++j) {
      acc = field->Add(acc, field->Mul(A.Raw(i, j), w.Raw(j)));
    }
    EXPECT_EQ(Aw.Raw(i), acc);
  }

  auto out_uv = linalg::outer(u, v);
  ASSERT_EQ(out_uv.Shape(0), 3);
  ASSERT_EQ(out_uv.Shape(1), 3);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(out_uv.Raw(i, j), field->Mul(u.Raw(i), v.Raw(j)));
    }
  }

  auto cr = linalg::cross(u, v);
  ASSERT_EQ(cr.Shape(0), 3);
  EXPECT_TRUE(linalg::vdot(u, cr).IsZero());
  EXPECT_TRUE(linalg::vdot(v, cr).IsZero());

  auto K2 = MakeExtMatrixFromInts(field, 2, 2, {1, 2, 3, 4});
  auto K3 = MakeExtMatrixFromInts(field, 2, 2, {2, 0, 5, 1});
  auto kron_res = linalg::kron(K2, K3);
  ASSERT_EQ(kron_res.Shape(0), 4);
  ASSERT_EQ(kron_res.Shape(1), 4);
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
      for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t c = 0; c < 2; ++c) {
          EXPECT_EQ(kron_res.Raw(i * 2 + r, j * 2 + c),
                    field->Mul(K2.Raw(i, j), K3.Raw(r, c)));
        }
      }
    }
  }

  // 3. Invertible 4x4 matrix: det, trace, inv, solve (1D & 2D), matrix_power
  GFArray<ExtField> M;
  for (int attempt = 0; attempt < 100; ++attempt) {
    auto candidate = linalg::random_matrix(field, {4, 4}, rng);
    if (!linalg::det(candidate).IsZero()) {
      M = std::move(candidate);
      break;
    }
  }
  ASSERT_TRUE(M.HasField());
  ASSERT_FALSE(linalg::det(M).IsZero());

  // Multiplicativity of determinant: det(M * N) == det(M) * det(N)
  auto N_sq = linalg::random_matrix(field, {4, 4}, rng);
  EXPECT_EQ(linalg::det(linalg::dot(M, N_sq)),
            linalg::det(M) * linalg::det(N_sq));

  // Zero-pivot row swap sign flip: [[0, a], [b, c]] has det = -a*b
  auto alpha = field->MultiplicativeGenerator();
  auto alpha2 = field->Mul(alpha, alpha);
  auto ZP = GFArray<ExtField>::Zeros(field, {2, 2});
  ZP.Raw(0, 1) = alpha;
  ZP.Raw(1, 0) = alpha2;
  ZP.Raw(1, 1) = one;
  EXPECT_EQ(linalg::det(ZP).RawValue(), field->Neg(field->Mul(alpha, alpha2)));

  // Inversion
  auto M_inv = linalg::inv(M);
  EXPECT_EQ(linalg::dot(M, M_inv), I4);
  EXPECT_EQ(linalg::dot(M_inv, M), I4);

  // Matrix powers (0, positive, negative)
  EXPECT_EQ(linalg::matrix_power(M, 0), I4);
  EXPECT_EQ(linalg::matrix_power(M, 1), M);
  EXPECT_EQ(linalg::matrix_power(M, -1), M_inv);
  EXPECT_EQ(
      linalg::dot(linalg::matrix_power(M, 5), linalg::matrix_power(M, -5)), I4);

  // Solve 1D & multi-RHS 2D
  auto rhs1 = linalg::random_matrix(field, {4}, rng);
  auto sol1 = linalg::solve(M, rhs1);
  EXPECT_EQ(linalg::dot(M, sol1), rhs1);

  auto rhs2 = linalg::random_matrix(field, {4, 3}, rng);
  auto sol2 = linalg::solve(M, rhs2);
  EXPECT_EQ(linalg::dot(M, sol2), rhs2);

  // Singular matrix throws
  auto Sing = GFArray<ExtField>::Zeros(field, {3, 3});
  Sing.Raw(0, 0) = alpha;
  Sing.Raw(1, 0) = alpha2;
  EXPECT_TRUE(linalg::det(Sing).IsZero());
  EXPECT_THROW(linalg::inv(Sing), std::runtime_error);
  EXPECT_THROW(linalg::solve(Sing, linalg::ones(field, {3})),
               std::runtime_error);

  // 4. RREF, rank, nullspace, left_nullspace, row_space, column_space,
  //    systematic_form
  auto Rect = linalg::random_matrix(field, {3, 5}, rng);
  const std::size_t rank = linalg::matrix_rank(Rect);
  auto E = linalg::row_echelon(Rect);
  auto R = linalg::rref(Rect);
  EXPECT_EQ(linalg::matrix_rank(E), rank);
  EXPECT_EQ(linalg::matrix_rank(R), rank);

  auto NS = linalg::nullspace(Rect);
  ASSERT_EQ(NS.Shape(0), 5 - rank);
  ASSERT_EQ(NS.Shape(1), 5);
  EXPECT_EQ(linalg::matrix_rank(NS), 5 - rank);
  EXPECT_EQ(linalg::dot(Rect, linalg::transpose(NS)),
            linalg::zeros(field, {3, 5 - rank}));

  auto LNS = linalg::left_nullspace(Rect);
  ASSERT_EQ(LNS.Shape(0), 3 - rank);
  EXPECT_EQ(linalg::dot(LNS, Rect), linalg::zeros(field, {3 - rank, 5}));

  EXPECT_EQ(linalg::matrix_rank(linalg::row_space(Rect)), rank);
  EXPECT_EQ(linalg::matrix_rank(linalg::column_space(Rect)), rank);

  auto sys = linalg::systematic_form(Rect);
  ASSERT_EQ(sys.pivot_cols.size(), rank);
  ASSERT_EQ(sys.G_sys.Shape(0), rank);
  ASSERT_EQ(sys.G_sys.Shape(1), 5);
  for (std::size_t i = 0; i < rank; ++i) {
    for (std::size_t j = 0; j < rank; ++j) {
      EXPECT_EQ(sys.G_sys.Raw(i, j), i == j ? one : zero);
    }
  }

  // 5. PLU decomposition: P * A == L * U
  auto plu_res = linalg::plu(M);
  EXPECT_EQ(linalg::dot(plu_res.P, M), linalg::dot(plu_res.L, plu_res.U));
  auto diag_prod = one;
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(plu_res.L.Raw(i, i), one);
    diag_prod = field->Mul(diag_prod, plu_res.U.Raw(i, i));
  }
  if (plu_res.num_swaps % 2 == 1) {
    diag_prod = field->Neg(diag_prod);
  }
  EXPECT_EQ(diag_prod, linalg::det(M).RawValue());

  auto plu_rect = linalg::plu(Rect);
  EXPECT_EQ(linalg::dot(plu_rect.P, Rect), linalg::dot(plu_rect.L, plu_rect.U));

  // 6. Characteristic polynomial & Cayley-Hamilton
  auto p = linalg::characteristic_poly(M);
  EXPECT_EQ(p.Degree(), 4);
  EXPECT_EQ(p[4].RawValue(), one);
  EXPECT_EQ(p[3], -linalg::trace(M));
  EXPECT_EQ(p[0], linalg::det(M));  // (-1)^4 == +1
  EXPECT_EQ(EvalPolyAtMatrixExt(p, M), linalg::zeros(field, {4, 4}));

  // 3x3 odd dimension constant term sign check: (-1)^3 * det == -det
  auto M3 = linalg::random_matrix(field, {3, 3}, rng);
  auto p3 = linalg::characteristic_poly(M3);
  EXPECT_EQ(p3.Degree(), 3);
  EXPECT_EQ(p3[3].RawValue(), one);
  EXPECT_EQ(p3[2], -linalg::trace(M3));
  EXPECT_EQ(p3[0], -linalg::det(M3));
  EXPECT_EQ(EvalPolyAtMatrixExt(p3, M3), linalg::zeros(field, {3, 3}));

  // 7. Legacy garray API compatibility over GaloisFieldExtension (Bug B-15)
  auto M_ga = M.ToElementArray();
  auto rhs1_ga = rhs1.ToElementArray();
  auto sol_ga = linalg::solve(M_ga, rhs1_ga);
  auto check_ga = linalg::dot(M_ga, sol_ga);
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(check_ga(i), rhs1_ga(i));
  }
  auto M_inv_ga = linalg::inv(M_ga);
  auto prod_ga = linalg::dot(M_ga, M_inv_ga);
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      EXPECT_EQ(prod_ga(i, j).RawValue(), i == j ? one : zero);
    }
  }
  auto p_ga = linalg::characteristic_poly(M_ga);
  EXPECT_EQ(p_ga, p);
}

}  // namespace

TEST(LinalgExtensionTest, GF3SquaredFullSuite) {
  auto field = std::make_shared<ExtField>(std::make_pair(3ULL, 2ULL));
  ASSERT_EQ(field->Characteristic(), 3);
  ASSERT_EQ(field->Order(), 9);
  RunFullExtensionLinalgSuite(field, 1001);
}

TEST(LinalgExtensionTest, GF5SquaredFullSuite) {
  auto field = std::make_shared<ExtField>(std::make_pair(5ULL, 2ULL));
  ASSERT_EQ(field->Characteristic(), 5);
  ASSERT_EQ(field->Order(), 25);
  RunFullExtensionLinalgSuite(field, 2002);
}

TEST(LinalgExtensionTest, ExplicitSmallGF9ArithmeticVerification) {
  auto field = std::make_shared<ExtField>(std::make_pair(3ULL, 2ULL));
  // Test explicit 2x2 matrix over GF(3^2) with non-base-field polynomial
  // elements (e.g., alpha, alpha + 1, 2*alpha + 1).
  auto A = MakeExtMatrixFromInts(field, 2, 2, {3, 4, 7, 5});
  auto b = MakeExtVectorFromInts(field, {5, 8});

  // det(A) = A(0,0)*A(1,1) - A(0,1)*A(1,0)
  auto expected_det = field->Sub(field->Mul(A.Raw(0, 0), A.Raw(1, 1)),
                                 field->Mul(A.Raw(0, 1), A.Raw(1, 0)));
  EXPECT_EQ(linalg::det(A).RawValue(), expected_det);
  ASSERT_FALSE(field->IsZero(expected_det));

  auto x = linalg::solve(A, b);
  EXPECT_EQ(linalg::dot(A, x), b);

  auto A_inv = linalg::inv(A);
  EXPECT_EQ(linalg::dot(A, A_inv), linalg::eye(field, 2));
}
