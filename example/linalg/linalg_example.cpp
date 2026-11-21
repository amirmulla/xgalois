#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

namespace {

template <typename Field>
void PrintMatrix(const std::string& label, const GFArray<Field>& mat) {
  std::cout << label << " (" << mat.Shape(0) << "x" << mat.Shape(1) << "):\n";
  for (size_t i = 0; i < mat.Shape(0); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < mat.Shape(1); ++j) {
      std::cout << mat.At(i, j);
      if (j + 1 < mat.Shape(1)) std::cout << ", ";
    }
    std::cout << " ]\n";
  }
}

template <typename Field>
void PrintVector(const std::string& label, const GFArray<Field>& vec) {
  std::cout << label << ": [ ";
  for (size_t i = 0; i < vec.Size(); ++i) {
    std::cout << vec.At(i);
    if (i + 1 < vec.Size()) std::cout << ", ";
  }
  std::cout << " ]\n";
}

template <typename Field>
GFArray<Field> EvalPolyAtMatrix(const PolynomialDense<Field>& p,
                                const GFArray<Field>& A) {
  const auto field = A.Field();
  const size_t n = A.Shape(0);
  auto acc = GFArray<Field>::Zeros(field, {n, n});
  if (p.IsZero()) {
    return acc;
  }
  const auto I = linalg::eye(field, n);
  for (int d = p.Degree(); d >= 0; --d) {
    acc = linalg::dot(acc, A);
    const auto c = p[static_cast<size_t>(d)];
    if (!c.IsZero()) {
      acc += I * c;
    }
  }
  return acc;
}

}  // namespace

int main() {
  std::cout << "XGalois Linear Algebra Examples\n";
  std::cout << "===============================\n";

  using PrimeField = GFP<uint8_t>;
  auto gf7 = std::make_shared<PrimeField>(7);

  std::cout << "\n=== Contiguous GFArray & xg::linalg over GF(7) ===\n";
  GFArray<PrimeField> A(gf7,
                        xt::xarray<uint8_t>{{1, 2, 3}, {2, 5, 1}, {3, 1, 4}});
  GFArray<PrimeField> B(gf7,
                        xt::xarray<uint8_t>{{2, 0, 1}, {1, 3, 2}, {4, 1, 5}});
  GFArray<PrimeField> b(gf7, xt::xarray<uint8_t>{4, 6, 2});

  PrintMatrix("A", A);
  PrintMatrix("B", B);

  const auto AB = linalg::dot(A, B);
  PrintMatrix("A * B (linalg::dot)", AB);

  std::cout << "\n--- Determinant, Trace, Rank, Inverse & Solve ---\n";
  const auto det_a = linalg::det(A);
  const auto tr_a = linalg::trace(A);
  const size_t rank_a = linalg::matrix_rank(A);
  const auto A_inv = linalg::inv(A);
  const auto x = linalg::solve(A, b);

  std::cout << "  det(A)   = " << det_a << "\n";
  std::cout << "  trace(A) = " << tr_a << "\n";
  std::cout << "  rank(A)  = " << rank_a << "\n";
  PrintVector("  Solution x to A*x = b", x);
  std::cout << "  Verify A * x == b:     " << std::boolalpha
            << (linalg::dot(A, x) == b) << "\n";
  std::cout << "  Verify A * inv(A) == I: "
            << (linalg::dot(A, A_inv) == linalg::eye(gf7, 3)) << "\n";

  std::cout << "\n--- RREF & Nullspace of a Rank-Deficient Matrix ---\n";
  // Row 2 = Row 0 + Row 1 (mod 7) -> rank 2, nullity 2
  GFArray<PrimeField> M(
      gf7, xt::xarray<uint8_t>{{1, 2, 3, 4}, {2, 4, 1, 3}, {3, 6, 4, 0}});
  const auto R = linalg::rref(M);
  const auto N = linalg::nullspace(M);
  const auto M_NT = linalg::dot(M, linalg::transpose(N));

  PrintMatrix("RREF(M)", R);
  PrintMatrix("Nullspace basis N", N);
  std::cout << "  rank(M) = " << linalg::matrix_rank(M)
            << ", nullity(M) = " << N.Shape(0) << "\n";
  std::cout << "  Verify M * N^T == 0: "
            << (M_NT == GFArray<PrimeField>::Zeros(gf7, {3, 2})) << "\n";

  std::cout << "\n--- PLU Decomposition ---\n";
  const auto plu_res = linalg::plu(A);
  const auto PA = linalg::dot(plu_res.P, A);
  const auto LU = linalg::dot(plu_res.L, plu_res.U);
  std::cout << "  Row swaps = " << plu_res.num_swaps
            << ", verify P * A == L * U: " << (PA == LU) << "\n";

  std::cout << "\n--- Characteristic Polynomial (Berkowitz Algorithm) ---\n";
  const auto char_poly = linalg::characteristic_poly(A);
  const auto poly_at_A = EvalPolyAtMatrix(char_poly, A);
  std::cout << "  chi_A(x) = det(xI - A) = " << char_poly.ToString() << "\n";
  std::cout << "  Verify Cayley-Hamilton chi_A(A) == 0: "
            << (poly_at_A == GFArray<PrimeField>::Zeros(gf7, {3, 3})) << "\n";

  std::cout << "\n=== Bit-Packed GF(2) M4RI Linear Algebra ===\n";
  auto gf2 = std::make_shared<GF2>();
  std::mt19937_64 rng(42);
  constexpr size_t kBinSize = 64;

  // Construct a guaranteed-invertible 64x64 GF(2) matrix G = L * U
  auto L_bin = linalg::eye(gf2, kBinSize);
  auto U_bin = linalg::eye(gf2, kBinSize);
  for (size_t i = 0; i < kBinSize; ++i) {
    for (size_t j = 0; j < i; ++j) {
      L_bin.Raw(i, j) = gf2->Random(rng);
    }
    for (size_t j = i + 1; j < kBinSize; ++j) {
      U_bin.Raw(i, j) = gf2->Random(rng);
    }
  }
  const auto G_bin = linalg::dot(L_bin, U_bin);
  const auto rhs_bin = GFArray<GF2>::Random(gf2, {kBinSize}, rng);

  const auto det_bin = linalg::det(G_bin);
  const auto G_bin_inv = linalg::inv(G_bin);
  const auto x_bin = linalg::solve(G_bin, rhs_bin);
  const auto plu_bin = linalg::plu(G_bin);

  std::cout << "  64x64 GF(2) matrix G = L * U:\n";
  std::cout << "    det(G)              = " << det_bin << "\n";
  std::cout << "    rank(G)             = " << linalg::matrix_rank(G_bin)
            << "\n";
  std::cout << "    G * inv(G) == I_64  : "
            << (linalg::dot(G_bin, G_bin_inv) == linalg::eye(gf2, kBinSize))
            << "\n";
  std::cout << "    G * solve(G,b) == b : "
            << (linalg::dot(G_bin, x_bin) == rhs_bin) << "\n";
  std::cout << "    P * G == L * U      : "
            << (linalg::dot(plu_bin.P, G_bin) ==
                linalg::dot(plu_bin.L, plu_bin.U))
            << "\n";

  std::cout << "\nLinear algebra examples completed successfully!\n";
  return 0;
}
