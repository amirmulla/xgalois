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
  const auto I = linalg::Eye(field, n);
  for (int d = p.Degree(); d >= 0; --d) {
    acc = linalg::Dot(acc, A);
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

  const auto AB = linalg::Dot(A, B);
  PrintMatrix("A * B (linalg::Dot)", AB);

  std::cout << "\n--- Determinant, Trace, Rank, Inverse, Power & Solve ---\n";
  const auto det_a = linalg::Det(A);
  const auto tr_a = linalg::Trace(A);
  const size_t rank_a = linalg::MatrixRank(A);
  const auto A_inv = linalg::Inv(A);
  const auto A_cube = linalg::MatrixPower(A, 3);
  const auto x = linalg::Solve(A, b);

  std::cout << "  Det(A)   = " << det_a << "\n";
  std::cout << "  Trace(A) = " << tr_a << "\n";
  std::cout << "  Rank(A)  = " << rank_a << "\n";
  PrintVector("  Solution x to A*x = b", x);
  std::cout << "  Verify A * x == b:             " << std::boolalpha
            << (linalg::Dot(A, x) == b) << "\n";
  std::cout << "  Verify A * Inv(A) == I:        "
            << (linalg::Dot(A, A_inv) == linalg::Eye(gf7, 3)) << "\n";
  std::cout << "  Verify MatrixPower(A,3)==A^3:  "
            << (A_cube == linalg::Dot(linalg::Dot(A, A), A)) << "\n";
  std::cout << "  Verify MatrixPower(A,-1)==A^-1:"
            << (linalg::MatrixPower(A, -1) == A_inv) << "\n";

  std::cout << "\n--- RREF & Nullspace of a Rank-Deficient Matrix ---\n";
  // Row 2 = Row 0 + Row 1 (mod 7) -> rank 2, nullity 2
  GFArray<PrimeField> M(
      gf7, xt::xarray<uint8_t>{{1, 2, 3, 4}, {2, 4, 1, 3}, {3, 6, 4, 0}});
  const auto R = linalg::Rref(M);
  const auto N = linalg::Nullspace(M);
  const auto M_NT = linalg::Dot(M, linalg::Transpose(N));

  PrintMatrix("Rref(M)", R);
  PrintMatrix("Nullspace basis N", N);
  std::cout << "  MatrixRank(M) = " << linalg::MatrixRank(M)
            << ", nullity(M) = " << N.Shape(0) << "\n";
  std::cout << "  Verify M * N^T == 0: "
            << (M_NT == GFArray<PrimeField>::Zeros(gf7, {3, 2})) << "\n";

  std::cout << "\n--- PLU Decomposition ---\n";
  const auto plu_res = linalg::Plu(A);
  const auto PA = linalg::Dot(plu_res.p, A);
  const auto LU = linalg::Dot(plu_res.l, plu_res.u);
  std::cout << "  Row swaps = " << plu_res.num_swaps
            << ", verify P * A == L * U: " << (PA == LU) << "\n";

  std::cout << "\n--- Characteristic Polynomial (Berkowitz Algorithm) ---\n";
  const auto char_poly = linalg::CharacteristicPoly(A);
  const auto poly_at_A = EvalPolyAtMatrix(char_poly, A);
  std::cout << "  chi_A(x) = det(xI - A) = " << char_poly.ToString() << "\n";
  std::cout << "  Verify Cayley-Hamilton chi_A(A) == 0: "
            << (poly_at_A == GFArray<PrimeField>::Zeros(gf7, {3, 3})) << "\n";

  std::cout << "\n=== Bit-Packed GF(2) M4RI Linear Algebra ===\n";
  auto gf2 = std::make_shared<GF2>();
  std::mt19937_64 rng(42);
  constexpr size_t kBinSize = 64;

  // Construct a guaranteed-invertible 64x64 GF(2) matrix G = L * U
  auto L_bin = linalg::Eye(gf2, kBinSize);
  auto U_bin = linalg::Eye(gf2, kBinSize);
  for (size_t i = 0; i < kBinSize; ++i) {
    for (size_t j = 0; j < i; ++j) {
      L_bin.Raw(i, j) = gf2->Random(rng);
    }
    for (size_t j = i + 1; j < kBinSize; ++j) {
      U_bin.Raw(i, j) = gf2->Random(rng);
    }
  }
  const auto G_bin = linalg::Dot(L_bin, U_bin);
  const auto rhs_bin = GFArray<GF2>::Random(gf2, {kBinSize}, rng);

  const auto det_bin = linalg::Det(G_bin);
  const auto G_bin_inv = linalg::Inv(G_bin);
  const auto x_bin = linalg::Solve(G_bin, rhs_bin);
  const auto plu_bin = linalg::Plu(G_bin);

  std::cout << "  64x64 GF(2) matrix G = L * U:\n";
  std::cout << "    Det(G)              = " << det_bin << "\n";
  std::cout << "    MatrixRank(G)       = " << linalg::MatrixRank(G_bin)
            << "\n";
  std::cout << "    G * Inv(G) == I_64  : "
            << (linalg::Dot(G_bin, G_bin_inv) == linalg::Eye(gf2, kBinSize))
            << "\n";
  std::cout << "    G * Solve(G,b) == b : "
            << (linalg::Dot(G_bin, x_bin) == rhs_bin) << "\n";
  std::cout << "    P * G == L * U      : "
            << (linalg::Dot(plu_bin.p, G_bin) ==
                linalg::Dot(plu_bin.l, plu_bin.u))
            << "\n";

  std::cout << "\nLinear algebra examples completed successfully!\n";
  return 0;
}
