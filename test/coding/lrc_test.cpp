#include "xgalois/coding/lrc.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace {

using namespace xg;
using namespace xg::coding;

template <typename GaloisField>
bool IsZeroMatrix(const GFMatrix<GaloisField> &mat) {
  const auto field = mat.Field();
  for (const auto &v : mat.RawSpan()) {
    if (!field->IsZero(v)) {
      return false;
    }
  }
  return true;
}

template <typename GaloisField>
typename GaloisField::element_type RandomNonZero(
    const std::shared_ptr<GaloisField> &field, std::mt19937_64 &rng) {
  auto v = field->Random(rng);
  while (field->IsZero(v)) {
    v = field->Random(rng);
  }
  return v;
}

template <typename GaloisField>
GFArray<GaloisField> ErasePositions(
    const GFArray<GaloisField> &cw,
    const std::vector<std::size_t> &erasure_positions) {
  auto field = cw.Field();
  auto received = cw;
  const auto zero = field->AdditiveIdentity();
  for (std::size_t pos : erasure_positions) {
    received.Raw(pos) = zero;
  }
  return received;
}

template <typename GaloisField>
GFArray<GaloisField> CorruptPositions(const GFArray<GaloisField> &cw,
                                      const std::vector<std::size_t> &positions,
                                      std::mt19937_64 &rng) {
  auto field = cw.Field();
  auto received = cw;
  for (std::size_t pos : positions) {
    const auto err = RandomNonZero(field, rng);
    received.Raw(pos) = field->Add(received.Raw(pos), err);
  }
  return received;
}

template <typename GaloisField>
std::size_t ExhaustiveMinimumDistance(
    const LocallyRecoverableCode<GaloisField> &code) {
  const auto codewords = code.GetCodewords(65536);
  std::size_t min_w = code.Length();
  for (const auto &cw : codewords) {
    const std::size_t w = code.Weight(cw);
    if (w > 0 && w < min_w) {
      min_w = w;
    }
  }
  return min_w;
}

template <typename GaloisField>
void VerifyLinearAlgebraInvariants(
    const LocallyRecoverableCode<GaloisField> &code, uint64_t seed) {
  const auto &G = code.GeneratorMatrix();
  const auto &H = code.ParityCheckMatrix();
  EXPECT_EQ(G.Shape(0), code.Dimension());
  EXPECT_EQ(G.Shape(1), code.Length());
  EXPECT_EQ(H.Shape(0), code.Redundancy());
  EXPECT_EQ(H.Shape(1), code.Length());
  EXPECT_EQ(linalg::MatrixRank(G), code.Dimension());
  if (code.Redundancy() > 0) {
    EXPECT_EQ(linalg::MatrixRank(H), code.Redundancy());
    EXPECT_TRUE(IsZeroMatrix(linalg::Dot(G, linalg::Transpose(H))));
  }

  std::mt19937_64 rng(seed);
  for (int trial = 0; trial < 5; ++trial) {
    auto msg =
        GFArray<GaloisField>::Random(code.Field(), {code.Dimension()}, rng);
    auto cw = code.Encode(msg);
    EXPECT_TRUE(code.Contains(cw));
    auto recovered_msg = code.Unencode(cw);
    EXPECT_EQ(recovered_msg, msg);
  }
}

template <typename GaloisField>
void VerifyStrictPairwiseDisjointRepairSets(
    const LocallyRecoverableCode<GaloisField> &code) {
  const std::size_t n = code.Length();
  const std::size_t t = code.Availability();
  for (std::size_t coord = 0; coord < n; ++coord) {
    const auto &sets = code.RepairSets(coord);
    ASSERT_EQ(sets.size(), t);
    std::set<std::size_t> seen;
    for (std::size_t u = 0; u < t; ++u) {
      EXPECT_EQ(sets[u], code.RepairSet(coord, u));
      EXPECT_EQ(sets[u].size(), code.BlockSize(u) - 1);
      for (std::size_t idx : sets[u]) {
        EXPECT_NE(idx, coord) << "Repair set " << u << " for coord " << coord
                              << " must not include coord itself";
        EXPECT_LT(idx, n);
        EXPECT_TRUE(seen.insert(idx).second)
            << "Coordinate " << idx << " appeared in multiple repair sets for "
            << coord;
      }
    }
  }
}

template <typename GaloisField>
void VerifyAllCoordinatesLocalRepair(
    const LocallyRecoverableCode<GaloisField> &code, uint64_t seed) {
  std::mt19937_64 rng(seed);
  auto msg =
      GFArray<GaloisField>::Random(code.Field(), {code.Dimension()}, rng);
  auto cw = code.Encode(msg);
  ASSERT_TRUE(code.Contains(cw));

  const std::size_t n = code.Length();
  const std::size_t t = code.Availability();
  for (std::size_t coord = 0; coord < n; ++coord) {
    auto received = cw;
    // Corrupt the erased symbol to a different value to prove LocalRepair never
    // reads the erased slot itself.
    received.Raw(coord) =
        code.Field()->Add(cw.Raw(coord), RandomNonZero(code.Field(), rng));

    for (std::size_t u = 0; u < t; ++u) {
      auto res = code.LocalRepair(received, coord, u);
      ASSERT_TRUE(res.Ok()) << res.error_message;
      EXPECT_EQ(res.erased_index, coord);
      EXPECT_EQ(res.repair_set_index, u);
      EXPECT_EQ(res.symbols_read, code.Locality(u));
      EXPECT_EQ(res.recovered_symbol.RawValue(), cw.Raw(coord));
    }
  }
}

// ============================================================================
// 1. Single-Set Tamo-Barg LRCs (Multiplicative & Additive, r | k and r !| k)
// ============================================================================

TEST(LRCTamoBargTest, MultiplicativeDivisibleAndNonDivisibleOverGF13) {
  auto gf13 = std::make_shared<GaloisFieldPrime<uint16_t>>(13);

  // (12, 6, r=3, delta=2): r | k, exact d = 12 - 6 - ceil(6/3) + 2 = 6.
  auto lrc_div =
      LocallyRecoverableCode<GaloisFieldPrime<uint16_t>>::CreateTamoBarg(
          gf13, 12, 6, 3);
  EXPECT_EQ(lrc_div.Length(), 12u);
  EXPECT_EQ(lrc_div.Dimension(), 6u);
  EXPECT_EQ(lrc_div.Locality(), 3u);
  EXPECT_EQ(lrc_div.LocalDistance(), 2u);
  EXPECT_EQ(lrc_div.Availability(), 1u);
  EXPECT_EQ(lrc_div.PartitionType(), LRCPartitionType::kMultiplicative);
  EXPECT_EQ(lrc_div.MinimumDistance(), 6u);
  EXPECT_EQ(lrc_div.GoodPolynomial().Degree(), 4);
  VerifyLinearAlgebraInvariants(lrc_div, 101);
  VerifyStrictPairwiseDisjointRepairSets(lrc_div);
  VerifyAllCoordinatesLocalRepair(lrc_div, 102);

  // (12, 5, r=3, delta=2): r !| k, exact d = 12 - 5 - ceil(5/3) + 2 = 7.
  auto lrc_nondiv =
      LocallyRecoverableCode<GaloisFieldPrime<uint16_t>>::CreateTamoBarg(
          gf13, 12, 5, 3);
  EXPECT_EQ(lrc_nondiv.Length(), 12u);
  EXPECT_EQ(lrc_nondiv.Dimension(), 5u);
  EXPECT_EQ(lrc_nondiv.Locality(), 3u);
  EXPECT_EQ(lrc_nondiv.LocalDistance(), 2u);
  EXPECT_EQ(lrc_nondiv.PartitionType(), LRCPartitionType::kMultiplicative);
  EXPECT_EQ(lrc_nondiv.MinimumDistance(), 7u);
  VerifyLinearAlgebraInvariants(lrc_nondiv, 103);
  VerifyAllCoordinatesLocalRepair(lrc_nondiv, 104);
}

TEST(LRCTamoBargTest, AdditiveOverGF16IncludesZeroAndSupportsFullFieldLength) {
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);

  // (16, 7, r=3, delta=2): s = 4 = 2^2 (auto-selects kAdditive), r !| k,
  // exact d = 16 - 7 - ceil(7/3) + 2 = 8.
  auto lrc = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf16, 16, 7, 3);
  EXPECT_EQ(lrc.Length(), 16u);
  EXPECT_EQ(lrc.Dimension(), 7u);
  EXPECT_EQ(lrc.Locality(), 3u);
  EXPECT_EQ(lrc.LocalDistance(), 2u);
  EXPECT_EQ(lrc.PartitionType(), LRCPartitionType::kAdditive);
  EXPECT_EQ(lrc.MinimumDistance(), 8u);
  EXPECT_EQ(lrc.GoodPolynomial().Degree(), 4);

  // Verify 0 is contained in fiber A_0.
  const auto evals = lrc.EvaluationPoints().RawSpan();
  const auto &fiber0 = lrc.FiberPartition(0)[0];
  bool found_zero_in_a0 = false;
  for (std::size_t idx : fiber0) {
    if (gf16->IsZero(evals[idx])) {
      found_zero_in_a0 = true;
    }
  }
  EXPECT_TRUE(found_zero_in_a0);

  // Verify good polynomial is constant on every fiber and distinct across
  // fibers.
  const auto &g = lrc.GoodPolynomial();
  std::set<uint8_t> fiber_vals;
  for (const auto &fiber : lrc.FiberPartition(0)) {
    const uint8_t c = g.EvalRaw(evals[fiber[0]]);
    for (std::size_t idx : fiber) {
      EXPECT_EQ(g.EvalRaw(evals[idx]), c);
    }
    EXPECT_TRUE(fiber_vals.insert(c).second);
  }

  VerifyLinearAlgebraInvariants(lrc, 201);
  VerifyStrictPairwiseDisjointRepairSets(lrc);
  VerifyAllCoordinatesLocalRepair(lrc, 202);

  // Also check r | k additive: (16, 6, 3) -> d = 16 - 6 - 2 + 2 = 10.
  auto lrc_div = LocallyRecoverableCode<GaloisFieldBinaryExtension<uint8_t>>::
      CreateTamoBarg(gf16, 16, 6, 3, 2, LRCPartitionType::kAdditive);
  EXPECT_EQ(lrc_div.MinimumDistance(), 10u);
  VerifyLinearAlgebraInvariants(lrc_div, 203);
  VerifyAllCoordinatesLocalRepair(lrc_div, 204);
}

TEST(LRCTamoBargTest, AdditiveWorksWithZechLogAndOddPrimeExtensionFields) {
  // GF(2^4) with Zech log internal representation.
  auto gf16_zech = std::make_shared<GFBEZechLogTables<uint8_t>>(4);
  auto lrc_zech =
      LocallyRecoverableCode<GFBEZechLogTables<uint8_t>>::CreateTamoBarg(
          gf16_zech, 16, 7, 3);
  EXPECT_EQ(lrc_zech.PartitionType(), LRCPartitionType::kAdditive);
  EXPECT_EQ(lrc_zech.MinimumDistance(), 8u);
  VerifyLinearAlgebraInvariants(lrc_zech, 210);
  VerifyAllCoordinatesLocalRepair(lrc_zech, 211);

  // GF(3^3) = GF(27) with GaloisFieldExtension (odd characteristic p = 3).
  auto gf27 =
      std::make_shared<GaloisFieldExtension<uint16_t>>(std::make_pair(3, 3));
  // s = 3 = 3^1 -> additive cosets of size 3, n = 15 (5 blocks), k = 7 (r=2,
  // r !| k): d = 15 - 7 - ceil(7/2) + 2 = 6.
  auto lrc_gf27 =
      LocallyRecoverableCode<GaloisFieldExtension<uint16_t>>::CreateTamoBarg(
          gf27, 15, 7, 2);
  EXPECT_EQ(lrc_gf27.PartitionType(), LRCPartitionType::kAdditive);
  EXPECT_EQ(lrc_gf27.MinimumDistance(), 6u);
  VerifyLinearAlgebraInvariants(lrc_gf27, 212);
  VerifyAllCoordinatesLocalRepair(lrc_gf27, 213);
}

TEST(LRCTamoBargTest, ExhaustiveExactMinimumDistanceVerification) {
  // 1. Multiplicative r | k over GF(7): (6, 4, r=2, delta=2) -> exact d = 2.
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto m_div =
      LocallyRecoverableCode<GaloisFieldPrime<uint8_t>>::CreateTamoBarg(gf7, 6,
                                                                        4, 2);
  EXPECT_EQ(m_div.MinimumDistance(), 2u);
  EXPECT_EQ(ExhaustiveMinimumDistance(m_div), 2u);

  // 2. Multiplicative r !| k over GF(7): (6, 3, r=2, delta=2) -> exact d = 3.
  auto m_nondiv =
      LocallyRecoverableCode<GaloisFieldPrime<uint8_t>>::CreateTamoBarg(gf7, 6,
                                                                        3, 2);
  EXPECT_EQ(m_nondiv.MinimumDistance(), 3u);
  EXPECT_EQ(ExhaustiveMinimumDistance(m_nondiv), 3u);

  // 3. Multiplicative r !| k over GF(13): (12, 3, r=2, delta=2) -> exact d = 9.
  auto gf13 = std::make_shared<GaloisFieldPrime<uint16_t>>(13);
  auto m13 = LocallyRecoverableCode<GaloisFieldPrime<uint16_t>>::CreateTamoBarg(
      gf13, 12, 3, 2);
  EXPECT_EQ(m13.MinimumDistance(), 9u);
  EXPECT_EQ(ExhaustiveMinimumDistance(m13), 9u);

  // 4. Additive r !| k over GF(8): (8, 4, r=3, delta=2) -> exact d = 4.
  auto gf8 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);
  auto a8_nondiv = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf8, 8, 4, 3);
  EXPECT_EQ(a8_nondiv.PartitionType(), LRCPartitionType::kAdditive);
  EXPECT_EQ(a8_nondiv.MinimumDistance(), 4u);
  EXPECT_EQ(ExhaustiveMinimumDistance(a8_nondiv), 4u);

  // 5. Additive r | k over GF(8): (8, 3, r=3, delta=2) -> exact d = 6 (MDS
  //    regime k = r).
  auto a8_eq = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf8, 8, 3, 3);
  EXPECT_EQ(a8_eq.MinimumDistance(), 6u);
  EXPECT_EQ(ExhaustiveMinimumDistance(a8_eq), 6u);

  // 6. Additive over GF(9) (p = 3): (9, 3, r=2, delta=2) -> exact d = 6,
  //    and (9, 4, r=2, delta=2) -> exact d = 5.
  auto gf9 =
      std::make_shared<GaloisFieldExtension<uint8_t>>(std::make_pair(3, 2));
  auto a9_nondiv =
      LocallyRecoverableCode<GaloisFieldExtension<uint8_t>>::CreateTamoBarg(
          gf9, 9, 3, 2);
  EXPECT_EQ(a9_nondiv.PartitionType(), LRCPartitionType::kAdditive);
  EXPECT_EQ(a9_nondiv.MinimumDistance(), 6u);
  EXPECT_EQ(ExhaustiveMinimumDistance(a9_nondiv), 6u);

  auto a9_div =
      LocallyRecoverableCode<GaloisFieldExtension<uint8_t>>::CreateTamoBarg(
          gf9, 9, 4, 2);
  EXPECT_EQ(a9_div.MinimumDistance(), 5u);
  EXPECT_EQ(ExhaustiveMinimumDistance(a9_div), 5u);

  // 7. Generalized (r, delta)-locality over GF(16): (15, 3, r=2, delta=4)
  //    s = 5, ceil(3/2) = 2 -> exact d = 15 - 3 + 1 - 1*3 = 10.
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto r_delta_small = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf16, 15, 3, 2, 4);
  EXPECT_EQ(r_delta_small.MinimumDistance(), 10u);
  EXPECT_EQ(ExhaustiveMinimumDistance(r_delta_small), 10u);
}

// ============================================================================
// 2. Generalized (r, delta)-Locality
// ============================================================================

TEST(LRCRDeltaTest, LocalRepairTwoSimultaneousErasuresInOneFiber) {
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);

  // (15, 6, r=3, delta=3) over GF(16): s = 5 | 15, m_B = 3,
  // exact d = 15 - 6 + 1 - (2 - 1) * (3 - 1) = 8.
  auto lrc = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf16, 15, 6, 3, 3);
  EXPECT_EQ(lrc.Length(), 15u);
  EXPECT_EQ(lrc.Dimension(), 6u);
  EXPECT_EQ(lrc.Locality(), 3u);
  EXPECT_EQ(lrc.LocalDistance(), 3u);
  EXPECT_EQ(lrc.BlockSize(), 5u);
  EXPECT_EQ(lrc.MinimumDistance(), 8u);
  EXPECT_EQ(lrc.PartitionType(), LRCPartitionType::kMultiplicative);
  VerifyLinearAlgebraInvariants(lrc, 301);

  std::mt19937_64 rng(302);
  auto msg =
      GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf16, {6}, rng);
  auto cw = lrc.Encode(msg);

  // Test EVERY pair of 2 simultaneous erasures inside fiber 0 (size 5).
  const auto &fiber0 = lrc.FiberPartition(0)[0];
  ASSERT_EQ(fiber0.size(), 5u);
  for (std::size_t a = 0; a < fiber0.size(); ++a) {
    for (std::size_t b = a + 1; b < fiber0.size(); ++b) {
      const std::size_t pos_a = fiber0[a];
      const std::size_t pos_b = fiber0[b];
      auto received = ErasePositions(cw, {pos_a, pos_b});

      // LocalRepair on pos_a with pos_b also erased must read strictly r = 3.
      auto res_a = lrc.LocalRepair(received, pos_a, 0, {pos_b});
      ASSERT_TRUE(res_a.Ok()) << res_a.error_message;
      EXPECT_EQ(res_a.symbols_read, 3u);
      EXPECT_EQ(res_a.recovered_symbol.RawValue(), cw.Raw(pos_a));

      // LocalRepair on pos_b with pos_a also erased must read strictly r = 3.
      auto res_b = lrc.LocalRepair(received, pos_b, 0, {pos_a});
      ASSERT_TRUE(res_b.Ok()) << res_b.error_message;
      EXPECT_EQ(res_b.symbols_read, 3u);
      EXPECT_EQ(res_b.recovered_symbol.RawValue(), cw.Raw(pos_b));

      // LocalRepairMultiple must recover BOTH erasures in one single local
      // interpolation reading only r = 3 symbols total!
      auto multi = lrc.LocalRepairMultiple(received, {pos_a, pos_b},
                                           /*allow_global_fallback=*/false);
      ASSERT_TRUE(multi.Ok()) << multi.error_message;
      EXPECT_EQ(multi.codeword, cw);
      EXPECT_EQ(multi.message, msg);
      EXPECT_EQ(multi.num_local_repairs, 2u);
      EXPECT_EQ(multi.num_global_repairs, 0u);
      EXPECT_EQ(multi.total_symbols_read, 3u);
      EXPECT_FALSE(multi.used_global_fallback);
    }
  }

  // 3 erasures (= delta) in one fiber must fail local repair.
  auto received3 = ErasePositions(cw, {fiber0[0], fiber0[1], fiber0[2]});
  auto fail_local =
      lrc.LocalRepair(received3, fiber0[0], 0, {fiber0[1], fiber0[2]});
  EXPECT_FALSE(fail_local.Ok());
}

// ============================================================================
// 3. Custom Good Polynomial Construction & Validation Errors
// ============================================================================

TEST(LRCCustomGoodPolyTest, ValidConstructionAndValidationRejections) {
  using Field = GaloisFieldPrime<uint16_t>;
  using Elem = GaloisFieldElement<Field>;
  using Poly = PolynomialDense<Field>;

  auto gf13 = std::make_shared<Field>(13);
  // Good polynomial g(x) = x^3 + 1 over GF(13).
  // Cubic residues mod 13:
  //   1^3 = 3^3 = 9^3 = 1  -> g = 2
  //   2^3 = 6^3 = 5^3 = 8  -> g = 9
  //   4^3 = 12^3 = 10^3 = 12 -> g = 0
  //   7^3 = 8^3 = 11^3 = 5 -> g = 6
  Poly g = Poly::FromRaw({1, 0, 0, 1}, gf13);
  std::vector<std::vector<Elem>> fibers = {
      {Elem(1, gf13), Elem(3, gf13), Elem(9, gf13)},
      {Elem(2, gf13), Elem(5, gf13), Elem(6, gf13)},
      {Elem(4, gf13), Elem(10, gf13), Elem(12, gf13)},
      {Elem(7, gf13), Elem(8, gf13), Elem(11, gf13)},
  };

  // n = 12, k = 5, r = 2, delta = 2 -> exact d = 12 - 5 - ceil(5/2) + 2 = 6.
  auto lrc = LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
      gf13, g, fibers, 5, 2);
  EXPECT_EQ(lrc.PartitionType(), LRCPartitionType::kCustom);
  EXPECT_EQ(lrc.Length(), 12u);
  EXPECT_EQ(lrc.Dimension(), 5u);
  EXPECT_EQ(lrc.Locality(), 2u);
  EXPECT_EQ(lrc.LocalDistance(), 2u);
  EXPECT_EQ(lrc.MinimumDistance(), 6u);
  EXPECT_EQ(lrc.GoodPolynomial(), g);
  VerifyLinearAlgebraInvariants(lrc, 401);
  VerifyAllCoordinatesLocalRepair(lrc, 402);

  // Also test the GFArray fiber overload.
  std::vector<GFArray<Field>> arr_fibers;
  for (const auto &f : fibers) {
    typename GFArray<Field>::raw_array_type raw =
        GFArray<Field>::raw_array_type::from_shape({f.size()});
    for (std::size_t i = 0; i < f.size(); ++i) {
      raw(i) = f[i].RawValue();
    }
    arr_fibers.push_back(GFArray<Field>::FromRaw(gf13, std::move(raw)));
  }
  auto lrc_arr = LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
      gf13, g, arr_fibers, 5, 2);
  EXPECT_EQ(lrc_arr.GeneratorMatrix(), lrc.GeneratorMatrix());

  // Validation error checks:
  // (a) Wrong fiber size.
  auto bad_size = fibers;
  bad_size[1].pop_back();
  EXPECT_THROW(LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
                   gf13, g, bad_size, 5, 2),
               std::invalid_argument);

  // (b) Duplicate evaluation point across fibers.
  auto dup_point = fibers;
  dup_point[1] = fibers[0];
  EXPECT_THROW(LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
                   gf13, g, dup_point, 5, 2),
               std::invalid_argument);

  // (c) Non-constant g on a fiber.
  auto non_const = fibers;
  non_const[0][2] = Elem(2, gf13);  // g(2) = 9 != g(1) = 2.
  EXPECT_THROW(LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
                   gf13, g, non_const, 5, 2),
               std::invalid_argument);

  // (d) Duplicate fiber constant g(A_0) == g(A_1) (using two fibers with
  //     distinct points that share the same g value — build via g(x) = x^2
  //     where g(1) = g(12) = 1, plus duplicate constant test using constant
  //     g(x) = x^3 + 1 with fibers of size 3: use a deg-2 poly where 2 fibers
  //     would collide only if points repeat; instead use g(x) = x^4 - x^2 on
  //     GF(13) or verify duplicate fiber constant with g(x) = x^6 + 1 where
  //     6 points map to 2: {1,3,9} and {4,10,12} don't have size 6. With
  //     g(x) = x^2: 1^2 = 12^2 = 1, 2^2 = 11^2 = 4 -> distinct. To get two
  //     disjoint fibers of size 2 with the same g-value, use deg(g) = 2 is
  //     at most 2 roots per value, so deg(g) = s means at most s points can
  //     share any value c_ell; thus distinct points with equal c_ell can only
  //     happen if fiber size != s or points duplicate. Still verify degree < 2
  //     and invalid locality/dimension!).
  Poly deg1 = Poly::FromRaw({1, 1}, gf13);
  EXPECT_THROW(LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
                   gf13, deg1, fibers, 2, 1),
               std::invalid_argument);
  EXPECT_THROW(LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
                   gf13, g, fibers, 5, 3),  // locality >= deg(g)
               std::invalid_argument);
  EXPECT_THROW(LocallyRecoverableCode<Field>::CreateFromGoodPolynomial(
                   gf13, g, fibers, 9, 2),  // k > m_B * r = 8
               std::invalid_argument);
}

// ============================================================================
// 4. Multi-Set Availability LRCs
// ============================================================================

TEST(LRCAvailabilityTest,
     MultiplicativeCRTTwoSetsOverGF13AndThreeSetsOverGF31) {
  // t = 2 over GF(13): s_1 = 3 (r_1=2, d_1=2), s_2 = 4 (r_2=3, d_2=2), S = 12.
  auto gf13 = std::make_shared<GaloisFieldPrime<uint16_t>>(13);
  auto lrc13 = LocallyRecoverableCode<
      GaloisFieldPrime<uint16_t>>::CreateMultiplicativeAvailability(gf13, 12, 4,
                                                                    {2, 3});
  EXPECT_EQ(lrc13.PartitionType(),
            LRCPartitionType::kAvailabilityMultiplicative);
  EXPECT_EQ(lrc13.Availability(), 2u);
  EXPECT_EQ(lrc13.Locality(0), 2u);
  EXPECT_EQ(lrc13.Locality(1), 3u);
  EXPECT_EQ(lrc13.LocalDistance(0), 2u);
  EXPECT_EQ(lrc13.LocalDistance(1), 2u);
  // CRT exponents E = {0, 1, 4, 6, 9, 10}; smallest 4 -> e_3 = 6 -> d >= 6.
  EXPECT_EQ(lrc13.MinimumDistance(), 6u);
  VerifyLinearAlgebraInvariants(lrc13, 501);
  VerifyStrictPairwiseDisjointRepairSets(lrc13);
  VerifyAllCoordinatesLocalRepair(lrc13, 502);

  // Verify repair via set 1 still succeeds even when ALL of set 0 is erased!
  {
    std::mt19937_64 rng(503);
    auto msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf13, {4}, rng);
    auto cw = lrc13.Encode(msg);
    const std::size_t target = 0;
    auto other_erased = lrc13.RepairSet(target, 0);
    auto all_erased = other_erased;
    all_erased.push_back(target);
    auto received = ErasePositions(cw, all_erased);

    EXPECT_FALSE(lrc13.LocalRepair(received, target, 0, other_erased).Ok());
    auto res1 = lrc13.LocalRepair(received, target, 1, other_erased);
    ASSERT_TRUE(res1.Ok()) << res1.error_message;
    EXPECT_EQ(res1.symbols_read, 3u);
    EXPECT_EQ(res1.recovered_symbol.RawValue(), cw.Raw(target));
  }

  // t = 3 over GF(31): s_1 = 2 (r_1=1), s_2 = 3 (r_2=2), s_3 = 5 (r_3=4),
  // S = 30 | 30.
  auto gf31 = std::make_shared<GaloisFieldPrime<uint16_t>>(31);
  auto lrc31 = LocallyRecoverableCode<
      GaloisFieldPrime<uint16_t>>::CreateMultiplicativeAvailability(gf31, 30, 6,
                                                                    {1, 2, 4});
  EXPECT_EQ(lrc31.Availability(), 3u);
  EXPECT_EQ(lrc31.Locality(0), 1u);
  EXPECT_EQ(lrc31.Locality(1), 2u);
  EXPECT_EQ(lrc31.Locality(2), 4u);
  // E = {0, 6, 10, 12, 16, 18, 22, 28}; smallest 6 -> e_5 = 18 -> d >= 12.
  EXPECT_EQ(lrc31.MinimumDistance(), 12u);
  VerifyLinearAlgebraInvariants(lrc31, 504);
  VerifyStrictPairwiseDisjointRepairSets(lrc31);
  VerifyAllCoordinatesLocalRepair(lrc31, 505);

  // Also test t = 3 with (r_3, delta_3) = (3, 3) (s_3 = 5) over GF(31).
  auto lrc31_rd = LocallyRecoverableCode<
      GaloisFieldPrime<uint16_t>>::CreateMultiplicativeAvailability(gf31, 30, 5,
                                                                    {1, 2, 3},
                                                                    {2, 2, 3});
  EXPECT_EQ(lrc31_rd.LocalDistance(2), 3u);
  VerifyStrictPairwiseDisjointRepairSets(lrc31_rd);
  VerifyAllCoordinatesLocalRepair(lrc31_rd, 506);

  // Non-coprime block sizes must throw.
  EXPECT_THROW(
      LocallyRecoverableCode<GaloisFieldPrime<uint16_t>>::
          CreateMultiplicativeAvailability(gf13, 12, 2, {1, 3}),  // s = 2, 4
      std::invalid_argument);
}

TEST(LRCAvailabilityTest, AdditiveDirectSumOverGF16GF64AndGF9) {
  // 1. H_1 (+) H_2 = GF(16): s_1 = s_2 = 4, m_B = 1, m_blk = 1.
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto lrc16 = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateAdditiveAvailability(gf16, 1,
                                                                       1,
                                                                       {3, 3});
  EXPECT_EQ(lrc16.PartitionType(), LRCPartitionType::kAvailabilityAdditive);
  EXPECT_EQ(lrc16.Availability(), 2u);
  EXPECT_EQ(lrc16.Length(), 16u);
  EXPECT_EQ(lrc16.Dimension(), 9u);
  EXPECT_EQ(lrc16.MinimumDistance(), 4u);  // (1 - 1 + 1) * 2 * 2 = 4.
  VerifyLinearAlgebraInvariants(lrc16, 601);
  VerifyStrictPairwiseDisjointRepairSets(lrc16);
  VerifyAllCoordinatesLocalRepair(lrc16, 602);

  // Also test higher local distance delta_1 = delta_2 = 3 (r_1 = r_2 = 2):
  // exact d = 1 * 3 * 3 = 9.
  auto lrc16_rd = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateAdditiveAvailability(gf16, 1,
                                                                       1,
                                                                       {2, 2},
                                                                       {3, 3});
  EXPECT_EQ(lrc16_rd.Dimension(), 4u);
  EXPECT_EQ(lrc16_rd.MinimumDistance(), 9u);
  VerifyLinearAlgebraInvariants(lrc16_rd, 603);
  VerifyStrictPairwiseDisjointRepairSets(lrc16_rd);
  VerifyAllCoordinatesLocalRepair(lrc16_rd, 604);

  // 2. Multi-superblock over GF(64): s_1 = s_2 = 4, m_B = 4, m_blk = 2:
  //    n = 64, k = 2 * 3 * 3 = 18, exact d = (4 - 2 + 1) * 2 * 2 = 12.
  auto gf64 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(6);
  auto lrc64 = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateAdditiveAvailability(gf64, 4,
                                                                       2,
                                                                       {3, 3});
  EXPECT_EQ(lrc64.Length(), 64u);
  EXPECT_EQ(lrc64.Dimension(), 18u);
  EXPECT_EQ(lrc64.MinimumDistance(), 12u);
  VerifyLinearAlgebraInvariants(lrc64, 605);
  VerifyStrictPairwiseDisjointRepairSets(lrc64);
  VerifyAllCoordinatesLocalRepair(lrc64, 606);

  // 3. Exhaustive exact minimum distance check over GF(9) (p = 3):
  //    s_1 = s_2 = 3, r_1 = r_2 = 2, m_B = 1, m_blk = 1 -> n = 9, k = 4, d = 4.
  auto gf9 =
      std::make_shared<GaloisFieldExtension<uint8_t>>(std::make_pair(3, 2));
  auto lrc9 = LocallyRecoverableCode<
      GaloisFieldExtension<uint8_t>>::CreateAdditiveAvailability(gf9, 1, 1,
                                                                 {2, 2});
  EXPECT_EQ(lrc9.Length(), 9u);
  EXPECT_EQ(lrc9.Dimension(), 4u);
  EXPECT_EQ(lrc9.MinimumDistance(), 4u);
  EXPECT_EQ(ExhaustiveMinimumDistance(lrc9), 4u);
  VerifyStrictPairwiseDisjointRepairSets(lrc9);
  VerifyAllCoordinatesLocalRepair(lrc9, 607);
}

TEST(LRCAvailabilityTest, HybridSemidirectProductOverGF16AndGF64) {
  // 1. H_1 |><| H_2 on GF(16) \ GF(4): b = 2, a = 1 (s_1 = 4), s_2 = 3,
  //    num_orbits = 1 -> n = 12, r_1 = 3, r_2 = 2, D = 8.
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto hybrid16 = LocallyRecoverableCode<GaloisFieldBinaryExtension<uint8_t>>::
      CreateHybridAvailability(gf16, /*num_orbits=*/1, /*subfield_degree_b=*/2,
                               /*additive_dim_a=*/1, /*mult_order_s2=*/3,
                               /*r1=*/3, /*r2=*/2, /*max_poly_degree=*/8);
  EXPECT_EQ(hybrid16.PartitionType(), LRCPartitionType::kAvailabilityHybrid);
  EXPECT_EQ(hybrid16.Availability(), 2u);
  EXPECT_EQ(hybrid16.Length(), 12u);
  EXPECT_EQ(hybrid16.Locality(0), 3u);
  EXPECT_EQ(hybrid16.Locality(1), 2u);
  EXPECT_EQ(hybrid16.LocalDistance(0), 2u);
  EXPECT_EQ(hybrid16.LocalDistance(1), 2u);
  EXPECT_GE(hybrid16.Dimension(), 4u);
  EXPECT_GE(hybrid16.MinimumDistance(), 5u);
  VerifyLinearAlgebraInvariants(hybrid16, 701);
  VerifyStrictPairwiseDisjointRepairSets(hybrid16);
  VerifyAllCoordinatesLocalRepair(hybrid16, 702);

  // Verify global decoder works on hybrid availability LRC up to 2e + s <= d-1.
  {
    std::mt19937_64 rng(703);
    auto msg = GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(
        gf16, {hybrid16.Dimension()}, rng);
    auto cw = hybrid16.Encode(msg);
    const std::size_t d = hybrid16.MinimumDistance();
    const std::size_t t_err = (d - 1) / 2;
    if (t_err >= 1) {
      std::vector<std::size_t> err_pos(t_err);
      std::iota(err_pos.begin(), err_pos.end(), std::size_t{0});
      auto noisy = CorruptPositions(cw, err_pos, rng);
      auto dec = hybrid16.DecodeWithResult(noisy);
      ASSERT_TRUE(dec.Ok()) << dec.error_message;
      EXPECT_EQ(dec.codeword, cw);
      EXPECT_EQ(dec.message, msg);
      EXPECT_EQ(dec.num_errors, t_err);
    }
  }

  // 2. Multi-orbit hybrid over GF(64): b = 2, a = 1 (s_1 = 4), s_2 = 3,
  //    num_orbits = 2 -> n = 24.
  auto gf64 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(6);
  auto hybrid64 = LocallyRecoverableCode<GaloisFieldBinaryExtension<uint8_t>>::
      CreateHybridAvailability(gf64, /*num_orbits=*/2, /*subfield_degree_b=*/2,
                               /*additive_dim_a=*/1, /*mult_order_s2=*/3,
                               /*r1=*/3, /*r2=*/2, /*max_poly_degree=*/14);
  EXPECT_EQ(hybrid64.Length(), 24u);
  EXPECT_EQ(hybrid64.Availability(), 2u);
  EXPECT_GE(hybrid64.Dimension(), 6u);
  VerifyLinearAlgebraInvariants(hybrid64, 704);
  VerifyStrictPairwiseDisjointRepairSets(hybrid64);
  VerifyAllCoordinatesLocalRepair(hybrid64, 705);
}

// ============================================================================
// 5. Cascade LocalRepairMultiple + Global Fallback + Global/Hybrid Decoders
// ============================================================================

TEST(LRCRepairAndDecodingTest, PeelingCascadeAndForcedGlobalFallback) {
  auto gf13 = std::make_shared<GaloisFieldPrime<uint16_t>>(13);

  // (12, 6, r=3, delta=2), 3 fibers of size 4, d = 6.
  auto lrc = LocallyRecoverableCode<GaloisFieldPrime<uint16_t>>::CreateTamoBarg(
      gf13, 12, 6, 3);
  std::mt19937_64 rng(801);
  auto msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf13, {6}, rng);
  auto cw = lrc.Encode(msg);

  // (a) Pure local cascade: 1 erasure in each of the 3 fibers.
  const auto &fibers = lrc.FiberPartition(0);
  const std::vector<std::size_t> one_per_fiber = {fibers[0][1], fibers[1][2],
                                                  fibers[2][0]};
  auto rec_local = ErasePositions(cw, one_per_fiber);
  auto res_local = lrc.LocalRepairMultiple(rec_local, one_per_fiber,
                                           /*allow_global_fallback=*/false);
  ASSERT_TRUE(res_local.Ok()) << res_local.error_message;
  EXPECT_EQ(res_local.codeword, cw);
  EXPECT_EQ(res_local.message, msg);
  EXPECT_EQ(res_local.num_local_repairs, 3u);
  EXPECT_EQ(res_local.num_global_repairs, 0u);
  EXPECT_EQ(res_local.total_symbols_read, 9u);  // 3 fibers * r=3.
  EXPECT_FALSE(res_local.used_global_fallback);

  // (b) Forced global fallback: 2 erasures (= delta) in fiber 0 + 1 in fiber 1
  //     + 2 erasures in fiber 2 -> total 5 = d - 1 erasures!
  //     Fiber 1 repairs locally (1 repair), fibers 0 & 2 require global solve
  //     (4 global repairs).
  const std::vector<std::size_t> heavy_erasures = {
      fibers[0][0], fibers[0][1], fibers[1][2], fibers[2][0], fibers[2][3]};
  auto rec_heavy = ErasePositions(cw, heavy_erasures);

  auto res_no_fb = lrc.LocalRepairMultiple(rec_heavy, heavy_erasures,
                                           /*allow_global_fallback=*/false);
  EXPECT_FALSE(res_no_fb.Ok());
  EXPECT_EQ(res_no_fb.num_local_repairs, 1u);

  auto res_fb = lrc.LocalRepairMultiple(rec_heavy, heavy_erasures,
                                        /*allow_global_fallback=*/true);
  ASSERT_TRUE(res_fb.Ok()) << res_fb.error_message;
  EXPECT_EQ(res_fb.codeword, cw);
  EXPECT_EQ(res_fb.message, msg);
  EXPECT_EQ(res_fb.num_local_repairs, 1u);
  EXPECT_EQ(res_fb.num_global_repairs, 4u);
  EXPECT_TRUE(res_fb.used_global_fallback);

  // (c) Multi-round peeling cascade on availability LRC (GF(16) H_1 (+) H_2):
  //     Erase 2 coordinates in the same H_1-fiber (blocks H_1 initially!) that
  //     lie in DISTINCT H_2-fibers -> H_2 peels them locally without fallback!
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto avail16 = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateAdditiveAvailability(gf16, 1,
                                                                       1,
                                                                       {3, 3});
  auto msg16 = GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(
      gf16, {avail16.Dimension()}, rng);
  auto cw16 = avail16.Encode(msg16);

  const auto &h1_fiber0 = avail16.FiberPartition(0)[0];
  const auto &h1_fiber1 = avail16.FiberPartition(0)[1];
  // Erase 2 symbols in h1_fiber0 (cannot repair via partition 0), plus 1 symbol
  // in h1_fiber1 sharing an H_2-fiber with h1_fiber0[0]:
  //   - Round 1: partition 0 repairs h1_fiber1[0] (only 1 erasure in h1_fiber1)
  //              AND partition 1 repairs h1_fiber0[1] (only 1 erasure in its
  //              H_2-fiber)!
  //   - Round 2: h1_fiber0[0] now has BOTH partition 0 and partition 1
  //   unblocked!
  const std::vector<std::size_t> cascade_erasures = {h1_fiber0[0], h1_fiber0[1],
                                                     h1_fiber1[0]};
  auto rec_cascade = ErasePositions(cw16, cascade_erasures);
  auto res_cascade = avail16.LocalRepairMultiple(
      rec_cascade, cascade_erasures, /*allow_global_fallback=*/false);
  ASSERT_TRUE(res_cascade.Ok()) << res_cascade.error_message;
  EXPECT_EQ(res_cascade.codeword, cw16);
  EXPECT_EQ(res_cascade.message, msg16);
  EXPECT_EQ(res_cascade.num_local_repairs, 3u);
  EXPECT_EQ(res_cascade.num_global_repairs, 0u);
  EXPECT_FALSE(res_cascade.used_global_fallback);
}

TEST(LRCRepairAndDecodingTest,
     GlobalDecoderCorrectsErrorsAndErasuresAndRejectsAmbientRSWords) {
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);

  // Additive (16, 7, r=3, delta=2) over GF(16): d = 8, includes 0 in A_0!
  // Max polynomial degree D_max = n - d = 8, ambient RS is [16, 9, 8].
  auto lrc = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf16, 16, 7, 3);
  ASSERT_EQ(lrc.MinimumDistance(), 8u);
  ASSERT_EQ(lrc.DefaultDecoder(), "LRCGlobal");

  std::mt19937_64 rng(901);
  auto msg =
      GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf16, {7}, rng);
  auto cw = lrc.Encode(msg);

  // Test all valid (e, s) combinations with 2e + s <= d - 1 = 7, including
  // index 0 (which is the zero evaluation point 0 in A_0!).
  for (std::size_t num_err = 0; num_err <= 3; ++num_err) {
    for (std::size_t num_era = 0; 2 * num_err + num_era <= 7; ++num_era) {
      std::vector<std::size_t> eras(num_era);
      std::iota(eras.begin(), eras.end(), std::size_t{0});
      std::vector<std::size_t> errs(num_err);
      std::iota(errs.begin(), errs.end(), num_era);

      auto received = CorruptPositions(cw, errs, rng);
      received = ErasePositions(received, eras);

      // Default LRCGlobal decoder:
      auto res_global = lrc.DecodeWithErasures(received, eras, "LRCGlobal");
      ASSERT_TRUE(res_global.Ok())
          << "LRCGlobal failed for e=" << num_err << ", s=" << num_era << ": "
          << res_global.error_message;
      EXPECT_EQ(res_global.codeword, cw);
      EXPECT_EQ(res_global.message, msg);

      // LRCLocalThenGlobal decoder:
      auto res_hybrid =
          lrc.DecodeWithErasures(received, eras, "LRCLocalThenGlobal");
      ASSERT_TRUE(res_hybrid.Ok())
          << "LRCLocalThenGlobal failed for e=" << num_err << ", s=" << num_era
          << ": " << res_hybrid.error_message;
      EXPECT_EQ(res_hybrid.codeword, cw);
      EXPECT_EQ(res_hybrid.message, msg);
    }
  }

  // Test LRCLocalThenGlobal on pure erasures BEYOND d - 1 = 7 when locally
  // repairable: over (15, 6, r=3, delta=3), d = 8, 2 erasures in each of 3
  // fibers = 6 erasures; for (16, 12, r=3, delta=2), d = 2, 1 erasure in each
  // of 4 fibers = 4 erasures > d - 1 = 1!
  auto high_rate = LocallyRecoverableCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateTamoBarg(gf16, 16, 12, 3);
  ASSERT_EQ(high_rate.MinimumDistance(), 2u);
  auto hr_msg =
      GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf16, {12}, rng);
  auto hr_cw = high_rate.Encode(hr_msg);
  const auto &hr_fibers = high_rate.FiberPartition(0);
  const std::vector<std::size_t> four_local_erasures = {
      hr_fibers[0][0], hr_fibers[1][1], hr_fibers[2][2], hr_fibers[3][3]};
  auto hr_recv = ErasePositions(hr_cw, four_local_erasures);
  // Ambient RS global decoder cannot decode 4 > d - 1 = 1 erasures,
  // but LRCLocalThenGlobal succeeds via local repair!
  EXPECT_FALSE(
      high_rate.DecodeWithErasures(hr_recv, four_local_erasures, "LRCGlobal")
          .Ok());
  auto hr_res = high_rate.DecodeWithErasures(hr_recv, four_local_erasures,
                                             "LRCLocalThenGlobal");
  ASSERT_TRUE(hr_res.Ok()) << hr_res.error_message;
  EXPECT_EQ(hr_res.codeword, hr_cw);
  EXPECT_EQ(hr_res.message, hr_msg);

  // CRITICAL: Reject non-LRC ambient RS codewords!
  // Polynomial f(x) = x^3 has degree 3 <= D_max = 8, so its evaluation lies in
  // the ambient [16, 9, 8] RS code, but NOT in the (16, 7, 3) LRC subcode
  // (whose monomials skip degree 3 mod 4).
  auto ambient_rs =
      GeneralizedReedSolomonCode<GaloisFieldBinaryExtension<uint8_t>>::
          StandardReedSolomon(gf16, lrc.EvaluationPoints(),
                              lrc.MaxPolynomialDegree() + 1);
  auto forbidden_poly =
      PolynomialDense<GaloisFieldBinaryExtension<uint8_t>>::Monomial(
          gf16, 3, gf16->MultiplicativeIdentity());
  auto non_lrc_rs_word = ambient_rs.EncodePolynomial(forbidden_poly);
  ASSERT_TRUE(ambient_rs.Contains(non_lrc_rs_word));
  ASSERT_FALSE(lrc.Contains(non_lrc_rs_word));

  // Even with 0 errors or 1 error on the ambient RS word, LRCGlobal MUST
  // reject it instead of returning a non-LRC word!
  auto reject0 = lrc.DecodeWithResult(non_lrc_rs_word, "LRCGlobal");
  EXPECT_FALSE(reject0.Ok());
  auto reject0_hybrid =
      lrc.DecodeWithResult(non_lrc_rs_word, "LRCLocalThenGlobal");
  EXPECT_FALSE(reject0_hybrid.Ok());

  auto noisy_non_lrc = CorruptPositions(non_lrc_rs_word, {2}, rng);
  auto reject1 = lrc.DecodeWithResult(noisy_non_lrc, "LRCGlobal");
  EXPECT_FALSE(reject1.Ok());
}

TEST(LRCRepairAndDecodingTest, CopyAndMovePreserveRuleOfFiveCacheSafety) {
  auto gf13 = std::make_shared<GaloisFieldPrime<uint16_t>>(13);
  using LRC = LocallyRecoverableCode<GaloisFieldPrime<uint16_t>>;

  std::mt19937_64 rng(1001);
  auto msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf13, {5}, rng);

  std::unique_ptr<LRC> copied;
  std::unique_ptr<LRC> moved;
  GFArray<GaloisFieldPrime<uint16_t>> expected_cw;

  {
    auto original = LRC::CreateTamoBarg(gf13, 12, 5, 3);
    // Warm up encoder and decoder caches on `original`.
    expected_cw = original.Encode(msg);
    auto dec = original.DecodeWithResult(expected_cw);
    ASSERT_TRUE(dec.Ok());
    auto dec2 = original.DecodeWithResult(expected_cw, "LRCLocalThenGlobal");
    ASSERT_TRUE(dec2.Ok());

    // Copy-construct and move-construct before destroying `original`.
    copied = std::make_unique<LRC>(original);
    LRC temp_copy = original;
    moved = std::make_unique<LRC>(std::move(temp_copy));
  }  // `original` and `temp_copy` are now destroyed.

  // Verify copy-constructed and move-constructed LRCs work cleanly.
  for (const LRC *inst : {copied.get(), moved.get()}) {
    auto cw = inst->Encode(msg);
    EXPECT_EQ(cw, expected_cw);
    auto noisy = CorruptPositions(cw, {0, 4, 8}, rng);
    auto res = inst->DecodeWithResult(noisy);
    ASSERT_TRUE(res.Ok()) << res.error_message;
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);

    auto erased = ErasePositions(cw, {3});
    auto loc = inst->LocalRepair(erased, 3);
    ASSERT_TRUE(loc.Ok());
    EXPECT_EQ(loc.recovered_symbol.RawValue(), cw.Raw(3));
  }

  // Also test copy-assignment and move-assignment.
  auto target_copy = LRC::CreateTamoBarg(gf13, 12, 6, 3);
  (void)target_copy.Encode(
      GFArray<GaloisFieldPrime<uint16_t>>::Zeros(gf13, {6}));
  target_copy = *copied;
  EXPECT_EQ(target_copy.Dimension(), 5u);
  EXPECT_EQ(target_copy.Encode(msg), expected_cw);

  auto target_move = LRC::CreateTamoBarg(gf13, 12, 6, 3);
  target_move = std::move(*moved);
  EXPECT_EQ(target_move.Dimension(), 5u);
  EXPECT_EQ(target_move.Encode(msg), expected_cw);
}

}  // namespace
