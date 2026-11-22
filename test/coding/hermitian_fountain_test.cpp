#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/coding/fountain.hpp"
#include "xgalois/coding/hermitian.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/linalg/linalg.hpp"

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
GFArray<GaloisField> CorruptRandomPositions(const GFArray<GaloisField> &cw,
                                            std::size_t num_errors,
                                            std::mt19937_64 &rng) {
  auto field = cw.Field();
  const std::size_t n = cw.Size();
  std::vector<std::size_t> indices(n);
  std::iota(indices.begin(), indices.end(), std::size_t{0});
  std::shuffle(indices.begin(), indices.end(), rng);

  auto received = cw;
  for (std::size_t i = 0; i < num_errors; ++i) {
    const auto err = RandomNonZero(field, rng);
    received.Raw(indices[i]) = field->Add(received.Raw(indices[i]), err);
  }
  return received;
}

template <typename GaloisField>
void VerifyHermitianCurveAndDual(const std::shared_ptr<GaloisField> &field,
                                 std::size_t expected_q,
                                 std::size_t expected_genus, std::size_t m) {
  HermitianCode<GaloisField> code(field, m);
  const std::size_t expected_n = expected_q * expected_q * expected_q;
  EXPECT_EQ(code.Q(), expected_q);
  EXPECT_EQ(code.Genus(), expected_genus);
  EXPECT_EQ(code.Length(), expected_n);
  EXPECT_EQ(code.PoleOrderM(), m);

  const auto &pts = code.RationalPoints();
  ASSERT_EQ(pts.size(), expected_n);

  // Verify every point satisfies y^q + y == x^(q+1) and all points are
  // distinct.
  std::set<std::pair<typename GaloisField::element_type,
                     typename GaloisField::element_type>>
      unique_pts;
  for (const auto &pt : pts) {
    const auto lhs = field->Add(field->Pow(pt.y, expected_q), pt.y);
    const auto rhs = field->Mul(field->Pow(pt.x, expected_q), pt.x);
    EXPECT_EQ(lhs, rhs);
    unique_pts.insert({pt.x, pt.y});
  }
  EXPECT_EQ(unique_pts.size(), expected_n);

  // Verify exact O(1) dual: C_L(D, m P_inf)^perp = C_L(D, (n + 2g - 2 - m)
  // P_inf) and G * H^T == 0.
  auto dual = code.DualCode();
  EXPECT_EQ(dual.PoleOrderM(), expected_n + 2 * expected_genus - 2 - m);
  EXPECT_EQ(code.Dimension() + dual.Dimension(), expected_n);

  const auto prod = linalg::Dot(code.GeneratorMatrix(),
                                linalg::Transpose(dual.GeneratorMatrix()));
  EXPECT_TRUE(IsZeroMatrix(prod));

  const auto prod_h = linalg::Dot(code.GeneratorMatrix(),
                                  linalg::Transpose(code.ParityCheckMatrix()));
  EXPECT_TRUE(IsZeroMatrix(prod_h));
}

// ============================================================================
// Hermitian Curve & Dual Tests over GF(4), GF(9), GF(16)
// ============================================================================

TEST(HermitianCodeTest, CurvePointsAndExactDualGF4) {
  auto gf4 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(2);
  VerifyHermitianCurveAndDual(gf4, /*expected_q=*/2, /*expected_genus=*/1,
                              /*m=*/4);
  // Also test m >= q^3 (linear dependence removal via RowSpace).
  VerifyHermitianCurveAndDual(gf4, /*expected_q=*/2, /*expected_genus=*/1,
                              /*m=*/8);
}

TEST(HermitianCodeTest, CurvePointsAndExactDualGF9) {
  auto gf9 = std::make_shared<GaloisFieldExtension<uint32_t>>(
      std::pair<uint64_t, uint64_t>{3, 2});
  VerifyHermitianCurveAndDual(gf9, /*expected_q=*/3, /*expected_genus=*/3,
                              /*m=*/15);
  // Also test m >= q^3.
  VerifyHermitianCurveAndDual(gf9, /*expected_q=*/3, /*expected_genus=*/3,
                              /*m=*/29);
}

TEST(HermitianCodeTest, CurvePointsAndExactDualGF16) {
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  VerifyHermitianCurveAndDual(gf16, /*expected_q=*/4, /*expected_genus=*/6,
                              /*m=*/36);
}

TEST(HermitianCodeTest, RejectsOddDegreeFieldAndOutOfRangeM) {
  auto gf8 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);
  EXPECT_THROW((HermitianCode<GaloisFieldBinaryExtension<uint8_t>>(gf8, 2)),
               std::invalid_argument);

  auto gf4 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(2);
  // max_m = 8 + 2 - 1 = 9, so m = 9 is out of range.
  EXPECT_THROW((HermitianCode<GaloisFieldBinaryExtension<uint8_t>>(gf4, 9)),
               std::invalid_argument);
}

TEST(HermitianCodeTest, RuleOfFiveCopyAndMoveSafety) {
  auto gf4 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(2);
  HermitianCode<GaloisFieldBinaryExtension<uint8_t>> orig(gf4, 4);
  std::mt19937_64 rng(42);
  auto msg = GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(
      gf4, {orig.Dimension()}, rng);
  auto cw = orig.Encode(msg);

  HermitianCode<GaloisFieldBinaryExtension<uint8_t>> copied(orig);
  EXPECT_EQ(copied.Encode(msg), cw);
  auto dec_copy = copied.DecodeWithResult(cw);
  ASSERT_TRUE(dec_copy.Ok());
  EXPECT_EQ(dec_copy.codeword, cw);

  HermitianCode<GaloisFieldBinaryExtension<uint8_t>> moved(std::move(copied));
  auto dec_move = moved.DecodeWithResult(cw);
  ASSERT_TRUE(dec_move.Ok());
  EXPECT_EQ(dec_move.codeword, cw);
}

// ============================================================================
// HermitianFengRaoDecoder Tests (Full floor((d* - 1) / 2) via Majority Voting)
// ============================================================================

TEST(HermitianFengRaoTest, H3DecodesFullRadiusT5BeyondBasicBound) {
  // H_3 over GF(9): q = 3, n = 27, g = 3, m = 15 => d* = 12, t = 5.
  // Basic AG decoder radius is floor((d* - 1 - g) / 2) = floor(8 / 2) = 4,
  // so correcting t = 5 proves true Feng-Rao majority voting!
  using GF9 = GaloisFieldExtension<uint32_t>;
  auto gf9 = std::make_shared<GF9>(std::pair<uint64_t, uint64_t>{3, 2});
  HermitianCode<GF9> code(gf9, 15);

  EXPECT_EQ(code.Length(), 27U);
  EXPECT_EQ(code.Genus(), 3U);
  EXPECT_EQ(code.DesignedDistance(), 12U);
  auto decoder = code.GetDecoder("HermitianFengRao");
  EXPECT_EQ(decoder->DecodingRadius(), 5U);

  std::mt19937_64 rng(202601);
  for (std::size_t weight : {0U, 1U, 3U, 4U, 5U}) {
    const std::size_t num_trials = (weight == 5U) ? 10U : 3U;
    for (std::size_t trial = 0; trial < num_trials; ++trial) {
      const auto cw = code.RandomCodeword(rng);
      const auto received = CorruptRandomPositions(cw, weight, rng);
      const auto res = code.DecodeWithResult(received);
      ASSERT_TRUE(res.Ok()) << "Failed at weight=" << weight
                            << " trial=" << trial << ": " << res.error_message;
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.num_errors, weight);
    }
  }
}

TEST(HermitianFengRaoTest, H4DecodesFullRadiusT13BeyondBasicBound) {
  // H_4 over GF(16): q = 4, n = 64, g = 6, m = 36 => d* = 28, t = 13.
  // Basic AG decoder radius is floor((28 - 1 - 6) / 2) = 10, so t = 13 is 3
  // errors strictly beyond the basic decoder!
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);
  HermitianCode<GF16> code(gf16, 36);

  EXPECT_EQ(code.Length(), 64U);
  EXPECT_EQ(code.Genus(), 6U);
  EXPECT_EQ(code.DesignedDistance(), 28U);
  auto decoder = code.GetDecoder("HermitianFengRao");
  EXPECT_EQ(decoder->DecodingRadius(), 13U);

  std::mt19937_64 rng(202602);
  for (std::size_t weight : {0U, 6U, 11U, 13U}) {
    const std::size_t num_trials = (weight == 13U) ? 5U : 2U;
    for (std::size_t trial = 0; trial < num_trials; ++trial) {
      const auto cw = code.RandomCodeword(rng);
      const auto received = CorruptRandomPositions(cw, weight, rng);
      const auto res = code.DecodeWithResult(received);
      ASSERT_TRUE(res.Ok()) << "Failed at weight=" << weight
                            << " trial=" << trial << ": " << res.error_message;
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.num_errors, weight);
    }
  }
}

// ============================================================================
// HermitianSudanListDecoder Tests (Beyond Unique Decoding Radius)
// ============================================================================

TEST(HermitianSudanTest, H3ListDecodesAt9ErrorsBeyondUniqueRadius8) {
  // H_3 over GF(9): q = 3, n = 27, g = 3, m = 9 => d* = 18, t_unique = 8.
  // List decoding with target radius 9 (mu = 3, L = 5) recovers transmitted
  // codewords at t = 9 errors!
  using GF9 = GaloisFieldExtension<uint32_t>;
  auto gf9 = std::make_shared<GF9>(std::pair<uint64_t, uint64_t>{3, 2});
  HermitianCode<GF9> code(gf9, 9);

  EXPECT_EQ(code.DesignedDistance(), 18U);
  EXPECT_EQ(code.GetDecoder("HermitianFengRao")->DecodingRadius(), 8U);

  auto list_decoder = HermitianSudanListDecoder<GF9>::FromTargetRadius(code, 9);
  EXPECT_GE(list_decoder.DecodingRadius(), 9U);

  std::mt19937_64 rng(202603);
  for (std::size_t trial = 0; trial < 5; ++trial) {
    const auto cw = code.RandomCodeword(rng);
    const auto received = CorruptRandomPositions(cw, 9, rng);
    const auto list = list_decoder.DecodeToList(received);
    ASSERT_FALSE(list.empty()) << "Empty list on trial " << trial;
    bool found_transmitted = false;
    for (const auto &cand : list) {
      EXPECT_TRUE(code.Contains(cand));
      EXPECT_LE(code.Distance(cand, received), list_decoder.DecodingRadius());
      if (cand == cw) {
        found_transmitted = true;
      }
    }
    EXPECT_TRUE(found_transmitted)
        << "Transmitted codeword missing from Sudan list on trial " << trial;
  }
}

// ============================================================================
// LTCode Tests over GF(2) and GF(2^8)
// ============================================================================

TEST(LTCodeTest, PeelingAndGaussianFallbackGF2) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  constexpr std::size_t kK = 64;
  constexpr std::size_t kSymSize = 16;
  LTCode<GaloisFieldBinary> lt(
      gf2, kK, kSymSize,
      LTCode<GaloisFieldBinary>::DegreeDistribution::kRobustSoliton, 0.1, 0.05);

  std::mt19937_64 rng(1001);
  const auto source =
      GFMatrix<GaloisFieldBinary>::Random(gf2, {kK, kSymSize}, rng);

  // Peeling succeeds cleanly with Robust Soliton overhead (K * Z ~ 108).
  const auto syms_peel = lt.GenerateSymbols(source, 140, 0);
  const auto peel_res = lt.DecodePeeling(syms_peel);
  ASSERT_TRUE(peel_res.Ok()) << peel_res.error_message;
  EXPECT_FALSE(peel_res.used_gaussian_fallback);
  EXPECT_EQ(peel_res.source_block, source);

  // Gaussian fallback succeeds at much tighter overhead where peeling stalls.
  const auto syms_tight = lt.GenerateSymbols(source, 76, 500);
  const auto gauss_res = lt.DecodeWithGaussianFallback(syms_tight);
  ASSERT_TRUE(gauss_res.Ok()) << gauss_res.error_message;
  EXPECT_TRUE(gauss_res.used_gaussian_fallback);
  EXPECT_EQ(gauss_res.source_block, source);
}

TEST(LTCodeTest, PeelingAndGaussianFallbackGF256) {
  using GF256 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf256 = std::make_shared<GF256>(8);
  constexpr std::size_t kK = 64;
  constexpr std::size_t kSymSize = 16;
  LTCode<GF256> lt(gf256, kK, kSymSize,
                   LTCode<GF256>::DegreeDistribution::kRobustSoliton, 0.1,
                   0.05);

  std::mt19937_64 rng(1002);
  const auto source = GFMatrix<GF256>::Random(gf256, {kK, kSymSize}, rng);

  // Peeling decoder.
  const auto syms_peel = lt.GenerateSymbols(source, 140, 0);
  const auto peel_res = lt.DecodePeeling(syms_peel);
  ASSERT_TRUE(peel_res.Ok()) << peel_res.error_message;
  EXPECT_FALSE(peel_res.used_gaussian_fallback);
  EXPECT_EQ(peel_res.source_block, source);

  // Gaussian fallback at low overhead (72 symbols for K = 64).
  const auto syms_low = lt.GenerateSymbols(source, 72, 1000);
  const auto gauss_res = lt.DecodeWithGaussianFallback(syms_low);
  ASSERT_TRUE(gauss_res.Ok()) << gauss_res.error_message;
  EXPECT_TRUE(gauss_res.used_gaussian_fallback);
  EXPECT_EQ(gauss_res.source_block, source);
}

// ============================================================================
// RaptorQCode Tests over GF(2) and GF(2^8)
// ============================================================================

TEST(RaptorQCodeTest, SystematicPureRepairAndMixedGF2) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  constexpr std::size_t kK = 64;
  constexpr std::size_t kSymSize = 8;
  RaptorQCode<GaloisFieldBinary> rq(gf2, kK, kSymSize);

  std::mt19937_64 rng(2001);
  const auto source =
      GFMatrix<GaloisFieldBinary>::Random(gf2, {kK, kSymSize}, rng);
  const auto intermediate = rq.EncodeIntermediate(source);

  // 1. Systematic symbols (0 <= esi < K) match source verbatim.
  for (uint32_t esi = 0; esi < kK; ++esi) {
    const auto sym = rq.GenerateSymbol(intermediate, esi);
    for (std::size_t j = 0; j < kSymSize; ++j) {
      EXPECT_EQ(sym.payload.Raw(j), source.Raw(esi, j));
    }
  }

  // 2. Pure repair decoding (esi >= K).
  const auto pure_repair = rq.GenerateSymbols(intermediate, kK + 12, kK);
  const auto repair_res = rq.Decode(pure_repair);
  ASSERT_TRUE(repair_res.Ok()) << repair_res.error_message;
  EXPECT_EQ(repair_res.source_block, source);

  // 3. Mixed systematic + repair decoding.
  std::vector<FountainSymbol<GaloisFieldBinary>> mixed;
  for (uint32_t esi = 0; esi < kK; esi += 2) {
    mixed.push_back(rq.GenerateSymbol(intermediate, esi));
  }
  for (uint32_t esi = kK; esi < kK + 44; ++esi) {
    mixed.push_back(rq.GenerateSymbol(intermediate, esi));
  }
  const auto mixed_res = rq.Decode(mixed);
  ASSERT_TRUE(mixed_res.Ok()) << mixed_res.error_message;
  EXPECT_EQ(mixed_res.source_block, source);
}

TEST(RaptorQCodeTest, SystematicAndNearZeroOverheadGF256K64) {
  using GF256 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf256 = std::make_shared<GF256>(8);
  constexpr std::size_t kK = 64;
  constexpr std::size_t kSymSize = 8;
  RaptorQCode<GF256> rq(gf256, kK, kSymSize);

  std::mt19937_64 rng(2002);
  const auto source = GFMatrix<GF256>::Random(gf256, {kK, kSymSize}, rng);
  const auto intermediate = rq.EncodeIntermediate(source);

  // 1. Systematic symbols match source verbatim.
  for (uint32_t esi = 0; esi < kK; ++esi) {
    const auto sym = rq.GenerateSymbol(intermediate, esi);
    for (std::size_t j = 0; j < kSymSize; ++j) {
      EXPECT_EQ(sym.payload.Raw(j), source.Raw(esi, j));
    }
  }

  // 2. Pure repair decoding at +2 overhead.
  const auto pure_repair = rq.GenerateSymbols(intermediate, kK + 2, kK);
  const auto repair_res = rq.Decode(pure_repair);
  ASSERT_TRUE(repair_res.Ok()) << repair_res.error_message;
  EXPECT_EQ(repair_res.source_block, source);

  // 3. Empirical > 98% success rate at +2 symbols overhead across 50 random
  // erasure trials (mixed systematic + repair pool).
  const auto pool = rq.GenerateSymbols(intermediate, 2 * kK, 0);
  std::vector<std::size_t> idxs(pool.size());
  std::iota(idxs.begin(), idxs.end(), std::size_t{0});

  constexpr std::size_t kTrials = 50;
  std::size_t successes = 0;
  for (std::size_t t = 0; t < kTrials; ++t) {
    std::shuffle(idxs.begin(), idxs.end(), rng);
    std::vector<FountainSymbol<GF256>> subset;
    subset.reserve(kK + 2);
    for (std::size_t i = 0; i < kK + 2; ++i) {
      subset.push_back(pool[idxs[i]]);
    }
    const auto dec = rq.Decode(subset);
    if (dec.Ok() && dec.source_block == source) {
      ++successes;
    }
  }
  EXPECT_GE(successes, 49U) << "Expected > 98% (>= 49/50) at +2 overhead, got "
                            << successes << "/" << kTrials;
}

TEST(RaptorQCodeTest, NearZeroOverheadGF256K128) {
  using GF256 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf256 = std::make_shared<GF256>(8);
  constexpr std::size_t kK = 128;
  constexpr std::size_t kSymSize = 8;
  RaptorQCode<GF256> rq(gf256, kK, kSymSize);

  std::mt19937_64 rng(2003);
  const auto source = GFMatrix<GF256>::Random(gf256, {kK, kSymSize}, rng);
  const auto intermediate = rq.EncodeIntermediate(source);

  const auto pool = rq.GenerateSymbols(intermediate, 2 * kK, 0);
  std::vector<std::size_t> idxs(pool.size());
  std::iota(idxs.begin(), idxs.end(), std::size_t{0});

  constexpr std::size_t kTrials = 50;
  std::size_t successes = 0;
  for (std::size_t t = 0; t < kTrials; ++t) {
    std::shuffle(idxs.begin(), idxs.end(), rng);
    std::vector<FountainSymbol<GF256>> subset;
    subset.reserve(kK + 2);
    for (std::size_t i = 0; i < kK + 2; ++i) {
      subset.push_back(pool[idxs[i]]);
    }
    const auto dec = rq.Decode(subset);
    if (dec.Ok() && dec.source_block == source) {
      ++successes;
    }
  }
  EXPECT_GE(successes, 49U) << "Expected > 98% (>= 49/50) at +2 overhead, got "
                            << successes << "/" << kTrials;
}

}  // namespace
