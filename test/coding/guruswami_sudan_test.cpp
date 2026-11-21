#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace {

using namespace xg;
using namespace xg::coding;

template <typename GaloisField>
GFArray<GaloisField> MakeGFArrayFromInts(
    const std::shared_ptr<GaloisField> &field,
    const std::vector<uint64_t> &vals) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;
  RawArray raw = RawArray::from_shape(ShapeType{vals.size()});
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if constexpr (std::is_integral_v<StorageType>) {
      raw(i) = field->SetElementValue(static_cast<StorageType>(vals[i]));
    } else {
      raw(i) = field->FromInt(vals[i]);
    }
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
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
GFArray<GaloisField> CorruptAtPositions(
    const GFArray<GaloisField> &cw, const std::vector<std::size_t> &positions,
    std::mt19937_64 &rng) {
  const auto field = cw.Field();
  auto received = cw;
  for (std::size_t pos : positions) {
    const auto err = RandomNonZero(field, rng);
    received.Raw(pos) = field->Add(received.Raw(pos), err);
  }
  return received;
}

template <typename GaloisField>
bool ContainsCodeword(const std::vector<GFArray<GaloisField>> &list,
                      const GFArray<GaloisField> &target) {
  return std::find(list.begin(), list.end(), target) != list.end();
}

template <typename GaloisField>
void VerifyListInvariants(const GeneralizedReedSolomonCode<GaloisField> &code,
                          const GuruswamiSudanDecoder<GaloisField> &decoder,
                          const GFArray<GaloisField> &received,
                          const GFArray<GaloisField> &expected_cw,
                          const GFArray<GaloisField> &expected_msg) {
  const auto cw_list = decoder.DecodeToList(received);
  const auto msg_list = decoder.DecodeToMessageList(received);
  ASSERT_EQ(cw_list.size(), msg_list.size());
  ASSERT_FALSE(cw_list.empty());
  EXPECT_TRUE(ContainsCodeword(cw_list, expected_cw));
  EXPECT_TRUE(ContainsCodeword(msg_list, expected_msg));

  std::size_t prev_dist = 0;
  for (std::size_t i = 0; i < cw_list.size(); ++i) {
    EXPECT_TRUE(code.Contains(cw_list[i]));
    EXPECT_EQ(code.Unencode(cw_list[i]), msg_list[i]);
    const std::size_t dist = code.Distance(cw_list[i], received);
    EXPECT_LE(dist, decoder.DecodingRadius());
    if (i > 0) {
      EXPECT_GE(dist, prev_dist);
    }
    prev_dist = dist;
  }

  const auto res = decoder.DecodeWithResult(received);
  ASSERT_TRUE(res.Ok());
  EXPECT_TRUE(code.Contains(res.codeword));
  EXPECT_EQ(res.codeword, cw_list[0]);
  EXPECT_EQ(res.message, msg_list[0]);
  EXPECT_EQ(res.num_errors, code.Distance(cw_list[0], received));
  EXPECT_EQ(res.num_erasures, 0U);
}

// Verify beyond-half-minimum-distance decoding for RS(15, 3) over GF(2^4):
// d = 13, unique radius t = 6, Johnson radius = 15 - 1 - floor(sqrt(30)) = 9.
TEST(GuruswamiSudanTest, RS15_3_GF16_BeyondUniqueRadius_7_8_9_Errors) {
  using GF = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF>(4);
  auto rs = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf16, 15, 3);

  ASSERT_EQ(rs.MinimumDistance(), 13U);
  const std::size_t unique_t = (rs.MinimumDistance() - 1) / 2;
  EXPECT_EQ(unique_t, 6U);

  auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 9);
  EXPECT_EQ(decoder.JohnsonRadius(), 9U);
  EXPECT_EQ(decoder.DecodingRadius(), 9U);
  EXPECT_EQ(decoder.Name(), "GuruswamiSudan");

  std::mt19937_64 rng(202601);
  const std::vector<GFArray<GF>> messages = {
      GFArray<GF>::Zeros(gf16, {3}),
      GFArray<GF>::Random(gf16, {3}, rng),
      GFArray<GF>::Random(gf16, {3}, rng),
  };

  std::vector<std::size_t> all_positions(15);
  std::iota(all_positions.begin(), all_positions.end(), std::size_t{0});

  for (const auto &msg : messages) {
    const auto cw = rs.Encode(msg);
    ASSERT_TRUE(rs.Contains(cw));
    EXPECT_EQ(rs.Unencode(cw), msg);

    for (std::size_t num_err :
         {std::size_t{7}, std::size_t{8}, std::size_t{9}}) {
      std::shuffle(all_positions.begin(), all_positions.end(), rng);
      const std::vector<std::size_t> err_pos(
          all_positions.begin(),
          all_positions.begin() + static_cast<std::ptrdiff_t>(num_err));
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      ASSERT_EQ(rs.Distance(cw, received), num_err);

      // Unique decoders (radius t = 6) cannot uniquely recover cw from >= 7
      // errors (either fail outright or misdecode to a different codeword).
      const auto bm_res = rs.DecodeWithResult(received, "BerlekampMassey");
      const auto sug_res = rs.DecodeWithResult(received, "Sugiyama");
      EXPECT_FALSE(bm_res.Ok() && bm_res.codeword == cw);
      EXPECT_FALSE(sug_res.Ok() && sug_res.codeword == cw);

      // Guruswami-Sudan recovers the exact transmitted codeword and message!
      VerifyListInvariants(rs, decoder, received, cw, msg);
    }
  }
}

// Verify RS(15, 3) over GF2XZECH<uint8_t> (non-polynomial internal Zech log
// representation) decodes 7, 8, and 9 errors cleanly.
TEST(GuruswamiSudanTest, RS15_3_GF16ZechRepresentation_9Errors) {
  using GF = GF2XZECH<uint8_t>;
  auto gf16 = std::make_shared<GF>(4);
  auto rs = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf16, 15, 3);

  auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 9);
  EXPECT_EQ(decoder.JohnsonRadius(), 9U);
  EXPECT_EQ(decoder.DecodingRadius(), 9U);

  std::mt19937_64 rng(202602);
  const auto msg = GFArray<GF>::Random(gf16, {3}, rng);
  const auto cw = rs.Encode(msg);

  const std::vector<std::size_t> err_pos = {0, 2, 4, 5, 7, 9, 11, 13, 14};
  const auto received = CorruptAtPositions(cw, err_pos, rng);
  ASSERT_EQ(rs.Distance(cw, received), 9U);

  VerifyListInvariants(rs, decoder, received, cw, msg);
}

// Verify extended RS(16, 4) over GF(2^4) with alpha_0 = 0:
// d = 13, unique radius t = 6, Johnson radius = 16 - 1 - floor(sqrt(48)) = 9.
// Verify recovery with 7 and 8 errors, explicitly corrupting alpha_0 = 0!
TEST(GuruswamiSudanTest, ExtendedRS16_4_GF16_WithZeroEvalPoint_7_8_Errors) {
  using GF = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF>(4);
  auto rs = GeneralizedReedSolomonCode<GF>::ExtendedReedSolomon(gf16, 16, 4);

  ASSERT_TRUE(rs.ZeroEvaluationIndex().has_value());
  ASSERT_EQ(*rs.ZeroEvaluationIndex(), 0U);
  ASSERT_EQ(rs.MinimumDistance(), 13U);

  auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 8);
  EXPECT_EQ(decoder.JohnsonRadius(), 9U);
  EXPECT_EQ(decoder.DecodingRadius(), 8U);

  std::mt19937_64 rng(202603);
  const std::vector<GFArray<GF>> messages = {
      GFArray<GF>::Zeros(gf16, {4}),
      GFArray<GF>::Random(gf16, {4}, rng),
  };

  for (const auto &msg : messages) {
    const auto cw = rs.Encode(msg);

    // 7 errors including index 0 (alpha_0 = 0).
    {
      const std::vector<std::size_t> err_pos = {0, 1, 3, 6, 9, 12, 15};
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      ASSERT_EQ(rs.Distance(cw, received), 7U);

      const auto bm_res = rs.DecodeWithResult(received, "BerlekampMassey");
      EXPECT_FALSE(bm_res.Ok() && bm_res.codeword == cw);

      VerifyListInvariants(rs, decoder, received, cw, msg);
    }

    // 8 errors including index 0 (alpha_0 = 0).
    {
      const std::vector<std::size_t> err_pos = {0, 2, 4, 5, 8, 10, 13, 14};
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      ASSERT_EQ(rs.Distance(cw, received), 8U);

      const auto sug_res = rs.DecodeWithResult(received, "Sugiyama");
      EXPECT_FALSE(sug_res.Ok() && sug_res.codeword == cw);

      VerifyListInvariants(rs, decoder, received, cw, msg);
    }
  }
}

// Verify RS(12, 2) over odd-characteristic prime field GF(13):
// d = 11, unique radius t = 5, Johnson radius = 12 - 1 - floor(sqrt(12)) = 8.
// Exercises exact Hasse derivatives modulo p = 13 with multiplicity m = 3.
TEST(GuruswamiSudanTest, RS12_2_GF13_OddPrimeField_6_7_8_Errors) {
  using GF = GaloisFieldPrime<uint8_t>;
  auto gf13 = std::make_shared<GF>(13);
  auto rs = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf13, 12, 2);

  ASSERT_EQ(rs.MinimumDistance(), 11U);
  auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 8);
  EXPECT_EQ(decoder.JohnsonRadius(), 8U);
  EXPECT_EQ(decoder.DecodingRadius(), 8U);
  EXPECT_EQ(decoder.Multiplicity(), 3U);

  std::mt19937_64 rng(202604);
  const std::vector<GFArray<GF>> messages = {
      GFArray<GF>::Zeros(gf13, {2}),
      GFArray<GF>::Random(gf13, {2}, rng),
      GFArray<GF>::Random(gf13, {2}, rng),
  };

  std::vector<std::size_t> all_positions(12);
  std::iota(all_positions.begin(), all_positions.end(), std::size_t{0});

  for (const auto &msg : messages) {
    const auto cw = rs.Encode(msg);
    for (std::size_t num_err :
         {std::size_t{6}, std::size_t{7}, std::size_t{8}}) {
      std::shuffle(all_positions.begin(), all_positions.end(), rng);
      const std::vector<std::size_t> err_pos(
          all_positions.begin(),
          all_positions.begin() + static_cast<std::ptrdiff_t>(num_err));
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      ASSERT_EQ(rs.Distance(cw, received), num_err);

      const auto bm_res = rs.DecodeWithResult(received, "BerlekampMassey");
      EXPECT_FALSE(bm_res.Ok() && bm_res.codeword == cw);

      VerifyListInvariants(rs, decoder, received, cw, msg);
    }
  }
}

// Verify GRS(15, 3) with non-trivial column multipliers v_j != 1 over GF(2^4),
// including both nonzero evaluation points and an extended variant with
// alpha_0 = 0 and v_0 != 1.
TEST(GuruswamiSudanTest, GRS15_3_GF16_NonTrivialColumnMultipliers) {
  using GF = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF>(4);
  auto alphas = MakeGFArrayFromInts(
      gf16, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
  auto vs = MakeGFArrayFromInts(
      gf16, {5, 11, 3, 14, 7, 2, 13, 9, 6, 15, 4, 8, 12, 10, 3});
  GeneralizedReedSolomonCode<GF> grs(gf16, alphas, vs, 3);

  ASSERT_TRUE(grs.ZeroEvaluationIndex().has_value());
  ASSERT_EQ(grs.MinimumDistance(), 13U);

  auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(grs, 9);
  EXPECT_EQ(decoder.JohnsonRadius(), 9U);
  EXPECT_EQ(decoder.DecodingRadius(), 9U);

  std::mt19937_64 rng(202605);
  const std::vector<GFArray<GF>> messages = {
      GFArray<GF>::Zeros(gf16, {3}),
      GFArray<GF>::Random(gf16, {3}, rng),
  };

  for (const auto &msg : messages) {
    const auto cw = grs.Encode(msg);
    ASSERT_TRUE(grs.Contains(cw));

    for (std::size_t num_err :
         {std::size_t{7}, std::size_t{8}, std::size_t{9}}) {
      // Always include position 0 (alpha_0 = 0, v_0 = 5) among the errors.
      std::vector<std::size_t> err_pos = {0};
      for (std::size_t j = 1; err_pos.size() < num_err; j += 1) {
        if (j < 15) {
          err_pos.push_back(j);
        }
      }
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      ASSERT_EQ(grs.Distance(cw, received), num_err);

      VerifyListInvariants(grs, decoder, received, cw, msg);
    }
  }
}

// Verify RS(9, 2) over general odd-characteristic extension field GF(3^2):
// d = 8, unique radius t = 3, Johnson radius = 9 - 1 - floor(sqrt(9)) = 5.
// Tests both m = 1 (Sudan, achieves radius 5) and m = 2 (exercises higher-order
// Hasse derivatives in characteristic 3).
TEST(GuruswamiSudanTest, RS9_2_GF9_OddExtensionField_4_5_Errors) {
  using GF = GaloisFieldExtension<uint8_t>;
  auto gf9 = std::make_shared<GF>(std::pair<uint64_t, uint64_t>{3, 2});
  auto rs = GeneralizedReedSolomonCode<GF>::ExtendedReedSolomon(gf9, 9, 2);

  ASSERT_TRUE(rs.ZeroEvaluationIndex().has_value());
  ASSERT_EQ(rs.MinimumDistance(), 8U);

  auto dec_auto = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 5);
  EXPECT_EQ(dec_auto.JohnsonRadius(), 5U);
  EXPECT_EQ(dec_auto.DecodingRadius(), 5U);

  // Explicit multiplicity m = 2 exercises characteristic-3 Hasse derivatives.
  GuruswamiSudanDecoder<GF> dec_m2(rs, 2);
  EXPECT_EQ(dec_m2.Multiplicity(), 2U);
  EXPECT_EQ(dec_m2.DecodingRadius(), 5U);

  std::mt19937_64 rng(202606);
  const std::vector<GFArray<GF>> messages = {
      GFArray<GF>::Zeros(gf9, {2}),
      GFArray<GF>::Random(gf9, {2}, rng),
  };

  for (const auto &msg : messages) {
    const auto cw = rs.Encode(msg);

    for (const std::vector<std::size_t> &err_pos :
         {std::vector<std::size_t>{0, 2, 5, 8},
          std::vector<std::size_t>{0, 1, 3, 6, 8}}) {
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      ASSERT_EQ(rs.Distance(cw, received), err_pos.size());

      const auto bm_res = rs.DecodeWithResult(received, "BerlekampMassey");
      EXPECT_FALSE(bm_res.Ok() && bm_res.codeword == cw);

      VerifyListInvariants(rs, dec_auto, received, cw, msg);
      VerifyListInvariants(rs, dec_m2, received, cw, msg);
    }
  }
}

// Genuine multi-codeword list tests:
// 1. Two codewords equidistant (distance 7) from received word in RS(15, 3).
// 2. Three codewords equidistant (distance 8) from received word in RS(12, 2).
// 3. Two codewords at distinct distances (6 and 7) verifying ascending sort.
TEST(GuruswamiSudanTest, GenuineMultiCodewordListsAndSorting) {
  // Part 1: Two equidistant codewords at distance 7 in RS(15, 3) over GF(2^4).
  {
    using GF = GaloisFieldBinaryExtension<uint8_t>;
    using Poly = PolynomialDense<GF>;
    auto gf16 = std::make_shared<GF>(4);
    auto rs = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf16, 15, 3);

    const Poly x_poly = Poly::Monomial(gf16, 1, gf16->MultiplicativeIdentity());
    const Poly a0 = Poly::FromRaw({rs.EvaluationPoints().Raw(0)}, gf16);
    const Poly a1 = Poly::FromRaw({rs.EvaluationPoints().Raw(1)}, gf16);

    const Poly f1 = Poly::Zero(gf16);
    const Poly f2 = (x_poly - a0) * (x_poly - a1);
    const auto c1 = rs.EncodePolynomial(f1);
    const auto c2 = rs.EncodePolynomial(f2);
    ASSERT_EQ(rs.Distance(c1, c2), 13U);

    // Pick a symbol at index 14 differing from both c1[14] and c2[14].
    typename GF::element_type third_sym = gf16->MultiplicativeIdentity();
    for (uint8_t val = 1; val < 16; ++val) {
      const auto cand = gf16->SetElementValue(val);
      if (cand != c1.Raw(14) && cand != c2.Raw(14)) {
        third_sym = cand;
        break;
      }
    }

    // Build r agreeing with both on {0, 1}, with c1 on {2..7}, with c2 on
    // {8..13}, and with neither on {14} -> dist(r, c1) = dist(r, c2) = 7.
    auto received = c1;
    for (std::size_t j = 8; j <= 13; ++j) {
      received.Raw(j) = c2.Raw(j);
    }
    received.Raw(14) = third_sym;
    ASSERT_EQ(rs.Distance(received, c1), 7U);
    ASSERT_EQ(rs.Distance(received, c2), 7U);

    auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 7);
    const auto list = decoder.DecodeToList(received);
    // Any codeword c3 != c1, c2 agrees with c1 on <= 2 positions, with c2 on
    // <= 2 positions, and with received[14] on <= 1 position -> dist >= 10.
    // Hence {c1, c2} is the exact complete list within radius 7.
    ASSERT_EQ(list.size(), 2U);
    EXPECT_TRUE(ContainsCodeword(list, c1));
    EXPECT_TRUE(ContainsCodeword(list, c2));

    const auto msg_list = decoder.DecodeToMessageList(received);
    ASSERT_EQ(msg_list.size(), 2U);
    EXPECT_TRUE(ContainsCodeword(msg_list, rs.Unencode(c1)));
    EXPECT_TRUE(ContainsCodeword(msg_list, rs.Unencode(c2)));

    const auto res = decoder.DecodeWithResult(received);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.num_errors, 7U);
    EXPECT_TRUE(res.codeword == c1 || res.codeword == c2);
    EXPECT_FALSE(res.error_message.empty());

    // Part 3: Two codewords at distinct distances 6 and 7 -> c1 must sort
    // strictly before c2 and be chosen uniquely by DecodeWithResult.
    auto received_asym = c1;
    for (std::size_t j = 9; j <= 14; ++j) {
      received_asym.Raw(j) = c2.Raw(j);
    }
    ASSERT_EQ(rs.Distance(received_asym, c1), 6U);
    ASSERT_EQ(rs.Distance(received_asym, c2), 7U);

    const auto list_asym = decoder.DecodeToList(received_asym);
    ASSERT_EQ(list_asym.size(), 2U);
    EXPECT_EQ(list_asym[0], c1);
    EXPECT_EQ(list_asym[1], c2);

    const auto res_asym = decoder.DecodeWithResult(received_asym);
    ASSERT_TRUE(res_asym.Ok());
    EXPECT_EQ(res_asym.codeword, c1);
    EXPECT_EQ(res_asym.num_errors, 6U);
    EXPECT_TRUE(res_asym.error_message.empty());
  }

  // Part 2: Three equidistant codewords at distance 8 in RS(12, 2) over GF(13).
  {
    using GF = GaloisFieldPrime<uint8_t>;
    using Poly = PolynomialDense<GF>;
    auto gf13 = std::make_shared<GF>(13);
    auto rs = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf13, 12, 2);

    const auto c0 = rs.EncodePolynomial(Poly::Zero(gf13));
    const auto c1 =
        rs.EncodePolynomial(Poly::FromRaw({gf13->SetElementValue(1)}, gf13));
    const auto c2 =
        rs.EncodePolynomial(Poly::FromRaw({gf13->SetElementValue(2)}, gf13));

    auto received = c0;
    for (std::size_t j = 4; j < 8; ++j) {
      received.Raw(j) = c1.Raw(j);
    }
    for (std::size_t j = 8; j < 12; ++j) {
      received.Raw(j) = c2.Raw(j);
    }
    ASSERT_EQ(rs.Distance(received, c0), 8U);
    ASSERT_EQ(rs.Distance(received, c1), 8U);
    ASSERT_EQ(rs.Distance(received, c2), 8U);

    auto decoder = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 8);
    const auto list = decoder.DecodeToList(received);
    ASSERT_GE(list.size(), 3U);
    EXPECT_TRUE(ContainsCodeword(list, c0));
    EXPECT_TRUE(ContainsCodeword(list, c1));
    EXPECT_TRUE(ContainsCodeword(list, c2));

    const auto poly_list = decoder.DecodeToPolynomialList(received);
    EXPECT_EQ(poly_list.size(), list.size());
  }
}

// Verify FromTargetRadius sizing, JohnsonRadius, DecodingRadius, standard
// 0..t recovery, uncorrectable Failure reporting, k = 1 repetition edge case,
// and invalid argument validation.
TEST(GuruswamiSudanTest, ParameterSizingBoundsStandardErrorsAndEdgeCases) {
  using GF = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF>(4);
  auto rs = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf16, 15, 3);

  // Verify exact JohnsonRadius static helper across edge parameters.
  EXPECT_EQ(GuruswamiSudanDecoder<GF>::ComputeJohnsonRadius(15, 3), 9U);
  EXPECT_EQ(GuruswamiSudanDecoder<GF>::ComputeJohnsonRadius(16, 4), 9U);
  EXPECT_EQ(GuruswamiSudanDecoder<GF>::ComputeJohnsonRadius(12, 2), 8U);
  EXPECT_EQ(GuruswamiSudanDecoder<GF>::ComputeJohnsonRadius(9, 2), 5U);
  EXPECT_EQ(GuruswamiSudanDecoder<GF>::ComputeJohnsonRadius(7, 1), 6U);
  EXPECT_EQ(GuruswamiSudanDecoder<GF>::ComputeJohnsonRadius(15, 15), 0U);

  // FromTargetRadius throws when target exceeds JohnsonRadius.
  EXPECT_THROW(GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, 10),
               std::invalid_argument);
  EXPECT_THROW(GuruswamiSudanDecoder<GF>(rs, 0), std::invalid_argument);
  EXPECT_THROW(GuruswamiSudanDecoder<GF>(rs, 1, 0), std::invalid_argument);

  // Verify minimal parameter selection for every target radius 0..9 on
  // RS(15, 3).
  for (std::size_t target = 0; target <= 9; ++target) {
    auto dec = GuruswamiSudanDecoder<GF>::FromTargetRadius(rs, target);
    EXPECT_EQ(dec.JohnsonRadius(), 9U);
    EXPECT_GE(dec.DecodingRadius(), target);
  }

  // Verify default (m = 1, Sudan) constructor on RS(15, 3): achieves radius 8.
  GuruswamiSudanDecoder<GF> sudan(rs);
  EXPECT_EQ(sudan.Multiplicity(), 1U);
  EXPECT_EQ(sudan.DecodingRadius(), 8U);

  // Verify Welch-Berlekamp specialization (m = 1, L = 1) achieves unique
  // radius t = 6 and recovers all 0..6 errors including 0 errors & all-zero
  // message, while returning explicit Failure when no codeword is within
  // radius 6.
  GuruswamiSudanDecoder<GF> wb(rs, 1, 1);
  EXPECT_EQ(wb.DecodingRadius(), 6U);

  std::mt19937_64 rng(202607);
  std::vector<std::size_t> all_pos(15);
  std::iota(all_pos.begin(), all_pos.end(), std::size_t{0});

  for (const auto &msg :
       {GFArray<GF>::Zeros(gf16, {3}), GFArray<GF>::Random(gf16, {3}, rng)}) {
    const auto cw = rs.Encode(msg);
    for (std::size_t num_err = 0; num_err <= 6; ++num_err) {
      std::shuffle(all_pos.begin(), all_pos.end(), rng);
      const std::vector<std::size_t> err_pos(
          all_pos.begin(),
          all_pos.begin() + static_cast<std::ptrdiff_t>(num_err));
      const auto received = CorruptAtPositions(cw, err_pos, rng);
      const auto res = wb.DecodeWithResult(received);
      ASSERT_TRUE(res.Ok()) << "Failed at num_err=" << num_err;
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, num_err);
    }
  }

  // Uncorrectable word test: construct a received word whose distance to c1
  // and c2 is 7 and whose distance to every other codeword in RS(15, 3) is
  // >= 10, guaranteeing zero codewords within DecodingRadius() == 6.
  {
    using Poly = PolynomialDense<GF>;
    const Poly x_poly = Poly::Monomial(gf16, 1, gf16->MultiplicativeIdentity());
    const Poly a0 = Poly::FromRaw({rs.EvaluationPoints().Raw(0)}, gf16);
    const Poly a1 = Poly::FromRaw({rs.EvaluationPoints().Raw(1)}, gf16);
    const auto c1 = rs.EncodePolynomial(Poly::Zero(gf16));
    const auto c2 = rs.EncodePolynomial((x_poly - a0) * (x_poly - a1));

    typename GF::element_type third_sym = gf16->MultiplicativeIdentity();
    for (uint8_t val = 1; val < 16; ++val) {
      const auto cand = gf16->SetElementValue(val);
      if (cand != c1.Raw(14) && cand != c2.Raw(14)) {
        third_sym = cand;
        break;
      }
    }

    auto received = c1;
    for (std::size_t j = 8; j <= 13; ++j) {
      received.Raw(j) = c2.Raw(j);
    }
    received.Raw(14) = third_sym;
    ASSERT_EQ(rs.Distance(received, c1), 7U);
    ASSERT_EQ(rs.Distance(received, c2), 7U);

    const auto wb_res = wb.DecodeWithResult(received);
    EXPECT_FALSE(wb_res.Ok());
    EXPECT_FALSE(wb_res.error_message.empty());
    EXPECT_TRUE(wb.DecodeToList(received).empty());
  }

  // k = 1 repetition GRS code edge case: RS(7, 1) over GF(2^3).
  {
    auto gf8 = std::make_shared<GF>(3);
    auto rep = GeneralizedReedSolomonCode<GF>::ReedSolomon(gf8, 7, 1);
    auto rep_dec = GuruswamiSudanDecoder<GF>::FromTargetRadius(rep, 6);
    EXPECT_EQ(rep_dec.JohnsonRadius(), 6U);
    EXPECT_EQ(rep_dec.DecodingRadius(), 6U);

    // Received word with symbol 3 appearing 4 times, symbol 5 twice, symbol 1
    // once -> all three constant codewords are within radius 6, sorted by
    // distance 3, 5, 6!
    auto received = MakeGFArrayFromInts(gf8, {3, 3, 5, 3, 1, 5, 3});
    const auto list = rep_dec.DecodeToList(received);
    ASSERT_EQ(list.size(), 3U);
    EXPECT_EQ(list[0], MakeGFArrayFromInts(gf8, {3, 3, 3, 3, 3, 3, 3}));
    EXPECT_EQ(list[1], MakeGFArrayFromInts(gf8, {5, 5, 5, 5, 5, 5, 5}));
    EXPECT_EQ(list[2], MakeGFArrayFromInts(gf8, {1, 1, 1, 1, 1, 1, 1}));

    const auto res = rep_dec.DecodeWithResult(received);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, list[0]);
    EXPECT_EQ(res.num_errors, 3U);

    // All 7 symbols distinct -> distance >= 6 from every codeword, so radius-3
    // decoder must return Failure.
    GuruswamiSudanDecoder<GF> rep_wb(rep, 1, 1);
    EXPECT_EQ(rep_wb.DecodingRadius(), 3U);
    auto all_distinct = MakeGFArrayFromInts(gf8, {0, 1, 2, 3, 4, 5, 6});
    const auto fail_res = rep_wb.DecodeWithResult(all_distinct);
    EXPECT_FALSE(fail_res.Ok());
    EXPECT_TRUE(rep_wb.DecodeToList(all_distinct).empty());
  }

  // Invalid received word dimensions / field mismatch throw invalid_argument.
  auto wrong_len = GFArray<GF>::Zeros(gf16, {10});
  EXPECT_THROW(sudan.DecodeToList(wrong_len), std::invalid_argument);
  auto gf8_other = std::make_shared<GF>(3);
  auto wrong_field = GFArray<GF>::Zeros(gf8_other, {15});
  EXPECT_THROW(sudan.DecodeToList(wrong_field), std::invalid_argument);
}

}  // namespace
