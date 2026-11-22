#include "xgalois/coding/interleaved_rs.hpp"

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
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/linear_code.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
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
GFArray<GaloisField> CorruptWordAtPositions(
    const GFArray<GaloisField> &cw, const std::vector<std::size_t> &positions,
    std::mt19937_64 &rng) {
  auto field = cw.Field();
  auto received = cw;
  for (std::size_t pos : positions) {
    const auto err = RandomNonZero(field, rng);
    received.Raw(pos) = field->Add(received.Raw(pos), err);
  }
  return received;
}

// Injects non-zero burst error rows (at least one non-zero entry per burst
// position; in fact every column receives a random field element with the row
// guaranteed non-zero) at the specified row indices.
template <typename GaloisField>
GFMatrix<GaloisField> InjectBurstErrors(
    const GFMatrix<GaloisField> &cw_mat,
    const std::vector<std::size_t> &burst_positions, std::mt19937_64 &rng) {
  auto field = cw_mat.Field();
  const std::size_t ell = cw_mat.Shape(1);
  auto rx = cw_mat;
  for (std::size_t row : burst_positions) {
    // Draw a uniformly random non-zero vector in F_q^ell \ {0}.
    std::vector<typename GaloisField::element_type> err_row(ell);
    bool any_nz = false;
    while (!any_nz) {
      for (std::size_t m = 0; m < ell; ++m) {
        err_row[m] = field->Random(rng);
        if (!field->IsZero(err_row[m])) {
          any_nz = true;
        }
      }
    }
    for (std::size_t m = 0; m < ell; ++m) {
      rx.Raw(row, m) = field->Add(rx.Raw(row, m), err_row[m]);
    }
  }
  return rx;
}

// Injects burst errors where EVERY column at every burst position is non-zero
// (ensuring each individual column suffers exactly t errors so independent BM
// strictly fails whenever t > t_u).
template <typename GaloisField>
GFMatrix<GaloisField> InjectFullColumnBurstErrors(
    const GFMatrix<GaloisField> &cw_mat,
    const std::vector<std::size_t> &burst_positions, std::mt19937_64 &rng) {
  auto field = cw_mat.Field();
  const std::size_t ell = cw_mat.Shape(1);
  auto rx = cw_mat;
  for (std::size_t row : burst_positions) {
    for (std::size_t m = 0; m < ell; ++m) {
      const auto err = RandomNonZero(field, rng);
      rx.Raw(row, m) = field->Add(rx.Raw(row, m), err);
    }
  }
  return rx;
}

TEST(InterleavedRSTest, HomogeneousIRS15_7_3BeyondUniqueRadius) {
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);

  auto irs = InterleavedGRSCode<GF16>::CreateHomogeneous(gf16, 15, 7, 3);
  EXPECT_EQ(irs.Field(), gf16);
  EXPECT_EQ(irs.Length(), 15U);
  EXPECT_EQ(irs.InterleavingOrder(), 3U);
  EXPECT_TRUE(irs.IsHomogeneous());
  EXPECT_EQ(irs.Dimensions(), (std::vector<std::size_t>{7, 7, 7}));
  EXPECT_EQ(irs.UniqueDecodingRadius(), 4U);
  EXPECT_EQ(irs.CollaborativeDecodingRadius(), 6U);
  EXPECT_EQ(irs.ComponentCode(0).Dimension(), 7U);
  EXPECT_EQ(irs.ComponentCode(2).Length(), 15U);
  EXPECT_THROW(static_cast<void>(irs.ComponentCode(3)), std::out_of_range);

  std::mt19937_64 rng(20260401ULL);

  // Verify Encode vs EncodeHomogeneous, Contains, Syndromes, and ToLinearCode.
  auto msg_mat = GFMatrix<GF16>::Random(gf16, {7, 3}, rng);
  std::vector<GFArray<GF16>> msgs;
  for (std::size_t m = 0; m < 3; ++m) {
    using RawArray = typename GFArray<GF16>::raw_array_type;
    using ShapeType = typename GFArray<GF16>::shape_type;
    RawArray r = RawArray::from_shape(ShapeType{7});
    for (std::size_t i = 0; i < 7; ++i) {
      r(i) = msg_mat.Raw(i, m);
    }
    msgs.push_back(GFArray<GF16>::FromRaw(gf16, std::move(r)));
  }

  const auto cw_mat = irs.Encode(msgs);
  const auto cw_mat_hom = irs.EncodeHomogeneous(msg_mat);
  EXPECT_EQ(cw_mat, cw_mat_hom);
  EXPECT_TRUE(irs.Contains(cw_mat));

  const auto zero_syns = irs.Syndromes(cw_mat);
  ASSERT_EQ(zero_syns.size(), 3U);
  for (std::size_t m = 0; m < 3; ++m) {
    ASSERT_EQ(zero_syns[m].Size(), 8U);
    for (const auto &v : zero_syns[m].RawSpan()) {
      EXPECT_TRUE(gf16->IsZero(v));
    }
  }

  // Verify ToLinearCode produces flattened [45, 21, 9] code matching row-major
  // cw_mat.
  const auto flat_code = irs.ToLinearCode();
  EXPECT_EQ(flat_code.Length(), 45U);
  EXPECT_EQ(flat_code.Dimension(), 21U);
  EXPECT_EQ(flat_code.MinimumDistance(), 9U);
  {
    using RawArray = typename GFArray<GF16>::raw_array_type;
    using ShapeType = typename GFArray<GF16>::shape_type;
    RawArray flat_cw_raw = RawArray::from_shape(ShapeType{45});
    const auto span = cw_mat.RawSpan();
    std::copy(span.begin(), span.end(), flat_cw_raw.data());
    auto flat_cw = GFArray<GF16>::FromRaw(gf16, std::move(flat_cw_raw));
    EXPECT_TRUE(flat_code.Contains(flat_cw));
  }

  // Zero-error fast path.
  {
    const auto res = irs.Decode(cw_mat, IRSDecoderMode::kSsb, false);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword_matrix, cw_mat);
    EXPECT_EQ(res.messages, msgs);
    EXPECT_TRUE(res.burst_error_positions.empty());
    EXPECT_EQ(res.num_burst_errors, 0U);
    EXPECT_EQ(res.num_erasures, 0U);
    EXPECT_FALSE(res.used_fallback);
  }

  // Inject t = 5 and t = 6 burst errors (strictly beyond t_u = 4):
  // verify kSsb, kSsbExact, and kBkyHankel ALL recover exact codewords while
  // kIndependent fails!
  const std::vector<IRSDecoderMode> collab_modes = {
      IRSDecoderMode::kSsb,
      IRSDecoderMode::kSsbExact,
      IRSDecoderMode::kBkyHankel,
  };

  for (std::size_t t : {std::size_t{5}, std::size_t{6}}) {
    std::vector<std::size_t> pos(15);
    std::iota(pos.begin(), pos.end(), std::size_t{0});
    std::size_t verified_trials = 0;

    for (std::size_t trial = 0; trial < 12; ++trial) {
      std::shuffle(pos.begin(), pos.end(), rng);
      std::vector<std::size_t> burst(
          pos.begin(), pos.begin() + static_cast<std::ptrdiff_t>(t));
      std::sort(burst.begin(), burst.end());

      const auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
      EXPECT_FALSE(irs.Contains(rx));

      // Verify independent BM fails (since every column has t > 4 errors).
      const auto ind_res =
          irs.Decode(rx, IRSDecoderMode::kIndependent, /*fallback=*/false);
      EXPECT_FALSE(ind_res.Ok());

      // Check BKY first; if the random draw over GF(16) at t = 6 happens to be
      // full-rank (happens > 93% of trials), ALL 3 collaborative modes MUST
      // agree and recover the exact codeword!
      const auto bky_res =
          irs.Decode(rx, IRSDecoderMode::kBkyHankel, /*fallback=*/false);
      if (!bky_res.Ok()) {
        // Rare (~6% at t = 6 over GF(16)) random rank-deficient draw.
        continue;
      }

      ++verified_trials;
      for (IRSDecoderMode mode : collab_modes) {
        const auto res =
            irs.Decode(rx, mode, /*fallback_to_independent=*/false);
        ASSERT_TRUE(res.Ok())
            << "Mode " << static_cast<int>(mode) << " failed at t=" << t
            << " trial=" << trial << ": " << res.error_message;
        EXPECT_EQ(res.codeword_matrix, cw_mat);
        EXPECT_EQ(res.messages, msgs);
        EXPECT_EQ(res.burst_error_positions, burst);
        EXPECT_EQ(res.num_burst_errors, t);
        EXPECT_EQ(res.num_erasures, 0U);
        EXPECT_FALSE(res.used_fallback);
        EXPECT_LT(res.failure_probability_bound, 1.0);
      }
    }
    EXPECT_GE(verified_trials, 8U)
        << "Expected high collaborative decoding success rate at t=" << t;
  }
}

TEST(InterleavedRSTest, ExtendedGRSWithZeroEvaluationPoint) {
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);

  // Extended IRS(n = 16, k = 8; ell = 3) over GF(16): alpha_0 = 0, r = 8,
  // t_u = 4, t_max = 6.
  auto irs =
      InterleavedGRSCode<GF16>::CreateExtendedHomogeneous(gf16, 16, 8, 3);
  ASSERT_TRUE(irs.ZeroEvaluationIndex().has_value());
  EXPECT_EQ(*irs.ZeroEvaluationIndex(), 0U);
  EXPECT_EQ(irs.UniqueDecodingRadius(), 4U);
  EXPECT_EQ(irs.CollaborativeDecodingRadius(), 6U);

  std::mt19937_64 rng(77777ULL);
  auto msg_mat = GFMatrix<GF16>::Random(gf16, {8, 3}, rng);
  const auto cw_mat = irs.EncodeHomogeneous(msg_mat);
  ASSERT_TRUE(irs.Contains(cw_mat));

  // Case 1: Single burst error hitting ONLY alpha_{j_0} = 0.
  {
    const std::vector<std::size_t> burst = {0};
    const auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res = irs.Decode(rx, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 1U);
    }
  }

  // Case 2: Burst error hitting alpha_{j_0} = 0 PLUS 5 non-zero positions
  // (6 total burst errors = t_max, strictly exceeding t_u = 4).
  {
    std::vector<std::size_t> nonzero_pos(15);
    std::iota(nonzero_pos.begin(), nonzero_pos.end(), std::size_t{1});
    std::size_t successes = 0;

    for (std::size_t trial = 0; trial < 10; ++trial) {
      std::shuffle(nonzero_pos.begin(), nonzero_pos.end(), rng);
      std::vector<std::size_t> burst = {0};
      burst.insert(burst.end(), nonzero_pos.begin(), nonzero_pos.begin() + 5);
      std::sort(burst.begin(), burst.end());

      const auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
      EXPECT_FALSE(irs.Decode(rx, IRSDecoderMode::kIndependent, false).Ok());

      const auto bky = irs.Decode(rx, IRSDecoderMode::kBkyHankel, false);
      if (!bky.Ok()) {
        continue;
      }
      ++successes;
      for (IRSDecoderMode mode :
           {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
            IRSDecoderMode::kBkyHankel}) {
        const auto res = irs.Decode(rx, mode, /*fallback=*/false);
        ASSERT_TRUE(res.Ok()) << "Extended alpha=0 pass 2 failed for mode "
                              << static_cast<int>(mode);
        EXPECT_EQ(res.codeword_matrix, cw_mat);
        EXPECT_EQ(res.burst_error_positions, burst);
        EXPECT_EQ(res.num_burst_errors, 6U);
        EXPECT_FALSE(res.used_fallback);
      }
    }
    EXPECT_GE(successes, 8U);
  }

  // Case 3: alpha_{j_0} = 0 placed in the middle (index 7) with non-trivial
  // column multipliers v over odd prime field GF(17).
  {
    using GF17 = GaloisFieldPrime<uint16_t>;
    auto gf17 = std::make_shared<GF17>(17);
    auto alphas = MakeGFArrayFromInts(
        gf17, {1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 12, 13, 14, 15});
    auto vs = MakeGFArrayFromInts(
        gf17, {3, 5, 7, 11, 13, 2, 6, 9, 14, 4, 8, 12, 15, 16, 1, 10});
    auto irs17 =
        InterleavedGRSCode<GF17>::CreateHomogeneous(gf17, alphas, vs, 8, 3);
    ASSERT_TRUE(irs17.ZeroEvaluationIndex().has_value());
    EXPECT_EQ(*irs17.ZeroEvaluationIndex(), 7U);

    auto m17 = GFMatrix<GF17>::Random(gf17, {8, 3}, rng);
    const auto cw17 = irs17.EncodeHomogeneous(m17);
    const std::vector<std::size_t> burst = {1, 4, 7, 10, 13, 15};
    const auto rx17 = InjectFullColumnBurstErrors(cw17, burst, rng);

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res = irs17.Decode(rx17, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword_matrix, cw17);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 6U);
    }
  }
}

TEST(InterleavedRSTest, HeterogeneousIRS15_5_7_9) {
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);

  // Heterogeneous IRS(n = 15; k = (5, 7, 9)) over GF(16):
  // Redundancies r = (10, 8, 6).
  // t_u = floor(6 / 2) = 3.
  // sum_{r_m > 5}(r_m - 5) = 5 + 3 + 1 = 9 >= 5.
  // sum_{r_m > 6}(r_m - 6) = 4 + 2 = 6 >= 6.
  auto irs = InterleavedGRSCode<GF16>::CreateHeterogeneous(gf16, 15, {5, 7, 9});
  EXPECT_FALSE(irs.IsHomogeneous());
  EXPECT_EQ(irs.Dimensions(), (std::vector<std::size_t>{5, 7, 9}));
  EXPECT_EQ(irs.UniqueDecodingRadius(), 3U);
  EXPECT_EQ(irs.CollaborativeDecodingRadius(), 6U);

  // EncodeHomogeneous must reject heterogeneous codes.
  auto dummy = GFMatrix<GF16>::Zeros(gf16, {5, 3});
  EXPECT_THROW(static_cast<void>(irs.EncodeHomogeneous(dummy)),
               std::invalid_argument);

  // Flattened linear code has length 45, dimension 5+7+9 = 21, d_min = 7.
  const auto flat = irs.ToLinearCode();
  EXPECT_EQ(flat.Length(), 45U);
  EXPECT_EQ(flat.Dimension(), 21U);
  EXPECT_EQ(flat.MinimumDistance(), 7U);

  std::mt19937_64 rng(90909ULL);
  std::vector<GFArray<GF16>> msgs = {
      GFArray<GF16>::Random(gf16, {5}, rng),
      GFArray<GF16>::Random(gf16, {7}, rng),
      GFArray<GF16>::Random(gf16, {9}, rng),
  };
  const auto cw_mat = irs.Encode(msgs);
  ASSERT_TRUE(irs.Contains(cw_mat));

  // Verify all 3 collaborative modes at t = 5 (strictly above t_u = 3 and
  // above stream 1's unique radius 4).
  std::vector<std::size_t> pos(15);
  std::iota(pos.begin(), pos.end(), std::size_t{0});
  for (std::size_t trial = 0; trial < 8; ++trial) {
    std::shuffle(pos.begin(), pos.end(), rng);
    std::vector<std::size_t> burst(pos.begin(), pos.begin() + 5);
    std::sort(burst.begin(), burst.end());

    const auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
    EXPECT_FALSE(irs.Decode(rx, IRSDecoderMode::kIndependent, false).Ok());

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res = irs.Decode(rx, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok()) << "Heterogeneous mode " << static_cast<int>(mode)
                            << " failed on trial " << trial;
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.messages, msgs);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 5U);
      EXPECT_FALSE(res.used_fallback);
    }
  }
}

TEST(InterleavedRSTest, ColumnErasuresAndBurstErrorsCombined) {
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);
  std::mt19937_64 rng(314159ULL);

  // 1. Primitive IRS(15, 7; 3): r = 8. With s = 2 erasures, r - s = 6:
  //    unique error budget = floor(6/2) = 3,
  //    collaborative error budget = floor(3/4 * 6) = 4!
  {
    auto irs = InterleavedGRSCode<GF16>::CreateHomogeneous(gf16, 15, 7, 3);
    auto msg_mat = GFMatrix<GF16>::Random(gf16, {7, 3}, rng);
    const auto cw_mat = irs.EncodeHomogeneous(msg_mat);

    const std::vector<std::size_t> erasures = {2, 11};
    const std::vector<std::size_t> burst = {0, 5, 8, 14};

    // Corrupt BOTH the burst positions and the erased positions (erased values
    // must be completely ignored by the decoder).
    auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
    rx = InjectFullColumnBurstErrors(rx, erasures, rng);

    // Independent BM fails since 2*4 + 2 = 10 > 8 = r.
    EXPECT_FALSE(irs.DecodeWithErasures(rx, erasures,
                                        IRSDecoderMode::kIndependent, false)
                     .Ok());

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res =
          irs.DecodeWithErasures(rx, erasures, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok())
          << "Erasures+errors failed for mode " << static_cast<int>(mode);
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 4U);
      EXPECT_EQ(res.num_erasures, 2U);
      EXPECT_FALSE(res.used_fallback);
    }

    // Invalid erasure checks.
    EXPECT_FALSE(
        irs.DecodeWithErasures(rx, {15}, IRSDecoderMode::kSsb, false).Ok());
    EXPECT_FALSE(
        irs.DecodeWithErasures(rx, {3, 3}, IRSDecoderMode::kSsb, false).Ok());
    EXPECT_FALSE(irs.DecodeWithErasures(rx, {0, 1, 2, 3, 4, 5, 6, 7, 8},
                                        IRSDecoderMode::kSsb, false)
                     .Ok());
  }

  // 2. Extended IRS(16, 8; 3) with alpha_0 = 0 ERASED + 1 non-zero erasure +
  //    4 burst errors (tests Case A shifted path!).
  {
    auto irs =
        InterleavedGRSCode<GF16>::CreateExtendedHomogeneous(gf16, 16, 8, 3);
    auto msg_mat = GFMatrix<GF16>::Random(gf16, {8, 3}, rng);
    const auto cw_mat = irs.EncodeHomogeneous(msg_mat);

    const std::vector<std::size_t> erasures = {0, 9};
    const std::vector<std::size_t> burst = {2, 6, 12, 15};
    auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
    rx = InjectFullColumnBurstErrors(rx, erasures, rng);

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res =
          irs.DecodeWithErasures(rx, erasures, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok()) << "Extended alpha_0 erased failed for mode "
                            << static_cast<int>(mode);
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 4U);
      EXPECT_EQ(res.num_erasures, 2U);
    }
  }

  // 3. Extended IRS(16, 8; 3) with 2 non-zero erasures + burst error hitting
  //    alpha_0 = 0 plus 3 non-zero burst positions (4 total burst errors,
  //    tests Case C Pass 2 with erasures!).
  {
    auto irs =
        InterleavedGRSCode<GF16>::CreateExtendedHomogeneous(gf16, 16, 8, 3);
    auto msg_mat = GFMatrix<GF16>::Random(gf16, {8, 3}, rng);
    const auto cw_mat = irs.EncodeHomogeneous(msg_mat);

    const std::vector<std::size_t> erasures = {3, 10};
    const std::vector<std::size_t> burst = {0, 5, 11, 14};
    auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);
    rx = InjectFullColumnBurstErrors(rx, erasures, rng);

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res =
          irs.DecodeWithErasures(rx, erasures, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok())
          << "Extended alpha_0 error + erasures failed for mode "
          << static_cast<int>(mode);
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 4U);
      EXPECT_EQ(res.num_erasures, 2U);
    }
  }
}

TEST(InterleavedRSTest, MultiSequenceLfsrEdgeCaseAndRankDeficientFallback) {
  // 1. Regression test for the multi-sequence LFSR edge case:
  //    S^(0) = (1, 0, 1, 0), S^(1) = (1, 1, 1, 1).
  //    Unique minimal LFSR over ANY field is Lambda(z) = 1 - z^2 (L = 2).
  //    Without the exact early-break invariant on !valid[m], SSB erroneously
  //    produces L = 3!
  {
    // Over GF(16) (characteristic 2: 1 - z^2 == 1 + z^2).
    using GF16 = GaloisFieldBinaryExtension<uint8_t>;
    auto gf16 = std::make_shared<GF16>(4);
    const uint8_t z = gf16->AdditiveIdentity();
    const uint8_t o = gf16->MultiplicativeIdentity();
    const std::vector<std::vector<uint8_t>> seqs = {
        {o, z, o, z},
        {o, o, o, o},
    };

    const auto ssb = xg::coding::detail::SolveMultiSequenceSsb(gf16, seqs, 3);
    EXPECT_EQ(ssb.length, 2U);
    EXPECT_EQ(ssb.locator.Degree(), 2);
    EXPECT_EQ(ssb.locator.RawCoeff(0), o);
    EXPECT_EQ(ssb.locator.RawCoeff(1), z);
    EXPECT_EQ(ssb.locator.RawCoeff(2), o);

    const auto popov =
        xg::coding::detail::SolveMultiSequencePopov(gf16, seqs, 3);
    ASSERT_TRUE(popov.has_value());
    EXPECT_EQ(popov->Degree(), 2);
    EXPECT_EQ(popov->RawCoeff(0), o);
    EXPECT_EQ(popov->RawCoeff(1), z);
    EXPECT_EQ(popov->RawCoeff(2), o);

    const auto bky =
        xg::coding::detail::SolveMultiSequenceBkyHankel(gf16, seqs, 3);
    ASSERT_TRUE(bky.has_value());
    EXPECT_EQ(bky->Degree(), 2);
    EXPECT_EQ(bky->RawCoeff(0), o);
    EXPECT_EQ(bky->RawCoeff(1), z);
    EXPECT_EQ(bky->RawCoeff(2), o);
  }

  {
    // Over GF(17) (odd characteristic: 1 - z^2 has coeff[2] = -1 = 16).
    using GF17 = GaloisFieldPrime<uint16_t>;
    auto gf17 = std::make_shared<GF17>(17);
    const uint16_t z = gf17->AdditiveIdentity();
    const uint16_t o = gf17->MultiplicativeIdentity();
    const uint16_t neg_o = gf17->Neg(o);
    const std::vector<std::vector<uint16_t>> seqs = {
        {o, z, o, z},
        {o, o, o, o},
    };

    const auto ssb = xg::coding::detail::SolveMultiSequenceSsb(gf17, seqs, 3);
    EXPECT_EQ(ssb.length, 2U);
    EXPECT_EQ(ssb.locator.Degree(), 2);
    EXPECT_EQ(ssb.locator.RawCoeff(0), o);
    EXPECT_EQ(ssb.locator.RawCoeff(1), z);
    EXPECT_EQ(ssb.locator.RawCoeff(2), neg_o);

    const auto popov =
        xg::coding::detail::SolveMultiSequencePopov(gf17, seqs, 3);
    ASSERT_TRUE(popov.has_value());
    EXPECT_EQ(popov->Degree(), 2);
    EXPECT_EQ(popov->RawCoeff(0), o);
    EXPECT_EQ(popov->RawCoeff(1), z);
    EXPECT_EQ(popov->RawCoeff(2), neg_o);

    const auto bky =
        xg::coding::detail::SolveMultiSequenceBkyHankel(gf17, seqs, 3);
    ASSERT_TRUE(bky.has_value());
    EXPECT_EQ(bky->Degree(), 2);
    EXPECT_EQ(bky->RawCoeff(0), o);
    EXPECT_EQ(bky->RawCoeff(1), z);
    EXPECT_EQ(bky->RawCoeff(2), neg_o);
  }

  // 2. Fallback test (used_fallback == true):
  //    IRS(15, 7; 3) over GF(16): t_u = 4, t_max = 6.
  //    (a) Rank-deficient t = 6 pattern: col 0 & col 1 share identical 4 errors
  //        at {0, 1, 2, 3}, col 2 has 2 errors at {4, 5}.
  //        rank(M_6) <= 2 + 0 + 2 = 4 < 6, while every column has <= 4 errors!
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);
  auto irs = InterleavedGRSCode<GF16>::CreateHomogeneous(gf16, 15, 7, 3);
  std::mt19937_64 rng(55555ULL);
  auto msg_mat = GFMatrix<GF16>::Random(gf16, {7, 3}, rng);
  const auto cw_mat = irs.EncodeHomogeneous(msg_mat);

  {
    auto rx = cw_mat;
    for (std::size_t row :
         {std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
      const auto e = RandomNonZero(gf16, rng);
      rx.Raw(row, 0) = gf16->Add(rx.Raw(row, 0), e);
      rx.Raw(row, 1) = gf16->Add(rx.Raw(row, 1), e);
    }
    for (std::size_t row : {std::size_t{4}, std::size_t{5}}) {
      const auto e = RandomNonZero(gf16, rng);
      rx.Raw(row, 2) = gf16->Add(rx.Raw(row, 2), e);
    }

    // Without fallback, BKY and exact Popov must detect rank deficiency and
    // decline to guess; with fallback, they recover exact codeword via BM!
    for (IRSDecoderMode mode :
         {IRSDecoderMode::kSsbExact, IRSDecoderMode::kBkyHankel}) {
      const auto no_fb = irs.Decode(rx, mode, /*fallback=*/false);
      EXPECT_FALSE(no_fb.Ok());

      const auto with_fb = irs.Decode(rx, mode, /*fallback=*/true);
      ASSERT_TRUE(with_fb.Ok());
      EXPECT_TRUE(with_fb.used_fallback);
      EXPECT_EQ(with_fb.codeword_matrix, cw_mat);
      EXPECT_EQ(with_fb.num_burst_errors, 6U);
    }
  }

  //    (b) Spread-across-columns pattern (total burst weight 7 > t_max = 6,
  //        each column has 3 <= t_u = 4 errors): ALL 3 collaborative modes
  //        must fail without fallback and succeed with used_fallback == true!
  {
    auto rx = cw_mat;
    for (std::size_t row : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
      rx.Raw(row, 0) = gf16->Add(rx.Raw(row, 0), RandomNonZero(gf16, rng));
    }
    for (std::size_t row : {std::size_t{2}, std::size_t{3}, std::size_t{4}}) {
      rx.Raw(row, 1) = gf16->Add(rx.Raw(row, 1), RandomNonZero(gf16, rng));
    }
    for (std::size_t row : {std::size_t{4}, std::size_t{5}, std::size_t{6}}) {
      rx.Raw(row, 2) = gf16->Add(rx.Raw(row, 2), RandomNonZero(gf16, rng));
    }

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      EXPECT_FALSE(irs.Decode(rx, mode, /*fallback=*/false).Ok());

      const auto with_fb = irs.Decode(rx, mode, /*fallback=*/true);
      ASSERT_TRUE(with_fb.Ok());
      EXPECT_TRUE(with_fb.used_fallback);
      EXPECT_EQ(with_fb.codeword_matrix, cw_mat);
      EXPECT_EQ(with_fb.num_burst_errors, 7U);
      EXPECT_EQ(with_fb.burst_error_positions,
                (std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6}));
    }
  }

  //    (c) Rank-1 uncorrectable t = 5 pattern (identical errors across all 3
  //        columns): all 3 collaborative modes detect rank deficiency and fail
  //        cleanly (no silent miscorrection).
  {
    auto rx = cw_mat;
    for (std::size_t row : {std::size_t{0}, std::size_t{2}, std::size_t{5},
                            std::size_t{9}, std::size_t{13}}) {
      const auto e = RandomNonZero(gf16, rng);
      for (std::size_t m = 0; m < 3; ++m) {
        rx.Raw(row, m) = gf16->Add(rx.Raw(row, m), e);
      }
    }
    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      EXPECT_FALSE(irs.Decode(rx, mode, /*fallback=*/true).Ok());
    }
  }
}

TEST(InterleavedRSTest, FailureProbabilityBoundSanityChecks) {
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto gf64 = std::make_shared<GF2XLOG<uint8_t>>(6);
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);

  auto irs16 = InterleavedGRSCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateHomogeneous(gf16, 15, 7, 3);
  auto irs64 =
      InterleavedGRSCode<GF2XLOG<uint8_t>>::CreateHomogeneous(gf64, 15, 7, 3);
  auto irs256 =
      InterleavedGRSCode<GF2XLOG<uint8_t>>::CreateHomogeneous(gf256, 15, 7, 3);

  // Within unique decoding radius (t <= 4): exact 0.0.
  for (std::size_t t = 0; t <= 4; ++t) {
    EXPECT_DOUBLE_EQ(irs16.FailureProbabilityBound(t), 0.0);
    EXPECT_DOUBLE_EQ(irs64.FailureProbabilityBound(t), 0.0);
    EXPECT_DOUBLE_EQ(irs256.FailureProbabilityBound(t), 0.0);
  }

  // Strictly below and at t_max = 6: strictly in (0, 1).
  const double p16_5 = irs16.FailureProbabilityBound(5);
  const double p16_6 = irs16.FailureProbabilityBound(6);
  const double p64_5 = irs64.FailureProbabilityBound(5);
  const double p64_6 = irs64.FailureProbabilityBound(6);
  const double p256_5 = irs256.FailureProbabilityBound(5);
  const double p256_6 = irs256.FailureProbabilityBound(6);

  EXPECT_GT(p16_5, 0.0);
  EXPECT_LT(p16_5, 1.0);
  EXPECT_GT(p16_6, p16_5);
  EXPECT_LT(p16_6, 1.0);

  // Monotonically decreases with field size q.
  EXPECT_GT(p16_5, p64_5);
  EXPECT_GT(p64_5, p256_5);
  EXPECT_GT(p16_6, p64_6);
  EXPECT_GT(p64_6, p256_6);

  // Beyond collaborative decoding radius (t > 6): 1.0.
  EXPECT_DOUBLE_EQ(irs16.FailureProbabilityBound(7), 1.0);
  EXPECT_DOUBLE_EQ(irs256.FailureProbabilityBound(8), 1.0);

  // Heterogeneous bound sanity check.
  auto het = InterleavedGRSCode<
      GaloisFieldBinaryExtension<uint8_t>>::CreateHeterogeneous(gf16, 15,
                                                                {5, 7, 9});
  EXPECT_DOUBLE_EQ(het.FailureProbabilityBound(3), 0.0);
  EXPECT_GT(het.FailureProbabilityBound(4), 0.0);
  EXPECT_LT(het.FailureProbabilityBound(4), 1.0);
  EXPECT_GT(het.FailureProbabilityBound(5), het.FailureProbabilityBound(4));
  EXPECT_LT(het.FailureProbabilityBound(5), 1.0);
  EXPECT_DOUBLE_EQ(het.FailureProbabilityBound(7), 1.0);
}

TEST(InterleavedRSTest, PowerDecodingRS63_10AndGRS63_10BeyondHalfDistance) {
  using GF64 = GF2XLOG<uint8_t>;
  auto gf64 = std::make_shared<GF64>(6);

  // Verify closed-form radius helpers on [63, 10]:
  // t_u = 26, ell_p = 2 -> t_power = 32, ell_p = 3 -> t_power = 33.
  EXPECT_EQ(PowerDecodingDecoder<GF64>::PowerDecodingRadius(63, 10, 1), 26U);
  EXPECT_EQ(PowerDecodingDecoder<GF64>::PowerDecodingRadius(63, 10, 2), 32U);
  EXPECT_EQ(PowerDecodingDecoder<GF64>::PowerDecodingRadius(63, 10, 3), 33U);
  EXPECT_EQ(PowerDecodingDecoder<GF64>::OptimalPoweringDegree(63, 10), 3U);

  std::mt19937_64 rng(63102026ULL);

  // 1. Standard RS(63, 10) over GF(64) with ell_p = 2 (t_power = 32) and
  //    auto optimal ell_p = 3 (t_power = 33).
  {
    auto rs = GeneralizedReedSolomonCode<GF64>::ReedSolomon(gf64, 63, 10);
    BerlekampMasseyDecoder<GF64> bm(rs);
    EXPECT_EQ(bm.DecodingRadius(), 26U);

    PowerDecodingDecoder<GF64> pd_ell2_ssb(rs, 2, IRSDecoderMode::kSsb, false);
    PowerDecodingDecoder<GF64> pd_ell2_exact(rs, 2, IRSDecoderMode::kSsbExact,
                                             false);
    PowerDecodingDecoder<GF64> pd_ell2_bky(rs, 2, IRSDecoderMode::kBkyHankel,
                                           false);
    PowerDecodingDecoder<GF64> pd_auto(rs);

    EXPECT_EQ(pd_ell2_ssb.Name(), "PowerDecoding");
    EXPECT_EQ(pd_ell2_ssb.PoweringDegree(), 2U);
    EXPECT_EQ(pd_ell2_ssb.DecodingRadius(), 32U);
    EXPECT_EQ(pd_auto.PoweringDegree(), 3U);
    EXPECT_EQ(pd_auto.DecodingRadius(), 33U);

    const auto msg = GFArray<GF64>::Random(gf64, {10}, rng);
    const auto cw = rs.Encode(msg);
    ASSERT_TRUE(rs.Contains(cw));

    std::vector<std::size_t> indices(63);
    std::iota(indices.begin(), indices.end(), std::size_t{0});

    for (std::size_t t :
         {std::size_t{29}, std::size_t{30}, std::size_t{31}, std::size_t{32}}) {
      std::shuffle(indices.begin(), indices.end(), rng);
      std::vector<std::size_t> err_pos(
          indices.begin(), indices.begin() + static_cast<std::ptrdiff_t>(t));
      const auto rx = CorruptWordAtPositions(cw, err_pos, rng);

      // Unique decoder strictly fails beyond t_u = 26.
      EXPECT_FALSE(bm.DecodeWithResult(rx).Ok());

      for (const PowerDecodingDecoder<GF64> *dec :
           {&pd_ell2_ssb, &pd_ell2_exact, &pd_ell2_bky, &pd_auto}) {
        const auto res = dec->DecodeWithResult(rx);
        ASSERT_TRUE(res.Ok()) << "PowerDecoding RS(63,10) failed at t=" << t
                              << " ell_p=" << dec->PoweringDegree()
                              << " mode=" << static_cast<int>(dec->Mode())
                              << ": " << res.error_message;
        EXPECT_EQ(res.codeword, cw);
        EXPECT_EQ(res.message, msg);
        EXPECT_EQ(res.num_errors, t);
        EXPECT_EQ(res.num_erasures, 0U);
      }
    }
  }

  // 2. GRS(63, 10) over GF(64) with non-trivial column multipliers v_j != 1.
  {
    auto evals = xg::coding::detail::BuildPrimitiveEvalPoints(gf64, 63);
    using RawArray = typename GFArray<GF64>::raw_array_type;
    using ShapeType = typename GFArray<GF64>::shape_type;
    RawArray v_raw = RawArray::from_shape(ShapeType{63});
    for (std::size_t j = 0; j < 63; ++j) {
      v_raw(j) = RandomNonZero(gf64, rng);
    }
    auto vs = GFArray<GF64>::FromRaw(gf64, std::move(v_raw));

    GeneralizedReedSolomonCode<GF64> grs(gf64, evals, vs, 10);
    BerlekampMasseyDecoder<GF64> bm(grs);
    PowerDecodingDecoder<GF64> pd_ssb(grs, 2, IRSDecoderMode::kSsb, false);
    PowerDecodingDecoder<GF64> pd_exact(grs, 2, IRSDecoderMode::kSsbExact,
                                        false);
    PowerDecodingDecoder<GF64> pd_bky(grs, 2, IRSDecoderMode::kBkyHankel,
                                      false);

    const auto msg = GFArray<GF64>::Random(gf64, {10}, rng);
    const auto cw = grs.Encode(msg);
    ASSERT_TRUE(grs.Contains(cw));

    std::vector<std::size_t> indices(63);
    std::iota(indices.begin(), indices.end(), std::size_t{0});

    for (std::size_t t :
         {std::size_t{29}, std::size_t{30}, std::size_t{31}, std::size_t{32}}) {
      std::shuffle(indices.begin(), indices.end(), rng);
      std::vector<std::size_t> err_pos(
          indices.begin(), indices.begin() + static_cast<std::ptrdiff_t>(t));
      const auto rx = CorruptWordAtPositions(cw, err_pos, rng);

      EXPECT_FALSE(bm.DecodeWithResult(rx).Ok());

      for (const PowerDecodingDecoder<GF64> *dec :
           {&pd_ssb, &pd_exact, &pd_bky}) {
        const auto res = dec->DecodeWithResult(rx);
        ASSERT_TRUE(res.Ok()) << "PowerDecoding GRS(63,10) failed at t=" << t
                              << " mode=" << static_cast<int>(dec->Mode());
        EXPECT_EQ(res.codeword, cw);
        EXPECT_EQ(res.message, msg);
        EXPECT_EQ(res.num_errors, t);
        EXPECT_EQ(res.num_erasures, 0U);
      }
    }
  }

  // 3. Extended GRS(64, 10) over GF(64) with alpha_0 = 0, non-trivial v,
  //    and an error hitting alpha_0 = 0 plus 30 non-zero positions (31 > 27).
  {
    auto evals = xg::coding::detail::BuildExtendedEvalPoints(gf64, 64);
    using RawArray = typename GFArray<GF64>::raw_array_type;
    using ShapeType = typename GFArray<GF64>::shape_type;
    RawArray v_raw = RawArray::from_shape(ShapeType{64});
    for (std::size_t j = 0; j < 64; ++j) {
      v_raw(j) = RandomNonZero(gf64, rng);
    }
    auto vs = GFArray<GF64>::FromRaw(gf64, std::move(v_raw));

    GeneralizedReedSolomonCode<GF64> ext_grs(gf64, evals, vs, 10);
    PowerDecodingDecoder<GF64> pd(ext_grs, 2, IRSDecoderMode::kSsb, false);

    const auto msg = GFArray<GF64>::Random(gf64, {10}, rng);
    const auto cw = ext_grs.Encode(msg);

    std::vector<std::size_t> nz_indices(63);
    std::iota(nz_indices.begin(), nz_indices.end(), std::size_t{1});
    std::shuffle(nz_indices.begin(), nz_indices.end(), rng);
    std::vector<std::size_t> err_pos = {0};
    err_pos.insert(err_pos.end(), nz_indices.begin(), nz_indices.begin() + 30);

    const auto rx = CorruptWordAtPositions(cw, err_pos, rng);
    const auto res = pd.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 31U);
  }

  // 4. Power Decoding with erasures: RS(63, 10), s = 4 erasures + t = 28
  //    random errors (2*28 + 4 = 60 > 53 = r -> BM fails, PowerDecoding wins!).
  {
    auto rs = GeneralizedReedSolomonCode<GF64>::ReedSolomon(gf64, 63, 10);
    BerlekampMasseyDecoder<GF64> bm(rs);
    PowerDecodingDecoder<GF64> pd(rs, 2, IRSDecoderMode::kSsb, false);

    const auto msg = GFArray<GF64>::Random(gf64, {10}, rng);
    const auto cw = rs.Encode(msg);

    std::vector<std::size_t> indices(63);
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    std::shuffle(indices.begin(), indices.end(), rng);
    std::vector<std::size_t> erasures(indices.begin(), indices.begin() + 4);
    std::vector<std::size_t> err_pos(indices.begin() + 4, indices.begin() + 32);

    auto rx = CorruptWordAtPositions(cw, err_pos, rng);
    rx = CorruptWordAtPositions(rx, erasures, rng);

    EXPECT_FALSE(bm.DecodeWithErasures(rx, erasures).Ok());

    const auto res = pd.DecodeWithErasures(rx, erasures);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 28U);
    EXPECT_EQ(res.num_erasures, 4U);
  }
}

TEST(InterleavedRSTest, MultiFieldZechAndOddExtensionSupport) {
  std::mt19937_64 rng(88888ULL);

  // GF(2^4) Zech log representation (internal value != integer value).
  {
    using GFZech = GF2XZECH<uint8_t>;
    auto gf16_zech = std::make_shared<GFZech>(4);
    auto irs =
        InterleavedGRSCode<GFZech>::CreateHomogeneous(gf16_zech, 15, 7, 3);
    auto msg_mat = GFMatrix<GFZech>::Random(gf16_zech, {7, 3}, rng);
    const auto cw_mat = irs.EncodeHomogeneous(msg_mat);
    const std::vector<std::size_t> burst = {1, 4, 7, 10, 14};
    const auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res = irs.Decode(rx, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.burst_error_positions, burst);
    }
  }

  // GF(3^3) odd-characteristic extension field (non-integral ExtensionElement).
  {
    using GF27 = GaloisFieldExtension<uint8_t>;
    auto gf27 = std::make_shared<GF27>(std::pair<uint64_t, uint64_t>{3, 3});
    auto irs = InterleavedGRSCode<GF27>::CreateHomogeneous(gf27, 20, 10, 3);
    EXPECT_EQ(irs.UniqueDecodingRadius(), 5U);
    EXPECT_EQ(irs.CollaborativeDecodingRadius(), 7U);

    auto msg_mat = GFMatrix<GF27>::Random(gf27, {10, 3}, rng);
    const auto cw_mat = irs.EncodeHomogeneous(msg_mat);
    const std::vector<std::size_t> burst = {0, 3, 6, 9, 12, 15, 19};
    const auto rx = InjectFullColumnBurstErrors(cw_mat, burst, rng);

    for (IRSDecoderMode mode : {IRSDecoderMode::kSsb, IRSDecoderMode::kSsbExact,
                                IRSDecoderMode::kBkyHankel}) {
      const auto res = irs.Decode(rx, mode, /*fallback=*/false);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword_matrix, cw_mat);
      EXPECT_EQ(res.burst_error_positions, burst);
      EXPECT_EQ(res.num_burst_errors, 7U);
    }
  }
}

TEST(InterleavedRSTest, ValidationAndErrorHandling) {
  using GF16 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf16 = std::make_shared<GF16>(4);

  EXPECT_THROW((InterleavedGRSCode<GF16>::CreateHomogeneous(nullptr, 15, 7, 3)),
               std::invalid_argument);
  EXPECT_THROW((InterleavedGRSCode<GF16>::CreateHomogeneous(gf16, 15, 7, 0)),
               std::invalid_argument);

  auto rs = GeneralizedReedSolomonCode<GF16>::ReedSolomon(gf16, 15, 7);
  EXPECT_THROW((PowerDecodingDecoder<GF16>(rs, 0)), std::invalid_argument);
  // ell_p = 3 violates 3*(7-1)+1 = 19 >= 15.
  EXPECT_THROW((PowerDecodingDecoder<GF16>(rs, 3)), std::invalid_argument);
  EXPECT_THROW(
      (PowerDecodingDecoder<GF16>(rs, 2, IRSDecoderMode::kIndependent)),
      std::invalid_argument);

  // PowerDecodingDecoder requires a GRS code.
  auto ham = HammingCode(3);
  EXPECT_THROW((PowerDecodingDecoder<GaloisFieldBinary>(ham)),
               std::invalid_argument);
}

}  // namespace
