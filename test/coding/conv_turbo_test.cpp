#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/convolutional.hpp"
#include "xgalois/coding/turbo.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace {

using namespace xg;
using namespace xg::coding;
using namespace xg::channel;

bool ArraysEqual(const GFArray<GaloisFieldBinary> &a,
                 const GFArray<GaloisFieldBinary> &b) {
  if (!a.HasField() || !b.HasField() || a.Size() != b.Size() ||
      a.Dimension() != b.Dimension()) {
    return false;
  }
  const auto sa = a.RawSpan();
  const auto sb = b.RawSpan();
  for (std::size_t i = 0; i < sa.size(); ++i) {
    if ((sa[i] & 1U) != (sb[i] & 1U)) {
      return false;
    }
  }
  return true;
}

std::size_t CountBitDiffs(const GFArray<GaloisFieldBinary> &a,
                          const GFArray<GaloisFieldBinary> &b) {
  const auto sa = a.RawSpan();
  const auto sb = b.RawSpan();
  std::size_t diffs = 0;
  for (std::size_t i = 0; i < sa.size(); ++i) {
    if ((sa[i] & 1U) != (sb[i] & 1U)) {
      ++diffs;
    }
  }
  return diffs;
}

double LogSumExp(const std::vector<double> &vals) {
  double m = -std::numeric_limits<double>::infinity();
  for (double v : vals) {
    if (v > m) {
      m = v;
    }
  }
  double sum = 0.0;
  for (double v : vals) {
    sum += std::exp(v - m);
  }
  return m + std::log(sum);
}

// ============================================================================
// 1. NASA (2, 1, 7) {0133, 0171} Feedforward Convolutional Code & Viterbi
// ============================================================================

TEST(ConvolutionalCodeTest, NasaRateHalfStructureAndLinearity) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 24;
  auto code = ConvolutionalCode::NasaRateHalf(gf2, k);

  EXPECT_EQ(code.Dimension(), k);
  EXPECT_EQ(code.Memory(), 6U);
  EXPECT_EQ(code.NumOutputBits(), 2U);
  EXPECT_EQ(code.GetTrellis().num_states, 64U);
  EXPECT_EQ(code.Length(), (k + 6) * 2);
  EXPECT_FALSE(code.IsRecursiveSystematic());

  std::mt19937_64 rng(101);
  auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
  auto cw = code.Encode(msg);
  EXPECT_EQ(cw.Size(), code.Length());

  // Lazy GeneratorMatrix and ParityCheckMatrix consistency.
  const auto &g = code.GeneratorMatrix();
  const auto &h = code.ParityCheckMatrix();
  EXPECT_EQ(g.Shape(0), k);
  EXPECT_EQ(g.Shape(1), code.Length());
  EXPECT_EQ(h.Shape(0), code.Length() - k);
  EXPECT_EQ(h.Shape(1), code.Length());

  auto cw_mat = linalg::Dot(msg, g);
  EXPECT_TRUE(ArraysEqual(cw, cw_mat));
  EXPECT_TRUE(code.Contains(cw));

  auto unenc = code.Unencode(cw);
  EXPECT_TRUE(ArraysEqual(unenc, msg));
}

TEST(ConvolutionalCodeTest, NasaRateHalfViterbiHardNoiselessAndIsolatedErrors) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 48;
  auto code = ConvolutionalCode::NasaRateHalf(gf2, k);
  ViterbiDecoder viterbi(code);

  std::mt19937_64 rng(202);
  auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
  auto cw = code.Encode(msg);

  // Noiseless via default registered decoder "ViterbiHard".
  auto res0 = code.DecodeWithResult(cw);
  ASSERT_TRUE(res0.Ok());
  EXPECT_EQ(res0.num_errors, 0U);
  EXPECT_TRUE(ArraysEqual(res0.message, msg));
  EXPECT_TRUE(ArraysEqual(res0.codeword, cw));

  // NASA K=7 R=1/2 has d_free = 10, so up to 4 arbitrary errors are always
  // uniquely correctable.
  auto corrupted = cw;
  const std::vector<std::size_t> flip_positions = {5, 29, 61, 88};
  auto span = corrupted.RawSpan();
  for (std::size_t pos : flip_positions) {
    span[pos] ^= 1U;
  }

  auto res4 = viterbi.DecodeHard(corrupted);
  ASSERT_TRUE(res4.Ok());
  EXPECT_EQ(res4.num_errors, flip_positions.size());
  EXPECT_TRUE(ArraysEqual(res4.message, msg));
  EXPECT_TRUE(ArraysEqual(res4.codeword, cw));
}

TEST(ConvolutionalCodeTest, NasaRateHalfViterbiSoftAwgnGain) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 100;
  auto code = ConvolutionalCode::NasaRateHalf(gf2, k);
  ViterbiDecoder viterbi(code);

  std::mt19937_64 msg_rng(303);
  BasicAWGNChannel<GaloisFieldBinary> ch =
      BasicAWGNChannel<GaloisFieldBinary>::FromEbN0(gf2, code.Length(), 3.5,
                                                    code.Rate(), 404);

  std::size_t hard_bit_errors = 0;
  std::size_t soft_bit_errors = 0;
  const std::size_t num_frames = 25;
  for (std::size_t f = 0; f < num_frames; ++f) {
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, msg_rng);
    auto cw = code.Encode(msg);
    auto tx = ch.TransmitDetailed(cw);

    auto res_hard = viterbi.DecodeHard(tx.hard_decision);
    auto res_soft = viterbi.DecodeSoft(tx.llr);
    ASSERT_TRUE(res_hard.Ok());
    ASSERT_TRUE(res_soft.Ok());

    hard_bit_errors += CountBitDiffs(res_hard.message, msg);
    soft_bit_errors += CountBitDiffs(res_soft.message, msg);
  }

  // Soft Viterbi should achieve near-zero errors at Eb/N0 = 3.5 dB and beat
  // hard-decision Viterbi.
  EXPECT_LT(soft_bit_errors, hard_bit_errors);
  EXPECT_LE(soft_bit_errors, 5U);
}

TEST(ConvolutionalCodeTest, TruncatedModeEncodeDecode) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 32;
  auto code = ConvolutionalCode::CreateFeedforward(gf2, k, {05U, 07U},
                                                   TerminationMode::kTruncated);
  EXPECT_EQ(code.Length(), 2 * k);

  std::mt19937_64 rng(505);
  auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
  auto cw = code.Encode(msg);

  ViterbiDecoder viterbi(code);
  auto res = viterbi.DecodeHard(cw);
  ASSERT_TRUE(res.Ok());
  EXPECT_EQ(res.num_errors, 0U);
  EXPECT_TRUE(ArraysEqual(res.message, msg));
}

// ============================================================================
// 2. Tail-Biting Convolutional Code
// ============================================================================

TEST(ConvolutionalCodeTest, TailBitingCircularStateAndViterbi) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // Test both feedforward (NASA K=7) and RSC (LTE 8-state) tail-biting codes.
  {
    const std::size_t k = 36;
    auto tb_ff =
        ConvolutionalCode::NasaRateHalf(gf2, k, TerminationMode::kTailBiting);
    EXPECT_EQ(tb_ff.Length(), 2 * k);

    std::mt19937_64 rng(606);
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
    const uint32_t s0 = tb_ff.ComputeTailBitingInitialState(msg.RawSpan());
    uint32_t state = s0;
    const auto &trellis = tb_ff.GetTrellis();
    for (std::size_t t = 0; t < k; ++t) {
      state = trellis.NextState(state, msg.RawSpan()[t] & 1U);
    }
    EXPECT_EQ(state, s0);

    auto cw = tb_ff.Encode(msg);
    EXPECT_TRUE(tb_ff.Contains(cw));

    // Introduce 3 isolated bit errors and verify Viterbi recovers the exact
    // message.
    auto corrupted = cw;
    corrupted.RawSpan()[3] ^= 1U;
    corrupted.RawSpan()[31] ^= 1U;
    corrupted.RawSpan()[64] ^= 1U;

    ViterbiDecoder viterbi(tb_ff);
    auto res = viterbi.DecodeHard(corrupted);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.num_errors, 3U);
    EXPECT_TRUE(ArraysEqual(res.message, msg));
  }

  // RSC tail-biting: LTE feedback 013 has period 7, so K=24 is invertible
  // while K=7 is singular.
  EXPECT_THROW(ConvolutionalCode::LteTurboConstituent(
                   gf2, 7, TerminationMode::kTailBiting),
               std::invalid_argument);

  {
    const std::size_t k = 24;
    auto tb_rsc = ConvolutionalCode::LteTurboConstituent(
        gf2, k, TerminationMode::kTailBiting);
    EXPECT_EQ(tb_rsc.Length(), 2 * k);

    std::mt19937_64 rng(707);
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
    const uint32_t s0 = tb_rsc.ComputeTailBitingInitialState(msg.RawSpan());
    uint32_t state = s0;
    const auto &trellis = tb_rsc.GetTrellis();
    for (std::size_t t = 0; t < k; ++t) {
      state = trellis.NextState(state, msg.RawSpan()[t] & 1U);
    }
    EXPECT_EQ(state, s0);

    auto cw = tb_rsc.Encode(msg);
    EXPECT_TRUE(tb_rsc.Contains(cw));

    BasicAWGNChannel<GaloisFieldBinary> ch =
        BasicAWGNChannel<GaloisFieldBinary>::FromEbN0(gf2, tb_rsc.Length(), 5.0,
                                                      tb_rsc.Rate(), 808);
    auto tx = ch.TransmitDetailed(cw);
    ViterbiDecoder viterbi(tb_rsc);
    auto res_soft = viterbi.DecodeSoft(tx.llr);
    ASSERT_TRUE(res_soft.Ok());
    EXPECT_TRUE(ArraysEqual(res_soft.message, msg));
  }
}

// ============================================================================
// 3. Exact BCJR Verification against Brute-Force Marginalization (K = 8 RSC)
// ============================================================================

TEST(BcjrDecoderTest, ExactLogMapMatchesBruteForceMarginalization) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 8;
  auto rsc = ConvolutionalCode::LteTurboConstituent(
      gf2, k, TerminationMode::kZeroTailed);
  const std::size_t n = rsc.Length();
  ASSERT_EQ(n, (k + 3) * 2);

  BcjrDecoder bcjr(rsc);

  // Generate deterministic noisy channel LLRs and non-zero a priori LLRs.
  std::mt19937_64 rng(909);
  std::normal_distribution<double> llr_dist(0.0, 1.8);
  std::vector<double> channel_llrs(n);
  for (double &v : channel_llrs) {
    v = llr_dist(rng);
  }
  std::vector<double> a_priori_llrs(k);
  for (double &v : a_priori_llrs) {
    v = 0.6 * llr_dist(rng);
  }

  // Test both with zero a priori and with non-trivial a priori LLRs.
  for (int trial = 0; trial < 2; ++trial) {
    const std::vector<double> prior =
        (trial == 0) ? std::vector<double>{} : a_priori_llrs;
    BcjrOutput out = bcjr.DecodeSoft(channel_llrs, prior, BcjrMode::kLogMap);
    ASSERT_EQ(out.app_llrs.size(), k);
    ASSERT_EQ(out.extrinsic_llrs.size(), k);
    ASSERT_EQ(out.hard_bits.Size(), k);

    // Exhaustively enumerate all 2^8 = 256 messages and compute exact posterior
    // log-weights.
    std::vector<std::vector<double>> log_w0(k);
    std::vector<std::vector<double>> log_w1(k);
    for (std::size_t t = 0; t < k; ++t) {
      log_w0[t].reserve(128);
      log_w1[t].reserve(128);
    }

    std::vector<uint8_t> msg_bits(k, 0U);
    for (uint32_t mask = 0; mask < (uint32_t{1} << k); ++mask) {
      double log_prob = 0.0;
      for (std::size_t t = 0; t < k; ++t) {
        const uint8_t bit = static_cast<uint8_t>((mask >> t) & 1U);
        msg_bits[t] = bit;
        if (!prior.empty()) {
          log_prob += (bit == 0U ? 0.5 : -0.5) * prior[t];
        }
      }
      const std::vector<uint8_t> cw = rsc.EncodeBits(msg_bits);
      for (std::size_t j = 0; j < n; ++j) {
        log_prob += (cw[j] == 0U ? 0.5 : -0.5) * channel_llrs[j];
      }
      for (std::size_t t = 0; t < k; ++t) {
        if (msg_bits[t] == 0U) {
          log_w0[t].push_back(log_prob);
        } else {
          log_w1[t].push_back(log_prob);
        }
      }
    }

    for (std::size_t t = 0; t < k; ++t) {
      const double exact_app = LogSumExp(log_w0[t]) - LogSumExp(log_w1[t]);
      EXPECT_NEAR(out.app_llrs[t], exact_app, 1e-11)
          << "Mismatch at trial=" << trial << ", bit=" << t;

      const double la = prior.empty() ? 0.0 : prior[t];
      const double sys_llr = channel_llrs[t * 2];
      EXPECT_NEAR(out.extrinsic_llrs[t], exact_app - la - sys_llr, 1e-11);

      const uint8_t expected_hard =
          (exact_app >= 0.0) ? uint8_t{0} : uint8_t{1};
      EXPECT_EQ(out.hard_bits.RawSpan()[t] & 1U, expected_hard);
    }
  }

  // Max-Log-MAP smoke check on the same LLRs (finite, reasonable agreement).
  BcjrOutput max_out =
      bcjr.DecodeSoft(channel_llrs, a_priori_llrs, BcjrMode::kMaxLogMap);
  ASSERT_EQ(max_out.app_llrs.size(), k);
  for (std::size_t t = 0; t < k; ++t) {
    EXPECT_TRUE(std::isfinite(max_out.app_llrs[t]));
    EXPECT_TRUE(std::isfinite(max_out.extrinsic_llrs[t]));
  }
}

// ============================================================================
// 4. Interleaver Tests (QPP, S-Random, Matrix, Explicit)
// ============================================================================

TEST(InterleaverTest, QppLtePresetsBijectionAndRoundTrip) {
  const std::vector<std::size_t> lte_sizes = {40, 48, 64, 128, 256, 512, 1024};
  for (std::size_t k : lte_sizes) {
    Interleaver il = Interleaver::CreateLteQpp(k);
    EXPECT_EQ(il.Size(), k);

    std::vector<int> data(k);
    std::iota(data.begin(), data.end(), 100);
    auto interleaved = il.Interleave(data);
    auto restored = il.Deinterleave(interleaved);
    EXPECT_EQ(restored, data);

    auto deint_first = il.Deinterleave(data);
    auto restored2 = il.Interleave(deint_first);
    EXPECT_EQ(restored2, data);
  }

  // Non-bijective QPP parameters (gcd(f1, K) != 1) must throw.
  EXPECT_THROW(Interleaver::QppInterleaver(40, 2, 10), std::invalid_argument);
}

TEST(InterleaverTest, SRandomSpreadAndRoundTrip) {
  const std::size_t k = 128;
  const std::size_t s = 6;
  Interleaver il = Interleaver::SRandomInterleaver(k, s, 42);
  EXPECT_EQ(il.Size(), k);

  // Verify S-spread property: |pi(j) - pi(i)| >= s for 0 < j - i < s.
  const auto &perm = il.Permutation();
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = i + 1; j < std::min(k, i + s); ++j) {
      const std::size_t diff =
          (perm[j] > perm[i]) ? (perm[j] - perm[i]) : (perm[i] - perm[j]);
      EXPECT_GE(diff, s) << "Spread violation at i=" << i << ", j=" << j;
    }
  }

  std::vector<double> vals(k);
  for (std::size_t i = 0; i < k; ++i) {
    vals[i] = static_cast<double>(i) * 0.25 - 10.0;
  }
  EXPECT_EQ(il.Deinterleave(il.Interleave(vals)), vals);
}

TEST(InterleaverTest, MatrixAndExplicitInterleaver) {
  Interleaver mat = Interleaver::MatrixInterleaver(4, 5);
  EXPECT_EQ(mat.Size(), 20U);

  std::vector<int> seq(20);
  std::iota(seq.begin(), seq.end(), 0);
  auto out = mat.Interleave(seq);
  // Column 0 of 4x5 row-major matrix is {0, 5, 10, 15}.
  EXPECT_EQ(out[0], 0);
  EXPECT_EQ(out[1], 5);
  EXPECT_EQ(out[2], 10);
  EXPECT_EQ(out[3], 15);
  EXPECT_EQ(out[4], 1);
  EXPECT_EQ(mat.Deinterleave(out), seq);

  Interleaver exp_il = Interleaver::ExplicitInterleaver({2, 0, 3, 1});
  std::vector<int> v4 = {10, 20, 30, 40};
  auto v4_int = exp_il.Interleave(v4);
  EXPECT_EQ(v4_int, (std::vector<int>{30, 10, 40, 20}));
  EXPECT_EQ(exp_il.Deinterleave(v4_int), v4);

  EXPECT_THROW(Interleaver::ExplicitInterleaver({0, 1, 1, 2}),
               std::invalid_argument);
}

// ============================================================================
// 5. 3GPP LTE Turbo Code (K = 40 and K = 128, Rate 1/3 and Punctured Rate 1/2)
// ============================================================================

TEST(TurboCodeTest, LteNoiselessRoundTripAndLinearStructure) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  for (std::size_t k : {std::size_t{40}, std::size_t{128}}) {
    for (bool rate_half : {false, true}) {
      auto turbo = TurboCode::CreateLte(gf2, k, rate_half);
      const std::size_t expected_len = rate_half ? (2 * k + 12) : (3 * k + 12);
      EXPECT_EQ(turbo.Dimension(), k);
      EXPECT_EQ(turbo.Length(), expected_len);

      std::mt19937_64 rng(1000 + k + (rate_half ? 1 : 0));
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
      auto cw = turbo.Encode(msg);
      EXPECT_EQ(cw.Size(), expected_len);

      auto unenc = turbo.Unencode(cw);
      EXPECT_TRUE(ArraysEqual(unenc, msg));

      // Strong noiseless LLRs should converge immediately in 1 iteration.
      std::vector<double> clean_llrs(expected_len);
      for (std::size_t i = 0; i < expected_len; ++i) {
        clean_llrs[i] = ((cw.RawSpan()[i] & 1U) == 0U) ? 8.0 : -8.0;
      }
      auto dec = turbo.DecodeSoft(clean_llrs, 6, BcjrMode::kLogMap, 1.0);
      EXPECT_TRUE(dec.success);
      EXPECT_TRUE(dec.converged);
      EXPECT_EQ(dec.iterations_used, 1U);
      EXPECT_TRUE(ArraysEqual(dec.message, msg));
      EXPECT_TRUE(ArraysEqual(dec.codeword, cw));

      // Verify exact linear code properties on K=40.
      if (k == 40) {
        const auto &g = turbo.GeneratorMatrix();
        const auto &h = turbo.ParityCheckMatrix();
        EXPECT_EQ(g.Shape(0), k);
        EXPECT_EQ(g.Shape(1), expected_len);
        EXPECT_EQ(h.Shape(0), expected_len - k);
        EXPECT_EQ(h.Shape(1), expected_len);
        EXPECT_TRUE(ArraysEqual(linalg::Dot(msg, g), cw));
        EXPECT_TRUE(turbo.Contains(cw));
      }
    }
  }
}

TEST(TurboCodeTest, LteRateOneThirdAwgnIterativeImprovement) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 128;
  auto turbo = TurboCode::CreateLte(gf2, k, /*rate_half_punctured=*/false);

  std::mt19937_64 msg_rng(2025);
  BasicAWGNChannel<GaloisFieldBinary> ch =
      BasicAWGNChannel<GaloisFieldBinary>::FromEbN0(gf2, turbo.Length(), 1.6,
                                                    turbo.Rate(), 2026);

  std::size_t errors_iter1 = 0;
  std::size_t errors_iter6 = 0;
  std::size_t converged_frames = 0;
  const std::size_t num_frames = 20;

  for (std::size_t f = 0; f < num_frames; ++f) {
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, msg_rng);
    auto cw = turbo.Encode(msg);
    auto tx = ch.TransmitDetailed(cw);

    auto res1 = turbo.DecodeSoft(tx.llr, 1, BcjrMode::kLogMap, 1.0);
    auto res6 = turbo.DecodeSoft(tx.llr, 6, BcjrMode::kLogMap, 1.0);

    errors_iter1 += CountBitDiffs(res1.message, msg);
    errors_iter6 += CountBitDiffs(res6.message, msg);
    if (res6.converged) {
      ++converged_frames;
    }
  }

  // Iterative turbo decoding must substantially reduce errors from iter 1 to 6
  // and converge on most frames at Eb/N0 = 1.2 dB.
  EXPECT_LT(errors_iter6, errors_iter1);
  EXPECT_LE(errors_iter6, 15U);
  EXPECT_GT(converged_frames, num_frames / 2);
}

TEST(TurboCodeTest, LteRateHalfPuncturedAwgnDecoding) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  const std::size_t k = 128;
  auto turbo = TurboCode::CreateLte(gf2, k, /*rate_half_punctured=*/true);
  EXPECT_EQ(turbo.Length(), 2 * k + 12);

  std::mt19937_64 msg_rng(3030);
  BasicAWGNChannel<GaloisFieldBinary> ch =
      BasicAWGNChannel<GaloisFieldBinary>::FromEbN0(gf2, turbo.Length(), 2.2,
                                                    turbo.Rate(), 3031);

  std::size_t errors_iter1 = 0;
  std::size_t errors_iter6 = 0;
  const std::size_t num_frames = 20;

  for (std::size_t f = 0; f < num_frames; ++f) {
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, msg_rng);
    auto cw = turbo.Encode(msg);
    auto tx = ch.TransmitDetailed(cw);

    auto res1 = turbo.DecodeSoft(tx.llr, 1, BcjrMode::kLogMap, 1.0);
    auto res6 = turbo.DecodeSoft(tx.llr, 6, BcjrMode::kLogMap, 1.0);

    errors_iter1 += CountBitDiffs(res1.message, msg);
    errors_iter6 += CountBitDiffs(res6.message, msg);
  }

  EXPECT_LT(errors_iter6, errors_iter1);
  EXPECT_LE(errors_iter6, 15U);

  // Also verify Max-Log-MAP with extrinsic scaling 0.75 on K=40 punctured code.
  auto turbo40 = TurboCode::CreateLte(gf2, 40, /*rate_half_punctured=*/true);
  BasicAWGNChannel<GaloisFieldBinary> ch40 =
      BasicAWGNChannel<GaloisFieldBinary>::FromEbN0(gf2, turbo40.Length(), 3.5,
                                                    turbo40.Rate(), 3032);
  auto msg40 = GFArray<GaloisFieldBinary>::Random(gf2, {40}, msg_rng);
  auto cw40 = turbo40.Encode(msg40);
  auto tx40 = ch40.TransmitDetailed(cw40);
  auto res40 = turbo40.DecodeSoft(tx40.llr, 6, BcjrMode::kMaxLogMap, 0.75);
  EXPECT_TRUE(res40.success);
  EXPECT_TRUE(ArraysEqual(res40.message, msg40));
}

}  // namespace
