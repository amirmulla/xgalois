#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <span>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// Co-include channel.hpp alongside all v0.4 channel headers to verify zero
// symbol or macro collisions across the entire xgalois/channel/ module.
#include "xgalois/channel/asymmetric_sync_channel.hpp"
#include "xgalois/channel/burst_channel.hpp"
#include "xgalois/channel/channel.hpp"
#include "xgalois/channel/fading_channel.hpp"
#include "xgalois/channel/modulation.hpp"
#include "xgalois/channel/rng.hpp"
#include "xgalois/channel/simulation.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"

namespace {

using xg::GaloisFieldBinary;
using xg::GaloisFieldExtension;
using xg::GaloisFieldPrime;
using xg::GFArray;
using xg::GFMatrix;
using xg::channel::BinaryAsymmetricChannel;
using xg::channel::BurstChannelMode;
using xg::channel::ComplexSymbol;
using xg::channel::CsiMode;
using xg::channel::DemapMode;
using xg::channel::Demodulator;
using xg::channel::DeriveStreamSeed;
using xg::channel::EditOp;
using xg::channel::ExportSweepCsv;
using xg::channel::ExportSweepJson;
using xg::channel::FrameTrialOutcome;
using xg::channel::GilbertElliottChannel;
using xg::channel::InsertionDeletionSubstitutionChannel;
using xg::channel::InverseNormalCdf;
using xg::channel::LevenshteinDistance;
using xg::channel::ModulationScheme;
using xg::channel::Modulator;
using xg::channel::ParallelSimulationConfig;
using xg::channel::Pcg64;
using xg::channel::RayleighFadingChannel;
using xg::channel::RicianFadingChannel;
using xg::channel::RunParameterSweep;
using xg::channel::SimulateCodePerformanceParallel;
using xg::channel::SplitMix64;
using xg::channel::Uint128;
using xg::channel::WilsonScoreInterval;

// Brute-force 2-D Bayes LLR reference for validating separable and 2-D
// demappers:
//   L(b_k) = ln( sum_{x: b_k=0} exp(-|y - h x|^2 / N_0)
//              / sum_{x: b_k=1} exp(-|y - h x|^2 / N_0) )
std::vector<double> BruteForceExactLogMap(
    const std::vector<ComplexSymbol>& constellation, std::size_t m,
    ComplexSymbol y, ComplexSymbol h, double n0) {
  const std::size_t num_pts = constellation.size();
  std::vector<double> metrics(num_pts, 0.0);
  for (std::size_t idx = 0; idx < num_pts; ++idx) {
    const ComplexSymbol diff = y - h * constellation[idx];
    metrics[idx] = -diff.NormSq() / n0;
  }

  auto log_sum_exp = [](const std::vector<double>& vals) {
    const double mx = *std::max_element(vals.begin(), vals.end());
    double s = 0.0;
    for (double v : vals) {
      s += std::exp(v - mx);
    }
    return mx + std::log(s);
  };

  std::vector<double> llrs(m, 0.0);
  for (std::size_t b = 0; b < m; ++b) {
    const std::size_t mask = std::size_t{1} << (m - 1 - b);
    std::vector<double> v0;
    std::vector<double> v1;
    for (std::size_t idx = 0; idx < num_pts; ++idx) {
      if ((idx & mask) == 0) {
        v0.push_back(metrics[idx]);
      } else {
        v1.push_back(metrics[idx]);
      }
    }
    llrs[b] = log_sum_exp(v0) - log_sum_exp(v1);
  }
  return llrs;
}

std::vector<double> BruteForceMaxLogMap(
    const std::vector<ComplexSymbol>& constellation, std::size_t m,
    ComplexSymbol y, ComplexSymbol h, double n0) {
  const std::size_t num_pts = constellation.size();
  std::vector<double> metrics(num_pts, 0.0);
  for (std::size_t idx = 0; idx < num_pts; ++idx) {
    const ComplexSymbol diff = y - h * constellation[idx];
    metrics[idx] = -diff.NormSq() / n0;
  }

  std::vector<double> llrs(m, 0.0);
  for (std::size_t b = 0; b < m; ++b) {
    const std::size_t mask = std::size_t{1} << (m - 1 - b);
    double max0 = -std::numeric_limits<double>::infinity();
    double max1 = -std::numeric_limits<double>::infinity();
    for (std::size_t idx = 0; idx < num_pts; ++idx) {
      if ((idx & mask) == 0) {
        max0 = std::max(max0, metrics[idx]);
      } else {
        max1 = std::max(max1, metrics[idx]);
      }
    }
    llrs[b] = max0 - max1;
  }
  return llrs;
}

// ============================================================================
// 1. PRNG Tests (SplitMix64, Pcg64, DeriveStreamSeed)
// ============================================================================

TEST(RngTest, SplitMix64DiscardAndUrbgCompliance) {
  static_assert(std::uniform_random_bit_generator<SplitMix64>);
  SplitMix64 a(123456789ULL);
  SplitMix64 b(123456789ULL);

  for (uint64_t i = 0; i < 1000; ++i) {
    (void)a();
  }
  b.Discard(1000);
  EXPECT_EQ(a, b);
  for (int i = 0; i < 32; ++i) {
    EXPECT_EQ(a(), b());
  }

  std::uniform_int_distribution<int> dist(0, 99);
  for (int i = 0; i < 100; ++i) {
    const int v = dist(a);
    EXPECT_GE(v, 0);
    EXPECT_LE(v, 99);
  }
}

TEST(RngTest, Pcg64DiscardExactBitIdentical) {
  static_assert(std::uniform_random_bit_generator<Pcg64>);

  const std::array<uint64_t, 7> deltas = {0, 1, 2, 3, 17, 1000, 65537};
  for (uint64_t delta : deltas) {
    Pcg64 seq(0xdeadbeefcafebabeULL, 0x123456789abcdef0ULL);
    Pcg64 jumped(0xdeadbeefcafebabeULL, 0x123456789abcdef0ULL);

    for (uint64_t i = 0; i < delta; ++i) {
      (void)seq();
    }
    jumped.Discard(delta);

    EXPECT_EQ(seq.State(), jumped.State()) << "Mismatch at delta=" << delta;
    for (int k = 0; k < 20; ++k) {
      EXPECT_EQ(seq(), jumped()) << "Output mismatch at delta=" << delta;
    }
  }
}

TEST(RngTest, Pcg64IndependentStreamsAnd128BitSeed) {
  Pcg64 s1(42ULL, 1ULL);
  Pcg64 s2(42ULL, 2ULL);
  EXPECT_NE(s1(), s2());

  Pcg64 hi1(Uint128{1ULL, 2ULL}, Uint128{0ULL, 10ULL});
  Pcg64 hi2(Uint128{1ULL, 2ULL}, Uint128{1ULL, 10ULL});
  EXPECT_NE(hi1(), hi2());
}

TEST(RngTest, DeriveStreamSeedUniquenessAndDeterminism) {
  std::unordered_set<uint64_t> seen;
  for (uint64_t sweep = 0; sweep < 10; ++sweep) {
    for (uint64_t chunk = 0; chunk < 200; ++chunk) {
      const uint64_t s = DeriveStreamSeed(42ULL, sweep, chunk);
      EXPECT_EQ(s, DeriveStreamSeed(42ULL, sweep, chunk));
      EXPECT_TRUE(seen.insert(s).second);
    }
  }
}

// ============================================================================
// 2. Gilbert-Elliott Burst Channel Tests
// ============================================================================

TEST(BurstChannelTest, AnalyticalAndEmpiricalStationaryStatistics) {
  using GF256 = GaloisFieldExtension<uint16_t>;
  auto gf256 = std::make_shared<GF256>(std::pair{2ULL, 8ULL});

  const double avg_burst = 10.0;
  const double avg_good = 90.0;
  const double err_good = 0.005;
  const double err_bad = 0.80;

  auto ch = GilbertElliottChannel<GF256>::FromBurstStats(
      gf256, avg_burst, avg_good, err_good, err_bad,
      BurstChannelMode::kSymbolError, 202601ULL);

  EXPECT_NEAR(ch.SteadyStateBadProb(), 0.1, 1e-12);
  EXPECT_NEAR(ch.SteadyStateGoodProb(), 0.9, 1e-12);
  EXPECT_NEAR(ch.AverageBurstLength(), 10.0, 1e-12);
  EXPECT_NEAR(ch.AverageGoodLength(), 90.0, 1e-12);
  EXPECT_NEAR(ch.SteadyStateErrorRate(), 0.9 * 0.005 + 0.1 * 0.80, 1e-12);

  constexpr std::size_t kLen = 200000;
  auto cw = GFArray<GF256>::Ones(gf256, {kLen});
  auto res = ch.Transmit(cw);

  ASSERT_EQ(res.state_trace.size(), kLen);
  EXPECT_TRUE(res.erasure_positions.empty());

  std::size_t bad_count = 0;
  std::size_t bad_runs = 0;
  for (std::size_t i = 0; i < kLen; ++i) {
    if (res.state_trace[i] == 1) {
      ++bad_count;
      if (i == 0 || res.state_trace[i - 1] == 0) {
        ++bad_runs;
      }
    }
  }

  const double empirical_bad_prob =
      static_cast<double>(bad_count) / static_cast<double>(kLen);
  const double empirical_avg_burst =
      static_cast<double>(bad_count) / static_cast<double>(bad_runs);
  const double empirical_err_rate =
      static_cast<double>(res.error_positions.size()) /
      static_cast<double>(kLen);

  EXPECT_NEAR(empirical_bad_prob, ch.SteadyStateBadProb(), 0.015);
  EXPECT_NEAR(empirical_avg_burst, ch.AverageBurstLength(), 1.2);
  EXPECT_NEAR(empirical_err_rate, ch.SteadyStateErrorRate(), 0.015);
  EXPECT_GT(res.num_burst_transitions, 100U);

  // Verify exact error_positions correspondence.
  std::vector<bool> is_err(kLen, false);
  for (std::size_t pos : res.error_positions) {
    is_err[pos] = true;
  }
  const auto tx_span = cw.RawSpan();
  const auto rx_span = res.received_word.RawSpan();
  for (std::size_t i = 0; i < kLen; ++i) {
    if (is_err[i]) {
      EXPECT_NE(rx_span[i], tx_span[i]);
    } else {
      EXPECT_EQ(rx_span[i], tx_span[i]);
    }
  }
}

TEST(BurstChannelTest, MatrixSynchronizedBurstsAndErasureMode) {
  using GF17 = GaloisFieldPrime<uint8_t>;
  auto gf17 = std::make_shared<GF17>(17);
  constexpr std::size_t kRows = 64;
  constexpr std::size_t kCols = 5;

  // Build non-zero n x ell matrix (shape {64, 5}).
  auto mat = GFMatrix<GF17>::Ones(gf17, {kRows, kCols});
  for (std::size_t r = 0; r < kRows; ++r) {
    for (std::size_t c = 0; c < kCols; ++c) {
      mat.SetExternal(r, c, static_cast<uint8_t>((r + c) % 16 + 1));
    }
  }

  // Error mode on 2-D IRS matrix.
  auto err_ch = GilbertElliottChannel<GF17>::FromBurstStats(
      gf17, 4.0, 12.0, 0.0, 1.0, BurstChannelMode::kSymbolError, 777ULL);
  auto err_res = err_ch.TransmitMatrix(mat);
  EXPECT_EQ(err_res.state_trace.size(), kRows);
  EXPECT_TRUE(err_res.erasure_positions.empty());
  EXPECT_FALSE(err_res.error_positions.empty());

  std::vector<bool> hit_row(kRows, false);
  for (std::size_t r : err_res.error_positions) {
    hit_row[r] = true;
  }
  for (std::size_t r = 0; r < kRows; ++r) {
    for (std::size_t c = 0; c < kCols; ++c) {
      if (hit_row[r]) {
        EXPECT_NE(err_res.received_word.Raw(r, c), mat.Raw(r, c));
      } else {
        EXPECT_EQ(err_res.received_word.Raw(r, c), mat.Raw(r, c));
      }
    }
  }

  // Erasure mode on 2-D matrix (err_good=0, err_bad=1 => erased iff Bad).
  auto era_ch = GilbertElliottChannel<GF17>::FromBurstStats(
      gf17, 4.0, 12.0, 0.0, 1.0, BurstChannelMode::kSymbolErasure, 777ULL);
  auto era_res = era_ch.TransmitMatrix(mat);
  EXPECT_EQ(era_res.state_trace.size(), kRows);
  EXPECT_TRUE(era_res.error_positions.empty());
  EXPECT_FALSE(era_res.erasure_positions.empty());

  std::vector<bool> era_hit_row(kRows, false);
  for (std::size_t r : era_res.erasure_positions) {
    era_hit_row[r] = true;
  }
  const auto zero = gf17->AdditiveIdentity();
  for (std::size_t r = 0; r < kRows; ++r) {
    EXPECT_EQ(era_hit_row[r], era_res.state_trace[r] == 1);
    for (std::size_t c = 0; c < kCols; ++c) {
      if (era_hit_row[r]) {
        EXPECT_EQ(era_res.received_word.Raw(r, c), zero);
      } else {
        EXPECT_EQ(era_res.received_word.Raw(r, c), mat.Raw(r, c));
      }
    }
  }

  // CloneWithSeed determinism.
  auto clone_a = era_ch.CloneWithSeed(999ULL);
  auto clone_b = era_ch.CloneWithSeed(999ULL);
  auto ra = clone_a.TransmitMatrix(mat);
  auto rb = clone_b.TransmitMatrix(mat);
  EXPECT_EQ(ra.received_word, rb.received_word);
  EXPECT_EQ(ra.erasure_positions, rb.erasure_positions);
  EXPECT_EQ(ra.state_trace, rb.state_trace);
}

// ============================================================================
// 3. Modulation & Demodulation Tests
// ============================================================================

TEST(ModulationTest, UnitAverageEnergyAndGrayCodingAllSchemes) {
  const std::array<ModulationScheme, 5> schemes = {
      ModulationScheme::kBpsk, ModulationScheme::kQpsk, ModulationScheme::k8Psk,
      ModulationScheme::k16Qam, ModulationScheme::k64Qam};
  const std::array<std::size_t, 5> expected_m = {1, 2, 3, 4, 6};

  for (std::size_t idx = 0; idx < schemes.size(); ++idx) {
    Modulator mod(schemes[idx]);
    EXPECT_EQ(mod.BitsPerSymbol(), expected_m[idx]);
    const auto& pts = mod.Constellation();
    ASSERT_EQ(pts.size(), std::size_t{1} << expected_m[idx]);

    // Verify exact unit average energy E_s = 1.0.
    double energy_sum = 0.0;
    for (const auto& s : pts) {
      energy_sum += s.NormSq();
    }
    const double avg_energy = energy_sum / static_cast<double>(pts.size());
    EXPECT_NEAR(avg_energy, 1.0, 1e-13) << "Scheme idx=" << idx;

    // Verify Gray coding: all minimum-distance neighbor pairs differ by 1 bit.
    double min_dist_sq = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < pts.size(); ++i) {
      for (std::size_t j = i + 1; j < pts.size(); ++j) {
        min_dist_sq = std::min(min_dist_sq, (pts[i] - pts[j]).NormSq());
      }
    }
    for (std::size_t i = 0; i < pts.size(); ++i) {
      for (std::size_t j = i + 1; j < pts.size(); ++j) {
        const double d2 = (pts[i] - pts[j]).NormSq();
        if (d2 <= min_dist_sq * 1.0001) {
          const unsigned int xor_bits = static_cast<unsigned int>(i ^ j);
          EXPECT_EQ(__builtin_popcount(xor_bits), 1)
              << "Non-Gray nearest neighbors " << i << " and " << j
              << " in scheme idx=" << idx;
        }
      }
    }
  }
}

TEST(ModulationTest, ZeroNoiseRoundTripAndExactBayesLlrVerification) {
  const std::array<ModulationScheme, 5> schemes = {
      ModulationScheme::kBpsk, ModulationScheme::kQpsk, ModulationScheme::k8Psk,
      ModulationScheme::k16Qam, ModulationScheme::k64Qam};

  Pcg64 rng(314159ULL);

  for (ModulationScheme scheme : schemes) {
    Modulator mod(scheme);
    Demodulator demod(scheme);
    const std::size_t m = mod.BitsPerSymbol();
    const std::size_t num_syms = std::size_t{1} << m;

    // Enumerate every constellation label once.
    std::vector<uint8_t> bits;
    bits.reserve(num_syms * m);
    for (std::size_t lbl = 0; lbl < num_syms; ++lbl) {
      for (std::size_t b = 0; b < m; ++b) {
        bits.push_back(static_cast<uint8_t>((lbl >> (m - 1 - b)) & 1U));
      }
    }

    const auto tx = mod.Modulate(bits);
    ASSERT_EQ(tx.size(), num_syms);

    // Zero-noise hard & soft round-trip with unit gains.
    const auto hard_unit = demod.DemodulateHard(tx);
    EXPECT_EQ(hard_unit, bits);

    const auto soft_zero =
        demod.DemodulateSoft(tx, 0.0, DemapMode::kExactLogMap);
    ASSERT_EQ(soft_zero.size(), bits.size());
    for (std::size_t i = 0; i < bits.size(); ++i) {
      if (bits[i] == 0) {
        EXPECT_GT(soft_zero[i], 0.0);
      } else {
        EXPECT_LT(soft_zero[i], 0.0);
      }
    }

    // Noisy reception with arbitrary complex fading gains h_s.
    const double n0 = 0.25;
    std::vector<ComplexSymbol> gains(num_syms);
    std::vector<ComplexSymbol> rx(num_syms);
    for (std::size_t s = 0; s < num_syms; ++s) {
      const auto [hi, hq] = rng.NextGaussianPair();
      const auto [ni, nq] = rng.NextGaussianPair();
      // Keep gain magnitude away from 0.
      gains[s] = ComplexSymbol{0.7 + 0.5 * hi, 0.3 + 0.5 * hq};
      const double sigma = std::sqrt(0.5 * n0);
      rx[s] = gains[s] * tx[s] + ComplexSymbol{sigma * ni, sigma * nq};
    }

    // Verify Exact Log-MAP matches brute-force 2-D Bayes formula to < 1e-11!
    const auto exact_llrs =
        demod.DemodulateSoft(rx, n0, DemapMode::kExactLogMap, gains);
    const auto maxlog_llrs =
        demod.DemodulateSoft(rx, n0, DemapMode::kMaxLogMap, gains);

    for (std::size_t s = 0; s < num_syms; ++s) {
      const auto ref_exact =
          BruteForceExactLogMap(mod.Constellation(), m, rx[s], gains[s], n0);
      const auto ref_max =
          BruteForceMaxLogMap(mod.Constellation(), m, rx[s], gains[s], n0);
      for (std::size_t b = 0; b < m; ++b) {
        EXPECT_NEAR(exact_llrs[s * m + b], ref_exact[b], 1e-11);
        EXPECT_NEAR(maxlog_llrs[s * m + b], ref_max[b], 1e-11);
      }
    }
  }
}

// ============================================================================
// 4. Fading Channels (Rayleigh & Rician) Tests
// ============================================================================

TEST(FadingChannelTest, RayleighEmpiricalBpskBerMatchesAnalyticalFormula) {
  constexpr std::size_t kNumBits = 120000;
  std::vector<uint8_t> tx_bits(kNumBits);
  Pcg64 bit_rng(1001ULL);
  for (std::size_t i = 0; i < kNumBits; ++i) {
    tx_bits[i] = static_cast<uint8_t>(bit_rng() & 1ULL);
  }

  for (double snr_db : {0.0, 4.0, 8.0}) {
    RayleighFadingChannel ch(snr_db, ModulationScheme::kBpsk,
                             /*coherence_symbols=*/1, CsiMode::kPerfect, 0.0,
                             2026ULL + static_cast<uint64_t>(snr_db * 10.0));
    const auto res = ch.Transmit(tx_bits);
    ASSERT_EQ(res.llrs.size(), kNumBits);

    std::size_t bit_errors = 0;
    for (std::size_t i = 0; i < kNumBits; ++i) {
      const uint8_t dec = (res.llrs[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
      if (dec != tx_bits[i]) {
        ++bit_errors;
      }
    }
    const double empirical_ber =
        static_cast<double>(bit_errors) / static_cast<double>(kNumBits);
    const double analytical_ber =
        RayleighFadingChannel::AnalyticalRayleighBpskBer(snr_db);

    // 5-sigma binomial tolerance.
    const double std_err =
        std::sqrt(analytical_ber * (1.0 - analytical_ber) / kNumBits);
    EXPECT_NEAR(empirical_ber, analytical_ber, 5.0 * std_err + 1e-4)
        << "SNR=" << snr_db << " dB";
  }
}

TEST(FadingChannelTest, BlockCoherenceCsiDegradationAndRicianHighK) {
  constexpr std::size_t kNumBits = 40000;
  std::vector<uint8_t> tx_bits(kNumBits);
  Pcg64 bit_rng(555ULL);
  for (std::size_t i = 0; i < kNumBits; ++i) {
    tx_bits[i] = static_cast<uint8_t>(bit_rng() & 1ULL);
  }

  // Block coherence check (coherence_symbols = 8).
  RayleighFadingChannel block_ch(6.0, ModulationScheme::kQpsk,
                                 /*coherence_symbols=*/8, CsiMode::kPerfect,
                                 0.0, 111ULL);
  const auto block_res =
      block_ch.Transmit(std::span<const uint8_t>(tx_bits).subspan(0, 64));
  ASSERT_EQ(block_res.channel_gains.size(), 32U);
  for (std::size_t blk = 0; blk < 4; ++blk) {
    const ComplexSymbol g0 = block_res.channel_gains[blk * 8];
    for (std::size_t k = 1; k < 8; ++k) {
      EXPECT_EQ(block_res.channel_gains[blk * 8 + k], g0);
    }
  }
  EXPECT_NE(block_res.channel_gains[0], block_res.channel_gains[8]);

  // Perfect CSI strictly beats Noisy CSI on Rayleigh BPSK.
  RayleighFadingChannel perfect_ch(8.0, ModulationScheme::kBpsk, 1,
                                   CsiMode::kPerfect, 0.0, 909ULL);
  RayleighFadingChannel noisy_ch(8.0, ModulationScheme::kBpsk, 1,
                                 CsiMode::kNoisyEstimate, 0.25, 909ULL);
  const auto res_perf = perfect_ch.Transmit(tx_bits);
  const auto res_noisy = noisy_ch.Transmit(tx_bits);

  std::size_t err_perf = 0;
  std::size_t err_noisy = 0;
  for (std::size_t i = 0; i < kNumBits; ++i) {
    if ((res_perf.llrs[i] >= 0.0 ? 0 : 1) != tx_bits[i]) {
      ++err_perf;
    }
    if ((res_noisy.llrs[i] >= 0.0 ? 0 : 1) != tx_bits[i]) {
      ++err_noisy;
    }
  }
  EXPECT_LT(err_perf, err_noisy);

  // Rician high-K (K = 30 dB = 1000) converges toward AWGN BPSK BER.
  const double snr_db = 4.0;
  auto rician_hi =
      RicianFadingChannel::FromKFactorDb(snr_db, 30.0, ModulationScheme::kBpsk,
                                         1, CsiMode::kPerfect, 0.0, 4242ULL);
  const auto rician_res = rician_hi.Transmit(tx_bits);
  std::size_t rician_errs = 0;
  for (std::size_t i = 0; i < kNumBits; ++i) {
    if ((rician_res.llrs[i] >= 0.0 ? 0 : 1) != tx_bits[i]) {
      ++rician_errs;
    }
  }
  const double rician_ber =
      static_cast<double>(rician_errs) / static_cast<double>(kNumBits);
  const double gamma = std::pow(10.0, snr_db / 10.0);
  const double awgn_ber = 0.5 * std::erfc(std::sqrt(gamma));
  const double rayleigh_ber =
      RayleighFadingChannel::AnalyticalRayleighBpskBer(snr_db);

  EXPECT_NEAR(rician_ber, awgn_ber, 0.004);
  EXPECT_LT(rician_ber, rayleigh_ber * 0.5);
}

// ============================================================================
// 5. Binary Asymmetric Channel (BAC / Z-Channel) & IDS Channel Tests
// ============================================================================

TEST(AsymmetricSyncChannelTest, BacSymmetryZChannelCapacityAndLlrs) {
  // Symmetric BAC reduces to BSC.
  const double p = 0.1;
  BinaryAsymmetricChannel sym(p, p, 42ULL);
  const double h2 = -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
  EXPECT_NEAR(sym.Capacity(), 1.0 - h2, 1e-12);

  const std::vector<uint8_t> sample_rx = {0, 1};
  const auto sym_llrs = sym.ComputeLlrs(sample_rx, 50.0);
  EXPECT_NEAR(sym_llrs[0], std::log((1.0 - p) / p), 1e-12);
  EXPECT_NEAR(sym_llrs[1], -sym_llrs[0], 1e-12);
  EXPECT_GT(sym_llrs[0], 0.0);
  EXPECT_LT(sym_llrs[1], 0.0);

  // Z-channel: 0 -> 1 never flips; exact closed-form capacity.
  const double p10 = 0.25;
  auto zch = BinaryAsymmetricChannel::ZChannel(p10, 777ULL);
  const double expected_z_cap =
      std::log2(1.0 + (1.0 - p10) * std::pow(p10, p10 / (1.0 - p10)));
  EXPECT_NEAR(zch.Capacity(), expected_z_cap, 1e-12);

  // Z-channel LLRs: receiving 1 implies X=1 with certainty (-max_llr).
  const auto z_llrs = zch.ComputeLlrs(sample_rx, 40.0);
  EXPECT_NEAR(z_llrs[0], std::log(1.0 / p10), 1e-12);
  EXPECT_DOUBLE_EQ(z_llrs[1], -40.0);

  // Empirical Z-channel transmission check.
  std::vector<uint8_t> zeros(50000, 0);
  auto rx_zeros = zch.Transmit(zeros);
  EXPECT_EQ(rx_zeros, zeros);  // Zero 0 -> 1 flips!

  std::vector<uint8_t> ones(50000, 1);
  auto rx_ones = zch.Transmit(ones);
  std::size_t flipped = 0;
  for (uint8_t b : rx_ones) {
    if (b == 0) {
      ++flipped;
    }
  }
  const double emp_p10 = static_cast<double>(flipped) / 50000.0;
  EXPECT_NEAR(emp_p10, p10, 0.01);

  // General asymmetric BAC capacity bounds and numerical verification.
  BinaryAsymmetricChannel asym(0.05, 0.30, 1ULL);
  const double cap = asym.Capacity();
  EXPECT_GT(cap, 0.0);
  EXPECT_LT(cap, 1.0);
  // Output-flip invariance: C(a, b) == C(1 - a, 1 - b).
  EXPECT_NEAR(cap, BinaryAsymmetricChannel(0.95, 0.70, 1ULL).Capacity(), 1e-12);
  // Degenerate zero-capacity channel: p01 + p10 == 1.
  EXPECT_DOUBLE_EQ(BinaryAsymmetricChannel(0.3, 0.7, 1ULL).Capacity(), 0.0);
}

TEST(AsymmetricSyncChannelTest, LevenshteinAndIdsChannelInvariants) {
  // Known edit distance cases.
  const std::string s1 = "kitten";
  const std::string s2 = "sitting";
  EXPECT_EQ(LevenshteinDistance<char>(s1, s2), 3U);
  EXPECT_EQ(LevenshteinDistance<char>(s1, s1), 0U);
  EXPECT_EQ(LevenshteinDistance<char>(std::span<const char>{}, s2), s2.size());

  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto tx = GFArray<GaloisFieldBinary>::Zeros(gf2, {128});
  Pcg64 rng(98765ULL);
  for (std::size_t i = 0; i < 128; ++i) {
    tx.Raw(i) = static_cast<uint8_t>(rng() & 1ULL);
  }

  // Test both single-insertion and geometric multi-insertion modes.
  for (std::size_t max_ins : {std::size_t{1}, std::size_t{4}}) {
    InsertionDeletionSubstitutionChannel<GaloisFieldBinary> ids(
        gf2, /*p_ins=*/0.08, /*p_del=*/0.07, /*p_sub=*/0.05, max_ins, 321ULL);

    for (int trial = 0; trial < 100; ++trial) {
      const auto res = ids.Transmit(tx);
      const std::size_t expected_len =
          tx.Size() + res.num_insertions - res.num_deletions;
      EXPECT_EQ(res.received_word.Size(), expected_len);

      // Edit trace must exactly account for all ops.
      std::size_t n_ins = 0;
      std::size_t n_del = 0;
      std::size_t n_sub = 0;
      std::size_t n_mat = 0;
      for (EditOp op : res.edit_trace) {
        if (op == EditOp::kInsert) ++n_ins;
        if (op == EditOp::kDelete) ++n_del;
        if (op == EditOp::kSubstitute) ++n_sub;
        if (op == EditOp::kMatch) ++n_mat;
      }
      EXPECT_EQ(n_ins, res.num_insertions);
      EXPECT_EQ(n_del, res.num_deletions);
      EXPECT_EQ(n_sub, res.num_substitutions);
      EXPECT_EQ(n_del + n_sub + n_mat, tx.Size());

      const std::size_t dist = LevenshteinDistance(tx, res.received_word);
      EXPECT_LE(dist,
                res.num_insertions + res.num_deletions + res.num_substitutions);
    }
  }

  // Edge case: all symbols deleted -> valid empty GFArray.
  InsertionDeletionSubstitutionChannel<GaloisFieldBinary> del_all(gf2, 0.0, 1.0,
                                                                  0.0, 1, 1ULL);
  const auto empty_res = del_all.Transmit(tx);
  EXPECT_EQ(empty_res.received_word.Size(), 0U);
  EXPECT_EQ(empty_res.num_deletions, tx.Size());
  EXPECT_EQ(LevenshteinDistance(tx, empty_res.received_word), tx.Size());
}

// ============================================================================
// 6. Parallel Simulation Engine, Statistics, & Exporters Tests
// ============================================================================

TEST(SimulationTest, InverseNormalCdfAndWilsonScoreInterval) {
  EXPECT_NEAR(InverseNormalCdf(0.5), 0.0, 1e-12);
  EXPECT_NEAR(InverseNormalCdf(0.975), 1.959963984540054, 1e-9);
  EXPECT_NEAR(InverseNormalCdf(0.995), 2.5758293035489004, 1e-9);
  EXPECT_NEAR(InverseNormalCdf(0.025), -1.959963984540054, 1e-9);
  EXPECT_NEAR(InverseNormalCdf(0.001), -3.090232306167813, 1e-9);

  // Round-trip Phi(Phi^{-1}(p)) == p across center and tails.
  for (double p : {1e-8, 1e-4, 0.01, 0.02424, 0.02426, 0.25, 0.75, 0.97574,
                   0.97576, 0.9999, 1.0 - 1e-8}) {
    const double z = InverseNormalCdf(p);
    const double cdf = 0.5 * std::erfc(-z / std::numbers::sqrt2);
    EXPECT_NEAR(cdf, p, 1e-9) << "Failed at p=" << p;
  }

  EXPECT_THROW((void)InverseNormalCdf(0.0), std::domain_error);
  EXPECT_THROW((void)InverseNormalCdf(1.0), std::domain_error);

  // Wilson score CI checks.
  const auto [lo0, hi0] = WilsonScoreInterval(0, 100, 0.95);
  EXPECT_DOUBLE_EQ(lo0, 0.0);
  EXPECT_GT(hi0, 0.03);
  EXPECT_LT(hi0, 0.05);

  const auto [lo50, hi50] = WilsonScoreInterval(50, 100, 0.95);
  EXPECT_NEAR(0.5 * (lo50 + hi50), 0.5, 1e-12);
  EXPECT_LT(lo50, 0.5);
  EXPECT_GT(hi50, 0.5);

  const auto [lo_big, hi_big] = WilsonScoreInterval(5000, 10000, 0.95);
  EXPECT_LT(hi_big - lo_big, hi50 - lo50);
}

TEST(SimulationTest, BitIdenticalAcrossThreadCountsAndStoppingModes) {
  auto make_worker = [](double err_prob) {
    return [err_prob](Pcg64& rng) -> FrameTrialOutcome {
      FrameTrialOutcome out;
      out.symbols = 32;
      out.iterations = 1.0 + static_cast<double>(rng.NextBounded(8));
      out.repair_symbols_read = 4.0 + static_cast<double>(rng.NextBounded(3));
      out.used_local_repair = (rng.NextUniform01() < 0.75);
      for (uint64_t i = 0; i < out.symbols; ++i) {
        if (rng.NextUniform01() < err_prob) {
          ++out.symbol_errors;
        }
      }
      out.frame_error = (out.symbol_errors > 0);
      return out;
    };
  };

  // Case A: Fixed max_frames (no early stop).
  {
    ParallelSimulationConfig cfg1;
    cfg1.max_frames = 1000;
    cfg1.target_frame_errors = 0;
    cfg1.chunk_size = 64;
    cfg1.base_seed = 424242ULL;
    cfg1.num_threads = 1;

    ParallelSimulationConfig cfg4 = cfg1;
    cfg4.num_threads = 4;

    const auto r1 = SimulateCodePerformanceParallel(cfg1, make_worker(0.01));
    const auto r4 = SimulateCodePerformanceParallel(cfg4, make_worker(0.01));

    EXPECT_EQ(r1.total_frames, 1000U);
    EXPECT_EQ(r1.total_frames, r4.total_frames);
    EXPECT_EQ(r1.frame_errors, r4.frame_errors);
    EXPECT_EQ(r1.total_symbols, r4.total_symbols);
    EXPECT_EQ(r1.symbol_errors, r4.symbol_errors);
    EXPECT_DOUBLE_EQ(r1.fer, r4.fer);
    EXPECT_DOUBLE_EQ(r1.ser, r4.ser);
    EXPECT_DOUBLE_EQ(r1.fer_ci95.first, r4.fer_ci95.first);
    EXPECT_DOUBLE_EQ(r1.fer_ci95.second, r4.fer_ci95.second);
    EXPECT_DOUBLE_EQ(r1.avg_decoder_iterations, r4.avg_decoder_iterations);
    EXPECT_DOUBLE_EQ(r1.avg_repair_symbols_read, r4.avg_repair_symbols_read);
    EXPECT_DOUBLE_EQ(r1.local_repair_fraction, r4.local_repair_fraction);
  }

  // Case B: Adaptive early stopping on target_frame_errors must ALSO be 100%
  // bit-identical across num_threads = 1, 2, 4, 8!
  {
    ParallelSimulationConfig base_cfg;
    base_cfg.max_frames = 50000;
    base_cfg.target_frame_errors = 75;
    base_cfg.chunk_size = 32;
    base_cfg.base_seed = 99999ULL;

    base_cfg.num_threads = 1;
    const auto ref =
        SimulateCodePerformanceParallel(base_cfg, make_worker(0.008));
    EXPECT_GE(ref.frame_errors, 75U);
    EXPECT_LT(ref.total_frames, 50000U);

    for (std::size_t nt : {std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
      ParallelSimulationConfig cfg = base_cfg;
      cfg.num_threads = nt;
      const auto res = SimulateCodePerformanceParallel(cfg, make_worker(0.008));
      EXPECT_EQ(res.total_frames, ref.total_frames) << "nt=" << nt;
      EXPECT_EQ(res.frame_errors, ref.frame_errors) << "nt=" << nt;
      EXPECT_EQ(res.symbol_errors, ref.symbol_errors) << "nt=" << nt;
      EXPECT_DOUBLE_EQ(res.fer, ref.fer) << "nt=" << nt;
      EXPECT_DOUBLE_EQ(res.ser, ref.ser) << "nt=" << nt;
      EXPECT_DOUBLE_EQ(res.avg_decoder_iterations, ref.avg_decoder_iterations)
          << "nt=" << nt;
    }
  }

  // Case C: Adaptive Wilson CI relative half-width stopping is bit-identical
  // for num_threads = 1 vs 4.
  {
    ParallelSimulationConfig cfg1;
    cfg1.max_frames = 50000;
    cfg1.target_frame_errors = 0;
    cfg1.min_ci_relative_half_width = 0.20;
    cfg1.chunk_size = 64;
    cfg1.base_seed = 7777ULL;
    cfg1.num_threads = 1;

    ParallelSimulationConfig cfg4 = cfg1;
    cfg4.num_threads = 4;

    const auto r1 = SimulateCodePerformanceParallel(cfg1, make_worker(0.02));
    const auto r4 = SimulateCodePerformanceParallel(cfg4, make_worker(0.02));
    EXPECT_EQ(r1.total_frames, r4.total_frames);
    EXPECT_EQ(r1.frame_errors, r4.frame_errors);
    EXPECT_LT(r1.total_frames, 50000U);
    const double rel_hw =
        0.5 * (r1.fer_ci95.second - r1.fer_ci95.first) / r1.fer;
    EXPECT_LE(rel_hw, 0.20);
  }
}

TEST(SimulationTest, ParameterSweepAndCsvJsonExporters) {
  ParallelSimulationConfig cfg;
  cfg.max_frames = 1024;
  cfg.target_frame_errors = 0;
  cfg.chunk_size = 64;
  cfg.base_seed = 12345ULL;
  cfg.num_threads = 4;

  auto sweep_factory = [](double p_err) {
    return [p_err](Pcg64& rng) -> FrameTrialOutcome {
      FrameTrialOutcome out;
      out.symbols = 16;
      out.iterations = 3.0;
      out.repair_symbols_read = 2.0;
      out.used_local_repair = true;
      for (uint64_t i = 0; i < out.symbols; ++i) {
        if (rng.NextUniform01() < p_err) {
          ++out.symbol_errors;
        }
      }
      out.frame_error = (out.symbol_errors > 0);
      return out;
    };
  };

  const auto report = RunParameterSweep(
      "symbol_error_prob", {0.05, 0.01, 0.001}, cfg, sweep_factory);
  ASSERT_EQ(report.schema_version, "1.0");
  ASSERT_EQ(report.param_name, "symbol_error_prob");
  ASSERT_EQ(report.points.size(), 3U);
  EXPECT_GT(report.points[0].fer, report.points[1].fer);
  EXPECT_GT(report.points[1].fer, report.points[2].fer);

  // Verify CSV export format and line count (1 header + 3 data rows).
  const std::string csv = ExportSweepCsv(report);
  std::istringstream csv_stream(csv);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(csv_stream, line)) {
    lines.push_back(line);
  }
  ASSERT_EQ(lines.size(), 4U);
  EXPECT_NE(lines[0].find("symbol_error_prob,total_frames,frame_errors"),
            std::string::npos);
  EXPECT_NE(lines[0].find("fer_ci95_low,fer_ci95_high"), std::string::npos);

  // Verify JSON export contains schema v1.0, param name, and balanced braces.
  const std::string json = ExportSweepJson(report);
  EXPECT_NE(json.find("\"schema_version\": \"1.0\""), std::string::npos);
  EXPECT_NE(json.find("\"param_name\": \"symbol_error_prob\""),
            std::string::npos);
  EXPECT_NE(json.find("\"fer_ci95\": ["), std::string::npos);
  EXPECT_EQ(std::count(json.begin(), json.end(), '{'),
            std::count(json.begin(), json.end(), '}'));
  EXPECT_EQ(std::count(json.begin(), json.end(), '['),
            std::count(json.begin(), json.end(), ']'));
}

}  // namespace
