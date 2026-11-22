#ifndef XGALOIS_CHANNEL_FADING_CHANNEL_HPP_
#define XGALOIS_CHANNEL_FADING_CHANNEL_HPP_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/channel/modulation.hpp"
#include "xgalois/channel/rng.hpp"

namespace xg {
namespace channel {

enum class CsiMode : uint8_t {
  kPerfect,
  kNoisyEstimate,
  kNone,
};

struct FadingTransmissionResult {
  std::vector<ComplexSymbol> received_symbols;
  std::vector<ComplexSymbol> channel_gains;
  std::vector<ComplexSymbol> estimated_gains;
  std::vector<double> llrs;
};

namespace detail {

inline double SnrDbToNoiseVarianceN0(double snr_db, const char* caller) {
  if (std::isnan(snr_db) || std::isinf(snr_db)) {
    throw std::invalid_argument(std::string(caller) +
                                ": snr_db must be finite");
  }
  // Constellation symbols are normalized to unit average energy E_s = 1.0,
  // so E_s / N_0 = 10^(snr_db / 10) implies N_0 = 10^(-snr_db / 10).
  return std::pow(10.0, -snr_db / 10.0);
}

inline ComplexSymbol SampleComplexGaussian(Pcg64& rng,
                                           double total_variance) noexcept {
  if (total_variance <= 0.0) {
    return ComplexSymbol{0.0, 0.0};
  }
  const double sigma_1d = std::sqrt(0.5 * total_variance);
  const auto [z_i, z_q] = rng.NextGaussianPair();
  return ComplexSymbol{sigma_1d * z_i, sigma_1d * z_q};
}

}  // namespace detail

// Frequency-flat Rayleigh fading channel with configurable block coherence,
// CSI quality (perfect, MMSE noisy estimate, or blind/none), and M-ary
// Gray-coded modulation/demodulation.
//
// Model: y_t = h_t * x_t + w_t,
//   where h_t ~ CN(0, 1) (E[|h_t|^2] = 1, constant over `coherence_symbols`),
//         w_t ~ CN(0, N_0), N_0 = 10^(-snr_db / 10) for E_s = 1.0.
// Under CsiMode::kNoisyEstimate, h_hat_t = h_t + e_t with e_t ~ CN(0,
// csi_error_variance), and soft demapping uses the MMSE effective noise
// variance sigma_eff^2 = N_0 + csi_error_variance * E_s.
class RayleighFadingChannel {
 public:
  explicit RayleighFadingChannel(
      double snr_db, ModulationScheme scheme = ModulationScheme::kBpsk,
      std::size_t coherence_symbols = 1, CsiMode csi_mode = CsiMode::kPerfect,
      double csi_error_variance = 0.0, uint64_t seed = 42)
      : snr_db_(snr_db),
        noise_variance_n0_(
            detail::SnrDbToNoiseVarianceN0(snr_db, "RayleighFadingChannel")),
        scheme_(scheme),
        coherence_symbols_(coherence_symbols),
        csi_mode_(csi_mode),
        csi_error_variance_(csi_error_variance),
        modulator_(scheme),
        demodulator_(scheme),
        rng_(seed) {
    if (coherence_symbols_ == 0) {
      throw std::invalid_argument(
          "RayleighFadingChannel: coherence_symbols must be >= 1");
    }
    if (std::isnan(csi_error_variance_) || std::isinf(csi_error_variance_) ||
        csi_error_variance_ < 0.0) {
      throw std::invalid_argument(
          "RayleighFadingChannel: csi_error_variance must be finite and >= 0");
    }
  }

  // Exact analytical BER of uncoded coherent BPSK over Rayleigh fading:
  //   P_b = 0.5 * (1 - sqrt(gamma / (1 + gamma))), where gamma = 10^(Eb/N0 /
  //   10).
  static double AnalyticalRayleighBpskBer(double eb_n0_db) {
    if (std::isnan(eb_n0_db)) {
      throw std::invalid_argument(
          "AnalyticalRayleighBpskBer: eb_n0_db cannot be NaN");
    }
    if (std::isinf(eb_n0_db)) {
      return (eb_n0_db > 0.0) ? 0.0 : 0.5;
    }
    const double gamma = std::pow(10.0, eb_n0_db / 10.0);
    return 0.5 * (1.0 - std::sqrt(gamma / (1.0 + gamma)));
  }

  double SnrDb() const noexcept { return snr_db_; }
  double NoiseVarianceN0() const noexcept { return noise_variance_n0_; }
  ModulationScheme Scheme() const noexcept { return scheme_; }
  std::size_t CoherenceSymbols() const noexcept { return coherence_symbols_; }
  CsiMode GetCsiMode() const noexcept { return csi_mode_; }
  double CsiErrorVariance() const noexcept { return csi_error_variance_; }

  void SetSeed(uint64_t seed) noexcept { rng_.Seed(seed); }

  RayleighFadingChannel CloneWithSeed(uint64_t seed) const {
    return RayleighFadingChannel(snr_db_, scheme_, coherence_symbols_,
                                 csi_mode_, csi_error_variance_, seed);
  }

  FadingTransmissionResult Transmit(
      std::span<const uint8_t> bits,
      DemapMode demap_mode = DemapMode::kExactLogMap) {
    const std::vector<ComplexSymbol> tx = modulator_.Modulate(bits);
    const std::size_t num_symbols = tx.size();

    FadingTransmissionResult result;
    result.received_symbols.resize(num_symbols);
    result.channel_gains.resize(num_symbols);
    result.estimated_gains.resize(num_symbols);

    ComplexSymbol block_h{1.0, 0.0};
    ComplexSymbol block_h_est{1.0, 0.0};

    for (std::size_t t = 0; t < num_symbols; ++t) {
      if (t % coherence_symbols_ == 0) {
        block_h = detail::SampleComplexGaussian(rng_, 1.0);
        if (csi_mode_ == CsiMode::kPerfect) {
          block_h_est = block_h;
        } else if (csi_mode_ == CsiMode::kNoisyEstimate) {
          const ComplexSymbol est_err =
              detail::SampleComplexGaussian(rng_, csi_error_variance_);
          block_h_est = block_h + est_err;
        } else {
          block_h_est = ComplexSymbol{1.0, 0.0};
        }
      }

      const ComplexSymbol noise =
          detail::SampleComplexGaussian(rng_, noise_variance_n0_);
      result.channel_gains[t] = block_h;
      result.estimated_gains[t] = block_h_est;
      result.received_symbols[t] = block_h * tx[t] + noise;
    }

    // Effective noise variance for LLR computation (E_s = 1.0):
    //   kPerfect:       sigma_eff^2 = N_0
    //   kNoisyEstimate: sigma_eff^2 = N_0 + sigma_csi^2 * E_s
    //   kNone:          sigma_eff^2 = N_0 + 1.0
    double effective_n0 = noise_variance_n0_;
    if (csi_mode_ == CsiMode::kNoisyEstimate) {
      effective_n0 = noise_variance_n0_ + csi_error_variance_;
    } else if (csi_mode_ == CsiMode::kNone) {
      effective_n0 = noise_variance_n0_ + 1.0;
    }

    result.llrs =
        demodulator_.DemodulateSoft(result.received_symbols, effective_n0,
                                    demap_mode, result.estimated_gains);
    return result;
  }

 private:
  double snr_db_;
  double noise_variance_n0_;
  ModulationScheme scheme_;
  std::size_t coherence_symbols_;
  CsiMode csi_mode_;
  double csi_error_variance_;
  Modulator modulator_;
  Demodulator demodulator_;
  Pcg64 rng_;
};

// Frequency-flat Rician fading channel with unit average fading power:
//   h_t = sqrt(K / (K + 1)) + sqrt(1 / (K + 1)) * g_t,  g_t ~ CN(0, 1).
// K = 0 recovers exact Rayleigh fading; K -> +inf converges to AWGN.
class RicianFadingChannel {
 public:
  RicianFadingChannel(double snr_db, double k_factor_linear,
                      ModulationScheme scheme = ModulationScheme::kBpsk,
                      std::size_t coherence_symbols = 1,
                      CsiMode csi_mode = CsiMode::kPerfect,
                      double csi_error_variance = 0.0, uint64_t seed = 42)
      : snr_db_(snr_db),
        noise_variance_n0_(
            detail::SnrDbToNoiseVarianceN0(snr_db, "RicianFadingChannel")),
        k_factor_linear_(k_factor_linear),
        scheme_(scheme),
        coherence_symbols_(coherence_symbols),
        csi_mode_(csi_mode),
        csi_error_variance_(csi_error_variance),
        modulator_(scheme),
        demodulator_(scheme),
        rng_(seed) {
    if (std::isnan(k_factor_linear_) || std::isinf(k_factor_linear_) ||
        k_factor_linear_ < 0.0) {
      throw std::invalid_argument(
          "RicianFadingChannel: k_factor_linear must be finite and >= 0");
    }
    if (coherence_symbols_ == 0) {
      throw std::invalid_argument(
          "RicianFadingChannel: coherence_symbols must be >= 1");
    }
    if (std::isnan(csi_error_variance_) || std::isinf(csi_error_variance_) ||
        csi_error_variance_ < 0.0) {
      throw std::invalid_argument(
          "RicianFadingChannel: csi_error_variance must be finite and >= 0");
    }
    los_amp_ = std::sqrt(k_factor_linear_ / (k_factor_linear_ + 1.0));
    nlos_var_ = 1.0 / (k_factor_linear_ + 1.0);
  }

  static RicianFadingChannel FromKFactorDb(
      double snr_db, double k_factor_db,
      ModulationScheme scheme = ModulationScheme::kBpsk,
      std::size_t coherence_symbols = 1, CsiMode csi_mode = CsiMode::kPerfect,
      double csi_error_variance = 0.0, uint64_t seed = 42) {
    if (std::isnan(k_factor_db) || std::isinf(k_factor_db)) {
      throw std::invalid_argument(
          "RicianFadingChannel::FromKFactorDb: k_factor_db must be finite");
    }
    const double k_linear = std::pow(10.0, k_factor_db / 10.0);
    return RicianFadingChannel(snr_db, k_linear, scheme, coherence_symbols,
                               csi_mode, csi_error_variance, seed);
  }

  double SnrDb() const noexcept { return snr_db_; }
  double NoiseVarianceN0() const noexcept { return noise_variance_n0_; }
  double KFactorLinear() const noexcept { return k_factor_linear_; }
  double KFactorDb() const noexcept {
    if (k_factor_linear_ == 0.0) {
      return -std::numeric_limits<double>::infinity();
    }
    return 10.0 * std::log10(k_factor_linear_);
  }
  ModulationScheme Scheme() const noexcept { return scheme_; }
  std::size_t CoherenceSymbols() const noexcept { return coherence_symbols_; }
  CsiMode GetCsiMode() const noexcept { return csi_mode_; }
  double CsiErrorVariance() const noexcept { return csi_error_variance_; }

  void SetSeed(uint64_t seed) noexcept { rng_.Seed(seed); }

  RicianFadingChannel CloneWithSeed(uint64_t seed) const {
    return RicianFadingChannel(snr_db_, k_factor_linear_, scheme_,
                               coherence_symbols_, csi_mode_,
                               csi_error_variance_, seed);
  }

  FadingTransmissionResult Transmit(
      std::span<const uint8_t> bits,
      DemapMode demap_mode = DemapMode::kExactLogMap) {
    const std::vector<ComplexSymbol> tx = modulator_.Modulate(bits);
    const std::size_t num_symbols = tx.size();

    FadingTransmissionResult result;
    result.received_symbols.resize(num_symbols);
    result.channel_gains.resize(num_symbols);
    result.estimated_gains.resize(num_symbols);

    ComplexSymbol block_h{1.0, 0.0};
    ComplexSymbol block_h_est{1.0, 0.0};

    for (std::size_t t = 0; t < num_symbols; ++t) {
      if (t % coherence_symbols_ == 0) {
        const ComplexSymbol scatter =
            detail::SampleComplexGaussian(rng_, nlos_var_);
        block_h = ComplexSymbol{los_amp_ + scatter.i, scatter.q};

        if (csi_mode_ == CsiMode::kPerfect) {
          block_h_est = block_h;
        } else if (csi_mode_ == CsiMode::kNoisyEstimate) {
          const ComplexSymbol est_err =
              detail::SampleComplexGaussian(rng_, csi_error_variance_);
          block_h_est = block_h + est_err;
        } else {
          // Blind receiver: use deterministic LOS mean gain when K > 0, or
          // nominal unit gain when K == 0.
          block_h_est = ComplexSymbol{los_amp_ > 0.0 ? los_amp_ : 1.0, 0.0};
        }
      }

      const ComplexSymbol noise =
          detail::SampleComplexGaussian(rng_, noise_variance_n0_);
      result.channel_gains[t] = block_h;
      result.estimated_gains[t] = block_h_est;
      result.received_symbols[t] = block_h * tx[t] + noise;
    }

    double effective_n0 = noise_variance_n0_;
    if (csi_mode_ == CsiMode::kNoisyEstimate) {
      effective_n0 = noise_variance_n0_ + csi_error_variance_;
    } else if (csi_mode_ == CsiMode::kNone) {
      effective_n0 = noise_variance_n0_ + nlos_var_;
    }

    result.llrs =
        demodulator_.DemodulateSoft(result.received_symbols, effective_n0,
                                    demap_mode, result.estimated_gains);
    return result;
  }

 private:
  double snr_db_;
  double noise_variance_n0_;
  double k_factor_linear_;
  double los_amp_ = 0.0;
  double nlos_var_ = 1.0;
  ModulationScheme scheme_;
  std::size_t coherence_symbols_;
  CsiMode csi_mode_;
  double csi_error_variance_;
  Modulator modulator_;
  Demodulator demodulator_;
  Pcg64 rng_;
};

}  // namespace channel
}  // namespace xg

#endif  // XGALOIS_CHANNEL_FADING_CHANNEL_HPP_
