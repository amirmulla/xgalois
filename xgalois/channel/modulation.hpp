#ifndef XGALOIS_CHANNEL_MODULATION_HPP_
#define XGALOIS_CHANNEL_MODULATION_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace xg {
namespace channel {

enum class ModulationScheme : uint8_t {
  kBpsk,
  kQpsk,
  k8Psk,
  k16Qam,
  k64Qam,
};

enum class DemapMode : uint8_t {
  kExactLogMap,
  kMaxLogMap,
};

struct ComplexSymbol {
  double i = 0.0;
  double q = 0.0;

  constexpr double NormSq() const noexcept { return i * i + q * q; }
  double Abs() const noexcept { return std::hypot(i, q); }
  constexpr ComplexSymbol Conj() const noexcept { return {i, -q}; }

  constexpr ComplexSymbol operator+(const ComplexSymbol& rhs) const noexcept {
    return {i + rhs.i, q + rhs.q};
  }
  constexpr ComplexSymbol operator-(const ComplexSymbol& rhs) const noexcept {
    return {i - rhs.i, q - rhs.q};
  }
  constexpr ComplexSymbol operator*(const ComplexSymbol& rhs) const noexcept {
    return {i * rhs.i - q * rhs.q, i * rhs.q + q * rhs.i};
  }
  constexpr ComplexSymbol operator*(double s) const noexcept {
    return {i * s, q * s};
  }
  constexpr bool operator==(const ComplexSymbol& rhs) const noexcept = default;
};

inline constexpr ComplexSymbol operator*(double s,
                                         const ComplexSymbol& c) noexcept {
  return {s * c.i, s * c.q};
}

namespace detail {

inline constexpr std::size_t SchemeBitsPerSymbol(
    ModulationScheme scheme) noexcept {
  switch (scheme) {
    case ModulationScheme::kBpsk:
      return 1;
    case ModulationScheme::kQpsk:
      return 2;
    case ModulationScheme::k8Psk:
      return 3;
    case ModulationScheme::k16Qam:
      return 4;
    case ModulationScheme::k64Qam:
      return 6;
  }
  return 1;
}

// 1-D Gray-coded PAM levels indexed by the m_1d-bit integer label (MSB first).
// Every pair of geometrically adjacent levels differs in exactly 1 bit, and
// MSB == 0 corresponds to the positive half-axis.
inline std::vector<double> BuildGrayPamLevels(std::size_t bits_per_axis,
                                              double scale) {
  const std::size_t num_levels = std::size_t{1} << bits_per_axis;
  std::vector<double> levels(num_levels, 0.0);
  for (std::size_t rank = 0; rank < num_levels; ++rank) {
    // Rank 0 -> largest positive level +(num_levels - 1);
    // Rank (num_levels - 1) -> largest negative level -(num_levels - 1).
    const std::size_t gray_label = rank ^ (rank >> 1);
    const double unnorm =
        static_cast<double>(num_levels - 1) - 2.0 * static_cast<double>(rank);
    levels[gray_label] = unnorm * scale;
  }
  return levels;
}

// Builds the unit-average-energy (E_s = 1.0) Gray-coded constellation table
// indexed by the big-endian m-bit label (bit 0 = MSB).
inline std::vector<ComplexSymbol> BuildConstellation(ModulationScheme scheme) {
  switch (scheme) {
    case ModulationScheme::kBpsk:
      // Matches XGalois AWGNChannel convention: bit 0 -> +1.0, bit 1 -> -1.0.
      return {ComplexSymbol{1.0, 0.0}, ComplexSymbol{-1.0, 0.0}};

    case ModulationScheme::kQpsk: {
      const double s = 1.0 / std::sqrt(2.0);
      const auto pam = BuildGrayPamLevels(1, s);
      std::vector<ComplexSymbol> table(4);
      for (std::size_t i_bits = 0; i_bits < 2; ++i_bits) {
        for (std::size_t q_bits = 0; q_bits < 2; ++q_bits) {
          table[(i_bits << 1) | q_bits] =
              ComplexSymbol{pam[i_bits], pam[q_bits]};
        }
      }
      return table;
    }

    case ModulationScheme::k8Psk: {
      // Binary-reflected Gray code around the unit circle: phase k has
      // 3-bit Gray label k ^ (k >> 1), so adjacent points (including wrap
      // 7 -> 0) differ in exactly 1 bit and every symbol has |s|^2 = 1.0.
      std::vector<ComplexSymbol> table(8);
      for (std::size_t k = 0; k < 8; ++k) {
        const std::size_t gray_label = k ^ (k >> 1);
        const double angle =
            (2.0 * std::numbers::pi * static_cast<double>(k)) / 8.0;
        table[gray_label] = ComplexSymbol{std::cos(angle), std::sin(angle)};
      }
      return table;
    }

    case ModulationScheme::k16Qam: {
      // Unnormalized levels {-3, -1, +1, +3} per axis have E[x_I^2] = 5,
      // so E[|x|^2] = 10 -> scale = 1 / sqrt(10).
      const double s = 1.0 / std::sqrt(10.0);
      const auto pam = BuildGrayPamLevels(2, s);
      std::vector<ComplexSymbol> table(16);
      for (std::size_t i_bits = 0; i_bits < 4; ++i_bits) {
        for (std::size_t q_bits = 0; q_bits < 4; ++q_bits) {
          table[(i_bits << 2) | q_bits] =
              ComplexSymbol{pam[i_bits], pam[q_bits]};
        }
      }
      return table;
    }

    case ModulationScheme::k64Qam: {
      // Unnormalized levels {-7, -5, -3, -1, +1, +3, +5, +7} per axis have
      // E[x_I^2] = 21, so E[|x|^2] = 42 -> scale = 1 / sqrt(42).
      const double s = 1.0 / std::sqrt(42.0);
      const auto pam = BuildGrayPamLevels(3, s);
      std::vector<ComplexSymbol> table(64);
      for (std::size_t i_bits = 0; i_bits < 8; ++i_bits) {
        for (std::size_t q_bits = 0; q_bits < 8; ++q_bits) {
          table[(i_bits << 3) | q_bits] =
              ComplexSymbol{pam[i_bits], pam[q_bits]};
        }
      }
      return table;
    }
  }
  throw std::invalid_argument("Unsupported ModulationScheme");
}

inline double LogSumExp(std::span<const double> metrics) noexcept {
  double max_val = metrics[0];
  for (std::size_t i = 1; i < metrics.size(); ++i) {
    if (metrics[i] > max_val) {
      max_val = metrics[i];
    }
  }
  if (!std::isfinite(max_val)) {
    return max_val;
  }
  double sum = 0.0;
  for (double v : metrics) {
    sum += std::exp(v - max_val);
  }
  return max_val + std::log(sum);
}

inline double MaxMetric(std::span<const double> metrics) noexcept {
  double max_val = metrics[0];
  for (std::size_t i = 1; i < metrics.size(); ++i) {
    if (metrics[i] > max_val) {
      max_val = metrics[i];
    }
  }
  return max_val;
}

}  // namespace detail

// Gray-coded baseband modulator normalized to exact unit average symbol energy
// E_s = 1.0 for BPSK, QPSK, 8-PSK, 16-QAM, and 64-QAM.
class Modulator {
 public:
  explicit Modulator(ModulationScheme scheme)
      : scheme_(scheme),
        bits_per_symbol_(detail::SchemeBitsPerSymbol(scheme)),
        constellation_(detail::BuildConstellation(scheme)) {}

  ModulationScheme Scheme() const noexcept { return scheme_; }
  std::size_t BitsPerSymbol() const noexcept { return bits_per_symbol_; }

  const std::vector<ComplexSymbol>& Constellation() const noexcept {
    return constellation_;
  }

  std::vector<ComplexSymbol> Modulate(std::span<const uint8_t> bits) const {
    if (bits.size() % bits_per_symbol_ != 0) {
      throw std::invalid_argument(
          "Modulator::Modulate: input bit length must be a multiple of "
          "BitsPerSymbol()");
    }
    const std::size_t num_symbols = bits.size() / bits_per_symbol_;
    std::vector<ComplexSymbol> symbols(num_symbols);
    for (std::size_t s = 0; s < num_symbols; ++s) {
      std::size_t label = 0;
      const std::size_t base = s * bits_per_symbol_;
      for (std::size_t b = 0; b < bits_per_symbol_; ++b) {
        const uint8_t bit = bits[base + b];
        if (bit > 1) {
          throw std::invalid_argument(
              "Modulator::Modulate: bit values must be 0 or 1");
        }
        label = (label << 1) | static_cast<std::size_t>(bit);
      }
      symbols[s] = constellation_[label];
    }
    return symbols;
  }

 private:
  ModulationScheme scheme_;
  std::size_t bits_per_symbol_;
  std::vector<ComplexSymbol> constellation_;
};

// Coherent hard-decision and soft-decision (Exact Log-MAP / Max-Log-MAP)
// baseband demodulator following XGalois's LLR sign convention:
//   L(b) = ln(P(b = 0 | y, h) / P(b = 1 | y, h))
//   (positive => bit 0, negative => bit 1).
//
// Noise convention: complex circularly symmetric AWGN w ~ CN(0, N_0), so
// E[|w|^2] = N_0 and each real dimension (I and Q) has variance N_0 / 2.
//
// Uses O(sqrt(M)) separable 1-D I/Q demapping for QPSK, 16-QAM, and 64-QAM
// (algebraically identical to full 2-D Bayes demapping) and exact log-sum-exp
// over all 8 symbols for 8-PSK.
class Demodulator {
 public:
  explicit Demodulator(ModulationScheme scheme)
      : scheme_(scheme),
        bits_per_symbol_(detail::SchemeBitsPerSymbol(scheme)),
        constellation_(detail::BuildConstellation(scheme)) {
    InitSeparablePam();
  }

  ModulationScheme Scheme() const noexcept { return scheme_; }
  std::size_t BitsPerSymbol() const noexcept { return bits_per_symbol_; }

  const std::vector<ComplexSymbol>& Constellation() const noexcept {
    return constellation_;
  }

  // Coherent minimum-Euclidean-distance hard decision.
  // If `channel_gains` is empty, unit gain h = (1, 0) is assumed.
  std::vector<uint8_t> DemodulateHard(
      std::span<const ComplexSymbol> rx,
      std::span<const ComplexSymbol> channel_gains = {}) const {
    ValidateGains(rx.size(), channel_gains);
    std::vector<uint8_t> bits(rx.size() * bits_per_symbol_, 0);
    const bool has_gains = !channel_gains.empty();

    for (std::size_t s = 0; s < rx.size(); ++s) {
      const ComplexSymbol h =
          has_gains ? channel_gains[s] : ComplexSymbol{1.0, 0.0};
      const double h_norm_sq = h.NormSq();
      std::size_t best_label = 0;

      if (h_norm_sq > 0.0) {
        const ComplexSymbol z = (h.Conj() * rx[s]) * (1.0 / h_norm_sq);
        if (is_separable_) {
          const std::size_t best_i = NearestPamLabel(z.i);
          if (scheme_ == ModulationScheme::kBpsk) {
            best_label = best_i;
          } else {
            const std::size_t best_q = NearestPamLabel(z.q);
            best_label = (best_i << bits_per_axis_) | best_q;
          }
        } else {
          double best_dist = std::numeric_limits<double>::infinity();
          for (std::size_t idx = 0; idx < constellation_.size(); ++idx) {
            const double d = (z - constellation_[idx]).NormSq();
            if (d < best_dist) {
              best_dist = d;
              best_label = idx;
            }
          }
        }
      }

      const std::size_t base = s * bits_per_symbol_;
      for (std::size_t b = 0; b < bits_per_symbol_; ++b) {
        const std::size_t shift = bits_per_symbol_ - 1 - b;
        bits[base + b] = static_cast<uint8_t>((best_label >> shift) & 1U);
      }
    }

    return bits;
  }

  // Coherent soft-decision LLR demapper.
  // `noise_variance_n0` is the complex noise variance N_0 = E[|w|^2] =
  // 2*sigma^2. If `channel_gains` is empty, unit gain h = (1, 0) is assumed.
  std::vector<double> DemodulateSoft(
      std::span<const ComplexSymbol> rx, double noise_variance_n0,
      DemapMode mode = DemapMode::kExactLogMap,
      std::span<const ComplexSymbol> channel_gains = {}) const {
    if (std::isnan(noise_variance_n0) || noise_variance_n0 < 0.0) {
      throw std::invalid_argument(
          "Demodulator::DemodulateSoft: noise_variance_n0 must be finite and "
          "non-negative");
    }
    ValidateGains(rx.size(), channel_gains);

    std::vector<double> llrs(rx.size() * bits_per_symbol_, 0.0);
    const bool has_gains = !channel_gains.empty();

    for (std::size_t s = 0; s < rx.size(); ++s) {
      const ComplexSymbol h =
          has_gains ? channel_gains[s] : ComplexSymbol{1.0, 0.0};
      const double h_norm_sq = h.NormSq();
      double* out = llrs.data() + s * bits_per_symbol_;

      if (h_norm_sq == 0.0 || std::isinf(noise_variance_n0)) {
        for (std::size_t b = 0; b < bits_per_symbol_; ++b) {
          out[b] = 0.0;
        }
        continue;
      }

      // Zero-noise limit: deterministic +/-inf LLRs from minimum-distance
      // decision (0.0 on exact ties).
      if (noise_variance_n0 == 0.0) {
        DemapZeroNoise(rx[s], h, out);
        continue;
      }

      // Sufficient statistic z = (h* y) / |h|^2 = x + w_tilde, where
      // |y - h x|^2 / N_0 = (|h|^2 / N_0) * |z - x|^2.
      const ComplexSymbol z = (h.Conj() * rx[s]) * (1.0 / h_norm_sq);
      const double inv_eff_n0 = h_norm_sq / noise_variance_n0;

      if (scheme_ == ModulationScheme::kBpsk) {
        // Exact BPSK LLR: (|z + 1|^2 - |z - 1|^2) * |h|^2 / N_0
        //               = 4 * Re(h* y) / N_0.
        out[0] = 4.0 * z.i * inv_eff_n0;
      } else if (is_separable_) {
        DemapAxis1D(z.i, inv_eff_n0, mode, out);
        DemapAxis1D(z.q, inv_eff_n0, mode, out + bits_per_axis_);
      } else {
        Demap2DExact(z, inv_eff_n0, mode, out);
      }
    }

    return llrs;
  }

  // Convenience overload accepting `channel_gains` directly after `rx`.
  std::vector<double> DemodulateSoft(
      std::span<const ComplexSymbol> rx,
      std::span<const ComplexSymbol> channel_gains, double noise_variance_n0,
      DemapMode mode = DemapMode::kExactLogMap) const {
    return DemodulateSoft(rx, noise_variance_n0, mode, channel_gains);
  }

 private:
  void InitSeparablePam() {
    switch (scheme_) {
      case ModulationScheme::kBpsk:
        is_separable_ = true;
        bits_per_axis_ = 1;
        pam_levels_ = {1.0, -1.0};
        break;
      case ModulationScheme::kQpsk:
        is_separable_ = true;
        bits_per_axis_ = 1;
        pam_levels_ = detail::BuildGrayPamLevels(1, 1.0 / std::sqrt(2.0));
        break;
      case ModulationScheme::k8Psk:
        is_separable_ = false;
        bits_per_axis_ = 0;
        break;
      case ModulationScheme::k16Qam:
        is_separable_ = true;
        bits_per_axis_ = 2;
        pam_levels_ = detail::BuildGrayPamLevels(2, 1.0 / std::sqrt(10.0));
        break;
      case ModulationScheme::k64Qam:
        is_separable_ = true;
        bits_per_axis_ = 3;
        pam_levels_ = detail::BuildGrayPamLevels(3, 1.0 / std::sqrt(42.0));
        break;
    }
  }

  static void ValidateGains(std::size_t num_rx,
                            std::span<const ComplexSymbol> channel_gains) {
    if (!channel_gains.empty() && channel_gains.size() != num_rx) {
      throw std::invalid_argument(
          "Demodulator: channel_gains size must match rx size");
    }
  }

  std::size_t NearestPamLabel(double coord) const noexcept {
    std::size_t best_label = 0;
    double best_dist = std::numeric_limits<double>::infinity();
    for (std::size_t idx = 0; idx < pam_levels_.size(); ++idx) {
      const double diff = coord - pam_levels_[idx];
      const double d2 = diff * diff;
      if (d2 < best_dist) {
        best_dist = d2;
        best_label = idx;
      }
    }
    return best_label;
  }

  // O(sqrt(M)) 1-D axis demapping: computes LLR for each of the
  // `bits_per_axis_` bits controlling a single real axis (I or Q).
  void DemapAxis1D(double coord, double inv_eff_n0, DemapMode mode,
                   double* out_llrs) const noexcept {
    // Maximum number of 1-D levels is 8 (for 64-QAM).
    std::array<double, 8> metrics{};
    const std::size_t num_levels = pam_levels_.size();
    for (std::size_t idx = 0; idx < num_levels; ++idx) {
      const double diff = coord - pam_levels_[idx];
      metrics[idx] = -(diff * diff) * inv_eff_n0;
    }

    const std::size_t half = num_levels / 2;
    std::array<double, 4> m0{};
    std::array<double, 4> m1{};

    for (std::size_t b = 0; b < bits_per_axis_; ++b) {
      const std::size_t bit_mask = std::size_t{1} << (bits_per_axis_ - 1 - b);
      std::size_t c0 = 0;
      std::size_t c1 = 0;
      for (std::size_t idx = 0; idx < num_levels; ++idx) {
        if ((idx & bit_mask) == 0) {
          m0[c0++] = metrics[idx];
        } else {
          m1[c1++] = metrics[idx];
        }
      }
      const std::span<const double> s0(m0.data(), half);
      const std::span<const double> s1(m1.data(), half);
      if (mode == DemapMode::kExactLogMap) {
        out_llrs[b] = detail::LogSumExp(s0) - detail::LogSumExp(s1);
      } else {
        out_llrs[b] = detail::MaxMetric(s0) - detail::MaxMetric(s1);
      }
    }
  }

  // Full 2-D demapping (used for non-separable constellations such as 8-PSK).
  void Demap2DExact(ComplexSymbol z, double inv_eff_n0, DemapMode mode,
                    double* out_llrs) const noexcept {
    const std::size_t m_size = constellation_.size();
    std::array<double, 64> metrics{};
    for (std::size_t idx = 0; idx < m_size; ++idx) {
      metrics[idx] = -(z - constellation_[idx]).NormSq() * inv_eff_n0;
    }

    const std::size_t half = m_size / 2;
    std::array<double, 32> m0{};
    std::array<double, 32> m1{};

    for (std::size_t b = 0; b < bits_per_symbol_; ++b) {
      const std::size_t bit_mask = std::size_t{1} << (bits_per_symbol_ - 1 - b);
      std::size_t c0 = 0;
      std::size_t c1 = 0;
      for (std::size_t idx = 0; idx < m_size; ++idx) {
        if ((idx & bit_mask) == 0) {
          m0[c0++] = metrics[idx];
        } else {
          m1[c1++] = metrics[idx];
        }
      }
      const std::span<const double> s0(m0.data(), half);
      const std::span<const double> s1(m1.data(), half);
      if (mode == DemapMode::kExactLogMap) {
        out_llrs[b] = detail::LogSumExp(s0) - detail::LogSumExp(s1);
      } else {
        out_llrs[b] = detail::MaxMetric(s0) - detail::MaxMetric(s1);
      }
    }
  }

  void DemapZeroNoise(ComplexSymbol y, ComplexSymbol h,
                      double* out_llrs) const noexcept {
    const double h_norm_sq = h.NormSq();
    const ComplexSymbol z = (h.Conj() * y) * (1.0 / h_norm_sq);
    const std::size_t m_size = constellation_.size();
    const double inf = std::numeric_limits<double>::infinity();

    for (std::size_t b = 0; b < bits_per_symbol_; ++b) {
      const std::size_t bit_mask = std::size_t{1} << (bits_per_symbol_ - 1 - b);
      double min_d0 = inf;
      double min_d1 = inf;
      for (std::size_t idx = 0; idx < m_size; ++idx) {
        const double d = (z - constellation_[idx]).NormSq();
        if ((idx & bit_mask) == 0) {
          min_d0 = std::min(min_d0, d);
        } else {
          min_d1 = std::min(min_d1, d);
        }
      }
      if (min_d0 < min_d1) {
        out_llrs[b] = inf;
      } else if (min_d1 < min_d0) {
        out_llrs[b] = -inf;
      } else {
        out_llrs[b] = 0.0;
      }
    }
  }

  ModulationScheme scheme_;
  std::size_t bits_per_symbol_;
  std::vector<ComplexSymbol> constellation_;
  bool is_separable_ = false;
  std::size_t bits_per_axis_ = 0;
  std::vector<double> pam_levels_;
};

}  // namespace channel
}  // namespace xg

#endif  // XGALOIS_CHANNEL_MODULATION_HPP_
