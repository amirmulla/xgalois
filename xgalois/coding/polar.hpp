#ifndef XGALOIS_CODING_POLAR_HPP_
#define XGALOIS_CODING_POLAR_HPP_

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

namespace polar_detail {

constexpr double kPolarMaxInputLlr = 50.0;
constexpr double kPolarMaxStageLlr = 1000.0;

// Exact Jacobian-logarithm box-plus for polar f(a, b).
inline double ExactBoxPlus(double a, double b) noexcept {
  const double abs_a = std::abs(a);
  const double abs_b = std::abs(b);
  const double min_mag = (abs_a < abs_b) ? abs_a : abs_b;
  const double sign = ((a >= 0.0) == (b >= 0.0)) ? 1.0 : -1.0;
  const double corr_pos = std::log1p(std::exp(-std::abs(a + b)));
  const double corr_neg = std::log1p(std::exp(-std::abs(a - b)));
  return std::clamp(sign * min_mag + corr_pos - corr_neg, -kPolarMaxStageLlr,
                    kPolarMaxStageLlr);
}

// Min-sum approximation for polar f(a, b) = sign(a)*sign(b)*min(|a|, |b|).
inline double MinSumF(double a, double b) noexcept {
  const double abs_a = std::abs(a);
  const double abs_b = std::abs(b);
  const double min_mag = (abs_a < abs_b) ? abs_a : abs_b;
  const double sign = ((a >= 0.0) == (b >= 0.0)) ? 1.0 : -1.0;
  return std::clamp(sign * min_mag, -kPolarMaxStageLlr, kPolarMaxStageLlr);
}

// Polar g(a, b, u) = (1 - 2u) * a + b.
inline double PolarG(double a, double b, uint8_t u) noexcept {
  const double val = ((u & 1U) == 0) ? (a + b) : (b - a);
  return std::clamp(val, -kPolarMaxStageLlr, kPolarMaxStageLlr);
}

// LLR path metric penalty: PM += |L| whenever decision disagrees with sign(L).
inline double PathMetricPenalty(double llr, uint8_t bit) noexcept {
  const uint8_t hard = (llr >= 0.0) ? uint8_t{0} : uint8_t{1};
  return ((bit & 1U) == hard) ? 0.0 : std::abs(llr);
}

// Continuous, strictly decreasing log(phi(x)) approximation for Gaussian
// Approximation (Chung et al. / Trifonov) of BPSK-AWGN channel polarization.
inline double LogPhi(double x) noexcept {
  if (x <= 0.0) {
    return 0.0;
  }
  if (x < 10.0) {
    return -0.4527 * std::pow(x, 0.86);
  }
  // Matched continuously at x = 10 to -0.4527 * 10^0.86.
  constexpr double kMatchOffset = -0.046582721911887;
  return -0.25 * x + 0.5 * std::log(std::numbers::pi / x) +
         std::log(1.0 - 10.0 / (7.0 * x)) + kMatchOffset;
}

// Monotone bisection inverse of LogPhi(x) for lp <= 0.
inline double InvLogPhi(double lp) noexcept {
  if (lp >= 0.0) {
    return 0.0;
  }
  double lo = 0.0;
  double hi = std::max(1.0, -4.0 * lp + 20.0);
  while (LogPhi(hi) > lp && hi < 1e12) {
    hi *= 2.0;
  }
  for (int iter = 0; iter < 60; ++iter) {
    const double mid = 0.5 * (lo + hi);
    if (LogPhi(mid) > lp) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

inline uint64_t DefaultCrcPolynomial(std::size_t crc_bits) {
  switch (crc_bits) {
    case 4:
      return 0x13ULL;  // x^4 + x + 1
    case 6:
      return 0x61ULL;  // 3GPP NR CRC-6: x^6 + x^5 + 1
    case 8:
      return 0x107ULL;  // CRC-8-CCITT: x^8 + x^2 + x + 1
    case 11:
      return 0xE21ULL;  // 3GPP NR CRC-11: x^11 + x^10 + x^9 + x^5 + 1
    case 16:
      return 0x11021ULL;  // CRC-16-CCITT: x^16 + x^12 + x^5 + 1
    case 24:
      return 0x1864CFBULL;  // CRC-24A
    default:
      return (uint64_t{1} << crc_bits) | 3ULL;
  }
}

inline uint64_t NormalizeCrcPolynomial(std::size_t crc_bits,
                                       uint64_t crc_poly) {
  if (crc_bits == 0) {
    if (crc_poly != 0) {
      throw std::invalid_argument(
          "PolarCode: crc_poly must be 0 when crc_bits == 0");
    }
    return 0;
  }
  if (crc_bits > 32) {
    throw std::invalid_argument("PolarCode: crc_bits must be between 0 and 32");
  }
  if (crc_poly == 0) {
    return DefaultCrcPolynomial(crc_bits);
  }
  const uint64_t top_bit = uint64_t{1} << crc_bits;
  const uint64_t max_allowed = (top_bit << 1) - 1ULL;
  if (crc_poly > max_allowed) {
    throw std::invalid_argument("PolarCode: crc_poly degree exceeds crc_bits");
  }
  return crc_poly | top_bit;
}

inline void ValidatePolarLlr(const std::vector<double> &channel_llrs,
                             std::size_t expected_length, const char *fn_name) {
  if (channel_llrs.size() != expected_length) {
    throw std::invalid_argument(std::string(fn_name) +
                                ": channel LLR length (" +
                                std::to_string(channel_llrs.size()) +
                                ") does not match block length (" +
                                std::to_string(expected_length) + ")");
  }
  for (std::size_t i = 0; i < expected_length; ++i) {
    if (std::isnan(channel_llrs[i])) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": channel LLRs must not contain NaN");
    }
  }
}

}  // namespace polar_detail

class PolarCode;
class PolarEncoder;
class SuccessiveCancellationDecoder;
class SuccessiveCancellationListDecoder;

// Binary Arikan Polar Code over GF(2) with block length N = 2^n, message
// dimension K, non-systematic or systematic encoding, and optional CRC-aided
// polar coding (CA-Polar) for CA-SCL decoding.
//
// Transform convention: G_N = F^{\otimes n} where F = [[1, 0], [1, 1]] (natural
// order, self-inverse G_N^2 = I_N).
class PolarCode : public AbstractLinearCode<GaloisFieldBinary> {
 public:
  using Base = AbstractLinearCode<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using matrix_type = GFMatrix<GaloisFieldBinary>;
  using syndrome_type = GFArray<GaloisFieldBinary>;

  // In-place O(N log N) Arikan butterfly polar transform: x = u * F^{\otimes
  // n}. Since F^2 = I_2 over GF(2), PolarTransform is an involution
  // (self-inverse).
  static void PolarTransform(std::span<uint8_t> data) {
    const std::size_t n = data.size();
    if (n == 0 || !std::has_single_bit(n)) {
      throw std::invalid_argument(
          "PolarCode::PolarTransform: input length must be a power of 2");
    }
    for (uint8_t &v : data) {
      v = static_cast<uint8_t>(v & 1U);
    }
    for (std::size_t step = 1; step < n; step <<= 1) {
      for (std::size_t i = 0; i < n; i += (step << 1)) {
        for (std::size_t j = 0; j < step; ++j) {
          data[i + j] = static_cast<uint8_t>(data[i + j] ^ data[i + j + step]);
        }
      }
    }
  }

  static GFArray<GaloisFieldBinary> PolarTransform(
      const GFArray<GaloisFieldBinary> &input) {
    if (!input.HasField()) {
      throw std::invalid_argument(
          "PolarCode::PolarTransform: input GFArray has null field");
    }
    if (input.Dimension() != 1 || input.Empty() ||
        !std::has_single_bit(input.Size())) {
      throw std::invalid_argument(
          "PolarCode::PolarTransform: input must be 1D with power-of-2 length");
    }
    GFArray<GaloisFieldBinary> out = input;
    PolarTransform(out.RawSpan());
    return out;
  }

  // Construct frozen bit mask (true = frozen, false = info) of length N = 2^n
  // with K info bits using log-domain Bhattacharyya parameter recursion:
  //   Z(W^-) = 2Z - Z^2,   Z(W^+) = Z^2.
  static std::vector<bool> ConstructBhattacharyya(
      std::size_t n, std::size_t k, double design_erasure_prob = 0.5) {
    ValidateBlockAndDimension(n, k, "PolarCode::ConstructBhattacharyya");
    if (!std::isfinite(design_erasure_prob) || design_erasure_prob <= 0.0 ||
        design_erasure_prob >= 1.0) {
      throw std::invalid_argument(
          "PolarCode::ConstructBhattacharyya: design_erasure_prob must be in "
          "(0, 1)");
    }

    const std::size_t num_stages =
        static_cast<std::size_t>(std::countr_zero(n));
    std::vector<double> log_z{std::log(design_erasure_prob)};
    std::vector<double> next;
    for (std::size_t stage = 0; stage < num_stages; ++stage) {
      next.resize(log_z.size() * 2);
      for (std::size_t i = 0; i < log_z.size(); ++i) {
        const double lz = log_z[i];
        // ln(2Z - Z^2) = ln(Z) + ln(2 - Z) = lz + log1p(1 - exp(lz))
        //              = lz + log1p(-expm1(lz)).
        next[2 * i] = std::min(0.0, lz + std::log1p(-std::expm1(lz)));
        next[2 * i + 1] = 2.0 * lz;
      }
      log_z.swap(next);
    }

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(
        order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
          if (log_z[a] != log_z[b]) {
            return log_z[a] <
                   log_z[b];  // smaller Bhattacharyya Z = more reliable
          }
          const int pa = std::popcount(a);
          const int pb = std::popcount(b);
          if (pa != pb) {
            return pa > pb;
          }
          return a > b;
        });

    std::vector<bool> frozen_mask(n, true);
    for (std::size_t r = 0; r < k; ++r) {
      frozen_mask[order[r]] = false;
    }
    return frozen_mask;
  }

  // Construct frozen bit mask (true = frozen, false = info) of length N = 2^n
  // with K info bits using Gaussian Approximation at `design_snr_db` (Eb/N0 in
  // dB, matching AWGNChannel's variance sigma^2 = 1 / (2 * R * 10^(snr/10))).
  static std::vector<bool> ConstructDesignSnr(std::size_t n, std::size_t k,
                                              double design_snr_db = 2.0) {
    ValidateBlockAndDimension(n, k, "PolarCode::ConstructDesignSnr");
    if (!std::isfinite(design_snr_db)) {
      throw std::invalid_argument(
          "PolarCode::ConstructDesignSnr: design_snr_db must be finite");
    }

    const double rate = static_cast<double>(k) / static_cast<double>(n);
    const double mu0 = 4.0 * rate * std::pow(10.0, design_snr_db / 10.0);
    const std::size_t num_stages =
        static_cast<std::size_t>(std::countr_zero(n));

    std::vector<double> mu{mu0};
    std::vector<double> next;
    for (std::size_t stage = 0; stage < num_stages; ++stage) {
      next.resize(mu.size() * 2);
      for (std::size_t i = 0; i < mu.size(); ++i) {
        const double lp = polar_detail::LogPhi(mu[i]);
        const double lp_minus = std::min(0.0, lp + std::log1p(-std::expm1(lp)));
        next[2 * i] = polar_detail::InvLogPhi(lp_minus);
        next[2 * i + 1] = 2.0 * mu[i];
      }
      mu.swap(next);
    }

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(
        order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
          if (mu[a] != mu[b]) {
            return mu[a] > mu[b];  // larger mean LLR = more reliable
          }
          const int pa = std::popcount(a);
          const int pb = std::popcount(b);
          if (pa != pb) {
            return pa > pb;
          }
          return a > b;
        });

    std::vector<bool> frozen_mask(n, true);
    for (std::size_t r = 0; r < k; ++r) {
      frozen_mask[order[r]] = false;
    }
    return frozen_mask;
  }

  // Construct a PolarCode from an explicit frozen mask (true = frozen bit,
  // false = non-frozen bit). When `crc_bits > 0`, the message dimension is
  // `num_non_frozen - crc_bits` (must be >= 1).
  PolarCode(std::shared_ptr<GaloisFieldBinary> field, std::size_t block_length,
            std::vector<bool> frozen_mask, bool systematic = false,
            std::size_t crc_bits = 0, uint64_t crc_poly = 0)
      : Base(std::move(field),
             ValidatePowerOfTwoLength(block_length, "PolarCode"),
             ComputeMessageDimension(block_length, frozen_mask, crc_bits),
             Metric::HAMMING, "Polar", "SC"),
        num_stages_(static_cast<std::size_t>(std::countr_zero(block_length))),
        frozen_mask_(std::move(frozen_mask)),
        systematic_(systematic),
        crc_bits_(crc_bits),
        crc_poly_(polar_detail::NormalizeCrcPolynomial(crc_bits, crc_poly)) {
    InitializeStructure();
  }

  // Convenience constructor using Gaussian Approximation at `design_snr_db`.
  PolarCode(std::shared_ptr<GaloisFieldBinary> field, std::size_t block_length,
            std::size_t message_dimension, double design_snr_db = 2.0,
            bool systematic = false, std::size_t crc_bits = 0,
            uint64_t crc_poly = 0)
      : PolarCode(std::move(field), block_length,
                  ConstructDesignSnr(
                      block_length,
                      ValidateMessageWithCrc(block_length, message_dimension,
                                             crc_bits, "PolarCode"),
                      design_snr_db),
                  systematic, crc_bits, crc_poly) {}

  static PolarCode FromBhattacharyya(std::shared_ptr<GaloisFieldBinary> field,
                                     std::size_t block_length,
                                     std::size_t message_dimension,
                                     double design_erasure_prob = 0.5,
                                     bool systematic = false,
                                     std::size_t crc_bits = 0,
                                     uint64_t crc_poly = 0) {
    const std::size_t k_polar =
        ValidateMessageWithCrc(block_length, message_dimension, crc_bits,
                               "PolarCode::FromBhattacharyya");
    auto mask =
        ConstructBhattacharyya(block_length, k_polar, design_erasure_prob);
    return PolarCode(std::move(field), block_length, std::move(mask),
                     systematic, crc_bits, crc_poly);
  }

  static PolarCode FromDesignSnr(std::shared_ptr<GaloisFieldBinary> field,
                                 std::size_t block_length,
                                 std::size_t message_dimension,
                                 double design_snr_db = 2.0,
                                 bool systematic = false,
                                 std::size_t crc_bits = 0,
                                 uint64_t crc_poly = 0) {
    const std::size_t k_polar = ValidateMessageWithCrc(
        block_length, message_dimension, crc_bits, "PolarCode::FromDesignSnr");
    auto mask = ConstructDesignSnr(block_length, k_polar, design_snr_db);
    return PolarCode(std::move(field), block_length, std::move(mask),
                     systematic, crc_bits, crc_poly);
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  std::size_t NumStages() const noexcept { return num_stages_; }
  bool IsSystematicPolar() const noexcept { return systematic_; }
  const std::vector<bool> &FrozenMask() const noexcept { return frozen_mask_; }
  const std::vector<std::size_t> &InfoPositions() const noexcept {
    return info_positions_;
  }
  const std::vector<std::size_t> &FrozenPositions() const noexcept {
    return frozen_positions_;
  }
  std::size_t CrcBits() const noexcept { return crc_bits_; }
  uint64_t CrcPoly() const noexcept { return crc_poly_; }

  // Linear binary CRC remainder m(x) * x^r mod g(x) with zero initial state.
  uint64_t ComputeCrc(std::span<const uint8_t> message_bits) const noexcept {
    if (crc_bits_ == 0) {
      return 0;
    }
    const uint64_t top_bit = uint64_t{1} << crc_bits_;
    const uint64_t mask = top_bit - 1ULL;
    uint64_t rem = 0;
    for (uint8_t bit : message_bits) {
      rem = (rem << 1) | static_cast<uint64_t>(bit & 1U);
      if ((rem & top_bit) != 0) {
        rem ^= crc_poly_;
      }
    }
    for (std::size_t i = 0; i < crc_bits_; ++i) {
      rem <<= 1;
      if ((rem & top_bit) != 0) {
        rem ^= crc_poly_;
      }
    }
    return rem & mask;
  }

  // Verify CRC on a payload vector [message (K bits) | crc (crc_bits bits)].
  bool VerifyCrc(std::span<const uint8_t> payload_bits) const noexcept {
    if (crc_bits_ == 0) {
      return true;
    }
    if (payload_bits.size() != this->dimension_ + crc_bits_) {
      return false;
    }
    const uint64_t expected =
        ComputeCrc(payload_bits.subspan(0, this->dimension_));
    uint64_t observed = 0;
    for (std::size_t i = 0; i < crc_bits_; ++i) {
      observed = (observed << 1) |
                 static_cast<uint64_t>(payload_bits[this->dimension_ + i] & 1U);
    }
    return observed == expected;
  }

  // Fast O(N log N) polar encoding (supports non-systematic, systematic, and
  // CRC-aided polar codes).
  codeword_type EncodeFast(const message_type &message) const {
    if (!message.HasField() || message.Field() != this->field_) {
      throw std::invalid_argument("PolarCode::EncodeFast: field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != this->dimension_) {
      throw std::invalid_argument(
          "PolarCode::EncodeFast: message must be 1D with length equal to "
          "Dimension()");
    }
    std::vector<uint8_t> cw_bits = EncodeRawBits(message.RawSpan());
    using RawArray = codeword_type::raw_array_type;
    using ShapeType = codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{this->length_});
    std::copy(cw_bits.begin(), cw_bits.end(), raw.data());
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  // Fast O(N log N) polar unencoding with frozen-bit and CRC validation.
  message_type UnencodeFast(const codeword_type &codeword) const {
    if (!codeword.HasField() || codeword.Field() != this->field_) {
      throw std::invalid_argument("PolarCode::UnencodeFast: field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != this->length_) {
      throw std::invalid_argument(
          "PolarCode::UnencodeFast: codeword must be 1D with length equal to "
          "Length()");
    }

    const uint8_t *cw_raw = codeword.RawData().data();
    std::vector<uint8_t> u(this->length_);
    for (std::size_t i = 0; i < this->length_; ++i) {
      u[i] = static_cast<uint8_t>(cw_raw[i] & 1U);
    }
    PolarTransform(u);

    for (std::size_t f : frozen_positions_) {
      if (u[f] != 0) {
        throw std::invalid_argument(
            "PolarCode::UnencodeFast: input is not a valid polar codeword");
      }
    }

    const std::size_t k_polar = info_positions_.size();
    std::vector<uint8_t> payload(k_polar);
    if (systematic_) {
      for (std::size_t r = 0; r < k_polar; ++r) {
        payload[r] = static_cast<uint8_t>(cw_raw[info_positions_[r]] & 1U);
      }
    } else {
      for (std::size_t r = 0; r < k_polar; ++r) {
        payload[r] = u[info_positions_[r]];
      }
    }

    if (crc_bits_ > 0 && !VerifyCrc(payload)) {
      throw std::invalid_argument(
          "PolarCode::UnencodeFast: input failed CRC validation");
    }

    using RawArray = message_type::raw_array_type;
    using ShapeType = message_type::shape_type;
    RawArray msg_raw = RawArray::from_shape(ShapeType{this->dimension_});
    std::copy(payload.begin(),
              payload.begin() + static_cast<std::ptrdiff_t>(this->dimension_),
              msg_raw.data());
    return message_type::FromRaw(this->field_, std::move(msg_raw));
  }

  // Extract the (K + crc_bits)-bit payload from decoded (u_hat, x_hat).
  std::vector<uint8_t> ExtractPayload(std::span<const uint8_t> u_hat,
                                      std::span<const uint8_t> x_hat) const {
    const std::size_t k_polar = info_positions_.size();
    std::vector<uint8_t> payload(k_polar);
    const std::span<const uint8_t> src = systematic_ ? x_hat : u_hat;
    for (std::size_t r = 0; r < k_polar; ++r) {
      payload[r] = static_cast<uint8_t>(src[info_positions_[r]] & 1U);
    }
    return payload;
  }

  std::string ToString() const override {
    std::string desc = "[" + std::to_string(this->length_) + ", " +
                       std::to_string(this->dimension_) + "] ";
    if (systematic_) {
      desc += "systematic ";
    }
    if (crc_bits_ > 0) {
      desc += "CRC-" + std::to_string(crc_bits_) + "-aided ";
    }
    desc += "polar code over GF(2)";
    return desc;
  }

 private:
  static std::size_t ValidatePowerOfTwoLength(std::size_t n,
                                              const char *fn_name) {
    if (n < 2 || !std::has_single_bit(n)) {
      throw std::invalid_argument(
          std::string(fn_name) +
          ": block length N must be a power of 2 with N >= 2");
    }
    return n;
  }

  static void ValidateBlockAndDimension(std::size_t n, std::size_t k,
                                        const char *fn_name) {
    ValidatePowerOfTwoLength(n, fn_name);
    if (k == 0 || k > n) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": dimension K must be in [1, N]");
    }
  }

  static std::size_t ValidateMessageWithCrc(std::size_t n, std::size_t k,
                                            std::size_t crc_bits,
                                            const char *fn_name) {
    ValidateBlockAndDimension(n, k, fn_name);
    if (crc_bits > 32) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": crc_bits must be between 0 and 32");
    }
    if (k + crc_bits > n) {
      throw std::invalid_argument(
          std::string(fn_name) +
          ": message dimension K + crc_bits cannot exceed block length N");
    }
    return k + crc_bits;
  }

  static std::size_t ComputeMessageDimension(
      std::size_t n, const std::vector<bool> &frozen_mask,
      std::size_t crc_bits) {
    ValidatePowerOfTwoLength(n, "PolarCode");
    if (frozen_mask.size() != n) {
      throw std::invalid_argument(
          "PolarCode: frozen_mask size must equal block length N");
    }
    std::size_t non_frozen = 0;
    for (bool is_frozen : frozen_mask) {
      if (!is_frozen) {
        ++non_frozen;
      }
    }
    if (non_frozen <= crc_bits) {
      throw std::invalid_argument(
          "PolarCode: number of non-frozen bits must exceed crc_bits");
    }
    return non_frozen - crc_bits;
  }

  std::vector<uint8_t> EncodePayloadBits(
      std::span<const uint8_t> payload) const {
    const std::size_t n = this->length_;
    const std::size_t k_polar = info_positions_.size();

    if (!systematic_) {
      std::vector<uint8_t> u(n, 0);
      for (std::size_t r = 0; r < k_polar; ++r) {
        u[info_positions_[r]] = static_cast<uint8_t>(payload[r] & 1U);
      }
      PolarTransform(u);
      return u;
    }

    if (is_upper_set_) {
      // Exact two-transform O(N log N) Arikan systematic encoder:
      // holds whenever info set A is an upper set in the Boolean lattice
      // (guaranteed for Bhattacharyya and Gaussian Approximation designs).
      std::vector<uint8_t> buf(n, 0);
      for (std::size_t r = 0; r < k_polar; ++r) {
        buf[info_positions_[r]] = static_cast<uint8_t>(payload[r] & 1U);
      }
      PolarTransform(buf);
      for (std::size_t f : frozen_positions_) {
        buf[f] = 0;
      }
      PolarTransform(buf);
      return buf;
    }

    // General exact systematic polar encoder for arbitrary custom frozen masks:
    // solve u_A * (G_N)_{A, A} = payload using the lower-triangular unit-diag
    // structure of G_N = F^{\otimes n}, then apply PolarTransform(u).
    std::vector<uint8_t> u(n, 0);
    for (std::size_t r = k_polar; r-- > 0;) {
      const std::size_t j = info_positions_[r];
      uint8_t acc = static_cast<uint8_t>(payload[r] & 1U);
      for (std::size_t q = r + 1; q < k_polar; ++q) {
        const std::size_t i = info_positions_[q];
        if ((i & j) == j) {
          acc = static_cast<uint8_t>(acc ^ u[i]);
        }
      }
      u[j] = acc;
    }
    PolarTransform(u);
    return u;
  }

  std::vector<uint8_t> EncodeRawBits(
      std::span<const uint8_t> message_bits) const {
    const std::size_t k = this->dimension_;
    if (crc_bits_ == 0) {
      return EncodePayloadBits(message_bits);
    }
    std::vector<uint8_t> payload(k + crc_bits_);
    for (std::size_t i = 0; i < k; ++i) {
      payload[i] = static_cast<uint8_t>(message_bits[i] & 1U);
    }
    const uint64_t crc = ComputeCrc(message_bits);
    for (std::size_t i = 0; i < crc_bits_; ++i) {
      payload[k + i] =
          static_cast<uint8_t>((crc >> (crc_bits_ - 1 - i)) & 1ULL);
    }
    return EncodePayloadBits(payload);
  }

  void InitializeStructure() {
    const std::size_t n = this->length_;
    const std::size_t k = this->dimension_;

    info_positions_.clear();
    frozen_positions_.clear();
    info_positions_.reserve(k + crc_bits_);
    frozen_positions_.reserve(n - (k + crc_bits_));
    for (std::size_t i = 0; i < n; ++i) {
      if (frozen_mask_[i]) {
        frozen_positions_.push_back(i);
      } else {
        info_positions_.push_back(i);
      }
    }

    // Check if info_positions_ is an upper set under bitwise subset order.
    is_upper_set_ = true;
    for (std::size_t idx : info_positions_) {
      for (std::size_t bit = 0; bit < num_stages_; ++bit) {
        const std::size_t sup = idx | (std::size_t{1} << bit);
        if (frozen_mask_[sup]) {
          is_upper_set_ = false;
          break;
        }
      }
      if (!is_upper_set_) {
        break;
      }
    }

    // Build exact K x N generator matrix G.
    using RawArray = matrix_type::raw_array_type;
    using ShapeType = matrix_type::shape_type;
    RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
    g_raw.fill(uint8_t{0});
    std::vector<uint8_t> basis_msg(k, 0);
    uint8_t *g_data = g_raw.data();
    for (std::size_t r = 0; r < k; ++r) {
      basis_msg[r] = 1U;
      std::vector<uint8_t> row_cw = EncodeRawBits(basis_msg);
      std::copy(row_cw.begin(), row_cw.end(), g_data + r * n);
      basis_msg[r] = 0U;
    }
    generator_matrix_ = matrix_type::FromRaw(this->field_, std::move(g_raw));

    // Build exact (N - K) x N parity-check matrix H satisfying G * H^T = 0.
    if (crc_bits_ == 0) {
      const std::size_t m = frozen_positions_.size();
      RawArray h_raw = RawArray::from_shape(ShapeType{m, n});
      h_raw.fill(uint8_t{0});
      uint8_t *h_data = h_raw.data();
      for (std::size_t s = 0; s < m; ++s) {
        const std::size_t f = frozen_positions_[s];
        uint8_t *h_row = h_data + s * n;
        for (std::size_t j = 0; j < n; ++j) {
          if ((j & f) == f) {
            h_row[j] = 1U;
          }
        }
      }
      parity_check_matrix_ =
          matrix_type::FromRaw(this->field_, std::move(h_raw));
    } else {
      parity_check_matrix_ =
          Base::ComputeParityCheckFromGenerator(generator_matrix_);
    }

    RegisterEncodersAndDecoders();
  }

  inline void RegisterEncodersAndDecoders();

  std::size_t num_stages_;
  std::vector<bool> frozen_mask_;
  bool systematic_;
  std::size_t crc_bits_;
  uint64_t crc_poly_;
  bool is_upper_set_ = false;
  std::vector<std::size_t> info_positions_;
  std::vector<std::size_t> frozen_positions_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;
};

// Fast O(N log N) Polar encoder (non-systematic, systematic, and CRC-aided).
class PolarEncoder : public Encoder<GaloisFieldBinary> {
 public:
  using Base = Encoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;

  explicit PolarEncoder(const PolarCode &code)
      : Base(code), polar_code_(&code) {}

  codeword_type Encode(const message_type &message) const override {
    return polar_code_->EncodeFast(message);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    return polar_code_->UnencodeFast(codeword);
  }

  std::string Name() const override { return "Polar"; }

 private:
  const PolarCode *polar_code_;
};

// O(N log N) recursive Successive Cancellation (SC) Polar decoder.
class SuccessiveCancellationDecoder : public Decoder<GaloisFieldBinary> {
 public:
  using Base = Decoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using result_type = DecodingResult<GaloisFieldBinary>;

  explicit SuccessiveCancellationDecoder(const PolarCode &code,
                                         bool use_min_sum = false,
                                         double hard_llr_magnitude = 4.0)
      : Base(code),
        polar_code_(&code),
        use_min_sum_(use_min_sum),
        hard_llr_magnitude_(hard_llr_magnitude) {
    if (!std::isfinite(hard_llr_magnitude_) || hard_llr_magnitude_ <= 0.0) {
      throw std::invalid_argument(
          "SuccessiveCancellationDecoder: hard_llr_magnitude must be positive "
          "and finite");
    }
  }

  result_type DecodeLlr(const std::vector<double> &channel_llrs) const {
    const std::size_t n = polar_code_->Length();
    polar_detail::ValidatePolarLlr(channel_llrs, n,
                                   "SuccessiveCancellationDecoder::DecodeLlr");

    std::vector<double> root_llr(n);
    std::vector<uint8_t> input_hard(n);
    for (std::size_t i = 0; i < n; ++i) {
      root_llr[i] =
          std::clamp(channel_llrs[i], -polar_detail::kPolarMaxInputLlr,
                     polar_detail::kPolarMaxInputLlr);
      input_hard[i] = (channel_llrs[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
    }

    const std::size_t stages = polar_code_->NumStages();
    std::vector<std::vector<double>> llr_levels(stages + 1);
    std::vector<std::vector<uint8_t>> beta_levels(stages + 1);
    for (std::size_t s = 0; s <= stages; ++s) {
      const std::size_t len = std::size_t{1} << s;
      llr_levels[s].resize(len);
      beta_levels[s].resize(len);
    }
    llr_levels[stages] = std::move(root_llr);

    std::vector<uint8_t> u_hat(n, 0);
    DecodeRecursive(stages, 0, llr_levels, beta_levels, u_hat);

    const std::vector<uint8_t> &x_hat = beta_levels[stages];
    const std::size_t num_errors = HammingDistance(x_hat, input_hard);
    codeword_type cw = BuildArray(this->code_->Field(), x_hat);
    std::vector<uint8_t> payload = polar_code_->ExtractPayload(u_hat, x_hat);

    if (polar_code_->CrcBits() > 0 && !polar_code_->VerifyCrc(payload)) {
      result_type res = result_type::Failure(
          "SuccessiveCancellationDecoder: decoded candidate failed CRC check");
      res.codeword = std::move(cw);
      res.num_errors = num_errors;
      return res;
    }

    const std::size_t k = polar_code_->Dimension();
    message_type msg = BuildArray(this->code_->Field(),
                                  std::span<const uint8_t>(payload.data(), k));
    return result_type::Success(std::move(cw), std::move(msg), num_errors);
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "SuccessiveCancellationDecoder::DecodeWithResult: field mismatch");
    }
    const std::size_t n = polar_code_->Length();
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "SuccessiveCancellationDecoder::DecodeWithResult: length mismatch");
    }
    const uint8_t *rx = received_word.RawData().data();
    std::vector<double> llrs(n);
    for (std::size_t i = 0; i < n; ++i) {
      llrs[i] =
          ((rx[i] & 1U) == 0) ? hard_llr_magnitude_ : -hard_llr_magnitude_;
    }
    return DecodeLlr(llrs);
  }

  std::string Name() const override { return "SC"; }
  std::size_t DecodingRadius() const override { return 0; }

 private:
  void DecodeRecursive(std::size_t stage, std::size_t offset,
                       std::vector<std::vector<double>> &llr_levels,
                       std::vector<std::vector<uint8_t>> &beta_levels,
                       std::vector<uint8_t> &u_hat) const {
    if (stage == 0) {
      const uint8_t bit =
          polar_code_->FrozenMask()[offset]
              ? uint8_t{0}
              : ((llr_levels[0][0] >= 0.0) ? uint8_t{0} : uint8_t{1});
      u_hat[offset] = bit;
      beta_levels[0][0] = bit;
      return;
    }

    const std::size_t half = std::size_t{1} << (stage - 1);
    const auto &cur_llr = llr_levels[stage];
    auto &child_llr = llr_levels[stage - 1];

    // Left branch: f(a, b).
    for (std::size_t i = 0; i < half; ++i) {
      child_llr[i] =
          use_min_sum_
              ? polar_detail::MinSumF(cur_llr[i], cur_llr[i + half])
              : polar_detail::ExactBoxPlus(cur_llr[i], cur_llr[i + half]);
    }
    DecodeRecursive(stage - 1, offset, llr_levels, beta_levels, u_hat);

    // Save left partial sum into the first half of beta_levels[stage].
    auto &cur_beta = beta_levels[stage];
    const auto &child_beta = beta_levels[stage - 1];
    for (std::size_t i = 0; i < half; ++i) {
      cur_beta[i] = child_beta[i];
    }

    // Right branch: g(a, b, beta_left).
    for (std::size_t i = 0; i < half; ++i) {
      child_llr[i] =
          polar_detail::PolarG(cur_llr[i], cur_llr[i + half], cur_beta[i]);
    }
    DecodeRecursive(stage - 1, offset + half, llr_levels, beta_levels, u_hat);

    // Combine partial sums: x_L = beta_left ^ beta_right, x_R = beta_right.
    for (std::size_t i = 0; i < half; ++i) {
      const uint8_t right_bit = child_beta[i];
      cur_beta[i] = static_cast<uint8_t>(cur_beta[i] ^ right_bit);
      cur_beta[i + half] = right_bit;
    }
  }

  static GFArray<GaloisFieldBinary> BuildArray(
      const std::shared_ptr<GaloisFieldBinary> &field,
      std::span<const uint8_t> bits) {
    using RawArray = GFArray<GaloisFieldBinary>::raw_array_type;
    using ShapeType = GFArray<GaloisFieldBinary>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{bits.size()});
    std::copy(bits.begin(), bits.end(), raw.data());
    return GFArray<GaloisFieldBinary>::FromRaw(field, std::move(raw));
  }

  static std::size_t HammingDistance(const std::vector<uint8_t> &a,
                                     const std::vector<uint8_t> &b) noexcept {
    std::size_t dist = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      if ((a[i] & 1U) != (b[i] & 1U)) {
        ++dist;
      }
    }
    return dist;
  }

  const PolarCode *polar_code_;
  bool use_min_sum_;
  double hard_llr_magnitude_;
};

// O(L * N log N) LLR-based Successive Cancellation List (SCL / CA-SCL) Polar
// decoder maintaining up to `list_size` candidate paths with LLR path metrics
// and optional CRC-aided path selection.
class SuccessiveCancellationListDecoder : public Decoder<GaloisFieldBinary> {
 public:
  using Base = Decoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using result_type = DecodingResult<GaloisFieldBinary>;

  explicit SuccessiveCancellationListDecoder(const PolarCode &code,
                                             std::size_t list_size = 8,
                                             bool use_min_sum = false,
                                             double hard_llr_magnitude = 4.0)
      : Base(code),
        polar_code_(&code),
        list_size_(list_size),
        use_min_sum_(use_min_sum),
        hard_llr_magnitude_(hard_llr_magnitude) {
    if (list_size_ == 0) {
      throw std::invalid_argument(
          "SuccessiveCancellationListDecoder: list_size must be at least 1");
    }
    if (!std::isfinite(hard_llr_magnitude_) || hard_llr_magnitude_ <= 0.0) {
      throw std::invalid_argument(
          "SuccessiveCancellationListDecoder: hard_llr_magnitude must be "
          "positive and finite");
    }
  }

  result_type DecodeLlr(const std::vector<double> &channel_llrs) const {
    const std::size_t n = polar_code_->Length();
    const std::size_t stages = polar_code_->NumStages();
    polar_detail::ValidatePolarLlr(
        channel_llrs, n, "SuccessiveCancellationListDecoder::DecodeLlr");

    std::vector<double> root_llr(n);
    std::vector<uint8_t> input_hard(n);
    for (std::size_t i = 0; i < n; ++i) {
      root_llr[i] =
          std::clamp(channel_llrs[i], -polar_detail::kPolarMaxInputLlr,
                     polar_detail::kPolarMaxInputLlr);
      input_hard[i] = (channel_llrs[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
    }

    std::vector<PathState> paths;
    paths.reserve(list_size_);
    paths.push_back(CreateInitialPath(n, stages, root_llr));

    const auto &frozen_mask = polar_code_->FrozenMask();
    std::vector<Candidate> candidates;
    candidates.reserve(2 * list_size_);

    for (std::size_t phi = 0; phi < n; ++phi) {
      for (PathState &p : paths) {
        CalcLlr(p, 0, phi, stages);
      }

      if (frozen_mask[phi]) {
        for (PathState &p : paths) {
          const double leaf_llr = p.llr[0][0];
          p.path_metric += polar_detail::PathMetricPenalty(leaf_llr, 0U);
          p.u_hat[phi] = 0U;
          UpdateBits(p, 0, phi, {0U}, stages);
        }
        continue;
      }

      // Non-frozen bit: expand every active path into bit 0 and bit 1, then
      // prune to top `list_size_` paths by ascending LLR path metric.
      candidates.clear();
      for (std::size_t idx = 0; idx < paths.size(); ++idx) {
        const double leaf_llr = paths[idx].llr[0][0];
        const double base_pm = paths[idx].path_metric;
        candidates.push_back(
            {base_pm + polar_detail::PathMetricPenalty(leaf_llr, 0U), idx,
             uint8_t{0}});
        candidates.push_back(
            {base_pm + polar_detail::PathMetricPenalty(leaf_llr, 1U), idx,
             uint8_t{1}});
      }

      const std::size_t keep = std::min(list_size_, candidates.size());
      auto cand_less = [](const Candidate &a, const Candidate &b) noexcept {
        if (a.metric != b.metric) {
          return a.metric < b.metric;
        }
        if (a.parent_idx != b.parent_idx) {
          return a.parent_idx < b.parent_idx;
        }
        return a.bit < b.bit;
      };
      if (candidates.size() > keep) {
        std::partial_sort(
            candidates.begin(),
            candidates.begin() + static_cast<std::ptrdiff_t>(keep),
            candidates.end(), cand_less);
        candidates.resize(keep);
      } else {
        std::sort(candidates.begin(), candidates.end(), cand_less);
      }

      std::vector<PathState> next_paths;
      next_paths.reserve(keep);
      for (const Candidate &cand : candidates) {
        PathState child = paths[cand.parent_idx];
        child.path_metric = cand.metric;
        child.u_hat[phi] = cand.bit;
        UpdateBits(child, 0, phi, {cand.bit}, stages);
        next_paths.push_back(std::move(child));
      }
      paths.swap(next_paths);
    }

    // Final sort by ascending path metric (accounts for any trailing frozen
    // bits after the last information position).
    std::stable_sort(paths.begin(), paths.end(),
                     [](const PathState &a, const PathState &b) noexcept {
                       return a.path_metric < b.path_metric;
                     });

    auto field = this->code_->Field();
    const std::size_t k = polar_code_->Dimension();

    if (polar_code_->CrcBits() > 0) {
      for (const PathState &p : paths) {
        std::vector<uint8_t> payload =
            polar_code_->ExtractPayload(p.u_hat, p.codeword_bits);
        if (polar_code_->VerifyCrc(payload)) {
          const std::size_t num_errors =
              HammingDistance(p.codeword_bits, input_hard);
          codeword_type cw = BuildArray(field, p.codeword_bits);
          message_type msg =
              BuildArray(field, std::span<const uint8_t>(payload.data(), k));
          return result_type::Success(std::move(cw), std::move(msg),
                                      num_errors);
        }
      }

      // No surviving path passed CRC validation.
      const PathState &best = paths[0];
      result_type res = result_type::Failure(
          "SuccessiveCancellationListDecoder: no candidate path passed CRC "
          "check");
      res.codeword = BuildArray(field, best.codeword_bits);
      res.num_errors = HammingDistance(best.codeword_bits, input_hard);
      return res;
    }

    const PathState &best = paths[0];
    std::vector<uint8_t> payload =
        polar_code_->ExtractPayload(best.u_hat, best.codeword_bits);
    const std::size_t num_errors =
        HammingDistance(best.codeword_bits, input_hard);
    codeword_type cw = BuildArray(field, best.codeword_bits);
    message_type msg =
        BuildArray(field, std::span<const uint8_t>(payload.data(), k));
    return result_type::Success(std::move(cw), std::move(msg), num_errors);
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "SuccessiveCancellationListDecoder::DecodeWithResult: field "
          "mismatch");
    }
    const std::size_t n = polar_code_->Length();
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "SuccessiveCancellationListDecoder::DecodeWithResult: length "
          "mismatch");
    }
    const uint8_t *rx = received_word.RawData().data();
    std::vector<double> llrs(n);
    for (std::size_t i = 0; i < n; ++i) {
      llrs[i] =
          ((rx[i] & 1U) == 0) ? hard_llr_magnitude_ : -hard_llr_magnitude_;
    }
    return DecodeLlr(llrs);
  }

  std::string Name() const override { return "SCL"; }
  std::size_t DecodingRadius() const override { return 0; }
  std::size_t ListSize() const noexcept { return list_size_; }

 private:
  struct PathState {
    double path_metric = 0.0;
    std::vector<uint8_t> u_hat;
    std::vector<uint8_t> codeword_bits;
    std::vector<std::vector<double>> llr;
    std::vector<std::vector<uint8_t>> bit_left;
  };

  struct Candidate {
    double metric = 0.0;
    std::size_t parent_idx = 0;
    uint8_t bit = 0;
  };

  static PathState CreateInitialPath(std::size_t n, std::size_t stages,
                                     const std::vector<double> &root_llr) {
    PathState p;
    p.path_metric = 0.0;
    p.u_hat.assign(n, 0);
    p.codeword_bits.assign(n, 0);
    p.llr.resize(stages + 1);
    p.bit_left.resize(stages);
    for (std::size_t s = 0; s <= stages; ++s) {
      p.llr[s].assign(std::size_t{1} << s, 0.0);
    }
    for (std::size_t s = 0; s < stages; ++s) {
      p.bit_left[s].assign(std::size_t{1} << s, 0);
    }
    p.llr[stages] = root_llr;
    return p;
  }

  void CalcLlr(PathState &p, std::size_t stage, std::size_t phi,
               std::size_t max_stage) const {
    if (stage == max_stage) {
      return;
    }
    const bool is_right = ((phi & 1U) != 0);
    if (!is_right) {
      CalcLlr(p, stage + 1, phi >> 1, max_stage);
    }
    const std::size_t len = std::size_t{1} << stage;
    const auto &parent_llr = p.llr[stage + 1];
    auto &cur_llr = p.llr[stage];
    if (!is_right) {
      for (std::size_t i = 0; i < len; ++i) {
        cur_llr[i] = use_min_sum_ ? polar_detail::MinSumF(parent_llr[i],
                                                          parent_llr[i + len])
                                  : polar_detail::ExactBoxPlus(
                                        parent_llr[i], parent_llr[i + len]);
      }
    } else {
      const auto &left_bits = p.bit_left[stage];
      for (std::size_t i = 0; i < len; ++i) {
        cur_llr[i] = polar_detail::PolarG(parent_llr[i], parent_llr[i + len],
                                          left_bits[i]);
      }
    }
  }

  void UpdateBits(PathState &p, std::size_t stage, std::size_t phi,
                  const std::vector<uint8_t> &cur_bits,
                  std::size_t max_stage) const {
    const bool is_right = ((phi & 1U) != 0);
    if (!is_right) {
      p.bit_left[stage] = cur_bits;
      return;
    }
    const std::size_t len = std::size_t{1} << stage;
    std::vector<uint8_t> combined(2 * len);
    const auto &left_bits = p.bit_left[stage];
    for (std::size_t i = 0; i < len; ++i) {
      combined[i] = static_cast<uint8_t>(left_bits[i] ^ cur_bits[i]);
      combined[i + len] = cur_bits[i];
    }
    if (stage + 1 < max_stage) {
      UpdateBits(p, stage + 1, phi >> 1, combined, max_stage);
    } else {
      p.codeword_bits = std::move(combined);
    }
  }

  static GFArray<GaloisFieldBinary> BuildArray(
      const std::shared_ptr<GaloisFieldBinary> &field,
      std::span<const uint8_t> bits) {
    using RawArray = GFArray<GaloisFieldBinary>::raw_array_type;
    using ShapeType = GFArray<GaloisFieldBinary>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{bits.size()});
    std::copy(bits.begin(), bits.end(), raw.data());
    return GFArray<GaloisFieldBinary>::FromRaw(field, std::move(raw));
  }

  static std::size_t HammingDistance(const std::vector<uint8_t> &a,
                                     const std::vector<uint8_t> &b) noexcept {
    std::size_t dist = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      if ((a[i] & 1U) != (b[i] & 1U)) {
        ++dist;
      }
    }
    return dist;
  }

  const PolarCode *polar_code_;
  std::size_t list_size_;
  bool use_min_sum_;
  double hard_llr_magnitude_;
};

inline void PolarCode::RegisterEncodersAndDecoders() {
  this->RegisterEncoder(
      "Polar", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &polar = dynamic_cast<const PolarCode &>(code);
        return std::make_shared<PolarEncoder>(polar);
      });
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<GaloisFieldBinary> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<GaloisFieldBinary>>(lin);
      });
  this->RegisterDecoder("SC", [](const AbstractCode<GaloisFieldBinary> &code) {
    const auto &polar = dynamic_cast<const PolarCode &>(code);
    return std::make_shared<SuccessiveCancellationDecoder>(polar);
  });
  this->RegisterDecoder("SCL", [](const AbstractCode<GaloisFieldBinary> &code) {
    const auto &polar = dynamic_cast<const PolarCode &>(code);
    return std::make_shared<SuccessiveCancellationListDecoder>(polar, 8);
  });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_POLAR_HPP_
