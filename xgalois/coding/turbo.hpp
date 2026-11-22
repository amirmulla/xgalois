#ifndef XGALOIS_CODING_TURBO_HPP_
#define XGALOIS_CODING_TURBO_HPP_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/convolutional.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

namespace turbo_detail {

// Overflow-safe standard C++20 modular multiplication (a * b) % mod.
inline std::size_t MulMod(std::size_t a, std::size_t b,
                          std::size_t mod) noexcept {
  uint64_t x = static_cast<uint64_t>(a % mod);
  uint64_t y = static_cast<uint64_t>(b % mod);
  const uint64_t m = static_cast<uint64_t>(mod);
  uint64_t res = 0;
  while (y > 0) {
    if ((y & 1ULL) != 0ULL) {
      res = (res >= m - x) ? (res - (m - x)) : (res + x);
    }
    x = (x >= m - x) ? (x - (m - x)) : (x + x);
    y >>= 1;
  }
  return static_cast<std::size_t>(res);
}

inline std::size_t AddMod(std::size_t a, std::size_t b,
                          std::size_t mod) noexcept {
  const std::size_t x = a % mod;
  const std::size_t y = b % mod;
  return (x >= mod - y) ? (x - (mod - y)) : (x + y);
}

}  // namespace turbo_detail

// Permutation interleaver supporting 3GPP LTE Quadratic Permutation Polynomial
// (QPP), S-random, block matrix, and explicit permutations.
//
// Convention (matching 3GPP TS 36.212):
//   interleaved[i] = input[perm[i]]
//   deinterleaved[perm[i]] = interleaved[i]
class Interleaver {
 public:
  // Quadratic Permutation Polynomial (QPP): pi(i) = (f1 * i + f2 * i^2) mod K.
  static Interleaver QppInterleaver(std::size_t K, std::size_t f1,
                                    std::size_t f2) {
    if (K == 0) {
      throw std::invalid_argument(
          "Interleaver::QppInterleaver: block length K must be positive");
    }
    std::vector<std::size_t> perm(K);
    for (std::size_t i = 0; i < K; ++i) {
      const std::size_t i_sq = turbo_detail::MulMod(i, i, K);
      const std::size_t term1 = turbo_detail::MulMod(f1, i, K);
      const std::size_t term2 = turbo_detail::MulMod(f2, i_sq, K);
      perm[i] = turbo_detail::AddMod(term1, term2, K);
    }
    return ExplicitInterleaver(std::move(perm));
  }

  // Standard 3GPP TS 36.212 Table 5.1.3-3 LTE QPP interleaver presets.
  static Interleaver CreateLteQpp(std::size_t K) {
    struct LteQppParams {
      std::size_t k;
      std::size_t f1;
      std::size_t f2;
    };
    static constexpr LteQppParams kTable[] = {
        {40, 3, 10},      {48, 7, 12},    {56, 19, 42},     {64, 7, 16},
        {72, 7, 18},      {80, 11, 20},   {96, 11, 24},     {104, 7, 26},
        {112, 41, 84},    {120, 103, 90}, {128, 15, 32},    {144, 17, 108},
        {160, 21, 120},   {192, 23, 48},  {224, 27, 56},    {256, 15, 32},
        {320, 21, 120},   {384, 23, 48},  {448, 29, 168},   {512, 31, 64},
        {640, 39, 80},    {768, 217, 48}, {896, 157, 56},   {1024, 31, 64},
        {1536, 43, 384},  {2048, 31, 64}, {3072, 157, 192}, {4096, 31, 64},
        {6144, 263, 480},
    };
    for (const auto &entry : kTable) {
      if (entry.k == K) {
        return QppInterleaver(entry.k, entry.f1, entry.f2);
      }
    }
    throw std::invalid_argument(
        "Interleaver::CreateLteQpp: unsupported LTE block length K=" +
        std::to_string(K));
  }

  // S-random interleaver: guarantees |pi(j) - pi(i)| >= s_spread whenever
  // 0 < j - i < s_spread.
  static Interleaver SRandomInterleaver(std::size_t K, std::size_t s_spread,
                                        uint64_t seed = 42) {
    if (K == 0) {
      throw std::invalid_argument(
          "Interleaver::SRandomInterleaver: block length K must be positive");
    }
    if (s_spread > K) {
      throw std::invalid_argument(
          "Interleaver::SRandomInterleaver: s_spread cannot exceed K");
    }
    std::mt19937_64 rng(seed);
    if (s_spread <= 1) {
      std::vector<std::size_t> perm(K);
      std::iota(perm.begin(), perm.end(), std::size_t{0});
      std::shuffle(perm.begin(), perm.end(), rng);
      return ExplicitInterleaver(std::move(perm));
    }

    constexpr std::size_t kMaxRestarts = 200;
    std::vector<std::size_t> perm(K, 0);
    std::vector<std::size_t> pool(K, 0);

    for (std::size_t attempt = 0; attempt < kMaxRestarts; ++attempt) {
      std::iota(pool.begin(), pool.end(), std::size_t{0});
      std::size_t remaining = K;
      bool failed = false;

      for (std::size_t i = 0; i < K; ++i) {
        const std::size_t window_start =
            (i >= s_spread) ? (i - s_spread + 1) : 0;
        const std::size_t max_tries =
            std::min<std::size_t>(remaining * 4 + 16, 256);
        bool placed = false;

        for (std::size_t tries = 0; tries < max_tries; ++tries) {
          std::uniform_int_distribution<std::size_t> dist(0, remaining - 1);
          const std::size_t pick_idx = dist(rng);
          const std::size_t candidate = pool[pick_idx];

          bool ok = true;
          for (std::size_t j = window_start; j < i; ++j) {
            const std::size_t diff = (candidate > perm[j])
                                         ? (candidate - perm[j])
                                         : (perm[j] - candidate);
            if (diff < s_spread) {
              ok = false;
              break;
            }
          }
          if (ok) {
            perm[i] = candidate;
            pool[pick_idx] = pool[remaining - 1];
            --remaining;
            placed = true;
            break;
          }
        }

        if (!placed) {
          // Deterministic scan over remaining pool before giving up on attempt.
          for (std::size_t pick_idx = 0; pick_idx < remaining; ++pick_idx) {
            const std::size_t candidate = pool[pick_idx];
            bool ok = true;
            for (std::size_t j = window_start; j < i; ++j) {
              const std::size_t diff = (candidate > perm[j])
                                           ? (candidate - perm[j])
                                           : (perm[j] - candidate);
              if (diff < s_spread) {
                ok = false;
                break;
              }
            }
            if (ok) {
              perm[i] = candidate;
              pool[pick_idx] = pool[remaining - 1];
              --remaining;
              placed = true;
              break;
            }
          }
        }

        if (!placed) {
          failed = true;
          break;
        }
      }

      if (!failed) {
        return ExplicitInterleaver(std::move(perm));
      }
    }

    throw std::invalid_argument(
        "Interleaver::SRandomInterleaver: failed to construct S-random "
        "interleaver for K=" +
        std::to_string(K) + ", s_spread=" + std::to_string(s_spread));
  }

  // Row-write, column-read block matrix interleaver of size rows * cols.
  static Interleaver MatrixInterleaver(std::size_t rows, std::size_t cols) {
    if (rows == 0 || cols == 0) {
      throw std::invalid_argument(
          "Interleaver::MatrixInterleaver: rows and cols must be positive");
    }
    const std::size_t k = rows * cols;
    std::vector<std::size_t> perm(k);
    for (std::size_t c = 0; c < cols; ++c) {
      for (std::size_t r = 0; r < rows; ++r) {
        perm[c * rows + r] = r * cols + c;
      }
    }
    return ExplicitInterleaver(std::move(perm));
  }

  // Explicit permutation interleaver with strict bijection validation.
  static Interleaver ExplicitInterleaver(std::vector<std::size_t> perm) {
    const std::size_t k = perm.size();
    if (k == 0) {
      throw std::invalid_argument(
          "Interleaver::ExplicitInterleaver: permutation cannot be empty");
    }
    std::vector<bool> seen(k, false);
    std::vector<std::size_t> inv(k, 0);
    for (std::size_t i = 0; i < k; ++i) {
      const std::size_t dst = perm[i];
      if (dst >= k || seen[dst]) {
        throw std::invalid_argument(
            "Interleaver: permutation is not a valid bijection on [0, K)");
      }
      seen[dst] = true;
      inv[dst] = i;
    }
    return Interleaver(std::move(perm), std::move(inv));
  }

  std::size_t Size() const noexcept { return perm_.size(); }
  const std::vector<std::size_t> &Permutation() const noexcept { return perm_; }
  const std::vector<std::size_t> &InversePermutation() const noexcept {
    return inv_perm_;
  }

  template <typename T>
  std::vector<T> Interleave(std::span<const T> input) const {
    if (input.size() != perm_.size()) {
      throw std::invalid_argument(
          "Interleaver::Interleave: input size mismatch");
    }
    std::vector<T> out(perm_.size());
    for (std::size_t i = 0; i < perm_.size(); ++i) {
      out[i] = input[perm_[i]];
    }
    return out;
  }

  template <typename T>
  std::vector<T> Interleave(const std::vector<T> &input) const {
    return Interleave(std::span<const T>(input));
  }

  template <typename T>
  std::vector<T> Deinterleave(std::span<const T> input) const {
    if (input.size() != perm_.size()) {
      throw std::invalid_argument(
          "Interleaver::Deinterleave: input size mismatch");
    }
    std::vector<T> out(perm_.size());
    for (std::size_t i = 0; i < perm_.size(); ++i) {
      out[perm_[i]] = input[i];
    }
    return out;
  }

  template <typename T>
  std::vector<T> Deinterleave(const std::vector<T> &input) const {
    return Deinterleave(std::span<const T>(input));
  }

  GFArray<GaloisFieldBinary> Interleave(
      const GFArray<GaloisFieldBinary> &input) const {
    if (!input.HasField()) {
      throw std::invalid_argument(
          "Interleaver::Interleave: GFArray has null field");
    }
    if (input.Dimension() != 1 || input.Size() != perm_.size()) {
      throw std::invalid_argument(
          "Interleaver::Interleave: GFArray size mismatch");
    }
    std::vector<uint8_t> out = Interleave(input.RawSpan());
    return conv_detail::BuildBinaryArray(input.Field(), out);
  }

  GFArray<GaloisFieldBinary> Deinterleave(
      const GFArray<GaloisFieldBinary> &input) const {
    if (!input.HasField()) {
      throw std::invalid_argument(
          "Interleaver::Deinterleave: GFArray has null field");
    }
    if (input.Dimension() != 1 || input.Size() != perm_.size()) {
      throw std::invalid_argument(
          "Interleaver::Deinterleave: GFArray size mismatch");
    }
    std::vector<uint8_t> out = Deinterleave(input.RawSpan());
    return conv_detail::BuildBinaryArray(input.Field(), out);
  }

 private:
  Interleaver(std::vector<std::size_t> perm, std::vector<std::size_t> inv_perm)
      : perm_(std::move(perm)), inv_perm_(std::move(inv_perm)) {}

  std::vector<std::size_t> perm_;
  std::vector<std::size_t> inv_perm_;
};

struct TurboDecodingResult {
  GFArray<GaloisFieldBinary> message;
  GFArray<GaloisFieldBinary> codeword;
  std::vector<double> app_llrs;
  std::size_t iterations_used = 0;
  bool converged = false;
  bool success = false;
};

class TurboCode;
class TurboEncoder;
class TurboHardDecoder;

// Parallel Concatenated Convolutional Code (PCCC / Turbo Code) over GF(2)
// combining two Recursive Systematic Convolutional (RSC) constituents, an
// Interleaver, and an optional periodic puncturing mask on the data section.
//
// Unpunctured codeword layout:
//   - Data section (K steps): [u_t, parity1_t..., parity2_t...] for t = 0..K-1
//   - Tail section (unpunctured when constituents are zero-tailed):
//       [RSC1 tail bits (m1 * n0_1)] followed by [RSC2 tail bits (m2 * n0_2)]
class TurboCode : public AbstractLinearCode<GaloisFieldBinary> {
 public:
  using Base = AbstractLinearCode<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using matrix_type = GFMatrix<GaloisFieldBinary>;
  using syndrome_type = GFArray<GaloisFieldBinary>;

  // Standard alternating parity puncturing mask for rate-1/3 -> rate-1/2 PCCC
  // over (sys, p1, p2) per step: keeps (sys, p1) on even steps and (sys, p2) on
  // odd steps.
  static std::vector<bool> RateHalfPunctureMask() {
    return {true, true, false, true, false, true};
  }

  // Construct a TurboCode from two RSC constituents, an interleaver, and an
  // optional periodic puncturing mask applied to the data section.
  TurboCode(std::shared_ptr<GaloisFieldBinary> field, ConvolutionalCode rsc1,
            ConvolutionalCode rsc2, Interleaver interleaver,
            std::vector<bool> puncture_mask = {})
      : Base(std::move(field),
             ComputePuncturedLength(rsc1, rsc2, interleaver, puncture_mask),
             interleaver.Size(), Metric::HAMMING, "Turbo", "TurboHard"),
        rsc1_(std::move(rsc1)),
        rsc2_(std::move(rsc2)),
        interleaver_(std::move(interleaver)),
        puncture_mask_(std::move(puncture_mask)) {
    InitializeLayout();
    RegisterEncodersAndDecoders();
  }

  // Convenience constructor using identical 3GPP LTE RSC constituents
  // (feedback=013, feedforward={015}) with a custom interleaver & puncturing.
  TurboCode(std::shared_ptr<GaloisFieldBinary> field,
            const Interleaver &interleaver,
            std::vector<bool> puncture_mask = {})
      : TurboCode(
            field,
            ConvolutionalCode::LteTurboConstituent(field, interleaver.Size()),
            ConvolutionalCode::LteTurboConstituent(field, interleaver.Size()),
            interleaver, std::move(puncture_mask)) {}

  // Factory: 3GPP LTE Turbo code with standard QPP interleaver for block length
  // K (rate ~1/3 unpunctured, or rate ~1/2 with alternating parity puncturing).
  static TurboCode CreateLte(std::shared_ptr<GaloisFieldBinary> field,
                             std::size_t K, bool rate_half_punctured = false) {
    Interleaver interleaver = Interleaver::CreateLteQpp(K);
    std::vector<bool> mask =
        rate_half_punctured ? RateHalfPunctureMask() : std::vector<bool>{};
    return TurboCode(std::move(field), std::move(interleaver), std::move(mask));
  }

  static TurboCode CreateLte(std::size_t K, bool rate_half_punctured = false) {
    return CreateLte(std::make_shared<GaloisFieldBinary>(), K,
                     rate_half_punctured);
  }

  // Rule of Five (mutex is not copyable/movable).
  TurboCode(const TurboCode &other)
      : Base(other),
        rsc1_(other.rsc1_),
        rsc2_(other.rsc2_),
        interleaver_(other.interleaver_),
        puncture_mask_(other.puncture_mask_),
        unpunctured_data_len_(other.unpunctured_data_len_),
        unpunctured_total_len_(other.unpunctured_total_len_),
        keep_mask_(other.keep_mask_) {
    std::lock_guard<std::mutex> lock(other.matrix_mutex_);
    generator_matrix_ = other.generator_matrix_;
    parity_check_matrix_ = other.parity_check_matrix_;
  }

  TurboCode &operator=(const TurboCode &other) {
    if (this != &other) {
      Base::operator=(other);
      rsc1_ = other.rsc1_;
      rsc2_ = other.rsc2_;
      interleaver_ = other.interleaver_;
      puncture_mask_ = other.puncture_mask_;
      unpunctured_data_len_ = other.unpunctured_data_len_;
      unpunctured_total_len_ = other.unpunctured_total_len_;
      keep_mask_ = other.keep_mask_;
      std::scoped_lock lock(matrix_mutex_, other.matrix_mutex_);
      generator_matrix_ = other.generator_matrix_;
      parity_check_matrix_ = other.parity_check_matrix_;
    }
    return *this;
  }

  TurboCode(TurboCode &&other) noexcept
      : Base(std::move(other)),
        rsc1_(std::move(other.rsc1_)),
        rsc2_(std::move(other.rsc2_)),
        interleaver_(std::move(other.interleaver_)),
        puncture_mask_(std::move(other.puncture_mask_)),
        unpunctured_data_len_(other.unpunctured_data_len_),
        unpunctured_total_len_(other.unpunctured_total_len_),
        keep_mask_(std::move(other.keep_mask_)),
        generator_matrix_(std::move(other.generator_matrix_)),
        parity_check_matrix_(std::move(other.parity_check_matrix_)) {}

  TurboCode &operator=(TurboCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      rsc1_ = std::move(other.rsc1_);
      rsc2_ = std::move(other.rsc2_);
      interleaver_ = std::move(other.interleaver_);
      puncture_mask_ = std::move(other.puncture_mask_);
      unpunctured_data_len_ = other.unpunctured_data_len_;
      unpunctured_total_len_ = other.unpunctured_total_len_;
      keep_mask_ = std::move(other.keep_mask_);
      generator_matrix_ = std::move(other.generator_matrix_);
      parity_check_matrix_ = std::move(other.parity_check_matrix_);
    }
    return *this;
  }

  ~TurboCode() override = default;

  const ConvolutionalCode &Constituent1() const noexcept { return rsc1_; }
  const ConvolutionalCode &Constituent2() const noexcept { return rsc2_; }
  const Interleaver &GetInterleaver() const noexcept { return interleaver_; }
  const std::vector<bool> &PunctureMask() const noexcept {
    return puncture_mask_;
  }
  std::size_t UnpuncturedLength() const noexcept {
    return unpunctured_total_len_;
  }

  // Fast O(K) systematic PCCC encoding on raw binary bits (0/1).
  std::vector<uint8_t> EncodeBits(std::span<const uint8_t> message_bits) const {
    const std::size_t k = this->dimension_;
    if (message_bits.size() != k) {
      throw std::invalid_argument(
          "TurboCode::EncodeBits: message length mismatch");
    }

    const std::vector<uint8_t> cw1 = rsc1_.EncodeBits(message_bits);
    const std::vector<uint8_t> interleaved_msg =
        interleaver_.Interleave(message_bits);
    const std::vector<uint8_t> cw2 =
        rsc2_.EncodeBits(std::span<const uint8_t>(interleaved_msg));

    const std::size_t n0_1 = rsc1_.NumOutputBits();
    const std::size_t n0_2 = rsc2_.NumOutputBits();
    const std::size_t p1 = n0_1 - 1;
    const std::size_t p2 = n0_2 - 1;

    std::vector<uint8_t> unpunctured;
    unpunctured.reserve(unpunctured_total_len_);

    // Data section: [sys, parity1..., parity2...] per step t.
    for (std::size_t t = 0; t < k; ++t) {
      unpunctured.push_back(static_cast<uint8_t>(message_bits[t] & 1U));
      for (std::size_t j = 0; j < p1; ++j) {
        unpunctured.push_back(cw1[t * n0_1 + 1 + j]);
      }
      for (std::size_t j = 0; j < p2; ++j) {
        unpunctured.push_back(cw2[t * n0_2 + 1 + j]);
      }
    }

    // Unpunctured termination tails.
    const std::size_t tail1_start = k * n0_1;
    for (std::size_t i = tail1_start; i < cw1.size(); ++i) {
      unpunctured.push_back(cw1[i]);
    }
    const std::size_t tail2_start = k * n0_2;
    for (std::size_t i = tail2_start; i < cw2.size(); ++i) {
      unpunctured.push_back(cw2[i]);
    }

    if (puncture_mask_.empty()) {
      return unpunctured;
    }

    std::vector<uint8_t> punctured;
    punctured.reserve(this->length_);
    for (std::size_t i = 0; i < unpunctured_total_len_; ++i) {
      if (keep_mask_[i]) {
        punctured.push_back(unpunctured[i]);
      }
    }
    return punctured;
  }

  // Fast O(K) systematic encoding returning a GFArray<GaloisFieldBinary>.
  codeword_type EncodeFast(const message_type &message) const {
    if (!message.HasField() || message.Field() != this->field_) {
      throw std::invalid_argument("TurboCode::EncodeFast: field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != this->dimension_) {
      throw std::invalid_argument(
          "TurboCode::EncodeFast: message length mismatch");
    }
    std::vector<uint8_t> bits = EncodeBits(message.RawSpan());
    return conv_detail::BuildBinaryArray(this->field_, bits);
  }

  // Fast unencoding for unpunctured (or systematic-preserving punctured) Turbo
  // codewords with full codeword verification.
  message_type UnencodeFast(const codeword_type &codeword) const {
    if (!codeword.HasField() || codeword.Field() != this->field_) {
      throw std::invalid_argument("TurboCode::UnencodeFast: field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != this->length_) {
      throw std::invalid_argument(
          "TurboCode::UnencodeFast: codeword length mismatch");
    }
    const auto span = codeword.RawSpan();
    const std::size_t k = this->dimension_;
    const std::size_t step_width =
        1 + (rsc1_.NumOutputBits() - 1) + (rsc2_.NumOutputBits() - 1);

    // Check whether all systematic positions are retained (not punctured).
    bool all_sys_kept = true;
    for (std::size_t t = 0; t < k; ++t) {
      if (!keep_mask_[t * step_width]) {
        all_sys_kept = false;
        break;
      }
    }

    std::vector<uint8_t> msg_bits(k, 0U);
    if (all_sys_kept) {
      std::size_t punc_idx = 0;
      for (std::size_t i = 0; i < unpunctured_data_len_; ++i) {
        if (keep_mask_[i]) {
          if (i % step_width == 0) {
            msg_bits[i / step_width] =
                static_cast<uint8_t>(span[punc_idx] & 1U);
          }
          ++punc_idx;
        }
      }
      const std::vector<uint8_t> reencoded = EncodeBits(msg_bits);
      for (std::size_t i = 0; i < this->length_; ++i) {
        if ((reencoded[i] & 1U) != (span[i] & 1U)) {
          throw std::invalid_argument(
              "TurboCode::UnencodeFast: input is not a valid Turbo codeword");
        }
      }
      return conv_detail::BuildBinaryArray(this->field_, msg_bits);
    }

    // Fallback for puncturing masks that puncture systematic bits: decode with
    // strong hard-decision LLRs and verify exact codeword match.
    std::vector<double> hard_llrs(this->length_);
    for (std::size_t i = 0; i < this->length_; ++i) {
      hard_llrs[i] = ((span[i] & 1U) == 0U) ? 8.0 : -8.0;
    }
    auto res = DecodeSoft(hard_llrs, 8, BcjrMode::kLogMap, 1.0);
    const auto cw_span = res.codeword.RawSpan();
    for (std::size_t i = 0; i < this->length_; ++i) {
      if ((cw_span[i] & 1U) != (span[i] & 1U)) {
        throw std::invalid_argument(
            "TurboCode::UnencodeFast: input is not a valid Turbo codeword");
      }
    }
    return std::move(res.message);
  }

  // Iterative Soft-Input Soft-Output (SISO) BCJR Turbo decoder with
  // depuncturing, extrinsic scaling, and consecutive half-iteration early
  // stopping.
  TurboDecodingResult DecodeSoft(const std::vector<double> &channel_llrs,
                                 std::size_t max_iterations = 8,
                                 BcjrMode mode = BcjrMode::kLogMap,
                                 double extrinsic_scale = 1.0) const {
    conv_detail::ValidateLlrVector(channel_llrs, this->length_,
                                   "TurboCode::DecodeSoft");
    if (max_iterations == 0) {
      throw std::invalid_argument(
          "TurboCode::DecodeSoft: max_iterations must be at least 1");
    }
    if (!std::isfinite(extrinsic_scale) || extrinsic_scale <= 0.0) {
      throw std::invalid_argument(
          "TurboCode::DecodeSoft: extrinsic_scale must be positive and finite");
    }

    const std::size_t k = this->dimension_;
    const std::size_t n0_1 = rsc1_.NumOutputBits();
    const std::size_t n0_2 = rsc2_.NumOutputBits();
    const std::size_t p1 = n0_1 - 1;
    const std::size_t p2 = n0_2 - 1;

    // Step 1: Depuncture channel LLRs (inserting 0.0 at punctured positions).
    std::vector<double> full_llrs(unpunctured_total_len_, 0.0);
    std::size_t src_idx = 0;
    for (std::size_t i = 0; i < unpunctured_total_len_; ++i) {
      if (keep_mask_[i]) {
        full_llrs[i] = conv_detail::ClampLlr(channel_llrs[src_idx++]);
      }
    }

    // Step 2: Demultiplex channel LLRs for RSC-1 and RSC-2.
    std::vector<double> sys_llr(k, 0.0);
    std::vector<double> llr_c1(rsc1_.Length(), 0.0);
    std::vector<double> llr_c2(rsc2_.Length(), 0.0);

    std::size_t pos = 0;
    for (std::size_t t = 0; t < k; ++t) {
      const double s_llr = full_llrs[pos++];
      sys_llr[t] = s_llr;
      llr_c1[t * n0_1] = s_llr;
      for (std::size_t j = 0; j < p1; ++j) {
        llr_c1[t * n0_1 + 1 + j] = full_llrs[pos++];
      }
      for (std::size_t j = 0; j < p2; ++j) {
        llr_c2[t * n0_2 + 1 + j] = full_llrs[pos++];
      }
    }

    // Interleaved systematic channel LLRs for RSC-2.
    const std::vector<double> sys_llr_int = interleaver_.Interleave(sys_llr);
    for (std::size_t t = 0; t < k; ++t) {
      llr_c2[t * n0_2] = sys_llr_int[t];
    }

    // Termination tail LLRs.
    const std::size_t tail1_len = rsc1_.Length() - k * n0_1;
    for (std::size_t i = 0; i < tail1_len; ++i) {
      llr_c1[k * n0_1 + i] = full_llrs[pos++];
    }
    const std::size_t tail2_len = rsc2_.Length() - k * n0_2;
    for (std::size_t i = 0; i < tail2_len; ++i) {
      llr_c2[k * n0_2 + i] = full_llrs[pos++];
    }

    // Step 3: Iterative SISO BCJR message passing.
    BcjrDecoder bcjr1(rsc1_);
    BcjrDecoder bcjr2(rsc2_);

    std::vector<double> a_priori_1(k, 0.0);
    std::vector<double> ext1_scaled(k, 0.0);
    std::vector<double> ext2_scaled(k, 0.0);
    std::vector<double> final_app(k, 0.0);
    std::vector<uint8_t> hard1(k, 0U);
    std::vector<uint8_t> hard2_deint(k, 0U);
    std::vector<uint8_t> prev_hard(k, 0U);

    bool converged = false;
    std::size_t iterations_used = 0;

    for (std::size_t iter = 1; iter <= max_iterations; ++iter) {
      iterations_used = iter;

      // Half-iteration 1: SISO BCJR on RSC-1 (natural order).
      BcjrOutput out1 = bcjr1.DecodeSoft(llr_c1, a_priori_1, mode);
      for (std::size_t i = 0; i < k; ++i) {
        hard1[i] = (out1.app_llrs[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
        ext1_scaled[i] = extrinsic_scale * out1.extrinsic_llrs[i];
      }
      if (iter > 1 && hard1 == prev_hard) {
        final_app = std::move(out1.app_llrs);
        converged = true;
        break;
      }

      // Interleave extrinsic LLRs from RSC-1 -> a priori LLRs for RSC-2.
      std::vector<double> a_priori_2 = interleaver_.Interleave(ext1_scaled);

      // Half-iteration 2: SISO BCJR on RSC-2 (interleaved order).
      BcjrOutput out2 = bcjr2.DecodeSoft(llr_c2, a_priori_2, mode);
      final_app = interleaver_.Deinterleave(out2.app_llrs);
      for (std::size_t i = 0; i < k; ++i) {
        hard2_deint[i] = (final_app[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
        ext2_scaled[i] = extrinsic_scale * out2.extrinsic_llrs[i];
      }

      if (hard2_deint == hard1) {
        converged = true;
        break;
      }

      prev_hard = hard2_deint;
      a_priori_1 = interleaver_.Deinterleave(ext2_scaled);
    }

    std::vector<uint8_t> final_hard(k, 0U);
    for (std::size_t i = 0; i < k; ++i) {
      final_hard[i] = (final_app[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
    }

    TurboDecodingResult result;
    result.message = conv_detail::BuildBinaryArray(this->field_, final_hard);
    result.codeword = EncodeFast(result.message);
    result.app_llrs = std::move(final_app);
    result.iterations_used = iterations_used;
    result.converged = converged;
    result.success = true;
    return result;
  }

  // Lazy O(K * N) generator matrix construction.
  const matrix_type &GeneratorMatrix() const override {
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!generator_matrix_.has_value()) {
      const std::size_t k = this->dimension_;
      const std::size_t n = this->length_;
      using RawArray = matrix_type::raw_array_type;
      using ShapeType = matrix_type::shape_type;
      RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
      uint8_t *g_data = g_raw.data();
      std::vector<uint8_t> basis(k, 0U);
      for (std::size_t r = 0; r < k; ++r) {
        basis[r] = 1U;
        std::vector<uint8_t> cw = EncodeBits(basis);
        std::copy(cw.begin(), cw.end(), g_data + r * n);
        basis[r] = 0U;
      }
      generator_matrix_ = matrix_type::FromRaw(this->field_, std::move(g_raw));
    }
    return *generator_matrix_;
  }

  // Lazy parity-check matrix construction via GF(2) nullspace of G.
  const matrix_type &ParityCheckMatrix() const override {
    const matrix_type &g = GeneratorMatrix();
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!parity_check_matrix_.has_value()) {
      parity_check_matrix_ = Base::ComputeParityCheckFromGenerator(g);
    }
    return *parity_check_matrix_;
  }

  // O(K) random codeword sampling without materializing GeneratorMatrix().
  codeword_type RandomCodeword() const override {
    auto msg = message_type::Random(this->field_, {this->dimension_});
    return EncodeFast(msg);
  }

  codeword_type RandomCodeword(std::mt19937_64 &rng) const override {
    auto msg = message_type::Random(this->field_, {this->dimension_}, rng);
    return EncodeFast(msg);
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) +
           "] PCCC Turbo code over GF(2) (unpunctured_len=" +
           std::to_string(unpunctured_total_len_) + ")";
  }

 private:
  static std::size_t ComputePuncturedLength(
      const ConvolutionalCode &rsc1, const ConvolutionalCode &rsc2,
      const Interleaver &interleaver, const std::vector<bool> &puncture_mask) {
    if (!rsc1.IsRecursiveSystematic() || !rsc2.IsRecursiveSystematic()) {
      throw std::invalid_argument(
          "TurboCode: constituent encoders must be recursive systematic (RSC)");
    }
    const std::size_t k = interleaver.Size();
    if (rsc1.InputBlockLength() != k || rsc2.InputBlockLength() != k) {
      throw std::invalid_argument(
          "TurboCode: constituent input block lengths must match interleaver "
          "size");
    }
    const std::size_t p1 = rsc1.NumOutputBits() - 1;
    const std::size_t p2 = rsc2.NumOutputBits() - 1;
    const std::size_t data_len = k * (1 + p1 + p2);
    const std::size_t tail1_len = rsc1.Length() - k * rsc1.NumOutputBits();
    const std::size_t tail2_len = rsc2.Length() - k * rsc2.NumOutputBits();

    if (puncture_mask.empty()) {
      return data_len + tail1_len + tail2_len;
    }

    std::size_t kept_data = 0;
    const std::size_t period = puncture_mask.size();
    for (std::size_t i = 0; i < data_len; ++i) {
      if (puncture_mask[i % period]) {
        ++kept_data;
      }
    }
    if (kept_data == 0) {
      throw std::invalid_argument(
          "TurboCode: puncture_mask must retain at least one bit");
    }
    const std::size_t total = kept_data + tail1_len + tail2_len;
    if (total < k) {
      throw std::invalid_argument(
          "TurboCode: punctured code length cannot be less than dimension K");
    }
    return total;
  }

  void InitializeLayout() {
    const std::size_t k = this->dimension_;
    const std::size_t p1 = rsc1_.NumOutputBits() - 1;
    const std::size_t p2 = rsc2_.NumOutputBits() - 1;
    unpunctured_data_len_ = k * (1 + p1 + p2);
    const std::size_t tail1_len = rsc1_.Length() - k * rsc1_.NumOutputBits();
    const std::size_t tail2_len = rsc2_.Length() - k * rsc2_.NumOutputBits();
    unpunctured_total_len_ = unpunctured_data_len_ + tail1_len + tail2_len;

    keep_mask_.assign(unpunctured_total_len_, true);
    if (!puncture_mask_.empty()) {
      const std::size_t period = puncture_mask_.size();
      for (std::size_t i = 0; i < unpunctured_data_len_; ++i) {
        keep_mask_[i] = puncture_mask_[i % period];
      }
    }
  }

  inline void RegisterEncodersAndDecoders();

  ConvolutionalCode rsc1_;
  ConvolutionalCode rsc2_;
  Interleaver interleaver_;
  std::vector<bool> puncture_mask_;

  std::size_t unpunctured_data_len_ = 0;
  std::size_t unpunctured_total_len_ = 0;
  std::vector<bool> keep_mask_;

  mutable std::mutex matrix_mutex_;
  mutable std::optional<matrix_type> generator_matrix_;
  mutable std::optional<matrix_type> parity_check_matrix_;
};

// Fast O(K) Turbo encoder.
class TurboEncoder : public Encoder<GaloisFieldBinary> {
 public:
  using Base = Encoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;

  explicit TurboEncoder(const TurboCode &code)
      : Base(code), turbo_code_(&code) {}

  codeword_type Encode(const message_type &message) const override {
    return turbo_code_->EncodeFast(message);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    return turbo_code_->UnencodeFast(codeword);
  }

  std::string Name() const override { return "Turbo"; }

 private:
  const TurboCode *turbo_code_;
};

// Hard-decision wrapper around TurboCode::DecodeSoft for the AbstractCode
// Decoder<GaloisFieldBinary> interface.
class TurboHardDecoder : public Decoder<GaloisFieldBinary> {
 public:
  using Base = Decoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using result_type = DecodingResult<GaloisFieldBinary>;

  explicit TurboHardDecoder(const TurboCode &code,
                            std::size_t max_iterations = 8,
                            BcjrMode mode = BcjrMode::kLogMap,
                            double hard_llr_magnitude = 4.0)
      : Base(code),
        turbo_code_(&code),
        max_iterations_(max_iterations),
        mode_(mode),
        hard_llr_magnitude_(hard_llr_magnitude) {
    if (max_iterations_ == 0) {
      throw std::invalid_argument(
          "TurboHardDecoder: max_iterations must be at least 1");
    }
    if (!std::isfinite(hard_llr_magnitude_) || hard_llr_magnitude_ <= 0.0) {
      throw std::invalid_argument(
          "TurboHardDecoder: hard_llr_magnitude must be positive and finite");
    }
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "TurboHardDecoder::DecodeWithResult: field mismatch");
    }
    const std::size_t n = turbo_code_->Length();
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "TurboHardDecoder::DecodeWithResult: length mismatch");
    }
    const auto rx = received_word.RawSpan();
    std::vector<double> llrs(n);
    for (std::size_t i = 0; i < n; ++i) {
      llrs[i] =
          ((rx[i] & 1U) == 0U) ? hard_llr_magnitude_ : -hard_llr_magnitude_;
    }
    auto turbo_res = turbo_code_->DecodeSoft(llrs, max_iterations_, mode_, 1.0);
    std::size_t num_errors = 0;
    const auto dec = turbo_res.codeword.RawSpan();
    for (std::size_t i = 0; i < n; ++i) {
      if ((dec[i] & 1U) != (rx[i] & 1U)) {
        ++num_errors;
      }
    }
    return result_type::Success(std::move(turbo_res.codeword),
                                std::move(turbo_res.message), num_errors);
  }

  std::string Name() const override { return "TurboHard"; }
  std::size_t DecodingRadius() const override { return 0; }

 private:
  const TurboCode *turbo_code_;
  std::size_t max_iterations_;
  BcjrMode mode_;
  double hard_llr_magnitude_;
};

inline void TurboCode::RegisterEncodersAndDecoders() {
  this->RegisterEncoder(
      "Turbo", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &turbo = dynamic_cast<const TurboCode &>(code);
        return std::make_shared<TurboEncoder>(turbo);
      });
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<GaloisFieldBinary> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<GaloisFieldBinary>>(lin);
      });
  this->RegisterDecoder(
      "TurboHard", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &turbo = dynamic_cast<const TurboCode &>(code);
        return std::make_shared<TurboHardDecoder>(turbo);
      });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_TURBO_HPP_
