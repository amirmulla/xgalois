#ifndef XGALOIS_CHANNEL_ASYMMETRIC_SYNC_CHANNEL_HPP_
#define XGALOIS_CHANNEL_ASYMMETRIC_SYNC_CHANNEL_HPP_

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/channel/rng.hpp"
#include "xgalois/field/gf_array.hpp"

namespace xg {
namespace channel {

namespace detail {

inline double BinaryEntropy(double p) noexcept {
  if (p <= 0.0 || p >= 1.0) {
    return 0.0;
  }
  return -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
}

}  // namespace detail

// Binary Asymmetric Channel (BAC) with crossover probabilities:
//   p01 = P(Y = 1 | X = 0),  p10 = P(Y = 0 | X = 1).
// Includes Z-channel factory (p01 = 0, p10 = p), exact clamped LLRs matching
// L = ln(P(X=0 | Y=y) / P(X=1 | Y=y)), and exact closed-form Shannon capacity.
class BinaryAsymmetricChannel {
 public:
  BinaryAsymmetricChannel(double p01, double p10, uint64_t seed = 42)
      : p01_(p01), p10_(p10), rng_(seed) {
    if (std::isnan(p01_) || p01_ < 0.0 || p01_ > 1.0) {
      throw std::invalid_argument(
          "BinaryAsymmetricChannel: p01 must be in [0, 1]");
    }
    if (std::isnan(p10_) || p10_ < 0.0 || p10_ > 1.0) {
      throw std::invalid_argument(
          "BinaryAsymmetricChannel: p10 must be in [0, 1]");
    }
  }

  // Z-channel: 0 -> 1 never occurs (p01 = 0), 1 -> 0 occurs with prob p10.
  static BinaryAsymmetricChannel ZChannel(double p10, uint64_t seed = 42) {
    return BinaryAsymmetricChannel(0.0, p10, seed);
  }

  double P01() const noexcept { return p01_; }
  double P10() const noexcept { return p10_; }

  void SetSeed(uint64_t seed) noexcept { rng_.Seed(seed); }

  BinaryAsymmetricChannel CloneWithSeed(uint64_t seed) const {
    return BinaryAsymmetricChannel(p01_, p10_, seed);
  }

  std::vector<uint8_t> Transmit(std::span<const uint8_t> bits) {
    std::vector<uint8_t> out(bits.size(), 0);
    for (std::size_t i = 0; i < bits.size(); ++i) {
      const uint8_t b = bits[i];
      if (b > 1) {
        throw std::invalid_argument(
            "BinaryAsymmetricChannel::Transmit: bit values must be 0 or 1");
      }
      if (b == 0) {
        const bool flip =
            (p01_ >= 1.0) || (p01_ > 0.0 && rng_.NextUniform01() < p01_);
        out[i] = flip ? uint8_t{1} : uint8_t{0};
      } else {
        const bool flip =
            (p10_ >= 1.0) || (p10_ > 0.0 && rng_.NextUniform01() < p10_);
        out[i] = flip ? uint8_t{0} : uint8_t{1};
      }
    }
    return out;
  }

  // Exact channel LLRs under uniform input prior:
  //   L(y) = ln(P(Y = y | X = 0) / P(Y = y | X = 1)),
  // clamped to [-max_llr, +max_llr]. Positive => bit 0, negative => bit 1.
  std::vector<double> ComputeLlrs(std::span<const uint8_t> rx_bits,
                                  double max_llr = 50.0) const {
    if (std::isnan(max_llr) || max_llr < 0.0) {
      throw std::invalid_argument(
          "BinaryAsymmetricChannel::ComputeLlrs: max_llr must be >= 0");
    }

    const double llr_y0 = ComputeSingleLlr(1.0 - p01_, p10_, max_llr);
    const double llr_y1 = ComputeSingleLlr(p01_, 1.0 - p10_, max_llr);

    std::vector<double> llrs(rx_bits.size(), 0.0);
    for (std::size_t i = 0; i < rx_bits.size(); ++i) {
      const uint8_t y = rx_bits[i];
      if (y > 1) {
        throw std::invalid_argument(
            "BinaryAsymmetricChannel::ComputeLlrs: rx bit values must be 0 or "
            "1");
      }
      llrs[i] = (y == 0) ? llr_y0 : llr_y1;
    }
    return llrs;
  }

  // Exact Shannon capacity C(p01, p10) in bits per channel use:
  //   C = max_{q in [0,1]} [ H_2(a + (1 - a - b) q) - (1 - q) H_2(a) - q H_2(b)
  //   ]
  // where (a, b) = (p01, p10) if p01 + p10 < 1, or (1 - p01, 1 - p10) if > 1.
  double Capacity() const noexcept {
    double a = p01_;
    double b = p10_;
    const double sum = a + b;
    if (std::abs(sum - 1.0) <= 1e-15) {
      return 0.0;
    }
    if (sum > 1.0) {
      a = 1.0 - a;
      b = 1.0 - b;
    }
    const double span = 1.0 - a - b;
    if (span <= 0.0) {
      return 0.0;
    }

    const double ha = detail::BinaryEntropy(a);
    const double hb = detail::BinaryEntropy(b);
    const double z = (hb - ha) / span;

    // y* = 1 / (1 + 2^z) computed without overflow for large |z|.
    double y_star = 0.5;
    if (z >= 0.0) {
      const double t = std::exp2(-z);
      y_star = t / (1.0 + t);
    } else {
      const double t = std::exp2(z);
      y_star = 1.0 / (1.0 + t);
    }

    const double q_star = std::clamp((y_star - a) / span, 0.0, 1.0);
    const double cap =
        detail::BinaryEntropy(y_star) - (1.0 - q_star) * ha - q_star * hb;
    return std::max(0.0, cap);
  }

 private:
  static double ComputeSingleLlr(double p_y_given_0, double p_y_given_1,
                                 double max_llr) noexcept {
    if (p_y_given_0 == 0.0 && p_y_given_1 == 0.0) {
      return 0.0;
    }
    if (p_y_given_1 == 0.0) {
      return max_llr;
    }
    if (p_y_given_0 == 0.0) {
      return -max_llr;
    }
    return std::clamp(std::log(p_y_given_0 / p_y_given_1), -max_llr, max_llr);
  }

  double p01_;
  double p10_;
  Pcg64 rng_;
};

enum class EditOp : uint8_t {
  kMatch,
  kSubstitute,
  kDelete,
  kInsert,
};

template <typename GaloisField>
struct IdsTransmissionResult {
  GFArray<GaloisField> received_word;
  std::vector<EditOp> edit_trace;
  std::size_t num_insertions = 0;
  std::size_t num_deletions = 0;
  std::size_t num_substitutions = 0;
};

// Davey-MacKay Insertion-Deletion-Substitution (IDS) synchronization channel
// over an arbitrary Galois field F_q.
//
// For each input symbol x_i:
//   1. Up to `max_insertions_per_pos` geometric insertions occur before x_i,
//      each with probability `p_ins`, inserting a uniform random element of
//      F_q. Setting `max_insertions_per_pos = 1` selects single-insertion mode.
//   2. With probability `p_del`, x_i is deleted.
//   3. Otherwise x_i is transmitted, and with probability `p_sub` it is
//      corrupted by adding a uniform non-zero error from F_q \ {0}.
template <typename GaloisField>
class InsertionDeletionSubstitutionChannel {
 public:
  using field_type = GaloisField;
  using storage_type = typename GaloisField::element_type;

  InsertionDeletionSubstitutionChannel(std::shared_ptr<GaloisField> field,
                                       double p_ins, double p_del, double p_sub,
                                       std::size_t max_insertions_per_pos = 1,
                                       uint64_t seed = 42)
      : field_(std::move(field)),
        p_ins_(p_ins),
        p_del_(p_del),
        p_sub_(p_sub),
        max_insertions_per_pos_(max_insertions_per_pos),
        rng_(seed) {
    if (!field_) {
      throw std::invalid_argument(
          "InsertionDeletionSubstitutionChannel: field pointer cannot be null");
    }
    ValidateProbability(p_ins_, "p_ins");
    ValidateProbability(p_del_, "p_del");
    ValidateProbability(p_sub_, "p_sub");
    if (p_ins_ == 1.0 && max_insertions_per_pos_ > 1024) {
      throw std::invalid_argument(
          "InsertionDeletionSubstitutionChannel: p_ins == 1.0 requires a small "
          "finite max_insertions_per_pos (<= 1024)");
    }
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }
  double PIns() const noexcept { return p_ins_; }
  double PDel() const noexcept { return p_del_; }
  double PSub() const noexcept { return p_sub_; }
  std::size_t MaxInsertionsPerPos() const noexcept {
    return max_insertions_per_pos_;
  }

  void SetSeed(uint64_t seed) noexcept { rng_.Seed(seed); }

  InsertionDeletionSubstitutionChannel CloneWithSeed(uint64_t seed) const {
    return InsertionDeletionSubstitutionChannel(field_, p_ins_, p_del_, p_sub_,
                                                max_insertions_per_pos_, seed);
  }

  IdsTransmissionResult<GaloisField> Transmit(
      const GFArray<GaloisField>& codeword) {
    if (!codeword.HasField() || codeword.Field() != field_) {
      throw std::invalid_argument(
          "InsertionDeletionSubstitutionChannel::Transmit: field mismatch");
    }
    const auto in_span = codeword.RawSpan();
    const std::size_t n = in_span.size();

    std::vector<storage_type> out_buf;
    out_buf.reserve(n);

    IdsTransmissionResult<GaloisField> result;
    result.edit_trace.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
      std::size_t inserted = 0;
      while (inserted < max_insertions_per_pos_ && DrawBernoulli(p_ins_)) {
        out_buf.push_back(SampleUniformElement());
        result.edit_trace.push_back(EditOp::kInsert);
        ++result.num_insertions;
        ++inserted;
      }

      if (DrawBernoulli(p_del_)) {
        result.edit_trace.push_back(EditOp::kDelete);
        ++result.num_deletions;
      } else if (DrawBernoulli(p_sub_)) {
        out_buf.push_back(field_->Add(in_span[i], SampleNonzeroError()));
        result.edit_trace.push_back(EditOp::kSubstitute);
        ++result.num_substitutions;
      } else {
        out_buf.push_back(in_span[i]);
        result.edit_trace.push_back(EditOp::kMatch);
      }
    }

    using RawArray = typename GFArray<GaloisField>::raw_array_type;
    using ShapeType = typename GFArray<GaloisField>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{out_buf.size()});
    for (std::size_t i = 0; i < out_buf.size(); ++i) {
      raw(i) = out_buf[i];
    }
    result.received_word =
        GFArray<GaloisField>::FromRaw(field_, std::move(raw));
    return result;
  }

 private:
  static void ValidateProbability(double p, const char* name) {
    if (std::isnan(p) || p < 0.0 || p > 1.0) {
      throw std::invalid_argument(
          std::string("InsertionDeletionSubstitutionChannel: ") + name +
          " must be in [0, 1]");
    }
  }

  bool DrawBernoulli(double p) noexcept {
    if (p <= 0.0) {
      return false;
    }
    if (p >= 1.0) {
      return true;
    }
    return rng_.NextUniform01() < p;
  }

  storage_type SampleUniformElement() {
    const uint64_t order = field_->Order();
    if constexpr (std::is_integral_v<storage_type>) {
      const uint64_t val = (order == 0) ? rng_() : rng_.NextBounded(order);
      return field_->SetElementValue(static_cast<storage_type>(val));
    } else if constexpr (requires(const GaloisField& f, uint64_t v) {
                           {
                             f.FromInt(v)
                           } -> std::convertible_to<storage_type>;
                         }) {
      const uint64_t val = (order == 0) ? rng_() : rng_.NextBounded(order);
      return field_->FromInt(val);
    } else {
      std::mt19937_64 fallback_rng(rng_());
      return field_->Random(fallback_rng);
    }
  }

  storage_type SampleNonzeroError() {
    const uint64_t order = field_->Order();
    const uint64_t max_k =
        (order == 0) ? std::numeric_limits<uint64_t>::max() : (order - 1);
    if (max_k == 0) {
      throw std::logic_error(
          "InsertionDeletionSubstitutionChannel: field order must be >= 2");
    }

    auto sample_k = [this, max_k]() -> uint64_t {
      if (max_k == 1) {
        return 1ULL;
      }
      if (max_k == std::numeric_limits<uint64_t>::max()) {
        uint64_t v = rng_();
        while (v == 0) {
          v = rng_();
        }
        return v;
      }
      return uint64_t{1} + rng_.NextBounded(max_k);
    };

    if constexpr (std::is_integral_v<storage_type>) {
      if (max_k == 1) {
        return field_->MultiplicativeIdentity();
      }
      return field_->SetElementValue(static_cast<storage_type>(sample_k()));
    } else if constexpr (requires(const GaloisField& f, uint64_t v) {
                           {
                             f.FromInt(v)
                           } -> std::convertible_to<storage_type>;
                         }) {
      if (max_k == 1) {
        return field_->MultiplicativeIdentity();
      }
      return field_->FromInt(sample_k());
    } else {
      std::mt19937_64 fallback_rng(rng_());
      storage_type e = field_->Random(fallback_rng);
      while (field_->IsZero(e)) {
        e = field_->Random(fallback_rng);
      }
      return e;
    }
  }

  std::shared_ptr<GaloisField> field_;
  double p_ins_;
  double p_del_;
  double p_sub_;
  std::size_t max_insertions_per_pos_;
  Pcg64 rng_;
};

// Exact Levenshtein edit distance in O(n * m) time and O(min(n, m)) space.
template <typename T>
std::size_t LevenshteinDistance(std::span<const T> a, std::span<const T> b) {
  if (a.size() < b.size()) {
    std::swap(a, b);
  }
  const std::size_t n = a.size();
  const std::size_t m = b.size();
  if (m == 0) {
    return n;
  }

  std::vector<std::size_t> dp(m + 1);
  for (std::size_t j = 0; j <= m; ++j) {
    dp[j] = j;
  }

  for (std::size_t i = 1; i <= n; ++i) {
    std::size_t prev_diag = dp[0];
    dp[0] = i;
    for (std::size_t j = 1; j <= m; ++j) {
      const std::size_t prev_up = dp[j];
      const std::size_t sub_cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
      dp[j] = std::min({prev_up + 1, dp[j - 1] + 1, prev_diag + sub_cost});
      prev_diag = prev_up;
    }
  }
  return dp[m];
}

template <typename Arr>
  requires requires(const Arr& x) {
    typename Arr::field_type;
    typename Arr::element_type;
    x.RawSpan();
  }
std::size_t LevenshteinDistance(const Arr& a, const Arr& b) {
  return LevenshteinDistance(
      std::span<const typename Arr::element_type>(a.RawSpan()),
      std::span<const typename Arr::element_type>(b.RawSpan()));
}

}  // namespace channel
}  // namespace xg

#endif  // XGALOIS_CHANNEL_ASYMMETRIC_SYNC_CHANNEL_HPP_
