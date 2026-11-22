#ifndef XGALOIS_CHANNEL_BURST_CHANNEL_HPP_
#define XGALOIS_CHANNEL_BURST_CHANNEL_HPP_

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/channel/rng.hpp"
#include "xgalois/field/gf_array.hpp"

namespace xg {
namespace channel {

enum class BurstChannelMode : uint8_t { kSymbolError, kSymbolErasure };

// Structured output from GilbertElliottChannel::Transmit and TransmitMatrix.
// `state_trace` records the per-symbol (or per-row for matrices) Markov state:
// 0 = Good, 1 = Bad. `num_burst_transitions` counts Good -> Bad (0 -> 1)
// transitions observed within this transmission.
template <typename GaloisField>
struct BurstTransmissionResult {
  GFArray<GaloisField> received_word;
  std::vector<std::size_t> error_positions;
  std::vector<std::size_t> erasure_positions;
  std::vector<uint8_t> state_trace;
  std::size_t num_burst_transitions = 0;
};

// 2-state Markov Gilbert-Elliott burst error / burst erasure channel over an
// arbitrary Galois field. Supports both 1-D codewords and 2-D n x ell
// Interleaved Reed-Solomon (IRS) codeword matrices with synchronized row
// bursts.
template <typename GaloisField>
class GilbertElliottChannel {
 public:
  using field_type = GaloisField;
  using storage_type = typename GaloisField::element_type;

  static constexpr uint8_t kGoodState = 0;
  static constexpr uint8_t kBadState = 1;

  GilbertElliottChannel(std::shared_ptr<GaloisField> field, double p_gb,
                        double p_bg, double err_good, double err_bad,
                        BurstChannelMode mode = BurstChannelMode::kSymbolError,
                        uint64_t seed = 42)
      : field_(std::move(field)),
        p_gb_(p_gb),
        p_bg_(p_bg),
        err_good_(err_good),
        err_bad_(err_bad),
        mode_(mode),
        rng_(seed) {
    if (!field_) {
      throw std::invalid_argument(
          "GilbertElliottChannel: field pointer cannot be null");
    }
    ValidateProbability(p_gb_, "p_gb");
    ValidateProbability(p_bg_, "p_bg");
    ValidateProbability(err_good_, "err_good");
    ValidateProbability(err_bad_, "err_bad");
    ResetToStationary();
  }

  // Constructs a Gilbert-Elliott channel from desired mean burst length
  // (1 / p_bg >= 1) and mean good-run length (1 / p_gb >= 1).
  static GilbertElliottChannel FromBurstStats(
      std::shared_ptr<GaloisField> field, double avg_burst_length,
      double avg_good_length, double err_good, double err_bad,
      BurstChannelMode mode = BurstChannelMode::kSymbolError,
      uint64_t seed = 42) {
    if (std::isnan(avg_burst_length) || std::isinf(avg_burst_length) ||
        avg_burst_length < 1.0) {
      throw std::invalid_argument(
          "GilbertElliottChannel::FromBurstStats: avg_burst_length must be "
          "finite and >= 1.0");
    }
    if (std::isnan(avg_good_length) || std::isinf(avg_good_length) ||
        avg_good_length < 1.0) {
      throw std::invalid_argument(
          "GilbertElliottChannel::FromBurstStats: avg_good_length must be "
          "finite and >= 1.0");
    }
    const double p_bg = 1.0 / avg_burst_length;
    const double p_gb = 1.0 / avg_good_length;
    return GilbertElliottChannel(std::move(field), p_gb, p_bg, err_good,
                                 err_bad, mode, seed);
  }

  // Analytical stationary probability of occupying the Bad state:
  // pi_B = p_gb / (p_gb + p_bg).
  double SteadyStateBadProb() const noexcept {
    const double denom = p_gb_ + p_bg_;
    if (denom == 0.0) {
      return 0.0;
    }
    return p_gb_ / denom;
  }

  double SteadyStateGoodProb() const noexcept {
    return 1.0 - SteadyStateBadProb();
  }

  // Analytical mean sojourn time in the Bad state: E[L_B] = 1 / p_bg.
  double AverageBurstLength() const noexcept {
    if (p_bg_ == 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    return 1.0 / p_bg_;
  }

  // Analytical mean sojourn time in the Good state: E[L_G] = 1 / p_gb.
  double AverageGoodLength() const noexcept {
    if (p_gb_ == 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    return 1.0 / p_gb_;
  }

  // Analytical stationary symbol error/erasure rate:
  // (1 - pi_B) * err_good + pi_B * err_bad.
  double SteadyStateErrorRate() const noexcept {
    const double pi_b = SteadyStateBadProb();
    return (1.0 - pi_b) * err_good_ + pi_b * err_bad_;
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }
  double PGb() const noexcept { return p_gb_; }
  double PBg() const noexcept { return p_bg_; }
  double ErrGood() const noexcept { return err_good_; }
  double ErrBad() const noexcept { return err_bad_; }
  BurstChannelMode Mode() const noexcept { return mode_; }
  uint8_t CurrentState() const noexcept { return state_; }

  void SetState(uint8_t state) {
    if (state != kGoodState && state != kBadState) {
      throw std::invalid_argument(
          "GilbertElliottChannel::SetState: state must be 0 (Good) or 1 (Bad)");
    }
    state_ = state;
  }

  void SetSeed(uint64_t seed) noexcept {
    rng_.Seed(seed);
    ResetToStationary();
  }

  void ResetToStationary() noexcept {
    const double pi_b = SteadyStateBadProb();
    state_ =
        (pi_b > 0.0 && rng_.NextUniform01() < pi_b) ? kBadState : kGoodState;
  }

  // Returns an independent channel instance with identical Markov parameters
  // initialized with `seed`.
  GilbertElliottChannel CloneWithSeed(uint64_t seed) const {
    return GilbertElliottChannel(field_, p_gb_, p_bg_, err_good_, err_bad_,
                                 mode_, seed);
  }

  // Transmits a 1-D codeword (or flattened symbol stream) through the
  // Gilbert-Elliott channel.
  BurstTransmissionResult<GaloisField> Transmit(
      const GFArray<GaloisField>& codeword) {
    ValidateField(codeword);
    GFArray<GaloisField> received = codeword;
    auto span = received.RawSpan();
    const std::size_t n = span.size();

    BurstTransmissionResult<GaloisField> result;
    result.state_trace.resize(n);
    const storage_type zero = field_->AdditiveIdentity();

    for (std::size_t i = 0; i < n; ++i) {
      const uint8_t cur_state = state_;
      result.state_trace[i] = cur_state;
      if (i > 0 && result.state_trace[i - 1] == kGoodState &&
          cur_state == kBadState) {
        ++result.num_burst_transitions;
      }

      const double err_prob = (cur_state == kGoodState) ? err_good_ : err_bad_;
      if (DrawBernoulli(err_prob)) {
        if (mode_ == BurstChannelMode::kSymbolError) {
          span[i] = field_->Add(span[i], SampleNonzeroError());
          result.error_positions.push_back(i);
        } else {
          span[i] = zero;
          result.erasure_positions.push_back(i);
        }
      }

      AdvanceState();
    }

    result.received_word = std::move(received);
    return result;
  }

  // Transmits an n x ell interleaved codeword matrix (e.g., for Interleaved
  // Reed-Solomon codes). The Markov chain steps once per row (0 .. n - 1), and
  // a burst error or erasure at row i affects the entire row across all ell
  // interleaved columns simultaneously.
  BurstTransmissionResult<GaloisField> TransmitMatrix(
      const GFMatrix<GaloisField>& codeword_matrix) {
    ValidateField(codeword_matrix);
    if (codeword_matrix.Dimension() != 2) {
      throw std::invalid_argument(
          "GilbertElliottChannel::TransmitMatrix: input must be a 2-D matrix");
    }

    GFMatrix<GaloisField> received = codeword_matrix;
    const std::size_t num_rows = received.Shape(0);
    const std::size_t num_cols = received.Shape(1);
    auto span = received.RawSpan();
    const storage_type zero = field_->AdditiveIdentity();

    BurstTransmissionResult<GaloisField> result;
    result.state_trace.resize(num_rows);

    for (std::size_t row = 0; row < num_rows; ++row) {
      const uint8_t cur_state = state_;
      result.state_trace[row] = cur_state;
      if (row > 0 && result.state_trace[row - 1] == kGoodState &&
          cur_state == kBadState) {
        ++result.num_burst_transitions;
      }

      const double err_prob = (cur_state == kGoodState) ? err_good_ : err_bad_;
      if (DrawBernoulli(err_prob)) {
        const std::size_t row_offset = row * num_cols;
        if (mode_ == BurstChannelMode::kSymbolError) {
          for (std::size_t col = 0; col < num_cols; ++col) {
            span[row_offset + col] =
                field_->Add(span[row_offset + col], SampleNonzeroError());
          }
          result.error_positions.push_back(row);
        } else {
          for (std::size_t col = 0; col < num_cols; ++col) {
            span[row_offset + col] = zero;
          }
          result.erasure_positions.push_back(row);
        }
      }

      AdvanceState();
    }

    result.received_word = std::move(received);
    return result;
  }

 private:
  static void ValidateProbability(double p, const char* name) {
    if (std::isnan(p) || p < 0.0 || p > 1.0) {
      throw std::invalid_argument(std::string("GilbertElliottChannel: ") +
                                  name + " must be in [0, 1]");
    }
  }

  void ValidateField(const GFArray<GaloisField>& arr) const {
    if (!arr.HasField() || arr.Field() != field_) {
      throw std::invalid_argument(
          "GilbertElliottChannel: input field does not match channel field");
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

  void AdvanceState() noexcept {
    if (state_ == kGoodState) {
      if (DrawBernoulli(p_gb_)) {
        state_ = kBadState;
      }
    } else {
      if (DrawBernoulli(p_bg_)) {
        state_ = kGoodState;
      }
    }
  }

  storage_type SampleNonzeroError() {
    const uint64_t order = field_->Order();
    const uint64_t max_k =
        (order == 0) ? std::numeric_limits<uint64_t>::max() : (order - 1);
    if (max_k == 0) {
      throw std::logic_error(
          "GilbertElliottChannel: field order must be at least 2");
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
  double p_gb_;
  double p_bg_;
  double err_good_;
  double err_bad_;
  BurstChannelMode mode_;
  Pcg64 rng_;
  uint8_t state_ = kGoodState;
};

}  // namespace channel
}  // namespace xg

#endif  // XGALOIS_CHANNEL_BURST_CHANNEL_HPP_
