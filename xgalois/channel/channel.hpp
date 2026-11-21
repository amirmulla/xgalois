#ifndef XGALOIS_CHANNEL_CHANNEL_HPP_
#define XGALOIS_CHANNEL_CHANNEL_HPP_

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xtensor/containers/xarray.hpp"

namespace xg {

class GaloisFieldBinary;

namespace channel {

namespace detail {

inline uint64_t GenerateRandomSeed() {
  std::random_device rd;
  return (static_cast<uint64_t>(rd()) << 32) ^ static_cast<uint64_t>(rd());
}

// Numerically stable binomial PMF P(X = k) for X ~ Binomial(n, p),
// computed in log-space via std::lgamma to avoid overflow / inf*0 = NaN.
inline double BinomialPMF(size_t n, size_t k, double p) {
  if (k > n) {
    return 0.0;
  }
  if (p == 0.0) {
    return (k == 0) ? 1.0 : 0.0;
  }
  if (p == 1.0) {
    return (k == n) ? 1.0 : 0.0;
  }
  const double dn = static_cast<double>(n);
  const double dk = static_cast<double>(k);
  const double dnk = static_cast<double>(n - k);
  const double log_binom =
      std::lgamma(dn + 1.0) - std::lgamma(dk + 1.0) - std::lgamma(dnk + 1.0);
  const double log_prob = log_binom + dk * std::log(p) + dnk * std::log1p(-p);
  return std::exp(log_prob);
}

inline double BinomialCDF(size_t n, size_t k, double p) {
  if (k >= n) {
    return 1.0;
  }
  double prob = 0.0;
  for (size_t i = 0; i <= k; ++i) {
    prob += BinomialPMF(n, i, p);
  }
  return std::min(prob, 1.0);
}

template <typename GaloisField>
size_t BitsPerSymbol(const GaloisField &field) {
  const uint64_t order = field.Order();
  if (order == 0) {
    return 64;  // GF(2^64)
  }
  if (order <= 2) {
    return 1;
  }
  return static_cast<size_t>(std::bit_width(order - 1));
}

template <typename GaloisField>
size_t CountSymbolBitErrors(const GaloisField &field,
                            const typename GaloisField::element_type &a,
                            const typename GaloisField::element_type &b,
                            size_t bits_per_symbol) {
  using Storage = typename GaloisField::element_type;
  if (a == b) {
    return 0;
  }
  if constexpr (std::is_integral_v<Storage>) {
    const uint64_t va = static_cast<uint64_t>(field.GetElementValue(a));
    const uint64_t vb = static_cast<uint64_t>(field.GetElementValue(b));
    return static_cast<size_t>(std::popcount(va ^ vb));
  } else if constexpr (requires(const GaloisField &f, const Storage &x) {
                         { f.ToInt(x) } -> std::convertible_to<uint64_t>;
                       }) {
    const uint64_t va = static_cast<uint64_t>(field.ToInt(a));
    const uint64_t vb = static_cast<uint64_t>(field.ToInt(b));
    return static_cast<size_t>(std::popcount(va ^ vb));
  } else {
    return bits_per_symbol;
  }
}

}  // namespace detail

// Structured result returned by erasure-capable channels for GFArray inputs.
template <typename GaloisField>
struct ErasureTransmission {
  GFArray<GaloisField> received;
  std::vector<size_t> erasure_positions;
  xt::xarray<bool> erasure_mask;
};

// Structured result returned by AWGNChannel::transmit_detailed().
template <typename BinaryField = GaloisFieldBinary>
struct AWGNTransmission {
  std::vector<double> soft_received;
  std::vector<double> llr;
  GFArray<BinaryField> hard_decision;
};

template <typename GaloisField>
class Channel {
 public:
  using field_type = GaloisField;
  using storage_type = typename GaloisField::element_type;
  using ArrayType = GFArray<GaloisField>;
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  Channel(std::shared_ptr<GaloisField> field, size_t input_size,
          size_t output_size = 0, std::optional<uint64_t> seed = std::nullopt)
      : field_(std::move(field)),
        input_size_(input_size),
        output_size_(output_size == 0 ? input_size : output_size),
        rng_(seed.has_value() ? *seed : detail::GenerateRandomSeed()) {
    if (!field_) {
      throw std::invalid_argument("Channel: field pointer cannot be null");
    }
  }

  virtual ~Channel() = default;

  std::shared_ptr<GaloisField> field() const { return field_; }
  std::shared_ptr<GaloisField> Field() const { return field_; }

  size_t input_size() const noexcept { return input_size_; }
  size_t InputSize() const noexcept { return input_size_; }

  size_t output_size() const noexcept { return output_size_; }
  size_t OutputSize() const noexcept { return output_size_; }

  void SetSeed(uint64_t seed) { rng_.seed(seed); }
  void set_seed(uint64_t seed) { rng_.seed(seed); }

  // Primary GFArray API (zero per-element shared_ptr overhead).
  ArrayType transmit(const ArrayType &message) {
    ValidateMessage(message);
    return transmit_unsafe(message);
  }

  ArrayType Transmit(const ArrayType &message) { return transmit(message); }

  virtual ArrayType transmit_unsafe(const ArrayType &message) = 0;

  ArrayType operator()(const ArrayType &message) { return transmit(message); }

  // Backwards-compatible legacy xt::xarray<GaloisFieldElement<GaloisField>> API
  // delegating cleanly to the GFArray core.
  VectorType transmit(const VectorType &message) {
    ValidateLegacyMessageSize(message);
    return transmit_unsafe(ToGFArray(message)).ToElementArray();
  }

  VectorType Transmit(const VectorType &message) { return transmit(message); }

  VectorType transmit_unsafe(const VectorType &message) {
    return transmit_unsafe(ToGFArray(message)).ToElementArray();
  }

  VectorType operator()(const VectorType &message) { return transmit(message); }

  virtual void print(std::ostream &os) const = 0;

  void Print(std::ostream &os) const { print(os); }

  virtual std::string to_string() const {
    std::ostringstream oss;
    print(oss);
    return oss.str();
  }

  std::string ToString() const { return to_string(); }

 protected:
  std::shared_ptr<GaloisField> field_;
  size_t input_size_;
  size_t output_size_;
  std::mt19937_64 rng_;

  void ValidateMessage(const ArrayType &message) const {
    if (!message.HasField() || message.Field() != field_) {
      throw std::invalid_argument(
          "Channel: message field does not match channel field");
    }
    if (message.Size() != input_size_) {
      throw std::invalid_argument("Message size (" +
                                  std::to_string(message.Size()) +
                                  ") does not match channel input size (" +
                                  std::to_string(input_size_) + ")");
    }
  }

  void ValidateLegacyMessageSize(const VectorType &message) const {
    if (message.size() != input_size_) {
      throw std::invalid_argument("Message size (" +
                                  std::to_string(message.size()) +
                                  ") does not match channel input size (" +
                                  std::to_string(input_size_) + ")");
    }
  }

  ArrayType ToGFArray(const VectorType &message) const {
    if (message.size() == 0) {
      return ArrayType::Zeros(field_, {0});
    }
    typename ArrayType::raw_array_type raw =
        ArrayType::raw_array_type::from_shape(message.shape());
    const auto *src = message.data();
    auto *dst = raw.data();
    const size_t n = message.size();
    for (size_t i = 0; i < n; ++i) {
      if (!src[i].HasField() || src[i].Field() != field_) {
        throw std::invalid_argument(
            "Channel: message element field does not match channel field");
      }
      dst[i] = src[i].RawValue();
    }
    return ArrayType::FromRaw(field_, std::move(raw));
  }

  // O(1) deterministic uniform non-zero field element sampler.
  // Adding the returned offset e to any symbol x produces a symbol uniformly
  // distributed over F_q \ {x} without rejection sampling.
  storage_type sample_nonzero_error() {
    const uint64_t order = field_->Order();
    // For GF(2^64), Order() returns 0; non-zero external values span
    // [1, UINT64_MAX].
    const uint64_t max_k =
        (order == 0) ? std::numeric_limits<uint64_t>::max() : (order - 1);
    if (max_k == 0) {
      throw std::logic_error("Channel: field order must be at least 2");
    }

    if constexpr (std::is_integral_v<storage_type>) {
      if (max_k == 1) {
        return field_->MultiplicativeIdentity();
      }
      std::uniform_int_distribution<uint64_t> dist(1, max_k);
      const uint64_t k = dist(rng_);
      return field_->SetElementValue(static_cast<storage_type>(k));
    } else if constexpr (requires(const GaloisField &f, uint64_t v) {
                           {
                             f.FromInt(v)
                           } -> std::convertible_to<storage_type>;
                         }) {
      if (max_k == 1) {
        return field_->MultiplicativeIdentity();
      }
      std::uniform_int_distribution<uint64_t> dist(1, max_k);
      const uint64_t k = dist(rng_);
      return field_->FromInt(k);
    } else {
      storage_type e = field_->Random(rng_);
      while (field_->IsZero(e)) {
        e = field_->Random(rng_);
      }
      return e;
    }
  }

  // Samples `count` distinct positions uniformly without replacement in
  // O(count) time and memory when count << input_size_ (Floyd's algorithm +
  // O(count) permutation shuffle), falling back to partial Fisher-Yates when
  // count is a substantial fraction of input_size_.
  std::vector<size_t> sample_distinct_positions(size_t count) {
    if (count > input_size_) {
      throw std::invalid_argument(
          "Number of positions cannot exceed input size");
    }
    std::vector<size_t> result;
    result.reserve(count);
    if (count == 0) {
      return result;
    }

    if (count * 4 >= input_size_) {
      std::vector<size_t> pool(input_size_);
      std::iota(pool.begin(), pool.end(), size_t{0});
      for (size_t i = 0; i < count; ++i) {
        std::uniform_int_distribution<size_t> dist(i, input_size_ - 1);
        const size_t j = dist(rng_);
        std::swap(pool[i], pool[j]);
        result.push_back(pool[i]);
      }
      return result;
    }

    std::unordered_set<size_t> chosen;
    chosen.reserve(count * 2);
    for (size_t j = input_size_ - count; j < input_size_; ++j) {
      std::uniform_int_distribution<size_t> dist(0, j);
      const size_t t = dist(rng_);
      if (chosen.insert(t).second) {
        result.push_back(t);
      } else {
        chosen.insert(j);
        result.push_back(j);
      }
    }
    // Shuffle the `count` chosen indices so order is also uniformly random
    // (critical when partitioning into error vs erasure subsets).
    for (size_t i = count; i > 1; --i) {
      std::uniform_int_distribution<size_t> dist(0, i - 1);
      std::swap(result[i - 1], result[dist(rng_)]);
    }
    return result;
  }

  std::vector<size_t> generate_error_positions(size_t num_errors) {
    return sample_distinct_positions(num_errors);
  }

  ArrayType generate_error_array(const std::vector<size_t> &error_positions) {
    ArrayType error_array = ArrayType::Zeros(field_, {input_size_});
    auto span = error_array.RawSpan();
    for (size_t pos : error_positions) {
      if (pos < input_size_) {
        span[pos] = sample_nonzero_error();
      }
    }
    return error_array;
  }

  VectorType generate_error_vector(const std::vector<size_t> &error_positions) {
    return generate_error_array(error_positions).ToElementArray();
  }

  std::string format_interval(size_t min_val, size_t max_val) const {
    if (min_val == max_val) {
      return std::to_string(min_val);
    }
    return "between " + std::to_string(min_val) + " and " +
           std::to_string(max_val);
  }

  void print_field(std::ostream &os) const {
    if constexpr (std::is_base_of_v<GaloisFieldBase<storage_type>,
                                    GaloisField>) {
      static_cast<const GaloisFieldBase<storage_type> &>(*field_).Print(os);
    } else {
      field_->Print(os);
    }
  }
};

template <typename GaloisField>
class StaticErrorRateChannel : public Channel<GaloisField> {
 public:
  using ArrayType = GFArray<GaloisField>;
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  using Channel<GaloisField>::operator();
  using Channel<GaloisField>::transmit;
  using Channel<GaloisField>::Transmit;
  using Channel<GaloisField>::transmit_unsafe;

  StaticErrorRateChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                         size_t num_errors,
                         std::optional<uint64_t> seed = std::nullopt)
      : StaticErrorRateChannel(std::move(field), input_size, num_errors,
                               num_errors, seed) {}

  StaticErrorRateChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                         size_t min_errors, size_t max_errors,
                         std::optional<uint64_t> seed = std::nullopt)
      : Channel<GaloisField>(std::move(field), input_size, 0, seed),
        min_errors_(min_errors),
        max_errors_(max_errors) {
    if (min_errors_ > max_errors_) {
      throw std::invalid_argument(
          "StaticErrorRateChannel: min_errors cannot be greater than "
          "max_errors");
    }
    if (max_errors_ > this->input_size_) {
      throw std::invalid_argument(
          "StaticErrorRateChannel: max_errors cannot exceed input_size");
    }
  }

  std::pair<size_t, size_t> number_errors() const noexcept {
    return {min_errors_, max_errors_};
  }

  ArrayType transmit_unsafe(const ArrayType &message) override {
    size_t num_errors = min_errors_;
    if (min_errors_ != max_errors_) {
      std::uniform_int_distribution<size_t> dist(min_errors_, max_errors_);
      num_errors = dist(this->rng_);
    }

    ArrayType result = message;
    if (num_errors == 0) {
      return result;
    }

    const size_t clamped_errors = std::min(num_errors, this->input_size_);
    const auto error_positions =
        this->sample_distinct_positions(clamped_errors);
    auto span = result.RawSpan();
    const size_t msg_len = span.size();
    for (size_t pos : error_positions) {
      if (pos < msg_len) {
        span[pos] = this->field_->Add(span[pos], this->sample_nonzero_error());
      }
    }
    return result;
  }

  void print(std::ostream &os) const override {
    os << "Static error rate channel creating "
       << this->format_interval(min_errors_, max_errors_)
       << " errors, of input and output space vector of dimension "
       << this->input_size_ << " over ";
    this->print_field(os);
  }

 private:
  size_t min_errors_;
  size_t max_errors_;
};

template <typename GaloisField>
class ErrorErasureChannel : public Channel<GaloisField> {
 public:
  using ArrayType = GFArray<GaloisField>;
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  using Channel<GaloisField>::operator();
  using Channel<GaloisField>::transmit;
  using Channel<GaloisField>::Transmit;
  using Channel<GaloisField>::transmit_unsafe;

  ErrorErasureChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                      size_t num_errors, size_t num_erasures,
                      std::optional<uint64_t> seed = std::nullopt)
      : ErrorErasureChannel(std::move(field), input_size, num_errors,
                            num_errors, num_erasures, num_erasures, seed) {}

  ErrorErasureChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                      size_t min_errors, size_t max_errors, size_t min_erasures,
                      size_t max_erasures,
                      std::optional<uint64_t> seed = std::nullopt)
      : Channel<GaloisField>(std::move(field), input_size, 0, seed),
        min_errors_(min_errors),
        max_errors_(max_errors),
        min_erasures_(min_erasures),
        max_erasures_(max_erasures) {
    if (min_errors_ > max_errors_ || min_erasures_ > max_erasures_) {
      throw std::invalid_argument(
          "ErrorErasureChannel: min values cannot be greater than max values");
    }
    if (max_errors_ > this->input_size_) {
      throw std::invalid_argument(
          "ErrorErasureChannel: max_errors cannot exceed input_size");
    }
    if (max_erasures_ > this->input_size_) {
      throw std::invalid_argument(
          "ErrorErasureChannel: max_erasures cannot exceed input_size");
    }
    // Overflow-safe check: max_errors_ + max_erasures_ <= input_size_.
    if (max_errors_ > this->input_size_ - max_erasures_) {
      throw std::invalid_argument(
          "ErrorErasureChannel: max_errors + max_erasures cannot exceed "
          "input_size");
    }
  }

  std::pair<size_t, size_t> number_errors() const noexcept {
    return {min_errors_, max_errors_};
  }

  std::pair<size_t, size_t> number_erasures() const noexcept {
    return {min_erasures_, max_erasures_};
  }

  // Primary safe GFArray API returning {received, erasure_positions,
  // erasure_mask}.
  ErasureTransmission<GaloisField> transmit_with_erasures(
      const ArrayType &message) {
    this->ValidateMessage(message);
    return transmit_unsafe_with_erasures(message);
  }

  ErasureTransmission<GaloisField> TransmitWithErasures(
      const ArrayType &message) {
    return transmit_with_erasures(message);
  }

  ErasureTransmission<GaloisField> transmit_unsafe_with_erasures(
      const ArrayType &message) {
    const size_t n = this->input_size_;
    size_t num_errors = min_errors_;
    if (min_errors_ != max_errors_) {
      std::uniform_int_distribution<size_t> dist(min_errors_, max_errors_);
      num_errors = dist(this->rng_);
    }
    size_t num_erasures = min_erasures_;
    if (min_erasures_ != max_erasures_) {
      std::uniform_int_distribution<size_t> dist(min_erasures_, max_erasures_);
      num_erasures = dist(this->rng_);
    }

    // Defensive clamping guarantees 0 <= first <= last <= n ALWAYS holds.
    const size_t first = std::min(num_errors, n);
    const size_t last = first + std::min(num_erasures, n - first);

    const auto chosen = this->sample_distinct_positions(last);

    ArrayType received = message;
    auto span = received.RawSpan();
    const size_t msg_len = span.size();

    for (size_t i = 0; i < first; ++i) {
      const size_t pos = chosen[i];
      if (pos < msg_len) {
        span[pos] = this->field_->Add(span[pos], this->sample_nonzero_error());
      }
    }

    xt::xarray<bool> erasure_mask = xt::xarray<bool>::from_shape({n});
    std::fill(erasure_mask.begin(), erasure_mask.end(), false);

    std::vector<size_t> erasure_positions;
    erasure_positions.reserve(last - first);
    const auto zero = this->field_->AdditiveIdentity();
    for (size_t i = first; i < last; ++i) {
      const size_t pos = chosen[i];
      erasure_positions.push_back(pos);
      if (pos < n) {
        erasure_mask[pos] = true;
      }
      if (pos < msg_len) {
        span[pos] = zero;
      }
    }
    std::sort(erasure_positions.begin(), erasure_positions.end());

    return ErasureTransmission<GaloisField>{std::move(received),
                                            std::move(erasure_positions),
                                            std::move(erasure_mask)};
  }

  // Legacy xt::xarray<GaloisFieldElement<GaloisField>> overloads.
  std::pair<VectorType, xt::xarray<bool>> transmit_with_erasures(
      const VectorType &message) {
    this->ValidateLegacyMessageSize(message);
    auto res = transmit_unsafe_with_erasures(this->ToGFArray(message));
    return {res.received.ToElementArray(), std::move(res.erasure_mask)};
  }

  std::pair<VectorType, xt::xarray<bool>> transmit_unsafe_with_erasures(
      const VectorType &message) {
    auto res = transmit_unsafe_with_erasures(this->ToGFArray(message));
    return {res.received.ToElementArray(), std::move(res.erasure_mask)};
  }

  ArrayType transmit_unsafe(const ArrayType &message) override {
    return transmit_unsafe_with_erasures(message).received;
  }

  void print(std::ostream &os) const override {
    os << "Error-and-erasure channel creating "
       << this->format_interval(min_errors_, max_errors_) << " errors and "
       << this->format_interval(min_erasures_, max_erasures_)
       << " erasures of input space vector of dimension " << this->input_size_
       << " over ";
    this->print_field(os);
  }

 private:
  size_t min_errors_;
  size_t max_errors_;
  size_t min_erasures_;
  size_t max_erasures_;
};

template <typename GaloisField>
class QarySymmetricChannel : public Channel<GaloisField> {
 public:
  using ArrayType = GFArray<GaloisField>;
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  using Channel<GaloisField>::operator();
  using Channel<GaloisField>::transmit;
  using Channel<GaloisField>::Transmit;
  using Channel<GaloisField>::transmit_unsafe;

  QarySymmetricChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                       double epsilon,
                       std::optional<uint64_t> seed = std::nullopt)
      : Channel<GaloisField>(std::move(field), input_size, 0, seed),
        epsilon_(epsilon) {
    if (std::isnan(epsilon_) || epsilon_ < 0.0 || epsilon_ > 1.0) {
      throw std::invalid_argument(
          "QarySymmetricChannel: error probability must be between 0 and 1");
    }
  }

  double error_probability() const noexcept { return epsilon_; }

  double probability_of_exactly_t_errors(size_t t) const {
    return detail::BinomialPMF(this->input_size_, t, epsilon_);
  }

  double probability_of_at_most_t_errors(size_t t) const {
    return detail::BinomialCDF(this->input_size_, t, epsilon_);
  }

  ArrayType transmit_unsafe(const ArrayType &message) override {
    ArrayType result = message;
    if (epsilon_ == 0.0) {
      return result;
    }

    auto span = result.RawSpan();
    if (epsilon_ == 1.0) {
      for (size_t i = 0; i < span.size(); ++i) {
        span[i] = this->field_->Add(span[i], this->sample_nonzero_error());
      }
      return result;
    }

    std::bernoulli_distribution error_dist(epsilon_);
    for (size_t i = 0; i < span.size(); ++i) {
      if (error_dist(this->rng_)) {
        span[i] = this->field_->Add(span[i], this->sample_nonzero_error());
      }
    }
    return result;
  }

  void print(std::ostream &os) const override {
    os << "q-ary symmetric channel with error probability " << epsilon_
       << ", of input and output space vector of dimension "
       << this->input_size_ << " over ";
    this->print_field(os);
  }

 private:
  double epsilon_;
};

// Specialized Binary Symmetric Channel (BSC) over GF(2).
template <typename BinaryField = GaloisFieldBinary>
class BasicBinarySymmetricChannel : public Channel<BinaryField> {
 public:
  using ArrayType = GFArray<BinaryField>;
  using ElementType = GaloisFieldElementBase<BinaryField>;
  using VectorType = xt::xarray<ElementType>;

  using Channel<BinaryField>::operator();
  using Channel<BinaryField>::transmit;
  using Channel<BinaryField>::Transmit;
  using Channel<BinaryField>::transmit_unsafe;

  BasicBinarySymmetricChannel(std::shared_ptr<BinaryField> field,
                              size_t input_size, double crossover_prob,
                              std::optional<uint64_t> seed = std::nullopt)
      : Channel<BinaryField>(std::move(field), input_size, 0, seed),
        p_(crossover_prob) {
    if (this->field_->Order() != 2) {
      throw std::invalid_argument(
          "BinarySymmetricChannel requires a binary field GF(2)");
    }
    if (std::isnan(p_) || p_ < 0.0 || p_ > 1.0) {
      throw std::invalid_argument(
          "BinarySymmetricChannel: crossover probability must be in [0, 1]");
    }
  }

  BasicBinarySymmetricChannel(size_t input_size, double crossover_prob,
                              std::optional<uint64_t> seed = std::nullopt)
      : BasicBinarySymmetricChannel(std::make_shared<BinaryField>(), input_size,
                                    crossover_prob, seed) {}

  double crossover_probability() const noexcept { return p_; }
  double error_probability() const noexcept { return p_; }

  // Exact Shannon capacity C(p) = 1 - H_2(p) in bits per channel use.
  double capacity() const {
    if (p_ == 0.0 || p_ == 1.0) {
      return 1.0;
    }
    const double h2 = -p_ * std::log2(p_) - (1.0 - p_) * std::log2(1.0 - p_);
    return std::max(0.0, 1.0 - h2);
  }

  double probability_of_exactly_t_errors(size_t t) const {
    return detail::BinomialPMF(this->input_size_, t, p_);
  }

  double probability_of_at_most_t_errors(size_t t) const {
    return detail::BinomialCDF(this->input_size_, t, p_);
  }

  ArrayType transmit_unsafe(const ArrayType &message) override {
    ArrayType result = message;
    if (p_ == 0.0) {
      return result;
    }

    auto span = result.RawSpan();
    const auto one = this->field_->MultiplicativeIdentity();
    if (p_ == 1.0) {
      for (size_t i = 0; i < span.size(); ++i) {
        span[i] = this->field_->Add(span[i], one);
      }
      return result;
    }

    std::bernoulli_distribution flip_dist(p_);
    for (size_t i = 0; i < span.size(); ++i) {
      if (flip_dist(this->rng_)) {
        span[i] = this->field_->Add(span[i], one);
      }
    }
    return result;
  }

  void print(std::ostream &os) const override {
    os << "Binary symmetric channel with crossover probability " << p_
       << ", of input and output space vector of dimension "
       << this->input_size_ << " over ";
    this->print_field(os);
  }

 private:
  double p_;
};

using BinarySymmetricChannel = BasicBinarySymmetricChannel<GaloisFieldBinary>;

// Memoryless q-ary Erasure Channel (BEC when GaloisField is GaloisFieldBinary).
template <typename GaloisField>
class QaryErasureChannel : public Channel<GaloisField> {
 public:
  using ArrayType = GFArray<GaloisField>;
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  using Channel<GaloisField>::operator();
  using Channel<GaloisField>::transmit;
  using Channel<GaloisField>::Transmit;
  using Channel<GaloisField>::transmit_unsafe;

  QaryErasureChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                     double epsilon,
                     std::optional<uint64_t> seed = std::nullopt)
      : Channel<GaloisField>(std::move(field), input_size, 0, seed),
        epsilon_(epsilon) {
    if (std::isnan(epsilon_) || epsilon_ < 0.0 || epsilon_ > 1.0) {
      throw std::invalid_argument(
          "QaryErasureChannel: erasure probability must be in [0, 1]");
    }
  }

  double erasure_probability() const noexcept { return epsilon_; }

  double probability_of_exactly_s_erasures(size_t s) const {
    return detail::BinomialPMF(this->input_size_, s, epsilon_);
  }

  double probability_of_at_most_s_erasures(size_t s) const {
    return detail::BinomialCDF(this->input_size_, s, epsilon_);
  }

  ErasureTransmission<GaloisField> transmit_with_erasures(
      const ArrayType &message) {
    this->ValidateMessage(message);
    return transmit_unsafe_with_erasures(message);
  }

  ErasureTransmission<GaloisField> TransmitWithErasures(
      const ArrayType &message) {
    return transmit_with_erasures(message);
  }

  ErasureTransmission<GaloisField> transmit_unsafe_with_erasures(
      const ArrayType &message) {
    ArrayType received = message;
    auto span = received.RawSpan();
    const size_t n = span.size();
    const auto zero = this->field_->AdditiveIdentity();

    xt::xarray<bool> erasure_mask = xt::xarray<bool>::from_shape({n});
    std::fill(erasure_mask.begin(), erasure_mask.end(), false);
    std::vector<size_t> erasure_positions;

    if (epsilon_ == 1.0) {
      erasure_positions.reserve(n);
      for (size_t i = 0; i < n; ++i) {
        span[i] = zero;
        erasure_mask[i] = true;
        erasure_positions.push_back(i);
      }
    } else if (epsilon_ > 0.0) {
      std::bernoulli_distribution erase_dist(epsilon_);
      for (size_t i = 0; i < n; ++i) {
        if (erase_dist(this->rng_)) {
          span[i] = zero;
          erasure_mask[i] = true;
          erasure_positions.push_back(i);
        }
      }
    }

    return ErasureTransmission<GaloisField>{std::move(received),
                                            std::move(erasure_positions),
                                            std::move(erasure_mask)};
  }

  std::pair<VectorType, xt::xarray<bool>> transmit_with_erasures(
      const VectorType &message) {
    this->ValidateLegacyMessageSize(message);
    auto res = transmit_unsafe_with_erasures(this->ToGFArray(message));
    return {res.received.ToElementArray(), std::move(res.erasure_mask)};
  }

  std::pair<VectorType, xt::xarray<bool>> transmit_unsafe_with_erasures(
      const VectorType &message) {
    auto res = transmit_unsafe_with_erasures(this->ToGFArray(message));
    return {res.received.ToElementArray(), std::move(res.erasure_mask)};
  }

  ArrayType transmit_unsafe(const ArrayType &message) override {
    return transmit_unsafe_with_erasures(message).received;
  }

  void print(std::ostream &os) const override {
    os << "q-ary erasure channel with erasure probability " << epsilon_
       << ", of input and output space vector of dimension "
       << this->input_size_ << " over ";
    this->print_field(os);
  }

 private:
  double epsilon_;
};

using BinaryErasureChannel = QaryErasureChannel<GaloisFieldBinary>;

// BPSK Additive White Gaussian Noise (AWGN) channel:
//   BPSK mapping: 0 -> +1.0, 1 -> -1.0 (unit symbol energy E_s = 1)
//   Noise: z ~ N(0, sigma^2), where sigma^2 = 1 / (2 * R * 10^(snr_db / 10))
//   Exact LLR: ln(P(bit=0 | y) / P(bit=1 | y)) = 2 * y / sigma^2
//   Hard decision: y >= 0 -> 0, y < 0 -> 1.
template <typename BinaryField = GaloisFieldBinary>
class BasicAWGNChannel : public Channel<BinaryField> {
 public:
  using ArrayType = GFArray<BinaryField>;
  using ElementType = GaloisFieldElementBase<BinaryField>;
  using VectorType = xt::xarray<ElementType>;

  using Channel<BinaryField>::operator();
  using Channel<BinaryField>::transmit;
  using Channel<BinaryField>::Transmit;
  using Channel<BinaryField>::transmit_unsafe;

  BasicAWGNChannel(std::shared_ptr<BinaryField> field, size_t input_size,
                   double snr_db, double code_rate = 1.0,
                   std::optional<uint64_t> seed = std::nullopt)
      : Channel<BinaryField>(std::move(field), input_size, 0, seed),
        snr_db_(snr_db),
        code_rate_(code_rate),
        sigma_(ComputeSigmaFromSNR(snr_db, code_rate)) {
    if (this->field_->Order() != 2) {
      throw std::invalid_argument("AWGNChannel requires a binary field GF(2)");
    }
  }

  BasicAWGNChannel(size_t input_size, double snr_db, double code_rate = 1.0,
                   std::optional<uint64_t> seed = std::nullopt)
      : BasicAWGNChannel(std::make_shared<BinaryField>(), input_size, snr_db,
                         code_rate, seed) {}

  static BasicAWGNChannel FromSNR(std::shared_ptr<BinaryField> field,
                                  size_t input_size, double snr_db,
                                  double code_rate = 1.0,
                                  std::optional<uint64_t> seed = std::nullopt) {
    return BasicAWGNChannel(std::move(field), input_size, snr_db, code_rate,
                            seed);
  }

  static BasicAWGNChannel FromSNR(size_t input_size, double snr_db,
                                  double code_rate = 1.0,
                                  std::optional<uint64_t> seed = std::nullopt) {
    return BasicAWGNChannel(input_size, snr_db, code_rate, seed);
  }

  static BasicAWGNChannel FromEsN0(
      std::shared_ptr<BinaryField> field, size_t input_size, double es_n0_db,
      std::optional<uint64_t> seed = std::nullopt) {
    return BasicAWGNChannel(std::move(field), input_size, es_n0_db, 1.0, seed);
  }

  static BasicAWGNChannel FromEbN0(
      std::shared_ptr<BinaryField> field, size_t input_size, double eb_n0_db,
      double code_rate, std::optional<uint64_t> seed = std::nullopt) {
    return BasicAWGNChannel(std::move(field), input_size, eb_n0_db, code_rate,
                            seed);
  }

  static BasicAWGNChannel FromSigma(
      std::shared_ptr<BinaryField> field, size_t input_size, double sigma,
      std::optional<uint64_t> seed = std::nullopt) {
    return BasicAWGNChannel(SigmaTag{}, std::move(field), input_size, sigma,
                            seed);
  }

  static BasicAWGNChannel FromSigma(
      size_t input_size, double sigma,
      std::optional<uint64_t> seed = std::nullopt) {
    return BasicAWGNChannel(SigmaTag{}, std::make_shared<BinaryField>(),
                            input_size, sigma, seed);
  }

  double sigma() const noexcept { return sigma_; }
  double noise_variance() const noexcept { return sigma_ * sigma_; }
  double snr_db() const noexcept { return snr_db_; }
  double code_rate() const noexcept { return code_rate_; }

  // Exact uncoded BPSK BER over AWGN: Q(1 / sigma) = 0.5 * erfc(1 / (sigma *
  // sqrt(2))).
  double theoretical_ber() const {
    if (sigma_ == 0.0) {
      return 0.0;
    }
    return 0.5 * std::erfc(1.0 / (sigma_ * std::sqrt(2.0)));
  }

  static std::vector<double> ModulateBPSK(const ArrayType &bits) {
    if (!bits.HasField()) {
      throw std::invalid_argument("ModulateBPSK: input GFArray has null field");
    }
    const auto field = bits.Field();
    const auto span = bits.RawSpan();
    std::vector<double> symbols(span.size());
    for (size_t i = 0; i < span.size(); ++i) {
      symbols[i] = field->IsZero(span[i]) ? 1.0 : -1.0;
    }
    return symbols;
  }

  std::vector<double> ComputeLLR(
      const std::vector<double> &soft_received) const {
    std::vector<double> llr(soft_received.size());
    const double var = sigma_ * sigma_;
    if (var == 0.0) {
      for (size_t i = 0; i < soft_received.size(); ++i) {
        if (soft_received[i] > 0.0) {
          llr[i] = std::numeric_limits<double>::infinity();
        } else if (soft_received[i] < 0.0) {
          llr[i] = -std::numeric_limits<double>::infinity();
        } else {
          llr[i] = 0.0;
        }
      }
      return llr;
    }
    const double scale = 2.0 / var;
    for (size_t i = 0; i < soft_received.size(); ++i) {
      llr[i] = scale * soft_received[i];
    }
    return llr;
  }

  ArrayType HardDecision(const std::vector<double> &soft_or_llr) const {
    ArrayType hard = ArrayType::Zeros(this->field_, {soft_or_llr.size()});
    auto span = hard.RawSpan();
    const auto zero = this->field_->AdditiveIdentity();
    const auto one = this->field_->MultiplicativeIdentity();
    for (size_t i = 0; i < soft_or_llr.size(); ++i) {
      span[i] = (soft_or_llr[i] >= 0.0) ? zero : one;
    }
    return hard;
  }

  std::vector<double> transmit_soft(const ArrayType &message) {
    this->ValidateMessage(message);
    return transmit_soft_unsafe(message);
  }

  std::vector<double> transmit_llr(const ArrayType &message) {
    return ComputeLLR(transmit_soft(message));
  }

  AWGNTransmission<BinaryField> transmit_detailed(const ArrayType &message) {
    this->ValidateMessage(message);
    auto soft = transmit_soft_unsafe(message);
    auto llr = ComputeLLR(soft);
    auto hard = HardDecision(soft);
    return AWGNTransmission<BinaryField>{std::move(soft), std::move(llr),
                                         std::move(hard)};
  }

  ArrayType transmit_unsafe(const ArrayType &message) override {
    return HardDecision(transmit_soft_unsafe(message));
  }

  void print(std::ostream &os) const override {
    os << "AWGN channel (BPSK, sigma=" << sigma_ << ", snr_db=" << snr_db_
       << ", rate=" << code_rate_ << ") of dimension " << this->input_size_
       << " over ";
    this->print_field(os);
  }

 private:
  struct SigmaTag {};

  BasicAWGNChannel(SigmaTag, std::shared_ptr<BinaryField> field,
                   size_t input_size, double sigma,
                   std::optional<uint64_t> seed)
      : Channel<BinaryField>(std::move(field), input_size, 0, seed),
        snr_db_(ComputeSNRFromSigma(sigma)),
        code_rate_(1.0),
        sigma_(sigma) {
    if (this->field_->Order() != 2) {
      throw std::invalid_argument("AWGNChannel requires a binary field GF(2)");
    }
  }

  static double ComputeSigmaFromSNR(double snr_db, double code_rate) {
    if (std::isnan(snr_db) || std::isinf(snr_db)) {
      throw std::invalid_argument("AWGNChannel: snr_db must be finite");
    }
    if (std::isnan(code_rate) || std::isinf(code_rate) || code_rate <= 0.0 ||
        code_rate > 1.0) {
      throw std::invalid_argument("AWGNChannel: code_rate must be in (0, 1]");
    }
    const double es_n0_linear = code_rate * std::pow(10.0, snr_db / 10.0);
    return std::sqrt(1.0 / (2.0 * es_n0_linear));
  }

  static double ComputeSNRFromSigma(double sigma) {
    if (std::isnan(sigma) || std::isinf(sigma) || sigma < 0.0) {
      throw std::invalid_argument(
          "AWGNChannel: sigma must be finite and non-negative");
    }
    if (sigma == 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    return 10.0 * std::log10(1.0 / (2.0 * sigma * sigma));
  }

  std::vector<double> transmit_soft_unsafe(const ArrayType &message) {
    auto soft = ModulateBPSK(message);
    if (sigma_ == 0.0) {
      return soft;
    }
    std::normal_distribution<double> noise_dist(0.0, sigma_);
    for (double &y : soft) {
      y += noise_dist(this->rng_);
    }
    return soft;
  }

  double snr_db_;
  double code_rate_;
  double sigma_;
};

using AWGNChannel = BasicAWGNChannel<GaloisFieldBinary>;

// Configuration for deterministic Monte Carlo BER/SER/FER simulations.
struct SimulationConfig {
  size_t max_frames = 10000;
  size_t target_frame_errors = 100;  // 0 disables early stopping.
  uint64_t seed = 42;
};

// Aggregated Monte Carlo error-rate statistics.
struct SimulationStats {
  size_t frames_sent = 0;
  size_t frame_errors = 0;
  size_t symbols_sent = 0;
  size_t symbol_errors = 0;
  size_t bits_sent = 0;
  size_t bit_errors = 0;
  size_t decoding_failures = 0;

  double FER() const noexcept {
    return frames_sent == 0 ? 0.0
                            : static_cast<double>(frame_errors) /
                                  static_cast<double>(frames_sent);
  }

  double SER() const noexcept {
    return symbols_sent == 0 ? 0.0
                             : static_cast<double>(symbol_errors) /
                                   static_cast<double>(symbols_sent);
  }

  double BER() const noexcept {
    return bits_sent == 0 ? 0.0
                          : static_cast<double>(bit_errors) /
                                static_cast<double>(bits_sent);
  }
};

namespace detail {

template <typename GaloisField>
void TallyMessageErrors(const GaloisField &field,
                        const GFArray<GaloisField> &expected,
                        const GFArray<GaloisField> &actual,
                        size_t bits_per_symbol, SimulationStats &stats) {
  const size_t k = expected.Size();
  if (!actual.HasField() || actual.Size() != k) {
    ++stats.frame_errors;
    stats.symbol_errors += k;
    stats.bit_errors += k * bits_per_symbol;
    return;
  }
  const auto exp_span = expected.RawSpan();
  const auto act_span = actual.RawSpan();
  size_t sym_err = 0;
  size_t bit_err = 0;
  for (size_t i = 0; i < k; ++i) {
    if (exp_span[i] != act_span[i]) {
      ++sym_err;
      bit_err += CountSymbolBitErrors(field, exp_span[i], act_span[i],
                                      bits_per_symbol);
    }
  }
  if (sym_err > 0) {
    ++stats.frame_errors;
    stats.symbol_errors += sym_err;
    stats.bit_errors += bit_err;
  }
}

}  // namespace detail

// Runs deterministic seeded Monte Carlo trials for any (encode_fn, channel_fn,
// decode_fn) pipeline. `decode_fn` may return:
//   - a DecodingResult-like struct with `.Ok()` and `.message`,
//   - a `std::optional<GFArray<GaloisField>>`, or
//   - a `GFArray<GaloisField>` directly.
template <typename GaloisField, typename EncodeFn, typename ChannelFn,
          typename DecodeFn>
SimulationStats SimulateCodePerformance(
    const std::shared_ptr<GaloisField> &field, size_t message_length,
    EncodeFn &&encode_fn, ChannelFn &&channel_fn, DecodeFn &&decode_fn,
    SimulationConfig config = {}) {
  if (!field) {
    throw std::invalid_argument(
        "SimulateCodePerformance: field pointer cannot be null");
  }

  SimulationStats stats;
  if (config.max_frames == 0) {
    return stats;
  }

  const size_t bits_per_symbol = detail::BitsPerSymbol(*field);
  std::mt19937_64 msg_rng(config.seed);

  for (size_t frame = 0; frame < config.max_frames; ++frame) {
    const GFArray<GaloisField> message =
        GFArray<GaloisField>::Random(field, {message_length}, msg_rng);
    auto codeword = encode_fn(message);
    auto received = channel_fn(codeword);
    auto decoded = decode_fn(received);
    using DecodedType = std::remove_cvref_t<decltype(decoded)>;

    ++stats.frames_sent;
    stats.symbols_sent += message_length;
    stats.bits_sent += message_length * bits_per_symbol;

    if constexpr (requires(const DecodedType &d) {
                    { d.Ok() } -> std::convertible_to<bool>;
                    d.message;
                  }) {
      if (!decoded.Ok()) {
        ++stats.decoding_failures;
        ++stats.frame_errors;
        stats.symbol_errors += message_length;
        stats.bit_errors += message_length * bits_per_symbol;
      } else {
        detail::TallyMessageErrors(*field, message, decoded.message,
                                   bits_per_symbol, stats);
      }
    } else if constexpr (requires(const DecodedType &d) {
                           { d.has_value() } -> std::convertible_to<bool>;
                           *d;
                         }) {
      if (!decoded.has_value()) {
        ++stats.decoding_failures;
        ++stats.frame_errors;
        stats.symbol_errors += message_length;
        stats.bit_errors += message_length * bits_per_symbol;
      } else {
        detail::TallyMessageErrors(*field, message, *decoded, bits_per_symbol,
                                   stats);
      }
    } else {
      detail::TallyMessageErrors(*field, message, decoded, bits_per_symbol,
                                 stats);
    }

    if (config.target_frame_errors > 0 &&
        stats.frame_errors >= config.target_frame_errors) {
      break;
    }
  }

  return stats;
}

template <typename GaloisField, typename EncodeFn, typename ChannelFn,
          typename DecodeFn>
SimulationStats SimulateCodePerformance(
    const std::shared_ptr<GaloisField> &field, size_t message_length,
    EncodeFn &&encode_fn, ChannelFn &&channel_fn, DecodeFn &&decode_fn,
    size_t max_frames, size_t target_frame_errors, uint64_t seed = 42) {
  return SimulateCodePerformance(
      field, message_length, std::forward<EncodeFn>(encode_fn),
      std::forward<ChannelFn>(channel_fn), std::forward<DecodeFn>(decode_fn),
      SimulationConfig{max_frames, target_frame_errors, seed});
}

// Convenience overload integrating any AbstractCode-like object and Channel.
template <typename Code, typename ChannelType>
SimulationStats SimulateCodePerformance(const Code &code, ChannelType &channel,
                                        SimulationConfig config = {}) {
  using FieldType = typename Code::field_type;
  const auto field = code.Field();
  const size_t k = code.Dimension();

  // Derive a deterministic, independent channel seed from config.seed.
  channel.SetSeed(config.seed ^ 0x9e3779b97f4a7c15ULL);

  auto encode_fn = [&code](const GFArray<FieldType> &msg) {
    return code.Encode(msg);
  };

  if constexpr (requires(ChannelType &ch, const GFArray<FieldType> &cw,
                         const Code &c, const std::vector<size_t> &pos) {
                  ch.transmit_with_erasures(cw);
                  c.DecodeWithErasures(cw, pos);
                }) {
    auto channel_fn = [&channel](const GFArray<FieldType> &cw) {
      return channel.transmit_with_erasures(cw);
    };
    auto decode_fn = [&code, k](const ErasureTransmission<FieldType> &tx) {
      auto res = code.DecodeWithErasures(tx.received, tx.erasure_positions);
      if constexpr (requires(const Code &c, const GFArray<FieldType> &cw) {
                      c.Unencode(cw);
                    }) {
        if (res.Ok() && (!res.message.HasField() || res.message.Size() != k)) {
          res.message = code.Unencode(res.codeword);
        }
      }
      return res;
    };
    return SimulateCodePerformance(field, k, encode_fn, channel_fn, decode_fn,
                                   config);
  } else {
    auto channel_fn = [&channel](const GFArray<FieldType> &cw) {
      return channel.transmit(cw);
    };
    auto decode_fn = [&code, k](const GFArray<FieldType> &rx) {
      auto res = code.DecodeWithResult(rx);
      if constexpr (requires(const Code &c, const GFArray<FieldType> &cw) {
                      c.Unencode(cw);
                    }) {
        if (res.Ok() && (!res.message.HasField() || res.message.Size() != k)) {
          res.message = code.Unencode(res.codeword);
        }
      }
      return res;
    };
    return SimulateCodePerformance(field, k, encode_fn, channel_fn, decode_fn,
                                   config);
  }
}

template <typename GaloisField>
std::ostream &operator<<(std::ostream &os,
                         const Channel<GaloisField> &channel) {
  channel.print(os);
  return os;
}

template <typename GaloisField>
std::shared_ptr<StaticErrorRateChannel<GaloisField>>
CreateStaticErrorRateChannel(const std::shared_ptr<GaloisField> &field,
                             size_t input_size, size_t num_errors,
                             std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<StaticErrorRateChannel<GaloisField>>(
      field, input_size, num_errors, seed);
}

template <typename GaloisField>
std::shared_ptr<StaticErrorRateChannel<GaloisField>>
CreateStaticErrorRateChannel(const std::shared_ptr<GaloisField> &field,
                             size_t input_size, size_t min_errors,
                             size_t max_errors,
                             std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<StaticErrorRateChannel<GaloisField>>(
      field, input_size, min_errors, max_errors, seed);
}

template <typename GaloisField>
std::shared_ptr<ErrorErasureChannel<GaloisField>> CreateErrorErasureChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    size_t num_errors, size_t num_erasures,
    std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<ErrorErasureChannel<GaloisField>>(
      field, input_size, num_errors, num_erasures, seed);
}

template <typename GaloisField>
std::shared_ptr<ErrorErasureChannel<GaloisField>> CreateErrorErasureChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    size_t min_errors, size_t max_errors, size_t min_erasures,
    size_t max_erasures, std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<ErrorErasureChannel<GaloisField>>(
      field, input_size, min_errors, max_errors, min_erasures, max_erasures,
      seed);
}

template <typename GaloisField>
std::shared_ptr<QarySymmetricChannel<GaloisField>> CreateQarySymmetricChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    double epsilon, std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<QarySymmetricChannel<GaloisField>>(field, input_size,
                                                             epsilon, seed);
}

template <typename BinaryField = GaloisFieldBinary>
std::shared_ptr<BasicBinarySymmetricChannel<BinaryField>>
CreateBinarySymmetricChannel(const std::shared_ptr<BinaryField> &field,
                             size_t input_size, double crossover_prob,
                             std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<BasicBinarySymmetricChannel<BinaryField>>(
      field, input_size, crossover_prob, seed);
}

template <typename BinaryField = GaloisFieldBinary>
std::shared_ptr<BasicBinarySymmetricChannel<BinaryField>>
CreateBinarySymmetricChannel(size_t input_size, double crossover_prob,
                             std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<BasicBinarySymmetricChannel<BinaryField>>(
      input_size, crossover_prob, seed);
}

template <typename GaloisField>
std::shared_ptr<QaryErasureChannel<GaloisField>> CreateQaryErasureChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    double epsilon, std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<QaryErasureChannel<GaloisField>>(field, input_size,
                                                           epsilon, seed);
}

template <typename BinaryField = GaloisFieldBinary>
std::shared_ptr<BasicAWGNChannel<BinaryField>> CreateAWGNChannel(
    const std::shared_ptr<BinaryField> &field, size_t input_size, double snr_db,
    double code_rate = 1.0, std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<BasicAWGNChannel<BinaryField>>(
      field, input_size, snr_db, code_rate, seed);
}

template <typename BinaryField = GaloisFieldBinary>
std::shared_ptr<BasicAWGNChannel<BinaryField>> CreateAWGNChannel(
    size_t input_size, double snr_db, double code_rate = 1.0,
    std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<BasicAWGNChannel<BinaryField>>(input_size, snr_db,
                                                         code_rate, seed);
}

template <typename BinaryField = GaloisFieldBinary>
std::shared_ptr<BasicAWGNChannel<BinaryField>> CreateAWGNChannelFromSigma(
    const std::shared_ptr<BinaryField> &field, size_t input_size, double sigma,
    std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<BasicAWGNChannel<BinaryField>>(
      BasicAWGNChannel<BinaryField>::FromSigma(field, input_size, sigma, seed));
}

template <typename BinaryField = GaloisFieldBinary>
std::shared_ptr<BasicAWGNChannel<BinaryField>> CreateAWGNChannelFromSigma(
    size_t input_size, double sigma,
    std::optional<uint64_t> seed = std::nullopt) {
  return std::make_shared<BasicAWGNChannel<BinaryField>>(
      BasicAWGNChannel<BinaryField>::FromSigma(input_size, sigma, seed));
}

}  // namespace channel

// Backwards-compatibility alias for existing `xg::channels::*` callers.
namespace channels = channel;

}  // namespace xg

#endif  // XGALOIS_CHANNEL_CHANNEL_HPP_
