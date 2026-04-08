#ifndef XGALOIS_CHANNEL_CHANNEL_HPP
#define XGALOIS_CHANNEL_CHANNEL_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/field/gf_element.hpp"
#include "xtensor/containers/xarray.hpp"

namespace xg {
namespace channels {

template <typename GaloisField>
class Channel {
 public:
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  Channel(std::shared_ptr<GaloisField> field, size_t input_size,
          size_t output_size = 0)
      : field_(std::move(field)),
        input_size_(input_size),
        output_size_(output_size == 0 ? input_size : output_size) {}

  virtual ~Channel() = default;

  std::shared_ptr<GaloisField> field() const { return field_; }

  size_t input_size() const { return input_size_; }

  size_t output_size() const { return output_size_; }

  VectorType transmit(const VectorType &message) {
    if (message.size() != input_size_) {
      throw std::invalid_argument("Message size (" +
                                  std::to_string(message.size()) +
                                  ") does not match channel input size (" +
                                  std::to_string(input_size_) + ")");
    }
    return transmit_unsafe(message);
  }

  virtual VectorType transmit_unsafe(const VectorType &message) = 0;

  VectorType operator()(const VectorType &message) { return transmit(message); }

  virtual void print(std::ostream &os) const = 0;

  virtual std::string to_string() const {
    std::ostringstream oss;
    print(oss);
    return oss.str();
  }

 protected:
  std::shared_ptr<GaloisField> field_;
  size_t input_size_;
  size_t output_size_;

  VectorType generate_error_vector(const std::vector<size_t> &error_positions) {
    VectorType error_vector =
        xt::xarray<ElementType>::from_shape({input_size_});

    ElementType zero_element(field_->AdditiveIdentity(), field_);
    for (size_t i = 0; i < input_size_; ++i) {
      error_vector[i] = zero_element;
    }

    for (size_t pos : error_positions) {
      if (pos < input_size_) {
        ElementType error_value;
        do {
          error_value = ElementType(field_->Random(), field_);
        } while (error_value == zero_element);
        error_vector[pos] = error_value;
      }
    }

    return error_vector;
  }

  std::vector<size_t> generate_error_positions(size_t num_errors) {
    if (num_errors > input_size_) {
      throw std::invalid_argument("Number of errors cannot exceed input size");
    }

    std::vector<size_t> all_positions(input_size_);
    std::iota(all_positions.begin(), all_positions.end(), 0);

    static thread_local std::mt19937 gen(std::random_device{}());
    std::shuffle(all_positions.begin(), all_positions.end(), gen);

    std::vector<size_t> error_positions(all_positions.begin(),
                                        all_positions.begin() + num_errors);
    return error_positions;
  }

  std::string format_interval(size_t min_val, size_t max_val) const {
    if (min_val == max_val) {
      return std::to_string(min_val);
    } else {
      return "between " + std::to_string(min_val) + " and " +
             std::to_string(max_val);
    }
  }
};

template <typename GaloisField>
class StaticErrorRateChannel : public Channel<GaloisField> {
 public:
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  StaticErrorRateChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                         size_t num_errors)
      : Channel<GaloisField>(std::move(field), input_size),
        min_errors_(num_errors),
        max_errors_(num_errors) {}

  StaticErrorRateChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                         size_t min_errors, size_t max_errors)
      : Channel<GaloisField>(std::move(field), input_size),
        min_errors_(min_errors),
        max_errors_(max_errors) {
    if (min_errors > max_errors) {
      throw std::invalid_argument(
          "min_errors cannot be greater than max_errors");
    }
  }

  std::pair<size_t, size_t> number_errors() const {
    return {min_errors_, max_errors_};
  }

  VectorType transmit_unsafe(const VectorType &message) override {
    size_t num_errors;
    if (min_errors_ == max_errors_) {
      num_errors = min_errors_;
    } else {
      static thread_local std::mt19937 gen(std::random_device{}());
      std::uniform_int_distribution<size_t> dist(min_errors_, max_errors_);
      num_errors = dist(gen);
    }

    if (num_errors == 0) {
      return message;
    }

    auto error_positions = this->generate_error_positions(num_errors);
    auto error_vector = this->generate_error_vector(error_positions);

    VectorType result = message;
    for (size_t i = 0; i < message.size(); ++i) {
      result[i] = message[i] + error_vector[i];
    }

    return result;
  }

  void print(std::ostream &os) const override {
    os << "Static error rate channel creating "
       << this->format_interval(min_errors_, max_errors_)
       << " errors, of input and output space vector of dimension "
       << this->input_size_ << " over ";
    this->field_->Print(os);
  }

 private:
  size_t min_errors_;
  size_t max_errors_;
};

template <typename GaloisField>
class ErrorErasureChannel : public Channel<GaloisField> {
 public:
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  ErrorErasureChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                      size_t num_errors, size_t num_erasures)
      : Channel<GaloisField>(std::move(field), input_size),
        min_errors_(num_errors),
        max_errors_(num_errors),
        min_erasures_(num_erasures),
        max_erasures_(num_erasures) {}

  ErrorErasureChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                      size_t min_errors, size_t max_errors, size_t min_erasures,
                      size_t max_erasures)
      : Channel<GaloisField>(std::move(field), input_size),
        min_errors_(min_errors),
        max_errors_(max_errors),
        min_erasures_(min_erasures),
        max_erasures_(max_erasures) {
    if (min_errors > max_errors || min_erasures > max_erasures) {
      throw std::invalid_argument(
          "min values cannot be greater than max values");
    }
  }

  std::pair<size_t, size_t> number_errors() const {
    return {min_errors_, max_errors_};
  }

  std::pair<size_t, size_t> number_erasures() const {
    return {min_erasures_, max_erasures_};
  }

  std::pair<VectorType, xt::xarray<bool>> transmit_unsafe_with_erasures(
      const VectorType &message) {
    size_t num_errors =
        (min_errors_ == max_errors_)
            ? min_errors_
            : std::uniform_int_distribution<size_t>(min_errors_, max_errors_)(
                  *reinterpret_cast<std::mt19937 *>(&thread_local_gen()));

    size_t num_erasures =
        (min_erasures_ == max_erasures_)
            ? min_erasures_
            : std::uniform_int_distribution<size_t>(min_erasures_,
                                                    max_erasures_)(
                  *reinterpret_cast<std::mt19937 *>(&thread_local_gen()));

    std::vector<size_t> all_positions(this->input_size_);
    std::iota(all_positions.begin(), all_positions.end(), 0);
    std::shuffle(all_positions.begin(), all_positions.end(),
                 *reinterpret_cast<std::mt19937 *>(&thread_local_gen()));

    std::vector<size_t> error_positions(
        all_positions.begin(),
        all_positions.begin() + std::min(num_errors, this->input_size_));
    std::vector<size_t> erasure_positions(
        all_positions.begin() + num_errors,
        all_positions.begin() +
            std::min(num_errors + num_erasures, this->input_size_));

    VectorType result = message;
    xt::xarray<bool> erasure_vector =
        xt::xarray<bool>::from_shape({this->input_size_});
    std::fill(erasure_vector.begin(), erasure_vector.end(), false);

    auto error_vector = this->generate_error_vector(error_positions);
    for (size_t i = 0; i < message.size(); ++i) {
      result[i] = message[i] + error_vector[i];
    }

    for (size_t pos : erasure_positions) {
      if (pos < this->input_size_) {
        result[pos] =
            ElementType(this->field_->AdditiveIdentity(), this->field_);
        erasure_vector[pos] = true;
      }
    }

    return {result, erasure_vector};
  }

  VectorType transmit_unsafe(const VectorType &message) override {
    return transmit_unsafe_with_erasures(message).first;
  }

  void print(std::ostream &os) const override {
    os << "Error-and-erasure channel creating "
       << this->format_interval(min_errors_, max_errors_) << " errors and "
       << this->format_interval(min_erasures_, max_erasures_)
       << " erasures of input space vector of dimension " << this->input_size_
       << " over ";
    this->field_->Print(os);
  }

 private:
  size_t min_errors_;
  size_t max_errors_;
  size_t min_erasures_;
  size_t max_erasures_;

  static std::mt19937 &thread_local_gen() {
    static thread_local std::mt19937 gen(std::random_device{}());
    return gen;
  }
};

template <typename GaloisField>
class QarySymmetricChannel : public Channel<GaloisField> {
 public:
  using ElementType = GaloisFieldElementBase<GaloisField>;
  using VectorType = xt::xarray<ElementType>;

  QarySymmetricChannel(std::shared_ptr<GaloisField> field, size_t input_size,
                       double epsilon)
      : Channel<GaloisField>(std::move(field), input_size), epsilon_(epsilon) {
    if (epsilon < 0.0 || epsilon > 1.0) {
      throw std::invalid_argument("Error probability must be between 0 and 1");
    }
  }

  double error_probability() const { return epsilon_; }

  double probability_of_exactly_t_errors(size_t t) const {
    if (t > this->input_size_) return 0.0;

    double binomial_coeff = 1.0;
    for (size_t i = 0; i < t; ++i) {
      binomial_coeff *= (this->input_size_ - i) / (i + 1.0);
    }

    return binomial_coeff * std::pow(epsilon_, t) *
           std::pow(1.0 - epsilon_, this->input_size_ - t);
  }

  double probability_of_at_most_t_errors(size_t t) const {
    double prob = 0.0;
    for (size_t i = 0; i <= std::min(t, this->input_size_); ++i) {
      prob += probability_of_exactly_t_errors(i);
    }
    return prob;
  }

  VectorType transmit_unsafe(const VectorType &message) override {
    VectorType result = message;
    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (size_t i = 0; i < message.size(); ++i) {
      if (dist(gen) < epsilon_) {
        ElementType new_symbol;
        do {
          new_symbol = ElementType(this->field_->Random(), this->field_);
        } while (new_symbol == message[i]);
        result[i] = new_symbol;
      }
    }

    return result;
  }

  void print(std::ostream &os) const override {
    os << "q-ary symmetric channel with error probability " << epsilon_
       << ", of input and output space vector of dimension "
       << this->input_size_ << " over ";
    this->field_->Print(os);
  }

 private:
  double epsilon_;
};

template <typename GaloisField>
std::ostream &operator<<(std::ostream &os,
                         const Channel<GaloisField> &channel) {
  channel.print(os);
  return os;
}

template <typename GaloisField>
std::shared_ptr<StaticErrorRateChannel<GaloisField>>
CreateStaticErrorRateChannel(const std::shared_ptr<GaloisField> &field,
                             size_t input_size, size_t num_errors) {
  return std::make_shared<StaticErrorRateChannel<GaloisField>>(
      field, input_size, num_errors);
}

template <typename GaloisField>
std::shared_ptr<StaticErrorRateChannel<GaloisField>>
CreateStaticErrorRateChannel(const std::shared_ptr<GaloisField> &field,
                             size_t input_size, size_t min_errors,
                             size_t max_errors) {
  return std::make_shared<StaticErrorRateChannel<GaloisField>>(
      field, input_size, min_errors, max_errors);
}

template <typename GaloisField>
std::shared_ptr<ErrorErasureChannel<GaloisField>> CreateErrorErasureChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    size_t num_errors, size_t num_erasures) {
  return std::make_shared<ErrorErasureChannel<GaloisField>>(
      field, input_size, num_errors, num_erasures);
}

template <typename GaloisField>
std::shared_ptr<ErrorErasureChannel<GaloisField>> CreateErrorErasureChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    size_t min_errors, size_t max_errors, size_t min_erasures,
    size_t max_erasures) {
  return std::make_shared<ErrorErasureChannel<GaloisField>>(
      field, input_size, min_errors, max_errors, min_erasures, max_erasures);
}

template <typename GaloisField>
std::shared_ptr<QarySymmetricChannel<GaloisField>> CreateQarySymmetricChannel(
    const std::shared_ptr<GaloisField> &field, size_t input_size,
    double epsilon) {
  return std::make_shared<QarySymmetricChannel<GaloisField>>(field, input_size,
                                                             epsilon);
}

}  // namespace channels
}  // namespace xg

#endif
