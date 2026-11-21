#ifndef XGALOIS_CODING_ENCODER_CYCLIC_ENCODER_HPP_
#define XGALOIS_CODING_ENCODER_CYCLIC_ENCODER_HPP_

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class CyclicEncoder : public Encoder<GaloisField> {
 public:
  using Base = Encoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using polynomial_type = PolynomialDense<GaloisField>;

  CyclicEncoder(const AbstractLinearCode<GaloisField> &code,
                polynomial_type generator_poly, bool systematic = false)
      : Base(code),
        generator_poly_(std::move(generator_poly)),
        n_(code.Length()),
        k_(code.Dimension()),
        systematic_(systematic) {
    if (!generator_poly_.HasField() ||
        generator_poly_.Field() != this->code_->Field()) {
      throw std::invalid_argument(
          "CyclicEncoder: generator polynomial field mismatch");
    }
    if (generator_poly_.IsZero() ||
        static_cast<std::size_t>(generator_poly_.Degree()) != n_ - k_) {
      throw std::invalid_argument(
          "CyclicEncoder: generator polynomial degree must equal n - k");
    }
  }

  template <typename CyclicCodeType>
    requires requires(const CyclicCodeType &c) { c.GeneratorPolynomial(); }
  explicit CyclicEncoder(const CyclicCodeType &code, bool systematic = false)
      : CyclicEncoder(code, code.GeneratorPolynomial(), systematic) {}

  codeword_type Encode(const message_type &message) const override {
    auto field = this->code_->Field();
    if (!message.HasField() || message.Field() != field) {
      throw std::invalid_argument(
          "CyclicEncoder::Encode: message field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != k_) {
      throw std::invalid_argument(
          "CyclicEncoder::Encode: message must be 1D with length equal to "
          "code dimension");
    }
    if (k_ == 0) {
      return codeword_type::Zeros(field, {n_});
    }

    std::vector<storage_type> m_raw(message.RawSpan().begin(),
                                    message.RawSpan().end());
    polynomial_type m_poly = polynomial_type::FromRaw(std::move(m_raw), field);

    polynomial_type c_poly;
    if (!systematic_) {
      // Non-systematic: c(x) = m(x) * g(x)
      c_poly = m_poly * generator_poly_;
    } else {
      // Systematic: c(x) = x^(n-k) * m(x) - (x^(n-k) * m(x) mod g(x))
      const std::size_t r = n_ - k_;
      polynomial_type shifted = m_poly.Shift(r);
      polynomial_type parity = shifted % generator_poly_;
      c_poly = shifted - parity;
    }

    return PolyToCodeword(c_poly);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    auto field = this->code_->Field();
    if (!codeword.HasField() || codeword.Field() != field) {
      throw std::invalid_argument(
          "CyclicEncoder::Unencode: codeword field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != n_) {
      throw std::invalid_argument(
          "CyclicEncoder::Unencode: codeword must be 1D with length equal to "
          "code length");
    }

    std::vector<storage_type> c_raw(codeword.RawSpan().begin(),
                                    codeword.RawSpan().end());
    polynomial_type c_poly = polynomial_type::FromRaw(std::move(c_raw), field);

    if (!systematic_) {
      auto [q, rem] = c_poly.DivRem(generator_poly_);
      if (!rem.IsZero()) {
        throw std::invalid_argument(
            "CyclicEncoder::Unencode: input is not a valid cyclic codeword");
      }
      if (k_ == 0) {
        return message_type::Zeros(field, {0});
      }
      return PolyToMessage(q);
    }

    // Systematic: verify c(x) % g(x) == 0 and extract top k coefficients.
    polynomial_type rem = c_poly % generator_poly_;
    if (!rem.IsZero()) {
      throw std::invalid_argument(
          "CyclicEncoder::Unencode: input is not a valid cyclic codeword");
    }
    if (k_ == 0) {
      return message_type::Zeros(field, {0});
    }

    using RawArray = typename message_type::raw_array_type;
    using ShapeType = typename message_type::shape_type;
    RawArray out = RawArray::from_shape(ShapeType{k_});
    const std::size_t r = n_ - k_;
    const storage_type *src = codeword.RawData().data() + r;
    std::copy(src, src + k_, out.data());
    return message_type::FromRaw(std::move(field), std::move(out));
  }

  std::string Name() const override {
    return systematic_ ? "CyclicSystematic" : "Cyclic";
  }

  bool IsSystematic() const noexcept { return systematic_; }

  std::size_t MessageLength() const noexcept { return k_; }

  std::string ToString() const override {
    return Name() + " encoder for [" + std::to_string(n_) + ", " +
           std::to_string(k_) + "] cyclic code";
  }

 private:
  codeword_type PolyToCodeword(const polynomial_type &poly) const {
    auto field = this->code_->Field();
    codeword_type out = codeword_type::Zeros(field, {n_});
    if (!poly.IsZero()) {
      const std::size_t limit = std::min<std::size_t>(
          n_, static_cast<std::size_t>(poly.Degree()) + 1);
      storage_type *dst = out.RawData().data();
      for (std::size_t i = 0; i < limit; ++i) {
        dst[i] = poly.RawCoeff(i);
      }
    }
    return out;
  }

  message_type PolyToMessage(const polynomial_type &poly) const {
    auto field = this->code_->Field();
    message_type out = message_type::Zeros(field, {k_});
    if (!poly.IsZero()) {
      const std::size_t limit = std::min<std::size_t>(
          k_, static_cast<std::size_t>(poly.Degree()) + 1);
      storage_type *dst = out.RawData().data();
      for (std::size_t i = 0; i < limit; ++i) {
        dst[i] = poly.RawCoeff(i);
      }
    }
    return out;
  }

  polynomial_type generator_poly_;
  std::size_t n_;
  std::size_t k_;
  bool systematic_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ENCODER_CYCLIC_ENCODER_HPP_
