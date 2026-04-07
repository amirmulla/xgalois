#ifndef XGALOIS_CODING_GENERATOR_MATRIX_ENCODER_HPP
#define XGALOIS_CODING_GENERATOR_MATRIX_ENCODER_HPP

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class GeneratorMatrixEncoder : public Encoder<GaloisField> {
 public:
  using element_type = typename Encoder<GaloisField>::element_type;
  using codeword_type = typename Encoder<GaloisField>::codeword_type;
  using message_type = typename Encoder<GaloisField>::message_type;
  using matrix_type = xg::garray<GaloisField>;
  using vector_type = xg::garray<GaloisField>;

  explicit GeneratorMatrixEncoder(const AbstractCode<GaloisField>* code)
      : Encoder<GaloisField>(code) {

    linear_code_ = dynamic_cast<const AbstractLinearCode<GaloisField>*>(code);
    if (!linear_code_) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder requires a linear code");
    }
  }

  codeword_type Encode(const message_type& message) const override {
    if (message.size() != MessageLength()) {
      throw std::invalid_argument("Message length must match code dimension");
    }

    auto generator = linear_code_->GeneratorMatrix();
    codeword_type codeword_vec = xg::linalg::dot(message, generator);

    return codeword_vec;
  }

  message_type Unencode(const codeword_type& codeword) const override {
    if (codeword.size() != this->code_->Length()) {
      throw std::invalid_argument("Codeword length must match code length");
    }

    if (!this->code_->Contains(codeword)) {
      throw std::invalid_argument("Input is not a valid codeword");
    }

    auto generator = linear_code_->GeneratorMatrix();

    if (IsSystematic(generator)) {
      message_type message = xg::linalg::zeros<GaloisField>({MessageLength()},
                                                            this->code_->Field());
      for (size_t i = 0; i < MessageLength(); ++i) {
        message(i) = codeword(i);
      }
      return message;
    } else {

      return SolveForMessage(codeword, generator);
    }
  }

  size_t MessageLength() const override { return linear_code_->Dimension(); }

  std::string ToString() const override {
    return "Generator Matrix Encoder for [" +
           std::to_string(this->code_->Length()) + ", " +
           std::to_string(linear_code_->Dimension()) + "] linear code";
  }

 private:
  const AbstractLinearCode<GaloisField>* linear_code_;

  codeword_type ToCodeword(const vector_type& vec) const {
    return vec;
  }

  bool IsSystematic(const matrix_type& generator) const {

    size_t k = generator.shape(0);

    for (size_t i = 0; i < k; ++i) {
      for (size_t j = 0; j < k; ++j) {
        if (i == j) {
          if (generator(i, j) != element_type(1, this->code_->Field())) {
            return false;
          }
        } else {
          if (generator(i, j) != element_type(0, this->code_->Field())) {
            return false;
          }
        }
      }
    }
    return true;
  }

  message_type SolveForMessage(const codeword_type& codeword,
                               const matrix_type& generator) const {

    auto field = this->code_->Field();
    size_t q = field->Order();
    size_t k = MessageLength();

    for (size_t i = 0; i < static_cast<size_t>(std::pow(q, k)); ++i) {
      message_type candidate = xg::linalg::zeros<GaloisField>({k}, field);
      size_t temp = i;

      for (size_t j = 0; j < k; ++j) {
        candidate(j) = element_type(temp % q, field);
        temp /= q;
      }

      if (Encode(candidate) == codeword) {
        return candidate;
      }
    }

    throw std::runtime_error("Failed to find message for codeword");
  }
};

}
}

#endif
