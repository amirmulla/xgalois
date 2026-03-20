#ifndef XGALOIS_CODING_GENERATOR_MATRIX_ENCODER_HPP
#define XGALOIS_CODING_GENERATOR_MATRIX_ENCODER_HPP

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename ElementType>
class GeneratorMatrixEncoder : public Encoder<ElementType> {
 public:
  using element_type = ElementType;
  using codeword_type = std::vector<ElementType>;
  using message_type = std::vector<ElementType>;
  using matrix_type = xg::linalg::Matrix<ElementType>;
  using vector_type = xg::linalg::Vector<ElementType>;

  // Constructor
  explicit GeneratorMatrixEncoder(const AbstractCode<ElementType>* code)
      : Encoder<ElementType>(code) {
    // Try to cast to linear code
    linear_code_ = dynamic_cast<const AbstractLinearCode<ElementType>*>(code);
    if (!linear_code_) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder requires a linear code");
    }
  }

  // Encode using generator matrix: c = m * G
  codeword_type Encode(const message_type& message) const override {
    if (message.size() != MessageLength()) {
      throw std::invalid_argument("Message length must match code dimension");
    }

    auto generator = linear_code_->GeneratorMatrix();
    vector_type msg_vec(message);
    vector_type codeword_vec = msg_vec * generator;

    return ToCodeword(codeword_vec);
  }

  // Unencode: find message from codeword
  message_type Unencode(const codeword_type& codeword) const override {
    if (codeword.size() != this->code_->Length()) {
      throw std::invalid_argument("Codeword length must match code length");
    }

    if (!this->code_->Contains(codeword)) {
      throw std::invalid_argument("Input is not a valid codeword");
    }

    // For systematic codes, the message is the first k positions
    // For general codes, we need to solve the system m * G = c
    auto generator = linear_code_->GeneratorMatrix();

    // Check if generator is in systematic form [I_k | P]
    if (IsSystematic(generator)) {
      message_type message(MessageLength());
      for (size_t i = 0; i < MessageLength(); ++i) {
        message[i] = codeword[i];
      }
      return message;
    } else {
      // Solve the linear system
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
  const AbstractLinearCode<ElementType>* linear_code_;

  // Helper functions
  codeword_type ToCodeword(const vector_type& vec) const {
    codeword_type result(vec.Size());
    for (size_t i = 0; i < vec.Size(); ++i) {
      result[i] = vec[i];
    }
    return result;
  }

  bool IsSystematic(const matrix_type& generator) const {
    // Check if the first k columns form an identity matrix
    size_t k = generator.Rows();

    for (size_t i = 0; i < k; ++i) {
      for (size_t j = 0; j < k; ++j) {
        if (i == j) {
          if (generator.Get(i, j) != ElementType(1)) {
            return false;
          }
        } else {
          if (generator.Get(i, j) != ElementType{}) {
            return false;
          }
        }
      }
    }
    return true;
  }

  message_type SolveForMessage(const codeword_type& codeword,
                               const matrix_type& generator) const {
    // This is a simplified implementation
    // In practice, you'd use proper linear algebra to solve m * G = c

    // For now, use brute force search (only feasible for small codes)
    auto field = this->code_->Field();
    size_t q = field->Order();
    size_t k = MessageLength();

    // Try all possible messages
    for (size_t i = 0; i < static_cast<size_t>(std::pow(q, k)); ++i) {
      message_type candidate(k);
      size_t temp = i;

      // Convert i to base-q representation
      for (size_t j = 0; j < k; ++j) {
        candidate[j] = ElementType(temp % q);
        temp /= q;
      }

      // Check if this message encodes to our codeword
      if (Encode(candidate) == codeword) {
        return candidate;
      }
    }

    throw std::runtime_error("Failed to find message for codeword");
  }
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_GENERATOR_MATRIX_ENCODER_HPP
