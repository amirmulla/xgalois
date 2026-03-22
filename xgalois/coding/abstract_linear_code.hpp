#ifndef XGALOIS_CODING_ABSTRACT_LINEAR_CODE_HPP
#define XGALOIS_CODING_ABSTRACT_LINEAR_CODE_HPP

#include <memory>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class AbstractLinearCode : public AbstractCode<GaloisField> {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xg::garray<GaloisField>;
  using message_type = xg::garray<GaloisField>;
  using matrix_type = xg::garray<GaloisField>;
  using vector_type = xg::garray<GaloisField>;

  // Constructor
  AbstractLinearCode(
      size_t length, size_t dimension,
      const std::string& default_encoder_name = "GeneratorMatrix",
      const std::string& default_decoder_name = "Syndrome",
      Metric metric = Metric::HAMMING)
      : AbstractCode<GaloisField>(length, default_encoder_name,
                                  default_decoder_name, metric),
        dimension_(dimension) {}

  virtual ~AbstractLinearCode() = default;

  // Override dimension
  size_t Dimension() const override { return dimension_; }

  // Linear code specific methods
  virtual matrix_type GeneratorMatrix() const = 0;
  virtual matrix_type ParityCheckMatrix() const = 0;

  // Syndrome computation
  virtual vector_type Syndrome(const codeword_type& word) const;

  // Dual code
  virtual std::unique_ptr<AbstractLinearCode<GaloisField>> DualCode() const = 0;

  // Check if word is in the code using syndrome
  bool Contains(const codeword_type& word) const override;

  // Get all codewords (enumeration)
  std::vector<codeword_type> GetCodewords() const override;

  // Random codeword generation
  codeword_type RandomCodeword() const override;

  // Information rate
  double Rate() const {
    return static_cast<double>(dimension_) / this->length_;
  }

  // Redundancy
  size_t Redundancy() const { return this->length_ - dimension_; }

  // Linear combination of codewords
  virtual codeword_type LinearCombination(
      const std::vector<codeword_type>& codewords,
      const xt::xarray<GaloisField>& coefficients) const;

 protected:
  size_t dimension_;

  // Convert between vector types
  vector_type ToVector(const codeword_type& word) const;
  codeword_type ToCodeword(const vector_type& vec) const;

  // Matrix operations helpers
  matrix_type ComputeParityCheckFromGenerator(
      const matrix_type& generator) const;
  matrix_type ComputeGeneratorFromParityCheck(
      const matrix_type& parity_check) const;
};

// Implementation of template methods

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::vector_type
AbstractLinearCode<GaloisField>::Syndrome(const codeword_type& word) const {
  auto H = ParityCheckMatrix();
  auto w = ToVector(word);
  return xg::linalg::dot(H, w);
}

template <typename GaloisField>
bool AbstractLinearCode<GaloisField>::Contains(
    const codeword_type& word) const {
  if (word.size() != this->length_) {
    return false;
  }

  auto syndrome = Syndrome(word);
  auto field = this->Field();

  // Check if syndrome is zero
  for (size_t i = 0; i < syndrome.size(); ++i) {
    if (syndrome(i) != element_type{}) {
      return false;
    }
  }
  return true;
}

template <typename GaloisField>
std::vector<typename AbstractLinearCode<GaloisField>::codeword_type>
AbstractLinearCode<GaloisField>::GetCodewords() const {
  std::vector<codeword_type> codewords;
  auto field = this->Field();
  auto generator = GeneratorMatrix();

  // Generate all possible messages
  size_t q = field->Order();
  size_t total_messages = 1;
  for (size_t i = 0; i < dimension_; ++i) {
    total_messages *= q;
  }

  for (size_t i = 0; i < total_messages; ++i) {
    message_type message({dimension_});
    size_t temp = i;

    // Convert i to base-q representation
    for (size_t j = 0; j < dimension_; ++j) {
      message[j] = element_type(temp % q);
      temp /= q;
    }

    // Encode message
    auto msg_vec = ToVector(message);
    msg_vec.reshape({1, dimension_});
    auto codeword_vec_2d = xg::linalg::dot(msg_vec, generator);
    codeword_vec_2d.reshape({this->length_});
    codewords.push_back(ToCodeword(codeword_vec_2d));
  }

  return codewords;
}

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::codeword_type
AbstractLinearCode<GaloisField>::RandomCodeword() const {
  auto field = this->Field();
  message_type message({dimension_});

  // Generate random message
  for (size_t i = 0; i < dimension_; ++i) {
    message[i] = field->Random();
  }

  // Encode using default encoder
  return this->Encode(message);
}

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::codeword_type
AbstractLinearCode<GaloisField>::LinearCombination(
    const std::vector<codeword_type>& codewords,
    const xt::xarray<GaloisField>& coefficients) const {
  if (codewords.size() != coefficients.size()) {
    throw std::invalid_argument(
        "Number of codewords must match number of coefficients");
  }

  if (codewords.empty()) {
    return xt::zeros<GaloisField>({this->length_});
  }

  auto field = this->Field();
  codeword_type result = xt::zeros<GaloisField>({this->length_});

  for (size_t i = 0; i < codewords.size(); ++i) {
    if (codewords[i].size() != this->length_) {
      throw std::invalid_argument("All codewords must have the same length");
    }

    for (size_t j = 0; j < this->length_; ++j) {
      auto term = field->Mul(coefficients(i), codewords[i](j));
      result(j) = field->Add(result(j), term);
    }
  }

  return result;
}

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::vector_type
AbstractLinearCode<GaloisField>::ToVector(const codeword_type& word) const {
  return vector_type(word);
}

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::codeword_type
AbstractLinearCode<GaloisField>::ToCodeword(const vector_type& vec) const {
  return vec;
}

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::matrix_type
AbstractLinearCode<GaloisField>::ComputeParityCheckFromGenerator(
    const matrix_type& generator) const {
  // For a generator matrix G of size k x n, the parity check matrix H is (n-k)
  // x n such that G * H^T = 0

  // This is a simplified implementation - in practice, you'd want to use
  // proper matrix operations to compute the null space
  size_t k = generator.shape(0);
  size_t n = generator.shape(1);

  if (k >= n) {
    throw std::invalid_argument(
        "Generator matrix must have more columns than rows");
  }

  // For systematic form [I_k | P], H = [-P^T | I_{n-k}]
  // This is a simplified approach - real implementation would handle general
  // case
  matrix_type H({n - k, n});

  // This needs proper implementation based on your matrix library
  // For now, return identity matrix as placeholder
  auto field = this->Field();
  for (size_t i = 0; i < n - k; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i == j - k) {
        H(i, j) = element_type(1);
      } else {
        H(i, j) = element_type{};
      }
    }
  }

  return H;
}

template <typename GaloisField>
typename AbstractLinearCode<GaloisField>::matrix_type
AbstractLinearCode<GaloisField>::ComputeGeneratorFromParityCheck(
    const matrix_type& parity_check) const {
  // This is the dual operation - compute generator from parity check
  // For a parity check matrix H of size (n-k) x n, find generator G of size k x
  // n such that G * H^T = 0

  // This is a placeholder implementation
  size_t n_minus_k = parity_check.shape(0);
  size_t n = parity_check.shape(1);
  size_t k = n - n_minus_k;

  matrix_type G({k, n});

  // TODO(amirmulla): This needs proper implementation based on your matrix
  // library. For now, return identity matrix as placeholder.
  for (size_t i = 0; i < k; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i == j) {
        G(i, j) = element_type(1);
      } else {
        G(i, j) = element_type{};
      }
    }
  }

  return G;
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ABSTRACT_LINEAR_CODE_HPP
