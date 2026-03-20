#ifndef XGALOIS_CODING_CYCLIC_CODE_HPP
#define XGALOIS_CODING_CYCLIC_CODE_HPP

#include <memory>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename ElementType>
class CyclicCode : public AbstractLinearCode<ElementType> {
 public:
  using element_type = ElementType;
  using codeword_type = std::vector<ElementType>;
  using message_type = std::vector<ElementType>;
  using matrix_type = xg::linalg::Matrix<ElementType>;
  using vector_type = xg::linalg::Vector<ElementType>;
  using polynomial_type = xg::poly::PolyDense<ElementType>;

  // Constructor from generator polynomial
  CyclicCode(std::shared_ptr<GaloisFieldBase<ElementType>> field, size_t length,
             const polynomial_type& generator_poly)
      : AbstractLinearCode<ElementType>(
            length, length - generator_poly.Degree(), "CyclicEncoder",
            "CyclicSyndrome", Metric::HAMMING),
        field_(field),
        generator_poly_(generator_poly) {
    if (generator_poly.Degree() >= length) {
      throw std::invalid_argument(
          "Generator polynomial degree must be less than code length");
    }

    // Compute parity check polynomial
    parity_check_poly_ = ComputeParityCheckPolynomial();

    // Build generator and parity check matrices
    BuildMatrices();

    // Register encoders and decoders
    RegisterEncodersAndDecoders();
  }

  // Constructor from parity check polynomial
  static std::unique_ptr<CyclicCode<ElementType>> FromParityCheckPolynomial(
      std::shared_ptr<GaloisFieldBase<ElementType>> field, size_t length,
      const polynomial_type& parity_check_poly) {
    auto generator_poly =
        ComputeGeneratorFromParityCheck(field, length, parity_check_poly);
    return std::make_unique<CyclicCode>(field, length, generator_poly);
  }

  // Override methods from AbstractLinearCode
  std::shared_ptr<GaloisFieldBase<ElementType>> Field() const override {
    return field_;
  }

  matrix_type GeneratorMatrix() const override { return generator_matrix_; }

  matrix_type ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  size_t MinimumDistance() const override {
    // For cyclic codes, minimum distance computation is complex
    // This is a placeholder - in practice, you'd use specific algorithms
    return ComputeMinimumDistance();
  }

  std::unique_ptr<AbstractLinearCode<ElementType>> DualCode() const override {
    // The dual of a cyclic code is also cyclic
    return std::make_unique<CyclicCode>(field_, this->length_,
                                        parity_check_poly_);
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + "] cyclic code over GF(" +
           std::to_string(field_->Order()) + ")";
  }

  // Cyclic code specific methods
  const polynomial_type& GeneratorPolynomial() const { return generator_poly_; }
  const polynomial_type& ParityCheckPolynomial() const {
    return parity_check_poly_;
  }

  // Check if code is a BCH code
  bool IsBCH() const {
    // Implementation would check if generator polynomial has consecutive roots
    return false;  // Placeholder
  }

  // Cyclic shift operation
  codeword_type CyclicShift(const codeword_type& codeword,
                            int shift = 1) const {
    if (codeword.size() != this->length_) {
      throw std::invalid_argument("Codeword length must match code length");
    }

    codeword_type shifted(this->length_);
    for (size_t i = 0; i < this->length_; ++i) {
      size_t new_pos = (i + shift + this->length_) % this->length_;
      shifted[new_pos] = codeword[i];
    }
    return shifted;
  }

  // Convert between polynomial and vector representations
  polynomial_type VectorToPolynomial(const codeword_type& vector) const {
    return polynomial_type(vector, field_);
  }

  codeword_type PolynomialToVector(const polynomial_type& poly) const {
    codeword_type result(this->length_, ElementType{});
    for (size_t i = 0; i <= poly.Degree() && i < this->length_; ++i) {
      result[i] = poly.GetCoefficient(i);
    }
    return result;
  }

 private:
  std::shared_ptr<GaloisFieldBase<ElementType>> field_;
  polynomial_type generator_poly_;
  polynomial_type parity_check_poly_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;

  polynomial_type ComputeParityCheckPolynomial() const {
    // h(x) = (x^n - 1) / g(x)
    // Create polynomial x^n - 1
    std::vector<ElementType> xn_minus_1_coeffs(this->length_ + 1,
                                               ElementType{});
    xn_minus_1_coeffs[0] = field_->Neg(ElementType(1));  // -1
    xn_minus_1_coeffs[this->length_] = ElementType(1);   // x^n

    polynomial_type xn_minus_1(xn_minus_1_coeffs, field_);

    // Compute polynomial division
    return xn_minus_1.Div(generator_poly_);
  }

  static polynomial_type ComputeGeneratorFromParityCheck(
      std::shared_ptr<GaloisFieldBase<ElementType>> field, size_t length,
      const polynomial_type& parity_check_poly) {
    // g(x) = (x^n - 1) / h(x)
    std::vector<ElementType> xn_minus_1_coeffs(length + 1, ElementType{});
    xn_minus_1_coeffs[0] = field->Neg(ElementType(1));  // -1
    xn_minus_1_coeffs[length] = ElementType(1);         // x^n

    polynomial_type xn_minus_1(xn_minus_1_coeffs, field);

    return xn_minus_1.Div(parity_check_poly);
  }

  void BuildMatrices() {
    BuildGeneratorMatrix();
    BuildParityCheckMatrix();
  }

  void BuildGeneratorMatrix() {
    size_t k = this->dimension_;
    size_t n = this->length_;

    generator_matrix_ = matrix_type(k, n);

    // Each row i is x^i * g(x) mod (x^n - 1)
    for (size_t i = 0; i < k; ++i) {
      // Create polynomial x^i
      std::vector<ElementType> xi_coeffs(i + 1, ElementType{});
      xi_coeffs[i] = ElementType(1);
      polynomial_type xi(xi_coeffs, field_);

      // Compute x^i * g(x)
      auto row_poly = xi.Mul(generator_poly_);

      // Reduce modulo x^n - 1
      row_poly = ReduceModuloXnMinus1(row_poly);

      // Fill matrix row
      for (size_t j = 0; j < n; ++j) {
        if (j <= row_poly.Degree()) {
          generator_matrix_.Set(i, j, row_poly.GetCoefficient(j));
        } else {
          generator_matrix_.Set(i, j, ElementType{});
        }
      }
    }
  }

  void BuildParityCheckMatrix() {
    size_t k = this->dimension_;
    size_t n = this->length_;
    size_t r = n - k;

    parity_check_matrix_ = matrix_type(r, n);

    // Each row i is x^i * h(x) mod (x^n - 1)
    for (size_t i = 0; i < r; ++i) {
      // Create polynomial x^i
      std::vector<ElementType> xi_coeffs(i + 1, ElementType{});
      xi_coeffs[i] = ElementType(1);
      polynomial_type xi(xi_coeffs, field_);

      // Compute x^i * h(x)
      auto row_poly = xi.Mul(parity_check_poly_);

      // Reduce modulo x^n - 1
      row_poly = ReduceModuloXnMinus1(row_poly);

      // Fill matrix row
      for (size_t j = 0; j < n; ++j) {
        if (j <= row_poly.Degree()) {
          parity_check_matrix_.Set(i, j, row_poly.GetCoefficient(j));
        } else {
          parity_check_matrix_.Set(i, j, ElementType{});
        }
      }
    }
  }

  polynomial_type ReduceModuloXnMinus1(const polynomial_type& poly) const {
    // Reduce polynomial modulo x^n - 1
    std::vector<ElementType> coeffs(this->length_, ElementType{});

    for (size_t i = 0; i <= poly.Degree(); ++i) {
      size_t pos = i % this->length_;
      coeffs[pos] = field_->Add(coeffs[pos], poly.GetCoefficient(i));
    }

    return polynomial_type(coeffs, field_);
  }

  size_t ComputeMinimumDistance() const {
    // This is a placeholder implementation
    // In practice, you'd use specific algorithms for cyclic codes

    // For now, compute by brute force (only feasible for small codes)
    size_t min_distance = this->length_;
    auto codewords = this->GetCodewords();

    for (size_t i = 1; i < codewords.size(); ++i) {
      size_t weight = this->Weight(codewords[i]);
      if (weight < min_distance) {
        min_distance = weight;
      }
    }

    return min_distance;
  }

  void RegisterEncodersAndDecoders() {
    // Register cyclic encoder
    this->RegisterEncoder(
        "CyclicEncoder", [](const AbstractCode<ElementType>* code) {
          return std::make_unique<CyclicEncoder<ElementType>>(code);
        });

    // Register cyclic syndrome decoder
    this->RegisterDecoder(
        "CyclicSyndrome", [](const AbstractCode<ElementType>* code) {
          return std::make_unique<CyclicSyndromeDecoder<ElementType>>(code);
        });
  }
};

// Forward declarations for cyclic-specific encoder/decoder
template <typename ElementType>
class CyclicEncoder;
template <typename ElementType>
class CyclicSyndromeDecoder;

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_CYCLIC_CODE_HPP
