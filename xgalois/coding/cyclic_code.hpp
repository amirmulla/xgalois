#ifndef XGALOIS_CODING_CYCLIC_CODE_HPP
#define XGALOIS_CODING_CYCLIC_CODE_HPP

#include <memory>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class CyclicEncoder;
template <typename GaloisField>
class CyclicSyndromeDecoder;

template <typename GaloisField>
class CyclicCode : public AbstractLinearCode<GaloisField> {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xg::garray<GaloisField>;
  using message_type = xg::garray<GaloisField>;
  using matrix_type = xg::garray<GaloisField>;
  using vector_type = xg::garray<GaloisField>;
  using polynomial_type = xg::PolynomialDense<GaloisField>;

  CyclicCode(std::shared_ptr<GaloisField> field, size_t length,
             const polynomial_type& generator_poly)
      : AbstractLinearCode<GaloisField>(
            length, length - generator_poly.Degree(), "CyclicEncoder",
            "CyclicSyndrome", Metric::HAMMING),
        field_(field),
        generator_poly_(generator_poly) {
    if (generator_poly.Degree() >= length) {
      throw std::invalid_argument(
          "Generator polynomial degree must be less than code length");
    }

    parity_check_poly_ = ComputeParityCheckPolynomial();

    BuildMatrices();

    RegisterEncodersAndDecoders();
  }

  static std::unique_ptr<CyclicCode<GaloisField>> FromParityCheckPolynomial(
      std::shared_ptr<GaloisField> field, size_t length,
      const polynomial_type& parity_check_poly) {
    auto generator_poly =
        ComputeGeneratorFromParityCheck(field, length, parity_check_poly);
    return std::make_unique<CyclicCode>(field, length, generator_poly);
  }

  std::shared_ptr<GaloisField> Field() const override {
    return field_;
  }

  matrix_type GeneratorMatrix() const override { return generator_matrix_; }

  matrix_type ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  size_t MinimumDistance() const override {

    return ComputeMinimumDistance();
  }

  std::unique_ptr<AbstractLinearCode<GaloisField>> DualCode() const override {

    return std::make_unique<CyclicCode>(field_, this->length_,
                                        parity_check_poly_);
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + "] cyclic code over GF(" +
           std::to_string(field_->Order()) + ")";
  }

  const polynomial_type& GeneratorPolynomial() const { return generator_poly_; }
  const polynomial_type& ParityCheckPolynomial() const {
    return parity_check_poly_;
  }

  bool IsBCH() const {

    return false;
  }

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

  polynomial_type VectorToPolynomial(const codeword_type& vector) const {
    return polynomial_type(vector, field_);
  }

  codeword_type PolynomialToVector(const polynomial_type& poly) const {
    codeword_type result = xt::zeros<GaloisField>({this->length_});
    for (size_t i = 0; i <= poly.Degree() && i < this->length_; ++i) {
      result[i] = poly.GetCoefficient(i);
    }
    return result;
  }

 private:
  std::shared_ptr<GaloisField> field_;
  polynomial_type generator_poly_;
  polynomial_type parity_check_poly_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;

  polynomial_type ComputeParityCheckPolynomial() const {

    std::vector<GaloisField> xn_minus_1_coeffs(this->length_ + 1);

    for (auto& c : xn_minus_1_coeffs) {
      c = element_type(0, field_);
    }

    xn_minus_1_coeffs[0] = field_->Neg(element_type(1, field_));
    xn_minus_1_coeffs[this->length_] = element_type(1, field_);

    polynomial_type xn_minus_1(xn_minus_1_coeffs, field_);

    return xn_minus_1.Div(generator_poly_);
  }

  static polynomial_type ComputeGeneratorFromParityCheck(
      std::shared_ptr<GaloisField> field, size_t length,
      const polynomial_type& parity_check_poly) {

    std::vector<GaloisField> xn_minus_1_coeffs(length + 1);
    for (auto& c : xn_minus_1_coeffs) {
      c = element_type(0, field);
    }

    xn_minus_1_coeffs[0] = field->Neg(element_type(1, field));
    xn_minus_1_coeffs[length] = element_type(1, field);

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

    for (size_t i = 0; i < k; ++i) {

      std::vector<GaloisField> xi_coeffs(i + 1);
      for (auto& c : xi_coeffs) {
        c = element_type(0, field_);
      }
      xi_coeffs[i] = element_type(1, field_);
      polynomial_type xi(xi_coeffs, field_);

      auto row_poly = xi.Mul(generator_poly_);

      row_poly = ReduceModuloXnMinus1(row_poly);

      for (size_t j = 0; j < n; ++j) {
        if (j <= row_poly.Degree()) {
          generator_matrix_.Set(i, j, row_poly.GetCoefficient(j));
        } else {
          generator_matrix_.Set(i, j, element_type{});
        }
      }
    }
  }

  void BuildParityCheckMatrix() {
    size_t k = this->dimension_;
    size_t n = this->length_;
    size_t r = n - k;

    parity_check_matrix_ = matrix_type(r, n);

    for (size_t i = 0; i < r; ++i) {

      std::vector<GaloisField> xi_coeffs(i + 1);
      for (auto& c : xi_coeffs) {
        c = element_type(0, field_);
      }
      xi_coeffs[i] = element_type(1, field_);
      polynomial_type xi(xi_coeffs, field_);

      auto row_poly = xi.Mul(parity_check_poly_);

      row_poly = ReduceModuloXnMinus1(row_poly);

      for (size_t j = 0; j < n; ++j) {
        if (j <= row_poly.Degree()) {
          parity_check_matrix_.Set(i, j, row_poly.GetCoefficient(j));
        } else {
          parity_check_matrix_.Set(i, j, element_type{});
        }
      }
    }
  }

  polynomial_type ReduceModuloXnMinus1(const polynomial_type& poly) const {

    std::vector<GaloisField> coeffs(this->length_, element_type{});

    for (size_t i = 0; i <= poly.Degree(); ++i) {
      size_t pos = i % this->length_;
      coeffs[pos] = field_->Add(coeffs[pos], poly.GetCoefficient(i));
    }

    return polynomial_type(coeffs, field_);
  }

  size_t ComputeMinimumDistance() const {

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

    this->RegisterEncoder(
        "CyclicEncoder", [](const AbstractCode<GaloisField>* code) {
          return std::make_unique<CyclicEncoder<GaloisField>>(code);
        });

    this->RegisterDecoder(
        "CyclicSyndrome", [](const AbstractCode<GaloisField>* code) {
          return std::make_unique<CyclicSyndromeDecoder<GaloisField>>(code);
        });
  }
};

}
}

#endif
