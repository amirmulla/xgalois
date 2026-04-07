#ifndef XGALOIS_CODING_GRS_HPP
#define XGALOIS_CODING_GRS_HPP

#include <memory>
#include <vector>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/views/xview.hpp>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class GRSEncoder;
template <typename GaloisField>
class GRSDecoder;
template <typename GaloisField>
class BerlekampMasseyDecoder;

template <typename GaloisField>
class GeneralizedReedSolomonCode : public AbstractLinearCode<GaloisField> {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<GaloisField>;
  using message_type = xt::xarray<GaloisField>;
  using matrix_type = xt::xarray<GaloisField>;
  using vector_type = xt::xarray<GaloisField>;
  using polynomial_type = xg::PolynomialDense<GaloisField>;

  GeneralizedReedSolomonCode(
    std::shared_ptr<GaloisField> field,
      const xt::xarray<GaloisField>& evaluation_points,
      const xt::xarray<GaloisField>& column_multipliers, size_t dimension)
      : AbstractLinearCode<GaloisField>(evaluation_points.size(), dimension,
                                        "GRSEncoder", "GRSDecoder",
                                        Metric::HAMMING),
        field_(field),
        evaluation_points_(evaluation_points),
        column_multipliers_(column_multipliers) {
    if (evaluation_points.size() != column_multipliers.size()) {
      throw std::invalid_argument(
          "Evaluation points and column multipliers must have same size");
    }

    if (dimension > evaluation_points.size()) {
      throw std::invalid_argument(
          "Dimension cannot exceed number of evaluation points");
    }

    for (size_t i = 0; i < evaluation_points.size(); ++i) {
      for (size_t j = i + 1; j < evaluation_points.size(); ++j) {
        if (evaluation_points(i) == evaluation_points(j)) {
          throw std::invalid_argument("Evaluation points must be distinct");
        }
      }
    }

    for (size_t i = 0; i < column_multipliers.size(); ++i) {
      if (column_multipliers(i) == element_type{}) {
        throw std::invalid_argument("Column multipliers must be non-zero");
      }
    }

    BuildMatrices();

    RegisterEncodersAndDecoders();
  }

  static std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>>
  StandardReedSolomon(std::shared_ptr<GaloisField> field,
                      const xt::xarray<GaloisField>& evaluation_points,
                      size_t dimension) {
    xt::xarray<GaloisField> multipliers =
        xt::ones<GaloisField>({evaluation_points.size()});
    return std::make_unique<GeneralizedReedSolomonCode>(
        field, evaluation_points, multipliers, dimension);
  }

  static std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>> ReedSolomon(
      std::shared_ptr<GaloisField> field, size_t length,
      size_t dimension, const element_type& primitive_element) {
    xt::xarray<GaloisField> evaluation_points =
        xt::zeros<GaloisField>({length});
    xt::xarray<GaloisField> multipliers = xt::ones<GaloisField>({length});

    element_type power = element_type(1);
    for (size_t i = 0; i < length; ++i) {
      evaluation_points(i) = power;
      power = field->Mul(power, primitive_element);
    }

    return std::make_unique<GeneralizedReedSolomonCode>(
        field, evaluation_points, multipliers, dimension);
  }

  std::shared_ptr<GaloisField> Field() const override {
    return field_;
  }

  matrix_type GeneratorMatrix() const override { return generator_matrix_; }

  matrix_type ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  size_t MinimumDistance() const override {

    return this->length_ - this->dimension_ + 1;
  }

  std::unique_ptr<AbstractLinearCode<GaloisField>> DualCode() const override {

    return ComputeDualCode();
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + ", " +
           std::to_string(MinimumDistance()) + "] GRS code over GF(" +
           std::to_string(field_->Order()) + ")";
  }

  const xt::xarray<GaloisField>& EvaluationPoints() const {
    return evaluation_points_;
  }

  const xt::xarray<GaloisField>& ColumnMultipliers() const {
    return column_multipliers_;
  }

  codeword_type EncodePolynomial(const polynomial_type& message_poly) const {
    if (message_poly.Degree() >= this->dimension_) {
      throw std::invalid_argument(
          "Message polynomial degree must be less than dimension");
    }

    codeword_type codeword = xt::zeros<GaloisField>({this->length_});
    for (size_t i = 0; i < this->length_; ++i) {
      element_type eval = message_poly.Evaluate(evaluation_points_(i));
      codeword(i) = field_->Mul(eval, column_multipliers_(i));
    }

    return codeword;
  }

  polynomial_type DecodeToPolynomial(const codeword_type& received_word) const {

    auto decoder = this->GetDecoder();
    auto corrected = decoder->DecodeToCode(received_word);

    return InterpolateFromCodeword(corrected);
  }

  bool IsMDS() const {
    return MinimumDistance() == this->length_ - this->dimension_ + 1;
  }

  xt::xarray<size_t> SystematicPositions() const {
    xt::xarray<size_t> positions = xt::zeros<size_t>({this->dimension_});
    for (size_t i = 0; i < this->dimension_; ++i) {
      positions(i) = i;
    }
    return positions;
  }

 private:
  std::shared_ptr<GaloisField> field_;
  xt::xarray<GaloisField> evaluation_points_;
  xt::xarray<GaloisField> column_multipliers_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;

  void BuildMatrices() {
    BuildGeneratorMatrix();
    BuildParityCheckMatrix();
  }

  void BuildGeneratorMatrix() {
    size_t k = this->dimension_;
    size_t n = this->length_;

    generator_matrix_ = matrix_type(k, n);

    for (size_t i = 0; i < k; ++i) {
      for (size_t j = 0; j < n; ++j) {
        element_type power = field_->Pow(evaluation_points_(j), i);
        element_type entry = field_->Mul(column_multipliers_(j), power);
        generator_matrix_.Set(i, j, entry);
      }
    }
  }

  void BuildParityCheckMatrix() {
    size_t k = this->dimension_;
    size_t n = this->length_;
    size_t r = n - k;

    parity_check_matrix_ = matrix_type(r, n);

    auto dual_multipliers = ComputeDualMultipliers();

    for (size_t i = 0; i < r; ++i) {
      for (size_t j = 0; j < n; ++j) {
        element_type power = field_->Pow(evaluation_points_(j), i);
        element_type entry = field_->Mul(dual_multipliers(j), power);
        parity_check_matrix_.Set(i, j, entry);
      }
    }
  }

  xt::xarray<GaloisField> ComputeDualMultipliers() const {

    xt::xarray<GaloisField> dual_multipliers =
        xt::zeros<GaloisField>({this->length_});

    for (size_t j = 0; j < this->length_; ++j) {

      element_type product = element_type(1);
      for (size_t i = 0; i < this->length_; ++i) {
        if (i != j) {
          element_type diff =
              field_->Sub(evaluation_points_(j), evaluation_points_(i));
          product = field_->Mul(product, diff);
        }
      }

      element_type denominator = field_->Mul(column_multipliers_(j), product);
      dual_multipliers(j) = field_->Inv(denominator);
    }

    return dual_multipliers;
  }

  std::unique_ptr<AbstractLinearCode<GaloisField>> ComputeDualCode() const {
    auto dual_multipliers = ComputeDualMultipliers();
    size_t dual_dimension = this->length_ - this->dimension_;

    return std::make_unique<GeneralizedReedSolomonCode>(
        field_, evaluation_points_, dual_multipliers, dual_dimension);
  }

  polynomial_type InterpolateFromCodeword(const codeword_type& codeword) const {

    xt::xarray<GaloisField> values = xt::zeros<GaloisField>({this->length_});

    for (size_t i = 0; i < this->length_; ++i) {
      values(i) = field_->Div(codeword(i), column_multipliers_(i));
    }

    xt::xarray<GaloisField> points = xt::zeros<GaloisField>({this->dimension_});
    xt::xarray<GaloisField> vals = xt::zeros<GaloisField>({this->dimension_});

    for (size_t i = 0; i < this->dimension_; ++i) {
      points(i) = evaluation_points_(i);
      vals(i) = values(i);
    }

    return LagrangeInterpolation(points, vals);
  }

  polynomial_type LagrangeInterpolation(
      const xt::xarray<GaloisField>& points,
      const xt::xarray<GaloisField>& values) const {
    if (points.size() != values.size()) {
      throw std::invalid_argument("Points and values must have same size");
    }

    size_t n = points.size();
    xt::xarray<GaloisField> result_coeffs = xt::zeros<GaloisField>({n});

    for (size_t i = 0; i < n; ++i) {

      xt::xarray<GaloisField> basis_coeffs = xt::ones<GaloisField>({1});

      for (size_t j = 0; j < n; ++j) {
        if (i != j) {

          element_type denominator = field_->Sub(points(i), points(j));
          element_type inv_denom = field_->Inv(denominator);

          xt::xarray<GaloisField> new_coeffs =
              xt::zeros<GaloisField>({basis_coeffs.size() + 1});
          for (size_t k = 0; k < basis_coeffs.size(); ++k) {
            new_coeffs(k) = field_->Add(
                new_coeffs(k),
                field_->Mul(basis_coeffs(k), field_->Neg(points(j))));
            new_coeffs(k + 1) = field_->Add(new_coeffs(k + 1), basis_coeffs(k));
          }

          for (size_t k = 0; k < new_coeffs.size(); ++k) {
            new_coeffs(k) = field_->Mul(new_coeffs(k), inv_denom);
          }

          basis_coeffs = new_coeffs;
        }
      }

      for (size_t k = 0; k < basis_coeffs.size() && k < result_coeffs.size();
           ++k) {
        element_type term = field_->Mul(values(i), basis_coeffs(k));
        result_coeffs(k) = field_->Add(result_coeffs(k), term);
      }
    }

    std::vector<GaloisField> coeffs_vec(result_coeffs.size());
    for (size_t i = 0; i < result_coeffs.size(); ++i) {
      coeffs_vec[i] = result_coeffs(i);
    }

    return polynomial_type(coeffs_vec, field_);
  }

  void RegisterEncodersAndDecoders() {

    this->RegisterEncoder(
        "GRSEncoder", [](const AbstractCode<GaloisField>* code) {
          return std::make_unique<GRSEncoder<GaloisField>>(code);
        });

    this->RegisterDecoder(
        "GRSDecoder", [](const AbstractCode<GaloisField>* code) {
          return std::make_unique<GRSDecoder<GaloisField>>(code);
        });

    this->RegisterDecoder(
        "BerlekampMassey", [](const AbstractCode<GaloisField>* code) {
          return std::make_unique<BerlekampMasseyDecoder<GaloisField>>(code);
        });
  }
};

}
}

#endif
