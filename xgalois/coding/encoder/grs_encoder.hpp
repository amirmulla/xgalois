#ifndef XGALOIS_CODING_GRS_ENCODER_HPP
#define XGALOIS_CODING_GRS_ENCODER_HPP

#include <xtensor/containers/xarray.hpp>

#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class GRSEncoder : public Encoder<GaloisField> {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<GaloisField>;
  using message_type = xt::xarray<GaloisField>;
  using polynomial_type = xg::PolynomialDense<GaloisField>;

  explicit GRSEncoder(const AbstractCode<GaloisField>* code)
      : Encoder<GaloisField>(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField>*>(code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSEncoder requires a GeneralizedReedSolomonCode");
    }
  }

  codeword_type Encode(const message_type& message) const override {
    if (message.size() != MessageLength()) {
      throw std::invalid_argument("Message length must match code dimension");
    }

    std::vector<GaloisField> message_vec(message.size());
    for (size_t i = 0; i < message.size(); ++i) {
      message_vec[i] = message(i);
    }
    polynomial_type message_poly(message_vec, grs_code_->Field());

    return grs_code_->EncodePolynomial(message_poly);
  }

  message_type Unencode(const codeword_type& codeword) const override {
    if (codeword.size() != this->code_->Length()) {
      throw std::invalid_argument("Codeword length must match code length");
    }

    if (!this->code_->Contains(codeword)) {
      throw std::invalid_argument("Input is not a valid codeword");
    }

    auto poly = grs_code_->DecodeToPolynomial(codeword);

    message_type message = xt::zeros<GaloisField>({MessageLength()});
    for (size_t i = 0; i < MessageLength(); ++i) {
      if (i <= poly.Degree()) {
        message(i) = poly.GetCoefficient(i);
      } else {
        message(i) = element_type{};
      }
    }

    return message;
  }

  size_t MessageLength() const override { return grs_code_->Dimension(); }

  std::string ToString() const override {
    return "GRS Encoder for [" + std::to_string(this->code_->Length()) + ", " +
           std::to_string(grs_code_->Dimension()) + "] GRS code";
  }

  codeword_type EncodePolynomial(const polynomial_type& message_poly) const {
    return grs_code_->EncodePolynomial(message_poly);
  }

  codeword_type SystematicEncode(const message_type& message) const {
    if (message.size() != MessageLength()) {
      throw std::invalid_argument("Message length must match code dimension");
    }

    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    auto multipliers = grs_code_->ColumnMultipliers();

    xt::xarray<GaloisField> points = xt::zeros<GaloisField>({MessageLength()});
    xt::xarray<GaloisField> values = xt::zeros<GaloisField>({MessageLength()});

    for (size_t i = 0; i < MessageLength(); ++i) {
      points(i) = eval_points(i);
      values(i) = field->Div(message(i), multipliers(i));
    }

    auto message_poly = LagrangeInterpolation(points, values);

    return EncodePolynomial(message_poly);
  }

 private:
  const GeneralizedReedSolomonCode<GaloisField>* grs_code_;

  polynomial_type LagrangeInterpolation(
      const xt::xarray<GaloisField>& points,
      const xt::xarray<GaloisField>& values) const {
    if (points.size() != values.size()) {
      throw std::invalid_argument("Points and values must have same size");
    }

    auto field = grs_code_->Field();
    size_t n = points.size();
    xt::xarray<GaloisField> result_coeffs = xt::zeros<GaloisField>({n});

    for (size_t i = 0; i < n; ++i) {
      xt::xarray<GaloisField> basis_coeffs = xt::ones<GaloisField>({1});

      for (size_t j = 0; j < n; ++j) {
        if (i != j) {
          element_type denominator = field->Sub(points(i), points(j));
          element_type inv_denom = field->Inv(denominator);

          xt::xarray<GaloisField> new_coeffs =
              xt::zeros<GaloisField>({basis_coeffs.size() + 1});
          for (size_t k = 0; k < basis_coeffs.size(); ++k) {
            new_coeffs(k) =
                field->Add(new_coeffs(k),
                           field->Mul(basis_coeffs(k), field->Neg(points(j))));
            new_coeffs(k + 1) = field->Add(new_coeffs(k + 1), basis_coeffs(k));
          }

          for (size_t k = 0; k < new_coeffs.size(); ++k) {
            new_coeffs(k) = field->Mul(new_coeffs(k), inv_denom);
          }

          basis_coeffs = new_coeffs;
        }
      }

      for (size_t k = 0; k < basis_coeffs.size() && k < result_coeffs.size();
           ++k) {
        element_type term = field->Mul(values(i), basis_coeffs(k));
        result_coeffs(k) = field->Add(result_coeffs(k), term);
      }
    }

    std::vector<GaloisField> coeffs_vec(result_coeffs.size());
    for (size_t i = 0; i < result_coeffs.size(); ++i) {
      coeffs_vec[i] = result_coeffs(i);
    }

    return polynomial_type(coeffs_vec, field);
  }
};

}  // namespace coding
}  // namespace xg

#endif
