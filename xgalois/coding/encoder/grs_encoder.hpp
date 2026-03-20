#ifndef XGALOIS_CODING_GRS_ENCODER_HPP
#define XGALOIS_CODING_GRS_ENCODER_HPP

#include <xtensor/containers/xarray.hpp>

#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename ElementType>
class GRSEncoder : public Encoder<ElementType> {
 public:
  using element_type = ElementType;
  using codeword_type = xt::xarray<ElementType>;
  using message_type = xt::xarray<ElementType>;
  using polynomial_type = xg::poly::PolyDense<ElementType>;

  // Constructor
  explicit GRSEncoder(const AbstractCode<ElementType>* code)
      : Encoder<ElementType>(code) {
    // Try to cast to GRS code
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<ElementType>*>(code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSEncoder requires a GeneralizedReedSolomonCode");
    }
  }

  // Encode message as polynomial evaluation
  codeword_type Encode(const message_type& message) const override {
    if (message.size() != MessageLength()) {
      throw std::invalid_argument("Message length must match code dimension");
    }

    // Create polynomial from message coefficients
    std::vector<ElementType> message_vec(message.size());
    for (size_t i = 0; i < message.size(); ++i) {
      message_vec[i] = message(i);
    }
    polynomial_type message_poly(message_vec, grs_code_->Field());

    // Encode using polynomial evaluation
    return grs_code_->EncodePolynomial(message_poly);
  }

  // Unencode: interpolate polynomial from codeword
  message_type Unencode(const codeword_type& codeword) const override {
    if (codeword.size() != this->code_->Length()) {
      throw std::invalid_argument("Codeword length must match code length");
    }

    if (!this->code_->Contains(codeword)) {
      throw std::invalid_argument("Input is not a valid codeword");
    }

    // Interpolate polynomial from codeword
    auto poly = grs_code_->DecodeToPolynomial(codeword);

    // Extract coefficients as message
    message_type message = xt::zeros<ElementType>({MessageLength()});
    for (size_t i = 0; i < MessageLength(); ++i) {
      if (i <= poly.Degree()) {
        message(i) = poly.GetCoefficient(i);
      } else {
        message(i) = ElementType{};
      }
    }

    return message;
  }

  size_t MessageLength() const override { return grs_code_->Dimension(); }

  std::string ToString() const override {
    return "GRS Encoder for [" + std::to_string(this->code_->Length()) + ", " +
           std::to_string(grs_code_->Dimension()) + "] GRS code";
  }

  // GRS-specific encoding methods

  // Encode polynomial directly
  codeword_type EncodePolynomial(const polynomial_type& message_poly) const {
    return grs_code_->EncodePolynomial(message_poly);
  }

  // Systematic encoding (message appears in first k positions)
  codeword_type SystematicEncode(const message_type& message) const {
    if (message.size() != MessageLength()) {
      throw std::invalid_argument("Message length must match code dimension");
    }

    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    auto multipliers = grs_code_->ColumnMultipliers();

    // For systematic encoding, we need to solve for the polynomial
    // such that f(α_0) = m_0/v_0, f(α_1) = m_1/v_1, ..., f(α_{k-1}) =
    // m_{k-1}/v_{k-1}

    xt::xarray<ElementType> points = xt::zeros<ElementType>({MessageLength()});
    xt::xarray<ElementType> values = xt::zeros<ElementType>({MessageLength()});

    for (size_t i = 0; i < MessageLength(); ++i) {
      points(i) = eval_points(i);
      values(i) = field->Div(message(i), multipliers(i));
    }

    // Interpolate polynomial
    auto message_poly = LagrangeInterpolation(points, values);

    // Encode polynomial
    return EncodePolynomial(message_poly);
  }

 private:
  const GeneralizedReedSolomonCode<ElementType>* grs_code_;

  // Lagrange interpolation helper
  polynomial_type LagrangeInterpolation(
      const xt::xarray<ElementType>& points,
      const xt::xarray<ElementType>& values) const {
    if (points.size() != values.size()) {
      throw std::invalid_argument("Points and values must have same size");
    }

    auto field = grs_code_->Field();
    size_t n = points.size();
    xt::xarray<ElementType> result_coeffs = xt::zeros<ElementType>({n});

    for (size_t i = 0; i < n; ++i) {
      // Compute Lagrange basis polynomial L_i(x)
      xt::xarray<ElementType> basis_coeffs = xt::ones<ElementType>({1});

      for (size_t j = 0; j < n; ++j) {
        if (i != j) {
          // Multiply by (x - points[j]) / (points[i] - points[j])
          ElementType denominator = field->Sub(points(i), points(j));
          ElementType inv_denom = field->Inv(denominator);

          // Multiply basis_coeffs by (x - points[j])
          xt::xarray<ElementType> new_coeffs =
              xt::zeros<ElementType>({basis_coeffs.size() + 1});
          for (size_t k = 0; k < basis_coeffs.size(); ++k) {
            new_coeffs(k) =
                field->Add(new_coeffs(k),
                           field->Mul(basis_coeffs(k), field->Neg(points(j))));
            new_coeffs(k + 1) = field->Add(new_coeffs(k + 1), basis_coeffs(k));
          }

          // Multiply by 1 / (points[i] - points[j])
          for (size_t k = 0; k < new_coeffs.size(); ++k) {
            new_coeffs(k) = field->Mul(new_coeffs(k), inv_denom);
          }

          basis_coeffs = new_coeffs;
        }
      }

      // Add values[i] * L_i(x) to result
      for (size_t k = 0; k < basis_coeffs.size() && k < result_coeffs.size();
           ++k) {
        ElementType term = field->Mul(values(i), basis_coeffs(k));
        result_coeffs(k) = field->Add(result_coeffs(k), term);
      }
    }

    // Convert xtensor array to std::vector for polynomial constructor
    std::vector<ElementType> coeffs_vec(result_coeffs.size());
    for (size_t i = 0; i < result_coeffs.size(); ++i) {
      coeffs_vec[i] = result_coeffs(i);
    }

    return polynomial_type(coeffs_vec, field);
  }
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_GRS_ENCODER_HPP
