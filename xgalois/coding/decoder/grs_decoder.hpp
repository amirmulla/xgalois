#ifndef XGALOIS_CODING_GRS_DECODER_HPP
#define XGALOIS_CODING_GRS_DECODER_HPP

#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class GRSDecoder : public Decoder<GaloisField> {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<element_type>;
  using message_type = xt::xarray<element_type>;
  using polynomial_type = xg::PolynomialDense<GaloisField>;

  explicit GRSDecoder(const AbstractCode<GaloisField>* code)
      : Decoder<GaloisField>(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField>*>(code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSDecoder requires a GeneralizedReedSolomonCode");
    }

    max_errors_ = (grs_code_->MinimumDistance() - 1) / 2;
  }

  codeword_type DecodeToCode(
      const codeword_type& received_word) const override {
    if (received_word.size() != this->code_->Length()) {
      throw std::invalid_argument(
          "Received word length must match code length");
    }

    auto syndrome = ComputeSyndrome(received_word);

    if (IsSyndromeZero(syndrome)) {
      return received_word;
    }

    auto error_locator = FindErrorLocator(syndrome);

    auto error_locations = FindErrorLocations(error_locator);

    auto error_values =
        FindErrorValues(syndrome, error_locator, error_locations);

    return CorrectErrors(received_word, error_locations, error_values);
  }

  message_type DecodeToMessage(
      const codeword_type& received_word) const override {
    auto corrected_codeword = DecodeToCode(received_word);

    auto encoder = this->code_->GetEncoder();
    return encoder->Unencode(corrected_codeword);
  }

  size_t MaxErrors() const { return max_errors_; }

  std::string ToString() const override {
    return "GRS Decoder (Peterson-Gorenstein-Zierler) for [" +
           std::to_string(this->code_->Length()) + ", " +
           std::to_string(grs_code_->Dimension()) + "] GRS code" +
           " (max errors: " + std::to_string(max_errors_) + ")";
  }

 private:
  const GeneralizedReedSolomonCode<GaloisField>* grs_code_;
  size_t max_errors_;

  std::vector<element_type> ComputeSyndrome(
      const codeword_type& received_word) const {
    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    auto multipliers = grs_code_->ColumnMultipliers();

    size_t syndrome_length = 2 * max_errors_;
    std::vector<element_type> syndrome(syndrome_length);

    for (size_t i = 0; i < syndrome_length; ++i) {
      syndrome[i] = element_type(0, field);
      for (size_t j = 0; j < received_word.size(); ++j) {
        element_type power = field->Pow(eval_points(j), i);
        element_type term = field->Mul(received_word(j), power);
        syndrome[i] = field->Add(syndrome[i], term);
      }
    }

    return syndrome;
  }

  bool IsSyndromeZero(const std::vector<element_type>& syndrome) const {
    for (const auto& s : syndrome) {
      if (s.Value() != 0) {
        return false;
      }
    }
    return true;
  }

  polynomial_type FindErrorLocator(
      const std::vector<element_type>& syndrome) const {
    auto field = grs_code_->Field();

    for (size_t t = 1; t <= max_errors_; ++t) {
      if (2 * t > syndrome.size()) continue;

      std::vector<std::vector<element_type>> matrix(
          t, std::vector<element_type>(t));
      std::vector<element_type> rhs(t);

      for (size_t i = 0; i < t; ++i) {
        for (size_t j = 0; j < t; ++j) {
          matrix[i][j] = syndrome[i + j];
        }
        rhs[i] = field->Neg(syndrome[i + t]);
      }

      auto solution = SolveLinearSystem(matrix, rhs);

      if (!solution.empty()) {
        std::vector<element_type> coeffs(t + 1);
        coeffs[0] = element_type(1, field);
        for (size_t i = 0; i < t; ++i) {
          coeffs[i + 1] = solution[i];
        }

        return polynomial_type(coeffs);
      }
    }

    return polynomial_type({element_type(1, field)});
  }

  std::vector<size_t> FindErrorLocations(
      const polynomial_type& error_locator) const {
    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    std::vector<size_t> error_locations;

    for (size_t i = 0; i < eval_points.size(); ++i) {
      element_type inv_alpha = field->Inv(eval_points[i]);
      if (error_locator.Evaluate(inv_alpha) == element_type{}) {
        error_locations.push_back(i);
      }
    }

    return error_locations;
  }

  std::vector<element_type> FindErrorValues(
      const std::vector<element_type>& syndrome,
      const polynomial_type& error_locator,
      const std::vector<size_t>& error_locations) const {
    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    auto multipliers = grs_code_->ColumnMultipliers();

    std::vector<element_type> error_values(error_locations.size());

    auto error_evaluator = ComputeErrorEvaluator(syndrome, error_locator);

    auto error_locator_derivative = error_locator.Derivative();

    for (size_t i = 0; i < error_locations.size(); ++i) {
      size_t loc = error_locations[i];
      element_type alpha_inv = field->Inv(eval_points[loc]);

      element_type numerator = error_evaluator.Evaluate(alpha_inv);
      element_type denominator = error_locator_derivative.Evaluate(alpha_inv);

      if (denominator == element_type{}) {
        throw std::runtime_error(
            "Error in Forney's algorithm: zero denominator");
      }

      element_type error_value = field->Div(numerator, denominator);
      error_value = field->Mul(error_value, field->Neg(eval_points[loc]));
      error_value = field->Div(error_value, multipliers[loc]);

      error_values[i] = error_value;
    }

    return error_values;
  }

  polynomial_type ComputeErrorEvaluator(
      const std::vector<element_type>& syndrome,
      const polynomial_type& error_locator) const {
    auto field = grs_code_->Field();

    std::vector<element_type> syndrome_coeffs(syndrome.size() + 1,
                                              element_type(0, field));
    for (size_t i = 0; i < syndrome.size(); ++i) {
      syndrome_coeffs[i + 1] = syndrome[i];
    }
    polynomial_type syndrome_poly(syndrome_coeffs);

    auto product = syndrome_poly.Mul(error_locator);

    size_t max_degree = 2 * max_errors_;
    std::vector<GaloisField> evaluator_coeffs(max_degree, element_type{});
    for (size_t i = 0; i < max_degree && i <= product.Degree(); ++i) {
      evaluator_coeffs[i] = product.GetCoefficient(i);
    }

    return polynomial_type(evaluator_coeffs, field);
  }

  codeword_type CorrectErrors(
      const codeword_type& received_word,
      const std::vector<size_t>& error_locations,
      const std::vector<element_type>& error_values) const {
    auto field = grs_code_->Field();
    codeword_type corrected = received_word;

    for (size_t i = 0; i < error_locations.size(); ++i) {
      size_t loc = error_locations[i];
      corrected(loc) = field->Sub(corrected(loc), error_values[i]);
    }

    return corrected;
  }

  std::vector<element_type> SolveLinearSystem(
      const std::vector<std::vector<element_type>>& A,
      const std::vector<element_type>& b) const {
    auto field = grs_code_->Field();
    size_t n = A.size();

    if (n == 0 || A[0].size() != n || b.size() != n) {
      return {};
    }

    std::vector<std::vector<element_type>> augmented(
        n, std::vector<element_type>(n + 1));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        augmented[i][j] = A[i][j];
      }
      augmented[i][n] = b[i];
    }

    for (size_t i = 0; i < n; ++i) {
      size_t pivot = i;
      for (size_t j = i + 1; j < n; ++j) {
        if (augmented[j][i].Value() != 0) {
          pivot = j;
          break;
        }
      }

      if (augmented[pivot][i].Value() == 0) {
        return {};
      }

      if (pivot != i) {
        std::swap(augmented[i], augmented[pivot]);
      }

      element_type pivot_inv = field->Inv(augmented[i][i]);
      for (size_t j = i + 1; j < n; ++j) {
        if (augmented[j][i] != element_type{}) {
          element_type factor = field->Mul(augmented[j][i], pivot_inv);
          for (size_t k = i; k <= n; ++k) {
            element_type term = field->Mul(factor, augmented[i][k]);
            augmented[j][k] = field->Sub(augmented[j][k], term);
          }
        }
      }
    }

    std::vector<GaloisField> solution(n);
    for (int i = n - 1; i >= 0; --i) {
      solution[i] = augmented[i][n];
      for (size_t j = i + 1; j < n; ++j) {
        element_type term = field->Mul(augmented[i][j], solution[j]);
        solution[i] = field->Sub(solution[i], term);
      }
      solution[i] = field->Div(solution[i], augmented[i][i]);
    }

    return solution;
  }
};

}  // namespace coding
}  // namespace xg

#endif
