#ifndef XGALOIS_CODING_GRS_DECODER_HPP
#define XGALOIS_CODING_GRS_DECODER_HPP

#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/grs.hpp"

namespace xg {
namespace coding {

template <typename ElementType>
class GRSDecoder : public Decoder<ElementType> {
 public:
  using element_type = ElementType;
  using codeword_type = std::vector<ElementType>;
  using message_type = std::vector<ElementType>;
  using polynomial_type = xg::poly::PolyDense<ElementType>;

  // Constructor
  explicit GRSDecoder(const AbstractCode<ElementType>* code)
      : Decoder<ElementType>(code) {
    // Try to cast to GRS code
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<ElementType>*>(code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSDecoder requires a GeneralizedReedSolomonCode");
    }

    // Calculate maximum correctable errors
    max_errors_ = (grs_code_->MinimumDistance() - 1) / 2;
  }

  // Decode using Peterson-Gorenstein-Zierler algorithm
  codeword_type DecodeToCode(
      const codeword_type& received_word) const override {
    if (received_word.size() != this->code_->Length()) {
      throw std::invalid_argument(
          "Received word length must match code length");
    }

    // Step 1: Compute syndrome
    auto syndrome = ComputeSyndrome(received_word);

    // Step 2: Check if syndrome is zero (no errors)
    if (IsSyndromeZero(syndrome)) {
      return received_word;
    }

    // Step 3: Find error locator polynomial using Peterson's algorithm
    auto error_locator = FindErrorLocator(syndrome);

    // Step 4: Find error locations
    auto error_locations = FindErrorLocations(error_locator);

    // Step 5: Find error values using Forney's algorithm
    auto error_values =
        FindErrorValues(syndrome, error_locator, error_locations);

    // Step 6: Correct errors
    return CorrectErrors(received_word, error_locations, error_values);
  }

  // Decode to message
  message_type DecodeToMessage(
      const codeword_type& received_word) const override {
    auto corrected_codeword = DecodeToCode(received_word);

    // Use the encoder to unencode
    auto encoder = this->code_->GetEncoder();
    return encoder->Unencode(corrected_codeword);
  }

  // Get maximum number of errors this decoder can handle
  size_t MaxErrors() const { return max_errors_; }

  std::string ToString() const override {
    return "GRS Decoder (Peterson-Gorenstein-Zierler) for [" +
           std::to_string(this->code_->Length()) + ", " +
           std::to_string(grs_code_->Dimension()) + "] GRS code" +
           " (max errors: " + std::to_string(max_errors_) + ")";
  }

 private:
  const GeneralizedReedSolomonCode<ElementType>* grs_code_;
  size_t max_errors_;

  // Compute syndrome polynomial
  std::vector<ElementType> ComputeSyndrome(
      const codeword_type& received_word) const {
    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    auto multipliers = grs_code_->ColumnMultipliers();

    size_t syndrome_length = 2 * max_errors_;
    std::vector<ElementType> syndrome(syndrome_length);

    // Syndrome S_i = sum_{j=0}^{n-1} r_j * alpha_j^i for i = 0, 1, ..., 2t-1
    for (size_t i = 0; i < syndrome_length; ++i) {
      syndrome[i] = ElementType{};
      for (size_t j = 0; j < received_word.size(); ++j) {
        ElementType power = field->Pow(eval_points[j], i);
        ElementType term = field->Mul(received_word[j], power);
        syndrome[i] = field->Add(syndrome[i], term);
      }
    }

    return syndrome;
  }

  bool IsSyndromeZero(const std::vector<ElementType>& syndrome) const {
    for (const auto& s : syndrome) {
      if (s != ElementType{}) {
        return false;
      }
    }
    return true;
  }

  // Find error locator polynomial using Peterson's algorithm
  polynomial_type FindErrorLocator(
      const std::vector<ElementType>& syndrome) const {
    auto field = grs_code_->Field();

    // Try different numbers of errors from 1 to max_errors_
    for (size_t t = 1; t <= max_errors_; ++t) {
      if (2 * t > syndrome.size()) continue;

      // Set up linear system for Peterson's algorithm
      // S * Lambda = -S_shifted
      std::vector<std::vector<ElementType>> matrix(t,
                                                   std::vector<ElementType>(t));
      std::vector<ElementType> rhs(t);

      for (size_t i = 0; i < t; ++i) {
        for (size_t j = 0; j < t; ++j) {
          matrix[i][j] = syndrome[i + j];
        }
        rhs[i] = field->Neg(syndrome[i + t]);
      }

      // Solve linear system
      auto solution = SolveLinearSystem(matrix, rhs);

      if (!solution.empty()) {
        // Found solution, construct error locator polynomial
        std::vector<ElementType> coeffs(t + 1);
        coeffs[0] = ElementType(1);  // Leading coefficient
        for (size_t i = 0; i < t; ++i) {
          coeffs[i + 1] = solution[i];
        }

        return polynomial_type(coeffs, field);
      }
    }

    // No solution found, return empty polynomial
    return polynomial_type({ElementType(1)}, field);
  }

  // Find error locations by finding roots of error locator polynomial
  std::vector<size_t> FindErrorLocations(
      const polynomial_type& error_locator) const {
    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    std::vector<size_t> error_locations;

    // Check each evaluation point
    for (size_t i = 0; i < eval_points.size(); ++i) {
      ElementType inv_alpha = field->Inv(eval_points[i]);
      if (error_locator.Evaluate(inv_alpha) == ElementType{}) {
        error_locations.push_back(i);
      }
    }

    return error_locations;
  }

  // Find error values using Forney's algorithm
  std::vector<ElementType> FindErrorValues(
      const std::vector<ElementType>& syndrome,
      const polynomial_type& error_locator,
      const std::vector<size_t>& error_locations) const {
    auto field = grs_code_->Field();
    auto eval_points = grs_code_->EvaluationPoints();
    auto multipliers = grs_code_->ColumnMultipliers();

    std::vector<ElementType> error_values(error_locations.size());

    // Compute error evaluator polynomial
    auto error_evaluator = ComputeErrorEvaluator(syndrome, error_locator);

    // Compute derivative of error locator
    auto error_locator_derivative = error_locator.Derivative();

    // Apply Forney's formula
    for (size_t i = 0; i < error_locations.size(); ++i) {
      size_t loc = error_locations[i];
      ElementType alpha_inv = field->Inv(eval_points[loc]);

      ElementType numerator = error_evaluator.Evaluate(alpha_inv);
      ElementType denominator = error_locator_derivative.Evaluate(alpha_inv);

      if (denominator == ElementType{}) {
        throw std::runtime_error(
            "Error in Forney's algorithm: zero denominator");
      }

      ElementType error_value = field->Div(numerator, denominator);
      error_value = field->Mul(error_value, field->Neg(eval_points[loc]));
      error_value = field->Div(error_value, multipliers[loc]);

      error_values[i] = error_value;
    }

    return error_values;
  }

  // Compute error evaluator polynomial
  polynomial_type ComputeErrorEvaluator(
      const std::vector<ElementType>& syndrome,
      const polynomial_type& error_locator) const {
    auto field = grs_code_->Field();

    // Create syndrome polynomial
    std::vector<ElementType> syndrome_coeffs(syndrome.size() + 1,
                                             ElementType{});
    for (size_t i = 0; i < syndrome.size(); ++i) {
      syndrome_coeffs[i + 1] = syndrome[i];
    }
    polynomial_type syndrome_poly(syndrome_coeffs, field);

    // Error evaluator = (syndrome_poly * error_locator) mod x^(2t)
    auto product = syndrome_poly.Mul(error_locator);

    // Take only first 2t coefficients
    size_t max_degree = 2 * max_errors_;
    std::vector<ElementType> evaluator_coeffs(max_degree, ElementType{});
    for (size_t i = 0; i < max_degree && i <= product.Degree(); ++i) {
      evaluator_coeffs[i] = product.GetCoefficient(i);
    }

    return polynomial_type(evaluator_coeffs, field);
  }

  // Correct errors in received word
  codeword_type CorrectErrors(
      const codeword_type& received_word,
      const std::vector<size_t>& error_locations,
      const std::vector<ElementType>& error_values) const {
    auto field = grs_code_->Field();
    codeword_type corrected = received_word;

    for (size_t i = 0; i < error_locations.size(); ++i) {
      size_t loc = error_locations[i];
      corrected[loc] = field->Sub(corrected[loc], error_values[i]);
    }

    return corrected;
  }

  // Solve linear system Ax = b using Gaussian elimination
  std::vector<ElementType> SolveLinearSystem(
      const std::vector<std::vector<ElementType>>& A,
      const std::vector<ElementType>& b) const {
    auto field = grs_code_->Field();
    size_t n = A.size();

    if (n == 0 || A[0].size() != n || b.size() != n) {
      return {};  // Invalid dimensions
    }

    // Create augmented matrix
    std::vector<std::vector<ElementType>> augmented(
        n, std::vector<ElementType>(n + 1));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        augmented[i][j] = A[i][j];
      }
      augmented[i][n] = b[i];
    }

    // Forward elimination
    for (size_t i = 0; i < n; ++i) {
      // Find pivot
      size_t pivot = i;
      for (size_t j = i + 1; j < n; ++j) {
        if (augmented[j][i] != ElementType{}) {
          pivot = j;
          break;
        }
      }

      if (augmented[pivot][i] == ElementType{}) {
        return {};  // No solution
      }

      // Swap rows
      if (pivot != i) {
        std::swap(augmented[i], augmented[pivot]);
      }

      // Eliminate column
      ElementType pivot_inv = field->Inv(augmented[i][i]);
      for (size_t j = i + 1; j < n; ++j) {
        if (augmented[j][i] != ElementType{}) {
          ElementType factor = field->Mul(augmented[j][i], pivot_inv);
          for (size_t k = i; k <= n; ++k) {
            ElementType term = field->Mul(factor, augmented[i][k]);
            augmented[j][k] = field->Sub(augmented[j][k], term);
          }
        }
      }
    }

    // Back substitution
    std::vector<ElementType> solution(n);
    for (int i = n - 1; i >= 0; --i) {
      solution[i] = augmented[i][n];
      for (size_t j = i + 1; j < n; ++j) {
        ElementType term = field->Mul(augmented[i][j], solution[j]);
        solution[i] = field->Sub(solution[i], term);
      }
      solution[i] = field->Div(solution[i], augmented[i][i]);
    }

    return solution;
  }
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_GRS_DECODER_HPP
