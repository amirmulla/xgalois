#ifndef XG_POLYNOMIAL_DENSE_HPP_
#define XG_POLYNOMIAL_DENSE_HPP_

// C system headers
#include <cassert>
#include <cmath>

// C++ standard library headers
#include <algorithm>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Project headers
#include "xgalois/field/gf_element.hpp"

namespace xg {

//------------------------------------------------------------------------------
// PolynomialDense Class
//------------------------------------------------------------------------------

// Represents a polynomial with coefficients from a Galois field,
// stored using a dense vector.
template <typename GaloisField>
class PolynomialDense {
 public:
  using ElementType = GaloisFieldElementBase<GaloisField>;

  // Constructs a polynomial from a vector of coefficients.
  // Coefficients are ordered from lowest degree (index 0) to highest.
  explicit PolynomialDense(const std::vector<ElementType> &coeffs,
                           const std::string &var = "x")
      : coefficients_(coeffs),
        field_(coeffs[0].Field()),
        kGfZero(field_->AdditiveIdentity(), field_),
        kGfOne(field_->MultiplicativeIdentity(), field_),
        variable_(var) {  // Initialize variable_ with var
    assert(!coeffs.empty() && "Coefficient vector cannot be empty.");
    if (!field_) {
      throw std::invalid_argument(
          "Cannot determine Galois field from coefficients.");
    }
    Trim();  // Ensures polynomial is in canonical form.
  }

  // Default copy constructor.
  PolynomialDense(const PolynomialDense &other) = default;
  // Default copy assignment operator.
  PolynomialDense &operator=(const PolynomialDense &other) = default;
  // Default move constructor.
  PolynomialDense(PolynomialDense &&other) = default;
  // Default move assignment operator.
  PolynomialDense &operator=(PolynomialDense &&other) = default;
  // Default destructor.
  ~PolynomialDense() = default;

  // Sets the variable string for polynomial representation.
  void SetVariable(const std::string &var) {
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    variable_ = var;
  }

  // Returns the variable string used for polynomial representation.
  std::string GetVariable() const { return variable_; }

  // Returns the degree of the polynomial.
  // Returns -1 for the zero polynomial.
  int Degree() const {
    if (coefficients_.size() == 1 && coefficients_[0] == kGfZero) {
      return -1;  // Zero polynomial case.
    }
    return coefficients_.size() - 1;
  }

  // Returns the number of coefficients stored (degree + 1 for non-zero, 1 for
  // zero).
  size_t Size() const { return coefficients_.size(); }

  // Accesses the coefficient at the given index (read-only).
  const ElementType &operator[](int index) const {
    if (index < 0 || index >= coefficients_.size()) {
      throw std::out_of_range("Coefficient index out of bounds.");
    }
    return coefficients_[index];
  }

  // Accesses the coefficient at the given index (read/write).
  // Note: Modifying coefficients may require a subsequent Trim() call
  // if the leading coefficient becomes zero.
  ElementType &operator[](int index) {
    if (index < 0 || index >= coefficients_.size()) {
      throw std::out_of_range("Coefficient index out of bounds.");
    }
    return coefficients_[index];
    // Note: Modifying coefficients directly might require re-normalization
    // if the highest coefficient becomes zero. This is a potential area
    // for refinement depending on intended usage.
  }

  // Adds two polynomials.
  PolynomialDense<GaloisField> operator+(
      const PolynomialDense<GaloisField> &other) const {
    size_t result_size = std::max(Size(), other.Size());
    std::vector<ElementType> result_coeffs(result_size, kGfZero);

    size_t min_size = std::min(Size(), other.Size());
    for (size_t i = 0; i < min_size; ++i) {
      result_coeffs[i] = coefficients_[i] + other[i];
    }

    if (Size() > other.Size()) {
      for (size_t i = min_size; i < Size(); ++i) {
        result_coeffs[i] = coefficients_[i];
      }
    } else if (other.Size() > Size()) {
      for (size_t i = min_size; i < other.Size(); ++i) {
        result_coeffs[i] = other[i];
      }
    }

    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();  // Ensure canonical form.
    return result;
  }

  // Subtracts one polynomial from another.
  PolynomialDense<GaloisField> operator-(
      const PolynomialDense<GaloisField> &other) const {
    size_t result_size = std::max(Size(), other.Size());
    std::vector<ElementType> result_coeffs(result_size, kGfZero);

    size_t min_size = std::min(Size(), other.Size());
    for (size_t i = 0; i < min_size; ++i) {
      result_coeffs[i] = coefficients_[i] - other[i];
    }

    if (Size() > other.Size()) {
      for (size_t i = min_size; i < Size(); ++i) {
        result_coeffs[i] = coefficients_[i];
      }
    } else if (other.Size() > Size()) {
      for (size_t i = min_size; i < other.Size(); ++i) {
        result_coeffs[i] = kGfZero - other[i];
      }
    }

    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();  // Ensure canonical form.
    return result;
  }

  // Multiplies two polynomials (convolution).
  PolynomialDense<GaloisField> operator*(
      const PolynomialDense<GaloisField> &other) const {
    if ((Size() == 1 && coefficients_[0] == kGfZero) ||
        (other.Size() == 1 && other[0] == kGfZero)) {
      return PolynomialDense<GaloisField>(std::vector<ElementType>{kGfZero},
                                          variable_);
    }

    size_t result_size = (Size() - 1) + (other.Size() - 1) + 1;
    std::vector<ElementType> result_coeffs(result_size, kGfZero);

    for (size_t i = 0; i < Size(); ++i) {
      for (size_t j = 0; j < other.Size(); ++j) {
        result_coeffs[i + j] += coefficients_[i] * other[j];
      }
    }

    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();  // Ensure canonical form.
    return result;
  }

  // Divides this polynomial by another, returning the quotient.
  PolynomialDense<GaloisField> operator/(
      const PolynomialDense<GaloisField> &other) const {
    auto result_pair = DivRem(other);
    return std::move(result_pair.first);
  }

  // Computes the remainder of this polynomial divided by another.
  PolynomialDense<GaloisField> operator%(
      const PolynomialDense<GaloisField> &other) const {
    auto result_pair = DivRem(other);
    return std::move(result_pair.second);
  }

  // Multiplies this polynomial by a scalar field element.
  PolynomialDense<GaloisField> operator*(const ElementType &scalar) const {
    std::vector<ElementType> result_coeffs(coefficients_.size());
    for (size_t i = 0; i < coefficients_.size(); ++i) {
      result_coeffs[i] = coefficients_[i] * scalar;
    }
    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();  // Ensure canonical form.
    return result;
  }

  // Negates this polynomial.
  PolynomialDense<GaloisField> operator-() const {
    std::vector<ElementType> result_coeffs(coefficients_.size());
    for (size_t i = 0; i < coefficients_.size(); ++i) {
      result_coeffs[i] = kGfZero - coefficients_[i];
    }
    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();  // Ensure canonical form.
    return result;
  }

  // Raises this polynomial to the given power using exponentiation by squaring.
  PolynomialDense<GaloisField> operator^(uint64_t exponent) const {
    if (exponent == 0) {
      // p(x)^0 = 1 (the constant polynomial).
      return PolynomialDense<GaloisField>(std::vector<ElementType>{kGfOne},
                                          variable_);
    }
    if (exponent == 1) {
      return *this;
    }

    PolynomialDense<GaloisField> base = *this;
    PolynomialDense<GaloisField> result(std::vector<ElementType>{kGfOne},
                                        variable_);

    uint64_t temp_exponent = exponent;

    while (temp_exponent > 0) {
      if (temp_exponent % 2 == 1) {
        result = result * base;  // Result is trimmed in operator*
      }
      base = base * base;  // Base is trimmed in operator*
      temp_exponent /= 2;
    }
    // Result is already trimmed by the final multiplication.
    return result;
  }

  // Compound assignment for addition.
  PolynomialDense<GaloisField> &operator+=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this + other;
    return *this;
  }

  // Compound assignment for subtraction.
  PolynomialDense<GaloisField> &operator-=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this - other;
    return *this;
  }

  // Compound assignment for multiplication.
  PolynomialDense<GaloisField> &operator*=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this * other;
    return *this;
  }

  // Compound assignment for scalar multiplication.
  PolynomialDense<GaloisField> &operator*=(const ElementType &scalar) {
    *this = *this * scalar;
    return *this;
  }

  // Compound assignment for division.
  PolynomialDense<GaloisField> &operator/=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this / other;
    return *this;
  }

  // Compound assignment for modulo.
  PolynomialDense<GaloisField> &operator%=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this % other;
    return *this;
  }

  // Checks for equality between two polynomials.
  // Relies on Trim() to ensure canonical representation for comparison.
  bool operator==(const PolynomialDense<GaloisField> &other) const {
    // Assumes both this and other are already trimmed or will be
    // implicitly by their operations. Direct coefficient comparison is fine.
    return coefficients_ == other.coefficients_;
  }

  // Checks for inequality between two polynomials.
  bool operator!=(const PolynomialDense<GaloisField> &other) const {
    return !(*this == other);
  }

  // Evaluates the polynomial at a given field element `x` using Horner's
  // method.
  ElementType operator()(const ElementType &x) const {
    ElementType result = kGfZero;
    for (int i = coefficients_.size() - 1; i >= 0; --i) {
      result = result * x + coefficients_[i];
    }
    return result;
  }

  // Returns a shared pointer to the Galois field of the coefficients.
  std::shared_ptr<GaloisField> Field() const { return field_; }

  // Prints the polynomial to an output stream.
  void Print(std::ostream &os) const {
    if (coefficients_.empty() ||
        (coefficients_.size() == 1 && coefficients_[0] == kGfZero)) {
      os << "0";  // Zero polynomial.
      return;
    }

    // Helper method to check if a coefficient needs brackets
    auto needs_brackets = [](const ElementType &coeff) {
      std::stringstream ss;
      ss << coeff;
      std::string coeff_str = ss.str();

      for (char c : coeff_str) {
        if (!std::isdigit(c)) {
          return true;
        }
      }
      return false;
    };

    bool first_term = true;
    for (int i = static_cast<int>(coefficients_.size()) - 1; i >= 0; --i) {
      const ElementType &coeff = coefficients_[i];

      if (coeff == kGfZero && coefficients_.size() > 1) {
        continue;  // Skip zero terms unless it's the zero polynomial itself.
      }

      if (!first_term) {
        os << " + ";
      }

      // Print coefficient if it's not 1 or if it's the constant term.
      if (i == 0 || coeff != kGfOne ||
          (coeff == kGfOne && i == 0 && coefficients_.size() == 1)) {
        if (needs_brackets(coeff)) {
          os << "(" << coeff << ")";
        } else {
          os << coeff;
        }
      } else if (coeff == kGfOne && i > 0 && first_term &&
                 coefficients_.size() - 1 == i) {
        // Special case for leading term x^n (coeff is 1)
      } else if (coeff != kGfOne) {
        os << coeff;
      }

      if (i > 0) {
        os << variable_;  // Use the variable_ member
        if (i > 1) {
          os << "^" << i;
        }
      }
      first_term = false;
    }
    if (first_term) {  // Handles case where all coefficients were zero but not
                       // caught by initial check (e.g. [0,0,0] after some ops)
      os << "0";
    }
  }

  // Performs polynomial long division, returning a pair {quotient, remainder}.
  std::pair<PolynomialDense<GaloisField>, PolynomialDense<GaloisField>> DivRem(
      const PolynomialDense<GaloisField> &divisor) const {
    if (divisor.Size() == 1 && divisor[0] == kGfZero) {
      throw std::invalid_argument("Division by zero polynomial.");
    }

    PolynomialDense<GaloisField> quotient(std::vector<ElementType>{kGfZero},
                                          variable_);
    PolynomialDense<GaloisField> remainder =
        *this;  // Inherits variable via copy

    int divisor_degree = divisor.Degree();

    // If dividend degree is less than divisor degree, quotient is 0, remainder
    // is dividend.
    if (remainder.Degree() < divisor_degree) {
      // Constructor already handles trimming, no need to call Trim() explicitly
      // here
      return {std::move(quotient), std::move(remainder)};
    }

    while (remainder.Degree() >= divisor_degree && remainder.Degree() != -1) {
      int remainder_degree = remainder.Degree();
      int degree_difference = remainder_degree - divisor_degree;

      const ElementType &remainder_leading_coeff = remainder[remainder_degree];
      const ElementType &divisor_leading_coeff = divisor[divisor_degree];

      if (divisor_leading_coeff == kGfZero) {
        throw std::runtime_error(
            "Internal error: Divisor leading coefficient "
            "is zero during division.");
      }

      ElementType term_coeff = remainder_leading_coeff / divisor_leading_coeff;

      std::vector<ElementType> term_coeffs(degree_difference + 1, kGfZero);
      term_coeffs[degree_difference] = term_coeff;
      // Constructor already calls Trim()
      PolynomialDense<GaloisField> term_poly(term_coeffs, variable_);

      // += operator already handles trimming
      quotient += term_poly;

      // * and -= operators already handle trimming
      PolynomialDense<GaloisField> product = term_poly * divisor;
      remainder -= product;

      // This Trim is needed because remainder might have lost its highest
      // degree term
      remainder.Trim();
    }

    // No need for explicit trimming here as operations already handle it
    return {std::move(quotient), std::move(remainder)};
  }

  // Computes the formal derivative of the polynomial.
  // For polynomial f(x) = a_n*x^n + ... + a_1*x + a_0,
  // the derivative is f'(x) = n*a_n*x^(n-1) + ... + 1*a_1
  PolynomialDense<GaloisField> Derivative() const {
    if (coefficients_.size() <= 1) {
      // Derivative of constant polynomial is zero
      return PolynomialDense<GaloisField>(std::vector<ElementType>{kGfZero},
                                          variable_);
    }

    std::vector<ElementType> derivative_coeffs;
    derivative_coeffs.reserve(coefficients_.size() - 1);

    // For each term a_i * x^i, derivative is i * a_i * x^(i-1)
    for (size_t i = 1; i < coefficients_.size(); ++i) {
      // Multiply coefficient by the power (i)
      ElementType coeff = coefficients_[i];

      // In finite fields, we need to compute i * coeff using repeated addition
      // This avoids the need for ElementFromInt which doesn't exist for all
      // field types
      // TODO(amirmulla): more efficient way to compute i * coeff
      ElementType derivative_coeff = kGfZero;
      for (size_t j = 0; j < i; ++j) {
        derivative_coeff += coeff;
      }

      derivative_coeffs.push_back(derivative_coeff);
    }

    // Handle case where all coefficients resulted in zero derivative
    if (derivative_coeffs.empty()) {
      derivative_coeffs.push_back(kGfZero);
    }

    PolynomialDense<GaloisField> result(derivative_coeffs, variable_);
    result.Trim();  // Ensure canonical form
    return result;
  }
  // Removes trailing zero coefficients to maintain a canonical representation.
  // Ensures the zero polynomial is represented as {0}.
  void Trim() {
    int last_nonzero = coefficients_.size() - 1;
    while (last_nonzero >= 0 && coefficients_[last_nonzero] == kGfZero) {
      last_nonzero--;
    }

    if (last_nonzero < 0) {
      coefficients_ = {kGfZero};
    } else {
      coefficients_.resize(last_nonzero + 1);
    }
  }

 private:
  // Coefficients of the polynomial, coefficients_[i] is for x^i.
  std::vector<ElementType> coefficients_;
  // Shared pointer to the Galois field for coefficient arithmetic.
  std::shared_ptr<GaloisField> field_;
  // Cached additive identity (0) of the field.
  ElementType kGfZero;
  // Cached multiplicative identity (1) of the field.
  ElementType kGfOne;
  // String representation of the polynomial variable.
  std::string variable_;
};

/// Overloads the << operator for printing PolynomialDense objects.
template <typename GaloisField>
std::ostream &operator<<(std::ostream &os,
                         const PolynomialDense<GaloisField> &poly) {
  poly.Print(os);
  return os;
}

}  // namespace xg

#endif  // XG_POLYNOMIAL_DENSE_HPP_
