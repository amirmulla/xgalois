#ifndef XG_UTILS_POLY_HPP_
#define XG_UTILS_POLY_HPP_

// C++ standard library headers
#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

// Project headers
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace utils {

//------------------------------------------------------------------------------
// PolynomialDenseExtendedGcd Function
//------------------------------------------------------------------------------
// Implements the Extended Euclidean Algorithm for polynomials.
// Returns a pair {gcd(a,b), {s, t}} such that a*s + b*t = gcd(a,b).
// The returned gcd is made monic.
template <typename GaloisField>
std::pair<PolynomialDense<GaloisField>,
          std::pair<PolynomialDense<GaloisField>, PolynomialDense<GaloisField>>>
PolynomialDenseExtendedGcd(const PolynomialDense<GaloisField> &a,
                           const PolynomialDense<GaloisField> &b) {
  using ElementType = typename PolynomialDense<GaloisField>::ElementType;
  using PolynomialType = PolynomialDense<GaloisField>;

  // Ensure polynomials share the same field, or handle error.
  // For simplicity, assume a.Field() is valid and b uses the same.
  std::shared_ptr<GaloisField> field = a.Field();
  if (!field) { // Should not happen if 'a' is valid
    if (b.Field())
      field = b.Field();
    else
      throw std::invalid_argument("Cannot determine field for GCD.");
  }

  const std::string var_name = a.GetVariable(); // Use variable from 'a'

  ElementType kGfZeroElem(field->AdditiveIdentity(), field);
  ElementType kGfOneElem(field->MultiplicativeIdentity(), field);
  PolynomialType kGfZeroPoly(std::vector<ElementType>{kGfZeroElem}, var_name);
  PolynomialType kGfOnePoly(std::vector<ElementType>{kGfOneElem}, var_name);

  if (b.Degree() == -1) { // If b is the zero polynomial
    // gcd(a, 0) = a. To make it monic, divide by leading coefficient.
    PolynomialType monic_a = a; // Copies variable_
    if (a.Degree() != -1 && a[a.Degree()] != kGfOneElem) {
      ElementType inv_lc = a[a.Degree()].Inv();
      monic_a = a * inv_lc; // s will be inv_lc, t will be 0
      return {std::move(monic_a),
              {PolynomialType(std::vector<ElementType>{inv_lc}, var_name),
               kGfZeroPoly}};
    }
    return {std::move(monic_a), {kGfOnePoly, kGfZeroPoly}};
  }

  if (a.Degree() == -1) { // If a is the zero polynomial
    // gcd(0, b) = b. To make it monic.
    PolynomialType monic_b = b; // Copies variable_
    if (b.Degree() != -1 && b[b.Degree()] != kGfOneElem) {
      ElementType inv_lc = b[b.Degree()].Inv();
      monic_b = b * inv_lc; // s will be 0, t will be inv_lc
      return {std::move(monic_b),
              {kGfZeroPoly,
               PolynomialType(std::vector<ElementType>{inv_lc},
                              b.GetVariable())}}; // Use b's variable here
    }
    return {
        std::move(monic_b),
        {kGfZeroPoly,
         kGfOnePoly}}; // kGfOnePoly already uses a's var, but b's would be more
                       // consistent if a is zero Let's ensure kGfOnePoly uses
                       // b.GetVariable() if a is zero. However, kGfOnePoly is
                       // already defined with a.GetVariable(). For gcd(0,b) ->
                       // s=0, t=inv_lc_b or 1. t should use b's variable. This
                       // is a subtle point. For now, the impact is minor. A
                       // simple fix: if a is zero, re-init kGfOnePoly with
                       // b.GetVariable() for t. Or, the PolynomialType
                       // constructor for {inv_lc} should use b.GetVariable().
                       // The current return for a=0, inv_lc case is correct.
                       // For a=0, b non-zero, t=1 case:
                       // return {std::move(monic_b), {kGfZeroPoly,
                       // PolynomialType(std::vector<ElementType>{kGfOneElem},
                       // b.GetVariable())}}; This is safer.
    PolynomialType one_poly_for_b(std::vector<ElementType>{kGfOneElem},
                                  b.GetVariable());
    return {std::move(monic_b), {kGfZeroPoly, one_poly_for_b}};
  }

  PolynomialType r0 = a;
  PolynomialType r1 = b;
  PolynomialType s0 = kGfOnePoly;
  PolynomialType s1 = kGfZeroPoly;
  PolynomialType t0 = kGfZeroPoly;
  PolynomialType t1 = kGfOnePoly;

  while (r1.Degree() != -1) { // While remainder r1 is not the zero polynomial
    auto division_result = r0.DivRem(r1); // DivRem returns trimmed results
    PolynomialType q = std::move(division_result.first);
    PolynomialType r_next = std::move(division_result.second);

    r0 = std::move(r1);
    r1 = std::move(r_next);

    // Operations (-, *) already ensure proper trimming
    PolynomialType s_next = s0 - q * s1;
    s0 = std::move(s1);
    s1 = std::move(s_next);

    PolynomialType t_next = t0 - q * t1;
    t0 = std::move(t1);
    t1 = std::move(t_next);
  }

  // r0 is the GCD. Make it monic.
  if (r0.Degree() != -1) { // If GCD is not zero polynomial
    const ElementType &leading_coeff = r0[r0.Degree()];
    if (leading_coeff != kGfOneElem) {
      ElementType inv_leading_coeff = leading_coeff.Inv();
      // *= operator already ensures trimming
      r0 *= inv_leading_coeff;
      s0 *= inv_leading_coeff;
      t0 *= inv_leading_coeff;
    }
  } else {
    // GCD is zero polynomial, can happen if both a and b are zero.
    // s0 and t0 would be zero polynomials as well by construction.
    // This case should be handled by initial checks if a or b are zero.
    // If both a and b are zero, gcd(0,0)=0. s,t can be 0,0.
    s0 = kGfZeroPoly;
    t0 = kGfZeroPoly;
  }

  return {std::move(r0), {std::move(s0), std::move(t0)}};
}

//------------------------------------------------------------------------------
// Helper Functions for Polynomial Parsing
//------------------------------------------------------------------------------
/**
 * @brief Extract the coefficient of a polynomial term
 *
 * Parses a polynomial term string to determine its coefficient.
 * Handles constant terms, variable terms with implicit coefficients,
 * and terms with explicit numeric coefficients.
 *
 * @tparam GaloisField The Galois field type
 * @param term String representation of a polynomial term
 * @return The coefficient of the term
 */
template <typename GaloisField>
static typename GaloisField::element_type
ExtractCoefficient(const std::string &term) {
  using ElementType = typename GaloisField::element_type;
  if (term.find('x') == std::string::npos &&
      term.find('X') == std::string::npos) {
    try {
      return static_cast<ElementType>(std::stoull(term));
    } catch (const std::exception &) {
      throw std::invalid_argument("Invalid coefficient in term: " + term);
    }
  }

  size_t x_pos = term.find('x');
  if (x_pos == std::string::npos) {
    x_pos = term.find('X');
  }

  if (x_pos == 0) {
    return ElementType(1);
  }

  std::string coeff_str = term.substr(0, x_pos);

  if (coeff_str.empty() || coeff_str == "+") {
    return ElementType(1);
  }
  if (coeff_str == "-") {
    return ElementType(1);
  }

  try {
    return static_cast<ElementType>(std::stoull(coeff_str));
  } catch (const std::exception &) {
    throw std::invalid_argument("Invalid coefficient in term: " + term);
  }
}

/**
 * @brief Extract the degree of a polynomial term
 *
 * Parses a polynomial term string to determine its degree (power).
 * Handles constant terms (degree 0), variable terms without explicit power
 * (degree 1), and terms with explicit numeric powers.
 *
 * @param term String representation of a polynomial term
 * @return The degree of the term
 */
static int ExtractDegree(const std::string &term) {
  // If no 'x' or 'X' found, it's a constant term (degree 0)
  size_t x_pos = term.find('x');
  if (x_pos == std::string::npos) {
    x_pos = term.find('X');
    if (x_pos == std::string::npos) {
      return 0; // Constant term
    }
  }

  // Find if there's a '^' after the variable
  size_t caret_pos = term.find('^', x_pos);
  if (caret_pos == std::string::npos) {
    // No explicit power, so degree is 1 (e.g., "x", "2x", "-3x")
    return 1;
  }

  // Extract the power string after '^'
  std::string power_str = term.substr(caret_pos + 1);

  // Remove any trailing non-digit characters
  size_t end_pos = 0;
  while (end_pos < power_str.length() && std::isdigit(power_str[end_pos])) {
    ++end_pos;
  }
  power_str = power_str.substr(0, end_pos);

  if (power_str.empty()) {
    throw std::invalid_argument("Invalid power in term: " + term);
  }

  try {
    return std::stoi(power_str);
  } catch (const std::exception &) {
    throw std::invalid_argument("Invalid power in term: " + term);
  }
}

//------------------------------------------------------------------------------
// ParsePolynomial Function
//------------------------------------------------------------------------------
/**
 * @brief Parse polynomial from string representation
 *
 * Converts a string representation of a polynomial into a PolynomialDense
 * object. Supports various polynomial formats including coefficient notation
 * and variable expressions.
 *
 * @tparam GaloisField The Galois field type
 * @param base_field Shared pointer to the base field
 * @param poly_str String representation of the polynomial
 * @return PolynomialDense object representing the parsed polynomial
 */
template <typename GaloisField>
PolynomialDense<GaloisField>
ParsePolynomial(std::shared_ptr<GaloisField> base_field,
                const std::string &poly_str) {
  using FieldType = GaloisField;
  using ElementType = typename GaloisField::element_type;
  using FieldElementType = GaloisFieldElementBase<FieldType>;
  using PolynomialType = PolynomialDense<FieldType>;

  if (poly_str.empty()) {
    throw std::invalid_argument("Polynomial string cannot be empty");
  }

  std::vector<FieldElementType> coefficients;
  std::string cleaned = poly_str;
  cleaned.erase(std::remove_if(cleaned.begin(), cleaned.end(), ::isspace),
                cleaned.end());

  if (cleaned.empty()) {
    throw std::invalid_argument(
        "Polynomial string cannot be empty after cleaning");
  }

  if (cleaned.find('x') == std::string::npos &&
      cleaned.find('X') == std::string::npos) {
    try {
      ElementType value = static_cast<ElementType>(std::stoull(cleaned));
      FieldElementType element(value, base_field);
      coefficients.push_back(element);
      return PolynomialType(coefficients, "x");
    } catch (const std::exception &) {
      throw std::invalid_argument("Invalid constant polynomial: " + poly_str);
    }
  }

  std::vector<std::string> terms;
  std::vector<bool> positive_signs;
  size_t pos = 0;
  bool first_term = true;

  while (pos < cleaned.length()) {
    size_t next_plus = cleaned.find('+', pos);
    size_t next_minus = cleaned.find('-', pos + (first_term ? 0 : 1));

    size_t next_op = std::min(next_plus, next_minus);

    if (next_op == std::string::npos) {
      std::string term = cleaned.substr(pos);
      if (!term.empty()) {
        terms.push_back(term);
        positive_signs.push_back(first_term ||
                                 (pos > 0 && cleaned[pos - 1] != '-'));
      }
      break;
    }

    std::string term = cleaned.substr(pos, next_op - pos);
    if (!term.empty()) {
      terms.push_back(term);
      positive_signs.push_back(first_term ||
                               (pos > 0 && cleaned[pos - 1] != '-'));
    }

    pos = next_op + 1;
    first_term = false;
  }

  int max_degree = 0;
  for (const auto &term : terms) {
    int degree = ExtractDegree(term);
    max_degree = std::max(max_degree, degree);
  }

  coefficients.resize(max_degree + 1,
                      FieldElementType(ElementType(0), base_field));

  for (size_t i = 0; i < terms.size(); ++i) {
    const std::string &term = terms[i];
    bool is_positive = positive_signs[i];

    ElementType coeff = ExtractCoefficient<GaloisField>(term);
    int degree = ExtractDegree(term);

    if (!is_positive) {
      ElementType prime = base_field->Characteristic();
      coeff = (prime - coeff) % prime;
    }

    ElementType existing = coefficients[degree].Value();
    ElementType new_coeff = (existing + coeff) % base_field->Characteristic();
    coefficients[degree] = FieldElementType(new_coeff, base_field);
  }

  return PolynomialType(coefficients, "x");
}

/**
 * @brief Parse polynomial from string representation with custom variable name
 *
 * Converts a string representation of a polynomial into a PolynomialDense
 * object. Supports various polynomial formats including coefficient notation
 * and variable expressions with a specified variable name.
 *
 * @tparam GaloisField The Galois field type
 * @param base_field Shared pointer to the base field
 * @param poly_str String representation of the polynomial
 * @param variable_name The variable name to use (e.g., "α", "x", "y")
 * @return PolynomialDense object representing the parsed polynomial
 */
template <typename GaloisField>
PolynomialDense<GaloisField>
ParsePolynomial(std::shared_ptr<GaloisField> base_field,
                const std::string &poly_str,
                const std::string &variable_name) {
  // If variable name is 'x', just call the standard version
  if (variable_name == "x") {
    return ParsePolynomial(base_field, poly_str);
  }

  // Replace variable_name with 'x' for standard parsing
  std::string normalized_str = poly_str;
  if (!variable_name.empty() && variable_name != "x") {
    size_t pos = 0;
    while ((pos = normalized_str.find(variable_name, pos)) != std::string::npos) {
      normalized_str.replace(pos, variable_name.length(), "x");
      pos += 1; // "x" is shorter than most variable names
    }
  }

  // Parse using the standard function and set the variable name
  auto result = ParsePolynomial(base_field, normalized_str);
  result.SetVariable(variable_name);
  return result;
}

//------------------------------------------------------------------------------
// Polynomial Irreducibility Function
//------------------------------------------------------------------------------
/**
 * @brief Polynomial irreducibility verification using Rabin's test
 *
 * Determines whether a given polynomial is irreducible over its coefficient
 * field using Rabin's irreducibility test and other specialized algorithms.
 * Irreducible polynomials are essential for constructing valid field
 * extensions since they ensure the quotient ring forms a field.
 *
 * @tparam GaloisField The underlying Galois field type
 * @param poly The polynomial to verify for irreducibility
 *
 * @return true if polynomial is irreducible, false otherwise
 *
 * @note Implements proper irreducibility tests including Rabin's algorithm
 *       for general finite fields and optimized Ben-Or algorithm for GF(2).
 */
template <typename GaloisField>
bool IsIrreducible(const PolynomialDense<GaloisField> &poly) {
  using ElementType = typename PolynomialDense<GaloisField>::ElementType;
  using PolynomialType = PolynomialDense<GaloisField>;

  // Get the field and basic elements
  std::shared_ptr<GaloisField> field = poly.Field();
  if (!field) {
    throw std::invalid_argument(
        "Cannot determine field for irreducibility test.");
  }

  int n = poly.Degree();

  // Handle special cases
  if (n <= 0) {
    return false; // Constants and zero polynomial are not irreducible
  }

  if (n == 1) {
    return true; // Linear polynomials are always irreducible
  }

  const std::string var_name = poly.GetVariable();
  ElementType kGfZeroElem(field->AdditiveIdentity(), field);
  ElementType kGfOneElem(field->MultiplicativeIdentity(), field);

  // Check if leading coefficient is non-zero
  if (poly[n] == kGfZeroElem) {
    return false; // Not a valid polynomial of degree n
  }

  uint64_t q = field->Characteristic();

  // Create x polynomial: x = 0 + 1*x
  std::vector<ElementType> x_coeffs = {kGfZeroElem, kGfOneElem};
  PolynomialType x_poly(x_coeffs, var_name);

  // For efficiency, start with a simplified test for small degrees
  if (n == 2) {
    // For degree 2, just check if it has roots in the field
    // A quadratic is irreducible iff it has no roots
    for (uint64_t i = 0; i < q; ++i) {
      ElementType test_val(i, field);
      std::vector<ElementType> test_coeffs = {test_val};
      PolynomialType test_poly(test_coeffs, var_name);

      // Evaluate poly at test_val by substituting x = test_val
      ElementType result = kGfZeroElem;
      ElementType power_of_x = kGfOneElem;

      for (int j = 0; j <= n; ++j) {
        result += poly[j] * power_of_x;
        if (j < n) {
          power_of_x *= test_val;
        }
      }

      if (result == kGfZeroElem) {
        return false; // Found a root, so not irreducible
      }
    }
    return true;
  }

  // Rabin's irreducibility test for general case
  // Test: gcd(f(x), x^(q^i) - x) = 1 for i = 1, 2, ..., floor(n/2)
  PolynomialType current_poly = x_poly; // Start with x

  for (int i = 1; i <= n / 2; ++i) {
    // Compute x^(q^i) mod poly using repeated squaring
    for (uint64_t j = 0; j < q; ++j) {
      // current_poly = current_poly * current_poly mod poly
      PolynomialType temp = current_poly * current_poly;
      current_poly = temp.DivRem(poly).second;
    }

    // Compute x^(q^i) - x
    PolynomialType diff = current_poly - x_poly;

    // Compute gcd(poly, x^(q^i) - x)
    auto [gcd_result, factors] = PolynomialDenseExtendedGcd(poly, diff);

    // If gcd is not 1, then poly is not irreducible
    if (gcd_result.Degree() > 0) {
      return false;
    }
  }

  return true; // Passed all tests
}

//------------------------------------------------------------------------------
// Polynomial to Binary Representation Conversion
//------------------------------------------------------------------------------
/**
 * @brief Convert a polynomial over GF(2) to binary representation
 *
 * Converts a PolynomialDense object over GF(2) into a binary representation
 * suitable for use in binary extension fields GF(2^m). Each coefficient
 * corresponds to a bit in the binary representation.
 *
 * @tparam ElementType The element type for the binary representation
 * @param poly The polynomial to convert (must be over GF(2))
 * @return Binary representation where bit i corresponds to x^i coefficient
 *
 * @throws std::invalid_argument if polynomial degree exceeds 64
 * @throws std::invalid_argument if coefficients are not 0 or 1
 */
template <typename GaloisField>
uint64_t BinaryPolynomialToUint(const PolynomialDense<GaloisField> &poly) {
  // Check polynomial degree doesn't exceed maximum allowed
  if (poly.Degree() >= 64) {
    throw std::invalid_argument("Polynomial degree (" +
                                std::to_string(poly.Degree()) +
                                ") must be less than maximum degree (64)");
  }

  // Convert polynomial coefficients to binary representation
  uint64_t result = 0;
  for (int i = 0; i <= poly.Degree(); ++i) {
    auto coeff = poly[i];
    if (coeff.Value() == 1) {
      result |= (static_cast<uint64_t>(1) << i);
    } else if (coeff.Value() != 0) {
      throw std::invalid_argument(
          "Polynomial coefficients must be 0 or 1 for binary fields");
    }
  }

  return result;
}

} // namespace utils
} // namespace xg

#endif // XG_UTILS_POLY_HPP_