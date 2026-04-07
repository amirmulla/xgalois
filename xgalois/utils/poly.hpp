#ifndef XG_UTILS_POLY_HPP_
#define XG_UTILS_POLY_HPP_

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace utils {

template <typename GaloisField>
std::pair<PolynomialDense<GaloisField>,
          std::pair<PolynomialDense<GaloisField>, PolynomialDense<GaloisField>>>
PolynomialDenseExtendedGcd(const PolynomialDense<GaloisField> &a,
                           const PolynomialDense<GaloisField> &b) {
  using ElementType = typename PolynomialDense<GaloisField>::ElementType;
  using PolynomialType = PolynomialDense<GaloisField>;

  std::shared_ptr<GaloisField> field = a.Field();
  if (!field) {
    if (b.Field()) {
      field = b.Field();
    } else {
      throw std::invalid_argument("Cannot determine field for GCD.");
    }
  }

  const std::string var_name = a.GetVariable();

  ElementType kGfZeroElem(field->AdditiveIdentity(), field);
  ElementType kGfOneElem(field->MultiplicativeIdentity(), field);
  PolynomialType kGfZeroPoly(std::vector<ElementType>{kGfZeroElem}, var_name);
  PolynomialType kGfOnePoly(std::vector<ElementType>{kGfOneElem}, var_name);

  if (b.Degree() == -1) {

    PolynomialType monic_a = a;
    if (a.Degree() != -1 && a[a.Degree()] != kGfOneElem) {
      ElementType inv_lc = a[a.Degree()].Inv();
      monic_a = a * inv_lc;
      return {std::move(monic_a),
              {PolynomialType(std::vector<ElementType>{inv_lc}, var_name),
               kGfZeroPoly}};
    }
    return {std::move(monic_a), {kGfOnePoly, kGfZeroPoly}};
  }

  if (a.Degree() == -1) {

    PolynomialType monic_b = b;
    if (b.Degree() != -1 && b[b.Degree()] != kGfOneElem) {
      ElementType inv_lc = b[b.Degree()].Inv();
      monic_b = b * inv_lc;
      return {std::move(monic_b),
              {kGfZeroPoly,
               PolynomialType(std::vector<ElementType>{inv_lc},
                              b.GetVariable())}};
    }
    return {
        std::move(monic_b),
        {kGfZeroPoly,
         kGfOnePoly}};

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

  while (r1.Degree() != -1) {
    auto division_result = r0.DivRem(r1);
    PolynomialType q = std::move(division_result.first);
    PolynomialType r_next = std::move(division_result.second);

    r0 = std::move(r1);
    r1 = std::move(r_next);

    PolynomialType s_next = s0 - q * s1;
    s0 = std::move(s1);
    s1 = std::move(s_next);

    PolynomialType t_next = t0 - q * t1;
    t0 = std::move(t1);
    t1 = std::move(t_next);
  }

  if (r0.Degree() != -1) {
    const ElementType &leading_coeff = r0[r0.Degree()];
    if (leading_coeff != kGfOneElem) {
      ElementType inv_leading_coeff = leading_coeff.Inv();

      r0 *= inv_leading_coeff;
      s0 *= inv_leading_coeff;
      t0 *= inv_leading_coeff;
    }
  } else {

    s0 = kGfZeroPoly;
    t0 = kGfZeroPoly;
  }

  return {std::move(r0), {std::move(s0), std::move(t0)}};
}

template <typename GaloisField>
static typename GaloisField::element_type ExtractCoefficient(
    const std::string &term) {
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

static int ExtractDegree(const std::string &term) {

  size_t x_pos = term.find('x');
  if (x_pos == std::string::npos) {
    x_pos = term.find('X');
    if (x_pos == std::string::npos) {
      return 0;
    }
  }

  size_t caret_pos = term.find('^', x_pos);
  if (caret_pos == std::string::npos) {

    return 1;
  }

  std::string power_str = term.substr(caret_pos + 1);

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

template <typename GaloisField>
PolynomialDense<GaloisField> ParsePolynomial(
    const std::shared_ptr<GaloisField> &base_field,
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

template <typename GaloisField>
PolynomialDense<GaloisField> ParsePolynomial(
    const std::shared_ptr<GaloisField> &base_field, const std::string &poly_str,
    const std::string &variable_name) {

  if (variable_name == "x") {
    return ParsePolynomial(base_field, poly_str);
  }

  std::string normalized_str = poly_str;
  if (!variable_name.empty() && variable_name != "x") {
    size_t pos = 0;
    while ((pos = normalized_str.find(variable_name, pos)) !=
           std::string::npos) {
      normalized_str.replace(pos, variable_name.length(), "x");
      pos += 1;
    }
  }

  auto result = ParsePolynomial(base_field, normalized_str);
  result.SetVariable(variable_name);
  return result;
}

template <typename GaloisField>
bool IsIrreducible(const PolynomialDense<GaloisField> &poly) {
  using ElementType = typename PolynomialDense<GaloisField>::ElementType;
  using PolynomialType = PolynomialDense<GaloisField>;

  std::shared_ptr<GaloisField> field = poly.Field();
  if (!field) {
    throw std::invalid_argument(
        "Cannot determine field for irreducibility test.");
  }

  int n = poly.Degree();

  if (n <= 0) {
    return false;
  }

  if (n == 1) {
    return true;
  }

  const std::string var_name = poly.GetVariable();
  ElementType kGfZeroElem(field->AdditiveIdentity(), field);
  ElementType kGfOneElem(field->MultiplicativeIdentity(), field);

  if (poly[n] == kGfZeroElem) {
    return false;
  }

  uint64_t q = field->Characteristic();

  std::vector<ElementType> x_coeffs = {kGfZeroElem, kGfOneElem};
  PolynomialType x_poly(x_coeffs, var_name);

  if (n == 2) {

    for (uint64_t i = 0; i < q; ++i) {
      ElementType test_val(i, field);
      std::vector<ElementType> test_coeffs = {test_val};
      PolynomialType test_poly(test_coeffs, var_name);

      ElementType result = kGfZeroElem;
      ElementType power_of_x = kGfOneElem;

      for (int j = 0; j <= n; ++j) {
        result += poly[j] * power_of_x;
        if (j < n) {
          power_of_x *= test_val;
        }
      }

      if (result == kGfZeroElem) {
        return false;
      }
    }
    return true;
  }

  PolynomialType current_poly = x_poly;

  for (int i = 1; i <= n / 2; ++i) {

    for (uint64_t j = 0; j < q; ++j) {

      PolynomialType temp = current_poly * current_poly;
      current_poly = temp.DivRem(poly).second;
    }

    PolynomialType diff = current_poly - x_poly;

    auto [gcd_result, factors] = PolynomialDenseExtendedGcd(poly, diff);

    if (gcd_result.Degree() > 0) {
      return false;
    }
  }

  return true;
}

template <typename GaloisField>
uint64_t BinaryPolynomialToUint(const PolynomialDense<GaloisField> &poly) {

  if (poly.Degree() >= 64) {
    throw std::invalid_argument("Polynomial degree (" +
                                std::to_string(poly.Degree()) +
                                ") must be less than maximum degree (64)");
  }

  uint64_t result = 0;
  for (int i = 0; i <= poly.Degree(); ++i) {
    const auto &coeff = poly[i];
    if (coeff.Value() == 1) {
      result |= (static_cast<uint64_t>(1) << i);
    } else if (coeff.Value() != 0) {
      throw std::invalid_argument(
          "Polynomial coefficients must be 0 or 1 for binary fields");
    }
  }

  return result;
}

}
}

#endif
