#pragma once

#include <stdexcept>
#include <string>

#include "xgalois/field/gf_base.hpp"

namespace xg {
namespace utils {

// Convert string representation to enum
FieldRepresentation ConvertRepresentation(const std::string &rep) {
  if (rep == "int") return FieldRepresentation::INT;
  if (rep == "hex") return FieldRepresentation::HEX;
  if (rep == "pow") return FieldRepresentation::POW;
  if (rep == "log") return FieldRepresentation::LOG;
  if (rep == "poly") return FieldRepresentation::POLY;
  throw std::invalid_argument("Unknown representation: " + rep);
}

// Parse power string like "g^5" or "g^-3" and return the field element
template <typename ElementType>
ElementType ParsePowerString(const std::string &pow_str,
                             const GaloisFieldBase<ElementType> &field) {
  if (pow_str.find('^') == std::string::npos) {
    throw std::invalid_argument("Invalid power format");
  }

  size_t pos = pow_str.find('^');
  std::string exp_str = pow_str.substr(pos + 1);

  // Handle negative exponents
  bool is_negative = false;
  if (!exp_str.empty() && exp_str[0] == '-') {
    is_negative = true;
    exp_str = exp_str.substr(1);  // Remove the minus sign
  }

  if (exp_str.empty()) {
    throw std::invalid_argument("Invalid exponent format");
  }

  uint32_t exp = std::stoul(exp_str);

  // Get the multiplicative group order (field order - 1)
  uint32_t group_order = field.Order() - 1;

  // Handle negative exponents: g^(-n) = g^(group_order - n)
  if (is_negative) {
    exp = exp % group_order;
    if (exp == 0) {
      exp = 0;  // g^0 = 1, so g^(-0) = g^0 = 1
    } else {
      exp = group_order - exp;
    }
  } else {
    exp = exp % group_order;
  }

  ElementType generator = field.MultiplicativeGenerator();
  return field.Pow(generator, exp);
}

}  // namespace utils
}  // namespace xg