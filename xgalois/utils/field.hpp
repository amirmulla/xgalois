#ifndef XGALOIS_UTILS_FIELD_HPP_
#define XGALOIS_UTILS_FIELD_HPP_

#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "xgalois/field/gf_base.hpp"

namespace xg {
namespace utils {

inline FieldRepresentation ConvertRepresentation(const std::string& rep) {
  if (rep == "int") return FieldRepresentation::INT;
  if (rep == "hex") return FieldRepresentation::HEX;
  if (rep == "pow") return FieldRepresentation::POW;
  if (rep == "log") return FieldRepresentation::LOG;
  if (rep == "poly") return FieldRepresentation::POLY;
  throw std::invalid_argument("Unknown representation: " + rep);
}

template <typename ElementType>
inline ElementType ParsePowerString(
    const std::string& pow_str, const GaloisFieldBase<ElementType>& field,
    const std::string& expected_generator = "") {
  size_t pos = pow_str.find('^');
  if (pos == std::string::npos) {
    throw std::invalid_argument("Invalid power format");
  }

  std::string base = pow_str.substr(0, pos);
  if (base.empty()) {
    throw std::invalid_argument("Invalid generator in power format");
  }

  if (!expected_generator.empty()) {
    if (base != expected_generator) {
      throw std::invalid_argument("Unexpected generator in power format: " +
                                  base);
    }
  } else {
    auto is_valid_start = [](unsigned char c) {
      return std::isalpha(c) != 0 || c == '_' || c >= 0x80;
    };
    auto is_valid_cont = [](unsigned char c) {
      return std::isalnum(c) != 0 || c == '_' || c >= 0x80;
    };

    if (!is_valid_start(static_cast<unsigned char>(base[0]))) {
      throw std::invalid_argument("Invalid generator in power format: " + base);
    }
    for (size_t i = 1; i < base.size(); ++i) {
      if (!is_valid_cont(static_cast<unsigned char>(base[i]))) {
        throw std::invalid_argument("Invalid generator in power format: " +
                                    base);
      }
    }
  }

  std::string exp_str = pow_str.substr(pos + 1);
  if (exp_str.empty()) {
    throw std::invalid_argument("Invalid exponent format");
  }

  bool is_negative = false;
  if (exp_str[0] == '-') {
    is_negative = true;
    exp_str = exp_str.substr(1);
    if (exp_str.empty()) {
      throw std::invalid_argument("Invalid exponent format");
    }
  }

  // Strictly require every remaining character to be an ASCII decimal digit.
  // Rejects leading '+', duplicate '-', whitespace, trailing garbage (e.g.,
  // "g^3abc"), and secondary '^'.
  for (char ch : exp_str) {
    if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
      throw std::invalid_argument("Invalid exponent format: " + pow_str);
    }
  }

  uint64_t raw_exp = 0;
  try {
    size_t parsed_len = 0;
    raw_exp = std::stoull(exp_str, &parsed_len, 10);
    if (parsed_len != exp_str.size()) {
      throw std::invalid_argument("Invalid exponent format: " + pow_str);
    }
  } catch (const std::out_of_range&) {
    throw std::invalid_argument("Exponent out of range: " + pow_str);
  }

  uint64_t field_order = static_cast<uint64_t>(field.Order());
  if (field_order <= 1) {
    throw std::invalid_argument("Invalid field order for power parsing");
  }
  uint64_t group_order = field_order - 1;

  // Reduce 64-bit exponent modulo group_order BEFORE narrowing to uint32_t.
  uint64_t reduced_exp = raw_exp % group_order;
  if (is_negative && reduced_exp != 0) {
    reduced_exp = group_order - reduced_exp;
  }

  ElementType generator = field.MultiplicativeGenerator();
  return field.Pow(generator, static_cast<uint32_t>(reduced_exp));
}

}  // namespace utils
}  // namespace xg

#endif  // XGALOIS_UTILS_FIELD_HPP_
