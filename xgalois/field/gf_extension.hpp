#ifndef XGALOIS_FIELD_EXTENSION_HPP_
#define XGALOIS_FIELD_EXTENSION_HPP_

#include <cassert>
#include <cstdint>

#include <sstream>
#include <stdexcept>
#include <vector>

#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {

template <typename BaseFieldElementType = uint32_t>
class GaloisFieldExtension
    : public GaloisFieldBase<
          PolynomialDense<GaloisFieldPrime<BaseFieldElementType>>> {
  static_assert(std::is_integral<BaseFieldElementType>::value,
                "BaseFieldElementType must be an integral type.");
  static_assert(sizeof(BaseFieldElementType) <= sizeof(uint32_t),
                "BaseFieldElementType must be at most 32 bits to prevent "
                "overflow issues.");

 public:
  using BaseField = GaloisFieldPrime<BaseFieldElementType>;
  using BaseFieldElement = GaloisFieldElementBase<BaseField>;
  using PolynomialType = PolynomialDense<BaseField>;

  explicit GaloisFieldExtension(std::pair<uint64_t, uint64_t> prime_exp,
                                const std::string &irreducible_poly = "",
                                const std::string &rep = "poly",
                                const std::string &poly_var = "α",
                                bool check_irreducible = false,
                                bool prime_testing = false)
      : base_field_(static_cast<BaseFieldElementType>(prime_exp.first), "int",
                    prime_testing),
        irreducible_poly_(CreateIrreduciblePolynomial(
            base_field_, irreducible_poly, prime_exp.second)),
        representation_(utils::ConvertRepresentation(rep)),
        variable_name_(poly_var) {

    assert(prime_exp.first >= 2 && "Prime must be at least 2.");
    assert(prime_exp.second >= 1 && "Exponent must be at least 1.");
    assert(irreducible_poly_.Degree() > 0 &&
           "Irreducible polynomial degree must be positive.");
    assert(utils::ConvertRepresentation(rep) == FieldRepresentation::POLY &&
           "Only polynomial representation is supported for extensions.");

    if (check_irreducible && !utils::IsIrreducible(irreducible_poly_)) {
      throw std::invalid_argument("Provided polynomial is not irreducible");
    }

    field_order_ = static_cast<uint32_t>(utils::SafeIntegerPower(
        base_field_.Order(), irreducible_poly_.Degree()));

    assert(field_order_ > 0 && "Field order must be positive.");
  }

  explicit GaloisFieldExtension(uint64_t order,
                                const std::string &irreducible_poly = "",
                                const std::string &rep = "poly",
                                const std::string &poly_var = "α",
                                bool check_irreducible = false,
                                bool prime_testing = false)
      : GaloisFieldExtension(DecomposeOrderToPrimeExp(order), irreducible_poly,
                             rep, poly_var, check_irreducible, prime_testing) {

  }

  uint32_t Characteristic() const override {
    return base_field_.Characteristic();
  }

  uint32_t Order() const override { return field_order_; }

  PolynomialType Modulus() const { return irreducible_poly_; }

  PolynomialType Add(const PolynomialType &a,
                     const PolynomialType &b) const override {
    return (a + b);
  }

  PolynomialType Sub(const PolynomialType &a,
                     const PolynomialType &b) const override {
    return (a - b);
  }

  PolynomialType Mul(const PolynomialType &a,
                     const PolynomialType &b) const override {
    return (a * b) % irreducible_poly_;
  }

  PolynomialType Div(const PolynomialType &a,
                     const PolynomialType &b) const override {
    if (b == AdditiveIdentity()) {
      throw std::domain_error("Division by zero in GaloisFieldExtension");
    }
    return Mul(a, Inv(b));
  }

  PolynomialType Neg(const PolynomialType &a) const override { return -a; }

  PolynomialType Inv(const PolynomialType &a) const override {
    if (a == AdditiveIdentity()) {
      throw std::domain_error(
          "Multiplicative inverse of zero is undefined in "
          "GaloisFieldExtension");
    }

    auto [gcd, factors] =
        utils::PolynomialDenseExtendedGcd(a, irreducible_poly_);

    if (gcd.Degree() > 0 || gcd == AdditiveIdentity()) {
      throw std::runtime_error(
          "ExtendedGcd did not return a constant for inverse calculation. "
          "Ensure the modulus is irreducible and input is non-zero.");
    }

    return factors.first % irreducible_poly_;
  }

  PolynomialType Pow(const PolynomialType &base, uint32_t exp) const override {

    if (exp == 0) {
      return MultiplicativeIdentity();
    }
    if (exp == 1) {
      return base % irreducible_poly_;
    }

    PolynomialType result = MultiplicativeIdentity();
    PolynomialType current_base = base % irreducible_poly_;
    uint64_t current_exp = exp;

    while (current_exp > 0) {
      if (current_exp % 2 == 1) {
        result = (result * current_base) % irreducible_poly_;
      }
      current_base = (current_base * current_base) % irreducible_poly_;
      current_exp /= 2;
    }

    return result;
  }

  PolynomialType Sqrt(const PolynomialType & ) const override {
    throw std::logic_error(
        "Sqrt not implemented for generic Galois field extensions");
  }

  uint32_t Log(const PolynomialType & ) const override {
    throw std::logic_error(
        "Log not implemented for generic Galois field extensions");
  }

  uint32_t Log(const PolynomialType & ,
               const PolynomialType & ) const override {
    throw std::logic_error(
        "Log with generator not implemented for generic "
        "Galois field extensions");
  }

  PolynomialType Random() const override {

    auto base_field_ptr = std::make_shared<BaseField>(base_field_);
    std::vector<BaseFieldElement> coeffs(irreducible_poly_.Degree());
    for (size_t i = 0; i < irreducible_poly_.Degree(); ++i) {
      coeffs[i] = BaseFieldElement(base_field_.Random(), base_field_ptr);
    }

    if (coeffs.empty() && irreducible_poly_.Degree() > 0) {

    }
    if (coeffs.empty()) {

      BaseFieldElement kGfZero(base_field_.AdditiveIdentity(), base_field_ptr);
      return PolynomialType({kGfZero}, variable_name_);
    }
    return PolynomialType(coeffs, variable_name_);
  }

  PolynomialType MultiplicativeIdentity() const override {
    auto base_field_ptr = std::make_shared<BaseField>(base_field_);
    BaseFieldElement kGfOne(base_field_.MultiplicativeIdentity(),
                            base_field_ptr);
    return PolynomialType({kGfOne}, variable_name_);
  }

  PolynomialType AdditiveIdentity() const override {
    auto base_field_ptr = std::make_shared<BaseField>(base_field_);
    BaseFieldElement kGfZero(base_field_.AdditiveIdentity(), base_field_ptr);
    return PolynomialType({kGfZero}, variable_name_);
  }

  PolynomialType GetElementValue(const PolynomialType &value) const override {
    return value;
  }

  PolynomialType SetElementValue(const PolynomialType &value) const override {

    return value % irreducible_poly_;
  }
  PolynomialType SetElementValue(const std::string &value_str) const override {

    if (value_str.find("g^") == 0) {
      try {
        PolynomialType value = utils::ParsePowerString(value_str, *this);
        return SetElementValue(value);
      } catch (const std::exception &) {
        throw std::invalid_argument("Invalid power format: " + value_str);
      }
    }

    if (value_str.find(variable_name_) != std::string::npos ||
        value_str.find('x') != std::string::npos ||
        value_str.find('+') != std::string::npos ||
        value_str.find('-') != std::string::npos || value_str == "0" ||
        value_str == "1") {
      try {
        auto base_field_ptr = std::make_shared<BaseField>(base_field_);
        PolynomialType poly =
            utils::ParsePolynomial(base_field_ptr, value_str, variable_name_);
        return SetElementValue(poly);
      } catch (const std::exception &) {
        throw std::invalid_argument("Invalid polynomial format: " + value_str);
      }
    }

    throw std::invalid_argument(
        "Unsupported string format for GF(q^m) element: " + value_str +
        ". Expected formats: g^k (power) or polynomial with " + variable_name_);
  }

  PolynomialType MultiplicativeGenerator() const override {
    throw std::logic_error(
        "MultiplicativeGenerator not implemented for "
        "generic Galois field extensions");
  }

  std::vector<PolynomialType> MultiplicativeGenerators() const override {
    throw std::logic_error(
        "MultiplicativeGenerators not implemented for "
        "generic Galois field extensions");
  }

  void Print(std::ostream &os) const override {
    os << "GaloisFieldExtension GF(" << base_field_.Order() << "^"
       << irreducible_poly_.Degree() << ")";
    os << " over GF(" << base_field_.Order() << ")";
    os << " with modulus " << irreducible_poly_;
    os << " (elements are polynomials in " << variable_name_
       << ")";

    switch (representation_) {
      case FieldRepresentation::POLY:
        os << " [Rep: POLY]";
        break;
      default:
        os << " [Rep: Other]";
        break;
    }
  }

  void Print(const PolynomialType &a, std::ostream &os) const override {
    os << a;
  }

  std::string ToString(const PolynomialType &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_;
  }

  void SetRepresentation(FieldRepresentation rep) override {

    representation_ = rep;
  }

 private:

  static std::pair<uint64_t, uint64_t> DecomposeOrderToPrimeExp(
      uint64_t order) {

    if (order < 2) {
      throw std::invalid_argument("Order of finite field must be at least 2");
    }

    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument(
          "Field order exceeds maximum supported size (uint32_t)");
    }

    auto [prime, exponent] = utils::DecomposePrimePower(order);
    if (prime == 0) {
      throw std::invalid_argument("Order must be a prime power");
    }

    return std::make_pair(prime, exponent);
  }

  static PolynomialType CreateIrreduciblePolynomial(
      const BaseField &base_field, const std::string &irreducible_poly_str,
      uint64_t degree) {
    if (!irreducible_poly_str.empty()) {

      return utils::ParsePolynomial<BaseField>(
          std::make_shared<BaseField>(base_field), irreducible_poly_str);
    }

    try {
      std::string conway_poly = databases::GetConwayPolynomial(
          static_cast<int>(base_field.Characteristic()),
          static_cast<int>(degree));
      return utils::ParsePolynomial<BaseField>(
          std::make_shared<BaseField>(base_field), conway_poly);
    } catch (const std::runtime_error &) {

      try {
        std::string irreducible_poly = databases::GetIrreduciblePolynomial(
            static_cast<int>(base_field.Characteristic()),
            static_cast<int>(degree));
        return utils::ParsePolynomial<BaseField>(
            std::make_shared<BaseField>(base_field), irreducible_poly);
      } catch (const std::runtime_error &) {
        throw std::invalid_argument(
            "No suitable polynomial found in databases for GF(" +
            std::to_string(base_field.Characteristic()) + "^" +
            std::to_string(degree) +
            "). Please provide modulus polynomial manually.");
      }
    }
  }

  BaseField base_field_;
  const PolynomialType irreducible_poly_;
  uint32_t field_order_;
  FieldRepresentation representation_;
  std::string variable_name_;
};

template <typename BaseFieldElementType = uint32_t>
using GFPX = GaloisFieldExtension<BaseFieldElementType>;

}

#endif
