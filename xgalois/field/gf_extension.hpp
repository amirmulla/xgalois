#ifndef XGALOIS_FIELD_EXTENSION_HPP_
#define XGALOIS_FIELD_EXTENSION_HPP_

// C system headers
#include <cassert>
#include <cmath>
#include <cstdint>

// C++ standard library headers
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

// Project headers
#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {

//------------------------------------------------------------------------------
// GaloisFieldExtension Class
//------------------------------------------------------------------------------

// Represents a Galois field extension GF(q^m) over a base field GF(q).
// Elements are polynomials modulo an irreducible polynomial.
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

  // Constructs a Galois field extension GF(q^m) with optional polynomial
  // string.
  GaloisFieldExtension(std::pair<uint64_t, uint64_t> prime_exp,
                       const std::string &irreducible_poly = "",
                       const std::string &rep = "poly",
                       const std::string &poly_var = "α",
                       bool check_irreducible = false,
                       bool prime_testing = false)
      : base_field_(static_cast<BaseFieldElementType>(prime_exp.first),
                    "int", prime_testing),
        irreducible_poly_(CreateIrreduciblePolynomial(
            base_field_, irreducible_poly, prime_exp.second)),
        representation_(utils::ConvertRepresentation(rep)), variable_name_(poly_var) {

    // Validate parameters first before using them
    assert(prime_exp.first >= 2 && "Prime must be at least 2.");
    assert(prime_exp.second >= 1 && "Exponent must be at least 1.");
    assert(irreducible_poly_.Degree() > 0 &&
           "Irreducible polynomial degree must be positive.");
    assert(utils::ConvertRepresentation(rep) == FieldRepresentation::POLY &&
           "Only polynomial representation is supported for extensions.");

    // Check irreducibility if requested
    if (check_irreducible && !utils::IsIrreducible(irreducible_poly_)) {
      throw std::invalid_argument("Provided polynomial is not irreducible");
    }

    // Calculate field order using safe integer power to avoid precision issues
    field_order_ = static_cast<uint32_t>(utils::SafeIntegerPower(base_field_.Order(),
                                           irreducible_poly_.Degree()));

    assert(field_order_ > 0 && "Field order must be positive.");
  }

  // Constructs a Galois field extension by decomposing the field order.
  // This constructor accepts a field order and automatically decomposes it
  // into prime and exponent components, similar to the factory pattern.
  GaloisFieldExtension(uint64_t order,
                       const std::string &irreducible_poly = "",
                       const std::string &rep = "poly",
                       const std::string &poly_var = "α",
                       bool check_irreducible = false,
                       bool prime_testing = false)
      : GaloisFieldExtension(DecomposeOrderToPrimeExp(order), irreducible_poly,
                             rep, poly_var, check_irreducible, prime_testing) {
    // Delegate to the main constructor after decomposing the order
  }

  // Returns the characteristic of the field, same as the base field.
  uint32_t Characteristic() const override {
    return base_field_.Characteristic();
  }

  // Returns the order (number of elements) of the field, |BaseField|^degree.
  uint32_t Order() const override { return field_order_; }

  // Returns the irreducible polynomial used as the modulus.
  PolynomialType Modulus() const { return irreducible_poly_; }

  // Adds two field elements (polynomial addition).
  PolynomialType Add(const PolynomialType &a,
                     const PolynomialType &b) const override {
    return (a + b);
  }

  // Subtracts one field element from another (polynomial subtraction).
  PolynomialType Sub(const PolynomialType &a,
                     const PolynomialType &b) const override {
    return (a - b);
  }

  // Multiplies two field elements (polynomial multiplication modulo modulus).
  PolynomialType Mul(const PolynomialType &a,
                     const PolynomialType &b) const override {
    return (a * b) % irreducible_poly_;
  }

  // Divides one field element by another (a * b_inv).
  PolynomialType Div(const PolynomialType &a,
                     const PolynomialType &b) const override {
    if (b == AdditiveIdentity()) {
      throw std::domain_error("Division by zero in GaloisFieldExtension");
    }
    return Mul(a, Inv(b));
  }

  // Returns the additive inverse (negation) of a field element.
  PolynomialType Neg(const PolynomialType &a) const override { return -a; }

  // Returns the multiplicative inverse of a non-zero field element.
  PolynomialType Inv(const PolynomialType &a) const override {
    if (a == AdditiveIdentity()) {
      throw std::domain_error("Multiplicative inverse of zero is undefined in "
                              "GaloisFieldExtension");
    }
    // Uses Extended Euclidean Algorithm: a * inv + k * modulus = gcd(a,
    // modulus). For irreducible modulus, gcd is a non-zero constant.
    auto [gcd, factors] =
        utils::PolynomialDenseExtendedGcd(a, irreducible_poly_);

    if (gcd.Degree() > 0 || gcd == AdditiveIdentity()) {
      throw std::runtime_error(
          "ExtendedGcd did not return a constant for inverse calculation. "
          "Ensure the modulus is irreducible and input is non-zero.");
    }
    // Assumes PolynomialDenseExtendedGcd provides factors.first as the inverse
    // or a form that becomes the inverse modulo irreducible_poly_ (e.g., if gcd
    // is 1 or factors.first is already scaled by gcd.Inv()).
    return factors.first % irreducible_poly_;
  }

  // Raises a field element to a power (polynomial exponentiation modulo
  // modulus).
  PolynomialType Pow(const PolynomialType &base, uint32_t exp) const override {
    // Use modular exponentiation to avoid computing large powers
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

  // Computes the square root of an element. Not implemented.
  PolynomialType Sqrt(const PolynomialType & /*a*/) const override {
    throw std::logic_error(
        "Sqrt not implemented for generic Galois field extensions");
  }

  // Computes the logarithm of an element. Not implemented.
  uint32_t Log(const PolynomialType & /*a*/) const override {
    throw std::logic_error(
        "Log not implemented for generic Galois field extensions");
  }

  // Computes the logarithm of an element with respect to a generator. Not
  // implemented.
  uint32_t Log(const PolynomialType & /*a*/,
               const PolynomialType & /*generator*/) const override {
    throw std::logic_error("Log with generator not implemented for generic "
                           "Galois field extensions");
  }

  // Returns a random element from the field.
  PolynomialType Random() const override {
    // Generates a random polynomial with degree less than
    // irreducible_poly_.Degree().
    auto base_field_ptr = std::make_shared<BaseField>(base_field_);
    std::vector<BaseFieldElement> coeffs(irreducible_poly_.Degree());
    for (size_t i = 0; i < irreducible_poly_.Degree(); ++i) {
      coeffs[i] = BaseFieldElement(base_field_.Random(), base_field_ptr);
    }
    // Ensure coeffs is not empty if degree is 0 (field_order_ = q^1,
    // irreducible_poly_.Degree() = 1)
    if (coeffs.empty() && irreducible_poly_.Degree() > 0) {
      // This case should ideally not happen if Degree() > 0 means at least one
      // coeff For degree 0 (constant), coeffs would have 1 element. If
      // irreducible_poly_.Degree() is 1, coeffs will have 1 element. If
      // irreducible_poly_.Degree() is 0, this loop doesn't run. However,
      // constructor asserts Degree > 0. If Degree is 1, coeffs has size 1. If
      // Degree is > 0 and coeffs becomes empty, it's an issue. For safety, if
      // coeffs ends up empty (e.g. degree 0 poly), create a zero poly. This
      // path should be reviewed based on how degree 0 polys are handled.
      // Assuming Degree() > 0 means coeffs will not be empty.
    }
    if (coeffs.empty()) { // Should only happen if irreducible_poly_.Degree()
                          // was 0, but asserted > 0
      BaseFieldElement kGfZero(base_field_.AdditiveIdentity(), base_field_ptr);
      return PolynomialType({kGfZero}, variable_name_);
    }
    return PolynomialType(coeffs, variable_name_); // Use variable_name_
  }

  // Returns the multiplicative identity (polynomial 1).
  PolynomialType MultiplicativeIdentity() const override {
    auto base_field_ptr = std::make_shared<BaseField>(base_field_);
    BaseFieldElement kGfOne(base_field_.MultiplicativeIdentity(),
                            base_field_ptr);
    return PolynomialType({kGfOne}, variable_name_); // Use variable_name_
  }

  // Returns the additive identity (polynomial 0).
  PolynomialType AdditiveIdentity() const override {
    auto base_field_ptr = std::make_shared<BaseField>(base_field_);
    BaseFieldElement kGfZero(base_field_.AdditiveIdentity(), base_field_ptr);
    return PolynomialType({kGfZero}, variable_name_); // Use variable_name_
  }

  PolynomialType GetElementValue(const PolynomialType &value) const override {
    return value;
  }

  PolynomialType SetElementValue(const PolynomialType &value) const override {
    // For extension fields, reduce the polynomial modulo the irreducible polynomial
    return value % irreducible_poly_;
  }  PolynomialType SetElementValue(const std::string &value_str) const override {
    // Handle power representation (e.g., "g^5" or "g^-3")
    if (value_str.find("g^") == 0) {
      try {
        PolynomialType value = utils::ParsePowerString(value_str, *this);
        return SetElementValue(value);
      } catch (const std::exception &) {
        throw std::invalid_argument("Invalid power format: " + value_str);
      }
    }

    // Handle polynomial representation (e.g., "α^2 + α + 1")
    if (value_str.find(variable_name_) != std::string::npos ||
        value_str.find("x") != std::string::npos ||
        value_str.find("+") != std::string::npos ||
        value_str.find("-") != std::string::npos ||
        value_str == "0" || value_str == "1") {
      try {
        auto base_field_ptr = std::make_shared<BaseField>(base_field_);
        PolynomialType poly = utils::ParsePolynomial(base_field_ptr, value_str, variable_name_);
        return SetElementValue(poly);
      } catch (const std::exception &) {
        throw std::invalid_argument("Invalid polynomial format: " + value_str);
      }
    }

    throw std::invalid_argument(
        "Unsupported string format for GF(q^m) element: " + value_str +
        ". Expected formats: g^k (power) or polynomial with " + variable_name_);
  }

  // Returns a multiplicative generator. Not implemented.
  PolynomialType MultiplicativeGenerator() const override {
    throw std::logic_error("MultiplicativeGenerator not implemented for "
                           "generic Galois field extensions");
  }

  // Returns all multiplicative generators. Not implemented.
  std::vector<PolynomialType> MultiplicativeGenerators() const override {
    throw std::logic_error("MultiplicativeGenerators not implemented for "
                           "generic Galois field extensions");
  }

  // Prints a description of the Galois field extension to the stream.
  void Print(std::ostream &os) const override {
    os << "GaloisFieldExtension GF(" << base_field_.Order() << "^"
       << irreducible_poly_.Degree() << ")";
    os << " over GF(" << base_field_.Order() << ")";
    os << " with modulus " << irreducible_poly_;
    os << " (elements are polynomials in " << variable_name_
       << ")"; // Added variable name info

    switch (representation_) {
    case FieldRepresentation::POLY:
      os << " [Rep: POLY]";
      break;
    default:
      os << " [Rep: Other]";
      break;
    }
  }

  // Prints a field element to the stream.
  void Print(const PolynomialType &a, std::ostream &os) const override {
    os << a;
  }

  // Converts a field element to its string representation.
  std::string ToString(const PolynomialType &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_;
  }

  void SetRepresentation(FieldRepresentation rep) override {
    // Currently, only POLY is asserted in constructor.
    // If other representations are added, ensure they are handled.
    representation_ = rep;
  }

private:
  // Helper function to decompose field order into prime and exponent
  static std::pair<uint64_t, uint64_t> DecomposeOrderToPrimeExp(uint64_t order) {
    // Validate field order constraints
    if (order < 2) {
      throw std::invalid_argument("Order of finite field must be at least 2");
    }

    // Validate field order is within supported range
    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument(
          "Field order exceeds maximum supported size (uint32_t)");
    }

    // Decompose order into prime power components using utility function
    auto [prime, exponent] = utils::DecomposePrimePower(order);
    if (prime == 0) {
      throw std::invalid_argument("Order must be a prime power");
    }

    return std::make_pair(prime, exponent);
  }

  // Helper function to create irreducible polynomial from string or database
  static PolynomialType
  CreateIrreduciblePolynomial(const BaseField &base_field,
                              const std::string &irreducible_poly_str,
                              uint64_t degree) {

    if (!irreducible_poly_str.empty()) {
      // Use provided polynomial string
      return utils::ParsePolynomial<BaseField>(
          std::make_shared<BaseField>(base_field), irreducible_poly_str);
    }

    // Get polynomial from database - try Conway polynomial first
    try {
      std::string conway_poly = databases::GetConwayPolynomial(
          static_cast<int>(base_field.Characteristic()),
          static_cast<int>(degree));
      return utils::ParsePolynomial<BaseField>(
          std::make_shared<BaseField>(base_field), conway_poly);
    } catch (const std::runtime_error &) {
      // Fall back to irreducible polynomial database
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

//------------------------------------------------------------------------------
// Type alias for convenience
template <typename BaseFieldElementType = uint32_t>
using GFPX = GaloisFieldExtension<BaseFieldElementType>;
//------------------------------------------------------------------------------

} // namespace xg

#endif // XGALOIS_FIELD_EXTENSION_HPP_
