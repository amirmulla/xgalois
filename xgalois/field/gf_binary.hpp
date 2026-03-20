#ifndef XGALOIS_FIELD_GF_BINARY_HPP_
#define XGALOIS_FIELD_GF_BINARY_HPP_

// C system headers
#include <cassert>
#include <cstdint>

// C++ standard library headers
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project headers
#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {

//------------------------------------------------------------------------------
// GaloisFieldBinary Class - Represents GF(2)
//------------------------------------------------------------------------------
class GaloisFieldBinary : public GaloisFieldBase<uint8_t> {
 public:
  GaloisFieldBinary(const std::string &rep = "int")
      : representation_(utils::ConvertRepresentation(rep)) {
    // Ensure the representation is valid for GF(2)
    if (representation_ != FieldRepresentation::INT &&
        representation_ != FieldRepresentation::HEX) {
      throw std::invalid_argument(
          "GF(2) supports only INT and HEX representations for its elements.");
    }
  }

  uint32_t Characteristic() const override { return 2; }
  uint32_t Order() const override { return 2; }
  uint64_t Modulus() const { return 2; }

  inline uint8_t Add(const uint8_t &a, const uint8_t &b) const override {
    return a ^ b;
  }
  inline uint8_t Sub(const uint8_t &a, const uint8_t &b) const override {
    return a ^ b;
  }
  inline uint8_t Mul(const uint8_t &a, const uint8_t &b) const override {
    return a & b;
  }
  inline uint8_t Div(const uint8_t &a, const uint8_t &b) const override {
    if (b == 0) throw std::domain_error("Division by zero in GF(2)");
    return a;  // a/1 = a
  }
  inline uint8_t Neg(const uint8_t &a) const override {
    return a;
  }  // -a = a in GF(2)
  inline uint8_t Inv(const uint8_t &a) const override {
    if (a == 0) throw std::domain_error("Inverse of zero in GF(2)");
    return 1;  // 1^-1 = 1
  }
  inline uint8_t Pow(const uint8_t &a, uint32_t exp) const override {
    if (a == 0) return (exp == 0) ? 1 : 0;  // 0^0 = 1, 0^i = 0 for i > 0
    return 1;                               // 1^exp = 1
  }
  inline uint8_t Sqrt(const uint8_t &a) const override {
    return a;
  }  // sqrt(a) = a in GF(2)

  uint32_t Log(const uint8_t &a) const override {
    if (a == 0) throw std::domain_error("Log of zero is undefined in GF(2)");
    if (a == 1) return 0;  // GF(2)^* = {1}, generator is 1, log_1(1) = 0.
    // Should not be reached if 'a' is a valid non-zero element (i.e., 1)
    throw std::invalid_argument("Log defined only for 1 in GF(2)");
  }

  uint32_t Log(const uint8_t &a, const uint8_t &generator) const override {
    if (generator != 1)
      throw std::invalid_argument("Generator in GF(2) must be 1");
    return Log(a);  // Delegate to single-argument Log
  }

  uint8_t Random() const override {
    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<uint8_t> distrib(0, 1);
    return distrib(gen);
  }

  uint8_t MultiplicativeIdentity() const override { return 1; }
  uint8_t AdditiveIdentity() const override { return 0; }

  uint8_t GetElementValue(const uint8_t &value) const override { return value; }
  uint8_t SetElementValue(const uint8_t &value) const override { return value; }
  uint8_t SetElementValue(const std::string &value_str) const override {
    throw std::invalid_argument(
        "GF(2) does not support string representation for elements.");
  }

  uint8_t MultiplicativeGenerator() const override { return 1; }
  std::vector<uint8_t> MultiplicativeGenerators() const override { return {1}; }

  void Print(std::ostream &os) const override { os << "GF(2)"; }
  void Print(const uint8_t &a, std::ostream &os) const override {
    if (representation_ == FieldRepresentation::INT) {
      os << static_cast<uint64_t>(a);
    } else if (representation_ == FieldRepresentation::HEX) {
      os << "0x" << std::hex << static_cast<uint64_t>(a) << std::dec;
    } else {
      throw std::invalid_argument("Unsupported representation for GF(2)");
    }
  }

  std::string ToString(const uint8_t &a) const override {
    if (representation_ == FieldRepresentation::INT) {
      return std::to_string(static_cast<uint64_t>(a));
    } else if (representation_ == FieldRepresentation::HEX) {
      std::ostringstream oss;
      oss << "0x" << std::hex << static_cast<uint64_t>(a) << std::dec;
      return oss.str();
    }
    throw std::invalid_argument("Unsupported representation for GF(2)");
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_;
  }
  void SetRepresentation(FieldRepresentation rep) override {
    // For GF(2), INT, HEX representations of 0 and 1 are trivial.
    // POLY, POW and LOG are not very meaningful.
    if (rep != FieldRepresentation::INT && rep != FieldRepresentation::HEX) {
      throw std::invalid_argument(
          "GF(2) supports INT, HEX representation for its elements.");
    }
    representation_ = rep;
  }

 private:
  FieldRepresentation representation_;
};  // End of GaloisFieldBinary class

//------------------------------------------------------------------------------
// GaloisFieldBinaryExtension Class - Represents GF(2^m)
//------------------------------------------------------------------------------
template <typename ElementType = uint32_t>
class GaloisFieldBinaryExtension : public GaloisFieldBase<ElementType> {
 public:
  // Ensure ElementType is an unsigned integral type
  static_assert(std::is_unsigned<ElementType>::value &&
                    std::is_integral<ElementType>::value,
                "ElementType must be an unsigned integral type");

  // Constructor with string parameters
  GaloisFieldBinaryExtension(uint8_t m, const std::string &rep = "int",
                             const std::string &irreducible_poly = "",
                             const std::string &variable_name = "α",
                             bool check_irreducible = false,
                             const std::string &generator_name = "g")
      : m_(m),
        representation_(utils::ConvertRepresentation(rep)),
        variable_name_(variable_name),
        generator_name_(generator_name) {
    // Validate parameters
    // Degree must be between 1 and (sizeof(ElementType) * 8) - 1, to prevent
    // overflow and ensure we can represent the field elements.
    if (m_ < 1 || m_ > (sizeof(ElementType) * 8) - 1) {
      throw std::invalid_argument(
          "Degree must be between 1 and " +
          std::to_string((sizeof(ElementType) * 8) - 1));
    }

    // If no irreducible polynomial is provided, retrieve from database
    std::string poly_str;
    if (irreducible_poly.empty()) {
      // Try to get Conway polynomial from database first
      try {
        poly_str = databases::GetConwayPolynomial(2, static_cast<uint64_t>(m_));
      } catch (const std::runtime_error &) {
        // Fall back to irreducible polynomial database
        try {
          poly_str =
              databases::GetIrreduciblePolynomial(2, static_cast<uint64_t>(m_));
        } catch (const std::runtime_error &) {
          throw std::invalid_argument(
              "No suitable polynomial found in databases for GF(2^" +
              std::to_string(m_) +
              "). Please provide irreducible_poly manually.");
        }
      }
    } else {
      poly_str = irreducible_poly;
    }

    // Parse the irreducible polynomial using ParsePolynomial
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto parsed_poly = utils::ParsePolynomial(gf2, poly_str);

    // Check polynomial degree matches the field extension
    if (parsed_poly.Degree() != m_) {
      throw std::invalid_argument("Irreducible polynomial degree (" +
                                  std::to_string(parsed_poly.Degree()) +
                                  ") must equal field extension degree (" +
                                  std::to_string(m_) + ")");
    }

    // Check irreducibility if requested
    if (check_irreducible && !utils::IsIrreducible(parsed_poly)) {
      throw std::invalid_argument("Provided polynomial is not irreducible");
    }

    // Convert PolyDense to uint64_t representation using utility function
    irreducible_poly_ =
        utils::BinaryPolynomialToUint<GaloisFieldBinary>(parsed_poly);

    // Basic validation of irreducible polynomial
    if (irreducible_poly_ <= 1) {
      throw std::invalid_argument("Invalid irreducible polynomial");
    }

    // Calculate field order
    order_ = static_cast<uint32_t>(1) << m_;
    group_order_ = order_ - 1;
  }

  uint32_t Characteristic() const override { return 2; }
  uint32_t Order() const override { return order_; }
  uint64_t Modulus() const { return irreducible_poly_; }
  uint8_t Degree() const { return m_; }

  // Addition is just XOR (same as in GF(2))
  inline ElementType Add(const ElementType &a,
                         const ElementType &b) const override {
    return a ^ b;
  }

  // Subtraction is also XOR in characteristic 2
  inline ElementType Sub(const ElementType &a,
                         const ElementType &b) const override {
    return a ^ b;
  }

  // Multiplication using polynomial multiplication mod irreducible polynomial
  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    if (a == 0 || b == 0) return 0;

    // Step 1: Polynomial multiplication (carry-free)
    uint64_t product = 0;
    for (uint8_t i = 0; i < m_; i++) {
      if (b & (ElementType(1) << i)) {
        product ^= (static_cast<uint64_t>(a) << i);
      }
    }

    // Step 2: Reduction modulo irreducible polynomial
    for (int i = 2 * m_ - 2; i >= m_; i--) {
      if (product & (1ULL << i)) {
        product ^= (static_cast<uint64_t>(irreducible_poly_) << (i - m_));
      }
    }

    return static_cast<ElementType>(product);
  }

  // Negation in characteristic 2 is identity
  inline ElementType Neg(const ElementType &a) const override { return a; }

  // Inverse using Extended Euclidean Algorithm
  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }
    return ExtendedEuclideanAlgorithm(irreducible_poly_, a);
  }

  // Division using multiplication by inverse
  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }
    return Mul(a, Inv(b));
  }

  // Exponentiation by squaring
  inline ElementType Pow(const ElementType &a, uint32_t exp) const override {
    if (a == 0) {
      return (exp == 0) ? 1 : 0;
    }

    ElementType result = 1;
    ElementType base = a;

    while (exp > 0) {
      if (exp & 1) {
        result = Mul(result, base);
      }
      base = Mul(base, base);
      exp >>= 1;
    }

    return result;
  }

  // Square root - for binary fields, can be implemented via repeated squaring
  // because (a^2)^(2^(m-1)) = a in GF(2^m)
  inline ElementType Sqrt(const ElementType &a) const override {
    if (a == 0) return 0;

    // Calculate 2^(m-1) - 1 which is the exponent needed
    uint64_t exp = (1ULL << (m_ - 1)) - 1;
    return Pow(a, exp);
  }

  // Logarithm base g (multiplicative generator)
  uint32_t Log(const ElementType &a) const override {
    return Log(a, MultiplicativeGenerator());
  }

  // Logarithm with specified generator
  uint32_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }

    if (generator == 0) {
      throw std::invalid_argument("Generator cannot be zero");
    }

    // Naive implementation: try all powers until we find a match
    ElementType power = 1;
    uint32_t result = 0;
    uint32_t group_order = group_order_;

    while (result < group_order) {
      if (power == a) {
        return result;
      }
      power = Mul(power, generator);
      result++;
    }

    throw std::runtime_error("Element not in the group or invalid generator");
  }

  // Random element generation
  ElementType Random() const override {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<uint64_t> distrib(
        0, static_cast<uint64_t>(group_order_));
    return static_cast<ElementType>(distrib(gen));
  }

  ElementType MultiplicativeIdentity() const override { return 1; }

  ElementType AdditiveIdentity() const override { return 0; }

  ElementType GetElementValue(const ElementType &value) const override {
    return value;
  }

  ElementType SetElementValue(const ElementType &value) const override {
    return value;
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    // Try power format using the configured generator name (e.g., "g^5")
    // Check for exact pattern: generator_name followed immediately by ^
    std::string power_pattern = generator_name_ + "^";
    size_t pos = value_str.find(power_pattern);
    if (pos != std::string::npos) {
      // Verify it's at word boundary (start of string or preceded by
      // non-alphanumeric)
      if (pos == 0 || !std::isalnum(value_str[pos - 1])) {
        // Verify there's a number after the ^
        size_t after_caret = pos + power_pattern.length();
        if (after_caret < value_str.length() &&
            (std::isdigit(value_str[after_caret]) ||
             value_str[after_caret] == '-')) {
          return utils::ParsePowerString(value_str, *this);
        }
      }
    }

    // Try polynomial format using the configured variable name (e.g., "α^2 + α
    // + 1")
    if (value_str.find(variable_name_) != std::string::npos ||
        value_str.find("x") != std::string::npos ||
        value_str.find("+") != std::string::npos ||
        value_str.find("-") != std::string::npos) {
      return ParsePolynomialString(value_str);
    }

    throw std::invalid_argument(
        "Unsupported string format for GF(2^m) element: " + value_str +
        ". Expected formats: " + generator_name_ + "^n, polynomial with " +
        variable_name_ + ".");
  }

  // Find a multiplicative generator
  ElementType MultiplicativeGenerator() const override {
    // A naive approach: try elements until one has order 2^m - 1
    for (ElementType candidate = 2; candidate < Order(); candidate++) {
      if (HasFullOrder(candidate)) {
        return candidate;
      }
    }
    throw std::runtime_error("Could not find a primitive element");
  }

  // Get all multiplicative generators (expensive for large fields)
  std::vector<ElementType> MultiplicativeGenerators() const override {
    std::vector<ElementType> generators;
    ElementType max_element = Order();

    for (ElementType candidate = 2; candidate < max_element; candidate++) {
      if (HasFullOrder(candidate)) {
        generators.push_back(candidate);
      }
    }

    return generators;
  }

  void Print(std::ostream &os) const override {
    os << "GF(2^" << static_cast<uint64_t>(m_) << ")";
  }

  void Print(const ElementType &a, std::ostream &os) const override {
    if (representation_ == FieldRepresentation::INT) {
      os << a;
    } else if (representation_ == FieldRepresentation::HEX) {
      os << "0x" << std::hex << a << std::dec;
    } else if (representation_ == FieldRepresentation::POLY) {
      os << BinaryPolyToString(a);
    } else if (representation_ == FieldRepresentation::POW) {
      os << "g^" << Log(a);
    } else if (representation_ == FieldRepresentation::LOG) {
      os << Log(a);
    } else {
      throw std::invalid_argument("Unsupported representation for GF(2^m)");
    }
  }

  std::string ToString(const ElementType &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_;
  }

  void SetRepresentation(FieldRepresentation rep) override {
    if (rep != FieldRepresentation::INT && rep != FieldRepresentation::HEX &&
        rep != FieldRepresentation::POLY && rep != FieldRepresentation::POW &&
        rep != FieldRepresentation::LOG) {
      throw std::invalid_argument(
          "GF(2^m) supports INT, HEX, POLY, POW, or LOG "
          "representation for its elements.");
    }
    representation_ = rep;
  }

 protected:
  uint8_t m_;                  // Extension degree
  uint64_t irreducible_poly_;  // Irreducible polynomial
  uint32_t order_;             // Field order = 2^m
  uint32_t group_order_;       // Multiplicative group order = 2^m - 1
  FieldRepresentation representation_;
  std::string variable_name_;   // Variable name for polynomial representation
  std::string generator_name_;  // Name for multiplicative generator

  // Computes polynomial inverse using Extended Euclidean Algorithm
  ElementType ExtendedEuclideanAlgorithm(const ElementType &poly,
                                         const ElementType &element) const {
    uint64_t r0 = static_cast<uint64_t>(poly);
    uint64_t r1 = static_cast<uint64_t>(element);
    uint64_t s0 = 1;
    uint64_t s1 = 0;
    uint64_t t0 = 0;
    uint64_t t1 = 1;

    while (r1 != 0) {
      // Calculate quotient and remainder
      uint64_t q = 0;
      int deg_diff = HighestBit64(r0) - HighestBit64(r1);

      if (deg_diff >= 0) {
        for (int i = deg_diff; i >= 0; i--) {
          if ((r0 & (1ULL << (HighestBit64(r1) + i)))) {
            r0 ^= (r1 << i);
            q ^= (1ULL << i);
          }
        }
      }

      // Update coefficients
      uint64_t temp_r = r0;
      r0 = r1;
      r1 = temp_r;

      uint64_t temp_s = s0 ^ Mul64(q, s1);
      s0 = s1;
      s1 = temp_s;

      uint64_t temp_t = t0 ^ Mul64(q, t1);
      t0 = t1;
      t1 = temp_t;
    }

    return static_cast<ElementType>(t0);
  }

  // Helper method to find the highest bit position
  int HighestBit(ElementType n) const {
    if (n == 0) return -1;

    int position = 0;
    while (n != 0) {
      n >>= 1;
      position++;
    }
    return position - 1;
  }

  // Helper method to find the highest bit position for uint64_t
  int HighestBit64(uint64_t n) const {
    if (n == 0) return -1;

    int position = 0;
    while (n != 0) {
      n >>= 1;
      position++;
    }
    return position - 1;
  }

  // 64-bit multiplication for Extended Euclidean Algorithm
  uint64_t Mul64(uint64_t a, uint64_t b) const {
    if (a == 0 || b == 0) return 0;

    uint64_t product = 0;
    for (uint8_t i = 0; i < m_; i++) {
      if (b & (1ULL << i)) {
        product ^= (a << i);
      }
    }

    for (int i = 2 * m_ - 2; i >= m_; i--) {
      if (product & (1ULL << i)) {
        product ^= (static_cast<uint64_t>(irreducible_poly_) << (i - m_));
      }
    }

    return product;
  }

  // Check if an element has the full multiplicative order
  bool HasFullOrder(ElementType element) const {
    ElementType power = element;

    // Check if element^k = 1 for any k less than the full order
    for (uint64_t k = 2; k < Order(); k++) {
      if (power == 1) {
        return false;
      }
      power = GaloisFieldBinaryExtension<ElementType>::Mul(power, element);
    }

    // element^(2^m - 1) should be 1
    return power == 1;
  }

  // Convert polynomial representation to string
  std::string BinaryPolyToString(ElementType poly) const {
    if (poly == 0) return "0";

    std::ostringstream oss;
    bool first = true;

    for (int i = m_ - 1; i >= 0; i--) {
      if (poly & (static_cast<ElementType>(1) << i)) {
        if (!first) oss << " + ";
        first = false;

        if (i == 0) {
          oss << "1";
        } else if (i == 1) {
          oss << variable_name_;
        } else {
          oss << variable_name_ << "^" << i;
        }
      }
    }

    return oss.str();
  }

  // Parse polynomial string like "α^2 + α + 1" using utils::ParsePolynomial
  ElementType ParsePolynomialString(const std::string &poly_str) const {
    if (poly_str == "0") return 0;

    // Create a GF(2) field for parsing
    auto gf2 = std::make_shared<GaloisFieldBinary>();

    // Parse polynomial using utils::ParsePolynomial with variable name
    auto parsed_poly = utils::ParsePolynomial(gf2, poly_str, variable_name_);

    // Check polynomial degree doesn't exceed field extension degree
    if (parsed_poly.Degree() >= m_) {
      throw std::invalid_argument(
          "Polynomial degree (" + std::to_string(parsed_poly.Degree()) +
          ") must be less than field extension degree (" + std::to_string(m_) +
          ")");
    }

    // Convert polynomial coefficients to binary representation using utility
    // function
    return utils::BinaryPolynomialToUint<GaloisFieldBinary>(parsed_poly);
  }
};

//------------------------------------------------------------------------------
// GFBELogTables Class - Implements GF(2^m) using logarithm tables
//------------------------------------------------------------------------------
template <typename ElementType = uint32_t>
class GFBELogTables : public GaloisFieldBinaryExtension<ElementType> {
 public:
  // Constructor
  GFBELogTables(uint8_t m, const std::string &rep = "int",
                const std::string &irreducible_poly = "",
                const std::string &variable_name = "α",
                bool check_irreducible = false,
                const std::string &generator_name = "g")
      : GaloisFieldBinaryExtension<ElementType>(
            m, rep, irreducible_poly, variable_name, check_irreducible,
            generator_name),
        is_generator_cached_(false) {  // Initialize the flag
    PrecomputeTables();
  }

  // Override arithmetic operations to use log tables
  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    if (a == 0 || b == 0)  // Handle zero multiplication
      return 0;

    ElementType log_sum = (log_table_[a] + log_table_[b]);
    if (log_sum >= this->group_order_) {
      log_sum -= this->group_order_;  // Handle wrap-around
    }
    return exp_table_[log_sum];
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0)  // Handle division by zero
      throw std::domain_error("Division by zero");
    if (a == 0)  // Handle zero division
      return 0;

    ElementType log_diff =
        (log_table_[a] > log_table_[b])
            ? (log_table_[a] - log_table_[b])
            : (log_table_[a] + this->group_order_ - log_table_[b]);

    return exp_table_[log_diff];
  }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0)  // Handle inverse of zero
      throw std::domain_error("Inverse of zero is undefined");

    return exp_table_[(this->group_order_ - log_table_[a])];
  }

  inline ElementType Pow(const ElementType &a, uint32_t exp) const override {
    if (a == 0) return (exp == 0) ? 1 : 0;  // 0^0 = 1, 0^k = 0 for k > 0

    return exp_table_[(static_cast<uint64_t>(log_table_[a]) * exp) %
                      this->group_order_];
  }

  // Logarithm using precomputed table
  uint32_t Log(const ElementType &a) const override {
    if (a == 0)  // Handle log of zero
      throw std::domain_error("Log of zero is undefined");

    return log_table_[a];
  }

  uint32_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0)  // Handle log of zero
      throw std::domain_error("Log of zero is undefined");

    if (generator == this->MultiplicativeGenerator())  // Use cached generator
      return log_table_[a];

    return GaloisFieldBinaryExtension<ElementType>::Log(a, generator);
  }

  // Override to return the cached generator
  ElementType MultiplicativeGenerator() const override {
    if (is_generator_cached_) {
      return cached_generator_;
    }
    // This path should ideally not be taken if PrecomputeTables always runs
    // first. However, to be safe, compute and cache it if not already done.
    cached_generator_ =
        GaloisFieldBinaryExtension<ElementType>::MultiplicativeGenerator();
    is_generator_cached_ = true;
    return cached_generator_;
  }

 protected:
  std::vector<uint32_t> log_table_;     // Maps field element to its log
  std::vector<ElementType> exp_table_;  // Maps log value to field element
  mutable ElementType
      cached_generator_;              // Stores the generator used for tables
  mutable bool is_generator_cached_;  // Flag to check if generator is cached

  static constexpr ElementType LOG_ZERO =
      std::numeric_limits<ElementType>::max();

  virtual void PrecomputeTables() {
    // Initialize tables
    log_table_.resize(this->order_, 0);
    exp_table_.resize(this->group_order_, 0);

    // Find a multiplicative generator
    // Use the base class method to find it the first time, which will be cached
    // by our overridden MultiplicativeGenerator()
    ElementType generator = this->MultiplicativeGenerator();
    // At this point, generator is stored in cached_generator_ and
    // is_generator_cached_ is true.

    // Compute exp table: g^i
    ElementType element = 1;  // g^0 = 1
    for (uint32_t i = 0; i < this->group_order_; i++) {
      exp_table_[i] = element;  // g^i
      log_table_[element] = i;  // log_g(g^i) = i
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }

    // Log of 0 is undefined, but we'll set it to a sentinel value
    log_table_[0] = LOG_ZERO;
  }
};

//------------------------------------------------------------------------------
// GFBELogTablesOpt Class - Optimized GF(2^m) using larger logarithm tables
//------------------------------------------------------------------------------
template <typename ElementType = uint32_t>
class GFBELogTablesOpt : public GFBELogTables<ElementType> {
 public:
  // Constructor
  GFBELogTablesOpt(uint8_t m, const std::string &rep = "int",
                   const std::string &irreducible_poly = "",
                   const std::string &variable_name = "α",
                   bool check_irreducible = false,
                   const std::string &generator_name = "g")
      : GFBELogTables<ElementType>(m, rep, irreducible_poly, variable_name,
                                   check_irreducible, generator_name) {}

  // Optimized arithmetic operations without zero checks
  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    return this->exp_table_[(this->log_table_[a]) + (this->log_table_[b])];
  }

  // Division with zero checks, with modulo handling
  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0)  // Handle division by zero
      throw std::domain_error("Division by zero");

    return this->exp_table_[this->group_order_ + this->log_table_[a] -
                            this->log_table_[b]];
  }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0)  // Handle inverse of zero
      throw std::domain_error("Inverse of zero is undefined");

    return this->exp_table_[this->group_order_ - this->log_table_[a]];
  }

 protected:
  void PrecomputeTables() override {
    // Initialize optimized tables
    this->log_table_.resize(this->order_, 0);
    this->exp_table_.resize(4 * this->group_order_ + 1, 0);

    // Find a multiplicative generator
    ElementType generator = this->MultiplicativeGenerator();

    // Compute exp table with extended wraparound
    ElementType element = 1;  // g^0 = 1
    for (uint32_t i = 0; i < this->group_order_; i++) {
      this->exp_table_[i] = element;                       // g^i
      this->exp_table_[i + this->group_order_] = element;  // wraparound copy
      this->log_table_[element] = i;                       // log_g(g^i) = i
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }

    // Set special values for zero handling
    this->log_table_[0] =
        2 * this->group_order_;  // Set log(0) to 2*group_order
  }
};

//------------------------------------------------------------------------------
// GFBEZechLogTables Class - Implements GF(2^m) using Zech logarithm tables
//------------------------------------------------------------------------------
template <typename ElementType = uint32_t>
class GFBEZechLogTables : public GaloisFieldBinaryExtension<ElementType> {
 public:
  // Constructor
  GFBEZechLogTables(uint8_t m, const std::string &rep = "log",
                    const std::string &irreducible_poly = "",
                    const std::string &variable_name = "α",
                    bool check_irreducible = false,
                    const std::string &generator_name = "g")
      : GaloisFieldBinaryExtension<ElementType>(
            m, rep, irreducible_poly, variable_name, check_irreducible,
            generator_name) {
    PrecomputeTables();
  }

  // Override arithmetic operations assuming inputs are in log form
  inline ElementType Add(const ElementType &log_a,
                         const uint32_t &log_b) const override {
    if (log_a == LOG_ZERO) return log_b;  // log(0) + log(b) = log(b)
    if (log_b == LOG_ZERO) return log_a;  // log(a) + log(0) = log(a)

    uint32_t diff = (log_b >= log_a) ? (log_b - log_a)
                                     : (this->group_order_ - log_a + log_b);

    return (zech_table_[diff] == LOG_ZERO)  // Check for infinity
               ? LOG_ZERO
               : ((log_a + zech_table_[diff]) <
                  this->group_order_)  // Check for wrap
                     ? (log_a + zech_table_[diff])
                     : (log_a + zech_table_[diff] - this->group_order_);
  }

  inline ElementType Sub(const ElementType &log_a,
                         const ElementType &log_b) const override {
    return Add(log_a, log_b);
  }

  inline ElementType Mul(const ElementType &log_a,
                         const ElementType &log_b) const override {
    if (log_a == LOG_ZERO || log_b == LOG_ZERO) return LOG_ZERO;  // log(0)

    return ((log_a + log_b) < this->group_order_)
               ? (log_a + log_b)
               : (log_a + log_b - this->group_order_);
  }

  inline ElementType Div(const ElementType &log_a,
                         const ElementType &log_b) const override {
    if (log_b == LOG_ZERO) throw std::domain_error("Division by zero");

    return (log_a < log_b)       ? (log_a + this->group_order_ - log_b)
           : (log_a == LOG_ZERO) ? LOG_ZERO
                                 : (log_a - log_b);
  }

  inline ElementType Inv(const ElementType &log_a) const override {
    if (log_a == LOG_ZERO)
      throw std::domain_error("Inverse of zero is undefined");

    return ((log_a == 0) ? 0 : (this->group_order_ - log_a));
  }

  inline ElementType Pow(const ElementType &log_a,
                         uint32_t exp) const override {
    if (log_a == LOG_ZERO) return (exp == 0) ? 0 : LOG_ZERO;  // log(1) : log(0)

    return (static_cast<uint64_t>(log_a) * exp) % this->group_order_;
  }

  // Logarithm is identity since input is already in log form
  uint32_t Log(const ElementType &log_a) const override {
    if (log_a == LOG_ZERO) {
      throw std::domain_error("Log of zero is undefined");
    }
    return static_cast<uint32_t>(log_a);
  }

  uint32_t Log(const ElementType &log_a,
               const ElementType &generator) const override {
    if (log_a == LOG_ZERO) {
      throw std::domain_error("Log of zero is undefined");
    }
    if (generator == this->MultiplicativeGenerator()) {
      return static_cast<uint32_t>(log_a);
    }

    ElementType a = exp_table_[log_a];
    return GaloisFieldBinaryExtension<ElementType>::Log(generator);
  }

  ElementType GetElementValue(const ElementType &log_a) const override {
    return (log_a == LOG_ZERO) ? 0 : exp_table_[log_a];
  }

  ElementType SetElementValue(const ElementType &value) const override {
    return log_table[value];
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    ElementType value =
        GaloisFieldBinaryExtension<ElementType>::SetElementValue(value_str);
    return log_table[value];
  }

 protected:
  std::vector<ElementType> zech_table_;
  std::vector<ElementType> exp_table_;
  std::vector<ElementType> log_table;

  static constexpr ElementType LOG_ZERO =
      std::numeric_limits<ElementType>::max();

  void PrecomputeTables() {
    zech_table_.resize(this->group_order_, LOG_ZERO);
    exp_table_.resize(this->group_order_, 0);
    log_table.resize(this->order_, LOG_ZERO);

    // Create log table first for O(1) lookups
    ElementType generator = this->MultiplicativeGenerator();
    ElementType element;

    // Build log table
    element = 1;  // Start with g^0 = 1
    for (uint32_t i = 0; i < this->group_order_; i++) {
      exp_table_[i] = element;  // g^i
      log_table[element] = i;   // log_g(g^i) = i
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }

    // Build Zech table
    element = 1;  // Start with g^0 = 1
    for (uint32_t i = 0; i < this->group_order_; i++) {
      // Compute 1 + g^i
      ElementType one_plus_element =
          GaloisFieldBinaryExtension<ElementType>::Add(1, element);
      if (one_plus_element == 0)  // if (1 + g^i) == 0 -> log(0) = infinity
        zech_table_[i] = LOG_ZERO;
      else  // Otherwise, store the log value
        zech_table_[i] = log_table[one_plus_element];
      // Compute next element, g^(i+1) = g^i * g
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }
  }
};

//------------------------------------------------------------------------------
// Type alias
using GF2 = GaloisFieldBinary;
template <typename ElementType = uint32_t>
using GF2X = GaloisFieldBinaryExtension<ElementType>;
template <typename ElementType = uint32_t>
using GF2XLOG = GFBELogTablesOpt<ElementType>;
using GF2XZECH = GFBEZechLogTables<uint32_t>;
//------------------------------------------------------------------------------

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_BINARY_HPP_
