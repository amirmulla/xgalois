#ifndef XGALOIS_FIELD_GF_PRIME_HPP_
#define XGALOIS_FIELD_GF_PRIME_HPP_

// C system headers
#include <cassert>
#include <cstdint>

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

// Project headers
#include "xgalois/field/gf_base.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

//------------------------------------------------------------------------------
// GaloisFieldPrime Class
//------------------------------------------------------------------------------

// Represents a Galois Field of prime order p, GF(p).
// Operations are performed modulo p.
template <typename ElementType = uint32_t>
class GaloisFieldPrime : public GaloisFieldBase<ElementType> {
  static_assert(std::is_integral<ElementType>::value,
                "ElementType must be an integral type.");
  static_assert(
      sizeof(ElementType) <= sizeof(uint32_t),
      "ElementType must be at most 32 bits to prevent overflow issues.");

 public:
  // Constructs a prime field with characteristic p.
  // Initializes the prime characteristic, representation, and finds a
  // multiplicative generator.
  explicit GaloisFieldPrime(ElementType prime, const std::string &rep = "int",
                            bool prime_testing = false)
      : p_(prime), representation_(utils::ConvertRepresentation(rep)) {
    if (p_ < 2) {
      throw std::invalid_argument("Prime must be >= 2.");
    }

    // Check if the provided number is actually prime when prime_testing is
    // enabled
    if (prime_testing && !xg::utils::IsPrime(static_cast<long long>(p_))) {
      throw std::invalid_argument("Provided value is not prime.");
    }

    generator_ = 0;  // Will be set by MultiplicativeGenerator().

    assert(p_ > 0 && "Prime must be positive.");
    assert((representation_ == FieldRepresentation::INT ||
            representation_ == FieldRepresentation::HEX ||
            representation_ == FieldRepresentation::POW ||
            representation_ == FieldRepresentation::LOG) &&
           "Invalid field representation.");

    if (this->p_ > 2) {
      this->generator_ = this->MultiplicativeGenerator();
    } else if (this->p_ == 2) {
      // For GF(2), 1 is the only non-zero element and thus the generator.
      this->generator_ = static_cast<ElementType>(1);
    }
  }

  uint32_t Characteristic() const override { return static_cast<uint32_t>(p_); }
  uint32_t Order() const override { return static_cast<uint32_t>(p_); }
  uint64_t Modulus() const { return static_cast<uint64_t>(p_); }

  ElementType MultiplicativeIdentity() const override {
    return static_cast<ElementType>(1);
  }

  ElementType AdditiveIdentity() const override {
    return static_cast<ElementType>(0);
  }

  ElementType GetElementValue(const ElementType &value) const override {
    return value;
  }

  ElementType SetElementValue(const ElementType &value) const override {
    // For prime fields, values should be in range [0, p-1]
    return static_cast<ElementType>(static_cast<uint64_t>(value) %
                                    static_cast<uint64_t>(p_));
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    // Handle power representation (e.g., "g^5" or "g^-3")
    if (value_str.find("g^") == 0) {
      try {
        ElementType value = utils::ParsePowerString(value_str, *this);
        return SetElementValue(value);
      } catch (const std::exception &) {
        throw std::invalid_argument("Invalid power format: " + value_str);
      }
    }

    throw std::invalid_argument(
        "Unsupported string format for GF(p) element: " + value_str +
        ". Expected format: g^k, where g is the generator.");
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_;
  }
  void SetRepresentation(FieldRepresentation rep) override {
    representation_ = rep;
  }

  // (a + b) mod p.
  inline ElementType Add(const ElementType &a,
                         const ElementType &b) const override {
    return static_cast<ElementType>(
        (static_cast<uint64_t>(a) + static_cast<uint64_t>(b)) %
        static_cast<uint64_t>(this->p_));
  }

  // (a - b) mod p.
  // Ensures a non-negative result before modulo.
  inline ElementType Sub(const ElementType &a,
                         const ElementType &b) const override {
    return static_cast<ElementType>((static_cast<uint64_t>(this->p_) +
                                     static_cast<uint64_t>(a) -
                                     static_cast<uint64_t>(b)) %
                                    static_cast<uint64_t>(this->p_));
  }

  // (a * b) mod p.
  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    return static_cast<ElementType>(
        (static_cast<uint64_t>(a) * static_cast<uint64_t>(b)) %
        static_cast<uint64_t>(this->p_));
  }

  // (a / b) mod p, equivalent to a * b^(-1) mod p.
  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero.");
    }
    return Mul(a, Inv(b));
  }

  // (-a) mod p.
  inline ElementType Neg(const ElementType &a) const override {
    return static_cast<ElementType>(
        (static_cast<uint64_t>(this->p_) - static_cast<uint64_t>(a)) %
        static_cast<uint64_t>(this->p_));
  }

  // Multiplicative inverse a^(-1) mod p using Fermat's Little Theorem.
  // a^(p-2) mod p = a^(-1) mod p for a != 0.
  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Element has no inverse (is zero).");
    }
    // For p=2, p-2=0, Pow(a,0)=1. inv(1)=1 in GF(2). Correct.
    return Pow(a, static_cast<uint64_t>(this->p_) - 2);
  }

  // base^exp mod p using square-and-multiply algorithm.
  inline ElementType Pow(const ElementType &base, uint32_t exp) const override {
    uint64_t result = 1;
    uint64_t current_base = static_cast<uint64_t>(base);
    uint64_t current_modulus = static_cast<uint64_t>(this->p_);
    uint32_t current_exp = exp;

    if (current_base == 0) {
      // By convention, 0^0 = 1.
      return (exp == 0) ? static_cast<ElementType>(1)
                        : static_cast<ElementType>(0);
    }

    while (current_exp > 0) {
      if (current_exp & 1) {  // If exp is odd.
        result = (result * current_base) % current_modulus;
      }
      current_base = (current_base * current_base) % current_modulus;
      current_exp >>= 1;  // exp = exp / 2.
    }
    return static_cast<ElementType>(result);
  }

  // Square root. Not yet implemented.
  // TODO(amirmulla): Implement square root for GF(p).
  inline ElementType Sqrt(const ElementType &a) const override {
    throw std::logic_error("sqrt not implemented for GF(p).");
  }

  // Computes discrete logarithm log_generator(a).
  // Finds k such that generator^k = a by iterative search.
  uint32_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Logarithm of zero is undefined.");
    }
    if (generator == 0) {
      throw std::invalid_argument("Generator for Log cannot be zero.");
    }
    if (this->p_ > 2 && generator == 1) {
      throw std::invalid_argument(
          "Generator for Log cannot be 1 for fields larger than GF(2).");
    }
    if (this->p_ == 2 && generator != 1) {
      throw std::invalid_argument("Generator for Log must be 1 for GF(2).");
    }

    if (a == this->MultiplicativeIdentity()) {  // log_g(1) is always 0.
      return 0;
    }

    uint64_t phi = this->Order() - 1;  // Order of the multiplicative group.
    if (phi == 0) {                    // Only for p=2 (phi=1) or invalid p=1.
      if (this->p_ == 2 && a == 1 && generator == 1) return 0;
      throw std::domain_error(
          "Logarithm undefined for this field configuration.");
    }

    ElementType current_power = generator;
    // Check powers from g^1 up to g^(phi-1).
    for (uint64_t k = 1; k < phi; ++k) {
      if (current_power == a) {
        return k;
      }
      current_power = Mul(current_power, generator);
      // Optimization: if g is not primitive, and we cycle back to 1 early.
      if (current_power == this->MultiplicativeIdentity()) {
        break;
      }
    }
    throw std::domain_error(
        "Element not found in the cyclic group generated by the "
        "generator, or generator is not primitive for this element.");
  }

  // Computes discrete logarithm using the field's default generator.
  uint32_t Log(const ElementType &a) const override {
    return Log(a, this->generator_);
  }

  // Generates a random element in [0, p-1].
  // Uses a thread-local RNG for thread safety.
  ElementType Random() const override {
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_int_distribution<uint64_t> dist(
        0, static_cast<uint64_t>(this->p_) - 1);
    return static_cast<ElementType>(dist(rng));
  }

  // Finds the first multiplicative generator (primitive root) of the field.
  // A generator g is an element whose powers generate all non-zero elements.
  // Algorithm: iterate g from 2 to p-1. For each g, check if
  // g^((p-1)/f) != 1 for all distinct prime factors f of p-1.
  ElementType MultiplicativeGenerator() const override {
    // If called after construction for p_ > 2, generator_ is already set.
    if (this->generator_ != 0 && this->p_ > 2) {
      return this->generator_;
    }
    if (this->p_ == 2) {  // Generator for GF(2) is 1.
      return static_cast<ElementType>(1);
    }
    if (this->p_ < 2) {  // Should be caught by constructor.
      throw std::runtime_error(
          "Multiplicative generator not applicable for "
          "fields smaller than GF(2).");
    }

    uint64_t phi = static_cast<uint64_t>(this->p_) - 1;
    std::vector<long long> prime_factors = utils::PrimeFactorize(phi);
    prime_factors.erase(std::unique(prime_factors.begin(), prime_factors.end()),
                        prime_factors.end());

    for (ElementType candidate = 2; candidate < this->p_; ++candidate) {
      bool is_generator = true;
      for (long long factor : prime_factors) {
        uint64_t exp = phi / static_cast<uint64_t>(factor);
        if (Pow(candidate, exp) == 1) {
          is_generator = false;
          break;
        }
      }
      if (is_generator) {
        return candidate;
      }
    }
    throw std::runtime_error("No multiplicative generator found.");
  }

  // Finds all multiplicative generators of the field.
  std::vector<ElementType> MultiplicativeGenerators() const override {
    if (this->p_ < 2) {
      throw std::runtime_error(
          "Multiplicative generators not applicable for "
          "fields smaller than GF(2).");
    }
    if (this->p_ == 2) {
      return {static_cast<ElementType>(1)};
    }

    uint64_t phi = static_cast<uint64_t>(this->p_) - 1;
    std::vector<long long> prime_factors = utils::PrimeFactorize(phi);
    prime_factors.erase(std::unique(prime_factors.begin(), prime_factors.end()),
                        prime_factors.end());

    std::vector<ElementType> generators;
    for (ElementType candidate = 2; candidate < this->p_; ++candidate) {
      bool is_generator = true;
      for (long long factor : prime_factors) {
        uint64_t exp = phi / static_cast<uint64_t>(factor);
        if (Pow(candidate, exp) == 1) {
          is_generator = false;
          break;
        }
      }
      if (is_generator) {
        generators.push_back(candidate);
      }
    }

    if (generators.empty()) {
      throw std::runtime_error("No multiplicative generators found.");
    }
    return generators;
  }

  // Prints a description of the field.
  void Print(std::ostream &os) const override {
    os << "Galois Prime Field GF(" << this->p_ << ")";
    switch (this->representation_) {
      case FieldRepresentation::INT:
        os << " [Rep: INT]";
        break;
      case FieldRepresentation::HEX:
        os << " [Rep: HEX]";
        break;
      case FieldRepresentation::POW:
        os << " [Rep: POW]";
        break;
      case FieldRepresentation::LOG:
        os << " [Rep: LOG]";
        break;
      default:
        os << " [Rep: UNKNOWN]";
        break;
    }
  }

  // Prints a field element according to the current representation.
  void Print(const ElementType &a, std::ostream &os) const override {
    switch (this->representation_) {
      case FieldRepresentation::HEX:
        os << "0x" << std::hex << static_cast<uint64_t>(a) << std::dec;
        break;
      case FieldRepresentation::POW:
        if (a == 0) {
          os << "0";
        } else if (a == 1) {
          os << "g^0";  // Or 1, depending on convention for g^0.
        } else {
          try {
            uint64_t power = Log(a, this->generator_);
            os << "g^" << power;
          } catch (const std::domain_error &) {
            // Fallback if log fails (e.g. generator not primitive for 'a').
            os << static_cast<uint64_t>(a);
          }
        }
        break;
      case FieldRepresentation::LOG:
        if (a == 0) {
          os << "undefined";  // Log of 0 is undefined.
        } else {
          try {
            uint32_t log_val = Log(a, this->generator_);
            os << log_val;
          } catch (const std::domain_error &) {
            os << static_cast<uint64_t>(a);  // Fallback.
          }
        }
        break;
      case FieldRepresentation::INT:
      default:
        os << static_cast<uint64_t>(a);
        break;
    }
  }

  // Converts a field element to its string representation.
  std::string ToString(const ElementType &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

 protected:
  ElementType p_;                       // Prime characteristic of the field.
  FieldRepresentation representation_;  // Representation for printing elements.
  ElementType generator_;  // Default multiplicative generator (primitive root).

};  // class GaloisFieldPrime

//------------------------------------------------------------------------------
// GaloisFieldPrimeTable Class
//------------------------------------------------------------------------------

// Represents a Galois Field of prime order p, GF(p), using lookup tables
// for accelerated Mul, Div, Inv, and Pow operations.
template <typename ElementType = uint32_t>
class GaloisFieldPrimeTable : public GaloisFieldPrime<ElementType> {
  static_assert(
      sizeof(ElementType) <= sizeof(uint32_t),
      "ElementType must be at most 32 bits to prevent overflow issues.");

 public:
  // Constructs a table-based prime field.
  // Initializes tables for log and anti-log.
  explicit GaloisFieldPrimeTable(ElementType prime,
                                 const std::string &rep = "int",
                                 bool prime_testing = false)
      : GaloisFieldPrime<ElementType>(prime, rep, prime_testing) {
    InitializeTables();
  }

 private:
  // Precomputes logarithm (log_alpha_) and anti-logarithm (alpha_pow_) tables.
  // generator_^i = alpha_pow_[i]
  // log_alpha_[alpha_pow_[i]] = i
  void InitializeTables() {
    // For p=2, tables are trivial and might not offer performance benefits.
    // The base class operations are simple enough for GF(2).
    if (this->p_ == 2) {
      alpha_pow_.resize(this->p_);  // Size 2
      log_alpha_.resize(this->p_);  // Size 2
      // this->generator_ is already 1 from base constructor.
      alpha_pow_[0] = static_cast<ElementType>(1);  // g^0 = 1
      // log_alpha_[0] is undefined, often set to a special value like p-1.
      // log_alpha_[1] = 0 since g^0 = 1.
      if (this->p_ > 0)
        log_alpha_[0] = this->p_ - 1;  // Special value for log(0)
      if (this->p_ > 1) log_alpha_[1] = static_cast<ElementType>(0);
      return;
    }

    // Ensure generator is found using the appropriate (potentially overridden)
    // method before table initialization.
    // If this class overrides MultiplicativeGenerator, it will be used.
    // Otherwise, GaloisFieldPrime::MultiplicativeGenerator is used.
    this->generator_ = this->MultiplicativeGenerator();

    alpha_pow_.resize(this->p_);
    log_alpha_.resize(this->p_);

    alpha_pow_[0] = static_cast<ElementType>(1);  // generator^0 = 1.
    for (uint64_t i = 1; i < static_cast<uint64_t>(this->p_); ++i) {
      // Use base class Mul for table construction to avoid dependency on
      // potentially uninitialized tables of this class.
      alpha_pow_[i] = GaloisFieldPrime<ElementType>::Mul(alpha_pow_[i - 1],
                                                         this->generator_);
    }

    // log_alpha_[0] is undefined; use p-1 as a placeholder.
    std::fill(log_alpha_.begin(), log_alpha_.end(), this->p_ - 1);
    for (uint64_t i = 0; i < static_cast<uint64_t>(this->p_) - 1;
         ++i) {                        // Powers from g^0 to g^(p-2).
      if (alpha_pow_[i] < this->p_) {  // Sanity check.
        log_alpha_[alpha_pow_[i]] = static_cast<ElementType>(i);
      } else {
        throw std::runtime_error(
            "Error during log_alpha_ table construction: value out of bounds.");
      }
    }
  }

 public:
  // (a * b) mod p using lookup tables.
  // a * b = g^(log(a) + log(b)).
  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    if (a == 0 || b == 0) {
      return 0;
    }
    if (this->p_ == 2)
      return GaloisFieldPrime<ElementType>::Mul(a, b);  // Base for GF(2).
    // Check bounds for table access.
    if (a >= this->p_ || b >= this->p_) {
      throw std::out_of_range("Operands out of field range in table Mul.");
    }
    uint64_t log_a = log_alpha_[a];
    uint64_t log_b = log_alpha_[b];
    // Exponent is modulo (p-1) because that's the order of the mult. group.
    uint64_t log_result =
        (log_a + log_b) % (static_cast<uint64_t>(this->p_) - 1);
    return alpha_pow_[log_result];
  }

  // (a / b) mod p using lookup tables.
  // a / b = g^(log(a) - log(b)).
  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero.");
    }
    if (a == 0) {
      return 0;
    }
    if (this->p_ == 2)
      return GaloisFieldPrime<ElementType>::Div(a, b);  // Base for GF(2).
    if (a >= this->p_ || b >= this->p_) {
      throw std::out_of_range("Operands out of field range in table Div.");
    }
    uint64_t log_a = log_alpha_[a];
    uint64_t log_b = log_alpha_[b];
    uint64_t log_result =
        (static_cast<uint64_t>(this->p_) - 1 + log_a - log_b) %
        (static_cast<uint64_t>(this->p_) - 1);
    return alpha_pow_[log_result];
  }

  // a^(-1) mod p using lookup tables.
  // a^(-1) = g^(-log(a)) = g^((p-1) - log(a)).
  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Element has no inverse (is zero).");
    }
    if (this->p_ == 2)
      return GaloisFieldPrime<ElementType>::Inv(a);  // Base for GF(2).
    if (a >= this->p_) {
      throw std::out_of_range("Operand out of field range in table Inv.");
    }
    uint64_t log_a = log_alpha_[a];
    // (p-1) - log_a handles the negative exponent modulo (p-1).
    uint64_t log_result = (static_cast<uint64_t>(this->p_) - 1 - log_a) %
                          (static_cast<uint64_t>(this->p_) - 1);
    return alpha_pow_[log_result];
  }

  // base^exp mod p using lookup tables.
  // base^exp = g^(exp * log(base)).
  inline ElementType Pow(const ElementType &base, uint32_t exp) const override {
    if (this->p_ == 2) {
      return GaloisFieldPrime<ElementType>::Pow(base, exp);
    }
    if (base == 0) {
      return (exp == 0) ? static_cast<ElementType>(1)
                        : static_cast<ElementType>(0);
    }
    if (base >= this->p_) {
      throw std::out_of_range("Base out of field range in table Pow.");
    }
    uint64_t log_base = log_alpha_[base];
    uint64_t log_result = (static_cast<uint64_t>(log_base) * exp) %
                          (static_cast<uint64_t>(this->p_) - 1);
    return alpha_pow_[log_result];
  }

  // Computes log_generator(a). Uses table if `generator` is the field's
  // default generator; otherwise, delegates to base class.
  uint32_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Logarithm of zero is undefined.");
    }
    if (this->p_ == 2) {
      return GaloisFieldPrime<ElementType>::Log(a, generator);
    }

    if (generator == this->generator_) {
      if (a >= this->p_) {
        throw std::out_of_range(
            "Element out of field range in table Log lookup.");
      }
      // log_alpha_[0] is a special value; a!=0 is guaranteed here.
      return static_cast<uint32_t>(this->log_alpha_[a]);
    } else {
      // Fallback to iterative search for non-default generator.
      return GaloisFieldPrime<ElementType>::Log(a, generator);
    }
  }

  // Computes log_g(a) using the field's table-specific (default) generator.
  uint32_t Log(const ElementType &a) const override {
    if (this->p_ == 2) {
      return GaloisFieldPrime<ElementType>::Log(a, this->generator_);
    }
    if (a == 0) {
      throw std::domain_error("Logarithm of zero is undefined.");
    }
    if (a >= this->p_) {
      throw std::out_of_range(
          "Element out of field range in table Log lookup.");
    }
    return static_cast<uint32_t>(this->log_alpha_[a]);
  }

  // Finds the first multiplicative generator.
  // Must use base class Pow to avoid using tables during their own setup.
  ElementType MultiplicativeGenerator() const override {
    if (this->p_ <= 2) {
      // GF(2) generator is 1, handled by base.
      // This override is primarily for ensuring base Pow is used.
      // Throwing here if p_ <= 2 to match base class behavior if called
      // directly. However, InitializeTables calls this, and for p_=2, it has
      // special path. If p_ > 2, this method proceeds. If p_ == 2,
      // InitializeTables handles it, this won't be called for table init. If
      // called externally for p_ == 2, base class would return 1. This override
      // is for p_ > 2 context during table init.
      if (this->p_ == 2)
        return static_cast<ElementType>(1);  // Should be consistent.
      throw std::runtime_error(
          "Multiplicative generator context error for p <= 2 in table class.");
    }

    uint64_t phi = static_cast<uint64_t>(this->p_) - 1;
    std::vector<long long> prime_factors = utils::PrimeFactorize(phi);
    prime_factors.erase(std::unique(prime_factors.begin(), prime_factors.end()),
                        prime_factors.end());

    for (ElementType candidate = 2; candidate < this->p_; ++candidate) {
      bool is_generator = true;
      for (long long factor : prime_factors) {
        uint64_t exp = phi / static_cast<uint64_t>(factor);
        // Critical: Use GaloisFieldPrime::Pow to avoid recursion if tables
        // are not yet initialized or if this method is called during init.
        if (GaloisFieldPrime<ElementType>::Pow(candidate, exp) == 1) {
          is_generator = false;
          break;
        }
      }
      if (is_generator) {
        return candidate;
      }
    }
    throw std::runtime_error("No multiplicative generator found.");
  }

  // Finds all multiplicative generators.
  // Must use base class Pow.
  std::vector<ElementType> MultiplicativeGenerators() const override {
    if (this->p_ < 2) {
      throw std::runtime_error(
          "Multiplicative generators not applicable for "
          "fields smaller than GF(2).");
    }
    if (this->p_ == 2) {
      return {static_cast<ElementType>(1)};
    }

    uint64_t phi = static_cast<uint64_t>(this->p_) - 1;
    std::vector<long long> prime_factors = utils::PrimeFactorize(phi);
    prime_factors.erase(std::unique(prime_factors.begin(), prime_factors.end()),
                        prime_factors.end());

    std::vector<ElementType> generators;
    for (ElementType candidate = 2; candidate < this->p_; ++candidate) {
      bool is_generator = true;
      for (long long factor : prime_factors) {
        uint64_t exp = phi / static_cast<uint64_t>(factor);
        if (GaloisFieldPrime<ElementType>::Pow(candidate, exp) == 1) {
          is_generator = false;
          break;
        }
      }
      if (is_generator) {
        generators.push_back(candidate);
      }
    }

    if (generators.empty()) {
      throw std::runtime_error("No multiplicative generators found.");
    }
    return generators;
  }

  // Prints a description of the table-based field.
  void Print(std::ostream &os) const override {
    os << "Galois Prime Field GF(" << this->p_ << ") [Table-Based]";
    switch (this->representation_) {
      case FieldRepresentation::INT:
        os << " [Rep: INT]";
        break;
      case FieldRepresentation::HEX:
        os << " [Rep: HEX]";
        break;
      case FieldRepresentation::POW:
        os << " [Rep: POW]";
        break;
      case FieldRepresentation::LOG:
        os << " [Rep: LOG]";
        break;
      default:
        os << " [Rep: UNKNOWN]";
        break;
    }
  }

 private:
  std::vector<ElementType> alpha_pow_;  // alpha_pow_[i] = g^i.
  std::vector<ElementType> log_alpha_;  // log_alpha_[val] = i where g^i = val.
};  // class GaloisFieldPrimeTable

//------------------------------------------------------------------------------
// Type aliases for convenience
template <typename ElementType = uint32_t>
using GFP = GaloisFieldPrime<ElementType>;
template <typename ElementType = uint32_t>
using GFPLOG = GaloisFieldPrimeTable<ElementType>;
//------------------------------------------------------------------------------

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_PRIME_HPP_
