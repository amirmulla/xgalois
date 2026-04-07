#ifndef XGALOIS_FIELD_GF_BINARY_HPP_
#define XGALOIS_FIELD_GF_BINARY_HPP_

#include <cassert>
#include <cstdint>

#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {

class GaloisFieldBinary : public GaloisFieldBase<uint8_t> {
 public:
  explicit GaloisFieldBinary(const std::string &rep = "int")
      : representation_(utils::ConvertRepresentation(rep)) {

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
    return a;
  }
  inline uint8_t Neg(const uint8_t &a) const override {
    return a;
  }
  inline uint8_t Inv(const uint8_t &a) const override {
    if (a == 0) throw std::domain_error("Inverse of zero in GF(2)");
    return 1;
  }
  inline uint8_t Pow(const uint8_t &a, uint32_t exp) const override {
    if (a == 0) return (exp == 0) ? 1 : 0;
    return 1;
  }
  inline uint8_t Sqrt(const uint8_t &a) const override {
    return a;
  }

  uint32_t Log(const uint8_t &a) const override {
    if (a == 0) throw std::domain_error("Log of zero is undefined in GF(2)");
    if (a == 1) return 0;

    throw std::invalid_argument("Log defined only for 1 in GF(2)");
  }

  uint32_t Log(const uint8_t &a, const uint8_t &generator) const override {
    if (generator != 1) {
      throw std::invalid_argument("Generator in GF(2) must be 1");
    }
    return Log(a);
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

    if (rep != FieldRepresentation::INT && rep != FieldRepresentation::HEX) {
      throw std::invalid_argument(
          "GF(2) supports INT, HEX representation for its elements.");
    }
    representation_ = rep;
  }

 private:
  FieldRepresentation representation_;
};

template <typename ElementType = uint32_t>
class GaloisFieldBinaryExtension : public GaloisFieldBase<ElementType> {
 public:

  static_assert(std::is_unsigned<ElementType>::value &&
                    std::is_integral<ElementType>::value,
                "ElementType must be an unsigned integral type");

  explicit GaloisFieldBinaryExtension(uint8_t m, const std::string &rep = "int",
                                      const std::string &irreducible_poly = "",
                                      const std::string &variable_name = "α",
                                      bool check_irreducible = false,
                                      const std::string &generator_name = "g")
      : m_(m),
        representation_(utils::ConvertRepresentation(rep)),
        variable_name_(variable_name),
        generator_name_(generator_name) {

    if (m_ < 1 || m_ > (sizeof(ElementType) * 8) - 1) {
      throw std::invalid_argument(
          "Degree must be between 1 and " +
          std::to_string((sizeof(ElementType) * 8) - 1));
    }

    std::string poly_str;
    if (irreducible_poly.empty()) {

      try {
        poly_str = databases::GetConwayPolynomial(2, static_cast<uint64_t>(m_));
      } catch (const std::runtime_error &) {

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

    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto parsed_poly = utils::ParsePolynomial(gf2, poly_str);

    if (parsed_poly.Degree() != m_) {
      throw std::invalid_argument("Irreducible polynomial degree (" +
                                  std::to_string(parsed_poly.Degree()) +
                                  ") must equal field extension degree (" +
                                  std::to_string(m_) + ")");
    }

    if (check_irreducible && !utils::IsIrreducible(parsed_poly)) {
      throw std::invalid_argument("Provided polynomial is not irreducible");
    }

    irreducible_poly_ =
        utils::BinaryPolynomialToUint<GaloisFieldBinary>(parsed_poly);

    if (irreducible_poly_ <= 1) {
      throw std::invalid_argument("Invalid irreducible polynomial");
    }

    order_ = static_cast<uint32_t>(1) << m_;
    group_order_ = order_ - 1;
  }

  uint32_t Characteristic() const override { return 2; }
  uint32_t Order() const override { return order_; }
  uint64_t Modulus() const { return irreducible_poly_; }
  uint8_t Degree() const { return m_; }

  inline ElementType Add(const ElementType &a,
                         const ElementType &b) const override {
    return a ^ b;
  }

  inline ElementType Sub(const ElementType &a,
                         const ElementType &b) const override {
    return a ^ b;
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    if (a == 0 || b == 0) return 0;

    uint64_t product = 0;
    for (uint8_t i = 0; i < m_; i++) {
      if (b & (ElementType(1) << i)) {
        product ^= (static_cast<uint64_t>(a) << i);
      }
    }

    for (int i = 2 * m_ - 2; i >= m_; i--) {
      if (product & (1ULL << i)) {
        product ^= (static_cast<uint64_t>(irreducible_poly_) << (i - m_));
      }
    }

    return static_cast<ElementType>(product);
  }

  inline ElementType Neg(const ElementType &a) const override { return a; }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }
    return ExtendedEuclideanAlgorithm(irreducible_poly_, a);
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }
    return Mul(a, Inv(b));
  }

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

  inline ElementType Sqrt(const ElementType &a) const override {
    if (a == 0) return 0;

    uint64_t exp = (1ULL << (m_ - 1)) - 1;
    return Pow(a, exp);
  }

  uint32_t Log(const ElementType &a) const override {
    return Log(a, MultiplicativeGenerator());
  }

  uint32_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }

    if (generator == 0) {
      throw std::invalid_argument("Generator cannot be zero");
    }

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

    std::string power_pattern = generator_name_ + "^";
    size_t pos = value_str.find(power_pattern);
    if (pos != std::string::npos) {

      if (pos == 0 || !std::isalnum(value_str[pos - 1])) {

        size_t after_caret = pos + power_pattern.length();
        if (after_caret < value_str.length() &&
            (std::isdigit(value_str[after_caret]) ||
             value_str[after_caret] == '-')) {
          return utils::ParsePowerString(value_str, *this);
        }
      }
    }

    if (value_str.find(variable_name_) != std::string::npos ||
        value_str.find('x') != std::string::npos ||
        value_str.find('+') != std::string::npos ||
        value_str.find('-') != std::string::npos) {
      return ParsePolynomialString(value_str);
    }

    throw std::invalid_argument(
        "Unsupported string format for GF(2^m) element: " + value_str +
        ". Expected formats: " + generator_name_ + "^n, polynomial with " +
        variable_name_ + ".");
  }

  ElementType MultiplicativeGenerator() const override {

    for (ElementType candidate = 2; candidate < Order(); candidate++) {
      if (HasFullOrder(candidate)) {
        return candidate;
      }
    }
    throw std::runtime_error("Could not find a primitive element");
  }

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
  uint8_t m_;
  uint64_t irreducible_poly_;
  uint32_t order_;
  uint32_t group_order_;
  FieldRepresentation representation_;
  std::string variable_name_;
  std::string generator_name_;

  ElementType ExtendedEuclideanAlgorithm(const ElementType &poly,
                                         const ElementType &element) const {
    uint64_t r0 = static_cast<uint64_t>(poly);
    uint64_t r1 = static_cast<uint64_t>(element);
    uint64_t s0 = 1;
    uint64_t s1 = 0;
    uint64_t t0 = 0;
    uint64_t t1 = 1;

    while (r1 != 0) {

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

  int HighestBit(ElementType n) const {
    if (n == 0) return -1;

    int position = 0;
    while (n != 0) {
      n >>= 1;
      position++;
    }
    return position - 1;
  }

  int HighestBit64(uint64_t n) const {
    if (n == 0) return -1;

    int position = 0;
    while (n != 0) {
      n >>= 1;
      position++;
    }
    return position - 1;
  }

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

  bool HasFullOrder(ElementType element) const {
    ElementType power = element;

    for (uint64_t k = 2; k < Order(); k++) {
      if (power == 1) {
        return false;
      }
      power = GaloisFieldBinaryExtension<ElementType>::Mul(power, element);
    }

    return power == 1;
  }

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

  ElementType ParsePolynomialString(const std::string &poly_str) const {
    if (poly_str == "0") return 0;

    auto gf2 = std::make_shared<GaloisFieldBinary>();

    auto parsed_poly = utils::ParsePolynomial(gf2, poly_str, variable_name_);

    if (parsed_poly.Degree() >= m_) {
      throw std::invalid_argument(
          "Polynomial degree (" + std::to_string(parsed_poly.Degree()) +
          ") must be less than field extension degree (" + std::to_string(m_) +
          ")");
    }

    return utils::BinaryPolynomialToUint<GaloisFieldBinary>(parsed_poly);
  }
};

template <typename ElementType = uint32_t>
class GFBELogTables : public GaloisFieldBinaryExtension<ElementType> {
 public:

  explicit GFBELogTables(uint8_t m, const std::string &rep = "int",
                         const std::string &irreducible_poly = "",
                         const std::string &variable_name = "α",
                         bool check_irreducible = false,
                         const std::string &generator_name = "g")
      : GaloisFieldBinaryExtension<ElementType>(
            m, rep, irreducible_poly, variable_name, check_irreducible,
            generator_name),
        is_generator_cached_(false) {
    PrecomputeTables();
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    if (a == 0 || b == 0) {
      return 0;
    }

    ElementType log_sum = (log_table_[a] + log_table_[b]);
    if (log_sum >= this->group_order_) {
      log_sum -= this->group_order_;
    }
    return exp_table_[log_sum];
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }
    if (a == 0) {
      return 0;
    }

    ElementType log_diff =
        (log_table_[a] > log_table_[b])
            ? (log_table_[a] - log_table_[b])
            : (log_table_[a] + this->group_order_ - log_table_[b]);

    return exp_table_[log_diff];
  }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }

    return exp_table_[(this->group_order_ - log_table_[a])];
  }

  inline ElementType Pow(const ElementType &a, uint32_t exp) const override {
    if (a == 0) return (exp == 0) ? 1 : 0;

    return exp_table_[(static_cast<uint64_t>(log_table_[a]) * exp) %
                      this->group_order_];
  }

  uint32_t Log(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }

    return log_table_[a];
  }

  uint32_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }

    if (generator == this->MultiplicativeGenerator()) {
      return log_table_[a];
    }

    return GaloisFieldBinaryExtension<ElementType>::Log(a, generator);
  }

  ElementType MultiplicativeGenerator() const override {
    if (is_generator_cached_) {
      return cached_generator_;
    }

    cached_generator_ =
        GaloisFieldBinaryExtension<ElementType>::MultiplicativeGenerator();
    is_generator_cached_ = true;
    return cached_generator_;
  }

 protected:
  std::vector<uint32_t> log_table_;
  std::vector<ElementType> exp_table_;
  mutable ElementType
      cached_generator_;
  mutable bool is_generator_cached_;

  static constexpr ElementType LOG_ZERO =
      std::numeric_limits<ElementType>::max();

  virtual void PrecomputeTables() {

    log_table_.resize(this->order_, 0);
    exp_table_.resize(this->group_order_, 0);

    ElementType generator = this->MultiplicativeGenerator();

    ElementType element = 1;
    for (uint32_t i = 0; i < this->group_order_; i++) {
      exp_table_[i] = element;
      log_table_[element] = i;
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }

    log_table_[0] = LOG_ZERO;
  }
};

template <typename ElementType = uint32_t>
class GFBELogTablesOpt : public GFBELogTables<ElementType> {
 public:

  explicit GFBELogTablesOpt(uint8_t m, const std::string &rep = "int",
                            const std::string &irreducible_poly = "",
                            const std::string &variable_name = "α",
                            bool check_irreducible = false,
                            const std::string &generator_name = "g")
      : GFBELogTables<ElementType>(m, rep, irreducible_poly, variable_name,
                                   check_irreducible, generator_name) {}

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    return this->exp_table_[(this->log_table_[a]) + (this->log_table_[b])];
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }

    return this->exp_table_[this->group_order_ + this->log_table_[a] -
                            this->log_table_[b]];
  }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }

    return this->exp_table_[this->group_order_ - this->log_table_[a]];
  }

 protected:
  void PrecomputeTables() override {

    this->log_table_.resize(this->order_, 0);
    this->exp_table_.resize(4 * this->group_order_ + 1, 0);

    ElementType generator = this->MultiplicativeGenerator();

    ElementType element = 1;
    for (uint32_t i = 0; i < this->group_order_; i++) {
      this->exp_table_[i] = element;
      this->exp_table_[i + this->group_order_] = element;
      this->log_table_[element] = i;
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }

    this->log_table_[0] =
        2 * this->group_order_;
  }
};

template <typename ElementType = uint32_t>
class GFBEZechLogTables : public GaloisFieldBinaryExtension<ElementType> {
 public:

  explicit GFBEZechLogTables(uint8_t m, const std::string &rep = "log",
                             const std::string &irreducible_poly = "",
                             const std::string &variable_name = "α",
                             bool check_irreducible = false,
                             const std::string &generator_name = "g")
      : GaloisFieldBinaryExtension<ElementType>(
            m, rep, irreducible_poly, variable_name, check_irreducible,
            generator_name) {
    PrecomputeTables();
  }

  inline ElementType Add(const ElementType &log_a,
                         const uint32_t &log_b) const override {
    if (log_a == LOG_ZERO) return log_b;
    if (log_b == LOG_ZERO) return log_a;

    uint32_t diff = (log_b >= log_a) ? (log_b - log_a)
                                     : (this->group_order_ - log_a + log_b);

    return (zech_table_[diff] == LOG_ZERO)
               ? LOG_ZERO
               : ((log_a + zech_table_[diff]) <
                  this->group_order_)
                     ? (log_a + zech_table_[diff])
                     : (log_a + zech_table_[diff] - this->group_order_);
  }

  inline ElementType Sub(const ElementType &log_a,
                         const ElementType &log_b) const override {
    return Add(log_a, log_b);
  }

  inline ElementType Mul(const ElementType &log_a,
                         const ElementType &log_b) const override {
    if (log_a == LOG_ZERO || log_b == LOG_ZERO) return LOG_ZERO;

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
    if (log_a == LOG_ZERO) {
      throw std::domain_error("Inverse of zero is undefined");
    }

    return ((log_a == 0) ? 0 : (this->group_order_ - log_a));
  }

  inline ElementType Pow(const ElementType &log_a,
                         uint32_t exp) const override {
    if (log_a == LOG_ZERO) return (exp == 0) ? 0 : LOG_ZERO;

    return (static_cast<uint64_t>(log_a) * exp) % this->group_order_;
  }

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

    ElementType generator = this->MultiplicativeGenerator();
    ElementType element;

    element = 1;
    for (uint32_t i = 0; i < this->group_order_; i++) {
      exp_table_[i] = element;
      log_table[element] = i;
      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }

    element = 1;
    for (uint32_t i = 0; i < this->group_order_; i++) {

      ElementType one_plus_element =
          GaloisFieldBinaryExtension<ElementType>::Add(1, element);
      if (one_plus_element == 0) {
        zech_table_[i] = LOG_ZERO;
      } else {
        zech_table_[i] = log_table[one_plus_element];
      }

      element =
          GaloisFieldBinaryExtension<ElementType>::Mul(element, generator);
    }
  }
};

using GF2 = GaloisFieldBinary;
template <typename ElementType = uint32_t>
using GF2X = GaloisFieldBinaryExtension<ElementType>;
template <typename ElementType = uint32_t>
using GF2XLOG = GFBELogTablesOpt<ElementType>;
using GF2XZECH = GFBEZechLogTables<uint32_t>;

}

#endif
