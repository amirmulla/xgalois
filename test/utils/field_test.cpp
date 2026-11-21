#include "xgalois/utils/field.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace xg::utils;

namespace {

// Minimal mock prime field GF(p) implementing GaloisFieldBase<uint32_t> for
// deterministic unit testing of ParsePowerString without pulling in gf_prime.
class MockPrimeField : public xg::GaloisFieldBase<uint32_t> {
 public:
  MockPrimeField(uint32_t prime, uint32_t generator)
      : p_(prime), gen_(generator) {}

  uint64_t Characteristic() const override { return p_; }
  uint64_t Order() const override { return p_; }

  uint32_t Add(const uint32_t& a, const uint32_t& b) const override {
    return (a + b) % p_;
  }
  uint32_t Sub(const uint32_t& a, const uint32_t& b) const override {
    return (a + p_ - (b % p_)) % p_;
  }
  uint32_t Mul(const uint32_t& a, const uint32_t& b) const override {
    return static_cast<uint32_t>(
        (static_cast<uint64_t>(a) * static_cast<uint64_t>(b)) % p_);
  }
  uint32_t Div(const uint32_t& a, const uint32_t& b) const override {
    return Mul(a, Inv(b));
  }
  uint32_t Neg(const uint32_t& a) const override {
    return (a % p_ == 0) ? 0 : (p_ - (a % p_));
  }
  uint32_t Inv(const uint32_t& a) const override {
    if (a % p_ == 0) {
      throw std::domain_error("Division by zero");
    }
    return Pow(a, p_ - 2);
  }
  uint32_t Pow(const uint32_t& a, uint64_t exp) const override {
    uint64_t result = 1;
    uint64_t base = a % p_;
    while (exp > 0) {
      if (exp & 1ULL) {
        result = (result * base) % p_;
      }
      base = (base * base) % p_;
      exp >>= 1;
    }
    return static_cast<uint32_t>(result);
  }
  uint32_t Sqrt(const uint32_t& a) const override { return a; }
  uint64_t Log(const uint32_t&) const override { return 0; }
  uint64_t Log(const uint32_t&, const uint32_t&) const override { return 0; }
  uint32_t Random() const override { return 1; }
  uint32_t MultiplicativeIdentity() const override { return 1; }
  uint32_t AdditiveIdentity() const override { return 0; }
  uint32_t MultiplicativeGenerator() const override { return gen_; }
  std::vector<uint32_t> MultiplicativeGenerators() const override {
    return {gen_};
  }
  uint32_t GetElementValue(const uint32_t& value) const override {
    return value;
  }
  uint32_t SetElementValue(const uint32_t& value) const override {
    return value % p_;
  }
  uint32_t SetElementValue(const std::string&) const override { return 0; }
  void Print(std::ostream& os) const override { os << "MockGF(" << p_ << ")"; }
  void Print(const uint32_t& a, std::ostream& os) const override { os << a; }
  std::string ToString(const uint32_t& a) const override {
    return std::to_string(a);
  }
  xg::FieldRepresentation GetRepresentation() const override { return rep_; }
  void SetRepresentation(xg::FieldRepresentation rep) override { rep_ = rep; }

 private:
  uint32_t p_;
  uint32_t gen_;
  xg::FieldRepresentation rep_ = xg::FieldRepresentation::INT;
};

}  // namespace

class FieldUtilsTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(FieldUtilsTest, ConvertRepresentationValidCases) {
  EXPECT_EQ(ConvertRepresentation("int"), xg::FieldRepresentation::INT);
  EXPECT_EQ(ConvertRepresentation("hex"), xg::FieldRepresentation::HEX);
  EXPECT_EQ(ConvertRepresentation("pow"), xg::FieldRepresentation::POW);
  EXPECT_EQ(ConvertRepresentation("log"), xg::FieldRepresentation::LOG);
  EXPECT_EQ(ConvertRepresentation("poly"), xg::FieldRepresentation::POLY);
}

TEST_F(FieldUtilsTest, ConvertRepresentationCaseSensitive) {
  EXPECT_THROW(ConvertRepresentation("INT"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("HEX"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("POW"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("LOG"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("POLY"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ConvertRepresentationInvalidCases) {
  EXPECT_THROW(ConvertRepresentation("invalid"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation(""), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("integer"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("hexadecimal"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("power"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("logarithm"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("polynomial"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ConvertRepresentationExceptionMessage) {
  try {
    ConvertRepresentation("unknown");
    FAIL() << "Expected std::invalid_argument to be thrown";
  } catch (const std::invalid_argument& e) {
    EXPECT_STREQ(e.what(), "Unknown representation: unknown");
  }
}

TEST_F(FieldUtilsTest, ConvertRepresentationWithWhitespace) {
  EXPECT_THROW(ConvertRepresentation(" int"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("int "), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation(" int "), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("i nt"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ConvertRepresentationSpecialCharacters) {
  EXPECT_THROW(ConvertRepresentation("int\n"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("int\t"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("int\r"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ParsePowerStringValidPowers) {
  // GF(7) with primitive root g = 3:
  // g^0=1, g^1=3, g^2=2, g^3=6, g^4=4, g^5=5, g^6=1
  MockPrimeField gf7(7, 3);

  EXPECT_EQ(ParsePowerString("g^0", gf7), 1U);
  EXPECT_EQ(ParsePowerString("g^1", gf7), 3U);
  EXPECT_EQ(ParsePowerString("g^2", gf7), 2U);
  EXPECT_EQ(ParsePowerString("g^3", gf7), 6U);
  EXPECT_EQ(ParsePowerString("g^4", gf7), 4U);
  EXPECT_EQ(ParsePowerString("g^5", gf7), 5U);
  EXPECT_EQ(ParsePowerString("g^6", gf7), 1U);
  EXPECT_EQ(ParsePowerString("g^7", gf7), 3U);

  // Valid custom identifiers and UTF-8 Greek symbols.
  EXPECT_EQ(ParsePowerString("alpha^2", gf7), 2U);
  EXPECT_EQ(ParsePowerString("g_1^3", gf7), 6U);
  EXPECT_EQ(ParsePowerString("α^4", gf7), 4U);

  // Optional expected_generator matching.
  EXPECT_EQ(ParsePowerString("g^2", gf7, "g"), 2U);
  EXPECT_THROW(ParsePowerString("h^2", gf7, "g"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ParsePowerString64BitExponentReduction) {
  // GF(7) has multiplicative group order 6.
  // 2^32 + 1 = 4294967297.
  // Truncating to uint32_t would yield 1 (giving g^1 = 3),
  // whereas exact 64-bit reduction yields 4294967297 % 6 = 5 (giving g^5 = 5).
  MockPrimeField gf7(7, 3);
  EXPECT_EQ(ParsePowerString("g^4294967297", gf7), 5U);

  // UINT64_MAX = 18446744073709551615; UINT64_MAX % 6 = 3 => g^3 = 6.
  EXPECT_EQ(ParsePowerString("g^18446744073709551615", gf7), 6U);

  // GF(2) edge case (group order = 1).
  MockPrimeField gf2(2, 1);
  EXPECT_EQ(ParsePowerString("g^0", gf2), 1U);
  EXPECT_EQ(ParsePowerString("g^123456789012345", gf2), 1U);
}

TEST_F(FieldUtilsTest, ParsePowerStringNegativeExponents) {
  // GF(7) with generator g = 3:
  // g^-1 = g^5 = 5, g^-2 = g^4 = 4, g^-6 = g^0 = 1
  MockPrimeField gf7(7, 3);
  EXPECT_EQ(ParsePowerString("g^-0", gf7), 1U);
  EXPECT_EQ(ParsePowerString("g^-1", gf7), gf7.Inv(3U));
  EXPECT_EQ(ParsePowerString("g^-2", gf7), gf7.Inv(gf7.Mul(3U, 3U)));
  EXPECT_EQ(ParsePowerString("g^-6", gf7), 1U);

  // Large negative 64-bit exponent: -(2^32 + 1) mod 6 = -5 mod 6 = 1 => g^1
  // = 3.
  EXPECT_EQ(ParsePowerString("g^-4294967297", gf7), 3U);

  // Malformed negative exponents must be rejected cleanly.
  EXPECT_THROW(ParsePowerString("g^-", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^--1", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^-+1", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^- 1", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^-1abc", gf7), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ParsePowerStringInvalidBasesAndTrailingGarbage) {
  MockPrimeField gf7(7, 3);

  // Missing caret or invalid/empty/numeric bases.
  EXPECT_THROW(ParsePowerString("", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("^2", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("0^2", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("123^2", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString(" g^2", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g ^2", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g+1^2", gf7), std::invalid_argument);

  // Malformed exponents and trailing garbage.
  EXPECT_THROW(ParsePowerString("g^", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^abc", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^3abc", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^ 3", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^3 ", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^+3", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^1^2", gf7), std::invalid_argument);
  EXPECT_THROW(ParsePowerString("g^999999999999999999999999999999", gf7),
               std::invalid_argument);
}
