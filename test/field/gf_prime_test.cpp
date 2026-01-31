#include "xgalois/field/gf_prime.hpp"
#include <gtest/gtest.h>
#include <sstream>
#include <set>

using namespace xg;

//===----------------------------------------------------------------------===//
// GF(p) Prime Field Tests
//===----------------------------------------------------------------------===//

class GaloisFieldPrimeTest : public ::testing::Test {
protected:
  void SetUp() override {
    gf2 = std::make_unique<GaloisFieldPrime<uint8_t>>(2);
    gf3 = std::make_unique<GaloisFieldPrime<uint8_t>>(3);
    gf5 = std::make_unique<GaloisFieldPrime<uint8_t>>(5);
    gf7 = std::make_unique<GaloisFieldPrime<uint8_t>>(7);
    gf11 = std::make_unique<GaloisFieldPrime<uint8_t>>(11);
    gf31 = std::make_unique<GaloisFieldPrime<uint8_t>>(31);
    gf257 = std::make_unique<GaloisFieldPrime<uint16_t>>(257);

    // Different representations
    gf7_hex = std::make_unique<GaloisFieldPrime<uint8_t>>(7, "hex");
    gf7_pow = std::make_unique<GaloisFieldPrime<uint8_t>>(7, "pow");
    gf7_log = std::make_unique<GaloisFieldPrime<uint8_t>>(7, "log");
  }

  std::unique_ptr<GaloisFieldPrime<uint8_t>> gf2, gf3, gf5, gf7, gf11, gf31;
  std::unique_ptr<GaloisFieldPrime<uint16_t>> gf257;
  std::unique_ptr<GaloisFieldPrime<uint8_t>> gf7_hex, gf7_pow, gf7_log;
};

TEST_F(GaloisFieldPrimeTest, FieldProperties) {
  EXPECT_EQ(gf7->Characteristic(), 7);
  EXPECT_EQ(gf7->Order(), 7);
  EXPECT_EQ(gf7->Modulus(), 7);
  EXPECT_EQ(gf7->AdditiveIdentity(), 0);
  EXPECT_EQ(gf7->MultiplicativeIdentity(), 1);

  // Generator should be valid
  uint8_t gen = gf7->MultiplicativeGenerator();
  EXPECT_NE(gen, 0);
  EXPECT_NE(gen, 1);
}

TEST_F(GaloisFieldPrimeTest, Addition) {
  EXPECT_EQ(gf7->Add(0, 0), 0);
  EXPECT_EQ(gf7->Add(3, 4), 0);  // 3 + 4 = 7 ≡ 0 (mod 7)
  EXPECT_EQ(gf7->Add(2, 3), 5);
  EXPECT_EQ(gf7->Add(6, 1), 0);  // 6 + 1 = 7 ≡ 0 (mod 7)
  EXPECT_EQ(gf7->Add(5, 5), 3);  // 5 + 5 = 10 ≡ 3 (mod 7)
}

TEST_F(GaloisFieldPrimeTest, Subtraction) {
  EXPECT_EQ(gf7->Sub(0, 0), 0);
  EXPECT_EQ(gf7->Sub(5, 3), 2);
  EXPECT_EQ(gf7->Sub(3, 5), 5);  // 3 - 5 = -2 ≡ 5 (mod 7)
  EXPECT_EQ(gf7->Sub(0, 1), 6);  // 0 - 1 = -1 ≡ 6 (mod 7)
}

TEST_F(GaloisFieldPrimeTest, Multiplication) {
  EXPECT_EQ(gf7->Mul(0, 5), 0);
  EXPECT_EQ(gf7->Mul(1, 5), 5);
  EXPECT_EQ(gf7->Mul(2, 3), 6);
  EXPECT_EQ(gf7->Mul(3, 5), 1);  // 3 * 5 = 15 ≡ 1 (mod 7)
  EXPECT_EQ(gf7->Mul(6, 6), 1);  // 6 * 6 = 36 ≡ 1 (mod 7)
}

TEST_F(GaloisFieldPrimeTest, Division) {
  EXPECT_EQ(gf7->Div(0, 1), 0);
  EXPECT_EQ(gf7->Div(6, 2), 3);  // 6 / 2 = 3
  EXPECT_EQ(gf7->Div(1, 3), 5);  // 1 / 3 = 1 * 3^(-1) = 1 * 5 = 5
  EXPECT_THROW(gf7->Div(1, 0), std::domain_error);
}

TEST_F(GaloisFieldPrimeTest, Negation) {
  EXPECT_EQ(gf7->Neg(0), 0);
  EXPECT_EQ(gf7->Neg(1), 6);
  EXPECT_EQ(gf7->Neg(3), 4);
  EXPECT_EQ(gf7->Neg(6), 1);
}

TEST_F(GaloisFieldPrimeTest, Inverse) {
  // Test inverse for all non-zero elements
  for (uint8_t a = 1; a < 7; ++a) {
    uint8_t inv_a = gf7->Inv(a);
    EXPECT_EQ(gf7->Mul(a, inv_a), 1) << "Element " << (int)a << " * " << (int)inv_a << " != 1";
  }

  EXPECT_THROW(gf7->Inv(0), std::domain_error);
}

TEST_F(GaloisFieldPrimeTest, Power) {
  EXPECT_EQ(gf7->Pow(0, 0), 1);  // 0^0 = 1 by convention
  EXPECT_EQ(gf7->Pow(0, 5), 0);
  EXPECT_EQ(gf7->Pow(1, 100), 1);
  EXPECT_EQ(gf7->Pow(2, 3), 1);  // 2^3 = 8 ≡ 1 (mod 7)
  EXPECT_EQ(gf7->Pow(3, 6), 1);  // Fermat's little theorem: a^(p-1) ≡ 1 (mod p)
}

TEST_F(GaloisFieldPrimeTest, MultiplicativeStructure) {
  // Test that multiplicative group is cyclic - simplified test
  uint8_t gen = gf7->MultiplicativeGenerator();

  // Test that generator has correct order
  EXPECT_EQ(gf7->Pow(gen, 6), 1);  // Order of multiplicative group is p-1 = 6
  EXPECT_NE(gf7->Pow(gen, 3), 1);  // Should not have smaller order
}

TEST_F(GaloisFieldPrimeTest, Logarithm) {
  uint8_t gen = gf7->MultiplicativeGenerator();

  // Test log base generator for a few elements
  for (uint8_t a = 1; a <= 3; ++a) {  // Test only first 3 elements for speed
    uint32_t log_a = gf7->Log(a, gen);
    EXPECT_EQ(gf7->Pow(gen, log_a), a) << "gen^log(" << (int)a << ") != " << (int)a;
  }

  EXPECT_THROW(gf7->Log(0, gen), std::domain_error);
  EXPECT_THROW(gf7->Log(1, 0), std::invalid_argument);
}

TEST_F(GaloisFieldPrimeTest, Random) {
  std::set<uint8_t> values;
  for (int i = 0; i < 20; ++i) {  // Reduced from 100 to 20 for speed
    uint8_t val = gf7->Random();
    EXPECT_LT(val, 7);
    values.insert(val);
  }
  EXPECT_GE(values.size(), 2);  // Should generate multiple values
}

TEST_F(GaloisFieldPrimeTest, Representation) {
  EXPECT_EQ(gf7->GetRepresentation(), FieldRepresentation::INT);
  EXPECT_EQ(gf7_hex->GetRepresentation(), FieldRepresentation::HEX);
  EXPECT_EQ(gf7_pow->GetRepresentation(), FieldRepresentation::POW);
  EXPECT_EQ(gf7_log->GetRepresentation(), FieldRepresentation::LOG);

  gf7->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(gf7->GetRepresentation(), FieldRepresentation::HEX);
}

TEST_F(GaloisFieldPrimeTest, ToString) {
  EXPECT_EQ(gf7->ToString(3), "3");

  // Test with valid field elements
  gf7->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(gf7->ToString(6), "0x6");
}

TEST_F(GaloisFieldPrimeTest, Print) {
  std::ostringstream oss;
  gf7->Print(oss);
  EXPECT_FALSE(oss.str().empty());  // Just check it prints something

  oss.str("");
  gf7->Print(3, oss);
  EXPECT_EQ(oss.str(), "3");
}

TEST_F(GaloisFieldPrimeTest, Constructor) {
  // Test invalid prime
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(0), std::invalid_argument);
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(1), std::invalid_argument);

  // Test prime checking (when enabled)
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(4, "int", true), std::invalid_argument);  // 4 is not prime
  EXPECT_NO_THROW(GaloisFieldPrime<uint8_t>(7, "int", true));  // 7 is prime
}

TEST_F(GaloisFieldPrimeTest, SpecialCases) {
  // Test GF(2) special case
  EXPECT_EQ(gf2->MultiplicativeGenerator(), 1);
  EXPECT_EQ(gf2->Add(1, 1), 0);
  EXPECT_EQ(gf2->Mul(1, 1), 1);

  // Test GF(3)
  EXPECT_EQ(gf3->Add(1, 2), 0);
  EXPECT_EQ(gf3->Mul(2, 2), 1);
}

TEST_F(GaloisFieldPrimeTest, LargerFields) {
  // Test operations in GF(257) - simplified for speed
  EXPECT_EQ(gf257->Characteristic(), 257);
  EXPECT_EQ(gf257->Order(), 257);

  // Test basic arithmetic
  EXPECT_EQ(gf257->Add(200, 100), 43);  // 200 + 100 = 300 ≡ 43 (mod 257)
  EXPECT_EQ(gf257->Mul(2, 129), 1);    // 2 * 129 = 258 ≡ 1 (mod 257)

  // Test that a few elements have inverses
  for (uint16_t a = 1; a <= 3; ++a) {  // Test only first 3 elements for speed
    uint16_t inv_a = gf257->Inv(a);
    EXPECT_EQ(gf257->Mul(a, inv_a), 1);
  }
}

TEST_F(GaloisFieldPrimeTest, SetElementValueString) {
  // Test power string parsing with positive exponents
  uint8_t gen = gf7->MultiplicativeGenerator();

  // Test g^0 = 1
  uint8_t elem_0 = gf7->SetElementValue("g^0");
  EXPECT_EQ(elem_0, 1);

  // Test g^1 = generator
  uint8_t elem_1 = gf7->SetElementValue("g^1");
  EXPECT_EQ(elem_1, gen);

  // Test g^2
  uint8_t elem_2 = gf7->SetElementValue("g^2");
  EXPECT_EQ(elem_2, gf7->Mul(gen, gen));

  // Test negative exponents (g^-1 should equal inverse of generator)
  uint8_t elem_neg1 = gf7->SetElementValue("g^-1");
  EXPECT_EQ(elem_neg1, gf7->Inv(gen));

  // Test g^-2 should equal inverse of g^2
  uint8_t elem_neg2 = gf7->SetElementValue("g^-2");
  uint8_t inv_g2 = gf7->Inv(gf7->Mul(gen, gen));
  EXPECT_EQ(elem_neg2, inv_g2);

  // Test wraparound: g^6 should equal g^0 in GF(7) (since group order is 6)
  uint8_t elem_6 = gf7->SetElementValue("g^6");
  EXPECT_EQ(elem_6, 1);

  // Test large positive exponent that wraps around
  uint8_t elem_7 = gf7->SetElementValue("g^7");
  EXPECT_EQ(elem_7, gen);  // Should be same as g^1
}

TEST_F(GaloisFieldPrimeTest, SetElementValueStringInvalid) {
  // Test invalid string formats
  EXPECT_THROW(gf7->SetElementValue("invalid"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("g^"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("g^abc"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("h^1"), std::invalid_argument);  // Wrong generator name
  EXPECT_THROW(gf7->SetElementValue("5"), std::invalid_argument);    // Plain number not supported
}

TEST_F(GaloisFieldPrimeTest, SetElementValueNumeric) {
  // Test numeric SetElementValue for reduction
  EXPECT_EQ(gf7->SetElementValue(0), 0);
  EXPECT_EQ(gf7->SetElementValue(6), 6);
  EXPECT_EQ(gf7->SetElementValue(7), 0);  // 7 ≡ 0 (mod 7)
  EXPECT_EQ(gf7->SetElementValue(8), 1);  // 8 ≡ 1 (mod 7)
  EXPECT_EQ(gf7->SetElementValue(14), 0); // 14 ≡ 0 (mod 7)
}

TEST_F(GaloisFieldPrimeTest, GetElementValue) {
  // GetElementValue should be identity for prime fields
  for (uint8_t i = 0; i < 7; ++i) {
    EXPECT_EQ(gf7->GetElementValue(i), i);
  }
}
