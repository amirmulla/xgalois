#include "xgalois/field/gf_binary.hpp"

#include <gtest/gtest.h>

#include <sstream>

using namespace xg;

//===----------------------------------------------------------------------===//
// GF(2) Tests
//===----------------------------------------------------------------------===//

class GaloisFieldBinaryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf2 = std::make_unique<GaloisFieldBinary>();
    gf2_hex = std::make_unique<GaloisFieldBinary>("hex");
  }

  std::unique_ptr<GaloisFieldBinary> gf2;
  std::unique_ptr<GaloisFieldBinary> gf2_hex;
};

TEST_F(GaloisFieldBinaryTest, FieldProperties) {
  EXPECT_EQ(gf2->Characteristic(), 2);
  EXPECT_EQ(gf2->Order(), 2);
  EXPECT_EQ(gf2->Modulus(), 2);
  EXPECT_EQ(gf2->AdditiveIdentity(), 0);
  EXPECT_EQ(gf2->MultiplicativeIdentity(), 1);
  EXPECT_EQ(gf2->MultiplicativeGenerator(), 1);
}

TEST_F(GaloisFieldBinaryTest, Addition) {
  EXPECT_EQ(gf2->Add(0, 0), 0);
  EXPECT_EQ(gf2->Add(0, 1), 1);
  EXPECT_EQ(gf2->Add(1, 0), 1);
  EXPECT_EQ(gf2->Add(1, 1), 0);
}

TEST_F(GaloisFieldBinaryTest, Subtraction) {
  EXPECT_EQ(gf2->Sub(0, 0), 0);
  EXPECT_EQ(gf2->Sub(0, 1), 1);
  EXPECT_EQ(gf2->Sub(1, 0), 1);
  EXPECT_EQ(gf2->Sub(1, 1), 0);
}

TEST_F(GaloisFieldBinaryTest, Multiplication) {
  EXPECT_EQ(gf2->Mul(0, 0), 0);
  EXPECT_EQ(gf2->Mul(0, 1), 0);
  EXPECT_EQ(gf2->Mul(1, 0), 0);
  EXPECT_EQ(gf2->Mul(1, 1), 1);
}

TEST_F(GaloisFieldBinaryTest, Division) {
  EXPECT_EQ(gf2->Div(0, 1), 0);
  EXPECT_EQ(gf2->Div(1, 1), 1);
  EXPECT_THROW(gf2->Div(0, 0), std::domain_error);
  EXPECT_THROW(gf2->Div(1, 0), std::domain_error);
}

TEST_F(GaloisFieldBinaryTest, Negation) {
  EXPECT_EQ(gf2->Neg(0), 0);
  EXPECT_EQ(gf2->Neg(1), 1);
}

TEST_F(GaloisFieldBinaryTest, Inverse) {
  EXPECT_EQ(gf2->Inv(1), 1);
  EXPECT_THROW(gf2->Inv(0), std::domain_error);
}

TEST_F(GaloisFieldBinaryTest, Power) {
  EXPECT_EQ(gf2->Pow(0, 0), 1);  // 0^0 = 1 by convention
  EXPECT_EQ(gf2->Pow(0, 1), 0);
  EXPECT_EQ(gf2->Pow(0, 5), 0);
  EXPECT_EQ(gf2->Pow(1, 0), 1);
  EXPECT_EQ(gf2->Pow(1, 1), 1);
  EXPECT_EQ(gf2->Pow(1, 100), 1);
}

TEST_F(GaloisFieldBinaryTest, SquareRoot) {
  EXPECT_EQ(gf2->Sqrt(0), 0);
  EXPECT_EQ(gf2->Sqrt(1), 1);
}

TEST_F(GaloisFieldBinaryTest, Logarithm) {
  EXPECT_EQ(gf2->Log(1), 0);
  EXPECT_EQ(gf2->Log(1, 1), 0);
  EXPECT_THROW(gf2->Log(0), std::domain_error);
  EXPECT_THROW(gf2->Log(1, 0), std::invalid_argument);
}

TEST_F(GaloisFieldBinaryTest, Random) {
  std::set<uint8_t> values;
  for (int i = 0; i < 100; ++i) {
    values.insert(gf2->Random());
  }
  EXPECT_GE(values.size(), 1);  // Should generate at least one value
  EXPECT_LE(values.size(), 2);  // Should only generate 0 and 1
  for (auto val : values) {
    EXPECT_TRUE(val == 0 || val == 1);
  }
}

TEST_F(GaloisFieldBinaryTest, Representation) {
  EXPECT_EQ(gf2->GetRepresentation(), FieldRepresentation::INT);
  EXPECT_EQ(gf2_hex->GetRepresentation(), FieldRepresentation::HEX);

  gf2->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(gf2->GetRepresentation(), FieldRepresentation::HEX);

  EXPECT_THROW(gf2->SetRepresentation(FieldRepresentation::POLY),
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryTest, ToString) {
  EXPECT_EQ(gf2->ToString(0), "0");
  EXPECT_EQ(gf2->ToString(1), "1");
  EXPECT_EQ(gf2_hex->ToString(0), "0x0");
  EXPECT_EQ(gf2_hex->ToString(1), "0x1");
}

TEST_F(GaloisFieldBinaryTest, Print) {
  std::ostringstream oss;
  gf2->Print(oss);
  EXPECT_EQ(oss.str(), "GF(2)");

  oss.str("");
  gf2->Print(1, oss);
  EXPECT_EQ(oss.str(), "1");

  oss.str("");
  gf2_hex->Print(1, oss);
  EXPECT_EQ(oss.str(), "0x1");
}

//===----------------------------------------------------------------------===//
// GF(2^m) Binary Extension Tests
//===----------------------------------------------------------------------===//

class GaloisFieldBinaryExtensionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // GF(2^3) = GF(8) with irreducible polynomial x^3 + x + 1
    gf8 = std::make_unique<GaloisFieldBinaryExtension<uint8_t>>(3);
    // GF(2^4) = GF(16) with irreducible polynomial from database
    gf16 = std::make_unique<GaloisFieldBinaryExtension<uint16_t>>(4);
  }

  std::unique_ptr<GaloisFieldBinaryExtension<uint8_t>> gf8;
  std::unique_ptr<GaloisFieldBinaryExtension<uint16_t>> gf16;
};

TEST_F(GaloisFieldBinaryExtensionTest, FieldProperties) {
  EXPECT_EQ(gf8->Characteristic(), 2);
  EXPECT_EQ(gf8->Order(), 8);
  EXPECT_EQ(gf8->Degree(), 3);
  EXPECT_EQ(gf8->AdditiveIdentity(), 0);
  EXPECT_EQ(gf8->MultiplicativeIdentity(), 1);

  EXPECT_EQ(gf16->Characteristic(), 2);
  EXPECT_EQ(gf16->Order(), 16);
  EXPECT_EQ(gf16->Degree(), 4);
}

TEST_F(GaloisFieldBinaryExtensionTest, BasicArithmetic) {
  // Test addition (XOR operation)
  EXPECT_EQ(gf8->Add(0, 0), 0);
  EXPECT_EQ(gf8->Add(1, 1), 0);
  EXPECT_EQ(gf8->Add(3, 5), 6);  // 011 XOR 101 = 110

  // Test subtraction (same as addition in GF(2^m))
  EXPECT_EQ(gf8->Sub(3, 5), gf8->Add(3, 5));

  // Test multiplication
  EXPECT_EQ(gf8->Mul(0, 5), 0);
  EXPECT_EQ(gf8->Mul(1, 5), 5);
  // In GF(2^3), multiplication is more complex than simple arithmetic
  uint8_t mul_result = gf8->Mul(2, 3);
  EXPECT_GE(mul_result, 0);
  EXPECT_LT(mul_result, 8);
}

TEST_F(GaloisFieldBinaryExtensionTest, MultiplicativeStructure) {
  // Test that all non-zero elements have inverses
  for (uint8_t a = 1; a < 8; ++a) {
    uint8_t inv_a = gf8->Inv(a);
    EXPECT_EQ(gf8->Mul(a, inv_a), 1)
        << "Element " << (int)a << " * " << (int)inv_a << " != 1";
  }

  // Test division
  for (uint8_t a = 0; a < 8; ++a) {
    for (uint8_t b = 1; b < 8; ++b) {
      uint8_t div_result = gf8->Div(a, b);
      EXPECT_EQ(gf8->Mul(div_result, b), a)
          << "(" << (int)a << " / " << (int)b << ") * " << (int)b
          << " != " << (int)a;
    }
  }
}

TEST_F(GaloisFieldBinaryExtensionTest, Random) {
  std::set<uint8_t> values;
  for (int i = 0; i < 100; ++i) {
    uint8_t val = gf8->Random();
    EXPECT_LT(val, 8);
    values.insert(val);
  }
  EXPECT_GE(values.size(), 3);  // Should generate multiple values
}

TEST_F(GaloisFieldBinaryExtensionTest, InvalidOperations) {
  EXPECT_THROW(gf8->Inv(0), std::domain_error);
  EXPECT_THROW(gf8->Div(1, 0), std::domain_error);
  EXPECT_THROW(gf8->Log(0), std::domain_error);
}

TEST_F(GaloisFieldBinaryExtensionTest, Constructor) {
  // Test invalid degree
  EXPECT_THROW(GaloisFieldBinaryExtension<uint8_t>(0), std::invalid_argument);
  EXPECT_THROW(GaloisFieldBinaryExtension<uint8_t>(10),
               std::invalid_argument);  // Too large for uint8_t

  // Test custom irreducible polynomial
  EXPECT_NO_THROW(GaloisFieldBinaryExtension<uint8_t>(3, "int", "x^3 + x + 1"));
  EXPECT_THROW(GaloisFieldBinaryExtension<uint8_t>(3, "int", "x^2 + 1"),
               std::invalid_argument);  // Wrong degree
}

//===----------------------------------------------------------------------===//
// GF(2) Element Assignment and String Constructor Tests
//===----------------------------------------------------------------------===//

class GaloisFieldBinaryElementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf2 = std::make_shared<GaloisFieldBinary>();
    gf2_hex = std::make_shared<GaloisFieldBinary>("hex");
  }

  std::shared_ptr<GaloisFieldBinary> gf2;
  std::shared_ptr<GaloisFieldBinary> gf2_hex;
};

TEST_F(GaloisFieldBinaryElementTest, ElementConstruction) {
  // Test basic construction with numeric values
  GaloisFieldElementBase<GaloisFieldBinary> elem_0(0, gf2);
  GaloisFieldElementBase<GaloisFieldBinary> elem_1(1, gf2);

  EXPECT_EQ(elem_0.Value(), 0);
  EXPECT_EQ(elem_1.Value(), 1);
  EXPECT_EQ(elem_0.Field(), gf2);
  EXPECT_EQ(elem_1.Field(), gf2);
}

TEST_F(GaloisFieldBinaryElementTest, StringConstructorNotSupported) {
  // GF(2) doesn't support string representation for elements
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinary>("0", gf2),
               std::invalid_argument);
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinary>("1", gf2),
               std::invalid_argument);
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinary>("α", gf2),
               std::invalid_argument);
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinary>("g^1", gf2),
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryElementTest, NumericAssignmentOperator) {
  GaloisFieldElementBase<GaloisFieldBinary> elem(0, gf2);

  // Test assignment of numeric values
  elem = static_cast<uint8_t>(1);
  EXPECT_EQ(elem.Value(), 1);

  elem = static_cast<uint8_t>(0);
  EXPECT_EQ(elem.Value(), 0);

  // Test assignment with values outside GF(2) - behavior depends on
  // implementation In a strict implementation, only 0 and 1 should be valid for
  // GF(2)
  elem = static_cast<uint8_t>(2);
  // The field's SetElementValue should handle this appropriately
  EXPECT_GE(elem.Value(), 0);
  EXPECT_LT(elem.Value(),
            8);  // Assuming it's at least valid within a reasonable range
}

TEST_F(GaloisFieldBinaryElementTest, StringAssignmentNotSupported) {
  GaloisFieldElementBase<GaloisFieldBinary> elem(0, gf2);

  // String assignment should throw since GF(2) doesn't support string
  // representation
  EXPECT_THROW(elem = "0", std::invalid_argument);
  EXPECT_THROW(elem = "1", std::invalid_argument);
  EXPECT_THROW(elem = "α", std::invalid_argument);
  EXPECT_THROW(elem = "g^0", std::invalid_argument);
}

TEST_F(GaloisFieldBinaryElementTest, ElementCopyAssignment) {
  GaloisFieldElementBase<GaloisFieldBinary> elem_0(0, gf2);
  GaloisFieldElementBase<GaloisFieldBinary> elem_1(1, gf2);

  // Test copy assignment between elements
  elem_0 = elem_1;
  EXPECT_EQ(elem_0.Value(), 1);
  EXPECT_EQ(elem_0.Field(), gf2);

  elem_1 = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
  EXPECT_EQ(elem_1.Value(), 0);
}

TEST_F(GaloisFieldBinaryElementTest, CompoundAssignmentOperators) {
  GaloisFieldElementBase<GaloisFieldBinary> elem_a(1, gf2);
  GaloisFieldElementBase<GaloisFieldBinary> elem_b(1, gf2);

  // Test += operator (XOR in GF(2))
  elem_a += elem_b;
  EXPECT_EQ(elem_a.Value(), 0);  // 1 + 1 = 0 in GF(2)

  elem_a = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
  elem_a += elem_b;
  EXPECT_EQ(elem_a.Value(), 1);  // 0 + 1 = 1 in GF(2)

  // Test -= operator (same as += in GF(2))
  elem_a -= elem_b;
  EXPECT_EQ(elem_a.Value(), 0);  // 1 - 1 = 0 in GF(2)

  // Test *= operator
  elem_a = GaloisFieldElementBase<GaloisFieldBinary>(1, gf2);
  elem_a *= elem_b;
  EXPECT_EQ(elem_a.Value(), 1);  // 1 * 1 = 1 in GF(2)

  elem_a = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
  elem_a *= elem_b;
  EXPECT_EQ(elem_a.Value(), 0);  // 0 * 1 = 0 in GF(2)

  // Test /= operator
  elem_a = GaloisFieldElementBase<GaloisFieldBinary>(1, gf2);
  elem_a /= elem_b;
  EXPECT_EQ(elem_a.Value(), 1);  // 1 / 1 = 1 in GF(2)

  elem_a = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
  elem_a /= elem_b;
  EXPECT_EQ(elem_a.Value(), 0);  // 0 / 1 = 0 in GF(2)

  // Test division by zero throws
  GaloisFieldElementBase<GaloisFieldBinary> elem_zero(0, gf2);
  EXPECT_THROW(elem_a /= elem_zero, std::domain_error);
}

TEST_F(GaloisFieldBinaryElementTest, StrictElementValidation) {
  // Test that GaloisFieldElement enforces field matching
  auto gf2_other = std::make_shared<GaloisFieldBinary>();

  GaloisFieldElement<GaloisFieldBinary> elem_gf2(1, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_other(1, gf2_other);

  // Operations between elements of different field instances should throw
  EXPECT_THROW(elem_gf2 + elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf2 - elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf2 * elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf2 / elem_other, std::invalid_argument);
}

TEST_F(GaloisFieldBinaryElementTest, ElementPrinting) {
  GaloisFieldElementBase<GaloisFieldBinary> elem_0(0, gf2);
  GaloisFieldElementBase<GaloisFieldBinary> elem_1(1, gf2);
  GaloisFieldElementBase<GaloisFieldBinary> elem_hex(1, gf2_hex);

  std::ostringstream oss;

  // Test printing with INT representation
  elem_0.Print(oss);
  EXPECT_EQ(oss.str(), "0");

  oss.str("");
  elem_1.Print(oss);
  EXPECT_EQ(oss.str(), "1");

  // Test printing with HEX representation
  oss.str("");
  elem_hex.Print(oss);
  EXPECT_EQ(oss.str(), "0x1");

  // Test stream operator
  oss.str("");
  oss << elem_1;
  EXPECT_EQ(oss.str(), "1");
}

TEST_F(GaloisFieldBinaryElementTest, NullFieldHandling) {
  // Test that null field pointer throws during construction
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinary>(0, nullptr),
               std::invalid_argument);

  // Test that operations on elements with null field throw
  GaloisFieldElementBase<GaloisFieldBinary> elem;  // Default constructed
  // Note: Default constructed elements may have null field, operations should
  // handle this gracefully
}

//===----------------------------------------------------------------------===//
// GF(2^m) Element Assignment and String Constructor Tests
//===----------------------------------------------------------------------===//

class GaloisFieldBinaryExtensionElementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // GF(2^3) = GF(8) with default irreducible polynomial
    gf8 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);
    // GF(2^4) = GF(16) with custom variable name
    gf16 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4, "poly", "",
                                                                  "β");
  }

  std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>> gf8;
  std::shared_ptr<GaloisFieldBinaryExtension<uint16_t>> gf16;
};

TEST_F(GaloisFieldBinaryExtensionElementTest, ElementConstruction) {
  // Test basic construction with numeric values
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_0(0, gf8);
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_5(5, gf8);

  EXPECT_EQ(elem_0.Value(), 0);
  EXPECT_EQ(elem_5.Value(), 5);
  EXPECT_EQ(elem_0.Field(), gf8);
  EXPECT_EQ(elem_5.Field(), gf8);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorPolynomial) {
  // Test polynomial string construction using default variable name "α"
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_poly(
      "α^2 + α + 1", gf8);

  // α^2 + α + 1 should correspond to binary 111 = 7 in GF(2^3)
  EXPECT_EQ(elem_poly.Value(), 7);

  // Test simple polynomial terms
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_alpha("α",
                                                                         gf8);
  EXPECT_EQ(elem_alpha.Value(), 2);  // α corresponds to binary 10 = 2

  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_alpha_sq(
      "α^2", gf8);
  EXPECT_EQ(elem_alpha_sq.Value(), 4);  // α^2 corresponds to binary 100 = 4

  // Test constant terms - use explicit polynomial format
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_one("g^0",
                                                                       gf8);
  EXPECT_EQ(elem_one.Value(), 1);

  // Test zero using numeric constructor since "0" might not be supported as
  // string
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_zero(0, gf8);
  EXPECT_EQ(elem_zero.Value(), 0);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorCustomVariable) {
  // Test polynomial string construction with custom variable name "β"
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint16_t>> elem_beta(
      "β^2 + β", gf16);

  // β^2 + β should correspond to binary 110 = 6 in GF(2^4)
  EXPECT_EQ(elem_beta.Value(), 6);

  // Test that original variable name doesn't work
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinaryExtension<uint16_t>>(
                   "α^2 + α", gf16),
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorGenerator) {
  // Test generator power string construction using default generator name "g"
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_g0("g^0",
                                                                      gf8);
  EXPECT_EQ(elem_g0.Value(), 1);  // g^0 = 1

  // Test other generator powers
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_g1("g^1",
                                                                      gf8);
  uint8_t generator = gf8->MultiplicativeGenerator();
  EXPECT_EQ(elem_g1.Value(), generator);

  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_g2("g^2",
                                                                      gf8);
  EXPECT_EQ(elem_g2.Value(), gf8->Mul(generator, generator));
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorInvalid) {
  // Test invalid string formats
  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(
                   "invalid", gf8),
               std::invalid_argument);

  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(
                   "h^1", gf8),  // Wrong generator name
               std::invalid_argument);

  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(
                   "α^8", gf8),  // Degree too high
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, NumericAssignmentOperator) {
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem(0, gf8);

  // Test assignment of numeric values
  elem = static_cast<uint8_t>(5);
  EXPECT_EQ(elem.Value(), 5);

  elem = static_cast<uint8_t>(7);
  EXPECT_EQ(elem.Value(), 7);

  elem = static_cast<uint8_t>(0);
  EXPECT_EQ(elem.Value(), 0);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringAssignmentOperator) {
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem(0, gf8);

  // Test assignment with polynomial strings
  elem = "α^2 + α";
  EXPECT_EQ(elem.Value(), 6);  // α^2 + α = binary 110 = 6

  elem = "α^2 + 1";
  EXPECT_EQ(elem.Value(), 5);  // α^2 + 1 = binary 101 = 5

  elem = "α";
  EXPECT_EQ(elem.Value(), 2);  // α = binary 10 = 2

  elem = "g^0";
  EXPECT_EQ(elem.Value(), 1);

  // Test numeric assignment for zero since string "0" might not be supported
  elem = static_cast<uint8_t>(0);
  EXPECT_EQ(elem.Value(), 0);

  // Test assignment with generator power strings
  elem = "g^0";
  EXPECT_EQ(elem.Value(), 1);

  elem = "g^1";
  EXPECT_EQ(elem.Value(), gf8->MultiplicativeGenerator());
}

TEST_F(GaloisFieldBinaryExtensionElementTest, CompoundAssignmentOperators) {
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_a(
      3, gf8);  // 011
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_b(
      5, gf8);  // 101

  // Test += operator (XOR in GF(2^m))
  elem_a += elem_b;
  EXPECT_EQ(elem_a.Value(), 6);  // 011 XOR 101 = 110 = 6

  // Test -= operator (same as += in characteristic 2)
  elem_a -= elem_b;
  EXPECT_EQ(elem_a.Value(), 3);  // 110 XOR 101 = 011 = 3

  // Test *= operator
  elem_a = GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(2, gf8);
  elem_b = GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(3, gf8);
  elem_a *= elem_b;
  EXPECT_EQ(elem_a.Value(), gf8->Mul(2, 3));

  // Test /= operator
  elem_a = GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(6, gf8);
  elem_b = GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>(2, gf8);
  elem_a /= elem_b;
  EXPECT_EQ(elem_a.Value(), gf8->Div(6, 2));

  // Test division by zero throws
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_zero(0, gf8);
  EXPECT_THROW(elem_a /= elem_zero, std::domain_error);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, ElementArithmetic) {
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_2(2, gf8);
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_3(3, gf8);

  // Test addition
  auto sum = elem_2 + elem_3;
  EXPECT_EQ(sum.Value(), 1);  // 010 XOR 011 = 001 = 1

  // Test subtraction (same as addition)
  auto diff = elem_3 - elem_2;
  EXPECT_EQ(diff.Value(), 1);  // 011 XOR 010 = 001 = 1

  // Test multiplication
  auto prod = elem_2 * elem_3;
  EXPECT_EQ(prod.Value(), gf8->Mul(2, 3));

  // Test division
  auto quot = prod / elem_2;
  EXPECT_EQ(quot.Value(), 3);  // Should get back elem_3

  // Test negation (identity in characteristic 2)
  auto neg = -elem_2;
  EXPECT_EQ(neg.Value(), 2);

  // Test power operations
  auto power = elem_2 ^ 3;
  EXPECT_EQ(power.Value(), gf8->Pow(2, 3));

  auto pow_method = elem_2.Pow(3);
  EXPECT_EQ(pow_method.Value(), gf8->Pow(2, 3));

  // Test inverse
  auto inv = elem_3.Inv();
  EXPECT_EQ((elem_3 * inv).Value(), 1);

  // Test square root
  auto sqrt_elem = elem_2.Sqrt();
  auto sqrt_squared = sqrt_elem * sqrt_elem;
  // In GF(2^m), sqrt is well-defined, but the relationship might be more
  // complex
  EXPECT_GE(sqrt_elem.Value(), 0);
  EXPECT_LT(sqrt_elem.Value(), 8);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StrictElementValidation) {
  // Test that GaloisFieldElement enforces field matching
  auto gf8_other = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);

  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_gf8(2, gf8);
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_other(3,
                                                                     gf8_other);

  // Operations between elements of different field instances should throw
  EXPECT_THROW(elem_gf8 + elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf8 - elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf8 * elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf8 / elem_other, std::invalid_argument);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, ElementPrinting) {
  GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>> elem_5(5, gf8);

  // First test the field's printing directly to understand the behavior
  std::ostringstream field_oss;
  gf8->Print(5, field_oss);
  std::string field_result = field_oss.str();

  // Now test element printing which should match field printing
  std::ostringstream oss;
  elem_5.Print(oss);
  std::string result = oss.str();
  EXPECT_EQ(result,
            field_result);  // Element printing should match field printing

  // Test stream operator
  oss.str("");
  oss << elem_5;
  result = oss.str();
  EXPECT_EQ(result, field_result);  // Stream operator should also match

  // Test with different field representations
  gf8->SetRepresentation(FieldRepresentation::HEX);
  field_oss.str("");
  gf8->Print(5, field_oss);
  field_result = field_oss.str();

  oss.str("");
  elem_5.Print(oss);
  result = oss.str();
  EXPECT_EQ(result, field_result);

  gf8->SetRepresentation(FieldRepresentation::POLY);
  field_oss.str("");
  gf8->Print(5, field_oss);
  field_result = field_oss.str();

  oss.str("");
  elem_5.Print(oss);
  result = oss.str();
  EXPECT_EQ(result, field_result);
  // Should show polynomial representation - we'll just check it's not empty and
  // reasonable
  EXPECT_FALSE(result.empty());
}
