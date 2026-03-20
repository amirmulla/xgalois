#include "xgalois/field/gf_extension.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <memory>
#include <set>

#include "xgalois/field/gf_prime.hpp"

using namespace xg;

//===----------------------------------------------------------------------===//
// GF(p^n) Extension Field Tests
//===----------------------------------------------------------------------===//

class GaloisFieldExtensionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create GF(3^2) = GF(9)
    gf9 = std::make_unique<GaloisFieldExtension<uint8_t>>(std::make_pair(3, 2));

    // Create GF(5^2) = GF(25)
    gf25 =
        std::make_unique<GaloisFieldExtension<uint8_t>>(std::make_pair(5, 2));
  }

  std::unique_ptr<GaloisFieldExtension<uint8_t>> gf9, gf25;
};

TEST_F(GaloisFieldExtensionTest, FieldProperties) {
  EXPECT_EQ(gf9->Characteristic(), 3);
  EXPECT_EQ(gf9->Order(), 9);

  EXPECT_EQ(gf25->Characteristic(), 5);
  EXPECT_EQ(gf25->Order(), 25);
}

TEST_F(GaloisFieldExtensionTest, IdentityElements) {
  auto zero_poly = gf9->AdditiveIdentity();
  auto one_poly = gf9->MultiplicativeIdentity();

  EXPECT_EQ(zero_poly.Degree(), -1);  // Zero polynomial has degree -1

  EXPECT_EQ(one_poly.Degree(), 0);
  EXPECT_EQ(one_poly[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, BasicArithmetic) {
  auto one = gf9->MultiplicativeIdentity();
  auto zero = gf9->AdditiveIdentity();

  // Test addition
  auto sum = gf9->Add(one, zero);
  EXPECT_EQ(sum.Degree(), 0);
  EXPECT_EQ(sum[0].Value(), 1);

  // Test subtraction
  auto diff = gf9->Sub(one, zero);
  EXPECT_EQ(diff.Degree(), 0);
  EXPECT_EQ(diff[0].Value(), 1);

  // Test multiplication
  auto prod = gf9->Mul(one, one);
  EXPECT_EQ(prod.Degree(), 0);
  EXPECT_EQ(prod[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, MultiplicativeStructure) {
  // Test that the one element has an inverse (itself)
  auto one = gf9->MultiplicativeIdentity();
  auto one_inv = gf9->Inv(one);
  auto product = gf9->Mul(one, one_inv);

  EXPECT_EQ(product.Degree(), 0);
  EXPECT_EQ(product[0].Value(), 1);

  // Test division by one
  auto div_result = gf9->Div(one, one);
  EXPECT_EQ(div_result.Degree(), 0);
  EXPECT_EQ(div_result[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, PowerOperations) {
  auto one = gf9->MultiplicativeIdentity();

  // Test power operations
  auto one_squared = gf9->Pow(one, 2);
  EXPECT_EQ(one_squared.Degree(), 0);
  EXPECT_EQ(one_squared[0].Value(), 1);

  auto one_to_zero = gf9->Pow(one, 0);
  EXPECT_EQ(one_to_zero.Degree(), 0);
  EXPECT_EQ(one_to_zero[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, Random) {
  std::set<std::string> values;
  for (int i = 0; i < 50; ++i) {
    auto rand_elem = gf9->Random();
    values.insert(gf9->ToString(rand_elem));
  }
  EXPECT_GE(values.size(), 3);  // Should generate multiple different values
}

TEST_F(GaloisFieldExtensionTest, InvalidOperations) {
  auto zero = gf9->AdditiveIdentity();
  EXPECT_THROW(gf9->Inv(zero), std::domain_error);
  EXPECT_THROW(gf9->Div(gf9->MultiplicativeIdentity(), zero),
               std::domain_error);
  EXPECT_THROW(gf9->Log(zero), std::logic_error);   // Not implemented
  EXPECT_THROW(gf9->Sqrt(zero), std::logic_error);  // Not implemented
}

TEST_F(GaloisFieldExtensionTest, Representation) {
  EXPECT_EQ(gf9->GetRepresentation(), FieldRepresentation::POLY);

  gf9->SetRepresentation(FieldRepresentation::POLY);
  EXPECT_EQ(gf9->GetRepresentation(), FieldRepresentation::POLY);
}

TEST_F(GaloisFieldExtensionTest, ToString) {
  auto one = gf9->MultiplicativeIdentity();

  std::string str = gf9->ToString(one);
  EXPECT_FALSE(str.empty());
}

TEST_F(GaloisFieldExtensionTest, Print) {
  std::ostringstream oss;
  gf9->Print(oss);
  std::string field_str = oss.str();
  EXPECT_FALSE(field_str.empty());

  oss.str("");
  auto one = gf9->MultiplicativeIdentity();
  gf9->Print(one, oss);
  std::string elem_str = oss.str();
  EXPECT_FALSE(elem_str.empty());
}

TEST_F(GaloisFieldExtensionTest, Constructor) {
  // Test construction with order
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(9));   // 3^2
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(25));  // 5^2
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(8));   // 2^3

  // Test construction with pair
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(std::make_pair(3, 2)));
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(std::make_pair(5, 2)));

  // Test invalid orders
  EXPECT_THROW(GaloisFieldExtension<uint8_t>(6),
               std::invalid_argument);  // Not a prime power
  EXPECT_THROW(GaloisFieldExtension<uint8_t>(10),
               std::invalid_argument);  // Not a prime power
}

TEST_F(GaloisFieldExtensionTest, MultipleFields) {
  // Test operations in different extension fields
  auto one_gf9 = gf9->MultiplicativeIdentity();
  auto one_gf25 = gf25->MultiplicativeIdentity();

  // Each field should maintain its own identity
  EXPECT_EQ(one_gf9.Degree(), 0);
  EXPECT_EQ(one_gf9[0].Value(), 1);

  EXPECT_EQ(one_gf25.Degree(), 0);
  EXPECT_EQ(one_gf25[0].Value(), 1);

  // Powers should behave identically for the identity element
  auto one_squared_gf9 = gf9->Pow(one_gf9, 2);
  auto one_squared_gf25 = gf25->Pow(one_gf25, 2);

  EXPECT_EQ(one_squared_gf9.Degree(), 0);
  EXPECT_EQ(one_squared_gf9[0].Value(), 1);

  EXPECT_EQ(one_squared_gf25.Degree(), 0);
  EXPECT_EQ(one_squared_gf25[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, Performance) {
  // Simple performance test - should complete in reasonable time
  auto start = std::chrono::high_resolution_clock::now();

  auto one = gf9->MultiplicativeIdentity();
  auto result = one;

  // Perform 1000 multiplications
  for (int i = 0; i < 1000; ++i) {
    result = gf9->Mul(result, one);
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  EXPECT_LT(duration.count(), 1000);  // Should complete in less than 1 second

  // Result should still be one
  EXPECT_EQ(result.Degree(), 0);
  EXPECT_EQ(result[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, SetElementValueString) {
  // Test polynomial string parsing
  auto base_field_ptr = std::make_shared<GaloisFieldPrime<uint8_t>>(3);

  // Test simple polynomial: "α + 1" should be polynomial with coefficients [1,
  // 1]
  auto elem1 = gf9->SetElementValue("α + 1");
  EXPECT_EQ(elem1.Degree(), 1);
  EXPECT_EQ(elem1[0].Value(), 1);  // Constant term
  EXPECT_EQ(elem1[1].Value(), 1);  // α term

  // Test "α" should be polynomial with coefficient [0, 1]
  auto elem2 = gf9->SetElementValue("α");
  EXPECT_EQ(elem2.Degree(), 1);
  EXPECT_EQ(elem2[0].Value(), 0);  // Constant term
  EXPECT_EQ(elem2[1].Value(), 1);  // α term

  // Test "1" should be polynomial with coefficient [1]
  auto elem3 = gf9->SetElementValue("1");
  EXPECT_EQ(elem3.Degree(), 0);
  EXPECT_EQ(elem3[0].Value(), 1);

  // Test "0" should be zero polynomial
  auto elem4 = gf9->SetElementValue("0");
  EXPECT_EQ(elem4.Degree(), -1);  // Zero polynomial has degree -1

  // Test "α^2" - this will be reduced modulo the irreducible polynomial
  // In GF(3^2), α^2 gets reduced to some linear combination
  auto elem5 = gf9->SetElementValue("α^2");
  EXPECT_GE(elem5.Degree(), -1);  // Valid polynomial
  EXPECT_LT(elem5.Degree(), 2);   // Should be reduced to degree < 2

  // Test complex polynomial: "2*α + 2" in GF(3^2)
  auto elem6 = gf9->SetElementValue("2*α + 2");
  EXPECT_EQ(elem6.Degree(), 1);
  EXPECT_EQ(elem6[0].Value(), 2);  // Constant term
  EXPECT_EQ(elem6[1].Value(), 2);  // α term
}

TEST_F(GaloisFieldExtensionTest, SetElementValueStringInvalid) {
  // Test invalid string formats
  EXPECT_THROW(gf9->SetElementValue("invalid"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("β + 1"),
               std::invalid_argument);  // Wrong variable
  EXPECT_THROW(gf9->SetElementValue("h^1"),
               std::invalid_argument);  // Wrong generator name for power
  EXPECT_THROW(gf9->SetElementValue("g^"),
               std::invalid_argument);  // Incomplete power format
}

TEST_F(GaloisFieldExtensionTest, SetElementValuePolynomial) {
  auto base_field_ptr = std::make_shared<GaloisFieldPrime<uint8_t>>(3);

  // Create a test polynomial directly
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
  coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
      2, base_field_ptr));  // 2
  coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
      1, base_field_ptr));  // α

  PolynomialDense<GaloisFieldPrime<uint8_t>> test_poly(coeffs, "α");

  // Test polynomial SetElementValue with reduction
  auto result = gf9->SetElementValue(test_poly);

  // Should be reduced modulo the irreducible polynomial
  EXPECT_GE(result.Degree(), -1);  // Valid polynomial
  EXPECT_LT(result.Degree(), 2);  // Should be reduced to degree < 2 for GF(3^2)
}

TEST_F(GaloisFieldExtensionTest, GetElementValue) {
  // GetElementValue should be identity for extension fields
  auto one = gf9->MultiplicativeIdentity();
  auto result = gf9->GetElementValue(one);

  EXPECT_EQ(result.Degree(), one.Degree());
  EXPECT_EQ(result[0].Value(), one[0].Value());
}

TEST_F(GaloisFieldExtensionTest, SetElementValueReduction) {
  auto base_field_ptr = std::make_shared<GaloisFieldPrime<uint8_t>>(3);

  // Create a polynomial that exceeds the field degree (degree >= 2 for GF(3^2))
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
  coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
      1, base_field_ptr));  // 1
  coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
      0, base_field_ptr));  // 0*α
  coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
      0, base_field_ptr));  // 0*α^2
  coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
      1, base_field_ptr));  // 1*α^3

  PolynomialDense<GaloisFieldPrime<uint8_t>> high_degree_poly(coeffs, "α");

  // This should be reduced modulo the irreducible polynomial
  auto result = gf9->SetElementValue(high_degree_poly);

  // Result should have degree < 2 (the extension degree)
  EXPECT_LT(result.Degree(), 2);
}

TEST_F(GaloisFieldExtensionTest, SetElementValuePowerStringInvalid) {
  // Test that power string format parsing works correctly for invalid formats
  EXPECT_THROW(gf9->SetElementValue("g^"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("g^abc"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("h^1"),
               std::invalid_argument);  // Wrong generator name
}
