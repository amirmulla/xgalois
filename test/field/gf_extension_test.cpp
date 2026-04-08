#include "xgalois/field/gf_extension.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <memory>
#include <set>

#include "xgalois/field/gf_prime.hpp"

using namespace xg;

class GaloisFieldExtensionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf9 = std::make_unique<GaloisFieldExtension<uint8_t>>(std::make_pair(3, 2));

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

  EXPECT_EQ(zero_poly.Degree(), -1);

  EXPECT_EQ(one_poly.Degree(), 0);
  EXPECT_EQ(one_poly[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, BasicArithmetic) {
  auto one = gf9->MultiplicativeIdentity();
  auto zero = gf9->AdditiveIdentity();

  auto sum = gf9->Add(one, zero);
  EXPECT_EQ(sum.Degree(), 0);
  EXPECT_EQ(sum[0].Value(), 1);

  auto diff = gf9->Sub(one, zero);
  EXPECT_EQ(diff.Degree(), 0);
  EXPECT_EQ(diff[0].Value(), 1);

  auto prod = gf9->Mul(one, one);
  EXPECT_EQ(prod.Degree(), 0);
  EXPECT_EQ(prod[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, MultiplicativeStructure) {
  auto one = gf9->MultiplicativeIdentity();
  auto one_inv = gf9->Inv(one);
  auto product = gf9->Mul(one, one_inv);

  EXPECT_EQ(product.Degree(), 0);
  EXPECT_EQ(product[0].Value(), 1);

  auto div_result = gf9->Div(one, one);
  EXPECT_EQ(div_result.Degree(), 0);
  EXPECT_EQ(div_result[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, PowerOperations) {
  auto one = gf9->MultiplicativeIdentity();

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
  EXPECT_GE(values.size(), 3);
}

TEST_F(GaloisFieldExtensionTest, InvalidOperations) {
  auto zero = gf9->AdditiveIdentity();
  EXPECT_THROW(gf9->Inv(zero), std::domain_error);
  EXPECT_THROW(gf9->Div(gf9->MultiplicativeIdentity(), zero),
               std::domain_error);
  EXPECT_THROW(gf9->Log(zero), std::logic_error);
  EXPECT_THROW(gf9->Sqrt(zero), std::logic_error);
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
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(9));
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(25));
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(8));

  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(std::make_pair(3, 2)));
  EXPECT_NO_THROW(GaloisFieldExtension<uint8_t>(std::make_pair(5, 2)));

  EXPECT_THROW(GaloisFieldExtension<uint8_t>(6), std::invalid_argument);
  EXPECT_THROW(GaloisFieldExtension<uint8_t>(10), std::invalid_argument);
}

TEST_F(GaloisFieldExtensionTest, MultipleFields) {
  auto one_gf9 = gf9->MultiplicativeIdentity();
  auto one_gf25 = gf25->MultiplicativeIdentity();

  EXPECT_EQ(one_gf9.Degree(), 0);
  EXPECT_EQ(one_gf9[0].Value(), 1);

  EXPECT_EQ(one_gf25.Degree(), 0);
  EXPECT_EQ(one_gf25[0].Value(), 1);

  auto one_squared_gf9 = gf9->Pow(one_gf9, 2);
  auto one_squared_gf25 = gf25->Pow(one_gf25, 2);

  EXPECT_EQ(one_squared_gf9.Degree(), 0);
  EXPECT_EQ(one_squared_gf9[0].Value(), 1);

  EXPECT_EQ(one_squared_gf25.Degree(), 0);
  EXPECT_EQ(one_squared_gf25[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, Performance) {
  auto start = std::chrono::high_resolution_clock::now();

  auto one = gf9->MultiplicativeIdentity();
  auto result = one;

  for (int i = 0; i < 1000; ++i) {
    result = gf9->Mul(result, one);
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  EXPECT_LT(duration.count(), 1000);

  EXPECT_EQ(result.Degree(), 0);
  EXPECT_EQ(result[0].Value(), 1);
}

TEST_F(GaloisFieldExtensionTest, SetElementValueString) {
  auto base_field_ptr = std::make_shared<GaloisFieldPrime<uint8_t>>(3);

  auto elem1 = gf9->SetElementValue("α + 1");
  EXPECT_EQ(elem1.Degree(), 1);
  EXPECT_EQ(elem1[0].Value(), 1);
  EXPECT_EQ(elem1[1].Value(), 1);

  auto elem2 = gf9->SetElementValue("α");
  EXPECT_EQ(elem2.Degree(), 1);
  EXPECT_EQ(elem2[0].Value(), 0);
  EXPECT_EQ(elem2[1].Value(), 1);

  auto elem3 = gf9->SetElementValue("1");
  EXPECT_EQ(elem3.Degree(), 0);
  EXPECT_EQ(elem3[0].Value(), 1);

  auto elem4 = gf9->SetElementValue("0");
  EXPECT_EQ(elem4.Degree(), -1);

  auto elem5 = gf9->SetElementValue("α^2");
  EXPECT_GE(elem5.Degree(), -1);
  EXPECT_LT(elem5.Degree(), 2);

  auto elem6 = gf9->SetElementValue("2*α + 2");
  EXPECT_EQ(elem6.Degree(), 1);
  EXPECT_EQ(elem6[0].Value(), 2);
  EXPECT_EQ(elem6[1].Value(), 2);
}

TEST_F(GaloisFieldExtensionTest, SetElementValueStringInvalid) {
  EXPECT_THROW(gf9->SetElementValue("invalid"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("β + 1"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("h^1"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("g^"), std::invalid_argument);
}

TEST_F(GaloisFieldExtensionTest, SetElementValuePolynomial) {
  auto base_field_ptr = std::make_shared<GaloisFieldPrime<uint8_t>>(3);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
  coeffs.push_back(
      GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(2, base_field_ptr));
  coeffs.push_back(
      GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, base_field_ptr));

  PolynomialDense<GaloisFieldPrime<uint8_t>> test_poly(coeffs, "α");

  auto result = gf9->SetElementValue(test_poly);

  EXPECT_GE(result.Degree(), -1);
  EXPECT_LT(result.Degree(), 2);
}

TEST_F(GaloisFieldExtensionTest, GetElementValue) {
  auto one = gf9->MultiplicativeIdentity();
  auto result = gf9->GetElementValue(one);

  EXPECT_EQ(result.Degree(), one.Degree());
  EXPECT_EQ(result[0].Value(), one[0].Value());
}

TEST_F(GaloisFieldExtensionTest, SetElementValueReduction) {
  auto base_field_ptr = std::make_shared<GaloisFieldPrime<uint8_t>>(3);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
  coeffs.push_back(
      GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, base_field_ptr));
  coeffs.push_back(
      GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(0, base_field_ptr));
  coeffs.push_back(
      GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(0, base_field_ptr));
  coeffs.push_back(
      GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, base_field_ptr));

  PolynomialDense<GaloisFieldPrime<uint8_t>> high_degree_poly(coeffs, "α");

  auto result = gf9->SetElementValue(high_degree_poly);

  EXPECT_LT(result.Degree(), 2);
}

TEST_F(GaloisFieldExtensionTest, SetElementValuePowerStringInvalid) {
  EXPECT_THROW(gf9->SetElementValue("g^"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("g^abc"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("h^1"), std::invalid_argument);
}
