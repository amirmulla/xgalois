#include "xgalois/field/gf_element.hpp"

#include <gtest/gtest.h>

#include <memory>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

class GaloisFieldElementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    gf2 = std::make_shared<GaloisFieldBinary>();

    elem_3 = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(3, gf7);
    elem_5 = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(5, gf7);
    elem_0 = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(0, gf7);
    elem_1 = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, gf7);

    bin_0 = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
    bin_1 = GaloisFieldElementBase<GaloisFieldBinary>(1, gf2);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7;
  std::shared_ptr<GaloisFieldBinary> gf2;

  GaloisFieldElementBase<GaloisFieldPrime<uint8_t>> elem_3, elem_5, elem_0,
      elem_1;
  GaloisFieldElementBase<GaloisFieldBinary> bin_0, bin_1;
};

TEST_F(GaloisFieldElementTest, Construction) {
  EXPECT_EQ(elem_3.Value(), 3);
  EXPECT_EQ(elem_3.Field(), gf7);

  EXPECT_THROW(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(3, nullptr),
               std::invalid_argument);
}

TEST_F(GaloisFieldElementTest, Addition) {
  auto result = elem_3 + elem_5;
  EXPECT_EQ(result.Value(), 1);
  EXPECT_EQ(result.Field(), gf7);

  auto zero_sum =
      elem_3 + GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(4, gf7);
  EXPECT_EQ(zero_sum.Value(), 0);
}

TEST_F(GaloisFieldElementTest, Subtraction) {
  auto result = elem_5 - elem_3;
  EXPECT_EQ(result.Value(), 2);

  auto neg_result = elem_3 - elem_5;
  EXPECT_EQ(neg_result.Value(), 5);
}

TEST_F(GaloisFieldElementTest, Multiplication) {
  auto result = elem_3 * elem_5;
  EXPECT_EQ(result.Value(), 1);

  auto zero_result = elem_0 * elem_5;
  EXPECT_EQ(zero_result.Value(), 0);

  auto identity_result = elem_1 * elem_5;
  EXPECT_EQ(identity_result.Value(), 5);
}

TEST_F(GaloisFieldElementTest, Division) {
  auto result = elem_1 / elem_3;
  EXPECT_EQ(result.Value(), 5);

  auto self_div = elem_5 / elem_5;
  EXPECT_EQ(self_div.Value(), 1);

  EXPECT_THROW(elem_3 / elem_0, std::domain_error);
}

TEST_F(GaloisFieldElementTest, Negation) {
  auto neg_3 = -elem_3;
  EXPECT_EQ(neg_3.Value(), 4);

  auto neg_0 = -elem_0;
  EXPECT_EQ(neg_0.Value(), 0);
}

TEST_F(GaloisFieldElementTest, Power) {
  auto power_result = elem_3 ^ 2;
  EXPECT_EQ(power_result.Value(), 2);

  auto pow_method = elem_3.Pow(3);
  EXPECT_EQ(pow_method.Value(), 6);

  auto pow_zero = elem_5 ^ 0;
  EXPECT_EQ(pow_zero.Value(), 1);
}

TEST_F(GaloisFieldElementTest, Inverse) {
  auto inv_3 = elem_3.Inv();
  EXPECT_EQ(inv_3.Value(), 5);

  auto check = elem_3 * inv_3;
  EXPECT_EQ(check.Value(), 1);

  EXPECT_THROW(elem_0.Inv(), std::domain_error);
}

TEST_F(GaloisFieldElementTest, SquareRoot) {
  auto sqrt_1 = bin_1.Sqrt();
  EXPECT_EQ(sqrt_1.Value(), 1);

  auto sqrt_0 = bin_0.Sqrt();
  EXPECT_EQ(sqrt_0.Value(), 0);

  EXPECT_THROW(elem_3.Sqrt(), std::logic_error);
}

TEST_F(GaloisFieldElementTest, CompoundAssignment) {
  auto elem_copy = elem_3;
  elem_copy += elem_5;
  EXPECT_EQ(elem_copy.Value(), 1);

  elem_copy = elem_3;
  elem_copy -= elem_5;
  EXPECT_EQ(elem_copy.Value(), 5);

  elem_copy = elem_3;
  elem_copy *= elem_5;
  EXPECT_EQ(elem_copy.Value(), 1);

  elem_copy = elem_1;
  elem_copy /= elem_3;
  EXPECT_EQ(elem_copy.Value(), 5);
}

TEST_F(GaloisFieldElementTest, Comparison) {
  auto elem_3_copy = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(3, gf7);
  auto elem_4 = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(4, gf7);

  EXPECT_TRUE(elem_3 == elem_3_copy);
  EXPECT_FALSE(elem_3 == elem_4);
  EXPECT_TRUE(elem_3 != elem_4);
  EXPECT_FALSE(elem_3 != elem_3_copy);
}

TEST_F(GaloisFieldElementTest, BinaryFieldSpecific) {
  auto sum = bin_0 + bin_1;
  EXPECT_EQ(sum.Value(), 1);

  auto sum2 = bin_1 + bin_1;
  EXPECT_EQ(sum2.Value(), 0);

  auto prod = bin_1 * bin_1;
  EXPECT_EQ(prod.Value(), 1);

  auto sqrt_bin = bin_1.Sqrt();
  EXPECT_EQ(sqrt_bin.Value(), 1);
}

TEST_F(GaloisFieldElementTest, StringOperations) {
  auto str_3 = gf7->ToString(elem_3.Value());
  EXPECT_EQ(str_3, "3");

  auto str_bin = gf2->ToString(bin_1.Value());
  EXPECT_EQ(str_bin, "1");
}

TEST_F(GaloisFieldElementTest, StreamOutput) {
  std::ostringstream oss;
  oss << elem_3;
  EXPECT_EQ(oss.str(), "3");

  oss.str("");
  oss << bin_1;
  EXPECT_EQ(oss.str(), "1");
}

TEST_F(GaloisFieldElementTest, LogarithmOperations) {
  uint8_t gen = gf7->MultiplicativeGenerator();
  auto gen_elem = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(gen, gf7);

  uint32_t log_val = gf7->Log(elem_3.Value());
  auto check = gen_elem.Pow(log_val);
  EXPECT_EQ(check.Value(), 3);

  uint32_t log_val_base = gf7->Log(elem_3.Value(), gen_elem.Value());
  auto check_base = gen_elem.Pow(log_val_base);
  EXPECT_EQ(check_base.Value(), 3);

  EXPECT_THROW(gf7->Log(elem_0.Value()), std::domain_error);
}

class GaloisFieldElementStrictTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);

    elem_3_gf7 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(3, gf7);
    elem_2_gf5 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(2, gf5);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7, gf5;
  GaloisFieldElement<GaloisFieldPrime<uint8_t>> elem_3_gf7, elem_2_gf5;
};

TEST_F(GaloisFieldElementStrictTest, SameFieldOperations) {
  auto elem_4_gf7 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(4, gf7);

  auto sum = elem_3_gf7 + elem_4_gf7;
  EXPECT_EQ(sum.Value(), 0);

  auto prod = elem_3_gf7 * elem_4_gf7;
  EXPECT_EQ(prod.Value(), 5);
}

TEST_F(GaloisFieldElementStrictTest, DifferentFieldOperations) {
  EXPECT_THROW(elem_3_gf7 + elem_2_gf5, std::invalid_argument);
  EXPECT_THROW(elem_3_gf7 * elem_2_gf5, std::invalid_argument);
  EXPECT_THROW(elem_3_gf7 - elem_2_gf5, std::invalid_argument);
  EXPECT_THROW(elem_3_gf7 / elem_2_gf5, std::invalid_argument);
}
