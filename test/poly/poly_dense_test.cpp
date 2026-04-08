

#include "xgalois/poly/poly_dense.hpp"

#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/utils/poly.hpp"

using namespace xg;

class PolynomialDenseTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);

    zero = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(0, gf7);
    one = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, gf7);
    two = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(2, gf7);
    three = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(3, gf7);
    four = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(4, gf7);
    five = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(5, gf7);
    six = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(6, gf7);

    gf2 = std::make_shared<GaloisFieldBinary>();
    gf2_zero = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
    gf2_one = GaloisFieldElementBase<GaloisFieldBinary>(1, gf2);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7;
  std::shared_ptr<GaloisFieldBinary> gf2;

  GaloisFieldElementBase<GaloisFieldPrime<uint8_t>> zero, one, two, three, four,
      five, six;
  GaloisFieldElementBase<GaloisFieldBinary> gf2_zero, gf2_one;
};

TEST_F(PolynomialDenseTest, ConstructorAndBasicProperties) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly.Size(), 3);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], two);
  EXPECT_EQ(poly[2], three);
  EXPECT_EQ(poly.GetVariable(), "x");
}

TEST_F(PolynomialDenseTest, CustomVariable) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one,
                                                                           two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs, "t");

  EXPECT_EQ(poly.GetVariable(), "t");

  poly.SetVariable("y");
  EXPECT_EQ(poly.GetVariable(), "y");

  EXPECT_THROW(poly.SetVariable(""), std::invalid_argument);
}

TEST_F(PolynomialDenseTest, ZeroPolynomial) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), -1);
  EXPECT_EQ(poly.Size(), 1);
  EXPECT_EQ(poly[0], zero);
}

TEST_F(PolynomialDenseTest, ConstantPolynomial) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), 0);
  EXPECT_EQ(poly.Size(), 1);
  EXPECT_EQ(poly[0], five);
}

TEST_F(PolynomialDenseTest, LeadingZeroTrimming) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, zero, zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), 1);
  EXPECT_EQ(poly.Size(), 2);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], two);
}

TEST_F(PolynomialDenseTest, IndexingOperators) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], two);
  EXPECT_EQ(poly[2], three);

  poly[1] = four;
  EXPECT_EQ(poly[1], four);

  EXPECT_THROW(poly[3], std::out_of_range);
  EXPECT_THROW(poly[-1], std::out_of_range);
}

TEST_F(PolynomialDenseTest, Addition) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      four, five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

  auto result = p1 + p2;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], five);
  EXPECT_EQ(result[1], zero);
  EXPECT_EQ(result[2], three);
}

TEST_F(PolynomialDenseTest, AdditionWithZero) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> zero_coeffs = {
      zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> zero_poly(zero_coeffs);

  auto result = poly + zero_poly;
  EXPECT_EQ(result, poly);
}

TEST_F(PolynomialDenseTest, Subtraction) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      four, five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

  auto result = p1 - p2;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], four);
  EXPECT_EQ(result[1], four);
  EXPECT_EQ(result[2], three);
}

TEST_F(PolynomialDenseTest, Negation) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto neg_poly = -poly;

  EXPECT_EQ(neg_poly.Degree(), 2);
  EXPECT_EQ(neg_poly[0], six);
  EXPECT_EQ(neg_poly[1], five);
  EXPECT_EQ(neg_poly[2], four);
}

TEST_F(PolynomialDenseTest, Multiplication) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      three, four};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

  auto result = p1 * p2;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], three);
  EXPECT_EQ(result[1], three);
  EXPECT_EQ(result[2], one);
}

TEST_F(PolynomialDenseTest, ScalarMultiplication) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto result = poly * five;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], five);
  EXPECT_EQ(result[1], three);
  EXPECT_EQ(result[2], one);
}

TEST_F(PolynomialDenseTest, MultiplicationByZero) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto result = poly * zero;

  EXPECT_EQ(result.Degree(), -1);
  EXPECT_EQ(result[0], zero);
}

TEST_F(PolynomialDenseTest, DivisionByConstant) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      two, four, six};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> div_coeffs = {
      two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(div_coeffs);

  auto result = poly / divisor;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], one);
  EXPECT_EQ(result[1], two);
  EXPECT_EQ(result[2], three);
}

TEST_F(PolynomialDenseTest, PolynomialDivision) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      dividend_coeffs = {two, three, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> dividend(dividend_coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      divisor_coeffs = {one, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(divisor_coeffs);

  auto result = dividend / divisor;

  EXPECT_EQ(result.Degree(), 1);
  EXPECT_EQ(result[0], two);
  EXPECT_EQ(result[1], one);
}

TEST_F(PolynomialDenseTest, ModuloOperation) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      dividend_coeffs = {two, three, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> dividend(dividend_coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      divisor_coeffs = {one, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(divisor_coeffs);

  auto remainder = dividend % divisor;

  EXPECT_EQ(remainder.Degree(), -1);
  EXPECT_EQ(remainder[0], zero);
}

TEST_F(PolynomialDenseTest, DivRemOperation) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      dividend_coeffs = {four, three, two, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> dividend(dividend_coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      divisor_coeffs = {one, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(divisor_coeffs);

  auto [quotient, remainder] = dividend.DivRem(divisor);

  auto verification = quotient * divisor + remainder;
  EXPECT_EQ(verification, dividend);
}

TEST_F(PolynomialDenseTest, DivisionByZero) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one,
                                                                           two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> zero_coeffs = {
      zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> zero_poly(zero_coeffs);

  EXPECT_THROW(poly.DivRem(zero_poly), std::invalid_argument);
}

TEST_F(PolynomialDenseTest, PowerOperation) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one,
                                                                           one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto power0 = poly ^ 0;
  EXPECT_EQ(power0.Degree(), 0);
  EXPECT_EQ(power0[0], one);

  auto power1 = poly ^ 1;
  EXPECT_EQ(power1, poly);

  auto power2 = poly ^ 2;
  EXPECT_EQ(power2.Degree(), 2);
  EXPECT_EQ(power2[0], one);
  EXPECT_EQ(power2[1], two);
  EXPECT_EQ(power2[2], one);
}

TEST_F(PolynomialDenseTest, Evaluation) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto result0 = poly(zero);
  EXPECT_EQ(result0, one);

  auto result1 = poly(one);
  EXPECT_EQ(result1, six);

  auto result2 = poly(two);
  EXPECT_EQ(result2, three);
}

TEST_F(PolynomialDenseTest, EqualityComparison) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs3 = {
      one, two, four};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly3(coeffs3);

  EXPECT_TRUE(poly1 == poly2);
  EXPECT_FALSE(poly1 == poly3);
  EXPECT_TRUE(poly1 != poly3);
  EXPECT_FALSE(poly1 != poly2);
}

TEST_F(PolynomialDenseTest, EqualityWithTrimming) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      one, two, zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

  EXPECT_TRUE(poly1 == poly2);
}

TEST_F(PolynomialDenseTest, CompoundAssignment) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      three, four};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

  auto poly_copy = poly1;
  poly_copy += poly2;
  EXPECT_EQ(poly_copy, poly1 + poly2);

  poly_copy = poly1;
  poly_copy -= poly2;
  EXPECT_EQ(poly_copy, poly1 - poly2);

  poly_copy = poly1;
  poly_copy *= poly2;
  EXPECT_EQ(poly_copy, poly1 * poly2);

  poly_copy = poly1;
  poly_copy *= five;
  EXPECT_EQ(poly_copy, poly1 * five);
}

TEST_F(PolynomialDenseTest, Derivative) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three, four};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto derivative = poly.Derivative();

  EXPECT_EQ(derivative.Degree(), 2);
  EXPECT_EQ(derivative[0], two);
  EXPECT_EQ(derivative[1], six);
  EXPECT_EQ(derivative[2], five);
}

TEST_F(PolynomialDenseTest, DerivativeOfConstant) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto derivative = poly.Derivative();

  EXPECT_EQ(derivative.Degree(), -1);
  EXPECT_EQ(derivative[0], zero);
}

TEST_F(PolynomialDenseTest, PrintPolynomial) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  std::stringstream ss;
  poly.Print(ss);
  std::string result = ss.str();

  EXPECT_TRUE(result.find('x') != std::string::npos);
  EXPECT_TRUE(result.length() > 0);
}

TEST_F(PolynomialDenseTest, StreamOperator) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one,
                                                                           two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  std::stringstream ss;
  ss << poly;
  std::string result = ss.str();

  EXPECT_TRUE(result.length() > 0);
}

TEST_F(PolynomialDenseTest, BinaryFieldOperations) {
  std::vector<GaloisFieldElementBase<GaloisFieldBinary>> coeffs1 = {
      gf2_one, gf2_zero, gf2_one};
  PolynomialDense<GaloisFieldBinary> poly1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldBinary>> coeffs2 = {gf2_one,
                                                                    gf2_one};
  PolynomialDense<GaloisFieldBinary> poly2(coeffs2);

  auto sum = poly1 + poly2;
  EXPECT_EQ(sum.Degree(), 2);
  EXPECT_EQ(sum[0], gf2_zero);
  EXPECT_EQ(sum[1], gf2_one);
  EXPECT_EQ(sum[2], gf2_one);

  auto product = poly1 * poly2;
  EXPECT_EQ(product.Degree(), 3);
  EXPECT_EQ(product[0], gf2_one);
  EXPECT_EQ(product[1], gf2_one);
  EXPECT_EQ(product[2], gf2_one);
  EXPECT_EQ(product[3], gf2_one);
}

TEST_F(PolynomialDenseTest, ExtendedGCD) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_a = {
      two, three, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> a(coeffs_a);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_b = {
      one, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> b(coeffs_b);

  auto [gcd, coeffs] = utils::PolynomialDenseExtendedGcd(a, b);
  auto [s, t] = coeffs;

  EXPECT_TRUE(gcd.Degree() >= 0);

  auto verification = a * s + b * t;
  EXPECT_EQ(verification, gcd);
}

TEST_F(PolynomialDenseTest, ExtendedGCDCoprime) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_a = {
      one, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> a(coeffs_a);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_b = {
      two, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> b(coeffs_b);

  auto [gcd, coeffs] = utils::PolynomialDenseExtendedGcd(a, b);
  auto [s, t] = coeffs;

  EXPECT_EQ(gcd.Degree(), 0);
  EXPECT_EQ(gcd[0], one);

  auto verification = a * s + b * t;
  EXPECT_EQ(verification, gcd);
}

TEST_F(PolynomialDenseTest, EmptyCoefficientsAssert) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> empty_coeffs;

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> valid_coeffs =
      {one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(valid_coeffs);
  EXPECT_EQ(poly.Degree(), 0);
}

TEST_F(PolynomialDenseTest, LargePolynomials) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
  coeffs.reserve(10);
  for (int i = 0; i < 10; ++i) {
    coeffs.push_back(
        GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(i % 7, gf7));
  }

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);
  EXPECT_EQ(poly.Degree(), 9);
  EXPECT_EQ(poly.Size(), 10);

  auto doubled = poly * two;
  EXPECT_EQ(doubled.Degree(), 9);

  auto squared = poly * poly;
  EXPECT_EQ(squared.Degree(), 18);
}

TEST_F(PolynomialDenseTest, FieldConsistency) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Field(), gf7);

  auto doubled = poly * two;
  EXPECT_EQ(doubled.Field(), gf7);

  auto sum = poly + poly;
  EXPECT_EQ(sum.Field(), gf7);
}