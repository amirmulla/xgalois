

#include "xgalois/utils/poly.hpp"

#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"

using namespace xg;
using namespace xg::utils;

class PolyUtilsTest : public ::testing::Test {
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

TEST_F(PolyUtilsTest, ExtractCoefficientConstantTerm) {
  using GF7Field = GaloisFieldPrime<uint8_t>;

  EXPECT_EQ(ExtractCoefficient<GF7Field>("5"), 5);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("0"), 0);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("123"), 123);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("1"), 1);
}

TEST_F(PolyUtilsTest, ExtractCoefficientVariableTerms) {
  using GF7Field = GaloisFieldPrime<uint8_t>;

  EXPECT_EQ(ExtractCoefficient<GF7Field>("3x"), 3);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("5x^2"), 5);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("12x^3"), 12);

  EXPECT_EQ(ExtractCoefficient<GF7Field>("x"), 1);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("x^2"), 1);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("x^10"), 1);

  EXPECT_EQ(ExtractCoefficient<GF7Field>("+x"), 1);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("+5x"), 5);

  EXPECT_EQ(ExtractCoefficient<GF7Field>("-x"), 1);
}

TEST_F(PolyUtilsTest, ExtractCoefficientUppercaseX) {
  using GF7Field = GaloisFieldPrime<uint8_t>;

  EXPECT_EQ(ExtractCoefficient<GF7Field>("3X"), 3);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("X"), 1);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("5X^2"), 5);
}

TEST_F(PolyUtilsTest, ExtractCoefficientInvalidTerms) {
  using GF7Field = GaloisFieldPrime<uint8_t>;

  EXPECT_THROW(ExtractCoefficient<GF7Field>("abc"), std::invalid_argument);

  EXPECT_EQ(ExtractCoefficient<GF7Field>("3.5x"), 3);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("ax"), std::invalid_argument);
}

TEST_F(PolyUtilsTest, ExtractDegreeConstantTerm) {
  EXPECT_EQ(ExtractDegree("5"), 0);
  EXPECT_EQ(ExtractDegree("0"), 0);
  EXPECT_EQ(ExtractDegree("123"), 0);
  EXPECT_EQ(ExtractDegree("1"), 0);
}

TEST_F(PolyUtilsTest, ExtractDegreeLinearTerms) {
  EXPECT_EQ(ExtractDegree("x"), 1);
  EXPECT_EQ(ExtractDegree("X"), 1);
  EXPECT_EQ(ExtractDegree("3x"), 1);
  EXPECT_EQ(ExtractDegree("123x"), 1);
  EXPECT_EQ(ExtractDegree("-x"), 1);
  EXPECT_EQ(ExtractDegree("+x"), 1);
}

TEST_F(PolyUtilsTest, ExtractDegreeHigherDegreeTerms) {
  EXPECT_EQ(ExtractDegree("x^2"), 2);
  EXPECT_EQ(ExtractDegree("X^2"), 2);
  EXPECT_EQ(ExtractDegree("3x^2"), 2);
  EXPECT_EQ(ExtractDegree("x^10"), 10);
  EXPECT_EQ(ExtractDegree("5x^100"), 100);
  EXPECT_EQ(ExtractDegree("-x^3"), 3);
  EXPECT_EQ(ExtractDegree("+2x^5"), 5);
}

TEST_F(PolyUtilsTest, ExtractDegreeZeroDegree) {
  EXPECT_EQ(ExtractDegree("x^0"), 0);
  EXPECT_EQ(ExtractDegree("5x^0"), 0);
}

TEST_F(PolyUtilsTest, ExtractDegreeInvalidTerms) {
  EXPECT_THROW(ExtractDegree("x^"), std::invalid_argument);
  EXPECT_THROW(ExtractDegree("x^abc"), std::invalid_argument);

  EXPECT_EQ(ExtractDegree("x^3.5"), 3);
}

TEST_F(PolyUtilsTest, ParsePolynomialConstant) {
  auto poly = ParsePolynomial(gf7, "5");
  EXPECT_EQ(poly.Degree(), 0);
  EXPECT_EQ(poly[0], five);

  auto zero_poly = ParsePolynomial(gf7, "0");
  EXPECT_EQ(zero_poly.Degree(), -1);
  EXPECT_EQ(zero_poly[0], zero);
}

TEST_F(PolyUtilsTest, ParsePolynomialLinear) {
  auto poly = ParsePolynomial(gf7, "3x + 2");
  EXPECT_EQ(poly.Degree(), 1);
  EXPECT_EQ(poly[0], two);
  EXPECT_EQ(poly[1], three);

  auto poly2 = ParsePolynomial(gf7, "x");
  EXPECT_EQ(poly2.Degree(), 1);
  EXPECT_EQ(poly2[0], zero);
  EXPECT_EQ(poly2[1], one);
}

TEST_F(PolyUtilsTest, ParsePolynomialQuadratic) {
  auto poly = ParsePolynomial(gf7, "2x^2 + 3x + 1");
  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], three);
  EXPECT_EQ(poly[2], two);
}

TEST_F(PolyUtilsTest, ParsePolynomialWithSubtraction) {
  auto poly = ParsePolynomial(gf7, "x^2 - 3x + 2");
  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly[0], two);
  EXPECT_EQ(poly[1], four);
  EXPECT_EQ(poly[2], one);
}

TEST_F(PolyUtilsTest, ParsePolynomialWithWhitespace) {
  auto poly1 = ParsePolynomial(gf7, "x^2 + 3x + 1");
  auto poly2 = ParsePolynomial(gf7, "x^2+3x+1");
  auto poly3 = ParsePolynomial(gf7, "  x^2  +  3x  +  1  ");

  EXPECT_EQ(poly1, poly2);
  EXPECT_EQ(poly2, poly3);
}

TEST_F(PolyUtilsTest, ParsePolynomialSparseTerms) {
  auto poly = ParsePolynomial(gf7, "x^5 + 3x^2 + 1");
  EXPECT_EQ(poly.Degree(), 5);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], zero);
  EXPECT_EQ(poly[2], three);
  EXPECT_EQ(poly[3], zero);
  EXPECT_EQ(poly[4], zero);
  EXPECT_EQ(poly[5], one);
}

TEST_F(PolyUtilsTest, ParsePolynomialBinaryField) {
  auto poly = ParsePolynomial(gf2, "x^2 + x + 1");
  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly[0], gf2_one);
  EXPECT_EQ(poly[1], gf2_one);
  EXPECT_EQ(poly[2], gf2_one);
}

TEST_F(PolyUtilsTest, ParsePolynomialInvalidInput) {
  EXPECT_THROW(ParsePolynomial(gf7, ""), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "   "), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "abc"), std::invalid_argument);

  auto poly = ParsePolynomial(gf7, "3.5x + 2");
  EXPECT_EQ(poly.Degree(), 1);
  EXPECT_EQ(poly[1], three);
}

TEST_F(PolyUtilsTest, PolynomialExtendedGcdBasicCase) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      six, zero, one};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      six, one};

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s, t] = bezout;

  auto left_side = poly1 * s + poly2 * t;
  EXPECT_EQ(left_side, gcd);

  if (gcd.Degree() >= 0) {
    EXPECT_EQ(gcd[gcd.Degree()], one);
  }
}

TEST_F(PolyUtilsTest, PolynomialExtendedGcdWithZero) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      two, three};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      zero};

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s, t] = bezout;

  EXPECT_EQ(gcd.Degree(), 1);
  EXPECT_EQ(gcd[1], one);

  auto left_side = poly1 * s + poly2 * t;
  EXPECT_EQ(left_side, gcd);
}

TEST_F(PolyUtilsTest, PolynomialExtendedGcdCoprimePolynomials) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, zero, one};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      one, one};

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s, t] = bezout;

  EXPECT_EQ(gcd.Degree(), 0);
  EXPECT_EQ(gcd[0], one);

  auto left_side = poly1 * s + poly2 * t;
  EXPECT_EQ(left_side, gcd);
}

TEST_F(PolyUtilsTest, IsIrreducibleConstantPolynomials) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_zero = {
      zero};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_one = {
      one};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_five = {
      five};

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly_zero(coeffs_zero);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly_one(coeffs_one);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly_five(coeffs_five);

  EXPECT_FALSE(IsIrreducible(poly_zero));
  EXPECT_FALSE(IsIrreducible(poly_one));
  EXPECT_FALSE(IsIrreducible(poly_five));
}

TEST_F(PolyUtilsTest, IsIrreducibleLinearPolynomials) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, one};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      two, three};
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs3 = {
      zero, one};

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly3(coeffs3);

  EXPECT_TRUE(IsIrreducible(poly1));
  EXPECT_TRUE(IsIrreducible(poly2));
  EXPECT_TRUE(IsIrreducible(poly3));
}

TEST_F(PolyUtilsTest, IsIrreducibleQuadraticPolynomials) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      coeffs_irreducible = {one, zero, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly_irreducible(
      coeffs_irreducible);
  EXPECT_TRUE(IsIrreducible(poly_irreducible));

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>
      coeffs_reducible = {six, zero, one};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly_reducible(coeffs_reducible);
  EXPECT_FALSE(IsIrreducible(poly_reducible));
}

TEST_F(PolyUtilsTest, IsIrreducibleBinaryField) {
  std::vector<GaloisFieldElementBase<GaloisFieldBinary>> coeffs_irreducible = {
      gf2_one, gf2_one, gf2_one};
  PolynomialDense<GaloisFieldBinary> poly_irreducible(coeffs_irreducible);
  EXPECT_TRUE(IsIrreducible(poly_irreducible));

  std::vector<GaloisFieldElementBase<GaloisFieldBinary>> coeffs_reducible = {
      gf2_one, gf2_zero, gf2_one};
  PolynomialDense<GaloisFieldBinary> poly_reducible(coeffs_reducible);
  EXPECT_FALSE(IsIrreducible(poly_reducible));
}

TEST_F(PolyUtilsTest, IsIrreducibleInvalidPolynomial) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_TRUE(IsIrreducible(poly));
}

TEST_F(PolyUtilsTest, ParseAndEvaluatePolynomial) {
  auto poly = ParsePolynomial(gf7, "2x^2 + 3x + 1");

  auto result = poly(two);
  EXPECT_EQ(result, one);

  auto result_zero = poly(zero);
  EXPECT_EQ(result_zero, one);
}

TEST_F(PolyUtilsTest, ParsePolynomialThenExtendedGcd) {
  auto poly1 = ParsePolynomial(gf7, "x^2 - 1");
  auto poly2 = ParsePolynomial(gf7, "x - 1");

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);

  EXPECT_EQ(gcd.Degree(), 1);
  EXPECT_EQ(gcd[1], one);
}

TEST_F(PolyUtilsTest, ParsePolynomialThenIrreducibilityTest) {
  auto irreducible_poly = ParsePolynomial(gf7, "x^2 + 1");
  EXPECT_TRUE(IsIrreducible(irreducible_poly));

  auto reducible_poly = ParsePolynomial(gf7, "x^2 - 1");
  EXPECT_FALSE(IsIrreducible(reducible_poly));
}
