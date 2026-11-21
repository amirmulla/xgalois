#include "xgalois/utils/poly.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"

using namespace xg;
using namespace xg::utils;

class PolyUtilsTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    zero = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(0, gf7);
    one = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(1, gf7);
    two = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(2, gf7);
    three = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(3, gf7);
    four = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(4, gf7);
    five = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(5, gf7);
    six = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(6, gf7);

    gf2 = std::make_shared<GaloisFieldBinary>();
    gf2_zero = GaloisFieldElement<GaloisFieldBinary>(0, gf2);
    gf2_one = GaloisFieldElement<GaloisFieldBinary>(1, gf2);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7;
  std::shared_ptr<GaloisFieldBinary> gf2;

  GaloisFieldElement<GaloisFieldPrime<uint8_t>> zero, one, two, three, four,
      five, six;
  GaloisFieldElement<GaloisFieldBinary> gf2_zero, gf2_one;
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
  EXPECT_EQ(ExtractCoefficient<GF7Field>("3*x"), 3);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("5x^2"), 5);
  EXPECT_EQ(ExtractCoefficient<GF7Field>("5*x^2"), 5);
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

  EXPECT_THROW(ExtractCoefficient<GF7Field>(""), std::invalid_argument);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("abc"), std::invalid_argument);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("3.5x"), std::invalid_argument);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("ax"), std::invalid_argument);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("2abcx"), std::invalid_argument);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("x2"), std::invalid_argument);
  EXPECT_THROW(ExtractCoefficient<GF7Field>("*x"), std::invalid_argument);
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
  EXPECT_EQ(ExtractDegree("3*x"), 1);
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
  EXPECT_THROW(ExtractDegree(""), std::invalid_argument);
  EXPECT_THROW(ExtractDegree("x^"), std::invalid_argument);
  EXPECT_THROW(ExtractDegree("x^abc"), std::invalid_argument);
  EXPECT_THROW(ExtractDegree("x^3.5"), std::invalid_argument);
  EXPECT_THROW(ExtractDegree("x^-2"), std::invalid_argument);
  EXPECT_THROW(ExtractDegree("x2"), std::invalid_argument);
}

TEST_F(PolyUtilsTest, ParsePolynomialConstantAndNegative) {
  auto poly = ParsePolynomial(gf7, "5");
  EXPECT_EQ(poly.Degree(), 0);
  EXPECT_EQ(poly[0], five);

  auto neg_const = ParsePolynomial(gf7, "-5");
  EXPECT_EQ(neg_const.Degree(), 0);
  EXPECT_EQ(neg_const[0], two);

  auto neg_large = ParsePolynomial(gf7, "-10");
  EXPECT_EQ(neg_large.Degree(), 0);
  EXPECT_EQ(neg_large[0], four);

  auto zero_poly = ParsePolynomial(gf7, "0");
  EXPECT_TRUE(zero_poly.IsZero());
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

  auto neg_x = ParsePolynomial(gf7, "-x");
  EXPECT_EQ(neg_x.Degree(), 1);
  EXPECT_EQ(neg_x[0], zero);
  EXPECT_EQ(neg_x[1], six);

  auto neg_10x = ParsePolynomial(gf7, "-10x");
  EXPECT_EQ(neg_10x.Degree(), 1);
  EXPECT_EQ(neg_10x[1], four);
}

TEST_F(PolyUtilsTest, ParsePolynomialQuadraticAndMultiplicationStar) {
  auto poly = ParsePolynomial(gf7, "2x^2 + 3x + 1");
  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], three);
  EXPECT_EQ(poly[2], two);

  auto star_poly = ParsePolynomial(gf7, "2*x^2 + 3*x + 1");
  EXPECT_EQ(star_poly, poly);

  auto leading_neg = ParsePolynomial(gf7, "-x^2 + 1");
  EXPECT_EQ(leading_neg.Degree(), 2);
  EXPECT_EQ(leading_neg[0], one);
  EXPECT_EQ(leading_neg[1], zero);
  EXPECT_EQ(leading_neg[2], six);
}

TEST_F(PolyUtilsTest, ParsePolynomialWithSubtractionAndModularReduction) {
  auto poly = ParsePolynomial(gf7, "x^2 - 3x + 2");
  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly[0], two);
  EXPECT_EQ(poly[1], four);
  EXPECT_EQ(poly[2], one);

  // Coefficient > p before negation: -10 === -(3) === 4 mod 7
  auto poly2 = ParsePolynomial(gf7, "x^2 - 10x + 1");
  EXPECT_EQ(poly2.Degree(), 2);
  EXPECT_EQ(poly2[0], one);
  EXPECT_EQ(poly2[1], four);
  EXPECT_EQ(poly2[2], one);

  // Duplicate degree terms accumulate mod p.
  auto dup = ParsePolynomial(gf7, "x^2 + 3*x^2 - 2x^2 + 5x + 2x");
  EXPECT_EQ(dup.Degree(), 2);
  EXPECT_EQ(dup[0], zero);
  EXPECT_EQ(dup[1], zero);
  EXPECT_EQ(dup[2], two);
}

TEST_F(PolyUtilsTest, ParsePolynomialCustomVariable) {
  auto poly = ParsePolynomial(gf7, "2*t^3 - t + 4", "t");
  EXPECT_EQ(poly.GetVariable(), "t");
  EXPECT_EQ(poly.Degree(), 3);
  EXPECT_EQ(poly[0], four);
  EXPECT_EQ(poly[1], six);
  EXPECT_EQ(poly[2], zero);
  EXPECT_EQ(poly[3], two);

  // Must NOT silently merge 'x' and 't'.
  EXPECT_THROW(ParsePolynomial(gf7, "x + t", "t"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "t^2 + 1", ""), std::invalid_argument);
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

TEST_F(PolyUtilsTest, ParsePolynomialBinaryAndExtensionFields) {
  auto poly = ParsePolynomial(gf2, "x^4 + x + 1");
  EXPECT_EQ(poly.Degree(), 4);
  EXPECT_EQ(poly[0], gf2_one);
  EXPECT_EQ(poly[1], gf2_one);
  EXPECT_EQ(poly[2], gf2_zero);
  EXPECT_EQ(poly[3], gf2_zero);
  EXPECT_EQ(poly[4], gf2_one);

  // GF(2^4) with hex and integer coefficients
  auto gf2_4 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4);
  auto p_bin_ext = ParsePolynomial(gf2_4, "(0x5)x^2 + 3*x + 0x7");
  EXPECT_EQ(p_bin_ext.Degree(), 2);
  EXPECT_EQ(p_bin_ext[0].Value(), 7u);
  EXPECT_EQ(p_bin_ext[1].Value(), 3u);
  EXPECT_EQ(p_bin_ext[2].Value(), 5u);

  // GF(3^2) with parenthesized extension coefficients
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  auto p_ext = ParsePolynomial(gf3_2, "(α + 1)x^2 - (2α)x + 2");
  EXPECT_EQ(p_ext.Degree(), 2);
  EXPECT_EQ(p_ext[0],
            GaloisFieldElement<GaloisFieldExtension<uint32_t>>("2", gf3_2));
  EXPECT_EQ(p_ext[1],
            GaloisFieldElement<GaloisFieldExtension<uint32_t>>("α", gf3_2));
  EXPECT_EQ(p_ext[2],
            GaloisFieldElement<GaloisFieldExtension<uint32_t>>("α + 1", gf3_2));
}

TEST_F(PolyUtilsTest, ParsePolynomialInvalidInput) {
  EXPECT_THROW(ParsePolynomial(gf7, ""), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "   "), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "abc"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "3.5x + 2"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "x^3.5"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "2abcx"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "x2"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "x + "), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "++x"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "x + -1"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "*x"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "()x"), std::invalid_argument);
  EXPECT_THROW(ParsePolynomial(gf7, "(2x"), std::invalid_argument);
}

TEST_F(PolyUtilsTest, BinaryPolynomialToUintTest) {
  auto z = PolynomialDense<GaloisFieldBinary>::Zero(gf2);
  EXPECT_EQ(BinaryPolynomialToUint(z), 0ULL);

  auto one_p = PolynomialDense<GaloisFieldBinary>::One(gf2);
  EXPECT_EQ(BinaryPolynomialToUint(one_p), 1ULL);

  auto conway4 = ParsePolynomial(gf2, "x^4 + x + 1");
  EXPECT_EQ(BinaryPolynomialToUint(conway4), 0b10011ULL);

  // Degree 63 is valid (bit 63 set); degree 64 must throw.
  auto deg63 = PolynomialDense<GaloisFieldBinary>::Monomial(gf2, 63, 1);
  EXPECT_EQ(BinaryPolynomialToUint(deg63), (uint64_t{1} << 63));
  auto deg64 = PolynomialDense<GaloisFieldBinary>::Monomial(gf2, 64, 1);
  EXPECT_THROW(BinaryPolynomialToUint(deg64), std::invalid_argument);

  // Non-characteristic-2 field must throw.
  auto p7 = ParsePolynomial(gf7, "x^2 + 1");
  EXPECT_THROW(BinaryPolynomialToUint(p7), std::invalid_argument);

  // GF(2^4) with {0, 1} coefficients succeeds; with coeff 2 throws.
  auto gf2_4 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4);
  auto bin_ok = ParsePolynomial(gf2_4, "x^3 + x + 1");
  EXPECT_EQ(BinaryPolynomialToUint(bin_ok), 0b1011ULL);
  auto bin_bad = ParsePolynomial(gf2_4, "2x + 1");
  EXPECT_THROW(BinaryPolynomialToUint(bin_bad), std::invalid_argument);
}

TEST_F(PolyUtilsTest, PolynomialGcdLcmAndExtendedGcd) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(
      {six, zero, one});                                         // x^2 - 1
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2({six, one});  // x - 1

  auto g = PolynomialGcd(poly1, poly2);
  EXPECT_EQ(g, poly2);
  EXPECT_TRUE(g.IsMonic());

  auto l = PolynomialLcm(poly1, poly2);
  EXPECT_EQ(l, poly1);

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s, t] = bezout;
  EXPECT_EQ(gcd, g);
  EXPECT_EQ(poly1 * s + poly2 * t, gcd);
  EXPECT_TRUE(gcd.IsMonic());

  // Field mismatch must throw.
  auto gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);
  auto p5 = ParsePolynomial(gf5, "x + 1");
  EXPECT_THROW(PolynomialGcd(poly1, p5), std::invalid_argument);
  EXPECT_THROW(PolynomialLcm(poly1, p5), std::invalid_argument);
  EXPECT_THROW(PolynomialDenseExtendedGcd(poly1, p5), std::invalid_argument);
}

TEST_F(PolyUtilsTest, PolynomialExtendedGcdWithZero) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1({two, three});
  auto poly2 = PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(gf7);

  auto [gcd1, bezout1] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s1, t1] = bezout1;
  EXPECT_EQ(gcd1.Degree(), 1);
  EXPECT_TRUE(gcd1.IsMonic());
  EXPECT_EQ(poly1 * s1 + poly2 * t1, gcd1);

  auto [gcd2, bezout2] = PolynomialDenseExtendedGcd(poly2, poly1);
  auto [s2, t2] = bezout2;
  EXPECT_EQ(gcd2, gcd1);
  EXPECT_EQ(poly2 * s2 + poly1 * t2, gcd2);
}

TEST_F(PolyUtilsTest, PolynomialExtendedGcdCoprimePolynomials) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1({one, zero, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2({one, one});

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s, t] = bezout;

  EXPECT_EQ(gcd.Degree(), 0);
  EXPECT_EQ(gcd[0], one);
  EXPECT_EQ(poly1 * s + poly2 * t, gcd);
}

TEST_F(PolyUtilsTest, IsIrreducibleConstantAndLinearPolynomials) {
  auto poly_zero = PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(gf7);
  auto poly_one = PolynomialDense<GaloisFieldPrime<uint8_t>>::One(gf7);
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly_five({five});

  EXPECT_FALSE(IsIrreducible(poly_zero));
  EXPECT_FALSE(IsIrreducible(poly_one));
  EXPECT_FALSE(IsIrreducible(poly_five));

  EXPECT_TRUE(
      IsIrreducible(PolynomialDense<GaloisFieldPrime<uint8_t>>({one, one})));
  EXPECT_TRUE(
      IsIrreducible(PolynomialDense<GaloisFieldPrime<uint8_t>>({two, three})));
  EXPECT_TRUE(
      IsIrreducible(PolynomialDense<GaloisFieldPrime<uint8_t>>({zero, one})));
}

TEST_F(PolyUtilsTest, IsIrreducibleDegrees1To8OverGF2AndGF3) {
  // Known Conway / irreducible polynomials over GF(2) for degrees 1..8
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^2 + x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^3 + x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^4 + x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^5 + x^2 + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^6 + x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^7 + x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf2, "x^8 + x^4 + x^3 + x^2 + 1")));

  // Reducible polynomials over GF(2) (including products of irreducibles with
  // no roots in GF(2)!)
  EXPECT_FALSE(IsIrreducible(ParsePolynomial(gf2, "x^2 + 1")));
  EXPECT_FALSE(
      IsIrreducible(ParsePolynomial(gf2, "x^4 + x^2 + 1")));  // (x^2+x+1)^2
  EXPECT_FALSE(IsIrreducible(
      ParsePolynomial(gf2, "x^5 + x^4 + 1")));  // (x^2+x+1)(x^3+x+1)
  EXPECT_FALSE(IsIrreducible(ParsePolynomial(
      gf2, "x^6 + x^5 + x^4 + x^3 + x^2 + x + 1")));  // (x^3+x+1)(x^3+x^2+1)

  // GF(3): reducible and irreducible cubics
  auto gf3 = std::make_shared<GaloisFieldPrime<uint8_t>>(3);
  EXPECT_FALSE(IsIrreducible(ParsePolynomial(gf3, "x^3 - x")));
  EXPECT_FALSE(IsIrreducible(
      ParsePolynomial(gf3, "x^3 + x^2 + x + 1")));  // (x+1)(x^2+1)
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf3, "x^3 + 2x + 1")));
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf3, "x^3 + 2x + 2")));

  // GF(7): quadratics
  EXPECT_TRUE(IsIrreducible(ParsePolynomial(gf7, "x^2 + 1")));
  EXPECT_FALSE(IsIrreducible(ParsePolynomial(gf7, "x^2 - 1")));
}

TEST_F(PolyUtilsTest, ParseAndEvaluatePolynomial) {
  auto poly = ParsePolynomial(gf7, "2x^2 + 3x + 1");

  EXPECT_EQ(poly(two), one);
  EXPECT_EQ(poly(zero), one);
}

TEST_F(PolyUtilsTest, ParsePolynomialThenExtendedGcd) {
  auto poly1 = ParsePolynomial(gf7, "x^2 - 1");
  auto poly2 = ParsePolynomial(gf7, "x - 1");

  auto [gcd, bezout] = PolynomialDenseExtendedGcd(poly1, poly2);
  auto [s, t] = bezout;

  EXPECT_EQ(gcd.Degree(), 1);
  EXPECT_EQ(gcd[1], one);
  EXPECT_EQ(poly1 * s + poly2 * t, gcd);
}
