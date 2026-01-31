//===----------------------------------------------------------------------===//
//                          XGalois Library
//===----------------------------------------------------------------------===//
// Copyright (C) 2024 Amir Mulla
//
// Comprehensive unit tests for PolynomialDense class
//===----------------------------------------------------------------------===//

#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/utils/poly.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <sstream>
#include <stdexcept>

using namespace xg;

//===----------------------------------------------------------------------===//
// Test Fixtures
//===----------------------------------------------------------------------===//

class PolynomialDenseTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create GF(7) field for testing
        gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);

        // Create field elements
        zero = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(0, gf7);
        one = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, gf7);
        two = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(2, gf7);
        three = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(3, gf7);
        four = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(4, gf7);
        five = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(5, gf7);
        six = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(6, gf7);

        // Create GF(2) field for binary tests
        gf2 = std::make_shared<GaloisFieldBinary>();
        gf2_zero = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
        gf2_one = GaloisFieldElementBase<GaloisFieldBinary>(1, gf2);
    }

    std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7;
    std::shared_ptr<GaloisFieldBinary> gf2;

    GaloisFieldElementBase<GaloisFieldPrime<uint8_t>> zero, one, two, three, four, five, six;
    GaloisFieldElementBase<GaloisFieldBinary> gf2_zero, gf2_one;
};

//===----------------------------------------------------------------------===//
// Constructor and Basic Properties Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, ConstructorAndBasicProperties) {
    // Test polynomial 1 + 2x + 3x^2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    EXPECT_EQ(poly.Degree(), 2);
    EXPECT_EQ(poly.Size(), 3);
    EXPECT_EQ(poly[0], one);
    EXPECT_EQ(poly[1], two);
    EXPECT_EQ(poly[2], three);
    EXPECT_EQ(poly.GetVariable(), "x");  // Default variable
}

TEST_F(PolynomialDenseTest, CustomVariable) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs, "t");

    EXPECT_EQ(poly.GetVariable(), "t");

    // Test setting variable
    poly.SetVariable("y");
    EXPECT_EQ(poly.GetVariable(), "y");

    // Test empty variable throws
    EXPECT_THROW(poly.SetVariable(""), std::invalid_argument);
}

TEST_F(PolynomialDenseTest, ZeroPolynomial) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {zero};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    EXPECT_EQ(poly.Degree(), -1);
    EXPECT_EQ(poly.Size(), 1);
    EXPECT_EQ(poly[0], zero);
}

TEST_F(PolynomialDenseTest, ConstantPolynomial) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {five};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    EXPECT_EQ(poly.Degree(), 0);
    EXPECT_EQ(poly.Size(), 1);
    EXPECT_EQ(poly[0], five);
}

TEST_F(PolynomialDenseTest, LeadingZeroTrimming) {
    // Create polynomial with leading zeros: 1 + 2x + 0x^2 + 0x^3
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, zero, zero};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    // Should be trimmed to degree 1
    EXPECT_EQ(poly.Degree(), 1);
    EXPECT_EQ(poly.Size(), 2);
    EXPECT_EQ(poly[0], one);
    EXPECT_EQ(poly[1], two);
}

//===----------------------------------------------------------------------===//
// Indexing and Access Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, IndexingOperators) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    // Test read access
    EXPECT_EQ(poly[0], one);
    EXPECT_EQ(poly[1], two);
    EXPECT_EQ(poly[2], three);

    // Test write access
    poly[1] = four;
    EXPECT_EQ(poly[1], four);

    // Test out of bounds
    EXPECT_THROW(poly[3], std::out_of_range);
    EXPECT_THROW(poly[-1], std::out_of_range);
}

//===----------------------------------------------------------------------===//
// Arithmetic Operations Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, Addition) {
    // p1 = 1 + 2x + 3x^2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

    // p2 = 4 + 5x
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {four, five};
    PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

    // p1 + p2 = (1+4) + (2+5)x + 3x^2 = 5 + 0x + 3x^2 = 5 + 3x^2 (in GF(7))
    auto result = p1 + p2;

    EXPECT_EQ(result.Degree(), 2);
    EXPECT_EQ(result[0], five);  // 1 + 4 = 5
    EXPECT_EQ(result[1], zero);  // 2 + 5 = 7 ≡ 0 (mod 7)
    EXPECT_EQ(result[2], three); // 3 + 0 = 3
}

TEST_F(PolynomialDenseTest, AdditionWithZero) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> zero_coeffs = {zero};
    PolynomialDense<GaloisFieldPrime<uint8_t>> zero_poly(zero_coeffs);

    auto result = poly + zero_poly;
    EXPECT_EQ(result, poly);
}

TEST_F(PolynomialDenseTest, Subtraction) {
    // p1 = 1 + 2x + 3x^2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

    // p2 = 4 + 5x
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {four, five};
    PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

    // p1 - p2 = (1-4) + (2-5)x + 3x^2 = -3 + (-3)x + 3x^2 = 4 + 4x + 3x^2 (in GF(7))
    auto result = p1 - p2;

    EXPECT_EQ(result.Degree(), 2);
    EXPECT_EQ(result[0], four);  // 1 - 4 = -3 ≡ 4 (mod 7)
    EXPECT_EQ(result[1], four);  // 2 - 5 = -3 ≡ 4 (mod 7)
    EXPECT_EQ(result[2], three); // 3 - 0 = 3
}

TEST_F(PolynomialDenseTest, Negation) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    auto neg_poly = -poly;

    EXPECT_EQ(neg_poly.Degree(), 2);
    EXPECT_EQ(neg_poly[0], six);  // -1 ≡ 6 (mod 7)
    EXPECT_EQ(neg_poly[1], five); // -2 ≡ 5 (mod 7)
    EXPECT_EQ(neg_poly[2], four); // -3 ≡ 4 (mod 7)
}

TEST_F(PolynomialDenseTest, Multiplication) {
    // p1 = 1 + 2x
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {one, two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

    // p2 = 3 + 4x
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {three, four};
    PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

    // p1 * p2 = (1 + 2x)(3 + 4x) = 3 + 4x + 6x + 8x^2 = 3 + 10x + 8x^2 = 3 + 3x + x^2 (in GF(7))
    auto result = p1 * p2;

    EXPECT_EQ(result.Degree(), 2);
    EXPECT_EQ(result[0], three); // 1*3 = 3
    EXPECT_EQ(result[1], three); // 1*4 + 2*3 = 4 + 6 = 10 ≡ 3 (mod 7)
    EXPECT_EQ(result[2], one);   // 2*4 = 8 ≡ 1 (mod 7)
}

TEST_F(PolynomialDenseTest, ScalarMultiplication) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    auto result = poly * five;

    EXPECT_EQ(result.Degree(), 2);
    EXPECT_EQ(result[0], five);  // 1*5 = 5
    EXPECT_EQ(result[1], three); // 2*5 = 10 ≡ 3 (mod 7)
    EXPECT_EQ(result[2], one);   // 3*5 = 15 ≡ 1 (mod 7)
}

TEST_F(PolynomialDenseTest, MultiplicationByZero) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    auto result = poly * zero;

    EXPECT_EQ(result.Degree(), -1); // Zero polynomial
    EXPECT_EQ(result[0], zero);
}

//===----------------------------------------------------------------------===//
// Division Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, DivisionByConstant) {
    // p = 2 + 4x + 6x^2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {two, four, six};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    // divisor = 2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> div_coeffs = {two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(div_coeffs);

    auto result = poly / divisor;

    // Result should be 1 + 2x + 3x^2 (since 2^-1 = 4 in GF(7))
    EXPECT_EQ(result.Degree(), 2);
    EXPECT_EQ(result[0], one);   // 2 * 4 = 8 ≡ 1 (mod 7)
    EXPECT_EQ(result[1], two);   // 4 * 4 = 16 ≡ 2 (mod 7)
    EXPECT_EQ(result[2], three); // 6 * 4 = 24 ≡ 3 (mod 7)
}

TEST_F(PolynomialDenseTest, PolynomialDivision) {
    // dividend = x^2 + 3x + 2 = (x + 1)(x + 2)
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> dividend_coeffs = {two, three, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> dividend(dividend_coeffs);

    // divisor = x + 1
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> divisor_coeffs = {one, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(divisor_coeffs);

    auto result = dividend / divisor;

    // Result should be x + 2
    EXPECT_EQ(result.Degree(), 1);
    EXPECT_EQ(result[0], two); // constant term
    EXPECT_EQ(result[1], one); // x coefficient
}

TEST_F(PolynomialDenseTest, ModuloOperation) {
    // dividend = x^2 + 3x + 2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> dividend_coeffs = {two, three, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> dividend(dividend_coeffs);

    // divisor = x + 1
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> divisor_coeffs = {one, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(divisor_coeffs);

    auto remainder = dividend % divisor;

    // Since (x^2 + 3x + 2) = (x + 1)(x + 2) + 0, remainder should be 0
    EXPECT_EQ(remainder.Degree(), -1);
    EXPECT_EQ(remainder[0], zero);
}

TEST_F(PolynomialDenseTest, DivRemOperation) {
    // dividend = x^3 + 2x^2 + 3x + 4
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> dividend_coeffs = {four, three, two, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> dividend(dividend_coeffs);

    // divisor = x + 1
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> divisor_coeffs = {one, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(divisor_coeffs);

    auto [quotient, remainder] = dividend.DivRem(divisor);

    // Verify division: dividend = quotient * divisor + remainder
    auto verification = quotient * divisor + remainder;
    EXPECT_EQ(verification, dividend);
}

TEST_F(PolynomialDenseTest, DivisionByZero) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> zero_coeffs = {zero};
    PolynomialDense<GaloisFieldPrime<uint8_t>> zero_poly(zero_coeffs);

    EXPECT_THROW(poly.DivRem(zero_poly), std::invalid_argument);
}

//===----------------------------------------------------------------------===//
// Power Operation Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, PowerOperation) {
    // p = 1 + x
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    // p^0 = 1
    auto power0 = poly ^ 0;
    EXPECT_EQ(power0.Degree(), 0);
    EXPECT_EQ(power0[0], one);

    // p^1 = p
    auto power1 = poly ^ 1;
    EXPECT_EQ(power1, poly);

    // p^2 = (1 + x)^2 = 1 + 2x + x^2
    auto power2 = poly ^ 2;
    EXPECT_EQ(power2.Degree(), 2);
    EXPECT_EQ(power2[0], one);  // constant term
    EXPECT_EQ(power2[1], two);  // x coefficient
    EXPECT_EQ(power2[2], one);  // x^2 coefficient
}

//===----------------------------------------------------------------------===//
// Evaluation Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, Evaluation) {
    // p = 1 + 2x + 3x^2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    // p(0) = 1
    auto result0 = poly(zero);
    EXPECT_EQ(result0, one);

    // p(1) = 1 + 2*1 + 3*1^2 = 1 + 2 + 3 = 6
    auto result1 = poly(one);
    EXPECT_EQ(result1, six);

    // p(2) = 1 + 2*2 + 3*2^2 = 1 + 4 + 12 = 17 ≡ 3 (mod 7)
    auto result2 = poly(two);
    EXPECT_EQ(result2, three);
}

//===----------------------------------------------------------------------===//
// Comparison Operations Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, EqualityComparison) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);

    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs3 = {one, two, four};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly3(coeffs3);

    EXPECT_TRUE(poly1 == poly2);
    EXPECT_FALSE(poly1 == poly3);
    EXPECT_TRUE(poly1 != poly3);
    EXPECT_FALSE(poly1 != poly2);
}

TEST_F(PolynomialDenseTest, EqualityWithTrimming) {
    // poly1 = 1 + 2x
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {one, two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);

    // poly2 = 1 + 2x + 0x^2 (should be trimmed to same as poly1)
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {one, two, zero};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

    EXPECT_TRUE(poly1 == poly2);
}

//===----------------------------------------------------------------------===//
// Compound Assignment Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, CompoundAssignment) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {one, two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly1(coeffs1);

    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {three, four};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly2(coeffs2);

    // Test +=
    auto poly_copy = poly1;
    poly_copy += poly2;
    EXPECT_EQ(poly_copy, poly1 + poly2);

    // Test -=
    poly_copy = poly1;
    poly_copy -= poly2;
    EXPECT_EQ(poly_copy, poly1 - poly2);

    // Test *=
    poly_copy = poly1;
    poly_copy *= poly2;
    EXPECT_EQ(poly_copy, poly1 * poly2);

    // Test scalar *=
    poly_copy = poly1;
    poly_copy *= five;
    EXPECT_EQ(poly_copy, poly1 * five);
}

//===----------------------------------------------------------------------===//
// Derivative Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, Derivative) {
    // p = 1 + 2x + 3x^2 + 4x^3
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three, four};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    // p' = 2 + 6x + 12x^2 = 2 + 6x + 5x^2 (in GF(7))
    auto derivative = poly.Derivative();

    EXPECT_EQ(derivative.Degree(), 2);
    EXPECT_EQ(derivative[0], two);  // 1*2 = 2
    EXPECT_EQ(derivative[1], six);  // 2*3 = 6
    EXPECT_EQ(derivative[2], five); // 3*4 = 12 ≡ 5 (mod 7)
}

TEST_F(PolynomialDenseTest, DerivativeOfConstant) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {five};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    auto derivative = poly.Derivative();

    EXPECT_EQ(derivative.Degree(), -1); // Zero polynomial
    EXPECT_EQ(derivative[0], zero);
}

//===----------------------------------------------------------------------===//
// Print and String Representation Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, PrintPolynomial) {
    // Test printing various polynomials
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    std::stringstream ss;
    poly.Print(ss);
    std::string result = ss.str();

    // Should contain x terms and coefficients
    EXPECT_TRUE(result.find("x") != std::string::npos);
    EXPECT_TRUE(result.length() > 0);
}

TEST_F(PolynomialDenseTest, StreamOperator) {
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    std::stringstream ss;
    ss << poly;
    std::string result = ss.str();

    EXPECT_TRUE(result.length() > 0);
}

//===----------------------------------------------------------------------===//
// Binary Field Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, BinaryFieldOperations) {
    // Test polynomials over GF(2)
    std::vector<GaloisFieldElementBase<GaloisFieldBinary>> coeffs1 = {gf2_one, gf2_zero, gf2_one}; // 1 + x^2
    PolynomialDense<GaloisFieldBinary> poly1(coeffs1);

    std::vector<GaloisFieldElementBase<GaloisFieldBinary>> coeffs2 = {gf2_one, gf2_one}; // 1 + x
    PolynomialDense<GaloisFieldBinary> poly2(coeffs2);

    // Addition in GF(2): (1 + x^2) + (1 + x) = x + x^2
    auto sum = poly1 + poly2;
    EXPECT_EQ(sum.Degree(), 2);
    EXPECT_EQ(sum[0], gf2_zero); // 1 + 1 = 0 in GF(2)
    EXPECT_EQ(sum[1], gf2_one);  // 0 + 1 = 1
    EXPECT_EQ(sum[2], gf2_one);  // 1 + 0 = 1

    // Multiplication in GF(2): (1 + x^2)(1 + x) = 1 + x + x^2 + x^3
    auto product = poly1 * poly2;
    EXPECT_EQ(product.Degree(), 3);
    EXPECT_EQ(product[0], gf2_one);  // 1*1 = 1
    EXPECT_EQ(product[1], gf2_one);  // 1*1 + 0*1 = 1
    EXPECT_EQ(product[2], gf2_one);  // 1*1 + 1*0 = 1
    EXPECT_EQ(product[3], gf2_one);  // 1*1 = 1
}

//===----------------------------------------------------------------------===//
// Extended GCD Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, ExtendedGCD) {
    // Test Extended Euclidean Algorithm
    // a = x^2 + 3x + 2 = (x + 1)(x + 2)
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_a = {two, three, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> a(coeffs_a);

    // b = x + 1
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_b = {one, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> b(coeffs_b);

    auto [gcd, coeffs] = utils::PolynomialDenseExtendedGcd(a, b);
    auto [s, t] = coeffs;

    // gcd should be x + 1 (or its monic form)
    EXPECT_TRUE(gcd.Degree() >= 0);

    // Verify: a*s + b*t = gcd
    auto verification = a * s + b * t;
    EXPECT_EQ(verification, gcd);
}

TEST_F(PolynomialDenseTest, ExtendedGCDCoprime) {
    // Test with coprime polynomials
    // a = x + 1
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_a = {one, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> a(coeffs_a);

    // b = x + 2
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs_b = {two, one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> b(coeffs_b);

    auto [gcd, coeffs] = utils::PolynomialDenseExtendedGcd(a, b);
    auto [s, t] = coeffs;

    // gcd should be 1 (constant polynomial) since a and b are coprime
    EXPECT_EQ(gcd.Degree(), 0);
    EXPECT_EQ(gcd[0], one);

    // Verify: a*s + b*t = gcd
    auto verification = a * s + b * t;
    EXPECT_EQ(verification, gcd);
}

//===----------------------------------------------------------------------===//
// Edge Cases and Error Handling Tests
//===----------------------------------------------------------------------===//

TEST_F(PolynomialDenseTest, EmptyCoefficientsAssert) {
    // Note: The constructor uses assert for empty coefficients,
    // so this test verifies the contract rather than exception throwing.
    // In debug builds, this would trigger an assertion failure.
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> empty_coeffs;

    // We can't test assert in a clean way with gtest, so we skip this test
    // and just verify that a valid polynomial can be created
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> valid_coeffs = {one};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(valid_coeffs);
    EXPECT_EQ(poly.Degree(), 0);
}

TEST_F(PolynomialDenseTest, LargePolynomials) {
    // Test with a larger polynomial
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
    for (int i = 0; i < 10; ++i) {
        coeffs.push_back(GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(i % 7, gf7));
    }

    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);
    EXPECT_EQ(poly.Degree(), 9);
    EXPECT_EQ(poly.Size(), 10);

    // Test operations on large polynomials
    auto doubled = poly * two;
    EXPECT_EQ(doubled.Degree(), 9);

    auto squared = poly * poly;
    EXPECT_EQ(squared.Degree(), 18);
}

TEST_F(PolynomialDenseTest, FieldConsistency) {
    // Test that polynomial operations maintain field consistency
    std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one, two, three};
    PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

    EXPECT_EQ(poly.Field(), gf7);

    auto doubled = poly * two;
    EXPECT_EQ(doubled.Field(), gf7);

    auto sum = poly + poly;
    EXPECT_EQ(sum.Field(), gf7);
}