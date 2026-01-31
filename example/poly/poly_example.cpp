/**
 * @file poly_example.cpp
 * @brief Examples for polynomial operations over Galois fields using direct constructors
 */

#include <iostream>
#include <memory>
#include <vector>
#include <sstream>
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/poly/poly_dense.hpp"

using namespace xg;

// Helper function to convert polynomial to string
template<typename GaloisField>
std::string PolynomialToString(const PolynomialDense<GaloisField>& poly) {
    std::stringstream ss;
    poly.Print(ss);
    return ss.str();
}

int main() {
    std::cout << "XGalois Polynomial Examples\n";
    std::cout << "===========================\n";
    std::cout << "This example uses direct field constructors (NO factory)\n";

    // Create fields using direct constructors
    std::cout << "\n=== Creating Fields for Polynomial Examples ===\n";

    auto gf7_field = std::make_shared<GFP<uint8_t>>(7);
    std::cout << "Created GF(7) using direct GFP constructor\n";

    auto gf2_field = std::make_shared<GF2>();
    std::cout << "Created GF(2) using direct GF2 constructor\n";

    // Create polynomials over GF(7)
    std::cout << "\n=== Polynomials over GF(7) ===\n";

    // Create polynomial 3x^2 + 2x + 1 over GF(7)
    using ElementType7 = GaloisFieldElementBase<GFP<uint8_t>>;
    std::vector<ElementType7> coeffs1 = {
        ElementType7(1, gf7_field),  // constant term
        ElementType7(2, gf7_field),  // x coefficient
        ElementType7(3, gf7_field)   // x^2 coefficient
    };
    PolynomialDense<GFP<uint8_t>> poly1(coeffs1);
    std::cout << "Created polynomial p1(x) = " << PolynomialToString(poly1) << "\n";

    // Create polynomial 5x + 4 over GF(7)
    std::vector<ElementType7> coeffs2 = {
        ElementType7(4, gf7_field),  // constant term
        ElementType7(5, gf7_field)   // x coefficient
    };
    PolynomialDense<GFP<uint8_t>> poly2(coeffs2);
    std::cout << "Created polynomial p2(x) = " << PolynomialToString(poly2) << "\n";

    // Polynomial arithmetic
    std::cout << "\n=== Polynomial Arithmetic over GF(7) ===\n";

    auto sum = poly1 + poly2;
    std::cout << "p1 + p2 = " << PolynomialToString(sum) << "\n";

    auto difference = poly1 - poly2;
    std::cout << "p1 - p2 = " << PolynomialToString(difference) << "\n";

    auto product = poly1 * poly2;
    std::cout << "p1 * p2 = " << PolynomialToString(product) << "\n";

    // Polynomial evaluation
    std::cout << "\n=== Polynomial Evaluation ===\n";
    ElementType7 eval_point(2, gf7_field);
    auto result1 = poly1(eval_point);  // Use operator() instead of Evaluate
    auto result2 = poly2(eval_point);

    std::stringstream ss1, ss2;
    result1.Print(ss1);
    result2.Print(ss2);
    std::cout << "p1(2) = " << ss1.str() << "\n";
    std::cout << "p2(2) = " << ss2.str() << "\n";

    // Polynomial properties
    std::cout << "\n=== Polynomial Properties ===\n";
    std::cout << "Degree of p1: " << poly1.Degree() << "\n";
    std::cout << "Degree of p2: " << poly2.Degree() << "\n";
    std::cout << "Degree of p1*p2: " << product.Degree() << "\n";

    // Simple leading coefficient check (manually access the highest degree coefficient)
    std::stringstream ss_lead;
    poly1[poly1.Degree()].Print(ss_lead);
    std::cout << "Leading coefficient of p1: " << ss_lead.str() << "\n";
    std::cout << "Is p1 monic? " << (poly1[poly1.Degree()].Value() == 1 ? "Yes" : "No") << "\n";

    // Polynomials over GF(2)
    std::cout << "\n=== Polynomials over GF(2) ===\n";

    // Create polynomial x^3 + x + 1 over GF(2)
    using ElementType2 = GaloisFieldElementBase<GF2>;
    std::vector<ElementType2> bin_coeffs1 = {
        ElementType2(1, gf2_field),  // constant term: 1
        ElementType2(1, gf2_field),  // x coefficient: 1
        ElementType2(0, gf2_field),  // x^2 coefficient: 0
        ElementType2(1, gf2_field)   // x^3 coefficient: 1
    };
    PolynomialDense<GF2> bin_poly1(bin_coeffs1);
    std::cout << "Created binary polynomial q1(x) = " << PolynomialToString(bin_poly1) << "\n";

    // Create polynomial x^2 + 1 over GF(2)
    std::vector<ElementType2> bin_coeffs2 = {
        ElementType2(1, gf2_field),  // constant term: 1
        ElementType2(0, gf2_field),  // x coefficient: 0
        ElementType2(1, gf2_field)   // x^2 coefficient: 1
    };
    PolynomialDense<GF2> bin_poly2(bin_coeffs2);
    std::cout << "Created binary polynomial q2(x) = " << PolynomialToString(bin_poly2) << "\n";

    // Binary polynomial arithmetic
    std::cout << "\n=== Binary Polynomial Arithmetic ===\n";

    auto bin_sum = bin_poly1 + bin_poly2;
    std::cout << "q1 + q2 = " << PolynomialToString(bin_sum) << "\n";

    auto bin_product = bin_poly1 * bin_poly2;
    std::cout << "q1 * q2 = " << PolynomialToString(bin_product) << "\n";

    // Polynomial division
    std::cout << "\n=== Polynomial Division ===\n";

    // Divide product by one of the factors
    auto quotient = bin_product / bin_poly2;
    auto remainder = bin_product % bin_poly2;

    std::cout << "(q1 * q2) / q2 = " << PolynomialToString(quotient) << "\n";
    std::cout << "(q1 * q2) % q2 = " << PolynomialToString(remainder) << "\n";
    std::cout << "Verification: quotient should equal q1, remainder should be 0\n";

    // Polynomial derivatives
    std::cout << "\n=== Polynomial Derivatives ===\n";

    auto derivative1 = poly1.Derivative();
    std::cout << "Derivative of p1 = " << PolynomialToString(derivative1) << "\n";

    auto bin_derivative = bin_poly1.Derivative();
    std::cout << "Derivative of q1 = " << PolynomialToString(bin_derivative) << "\n";

    // Simple GCD demonstration using repeated subtraction approach
    std::cout << "\n=== Simple Polynomial Operations ===\n";

    // Create some polynomials for demonstration
    std::vector<ElementType2> simple_coeffs1 = {
        ElementType2(0, gf2_field),  // constant: 0
        ElementType2(1, gf2_field),  // x: 1
        ElementType2(1, gf2_field)   // x^2: 1
    };  // x + x^2 = x(1 + x)

    std::vector<ElementType2> simple_coeffs2 = {
        ElementType2(0, gf2_field),  // constant: 0
        ElementType2(0, gf2_field),  // x: 0
        ElementType2(1, gf2_field)   // x^2: 1
    };  // x^2

    PolynomialDense<GF2> simple_poly1(simple_coeffs1);
    PolynomialDense<GF2> simple_poly2(simple_coeffs2);

    std::cout << "Polynomial r1(x) = " << PolynomialToString(simple_poly1) << "\n";
    std::cout << "Polynomial r2(x) = " << PolynomialToString(simple_poly2) << "\n";

    // Note: Since GCD is not implemented, we'll show the polynomials and explain
    std::cout << "Note: GCD computation would show that gcd(r1, r2) = x\n";
    std::cout << "      since r1 = x(1 + x) and r2 = x^2 = x*x, so gcd = x\n";

    // Demonstrate some useful polynomial operations
    std::cout << "\n=== Additional Polynomial Operations ===\n";

    // Check if polynomial is zero
    std::vector<ElementType7> zero_coeffs = {ElementType7(0, gf7_field)};
    PolynomialDense<GFP<uint8_t>> zero_poly(zero_coeffs);
    std::cout << "Zero polynomial: " << PolynomialToString(zero_poly) << "\n";
    std::cout << "Is zero? " << (zero_poly.Degree() == -1 ? "Yes" : "No") << "\n";

    // Check polynomial equality
    auto poly2_copy = poly2;
    std::cout << "p2 == copy of p2? " << (poly2 == poly2_copy ? "Yes" : "No") << "\n";
    std::cout << "p1 == p2? " << (poly1 == poly2 ? "Yes" : "No") << "\n";

    std::cout << "\nPolynomial examples completed successfully!\n";
    std::cout << "Note: This example uses direct field constructors (NO factory).\n";
    std::cout << "Polynomials are represented with coefficients as field elements.\n";

    return 0;
}
