/**
 * @file gf_binary_example.cpp
 * @brief Examples for binary field operations using GF2, GF2X aliases and
 * direct constructors
 */

#include <bitset>
#include <iostream>
#include <memory>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Binary Field Examples\n";
  std::cout << "=============================\n";
  std::cout << "This example uses GF2, GF2X aliases and direct constructors "
               "(NO factory)\n";

  // Create the base binary field GF(2)
  std::cout << "\n=== Creating Base Binary Field GF(2) ===\n";
  auto gf2_field = std::make_shared<GF2>();
  std::cout << "Created GF(2) using GF2 alias\n";
  std::cout << "  Order: " << gf2_field->Order()
            << ", Characteristic: " << gf2_field->Characteristic() << "\n";

  // Demonstrate basic GF(2) operations
  std::cout << "\n=== Operations in GF(2) ===\n";
  uint8_t a = 1, b = 1, c = 0;
  std::cout << "a = " << (int)a << ", b = " << (int)b << ", c = " << (int)c
            << "\n";
  std::cout << "a + b = " << (int)gf2_field->Add(a, b) << " (XOR in GF(2))\n";
  std::cout << "a * b = " << (int)gf2_field->Mul(a, b) << "\n";
  std::cout << "a + c = " << (int)gf2_field->Add(a, c) << "\n";

  // Create binary extension fields using GF2X alias
  std::cout << "\n=== Creating Binary Extension Fields with GF2X Alias ===\n";

  // GF(2^2) using GF2X alias
  auto gf4_field = std::make_shared<GF2X<uint8_t>>(2);  // Degree 2
  std::cout << "Created GF(2^2) using GF2X<uint8_t> alias\n";
  std::cout << "  Order: " << gf4_field->Order()
            << ", Characteristic: " << gf4_field->Characteristic() << "\n";

  // GF(2^3) using GF2X alias
  auto gf8_field = std::make_shared<GF2X<uint8_t>>(3);  // Degree 3
  std::cout << "Created GF(2^3) using GF2X<uint8_t> alias\n";
  std::cout << "  Order: " << gf8_field->Order()
            << ", Characteristic: " << gf8_field->Characteristic() << "\n";

  // GF(2^4) using GF2X alias
  auto gf16_field = std::make_shared<GF2X<uint8_t>>(4);  // Degree 4
  std::cout << "Created GF(2^4) using GF2X<uint8_t> alias\n";
  std::cout << "  Order: " << gf16_field->Order()
            << ", Characteristic: " << gf16_field->Characteristic() << "\n";

  // Demonstrate operations in GF(2^4)
  std::cout << "\n=== Operations in GF(2^4) ===\n";

  // Elements in GF(2^4) represented as polynomials over GF(2)
  // For example: 3 = 0011₂ represents polynomial x + 1
  //              5 = 0101₂ represents polynomial x² + 1
  //              9 = 1001₂ represents polynomial x³ + 1

  uint8_t x = 3;  // x + 1
  uint8_t y = 5;  // x² + 1
  uint8_t z = 9;  // x³ + 1

  std::cout << "Elements in polynomial representation:\n";
  std::cout << "  x = " << (int)x << " (binary: " << std::bitset<4>(x)
            << ", polynomial: x + 1)\n";
  std::cout << "  y = " << (int)y << " (binary: " << std::bitset<4>(y)
            << ", polynomial: x² + 1)\n";
  std::cout << "  z = " << (int)z << " (binary: " << std::bitset<4>(z)
            << ", polynomial: x³ + 1)\n";

  auto sum = gf16_field->Add(x, y);
  auto product = gf16_field->Mul(x, y);
  auto quotient = gf16_field->Div(x, y);
  auto x_inv = gf16_field->Inv(x);

  std::cout << "\nArithmetic operations:\n";
  std::cout << "  x + y = " << (int)sum << " (binary: " << std::bitset<4>(sum)
            << ")\n";
  std::cout << "  x * y = " << (int)product
            << " (binary: " << std::bitset<4>(product) << ")\n";
  std::cout << "  x / y = " << (int)quotient
            << " (binary: " << std::bitset<4>(quotient) << ")\n";
  std::cout << "  x^(-1) = " << (int)x_inv
            << " (binary: " << std::bitset<4>(x_inv) << ")\n";

  // Verify multiplicative inverse
  auto verification = gf16_field->Mul(x, x_inv);
  std::cout << "  Verification: x * x^(-1) = " << (int)verification
            << " (should be 1)\n";

  // Demonstrate field properties
  std::cout << "\n=== Field Properties ===\n";
  std::cout << "GF(2^4) characteristics:\n";
  std::cout << "  Extension degree: " << (int)gf16_field->Degree() << "\n";
  std::cout << "  Primitive element: "
            << (int)gf16_field->MultiplicativeGenerator() << "\n";

  // Generate all elements using primitive element
  std::cout << "\n=== All Elements using Primitive Element ===\n";
  auto primitive = gf16_field->MultiplicativeGenerator();
  uint8_t current = 1;  // Start with 1

  std::cout << "Powers of primitive element " << (int)primitive << ":\n";
  for (int i = 0; i < 15; ++i) {  // 2^4 - 1 = 15 non-zero elements
    std::cout << "  " << (int)primitive << "^" << i << " = " << (int)current
              << " (binary: " << std::bitset<4>(current) << ")\n";
    current = gf16_field->Mul(current, primitive);
  }

  std::cout << "\nBinary field examples completed successfully!\n";
  std::cout << "Note: This example uses direct constructors and type aliases "
               "(NO factory).\n";

  return 0;
}
