

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

  std::cout << "\n=== Creating Base Binary Field GF(2) ===\n";
  auto gf2_field = std::make_shared<GF2>();
  std::cout << "Created GF(2) using GF2 alias\n";
  std::cout << "  Order: " << gf2_field->Order()
            << ", Characteristic: " << gf2_field->Characteristic() << "\n";

  std::cout << "\n=== Operations in GF(2) ===\n";
  uint8_t a = 1, b = 1, c = 0;
  std::cout << "a = " << static_cast<int>(a) << ", b = " << static_cast<int>(b)
            << ", c = " << static_cast<int>(c) << "\n";
  std::cout << "a + b = " << static_cast<int>(gf2_field->Add(a, b))
            << " (XOR in GF(2))\n";
  std::cout << "a * b = " << static_cast<int>(gf2_field->Mul(a, b)) << "\n";
  std::cout << "a + c = " << static_cast<int>(gf2_field->Add(a, c)) << "\n";

  std::cout << "\n=== Creating Binary Extension Fields with GF2X Alias ===\n";

  auto gf4_field = std::make_shared<GF2X<uint8_t>>(2);
  std::cout << "Created GF(2^2) using GF2X<uint8_t> alias\n";
  std::cout << "  Order: " << gf4_field->Order()
            << ", Characteristic: " << gf4_field->Characteristic() << "\n";

  auto gf8_field = std::make_shared<GF2X<uint8_t>>(3);
  std::cout << "Created GF(2^3) using GF2X<uint8_t> alias\n";
  std::cout << "  Order: " << gf8_field->Order()
            << ", Characteristic: " << gf8_field->Characteristic() << "\n";

  auto gf16_field = std::make_shared<GF2X<uint8_t>>(4);
  std::cout << "Created GF(2^4) using GF2X<uint8_t> alias\n";
  std::cout << "  Order: " << gf16_field->Order()
            << ", Characteristic: " << gf16_field->Characteristic() << "\n";

  std::cout << "\n=== Operations in GF(2^4) ===\n";

  uint8_t x = 3;
  uint8_t y = 5;
  uint8_t z = 9;

  std::cout << "Elements in polynomial representation:\n";
  std::cout << "  x = " << static_cast<int>(x)
            << " (binary: " << std::bitset<4>(x) << ", polynomial: x + 1)\n";
  std::cout << "  y = " << static_cast<int>(y)
            << " (binary: " << std::bitset<4>(y) << ", polynomial: x² + 1)\n";
  std::cout << "  z = " << static_cast<int>(z)
            << " (binary: " << std::bitset<4>(z) << ", polynomial: x³ + 1)\n";

  auto sum = gf16_field->Add(x, y);
  auto product = gf16_field->Mul(x, y);
  auto quotient = gf16_field->Div(x, y);
  auto x_inv = gf16_field->Inv(x);

  std::cout << "\nArithmetic operations:\n";
  std::cout << "  x + y = " << static_cast<int>(sum)
            << " (binary: " << std::bitset<4>(sum) << ")\n";
  std::cout << "  x * y = " << static_cast<int>(product)
            << " (binary: " << std::bitset<4>(product) << ")\n";
  std::cout << "  x / y = " << static_cast<int>(quotient)
            << " (binary: " << std::bitset<4>(quotient) << ")\n";
  std::cout << "  x^(-1) = " << static_cast<int>(x_inv)
            << " (binary: " << std::bitset<4>(x_inv) << ")\n";

  auto verification = gf16_field->Mul(x, x_inv);
  std::cout << "  Verification: x * x^(-1) = " << static_cast<int>(verification)
            << " (should be 1)\n";

  std::cout << "\n=== Field Properties ===\n";
  std::cout << "GF(2^4) characteristics:\n";
  std::cout << "  Extension degree: " << static_cast<int>(gf16_field->Degree())
            << "\n";
  std::cout << "  Primitive element: "
            << static_cast<int>(gf16_field->MultiplicativeGenerator()) << "\n";

  std::cout << "\n=== All Elements using Primitive Element ===\n";
  auto primitive = gf16_field->MultiplicativeGenerator();
  uint8_t current = 1;

  std::cout << "Powers of primitive element " << static_cast<int>(primitive)
            << ":\n";
  for (int i = 0; i < 15; ++i) {
    std::cout << "  " << static_cast<int>(primitive) << "^" << i << " = "
              << static_cast<int>(current)
              << " (binary: " << std::bitset<4>(current) << ")\n";
    current = gf16_field->Mul(current, primitive);
  }

  std::cout << "\nBinary field examples completed successfully!\n";
  std::cout << "Note: This example uses direct constructors and type aliases "
               "(NO factory).\n";

  return 0;
}
