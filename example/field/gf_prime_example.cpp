

#include <iostream>
#include <memory>

#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Prime Field Examples\n";
  std::cout << "============================\n";
  std::cout
      << "This example uses GFP alias and direct constructors (NO factory)\n";

  std::cout << "\n=== Creating Prime Fields with GFP Alias ===\n";

  auto gf7_field = std::make_shared<GFP<uint8_t>>(7);
  std::cout << "Created GF(7) using GFP<uint8_t> alias\n";
  std::cout << "  Order: " << gf7_field->Order()
            << ", Characteristic: " << gf7_field->Characteristic() << "\n";

  auto gf11_field = std::make_shared<GFP<uint8_t>>(11);
  std::cout << "Created GF(11) using GFP<uint8_t> alias\n";
  std::cout << "  Order: " << gf11_field->Order()
            << ", Characteristic: " << gf11_field->Characteristic() << "\n";

  auto gf101_field = std::make_shared<GFP<uint8_t>>(101);
  std::cout << "Created GF(101) using GFP<uint8_t> alias\n";
  std::cout << "  Order: " << gf101_field->Order()
            << ", Characteristic: " << gf101_field->Characteristic() << "\n";

  std::cout << "\n=== Operations in GF(11) ===\n";

  using ElementType = GaloisFieldElementBase<GFP<uint8_t>>;
  ElementType a(3, gf11_field);
  ElementType b(5, gf11_field);
  ElementType c(7, gf11_field);

  std::cout << "Elements: a = " << static_cast<int>(a.Value())
            << ", b = " << static_cast<int>(b.Value())
            << ", c = " << static_cast<int>(c.Value()) << "\n";

  auto sum = a + b;
  auto diff = a - b;
  auto product = a * b;
  auto quotient = a / b;

  std::cout << "a + b = " << static_cast<int>(sum.Value()) << "\n";
  std::cout << "a - b = " << static_cast<int>(diff.Value()) << "\n";
  std::cout << "a * b = " << static_cast<int>(product.Value()) << "\n";
  std::cout << "a / b = " << static_cast<int>(quotient.Value()) << "\n";

  auto a_inv = a.Inv();
  std::cout << "a^(-1) = " << static_cast<int>(a_inv.Value()) << "\n";
  std::cout << "Verification: a * a^(-1) = "
            << static_cast<int>((a * a_inv).Value()) << " (should be 1)\n";

  auto a_squared = a * a;
  auto a_cubed = a_squared * a;
  std::cout << "a^2 = " << static_cast<int>(a_squared.Value()) << "\n";
  std::cout << "a^3 = " << static_cast<int>(a_cubed.Value()) << "\n";

  auto complex_expr = (a + b) * c - a * b;
  std::cout << "(a + b) * c - a * b = "
            << static_cast<int>(complex_expr.Value()) << "\n";

  std::cout << "\n=== Field Properties ===\n";
  std::cout << "GF(11) characteristics:\n";
  std::cout << "  Order: " << gf11_field->Order() << "\n";
  std::cout << "  Characteristic: " << gf11_field->Characteristic() << "\n";
  std::cout << "  Modulus: " << gf11_field->Modulus() << "\n";

  std::cout << "\n=== All Non-Zero Elements and Their Inverses ===\n";
  for (int i = 1; i < 11; ++i) {
    ElementType elem(i, gf11_field);
    auto inv = elem.Inv();
    std::cout << i << "^(-1) = " << static_cast<int>(inv.Value())
              << " (check: " << i << " * " << static_cast<int>(inv.Value())
              << " = " << static_cast<int>((elem * inv).Value()) << ")\n";
  }

  std::cout << "\nPrime field examples completed successfully!\n";
  return 0;
}
