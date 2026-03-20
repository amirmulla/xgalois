/**
 * @file gf_extension_example.cpp
 * @brief Examples for extension field operations using GFPX alias and direct
 * constructors
 */

#include <iostream>
#include <memory>

#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Extension Field Examples\n";
  std::cout << "================================\n";
  std::cout
      << "This example uses GFPX alias and direct constructors (NO factory)\n";

  // Create extension fields using GFPX alias
  std::cout << "\n=== Creating Extension Fields with GFPX Alias ===\n";

  // Create GF(3^2) = GF(9) using GFPX alias
  auto gf9_field = std::make_shared<GFPX<uint8_t>>(std::make_pair(3, 2));
  std::cout << "Created GF(3^2) = GF(9) using GFPX<uint8_t> alias\n";
  std::cout << "  Order: " << gf9_field->Order()
            << ", Characteristic: " << gf9_field->Characteristic() << "\n";

  // Create GF(5^2) = GF(25) using GFPX alias
  auto gf25_field = std::make_shared<GFPX<uint8_t>>(std::make_pair(5, 2));
  std::cout << "Created GF(5^2) = GF(25) using GFPX<uint8_t> alias\n";
  std::cout << "  Order: " << gf25_field->Order()
            << ", Characteristic: " << gf25_field->Characteristic() << "\n";

  // Create GF(7^2) = GF(49) using GFPX alias
  auto gf49_field = std::make_shared<GFPX<uint8_t>>(std::make_pair(7, 2));
  std::cout << "Created GF(7^2) = GF(49) using GFPX<uint8_t> alias\n";
  std::cout << "  Order: " << gf49_field->Order()
            << ", Characteristic: " << gf49_field->Characteristic() << "\n";

  // Demonstrate working with polynomial elements
  std::cout << "\n=== Extension Field Elements (Polynomials) ===\n";
  std::cout
      << "Elements in extension fields are polynomials over the base field.\n";
  std::cout << "For example, in GF(3^2), elements are polynomials a₀ + a₁x "
               "where a₀, a₁ ∈ GF(3).\n";

  // Create polynomial elements using field methods
  std::cout << "\n=== Creating Polynomial Elements ===\n";
  std::cout << "In extension fields, we work with polynomial elements created "
               "through field operations.\n";

  // Get basic field elements
  auto poly_zero = gf9_field->AdditiveIdentity();
  auto poly_one = gf9_field->MultiplicativeIdentity();

  std::cout << "Zero polynomial: " << gf9_field->ToString(poly_zero) << "\n";
  std::cout << "Unit polynomial: " << gf9_field->ToString(poly_one) << "\n";

  // Generate some random elements to demonstrate
  auto random_elem1 = gf9_field->Random();
  auto random_elem2 = gf9_field->Random();

  std::cout << "Random element 1: " << gf9_field->ToString(random_elem1)
            << "\n";
  std::cout << "Random element 2: " << gf9_field->ToString(random_elem2)
            << "\n";

  // Demonstrate field operations
  std::cout << "\n=== Field Operations ===\n";

  // Addition
  auto sum = gf9_field->Add(random_elem1, poly_one);
  std::cout << "random1 + 1 = " << gf9_field->ToString(sum) << "\n";

  // Multiplication
  auto product = gf9_field->Mul(random_elem1, random_elem2);
  std::cout << "random1 * random2 = " << gf9_field->ToString(product) << "\n";

  // Demonstrate field properties
  std::cout << "\n=== Field Properties ===\n";
  std::cout << "GF(3^2) characteristics:\n";
  std::cout << "  Extension degree: " << gf9_field->Modulus().Degree() << "\n";
  std::cout << "  Base field characteristic: " << gf9_field->Characteristic()
            << "\n";

  // Show irreducible polynomial (modulus)
  std::cout << "  Irreducible polynomial: "
            << gf9_field->ToString(gf9_field->Modulus()) << "\n";

  // Random element
  auto random_elem = gf9_field->Random();
  std::cout << "  Random element: " << gf9_field->ToString(random_elem) << "\n";

  // Demonstrate multiplicative group
  std::cout << "\n=== Multiplicative Group Structure ===\n";
  std::cout << "The multiplicative group of GF(9) has order 8 (9-1).\n";
  std::cout << "Note: Primitive element computation is not implemented in this "
               "version.\n";

  // Instead, demonstrate powers of a random non-zero element
  auto non_zero_elem = gf9_field->Random();
  // Make sure it's not zero
  while (non_zero_elem == gf9_field->AdditiveIdentity()) {
    non_zero_elem = gf9_field->Random();
  }

  std::cout << "Powers of element: " << gf9_field->ToString(non_zero_elem)
            << "\n";
  auto current = gf9_field->MultiplicativeIdentity();  // Start with 1
  for (int i = 0; i < 8; ++i) {
    std::cout << "  elem^" << i << " = " << gf9_field->ToString(current)
              << "\n";
    if (i < 7) {  // Don't multiply after the last iteration
      current = gf9_field->Mul(current, non_zero_elem);
    }
  }

  // Demonstrate another extension field
  std::cout << "\n=== Operations in GF(5^2) = GF(25) ===\n";

  // Create random elements in GF(25) to demonstrate operations
  auto elem_2 = gf25_field->Random();
  auto elem_3x_plus_1 = gf25_field->Random();

  std::cout << "Working with random elements in GF(25):\n";
  std::cout << "Element 1: " << gf25_field->ToString(elem_2) << "\n";
  std::cout << "Element 2: " << gf25_field->ToString(elem_3x_plus_1) << "\n";

  auto sum25 = gf25_field->Add(elem_2, elem_3x_plus_1);
  auto product25 = gf25_field->Mul(elem_2, elem_3x_plus_1);

  std::cout << "Sum: " << gf25_field->ToString(sum25) << "\n";
  std::cout << "Product: " << gf25_field->ToString(product25) << "\n";

  // Show field properties
  std::cout << "GF(25) characteristics:\n";
  std::cout << "  Extension degree: " << gf25_field->Modulus().Degree() << "\n";
  std::cout << "  Base field characteristic: " << gf25_field->Characteristic()
            << "\n";
  std::cout << "  Irreducible polynomial: "
            << gf25_field->ToString(gf25_field->Modulus()) << "\n";

  std::cout << "\nExtension field examples completed successfully!\n";
  std::cout << "Note: This example uses direct constructors and type aliases "
               "(NO factory).\n";
  std::cout
      << "Extension fields use polynomial arithmetic over the base field.\n";

  return 0;
}
