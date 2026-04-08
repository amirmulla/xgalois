

#include <iostream>
#include <memory>
#include <variant>

#include "xgalois/field/gf_factory.hpp"

using namespace xg;

auto getFieldProperties = [](const GaloisFieldVariant& field) {
  return std::visit(
      [](const auto& field_ptr) -> std::pair<uint64_t, uint64_t> {
        return {field_ptr->Order(), field_ptr->Characteristic()};
      },
      field);
};

auto performSimpleOperation = [](const GaloisFieldVariant& field, uint64_t a,
                                 uint64_t b,
                                 const std::string& op) -> uint64_t {
  return std::visit(
      [a, b, &op](const auto& field_ptr) -> uint64_t {
        using FieldType = typename std::decay_t<decltype(*field_ptr)>;
        using ElementType = typename FieldType::element_type;

        if constexpr (std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint8_t>> ||
                      std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint16_t>> ||
                      std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint32_t>>) {
          return 0;
        } else {
          ElementType elem_a = static_cast<ElementType>(a);
          ElementType elem_b = static_cast<ElementType>(b);
          ElementType result;

          if (op == "add") {
            result = field_ptr->Add(elem_a, elem_b);
          } else if (op == "mul") {
            result = field_ptr->Mul(elem_a, elem_b);
          } else if (op == "div") {
            result = field_ptr->Div(elem_a, elem_b);
          } else if (op == "inv") {
            result = field_ptr->Inv(elem_a);
          } else {
            result = elem_a;
          }

          return static_cast<uint64_t>(result);
        }
      },
      field);
};

int main() {
  std::cout << "XGalois Factory Examples\n";
  std::cout << "========================\n";

  std::cout << "\n=== Using GF() Factory to Create Different Field Types ===\n";

  auto gf7 = GF(7);
  auto gf11 = GF(11);
  auto gf101 = GF(101);

  std::cout << "Created prime fields:\n";
  auto [order7, char7] = getFieldProperties(gf7);
  std::cout << "  GF(7) - Order: " << order7 << ", Characteristic: " << char7
            << "\n";

  auto [order11, char11] = getFieldProperties(gf11);
  std::cout << "  GF(11) - Order: " << order11 << ", Characteristic: " << char11
            << "\n";

  auto [order101, char101] = getFieldProperties(gf101);
  std::cout << "  GF(101) - Order: " << order101
            << ", Characteristic: " << char101 << "\n";

  auto gf4 = GF(4);
  auto gf8 = GF(8);
  auto gf16 = GF(16);

  std::cout << "\nCreated binary extension fields:\n";
  auto [order4, char4] = getFieldProperties(gf4);
  std::cout << "  GF(2^2) - Order: " << order4 << ", Characteristic: " << char4
            << "\n";

  auto [order8, char8] = getFieldProperties(gf8);
  std::cout << "  GF(2^3) - Order: " << order8 << ", Characteristic: " << char8
            << "\n";

  auto [order16, char16] = getFieldProperties(gf16);
  std::cout << "  GF(2^4) - Order: " << order16
            << ", Characteristic: " << char16 << "\n";

  auto gf9 = GF(9, "poly");
  auto gf25 = GF(25, "poly");

  std::cout << "\nCreated general extension fields:\n";
  auto [order9, char9] = getFieldProperties(gf9);
  std::cout << "  GF(3^2) - Order: " << order9 << ", Characteristic: " << char9
            << "\n";

  auto [order25, char25] = getFieldProperties(gf25);
  std::cout << "  GF(5^2) - Order: " << order25
            << ", Characteristic: " << char25 << "\n";

  std::cout << "\n=== Field Operations Using Factory-Created Fields ===\n";

  std::cout << "Operations in GF(11):\n";
  uint64_t a = 3, b = 5, c = 7;

  std::cout << "  a = " << a << ", b = " << b << ", c = " << c << "\n";
  std::cout << "  a + b = " << performSimpleOperation(gf11, a, b, "add")
            << "\n";
  std::cout << "  a * b = " << performSimpleOperation(gf11, a, b, "mul")
            << "\n";
  std::cout << "  a / b = " << performSimpleOperation(gf11, a, b, "div")
            << "\n";
  std::cout << "  a^(-1) = " << performSimpleOperation(gf11, a, 0, "inv")
            << "\n";

  std::cout << "\nOperations in GF(2^4):\n";
  uint64_t x = 3, y = 5;

  std::cout << "  x = " << x << ", y = " << y << "\n";
  std::cout << "  x + y = " << performSimpleOperation(gf16, x, y, "add")
            << "\n";
  std::cout << "  x * y = " << performSimpleOperation(gf16, x, y, "mul")
            << "\n";
  std::cout << "  x / y = " << performSimpleOperation(gf16, x, y, "div")
            << "\n";

  std::cout << "\n=== Different Field Representations ===\n";

  auto gf7_int = GF(7, "int");
  auto gf7_hex = GF(7, "hex");

  std::cout << "Same field with different representations:\n";
  std::cout << "  GF(7) with INT representation created\n";
  std::cout << "  GF(7) with HEX representation created\n";
  std::cout << "  Element operations work the same way regardless of "
               "representation\n";

  std::cout << "\n=== Using Factory Helper Functions ===\n";

  try {
    auto elem_variant = FetchElement(gf11, 5);
    std::cout << "Successfully created element with value 5 in GF(11) using "
                 "FetchElement\n";

    std::cout << "Note: Extension fields (like GF(9), GF(25)) require "
                 "polynomial element creation\n";
    std::cout << "      and should use field methods directly rather than "
                 "FetchElement\n";
  } catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << "\n";
  }

  std::cout << "\n=== Factory Pattern Benefits ===\n";
  std::cout << "The GF() factory provides several advantages:\n";
  std::cout << "1. Unified interface for creating different field types\n";
  std::cout << "2. Automatic optimal implementation selection\n";
  std::cout << "3. Database integration for irreducible polynomials\n";
  std::cout << "4. Type-safe variant container for runtime polymorphism\n";
  std::cout << "5. Memory-efficient element type selection\n";

  std::cout << "\nFactory examples completed successfully!\n";
  std::cout << "Note: This is the ONLY example that uses GF() factory.\n";
  std::cout << "All other examples use direct constructors and type aliases.\n";
  std::cout << "The factory returns a std::variant that requires std::visit "
               "for operations.\n";

  return 0;
}
