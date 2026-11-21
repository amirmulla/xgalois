#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include "xgalois/field/gf_factory.hpp"

using namespace xg;

namespace {

std::pair<uint64_t, uint64_t> GetFieldProperties(
    const GaloisFieldVariant& field) {
  return std::visit(
      [](const auto& field_ptr) -> std::pair<uint64_t, uint64_t> {
        return {field_ptr->Order(), field_ptr->Characteristic()};
      },
      field);
}

uint64_t PerformOperation(const GaloisFieldVariant& field, uint64_t a,
                          uint64_t b, const std::string& op) {
  return std::visit(
      [a, b, &op](const auto& field_ptr) -> uint64_t {
        using FieldType = std::decay_t<decltype(*field_ptr)>;
        using ElementType = typename FieldType::element_type;

        auto to_internal = [&](uint64_t v) -> ElementType {
          if constexpr (detail::kIsGaloisFieldExtension<FieldType>) {
            return field_ptr->FromInt(v);
          } else {
            return field_ptr->SetElementValue(static_cast<ElementType>(v));
          }
        };

        auto to_external = [&](const ElementType& v) -> uint64_t {
          if constexpr (detail::kIsGaloisFieldExtension<FieldType>) {
            return field_ptr->ToInt(v);
          } else {
            return static_cast<uint64_t>(field_ptr->GetElementValue(v));
          }
        };

        const ElementType elem_a = to_internal(a);
        const ElementType elem_b = to_internal(b);
        ElementType result{};

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

        return to_external(result);
      },
      field);
}

void PrintElementVariant(const std::string& label,
                         const GaloisFieldElementVariant& elem) {
  std::visit([&label](const auto& e) { std::cout << label << e << "\n"; },
             elem);
}

}  // namespace

int main() {
  std::cout << "XGalois Factory Examples\n";
  std::cout << "========================\n";

  std::cout << "\n=== Creating Prime Fields GF(p) ===\n";
  auto gf7 = GF(7);
  auto gf11_table = GF(11, "int", "", "α", "table");
  auto gf101 = GF(101);

  for (const auto& [name, field] :
       std::initializer_list<std::pair<const char*, const GaloisFieldVariant*>>{
           {"GF(7) [auto]", &gf7},
           {"GF(11) [table]", &gf11_table},
           {"GF(101) [auto]", &gf101}}) {
    auto [order, characteristic] = GetFieldProperties(*field);
    std::cout << "  " << name << " -> Order: " << order
              << ", Characteristic: " << characteristic << "\n";
  }

  std::cout << "\n=== Creating Binary Extension Fields GF(2^m) ===\n";
  auto gf16_log_opt =
      GF(std::pair<uint64_t, uint64_t>{2, 4}, "int", "", "α", "log_opt");
  auto gf16_zech =
      GF(std::pair<uint64_t, uint64_t>{2, 4}, "int", "", "α", "zech");
  auto gf256_auto = GF(std::pair<uint64_t, uint64_t>{2, 8});

  for (const auto& [name, field] :
       std::initializer_list<std::pair<const char*, const GaloisFieldVariant*>>{
           {"GF(2^4) [log_opt]", &gf16_log_opt},
           {"GF(2^4) [zech]", &gf16_zech},
           {"GF(2^8) [auto/log_opt]", &gf256_auto}}) {
    auto [order, characteristic] = GetFieldProperties(*field);
    std::cout << "  " << name << " -> Order: " << order
              << ", Characteristic: " << characteristic << "\n";
  }

  std::cout << "\n=== Creating General Extension Fields GF(p^m) ===\n";
  auto gf9 = GF(9, "poly");
  auto gf25 = GF(25, "poly");

  for (const auto& [name, field] :
       std::initializer_list<std::pair<const char*, const GaloisFieldVariant*>>{
           {"GF(3^2) [poly]", &gf9}, {"GF(5^2) [poly]", &gf25}}) {
    auto [order, characteristic] = GetFieldProperties(*field);
    std::cout << "  " << name << " -> Order: " << order
              << ", Characteristic: " << characteristic << "\n";
  }

  std::cout << "\n=== Arithmetic Across Factory-Created Fields ===\n";
  const uint64_t a = 3;
  const uint64_t b = 5;

  std::cout << "In GF(11) [table] (a=" << a << ", b=" << b << "):\n";
  std::cout << "  a + b  = " << PerformOperation(gf11_table, a, b, "add")
            << "\n";
  std::cout << "  a * b  = " << PerformOperation(gf11_table, a, b, "mul")
            << "\n";
  std::cout << "  a / b  = " << PerformOperation(gf11_table, a, b, "div")
            << "\n";
  std::cout << "  a^(-1) = " << PerformOperation(gf11_table, a, 0, "inv")
            << "\n";

  std::cout << "\nComparing GF(2^4) [log_opt] vs [zech] (a=" << a << ", b=" << b
            << "):\n";
  std::cout << "  log_opt: a + b = "
            << PerformOperation(gf16_log_opt, a, b, "add")
            << ", a * b = " << PerformOperation(gf16_log_opt, a, b, "mul")
            << ", a / b = " << PerformOperation(gf16_log_opt, a, b, "div")
            << "\n";
  std::cout << "  zech:    a + b = " << PerformOperation(gf16_zech, a, b, "add")
            << ", a * b = " << PerformOperation(gf16_zech, a, b, "mul")
            << ", a / b = " << PerformOperation(gf16_zech, a, b, "div") << "\n";

  std::cout << "\nIn GF(3^2) (base-3 int encodings a=4 -> α+1, b=5 -> α+2):\n";
  std::cout << "  a + b  = " << PerformOperation(gf9, 4, 5, "add") << "\n";
  std::cout << "  a * b  = " << PerformOperation(gf9, 4, 5, "mul") << "\n";
  std::cout << "  a^(-1) = " << PerformOperation(gf9, 4, 0, "inv") << "\n";

  std::cout << "\n=== Field Representations & FetchElement ===\n";
  auto gf7_int = GF(7, "int");
  auto gf7_hex = GF(7, "hex");

  PrintElementVariant("  GF(7)  [INT]  FetchElement(6) = ",
                      FetchElement(gf7_int, 6));
  PrintElementVariant("  GF(7)  [HEX]  FetchElement(6) = ",
                      FetchElement(gf7_hex, 6));
  PrintElementVariant("  GF(16) [ZECH] FetchElement(5) = ",
                      FetchElement(gf16_zech, 5));
  PrintElementVariant("  GF(9)  [POLY] FetchElement(5) = ",
                      FetchElement(gf9, 5));

  std::cout << "\nFactory examples completed successfully!\n";
  return 0;
}
