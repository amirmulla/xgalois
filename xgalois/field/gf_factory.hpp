

#ifndef XGALOIS_FIELD_GF_FACTORY_HPP_
#define XGALOIS_FIELD_GF_FACTORY_HPP_

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

using GaloisFieldVariant =
    std::variant<std::shared_ptr<GaloisFieldBinary>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint16_t>>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint32_t>>,
                 std::shared_ptr<GFBELogTables<uint8_t>>,
                 std::shared_ptr<GFBELogTables<uint16_t>>,
                 std::shared_ptr<GFBELogTables<uint32_t>>,
                 std::shared_ptr<GFBELogTablesOpt<uint8_t>>,
                 std::shared_ptr<GFBELogTablesOpt<uint16_t>>,
                 std::shared_ptr<GFBELogTablesOpt<uint32_t>>,
                 std::shared_ptr<GFBEZechLogTables<uint32_t>>,
                 std::shared_ptr<GaloisFieldPrime<uint8_t>>,
                 std::shared_ptr<GaloisFieldPrime<uint16_t>>,
                 std::shared_ptr<GaloisFieldPrime<uint32_t>>,
                 std::shared_ptr<GaloisFieldExtension<uint8_t>>,
                 std::shared_ptr<GaloisFieldExtension<uint16_t>>,
                 std::shared_ptr<GaloisFieldExtension<uint32_t>>>;
;

using GaloisFieldElementVariant =
    std::variant<GaloisFieldElementBase<GaloisFieldBinary>,
                 GaloisFieldElementBase<GaloisFieldBinaryExtension<uint8_t>>,
                 GaloisFieldElementBase<GaloisFieldBinaryExtension<uint16_t>>,
                 GaloisFieldElementBase<GaloisFieldBinaryExtension<uint32_t>>,
                 GaloisFieldElementBase<GFBELogTables<uint8_t>>,
                 GaloisFieldElementBase<GFBELogTables<uint16_t>>,
                 GaloisFieldElementBase<GFBELogTables<uint32_t>>,
                 GaloisFieldElementBase<GFBELogTablesOpt<uint8_t>>,
                 GaloisFieldElementBase<GFBELogTablesOpt<uint16_t>>,
                 GaloisFieldElementBase<GFBELogTablesOpt<uint32_t>>,
                 GaloisFieldElementBase<GFBEZechLogTables<uint32_t>>,
                 GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>,
                 GaloisFieldElementBase<GaloisFieldPrime<uint16_t>>,
                 GaloisFieldElementBase<GaloisFieldPrime<uint32_t>>>;

GaloisFieldElementVariant FetchElement(const GaloisFieldVariant &field_variant,
                                       uint64_t value) {
  return std::visit(
      [value](const auto &field_ptr) -> GaloisFieldElementVariant {
        if (!field_ptr) {
          throw std::invalid_argument("Field pointer cannot be null");
        }

        using FieldType = typename std::decay_t<decltype(*field_ptr)>;

        if constexpr (std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint8_t>> ||
                      std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint16_t>> ||
                      std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint32_t>>) {
          throw std::invalid_argument(
              "Extension fields require polynomial element creation - use "
              "field methods directly");
        } else {
          using ElementType = typename FieldType::element_type;
          ElementType element_value = static_cast<ElementType>(value);
          return GaloisFieldElementBase<FieldType>(element_value, field_ptr);
        }
      },
      field_variant);
}

using GaloisFieldElementStrictVariant =
    std::variant<GaloisFieldElement<GaloisFieldBinary>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint16_t>>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint32_t>>,
                 GaloisFieldElement<GFBELogTables<uint8_t>>,
                 GaloisFieldElement<GFBELogTables<uint16_t>>,
                 GaloisFieldElement<GFBELogTables<uint32_t>>,
                 GaloisFieldElement<GFBELogTablesOpt<uint8_t>>,
                 GaloisFieldElement<GFBELogTablesOpt<uint16_t>>,
                 GaloisFieldElement<GFBELogTablesOpt<uint32_t>>,
                 GaloisFieldElement<GFBEZechLogTables<uint32_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint8_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint16_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint32_t>>>;

GaloisFieldElementStrictVariant FetchElementStrict(
    const GaloisFieldVariant &field_variant, uint64_t value) {
  return std::visit(
      [value](const auto &field_ptr) -> GaloisFieldElementStrictVariant {
        if (!field_ptr) {
          throw std::invalid_argument("Field pointer cannot be null");
        }

        using FieldType = typename std::decay_t<decltype(*field_ptr)>;

        if constexpr (std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint8_t>> ||
                      std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint16_t>> ||
                      std::is_same_v<FieldType,
                                     GaloisFieldExtension<uint32_t>>) {
          throw std::invalid_argument(
              "Extension fields require polynomial element creation - use "
              "field methods directly");
        } else {
          using ElementType = typename FieldType::element_type;
          ElementType element_value = static_cast<ElementType>(value);
          return GaloisFieldElement<FieldType>(element_value, field_ptr);
        }
      },
      field_variant);
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template <typename FieldType>
GaloisFieldElement<FieldType> FetchElement(
    const GaloisFieldVariant &field_variant, uint64_t value) {
  auto field_ptr = std::get<std::shared_ptr<FieldType>>(field_variant);
  if (!field_ptr) {
    throw std::invalid_argument("Field pointer cannot be null");
  }

  using ElementType = typename FieldType::element_type;
  ElementType element_value = static_cast<ElementType>(value);
  return GaloisFieldElement<FieldType>(element_value, field_ptr);
}
#endif

class GaloisFieldFactory {
 public:
  static GaloisFieldVariant Create(std::pair<uint64_t, uint64_t> prime_exp,
                                   const std::string &representation = "int",
                                   const std::string &modulus = "",
                                   const std::string &variable_name = "α",
                                   const std::string &impl = "auto",
                                   bool check_irreducible = false,
                                   bool prime_testing = false) {
    uint64_t prime = prime_exp.first;
    uint64_t exponent = prime_exp.second;

    if (prime < 2) {
      throw std::invalid_argument("Prime must be at least 2");
    }
    if (exponent < 1) {
      throw std::invalid_argument("Exponent must be at least 1");
    }

    if (prime == 2 && exponent == 1) {
      return CreateBinaryField(representation, impl);
    }

    if (prime == 2 && exponent > 1) {
      return CreateBinaryExtensionField(exponent, representation, modulus,
                                        variable_name, impl, check_irreducible);
    }

    if (exponent == 1) {
      return CreatePrimeField(prime, representation, prime_testing);
    }

    return CreatePrimeExtensionField(prime, exponent, modulus, variable_name,
                                     representation, impl, check_irreducible,
                                     prime_testing);
  }

  static GaloisFieldVariant Create(uint64_t order,
                                   const std::string &representation = "int",
                                   const std::string &modulus = "",
                                   const std::string &variable_name = "α",
                                   const std::string &impl = "auto",
                                   bool check_irreducible = false,
                                   bool prime_testing = false) {
    if (order < 2) {
      throw std::invalid_argument("Order of finite field must be at least 2");
    }

    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument(
          "Field order exceeds maximum supported size (uint32_t)");
    }

    auto [prime, exponent] = utils::DecomposePrimePower(order);
    if (prime == 0) {
      throw std::invalid_argument("Order must be a prime power");
    }

    return Create(std::make_pair(prime, exponent), representation, modulus,
                  variable_name, impl, check_irreducible, prime_testing);
  }

 private:
  static GaloisFieldVariant CreatePrimeField(
      uint64_t prime, const std::string &representation = "int",
      bool prime_testing = false) {
    if (prime <= std::numeric_limits<uint8_t>::max()) {
      return std::make_shared<GaloisFieldPrime<uint8_t>>(
          static_cast<uint8_t>(prime), representation, prime_testing);
    } else if (prime <= std::numeric_limits<uint16_t>::max()) {
      return std::make_shared<GaloisFieldPrime<uint16_t>>(
          static_cast<uint16_t>(prime), representation, prime_testing);
    } else if (prime <= std::numeric_limits<uint32_t>::max()) {
      return std::make_shared<GaloisFieldPrime<uint32_t>>(
          static_cast<uint32_t>(prime), representation, prime_testing);
    } else {
      throw std::invalid_argument(
          "Prime exceeds maximum supported size (uint32_t)");
    }
  }

  static GaloisFieldVariant CreatePrimeExtensionField(
      uint64_t prime, uint64_t exponent, const std::string &modulus = "",
      const std::string &variable_name = "α",
      const std::string &representation = "poly",
      const std::string &impl = "auto", bool check_irreducible = false,
      bool prime_testing = false) {
    if (exponent < 2) {
      throw std::invalid_argument("Prime extension degree must be at least 2");
    }

    uint64_t order = static_cast<uint64_t>(std::pow(prime, exponent));

    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument(
          "Prime extension field order exceeds maximum "
          "supported size (uint32_t)");
    }

    if (order <= std::numeric_limits<uint8_t>::max()) {
      return CreatePrimeExtensionFieldType<uint8_t>(
          prime, exponent, modulus, variable_name, representation, impl,
          check_irreducible, prime_testing);
    } else if (order <= std::numeric_limits<uint16_t>::max()) {
      return CreatePrimeExtensionFieldType<uint16_t>(
          prime, exponent, modulus, variable_name, representation, impl,
          check_irreducible, prime_testing);
    } else {
      return CreatePrimeExtensionFieldType<uint32_t>(
          prime, exponent, modulus, variable_name, representation, impl,
          check_irreducible, prime_testing);
    }
  }

  template <typename ElementType>
  static GaloisFieldVariant CreatePrimeExtensionFieldType(
      uint64_t prime, uint64_t exponent, const std::string &modulus,
      const std::string &variable_name, const std::string &representation,
      const std::string &impl, bool check_irreducible, bool prime_testing) {
    if (exponent < 2) {
      throw std::invalid_argument("Extension degree must be at least 2");
    }

    return std::make_shared<GaloisFieldExtension<ElementType>>(
        std::make_pair(prime, exponent), modulus, representation, variable_name,
        check_irreducible, prime_testing);
  }

  static GaloisFieldVariant CreateBinaryField(
      const std::string &representation = "int",
      const std::string &impl = "auto") {
    auto field = std::make_shared<GaloisFieldBinary>();
    field->SetRepresentation(utils::ConvertRepresentation(representation));
    return field;
  }

  static GaloisFieldVariant CreateBinaryExtensionField(
      uint64_t exponent, const std::string &representation = "int",
      const std::string &modulus = "", const std::string &variable_name = "α",
      const std::string &impl = "auto", bool check_irreducible = false) {
    if (exponent < 2) {
      throw std::invalid_argument("Binary extension degree must be at least 2");
    }

    uint64_t order = 1ULL << exponent;

    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument(
          "Binary extension field order exceeds "
          "maximum supported size (uint32_t)");
    }

    if (impl == "auto") {
      if (exponent <= 7) {
        return CreateBinaryExtensionFieldTyped<uint8_t>(
            exponent, representation, modulus, variable_name, "log_opt",
            check_irreducible);
      } else if (exponent <= 15) {
        return CreateBinaryExtensionFieldTyped<uint16_t>(
            exponent, representation, modulus, variable_name, "log",
            check_irreducible);
      } else {
        return CreateBinaryExtensionFieldTyped<uint32_t>(
            exponent, representation, modulus, variable_name, "standard",
            check_irreducible);
      }
    }

    if (order <= std::numeric_limits<uint8_t>::max()) {
      return CreateBinaryExtensionFieldTyped<uint8_t>(exponent, representation,
                                                      modulus, variable_name,
                                                      impl, check_irreducible);
    } else if (order <= std::numeric_limits<uint16_t>::max()) {
      return CreateBinaryExtensionFieldTyped<uint16_t>(exponent, representation,
                                                       modulus, variable_name,
                                                       impl, check_irreducible);
    } else {
      return CreateBinaryExtensionFieldTyped<uint32_t>(exponent, representation,
                                                       modulus, variable_name,
                                                       impl, check_irreducible);
    }
  }

  template <typename ElementType>
  static GaloisFieldVariant CreateBinaryExtensionFieldTyped(
      uint64_t exponent, const std::string &representation,
      const std::string &modulus, const std::string &variable_name,
      const std::string &impl, bool check_irreducible) {
    uint8_t m = static_cast<uint8_t>(exponent);

    if (impl == "standard" || impl == "poly") {
      auto field = std::make_shared<GaloisFieldBinaryExtension<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
      return field;

    } else if (impl == "log" || impl == "log_tables") {
      auto field = std::make_shared<GFBELogTables<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
      return field;

    } else if (impl == "log_opt" || impl == "log_optimized") {
      auto field = std::make_shared<GFBELogTablesOpt<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
      return field;

    } else if (impl == "zech" || impl == "zech_log") {
      if constexpr (std::is_same_v<ElementType, uint32_t>) {
        auto field = std::make_shared<GFBEZechLogTables<uint32_t>>(
            m, representation, modulus, variable_name, check_irreducible);
        return field;
      } else {
        throw std::invalid_argument(
            "Zech logarithm implementation only "
            "supports uint32_t element type");
      }

    } else {
      throw std::invalid_argument(
          "Unknown binary field implementation: " + impl +
          ". Available options: standard, log, log_opt, zech");
    }
  }
};

GaloisFieldVariant GF(std::pair<uint64_t, uint64_t> prime_exp,
                      const std::string &representation = "int",
                      const std::string &modulus = "",
                      const std::string &variable_name = "",
                      const std::string &impl = "auto",
                      bool check_irreducible = false,
                      bool prime_testing = false) {
  return GaloisFieldFactory::Create(prime_exp, representation, modulus,
                                    variable_name, impl, check_irreducible,
                                    prime_testing);
}

GaloisFieldVariant GF(uint64_t order, const std::string &representation = "int",
                      const std::string &modulus = "",
                      const std::string &variable_name = "",
                      const std::string &impl = "auto",
                      bool check_irreducible = false,
                      bool prime_testing = false) {
  return GaloisFieldFactory::Create(order, representation, modulus,
                                    variable_name, impl, check_irreducible,
                                    prime_testing);
}

}  // namespace xg

#endif
