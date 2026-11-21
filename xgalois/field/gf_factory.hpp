#ifndef XGALOIS_FIELD_GF_FACTORY_HPP_
#define XGALOIS_FIELD_GF_FACTORY_HPP_

#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

namespace detail {

template <typename T>
inline constexpr bool kIsGaloisFieldExtension = false;

template <UnsignedFieldWord T>
inline constexpr bool kIsGaloisFieldExtension<GaloisFieldExtension<T>> = true;

}  // namespace detail

using GaloisFieldVariant =
    std::variant<std::shared_ptr<GaloisFieldBinary>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint16_t>>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint32_t>>,
                 std::shared_ptr<GaloisFieldBinaryExtension<uint64_t>>,
                 std::shared_ptr<GFBELogTables<uint8_t>>,
                 std::shared_ptr<GFBELogTables<uint16_t>>,
                 std::shared_ptr<GFBELogTables<uint32_t>>,
                 std::shared_ptr<GFBELogTablesOpt<uint8_t>>,
                 std::shared_ptr<GFBELogTablesOpt<uint16_t>>,
                 std::shared_ptr<GFBELogTablesOpt<uint32_t>>,
                 std::shared_ptr<GFBEZechLogTables<uint8_t>>,
                 std::shared_ptr<GFBEZechLogTables<uint16_t>>,
                 std::shared_ptr<GFBEZechLogTables<uint32_t>>,
                 std::shared_ptr<GaloisFieldPrime<uint8_t>>,
                 std::shared_ptr<GaloisFieldPrime<uint16_t>>,
                 std::shared_ptr<GaloisFieldPrime<uint32_t>>,
                 std::shared_ptr<GaloisFieldPrimeTable<uint8_t>>,
                 std::shared_ptr<GaloisFieldPrimeTable<uint16_t>>,
                 std::shared_ptr<GaloisFieldPrimeTable<uint32_t>>,
                 std::shared_ptr<GaloisFieldExtension<uint8_t>>,
                 std::shared_ptr<GaloisFieldExtension<uint16_t>>,
                 std::shared_ptr<GaloisFieldExtension<uint32_t>>>;

using GaloisFieldElementVariant =
    std::variant<GaloisFieldElement<GaloisFieldBinary>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint16_t>>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint32_t>>,
                 GaloisFieldElement<GaloisFieldBinaryExtension<uint64_t>>,
                 GaloisFieldElement<GFBELogTables<uint8_t>>,
                 GaloisFieldElement<GFBELogTables<uint16_t>>,
                 GaloisFieldElement<GFBELogTables<uint32_t>>,
                 GaloisFieldElement<GFBELogTablesOpt<uint8_t>>,
                 GaloisFieldElement<GFBELogTablesOpt<uint16_t>>,
                 GaloisFieldElement<GFBELogTablesOpt<uint32_t>>,
                 GaloisFieldElement<GFBEZechLogTables<uint8_t>>,
                 GaloisFieldElement<GFBEZechLogTables<uint16_t>>,
                 GaloisFieldElement<GFBEZechLogTables<uint32_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint8_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint16_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint32_t>>,
                 GaloisFieldElement<GaloisFieldPrimeTable<uint8_t>>,
                 GaloisFieldElement<GaloisFieldPrimeTable<uint16_t>>,
                 GaloisFieldElement<GaloisFieldPrimeTable<uint32_t>>,
                 GaloisFieldElement<GaloisFieldExtension<uint8_t>>,
                 GaloisFieldElement<GaloisFieldExtension<uint16_t>>,
                 GaloisFieldElement<GaloisFieldExtension<uint32_t>>>;

using GaloisFieldElementStrictVariant = GaloisFieldElementVariant;

inline GaloisFieldElementVariant FetchElement(
    const GaloisFieldVariant &field_variant, uint64_t value) {
  return std::visit(
      [value](const auto &field_ptr) -> GaloisFieldElementVariant {
        if (!field_ptr) {
          throw std::invalid_argument("Field pointer cannot be null");
        }

        using FieldType = std::decay_t<decltype(*field_ptr)>;

        if constexpr (detail::kIsGaloisFieldExtension<FieldType>) {
          return GaloisFieldElement<FieldType>(field_ptr->FromInt(value),
                                               field_ptr);
        } else {
          using ElementType = typename FieldType::element_type;
          const uint64_t order = field_ptr->Order();
          const uint64_t reduced = (order == 0) ? value : (value % order);
          return GaloisFieldElement<FieldType>(
              static_cast<ElementType>(reduced), field_ptr);
        }
      },
      field_variant);
}

inline GaloisFieldElementStrictVariant FetchElementStrict(
    const GaloisFieldVariant &field_variant, uint64_t value) {
  return FetchElement(field_variant, value);
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template <typename FieldType>
GaloisFieldElement<FieldType> FetchElement(
    const GaloisFieldVariant &field_variant, uint64_t value) {
  auto field_ptr = std::get<std::shared_ptr<FieldType>>(field_variant);
  if (!field_ptr) {
    throw std::invalid_argument("Field pointer cannot be null");
  }

  if constexpr (detail::kIsGaloisFieldExtension<FieldType>) {
    return GaloisFieldElement<FieldType>(field_ptr->FromInt(value), field_ptr);
  } else {
    using ElementType = typename FieldType::element_type;
    const uint64_t order = field_ptr->Order();
    const uint64_t reduced = (order == 0) ? value : (value % order);
    return GaloisFieldElement<FieldType>(static_cast<ElementType>(reduced),
                                         field_ptr);
  }
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
                                   bool prime_testing = true) {
    const uint64_t prime = prime_exp.first;
    const uint64_t exponent = prime_exp.second;

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
      return CreatePrimeField(prime, representation, impl, prime_testing);
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
                                   bool prime_testing = true) {
    if (order < 2) {
      throw std::invalid_argument("Order of finite field must be at least 2");
    }

    auto [prime, exponent] = utils::DecomposePrimePower(order);
    if (prime == 0) {
      throw std::invalid_argument("Order must be a prime power");
    }

    return Create(std::make_pair(prime, exponent), representation, modulus,
                  variable_name, impl, check_irreducible, prime_testing);
  }

 private:
  template <UnsignedFieldWord ElementType>
  static GaloisFieldVariant CreatePrimeFieldTyped(
      uint64_t prime, const std::string &representation,
      const std::string &impl, bool prime_testing) {
    const ElementType p = static_cast<ElementType>(prime);
    if (impl == "auto" || impl == "standard") {
      return std::make_shared<GaloisFieldPrime<ElementType>>(p, representation,
                                                             prime_testing);
    }
    if (impl == "table" || impl == "log") {
      return std::make_shared<GaloisFieldPrimeTable<ElementType>>(
          p, representation, prime_testing);
    }
    throw std::invalid_argument(
        "Unknown prime field implementation: " + impl +
        ". Available options: auto, standard, table, log");
  }

  static GaloisFieldVariant CreatePrimeField(
      uint64_t prime, const std::string &representation = "int",
      const std::string &impl = "auto", bool prime_testing = true) {
    if (prime <= std::numeric_limits<uint8_t>::max()) {
      return CreatePrimeFieldTyped<uint8_t>(prime, representation, impl,
                                            prime_testing);
    }
    if (prime <= std::numeric_limits<uint16_t>::max()) {
      return CreatePrimeFieldTyped<uint16_t>(prime, representation, impl,
                                             prime_testing);
    }
    if (prime <= std::numeric_limits<uint32_t>::max()) {
      return CreatePrimeFieldTyped<uint32_t>(prime, representation, impl,
                                             prime_testing);
    }
    throw std::invalid_argument(
        "Prime exceeds maximum supported size (uint32_t)");
  }

  static GaloisFieldVariant CreatePrimeExtensionField(
      uint64_t prime, uint64_t exponent, const std::string &modulus = "",
      const std::string &variable_name = "α",
      const std::string &representation = "poly",
      const std::string &impl = "auto", bool check_irreducible = false,
      bool prime_testing = true) {
    if (exponent < 2) {
      throw std::invalid_argument("Prime extension degree must be at least 2");
    }

    if (impl != "auto" && impl != "standard" && impl != "poly") {
      throw std::invalid_argument(
          "Unknown prime extension field implementation: " + impl +
          ". Available options: auto, standard, poly");
    }

    try {
      (void)utils::SafeIntegerPower(prime, exponent);
    } catch (const std::overflow_error &) {
      throw std::invalid_argument(
          "Prime extension field order exceeds maximum supported size "
          "(uint64_t)");
    }

    if (prime <= std::numeric_limits<uint8_t>::max()) {
      return std::make_shared<GaloisFieldExtension<uint8_t>>(
          std::make_pair(prime, exponent), modulus, representation,
          variable_name, check_irreducible, prime_testing);
    }
    if (prime <= std::numeric_limits<uint16_t>::max()) {
      return std::make_shared<GaloisFieldExtension<uint16_t>>(
          std::make_pair(prime, exponent), modulus, representation,
          variable_name, check_irreducible, prime_testing);
    }
    if (prime <= std::numeric_limits<uint32_t>::max()) {
      return std::make_shared<GaloisFieldExtension<uint32_t>>(
          std::make_pair(prime, exponent), modulus, representation,
          variable_name, check_irreducible, prime_testing);
    }
    throw std::invalid_argument(
        "Prime exceeds maximum supported size (uint32_t)");
  }

  static GaloisFieldVariant CreateBinaryField(
      const std::string &representation = "int",
      const std::string &impl = "auto") {
    if (impl != "auto" && impl != "standard" && impl != "binary") {
      throw std::invalid_argument("Unknown binary field implementation: " +
                                  impl + ". Available options: auto, standard");
    }
    return std::make_shared<GaloisFieldBinary>(representation);
  }

  static GaloisFieldVariant CreateBinaryExtensionField(
      uint64_t exponent, const std::string &representation = "int",
      const std::string &modulus = "", const std::string &variable_name = "α",
      const std::string &impl = "auto", bool check_irreducible = false) {
    if (exponent < 2 || exponent > 64) {
      throw std::invalid_argument(
          "Binary extension degree must be between 2 and 64");
    }

    if (impl == "auto") {
      if (exponent <= 8) {
        return CreateBinaryExtensionFieldTyped<uint8_t>(
            exponent, representation, modulus, variable_name, "log_opt",
            check_irreducible);
      }
      if (exponent <= 16) {
        return CreateBinaryExtensionFieldTyped<uint16_t>(
            exponent, representation, modulus, variable_name, "log",
            check_irreducible);
      }
      if (exponent <= 32) {
        return CreateBinaryExtensionFieldTyped<uint32_t>(
            exponent, representation, modulus, variable_name, "standard",
            check_irreducible);
      }
      return CreateBinaryExtensionFieldTyped<uint64_t>(
          exponent, representation, modulus, variable_name, "standard",
          check_irreducible);
    }

    if (exponent <= 8) {
      return CreateBinaryExtensionFieldTyped<uint8_t>(exponent, representation,
                                                      modulus, variable_name,
                                                      impl, check_irreducible);
    }
    if (exponent <= 16) {
      return CreateBinaryExtensionFieldTyped<uint16_t>(exponent, representation,
                                                       modulus, variable_name,
                                                       impl, check_irreducible);
    }
    if (exponent <= 32) {
      return CreateBinaryExtensionFieldTyped<uint32_t>(exponent, representation,
                                                       modulus, variable_name,
                                                       impl, check_irreducible);
    }
    return CreateBinaryExtensionFieldTyped<uint64_t>(exponent, representation,
                                                     modulus, variable_name,
                                                     impl, check_irreducible);
  }

  template <UnsignedFieldWord ElementType>
  static GaloisFieldVariant CreateBinaryExtensionFieldTyped(
      uint64_t exponent, const std::string &representation,
      const std::string &modulus, const std::string &variable_name,
      const std::string &impl, bool check_irreducible) {
    const uint8_t m = static_cast<uint8_t>(exponent);

    if (impl == "standard" || impl == "poly") {
      return std::make_shared<GaloisFieldBinaryExtension<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
    }

    if (impl == "log" || impl == "log_tables") {
      if constexpr (!std::is_same_v<ElementType, uint64_t>) {
        return std::make_shared<GFBELogTables<ElementType>>(
            m, representation, modulus, variable_name, check_irreducible);
      } else {
        throw std::invalid_argument(
            "Log table implementation is not supported for degree > 32");
      }
    }

    if (impl == "log_opt" || impl == "log_optimized") {
      if constexpr (!std::is_same_v<ElementType, uint64_t>) {
        return std::make_shared<GFBELogTablesOpt<ElementType>>(
            m, representation, modulus, variable_name, check_irreducible);
      } else {
        throw std::invalid_argument(
            "Optimized log table implementation is not supported for "
            "degree > 32");
      }
    }

    if (impl == "zech" || impl == "zech_log") {
      if constexpr (!std::is_same_v<ElementType, uint64_t>) {
        if (exponent > 20) {
          throw std::invalid_argument(
              "Zech logarithm implementation supports degree at most 20");
        }
        return std::make_shared<GFBEZechLogTables<ElementType>>(
            m, representation, modulus, variable_name, check_irreducible);
      } else {
        throw std::invalid_argument(
            "Zech logarithm implementation supports degree at most 20");
      }
    }

    throw std::invalid_argument(
        "Unknown binary field implementation: " + impl +
        ". Available options: standard, log, log_opt, zech");
  }
};

inline GaloisFieldVariant GF(std::pair<uint64_t, uint64_t> prime_exp,
                             const std::string &representation = "int",
                             const std::string &modulus = "",
                             const std::string &variable_name = "α",
                             const std::string &impl = "auto",
                             bool check_irreducible = false,
                             bool prime_testing = true) {
  return GaloisFieldFactory::Create(prime_exp, representation, modulus,
                                    variable_name, impl, check_irreducible,
                                    prime_testing);
}

inline GaloisFieldVariant GF(uint64_t order,
                             const std::string &representation = "int",
                             const std::string &modulus = "",
                             const std::string &variable_name = "α",
                             const std::string &impl = "auto",
                             bool check_irreducible = false,
                             bool prime_testing = true) {
  return GaloisFieldFactory::Create(order, representation, modulus,
                                    variable_name, impl, check_irreducible,
                                    prime_testing);
}

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_FACTORY_HPP_
