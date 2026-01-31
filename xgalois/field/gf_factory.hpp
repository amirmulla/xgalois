/*
 * Galois Field Factory Implementation
 *
 * This header provides a comprehensive factory system for creating and managing
 * various types of finite fields (Galois fields). The implementation supports
 * prime fields GF(p), binary fields GF(2), binary extension fields GF(2^m),
 * and general extension fields GF(p^n) with multiple optimized implementations.
 */

#ifndef XGALOIS_FIELD_GF_FACTORY_HPP_
#define XGALOIS_FIELD_GF_FACTORY_HPP_

/* Standard C library headers for mathematical operations and type definitions
 */
#include <cmath>
#include <cstdint>

/* C++ standard library headers for containers, memory management, and utilities
 */
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

/* Project-specific headers providing field implementations and utilities */
#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {

/**
 * @brief Type-safe variant container for all supported Galois field
 * implementations
 *
 * This std::variant encapsulates all possible Galois field types that can be
 * created by the factory system. It provides runtime polymorphism while
 * maintaining type safety and avoiding virtual function overhead. The variant
 * includes:
 *
 * - GaloisFieldBinary: The base binary field GF(2)
 * - GaloisFieldBinaryExtension<T>: Standard binary extension fields GF(2^m)
 * - GFBELogTables<T>: Binary extensions with logarithm table optimizations
 * - GFBELogTablesOpt<T>: Binary extensions with enhanced logarithm tables
 * - GFBEZechLogTables: Binary extensions using Zech logarithm tables
 * - GaloisFieldPrime<T>: Prime fields GF(p) for various element types
 * - GaloisFieldExtension<T>: General extension fields GF(p^n)
 *
 * Each template parameter T represents the underlying integer type (uint8_t,
 * uint16_t, uint32_t) used to store field elements, automatically
 * selected based on the field order for optimal memory usage and performance.
 */
using GaloisFieldVariant = std::variant<
    std::shared_ptr<GaloisFieldBinary>,
    std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>>,
    std::shared_ptr<GaloisFieldBinaryExtension<uint16_t>>,
    std::shared_ptr<GaloisFieldBinaryExtension<uint32_t>>,
    std::shared_ptr<GFBELogTables<uint8_t>>,
    std::shared_ptr<GFBELogTables<uint16_t>>,
    std::shared_ptr<GFBELogTables<uint32_t>>,
    std::shared_ptr<GFBELogTablesOpt<uint8_t>>,
    std::shared_ptr<GFBELogTablesOpt<uint16_t>>,
    std::shared_ptr<GFBELogTablesOpt<uint32_t>>,
    std::shared_ptr<GFBEZechLogTables>,
    std::shared_ptr<GaloisFieldPrime<uint8_t>>,
    std::shared_ptr<GaloisFieldPrime<uint16_t>>,
    std::shared_ptr<GaloisFieldPrime<uint32_t>>,
    std::shared_ptr<GaloisFieldExtension<uint8_t>>,
    std::shared_ptr<GaloisFieldExtension<uint16_t>>,
    std::shared_ptr<GaloisFieldExtension<uint32_t>>>;;

/**
 * @brief Type-safe variant container for Galois field elements with simple
 * representations
 *
 * This variant supports field elements that can be represented as simple
 * integer values. Extension fields are not included as they require polynomial
 * representations.
 */
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
                 GaloisFieldElementBase<GFBEZechLogTables>,
                 GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>,
                 GaloisFieldElementBase<GaloisFieldPrime<uint16_t>>,
                 GaloisFieldElementBase<GaloisFieldPrime<uint32_t>>>;

/**
 * @brief Factory function to create field elements from integer values
 *
 * Creates field elements for fields that can represent elements as simple
 * integers. This excludes extension fields which require polynomial
 * representations.
 *
 * @param field_variant The GaloisFieldVariant containing the field
 * implementation
 * @param value The integer value for the element
 *
 * @return GaloisFieldElementVariant containing the created element with proper
 * type
 *
 * @throws std::invalid_argument if field_variant contains null pointer or
 * unsupported field type
 *
 * @example
 *   auto gf8 = GF(8);     // GF(2^3) - binary extension field
 *   auto elem = FetchElement(gf8, 5);  // Creates element with value 5 in GF(8)
 *
 *   auto gf7 = GF(7);     // GF(7) - prime field
 *   auto elem2 = FetchElement(gf7, 3); // Creates element with value 3 in GF(7)
 */
GaloisFieldElementVariant FetchElement(const GaloisFieldVariant &field_variant,
                                       uint64_t value) {
  return std::visit(
      [value](const auto &field_ptr) -> GaloisFieldElementVariant {
        if (!field_ptr) {
          throw std::invalid_argument("Field pointer cannot be null");
        }

        using FieldType = typename std::decay_t<decltype(*field_ptr)>;

        // Check if this is an extension field (which we don't support for
        // simple element creation)
        if constexpr (
            std::is_same_v<FieldType,
                           GaloisFieldExtension<uint8_t>> ||
            std::is_same_v<FieldType,
                           GaloisFieldExtension<uint16_t>> ||
            std::is_same_v<FieldType,
                           GaloisFieldExtension<uint32_t>>) {
          throw std::invalid_argument(
              "Extension fields require polynomial element creation - use "
              "field methods directly");
        } else {
          // For supported field types, create element directly
          using ElementType = typename FieldType::element_type;
          ElementType element_value = static_cast<ElementType>(value);
          return GaloisFieldElementBase<FieldType>(element_value, field_ptr);
        }
      },
      field_variant);
}

/**
 * @brief Type-safe variant container for Galois field elements with strict
 * validation
 *
 * This variant supports field elements with enhanced validation that ensures
 * arithmetic operations only occur between elements of the same field.
 * Extension fields are not included as they require polynomial representations.
 */
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
                 GaloisFieldElement<GFBEZechLogTables>,
                 GaloisFieldElement<GaloisFieldPrime<uint8_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint16_t>>,
                 GaloisFieldElement<GaloisFieldPrime<uint32_t>>>;

/**
 * @brief Factory function to create field elements with strict validation from
 * integer values
 *
 * Creates field elements for fields that can represent elements as simple
 * integers. This excludes extension fields which require polynomial
 * representations. Returns GaloisFieldElement types which provide enhanced
 * validation ensuring arithmetic operations only occur between elements of the
 * same field.
 *
 * @param field_variant The GaloisFieldVariant containing the field
 * implementation
 * @param value The integer value for the element
 *
 * @return GaloisFieldElementStrictVariant containing the created element with
 * proper type and validation
 *
 * @throws std::invalid_argument if field_variant contains null pointer or
 * unsupported field type
 *
 * @example
 *   auto gf8 = GF(8);     // GF(2^3) - binary extension field
 *   auto elem = FetchElementStrict(gf8, 5);  // Creates element with value 5 in
 * GF(8) with strict validation
 *
 *   auto gf7 = GF(7);     // GF(7) - prime field
 *   auto elem2 = FetchElementStrict(gf7, 3); // Creates element with value 3 in
 * GF(7) with strict validation
 */
GaloisFieldElementStrictVariant
FetchElementStrict(const GaloisFieldVariant &field_variant, uint64_t value) {
  return std::visit(
      [value](const auto &field_ptr) -> GaloisFieldElementStrictVariant {
        if (!field_ptr) {
          throw std::invalid_argument("Field pointer cannot be null");
        }

        using FieldType = typename std::decay_t<decltype(*field_ptr)>;

        // Check if this is an extension field (which we don't support for
        // simple element creation)
        if constexpr (
            std::is_same_v<FieldType,
                           GaloisFieldExtension<uint8_t>> ||
            std::is_same_v<FieldType,
                           GaloisFieldExtension<uint16_t>> ||
            std::is_same_v<FieldType,
                           GaloisFieldExtension<uint32_t>>) {
          throw std::invalid_argument(
              "Extension fields require polynomial element creation - use "
              "field methods directly");
        } else {
          // For supported field types, create element directly
          using ElementType = typename FieldType::element_type;
          ElementType element_value = static_cast<ElementType>(value);
          return GaloisFieldElement<FieldType>(element_value, field_ptr);
        }
      },
      field_variant);
}

/**
 * @brief Template factory function to create a specific GaloisFieldElement type
 * from integer values
 *
 * Creates a field element of a specific GaloisFieldElement<FieldType> from a
 * field variant. This function requires explicit template specification and
 * validates that the variant contains the expected field type.
 *
 * @tparam FieldType The specific field type (e.g., GaloisFieldPrime<uint32_t>)
 * @param field_variant The GaloisFieldVariant containing the field
 * implementation
 * @param value The integer value for the element
 *
 * @return GaloisFieldElement<FieldType> containing the created element
 *
 * @throws std::invalid_argument if field_variant doesn't contain the expected
 * field type
 * @throws std::bad_variant_access if the variant doesn't hold the expected type
 *
 * @example
 *   auto gf7 = GF(7);     // GF(7) - prime field
 *   auto elem = FetchElement<GaloisFieldPrime<uint32_t>>(gf7, 3); // Creates
 * element with value 3
 */
template <typename FieldType>
GaloisFieldElement<FieldType>
FetchElement(const GaloisFieldVariant &field_variant, uint64_t value) {
  auto field_ptr = std::get<std::shared_ptr<FieldType>>(field_variant);
  if (!field_ptr) {
    throw std::invalid_argument("Field pointer cannot be null");
  }

  using ElementType = typename FieldType::element_type;
  ElementType element_value = static_cast<ElementType>(value);
  return GaloisFieldElement<FieldType>(element_value, field_ptr);
}

/**
 * @brief Advanced factory class for creating optimized finite field
 * implementations
 *
 * The GaloisFieldFactory provides a sophisticated interface
 * for constructing finite fields with automatic optimization selection. The
 * factory intelligently chooses the most efficient implementation based on
 * field characteristics such as order, prime base, and extension degree.
 *
 * Key capabilities:
 * - Automatic element type selection based on field order for memory efficiency
 *   (supports up to uint32_t maximum field size)
 * - Multiple implementation strategies for binary extension fields (standard
 *   polynomial arithmetic, logarithm tables, optimized tables, Zech logarithms)
 * - Prime power decomposition using advanced factorization algorithms
 * - String-based configuration for representations and implementation
 * preferences
 * - Comprehensive error checking with meaningful diagnostic messages
 *
 * The factory supports the mathematical notation GF(p^n) where p is prime and
 * n is the extension degree, automatically handling special cases like GF(2)
 * and binary extension fields GF(2^m) with specialized optimizations.
 */
class GaloisFieldFactory {
public:
  /**
   * @brief Primary factory method for creating finite fields using (prime,
   * exponent) specification
   *
   * This is the principal interface for field creation, accepting a (prime,
   * exponent) pair to define the field GF(p^n). The method automatically
   * selects the optimal element type and implementation strategy based on the
   * resulting field order and characteristics.
   *
   * @param prime_exp Pair containing (prime_base, extension_degree) defining
   * GF(p^n)
   * @param representation Display format for field elements: "int", "hex",
   * "pow", "log", "poly"
   * @param modulus String representation of irreducible polynomial for
   * extensions (optional - uses Conway polynomial from database if available)
   * @param variable_name Symbol used for the extension field generator
   * (default: empty)
   * @param impl Implementation preference: "auto", "standard", "log",
   * "log_opt", "zech"
   * @param check_irreducible Whether to verify polynomial irreducibility
   * (disabled by default for performance)
   * @param prime_testing Whether to verify that the base prime is actually prime
   *
   * @return GaloisFieldVariant containing the optimally configured field
   * implementation
   *
   * @throws std::invalid_argument if prime < 2, exponent < 1, or no suitable
   * polynomial found
   *
   * Special handling:
   * - GF(2): Returns specialized binary field implementation
   * - GF(2^m): Uses optimized binary extension field with automatic algorithm
   * selection
   * - GF(p^n): Creates general extension field over prime base field
   */
  static GaloisFieldVariant Create(std::pair<uint64_t, uint64_t> prime_exp,
                                   const std::string &representation = "int",
                                   const std::string &modulus = "",
                                   const std::string &variable_name = "α",
                                   const std::string &impl = "auto",
                                   bool check_irreducible = false,
                                   bool prime_testing = false) {

    uint64_t prime = prime_exp.first;
    uint64_t exponent = prime_exp.second;

    /* Input validation with comprehensive error messages */
    if (prime < 2) {
      throw std::invalid_argument("Prime must be at least 2");
    }
    if (exponent < 1) {
      throw std::invalid_argument("Exponent must be at least 1");
    }

    /* Specialized optimized path for the fundamental binary field GF(2) */
    if (prime == 2 && exponent == 1) {
      return CreateBinaryField(representation, impl);
    }

    /* Specialized optimized path for binary extension fields GF(2^m) with multiple algorithms */
    if (prime == 2 && exponent > 1) {
      return CreateBinaryExtensionField(exponent, representation, modulus,
                                        variable_name, impl, check_irreducible);
    }

    /* Handle prime fields GF(p) */
    if (exponent == 1) {
      return CreatePrimeField(prime, representation, prime_testing);
    }

    /* Handle prime extension fields GF(p^n) */
    return CreatePrimeExtensionField(prime, exponent, modulus, variable_name,
                                     representation, impl, check_irreducible, prime_testing);
  }

  /**
   * @brief Factory method for creating finite fields using field order
   * specification
   *
   * This method allows creation of finite fields by directly specifying the
   * field order as a single integer. It automatically decomposes the order into
   * its prime power components and selects the appropriate field type based on
   * the order characteristics.
   *
   * @param order The total number of elements in the field (must be ≥ 2)
   * @param representation Display format for field elements: "int", "hex",
   * "pow", "log", "poly"
   * @param modulus String representation of irreducible polynomial for
   * extensions (optional - uses Conway polynomial from database if available)
   * @param variable_name Symbol used for the extension field generator
   * (default: "α")
   * @param impl Implementation preference: "auto", "standard", "log",
   * "log_opt", "zech"
   * @param check_irreducible Whether to verify polynomial irreducibility
   * @param prime_testing Whether to verify that the base prime is actually prime
   *
   * @return GaloisFieldVariant containing the constructed finite field
   *
   * @throws std::invalid_argument if order < 2, decomposition fails, or no
   * suitable polynomial found
   */
  static GaloisFieldVariant Create(uint64_t order,
                                   const std::string &representation = "int",
                                   const std::string &modulus = "",
                                   const std::string &variable_name = "α",
                                   const std::string &impl = "auto",
                                   bool check_irreducible = false,
                                   bool prime_testing = false) {

    /* Validate field order constraints */
    if (order < 2) {
      throw std::invalid_argument("Order of finite field must be at least 2");
    }

    /* Validate field order is within supported range */
    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument(
          "Field order exceeds maximum supported size (uint32_t)");
    }

    /* Decompose order into prime power components using advanced factorization
     */
    auto [prime, exponent] = utils::DecomposePrimePower(order);
    if (prime == 0) {
      throw std::invalid_argument("Order must be a prime power");
    }

    /* Delegate to primary factory method with decomposed parameters */
    return Create(std::make_pair(prime, exponent), representation, modulus,
                  variable_name, impl, check_irreducible, prime_testing);
  }

private:
  /**
   * @brief Factory method for prime fields with automatic element type
   * selection
   *
   * Creates a prime field GF(p) with the smallest suitable integer type based
   * on the prime value. This optimization reduces memory usage and can improve
   * performance for smaller primes.
   *
   * @param prime The prime number defining the field
   * @param representation Display format for field elements
   * @param prime_testing Whether to verify that the input is actually prime
   *
   * @return GaloisFieldVariant containing optimally typed prime field
   */
  static GaloisFieldVariant
  CreatePrimeField(uint64_t prime, const std::string &representation = "int",
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

  /**
   * @brief Advanced factory for prime extension fields GF(p^n) with automatic
   * type selection
   *
   * Creates prime extension fields over any prime base field with automatic
   * element type selection based on field order. This method provides the same
   * interface consistency as CreateBinaryExtensionField for prime-based
   * extensions.
   *
   * @param prime Base prime p for the extension GF(p^n)
   * @param exponent Extension degree n (must be ≥ 2)
   * @param modulus Irreducible polynomial specification (optional - uses Conway
   * polynomial from database if available)
   * @param variable_name Symbol for the extension generator
   * @param representation Element display format (typically "poly" for
   * extensions)
   * @param impl Implementation preference (may be limited for general
   * extensions)
   * @param check_irreducible Polynomial verification flag
   * @param prime_testing Whether to verify that the base prime is actually prime
   *
   * @return GaloisFieldVariant containing the constructed extension field
   *
   * @throws std::invalid_argument if exponent < 2, no suitable polynomial
   * found, or field order exceeds limits
   */
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

    /* Validate field order is within supported range */
    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument("Prime extension field order exceeds maximum "
                                  "supported size (uint32_t)");
    }

    /* Automatic element type selection based on field order for optimal memory usage */
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

  // Helper method to create prime extension fields with specific element type
  template <typename ElementType>
  static GaloisFieldVariant CreatePrimeExtensionFieldType(
      uint64_t prime, uint64_t exponent, const std::string &modulus,
      const std::string &variable_name, const std::string &representation,
      const std::string &impl, bool check_irreducible, bool prime_testing) {

    if (exponent < 2) {
      throw std::invalid_argument("Extension degree must be at least 2");
    }

    /* Construct the extension field with all validated parameters - let the constructor
     * handle database fetching, polynomial parsing, and irreducibility checking */
    return std::make_shared<GaloisFieldExtension<ElementType>>(
        std::make_pair(prime, exponent), modulus,
        representation, variable_name, check_irreducible, prime_testing);
  }

  /**
   * @brief Specialized factory for the fundamental binary field GF(2)
   *
   * Creates the base binary field with two elements {0, 1}. This field serves
   * as the foundation for all binary extension fields and has unique properties
   * in cryptographic and coding theory applications.
   *
   * @param representation Display format (typically "int" for binary)
   * @param impl Implementation preference (ignored for GF(2) as only one
   * implementation exists)
   *
   * @return GaloisFieldVariant containing the binary field GF(2)
   *
   * @note The impl parameter is accepted for interface consistency but has no
   *       effect since GF(2) has only one possible implementation
   */
  static GaloisFieldVariant
  CreateBinaryField(const std::string &representation = "int",
                    const std::string &impl = "auto") {
    /* GF(2) has only one standard implementation regardless of impl parameter
     */
    auto field = std::make_shared<GaloisFieldBinary>();
    field->SetRepresentation(utils::ConvertRepresentation(representation));
    return field;
  }

  /**
   * @brief Advanced factory for binary extension fields GF(2^m) with algorithm
   * selection
   *
   * Creates binary extension fields with multiple optimized implementation
   * strategies. The method provides sophisticated algorithm selection based on
   * field size and performance characteristics, enabling optimal performance
   * for different applications.
   *
   * Implementation algorithms available:
   * - "standard"/"poly": Classical polynomial arithmetic with carry-less
   * multiplication
   * - "log"/"log_tables": Logarithm-based multiplication using discrete log
   * tables
   * - "log_opt"/"log_optimized": Enhanced logarithm tables with extended lookup
   * optimization
   * - "zech"/"zech_log": Zech logarithm tables for accelerated addition and
   * multiplication
   * - "auto": Intelligent algorithm selection based on extension degree and
   * performance profiles
   *
   * @param exponent Extension degree m for GF(2^m) (must be ≥ 2)
   * @param representation Element display format
   * @param modulus Irreducible polynomial specification (optional - uses Conway
   * polynomial from database if available)
   * @param variable_name Symbol for the primitive element (default: α)
   * @param impl Implementation algorithm preference
   * @param check_irreducible Whether to verify polynomial irreducibility
   *
   * @return GaloisFieldVariant with optimally configured binary extension field
   *
   * @throws std::invalid_argument if exponent < 2, no suitable polynomial
   * found, or polynomial verification fails
   *
   * @note Automatic algorithm selection uses heuristics:
   *       - Small fields (m ≤ 8): Optimized logarithm tables for fastest
   * operations
   *       - Medium fields (8 < m ≤ 16): Standard logarithm tables balancing
   * speed/memory
   *       - Large fields (m > 16): Polynomial arithmetic to avoid memory
   * overhead
   *       Maximum supported field size: 2^32 elements (uint32_t limit)
   */
  static GaloisFieldVariant CreateBinaryExtensionField(
      uint64_t exponent, const std::string &representation = "int",
      const std::string &modulus = "", const std::string &variable_name = "α",
      const std::string &impl = "auto", bool check_irreducible = false) {
    if (exponent < 2) {
      throw std::invalid_argument("Binary extension degree must be at least 2");
    }

    /* Calculate field order for element type selection */
    uint64_t order = 1ULL << exponent;

    /* Validate field order is within supported range */
    if (order > std::numeric_limits<uint32_t>::max()) {
      throw std::invalid_argument("Binary extension field order exceeds "
                                  "maximum supported size (uint32_t)");
    }

    /* Intelligent implementation selection for optimal performance */
    if (impl == "auto") {
      /* Heuristic-based algorithm selection optimized for different field sizes
       */
      if (exponent <= 8) {
        return CreateBinaryExtensionFieldTyped<uint8_t>(
            exponent, representation, modulus, variable_name,
            "log_opt", check_irreducible);
      } else if (exponent <= 16) {
        return CreateBinaryExtensionFieldTyped<uint16_t>(
            exponent, representation, modulus, variable_name, "log", check_irreducible);
      } else {
        return CreateBinaryExtensionFieldTyped<uint32_t>(
            exponent, representation, modulus, variable_name,
            "standard", check_irreducible);
      }
    }
    // Manual implementation selection with automatic element type optimization
    if (order <= std::numeric_limits<uint8_t>::max()) {
      return CreateBinaryExtensionFieldTyped<uint8_t>(
          exponent, representation, modulus, variable_name, impl, check_irreducible);
    } else if (order <= std::numeric_limits<uint16_t>::max()) {
      return CreateBinaryExtensionFieldTyped<uint16_t>(
          exponent, representation, modulus, variable_name, impl, check_irreducible);
    } else {
      return CreateBinaryExtensionFieldTyped<uint32_t>(
          exponent, representation, modulus, variable_name, impl, check_irreducible);
    }
  }
  // Helper method to create binary extension fields with specific element type
  template <typename ElementType>
  static GaloisFieldVariant CreateBinaryExtensionFieldTyped(
      uint64_t exponent, const std::string &representation,
      const std::string &modulus, const std::string &variable_name,
      const std::string &impl, bool check_irreducible) {
    uint8_t m = static_cast<uint8_t>(exponent);

    if (impl == "standard" || impl == "poly") {
      // Standard polynomial arithmetic implementation
      auto field = std::make_shared<GaloisFieldBinaryExtension<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
      return field;

    } else if (impl == "log" || impl == "log_tables") {
      // Logarithm tables implementation
      auto field = std::make_shared<GFBELogTables<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
      return field;

    } else if (impl == "log_opt" || impl == "log_optimized") {
      // Optimized logarithm tables implementation
      auto field = std::make_shared<GFBELogTablesOpt<ElementType>>(
          m, representation, modulus, variable_name, check_irreducible);
      return field;

    } else if (impl == "zech" || impl == "zech_log") {
      // Zech logarithm tables implementation (only for specific types)
      if constexpr (std::is_same_v<ElementType, uint32_t>) {
        auto field = std::make_shared<GFBEZechLogTables>(m, representation,
                                                         modulus, variable_name, check_irreducible);
        return field;
      } else {
        throw std::invalid_argument("Zech logarithm implementation only "
                                    "supports uint32_t element type");
      }

    } else {
      throw std::invalid_argument(
          "Unknown binary field implementation: " + impl +
          ". Available options: standard, log, log_opt, zech");
    }
  }
}; // class GaloisFieldFactory

// -------------------------------------------------------------
// Interface methods
// -------------------------------------------------------------

/**
 * @brief Creates a Galois finite field GF(p^n) using prime-exponent
 * specification
 *
 * Constructs a finite field with prime base p and extension degree n,
 * automatically selecting the optimal implementation based on the field
 * characteristics. This interface is convenient for working with prime power
 * specifications.
 *
 * @param prime_exp Pair (prime_base, extension_degree) defining GF(p^n)
 * @param representation Display format for field elements:
 *                       - "int": Integer representation (default)
 *                       - "hex": Hexadecimal format
 *                       - "pow": Power notation (α^i)
 *                       - "log": Logarithmic representation
 *                       - "poly": Polynomial notation
 * @param modulus Irreducible polynomial for field extensions as string:
 *                - Required when n > 1
 *                - Format: decimal (e.g., "19") or hex (e.g., "0x13") for
 * binary fields
 *                - Polynomial representation for general fields
 * @param variable_name Symbol for the primitive element in extension fields
 * (e.g., "α", "x")
 *                      - Used only when n > 1 (extension fields)
 *                      - Default: empty string
 * @param impl Implementation algorithm preference:
 *             - "auto": Automatic selection based on field characteristics
 * (recommended)
 *             - "standard"/"poly": Classical polynomial arithmetic
 *             - "log"/"log_tables": Logarithm table optimization
 *             - "log_opt": Enhanced logarithm tables
 *             - "zech": Zech logarithm tables (binary fields only)
 * @param check_irreducible Verify that the modulus polynomial is irreducible:
 *                          - false: Skip verification for performance
 * (default)
 *                          - true: Verify polynomial irreducibility (slower
 * but safer)
 * @param prime_testing Whether to verify that the base prime is actually prime:
 *                      - false: Skip prime verification for performance (default)
 *                      - true: Verify primality of base field characteristic
 *
 * @return GaloisFieldVariant containing the optimized field implementation
 *
 * @throws std::invalid_argument If:
 *         - prime < 2 or exponent < 1
 *         - Field order p^n exceeds uint32_t maximum
 *         - No suitable polynomial found in databases (when modulus not
 * provided)
 *         - Modulus polynomial is not irreducible (when check_irreducible =
 * true)
 *
 * @note Special optimizations:
 *       - GF(2): Uses specialized binary field implementation
 *       - GF(2^m): Automatic algorithm selection for binary extensions
 *       - Element type (uint8_t/uint16_t/uint32_t) chosen automatically for
 * memory efficiency
 *       - Database integration: Conway polynomials preferred, irreducible
 * polynomials as fallback
 *
 * @example
 *   auto gf4 = GF({2, 2});                           // GF(2^2) with Conway
 * polynomial from database auto gf4_manual = GF({2, 2}, "int", "0x7", "α"); //
 * GF(2^2) with manual modulus x^2+x+1 auto gf5 = GF({5, 1}); // GF(5) - prime
 * field auto gf8 = GF({2, 3}, "pow", "", "β");           // GF(2^3) with power
 * representation, database polynomial
 */
GaloisFieldVariant GF(std::pair<uint64_t, uint64_t> prime_exp,
                      const std::string &representation = "int",
                      const std::string &modulus = "",
                      const std::string &variable_name = "",
                      const std::string &impl = "auto",
                      bool check_irreducible = false,
                      bool prime_testing = false) {

  return GaloisFieldFactory::Create(prime_exp, representation, modulus,
                                    variable_name, impl, check_irreducible, prime_testing);
}
/**
 * @brief Creates a Galois finite field GF(q) using total field order
 *
 * Constructs a finite field with exactly q elements by automatically
 * decomposing the order into prime power form (p^n). This interface is
 * convenient when working directly with field sizes rather than
 * prime-exponent specifications.
 *
 * @param order Total number of elements in the field (must be a prime power
 * p^n):
 *              - Must be ≥ 2 and ≤ 2^32
 *              - Examples: 2, 3, 4, 5, 7, 8, 9, 11, 16, 25, 27, 32, ...
 *              - Invalid: 6, 10, 12, 14, 15, 18, 20, ... (not prime powers)
 * @param representation Display format for field elements:
 *                       - "int": Integer representation (default)
 *                       - "hex": Hexadecimal format
 *                       - "pow": Power notation (α^i)
 *                       - "log": Logarithmic representation
 *                       - "poly": Polynomial notation
 * @param modulus Irreducible polynomial for field extensions as string:
 *                - Required when n > 1
 *                - Format: decimal (e.g., "19") or hex (e.g., "0x13") for
 * binary fields
 *                - Polynomial representation for general fields
 * @param variable_name Symbol for the primitive element in extension fields
 * (e.g., "α", "x")
 *                      - Used only when n > 1 (extension fields)
 *                      - Default: empty string
 * @param impl Implementation algorithm preference:
 *             - "auto": Automatic selection based on field characteristics
 * (recommended)
 *             - "standard"/"poly": Classical polynomial arithmetic
 *             - "log"/"log_tables": Logarithm table optimization
 *             - "log_opt": Enhanced logarithm tables
 *             - "zech": Zech logarithm tables (binary fields only)
 * @param check_irreducible Verify that the modulus polynomial is irreducible:
 *                          - false: Skip verification for performance
 * (default)
 *                          - true: Verify polynomial irreducibility (slower
 * but safer)
 * @param prime_testing Whether to verify that the base prime is actually prime:
 *                      - false: Skip prime verification for performance (default)
 *                      - true: Verify primality of base field characteristic
 *
 * @return GaloisFieldVariant containing the optimized field implementation
 *         Can be cast to specific field types or used with std::visit
 *
 * @throws std::invalid_argument If:
 *         - prime < 2 or exponent < 1
 *         - Field order p^n exceeds uint32_t maximum
 *         - No suitable polynomial found in databases (when modulus not
 * provided)
 *         - Modulus polynomial is not irreducible (when check_irreducible =
 * true)
 *
 * @note Special optimizations:
 *       - GF(2): Uses specialized binary field implementation
 *       - GF(2^m): Automatic algorithm selection for binary extensions
 *       - Element type (uint8_t/uint16_t/uint32_t) chosen automatically for
 * memory efficiency
 *       - Database integration: Conway polynomials preferred, irreducible
 * polynomials as fallback
 *
 * @example
 *   auto gf4 = GF(4);                                // GF(2^2) with database
 * polynomial auto gf4_manual = GF(4, "int", "0x7", "α");      // GF(2^2) with
 * manual modulus x^2+x+1 auto gf5 = GF(5);                                //
 * GF(5) - prime field auto gf8 = GF(8, "pow", "", "β");                //
 * GF(2^3) with power representation, database polynomial representation
 */
GaloisFieldVariant GF(uint64_t order, const std::string &representation = "int",
                      const std::string &modulus = "",
                      const std::string &variable_name = "",
                      const std::string &impl = "auto",
                      bool check_irreducible = false,
                      bool prime_testing = false) {

  return GaloisFieldFactory::Create(order, representation, modulus,
                                    variable_name, impl, check_irreducible, prime_testing);
}

} // namespace xg

#endif // XGALOIS_FIELD_GF_FACTORY_HPP_