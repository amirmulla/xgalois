#include "xgalois/field/gf_factory.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <memory>
#include <set>
#include <variant>
#include <vector>

using namespace xg;

//===----------------------------------------------------------------------===//
// GaloisField Factory Tests
//===----------------------------------------------------------------------===//

class GaloisFieldFactoryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Test various field orders that should trigger different implementations
  }
};

TEST_F(GaloisFieldFactoryTest, PrimeFields) {
  // Test creation of various prime fields
  auto gf2 = GF(2);
  auto gf3 = GF(3);
  auto gf5 = GF(5);
  auto gf7 = GF(7);
  auto gf11 = GF(11);
  auto gf257 = GF(257);

  // Verify they are the correct field types using std::visit
  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 2);
        EXPECT_EQ(field->Order(), 2);
      },
      gf2);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 7);
        EXPECT_EQ(field->Order(), 7);
      },
      gf7);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 257);
        EXPECT_EQ(field->Order(), 257);
      },
      gf257);
}

TEST_F(GaloisFieldFactoryTest, BinaryExtensionFields) {
  // Test creation of binary extension fields
  auto gf4 = GF(4);      // GF(2^2)
  auto gf8 = GF(8);      // GF(2^3)
  auto gf16 = GF(16);    // GF(2^4)
  auto gf32 = GF(32);    // GF(2^5)
  auto gf256 = GF(256);  // GF(2^8)

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 2);
        EXPECT_EQ(field->Order(), 4);
      },
      gf4);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 2);
        EXPECT_EQ(field->Order(), 8);
      },
      gf8);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 2);
        EXPECT_EQ(field->Order(), 256);
      },
      gf256);
}

TEST_F(GaloisFieldFactoryTest, CustomRepresentation) {
  // Test fields with different representations
  auto gf7_int = GF(7, "int");
  auto gf7_hex = GF(7, "hex");
  auto gf7_pow = GF(7, "pow");
  auto gf7_log = GF(7, "log");

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::INT);
      },
      gf7_int);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::HEX);
      },
      gf7_hex);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::POW);
      },
      gf7_pow);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::LOG);
      },
      gf7_log);
}

TEST_F(GaloisFieldFactoryTest, ElementCreation) {
  auto gf7 = GF(7);
  auto gf8 = GF(8);

  // Test element creation for supported field types
  auto elem_gf7 = FetchElement(gf7, 3);
  auto elem_gf8 = FetchElement(gf8, 5);

  // Verify elements have correct values
  std::visit([](const auto& elem) { EXPECT_EQ(elem.Value(), 3); }, elem_gf7);

  std::visit([](const auto& elem) { EXPECT_EQ(elem.Value(), 5); }, elem_gf8);
}

TEST_F(GaloisFieldFactoryTest, StrictElementCreation) {
  auto gf7 = GF(7);
  auto gf8 = GF(8);

  // Test strict element creation
  auto elem_strict_gf7 = FetchElementStrict(gf7, 3);
  auto elem_strict_gf8 = FetchElementStrict(gf8, 5);

  // Verify elements have correct values
  std::visit([](const auto& elem) { EXPECT_EQ(elem.Value(), 3); },
             elem_strict_gf7);

  std::visit([](const auto& elem) { EXPECT_EQ(elem.Value(), 5); },
             elem_strict_gf8);
}

TEST_F(GaloisFieldFactoryTest, InvalidFieldOrders) {
  // Test invalid field orders
  EXPECT_THROW(GF(0), std::invalid_argument);
  EXPECT_THROW(GF(1), std::invalid_argument);
  EXPECT_THROW(GF(6), std::invalid_argument);   // 6 = 2*3, not a prime power
  EXPECT_THROW(GF(10), std::invalid_argument);  // 10 = 2*5, not a prime power
  EXPECT_THROW(GF(12), std::invalid_argument);  // 12 = 2^2*3, not a prime power
}

TEST_F(GaloisFieldFactoryTest, LargeFields) {
  // Test creation of larger fields (but not too large for quick tests)
  auto gf1024 = GF(1024);  // GF(2^10)
  auto gf1009 = GF(1009);  // Large prime

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 2);
        EXPECT_EQ(field->Order(), 1024);
      },
      gf1024);

  std::visit(
      [](const auto& field) {
        EXPECT_EQ(field->Characteristic(), 1009);
        EXPECT_EQ(field->Order(), 1009);
      },
      gf1009);
}

TEST_F(GaloisFieldFactoryTest, ElementTypeOptimization) {
  // Test that appropriate element types are chosen based on field size
  auto gf2 = GF(2);          // Should use uint8_t
  auto gf256 = GF(256);      // Should use uint8_t
  auto gf257 = GF(257);      // Should use uint16_t
  auto gf65536 = GF(65536);  // Should use uint16_t or uint32_t

  // All should be created successfully with appropriate types
  std::visit([](const auto& field) { EXPECT_NE(field, nullptr); }, gf2);

  std::visit([](const auto& field) { EXPECT_NE(field, nullptr); }, gf256);

  std::visit([](const auto& field) { EXPECT_NE(field, nullptr); }, gf257);

  std::visit([](const auto& field) { EXPECT_NE(field, nullptr); }, gf65536);
}

TEST_F(GaloisFieldFactoryTest, FieldOperations) {
  auto gf7 = GF(7);
  auto gf8 = GF(8);

  // Test basic operations through variant interface
  std::visit(
      [](const auto& field) {
        auto a = field->MultiplicativeIdentity();
        auto b = field->AdditiveIdentity();
        auto sum = field->Add(a, b);
        EXPECT_EQ(sum, a);  // 1 + 0 = 1

        auto prod = field->Mul(a, a);
        EXPECT_EQ(prod, a);  // 1 * 1 = 1
      },
      gf7);

  std::visit(
      [](const auto& field) {
        auto a = field->MultiplicativeIdentity();
        auto b = field->AdditiveIdentity();
        auto sum = field->Add(a, b);
        EXPECT_EQ(sum, a);  // 1 + 0 = 1
      },
      gf8);
}

TEST_F(GaloisFieldFactoryTest, PerformanceTest) {
  // Quick performance test to ensure factory creation is efficient
  auto start = std::chrono::high_resolution_clock::now();

  // Create multiple fields
  std::vector<GaloisFieldVariant> fields;
  for (int i = 0; i < 10; ++i) {
    fields.push_back(GF(7));
    fields.push_back(GF(8));
    fields.push_back(GF(16));
    fields.push_back(GF(31));
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  EXPECT_LT(duration.count(),
            30000);  // Should complete in less than 30 seconds
  EXPECT_EQ(fields.size(), 40);
}

TEST_F(GaloisFieldFactoryTest, PrimeFactorization) {
  // Test that the factory correctly identifies prime powers
  // These should work (prime powers)
  EXPECT_NO_THROW(GF(4));  // 2^2 - binary extension field
  EXPECT_NO_THROW(GF(
      9,
      "poly"));  // 3^2 - general extension field (requires poly representation)
  EXPECT_NO_THROW(GF(
      25,
      "poly"));  // 5^2 - general extension field (requires poly representation)
  EXPECT_NO_THROW(GF(
      27,
      "poly"));  // 3^3 - general extension field (requires poly representation)
  EXPECT_NO_THROW(GF(
      49,
      "poly"));  // 7^2 - general extension field (requires poly representation)
  EXPECT_NO_THROW(GF(
      125,
      "poly"));  // 5^3 - general extension field (requires poly representation)

  // These should fail (not prime powers)
  EXPECT_THROW(GF(6), std::invalid_argument);   // 2*3
  EXPECT_THROW(GF(15), std::invalid_argument);  // 3*5
  EXPECT_THROW(GF(21), std::invalid_argument);  // 3*7
  EXPECT_THROW(GF(30), std::invalid_argument);  // 2*3*5
}

TEST_F(GaloisFieldFactoryTest, VariantTypeHandling) {
  auto gf7 = GF(7);
  auto gf8 = GF(8);

  // Test that we can distinguish between different field types
  bool is_prime_field = std::visit(
      [](const auto& field) -> bool {
        using FieldType = std::decay_t<decltype(*field)>;
        return std::is_same_v<FieldType, GaloisFieldPrime<uint8_t>> ||
               std::is_same_v<FieldType, GaloisFieldPrime<uint16_t>> ||
               std::is_same_v<FieldType, GaloisFieldPrime<uint32_t>>;
      },
      gf7);

  bool is_binary_extension = std::visit(
      [](const auto& field) -> bool {
        using FieldType = std::decay_t<decltype(*field)>;
        // Check for all possible binary extension field implementations
        return std::is_same_v<FieldType, GaloisFieldBinaryExtension<uint8_t>> ||
               std::is_same_v<FieldType,
                              GaloisFieldBinaryExtension<uint16_t>> ||
               std::is_same_v<FieldType,
                              GaloisFieldBinaryExtension<uint32_t>> ||
               std::is_same_v<FieldType, GFBELogTables<uint8_t>> ||
               std::is_same_v<FieldType, GFBELogTables<uint16_t>> ||
               std::is_same_v<FieldType, GFBELogTables<uint32_t>> ||
               std::is_same_v<FieldType, GFBELogTablesOpt<uint8_t>> ||
               std::is_same_v<FieldType, GFBELogTablesOpt<uint16_t>> ||
               std::is_same_v<FieldType, GFBELogTablesOpt<uint32_t>> ||
               std::is_same_v<FieldType, GFBEZechLogTables<uint32_t>>;
      },
      gf8);

  EXPECT_TRUE(is_prime_field);
  EXPECT_TRUE(is_binary_extension);
}
