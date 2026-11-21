#include "xgalois/field/gf_factory.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

using namespace xg;

namespace {

uint64_t ExtractElementInt(const GaloisFieldElementVariant &elem_var) {
  return std::visit(
      [](const auto &elem) -> uint64_t {
        using ValueType = decltype(elem.Value());
        if constexpr (std::is_integral_v<ValueType>) {
          return static_cast<uint64_t>(elem.Value());
        } else {
          return elem.Field()->ToInt(elem.Value());
        }
      },
      elem_var);
}

}  // namespace

class GaloisFieldFactoryTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(GaloisFieldFactoryTest, PrimeFields) {
  auto gf2 = GF(2);
  auto gf3 = GF(3);
  auto gf7 = GF(7);
  auto gf257 = GF(257);
  auto gf65537 = GF(65537);

  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->Characteristic(), 2U);
        EXPECT_EQ(field->Order(), 2U);
        EXPECT_EQ(field->Degree(), 1U);
      },
      gf2);

  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->Characteristic(), 3U);
        EXPECT_EQ(field->Order(), 3U);
      },
      gf3);

  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->Characteristic(), 7U);
        EXPECT_EQ(field->Order(), 7U);
      },
      gf7);

  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->Characteristic(), 257U);
        EXPECT_EQ(field->Order(), 257U);
      },
      gf257);

  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->Characteristic(), 65537U);
        EXPECT_EQ(field->Order(), 65537U);
      },
      gf65537);

  EXPECT_TRUE(std::holds_alternative<std::shared_ptr<GaloisFieldBinary>>(gf2));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrime<uint8_t>>>(gf7));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrime<uint16_t>>>(
          gf257));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrime<uint32_t>>>(
          gf65537));
}

TEST_F(GaloisFieldFactoryTest, PrimeFieldImplementations) {
  auto gf7_std = GF(7, "int", "", "α", "standard");
  auto gf7_tbl = GF(7, "int", "", "α", "table");
  auto gf7_log = GF(7, "int", "", "α", "log");
  auto gf257_tbl = GF(257, "int", "", "α", "table");

  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrime<uint8_t>>>(
          gf7_std));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrimeTable<uint8_t>>>(
          gf7_tbl));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrimeTable<uint8_t>>>(
          gf7_log));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldPrimeTable<uint16_t>>>(
          gf257_tbl));

  auto f_tbl =
      std::get<std::shared_ptr<GaloisFieldPrimeTable<uint8_t>>>(gf7_tbl);
  EXPECT_EQ(f_tbl->Mul(3, 5), 1U);
  EXPECT_EQ(f_tbl->Inv(3), 5U);

  EXPECT_THROW(GF(7, "int", "", "α", "unknown_impl"), std::invalid_argument);
}

TEST_F(GaloisFieldFactoryTest, BinaryExtensionDegrees1To20AndLarge) {
  for (uint64_t m = 1; m <= 20; ++m) {
    const uint64_t expected_order = 1ULL << m;
    auto field_by_order = GF(expected_order);
    auto field_by_pair = GF(std::make_pair(uint64_t{2}, m));

    std::visit(
        [m, expected_order](const auto &field) {
          EXPECT_EQ(field->Characteristic(), 2U);
          EXPECT_EQ(field->Degree(), static_cast<uint32_t>(m));
          EXPECT_EQ(field->Order(), expected_order);
        },
        field_by_order);

    std::visit(
        [m, expected_order](const auto &field) {
          EXPECT_EQ(field->Characteristic(), 2U);
          EXPECT_EQ(field->Degree(), static_cast<uint32_t>(m));
          EXPECT_EQ(field->Order(), expected_order);
        },
        field_by_pair);
  }

  // Exact word-width selection at boundaries m = 8, 16, 32, 64.
  auto gf256 = GF(std::make_pair(uint64_t{2}, uint64_t{8}));
  auto gf65536 = GF(std::make_pair(uint64_t{2}, uint64_t{16}));
  auto gf2_32 = GF(std::make_pair(uint64_t{2}, uint64_t{32}));
  auto gf2_64 = GF(std::make_pair(uint64_t{2}, uint64_t{64}));

  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBELogTablesOpt<uint8_t>>>(
          gf256));
  EXPECT_TRUE(std::holds_alternative<std::shared_ptr<GFBELogTables<uint16_t>>>(
      gf65536));
  EXPECT_TRUE(std::holds_alternative<
              std::shared_ptr<GaloisFieldBinaryExtension<uint32_t>>>(gf2_32));
  EXPECT_TRUE(std::holds_alternative<
              std::shared_ptr<GaloisFieldBinaryExtension<uint64_t>>>(gf2_64));

  auto f32 =
      std::get<std::shared_ptr<GaloisFieldBinaryExtension<uint32_t>>>(gf2_32);
  EXPECT_EQ(f32->Degree(), 32U);
  EXPECT_EQ(f32->Order(), uint64_t{1} << 32);
  EXPECT_EQ(f32->Mul(0x12345678U, f32->Inv(0x12345678U)), 1U);

  auto f64 =
      std::get<std::shared_ptr<GaloisFieldBinaryExtension<uint64_t>>>(gf2_64);
  EXPECT_EQ(f64->Characteristic(), 2U);
  EXPECT_EQ(f64->Degree(), 64U);
  EXPECT_EQ(f64->Order(), 0U);
  constexpr uint64_t kSample64 = 0x123456789abcdef1ULL;
  EXPECT_EQ(f64->Mul(kSample64, f64->Inv(kSample64)), 1ULL);
}

TEST_F(GaloisFieldFactoryTest, BinaryExtensionAllImplModes) {
  // Test standard, log, log_opt, and zech across uint8_t, uint16_t, uint32_t.
  auto gf8_std = GF(8, "int", "", "α", "standard");
  auto gf8_log = GF(8, "int", "", "α", "log");
  auto gf8_opt = GF(8, "int", "", "α", "log_opt");
  auto gf8_zech = GF(8, "int", "", "α", "zech");

  EXPECT_TRUE(std::holds_alternative<
              std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>>>(gf8_std));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBELogTables<uint8_t>>>(gf8_log));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBELogTablesOpt<uint8_t>>>(
          gf8_opt));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBEZechLogTables<uint8_t>>>(
          gf8_zech));

  // Verify Zech works for m <= 15 (uint8_t and uint16_t) AND m = 17 (uint32_t).
  auto gf256_zech = GF(256, "int", "", "α", "zech");
  auto gf4096_zech = GF(4096, "int", "", "α", "zech");
  auto gf2_17_zech =
      GF(std::make_pair(uint64_t{2}, uint64_t{17}), "int", "", "α", "zech");

  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBEZechLogTables<uint8_t>>>(
          gf256_zech));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBEZechLogTables<uint16_t>>>(
          gf4096_zech));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GFBEZechLogTables<uint32_t>>>(
          gf2_17_zech));

  // Verify arithmetic consistency between standard and zech in GF(8).
  auto f_std =
      std::get<std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>>>(gf8_std);
  auto f_zech = std::get<std::shared_ptr<GFBEZechLogTables<uint8_t>>>(gf8_zech);
  for (uint8_t a = 0; a < 8; ++a) {
    for (uint8_t b = 0; b < 8; ++b) {
      const uint8_t za = f_zech->SetElementValue(a);
      const uint8_t zb = f_zech->SetElementValue(b);
      EXPECT_EQ(f_zech->GetElementValue(f_zech->Add(za, zb)), f_std->Add(a, b));
      EXPECT_EQ(f_zech->GetElementValue(f_zech->Mul(za, zb)), f_std->Mul(a, b));
    }
  }

  // Zech > 20 rejected; unknown impl rejected.
  EXPECT_THROW(
      GF(std::make_pair(uint64_t{2}, uint64_t{21}), "int", "", "α", "zech"),
      std::invalid_argument);
  EXPECT_THROW(GF(8, "int", "", "α", "bogus_impl"), std::invalid_argument);
}

TEST_F(GaloisFieldFactoryTest, DefaultVariableNameAndStringParsing) {
  // Bug B-18 regression: GF(8) must default variable_name to "α" (not "").
  auto gf8 = GF(8);
  std::visit(
      [](const auto &field) {
        using FieldType = std::decay_t<decltype(*field)>;
        if constexpr (!detail::kIsGaloisFieldExtension<FieldType> &&
                      !std::is_same_v<FieldType, GaloisFieldBinary>) {
          auto elem = field->SetElementValue("α + 1");
          EXPECT_EQ(field->GetElementValue(elem), 3U);
        }
      },
      gf8);
}

TEST_F(GaloisFieldFactoryTest, PrimeExtensionFields) {
  // Default "int" representation works without aborting, as does "poly".
  auto gf9_default = GF(9);
  auto gf9_poly = GF(9, "poly");
  auto gf25 = GF(25);
  auto gf27 = GF(std::make_pair(uint64_t{3}, uint64_t{3}));
  auto gf729 = GF(std::make_pair(uint64_t{3}, uint64_t{6}));
  auto gf257_sq =
      GF(std::make_pair(uint64_t{257}, uint64_t{2}), "poly", "x^2 + 3");

  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(
          gf9_default));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(
          gf9_poly));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(
          gf25));
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(
          gf27));
  // GF(3^6) has order 729 > 255, but prime 3 <= 255 -> selects uint8_t!
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(
          gf729));
  // Prime 257 > 255 -> selects uint16_t.
  EXPECT_TRUE(
      std::holds_alternative<std::shared_ptr<GaloisFieldExtension<uint16_t>>>(
          gf257_sq));

  auto f9 =
      std::get<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(gf9_default);
  EXPECT_EQ(f9->Characteristic(), 3U);
  EXPECT_EQ(f9->Order(), 9U);
  EXPECT_EQ(f9->Degree(), 2U);
  EXPECT_EQ(f9->GetRepresentation(), FieldRepresentation::INT);

  auto f27 = std::get<std::shared_ptr<GaloisFieldExtension<uint8_t>>>(gf27);
  EXPECT_EQ(f27->Characteristic(), 3U);
  EXPECT_EQ(f27->Order(), 27U);
  EXPECT_EQ(f27->Degree(), 3U);
}

TEST_F(GaloisFieldFactoryTest, CustomRepresentations) {
  auto gf7_int = GF(7, "int");
  auto gf7_hex = GF(7, "hex");
  auto gf7_pow = GF(7, "pow");
  auto gf7_log = GF(7, "log");

  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::INT);
      },
      gf7_int);
  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::HEX);
      },
      gf7_hex);
  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::POW);
      },
      gf7_pow);
  std::visit(
      [](const auto &field) {
        EXPECT_EQ(field->GetRepresentation(), FieldRepresentation::LOG);
      },
      gf7_log);
}

TEST_F(GaloisFieldFactoryTest, FetchElementReductionAndExtensionSupport) {
  auto gf7 = GF(7);
  auto gf8 = GF(8);
  auto gf9 = GF(9);
  auto gf2_64 = GF(std::make_pair(uint64_t{2}, uint64_t{64}));

  EXPECT_EQ(ExtractElementInt(FetchElement(gf7, 3)), 3U);
  EXPECT_EQ(ExtractElementInt(FetchElement(gf8, 5)), 5U);
  EXPECT_EQ(ExtractElementInt(FetchElementStrict(gf7, 3)), 3U);
  EXPECT_EQ(ExtractElementInt(FetchElementStrict(gf8, 5)), 5U);

  // Reduction modulo field order BEFORE narrowing to uint8_t: 256 % 7 == 4.
  EXPECT_EQ(ExtractElementInt(FetchElement(gf7, 256)), 4U);
  EXPECT_EQ(ExtractElementInt(FetchElementStrict(gf7, 256)), 4U);

  // Extension field FetchElement via FromInt(val).
  auto elem_gf9 = FetchElement(gf9, 5);
  EXPECT_EQ(ExtractElementInt(elem_gf9), 5U);
  auto elem_gf9_wrapped = FetchElement(gf9, 14);  // 14 % 9 == 5
  EXPECT_EQ(ExtractElementInt(elem_gf9_wrapped), 5U);

  // GF(2^64) FetchElement (Order() == 0) must not divide by zero.
  constexpr uint64_t kVal64 = 0x123456789abcdef0ULL;
  EXPECT_EQ(ExtractElementInt(FetchElement(gf2_64, kVal64)), kVal64);
}

TEST_F(GaloisFieldFactoryTest, InvalidOrdersAndExponents) {
  EXPECT_THROW(GF(0), std::invalid_argument);
  EXPECT_THROW(GF(1), std::invalid_argument);
  EXPECT_THROW(GF(6), std::invalid_argument);
  EXPECT_THROW(GF(10), std::invalid_argument);
  EXPECT_THROW(GF(12), std::invalid_argument);
  EXPECT_THROW(GF(15), std::invalid_argument);
  EXPECT_THROW(GF(std::make_pair(uint64_t{2}, uint64_t{0})),
               std::invalid_argument);
  EXPECT_THROW(GF(std::make_pair(uint64_t{2}, uint64_t{65})),
               std::invalid_argument);
  EXPECT_THROW(GF(std::make_pair(uint64_t{4}, uint64_t{2})),
               std::invalid_argument);
  EXPECT_THROW(GF(std::make_pair(uint64_t{1}, uint64_t{3})),
               std::invalid_argument);
}
