//===----------------------------------------------------------------------===//
//                          XGalois Library
//===----------------------------------------------------------------------===//
// Copyright (C) 2024 Amir Mulla
//
// Comprehensive unit tests for field utility functions
//===----------------------------------------------------------------------===//

#include "xgalois/utils/field.hpp"
#include <gtest/gtest.h>
#include <stdexcept>

using namespace xg::utils;

//===----------------------------------------------------------------------===//
// Test Fixtures
//===----------------------------------------------------------------------===//

class FieldUtilsTest : public ::testing::Test {
protected:
  void SetUp() override {
    // No special setup needed for field utilities
  }
};

//===----------------------------------------------------------------------===//
// ConvertRepresentation Tests
//===----------------------------------------------------------------------===//

TEST_F(FieldUtilsTest, ConvertRepresentationValidCases) {
  EXPECT_EQ(ConvertRepresentation("int"), xg::FieldRepresentation::INT);
  EXPECT_EQ(ConvertRepresentation("hex"), xg::FieldRepresentation::HEX);
  EXPECT_EQ(ConvertRepresentation("pow"), xg::FieldRepresentation::POW);
  EXPECT_EQ(ConvertRepresentation("log"), xg::FieldRepresentation::LOG);
  EXPECT_EQ(ConvertRepresentation("poly"), xg::FieldRepresentation::POLY);
}

TEST_F(FieldUtilsTest, ConvertRepresentationCaseSensitive) {
  // Test that the function is case-sensitive
  EXPECT_THROW(ConvertRepresentation("INT"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("HEX"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("POW"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("LOG"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("POLY"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ConvertRepresentationInvalidCases) {
  EXPECT_THROW(ConvertRepresentation("invalid"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation(""), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("integer"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("hexadecimal"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("power"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("logarithm"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("polynomial"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ConvertRepresentationExceptionMessage) {
  try {
    ConvertRepresentation("unknown");
    FAIL() << "Expected std::invalid_argument to be thrown";
  } catch (const std::invalid_argument& e) {
    EXPECT_STREQ(e.what(), "Unknown representation: unknown");
  }
}

TEST_F(FieldUtilsTest, ConvertRepresentationWithWhitespace) {
  // Test that whitespace is not trimmed (should throw)
  EXPECT_THROW(ConvertRepresentation(" int"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("int "), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation(" int "), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("i nt"), std::invalid_argument);
}

TEST_F(FieldUtilsTest, ConvertRepresentationSpecialCharacters) {
  // Test with special characters
  EXPECT_THROW(ConvertRepresentation("int\n"), std::invalid_argument);
  EXPECT_THROW(ConvertRepresentation("int\t"), std::invalid_argument);
  // Note: Testing with null character is tricky as it would terminate the string
  // Instead test with other invalid variations
  EXPECT_THROW(ConvertRepresentation("int\r"), std::invalid_argument);
}
