#include <gtest/gtest.h>
#include "coding/coding.hpp"
#include "field/gf_prime.hpp"
#include "field/gf_extension.hpp"

using namespace xg;
using namespace xg::coding;

class CodingTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create various fields for testing
        gf2 = std::make_shared<GaloisFieldPrime<uint8_t>>(2);
        gf3 = std::make_shared<GaloisFieldPrime<uint8_t>>(3);
        gf4 = std::make_shared<GaloisFieldBinary<uint8_t>>(2);
        gf8 = std::make_shared<GaloisFieldBinary<uint8_t>>(3);
        gf16 = std::make_shared<GaloisFieldBinary<uint8_t>>(4);
    }

    std::shared_ptr<GaloisFieldBase<uint8_t>> gf2, gf3, gf4, gf8, gf16;
};

// Test Reed-Solomon code creation and basic properties
TEST_F(CodingTest, ReedSolomonCreation) {
    auto rs_code = CreateReedSolomonCode(gf16, 15, 11);

    EXPECT_EQ(rs_code->Length(), 15);
    EXPECT_EQ(rs_code->Dimension(), 11);
    EXPECT_EQ(rs_code->MinimumDistance(), 5);
    EXPECT_TRUE(IsMDS(*rs_code));
    EXPECT_EQ(ErrorCorrectionCapability(*rs_code), 2);
}

// Test Generalized Reed-Solomon code
TEST_F(CodingTest, GeneralizedReedSolomon) {
    std::vector<uint8_t> eval_points = {1, 2, 3, 4, 5};
    std::vector<uint8_t> multipliers = {1, 1, 2, 2, 3};

    auto grs_code = CreateGeneralizedReedSolomonCode(gf8, eval_points, multipliers, 3);

    EXPECT_EQ(grs_code->Length(), 5);
    EXPECT_EQ(grs_code->Dimension(), 3);
    EXPECT_EQ(grs_code->MinimumDistance(), 3);
    EXPECT_TRUE(IsMDS(*grs_code));
}

// Test cyclic code creation
TEST_F(CodingTest, CyclicCode) {
    // Create generator polynomial g(x) = x^3 + x + 1
    std::vector<uint8_t> gen_coeffs = {1, 1, 0, 1};
    auto gen_poly = xg::poly::PolyDense<uint8_t>(gen_coeffs, gf2);

    auto cyclic_code = CreateCyclicCode(gf2, 7, gen_poly);

    EXPECT_EQ(cyclic_code->Length(), 7);
    EXPECT_EQ(cyclic_code->Dimension(), 4);
    EXPECT_EQ(cyclic_code->GeneratorPolynomial().Degree(), 3);
}

// Test encoding and decoding
TEST_F(CodingTest, EncodingDecoding) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

    std::vector<uint8_t> message = {1, 2, 3, 4};

    // Test encoding
    auto codeword = rs_code->Encode(message);
    EXPECT_EQ(codeword.size(), 7);

    // Test that encoded message is a valid codeword
    EXPECT_TRUE(rs_code->Contains(codeword));

    // Test decoding
    auto decoded = rs_code->DecodeToMessage(codeword);
    EXPECT_EQ(decoded, message);

    // Test unencode
    auto unencoded = rs_code->Unencode(codeword);
    EXPECT_EQ(unencoded, message);
}

// Test error correction
TEST_F(CodingTest, ErrorCorrection) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 3);

    std::vector<uint8_t> message = {1, 2, 3};
    auto codeword = rs_code->Encode(message);

    // Introduce single error
    auto received = codeword;
    received[0] = (received[0] + 1) % 8;

    // Should be able to correct single error
    auto corrected = rs_code->DecodeToCode(received);
    EXPECT_EQ(corrected, codeword);

    auto decoded_message = rs_code->DecodeToMessage(received);
    EXPECT_EQ(decoded_message, message);
}

// Test cyclic shift property
TEST_F(CodingTest, CyclicShift) {
    std::vector<uint8_t> gen_coeffs = {1, 1, 0, 1}; // x^3 + x + 1
    auto gen_poly = xg::poly::PolyDense<uint8_t>(gen_coeffs, gf2);
    auto cyclic_code = CreateCyclicCode(gf2, 7, gen_poly);

    std::vector<uint8_t> message = {1, 0, 1, 1};
    auto codeword = cyclic_code->Encode(message);

    // Test cyclic shift
    auto shifted = cyclic_code->CyclicShift(codeword, 1);
    EXPECT_TRUE(cyclic_code->Contains(shifted));

    // Test multiple shifts
    auto shifted2 = cyclic_code->CyclicShift(shifted, 1);
    EXPECT_TRUE(cyclic_code->Contains(shifted2));
}

// Test code bounds
TEST_F(CodingTest, CodeBounds) {
    auto rs_code = CreateReedSolomonCode(gf16, 15, 11);

    // Reed-Solomon codes are MDS, so they meet the Singleton bound with equality
    EXPECT_TRUE(SatisfiesSingletonBound(*rs_code));
    EXPECT_TRUE(IsMDS(*rs_code));

    // Should also satisfy Hamming bound
    EXPECT_TRUE(SatisfiesHammingBound(*rs_code));

    // RS codes are not perfect (except for trivial cases)
    EXPECT_FALSE(IsPerfectCode(*rs_code));
}

// Test encoder/decoder registration
TEST_F(CodingTest, EncoderDecoderRegistration) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

    // Check that encoders are registered
    auto encoders = rs_code->AvailableEncoders();
    EXPECT_GT(encoders.size(), 0);

    // Check that decoders are registered
    auto decoders = rs_code->AvailableDecoders();
    EXPECT_GT(decoders.size(), 0);

    // Test getting specific encoder
    auto encoder = rs_code->GetEncoder("GRSEncoder");
    EXPECT_NE(encoder, nullptr);

    // Test getting specific decoder
    auto decoder = rs_code->GetDecoder("GRSDecoder");
    EXPECT_NE(decoder, nullptr);
}

// Test polynomial operations in GRS codes
TEST_F(CodingTest, PolynomialOperations) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

    // Create a polynomial
    std::vector<uint8_t> coeffs = {1, 2, 3, 4}; // 1 + 2x + 3x^2 + 4x^3
    auto poly = xg::poly::PolyDense<uint8_t>(coeffs, gf8);

    // Test polynomial encoding
    auto codeword = rs_code->EncodePolynomial(poly);
    EXPECT_EQ(codeword.size(), 7);
    EXPECT_TRUE(rs_code->Contains(codeword));

    // Test polynomial decoding
    auto decoded_poly = rs_code->DecodeToPolynomial(codeword);
    EXPECT_EQ(decoded_poly.Degree(), poly.Degree());

    // Coefficients should match
    for (size_t i = 0; i <= poly.Degree(); ++i) {
        EXPECT_EQ(decoded_poly.GetCoefficient(i), poly.GetCoefficient(i));
    }
}

// Test distance calculations
TEST_F(CodingTest, DistanceCalculations) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

    std::vector<uint8_t> word1 = {1, 2, 3, 4, 5, 6, 7};
    std::vector<uint8_t> word2 = {1, 2, 3, 0, 0, 6, 7};

    // Should have Hamming distance 2
    size_t distance = rs_code->Distance(word1, word2);
    EXPECT_EQ(distance, 2);

    // Test weight calculation
    std::vector<uint8_t> word3 = {0, 2, 0, 4, 0, 6, 0};
    size_t weight = rs_code->Weight(word3);
    EXPECT_EQ(weight, 3);
}

// Test dual codes
TEST_F(CodingTest, DualCodes) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 4);
    auto dual_code = rs_code->DualCode();

    // Dual of [n, k] code is [n, n-k]
    EXPECT_EQ(dual_code->Length(), 7);
    EXPECT_EQ(dual_code->Dimension(), 3);

    // For Reed-Solomon codes, dual is also Reed-Solomon
    EXPECT_TRUE(IsMDS(*dual_code));
}

// Test systematic encoding
TEST_F(CodingTest, SystematicEncoding) {
    auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

    std::vector<uint8_t> message = {1, 2, 3, 4};

    // Get systematic encoder
    auto encoder = rs_code->GetEncoder("GRSEncoder");
    auto sys_encoder = dynamic_cast<GRSEncoder<uint8_t>*>(encoder.get());

    if (sys_encoder) {
        auto sys_codeword = sys_encoder->SystematicEncode(message);
        EXPECT_EQ(sys_codeword.size(), 7);
        EXPECT_TRUE(rs_code->Contains(sys_codeword));

        // For systematic encoding, message should appear in specific positions
        auto sys_positions = rs_code->SystematicPositions();
        for (size_t i = 0; i < message.size(); ++i) {
            EXPECT_EQ(sys_codeword[sys_positions[i]], message[i]);
        }
    }
}

// Run all tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
