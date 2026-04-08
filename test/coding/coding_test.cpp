#include "xgalois/coding/coding.hpp"

#include <gtest/gtest.h>

using namespace xg;
using namespace xg::coding;

class CodingTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf2 = std::make_shared<GaloisFieldPrime<uint8_t>>(2);
    gf3 = std::make_shared<GaloisFieldPrime<uint8_t>>(3);
    gf4 = std::make_shared<GaloisFieldBinary<uint8_t>>(2);
    gf8 = std::make_shared<GaloisFieldBinary<uint8_t>>(3);
    gf16 = std::make_shared<GaloisFieldBinary<uint8_t>>(4);
  }

  std::shared_ptr<GaloisFieldBase<uint8_t>> gf2, gf3, gf4, gf8, gf16;
};

TEST_F(CodingTest, ReedSolomonCreation) {
  auto rs_code = CreateReedSolomonCode(gf16, 15, 11);

  EXPECT_EQ(rs_code->Length(), 15);
  EXPECT_EQ(rs_code->Dimension(), 11);
  EXPECT_EQ(rs_code->MinimumDistance(), 5);
  EXPECT_TRUE(IsMDS(*rs_code));
  EXPECT_EQ(ErrorCorrectionCapability(*rs_code), 2);
}

TEST_F(CodingTest, GeneralizedReedSolomon) {
  std::vector<uint8_t> eval_points = {1, 2, 3, 4, 5};
  std::vector<uint8_t> multipliers = {1, 1, 2, 2, 3};

  auto grs_code =
      CreateGeneralizedReedSolomonCode(gf8, eval_points, multipliers, 3);

  EXPECT_EQ(grs_code->Length(), 5);
  EXPECT_EQ(grs_code->Dimension(), 3);
  EXPECT_EQ(grs_code->MinimumDistance(), 3);
  EXPECT_TRUE(IsMDS(*grs_code));
}

TEST_F(CodingTest, CyclicCode) {
  std::vector<uint8_t> gen_coeffs = {1, 1, 0, 1};
  auto gen_poly = xg::poly::PolyDense<uint8_t>(gen_coeffs, gf2);

  auto cyclic_code = CreateCyclicCode(gf2, 7, gen_poly);

  EXPECT_EQ(cyclic_code->Length(), 7);
  EXPECT_EQ(cyclic_code->Dimension(), 4);
  EXPECT_EQ(cyclic_code->GeneratorPolynomial().Degree(), 3);
}

TEST_F(CodingTest, EncodingDecoding) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

  std::vector<uint8_t> message = {1, 2, 3, 4};

  auto codeword = rs_code->Encode(message);
  EXPECT_EQ(codeword.size(), 7);

  EXPECT_TRUE(rs_code->Contains(codeword));

  auto decoded = rs_code->DecodeToMessage(codeword);
  EXPECT_EQ(decoded, message);

  auto unencoded = rs_code->Unencode(codeword);
  EXPECT_EQ(unencoded, message);
}

TEST_F(CodingTest, ErrorCorrection) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 3);

  std::vector<uint8_t> message = {1, 2, 3};
  auto codeword = rs_code->Encode(message);

  auto received = codeword;
  received[0] = (received[0] + 1) % 8;

  auto corrected = rs_code->DecodeToCode(received);
  EXPECT_EQ(corrected, codeword);

  auto decoded_message = rs_code->DecodeToMessage(received);
  EXPECT_EQ(decoded_message, message);
}

TEST_F(CodingTest, CyclicShift) {
  std::vector<uint8_t> gen_coeffs = {1, 1, 0, 1};
  auto gen_poly = xg::poly::PolyDense<uint8_t>(gen_coeffs, gf2);
  auto cyclic_code = CreateCyclicCode(gf2, 7, gen_poly);

  std::vector<uint8_t> message = {1, 0, 1, 1};
  auto codeword = cyclic_code->Encode(message);

  auto shifted = cyclic_code->CyclicShift(codeword, 1);
  EXPECT_TRUE(cyclic_code->Contains(shifted));

  auto shifted2 = cyclic_code->CyclicShift(shifted, 1);
  EXPECT_TRUE(cyclic_code->Contains(shifted2));
}

TEST_F(CodingTest, CodeBounds) {
  auto rs_code = CreateReedSolomonCode(gf16, 15, 11);

  EXPECT_TRUE(SatisfiesSingletonBound(*rs_code));
  EXPECT_TRUE(IsMDS(*rs_code));

  EXPECT_TRUE(SatisfiesHammingBound(*rs_code));

  EXPECT_FALSE(IsPerfectCode(*rs_code));
}

TEST_F(CodingTest, EncoderDecoderRegistration) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

  auto encoders = rs_code->AvailableEncoders();
  EXPECT_GT(encoders.size(), 0);

  auto decoders = rs_code->AvailableDecoders();
  EXPECT_GT(decoders.size(), 0);

  auto encoder = rs_code->GetEncoder("GRSEncoder");
  EXPECT_NE(encoder, nullptr);

  auto decoder = rs_code->GetDecoder("GRSDecoder");
  EXPECT_NE(decoder, nullptr);
}

TEST_F(CodingTest, PolynomialOperations) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

  std::vector<uint8_t> coeffs = {1, 2, 3, 4};
  auto poly = xg::poly::PolyDense<uint8_t>(coeffs, gf8);

  auto codeword = rs_code->EncodePolynomial(poly);
  EXPECT_EQ(codeword.size(), 7);
  EXPECT_TRUE(rs_code->Contains(codeword));

  auto decoded_poly = rs_code->DecodeToPolynomial(codeword);
  EXPECT_EQ(decoded_poly.Degree(), poly.Degree());

  for (size_t i = 0; i <= poly.Degree(); ++i) {
    EXPECT_EQ(decoded_poly.GetCoefficient(i), poly.GetCoefficient(i));
  }
}

TEST_F(CodingTest, DistanceCalculations) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

  std::vector<uint8_t> word1 = {1, 2, 3, 4, 5, 6, 7};
  std::vector<uint8_t> word2 = {1, 2, 3, 0, 0, 6, 7};

  size_t distance = rs_code->Distance(word1, word2);
  EXPECT_EQ(distance, 2);

  std::vector<uint8_t> word3 = {0, 2, 0, 4, 0, 6, 0};
  size_t weight = rs_code->Weight(word3);
  EXPECT_EQ(weight, 3);
}

TEST_F(CodingTest, DualCodes) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 4);
  auto dual_code = rs_code->DualCode();

  EXPECT_EQ(dual_code->Length(), 7);
  EXPECT_EQ(dual_code->Dimension(), 3);

  EXPECT_TRUE(IsMDS(*dual_code));
}

TEST_F(CodingTest, SystematicEncoding) {
  auto rs_code = CreateReedSolomonCode(gf8, 7, 4);

  std::vector<uint8_t> message = {1, 2, 3, 4};

  auto encoder = rs_code->GetEncoder("GRSEncoder");
  auto sys_encoder = dynamic_cast<GRSEncoder<uint8_t> *>(encoder.get());

  if (sys_encoder) {
    auto sys_codeword = sys_encoder->SystematicEncode(message);
    EXPECT_EQ(sys_codeword.size(), 7);
    EXPECT_TRUE(rs_code->Contains(sys_codeword));

    auto sys_positions = rs_code->SystematicPositions();
    for (size_t i = 0; i < message.size(); ++i) {
      EXPECT_EQ(sys_codeword[sys_positions[i]], message[i]);
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
