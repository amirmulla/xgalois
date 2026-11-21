#include "xgalois/coding/coding.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"

namespace {

using namespace xg;
using namespace xg::coding;

template <typename GaloisField>
GFArray<GaloisField> MakeWord(const std::shared_ptr<GaloisField> &field,
                              const std::vector<uint64_t> &vals) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;
  RawArray raw = RawArray::from_shape(ShapeType{vals.size()});
  for (std::size_t i = 0; i < vals.size(); ++i) {
    if constexpr (std::is_integral_v<StorageType>) {
      raw(i) = field->SetElementValue(static_cast<StorageType>(vals[i]));
    } else {
      raw(i) = field->FromInt(vals[i]);
    }
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
}

TEST(CodingTest, CreateReedSolomonOverPrimeAndBinaryFields) {
  // GF(7): 2 is NOT a primitive root mod 7 (order 3); MultiplicativeGenerator()
  // returns a true primitive root (e.g., 3), so length 6 succeeds cleanly!
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto rs7 = CreateReedSolomonCode(gf7, 6, 4);
  EXPECT_EQ(rs7.Length(), 6U);
  EXPECT_EQ(rs7.Dimension(), 4U);
  EXPECT_EQ(rs7.MinimumDistance(), 3U);
  EXPECT_EQ(Redundancy(rs7), 2U);
  EXPECT_EQ(ErrorCorrectionCapability(rs7), 1U);
  EXPECT_TRUE(IsMDSCode(rs7));
  EXPECT_TRUE(SatisfiesSingletonBound(rs7));
  EXPECT_TRUE(SatisfiesHammingBound(rs7));

  // GF(17): 2 has order 8 < 16, whereas MultiplicativeGenerator() has order 16.
  auto gf17 = std::make_shared<GaloisFieldPrime<uint16_t>>(17);
  auto rs17 = CreateReedSolomonCode(gf17, 16, 12);
  EXPECT_EQ(rs17.Length(), 16U);
  EXPECT_EQ(rs17.Dimension(), 12U);
  EXPECT_EQ(rs17.MinimumDistance(), 5U);
  EXPECT_EQ(ErrorCorrectionCapability(rs17), 2U);
  EXPECT_TRUE(IsMDSCode(rs17));

  // GF(2^4)
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto rs16 = CreateReedSolomonCode(gf16, 15, 11);
  EXPECT_EQ(rs16.Length(), 15U);
  EXPECT_EQ(rs16.Dimension(), 11U);
  EXPECT_EQ(rs16.MinimumDistance(), 5U);
  EXPECT_EQ(ErrorCorrectionCapability(rs16), 2U);
  EXPECT_TRUE(IsMDSCode(rs16));
  EXPECT_TRUE(SatisfiesHammingBound(rs16));
  EXPECT_FALSE(IsPerfectCode(rs16));
}

TEST(CodingTest, HammingBoundAndPerfectCodes) {
  // Binary Hamming [7, 4, 3]: V_2(7, 1) = 1 + 7 = 8 = 2^3 => perfect!
  auto hamming7 = HammingCode(3);
  EXPECT_EQ(hamming7.Length(), 7U);
  EXPECT_EQ(hamming7.Dimension(), 4U);
  EXPECT_EQ(hamming7.MinimumDistance(), 3U);
  EXPECT_EQ(ErrorCorrectionCapability(hamming7), 1U);
  EXPECT_TRUE(SatisfiesSingletonBound(hamming7));
  EXPECT_FALSE(IsMDSCode(hamming7));
  EXPECT_TRUE(SatisfiesHammingBound(hamming7));
  EXPECT_TRUE(IsPerfectCode(hamming7));

  // Extended binary Hamming [8, 4, 4]: t = 1, V_2(8, 1) = 9 < 16 => not
  // perfect.
  auto ext_hamming8 = ExtendedHammingCode(3);
  EXPECT_TRUE(SatisfiesHammingBound(ext_hamming8));
  EXPECT_FALSE(IsPerfectCode(ext_hamming8));

  // Binary Golay [23, 12, 7]: V_2(23, 3) = 1 + 23 + 253 + 1771 = 2048 = 2^11 =>
  // perfect!
  auto golay23 = BinaryGolayCode();
  EXPECT_EQ(golay23.Length(), 23U);
  EXPECT_EQ(golay23.Dimension(), 12U);
  EXPECT_EQ(golay23.MinimumDistance(), 7U);
  EXPECT_EQ(ErrorCorrectionCapability(golay23), 3U);
  EXPECT_TRUE(SatisfiesSingletonBound(golay23));
  EXPECT_TRUE(SatisfiesHammingBound(golay23));
  EXPECT_TRUE(IsPerfectCode(golay23));

  // Extended binary Golay [24, 12, 8]: not perfect.
  auto golay24 = ExtendedBinaryGolayCode();
  EXPECT_TRUE(SatisfiesHammingBound(golay24));
  EXPECT_FALSE(IsPerfectCode(golay24));

  // Ternary Golay [11, 6, 5]: V_3(11, 2) = 1 + 22 + 220 = 243 = 3^5 => perfect!
  auto golay11 = TernaryGolayCode();
  EXPECT_EQ(golay11.MinimumDistance(), 5U);
  EXPECT_EQ(ErrorCorrectionCapability(golay11), 2U);
  EXPECT_TRUE(SatisfiesHammingBound(golay11));
  EXPECT_TRUE(IsPerfectCode(golay11));

  // Large RS(255, 223) over GF(2^8): 256^32 = 2^256 overflows uint64/uint128;
  // verify overflow-free exact Hamming bound check succeeds.
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
  auto rs255 = CreateReedSolomonCode(gf256, 255, 223);
  EXPECT_EQ(SingletonBound(255, 223), 33U);
  EXPECT_TRUE(SatisfiesSingletonBound(rs255));
  EXPECT_TRUE(IsMDSCode(rs255));
  EXPECT_TRUE(SatisfiesHammingBound(rs255));
  EXPECT_FALSE(IsPerfectCode(rs255));
}

TEST(CodingTest, EndToEndPrimeAndBinaryExtensionCodes) {
  std::mt19937_64 rng(12345);

  // Prime field GF(17) RS(16, 10) -> corrects up to 3 errors
  {
    auto gf17 = std::make_shared<GaloisFieldPrime<uint16_t>>(17);
    auto rs = CreateReedSolomonCode(gf17, 16, 10);
    auto msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf17, {10}, rng);
    auto cw = rs.Encode(msg);
    EXPECT_TRUE(rs.Contains(cw));
    EXPECT_EQ(rs.Unencode(cw), msg);

    auto rx = cw;
    rx.Raw(1) = gf17->Add(rx.Raw(1), gf17->SetElementValue(5));
    rx.Raw(7) = gf17->Add(rx.Raw(7), gf17->SetElementValue(11));
    rx.Raw(14) = gf17->Add(rx.Raw(14), gf17->SetElementValue(3));

    auto res = rs.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 3U);
  }

  // Binary extension field GF(2^4) Extended RS(16, 12) -> corrects up to 2
  // errors including at alpha = 0
  {
    auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
    auto ext_rs = CreateExtendedReedSolomonCode(gf16, 16, 12);
    auto msg =
        GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf16, {12}, rng);
    auto cw = ext_rs.Encode(msg);
    EXPECT_TRUE(ext_rs.Contains(cw));

    auto rx = cw;
    rx.Raw(0) = gf16->Add(rx.Raw(0), gf16->SetElementValue(9));
    rx.Raw(11) = gf16->Add(rx.Raw(11), gf16->SetElementValue(4));

    auto res = ext_rs.DecodeWithResult(rx, "Sugiyama");
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 2U);
  }

  // General extension field GF(3^2) GRS(9, 5) with alpha_0 = 0
  {
    auto gf9 = std::make_shared<GaloisFieldExtension<uint8_t>>(
        std::pair<uint64_t, uint64_t>{3, 2});
    auto ext_rs = CreateExtendedReedSolomonCode(gf9, 9, 5);
    auto msg = GFArray<GaloisFieldExtension<uint8_t>>::Random(gf9, {5}, rng);
    auto cw = ext_rs.Encode(msg);
    EXPECT_TRUE(ext_rs.Contains(cw));

    auto rx = cw;
    rx.Raw(0) = gf9->Add(rx.Raw(0), gf9->FromInt(4));
    rx.Raw(6) = gf9->Add(rx.Raw(6), gf9->FromInt(7));

    auto res = ext_rs.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 2U);
  }
}

TEST(CodingTest, EncoderDecoderRegistrationAndUtilities) {
  auto gf8 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);
  auto rs = CreateReedSolomonCode(gf8, 7, 4);

  EXPECT_EQ(rs.DefaultEncoder(), "GRS");
  EXPECT_EQ(rs.DefaultDecoder(), "BerlekampMassey");
  EXPECT_NE(rs.GetEncoder("GRS"), nullptr);
  EXPECT_NE(rs.GetEncoder("GRSSystematic"), nullptr);
  EXPECT_NE(rs.GetEncoder("GeneratorMatrix"), nullptr);
  EXPECT_NE(rs.GetDecoder("BerlekampMassey"), nullptr);
  EXPECT_NE(rs.GetDecoder("Sugiyama"), nullptr);
  EXPECT_NE(rs.GetDecoder("PGZ"), nullptr);
  EXPECT_NE(rs.GetDecoder("GRS"), nullptr);

  // Distance & Weight
  auto w1 = MakeWord(gf8, {1, 2, 3, 4, 5, 6, 7});
  auto w2 = MakeWord(gf8, {1, 2, 3, 0, 0, 6, 7});
  EXPECT_EQ(rs.Distance(w1, w2), 2U);
  auto w3 = MakeWord(gf8, {0, 2, 0, 4, 0, 6, 0});
  EXPECT_EQ(rs.Weight(w3), 3U);

  // Polynomial encode / decode round-trip
  auto poly = PolynomialDense<GaloisFieldBinaryExtension<uint8_t>>::FromRaw(
      {1, 2, 3, 4}, gf8);
  auto cw = rs.EncodePolynomial(poly);
  EXPECT_TRUE(rs.Contains(cw));
  auto decoded_poly = rs.DecodeToPolynomial(cw);
  EXPECT_EQ(decoded_poly, poly);

  // Dual code
  auto dual = rs.DualCode();
  EXPECT_EQ(dual.Length(), 7U);
  EXPECT_EQ(dual.Dimension(), 3U);
  EXPECT_TRUE(IsMDSCode(dual));
}

TEST(CodingTest, CyclicCodeFactory) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  // g(x) = 1 + x + x^3 divides x^7 - 1 => [7, 4, 3] cyclic Hamming code.
  auto gen_poly =
      PolynomialDense<GaloisFieldBinary>::FromRaw({1, 1, 0, 1}, gf2);
  auto cyclic = CreateCyclicCode(gf2, 7, gen_poly);

  EXPECT_EQ(cyclic.Length(), 7U);
  EXPECT_EQ(cyclic.Dimension(), 4U);
  EXPECT_EQ(cyclic.GeneratorPolynomial().Degree(), 3);

  auto msg = MakeWord(gf2, {1, 0, 1, 1});
  auto cw = cyclic.Encode(msg);
  EXPECT_TRUE(cyclic.Contains(cw));
  EXPECT_EQ(cyclic.Unencode(cw), msg);

  auto shifted = cyclic.CyclicShift(cw, 1);
  EXPECT_TRUE(cyclic.Contains(shifted));
}

}  // namespace
