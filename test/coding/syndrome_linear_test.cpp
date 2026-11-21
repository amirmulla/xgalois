#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/coding/decoder/syndrome_decoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/linear_code.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace {

using namespace xg;
using namespace xg::coding;

template <typename Field>
void ExpectOrthogonal(const GFMatrix<Field> &G, const GFMatrix<Field> &H) {
  ASSERT_EQ(G.Shape(1), H.Shape(1));
  if (G.Shape(0) == 0 || H.Shape(0) == 0) {
    return;
  }
  auto prod = linalg::Dot(G, linalg::Transpose(H));
  ASSERT_EQ(prod.Shape(0), G.Shape(0));
  ASSERT_EQ(prod.Shape(1), H.Shape(0));
  auto field = G.Field();
  for (const auto &v : prod.RawSpan()) {
    EXPECT_TRUE(field->IsZero(v));
  }
}

template <typename Field>
GFArray<Field> MakeArray1D(const std::shared_ptr<Field> &field,
                           const std::vector<typename Field::element_type> &v) {
  using RawArray = typename GFArray<Field>::raw_array_type;
  using ShapeType = typename GFArray<Field>::shape_type;
  RawArray raw = RawArray::from_shape(ShapeType{v.size()});
  for (std::size_t i = 0; i < v.size(); ++i) {
    raw(i) = field->SetElementValue(v[i]);
  }
  return GFArray<Field>::FromRaw(field, std::move(raw));
}

template <typename Field>
GFMatrix<Field> MakeMatrix2D(
    const std::shared_ptr<Field> &field,
    const std::vector<std::vector<typename Field::element_type>> &rows) {
  using RawArray = typename GFMatrix<Field>::raw_array_type;
  using ShapeType = typename GFMatrix<Field>::shape_type;
  const std::size_t r = rows.size();
  const std::size_t c = r == 0 ? 0 : rows[0].size();
  RawArray raw = RawArray::from_shape(ShapeType{r, c});
  for (std::size_t i = 0; i < r; ++i) {
    for (std::size_t j = 0; j < c; ++j) {
      raw(i, j) = field->SetElementValue(rows[i][j]);
    }
  }
  return GFMatrix<Field>::FromRaw(field, std::move(raw));
}

TEST(LinearCodeTest, FromGeneratorAndParityCheckMatrices) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  // Non-systematic [7, 4, 3] Hamming generator matrix (scrambled columns/rows).
  auto G = MakeMatrix2D(gf2, {
                                 {1, 1, 0, 1, 0, 0, 0},
                                 {0, 1, 1, 0, 1, 0, 0},
                                 {0, 0, 1, 1, 0, 1, 0},
                                 {0, 0, 0, 1, 1, 0, 1},
                             });

  auto code = LinearCode<GaloisFieldBinary>::FromGeneratorMatrix(gf2, G);
  EXPECT_EQ(code.Length(), 7u);
  EXPECT_EQ(code.Dimension(), 4u);
  EXPECT_EQ(code.Redundancy(), 3u);
  EXPECT_EQ(code.MinimumDistance(), 3u);
  EXPECT_FALSE(code.IsSystematic());
  ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

  // Reconstruct from H (with an added redundant row = row 0 + row 1).
  const auto &H = code.ParityCheckMatrix();
  auto H_red = linalg::VStack(
      H,
      MakeMatrix2D(gf2, {{static_cast<uint8_t>(H.Raw(0, 0) ^ H.Raw(1, 0)),
                          static_cast<uint8_t>(H.Raw(0, 1) ^ H.Raw(1, 1)),
                          static_cast<uint8_t>(H.Raw(0, 2) ^ H.Raw(1, 2)),
                          static_cast<uint8_t>(H.Raw(0, 3) ^ H.Raw(1, 3)),
                          static_cast<uint8_t>(H.Raw(0, 4) ^ H.Raw(1, 4)),
                          static_cast<uint8_t>(H.Raw(0, 5) ^ H.Raw(1, 5)),
                          static_cast<uint8_t>(H.Raw(0, 6) ^ H.Raw(1, 6))}}));
  auto code2 = LinearCode<GaloisFieldBinary>::FromParityCheckMatrix(gf2, H_red);
  EXPECT_EQ(code2.Length(), 7u);
  EXPECT_EQ(code2.Dimension(), 4u);
  EXPECT_EQ(code2.ParityCheckMatrix().Shape(0), 3u);
  ExpectOrthogonal(code2.GeneratorMatrix(), code2.ParityCheckMatrix());

  // Every codeword of code is also in code2.
  for (const auto &cw : code.GetCodewords()) {
    EXPECT_TRUE(code2.Contains(cw));
  }
}

TEST(LinearCodeTest, DualPuncturedShortened) {
  auto code = ExtendedHammingCode(3);  // [8, 4, 4] self-dual
  ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

  auto dual = code.DualCode();
  EXPECT_EQ(dual.Length(), 8u);
  EXPECT_EQ(dual.Dimension(), 4u);
  EXPECT_EQ(dual.MinimumDistance(), 4u);
  ExpectOrthogonal(dual.GeneratorMatrix(), dual.ParityCheckMatrix());
  ExpectOrthogonal(code.GeneratorMatrix(), dual.GeneratorMatrix());

  // Puncturing [8, 4, 4] at the overall parity bit yields [7, 4, 3].
  auto punctured = code.Punctured(7);
  EXPECT_EQ(punctured.Length(), 7u);
  EXPECT_EQ(punctured.Dimension(), 4u);
  EXPECT_EQ(punctured.MinimumDistance(), 3u);
  ExpectOrthogonal(punctured.GeneratorMatrix(), punctured.ParityCheckMatrix());

  // Shortening [8, 4, 4] at coordinate 0 yields [7, 3, 4].
  auto shortened = code.Shortened(0);
  EXPECT_EQ(shortened.Length(), 7u);
  EXPECT_EQ(shortened.Dimension(), 3u);
  EXPECT_EQ(shortened.MinimumDistance(), 4u);
  ExpectOrthogonal(shortened.GeneratorMatrix(), shortened.ParityCheckMatrix());

  // Verify every codeword of `code` with c[0] == 0 maps into `shortened`.
  for (const auto &cw : code.GetCodewords()) {
    if (cw.Raw(0) == 0) {
      std::vector<uint8_t> dropped(7);
      for (std::size_t i = 0; i < 7; ++i) {
        dropped[i] = cw.Raw(i + 1);
      }
      EXPECT_TRUE(shortened.Contains(MakeArray1D(code.Field(), dropped)));
    }
  }
}

TEST(GeneratorMatrixEncoderTest, NonSystematicAcrossAllFieldTypes) {
  std::mt19937_64 rng(42);

  // 1. GF(2) non-systematic G
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    auto G = MakeMatrix2D(gf2, {
                                   {0, 1, 1, 1, 0, 0, 1},
                                   {1, 0, 1, 0, 1, 1, 0},
                                   {1, 1, 0, 1, 1, 0, 1},
                               });
    auto code = LinearCode<GaloisFieldBinary>::FromGeneratorMatrix(gf2, G);
    EXPECT_FALSE(code.IsSystematic());
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

    auto zero_msg = GFArray<GaloisFieldBinary>::Zeros(gf2, {3});
    auto zero_cw = code.Encode(zero_msg);
    EXPECT_TRUE(code.Contains(zero_cw));
    EXPECT_EQ(code.Unencode(zero_cw), zero_msg);

    for (int trial = 0; trial < 20; ++trial) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {3}, rng);
      auto cw = code.Encode(msg);
      EXPECT_TRUE(code.Contains(cw));
      EXPECT_EQ(code.Unencode(cw), msg);
    }

    // Invalid codeword must throw on Unencode.
    auto bad = zero_cw;
    bad.Raw(0) = 1;
    if (!code.Contains(bad)) {
      EXPECT_THROW(code.Unencode(bad), std::invalid_argument);
    }
  }

  // 2. GF(7) prime field non-systematic G
  {
    auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    auto G = MakeMatrix2D(gf7, {
                                   {2, 3, 1, 4, 5, 6},
                                   {0, 5, 2, 1, 3, 2},
                                   {4, 1, 6, 0, 2, 5},
                               });
    auto code =
        LinearCode<GaloisFieldPrime<uint8_t>>::FromGeneratorMatrix(gf7, G);
    EXPECT_FALSE(code.IsSystematic());
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

    auto zero_msg = GFArray<GaloisFieldPrime<uint8_t>>::Zeros(gf7, {3});
    EXPECT_EQ(code.Unencode(code.Encode(zero_msg)), zero_msg);

    for (int trial = 0; trial < 25; ++trial) {
      auto msg = GFArray<GaloisFieldPrime<uint8_t>>::Random(gf7, {3}, rng);
      auto cw = code.Encode(msg);
      EXPECT_TRUE(code.Contains(cw));
      EXPECT_EQ(code.Unencode(cw), msg);
    }
  }

  // 3. GF(2^4) binary extension field (standard and Zech log!)
  {
    auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
    auto G = MakeMatrix2D(gf16, {
                                    {3, 7, 1, 9, 12, 5},
                                    {0, 4, 11, 2, 6, 15},
                                    {5, 13, 8, 3, 1, 10},
                                });
    auto code = LinearCode<GF2X<uint8_t>>::FromGeneratorMatrix(gf16, G);
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

    auto zero_msg = GFArray<GF2X<uint8_t>>::Zeros(gf16, {3});
    EXPECT_EQ(code.Unencode(code.Encode(zero_msg)), zero_msg);

    for (int trial = 0; trial < 25; ++trial) {
      auto msg = GFArray<GF2X<uint8_t>>::Random(gf16, {3}, rng);
      auto cw = code.Encode(msg);
      EXPECT_TRUE(code.Contains(cw));
      EXPECT_EQ(code.Unencode(cw), msg);
    }
  }
  {
    auto gf16z = std::make_shared<GF2XZECH<uint8_t>>(4);
    auto G = MakeMatrix2D(gf16z, {
                                     {3, 7, 1, 9, 12},
                                     {2, 0, 11, 4, 6},
                                 });
    auto code = LinearCode<GF2XZECH<uint8_t>>::FromGeneratorMatrix(gf16z, G);
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

    auto zero_msg = GFArray<GF2XZECH<uint8_t>>::Zeros(gf16z, {2});
    EXPECT_EQ(code.Unencode(code.Encode(zero_msg)), zero_msg);

    for (int trial = 0; trial < 20; ++trial) {
      auto msg = GFArray<GF2XZECH<uint8_t>>::Random(gf16z, {2}, rng);
      auto cw = code.Encode(msg);
      EXPECT_TRUE(code.Contains(cw));
      EXPECT_EQ(code.Unencode(cw), msg);
    }
  }

  // 4. GF(3^2) general extension field (ExtensionElement)
  {
    using Ext9 = GaloisFieldExtension<uint8_t>;
    auto gf9 = std::make_shared<Ext9>(std::make_pair(3ULL, 2ULL));
    auto e = [&](uint64_t v) { return gf9->FromInt(v); };
    auto G = MakeMatrix2D(gf9, {
                                   {e(2), e(5), e(1), e(7), e(4)},
                                   {e(0), e(3), e(8), e(2), e(6)},
                               });
    auto code = LinearCode<Ext9>::FromGeneratorMatrix(gf9, G);
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

    auto zero_msg = GFArray<Ext9>::Zeros(gf9, {2});
    EXPECT_EQ(code.Unencode(code.Encode(zero_msg)), zero_msg);

    for (int trial = 0; trial < 20; ++trial) {
      auto msg = GFArray<Ext9>::Random(gf9, {2}, rng);
      auto cw = code.Encode(msg);
      EXPECT_TRUE(code.Contains(cw));
      EXPECT_EQ(code.Unencode(cw), msg);
    }

    // Also test SyndromeDecoder on GF(3^2) (ExtensionElement hash!).
    SyndromeDecoder<Ext9> dec(code, 1);
    auto msg = GFArray<Ext9>::Random(gf9, {2}, rng);
    auto cw = code.Encode(msg);
    auto rx = cw;
    rx.Raw(3) = gf9->Add(rx.Raw(3), gf9->FromInt(5));
    auto res = dec.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 1u);
  }
}

TEST(SyndromeDecoderTest, EncoderDecoderCachingAndCopySafety) {
  auto code = HammingCode(3);
  auto enc1 = code.GetEncoder();
  auto enc2 = code.GetEncoder();
  EXPECT_EQ(enc1.get(), enc2.get());

  auto dec1 = code.GetDecoder();
  auto dec2 = code.GetDecoder();
  EXPECT_EQ(dec1.get(), dec2.get());

  // Copy construction must not reuse encoder/decoder pointers bound to
  // original.
  auto copy = code;
  auto copy_enc = copy.GetEncoder();
  EXPECT_NE(copy_enc.get(), enc1.get());
  EXPECT_EQ(&copy_enc->Code(), &copy);
}

TEST(SyndromeDecoderTest, MemoryCapThrowsClearly) {
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto code = HammingCode(gf7, 3);  // n = (343-1)/6 = 57
  // Requesting max_errors = 2 with a tiny cap must throw std::runtime_error.
  EXPECT_THROW((SyndromeDecoder<GaloisFieldPrime<uint8_t>>(code, 2, 50)),
               std::runtime_error);
}

TEST(ClassicCodesTest, BinaryHamming7And15) {
  auto ham7 = HammingCode(3);
  EXPECT_EQ(ham7.Length(), 7u);
  EXPECT_EQ(ham7.Dimension(), 4u);
  EXPECT_EQ(ham7.MinimumDistance(), 3u);
  ExpectOrthogonal(ham7.GeneratorMatrix(), ham7.ParityCheckMatrix());

  auto msg = MakeArray1D(ham7.Field(), {1, 0, 1, 1});
  auto cw = ham7.Encode(msg);
  EXPECT_TRUE(ham7.Contains(cw));
  EXPECT_EQ(ham7.Unencode(cw), msg);

  // Zero-error decode
  auto res0 = ham7.DecodeWithResult(cw);
  ASSERT_TRUE(res0.Ok());
  EXPECT_EQ(res0.codeword, cw);
  EXPECT_EQ(res0.message, msg);
  EXPECT_EQ(res0.num_errors, 0u);

  // Exhaustive single-bit error correction
  for (std::size_t pos = 0; pos < 7; ++pos) {
    auto rx = cw;
    rx.Raw(pos) ^= 1u;
    auto res = ham7.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 1u);
  }

  auto ham15 = HammingCode(4);
  EXPECT_EQ(ham15.Length(), 15u);
  EXPECT_EQ(ham15.Dimension(), 11u);
  EXPECT_EQ(ham15.MinimumDistance(), 3u);
  ExpectOrthogonal(ham15.GeneratorMatrix(), ham15.ParityCheckMatrix());

  std::mt19937_64 rng(123);
  auto msg15 = GFArray<GaloisFieldBinary>::Random(ham15.Field(), {11}, rng);
  auto cw15 = ham15.Encode(msg15);
  for (std::size_t pos = 0; pos < 15; ++pos) {
    auto rx = cw15;
    rx.Raw(pos) ^= 1u;
    auto res = ham15.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw15);
    EXPECT_EQ(res.message, msg15);
    EXPECT_EQ(res.num_errors, 1u);
  }
}

TEST(ClassicCodesTest, QAryHammingOverGF3) {
  auto gf3 = std::make_shared<GaloisFieldPrime<uint8_t>>(3);
  auto ham = HammingCode(gf3, 3);  // [(27-1)/2 = 13, 10, 3]
  EXPECT_EQ(ham.Length(), 13u);
  EXPECT_EQ(ham.Dimension(), 10u);
  EXPECT_EQ(ham.MinimumDistance(), 3u);
  ExpectOrthogonal(ham.GeneratorMatrix(), ham.ParityCheckMatrix());

  std::mt19937_64 rng(99);
  auto msg = GFArray<GaloisFieldPrime<uint8_t>>::Random(gf3, {10}, rng);
  auto cw = ham.Encode(msg);
  EXPECT_TRUE(ham.Contains(cw));

  // Exhaustive single-symbol error correction (all 13 positions x 2 nonzero
  // vals)
  for (std::size_t pos = 0; pos < 13; ++pos) {
    for (uint8_t err : {uint8_t{1}, uint8_t{2}}) {
      auto rx = cw;
      rx.Raw(pos) = gf3->Add(rx.Raw(pos), err);
      auto res = ham.DecodeWithResult(rx);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 1u);
    }
  }
}

TEST(ClassicCodesTest, ExtendedHammingDetectsDoubleErrors) {
  auto ext8 = ExtendedHammingCode(3);
  EXPECT_EQ(ext8.Length(), 8u);
  EXPECT_EQ(ext8.Dimension(), 4u);
  EXPECT_EQ(ext8.MinimumDistance(), 4u);
  ExpectOrthogonal(ext8.GeneratorMatrix(), ext8.ParityCheckMatrix());

  auto msg = MakeArray1D(ext8.Field(), {1, 1, 0, 1});
  auto cw = ext8.Encode(msg);

  // Corrects all 1-bit errors
  for (std::size_t i = 0; i < 8; ++i) {
    auto rx = cw;
    rx.Raw(i) ^= 1u;
    auto res = ext8.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 1u);
  }

  // Explicitly detects ALL 2-bit errors as uncorrectable (never silently
  // corrupt)
  for (std::size_t i = 0; i < 8; ++i) {
    for (std::size_t j = i + 1; j < 8; ++j) {
      auto rx = cw;
      rx.Raw(i) ^= 1u;
      rx.Raw(j) ^= 1u;
      auto res = ext8.DecodeWithResult(rx);
      EXPECT_FALSE(res.Ok());
      EXPECT_FALSE(res.error_message.empty());
    }
  }
}

TEST(ClassicCodesTest, BinaryGolay23Exhaustive3BitCorrection) {
  auto golay23 = BinaryGolayCode();
  EXPECT_EQ(golay23.Length(), 23u);
  EXPECT_EQ(golay23.Dimension(), 12u);
  EXPECT_EQ(golay23.MinimumDistance(), 7u);
  ExpectOrthogonal(golay23.GeneratorMatrix(), golay23.ParityCheckMatrix());

  std::mt19937_64 rng(2025);
  auto msg = GFArray<GaloisFieldBinary>::Random(golay23.Field(), {12}, rng);
  auto cw = golay23.Encode(msg);
  EXPECT_TRUE(golay23.Contains(cw));
  EXPECT_EQ(golay23.Unencode(cw), msg);

  // Exhaustively verify ALL weight 1, 2, and 3 error patterns (2047 total!)
  for (std::size_t i = 0; i < 23; ++i) {
    auto rx1 = cw;
    rx1.Raw(i) ^= 1u;
    auto res1 = golay23.DecodeWithResult(rx1);
    ASSERT_TRUE(res1.Ok());
    EXPECT_EQ(res1.codeword, cw);
    EXPECT_EQ(res1.message, msg);
    EXPECT_EQ(res1.num_errors, 1u);

    for (std::size_t j = i + 1; j < 23; ++j) {
      auto rx2 = rx1;
      rx2.Raw(j) ^= 1u;
      auto res2 = golay23.DecodeWithResult(rx2);
      ASSERT_TRUE(res2.Ok());
      EXPECT_EQ(res2.codeword, cw);
      EXPECT_EQ(res2.num_errors, 2u);

      for (std::size_t l = j + 1; l < 23; ++l) {
        auto rx3 = rx2;
        rx3.Raw(l) ^= 1u;
        auto res3 = golay23.DecodeWithResult(rx3);
        ASSERT_TRUE(res3.Ok());
        EXPECT_EQ(res3.codeword, cw);
        EXPECT_EQ(res3.num_errors, 3u);
      }
    }
  }
}

TEST(ClassicCodesTest, ExtendedBinaryGolay24SelfDualAnd4BitDetection) {
  auto golay24 = ExtendedBinaryGolayCode();
  EXPECT_EQ(golay24.Length(), 24u);
  EXPECT_EQ(golay24.Dimension(), 12u);
  EXPECT_EQ(golay24.MinimumDistance(), 8u);
  ExpectOrthogonal(golay24.GeneratorMatrix(), golay24.ParityCheckMatrix());
  // Self-dual: G * G^T = 0.
  ExpectOrthogonal(golay24.GeneratorMatrix(), golay24.GeneratorMatrix());

  auto msg = GFArray<GaloisFieldBinary>::Ones(golay24.Field(), {12});
  auto cw = golay24.Encode(msg);
  EXPECT_EQ(golay24.Weight(cw) % 4, 0u);  // Doubly-even

  // 3-bit error corrected
  auto rx3 = cw;
  rx3.Raw(0) ^= 1u;
  rx3.Raw(11) ^= 1u;
  rx3.Raw(23) ^= 1u;
  auto res3 = golay24.DecodeWithResult(rx3);
  ASSERT_TRUE(res3.Ok());
  EXPECT_EQ(res3.codeword, cw);
  EXPECT_EQ(res3.message, msg);
  EXPECT_EQ(res3.num_errors, 3u);

  // 4-bit error detected as uncorrectable
  auto rx4 = rx3;
  rx4.Raw(5) ^= 1u;
  auto res4 = golay24.DecodeWithResult(rx4);
  EXPECT_FALSE(res4.Ok());
}

TEST(ClassicCodesTest, TernaryGolay11Exhaustive2SymbolCorrection) {
  auto tg = TernaryGolayCode();
  auto gf3 = tg.Field();
  EXPECT_EQ(tg.Length(), 11u);
  EXPECT_EQ(tg.Dimension(), 6u);
  EXPECT_EQ(tg.MinimumDistance(), 5u);
  ExpectOrthogonal(tg.GeneratorMatrix(), tg.ParityCheckMatrix());

  auto msg = MakeArray1D(gf3, {1, 2, 0, 2, 1, 1});
  auto cw = tg.Encode(msg);
  EXPECT_TRUE(tg.Contains(cw));
  EXPECT_EQ(tg.Unencode(cw), msg);

  // Exhaustively verify all weight-1 and weight-2 errors over GF(3)
  for (std::size_t i = 0; i < 11; ++i) {
    for (uint8_t e1 : {uint8_t{1}, uint8_t{2}}) {
      auto rx1 = cw;
      rx1.Raw(i) = gf3->Add(rx1.Raw(i), e1);
      auto res1 = tg.DecodeWithResult(rx1);
      ASSERT_TRUE(res1.Ok());
      EXPECT_EQ(res1.codeword, cw);
      EXPECT_EQ(res1.message, msg);
      EXPECT_EQ(res1.num_errors, 1u);

      for (std::size_t j = i + 1; j < 11; ++j) {
        for (uint8_t e2 : {uint8_t{1}, uint8_t{2}}) {
          auto rx2 = rx1;
          rx2.Raw(j) = gf3->Add(rx2.Raw(j), e2);
          auto res2 = tg.DecodeWithResult(rx2);
          ASSERT_TRUE(res2.Ok());
          EXPECT_EQ(res2.codeword, cw);
          EXPECT_EQ(res2.num_errors, 2u);
        }
      }
    }
  }
}

TEST(ClassicCodesTest, ReedMullerCodes) {
  auto rm13 = ReedMullerCode(1, 3);
  EXPECT_EQ(rm13.Length(), 8u);
  EXPECT_EQ(rm13.Dimension(), 4u);
  EXPECT_EQ(rm13.MinimumDistance(), 4u);
  ExpectOrthogonal(rm13.GeneratorMatrix(), rm13.ParityCheckMatrix());

  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto rm14 = ReedMullerCode(gf2, 1, 4);
  EXPECT_EQ(rm14.Length(), 16u);
  EXPECT_EQ(rm14.Dimension(), 5u);
  EXPECT_EQ(rm14.MinimumDistance(), 8u);
  ExpectOrthogonal(rm14.GeneratorMatrix(), rm14.ParityCheckMatrix());

  auto rm24 = ReedMullerCode(gf2, 2, 4);
  EXPECT_EQ(rm24.Length(), 16u);
  EXPECT_EQ(rm24.Dimension(), 11u);
  EXPECT_EQ(rm24.MinimumDistance(), 4u);
  ExpectOrthogonal(rm24.GeneratorMatrix(), rm24.ParityCheckMatrix());

  // Duality theorem: RM(1, 4)^perp == RM(2, 4) (G_{1,4} * G_{2,4}^T == 0).
  ExpectOrthogonal(rm14.GeneratorMatrix(), rm24.GeneratorMatrix());

  // RM(1, 4) has d = 8, so SyndromeDecoder corrects up to 3 errors.
  auto msg = MakeArray1D(rm14.Field(), {1, 0, 1, 1, 0});
  auto cw = rm14.Encode(msg);
  auto rx = cw;
  rx.Raw(1) ^= 1u;
  rx.Raw(7) ^= 1u;
  rx.Raw(14) ^= 1u;
  auto res = rm14.DecodeWithResult(rx);
  ASSERT_TRUE(res.Ok());
  EXPECT_EQ(res.codeword, cw);
  EXPECT_EQ(res.message, msg);
  EXPECT_EQ(res.num_errors, 3u);
}

TEST(MetricTest, LeeWeightAndDistanceOverPrimeFields) {
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto G = MakeMatrix2D(gf7, {{1, 1, 1, 1, 1}});
  auto code = LinearCode<GaloisFieldPrime<uint8_t>>::FromGeneratorMatrix(
      gf7, G, Metric::LEE);

  // In GF(7), Lee weights are: 0->0, 1->1, 2->2, 3->3, 4->3, 5->2, 6->1.
  auto w1 = MakeArray1D(gf7, {0, 1, 2, 3, 6});
  EXPECT_EQ(code.Weight(w1), 0u + 1u + 2u + 3u + 1u);

  auto w2 = MakeArray1D(gf7, {1, 6, 4, 3, 2});
  // w1 - w2 mod 7 = {6, 2, 5, 0, 4} -> Lee weights {1, 2, 2, 0, 3} -> sum 8.
  EXPECT_EQ(code.Distance(w1, w2), 8u);
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
