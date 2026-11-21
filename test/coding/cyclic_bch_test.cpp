#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/coding/bch.hpp"
#include "xgalois/coding/cyclic_code.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"

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
PolynomialDense<Field> MakePoly(
    const std::shared_ptr<Field> &field,
    const std::vector<typename Field::element_type> &coeffs) {
  std::vector<typename Field::element_type> raw;
  raw.reserve(coeffs.size());
  for (const auto &c : coeffs) {
    raw.push_back(field->SetElementValue(c));
  }
  return PolynomialDense<Field>::FromRaw(std::move(raw), field);
}

TEST(CyclicCodeTest, ValidationRejectsInvalidGeneratorPolynomials) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // Zero polynomial
  auto zero_p = PolynomialDense<GaloisFieldBinary>::Zero(gf2);
  EXPECT_THROW((CyclicCode<GaloisFieldBinary>(gf2, 7, zero_p)),
               std::invalid_argument);

  // Constant polynomial
  auto one_p = PolynomialDense<GaloisFieldBinary>::One(gf2);
  EXPECT_THROW((CyclicCode<GaloisFieldBinary>(gf2, 7, one_p)),
               std::invalid_argument);

  // Does not divide x^7 - 1 (x^2 + x + 1 divides x^3 - 1, not x^7 - 1)
  auto bad_g = MakePoly(gf2, {1, 1, 1});
  EXPECT_THROW((CyclicCode<GaloisFieldBinary>(gf2, 7, bad_g)),
               std::invalid_argument);

  // Degree >= length
  auto too_big = MakePoly(gf2, {1, 0, 0, 0, 0, 0, 0, 1});  // x^7 + 1
  EXPECT_THROW((CyclicCode<GaloisFieldBinary>(gf2, 7, too_big)),
               std::invalid_argument);
}

TEST(CyclicCodeTest, BinaryHamming7AndDualAndShifts) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  // g(x) = 1 + x + x^3 divides x^7 - 1
  auto g = MakePoly(gf2, {1, 1, 0, 1});
  CyclicCode<GaloisFieldBinary> code(gf2, 7, g);

  EXPECT_EQ(code.Length(), 7u);
  EXPECT_EQ(code.Dimension(), 4u);
  EXPECT_EQ(code.MinimumDistance(), 3u);
  EXPECT_TRUE(code.IsBCH());
  ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

  // FromParityCheckPolynomial reproduces exact same code.
  auto from_h = CyclicCode<GaloisFieldBinary>::FromParityCheckPolynomial(
      gf2, 7, code.ParityCheckPolynomial());
  EXPECT_EQ(from_h.GeneratorPolynomial(), code.GeneratorPolynomial());

  // Dual code orthogonality & involution (C^perp)^perp == C.
  auto dual = code.DualCode();
  EXPECT_EQ(dual.Length(), 7u);
  EXPECT_EQ(dual.Dimension(), 3u);
  EXPECT_EQ(dual.MinimumDistance(), 4u);
  ExpectOrthogonal(dual.GeneratorMatrix(), dual.ParityCheckMatrix());
  ExpectOrthogonal(code.GeneratorMatrix(), dual.GeneratorMatrix());
  EXPECT_EQ(dual.DualCode().GeneratorPolynomial(), code.GeneratorPolynomial());

  // CyclicEncoder: non-systematic and systematic round-trip (including 0 msg)
  auto zero_msg = GFArray<GaloisFieldBinary>::Zeros(gf2, {4});
  auto cw0 = code.Encode(zero_msg, "Cyclic");
  EXPECT_TRUE(code.Contains(cw0));
  EXPECT_EQ(code.Unencode(cw0, "Cyclic"), zero_msg);

  auto cw0_sys = code.Encode(zero_msg, "CyclicSystematic");
  EXPECT_TRUE(code.Contains(cw0_sys));
  EXPECT_EQ(code.Unencode(cw0_sys, "CyclicSystematic"), zero_msg);

  auto msg = MakeArray1D(gf2, {1, 0, 1, 1});
  auto cw = code.Encode(msg, "Cyclic");
  EXPECT_TRUE(code.Contains(cw));
  EXPECT_EQ(code.Unencode(cw, "Cyclic"), msg);

  auto cw_sys = code.Encode(msg, "CyclicSystematic");
  EXPECT_TRUE(code.Contains(cw_sys));
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(cw_sys.Raw(3 + i), msg.Raw(i));
  }
  EXPECT_EQ(code.Unencode(cw_sys, "CyclicSystematic"), msg);

  // Non-codeword throws on Unencode in both modes
  auto corrupted = cw;
  corrupted.Raw(2) ^= 1u;
  EXPECT_THROW(code.Unencode(corrupted, "Cyclic"), std::invalid_argument);
  EXPECT_THROW(code.Unencode(corrupted, "CyclicSystematic"),
               std::invalid_argument);

  // CyclicShift positive, negative, large shifts preserve codeword membership.
  for (int s : {-15, -7, -3, -1, 0, 1, 3, 7, 15}) {
    auto shifted = code.CyclicShift(cw, s);
    EXPECT_TRUE(code.Contains(shifted));
    EXPECT_EQ(code.CyclicShift(shifted, -s), cw);
  }

  // CyclicSyndromeDecoder corrects all single-bit errors.
  for (std::size_t pos = 0; pos < 7; ++pos) {
    auto rx = cw;
    rx.Raw(pos) ^= 1u;
    auto res = code.DecodeWithResult(rx, "CyclicSyndrome");
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 1u);
  }
}

TEST(CyclicCodeTest, NonBchCyclicCodeIsDetected) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  // Over GF(2), n = 15: g(x) = (x + 1)(x^2 + x + 1) = x^3 + 1 divides x^15 + 1.
  // Its defining set is {0, 5, 10}, which cannot be generated by any
  // consecutive run for any primitive 15-th root of unity.
  auto g = MakePoly(gf2, {1, 0, 0, 1});
  CyclicCode<GaloisFieldBinary> code(gf2, 15, g);
  ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());
  EXPECT_FALSE(code.IsBCH());
}

TEST(CyclicCodeTest, PrimeFieldCyclicCodesOverGF5AndGF7) {
  // GF(5), n = 4: x^4 - 1 = (x - 1)(x - 2)(x - 3)(x - 4)
  // Take g(x) = (x - 1)(x - 2) = x^2 - 3x + 2 = x^2 + 2x + 2 (mod 5).
  {
    auto gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);
    auto g = MakePoly(gf5, {2, 2, 1});
    CyclicCode<GaloisFieldPrime<uint8_t>> code(gf5, 4, g);
    EXPECT_EQ(code.Length(), 4u);
    EXPECT_EQ(code.Dimension(), 2u);
    EXPECT_EQ(code.MinimumDistance(), 3u);
    EXPECT_TRUE(code.IsBCH());
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());
    ExpectOrthogonal(code.GeneratorMatrix(), code.DualCode().GeneratorMatrix());

    auto msg = MakeArray1D(gf5, {3, 4});
    auto cw = code.Encode(msg, "Cyclic");
    EXPECT_EQ(code.Unencode(cw, "Cyclic"), msg);

    auto cw_sys = code.Encode(msg, "CyclicSystematic");
    EXPECT_TRUE(code.Contains(cw_sys));
    EXPECT_EQ(code.Unencode(cw_sys, "CyclicSystematic"), msg);

    // Exhaustive single-symbol error correction via CyclicSyndromeDecoder
    for (std::size_t pos = 0; pos < 4; ++pos) {
      for (uint8_t err = 1; err < 5; ++err) {
        auto rx = cw;
        rx.Raw(pos) = gf5->Add(rx.Raw(pos), err);
        auto res = code.DecodeWithResult(rx, "CyclicSyndrome");
        ASSERT_TRUE(res.Ok());
        EXPECT_EQ(res.codeword, cw);
        EXPECT_EQ(res.message, msg);
        EXPECT_EQ(res.num_errors, 1u);
      }
    }
  }

  // GF(7), n = 6: primitive root 3; x^6 - 1 = prod_{k=0}^5 (x - 3^k).
  // Take g(x) = (x - 3)(x - 2)(x - 6)(x - 4) of degree 4 => [6, 2, 5] MDS!
  {
    auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    auto f1 = MakePoly(gf7, {4, 1});  // x - 3
    auto f2 = MakePoly(gf7, {5, 1});  // x - 2 (3^2 = 2)
    auto f3 = MakePoly(gf7, {1, 1});  // x - 6 (3^3 = 6)
    auto f4 = MakePoly(gf7, {3, 1});  // x - 4 (3^4 = 4)
    auto g = f1 * f2 * f3 * f4;

    CyclicCode<GaloisFieldPrime<uint8_t>> code(gf7, 6, g);
    EXPECT_EQ(code.Length(), 6u);
    EXPECT_EQ(code.Dimension(), 2u);
    EXPECT_EQ(code.MinimumDistance(), 5u);
    EXPECT_TRUE(code.IsBCH());
    ExpectOrthogonal(code.GeneratorMatrix(), code.ParityCheckMatrix());

    auto msg = MakeArray1D(gf7, {5, 2});
    auto cw = code.Encode(msg);
    EXPECT_TRUE(code.Contains(cw));

    // Correct a weight-2 error with CyclicSyndromeDecoder
    auto rx = cw;
    rx.Raw(1) = gf7->Add(rx.Raw(1), 4);
    rx.Raw(4) = gf7->Add(rx.Raw(4), 6);
    auto res = code.DecodeWithResult(rx, "CyclicSyndrome");
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 2u);
  }
}

TEST(BCHCodeTest, BinaryBCH15_7_5Exhaustive2BitErrors) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);

  auto bch =
      BCHCode<GaloisFieldBinary, GF2X<uint8_t>>::FromErrorCorrectionCapability(
          gf2, gf16, 2);
  EXPECT_EQ(bch.Length(), 15u);
  EXPECT_EQ(bch.Dimension(), 7u);
  EXPECT_EQ(bch.DesignedDistance(), 5u);
  EXPECT_EQ(bch.MinimumDistance(), 5u);
  EXPECT_EQ(bch.ErrorCorrectionCapability(), 2u);
  EXPECT_TRUE(bch.IsBCH());
  ExpectOrthogonal(bch.GeneratorMatrix(), bch.ParityCheckMatrix());

  // All-zero message round-trip
  auto zero_msg = GFArray<GaloisFieldBinary>::Zeros(gf2, {7});
  auto zero_cw = bch.Encode(zero_msg);
  auto res0 = bch.DecodeWithResult(zero_cw);
  ASSERT_TRUE(res0.Ok());
  EXPECT_EQ(res0.codeword, zero_cw);
  EXPECT_EQ(res0.message, zero_msg);
  EXPECT_EQ(res0.num_errors, 0u);

  auto msg = MakeArray1D(gf2, {1, 0, 1, 1, 0, 1, 0});
  auto cw = bch.Encode(msg);
  EXPECT_TRUE(bch.Contains(cw));

  // Exhaustive weight-1 and weight-2 error correction (120 error patterns)
  for (std::size_t i = 0; i < 15; ++i) {
    auto rx1 = cw;
    rx1.Raw(i) ^= 1u;
    auto res1 = bch.DecodeWithResult(rx1);
    ASSERT_TRUE(res1.Ok());
    EXPECT_EQ(res1.codeword, cw);
    EXPECT_EQ(res1.message, msg);
    EXPECT_EQ(res1.num_errors, 1u);

    for (std::size_t j = i + 1; j < 15; ++j) {
      auto rx2 = rx1;
      rx2.Raw(j) ^= 1u;
      auto res2 = bch.DecodeWithResult(rx2);
      ASSERT_TRUE(res2.Ok());
      EXPECT_EQ(res2.codeword, cw);
      EXPECT_EQ(res2.message, msg);
      EXPECT_EQ(res2.num_errors, 2u);
    }
  }

  // Also cross-check CyclicSyndromeDecoder on the same [15, 7, 5] BCH code.
  auto rx_cyc = cw;
  rx_cyc.Raw(2) ^= 1u;
  rx_cyc.Raw(11) ^= 1u;
  auto res_cyc = bch.DecodeWithResult(rx_cyc, "CyclicSyndrome");
  ASSERT_TRUE(res_cyc.Ok());
  EXPECT_EQ(res_cyc.codeword, cw);
  EXPECT_EQ(res_cyc.message, msg);
  EXPECT_EQ(res_cyc.num_errors, 2u);
}

TEST(BCHCodeTest, BinaryBCH15_5_7Exhaustive3BitErrors) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);

  BCHCode<GaloisFieldBinary, GF2X<uint8_t>> bch(gf2, gf16, 7);
  EXPECT_EQ(bch.Length(), 15u);
  EXPECT_EQ(bch.Dimension(), 5u);
  EXPECT_EQ(bch.MinimumDistance(), 7u);
  ExpectOrthogonal(bch.GeneratorMatrix(), bch.ParityCheckMatrix());

  auto msg = MakeArray1D(gf2, {1, 1, 0, 1, 1});
  auto cw = bch.Encode(msg);

  // Exhaustive weight-1, weight-2, and weight-3 error correction (575 patterns)
  for (std::size_t i = 0; i < 15; ++i) {
    auto rx1 = cw;
    rx1.Raw(i) ^= 1u;
    auto res1 = bch.DecodeWithResult(rx1);
    ASSERT_TRUE(res1.Ok());
    EXPECT_EQ(res1.codeword, cw);
    EXPECT_EQ(res1.num_errors, 1u);

    for (std::size_t j = i + 1; j < 15; ++j) {
      auto rx2 = rx1;
      rx2.Raw(j) ^= 1u;
      auto res2 = bch.DecodeWithResult(rx2);
      ASSERT_TRUE(res2.Ok());
      EXPECT_EQ(res2.codeword, cw);
      EXPECT_EQ(res2.num_errors, 2u);

      for (std::size_t l = j + 1; l < 15; ++l) {
        auto rx3 = rx2;
        rx3.Raw(l) ^= 1u;
        auto res3 = bch.DecodeWithResult(rx3);
        ASSERT_TRUE(res3.Ok());
        EXPECT_EQ(res3.codeword, cw);
        EXPECT_EQ(res3.message, msg);
        EXPECT_EQ(res3.num_errors, 3u);
      }
    }
  }
}

TEST(BCHCodeTest, BinaryBCH31_16_7And63_45_7Random3BitAndUncorrectable) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  std::mt19937_64 rng(777);

  // BCH [31, 16, 7] (t = 3)
  {
    auto gf32 = std::make_shared<GF2X<uint8_t>>(5);
    BCHCode<GaloisFieldBinary, GF2X<uint8_t>> bch31(gf2, gf32, 7);
    EXPECT_EQ(bch31.Length(), 31u);
    EXPECT_EQ(bch31.Dimension(), 16u);
    EXPECT_EQ(bch31.MinimumDistance(), 7u);
    ExpectOrthogonal(bch31.GeneratorMatrix(), bch31.ParityCheckMatrix());

    std::vector<std::size_t> indices(31);
    for (std::size_t i = 0; i < 31; ++i) indices[i] = i;

    for (int trial = 0; trial < 50; ++trial) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {16}, rng);
      auto cw = bch31.Encode(msg, "CyclicSystematic");
      EXPECT_TRUE(bch31.Contains(cw));

      std::shuffle(indices.begin(), indices.end(), rng);
      const std::size_t weight = static_cast<std::size_t>(trial % 4);  // 0..3
      auto rx = cw;
      for (std::size_t e = 0; e < weight; ++e) {
        rx.Raw(indices[e]) ^= 1u;
      }
      auto res = bch31.DecodeWithResult(rx);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.num_errors, weight);
      EXPECT_EQ(bch31.Unencode(res.codeword, "CyclicSystematic"), msg);
    }

    // Verify uncorrectable errors (> 3 bits) are explicitly detected.
    auto cw = bch31.Encode(GFArray<GaloisFieldBinary>::Zeros(gf2, {16}));
    std::size_t detected_uncorrectable = 0;
    for (int trial = 0; trial < 40; ++trial) {
      std::shuffle(indices.begin(), indices.end(), rng);
      auto rx = cw;
      for (std::size_t e = 0; e < 5; ++e) {
        rx.Raw(indices[e]) ^= 1u;
      }
      auto res = bch31.DecodeWithResult(rx);
      if (!res.Ok()) {
        ++detected_uncorrectable;
      } else {
        // If it did decode (bounded-distance hit another codeword), result MUST
        // be a genuine codeword at distance <= 3 from rx!
        EXPECT_TRUE(bch31.Contains(res.codeword));
        EXPECT_LE(bch31.Distance(rx, res.codeword), 3u);
      }
    }
    EXPECT_GT(detected_uncorrectable, 20u);
  }

  // BCH [63, 45, 7] (t = 3)
  {
    auto gf64 = std::make_shared<GF2X<uint16_t>>(6);
    BCHCode<GaloisFieldBinary, GF2X<uint16_t>> bch63(gf2, gf64, 7);
    EXPECT_EQ(bch63.Length(), 63u);
    EXPECT_EQ(bch63.Dimension(), 45u);
    EXPECT_EQ(bch63.MinimumDistance(), 7u);
    ExpectOrthogonal(bch63.GeneratorMatrix(), bch63.ParityCheckMatrix());

    std::vector<std::size_t> indices(63);
    for (std::size_t i = 0; i < 63; ++i) indices[i] = i;

    for (int trial = 0; trial < 40; ++trial) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {45}, rng);
      auto cw = bch63.Encode(msg);
      EXPECT_TRUE(bch63.Contains(cw));

      std::shuffle(indices.begin(), indices.end(), rng);
      const std::size_t weight = static_cast<std::size_t>(trial % 4);  // 0..3
      auto rx = cw;
      for (std::size_t e = 0; e < weight; ++e) {
        rx.Raw(indices[e]) ^= 1u;
      }
      auto res = bch63.DecodeWithResult(rx);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, weight);
    }

    // Verify uncorrectable errors on BCH [63, 45, 7] are detected.
    auto cw = bch63.Encode(GFArray<GaloisFieldBinary>::Zeros(gf2, {45}));
    std::size_t detected = 0;
    for (int trial = 0; trial < 30; ++trial) {
      std::shuffle(indices.begin(), indices.end(), rng);
      auto rx = cw;
      for (std::size_t e = 0; e < 5; ++e) {
        rx.Raw(indices[e]) ^= 1u;
      }
      auto res = bch63.DecodeWithResult(rx);
      if (!res.Ok()) {
        ++detected;
      } else {
        EXPECT_TRUE(bch63.Contains(res.codeword));
        EXPECT_LE(bch63.Distance(rx, res.codeword), 3u);
      }
    }
    EXPECT_GT(detected, 20u);
  }
}

TEST(BCHCodeTest, NonBinaryBCHWithForney) {
  // Ternary BCH of length n = 8 over GF(3) with splitting field GF(3^2),
  // designed distance delta = 5 (t = 2), b = 1.
  using GF3 = GaloisFieldPrime<uint8_t>;
  using GF9 = GaloisFieldExtension<uint8_t>;
  auto gf3 = std::make_shared<GF3>(3);
  auto gf9 = std::make_shared<GF9>(std::make_pair(3ULL, 2ULL));

  BCHCode<GF3, GF9> tbch(gf3, gf9, 5, 1);
  EXPECT_EQ(tbch.Length(), 8u);
  EXPECT_EQ(tbch.ErrorCorrectionCapability(), 2u);
  ExpectOrthogonal(tbch.GeneratorMatrix(), tbch.ParityCheckMatrix());

  std::mt19937_64 rng(31415);
  auto msg = GFArray<GF3>::Random(gf3, {tbch.Dimension()}, rng);
  auto cw = tbch.Encode(msg);
  EXPECT_TRUE(tbch.Contains(cw));

  // Exhaustive weight-1 and weight-2 non-binary errors over GF(3).
  for (std::size_t i = 0; i < 8; ++i) {
    for (uint8_t e1 : {uint8_t{1}, uint8_t{2}}) {
      auto rx1 = cw;
      rx1.Raw(i) = gf3->Add(rx1.Raw(i), e1);
      auto res1 = tbch.DecodeWithResult(rx1);
      ASSERT_TRUE(res1.Ok());
      EXPECT_EQ(res1.codeword, cw);
      EXPECT_EQ(res1.message, msg);
      EXPECT_EQ(res1.num_errors, 1u);

      for (std::size_t j = i + 1; j < 8; ++j) {
        for (uint8_t e2 : {uint8_t{1}, uint8_t{2}}) {
          auto rx2 = rx1;
          rx2.Raw(j) = gf3->Add(rx2.Raw(j), e2);
          auto res2 = tbch.DecodeWithResult(rx2);
          ASSERT_TRUE(res2.Ok());
          EXPECT_EQ(res2.codeword, cw);
          EXPECT_EQ(res2.message, msg);
          EXPECT_EQ(res2.num_errors, 2u);
        }
      }
    }
  }

  // Also test b = 0 (non-narrow-sense) RS-as-BCH over GF(2^4): [15, 11, 5].
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  BCHCode<GF2X<uint8_t>, GF2X<uint8_t>> rs_bch(gf16, gf16, 5, 0);
  EXPECT_EQ(rs_bch.Length(), 15u);
  EXPECT_EQ(rs_bch.Dimension(), 11u);
  ExpectOrthogonal(rs_bch.GeneratorMatrix(), rs_bch.ParityCheckMatrix());

  auto rs_msg = GFArray<GF2X<uint8_t>>::Random(gf16, {11}, rng);
  auto rs_cw = rs_bch.Encode(rs_msg);
  auto rs_rx = rs_cw;
  rs_rx.Raw(3) = gf16->Add(rs_rx.Raw(3), 11);
  rs_rx.Raw(12) = gf16->Add(rs_rx.Raw(12), 7);
  auto rs_res = rs_bch.DecodeWithResult(rs_rx);
  ASSERT_TRUE(rs_res.Ok());
  EXPECT_EQ(rs_res.codeword, rs_cw);
  EXPECT_EQ(rs_res.message, rs_msg);
  EXPECT_EQ(rs_res.num_errors, 2u);
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
