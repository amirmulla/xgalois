#include "xgalois/coding/grs.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace {

using namespace xg;
using namespace xg::coding;

template <typename GaloisField>
bool IsZeroMatrix(const GFMatrix<GaloisField> &mat) {
  const auto field = mat.Field();
  for (const auto &v : mat.RawSpan()) {
    if (!field->IsZero(v)) {
      return false;
    }
  }
  return true;
}

template <typename GaloisField>
GFArray<GaloisField> MakeGFArrayFromInts(
    const std::shared_ptr<GaloisField> &field,
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

template <typename GaloisField>
typename GaloisField::element_type RandomNonZero(
    const std::shared_ptr<GaloisField> &field, std::mt19937_64 &rng) {
  auto v = field->Random(rng);
  while (field->IsZero(v)) {
    v = field->Random(rng);
  }
  return v;
}

template <typename GaloisField>
GFArray<GaloisField> CorruptAtPositions(
    const GFArray<GaloisField> &cw, const std::vector<std::size_t> &positions,
    std::mt19937_64 &rng) {
  auto field = cw.Field();
  auto received = cw;
  for (std::size_t pos : positions) {
    const auto err = RandomNonZero(field, rng);
    received.Raw(pos) = field->Add(received.Raw(pos), err);
  }
  return received;
}

template <typename GaloisField>
void VerifyOrthogonalityAndDual(
    const GeneralizedReedSolomonCode<GaloisField> &code) {
  const auto &G = code.GeneratorMatrix();
  const auto &H = code.ParityCheckMatrix();
  EXPECT_EQ(G.Shape(0), code.Dimension());
  EXPECT_EQ(G.Shape(1), code.Length());
  EXPECT_EQ(H.Shape(0), code.Redundancy());
  EXPECT_EQ(H.Shape(1), code.Length());

  if (code.Redundancy() > 0) {
    const auto prod = linalg::Dot(G, linalg::Transpose(H));
    EXPECT_TRUE(IsZeroMatrix(prod));

    const auto dual = code.DualCode();
    EXPECT_EQ(dual.Length(), code.Length());
    EXPECT_EQ(dual.Dimension(), code.Redundancy());
    EXPECT_EQ(dual.MinimumDistance(), code.Dimension() + 1);
    EXPECT_EQ(dual.GeneratorMatrix(), H);
    const auto dual_prod =
        linalg::Dot(dual.GeneratorMatrix(), linalg::Transpose(G));
    EXPECT_TRUE(IsZeroMatrix(dual_prod));
    const auto dual_self_prod = linalg::Dot(
        dual.GeneratorMatrix(), linalg::Transpose(dual.ParityCheckMatrix()));
    EXPECT_TRUE(IsZeroMatrix(dual_self_prod));
  }
}

template <typename GaloisField>
void VerifyAllDecodersUpToRadius(
    const GeneralizedReedSolomonCode<GaloisField> &code, uint64_t seed) {
  const auto field = code.Field();
  const std::size_t n = code.Length();
  const std::size_t k = code.Dimension();
  const std::size_t t = (code.MinimumDistance() - 1) / 2;
  std::mt19937_64 rng(seed);

  const std::vector<std::string> decoders = {"PGZ", "BerlekampMassey",
                                             "Sugiyama"};
  std::vector<std::size_t> all_indices(n);
  std::iota(all_indices.begin(), all_indices.end(), std::size_t{0});

  // Test both non-zero random messages and the all-zero message.
  const std::vector<GFArray<GaloisField>> messages = {
      GFArray<GaloisField>::Zeros(field, {k}),
      GFArray<GaloisField>::Random(field, {k}, rng),
  };

  for (const auto &msg : messages) {
    const auto cw = code.Encode(msg);
    ASSERT_TRUE(code.Contains(cw));
    EXPECT_EQ(code.Unencode(cw), msg);

    for (std::size_t num_err = 0; num_err <= t; ++num_err) {
      std::shuffle(all_indices.begin(), all_indices.end(), rng);
      std::vector<std::size_t> err_pos(
          all_indices.begin(),
          all_indices.begin() + static_cast<std::ptrdiff_t>(num_err));
      const auto received = CorruptAtPositions(cw, err_pos, rng);

      for (const auto &dec_name : decoders) {
        const auto res = code.DecodeWithResult(received, dec_name);
        ASSERT_TRUE(res.Ok())
            << "Decoder " << dec_name << " failed on " << num_err << " errors";
        EXPECT_EQ(res.codeword, cw) << "Decoder " << dec_name;
        EXPECT_EQ(res.message, msg) << "Decoder " << dec_name;
        EXPECT_EQ(res.num_errors, num_err) << "Decoder " << dec_name;
        EXPECT_EQ(res.num_erasures, 0U);
      }
    }
  }
}

TEST(GRSTest, OrthogonalityAndDualMultiField) {
  // GF(7) with non-trivial multipliers v_j != 1
  {
    auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    auto alphas = MakeGFArrayFromInts(gf7, {1, 2, 3, 4, 5, 6});
    auto vs = MakeGFArrayFromInts(gf7, {2, 5, 3, 6, 4, 1});
    GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>> grs(gf7, alphas, vs,
                                                              2);
    EXPECT_EQ(grs.MinimumDistance(), 5U);
    VerifyOrthogonalityAndDual(grs);
    VerifyAllDecodersUpToRadius(grs, 101);
  }

  // GF(17) with non-trivial multipliers
  {
    auto gf17 = std::make_shared<GaloisFieldPrime<uint16_t>>(17);
    auto alphas =
        MakeGFArrayFromInts(gf17, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
    auto vs =
        MakeGFArrayFromInts(gf17, {3, 14, 7, 2, 9, 11, 5, 16, 4, 8, 13, 6});
    GeneralizedReedSolomonCode<GaloisFieldPrime<uint16_t>> grs(gf17, alphas, vs,
                                                               5);
    EXPECT_EQ(grs.MinimumDistance(), 8U);
    VerifyOrthogonalityAndDual(grs);
    VerifyAllDecodersUpToRadius(grs, 102);
  }

  // GF(2^4) with non-trivial multipliers
  {
    auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
    auto alphas = MakeGFArrayFromInts(
        gf16, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
    auto vs = MakeGFArrayFromInts(
        gf16, {3, 7, 11, 2, 9, 5, 14, 6, 13, 4, 8, 12, 15, 1, 10});
    GeneralizedReedSolomonCode<GaloisFieldBinaryExtension<uint8_t>> grs(
        gf16, alphas, vs, 9);
    EXPECT_EQ(grs.MinimumDistance(), 7U);
    VerifyOrthogonalityAndDual(grs);
    VerifyAllDecodersUpToRadius(grs, 103);
  }

  // GF(2^4) Zech log representation (internal representation differs from int!)
  {
    auto gf16_zech = std::make_shared<GF2XZECH<uint8_t>>(4);
    auto rs = GeneralizedReedSolomonCode<GF2XZECH<uint8_t>>::ReedSolomon(
        gf16_zech, 15, 9);
    VerifyOrthogonalityAndDual(rs);
    VerifyAllDecodersUpToRadius(rs, 104);
  }

  // GF(2^8) log-table field with non-trivial multipliers
  {
    auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
    std::vector<uint64_t> a_vals;
    std::vector<uint64_t> v_vals;
    for (uint64_t i = 1; i <= 20; ++i) {
      a_vals.push_back(i);
      v_vals.push_back((i * 13 + 5) % 255 + 1);
    }
    auto alphas = MakeGFArrayFromInts(gf256, a_vals);
    auto vs = MakeGFArrayFromInts(gf256, v_vals);
    GeneralizedReedSolomonCode<GF2XLOG<uint8_t>> grs(gf256, alphas, vs, 12);
    EXPECT_EQ(grs.MinimumDistance(), 9U);
    VerifyOrthogonalityAndDual(grs);
    VerifyAllDecodersUpToRadius(grs, 105);
  }

  // GF(3^2) general odd-characteristic extension field
  {
    auto gf9 = std::make_shared<GaloisFieldExtension<uint8_t>>(
        std::pair<uint64_t, uint64_t>{3, 2});
    auto alphas = MakeGFArrayFromInts(gf9, {1, 2, 3, 4, 5, 6, 7, 8});
    auto vs = MakeGFArrayFromInts(gf9, {2, 5, 7, 3, 8, 4, 6, 1});
    GeneralizedReedSolomonCode<GaloisFieldExtension<uint8_t>> grs(gf9, alphas,
                                                                  vs, 4);
    EXPECT_EQ(grs.MinimumDistance(), 5U);
    VerifyOrthogonalityAndDual(grs);
    VerifyAllDecodersUpToRadius(grs, 106);
  }
}

TEST(GRSTest, EncodersAllZeroAndSystematicConsistency) {
  auto gf17 = std::make_shared<GaloisFieldPrime<uint16_t>>(17);
  auto alphas = MakeGFArrayFromInts(gf17, {0, 1, 2, 3, 5, 7, 11, 13, 15, 16});
  auto vs = MakeGFArrayFromInts(gf17, {4, 2, 9, 1, 13, 6, 8, 12, 3, 15});
  GeneralizedReedSolomonCode<GaloisFieldPrime<uint16_t>> grs(gf17, alphas, vs,
                                                             5);

  // Zero polynomial encodes cleanly to all-zero codeword (regression test).
  auto zero_poly = PolynomialDense<GaloisFieldPrime<uint16_t>>::Zero(gf17);
  auto zero_cw = grs.EncodePolynomial(zero_poly);
  EXPECT_EQ(zero_cw, GFArray<GaloisFieldPrime<uint16_t>>::Zeros(gf17, {10}));
  EXPECT_TRUE(grs.Contains(zero_cw));

  std::mt19937_64 rng(42);
  auto zero_msg = GFArray<GaloisFieldPrime<uint16_t>>::Zeros(gf17, {5});
  auto rand_msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf17, {5}, rng);

  for (const auto &msg : {zero_msg, rand_msg}) {
    // Non-systematic GRS encoder and GeneratorMatrix encoder must agree.
    auto cw_grs = grs.Encode(msg, "GRS");
    auto cw_gm = grs.Encode(msg, "GeneratorMatrix");
    EXPECT_EQ(cw_grs, cw_gm);
    EXPECT_TRUE(grs.Contains(cw_grs));
    EXPECT_EQ(grs.Unencode(cw_grs, "GRS"), msg);
    EXPECT_EQ(grs.Unencode(cw_gm, "GeneratorMatrix"), msg);

    // Systematic GRS encoder places message verbatim at 0..k-1.
    auto cw_sys = grs.Encode(msg, "GRSSystematic");
    EXPECT_TRUE(grs.Contains(cw_sys));
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_EQ(cw_sys.Raw(i), msg.Raw(i));
    }
    EXPECT_EQ(grs.Unencode(cw_sys, "GRSSystematic"), msg);
  }

  // Non-codeword must be rejected by Unencode.
  auto bad = grs.Encode(rand_msg, "GRS");
  bad.Raw(0) = gf17->Add(bad.Raw(0), gf17->MultiplicativeIdentity());
  EXPECT_FALSE(grs.Contains(bad));
  EXPECT_THROW(grs.Unencode(bad, "GRS"), std::invalid_argument);
  EXPECT_THROW(grs.Unencode(bad, "GRSSystematic"), std::invalid_argument);
}

TEST(GRSTest, UncorrectableErrorsReturnFailure) {
  // Odd redundancy r = 7 => d = 8, t = 3. Every weight-4 error is strictly at
  // distance >= 4 from EVERY codeword, so ALL decoders must return !Ok().
  auto gf17 = std::make_shared<GaloisFieldPrime<uint16_t>>(17);
  auto rs = GeneralizedReedSolomonCode<GaloisFieldPrime<uint16_t>>::ReedSolomon(
      gf17, 14, 7);
  ASSERT_EQ(rs.MinimumDistance(), 8U);

  std::mt19937_64 rng(777);
  const auto msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf17, {7}, rng);
  const auto cw = rs.Encode(msg);
  const auto received = CorruptAtPositions(cw, {0, 3, 7, 12}, rng);

  for (const char *dec : {"PGZ", "BerlekampMassey", "Sugiyama"}) {
    const auto res = rs.DecodeWithResult(received, dec);
    EXPECT_FALSE(res.Ok()) << "Decoder " << dec
                           << " must reject weight-(t+1) error";
    EXPECT_FALSE(res.error_message.empty());
    EXPECT_THROW(rs.DecodeToCode(received, dec), std::runtime_error);
  }
}

TEST(GRSTest, ExtendedRSAndZeroEvaluationPointErrors) {
  // Extended RS over GF(16) of full length q = 16, k = 10 => r = 6, d = 7, t =
  // 3, with alpha_0 = 0.
  auto gf16 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(4);
  auto ext_rs = GeneralizedReedSolomonCode<
      GaloisFieldBinaryExtension<uint8_t>>::ExtendedReedSolomon(gf16, 16, 10);
  ASSERT_TRUE(ext_rs.ZeroEvaluationIndex().has_value());
  EXPECT_EQ(*ext_rs.ZeroEvaluationIndex(), 0U);
  VerifyOrthogonalityAndDual(ext_rs);

  std::mt19937_64 rng(2025);
  const auto msg =
      GFArray<GaloisFieldBinaryExtension<uint8_t>>::Random(gf16, {10}, rng);
  const auto cw = ext_rs.Encode(msg);
  ASSERT_TRUE(ext_rs.Contains(cw));

  const std::vector<std::string> decoders = {"PGZ", "BerlekampMassey",
                                             "Sugiyama"};

  // 1) Single error strictly at alpha_0 = 0.
  {
    const auto rx = CorruptAtPositions(cw, {0}, rng);
    for (const auto &dec : decoders) {
      const auto res = ext_rs.DecodeWithResult(rx, dec);
      ASSERT_TRUE(res.Ok()) << dec << " failed on single error at alpha=0";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 1U);
    }
  }

  // 2) Error at alpha_0 = 0 combined with 1 and 2 non-zero-point errors (total
  // weight 2 and 3 = t).
  for (const auto &pattern :
       std::vector<std::vector<std::size_t>>{{0, 5}, {0, 2, 14}}) {
    const auto rx = CorruptAtPositions(cw, pattern, rng);
    for (const auto &dec : decoders) {
      const auto res = ext_rs.DecodeWithResult(rx, dec);
      ASSERT_TRUE(res.Ok()) << dec << " failed on weight-" << pattern.size()
                            << " pattern including alpha=0";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, pattern.size());
    }
  }

  // 3) GRS over GF(7) with alpha_j = 0 at interior position j_0 = 3 and
  // non-trivial column multipliers.
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto alphas = MakeGFArrayFromInts(gf7, {2, 5, 1, 0, 6, 3, 4});
  auto vs = MakeGFArrayFromInts(gf7, {3, 1, 6, 4, 2, 5, 3});
  GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>> grs7(gf7, alphas, vs,
                                                             3);
  ASSERT_TRUE(grs7.ZeroEvaluationIndex().has_value());
  EXPECT_EQ(*grs7.ZeroEvaluationIndex(), 3U);
  VerifyOrthogonalityAndDual(grs7);

  const auto msg7 = GFArray<GaloisFieldPrime<uint8_t>>::Random(gf7, {3}, rng);
  const auto cw7 = grs7.Encode(msg7);
  for (const auto &pattern :
       std::vector<std::vector<std::size_t>>{{3}, {3, 0}, {3, 6}, {1, 5}}) {
    const auto rx = CorruptAtPositions(cw7, pattern, rng);
    for (const auto &dec : decoders) {
      const auto res = grs7.DecodeWithResult(rx, dec);
      ASSERT_TRUE(res.Ok()) << dec << " failed on GF(7) interior alpha=0";
      EXPECT_EQ(res.codeword, cw7);
      EXPECT_EQ(res.message, msg7);
      EXPECT_EQ(res.num_errors, pattern.size());
    }
  }
}

TEST(GRSTest, ErrorsAndErasuresAllCombinationsIncludingZeroPoint) {
  // Extended RS over GF(17) of length 17, k = 11 => r = 6, alpha_0 = 0.
  auto gf17 = std::make_shared<GaloisFieldPrime<uint16_t>>(17);
  auto ext_rs = GeneralizedReedSolomonCode<
      GaloisFieldPrime<uint16_t>>::ExtendedReedSolomon(gf17, 17, 11);
  const std::size_t n = ext_rs.Length();
  const std::size_t r = ext_rs.Redundancy();  // 6

  std::mt19937_64 rng(9001);
  const auto msg = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf17, {11}, rng);
  const auto cw = ext_rs.Encode(msg);

  const std::vector<std::string> decoders = {"PGZ", "BerlekampMassey",
                                             "Sugiyama"};
  std::vector<std::size_t> indices(n);
  std::iota(indices.begin(), indices.end(), std::size_t{0});

  // Exhaustively test every valid (e, s) with 2e + s <= r = 6.
  for (std::size_t s = 0; s <= r; ++s) {
    for (std::size_t e = 0; 2 * e + s <= r; ++e) {
      std::shuffle(indices.begin(), indices.end(), rng);
      std::vector<std::size_t> erasures(
          indices.begin(), indices.begin() + static_cast<std::ptrdiff_t>(s));
      std::vector<std::size_t> errors(
          indices.begin() + static_cast<std::ptrdiff_t>(s),
          indices.begin() + static_cast<std::ptrdiff_t>(s + e));

      // Corrupt both erased positions (arbitrary garbage) and unknown errors.
      auto rx = CorruptAtPositions(cw, errors, rng);
      for (std::size_t p : erasures) {
        rx.Raw(p) = gf17->Random(rng);
      }

      for (const auto &dec : decoders) {
        const auto res = ext_rs.DecodeWithErasures(rx, erasures, dec);
        ASSERT_TRUE(res.Ok())
            << dec << " failed for (e=" << e << ", s=" << s << ")";
        EXPECT_EQ(res.codeword, cw);
        EXPECT_EQ(res.message, msg);
        EXPECT_EQ(res.num_errors, e);
        EXPECT_EQ(res.num_erasures, s);
      }
    }
  }

  // Explicit test: erasure AT alpha_0 = 0 combined with another erasure and 2
  // unknown errors (2*2 + 2 = 6 = r).
  {
    const std::vector<std::size_t> erasures = {0, 4};
    const std::vector<std::size_t> errors = {2, 13};
    auto rx = CorruptAtPositions(cw, errors, rng);
    rx.Raw(0) = gf17->Random(rng);
    rx.Raw(4) = gf17->Random(rng);
    for (const auto &dec : decoders) {
      const auto res = ext_rs.DecodeWithErasures(rx, erasures, dec);
      ASSERT_TRUE(res.Ok()) << dec << " failed with erasure at alpha=0";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 2U);
      EXPECT_EQ(res.num_erasures, 2U);
    }
  }

  // Explicit test: unknown error AT alpha_0 = 0 while positions {1, 5} are
  // erased plus 1 more error at position 10 (2*2 + 2 = 6 = r).
  {
    const std::vector<std::size_t> erasures = {1, 5};
    const std::vector<std::size_t> errors = {0, 10};
    auto rx = CorruptAtPositions(cw, errors, rng);
    rx.Raw(1) = gf17->Random(rng);
    rx.Raw(5) = gf17->Random(rng);
    for (const auto &dec : decoders) {
      const auto res = ext_rs.DecodeWithErasures(rx, erasures, dec);
      ASSERT_TRUE(res.Ok())
          << dec << " failed with error at alpha=0 + erasures";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 2U);
      EXPECT_EQ(res.num_erasures, 2U);
    }
  }

  // Over-budget erasures (s = r + 1), duplicate erasures, out-of-range
  // erasures, and 2e + s > r must all return Failure.
  for (const auto &dec : decoders) {
    EXPECT_FALSE(
        ext_rs.DecodeWithErasures(cw, {0, 1, 2, 3, 4, 5, 6}, dec).Ok());
    EXPECT_FALSE(ext_rs.DecodeWithErasures(cw, {2, 2}, dec).Ok());
    EXPECT_FALSE(ext_rs.DecodeWithErasures(cw, {n}, dec).Ok());

    // Odd remaining redundancy check: s = 1 erasure, e = 3 errors => 2e+s = 7 >
    // 6. Use r = 7 code so 2*3 + 2 = 8 > 7 is strictly uncorrectable.
    auto rs_odd = GeneralizedReedSolomonCode<
        GaloisFieldPrime<uint16_t>>::ExtendedReedSolomon(gf17, 17, 10);
    auto m_odd = GFArray<GaloisFieldPrime<uint16_t>>::Random(gf17, {10}, rng);
    auto c_odd = rs_odd.Encode(m_odd);
    // s = 2 erasures => max errors = (7 - 2) / 2 = 2; inject 3 errors (2*3+2=8
    // > 7, and remaining parity distance after puncturing 2 erasures is 8-2=6,
    // so every weight-3 error is at distance >= 3 from every other codeword!).
    auto rx_bad = CorruptAtPositions(c_odd, {3, 7, 12}, rng);
    EXPECT_FALSE(rs_odd.DecodeWithErasures(rx_bad, {1, 5}, dec).Ok());
  }
}

TEST(GRSTest, CCSDSReedSolomon255_223) {
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
  auto rs = GeneralizedReedSolomonCode<GF2XLOG<uint8_t>>::ReedSolomon(gf256,
                                                                      255, 223);

  EXPECT_EQ(rs.Length(), 255U);
  EXPECT_EQ(rs.Dimension(), 223U);
  EXPECT_EQ(rs.Redundancy(), 32U);
  EXPECT_EQ(rs.MinimumDistance(), 33U);
  EXPECT_EQ(rs.GetDecoder("BerlekampMassey")->DecodingRadius(), 16U);
  EXPECT_EQ(rs.GetDecoder("Sugiyama")->DecodingRadius(), 16U);

  // Verify exact G * H^T = 0 on the full 223 x 255 and 32 x 255 matrices.
  const auto ght = linalg::Dot(rs.GeneratorMatrix(),
                               linalg::Transpose(rs.ParityCheckMatrix()));
  EXPECT_TRUE(IsZeroMatrix(ght));

  std::mt19937_64 rng(0xCC5D5255ULL);
  std::vector<std::size_t> all_pos(255);
  std::iota(all_pos.begin(), all_pos.end(), std::size_t{0});

  const std::vector<GFArray<GF2XLOG<uint8_t>>> test_messages = {
      GFArray<GF2XLOG<uint8_t>>::Zeros(gf256, {223}),
      GFArray<GF2XLOG<uint8_t>>::Random(gf256, {223}, rng),
  };

  for (const auto &msg : test_messages) {
    const auto cw = rs.Encode(msg);
    ASSERT_TRUE(rs.Contains(cw));
    EXPECT_EQ(rs.Unencode(cw), msg);

    // Systematic encoding round-trip.
    const auto sys_cw = rs.Encode(msg, "GRSSystematic");
    ASSERT_TRUE(rs.Contains(sys_cw));
    EXPECT_EQ(rs.Unencode(sys_cw, "GRSSystematic"), msg);

    // 1) Correct 16 random errors with BOTH Berlekamp-Massey and Sugiyama.
    std::shuffle(all_pos.begin(), all_pos.end(), rng);
    std::vector<std::size_t> err16(all_pos.begin(), all_pos.begin() + 16);
    const auto rx16 = CorruptAtPositions(cw, err16, rng);

    for (const char *dec : {"BerlekampMassey", "Sugiyama"}) {
      const auto res = rs.DecodeWithResult(rx16, dec);
      ASSERT_TRUE(res.Ok()) << dec << " failed on 16 errors in RS(255,223)";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 16U);
      EXPECT_EQ(res.num_erasures, 0U);
    }

    // 2) Correct 32 random erasures with BOTH Berlekamp-Massey and Sugiyama.
    std::shuffle(all_pos.begin(), all_pos.end(), rng);
    std::vector<std::size_t> era32(all_pos.begin(), all_pos.begin() + 32);
    auto rx32 = cw;
    for (std::size_t p : era32) {
      rx32.Raw(p) = gf256->Random(rng);
    }
    for (const char *dec : {"BerlekampMassey", "Sugiyama"}) {
      const auto res = rs.DecodeWithErasures(rx32, era32, dec);
      ASSERT_TRUE(res.Ok()) << dec << " failed on 32 erasures in RS(255,223)";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 0U);
      EXPECT_EQ(res.num_erasures, 32U);
    }

    // 3) Correct 10 random errors + 12 random erasures (2*10 + 12 = 32).
    std::shuffle(all_pos.begin(), all_pos.end(), rng);
    std::vector<std::size_t> era12(all_pos.begin(), all_pos.begin() + 12);
    std::vector<std::size_t> err10(all_pos.begin() + 12, all_pos.begin() + 22);
    auto rx_mix = CorruptAtPositions(cw, err10, rng);
    for (std::size_t p : era12) {
      rx_mix.Raw(p) = gf256->Random(rng);
    }
    for (const char *dec : {"BerlekampMassey", "Sugiyama"}) {
      const auto res = rs.DecodeWithErasures(rx_mix, era12, dec);
      ASSERT_TRUE(res.Ok())
          << dec << " failed on 10 errors + 12 erasures in RS(255,223)";
      EXPECT_EQ(res.codeword, cw);
      EXPECT_EQ(res.message, msg);
      EXPECT_EQ(res.num_errors, 10U);
      EXPECT_EQ(res.num_erasures, 12U);
    }
  }

  // 4) Uncorrectable 17 random errors must fail cleanly.
  const auto msg = GFArray<GF2XLOG<uint8_t>>::Random(gf256, {223}, rng);
  const auto cw = rs.Encode(msg);
  std::shuffle(all_pos.begin(), all_pos.end(), rng);
  std::vector<std::size_t> err17(all_pos.begin(), all_pos.begin() + 17);
  const auto rx17 = CorruptAtPositions(cw, err17, rng);
  for (const char *dec : {"BerlekampMassey", "Sugiyama"}) {
    const auto res = rs.DecodeWithResult(rx17, dec);
    EXPECT_FALSE(res.Ok()) << dec << " must fail on 17 errors in RS(255,223)";
  }
}

TEST(GRSTest, ValidationAndErrorHandling) {
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto dup_alphas = MakeGFArrayFromInts(gf7, {1, 2, 2, 4});
  auto ok_alphas = MakeGFArrayFromInts(gf7, {1, 2, 3, 4});
  auto ok_vs = MakeGFArrayFromInts(gf7, {1, 2, 3, 4});
  auto zero_vs = MakeGFArrayFromInts(gf7, {1, 0, 3, 4});

  EXPECT_THROW((GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>>(
                   gf7, dup_alphas, ok_vs, 2)),
               std::invalid_argument);
  EXPECT_THROW((GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>>(
                   gf7, ok_alphas, zero_vs, 2)),
               std::invalid_argument);
  EXPECT_THROW((GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>>(
                   gf7, ok_alphas, ok_vs, 0)),
               std::invalid_argument);
  EXPECT_THROW((GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>>(
                   gf7, ok_alphas, ok_vs, 5)),
               std::invalid_argument);
  EXPECT_THROW(
      GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>>::ReedSolomon(gf7, 7,
                                                                         3),
      std::invalid_argument);
  EXPECT_THROW(GeneralizedReedSolomonCode<
                   GaloisFieldPrime<uint8_t>>::ExtendedReedSolomon(gf7, 8, 3),
               std::invalid_argument);

  auto rs = GeneralizedReedSolomonCode<GaloisFieldPrime<uint8_t>>::ReedSolomon(
      gf7, 6, 3);
  auto high_deg = PolynomialDense<GaloisFieldPrime<uint8_t>>::Monomial(
      gf7, 3, gf7->MultiplicativeIdentity());
  EXPECT_THROW(rs.EncodePolynomial(high_deg), std::invalid_argument);
}

}  // namespace
