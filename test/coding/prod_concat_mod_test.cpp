#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "xgalois/coding/concatenated.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/linear_code.hpp"
#include "xgalois/coding/modifiers.hpp"
#include "xgalois/coding/product.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/utils/bigint.hpp"

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
  auto field = G.Field();
  for (const auto &v : prod.RawSpan()) {
    EXPECT_TRUE(field->IsZero(v));
  }
}

BigInt BinomialBigInt(std::size_t n, std::size_t k) {
  if (k > n) {
    return BigInt(0);
  }
  if (k == 0 || k == n) {
    return BigInt(1);
  }
  k = std::min(k, n - k);
  BigInt res(1);
  for (std::size_t i = 1; i <= k; ++i) {
    res = (res * BigInt::FromUint64(n - k + i)) / BigInt::FromUint64(i);
  }
  return res;
}

// Closed-form MDS weight enumerator for an [n, k, d = n - k + 1] MDS code
// over GF(q):
//   A_0 = 1, A_w = 0 for 1 <= w < d,
//   A_w = C(n, w) * (q - 1) * sum_{j=0}^{w - d} (-1)^j C(w - 1, j) q^{w - d -
//   j}.
BigInt ExactMdsWeightCount(std::size_t n, std::size_t k, uint64_t q,
                           std::size_t w) {
  const std::size_t d = n - k + 1;
  if (w == 0) {
    return BigInt(1);
  }
  if (w < d || w > n) {
    return BigInt(0);
  }
  BigInt sum(0);
  const BigInt q_big = BigInt::FromUint64(q);
  for (std::size_t j = 0; j <= w - d; ++j) {
    BigInt term = BinomialBigInt(w - 1, j) * q_big.Pow(w - d - j);
    if ((j & std::size_t{1}) != 0) {
      sum -= term;
    } else {
      sum += term;
    }
  }
  return BinomialBigInt(n, w) * BigInt::FromUint64(q - 1) * sum;
}

// ============================================================================
// 1. ProductCode Tests
// ============================================================================

TEST(ProductCodeTest, Rs75xRs75ParametersEncodingAndIterativeBurstRecovery) {
  using GF8 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf8 = std::make_shared<GF8>(3);

  auto rs_row = std::make_shared<GeneralizedReedSolomonCode<GF8>>(
      GeneralizedReedSolomonCode<GF8>::ReedSolomon(gf8, 7, 5));
  auto rs_col = std::make_shared<GeneralizedReedSolomonCode<GF8>>(
      GeneralizedReedSolomonCode<GF8>::ReedSolomon(gf8, 7, 5));

  ProductCode<GF8> prod(rs_row, rs_col);
  EXPECT_EQ(prod.Length(), 49u);
  EXPECT_EQ(prod.Dimension(), 25u);
  EXPECT_EQ(prod.MinimumDistance(), 9u);
  EXPECT_EQ(prod.RowLength(), 7u);
  EXPECT_EQ(prod.ColLength(), 7u);

  // Verify lazy Kronecker G and parity-check matrix H.
  const auto &G = prod.GeneratorMatrix();
  const auto &H = prod.ParityCheckMatrix();
  EXPECT_EQ(G.Shape(0), 25u);
  EXPECT_EQ(G.Shape(1), 49u);
  EXPECT_EQ(H.Shape(0), 24u);
  EXPECT_EQ(H.Shape(1), 49u);
  ExpectOrthogonal(G, H);

  std::mt19937_64 rng(42);
  auto msg = GFArray<GF8>::Random(gf8, {25}, rng);
  auto cw_fast = prod.Encode(msg);
  auto cw_mat_enc = prod.Encode(msg, "GeneratorMatrix");
  EXPECT_EQ(cw_fast, cw_mat_enc);
  EXPECT_TRUE(prod.Contains(cw_fast));
  EXPECT_EQ(prod.Unencode(cw_fast), msg);

  // Inject a 2x2 burst (rows 1 & 4, cols 2 & 5): choose error values e2, e5
  // satisfying w_2 * e2 + w_5 * e5 = 0 so S_0 = 0, S_1 != 0 in RS(7, 5, 3).
  // That syndrome has NO weight-1 coset leader, so naive row decoding fails
  // on rows 1 & 4; ProductCode flags rows 1 & 4 as erasures for the column
  // pass, and col_code RS(7, 5, 3) corrects 2 erasures (2*0 + 2 < 3)!
  auto cw_mat = prod.UnflattenMatrix(cw_fast, 7, 7);
  auto corrupted_mat = cw_mat;
  const auto &w_dual = rs_row->DualMultipliers();
  const uint8_t e2 = gf8->SetElementValue(3);
  const uint8_t e5 =
      gf8->Div(gf8->Neg(gf8->Mul(w_dual.Raw(2), e2)), w_dual.Raw(5));
  for (std::size_t r : {std::size_t{1}, std::size_t{4}}) {
    corrupted_mat.Raw(r, 2) = gf8->Add(corrupted_mat.Raw(r, 2), e2);
    corrupted_mat.Raw(r, 5) = gf8->Add(corrupted_mat.Raw(r, 5), e5);
  }
  EXPECT_FALSE(prod.ContainsMatrix(corrupted_mat));

  // Verify naive row-only decoding genuinely fails on row 1.
  {
    using RawArray = typename GFArray<GF8>::raw_array_type;
    using ShapeType = typename GFArray<GF8>::shape_type;
    RawArray row1_raw = RawArray::from_shape(ShapeType{7});
    for (std::size_t c = 0; c < 7; ++c) {
      row1_raw(c) = corrupted_mat.Raw(1, c);
    }
    auto row1 = GFArray<GF8>::FromRaw(gf8, std::move(row1_raw));
    EXPECT_FALSE(rs_row->DecodeWithResult(row1).Ok());
  }

  // Iterative hard matrix decoder recovers the exact 2D codeword & message.
  auto mat_res = prod.DecodeMatrixHard(corrupted_mat);
  ASSERT_TRUE(mat_res.Ok()) << mat_res.error_message;
  EXPECT_EQ(mat_res.codeword_matrix, cw_mat);
  EXPECT_EQ(mat_res.num_errors, 4u);
  EXPECT_EQ(prod.FlattenMatrix(mat_res.message_matrix), msg);

  // Registered "ProductHard" 1D decoder also succeeds.
  auto flat_recv = prod.FlattenMatrix(corrupted_mat);
  auto flat_res = prod.DecodeWithResult(flat_recv, "ProductHard");
  ASSERT_TRUE(flat_res.Ok()) << flat_res.error_message;
  EXPECT_EQ(flat_res.codeword, cw_fast);
  EXPECT_EQ(flat_res.message, msg);
  EXPECT_EQ(flat_res.num_errors, 4u);
}

TEST(ProductCodeTest, NonSquareProductCodeIndexing) {
  using GF8 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf8 = std::make_shared<GF8>(3);

  auto rs_row = std::make_shared<GeneralizedReedSolomonCode<GF8>>(
      GeneralizedReedSolomonCode<GF8>::ReedSolomon(gf8, 7, 5));
  auto rs_col = std::make_shared<GeneralizedReedSolomonCode<GF8>>(
      GeneralizedReedSolomonCode<GF8>::ReedSolomon(gf8, 6, 4));

  ProductCode<GF8> prod(rs_row, rs_col);
  EXPECT_EQ(prod.Length(), 42u);
  EXPECT_EQ(prod.Dimension(), 20u);
  EXPECT_EQ(prod.MinimumDistance(), 9u);

  std::mt19937_64 rng(77);
  auto msg = GFArray<GF8>::Random(gf8, {20}, rng);
  auto cw_fast = prod.Encode(msg);
  auto cw_gen = prod.Encode(msg, "GeneratorMatrix");
  EXPECT_EQ(cw_fast, cw_gen);
  EXPECT_EQ(prod.Unencode(cw_fast), msg);
}

TEST(ProductCodeTest, ExtendedHamming84ChasePyndiahSoftTpcAwgn) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto ext_ham = std::make_shared<LinearCode<GaloisFieldBinary>>(
      ExtendedHammingCode(gf2, 3));
  ProductCode<GaloisFieldBinary> tpc(ext_ham, ext_ham);

  EXPECT_EQ(tpc.Length(), 64u);
  EXPECT_EQ(tpc.Dimension(), 16u);
  EXPECT_EQ(tpc.MinimumDistance(), 16u);

  std::mt19937_64 rng(2026);
  constexpr double kSigma = 0.85;
  std::normal_distribution<double> noise(0.0, kSigma);

  std::size_t soft_successes = 0;
  std::size_t hard_successes = 0;
  std::size_t hard_failures_recovered_by_soft = 0;
  constexpr std::size_t kTrials = 30;

  for (std::size_t trial = 0; trial < kTrials; ++trial) {
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {16}, rng);
    auto cw = tpc.Encode(msg);
    const auto cw_span = cw.RawSpan();

    std::vector<double> llrs(64);
    auto hard_recv = GFArray<GaloisFieldBinary>::Zeros(gf2, {64});
    const double sigma2 = kSigma * kSigma;
    for (std::size_t i = 0; i < 64; ++i) {
      const double tx = (cw_span[i] == 0) ? 1.0 : -1.0;
      const double y = tx + noise(rng);
      llrs[i] = 2.0 * y / sigma2;
      hard_recv.Raw(i) = (llrs[i] < 0.0) ? uint8_t{1} : uint8_t{0};
    }

    auto hard_res = tpc.DecodeWithResult(hard_recv, "ProductHard");
    const bool hard_ok = hard_res.Ok() && hard_res.codeword == cw;
    if (hard_ok) {
      ++hard_successes;
    }

    auto soft_res = tpc.DecodeSoftChasePyndiah(llrs, 4, 4);
    const bool soft_ok =
        soft_res.Ok() && tpc.FlattenMatrix(soft_res.codeword_matrix) == cw;
    if (soft_ok) {
      ++soft_successes;
      if (!hard_ok) {
        ++hard_failures_recovered_by_soft;
      }
    }
  }

  // Chase-Pyndiah soft TPC should substantially outperform hard decoding.
  EXPECT_GE(soft_successes, 24u);
  EXPECT_GT(soft_successes, hard_successes);
  EXPECT_GE(hard_failures_recovered_by_soft, 3u);

  // Deterministic multi-error pattern (6 bit errors in a 2x3 block where
  // hard decoding fails because Extended Hamming [8,4,4] only corrects 1 bit).
  auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {16}, rng);
  auto cw = tpc.Encode(msg);
  const auto cw_span = cw.RawSpan();
  std::vector<double> llrs(64);
  for (std::size_t i = 0; i < 64; ++i) {
    llrs[i] = (cw_span[i] == 0) ? 3.0 : -3.0;
  }
  // Flip 6 positions in rows {0, 1} x cols {0, 1, 2} with low reliability.
  for (std::size_t r : {std::size_t{0}, std::size_t{1}}) {
    for (std::size_t c : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
      const std::size_t idx = r * 8 + c;
      llrs[idx] = (cw_span[idx] == 0) ? -0.4 : 0.4;
    }
  }
  auto soft_res = tpc.DecodeSoftChasePyndiah(llrs, 3, 4);
  ASSERT_TRUE(soft_res.Ok());
  EXPECT_EQ(tpc.FlattenMatrix(soft_res.codeword_matrix), cw);
  EXPECT_EQ(tpc.FlattenMatrix(soft_res.message_matrix), msg);
}

// ============================================================================
// 2. ConcatenatedCode & Forney GMD Tests
// ============================================================================

TEST(ConcatenatedCodeTest, Rs73OverGf8WithSimplex734StandardAndForneyGmd) {
  using GF2 = GaloisFieldBinary;
  using GF8 = GaloisFieldBinaryExtension<uint8_t>;

  auto gf2 = std::make_shared<GF2>();
  auto gf8 = std::make_shared<GF8>(3);

  // Outer code: RS(7, 3, 5) over GF(8).
  auto outer = std::make_shared<GeneralizedReedSolomonCode<GF8>>(
      GeneralizedReedSolomonCode<GF8>::ReedSolomon(gf8, 7, 3));
  // Inner code: Simplex [7, 3, 4] over GF(2) (dual of Hamming [7, 4, 3]).
  auto hamming74 = HammingCode(gf2, 3);
  auto inner = std::make_shared<LinearCode<GF2>>(hamming74.DualCode());
  ASSERT_EQ(inner->Length(), 7u);
  ASSERT_EQ(inner->Dimension(), 3u);
  ASSERT_EQ(inner->MinimumDistance(), 4u);

  ConcatenatedCode<GF2, GF8> concat(outer, inner);
  EXPECT_EQ(concat.Length(), 49u);
  EXPECT_EQ(concat.Dimension(), 9u);
  EXPECT_EQ(concat.DesignedDistance(), 20u);
  EXPECT_EQ(concat.MinimumDistance(), 20u);

  // Verify GeneratorMatrix / ParityCheckMatrix consistency.
  const auto &G = concat.GeneratorMatrix();
  const auto &H = concat.ParityCheckMatrix();
  EXPECT_EQ(G.Shape(0), 9u);
  EXPECT_EQ(G.Shape(1), 49u);
  EXPECT_EQ(H.Shape(0), 40u);
  EXPECT_EQ(H.Shape(1), 49u);
  ExpectOrthogonal(G, H);

  std::mt19937_64 rng(1234);
  auto msg = GFArray<GF2>::Random(gf2, {9}, rng);
  auto cw = concat.Encode(msg);
  EXPECT_EQ(cw, concat.Encode(msg, "GeneratorMatrix"));
  EXPECT_TRUE(concat.Contains(cw));
  EXPECT_EQ(concat.Unencode(cw), msg);

  // Find a non-zero weight-4 codeword of Simplex [7, 3, 4] so that flipping 3
  // of its 4 support positions causes inner decoding to genuinely MISDECODE
  // (at distance 1 from that wrong codeword) rather than fail!
  std::vector<std::size_t> wt4_support;
  for (const auto &in_cw : inner->GetCodewords()) {
    if (inner->Weight(in_cw) == 4u) {
      for (std::size_t i = 0; i < 7; ++i) {
        if (in_cw.Raw(i) == 1) {
          wt4_support.push_back(i);
        }
      }
      break;
    }
  }
  ASSERT_EQ(wt4_support.size(), 4u);

  // Helper to inject a given number of bit errors in each of the 7 blocks:
  //   - 3 errors: flip first 3 bits of wt4_support -> inner MISDECODES (e_j=1).
  //   - 2 errors: flip first 2 bits of wt4_support -> inner FAILS (erasure).
  //   - 1 error : flip first bit -> inner CORRECTS (e_j=1).
  auto corrupt_blocks = [&](const std::vector<std::size_t> &per_block_errors) {
    auto recv = cw;
    for (std::size_t b = 0; b < 7; ++b) {
      for (std::size_t e = 0; e < per_block_errors[b]; ++e) {
        recv.Raw(b * 7 + wt4_support[e]) ^= uint8_t{1};
      }
    }
    return recv;
  };

  // Standard Inner -> Outer succeeds when 2*e_out + s_out < d_out = 5:
  // Blocks 0, 1 have 2 bit errors (inner fails -> 2 erasures);
  // Block 2 has 3 bit errors (inner misdecodes -> 1 outer error);
  // Remaining blocks have 0 or 1 bit errors (inner corrects).
  // Total: 2*1 + 2 = 4 < 5 -> Standard succeeds!
  {
    auto recv = corrupt_blocks({2, 2, 3, 1, 1, 0, 0});
    auto std_res = concat.DecodeWithResult(recv, "ConcatenatedStandard");
    ASSERT_TRUE(std_res.Ok()) << std_res.error_message;
    EXPECT_EQ(std_res.codeword, cw);
    EXPECT_EQ(std_res.message, msg);
    EXPECT_EQ(std_res.num_errors, 9u);
  }

  // Forney GMD corrects up to floor((d_out * d_in - 1) / 2) = 9 bit errors
  // even when naive hard inner -> outer FAILS!
  // Pattern (3, 3, 3, 0, 0, 0, 0): 9 bit errors total.
  //   - Blocks 0, 1, 2 each have 3 bit errors -> inner Simplex [7,3,4]
  //     misdecodes each to a wrong codeword at distance 1 (e_j = 1, theta =
  //     0.5).
  //   - That creates 3 outer errors (2*3 = 6 >= 5) -> Standard FAILS!
  //   - Forney GMD erases the 4 least-reliable blocks (including 0, 1, 2) ->
  //     0 outer errors + 4 erasures < 5 -> GMD RECOVERS EXACT CODEWORD!
  {
    auto recv = corrupt_blocks({3, 3, 3, 0, 0, 0, 0});
    auto std_res = concat.DecodeWithResult(recv, "ConcatenatedStandard");
    EXPECT_FALSE(std_res.Ok() && std_res.codeword == cw);

    auto gmd_res = concat.DecodeGmdHard(recv);
    ASSERT_TRUE(gmd_res.Ok()) << gmd_res.error_message;
    EXPECT_EQ(gmd_res.codeword, cw);
    EXPECT_EQ(gmd_res.message, msg);
    EXPECT_EQ(gmd_res.num_errors, 9u);

    auto gmd_reg = concat.DecodeWithResult(recv, "ConcatenatedGMD");
    ASSERT_TRUE(gmd_reg.Ok());
    EXPECT_EQ(gmd_reg.codeword, cw);
  }

  // Pattern (3, 3, 2, 1, 0, 0, 0): 9 bit errors total.
  //   - 2 outer errors + 1 outer erasure => 2*2 + 1 = 5 >= 5 -> Standard FAILS!
  //   - GMD trial s = 4 erases blocks {2 (failed), 0, 1, 3 (theta=0.5)} ->
  //     0 outer errors + 4 erasures < 5 -> GMD SUCCEEDS!
  {
    auto recv = corrupt_blocks({3, 3, 2, 1, 0, 0, 0});
    auto std_res = concat.DecodeWithResult(recv, "ConcatenatedStandard");
    EXPECT_FALSE(std_res.Ok() && std_res.codeword == cw);

    auto gmd_res = concat.DecodeGmdHard(recv);
    ASSERT_TRUE(gmd_res.Ok()) << gmd_res.error_message;
    EXPECT_EQ(gmd_res.codeword, cw);
    EXPECT_EQ(gmd_res.message, msg);
    EXPECT_EQ(gmd_res.num_errors, 9u);
  }

  // Explicit reliability-weighted DecodeGmd call.
  {
    auto recv = corrupt_blocks({3, 3, 3, 0, 0, 0, 0});
    std::vector<double> rels = {0.1, 0.2, 0.3, 0.9, 0.95, 1.0, 1.0};
    auto gmd_res = concat.DecodeGmd(recv, rels);
    ASSERT_TRUE(gmd_res.Ok());
    EXPECT_EQ(gmd_res.codeword, cw);
    EXPECT_EQ(gmd_res.message, msg);
  }
}

TEST(ConcatenatedCodeTest, MultiStreamConcatenationWhenMDividesKIn) {
  using GF2 = GaloisFieldBinary;
  using GF4 = GaloisFieldBinaryExtension<uint8_t>;

  auto gf2 = std::make_shared<GF2>();
  auto gf4 = std::make_shared<GF4>(2);  // relative degree m = 2 over GF(2)

  // Outer: RS(3, 1, 3) over GF(4); Inner: Hamming [7, 4, 3] over GF(2) (k_in =
  // 4 = 2 * 2).
  auto outer = std::make_shared<GeneralizedReedSolomonCode<GF4>>(
      GeneralizedReedSolomonCode<GF4>::ReedSolomon(gf4, 3, 1));
  auto inner = std::make_shared<LinearCode<GF2>>(HammingCode(gf2, 3));

  ConcatenatedCode<GF2, GF4> concat(outer, inner);
  EXPECT_EQ(concat.Length(), 21u);
  EXPECT_EQ(concat.Dimension(), 4u);
  EXPECT_EQ(concat.InterleavedStreams(), 2u);
  EXPECT_EQ(concat.DesignedDistance(), 9u);

  std::mt19937_64 rng(99);
  auto msg = GFArray<GF2>::Random(gf2, {4}, rng);
  auto cw = concat.Encode(msg);
  EXPECT_EQ(cw, concat.Encode(msg, "GeneratorMatrix"));
  EXPECT_TRUE(concat.Contains(cw));
  EXPECT_EQ(concat.Unencode(cw), msg);

  // Inject 3 bit errors across blocks (2 in block 0, 1 in block 1).
  auto recv = cw;
  recv.Raw(0) ^= 1;
  recv.Raw(1) ^= 1;
  recv.Raw(7) ^= 1;
  auto res = concat.DecodeGmdHard(recv);
  ASSERT_TRUE(res.Ok());
  EXPECT_EQ(res.codeword, cw);
  EXPECT_EQ(res.message, msg);
}

// ============================================================================
// 3. Code Modifiers & PlotkinCode Tests
// ============================================================================

TEST(ModifiersTest, PunctureShortenExtendExpurgateAugment) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto ham = HammingCode(gf2, 3);  // [7, 4, 3]

  // ExtendedCode: Hamming [7, 4, 3] -> [8, 4, 4] with even-weight codewords.
  auto ext = ExtendedCode(ham);
  EXPECT_EQ(ext.Length(), 8u);
  EXPECT_EQ(ext.Dimension(), 4u);
  EXPECT_EQ(ext.MinimumDistance(), 4u);
  for (const auto &cw : ext.GetCodewords()) {
    EXPECT_EQ(ext.Weight(cw) % 2u, 0u);
  }

  // PuncturedCode: puncturing the parity bit of Extended Hamming [8, 4, 4]
  // restores a [7, 4, 3] code.
  auto punc = PuncturedCode(ext, {7});
  EXPECT_EQ(punc.Length(), 7u);
  EXPECT_EQ(punc.Dimension(), 4u);
  EXPECT_EQ(punc.MinimumDistance(), 3u);

  // Multi-coordinate puncturing.
  auto punc2 = PuncturedCode(ext, {0, 7});
  EXPECT_EQ(punc2.Length(), 6u);
  EXPECT_EQ(punc2.Dimension(), 4u);
  EXPECT_EQ(punc2.MinimumDistance(), 2u);

  // ShortenedCode: shortening Hamming [7, 4, 3] at {0, 6} -> [5, 2, >=3].
  auto shrt = ShortenedCode(ham, {0, 6});
  EXPECT_EQ(shrt.Length(), 5u);
  EXPECT_EQ(shrt.Dimension(), 2u);
  EXPECT_GE(shrt.MinimumDistance(), 3u);
  // Every shortened codeword padded with 0 at {0, 6} must lie in Hamming
  // [7,4,3].
  for (const auto &cw : shrt.GetCodewords()) {
    auto padded = GFArray<GaloisFieldBinary>::Zeros(gf2, {7});
    for (std::size_t i = 0; i < 5; ++i) {
      padded.Raw(i + 1) = cw.Raw(i);
    }
    EXPECT_TRUE(ham.Contains(padded));
  }

  // ExpurgatedCode: default all-ones parity check removes all odd-weight words
  // from Hamming [7, 4, 3] -> [7, 3, 4].
  auto expurg = ExpurgatedCode(ham);
  EXPECT_EQ(expurg.Length(), 7u);
  EXPECT_EQ(expurg.Dimension(), 3u);
  EXPECT_EQ(expurg.MinimumDistance(), 4u);
  for (const auto &cw : expurg.GetCodewords()) {
    EXPECT_EQ(expurg.Weight(cw) % 2u, 0u);
    EXPECT_TRUE(ham.Contains(cw));
  }

  // AugmentedCode: adjoining all-ones row to Simplex [7, 3, 4] -> [7, 4, 3].
  auto simplex = ham.DualCode();
  auto aug = AugmentedCode(simplex);
  EXPECT_EQ(aug.Length(), 7u);
  EXPECT_EQ(aug.Dimension(), 4u);
  EXPECT_EQ(aug.MinimumDistance(), 3u);
  EXPECT_TRUE(aug.Contains(GFArray<GaloisFieldBinary>::Ones(gf2, {7})));

  // Non-binary ExtendedCode over GF(7): every codeword sums to 0 mod 7.
  using GF7 = GaloisFieldPrime<uint8_t>;
  auto gf7 = std::make_shared<GF7>(7);
  auto rs53 = GeneralizedReedSolomonCode<GF7>::ReedSolomon(gf7, 5, 3);
  auto rs53_ext = ExtendedCode(rs53);
  EXPECT_EQ(rs53_ext.Length(), 6u);
  EXPECT_EQ(rs53_ext.Dimension(), 3u);
  for (const auto &cw : rs53_ext.GetCodewords()) {
    uint8_t sum = 0;
    for (uint8_t v : cw.RawSpan()) {
      sum = gf7->Add(sum, v);
    }
    EXPECT_TRUE(gf7->IsZero(sum));
  }
}

TEST(PlotkinCodeTest, Rm13FromRm12AndRm02HardAndSoftDecoders) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto rm12 = std::make_shared<LinearCode<GaloisFieldBinary>>(
      ReedMullerCode(gf2, 1, 2));
  auto rm02 = std::make_shared<LinearCode<GaloisFieldBinary>>(
      ReedMullerCode(gf2, 0, 2));

  PlotkinCode<GaloisFieldBinary> rm13_plotkin(rm12, rm02);
  EXPECT_EQ(rm13_plotkin.Length(), 8u);
  EXPECT_EQ(rm13_plotkin.Dimension(), 4u);
  EXPECT_EQ(rm13_plotkin.MinimumDistance(), 4u);

  ExpectOrthogonal(rm13_plotkin.GeneratorMatrix(),
                   rm13_plotkin.ParityCheckMatrix());

  // Verify exact codeword set matches ReedMullerCode(1, 3).
  auto rm13_ref = ReedMullerCode(gf2, 1, 3);
  const auto cws = rm13_plotkin.GetCodewords();
  EXPECT_EQ(cws.size(), 16u);
  for (const auto &cw : cws) {
    EXPECT_TRUE(rm13_ref.Contains(cw));
  }

  // Verify Plotkin encoder matches GeneratorMatrix encoder & round-trips.
  std::mt19937_64 rng(555);
  auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {4}, rng);
  auto cw = rm13_plotkin.Encode(msg);
  EXPECT_EQ(cw, rm13_plotkin.Encode(msg, "GeneratorMatrix"));
  EXPECT_EQ(rm13_plotkin.Unencode(cw), msg);

  // Hard Plotkin decoder corrects every single-bit error (t = 1).
  for (std::size_t pos = 0; pos < 8; ++pos) {
    auto recv = cw;
    recv.Raw(pos) ^= uint8_t{1};
    auto res = rm13_plotkin.DecodeWithResult(recv, "PlotkinHard");
    ASSERT_TRUE(res.Ok()) << "Failed at error pos " << pos;
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 1u);
  }

  // Soft Plotkin LLR decoder (both exact boxplus and min-sum) corrects 2 low-
  // reliability bit flips where hard decoding fails!
  std::vector<double> llrs(8);
  const auto cw_span = cw.RawSpan();
  for (std::size_t i = 0; i < 8; ++i) {
    llrs[i] = (cw_span[i] == 0) ? 4.0 : -4.0;
  }
  // Flip bit 1 and bit 6 with small magnitude.
  llrs[1] = (cw_span[1] == 0) ? -0.5 : 0.5;
  llrs[6] = (cw_span[6] == 0) ? -0.5 : 0.5;

  auto soft_exact =
      rm13_plotkin.DecodeSoftPlotkin(llrs, /*exact_boxplus=*/true);
  ASSERT_TRUE(soft_exact.Ok());
  EXPECT_EQ(soft_exact.codeword, cw);
  EXPECT_EQ(soft_exact.message, msg);

  auto soft_minsum =
      rm13_plotkin.DecodeSoftPlotkin(llrs, /*exact_boxplus=*/false);
  ASSERT_TRUE(soft_minsum.Ok());
  EXPECT_EQ(soft_minsum.codeword, cw);
  EXPECT_EQ(soft_minsum.message, msg);

  // Nested 2-level PlotkinCode: RM(1, 4) = [16, 5, 8] from RM(1, 3) & RM(0, 3).
  auto rm13_ptr = std::make_shared<PlotkinCode<GaloisFieldBinary>>(rm12, rm02);
  auto rm03_ptr = std::make_shared<LinearCode<GaloisFieldBinary>>(
      ReedMullerCode(gf2, 0, 3));
  PlotkinCode<GaloisFieldBinary> rm14(rm13_ptr, rm03_ptr);
  EXPECT_EQ(rm14.Length(), 16u);
  EXPECT_EQ(rm14.Dimension(), 5u);
  EXPECT_EQ(rm14.MinimumDistance(), 8u);

  auto msg14 = GFArray<GaloisFieldBinary>::Random(gf2, {5}, rng);
  auto cw14 = rm14.Encode(msg14);
  std::vector<double> llrs14(16);
  for (std::size_t i = 0; i < 16; ++i) {
    llrs14[i] = (cw14.Raw(i) == 0) ? 3.5 : -3.5;
  }
  // Flip 3 low-reliability bits across both halves.
  for (std::size_t pos : {std::size_t{2}, std::size_t{7}, std::size_t{11}}) {
    llrs14[pos] = (cw14.Raw(pos) == 0) ? -0.6 : 0.6;
  }
  auto res14 = rm14.DecodeSoftPlotkin(llrs14);
  ASSERT_TRUE(res14.Ok());
  EXPECT_EQ(res14.codeword, cw14);
  EXPECT_EQ(res14.message, msg14);
}

TEST(PlotkinCodeTest, NonBinaryPlotkinHardDecoder) {
  using GF5 = GaloisFieldPrime<uint8_t>;
  auto gf5 = std::make_shared<GF5>(5);

  auto u_code = std::make_shared<GeneralizedReedSolomonCode<GF5>>(
      GeneralizedReedSolomonCode<GF5>::ReedSolomon(gf5, 4, 2));  // [4, 2, 3]
  auto v_code = std::make_shared<GeneralizedReedSolomonCode<GF5>>(
      GeneralizedReedSolomonCode<GF5>::ReedSolomon(gf5, 4, 1));  // [4, 1, 4]

  PlotkinCode<GF5> plotkin(u_code, v_code);
  EXPECT_EQ(plotkin.Length(), 8u);
  EXPECT_EQ(plotkin.Dimension(), 3u);
  EXPECT_EQ(plotkin.MinimumDistance(), 4u);
  ExpectOrthogonal(plotkin.GeneratorMatrix(), plotkin.ParityCheckMatrix());

  std::mt19937_64 rng(321);
  auto msg = GFArray<GF5>::Random(gf5, {3}, rng);
  auto cw = plotkin.Encode(msg);
  EXPECT_EQ(cw, plotkin.Encode(msg, "GeneratorMatrix"));

  // Single symbol error at any position is corrected by PlotkinHard.
  for (std::size_t pos = 0; pos < 8; ++pos) {
    auto recv = cw;
    recv.Raw(pos) = gf5->Add(recv.Raw(pos), gf5->SetElementValue(2));
    auto res = plotkin.DecodeWithResult(recv, "PlotkinHard");
    ASSERT_TRUE(res.Ok()) << "Failed at pos " << pos;
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 1u);
  }
}

// ============================================================================
// 4. WeightEnumerator & Exact MacWilliams Identity Tests
// ============================================================================

TEST(WeightEnumeratorTest, ExtendedBinaryGolay24128SelfDuality) {
  auto golay = ExtendedBinaryGolayCode();
  auto dist = ComputeWeightEnumerator(golay);

  EXPECT_EQ(dist.length, 24u);
  EXPECT_EQ(dist.dimension, 12u);
  EXPECT_EQ(dist.field_order, 2u);
  EXPECT_EQ(dist.MinimumDistance(), 8u);

  // Exact weight enumerator: 1 + 759*z^8 + 2576*z^12 + 759*z^16 + z^24.
  for (std::size_t w = 0; w <= 24; ++w) {
    if (w == 0 || w == 24) {
      EXPECT_EQ(dist.counts[w], BigInt(1)) << "w=" << w;
    } else if (w == 8 || w == 16) {
      EXPECT_EQ(dist.counts[w], BigInt(759)) << "w=" << w;
    } else if (w == 12) {
      EXPECT_EQ(dist.counts[w], BigInt(2576)) << "w=" << w;
    } else {
      EXPECT_EQ(dist.counts[w], BigInt(0)) << "w=" << w;
    }
  }
  EXPECT_EQ(dist.ToPolynomialString("z"),
            "1 + 759*z^8 + 2576*z^12 + 759*z^16 + z^24");

  // Self-duality: MacWilliamsTransform(A) == A exactly!
  auto dual_dist = MacWilliamsTransform(dist);
  EXPECT_EQ(dual_dist, dist);
}

TEST(WeightEnumeratorTest, Hamming743AndSimplex734MacWilliamsRoundTrip) {
  auto ham = HammingCode(3);  // [7, 4, 3]: k > n - k -> uses dual + MacWilliams
  auto simplex = ham.DualCode();  // [7, 3, 4]: k <= n - k -> direct enumeration

  auto ham_dist = ComputeWeightEnumerator(ham);
  auto simp_dist = ComputeWeightEnumerator(simplex);

  EXPECT_EQ(ham_dist.MinimumDistance(), 3u);
  EXPECT_EQ(ham_dist.ToPolynomialString("z"), "1 + 7*z^3 + 7*z^4 + z^7");

  EXPECT_EQ(simp_dist.MinimumDistance(), 4u);
  EXPECT_EQ(simp_dist.ToPolynomialString("z"), "1 + 7*z^4");

  EXPECT_EQ(MacWilliamsTransform(simp_dist), ham_dist);
  EXPECT_EQ(MacWilliamsTransform(ham_dist), simp_dist);
}

TEST(WeightEnumeratorTest, NonBinaryMdsRs744OverGf8ClosedFormAndDual) {
  using GF8 = GaloisFieldBinaryExtension<uint8_t>;
  auto gf8 = std::make_shared<GF8>(3);

  // RS(7, 4, 4) over GF(8): k = 4 > n - k = 3, so ComputeWeightEnumerator
  // enumerates the dual RS(7, 3, 5) (8^3 = 512 words) and applies MacWilliams.
  auto rs74 = GeneralizedReedSolomonCode<GF8>::ReedSolomon(gf8, 7, 4);
  auto dist = ComputeWeightEnumerator(rs74);

  EXPECT_EQ(dist.length, 7u);
  EXPECT_EQ(dist.dimension, 4u);
  EXPECT_EQ(dist.field_order, 8u);
  EXPECT_EQ(dist.MinimumDistance(), 4u);

  BigInt total(0);
  for (std::size_t w = 0; w <= 7; ++w) {
    const BigInt expected = ExactMdsWeightCount(7, 4, 8, w);
    EXPECT_EQ(dist.counts[w], expected) << "Mismatch at weight w=" << w;
    total += dist.counts[w];
  }
  EXPECT_EQ(total, BigInt::FromUint64(8).Pow(4));

  // Verify MacWilliamsTransform of RS(7, 4) matches RS(7, 4).DualCode().
  auto dual_rs = rs74.DualCode();
  auto dual_dist_direct = ComputeWeightEnumerator(dual_rs);
  auto dual_dist_mw = MacWilliamsTransform(dist);
  EXPECT_EQ(dual_dist_mw, dual_dist_direct);
  EXPECT_EQ(MacWilliamsTransform(dual_dist_direct), dist);
}

TEST(WeightEnumeratorTest, MacWilliamsRejectsInvalidDistribution) {
  WeightDistribution bogus;
  bogus.length = 7;
  bogus.dimension = 3;
  bogus.field_order = 2;
  bogus.counts = {BigInt(1), BigInt(1), BigInt(0), BigInt(0),
                  BigInt(0), BigInt(0), BigInt(0), BigInt(0)};
  EXPECT_THROW((void)MacWilliamsTransform(bogus), std::runtime_error);
}

}  // namespace
