#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/ldpc.hpp"
#include "xgalois/coding/linear_code.hpp"
#include "xgalois/coding/polar.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace {

using namespace xg;
using namespace xg::coding;
using namespace xg::channel;

GFArray<GaloisFieldBinary> MakeBinaryWord(
    const std::shared_ptr<GaloisFieldBinary> &gf2,
    const std::vector<uint8_t> &bits) {
  using RawArray = GFArray<GaloisFieldBinary>::raw_array_type;
  using ShapeType = GFArray<GaloisFieldBinary>::shape_type;
  RawArray raw = RawArray::from_shape(ShapeType{bits.size()});
  for (std::size_t i = 0; i < bits.size(); ++i) {
    raw(i) = static_cast<uint8_t>(bits[i] & 1U);
  }
  return GFArray<GaloisFieldBinary>::FromRaw(gf2, std::move(raw));
}

bool IsAllZero(const GFArray<GaloisFieldBinary> &arr) {
  for (uint8_t v : arr.RawSpan()) {
    if ((v & 1U) != 0) {
      return false;
    }
  }
  return true;
}

bool ArraysEqual(const GFArray<GaloisFieldBinary> &a,
                 const GFArray<GaloisFieldBinary> &b) {
  if (!a.HasField() || !b.HasField() || a.Size() != b.Size() ||
      a.Dimension() != b.Dimension()) {
    return false;
  }
  const auto sa = a.RawSpan();
  const auto sb = b.RawSpan();
  for (std::size_t i = 0; i < sa.size(); ++i) {
    if ((sa[i] & 1U) != (sb[i] & 1U)) {
      return false;
    }
  }
  return true;
}

std::size_t CountBitDiffs(const GFArray<GaloisFieldBinary> &a,
                          const GFArray<GaloisFieldBinary> &b) {
  std::size_t diffs = 0;
  const auto sa = a.RawSpan();
  const auto sb = b.RawSpan();
  for (std::size_t i = 0; i < sa.size(); ++i) {
    if ((sa[i] & 1U) != (sb[i] & 1U)) {
      ++diffs;
    }
  }
  return diffs;
}

// Construct a girth >= 6 regular (3, 6) QC-LDPC base matrix with prime lifting
// factor Z = 17: B(i, j) = (i * j) mod 17 has zero 4-cycles because 17 is prime
// and (i1 - i2) * (j1 - j2) != 0 (mod 17) for distinct i1,i2 in [0,2] and
// j1,j2 in [0,5].
std::vector<std::vector<int>> MakeGirth6QcBaseMatrix() {
  return {
      {0, 0, 0, 0, 0, 0},
      {0, 1, 2, 3, 4, 5},
      {0, 2, 4, 6, 8, 10},
  };
}

// ============================================================================
// LDPC Tests
// ============================================================================

TEST(LdpcTest, SparseBinaryMatrixRepresentationsAndSyndrome) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // Base matrix with -1 (zero block) and non-negative shifts.
  const std::vector<std::vector<int>> base = {
      {0, 1, -1, 2},
      {-1, 0, 3, 1},
  };
  const std::size_t z = 4;
  SparseBinaryMatrix sparse = SparseBinaryMatrix::FromQuasiCyclic(base, z);
  EXPECT_EQ(sparse.NumRows(), 8U);
  EXPECT_EQ(sparse.NumCols(), 16U);
  EXPECT_EQ(sparse.NumNonZeros(), 6U * z);

  for (std::size_t r = 0; r < sparse.NumRows(); ++r) {
    EXPECT_EQ(sparse.RowDegree(r), 3U);
  }
  // Col blocks 0 and 2 have degree 1; col blocks 1 and 3 have degree 2.
  for (std::size_t c = 0; c < 4; ++c) {
    EXPECT_EQ(sparse.ColDegree(c), 1U);
    EXPECT_EQ(sparse.ColDegree(4 + c), 2U);
    EXPECT_EQ(sparse.ColDegree(8 + c), 1U);
    EXPECT_EQ(sparse.ColDegree(12 + c), 2U);
  }

  auto dense = sparse.ToDense(gf2);
  SparseBinaryMatrix roundtrip = SparseBinaryMatrix::FromDense(dense);
  EXPECT_EQ(roundtrip.NumRows(), sparse.NumRows());
  EXPECT_EQ(roundtrip.NumCols(), sparse.NumCols());
  EXPECT_EQ(roundtrip.NumNonZeros(), sparse.NumNonZeros());

  std::mt19937_64 rng(123);
  for (int trial = 0; trial < 20; ++trial) {
    auto word = GFArray<GaloisFieldBinary>::Random(gf2, {16}, rng);
    auto syn_sparse = sparse.Syndrome(word);
    auto syn_dense = linalg::Dot(dense, word);
    EXPECT_TRUE(ArraysEqual(syn_sparse, syn_dense));
    EXPECT_EQ(sparse.CheckParity(word), IsAllZero(syn_dense));
  }
}

TEST(LdpcTest, AlgebraicInvariantsHammingRegularAndQuasiCyclic) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // 1. [7, 4] Hamming code represented as an LDPC code.
  auto hamming = HammingCode(gf2, 3);
  LDPCCode hamming_ldpc(gf2, hamming.ParityCheckMatrix());
  EXPECT_EQ(hamming_ldpc.Length(), 7U);
  EXPECT_EQ(hamming_ldpc.Dimension(), 4U);

  // 2. Regular Gallager (48, 3, 6) LDPC code.
  auto reg_ldpc = LDPCCode::MakeRegularLDPC(gf2, 48, 3, 6, 42);
  EXPECT_EQ(reg_ldpc.Length(), 48U);
  EXPECT_GT(reg_ldpc.Dimension(), 20U);
  for (std::size_t r = 0; r < reg_ldpc.SparseParityCheckMatrix().NumRows();
       ++r) {
    EXPECT_EQ(reg_ldpc.SparseParityCheckMatrix().RowDegree(r), 6U);
  }
  for (std::size_t c = 0; c < reg_ldpc.SparseParityCheckMatrix().NumCols();
       ++c) {
    EXPECT_EQ(reg_ldpc.SparseParityCheckMatrix().ColDegree(c), 3U);
  }

  // 3. Girth-6 Quasi-Cyclic (102, 53) LDPC code.
  auto qc_ldpc = LDPCCode::FromQuasiCyclic(gf2, MakeGirth6QcBaseMatrix(), 17);
  EXPECT_EQ(qc_ldpc.Length(), 102U);
  EXPECT_EQ(qc_ldpc.Dimension(), 53U);

  std::mt19937_64 rng(2026);
  const std::vector<const LDPCCode *> codes = {&hamming_ldpc, &reg_ldpc,
                                               &qc_ldpc};
  for (const LDPCCode *code : codes) {
    const auto &G = code->GeneratorMatrix();
    const auto &H = code->ParityCheckMatrix();
    EXPECT_EQ(G.Shape(0), code->Dimension());
    EXPECT_EQ(G.Shape(1), code->Length());
    EXPECT_EQ(H.Shape(1), code->Length());

    // Exact orthogonality G * H^T = 0.
    auto ght = linalg::Dot(G, linalg::Transpose(H));
    EXPECT_TRUE(IsAllZero(ght));

    // All-zero message round-trip.
    auto zero_msg = GFArray<GaloisFieldBinary>::Zeros(gf2, {code->Dimension()});
    auto zero_cw = code->Encode(zero_msg);
    EXPECT_TRUE(IsAllZero(zero_cw));
    EXPECT_TRUE(code->Contains(zero_cw));
    EXPECT_TRUE(code->CheckParity(zero_cw));
    EXPECT_TRUE(ArraysEqual(code->Unencode(zero_cw), zero_msg));

    // Random messages: Encode/Unencode, Contains, Sparse Syndrome == Dense.
    for (int t = 0; t < 15; ++t) {
      auto msg =
          GFArray<GaloisFieldBinary>::Random(gf2, {code->Dimension()}, rng);
      auto cw = code->Encode(msg);
      EXPECT_TRUE(code->Contains(cw));
      EXPECT_TRUE(code->CheckParity(cw));
      EXPECT_TRUE(IsAllZero(code->SyndromeSparse(cw)));
      EXPECT_TRUE(IsAllZero(code->Syndrome(cw)));
      EXPECT_TRUE(ArraysEqual(code->Unencode(cw), msg));

      // Perturb 1 bit and check sparse syndrome matches dense syndrome.
      auto perturbed = cw;
      const std::size_t flip_pos = static_cast<std::size_t>(t) % code->Length();
      perturbed.Raw(flip_pos) ^= 1U;
      EXPECT_TRUE(ArraysEqual(code->SyndromeSparse(perturbed),
                              code->Syndrome(perturbed)));
    }
  }
}

TEST(LdpcTest, BitFlippingDecoderCorrectsErrorsAndDetectsFailure) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // [7, 4] Hamming-as-LDPC: must correct every single-bit error.
  auto hamming = HammingCode(gf2, 3);
  LDPCCode hamming_ldpc(gf2, hamming.ParityCheckMatrix());
  auto msg = MakeBinaryWord(gf2, {1, 0, 1, 1});
  auto cw = hamming_ldpc.Encode(msg);
  for (std::size_t pos = 0; pos < 7; ++pos) {
    auto rx = cw;
    rx.Raw(pos) ^= 1U;
    auto res = hamming_ldpc.DecodeWithResult(rx, "BitFlipping");
    ASSERT_TRUE(res.Ok()) << "Failed on Hamming error position " << pos;
    EXPECT_EQ(res.num_errors, 1U);
    EXPECT_TRUE(ArraysEqual(res.codeword, cw));
    EXPECT_TRUE(ArraysEqual(res.message, msg));
  }

  // Girth-6 QC-LDPC (102, 53): corrects 1-bit and 2-bit hard errors.
  auto qc_ldpc = LDPCCode::FromQuasiCyclic(gf2, MakeGirth6QcBaseMatrix(), 17);
  BitFlippingDecoder bf(qc_ldpc, 50);
  std::mt19937_64 rng(777);
  for (int trial = 0; trial < 25; ++trial) {
    auto m =
        GFArray<GaloisFieldBinary>::Random(gf2, {qc_ldpc.Dimension()}, rng);
    auto c = qc_ldpc.Encode(m);
    auto rx = c;
    const std::size_t p1 = rng() % qc_ldpc.Length();
    std::size_t p2 = rng() % qc_ldpc.Length();
    while (p2 == p1) {
      p2 = rng() % qc_ldpc.Length();
    }
    rx.Raw(p1) ^= 1U;
    rx.Raw(p2) ^= 1U;
    auto res = bf.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.num_errors, 2U);
    EXPECT_TRUE(ArraysEqual(res.codeword, c));
    EXPECT_TRUE(ArraysEqual(res.message, m));
  }

  // Uncorrectable word with massive corruption should fail cleanly (not claim
  // success on a non-codeword).
  auto zero_cw = GFArray<GaloisFieldBinary>::Zeros(gf2, {qc_ldpc.Length()});
  auto heavily_corrupted = zero_cw;
  for (std::size_t i = 0; i < 40; ++i) {
    heavily_corrupted.Raw(i) = 1U;
  }
  auto fail_res = bf.DecodeWithResult(heavily_corrupted);
  if (fail_res.Ok()) {
    EXPECT_TRUE(qc_ldpc.CheckParity(fail_res.codeword));
  } else {
    EXPECT_FALSE(fail_res.error_message.empty());
  }
}

TEST(LdpcTest, BeliefPropagationAwgnSumProductAndMinSumVariants) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto qc_ldpc = LDPCCode::FromQuasiCyclic(gf2, MakeGirth6QcBaseMatrix(), 17);
  const double rate = qc_ldpc.Rate();

  BeliefPropagationDecoder sp_dec(qc_ldpc, BpAlgorithm::kSumProduct, 50);
  BeliefPropagationDecoder nms_dec(qc_ldpc, BpAlgorithm::kNormalizedMinSum, 50,
                                   0.75);
  BeliefPropagationDecoder oms_dec(qc_ldpc, BpAlgorithm::kOffsetMinSum, 50, 1.0,
                                   0.4);

  // 1. Noiseless channel (infinite LLRs) must decode cleanly without NaN.
  {
    auto ch_noiseless = AWGNChannel::FromSigma(gf2, qc_ldpc.Length(), 0.0, 11);
    std::mt19937_64 rng(11);
    auto msg =
        GFArray<GaloisFieldBinary>::Random(gf2, {qc_ldpc.Dimension()}, rng);
    auto cw = qc_ldpc.Encode(msg);
    auto llrs = ch_noiseless.TransmitLlr(cw);
    for (const auto *dec : {&sp_dec, &nms_dec, &oms_dec}) {
      auto res = dec->DecodeLlr(llrs);
      ASSERT_TRUE(res.Ok());
      EXPECT_EQ(res.num_errors, 0U);
      EXPECT_TRUE(ArraysEqual(res.codeword, cw));
      EXPECT_TRUE(ArraysEqual(res.message, msg));
    }
  }

  // 2. High SNR (7.0 dB): near-zero FER across all three BP variants.
  {
    auto ch_high =
        AWGNChannel::FromEbN0(gf2, qc_ldpc.Length(), 7.0, rate, 1001);
    std::mt19937_64 rng(1002);
    constexpr int kHighFrames = 40;
    int sp_err = 0;
    int nms_err = 0;
    int oms_err = 0;
    for (int f = 0; f < kHighFrames; ++f) {
      auto msg =
          GFArray<GaloisFieldBinary>::Random(gf2, {qc_ldpc.Dimension()}, rng);
      auto cw = qc_ldpc.Encode(msg);
      auto tx = ch_high.TransmitDetailed(cw);

      auto r_sp = sp_dec.DecodeLlr(tx.llr);
      auto r_nms = nms_dec.DecodeLlr(tx.llr);
      auto r_oms = oms_dec.DecodeLlr(tx.llr);

      if (!r_sp.Ok() || !ArraysEqual(r_sp.message, msg)) ++sp_err;
      if (!r_nms.Ok() || !ArraysEqual(r_nms.message, msg)) ++nms_err;
      if (!r_oms.Ok() || !ArraysEqual(r_oms.message, msg)) ++oms_err;
    }
    EXPECT_LE(sp_err, 1);
    EXPECT_LE(nms_err, 1);
    EXPECT_LE(oms_err, 1);
  }

  // 3. Moderate SNR (3.5 dB): Sum-Product, Normalized Min-Sum, and Offset
  // Min-Sum must all substantially outperform raw uncoded hard decisions.
  {
    auto ch_mod = AWGNChannel::FromEbN0(gf2, qc_ldpc.Length(), 3.5, rate, 2001);
    std::mt19937_64 rng(2002);
    constexpr int kModFrames = 80;
    int raw_frame_errors = 0;
    int sp_frame_errors = 0;
    int nms_frame_errors = 0;
    int oms_frame_errors = 0;

    for (int f = 0; f < kModFrames; ++f) {
      auto msg =
          GFArray<GaloisFieldBinary>::Random(gf2, {qc_ldpc.Dimension()}, rng);
      auto cw = qc_ldpc.Encode(msg);
      auto tx = ch_mod.TransmitDetailed(cw);

      if (!ArraysEqual(tx.hard_decision, cw)) {
        ++raw_frame_errors;
      }

      auto r_sp = sp_dec.DecodeLlr(tx.llr);
      if (r_sp.Ok()) {
        EXPECT_TRUE(qc_ldpc.CheckParity(r_sp.codeword));
        EXPECT_EQ(r_sp.num_errors,
                  CountBitDiffs(r_sp.codeword, tx.hard_decision));
      }
      if (!r_sp.Ok() || !ArraysEqual(r_sp.message, msg)) {
        ++sp_frame_errors;
      }

      auto r_nms = nms_dec.DecodeLlr(tx.llr);
      if (!r_nms.Ok() || !ArraysEqual(r_nms.message, msg)) {
        ++nms_frame_errors;
      }

      auto r_oms = oms_dec.DecodeLlr(tx.llr);
      if (!r_oms.Ok() || !ArraysEqual(r_oms.message, msg)) {
        ++oms_frame_errors;
      }
    }

    EXPECT_GT(raw_frame_errors, 25);
    EXPECT_LT(sp_frame_errors, raw_frame_errors / 2);
    EXPECT_LT(nms_frame_errors, raw_frame_errors / 2);
    EXPECT_LT(oms_frame_errors, raw_frame_errors / 2);
  }

  // 4. Hard-decision Decoder interface on BeliefPropagationDecoder.
  {
    std::mt19937_64 rng(3003);
    auto msg =
        GFArray<GaloisFieldBinary>::Random(gf2, {qc_ldpc.Dimension()}, rng);
    auto cw = qc_ldpc.Encode(msg);
    auto rx = cw;
    rx.Raw(5) ^= 1U;
    rx.Raw(42) ^= 1U;
    auto res = qc_ldpc.DecodeWithResult(rx, "BeliefPropagation");
    ASSERT_TRUE(res.Ok());
    EXPECT_EQ(res.num_errors, 2U);
    EXPECT_TRUE(ArraysEqual(res.message, msg));
  }
}

// ============================================================================
// Polar Code Tests
// ============================================================================

TEST(PolarTest, TransformMatchesDenseGeneratorAndOrthogonality) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  std::mt19937_64 rng(99);

  for (std::size_t n :
       {std::size_t{8}, std::size_t{16}, std::size_t{32}, std::size_t{64}}) {
    const std::size_t k = n / 2;

    // Rate-1 full transform check: PolarTransform(u) == u * F^{\otimes
    // log2(n)}.
    std::vector<bool> no_frozen(n, false);
    PolarCode full_polar(gf2, n, no_frozen, /*systematic=*/false);
    const auto &G_full = full_polar.GeneratorMatrix();
    for (int t = 0; t < 10; ++t) {
      auto u = GFArray<GaloisFieldBinary>::Random(gf2, {n}, rng);
      auto x_fast = PolarCode::PolarTransform(u);
      auto x_dense = linalg::Dot(u, G_full);
      EXPECT_TRUE(ArraysEqual(x_fast, x_dense));
      // Involution: PolarTransform(PolarTransform(u)) == u.
      EXPECT_TRUE(ArraysEqual(PolarCode::PolarTransform(x_fast), u));
    }

    // Half-rate polar codes (both Bhattacharyya and Design-SNR constructions).
    auto polar_bhatt = PolarCode::FromBhattacharyya(gf2, n, k, 0.5, false);
    auto polar_ga = PolarCode::FromDesignSnr(gf2, n, k, 2.0, false);

    for (const PolarCode *code : {&polar_bhatt, &polar_ga}) {
      EXPECT_EQ(code->Length(), n);
      EXPECT_EQ(code->Dimension(), k);
      EXPECT_EQ(code->InfoPositions().size(), k);
      EXPECT_EQ(code->FrozenPositions().size(), n - k);
      EXPECT_TRUE(code->FrozenMask()[0]);
      EXPECT_FALSE(code->FrozenMask()[n - 1]);

      // Exact G * H^T = 0.
      auto ght = linalg::Dot(code->GeneratorMatrix(),
                             linalg::Transpose(code->ParityCheckMatrix()));
      EXPECT_TRUE(IsAllZero(ght));

      // Verify O(N log N) Polar encoder matches dense GeneratorMatrix encoder.
      for (int t = 0; t < 10; ++t) {
        auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {k}, rng);
        auto cw_polar = code->Encode(msg, "Polar");
        auto cw_dense = code->Encode(msg, "GeneratorMatrix");
        EXPECT_TRUE(ArraysEqual(cw_polar, cw_dense));
        EXPECT_TRUE(code->Contains(cw_polar));
        EXPECT_TRUE(ArraysEqual(code->Unencode(cw_polar, "Polar"), msg));
        EXPECT_TRUE(
            ArraysEqual(code->Unencode(cw_polar, "GeneratorMatrix"), msg));
      }
    }
  }
}

TEST(PolarTest, NonSystematicAndSystematicRoundTrip) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  std::mt19937_64 rng(555);

  // Standard Bhattacharyya & Design-SNR systematic codes (upper-set fast path).
  auto sys_bhatt = PolarCode::FromBhattacharyya(gf2, 32, 16, 0.5, true);
  auto sys_ga = PolarCode::FromDesignSnr(gf2, 64, 32, 2.0, true);

  // Custom non-upper-set frozen mask (exercises general lower-triangular
  // solve). N = 8, info = {1, 5, 7}, frozen = {0, 2, 3, 4, 6} (3 is between 1
  // and 7).
  std::vector<bool> non_upper_mask = {true, false, true, true,
                                      true, false, true, false};
  PolarCode sys_custom(gf2, 8, non_upper_mask, /*systematic=*/true);

  for (const PolarCode *code : {&sys_bhatt, &sys_ga, &sys_custom}) {
    EXPECT_TRUE(code->IsSystematicPolar());
    auto ght = linalg::Dot(code->GeneratorMatrix(),
                           linalg::Transpose(code->ParityCheckMatrix()));
    EXPECT_TRUE(IsAllZero(ght));

    // All-zero message.
    auto zero_msg = GFArray<GaloisFieldBinary>::Zeros(gf2, {code->Dimension()});
    auto zero_cw = code->Encode(zero_msg);
    EXPECT_TRUE(IsAllZero(zero_cw));
    EXPECT_TRUE(code->Contains(zero_cw));
    EXPECT_TRUE(ArraysEqual(code->Unencode(zero_cw), zero_msg));

    // Random messages: systematic bits sit directly at InfoPositions().
    for (int t = 0; t < 20; ++t) {
      auto msg =
          GFArray<GaloisFieldBinary>::Random(gf2, {code->Dimension()}, rng);
      auto cw = code->Encode(msg, "Polar");
      auto cw_dense = code->Encode(msg, "GeneratorMatrix");
      EXPECT_TRUE(ArraysEqual(cw, cw_dense));
      EXPECT_TRUE(code->Contains(cw));

      const auto &info_pos = code->InfoPositions();
      for (std::size_t r = 0; r < code->Dimension(); ++r) {
        EXPECT_EQ(cw.Raw(info_pos[r]), msg.Raw(r));
      }
      EXPECT_TRUE(ArraysEqual(code->Unencode(cw), msg));
    }
  }
}

TEST(PolarTest, ScAndSclOverAwgnPolar64And128) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // 1. Polar(64, 32) non-systematic and systematic: noiseless + high-SNR +
  // moderate-SNR comparison of SC vs SCL(L=1, 4, 8).
  for (bool systematic : {false, true}) {
    auto polar64 = PolarCode::FromDesignSnr(gf2, 64, 32, 2.5, systematic);
    SuccessiveCancellationDecoder sc(polar64);
    SuccessiveCancellationListDecoder scl1(polar64, 1);
    SuccessiveCancellationListDecoder scl4(polar64, 4);
    SuccessiveCancellationListDecoder scl8(polar64, 8);

    // Noiseless (sigma = 0, infinite LLRs) & high-SNR exact recovery.
    auto ch_zero = AWGNChannel::FromSigma(gf2, 64, 0.0, 401);
    std::mt19937_64 rng(402);
    for (int f = 0; f < 10; ++f) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {32}, rng);
      auto cw = polar64.Encode(msg);
      auto llrs = ch_zero.TransmitLlr(cw);
      auto r_sc = sc.DecodeLlr(llrs);
      auto r_scl = scl8.DecodeLlr(llrs);
      ASSERT_TRUE(r_sc.Ok());
      ASSERT_TRUE(r_scl.Ok());
      EXPECT_TRUE(ArraysEqual(r_sc.message, msg));
      EXPECT_TRUE(ArraysEqual(r_scl.message, msg));
    }

    // Moderate SNR (3.0 dB): SCL(L=1) matches SC identically, and SCL(L=8)
    // achieves FER <= SC on identical seeded AWGN frames.
    auto ch_mod = AWGNChannel::FromEbN0(gf2, 64, 3.0, 0.5, 501);
    int sc_errors = 0;
    int scl4_errors = 0;
    int scl8_errors = 0;
    int raw_errors = 0;
    constexpr int kFrames = 80;
    for (int f = 0; f < kFrames; ++f) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {32}, rng);
      auto cw = polar64.Encode(msg);
      auto tx = ch_mod.TransmitDetailed(cw);

      if (!ArraysEqual(tx.hard_decision, cw)) {
        ++raw_errors;
      }

      auto r_sc = sc.DecodeLlr(tx.llr);
      auto r_scl1 = scl1.DecodeLlr(tx.llr);
      auto r_scl4 = scl4.DecodeLlr(tx.llr);
      auto r_scl8 = scl8.DecodeLlr(tx.llr);

      ASSERT_TRUE(r_sc.Ok());
      ASSERT_TRUE(r_scl1.Ok());
      ASSERT_TRUE(r_scl4.Ok());
      ASSERT_TRUE(r_scl8.Ok());
      EXPECT_TRUE(ArraysEqual(r_sc.codeword, r_scl1.codeword));
      EXPECT_TRUE(polar64.Contains(r_scl8.codeword));

      if (!ArraysEqual(r_sc.message, msg)) ++sc_errors;
      if (!ArraysEqual(r_scl4.message, msg)) ++scl4_errors;
      if (!ArraysEqual(r_scl8.message, msg)) ++scl8_errors;
    }
    EXPECT_GT(raw_errors, 15);
    EXPECT_LE(scl8_errors, sc_errors);
    EXPECT_LE(scl4_errors, raw_errors);
  }

  // 2. Polar(128, 64) over AWGN at high SNR (5.5 dB) and moderate SNR (2.5 dB).
  {
    auto polar128 = PolarCode::FromDesignSnr(gf2, 128, 64, 2.5, false);
    SuccessiveCancellationDecoder sc128(polar128);
    SuccessiveCancellationListDecoder scl128(polar128, 8);

    auto ch_high = AWGNChannel::FromEbN0(gf2, 128, 5.5, 0.5, 601);
    std::mt19937_64 rng(602);
    int high_err = 0;
    for (int f = 0; f < 30; ++f) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {64}, rng);
      auto cw = polar128.Encode(msg);
      auto llrs = ch_high.TransmitLlr(cw);
      auto r_scl = scl128.DecodeLlr(llrs);
      if (!r_scl.Ok() || !ArraysEqual(r_scl.message, msg)) {
        ++high_err;
      }
    }
    EXPECT_EQ(high_err, 0);

    auto ch_mod = AWGNChannel::FromEbN0(gf2, 128, 2.5, 0.5, 603);
    int sc_err = 0;
    int scl_err = 0;
    for (int f = 0; f < 50; ++f) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {64}, rng);
      auto cw = polar128.Encode(msg);
      auto llrs = ch_mod.TransmitLlr(cw);
      auto r_sc = sc128.DecodeLlr(llrs);
      auto r_scl = scl128.DecodeLlr(llrs);
      if (!r_sc.Ok() || !ArraysEqual(r_sc.message, msg)) ++sc_err;
      if (!r_scl.Ok() || !ArraysEqual(r_scl.message, msg)) ++scl_err;
    }
    EXPECT_LE(scl_err, sc_err);
  }
}

TEST(PolarTest, CrcAidedPolarCaSclImprovesOverPlainSc) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // Polar(64, 32) plain vs Polar(64, 32) with CRC-8 (40 non-frozen bits: 32
  // message bits + 8 CRC bits, true linear dimension K = 32).
  auto plain_polar = PolarCode::FromDesignSnr(gf2, 64, 32, 2.5, false, 0);
  auto ca_polar = PolarCode::FromDesignSnr(gf2, 64, 32, 2.5, false, 8, 0x107);
  auto ca_polar_sys =
      PolarCode::FromDesignSnr(gf2, 64, 32, 2.5, true, 8, 0x107);

  EXPECT_EQ(ca_polar.Length(), 64U);
  EXPECT_EQ(ca_polar.Dimension(), 32U);
  EXPECT_EQ(ca_polar.CrcBits(), 8U);
  EXPECT_EQ(ca_polar.InfoPositions().size(), 40U);

  // Verify G * H^T = 0 and fast vs dense encoding on CA-Polar (both modes).
  std::mt19937_64 rng(808);
  for (const PolarCode *code : {&ca_polar, &ca_polar_sys}) {
    auto ght = linalg::Dot(code->GeneratorMatrix(),
                           linalg::Transpose(code->ParityCheckMatrix()));
    EXPECT_TRUE(IsAllZero(ght));
    for (int t = 0; t < 15; ++t) {
      auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {32}, rng);
      auto cw_fast = code->Encode(msg, "Polar");
      auto cw_dense = code->Encode(msg, "GeneratorMatrix");
      EXPECT_TRUE(ArraysEqual(cw_fast, cw_dense));
      EXPECT_TRUE(code->Contains(cw_fast));
      EXPECT_TRUE(ArraysEqual(code->Unencode(cw_fast), msg));
    }
  }

  // Compare plain SC vs CA-SCL(L=8) over identical AWGN noise realizations at
  // moderate SNR (2.2 dB) where SC experiences noticeable frame errors while
  // CA-SCL recovers substantially more frames via CRC path selection.
  SuccessiveCancellationDecoder plain_sc(plain_polar);
  SuccessiveCancellationListDecoder ca_scl(ca_polar, 8);

  constexpr int kFrames = 120;
  int plain_sc_errors = 0;
  int ca_scl_errors = 0;
  int ca_scl_crc_passes = 0;

  // Use identical channel seed for fair comparison at the same Eb/N0.
  auto ch_plain = AWGNChannel::FromEbN0(gf2, 64, 2.2, 0.5, 9001);
  auto ch_ca = AWGNChannel::FromEbN0(gf2, 64, 2.2, 0.5, 9001);

  for (int f = 0; f < kFrames; ++f) {
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {32}, rng);

    auto cw_plain = plain_polar.Encode(msg);
    auto llr_plain = ch_plain.TransmitLlr(cw_plain);
    auto res_plain = plain_sc.DecodeLlr(llr_plain);
    if (!res_plain.Ok() || !ArraysEqual(res_plain.message, msg)) {
      ++plain_sc_errors;
    }

    auto cw_ca = ca_polar.Encode(msg);
    auto llr_ca = ch_ca.TransmitLlr(cw_ca);
    auto res_ca = ca_scl.DecodeLlr(llr_ca);
    if (res_ca.Ok()) {
      ++ca_scl_crc_passes;
      EXPECT_TRUE(ca_polar.Contains(res_ca.codeword));
    }
    if (!res_ca.Ok() || !ArraysEqual(res_ca.message, msg)) {
      ++ca_scl_errors;
    }
  }

  EXPECT_GT(plain_sc_errors, 10);
  EXPECT_GT(ca_scl_crc_passes, kFrames / 2);
  EXPECT_LT(ca_scl_errors, plain_sc_errors);
}

TEST(LdpcPolarValidationTest, InvalidParametersThrowCleanly) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // LDPC invalid constructions.
  EXPECT_THROW(SparseBinaryMatrix::FromQuasiCyclic({{0, 1}}, 0),
               std::invalid_argument);
  EXPECT_THROW(SparseBinaryMatrix::FromQuasiCyclic({}, 4),
               std::invalid_argument);
  EXPECT_THROW(LDPCCode::MakeRegularLDPC(gf2, 10, 3, 4), std::invalid_argument);
  EXPECT_THROW(LDPCCode::MakeRegularLDPC(gf2, 12, 4, 4), std::invalid_argument);

  auto qc = LDPCCode::FromQuasiCyclic(gf2, {{0, 1, 2}}, 4);
  EXPECT_THROW(
      BeliefPropagationDecoder(qc, BpAlgorithm::kNormalizedMinSum, 10, 1.5),
      std::invalid_argument);
  EXPECT_THROW(
      BeliefPropagationDecoder(qc, BpAlgorithm::kOffsetMinSum, 10, 1.0, -0.1),
      std::invalid_argument);
  BeliefPropagationDecoder bp(qc);
  EXPECT_THROW(bp.DecodeLlr(std::vector<double>(5, 1.0)),
               std::invalid_argument);
  std::vector<double> nan_llrs(qc.Length(),
                               std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(bp.DecodeLlr(nan_llrs), std::invalid_argument);

  // Polar invalid constructions.
  EXPECT_THROW(PolarCode::ConstructBhattacharyya(6, 3), std::invalid_argument);
  EXPECT_THROW(PolarCode::ConstructBhattacharyya(8, 0), std::invalid_argument);
  EXPECT_THROW(PolarCode::ConstructBhattacharyya(8, 9), std::invalid_argument);
  EXPECT_THROW(PolarCode::ConstructBhattacharyya(8, 4, 0.0),
               std::invalid_argument);
  EXPECT_THROW(PolarCode::FromDesignSnr(gf2, 16, 12, 2.0, false, 8),
               std::invalid_argument);

  auto p16 = PolarCode::FromDesignSnr(gf2, 16, 8);
  EXPECT_THROW(SuccessiveCancellationListDecoder(p16, 0),
               std::invalid_argument);
  SuccessiveCancellationDecoder sc(p16);
  EXPECT_THROW(sc.DecodeLlr(std::vector<double>(10, 1.0)),
               std::invalid_argument);
  std::vector<double> polar_nan(16, std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(sc.DecodeLlr(polar_nan), std::invalid_argument);
}

}  // namespace
