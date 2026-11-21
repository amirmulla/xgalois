#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
#include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
#include "xgalois/coding/ldpc.hpp"
#include "xgalois/coding/polar.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"

using namespace xg;
using namespace xg::coding;

namespace {

template <typename Field>
void PrintArray(const std::string& label, const GFArray<Field>& arr) {
  std::cout << label << ": [";
  for (size_t i = 0; i < arr.Size(); ++i) {
    std::cout << arr.At(i);
    if (i + 1 < arr.Size()) std::cout << ", ";
  }
  std::cout << "]\n";
}

void ExampleRS15_11() {
  std::cout << "=== 1. Reed-Solomon RS(15, 11) over GF(2^4) ===\n";

  using Field16 = GF2XLOG<uint8_t>;
  auto gf16 = std::make_shared<Field16>(4);
  auto rs = GeneralizedReedSolomonCode<Field16>::ReedSolomon(
      gf16, /*length=*/15, /*dimension=*/11);

  std::cout << "  Code: " << rs.ToString() << "\n";
  std::cout << "  Minimum distance d = " << rs.MinimumDistance()
            << ", t = " << ErrorCorrectionCapability(rs)
            << ", rate = " << std::fixed << std::setprecision(4) << CodeRate(rs)
            << "\n";
  std::cout << "  IsMDS=" << std::boolalpha << IsMDS(rs)
            << ", SatisfiesSingleton=" << SatisfiesSingletonBound(rs)
            << ", SatisfiesHamming=" << SatisfiesHammingBound(rs)
            << ", IsPerfect=" << IsPerfectCode(rs) << "\n";

  GFArray<Field16> message(
      gf16, xt::xarray<uint8_t>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11});
  const auto codeword = rs.Encode(message);
  PrintArray("  Message ", message);
  PrintArray("  Codeword", codeword);

  // Compare Berlekamp-Massey, Sugiyama, and PGZ on 2 symbol errors
  auto err_ch = channel::StaticErrorRateChannel<Field16>::Fixed(
      gf16, 15, /*num_errors=*/2, /*seed=*/101ULL);
  const auto rx_err = err_ch.Transmit(codeword);
  for (const char* dec_name : {"BerlekampMassey", "Sugiyama", "PGZ"}) {
    const auto res = rs.DecodeWithResult(rx_err, dec_name);
    std::cout << "  Decoder " << std::setw(16) << std::left << dec_name
              << " -> Ok=" << res.Ok() << ", errors=" << res.num_errors
              << ", exact_match=" << (res.message == message) << "\n";
  }

  // Errors-and-erasures decoding: 1 error + 2 erasures (2*1 + 2 = 4 <= n - k)
  channel::ErrorErasureChannel<Field16> ee_ch(
      gf16, 15, /*num_errors=*/1, /*num_erasures=*/2, /*seed=*/202ULL);
  const auto tx = ee_ch.TransmitWithErasures(codeword);
  BerlekampMasseyDecoder<Field16> bm_decoder(rs);
  const auto ee_res =
      bm_decoder.DecodeWithErasures(tx.received, tx.erasure_positions);
  std::cout << "  BerlekampMassey Errors+Erasures (e=1, s=2) -> Ok="
            << ee_res.Ok() << ", errors=" << ee_res.num_errors
            << ", erasures=" << ee_res.num_erasures
            << ", exact_match=" << (ee_res.message == message) << "\n\n";
}

void ExampleRS255_223() {
  std::cout << "=== 2. CCSDS Reed-Solomon RS(255, 223) over GF(2^8) ===\n";

  using Field256 = GF2XLOG<uint8_t>;
  auto gf256 = std::make_shared<Field256>(8);
  auto rs255 = GeneralizedReedSolomonCode<Field256>::ReedSolomon(
      gf256, /*length=*/255, /*dimension=*/223);

  std::cout << "  Code: " << rs255.ToString()
            << " (t = " << ErrorCorrectionCapability(rs255) << ")\n";

  std::mt19937_64 rng(42);
  const auto message = GFArray<Field256>::Random(gf256, {223}, rng);
  const auto codeword = rs255.Encode(message);

  // 10 symbol errors + 12 erasures => 2*10 + 12 = 32 == n - k (capacity limit)
  channel::ErrorErasureChannel<Field256> ee_ch(
      gf256, 255, /*num_errors=*/10, /*num_erasures=*/12, /*seed=*/777ULL);
  const auto tx = ee_ch.TransmitWithErasures(codeword);

  BerlekampMasseyDecoder<Field256> bm_decoder(rs255);
  const auto res =
      bm_decoder.DecodeWithErasures(tx.received, tx.erasure_positions);
  std::cout << "  Decoded 10 errors + 12 erasures -> Ok=" << std::boolalpha
            << res.Ok() << ", errors=" << res.num_errors
            << ", erasures=" << res.num_erasures
            << ", exact_match=" << (res.message == message) << "\n\n";
}

void ExampleBCHCode() {
  std::cout << "=== 3. Binary Primitive BCH(15, 7, t=2) Code ===\n";

  auto gf2 = std::make_shared<GF2>();
  auto gf16 = std::make_shared<GF2XLOG<uint8_t>>(4);
  auto bch = BCHCode<GF2, GF2XLOG<uint8_t>>::FromErrorCorrectionCapability(
      gf2, gf16, /*t=*/2);

  std::cout << "  Code: " << bch.ToString()
            << ", g(x) = " << bch.GeneratorPolynomial().ToString() << "\n";

  std::mt19937_64 rng(123);
  const auto message = GFArray<GF2>::Random(gf2, {bch.Dimension()}, rng);
  const auto codeword = bch.Encode(message);

  auto ch = channel::StaticErrorRateChannel<GF2>::Fixed(
      gf2, bch.Length(), /*num_errors=*/2, /*seed=*/321ULL);
  const auto rx = ch.Transmit(codeword);
  const auto res = bch.DecodeWithResult(rx);

  std::cout << "  Decoded 2 bit errors -> Ok=" << std::boolalpha << res.Ok()
            << ", errors=" << res.num_errors
            << ", exact_match=" << (res.message == message) << "\n\n";
}

void ExampleBinaryGolayCode() {
  std::cout << "=== 4. Perfect Binary Golay [23, 12, 7] Code ===\n";

  auto gf2 = std::make_shared<GF2>();
  auto golay = BinaryGolayCode(gf2);

  std::cout << "  Code: " << golay.ToString()
            << ", d = " << golay.MinimumDistance()
            << ", t = " << ErrorCorrectionCapability(golay)
            << ", IsPerfectCode=" << std::boolalpha << IsPerfectCode(golay)
            << "\n";

  std::mt19937_64 rng(456);
  const auto message = GFArray<GF2>::Random(gf2, {golay.Dimension()}, rng);
  const auto codeword = golay.Encode(message);

  auto ch = channel::StaticErrorRateChannel<GF2>::Fixed(
      gf2, golay.Length(), /*num_errors=*/3, /*seed=*/654ULL);
  const auto rx = ch.Transmit(codeword);
  const auto res = golay.DecodeWithResult(rx);

  std::cout << "  Decoded 3 bit errors -> Ok=" << res.Ok()
            << ", errors=" << res.num_errors
            << ", exact_match=" << (res.message == message) << "\n\n";
}

void ExampleSimulateCodePerformance() {
  std::cout << "=== 5. Monte Carlo Simulation (SimulateCodePerformance) ===\n";

  using Field16 = GF2XLOG<uint8_t>;
  auto gf16 = std::make_shared<Field16>(4);
  auto rs = GeneralizedReedSolomonCode<Field16>::ReedSolomon(gf16, 15, 11);

  channel::QarySymmetricChannel<Field16> qsc(gf16, 15, /*epsilon=*/0.06,
                                             /*seed=*/99ULL);
  const channel::SimulationConfig qsc_cfg{
      .max_frames = 400, .target_frame_errors = 50, .seed = 2026};
  const auto qsc_stats = channel::SimulateCodePerformance(rs, qsc, qsc_cfg);

  std::cout << "  RS(15,11) over QSC(eps=0.06): frames="
            << qsc_stats.frames_sent
            << ", frame_errors=" << qsc_stats.frame_errors << std::fixed
            << std::setprecision(4) << ", FER=" << qsc_stats.Fer()
            << ", SER=" << qsc_stats.Ser() << ", BER=" << qsc_stats.Ber()
            << "\n";

  channel::ErrorErasureChannel<Field16> ee_ch(
      gf16, 15, /*num_errors=*/1, /*num_erasures=*/2, /*seed=*/88ULL);
  const channel::SimulationConfig ee_cfg{
      .max_frames = 100, .target_frame_errors = 0, .seed = 2026};
  const auto ee_stats = channel::SimulateCodePerformance(rs, ee_ch, ee_cfg);

  std::cout << "  RS(15,11) over ErrorErasure(e=1,s=2): frames="
            << ee_stats.frames_sent
            << ", frame_errors=" << ee_stats.frame_errors
            << ", FER=" << ee_stats.Fer() << "\n\n";
}

void ExampleGuruswamiSudanListDecoding() {
  std::cout << "=== 6. Guruswami-Sudan List Decoding RS(15, 3) (e=8 > t=6) "
               "===\n";

  using Field16 = GF2XLOG<uint8_t>;
  auto gf16 = std::make_shared<Field16>(4);
  auto rs = GeneralizedReedSolomonCode<Field16>::ReedSolomon(
      gf16, /*length=*/15, /*dimension=*/3);

  // d = 13 -> classical unique radius t = 6, Johnson radius tau_J = 9.
  auto gs = GuruswamiSudanDecoder<Field16>::FromTargetRadius(rs, /*target=*/8);
  std::cout << "  Code: " << rs.ToString()
            << ", unique t = " << ErrorCorrectionCapability(rs)
            << ", Johnson radius = " << gs.JohnsonRadius()
            << ", m = " << gs.Multiplicity() << ", L = " << gs.ListSize()
            << "\n";

  std::mt19937_64 rng(2026);
  const auto message = GFArray<Field16>::Random(gf16, {3}, rng);
  const auto codeword = rs.Encode(message);

  // Corrupt 8 symbols (strictly beyond t = 6 where BM cannot decode).
  auto ch = channel::StaticErrorRateChannel<Field16>::Fixed(
      gf16, 15, /*num_errors=*/8, /*seed=*/77ULL);
  const auto received = ch.Transmit(codeword);

  const auto cw_list = gs.DecodeToList(received);
  const bool found_original =
      std::find(cw_list.begin(), cw_list.end(), codeword) != cw_list.end();
  const auto best = gs.DecodeWithResult(received);

  std::cout << "  List candidates within radius 8: " << cw_list.size()
            << ", contains transmitted codeword: " << std::boolalpha
            << found_original << "\n";
  std::cout << "  Closest codeword DecodeWithResult -> Ok=" << best.Ok()
            << ", distance=" << best.num_errors
            << ", exact_match=" << (best.message == message) << "\n\n";
}

void ExampleLdpcMinSumAwgn() {
  std::cout << "=== 7. Quasi-Cyclic LDPC Normalized Min-Sum over AWGN ===\n";

  auto gf2 = std::make_shared<GF2>();
  // Girth >= 6 regular (3, 6) QC-LDPC base matrix lifted by prime Z = 17.
  const std::vector<std::vector<int>> qc_base = {
      {0, 0, 0, 0, 0, 0},
      {0, 1, 2, 3, 4, 5},
      {0, 2, 4, 6, 8, 10},
  };
  auto ldpc = LDPCCode::FromQuasiCyclic(gf2, qc_base, /*lifting_factor_z=*/17);
  std::cout << "  Code: " << ldpc.ToString() << ", rate = " << std::fixed
            << std::setprecision(4) << ldpc.Rate() << "\n";

  BeliefPropagationDecoder nms_decoder(ldpc, BpAlgorithm::kNormalizedMinSum,
                                       /*max_iterations=*/50, /*alpha=*/0.75);
  auto awgn = channel::AWGNChannel::FromEbN0(
      gf2, ldpc.Length(), /*eb_n0_db=*/5.0, ldpc.Rate(), /*seed=*/1001ULL);

  std::mt19937_64 rng(1002);
  const auto message = GFArray<GF2>::Random(gf2, {ldpc.Dimension()}, rng);
  const auto codeword = ldpc.Encode(message);
  const auto llrs = awgn.TransmitLlr(codeword);
  const auto res = nms_decoder.DecodeLlr(llrs);

  std::cout << "  Decoder: " << nms_decoder.Name()
            << " at Eb/N0=" << awgn.SnrDb() << " dB -> Ok=" << std::boolalpha
            << res.Ok() << ", flipped_bits=" << res.num_errors
            << ", parity_ok=" << ldpc.CheckParity(res.codeword)
            << ", exact_match=" << (res.message == message) << "\n\n";
}

void ExamplePolarSclAwgn() {
  std::cout << "=== 8. Polar(64, 32) CRC-Aided SCL(L=8) over AWGN ===\n";

  auto gf2 = std::make_shared<GF2>();
  auto polar = PolarCode::FromDesignSnr(gf2, /*block_length=*/64,
                                        /*message_dimension=*/32,
                                        /*design_snr_db=*/2.5,
                                        /*systematic=*/false, /*crc_bits=*/8);
  std::cout << "  Code: " << polar.ToString() << "\n";

  SuccessiveCancellationDecoder sc_decoder(polar);
  SuccessiveCancellationListDecoder scl_decoder(polar, /*list_size=*/8);
  auto awgn = channel::AWGNChannel::FromEbN0(
      gf2, polar.Length(), /*eb_n0_db=*/4.5, polar.Rate(), /*seed=*/2026ULL);

  std::mt19937_64 rng(2027);
  const auto message = GFArray<GF2>::Random(gf2, {polar.Dimension()}, rng);
  const auto codeword = polar.Encode(message);
  const auto llrs = awgn.TransmitLlr(codeword);

  const auto sc_res = sc_decoder.DecodeLlr(llrs);
  const auto scl_res = scl_decoder.DecodeLlr(llrs);

  std::cout << "  SC         at Eb/N0=" << awgn.SnrDb()
            << " dB -> Ok=" << std::boolalpha << sc_res.Ok()
            << ", exact_match=" << (sc_res.Ok() && sc_res.message == message)
            << "\n";
  std::cout << "  CA-SCL(L=8) at Eb/N0=" << awgn.SnrDb()
            << " dB -> Ok=" << scl_res.Ok()
            << ", exact_match=" << (scl_res.Ok() && scl_res.message == message)
            << "\n\n";
}

}  // namespace

int main() {
  std::cout << "XGalois Coding Theory Examples\n";
  std::cout << "==============================\n\n";

  ExampleRS15_11();
  ExampleRS255_223();
  ExampleBCHCode();
  ExampleBinaryGolayCode();
  ExampleSimulateCodePerformance();
  ExampleGuruswamiSudanListDecoding();
  ExampleLdpcMinSumAwgn();
  ExamplePolarSclAwgn();

  std::cout << "All coding theory examples completed successfully!\n";
  return 0;
}
