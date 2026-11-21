#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
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
  channel::StaticErrorRateChannel<Field16> err_ch(
      gf16, 15, /*num_errors=*/2, std::optional<uint64_t>{101ULL});
  const auto rx_err = err_ch.transmit(codeword);
  for (const char* dec_name : {"BerlekampMassey", "Sugiyama", "PGZ"}) {
    const auto res = rs.DecodeWithResult(rx_err, dec_name);
    std::cout << "  Decoder " << std::setw(16) << std::left << dec_name
              << " -> Ok=" << res.Ok() << ", errors=" << res.num_errors
              << ", exact_match=" << (res.message == message) << "\n";
  }

  // Errors-and-erasures decoding: 1 error + 2 erasures (2*1 + 2 = 4 <= n - k)
  channel::ErrorErasureChannel<Field16> ee_ch(gf16, 15, /*num_errors=*/1,
                                              /*num_erasures=*/2, /*seed=*/202);
  const auto tx = ee_ch.transmit_with_erasures(codeword);
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
      gf256, 255, /*num_errors=*/10, /*num_erasures=*/12, /*seed=*/777);
  const auto tx = ee_ch.transmit_with_erasures(codeword);

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

  channel::StaticErrorRateChannel<GF2> ch(gf2, bch.Length(), /*num_errors=*/2,
                                          std::optional<uint64_t>{321ULL});
  const auto rx = ch.transmit(codeword);
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

  channel::StaticErrorRateChannel<GF2> ch(gf2, golay.Length(), /*num_errors=*/3,
                                          std::optional<uint64_t>{654ULL});
  const auto rx = ch.transmit(codeword);
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
                                             /*seed=*/99);
  const channel::SimulationConfig qsc_cfg{
      .max_frames = 400, .target_frame_errors = 50, .seed = 2026};
  const auto qsc_stats = channel::SimulateCodePerformance(rs, qsc, qsc_cfg);

  std::cout << "  RS(15,11) over QSC(eps=0.06): frames="
            << qsc_stats.frames_sent
            << ", frame_errors=" << qsc_stats.frame_errors << std::fixed
            << std::setprecision(4) << ", FER=" << qsc_stats.FER()
            << ", SER=" << qsc_stats.SER() << ", BER=" << qsc_stats.BER()
            << "\n";

  channel::ErrorErasureChannel<Field16> ee_ch(gf16, 15, /*num_errors=*/1,
                                              /*num_erasures=*/2, /*seed=*/88);
  const channel::SimulationConfig ee_cfg{
      .max_frames = 100, .target_frame_errors = 0, .seed = 2026};
  const auto ee_stats = channel::SimulateCodePerformance(rs, ee_ch, ee_cfg);

  std::cout << "  RS(15,11) over ErrorErasure(e=1,s=2): frames="
            << ee_stats.frames_sent
            << ", frame_errors=" << ee_stats.frame_errors
            << ", FER=" << ee_stats.FER() << "\n\n";
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

  std::cout << "All coding theory examples completed successfully!\n";
  return 0;
}
