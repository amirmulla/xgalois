#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "xgalois/channel/channel.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

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

template <typename Field>
size_t CountDifferences(const GFArray<Field>& a, const GFArray<Field>& b) {
  size_t diffs = 0;
  for (size_t i = 0; i < a.Size(); ++i) {
    if (a.Raw(i) != b.Raw(i)) {
      ++diffs;
    }
  }
  return diffs;
}

}  // namespace

int main() {
  std::cout << "XGalois Channel Examples\n";
  std::cout << "========================\n";

  constexpr size_t kMsgLen = 8;

  std::cout << "\n=== GF(2^4) Channels (GFArray + Seeded RNGs) ===\n";
  using Field16 = GF2XLOG<uint8_t>;
  auto gf16 = std::make_shared<Field16>(4);
  GFArray<Field16> msg16(gf16, xt::xarray<uint8_t>{1, 2, 3, 4, 5, 6, 7, 8});
  PrintArray("Original GF(2^4) message", msg16);

  std::cout << "\n--- StaticErrorRateChannel::Fixed (2 errors, seed=42) ---\n";
  auto static_ch = channel::StaticErrorRateChannel<Field16>::Fixed(
      gf16, kMsgLen, /*num_errors=*/2, /*seed=*/42ULL);
  std::cout << "Channel: " << static_ch << "\n";
  const auto rx_static = static_ch.Transmit(msg16);
  PrintArray("Received", rx_static);
  std::cout << "Symbol errors introduced: "
            << CountDifferences(msg16, rx_static) << "\n";

  std::cout << "\n--- StaticErrorRateChannel (1..3 errors, seed=100) ---\n";
  channel::StaticErrorRateChannel<Field16> var_ch(gf16, kMsgLen, 1, 3,
                                                  /*seed=*/100ULL);
  for (int trial = 1; trial <= 3; ++trial) {
    const auto rx = var_ch.Transmit(msg16);
    std::cout << "  Trial " << trial
              << " -> errors: " << CountDifferences(msg16, rx) << "\n";
  }

  std::cout
      << "\n--- ErrorErasureChannel (2 errors + 1 erasure, seed=77) ---\n";
  channel::ErrorErasureChannel<Field16> ee_ch(gf16, kMsgLen, 2, 1,
                                              /*seed=*/77ULL);
  std::cout << "Channel: " << ee_ch << "\n";
  const auto ee_tx = ee_ch.TransmitWithErasures(msg16);
  PrintArray("Received", ee_tx.received);
  std::cout << "Erasure positions: [";
  for (size_t i = 0; i < ee_tx.erasure_positions.size(); ++i) {
    std::cout << ee_tx.erasure_positions[i];
    if (i + 1 < ee_tx.erasure_positions.size()) std::cout << ", ";
  }
  std::cout << "]\n";
  std::cout << "Erasure mask:      [";
  for (size_t i = 0; i < ee_tx.erasure_mask.size(); ++i) {
    std::cout << (ee_tx.erasure_mask[i] ? "1" : "0");
    if (i + 1 < ee_tx.erasure_mask.size()) std::cout << ", ";
  }
  std::cout << "]\n";

  std::cout << "\n=== GF(31) QarySymmetricChannel (seed=2026) ===\n";
  using Field31 = GFP<uint32_t>;
  auto gf31 = std::make_shared<Field31>(31);
  GFArray<Field31> msg31(gf31, xt::xarray<uint32_t>{1, 2, 3, 4, 5, 6, 7, 8});
  PrintArray("Original GF(31) message", msg31);

  channel::QarySymmetricChannel<Field31> qsc(gf31, kMsgLen, /*epsilon=*/0.2,
                                             /*seed=*/2026ULL);
  std::cout << "Channel: " << qsc << "\n";
  std::cout << std::fixed << std::setprecision(4)
            << "P(== 2 errors) = " << qsc.ProbabilityOfExactlyTErrors(2)
            << ", P(<= 2 errors) = " << qsc.ProbabilityOfAtMostTErrors(2)
            << "\n";
  const auto rx_qsc = qsc.Transmit(msg31);
  PrintArray("Received", rx_qsc);
  std::cout << "Symbol errors introduced: " << CountDifferences(msg31, rx_qsc)
            << "\n";

  std::cout << "\n=== GF(2) BinarySymmetricChannel & BPSK AWGNChannel ===\n";
  auto gf2 = std::make_shared<GF2>();
  constexpr size_t kBitLen = 12;
  std::mt19937_64 rng(99);
  const auto bits = GFArray<GF2>::Random(gf2, {kBitLen}, rng);
  PrintArray("Original bits", bits);

  channel::BinarySymmetricChannel bsc(gf2, kBitLen, /*crossover_prob=*/0.1,
                                      /*seed=*/55ULL);
  const auto rx_bsc = bsc.Transmit(bits);
  std::cout << "BSC (p=0.1) capacity = " << bsc.Capacity()
            << " bits/use, bit flips = " << CountDifferences(bits, rx_bsc)
            << "\n";

  auto awgn = channel::AWGNChannel::FromSNR(gf2, kBitLen, /*snr_db=*/4.0,
                                            /*code_rate=*/1.0, /*seed=*/314ULL);
  std::cout << "Channel: " << awgn << "\n";
  std::cout << "  snr_db = " << awgn.SnrDb() << ", sigma = " << awgn.Sigma()
            << ", theoretical uncoded BER = " << awgn.TheoreticalBer() << "\n";
  const auto awgn_tx = awgn.TransmitDetailed(bits);
  PrintArray("AWGN hard decision", awgn_tx.hard_decision);
  std::cout << "  First 4 soft samples: [";
  for (size_t i = 0; i < 4; ++i) {
    std::cout << std::setprecision(3) << awgn_tx.soft_received[i];
    if (i < 3) std::cout << ", ";
  }
  std::cout << "]\n  First 4 LLRs:         [";
  for (size_t i = 0; i < 4; ++i) {
    std::cout << std::setprecision(3) << awgn_tx.llr[i];
    if (i < 3) std::cout << ", ";
  }
  std::cout << "]\n  Hard-decision bit errors: "
            << CountDifferences(bits, awgn_tx.hard_decision) << "\n";

  std::cout << "\nChannel examples completed successfully!\n";
  return 0;
}
