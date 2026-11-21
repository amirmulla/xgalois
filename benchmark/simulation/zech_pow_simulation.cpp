#include <benchmark/benchmark.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "xgalois/field/gf_binary.hpp"

using xg::GFBEZechLogTables;

namespace {

std::vector<uint32_t> GenerateRandomZechElements(
    const GFBEZechLogTables<uint32_t>& field, size_t count,
    std::mt19937_64& gen) {
  std::vector<uint32_t> elements;
  elements.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    elements.push_back(field.Random(gen));
  }
  return elements;
}

std::vector<uint64_t> GenerateRandomExponents(size_t count, uint64_t min_exp,
                                              uint64_t max_exp,
                                              std::mt19937_64& gen) {
  std::uniform_int_distribution<uint64_t> exp_dis(min_exp, max_exp);
  std::vector<uint64_t> exps;
  exps.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    exps.push_back(exp_dis(gen));
  }
  return exps;
}

}  // namespace

int main() {
  constexpr uint8_t m = 20;
  constexpr size_t kNumElements = size_t{1} << 20;
  constexpr size_t kSmallBufSize = 1024;
  constexpr size_t kSmallBufMask = kSmallBufSize - 1;

  std::cout << "Creating GF(2^" << static_cast<int>(m)
            << ") using Zech Log Tables...\n";

  auto field = std::make_shared<GFBEZechLogTables<uint32_t>>(m);

  std::cout << "Field order: " << field->Order() << '\n';
  std::cout << "Field characteristic: " << field->Characteristic() << '\n';
  std::cout << "Field memory: " << field->MemoryBytes() << " bytes\n";

  // Sanity check Wave 1 SetElementValue / GetElementValue / Pow round-trip.
  const uint32_t sample = field->SetElementValue(2U);
  if (field->GetElementValue(field->Pow(sample, 3)) != 8U) {
    std::cerr << "Sanity check failed for Zech Pow!\n";
    return 1;
  }

  std::mt19937_64 gen(42);
  std::cout << "Generating " << kNumElements
            << " random field elements and exponents...\n";
  auto log_elements = GenerateRandomZechElements(*field, kNumElements, gen);
  auto exponents = GenerateRandomExponents(kNumElements, 0, 1000, gen);

  std::cout << "Measuring exponentiation of " << kNumElements
            << " random pairs...\n";

  auto start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < kNumElements; ++i) {
    uint32_t result = field->Pow(log_elements[i], exponents[i]);
    benchmark::DoNotOptimize(result);
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      end_time - start_time);

  double total_time_ms = static_cast<double>(duration.count()) / 1e6;
  double avg_time_ns =
      static_cast<double>(duration.count()) / static_cast<double>(kNumElements);
  double operations_per_second = 1e9 / avg_time_ns;

  std::cout << "\n=== Power Performance Results ===\n";
  std::cout << "Total operations: " << kNumElements << '\n';
  std::cout << "Total time: " << total_time_ms << " ms\n";
  std::cout << "Average time per exponentiation: " << avg_time_ns << " ns\n";
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  std::cout << "\n=== Cache-Resident Batch Power Test ===\n";
  auto batch_bases = GenerateRandomZechElements(*field, kSmallBufSize, gen);
  auto batch_exps = GenerateRandomExponents(kSmallBufSize, 0, 1000, gen);
  constexpr uint64_t kBatchOps = 1000000;

  start_time = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < kBatchOps; ++i) {
    const size_t idx = static_cast<size_t>(i) & kSmallBufMask;
    uint32_t result = field->Pow(batch_bases[idx], batch_exps[idx]);
    benchmark::DoNotOptimize(result);
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                  start_time);

  total_time_ms = static_cast<double>(duration.count()) / 1e6;
  avg_time_ns =
      static_cast<double>(duration.count()) / static_cast<double>(kBatchOps);
  operations_per_second = 1e9 / avg_time_ns;

  std::cout << "Batch operations: " << kBatchOps << '\n';
  std::cout << "Total time: " << total_time_ms << " ms\n";
  std::cout << "Average time per exponentiation: " << avg_time_ns << " ns\n";
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  std::cout << "\n=== Performance Comparison Across Field Sizes ===\n";
  for (uint8_t test_m : {4, 6, 8, 10, 12}) {
    auto test_field = std::make_shared<GFBEZechLogTables<uint32_t>>(test_m);
    auto test_bases =
        GenerateRandomZechElements(*test_field, kSmallBufSize, gen);
    auto test_exps = GenerateRandomExponents(kSmallBufSize, 0, 1000, gen);

    constexpr uint64_t kTestOps = 100000;
    start_time = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < kTestOps; ++i) {
      const size_t idx = static_cast<size_t>(i) & kSmallBufMask;
      uint32_t result = test_field->Pow(test_bases[idx], test_exps[idx]);
      benchmark::DoNotOptimize(result);
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                    start_time);

    avg_time_ns =
        static_cast<double>(duration.count()) / static_cast<double>(kTestOps);
    operations_per_second = 1e9 / avg_time_ns;

    std::cout << "GF(2^" << static_cast<int>(test_m)
              << ") [Order: " << test_field->Order()
              << ", Bytes: " << test_field->MemoryBytes()
              << "]: " << avg_time_ns << " ns/op, "
              << static_cast<uint64_t>(operations_per_second) << " ops/sec\n";
  }

  std::cout << "\n=== Performance with Different Exponent Ranges ===\n";
  const std::vector<std::pair<uint64_t, uint64_t>> exp_ranges = {
      {0, 10}, {0, 100}, {0, 1000}, {0, 10000}};

  for (const auto& range : exp_ranges) {
    // Pre-generate bases and exponents OUTSIDE the timed loop.
    auto range_bases = GenerateRandomZechElements(*field, kSmallBufSize, gen);
    auto test_exps =
        GenerateRandomExponents(kSmallBufSize, range.first, range.second, gen);

    constexpr uint64_t kTestOps = 100000;
    start_time = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < kTestOps; ++i) {
      const size_t idx = static_cast<size_t>(i) & kSmallBufMask;
      uint32_t result = field->Pow(range_bases[idx], test_exps[idx]);
      benchmark::DoNotOptimize(result);
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                    start_time);

    avg_time_ns =
        static_cast<double>(duration.count()) / static_cast<double>(kTestOps);
    operations_per_second = 1e9 / avg_time_ns;

    std::cout << "Exponent range [" << range.first << ", " << range.second
              << "]: " << avg_time_ns << " ns/op, "
              << static_cast<uint64_t>(operations_per_second) << " ops/sec\n";
  }

  std::cout << "\nPower simulation completed successfully!\n";
  return 0;
}
