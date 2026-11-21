#include <benchmark/benchmark.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
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

}  // namespace

int main() {
  constexpr uint8_t m = 20;
  constexpr size_t kNumElements = size_t{1} << 20;
  constexpr size_t kNumMask = kNumElements - 1;
  constexpr size_t kSmallBufSize = 1024;
  constexpr size_t kSmallBufMask = kSmallBufSize - 1;

  std::cout << "Creating GF(2^" << static_cast<int>(m)
            << ") using Zech Log Tables...\n";

  auto field = std::make_shared<GFBEZechLogTables<uint32_t>>(m);

  std::cout << "Field order: " << field->Order() << '\n';
  std::cout << "Field characteristic: " << field->Characteristic() << '\n';
  std::cout << "Field memory: " << field->MemoryBytes() << " bytes\n";

  // Sanity check Wave 1 SetElementValue / GetElementValue / Add round-trip.
  const uint32_t sample_a = field->SetElementValue(3U);
  const uint32_t sample_b = field->SetElementValue(5U);
  const uint32_t sample_sum =
      field->GetElementValue(field->Add(sample_a, sample_b));
  if (sample_sum != (3U ^ 5U)) {
    std::cerr << "Sanity check failed for Zech Add!\n";
    return 1;
  }

  std::mt19937_64 gen(42);
  std::cout << "Generating " << kNumElements << " random field elements...\n";
  auto log_elements = GenerateRandomZechElements(*field, kNumElements, gen);

  std::cout << "Measuring addition of " << kNumElements << " random pairs...\n";

  auto start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < kNumElements; ++i) {
    uint32_t result =
        field->Add(log_elements[i], log_elements[(i + 1) & kNumMask]);
    benchmark::DoNotOptimize(result);
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      end_time - start_time);

  double total_time_ms = static_cast<double>(duration.count()) / 1e6;
  double avg_time_ns =
      static_cast<double>(duration.count()) / static_cast<double>(kNumElements);
  double operations_per_second = 1e9 / avg_time_ns;

  std::cout << "\n=== Addition Performance Results ===\n";
  std::cout << "Total operations: " << kNumElements << '\n';
  std::cout << "Total time: " << total_time_ms << " ms\n";
  std::cout << "Average time per addition: " << avg_time_ns << " ns\n";
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  std::cout << "\n=== Cache-Resident Batch Addition Test ===\n";
  auto batch_elements = GenerateRandomZechElements(*field, kSmallBufSize, gen);
  constexpr uint64_t kBatchOps = 1000000;

  start_time = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < kBatchOps; ++i) {
    const size_t idx = static_cast<size_t>(i) & kSmallBufMask;
    uint32_t result = field->Add(batch_elements[idx],
                                 batch_elements[(idx + 1) & kSmallBufMask]);
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
  std::cout << "Average time per addition: " << avg_time_ns << " ns\n";
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  std::cout << "\n=== Performance Comparison Across Field Sizes ===\n";
  for (uint8_t test_m : {4, 6, 8, 10, 12}) {
    auto test_field = std::make_shared<GFBEZechLogTables<uint32_t>>(test_m);
    auto test_elems =
        GenerateRandomZechElements(*test_field, kSmallBufSize, gen);

    constexpr uint64_t kTestOps = 100000;
    start_time = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < kTestOps; ++i) {
      const size_t idx = static_cast<size_t>(i) & kSmallBufMask;
      uint32_t result = test_field->Add(test_elems[idx],
                                        test_elems[(idx + 1) & kSmallBufMask]);
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

  std::cout << "\nAddition simulation completed successfully!\n";
  return 0;
}
