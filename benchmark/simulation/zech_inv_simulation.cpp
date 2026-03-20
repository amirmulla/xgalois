#include <benchmark/benchmark.h>

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "xgalois/field/gf_binary.hpp"

using namespace xg;

int main() {
  // Use GF(2^20) as an example field
  constexpr uint8_t m = 20;

  std::cout << "Creating GF(2^" << static_cast<int>(m)
            << ") using Zech Log Tables..." << '\n';

  // Create the field using Zech Log Tables implementation
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::cout << "Field order: " << field->Order() << std::endl;
  std::cout << "Field characteristic: " << field->Characteristic() << std::endl;

  // Generate 1000000 random field elements (in log representation)
  std::mt19937 gen(42);  // Fixed seed for reproducibility
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  uint64_t num_elements = 1e6;
  std::vector<uint32_t> log_elements;
  log_elements.reserve(num_elements);

  std::cout << "Generating " << num_elements << " random field elements..."
            << '\n';
  for (size_t i = 0; i < num_elements; ++i) {
    // Generate a random element in the field, avoiding zero (ZECH_INFINITY)
    uint32_t elm = static_cast<uint32_t>(dis(gen));
    while (elm == static_cast<uint32_t>(-1)) {  // -1 represents ZECH_INFINITY
      elm = static_cast<uint32_t>(dis(gen));
    }
    log_elements.push_back(elm);
  }

  // Measure inversion performance
  std::cout << "Measuring inversion of " << num_elements << " random pairs..."
            << '\n';

  auto start_time = std::chrono::high_resolution_clock::now();

  // Perform inversions between adjacent pairs (avoiding inversion by zero)
  for (size_t i = 0; i < log_elements.size() - 1; ++i) {
    uint32_t result = field->Inv(log_elements[i]);
    benchmark::DoNotOptimize(result);  // Prevent compiler optimization
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      end_time - start_time);

  // Calculate performance metrics
  double total_time_ms = duration.count() / 1e6;
  double avg_time_ns =
      static_cast<double>(duration.count()) / (log_elements.size() - 1);
  double operations_per_second = 1e9 / avg_time_ns;

  // Display results
  std::cout << "\n=== inversion Performance Results ===" << '\n';
  std::cout << "Total operations: " << (log_elements.size() - 1) << '\n';
  std::cout << "Total time: " << total_time_ms << " ms" << '\n';
  std::cout << "Average time per inversion: " << avg_time_ns << " ns" << '\n';
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  // Additional test: measure a batch of inversions with the same operands
  std::cout << "\n=== Batch inversion Test ===" << '\n';
  // Pick two random elements for repeated inversion
  std::uniform_int_distribution<size_t> index_dis(0, log_elements.size() - 1);

  // Ensure b is not zero (ZECH_INFINITY)
  auto elm = log_elements[index_dis(gen)];

  constexpr uint64_t batch_size = 1000000;  // 1 million operations

  start_time = std::chrono::high_resolution_clock::now();

  for (uint64_t i = 0; i < batch_size; ++i) {
    uint32_t result = field->Inv(elm);
    benchmark::DoNotOptimize(result);
  }

  end_time = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                  start_time);

  total_time_ms = duration.count() / 1e6;
  avg_time_ns = static_cast<double>(duration.count()) / batch_size;
  operations_per_second = 1e9 / avg_time_ns;

  std::cout << "Batch operations: " << batch_size << '\n';
  std::cout << "Total time: " << total_time_ms << " ms" << '\n';
  std::cout << "Average time per inversion: " << avg_time_ns << " ns" << '\n';
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  // Test with different field sizes
  std::cout << "\n=== Performance Comparison Across Field Sizes ===" << '\n';

  for (uint8_t test_m : {4, 6, 8, 10, 12}) {
    auto test_field = std::make_shared<GFBEZechLogTables>(test_m);

    // Generate some test elements
    std::uniform_int_distribution<uint64_t> test_dis(0,
                                                     test_field->Order() - 2);
    uint32_t test_elm = static_cast<uint32_t>(test_dis(gen));

    constexpr uint64_t test_operations = 100000;

    start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < test_operations; ++i) {
      uint32_t result = test_field->Inv(test_elm);
      benchmark::DoNotOptimize(result);
    }

    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                    start_time);

    avg_time_ns = static_cast<double>(duration.count()) / test_operations;
    operations_per_second = 1e9 / avg_time_ns;

    std::cout << "GF(2^" << static_cast<int>(test_m)
              << ") [Order: " << test_field->Order() << "]: " << avg_time_ns
              << " ns/op, " << static_cast<uint64_t>(operations_per_second)
              << " ops/sec" << std::endl;
  }

  std::cout << "\ninversion simulation completed successfully!" << '\n';
  return 0;
}
