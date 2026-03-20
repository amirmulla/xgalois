#include <benchmark/benchmark.h>

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "xgalois/field/gf_binary.hpp"

using namespace xg;

int main() {
  // Use GF(2^8) as an example field
  constexpr uint8_t m = 20;

  std::cout << "Creating GF(2^" << static_cast<int>(m)
            << ") using Zech Log Tables..." << std::endl;

  // Create the field using Zech Log Tables implementation
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::cout << "Field order: " << field->Order() << std::endl;
  std::cout << "Field characteristic: " << field->Characteristic() << std::endl;

  // Generate 1000 random field elements (in log representation)
  std::mt19937 gen(42);  // Fixed seed for reproducibility
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  uint64_t num_elements = 1e6;
  std::vector<uint32_t> log_elements;
  log_elements.reserve(num_elements);

  std::cout << "Generating " << num_elements << " random field elements..."
            << std::endl;
  for (size_t i = 0; i < num_elements; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  // Measure multiplication performance
  std::cout << "Measuring multiplication of " << num_elements
            << " random pairs..." << std::endl;

  auto start_time = std::chrono::high_resolution_clock::now();

  // Perform multiplications between adjacent pairs
  for (size_t i = 0; i < log_elements.size(); ++i) {
    uint32_t result = field->Mul(log_elements[i], log_elements[i + 1]);
    benchmark::DoNotOptimize(result);  // Prevent compiler optimization
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      end_time - start_time);

  // Calculate performance metrics
  double total_time_ms = duration.count() / 1e6;
  double avg_time_ns =
      static_cast<double>(duration.count()) / (log_elements.size());
  double operations_per_second = 1e9 / avg_time_ns;

  // Display results
  std::cout << "\n=== Multiplication Performance Results ===" << std::endl;
  std::cout << "Total operations: " << (log_elements.size()) << std::endl;
  std::cout << "Total time: " << total_time_ms << " ms" << std::endl;
  std::cout << "Average time per multiplication: " << avg_time_ns << " ns"
            << std::endl;
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << std::endl;

  // Additional test: measure a batch of multiplications with the same operands
  std::cout << "\n=== Batch Multiplication Test ===" << std::endl;
  // Pick two random elements for repeated multiplication
  std::uniform_int_distribution<size_t> index_dis(0, log_elements.size() - 1);
  uint32_t a = log_elements[index_dis(gen)];
  uint32_t b = log_elements[index_dis(gen)];

  constexpr uint64_t batch_size = 1000000;  // 1 million operations

  start_time = std::chrono::high_resolution_clock::now();

  for (uint64_t i = 0; i < batch_size; ++i) {
    uint32_t result = field->Mul(a, b);
    benchmark::DoNotOptimize(result);
  }

  end_time = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                  start_time);

  total_time_ms = duration.count() / 1e6;
  avg_time_ns = static_cast<double>(duration.count()) / batch_size;
  operations_per_second = 1e9 / avg_time_ns;

  std::cout << "Batch operations: " << batch_size << std::endl;
  std::cout << "Total time: " << total_time_ms << " ms" << std::endl;
  std::cout << "Average time per multiplication: " << avg_time_ns << " ns"
            << std::endl;
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << std::endl;

  // Test with different field sizes
  std::cout << "\n=== Performance Comparison Across Field Sizes ==="
            << std::endl;

  for (uint8_t test_m : {4, 6, 8, 10, 12}) {
    auto test_field = std::make_shared<GFBEZechLogTables>(test_m);

    // Generate some test elements
    std::uniform_int_distribution<uint64_t> test_dis(0,
                                                     test_field->Order() - 2);
    uint32_t test_a = static_cast<uint32_t>(test_dis(gen));
    uint32_t test_b = static_cast<uint32_t>(test_dis(gen));

    constexpr uint64_t test_operations = 100000;

    start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < test_operations; ++i) {
      uint32_t result = test_field->Mul(test_a, test_b);
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

  std::cout << "\nSimulation completed successfully!" << std::endl;
  return 0;
}
