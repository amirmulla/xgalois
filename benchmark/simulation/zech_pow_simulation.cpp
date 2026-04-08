#include <benchmark/benchmark.h>

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "xgalois/field/gf_binary.hpp"

using xg::GFBEZechLogTables;

int main() {
  constexpr uint8_t m = 20;

  std::cout << "Creating GF(2^" << static_cast<int>(m)
            << ") using Zech Log Tables..." << '\n';

  auto field = std::make_shared<GFBEZechLogTables<uint32_t>>(m);

  std::cout << "Field order: " << field->Order() << '\n';
  std::cout << "Field characteristic: " << field->Characteristic() << '\n';

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);
  std::uniform_int_distribution<uint32_t> exp_dis(0, 1000);

  uint64_t num_elements = 1e6;
  std::vector<uint32_t> log_elements;
  std::vector<uint32_t> exponents;
  log_elements.reserve(num_elements);
  exponents.reserve(num_elements);

  std::cout << "Generating " << num_elements
            << " random field elements and exponents..." << '\n';
  for (size_t i = 0; i < num_elements; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
    exponents.push_back(exp_dis(gen));
  }

  std::cout << "Measuring exponentiation of " << num_elements
            << " random pairs..." << '\n';

  auto start_time = std::chrono::high_resolution_clock::now();

  for (size_t i = 0; i < log_elements.size(); ++i) {
    uint32_t result = field->Pow(log_elements[i], exponents[i]);
    benchmark::DoNotOptimize(result);
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
      end_time - start_time);

  double total_time_ms = duration.count() / 1e6;
  double avg_time_ns =
      static_cast<double>(duration.count()) / log_elements.size();
  double operations_per_second = 1e9 / avg_time_ns;

  std::cout << "\n=== Power Performance Results ===" << '\n';
  std::cout << "Total operations: " << log_elements.size() << '\n';
  std::cout << "Total time: " << total_time_ms << " ms" << '\n';
  std::cout << "Average time per exponentiation: " << avg_time_ns << " ns"
            << '\n';
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  std::cout << "\n=== Batch Power Test ===" << '\n';

  std::uniform_int_distribution<size_t> index_dis(0, log_elements.size() - 1);
  uint32_t base = log_elements[index_dis(gen)];
  uint32_t exp = exponents[index_dis(gen)];

  constexpr uint64_t batch_size = 1000000;

  start_time = std::chrono::high_resolution_clock::now();

  for (uint64_t i = 0; i < batch_size; ++i) {
    uint32_t result = field->Pow(base, exp);
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
  std::cout << "Average time per exponentiation: " << avg_time_ns << " ns"
            << '\n';
  std::cout << "Operations per second: "
            << static_cast<uint64_t>(operations_per_second) << '\n';

  std::cout << "\n=== Performance Comparison Across Field Sizes ===" << '\n';

  for (uint8_t test_m : {4, 6, 8, 10, 12}) {
    auto test_field = std::make_shared<GFBEZechLogTables<uint32_t>>(test_m);

    std::uniform_int_distribution<uint64_t> test_dis(0,
                                                     test_field->Order() - 2);
    uint32_t test_base = static_cast<uint32_t>(test_dis(gen));
    uint32_t test_exp = exp_dis(gen);

    constexpr uint64_t test_operations = 100000;

    start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < test_operations; ++i) {
      uint32_t result = test_field->Pow(test_base, test_exp);
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
              << " ops/sec" << '\n';
  }

  std::cout << "\n=== Performance with Different Exponent Ranges ===" << '\n';

  std::vector<std::pair<uint32_t, uint32_t>> exp_ranges = {
      {0, 10}, {0, 100}, {0, 1000}, {0, 10000}};

  for (const auto& range : exp_ranges) {
    std::uniform_int_distribution<uint32_t> range_exp_dis(range.first,
                                                          range.second);
    uint32_t test_base = log_elements[index_dis(gen)];

    constexpr uint64_t test_operations = 100000;

    start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < test_operations; ++i) {
      uint32_t test_exp = range_exp_dis(gen);
      uint32_t result = field->Pow(test_base, test_exp);
      benchmark::DoNotOptimize(result);
    }

    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                    start_time);

    avg_time_ns = static_cast<double>(duration.count()) / test_operations;
    operations_per_second = 1e9 / avg_time_ns;

    std::cout << "Exponent range [" << range.first << ", " << range.second
              << "]: " << avg_time_ns << " ns/op, "
              << static_cast<uint64_t>(operations_per_second) << " ops/sec"
              << '\n';
  }

  std::cout << "\nPower simulation completed successfully!" << '\n';
  return 0;
}
