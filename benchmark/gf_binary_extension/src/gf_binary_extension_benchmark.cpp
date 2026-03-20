/**
 * @file gf_binary_extension_benchmark.cpp
 * @brief Comprehensive benchmark comparison between binary extension field
 * implementations Compares 4 different implementation classes across small,
 * medium, and large field sizes
 */

#include <benchmark/benchmark.h>
#include <sys/resource.h>
#include <unistd.h>

#include <chrono>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"

using namespace xg;

//------------------------------------------------------------------------------
// Memory Usage Utilities
//------------------------------------------------------------------------------

struct MemoryUsage {
  size_t peak_rss_kb;
  size_t current_rss_kb;

  MemoryUsage() : peak_rss_kb(0), current_rss_kb(0) {}
};

MemoryUsage GetMemoryUsage() {
  MemoryUsage usage;

  // Get peak memory usage
  struct rusage rusage_data;
  if (getrusage(RUSAGE_SELF, &rusage_data) == 0) {
    usage.peak_rss_kb = rusage_data.ru_maxrss / 1024;  // Convert to KB on macOS
  }

  // Get current memory usage using ps command (macOS compatible)
  pid_t pid = getpid();
  char command[256];
  snprintf(command, sizeof(command), "ps -o rss= -p %d", pid);

  FILE* pipe = popen(command, "r");
  if (pipe) {
    char buffer[128];
    if (fgets(buffer, sizeof(buffer), pipe)) {
      usage.current_rss_kb = atol(buffer);
    }
    pclose(pipe);
  }

  return usage;
}

//------------------------------------------------------------------------------
// Field Sizes for Testing
//------------------------------------------------------------------------------

// Small fields - Common small extension degrees
const std::vector<uint8_t> SMALL_FIELD_DEGREES = {2, 3, 4, 5, 6, 7, 8};

// Medium fields - Practical sizes for applications
const std::vector<uint8_t> MEDIUM_FIELD_DEGREES = {9,  10, 11, 12,
                                                   13, 14, 15, 16};

// Large fields - Note: Some implementations may be memory intensive
const std::vector<uint8_t> LARGE_FIELD_DEGREES = {17, 18, 19, 20};

//------------------------------------------------------------------------------
// Helper Functions
//------------------------------------------------------------------------------

template <typename FieldType>
std::vector<typename FieldType::element_type> GenerateRandomElements(
    std::shared_ptr<FieldType> field, size_t count, uint32_t seed = 42) {
  std::mt19937 gen(seed);
  std::vector<typename FieldType::element_type> elements;
  elements.reserve(count);

  for (size_t i = 0; i < count; ++i) {
    auto elem = field->Random();
    // Ensure we don't have zero for operations that don't allow it
    while (elem == 0) {
      elem = field->Random();
    }
    elements.push_back(elem);
  }

  return elements;
}

//------------------------------------------------------------------------------
// Small Field Benchmarks (GF(2^m) where m = 2-8)
//------------------------------------------------------------------------------

// GaloisFieldBinaryExtension - Base Implementation
static void BM_GF2X_Small_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Small_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Small_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Small_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Small_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// GFBELogTables - Log Table Implementation
static void BM_GF2XLOG_Small_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Small_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Small_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Small_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Small_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// GFBELogTablesOpt - Optimized Log Table Implementation
static void BM_GF2XLOGOPT_Small_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Small_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Small_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Small_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Small_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// GFBEZechLogTables - Zech Log Table Implementation
static void BM_GF2XZECH_Small_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  // For Zech logs, we work with logarithmic representations
  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 1000; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Small_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 1000; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Small_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 1000; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Small_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 1000; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(log_elements[idx % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Small_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 1000; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(log_elements[idx % log_elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

//------------------------------------------------------------------------------
// Medium Field Benchmarks (GF(2^m) where m = 9-16)
//------------------------------------------------------------------------------

// Base Implementation - Medium Fields
static void BM_GF2X_Medium_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Medium_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Medium_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Medium_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Medium_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 1024 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// Log Table Implementation - Medium Fields
static void BM_GF2XLOG_Medium_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Medium_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Medium_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Medium_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Medium_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 1024 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// Optimized Log Table Implementation - Medium Fields
static void BM_GF2XLOGOPT_Medium_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Medium_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Medium_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Medium_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Medium_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint32_t>>(m);

  auto elements = GenerateRandomElements(field, 500);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 1024 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// Zech Log Table Implementation - Medium Fields
static void BM_GF2XZECH_Medium_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 500; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Medium_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 500; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Medium_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 500; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Medium_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 500; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(log_elements[idx % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Medium_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 500; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(log_elements[idx % log_elements.size()],
                             static_cast<uint64_t>(idx % 1024 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

//------------------------------------------------------------------------------
// Large Field Benchmarks (GF(2^m)
//------------------------------------------------------------------------------

// Base Implementation - Large Fields
static void BM_GF2X_Large_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Large_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Large_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Large_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2X_Large_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldBinaryExtension<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 4096 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// Log Table Implementation - Large Fields
static void BM_GF2XLOG_Large_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Large_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Large_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Large_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOG_Large_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTables<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 4096 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// Optimized Log Table Implementation - Large Fields
static void BM_GF2XLOGOPT_Large_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Large_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Large_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Large_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XLOGOPT_Large_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBELogTablesOpt<uint64_t>>(m);

  auto elements = GenerateRandomElements(field, 200);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 4096 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

// Zech Log Table Implementation - Large Fields
static void BM_GF2XZECH_Large_Addition(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 200; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Large_Multiplication(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 200; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Large_Division(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 200; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(log_elements[idx % log_elements.size()],
                             log_elements[(idx + 1) % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Large_Inversion(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 200; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(log_elements[idx % log_elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

static void BM_GF2XZECH_Large_Exponentiation(benchmark::State& state) {
  uint8_t m = static_cast<uint8_t>(state.range(0));
  auto field = std::make_shared<GFBEZechLogTables>(m);

  std::mt19937 gen(42);
  std::uniform_int_distribution<uint64_t> dis(0, field->Order() - 2);

  std::vector<uint32_t> log_elements;
  for (size_t i = 0; i < 200; ++i) {
    log_elements.push_back(static_cast<uint32_t>(dis(gen)));
  }

  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(log_elements[idx % log_elements.size()],
                             static_cast<uint64_t>(idx % 4096 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  MemoryUsage mem_end = GetMemoryUsage();
  state.counters["MemoryPeak_KB"] = mem_end.peak_rss_kb;
  state.counters["FieldOrder"] = field->Order();
}

//------------------------------------------------------------------------------
// Benchmark Registration
//------------------------------------------------------------------------------

// Small Field Benchmarks (GF(2^m)
BENCHMARK(BM_GF2X_Small_Addition)->DenseRange(2, 8);
BENCHMARK(BM_GF2X_Small_Multiplication)->DenseRange(2, 8);
BENCHMARK(BM_GF2X_Small_Division)->DenseRange(2, 8);
BENCHMARK(BM_GF2X_Small_Inversion)->DenseRange(2, 8);
BENCHMARK(BM_GF2X_Small_Exponentiation)->DenseRange(2, 8);

BENCHMARK(BM_GF2XLOG_Small_Addition)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOG_Small_Multiplication)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOG_Small_Division)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOG_Small_Inversion)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOG_Small_Exponentiation)->DenseRange(2, 8);

BENCHMARK(BM_GF2XLOGOPT_Small_Addition)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOGOPT_Small_Multiplication)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOGOPT_Small_Division)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOGOPT_Small_Inversion)->DenseRange(2, 8);
BENCHMARK(BM_GF2XLOGOPT_Small_Exponentiation)->DenseRange(2, 8);

BENCHMARK(BM_GF2XZECH_Small_Addition)->DenseRange(2, 8);
BENCHMARK(BM_GF2XZECH_Small_Multiplication)->DenseRange(2, 8);
BENCHMARK(BM_GF2XZECH_Small_Division)->DenseRange(2, 8);
BENCHMARK(BM_GF2XZECH_Small_Inversion)->DenseRange(2, 8);
BENCHMARK(BM_GF2XZECH_Small_Exponentiation)->DenseRange(2, 8);

// Medium Field Benchmarks (GF(2^m) where m = 9-16)
BENCHMARK(BM_GF2X_Medium_Addition)->DenseRange(9, 16);
BENCHMARK(BM_GF2X_Medium_Multiplication)->DenseRange(9, 16);
BENCHMARK(BM_GF2X_Medium_Division)->DenseRange(9, 16);
BENCHMARK(BM_GF2X_Medium_Inversion)->DenseRange(9, 16);
BENCHMARK(BM_GF2X_Medium_Exponentiation)->DenseRange(9, 16);

BENCHMARK(BM_GF2XLOG_Medium_Addition)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOG_Medium_Multiplication)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOG_Medium_Division)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOG_Medium_Inversion)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOG_Medium_Exponentiation)->DenseRange(9, 16);

BENCHMARK(BM_GF2XLOGOPT_Medium_Addition)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOGOPT_Medium_Multiplication)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOGOPT_Medium_Division)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOGOPT_Medium_Inversion)->DenseRange(9, 16);
BENCHMARK(BM_GF2XLOGOPT_Medium_Exponentiation)->DenseRange(9, 16);

BENCHMARK(BM_GF2XZECH_Medium_Addition)->DenseRange(9, 16);
BENCHMARK(BM_GF2XZECH_Medium_Multiplication)->DenseRange(9, 16);
BENCHMARK(BM_GF2XZECH_Medium_Division)->DenseRange(9, 16);
BENCHMARK(BM_GF2XZECH_Medium_Inversion)->DenseRange(9, 16);
BENCHMARK(BM_GF2XZECH_Medium_Exponentiation)->DenseRange(9, 16);

// Large Field Benchmarks (GF(2^m) where m = 17-20)
BENCHMARK(BM_GF2X_Large_Addition)->DenseRange(17, 20);
BENCHMARK(BM_GF2X_Large_Multiplication)->DenseRange(17, 20);
BENCHMARK(BM_GF2X_Large_Division)->DenseRange(17, 20);
BENCHMARK(BM_GF2X_Large_Inversion)->DenseRange(17, 20);
BENCHMARK(BM_GF2X_Large_Exponentiation)->DenseRange(17, 20);

BENCHMARK(BM_GF2XLOG_Large_Addition)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOG_Large_Multiplication)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOG_Large_Division)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOG_Large_Inversion)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOG_Large_Exponentiation)->DenseRange(17, 20);

BENCHMARK(BM_GF2XLOGOPT_Large_Addition)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOGOPT_Large_Multiplication)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOGOPT_Large_Division)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOGOPT_Large_Inversion)->DenseRange(17, 20);
BENCHMARK(BM_GF2XLOGOPT_Large_Exponentiation)->DenseRange(17, 20);

BENCHMARK(BM_GF2XZECH_Large_Addition)->DenseRange(17, 20);
BENCHMARK(BM_GF2XZECH_Large_Multiplication)->DenseRange(17, 20);
BENCHMARK(BM_GF2XZECH_Large_Division)->DenseRange(17, 20);
BENCHMARK(BM_GF2XZECH_Large_Inversion)->DenseRange(17, 20);
BENCHMARK(BM_GF2XZECH_Large_Exponentiation)->DenseRange(17, 20);

BENCHMARK_MAIN();