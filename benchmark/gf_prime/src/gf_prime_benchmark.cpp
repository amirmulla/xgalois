/**
 * @file gf_prime_benchmark.cpp
 * @brief Comprehensive benchmark comparison between prime field implementations
 * Compares 2 different implementation classes (standard and logarithm) across
 * small, medium, and large field sizes
 */

#include <benchmark/benchmark.h>
#include <sys/resource.h>
#include <unistd.h>

#include <chrono>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"

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
// Prime Field Sizes for Testing
//------------------------------------------------------------------------------

// Small prime fields - Common small primes (reduced sample)
const std::vector<uint32_t> SMALL_PRIMES = {2, 3, 7, 13, 23, 31, 47};

// Medium prime fields - Practical sizes for applications (reduced sample)
const std::vector<uint32_t> MEDIUM_PRIMES = {53, 71, 89, 103, 127};

// Large prime fields - Larger primes for testing scalability (reduced sample)
const std::vector<uint32_t> LARGE_PRIMES = {131, 179, 241, 307, 389, 463,
                                            577, 683, 787, 907, 997};

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
// Small Prime Field Benchmarks (Small primes: 2, 3, 5, 7, 11, 13, 17, 19, 23,
// 29, 31, 37, 41, 43, 47)
//------------------------------------------------------------------------------

// GaloisFieldPrime - Standard Implementation
static void BM_GFP_Small_Addition(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

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

static void BM_GFP_Small_Multiplication(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Small_Division(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Small_Inversion(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Small_Exponentiation(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

// GaloisFieldPrimeTable - Logarithm Table Implementation
static void BM_GFPTABLE_Small_Addition(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

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

static void BM_GFPTABLE_Small_Multiplication(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Small_Division(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Small_Inversion(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Small_Exponentiation(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

//------------------------------------------------------------------------------
// Medium Prime Field Benchmarks
//------------------------------------------------------------------------------

// GaloisFieldPrime - Standard Implementation
static void BM_GFP_Medium_Addition(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Medium_Multiplication(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Medium_Division(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Medium_Inversion(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Medium_Exponentiation(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

// GaloisFieldPrimeTable - Logarithm Table Implementation
static void BM_GFPTABLE_Medium_Addition(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Medium_Multiplication(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Medium_Division(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Medium_Inversion(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Medium_Exponentiation(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

//------------------------------------------------------------------------------
// Large Prime Field Benchmarks
//------------------------------------------------------------------------------

// GaloisFieldPrime - Standard Implementation
static void BM_GFP_Large_Addition(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Large_Multiplication(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Large_Division(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Large_Inversion(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFP_Large_Exponentiation(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrime<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

// GaloisFieldPrimeTable - Logarithm Table Implementation
static void BM_GFPTABLE_Large_Addition(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Large_Multiplication(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Large_Division(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(elements[idx % elements.size()],
                             elements[(idx + 1) % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Large_Inversion(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx % elements.size()]);
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

static void BM_GFPTABLE_Large_Exponentiation(benchmark::State& state) {
  uint32_t p = static_cast<uint32_t>(state.range(0));
  auto field = std::make_shared<GaloisFieldPrimeTable<uint32_t>>(p);

  auto elements = GenerateRandomElements(field, 1000);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Pow(elements[idx % elements.size()],
                             static_cast<uint64_t>(idx % 256 + 1));
    benchmark::DoNotOptimize(result);
    idx++;
  }

  state.counters["FieldOrder"] = field->Order();
}

//------------------------------------------------------------------------------
// Benchmark Registration
//------------------------------------------------------------------------------

// Small prime field benchmarks - using specific small primes
BENCHMARK(BM_GFP_Small_Addition)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Small_Multiplication)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Small_Division)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Small_Inversion)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Small_Exponentiation)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GFPTABLE_Small_Addition)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Small_Multiplication)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Small_Division)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Small_Inversion)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Small_Exponentiation)
    ->Args({2})
    ->Args({3})
    ->Args({7})
    ->Args({13})
    ->Args({23})
    ->Args({31})
    ->Args({47})
    ->Unit(benchmark::kNanosecond);

// Medium prime field benchmarks - using specific medium primes
BENCHMARK(BM_GFP_Medium_Addition)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Medium_Multiplication)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Medium_Division)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Medium_Inversion)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Medium_Exponentiation)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GFPTABLE_Medium_Addition)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Medium_Multiplication)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Medium_Division)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Medium_Inversion)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Medium_Exponentiation)
    ->Args({53})
    ->Args({71})
    ->Args({89})
    ->Args({103})
    ->Args({127})
    ->Unit(benchmark::kNanosecond);

// Large prime field benchmarks - using specific large primes
BENCHMARK(BM_GFP_Large_Addition)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Large_Multiplication)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Large_Division)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Large_Inversion)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Large_Exponentiation)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GFPTABLE_Large_Addition)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Large_Multiplication)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Large_Division)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Large_Inversion)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFPTABLE_Large_Exponentiation)
    ->Args({131})
    ->Args({179})
    ->Args({241})
    ->Args({307})
    ->Args({389})
    ->Args({463})
    ->Args({577})
    ->Args({683})
    ->Args({787})
    ->Args({907})
    ->Args({997})
    ->Unit(benchmark::kNanosecond);

BENCHMARK_MAIN();
