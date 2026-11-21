#include <benchmark/benchmark.h>
#include <sys/resource.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "xgalois/field/gf_prime.hpp"

using xg::GaloisFieldPrime;
using xg::GaloisFieldPrimeTable;

namespace {

constexpr size_t kNumElements = 1024;
constexpr size_t kElementMask = kNumElements - 1;
constexpr uint64_t kDefaultSeed = 42;
constexpr uint64_t kDenomSeed = 1337;
constexpr uint64_t kExpSeed = 9001;

constexpr uint32_t kBarrett32Prime = 4294967291U;
constexpr uint64_t kBarrett64Prime = 2305843009213693951ULL;

struct MemoryUsage {
  size_t peak_rss_kb = 0;
};

[[maybe_unused]] MemoryUsage GetMemoryUsage() {
  MemoryUsage usage;
  struct rusage rusage_data {};
  if (getrusage(RUSAGE_SELF, &rusage_data) == 0) {
#ifdef __APPLE__
    usage.peak_rss_kb = static_cast<size_t>(rusage_data.ru_maxrss) / 1024;
#else
    usage.peak_rss_kb = static_cast<size_t>(rusage_data.ru_maxrss);
#endif
  }
  return usage;
}

// Includes zero operands naturally via seeded field->Random(gen).
template <typename FieldType>
std::vector<typename FieldType::element_type> GenerateRandomElements(
    const FieldType& field, size_t count, uint64_t seed = kDefaultSeed) {
  std::mt19937_64 gen(seed);
  std::vector<typename FieldType::element_type> elements;
  elements.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    elements.push_back(field.Random(gen));
  }
  return elements;
}

// Non-zero operands required for division denominators and multiplicative
// inversion.
template <typename FieldType>
std::vector<typename FieldType::element_type> GenerateNonZeroElements(
    const FieldType& field, size_t count, uint64_t seed = kDenomSeed) {
  std::mt19937_64 gen(seed);
  std::vector<typename FieldType::element_type> elements;
  elements.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    auto elem = field.Random(gen);
    while (field.IsZero(elem)) {
      elem = field.Random(gen);
    }
    elements.push_back(elem);
  }
  return elements;
}

std::vector<uint64_t> GenerateExponents(size_t count, uint64_t max_exp,
                                        uint64_t seed = kExpSeed) {
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<uint64_t> dist(0, max_exp);
  std::vector<uint64_t> exponents;
  exponents.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    exponents.push_back(dist(gen));
  }
  return exponents;
}

template <typename FieldType>
void RecordFieldCounters(benchmark::State& state, const FieldType& field) {
  state.counters["FieldBytes"] = static_cast<double>(field.MemoryBytes());
  state.counters["FieldOrder"] = static_cast<double>(field.Order());
}

template <typename FieldType>
void BenchAdd(benchmark::State& state, typename FieldType::element_type p) {
  auto field = std::make_shared<FieldType>(p);
  const auto elements = GenerateRandomElements(*field, kNumElements);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Add(elements[idx & kElementMask],
                             elements[(idx + 1) & kElementMask]);
    benchmark::DoNotOptimize(result);
    ++idx;
  }

  RecordFieldCounters(state, *field);
}

template <typename FieldType>
void BenchMul(benchmark::State& state, typename FieldType::element_type p) {
  auto field = std::make_shared<FieldType>(p);
  const auto elements = GenerateRandomElements(*field, kNumElements);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Mul(elements[idx & kElementMask],
                             elements[(idx + 1) & kElementMask]);
    benchmark::DoNotOptimize(result);
    ++idx;
  }

  RecordFieldCounters(state, *field);
}

template <typename FieldType>
void BenchDiv(benchmark::State& state, typename FieldType::element_type p) {
  auto field = std::make_shared<FieldType>(p);
  const auto numerators =
      GenerateRandomElements(*field, kNumElements, kDefaultSeed);
  const auto denominators =
      GenerateNonZeroElements(*field, kNumElements, kDenomSeed);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Div(numerators[idx & kElementMask],
                             denominators[idx & kElementMask]);
    benchmark::DoNotOptimize(result);
    ++idx;
  }

  RecordFieldCounters(state, *field);
}

template <typename FieldType>
void BenchInv(benchmark::State& state, typename FieldType::element_type p) {
  auto field = std::make_shared<FieldType>(p);
  const auto elements =
      GenerateNonZeroElements(*field, kNumElements, kDefaultSeed);
  size_t idx = 0;

  for (auto _ : state) {
    auto result = field->Inv(elements[idx & kElementMask]);
    benchmark::DoNotOptimize(result);
    ++idx;
  }

  RecordFieldCounters(state, *field);
}

template <typename FieldType>
void BenchPow(benchmark::State& state, typename FieldType::element_type p,
              uint64_t max_exp = 256) {
  auto field = std::make_shared<FieldType>(p);
  const auto elements = GenerateRandomElements(*field, kNumElements);
  const auto exponents = GenerateExponents(kNumElements, max_exp);
  size_t idx = 0;

  for (auto _ : state) {
    auto result =
        field->Pow(elements[idx & kElementMask], exponents[idx & kElementMask]);
    benchmark::DoNotOptimize(result);
    ++idx;
  }

  RecordFieldCounters(state, *field);
}

// Small primes: GF(2) .. GF(47)
void BM_GFP_Small_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Small_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Small_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Small_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Small_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrime<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)), 256);
}

void BM_GFPTABLE_Small_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Small_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Small_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Small_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Small_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)), 256);
}

// Medium primes: GF(53) .. GF(127)
void BM_GFP_Medium_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Medium_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Medium_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Medium_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Medium_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrime<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)), 256);
}

void BM_GFPTABLE_Medium_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Medium_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Medium_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Medium_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Medium_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)), 256);
}

// Large primes: GF(131) .. GF(997)
void BM_GFP_Large_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Large_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Large_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Large_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrime<uint32_t>>(state,
                                       static_cast<uint32_t>(state.range(0)));
}
void BM_GFP_Large_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrime<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)), 256);
}

void BM_GFPTABLE_Large_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Large_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Large_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Large_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)));
}
void BM_GFPTABLE_Large_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrimeTable<uint32_t>>(
      state, static_cast<uint32_t>(state.range(0)), 256);
}

// 32-bit Barrett prime: p = 4294967291 (2^32 - 5)
void BM_GFP_Barrett32_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrime<uint32_t>>(state, kBarrett32Prime);
}
void BM_GFP_Barrett32_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrime<uint32_t>>(state, kBarrett32Prime);
}
void BM_GFP_Barrett32_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrime<uint32_t>>(state, kBarrett32Prime);
}
void BM_GFP_Barrett32_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrime<uint32_t>>(state, kBarrett32Prime);
}
void BM_GFP_Barrett32_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrime<uint32_t>>(state, kBarrett32Prime, 4096);
}

// 64-bit Barrett prime: p = 2305843009213693951 (2^61 - 1)
void BM_GFP_Barrett64_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldPrime<uint64_t>>(state, kBarrett64Prime);
}
void BM_GFP_Barrett64_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldPrime<uint64_t>>(state, kBarrett64Prime);
}
void BM_GFP_Barrett64_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldPrime<uint64_t>>(state, kBarrett64Prime);
}
void BM_GFP_Barrett64_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldPrime<uint64_t>>(state, kBarrett64Prime);
}
void BM_GFP_Barrett64_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldPrime<uint64_t>>(state, kBarrett64Prime, 4096);
}

void ApplySmallPrimes(benchmark::internal::Benchmark* b) {
  for (int64_t p : {2, 3, 7, 13, 23, 31, 47}) {
    b->Arg(p);
  }
  b->Unit(benchmark::kNanosecond);
}

void ApplyMediumPrimes(benchmark::internal::Benchmark* b) {
  for (int64_t p : {53, 71, 89, 103, 127}) {
    b->Arg(p);
  }
  b->Unit(benchmark::kNanosecond);
}

void ApplyLargePrimes(benchmark::internal::Benchmark* b) {
  for (int64_t p : {131, 179, 241, 307, 389, 463, 577, 683, 787, 907, 997}) {
    b->Arg(p);
  }
  b->Unit(benchmark::kNanosecond);
}

}  // namespace

BENCHMARK(BM_GFP_Small_Addition)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFP_Small_Multiplication)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFP_Small_Division)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFP_Small_Inversion)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFP_Small_Exponentiation)->Apply(ApplySmallPrimes);

BENCHMARK(BM_GFPTABLE_Small_Addition)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFPTABLE_Small_Multiplication)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFPTABLE_Small_Division)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFPTABLE_Small_Inversion)->Apply(ApplySmallPrimes);
BENCHMARK(BM_GFPTABLE_Small_Exponentiation)->Apply(ApplySmallPrimes);

BENCHMARK(BM_GFP_Medium_Addition)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFP_Medium_Multiplication)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFP_Medium_Division)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFP_Medium_Inversion)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFP_Medium_Exponentiation)->Apply(ApplyMediumPrimes);

BENCHMARK(BM_GFPTABLE_Medium_Addition)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFPTABLE_Medium_Multiplication)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFPTABLE_Medium_Division)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFPTABLE_Medium_Inversion)->Apply(ApplyMediumPrimes);
BENCHMARK(BM_GFPTABLE_Medium_Exponentiation)->Apply(ApplyMediumPrimes);

BENCHMARK(BM_GFP_Large_Addition)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFP_Large_Multiplication)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFP_Large_Division)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFP_Large_Inversion)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFP_Large_Exponentiation)->Apply(ApplyLargePrimes);

BENCHMARK(BM_GFPTABLE_Large_Addition)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFPTABLE_Large_Multiplication)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFPTABLE_Large_Division)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFPTABLE_Large_Inversion)->Apply(ApplyLargePrimes);
BENCHMARK(BM_GFPTABLE_Large_Exponentiation)->Apply(ApplyLargePrimes);

BENCHMARK(BM_GFP_Barrett32_Addition)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett32_Multiplication)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett32_Division)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett32_Inversion)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett32_Exponentiation)->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GFP_Barrett64_Addition)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett64_Multiplication)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett64_Division)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett64_Inversion)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFP_Barrett64_Exponentiation)->Unit(benchmark::kNanosecond);
