#include <benchmark/benchmark.h>
#include <sys/resource.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "xgalois/field/gf_binary.hpp"

using xg::GaloisFieldBinaryExtension;
using xg::GF2X;
using xg::GF2XLOG;
using xg::GF2XZECH;
using xg::GFBELogTables;
using xg::GFBELogTablesOpt;
using xg::GFBEZechLogTables;

namespace {

constexpr size_t kNumElements = 1024;
constexpr size_t kElementMask = kNumElements - 1;
constexpr uint64_t kDefaultSeed = 42;
constexpr uint64_t kDenomSeed = 1337;
constexpr uint64_t kExpSeed = 9001;

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

// Includes zero operands naturally via seeded field.Random(gen) (works
// uniformly for standard polynomial representations and Zech-log kLogZero).
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
void BenchAdd(benchmark::State& state, uint8_t m) {
  auto field = std::make_shared<FieldType>(m);
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
void BenchMul(benchmark::State& state, uint8_t m) {
  auto field = std::make_shared<FieldType>(m);
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
void BenchDiv(benchmark::State& state, uint8_t m) {
  auto field = std::make_shared<FieldType>(m);
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
void BenchInv(benchmark::State& state, uint8_t m) {
  auto field = std::make_shared<FieldType>(m);
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
void BenchPow(benchmark::State& state, uint8_t m, uint64_t max_exp) {
  auto field = std::make_shared<FieldType>(m);
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

// ============================================================================
// Small fields: GF(2^2) .. GF(2^8) (uint32_t)
// ============================================================================

void BM_GF2X_Small_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Small_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Small_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Small_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Small_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 256);
}

void BM_GF2XLOG_Small_Addition(benchmark::State& state) {
  BenchAdd<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Small_Multiplication(benchmark::State& state) {
  BenchMul<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Small_Division(benchmark::State& state) {
  BenchDiv<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Small_Inversion(benchmark::State& state) {
  BenchInv<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Small_Exponentiation(benchmark::State& state) {
  BenchPow<GFBELogTables<uint32_t>>(state, static_cast<uint8_t>(state.range(0)),
                                    256);
}

void BM_GF2XLOGOPT_Small_Addition(benchmark::State& state) {
  BenchAdd<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Small_Multiplication(benchmark::State& state) {
  BenchMul<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Small_Division(benchmark::State& state) {
  BenchDiv<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Small_Inversion(benchmark::State& state) {
  BenchInv<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Small_Exponentiation(benchmark::State& state) {
  BenchPow<GFBELogTablesOpt<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 256);
}

void BM_GF2XZECH_Small_Addition(benchmark::State& state) {
  BenchAdd<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Small_Multiplication(benchmark::State& state) {
  BenchMul<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Small_Division(benchmark::State& state) {
  BenchDiv<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Small_Inversion(benchmark::State& state) {
  BenchInv<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Small_Exponentiation(benchmark::State& state) {
  BenchPow<GFBEZechLogTables<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 256);
}

// ============================================================================
// Medium fields: GF(2^9) .. GF(2^16) (uint32_t)
// ============================================================================

void BM_GF2X_Medium_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Medium_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Medium_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Medium_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Medium_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 1024);
}

void BM_GF2XLOG_Medium_Addition(benchmark::State& state) {
  BenchAdd<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Medium_Multiplication(benchmark::State& state) {
  BenchMul<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Medium_Division(benchmark::State& state) {
  BenchDiv<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Medium_Inversion(benchmark::State& state) {
  BenchInv<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Medium_Exponentiation(benchmark::State& state) {
  BenchPow<GFBELogTables<uint32_t>>(state, static_cast<uint8_t>(state.range(0)),
                                    1024);
}

void BM_GF2XLOGOPT_Medium_Addition(benchmark::State& state) {
  BenchAdd<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Medium_Multiplication(benchmark::State& state) {
  BenchMul<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Medium_Division(benchmark::State& state) {
  BenchDiv<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Medium_Inversion(benchmark::State& state) {
  BenchInv<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Medium_Exponentiation(benchmark::State& state) {
  BenchPow<GFBELogTablesOpt<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 1024);
}

void BM_GF2XZECH_Medium_Addition(benchmark::State& state) {
  BenchAdd<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Medium_Multiplication(benchmark::State& state) {
  BenchMul<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Medium_Division(benchmark::State& state) {
  BenchDiv<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Medium_Inversion(benchmark::State& state) {
  BenchInv<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Medium_Exponentiation(benchmark::State& state) {
  BenchPow<GFBEZechLogTables<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 1024);
}

// ============================================================================
// Large fields: GF(2^17) .. GF(2^20) (matching uint32_t across all 4)
// ============================================================================

void BM_GF2X_Large_Addition(benchmark::State& state) {
  BenchAdd<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Large_Multiplication(benchmark::State& state) {
  BenchMul<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Large_Division(benchmark::State& state) {
  BenchDiv<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Large_Inversion(benchmark::State& state) {
  BenchInv<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)));
}
void BM_GF2X_Large_Exponentiation(benchmark::State& state) {
  BenchPow<GaloisFieldBinaryExtension<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 4096);
}

void BM_GF2XLOG_Large_Addition(benchmark::State& state) {
  BenchAdd<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Large_Multiplication(benchmark::State& state) {
  BenchMul<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Large_Division(benchmark::State& state) {
  BenchDiv<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Large_Inversion(benchmark::State& state) {
  BenchInv<GFBELogTables<uint32_t>>(state,
                                    static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOG_Large_Exponentiation(benchmark::State& state) {
  BenchPow<GFBELogTables<uint32_t>>(state, static_cast<uint8_t>(state.range(0)),
                                    4096);
}

void BM_GF2XLOGOPT_Large_Addition(benchmark::State& state) {
  BenchAdd<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Large_Multiplication(benchmark::State& state) {
  BenchMul<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Large_Division(benchmark::State& state) {
  BenchDiv<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Large_Inversion(benchmark::State& state) {
  BenchInv<GFBELogTablesOpt<uint32_t>>(state,
                                       static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XLOGOPT_Large_Exponentiation(benchmark::State& state) {
  BenchPow<GFBELogTablesOpt<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 4096);
}

void BM_GF2XZECH_Large_Addition(benchmark::State& state) {
  BenchAdd<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Large_Multiplication(benchmark::State& state) {
  BenchMul<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Large_Division(benchmark::State& state) {
  BenchDiv<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Large_Inversion(benchmark::State& state) {
  BenchInv<GFBEZechLogTables<uint32_t>>(state,
                                        static_cast<uint8_t>(state.range(0)));
}
void BM_GF2XZECH_Large_Exponentiation(benchmark::State& state) {
  BenchPow<GFBEZechLogTables<uint32_t>>(
      state, static_cast<uint8_t>(state.range(0)), 4096);
}

// ============================================================================
// GF(2^8) uint8_t comparison across GF2X, GF2XLOG, GF2XZECH
// ============================================================================

void BM_GF2X_U8_GF256_Addition(benchmark::State& state) {
  BenchAdd<GF2X<uint8_t>>(state, 8);
}
void BM_GF2X_U8_GF256_Multiplication(benchmark::State& state) {
  BenchMul<GF2X<uint8_t>>(state, 8);
}
void BM_GF2X_U8_GF256_Division(benchmark::State& state) {
  BenchDiv<GF2X<uint8_t>>(state, 8);
}
void BM_GF2X_U8_GF256_Inversion(benchmark::State& state) {
  BenchInv<GF2X<uint8_t>>(state, 8);
}
void BM_GF2X_U8_GF256_Exponentiation(benchmark::State& state) {
  BenchPow<GF2X<uint8_t>>(state, 8, 256);
}

void BM_GF2XLOG_U8_GF256_Addition(benchmark::State& state) {
  BenchAdd<GF2XLOG<uint8_t>>(state, 8);
}
void BM_GF2XLOG_U8_GF256_Multiplication(benchmark::State& state) {
  BenchMul<GF2XLOG<uint8_t>>(state, 8);
}
void BM_GF2XLOG_U8_GF256_Division(benchmark::State& state) {
  BenchDiv<GF2XLOG<uint8_t>>(state, 8);
}
void BM_GF2XLOG_U8_GF256_Inversion(benchmark::State& state) {
  BenchInv<GF2XLOG<uint8_t>>(state, 8);
}
void BM_GF2XLOG_U8_GF256_Exponentiation(benchmark::State& state) {
  BenchPow<GF2XLOG<uint8_t>>(state, 8, 256);
}

void BM_GF2XZECH_U8_GF256_Addition(benchmark::State& state) {
  BenchAdd<GF2XZECH<uint8_t>>(state, 8);
}
void BM_GF2XZECH_U8_GF256_Multiplication(benchmark::State& state) {
  BenchMul<GF2XZECH<uint8_t>>(state, 8);
}
void BM_GF2XZECH_U8_GF256_Division(benchmark::State& state) {
  BenchDiv<GF2XZECH<uint8_t>>(state, 8);
}
void BM_GF2XZECH_U8_GF256_Inversion(benchmark::State& state) {
  BenchInv<GF2XZECH<uint8_t>>(state, 8);
}
void BM_GF2XZECH_U8_GF256_Exponentiation(benchmark::State& state) {
  BenchPow<GF2XZECH<uint8_t>>(state, 8, 256);
}

}  // namespace

BENCHMARK(BM_GF2X_Small_Addition)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Small_Multiplication)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Small_Division)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Small_Inversion)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Small_Exponentiation)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOG_Small_Addition)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Small_Multiplication)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Small_Division)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Small_Inversion)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Small_Exponentiation)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOGOPT_Small_Addition)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Small_Multiplication)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Small_Division)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Small_Inversion)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Small_Exponentiation)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XZECH_Small_Addition)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Small_Multiplication)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Small_Division)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Small_Inversion)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Small_Exponentiation)
    ->DenseRange(2, 8)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2X_Medium_Addition)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Medium_Multiplication)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Medium_Division)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Medium_Inversion)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Medium_Exponentiation)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOG_Medium_Addition)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Medium_Multiplication)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Medium_Division)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Medium_Inversion)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Medium_Exponentiation)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOGOPT_Medium_Addition)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Medium_Multiplication)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Medium_Division)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Medium_Inversion)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Medium_Exponentiation)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XZECH_Medium_Addition)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Medium_Multiplication)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Medium_Division)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Medium_Inversion)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Medium_Exponentiation)
    ->DenseRange(9, 16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2X_Large_Addition)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Large_Multiplication)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Large_Division)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Large_Inversion)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_Large_Exponentiation)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOG_Large_Addition)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Large_Multiplication)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Large_Division)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Large_Inversion)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_Large_Exponentiation)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOGOPT_Large_Addition)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Large_Multiplication)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Large_Division)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Large_Inversion)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOGOPT_Large_Exponentiation)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XZECH_Large_Addition)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Large_Multiplication)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Large_Division)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Large_Inversion)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_Large_Exponentiation)
    ->DenseRange(17, 20)
    ->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2X_U8_GF256_Addition)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_U8_GF256_Multiplication)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_U8_GF256_Division)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_U8_GF256_Inversion)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2X_U8_GF256_Exponentiation)->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XLOG_U8_GF256_Addition)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_U8_GF256_Multiplication)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_U8_GF256_Division)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_U8_GF256_Inversion)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XLOG_U8_GF256_Exponentiation)->Unit(benchmark::kNanosecond);

BENCHMARK(BM_GF2XZECH_U8_GF256_Addition)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_U8_GF256_Multiplication)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_U8_GF256_Division)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_U8_GF256_Inversion)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF2XZECH_U8_GF256_Exponentiation)->Unit(benchmark::kNanosecond);