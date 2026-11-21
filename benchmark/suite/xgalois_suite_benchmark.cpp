#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/decoder/sugiyama_decoder.hpp"
#include "xgalois/coding/encoder/grs_encoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"

using xg::GaloisFieldElement;
using xg::GaloisFieldPrime;
using xg::GaloisFieldPrimeTable;
using xg::GF2;
using xg::GF2X;
using xg::GF2XLOG;
using xg::GF2XZECH;
using xg::GFArray;
using xg::PolynomialDense;

namespace {

constexpr uint64_t kSeed = 42;
constexpr size_t kScalarBufSize = 1024;
constexpr size_t kScalarBufMask = kScalarBufSize - 1;
constexpr uint32_t kPrime65521 = 65521U;

template <typename FieldType>
PolynomialDense<FieldType> RandomPolynomial(
    const std::shared_ptr<FieldType>& field, size_t degree,
    std::mt19937_64& rng) {
  using RawType = typename FieldType::element_type;
  std::vector<RawType> coeffs(degree + 1);
  for (size_t i = 0; i <= degree; ++i) {
    coeffs[i] = field->Random(rng);
  }
  while (field->IsZero(coeffs[degree])) {
    coeffs[degree] = field->Random(rng);
  }
  return PolynomialDense<FieldType>::FromRaw(std::move(coeffs), field);
}

// ============================================================================
// 1. GFArray Bulk Kernels: AddArray, MulArray, FmaArray
// ============================================================================

template <typename FieldType, typename FieldFactory>
void BenchBulkAddArray(benchmark::State& state, FieldFactory make_field) {
  using T = typename FieldType::element_type;
  const size_t n = static_cast<size_t>(state.range(0));
  auto field = make_field();
  std::mt19937_64 rng(kSeed);

  const auto a = GFArray<FieldType>::Random(field, {n}, rng);
  const auto b = GFArray<FieldType>::Random(field, {n}, rng);
  auto dst = GFArray<FieldType>::Zeros(field, {n});

  for (auto _ : state) {
    field->AddArray(dst.RawSpan(), a.RawSpan(), b.RawSpan());
    benchmark::DoNotOptimize(dst.RawSpan().data());
    benchmark::ClobberMemory();
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n * sizeof(T)));
}

template <typename FieldType, typename FieldFactory>
void BenchBulkMulArray(benchmark::State& state, FieldFactory make_field) {
  using T = typename FieldType::element_type;
  const size_t n = static_cast<size_t>(state.range(0));
  auto field = make_field();
  std::mt19937_64 rng(kSeed);

  const auto a = GFArray<FieldType>::Random(field, {n}, rng);
  const auto b = GFArray<FieldType>::Random(field, {n}, rng);
  auto dst = GFArray<FieldType>::Zeros(field, {n});

  for (auto _ : state) {
    field->MulArray(dst.RawSpan(), a.RawSpan(), b.RawSpan());
    benchmark::DoNotOptimize(dst.RawSpan().data());
    benchmark::ClobberMemory();
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n * sizeof(T)));
}

template <typename FieldType, typename FieldFactory>
void BenchBulkFmaArray(benchmark::State& state, FieldFactory make_field,
                       typename FieldType::element_type scalar) {
  using T = typename FieldType::element_type;
  const size_t n = static_cast<size_t>(state.range(0));
  auto field = make_field();
  std::mt19937_64 rng(kSeed);

  const auto a = GFArray<FieldType>::Random(field, {n}, rng);
  auto dst = GFArray<FieldType>::Random(field, {n}, rng);

  for (auto _ : state) {
    field->FmaArray(dst.RawSpan(), a.RawSpan(), scalar);
    benchmark::DoNotOptimize(dst.RawSpan().data());
    benchmark::ClobberMemory();
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n * sizeof(T)));
}

void BM_GFArray_GF256_U8_AddArray(benchmark::State& state) {
  BenchBulkAddArray<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); });
}

void BM_GFArray_GF256_U8_MulArray(benchmark::State& state) {
  BenchBulkMulArray<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); });
}

void BM_GFArray_GF256_U8_FmaArray(benchmark::State& state) {
  BenchBulkFmaArray<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); },
      uint8_t{0x53});
}

void BM_GFArray_Prime65521_U32_AddArray(benchmark::State& state) {
  BenchBulkAddArray<GaloisFieldPrime<uint32_t>>(state, [] {
    return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521);
  });
}

void BM_GFArray_Prime65521_U32_MulArray(benchmark::State& state) {
  BenchBulkMulArray<GaloisFieldPrime<uint32_t>>(state, [] {
    return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521);
  });
}

void BM_GFArray_Prime65521_U32_FmaArray(benchmark::State& state) {
  BenchBulkFmaArray<GaloisFieldPrime<uint32_t>>(
      state,
      [] { return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521); },
      uint32_t{12345U});
}

// ============================================================================
// 2. GaloisFieldElement vs Raw Field Op Overhead Comparison
// ============================================================================

template <typename FieldType, typename FieldFactory>
void BenchRawMul(benchmark::State& state, FieldFactory make_field) {
  auto field = make_field();
  std::mt19937_64 rng(kSeed);
  std::vector<typename FieldType::element_type> vals(kScalarBufSize);
  for (auto& v : vals) {
    v = field->Random(rng);
  }
  size_t idx = 0;

  for (auto _ : state) {
    auto res = field->Mul(vals[idx & kScalarBufMask],
                          vals[(idx + 1) & kScalarBufMask]);
    benchmark::DoNotOptimize(res);
    ++idx;
  }
}

template <typename FieldType, typename FieldFactory>
void BenchElementMul(benchmark::State& state, FieldFactory make_field) {
  auto field = make_field();
  std::mt19937_64 rng(kSeed);
  std::vector<GaloisFieldElement<FieldType>> elems;
  elems.reserve(kScalarBufSize);
  for (size_t i = 0; i < kScalarBufSize; ++i) {
    elems.push_back(
        GaloisFieldElement<FieldType>::FromInternal(field->Random(rng), field));
  }
  size_t idx = 0;

  for (auto _ : state) {
    auto res = elems[idx & kScalarBufMask] * elems[(idx + 1) & kScalarBufMask];
    auto raw = res.RawValue();
    benchmark::DoNotOptimize(raw);
    ++idx;
  }
}

template <typename FieldType, typename FieldFactory>
void BenchElementMulAssign(benchmark::State& state, FieldFactory make_field) {
  auto field = make_field();
  std::mt19937_64 rng(kSeed);
  std::vector<GaloisFieldElement<FieldType>> elems;
  elems.reserve(kScalarBufSize);
  for (size_t i = 0; i < kScalarBufSize; ++i) {
    auto v = field->Random(rng);
    while (field->IsZero(v)) {
      v = field->Random(rng);
    }
    elems.push_back(GaloisFieldElement<FieldType>::FromInternal(v, field));
  }
  auto acc = GaloisFieldElement<FieldType>::FromInternal(
      field->MultiplicativeIdentity(), field);
  size_t idx = 0;

  for (auto _ : state) {
    acc *= elems[idx & kScalarBufMask];
    auto raw = acc.RawValue();
    benchmark::DoNotOptimize(raw);
    ++idx;
  }
}

void BM_Overhead_GF256_RawMul(benchmark::State& state) {
  BenchRawMul<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); });
}

void BM_Overhead_GF256_ElementMul(benchmark::State& state) {
  BenchElementMul<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); });
}

void BM_Overhead_GF256_ElementMulAssign(benchmark::State& state) {
  BenchElementMulAssign<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); });
}

void BM_Overhead_Prime65521_RawMul(benchmark::State& state) {
  BenchRawMul<GaloisFieldPrime<uint32_t>>(state, [] {
    return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521);
  });
}

void BM_Overhead_Prime65521_ElementMul(benchmark::State& state) {
  BenchElementMul<GaloisFieldPrime<uint32_t>>(state, [] {
    return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521);
  });
}

void BM_Overhead_Prime65521_ElementMulAssign(benchmark::State& state) {
  BenchElementMulAssign<GaloisFieldPrime<uint32_t>>(state, [] {
    return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521);
  });
}

// ============================================================================
// 3. PolynomialDense: Schoolbook (deg 16) vs Karatsuba (deg 64, 256) & DivRem
// ============================================================================

template <typename FieldType, typename FieldFactory>
void BenchPolyMul(benchmark::State& state, FieldFactory make_field) {
  const size_t deg = static_cast<size_t>(state.range(0));
  auto field = make_field();
  std::mt19937_64 rng(kSeed);

  const auto p1 = RandomPolynomial(field, deg, rng);
  const auto p2 = RandomPolynomial(field, deg, rng);

  for (auto _ : state) {
    auto prod = p1 * p2;
    benchmark::DoNotOptimize(prod.RawCoeffs().data());
    benchmark::ClobberMemory();
  }
}

template <typename FieldType, typename FieldFactory>
void BenchPolyDivRem(benchmark::State& state, FieldFactory make_field,
                     size_t deg_num, size_t deg_den) {
  auto field = make_field();
  std::mt19937_64 rng(kSeed);

  const auto num = RandomPolynomial(field, deg_num, rng);
  const auto den = RandomPolynomial(field, deg_den, rng);

  for (auto _ : state) {
    auto qr = num.DivRem(den);
    benchmark::DoNotOptimize(qr.first.RawCoeffs().data());
    benchmark::DoNotOptimize(qr.second.RawCoeffs().data());
    benchmark::ClobberMemory();
  }
}

void BM_PolyDense_GF256_Mul(benchmark::State& state) {
  BenchPolyMul<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); });
}

void BM_PolyDense_Prime65521_Mul(benchmark::State& state) {
  BenchPolyMul<GaloisFieldPrime<uint32_t>>(state, [] {
    return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521);
  });
}

void BM_PolyDense_GF256_DivRem_256_128(benchmark::State& state) {
  BenchPolyDivRem<GF2XLOG<uint8_t>>(
      state, [] { return std::make_shared<GF2XLOG<uint8_t>>(8); }, 256, 128);
}

void BM_PolyDense_Prime65521_DivRem_256_128(benchmark::State& state) {
  BenchPolyDivRem<GaloisFieldPrime<uint32_t>>(
      state,
      [] { return std::make_shared<GaloisFieldPrime<uint32_t>>(kPrime65521); },
      256, 128);
}

// ============================================================================
// 4. Linear Algebra: linalg::dot & linalg::rref on GF(2) M4RI and GF(2^8)
// ============================================================================

void BM_Linalg_GF2_M4RI_Dot(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf2 = std::make_shared<GF2>();
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::random_matrix(gf2, {n, n}, rng);
  const auto b = xg::linalg::random_matrix(gf2, {n, n}, rng);

  for (auto _ : state) {
    auto c = xg::linalg::dot(a, b);
    benchmark::DoNotOptimize(c.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

void BM_Linalg_GF2_M4RI_Rref(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf2 = std::make_shared<GF2>();
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::random_matrix(gf2, {n, n}, rng);

  for (auto _ : state) {
    auto r = xg::linalg::rref(a);
    benchmark::DoNotOptimize(r.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

void BM_Linalg_GF256_Dot(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::random_matrix(gf256, {n, n}, rng);
  const auto b = xg::linalg::random_matrix(gf256, {n, n}, rng);

  for (auto _ : state) {
    auto c = xg::linalg::dot(a, b);
    benchmark::DoNotOptimize(c.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

void BM_Linalg_GF256_Rref(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::random_matrix(gf256, {n, n}, rng);

  for (auto _ : state) {
    auto r = xg::linalg::rref(a);
    benchmark::DoNotOptimize(r.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

// ============================================================================
// 5. RS(255, 223) over GF2XLOG<uint8_t>(8): Encode, BM (16 err), Sugiyama (16
//    err)
// ============================================================================

struct Rs255Fixture {
  using Field = GF2XLOG<uint8_t>;
  using Code = xg::coding::GeneralizedReedSolomonCode<Field>;

  static constexpr size_t kCodeLen = 255;
  static constexpr size_t kMessageLen = 223;
  static constexpr size_t kNumErrors = 16;
  static constexpr size_t kBatchSize = 16;
  static constexpr size_t kBatchMask = kBatchSize - 1;

  std::shared_ptr<Field> field;
  Code rs;
  std::vector<GFArray<Field>> messages;
  std::vector<GFArray<Field>> codewords;
  std::vector<GFArray<Field>> received_16err;

  Rs255Fixture()
      : field(std::make_shared<Field>(8)),
        rs(Code::ReedSolomon(field, kCodeLen, kMessageLen)) {
    std::mt19937_64 rng(kSeed);
    xg::coding::GRSEncoder<Field> encoder(rs);

    std::vector<size_t> positions(kCodeLen);
    std::iota(positions.begin(), positions.end(), size_t{0});
    std::uniform_int_distribution<unsigned int> err_dist(1U, 255U);

    messages.reserve(kBatchSize);
    codewords.reserve(kBatchSize);
    received_16err.reserve(kBatchSize);

    for (size_t b = 0; b < kBatchSize; ++b) {
      auto msg = GFArray<Field>::Random(field, {kMessageLen}, rng);
      auto cw = encoder.Encode(msg);
      auto rx = cw;

      std::shuffle(positions.begin(), positions.end(), rng);
      for (size_t e = 0; e < kNumErrors; ++e) {
        const size_t pos = positions[e];
        const uint8_t err = static_cast<uint8_t>(err_dist(rng));
        rx.Raw(pos) = field->Add(rx.Raw(pos), err);
      }

      messages.push_back(std::move(msg));
      codewords.push_back(std::move(cw));
      received_16err.push_back(std::move(rx));
    }
  }
};

const Rs255Fixture& GetRs255Fixture() {
  static const Rs255Fixture fixture;
  return fixture;
}

void BM_RS255_223_Encode(benchmark::State& state) {
  const auto& fix = GetRs255Fixture();
  xg::coding::GRSEncoder<Rs255Fixture::Field> encoder(fix.rs);
  size_t idx = 0;

  for (auto _ : state) {
    auto cw = encoder.Encode(fix.messages[idx & Rs255Fixture::kBatchMask]);
    benchmark::DoNotOptimize(cw.RawSpan().data());
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(Rs255Fixture::kCodeLen));
}

void BM_RS255_223_BerlekampMassey_16Errors(benchmark::State& state) {
  const auto& fix = GetRs255Fixture();
  xg::coding::BerlekampMasseyDecoder<Rs255Fixture::Field> decoder(fix.rs);

  const auto check = decoder.DecodeWithResult(fix.received_16err[0]);
  if (!check.Ok() || check.num_errors != Rs255Fixture::kNumErrors ||
      check.codeword != fix.codewords[0]) {
    state.SkipWithError("BerlekampMasseyDecoder 16-error self-check failed");
    return;
  }

  size_t idx = 0;
  for (auto _ : state) {
    auto res = decoder.DecodeWithResult(
        fix.received_16err[idx & Rs255Fixture::kBatchMask]);
    benchmark::DoNotOptimize(res.codeword.RawSpan().data());
    benchmark::DoNotOptimize(res.num_errors);
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(Rs255Fixture::kCodeLen));
}

void BM_RS255_223_Sugiyama_16Errors(benchmark::State& state) {
  const auto& fix = GetRs255Fixture();
  xg::coding::SugiyamaDecoder<Rs255Fixture::Field> decoder(fix.rs);

  const auto check = decoder.DecodeWithResult(fix.received_16err[0]);
  if (!check.Ok() || check.num_errors != Rs255Fixture::kNumErrors ||
      check.codeword != fix.codewords[0]) {
    state.SkipWithError("SugiyamaDecoder 16-error self-check failed");
    return;
  }

  size_t idx = 0;
  for (auto _ : state) {
    auto res = decoder.DecodeWithResult(
        fix.received_16err[idx & Rs255Fixture::kBatchMask]);
    benchmark::DoNotOptimize(res.codeword.RawSpan().data());
    benchmark::DoNotOptimize(res.num_errors);
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(Rs255Fixture::kCodeLen));
}

// ============================================================================
// 6. Field Construction Time: GF(2^8)..GF(2^16) and GF(65521)
// ============================================================================

template <typename FieldType>
void BenchBinaryFieldConstruct(benchmark::State& state) {
  const uint8_t m = static_cast<uint8_t>(state.range(0));
  size_t field_bytes = 0;
  for (auto _ : state) {
    FieldType field(m);
    field_bytes = field.MemoryBytes();
    auto* ptr = &field;
    benchmark::DoNotOptimize(ptr);
    benchmark::ClobberMemory();
  }
  state.counters["FieldBytes"] = static_cast<double>(field_bytes);
}

void BM_FieldConstruct_GF2X(benchmark::State& state) {
  BenchBinaryFieldConstruct<GF2X<uint32_t>>(state);
}

void BM_FieldConstruct_GF2XLOG(benchmark::State& state) {
  BenchBinaryFieldConstruct<GF2XLOG<uint32_t>>(state);
}

void BM_FieldConstruct_GF2XZECH(benchmark::State& state) {
  BenchBinaryFieldConstruct<GF2XZECH<uint32_t>>(state);
}

void BM_FieldConstruct_GFP_65521(benchmark::State& state) {
  size_t field_bytes = 0;
  for (auto _ : state) {
    GaloisFieldPrime<uint32_t> field(kPrime65521);
    field_bytes = field.MemoryBytes();
    auto* ptr = &field;
    benchmark::DoNotOptimize(ptr);
    benchmark::ClobberMemory();
  }
  state.counters["FieldBytes"] = static_cast<double>(field_bytes);
}

void BM_FieldConstruct_GFPTABLE_65521(benchmark::State& state) {
  size_t field_bytes = 0;
  for (auto _ : state) {
    GaloisFieldPrimeTable<uint32_t> field(kPrime65521);
    field_bytes = field.MemoryBytes();
    auto* ptr = &field;
    benchmark::DoNotOptimize(ptr);
    benchmark::ClobberMemory();
  }
  state.counters["FieldBytes"] = static_cast<double>(field_bytes);
}

}  // namespace

// 1. GFArray bulk kernels
BENCHMARK(BM_GFArray_GF256_U8_AddArray)->Arg(256)->Arg(4096)->Arg(65536);
BENCHMARK(BM_GFArray_GF256_U8_MulArray)->Arg(256)->Arg(4096)->Arg(65536);
BENCHMARK(BM_GFArray_GF256_U8_FmaArray)->Arg(256)->Arg(4096)->Arg(65536);

BENCHMARK(BM_GFArray_Prime65521_U32_AddArray)->Arg(256)->Arg(4096)->Arg(65536);
BENCHMARK(BM_GFArray_Prime65521_U32_MulArray)->Arg(256)->Arg(4096)->Arg(65536);
BENCHMARK(BM_GFArray_Prime65521_U32_FmaArray)->Arg(256)->Arg(4096)->Arg(65536);

// 2. GaloisFieldElement vs raw field op overhead
BENCHMARK(BM_Overhead_GF256_RawMul)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_Overhead_GF256_ElementMul)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_Overhead_GF256_ElementMulAssign)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_Overhead_Prime65521_RawMul)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_Overhead_Prime65521_ElementMul)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_Overhead_Prime65521_ElementMulAssign)
    ->Unit(benchmark::kNanosecond);

// 3. PolynomialDense multiplication (schoolbook vs Karatsuba) & DivRem
BENCHMARK(BM_PolyDense_GF256_Mul)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_PolyDense_Prime65521_Mul)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_PolyDense_GF256_DivRem_256_128)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_PolyDense_Prime65521_DivRem_256_128)
    ->Unit(benchmark::kMicrosecond);

// 4. Linear algebra: GF(2) bit-packed M4RI and GF(2^8) GFArray
BENCHMARK(BM_Linalg_GF2_M4RI_Dot)
    ->Arg(64)
    ->Arg(256)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Linalg_GF2_M4RI_Rref)
    ->Arg(64)
    ->Arg(256)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Linalg_GF256_Dot)
    ->Arg(64)
    ->Arg(128)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Linalg_GF256_Rref)
    ->Arg(64)
    ->Arg(128)
    ->Unit(benchmark::kMicrosecond);

// 5. RS(255, 223) throughput in MB/s
BENCHMARK(BM_RS255_223_Encode)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_RS255_223_BerlekampMassey_16Errors)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_RS255_223_Sugiyama_16Errors)->Unit(benchmark::kMicrosecond);

// 6. Field construction latency
BENCHMARK(BM_FieldConstruct_GF2X)
    ->DenseRange(8, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_FieldConstruct_GF2XLOG)
    ->DenseRange(8, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_FieldConstruct_GF2XZECH)
    ->DenseRange(8, 16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_FieldConstruct_GFP_65521)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_FieldConstruct_GFPTABLE_65521)->Unit(benchmark::kMicrosecond);
