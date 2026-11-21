#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
#include "xgalois/coding/decoder/sugiyama_decoder.hpp"
#include "xgalois/coding/encoder/grs_encoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/ldpc.hpp"
#include "xgalois/coding/polar.hpp"
#include "xgalois/field/gf128.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"

using xg::GaloisFieldElement;
using xg::GaloisFieldPrime;
using xg::GaloisFieldPrimeTable;
using xg::GF128;
using xg::GF2;
using xg::GF2X;
using xg::GF2XLOG;
using xg::GF2XZECH;
using xg::GFArray;
using xg::PolynomialDense;
using xg::UInt128;

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
// 4. Linear Algebra: linalg::Dot & linalg::Rref on GF(2) M4RI and GF(2^8)
// ============================================================================

void BM_Linalg_GF2_M4RI_Dot(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf2 = std::make_shared<GF2>();
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::RandomMatrix(gf2, {n, n}, rng);
  const auto b = xg::linalg::RandomMatrix(gf2, {n, n}, rng);

  for (auto _ : state) {
    auto c = xg::linalg::Dot(a, b);
    benchmark::DoNotOptimize(c.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

void BM_Linalg_GF2_M4RI_Rref(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf2 = std::make_shared<GF2>();
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::RandomMatrix(gf2, {n, n}, rng);

  for (auto _ : state) {
    auto r = xg::linalg::Rref(a);
    benchmark::DoNotOptimize(r.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

void BM_Linalg_GF256_Dot(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::RandomMatrix(gf256, {n, n}, rng);
  const auto b = xg::linalg::RandomMatrix(gf256, {n, n}, rng);

  for (auto _ : state) {
    auto c = xg::linalg::Dot(a, b);
    benchmark::DoNotOptimize(c.RawSpan().data());
    benchmark::ClobberMemory();
  }
}

void BM_Linalg_GF256_Rref(benchmark::State& state) {
  const size_t n = static_cast<size_t>(state.range(0));
  auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
  std::mt19937_64 rng(kSeed);

  const auto a = xg::linalg::RandomMatrix(gf256, {n, n}, rng);

  for (auto _ : state) {
    auto r = xg::linalg::Rref(a);
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

// ============================================================================
// 7. GF(2^128) & GHASH (PCLMULQDQ / Portable Carry-Less)
// ============================================================================

void BM_GF128_Mul(benchmark::State& state) {
  GF128 field;
  std::mt19937_64 rng(kSeed);
  std::vector<UInt128> vals(kScalarBufSize);
  for (auto& v : vals) {
    v = field.Random(rng);
  }
  size_t idx = 0;

  for (auto _ : state) {
    const UInt128 res =
        field.Mul(vals[idx & kScalarBufMask], vals[(idx + 1) & kScalarBufMask]);
    benchmark::DoNotOptimize(res);
    ++idx;
  }
}

void BM_GF128_GhashMul(benchmark::State& state) {
  GF128 field;
  std::mt19937_64 rng(kSeed);
  std::vector<UInt128> vals(kScalarBufSize);
  for (auto& v : vals) {
    v = field.Random(rng);
  }
  size_t idx = 0;

  for (auto _ : state) {
    const UInt128 res = field.GhashMul(vals[idx & kScalarBufMask],
                                       vals[(idx + 1) & kScalarBufMask]);
    benchmark::DoNotOptimize(res);
    ++idx;
  }
}

void BM_GF128_Inv(benchmark::State& state) {
  GF128 field;
  std::mt19937_64 rng(kSeed);
  std::vector<UInt128> vals(kScalarBufSize);
  for (auto& v : vals) {
    v = field.Random(rng);
    while (field.IsZero(v)) {
      v = field.Random(rng);
    }
  }
  size_t idx = 0;

  for (auto _ : state) {
    const UInt128 res = field.Inv(vals[idx & kScalarBufMask]);
    benchmark::DoNotOptimize(res);
    ++idx;
  }
}

void BM_GFArray_GF128_MulArray(benchmark::State& state) {
  BenchBulkMulArray<GF128>(state, [] { return std::make_shared<GF128>(); });
}

// ============================================================================
// 8. QC-LDPC Belief Propagation (Normalized Min-Sum over AWGN LLRs)
// ============================================================================

struct LdpcQcFixture {
  static constexpr size_t kLiftingZ = 17;
  static constexpr size_t kBatchSize = 16;
  static constexpr size_t kBatchMask = kBatchSize - 1;
  static constexpr size_t kMaxIterations = 50;
  static constexpr double kAlpha = 0.75;
  static constexpr double kSnrDb = 4.5;

  std::shared_ptr<GF2> field;
  xg::coding::LDPCCode code;
  std::vector<GFArray<GF2>> codewords;
  std::vector<std::vector<double>> llrs;

  static std::vector<std::vector<int>> BaseMatrix() {
    return {
        {0, 0, 0, 0, 0, 0},
        {0, 1, 2, 3, 4, 5},
        {0, 2, 4, 6, 8, 10},
    };
  }

  LdpcQcFixture()
      : field(std::make_shared<GF2>()),
        code(xg::coding::LDPCCode::FromQuasiCyclic(field, BaseMatrix(),
                                                   kLiftingZ)) {
    std::mt19937_64 rng(kSeed);
    auto channel = xg::channel::AWGNChannel::FromEbN0(
        field, code.Length(), kSnrDb, code.Rate(), kSeed);
    xg::coding::BeliefPropagationDecoder decoder(
        code, xg::coding::BpAlgorithm::kNormalizedMinSum, kMaxIterations,
        kAlpha);

    codewords.reserve(kBatchSize);
    llrs.reserve(kBatchSize);

    for (size_t attempt = 0; attempt < 512 && llrs.size() < kBatchSize;
         ++attempt) {
      auto msg = GFArray<GF2>::Random(field, {code.Dimension()}, rng);
      auto cw = code.Encode(msg);
      auto tx = channel.TransmitDetailed(cw);
      if (code.CheckParity(tx.hard_decision)) {
        continue;
      }
      const auto res = decoder.DecodeLlr(tx.llr);
      if (res.Ok() && res.codeword == cw) {
        codewords.push_back(std::move(cw));
        llrs.push_back(std::move(tx.llr));
      }
    }
    while (llrs.size() < kBatchSize) {
      auto msg = GFArray<GF2>::Random(field, {code.Dimension()}, rng);
      auto cw = code.Encode(msg);
      auto tx_llr = channel.TransmitLlr(cw);
      codewords.push_back(std::move(cw));
      llrs.push_back(std::move(tx_llr));
    }
  }
};

const LdpcQcFixture& GetLdpcQcFixture() {
  static const LdpcQcFixture fixture;
  return fixture;
}

void BM_LDPC_NormalizedMinSum_Decode(benchmark::State& state) {
  const auto& fix = GetLdpcQcFixture();
  xg::coding::BeliefPropagationDecoder decoder(
      fix.code, xg::coding::BpAlgorithm::kNormalizedMinSum,
      LdpcQcFixture::kMaxIterations, LdpcQcFixture::kAlpha);

  const auto check = decoder.DecodeLlr(fix.llrs[0]);
  if (!check.Ok() || check.codeword != fix.codewords[0]) {
    state.SkipWithError(
        "BeliefPropagationDecoder NormalizedMinSum self-check failed");
    return;
  }

  size_t idx = 0;
  for (auto _ : state) {
    auto res = decoder.DecodeLlr(fix.llrs[idx & LdpcQcFixture::kBatchMask]);
    benchmark::DoNotOptimize(res.codeword.RawSpan().data());
    benchmark::DoNotOptimize(res.num_errors);
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

// ============================================================================
// 9. Polar(64, 32): Encode, SC Decode, SCL(L=4) Decode over AWGN LLRs
// ============================================================================

struct Polar64Fixture {
  static constexpr size_t kCodeLen = 64;
  static constexpr size_t kMessageLen = 32;
  static constexpr size_t kBatchSize = 16;
  static constexpr size_t kBatchMask = kBatchSize - 1;
  static constexpr double kDesignSnrDb = 2.5;
  static constexpr double kChannelSnrDb = 4.0;

  std::shared_ptr<GF2> field;
  xg::coding::PolarCode code;
  std::vector<GFArray<GF2>> messages;
  std::vector<GFArray<GF2>> codewords;
  std::vector<std::vector<double>> llrs;

  Polar64Fixture()
      : field(std::make_shared<GF2>()),
        code(xg::coding::PolarCode::FromDesignSnr(
            field, kCodeLen, kMessageLen, kDesignSnrDb, /*systematic=*/false)) {
    std::mt19937_64 rng(kSeed);
    xg::coding::PolarEncoder encoder(code);
    xg::coding::SuccessiveCancellationDecoder sc(code);
    xg::coding::SuccessiveCancellationListDecoder scl4(code, /*list_size=*/4);
    auto channel = xg::channel::AWGNChannel::FromEbN0(
        field, kCodeLen, kChannelSnrDb, code.Rate(), kSeed);

    messages.reserve(kBatchSize);
    codewords.reserve(kBatchSize);
    llrs.reserve(kBatchSize);

    for (size_t attempt = 0; attempt < 512 && llrs.size() < kBatchSize;
         ++attempt) {
      auto msg = GFArray<GF2>::Random(field, {kMessageLen}, rng);
      auto cw = encoder.Encode(msg);
      auto rx_llr = channel.TransmitLlr(cw);
      const auto r_sc = sc.DecodeLlr(rx_llr);
      const auto r_scl = scl4.DecodeLlr(rx_llr);
      if (r_sc.Ok() && r_sc.codeword == cw && r_scl.Ok() &&
          r_scl.codeword == cw) {
        messages.push_back(std::move(msg));
        codewords.push_back(std::move(cw));
        llrs.push_back(std::move(rx_llr));
      }
    }
    while (llrs.size() < kBatchSize) {
      auto msg = GFArray<GF2>::Random(field, {kMessageLen}, rng);
      auto cw = encoder.Encode(msg);
      auto rx_llr = channel.TransmitLlr(cw);
      messages.push_back(std::move(msg));
      codewords.push_back(std::move(cw));
      llrs.push_back(std::move(rx_llr));
    }
  }
};

const Polar64Fixture& GetPolar64Fixture() {
  static const Polar64Fixture fixture;
  return fixture;
}

void BM_Polar64_32_Encode(benchmark::State& state) {
  const auto& fix = GetPolar64Fixture();
  xg::coding::PolarEncoder encoder(fix.code);
  size_t idx = 0;

  for (auto _ : state) {
    auto cw = encoder.Encode(fix.messages[idx & Polar64Fixture::kBatchMask]);
    benchmark::DoNotOptimize(cw.RawSpan().data());
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

void BM_Polar64_32_SC_Decode(benchmark::State& state) {
  const auto& fix = GetPolar64Fixture();
  xg::coding::SuccessiveCancellationDecoder decoder(fix.code);

  const auto check = decoder.DecodeLlr(fix.llrs[0]);
  if (!check.Ok() || check.codeword != fix.codewords[0] ||
      check.message != fix.messages[0]) {
    state.SkipWithError("Polar SC decoder self-check failed");
    return;
  }

  size_t idx = 0;
  for (auto _ : state) {
    auto res = decoder.DecodeLlr(fix.llrs[idx & Polar64Fixture::kBatchMask]);
    benchmark::DoNotOptimize(res.codeword.RawSpan().data());
    benchmark::DoNotOptimize(res.num_errors);
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

void BM_Polar64_32_SCL4_Decode(benchmark::State& state) {
  const auto& fix = GetPolar64Fixture();
  xg::coding::SuccessiveCancellationListDecoder decoder(fix.code,
                                                        /*list_size=*/4);

  const auto check = decoder.DecodeLlr(fix.llrs[0]);
  if (!check.Ok() || check.codeword != fix.codewords[0] ||
      check.message != fix.messages[0]) {
    state.SkipWithError("Polar SCL(L=4) decoder self-check failed");
    return;
  }

  size_t idx = 0;
  for (auto _ : state) {
    auto res = decoder.DecodeLlr(fix.llrs[idx & Polar64Fixture::kBatchMask]);
    benchmark::DoNotOptimize(res.codeword.RawSpan().data());
    benchmark::DoNotOptimize(res.num_errors);
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
}

// ============================================================================
// 10. Guruswami-Sudan List Decoding: RS(15, 3) over GF(2^4) with 8 Errors
//     (unique radius t = 6, Johnson radius tau_J = 9)
// ============================================================================

struct GsRs15Fixture {
  using Field = GF2XLOG<uint8_t>;
  using Code = xg::coding::GeneralizedReedSolomonCode<Field>;

  static constexpr size_t kCodeLen = 15;
  static constexpr size_t kMessageLen = 3;
  static constexpr size_t kNumErrors = 8;
  static constexpr size_t kBatchSize = 16;
  static constexpr size_t kBatchMask = kBatchSize - 1;

  std::shared_ptr<Field> field;
  Code rs;
  std::vector<GFArray<Field>> codewords;
  std::vector<GFArray<Field>> received_8err;

  GsRs15Fixture()
      : field(std::make_shared<Field>(4)),
        rs(Code::ReedSolomon(field, kCodeLen, kMessageLen)) {
    std::mt19937_64 rng(kSeed);
    xg::coding::GRSEncoder<Field> encoder(rs);

    std::vector<size_t> positions(kCodeLen);
    std::iota(positions.begin(), positions.end(), size_t{0});

    codewords.reserve(kBatchSize);
    received_8err.reserve(kBatchSize);

    for (size_t b = 0; b < kBatchSize; ++b) {
      auto msg = GFArray<Field>::Random(field, {kMessageLen}, rng);
      auto cw = encoder.Encode(msg);
      auto rx = cw;

      std::shuffle(positions.begin(), positions.end(), rng);
      for (size_t e = 0; e < kNumErrors; ++e) {
        const size_t pos = positions[e];
        uint8_t err = field->Random(rng);
        while (field->IsZero(err)) {
          err = field->Random(rng);
        }
        rx.Raw(pos) = field->Add(rx.Raw(pos), err);
      }

      codewords.push_back(std::move(cw));
      received_8err.push_back(std::move(rx));
    }
  }
};

const GsRs15Fixture& GetGsRs15Fixture() {
  static const GsRs15Fixture fixture;
  return fixture;
}

void BM_GuruswamiSudan_RS15_3_8Errors(benchmark::State& state) {
  const auto& fix = GetGsRs15Fixture();
  const auto decoder =
      xg::coding::GuruswamiSudanDecoder<GsRs15Fixture::Field>::FromTargetRadius(
          fix.rs, GsRs15Fixture::kNumErrors);

  const auto check_list = decoder.DecodeToList(fix.received_8err[0]);
  if (std::find(check_list.begin(), check_list.end(), fix.codewords[0]) ==
      check_list.end()) {
    state.SkipWithError("GuruswamiSudanDecoder 8-error self-check failed");
    return;
  }

  size_t idx = 0;
  for (auto _ : state) {
    auto list = decoder.DecodeToList(
        fix.received_8err[idx & GsRs15Fixture::kBatchMask]);
    const size_t list_len = list.size();
    benchmark::DoNotOptimize(list.data());
    benchmark::DoNotOptimize(list_len);
    benchmark::ClobberMemory();
    ++idx;
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(GsRs15Fixture::kCodeLen));
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

// 7. GF(2^128) & GHASH arithmetic (PCLMULQDQ / portable carry-less)
BENCHMARK(BM_GF128_Mul)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF128_GhashMul)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GF128_Inv)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_GFArray_GF128_MulArray)->Arg(1024)->Unit(benchmark::kMicrosecond);

// 8. QC-LDPC Belief Propagation (Normalized Min-Sum)
BENCHMARK(BM_LDPC_NormalizedMinSum_Decode)->Unit(benchmark::kMicrosecond);

// 9. Polar(64, 32) encode and SC / SCL(L=4) decoding
BENCHMARK(BM_Polar64_32_Encode)->Unit(benchmark::kNanosecond);
BENCHMARK(BM_Polar64_32_SC_Decode)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Polar64_32_SCL4_Decode)->Unit(benchmark::kMicrosecond);

// 10. Guruswami-Sudan RS(15, 3) list decoding beyond unique radius (8 errors)
BENCHMARK(BM_GuruswamiSudan_RS15_3_8Errors)->Unit(benchmark::kMicrosecond);
