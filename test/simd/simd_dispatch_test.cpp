#include <gtest/gtest.h>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "xgalois/simd/dispatch.hpp"
#include "xgalois/simd/kernels_scalar.hpp"
#include "xgalois/simd/thread_pool.hpp"

namespace {

using xg::simd::ComputeShoupFactor;
using xg::simd::DetectIsaLevel;
using xg::simd::GetActiveIsaLevel;
using xg::simd::GetForcedIsaLevel;
using xg::simd::Gf128ClmulScalar;
using xg::simd::Gf128MulArray;
using xg::simd::Gf128MulScalar;
using xg::simd::Gf28AddArray;
using xg::simd::Gf28FmaArray;
using xg::simd::Gf28MulArray;
using xg::simd::Gf28MulScalar;
using xg::simd::Gf28ScalarMulArray;
using xg::simd::Gf2XorArray;
using xg::simd::GfPrimeAddArray;
using xg::simd::GfPrimeMulShoupArray;
using xg::simd::GfPrimeSubArray;
using xg::simd::GlobalThreadPool;
using xg::simd::IsaLevel;
using xg::simd::IsaLevelName;
using xg::simd::ScopedIsaOverride;
using xg::simd::SetForcedIsaLevel;
using xg::simd::ShoupMulMod;
using xg::simd::ThreadPool;

// Slow bit-by-bit reference GF(2^128) multiplication modulo x^128 + R(x).
std::pair<uint64_t, uint64_t> ReferenceGf128Mul(std::pair<uint64_t, uint64_t> a,
                                                std::pair<uint64_t, uint64_t> b,
                                                uint64_t reduction_poly_lo) {
  uint64_t res_lo = 0;
  uint64_t res_hi = 0;
  uint64_t cur_lo = a.first;
  uint64_t cur_hi = a.second;

  for (int bit = 0; bit < 128; ++bit) {
    const uint64_t b_word = (bit < 64) ? b.first : b.second;
    const unsigned b_shift = static_cast<unsigned>(bit & 63);
    if (((b_word >> b_shift) & 1ULL) != 0ULL) {
      res_lo ^= cur_lo;
      res_hi ^= cur_hi;
    }

    const bool top_bit = (cur_hi >> 63) != 0ULL;
    cur_hi = (cur_hi << 1) | (cur_lo >> 63);
    cur_lo = cur_lo << 1;
    if (top_bit) {
      cur_lo ^= reduction_poly_lo;
    }
  }
  return {res_lo, res_hi};
}

// Reference modular multiplication using 128-bit arithmetic when available,
// or Russian peasant modular doubling as a portable fallback.
uint64_t ReferenceMulMod(uint64_t a, uint64_t b, uint64_t mod) {
#if defined(__SIZEOF_INT128__)
  __extension__ using Uint128 = unsigned __int128;
  return static_cast<uint64_t>((static_cast<Uint128>(a) * b) % mod);
#else
  uint64_t res = 0;
  uint64_t lhs = a % mod;
  uint64_t rhs = b;
  while (rhs > 0) {
    if ((rhs & 1ULL) != 0ULL) {
      const uint64_t diff = mod - lhs;
      res = (res >= diff) ? (res - diff) : (res + lhs);
    }
    const uint64_t diff = mod - lhs;
    lhs = (lhs >= diff) ? (lhs - diff) : (lhs + lhs);
    rhs >>= 1;
  }
  return res;
#endif
}

// ============================================================================
// ISA Dispatch Tests
// ============================================================================

TEST(SimdDispatchTest, DetectIsaLevelIsConsistent) {
  const IsaLevel first = DetectIsaLevel();
  const IsaLevel second = DetectIsaLevel();
  EXPECT_EQ(first, second);
  EXPECT_FALSE(IsaLevelName(first).empty());
  EXPECT_NE(IsaLevelName(first), "Unknown");
}

TEST(SimdDispatchTest, IsaLevelNamesAreValidAndDistinct) {
  const std::vector<IsaLevel> levels = {
      IsaLevel::kScalar,     IsaLevel::kSse42Pclmul, IsaLevel::kAvx2,
      IsaLevel::kAvx512Gfni, IsaLevel::kNeon,        IsaLevel::kNeonPmull,
  };
  for (std::size_t i = 0; i < levels.size(); ++i) {
    const std::string_view name = IsaLevelName(levels[i]);
    EXPECT_FALSE(name.empty());
    EXPECT_NE(name, "Unknown");
    for (std::size_t j = i + 1; j < levels.size(); ++j) {
      EXPECT_NE(name, IsaLevelName(levels[j]));
    }
  }
  EXPECT_EQ(IsaLevelName(static_cast<IsaLevel>(255)), "Unknown");
}

TEST(SimdDispatchTest, ScopedIsaOverrideRestoresState) {
  SetForcedIsaLevel(std::nullopt);
  const IsaLevel detected = DetectIsaLevel();
  EXPECT_EQ(GetForcedIsaLevel(), std::nullopt);
  EXPECT_EQ(GetActiveIsaLevel(), detected);

  {
    ScopedIsaOverride outer(IsaLevel::kScalar);
    EXPECT_EQ(GetForcedIsaLevel(), std::optional<IsaLevel>(IsaLevel::kScalar));
    EXPECT_EQ(GetActiveIsaLevel(), IsaLevel::kScalar);

    {
      ScopedIsaOverride inner(IsaLevel::kAvx2);
      EXPECT_EQ(GetForcedIsaLevel(), std::optional<IsaLevel>(IsaLevel::kAvx2));
      EXPECT_EQ(GetActiveIsaLevel(), IsaLevel::kAvx2);

      {
        ScopedIsaOverride auto_mode(std::nullopt);
        EXPECT_EQ(GetForcedIsaLevel(), std::nullopt);
        EXPECT_EQ(GetActiveIsaLevel(), detected);
      }

      EXPECT_EQ(GetActiveIsaLevel(), IsaLevel::kAvx2);
    }

    EXPECT_EQ(GetActiveIsaLevel(), IsaLevel::kScalar);
  }

  EXPECT_EQ(GetForcedIsaLevel(), std::nullopt);
  EXPECT_EQ(GetActiveIsaLevel(), detected);
}

// ============================================================================
// GF(2) and GF(2^8) Scalar Kernel Tests
// ============================================================================

TEST(SimdKernelsScalarTest, Gf2XorArrayBasicAndInPlace) {
  std::vector<uint64_t> a = {0ULL, 0x0123456789ABCDEFULL, ~0ULL,
                             0xDEADBEEFCAFEBABEULL};
  std::vector<uint64_t> b = {~0ULL, 0xFEDCBA9876543210ULL,
                             0x5555555555555555ULL, 0xDEADBEEFCAFEBABEULL};
  std::vector<uint64_t> dst(a.size(), 0);

  Gf2XorArray(dst, a, b);
  for (std::size_t i = 0; i < a.size(); ++i) {
    EXPECT_EQ(dst[i], a[i] ^ b[i]);
  }

  // In-place aliasing dst == a.
  std::vector<uint64_t> in_place = a;
  Gf2XorArray(in_place, in_place, b);
  EXPECT_EQ(in_place, dst);

  // Empty spans should succeed.
  std::vector<uint64_t> empty;
  EXPECT_NO_THROW(Gf2XorArray(empty, empty, empty));

  // Size mismatch throws.
  std::vector<uint64_t> short_vec(2, 0);
  EXPECT_THROW(Gf2XorArray(short_vec, a, b), std::invalid_argument);
}

TEST(SimdKernelsScalarTest, Gf28AddArray) {
  std::vector<uint8_t> a(256);
  std::vector<uint8_t> b(256);
  std::vector<uint8_t> dst(256);
  for (std::size_t i = 0; i < 256; ++i) {
    a[i] = static_cast<uint8_t>(i);
    b[i] = static_cast<uint8_t>(255 - i);
  }
  Gf28AddArray(dst, a, b);
  for (std::size_t i = 0; i < 256; ++i) {
    EXPECT_EQ(dst[i], static_cast<uint8_t>(a[i] ^ b[i]));
  }

  std::vector<uint8_t> bad(10);
  EXPECT_THROW(Gf28AddArray(bad, a, b), std::invalid_argument);
}

TEST(SimdKernelsScalarTest, Gf28MulScalarAesAndReedSolomon) {
  constexpr uint16_t kAesPoly = 0x11B;
  constexpr uint16_t kRsPoly = 0x11D;

  // Known FIPS-197 AES GF(2^8) vectors.
  EXPECT_EQ(Gf28MulScalar(0x57, 0x83, kAesPoly), 0xC1);
  EXPECT_EQ(Gf28MulScalar(0x57, 0x13, kAesPoly), 0xFE);
  EXPECT_EQ(Gf28MulScalar(0x53, 0xCA, kAesPoly), 0x01);

  // Exhaustive field-axiom check on both 0x11D and 0x11B.
  for (uint16_t poly : {kRsPoly, kAesPoly}) {
    for (unsigned a = 0; a < 256; ++a) {
      EXPECT_EQ(Gf28MulScalar(static_cast<uint8_t>(a), 0, poly), 0);
      EXPECT_EQ(Gf28MulScalar(0, static_cast<uint8_t>(a), poly), 0);
      EXPECT_EQ(Gf28MulScalar(static_cast<uint8_t>(a), 1, poly), a);
      EXPECT_EQ(Gf28MulScalar(1, static_cast<uint8_t>(a), poly), a);

      if (a != 0) {
        int inverse_count = 0;
        for (unsigned b = 1; b < 256; ++b) {
          const uint8_t ab = Gf28MulScalar(static_cast<uint8_t>(a),
                                           static_cast<uint8_t>(b), poly);
          const uint8_t ba = Gf28MulScalar(static_cast<uint8_t>(b),
                                           static_cast<uint8_t>(a), poly);
          EXPECT_EQ(ab, ba);
          EXPECT_NE(ab, 0);
          if (ab == 1) {
            ++inverse_count;
          }
        }
        EXPECT_EQ(inverse_count, 1);
      }
    }
  }
}

TEST(SimdKernelsScalarTest, Gf28ArrayOpsMatchScalarOn0x11DAnd0x11B) {
  std::vector<uint8_t> src(256);
  std::iota(src.begin(), src.end(), uint8_t{0});

  for (uint16_t poly : {uint16_t{0x11D}, uint16_t{0x11B}}) {
    for (uint8_t scalar : {uint8_t{0x00}, uint8_t{0x01}, uint8_t{0x02},
                           uint8_t{0x57}, uint8_t{0xC3}, uint8_t{0xFF}}) {
      std::vector<uint8_t> dst(256, 0xAA);
      Gf28ScalarMulArray(dst, src, scalar, poly);
      for (std::size_t i = 0; i < 256; ++i) {
        EXPECT_EQ(dst[i], Gf28MulScalar(src[i], scalar, poly));
      }

      // In-place Gf28ScalarMulArray.
      std::vector<uint8_t> in_place = src;
      Gf28ScalarMulArray(in_place, in_place, scalar, poly);
      EXPECT_EQ(in_place, dst);

      // Gf28FmaArray: dst_new[i] == init[i] ^ (src[i] * scalar).
      std::vector<uint8_t> acc(256);
      for (std::size_t i = 0; i < 256; ++i) {
        acc[i] = static_cast<uint8_t>((i * 37u + 13u) & 0xFFu);
      }
      const std::vector<uint8_t> init = acc;
      Gf28FmaArray(acc, src, scalar, poly);
      for (std::size_t i = 0; i < 256; ++i) {
        EXPECT_EQ(acc[i], static_cast<uint8_t>(
                              init[i] ^ Gf28MulScalar(src[i], scalar, poly)));
      }
    }

    // Element-wise Gf28MulArray.
    std::vector<uint8_t> rhs(256);
    for (std::size_t i = 0; i < 256; ++i) {
      rhs[i] = static_cast<uint8_t>((i * 53u + 7u) & 0xFFu);
    }
    std::vector<uint8_t> prod(256);
    Gf28MulArray(prod, src, rhs, poly);
    for (std::size_t i = 0; i < 256; ++i) {
      EXPECT_EQ(prod[i], Gf28MulScalar(src[i], rhs[i], poly));
    }
  }

  std::vector<uint8_t> bad(16);
  EXPECT_THROW(Gf28ScalarMulArray(bad, src, 0x05, 0x11D),
               std::invalid_argument);
  EXPECT_THROW(Gf28FmaArray(bad, src, 0x05, 0x11D), std::invalid_argument);
  EXPECT_THROW(Gf28MulArray(bad, src, src, 0x11D), std::invalid_argument);
}

// ============================================================================
// GF(2^128) Carry-Less Multiplication & Reduction Tests
// ============================================================================

TEST(SimdKernelsScalarTest, Gf128ClmulScalarEdgeCasesAndLinearity) {
  EXPECT_EQ(Gf128ClmulScalar(0ULL, 0ULL),
            (std::pair<uint64_t, uint64_t>{0, 0}));
  EXPECT_EQ(Gf128ClmulScalar(0ULL, ~0ULL),
            (std::pair<uint64_t, uint64_t>{0, 0}));
  EXPECT_EQ(Gf128ClmulScalar(1ULL, 0x123456789ABCDEF0ULL),
            (std::pair<uint64_t, uint64_t>{0x123456789ABCDEF0ULL, 0ULL}));
  // x^63 * x = x^64 -> {0, 1}
  EXPECT_EQ(Gf128ClmulScalar(1ULL << 63, 2ULL),
            (std::pair<uint64_t, uint64_t>{0ULL, 1ULL}));
  // x^63 * x^63 = x^126 -> {0, 1ULL << 62}
  EXPECT_EQ(Gf128ClmulScalar(1ULL << 63, 1ULL << 63),
            (std::pair<uint64_t, uint64_t>{0ULL, 1ULL << 62}));

  std::mt19937_64 rng(0x123456789ABCDEF0ULL);
  for (int i = 0; i < 200; ++i) {
    const uint64_t a = rng();
    const uint64_t b = rng();
    const uint64_t c = rng();

    EXPECT_EQ(Gf128ClmulScalar(a, b), Gf128ClmulScalar(b, a));

    const auto ab = Gf128ClmulScalar(a ^ b, c);
    const auto ac = Gf128ClmulScalar(a, c);
    const auto bc = Gf128ClmulScalar(b, c);
    EXPECT_EQ(ab.first, ac.first ^ bc.first);
    EXPECT_EQ(ab.second, ac.second ^ bc.second);
  }
}

TEST(SimdKernelsScalarTest, Gf128MulArrayGcm0x87) {
  constexpr uint64_t kGcmPolyLo = 0x87ULL;  // x^128 + x^7 + x^2 + x + 1

  // Known exact reductions.
  // x^64 * x^64 = x^128 == 0x87
  EXPECT_EQ(Gf128MulScalar({0ULL, 1ULL}, {0ULL, 1ULL}, kGcmPolyLo),
            (std::pair<uint64_t, uint64_t>{kGcmPolyLo, 0ULL}));
  // x^127 * x = x^128 == 0x87
  EXPECT_EQ(Gf128MulScalar({0ULL, 1ULL << 63}, {2ULL, 0ULL}, kGcmPolyLo),
            (std::pair<uint64_t, uint64_t>{kGcmPolyLo, 0ULL}));
  // x^127 * x^127 = x^127 + x^126 + x^12 + x^6 + x^5 + x^2 + x + 1
  EXPECT_EQ(Gf128MulScalar({0ULL, 1ULL << 63}, {0ULL, 1ULL << 63}, kGcmPolyLo),
            (std::pair<uint64_t, uint64_t>{0x1067ULL, 0xC000000000000000ULL}));

  std::mt19937_64 rng(0xCAFEBABE12345678ULL);
  constexpr std::size_t kCount = 128;
  std::vector<std::pair<uint64_t, uint64_t>> a(kCount);
  std::vector<std::pair<uint64_t, uint64_t>> b(kCount);
  std::vector<std::pair<uint64_t, uint64_t>> dst(kCount);

  a[0] = {0ULL, 0ULL};
  b[0] = {~0ULL, ~0ULL};
  a[1] = {1ULL, 0ULL};
  b[1] = {0x0123456789ABCDEFULL, 0xFEDCBA9876543210ULL};
  a[2] = {~0ULL, ~0ULL};
  b[2] = {~0ULL, ~0ULL};
  for (std::size_t i = 3; i < kCount; ++i) {
    a[i] = {rng(), rng()};
    b[i] = {rng(), rng()};
  }

  Gf128MulArray(dst, a, b, kGcmPolyLo);
  for (std::size_t i = 0; i < kCount; ++i) {
    EXPECT_EQ(dst[i], ReferenceGf128Mul(a[i], b[i], kGcmPolyLo));
    EXPECT_EQ(dst[i], Gf128MulScalar(b[i], a[i], kGcmPolyLo));
  }

  // Also verify reduction when bit 63 of reduction_poly_lo is set.
  constexpr uint64_t kHighBitPolyLo = (1ULL << 63) | 0x87ULL;
  Gf128MulArray(dst, a, b, kHighBitPolyLo);
  for (std::size_t i = 0; i < kCount; ++i) {
    EXPECT_EQ(dst[i], ReferenceGf128Mul(a[i], b[i], kHighBitPolyLo));
  }

  std::vector<std::pair<uint64_t, uint64_t>> bad(5);
  EXPECT_THROW(Gf128MulArray(bad, a, b, kGcmPolyLo), std::invalid_argument);
}

// ============================================================================
// Prime Field & Shoup Modular Multiplication Tests
// ============================================================================

TEST(SimdKernelsScalarTest, GfPrimeAddSubArray) {
  const std::vector<uint64_t> moduli = {
      2013265921ULL,           // 31-bit NTT prime (15 * 2^27 + 1)
      1945555039024054273ULL,  // 61-bit NTT prime (27 * 2^56 + 1)
      UINT64_MAX - 58ULL,  // Largest 64-bit prime (2^64 - 59), tests overflow
  };

  std::mt19937_64 rng(0x9876543210FEDCBAULL);
  constexpr std::size_t kCount = 256;

  for (uint64_t mod : moduli) {
    std::vector<uint64_t> a(kCount);
    std::vector<uint64_t> b(kCount);
    a[0] = 0;
    b[0] = 0;
    a[1] = mod - 1;
    b[1] = mod - 1;
    a[2] = 0;
    b[2] = mod - 1;
    a[3] = mod - 1;
    b[3] = 1;
    for (std::size_t i = 4; i < kCount; ++i) {
      a[i] = rng() % mod;
      b[i] = rng() % mod;
    }

    std::vector<uint64_t> sum(kCount);
    std::vector<uint64_t> diff(kCount);
    std::vector<uint64_t> roundtrip(kCount);

    GfPrimeAddArray(sum, a, b, mod);
    GfPrimeSubArray(diff, a, b, mod);

    for (std::size_t i = 0; i < kCount; ++i) {
      EXPECT_LT(sum[i], mod);
      EXPECT_LT(diff[i], mod);
    }

    // (a + b) - b == a and (a - b) + b == a.
    GfPrimeSubArray(roundtrip, sum, b, mod);
    EXPECT_EQ(roundtrip, a);
    GfPrimeAddArray(roundtrip, diff, b, mod);
    EXPECT_EQ(roundtrip, a);
  }

  std::vector<uint64_t> v = {1, 2};
  EXPECT_THROW(GfPrimeAddArray(v, v, v, 0), std::invalid_argument);
  EXPECT_THROW(GfPrimeSubArray(v, v, v, 0), std::invalid_argument);
}

TEST(SimdKernelsScalarTest, ShoupMulModAndArrayOn31BitAnd61BitPrimes) {
  const std::vector<uint64_t> primes = {
      998244353ULL,            // 30-bit NTT prime (119 * 2^23 + 1)
      2013265921ULL,           // 31-bit NTT prime (15 * 2^27 + 1)
      1945555039024054273ULL,  // 61-bit NTT prime (27 * 2^56 + 1)
      2305843009213693951ULL,  // 61-bit Mersenne prime (2^61 - 1)
  };

  std::mt19937_64 rng(0x5A5A5A5A5A5A5A5AULL);
  constexpr std::size_t kCount = 512;

  for (uint64_t mod : primes) {
    const std::vector<uint64_t> multipliers = {
        0ULL, 1ULL, 2ULL, mod - 2ULL, mod - 1ULL, rng() % mod, rng() % mod,
    };

    std::vector<uint64_t> a(kCount);
    a[0] = 0ULL;
    a[1] = 1ULL;
    a[2] = mod - 1ULL;
    a[3] = mod - 2ULL;
    for (std::size_t i = 4; i < kCount; ++i) {
      a[i] = rng() % mod;
    }

    for (uint64_t b : multipliers) {
      const uint64_t b_shoup = ComputeShoupFactor(b, mod);
      std::vector<uint64_t> dst(kCount);
      GfPrimeMulShoupArray(dst, a, b, b_shoup, mod);

      for (std::size_t i = 0; i < kCount; ++i) {
        const uint64_t expected = ReferenceMulMod(a[i], b, mod);
        EXPECT_EQ(dst[i], expected)
            << "mod=" << mod << " a=" << a[i] << " b=" << b;
        EXPECT_EQ(ShoupMulMod(a[i], b, b_shoup, mod), expected);
      }

      // In-place aliasing dst == a.
      std::vector<uint64_t> in_place = a;
      GfPrimeMulShoupArray(in_place, in_place, b, b_shoup, mod);
      EXPECT_EQ(in_place, dst);
    }
  }

  // Invalid preconditions throw.
  EXPECT_THROW(ComputeShoupFactor(0, 0), std::invalid_argument);
  EXPECT_THROW(ComputeShoupFactor(1, 1ULL << 62), std::invalid_argument);
  EXPECT_THROW(ComputeShoupFactor(17, 17), std::invalid_argument);

  std::vector<uint64_t> v(4, 1);
  std::vector<uint64_t> bad(3, 1);
  EXPECT_THROW(GfPrimeMulShoupArray(bad, v, 1, 0, 17), std::invalid_argument);
  EXPECT_THROW(GfPrimeMulShoupArray(v, v, 1, 0, 0), std::invalid_argument);
  EXPECT_THROW(GfPrimeMulShoupArray(v, v, 1, 0, 1ULL << 62),
               std::invalid_argument);
}

// ============================================================================
// ThreadPool Tests
// ============================================================================

TEST(ThreadPoolTest, NumWorkersAndSubmitValuesAndExceptions) {
  ThreadPool pool(4);
  EXPECT_EQ(pool.NumWorkers(), 4u);

  ThreadPool min_pool(0);
  EXPECT_EQ(min_pool.NumWorkers(), 1u);

  auto f1 = pool.Submit([](int x, int y) { return x + y; }, 19, 23);
  EXPECT_EQ(f1.get(), 42);

  // Move-only argument.
  auto ptr = std::make_unique<int>(99);
  auto f2 = pool.Submit([](std::unique_ptr<int> p) { return *p * 2; },
                        std::move(ptr));
  EXPECT_EQ(f2.get(), 198);

  // Void task.
  std::atomic<int> flag{0};
  auto f3 =
      pool.Submit([&flag]() { flag.store(7, std::memory_order_release); });
  f3.get();
  EXPECT_EQ(flag.load(std::memory_order_acquire), 7);

  // Exception propagation through std::future.
  auto f4 =
      pool.Submit([]() -> int { throw std::runtime_error("submit failure"); });
  EXPECT_THROW(f4.get(), std::runtime_error);
}

TEST(ThreadPoolTest, ParallelForEdgeCasesAndLargeSum) {
  ThreadPool pool(4);

  // Empty ranges: begin == end and begin > end.
  std::atomic<std::size_t> calls{0};
  pool.ParallelFor(5, 5, [&](std::size_t) { calls.fetch_add(1); });
  pool.ParallelFor(10, 5, [&](std::size_t) { calls.fetch_add(1); });
  EXPECT_EQ(calls.load(), 0u);

  // Single element with both index and range callbacks.
  std::size_t single_idx = 0;
  pool.ParallelFor(7, 8, [&](std::size_t i) { single_idx = i; });
  EXPECT_EQ(single_idx, 7u);

  std::pair<std::size_t, std::size_t> single_range{0, 0};
  pool.ParallelFor(
      7, 8, [&](std::size_t lo, std::size_t hi) { single_range = {lo, hi}; });
  EXPECT_EQ(single_range, (std::pair<std::size_t, std::size_t>{7, 8}));

  // Large 100k sum using range callback.
  constexpr std::size_t kN = 100000;
  constexpr uint64_t kExpectedSum = static_cast<uint64_t>(kN) * (kN - 1) / 2ULL;

  std::atomic<uint64_t> range_sum{0};
  pool.ParallelFor(0, kN, [&](std::size_t lo, std::size_t hi) {
    uint64_t local = 0;
    for (std::size_t i = lo; i < hi; ++i) {
      local += i;
    }
    range_sum.fetch_add(local, std::memory_order_relaxed);
  });
  EXPECT_EQ(range_sum.load(), kExpectedSum);

  // Large 100k sum with per-index callback and explicit grain size.
  std::vector<uint64_t> squares(kN, 0);
  pool.ParallelFor(
      0, kN, [&](std::size_t i) { squares[i] = static_cast<uint64_t>(i); },
      1024);
  const uint64_t idx_sum =
      std::accumulate(squares.begin(), squares.end(), uint64_t{0});
  EXPECT_EQ(idx_sum, kExpectedSum);
}

TEST(ThreadPoolTest, NestedParallelForDoesNotDeadlock) {
  // Deliberately small pool (2 workers) with 3-level nesting so every worker is
  // saturated simultaneously.
  ThreadPool pool(2);
  constexpr std::size_t kOuter = 16;
  constexpr std::size_t kMiddle = 16;
  constexpr std::size_t kInner = 32;

  std::atomic<uint64_t> total{0};
  pool.ParallelFor(
      0, kOuter,
      [&](std::size_t i) {
        pool.ParallelFor(
            0, kMiddle,
            [&](std::size_t j) {
              pool.ParallelFor(
                  0, kInner,
                  [&](std::size_t k) {
                    total.fetch_add(1 + i + j + k, std::memory_order_relaxed);
                  },
                  4);
            },
            2);
      },
      1);

  uint64_t expected = 0;
  for (std::size_t i = 0; i < kOuter; ++i) {
    for (std::size_t j = 0; j < kMiddle; ++j) {
      for (std::size_t k = 0; k < kInner; ++k) {
        expected += 1 + i + j + k;
      }
    }
  }
  EXPECT_EQ(total.load(), expected);
}

TEST(ThreadPoolTest, ParallelForExceptionPropagation) {
  ThreadPool pool(4);
  std::atomic<std::size_t> visited{0};

  EXPECT_THROW(
      pool.ParallelFor(
          0, 64,
          [&](std::size_t i) {
            visited.fetch_add(1, std::memory_order_relaxed);
            if (i == 17 || i == 49) {
              throw std::runtime_error("chunk error at " + std::to_string(i));
            }
          },
          4),
      std::runtime_error);

  // Pool remains healthy and usable after an exception.
  std::atomic<int> after_sum{0};
  pool.ParallelFor(0, 100, [&](std::size_t i) {
    after_sum.fetch_add(static_cast<int>(i), std::memory_order_relaxed);
  });
  EXPECT_EQ(after_sum.load(), 4950);
}

TEST(ThreadPoolTest, TransformBatchAndGlobalWrappers) {
  constexpr std::size_t kCount = 4096;
  std::vector<uint32_t> input(kCount);
  std::iota(input.begin(), input.end(), uint32_t{1});
  std::vector<uint64_t> output(kCount, 0);

  ThreadPool pool(4);
  pool.TransformBatch(
      std::span<const uint32_t>(input), std::span<uint64_t>(output),
      [](uint32_t x) { return static_cast<uint64_t>(x) * x; }, 128);

  for (std::size_t i = 0; i < kCount; ++i) {
    const uint64_t v = i + 1;
    EXPECT_EQ(output[i], v * v);
  }

  // Size mismatch throws.
  std::vector<uint64_t> bad_out(10, 0);
  EXPECT_THROW(
      pool.TransformBatch(std::span<const uint32_t>(input),
                          std::span<uint64_t>(bad_out),
                          [](uint32_t x) { return static_cast<uint64_t>(x); }),
      std::invalid_argument);

  // Empty batch succeeds.
  std::vector<uint32_t> empty_in;
  std::vector<uint64_t> empty_out;
  EXPECT_NO_THROW(pool.TransformBatch(
      std::span<const uint32_t>(empty_in), std::span<uint64_t>(empty_out),
      [](uint32_t x) { return static_cast<uint64_t>(x); }));

  // GlobalThreadPool & free convenience wrappers.
  EXPECT_GE(GlobalThreadPool().NumWorkers(), 1u);
  std::vector<uint64_t> global_out(kCount, 0);
  xg::simd::TransformBatch(
      std::span<const uint32_t>(input), std::span<uint64_t>(global_out),
      [](uint32_t x) { return static_cast<uint64_t>(x) * 3ULL; });
  for (std::size_t i = 0; i < kCount; ++i) {
    EXPECT_EQ(global_out[i], static_cast<uint64_t>(input[i]) * 3ULL);
  }

  std::atomic<uint64_t> global_sum{0};
  xg::simd::ParallelFor(0, 1000, [&](std::size_t i) {
    global_sum.fetch_add(i, std::memory_order_relaxed);
  });
  EXPECT_EQ(global_sum.load(), 499500ULL);
}

}  // namespace
