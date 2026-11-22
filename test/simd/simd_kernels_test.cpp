#include <gtest/gtest.h>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/simd/dispatch.hpp"
#include "xgalois/simd/kernels_avx2.hpp"
#include "xgalois/simd/kernels_avx512.hpp"
#include "xgalois/simd/kernels_neon.hpp"
#include "xgalois/simd/kernels_scalar.hpp"

namespace {

using xg::simd::DetectIsaLevel;
using xg::simd::Gf128MulArray;
using xg::simd::Gf28AddArray;
using xg::simd::Gf28FmaArray;
using xg::simd::Gf28MulArray;
using xg::simd::Gf28MulScalar;
using xg::simd::Gf28ScalarMulArray;
using xg::simd::Gf2XorArray;
using xg::simd::GfPrimeAddArray;
using xg::simd::GfPrimeSubArray;
using xg::simd::IsaLevel;
using xg::simd::IsaLevelName;
using xg::simd::ScopedIsaOverride;

constexpr std::size_t kTestSizes[] = {
    0, 1, 7, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 1025,
};

// Returns all ISA levels safely executable on the current CPU.
std::vector<IsaLevel> SupportedIsaLevels() {
  std::vector<IsaLevel> levels = {IsaLevel::kScalar};
  const IsaLevel detected = DetectIsaLevel();
  switch (detected) {
    case IsaLevel::kAvx512Gfni:
      levels.push_back(IsaLevel::kAvx2);
      levels.push_back(IsaLevel::kAvx512Gfni);
      break;
    case IsaLevel::kAvx2:
      levels.push_back(IsaLevel::kAvx2);
      break;
    case IsaLevel::kNeonPmull:
      levels.push_back(IsaLevel::kNeon);
      levels.push_back(IsaLevel::kNeonPmull);
      break;
    case IsaLevel::kNeon:
      levels.push_back(IsaLevel::kNeon);
      break;
    case IsaLevel::kSse42Pclmul:
      levels.push_back(IsaLevel::kSse42Pclmul);
      break;
    case IsaLevel::kScalar:
      break;
  }
  return levels;
}

// Software reference model of Intel VGF2P8AFFINEQB with imm8 = 0:
//   out.bit[i] = parity(M.byte[7 - i] & x) for i in 0..7.
uint8_t EmulateGf2p8Affine(uint8_t x, uint64_t matrix) noexcept {
  uint8_t out = 0u;
  for (int i = 0; i < 8; ++i) {
    const uint8_t row =
        static_cast<uint8_t>((matrix >> (8 * (7 - i))) & 0xFFULL);
    const uint8_t parity =
        static_cast<uint8_t>(std::popcount(static_cast<unsigned>(row & x)) & 1);
    out |= static_cast<uint8_t>(parity << i);
  }
  return out;
}

TEST(SimdKernelsTest, Gf2XorArrayBitIdenticalAcrossSupportedIsas) {
  std::mt19937_64 rng(0x1122334455667788ULL);

  for (IsaLevel level : SupportedIsaLevels()) {
    SCOPED_TRACE(std::string(IsaLevelName(level)));
    ScopedIsaOverride guard(level);

    for (std::size_t n : kTestSizes) {
      // Allocate +1 element and slice at offset 1 to exercise unaligned tails.
      std::vector<uint64_t> a_buf(n + 1);
      std::vector<uint64_t> b_buf(n + 1);
      for (std::size_t i = 0; i < n + 1; ++i) {
        a_buf[i] = rng();
        b_buf[i] = rng();
      }
      const std::span<const uint64_t> a(a_buf.data() + 1, n);
      const std::span<const uint64_t> b(b_buf.data() + 1, n);

      std::vector<uint64_t> expected(n, 0ULL);
      xg::simd::scalar::Gf2XorArray(expected, a, b);

      std::vector<uint64_t> actual_buf(n + 1, 0xDEADBEEFULL);
      const std::span<uint64_t> actual(actual_buf.data() + 1, n);
      Gf2XorArray(actual, a, b);

      EXPECT_EQ(actual_buf[0], 0xDEADBEEFULL);
      for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQ(actual[i], expected[i]) << "n=" << n << " i=" << i;
      }

      // In-place aliasing dst == a.
      std::vector<uint64_t> in_place(a.begin(), a.end());
      Gf2XorArray(in_place, in_place, b);
      EXPECT_EQ(in_place, expected);
    }

    std::vector<uint64_t> v1(4, 1ULL);
    std::vector<uint64_t> v2(3, 1ULL);
    EXPECT_THROW(Gf2XorArray(v2, v1, v1), std::invalid_argument);
  }
}

TEST(SimdKernelsTest, Gf28AddArrayBitIdenticalAcrossSupportedIsas) {
  std::mt19937_64 rng(0x2233445566778899ULL);

  for (IsaLevel level : SupportedIsaLevels()) {
    SCOPED_TRACE(std::string(IsaLevelName(level)));
    ScopedIsaOverride guard(level);

    for (std::size_t n : kTestSizes) {
      std::vector<uint8_t> a_buf(n + 3);
      std::vector<uint8_t> b_buf(n + 3);
      for (std::size_t i = 0; i < n + 3; ++i) {
        a_buf[i] = static_cast<uint8_t>(rng() & 0xFFu);
        b_buf[i] = static_cast<uint8_t>(rng() & 0xFFu);
      }
      const std::span<const uint8_t> a(a_buf.data() + 1, n);
      const std::span<const uint8_t> b(b_buf.data() + 2, n);

      std::vector<uint8_t> expected(n, 0);
      xg::simd::scalar::Gf28AddArray(expected, a, b);

      std::vector<uint8_t> actual_buf(n + 3, 0xA5);
      const std::span<uint8_t> actual(actual_buf.data() + 1, n);
      Gf28AddArray(actual, a, b);

      EXPECT_EQ(actual_buf[0], 0xA5);
      EXPECT_EQ(actual_buf[n + 1], 0xA5);
      for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQ(actual[i], expected[i]) << "n=" << n << " i=" << i;
      }
    }

    std::vector<uint8_t> v1(16, 1);
    std::vector<uint8_t> v2(15, 1);
    EXPECT_THROW(Gf28AddArray(v2, v1, v1), std::invalid_argument);
  }
}

TEST(SimdKernelsTest, Gf28ScalarMulAndFmaExhaustiveOn0x11DAnd0x11B) {
  std::vector<uint8_t> all_bytes(256);
  std::iota(all_bytes.begin(), all_bytes.end(), uint8_t{0});

  std::vector<uint8_t> init_acc(256);
  for (std::size_t i = 0; i < 256; ++i) {
    init_acc[i] = static_cast<uint8_t>((i * 73u + 19u) & 0xFFu);
  }

  for (IsaLevel level : SupportedIsaLevels()) {
    SCOPED_TRACE(std::string(IsaLevelName(level)));
    ScopedIsaOverride guard(level);

    for (uint16_t poly : {uint16_t{0x11D}, uint16_t{0x11B}}) {
      // Exhaustive: all 256 scalars x all 256 byte values.
      for (unsigned s = 0; s < 256; ++s) {
        const uint8_t scalar = static_cast<uint8_t>(s);

        std::vector<uint8_t> expected_mul(256);
        xg::simd::scalar::Gf28ScalarMulArray(expected_mul, all_bytes, scalar,
                                             poly);
        std::vector<uint8_t> actual_mul(256, 0xCC);
        Gf28ScalarMulArray(actual_mul, all_bytes, scalar, poly);
        EXPECT_EQ(actual_mul, expected_mul)
            << "poly=0x" << std::hex << poly << " scalar=" << s;

        std::vector<uint8_t> expected_fma = init_acc;
        xg::simd::scalar::Gf28FmaArray(expected_fma, all_bytes, scalar, poly);
        std::vector<uint8_t> actual_fma = init_acc;
        Gf28FmaArray(actual_fma, all_bytes, scalar, poly);
        EXPECT_EQ(actual_fma, expected_fma)
            << "poly=0x" << std::hex << poly << " scalar=" << s;
      }

      // All odd/unaligned/tail sizes with representative scalars.
      std::mt19937_64 rng(0x33445566778899AAULL ^ poly);
      for (std::size_t n : kTestSizes) {
        std::vector<uint8_t> src_buf(n + 3);
        std::vector<uint8_t> acc_init(n + 3);
        for (std::size_t i = 0; i < n + 3; ++i) {
          src_buf[i] = static_cast<uint8_t>(rng() & 0xFFu);
          acc_init[i] = static_cast<uint8_t>(rng() & 0xFFu);
        }
        const std::span<const uint8_t> src(src_buf.data() + 1, n);

        for (uint8_t scalar : {uint8_t{0x00}, uint8_t{0x01}, uint8_t{0x02},
                               uint8_t{0x57}, uint8_t{0xA9}, uint8_t{0xFF}}) {
          std::vector<uint8_t> expected_mul(n);
          xg::simd::scalar::Gf28ScalarMulArray(expected_mul, src, scalar, poly);

          std::vector<uint8_t> actual_mul_buf(n + 3, 0xEE);
          const std::span<uint8_t> actual_mul(actual_mul_buf.data() + 1, n);
          Gf28ScalarMulArray(actual_mul, src, scalar, poly);
          EXPECT_EQ(actual_mul_buf[0], 0xEE);
          EXPECT_EQ(actual_mul_buf[n + 1], 0xEE);
          for (std::size_t i = 0; i < n; ++i) {
            EXPECT_EQ(actual_mul[i], expected_mul[i]);
          }

          // In-place Gf28ScalarMulArray.
          std::vector<uint8_t> in_place(src.begin(), src.end());
          Gf28ScalarMulArray(in_place, in_place, scalar, poly);
          EXPECT_EQ(in_place, expected_mul);

          // Unaligned Gf28FmaArray.
          std::vector<uint8_t> expected_fma(acc_init.begin() + 2,
                                            acc_init.begin() + 2 + n);
          xg::simd::scalar::Gf28FmaArray(expected_fma, src, scalar, poly);

          std::vector<uint8_t> actual_fma_buf = acc_init;
          const std::span<uint8_t> actual_fma(actual_fma_buf.data() + 2, n);
          Gf28FmaArray(actual_fma, src, scalar, poly);
          EXPECT_EQ(actual_fma_buf[1], acc_init[1]);
          EXPECT_EQ(actual_fma_buf[n + 2], acc_init[n + 2]);
          for (std::size_t i = 0; i < n; ++i) {
            EXPECT_EQ(actual_fma[i], expected_fma[i]);
          }
        }
      }
    }

    std::vector<uint8_t> v1(16, 1);
    std::vector<uint8_t> v2(15, 1);
    EXPECT_THROW(Gf28ScalarMulArray(v2, v1, 0x57, 0x11D),
                 std::invalid_argument);
    EXPECT_THROW(Gf28FmaArray(v2, v1, 0x57, 0x11D), std::invalid_argument);
  }
}

TEST(SimdKernelsTest, Gf28MulArrayExhaustiveAndTailSizesOn0x11DAnd0x11B) {
  for (IsaLevel level : SupportedIsaLevels()) {
    SCOPED_TRACE(std::string(IsaLevelName(level)));
    ScopedIsaOverride guard(level);

    for (uint16_t poly : {uint16_t{0x11D}, uint16_t{0x11B}}) {
      // Exhaustive 256x256 check batched in 256-byte rows.
      std::vector<uint8_t> lhs(256);
      std::iota(lhs.begin(), lhs.end(), uint8_t{0});
      std::vector<uint8_t> rhs(256);
      std::vector<uint8_t> actual(256);
      std::vector<uint8_t> expected(256);

      for (unsigned b = 0; b < 256; ++b) {
        std::fill(rhs.begin(), rhs.end(), static_cast<uint8_t>(b));
        xg::simd::scalar::Gf28MulArray(expected, lhs, rhs, poly);
        Gf28MulArray(actual, lhs, rhs, poly);
        EXPECT_EQ(actual, expected)
            << "poly=0x" << std::hex << poly << " b=" << b;
      }

      // Tail and unaligned sizes.
      std::mt19937_64 rng(0x445566778899AABBULL ^ poly);
      for (std::size_t n : kTestSizes) {
        std::vector<uint8_t> a_buf(n + 2);
        std::vector<uint8_t> b_buf(n + 2);
        for (std::size_t i = 0; i < n + 2; ++i) {
          a_buf[i] = static_cast<uint8_t>(rng() & 0xFFu);
          b_buf[i] = static_cast<uint8_t>(rng() & 0xFFu);
        }
        const std::span<const uint8_t> a(a_buf.data() + 1, n);
        const std::span<const uint8_t> b(b_buf.data() + 1, n);

        std::vector<uint8_t> exp(n);
        xg::simd::scalar::Gf28MulArray(exp, a, b, poly);

        std::vector<uint8_t> act(n, 0xFF);
        Gf28MulArray(act, a, b, poly);
        EXPECT_EQ(act, exp);
      }
    }

    std::vector<uint8_t> v1(16, 1);
    std::vector<uint8_t> v2(15, 1);
    EXPECT_THROW(Gf28MulArray(v2, v1, v1, 0x11B), std::invalid_argument);
  }
}

TEST(SimdKernelsTest, Avx512GfniMatrixMatchesScalarExhaustively) {
  // Verifies the GFNI 8x8 bit matrix bit-ordering against scalar multiplication
  // for all 256 scalars x all 256 byte values on both 0x11D and 0x11B, even on
  // CPUs without AVX-512/GFNI hardware.
  for (uint16_t poly : {uint16_t{0x11D}, uint16_t{0x11B}}) {
    for (unsigned s = 0; s < 256; ++s) {
      const uint8_t scalar = static_cast<uint8_t>(s);
      const uint64_t mat =
          xg::simd::avx512::detail::BuildGfniMulMatrix(scalar, poly);
      for (unsigned x = 0; x < 256; ++x) {
        const uint8_t actual = EmulateGf2p8Affine(static_cast<uint8_t>(x), mat);
        const uint8_t expected =
            Gf28MulScalar(static_cast<uint8_t>(x), scalar, poly);
        ASSERT_EQ(actual, expected)
            << "poly=0x" << std::hex << poly << " scalar=" << s << " x=" << x;
      }
    }
  }
}

TEST(SimdKernelsTest, Gf128MulArrayBitIdenticalAcrossSupportedIsas) {
  constexpr uint64_t kReductionPolys[] = {
      0x87ULL,                                 // GCM: x^128 + x^7 + x^2 + x + 1
      (1ULL << 63) | 0x87ULL,                  // High-bit reduction polynomial
      0x8000000000000000ULL | 0x123456789ULL,  // Arbitrary 64-bit high-bit poly
  };

  std::mt19937_64 rng(0x5566778899AABBCCULL);

  for (IsaLevel level : SupportedIsaLevels()) {
    SCOPED_TRACE(std::string(IsaLevelName(level)));
    ScopedIsaOverride guard(level);

    for (uint64_t poly_lo : kReductionPolys) {
      for (std::size_t n : kTestSizes) {
        std::vector<std::pair<uint64_t, uint64_t>> a(n);
        std::vector<std::pair<uint64_t, uint64_t>> b(n);
        for (std::size_t i = 0; i < n; ++i) {
          a[i] = {rng(), rng()};
          b[i] = {rng(), rng()};
        }
        if (n >= 4) {
          a[0] = {0ULL, 0ULL};
          b[0] = {~0ULL, ~0ULL};
          a[1] = {1ULL, 0ULL};
          b[1] = {0x0123456789ABCDEFULL, 0xFEDCBA9876543210ULL};
          a[2] = {0ULL, 1ULL << 63};
          b[2] = {0ULL, 1ULL << 63};
          a[3] = {~0ULL, ~0ULL};
          b[3] = {~0ULL, ~0ULL};
        }

        std::vector<std::pair<uint64_t, uint64_t>> expected(n);
        xg::simd::scalar::Gf128MulArray(expected, a, b, poly_lo);

        std::vector<std::pair<uint64_t, uint64_t>> actual(n);
        Gf128MulArray(actual, a, b, poly_lo);
        EXPECT_EQ(actual, expected)
            << "n=" << n << " poly=0x" << std::hex << poly_lo;

        // In-place aliasing dst == a.
        std::vector<std::pair<uint64_t, uint64_t>> in_place = a;
        Gf128MulArray(in_place, in_place, b, poly_lo);
        EXPECT_EQ(in_place, expected);
      }
    }

    std::vector<std::pair<uint64_t, uint64_t>> v1(4);
    std::vector<std::pair<uint64_t, uint64_t>> v2(3);
    EXPECT_THROW(Gf128MulArray(v2, v1, v1, 0x87ULL), std::invalid_argument);
  }
}

TEST(SimdKernelsTest, GfPrimeAddSubArrayBitIdenticalAcrossSupportedIsas) {
  constexpr uint64_t kModuli[] = {
      257ULL, 998244353ULL, (1ULL << 61) - 1ULL,
      UINT64_MAX - 58ULL,  // 2^64 - 59: tests 64-bit unsigned overflow handling
  };

  std::mt19937_64 rng(0x66778899AABBCCDDULL);

  for (IsaLevel level : SupportedIsaLevels()) {
    SCOPED_TRACE(std::string(IsaLevelName(level)));
    ScopedIsaOverride guard(level);

    for (uint64_t mod : kModuli) {
      for (std::size_t n : kTestSizes) {
        std::vector<uint64_t> a_buf(n + 1);
        std::vector<uint64_t> b_buf(n + 1);
        for (std::size_t i = 0; i < n + 1; ++i) {
          a_buf[i] = rng() % mod;
          b_buf[i] = rng() % mod;
        }
        const std::span<uint64_t> a_mut(a_buf.data() + 1, n);
        const std::span<uint64_t> b_mut(b_buf.data() + 1, n);
        if (n >= 6) {
          a_mut[0] = 0ULL;
          b_mut[0] = 0ULL;
          a_mut[1] = mod - 1ULL;
          b_mut[1] = mod - 1ULL;
          a_mut[2] = 0ULL;
          b_mut[2] = mod - 1ULL;
          a_mut[3] = mod - 1ULL;
          b_mut[3] = 0ULL;
          a_mut[4] = mod - 1ULL;
          b_mut[4] = 1ULL;
          a_mut[5] = 1ULL;
          b_mut[5] = mod - 1ULL;
        }
        const std::span<const uint64_t> a = a_mut;
        const std::span<const uint64_t> b = b_mut;

        std::vector<uint64_t> expected_add(n);
        std::vector<uint64_t> expected_sub(n);
        xg::simd::scalar::GfPrimeAddArray(expected_add, a, b, mod);
        xg::simd::scalar::GfPrimeSubArray(expected_sub, a, b, mod);

        std::vector<uint64_t> actual_add_buf(n + 1, 0xCAFEBABEULL);
        std::vector<uint64_t> actual_sub_buf(n + 1, 0xCAFEBABEULL);
        const std::span<uint64_t> actual_add(actual_add_buf.data() + 1, n);
        const std::span<uint64_t> actual_sub(actual_sub_buf.data() + 1, n);

        GfPrimeAddArray(actual_add, a, b, mod);
        GfPrimeSubArray(actual_sub, a, b, mod);

        EXPECT_EQ(actual_add_buf[0], 0xCAFEBABEULL);
        EXPECT_EQ(actual_sub_buf[0], 0xCAFEBABEULL);
        for (std::size_t i = 0; i < n; ++i) {
          EXPECT_EQ(actual_add[i], expected_add[i])
              << "mod=" << mod << " n=" << n << " i=" << i;
          EXPECT_EQ(actual_sub[i], expected_sub[i])
              << "mod=" << mod << " n=" << n << " i=" << i;
        }

        // In-place aliasing dst == a.
        std::vector<uint64_t> in_place_add(a.begin(), a.end());
        GfPrimeAddArray(in_place_add, in_place_add, b, mod);
        EXPECT_EQ(in_place_add, expected_add);

        std::vector<uint64_t> in_place_sub(a.begin(), a.end());
        GfPrimeSubArray(in_place_sub, in_place_sub, b, mod);
        EXPECT_EQ(in_place_sub, expected_sub);
      }
    }

    std::vector<uint64_t> v1(8, 1ULL);
    std::vector<uint64_t> v2(7, 1ULL);
    EXPECT_THROW(GfPrimeAddArray(v1, v1, v1, 0ULL), std::invalid_argument);
    EXPECT_THROW(GfPrimeSubArray(v1, v1, v1, 0ULL), std::invalid_argument);
    EXPECT_THROW(GfPrimeAddArray(v2, v1, v1, 257ULL), std::invalid_argument);
    EXPECT_THROW(GfPrimeSubArray(v2, v1, v1, 257ULL), std::invalid_argument);
  }
}

}  // namespace
