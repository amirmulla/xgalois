#ifndef XGALOIS_SIMD_KERNELS_NEON_HPP_
#define XGALOIS_SIMD_KERNELS_NEON_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

#if defined(__aarch64__)
#include <arm_neon.h>
#define XGALOIS_SIMD_HAS_NEON 1
#if defined(__ARM_FEATURE_CRYPTO) || defined(__GNUC__) || defined(__clang__)
#define XGALOIS_SIMD_HAS_NEON_PMULL 1
#else
#define XGALOIS_SIMD_HAS_NEON_PMULL 0
#endif
#else
#define XGALOIS_SIMD_HAS_NEON 0
#define XGALOIS_SIMD_HAS_NEON_PMULL 0
#endif

namespace xg {
namespace simd {
namespace neon {

#if XGALOIS_SIMD_HAS_NEON

namespace detail {

struct Gf28NibbleLuts {
  alignas(16) std::array<uint8_t, 16> lo;
  alignas(16) std::array<uint8_t, 16> hi;
};

constexpr uint8_t Gf28MulByX(uint8_t value, uint8_t poly_lo) noexcept {
  const uint8_t mask =
      static_cast<uint8_t>(0u - static_cast<unsigned>(value >> 7));
  return static_cast<uint8_t>((static_cast<unsigned>(value) << 1) ^
                              (poly_lo & mask));
}

constexpr Gf28NibbleLuts BuildGf28NibbleLuts(uint8_t scalar,
                                             uint16_t irred_poly) noexcept {
  const uint8_t poly_lo = static_cast<uint8_t>(irred_poly & 0xFFu);
  const uint8_t b0 = scalar;
  const uint8_t b1 = Gf28MulByX(b0, poly_lo);
  const uint8_t b2 = Gf28MulByX(b1, poly_lo);
  const uint8_t b3 = Gf28MulByX(b2, poly_lo);
  const uint8_t b4 = Gf28MulByX(b3, poly_lo);
  const uint8_t b5 = Gf28MulByX(b4, poly_lo);
  const uint8_t b6 = Gf28MulByX(b5, poly_lo);
  const uint8_t b7 = Gf28MulByX(b6, poly_lo);

  Gf28NibbleLuts luts{};
  luts.lo[0] = 0;
  luts.hi[0] = 0;
  for (unsigned idx = 1; idx < 16u; ++idx) {
    uint8_t lo_val = 0;
    uint8_t hi_val = 0;
    if ((idx & 1u) != 0u) {
      lo_val ^= b0;
      hi_val ^= b4;
    }
    if ((idx & 2u) != 0u) {
      lo_val ^= b1;
      hi_val ^= b5;
    }
    if ((idx & 4u) != 0u) {
      lo_val ^= b2;
      hi_val ^= b6;
    }
    if ((idx & 8u) != 0u) {
      lo_val ^= b3;
      hi_val ^= b7;
    }
    luts.lo[idx] = lo_val;
    luts.hi[idx] = hi_val;
  }
  return luts;
}

#if XGALOIS_SIMD_HAS_NEON_PMULL
#if !defined(__ARM_FEATURE_CRYPTO) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("+crypto")))
#endif
inline std::pair<uint64_t, uint64_t>
Clmul64(uint64_t a, uint64_t b) noexcept {
  const poly128_t prod =
      vmull_p64(static_cast<poly64_t>(a), static_cast<poly64_t>(b));
  const uint64x2_t vec = vreinterpretq_u64_p128(prod);
  return {vgetq_lane_u64(vec, 0), vgetq_lane_u64(vec, 1)};
}
#endif

}  // namespace detail

inline void Gf2XorArray(std::span<uint64_t> dst, std::span<const uint64_t> a,
                        std::span<const uint64_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf2XorArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint64_t* dst_ptr = dst.data();
  const uint64_t* a_ptr = a.data();
  const uint64_t* b_ptr = b.data();

  std::size_t i = 0;
  for (; i + 2 <= n; i += 2) {
    const uint64x2_t va = vld1q_u64(a_ptr + i);
    const uint64x2_t vb = vld1q_u64(b_ptr + i);
    vst1q_u64(dst_ptr + i, veorq_u64(va, vb));
  }
  for (; i < n; ++i) {
    dst_ptr[i] = a_ptr[i] ^ b_ptr[i];
  }
}

inline void Gf28AddArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                         std::span<const uint8_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf28AddArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* a_ptr = a.data();
  const uint8_t* b_ptr = b.data();

  std::size_t i = 0;
  for (; i + 16 <= n; i += 16) {
    const uint8x16_t va = vld1q_u8(a_ptr + i);
    const uint8x16_t vb = vld1q_u8(b_ptr + i);
    vst1q_u8(dst_ptr + i, veorq_u8(va, vb));
  }
  for (; i < n; ++i) {
    dst_ptr[i] = static_cast<uint8_t>(a_ptr[i] ^ b_ptr[i]);
  }
}

inline void Gf28ScalarMulArray(std::span<uint8_t> dst,
                               std::span<const uint8_t> src, uint8_t scalar,
                               uint16_t irred_poly) {
  if (dst.size() != src.size()) {
    throw std::invalid_argument("Gf28ScalarMulArray: span sizes must match");
  }
  if (scalar == 0u) {
    std::fill(dst.begin(), dst.end(), uint8_t{0});
    return;
  }
  if (scalar == 1u) {
    if (dst.data() != src.data()) {
      std::copy(src.begin(), src.end(), dst.begin());
    }
    return;
  }

  const detail::Gf28NibbleLuts luts =
      detail::BuildGf28NibbleLuts(scalar, irred_poly);
  const uint8x16_t v_lo_lut = vld1q_u8(luts.lo.data());
  const uint8x16_t v_hi_lut = vld1q_u8(luts.hi.data());
  const uint8x16_t v_mask = vdupq_n_u8(0x0F);

  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* src_ptr = src.data();

  std::size_t i = 0;
  for (; i + 16 <= n; i += 16) {
    const uint8x16_t v = vld1q_u8(src_ptr + i);
    const uint8x16_t lo_nib = vandq_u8(v, v_mask);
    const uint8x16_t hi_nib = vshrq_n_u8(v, 4);
    const uint8x16_t res =
        veorq_u8(vqtbl1q_u8(v_lo_lut, lo_nib), vqtbl1q_u8(v_hi_lut, hi_nib));
    vst1q_u8(dst_ptr + i, res);
  }
  for (; i < n; ++i) {
    const uint8_t v = src_ptr[i];
    dst_ptr[i] = static_cast<uint8_t>(luts.lo[v & 0x0Fu] ^ luts.hi[v >> 4]);
  }
}

inline void Gf28FmaArray(std::span<uint8_t> dst, std::span<const uint8_t> src,
                         uint8_t scalar, uint16_t irred_poly) {
  if (dst.size() != src.size()) {
    throw std::invalid_argument("Gf28FmaArray: span sizes must match");
  }
  if (scalar == 0u) {
    return;
  }
  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* src_ptr = src.data();

  if (scalar == 1u) {
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
      const uint8x16_t vd = vld1q_u8(dst_ptr + i);
      const uint8x16_t vs = vld1q_u8(src_ptr + i);
      vst1q_u8(dst_ptr + i, veorq_u8(vd, vs));
    }
    for (; i < n; ++i) {
      dst_ptr[i] = static_cast<uint8_t>(dst_ptr[i] ^ src_ptr[i]);
    }
    return;
  }

  const detail::Gf28NibbleLuts luts =
      detail::BuildGf28NibbleLuts(scalar, irred_poly);
  const uint8x16_t v_lo_lut = vld1q_u8(luts.lo.data());
  const uint8x16_t v_hi_lut = vld1q_u8(luts.hi.data());
  const uint8x16_t v_mask = vdupq_n_u8(0x0F);

  std::size_t i = 0;
  for (; i + 16 <= n; i += 16) {
    const uint8x16_t v = vld1q_u8(src_ptr + i);
    const uint8x16_t vd = vld1q_u8(dst_ptr + i);
    const uint8x16_t lo_nib = vandq_u8(v, v_mask);
    const uint8x16_t hi_nib = vshrq_n_u8(v, 4);
    const uint8x16_t prod =
        veorq_u8(vqtbl1q_u8(v_lo_lut, lo_nib), vqtbl1q_u8(v_hi_lut, hi_nib));
    vst1q_u8(dst_ptr + i, veorq_u8(vd, prod));
  }
  for (; i < n; ++i) {
    const uint8_t v = src_ptr[i];
    dst_ptr[i] =
        static_cast<uint8_t>(dst_ptr[i] ^ luts.lo[v & 0x0Fu] ^ luts.hi[v >> 4]);
  }
}

#if XGALOIS_SIMD_HAS_NEON_PMULL
#if !defined(__ARM_FEATURE_CRYPTO) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("+crypto")))
#endif
inline void
Gf128MulArray(std::span<std::pair<uint64_t, uint64_t>> dst,
              std::span<const std::pair<uint64_t, uint64_t>> a,
              std::span<const std::pair<uint64_t, uint64_t>> b,
              uint64_t reduction_poly_lo) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf128MulArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    const uint64_t a0 = a[i].first;
    const uint64_t a1 = a[i].second;
    const uint64_t b0 = b[i].first;
    const uint64_t b1 = b[i].second;

    const auto z0 = detail::Clmul64(a0, b0);
    const auto z2 = detail::Clmul64(a1, b1);
    const auto zm = detail::Clmul64(a0 ^ a1, b0 ^ b1);

    const uint64_t z1_lo = zm.first ^ z0.first ^ z2.first;
    const uint64_t z1_hi = zm.second ^ z0.second ^ z2.second;

    uint64_t c0 = z0.first;
    uint64_t c1 = z0.second ^ z1_lo;
    uint64_t c2 = z1_hi ^ z2.first;
    const uint64_t c3 = z2.second;

    const auto t3 = detail::Clmul64(c3, reduction_poly_lo);
    c1 ^= t3.first;
    c2 ^= t3.second;

    const auto t2 = detail::Clmul64(c2, reduction_poly_lo);
    c0 ^= t2.first;
    c1 ^= t2.second;

    dst[i] = {c0, c1};
  }
}
#endif  // XGALOIS_SIMD_HAS_NEON_PMULL

inline void GfPrimeAddArray(std::span<uint64_t> dst,
                            std::span<const uint64_t> a,
                            std::span<const uint64_t> b, uint64_t mod) {
  if (mod == 0ULL) {
    throw std::invalid_argument("GfPrimeAddArray: mod must be positive");
  }
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("GfPrimeAddArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint64_t* dst_ptr = dst.data();
  const uint64_t* a_ptr = a.data();
  const uint64_t* b_ptr = b.data();

  const uint64x2_t v_mod = vdupq_n_u64(mod);
  std::size_t i = 0;
  for (; i + 2 <= n; i += 2) {
    const uint64x2_t va = vld1q_u64(a_ptr + i);
    const uint64x2_t vb = vld1q_u64(b_ptr + i);
    const uint64x2_t v_diff = vsubq_u64(v_mod, vb);
    const uint64x2_t v_ge = vcgeq_u64(va, v_diff);
    const uint64x2_t v_sum = vaddq_u64(va, vb);
    const uint64x2_t v_res = vsubq_u64(v_sum, vandq_u64(v_ge, v_mod));
    vst1q_u64(dst_ptr + i, v_res);
  }
  for (; i < n; ++i) {
    const uint64_t av = a_ptr[i];
    const uint64_t bv = b_ptr[i];
    const uint64_t diff = mod - bv;
    dst_ptr[i] = (av >= diff) ? (av - diff) : (av + bv);
  }
}

inline void GfPrimeSubArray(std::span<uint64_t> dst,
                            std::span<const uint64_t> a,
                            std::span<const uint64_t> b, uint64_t mod) {
  if (mod == 0ULL) {
    throw std::invalid_argument("GfPrimeSubArray: mod must be positive");
  }
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("GfPrimeSubArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint64_t* dst_ptr = dst.data();
  const uint64_t* a_ptr = a.data();
  const uint64_t* b_ptr = b.data();

  const uint64x2_t v_mod = vdupq_n_u64(mod);
  std::size_t i = 0;
  for (; i + 2 <= n; i += 2) {
    const uint64x2_t va = vld1q_u64(a_ptr + i);
    const uint64x2_t vb = vld1q_u64(b_ptr + i);
    const uint64x2_t v_lt = vcltq_u64(va, vb);
    const uint64x2_t v_diff = vsubq_u64(va, vb);
    const uint64x2_t v_res = vaddq_u64(v_diff, vandq_u64(v_lt, v_mod));
    vst1q_u64(dst_ptr + i, v_res);
  }
  for (; i < n; ++i) {
    const uint64_t av = a_ptr[i];
    const uint64_t bv = b_ptr[i];
    dst_ptr[i] = (av >= bv) ? (av - bv) : (av + (mod - bv));
  }
}

#endif  // XGALOIS_SIMD_HAS_NEON

}  // namespace neon
}  // namespace simd
}  // namespace xg

#endif  // XGALOIS_SIMD_KERNELS_NEON_HPP_
