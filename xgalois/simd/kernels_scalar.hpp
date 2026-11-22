#ifndef XGALOIS_SIMD_KERNELS_SCALAR_HPP_
#define XGALOIS_SIMD_KERNELS_SCALAR_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

#include "xgalois/simd/dispatch.hpp"

namespace xg {
namespace simd {

namespace detail {

#if defined(__SIZEOF_INT128__)
__extension__ using Uint128 = unsigned __int128;
#endif

// Computes the high 64 bits of the unsigned 64x64 -> 128-bit product a * b.
constexpr uint64_t MulHi64(uint64_t a, uint64_t b) noexcept {
#if defined(__SIZEOF_INT128__)
  return static_cast<uint64_t>((static_cast<Uint128>(a) * b) >> 64);
#else
  const uint64_t a_lo = a & 0xFFFFFFFFULL;
  const uint64_t a_hi = a >> 32;
  const uint64_t b_lo = b & 0xFFFFFFFFULL;
  const uint64_t b_hi = b >> 32;

  const uint64_t p0 = a_lo * b_lo;
  const uint64_t p1 = a_lo * b_hi;
  const uint64_t p2 = a_hi * b_lo;
  const uint64_t p3 = a_hi * b_hi;

  const uint64_t mid = (p0 >> 32) + (p1 & 0xFFFFFFFFULL) + (p2 & 0xFFFFFFFFULL);
  return p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
#endif
}

// Computes floor(((hi << 64) | lo) / divisor) assuming hi < divisor < 2^62.
constexpr uint64_t Div128By64(uint64_t hi, uint64_t lo,
                              uint64_t divisor) noexcept {
#if defined(__SIZEOF_INT128__)
  const Uint128 dividend = (static_cast<Uint128>(hi) << 64) | lo;
  return static_cast<uint64_t>(dividend / divisor);
#else
  uint64_t rem = hi;
  uint64_t quot = 0;
  for (int bit = 63; bit >= 0; --bit) {
    rem = (rem << 1) | ((lo >> static_cast<unsigned>(bit)) & 1ULL);
    if (rem >= divisor) {
      rem -= divisor;
      quot |= (1ULL << static_cast<unsigned>(bit));
    }
  }
  return quot;
#endif
}

struct Gf28NibbleLuts {
  std::array<uint8_t, 16> lo;
  std::array<uint8_t, 16> hi;
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

}  // namespace detail

// Multiplies two GF(2^8) elements modulo an arbitrary degree-8 irreducible
// polynomial (with bit 8 set, e.g. 0x11D or AES 0x11B).
constexpr uint8_t Gf28MulScalar(uint8_t a, uint8_t b,
                                uint16_t irred_poly) noexcept {
  const uint8_t poly_lo = static_cast<uint8_t>(irred_poly & 0xFFu);
  uint8_t result = 0;
  uint8_t lhs = a;
  uint8_t rhs = b;
  for (int bit = 0; bit < 8; ++bit) {
    const uint8_t add_mask =
        static_cast<uint8_t>(0u - static_cast<unsigned>(rhs & 1u));
    result ^= static_cast<uint8_t>(lhs & add_mask);
    lhs = detail::Gf28MulByX(lhs, poly_lo);
    rhs = static_cast<uint8_t>(rhs >> 1);
  }
  return result;
}

// Carry-less 64x64 -> 128-bit polynomial multiplication over GF(2).
// Returns {lo, hi} representing lo + hi * x^64.
constexpr std::pair<uint64_t, uint64_t> Gf128ClmulScalar(uint64_t a,
                                                         uint64_t b) noexcept {
  uint64_t lo = ((b & 1ULL) != 0ULL) ? a : 0ULL;
  uint64_t hi = 0ULL;
  for (unsigned shift = 1u; shift < 64u; ++shift) {
    const uint64_t mask = 0ULL - ((b >> shift) & 1ULL);
    lo ^= (a << shift) & mask;
    hi ^= (a >> (64u - shift)) & mask;
  }
  return {lo, hi};
}

// Multiplies two GF(2^128) elements represented as {lo, hi} modulo
// x^128 + reduction_poly_lo (e.g. reduction_poly_lo = 0x87 for GCM).
constexpr std::pair<uint64_t, uint64_t> Gf128MulScalar(
    std::pair<uint64_t, uint64_t> a, std::pair<uint64_t, uint64_t> b,
    uint64_t reduction_poly_lo) noexcept {
  const uint64_t a0 = a.first;
  const uint64_t a1 = a.second;
  const uint64_t b0 = b.first;
  const uint64_t b1 = b.second;

  // Karatsuba 128x128 carry-less product (3 clmuls).
  const auto z0 = Gf128ClmulScalar(a0, b0);
  const auto z2 = Gf128ClmulScalar(a1, b1);
  const auto zm = Gf128ClmulScalar(a0 ^ a1, b0 ^ b1);

  const uint64_t z1_lo = zm.first ^ z0.first ^ z2.first;
  const uint64_t z1_hi = zm.second ^ z0.second ^ z2.second;

  uint64_t c0 = z0.first;
  uint64_t c1 = z0.second ^ z1_lo;
  uint64_t c2 = z1_hi ^ z2.first;
  const uint64_t c3 = z2.second;

  // Reduce c3 * x^192 = (c3 * R) * x^64, then c2 * x^128 = c2 * R.
  // Since deg(c3) <= 62 and deg(R) <= 63, deg(c3 * R) <= 125 so t3_hi fits in
  // limb 2; then deg(c2 * R) <= 126 < 128 so t2 fits completely in limbs 0..1.
  const auto t3 = Gf128ClmulScalar(c3, reduction_poly_lo);
  c1 ^= t3.first;
  c2 ^= t3.second;

  const auto t2 = Gf128ClmulScalar(c2, reduction_poly_lo);
  c0 ^= t2.first;
  c1 ^= t2.second;

  return {c0, c1};
}

// Precomputes Shoup's modular multiplication factor:
//   floor((b * 2^64) / mod)
// Requires 0 < mod < 2^62 and b < mod.
constexpr uint64_t ComputeShoupFactor(uint64_t b, uint64_t mod) {
  if (mod == 0ULL || mod >= (1ULL << 62)) {
    throw std::invalid_argument(
        "ComputeShoupFactor: mod must satisfy 0 < mod < 2^62");
  }
  if (b >= mod) {
    throw std::invalid_argument("ComputeShoupFactor: b must satisfy b < mod");
  }
  return detail::Div128By64(b, 0ULL, mod);
}

// Computes exact (a * b) mod mod using precomputed Shoup factor
// b_shoup = ComputeShoupFactor(b, mod), for mod < 2^62 and b < mod.
constexpr uint64_t ShoupMulMod(uint64_t a, uint64_t b, uint64_t b_shoup,
                               uint64_t mod) noexcept {
  const uint64_t q = detail::MulHi64(a, b_shoup);
  const uint64_t r = a * b - q * mod;
  return (r >= mod) ? (r - mod) : r;
}

// Portable, zero-UB scalar reference implementations.
namespace scalar {

inline void Gf2XorArray(std::span<uint64_t> dst, std::span<const uint64_t> a,
                        std::span<const uint64_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf2XorArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    dst[i] = a[i] ^ b[i];
  }
}

inline void Gf28AddArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                         std::span<const uint8_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf28AddArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    dst[i] = static_cast<uint8_t>(a[i] ^ b[i]);
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
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    const uint8_t v = src[i];
    dst[i] = static_cast<uint8_t>(luts.lo[v & 0x0Fu] ^ luts.hi[v >> 4]);
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
  if (scalar == 1u) {
    for (std::size_t i = 0; i < n; ++i) {
      dst[i] = static_cast<uint8_t>(dst[i] ^ src[i]);
    }
    return;
  }
  const detail::Gf28NibbleLuts luts =
      detail::BuildGf28NibbleLuts(scalar, irred_poly);
  for (std::size_t i = 0; i < n; ++i) {
    const uint8_t v = src[i];
    dst[i] =
        static_cast<uint8_t>(dst[i] ^ luts.lo[v & 0x0Fu] ^ luts.hi[v >> 4]);
  }
}

inline void Gf28MulArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                         std::span<const uint8_t> b, uint16_t irred_poly) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf28MulArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    dst[i] = Gf28MulScalar(a[i], b[i], irred_poly);
  }
}

inline void Gf128MulArray(std::span<std::pair<uint64_t, uint64_t>> dst,
                          std::span<const std::pair<uint64_t, uint64_t>> a,
                          std::span<const std::pair<uint64_t, uint64_t>> b,
                          uint64_t reduction_poly_lo) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf128MulArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    dst[i] = Gf128MulScalar(a[i], b[i], reduction_poly_lo);
  }
}

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
  for (std::size_t i = 0; i < n; ++i) {
    const uint64_t av = a[i];
    const uint64_t bv = b[i];
    const uint64_t diff = mod - bv;
    dst[i] = (av >= diff) ? (av - diff) : (av + bv);
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
  for (std::size_t i = 0; i < n; ++i) {
    const uint64_t av = a[i];
    const uint64_t bv = b[i];
    dst[i] = (av >= bv) ? (av - bv) : (av + (mod - bv));
  }
}

inline void GfPrimeMulShoupArray(std::span<uint64_t> dst,
                                 std::span<const uint64_t> a, uint64_t b,
                                 uint64_t b_shoup, uint64_t mod) {
  if (mod == 0ULL || mod >= (1ULL << 62)) {
    throw std::invalid_argument(
        "GfPrimeMulShoupArray: mod must satisfy 0 < mod < 2^62");
  }
  if (dst.size() != a.size()) {
    throw std::invalid_argument("GfPrimeMulShoupArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  for (std::size_t i = 0; i < n; ++i) {
    dst[i] = ShoupMulMod(a[i], b, b_shoup, mod);
  }
}

}  // namespace scalar

// Top-level runtime-dispatched entrypoints in xg::simd.
// Currently route all ISA levels to the portable scalar reference kernels;
// vectorized AVX2/AVX-512/NEON backends plug into these dispatch points in
// v0.5.

inline void Gf2XorArray(std::span<uint64_t> dst, std::span<const uint64_t> a,
                        std::span<const uint64_t> b) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::Gf2XorArray(dst, a, b);
      break;
  }
}

inline void Gf28AddArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                         std::span<const uint8_t> b) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::Gf28AddArray(dst, a, b);
      break;
  }
}

inline void Gf28ScalarMulArray(std::span<uint8_t> dst,
                               std::span<const uint8_t> src, uint8_t scalar,
                               uint16_t irred_poly) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::Gf28ScalarMulArray(dst, src, scalar, irred_poly);
      break;
  }
}

inline void Gf28FmaArray(std::span<uint8_t> dst, std::span<const uint8_t> src,
                         uint8_t scalar, uint16_t irred_poly) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::Gf28FmaArray(dst, src, scalar, irred_poly);
      break;
  }
}

inline void Gf28MulArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                         std::span<const uint8_t> b, uint16_t irred_poly) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::Gf28MulArray(dst, a, b, irred_poly);
      break;
  }
}

inline void Gf128MulArray(std::span<std::pair<uint64_t, uint64_t>> dst,
                          std::span<const std::pair<uint64_t, uint64_t>> a,
                          std::span<const std::pair<uint64_t, uint64_t>> b,
                          uint64_t reduction_poly_lo) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::Gf128MulArray(dst, a, b, reduction_poly_lo);
      break;
  }
}

inline void GfPrimeAddArray(std::span<uint64_t> dst,
                            std::span<const uint64_t> a,
                            std::span<const uint64_t> b, uint64_t mod) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::GfPrimeAddArray(dst, a, b, mod);
      break;
  }
}

inline void GfPrimeSubArray(std::span<uint64_t> dst,
                            std::span<const uint64_t> a,
                            std::span<const uint64_t> b, uint64_t mod) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::GfPrimeSubArray(dst, a, b, mod);
      break;
  }
}

inline void GfPrimeMulShoupArray(std::span<uint64_t> dst,
                                 std::span<const uint64_t> a, uint64_t b,
                                 uint64_t b_shoup, uint64_t mod) {
  switch (GetActiveIsaLevel()) {
    case IsaLevel::kAvx512Gfni:
    case IsaLevel::kAvx2:
    case IsaLevel::kSse42Pclmul:
    case IsaLevel::kNeonPmull:
    case IsaLevel::kNeon:
    case IsaLevel::kScalar:
    default:
      scalar::GfPrimeMulShoupArray(dst, a, b, b_shoup, mod);
      break;
  }
}

}  // namespace simd
}  // namespace xg

#endif  // XGALOIS_SIMD_KERNELS_SCALAR_HPP_
