#ifndef XGALOIS_SIMD_KERNELS_AVX2_HPP_
#define XGALOIS_SIMD_KERNELS_AVX2_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define XGALOIS_SIMD_HAS_AVX2 1
#else
#define XGALOIS_SIMD_HAS_AVX2 0
#endif

namespace xg {
namespace simd {
namespace avx2 {

#if XGALOIS_SIMD_HAS_AVX2

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

}  // namespace detail

__attribute__((target("avx2"))) inline void Gf2XorArray(
    std::span<uint64_t> dst, std::span<const uint64_t> a,
    std::span<const uint64_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf2XorArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint64_t* dst_ptr = dst.data();
  const uint64_t* a_ptr = a.data();
  const uint64_t* b_ptr = b.data();

  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const __m256i va =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a_ptr + i));
    const __m256i vb =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b_ptr + i));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i),
                        _mm256_xor_si256(va, vb));
  }
  for (; i < n; ++i) {
    dst_ptr[i] = a_ptr[i] ^ b_ptr[i];
  }
}

__attribute__((target("avx2"))) inline void Gf28AddArray(
    std::span<uint8_t> dst, std::span<const uint8_t> a,
    std::span<const uint8_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf28AddArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* a_ptr = a.data();
  const uint8_t* b_ptr = b.data();

  std::size_t i = 0;
  for (; i + 32 <= n; i += 32) {
    const __m256i va =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a_ptr + i));
    const __m256i vb =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b_ptr + i));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i),
                        _mm256_xor_si256(va, vb));
  }
  for (; i < n; ++i) {
    dst_ptr[i] = static_cast<uint8_t>(a_ptr[i] ^ b_ptr[i]);
  }
}

__attribute__((target("avx2"))) inline void Gf28ScalarMulArray(
    std::span<uint8_t> dst, std::span<const uint8_t> src, uint8_t scalar,
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
  const __m128i lo128 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(luts.lo.data()));
  const __m128i hi128 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(luts.hi.data()));
  const __m256i v_lo_lut = _mm256_broadcastsi128_si256(lo128);
  const __m256i v_hi_lut = _mm256_broadcastsi128_si256(hi128);
  const __m256i v_mask = _mm256_set1_epi8(0x0F);

  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* src_ptr = src.data();

  std::size_t i = 0;
  for (; i + 32 <= n; i += 32) {
    const __m256i v =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + i));
    const __m256i lo_nib = _mm256_and_si256(v, v_mask);
    const __m256i hi_nib = _mm256_and_si256(_mm256_srli_epi16(v, 4), v_mask);
    const __m256i res = _mm256_xor_si256(_mm256_shuffle_epi8(v_lo_lut, lo_nib),
                                         _mm256_shuffle_epi8(v_hi_lut, hi_nib));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i), res);
  }
  for (; i < n; ++i) {
    const uint8_t v = src_ptr[i];
    dst_ptr[i] = static_cast<uint8_t>(luts.lo[v & 0x0Fu] ^ luts.hi[v >> 4]);
  }
}

__attribute__((target("avx2"))) inline void Gf28FmaArray(
    std::span<uint8_t> dst, std::span<const uint8_t> src, uint8_t scalar,
    uint16_t irred_poly) {
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
    for (; i + 32 <= n; i += 32) {
      const __m256i vd =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst_ptr + i));
      const __m256i vs =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + i));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i),
                          _mm256_xor_si256(vd, vs));
    }
    for (; i < n; ++i) {
      dst_ptr[i] = static_cast<uint8_t>(dst_ptr[i] ^ src_ptr[i]);
    }
    return;
  }

  const detail::Gf28NibbleLuts luts =
      detail::BuildGf28NibbleLuts(scalar, irred_poly);
  const __m128i lo128 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(luts.lo.data()));
  const __m128i hi128 =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(luts.hi.data()));
  const __m256i v_lo_lut = _mm256_broadcastsi128_si256(lo128);
  const __m256i v_hi_lut = _mm256_broadcastsi128_si256(hi128);
  const __m256i v_mask = _mm256_set1_epi8(0x0F);

  std::size_t i = 0;
  for (; i + 32 <= n; i += 32) {
    const __m256i v =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + i));
    const __m256i vd =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst_ptr + i));
    const __m256i lo_nib = _mm256_and_si256(v, v_mask);
    const __m256i hi_nib = _mm256_and_si256(_mm256_srli_epi16(v, 4), v_mask);
    const __m256i prod =
        _mm256_xor_si256(_mm256_shuffle_epi8(v_lo_lut, lo_nib),
                         _mm256_shuffle_epi8(v_hi_lut, hi_nib));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i),
                        _mm256_xor_si256(vd, prod));
  }
  for (; i < n; ++i) {
    const uint8_t v = src_ptr[i];
    dst_ptr[i] =
        static_cast<uint8_t>(dst_ptr[i] ^ luts.lo[v & 0x0Fu] ^ luts.hi[v >> 4]);
  }
}

__attribute__((target("avx2,pclmul"))) inline void Gf128MulArray(
    std::span<std::pair<uint64_t, uint64_t>> dst,
    std::span<const std::pair<uint64_t, uint64_t>> a,
    std::span<const std::pair<uint64_t, uint64_t>> b,
    uint64_t reduction_poly_lo) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf128MulArray: span sizes must match");
  }

  const __m128i v_poly =
      _mm_cvtsi64_si128(static_cast<long long>(reduction_poly_lo));
  const std::size_t n = dst.size();

  for (std::size_t i = 0; i < n; ++i) {
    const __m128i va = _mm_set_epi64x(static_cast<long long>(a[i].second),
                                      static_cast<long long>(a[i].first));
    const __m128i vb = _mm_set_epi64x(static_cast<long long>(b[i].second),
                                      static_cast<long long>(b[i].first));

    // Karatsuba 128x128 carry-less product (3 pclmulqdq).
    const __m128i z0 = _mm_clmulepi64_si128(va, vb, 0x00);
    const __m128i z2 = _mm_clmulepi64_si128(va, vb, 0x11);
    const __m128i ax = _mm_xor_si128(va, _mm_srli_si128(va, 8));
    const __m128i bx = _mm_xor_si128(vb, _mm_srli_si128(vb, 8));
    const __m128i zm = _mm_clmulepi64_si128(ax, bx, 0x00);
    const __m128i z1 = _mm_xor_si128(_mm_xor_si128(zm, z0), z2);

    __m128i lo128 = _mm_xor_si128(z0, _mm_slli_si128(z1, 8));
    __m128i hi128 = _mm_xor_si128(z2, _mm_srli_si128(z1, 8));

    // Reduce c3 * x^192 = (c3 * R) * x^64, then c2' * x^128 = c2' * R.
    const __m128i t3 = _mm_clmulepi64_si128(hi128, v_poly, 0x01);
    lo128 = _mm_xor_si128(lo128, _mm_slli_si128(t3, 8));
    hi128 = _mm_xor_si128(hi128, _mm_srli_si128(t3, 8));

    const __m128i t2 = _mm_clmulepi64_si128(hi128, v_poly, 0x00);
    const __m128i res = _mm_xor_si128(lo128, t2);

    const uint64_t r0 = static_cast<uint64_t>(_mm_cvtsi128_si64(res));
    const uint64_t r1 =
        static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(res, 8)));
    dst[i] = {r0, r1};
  }
}

__attribute__((target("avx2"))) inline void GfPrimeAddArray(
    std::span<uint64_t> dst, std::span<const uint64_t> a,
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

  const __m256i v_mod = _mm256_set1_epi64x(static_cast<long long>(mod));
  const __m256i v_sign =
      _mm256_set1_epi64x(static_cast<long long>(0x8000000000000000ULL));

  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const __m256i va =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a_ptr + i));
    const __m256i vb =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b_ptr + i));
    const __m256i v_diff = _mm256_sub_epi64(v_mod, vb);
    // Unsigned (diff > a) <=> unsigned (a < mod - b), using sign-bit flip.
    const __m256i v_lt = _mm256_cmpgt_epi64(_mm256_xor_si256(v_diff, v_sign),
                                            _mm256_xor_si256(va, v_sign));
    const __m256i v_sum = _mm256_add_epi64(va, vb);
    const __m256i v_sub_mod = _mm256_andnot_si256(v_lt, v_mod);
    const __m256i v_res = _mm256_sub_epi64(v_sum, v_sub_mod);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i), v_res);
  }
  for (; i < n; ++i) {
    const uint64_t av = a_ptr[i];
    const uint64_t bv = b_ptr[i];
    const uint64_t diff = mod - bv;
    dst_ptr[i] = (av >= diff) ? (av - diff) : (av + bv);
  }
}

__attribute__((target("avx2"))) inline void GfPrimeSubArray(
    std::span<uint64_t> dst, std::span<const uint64_t> a,
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

  const __m256i v_mod = _mm256_set1_epi64x(static_cast<long long>(mod));
  const __m256i v_sign =
      _mm256_set1_epi64x(static_cast<long long>(0x8000000000000000ULL));

  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const __m256i va =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a_ptr + i));
    const __m256i vb =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b_ptr + i));
    // Unsigned (b > a) <=> unsigned (a < b), using sign-bit flip.
    const __m256i v_lt = _mm256_cmpgt_epi64(_mm256_xor_si256(vb, v_sign),
                                            _mm256_xor_si256(va, v_sign));
    const __m256i v_diff = _mm256_sub_epi64(va, vb);
    const __m256i v_add_mod = _mm256_and_si256(v_lt, v_mod);
    const __m256i v_res = _mm256_add_epi64(v_diff, v_add_mod);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i), v_res);
  }
  for (; i < n; ++i) {
    const uint64_t av = a_ptr[i];
    const uint64_t bv = b_ptr[i];
    dst_ptr[i] = (av >= bv) ? (av - bv) : (av + (mod - bv));
  }
}

#endif  // XGALOIS_SIMD_HAS_AVX2

}  // namespace avx2
}  // namespace simd
}  // namespace xg

#endif  // XGALOIS_SIMD_KERNELS_AVX2_HPP_
