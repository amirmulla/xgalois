#ifndef XGALOIS_SIMD_KERNELS_AVX512_HPP_
#define XGALOIS_SIMD_KERNELS_AVX512_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define XGALOIS_SIMD_HAS_AVX512 1
#else
#define XGALOIS_SIMD_HAS_AVX512 0
#endif

namespace xg {
namespace simd {
namespace avx512 {

namespace detail {

constexpr uint8_t Gf28MulByX(uint8_t value, uint8_t poly_lo) noexcept {
  const uint8_t mask =
      static_cast<uint8_t>(0u - static_cast<unsigned>(value >> 7));
  return static_cast<uint8_t>((static_cast<unsigned>(value) << 1) ^
                              (poly_lo & mask));
}

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
    lhs = Gf28MulByX(lhs, poly_lo);
    rhs = static_cast<uint8_t>(rhs >> 1);
  }
  return result;
}

// Precomputes the 8x8 GF(2) multiplication matrix M(scalar, irred_poly) packed
// into a uint64_t for VGF2P8AFFINEQB (_mm512_gf2p8affine_epi64_epi8(..., 0)).
//
// Intel VGF2P8AFFINEQB computes for each byte x and output bit i in 0..7:
//   out.bit[i] = parity(M.byte[7 - i] & x)
// Since (x * scalar) mod poly = XOR_{k=0..7} x.bit[k] * (scalar * x^k mod
// poly), we set M.byte[7 - i].bit[k] = (scalar * x^k mod poly).bit[i].
constexpr uint64_t BuildGfniMulMatrix(uint8_t scalar,
                                      uint16_t irred_poly) noexcept {
  const uint8_t poly_lo = static_cast<uint8_t>(irred_poly & 0xFFu);
  uint8_t cols[8] = {};
  cols[0] = scalar;
  for (int k = 1; k < 8; ++k) {
    cols[k] = Gf28MulByX(cols[k - 1], poly_lo);
  }

  uint64_t matrix = 0ULL;
  for (int i = 0; i < 8; ++i) {
    uint8_t row = 0u;
    for (int k = 0; k < 8; ++k) {
      const uint8_t bit = static_cast<uint8_t>((cols[k] >> i) & 1u);
      row |= static_cast<uint8_t>(bit << k);
    }
    matrix |= static_cast<uint64_t>(row) << (8 * (7 - i));
  }
  return matrix;
}

}  // namespace detail

#if XGALOIS_SIMD_HAS_AVX512

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
Gf2XorArray(std::span<uint64_t> dst, std::span<const uint64_t> a,
            std::span<const uint64_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf2XorArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint64_t* dst_ptr = dst.data();
  const uint64_t* a_ptr = a.data();
  const uint64_t* b_ptr = b.data();

  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    const __m512i va = _mm512_loadu_si512(a_ptr + i);
    const __m512i vb = _mm512_loadu_si512(b_ptr + i);
    _mm512_storeu_si512(dst_ptr + i, _mm512_xor_si512(va, vb));
  }
  for (; i < n; ++i) {
    dst_ptr[i] = a_ptr[i] ^ b_ptr[i];
  }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
Gf28AddArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
             std::span<const uint8_t> b) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf28AddArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* a_ptr = a.data();
  const uint8_t* b_ptr = b.data();

  std::size_t i = 0;
  for (; i + 64 <= n; i += 64) {
    const __m512i va = _mm512_loadu_si512(a_ptr + i);
    const __m512i vb = _mm512_loadu_si512(b_ptr + i);
    _mm512_storeu_si512(dst_ptr + i, _mm512_xor_si512(va, vb));
  }
  for (; i < n; ++i) {
    dst_ptr[i] = static_cast<uint8_t>(a_ptr[i] ^ b_ptr[i]);
  }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
Gf28ScalarMulArray(std::span<uint8_t> dst, std::span<const uint8_t> src,
                   uint8_t scalar, uint16_t irred_poly) {
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

  const uint64_t matrix = detail::BuildGfniMulMatrix(scalar, irred_poly);
  const __m512i v_mat = _mm512_set1_epi64(static_cast<long long>(matrix));

  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* src_ptr = src.data();

  std::size_t i = 0;
  for (; i + 64 <= n; i += 64) {
    const __m512i v_src = _mm512_loadu_si512(src_ptr + i);
    const __m512i v_res = _mm512_gf2p8affine_epi64_epi8(v_src, v_mat, 0);
    _mm512_storeu_si512(dst_ptr + i, v_res);
  }
  if (i < n) {
    const __mmask64 tail_mask = (uint64_t{1} << (n - i)) - 1ULL;
    const __m512i v_src = _mm512_maskz_loadu_epi8(tail_mask, src_ptr + i);
    const __m512i v_res = _mm512_gf2p8affine_epi64_epi8(v_src, v_mat, 0);
    _mm512_mask_storeu_epi8(dst_ptr + i, tail_mask, v_res);
  }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
Gf28FmaArray(std::span<uint8_t> dst, std::span<const uint8_t> src,
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
    for (; i + 64 <= n; i += 64) {
      const __m512i vd = _mm512_loadu_si512(dst_ptr + i);
      const __m512i vs = _mm512_loadu_si512(src_ptr + i);
      _mm512_storeu_si512(dst_ptr + i, _mm512_xor_si512(vd, vs));
    }
    for (; i < n; ++i) {
      dst_ptr[i] = static_cast<uint8_t>(dst_ptr[i] ^ src_ptr[i]);
    }
    return;
  }

  const uint64_t matrix = detail::BuildGfniMulMatrix(scalar, irred_poly);
  const __m512i v_mat = _mm512_set1_epi64(static_cast<long long>(matrix));

  std::size_t i = 0;
  for (; i + 64 <= n; i += 64) {
    const __m512i v_src = _mm512_loadu_si512(src_ptr + i);
    const __m512i v_dst = _mm512_loadu_si512(dst_ptr + i);
    const __m512i v_prod = _mm512_gf2p8affine_epi64_epi8(v_src, v_mat, 0);
    _mm512_storeu_si512(dst_ptr + i, _mm512_xor_si512(v_dst, v_prod));
  }
  if (i < n) {
    const __mmask64 tail_mask = (uint64_t{1} << (n - i)) - 1ULL;
    const __m512i v_src = _mm512_maskz_loadu_epi8(tail_mask, src_ptr + i);
    const __m512i v_dst = _mm512_maskz_loadu_epi8(tail_mask, dst_ptr + i);
    const __m512i v_prod = _mm512_gf2p8affine_epi64_epi8(v_src, v_mat, 0);
    _mm512_mask_storeu_epi8(dst_ptr + i, tail_mask,
                            _mm512_xor_si512(v_dst, v_prod));
  }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
Gf28MulArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
             std::span<const uint8_t> b, uint16_t irred_poly) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf28MulArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  uint8_t* dst_ptr = dst.data();
  const uint8_t* a_ptr = a.data();
  const uint8_t* b_ptr = b.data();

  if (irred_poly == 0x11Bu) {
    std::size_t i = 0;
    for (; i + 64 <= n; i += 64) {
      const __m512i va = _mm512_loadu_si512(a_ptr + i);
      const __m512i vb = _mm512_loadu_si512(b_ptr + i);
      _mm512_storeu_si512(dst_ptr + i, _mm512_gf2p8mul_epi8(va, vb));
    }
    if (i < n) {
      const __mmask64 tail_mask = (uint64_t{1} << (n - i)) - 1ULL;
      const __m512i va = _mm512_maskz_loadu_epi8(tail_mask, a_ptr + i);
      const __m512i vb = _mm512_maskz_loadu_epi8(tail_mask, b_ptr + i);
      _mm512_mask_storeu_epi8(dst_ptr + i, tail_mask,
                              _mm512_gf2p8mul_epi8(va, vb));
    }
    return;
  }

  for (std::size_t i = 0; i < n; ++i) {
    dst_ptr[i] = detail::Gf28MulScalar(a_ptr[i], b_ptr[i], irred_poly);
  }
}

__attribute__((
    target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline __m512i
Gf128Mul4(__m512i va, __m512i vb, __m512i v_poly) noexcept {
  const __m512i z0 = _mm512_clmulepi64_epi128(va, vb, 0x00);
  const __m512i z2 = _mm512_clmulepi64_epi128(va, vb, 0x11);
  const __m512i ax = _mm512_xor_si512(va, _mm512_bsrli_epi128(va, 8));
  const __m512i bx = _mm512_xor_si512(vb, _mm512_bsrli_epi128(vb, 8));
  const __m512i zm = _mm512_clmulepi64_epi128(ax, bx, 0x00);
  const __m512i z1 = _mm512_ternarylogic_epi64(zm, z0, z2, 0x96);

  __m512i lo128 = _mm512_xor_si512(z0, _mm512_bslli_epi128(z1, 8));
  __m512i hi128 = _mm512_xor_si512(z2, _mm512_bsrli_epi128(z1, 8));

  const __m512i t3 = _mm512_clmulepi64_epi128(hi128, v_poly, 0x01);
  lo128 = _mm512_xor_si512(lo128, _mm512_bslli_epi128(t3, 8));
  hi128 = _mm512_xor_si512(hi128, _mm512_bsrli_epi128(t3, 8));

  const __m512i t2 = _mm512_clmulepi64_epi128(hi128, v_poly, 0x00);
  return _mm512_xor_si512(lo128, t2);
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
Gf128MulArray(std::span<std::pair<uint64_t, uint64_t>> dst,
              std::span<const std::pair<uint64_t, uint64_t>> a,
              std::span<const std::pair<uint64_t, uint64_t>> b,
              uint64_t reduction_poly_lo) {
  if (dst.size() != a.size() || dst.size() != b.size()) {
    throw std::invalid_argument("Gf128MulArray: span sizes must match");
  }
  const std::size_t n = dst.size();
  if (n == 0) {
    return;
  }

  const long long poly_ll = static_cast<long long>(reduction_poly_lo);
  const __m512i v_poly =
      _mm512_set_epi64(0, poly_ll, 0, poly_ll, 0, poly_ll, 0, poly_ll);

  alignas(64) uint64_t a_buf[8];
  alignas(64) uint64_t b_buf[8];
  alignas(64) uint64_t dst_buf[8];

  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    for (std::size_t k = 0; k < 4; ++k) {
      a_buf[2 * k] = a[i + k].first;
      a_buf[2 * k + 1] = a[i + k].second;
      b_buf[2 * k] = b[i + k].first;
      b_buf[2 * k + 1] = b[i + k].second;
    }
    const __m512i va = _mm512_load_si512(a_buf);
    const __m512i vb = _mm512_load_si512(b_buf);
    _mm512_store_si512(dst_buf, Gf128Mul4(va, vb, v_poly));
    for (std::size_t k = 0; k < 4; ++k) {
      dst[i + k] = {dst_buf[2 * k], dst_buf[2 * k + 1]};
    }
  }

  if (i < n) {
    const std::size_t rem = n - i;
    for (std::size_t k = 0; k < 8; ++k) {
      a_buf[k] = 0ULL;
      b_buf[k] = 0ULL;
    }
    for (std::size_t k = 0; k < rem; ++k) {
      a_buf[2 * k] = a[i + k].first;
      a_buf[2 * k + 1] = a[i + k].second;
      b_buf[2 * k] = b[i + k].first;
      b_buf[2 * k + 1] = b[i + k].second;
    }
    const __m512i va = _mm512_load_si512(a_buf);
    const __m512i vb = _mm512_load_si512(b_buf);
    _mm512_store_si512(dst_buf, Gf128Mul4(va, vb, v_poly));
    for (std::size_t k = 0; k < rem; ++k) {
      dst[i + k] = {dst_buf[2 * k], dst_buf[2 * k + 1]};
    }
  }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
GfPrimeAddArray(std::span<uint64_t> dst, std::span<const uint64_t> a,
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

  const __m512i v_mod = _mm512_set1_epi64(static_cast<long long>(mod));

  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    const __m512i va = _mm512_loadu_si512(a_ptr + i);
    const __m512i vb = _mm512_loadu_si512(b_ptr + i);
    const __m512i v_diff = _mm512_sub_epi64(v_mod, vb);
    const __mmask8 ge_mask = _mm512_cmpge_epu64_mask(va, v_diff);
    const __m512i v_sum = _mm512_add_epi64(va, vb);
    const __m512i v_res = _mm512_mask_sub_epi64(v_sum, ge_mask, v_sum, v_mod);
    _mm512_storeu_si512(dst_ptr + i, v_res);
  }
  for (; i < n; ++i) {
    const uint64_t av = a_ptr[i];
    const uint64_t bv = b_ptr[i];
    const uint64_t diff = mod - bv;
    dst_ptr[i] = (av >= diff) ? (av - diff) : (av + bv);
  }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni,vpclmulqdq"))) inline void
GfPrimeSubArray(std::span<uint64_t> dst, std::span<const uint64_t> a,
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

  const __m512i v_mod = _mm512_set1_epi64(static_cast<long long>(mod));

  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    const __m512i va = _mm512_loadu_si512(a_ptr + i);
    const __m512i vb = _mm512_loadu_si512(b_ptr + i);
    const __mmask8 lt_mask = _mm512_cmplt_epu64_mask(va, vb);
    const __m512i v_diff = _mm512_sub_epi64(va, vb);
    const __m512i v_res = _mm512_mask_add_epi64(v_diff, lt_mask, v_diff, v_mod);
    _mm512_storeu_si512(dst_ptr + i, v_res);
  }
  for (; i < n; ++i) {
    const uint64_t av = a_ptr[i];
    const uint64_t bv = b_ptr[i];
    dst_ptr[i] = (av >= bv) ? (av - bv) : (av + (mod - bv));
  }
}

#endif  // XGALOIS_SIMD_HAS_AVX512

}  // namespace avx512
}  // namespace simd
}  // namespace xg

#endif  // XGALOIS_SIMD_KERNELS_AVX512_HPP_
