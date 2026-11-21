#ifndef XGALOIS_FIELD_GF128_HPP_
#define XGALOIS_FIELD_GF128_HPP_

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cctype>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__PCLMUL__) && defined(__x86_64__)
#include <emmintrin.h>
#include <wmmintrin.h>
#endif

#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/utils/field.hpp"

namespace xg {

namespace detail {

class StreamFlagsGuard128 {
 public:
  explicit StreamFlagsGuard128(std::ostream& os)
      : os_(os), flags_(os.flags()) {}
  ~StreamFlagsGuard128() { os_.flags(flags_); }

  StreamFlagsGuard128(const StreamFlagsGuard128&) = delete;
  StreamFlagsGuard128& operator=(const StreamFlagsGuard128&) = delete;

 private:
  std::ostream& os_;
  std::ios_base::fmtflags flags_;
};

inline std::string_view TrimAscii128(std::string_view sv) noexcept {
  size_t start = 0;
  while (start < sv.size() &&
         std::isspace(static_cast<unsigned char>(sv[start])) != 0) {
    ++start;
  }
  size_t end = sv.size();
  while (end > start &&
         std::isspace(static_cast<unsigned char>(sv[end - 1])) != 0) {
    --end;
  }
  return sv.substr(start, end - start);
}

constexpr uint64_t BitReverse64(uint64_t x) noexcept {
  x = ((x & 0x5555555555555555ULL) << 1) | ((x >> 1) & 0x5555555555555555ULL);
  x = ((x & 0x3333333333333333ULL) << 2) | ((x >> 2) & 0x3333333333333333ULL);
  x = ((x & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
  x = ((x & 0x00FF00FF00FF00FFULL) << 8) | ((x >> 8) & 0x00FF00FF00FF00FFULL);
  x = ((x & 0x0000FFFF0000FFFFULL) << 16) | ((x >> 16) & 0x0000FFFF0000FFFFULL);
  return (x << 32) | (x >> 32);
}

}  // namespace detail

// Portable, trivially copyable 128-bit unsigned value type for GF(2^128).
// Bit 0..63 reside in `lo`, bits 64..127 reside in `hi`.
struct UInt128 {
  uint64_t lo = 0;
  uint64_t hi = 0;

  constexpr bool IsZero() const noexcept { return lo == 0 && hi == 0; }
  constexpr bool IsOne() const noexcept { return lo == 1 && hi == 0; }

  constexpr bool Bit(int i) const noexcept {
    if (i < 0 || i >= 128) {
      return false;
    }
    if (i < 64) {
      return ((lo >> i) & 1ULL) != 0;
    }
    return ((hi >> (i - 64)) & 1ULL) != 0;
  }

  constexpr int BitWidth() const noexcept {
    if (hi != 0) {
      return 64 + static_cast<int>(std::bit_width(hi));
    }
    return static_cast<int>(std::bit_width(lo));
  }

  constexpr int PopCount() const noexcept {
    return std::popcount(lo) + std::popcount(hi);
  }

  constexpr UInt128 BitReverse() const noexcept {
    return UInt128{detail::BitReverse64(hi), detail::BitReverse64(lo)};
  }

  constexpr UInt128 operator^(const UInt128& other) const noexcept {
    return UInt128{lo ^ other.lo, hi ^ other.hi};
  }

  constexpr UInt128& operator^=(const UInt128& other) noexcept {
    lo ^= other.lo;
    hi ^= other.hi;
    return *this;
  }

  constexpr UInt128 operator&(const UInt128& other) const noexcept {
    return UInt128{lo & other.lo, hi & other.hi};
  }

  constexpr UInt128& operator&=(const UInt128& other) noexcept {
    lo &= other.lo;
    hi &= other.hi;
    return *this;
  }

  constexpr UInt128 operator|(const UInt128& other) const noexcept {
    return UInt128{lo | other.lo, hi | other.hi};
  }

  constexpr UInt128& operator|=(const UInt128& other) noexcept {
    lo |= other.lo;
    hi |= other.hi;
    return *this;
  }

  constexpr UInt128 operator~() const noexcept { return UInt128{~lo, ~hi}; }

  constexpr UInt128 operator<<(int shift) const noexcept {
    if (shift <= 0) {
      return *this;
    }
    if (shift >= 128) {
      return UInt128{0, 0};
    }
    if (shift >= 64) {
      return UInt128{0, lo << (shift - 64)};
    }
    return UInt128{lo << shift, (hi << shift) | (lo >> (64 - shift))};
  }

  constexpr UInt128& operator<<=(int shift) noexcept {
    *this = *this << shift;
    return *this;
  }

  constexpr UInt128 operator>>(int shift) const noexcept {
    if (shift <= 0) {
      return *this;
    }
    if (shift >= 128) {
      return UInt128{0, 0};
    }
    if (shift >= 64) {
      return UInt128{hi >> (shift - 64), 0};
    }
    return UInt128{(lo >> shift) | (hi << (64 - shift)), hi >> shift};
  }

  constexpr UInt128& operator>>=(int shift) noexcept {
    *this = *this >> shift;
    return *this;
  }

  constexpr bool operator==(const UInt128& other) const noexcept {
    return lo == other.lo && hi == other.hi;
  }

  constexpr std::strong_ordering operator<=>(
      const UInt128& other) const noexcept {
    if (hi != other.hi) {
      return hi <=> other.hi;
    }
    return lo <=> other.lo;
  }

  static UInt128 FromHex(std::string_view sv) {
    const std::string_view trimmed = detail::TrimAscii128(sv);
    if (trimmed.empty()) {
      throw std::invalid_argument("UInt128::FromHex: empty string");
    }
    std::string_view digits = trimmed;
    if (digits.size() >= 2 && digits[0] == '0' &&
        (digits[1] == 'x' || digits[1] == 'X')) {
      digits.remove_prefix(2);
    }
    if (digits.empty()) {
      throw std::invalid_argument("UInt128::FromHex: missing hex digits");
    }

    UInt128 value{0, 0};
    for (char ch : digits) {
      uint64_t nibble = 0;
      if (ch >= '0' && ch <= '9') {
        nibble = static_cast<uint64_t>(ch - '0');
      } else if (ch >= 'a' && ch <= 'f') {
        nibble = static_cast<uint64_t>(10 + (ch - 'a'));
      } else if (ch >= 'A' && ch <= 'F') {
        nibble = static_cast<uint64_t>(10 + (ch - 'A'));
      } else {
        throw std::invalid_argument("UInt128::FromHex: invalid hex digit");
      }
      if ((value.hi >> 60) != 0) {
        throw std::out_of_range("UInt128::FromHex: value exceeds 128 bits");
      }
      value <<= 4;
      value.lo |= nibble;
    }
    return value;
  }

  std::string ToHex(bool zero_pad = false) const {
    static constexpr char kHexDigits[] = "0123456789abcdef";
    std::string out;
    out.reserve(34);
    out.append("0x");
    if (zero_pad) {
      for (int i = 60; i >= 0; i -= 4) {
        out.push_back(kHexDigits[(hi >> i) & 0xFULL]);
      }
      for (int i = 60; i >= 0; i -= 4) {
        out.push_back(kHexDigits[(lo >> i) & 0xFULL]);
      }
      return out;
    }
    if (IsZero()) {
      out.push_back('0');
      return out;
    }
    if (hi != 0) {
      bool started = false;
      for (int i = 60; i >= 0; i -= 4) {
        const unsigned int nib = static_cast<unsigned int>((hi >> i) & 0xFULL);
        if (nib != 0 || started) {
          out.push_back(kHexDigits[nib]);
          started = true;
        }
      }
      for (int i = 60; i >= 0; i -= 4) {
        out.push_back(kHexDigits[(lo >> i) & 0xFULL]);
      }
      return out;
    }
    bool started = false;
    for (int i = 60; i >= 0; i -= 4) {
      const unsigned int nib = static_cast<unsigned int>((lo >> i) & 0xFULL);
      if (nib != 0 || started) {
        out.push_back(kHexDigits[nib]);
        started = true;
      }
    }
    return out;
  }

  static UInt128 FromDecimal(std::string_view sv) {
    const std::string_view trimmed = detail::TrimAscii128(sv);
    if (trimmed.empty()) {
      throw std::invalid_argument("UInt128::FromDecimal: empty string");
    }

    auto add_checked = [](UInt128& acc, UInt128 rhs) noexcept -> bool {
      const uint64_t new_lo = acc.lo + rhs.lo;
      const uint64_t carry = (new_lo < acc.lo) ? 1ULL : 0ULL;
      const uint64_t new_hi = acc.hi + rhs.hi;
      if (new_hi < acc.hi) {
        return false;
      }
      const uint64_t final_hi = new_hi + carry;
      if (final_hi < new_hi) {
        return false;
      }
      acc.lo = new_lo;
      acc.hi = final_hi;
      return true;
    };

    UInt128 value{0, 0};
    for (char ch : trimmed) {
      if (ch < '0' || ch > '9') {
        throw std::invalid_argument("UInt128::FromDecimal: invalid digit");
      }
      if ((value.hi >> 61) != 0) {
        throw std::out_of_range("UInt128::FromDecimal: value exceeds 128 bits");
      }
      UInt128 times8 = value << 3;
      const UInt128 times2 = value << 1;
      if (!add_checked(times8, times2) ||
          !add_checked(times8, UInt128{static_cast<uint64_t>(ch - '0'), 0})) {
        throw std::out_of_range("UInt128::FromDecimal: value exceeds 128 bits");
      }
      value = times8;
    }
    return value;
  }

  std::string ToDecimal() const {
    if (IsZero()) {
      return "0";
    }
    static constexpr uint64_t kChunkBase = 1000000000ULL;  // 10^9
    std::array<uint32_t, 4> limbs = {
        static_cast<uint32_t>(hi >> 32),
        static_cast<uint32_t>(hi & 0xFFFFFFFFULL),
        static_cast<uint32_t>(lo >> 32),
        static_cast<uint32_t>(lo & 0xFFFFFFFFULL),
    };

    std::vector<uint32_t> chunks;
    chunks.reserve(5);
    while ((limbs[0] | limbs[1] | limbs[2] | limbs[3]) != 0) {
      uint64_t rem = 0;
      for (uint32_t& limb : limbs) {
        const uint64_t cur = (rem << 32) | static_cast<uint64_t>(limb);
        limb = static_cast<uint32_t>(cur / kChunkBase);
        rem = cur % kChunkBase;
      }
      chunks.push_back(static_cast<uint32_t>(rem));
    }

    std::ostringstream oss;
    oss << chunks.back();
    for (size_t i = chunks.size() - 1; i > 0; --i) {
      oss << std::setw(9) << std::setfill('0') << chunks[i - 1];
    }
    return oss.str();
  }

  constexpr std::array<uint8_t, 16> ToBytesBigEndian() const noexcept {
    std::array<uint8_t, 16> bytes{};
    for (size_t i = 0; i < 8; ++i) {
      bytes[i] = static_cast<uint8_t>((hi >> (56 - 8 * i)) & 0xFFULL);
      bytes[8 + i] = static_cast<uint8_t>((lo >> (56 - 8 * i)) & 0xFFULL);
    }
    return bytes;
  }

  static constexpr UInt128 FromBytesBigEndian(
      const std::array<uint8_t, 16>& bytes) noexcept {
    uint64_t high = 0;
    uint64_t low = 0;
    for (size_t i = 0; i < 8; ++i) {
      high = (high << 8) | static_cast<uint64_t>(bytes[i]);
      low = (low << 8) | static_cast<uint64_t>(bytes[8 + i]);
    }
    return UInt128{low, high};
  }

  static UInt128 FromBytesBigEndian(std::span<const uint8_t> bytes) {
    if (bytes.size() != 16) {
      throw std::invalid_argument(
          "UInt128::FromBytesBigEndian: expected exactly 16 bytes");
    }
    std::array<uint8_t, 16> arr{};
    std::copy_n(bytes.data(), 16, arr.data());
    return FromBytesBigEndian(arr);
  }
};

inline std::ostream& operator<<(std::ostream& os, const UInt128& value) {
  os << value.ToHex();
  return os;
}

inline std::istream& operator>>(std::istream& is, UInt128& value) {
  std::string token;
  if (is >> token) {
    try {
      if (token.size() >= 2 && token[0] == '0' &&
          (token[1] == 'x' || token[1] == 'X')) {
        value = UInt128::FromHex(token);
      } else {
        value = UInt128::FromDecimal(token);
      }
    } catch (...) {
      is.setstate(std::ios::failbit);
    }
  }
  return is;
}

namespace detail {

// Portable 64x64 -> 128-bit carry-less multiplication over GF(2).
inline UInt128 ClMul64Portable128(uint64_t a, uint64_t b) noexcept {
  uint64_t lo = 0;
  uint64_t hi = 0;
  uint64_t tmp = b;
  while (tmp != 0) {
    const int bit = std::countr_zero(tmp);
    lo ^= (a << bit);
    if (bit > 0) {
      hi ^= (a >> (64 - bit));
    }
    tmp &= (tmp - 1);
  }
  return UInt128{lo, hi};
}

// 128x128 -> 256-bit carry-less multiplication over GF(2).
// Returns {hi_128, lo_128} representing lo_128 + hi_128 * x^128.
// Uses Karatsuba 3-PCLMULQDQ when __PCLMUL__ is enabled on x86-64, and a
// bit-identical 3-mul Karatsuba shift-XOR fallback otherwise.
inline std::pair<UInt128, UInt128> ClMul128(UInt128 a, UInt128 b) noexcept {
#if defined(__PCLMUL__) && defined(__x86_64__)
  const __m128i va = _mm_set_epi64x(static_cast<long long>(a.hi),
                                    static_cast<long long>(a.lo));
  const __m128i vb = _mm_set_epi64x(static_cast<long long>(b.hi),
                                    static_cast<long long>(b.lo));
  const __m128i lo_prod = _mm_clmulepi64_si128(va, vb, 0x00);
  const __m128i hi_prod = _mm_clmulepi64_si128(va, vb, 0x11);
  const __m128i a_xor = _mm_xor_si128(va, _mm_srli_si128(va, 8));
  const __m128i b_xor = _mm_xor_si128(vb, _mm_srli_si128(vb, 8));
  const __m128i mid_prod = _mm_clmulepi64_si128(a_xor, b_xor, 0x00);
  const __m128i mid = _mm_xor_si128(mid_prod, _mm_xor_si128(lo_prod, hi_prod));
  const __m128i res_lo = _mm_xor_si128(lo_prod, _mm_slli_si128(mid, 8));
  const __m128i res_hi = _mm_xor_si128(hi_prod, _mm_srli_si128(mid, 8));

  UInt128 lo_out{
      static_cast<uint64_t>(_mm_cvtsi128_si64(res_lo)),
      static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(res_lo, 8))),
  };
  UInt128 hi_out{
      static_cast<uint64_t>(_mm_cvtsi128_si64(res_hi)),
      static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(res_hi, 8))),
  };
  return {hi_out, lo_out};
#else
  const UInt128 lo_prod = ClMul64Portable128(a.lo, b.lo);
  const UInt128 hi_prod = ClMul64Portable128(a.hi, b.hi);
  const UInt128 mid_prod = ClMul64Portable128(a.lo ^ a.hi, b.lo ^ b.hi);
  const UInt128 mid = mid_prod ^ lo_prod ^ hi_prod;
  const UInt128 lo_out{lo_prod.lo, lo_prod.hi ^ mid.lo};
  const UInt128 hi_out{hi_prod.lo ^ mid.hi, hi_prod.hi};
  return {hi_out, lo_out};
#endif
}

// Fast 256 -> 128-bit reduction modulo AES-GCM/GHASH polynomial
// P(x) = x^128 + x^7 + x^2 + x + 1 (modulus_low = 0x87).
inline UInt128 ReduceGhashPoly128(UInt128 hi, UInt128 lo) noexcept {
  const uint64_t overflow = (hi.hi >> 57) ^ (hi.hi >> 62) ^ (hi.hi >> 63);
  const UInt128 h{hi.lo ^ overflow, hi.hi};
  const UInt128 red = h ^ (h << 1) ^ (h << 2) ^ (h << 7);
  return lo ^ red;
}

// Specialized fast multiply modulo x^128 + x^7 + x^2 + x + 1.
// Keeps the entire Karatsuba + 2-shift reduction in XMM registers on PCLMUL.
inline UInt128 MulGhashPoly128(UInt128 a, UInt128 b) noexcept {
#if defined(__PCLMUL__) && defined(__x86_64__)
  const __m128i va = _mm_set_epi64x(static_cast<long long>(a.hi),
                                    static_cast<long long>(a.lo));
  const __m128i vb = _mm_set_epi64x(static_cast<long long>(b.hi),
                                    static_cast<long long>(b.lo));
  const __m128i lo_prod = _mm_clmulepi64_si128(va, vb, 0x00);
  const __m128i hi_prod = _mm_clmulepi64_si128(va, vb, 0x11);
  const __m128i a_xor = _mm_xor_si128(va, _mm_srli_si128(va, 8));
  const __m128i b_xor = _mm_xor_si128(vb, _mm_srli_si128(vb, 8));
  const __m128i mid_prod = _mm_clmulepi64_si128(a_xor, b_xor, 0x00);
  const __m128i mid = _mm_xor_si128(mid_prod, _mm_xor_si128(lo_prod, hi_prod));
  const __m128i res_lo = _mm_xor_si128(lo_prod, _mm_slli_si128(mid, 8));
  const __m128i res_hi = _mm_xor_si128(hi_prod, _mm_srli_si128(mid, 8));

  // Two-step shift-XOR reduction modulo x^128 + x^7 + x^2 + x + 1.
  const __m128i c = _mm_xor_si128(
      _mm_srli_epi64(res_hi, 63),
      _mm_xor_si128(_mm_srli_epi64(res_hi, 62), _mm_srli_epi64(res_hi, 57)));
  const __m128i h = _mm_xor_si128(res_hi, _mm_srli_si128(c, 8));
  const __m128i shifted =
      _mm_xor_si128(_mm_slli_epi64(h, 1),
                    _mm_xor_si128(_mm_slli_epi64(h, 2), _mm_slli_epi64(h, 7)));
  const __m128i red =
      _mm_xor_si128(h, _mm_xor_si128(shifted, _mm_slli_si128(c, 8)));
  const __m128i out = _mm_xor_si128(res_lo, red);

  return UInt128{
      static_cast<uint64_t>(_mm_cvtsi128_si64(out)),
      static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(out, 8))),
  };
#else
  const auto [hi, lo] = ClMul128(a, b);
  return ReduceGhashPoly128(hi, lo);
#endif
}

// Computes the low 128 bits of mu(x) = floor(x^256 / P(x)) for a monic
// degree-128 binary polynomial P(x) = x^128 + modulus_low.
constexpr UInt128 ComputeBarrettMuLow128(UInt128 modulus_low) noexcept {
  UInt128 rem = modulus_low;
  UInt128 q{0, 0};
  for (int i = 127; i >= 0; --i) {
    const bool top = rem.Bit(127);
    rem <<= 1;
    if (top) {
      rem ^= modulus_low;
      q |= (UInt128{1, 0} << i);
    }
  }
  return q;
}

// Exact carry-less Barrett reduction of 256-bit (hi, lo) modulo
// P(x) = x^128 + modulus_low for any degree-128 binary polynomial.
inline UInt128 ReduceBarrett128(UInt128 hi, UInt128 lo, UInt128 modulus_low,
                                UInt128 barrett_mu_low) noexcept {
  if (hi.IsZero()) {
    return lo;
  }
  const auto [t1_hi, t1_lo] = ClMul128(hi, barrett_mu_low);
  (void)t1_lo;
  const UInt128 q = hi ^ t1_hi;
  const auto [t2_hi, t2_lo] = ClMul128(q, modulus_low);
  (void)t2_hi;
  return lo ^ t2_lo;
}

// Exact GCD of P(x) = x^128 + modulus_low and g(x) (deg(g) < 128).
inline UInt128 BinaryPoly128GcdWithModulus(UInt128 modulus_low,
                                           UInt128 g) noexcept {
  if (g.IsZero()) {
    return UInt128{2, 0};  // Non-1 sentinel since deg(P) == 128.
  }
  if (g.IsOne()) {
    return UInt128{1, 0};
  }

  const int deg_g = g.BitWidth() - 1;  // 1 <= deg_g <= 127
  const int shift0 = 128 - deg_g;      // 1 <= shift0 <= 127
  const UInt128 g_no_msb = g ^ (UInt128{1, 0} << deg_g);
  UInt128 r = modulus_low ^ (g_no_msb << shift0);

  while (!r.IsZero()) {
    const int deg_r = r.BitWidth() - 1;
    if (deg_r < deg_g) {
      break;
    }
    r ^= (g << (deg_r - deg_g));
  }

  UInt128 r0 = g;
  UInt128 r1 = r;
  while (!r1.IsZero()) {
    const int deg_r1 = r1.BitWidth() - 1;
    while (!r0.IsZero()) {
      const int deg_r0 = r0.BitWidth() - 1;
      if (deg_r0 < deg_r1) {
        break;
      }
      r0 ^= (r1 << (deg_r0 - deg_r1));
    }
    std::swap(r0, r1);
  }
  return r0;
}

// Rabin's exact irreducibility test for a monic degree-128 binary polynomial
// P(x) = x^128 + modulus_low over GF(2).
// Since 128 = 2^7 has sole prime divisor q = 2, P is irreducible iff:
//   (1) gcd(P(x), x^(2^64) - x) == 1, and
//   (2) x^(2^128) == x (mod P(x)).
inline bool IsBinaryPoly128Irreducible(UInt128 modulus_low) noexcept {
  // Divisible by x?
  if ((modulus_low.lo & 1ULL) == 0) {
    return false;
  }
  // Divisible by x + 1? P(1) == 0 iff popcount(modulus_low) is odd.
  if ((modulus_low.PopCount() & 1) != 0) {
    return false;
  }

  const UInt128 mu_low = ComputeBarrettMuLow128(modulus_low);
  auto square_mod = [&](UInt128 v) noexcept -> UInt128 {
    const auto [hi, lo] = ClMul128(v, v);
    return ReduceBarrett128(hi, lo, modulus_low, mu_low);
  };

  UInt128 h{2, 0};  // x mod P(x)
  for (int step = 0; step < 64; ++step) {
    h = square_mod(h);
  }
  const UInt128 diff = h ^ UInt128 { 2, 0 };  // x^(2^64) - x mod P(x)
  if (diff.IsZero() ||
      !BinaryPoly128GcdWithModulus(modulus_low, diff).IsOne()) {
    return false;
  }
  for (int step = 0; step < 64; ++step) {
    h = square_mod(h);
  }
  return h == UInt128{2, 0};
}

// Exact division of a 128-bit unsigned integer by a non-zero 64-bit divisor.
constexpr UInt128 Div128By64(UInt128 num, uint64_t den) noexcept {
  UInt128 quotient{0, 0};
  uint64_t rem = 0;
  for (int i = 127; i >= 0; --i) {
    const bool carry = (rem >> 63) != 0;
    rem = (rem << 1) | (num.Bit(i) ? 1ULL : 0ULL);
    if (carry || rem >= den) {
      rem -= den;
      quotient |= (UInt128{1, 0} << i);
    }
  }
  return quotient;
}

// Parses a binary polynomial of degree <= 127 in `variable_name` or x/X,
// e.g. "x^127 + x^7 + x^2 + x + 1", "α^3 + α + 1", "1", "0".
inline UInt128 ParseBinaryPoly128(std::string_view input,
                                  std::string_view variable_name) {
  std::string cleaned;
  cleaned.reserve(input.size());
  for (char ch : input) {
    if (std::isspace(static_cast<unsigned char>(ch)) == 0) {
      cleaned.push_back(ch);
    }
  }
  if (cleaned.empty()) {
    throw std::invalid_argument("Polynomial string cannot be empty");
  }

  std::vector<std::string_view> terms;
  size_t pos = 0;
  if (cleaned[0] == '+' || cleaned[0] == '-') {
    pos = 1;
    if (pos >= cleaned.size()) {
      throw std::invalid_argument("Invalid polynomial string: " +
                                  std::string(input));
    }
  }

  size_t term_start = pos;
  for (size_t i = pos; i < cleaned.size(); ++i) {
    const char c = cleaned[i];
    if (c == '+' || c == '-') {
      if (i > 0 && cleaned[i - 1] == '^') {
        throw std::invalid_argument("Negative exponent in polynomial: " +
                                    std::string(input));
      }
      if (i == term_start) {
        throw std::invalid_argument("Empty term in polynomial: " +
                                    std::string(input));
      }
      terms.push_back(
          std::string_view(cleaned).substr(term_start, i - term_start));
      term_start = i + 1;
    }
  }
  if (term_start >= cleaned.size()) {
    throw std::invalid_argument("Trailing operator in polynomial: " +
                                std::string(input));
  }
  terms.push_back(std::string_view(cleaned).substr(term_start));

  UInt128 poly{0, 0};
  for (std::string_view term : terms) {
    size_t var_pos = std::string_view::npos;
    size_t var_len = 0;

    if (!variable_name.empty()) {
      var_pos = term.find(variable_name);
      if (var_pos != std::string_view::npos) {
        var_len = variable_name.size();
      }
    }
    if (var_pos == std::string_view::npos) {
      var_pos = term.find('x');
      if (var_pos == std::string_view::npos) {
        var_pos = term.find('X');
      }
      if (var_pos != std::string_view::npos) {
        var_len = 1;
      }
    }

    uint64_t coeff_parity = 1;
    int deg = 0;

    if (var_pos == std::string_view::npos) {
      for (char ch : term) {
        if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
          throw std::invalid_argument("Invalid term in binary polynomial: " +
                                      std::string(term));
        }
      }
      size_t parsed = 0;
      uint64_t val = 0;
      try {
        val = std::stoull(std::string(term), &parsed, 10);
      } catch (const std::out_of_range&) {
        throw std::invalid_argument("Constant out of range in polynomial: " +
                                    std::string(term));
      }
      if (parsed != term.size()) {
        throw std::invalid_argument("Invalid constant term: " +
                                    std::string(term));
      }
      coeff_parity = val & 1ULL;
      deg = 0;
    } else {
      std::string_view prefix = term.substr(0, var_pos);
      if (!prefix.empty()) {
        if (prefix.back() == '*') {
          prefix.remove_suffix(1);
          if (prefix.empty()) {
            throw std::invalid_argument("Invalid coefficient in term: " +
                                        std::string(term));
          }
        }
        for (char ch : prefix) {
          if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
            throw std::invalid_argument("Invalid coefficient in term: " +
                                        std::string(term));
          }
        }
        size_t parsed = 0;
        uint64_t val = 0;
        try {
          val = std::stoull(std::string(prefix), &parsed, 10);
        } catch (const std::out_of_range&) {
          throw std::invalid_argument("Coefficient out of range in term: " +
                                      std::string(term));
        }
        if (parsed != prefix.size()) {
          throw std::invalid_argument("Invalid coefficient in term: " +
                                      std::string(term));
        }
        coeff_parity = val & 1ULL;
      }

      const std::string_view suffix = term.substr(var_pos + var_len);
      if (suffix.empty()) {
        deg = 1;
      } else if (suffix[0] == '^') {
        const std::string_view exp_sv = suffix.substr(1);
        if (exp_sv.empty()) {
          throw std::invalid_argument("Missing exponent in term: " +
                                      std::string(term));
        }
        for (char ch : exp_sv) {
          if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
            throw std::invalid_argument("Invalid exponent in term: " +
                                        std::string(term));
          }
        }
        size_t parsed = 0;
        uint64_t raw_deg = 0;
        try {
          raw_deg = std::stoull(std::string(exp_sv), &parsed, 10);
        } catch (const std::out_of_range&) {
          throw std::invalid_argument("Exponent out of range in term: " +
                                      std::string(term));
        }
        if (parsed != exp_sv.size()) {
          throw std::invalid_argument("Invalid exponent in term: " +
                                      std::string(term));
        }
        if (raw_deg > 127) {
          throw std::invalid_argument(
              "Polynomial degree exceeds 127 for GF(2^128): " +
              std::string(term));
        }
        deg = static_cast<int>(raw_deg);
      } else {
        throw std::invalid_argument("Unexpected characters in term: " +
                                    std::string(term));
      }
    }

    if (coeff_parity == 1) {
      poly ^= (UInt128{1, 0} << deg);
    }
  }
  return poly;
}

}  // namespace detail

// Binary extension field GF(2^128) with hardware PCLMULQDQ acceleration,
// specialized fast reduction for the AES-GCM / GHASH polynomial
// x^128 + x^7 + x^2 + x + 1, general carry-less Barrett reduction for custom
// irreducible polynomials, and built-in NIST SP 800-38D GHASH multiplication.
class GaloisFieldBinary128 final : public GaloisFieldBase<UInt128> {
 public:
  // Default AES-GCM / GHASH irreducible polynomial:
  // P(x) = x^128 + x^7 + x^2 + x + 1 (low 128 bits = 0x87).
  static constexpr UInt128 kGhashModulusLow{0x87ULL, 0ULL};

  explicit GaloisFieldBinary128(UInt128 modulus_low = kGhashModulusLow,
                                bool check_irreducible = false,
                                const std::string& rep = "hex",
                                const std::string& variable_name = "x")
      : modulus_low_(modulus_low),
        barrett_mu_low_(detail::ComputeBarrettMuLow128(modulus_low)),
        is_ghash_poly_(modulus_low == kGhashModulusLow),
        representation_(ValidateRep(utils::ConvertRepresentation(rep))),
        variable_name_(variable_name) {
    if (modulus_low_.IsZero()) {
      throw std::invalid_argument(
          "GF(2^128) modulus_low cannot be zero (x^128 is reducible)");
    }
    if (check_irreducible &&
        !detail::IsBinaryPoly128Irreducible(modulus_low_)) {
      throw std::invalid_argument(
          "Provided degree-128 polynomial is not irreducible over GF(2)");
    }
  }

  GaloisFieldBinary128(const GaloisFieldBinary128& other)
      : modulus_low_(other.modulus_low_),
        barrett_mu_low_(other.barrett_mu_low_),
        is_ghash_poly_(other.is_ghash_poly_),
        representation_(other.representation_.load(std::memory_order_relaxed)),
        variable_name_(other.variable_name_) {}

  GaloisFieldBinary128& operator=(const GaloisFieldBinary128& other) {
    if (this != &other) {
      modulus_low_ = other.modulus_low_;
      barrett_mu_low_ = other.barrett_mu_low_;
      is_ghash_poly_ = other.is_ghash_poly_;
      representation_.store(
          other.representation_.load(std::memory_order_relaxed),
          std::memory_order_relaxed);
      variable_name_ = other.variable_name_;
    }
    return *this;
  }

  static bool IsIrreducible(UInt128 modulus_low) noexcept {
    return detail::IsBinaryPoly128Irreducible(modulus_low);
  }

  uint64_t Characteristic() const override { return 2; }

  // Returns 0 because 2^128 overflows uint64_t (inspect Degree() == 128).
  uint64_t Order() const override { return 0; }

  uint32_t Degree() const override { return 128; }

  UInt128 ModulusLow() const noexcept { return modulus_low_; }

  UInt128 Add(const UInt128& a, const UInt128& b) const override {
    return a ^ b;
  }

  UInt128 Sub(const UInt128& a, const UInt128& b) const override {
    return a ^ b;
  }

  UInt128 Mul(const UInt128& a, const UInt128& b) const override {
    return MulFast(a, b);
  }

  UInt128 Neg(const UInt128& a) const override { return a; }

  UInt128 Inv(const UInt128& a) const override {
    if (a.IsZero()) {
      throw std::domain_error("Inverse of zero is undefined in GF(2^128)");
    }
    if (a.IsOne()) {
      return UInt128{1, 0};
    }

    const int deg_a = a.BitWidth() - 1;  // 1 <= deg_a <= 127
    const int shift0 = 128 - deg_a;      // 1 <= shift0 <= 127
    const UInt128 a_no_msb = a ^ (UInt128{1, 0} << deg_a);

    // First division step: P(x) = q * a + r with deg(r) < deg(a) <= 127.
    UInt128 r = modulus_low_ ^ (a_no_msb << shift0);
    UInt128 q = UInt128{1, 0} << shift0;

    while (!r.IsZero()) {
      const int deg_r = r.BitWidth() - 1;
      if (deg_r < deg_a) {
        break;
      }
      const int s = deg_r - deg_a;
      r ^= (a << s);
      q ^= (UInt128{1, 0} << s);
    }

    UInt128 r0 = a;
    UInt128 r1 = r;
    UInt128 t0{1, 0};
    UInt128 t1 = q;

    // Loop invariant: t0 * a == r0 (mod P), t1 * a == r1 (mod P),
    // and deg(t1) + deg(r0) == 128. Stopping once r1 <= 1 guarantees
    // deg(t1) + s <= 127 (zero 128-bit overflow).
    while (r1.hi != 0 || r1.lo > 1) {
      const int deg_r1 = r1.BitWidth() - 1;
      while (!r0.IsZero()) {
        const int deg_r0 = r0.BitWidth() - 1;
        if (deg_r0 < deg_r1) {
          break;
        }
        const int s = deg_r0 - deg_r1;
        r0 ^= (r1 << s);
        t0 ^= (t1 << s);
      }
      std::swap(r0, r1);
      std::swap(t0, t1);
    }

    if (r1.IsOne()) {
      return t1;
    }
    throw std::domain_error(
        "Element is not invertible in GF(2^128) (modulus may be reducible)");
  }

  UInt128 Div(const UInt128& a, const UInt128& b) const override {
    if (b.IsZero()) {
      throw std::domain_error("Division by zero in GF(2^128)");
    }
    return MulFast(a, Inv(b));
  }

  UInt128 Pow(const UInt128& a, uint64_t exp) const override {
    if (a.IsZero()) {
      return (exp == 0) ? UInt128{1, 0} : UInt128{0, 0};
    }
    if (a.IsOne() || exp == 0) {
      return UInt128{1, 0};
    }
    UInt128 result{1, 0};
    UInt128 base = a;
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        result = MulFast(result, base);
      }
      exp >>= 1;
      if (exp == 0) {
        break;
      }
      base = MulFast(base, base);
    }
    return result;
  }

  UInt128 Pow128(const UInt128& a, UInt128 exp) const noexcept {
    if (a.IsZero()) {
      return exp.IsZero() ? UInt128{1, 0} : UInt128{0, 0};
    }
    if (a.IsOne() || exp.IsZero()) {
      return UInt128{1, 0};
    }
    UInt128 result{1, 0};
    UInt128 base = a;
    while (!exp.IsZero()) {
      if ((exp.lo & 1ULL) != 0) {
        result = MulFast(result, base);
      }
      exp >>= 1;
      if (exp.IsZero()) {
        break;
      }
      base = MulFast(base, base);
    }
    return result;
  }

  // In GF(2^128), sqrt(a) = a^(2^127) via 127 repeated squarings.
  UInt128 Sqrt(const UInt128& a) const override {
    if (a.IsZero() || a.IsOne()) {
      return a;
    }
    UInt128 res = a;
    for (int i = 0; i < 127; ++i) {
      res = MulFast(res, res);
    }
    return res;
  }

  // Frobenius endomorphism: a -> a^2.
  UInt128 Frobenius(const UInt128& a) const override { return MulFast(a, a); }

  // Absolute field trace: Tr_{GF(2^128)/GF(2)}(a) = sum_{i=0}^{127} a^(2^i).
  UInt128 Trace(const UInt128& a) const override {
    UInt128 acc = a;
    UInt128 cur = a;
    for (int i = 1; i < 128; ++i) {
      cur = MulFast(cur, cur);
      acc ^= cur;
    }
    return acc;
  }

  // Absolute field norm: N_{GF(2^128)/GF(2)}(a) = a^(2^128 - 1) = (a != 0).
  UInt128 Norm(const UInt128& a) const override {
    return a.IsZero() ? UInt128{0, 0} : UInt128{1, 0};
  }

  // NIST SP 800-38D bit-reflected GHASH multiplication directly on top of Mul.
  // For the default modulus x^128 + x^7 + x^2 + x + 1, computes X * Y in the
  // standard AES-GCM block representation where bit 127 (MSB of byte 0) is t^0.
  UInt128 GhashMul(const UInt128& a, const UInt128& b) const noexcept {
    return MulFast(a.BitReverse(), b.BitReverse()).BitReverse();
  }

  uint64_t Log(const UInt128& a) const override {
    return Log(a, MultiplicativeGenerator());
  }

  uint64_t Log(const UInt128& a, const UInt128& generator) const override {
    if (a.IsZero()) {
      throw std::domain_error("Log of zero is undefined in GF(2^128)");
    }
    if (generator.IsZero()) {
      throw std::invalid_argument("Generator cannot be zero in GF(2^128)");
    }
    if (a.IsOne()) {
      return 0;
    }
    if (a == generator) {
      return 1;
    }
    throw std::runtime_error(
        "Discrete logarithm is not supported in GF(2^128) "
        "(group order 2^128 - 1 exceeds uint64_t)");
  }

  UInt128 MultiplicativeOrder128(const UInt128& a) const {
    if (a.IsZero()) {
      throw std::domain_error("Multiplicative order of zero is undefined");
    }
    if (a.IsOne()) {
      return UInt128{1, 0};
    }
    // Complete prime factorization of 2^128 - 1 (all 9 prime factors are
    // simple): 3 * 5 * 17 * 257 * 641 * 65537 * 274177 * 6700417 *
    // 67280421310721.
    UInt128 ord{~uint64_t{0}, ~uint64_t{0}};
    for (uint64_t p : kGroupOrderPrimeFactors) {
      const UInt128 candidate = detail::Div128By64(ord, p);
      if (Pow128(a, candidate).IsOne()) {
        ord = candidate;
      }
    }
    return ord;
  }

  uint64_t MultiplicativeOrder(const UInt128& a) const override {
    const UInt128 ord = MultiplicativeOrder128(a);
    if (ord.hi != 0) {
      throw std::overflow_error(
          "Multiplicative order in GF(2^128) exceeds uint64_t");
    }
    return ord.lo;
  }

  UInt128 MultiplicativeGenerator() const override {
    // x^128 + x^7 + x^2 + x + 1 is primitive, so x = UInt128{2, 0} generates
    // GF(2^128)^*.
    if (is_ghash_poly_) {
      return UInt128{2, 0};
    }
    for (uint64_t candidate = 2; candidate < ~uint64_t{0}; ++candidate) {
      const UInt128 elem{candidate, 0};
      if (IsPrimitiveElement(elem)) {
        return elem;
      }
    }
    throw std::runtime_error(
        "Could not find a multiplicative generator in GF(2^128)");
  }

  std::vector<UInt128> MultiplicativeGenerators() const override {
    throw std::length_error(
        "Field order 2^128 is too large to enumerate all multiplicative "
        "generators");
  }

  UInt128 Random() const override {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    return Random(gen);
  }

  UInt128 Random(std::mt19937_64& rng) const override {
    const uint64_t lo = rng();
    const uint64_t hi = rng();
    return UInt128{lo, hi};
  }

  UInt128 MultiplicativeIdentity() const override { return UInt128{1, 0}; }
  UInt128 AdditiveIdentity() const override { return UInt128{0, 0}; }

  bool IsZero(const UInt128& a) const override { return a.IsZero(); }
  bool IsOne(const UInt128& a) const override { return a.IsOne(); }

  UInt128 GetElementValue(const UInt128& value) const override { return value; }

  UInt128 SetElementValue(const UInt128& value) const override { return value; }

  UInt128 SetElementValue(const std::string& value_str) const override {
    const std::string_view trimmed = detail::TrimAscii128(value_str);
    if (trimmed.empty()) {
      throw std::invalid_argument("Empty string for GF(2^128) element");
    }

    // 1. Hexadecimal literal: "0x..." / "0X...".
    if (trimmed.size() >= 2 && trimmed[0] == '0' &&
        (trimmed[1] == 'x' || trimmed[1] == 'X')) {
      return UInt128::FromHex(trimmed);
    }

    // 2. Pure decimal integer literal.
    bool all_digits = true;
    for (char ch : trimmed) {
      if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
        all_digits = false;
        break;
      }
    }
    if (all_digits) {
      return UInt128::FromDecimal(trimmed);
    }

    // 3. Binary polynomial string (e.g. "x^7 + x^2 + x + 1").
    return detail::ParseBinaryPoly128(trimmed, variable_name_);
  }

  void Print(std::ostream& os) const override { os << "GF(2^128)"; }

  void Print(const UInt128& a, std::ostream& os) const override {
    detail::StreamFlagsGuard128 guard(os);
    const FieldRepresentation rep = GetRepresentation();
    if (rep == FieldRepresentation::HEX) {
      os << a.ToHex();
    } else if (rep == FieldRepresentation::INT) {
      os << a.ToDecimal();
    } else if (rep == FieldRepresentation::POLY) {
      os << PolyToString(a);
    } else {
      throw std::invalid_argument("Unsupported representation for GF(2^128)");
    }
  }

  std::string ToString(const UInt128& a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_.load(std::memory_order_relaxed);
  }

  void SetRepresentation(FieldRepresentation rep) override {
    representation_.store(ValidateRep(rep), std::memory_order_relaxed);
  }

  size_t MemoryBytes() const override {
    return sizeof(*this) + variable_name_.capacity();
  }

  void AddArray(std::span<UInt128> dst, std::span<const UInt128> a,
                std::span<const UInt128> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("AddArray: span sizes must match");
    }
    UInt128* out = dst.data();
    const UInt128* pa = a.data();
    const UInt128* pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = pa[i] ^ pb[i];
    }
  }

  void SubArray(std::span<UInt128> dst, std::span<const UInt128> a,
                std::span<const UInt128> b) const override {
    AddArray(dst, a, b);
  }

  void MulArray(std::span<UInt128> dst, std::span<const UInt128> a,
                std::span<const UInt128> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match");
    }
    UInt128* out = dst.data();
    const UInt128* pa = a.data();
    const UInt128* pb = b.data();
    if (is_ghash_poly_) {
      for (size_t i = 0; i < n; ++i) {
        out[i] = detail::MulGhashPoly128(pa[i], pb[i]);
      }
    } else {
      for (size_t i = 0; i < n; ++i) {
        out[i] = MulBarrett(pa[i], pb[i]);
      }
    }
  }

  void ScalarMulArray(std::span<UInt128> dst, std::span<const UInt128> a,
                      const UInt128& scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match");
    }
    if (scalar.IsZero()) {
      std::fill(dst.begin(), dst.end(), UInt128{0, 0});
      return;
    }
    if (scalar.IsOne()) {
      if (dst.data() != a.data()) {
        std::copy(a.begin(), a.end(), dst.begin());
      }
      return;
    }
    UInt128* out = dst.data();
    const UInt128* pa = a.data();
    if (is_ghash_poly_) {
      for (size_t i = 0; i < n; ++i) {
        out[i] = detail::MulGhashPoly128(pa[i], scalar);
      }
    } else {
      for (size_t i = 0; i < n; ++i) {
        out[i] = MulBarrett(pa[i], scalar);
      }
    }
  }

  void FmaArray(std::span<UInt128> dst, std::span<const UInt128> a,
                const UInt128& scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match");
    }
    if (scalar.IsZero()) {
      return;
    }
    UInt128* out = dst.data();
    const UInt128* pa = a.data();
    if (scalar.IsOne()) {
      for (size_t i = 0; i < n; ++i) {
        out[i] ^= pa[i];
      }
      return;
    }
    if (is_ghash_poly_) {
      for (size_t i = 0; i < n; ++i) {
        out[i] ^= detail::MulGhashPoly128(pa[i], scalar);
      }
    } else {
      for (size_t i = 0; i < n; ++i) {
        out[i] ^= MulBarrett(pa[i], scalar);
      }
    }
  }

  // Lazy-reduction dot product: accumulates exact unreduced 256-bit carry-less
  // products across all terms and performs a single modular reduction at the
  // end.
  UInt128 DotProduct(std::span<const UInt128> a,
                     std::span<const UInt128> b) const override {
    const size_t n = a.size();
    if (b.size() != n) {
      throw std::invalid_argument("DotProduct: span sizes must match");
    }
    UInt128 acc_hi{0, 0};
    UInt128 acc_lo{0, 0};
    const UInt128* pa = a.data();
    const UInt128* pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      const auto [hi, lo] = detail::ClMul128(pa[i], pb[i]);
      acc_hi ^= hi;
      acc_lo ^= lo;
    }
    if (is_ghash_poly_) {
      return detail::ReduceGhashPoly128(acc_hi, acc_lo);
    }
    return detail::ReduceBarrett128(acc_hi, acc_lo, modulus_low_,
                                    barrett_mu_low_);
  }

  // Montgomery's batch inversion (safe for both out-of-place and in-place).
  void BatchInv(std::span<UInt128> dst,
                std::span<const UInt128> src) const override {
    const size_t n = src.size();
    if (dst.size() != n) {
      throw std::invalid_argument("BatchInv: dst and src sizes must match");
    }
    if (n == 0) {
      return;
    }
    std::vector<UInt128> prefix(n);
    prefix[0] = src[0];
    if (prefix[0].IsZero()) {
      throw std::domain_error("BatchInv: inverse of zero is undefined");
    }
    for (size_t i = 1; i < n; ++i) {
      if (src[i].IsZero()) {
        throw std::domain_error("BatchInv: inverse of zero is undefined");
      }
      prefix[i] = MulFast(prefix[i - 1], src[i]);
    }
    UInt128 inv_acc = Inv(prefix[n - 1]);
    for (size_t i = n - 1; i > 0; --i) {
      const UInt128 cur = src[i];
      dst[i] = MulFast(inv_acc, prefix[i - 1]);
      inv_acc = MulFast(inv_acc, cur);
    }
    dst[0] = inv_acc;
  }

 private:
  static constexpr std::array<uint64_t, 9> kGroupOrderPrimeFactors = {
      3ULL,     5ULL,      17ULL,      257ULL,           641ULL,
      65537ULL, 274177ULL, 6700417ULL, 67280421310721ULL};

  static FieldRepresentation ValidateRep(FieldRepresentation rep) {
    if (rep != FieldRepresentation::HEX && rep != FieldRepresentation::INT &&
        rep != FieldRepresentation::POLY) {
      throw std::invalid_argument(
          "GF(2^128) supports HEX, INT, or POLY representations");
    }
    return rep;
  }

  inline UInt128 MulBarrett(UInt128 a, UInt128 b) const noexcept {
    const auto [hi, lo] = detail::ClMul128(a, b);
    return detail::ReduceBarrett128(hi, lo, modulus_low_, barrett_mu_low_);
  }

  inline UInt128 MulFast(UInt128 a, UInt128 b) const noexcept {
    if (is_ghash_poly_) {
      return detail::MulGhashPoly128(a, b);
    }
    return MulBarrett(a, b);
  }

  bool IsPrimitiveElement(UInt128 elem) const noexcept {
    if (elem.IsZero() || elem.IsOne()) {
      return false;
    }
    const UInt128 group_order{~uint64_t{0}, ~uint64_t{0}};
    if (!Pow128(elem, group_order).IsOne()) {
      return false;
    }
    for (uint64_t p : kGroupOrderPrimeFactors) {
      const UInt128 exp = detail::Div128By64(group_order, p);
      if (Pow128(elem, exp).IsOne()) {
        return false;
      }
    }
    return true;
  }

  std::string PolyToString(UInt128 poly) const {
    if (poly.IsZero()) {
      return "0";
    }
    std::ostringstream oss;
    bool first = true;
    for (int i = poly.BitWidth() - 1; i >= 0; --i) {
      if (poly.Bit(i)) {
        if (!first) {
          oss << " + ";
        }
        first = false;
        if (i == 0) {
          oss << "1";
        } else if (i == 1) {
          oss << variable_name_;
        } else {
          oss << variable_name_ << "^" << i;
        }
      }
    }
    return oss.str();
  }

  UInt128 modulus_low_;
  UInt128 barrett_mu_low_;
  bool is_ghash_poly_;
  std::atomic<FieldRepresentation> representation_;
  std::string variable_name_;
};

using GF128 = GaloisFieldBinary128;

static_assert(std::is_trivially_copyable_v<UInt128>);
static_assert(std::is_standard_layout_v<UInt128>);
static_assert(sizeof(UInt128) == 16);
static_assert(FiniteField<GaloisFieldBinary128>);

}  // namespace xg

namespace std {

template <>
struct hash<xg::UInt128> {
  size_t operator()(const xg::UInt128& v) const noexcept {
    const size_t h1 = std::hash<uint64_t>{}(v.lo);
    const size_t h2 = std::hash<uint64_t>{}(v.hi);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

}  // namespace std

#endif  // XGALOIS_FIELD_GF128_HPP_
