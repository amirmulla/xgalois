#ifndef XGALOIS_FIELD_GF_BINARY_LARGE_HPP_
#define XGALOIS_FIELD_GF_BINARY_LARGE_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {

namespace detail {

// Inserts a 0 bit between every pair of bits in a 32-bit word (squaring over
// GF(2)[x]).
constexpr uint64_t SpreadBits32(uint32_t w) noexcept {
  uint64_t x = static_cast<uint64_t>(w);
  x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
  x = (x | (x << 8)) & 0x00FF00FF00FF00FFULL;
  x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0FULL;
  x = (x | (x << 2)) & 0x3333333333333333ULL;
  x = (x | (x << 1)) & 0x5555555555555555ULL;
  return x;
}

// Extracts even-indexed bits (0, 2, ..., 62) of a 64-bit word into bits 0..31.
constexpr uint32_t CompactEvenBits64(uint64_t x) noexcept {
  x &= 0x5555555555555555ULL;
  x = (x | (x >> 1)) & 0x3333333333333333ULL;
  x = (x | (x >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
  x = (x | (x >> 4)) & 0x00FF00FF00FF00FFULL;
  x = (x | (x >> 8)) & 0x0000FFFF0000FFFFULL;
  x = (x | (x >> 16)) & 0x00000000FFFFFFFFULL;
  return static_cast<uint32_t>(x);
}

// Recursive zero-allocation Karatsuba carry-less multiplication of n-limb
// binary polynomials (1 <= n <= 9) into a 2*n-limb output buffer.
inline void ClMulKaratsubaImpl(const uint64_t* a, const uint64_t* b,
                               std::size_t n, uint64_t* out) noexcept {
  for (std::size_t i = 0; i < 2 * n; ++i) {
    out[i] = 0ULL;
  }
  if (n <= 2) {
    for (std::size_t i = 0; i < n; ++i) {
      if (a[i] == 0ULL) {
        continue;
      }
      for (std::size_t j = 0; j < n; ++j) {
        if (b[j] == 0ULL) {
          continue;
        }
        const auto [hi, lo] = ClMul64(a[i], b[j]);
        out[i + j] ^= lo;
        out[i + j + 1] ^= hi;
      }
    }
    return;
  }

  constexpr std::size_t kMaxHalf = 5;  // ceil(9 / 2) == 5.
  const std::size_t lo = n / 2;
  const std::size_t hi = n - lo;

  std::array<uint64_t, 2 * kMaxHalf> z0{};
  std::array<uint64_t, 2 * kMaxHalf> z2{};
  std::array<uint64_t, 2 * kMaxHalf> z1{};
  std::array<uint64_t, kMaxHalf> a_sum{};
  std::array<uint64_t, kMaxHalf> b_sum{};

  ClMulKaratsubaImpl(a, b, lo, z0.data());
  ClMulKaratsubaImpl(a + lo, b + lo, hi, z2.data());

  for (std::size_t i = 0; i < hi; ++i) {
    a_sum[i] = ((i < lo) ? a[i] : 0ULL) ^ a[lo + i];
    b_sum[i] = ((i < lo) ? b[i] : 0ULL) ^ b[lo + i];
  }
  ClMulKaratsubaImpl(a_sum.data(), b_sum.data(), hi, z1.data());

  for (std::size_t i = 0; i < 2 * lo; ++i) {
    z1[i] ^= z0[i];
  }
  for (std::size_t i = 0; i < 2 * hi; ++i) {
    z1[i] ^= z2[i];
  }

  for (std::size_t i = 0; i < 2 * lo; ++i) {
    out[i] ^= z0[i];
  }
  for (std::size_t i = 0; i < 2 * hi; ++i) {
    out[lo + i] ^= z1[i];
  }
  for (std::size_t i = 0; i < 2 * hi; ++i) {
    out[2 * lo + i] ^= z2[i];
  }
}

}  // namespace detail

// ============================================================================
// GaloisFieldBinaryLarge<M>: Large Binary Extension Field GF(2^M)
// (65 <= M <= 576)
// ============================================================================
//
// Features:
//   - BigUInt<kBits> element storage (kBits = ceil(M/64)*64)
//   - Zero-allocation Karatsuba carry-less multiplication + fast sparse modular
//     reduction for trinomials/pentanomials (and arbitrary sparse moduli)
//   - Bit-interleaving squaring Sqr()
//   - Itoh-Tsujii inversion Inv()
//   - Precomputed sqrt(x) square root Sqrt() in O(1) Mul
//   - Precomputed O(1) absolute trace mask Trace()
//   - Half-trace HalfTrace() for odd M solving z^2 + z = a when Tr(a) == 0
template <std::size_t M>
  requires(M >= 65 && M <= 576)
class GaloisFieldBinaryLarge {
 public:
  static constexpr std::size_t kDegree = M;
  static constexpr std::size_t kBits = ((M + 63) / 64) * 64;
  static constexpr std::size_t kLimbs = kBits / 64;

  using element_type = BigUInt<kBits>;
  using storage_type = BigUInt<kBits>;

  explicit GaloisFieldBinaryLarge(std::vector<std::size_t> exponents) {
    if (exponents.empty()) {
      throw std::invalid_argument(
          "GaloisFieldBinaryLarge: reduction polynomial exponents cannot be "
          "empty");
    }
    std::sort(exponents.begin(), exponents.end(), std::greater<std::size_t>{});
    exponents.erase(std::unique(exponents.begin(), exponents.end()),
                    exponents.end());

    if (exponents.front() != M) {
      throw std::invalid_argument(
          "GaloisFieldBinaryLarge: maximum exponent must equal M (" +
          std::to_string(M) + ")");
    }
    if (exponents.back() != 0) {
      throw std::invalid_argument(
          "GaloisFieldBinaryLarge: irreducible binary polynomial must include "
          "constant term 0");
    }

    low_exponents_.assign(exponents.begin() + 1, exponents.end());

    top_mask_ = (M % 64 == 0) ? ~uint64_t{0} : ((uint64_t{1} << (M % 64)) - 1);

    // Precompute sqrt(x) = x^(2^(M-1)) mod f(x).
    element_type sx(2ULL);
    for (std::size_t i = 0; i < M - 1; ++i) {
      sx = Sqr(sx);
    }
    sqrt_x_ = sx;

    // Precompute O(1) absolute trace mask: bit i is Tr(x^i) in {0, 1}.
    trace_mask_ = BigUInt<kBits>::Zero();
    std::vector<uint8_t> tr_bits(M, 0);
    tr_bits[0] = static_cast<uint8_t>(M & 1U);
    if (tr_bits[0] != 0U) {
      trace_mask_.SetBit(0);
    }
    for (std::size_t i = 1; i < M; ++i) {
      if ((i & 1U) == 0U) {
        // Tr(x^(2k)) == Tr((x^k)^2) == Tr(x^k).
        tr_bits[i] = tr_bits[i >> 1];
      } else {
        element_type xi = BigUInt<kBits>::Zero();
        xi.SetBit(i);
        element_type acc = xi;
        element_type cur = xi;
        for (std::size_t step = 1; step < M; ++step) {
          cur = Sqr(cur);
          acc ^= cur;
        }
        tr_bits[i] = static_cast<uint8_t>(acc.Limb(0) & 1ULL);
      }
      if (tr_bits[i] != 0U) {
        trace_mask_.SetBit(i);
      }
    }
  }

  uint64_t Characteristic() const noexcept { return 2; }
  uint32_t Degree() const noexcept { return static_cast<uint32_t>(M); }
  BigInt BigCharacteristic() const { return BigInt(2); }
  BigInt BigOrder() const { return BigInt(1) << M; }

  const std::vector<std::size_t>& LowExponents() const noexcept {
    return low_exponents_;
  }

  element_type AdditiveIdentity() const noexcept {
    return BigUInt<kBits>::Zero();
  }

  element_type MultiplicativeIdentity() const noexcept {
    return BigUInt<kBits>::One();
  }

  element_type Reduce(const element_type& a) const noexcept {
    if constexpr (M % 64 == 0) {
      return a;
    } else {
      std::array<uint64_t, kLimbs> limbs = a.Limbs();
      limbs[kLimbs - 1] &= top_mask_;
      return BigUInt<kBits>(limbs);
    }
  }

  bool IsZero(const element_type& a) const noexcept {
    return Reduce(a).IsZero();
  }

  bool IsOne(const element_type& a) const noexcept { return Reduce(a).IsOne(); }

  bool Equals(const element_type& a, const element_type& b) const noexcept {
    return Reduce(a) == Reduce(b);
  }

  element_type Add(const element_type& a,
                   const element_type& b) const noexcept {
    return Reduce(a ^ b);
  }

  element_type Sub(const element_type& a,
                   const element_type& b) const noexcept {
    return Reduce(a ^ b);
  }

  element_type Neg(const element_type& a) const noexcept { return Reduce(a); }

  // Karatsuba carry-less multiplication followed by fast sparse reduction.
  element_type Mul(const element_type& a,
                   const element_type& b) const noexcept {
    const element_type ar = Reduce(a);
    const element_type br = Reduce(b);
    if (ar.IsZero() || br.IsZero()) {
      return BigUInt<kBits>::Zero();
    }
    std::array<uint64_t, 2 * kLimbs> wide{};
    detail::ClMulKaratsubaImpl(ar.Limbs().data(), br.Limbs().data(), kLimbs,
                               wide.data());
    return ReduceWide(wide);
  }

  // Fast bit-interleaving squaring followed by sparse reduction.
  element_type Sqr(const element_type& a) const noexcept {
    const element_type ar = Reduce(a);
    if (ar.IsZero()) {
      return BigUInt<kBits>::Zero();
    }
    std::array<uint64_t, 2 * kLimbs> wide{};
    const auto& limbs = ar.Limbs();
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const uint32_t lo32 = static_cast<uint32_t>(limbs[i]);
      const uint32_t hi32 = static_cast<uint32_t>(limbs[i] >> 32);
      wide[2 * i] = detail::SpreadBits32(lo32);
      wide[2 * i + 1] = detail::SpreadBits32(hi32);
    }
    return ReduceWide(wide);
  }

  element_type Pow(const element_type& a, uint64_t exp) const noexcept {
    element_type result = BigUInt<kBits>::One();
    element_type base = Reduce(a);
    while (exp > 0ULL) {
      if ((exp & 1ULL) != 0ULL) {
        result = Mul(result, base);
      }
      exp >>= 1;
      if (exp > 0ULL) {
        base = Sqr(base);
      }
    }
    return result;
  }

  element_type PowBig(const element_type& a, const BigInt& exp) const {
    if (exp.IsZero()) {
      return BigUInt<kBits>::One();
    }
    if (exp.IsNegative()) {
      return PowBig(Inv(a), -exp);
    }
    element_type result = BigUInt<kBits>::One();
    element_type base = Reduce(a);
    const std::size_t bits = exp.BitLength();
    for (std::size_t i = 0; i < bits; ++i) {
      if (exp.GetBit(i)) {
        result = Mul(result, base);
      }
      if (i + 1 < bits) {
        base = Sqr(base);
      }
    }
    return result;
  }

  // Itoh-Tsujii inversion: a^(-1) = (a^(2^(M-1) - 1))^2 using O(log M)
  // multiplications and M - 1 fast bit-interleaving squarings.
  element_type Inv(const element_type& a) const {
    const element_type ar = Reduce(a);
    if (ar.IsZero()) {
      throw std::domain_error(
          "GaloisFieldBinaryLarge::Inv: inverse of zero is undefined");
    }
    if (ar.IsOne()) {
      return BigUInt<kBits>::One();
    }

    constexpr std::size_t kExp = M - 1;
    const int msb = static_cast<int>(std::bit_width(kExp)) - 1;
    element_type res = ar;  // res = a^(2^k - 1), initially k = 1.
    std::size_t k = 1;
    for (int bit = msb - 1; bit >= 0; --bit) {
      element_type shifted = res;
      for (std::size_t s = 0; s < k; ++s) {
        shifted = Sqr(shifted);
      }
      res = Mul(shifted, res);  // a^(2^(2k) - 1)
      k <<= 1;
      if (((kExp >> bit) & 1U) != 0U) {
        res = Mul(Sqr(res), ar);  // a^(2^(k+1) - 1)
        ++k;
      }
    }
    return Sqr(res);
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return Mul(a, Inv(b));
  }

  // Fast square root using precomputed sqrt(x) = x^(2^(M-1)):
  //   a(x) = a_even(x)^2 + x * a_odd(x)^2 => sqrt(a) = a_even + sqrt(x)*a_odd.
  element_type Sqrt(const element_type& a) const noexcept {
    const element_type ar = Reduce(a);
    if (ar.IsZero() || ar.IsOne()) {
      return ar;
    }
    std::array<uint64_t, kLimbs> even_limbs{};
    std::array<uint64_t, kLimbs> odd_limbs{};
    const auto& limbs = ar.Limbs();
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const uint64_t ev =
          static_cast<uint64_t>(detail::CompactEvenBits64(limbs[i]));
      const uint64_t od =
          static_cast<uint64_t>(detail::CompactEvenBits64(limbs[i] >> 1));
      const unsigned shift = static_cast<unsigned>((i & 1U) * 32U);
      even_limbs[i >> 1] |= (ev << shift);
      odd_limbs[i >> 1] |= (od << shift);
    }
    const element_type a_even(even_limbs);
    const element_type a_odd(odd_limbs);
    return Add(a_even, Mul(sqrt_x_, a_odd));
  }

  // O(1) absolute trace Tr_{GF(2^M)/GF(2)}(a) via precomputed trace_mask_.
  uint8_t Trace(const element_type& a) const noexcept {
    const element_type masked = Reduce(a) & trace_mask_;
    int parity = 0;
    for (uint64_t w : masked.Limbs()) {
      parity ^= (std::popcount(w) & 1);
    }
    return static_cast<uint8_t>(parity & 1);
  }

  // Half-trace H(a) = sum_{i=0}^{(M-1)/2} a^(2^(2i)) for odd M.
  // Satisfies H(a)^2 + H(a) = a whenever Trace(a) == 0.
  element_type HalfTrace(const element_type& a) const {
    if constexpr (M % 2 == 0) {
      throw std::domain_error(
          "GaloisFieldBinaryLarge::HalfTrace is only defined for odd degree M");
    } else {
      element_type acc = Reduce(a);
      element_type cur = acc;
      for (std::size_t i = 0; i < (M - 1) / 2; ++i) {
        cur = Sqr(Sqr(cur));
        acc ^= cur;
      }
      return acc;
    }
  }

  std::string ToString(const element_type& a) const {
    return Reduce(a).ToHex(/*prefix=*/true);
  }

 private:
  element_type ReduceWide(std::array<uint64_t, 2 * kLimbs>& t) const noexcept {
    auto xor_word_at = [&t](uint64_t w, std::size_t bit_pos) noexcept {
      if (w == 0ULL) {
        return;
      }
      const std::size_t limb = bit_pos / 64;
      const unsigned shift = static_cast<unsigned>(bit_pos % 64);
      t[limb] ^= (w << shift);
      if (shift != 0U && limb + 1 < t.size()) {
        t[limb + 1] ^= (w >> (64U - shift));
      }
    };

    auto has_high_bits = [this, &t]() noexcept -> bool {
      for (std::size_t i = kLimbs; i < 2 * kLimbs; ++i) {
        if (t[i] != 0ULL) {
          return true;
        }
      }
      if constexpr (M % 64 != 0) {
        if ((t[kLimbs - 1] & ~top_mask_) != 0ULL) {
          return true;
        }
      }
      return false;
    };

    while (has_high_bits()) {
      for (std::size_t i = 2 * kLimbs; i-- > kLimbs;) {
        const uint64_t w = t[i];
        if (w == 0ULL) {
          continue;
        }
        t[i] = 0ULL;
        const std::size_t base_shift = i * 64 - M;
        for (std::size_t e : low_exponents_) {
          xor_word_at(w, base_shift + e);
        }
      }
      if constexpr (M % 64 != 0) {
        const uint64_t hi = t[kLimbs - 1] >> (M % 64);
        if (hi != 0ULL) {
          t[kLimbs - 1] &= top_mask_;
          for (std::size_t e : low_exponents_) {
            xor_word_at(hi, e);
          }
        }
      }
    }

    std::array<uint64_t, kLimbs> out{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      out[i] = t[i];
    }
    return BigUInt<kBits>(out);
  }

  std::vector<std::size_t> low_exponents_;
  uint64_t top_mask_ = ~uint64_t{0};
  element_type sqrt_x_;
  element_type trace_mask_;
};

template <std::size_t M>
  requires(M >= 65 && M <= 576)
struct AlgebraTraits<GaloisFieldBinaryLarge<M>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = true;
};

// ============================================================================
// Standard Large Binary Field Presets (NIST / SECG Curves & Cryptography)
// ============================================================================

// NIST B-163 / K-163: GF(2^163) with modulus x^163 + x^7 + x^6 + x^3 + 1.
inline GaloisFieldBinaryLarge<163> CreateGF2_163() {
  return GaloisFieldBinaryLarge<163>({163, 7, 6, 3, 0});
}

// NIST B-233 / K-233: GF(2^233) with modulus x^233 + x^74 + 1.
inline GaloisFieldBinaryLarge<233> CreateGF2_233() {
  return GaloisFieldBinaryLarge<233>({233, 74, 0});
}

// GF(2^256) with pentanomial modulus x^256 + x^10 + x^5 + x^2 + 1.
inline GaloisFieldBinaryLarge<256> CreateGF2_256() {
  return GaloisFieldBinaryLarge<256>({256, 10, 5, 2, 0});
}

// NIST B-283 / K-283: GF(2^283) with modulus x^283 + x^12 + x^7 + x^5 + 1.
inline GaloisFieldBinaryLarge<283> CreateGF2_283() {
  return GaloisFieldBinaryLarge<283>({283, 12, 7, 5, 0});
}

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_BINARY_LARGE_HPP_
