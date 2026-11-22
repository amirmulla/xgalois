#ifndef XGALOIS_RING_MODULAR_RING_HPP_
#define XGALOIS_RING_MODULAR_RING_HPP_

#include <algorithm>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"

namespace xg {
namespace ring {

namespace detail {

__extension__ using u128 = unsigned __int128;
__extension__ using i128 = __int128;

inline constexpr uint64_t ModRingGcdU64(uint64_t a, uint64_t b) noexcept {
  while (b != 0) {
    const uint64_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

inline uint64_t ModRingMulU64(uint64_t a, uint64_t b, uint64_t mod) noexcept {
  return static_cast<uint64_t>((static_cast<u128>(a) * b) % mod);
}

inline uint64_t ModRingPowU64(uint64_t base, uint64_t exp,
                              uint64_t mod) noexcept {
  uint64_t res = 1 % mod;
  base %= mod;
  while (exp > 0) {
    if ((exp & 1ULL) != 0) {
      res = ModRingMulU64(res, base, mod);
    }
    base = ModRingMulU64(base, base, mod);
    exp >>= 1;
  }
  return res;
}

inline bool ModRingIsPrimeU64(uint64_t n) noexcept {
  if (n < 2) {
    return false;
  }
  static constexpr uint64_t kSmallPrimes[] = {2,  3,  5,  7,  11, 13,
                                              17, 19, 23, 29, 31, 37};
  for (uint64_t p : kSmallPrimes) {
    if (n == p) {
      return true;
    }
    if (n % p == 0) {
      return false;
    }
  }
  if (n < 41ULL * 41ULL) {
    return true;
  }

  static constexpr uint64_t kWitnesses[] = {
      2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};

  uint64_t d = n - 1;
  const int s = std::countr_zero(d);
  d >>= s;

  for (uint64_t a : kWitnesses) {
    const uint64_t a_mod = a % n;
    if (a_mod == 0) {
      continue;
    }
    uint64_t x = ModRingPowU64(a_mod, d, n);
    if (x == 1 || x == n - 1) {
      continue;
    }
    bool composite = true;
    for (int r = 1; r < s; ++r) {
      x = ModRingMulU64(x, x, n);
      if (x == n - 1) {
        composite = false;
        break;
      }
      if (x == 1) {
        return false;
      }
    }
    if (composite) {
      return false;
    }
  }
  return true;
}

inline uint64_t ModRingPollardRhoU64(uint64_t n) noexcept {
  if ((n & 1ULL) == 0) {
    return 2;
  }
  constexpr uint64_t kBatch = 128;
  for (uint64_t attempt = 1; attempt <= 64; ++attempt) {
    const uint64_t c = attempt;
    uint64_t y = attempt + 1;
    uint64_t q = 1;
    uint64_t x = y;
    uint64_t ys = y;
    uint64_t r = 1;
    uint64_t g = 1;

    auto step = [n, c](uint64_t v) noexcept -> uint64_t {
      return static_cast<uint64_t>((static_cast<u128>(v) * v + c) % n);
    };

    while (g == 1) {
      x = y;
      for (uint64_t i = 0; i < r; ++i) {
        y = step(y);
      }
      uint64_t k = 0;
      while (k < r && g == 1) {
        ys = y;
        const uint64_t limit = std::min(kBatch, r - k);
        for (uint64_t i = 0; i < limit; ++i) {
          y = step(y);
          const uint64_t diff = (x >= y) ? (x - y) : (y - x);
          if (diff == 0) {
            g = n;
            break;
          }
          const uint64_t next_q = ModRingMulU64(q, diff, n);
          if (next_q == 0) {
            g = ModRingGcdU64(diff, n);
            break;
          }
          q = next_q;
        }
        if (g == 1) {
          g = ModRingGcdU64(q, n);
        }
        k += limit;
      }
      r <<= 1;
    }

    if (g == n) {
      for (uint64_t cnt = 0; cnt <= kBatch + 2; ++cnt) {
        ys = step(ys);
        const uint64_t diff = (x >= ys) ? (x - ys) : (ys - x);
        g = ModRingGcdU64(diff, n);
        if (g > 1) {
          break;
        }
      }
    }

    if (g > 1 && g < n) {
      return g;
    }
  }
  return n;
}

inline void CollectPrimeFactorsU64(uint64_t n, std::vector<uint64_t>& factors) {
  if (n <= 1) {
    return;
  }
  if (ModRingIsPrimeU64(n)) {
    factors.push_back(n);
    return;
  }
  uint64_t d = ModRingPollardRhoU64(n);
  if (d == n) {
    // Fallback trial division if Pollard-Rho fails to split.
    for (uint64_t p = 3; p <= n / p; p += 2) {
      if (n % p == 0) {
        d = p;
        break;
      }
    }
  }
  if (d == n) {
    factors.push_back(n);
    return;
  }
  CollectPrimeFactorsU64(d, factors);
  CollectPrimeFactorsU64(n / d, factors);
}

inline uint64_t EulerTotientU64(uint64_t n) {
  if (n <= 1) {
    return n;
  }
  uint64_t result = n;
  uint64_t rem = n;

  static constexpr uint64_t kSmallPrimes[] = {2,  3,  5,  7,  11, 13, 17, 19,
                                              23, 29, 31, 37, 41, 43, 47, 53};
  for (uint64_t p : kSmallPrimes) {
    if (rem % p == 0) {
      while (rem % p == 0) {
        rem /= p;
      }
      result -= result / p;
    }
  }

  if (rem > 1) {
    std::vector<uint64_t> prime_factors;
    CollectPrimeFactorsU64(rem, prime_factors);
    std::sort(prime_factors.begin(), prime_factors.end());
    prime_factors.erase(std::unique(prime_factors.begin(), prime_factors.end()),
                        prime_factors.end());
    for (uint64_t p : prime_factors) {
      result -= result / p;
    }
  }
  return result;
}

}  // namespace detail

// Residue class ring Z/nZ for any modulus n >= 2 (prime or composite).
// Satisfies xg::CommutativeRing.
template <typename UIntT = uint64_t>
  requires(std::unsigned_integral<UIntT> &&
           !std::same_as<std::remove_cv_t<UIntT>, bool> &&
           sizeof(UIntT) <= sizeof(uint64_t))
class ModularRing {
 public:
  using element_type = UIntT;

  explicit ModularRing(UIntT modulus) : modulus_(modulus) {
    if (modulus_ < UIntT{2}) {
      throw std::invalid_argument("ModularRing: modulus must be >= 2");
    }
  }

  UIntT Modulus() const noexcept { return modulus_; }

  element_type Reduce(element_type a) const noexcept {
    return static_cast<element_type>(a % modulus_);
  }

  // Overflow-safe modular addition even when modulus_ == max(UIntT).
  element_type Add(const element_type& a,
                   const element_type& b) const noexcept {
    const UIntT x = static_cast<UIntT>(a % modulus_);
    const UIntT y = static_cast<UIntT>(b % modulus_);
    const UIntT rem = static_cast<UIntT>(modulus_ - x);
    return (y >= rem) ? static_cast<UIntT>(y - rem) : static_cast<UIntT>(x + y);
  }

  element_type Sub(const element_type& a,
                   const element_type& b) const noexcept {
    const UIntT x = static_cast<UIntT>(a % modulus_);
    const UIntT y = static_cast<UIntT>(b % modulus_);
    return (x >= y) ? static_cast<UIntT>(x - y)
                    : static_cast<UIntT>(modulus_ - (y - x));
  }

  element_type Neg(const element_type& a) const noexcept {
    const UIntT x = static_cast<UIntT>(a % modulus_);
    return (x == UIntT{0}) ? UIntT{0} : static_cast<UIntT>(modulus_ - x);
  }

  // Overflow-safe modular multiplication using __uint128_t for 64-bit UIntT.
  element_type Mul(const element_type& a,
                   const element_type& b) const noexcept {
    const UIntT x = static_cast<UIntT>(a % modulus_);
    const UIntT y = static_cast<UIntT>(b % modulus_);
    if constexpr (sizeof(UIntT) <= 4) {
      return static_cast<UIntT>(
          (static_cast<uint64_t>(x) * static_cast<uint64_t>(y)) %
          static_cast<uint64_t>(modulus_));
    } else {
      return static_cast<UIntT>(
          (static_cast<detail::u128>(x) * static_cast<detail::u128>(y)) %
          static_cast<detail::u128>(modulus_));
    }
  }

  element_type Pow(const element_type& a, uint64_t exp) const noexcept {
    element_type acc = MultiplicativeIdentity();
    element_type base = Reduce(a);
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  element_type AdditiveIdentity() const noexcept { return UIntT{0}; }

  element_type MultiplicativeIdentity() const noexcept { return UIntT{1}; }

  bool IsZero(const element_type& a) const noexcept {
    return (a % modulus_) == UIntT{0};
  }

  bool IsOne(const element_type& a) const noexcept {
    return (a % modulus_) == UIntT{1};
  }

  bool Equals(const element_type& a, const element_type& b) const noexcept {
    return (a % modulus_) == (b % modulus_);
  }

  std::string ToString(const element_type& a) const {
    return std::to_string(static_cast<uint64_t>(a % modulus_));
  }

  bool IsUnit(const element_type& a) const noexcept {
    const uint64_t x = static_cast<uint64_t>(a % modulus_);
    return x != 0 &&
           detail::ModRingGcdU64(x, static_cast<uint64_t>(modulus_)) == 1;
  }

  // Computes multiplicative inverse modulo modulus_ using signed 128-bit EEA
  // (overflow-free for all uint64_t moduli up to UINT64_MAX).
  std::optional<element_type> TryInv(const element_type& a) const noexcept {
    const uint64_t x = static_cast<uint64_t>(a % modulus_);
    if (x == 0) {
      return std::nullopt;
    }
    detail::i128 old_r = static_cast<detail::i128>(x);
    detail::i128 r = static_cast<detail::i128>(modulus_);
    detail::i128 old_s = 1;
    detail::i128 s = 0;

    while (r != 0) {
      const detail::i128 q = old_r / r;
      const detail::i128 next_r = old_r - q * r;
      old_r = r;
      r = next_r;

      const detail::i128 next_s = old_s - q * s;
      old_s = s;
      s = next_s;
    }

    if (old_r != 1) {
      return std::nullopt;
    }
    const detail::i128 mod128 = static_cast<detail::i128>(modulus_);
    detail::i128 inv = old_s % mod128;
    if (inv < 0) {
      inv += mod128;
    }
    return static_cast<element_type>(inv);
  }

  // Returns true iff a != 0 (mod n) and gcd(a, n) > 1.
  bool IsZeroDivisor(const element_type& a) const noexcept {
    const uint64_t x = static_cast<uint64_t>(a % modulus_);
    return x != 0 &&
           detail::ModRingGcdU64(x, static_cast<uint64_t>(modulus_)) > 1;
  }

  // Enumerates all elements of the unit group (Z/nZ)^x in ascending order.
  std::vector<element_type> UnitGroupElements() const {
    std::vector<element_type> units;
    const uint64_t n = static_cast<uint64_t>(modulus_);
    for (uint64_t u = 1; u < n; ++u) {
      if (detail::ModRingGcdU64(u, n) == 1) {
        units.push_back(static_cast<element_type>(u));
      }
    }
    return units;
  }

  // Computes Euler's totient phi(n) = |(Z/nZ)^x|.
  uint64_t EulerTotient() const {
    return detail::EulerTotientU64(static_cast<uint64_t>(modulus_));
  }

 private:
  UIntT modulus_;
};

}  // namespace ring

template <typename UIntT>
struct AlgebraTraits<ring::ModularRing<UIntT>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_MODULAR_RING_HPP_
