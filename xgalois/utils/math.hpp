#ifndef XGALOIS_UTILS_MATH_HPP_
#define XGALOIS_UTILS_MATH_HPP_

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace xg {
namespace utils {

namespace detail {

__extension__ using u128 = unsigned __int128;
__extension__ using i128 = __int128;

}  // namespace detail

template <typename T>
inline T Gcd(T a, T b) {
  return std::gcd(a, b);
}

template <typename T>
inline T Lcm(T a, T b) {
  return std::lcm(a, b);
}

template <typename T>
inline std::pair<T, std::pair<T, T>> ExtendedGcd(T a, T b) {
  if constexpr (std::is_integral_v<T> && sizeof(T) <= 8) {
    detail::i128 old_r = static_cast<detail::i128>(a);
    detail::i128 r = static_cast<detail::i128>(b);
    detail::i128 old_s = 1;
    detail::i128 s = 0;
    detail::i128 old_t = 0;
    detail::i128 t = 1;

    while (r != 0) {
      detail::i128 q = old_r / r;

      detail::i128 temp_r = r;
      r = old_r - q * r;
      old_r = temp_r;

      detail::i128 temp_s = s;
      s = old_s - q * s;
      old_s = temp_s;

      detail::i128 temp_t = t;
      t = old_t - q * t;
      old_t = temp_t;
    }

    if constexpr (std::is_signed_v<T>) {
      if (old_r < 0) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
      }
    }

    return {static_cast<T>(old_r),
            {static_cast<T>(old_s), static_cast<T>(old_t)}};
  } else {
    T s = 0;
    T old_s = 1;
    T t = 1;
    T old_t = 0;
    T r = b;
    T old_r = a;

    while (r != 0) {
      T q = old_r / r;

      T temp_r = r;
      r = old_r - q * r;
      old_r = temp_r;

      T temp_s = s;
      s = old_s - q * s;
      old_s = temp_s;

      T temp_t = t;
      t = old_t - q * t;
      old_t = temp_t;
    }

    if constexpr (std::is_signed_v<T>) {
      if (old_r < 0) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
      }
    }

    return {old_r, {old_s, old_t}};
  }
}

inline uint64_t ModMul(uint64_t a, uint64_t b, uint64_t mod) {
  if (mod == 0) {
    throw std::domain_error("Modulus must be nonzero in ModMul");
  }
  if (mod == 1) {
    return 0;
  }
  return static_cast<uint64_t>(
      (static_cast<detail::u128>(a) * static_cast<detail::u128>(b)) %
      static_cast<detail::u128>(mod));
}

inline uint64_t ModPow(uint64_t base, uint64_t exp, uint64_t mod) {
  if (mod == 0) {
    throw std::domain_error("Modulus must be nonzero in ModPow");
  }
  if (mod == 1) {
    return 0;
  }
  uint64_t result = 1;
  base %= mod;
  while (exp > 0) {
    if (exp & 1ULL) {
      result = ModMul(result, base, mod);
    }
    base = ModMul(base, base, mod);
    exp >>= 1;
  }
  return result;
}

inline uint64_t ModInv(uint64_t a, uint64_t mod) {
  if (mod <= 1) {
    throw std::domain_error("Modulus must be greater than 1 in ModInv");
  }
  a %= mod;
  if (a == 0) {
    throw std::domain_error("Modular inverse does not exist (gcd != 1)");
  }

  detail::i128 old_r = static_cast<detail::i128>(a);
  detail::i128 r = static_cast<detail::i128>(mod);
  detail::i128 old_s = 1;
  detail::i128 s = 0;

  while (r != 0) {
    detail::i128 q = old_r / r;

    detail::i128 temp_r = r;
    r = old_r - q * r;
    old_r = temp_r;

    detail::i128 temp_s = s;
    s = old_s - q * s;
    old_s = temp_s;
  }

  if (old_r != 1) {
    throw std::domain_error("Modular inverse does not exist (gcd != 1)");
  }

  detail::i128 mod_128 = static_cast<detail::i128>(mod);
  old_s %= mod_128;
  if (old_s < 0) {
    old_s += mod_128;
  }
  return static_cast<uint64_t>(old_s);
}

inline uint64_t ISqrt(uint64_t n) {
  if (n < 2) {
    return n;
  }
  unsigned shift = (static_cast<unsigned>(std::bit_width(n)) + 1U) / 2U;
  uint64_t x = 1ULL << shift;
  while (true) {
    uint64_t y = (x + n / x) >> 1;
    if (y >= x) {
      return x;
    }
    x = y;
  }
}

namespace detail {

inline uint64_t ISqrt128(u128 n) {
  if (n <= std::numeric_limits<uint64_t>::max()) {
    return ISqrt(static_cast<uint64_t>(n));
  }
  uint64_t hi = static_cast<uint64_t>(n >> 64);
  unsigned bw = 64U + static_cast<unsigned>(std::bit_width(hi));
  unsigned shift = (bw + 1U) / 2U;
  u128 x = static_cast<u128>(1) << shift;
  while (true) {
    u128 y = (x + n / x) >> 1;
    if (y >= x) {
      return static_cast<uint64_t>(x);
    }
    x = y;
  }
}

}  // namespace detail

inline bool IsPrime(uint64_t n) {
  if (n < 2) {
    return false;
  }

  constexpr uint64_t kSmallPrimes[] = {2,  3,  5,  7,  11, 13,
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

  // Deterministic 64-bit Miller-Rabin witness set (Sinclair 2011).
  // Proven 100% exact for all uint64_t values (n < 2^64).
  constexpr uint64_t kWitnesses[] = {
      2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};

  uint64_t d = n - 1;
  int s = std::countr_zero(d);
  d >>= s;

  for (uint64_t a : kWitnesses) {
    uint64_t a_mod = a % n;
    if (a_mod == 0) {
      continue;
    }
    uint64_t x = ModPow(a_mod, d, n);
    if (x == 1 || x == n - 1) {
      continue;
    }
    bool composite = true;
    for (int r = 1; r < s; ++r) {
      x = ModMul(x, x, n);
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

inline std::vector<uint64_t> TrialDivision(uint64_t n) {
  std::vector<uint64_t> factors;
  if (n <= 1) {
    return factors;
  }

  while ((n & 1ULL) == 0) {
    factors.push_back(2);
    n >>= 1;
  }

  for (uint64_t i = 3; i <= n / i; i += 2) {
    while (n % i == 0) {
      factors.push_back(i);
      n /= i;
    }
  }

  if (n > 1) {
    factors.push_back(n);
  }

  return factors;
}

inline uint64_t PollardsRho(uint64_t n) {
  if (n <= 1) {
    return n;
  }
  if ((n & 1ULL) == 0) {
    return 2;
  }

  constexpr uint64_t kSmallOddPrimes[] = {3,  5,  7,  11, 13, 17,
                                          19, 23, 29, 31, 37};
  for (uint64_t p : kSmallOddPrimes) {
    if (n == p) {
      return n;
    }
    if (n % p == 0) {
      return p;
    }
  }

  if (IsPrime(n)) {
    return n;
  }

  auto f = [n](uint64_t val, uint64_t c) -> uint64_t {
    return static_cast<uint64_t>(
        (static_cast<detail::u128>(val) * static_cast<detail::u128>(val) +
         static_cast<detail::u128>(c)) %
        static_cast<detail::u128>(n));
  };

  constexpr uint64_t kBatchSize = 128;

  // Deterministic Brent's cycle-finding Pollard's Rho with batched GCD.
  for (uint64_t c = 1; c < n; ++c) {
    uint64_t y = (c + 1) % n;
    if (y == 0) {
      y = 1;
    }
    uint64_t r = 1;
    uint64_t q = 1;
    uint64_t g = 1;
    uint64_t x = y;
    uint64_t ys = y;

    while (g == 1 && r <= (1ULL << 32)) {
      x = y;
      for (uint64_t i = 0; i < r; ++i) {
        y = f(y, c);
      }
      uint64_t k = 0;
      while (k < r && g == 1) {
        ys = y;
        uint64_t limit = std::min(kBatchSize, r - k);
        for (uint64_t i = 0; i < limit; ++i) {
          y = f(y, c);
          uint64_t diff = x > y ? x - y : y - x;
          q = ModMul(q, diff, n);
        }
        g = std::gcd(q, n);
        k += limit;
      }
      r <<= 1;
    }

    if (g == n) {
      do {
        ys = f(ys, c);
        uint64_t diff = x > ys ? x - ys : ys - x;
        g = std::gcd(diff, n);
      } while (g == 1);
    }

    if (g > 1 && g < n) {
      return g;
    }
  }

  for (uint64_t i = 41; i <= n / i; i += 2) {
    if (n % i == 0) {
      return i;
    }
  }
  return n;
}

inline std::vector<uint64_t> FactorizePollardsRho(uint64_t n) {
  std::vector<uint64_t> factors;
  if (n <= 1) {
    return factors;
  }

  constexpr uint64_t kSmallPrimes[] = {2,  3,  5,  7,  11, 13,
                                       17, 19, 23, 29, 31, 37};
  for (uint64_t p : kSmallPrimes) {
    while (n % p == 0) {
      factors.push_back(p);
      n /= p;
    }
  }

  if (n == 1) {
    return factors;
  }

  std::vector<uint64_t> stack;
  stack.push_back(n);
  while (!stack.empty()) {
    uint64_t cur = stack.back();
    stack.pop_back();
    if (cur <= 1) {
      continue;
    }
    if (IsPrime(cur)) {
      factors.push_back(cur);
      continue;
    }
    uint64_t d = PollardsRho(cur);
    stack.push_back(d);
    stack.push_back(cur / d);
  }

  std::sort(factors.begin(), factors.end());
  return factors;
}

namespace detail {

inline void FermatFactorizeOdd(uint64_t m, std::vector<uint64_t>& factors) {
  if (m <= 1) {
    return;
  }
  if (IsPrime(m)) {
    factors.push_back(m);
    return;
  }

  uint64_t a = ISqrt(m);
  if (static_cast<u128>(a) * static_cast<u128>(a) < static_cast<u128>(m)) {
    ++a;
  }

  if (static_cast<u128>(a) * static_cast<u128>(a) == static_cast<u128>(m)) {
    FermatFactorizeOdd(a, factors);
    FermatFactorizeOdd(a, factors);
    return;
  }

  constexpr uint64_t kMaxFermatIterations = 1000000;
  for (uint64_t iter = 0; iter < kMaxFermatIterations; ++iter) {
    u128 b_squared =
        static_cast<u128>(a) * static_cast<u128>(a) - static_cast<u128>(m);
    uint64_t b = ISqrt128(b_squared);
    if (static_cast<u128>(b) * static_cast<u128>(b) == b_squared) {
      uint64_t factor1 = a - b;
      uint64_t factor2 = a + b;
      if (factor1 > 1 && factor2 > 1) {
        FermatFactorizeOdd(factor1, factors);
        FermatFactorizeOdd(factor2, factors);
        return;
      }
      break;
    }
    ++a;
  }

  // Fallback cleanly to Pollard's Rho if Fermat iteration budget exceeded.
  std::vector<uint64_t> fallback = FactorizePollardsRho(m);
  factors.insert(factors.end(), fallback.begin(), fallback.end());
}

}  // namespace detail

inline std::vector<uint64_t> FermatFactorization(uint64_t n) {
  std::vector<uint64_t> factors;
  if (n <= 1) {
    return factors;
  }

  while ((n & 1ULL) == 0) {
    factors.push_back(2);
    n >>= 1;
  }

  if (n > 1) {
    detail::FermatFactorizeOdd(n, factors);
  }

  std::sort(factors.begin(), factors.end());
  return factors;
}

inline std::vector<uint64_t> PrimeFactorize(uint64_t n) {
  std::vector<uint64_t> factors;
  if (n <= 1) {
    return factors;
  }

  // Fast small-prime trial division up to 97.
  constexpr uint64_t kSmallPrimes[] = {2,  3,  5,  7,  11, 13, 17, 19, 23,
                                       29, 31, 37, 41, 43, 47, 53, 59, 61,
                                       67, 71, 73, 79, 83, 89, 97};
  for (uint64_t p : kSmallPrimes) {
    while (n % p == 0) {
      factors.push_back(p);
      n /= p;
    }
    if (n == 1) {
      return factors;
    }
  }

  if (n < 101ULL * 101ULL || IsPrime(n)) {
    factors.push_back(n);
    return factors;
  }

  std::vector<uint64_t> remaining = FactorizePollardsRho(n);
  factors.insert(factors.end(), remaining.begin(), remaining.end());
  std::sort(factors.begin(), factors.end());
  return factors;
}

inline std::vector<std::pair<uint64_t, uint64_t>>
PrimeFactorizeWithMultiplicity(uint64_t n) {
  std::vector<std::pair<uint64_t, uint64_t>> result;
  std::vector<uint64_t> factors = PrimeFactorize(n);
  for (uint64_t f : factors) {
    if (!result.empty() && result.back().first == f) {
      ++result.back().second;
    } else {
      result.emplace_back(f, 1ULL);
    }
  }
  return result;
}

inline std::vector<uint64_t> DistinctPrimeFactors(uint64_t n) {
  std::vector<uint64_t> distinct;
  std::vector<uint64_t> factors = PrimeFactorize(n);
  for (uint64_t f : factors) {
    if (distinct.empty() || distinct.back() != f) {
      distinct.push_back(f);
    }
  }
  return distinct;
}

inline uint64_t EulerTotient(uint64_t n) {
  if (n == 0) {
    return 0;
  }
  if (n == 1) {
    return 1;
  }
  uint64_t result = n;
  for (uint64_t p : DistinctPrimeFactors(n)) {
    result = (result / p) * (p - 1);
  }
  return result;
}

inline std::pair<uint64_t, uint64_t> DecomposePrimePower(uint64_t n) {
  if (n < 2) {
    return {0, 0};
  }
  if ((n & (n - 1)) == 0) {
    return {2ULL, static_cast<uint64_t>(std::countr_zero(n))};
  }
  if (IsPrime(n)) {
    return {n, 1ULL};
  }

  auto factor_pairs = PrimeFactorizeWithMultiplicity(n);
  if (factor_pairs.size() != 1) {
    return {0, 0};
  }
  return factor_pairs[0];
}

inline uint64_t SafeIntegerPower(uint64_t base, uint64_t exponent) {
  if (exponent == 0) {
    return 1;
  }

  if (base == 0) {
    return 0;
  }

  if (base == 1) {
    return 1;
  }

  uint64_t result = 1;
  uint64_t current_base = base;
  uint64_t current_exp = exponent;

  while (current_exp > 0) {
    if (current_exp & 1ULL) {
      if (result > std::numeric_limits<uint64_t>::max() / current_base) {
        throw std::overflow_error(
            "Integer power computation would overflow uint64_t");
      }
      result *= current_base;
    }

    current_exp >>= 1;
    if (current_exp == 0) {
      break;
    }

    if (current_base > std::numeric_limits<uint64_t>::max() / current_base) {
      throw std::overflow_error(
          "Integer power computation would overflow uint64_t");
    }
    current_base *= current_base;
  }

  return result;
}

}  // namespace utils
}  // namespace xg

#endif  // XGALOIS_UTILS_MATH_HPP_
