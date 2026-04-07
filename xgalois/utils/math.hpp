#ifndef XGALOIS_INCLUDE_UTILS_MATH_HPP_
#define XGALOIS_INCLUDE_UTILS_MATH_HPP_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "xgalois/databases/interface.hpp"

namespace xg {
namespace utils {

template <typename T>
inline T Gcd(T a, T b) {
  return std::gcd(a, b);
}

template <typename T>
inline T Lcm(T a, T b) {
  return std::lcm(a, b);
}

template <typename T>
std::pair<T, std::pair<T, T>> ExtendedGcd(T a, T b) {

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

  return {old_r, {old_s, old_t}};
}

std::vector<uint64_t> TrialDivision(uint64_t n) {
  std::vector<uint64_t> factors;
  if (n <= 1) {
    return factors;
  }

  while (n % 2 == 0) {
    factors.push_back(2);
    n /= 2;
  }

  for (uint64_t i = 3; i * i <= n; i += 2) {
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

std::vector<uint64_t> FermatFactorization(uint64_t n) {
  std::vector<uint64_t> factors;

  if (n <= 1 || n % 2 == 0) {
    if (n == 2) {
      factors.push_back(2);
    }

    return factors;
  }

  uint64_t a = static_cast<uint64_t>(std::sqrt(n));

  if (a * a < n) {
    a++;
  }

  if (a * a == n) {
    factors.push_back(a);
    factors.push_back(a);
    return factors;
  }

  const uint64_t kMaxFermatIterations = 1000000;
  uint64_t iterations = 0;

  while (iterations < kMaxFermatIterations) {
    uint64_t bSquared = a * a - n;

    if (bSquared < 0) {
      a++;
      iterations++;
      continue;
    }
    uint64_t b = static_cast<uint64_t>(std::sqrt(bSquared));

    if (b * b == bSquared) {
      uint64_t factor1 = a - b;
      uint64_t factor2 = a + b;

      std::vector<uint64_t> factors1 = TrialDivision(factor1);
      std::vector<uint64_t> factors2 = TrialDivision(factor2);
      factors.insert(factors.end(), factors1.begin(), factors1.end());
      factors.insert(factors.end(), factors2.begin(), factors2.end());

      std::sort(factors.begin(), factors.end());
      return factors;
    }
    a++;
    iterations++;
  }

  factors.push_back(n);
  return factors;
}

uint64_t PollardsRho(uint64_t n) {
  if (n <= 1) {
    return n;
  }
  if (n % 2 == 0) {
    return 2;
  }

  std::mt19937_64 rng(
      std::chrono::steady_clock::now().time_since_epoch().count());

  std::uniform_int_distribution<uint64_t> dist(1, n - 1 > 0 ? n - 1 : 1);

  uint64_t x = dist(rng);
  uint64_t y = x;
  uint64_t c = dist(rng);
  uint64_t d = 1;

  auto f = [&](uint64_t val) {
    return (static_cast<__int128>(val) * val + c) % n;
  };

  const int kMaxRhoIterations = 1000000;
  int iterations = 0;

  while (d == 1 && iterations < kMaxRhoIterations) {
    x = f(x);
    y = f(f(y));

    uint64_t diff = x > y ? x - y : y - x;
    d = std::gcd(diff, n);
    iterations++;
  }

  if (d == n || d == 1) {

    return n;
  }

  return d;
}

std::vector<uint64_t> FactorizePollardsRho(uint64_t n) {
  std::vector<uint64_t> factors;
  if (n <= 1) {
    return factors;
  }

  while (n % 2 == 0) {
    factors.push_back(2);
    n /= 2;
  }

  if (n == 1) {
    return factors;
  }

  const uint64_t kTrialDivisionThreshold = 1000000;
  if (n < kTrialDivisionThreshold) {
    std::vector<uint64_t> trialFactors = TrialDivision(n);
    factors.insert(factors.end(), trialFactors.begin(), trialFactors.end());
    return factors;
  }

  uint64_t factor = PollardsRho(n);

  if (factor == n) {

    factors.push_back(n);
  } else {

    std::vector<uint64_t> factors1 = FactorizePollardsRho(factor);
    std::vector<uint64_t> factors2 = FactorizePollardsRho(n / factor);
    factors.insert(factors.end(), factors1.begin(), factors1.end());
    factors.insert(factors.end(), factors2.begin(), factors2.end());
  }

  std::sort(factors.begin(), factors.end());
  return factors;
}

std::vector<uint64_t> PrimeFactorize(uint64_t n) {
  if (n <= 1) {
    return {};
  }

  try {
    xg::databases::PrimeFactorsDatabase db;
    xg::databases::PrimeFactorsResult result = db.fetch(n);

    std::vector<uint64_t> factors;
    for (size_t i = 0; i < result.factors.size(); ++i) {

      for (int j = 0; j < result.multiplicities[i]; ++j) {
        factors.push_back(result.factors[i]);
      }
    }

    std::sort(factors.begin(), factors.end());
    return factors;

  } catch (const std::runtime_error &e) {

  }

  const uint64_t kTrialDivisionLimit =
      10000000;

  if (n < kTrialDivisionLimit) {
    return TrialDivision(n);
  } else {

    return FactorizePollardsRho(n);
  }

}

bool IsPrime(uint64_t n) {

  if (n < 2) return false;
  if (n == 2) return true;
  if (n % 2 == 0) return false;

  const uint64_t kDirectTestLimit = 1000000;
  if (n <= kDirectTestLimit) {

    for (uint64_t i = 3; i * i <= n; i += 2) {
      if (n % i == 0) {
        return false;
      }
    }
    return true;
  }

  std::vector<uint64_t> factors = PrimeFactorize(n);
  return factors.size() == 1 && factors[0] == n;
}

std::pair<uint64_t, uint64_t> DecomposePrimePower(uint64_t n) {
  if (n < 2) return {0, 0};

  try {
    databases::PrimeFactorsDatabase prime_db;
    auto result = prime_db.fetch(static_cast<uint64_t>(n));

    if (!result.factors.empty()) {
      uint64_t prime = result.factors[0];
      uint64_t total_exponent = 0;

      for (size_t i = 0; i < result.factors.size(); ++i) {
        if (result.factors[i] != prime) {
          return {0, 0};
        }
        total_exponent += result.multiplicities[i];
      }

      return {static_cast<uint64_t>(prime), total_exponent};
    }
  } catch (const std::runtime_error &) {

  }

  std::vector<uint64_t> prime_factors =
      PrimeFactorize(static_cast<uint64_t>(n));

  if (prime_factors.empty()) {
    return {0, 0};
  }

  uint64_t prime = prime_factors[0];
  for (uint64_t factor : prime_factors) {
    if (factor != prime) {
      return {0, 0};
    }
  }

  uint64_t exponent = prime_factors.size();
  return {static_cast<uint64_t>(prime), exponent};
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

    if (current_exp & 1) {

      if (result > UINT64_MAX / current_base) {
        throw std::overflow_error(
            "Integer power computation would overflow uint64_t");
      }
      result *= current_base;
    }

    current_exp >>= 1;
    if (current_exp == 0) {
      break;
    }

    if (current_base > UINT64_MAX / current_base) {
      throw std::overflow_error(
          "Integer power computation would overflow uint64_t");
    }
    current_base *= current_base;
  }

  return result;
}

}
}

#endif
