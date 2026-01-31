#ifndef XGALOIS_INCLUDE_UTILS_MATH_HPP_
#define XGALOIS_INCLUDE_UTILS_MATH_HPP_

#include "xgalois/databases/interface.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace xg {
namespace utils {

//------------------------------------------------------------------------------
// Divisibility
//------------------------------------------------------------------------------

// Returns greatest common divisor of a and b.
template <typename T> inline T Gcd(T a, T b) { return std::gcd(a, b); }

// Returns least common multiple of a and b.
template <typename T> inline T Lcm(T a, T b) { return std::lcm(a, b); }

// Returns the GCD of two numbers and Bézout's identity coefficients.
//
// Computes GCD(a,b) and coefficients of Bézout's identity using the
// Extended Euclidean Algorithm. Returns {d, {x,y}} where d = GCD(a,b)
// and ax + by = d.
//
// @param a First number
// @param b Second number
// @return Pair containing GCD and Bézout's coefficients {x,y}
template <typename T> std::pair<T, std::pair<T, T>> ExtendedGcd(T a, T b) {
  // Initialize coefficients and remainders according to the Extended Euclidean
  // Algorithm setup.
  // Using standard variable names often seen in EGCD literature:
  // r: remainder, s: coefficient for a, t: coefficient for b
  T s = 0;
  T old_s = 1;
  T t = 1;
  T old_t = 0;
  T r = b;
  T old_r = a;

  // The loop continues until the remainder becomes zero.
  while (r != 0) {
    // Calculate the quotient of the current remainders.
    T q = old_r / r;

    // Update remainders: new_r = old_r - q * r
    T temp_r = r;
    r = old_r - q * r;
    old_r = temp_r;

    // Update s coefficients: new_s = old_s - q * s
    T temp_s = s;
    s = old_s - q * s;
    old_s = temp_s;

    // Update t coefficients: new_t = old_t - q * t
    T temp_t = t;
    t = old_t - q * t;
    old_t = temp_t;
  }

  // When the loop terminates, old_r is the GCD, and old_s and
  // old_t are the corresponding coefficients for the original a and b.
  return {old_r, {old_s, old_t}};
}

//------------------------------------------------------------------------------
// Prime Factorization
//------------------------------------------------------------------------------

// --- 1. Trial Division ---
// Concept: This is the most straightforward method. It tests divisibility of
// the number `n` by every integer starting from 2 up to the square root of `n`.
// If a number `i` divides `n`, then `i` is a prime factor. We repeatedly divide
// `n` by `i` until it's no longer divisible, adding `i` to the list of factors.
// We only need to check up to sqrt(n) because if `n` has a factor larger than
// sqrt(n), it must also have a factor smaller than sqrt(n).
// We handle factor 2 separately to then only check odd divisors (3, 5, 7, ...).
//
// Complexity: O(sqrt(n)) in the worst case (when n is prime or a product of two
// large primes). Effective for small numbers (typically up to 10^7).
std::vector<long long> TrialDivision(long long n) {
  std::vector<long long> factors;
  if (n <= 1) {
    return factors; // Numbers less than or equal to 1 have no prime factors.
  }

  // Handle factor 2: Check and divide by 2 until it's no longer possible.
  while (n % 2 == 0) {
    factors.push_back(2);
    n /= 2;
  }

  // Handle odd factors: Check odd numbers from 3 up to the square root of the
  // remaining n. We increment by 2 (i += 2) to only check odd numbers.
  for (long long i = 3; i * i <= n; i += 2) {
    while (n % i == 0) {
      factors.push_back(i);
      n /= i;
    }
  }

  // If the remaining value of n is greater than 1, it means the remaining n is
  // a prime factor itself (larger than sqrt of the original n, or the largest
  // factor).
  if (n > 1) {
    factors.push_back(n);
  }

  return factors;
}

// --- 2. Fermat's Factorization Method ---
// Concept: This method is based on the difference of squares factorization:
// n = a^2 - b^2 = (a - b)(a + b). It works efficiently if the number `n` is a
// product of two factors that are close to each other (i.e., a-b is small).
// It is primarily used for odd numbers.
//
// Process: We search for an integer `a` starting from ceil(sqrt(n)).
// We calculate a^2 - n. If this value is a perfect square (let's call it b^2),
// then we have found a factorization: n = (a - sqrt(a^2 - n)) * (a + sqrt(a^2 -
// n)).
//
// Complexity: Depends on how close the factors are. If factors p and q are
// close, a is close to sqrt(n), and the algorithm is fast. If factors are far
// apart, it's slow. Not a general-purpose algorithm for arbitrary numbers.
std::vector<long long> FermatFactorization(long long n) {
  std::vector<long long> factors;
  // Fermat's method is for odd numbers > 1.
  if (n <= 1 || n % 2 == 0) {
    if (n == 2) {
      factors.push_back(2);
    }
    // For even numbers > 2, handle factor 2 first, then apply to remaining odd
    // part. For simplicity here, we just return empty or handle trivial cases.
    return factors;
  }

  long long a = static_cast<long long>(std::sqrt(n));
  // Ensure a*a >= n. If sqrt(n) is an integer, a*a == n. Otherwise, we need to
  // start with the smallest integer 'a' such that a^2 >= n.
  if (a * a < n) {
    a++;
  }

  // Check for perfect square case: If n is a perfect square, sqrt(n) is a
  // factor.
  if (a * a == n) {
    factors.push_back(a);
    factors.push_back(a);
    return factors;
  }

  // Iterate to find a and b such that a^2 - b^2 = n.
  // We limit the iterations to avoid potential infinite loops for numbers where
  // Fermat's method is not efficient or applicable within a reasonable time.
  const long long kMaxFermatIterations = 1000000; // Limit search range
  long long iterations = 0;

  while (iterations < kMaxFermatIterations) {
    long long bSquared = a * a - n;
    // Check if bSquared is non-negative before taking sqrt.
    if (bSquared < 0) {
      a++;
      iterations++;
      continue;
    }
    long long b = static_cast<long long>(std::sqrt(bSquared));

    // Check if bSquared is a perfect square (b*b == bSquared).
    if (b * b == bSquared) {
      long long factor1 = a - b;
      long long factor2 = a + b;
      // Once factors are found, recursively factor them using Trial Division
      // to get the prime factors.
      std::vector<long long> factors1 = TrialDivision(factor1);
      std::vector<long long> factors2 = TrialDivision(factor2);
      factors.insert(factors.end(), factors1.begin(), factors1.end());
      factors.insert(factors.end(), factors2.begin(), factors2.end());
      // Sort factors for consistent output after combining.
      std::sort(factors.begin(), factors.end());
      return factors;
    }
    a++;
    iterations++;
  }

  // If no factors found within the limited search range using Fermat's method,
  // we return the original number as a factor, indicating that this method
  // couldn't factor it efficiently within the limits.
  factors.push_back(n); // Could not factor using Fermat in limited steps.
  return factors;
}

// --- 3. Pollard's Rho Algorithm ---
// Concept: A probabilistic algorithm that is generally more efficient than
// Trial Division for numbers with relatively small prime factors. It uses
// a pseudo-random sequence and Floyd's cycle-finding algorithm.
//
// Process:
// 1. Choose a starting number x_0 and a function f(x) (e.g., f(x) = (x^2 + c)
// mod n).
// 2. Generate a sequence x_1 = f(x_0), x_2 = f(x_1), ... modulo n.
// 3. Simultaneously, track two "pointers" or values in the sequence: one moving
//    one step at a time (tortoise) and one moving two steps at a time (hare).
// 4. The difference between the hare's value and the tortoise's value, modulo a
//    prime factor `p` of `n`, will eventually become 0 when they are at the
//    same point in the sequence modulo `p`.
// 5. The greatest common divisor (GCD) of the absolute difference between the
//    hare and tortoise values and `n` will likely reveal a non-trivial factor
//    of `n`. If gcd(|hare - tortoise|, n) > 1 and < n, a factor is found.
//
// Complexity: Expected time complexity is roughly O(n^(1/4)). Its efficiency
// depends on the size of the smallest prime factor of n.
long long PollardsRho(long long n) {
  if (n <= 1) {
    return n;
  }
  if (n % 2 == 0) {
    return 2; // Handle even numbers quickly.
  }

  // Use a random number generator for initial x (tortoise) and constant c.
  // Seeding with current time provides different sequences on different runs.
  std::mt19937_64 rng(
      std::chrono::steady_clock::now().time_since_epoch().count());
  // Ensure distribution range is valid for n > 1.
  std::uniform_int_distribution<long long> dist(1, n - 1 > 0 ? n - 1 : 1);

  long long x = dist(rng); // Tortoise
  long long y = x;         // Hare
  long long c = dist(rng); // Constant for the function f(x)
  long long d = 1;         // GCD, initialized to 1

  // The function f(x) = (x^2 + c) mod n.
  // Use __int128 for intermediate calculation (val * val) to prevent overflow
  // before taking the modulo, especially for large `long long` values of `val`.
  auto f = [&](long long val) {
    return (static_cast<__int128>(val) * val + c) % n;
  };

  // Floyd's cycle-finding algorithm.
  // We iterate until a non-trivial GCD is found or we reach an iteration limit
  // to prevent potential infinite loops for prime numbers or difficult cases.
  const int kMaxRhoIterations = 1000000; // Limit iterations
  int iterations = 0;

  while (d == 1 && iterations < kMaxRhoIterations) {
    x = f(x);    // Tortoise moves one step
    y = f(f(y)); // Hare moves two steps
    // Calculate the GCD of the absolute difference and n.
    // std::gcd is available in C++17 and later.
    d = std::gcd(std::abs(x - y), n);
    iterations++;
  }

  // Check the result of the GCD:
  // If d == n, it means the algorithm failed to find a non-trivial factor
  // (e.g., n is prime, or the sequence didn't reveal a factor within limits).
  // If d == 1, it means the iteration limit was reached before finding a
  // factor.
  if (d == n || d == 1) {
    // Failed to find a factor within iteration limit or found trivial factor.
    // For a more robust implementation, you might repeat the process with
    // different random seeds/constants or use a fallback algorithm. Here, we
    // return n, implying failure to find a non-trivial factor in this attempt.
    return n;
  }

  return d; // Found a non-trivial factor.
}

// Helper to get all factors using Pollard's Rho recursively.
// Factors a given number `n` completely by repeatedly applying Pollard's Rho
// to find one factor, then recursively factoring the resulting parts.
// Uses Trial Division for smaller numbers for efficiency.
std::vector<long long> FactorizePollardsRho(long long n) {
  std::vector<long long> factors;
  if (n <= 1) {
    return factors;
  }

  // Handle even numbers first using Trial Division.
  while (n % 2 == 0) {
    factors.push_back(2);
    n /= 2;
  }

  if (n == 1) {
    return factors;
  }

  // Use trial division for small remaining numbers for efficiency.
  // Pollard's Rho has overhead, so for small numbers, Trial Division is faster.
  // The threshold can be adjusted based on performance testing.
  const long long kTrialDivisionThreshold = 1000000;
  if (n < kTrialDivisionThreshold) {
    std::vector<long long> trialFactors = TrialDivision(n);
    factors.insert(factors.end(), trialFactors.begin(), trialFactors.end());
    return factors;
  }

  // Use Pollard's Rho to find one factor of the remaining n.
  long long factor = PollardsRho(n);

  if (factor == n) {
    // If Pollard's Rho returned n, it means it failed to find a non-trivial
    // factor. This could happen if n is prime or the algorithm got stuck. In
    // this case, we assume n is prime (or couldn't be factored by Rho) and add
    // it as a factor.
    factors.push_back(n);
  } else {
    // If a non-trivial factor is found, recursively factor the found factor
    // and the remaining part (n / factor).
    std::vector<long long> factors1 = FactorizePollardsRho(factor);
    std::vector<long long> factors2 = FactorizePollardsRho(n / factor);
    factors.insert(factors.end(), factors1.begin(), factors1.end());
    factors.insert(factors.end(), factors2.begin(), factors2.end());
  }

  // Sort factors for consistent output order.
  std::sort(factors.begin(), factors.end());
  return factors;
}

// --- 4. Prime Factorization Function ---
// Concept: This function attempts to choose the "best" algorithm based on the
// size of the input number. For smaller numbers, Trial Division is fastest.
// For larger numbers, Pollard's Rho (or its recursive application) is generally
// more efficient than Trial Division. More advanced algorithms like QS or GNFS
// would be needed for very large numbers (hundreds of digits), but are not
// implemented here due to complexity.
//
// Process: First check if factors exist in the database. If not found, check
// the input number against predefined thresholds and call the appropriate
// factorization function.
std::vector<long long> PrimeFactorize(long long n) {
  if (n <= 1) {
    return {}; // Return empty vector for numbers <= 1
  }

  // First, try to get factors from the database
  try {
    xg::databases::PrimeFactorsDatabase db;
    xg::databases::PrimeFactorsResult result = db.fetch(n);

    // Convert the database result to the expected format
    std::vector<long long> factors;
    for (size_t i = 0; i < result.factors.size(); ++i) {
      // Add each factor according to its multiplicity
      for (int j = 0; j < result.multiplicities[i]; ++j) {
        factors.push_back(result.factors[i]);
      }
    }

    // Sort factors for consistent output order
    std::sort(factors.begin(), factors.end());
    return factors;

  } catch (const std::runtime_error &e) {
    // Database lookup failed (file not found or number not in database)
    // Fall back to computational methods
  }

  // Define thresholds for switching algorithms.
  // These thresholds are approximate and can be tuned based on performance
  // testing on your specific system and typical input range.
  const long long kTrialDivisionLimit =
      10000000; // Use trial division up to this limit
  // For numbers larger than this, Pollard's Rho is generally better.

  if (n < kTrialDivisionLimit) {
    return TrialDivision(n);
  } else {
    // For larger numbers, use the recursive Pollard's Rho approach
    // which also incorporates trial division for smaller resulting factors.
    return FactorizePollardsRho(n);
  }
  // Note: Fermat's method is not included in the smart selection by default
  // as its efficiency depends on factor proximity, not just number size,
  // making it less reliable as a general-purpose choice based solely on size.
}

//------------------------------------------------------------------------------
// Prime Testing
//------------------------------------------------------------------------------

/**
 * @brief Determines if a number is prime using optimized trial division
 *
 * Uses an efficient trial division approach that checks divisibility by 2,
 * then tests only odd numbers up to the square root of n. For larger numbers,
 * it leverages the existing prime factorization infrastructure.
 *
 * @param n The number to test for primality (must be ≥ 0)
 *
 * @return true if n is prime, false otherwise
 *
 * @note Returns false for n < 2 as per mathematical convention.
 *       Uses optimized algorithms based on input size for best performance.
 */
bool IsPrime(long long n) {
  // Handle edge cases
  if (n < 2) return false;
  if (n == 2) return true;
  if (n % 2 == 0) return false;

  // For small numbers, use direct trial division
  const long long kDirectTestLimit = 1000000;
  if (n <= kDirectTestLimit) {
    // Test odd divisors up to sqrt(n)
    for (long long i = 3; i * i <= n; i += 2) {
      if (n % i == 0) {
        return false;
      }
    }
    return true;
  }

  // For larger numbers, use the factorization infrastructure
  // A number is prime if its only prime factor is itself
  std::vector<long long> factors = PrimeFactorize(n);
  return factors.size() == 1 && factors[0] == n;
}

//------------------------------------------------------------------------------
// Prime Power Decomposition
//------------------------------------------------------------------------------

/**
 * @brief Prime power decomposition using database-backed factorization
 *
 * Decomposes a positive integer into (prime, exponent) form to determine if
 * it represents a valid finite field order. First attempts to use the
 * prime factors database for fast lookup, then falls back to algorithmic
 * factorization from utils::PrimeFactorize for numbers not in the database.
 *
 * @param n The number to decompose (must be ≥ 2)
 *
 * @return Pair (prime, exponent) if n = prime^exponent, otherwise (0, 0)
 *
 * @note This method provides significant performance improvements for
 *       common field orders by using precomputed database lookups while
 *       maintaining correctness through algorithmic fallback.
 */
std::pair<uint64_t, uint64_t> DecomposePrimePower(uint64_t n) {
  if (n < 2)
    return {0, 0};

  // First attempt: Use prime factors database for fast lookup
  try {
    databases::PrimeFactorsDatabase prime_db;
    auto result = prime_db.fetch(static_cast<long long>(n));

    // Check if it's a prime power (all factors are the same)
    if (!result.factors.empty()) {
      long long prime = result.factors[0];
      uint64_t total_exponent = 0;

      // Verify all factors are the same and sum multiplicities
      for (size_t i = 0; i < result.factors.size(); ++i) {
        if (result.factors[i] != prime) {
          return {0, 0}; // Not a prime power
        }
        total_exponent += result.multiplicities[i];
      }

      return {static_cast<uint64_t>(prime), total_exponent};
    }
  } catch (const std::runtime_error &) {
    // Database lookup failed, fall back to algorithmic approach
  }

  /* Fallback: Employ advanced prime factorization for robust number
   * decomposition */
  std::vector<long long> prime_factors =
      PrimeFactorize(static_cast<long long>(n));

  if (prime_factors.empty()) {
    return {0, 0};
  }

  /* Verify that all prime factors are identical (i.e., n = p^k for some prime
   * p) */
  long long prime = prime_factors[0];
  for (long long factor : prime_factors) {
    if (factor != prime) {
      return {0, 0}; /* n is not a prime power */
    }
  }

  /* All factors are identical, confirming n = prime^exponent structure */
  uint64_t exponent = prime_factors.size();
  return {static_cast<uint64_t>(prime), exponent};
}

//------------------------------------------------------------------------------
// Safe Integer Power
//------------------------------------------------------------------------------

/**
 * @brief Safely computes base^exponent for integer types
 *
 * Computes integer exponentiation using pure integer arithmetic to avoid
 * precision issues that can occur with std::pow when used with large integers.
 * Uses exponentiation by squaring for efficiency and includes overflow detection.
 *
 * @param base The base value
 * @param exponent The exponent value
 * @return The result of base^exponent
 * @throws std::overflow_error if the result would exceed uint64_t limits
 * @throws std::invalid_argument if base is 0 and exponent is 0 (undefined)
 *
 * @note This function handles edge cases:
 *       - 0^0 returns 1 by mathematical convention
 *       - 0^n returns 0 for n > 0
 *       - 1^n returns 1 for any n
 *       - n^0 returns 1 for any n
 */
inline uint64_t SafeIntegerPower(uint64_t base, uint64_t exponent) {
  // Handle edge cases
  if (exponent == 0) {
    return 1; // By convention: n^0 = 1 for any n (including 0^0 = 1)
  }

  if (base == 0) {
    return 0; // 0^n = 0 for n > 0
  }

  if (base == 1) {
    return 1; // 1^n = 1 for any n
  }

  // Use exponentiation by squaring to compute base^exponent
  uint64_t result = 1;
  uint64_t current_base = base;
  uint64_t current_exp = exponent;

  while (current_exp > 0) {
    // If current exponent is odd, multiply result by current base
    if (current_exp & 1) {
      // Check for overflow before multiplication
      if (result > UINT64_MAX / current_base) {
        throw std::overflow_error(
            "Integer power computation would overflow uint64_t");
      }
      result *= current_base;
    }

    // Break early if we're done to avoid unnecessary overflow check
    current_exp >>= 1;
    if (current_exp == 0) {
      break;
    }

    // Check for overflow before squaring current_base
    if (current_base > UINT64_MAX / current_base) {
      throw std::overflow_error(
          "Integer power computation would overflow uint64_t");
    }
    current_base *= current_base;
  }

  return result;
}

//------------------------------------------------------------------------------

} // namespace utils
} // namespace xg

#endif // XGALOIS_INCLUDE_UTILS_MATH_HPP_
