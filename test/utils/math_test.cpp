#include "xgalois/utils/math.hpp"
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

using namespace xg::utils;

//===----------------------------------------------------------------------===//
// Test Fixtures
//===----------------------------------------------------------------------===//

class MathUtilsTest : public ::testing::Test {
protected:
  void SetUp() override {
    // No special setup needed for math utilities
  }
};

//===----------------------------------------------------------------------===//
// Divisibility Tests
//===----------------------------------------------------------------------===//

TEST_F(MathUtilsTest, GcdBasicCases) {
  EXPECT_EQ(Gcd(12, 8), 4);
  EXPECT_EQ(Gcd(21, 14), 7);
  EXPECT_EQ(Gcd(17, 13), 1); // Coprime numbers
  EXPECT_EQ(Gcd(0, 5), 5);
  EXPECT_EQ(Gcd(5, 0), 5);
  EXPECT_EQ(Gcd(0, 0), 0);
}

TEST_F(MathUtilsTest, GcdLargeNumbers) {
  EXPECT_EQ(Gcd(1071, 462), 21);
  EXPECT_EQ(Gcd(48, 18), 6);
  EXPECT_EQ(Gcd(252, 105), 21);
}

TEST_F(MathUtilsTest, LcmBasicCases) {
  EXPECT_EQ(Lcm(4, 6), 12);
  EXPECT_EQ(Lcm(21, 14), 42);
  EXPECT_EQ(Lcm(17, 13), 221); // Coprime numbers
  EXPECT_EQ(Lcm(12, 8), 24);
}

TEST_F(MathUtilsTest, ExtendedGcdBasicCases) {
  auto result = ExtendedGcd(240, 46);
  EXPECT_EQ(result.first, 2); // GCD(240, 46) = 2

  // Verify Bézout's identity: 240 * x + 46 * y = gcd
  long long x = result.second.first;
  long long y = result.second.second;
  EXPECT_EQ(240 * x + 46 * y, result.first);
}

TEST_F(MathUtilsTest, ExtendedGcdCoprimeNumbers) {
  auto result = ExtendedGcd(17, 13);
  EXPECT_EQ(result.first, 1); // GCD(17, 13) = 1

  // Verify Bézout's identity
  long long x = result.second.first;
  long long y = result.second.second;
  EXPECT_EQ(17 * x + 13 * y, 1);
}

//===----------------------------------------------------------------------===//
// Prime Testing Tests
//===----------------------------------------------------------------------===//

TEST_F(MathUtilsTest, IsPrimeSmallNumbers) {
  // Test edge cases
  EXPECT_FALSE(IsPrime(-1));
  EXPECT_FALSE(IsPrime(0));
  EXPECT_FALSE(IsPrime(1));

  // Test first few primes
  EXPECT_TRUE(IsPrime(2));
  EXPECT_TRUE(IsPrime(3));
  EXPECT_FALSE(IsPrime(4));
  EXPECT_TRUE(IsPrime(5));
  EXPECT_FALSE(IsPrime(6));
  EXPECT_TRUE(IsPrime(7));
  EXPECT_FALSE(IsPrime(8));
  EXPECT_FALSE(IsPrime(9));
  EXPECT_FALSE(IsPrime(10));
  EXPECT_TRUE(IsPrime(11));
  EXPECT_FALSE(IsPrime(12));
  EXPECT_TRUE(IsPrime(13));
}

TEST_F(MathUtilsTest, IsPrimeMediumNumbers) {
  // Test some medium-sized primes
  EXPECT_TRUE(IsPrime(97));
  EXPECT_TRUE(IsPrime(101));
  EXPECT_TRUE(IsPrime(103));
  EXPECT_TRUE(IsPrime(107));
  EXPECT_TRUE(IsPrime(109));
  EXPECT_TRUE(IsPrime(113));

  // Test some composite numbers
  EXPECT_FALSE(IsPrime(99));  // 9 * 11
  EXPECT_FALSE(IsPrime(100)); // 10 * 10
  EXPECT_FALSE(IsPrime(102)); // 2 * 51
  EXPECT_FALSE(IsPrime(104)); // 8 * 13
  EXPECT_FALSE(IsPrime(105)); // 3 * 5 * 7
}

TEST_F(MathUtilsTest, IsPrimeLargerNumbers) {
  // Test some larger primes
  EXPECT_TRUE(IsPrime(997));
  EXPECT_TRUE(IsPrime(1009));
  EXPECT_TRUE(IsPrime(1013));
  EXPECT_TRUE(IsPrime(1019));

  // Test some larger composite numbers
  EXPECT_FALSE(IsPrime(1001)); // 7 * 11 * 13
  EXPECT_FALSE(IsPrime(1003)); // 17 * 59
  EXPECT_FALSE(IsPrime(1005)); // 3 * 5 * 67
}

TEST_F(MathUtilsTest, IsPrimePerfectSquares) {
  // Test perfect squares (should all be composite except for powers of primes)
  EXPECT_FALSE(IsPrime(4));   // 2^2
  EXPECT_FALSE(IsPrime(9));   // 3^2
  EXPECT_FALSE(IsPrime(16));  // 4^2
  EXPECT_FALSE(IsPrime(25));  // 5^2
  EXPECT_FALSE(IsPrime(36));  // 6^2
  EXPECT_FALSE(IsPrime(49));  // 7^2
  EXPECT_FALSE(IsPrime(64));  // 8^2
  EXPECT_FALSE(IsPrime(81));  // 9^2
  EXPECT_FALSE(IsPrime(100)); // 10^2
}

//===----------------------------------------------------------------------===//
// Prime Factorization Tests
//===----------------------------------------------------------------------===//

TEST_F(MathUtilsTest, TrialDivisionSmallNumbers) {
  auto factors = TrialDivision(1);
  EXPECT_TRUE(factors.empty());

  factors = TrialDivision(2);
  EXPECT_EQ(factors, std::vector<long long>({2}));

  factors = TrialDivision(12);
  std::vector<long long> expected = {2, 2, 3};
  EXPECT_EQ(factors, expected);

  factors = TrialDivision(17);
  EXPECT_EQ(factors, std::vector<long long>({17}));
}

TEST_F(MathUtilsTest, TrialDivisionMediumNumbers) {
  auto factors = TrialDivision(60);
  std::vector<long long> expected = {2, 2, 3, 5};
  EXPECT_EQ(factors, expected);

  factors = TrialDivision(100);
  expected = {2, 2, 5, 5};
  EXPECT_EQ(factors, expected);

  factors = TrialDivision(210);
  expected = {2, 3, 5, 7};
  EXPECT_EQ(factors, expected);
}

TEST_F(MathUtilsTest, PrimeFactorizeComprehensive) {
  // Test edge cases
  auto factors = PrimeFactorize(1);
  EXPECT_TRUE(factors.empty());

  // Test prime numbers
  factors = PrimeFactorize(17);
  EXPECT_EQ(factors, std::vector<long long>({17}));

  // Test composite numbers
  factors = PrimeFactorize(12);
  std::vector<long long> expected = {2, 2, 3};
  EXPECT_EQ(factors, expected);

  factors = PrimeFactorize(360);
  expected = {2, 2, 2, 3, 3, 5};
  EXPECT_EQ(factors, expected);
}

TEST_F(MathUtilsTest, FermatFactorizationBasicCases) {
  // Test edge cases
  auto factors = FermatFactorization(1);
  EXPECT_TRUE(factors.empty());

  factors = FermatFactorization(2);
  EXPECT_EQ(factors, std::vector<long long>({2}));

  // Test perfect squares
  factors = FermatFactorization(9);
  std::vector<long long> expected = {3, 3};
  EXPECT_EQ(factors, expected);

  factors = FermatFactorization(25);
  expected = {5, 5};
  EXPECT_EQ(factors, expected);
}

//===----------------------------------------------------------------------===//
// Prime Power Decomposition Tests
//===----------------------------------------------------------------------===//

TEST_F(MathUtilsTest, DecomposePrimePowerBasicCases) {
  // Test edge cases
  auto result = DecomposePrimePower(0);
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));

  result = DecomposePrimePower(1);
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));

  // Test prime powers
  result = DecomposePrimePower(2);
  EXPECT_EQ(result, std::make_pair(2ULL, 1ULL));

  result = DecomposePrimePower(4);
  EXPECT_EQ(result, std::make_pair(2ULL, 2ULL));

  result = DecomposePrimePower(8);
  EXPECT_EQ(result, std::make_pair(2ULL, 3ULL));

  result = DecomposePrimePower(16);
  EXPECT_EQ(result, std::make_pair(2ULL, 4ULL));
}

TEST_F(MathUtilsTest, DecomposePrimePowerOddPrimes) {
  auto result = DecomposePrimePower(3);
  EXPECT_EQ(result, std::make_pair(3ULL, 1ULL));

  result = DecomposePrimePower(9);
  EXPECT_EQ(result, std::make_pair(3ULL, 2ULL));

  result = DecomposePrimePower(27);
  EXPECT_EQ(result, std::make_pair(3ULL, 3ULL));

  result = DecomposePrimePower(5);
  EXPECT_EQ(result, std::make_pair(5ULL, 1ULL));

  result = DecomposePrimePower(25);
  EXPECT_EQ(result, std::make_pair(5ULL, 2ULL));

  result = DecomposePrimePower(125);
  EXPECT_EQ(result, std::make_pair(5ULL, 3ULL));
}

TEST_F(MathUtilsTest, DecomposePrimePowerNonPrimePowers) {
  // Test composite numbers that are not prime powers
  auto result = DecomposePrimePower(6);  // 2 * 3
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));

  result = DecomposePrimePower(10); // 2 * 5
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));

  result = DecomposePrimePower(12); // 2^2 * 3
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));

  result = DecomposePrimePower(15); // 3 * 5
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));

  result = DecomposePrimePower(30); // 2 * 3 * 5
  EXPECT_EQ(result, std::make_pair(0ULL, 0ULL));
}

//===----------------------------------------------------------------------===//
// Pollard's Rho Algorithm Tests
//===----------------------------------------------------------------------===//

TEST_F(MathUtilsTest, PollardsRhoBasicCases) {
  // Test small composite numbers
  long long factor = PollardsRho(15);
  EXPECT_TRUE(factor == 3 || factor == 5 || factor == 15);  // Allow for algorithm failure

  factor = PollardsRho(21);
  EXPECT_TRUE(factor == 3 || factor == 7 || factor == 21);  // Allow for algorithm failure

  // Test even numbers (should return 2)
  factor = PollardsRho(14);
  EXPECT_EQ(factor, 2);

  factor = PollardsRho(22);
  EXPECT_EQ(factor, 2);
}

TEST_F(MathUtilsTest, FactorizePollardsRhoComprehensive) {
  // Test factorization with recursive Pollard's Rho
  auto factors = FactorizePollardsRho(12);
  std::sort(factors.begin(), factors.end());
  std::vector<long long> expected = {2, 2, 3};
  EXPECT_EQ(factors, expected);

  factors = FactorizePollardsRho(15);
  std::sort(factors.begin(), factors.end());
  expected = {3, 5};
  EXPECT_EQ(factors, expected);

  // Test prime numbers (should return the number itself)
  factors = FactorizePollardsRho(17);
  EXPECT_EQ(factors, std::vector<long long>({17}));
}

//===----------------------------------------------------------------------===//
// Performance and Stress Tests
//===----------------------------------------------------------------------===//

TEST_F(MathUtilsTest, PrimeFactorizeConsistency) {
  // Test that different factorization methods give consistent results
  // for numbers where both trial division and Pollard's Rho are applicable
  std::vector<long long> test_numbers = {60, 84, 120, 180, 210, 300, 360, 420};

  for (long long n : test_numbers) {
    auto trial_factors = TrialDivision(n);
    auto smart_factors = PrimeFactorize(n);

    // Sort both results for comparison
    std::sort(trial_factors.begin(), trial_factors.end());
    std::sort(smart_factors.begin(), smart_factors.end());

    // Verify that the product of factors equals the original number
    long long trial_product = 1;
    for (long long factor : trial_factors) {
      trial_product *= factor;
    }

    long long smart_product = 1;
    for (long long factor : smart_factors) {
      smart_product *= factor;
    }

    EXPECT_EQ(trial_product, n) << "Trial division factors don't multiply to " << n;
    EXPECT_EQ(smart_product, n) << "Smart factorization factors don't multiply to " << n;

    // Both should produce the same product (the original number)
    EXPECT_EQ(trial_product, smart_product) << "Inconsistent factorization products for " << n;
  }
}

TEST_F(MathUtilsTest, IsPrimeConsistencyWithFactorization) {
  // Test that IsPrime is consistent with prime factorization
  std::vector<long long> test_numbers = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 97, 98, 99, 100, 101};

  for (long long n : test_numbers) {
    bool is_prime_result = IsPrime(n);
    auto factors = PrimeFactorize(n);
    bool is_prime_from_factors = (factors.size() == 1 && factors[0] == n);

    EXPECT_EQ(is_prime_result, is_prime_from_factors) << "Inconsistent prime test for " << n;
  }
}

//------------------------------------------------------------------------------
// SafeIntegerPower Tests
//------------------------------------------------------------------------------

TEST_F(MathUtilsTest, SafeIntegerPowerBasicCases) {
  // Test edge cases
  EXPECT_EQ(SafeIntegerPower(0, 0), 1); // 0^0 = 1 by convention
  EXPECT_EQ(SafeIntegerPower(0, 1), 0); // 0^n = 0 for n > 0
  EXPECT_EQ(SafeIntegerPower(0, 100), 0);

  EXPECT_EQ(SafeIntegerPower(1, 0), 1); // 1^0 = 1
  EXPECT_EQ(SafeIntegerPower(1, 1), 1); // 1^n = 1 for any n
  EXPECT_EQ(SafeIntegerPower(1, 1000), 1);

  EXPECT_EQ(SafeIntegerPower(5, 0), 1); // n^0 = 1 for any n
  EXPECT_EQ(SafeIntegerPower(100, 0), 1);
}

TEST_F(MathUtilsTest, SafeIntegerPowerSmallNumbers) {
  // Test basic powers
  EXPECT_EQ(SafeIntegerPower(2, 1), 2);
  EXPECT_EQ(SafeIntegerPower(2, 2), 4);
  EXPECT_EQ(SafeIntegerPower(2, 3), 8);
  EXPECT_EQ(SafeIntegerPower(2, 4), 16);
  EXPECT_EQ(SafeIntegerPower(2, 10), 1024);

  EXPECT_EQ(SafeIntegerPower(3, 1), 3);
  EXPECT_EQ(SafeIntegerPower(3, 2), 9);
  EXPECT_EQ(SafeIntegerPower(3, 3), 27);
  EXPECT_EQ(SafeIntegerPower(3, 4), 81);

  EXPECT_EQ(SafeIntegerPower(5, 2), 25);
  EXPECT_EQ(SafeIntegerPower(5, 3), 125);

  EXPECT_EQ(SafeIntegerPower(10, 3), 1000);
}

TEST_F(MathUtilsTest, SafeIntegerPowerLargerNumbers) {
  // Test some larger but safe powers
  EXPECT_EQ(SafeIntegerPower(2, 20), 1048576ULL);
  EXPECT_EQ(SafeIntegerPower(3, 10), 59049ULL);
  EXPECT_EQ(SafeIntegerPower(7, 8), 5764801ULL);

  // Test field order calculations similar to what we use in extensions
  EXPECT_EQ(SafeIntegerPower(2, 8), 256ULL);    // GF(2^8)
  EXPECT_EQ(SafeIntegerPower(2, 16), 65536ULL); // GF(2^16)
  EXPECT_EQ(SafeIntegerPower(3, 5), 243ULL);    // GF(3^5)
  EXPECT_EQ(SafeIntegerPower(5, 3), 125ULL);    // GF(5^3)
}

TEST_F(MathUtilsTest, SafeIntegerPowerOverflowDetection) {
  // Test overflow detection
  // 2^64 should overflow
  EXPECT_THROW(SafeIntegerPower(2, 64), std::overflow_error);

  // Large base with moderate exponent should overflow
  EXPECT_THROW(SafeIntegerPower(UINT64_MAX, 2), std::overflow_error);

  // Numbers that should be close to but not exceed the limit
  EXPECT_NO_THROW(SafeIntegerPower(2, 63)); // 2^63 should be fine
}

TEST_F(MathUtilsTest, SafeIntegerPowerComparisonWithStdPow) {
  // Compare with std::pow for small numbers where precision isn't an issue
  for (uint64_t base = 2; base <= 10; ++base) {
    for (uint64_t exp = 0; exp <= 10; ++exp) {
      uint64_t safe_result = SafeIntegerPower(base, exp);
      uint64_t std_result = static_cast<uint64_t>(std::pow(base, exp));
      EXPECT_EQ(safe_result, std_result)
          << "Mismatch for " << base << "^" << exp;
    }
  }
}
