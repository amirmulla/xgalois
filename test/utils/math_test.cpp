#include "xgalois/utils/math.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace xg::utils;

namespace {

uint64_t MultiplyFactors(const std::vector<uint64_t>& factors) {
  uint64_t prod = 1;
  for (uint64_t f : factors) {
    prod *= f;
  }
  return prod;
}

}  // namespace

class MathUtilsTest : public ::testing::Test {
 protected:
  void SetUp() override {}
};

TEST_F(MathUtilsTest, GcdBasicCases) {
  EXPECT_EQ(Gcd(12, 8), 4);
  EXPECT_EQ(Gcd(21, 14), 7);
  EXPECT_EQ(Gcd(17, 13), 1);
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
  EXPECT_EQ(Lcm(17, 13), 221);
  EXPECT_EQ(Lcm(12, 8), 24);
}

TEST_F(MathUtilsTest, ExtendedGcdBasicCases) {
  auto result = ExtendedGcd(240, 46);
  EXPECT_EQ(result.first, 2);

  int x = result.second.first;
  int y = result.second.second;
  EXPECT_EQ(240 * x + 46 * y, result.first);
}

TEST_F(MathUtilsTest, ExtendedGcdCoprimeNumbers) {
  auto result = ExtendedGcd(17, 13);
  EXPECT_EQ(result.first, 1);

  int x = result.second.first;
  int y = result.second.second;
  EXPECT_EQ(17 * x + 13 * y, 1);
}

TEST_F(MathUtilsTest, ExtendedGcdSignedNegativeInputs) {
  auto [g1, bez1] = ExtendedGcd(-240, 46);
  EXPECT_EQ(g1, 2);
  EXPECT_EQ(-240 * bez1.first + 46 * bez1.second, g1);

  auto [g2, bez2] = ExtendedGcd(240, -46);
  EXPECT_EQ(g2, 2);
  EXPECT_EQ(240 * bez2.first + (-46) * bez2.second, g2);

  auto [g3, bez3] = ExtendedGcd(-240, -46);
  EXPECT_EQ(g3, 2);
  EXPECT_EQ(-240 * bez3.first + (-46) * bez3.second, g3);

  auto [g4, bez4] = ExtendedGcd<int64_t>(-1000000007LL, -1000000009LL);
  EXPECT_EQ(g4, 1);
  __extension__ using i128 = __int128;
  i128 lhs = static_cast<i128>(-1000000007LL) * bez4.first +
             static_cast<i128>(-1000000009LL) * bez4.second;
  EXPECT_EQ(lhs, 1);
}

TEST_F(MathUtilsTest, ModMulAndModPow) {
  EXPECT_EQ(ModMul(0, 12345, 97), 0ULL);
  EXPECT_EQ(ModMul(12345, 67890, 1), 0ULL);
  EXPECT_THROW(ModMul(3, 4, 0), std::domain_error);

  // Near UINT64_MAX modular multiplication (no overflow).
  constexpr uint64_t kLargePrime = 18446744073709551557ULL;  // 2^64 - 59
  EXPECT_EQ(ModMul(kLargePrime - 1, kLargePrime - 1, kLargePrime), 1ULL);
  EXPECT_EQ(ModMul(kLargePrime - 2, kLargePrime - 3, kLargePrime), 6ULL);

  EXPECT_EQ(ModPow(2, 10, 1000), 24ULL);
  EXPECT_EQ(ModPow(3, 0, 7), 1ULL);
  EXPECT_EQ(ModPow(5, 10, 1), 0ULL);
  EXPECT_THROW(ModPow(2, 3, 0), std::domain_error);

  // Fermat's Little Theorem for largest 64-bit prime.
  EXPECT_EQ(ModPow(2, kLargePrime - 1, kLargePrime), 1ULL);
  EXPECT_EQ(ModPow(123456789ULL, kLargePrime - 1, kLargePrime), 1ULL);
}

TEST_F(MathUtilsTest, ModInvCases) {
  EXPECT_EQ(ModInv(3, 7), 5ULL);
  EXPECT_EQ(ModMul(3, ModInv(3, 7), 7), 1ULL);

  constexpr uint64_t kLargePrime = 18446744073709551557ULL;
  uint64_t inv = ModInv(1234567890123456789ULL, kLargePrime);
  EXPECT_EQ(ModMul(1234567890123456789ULL, inv, kLargePrime), 1ULL);

  EXPECT_THROW(ModInv(0, 7), std::domain_error);
  EXPECT_THROW(ModInv(2, 6), std::domain_error);
  EXPECT_THROW(ModInv(1, 1), std::domain_error);
  EXPECT_THROW(ModInv(1, 0), std::domain_error);
}

TEST_F(MathUtilsTest, ISqrtExactFullRange) {
  EXPECT_EQ(ISqrt(0), 0ULL);
  EXPECT_EQ(ISqrt(1), 1ULL);
  EXPECT_EQ(ISqrt(2), 1ULL);
  EXPECT_EQ(ISqrt(3), 1ULL);
  EXPECT_EQ(ISqrt(4), 2ULL);
  EXPECT_EQ(ISqrt(15), 3ULL);
  EXPECT_EQ(ISqrt(16), 4ULL);

  // Values where double (53-bit mantissa) loses precision.
  constexpr uint64_t kBigRoot = 4294967295ULL;  // 2^32 - 1
  constexpr uint64_t kBigSq = kBigRoot * kBigRoot;
  EXPECT_EQ(ISqrt(kBigSq - 1), kBigRoot - 1);
  EXPECT_EQ(ISqrt(kBigSq), kBigRoot);
  EXPECT_EQ(ISqrt(kBigSq + 1), kBigRoot);
  EXPECT_EQ(ISqrt(std::numeric_limits<uint64_t>::max()), kBigRoot);

  constexpr uint64_t kMidRoot = (1ULL << 31) + 12345ULL;
  constexpr uint64_t kMidSq = kMidRoot * kMidRoot;
  EXPECT_EQ(ISqrt(kMidSq - 1), kMidRoot - 1);
  EXPECT_EQ(ISqrt(kMidSq), kMidRoot);
  EXPECT_EQ(ISqrt(kMidSq + 1), kMidRoot);
}

TEST_F(MathUtilsTest, IsPrimeSmallNumbers) {
  EXPECT_FALSE(IsPrime(0));
  EXPECT_FALSE(IsPrime(1));

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
  EXPECT_TRUE(IsPrime(97));
  EXPECT_TRUE(IsPrime(101));
  EXPECT_TRUE(IsPrime(103));
  EXPECT_TRUE(IsPrime(107));
  EXPECT_TRUE(IsPrime(109));
  EXPECT_TRUE(IsPrime(113));

  EXPECT_FALSE(IsPrime(99));
  EXPECT_FALSE(IsPrime(100));
  EXPECT_FALSE(IsPrime(102));
  EXPECT_FALSE(IsPrime(104));
  EXPECT_FALSE(IsPrime(105));
}

TEST_F(MathUtilsTest, IsPrimeLargerNumbers) {
  EXPECT_TRUE(IsPrime(997));
  EXPECT_TRUE(IsPrime(1009));
  EXPECT_TRUE(IsPrime(1013));
  EXPECT_TRUE(IsPrime(1019));

  EXPECT_FALSE(IsPrime(1001));
  EXPECT_FALSE(IsPrime(1003));
  EXPECT_FALSE(IsPrime(1005));
}

TEST_F(MathUtilsTest, IsPrimePerfectSquares) {
  EXPECT_FALSE(IsPrime(4));
  EXPECT_FALSE(IsPrime(9));
  EXPECT_FALSE(IsPrime(16));
  EXPECT_FALSE(IsPrime(25));
  EXPECT_FALSE(IsPrime(36));
  EXPECT_FALSE(IsPrime(49));
  EXPECT_FALSE(IsPrime(64));
  EXPECT_FALSE(IsPrime(81));
  EXPECT_FALSE(IsPrime(100));
}

TEST_F(MathUtilsTest, IsPrimeCarmichaelAndLarge64BitCases) {
  // Carmichael numbers & famous strong pseudoprimes must all return false.
  const std::vector<uint64_t> carmichael_and_pseudoprimes = {
      561ULL,  1105ULL, 1729ULL, 2465ULL,
      2821ULL, 6601ULL, 8911ULL, 341550071728321ULL};
  for (uint64_t c : carmichael_and_pseudoprimes) {
    EXPECT_FALSE(IsPrime(c)) << "False positive on pseudoprime " << c;
  }

  // Known large primes up to largest 64-bit prime.
  EXPECT_TRUE(IsPrime(4294967291ULL));
  EXPECT_TRUE(IsPrime(2305843009213693951ULL));   // 2^61 - 1 (Mersenne prime)
  EXPECT_TRUE(IsPrime(18446744073709551557ULL));  // Largest uint64_t prime

  // Large semiprimes and prime squares > 10^7.
  EXPECT_FALSE(IsPrime(10000019ULL * 10000019ULL));
  EXPECT_FALSE(IsPrime(4294967291ULL * 4294967291ULL));
  EXPECT_FALSE(IsPrime(1000000007ULL * 1000000009ULL));
  EXPECT_FALSE(IsPrime(4294967291ULL * 4294967279ULL));
  EXPECT_FALSE(IsPrime(std::numeric_limits<uint64_t>::max()));
}

TEST_F(MathUtilsTest, TrialDivisionSmallNumbers) {
  auto factors = TrialDivision(1);
  EXPECT_TRUE(factors.empty());

  factors = TrialDivision(2);
  EXPECT_EQ(factors, std::vector<uint64_t>({2}));

  factors = TrialDivision(12);
  std::vector<uint64_t> expected = {2, 2, 3};
  EXPECT_EQ(factors, expected);

  factors = TrialDivision(17);
  EXPECT_EQ(factors, std::vector<uint64_t>({17}));
}

TEST_F(MathUtilsTest, TrialDivisionMediumNumbers) {
  auto factors = TrialDivision(60);
  std::vector<uint64_t> expected = {2, 2, 3, 5};
  EXPECT_EQ(factors, expected);

  factors = TrialDivision(100);
  expected = {2, 2, 5, 5};
  EXPECT_EQ(factors, expected);

  factors = TrialDivision(210);
  expected = {2, 3, 5, 7};
  EXPECT_EQ(factors, expected);
}

TEST_F(MathUtilsTest, PrimeFactorizeComprehensive) {
  EXPECT_TRUE(PrimeFactorize(0).empty());
  EXPECT_TRUE(PrimeFactorize(1).empty());
  EXPECT_EQ(PrimeFactorize(17), std::vector<uint64_t>({17}));
  EXPECT_EQ(PrimeFactorize(12), std::vector<uint64_t>({2, 2, 3}));
  EXPECT_EQ(PrimeFactorize(360), std::vector<uint64_t>({2, 2, 2, 3, 3, 5}));
  EXPECT_EQ(PrimeFactorize(121), std::vector<uint64_t>({11, 11}));
  EXPECT_EQ(PrimeFactorize(21218), std::vector<uint64_t>({2, 103, 103}));
  EXPECT_EQ(PrimeFactorize(643482), std::vector<uint64_t>({2, 3, 3, 7, 5107}));
  EXPECT_EQ(PrimeFactorize(309550), std::vector<uint64_t>({2, 5, 5, 41, 151}));
  EXPECT_EQ(PrimeFactorize(823543),
            std::vector<uint64_t>({7, 7, 7, 7, 7, 7, 7}));
}

TEST_F(MathUtilsTest, PrimeFactorizeLarge64BitNumbers) {
  const std::vector<uint64_t> large_targets = {
      10000019ULL * 10000019ULL,
      4294967291ULL * 4294967291ULL,
      1000000007ULL * 1000000009ULL,
      4294967291ULL * 4294967279ULL,
      341550071728321ULL,
      18446744073709551557ULL,
      std::numeric_limits<uint64_t>::max(),
  };

  for (uint64_t n : large_targets) {
    auto factors = PrimeFactorize(n);
    ASSERT_FALSE(factors.empty()) << "Empty factorization for " << n;
    EXPECT_TRUE(std::is_sorted(factors.begin(), factors.end()));
    EXPECT_EQ(MultiplyFactors(factors), n) << "Product mismatch for " << n;
    for (uint64_t p : factors) {
      EXPECT_TRUE(IsPrime(p)) << "Non-prime factor " << p << " for " << n;
    }
  }

  // Exact factorization of UINT64_MAX = 2^64 - 1.
  EXPECT_EQ(PrimeFactorize(std::numeric_limits<uint64_t>::max()),
            std::vector<uint64_t>({3, 5, 17, 257, 641, 65537, 6700417}));
}

TEST_F(MathUtilsTest, PrimeFactorizeHelpersAndEulerTotient) {
  auto with_mult = PrimeFactorizeWithMultiplicity(360);
  std::vector<std::pair<uint64_t, uint64_t>> expected_mult = {
      {2, 3}, {3, 2}, {5, 1}};
  EXPECT_EQ(with_mult, expected_mult);

  EXPECT_EQ(DistinctPrimeFactors(360), std::vector<uint64_t>({2, 3, 5}));
  EXPECT_TRUE(DistinctPrimeFactors(1).empty());

  EXPECT_EQ(EulerTotient(0), 0ULL);
  EXPECT_EQ(EulerTotient(1), 1ULL);
  EXPECT_EQ(EulerTotient(2), 1ULL);
  EXPECT_EQ(EulerTotient(9), 6ULL);
  EXPECT_EQ(EulerTotient(12), 4ULL);
  EXPECT_EQ(EulerTotient(36), 12ULL);
  EXPECT_EQ(EulerTotient(97), 96ULL);
  EXPECT_EQ(EulerTotient(1000000007ULL), 1000000006ULL);
  EXPECT_EQ(EulerTotient(18446744073709551557ULL), 18446744073709551556ULL);
}

TEST_F(MathUtilsTest, FermatFactorizationBasicAndCompositeCases) {
  EXPECT_TRUE(FermatFactorization(1).empty());
  EXPECT_EQ(FermatFactorization(2), std::vector<uint64_t>({2}));
  EXPECT_EQ(FermatFactorization(9), std::vector<uint64_t>({3, 3}));
  EXPECT_EQ(FermatFactorization(25), std::vector<uint64_t>({5, 5}));

  // Even numbers and composite squares (recursively factored).
  EXPECT_EQ(FermatFactorization(12), std::vector<uint64_t>({2, 2, 3}));
  EXPECT_EQ(FermatFactorization(100), std::vector<uint64_t>({2, 2, 5, 5}));
  EXPECT_EQ(FermatFactorization(81), std::vector<uint64_t>({3, 3, 3, 3}));
  EXPECT_EQ(FermatFactorization(225), std::vector<uint64_t>({3, 3, 5, 5}));

  // Fallback when factors are far apart.
  EXPECT_EQ(FermatFactorization(3ULL * 1000000007ULL),
            std::vector<uint64_t>({3ULL, 1000000007ULL}));
}

TEST_F(MathUtilsTest, DecomposePrimePowerBasicCases) {
  using U64Pair = std::pair<uint64_t, uint64_t>;
  EXPECT_EQ(DecomposePrimePower(0), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(1), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(2), (U64Pair{2, 1}));
  EXPECT_EQ(DecomposePrimePower(4), (U64Pair{2, 2}));
  EXPECT_EQ(DecomposePrimePower(8), (U64Pair{2, 3}));
  EXPECT_EQ(DecomposePrimePower(16), (U64Pair{2, 4}));
  EXPECT_EQ(DecomposePrimePower(1ULL << 63), (U64Pair{2, 63}));
}

TEST_F(MathUtilsTest, DecomposePrimePowerOddAndLargePrimes) {
  using U64Pair = std::pair<uint64_t, uint64_t>;
  EXPECT_EQ(DecomposePrimePower(3), (U64Pair{3, 1}));
  EXPECT_EQ(DecomposePrimePower(9), (U64Pair{3, 2}));
  EXPECT_EQ(DecomposePrimePower(27), (U64Pair{3, 3}));
  EXPECT_EQ(DecomposePrimePower(5), (U64Pair{5, 1}));
  EXPECT_EQ(DecomposePrimePower(25), (U64Pair{5, 2}));
  EXPECT_EQ(DecomposePrimePower(125), (U64Pair{5, 3}));
  EXPECT_EQ(DecomposePrimePower(121), (U64Pair{11, 2}));
  EXPECT_EQ(DecomposePrimePower(823543), (U64Pair{7, 7}));

  // Large prime squares and largest 64-bit prime.
  EXPECT_EQ(DecomposePrimePower(4294967291ULL * 4294967291ULL),
            (U64Pair{4294967291ULL, 2}));
  EXPECT_EQ(DecomposePrimePower(18446744073709551557ULL),
            (U64Pair{18446744073709551557ULL, 1}));
}

TEST_F(MathUtilsTest, DecomposePrimePowerNonPrimePowers) {
  using U64Pair = std::pair<uint64_t, uint64_t>;
  EXPECT_EQ(DecomposePrimePower(6), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(10), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(12), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(15), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(30), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(21218), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(643482), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(309550), (U64Pair{0, 0}));
  EXPECT_EQ(DecomposePrimePower(std::numeric_limits<uint64_t>::max()),
            (U64Pair{0, 0}));
}

TEST_F(MathUtilsTest, PollardsRhoBasicCases) {
  uint64_t factor = PollardsRho(15);
  EXPECT_TRUE(factor == 3 || factor == 5);

  factor = PollardsRho(21);
  EXPECT_TRUE(factor == 3 || factor == 7);

  EXPECT_EQ(PollardsRho(14), 2ULL);
  EXPECT_EQ(PollardsRho(22), 2ULL);

  uint64_t large_semi = 1000000007ULL * 1000000009ULL;
  uint64_t large_factor = PollardsRho(large_semi);
  EXPECT_TRUE(large_factor == 1000000007ULL || large_factor == 1000000009ULL);
}

TEST_F(MathUtilsTest, FactorizePollardsRhoComprehensive) {
  EXPECT_EQ(FactorizePollardsRho(12), std::vector<uint64_t>({2, 2, 3}));
  EXPECT_EQ(FactorizePollardsRho(15), std::vector<uint64_t>({3, 5}));
  EXPECT_EQ(FactorizePollardsRho(17), std::vector<uint64_t>({17}));
  EXPECT_EQ(FactorizePollardsRho(4294967291ULL * 4294967291ULL),
            std::vector<uint64_t>({4294967291ULL, 4294967291ULL}));
}

TEST_F(MathUtilsTest, PrimeFactorizeConsistency) {
  std::vector<uint64_t> test_numbers = {60, 84, 120, 180, 210, 300, 360, 420};

  for (uint64_t n : test_numbers) {
    auto trial_factors = TrialDivision(n);
    auto smart_factors = PrimeFactorize(n);

    std::sort(trial_factors.begin(), trial_factors.end());
    std::sort(smart_factors.begin(), smart_factors.end());

    EXPECT_EQ(MultiplyFactors(trial_factors), n);
    EXPECT_EQ(MultiplyFactors(smart_factors), n);
    EXPECT_EQ(trial_factors, smart_factors);
  }
}

TEST_F(MathUtilsTest, IsPrimeConsistencyWithFactorization) {
  std::vector<uint64_t> test_numbers = {2,  3,  4,  5,  6,  7,  8,   9,
                                        10, 11, 12, 13, 14, 15, 16,  17,
                                        18, 19, 20, 97, 98, 99, 100, 101};

  for (uint64_t n : test_numbers) {
    bool is_prime_result = IsPrime(n);
    auto factors = PrimeFactorize(n);
    bool is_prime_from_factors = (factors.size() == 1 && factors[0] == n);

    EXPECT_EQ(is_prime_result, is_prime_from_factors)
        << "Inconsistent prime test for " << n;
  }
}

TEST_F(MathUtilsTest, SafeIntegerPowerBasicCases) {
  EXPECT_EQ(SafeIntegerPower(0, 0), 1ULL);
  EXPECT_EQ(SafeIntegerPower(0, 1), 0ULL);
  EXPECT_EQ(SafeIntegerPower(0, 100), 0ULL);

  EXPECT_EQ(SafeIntegerPower(1, 0), 1ULL);
  EXPECT_EQ(SafeIntegerPower(1, 1), 1ULL);
  EXPECT_EQ(SafeIntegerPower(1, 1000), 1ULL);

  EXPECT_EQ(SafeIntegerPower(5, 0), 1ULL);
  EXPECT_EQ(SafeIntegerPower(100, 0), 1ULL);
}

TEST_F(MathUtilsTest, SafeIntegerPowerSmallNumbers) {
  EXPECT_EQ(SafeIntegerPower(2, 1), 2ULL);
  EXPECT_EQ(SafeIntegerPower(2, 2), 4ULL);
  EXPECT_EQ(SafeIntegerPower(2, 3), 8ULL);
  EXPECT_EQ(SafeIntegerPower(2, 4), 16ULL);
  EXPECT_EQ(SafeIntegerPower(2, 10), 1024ULL);

  EXPECT_EQ(SafeIntegerPower(3, 1), 3ULL);
  EXPECT_EQ(SafeIntegerPower(3, 2), 9ULL);
  EXPECT_EQ(SafeIntegerPower(3, 3), 27ULL);
  EXPECT_EQ(SafeIntegerPower(3, 4), 81ULL);

  EXPECT_EQ(SafeIntegerPower(5, 2), 25ULL);
  EXPECT_EQ(SafeIntegerPower(5, 3), 125ULL);

  EXPECT_EQ(SafeIntegerPower(10, 3), 1000ULL);
}

TEST_F(MathUtilsTest, SafeIntegerPowerLargerNumbers) {
  EXPECT_EQ(SafeIntegerPower(2, 20), 1048576ULL);
  EXPECT_EQ(SafeIntegerPower(3, 10), 59049ULL);
  EXPECT_EQ(SafeIntegerPower(7, 8), 5764801ULL);

  EXPECT_EQ(SafeIntegerPower(2, 8), 256ULL);
  EXPECT_EQ(SafeIntegerPower(2, 16), 65536ULL);
  EXPECT_EQ(SafeIntegerPower(3, 5), 243ULL);
  EXPECT_EQ(SafeIntegerPower(5, 3), 125ULL);
}

TEST_F(MathUtilsTest, SafeIntegerPowerOverflowDetection) {
  EXPECT_THROW(SafeIntegerPower(2, 64), std::overflow_error);
  EXPECT_THROW(SafeIntegerPower(std::numeric_limits<uint64_t>::max(), 2),
               std::overflow_error);
  EXPECT_NO_THROW(SafeIntegerPower(2, 63));
}

TEST_F(MathUtilsTest, SafeIntegerPowerComparisonWithStdPow) {
  for (uint64_t base = 2; base <= 10; ++base) {
    for (uint64_t exp = 0; exp <= 10; ++exp) {
      uint64_t safe_result = SafeIntegerPower(base, exp);
      uint64_t std_result = static_cast<uint64_t>(std::pow(base, exp));
      EXPECT_EQ(safe_result, std_result)
          << "Mismatch for " << base << "^" << exp;
    }
  }
}
