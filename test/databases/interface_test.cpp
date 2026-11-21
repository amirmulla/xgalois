#include "xgalois/databases/interface.hpp"

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "xgalois/utils/math.hpp"

namespace xg {
namespace databases {

struct PrimeTestCase {
  uint64_t number;
  std::vector<uint64_t> expected_factors;
  std::vector<int> expected_multiplicities;
  int expected_composite;
};

struct PolyTestCase {
  int characteristic;
  int degree;
  std::vector<int> expected_degrees;
  std::vector<int> expected_coeffs;
};

static const std::vector<PrimeTestCase> kPrimeTestCases = {
    {273323, {273323}, {1}, 0},
    {142903, {142903}, {1}, 0},
    {21218, {2, 103}, {1, 2}, 1},
    {643482, {2, 3, 7, 5107}, {1, 2, 1, 1}, 1},
    {309550, {2, 5, 41, 151}, {1, 2, 1, 1}, 1},
    {823543, {7}, {7}, 1},
    {121, {11}, {2}, 1},
};

static const std::vector<PolyTestCase> kIrreducibleTestCases = {
    {2, 2, {2, 1, 0}, {1, 1, 1}},
    {2, 3, {3, 1, 0}, {1, 1, 1}},
    {2, 4, {4, 1, 0}, {1, 1, 1}},
};

static const std::vector<PolyTestCase> kConwayTestCases = {
    {2, 1, {0, 1}, {1, 1}},
    {2, 2, {0, 1, 2}, {1, 1, 1}},
    {2, 3, {0, 1, 3}, {1, 1, 1}},
};

class DatabaseInterfaceTest : public ::testing::Test {
 protected:
  void SetUp() override {
    start_time_ = std::chrono::high_resolution_clock::now();
    VerifyPolynomialDatabaseFiles();
  }

  void TearDown() override {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time_);
    (void)duration;
  }

  void VerifyFactorizationInvariant(uint64_t n,
                                    const PrimeFactorsResult& result) {
    ASSERT_EQ(result.factors.size(), result.multiplicities.size());
    ASSERT_FALSE(result.factors.empty());

    uint64_t product = 1;
    for (size_t i = 0; i < result.factors.size(); ++i) {
      uint64_t p = result.factors[i];
      int e = result.multiplicities[i];
      EXPECT_TRUE(xg::utils::IsPrime(p)) << "Non-prime factor " << p;
      EXPECT_GT(e, 0);
      if (i > 0) {
        EXPECT_LT(result.factors[i - 1], p);
      }
      for (int k = 0; k < e; ++k) {
        product *= p;
      }
    }
    EXPECT_EQ(product, n) << "Factor product mismatch for " << n;
  }

  void TestPrimeFactors(const PrimeTestCase& test_case) {
    PrimeFactorsDatabase db;
    PrimeFactorsResult result = db.fetch(test_case.number);

    VerifyFactorizationInvariant(test_case.number, result);

    EXPECT_EQ(result.factors, test_case.expected_factors);
    EXPECT_EQ(result.multiplicities, test_case.expected_multiplicities);
    EXPECT_EQ(result.composite, test_case.expected_composite);
  }

  void TestPolynomialFetch(const PolyTestCase& test_case,
                           bool is_irreducible = true) {
    if (is_irreducible) {
      IrreduciblePolyDatabase db;
      IrreduciblePolyResult result =
          db.fetch(test_case.characteristic, test_case.degree);

      EXPECT_EQ(result.nonzero_degrees, test_case.expected_degrees);
      EXPECT_EQ(result.nonzero_coeffs, test_case.expected_coeffs);
    } else {
      ConwayPolyDatabase db;
      ConwayPolyResult result =
          db.fetch(test_case.characteristic, test_case.degree);

      EXPECT_EQ(result.nonzero_degrees, test_case.expected_degrees);
      EXPECT_EQ(result.nonzero_coeffs, test_case.expected_coeffs);
    }
  }

 private:
  std::chrono::high_resolution_clock::time_point start_time_;

  void VerifyPolynomialDatabaseFiles() {
    std::vector<std::string> db_files = {"irreducible_polys.db",
                                         "conway_polys.db"};
    for (const auto& file : db_files) {
      std::string resolved = detail::ResolveDbPath(file);
      if (!std::filesystem::exists(resolved)) {
        GTEST_SKIP() << "Database file not found: " << resolved;
      }
    }
  }
};

TEST_F(DatabaseInterfaceTest, PrimeFactorsFetchPrime) {
  TestPrimeFactors(kPrimeTestCases[0]);
}

TEST_F(DatabaseInterfaceTest, PrimeFactorsFetchComposite) {
  TestPrimeFactors(kPrimeTestCases[2]);
}

TEST_F(DatabaseInterfaceTest, PrimeFactorsFetchMultipleFactors) {
  TestPrimeFactors(kPrimeTestCases[3]);
}

TEST_F(DatabaseInterfaceTest,
       PrimeFactorsMultipleFactorsWithHigherMultiplicity) {
  TestPrimeFactors(kPrimeTestCases[4]);
}

TEST_F(DatabaseInterfaceTest, PrimeFactorsOnTheFlyLargeNumbers) {
  PrimeFactorsDatabase db;

  const std::vector<uint64_t> targets = {
      999999999ULL,
      4294967291ULL,
      4294967291ULL * 4294967291ULL,
      18446744073709551557ULL,
      std::numeric_limits<uint64_t>::max(),
  };

  for (uint64_t n : targets) {
    PrimeFactorsResult result = db.fetch(n);
    VerifyFactorizationInvariant(n, result);
    EXPECT_EQ(result.composite, xg::utils::IsPrime(n) ? 0 : 1);
  }
}

TEST_F(DatabaseInterfaceTest, PrimeFactorsFetchNotFound) {
  PrimeFactorsDatabase db;
  EXPECT_THROW(db.fetch(0), std::runtime_error);
  EXPECT_THROW(db.fetch(1), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, IrreduciblePolyFetchGF2Degree2) {
  TestPolynomialFetch(kIrreducibleTestCases[0]);
}

TEST_F(DatabaseInterfaceTest, IrreduciblePolyFetchGF2Degree3) {
  TestPolynomialFetch(kIrreducibleTestCases[1]);
}

TEST_F(DatabaseInterfaceTest, IrreduciblePolyFetchGF2Degree4) {
  TestPolynomialFetch(kIrreducibleTestCases[2]);
}

TEST_F(DatabaseInterfaceTest, IrreduciblePolyFetchAll) {
  IrreduciblePolyDatabase db;

  auto deg2 = db.fetch_all(2, 2);
  ASSERT_EQ(deg2.size(), 1U);
  EXPECT_EQ(deg2[0].nonzero_degrees, std::vector<int>({2, 1, 0}));

  // irreducible_polys.db has 3 primitive polynomials for GF(2^5).
  auto deg5 = db.fetch_all(2, 5);
  EXPECT_EQ(deg5.size(), 3U);
  for (const auto& poly : deg5) {
    EXPECT_EQ(poly.nonzero_degrees.size(), poly.nonzero_coeffs.size());
    EXPECT_FALSE(poly.nonzero_degrees.empty());
    EXPECT_EQ(poly.nonzero_degrees.front(), 5);
  }

  auto missing = db.fetch_all(999, 999);
  EXPECT_TRUE(missing.empty());
}

TEST_F(DatabaseInterfaceTest, IrreduciblePolyFetchNotFound) {
  IrreduciblePolyDatabase db;
  EXPECT_THROW(db.fetch(999, 999), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, ConwayPolyFetchGF2Degree1) {
  TestPolynomialFetch(kConwayTestCases[0], false);
}

TEST_F(DatabaseInterfaceTest, ConwayPolyFetchGF2Degree2) {
  TestPolynomialFetch(kConwayTestCases[1], false);
}

TEST_F(DatabaseInterfaceTest, ConwayPolyFetchGF2Degree3) {
  TestPolynomialFetch(kConwayTestCases[2], false);
}

TEST_F(DatabaseInterfaceTest, EmbeddedConwayPolysZeroDiskIO) {
  ConwayPolyDatabase db;

  // All GF(2^m) for m = 1..64 must be served from the embedded constexpr table.
  for (int m = 1; m <= 64; ++m) {
    auto embedded =
        detail::LookupEmbeddedConwayPoly(2U, static_cast<uint32_t>(m));
    ASSERT_TRUE(embedded.has_value()) << "Missing embedded GF(2^" << m << ")";
    EXPECT_FALSE(embedded->empty());
    EXPECT_EQ(embedded->back().degree, static_cast<uint16_t>(m));
    EXPECT_EQ(embedded->back().coeff, 1U);

    ConwayPolyResult result = db.fetch(2, m);
    ASSERT_EQ(result.nonzero_degrees.size(), embedded->size());
    EXPECT_EQ(result.nonzero_degrees.back(), m);
  }

  // Common odd primes embedded entries.
  const std::vector<std::pair<int, int>> odd_embedded = {
      {3, 2},  {3, 4},  {5, 2},  {5, 3},  {7, 2},  {11, 2},
      {13, 2}, {17, 2}, {19, 2}, {23, 2}, {29, 2}, {31, 2}};
  for (const auto& [p, d] : odd_embedded) {
    EXPECT_TRUE(detail::LookupEmbeddedConwayPoly(static_cast<uint32_t>(p),
                                                 static_cast<uint32_t>(d))
                    .has_value())
        << "Missing embedded GF(" << p << "^" << d << ")";
    ConwayPolyResult result = db.fetch(p, d);
    EXPECT_EQ(result.nonzero_degrees.back(), d);
    EXPECT_EQ(result.nonzero_coeffs.back(), 1);
  }

  // Extended entry loaded from conway_polys.db (not embedded).
  EXPECT_FALSE(detail::LookupEmbeddedConwayPoly(2U, 65U).has_value());
  ConwayPolyResult ext = db.fetch(2, 65);
  EXPECT_EQ(ext.nonzero_degrees.back(), 65);
}

TEST_F(DatabaseInterfaceTest, ConwayPolyFetchNotFound) {
  ConwayPolyDatabase db;
  EXPECT_THROW(db.fetch(999, 999), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, MissingFileErrorHandling) {
  IrreduciblePolyDatabase missing_irred("/nonexistent/path/irreducible.db");
  ConwayPolyDatabase missing_conway("/nonexistent/path/conway.db");

  EXPECT_THROW(missing_irred.fetch(2, 2), std::runtime_error);
  EXPECT_THROW(missing_irred.fetch_all(2, 2), std::runtime_error);
  EXPECT_THROW(missing_conway.fetch(2, 2), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, DatabasePerformance) {
  PrimeFactorsDatabase db;

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<uint64_t> test_numbers = {273323, 142903, 21218, 643482, 309550};
  for (uint64_t num : test_numbers) {
    PrimeFactorsResult result = db.fetch(num);
    EXPECT_GT(result.factors.size(), 0U);
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  EXPECT_LT(duration.count(), 1000);
}

TEST_F(DatabaseInterfaceTest, DatabaseDataConsistency) {
  PrimeFactorsDatabase prime_db;
  IrreduciblePolyDatabase irreducible_db;
  ConwayPolyDatabase conway_db;

  PrimeFactorsResult prime_result = prime_db.fetch(273323);
  EXPECT_EQ(prime_result.factors.size(), prime_result.multiplicities.size());
  EXPECT_FALSE(prime_result.factors.empty());

  IrreduciblePolyResult irreducible_result = irreducible_db.fetch(2, 2);
  EXPECT_EQ(irreducible_result.nonzero_degrees.size(),
            irreducible_result.nonzero_coeffs.size());
  EXPECT_FALSE(irreducible_result.nonzero_degrees.empty());

  ConwayPolyResult conway_result = conway_db.fetch(2, 2);
  EXPECT_EQ(conway_result.nonzero_degrees.size(),
            conway_result.nonzero_coeffs.size());
  EXPECT_FALSE(conway_result.nonzero_degrees.empty());
}

TEST_F(DatabaseInterfaceTest, MultipleDatabaseOperations) {
  PrimeFactorsDatabase prime_db;
  IrreduciblePolyDatabase irreducible_db;
  ConwayPolyDatabase conway_db;

  PrimeFactorsResult prime1 = prime_db.fetch(273323);
  PrimeFactorsResult prime2 = prime_db.fetch(21218);
  EXPECT_NE(prime1.composite, prime2.composite);

  IrreduciblePolyResult irreducible1 = irreducible_db.fetch(2, 2);
  IrreduciblePolyResult irreducible2 = irreducible_db.fetch(2, 3);
  EXPECT_NE(irreducible1.nonzero_degrees[0], irreducible2.nonzero_degrees[0]);

  ConwayPolyResult conway1 = conway_db.fetch(2, 1);
  ConwayPolyResult conway2 = conway_db.fetch(2, 2);
  EXPECT_NE(conway1.nonzero_degrees.size(), conway2.nonzero_degrees.size());
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringBinary) {
  std::vector<int> degrees = {2, 1, 0};
  std::vector<int> coeffs = {1, 1, 1};
  EXPECT_EQ(ConvertPolyResultToString(degrees, coeffs, 2), "x^2+x+1");

  // Even coefficients must reduce to 0 mod 2 (not 1).
  EXPECT_EQ(ConvertPolyResultToString({2, 1, 0}, {2, 3, 4}, 2), "x");
  EXPECT_EQ(ConvertPolyResultToString({2, 1, 0}, {2, 4, 6}, 2), "0");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringPrime) {
  std::vector<int> degrees = {3, 1, 0};
  std::vector<int> coeffs = {2, 3, 1};
  EXPECT_EQ(ConvertPolyResultToString(degrees, coeffs, 5), "2*x^3+3*x+1");

  // Coefficient reduction modulo characteristic (including negative values).
  EXPECT_EQ(ConvertPolyResultToString({2, 1, 0}, {5, 6, -1}, 5), "x+4");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringSingleTerm) {
  EXPECT_EQ(ConvertPolyResultToString({5}, {1}, 7), "x^5");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringConstant) {
  EXPECT_EQ(ConvertPolyResultToString({0}, {4}, 5), "4");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringEmpty) {
  EXPECT_EQ(ConvertPolyResultToString({}, {}, 2), "0");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringMismatchAndInvalidChar) {
  EXPECT_THROW(ConvertPolyResultToString({1}, {1, 1}, 2),
               std::invalid_argument);
  EXPECT_THROW(ConvertPolyResultToString({1}, {1}, 1), std::invalid_argument);
  EXPECT_THROW(ConvertPolyResultToString({1}, {1}, 0), std::invalid_argument);
}

TEST_F(DatabaseInterfaceTest, GetIrreduciblePolynomialGF2Degree2) {
  EXPECT_EQ(GetIrreduciblePolynomial(2, 2), "x^2+x+1");
}

TEST_F(DatabaseInterfaceTest, GetIrreduciblePolynomialGF2Degree3) {
  EXPECT_EQ(GetIrreduciblePolynomial(2, 3), "x^3+x+1");
}

TEST_F(DatabaseInterfaceTest, GetIrreduciblePolynomialNotFound) {
  EXPECT_THROW(GetIrreduciblePolynomial(999, 999), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, GetConwayPolynomialGF2Degree2) {
  EXPECT_EQ(GetConwayPolynomial(2, 2), "x^2+x+1");
}

TEST_F(DatabaseInterfaceTest, GetConwayPolynomialGF3Degree2) {
  EXPECT_EQ(GetConwayPolynomial(3, 2), "x^2+2*x+2");
}

TEST_F(DatabaseInterfaceTest, GetConwayPolynomialNotFound) {
  EXPECT_THROW(GetConwayPolynomial(999, 999), std::runtime_error);
}

class PrimeFactorsParameterizedTest
    : public DatabaseInterfaceTest,
      public ::testing::WithParamInterface<PrimeTestCase> {};

TEST_P(PrimeFactorsParameterizedTest, FactorizationTest) {
  TestPrimeFactors(GetParam());
}

INSTANTIATE_TEST_SUITE_P(AllPrimeFactorizations, PrimeFactorsParameterizedTest,
                         ::testing::ValuesIn(kPrimeTestCases));

class IrreduciblePolyParameterizedTest
    : public DatabaseInterfaceTest,
      public ::testing::WithParamInterface<PolyTestCase> {};

TEST_P(IrreduciblePolyParameterizedTest, IrreduciblePolyTest) {
  TestPolynomialFetch(GetParam(), true);
}

INSTANTIATE_TEST_SUITE_P(AllIrreduciblePolynomials,
                         IrreduciblePolyParameterizedTest,
                         ::testing::ValuesIn(kIrreducibleTestCases));

class ConwayPolyParameterizedTest
    : public DatabaseInterfaceTest,
      public ::testing::WithParamInterface<PolyTestCase> {};

TEST_P(ConwayPolyParameterizedTest, ConwayPolyTest) {
  TestPolynomialFetch(GetParam(), false);
}

INSTANTIATE_TEST_SUITE_P(AllConwayPolynomials, ConwayPolyParameterizedTest,
                         ::testing::ValuesIn(kConwayTestCases));

TEST_F(DatabaseInterfaceTest, PrimeFactorsEdgeCases) {
  PrimeFactorsDatabase db;

  for (uint64_t invalid_num : {0ULL, 1ULL}) {
    EXPECT_THROW(db.fetch(invalid_num), std::runtime_error)
        << "Expected exception for number: " << invalid_num;
  }
}

TEST_F(DatabaseInterfaceTest, PolynomialEdgeCases) {
  IrreduciblePolyDatabase irreducible_db;
  ConwayPolyDatabase conway_db;

  std::vector<std::pair<int, int>> invalid_params = {
      {0, 1}, {1, 1}, {-1, 2}, {2, 0}, {2, -1}};

  for (const auto& [char_val, deg_val] : invalid_params) {
    EXPECT_THROW(irreducible_db.fetch(char_val, deg_val), std::runtime_error)
        << "Expected exception for characteristic=" << char_val
        << ", degree=" << deg_val;
    EXPECT_THROW(conway_db.fetch(char_val, deg_val), std::runtime_error)
        << "Expected exception for characteristic=" << char_val
        << ", degree=" << deg_val;
  }
}

TEST_F(DatabaseInterfaceTest, StressTestMultipleLookups) {
  PrimeFactorsDatabase prime_db;
  IrreduciblePolyDatabase irreducible_db;
  ConwayPolyDatabase conway_db;

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < 50; ++i) {
    for (const auto& test_case : kPrimeTestCases) {
      PrimeFactorsResult result = prime_db.fetch(test_case.number);
      EXPECT_FALSE(result.factors.empty());
    }

    for (const auto& test_case : kIrreducibleTestCases) {
      IrreduciblePolyResult result =
          irreducible_db.fetch(test_case.characteristic, test_case.degree);
      EXPECT_FALSE(result.nonzero_degrees.empty());
    }

    for (const auto& test_case : kConwayTestCases) {
      ConwayPolyResult result =
          conway_db.fetch(test_case.characteristic, test_case.degree);
      EXPECT_FALSE(result.nonzero_degrees.empty());
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  EXPECT_LT(duration.count(), 2000)
      << "Stress test took too long: " << duration.count() << "ms";
}

TEST_F(DatabaseInterfaceTest, ConcurrentDatabaseAccess) {
  constexpr int kNumThreads = 8;
  constexpr int kIterations = 50;
  std::atomic<int> failures{0};

  std::vector<std::thread> threads;
  threads.reserve(kNumThreads);

  for (int t = 0; t < kNumThreads; ++t) {
    threads.emplace_back([&failures, t]() {
      PrimeFactorsDatabase prime_db;
      IrreduciblePolyDatabase irred_db;
      ConwayPolyDatabase conway_db;

      for (int i = 0; i < kIterations; ++i) {
        try {
          auto p_res = prime_db.fetch(21218);
          if (p_res.factors != std::vector<uint64_t>({2, 103}) ||
              p_res.multiplicities != std::vector<int>({1, 2}) ||
              p_res.composite != 1) {
            ++failures;
          }

          auto i_res = irred_db.fetch(2, 2 + (t + i) % 3);
          if (i_res.nonzero_degrees.empty() ||
              i_res.nonzero_degrees.size() != i_res.nonzero_coeffs.size()) {
            ++failures;
          }

          auto c_emb = conway_db.fetch(2, 1 + (t + i) % 64);
          if (c_emb.nonzero_degrees.empty() ||
              c_emb.nonzero_degrees.back() != 1 + (t + i) % 64) {
            ++failures;
          }

          // Extended file-backed entry accessed concurrently across threads.
          auto c_ext = conway_db.fetch(2, 65);
          if (c_ext.nonzero_degrees.empty() ||
              c_ext.nonzero_degrees.back() != 65) {
            ++failures;
          }
        } catch (...) {
          ++failures;
        }
      }
    });
  }

  for (auto& th : threads) {
    th.join();
  }

  EXPECT_EQ(failures.load(), 0);
}

}  // namespace databases
}  // namespace xg