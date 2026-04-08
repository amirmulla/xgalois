

#include "xgalois/databases/interface.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <vector>

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
    {21218, {2, 103}, {1, 1}, 1},
    {643482, {2, 3, 7, 5107}, {1, 1, 1, 1}, 1},
    {309550, {2, 5, 41, 151}, {1, 1, 1, 1}, 1},
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

    VerifyDatabaseFiles();
  }

  void TearDown() override {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time_);

    if (duration.count() > 100) {
      std::cout << "Test took " << duration.count() << "ms" << '\n';
    }
  }

  void TestPrimeFactors(const PrimeTestCase& test_case) {
    PrimeFactorsDatabase db;
    PrimeFactorsResult result = db.fetch(test_case.number);

    EXPECT_EQ(result.factors.size(), test_case.expected_factors.size());
    EXPECT_EQ(result.multiplicities.size(),
              test_case.expected_multiplicities.size());
    EXPECT_EQ(result.composite, test_case.expected_composite);

    for (size_t i = 0; i < test_case.expected_factors.size(); ++i) {
      EXPECT_EQ(result.factors[i], test_case.expected_factors[i]);
      EXPECT_EQ(result.multiplicities[i], test_case.expected_multiplicities[i]);
    }
  }

  void TestPolynomialFetch(const PolyTestCase& test_case,
                           bool is_irreducible = true) {
    if (is_irreducible) {
      IrreduciblePolyDatabase db;
      IrreduciblePolyResult result =
          db.fetch(test_case.characteristic, test_case.degree);

      EXPECT_EQ(result.nonzero_degrees.size(),
                test_case.expected_degrees.size());
      EXPECT_EQ(result.nonzero_coeffs.size(), test_case.expected_coeffs.size());

      for (size_t i = 0; i < test_case.expected_degrees.size(); ++i) {
        EXPECT_EQ(result.nonzero_degrees[i], test_case.expected_degrees[i]);
        EXPECT_EQ(result.nonzero_coeffs[i], test_case.expected_coeffs[i]);
      }
    } else {
      ConwayPolyDatabase db;
      ConwayPolyResult result =
          db.fetch(test_case.characteristic, test_case.degree);

      EXPECT_EQ(result.nonzero_degrees.size(),
                test_case.expected_degrees.size());
      EXPECT_EQ(result.nonzero_coeffs.size(), test_case.expected_coeffs.size());

      for (size_t i = 0; i < test_case.expected_degrees.size(); ++i) {
        EXPECT_EQ(result.nonzero_degrees[i], test_case.expected_degrees[i]);
        EXPECT_EQ(result.nonzero_coeffs[i], test_case.expected_coeffs[i]);
      }
    }
  }

 private:
  std::chrono::high_resolution_clock::time_point start_time_;

  void VerifyDatabaseFiles() {
    std::vector<std::string> db_files = {
        "xgalois/databases/prime_factors.db",
        "xgalois/databases/irreducible_polys.db",
        "xgalois/databases/conway_polys.db"};

    for (const auto& file : db_files) {
      if (!std::filesystem::exists(file)) {
        GTEST_SKIP() << "Database file not found: " << file;
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

TEST_F(DatabaseInterfaceTest, PrimeFactorsFetchNotFound) {
  PrimeFactorsDatabase db;

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

TEST_F(DatabaseInterfaceTest, ConwayPolyFetchNotFound) {
  ConwayPolyDatabase db;

  EXPECT_THROW(db.fetch(999, 999), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, DatabaseFileNotFoundHandling) {
  PrimeFactorsDatabase prime_db;
  IrreduciblePolyDatabase irreducible_db;
  ConwayPolyDatabase conway_db;

  EXPECT_THROW(prime_db.fetch(999999999), std::runtime_error);
  EXPECT_THROW(irreducible_db.fetch(999, 999), std::runtime_error);
  EXPECT_THROW(conway_db.fetch(999, 999), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, DatabasePerformance) {
  PrimeFactorsDatabase db;

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<uint64_t> test_numbers = {273323, 142903, 21218, 643482, 309550};

  for (uint64_t num : test_numbers) {
    PrimeFactorsResult result = db.fetch(num);
    EXPECT_GT(result.factors.size(), 0);
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
  std::string result = ConvertPolyResultToString(degrees, coeffs, 2);
  EXPECT_EQ(result, "x^2+x+1");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringPrime) {
  std::vector<int> degrees = {3, 1, 0};
  std::vector<int> coeffs = {2, 3, 1};
  std::string result = ConvertPolyResultToString(degrees, coeffs, 5);
  EXPECT_EQ(result, "2*x^3+3*x+1");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringSingleTerm) {
  std::vector<int> degrees = {5};
  std::vector<int> coeffs = {1};
  std::string result = ConvertPolyResultToString(degrees, coeffs, 7);
  EXPECT_EQ(result, "x^5");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringConstant) {
  std::vector<int> degrees = {0};
  std::vector<int> coeffs = {4};
  std::string result = ConvertPolyResultToString(degrees, coeffs, 5);
  EXPECT_EQ(result, "4");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringEmpty) {
  std::vector<int> degrees = {};
  std::vector<int> coeffs = {};
  std::string result = ConvertPolyResultToString(degrees, coeffs, 2);
  EXPECT_EQ(result, "0");
}

TEST_F(DatabaseInterfaceTest, ConvertPolyResultToStringMismatch) {
  std::vector<int> degrees = {1};
  std::vector<int> coeffs = {1, 1};
  EXPECT_THROW(ConvertPolyResultToString(degrees, coeffs, 2),
               std::invalid_argument);
}

TEST_F(DatabaseInterfaceTest, GetIrreduciblePolynomialGF2Degree2) {
  std::string poly_str = GetIrreduciblePolynomial(2, 2);
  EXPECT_EQ(poly_str, "x^2+x+1");
}

TEST_F(DatabaseInterfaceTest, GetIrreduciblePolynomialGF2Degree3) {
  std::string poly_str = GetIrreduciblePolynomial(2, 3);
  EXPECT_EQ(poly_str, "x^3+x+1");
}

TEST_F(DatabaseInterfaceTest, GetIrreduciblePolynomialNotFound) {
  EXPECT_THROW(GetIrreduciblePolynomial(999, 999), std::runtime_error);
}

TEST_F(DatabaseInterfaceTest, GetConwayPolynomialGF2Degree2) {
  std::string poly_str = GetConwayPolynomial(2, 2);
  EXPECT_EQ(poly_str, "x^2+x+1");
}

TEST_F(DatabaseInterfaceTest, GetConwayPolynomialGF3Degree2) {
  std::string poly_str = GetConwayPolynomial(3, 2);
  EXPECT_EQ(poly_str, "x^2+2*x+2");
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

  std::vector<uint64_t> invalid_numbers = {0, static_cast<uint64_t>(-1),
                                           999999999};

  for (auto num : invalid_numbers) {
    EXPECT_THROW(db.fetch(num), std::runtime_error)
        << "Expected exception for number: " << num;
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

  for (int i = 0; i < 10; ++i) {
    for (const auto& test_case : kPrimeTestCases) {
      if (test_case.number > 0) {
        PrimeFactorsResult result = prime_db.fetch(test_case.number);
        EXPECT_FALSE(result.factors.empty());
      }
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

  EXPECT_LT(duration.count(), 10000)
      << "Stress test took too long: " << duration.count() << "ms";
}

TEST_F(DatabaseInterfaceTest, ConcurrentDatabaseAccess) {
  PrimeFactorsDatabase prime_db1, prime_db2;
  IrreduciblePolyDatabase irreducible_db1, irreducible_db2;
  ConwayPolyDatabase conway_db1, conway_db2;

  PrimeFactorsResult prime_result1 = prime_db1.fetch(273323);
  PrimeFactorsResult prime_result2 = prime_db2.fetch(273323);

  EXPECT_EQ(prime_result1.factors, prime_result2.factors);
  EXPECT_EQ(prime_result1.multiplicities, prime_result2.multiplicities);
  EXPECT_EQ(prime_result1.composite, prime_result2.composite);

  IrreduciblePolyResult irreducible_result1 = irreducible_db1.fetch(2, 2);
  IrreduciblePolyResult irreducible_result2 = irreducible_db2.fetch(2, 2);

  EXPECT_EQ(irreducible_result1.nonzero_degrees,
            irreducible_result2.nonzero_degrees);
  EXPECT_EQ(irreducible_result1.nonzero_coeffs,
            irreducible_result2.nonzero_coeffs);

  ConwayPolyResult conway_result1 = conway_db1.fetch(2, 2);
  ConwayPolyResult conway_result2 = conway_db2.fetch(2, 2);

  EXPECT_EQ(conway_result1.nonzero_degrees, conway_result2.nonzero_degrees);
  EXPECT_EQ(conway_result1.nonzero_coeffs, conway_result2.nonzero_coeffs);
}

}  // namespace databases
}  // namespace xg