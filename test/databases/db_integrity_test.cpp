#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/databases/embedded_conway_polys.hpp"
#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace {

PolynomialDense<GF2> BuildBinaryPoly(
    const std::shared_ptr<GF2> &gf2, uint32_t expected_degree,
    std::span<const databases::detail::EmbeddedConwayTerm> terms) {
  std::vector<uint8_t> coeffs(expected_degree + 1, uint8_t{0});
  for (const auto &term : terms) {
    EXPECT_LE(term.degree, expected_degree);
    EXPECT_EQ(term.coeff, uint16_t{1});
    if (term.degree <= expected_degree) {
      coeffs[term.degree] = static_cast<uint8_t>(term.coeff & 1U);
    }
  }
  return PolynomialDense<GF2>::FromRaw(std::move(coeffs), gf2);
}

PolynomialDense<GFP<uint16_t>> BuildPrimePoly(
    const std::shared_ptr<GFP<uint16_t>> &gfp, uint32_t expected_degree,
    std::span<const databases::detail::EmbeddedConwayTerm> terms) {
  const uint16_t p = static_cast<uint16_t>(gfp->Characteristic());
  std::vector<uint16_t> coeffs(expected_degree + 1, uint16_t{0});
  for (const auto &term : terms) {
    EXPECT_LE(term.degree, expected_degree);
    EXPECT_LT(term.coeff, p);
    if (term.degree <= expected_degree) {
      coeffs[term.degree] = static_cast<uint16_t>(term.coeff % p);
    }
  }
  return PolynomialDense<GFP<uint16_t>>::FromRaw(std::move(coeffs), gfp);
}

PolynomialDense<GF2> BuildBinaryPolyFromResult(
    const std::shared_ptr<GF2> &gf2, int expected_degree,
    const databases::IrreduciblePolyResult &res) {
  EXPECT_EQ(res.nonzero_degrees.size(), res.nonzero_coeffs.size());
  std::vector<uint8_t> coeffs(static_cast<size_t>(expected_degree) + 1,
                              uint8_t{0});
  for (size_t i = 0; i < res.nonzero_degrees.size(); ++i) {
    const int d = res.nonzero_degrees[i];
    EXPECT_GE(d, 0);
    EXPECT_LE(d, expected_degree);
    if (d >= 0 && d <= expected_degree) {
      coeffs[static_cast<size_t>(d)] =
          static_cast<uint8_t>(res.nonzero_coeffs[i] & 1);
    }
  }
  return PolynomialDense<GF2>::FromRaw(std::move(coeffs), gf2);
}

PolynomialDense<GFP<uint16_t>> BuildPrimePolyFromResult(
    const std::shared_ptr<GFP<uint16_t>> &gfp, int expected_degree,
    const databases::IrreduciblePolyResult &res) {
  EXPECT_EQ(res.nonzero_degrees.size(), res.nonzero_coeffs.size());
  const int p = static_cast<int>(gfp->Characteristic());
  std::vector<uint16_t> coeffs(static_cast<size_t>(expected_degree) + 1,
                               uint16_t{0});
  for (size_t i = 0; i < res.nonzero_degrees.size(); ++i) {
    const int d = res.nonzero_degrees[i];
    EXPECT_GE(d, 0);
    EXPECT_LE(d, expected_degree);
    if (d >= 0 && d <= expected_degree) {
      int c = res.nonzero_coeffs[i] % p;
      if (c < 0) {
        c += p;
      }
      coeffs[static_cast<size_t>(d)] = static_cast<uint16_t>(c);
    }
  }
  return PolynomialDense<GFP<uint16_t>>::FromRaw(std::move(coeffs), gfp);
}

void VerifyPrimeFactorization(const databases::PrimeFactorsDatabase &db,
                              uint64_t n) {
  const auto res = db.Fetch(n);
  ASSERT_EQ(res.factors.size(), res.multiplicities.size()) << "n=" << n;
  ASSERT_FALSE(res.factors.empty()) << "n=" << n;

  uint64_t product = 1;
  for (size_t i = 0; i < res.factors.size(); ++i) {
    const uint64_t p = res.factors[i];
    const int e = res.multiplicities[i];
    EXPECT_TRUE(utils::IsPrime(p))
        << "Non-prime factor " << p << " for n=" << n;
    EXPECT_GT(e, 0) << "Non-positive multiplicity for n=" << n;
    if (i > 0) {
      EXPECT_LT(res.factors[i - 1], p)
          << "Factors not strictly sorted for n=" << n;
    }
    for (int k = 0; k < e; ++k) {
      product *= p;
    }
  }
  EXPECT_EQ(product, n) << "Factor product mismatch for n=" << n;

  const int expected_composite =
      (res.factors.size() == 1 && res.multiplicities[0] == 1) ? 0 : 1;
  EXPECT_EQ(res.composite, expected_composite)
      << "Composite flag mismatch for n=" << n;
  EXPECT_EQ(res.composite == 0, utils::IsPrime(n))
      << "Primality mismatch for n=" << n;
}

// ============================================================================
// Embedded Conway Polynomial Table Integrity
// ============================================================================

TEST(DbIntegrityTest, EmbeddedConwayTableStructuralInvariants) {
  using namespace databases::detail;
  ASSERT_FALSE(kEmbeddedConwayEntries.empty());

  for (size_t i = 0; i < kEmbeddedConwayEntries.size(); ++i) {
    const auto &entry = kEmbeddedConwayEntries[i];
    EXPECT_TRUE(utils::IsPrime(entry.characteristic))
        << "Non-prime characteristic " << entry.characteristic;
    EXPECT_GE(entry.degree, 1U);
    EXPECT_GE(entry.count, 2U);
    EXPECT_LE(static_cast<size_t>(entry.offset) + entry.count,
              kEmbeddedConwayTerms.size());

    if (i > 0) {
      const auto &prev = kEmbeddedConwayEntries[i - 1];
      const bool strictly_sorted =
          (prev.characteristic < entry.characteristic) ||
          (prev.characteristic == entry.characteristic &&
           prev.degree < entry.degree);
      EXPECT_TRUE(strictly_sorted)
          << "Embedded table not sorted at index " << i;
    }

    auto span_opt =
        LookupEmbeddedConwayPoly(entry.characteristic, entry.degree);
    ASSERT_TRUE(span_opt.has_value());
    EXPECT_EQ(span_opt->size(), entry.count);
    EXPECT_EQ(span_opt->front().degree, 0U)
        << "Zero constant term for GF(" << entry.characteristic << "^"
        << entry.degree << ")";
    EXPECT_EQ(span_opt->back().degree, entry.degree);
    EXPECT_EQ(span_opt->back().coeff, 1U)
        << "Non-monic Conway poly for GF(" << entry.characteristic << "^"
        << entry.degree << ")";
  }
}

TEST(DbIntegrityTest, EmbeddedConwayBinaryIrreducibleAndPrimitive) {
  auto gf2 = std::make_shared<GF2>();

  // Exhaustive for d = 1..16.
  for (uint32_t d = 1; d <= 16; ++d) {
    const auto terms = databases::detail::LookupEmbeddedConwayPoly(2, d);
    ASSERT_TRUE(terms.has_value()) << "Missing embedded C_{2," << d << "}";
    const auto poly = BuildBinaryPoly(gf2, d, *terms);
    EXPECT_EQ(poly.Degree(), static_cast<int>(d));
    EXPECT_TRUE(poly.IsMonic());
    EXPECT_TRUE(utils::IsIrreducible(poly))
        << "C_{2," << d << "} failed IsIrreducible";
    EXPECT_TRUE(utils::IsPrimitive(poly))
        << "C_{2," << d << "} failed IsPrimitive";
  }

  // Sampled higher degrees up to 32.
  const std::vector<uint32_t> sampled_degrees = {17, 18, 20, 23,
                                                 24, 28, 31, 32};
  for (uint32_t d : sampled_degrees) {
    const auto terms = databases::detail::LookupEmbeddedConwayPoly(2, d);
    ASSERT_TRUE(terms.has_value()) << "Missing embedded C_{2," << d << "}";
    const auto poly = BuildBinaryPoly(gf2, d, *terms);
    EXPECT_EQ(poly.Degree(), static_cast<int>(d));
    EXPECT_TRUE(poly.IsMonic());
    EXPECT_TRUE(utils::IsIrreducible(poly))
        << "C_{2," << d << "} failed IsIrreducible";
    EXPECT_TRUE(utils::IsPrimitive(poly))
        << "C_{2," << d << "} failed IsPrimitive";
  }
}

TEST(DbIntegrityTest, EmbeddedConwayOddPrimesIrreducibleAndPrimitive) {
  std::unordered_map<uint32_t, std::shared_ptr<GFP<uint16_t>>> field_cache;
  size_t verified_count = 0;

  for (const auto &entry : databases::detail::kEmbeddedConwayEntries) {
    if (entry.characteristic == 2) {
      continue;
    }
    uint64_t order = 1;
    bool exceeds_limit = false;
    for (uint32_t i = 0; i < entry.degree; ++i) {
      order *= entry.characteristic;
      if (order > 10000ULL) {
        exceeds_limit = true;
        break;
      }
    }
    if (exceeds_limit) {
      continue;
    }

    auto &gfp = field_cache[entry.characteristic];
    if (!gfp) {
      gfp = std::make_shared<GFP<uint16_t>>(
          static_cast<uint16_t>(entry.characteristic));
    }

    const auto terms = databases::detail::LookupEmbeddedConwayPoly(
        entry.characteristic, entry.degree);
    ASSERT_TRUE(terms.has_value());
    const auto poly = BuildPrimePoly(gfp, entry.degree, *terms);
    EXPECT_EQ(poly.Degree(), static_cast<int>(entry.degree));
    EXPECT_TRUE(poly.IsMonic());
    EXPECT_TRUE(utils::IsIrreducible(poly))
        << "C_{" << entry.characteristic << "," << entry.degree
        << "} failed IsIrreducible";
    EXPECT_TRUE(utils::IsPrimitive(poly))
        << "C_{" << entry.characteristic << "," << entry.degree
        << "} failed IsPrimitive";
    ++verified_count;
  }

  EXPECT_GT(verified_count, 25U);
}

// ============================================================================
// Disk-Backed Conway & Irreducible Polynomial Databases
// ============================================================================

TEST(DbIntegrityTest, ConwayPolyDatabaseSampledEntriesAreIrreducible) {
  databases::ConwayPolyDatabase db;
  auto gf2 = std::make_shared<GF2>();

  // Binary Conway entries via ConwayPolyDatabase.
  for (int d : {2, 4, 8, 12, 16, 20, 24}) {
    const auto res = db.Fetch(2, d);
    const auto poly = BuildBinaryPolyFromResult(gf2, d, res);
    EXPECT_EQ(poly.Degree(), d);
    EXPECT_TRUE(utils::IsIrreducible(poly))
        << "ConwayPolyDatabase(2, " << d << ") not irreducible";
    EXPECT_TRUE(utils::IsPrimitive(poly))
        << "ConwayPolyDatabase(2, " << d << ") not primitive";
  }

  // Disk-backed odd primes outside the embedded table (p > 31) and in-table
  // primes.
  const std::vector<std::pair<int, int>> prime_samples = {
      {3, 3},  {5, 3},  {7, 3},  {11, 2}, {13, 3}, {37, 2},  {37, 3},
      {41, 2}, {43, 2}, {47, 2}, {53, 2}, {97, 2}, {101, 2},
  };
  for (const auto &[p, d] : prime_samples) {
    auto gfp = std::make_shared<GFP<uint16_t>>(static_cast<uint16_t>(p));
    const auto res = db.Fetch(p, d);
    const auto poly = BuildPrimePolyFromResult(gfp, d, res);
    EXPECT_EQ(poly.Degree(), d);
    EXPECT_TRUE(poly.IsMonic());
    EXPECT_TRUE(utils::IsIrreducible(poly))
        << "ConwayPolyDatabase(" << p << ", " << d << ") not irreducible";
    EXPECT_TRUE(utils::IsPrimitive(poly))
        << "ConwayPolyDatabase(" << p << ", " << d << ") not primitive";
  }
}

TEST(DbIntegrityTest, IrreduciblePolyDatabaseSampledEntriesAreIrreducible) {
  databases::IrreduciblePolyDatabase db;
  auto gf2 = std::make_shared<GF2>();

  const std::vector<int> degrees = {2,  3,  4,  5,  6,  7,  8,
                                    10, 12, 16, 20, 24, 31, 32};
  size_t total_checked = 0;
  for (int d : degrees) {
    const auto first = db.Fetch(2, d);
    const auto first_poly = BuildBinaryPolyFromResult(gf2, d, first);
    EXPECT_EQ(first_poly.Degree(), d);
    EXPECT_TRUE(utils::IsIrreducible(first_poly))
        << "IrreduciblePolyDatabase::Fetch(2, " << d << ") not irreducible";

    const auto all = db.FetchAll(2, d);
    ASSERT_FALSE(all.empty());
    const size_t limit = std::min<size_t>(all.size(), 5);
    for (size_t i = 0; i < limit; ++i) {
      const auto poly = BuildBinaryPolyFromResult(gf2, d, all[i]);
      EXPECT_EQ(poly.Degree(), d);
      EXPECT_TRUE(utils::IsIrreducible(poly))
          << "IrreduciblePolyDatabase::FetchAll(2, " << d << ")[" << i
          << "] not irreducible";
      ++total_checked;
    }
  }
  EXPECT_GE(total_checked, 40U);
}

// ============================================================================
// PrimeFactorsDatabase Exhaustive & Regression Integrity
// ============================================================================

TEST(DbIntegrityTest, PrimeFactorsDatabase500NumbersAndRegressions) {
  databases::PrimeFactorsDatabase db;

  // 1. Explicitly verify previously corrupted regression entries.
  struct RegressionCase {
    uint64_t n;
    std::vector<uint64_t> factors;
    std::vector<int> mults;
  };
  const std::vector<RegressionCase> regressions = {
      {121, {11}, {2}},
      {21218, {2, 103}, {1, 2}},
      {643482, {2, 3, 7, 5107}, {1, 2, 1, 1}},
      {309550, {2, 5, 41, 151}, {1, 2, 1, 1}},
      {823543, {7}, {7}},
  };
  for (const auto &rc : regressions) {
    VerifyPrimeFactorization(db, rc.n);
    const auto res = db.Fetch(rc.n);
    EXPECT_EQ(res.factors, rc.factors) << "n=" << rc.n;
    EXPECT_EQ(res.multiplicities, rc.mults) << "n=" << rc.n;
    EXPECT_EQ(res.composite, 1) << "n=" << rc.n;
  }

  // 2. Collect >= 500 diverse numbers: consecutive small integers, squares,
  //    prime powers, and deterministic random large integers.
  std::vector<uint64_t> numbers;
  numbers.reserve(550);

  // Consecutive integers 2..201 (200 numbers).
  for (uint64_t n = 2; n <= 201; ++n) {
    numbers.push_back(n);
  }

  // Squares of primes and composites (50 numbers).
  for (uint64_t k = 15; k < 65; ++k) {
    numbers.push_back(k * k);
  }

  // Prime powers p^e (50 numbers).
  const std::vector<uint64_t> base_primes = {2,  3,  5,  7,  11,
                                             13, 17, 19, 23, 29};
  for (uint64_t p : base_primes) {
    uint64_t val = p;
    for (int e = 1; e <= 5; ++e) {
      numbers.push_back(val);
      val *= p;
    }
  }

  // Deterministic random integers across wide range (205 numbers).
  std::mt19937_64 rng(0xDB123456789ULL);
  std::uniform_int_distribution<uint64_t> dist(202ULL, 50000000ULL);
  for (size_t i = 0; i < 205; ++i) {
    numbers.push_back(dist(rng));
  }

  ASSERT_GE(numbers.size(), 500U);
  for (uint64_t n : numbers) {
    VerifyPrimeFactorization(db, n);
  }
}

}  // namespace
}  // namespace xg
