#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

#include "xgalois/databases/embedded_conway_polys.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

using namespace xg;
using namespace xg::utils;

namespace {

int Mobius(uint64_t n) {
  auto factors = PrimeFactorizeWithMultiplicity(n);
  for (const auto &[p, mult] : factors) {
    (void)p;
    if (mult > 1) {
      return 0;
    }
  }
  return (factors.size() % 2 == 0) ? 1 : -1;
}

// Exact necklace formula N_q(n) = (1/n) * sum_{d | n} mu(n/d) * q^d.
uint64_t NecklaceIrreducibleCount(uint64_t q, uint32_t n) {
  int64_t sum = 0;
  for (uint32_t d = 1; d <= n; ++d) {
    if (n % d == 0) {
      sum += static_cast<int64_t>(Mobius(n / d)) *
             static_cast<int64_t>(SafeIntegerPower(q, d));
    }
  }
  return static_cast<uint64_t>(sum / static_cast<int64_t>(n));
}

// Exact number of monic primitive polynomials of degree n over F_q:
// phi(q^n - 1) / n.
uint64_t PrimitivePolyCount(uint64_t q, uint32_t n) {
  const uint64_t group_order = SafeIntegerPower(q, n) - 1;
  return EulerTotient(group_order) / static_cast<uint64_t>(n);
}

// Enumerates all elements of a small finite field in raw internal storage.
template <typename Field>
std::vector<typename Field::element_type> AllFieldElementsRaw(
    const std::shared_ptr<Field> &field) {
  using Raw = typename Field::element_type;
  const uint64_t q = field->Order();
  std::vector<Raw> elems;
  elems.reserve(static_cast<size_t>(q));
  if constexpr (std::is_integral_v<Raw>) {
    for (uint64_t v = 0; v < q; ++v) {
      elems.push_back(field->SetElementValue(static_cast<Raw>(v)));
    }
  } else {
    for (uint64_t v = 0; v < q; ++v) {
      elems.push_back(field->SetElementValue(std::to_string(v)));
    }
  }
  return elems;
}

template <typename Field>
void VerifyExhaustiveCounts(const std::shared_ptr<Field> &field, uint32_t deg) {
  using Raw = typename Field::element_type;
  const auto all_elems = AllFieldElementsRaw(field);
  const uint64_t q = field->Order();
  const uint64_t total_monic = SafeIntegerPower(q, deg);

  uint64_t irred_count = 0;
  uint64_t prim_count = 0;
  std::vector<Raw> coeffs(static_cast<size_t>(deg) + 1);
  coeffs[deg] = field->MultiplicativeIdentity();

  for (uint64_t idx = 0; idx < total_monic; ++idx) {
    uint64_t tmp = idx;
    for (uint32_t j = 0; j < deg; ++j) {
      coeffs[j] = all_elems[static_cast<size_t>(tmp % q)];
      tmp /= q;
    }
    auto poly = PolynomialDense<Field>::FromRaw(coeffs, field);
    const bool irred = IsIrreducible(poly);
    const bool prim = IsPrimitive(poly);
    if (prim) {
      EXPECT_TRUE(irred) << "Primitive polynomial must be irreducible: "
                         << poly.ToString();
    }
    if (irred) {
      ++irred_count;
    }
    if (prim) {
      ++prim_count;
    }
  }

  EXPECT_EQ(irred_count, NecklaceIrreducibleCount(q, deg))
      << "Irreducible count mismatch for q=" << q << ", deg=" << deg;
  EXPECT_EQ(prim_count, PrimitivePolyCount(q, deg))
      << "Primitive count mismatch for q=" << q << ", deg=" << deg;
}

template <typename Field>
void VerifyFactorization(const PolynomialDense<Field> &poly) {
  auto [lc, factors] = Factor(poly);
  auto field = poly.Field();

  EXPECT_EQ(lc, poly.LeadingCoeff());

  auto reconstructed = PolynomialDense<Field>::FromRaw({lc.RawValue()}, field,
                                                       poly.GetVariable());
  for (size_t i = 0; i < factors.size(); ++i) {
    const auto &[f, mult] = factors[i];
    EXPECT_GE(mult, 1u);
    EXPECT_GE(f.Degree(), 1);
    EXPECT_TRUE(f.IsMonic()) << "Factor not monic: " << f.ToString();
    EXPECT_TRUE(IsIrreducible(f)) << "Factor not irreducible: " << f.ToString();
    for (size_t j = i + 1; j < factors.size(); ++j) {
      EXPECT_NE(f, factors[j].first)
          << "Duplicate irreducible factor emitted: " << f.ToString();
    }
    reconstructed *= f.Pow(mult);
  }

  EXPECT_EQ(reconstructed, poly)
      << "Factorization product did not reconstruct original polynomial: "
      << poly.ToString();
}

}  // namespace

TEST(PolyAlgoTest, ExhaustiveIrreducibleAndPrimitiveCountsGF2) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  for (uint32_t deg = 1; deg <= 8; ++deg) {
    VerifyExhaustiveCounts(gf2, deg);
  }
}

TEST(PolyAlgoTest, ExhaustiveIrreducibleAndPrimitiveCountsGF3AndGF5) {
  auto gf3 = std::make_shared<GaloisFieldPrime<uint8_t>>(3);
  for (uint32_t deg = 1; deg <= 4; ++deg) {
    VerifyExhaustiveCounts(gf3, deg);
  }

  auto gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);
  for (uint32_t deg = 1; deg <= 3; ++deg) {
    VerifyExhaustiveCounts(gf5, deg);
  }
}

TEST(PolyAlgoTest, ExhaustiveIrreducibleAndPrimitiveCountsExtensionFields) {
  // GF(2^2) = F_4
  auto gf4 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(2);
  for (uint32_t deg = 1; deg <= 3; ++deg) {
    VerifyExhaustiveCounts(gf4, deg);
  }

  // GF(2^2) Zech-log table
  auto gf4_tbl = std::make_shared<GFBEZechLogTables<uint8_t>>(2);
  VerifyExhaustiveCounts(gf4_tbl, 1);
  VerifyExhaustiveCounts(gf4_tbl, 2);

  // GF(3^2) = F_9 (GaloisFieldExtension)
  auto gf9 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  VerifyExhaustiveCounts(gf9, 1);
  VerifyExhaustiveCounts(gf9, 2);
}

TEST(PolyAlgoTest, EmbeddedConwayPolynomialsAreIrreducibleAndPrimitive) {
  // GF(2) degrees 1..16
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  for (uint32_t deg = 1; deg <= 16; ++deg) {
    auto terms = databases::detail::LookupEmbeddedConwayPoly(2, deg);
    ASSERT_TRUE(terms.has_value());
    std::vector<uint8_t> raw(static_cast<size_t>(deg) + 1, 0);
    for (const auto &t : *terms) {
      raw[t.degree] = static_cast<uint8_t>(t.coeff);
    }
    auto poly = PolynomialDense<GaloisFieldBinary>::FromRaw(raw, gf2);
    EXPECT_TRUE(IsIrreducible(poly)) << "GF(2) Conway deg=" << deg;
    EXPECT_TRUE(IsPrimitive(poly)) << "GF(2) Conway deg=" << deg;
  }

  // Odd prime fields p in {3, 5, 7}
  for (uint32_t p : {3u, 5u, 7u}) {
    auto gfp =
        std::make_shared<GaloisFieldPrime<uint8_t>>(static_cast<uint8_t>(p));
    const uint32_t max_deg = (p == 3) ? 6u : ((p == 5) ? 5u : 4u);
    for (uint32_t deg = 1; deg <= max_deg; ++deg) {
      auto terms = databases::detail::LookupEmbeddedConwayPoly(p, deg);
      ASSERT_TRUE(terms.has_value());
      std::vector<uint8_t> raw(static_cast<size_t>(deg) + 1, 0);
      for (const auto &t : *terms) {
        raw[t.degree] = static_cast<uint8_t>(t.coeff);
      }
      auto poly = PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw(raw, gfp);
      EXPECT_TRUE(IsIrreducible(poly)) << "GF(" << p << ") Conway deg=" << deg;
      EXPECT_TRUE(IsPrimitive(poly)) << "GF(" << p << ") Conway deg=" << deg;
    }
  }
}

TEST(PolyAlgoTest, FindIrreducibleAndPrimitiveAcrossAllFieldTypes) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  for (uint32_t deg : {1u, 2u, 4u, 7u, 11u}) {
    auto irred = FindIrreduciblePolynomial(gf2, deg);
    EXPECT_EQ(irred.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(irred.IsMonic());
    EXPECT_TRUE(IsIrreducible(irred));

    auto prim = FindPrimitivePolynomial(gf2, deg);
    EXPECT_EQ(prim.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(prim.IsMonic());
    EXPECT_TRUE(IsPrimitive(prim));
  }

  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  for (uint32_t deg : {1u, 2u, 3u, 5u}) {
    auto irred = FindIrreduciblePolynomial(gf7, deg);
    EXPECT_EQ(irred.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(IsIrreducible(irred));

    auto prim = FindPrimitivePolynomial(gf7, deg);
    EXPECT_EQ(prim.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(IsPrimitive(prim));
  }

  auto gf2_4 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4);
  for (uint32_t deg : {1u, 2u, 3u, 4u}) {
    auto irred = FindIrreduciblePolynomial(gf2_4, deg);
    EXPECT_EQ(irred.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(IsIrreducible(irred));

    auto prim = FindPrimitivePolynomial(gf2_4, deg);
    EXPECT_EQ(prim.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(IsPrimitive(prim));
  }

  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  for (uint32_t deg : {1u, 2u, 3u}) {
    auto irred = FindIrreduciblePolynomial(gf3_2, deg);
    EXPECT_EQ(irred.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(IsIrreducible(irred));

    auto prim = FindPrimitivePolynomial(gf3_2, deg);
    EXPECT_EQ(prim.Degree(), static_cast<int>(deg));
    EXPECT_TRUE(IsPrimitive(prim));
  }
}

TEST(PolyAlgoTest, FactorOverGF2IncludingChar2RepeatedMultiplicities) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  EXPECT_THROW(Factor(PolynomialDense<GaloisFieldBinary>::Zero(gf2)),
               std::invalid_argument);

  auto [lc_one, empty_factors] =
      Factor(PolynomialDense<GaloisFieldBinary>::One(gf2));
  EXPECT_TRUE(lc_one.IsOne());
  EXPECT_TRUE(empty_factors.empty());

  auto f1 = ParsePolynomial(gf2, "x + 1");
  auto f2 = ParsePolynomial(gf2, "x^2 + x + 1");
  auto f3a = ParsePolynomial(gf2, "x^3 + x + 1");
  auto f3b = ParsePolynomial(gf2, "x^3 + x^2 + 1");

  // Product with multiplicities 1, 2 (= p), 3, 4 (= p^2) and equal-degree pair.
  auto composite = f1.Pow(4) * f2.Pow(2) * f3a * f3b.Pow(2);
  VerifyFactorization(composite);

  // Pure p^2 power: (x^3 + x + 1)^4
  VerifyFactorization(f3a.Pow(4));

  // Random degree 1..12 polynomials over GF(2)
  std::mt19937_64 rng(0xCAFEBABE1ULL);
  for (int trial = 0; trial < 30; ++trial) {
    const size_t deg = 1 + (trial % 12);
    std::vector<uint8_t> raw(deg + 1);
    for (size_t i = 0; i < deg; ++i) raw[i] = gf2->Random(rng);
    raw[deg] = 1;
    VerifyFactorization(PolynomialDense<GaloisFieldBinary>::FromRaw(raw, gf2));
  }
}

TEST(PolyAlgoTest, FactorOverOddPrimeFieldsGF3GF5GF7) {
  std::mt19937_64 rng(0xDEADBEEF2ULL);
  for (uint8_t p : {uint8_t{3}, uint8_t{5}, uint8_t{7}}) {
    auto gfp = std::make_shared<GaloisFieldPrime<uint8_t>>(p);
    auto x = PolynomialDense<GaloisFieldPrime<uint8_t>>::Monomial(gfp, 1, 1);

    // Explicit multiplicity p and 2p factors + non-monic leading coefficient.
    auto irred1 = FindIrreduciblePolynomial(gfp, 1);
    auto irred2 = FindIrreduciblePolynomial(gfp, 2);
    auto lead = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(p - 1, gfp);
    auto multi = (irred1.Pow(p) * irred2.Pow(p) *
                  (x + PolynomialDense<GaloisFieldPrime<uint8_t>>::One(gfp))) *
                 lead;
    VerifyFactorization(multi);

    // x^p - x splits into all p linear factors.
    auto split_all = x.Pow(p) - x;
    auto [lc, factors] = Factor(split_all);
    EXPECT_TRUE(lc.IsOne());
    EXPECT_EQ(factors.size(), static_cast<size_t>(p));
    VerifyFactorization(split_all);

    // Random polynomials of degree 1..10.
    for (int trial = 0; trial < 20; ++trial) {
      const size_t deg = 1 + (trial % 10);
      std::vector<uint8_t> raw(deg + 1);
      for (size_t i = 0; i < deg; ++i) raw[i] = gfp->Random(rng);
      raw[deg] = static_cast<uint8_t>(1 + (rng() % (p - 1)));
      VerifyFactorization(
          PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw(raw, gfp));
    }
  }
}

TEST(PolyAlgoTest, FactorOverExtensionFieldsGF2_4AndGF3_2) {
  std::mt19937_64 rng(0xF00DBABE3ULL);

  // GF(2^4) direct and Zech-log table
  auto gf2_4 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4);
  auto gf2_4_tbl = std::make_shared<GFBEZechLogTables<uint16_t>>(4);

  auto check_bin_ext = [&](auto field) {
    using Field = typename decltype(field)::element_type;
    using Elem = GaloisFieldElement<Field>;
    auto x = PolynomialDense<Field>::Monomial(field, 1,
                                              field->MultiplicativeIdentity());
    auto c1 = Elem(3, field);
    auto c2 = Elem(7, field);
    auto c3 = Elem(12, field);
    // (x + c1)^2 * (x + c2)^4 * (x + c3) tests inverse Frobenius on p-th roots!
    auto l1 = x + PolynomialDense<Field>::FromRaw({c1.RawValue()}, field);
    auto l2 = x + PolynomialDense<Field>::FromRaw({c2.RawValue()}, field);
    auto l3 = x + PolynomialDense<Field>::FromRaw({c3.RawValue()}, field);
    auto irred2 = FindIrreduciblePolynomial(field, 2);
    auto poly = (l1.Pow(2) * l2.Pow(4) * l3 * irred2.Pow(2)) * c2;
    VerifyFactorization(poly);

    for (int trial = 0; trial < 12; ++trial) {
      const size_t deg = 2 + (trial % 6);
      std::vector<uint16_t> raw(deg + 1);
      for (size_t i = 0; i < deg; ++i) raw[i] = field->Random(rng);
      raw[deg] = field->MultiplicativeIdentity();
      VerifyFactorization(PolynomialDense<Field>::FromRaw(raw, field));
    }
  };

  check_bin_ext(gf2_4);
  check_bin_ext(gf2_4_tbl);

  // GF(3^2) GaloisFieldExtension
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  using ExtField = GaloisFieldExtension<uint32_t>;
  using ExtElem = GaloisFieldElement<ExtField>;
  auto x = PolynomialDense<ExtField>::Monomial(gf3_2, 1,
                                               gf3_2->MultiplicativeIdentity());
  ExtElem a("α + 1", gf3_2);
  ExtElem b("2α + 2", gf3_2);
  auto l1 = x - PolynomialDense<ExtField>::FromRaw({a.RawValue()}, gf3_2);
  auto l2 = x - PolynomialDense<ExtField>::FromRaw({b.RawValue()}, gf3_2);
  auto irred2 = FindIrreduciblePolynomial(gf3_2, 2);
  // Multiplicity 3 (= p) tests inverse Frobenius a^(p^(m-1)) = a^3 in GF(3^2)!
  auto ext_poly = (l1.Pow(3) * l2 * irred2.Pow(3)) * b;
  VerifyFactorization(ext_poly);

  for (int trial = 0; trial < 10; ++trial) {
    const size_t deg = 2 + (trial % 5);
    std::vector<ExtensionElement<uint32_t>> raw(deg + 1);
    for (size_t i = 0; i < deg; ++i) raw[i] = gf3_2->Random(rng);
    raw[deg] = gf3_2->MultiplicativeIdentity();
    VerifyFactorization(PolynomialDense<ExtField>::FromRaw(raw, gf3_2));
  }
}

TEST(PolyAlgoTest, FindRootsAndChienSearchAcrossFields) {
  // Small prime field GF(7): compare FindRoots, ChienSearch, and brute force.
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  auto p7 = ParsePolynomial(gf7, "3x^4 + 2x^3 + 6x^2 + x + 2");
  auto roots7 = FindRoots(p7);
  std::vector<GaloisFieldElement<GaloisFieldPrime<uint8_t>>> all7;
  for (uint8_t v = 0; v < 7; ++v) all7.emplace_back(v, gf7);
  auto chien7 = ChienSearch(p7, all7);
  EXPECT_EQ(roots7, chien7);
  for (const auto &r : roots7) {
    EXPECT_TRUE(p7(r).IsZero());
  }

  // Large prime field GF(1000003): Cantor-Zassenhaus finds exact roots quickly.
  auto gf_big = std::make_shared<GaloisFieldPrime<uint32_t>>(1000003);
  using BigElem = GaloisFieldElement<GaloisFieldPrime<uint32_t>>;
  auto x_big =
      PolynomialDense<GaloisFieldPrime<uint32_t>>::Monomial(gf_big, 1, 1);
  std::vector<uint32_t> expected_vals = {0, 42, 12345, 500001, 999999};
  auto prod_big = PolynomialDense<GaloisFieldPrime<uint32_t>>::One(gf_big);
  for (uint32_t r : expected_vals) {
    prod_big *= (x_big - PolynomialDense<GaloisFieldPrime<uint32_t>>::FromRaw(
                             {r}, gf_big));
  }
  // Multiply by an irreducible quadratic so polynomial doesn't split
  // completely.
  auto irred_big = FindIrreduciblePolynomial(gf_big, 2);
  prod_big *= irred_big;

  auto found_big = FindRoots(prod_big);
  ASSERT_EQ(found_big.size(), expected_vals.size());
  for (size_t i = 0; i < expected_vals.size(); ++i) {
    EXPECT_EQ(found_big[i], BigElem(expected_vals[i], gf_big));
    EXPECT_TRUE(prod_big(found_big[i]).IsZero());
  }

  // GF(2^4) Zech-log table
  auto gf2_4 = std::make_shared<GFBEZechLogTables<uint16_t>>(4);
  using BinElem = GaloisFieldElement<GFBEZechLogTables<uint16_t>>;
  auto x_bin = PolynomialDense<GFBEZechLogTables<uint16_t>>::Monomial(
      gf2_4, 1, gf2_4->MultiplicativeIdentity());
  std::vector<uint16_t> bin_roots = {1, 5, 9, 14};
  auto prod_bin = PolynomialDense<GFBEZechLogTables<uint16_t>>::One(gf2_4);
  for (uint16_t r : bin_roots) {
    prod_bin *= (x_bin - PolynomialDense<GFBEZechLogTables<uint16_t>>::FromRaw(
                             {gf2_4->SetElementValue(r)}, gf2_4));
  }
  prod_bin *= FindIrreduciblePolynomial(gf2_4, 2);
  auto found_bin = FindRoots(prod_bin);
  ASSERT_EQ(found_bin.size(), bin_roots.size());
  for (size_t i = 0; i < bin_roots.size(); ++i) {
    EXPECT_EQ(found_bin[i], BinElem(bin_roots[i], gf2_4));
  }
  auto all_raw_bin = AllFieldElementsRaw(gf2_4);
  auto chien_raw =
      ChienSearch(prod_bin, std::span<const uint16_t>(all_raw_bin));
  EXPECT_EQ(chien_raw.size(), bin_roots.size());

  // GF(3^2) GaloisFieldExtension
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  using ExtElem = GaloisFieldElement<GaloisFieldExtension<uint32_t>>;
  auto x_ext = PolynomialDense<GaloisFieldExtension<uint32_t>>::Monomial(
      gf3_2, 1, gf3_2->MultiplicativeIdentity());
  ExtElem r1("0", gf3_2);
  ExtElem r2("α + 1", gf3_2);
  ExtElem r3("2α + 2", gf3_2);
  auto prod_ext =
      (x_ext - PolynomialDense<GaloisFieldExtension<uint32_t>>::FromRaw(
                   {r1.RawValue()}, gf3_2)) *
      (x_ext - PolynomialDense<GaloisFieldExtension<uint32_t>>::FromRaw(
                   {r2.RawValue()}, gf3_2)) *
      (x_ext - PolynomialDense<GaloisFieldExtension<uint32_t>>::FromRaw(
                   {r3.RawValue()}, gf3_2)) *
      FindIrreduciblePolynomial(gf3_2, 2);
  auto found_ext = FindRoots(prod_ext);
  ASSERT_EQ(found_ext.size(), 3u);
  for (const auto &r : found_ext) {
    EXPECT_TRUE(prod_ext(r).IsZero());
  }
}

TEST(PolyAlgoTest, InterpolatePropertyTestsAcrossFields) {
  std::mt19937_64 rng(0x55AA55AA4ULL);

  // GF(7)
  auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
  using Elem7 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>;
  {
    std::vector<Elem7> xs = {Elem7(0, gf7), Elem7(2, gf7), Elem7(4, gf7),
                             Elem7(6, gf7), Elem7(1, gf7)};
    std::vector<Elem7> ys = {Elem7(3, gf7), Elem7(1, gf7), Elem7(5, gf7),
                             Elem7(0, gf7), Elem7(6, gf7)};
    auto poly = Interpolate(xs, ys);
    EXPECT_LT(poly.Degree(), static_cast<int>(xs.size()));
    for (size_t i = 0; i < xs.size(); ++i) {
      EXPECT_EQ(poly(xs[i]), ys[i]);
    }

    // Error cases
    EXPECT_THROW(Interpolate(std::vector<Elem7>{}, std::vector<Elem7>{}),
                 std::invalid_argument);
    EXPECT_THROW(Interpolate(xs, std::vector<Elem7>{Elem7(1, gf7)}),
                 std::invalid_argument);
    std::vector<Elem7> dup_x = {Elem7(1, gf7), Elem7(2, gf7), Elem7(1, gf7)};
    std::vector<Elem7> dup_y = {Elem7(0, gf7), Elem7(1, gf7), Elem7(2, gf7)};
    EXPECT_THROW(Interpolate(dup_x, dup_y), std::invalid_argument);
  }

  // GF(2^4) exact recovery of degree-6 polynomial from 7 evaluation points
  auto gf2_4 = std::make_shared<GFBEZechLogTables<uint16_t>>(4);
  {
    std::vector<uint16_t> orig_raw(7);
    for (size_t i = 0; i < 7; ++i) orig_raw[i] = gf2_4->Random(rng);
    orig_raw[6] = gf2_4->MultiplicativeIdentity();
    auto orig =
        PolynomialDense<GFBEZechLogTables<uint16_t>>::FromRaw(orig_raw, gf2_4);

    std::vector<GaloisFieldElement<GFBEZechLogTables<uint16_t>>> xs, ys;
    for (uint16_t v = 0; v < 7; ++v) {
      GaloisFieldElement<GFBEZechLogTables<uint16_t>> pt(v, gf2_4);
      xs.push_back(pt);
      ys.push_back(orig(pt));
    }
    auto recovered = Interpolate(xs, ys);
    EXPECT_EQ(recovered, orig);
  }

  // GF(3^2) GaloisFieldExtension exact recovery of degree-4 polynomial
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  {
    auto all_raw = AllFieldElementsRaw(gf3_2);
    std::vector<ExtensionElement<uint32_t>> orig_raw(5);
    for (size_t i = 0; i < 5; ++i) orig_raw[i] = gf3_2->Random(rng);
    orig_raw[4] = gf3_2->MultiplicativeIdentity();
    auto orig = PolynomialDense<GaloisFieldExtension<uint32_t>>::FromRaw(
        orig_raw, gf3_2);

    std::vector<GaloisFieldElement<GaloisFieldExtension<uint32_t>>> xs, ys;
    for (size_t i = 0; i < 5; ++i) {
      auto pt =
          GaloisFieldElement<GaloisFieldExtension<uint32_t>>::FromInternal(
              all_raw[i], gf3_2);
      xs.push_back(pt);
      ys.push_back(orig(pt));
    }
    auto recovered = Interpolate(xs, ys);
    EXPECT_EQ(recovered, orig);
    EXPECT_LT(recovered.Degree(), 5);
    for (size_t i = 0; i < 5; ++i) {
      EXPECT_EQ(recovered(xs[i]), ys[i]);
    }
  }
}
