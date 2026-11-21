#include "xgalois/field/gf_extension.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/utils/math.hpp"

using namespace xg;

static_assert(FiniteField<GaloisFieldExtension<uint8_t>>);
static_assert(FiniteField<GaloisFieldExtension<uint16_t>>);
static_assert(FiniteField<GFPX<uint32_t>>);

namespace {

template <typename Field>
void VerifyFieldAxiomsExhaustive(const Field &field) {
  using Elem = typename Field::element_type;
  const uint64_t q = field.Order();
  const uint64_t p = field.Characteristic();
  const uint32_t m = field.Degree();
  const Elem zero = field.AdditiveIdentity();
  const Elem one = field.MultiplicativeIdentity();

  std::vector<Elem> elems;
  elems.reserve(static_cast<size_t>(q));
  for (uint64_t i = 0; i < q; ++i) {
    Elem e = field.FromInt(i);
    EXPECT_EQ(field.ToInt(e), i);
    EXPECT_LT(e.Degree(), static_cast<int>(m));
    elems.push_back(e);
  }

  // Distinctness and strict ordering.
  for (size_t i = 0; i + 1 < elems.size(); ++i) {
    EXPECT_LT(elems[i], elems[i + 1]);
    EXPECT_NE(elems[i], elems[i + 1]);
  }

  for (const Elem &a : elems) {
    EXPECT_EQ(field.Add(a, zero), a);
    EXPECT_EQ(field.Add(a, field.Neg(a)), zero);
    EXPECT_EQ(field.Mul(a, one), a);
    EXPECT_EQ(field.Mul(a, zero), zero);
    EXPECT_EQ(field.Pow(a, q), a);
    EXPECT_EQ(field.Pow(a, 0), one);
    EXPECT_EQ(field.Frobenius(a), field.Pow(a, p));

    // Trace and Norm lie in the prime subfield F_p (degree <= 0).
    const Elem tr = field.Trace(a);
    const Elem nm = field.Norm(a);
    EXPECT_LE(tr.Degree(), 0);
    EXPECT_LE(nm.Degree(), 0);

    // Minimal polynomial is monic, degree divides m, and annihilates a.
    const Elem min_poly = field.MinimalPolynomial(a);
    EXPECT_GE(min_poly.Degree(), 1);
    EXPECT_EQ(m % static_cast<uint32_t>(min_poly.Degree()), 0U);
    EXPECT_EQ(min_poly[static_cast<size_t>(min_poly.Degree())], 1);

    // Evaluate min_poly(a) in the extension field via Horner's method.
    Elem eval = zero;
    for (int d = min_poly.Degree(); d >= 0; --d) {
      eval = field.Add(field.Mul(eval, a),
                       field.FromCoeffs({min_poly[static_cast<size_t>(d)]}));
    }
    EXPECT_TRUE(field.IsZero(eval));

    if (!field.IsZero(a)) {
      const Elem inv_a = field.Inv(a);
      EXPECT_EQ(field.Mul(a, inv_a), one);
      EXPECT_EQ(field.Mul(inv_a, a), one);
      EXPECT_EQ(field.Pow(a, q - 1), one);
      EXPECT_EQ(field.PowSigned(a, -1), inv_a);
      EXPECT_EQ(field.PowSigned(a, -2), field.Mul(inv_a, inv_a));
    }

    for (const Elem &b : elems) {
      EXPECT_EQ(field.Add(a, b), field.Add(b, a));
      EXPECT_EQ(field.Sub(a, b), field.Add(a, field.Neg(b)));
      EXPECT_EQ(field.Mul(a, b), field.Mul(b, a));
      EXPECT_EQ(field.Frobenius(field.Add(a, b)),
                field.Add(field.Frobenius(a), field.Frobenius(b)));
      EXPECT_EQ(field.Frobenius(field.Mul(a, b)),
                field.Mul(field.Frobenius(a), field.Frobenius(b)));
      EXPECT_EQ(field.Trace(field.Add(a, b)),
                field.Add(field.Trace(a), field.Trace(b)));
      EXPECT_EQ(field.Norm(field.Mul(a, b)),
                field.Mul(field.Norm(a), field.Norm(b)));
      if (!field.IsZero(b)) {
        EXPECT_EQ(field.Div(a, b), field.Mul(a, field.Inv(b)));
      }
    }
  }

  // Associativity and distributivity on sampled/all triples (full for q <= 27,
  // strided for q == 49 to keep runtime brisk).
  const size_t step = (q <= 27) ? 1 : 3;
  for (size_t i = 0; i < elems.size(); i += step) {
    for (size_t j = 0; j < elems.size(); j += step) {
      for (size_t k = 0; k < elems.size(); k += step) {
        const Elem &a = elems[i];
        const Elem &b = elems[j];
        const Elem &c = elems[k];
        EXPECT_EQ(field.Add(field.Add(a, b), c), field.Add(a, field.Add(b, c)));
        EXPECT_EQ(field.Mul(field.Mul(a, b), c), field.Mul(a, field.Mul(b, c)));
        EXPECT_EQ(field.Mul(a, field.Add(b, c)),
                  field.Add(field.Mul(a, b), field.Mul(a, c)));
      }
    }
  }
}

template <typename Field>
void VerifySqrtAndLogExhaustive(const Field &field) {
  using Elem = typename Field::element_type;
  const uint64_t q = field.Order();
  const uint64_t p = field.Characteristic();
  const Elem zero = field.AdditiveIdentity();
  const Elem one = field.MultiplicativeIdentity();

  EXPECT_EQ(field.Sqrt(zero), zero);
  EXPECT_EQ(field.Sqrt(one), one);
  EXPECT_THROW(field.Log(zero), std::domain_error);

  const Elem gen = field.MultiplicativeGenerator();
  EXPECT_EQ(field.MultiplicativeOrder(gen), q - 1);

  uint64_t qr_count = 0;
  for (uint64_t i = 1; i < q; ++i) {
    const Elem a = field.FromInt(i);
    const uint64_t lg = field.Log(a);
    EXPECT_LT(lg, q - 1);
    EXPECT_EQ(field.Pow(gen, lg), a);
    EXPECT_EQ(field.Log(a, gen), lg);

    if (p == 2) {
      const Elem root = field.Sqrt(a);
      EXPECT_EQ(field.Mul(root, root), a);
      ++qr_count;
    } else {
      const bool is_qr = field.IsOne(field.Pow(a, (q - 1) / 2));
      if (is_qr) {
        const Elem root = field.Sqrt(a);
        EXPECT_EQ(field.Mul(root, root), a);
        ++qr_count;
      } else {
        EXPECT_THROW(field.Sqrt(a), std::domain_error);
      }
    }
  }

  if (p == 2) {
    EXPECT_EQ(qr_count, q - 1);
  } else {
    EXPECT_EQ(qr_count, (q - 1) / 2);
  }
}

}  // namespace

class GaloisFieldExtensionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf9 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3, 2));
    gf25 = std::make_shared<GFPX<uint8_t>>(std::make_pair(5, 2));
    gf27 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3, 3));
    gf49 = std::make_shared<GFPX<uint8_t>>(std::make_pair(7, 2));
  }

  std::shared_ptr<GFPX<uint8_t>> gf9, gf25, gf27, gf49;
};

TEST_F(GaloisFieldExtensionTest, FieldProperties) {
  EXPECT_EQ(gf9->Characteristic(), 3ULL);
  EXPECT_EQ(gf9->Order(), 9ULL);
  EXPECT_EQ(gf9->Degree(), 2U);

  EXPECT_EQ(gf25->Characteristic(), 5ULL);
  EXPECT_EQ(gf25->Order(), 25ULL);
  EXPECT_EQ(gf25->Degree(), 2U);

  EXPECT_EQ(gf27->Characteristic(), 3ULL);
  EXPECT_EQ(gf27->Order(), 27ULL);
  EXPECT_EQ(gf27->Degree(), 3U);

  EXPECT_EQ(gf49->Characteristic(), 7ULL);
  EXPECT_EQ(gf49->Order(), 49ULL);
  EXPECT_EQ(gf49->Degree(), 2U);
}

TEST_F(GaloisFieldExtensionTest, IdentityElements) {
  auto zero = gf9->AdditiveIdentity();
  auto one = gf9->MultiplicativeIdentity();

  EXPECT_TRUE(zero.IsZero());
  EXPECT_EQ(zero.Degree(), -1);
  EXPECT_EQ(zero[0], 0);
  EXPECT_TRUE(gf9->IsZero(zero));
  EXPECT_FALSE(gf9->IsOne(zero));

  EXPECT_FALSE(one.IsZero());
  EXPECT_EQ(one.Degree(), 0);
  EXPECT_EQ(one[0], 1);
  EXPECT_TRUE(gf9->IsOne(one));
  EXPECT_FALSE(gf9->IsZero(one));
}

TEST_F(GaloisFieldExtensionTest, ExhaustiveFieldAxiomsGF9) {
  VerifyFieldAxiomsExhaustive(*gf9);
  VerifySqrtAndLogExhaustive(*gf9);
}

TEST_F(GaloisFieldExtensionTest, ExhaustiveFieldAxiomsGF25) {
  VerifyFieldAxiomsExhaustive(*gf25);
  VerifySqrtAndLogExhaustive(*gf25);
}

TEST_F(GaloisFieldExtensionTest, ExhaustiveFieldAxiomsGF27) {
  VerifyFieldAxiomsExhaustive(*gf27);
  VerifySqrtAndLogExhaustive(*gf27);
}

TEST_F(GaloisFieldExtensionTest, ExhaustiveFieldAxiomsGF49) {
  VerifyFieldAxiomsExhaustive(*gf49);
  VerifySqrtAndLogExhaustive(*gf49);
}

TEST_F(GaloisFieldExtensionTest, CharacteristicTwoExtension) {
  GFPX<uint8_t> gf8(std::make_pair(2, 3));
  EXPECT_EQ(gf8.Characteristic(), 2ULL);
  EXPECT_EQ(gf8.Order(), 8ULL);
  EXPECT_EQ(gf8.Degree(), 3U);
  VerifyFieldAxiomsExhaustive(gf8);
  VerifySqrtAndLogExhaustive(gf8);
}

TEST_F(GaloisFieldExtensionTest, NonConstantPolynomialArithmetic) {
  // Explicit GF(3^2) = F_3[x] / (x^2 + 1), where alpha^2 = 2 = -1.
  GFPX<uint8_t> gf9_custom(std::make_pair(3, 2), "x^2 + 1", "poly", "α", true);

  auto alpha = gf9_custom.FromCoeffs({0, 1});
  auto a = gf9_custom.FromCoeffs({1, 2});  // 2*alpha + 1
  auto b = gf9_custom.FromCoeffs({2, 1});  // alpha + 2

  // (2a + 1) + (a + 2) = 3a + 3 = 0 mod 3
  EXPECT_TRUE(gf9_custom.IsZero(gf9_custom.Add(a, b)));

  // alpha^2 = -1 = 2
  auto alpha_sq = gf9_custom.Mul(alpha, alpha);
  EXPECT_EQ(alpha_sq.Degree(), 0);
  EXPECT_EQ(alpha_sq[0], 2);

  // (2a + 1) * (a + 2) = 2a^2 + 5a + 2 = 2(2) + 2a + 2 = 2a + 6 = 2a mod 3
  auto prod = gf9_custom.Mul(a, b);
  EXPECT_EQ(prod, gf9_custom.FromCoeffs({0, 2}));

  // Inverse of alpha: alpha * (2*alpha) = 2*alpha^2 = 2*2 = 1 mod 3
  EXPECT_EQ(gf9_custom.Inv(alpha), gf9_custom.FromCoeffs({0, 2}));

  // Minimal polynomial of alpha is x^2 + 1.
  EXPECT_EQ(gf9_custom.MinimalPolynomial(alpha), gf9_custom.Modulus());
}

TEST_F(GaloisFieldExtensionTest,
       MultiplicativeGeneratorsAndNonPrimitiveModulus) {
  // Conway polynomial: alpha = {0, 1} is primitive.
  auto conway_gen = gf9->MultiplicativeGenerator();
  EXPECT_EQ(conway_gen, ExtensionElement<uint8_t>({0, 1}));
  EXPECT_EQ(gf9->MultiplicativeOrder(conway_gen), 8U);

  auto all_gens_gf9 = gf9->MultiplicativeGenerators();
  EXPECT_EQ(all_gens_gf9.size(), utils::EulerTotient(8));
  for (const auto &g : all_gens_gf9) {
    EXPECT_EQ(gf9->MultiplicativeOrder(g), 8U);
  }

  auto all_gens_gf27 = gf27->MultiplicativeGenerators();
  EXPECT_EQ(all_gens_gf27.size(), utils::EulerTotient(26));
  for (const auto &g : all_gens_gf27) {
    EXPECT_EQ(gf27->MultiplicativeOrder(g), 26U);
  }

  // Non-primitive irreducible polynomial x^2 + 1 over GF(3): alpha has order 4,
  // so MultiplicativeGenerator() must find an element of full order 8.
  GFPX<uint8_t> gf9_non_prim(std::make_pair(3, 2), "x^2 + 1", "poly", "α",
                             true);
  auto alpha = gf9_non_prim.FromCoeffs({0, 1});
  EXPECT_EQ(gf9_non_prim.MultiplicativeOrder(alpha), 4U);

  auto gen = gf9_non_prim.MultiplicativeGenerator();
  EXPECT_EQ(gf9_non_prim.MultiplicativeOrder(gen), 8U);

  // Log with non-primitive base alpha: alpha^2 == 2 is in subgroup, but gen is
  // not.
  EXPECT_EQ(gf9_non_prim.Log(gf9_non_prim.FromCoeffs({2}), alpha), 2U);
  EXPECT_THROW(gf9_non_prim.Log(gen, alpha), std::domain_error);
}

TEST_F(GaloisFieldExtensionTest, IrreducibilityCheckAndValidation) {
  // Irreducible polynomials accepted.
  EXPECT_NO_THROW(
      GFPX<uint8_t>(std::make_pair(3, 2), "x^2 + 1", "poly", "α", true));
  EXPECT_NO_THROW(
      GFPX<uint8_t>(std::make_pair(3, 3), "x^3 + 2*x + 1", "poly", "α", true));
  EXPECT_NO_THROW(
      GFPX<uint8_t>(std::make_pair(5, 2), "x^2 + 2", "poly", "α", true));

  // Reducible polynomials rejected when check_irreducible == true.
  EXPECT_THROW(
      GFPX<uint8_t>(std::make_pair(3, 2), "x^2 + 2", "poly", "α", true),
      std::invalid_argument);
  EXPECT_THROW(GFPX<uint8_t>(std::make_pair(3, 2), "x^2", "poly", "α", true),
               std::invalid_argument);
  EXPECT_THROW(
      GFPX<uint8_t>(std::make_pair(5, 2), "x^2 + 1", "poly", "α", true),
      std::invalid_argument);
  EXPECT_THROW(
      GFPX<uint8_t>(std::make_pair(3, 3), "x^3 + 1", "poly", "α", true),
      std::invalid_argument);

  // Degree mismatch rejected even when check_irreducible == false.
  EXPECT_THROW(
      GFPX<uint8_t>(std::make_pair(3, 2), "x^3 + 2*x + 1", "poly", "α", false),
      std::invalid_argument);
  EXPECT_THROW(
      GFPX<uint8_t>(std::make_pair(3, 3), "x^2 + 1", "poly", "α", false),
      std::invalid_argument);

  // Composite prime, overflowed prime, zero exponent, non-prime-power order.
  EXPECT_THROW(GFPX<uint8_t>(std::make_pair(4, 2)), std::invalid_argument);
  EXPECT_THROW(GFPX<uint8_t>(std::make_pair(257, 2)), std::invalid_argument);
  EXPECT_THROW(GFPX<uint8_t>(std::make_pair(3, 0)), std::invalid_argument);
  EXPECT_THROW(GFPX<uint8_t>(6), std::invalid_argument);
  EXPECT_THROW(GFPX<uint8_t>(10), std::invalid_argument);
}

TEST_F(GaloisFieldExtensionTest, StringParsingAndRepresentations) {
  auto e1 = gf9->SetElementValue("2*α + 1");
  EXPECT_EQ(e1, ExtensionElement<uint8_t>({1, 2}));

  auto e2 = gf9->SetElementValue("-α + 2");
  EXPECT_EQ(e2, ExtensionElement<uint8_t>({2, 2}));

  auto e3 = gf9->SetElementValue("2");
  EXPECT_EQ(e3, ExtensionElement<uint8_t>({2}));

  auto e4 = gf9->SetElementValue("0");
  EXPECT_TRUE(e4.IsZero());

  auto e5 = gf9->SetElementValue("α^2");
  EXPECT_LT(e5.Degree(), 2);

  auto g = gf9->MultiplicativeGenerator();
  EXPECT_EQ(gf9->SetElementValue("g^3"), gf9->Pow(g, 3));
  EXPECT_EQ(gf9->SetElementValue("g^-1"), gf9->Inv(g));
  EXPECT_EQ(gf9->SetElementValue("g^0"), gf9->MultiplicativeIdentity());

  // Hex and base-p integer parsing.
  EXPECT_EQ(gf9->SetElementValue("5"), gf9->FromInt(5));
  EXPECT_EQ(gf9->SetElementValue("0x5"), gf9->FromInt(5));

  // Invalid strings throw std::invalid_argument.
  EXPECT_THROW(gf9->SetElementValue(""), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("invalid"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("β + 1"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("h^1"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("g^"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("g^abc"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("2*"), std::invalid_argument);
  EXPECT_THROW(gf9->SetElementValue("α + "), std::invalid_argument);

  // All 5 representations work without aborting.
  EXPECT_NO_THROW(GFPX<uint8_t>(9, "", "int"));
  EXPECT_NO_THROW(GFPX<uint8_t>(9, "", "hex"));
  EXPECT_NO_THROW(GFPX<uint8_t>(9, "", "pow"));
  EXPECT_NO_THROW(GFPX<uint8_t>(9, "", "log"));
  EXPECT_NO_THROW(GFPX<uint8_t>(9, "", "poly"));

  gf9->SetRepresentation(FieldRepresentation::POLY);
  EXPECT_EQ(gf9->ToString(e1), "2*α + 1");

  gf9->SetRepresentation(FieldRepresentation::INT);
  EXPECT_EQ(gf9->ToString(e1), std::to_string(gf9->ToInt(e1)));

  gf9->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(gf9->ToString(e1), "0x7");

  gf9->SetRepresentation(FieldRepresentation::POW);
  EXPECT_EQ(gf9->ToString(gf9->Pow(g, 3)), "g^3");
  EXPECT_EQ(gf9->ToString(gf9->AdditiveIdentity()), "0");

  gf9->SetRepresentation(FieldRepresentation::LOG);
  EXPECT_EQ(gf9->ToString(gf9->Pow(g, 3)), "3");
  EXPECT_EQ(gf9->ToString(gf9->AdditiveIdentity()), "undefined");

  gf9->SetRepresentation(FieldRepresentation::POLY);

  std::ostringstream oss;
  gf9->Print(oss);
  EXPECT_FALSE(oss.str().empty());
}

TEST_F(GaloisFieldExtensionTest, UnreducedInputsAndReduction) {
  // Unreduced coefficients >= p and degree >= m.
  ExtensionElement<uint8_t> unreduced({4, 5, 1, 2});
  auto reduced = gf9->SetElementValue(unreduced);
  EXPECT_LT(reduced.Degree(), 2);
  for (uint8_t c : reduced.Coeffs()) {
    EXPECT_LT(c, 3U);
  }
  EXPECT_EQ(gf9->Add(unreduced, gf9->AdditiveIdentity()), reduced);
  EXPECT_EQ(gf9->Mul(unreduced, gf9->MultiplicativeIdentity()), reduced);
}

TEST_F(GaloisFieldExtensionTest, RandomAndMemoryBytes) {
  std::mt19937_64 rng1(42);
  std::mt19937_64 rng2(42);
  std::set<ExtensionElement<uint8_t>> distinct;
  for (int i = 0; i < 64; ++i) {
    auto a = gf9->Random(rng1);
    auto b = gf9->Random(rng2);
    EXPECT_EQ(a, b);
    EXPECT_LT(a.Degree(), 2);
    distinct.insert(a);
  }
  EXPECT_GE(distinct.size(), 5U);
  EXPECT_GT(gf9->MemoryBytes(), 0U);
}

TEST_F(GaloisFieldExtensionTest, ElementContainersAndHashing) {
  std::unordered_set<ExtensionElement<uint8_t>> hash_set;
  std::set<ExtensionElement<uint8_t>> ordered_set;
  for (uint64_t i = 0; i < gf9->Order(); ++i) {
    auto e = gf9->FromInt(i);
    hash_set.insert(e);
    ordered_set.insert(e);
  }
  EXPECT_EQ(hash_set.size(), 9U);
  EXPECT_EQ(ordered_set.size(), 9U);
}

TEST_F(GaloisFieldExtensionTest, GaloisFieldElementWrapper) {
  using Elem = GaloisFieldElement<GFPX<uint8_t>>;
  Elem a("2*α + 1", gf9);
  Elem b("α + 2", gf9);
  Elem zero("0", gf9);
  Elem one("1", gf9);

  EXPECT_EQ((a + b).Value(), gf9->Add(a.Value(), b.Value()));
  EXPECT_EQ((a - b).Value(), gf9->Sub(a.Value(), b.Value()));
  EXPECT_EQ((a * b).Value(), gf9->Mul(a.Value(), b.Value()));
  EXPECT_EQ((a / b).Value(), gf9->Div(a.Value(), b.Value()));
  EXPECT_EQ((-a).Value(), gf9->Neg(a.Value()));
  EXPECT_EQ((a * a.Inv()), one);
  EXPECT_EQ(a.PowSigned(-1), a.Inv());
  EXPECT_EQ(a.Frobenius().Value(), gf9->Frobenius(a.Value()));
  EXPECT_EQ(a.Trace().Value(), gf9->Trace(a.Value()));
  EXPECT_EQ(a.Norm().Value(), gf9->Norm(a.Value()));

  auto sq = a * a;
  auto r = sq.Sqrt();
  EXPECT_EQ(r * r, sq);

  EXPECT_THROW(a / zero, std::domain_error);
  EXPECT_THROW(zero.Inv(), std::domain_error);

  std::unordered_set<Elem> elem_set{a, b, zero, one};
  EXPECT_EQ(elem_set.size(), 4U);

  std::ostringstream oss;
  oss << a;
  EXPECT_EQ(oss.str(), "2*α + 1");
}

TEST_F(GaloisFieldExtensionTest, GFArrayOverExtensionField) {
  using Arr = GFArray<GFPX<uint8_t>>;
  using Elem = GaloisFieldElement<GFPX<uint8_t>>;

  auto zeros = Arr::Zeros(gf9, {2, 3});
  auto ones = Arr::Ones(gf9, {2, 3});
  EXPECT_EQ(zeros.Size(), 6U);
  EXPECT_EQ(ones.Size(), 6U);
  EXPECT_TRUE(zeros.At(0, 0).IsZero());
  EXPECT_TRUE(ones.At(1, 2).IsOne());

  std::mt19937_64 rng(12345);
  auto a = Arr::Random(gf9, {2, 3}, rng);
  auto b = Arr::Random(gf9, {2, 3}, rng);

  auto sum = a + b;
  auto diff = sum - b;
  EXPECT_EQ(diff, a);

  auto prod = a * b;
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(prod.At(i, j), a.At(i, j) * b.At(i, j));
    }
  }

  Elem scalar("2*α + 1", gf9);
  auto scaled_r = a * scalar;
  auto scaled_l = scalar * a;
  EXPECT_EQ(scaled_r, scaled_l);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(scaled_r.At(i, j), a.At(i, j) * scalar);
    }
  }

  auto elem_arr = a.ToElementArray();
  auto round_trip = Arr::FromElementArray(elem_arr);
  EXPECT_EQ(round_trip, a);
}
