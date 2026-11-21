#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

namespace {

// Exhaustive differential oracle comparing GF2X, GFBELogTables,
// GFBELogTablesOpt, and GFBEZechLogTables across every pair (a, b) in GF(2^m).
template <typename ElementType>
void VerifyBinaryTablesExhaustive(uint8_t m) {
  const GF2X<ElementType> ref(m);
  const GFBELogTables<ElementType> tbl(m);
  const GFBELogTablesOpt<ElementType> opt(m);
  const GFBEZechLogTables<ElementType> zech(m);

  const uint64_t order = ref.Order();
  const uint64_t g_ord = ref.GroupOrder();
  ASSERT_EQ(order, uint64_t{1} << m);
  ASSERT_EQ(tbl.Order(), order);
  ASSERT_EQ(opt.Order(), order);
  ASSERT_EQ(zech.Order(), order);

  // Verify table memory accounting.
  EXPECT_GT(tbl.MemoryBytes(), ref.MemoryBytes());
  EXPECT_GT(opt.MemoryBytes(), tbl.MemoryBytes());
  EXPECT_GT(zech.MemoryBytes(), ref.MemoryBytes());

  // Verify identities & zero/one predicates (LSP compliance for Zech logs).
  EXPECT_TRUE(zech.IsZero(zech.AdditiveIdentity()));
  EXPECT_FALSE(zech.IsOne(zech.AdditiveIdentity()));
  EXPECT_TRUE(zech.IsOne(zech.MultiplicativeIdentity()));
  EXPECT_FALSE(zech.IsZero(zech.MultiplicativeIdentity()));
  EXPECT_EQ(zech.GetElementValue(zech.AdditiveIdentity()), ElementType{0});
  EXPECT_EQ(zech.GetElementValue(zech.MultiplicativeIdentity()),
            ElementType{1});

  // Verify zero error handling across all implementations.
  EXPECT_THROW(tbl.Inv(0), std::domain_error);
  EXPECT_THROW(opt.Inv(0), std::domain_error);
  EXPECT_THROW(zech.Inv(zech.AdditiveIdentity()), std::domain_error);
  EXPECT_THROW(tbl.Log(0), std::domain_error);
  EXPECT_THROW(opt.Log(0), std::domain_error);
  EXPECT_THROW(zech.Log(zech.AdditiveIdentity()), std::domain_error);
  EXPECT_THROW(tbl.Div(1, 0), std::domain_error);
  EXPECT_THROW(opt.Div(1, 0), std::domain_error);
  EXPECT_THROW(zech.Div(zech.MultiplicativeIdentity(), zech.AdditiveIdentity()),
               std::domain_error);

  const std::vector<uint64_t> test_exps = {
      0ULL, 1ULL, 2ULL, 3ULL, g_ord, g_ord + 1ULL, 1000ULL, ~uint64_t{0}};

  const auto gens = ref.MultiplicativeGenerators();
  const ElementType alt_gen = gens.back();
  const ElementType z_alt_gen = zech.SetElementValue(alt_gen);

  for (uint64_t ai = 0; ai < order; ++ai) {
    const auto a = static_cast<ElementType>(ai);
    const auto za = zech.SetElementValue(a);
    ASSERT_EQ(zech.GetElementValue(za), a);

    // Unary ops: Neg, Sqrt, Inv, Log, Frobenius, Trace, Norm.
    EXPECT_EQ(tbl.Neg(a), ref.Neg(a));
    EXPECT_EQ(opt.Neg(a), ref.Neg(a));
    EXPECT_EQ(zech.GetElementValue(zech.Neg(za)), ref.Neg(a));

    const ElementType ref_sqrt = ref.Sqrt(a);
    EXPECT_EQ(ref.Mul(ref_sqrt, ref_sqrt), a);
    EXPECT_EQ(tbl.Sqrt(a), ref_sqrt);
    EXPECT_EQ(opt.Sqrt(a), ref_sqrt);
    EXPECT_EQ(zech.GetElementValue(zech.Sqrt(za)), ref_sqrt);

    EXPECT_EQ(tbl.Frobenius(a), ref.Frobenius(a));
    EXPECT_EQ(opt.Frobenius(a), ref.Frobenius(a));
    EXPECT_EQ(zech.GetElementValue(zech.Frobenius(za)), ref.Frobenius(a));

    EXPECT_EQ(tbl.Trace(a), ref.Trace(a));
    EXPECT_EQ(opt.Trace(a), ref.Trace(a));
    EXPECT_EQ(zech.GetElementValue(zech.Trace(za)), ref.Trace(a));

    EXPECT_EQ(tbl.Norm(a), ref.Norm(a));
    EXPECT_EQ(opt.Norm(a), ref.Norm(a));
    EXPECT_EQ(zech.GetElementValue(zech.Norm(za)), ref.Norm(a));

    if (a != 0) {
      const ElementType ref_inv = ref.Inv(a);
      EXPECT_EQ(ref.Mul(a, ref_inv), ElementType{1});
      EXPECT_EQ(tbl.Inv(a), ref_inv);
      EXPECT_EQ(opt.Inv(a), ref_inv);
      EXPECT_EQ(zech.GetElementValue(zech.Inv(za)), ref_inv);

      const uint64_t ref_log = ref.Log(a);
      EXPECT_EQ(tbl.Log(a), ref_log);
      EXPECT_EQ(opt.Log(a), ref_log);
      EXPECT_EQ(zech.Log(za), ref_log);

      const uint64_t ref_log_alt = ref.Log(a, alt_gen);
      EXPECT_EQ(tbl.Log(a, alt_gen), ref_log_alt);
      EXPECT_EQ(opt.Log(a, alt_gen), ref_log_alt);
      EXPECT_EQ(zech.Log(za, z_alt_gen), ref_log_alt);

      const uint64_t ref_ord = ref.MultiplicativeOrder(a);
      EXPECT_EQ(tbl.MultiplicativeOrder(a), ref_ord);
      EXPECT_EQ(opt.MultiplicativeOrder(a), ref_ord);
      EXPECT_EQ(zech.MultiplicativeOrder(za), ref_ord);
    }

    for (uint64_t exp : test_exps) {
      const ElementType ref_pow = ref.Pow(a, exp);
      EXPECT_EQ(tbl.Pow(a, exp), ref_pow);
      EXPECT_EQ(opt.Pow(a, exp), ref_pow);
      EXPECT_EQ(zech.GetElementValue(zech.Pow(za, exp)), ref_pow);
    }

    // Binary ops across ALL pairs (a, b): Add, Sub, Mul, Div.
    for (uint64_t bi = 0; bi < order; ++bi) {
      const auto b = static_cast<ElementType>(bi);
      const auto zb = zech.SetElementValue(b);

      const ElementType ref_add = ref.Add(a, b);
      EXPECT_EQ(tbl.Add(a, b), ref_add);
      EXPECT_EQ(opt.Add(a, b), ref_add);
      EXPECT_EQ(zech.GetElementValue(zech.Add(za, zb)), ref_add);

      const ElementType ref_sub = ref.Sub(a, b);
      EXPECT_EQ(tbl.Sub(a, b), ref_sub);
      EXPECT_EQ(opt.Sub(a, b), ref_sub);
      EXPECT_EQ(zech.GetElementValue(zech.Sub(za, zb)), ref_sub);

      const ElementType ref_mul = ref.Mul(a, b);
      EXPECT_EQ(tbl.Mul(a, b), ref_mul);
      EXPECT_EQ(opt.Mul(a, b), ref_mul);
      EXPECT_EQ(zech.GetElementValue(zech.Mul(za, zb)), ref_mul);

      if (b != 0) {
        const ElementType ref_div = ref.Div(a, b);
        EXPECT_EQ(tbl.Div(a, b), ref_div);
        EXPECT_EQ(opt.Div(a, b), ref_div);
        EXPECT_EQ(zech.GetElementValue(zech.Div(za, zb)), ref_div);
      }
    }
  }
}

template <typename ElementType>
void VerifyBinaryTablesRandomized(uint8_t m, uint64_t seed) {
  const GF2X<ElementType> ref(m);
  const GFBELogTables<ElementType> tbl(m);
  const GFBELogTablesOpt<ElementType> opt(m);
  const GFBEZechLogTables<ElementType> zech(m);

  std::mt19937_64 rng(seed);
  constexpr size_t kNumSamples = 500;
  std::vector<ElementType> a_vec(kNumSamples), b_vec(kNumSamples);
  std::vector<ElementType> za_vec(kNumSamples), zb_vec(kNumSamples);

  for (size_t i = 0; i < kNumSamples; ++i) {
    a_vec[i] = ref.Random(rng);
    b_vec[i] = ref.Random(rng);
    // Force corner cases in the first few indices.
    if (i == 0) {
      a_vec[i] = 0;
      b_vec[i] = 0;
    } else if (i == 1) {
      a_vec[i] = 1;
      b_vec[i] = 1;
    } else if (i == 2) {
      a_vec[i] = 0;
      b_vec[i] = 1;
    } else if (i == 3) {
      a_vec[i] = 1;
      b_vec[i] = 0;
    } else if (i == 4) {
      b_vec[i] = a_vec[i];
    }
    za_vec[i] = zech.SetElementValue(a_vec[i]);
    zb_vec[i] = zech.SetElementValue(b_vec[i]);
  }

  for (size_t i = 0; i < kNumSamples; ++i) {
    const ElementType a = a_vec[i];
    const ElementType b = b_vec[i];
    const ElementType za = za_vec[i];
    const ElementType zb = zb_vec[i];

    EXPECT_EQ(tbl.Add(a, b), ref.Add(a, b));
    EXPECT_EQ(opt.Add(a, b), ref.Add(a, b));
    EXPECT_EQ(zech.GetElementValue(zech.Add(za, zb)), ref.Add(a, b));

    EXPECT_EQ(tbl.Mul(a, b), ref.Mul(a, b));
    EXPECT_EQ(opt.Mul(a, b), ref.Mul(a, b));
    EXPECT_EQ(zech.GetElementValue(zech.Mul(za, zb)), ref.Mul(a, b));

    EXPECT_EQ(tbl.Sqrt(a), ref.Sqrt(a));
    EXPECT_EQ(opt.Sqrt(a), ref.Sqrt(a));
    EXPECT_EQ(zech.GetElementValue(zech.Sqrt(za)), ref.Sqrt(a));

    if (b != 0) {
      EXPECT_EQ(tbl.Div(a, b), ref.Div(a, b));
      EXPECT_EQ(opt.Div(a, b), ref.Div(a, b));
      EXPECT_EQ(zech.GetElementValue(zech.Div(za, zb)), ref.Div(a, b));
      EXPECT_EQ(tbl.Inv(b), ref.Inv(b));
      EXPECT_EQ(opt.Inv(b), ref.Inv(b));
      EXPECT_EQ(zech.GetElementValue(zech.Inv(zb)), ref.Inv(b));
      EXPECT_EQ(tbl.Log(b), ref.Log(b));
      EXPECT_EQ(opt.Log(b), ref.Log(b));
      EXPECT_EQ(zech.Log(zb), ref.Log(b));
    }
  }

  // Verify bulk span kernels agree across all 4 implementations!
  std::vector<ElementType> d_ref(kNumSamples), d_tbl(kNumSamples),
      d_opt(kNumSamples), d_zech(kNumSamples);

  ref.AddArray(d_ref, a_vec, b_vec);
  tbl.AddArray(d_tbl, a_vec, b_vec);
  opt.AddArray(d_opt, a_vec, b_vec);
  zech.AddArray(d_zech, za_vec, zb_vec);
  EXPECT_EQ(d_tbl, d_ref);
  EXPECT_EQ(d_opt, d_ref);
  for (size_t i = 0; i < kNumSamples; ++i) {
    EXPECT_EQ(zech.GetElementValue(d_zech[i]), d_ref[i]);
  }

  ref.MulArray(d_ref, a_vec, b_vec);
  tbl.MulArray(d_tbl, a_vec, b_vec);
  opt.MulArray(d_opt, a_vec, b_vec);
  zech.MulArray(d_zech, za_vec, zb_vec);
  EXPECT_EQ(d_tbl, d_ref);
  EXPECT_EQ(d_opt, d_ref);
  for (size_t i = 0; i < kNumSamples; ++i) {
    EXPECT_EQ(zech.GetElementValue(d_zech[i]), d_ref[i]);
  }

  const ElementType scalar = ref.Random(rng) | ElementType{1};
  const ElementType z_scalar = zech.SetElementValue(scalar);

  ref.ScalarMulArray(d_ref, a_vec, scalar);
  tbl.ScalarMulArray(d_tbl, a_vec, scalar);
  opt.ScalarMulArray(d_opt, a_vec, scalar);
  zech.ScalarMulArray(d_zech, za_vec, z_scalar);
  EXPECT_EQ(d_tbl, d_ref);
  EXPECT_EQ(d_opt, d_ref);
  for (size_t i = 0; i < kNumSamples; ++i) {
    EXPECT_EQ(zech.GetElementValue(d_zech[i]), d_ref[i]);
  }

  d_ref = b_vec;
  d_tbl = b_vec;
  d_opt = b_vec;
  d_zech = zb_vec;
  ref.FmaArray(d_ref, a_vec, scalar);
  tbl.FmaArray(d_tbl, a_vec, scalar);
  opt.FmaArray(d_opt, a_vec, scalar);
  zech.FmaArray(d_zech, za_vec, z_scalar);
  EXPECT_EQ(d_tbl, d_ref);
  EXPECT_EQ(d_opt, d_ref);
  for (size_t i = 0; i < kNumSamples; ++i) {
    EXPECT_EQ(zech.GetElementValue(d_zech[i]), d_ref[i]);
  }

  const ElementType dot_ref = ref.DotProduct(a_vec, b_vec);
  EXPECT_EQ(tbl.DotProduct(a_vec, b_vec), dot_ref);
  EXPECT_EQ(opt.DotProduct(a_vec, b_vec), dot_ref);
  EXPECT_EQ(zech.GetElementValue(zech.DotProduct(za_vec, zb_vec)), dot_ref);
}

template <typename ElementType>
void VerifyPrimeTableExhaustive(ElementType prime) {
  const GaloisFieldPrime<ElementType> ref(prime);
  const GaloisFieldPrimeTable<ElementType> tbl(prime);

  ASSERT_EQ(ref.Order(), static_cast<uint64_t>(prime));
  ASSERT_EQ(tbl.Order(), static_cast<uint64_t>(prime));
  EXPECT_GT(tbl.MemoryBytes(), ref.MemoryBytes());

  EXPECT_THROW(tbl.Inv(0), std::domain_error);
  EXPECT_THROW(tbl.Div(1, 0), std::domain_error);
  EXPECT_THROW(tbl.Log(0), std::domain_error);

  const auto gens = ref.MultiplicativeGenerators();
  const ElementType alt_gen = gens.back();
  const uint64_t phi = static_cast<uint64_t>(prime) - 1;
  const std::vector<uint64_t> exps = {0ULL, 1ULL,       2ULL,   5ULL,
                                      phi,  phi + 1ULL, 1000ULL};

  for (uint64_t ai = 0; ai < static_cast<uint64_t>(prime); ++ai) {
    const auto a = static_cast<ElementType>(ai);
    EXPECT_EQ(tbl.Neg(a), ref.Neg(a));
    EXPECT_EQ(tbl.IsQuadraticResidue(a), ref.IsQuadraticResidue(a));
    if (ref.IsQuadraticResidue(a)) {
      const ElementType s = tbl.Sqrt(a);
      EXPECT_EQ(s, ref.Sqrt(a));
      EXPECT_EQ(tbl.Mul(s, s), a);
    } else {
      EXPECT_THROW(tbl.Sqrt(a), std::domain_error);
    }

    if (a != 0) {
      EXPECT_EQ(tbl.Inv(a), ref.Inv(a));
      EXPECT_EQ(tbl.Mul(a, tbl.Inv(a)), ElementType{1});
      EXPECT_EQ(tbl.Log(a), ref.Log(a));
      EXPECT_EQ(tbl.Log(a, alt_gen), ref.Log(a, alt_gen));
    }

    for (uint64_t exp : exps) {
      EXPECT_EQ(tbl.Pow(a, exp), ref.Pow(a, exp));
    }

    for (uint64_t bi = 0; bi < static_cast<uint64_t>(prime); ++bi) {
      const auto b = static_cast<ElementType>(bi);
      EXPECT_EQ(tbl.Add(a, b), ref.Add(a, b));
      EXPECT_EQ(tbl.Sub(a, b), ref.Sub(a, b));
      EXPECT_EQ(tbl.Mul(a, b), ref.Mul(a, b));
      if (b != 0) {
        EXPECT_EQ(tbl.Div(a, b), ref.Div(a, b));
      }
    }
  }
}

}  // namespace

TEST(GFTablesDifferentialTest, BinaryExhaustiveM1ToM8InUint8) {
  for (uint8_t m = 1; m <= 8; ++m) {
    SCOPED_TRACE("GF(2^" + std::to_string(m) + ") in uint8_t");
    VerifyBinaryTablesExhaustive<uint8_t>(m);
  }
}

TEST(GFTablesDifferentialTest, BinaryZechCompilesForAllWordWidths) {
  // Regression test for Bug B-05: GFBEZechLogTables::Add previously hardcoded
  // uint32_t and failed to compile for uint8_t and uint16_t.
  VerifyBinaryTablesExhaustive<uint8_t>(4);
  VerifyBinaryTablesExhaustive<uint16_t>(4);
  VerifyBinaryTablesExhaustive<uint32_t>(4);
}

TEST(GFTablesDifferentialTest, BinaryRandomizedM12AndM16) {
  VerifyBinaryTablesRandomized<uint16_t>(12, 202612);
  VerifyBinaryTablesRandomized<uint16_t>(16, 202616);
}

TEST(GFTablesDifferentialTest, ElementOperatorsZechAndLogOpt) {
  for (uint8_t m : {uint8_t{1}, uint8_t{3}, uint8_t{4}, uint8_t{8}}) {
    auto ref = std::make_shared<GF2X<uint8_t>>(m);
    auto opt = std::make_shared<GFBELogTablesOpt<uint8_t>>(m);
    auto zech = std::make_shared<GFBEZechLogTables<uint8_t>>(m);

    const uint64_t order = ref->Order();
    for (uint64_t ai = 0; ai < order; ++ai) {
      const auto a = static_cast<uint8_t>(ai);
      GaloisFieldElement<GF2X<uint8_t>> ea(a, ref);
      GaloisFieldElement<GFBELogTablesOpt<uint8_t>> oa(a, opt);
      GaloisFieldElement<GFBEZechLogTables<uint8_t>> za(a, zech);

      EXPECT_EQ(oa.Value(), a);
      EXPECT_EQ(za.Value(), a);
      EXPECT_EQ(oa.IsZero(), a == 0);
      EXPECT_EQ(za.IsZero(), a == 0);
      EXPECT_EQ(oa.IsOne(), a == 1);
      EXPECT_EQ(za.IsOne(), a == 1);

      // a + a == 0, 0 * a == 0
      EXPECT_TRUE((oa + oa).IsZero());
      EXPECT_TRUE((za + za).IsZero());
      GaloisFieldElement<GFBELogTablesOpt<uint8_t>> ozero(0, opt);
      GaloisFieldElement<GFBEZechLogTables<uint8_t>> zzero(0, zech);
      EXPECT_TRUE((ozero * oa).IsZero());
      EXPECT_TRUE((zzero * za).IsZero());

      // Sqrt(a)^2 == a
      EXPECT_EQ((oa.Sqrt() * oa.Sqrt()).Value(), a);
      EXPECT_EQ((za.Sqrt() * za.Sqrt()).Value(), a);

      if (a != 0) {
        // a / a == 1, Inv(1) == 1, a * Inv(a) == 1
        EXPECT_TRUE((oa / oa).IsOne());
        EXPECT_TRUE((za / za).IsOne());
        EXPECT_EQ((oa * oa.Inv()).Value(), 1);
        EXPECT_EQ((za * za.Inv()).Value(), 1);
        EXPECT_EQ(oa.PowSigned(-2).Value(), ea.PowSigned(-2).Value());
        EXPECT_EQ(za.PowSigned(-2).Value(), ea.PowSigned(-2).Value());
      }
    }

    // Inv(1) == 1 explicitly (Bug B-02 regression).
    GaloisFieldElement<GFBELogTablesOpt<uint8_t>> o1(1, opt);
    GaloisFieldElement<GFBEZechLogTables<uint8_t>> z1(1, zech);
    EXPECT_EQ(o1.Inv().Value(), 1);
    EXPECT_EQ(z1.Inv().Value(), 1);

    // Pairwise +, -, *, / for small fields (all pairs for m <= 4, sampled for
    // m == 8).
    const uint64_t step = (m <= 4) ? 1 : 17;
    for (uint64_t ai = 0; ai < order; ai += step) {
      for (uint64_t bi = 0; bi < order; bi += step) {
        const auto a = static_cast<uint8_t>(ai);
        const auto b = static_cast<uint8_t>(bi);
        GaloisFieldElement<GF2X<uint8_t>> ea(a, ref), eb(b, ref);
        GaloisFieldElement<GFBELogTablesOpt<uint8_t>> oa(a, opt), ob(b, opt);
        GaloisFieldElement<GFBEZechLogTables<uint8_t>> za(a, zech), zb(b, zech);

        EXPECT_EQ((oa + ob).Value(), (ea + eb).Value());
        EXPECT_EQ((za + zb).Value(), (ea + eb).Value());
        EXPECT_EQ((oa - ob).Value(), (ea - eb).Value());
        EXPECT_EQ((za - zb).Value(), (ea - eb).Value());
        EXPECT_EQ((oa * ob).Value(), (ea * eb).Value());
        EXPECT_EQ((za * zb).Value(), (ea * eb).Value());
        if (b != 0) {
          EXPECT_EQ((oa / ob).Value(), (ea / eb).Value());
          EXPECT_EQ((za / zb).Value(), (ea / eb).Value());
        }
      }
    }
  }
}

TEST(GFTablesDifferentialTest, ZechStringAndPrintAndRandomCompliance) {
  auto ref = std::make_shared<GF2X<uint8_t>>(4, "int", "", "α", false, "g");
  auto zech = std::make_shared<GFBEZechLogTables<uint8_t>>(4, "int", "", "α",
                                                           false, "g");

  // Verify SetElementValue(string) matches across all supported formats without
  // double-logging!
  for (const char *raw :
       {"0", "1", "5", "0xA", "g^0", "g^1", "g^5", "g^-1", "α^3 + α + 1"}) {
    const std::string s(raw);
    GaloisFieldElement<GF2X<uint8_t>> er(s, ref);
    GaloisFieldElement<GFBEZechLogTables<uint8_t>> ez(s, zech);
    EXPECT_EQ(ez.Value(), er.Value()) << "Mismatch for input string: " << s;
  }

  // Verify Print / ToString matches across all 5 representations!
  for (FieldRepresentation rep :
       {FieldRepresentation::INT, FieldRepresentation::HEX,
        FieldRepresentation::POLY, FieldRepresentation::POW,
        FieldRepresentation::LOG}) {
    ref->SetRepresentation(rep);
    zech->SetRepresentation(rep);
    for (uint8_t a = 0; a < 16; ++a) {
      const uint8_t za = zech->SetElementValue(a);
      EXPECT_EQ(zech->ToString(za), ref->ToString(a));
    }
  }

  // Verify Zech Random() only returns valid internal values {kLogZero, 0..14}
  // and never returns OOB group_order_ (15).
  std::mt19937_64 rng(777);
  std::set<uint8_t> external_seen;
  for (int i = 0; i < 400; ++i) {
    const uint8_t z = zech->Random(rng);
    EXPECT_TRUE(z == GFBEZechLogTables<uint8_t>::kLogZero || z < 15);
    const uint8_t ext = zech->GetElementValue(z);
    EXPECT_LT(ext, 16);
    external_seen.insert(ext);
  }
  EXPECT_EQ(external_seen.size(), 16U);
}

TEST(GFTablesDifferentialTest, BinaryTableMemoryGuards) {
  EXPECT_THROW(GFBELogTables<uint32_t>(23), std::invalid_argument);
  EXPECT_THROW(GFBELogTablesOpt<uint32_t>(23), std::invalid_argument);
  EXPECT_THROW(GFBEZechLogTables<uint32_t>(23), std::invalid_argument);
}

TEST(GFTablesDifferentialTest, PrimeTableExhaustiveDifferential) {
  for (uint8_t p : {uint8_t{2}, uint8_t{3}, uint8_t{5}, uint8_t{7}, uint8_t{13},
                    uint8_t{17}, uint8_t{251}}) {
    SCOPED_TRACE("GF(" + std::to_string(p) + ")");
    VerifyPrimeTableExhaustive<uint8_t>(p);
  }
  VerifyPrimeTableExhaustive<uint16_t>(257);
}

TEST(GFTablesDifferentialTest, PrimeTableBulkKernelsAndElementOperators) {
  auto ref = std::make_shared<GaloisFieldPrime<uint16_t>>(65521);
  auto tbl = std::make_shared<GaloisFieldPrimeTable<uint16_t>>(65521);

  std::mt19937_64 rng(2026);
  constexpr size_t kN = 256;
  std::vector<uint16_t> a(kN), b(kN), d_ref(kN), d_tbl(kN);
  for (size_t i = 0; i < kN; ++i) {
    a[i] = ref->Random(rng);
    b[i] = ref->Random(rng);
  }

  ref->AddArray(d_ref, a, b);
  tbl->AddArray(d_tbl, a, b);
  EXPECT_EQ(d_tbl, d_ref);

  ref->SubArray(d_ref, a, b);
  tbl->SubArray(d_tbl, a, b);
  EXPECT_EQ(d_tbl, d_ref);

  ref->MulArray(d_ref, a, b);
  tbl->MulArray(d_tbl, a, b);
  EXPECT_EQ(d_tbl, d_ref);

  const uint16_t scalar = 12345;
  ref->ScalarMulArray(d_ref, a, scalar);
  tbl->ScalarMulArray(d_tbl, a, scalar);
  EXPECT_EQ(d_tbl, d_ref);

  d_ref = b;
  d_tbl = b;
  ref->FmaArray(d_ref, a, scalar);
  tbl->FmaArray(d_tbl, a, scalar);
  EXPECT_EQ(d_tbl, d_ref);

  EXPECT_EQ(tbl->DotProduct(a, b), ref->DotProduct(a, b));

  // Element operator sanity check with GaloisFieldPrimeTable.
  GaloisFieldElement<GaloisFieldPrimeTable<uint16_t>> ea(12345, tbl);
  GaloisFieldElement<GaloisFieldPrimeTable<uint16_t>> eb(54321, tbl);
  EXPECT_EQ(((ea + eb) - eb).Value(), ea.Value());
  EXPECT_EQ(((ea * eb) / eb).Value(), ea.Value());
  EXPECT_EQ((ea * ea.Inv()).Value(), 1);

  // Memory cap guard for prime tables (> 2^22).
  EXPECT_THROW(GaloisFieldPrimeTable<uint32_t>(4194319), std::invalid_argument);
}
