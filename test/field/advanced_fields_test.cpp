#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/field/gf_bigint.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_binary_large.hpp"
#include "xgalois/field/gf_dual_basis.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_galois_group.hpp"
#include "xgalois/field/gf_normal_basis.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/field/gf_tower.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/bigint.hpp"
#include "xgalois/utils/poly.hpp"

namespace {

using xg::AnalyzeSplittingField;
using xg::BigInt;
using xg::BigUInt;
using xg::CompositeGF256Tower;
using xg::CreateGF2_163;
using xg::CreateGF2_233;
using xg::CreateGF2_256;
using xg::CreateGF2_283;
using xg::CyclotomicCosets;
using xg::DetectOnbType;
using xg::Field;
using xg::FieldIsomorphism;
using xg::FieldTower;
using xg::GaloisFieldBinary;
using xg::GaloisFieldBinaryExtension;
using xg::GaloisFieldBinaryLarge;
using xg::GaloisFieldDualBasis;
using xg::GaloisFieldExtension;
using xg::GaloisFieldNormalBinary;
using xg::GaloisFieldPrime;
using xg::GaloisFieldPrimeBig;
using xg::GaloisFieldPrimeDyn;
using xg::GaloisGroup;
using xg::NormalToPolyConverter;
using xg::OnbType;
using xg::PolynomialDense;
using xg::SplittingFieldInfo;
using xg::TowerElement;

// ============================================================================
// 1. Static Concept Verification (xg::Field)
// ============================================================================

static_assert(Field<FieldTower<GaloisFieldPrime<uint32_t>, 2>>);
static_assert(Field<FieldTower<FieldTower<GaloisFieldPrime<uint32_t>, 2>, 3>>);
static_assert(Field<FieldTower<GaloisFieldPrime<uint32_t>>>);
static_assert(Field<FieldTower<GaloisFieldPrimeBig<256>, 2>>);
static_assert(Field<GaloisFieldNormalBinary>);
static_assert(Field<GaloisFieldPrimeBig<64>>);
static_assert(Field<GaloisFieldPrimeBig<256>>);
static_assert(Field<GaloisFieldPrimeBig<384>>);
static_assert(Field<GaloisFieldPrimeDyn>);
static_assert(Field<GaloisFieldBinaryLarge<163>>);
static_assert(Field<GaloisFieldBinaryLarge<233>>);
static_assert(Field<GaloisFieldBinaryLarge<256>>);
static_assert(Field<GaloisFieldBinaryLarge<283>>);

// Reference FIPS-197 AES S-box table (all 256 entries).
constexpr std::array<uint8_t, 256> kFips197AesSbox = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b,
    0xfe, 0xd7, 0xab, 0x76, 0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, 0xb7, 0xfd, 0x93, 0x26,
    0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2,
    0xeb, 0x27, 0xb2, 0x75, 0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, 0x53, 0xd1, 0x00, 0xed,
    0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f,
    0x50, 0x3c, 0x9f, 0xa8, 0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, 0xcd, 0x0c, 0x13, 0xec,
    0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14,
    0xde, 0x5e, 0x0b, 0xdb, 0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, 0xe7, 0xc8, 0x37, 0x6d,
    0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f,
    0x4b, 0xbd, 0x8b, 0x8a, 0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, 0xe1, 0xf8, 0x98, 0x11,
    0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f,
    0xb0, 0x54, 0xbb, 0x16};

// ============================================================================
// 2. FieldTower & CompositeGF256Tower Tests
// ============================================================================

TEST(FieldTowerTest, QuadraticExtensionGF7Squared) {
  // F_{7^2} = F_7[u] / (u^2 + 1), irreducible since -1 == 6 is QNR mod 7.
  GaloisFieldPrime<uint32_t> f7(7);
  FieldTower<GaloisFieldPrime<uint32_t>, 2> f49(f7, {1U, 0U, 1U}, "u");

  EXPECT_TRUE(f49.IsBinomial());
  EXPECT_EQ(f49.RelativeDegree(), 2U);
  EXPECT_EQ(f49.AbsoluteDegree(), 2U);
  EXPECT_EQ(f49.Characteristic(), 7ULL);
  EXPECT_EQ(f49.BigCharacteristic(), BigInt(7));
  EXPECT_EQ(f49.BigOrder(), BigInt(49));

  const auto zero = f49.AdditiveIdentity();
  const auto one = f49.MultiplicativeIdentity();
  const auto u = f49.Generator();
  EXPECT_TRUE(f49.IsZero(zero));
  EXPECT_TRUE(f49.IsOne(one));
  EXPECT_EQ(f49.Mul(u, u), f49.EmbedFromBase(6U));  // u^2 == -1 == 6 mod 7.

  // EmbedFromBase / ProjectToBase round-trip.
  for (uint32_t c = 0; c < 7; ++c) {
    const auto emb = f49.EmbedFromBase(c);
    const auto proj = f49.ProjectToBase(emb);
    ASSERT_TRUE(proj.has_value());
    EXPECT_EQ(*proj, c);
  }
  EXPECT_FALSE(f49.ProjectToBase(u).has_value());

  // Exhaustively verify multiplicative inverses, Frobenius, Trace, and Norm on
  // all 49 elements of F_{7^2}.
  for (uint32_t a0 = 0; a0 < 7; ++a0) {
    for (uint32_t a1 = 0; a1 < 7; ++a1) {
      const TowerElement<uint32_t, 2> a{a0, a1};
      EXPECT_EQ(f49.Add(a, f49.Neg(a)), zero);

      // Relative Frobenius: (a0 + a1*u)^7 == a0 - a1*u.
      const auto frob = f49.RelativeFrobenius(a);
      EXPECT_EQ(frob, f49.Pow(a, 7));
      EXPECT_EQ(frob[0], a0);
      EXPECT_EQ(frob[1], f7.Neg(a1));
      EXPECT_EQ(f49.RelativeFrobenius(a, 2), a);

      // RelativeTrace = 2*a0 mod 7; RelativeNorm = a0^2 + a1^2 mod 7.
      EXPECT_EQ(f49.RelativeTrace(a), f7.Add(a0, a0));
      EXPECT_EQ(f49.RelativeNorm(a), f7.Add(f7.Mul(a0, a0), f7.Mul(a1, a1)));

      if (!f49.IsZero(a)) {
        const auto inv_a = f49.Inv(a);
        EXPECT_TRUE(f49.IsOne(f49.Mul(a, inv_a)));
        EXPECT_EQ(f49.Pow(a, 48), one);
        EXPECT_EQ(f49.PowBig(a, BigInt(-1)), inv_a);
      }
    }
  }
  EXPECT_THROW((void)f49.Inv(zero), std::domain_error);
}

TEST(FieldTowerTest, RecursiveTowerGF7ToGF7SquaredToGF7Sixth) {
  // F_7 -> K = F_{7^2} = F_7[u]/(u^2 + 1) -> L = F_{7^6} = K[v]/(v^3 - (1+u)).
  // Since (1 + u)^((49 - 1)/3) = (1 + u)^16 = 4 != 1 in K, v^3 - (1 + u) is
  // irreducible over K.
  using Fp = GaloisFieldPrime<uint32_t>;
  using Fp2 = FieldTower<Fp, 2>;
  using Fp6 = FieldTower<Fp2, 3>;

  Fp f7(7);
  Fp2 f49(f7, {1U, 0U, 1U}, "u");
  const TowerElement<uint32_t, 2> xi{1U, 1U};  // 1 + u
  EXPECT_FALSE(f49.IsOne(f49.Pow(xi, 16)));

  Fp6 f7_6 = Fp6::Binomial(f49, 3, xi, "v");

  EXPECT_TRUE(f7_6.IsBinomial());
  EXPECT_EQ(f7_6.RelativeDegree(), 3U);
  EXPECT_EQ(f7_6.AbsoluteDegree(), 6U);
  EXPECT_EQ(f7_6.Characteristic(), 7ULL);
  EXPECT_EQ(f7_6.BigOrder(), BigInt(7).Pow(6));

  std::mt19937_64 rng(20261122ULL);
  std::uniform_int_distribution<uint32_t> dist(0, 6);
  auto rand_fp2 = [&]() {
    return TowerElement<uint32_t, 2>{dist(rng), dist(rng)};
  };
  auto rand_fp6 = [&]() {
    return TowerElement<TowerElement<uint32_t, 2>, 3>{rand_fp2(), rand_fp2(),
                                                      rand_fp2()};
  };

  for (int iter = 0; iter < 40; ++iter) {
    const auto a = rand_fp6();
    const auto b = rand_fp6();
    const auto c = rand_fp6();

    // Distributivity & associativity.
    EXPECT_EQ(f7_6.Mul(a, f7_6.Add(b, c)),
              f7_6.Add(f7_6.Mul(a, b), f7_6.Mul(a, c)));
    EXPECT_EQ(f7_6.Mul(f7_6.Mul(a, b), c), f7_6.Mul(a, f7_6.Mul(b, c)));

    if (!f7_6.IsZero(a)) {
      const auto inv_a = f7_6.Inv(a);
      EXPECT_TRUE(f7_6.IsOne(f7_6.Mul(a, inv_a)));
      EXPECT_EQ(f7_6.Div(f7_6.Mul(b, a), a), b);
    }

    // Relative Frobenius: a -> a^(7^2) = a^49.
    EXPECT_EQ(f7_6.RelativeFrobenius(a, 1), f7_6.Pow(a, 49));
    EXPECT_EQ(f7_6.RelativeFrobenius(a, 2), f7_6.Pow(a, 49ULL * 49ULL));
    EXPECT_EQ(f7_6.RelativeFrobenius(a, 3), a);

    // Transitivity of Trace & Norm:
    //   Tr_{L/F}(a) = Tr_{K/F}(Tr_{L/K}(a)) = sum_{i=0}^{5} a^(7^i)
    //   N_{L/F}(a)  = N_{K/F}(N_{L/K}(a))   = prod_{i=0}^{5} a^(7^i).
    auto abs_trace_elem = f7_6.AdditiveIdentity();
    auto abs_norm_elem = f7_6.MultiplicativeIdentity();
    auto cur = a;
    for (int i = 0; i < 6; ++i) {
      abs_trace_elem = f7_6.Add(abs_trace_elem, cur);
      abs_norm_elem = f7_6.Mul(abs_norm_elem, cur);
      cur = f7_6.Pow(cur, 7);
    }

    const auto tr_l_to_k = f7_6.RelativeTrace(a);
    const uint32_t tr_trans = f49.RelativeTrace(tr_l_to_k);
    const auto proj_tr_k = f7_6.ProjectToBase(abs_trace_elem);
    ASSERT_TRUE(proj_tr_k.has_value());
    const auto proj_tr_f = f49.ProjectToBase(*proj_tr_k);
    ASSERT_TRUE(proj_tr_f.has_value());
    EXPECT_EQ(tr_trans, *proj_tr_f);

    const auto nm_l_to_k = f7_6.RelativeNorm(a);
    const uint32_t nm_trans = f49.RelativeNorm(nm_l_to_k);
    const auto proj_nm_k = f7_6.ProjectToBase(abs_norm_elem);
    ASSERT_TRUE(proj_nm_k.has_value());
    const auto proj_nm_f = f49.ProjectToBase(*proj_nm_k);
    ASSERT_TRUE(proj_nm_f.has_value());
    EXPECT_EQ(nm_trans, *proj_nm_f);
  }
}

TEST(FieldTowerTest, NonBinomialDynamicDegreeAndBigPrimeTower) {
  // General non-binomial irreducible: x^3 + x + 1 over GF(5) (no roots in F_5:
  // 0->1, 1->3, 2->1, 3->1, 4->4).
  GaloisFieldPrime<uint32_t> f5(5);
  FieldTower<GaloisFieldPrime<uint32_t>> f125(f5, {1U, 1U, 0U, 1U}, "x");
  EXPECT_FALSE(f125.IsBinomial());
  EXPECT_EQ(f125.RelativeDegree(), 3U);
  EXPECT_EQ(f125.BigOrder(), BigInt(125));

  const TowerElement<uint32_t> a{2U, 3U, 4U};
  const TowerElement<uint32_t> b{1U, 4U, 2U};
  EXPECT_TRUE(f125.IsOne(f125.Mul(a, f125.Inv(a))));
  EXPECT_EQ(f125.Div(f125.Mul(a, b), b), a);
  EXPECT_EQ(f125.RelativeFrobenius(a), f125.Pow(a, 5));
  EXPECT_EQ(f125.Pow(a, 124), f125.MultiplicativeIdentity());
  EXPECT_FALSE(f125.ToString(a).empty());

  // BN254 F_{p^2} = F_p[u] / (u^2 + 1) over a 254-bit prime!
  auto bn254 = xg::primes::Bn254();
  const auto one254 = bn254.MultiplicativeIdentity();
  const auto zero254 = bn254.AdditiveIdentity();
  FieldTower<GaloisFieldPrimeBig<256>, 2> bn254_fp2(
      bn254, {one254, zero254, one254}, "u");
  EXPECT_TRUE(bn254_fp2.IsBinomial());
  const TowerElement<BigUInt<256>, 2> z{BigUInt<256>(123456789ULL),
                                        BigUInt<256>(987654321ULL)};
  const auto z_inv = bn254_fp2.Inv(z);
  EXPECT_TRUE(bn254_fp2.IsOne(bn254_fp2.Mul(z, z_inv)));
  // Frobenius in F_{p^2} over F_p(i) is complex conjugation: (a + bi)^p = a -
  // bi.
  const auto z_conj = bn254_fp2.RelativeFrobenius(z);
  EXPECT_EQ(z_conj[0], z[0]);
  EXPECT_EQ(z_conj[1], bn254.Neg(z[1]));
}

TEST(CompositeGF256TowerTest, MatchesAll256Fips197AesSboxAndArithmetic) {
  GaloisFieldBinaryExtension<uint8_t> aes_gf(8, "int",
                                             "x^8 + x^4 + x^3 + x + 1");

  // Verify every single FIPS-197 AES S-box and InvSbox entry (all 256 bytes).
  for (uint32_t i = 0; i < 256; ++i) {
    const uint8_t x = static_cast<uint8_t>(i);
    const uint8_t s = CompositeGF256Tower::AesSbox(x);
    EXPECT_EQ(s, kFips197AesSbox[i]) << "AES S-box mismatch at x=" << i;
    EXPECT_EQ(CompositeGF256Tower::AesInvSbox(s), x)
        << "AES InvSbox mismatch at x=" << i;

    // Verify ToTower / FromTower round-trip and Inv.
    const uint8_t t = CompositeGF256Tower::ToTower(x);
    EXPECT_EQ(CompositeGF256Tower::FromTower(t), x);
    if (x == 0) {
      EXPECT_EQ(CompositeGF256Tower::Inv(0), 0);
    } else {
      EXPECT_EQ(CompositeGF256Tower::Inv(x), aes_gf.Inv(x));
    }
  }

  // Verify multiplication across all 256 * 256 = 65536 pairs.
  for (uint32_t a = 0; a < 256; ++a) {
    for (uint32_t b = 0; b < 256; ++b) {
      const uint8_t ua = static_cast<uint8_t>(a);
      const uint8_t ub = static_cast<uint8_t>(b);
      EXPECT_EQ(CompositeGF256Tower::Mul(ua, ub), aes_gf.Mul(ua, ub));
    }
  }
}

// ============================================================================
// 3. GaloisGroup, CyclotomicCosets, FieldIsomorphism, AnalyzeSplittingField
// ============================================================================

TEST(GaloisGroupTest, FrobeniusOrbitsAndMinimalPolynomials) {
  GaloisFieldBinaryExtension<uint8_t> gf256(8, "int",
                                            "x^8 + x^4 + x^3 + x + 1");
  GaloisGroup<GaloisFieldBinaryExtension<uint8_t>> gal(gf256);
  EXPECT_EQ(gal.Order(), 8U);

  // Orbit of x (= 2) has size 8; its minimal polynomial is x^8+x^4+x^3+x+1.
  const auto orbit_x = gal.FrobeniusOrbit(uint8_t{2});
  EXPECT_EQ(orbit_x.size(), 8U);
  EXPECT_EQ(gal.Apply(uint8_t{2}, 8), uint8_t{2});

  const auto minpoly_x = gal.MinimalPolynomialOverPrime(uint8_t{2});
  const std::vector<uint64_t> expected_aes = {1, 1, 0, 1, 1, 0, 0, 0, 1};
  EXPECT_EQ(minpoly_x, expected_aes);

  // Relative Galois group Gal(GF(2^8) / GF(2^2)) has order 4.
  GaloisGroup<GaloisFieldBinaryExtension<uint8_t>> rel_gal(gf256, 2);
  EXPECT_EQ(rel_gal.Order(), 4U);
  EXPECT_EQ(rel_gal.Apply(uint8_t{0x57}, 1), gf256.Pow(uint8_t{0x57}, 4));
  EXPECT_EQ(rel_gal.Apply(uint8_t{0x57}, 4), uint8_t{0x57});

  // Odd-characteristic GaloisFieldExtension GF(3^4).
  GaloisFieldExtension<uint32_t> gf81({3, 4});
  GaloisGroup<GaloisFieldExtension<uint32_t>> gal81(gf81);
  EXPECT_EQ(gal81.Order(), 4U);
  const auto gen = gf81.MultiplicativeGenerator();
  const auto minp = gal81.MinimalPolynomialOverPrime(gen);
  EXPECT_EQ(minp.size(), 5U);
  EXPECT_EQ(minp.back(), 1ULL);
  // Evaluate minp(gen) == 0 in GF(3^4).
  auto eval = gf81.AdditiveIdentity();
  for (std::size_t i = minp.size(); i-- > 0;) {
    eval = gf81.Add(gf81.Mul(eval, gen),
                    gf81.FromCoeffs({static_cast<uint32_t>(minp[i])}));
  }
  EXPECT_TRUE(gf81.IsZero(eval));
}

TEST(GaloisGroupTest, CyclotomicCosets) {
  const auto cosets = CyclotomicCosets(15, 2);
  const std::vector<std::vector<uint64_t>> expected = {
      {0}, {1, 2, 4, 8}, {3, 6, 12, 9}, {5, 10}, {7, 14, 13, 11}};
  EXPECT_EQ(cosets, expected);

  EXPECT_THROW((void)CyclotomicCosets(14, 2), std::invalid_argument);
  EXPECT_THROW((void)CyclotomicCosets(0, 2), std::invalid_argument);
}

TEST(FieldIsomorphismTest, BinaryGF256AndGF2_16AndOddExtension) {
  // GF(2^8, 0x11D) <-> GF(2^8, 0x11B)
  GaloisFieldBinaryExtension<uint8_t> gf_11d(8, "int",
                                             "x^8 + x^4 + x^3 + x^2 + 1");
  GaloisFieldBinaryExtension<uint8_t> gf_11b(8, "int",
                                             "x^8 + x^4 + x^3 + x + 1");

  FieldIsomorphism iso(gf_11d, gf_11b);
  EXPECT_EQ(iso.Forward(0), 0);
  EXPECT_EQ(iso.Forward(1), 1);

  for (uint32_t i = 0; i < 256; ++i) {
    const uint8_t a = static_cast<uint8_t>(i);
    const uint8_t fa = iso.Forward(a);
    EXPECT_EQ(iso.Inverse(fa), a);
    EXPECT_EQ(iso.Forward(iso.Inverse(a)), a);
    if (a != 0) {
      EXPECT_EQ(iso.Forward(gf_11d.Inv(a)), gf_11b.Inv(fa));
    }
  }
  for (uint32_t i = 0; i < 256; i += 5) {
    for (uint32_t j = 0; j < 256; j += 7) {
      const uint8_t a = static_cast<uint8_t>(i);
      const uint8_t b = static_cast<uint8_t>(j);
      EXPECT_EQ(iso.Forward(gf_11d.Add(a, b)),
                gf_11b.Add(iso.Forward(a), iso.Forward(b)));
      EXPECT_EQ(iso.Forward(gf_11d.Mul(a, b)),
                gf_11b.Mul(iso.Forward(a), iso.Forward(b)));
    }
  }

  // Larger binary field GF(2^16) exercising Cantor-Zassenhaus root finding.
  GaloisFieldBinaryExtension<uint16_t> gf16_a(16, "int",
                                              "x^16 + x^5 + x^3 + x^2 + 1");
  GaloisFieldBinaryExtension<uint16_t> gf16_b(16, "int",
                                              "x^16 + x^12 + x^3 + x + 1");
  FieldIsomorphism iso16(gf16_a, gf16_b);
  std::mt19937_64 rng(42);
  for (int t = 0; t < 50; ++t) {
    const uint16_t a = gf16_a.Random(rng);
    const uint16_t b = gf16_a.Random(rng);
    EXPECT_EQ(iso16.Inverse(iso16.Forward(a)), a);
    EXPECT_EQ(iso16.Forward(gf16_a.Mul(a, b)),
              gf16_b.Mul(iso16.Forward(a), iso16.Forward(b)));
  }

  // Odd-characteristic GF(5^3) with two distinct irreducibles.
  GaloisFieldExtension<uint32_t> gf125_a({5, 3}, "x^3 + x + 1", "poly", "x",
                                         /*check_irreducible=*/true);
  GaloisFieldExtension<uint32_t> gf125_b({5, 3}, "x^3 + 3*x + 2", "poly", "x",
                                         /*check_irreducible=*/true);
  FieldIsomorphism iso125(gf125_a, gf125_b);
  for (int t = 0; t < 30; ++t) {
    const auto a = gf125_a.Random(rng);
    const auto b = gf125_a.Random(rng);
    EXPECT_EQ(iso125.Inverse(iso125.Forward(a)), a);
    EXPECT_EQ(iso125.Forward(gf125_a.Add(a, b)),
              gf125_b.Add(iso125.Forward(a), iso125.Forward(b)));
    EXPECT_EQ(iso125.Forward(gf125_a.Mul(a, b)),
              gf125_b.Mul(iso125.Forward(a), iso125.Forward(b)));
  }
}

TEST(SplittingFieldTest, AnalyzeSplittingFieldDegreesAndLcm) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  // f(x) = (x^2 + x + 1) * (x^3 + x + 1) over GF(2): irreducible degrees {2,3},
  // splitting field extension degree = lcm(2, 3) = 6.
  auto p2 = xg::utils::ParsePolynomial(gf2, "x^2 + x + 1");
  auto p3 = xg::utils::ParsePolynomial(gf2, "x^3 + x + 1");
  auto f = p2 * p3;

  const SplittingFieldInfo info = AnalyzeSplittingField(f);
  EXPECT_EQ(info.base_order, 2ULL);
  EXPECT_EQ(info.extension_degree, 6U);
  EXPECT_EQ(info.factor_degrees, (std::vector<uint32_t>{2, 3}));

  // With repeated factor: (x + 1)^2 * (x^2 + x + 1) * (x^4 + x + 1):
  // degrees {1, 1, 2, 4}, extension_degree = lcm(1, 2, 4) = 4.
  auto p1 = xg::utils::ParsePolynomial(gf2, "x + 1");
  auto p4 = xg::utils::ParsePolynomial(gf2, "x^4 + x + 1");
  auto g = p1 * p1 * p2 * p4;
  const SplittingFieldInfo info_g = AnalyzeSplittingField(g);
  EXPECT_EQ(info_g.base_order, 2ULL);
  EXPECT_EQ(info_g.extension_degree, 4U);
  EXPECT_EQ(info_g.factor_degrees, (std::vector<uint32_t>{1, 1, 2, 4}));
}

// ============================================================================
// 4. GaloisFieldNormalBinary & NormalToPolyConverter Tests
// ============================================================================

TEST(NormalBasisTest, DetectOnbType) {
  EXPECT_EQ(DetectOnbType(1), OnbType::kNone);
  EXPECT_EQ(DetectOnbType(2), OnbType::kType1);
  EXPECT_EQ(DetectOnbType(3), OnbType::kType2);
  EXPECT_EQ(DetectOnbType(4), OnbType::kType1);
  EXPECT_EQ(DetectOnbType(5), OnbType::kType2);
  EXPECT_EQ(DetectOnbType(6), OnbType::kType2);
  EXPECT_EQ(DetectOnbType(7), OnbType::kNone);
  EXPECT_EQ(DetectOnbType(8), OnbType::kNone);
  EXPECT_EQ(DetectOnbType(9), OnbType::kType2);
  EXPECT_EQ(DetectOnbType(10), OnbType::kType1);
  EXPECT_EQ(DetectOnbType(11), OnbType::kType2);
  EXPECT_EQ(DetectOnbType(12), OnbType::kType1);
  EXPECT_EQ(DetectOnbType(23), OnbType::kType2);
  EXPECT_EQ(DetectOnbType(28), OnbType::kType1);
  EXPECT_THROW((void)GaloisFieldNormalBinary(8), std::invalid_argument);
}

TEST(NormalBasisTest, FieldAxiomsAndIsomorphismType1AndType2) {
  // Test Type I (m = 4, 10) and Type II (m = 5, 6, 9, 11) — note m = 11 has
  // p = 23 == 3 (mod 4) with ord_23(2) = 11 (p | 2^m - 1), while m = 5, 6, 9
  // have ord_p(2) = 2m (p | 2^m + 1).
  const std::vector<uint32_t> degrees = {4, 5, 6, 9, 10, 11};
  std::mt19937_64 rng(1337);

  for (uint32_t m : degrees) {
    GaloisFieldNormalBinary nf(m);
    GaloisFieldBinaryExtension<uint32_t> pf(static_cast<uint8_t>(m));
    NormalToPolyConverter conv(nf, pf);

    const uint64_t mask = (uint64_t{1} << m) - 1ULL;
    EXPECT_EQ(nf.MultiplicativeIdentity(), mask);
    EXPECT_EQ(nf.AdditiveIdentity(), 0ULL);
    EXPECT_EQ(conv.ToPoly(0ULL), 0U);
    EXPECT_EQ(conv.ToPoly(mask), 1U);
    EXPECT_EQ(conv.ToNormal(1U), mask);

    // Count 1s in lambda table == 2*m - 1.
    uint32_t ones = 0;
    for (uint32_t i = 0; i < m; ++i) {
      for (uint32_t j = 0; j < m; ++j) {
        if (nf.LambdaEntry(i, j)) {
          ++ones;
        }
      }
    }
    EXPECT_EQ(ones, 2U * m - 1U) << "ONB complexity mismatch for m=" << m;

    std::uniform_int_distribution<uint64_t> dist(0, mask);
    for (int t = 0; t < 80; ++t) {
      const uint64_t a = dist(rng);
      const uint64_t b = dist(rng);
      const uint64_t c = dist(rng);

      // O(1) Frobenius, Sqrt, Trace.
      EXPECT_EQ(nf.Frobenius(a), nf.Mul(a, a));
      EXPECT_EQ(nf.Frobenius(nf.Sqrt(a)), a);
      EXPECT_EQ(nf.Sqrt(nf.Frobenius(a)), a);
      EXPECT_EQ(nf.Trace(nf.Add(a, b)), nf.Trace(a) ^ nf.Trace(b));
      EXPECT_EQ(nf.Trace(nf.Frobenius(a)), nf.Trace(a));

      // Field axioms.
      EXPECT_EQ(nf.Mul(a, b), nf.Mul(b, a));
      EXPECT_EQ(nf.Mul(nf.Mul(a, b), c), nf.Mul(a, nf.Mul(b, c)));
      EXPECT_EQ(nf.Mul(a, nf.Add(b, c)), nf.Add(nf.Mul(a, b), nf.Mul(a, c)));
      EXPECT_EQ(nf.Mul(a, mask), a);

      if (a != 0) {
        const uint64_t inv_a = nf.Inv(a);
        EXPECT_EQ(nf.Mul(a, inv_a), mask);
        EXPECT_EQ(nf.Div(nf.Mul(b, a), a), b);
        EXPECT_EQ(nf.Pow(a, (uint64_t{1} << m) - 1ULL), mask);
      }

      // Exact ring isomorphism via NormalToPolyConverter.
      const uint32_t pa = conv.ToPoly(a);
      const uint32_t pb = conv.ToPoly(b);
      EXPECT_EQ(conv.ToNormal(pa), a);
      EXPECT_EQ(conv.ToPoly(nf.Add(a, b)), pf.Add(pa, pb));
      EXPECT_EQ(conv.ToPoly(nf.Mul(a, b)), pf.Mul(pa, pb));
      EXPECT_EQ(conv.ToPoly(nf.Frobenius(a)), pf.Mul(pa, pa));
    }
  }
}

// ============================================================================
// 5. GaloisFieldDualBasis Tests
// ============================================================================

TEST(DualBasisTest, OrthogonalityAndBerlekampMultiplication) {
  GaloisFieldBinaryExtension<uint8_t> gf256(8, "int",
                                            "x^8 + x^4 + x^3 + x + 1");

  // Standard dual basis (gamma = 1) and weakly-dual basis (gamma = 0x57).
  for (uint8_t gamma : {uint8_t{1}, uint8_t{0x57}}) {
    GaloisFieldDualBasis dual(gf256, gamma);
    EXPECT_TRUE(dual.VerifyDualOrthogonality());

    for (uint32_t a = 0; a < 256; ++a) {
      const uint8_t ua = static_cast<uint8_t>(a);
      const uint8_t da = dual.ToDual(ua);
      EXPECT_EQ(dual.FromDual(da), ua);
    }

    // Verify MulPrimalDual across ALL 256 * 256 pairs!
    for (uint32_t a = 0; a < 256; ++a) {
      const uint8_t ua = static_cast<uint8_t>(a);
      for (uint32_t b = 0; b < 256; ++b) {
        const uint8_t ub = static_cast<uint8_t>(b);
        const uint8_t dual_prod = dual.MulPrimalDual(ua, dual.ToDual(ub));
        EXPECT_EQ(dual.FromDual(dual_prod), gf256.Mul(ua, ub));
      }
    }
  }

  // Also test GF(2^16).
  GaloisFieldBinaryExtension<uint16_t> gf16(16, "int",
                                            "x^16 + x^5 + x^3 + x^2 + 1");
  GaloisFieldDualBasis dual16(gf16, uint16_t{0x1234});
  EXPECT_TRUE(dual16.VerifyDualOrthogonality());
  std::mt19937_64 rng(99);
  for (int i = 0; i < 100; ++i) {
    const uint16_t a = gf16.Random(rng);
    const uint16_t b = gf16.Random(rng);
    EXPECT_EQ(dual16.FromDual(dual16.MulPrimalDual(a, dual16.ToDual(b))),
              gf16.Mul(a, b));
  }
}

// ============================================================================
// 6. GaloisFieldPrimeBig & GaloisFieldPrimeDyn Tests
// ============================================================================

template <std::size_t Bits>
void VerifyBigPrimeFieldAxiomsAndSqrt(const GaloisFieldPrimeBig<Bits>& field,
                                      uint64_t seed) {
  using E = BigUInt<Bits>;
  EXPECT_TRUE(xg::IsProbablePrime(field.BigCharacteristic()));
  EXPECT_EQ(field.Degree(), 1U);

  std::mt19937_64 rng(seed);
  auto rand_elem = [&]() -> E {
    std::array<uint64_t, E::kLimbs> limbs{};
    for (uint64_t& w : limbs) {
      w = rng();
    }
    return field.Reduce(E(limbs));
  };

  for (int t = 0; t < 20; ++t) {
    const E a = rand_elem();
    const E b = rand_elem();
    const E c = rand_elem();

    EXPECT_TRUE(field.IsZero(field.Add(a, field.Neg(a))));
    EXPECT_EQ(field.Add(a, b), field.Add(b, a));
    EXPECT_EQ(field.Mul(a, b), field.Mul(b, a));
    EXPECT_EQ(field.Sqr(a), field.Mul(a, a));
    EXPECT_EQ(field.Mul(a, field.Add(b, c)),
              field.Add(field.Mul(a, b), field.Mul(a, c)));
    EXPECT_EQ(field.Mul(field.Mul(a, b), c), field.Mul(a, field.Mul(b, c)));

    // Montgomery domain round-trip.
    EXPECT_EQ(field.FromMont(field.MulMont(field.ToMont(a), field.ToMont(b))),
              field.Mul(a, b));

    if (!field.IsZero(a)) {
      const E inv_a = field.Inv(a);
      EXPECT_TRUE(field.IsOne(field.Mul(a, inv_a)));
      EXPECT_EQ(field.Div(field.Mul(b, a), a), b);
      EXPECT_EQ(field.PowBig(a, field.BigOrder() - BigInt(1)),
                field.MultiplicativeIdentity());
    }

    // Tonelli-Shanks square root on guaranteed square a^2.
    const E sq = field.Sqr(a);
    EXPECT_TRUE(field.HasSqrt(sq));
    const auto root = field.Sqrt(sq);
    ASSERT_TRUE(root.has_value());
    EXPECT_EQ(field.Sqr(*root), sq);

    // Check HasSqrt consistency on random element b.
    const auto root_b = field.Sqrt(b);
    EXPECT_EQ(root_b.has_value(), field.HasSqrt(b));
    if (root_b.has_value()) {
      EXPECT_EQ(field.Sqr(*root_b), b);
    }
  }
}

TEST(BigPrimeFieldTest, StandardCryptoPrimesAndTonelliShanks) {
  const auto c25519 = xg::primes::Curve25519();
  EXPECT_EQ(c25519.BigCharacteristic(), (BigInt(1) << 255) - BigInt(19));

  const auto secp = xg::primes::Secp256k1();
  EXPECT_EQ(secp.BigCharacteristic(),
            (BigInt(1) << 256) - (BigInt(1) << 32) - BigInt(977));

  const auto p256 = xg::primes::NistP256();
  EXPECT_EQ(p256.BigCharacteristic(), (BigInt(1) << 256) - (BigInt(1) << 224) +
                                          (BigInt(1) << 192) +
                                          (BigInt(1) << 96) - BigInt(1));

  const auto p384 = xg::primes::NistP384();
  EXPECT_EQ(p384.BigCharacteristic(), (BigInt(1) << 384) - (BigInt(1) << 128) -
                                          (BigInt(1) << 96) +
                                          (BigInt(1) << 32) - BigInt(1));

  const auto bn254 = xg::primes::Bn254();
  const auto bls381 = xg::primes::Bls12381();
  const auto goldilocks = xg::primes::Goldilocks();
  EXPECT_EQ(goldilocks.BigCharacteristic(),
            (BigInt(1) << 64) - (BigInt(1) << 32) + BigInt(1));

  VerifyBigPrimeFieldAxiomsAndSqrt(c25519, 101);
  VerifyBigPrimeFieldAxiomsAndSqrt(secp, 102);
  VerifyBigPrimeFieldAxiomsAndSqrt(p256, 103);
  VerifyBigPrimeFieldAxiomsAndSqrt(p384, 104);
  VerifyBigPrimeFieldAxiomsAndSqrt(bn254, 105);
  VerifyBigPrimeFieldAxiomsAndSqrt(bls381, 106);
  VerifyBigPrimeFieldAxiomsAndSqrt(goldilocks, 107);

  // Confirm quadratic non-residue returns std::nullopt on Curve25519 (2 is QNR
  // mod 2^255 - 19 since p == 5 mod 8).
  EXPECT_FALSE(c25519.HasSqrt(BigUInt<256>(2ULL)));
  EXPECT_FALSE(c25519.Sqrt(BigUInt<256>(2ULL)).has_value());
}

TEST(BigPrimeFieldTest, DynamicPrecisionGaloisFieldPrimeDyn) {
  // 255-bit Curve25519 prime and a 521-bit Mersenne prime 2^521 - 1.
  GaloisFieldPrimeDyn f_dyn((BigInt(1) << 255) - BigInt(19));
  GaloisFieldPrimeDyn f_m521((BigInt(1) << 521) - BigInt(1),
                             /*check_prime=*/false);

  const BigInt a = BigInt::FromHex("0x123456789abcdef0123456789abcdef012345");
  const BigInt b = BigInt::FromHex("0xfedcba9876543210fedcba9876543210fedcb");

  EXPECT_TRUE(f_dyn.IsOne(f_dyn.Mul(a, f_dyn.Inv(a))));
  EXPECT_EQ(f_dyn.Div(f_dyn.Mul(a, b), b), f_dyn.Reduce(a));
  EXPECT_EQ(f_dyn.PowBig(a, BigInt(-2)), f_dyn.Inv(f_dyn.Mul(a, a)));

  const BigInt a_sq = f_dyn.Mul(a, a);
  const auto r_dyn = f_dyn.Sqrt(a_sq);
  ASSERT_TRUE(r_dyn.has_value());
  EXPECT_EQ(f_dyn.Mul(*r_dyn, *r_dyn), a_sq);
  EXPECT_FALSE(f_dyn.Sqrt(BigInt(2)).has_value());

  // 521-bit Mersenne prime field axioms & square root.
  const BigInt x = (BigInt(1) << 400) + BigInt(123456789);
  const BigInt y = (BigInt(1) << 300) + BigInt(987654321);
  EXPECT_TRUE(f_m521.IsOne(f_m521.Mul(x, f_m521.Inv(x))));
  EXPECT_EQ(f_m521.Div(f_m521.Mul(x, y), y), x);
  const BigInt x_sq = f_m521.Mul(x, x);
  const auto r_521 = f_m521.Sqrt(x_sq);
  ASSERT_TRUE(r_521.has_value());
  EXPECT_EQ(f_m521.Mul(*r_521, *r_521), x_sq);
}

// ============================================================================
// 7. GaloisFieldBinaryLarge Tests (GF(2^163), GF(2^233), GF(2^256), GF(2^283))
// ============================================================================

template <std::size_t M>
void VerifyLargeBinaryField(const GaloisFieldBinaryLarge<M>& field,
                            uint64_t seed) {
  using E = typename GaloisFieldBinaryLarge<M>::element_type;
  EXPECT_EQ(field.Characteristic(), 2ULL);
  EXPECT_EQ(field.Degree(), static_cast<uint32_t>(M));
  EXPECT_EQ(field.BigOrder(), BigInt(1) << M);

  std::mt19937_64 rng(seed);
  auto rand_elem = [&]() -> E {
    std::array<uint64_t, E::kLimbs> limbs{};
    for (uint64_t& w : limbs) {
      w = rng();
    }
    return field.Reduce(E(limbs));
  };

  for (int t = 0; t < 20; ++t) {
    const E a = rand_elem();
    const E b = rand_elem();
    const E c = rand_elem();

    // Additive & multiplicative axioms.
    EXPECT_TRUE(field.IsZero(field.Add(a, a)));
    EXPECT_EQ(field.Mul(a, b), field.Mul(b, a));
    EXPECT_EQ(field.Mul(field.Mul(a, b), c), field.Mul(a, field.Mul(b, c)));
    EXPECT_EQ(field.Mul(a, field.Add(b, c)),
              field.Add(field.Mul(a, b), field.Mul(a, c)));
    EXPECT_EQ(field.Sqr(a), field.Mul(a, a));

    // Inv & Div.
    if (!field.IsZero(a)) {
      const E inv_a = field.Inv(a);
      EXPECT_TRUE(field.IsOne(field.Mul(a, inv_a)));
      EXPECT_EQ(field.Div(field.Mul(b, a), a), b);
    }

    // Sqrt: Sqrt(a)^2 == a and Sqrt(a^2) == a.
    const E sa = field.Sqrt(a);
    EXPECT_EQ(field.Sqr(sa), a);
    EXPECT_EQ(field.Sqrt(field.Sqr(a)), a);

    // Trace linearity & Frobenius invariance.
    EXPECT_EQ(field.Trace(field.Add(a, b)),
              static_cast<uint8_t>(field.Trace(a) ^ field.Trace(b)));
    EXPECT_EQ(field.Trace(field.Sqr(a)), field.Trace(a));

    // HalfTrace for odd M: for any element z = a^2 + a (guaranteed Tr(z) == 0),
    // HalfTrace(z)^2 + HalfTrace(z) == z.
    if constexpr (M % 2 == 1) {
      const E z = field.Add(field.Sqr(a), a);
      EXPECT_EQ(field.Trace(z), 0U);
      const E hz = field.HalfTrace(z);
      EXPECT_EQ(field.Add(field.Sqr(hz), hz), z);
    }
  }
}

TEST(LargeBinaryFieldTest, GF163GF233GF256GF283) {
  const auto gf163 = CreateGF2_163();
  const auto gf233 = CreateGF2_233();
  const auto gf256 = CreateGF2_256();
  const auto gf283 = CreateGF2_283();

  VerifyLargeBinaryField(gf163, 201);
  VerifyLargeBinaryField(gf233, 202);
  VerifyLargeBinaryField(gf256, 203);
  VerifyLargeBinaryField(gf283, 204);

  // Even-degree HalfTrace throws std::domain_error.
  EXPECT_THROW((void)gf256.HalfTrace(BigUInt<256>(1ULL)), std::domain_error);
}

}  // namespace
