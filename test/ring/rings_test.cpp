#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/ring_element.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/ring/ideal.hpp"
#include "xgalois/ring/integer_ring.hpp"
#include "xgalois/ring/matrix_ring.hpp"
#include "xgalois/ring/modular_ring.hpp"
#include "xgalois/ring/polynomial_ring.hpp"
#include "xgalois/ring/quadratic_integers.hpp"
#include "xgalois/ring/quotient_ring.hpp"
#include "xgalois/ring/skew_polynomial_ring.hpp"

namespace {

using xg::BigInt;
using xg::CommutativeRing;
using xg::EuclideanDomain;
using xg::Field;
using xg::GaloisFieldBinaryExtension;
using xg::GaloisFieldPrime;
using xg::IntegralDomain;
using xg::PrincipalIdealDomain;
using xg::Ring;
using xg::RingElement;
using xg::ring::ChineseRemainderTheorem;
using xg::ring::CyclicRing;
using xg::ring::EisensteinInt;
using xg::ring::EisensteinIntegerRing;
using xg::ring::FrobeniusSkewRing;
using xg::ring::GaussianInt;
using xg::ring::GaussianIntegerRing;
using xg::ring::IdealIntersection;
using xg::ring::IdealProduct;
using xg::ring::IdealSum;
using xg::ring::IntegerRing;
using xg::ring::MatrixRing;
using xg::ring::ModularRing;
using xg::ring::NegacyclicRing;
using xg::ring::PolynomialQuotientRing;
using xg::ring::PolynomialRing;
using xg::ring::PrincipalIdeal;
using xg::ring::QuotientRing;
using xg::ring::RingPoly;
using xg::ring::SkewPolynomialRing;

// ============================================================================
// Compile-Time Concept Hierarchy Verification
// ============================================================================

using ZRing = IntegerRing<BigInt>;
using Z64Ring = IntegerRing<int64_t>;
using Mod64Ring = ModularRing<uint64_t>;
using GF7 = GaloisFieldPrime<uint32_t>;
using GF256 = GaloisFieldBinaryExtension<uint8_t>;

static_assert(Ring<ZRing>);
static_assert(CommutativeRing<ZRing>);
static_assert(IntegralDomain<ZRing>);
static_assert(PrincipalIdealDomain<ZRing>);
static_assert(EuclideanDomain<ZRing>);
static_assert(!Field<ZRing>);

static_assert(EuclideanDomain<Z64Ring>);
static_assert(!Field<Z64Ring>);

static_assert(Ring<Mod64Ring>);
static_assert(CommutativeRing<Mod64Ring>);
static_assert(!IntegralDomain<Mod64Ring>);
static_assert(!EuclideanDomain<Mod64Ring>);

static_assert(EuclideanDomain<GaussianIntegerRing>);
static_assert(!Field<GaussianIntegerRing>);

static_assert(EuclideanDomain<EisensteinIntegerRing>);
static_assert(!Field<EisensteinIntegerRing>);

static_assert(CommutativeRing<PolynomialRing<Mod64Ring>>);
static_assert(!IntegralDomain<PolynomialRing<Mod64Ring>>);
static_assert(!EuclideanDomain<PolynomialRing<Mod64Ring>>);

static_assert(IntegralDomain<PolynomialRing<ZRing>>);
static_assert(!PrincipalIdealDomain<PolynomialRing<ZRing>>);
static_assert(!EuclideanDomain<PolynomialRing<ZRing>>);

static_assert(EuclideanDomain<PolynomialRing<GF7>>);
static_assert(!Field<PolynomialRing<GF7>>);

static_assert(Ring<PolynomialRing<MatrixRing<ZRing>>>);
static_assert(!CommutativeRing<PolynomialRing<MatrixRing<ZRing>>>);

static_assert(CommutativeRing<QuotientRing<ZRing>>);
static_assert(!IntegralDomain<QuotientRing<ZRing>>);

static_assert(CommutativeRing<PolynomialQuotientRing<Mod64Ring>>);
static_assert(CommutativeRing<NegacyclicRing<Mod64Ring>>);
static_assert(CommutativeRing<CyclicRing<Mod64Ring>>);

static_assert(Ring<MatrixRing<Mod64Ring>>);
static_assert(!CommutativeRing<MatrixRing<Mod64Ring>>);

static_assert(Ring<FrobeniusSkewRing<GF256>>);
static_assert(!CommutativeRing<FrobeniusSkewRing<GF256>>);

// ============================================================================
// IntegerRing Tests
// ============================================================================

TEST(IntegerRingTest, EuclideanDivisionAndCanonicalization) {
  const ZRing z;
  const std::vector<int64_t> dividends = {-25, -7, -1, 0, 1, 7, 25};
  const std::vector<int64_t> divisors = {-7, -3, -1, 1, 3, 7};

  for (int64_t a_val : dividends) {
    for (int64_t b_val : divisors) {
      const BigInt a(a_val);
      const BigInt b(b_val);
      auto [q, r] = z.DivRem(a, b);
      EXPECT_EQ(z.Add(z.Mul(q, b), r), a);
      EXPECT_GE(r, BigInt(0));
      EXPECT_LT(r, z.EuclideanNorm(b));
    }
  }

  EXPECT_THROW((void)z.DivRem(BigInt(5), BigInt(0)), std::domain_error);

  EXPECT_TRUE(z.IsUnit(BigInt(1)));
  EXPECT_TRUE(z.IsUnit(BigInt(-1)));
  EXPECT_FALSE(z.IsUnit(BigInt(0)));
  EXPECT_FALSE(z.IsUnit(BigInt(2)));
  EXPECT_EQ(z.TryInv(BigInt(-1)), std::optional<BigInt>(BigInt(-1)));
  EXPECT_EQ(z.TryInv(BigInt(3)), std::nullopt);

  auto [canon_neg, unit_neg] = z.AssociateCanonical(BigInt(-42));
  EXPECT_EQ(canon_neg, BigInt(42));
  EXPECT_EQ(unit_neg, BigInt(-1));
  EXPECT_EQ(z.Mul(unit_neg, BigInt(-42)), canon_neg);

  auto [g, s, t] = z.ExtendedGcd(BigInt(240), BigInt(-46));
  EXPECT_EQ(g, BigInt(2));
  EXPECT_EQ(z.Add(z.Mul(s, BigInt(240)), z.Mul(t, BigInt(-46))), g);

  // Interop with RingElement<IntegerRing<BigInt>>
  auto z_ptr = std::make_shared<const ZRing>();
  RingElement<ZRing> e1(BigInt(15), z_ptr);
  RingElement<ZRing> e2(BigInt(-4), z_ptr);
  EXPECT_EQ((e1 + e2).RawValue(), BigInt(11));
  EXPECT_EQ((e1 * e2).RawValue(), BigInt(-60));
  EXPECT_FALSE(e1.IsUnit());
}

// ============================================================================
// ModularRing Tests
// ============================================================================

TEST(ModularRingTest, ZMod12UnitsAndZeroDivisors) {
  const Mod64Ring z12(12);
  EXPECT_EQ(z12.Modulus(), 12ULL);
  EXPECT_EQ(z12.EulerTotient(), 4ULL);

  const std::vector<uint64_t> expected_units = {1, 5, 7, 11};
  EXPECT_EQ(z12.UnitGroupElements(), expected_units);

  for (uint64_t x = 0; x < 12; ++x) {
    const bool is_u = z12.IsUnit(x);
    const bool is_zd = z12.IsZeroDivisor(x);
    if (x == 0) {
      EXPECT_FALSE(is_u);
      EXPECT_FALSE(is_zd);
      EXPECT_EQ(z12.TryInv(x), std::nullopt);
    } else if (x == 1 || x == 5 || x == 7 || x == 11) {
      EXPECT_TRUE(is_u);
      EXPECT_FALSE(is_zd);
      auto inv = z12.TryInv(x);
      ASSERT_TRUE(inv.has_value());
      EXPECT_EQ(z12.Mul(x, *inv), 1ULL);
    } else {
      EXPECT_FALSE(is_u);
      EXPECT_TRUE(is_zd);
      EXPECT_EQ(z12.TryInv(x), std::nullopt);
    }
  }
}

TEST(ModularRingTest, ZMod17AndOverflowSafety) {
  const Mod64Ring z17(17);
  EXPECT_EQ(z17.EulerTotient(), 16ULL);
  for (uint64_t x = 1; x < 17; ++x) {
    EXPECT_TRUE(z17.IsUnit(x));
    EXPECT_FALSE(z17.IsZeroDivisor(x));
    auto inv = z17.TryInv(x);
    ASSERT_TRUE(inv.has_value());
    EXPECT_EQ(z17.Mul(x, *inv), 1ULL);
    EXPECT_EQ(z17.Pow(x, 16), 1ULL);
  }

  // Largest 64-bit prime: 2^64 - 59.
  const uint64_t big_p = std::numeric_limits<uint64_t>::max() - 58ULL;
  const Mod64Ring zbig(big_p);
  const uint64_t a = big_p - 1;
  const uint64_t b = big_p - 2;
  EXPECT_EQ(zbig.Add(a, b), big_p - 3);
  EXPECT_EQ(zbig.Sub(b, a), big_p - 1);
  EXPECT_EQ(zbig.Mul(a, b), 2ULL);
  auto inv_b = zbig.TryInv(b);
  ASSERT_TRUE(inv_b.has_value());
  EXPECT_EQ(zbig.Mul(b, *inv_b), 1ULL);
  EXPECT_EQ(zbig.EulerTotient(), big_p - 1);
}

// ============================================================================
// Gaussian & Eisenstein Integer Ring Tests
// ============================================================================

TEST(QuadraticIntegersTest, GaussianStressGridAndBezout) {
  const GaussianIntegerRing zi;

  // Verify units and canonical associates.
  const std::vector<GaussianInt> units = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
  for (const auto& u : units) {
    EXPECT_TRUE(zi.IsUnit(u));
    auto inv = zi.TryInv(u);
    ASSERT_TRUE(inv.has_value());
    EXPECT_EQ(zi.Mul(u, *inv), zi.MultiplicativeIdentity());
  }
  EXPECT_FALSE(zi.IsUnit({1, 1}));
  EXPECT_EQ(zi.TryInv({1, 1}), std::nullopt);

  // Stress-test DivRem norm decrease across a grid.
  for (int64_t yr = -4; yr <= 4; ++yr) {
    for (int64_t yi = -4; yi <= 4; ++yi) {
      if (yr == 0 && yi == 0) {
        continue;
      }
      const GaussianInt y{yr, yi};
      const BigInt ny = zi.Norm(y);
      for (int64_t xr = -6; xr <= 6; ++xr) {
        for (int64_t xi = -6; xi <= 6; ++xi) {
          const GaussianInt x{xr, xi};
          auto [q, r] = zi.DivRem(x, y);
          EXPECT_EQ(zi.Add(zi.Mul(q, y), r), x);
          EXPECT_LT(zi.Norm(r), ny);
        }
      }
    }
  }

  // ExtendedGcd + Bezout identity + first-quadrant canonicalization.
  const GaussianInt a{11, 3};
  const GaussianInt b{1, 8};
  auto [g, s, t] = zi.ExtendedGcd(a, b);
  EXPECT_EQ(zi.Add(zi.Mul(s, a), zi.Mul(t, b)), g);
  EXPECT_GT(g.re, BigInt(0));
  EXPECT_GE(g.im, BigInt(0));
  EXPECT_TRUE(zi.IsZero(zi.DivRem(a, g).second));
  EXPECT_TRUE(zi.IsZero(zi.DivRem(b, g).second));
  // gcd(11 + 3i, 1 + 8i) is associate to 1 - 2i -> canonical 2 + i (norm 5).
  EXPECT_EQ(g, (GaussianInt{2, 1}));
}

TEST(QuadraticIntegersTest, EisensteinStressGridUnitsAndBezout) {
  const EisensteinIntegerRing zw;

  // All 6 units: +/-1, +/-w, +/-(1 + w).
  const std::vector<EisensteinInt> units = {{1, 0},  {-1, 0}, {0, 1},
                                            {0, -1}, {1, 1},  {-1, -1}};
  for (const auto& u : units) {
    EXPECT_TRUE(zw.IsUnit(u));
    auto inv = zw.TryInv(u);
    ASSERT_TRUE(inv.has_value());
    EXPECT_EQ(zw.Mul(u, *inv), zw.MultiplicativeIdentity());
  }
  EXPECT_FALSE(zw.IsUnit({2, 1}));

  // All 6 associates of an element map to the same canonical representative.
  const EisensteinInt sample{5, 2};
  auto [canon0, u0] = zw.AssociateCanonical(sample);
  EXPECT_EQ(zw.Mul(u0, sample), canon0);
  for (const auto& u : units) {
    auto [c, cu] = zw.AssociateCanonical(zw.Mul(u, sample));
    EXPECT_EQ(c, canon0);
    EXPECT_EQ(zw.Mul(cu, zw.Mul(u, sample)), c);
  }

  // Stress-test DivRem exact 6-neighbor Voronoi guarantee: N(r) < N(y).
  for (int64_t ya = -4; ya <= 4; ++ya) {
    for (int64_t yb = -4; yb <= 4; ++yb) {
      if (ya == 0 && yb == 0) {
        continue;
      }
      const EisensteinInt y{ya, yb};
      const BigInt ny = zw.Norm(y);
      for (int64_t xa = -6; xa <= 6; ++xa) {
        for (int64_t xb = -6; xb <= 6; ++xb) {
          const EisensteinInt x{xa, xb};
          auto [q, r] = zw.DivRem(x, y);
          EXPECT_EQ(zw.Add(zw.Mul(q, y), r), x);
          EXPECT_LT(zw.Norm(r), ny);
        }
      }
    }
  }

  // Bezout identity over Z[w].
  const EisensteinInt a{7, 3};
  const EisensteinInt b{4, -2};
  auto [g, s, t] = zw.ExtendedGcd(a, b);
  EXPECT_EQ(zw.Add(zw.Mul(s, a), zw.Mul(t, b)), g);
  EXPECT_TRUE(zw.IsZero(zw.DivRem(a, g).second));
  EXPECT_TRUE(zw.IsZero(zw.DivRem(b, g).second));
}

// ============================================================================
// PolynomialRing Tests
// ============================================================================

TEST(PolynomialRingTest, OverZMod6ZeroDivisorCancellationAndMonicDivRem) {
  const Mod64Ring z6(6);
  const PolynomialRing<Mod64Ring> poly(z6);

  // (2x + 1)(3x + 1) = 6x^2 + 5x + 1 = 5x + 1 in (Z/6Z)[x].
  const RingPoly<Mod64Ring> p({1, 2});
  const RingPoly<Mod64Ring> q({1, 3});
  const auto prod = poly.Mul(p, q);
  EXPECT_EQ(poly.Degree(prod), 1);
  EXPECT_EQ(prod.Size(), 2U);  // trailing zero trimmed!
  EXPECT_EQ(prod[0], 1ULL);
  EXPECT_EQ(prod[1], 5ULL);

  // MonicDivRem over Z/6Z: divide (x^3 + 4x^2 + 2x + 5) by (x^2 + 2x + 3).
  const RingPoly<Mod64Ring> dividend({5, 2, 4, 1});
  const RingPoly<Mod64Ring> monic_divisor({3, 2, 1});
  auto [quot, rem] = poly.MonicDivRem(dividend, monic_divisor);
  EXPECT_TRUE(
      poly.Equals(poly.Add(poly.Mul(quot, monic_divisor), rem), dividend));
  EXPECT_LT(poly.Degree(rem), poly.Degree(monic_divisor));

  // Derivative and Eval check.
  // d/dx(x^3 + 4x^2 + 2x + 5) = 3x^2 + 2x + 2 mod 6.
  const auto deriv = poly.Derivative(dividend);
  EXPECT_TRUE(poly.Equals(deriv, RingPoly<Mod64Ring>({2, 2, 3})));
  EXPECT_EQ(poly.Eval(dividend, 2ULL), (5 + 4 + 16 + 8) % 6ULL);
}

TEST(PolynomialRingTest, PseudoDivRemOverZAndEgcdOverGF7) {
  // PseudoDivRem over Z[x].
  const PolynomialRing<ZRing> zx;
  const RingPoly<ZRing> a({BigInt(1), BigInt(2), BigInt(3)});  // 3x^2 + 2x + 1
  const RingPoly<ZRing> b({BigInt(1), BigInt(2)});             // 2x + 1
  auto [pq, pr] = zx.PseudoDivRem(a, b);
  // lc(b)^(2 - 1 + 1) = 2^2 = 4.
  EXPECT_TRUE(zx.Equals(zx.ScalarMul(BigInt(4), a), zx.Add(zx.Mul(pq, b), pr)));
  EXPECT_LT(zx.Degree(pr), zx.Degree(b));

  // ExtendedGcd over GF(7)[x].
  const GF7 gf7(7);
  const PolynomialRing<GF7> f7x(gf7);
  // f = (x - 1)(x - 2) = x^2 - 3x + 2 = x^2 + 4x + 2 mod 7.
  // g = (x - 2)(x - 3) = x^2 - 5x + 6 = x^2 + 2x + 6 mod 7.
  const RingPoly<GF7> f({2, 4, 1});
  const RingPoly<GF7> g({6, 2, 1});
  auto [gcd_fg, s, t] = f7x.ExtendedGcd(f, g);
  // Monic gcd should be x - 2 = x + 5 mod 7.
  EXPECT_TRUE(f7x.Equals(gcd_fg, RingPoly<GF7>({5, 1})));
  EXPECT_TRUE(f7x.Equals(f7x.Add(f7x.Mul(s, f), f7x.Mul(t, g)), gcd_fg));
}

// ============================================================================
// QuotientRing, PolynomialQuotientRing, NegacyclicRing, CyclicRing Tests
// ============================================================================

TEST(QuotientRingTest, QuotientNegacyclicAndCyclicRings) {
  // QuotientRing over Z[i]: Z[i] / (2 + i) is a field of size 5!
  const GaussianIntegerRing zi;
  const QuotientRing<GaussianIntegerRing> zi_mod(zi, GaussianInt{2, 1});
  for (int64_t r = 1; r <= 4; ++r) {
    const GaussianInt elem{r, 0};
    EXPECT_TRUE(zi_mod.IsUnit(elem));
    auto inv = zi_mod.TryInv(elem);
    ASSERT_TRUE(inv.has_value());
    EXPECT_TRUE(zi_mod.IsOne(zi_mod.Mul(elem, *inv)));
  }

  // NegacyclicRing: Z_17[x] / (x^4 + 1).
  const Mod64Ring z17(17);
  const NegacyclicRing<Mod64Ring> neg_ring(z17, 4);
  const std::vector<uint64_t> x = {0, 1, 0, 0};
  const auto x4 = neg_ring.Pow(x, 4);
  EXPECT_TRUE(neg_ring.Equals(x4, {16, 0, 0, 0}));  // x^4 == -1 mod 17
  EXPECT_TRUE(neg_ring.IsOne(neg_ring.Pow(x, 8)));

  // Compare NegacyclicRing::Mul against PolynomialQuotientRing mod (x^4 + 1).
  const PolynomialQuotientRing<Mod64Ring> pq_neg(
      z17, RingPoly<Mod64Ring>({1, 0, 0, 0, 1}));
  const std::vector<uint64_t> av = {3, 5, 7, 11};
  const std::vector<uint64_t> bv = {2, 13, 4, 9};
  const auto neg_prod = neg_ring.Mul(av, bv);
  const auto pq_prod =
      pq_neg.Mul(RingPoly<Mod64Ring>(av), RingPoly<Mod64Ring>(bv));
  for (std::size_t i = 0; i < 4; ++i) {
    const uint64_t expected = (i < pq_prod.Size()) ? pq_prod[i] : 0ULL;
    EXPECT_EQ(neg_prod[i], expected);
  }

  // CyclicRing: Z_17[x] / (x^4 - 1).
  const CyclicRing<Mod64Ring> cyc_ring(z17, 4);
  EXPECT_TRUE(cyc_ring.IsOne(cyc_ring.Pow(x, 4)));
  const PolynomialQuotientRing<Mod64Ring> pq_cyc(
      z17, RingPoly<Mod64Ring>({16, 0, 0, 0, 1}));
  const auto cyc_prod = cyc_ring.Mul(av, bv);
  const auto pq_cyc_prod =
      pq_cyc.Mul(RingPoly<Mod64Ring>(av), RingPoly<Mod64Ring>(bv));
  for (std::size_t i = 0; i < 4; ++i) {
    const uint64_t expected = (i < pq_cyc_prod.Size()) ? pq_cyc_prod[i] : 0ULL;
    EXPECT_EQ(cyc_prod[i], expected);
  }
}

// ============================================================================
// MatrixRing & Berkowitz Determinant Tests
// ============================================================================

TEST(MatrixRingTest, BerkowitzDetOverZMod6MatchesSarrus) {
  const Mod64Ring z6(6);
  const MatrixRing<Mod64Ring> m3(z6, 3);

  auto sarrus_det = [&](const std::vector<uint64_t>& m) -> uint64_t {
    // Explicit 3x3 Sarrus rule over Z/6Z.
    const uint64_t pos1 = z6.Mul(z6.Mul(m[0], m[4]), m[8]);
    const uint64_t pos2 = z6.Mul(z6.Mul(m[1], m[5]), m[6]);
    const uint64_t pos3 = z6.Mul(z6.Mul(m[2], m[3]), m[7]);
    const uint64_t neg1 = z6.Mul(z6.Mul(m[2], m[4]), m[6]);
    const uint64_t neg2 = z6.Mul(z6.Mul(m[0], m[5]), m[7]);
    const uint64_t neg3 = z6.Mul(z6.Mul(m[1], m[3]), m[8]);
    const uint64_t pos = z6.Add(z6.Add(pos1, pos2), pos3);
    const uint64_t neg = z6.Add(z6.Add(neg1, neg2), neg3);
    return z6.Sub(pos, neg);
  };

  // Matrix whose entries are ALL zero divisors in Z/6Z ({2, 3, 4}).
  const std::vector<uint64_t> a = {2, 3, 4, 3, 2, 3, 4, 3, 2};
  const std::vector<uint64_t> b = {1, 2, 3, 4, 5, 2, 3, 1, 4};
  const std::vector<uint64_t> c = {5, 3, 2, 2, 4, 1, 3, 5, 3};

  EXPECT_EQ(m3.Det(a), sarrus_det(a));
  EXPECT_EQ(m3.Det(b), sarrus_det(b));
  EXPECT_EQ(m3.Det(c), sarrus_det(c));

  // Multiplicativity: det(B * C) == det(B) * det(C) in Z/6Z.
  EXPECT_EQ(m3.Det(m3.Mul(b, c)), z6.Mul(m3.Det(b), m3.Det(c)));
  // Non-commutativity: B * C != C * B.
  EXPECT_FALSE(m3.Equals(m3.Mul(b, c), m3.Mul(c, b)));
  EXPECT_EQ(m3.Trace(b), z6.Add(z6.Add(1, 5), 4));
  EXPECT_TRUE(m3.Equals(m3.Transpose(m3.Transpose(b)), b));
  EXPECT_TRUE(m3.Equals(m3.Pow(b, 3), m3.Mul(m3.Mul(b, b), b)));

  // 4x4 determinant over Z.
  const MatrixRing<ZRing> m4(ZRing{}, 4);
  const std::vector<BigInt> d4 = {
      BigInt(1), BigInt(2),  BigInt(0), BigInt(1),   //
      BigInt(3), BigInt(-1), BigInt(2), BigInt(0),   //
      BigInt(0), BigInt(4),  BigInt(1), BigInt(-2),  //
      BigInt(2), BigInt(0),  BigInt(3), BigInt(5)};
  // Reference det computed analytically: -129.
  // Check via Laplace expansion on row 0:
  // 1 * det([-1,2,0; 4,1,-2; 0,3,5]) = 1*(-1*(5+6) - 2*(20-0)) = -11 - 40 = -51
  // -2 * det([3,2,0; 0,1,-2; 2,3,5]) = -2*(3*(5+6) - 2*(0+4)) = -2*(33 - 8) =
  // -50 -1 * det([3,-1,2; 0,4,1; 2,0,3]) = -1*(3*(12-0) + 1*(0-2) + 2*(0-8))
  // = -1*(36 - 2 - 16) = -18. Total = -51 - 50 - 18 = -119.
  EXPECT_EQ(m4.Det(d4), BigInt(-119));
}

// ============================================================================
// FrobeniusSkewRing & SkewPolynomialRing Tests
// ============================================================================

TEST(SkewPolynomialRingTest, FrobeniusSkewOverGF256) {
  // GF(2^8) with AES irreducible polynomial x^8 + x^4 + x^3 + x + 1.
  const GF256 gf256(8, "int", "x^8 + x^4 + x^3 + x + 1");
  const FrobeniusSkewRing<GF256> skew(gf256, /*frobenius_step=*/1);

  using SkewPoly = FrobeniusSkewRing<GF256>::element_type;

  // Non-commutativity: x * alpha == alpha^2 * x != alpha * x.
  const uint8_t alpha = 0x02;
  const SkewPoly x = skew.X();
  const SkewPoly c_alpha = skew.Constant(alpha);
  const SkewPoly xa = skew.Mul(x, c_alpha);
  const SkewPoly ax = skew.Mul(c_alpha, x);
  EXPECT_FALSE(skew.Equals(xa, ax));
  EXPECT_TRUE(skew.Equals(xa, SkewPoly{0, gf256.Mul(alpha, alpha)}));

  // RightDivRem and LeftDivRem.
  const SkewPoly a = {0x13, 0x07, 0x25, 0x41, 0x09};
  const SkewPoly b = {0x05, 0x11, 0x03};

  auto [qr, rr] = skew.RightDivRem(a, b);
  EXPECT_TRUE(skew.Equals(skew.Add(skew.Mul(qr, b), rr), a));
  EXPECT_LT(skew.Degree(rr), skew.Degree(b));

  auto [ql, rl] = skew.LeftDivRem(a, b);
  EXPECT_TRUE(skew.Equals(skew.Add(skew.Mul(b, ql), rl), a));
  EXPECT_LT(skew.Degree(rl), skew.Degree(b));

  // RightExtendedGcd (GCRD) with known common right factor.
  const SkewPoly common_r = {0x0b, 0x17, 0x01};  // monic deg 2
  const SkewPoly p1 = skew.Mul(SkewPoly{0x03, 0x05, 0x01}, common_r);
  const SkewPoly p2 = skew.Mul(SkewPoly{0x19, 0x02}, common_r);
  auto [gcrd, ur, vr] = skew.RightExtendedGcd(p1, p2);
  EXPECT_TRUE(skew.Equals(skew.Add(skew.Mul(ur, p1), skew.Mul(vr, p2)), gcrd));
  EXPECT_TRUE(gf256.IsOne(skew.LeadingCoeff(gcrd)));
  EXPECT_TRUE(skew.IsZero(skew.RightDivRem(p1, gcrd).second));
  EXPECT_TRUE(skew.IsZero(skew.RightDivRem(p2, gcrd).second));
  EXPECT_GE(skew.Degree(gcrd), skew.Degree(common_r));

  // Partial Skew-EEA (stop_deg = 2).
  auto [partial_r, partial_u, partial_v] =
      skew.RightExtendedGcd(p1, p2, /*stop_deg=*/2);
  EXPECT_LE(skew.Degree(partial_r), 2);
  EXPECT_TRUE(skew.Equals(
      skew.Add(skew.Mul(partial_u, p1), skew.Mul(partial_v, p2)), partial_r));

  // LeftExtendedGcd (GCLD) with known common left factor.
  const SkewPoly common_l = {0x0d, 0x21, 0x01};  // monic deg 2
  const SkewPoly q1 = skew.Mul(common_l, SkewPoly{0x07, 0x04, 0x01});
  const SkewPoly q2 = skew.Mul(common_l, SkewPoly{0x15, 0x03});
  auto [gcld, ul, vl] = skew.LeftExtendedGcd(q1, q2);
  EXPECT_TRUE(skew.Equals(skew.Add(skew.Mul(q1, ul), skew.Mul(q2, vl)), gcld));
  EXPECT_TRUE(gf256.IsOne(skew.LeadingCoeff(gcld)));
  EXPECT_TRUE(skew.IsZero(skew.LeftDivRem(q1, gcld).second));
  EXPECT_TRUE(skew.IsZero(skew.LeftDivRem(q2, gcld).second));
  EXPECT_GE(skew.Degree(gcld), skew.Degree(common_l));

  // LCLM and LCRM.
  const SkewPoly lclm = skew.Lclm(p1, p2);
  EXPECT_TRUE(gf256.IsOne(skew.LeadingCoeff(lclm)));
  EXPECT_TRUE(skew.IsZero(skew.RightDivRem(lclm, p1).second));
  EXPECT_TRUE(skew.IsZero(skew.RightDivRem(lclm, p2).second));
  EXPECT_EQ(skew.Degree(lclm),
            skew.Degree(p1) + skew.Degree(p2) - skew.Degree(gcrd));

  const SkewPoly lcrm = skew.Lcrm(q1, q2);
  EXPECT_TRUE(gf256.IsOne(skew.LeadingCoeff(lcrm)));
  EXPECT_TRUE(skew.IsZero(skew.LeftDivRem(lcrm, q1).second));
  EXPECT_TRUE(skew.IsZero(skew.LeftDivRem(lcrm, q2).second));
  EXPECT_EQ(skew.Degree(lcrm),
            skew.Degree(q1) + skew.Degree(q2) - skew.Degree(gcld));

  // Operator composition law: (a * b)(z) == a(b(z)).
  for (uint8_t z :
       {uint8_t{0x01}, uint8_t{0x02}, uint8_t{0x53}, uint8_t{0xa7}}) {
    EXPECT_EQ(skew.EvalOperator(skew.Mul(a, b), z),
              skew.EvalOperator(a, skew.EvalOperator(b, z)));
  }

  // MinimalPolynomial over 4 F_2-linearly independent points in GF(2^8):
  // {1, alpha, alpha^2, alpha^3} = {0x01, 0x02, 0x04, 0x08}.
  const std::vector<uint8_t> basis = {0x01, 0x02, 0x04, 0x08};
  const SkewPoly min_poly = skew.MinimalPolynomial(basis);
  EXPECT_EQ(skew.Degree(min_poly), 4);
  EXPECT_TRUE(gf256.IsOne(skew.LeadingCoeff(min_poly)));

  // Must vanish on ALL 2^4 = 16 F_2-linear combinations of basis.
  for (uint32_t mask = 0; mask < 16; ++mask) {
    uint8_t pt = 0;
    for (std::size_t bit = 0; bit < 4; ++bit) {
      if (((mask >> bit) & 1U) != 0) {
        pt = gf256.Add(pt, basis[bit]);
      }
    }
    EXPECT_EQ(skew.EvalOperator(min_poly, pt), 0U);
  }
  // Must NOT vanish on a point outside the span (e.g. alpha^4 = 0x10).
  EXPECT_NE(skew.EvalOperator(min_poly, 0x10), 0U);

  // InterpolateOperator on the 4 independent points.
  const std::vector<uint8_t> values = {0x31, 0x77, 0xc4, 0x09};
  const SkewPoly interp = skew.InterpolateOperator(basis, values);
  EXPECT_LT(skew.Degree(interp), 4);
  for (std::size_t i = 0; i < basis.size(); ++i) {
    EXPECT_EQ(skew.EvalOperator(interp, basis[i]), values[i]);
  }
}

TEST(SkewPolynomialRingTest, GeneralOreRingWithNonZeroDerivation) {
  const GF256 gf256(8, "int", "x^8 + x^4 + x^3 + x + 1");
  const uint8_t c = 0x53;

  // Inner sigma-derivation: delta(a) = sigma(a)*c - c*a.
  auto sigma = [](const GF256& f, uint8_t a) { return f.Frobenius(a); };
  auto delta = [c](const GF256& f, uint8_t a) {
    return f.Sub(f.Mul(f.Frobenius(a), c), f.Mul(c, a));
  };
  auto sigma_inv = [](const GF256& f, uint8_t a) {
    uint8_t cur = a;
    for (int i = 0; i < 7; ++i) {
      cur = f.Frobenius(cur);
    }
    return cur;
  };

  const SkewPolynomialRing<GF256, decltype(sigma), decltype(delta),
                           decltype(sigma_inv)>
      ore(gf256, sigma, delta, sigma_inv);

  using Poly = std::vector<uint8_t>;
  const Poly f = {0x11, 0x22, 0x03};
  const Poly g = {0x05, 0x09, 0x07};
  const Poly h = {0x1a, 0x04};

  // Associativity (f * g) * h == f * (g * h).
  EXPECT_TRUE(ore.Equals(ore.Mul(ore.Mul(f, g), h), ore.Mul(f, ore.Mul(g, h))));

  // Right and Left division with non-zero derivation.
  const Poly prod_a = ore.Add(ore.Mul(f, g), Poly{0x12, 0x34});
  auto [qr, rr] = ore.RightDivRem(prod_a, g);
  EXPECT_TRUE(ore.Equals(ore.Add(ore.Mul(qr, g), rr), prod_a));
  EXPECT_LT(ore.Degree(rr), ore.Degree(g));

  auto [ql, rl] = ore.LeftDivRem(prod_a, g);
  EXPECT_TRUE(ore.Equals(ore.Add(ore.Mul(g, ql), rl), prod_a));
  EXPECT_LT(ore.Degree(rl), ore.Degree(g));
}

// ============================================================================
// PrincipalIdeal & ChineseRemainderTheorem Tests
// ============================================================================

TEST(IdealAndCrtTest, PrincipalIdealsAndCrtOverZ) {
  const ZRing z;

  const PrincipalIdeal<ZRing> i12(z, BigInt(-12));
  const PrincipalIdeal<ZRing> i12_pos(z, BigInt(12));
  const PrincipalIdeal<ZRing> i18(z, BigInt(18));
  const PrincipalIdeal<ZRing> i0(z, BigInt(0));
  const PrincipalIdeal<ZRing> i_whole(z, BigInt(-1));

  EXPECT_EQ(i12, i12_pos);
  EXPECT_EQ(i12.Generator(), BigInt(12));
  EXPECT_TRUE(i12.Contains(BigInt(36)));
  EXPECT_FALSE(i12.Contains(BigInt(18)));
  EXPECT_TRUE(i0.IsZero());
  EXPECT_TRUE(i_whole.IsWholeRing());
  EXPECT_FALSE(i12.IsWholeRing());

  EXPECT_EQ(IdealSum(i12, i18), PrincipalIdeal<ZRing>(z, BigInt(6)));
  EXPECT_EQ(IdealProduct(i12, i18), PrincipalIdeal<ZRing>(z, BigInt(216)));
  EXPECT_EQ(IdealIntersection(i12, i18), PrincipalIdeal<ZRing>(z, BigInt(36)));

  // CRT: pairwise-coprime moduli (3, 5, 7).
  auto res1 = ChineseRemainderTheorem(z, {BigInt(2), BigInt(3), BigInt(2)},
                                      {BigInt(3), BigInt(5), BigInt(7)});
  ASSERT_TRUE(res1.has_value());
  EXPECT_EQ(res1->first, BigInt(23));
  EXPECT_EQ(res1->second, BigInt(105));

  // CRT: non-coprime consistent (x == 3 mod 6, x == 5 mod 8 => x == 21 mod 24).
  auto res2 = ChineseRemainderTheorem(z, {BigInt(3), BigInt(5)},
                                      {BigInt(6), BigInt(8)});
  ASSERT_TRUE(res2.has_value());
  EXPECT_EQ(res2->first, BigInt(21));
  EXPECT_EQ(res2->second, BigInt(24));

  // CRT: non-coprime inconsistent (x == 1 mod 6, x == 2 mod 9).
  auto res3 = ChineseRemainderTheorem(z, {BigInt(1), BigInt(2)},
                                      {BigInt(6), BigInt(9)});
  EXPECT_EQ(res3, std::nullopt);
}

TEST(IdealAndCrtTest, PrincipalIdealsAndCrtOverGF7Poly) {
  const GF7 gf7(7);
  const PolynomialRing<GF7> f7x(gf7);
  using Poly = RingPoly<GF7>;

  // m1 = (x - 1)(x - 2) = x^2 + 4x + 2
  // m2 = (x - 2)(x - 3) = x^2 + 2x + 6
  const Poly m1({2, 4, 1});
  const Poly m2({6, 2, 1});
  const Poly m3({4, 1});  // x - 3 = x + 4 (coprime to m1)

  const PrincipalIdeal<PolynomialRing<GF7>> id1(f7x, m1);
  const PrincipalIdeal<PolynomialRing<GF7>> id2(f7x, m2);
  // Sum = (x - 2) = (x + 5)
  EXPECT_EQ(IdealSum(id1, id2),
            PrincipalIdeal<PolynomialRing<GF7>>(f7x, Poly({5, 1})));
  // Intersection = (x - 1)(x - 2)(x - 3)
  const Poly lcm_expected = f7x.Mul(m1, m3);
  EXPECT_EQ(IdealIntersection(id1, id2),
            PrincipalIdeal<PolynomialRing<GF7>>(f7x, lcm_expected));

  // Coprime polynomial CRT: mod m1 and mod m3.
  const Poly target({3, 1, 5});  // 5x^2 + x + 3
  const Poly r1 = f7x.DivRem(target, m1).second;
  const Poly r3 = f7x.DivRem(target, m3).second;
  auto crt_coprime = ChineseRemainderTheorem(f7x, {r1, r3}, {m1, m3});
  ASSERT_TRUE(crt_coprime.has_value());
  EXPECT_TRUE(f7x.Equals(crt_coprime->first, target));
  EXPECT_TRUE(f7x.Equals(crt_coprime->second, lcm_expected));

  // Non-coprime consistent polynomial CRT: mod m1 and mod m2 (share x - 2).
  const Poly r2_ok = f7x.DivRem(target, m2).second;
  auto crt_consistent = ChineseRemainderTheorem(f7x, {r1, r2_ok}, {m1, m2});
  ASSERT_TRUE(crt_consistent.has_value());
  EXPECT_TRUE(f7x.Equals(crt_consistent->first, target));
  EXPECT_TRUE(f7x.Equals(crt_consistent->second, lcm_expected));

  // Non-coprime inconsistent: perturb r2 by +1 so values at x=2 disagree.
  const Poly r2_bad = f7x.Add(r2_ok, Poly({1}));
  auto crt_bad = ChineseRemainderTheorem(f7x, {r1, r2_bad}, {m1, m2});
  EXPECT_EQ(crt_bad, std::nullopt);
}

}  // namespace
