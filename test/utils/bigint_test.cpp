#include "xgalois/utils/bigint.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using xg::BigInt;
using xg::BigUInt;
using xg::DivRem;
using xg::EuclideanDivRem;
using xg::ExtendedGcd;
using xg::Factorial;
using xg::FactorInteger;
using xg::Gcd;
using xg::IsProbablePrime;
using xg::ISqrt;
using xg::JacobiSymbol;
using xg::ModAdd;
using xg::ModInv;
using xg::ModInverse;
using xg::ModMul;
using xg::ModPow;
using xg::ModSub;
using xg::MontgomeryCtx;
using xg::PowMod;

namespace {

uint64_t NextRng(uint64_t& state) {
  uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

template <std::size_t Bits>
BigUInt<Bits> RandomBigUInt(uint64_t& state) {
  std::array<uint64_t, BigUInt<Bits>::kLimbs> limbs{};
  for (std::size_t i = 0; i < BigUInt<Bits>::kLimbs; ++i) {
    limbs[i] = NextRng(state);
  }
  return BigUInt<Bits>(limbs);
}

BigInt RandomBigIntLimbs(std::size_t num_limbs, uint64_t& state) {
  BigInt res(0);
  for (std::size_t i = 0; i < num_limbs; ++i) {
    uint64_t w = NextRng(state);
    if (i == 0 && w == 0) {
      w = 1;
    }
    res = (res << 64) + BigInt::FromUint64(w);
  }
  return res;
}

template <std::size_t Bits>
void VerifyMontgomerySuite(const BigUInt<Bits>& p, uint64_t seed) {
  MontgomeryCtx<Bits> ctx(p);
  EXPECT_EQ(ctx.Modulus(), p);

  const BigUInt<Bits> zero = BigUInt<Bits>::Zero();
  const BigUInt<Bits> one = BigUInt<Bits>::One();
  const BigUInt<Bits> p_minus_1 = p - one;
  const BigUInt<Bits> p_minus_2 = p - BigUInt<Bits>(2);

  // Round-trip boundary values.
  EXPECT_EQ(ctx.FromMont(ctx.ToMont(zero)), zero);
  EXPECT_EQ(ctx.FromMont(ctx.ToMont(one)), one);
  EXPECT_EQ(ctx.FromMont(ctx.ToMont(p_minus_1)), p_minus_1);
  EXPECT_EQ(ctx.FromMont(ctx.ToMont(p)), zero);

  // Negation and additive identities.
  EXPECT_EQ(ctx.FromMont(ctx.Neg(ctx.ToMont(zero))), zero);
  EXPECT_EQ(ctx.FromMont(ctx.Neg(ctx.ToMont(one))), p_minus_1);
  EXPECT_EQ(ctx.FromMont(ctx.Add(ctx.ToMont(p_minus_1), ctx.ToMont(one))),
            zero);
  EXPECT_EQ(ctx.FromMont(ctx.Sub(ctx.ToMont(zero), ctx.ToMont(one))),
            p_minus_1);

  // Fermat's Little Theorem and inverse on boundary values.
  EXPECT_EQ(ctx.Pow(ctx.ToMont(p_minus_1), p_minus_1), ctx.ToMont(one));
  EXPECT_EQ(ctx.Mul(ctx.ToMont(p_minus_1), ctx.Inv(ctx.ToMont(p_minus_1))),
            ctx.ToMont(one));

  uint64_t rng = seed;
  for (int i = 0; i < 20; ++i) {
    const BigUInt<Bits> a = RandomBigUInt<Bits>(rng) % p;
    const BigUInt<Bits> b = RandomBigUInt<Bits>(rng) % p;
    const BigUInt<Bits> a_m = ctx.ToMont(a);
    const BigUInt<Bits> b_m = ctx.ToMont(b);

    EXPECT_EQ(ctx.FromMont(a_m), a);
    EXPECT_EQ(ctx.FromMont(ctx.Add(a_m, b_m)), ModAdd(a, b, p));
    EXPECT_EQ(ctx.FromMont(ctx.Sub(a_m, b_m)), ModSub(a, b, p));
    EXPECT_EQ(ctx.FromMont(ctx.Mul(a_m, b_m)), ModMul(a, b, p));
    EXPECT_EQ(ctx.FromMont(ctx.Sqr(a_m)), ModMul(a, a, p));

    if (!a.IsZero()) {
      const BigUInt<Bits> a_inv_m = ctx.Inv(a_m);
      EXPECT_EQ(ctx.Mul(a_m, a_inv_m), ctx.ToMont(one));
      EXPECT_EQ(ctx.FromMont(a_inv_m), ModInv(a, p));
    }
  }

  // Verify Pow against Fermat's Little Theorem (a^(p-2) == a^{-1} mod p).
  const BigUInt<Bits> sample = RandomBigUInt<Bits>(rng) % p;
  if (!sample.IsZero()) {
    const BigUInt<Bits> s_m = ctx.ToMont(sample);
    EXPECT_EQ(ctx.Pow(s_m, p_minus_1), ctx.ToMont(one));
    EXPECT_EQ(ctx.Pow(s_m, p_minus_2), ctx.Inv(s_m));
  }
}

bool IsPrimeTrialSmall(uint64_t n) {
  if (n < 2) return false;
  for (uint64_t d = 2; d * d <= n; ++d) {
    if (n % d == 0) return false;
  }
  return true;
}

}  // namespace

// ============================================================================
// BigUInt<Bits> tests
// ============================================================================

TEST(BigUIntTest, ConstructionParsingAndFormatting) {
  const BigUInt<128> zero;
  EXPECT_TRUE(zero.IsZero());
  EXPECT_FALSE(zero.IsOne());
  EXPECT_TRUE(zero.IsEven());
  EXPECT_FALSE(zero.IsOdd());
  EXPECT_EQ(zero.BitLength(), 0U);
  EXPECT_EQ(zero.ToHex(), "0");
  EXPECT_EQ(zero.ToHex(true), "0x0");
  EXPECT_EQ(zero.ToDecimal(), "0");
  EXPECT_EQ(zero.ToString(), "0");

  const BigUInt<128> one(1);
  EXPECT_FALSE(one.IsZero());
  EXPECT_TRUE(one.IsOne());
  EXPECT_TRUE(one.IsOdd());
  EXPECT_EQ(one.BitLength(), 1U);
  EXPECT_TRUE(one.GetBit(0));
  EXPECT_FALSE(one.GetBit(1));
  EXPECT_FALSE(one.GetBit(128));

  const auto max128 =
      BigUInt<128>::FromHex("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
  EXPECT_EQ(max128.BitLength(), 128U);
  EXPECT_EQ(max128.ToDecimal(), "340282366920938463463374607431768211455");
  EXPECT_EQ(BigUInt<128>::FromDecimal(max128.ToDecimal()), max128);
  EXPECT_EQ(BigUInt<128>::FromHex("0x0000000000000000000000000000000001"), one);

  EXPECT_THROW(BigUInt<128>::FromHex(""), std::invalid_argument);
  EXPECT_THROW(BigUInt<128>::FromHex("0x"), std::invalid_argument);
  EXPECT_THROW(BigUInt<128>::FromHex("0x123g"), std::invalid_argument);
  EXPECT_THROW(BigUInt<128>::FromHex("0x100000000000000000000000000000000"),
               std::overflow_error);

  EXPECT_THROW(BigUInt<128>::FromDecimal(""), std::invalid_argument);
  EXPECT_THROW(BigUInt<128>::FromDecimal("-5"), std::invalid_argument);
  EXPECT_THROW(
      BigUInt<128>::FromDecimal("340282366920938463463374607431768211456"),
      std::overflow_error);
}

TEST(BigUIntTest, BitwiseShiftsAndComparisons) {
  BigUInt<256> x;
  x.SetBit(0);
  x.SetBit(64);
  x.SetBit(191);
  x.SetBit(255);
  EXPECT_EQ(x.BitLength(), 256U);
  EXPECT_TRUE(x.GetBit(0));
  EXPECT_TRUE(x.GetBit(64));
  EXPECT_TRUE(x.GetBit(191));
  EXPECT_TRUE(x.GetBit(255));
  EXPECT_FALSE(x.GetBit(100));
  EXPECT_THROW(x.SetBit(256), std::out_of_range);
  EXPECT_THROW((void)x.Limb(4), std::out_of_range);

  // Ordering checks MSB-first (not lexicographical LSW-first).
  const BigUInt<128> a({0xFFFFFFFFFFFFFFFFULL, 0ULL});
  const BigUInt<128> b({0ULL, 1ULL});
  EXPECT_LT(a, b);
  EXPECT_GT(b, a);
  EXPECT_NE(a, b);

  // Shifts across limb boundaries.
  const BigUInt<256> one(1);
  EXPECT_EQ((one << 64).Limb(1), 1ULL);
  EXPECT_EQ((one << 127).BitLength(), 128U);
  EXPECT_EQ(((one << 200) >> 136), (one << 64));
  EXPECT_EQ((one << 256), BigUInt<256>::Zero());
  EXPECT_EQ((x >> 256), BigUInt<256>::Zero());
  EXPECT_EQ((x & ~x), BigUInt<256>::Zero());
  EXPECT_EQ((x | ~x), ~BigUInt<256>::Zero());
  EXPECT_EQ((x ^ x), BigUInt<256>::Zero());
}

TEST(BigUIntTest, ArithmeticWideMulSqrAndDivRem) {
  // Carry/borrow across all limbs in 128, 256, 384, 512 bits.
  const BigUInt<256> max256 = ~BigUInt<256>::Zero();
  const auto [sum_wrap, carry] =
      BigUInt<256>::AddWithCarry(max256, BigUInt<256>::One());
  EXPECT_EQ(sum_wrap, BigUInt<256>::Zero());
  EXPECT_EQ(carry, 1ULL);

  const auto [sub_wrap, borrow] =
      BigUInt<256>::SubWithBorrow(BigUInt<256>::Zero(), BigUInt<256>::One());
  EXPECT_EQ(sub_wrap, max256);
  EXPECT_EQ(borrow, 1ULL);

  // MulWide and SqrWide agreement across 128, 256, 384, 512 bits.
  uint64_t rng = 0x123456789abcdef0ULL;
  auto check_width = [&rng]<std::size_t Bits>() {
    const BigUInt<Bits> max_v = ~BigUInt<Bits>::Zero();
    EXPECT_EQ(max_v.SqrWide(), max_v.MulWide(max_v));

    for (int i = 0; i < 15; ++i) {
      const BigUInt<Bits> a = RandomBigUInt<Bits>(rng);
      const BigUInt<Bits> b = RandomBigUInt<Bits>(rng);

      // SqrWide == MulWide(a, a)
      EXPECT_EQ(a.SqrWide(), a.MulWide(a));

      // Truncated multiplication matches low limbs of MulWide.
      const BigUInt<Bits> trunc = a * b;
      const BigUInt<2 * Bits> wide = a.MulWide(b);
      for (std::size_t k = 0; k < BigUInt<Bits>::kLimbs; ++k) {
        EXPECT_EQ(trunc.Limb(k), wide.Limb(k));
      }

      // DivRem invariant: a == q * b + r and r < b.
      if (!b.IsZero()) {
        auto [q, r] = BigUInt<Bits>::DivRem(a, b);
        EXPECT_LT(r, b);
        EXPECT_EQ(q * b + r, a);
        EXPECT_EQ(a / b, q);
        EXPECT_EQ(a % b, r);
      }

      // Also test division with a smaller multi-limb divisor (~half width).
      const BigUInt<Bits> small_b = b >> (Bits / 2);
      if (!small_b.IsZero()) {
        auto [q, r] = BigUInt<Bits>::DivRem(a, small_b);
        EXPECT_LT(r, small_b);
        EXPECT_EQ(q * small_b + r, a);
      }
    }
  };

  check_width.template operator()<128>();
  check_width.template operator()<256>();
  check_width.template operator()<384>();
  check_width.template operator()<512>();

  EXPECT_THROW((void)(BigUInt<256>(10) / BigUInt<256>::Zero()),
               std::domain_error);
}

TEST(BigUIntTest, ConstantTimeSelectAndSwap) {
  uint64_t rng = 0xcafebabe01234567ULL;
  BigUInt<256> a = RandomBigUInt<256>(rng);
  BigUInt<256> b = RandomBigUInt<256>(rng);
  const BigUInt<256> orig_a = a;
  const BigUInt<256> orig_b = b;

  EXPECT_EQ(BigUInt<256>::CtSelect(true, a, b), orig_a);
  EXPECT_EQ(BigUInt<256>::CtSelect(false, a, b), orig_b);

  BigUInt<256>::CtSwap(false, a, b);
  EXPECT_EQ(a, orig_a);
  EXPECT_EQ(b, orig_b);

  BigUInt<256>::CtSwap(true, a, b);
  EXPECT_EQ(a, orig_b);
  EXPECT_EQ(b, orig_a);

  // Self-swap safety.
  BigUInt<256>::CtSwap(true, a, a);
  EXPECT_EQ(a, orig_b);
}

TEST(BigUIntTest, FreeModularHelpers) {
  const BigUInt<256> even_mod =
      BigUInt<256>::FromHex("0x100000000000000000000000000000000");  // 2^128
  const BigUInt<256> a(17);
  const BigUInt<256> inv = ModInv(a, even_mod);
  EXPECT_EQ(ModMul(a, inv, even_mod), BigUInt<256>::One());
  EXPECT_EQ(ModPow(BigUInt<256>(3), BigUInt<256>(20), even_mod),
            BigUInt<256>(3486784401ULL));

  EXPECT_THROW((void)ModInv(BigUInt<256>(2), even_mod), std::domain_error);
  EXPECT_THROW((void)ModInv(BigUInt<256>(1), BigUInt<256>::One()),
               std::domain_error);
  EXPECT_THROW((void)ModMul(a, a, BigUInt<256>::Zero()), std::domain_error);
}

// ============================================================================
// MontgomeryCtx<256> and MontgomeryCtx<384> cryptographic prime tests
// ============================================================================

TEST(MontgomeryTest, InvalidModuliThrow) {
  EXPECT_THROW(static_cast<void>(MontgomeryCtx<256>(BigUInt<256>::Zero())),
               std::invalid_argument);
  EXPECT_THROW(static_cast<void>(MontgomeryCtx<256>(BigUInt<256>::One())),
               std::invalid_argument);
  EXPECT_THROW(static_cast<void>(MontgomeryCtx<256>(BigUInt<256>(2))),
               std::invalid_argument);
  EXPECT_THROW(static_cast<void>(MontgomeryCtx<256>(BigUInt<256>(1000))),
               std::invalid_argument);
}

TEST(MontgomeryTest, Curve25519Prime) {
  // p = 2^255 - 19
  const auto p = BigUInt<256>::FromHex(
      "0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed");
  VerifyMontgomerySuite<256>(p, 0x25519ULL);
}

TEST(MontgomeryTest, Secp256k1Prime) {
  // p = 2^256 - 2^32 - 977
  const auto p = BigUInt<256>::FromHex(
      "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f");
  VerifyMontgomerySuite<256>(p, 0x256a1ULL);
}

TEST(MontgomeryTest, NistP256Prime) {
  // p = 2^256 - 2^224 + 2^192 + 2^96 - 1
  const auto p = BigUInt<256>::FromHex(
      "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
  VerifyMontgomerySuite<256>(p, 0x256b2ULL);
}

TEST(MontgomeryTest, Bls12_381BasePrime) {
  // BLS12-381 base field prime (381 bits).
  const auto p = BigUInt<384>::FromHex(
      "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabff"
      "feb153ffffb9feffffffffaaab");
  VerifyMontgomerySuite<384>(p, 0x381c3ULL);
}

TEST(MontgomeryTest, CompositeOddModulus) {
  // Product of two 64-bit primes (128-bit odd composite modulus).
  const auto p1 = BigUInt<128>(18446744073709551557ULL);  // 2^64 - 59
  const auto p2 = BigUInt<128>(18446744073709551533ULL);  // 2^64 - 83
  const auto n = p1 * p2;
  MontgomeryCtx<128> ctx(n);

  const BigUInt<128> a(123456789012345ULL);
  const BigUInt<128> b(987654321098765ULL);
  EXPECT_EQ(ctx.FromMont(ctx.Mul(ctx.ToMont(a), ctx.ToMont(b))),
            ModMul(a, b, n));
  EXPECT_EQ(ctx.Mul(ctx.ToMont(a), ctx.Inv(ctx.ToMont(a))),
            ctx.ToMont(BigUInt<128>::One()));
  EXPECT_THROW((void)ctx.Inv(ctx.ToMont(p1)), std::domain_error);
}

// ============================================================================
// BigInt tests
// ============================================================================

TEST(BigIntTest, ConstructionConversionsAndFormatting) {
  EXPECT_TRUE(BigInt().IsZero());
  EXPECT_EQ(BigInt().Sign(), 0);
  EXPECT_EQ(BigInt(0).ToString(), "0");
  EXPECT_EQ(BigInt::FromDecimal("-0"), BigInt(0));
  EXPECT_FALSE(BigInt::FromDecimal("-0").IsNegative());

  const int64_t kMin64 = std::numeric_limits<int64_t>::min();
  const int64_t kMax64 = std::numeric_limits<int64_t>::max();
  const uint64_t kMaxU64 = std::numeric_limits<uint64_t>::max();

  EXPECT_EQ(BigInt(kMin64).ToInt64(), kMin64);
  EXPECT_EQ(BigInt(kMin64).ToDecimal(), "-9223372036854775808");
  EXPECT_EQ(BigInt(kMax64).ToInt64(), kMax64);
  EXPECT_EQ(BigInt(kMaxU64).ToUint64(), kMaxU64);

  EXPECT_THROW((void)(BigInt(kMax64) + BigInt(1)).ToInt64(),
               std::overflow_error);
  EXPECT_THROW((void)(BigInt(kMin64) - BigInt(1)).ToInt64(),
               std::overflow_error);
  EXPECT_THROW((void)BigInt(-1).ToUint64(), std::overflow_error);
  EXPECT_THROW((void)(BigInt(kMaxU64) + BigInt(1)).ToUint64(),
               std::overflow_error);

  // Hex round-trip with negative values and prefix.
  const BigInt neg_hex = BigInt::FromHex("-0xdeadbeef0123456789abcdef");
  EXPECT_TRUE(neg_hex.IsNegative());
  EXPECT_EQ(neg_hex.ToHex(), "-deadbeef0123456789abcdef");
  EXPECT_EQ(neg_hex.ToHex(true), "-0xdeadbeef0123456789abcdef");
  EXPECT_EQ(BigInt::FromDecimal(neg_hex.ToDecimal()), neg_hex);

  // BigUInt <-> BigInt conversion.
  const auto u = BigUInt<256>::FromHex(
      "0x123456789abcdef0fedcba987654321000112233445566778899aabbccddeeff");
  const BigInt bu = BigInt::FromBigUInt(u);
  EXPECT_EQ(bu.ToBigUInt<256>(), u);
  EXPECT_THROW((void)(-bu).ToBigUInt<256>(), std::overflow_error);
  EXPECT_THROW((void)(BigInt(1) << 256).ToBigUInt<256>(), std::overflow_error);

  // Stream operator and std::hash.
  std::ostringstream oss;
  oss << BigInt(-123456789);
  EXPECT_EQ(oss.str(), "-123456789");

  std::unordered_set<BigInt> set;
  set.insert(BigInt::FromDecimal("123456789012345678901234567890"));
  set.insert(BigInt::FromDecimal("-123456789012345678901234567890"));
  set.insert(BigInt(0));
  EXPECT_EQ(set.size(), 3U);
  EXPECT_TRUE(
      set.contains(BigInt::FromDecimal("123456789012345678901234567890")));
}

TEST(BigIntTest, BitwiseAndShiftSemantics) {
  EXPECT_EQ(BigInt(1) << 130,
            BigInt::FromHex("0x400000000000000000000000000000000"));
  EXPECT_EQ(BigInt(-1) << 64, -BigInt::FromHex("0x10000000000000000"));
  EXPECT_EQ(BigInt(-5) >> 1, BigInt(-3));
  EXPECT_EQ(BigInt(-4) >> 1, BigInt(-2));
  EXPECT_EQ(BigInt(-1) >> 200, BigInt(-1));
  EXPECT_EQ(BigInt(5) & BigInt(3), BigInt(1));
  EXPECT_EQ(BigInt(-5) & BigInt(3), BigInt(3));
  EXPECT_EQ(BigInt(-5) | BigInt(3), BigInt(-5));
  EXPECT_EQ(BigInt(-1) ^ BigInt(42), BigInt(-43));
  EXPECT_EQ(~BigInt(0), BigInt(-1));
  EXPECT_EQ(~BigInt(-1), BigInt(0));
}

TEST(BigIntTest, KaratsubaVsSchoolbookLargeOperands) {
  uint64_t rng = 0xfeedfacecafebeefULL;

  // Directly verify Karatsuba matches schoolbook on 32, 50, 100 limbs and
  // unbalanced operand shapes.
  const std::vector<std::pair<std::size_t, std::size_t>> shapes = {
      {32, 32}, {50, 50}, {64, 64}, {100, 100}, {100, 40}, {40, 100}, {75, 33}};
  for (auto [na, nb] : shapes) {
    const BigInt a = RandomBigIntLimbs(na, rng);
    const BigInt b = RandomBigIntLimbs(nb, rng);
    const auto sb = xg::detail::MulLimbsSchoolbook(a.RawLimbs(), b.RawLimbs());
    const auto ka = xg::detail::MulLimbsKaratsuba(a.RawLimbs(), b.RawLimbs());
    EXPECT_EQ(ka, sb) << "Mismatch for shape (" << na << ", " << nb << ")";

    // Sign rules on large Karatsuba products.
    const BigInt pos = a * b;
    EXPECT_EQ((-a) * b, -pos);
    EXPECT_EQ(a * (-b), -pos);
    EXPECT_EQ((-a) * (-b), pos);
  }

  // Exact 100! value.
  const BigInt f100 = Factorial(100);
  EXPECT_EQ(
      f100.ToDecimal(),
      "933262154439441526816992388562667004907159682643816214685929638952175999"
      "932299156089414639761565182862536979208272237582511852109168640000000000"
      "00000000000000");

  // (2^2048 - 1)^2 == 2^4096 - 2^2049 + 1 (32 limbs * 32 limbs -> Karatsuba).
  const BigInt m2048 = (BigInt(1) << 2048) - BigInt(1);
  const BigInt expected = (BigInt(1) << 4096) - (BigInt(1) << 2049) + BigInt(1);
  EXPECT_EQ(m2048 * m2048, expected);
  EXPECT_EQ((-m2048) * m2048, -expected);
}

TEST(BigIntTest, KnuthAlgorithmDAndEuclideanDivRem) {
  // All 4 sign combinations for DivRem (truncated) vs EuclideanDivRem.
  auto check_signs = [](const BigInt& a, const BigInt& b) {
    auto [qt, rt] = DivRem(a, b);
    EXPECT_EQ(qt * b + rt, a);
    EXPECT_LT(rt.Abs(), b.Abs());
    if (!rt.IsZero()) {
      EXPECT_EQ(rt.Sign(), a.Sign());
    }

    auto [qe, re] = EuclideanDivRem(a, b);
    EXPECT_EQ(qe * b + re, a);
    EXPECT_GE(re, BigInt(0));
    EXPECT_LT(re, b.Abs());
  };

  check_signs(BigInt(17), BigInt(5));
  check_signs(BigInt(-17), BigInt(5));
  check_signs(BigInt(17), BigInt(-5));
  check_signs(BigInt(-17), BigInt(-5));
  check_signs(BigInt(-20), BigInt(5));
  check_signs(BigInt(-20), BigInt(-5));

  // Targeted Knuth Algorithm D Step D6 (add-back) trigger:
  // v = 2^63 * B^2 + 0 * B + (B - 1), already normalized (s = 0).
  // u = (2^63 - 1) * B^3 + 2^63 * B^2 + 0 * B + 0.
  // Then qhat = B - 1, rhat = 0, qhat * v1 = 0 <= 0 passes D3, yet qhat * v > u
  // by (B - 1)^2, forcing negative remainder in D4 and add-back in D6!
  const BigInt v_d6 = (BigInt(1) << 191) +
                      BigInt::FromUint64(std::numeric_limits<uint64_t>::max());
  const BigInt u_d6 =
      (BigInt::FromUint64((uint64_t{1} << 63) - 1) << 192) + (BigInt(1) << 191);
  {
    auto [q, r] = DivRem(u_d6, v_d6);
    EXPECT_EQ(q, BigInt::FromUint64(std::numeric_limits<uint64_t>::max() - 1));
    EXPECT_EQ(q * v_d6 + r, u_d6);
    EXPECT_GE(r, BigInt(0));
    EXPECT_LT(r, v_d6);
  }

  // Trigger D6 with non-zero normalization shift s (divide both by 2^k).
  for (std::size_t shift : {1U, 7U, 31U, 63U}) {
    const BigInt v_s = v_d6 >> shift;
    const BigInt u_s = u_d6 >> shift;
    auto [q, r] = DivRem(u_s, v_s);
    EXPECT_EQ(q * v_s + r, u_s);
    EXPECT_GE(r, BigInt(0));
    EXPECT_LT(r, v_s);
  }

  // Randomized multi-limb stress test for DivRem & EuclideanDivRem.
  uint64_t rng = 0xd6d6d6d612345678ULL;
  for (int i = 0; i < 60; ++i) {
    const std::size_t nb = 1 + (NextRng(rng) % 12);
    const std::size_t na = nb + (NextRng(rng) % 15);
    const BigInt a = RandomBigIntLimbs(na, rng);
    const BigInt b = RandomBigIntLimbs(nb, rng);
    check_signs(a, b);
    check_signs(-a, b);
    check_signs(a, -b);
    check_signs(-a, -b);

    // Exact product + remainder identity.
    const BigInt q_true = RandomBigIntLimbs(1 + (NextRng(rng) % 10), rng);
    const BigInt r_true = RandomBigIntLimbs(nb, rng) % b;
    const BigInt dividend = q_true * b + r_true;
    auto [q_got, r_got] = DivRem(dividend, b);
    EXPECT_EQ(q_got, q_true);
    EXPECT_EQ(r_got, r_true);
  }

  EXPECT_THROW((void)DivRem(BigInt(1), BigInt(0)), std::domain_error);
}

TEST(BigIntTest, GcdExtendedGcdModInversePowModISqrt) {
  // ExtendedGcd across signs and large integers.
  auto verify_egcd = [](const BigInt& a, const BigInt& b) {
    auto [g, x, y] = ExtendedGcd(a, b);
    EXPECT_GE(g, BigInt(0));
    EXPECT_EQ(g, Gcd(a, b));
    EXPECT_EQ(a * x + b * y, g);
  };

  verify_egcd(BigInt(0), BigInt(0));
  verify_egcd(BigInt(-42), BigInt(0));
  verify_egcd(BigInt(0), BigInt(-42));
  verify_egcd(BigInt(240), BigInt(46));
  verify_egcd(BigInt(-240), BigInt(46));
  verify_egcd(BigInt(240), BigInt(-46));
  verify_egcd(BigInt(-240), BigInt(-46));

  const BigInt big_a = Factorial(40) + BigInt(17);
  const BigInt big_b = Factorial(35) + BigInt(29);
  verify_egcd(big_a, big_b);
  verify_egcd(-big_a, big_b);

  // ModInverse.
  const BigInt mod = (BigInt(1) << 127) - BigInt(1);  // M127 prime
  const BigInt val = BigInt::FromDecimal("98765432109876543210987654321");
  const BigInt inv = ModInverse(val, mod);
  EXPECT_EQ((val * inv) % mod, BigInt(1));
  EXPECT_EQ(ModInverse(-val, mod), mod - inv);
  EXPECT_THROW((void)ModInverse(BigInt(6), BigInt(9)), std::domain_error);
  EXPECT_THROW((void)ModInverse(BigInt(1), BigInt(1)), std::domain_error);

  // PowMod including negative exponent and Fermat / RSA round-trip.
  EXPECT_EQ(PowMod(val, BigInt(-1), mod), inv);
  EXPECT_EQ(PowMod(val, mod - BigInt(1), mod), BigInt(1));
  EXPECT_EQ(PowMod(BigInt(-3), BigInt(3), BigInt(7)), BigInt(1));
  EXPECT_EQ(PowMod(BigInt(123), BigInt(456), BigInt(1)), BigInt(0));
  EXPECT_THROW((void)PowMod(BigInt(2), BigInt(3), BigInt(0)),
               std::domain_error);

  // ISqrt on small and large values (k^2 - 1, k^2, k^2 + 1).
  EXPECT_EQ(ISqrt(BigInt(0)), BigInt(0));
  EXPECT_EQ(ISqrt(BigInt(1)), BigInt(1));
  EXPECT_EQ(ISqrt(BigInt(2)), BigInt(1));
  EXPECT_EQ(ISqrt(BigInt(3)), BigInt(1));
  EXPECT_EQ(ISqrt(BigInt(4)), BigInt(2));
  EXPECT_THROW((void)ISqrt(BigInt(-1)), std::domain_error);

  const BigInt k = (BigInt(1) << 300) + BigInt::FromDecimal("123456789012345");
  const BigInt k_sq = k * k;
  EXPECT_EQ(ISqrt(k_sq - BigInt(1)), k - BigInt(1));
  EXPECT_EQ(ISqrt(k_sq), k);
  EXPECT_EQ(ISqrt(k_sq + BigInt(1)), k);
  EXPECT_EQ(ISqrt(k_sq + (k << 1)), k);
  EXPECT_EQ(ISqrt(k_sq + (k << 1) + BigInt(1)), k + BigInt(1));
}

TEST(BigIntTest, JacobiSymbolProperties) {
  // Invalid moduli.
  EXPECT_THROW((void)JacobiSymbol(BigInt(2), BigInt(0)), std::domain_error);
  EXPECT_THROW((void)JacobiSymbol(BigInt(2), BigInt(4)), std::domain_error);
  EXPECT_THROW((void)JacobiSymbol(BigInt(2), BigInt(-7)), std::domain_error);

  // Verify against Euler's criterion for odd primes p in [3..47].
  for (int p = 3; p <= 47; p += 2) {
    if (!IsPrimeTrialSmall(static_cast<uint64_t>(p))) {
      continue;
    }
    const BigInt bp(p);
    for (int a = -p; a <= 2 * p; ++a) {
      const int j = JacobiSymbol(BigInt(a), bp);
      if (a % p == 0) {
        EXPECT_EQ(j, 0);
      } else {
        const BigInt e = PowMod(BigInt(a), BigInt((p - 1) / 2), bp);
        const int expected = e.IsOne() ? 1 : -1;
        EXPECT_EQ(j, expected) << "a=" << a << ", p=" << p;
      }
    }
  }

  // Multiplicativity on odd composites: (a / (m*n)) == (a/m) * (a/n).
  const BigInt m(1000003);
  const BigInt n(1000033);
  const BigInt mn = m * n;
  for (int a : {-17, -1, 0, 2, 5, 19, 1000003}) {
    EXPECT_EQ(JacobiSymbol(BigInt(a), mn),
              JacobiSymbol(BigInt(a), m) * JacobiSymbol(BigInt(a), n));
  }
}

TEST(BigIntTest, IsProbablePrimePrimesCarmichaelsAndMersennes) {
  // Exact verification for 0..1000.
  for (uint64_t n = 0; n <= 1000; ++n) {
    EXPECT_EQ(IsProbablePrime(BigInt(n)), IsPrimeTrialSmall(n)) << "n=" << n;
  }

  // Negative numbers are not prime.
  EXPECT_FALSE(IsProbablePrime(BigInt(-7)));

  // Carmichael numbers must ALL be rejected.
  const std::vector<uint64_t> carmichaels = {
      561, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 341550071728321ULL};
  for (uint64_t c : carmichaels) {
    EXPECT_FALSE(IsProbablePrime(BigInt::FromUint64(c))) << "c=" << c;
  }

  // Large (75-bit > 64-bit) Carmichael number: Chernick's U_3(k) for k=3000251:
  // (6k+1)(12k+1)(18k+1) = 18001507 * 36003013 * 54004519.
  const BigInt p_c1(18001507ULL);
  const BigInt p_c2(36003013ULL);
  const BigInt p_c3(54004519ULL);
  ASSERT_TRUE(IsProbablePrime(p_c1));
  ASSERT_TRUE(IsProbablePrime(p_c2));
  ASSERT_TRUE(IsProbablePrime(p_c3));
  const BigInt large_carmichael = p_c1 * p_c2 * p_c3;
  ASSERT_GT(large_carmichael.BitLength(), 64U);
  EXPECT_FALSE(IsProbablePrime(large_carmichael));

  // Perfect square of large prime (>64-bit) must be rejected.
  const BigInt m61 = (BigInt(1) << 61) - BigInt(1);
  EXPECT_FALSE(IsProbablePrime(m61 * m61));

  // Mersenne primes: 2^61 - 1, 2^89 - 1, 2^107 - 1, 2^127 - 1.
  for (std::size_t p : {61U, 89U, 107U, 127U}) {
    const BigInt mp = (BigInt(1) << p) - BigInt(1);
    EXPECT_TRUE(IsProbablePrime(mp)) << "Mersenne 2^" << p << "-1";
  }

  // Composite Mersenne numbers: 2^67 - 1, 2^101 - 1, 2^128 - 1.
  for (std::size_t p : {67U, 101U, 128U}) {
    const BigInt mc = (BigInt(1) << p) - BigInt(1);
    EXPECT_FALSE(IsProbablePrime(mc)) << "Composite 2^" << p << "-1";
  }

  // Cryptographic primes (>250 bits).
  EXPECT_TRUE(IsProbablePrime((BigInt(1) << 255) - BigInt(19)));
  EXPECT_TRUE(
      IsProbablePrime((BigInt(1) << 256) - (BigInt(1) << 32) - BigInt(977)));
  EXPECT_TRUE(IsProbablePrime(BigInt::FromHex(
      "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabff"
      "feb153ffffb9feffffffffaaab")));
}

TEST(BigIntTest, FactorIntegerSmall64BitAndMediumComposites) {
  auto verify_factorization = [](const BigInt& n) {
    const auto factors = FactorInteger(n);
    BigInt product(1);
    for (std::size_t i = 0; i < factors.size(); ++i) {
      const auto& [p, e] = factors[i];
      EXPECT_TRUE(IsProbablePrime(p)) << "Factor not prime: " << p;
      EXPECT_GE(e, 1ULL);
      if (i > 0) {
        EXPECT_LT(factors[i - 1].first, p);
      }
      product *= p.Pow(e);
    }
    EXPECT_EQ(product, n);
    return factors;
  };

  EXPECT_THROW((void)FactorInteger(BigInt(0)), std::domain_error);
  EXPECT_THROW((void)FactorInteger(BigInt(-10)), std::domain_error);
  EXPECT_TRUE(FactorInteger(BigInt(1)).empty());

  // Small numbers and prime powers.
  verify_factorization(BigInt(2));
  verify_factorization(BigInt(360));
  verify_factorization(BigInt(1) << 100);
  verify_factorization(BigInt(3).Pow(40) * BigInt(7).Pow(12));

  // 64-bit semiprimes.
  const BigInt sp64 = BigInt(4294967291ULL) * BigInt(4294967279ULL);
  const auto f_sp64 = verify_factorization(sp64);
  ASSERT_EQ(f_sp64.size(), 2U);
  EXPECT_EQ(f_sp64[0].first, BigInt(4294967279ULL));
  EXPECT_EQ(f_sp64[1].first, BigInt(4294967291ULL));

  // Fermat number F6 = 2^64 + 1 = 274177 * 67280421310721.
  const BigInt f6 = (BigInt(1) << 64) + BigInt(1);
  const auto f_f6 = verify_factorization(f6);
  ASSERT_EQ(f_f6.size(), 2U);
  EXPECT_EQ(f_f6[0].first, BigInt(274177ULL));
  EXPECT_EQ(f_f6[0].second, 1ULL);
  EXPECT_EQ(f_f6[1].first, BigInt(67280421310721ULL));
  EXPECT_EQ(f_f6[1].second, 1ULL);

  // Perfect square of a 48-bit prime.
  const BigInt p48_a = BigInt::FromUint64(281474976710597ULL);  // 2^48 - 59
  ASSERT_TRUE(IsProbablePrime(p48_a));
  const auto f_sq = verify_factorization(p48_a * p48_a);
  ASSERT_EQ(f_sq.size(), 1U);
  EXPECT_EQ(f_sq[0].first, p48_a);
  EXPECT_EQ(f_sq[0].second, 2ULL);

  // 96-bit semiprime (product of 44-bit and 52-bit primes).
  const BigInt p44 = BigInt::FromUint64(17592186044399ULL);    // 2^44 - 17
  const BigInt p52 = BigInt::FromUint64(4503599627370449ULL);  // 2^52 - 47
  ASSERT_TRUE(IsProbablePrime(p44));
  ASSERT_TRUE(IsProbablePrime(p52));
  const BigInt semi96 = p44 * p52;
  ASSERT_GE(semi96.BitLength(), 95U);
  const auto f_96 = verify_factorization(semi96);
  ASSERT_EQ(f_96.size(), 2U);
  EXPECT_EQ(f_96[0].first, p44);
  EXPECT_EQ(f_96[1].first, p52);

  // 104-bit semiprime (product of 40-bit and 64-bit primes).
  const BigInt p40_a = BigInt::FromUint64(1099511627689ULL);       // 2^40 - 87
  const BigInt p64 = BigInt::FromUint64(18446744073709551557ULL);  // 2^64 - 59
  ASSERT_TRUE(IsProbablePrime(p40_a));
  ASSERT_TRUE(IsProbablePrime(p64));
  const BigInt semi104 = p40_a * p64;
  ASSERT_GE(semi104.BitLength(), 103U);
  const auto f_104 = verify_factorization(semi104);
  ASSERT_EQ(f_104.size(), 2U);
  EXPECT_EQ(f_104[0].first, p40_a);
  EXPECT_EQ(f_104[1].first, p64);

  // 112-bit composite (product of 32-bit, 40-bit, and 40-bit primes).
  const BigInt p32 = BigInt::FromUint64(4294967291ULL);       // 2^32 - 5
  const BigInt p40_b = BigInt::FromUint64(1099511627609ULL);  // 2^40 - 167
  ASSERT_TRUE(IsProbablePrime(p32));
  ASSERT_TRUE(IsProbablePrime(p40_b));
  const BigInt comp112 = p32 * p40_a * p40_b;
  ASSERT_GE(comp112.BitLength(), 110U);
  const auto f_112 = verify_factorization(comp112);
  ASSERT_EQ(f_112.size(), 3U);
  EXPECT_EQ(f_112[0].first, p32);
  EXPECT_EQ(f_112[1].first, p40_b);
  EXPECT_EQ(f_112[2].first, p40_a);
}
