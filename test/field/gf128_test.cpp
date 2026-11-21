#include "xgalois/field/gf128.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"

namespace {

using xg::FiniteField;
using xg::GaloisFieldBinary128;
using xg::GaloisFieldElement;
using xg::GF128;
using xg::GFArray;
using xg::UInt128;

constexpr UInt128 U128(uint64_t lo, uint64_t hi = 0) noexcept {
  return UInt128{lo, hi};
}

// Independent reference NIST SP 800-38D Algorithm 1 (block multiplication)
// to cross-check GhashMul on arbitrary random 128-bit blocks.
UInt128 ReferenceNistGhashMul(UInt128 x, UInt128 y) {
  const UInt128 kR{0ULL, 0xe100000000000000ULL};
  UInt128 z{0, 0};
  UInt128 v = y;
  for (int i = 0; i < 128; ++i) {
    if (x.Bit(127 - i)) {
      z ^= v;
    }
    const bool lsb = (v.lo & 1ULL) != 0;
    v >>= 1;
    if (lsb) {
      v ^= kR;
    }
  }
  return z;
}

TEST(UInt128Test, BasicPropertiesAndBitwiseOps) {
  static_assert(std::is_trivially_copyable_v<UInt128>);
  static_assert(std::is_standard_layout_v<UInt128>);
  static_assert(sizeof(UInt128) == 16);

  const UInt128 zero{};
  const UInt128 one{1ULL, 0ULL};
  const UInt128 high_bit{0ULL, 1ULL << 63};

  EXPECT_TRUE(zero.IsZero());
  EXPECT_FALSE(zero.IsOne());
  EXPECT_EQ(zero.BitWidth(), 0);
  EXPECT_EQ(zero.PopCount(), 0);

  EXPECT_FALSE(one.IsZero());
  EXPECT_TRUE(one.IsOne());
  EXPECT_EQ(one.BitWidth(), 1);
  EXPECT_TRUE(one.Bit(0));
  EXPECT_FALSE(one.Bit(1));
  EXPECT_FALSE(one.Bit(-1));
  EXPECT_FALSE(one.Bit(128));

  EXPECT_EQ(high_bit.BitWidth(), 128);
  EXPECT_TRUE(high_bit.Bit(127));
  EXPECT_FALSE(high_bit.Bit(126));

  const UInt128 a{0x0123456789abcdefULL, 0xfedcba9876543210ULL};
  const UInt128 b{0xf0f0f0f0f0f0f0f0ULL, 0x0f0f0f0f0f0f0f0fULL};
  EXPECT_EQ((a ^ b), U128(a.lo ^ b.lo, a.hi ^ b.hi));
  EXPECT_EQ((a & b), U128(a.lo & b.lo, a.hi & b.hi));
  EXPECT_EQ((a | b), U128(a.lo | b.lo, a.hi | b.hi));
  EXPECT_EQ(~a, U128(~a.lo, ~a.hi));

  // Boundary shifts (must have zero UB for 0, 63, 64, 127, 128).
  EXPECT_EQ(a << 0, a);
  EXPECT_EQ(a >> 0, a);
  EXPECT_EQ(a << 128, zero);
  EXPECT_EQ(a >> 128, zero);
  EXPECT_EQ(one << 63, U128(1ULL << 63, 0ULL));
  EXPECT_EQ(one << 64, U128(0ULL, 1ULL));
  EXPECT_EQ(one << 127, high_bit);
  EXPECT_EQ(high_bit >> 64, U128(1ULL << 63, 0ULL));
  EXPECT_EQ(high_bit >> 127, one);
  EXPECT_EQ((a << 13) >> 13,
            U128(a.lo, a.hi & ((uint64_t{1} << (64 - 13)) - 1ULL)));

  // Ordering compares hi before lo.
  EXPECT_LT(U128(~uint64_t{0}, 0ULL), U128(0ULL, 1ULL));
  EXPECT_GT(U128(0ULL, 1ULL), U128(~uint64_t{0}, 0ULL));
  EXPECT_EQ(a <=> a, std::strong_ordering::equal);
  EXPECT_NE(a, b);

  // Bit reversal is an involution.
  EXPECT_EQ(one.BitReverse(), high_bit);
  EXPECT_EQ(a.BitReverse().BitReverse(), a);
}

TEST(UInt128Test, HexDecimalBytesAndHashRoundTrips) {
  const UInt128 val = UInt128::FromHex("0x66e94bd4ef8a2c3b884cfa59ca342b2e");
  EXPECT_EQ(val.hi, 0x66e94bd4ef8a2c3bULL);
  EXPECT_EQ(val.lo, 0x884cfa59ca342b2eULL);
  EXPECT_EQ(val.ToHex(), "0x66e94bd4ef8a2c3b884cfa59ca342b2e");
  EXPECT_EQ(UInt128::FromHex("66E94BD4EF8A2C3B884CFA59CA342B2E"), val);

  const UInt128 small{0x87ULL, 0ULL};
  EXPECT_EQ(small.ToHex(), "0x87");
  EXPECT_EQ(small.ToHex(/*zero_pad=*/true),
            "0x00000000000000000000000000000087");
  EXPECT_EQ(UInt128::FromHex(small.ToHex()), small);
  EXPECT_EQ(UInt128::FromHex(small.ToHex(/*zero_pad=*/true)), small);
  EXPECT_EQ(U128(0, 0).ToHex(), "0x0");
  EXPECT_EQ(U128(1, 1).ToHex(), "0x10000000000000001");

  // Decimal round-trip including 0 and 2^128 - 1.
  const UInt128 max_val{~uint64_t{0}, ~uint64_t{0}};
  EXPECT_EQ(U128(0, 0).ToDecimal(), "0");
  EXPECT_EQ(max_val.ToDecimal(), "340282366920938463463374607431768211455");
  EXPECT_EQ(UInt128::FromDecimal("340282366920938463463374607431768211455"),
            max_val);
  EXPECT_EQ(UInt128::FromDecimal(val.ToDecimal()), val);
  EXPECT_THROW(UInt128::FromDecimal("340282366920938463463374607431768211456"),
               std::out_of_range);
  EXPECT_THROW(UInt128::FromHex("0x100000000000000000000000000000000"),
               std::out_of_range);
  EXPECT_THROW(UInt128::FromHex("0x"), std::invalid_argument);
  EXPECT_THROW(UInt128::FromHex("0xzzzz"), std::invalid_argument);

  // Big-endian byte conversion.
  const auto bytes = val.ToBytesBigEndian();
  EXPECT_EQ(bytes[0], 0x66);
  EXPECT_EQ(bytes[1], 0xe9);
  EXPECT_EQ(bytes[15], 0x2e);
  EXPECT_EQ(UInt128::FromBytesBigEndian(bytes), val);
  const std::vector<uint8_t> vec(bytes.begin(), bytes.end());
  EXPECT_EQ(UInt128::FromBytesBigEndian(std::span<const uint8_t>(vec)), val);
  const std::vector<uint8_t> bad_len(15, 0);
  EXPECT_THROW(UInt128::FromBytesBigEndian(std::span<const uint8_t>(bad_len)),
               std::invalid_argument);

  // Stream output and std::hash.
  std::ostringstream oss;
  oss << val;
  EXPECT_EQ(oss.str(), val.ToHex());

  std::unordered_set<UInt128> set;
  set.insert(UInt128{0, 0});
  set.insert(small);
  set.insert(val);
  set.insert(max_val);
  EXPECT_EQ(set.size(), 4U);
  EXPECT_TRUE(set.contains(val));
}

TEST(GaloisFieldBinary128Test, ConceptAndMetadata) {
  static_assert(FiniteField<GaloisFieldBinary128>);
  static_assert(FiniteField<GF128>);

  GF128 gf;
  EXPECT_EQ(gf.Characteristic(), 2ULL);
  EXPECT_EQ(gf.Order(), 0ULL);
  EXPECT_EQ(gf.Degree(), 128U);
  EXPECT_EQ(gf.ModulusLow(), GF128::kGhashModulusLow);
  EXPECT_EQ(gf.AdditiveIdentity(), U128(0, 0));
  EXPECT_EQ(gf.MultiplicativeIdentity(), U128(1, 0));
  EXPECT_TRUE(gf.IsZero(U128(0, 0)));
  EXPECT_TRUE(gf.IsOne(U128(1, 0)));
  EXPECT_EQ(gf.MultiplicativeGenerator(), U128(2, 0));
  EXPECT_EQ(gf.MultiplicativeOrder128(U128(2, 0)),
            U128(~uint64_t{0}, ~uint64_t{0}));
  EXPECT_THROW(gf.MultiplicativeGenerators(), std::length_error);
  EXPECT_GE(gf.MemoryBytes(), sizeof(GF128));

  // Verify copy construction and representation atomics.
  gf.SetRepresentation(xg::FieldRepresentation::POLY);
  GF128 gf_copy = gf;
  EXPECT_EQ(gf_copy.GetRepresentation(), xg::FieldRepresentation::POLY);

  std::ostringstream oss;
  oss << gf;
  EXPECT_EQ(oss.str(), "GF(2^128)");
}

TEST(GaloisFieldBinary128Test, FieldAxiomsAndAlgebraicIdentities) {
  GF128 gf;
  std::mt19937_64 rng(0x123456789abcdef0ULL);
  const UInt128 zero = gf.AdditiveIdentity();
  const UInt128 one = gf.MultiplicativeIdentity();
  const UInt128 group_order{~uint64_t{0}, ~uint64_t{0}};

  // Edge cases: 0 and 1.
  EXPECT_EQ(gf.Add(zero, one), one);
  EXPECT_EQ(gf.Add(one, one), zero);
  EXPECT_EQ(gf.Mul(zero, one), zero);
  EXPECT_EQ(gf.Mul(one, one), one);
  EXPECT_EQ(gf.Inv(one), one);
  EXPECT_THROW(gf.Inv(zero), std::domain_error);
  EXPECT_THROW(gf.Div(one, zero), std::domain_error);
  EXPECT_EQ(gf.Pow(zero, 0), one);
  EXPECT_EQ(gf.Pow(zero, 1), zero);
  EXPECT_EQ(gf.Pow(zero, 42), zero);
  EXPECT_EQ(gf.Pow128(zero, zero), one);
  EXPECT_EQ(gf.Pow128(zero, one), zero);
  EXPECT_EQ(gf.Sqrt(zero), zero);
  EXPECT_EQ(gf.Sqrt(one), one);
  EXPECT_EQ(gf.Norm(zero), zero);
  EXPECT_EQ(gf.Norm(one), one);
  EXPECT_EQ(gf.Trace(zero), zero);
  EXPECT_EQ(gf.Trace(one), zero);  // Tr_{GF(2^128)/GF(2)}(1) = 128 mod 2 = 0.
  EXPECT_EQ(gf.Log(one), 0ULL);
  EXPECT_EQ(gf.Log(gf.MultiplicativeGenerator()), 1ULL);
  EXPECT_THROW(gf.Log(zero), std::domain_error);

  // Verify x^128 == x^7 + x^2 + x + 1 in polynomial basis.
  const UInt128 x{2ULL, 0ULL};
  const UInt128 x127{0ULL, 1ULL << 63};
  EXPECT_EQ(gf.Mul(x, x127), U128(0x87ULL, 0ULL));

  // Cross-check fast GHASH reduction vs general Barrett reduction on 0x87.
  const UInt128 mu_ghash =
      xg::detail::ComputeBarrettMuLow128(GF128::kGhashModulusLow);

  for (int trial = 0; trial < 32; ++trial) {
    const UInt128 a = gf.Random(rng);
    const UInt128 b = gf.Random(rng);
    const UInt128 c = gf.Random(rng);

    // Commutativity, associativity, distributivity, characteristic 2.
    EXPECT_EQ(gf.Add(a, b), gf.Add(b, a));
    EXPECT_EQ(gf.Add(a, a), zero);
    EXPECT_EQ(gf.Sub(a, b), gf.Add(a, b));
    EXPECT_EQ(gf.Neg(a), a);
    EXPECT_EQ(gf.Mul(a, b), gf.Mul(b, a));
    EXPECT_EQ(gf.Mul(gf.Mul(a, b), c), gf.Mul(a, gf.Mul(b, c)));
    EXPECT_EQ(gf.Mul(a, gf.Add(b, c)), gf.Add(gf.Mul(a, b), gf.Mul(a, c)));
    EXPECT_EQ(gf.Mul(a, one), a);
    EXPECT_EQ(gf.Mul(a, zero), zero);

    // Verify fast GHASH reduction matches Barrett reduction bit-for-bit.
    const auto [p_hi, p_lo] = xg::detail::ClMul128(a, b);
    EXPECT_EQ(gf.Mul(a, b), xg::detail::ReduceBarrett128(
                                p_hi, p_lo, GF128::kGhashModulusLow, mu_ghash));

    // Frobenius additivity: (a + b)^2 == a^2 + b^2.
    EXPECT_EQ(gf.Frobenius(gf.Add(a, b)),
              gf.Add(gf.Frobenius(a), gf.Frobenius(b)));
    EXPECT_EQ(gf.Frobenius(a), gf.Mul(a, a));

    // Trace additivity & GF(2) codomain: Tr(a + b) == Tr(a) + Tr(b), Tr(a) in
    // {0, 1}.
    const UInt128 tr_a = gf.Trace(a);
    const UInt128 tr_b = gf.Trace(b);
    EXPECT_TRUE(tr_a == zero || tr_a == one);
    EXPECT_TRUE(tr_b == zero || tr_b == one);
    EXPECT_EQ(gf.Trace(gf.Add(a, b)), gf.Add(tr_a, tr_b));
    EXPECT_EQ(gf.Trace(gf.Frobenius(a)), tr_a);

    // Sqrt: Sqrt(a)^2 == a.
    const UInt128 s = gf.Sqrt(a);
    EXPECT_EQ(gf.Mul(s, s), a);

    // Inverse, Division, Norm, Signed Powers, Fermat's Little Theorem.
    if (!a.IsZero()) {
      const UInt128 inv_a = gf.Inv(a);
      EXPECT_EQ(gf.Mul(a, inv_a), one);
      EXPECT_EQ(gf.Norm(a), one);
      EXPECT_EQ(gf.Pow128(a, group_order), one);
      EXPECT_EQ(gf.PowSigned(a, -5), gf.Inv(gf.Pow(a, 5)));
    }
    if (!b.IsZero()) {
      const UInt128 q = gf.Div(a, b);
      EXPECT_EQ(gf.Mul(q, b), a);
    }
  }
}

TEST(GaloisFieldBinary128Test, NistSp80038DGhashOfficialVectors) {
  GF128 gf;

  // NIST SP 800-38D Test Case 1:
  // K = 0^128, IV = 0^96, P = empty, A = empty.
  // H = AES_K(0^128) = 66e94bd4ef8a2c3b884cfa59ca342b2e
  // len(A) || len(C) = 0^128 -> GHASH = 0 * H = 0^128.
  const UInt128 h1 = UInt128::FromHex("0x66e94bd4ef8a2c3b884cfa59ca342b2e");
  const UInt128 ghash1 = gf.GhashMul(UInt128{0, 0}, h1);
  EXPECT_EQ(ghash1, U128(0, 0));
  // Verify GHASH multiplicative identity is 0x8000...0000 (bit 127 = t^0).
  const UInt128 ghash_one =
      UInt128::FromHex("0x80000000000000000000000000000000");
  EXPECT_EQ(gf.GhashMul(ghash_one, h1), h1);
  EXPECT_EQ(gf.GhashMul(h1, ghash_one), h1);

  // NIST SP 800-38D Test Case 2:
  // K = 0^128, IV = 0^96, P = 0^128 (16 bytes), A = empty.
  // H    = 66e94bd4ef8a2c3b884cfa59ca342b2e
  // C1   = 0388dace60b6a392f328c2b971b2fe78
  // X1   = C1 * H = 5e2ec746917062882c85b0685353deb7
  // L    = 00000000000000000000000000000080 (len(A)=0, len(C)=128 bits)
  // X2   = (X1 ^ L) * H = f38cbb1ad69223dcc3457ae5b6b0f885
  // J0   = E_K(Y0) = 58e2fccefa7e3061367f1d57a4e7455a
  // Tag  = X2 ^ J0 = ab6e47d42cec13bdf53a67b21257bddf
  const UInt128 c1 = UInt128::FromHex("0x0388dace60b6a392f328c2b971b2fe78");
  const UInt128 expected_x1 =
      UInt128::FromHex("0x5e2ec746917062882c85b0685353deb7");
  const UInt128 x1 = gf.GhashMul(c1, h1);
  EXPECT_EQ(x1, expected_x1);

  const UInt128 len_block2 =
      UInt128::FromHex("0x00000000000000000000000000000080");
  const UInt128 expected_x2 =
      UInt128::FromHex("0xf38cbb1ad69223dcc3457ae5b6b0f885");
  const UInt128 x2 = gf.GhashMul(gf.Add(x1, len_block2), h1);
  EXPECT_EQ(x2, expected_x2);

  const UInt128 ek_y0_tc2 =
      UInt128::FromHex("0x58e2fccefa7e3061367f1d57a4e7455a");
  const UInt128 expected_tag2 =
      UInt128::FromHex("0xab6e47d42cec13bdf53a67b21257bddf");
  EXPECT_EQ(gf.Add(x2, ek_y0_tc2), expected_tag2);

  // NIST SP 800-38D Test Case 3 (4 ciphertext blocks, 512 bits):
  const UInt128 h3 = UInt128::FromHex("0xb83b533708bf535d0aa6e52980d53b78");
  const std::array<UInt128, 5> tc3_blocks = {
      UInt128::FromHex("0x42831ec2217774244b7221b784d0d49c"),
      UInt128::FromHex("0xe3aa212f2c02a4e035c17e2329aca12e"),
      UInt128::FromHex("0x21d514b25466931c7d8f6a5aac84aa05"),
      UInt128::FromHex("0x1ba30b396a0aac973d58e091473f5985"),
      UInt128::FromHex("0x00000000000000000000000000000200"),
  };
  UInt128 acc{0, 0};
  for (const UInt128& blk : tc3_blocks) {
    acc = gf.GhashMul(gf.Add(acc, blk), h3);
  }
  const UInt128 expected_ghash3 =
      UInt128::FromHex("0x7f1b32b81b820d02614f8895ac1d4eac");
  EXPECT_EQ(acc, expected_ghash3);
  const UInt128 ek_y0_tc3 =
      UInt128::FromHex("0x3247184b3c4f69a44dbcd22887bbb418");
  const UInt128 expected_tag3 =
      UInt128::FromHex("0x4d5c2af327cd64a62cf35abd2ba6fab4");
  EXPECT_EQ(gf.Add(acc, ek_y0_tc3), expected_tag3);

  // Cross-check GhashMul against independent NIST Algorithm 1 on random inputs.
  std::mt19937_64 rng(0xfeedfacecafebeefULL);
  for (int i = 0; i < 32; ++i) {
    const UInt128 a = gf.Random(rng);
    const UInt128 b = gf.Random(rng);
    EXPECT_EQ(gf.GhashMul(a, b), ReferenceNistGhashMul(a, b));
  }
}

TEST(GaloisFieldBinary128Test, CustomIrreducibleAndRabinRejection) {
  // Custom irreducible polynomial: x^128 + x^29 + x^27 + x^2 + 1.
  const UInt128 custom_mod_low{(1ULL << 29) | (1ULL << 27) | (1ULL << 2) | 1ULL,
                               0ULL};
  EXPECT_TRUE(GF128::IsIrreducible(GF128::kGhashModulusLow));
  EXPECT_TRUE(GF128::IsIrreducible(custom_mod_low));

  GF128 custom_gf(custom_mod_low, /*check_irreducible=*/true);
  EXPECT_EQ(custom_gf.ModulusLow(), custom_mod_low);

  // Verify x * x^127 == custom_mod_low.
  EXPECT_EQ(custom_gf.Mul(U128(2, 0), U128(0, 1ULL << 63)), custom_mod_low);

  // Verify field axioms & generator on custom polynomial.
  const UInt128 gen = custom_gf.MultiplicativeGenerator();
  EXPECT_EQ(custom_gf.MultiplicativeOrder128(gen),
            U128(~uint64_t{0}, ~uint64_t{0}));

  std::mt19937_64 rng(0x9876543210fedcbaULL);
  for (int i = 0; i < 16; ++i) {
    const UInt128 a = custom_gf.Random(rng);
    const UInt128 b = custom_gf.Random(rng);
    const UInt128 c = custom_gf.Random(rng);
    EXPECT_EQ(custom_gf.Mul(custom_gf.Mul(a, b), c),
              custom_gf.Mul(a, custom_gf.Mul(b, c)));
    EXPECT_EQ(custom_gf.Mul(a, custom_gf.Add(b, c)),
              custom_gf.Add(custom_gf.Mul(a, b), custom_gf.Mul(a, c)));
    if (!a.IsZero()) {
      EXPECT_EQ(custom_gf.Mul(a, custom_gf.Inv(a)), U128(1, 0));
    }
    const UInt128 s = custom_gf.Sqrt(a);
    EXPECT_EQ(custom_gf.Mul(s, s), a);
  }

  // Reject reducible polynomials when check_irreducible = true:
  // 1. x^128 (zero modulus_low) is always rejected.
  EXPECT_THROW(GF128(U128(0, 0), false), std::invalid_argument);
  // 2. Divisible by x: x^128 + x^7 + x^2 + x.
  EXPECT_FALSE(GF128::IsIrreducible(U128(0x86ULL, 0ULL)));
  EXPECT_THROW(GF128(U128(0x86ULL, 0ULL), true), std::invalid_argument);
  // 3. Divisible by x + 1: x^128 + 1 = (x + 1)^128.
  EXPECT_FALSE(GF128::IsIrreducible(U128(1ULL, 0ULL)));
  EXPECT_THROW(GF128(U128(1ULL, 0ULL), true), std::invalid_argument);
  // 4. Square: x^128 + x^2 + 1 = (x^64 + x + 1)^2.
  EXPECT_FALSE(GF128::IsIrreducible(U128(5ULL, 0ULL)));
  EXPECT_THROW(GF128(U128(5ULL, 0ULL), true), std::invalid_argument);
  // 5. Non-trivial product of two degree-64 polynomials with P(0) = P(1) = 1:
  //    (x^64 + x^4 + x^3 + x + 1) * (x^64 + x^11 + x^2 + x + 1).
  const UInt128 product_reducible =
      UInt128::FromHex("0x000000000000081c000000000000d841");
  EXPECT_FALSE(GF128::IsIrreducible(product_reducible));
  EXPECT_THROW(GF128(product_reducible, true), std::invalid_argument);
}

TEST(GaloisFieldBinary128Test, StringParsingAndFormatting) {
  GF128 gf;

  // Hex round-trip.
  const UInt128 v = gf.SetElementValue("0xdeadbeef0123456789abcdef11223344");
  EXPECT_EQ(v.hi, 0xdeadbeef01234567ULL);
  EXPECT_EQ(v.lo, 0x89abcdef11223344ULL);
  gf.SetRepresentation(xg::FieldRepresentation::HEX);
  EXPECT_EQ(gf.ToString(v), "0xdeadbeef0123456789abcdef11223344");
  EXPECT_EQ(gf.SetElementValue(gf.ToString(v)), v);

  // Decimal round-trip.
  gf.SetRepresentation(xg::FieldRepresentation::INT);
  const std::string dec_str = gf.ToString(v);
  EXPECT_EQ(gf.SetElementValue(dec_str), v);
  EXPECT_EQ(gf.SetElementValue("135"), U128(0x87ULL, 0ULL));

  // Binary polynomial round-trip.
  gf.SetRepresentation(xg::FieldRepresentation::POLY);
  const UInt128 p = gf.SetElementValue("x^127 + x^7 + x^2 + x + 1");
  EXPECT_EQ(p, U128(0x87ULL, 1ULL << 63));
  EXPECT_EQ(gf.ToString(p), "x^127 + x^7 + x^2 + x + 1");
  EXPECT_EQ(gf.SetElementValue(gf.ToString(p)), p);
  EXPECT_EQ(gf.ToString(U128(0, 0)), "0");

  // Invalid inputs & unsupported representations.
  EXPECT_THROW(gf.SetElementValue(""), std::invalid_argument);
  EXPECT_THROW(gf.SetElementValue("x^128 + 1"), std::invalid_argument);
  EXPECT_THROW(gf.SetRepresentation(xg::FieldRepresentation::POW),
               std::invalid_argument);
  EXPECT_THROW(gf.SetRepresentation(xg::FieldRepresentation::LOG),
               std::invalid_argument);
}

TEST(GaloisFieldBinary128Test, BulkSpanKernelsAndBatchInv) {
  const auto test_kernels = [](const GF128& gf) {
    std::mt19937_64 rng(0x42424242ULL);
    constexpr size_t kN = 37;
    std::vector<UInt128> a(kN), b(kN), dst(kN);
    for (size_t i = 0; i < kN; ++i) {
      a[i] = gf.Random(rng);
      b[i] = gf.Random(rng);
      if (a[i].IsZero()) {
        a[i] = UInt128{1, 0};
      }
    }

    // AddArray & SubArray.
    gf.AddArray(dst, a, b);
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(dst[i], gf.Add(a[i], b[i]));
    }
    gf.SubArray(dst, a, b);
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(dst[i], gf.Sub(a[i], b[i]));
    }

    // MulArray.
    gf.MulArray(dst, a, b);
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(dst[i], gf.Mul(a[i], b[i]));
    }

    // ScalarMulArray (zero, one, and random scalar).
    const UInt128 scalar = gf.Random(rng);
    gf.ScalarMulArray(dst, a, UInt128{0, 0});
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_TRUE(dst[i].IsZero());
    }
    gf.ScalarMulArray(dst, a, UInt128{1, 0});
    EXPECT_EQ(dst, a);
    gf.ScalarMulArray(dst, a, scalar);
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(dst[i], gf.Mul(a[i], scalar));
    }

    // FmaArray (zero, one, and random scalar).
    std::vector<UInt128> acc = b;
    gf.FmaArray(acc, a, UInt128{0, 0});
    EXPECT_EQ(acc, b);
    gf.FmaArray(acc, a, UInt128{1, 0});
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(acc[i], gf.Add(b[i], a[i]));
    }
    acc = b;
    gf.FmaArray(acc, a, scalar);
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(acc[i], gf.Add(b[i], gf.Mul(a[i], scalar)));
    }

    // DotProduct (lazy reduction matches naive accumulation).
    UInt128 expected_dot{0, 0};
    for (size_t i = 0; i < kN; ++i) {
      expected_dot = gf.Add(expected_dot, gf.Mul(a[i], b[i]));
    }
    EXPECT_EQ(gf.DotProduct(a, b), expected_dot);

    // BatchInv out-of-place and in-place.
    std::vector<UInt128> invs(kN);
    gf.BatchInv(invs, a);
    for (size_t i = 0; i < kN; ++i) {
      EXPECT_EQ(gf.Mul(a[i], invs[i]), U128(1, 0));
    }
    std::vector<UInt128> inplace = a;
    gf.BatchInv(inplace, inplace);
    EXPECT_EQ(inplace, invs);

    // BatchInv edge cases: empty & zero entry.
    std::vector<UInt128> empty;
    EXPECT_NO_THROW(gf.BatchInv(empty, empty));
    std::vector<UInt128> with_zero = a;
    with_zero[kN / 2] = UInt128{0, 0};
    EXPECT_THROW(gf.BatchInv(invs, with_zero), std::domain_error);
  };

  // Run on both default GHASH polynomial and custom irreducible polynomial.
  test_kernels(GF128());
  const UInt128 custom_mod{(1ULL << 29) | (1ULL << 27) | (1ULL << 2) | 1ULL,
                           0ULL};
  test_kernels(GF128(custom_mod, /*check_irreducible=*/true));
}

TEST(GaloisFieldBinary128Test, ElementAndGFArrayIntegration) {
  auto gf = std::make_shared<GF128>();

  GaloisFieldElement<GF128> e1(UInt128::FromHex("0x123456789abcdef011223344"),
                               gf);
  GaloisFieldElement<GF128> e2("x^7 + x^2 + x + 1", gf);
  EXPECT_EQ(e2.Value(), U128(0x87ULL, 0ULL));
  EXPECT_EQ((e1 + e2) - e2, e1);
  EXPECT_EQ((e1 * e2) / e2, e1);
  EXPECT_EQ(e1 * e1.Inv(), GaloisFieldElement<GF128>(U128(1, 0), gf));
  EXPECT_EQ(e1.Sqrt() * e1.Sqrt(), e1);
  EXPECT_EQ(e1.Frobenius(), e1 * e1);
  EXPECT_EQ(e1.Norm().Value(), U128(1, 0));
  EXPECT_TRUE(e1.Trace().IsZero() || e1.Trace().IsOne());

  // GFArray<GF128> integration.
  std::mt19937_64 rng(0x5555aaaaULL);
  auto arr_a = GFArray<GF128>::Random(gf, {2, 3}, rng);
  auto arr_b = GFArray<GF128>::Random(gf, {2, 3}, rng);
  auto arr_sum = arr_a + arr_b;
  auto arr_diff = arr_sum - arr_b;
  EXPECT_EQ(arr_diff, arr_a);

  auto arr_prod = arr_a * arr_b;
  auto arr_scaled = arr_a * e2;
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(arr_prod.At(i, j), arr_a.At(i, j) * arr_b.At(i, j));
      EXPECT_EQ(arr_scaled.At(i, j), arr_a.At(i, j) * e2);
    }
  }

  auto zeros = GFArray<GF128>::Zeros(gf, {4});
  auto ones = GFArray<GF128>::Ones(gf, {4});
  EXPECT_EQ(zeros + ones, ones);
  EXPECT_EQ(ones * ones, ones);
}

}  // namespace
