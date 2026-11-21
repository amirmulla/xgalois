#include "xgalois/field/gf_prime.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace xg;

// Verify FiniteField concept compliance and UnsignedFieldWord constraints.
static_assert(FiniteField<GaloisFieldPrime<uint8_t>>);
static_assert(FiniteField<GaloisFieldPrime<uint16_t>>);
static_assert(FiniteField<GaloisFieldPrime<uint32_t>>);
static_assert(FiniteField<GaloisFieldPrime<uint64_t>>);
static_assert(FiniteField<GaloisFieldPrimeTable<uint8_t>>);
static_assert(FiniteField<GaloisFieldPrimeTable<uint16_t>>);
static_assert(FiniteField<GaloisFieldPrimeTable<uint32_t>>);
static_assert(FiniteField<GaloisFieldPrimeTable<uint64_t>>);
static_assert(FiniteField<GFP<uint32_t>>);
static_assert(FiniteField<GFPLOG<uint32_t>>);
static_assert(!UnsignedFieldWord<int>);
static_assert(!UnsignedFieldWord<int32_t>);
static_assert(!UnsignedFieldWord<int64_t>);
static_assert(!UnsignedFieldWord<bool>);

namespace {

__extension__ using u128 = unsigned __int128;

uint64_t ExactModMul(uint64_t a, uint64_t b, uint64_t p) {
  return static_cast<uint64_t>(
      (static_cast<u128>(a % p) * static_cast<u128>(b % p)) % p);
}

template <typename ElementType>
void RunExhaustiveDifferentialTest(ElementType p) {
  GaloisFieldPrime<ElementType> prime_field(p);
  GaloisFieldPrimeTable<ElementType> table_field(p);

  const uint64_t p64 = static_cast<uint64_t>(p);
  EXPECT_EQ(prime_field.Characteristic(), p64);
  EXPECT_EQ(table_field.Characteristic(), p64);
  EXPECT_EQ(prime_field.Order(), p64);
  EXPECT_EQ(table_field.Order(), p64);
  EXPECT_EQ(prime_field.Degree(), 1U);
  EXPECT_EQ(table_field.Degree(), 1U);

  const ElementType gen = prime_field.MultiplicativeGenerator();
  EXPECT_EQ(table_field.MultiplicativeGenerator(), gen);

  // Check all pairs (a, b).
  for (uint64_t a = 0; a < p64; ++a) {
    const ElementType ea = static_cast<ElementType>(a);
    EXPECT_EQ(prime_field.Neg(ea), table_field.Neg(ea));
    EXPECT_EQ(prime_field.Add(ea, prime_field.Neg(ea)), 0);
    EXPECT_EQ(prime_field.Frobenius(ea), ea);
    EXPECT_EQ(table_field.Frobenius(ea), ea);
    EXPECT_EQ(prime_field.Trace(ea), ea);
    EXPECT_EQ(table_field.Trace(ea), ea);
    EXPECT_EQ(prime_field.Norm(ea), ea);
    EXPECT_EQ(table_field.Norm(ea), ea);

    if (a == 0) {
      EXPECT_THROW(prime_field.Inv(ea), std::domain_error);
      EXPECT_THROW(table_field.Inv(ea), std::domain_error);
      EXPECT_THROW(prime_field.Log(ea), std::domain_error);
      EXPECT_THROW(table_field.Log(ea), std::domain_error);
      EXPECT_EQ(prime_field.Sqrt(ea), 0);
      EXPECT_EQ(table_field.Sqrt(ea), 0);
      EXPECT_TRUE(prime_field.IsQuadraticResidue(ea));
      EXPECT_TRUE(table_field.IsQuadraticResidue(ea));
    } else {
      const ElementType inv_p = prime_field.Inv(ea);
      const ElementType inv_t = table_field.Inv(ea);
      EXPECT_EQ(inv_p, inv_t);
      EXPECT_EQ(prime_field.Mul(ea, inv_p), 1);
      EXPECT_EQ(table_field.Mul(ea, inv_t), 1);

      const uint64_t log_p = prime_field.Log(ea);
      const uint64_t log_t = table_field.Log(ea);
      EXPECT_EQ(log_p, log_t);
      EXPECT_EQ(prime_field.Pow(gen, log_p), ea);
      EXPECT_EQ(table_field.Pow(gen, log_t), ea);

      const bool qr_p = prime_field.IsQuadraticResidue(ea);
      const bool qr_t = table_field.IsQuadraticResidue(ea);
      EXPECT_EQ(qr_p, qr_t);
      if (qr_p) {
        const ElementType r_p = prime_field.Sqrt(ea);
        const ElementType r_t = table_field.Sqrt(ea);
        EXPECT_EQ(r_p, r_t);
        EXPECT_EQ(prime_field.Mul(r_p, r_p), ea);
        EXPECT_LE(static_cast<uint64_t>(r_p), p64 - static_cast<uint64_t>(r_p));
      } else {
        EXPECT_THROW(prime_field.Sqrt(ea), std::domain_error);
        EXPECT_THROW(table_field.Sqrt(ea), std::domain_error);
      }

      const uint64_t ord_p = prime_field.MultiplicativeOrder(ea);
      const uint64_t ord_t = table_field.MultiplicativeOrder(ea);
      EXPECT_EQ(ord_p, ord_t);
      EXPECT_EQ((p64 - 1) % ord_p, 0ULL);
      EXPECT_EQ(prime_field.Pow(ea, ord_p), 1);
    }

    for (uint64_t b = 0; b < p64; ++b) {
      const ElementType eb = static_cast<ElementType>(b);
      const uint64_t ref_add = (a + b) % p64;
      const uint64_t ref_sub = (a + p64 - b) % p64;
      const uint64_t ref_mul = (a * b) % p64;

      EXPECT_EQ(prime_field.Add(ea, eb), ref_add);
      EXPECT_EQ(table_field.Add(ea, eb), ref_add);
      EXPECT_EQ(prime_field.Sub(ea, eb), ref_sub);
      EXPECT_EQ(table_field.Sub(ea, eb), ref_sub);
      EXPECT_EQ(prime_field.Mul(ea, eb), ref_mul);
      EXPECT_EQ(table_field.Mul(ea, eb), ref_mul);

      if (b == 0) {
        EXPECT_THROW(prime_field.Div(ea, eb), std::domain_error);
        EXPECT_THROW(table_field.Div(ea, eb), std::domain_error);
      } else {
        const ElementType d_p = prime_field.Div(ea, eb);
        const ElementType d_t = table_field.Div(ea, eb);
        EXPECT_EQ(d_p, d_t);
        EXPECT_EQ(prime_field.Mul(d_p, eb), ea);
      }
    }
  }
}

}  // namespace

class GaloisFieldPrimeTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf2 = std::make_unique<GaloisFieldPrime<uint8_t>>(2);
    gf3 = std::make_unique<GaloisFieldPrime<uint8_t>>(3);
    gf5 = std::make_unique<GaloisFieldPrime<uint8_t>>(5);
    gf7 = std::make_unique<GaloisFieldPrime<uint8_t>>(7);
    gf13 = std::make_unique<GaloisFieldPrime<uint8_t>>(13);
    gf17 = std::make_unique<GaloisFieldPrime<uint8_t>>(17);
    gf251 = std::make_unique<GaloisFieldPrime<uint8_t>>(251);
    gf65521 = std::make_unique<GaloisFieldPrime<uint16_t>>(65521);

    gf7_hex = std::make_unique<GaloisFieldPrime<uint8_t>>(7, "hex");
    gf7_pow = std::make_unique<GaloisFieldPrime<uint8_t>>(7, "pow");
    gf7_log = std::make_unique<GaloisFieldPrime<uint8_t>>(7, "log");
  }

  std::unique_ptr<GaloisFieldPrime<uint8_t>> gf2, gf3, gf5, gf7, gf13, gf17,
      gf251;
  std::unique_ptr<GaloisFieldPrime<uint16_t>> gf65521;
  std::unique_ptr<GaloisFieldPrime<uint8_t>> gf7_hex, gf7_pow, gf7_log;
};

TEST_F(GaloisFieldPrimeTest, FieldProperties) {
  EXPECT_EQ(gf7->Characteristic(), 7ULL);
  EXPECT_EQ(gf7->Order(), 7ULL);
  EXPECT_EQ(gf7->Degree(), 1U);
  EXPECT_EQ(gf7->Modulus(), 7ULL);
  EXPECT_EQ(gf7->AdditiveIdentity(), 0);
  EXPECT_EQ(gf7->MultiplicativeIdentity(), 1);
  EXPECT_GT(gf7->MemoryBytes(), 0U);

  uint8_t gen = gf7->MultiplicativeGenerator();
  EXPECT_NE(gen, 0);
  EXPECT_NE(gen, 1);
  EXPECT_EQ(gf7->MultiplicativeOrder(gen), 6ULL);
}

TEST_F(GaloisFieldPrimeTest, ExhaustiveDifferentialPrimeVsTable) {
  RunExhaustiveDifferentialTest<uint8_t>(2);
  RunExhaustiveDifferentialTest<uint8_t>(3);
  RunExhaustiveDifferentialTest<uint8_t>(5);
  RunExhaustiveDifferentialTest<uint8_t>(7);
  RunExhaustiveDifferentialTest<uint8_t>(13);
  RunExhaustiveDifferentialTest<uint8_t>(17);
  RunExhaustiveDifferentialTest<uint8_t>(251);
}

TEST_F(GaloisFieldPrimeTest, DifferentialGF65521) {
  GaloisFieldPrime<uint16_t> prime_field(65521);
  GaloisFieldPrimeTable<uint16_t> table_field(65521);

  EXPECT_GT(table_field.MemoryBytes(), prime_field.MemoryBytes());
  EXPECT_EQ(prime_field.MultiplicativeGenerator(),
            table_field.MultiplicativeGenerator());

  std::mt19937_64 rng(20261121ULL);
  std::vector<uint16_t> samples = {0, 1, 2, 3, 17, 255, 32768, 65519, 65520};
  for (int i = 0; i < 64; ++i) {
    samples.push_back(prime_field.Random(rng));
  }

  for (uint16_t a : samples) {
    EXPECT_EQ(prime_field.Neg(a), table_field.Neg(a));
    if (a != 0) {
      EXPECT_EQ(prime_field.Inv(a), table_field.Inv(a));
      EXPECT_EQ(prime_field.Mul(a, prime_field.Inv(a)), 1);
      const uint64_t l_p = prime_field.Log(a);
      const uint64_t l_t = table_field.Log(a);
      EXPECT_EQ(l_p, l_t);
      EXPECT_EQ(prime_field.Pow(prime_field.MultiplicativeGenerator(), l_p), a);
    }
    EXPECT_EQ(prime_field.IsQuadraticResidue(a),
              table_field.IsQuadraticResidue(a));
    if (prime_field.IsQuadraticResidue(a)) {
      const uint16_t r_p = prime_field.Sqrt(a);
      const uint16_t r_t = table_field.Sqrt(a);
      EXPECT_EQ(r_p, r_t);
      EXPECT_EQ(prime_field.Mul(r_p, r_p), a);
    } else {
      EXPECT_THROW(prime_field.Sqrt(a), std::domain_error);
      EXPECT_THROW(table_field.Sqrt(a), std::domain_error);
    }

    for (uint16_t b : samples) {
      EXPECT_EQ(prime_field.Add(a, b), table_field.Add(a, b));
      EXPECT_EQ(prime_field.Sub(a, b), table_field.Sub(a, b));
      EXPECT_EQ(prime_field.Mul(a, b), table_field.Mul(a, b));
      if (b != 0) {
        EXPECT_EQ(prime_field.Div(a, b), table_field.Div(a, b));
      }
    }
  }
}

TEST_F(GaloisFieldPrimeTest, Large32BitPrime4294967291) {
  constexpr uint32_t kP = 4294967291U;  // 2^32 - 5 (prime, 3 mod 4)
  GaloisFieldPrime<uint32_t> gf(kP);

  EXPECT_EQ(gf.Characteristic(), kP);
  EXPECT_EQ(gf.Order(), kP);

  std::mt19937_64 rng(42);
  std::vector<uint32_t> values = {0,      1,          2,      3,
                                  kP / 2, kP / 2 + 1, kP - 2, kP - 1};
  for (int i = 0; i < 64; ++i) {
    values.push_back(gf.Random(rng));
  }

  for (uint32_t a : values) {
    EXPECT_EQ(gf.Add(a, gf.Neg(a)), 0U);
    if (a != 0) {
      const uint32_t inv = gf.Inv(a);
      EXPECT_EQ(gf.Mul(a, inv), 1U);
      EXPECT_EQ(ExactModMul(a, inv, kP), 1ULL);
    }
    for (uint32_t b : values) {
      const uint64_t expected_mul = ExactModMul(a, b, kP);
      EXPECT_EQ(gf.Mul(a, b), expected_mul);
      EXPECT_EQ(gf.Add(a, b), (uint64_t{a} + uint64_t{b}) % kP);
      EXPECT_EQ(gf.Sub(a, b), (uint64_t{a} + kP - uint64_t{b}) % kP);
    }
  }

  // Test Sqrt on squares and non-residues.
  for (uint32_t x : {0U, 1U, 2U, 123456789U, kP - 1}) {
    const uint32_t sq = gf.Mul(x, x);
    EXPECT_TRUE(gf.IsQuadraticResidue(sq));
    const uint32_t r = gf.Sqrt(sq);
    EXPECT_EQ(gf.Mul(r, r), sq);
    EXPECT_LE(r, kP - r);
  }

  // Test BSGS + Pohlig-Hellman Log on 32-bit prime.
  const uint32_t gen = gf.MultiplicativeGenerator();
  EXPECT_EQ(gf.MultiplicativeOrder(gen), uint64_t{kP} - 1);
  for (uint64_t exp : std::initializer_list<uint64_t>{
           0, 1, 2, 1234567, uint64_t{kP} - 2, uint64_t{kP} / 3}) {
    const uint32_t target = gf.Pow(gen, exp);
    const uint64_t found = gf.Log(target);
    EXPECT_EQ(found, exp);
    EXPECT_EQ(gf.Pow(gen, found), target);
  }
}

TEST_F(GaloisFieldPrimeTest, Large64BitMersenneAndNear2Pow63Primes) {
  // Mersenne prime 2^61 - 1.
  constexpr uint64_t kM61 = (uint64_t{1} << 61) - 1;
  GaloisFieldPrime<uint64_t> gf_m61(kM61);
  EXPECT_EQ(gf_m61.Characteristic(), kM61);
  EXPECT_EQ(gf_m61.Order(), kM61);

  std::mt19937_64 rng(61);
  std::vector<uint64_t> vals = {0,        1,       2, 3, kM61 / 2, kM61 / 2 + 1,
                                kM61 - 2, kM61 - 1};
  for (int i = 0; i < 50; ++i) {
    vals.push_back(gf_m61.Random(rng));
  }

  for (uint64_t a : vals) {
    EXPECT_EQ(gf_m61.Add(a, gf_m61.Neg(a)), 0ULL);
    if (a != 0) {
      const uint64_t inv = gf_m61.Inv(a);
      EXPECT_EQ(gf_m61.Mul(a, inv), 1ULL);
      EXPECT_EQ(ExactModMul(a, inv, kM61), 1ULL);
    }
    for (uint64_t b : vals) {
      EXPECT_EQ(gf_m61.Mul(a, b), ExactModMul(a, b, kM61));
      EXPECT_EQ(gf_m61.Add(a, b),
                static_cast<uint64_t>((static_cast<u128>(a) + b) % kM61));
      EXPECT_EQ(
          gf_m61.Sub(a, b),
          static_cast<uint64_t>((static_cast<u128>(a) + kM61 - b) % kM61));
    }
  }

  // Verify Sqrt on 64-bit Mersenne prime.
  for (uint64_t x :
       std::initializer_list<uint64_t>{0, 1, 42, 123456789012345, kM61 - 1}) {
    const uint64_t sq = gf_m61.Mul(x, x);
    EXPECT_TRUE(gf_m61.IsQuadraticResidue(sq));
    const uint64_t r = gf_m61.Sqrt(sq);
    EXPECT_EQ(gf_m61.Mul(r, r), sq);
    EXPECT_LE(r, kM61 - r);
  }

  // Verify Log (Pohlig-Hellman + BSGS) on 2^61 - 1.
  const uint64_t gen = gf_m61.MultiplicativeGenerator();
  for (uint64_t exp : std::initializer_list<uint64_t>{
           0, 1, 17, 987654321012345, kM61 - 2, (kM61 - 1) / 2}) {
    const uint64_t elem = gf_m61.Pow(gen, exp);
    const uint64_t logged = gf_m61.Log(elem);
    EXPECT_EQ(logged, exp);
    EXPECT_EQ(gf_m61.Pow(gen, logged), elem);
  }

  // Largest prime <= 2^63 - 1: 2^63 - 25 = 9223372036854775783 (7 mod 8 -> 3
  // mod 4).
  constexpr uint64_t kMax63Prime = (uint64_t{1} << 63) - 25;
  GaloisFieldPrime<uint64_t> gf_max(kMax63Prime);
  std::vector<uint64_t> max_vals = {
      0, 1, 2, 3, kMax63Prime / 2, kMax63Prime - 2, kMax63Prime - 1};
  for (int i = 0; i < 40; ++i) {
    max_vals.push_back(gf_max.Random(rng));
  }
  for (uint64_t a : max_vals) {
    if (a != 0) {
      const uint64_t inv = gf_max.Inv(a);
      EXPECT_EQ(gf_max.Mul(a, inv), 1ULL);
      EXPECT_EQ(ExactModMul(a, inv, kMax63Prime), 1ULL);
    }
    for (uint64_t b : max_vals) {
      EXPECT_EQ(gf_max.Mul(a, b), ExactModMul(a, b, kMax63Prime));
    }
  }

  // Verify Tonelli-Shanks (p == 1 mod 4) on a 64-bit prime: 2^63 - 165.
  // Check 2^63 - 165 = 9223372036854775643 (3 mod 4); use 2^63 - 555?
  // 2^61 + 15 = 2305843009213693967? Let's use 4294967311?
  // 1000000000000000009 (1 mod 4): check 2^62 + 135 = 4611686018427388039 (3
  // mod 4) 9223372036854775337 = 2^63 - 471 (1 mod 4) — test known prime
  // 2305843009213693969? Use 65537 (Fermat prime, 1 mod 4, s=16) and 998244353
  // (NTT prime, 1 mod 4, s=23).
  constexpr uint32_t kNttPrime = 998244353U;  // 119 * 2^23 + 1
  GaloisFieldPrime<uint32_t> gf_ntt(kNttPrime);
  for (uint32_t x : {0U, 1U, 2U, 3U, 999999U, kNttPrime - 1}) {
    const uint32_t sq = gf_ntt.Mul(x, x);
    EXPECT_TRUE(gf_ntt.IsQuadraticResidue(sq));
    const uint32_t r = gf_ntt.Sqrt(sq);
    EXPECT_EQ(gf_ntt.Mul(r, r), sq);
    EXPECT_LE(r, kNttPrime - r);
  }
  // 3 is a primitive root mod 998244353, hence a quadratic non-residue.
  EXPECT_FALSE(gf_ntt.IsQuadraticResidue(3));
  EXPECT_THROW(gf_ntt.Sqrt(3), std::domain_error);

  // Prime > 2^63 - 1 must be rejected.
  EXPECT_THROW(GaloisFieldPrime<uint64_t>((uint64_t{1} << 63) + 29),
               std::invalid_argument);
}

TEST_F(GaloisFieldPrimeTest, PowAndPowSigned) {
  EXPECT_EQ(gf7->Pow(0, 0), 1);
  EXPECT_EQ(gf7->Pow(0, 5), 0);
  EXPECT_EQ(gf7->Pow(1, 100), 1);
  EXPECT_EQ(gf7->Pow(2, 3), 1);
  EXPECT_EQ(gf7->Pow(3, 6), 1);

  EXPECT_EQ(gf7->PowSigned(3, 0), 1);
  EXPECT_EQ(gf7->PowSigned(3, 2), 2);
  EXPECT_EQ(gf7->PowSigned(3, -1), 5);
  EXPECT_EQ(gf7->PowSigned(3, -2), 4);
  EXPECT_EQ(gf7->Mul(gf7->PowSigned(3, 5), gf7->PowSigned(3, -5)), 1);
  EXPECT_THROW(gf7->PowSigned(0, -1), std::domain_error);

  // Verify INT64_MIN does not trigger signed overflow UB.
  const uint8_t r_min = gf7->PowSigned(3, std::numeric_limits<int64_t>::min());
  const uint64_t mag = static_cast<uint64_t>(1ULL << 63);
  EXPECT_EQ(r_min, gf7->Pow(gf7->Inv(3), mag));

  // Verify Table Pow with huge uint64_t exponent does not overflow.
  GaloisFieldPrimeTable<uint8_t> gf7_tbl(7);
  EXPECT_EQ(gf7_tbl.Pow(3, std::numeric_limits<uint64_t>::max()),
            gf7->Pow(3, std::numeric_limits<uint64_t>::max()));
}

TEST_F(GaloisFieldPrimeTest, LogEdgeCasesAndNonGenerators) {
  uint8_t gen = gf7->MultiplicativeGenerator();
  EXPECT_THROW(gf7->Log(0, gen), std::domain_error);
  EXPECT_THROW(gf7->Log(1, 0), std::invalid_argument);
  EXPECT_THROW(gf7->Log(1, 1), std::invalid_argument);

  // 2 has order 3 in GF(7): <2> = {1, 2, 4}.
  EXPECT_EQ(gf7->Log(1, 2), 0ULL);
  EXPECT_EQ(gf7->Log(2, 2), 1ULL);
  EXPECT_EQ(gf7->Log(4, 2), 2ULL);
  EXPECT_THROW(gf7->Log(3, 2), std::domain_error);
  EXPECT_THROW(gf7->Log(5, 2), std::domain_error);
  EXPECT_THROW(gf7->Log(6, 2), std::domain_error);

  // Same for table-backed field.
  GaloisFieldPrimeTable<uint8_t> gf7_tbl(7);
  EXPECT_EQ(gf7_tbl.Log(1, 2), 0ULL);
  EXPECT_EQ(gf7_tbl.Log(2, 2), 1ULL);
  EXPECT_EQ(gf7_tbl.Log(4, 2), 2ULL);
  EXPECT_THROW(gf7_tbl.Log(3, 2), std::domain_error);
  EXPECT_THROW(gf7_tbl.Log(1, 0), std::invalid_argument);
  EXPECT_THROW(gf7_tbl.Log(1, 1), std::invalid_argument);

  // Alternative primitive root (5 in GF(7)).
  for (uint8_t a = 1; a < 7; ++a) {
    const uint64_t lp = gf7->Log(a, 5);
    const uint64_t lt = gf7_tbl.Log(a, 5);
    EXPECT_EQ(lp, lt);
    EXPECT_EQ(gf7->Pow(5, lp), a);
  }

  // GF(2) log edge cases.
  GaloisFieldPrimeTable<uint8_t> gf2_tbl(2);
  EXPECT_EQ(gf2->Log(1), 0ULL);
  EXPECT_EQ(gf2->Log(1, 1), 0ULL);
  EXPECT_EQ(gf2_tbl.Log(1), 0ULL);
  EXPECT_EQ(gf2_tbl.Log(1, 1), 0ULL);
  EXPECT_THROW(gf2->Log(0), std::domain_error);
  EXPECT_THROW(gf2_tbl.Log(0), std::domain_error);
}

TEST_F(GaloisFieldPrimeTest, MultiplicativeGenerators) {
  auto gens7 = gf7->MultiplicativeGenerators();
  EXPECT_EQ(gens7, (std::vector<uint8_t>{3, 5}));

  auto gens2 = gf2->MultiplicativeGenerators();
  EXPECT_EQ(gens2, (std::vector<uint8_t>{1}));

  GaloisFieldPrimeTable<uint8_t> gf17_tbl(17);
  EXPECT_EQ(gf17->MultiplicativeGenerators(),
            gf17_tbl.MultiplicativeGenerators());
  EXPECT_EQ(gf17->MultiplicativeGenerators().size(), 8U);  // phi(16) = 8
}

TEST_F(GaloisFieldPrimeTest, BatchInv) {
  std::vector<uint8_t> src = {1, 2, 3, 4, 5, 6};
  std::vector<uint8_t> dst(src.size());
  gf7->BatchInv(dst, src);
  for (size_t i = 0; i < src.size(); ++i) {
    EXPECT_EQ(gf7->Mul(src[i], dst[i]), 1);
  }

  GaloisFieldPrimeTable<uint8_t> gf251_tbl(251);
  std::vector<uint8_t> src251(250);
  for (size_t i = 0; i < 250; ++i) {
    src251[i] = static_cast<uint8_t>(i + 1);
  }
  std::vector<uint8_t> dst251(250);
  gf251_tbl.BatchInv(dst251, src251);
  for (size_t i = 0; i < 250; ++i) {
    EXPECT_EQ(gf251_tbl.Mul(src251[i], dst251[i]), 1);
  }

  std::vector<uint8_t> with_zero = {1, 0, 3};
  std::vector<uint8_t> out_zero(3);
  EXPECT_THROW(gf7->BatchInv(out_zero, with_zero), std::domain_error);

  std::vector<uint8_t> wrong_size(2);
  EXPECT_THROW(gf7->BatchInv(wrong_size, with_zero), std::invalid_argument);
}

TEST_F(GaloisFieldPrimeTest, BulkSpanKernelsPrimeAndTable) {
  auto verify_kernels = [](const auto &field) {
    using T = typename std::decay_t<decltype(field)>::element_type;
    std::mt19937_64 rng(777);
    const size_t n = 128;
    std::vector<T> a(n), b(n), dst(n);
    for (size_t i = 0; i < n; ++i) {
      a[i] = field.Random(rng);
      b[i] = field.Random(rng);
    }

    // AddArray
    field.AddArray(dst, a, b);
    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(dst[i], field.Add(a[i], b[i]));
    }

    // In-place aliasing AddArray
    std::vector<T> aliased = a;
    field.AddArray(aliased, aliased, b);
    EXPECT_EQ(aliased, dst);

    // SubArray
    field.SubArray(dst, a, b);
    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(dst[i], field.Sub(a[i], b[i]));
    }

    // MulArray
    field.MulArray(dst, a, b);
    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(dst[i], field.Mul(a[i], b[i]));
    }

    // ScalarMulArray (0, 1, random)
    for (T scalar : {T{0}, T{1}, field.Random(rng)}) {
      field.ScalarMulArray(dst, a, scalar);
      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(dst[i], field.Mul(a[i], scalar));
      }
    }

    // FmaArray (0, 1, random)
    for (T scalar : {T{0}, T{1}, field.Random(rng)}) {
      std::vector<T> acc = b;
      field.FmaArray(acc, a, scalar);
      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(acc[i], field.Add(b[i], field.Mul(a[i], scalar)));
      }
    }

    // DotProduct
    T expected_dot = 0;
    for (size_t i = 0; i < n; ++i) {
      expected_dot = field.Add(expected_dot, field.Mul(a[i], b[i]));
    }
    EXPECT_EQ(field.DotProduct(a, b), expected_dot);

    // Empty spans
    std::vector<T> empty;
    field.AddArray(empty, empty, empty);
    field.MulArray(empty, empty, empty);
    EXPECT_EQ(field.DotProduct(empty, empty), T{0});

    // Size mismatch exceptions
    std::vector<T> short_vec(n - 1);
    EXPECT_THROW(field.AddArray(short_vec, a, b), std::invalid_argument);
    EXPECT_THROW(field.SubArray(short_vec, a, b), std::invalid_argument);
    EXPECT_THROW(field.MulArray(short_vec, a, b), std::invalid_argument);
    EXPECT_THROW(field.ScalarMulArray(short_vec, a, T{1}),
                 std::invalid_argument);
    EXPECT_THROW(field.FmaArray(short_vec, a, T{1}), std::invalid_argument);
    EXPECT_THROW(field.DotProduct(short_vec, b), std::invalid_argument);
  };

  verify_kernels(*gf251);
  verify_kernels(GaloisFieldPrimeTable<uint8_t>(251));
  verify_kernels(GaloisFieldPrime<uint32_t>(4294967291U));
  verify_kernels(GaloisFieldPrime<uint64_t>((uint64_t{1} << 61) - 1));
  verify_kernels(GaloisFieldPrimeTable<uint16_t>(65521));
}

TEST_F(GaloisFieldPrimeTest, Uint8PrintingAndStreamFlagPreservation) {
  // Verify uint8_t prints as digits ("7"), NOT raw ASCII control char "\x07"!
  std::ostringstream oss;
  gf7->Print(oss);
  EXPECT_EQ(oss.str(), "Galois Prime Field GF(7) [Rep: INT]");
  EXPECT_EQ(oss.str().find('\x07'), std::string::npos);

  EXPECT_EQ(gf251->ToString(7), "7");
  EXPECT_NE(gf251->ToString(7), std::string(1, '\x07'));
  EXPECT_EQ(gf251->ToString(250), "250");

  // Verify HEX printing saves & restores caller stream flags.
  std::ostringstream hex_oss;
  gf7_hex->Print(6, hex_oss);
  hex_oss << " " << 16;
  EXPECT_EQ(hex_oss.str(), "0x6 16");

  // Verify POW and LOG representations and custom generator name.
  EXPECT_EQ(gf7_pow->ToString(0), "0");
  EXPECT_EQ(gf7_pow->ToString(1), "g^0");
  EXPECT_EQ(gf7_pow->ToString(3), "g^1");

  EXPECT_EQ(gf7_log->ToString(0), "undefined");
  EXPECT_EQ(gf7_log->ToString(1), "0");
  EXPECT_EQ(gf7_log->ToString(3), "1");

  GaloisFieldPrime<uint8_t> custom_gen(7, "pow", true, "alpha");
  EXPECT_EQ(custom_gen.GeneratorName(), "alpha");
  EXPECT_EQ(custom_gen.ToString(1), "alpha^0");
  EXPECT_EQ(custom_gen.ToString(3), "alpha^1");
  EXPECT_EQ(custom_gen.SetElementValue("alpha^2"), 2);
  EXPECT_THROW(custom_gen.SetElementValue("g^2"), std::invalid_argument);
}

TEST_F(GaloisFieldPrimeTest, RepresentationValidation) {
  EXPECT_EQ(gf7->GetRepresentation(), FieldRepresentation::INT);
  EXPECT_EQ(gf7_hex->GetRepresentation(), FieldRepresentation::HEX);
  EXPECT_EQ(gf7_pow->GetRepresentation(), FieldRepresentation::POW);
  EXPECT_EQ(gf7_log->GetRepresentation(), FieldRepresentation::LOG);

  gf7->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(gf7->GetRepresentation(), FieldRepresentation::HEX);
  EXPECT_THROW(gf7->SetRepresentation(FieldRepresentation::POLY),
               std::invalid_argument);

  EXPECT_THROW(GaloisFieldPrime<uint8_t>(7, "poly"), std::invalid_argument);
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(7, "unknown"), std::invalid_argument);
}

TEST_F(GaloisFieldPrimeTest, SetElementValueStringFormats) {
  const uint8_t gen = gf7->MultiplicativeGenerator();

  // Power notation.
  EXPECT_EQ(gf7->SetElementValue("g^0"), 1);
  EXPECT_EQ(gf7->SetElementValue("g^1"), gen);
  EXPECT_EQ(gf7->SetElementValue("g^2"), gf7->Mul(gen, gen));
  EXPECT_EQ(gf7->SetElementValue("g^-1"), gf7->Inv(gen));
  EXPECT_EQ(gf7->SetElementValue("g^-2"), gf7->Inv(gf7->Mul(gen, gen)));
  EXPECT_EQ(gf7->SetElementValue("g^6"), 1);
  EXPECT_EQ(gf7->SetElementValue("g^7"), gen);

  // Decimal integers (including 0, >p, and negative).
  EXPECT_EQ(gf7->SetElementValue("0"), 0);
  EXPECT_EQ(gf7->SetElementValue("-0"), 0);
  EXPECT_EQ(gf7->SetElementValue("5"), 5);
  EXPECT_EQ(gf7->SetElementValue("14"), 0);
  EXPECT_EQ(gf7->SetElementValue("15"), 1);
  EXPECT_EQ(gf7->SetElementValue("-3"), 4);
  EXPECT_EQ(gf7->SetElementValue("-7"), 0);
  EXPECT_EQ(gf7->SetElementValue("-10"), 4);

  // Hexadecimal integers.
  EXPECT_EQ(gf7->SetElementValue("0x0"), 0);
  EXPECT_EQ(gf7->SetElementValue("0x6"), 6);
  EXPECT_EQ(gf7->SetElementValue("0x10"), 2);
  EXPECT_EQ(gf7->SetElementValue("0X1A"), 5);
  EXPECT_EQ(gf7->SetElementValue("-0x3"), 4);

  // Invalid strings.
  EXPECT_THROW(gf7->SetElementValue(""), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("-"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("invalid"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("g^"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("g^-"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("g^abc"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("h^1"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("0x"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("0xzz"), std::invalid_argument);
  EXPECT_THROW(gf7->SetElementValue("5abc"), std::invalid_argument);
}

TEST_F(GaloisFieldPrimeTest, UnreducedInputsSafe) {
  // Verify unreduced inputs never underflow or read OOB in either Prime or
  // Table.
  GaloisFieldPrimeTable<uint8_t> gf7_tbl(7);

  for (const GaloisFieldPrime<uint8_t> *f :
       std::vector<const GaloisFieldPrime<uint8_t> *>{gf7.get(), &gf7_tbl}) {
    EXPECT_EQ(f->Add(255, 254), (255 + 254) % 7);
    EXPECT_EQ(f->Sub(1, 255), (1 + 7 - (255 % 7)) % 7);
    EXPECT_EQ(f->Sub(255, 254), 1);
    EXPECT_EQ(f->Neg(255), (7 - (255 % 7)) % 7);
    EXPECT_EQ(f->Neg(14), 0);
    EXPECT_EQ(f->Mul(255, 254), (255 * 254) % 7);
    EXPECT_EQ(f->Div(255, 10), f->Div(255 % 7, 10 % 7));
    EXPECT_THROW(f->Div(3, 14), std::domain_error);
    EXPECT_THROW(f->Inv(14), std::domain_error);
    EXPECT_THROW(f->Inv(0), std::domain_error);
    EXPECT_EQ(f->Inv(10), f->Inv(3));
    EXPECT_TRUE(f->IsZero(14));
    EXPECT_TRUE(f->IsOne(15));
    EXPECT_EQ(f->GetElementValue(15), 1);
    EXPECT_EQ(f->SetElementValue(15), 1);
  }
}

TEST_F(GaloisFieldPrimeTest, ConstructorValidationAndCompositeRejection) {
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(0), std::invalid_argument);
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(1), std::invalid_argument);

  // Default prime_testing == true must reject composites!
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(4), std::invalid_argument);
  EXPECT_THROW(GaloisFieldPrime<uint8_t>(9), std::invalid_argument);
  EXPECT_THROW(GaloisFieldPrime<uint16_t>(561), std::invalid_argument);
  EXPECT_THROW(GaloisFieldPrime<uint64_t>(341550071728321ULL),
               std::invalid_argument);

  EXPECT_NO_THROW(GaloisFieldPrime<uint8_t>(7));
  EXPECT_NO_THROW(GaloisFieldPrime<uint8_t>(4, "int", false));

  // Table memory cap (> 2^22 rejected).
  EXPECT_THROW(GaloisFieldPrimeTable<uint32_t>(4194319U),
               std::invalid_argument);
  EXPECT_NO_THROW(GaloisFieldPrimeTable<uint32_t>(65521U));
}

TEST_F(GaloisFieldPrimeTest, RandomAndReproducibility) {
  std::mt19937_64 rng1(12345);
  std::mt19937_64 rng2(12345);

  std::set<uint8_t> seen;
  for (int i = 0; i < 64; ++i) {
    const uint8_t v1 = gf251->Random(rng1);
    const uint8_t v2 = gf251->Random(rng2);
    EXPECT_EQ(v1, v2);
    EXPECT_LT(v1, 251);
    seen.insert(v1);
  }
  EXPECT_GE(seen.size(), 10U);

  for (int i = 0; i < 32; ++i) {
    EXPECT_LT(gf7->Random(), 7);
  }
}

TEST_F(GaloisFieldPrimeTest, CopyAndThreadSafeLazyGenerator) {
  GaloisFieldPrime<uint32_t> field(4294967291U);
  std::vector<uint32_t> results(4, 0);
  std::vector<std::thread> threads;
  for (size_t i = 0; i < 4; ++i) {
    threads.emplace_back([&field, &results, i]() {
      results[i] = field.MultiplicativeGenerator();
    });
  }
  for (auto &t : threads) {
    t.join();
  }
  for (size_t i = 1; i < results.size(); ++i) {
    EXPECT_EQ(results[i], results[0]);
  }

  GaloisFieldPrime<uint32_t> copy(field);
  EXPECT_EQ(copy.MultiplicativeGenerator(), results[0]);
  EXPECT_EQ(copy.Mul(12345, 67890), field.Mul(12345, 67890));
}
