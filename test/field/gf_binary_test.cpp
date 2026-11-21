#include "xgalois/field/gf_binary.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <vector>

#include "xgalois/field/gf_element.hpp"

using namespace xg;

class GaloisFieldBinaryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf2 = std::make_unique<GaloisFieldBinary>();
    gf2_hex = std::make_unique<GaloisFieldBinary>("hex");
  }

  std::unique_ptr<GaloisFieldBinary> gf2;
  std::unique_ptr<GaloisFieldBinary> gf2_hex;
};

TEST_F(GaloisFieldBinaryTest, FieldProperties) {
  EXPECT_EQ(gf2->Characteristic(), 2ULL);
  EXPECT_EQ(gf2->Order(), 2ULL);
  EXPECT_EQ(gf2->Degree(), 1U);
  EXPECT_EQ(gf2->Modulus(), 2ULL);
  EXPECT_EQ(gf2->AdditiveIdentity(), 0);
  EXPECT_EQ(gf2->MultiplicativeIdentity(), 1);
  EXPECT_EQ(gf2->MultiplicativeGenerator(), 1);
  EXPECT_EQ(gf2->MultiplicativeGenerators(), std::vector<uint8_t>({1}));
  EXPECT_GE(gf2->MemoryBytes(), sizeof(GaloisFieldBinary));
}

TEST_F(GaloisFieldBinaryTest, Addition) {
  EXPECT_EQ(gf2->Add(0, 0), 0);
  EXPECT_EQ(gf2->Add(0, 1), 1);
  EXPECT_EQ(gf2->Add(1, 0), 1);
  EXPECT_EQ(gf2->Add(1, 1), 0);
}

TEST_F(GaloisFieldBinaryTest, Subtraction) {
  EXPECT_EQ(gf2->Sub(0, 0), 0);
  EXPECT_EQ(gf2->Sub(0, 1), 1);
  EXPECT_EQ(gf2->Sub(1, 0), 1);
  EXPECT_EQ(gf2->Sub(1, 1), 0);
}

TEST_F(GaloisFieldBinaryTest, Multiplication) {
  EXPECT_EQ(gf2->Mul(0, 0), 0);
  EXPECT_EQ(gf2->Mul(0, 1), 0);
  EXPECT_EQ(gf2->Mul(1, 0), 0);
  EXPECT_EQ(gf2->Mul(1, 1), 1);
}

TEST_F(GaloisFieldBinaryTest, Division) {
  EXPECT_EQ(gf2->Div(0, 1), 0);
  EXPECT_EQ(gf2->Div(1, 1), 1);
  EXPECT_THROW(gf2->Div(0, 0), std::domain_error);
  EXPECT_THROW(gf2->Div(1, 0), std::domain_error);
}

TEST_F(GaloisFieldBinaryTest, Negation) {
  EXPECT_EQ(gf2->Neg(0), 0);
  EXPECT_EQ(gf2->Neg(1), 1);
}

TEST_F(GaloisFieldBinaryTest, Inverse) {
  EXPECT_EQ(gf2->Inv(1), 1);
  EXPECT_THROW(gf2->Inv(0), std::domain_error);
}

TEST_F(GaloisFieldBinaryTest, Power) {
  EXPECT_EQ(gf2->Pow(0, 0), 1);
  EXPECT_EQ(gf2->Pow(0, 1), 0);
  EXPECT_EQ(gf2->Pow(0, 5), 0);
  EXPECT_EQ(gf2->Pow(1, 0), 1);
  EXPECT_EQ(gf2->Pow(1, 1), 1);
  EXPECT_EQ(gf2->Pow(1, 100), 1);
  EXPECT_EQ(gf2->Pow(1, ~uint64_t{0}), 1);
  EXPECT_EQ(gf2->PowSigned(1, -5), 1);
  EXPECT_THROW(gf2->PowSigned(0, -1), std::domain_error);
}

TEST_F(GaloisFieldBinaryTest, SquareRoot) {
  EXPECT_EQ(gf2->Sqrt(0), 0);
  EXPECT_EQ(gf2->Sqrt(1), 1);
}

TEST_F(GaloisFieldBinaryTest, Logarithm) {
  EXPECT_EQ(gf2->Log(1), 0ULL);
  EXPECT_EQ(gf2->Log(1, 1), 0ULL);
  EXPECT_THROW(gf2->Log(0), std::domain_error);
  EXPECT_THROW(gf2->Log(1, 0), std::invalid_argument);
}

TEST_F(GaloisFieldBinaryTest, Random) {
  std::set<uint8_t> values;
  for (int i = 0; i < 100; ++i) {
    values.insert(gf2->Random());
  }
  EXPECT_GE(values.size(), 1U);
  EXPECT_LE(values.size(), 2U);
  for (auto val : values) {
    EXPECT_TRUE(val == 0 || val == 1);
  }

  std::mt19937_64 rng1(42);
  std::mt19937_64 rng2(42);
  for (int i = 0; i < 64; ++i) {
    uint8_t v1 = gf2->Random(rng1);
    uint8_t v2 = gf2->Random(rng2);
    EXPECT_EQ(v1, v2);
    EXPECT_TRUE(v1 == 0 || v1 == 1);
  }
}

TEST_F(GaloisFieldBinaryTest, Representation) {
  EXPECT_EQ(gf2->GetRepresentation(), FieldRepresentation::INT);
  EXPECT_EQ(gf2_hex->GetRepresentation(), FieldRepresentation::HEX);

  gf2->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(gf2->GetRepresentation(), FieldRepresentation::HEX);

  EXPECT_THROW(gf2->SetRepresentation(FieldRepresentation::POLY),
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryTest, ToString) {
  EXPECT_EQ(gf2->ToString(0), "0");
  EXPECT_EQ(gf2->ToString(1), "1");
  EXPECT_EQ(gf2_hex->ToString(0), "0x0");
  EXPECT_EQ(gf2_hex->ToString(1), "0x1");
}

TEST_F(GaloisFieldBinaryTest, PrintAndStreamFlagPreservation) {
  std::ostringstream oss;
  gf2->Print(oss);
  EXPECT_EQ(oss.str(), "GF(2)");

  oss.str("");
  gf2->Print(1, oss);
  EXPECT_EQ(oss.str(), "1");

  oss.str("");
  gf2_hex->Print(1, oss);
  oss << " " << 16;
  EXPECT_EQ(oss.str(), "0x1 16");
}

TEST_F(GaloisFieldBinaryTest, BulkKernels) {
  std::vector<uint8_t> a = {0, 1, 0, 1, 1, 0, 1, 1, 0};
  std::vector<uint8_t> b = {0, 0, 1, 1, 0, 1, 1, 0, 1};
  std::vector<uint8_t> dst(a.size(), 0);

  gf2->AddArray(dst, a, b);
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_EQ(dst[i], gf2->Add(a[i], b[i]));
  }

  gf2->SubArray(dst, a, b);
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_EQ(dst[i], gf2->Sub(a[i], b[i]));
  }

  gf2->MulArray(dst, a, b);
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_EQ(dst[i], gf2->Mul(a[i], b[i]));
  }

  gf2->ScalarMulArray(dst, a, 0);
  for (uint8_t v : dst) {
    EXPECT_EQ(v, 0);
  }
  gf2->ScalarMulArray(dst, a, 1);
  EXPECT_EQ(dst, a);

  std::vector<uint8_t> acc = b;
  gf2->FmaArray(acc, a, 1);
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_EQ(acc[i], gf2->Add(b[i], a[i]));
  }

  uint8_t expected_dot = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    expected_dot = gf2->Add(expected_dot, gf2->Mul(a[i], b[i]));
  }
  EXPECT_EQ(gf2->DotProduct(a, b), expected_dot);
}

class GaloisFieldBinaryExtensionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf8 = std::make_unique<GaloisFieldBinaryExtension<uint8_t>>(3);
    gf16 = std::make_unique<GaloisFieldBinaryExtension<uint16_t>>(4);
  }

  std::unique_ptr<GaloisFieldBinaryExtension<uint8_t>> gf8;
  std::unique_ptr<GaloisFieldBinaryExtension<uint16_t>> gf16;
};

TEST_F(GaloisFieldBinaryExtensionTest, FieldProperties) {
  EXPECT_EQ(gf8->Characteristic(), 2ULL);
  EXPECT_EQ(gf8->Order(), 8ULL);
  EXPECT_EQ(gf8->Degree(), 3U);
  EXPECT_EQ(gf8->Modulus(), 0b1011ULL);  // Conway poly: x^3 + x + 1
  EXPECT_EQ(gf8->AdditiveIdentity(), 0);
  EXPECT_EQ(gf8->MultiplicativeIdentity(), 1);
  EXPECT_EQ(gf8->MultiplicativeGenerator(), 2);

  EXPECT_EQ(gf16->Characteristic(), 2ULL);
  EXPECT_EQ(gf16->Order(), 16ULL);
  EXPECT_EQ(gf16->Degree(), 4U);
  EXPECT_EQ(gf16->Modulus(), 0b10011ULL);  // Conway poly: x^4 + x + 1
}

TEST_F(GaloisFieldBinaryExtensionTest, BasicArithmeticExact) {
  EXPECT_EQ(gf8->Add(0, 0), 0);
  EXPECT_EQ(gf8->Add(1, 1), 0);
  EXPECT_EQ(gf8->Add(3, 5), 6);

  EXPECT_EQ(gf8->Sub(3, 5), gf8->Add(3, 5));

  EXPECT_EQ(gf8->Mul(0, 5), 0);
  EXPECT_EQ(gf8->Mul(1, 5), 5);

  // Over x^3 + x + 1:
  // 2 * 3 = x * (x + 1) = x^2 + x = 6
  EXPECT_EQ(gf8->Mul(2, 3), 6);
  // 4 * 2 = x^2 * x = x^3 = x + 1 = 3
  EXPECT_EQ(gf8->Mul(4, 2), 3);
  // 7 * 7 = (x^2 + x + 1)^2 = x^4 + x^2 + 1 = (x^2 + x) + x^2 + 1 = x + 1 = 3
  EXPECT_EQ(gf8->Mul(7, 7), 3);
}

TEST_F(GaloisFieldBinaryExtensionTest, MultiplicativeStructureAndSqrt) {
  for (uint8_t a = 0; a < 8; ++a) {
    uint8_t s = gf8->Sqrt(a);
    EXPECT_EQ(gf8->Mul(s, s), a)
        << "Sqrt(" << static_cast<int>(a) << ")^2 != " << static_cast<int>(a);
  }

  for (uint8_t a = 1; a < 8; ++a) {
    uint8_t inv_a = gf8->Inv(a);
    EXPECT_EQ(gf8->Mul(a, inv_a), 1)
        << "Element " << static_cast<int>(a) << " * " << static_cast<int>(inv_a)
        << " != 1";

    uint64_t log_a = gf8->Log(a);
    EXPECT_EQ(gf8->Pow(gf8->MultiplicativeGenerator(), log_a), a);
  }

  for (uint8_t a = 0; a < 8; ++a) {
    for (uint8_t b = 1; b < 8; ++b) {
      uint8_t div_result = gf8->Div(a, b);
      EXPECT_EQ(gf8->Mul(div_result, b), a);
    }
  }
}

TEST_F(GaloisFieldBinaryExtensionTest, DegreeOneField) {
  GF2X<uint8_t> gf2_ext(1);
  EXPECT_EQ(gf2_ext.Characteristic(), 2ULL);
  EXPECT_EQ(gf2_ext.Order(), 2ULL);
  EXPECT_EQ(gf2_ext.Degree(), 1U);
  EXPECT_EQ(gf2_ext.GroupOrder(), 1ULL);
  EXPECT_EQ(gf2_ext.MultiplicativeGenerator(), 1);
  EXPECT_EQ(gf2_ext.MultiplicativeGenerators(), std::vector<uint8_t>({1}));
  EXPECT_EQ(gf2_ext.Log(1), 0ULL);
  EXPECT_EQ(gf2_ext.Log(1, 1), 0ULL);
  EXPECT_EQ(gf2_ext.Inv(1), 1);
  EXPECT_EQ(gf2_ext.Sqrt(0), 0);
  EXPECT_EQ(gf2_ext.Sqrt(1), 1);
  EXPECT_EQ(gf2_ext.Add(1, 1), 0);
  EXPECT_EQ(gf2_ext.Mul(1, 1), 1);
  EXPECT_EQ(gf2_ext.Div(1, 1), 1);
  EXPECT_THROW(gf2_ext.Log(0), std::domain_error);
}

TEST_F(GaloisFieldBinaryExtensionTest, FullWidthUint8GF256AndAesPoly) {
  // GF(2^8) in uint8_t with default Conway polynomial x^8 + x^4 + x^3 + x^2 + 1
  GF2X<uint8_t> gf256(8);
  EXPECT_EQ(gf256.Order(), 256ULL);
  EXPECT_EQ(gf256.Degree(), 8U);
  EXPECT_EQ(gf256.Modulus(), 0x11dULL);
  EXPECT_EQ(gf256.MultiplicativeGenerator(), 2);
  EXPECT_EQ(gf256.MultiplicativeOrder(2), 255ULL);

  for (uint16_t a = 0; a < 256; ++a) {
    uint8_t elem = static_cast<uint8_t>(a);
    uint8_t s = gf256.Sqrt(elem);
    EXPECT_EQ(gf256.Mul(s, s), elem);
    if (elem != 0) {
      uint8_t inv = gf256.Inv(elem);
      EXPECT_EQ(gf256.Mul(elem, inv), 1);
      uint64_t l = gf256.Log(elem);
      EXPECT_EQ(gf256.Pow(2, l), elem);
    }
  }

  // AES field: GF(2^8) with x^8 + x^4 + x^3 + x + 1 (0x11b)
  GF2X<uint8_t> aes_gf(8, "hex", "x^8 + x^4 + x^3 + x + 1", "α", true);
  EXPECT_EQ(aes_gf.Modulus(), 0x11bULL);
  // Well-known AES field test vectors:
  EXPECT_EQ(aes_gf.Mul(0x57, 0x83), 0xc1);
  EXPECT_EQ(aes_gf.Inv(0x53), 0xca);
  EXPECT_EQ(aes_gf.Mul(0x53, 0xca), 0x01);
  // In AES field, x (2) is NOT primitive (order 51), x + 1 (3) IS primitive!
  EXPECT_EQ(aes_gf.MultiplicativeOrder(2), 51ULL);
  EXPECT_EQ(aes_gf.MultiplicativeGenerator(), 3);
  EXPECT_EQ(aes_gf.MultiplicativeOrder(3), 255ULL);
}

TEST_F(GaloisFieldBinaryExtensionTest, FullWidthUint16GF65536) {
  GF2X<uint16_t> gf65536(16);
  EXPECT_EQ(gf65536.Order(), 65536ULL);
  EXPECT_EQ(gf65536.Degree(), 16U);
  EXPECT_EQ(gf65536.MultiplicativeGenerator(), 2);
  EXPECT_EQ(gf65536.MultiplicativeOrder(2), 65535ULL);

  std::mt19937_64 rng(12345);
  for (int i = 0; i < 200; ++i) {
    uint16_t a = gf65536.Random(rng);
    uint16_t b = gf65536.Random(rng);
    EXPECT_EQ(gf65536.Mul(a, b), gf65536.Mul(b, a));
    uint16_t s = gf65536.Sqrt(a);
    EXPECT_EQ(gf65536.Mul(s, s), a);
    if (a != 0) {
      uint16_t inv = gf65536.Inv(a);
      EXPECT_EQ(gf65536.Mul(a, inv), 1);
      uint64_t l = gf65536.Log(a);
      EXPECT_EQ(gf65536.Pow(2, l), a);
    }
  }
}

TEST_F(GaloisFieldBinaryExtensionTest, FullWidthUint32GF2Pow32) {
  GF2X<uint32_t> gf32(32);
  EXPECT_EQ(gf32.Order(), uint64_t{1} << 32);
  EXPECT_EQ(gf32.Degree(), 32U);
  EXPECT_EQ(gf32.GroupOrder(), 0xFFFFFFFFULL);
  EXPECT_EQ(gf32.MultiplicativeGenerator(), 2U);
  EXPECT_EQ(gf32.MultiplicativeOrder(2U), 0xFFFFFFFFULL);

  std::mt19937_64 rng(999);
  for (int i = 0; i < 64; ++i) {
    uint32_t a = gf32.Random(rng);
    uint32_t b = gf32.Random(rng);
    uint32_t c = gf32.Random(rng);
    EXPECT_EQ(gf32.Mul(a, gf32.Add(b, c)),
              gf32.Add(gf32.Mul(a, b), gf32.Mul(a, c)));
    uint32_t s = gf32.Sqrt(a);
    EXPECT_EQ(gf32.Mul(s, s), a);
    if (a != 0) {
      EXPECT_EQ(gf32.Mul(a, gf32.Inv(a)), 1U);
      EXPECT_EQ(gf32.PowSigned(a, -3), gf32.Inv(gf32.Pow(a, 3)));
      uint32_t tr = gf32.Trace(a);
      EXPECT_TRUE(tr == 0U || tr == 1U);
      EXPECT_EQ(gf32.Norm(a), 1U);
    }
  }

  // Verify Pohlig-Hellman + BSGS Log in GF(2^32).
  uint32_t elem = gf32.Pow(2U, 123456789ULL);
  EXPECT_EQ(gf32.Log(elem), 123456789ULL);
}

TEST_F(GaloisFieldBinaryExtensionTest, FullWidthUint64GF2Pow64) {
  GF2X<uint64_t> gf64(64);
  EXPECT_EQ(gf64.Order(), 0ULL);  // 2^64 overflows uint64_t; contract returns 0
  EXPECT_EQ(gf64.Degree(), 64U);
  EXPECT_EQ(gf64.GroupOrder(), ~uint64_t{0});
  EXPECT_EQ(gf64.MultiplicativeGenerator(), 2ULL);
  EXPECT_EQ(gf64.MultiplicativeOrder(2ULL), ~uint64_t{0});

  std::mt19937_64 rng(77777);
  for (int i = 0; i < 64; ++i) {
    uint64_t a = gf64.Random(rng);
    uint64_t b = gf64.Random(rng);
    uint64_t c = gf64.Random(rng);
    // Associativity & distributivity
    EXPECT_EQ(gf64.Mul(gf64.Mul(a, b), c), gf64.Mul(a, gf64.Mul(b, c)));
    EXPECT_EQ(gf64.Mul(a, gf64.Add(b, c)),
              gf64.Add(gf64.Mul(a, b), gf64.Mul(a, c)));
    uint64_t s = gf64.Sqrt(a);
    EXPECT_EQ(gf64.Mul(s, s), a);
    if (a != 0) {
      uint64_t inv = gf64.Inv(a);
      EXPECT_EQ(gf64.Mul(a, inv), 1ULL);
      EXPECT_EQ(gf64.Div(b, a), gf64.Mul(b, inv));
      EXPECT_EQ(gf64.Frobenius(a), gf64.Mul(a, a));
      uint64_t tr = gf64.Trace(a);
      EXPECT_TRUE(tr == 0ULL || tr == 1ULL);
      EXPECT_EQ(gf64.Norm(a), 1ULL);
    }
  }

  // Verify Pohlig-Hellman + BSGS Log in GF(2^64)!
  const uint64_t secret_exp = 0x123456789ABCDEF0ULL;
  const uint64_t target = gf64.Pow(2ULL, secret_exp);
  EXPECT_EQ(gf64.Log(target), secret_exp);

  // Custom irreducible degree-64 polynomial: x^64 + x^4 + x^3 + x + 1
  EXPECT_NO_THROW(GF2X<uint64_t>(64, "hex", "x^64 + x^4 + x^3 + x + 1", "α",
                                 /*check_irreducible=*/true));
  // Reducible degree-64 polynomial: x^64 + 1 = (x + 1)^64
  EXPECT_THROW(GF2X<uint64_t>(64, "hex", "x^64 + 1", "α",
                              /*check_irreducible=*/true),
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryExtensionTest, TraceNormFrobeniusAndBatchInv) {
  GF2X<uint8_t> gf256(8);

  // Trace is F_2-linear and maps into {0, 1}.
  for (uint16_t a = 0; a < 256; ++a) {
    uint8_t tr = gf256.Trace(static_cast<uint8_t>(a));
    EXPECT_TRUE(tr == 0 || tr == 1);
    EXPECT_EQ(gf256.Frobenius(static_cast<uint8_t>(a)),
              gf256.Mul(static_cast<uint8_t>(a), static_cast<uint8_t>(a)));
    if (a != 0) {
      EXPECT_EQ(gf256.Norm(static_cast<uint8_t>(a)), 1);
    } else {
      EXPECT_EQ(gf256.Norm(0), 0);
    }
  }
  EXPECT_EQ(gf256.Trace(gf256.Add(0x37, 0x91)),
            gf256.Add(gf256.Trace(0x37), gf256.Trace(0x91)));

  // BatchInv across all 255 nonzero elements.
  std::vector<uint8_t> src(255), dst(255);
  for (size_t i = 0; i < 255; ++i) {
    src[i] = static_cast<uint8_t>(i + 1);
  }
  gf256.BatchInv(dst, src);
  for (size_t i = 0; i < 255; ++i) {
    EXPECT_EQ(gf256.Mul(src[i], dst[i]), 1);
  }
}

TEST_F(GaloisFieldBinaryExtensionTest, BulkSpanKernelsSimdAndPortable) {
  GF2X<uint8_t> gf256(8);
  std::mt19937_64 rng(314159);

  // Test lengths covering 0, sub-SIMD, exact SIMD, SIMD+tail, and large arrays.
  for (size_t n :
       {size_t{0}, size_t{1}, size_t{7}, size_t{15}, size_t{16}, size_t{17},
        size_t{31}, size_t{32}, size_t{33}, size_t{100}}) {
    std::vector<uint8_t> a(n), b(n), dst(n);
    for (size_t i = 0; i < n; ++i) {
      a[i] = gf256.Random(rng);
      b[i] = gf256.Random(rng);
    }

    gf256.AddArray(dst, a, b);
    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(dst[i], gf256.Add(a[i], b[i]));
    }

    gf256.SubArray(dst, a, b);
    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(dst[i], gf256.Sub(a[i], b[i]));
    }

    gf256.MulArray(dst, a, b);
    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(dst[i], gf256.Mul(a[i], b[i]));
    }

    for (uint8_t scalar :
         {uint8_t{0}, uint8_t{1}, uint8_t{2}, uint8_t{0x57}, uint8_t{0xFF}}) {
      gf256.ScalarMulArray(dst, a, scalar);
      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(dst[i], gf256.Mul(a[i], scalar));
      }

      // In-place ScalarMulArray (dst == a)
      std::vector<uint8_t> inplace = a;
      gf256.ScalarMulArray(inplace, inplace, scalar);
      EXPECT_EQ(inplace, dst);

      std::vector<uint8_t> fma_dst = b;
      gf256.FmaArray(fma_dst, a, scalar);
      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(fma_dst[i], gf256.Add(b[i], gf256.Mul(a[i], scalar)));
      }
    }

    uint8_t expected_dot = 0;
    for (size_t i = 0; i < n; ++i) {
      expected_dot = gf256.Add(expected_dot, gf256.Mul(a[i], b[i]));
    }
    EXPECT_EQ(gf256.DotProduct(a, b), expected_dot);
  }

  // Also test uint16_t nibble-split ScalarMulArray / FmaArray (>= 64 elements).
  GF2X<uint16_t> gf16_full(16);
  constexpr size_t kN16 = 85;
  std::vector<uint16_t> a16(kN16), b16(kN16), dst16(kN16);
  for (size_t i = 0; i < kN16; ++i) {
    a16[i] = gf16_full.Random(rng);
    b16[i] = gf16_full.Random(rng);
  }
  const uint16_t s16 = 0xBEEF;
  gf16_full.ScalarMulArray(dst16, a16, s16);
  for (size_t i = 0; i < kN16; ++i) {
    EXPECT_EQ(dst16[i], gf16_full.Mul(a16[i], s16));
  }
  std::vector<uint16_t> fma16 = b16;
  gf16_full.FmaArray(fma16, a16, s16);
  for (size_t i = 0; i < kN16; ++i) {
    EXPECT_EQ(fma16[i], gf16_full.Add(b16[i], gf16_full.Mul(a16[i], s16)));
  }
}

TEST_F(GaloisFieldBinaryExtensionTest, Random) {
  std::set<uint8_t> values;
  for (int i = 0; i < 100; ++i) {
    uint8_t val = gf8->Random();
    EXPECT_LT(val, 8);
    values.insert(val);
  }
  EXPECT_GE(values.size(), 3U);
}

TEST_F(GaloisFieldBinaryExtensionTest, InvalidOperations) {
  EXPECT_THROW(gf8->Inv(0), std::domain_error);
  EXPECT_THROW(gf8->Div(1, 0), std::domain_error);
  EXPECT_THROW(gf8->Log(0), std::domain_error);
  EXPECT_THROW(gf8->SetElementValue(static_cast<uint8_t>(8)),
               std::out_of_range);
}

TEST_F(GaloisFieldBinaryExtensionTest, Constructor) {
  EXPECT_THROW(GaloisFieldBinaryExtension<uint8_t>(0), std::invalid_argument);
  EXPECT_THROW(GaloisFieldBinaryExtension<uint8_t>(9), std::invalid_argument);
  EXPECT_NO_THROW(GaloisFieldBinaryExtension<uint8_t>(8));

  EXPECT_NO_THROW(
      GaloisFieldBinaryExtension<uint8_t>(3, "int", "x^3 + x + 1", "α", true));
  EXPECT_THROW(
      GaloisFieldBinaryExtension<uint8_t>(3, "int", "x^3 + 1", "α", true),
      std::invalid_argument);
  EXPECT_THROW(GaloisFieldBinaryExtension<uint8_t>(3, "int", "x^2 + 1"),
               std::invalid_argument);
}

class GaloisFieldBinaryElementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf2 = std::make_shared<GaloisFieldBinary>();
    gf2_hex = std::make_shared<GaloisFieldBinary>("hex");
  }

  std::shared_ptr<GaloisFieldBinary> gf2;
  std::shared_ptr<GaloisFieldBinary> gf2_hex;
};

TEST_F(GaloisFieldBinaryElementTest, ElementConstruction) {
  GaloisFieldElement<GaloisFieldBinary> elem_0(0, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_1(1, gf2);

  EXPECT_EQ(elem_0.Value(), 0);
  EXPECT_EQ(elem_1.Value(), 1);
  EXPECT_EQ(elem_0.Field(), gf2);
  EXPECT_EQ(elem_1.Field(), gf2);
}

TEST_F(GaloisFieldBinaryElementTest, StringConstructorNotSupported) {
  EXPECT_THROW(GaloisFieldElement<GaloisFieldBinary>("0", gf2),
               std::invalid_argument);
  EXPECT_THROW(GaloisFieldElement<GaloisFieldBinary>("1", gf2),
               std::invalid_argument);
  EXPECT_THROW(GaloisFieldElement<GaloisFieldBinary>("α", gf2),
               std::invalid_argument);
  EXPECT_THROW(GaloisFieldElement<GaloisFieldBinary>("g^1", gf2),
               std::invalid_argument);
}

TEST_F(GaloisFieldBinaryElementTest, NumericAssignmentOperatorReducesMod2) {
  GaloisFieldElement<GaloisFieldBinary> elem(0, gf2);

  elem = static_cast<uint8_t>(1);
  EXPECT_EQ(elem.Value(), 1);

  elem = static_cast<uint8_t>(0);
  EXPECT_EQ(elem.Value(), 0);

  // Bug B-19 fix: SetElementValue(uint8_t) reduces modulo 2.
  elem = static_cast<uint8_t>(2);
  EXPECT_EQ(elem.Value(), 0);

  elem = static_cast<uint8_t>(3);
  EXPECT_EQ(elem.Value(), 1);

  elem = static_cast<uint8_t>(255);
  EXPECT_EQ(elem.Value(), 1);
}

TEST_F(GaloisFieldBinaryElementTest, StringAssignmentNotSupported) {
  GaloisFieldElement<GaloisFieldBinary> elem(0, gf2);

  EXPECT_THROW(elem = "0", std::invalid_argument);
  EXPECT_THROW(elem = "1", std::invalid_argument);
  EXPECT_THROW(elem = "α", std::invalid_argument);
  EXPECT_THROW(elem = "g^0", std::invalid_argument);
}

TEST_F(GaloisFieldBinaryElementTest, ElementCopyAssignment) {
  GaloisFieldElement<GaloisFieldBinary> elem_0(0, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_1(1, gf2);

  elem_0 = elem_1;
  EXPECT_EQ(elem_0.Value(), 1);
  EXPECT_EQ(elem_0.Field(), gf2);

  elem_1 = GaloisFieldElement<GaloisFieldBinary>(0, gf2);
  EXPECT_EQ(elem_1.Value(), 0);
}

TEST_F(GaloisFieldBinaryElementTest, CompoundAssignmentOperators) {
  GaloisFieldElement<GaloisFieldBinary> elem_a(1, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_b(1, gf2);

  elem_a += elem_b;
  EXPECT_EQ(elem_a.Value(), 0);

  elem_a = GaloisFieldElement<GaloisFieldBinary>(0, gf2);
  elem_a += elem_b;
  EXPECT_EQ(elem_a.Value(), 1);

  elem_a -= elem_b;
  EXPECT_EQ(elem_a.Value(), 0);

  elem_a = GaloisFieldElement<GaloisFieldBinary>(1, gf2);
  elem_a *= elem_b;
  EXPECT_EQ(elem_a.Value(), 1);

  elem_a = GaloisFieldElement<GaloisFieldBinary>(0, gf2);
  elem_a *= elem_b;
  EXPECT_EQ(elem_a.Value(), 0);

  elem_a = GaloisFieldElement<GaloisFieldBinary>(1, gf2);
  elem_a /= elem_b;
  EXPECT_EQ(elem_a.Value(), 1);

  elem_a = GaloisFieldElement<GaloisFieldBinary>(0, gf2);
  elem_a /= elem_b;
  EXPECT_EQ(elem_a.Value(), 0);

  GaloisFieldElement<GaloisFieldBinary> elem_zero(0, gf2);
  EXPECT_THROW(elem_a /= elem_zero, std::domain_error);
}

TEST_F(GaloisFieldBinaryElementTest, StrictElementValidation) {
  auto gf2_other = std::make_shared<GaloisFieldBinary>();

  GaloisFieldElement<GaloisFieldBinary> elem_gf2(1, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_other(1, gf2_other);

  EXPECT_THROW(elem_gf2 + elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf2 - elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf2 * elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf2 / elem_other, std::invalid_argument);
}

TEST_F(GaloisFieldBinaryElementTest, ElementPrinting) {
  GaloisFieldElement<GaloisFieldBinary> elem_0(0, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_1(1, gf2);
  GaloisFieldElement<GaloisFieldBinary> elem_hex(1, gf2_hex);

  std::ostringstream oss;

  elem_0.Print(oss);
  EXPECT_EQ(oss.str(), "0");

  oss.str("");
  elem_1.Print(oss);
  EXPECT_EQ(oss.str(), "1");

  oss.str("");
  elem_hex.Print(oss);
  EXPECT_EQ(oss.str(), "0x1");

  oss.str("");
  oss << elem_1;
  EXPECT_EQ(oss.str(), "1");
}

TEST_F(GaloisFieldBinaryElementTest, NullFieldHandling) {
  EXPECT_THROW(GaloisFieldElement<GaloisFieldBinary>(0, nullptr),
               std::invalid_argument);

  GaloisFieldElement<GaloisFieldBinary> elem;
  EXPECT_FALSE(elem.HasField());
}

class GaloisFieldBinaryExtensionElementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf8 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);
    gf16 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4, "poly", "",
                                                                  "β");
  }

  std::shared_ptr<GaloisFieldBinaryExtension<uint8_t>> gf8;
  std::shared_ptr<GaloisFieldBinaryExtension<uint16_t>> gf16;
};

TEST_F(GaloisFieldBinaryExtensionElementTest, ElementConstruction) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_0(0, gf8);
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_5(5, gf8);

  EXPECT_EQ(elem_0.Value(), 0);
  EXPECT_EQ(elem_5.Value(), 5);
  EXPECT_EQ(elem_0.Field(), gf8);
  EXPECT_EQ(elem_5.Field(), gf8);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorAllFormats) {
  GaloisFieldElement<GF2X<uint8_t>> elem_poly("α^2 + α + 1", gf8);
  EXPECT_EQ(elem_poly.Value(), 7);

  GaloisFieldElement<GF2X<uint8_t>> elem_alpha("α", gf8);
  EXPECT_EQ(elem_alpha.Value(), 2);

  GaloisFieldElement<GF2X<uint8_t>> elem_alpha_sq("α^2", gf8);
  EXPECT_EQ(elem_alpha_sq.Value(), 4);

  GaloisFieldElement<GF2X<uint8_t>> elem_zero("0", gf8);
  EXPECT_EQ(elem_zero.Value(), 0);

  GaloisFieldElement<GF2X<uint8_t>> elem_one_dec("1", gf8);
  EXPECT_EQ(elem_one_dec.Value(), 1);

  GaloisFieldElement<GF2X<uint8_t>> elem_five("5", gf8);
  EXPECT_EQ(elem_five.Value(), 5);

  GaloisFieldElement<GF2X<uint8_t>> elem_hex("0x6", gf8);
  EXPECT_EQ(elem_hex.Value(), 6);

  GaloisFieldElement<GF2X<uint8_t>> elem_g0("g^0", gf8);
  EXPECT_EQ(elem_g0.Value(), 1);

  GaloisFieldElement<GF2X<uint8_t>> elem_neg("g^-1", gf8);
  EXPECT_EQ(elem_neg.Value(), gf8->Inv(gf8->MultiplicativeGenerator()));
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorCustomVariable) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint16_t>> elem_beta("β^2 + β",
                                                                     gf16);
  EXPECT_EQ(elem_beta.Value(), 6);

  EXPECT_THROW(
      GaloisFieldElement<GaloisFieldBinaryExtension<uint16_t>>("α^2 + α", gf16),
      std::invalid_argument);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorGenerator) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_g0("g^0", gf8);
  EXPECT_EQ(elem_g0.Value(), 1);

  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_g1("g^1", gf8);
  uint8_t generator = gf8->MultiplicativeGenerator();
  EXPECT_EQ(elem_g1.Value(), generator);

  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_g2("g^2", gf8);
  EXPECT_EQ(elem_g2.Value(), gf8->Mul(generator, generator));
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringConstructorInvalid) {
  EXPECT_THROW(
      GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>("invalid", gf8),
      std::invalid_argument);

  EXPECT_THROW(
      GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>("h^1", gf8),
      std::invalid_argument);

  EXPECT_THROW(
      GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>("α^8", gf8),
      std::invalid_argument);

  EXPECT_THROW(
      GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>("0x8", gf8),
      std::out_of_range);

  EXPECT_THROW(
      GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>("8", gf8),
      std::out_of_range);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, NumericAssignmentOperator) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem(0, gf8);

  elem = static_cast<uint8_t>(5);
  EXPECT_EQ(elem.Value(), 5);

  elem = static_cast<uint8_t>(7);
  EXPECT_EQ(elem.Value(), 7);

  elem = static_cast<uint8_t>(0);
  EXPECT_EQ(elem.Value(), 0);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StringAssignmentOperator) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem(0, gf8);

  elem = "α^2 + α";
  EXPECT_EQ(elem.Value(), 6);

  elem = "α^2 + 1";
  EXPECT_EQ(elem.Value(), 5);

  elem = "α";
  EXPECT_EQ(elem.Value(), 2);

  elem = "g^0";
  EXPECT_EQ(elem.Value(), 1);

  elem = static_cast<uint8_t>(0);
  EXPECT_EQ(elem.Value(), 0);

  elem = "g^1";
  EXPECT_EQ(elem.Value(), gf8->MultiplicativeGenerator());
}

TEST_F(GaloisFieldBinaryExtensionElementTest, CompoundAssignmentOperators) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_a(3, gf8);
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_b(5, gf8);

  elem_a += elem_b;
  EXPECT_EQ(elem_a.Value(), 6);

  elem_a -= elem_b;
  EXPECT_EQ(elem_a.Value(), 3);

  elem_a = GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>(2, gf8);
  elem_b = GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>(3, gf8);
  elem_a *= elem_b;
  EXPECT_EQ(elem_a.Value(), 6);

  elem_a = GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>(6, gf8);
  elem_b = GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>>(2, gf8);
  elem_a /= elem_b;
  EXPECT_EQ(elem_a.Value(), 3);

  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_zero(0, gf8);
  EXPECT_THROW(elem_a /= elem_zero, std::domain_error);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, ElementArithmetic) {
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_2(2, gf8);
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_3(3, gf8);

  auto sum = elem_2 + elem_3;
  EXPECT_EQ(sum.Value(), 1);

  auto diff = elem_3 - elem_2;
  EXPECT_EQ(diff.Value(), 1);

  auto prod = elem_2 * elem_3;
  EXPECT_EQ(prod.Value(), 6);

  auto quot = prod / elem_2;
  EXPECT_EQ(quot.Value(), 3);

  auto neg = -elem_2;
  EXPECT_EQ(neg.Value(), 2);

  auto power = elem_2 ^ 3;
  EXPECT_EQ(power.Value(), 3);  // x^3 = x + 1 = 3

  auto pow_method = elem_2.Pow(3);
  EXPECT_EQ(pow_method.Value(), 3);

  auto inv = elem_3.Inv();
  EXPECT_EQ((elem_3 * inv).Value(), 1);

  // Bug B-03 fix: Sqrt(a)^2 MUST equal a!
  auto sqrt_elem = elem_2.Sqrt();
  auto sqrt_squared = sqrt_elem * sqrt_elem;
  EXPECT_EQ(sqrt_squared.Value(), elem_2.Value());
}

TEST_F(GaloisFieldBinaryExtensionElementTest, StrictElementValidation) {
  auto gf8_other = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(3);

  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_gf8(2, gf8);
  GaloisFieldElement<GaloisFieldBinaryExtension<uint8_t>> elem_other(3,
                                                                     gf8_other);

  EXPECT_THROW(elem_gf8 + elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf8 - elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf8 * elem_other, std::invalid_argument);
  EXPECT_THROW(elem_gf8 / elem_other, std::invalid_argument);
}

TEST_F(GaloisFieldBinaryExtensionElementTest, ElementPrintingExact) {
  auto gf8_custom =
      std::make_shared<GF2X<uint8_t>>(3, "int", "", "α", false, "gen");
  GaloisFieldElement<GF2X<uint8_t>> elem_5(5, gf8_custom);
  GaloisFieldElement<GF2X<uint8_t>> elem_0(0, gf8_custom);

  // INT: must print "5", NOT raw ASCII 0x05!
  EXPECT_EQ(elem_5.ToString(), "5");
  EXPECT_EQ(elem_0.ToString(), "0");

  std::ostringstream oss;
  oss << elem_5;
  EXPECT_EQ(oss.str(), "5");

  // HEX: must print "0x5" and restore stream flags!
  gf8_custom->SetRepresentation(FieldRepresentation::HEX);
  EXPECT_EQ(elem_5.ToString(), "0x5");
  oss.str("");
  oss << elem_5 << " " << 16;
  EXPECT_EQ(oss.str(), "0x5 16");

  // POLY: 5 = 0b101 = α^2 + 1
  gf8_custom->SetRepresentation(FieldRepresentation::POLY);
  EXPECT_EQ(elem_5.ToString(), "α^2 + 1");
  EXPECT_EQ(elem_0.ToString(), "0");

  // POW: uses generator_name_ ("gen^...") and handles 0 gracefully!
  gf8_custom->SetRepresentation(FieldRepresentation::POW);
  EXPECT_EQ(elem_5.ToString(), "gen^" + std::to_string(gf8_custom->Log(5)));
  EXPECT_EQ(elem_0.ToString(), "0");

  // LOG: handles 0 gracefully without throwing!
  gf8_custom->SetRepresentation(FieldRepresentation::LOG);
  EXPECT_EQ(elem_5.ToString(), std::to_string(gf8_custom->Log(5)));
  EXPECT_EQ(elem_0.ToString(), "-inf");
}
