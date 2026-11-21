#include "xgalois/field/gf_element.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

static_assert(!noexcept(
    std::declval<std::ostream &>()
    << std::declval<const GaloisFieldElement<GaloisFieldPrime<uint8_t>> &>()));

class GaloisFieldElementTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    gf2 = std::make_shared<GaloisFieldBinary>();

    elem_0 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(0, gf7);
    elem_1 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(1, gf7);
    elem_2 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(2, gf7);
    elem_3 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(3, gf7);
    elem_4 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(4, gf7);
    elem_5 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(5, gf7);
    elem_6 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(6, gf7);

    bin_0 = GaloisFieldElement<GaloisFieldBinary>(0, gf2);
    bin_1 = GaloisFieldElement<GaloisFieldBinary>(1, gf2);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7;
  std::shared_ptr<GaloisFieldBinary> gf2;

  GaloisFieldElement<GaloisFieldPrime<uint8_t>> elem_0, elem_1, elem_2, elem_3,
      elem_4, elem_5, elem_6;
  GaloisFieldElement<GaloisFieldBinary> bin_0, bin_1;
};

TEST_F(GaloisFieldElementTest, Construction) {
  EXPECT_EQ(elem_3.Value(), 3);
  EXPECT_EQ(elem_3.Field(), gf7);
  EXPECT_TRUE(elem_3.HasField());

  EXPECT_THROW(GaloisFieldElement<GaloisFieldPrime<uint8_t>>(3, nullptr),
               std::invalid_argument);
}

TEST_F(GaloisFieldElementTest, Addition) {
  auto result = elem_3 + elem_5;
  EXPECT_EQ(result.Value(), 1);
  EXPECT_EQ(result.Field(), gf7);

  auto zero_sum = elem_3 + elem_4;
  EXPECT_EQ(zero_sum.Value(), 0);
}

TEST_F(GaloisFieldElementTest, Subtraction) {
  auto result = elem_5 - elem_3;
  EXPECT_EQ(result.Value(), 2);

  auto neg_result = elem_3 - elem_5;
  EXPECT_EQ(neg_result.Value(), 5);
}

TEST_F(GaloisFieldElementTest, Multiplication) {
  auto result = elem_3 * elem_5;
  EXPECT_EQ(result.Value(), 1);

  auto zero_result = elem_0 * elem_5;
  EXPECT_EQ(zero_result.Value(), 0);

  auto identity_result = elem_1 * elem_5;
  EXPECT_EQ(identity_result.Value(), 5);
}

TEST_F(GaloisFieldElementTest, Division) {
  auto result = elem_1 / elem_3;
  EXPECT_EQ(result.Value(), 5);

  auto self_div = elem_5 / elem_5;
  EXPECT_EQ(self_div.Value(), 1);

  EXPECT_THROW(elem_3 / elem_0, std::domain_error);
}

TEST_F(GaloisFieldElementTest, Negation) {
  auto neg_3 = -elem_3;
  EXPECT_EQ(neg_3.Value(), 4);

  auto neg_0 = -elem_0;
  EXPECT_EQ(neg_0.Value(), 0);
}

TEST_F(GaloisFieldElementTest, PowerAndPowSigned) {
  auto power_result = elem_3 ^ 2;
  EXPECT_EQ(power_result.Value(), 2);

  auto pow_method = elem_3.Pow(3);
  EXPECT_EQ(pow_method.Value(), 6);

  auto pow_zero = elem_5 ^ 0;
  EXPECT_EQ(pow_zero.Value(), 1);

  EXPECT_EQ(elem_3.PowSigned(0).Value(), 1);
  EXPECT_EQ(elem_3.PowSigned(2).Value(), 2);
  EXPECT_EQ(elem_3.PowSigned(-1), elem_3.Inv());
  EXPECT_EQ(elem_3.PowSigned(-2), elem_3.Inv().Pow(2));
  EXPECT_EQ((elem_3.PowSigned(4) * elem_3.PowSigned(-4)).Value(), 1);
  EXPECT_THROW(elem_0.PowSigned(-1), std::domain_error);
}

TEST_F(GaloisFieldElementTest, Inverse) {
  auto inv_3 = elem_3.Inv();
  EXPECT_EQ(inv_3.Value(), 5);

  auto check = elem_3 * inv_3;
  EXPECT_EQ(check.Value(), 1);

  EXPECT_THROW(elem_0.Inv(), std::domain_error);
}

TEST_F(GaloisFieldElementTest, SquareRootTonelliShanks) {
  EXPECT_EQ(bin_1.Sqrt().Value(), 1);
  EXPECT_EQ(bin_0.Sqrt().Value(), 0);

  EXPECT_EQ(elem_0.Sqrt().Value(), 0);
  EXPECT_EQ(elem_1.Sqrt().Value(), 1);

  auto sqrt_4 = elem_4.Sqrt();
  EXPECT_EQ((sqrt_4 * sqrt_4).Value(), 4);
  EXPECT_EQ(sqrt_4.Value(), 2);

  auto sqrt_2 = elem_2.Sqrt();
  EXPECT_EQ((sqrt_2 * sqrt_2).Value(), 2);
  EXPECT_TRUE(sqrt_2.Value() == 3 || sqrt_2.Value() == 4);

  // 3, 5, 6 are quadratic non-residues mod 7.
  EXPECT_THROW(elem_3.Sqrt(), std::domain_error);
  EXPECT_THROW(elem_5.Sqrt(), std::domain_error);
  EXPECT_THROW(elem_6.Sqrt(), std::domain_error);
}

TEST_F(GaloisFieldElementTest, FrobeniusTraceNormAndMultiplicativeOrder) {
  // Over prime field GF(7), Frobenius, Trace, and Norm are the identity map.
  for (uint8_t v = 0; v < 7; ++v) {
    GaloisFieldElement<GaloisFieldPrime<uint8_t>> e(v, gf7);
    EXPECT_EQ(e.Frobenius(), e);
    EXPECT_EQ(e.Trace(), e);
    EXPECT_EQ(e.Norm(), e);
  }

  EXPECT_EQ(elem_1.MultiplicativeOrder(), 1U);
  EXPECT_EQ(elem_6.MultiplicativeOrder(), 2U);
  EXPECT_EQ(elem_2.MultiplicativeOrder(), 3U);
  EXPECT_EQ(elem_4.MultiplicativeOrder(), 3U);
  EXPECT_EQ(elem_3.MultiplicativeOrder(), 6U);
  EXPECT_EQ(elem_5.MultiplicativeOrder(), 6U);
  EXPECT_THROW(elem_0.MultiplicativeOrder(), std::domain_error);

  // Non-trivial extension field GF(3^2): Frobenius(a) == a^3, Trace/Norm in F3.
  auto gf9 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3, 2));
  GaloisFieldElement<GFPX<uint8_t>> a("α + 1", gf9);
  GaloisFieldElement<GFPX<uint8_t>> b("2*α + 1", gf9);
  EXPECT_EQ(a.Frobenius(), a.Pow(3));
  EXPECT_EQ((a + b).Trace(), a.Trace() + b.Trace());
  EXPECT_EQ((a * b).Norm(), a.Norm() * b.Norm());
  EXPECT_LE(a.Trace().Value().Degree(), 0);
  EXPECT_LE(a.Norm().Value().Degree(), 0);
}

TEST_F(GaloisFieldElementTest, CompoundAssignment) {
  auto elem_copy = elem_3;
  elem_copy += elem_5;
  EXPECT_EQ(elem_copy.Value(), 1);

  elem_copy = elem_3;
  elem_copy -= elem_5;
  EXPECT_EQ(elem_copy.Value(), 5);

  elem_copy = elem_3;
  elem_copy *= elem_5;
  EXPECT_EQ(elem_copy.Value(), 1);

  elem_copy = elem_1;
  elem_copy /= elem_3;
  EXPECT_EQ(elem_copy.Value(), 5);
}

TEST_F(GaloisFieldElementTest, Comparison) {
  auto elem_3_copy = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(3, gf7);

  EXPECT_TRUE(elem_3 == elem_3_copy);
  EXPECT_FALSE(elem_3 == elem_4);
  EXPECT_TRUE(elem_3 != elem_4);
  EXPECT_FALSE(elem_3 != elem_3_copy);
  EXPECT_TRUE(elem_3 < elem_4);
  EXPECT_TRUE(elem_3 <= elem_3_copy);
  EXPECT_TRUE(elem_4 > elem_3);
  EXPECT_TRUE(elem_4 >= elem_3);
}

TEST_F(GaloisFieldElementTest, HashSupport) {
  using Elem = GaloisFieldElement<GaloisFieldPrime<uint8_t>>;
  std::hash<Elem> hasher;

  Elem elem_3_copy(3, gf7);
  EXPECT_EQ(hasher(elem_3), hasher(elem_3_copy));

  Elem null_elem;
  EXPECT_EQ(hasher(null_elem), 0U);

  std::unordered_set<Elem> set{elem_0, elem_1, elem_2, elem_3,
                               elem_4, elem_5, elem_6, elem_3_copy};
  EXPECT_EQ(set.size(), 7U);
}

TEST_F(GaloisFieldElementTest, NullFieldExceptionsAndNonNoexceptStream) {
  using Elem = GaloisFieldElement<GaloisFieldPrime<uint8_t>>;
  Elem null_elem;

  EXPECT_FALSE(null_elem.HasField());
  EXPECT_EQ(null_elem.Field(), nullptr);
  EXPECT_THROW(null_elem.Value(), std::runtime_error);
  EXPECT_THROW(null_elem.RawValue(), std::runtime_error);
  EXPECT_THROW(null_elem.IsZero(), std::runtime_error);
  EXPECT_THROW(null_elem.IsOne(), std::runtime_error);
  EXPECT_THROW(null_elem.Inv(), std::runtime_error);
  EXPECT_THROW(-null_elem, std::runtime_error);
  EXPECT_THROW(null_elem.ToString(), std::runtime_error);
  EXPECT_THROW(null_elem + elem_3, std::invalid_argument);
  EXPECT_THROW(elem_3 + null_elem, std::invalid_argument);

  std::ostringstream oss;
  EXPECT_THROW(oss << null_elem, std::runtime_error);
}

TEST_F(GaloisFieldElementTest, BinaryFieldSpecific) {
  auto sum = bin_0 + bin_1;
  EXPECT_EQ(sum.Value(), 1);

  auto sum2 = bin_1 + bin_1;
  EXPECT_EQ(sum2.Value(), 0);

  auto prod = bin_1 * bin_1;
  EXPECT_EQ(prod.Value(), 1);

  auto sqrt_bin = bin_1.Sqrt();
  EXPECT_EQ(sqrt_bin.Value(), 1);
}

TEST_F(GaloisFieldElementTest, StringAndStreamOperations) {
  EXPECT_EQ(elem_3.ToString(), "3");
  EXPECT_EQ(bin_1.ToString(), "1");

  std::ostringstream oss;
  oss << elem_3;
  EXPECT_EQ(oss.str(), "3");

  oss.str("");
  oss << bin_1;
  EXPECT_EQ(oss.str(), "1");
}

TEST_F(GaloisFieldElementTest, LogarithmOperations) {
  uint8_t gen = gf7->MultiplicativeGenerator();
  auto gen_elem = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(gen, gf7);

  uint64_t log_val = elem_3.Log();
  EXPECT_EQ(gen_elem.Pow(log_val).Value(), 3);

  uint64_t log_val_base = gf7->Log(elem_3.Value(), gen_elem.Value());
  EXPECT_EQ(gen_elem.Pow(log_val_base).Value(), 3);

  EXPECT_THROW(elem_0.Log(), std::domain_error);
}

TEST_F(GaloisFieldElementTest, BatchInv) {
  std::vector<uint8_t> src = {1, 2, 3, 4, 5, 6};
  std::vector<uint8_t> dst(src.size(), 0);

  gf7->BatchInv(std::span<uint8_t>(dst), std::span<const uint8_t>(src));
  for (size_t i = 0; i < src.size(); ++i) {
    EXPECT_EQ(gf7->Mul(src[i], dst[i]), 1U);
    EXPECT_EQ(dst[i], gf7->Inv(src[i]));
  }

  // Empty span is a valid no-op.
  std::vector<uint8_t> empty_src;
  std::vector<uint8_t> empty_dst;
  EXPECT_NO_THROW(gf7->BatchInv(std::span<uint8_t>(empty_dst),
                                std::span<const uint8_t>(empty_src)));

  // Single element.
  std::vector<uint8_t> one_src = {3};
  std::vector<uint8_t> one_dst = {0};
  gf7->BatchInv(std::span<uint8_t>(one_dst), std::span<const uint8_t>(one_src));
  EXPECT_EQ(one_dst[0], 5U);

  // Size mismatch throws invalid_argument.
  std::vector<uint8_t> wrong_dst(2, 0);
  EXPECT_THROW(gf7->BatchInv(std::span<uint8_t>(wrong_dst),
                             std::span<const uint8_t>(src)),
               std::invalid_argument);

  // Zero entry throws domain_error.
  std::vector<uint8_t> with_zero = {1, 2, 0, 4};
  std::vector<uint8_t> out_zero(with_zero.size(), 0);
  EXPECT_THROW(gf7->BatchInv(std::span<uint8_t>(out_zero),
                             std::span<const uint8_t>(with_zero)),
               std::domain_error);
}

TEST_F(GaloisFieldElementTest, GFArrayOperations) {
  using Arr = GFArray<GaloisFieldPrime<uint8_t>>;

  auto zeros = Arr::Zeros(gf7, {2, 3});
  auto ones = Arr::Ones(gf7, {2, 3});
  EXPECT_EQ(zeros.Size(), 6U);
  EXPECT_EQ(ones.Size(), 6U);
  EXPECT_EQ(zeros.Shape(0), 2U);
  EXPECT_EQ(zeros.Shape(1), 3U);

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(zeros.At(i, j).Value(), 0U);
      EXPECT_EQ(ones.At(i, j).Value(), 1U);
    }
  }

  std::mt19937_64 rng1(777);
  std::mt19937_64 rng2(777);
  auto a = Arr::Random(gf7, {2, 3}, rng1);
  auto a_same = Arr::Random(gf7, {2, 3}, rng2);
  EXPECT_EQ(a, a_same);

  auto b = Arr::Random(gf7, {2, 3}, rng1);

  auto sum = a + b;
  auto diff = sum - b;
  EXPECT_EQ(diff, a);

  auto neg_a = -a;
  EXPECT_EQ(a + neg_a, zeros);

  auto prod = a * b;
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(prod.At(i, j), a.At(i, j) * b.At(i, j));
    }
  }

  auto scaled_r = a * elem_3;
  auto scaled_l = elem_3 * a;
  EXPECT_EQ(scaled_r, scaled_l);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      EXPECT_EQ(scaled_r.At(i, j), a.At(i, j) * elem_3);
    }
  }

  auto accum = a;
  accum += b;
  EXPECT_EQ(accum, sum);
  accum -= b;
  EXPECT_EQ(accum, a);

  // Round-trip ToElementArray / FromElementArray.
  auto elem_arr = a.ToElementArray();
  auto from_elem = Arr::FromElementArray(elem_arr);
  EXPECT_EQ(from_elem, a);

  // Error handling.
  auto wrong_shape = Arr::Zeros(gf7, {3, 2});
  EXPECT_THROW(a + wrong_shape, std::invalid_argument);

  auto gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);
  auto other_field_arr = Arr::Zeros(gf5, {2, 3});
  EXPECT_THROW(a + other_field_arr, std::invalid_argument);

  EXPECT_THROW(Arr::Zeros(nullptr, {2, 2}), std::invalid_argument);
}

class GaloisFieldElementStrictTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);

    elem_3_gf7 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(3, gf7);
    elem_2_gf5 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(2, gf5);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7, gf5;
  GaloisFieldElement<GaloisFieldPrime<uint8_t>> elem_3_gf7, elem_2_gf5;
};

TEST_F(GaloisFieldElementStrictTest, SameFieldOperations) {
  auto elem_4_gf7 = GaloisFieldElement<GaloisFieldPrime<uint8_t>>(4, gf7);

  auto sum = elem_3_gf7 + elem_4_gf7;
  EXPECT_EQ(sum.Value(), 0);

  auto prod = elem_3_gf7 * elem_4_gf7;
  EXPECT_EQ(prod.Value(), 5);
}

TEST_F(GaloisFieldElementStrictTest, DifferentFieldOperations) {
  EXPECT_THROW(elem_3_gf7 + elem_2_gf5, std::invalid_argument);
  EXPECT_THROW(elem_3_gf7 * elem_2_gf5, std::invalid_argument);
  EXPECT_THROW(elem_3_gf7 - elem_2_gf5, std::invalid_argument);
  EXPECT_THROW(elem_3_gf7 / elem_2_gf5, std::invalid_argument);
}
