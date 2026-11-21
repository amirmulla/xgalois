#include "xgalois/poly/poly_dense.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/utils/poly.hpp"

using namespace xg;

namespace {

template <typename Field>
PolynomialDense<Field> NaiveSchoolbookMul(const PolynomialDense<Field> &a,
                                          const PolynomialDense<Field> &b) {
  auto field = a.Field();
  if (a.IsZero() || b.IsZero()) {
    return PolynomialDense<Field>::Zero(field, a.GetVariable());
  }
  using Raw = typename Field::element_type;
  const Raw zero = field->AdditiveIdentity();
  std::vector<Raw> out(a.Size() + b.Size() - 1, zero);
  for (size_t i = 0; i < a.Size(); ++i) {
    for (size_t j = 0; j < b.Size(); ++j) {
      out[i + j] =
          field->Add(out[i + j], field->Mul(a.RawCoeff(i), b.RawCoeff(j)));
    }
  }
  return PolynomialDense<Field>::FromRaw(std::move(out), field,
                                         a.GetVariable());
}

}  // namespace

class PolynomialDenseTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);

    zero = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(0, gf7);
    one = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(1, gf7);
    two = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(2, gf7);
    three = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(3, gf7);
    four = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(4, gf7);
    five = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(5, gf7);
    six = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(6, gf7);

    gf2 = std::make_shared<GaloisFieldBinary>();
    gf2_zero = GaloisFieldElementBase<GaloisFieldBinary>(0, gf2);
    gf2_one = GaloisFieldElementBase<GaloisFieldBinary>(1, gf2);
  }

  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf7;
  std::shared_ptr<GaloisFieldBinary> gf2;

  GaloisFieldElementBase<GaloisFieldPrime<uint8_t>> zero, one, two, three, four,
      five, six;
  GaloisFieldElementBase<GaloisFieldBinary> gf2_zero, gf2_one;
};

TEST_F(PolynomialDenseTest, DefaultConstructorAndStaticFactories) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> unattached;
  EXPECT_FALSE(unattached.HasField());
  EXPECT_TRUE(unattached.IsZero());
  EXPECT_EQ(unattached.Degree(), -1);
  EXPECT_EQ(unattached.Size(), 0u);
  EXPECT_EQ(unattached.ToString(), "0");

  auto z = PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(gf7);
  EXPECT_TRUE(z.HasField());
  EXPECT_TRUE(z.IsZero());
  EXPECT_FALSE(z.IsOne());
  EXPECT_FALSE(z.IsMonic());
  EXPECT_EQ(z.Degree(), -1);
  EXPECT_EQ(z.Size(), 0u);
  EXPECT_EQ(z[0], zero);
  EXPECT_EQ(z.LeadingCoeff(), zero);
  EXPECT_THROW(z.Monic(), std::domain_error);

  auto o = PolynomialDense<GaloisFieldPrime<uint8_t>>::One(gf7);
  EXPECT_FALSE(o.IsZero());
  EXPECT_TRUE(o.IsOne());
  EXPECT_TRUE(o.IsMonic());
  EXPECT_EQ(o.Degree(), 0);
  EXPECT_EQ(o.Size(), 1u);
  EXPECT_EQ(o[0], one);

  auto mono =
      PolynomialDense<GaloisFieldPrime<uint8_t>>::Monomial(gf7, 4, 5, "t");
  EXPECT_EQ(mono.Degree(), 4);
  EXPECT_EQ(mono.GetVariable(), "t");
  EXPECT_EQ(mono.RawCoeff(4), 5);
  EXPECT_EQ(mono.RawCoeff(3), 0);
  EXPECT_EQ(mono.RawCoeff(10), 0);

  auto zero_mono =
      PolynomialDense<GaloisFieldPrime<uint8_t>>::Monomial(gf7, 4, 0);
  EXPECT_TRUE(zero_mono.IsZero());
}

TEST_F(PolynomialDenseTest, ConstructorAndBasicProperties) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), 2);
  EXPECT_EQ(poly.Size(), 3u);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], two);
  EXPECT_EQ(poly[2], three);
  EXPECT_EQ(poly.LeadingCoeff(), three);
  EXPECT_FALSE(poly.IsMonic());
  EXPECT_EQ(poly.GetVariable(), "x");

  auto monic = poly.Monic();
  EXPECT_TRUE(monic.IsMonic());
  EXPECT_EQ(monic * three, poly);
}

TEST_F(PolynomialDenseTest, CustomVariable) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {one,
                                                                           two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs, "t");

  EXPECT_EQ(poly.GetVariable(), "t");

  poly.SetVariable("y");
  EXPECT_EQ(poly.GetVariable(), "y");

  EXPECT_THROW(poly.SetVariable(""), std::invalid_argument);
  EXPECT_THROW(PolynomialDense<GaloisFieldPrime<uint8_t>>(coeffs, ""),
               std::invalid_argument);
}

TEST_F(PolynomialDenseTest, ZeroPolynomial) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_TRUE(poly.IsZero());
  EXPECT_EQ(poly.Degree(), -1);
  EXPECT_EQ(poly.Size(), 0u);
  EXPECT_EQ(poly[0], zero);
}

TEST_F(PolynomialDenseTest, ConstantPolynomial) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), 0);
  EXPECT_EQ(poly.Size(), 1u);
  EXPECT_EQ(poly[0], five);
}

TEST_F(PolynomialDenseTest, LeadingZeroTrimming) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, zero, zero};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly.Degree(), 1);
  EXPECT_EQ(poly.Size(), 2u);
  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], two);
}

TEST_F(PolynomialDenseTest, IndexingAndSetCoeffPreservesInvariant) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  EXPECT_EQ(poly[0], one);
  EXPECT_EQ(poly[1], two);
  EXPECT_EQ(poly[2], three);
  EXPECT_EQ(poly[3], zero);
  EXPECT_EQ(poly[100], zero);
  EXPECT_THROW(poly.At(3), std::out_of_range);

  poly.SetCoeff(1, four);
  EXPECT_EQ(poly[1], four);

  // Setting leading coefficient to zero must automatically trim.
  poly.SetCoeff(2, zero);
  EXPECT_EQ(poly.Degree(), 1);
  EXPECT_EQ(poly.Size(), 2u);

  // Expanding beyond current degree via SetRawCoeff.
  poly.SetRawCoeff(5, 6);
  EXPECT_EQ(poly.Degree(), 5);
  EXPECT_EQ(poly[5], six);
  EXPECT_EQ(poly[3], zero);

  // Setting beyond degree to zero is a no-op.
  poly.SetRawCoeff(10, 0);
  EXPECT_EQ(poly.Degree(), 5);
}

TEST_F(PolynomialDenseTest, Addition) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      four, five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

  auto result = p1 + p2;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], five);
  EXPECT_EQ(result[1], zero);
  EXPECT_EQ(result[2], three);
}

TEST_F(PolynomialDenseTest, AdditionWithZero) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto zero_poly = PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(gf7);

  auto result = poly + zero_poly;
  EXPECT_EQ(result, poly);
}

TEST_F(PolynomialDenseTest, Subtraction) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      four, five};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

  auto result = p1 - p2;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], four);
  EXPECT_EQ(result[1], four);
  EXPECT_EQ(result[2], three);
}

TEST_F(PolynomialDenseTest, Negation) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto neg_poly = -poly;

  EXPECT_EQ(neg_poly.Degree(), 2);
  EXPECT_EQ(neg_poly[0], six);
  EXPECT_EQ(neg_poly[1], five);
  EXPECT_EQ(neg_poly[2], four);
}

TEST_F(PolynomialDenseTest, Multiplication) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs1 = {
      one, two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1(coeffs1);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs2 = {
      three, four};
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2(coeffs2);

  auto result = p1 * p2;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], three);
  EXPECT_EQ(result[1], three);
  EXPECT_EQ(result[2], one);
}

TEST_F(PolynomialDenseTest, ScalarMultiplicationAndDivision) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto result = poly * five;
  auto left_result = five * poly;

  EXPECT_EQ(result, left_result);
  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], five);
  EXPECT_EQ(result[1], three);
  EXPECT_EQ(result[2], one);

  auto back = result / five;
  EXPECT_EQ(back, poly);

  auto inplace = result;
  inplace /= five;
  EXPECT_EQ(inplace, poly);

  EXPECT_THROW(poly / zero, std::domain_error);
  EXPECT_THROW(inplace /= zero, std::domain_error);
}

TEST_F(PolynomialDenseTest, MultiplicationByZero) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      one, two, three};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  auto result = poly * zero;

  EXPECT_TRUE(result.IsZero());
  EXPECT_EQ(result.Degree(), -1);
  EXPECT_EQ(result[0], zero);
}

TEST_F(PolynomialDenseTest, DivisionByConstant) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs = {
      two, four, six};
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);

  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> div_coeffs = {
      two};
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor(div_coeffs);

  auto result = poly / divisor;

  EXPECT_EQ(result.Degree(), 2);
  EXPECT_EQ(result[0], one);
  EXPECT_EQ(result[1], two);
  EXPECT_EQ(result[2], three);
}

TEST_F(PolynomialDenseTest, PolynomialDivisionAndNonMonicDivRem) {
  // (x^2 + 3x + 2) / (x + 1) = x + 2
  PolynomialDense<GaloisFieldPrime<uint8_t>> dividend({two, three, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor({one, one});

  auto result = dividend / divisor;
  EXPECT_EQ(result.Degree(), 1);
  EXPECT_EQ(result[0], two);
  EXPECT_EQ(result[1], one);

  // Non-monic divisor: (5x^5 + 3x^4 + 2x^3 + 6x^2 + x + 4) / (3x^2 + 2x + 5)
  auto non_monic_a = PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw(
      {4, 1, 6, 2, 3, 5}, gf7);
  auto non_monic_b =
      PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({5, 2, 3}, gf7);

  auto [q, r] = non_monic_a.DivRem(non_monic_b);
  EXPECT_LT(r.Degree(), non_monic_b.Degree());
  EXPECT_EQ(q * non_monic_b + r, non_monic_a);
  EXPECT_EQ(non_monic_a / non_monic_b, q);
  EXPECT_EQ(non_monic_a % non_monic_b, r);
}

TEST_F(PolynomialDenseTest, ModuloOperation) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> dividend({two, three, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor({one, one});

  auto remainder = dividend % divisor;

  EXPECT_TRUE(remainder.IsZero());
  EXPECT_EQ(remainder.Degree(), -1);
  EXPECT_EQ(remainder[0], zero);
}

TEST_F(PolynomialDenseTest, DivRemOperation) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> dividend({four, three, two, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> divisor({one, one});

  auto [quotient, remainder] = dividend.DivRem(divisor);

  auto verification = quotient * divisor + remainder;
  EXPECT_EQ(verification, dividend);
}

TEST_F(PolynomialDenseTest, DivisionByZero) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly({one, two});
  auto zero_poly = PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(gf7);

  EXPECT_THROW(poly.DivRem(zero_poly), std::invalid_argument);
  EXPECT_THROW(poly / zero_poly, std::invalid_argument);
  EXPECT_THROW(poly % zero_poly, std::invalid_argument);
  EXPECT_THROW(poly /= zero_poly, std::invalid_argument);
  EXPECT_THROW(poly %= zero_poly, std::invalid_argument);
}

TEST_F(PolynomialDenseTest, PowerAndPowModOperations) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly({one, one});

  auto power0 = poly ^ 0;
  EXPECT_EQ(power0.Degree(), 0);
  EXPECT_EQ(power0[0], one);
  EXPECT_EQ(poly.Pow(0), power0);

  auto power1 = poly ^ 1;
  EXPECT_EQ(power1, poly);

  auto power2 = poly ^ 2;
  EXPECT_EQ(power2.Degree(), 2);
  EXPECT_EQ(power2[0], one);
  EXPECT_EQ(power2[1], two);
  EXPECT_EQ(power2[2], one);

  // Verify PowMod against Pow() % mod.
  auto mod =
      PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({1, 0, 1}, gf7);
  for (uint64_t e : {0ULL, 1ULL, 2ULL, 5ULL, 7ULL, 13ULL, 48ULL}) {
    EXPECT_EQ(poly.PowMod(e, mod), poly.Pow(e) % mod);
  }
  auto zero_poly = PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(gf7);
  EXPECT_THROW(poly.PowMod(3, zero_poly), std::invalid_argument);
  auto const_mod = PolynomialDense<GaloisFieldPrime<uint8_t>>::One(gf7);
  EXPECT_TRUE(poly.PowMod(5, const_mod).IsZero());
}

TEST_F(PolynomialDenseTest, EvaluationAndEvalBatch) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly({one, two, three});

  EXPECT_EQ(poly(zero), one);
  EXPECT_EQ(poly(one), six);
  EXPECT_EQ(poly(two), three);

  std::vector<uint8_t> pts = {0, 1, 2, 3, 4, 5, 6};
  auto batch = poly.EvalBatch(pts);
  ASSERT_EQ(batch.size(), pts.size());
  for (size_t i = 0; i < pts.size(); ++i) {
    EXPECT_EQ(batch[i], poly.EvalRaw(pts[i]));
    EXPECT_EQ(batch[i],
              poly(GaloisFieldElement<GaloisFieldPrime<uint8_t>>(pts[i], gf7))
                  .RawValue());
  }
}

TEST_F(PolynomialDenseTest, EqualityComparison) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1({one, two, three});
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2({one, two, three});
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly3({one, two, four});

  EXPECT_TRUE(poly1 == poly2);
  EXPECT_FALSE(poly1 == poly3);
  EXPECT_TRUE(poly1 != poly3);
  EXPECT_FALSE(poly1 != poly2);
}

TEST_F(PolynomialDenseTest, EqualityWithTrimming) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly1({one, two});
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly2({one, two, zero});

  EXPECT_TRUE(poly1 == poly2);
}

TEST_F(PolynomialDenseTest, CompoundAssignmentAllOperators) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> p1({four, three, two, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> p2({five, two, three});

  auto copy = p1;
  copy += p2;
  EXPECT_EQ(copy, p1 + p2);

  copy = p1;
  copy -= p2;
  EXPECT_EQ(copy, p1 - p2);

  // Self-subtraction.
  copy = p1;
  copy -= copy;
  EXPECT_TRUE(copy.IsZero());

  copy = p1;
  copy *= p2;
  EXPECT_EQ(copy, p1 * p2);

  copy = p1;
  copy *= five;
  EXPECT_EQ(copy, p1 * five);

  copy = p1;
  copy /= p2;
  EXPECT_EQ(copy, p1 / p2);

  copy = p1;
  copy %= p2;
  EXPECT_EQ(copy, p1 % p2);

  // Self-modulo.
  copy = p1;
  copy %= copy;
  EXPECT_TRUE(copy.IsZero());
}

TEST_F(PolynomialDenseTest, ShiftReverseCompose) {
  auto p = PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({1, 2, 3}, gf7);
  auto shifted = p.Shift(3);
  EXPECT_EQ(shifted.Degree(), 5);
  EXPECT_EQ(shifted.RawCoeff(0), 0);
  EXPECT_EQ(shifted.RawCoeff(2), 0);
  EXPECT_EQ(shifted.RawCoeff(3), 1);
  EXPECT_EQ(shifted.RawCoeff(4), 2);
  EXPECT_EQ(shifted.RawCoeff(5), 3);
  EXPECT_EQ(p.Shift(0), p);

  auto rev2 = p.Reverse(2);
  EXPECT_EQ(rev2.RawCoeff(0), 3);
  EXPECT_EQ(rev2.RawCoeff(1), 2);
  EXPECT_EQ(rev2.RawCoeff(2), 1);

  auto rev4 = p.Reverse(4);
  EXPECT_EQ(rev4.Degree(), 4);
  EXPECT_EQ(rev4.RawCoeff(0), 0);
  EXPECT_EQ(rev4.RawCoeff(1), 0);
  EXPECT_EQ(rev4.RawCoeff(2), 3);
  EXPECT_EQ(rev4.RawCoeff(3), 2);
  EXPECT_EQ(rev4.RawCoeff(4), 1);
  EXPECT_THROW(p.Reverse(1), std::invalid_argument);

  // Compose p(q(t)) and verify evaluation identity p(q(x)) == (p o q)(x).
  auto q =
      PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({5, 1, 2}, gf7, "t");
  auto comp = p.Compose(q);
  EXPECT_EQ(comp.GetVariable(), "t");
  for (uint8_t x = 0; x < 7; ++x) {
    EXPECT_EQ(comp.EvalRaw(x), p.EvalRaw(q.EvalRaw(x)));
  }
}

TEST_F(PolynomialDenseTest, DerivativeCharPrimeAndChar2) {
  // GF(7): d/dx (1 + 2x + 3x^2 + 4x^3 + 5x^7 + 6x^8)
  // = 2 + 6x + 5x^2 + 0*x^6 + 6x^7
  auto poly = PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw(
      {1, 2, 3, 4, 0, 0, 0, 5, 6}, gf7);
  auto deriv = poly.Derivative();
  EXPECT_EQ(deriv.Degree(), 7);
  EXPECT_EQ(deriv[0], two);
  EXPECT_EQ(deriv[1], six);
  EXPECT_EQ(deriv[2], five);
  EXPECT_EQ(deriv[6], zero);
  EXPECT_EQ(deriv[7], six);

  // Characteristic 2: d/dx (1 + x + x^2 + x^3 + x^4 + x^5) = 1 + x^2 + x^4
  auto bpoly =
      PolynomialDense<GaloisFieldBinary>::FromRaw({1, 1, 1, 1, 1, 1}, gf2);
  auto bderiv = bpoly.Derivative();
  EXPECT_EQ(bderiv.Degree(), 4);
  EXPECT_EQ(bderiv.RawCoeff(0), 1);
  EXPECT_EQ(bderiv.RawCoeff(1), 0);
  EXPECT_EQ(bderiv.RawCoeff(2), 1);
  EXPECT_EQ(bderiv.RawCoeff(3), 0);
  EXPECT_EQ(bderiv.RawCoeff(4), 1);

  // Derivative of pure square over GF(2) is zero.
  auto bsq = PolynomialDense<GaloisFieldBinary>::FromRaw({1, 0, 1, 0, 1}, gf2);
  EXPECT_TRUE(bsq.Derivative().IsZero());
}

TEST_F(PolynomialDenseTest, DerivativeOfConstant) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly({five});
  auto derivative = poly.Derivative();

  EXPECT_TRUE(derivative.IsZero());
  EXPECT_EQ(derivative.Degree(), -1);
  EXPECT_EQ(derivative[0], zero);
}

TEST_F(PolynomialDenseTest, KaratsubaMatchesSchoolbook) {
  std::mt19937_64 rng(0x123456789ABCDEF0ULL);

  // Test across GF(7), GF(2), GF(2^4) (direct & Zech-log), and GF(3^2).
  auto gf2_4 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4);
  auto gf2_4_tbl = std::make_shared<GFBEZechLogTables<uint16_t>>(4);
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));

  auto check_field = [&](auto field,
                         const std::vector<std::pair<size_t, size_t>> &sizes) {
    using Field = typename decltype(field)::element_type;
    using Raw = typename Field::element_type;
    for (auto [na, nb] : sizes) {
      std::vector<Raw> ra(na);
      std::vector<Raw> rb(nb);
      for (size_t i = 0; i < na; ++i) ra[i] = field->Random(rng);
      for (size_t i = 0; i < nb; ++i) rb[i] = field->Random(rng);
      auto pa = PolynomialDense<Field>::FromRaw(ra, field);
      auto pb = PolynomialDense<Field>::FromRaw(rb, field);
      auto fast = pa * pb;
      auto ref = NaiveSchoolbookMul(pa, pb);
      EXPECT_EQ(fast, ref) << "Mismatch for sizes " << na << " x " << nb;
    }
  };

  const std::vector<std::pair<size_t, size_t>> test_sizes = {
      {33, 33}, {64, 64}, {65, 65}, {70, 60}, {100, 35}, {35, 100}};

  check_field(gf7, test_sizes);
  check_field(gf2, test_sizes);
  check_field(gf2_4, {{33, 33}, {65, 60}, {90, 40}});
  check_field(gf2_4_tbl, {{33, 33}, {65, 60}, {90, 40}});
  check_field(gf3_2, {{33, 33}, {64, 50}});
}

TEST_F(PolynomialDenseTest, PrintAndToStringAcrossAllFieldTypes) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly({one, two, three});
  std::stringstream ss;
  poly.Print(ss);
  EXPECT_EQ(ss.str(), "3x^2 + 2x + 1");
  EXPECT_EQ(poly.ToString(), "3x^2 + 2x + 1");

  std::stringstream ss2;
  ss2 << poly;
  EXPECT_EQ(ss2.str(), "3x^2 + 2x + 1");

  // GF(2)
  auto bpoly = PolynomialDense<GaloisFieldBinary>::FromRaw({1, 1, 0, 1}, gf2);
  EXPECT_EQ(bpoly.ToString(), "x^3 + x + 1");

  // GF(2^4) in HEX rep
  auto gf2_4 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(4, "hex");
  auto hex_poly =
      PolynomialDense<GaloisFieldBinaryExtension<uint16_t>>::FromRaw(
          {gf2_4->SetElementValue(3), gf2_4->SetElementValue(5),
           gf2_4->SetElementValue(1)},
          gf2_4);
  std::string hex_str = hex_poly.ToString();
  EXPECT_NE(hex_str.find("x^2"), std::string::npos);
  EXPECT_NE(hex_str.find("(0x5)x"), std::string::npos);

  // GF(3^2) in POLY rep (UTF-8 alpha)
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  auto c0 = gf3_2->SetElementValue("α + 1");
  auto c1 = gf3_2->SetElementValue("2α + 1");
  auto c2 = gf3_2->MultiplicativeIdentity();
  auto ext_poly = PolynomialDense<GaloisFieldExtension<uint32_t>>::FromRaw(
      {c0, c1, c2}, gf3_2);
  std::string ext_str = ext_poly.ToString();
  EXPECT_NE(ext_str.find("x^2"), std::string::npos);
  EXPECT_NE(ext_str.find("α"), std::string::npos);
}

TEST_F(PolynomialDenseTest, BinaryFieldOperations) {
  PolynomialDense<GaloisFieldBinary> poly1({gf2_one, gf2_zero, gf2_one});
  PolynomialDense<GaloisFieldBinary> poly2({gf2_one, gf2_one});

  auto sum = poly1 + poly2;
  EXPECT_EQ(sum.Degree(), 2);
  EXPECT_EQ(sum[0], gf2_zero);
  EXPECT_EQ(sum[1], gf2_one);
  EXPECT_EQ(sum[2], gf2_one);

  auto product = poly1 * poly2;
  EXPECT_EQ(product.Degree(), 3);
  EXPECT_EQ(product[0], gf2_one);
  EXPECT_EQ(product[1], gf2_one);
  EXPECT_EQ(product[2], gf2_one);
  EXPECT_EQ(product[3], gf2_one);
}

TEST_F(PolynomialDenseTest, ExtensionFieldPolynomialOperations) {
  // GF(2^4) Zech-log table
  auto gf2_4 = std::make_shared<GFBEZechLogTables<uint16_t>>(4);
  auto a = GaloisFieldElement<GFBEZechLogTables<uint16_t>>(5, gf2_4);
  auto b = GaloisFieldElement<GFBEZechLogTables<uint16_t>>(11, gf2_4);
  auto c = GaloisFieldElement<GFBEZechLogTables<uint16_t>>(1, gf2_4);
  PolynomialDense<GFBEZechLogTables<uint16_t>> p1({a, b, c});
  PolynomialDense<GFBEZechLogTables<uint16_t>> p2({b, a});
  auto prod = p1 * p2;
  auto [q, r] = prod.DivRem(p2);
  EXPECT_EQ(q, p1);
  EXPECT_TRUE(r.IsZero());
  EXPECT_EQ(prod(a), p1(a) * p2(a));

  // GF(3^2) GaloisFieldExtension
  auto gf3_2 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::make_pair(3, 2));
  using ExtElem = GaloisFieldElement<GaloisFieldExtension<uint32_t>>;
  ExtElem e0("α + 2", gf3_2);
  ExtElem e1("2α + 1", gf3_2);
  ExtElem e2("1", gf3_2);
  PolynomialDense<GaloisFieldExtension<uint32_t>> ep1({e0, e1, e2});
  PolynomialDense<GaloisFieldExtension<uint32_t>> ep2({e1, e0});
  auto eprod = ep1 * ep2 + ep2;
  auto [eq, er] = eprod.DivRem(ep1);
  EXPECT_EQ(eq * ep1 + er, eprod);
  EXPECT_LT(er.Degree(), ep1.Degree());

  // Derivative in GF(3^2): d/dx (e0 + e1*x + e2*x^2 + e1*x^3) = e1 + 2*e2*x
  PolynomialDense<GaloisFieldExtension<uint32_t>> ep3({e0, e1, e2, e1});
  auto ed = ep3.Derivative();
  EXPECT_EQ(ed.Degree(), 1);
  EXPECT_EQ(ed[0], e1);
  EXPECT_EQ(ed[1], e2 + e2);
}

TEST_F(PolynomialDenseTest, ExtendedGCD) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> a({two, three, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> b({one, one});

  auto [gcd, coeffs] = utils::PolynomialDenseExtendedGcd(a, b);
  auto [s, t] = coeffs;

  EXPECT_TRUE(gcd.Degree() >= 0);
  EXPECT_TRUE(gcd.IsMonic());
  EXPECT_EQ(a * s + b * t, gcd);
}

TEST_F(PolynomialDenseTest, ExtendedGCDCoprime) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> a({one, one});
  PolynomialDense<GaloisFieldPrime<uint8_t>> b({two, one});

  auto [gcd, coeffs] = utils::PolynomialDenseExtendedGcd(a, b);
  auto [s, t] = coeffs;

  EXPECT_EQ(gcd.Degree(), 0);
  EXPECT_EQ(gcd[0], one);
  EXPECT_EQ(a * s + b * t, gcd);
}

TEST_F(PolynomialDenseTest, EmptyAndInvalidConstructorsThrow) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> empty_coeffs;
  EXPECT_THROW((void)PolynomialDense<GaloisFieldPrime<uint8_t>>(empty_coeffs),
               std::invalid_argument);

  // Unattached default element throws.
  GaloisFieldElement<GaloisFieldPrime<uint8_t>> unattached_elem;
  EXPECT_THROW(PolynomialDense<GaloisFieldPrime<uint8_t>>({unattached_elem}),
               std::invalid_argument);

  // Mixed-field elements throw.
  auto gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);
  GaloisFieldElement<GaloisFieldPrime<uint8_t>> five_one(1, gf5);
  EXPECT_THROW(PolynomialDense<GaloisFieldPrime<uint8_t>>({one, five_one}),
               std::invalid_argument);

  // Null field pointer in static factories throws.
  std::shared_ptr<GaloisFieldPrime<uint8_t>> null_field;
  EXPECT_THROW(PolynomialDense<GaloisFieldPrime<uint8_t>>::Zero(null_field),
               std::invalid_argument);
  EXPECT_THROW(PolynomialDense<GaloisFieldPrime<uint8_t>>::One(null_field),
               std::invalid_argument);
  EXPECT_THROW(
      PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({1}, null_field),
      std::invalid_argument);
}

TEST_F(PolynomialDenseTest, FieldMismatchThrowsAndInequality) {
  auto gf5 = std::make_shared<GaloisFieldPrime<uint8_t>>(5);
  auto gf7_other = std::make_shared<GaloisFieldPrime<uint8_t>>(7);

  auto p7 = PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({1, 2}, gf7);
  auto p5 = PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({1, 2}, gf5);
  auto p7_other =
      PolynomialDense<GaloisFieldPrime<uint8_t>>::FromRaw({1, 2}, gf7_other);
  GaloisFieldElement<GaloisFieldPrime<uint8_t>> s5(2, gf5);

  EXPECT_FALSE(p7 == p5);
  EXPECT_TRUE(p7 != p5);
  EXPECT_FALSE(p7 == p7_other);

  EXPECT_THROW(p7 + p5, std::invalid_argument);
  EXPECT_THROW(p7 - p5, std::invalid_argument);
  EXPECT_THROW(p7 * p5, std::invalid_argument);
  EXPECT_THROW(p7 / p5, std::invalid_argument);
  EXPECT_THROW(p7 % p5, std::invalid_argument);
  EXPECT_THROW(p7.DivRem(p5), std::invalid_argument);
  EXPECT_THROW(p7.PowMod(2, p5), std::invalid_argument);
  EXPECT_THROW(p7.Compose(p5), std::invalid_argument);
  EXPECT_THROW(p7 * s5, std::invalid_argument);
  EXPECT_THROW(p7 / s5, std::invalid_argument);
  EXPECT_THROW(p7(s5), std::invalid_argument);

  auto copy = p7;
  EXPECT_THROW(copy += p5, std::invalid_argument);
  EXPECT_THROW(copy -= p5, std::invalid_argument);
  EXPECT_THROW(copy *= p5, std::invalid_argument);
  EXPECT_THROW(copy /= p5, std::invalid_argument);
  EXPECT_THROW(copy %= p5, std::invalid_argument);
  EXPECT_THROW(copy.SetCoeff(0, s5), std::invalid_argument);
}

TEST_F(PolynomialDenseTest, LargePolynomials) {
  std::vector<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> coeffs;
  coeffs.reserve(10);
  for (int i = 0; i < 10; ++i) {
    coeffs.push_back(
        GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(i % 7, gf7));
  }

  PolynomialDense<GaloisFieldPrime<uint8_t>> poly(coeffs);
  EXPECT_EQ(poly.Degree(), 9);
  EXPECT_EQ(poly.Size(), 10u);

  auto doubled = poly * two;
  EXPECT_EQ(doubled.Degree(), 9);

  auto squared = poly * poly;
  EXPECT_EQ(squared.Degree(), 18);
}

TEST_F(PolynomialDenseTest, FieldConsistency) {
  PolynomialDense<GaloisFieldPrime<uint8_t>> poly({one, two, three});

  EXPECT_EQ(poly.Field(), gf7);
  EXPECT_EQ((poly * two).Field(), gf7);
  EXPECT_EQ((poly + poly).Field(), gf7);
}