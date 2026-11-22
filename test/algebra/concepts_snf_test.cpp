#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/ring_array.hpp"
#include "xgalois/algebra/ring_element.hpp"
#include "xgalois/algebra/snf.hpp"
#include "xgalois/field/gf128.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"

namespace {

using xg::AbelianGroup;
using xg::AdditiveGroupView;
using xg::Algebra;
using xg::AlgebraTraits;
using xg::CommutativeRing;
using xg::ComputeSmithNormalForm;
using xg::DecomposeFinitelyGeneratedAbelianGroup;
using xg::ElementwiseMul;
using xg::EuclideanDomain;
using xg::Field;
using xg::GaloisFieldBinary;
using xg::GaloisFieldBinary128;
using xg::GaloisFieldBinaryExtension;
using xg::GaloisFieldExtension;
using xg::GaloisFieldPrime;
using xg::GF128;
using xg::GFArray;
using xg::Group;
using xg::IntegralDomain;
using xg::LeftScalarMul;
using xg::Magma;
using xg::MatMul;
using xg::Module;
using xg::Monoid;
using xg::MultiplicativeGroupView;
using xg::MultiplicativeMonoidView;
using xg::PrincipalIdealDomain;
using xg::RightScalarMul;
using xg::Ring;
using xg::RingArray;
using xg::RingElement;
using xg::Semigroup;
using xg::ToGFArray;
using xg::ToRingArray;
using xg::VectorSpace;

// ============================================================================
// Mock Non-Commutative Ring: 2x2 Integer Matrices M_2(Z)
// ============================================================================

struct Mat2 {
  int64_t a00 = 0;
  int64_t a01 = 0;
  int64_t a10 = 0;
  int64_t a11 = 0;

  constexpr bool operator==(const Mat2& other) const = default;
};

class Mat2IntRing {
 public:
  using element_type = Mat2;
  using scalar_type = int64_t;

  Mat2 Add(const Mat2& x, const Mat2& y) const {
    return {x.a00 + y.a00, x.a01 + y.a01, x.a10 + y.a10, x.a11 + y.a11};
  }

  Mat2 Sub(const Mat2& x, const Mat2& y) const {
    return {x.a00 - y.a00, x.a01 - y.a01, x.a10 - y.a10, x.a11 - y.a11};
  }

  Mat2 Neg(const Mat2& x) const { return {-x.a00, -x.a01, -x.a10, -x.a11}; }

  Mat2 Scale(int64_t s, const Mat2& x) const {
    return {s * x.a00, s * x.a01, s * x.a10, s * x.a11};
  }

  Mat2 Mul(const Mat2& x, const Mat2& y) const {
    return {x.a00 * y.a00 + x.a01 * y.a10, x.a00 * y.a01 + x.a01 * y.a11,
            x.a10 * y.a00 + x.a11 * y.a10, x.a10 * y.a01 + x.a11 * y.a11};
  }

  Mat2 AdditiveIdentity() const { return {0, 0, 0, 0}; }
  Mat2 MultiplicativeIdentity() const { return {1, 0, 0, 1}; }

  bool IsZero(const Mat2& x) const { return x == AdditiveIdentity(); }
  bool IsOne(const Mat2& x) const { return x == MultiplicativeIdentity(); }

  bool IsUnit(const Mat2& x) const {
    const int64_t det = x.a00 * x.a11 - x.a01 * x.a10;
    return det == 1 || det == -1;
  }

  std::optional<Mat2> TryInv(const Mat2& x) const {
    const int64_t det = x.a00 * x.a11 - x.a01 * x.a10;
    if (det == 1) {
      return Mat2{x.a11, -x.a01, -x.a10, x.a00};
    }
    if (det == -1) {
      return Mat2{-x.a11, x.a01, x.a10, -x.a00};
    }
    return std::nullopt;
  }

  std::string ToString(const Mat2& x) const {
    std::ostringstream oss;
    oss << "[[" << x.a00 << "," << x.a01 << "],[" << x.a10 << "," << x.a11
        << "]]";
    return oss.str();
  }
};

// ============================================================================
// Mock Euclidean Domain: Integers Z (not a Field)
// ============================================================================

class Int64EuclideanRing {
 public:
  using element_type = int64_t;

  int64_t Add(int64_t a, int64_t b) const { return a + b; }
  int64_t Sub(int64_t a, int64_t b) const { return a - b; }
  int64_t Neg(int64_t a) const { return -a; }
  int64_t Mul(int64_t a, int64_t b) const { return a * b; }
  int64_t AdditiveIdentity() const { return 0; }
  int64_t MultiplicativeIdentity() const { return 1; }
  bool IsZero(int64_t a) const { return a == 0; }
  bool IsOne(int64_t a) const { return a == 1; }
  std::string ToString(int64_t a) const { return std::to_string(a); }

  std::pair<int64_t, int64_t> DivRem(int64_t a, int64_t b) const {
    if (b == 0) {
      throw std::domain_error("division by zero");
    }
    return {a / b, a % b};
  }

  uint64_t EuclideanNorm(int64_t a) const {
    return a < 0 ? static_cast<uint64_t>(-a) : static_cast<uint64_t>(a);
  }

  std::pair<int64_t, int64_t> AssociateCanonical(int64_t a) const {
    if (a < 0) {
      return {-a, -1};
    }
    return {a, 1};
  }

  bool IsUnit(int64_t a) const { return a == 1 || a == -1; }

  std::optional<int64_t> TryInv(int64_t a) const {
    if (a == 1 || a == -1) {
      return a;
    }
    return std::nullopt;
  }
};

// ============================================================================
// Mock Vector Space / Algebra over GF(p)
// ============================================================================

struct Vec2Prime {
  uint32_t x = 0;
  uint32_t y = 0;
  constexpr bool operator==(const Vec2Prime&) const = default;
};

class PrimeVec2Space {
 public:
  using element_type = Vec2Prime;
  using scalar_type = uint32_t;

  explicit PrimeVec2Space(const GaloisFieldPrime<uint32_t>& f) : f_(f) {}

  Vec2Prime Add(const Vec2Prime& a, const Vec2Prime& b) const {
    return {f_.Add(a.x, b.x), f_.Add(a.y, b.y)};
  }
  Vec2Prime Sub(const Vec2Prime& a, const Vec2Prime& b) const {
    return {f_.Sub(a.x, b.x), f_.Sub(a.y, b.y)};
  }
  Vec2Prime Neg(const Vec2Prime& a) const { return {f_.Neg(a.x), f_.Neg(a.y)}; }
  Vec2Prime Scale(uint32_t s, const Vec2Prime& a) const {
    return {f_.Mul(s, a.x), f_.Mul(s, a.y)};
  }
  Vec2Prime AdditiveIdentity() const { return {0, 0}; }
  bool IsZero(const Vec2Prime& a) const { return a.x == 0 && a.y == 0; }

 private:
  const GaloisFieldPrime<uint32_t>& f_;
};

// Minimal custom signed integer type (mimics xg::BigInt interface without
// specializing std::is_signed or std::abs) to verify SNF genericity.
class MockBigInt {
 public:
  MockBigInt() = default;
  MockBigInt(int64_t v) : v_(v) {}  // NOLINT(runtime/explicit)

  int64_t Value() const { return v_; }

  MockBigInt operator-() const { return MockBigInt(-v_); }
  MockBigInt operator+(const MockBigInt& o) const {
    return MockBigInt(v_ + o.v_);
  }
  MockBigInt operator-(const MockBigInt& o) const {
    return MockBigInt(v_ - o.v_);
  }
  MockBigInt operator*(const MockBigInt& o) const {
    return MockBigInt(v_ * o.v_);
  }
  MockBigInt operator/(const MockBigInt& o) const {
    return MockBigInt(v_ / o.v_);
  }
  MockBigInt operator%(const MockBigInt& o) const {
    return MockBigInt(v_ % o.v_);
  }

  bool operator==(const MockBigInt& o) const { return v_ == o.v_; }
  bool operator!=(const MockBigInt& o) const { return v_ != o.v_; }
  bool operator<(const MockBigInt& o) const { return v_ < o.v_; }
  bool operator<=(const MockBigInt& o) const { return v_ <= o.v_; }
  bool operator>(const MockBigInt& o) const { return v_ > o.v_; }
  bool operator>=(const MockBigInt& o) const { return v_ >= o.v_; }

 private:
  int64_t v_ = 0;
};

// ============================================================================
// Exact matrix helpers for SNF verification
// ============================================================================

template <typename IntT>
std::vector<std::vector<IntT>> MatMulInt(
    const std::vector<std::vector<IntT>>& a,
    const std::vector<std::vector<IntT>>& b) {
  const std::size_t m = a.size();
  const std::size_t k_dim = m == 0 ? 0 : a[0].size();
  const std::size_t n = b.empty() ? 0 : b[0].size();
  std::vector<std::vector<IntT>> out(m, std::vector<IntT>(n, IntT{0}));
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      IntT acc{0};
      for (std::size_t k = 0; k < k_dim; ++k) {
        acc = acc + a[i][k] * b[k][j];
      }
      out[i][j] = acc;
    }
  }
  return out;
}

// Exact integer determinant via Bareiss fraction-free algorithm.
template <typename IntT>
IntT ExactDeterminant(std::vector<std::vector<IntT>> mat) {
  const std::size_t n = mat.size();
  if (n == 0) {
    return IntT{1};
  }
  IntT sign{1};
  IntT prev_pivot{1};
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot_row = k;
    while (pivot_row < n && mat[pivot_row][k] == IntT{0}) {
      ++pivot_row;
    }
    if (pivot_row == n) {
      return IntT{0};
    }
    if (pivot_row != k) {
      std::swap(mat[k], mat[pivot_row]);
      sign = -sign;
    }
    const IntT pivot = mat[k][k];
    for (std::size_t i = k + 1; i < n; ++i) {
      for (std::size_t j = k + 1; j < n; ++j) {
        mat[i][j] = (mat[i][j] * pivot - mat[i][k] * mat[k][j]) / prev_pivot;
      }
      mat[i][k] = IntT{0};
    }
    prev_pivot = pivot;
  }
  return sign * mat[n - 1][n - 1];
}

template <typename IntT>
void VerifySnfInvariants(const std::vector<std::vector<IntT>>& a,
                         const xg::SmithNormalFormResult<IntT>& snf) {
  const std::size_t m = a.size();
  const std::size_t n = m == 0 ? 0 : a[0].size();

  ASSERT_EQ(snf.S.size(), m);
  ASSERT_EQ(snf.U.size(), m);
  ASSERT_EQ(snf.V.size(), n);

  // Verify unimodularity: |det(U)| == 1 and |det(V)| == 1.
  const IntT det_u = ExactDeterminant(snf.U);
  const IntT det_v = ExactDeterminant(snf.V);
  EXPECT_TRUE(det_u == IntT{1} || det_u == IntT{-1});
  EXPECT_TRUE(det_v == IntT{1} || det_v == IntT{-1});

  // Verify U * A * V == S (when both dimensions are non-empty).
  if (m > 0 && n > 0) {
    const auto ua = MatMulInt(snf.U, a);
    const auto uav = MatMulInt(ua, snf.V);
    EXPECT_EQ(uav, snf.S);
  }

  // Verify diagonal structure, positivity, and divisibility chain.
  EXPECT_EQ(snf.invariant_factors.size(), snf.rank);
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      if (i == j && i < snf.rank) {
        EXPECT_GT(snf.S[i][j], IntT{0});
        EXPECT_EQ(snf.S[i][j], snf.invariant_factors[i]);
      } else {
        EXPECT_EQ(snf.S[i][j], IntT{0});
      }
    }
  }
  for (std::size_t i = 0; i + 1 < snf.invariant_factors.size(); ++i) {
    EXPECT_EQ(snf.invariant_factors[i + 1] % snf.invariant_factors[i], IntT{0});
  }
}

}  // namespace

namespace xg {

template <>
struct AlgebraTraits<Int64EuclideanRing> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = false;
};

}  // namespace xg

// ============================================================================
// Compile-Time Concept Verification (static_assert)
// ============================================================================

template <typename F>
consteval bool SatisfiesFullFieldHierarchy() {
  return Ring<F> && CommutativeRing<F> && IntegralDomain<F> &&
         PrincipalIdealDomain<F> && EuclideanDomain<F> && Field<F>;
}

static_assert(SatisfiesFullFieldHierarchy<GaloisFieldPrime<uint32_t>>());
static_assert(SatisfiesFullFieldHierarchy<GaloisFieldBinary>());
static_assert(
    SatisfiesFullFieldHierarchy<GaloisFieldBinaryExtension<uint8_t>>());
static_assert(SatisfiesFullFieldHierarchy<GaloisFieldBinary128>());
static_assert(SatisfiesFullFieldHierarchy<GF128>());
static_assert(SatisfiesFullFieldHierarchy<GaloisFieldExtension<>>());

// Bridge views satisfy exact group/monoid concepts.
static_assert(AbelianGroup<AdditiveGroupView<GaloisFieldPrime<uint32_t>>>);
static_assert(AbelianGroup<AdditiveGroupView<GaloisFieldBinary>>);
static_assert(
    AbelianGroup<AdditiveGroupView<GaloisFieldBinaryExtension<uint8_t>>>);
static_assert(AbelianGroup<AdditiveGroupView<GaloisFieldBinary128>>);
static_assert(AbelianGroup<AdditiveGroupView<GaloisFieldExtension<>>>);
static_assert(AbelianGroup<AdditiveGroupView<Mat2IntRing>>);

static_assert(Monoid<MultiplicativeMonoidView<GaloisFieldPrime<uint32_t>>>);
static_assert(!Group<MultiplicativeMonoidView<GaloisFieldPrime<uint32_t>>>);
static_assert(Monoid<MultiplicativeMonoidView<Mat2IntRing>>);
static_assert(
    !AlgebraTraits<MultiplicativeMonoidView<Mat2IntRing>>::kIsCommutativeOp);

static_assert(
    AbelianGroup<MultiplicativeGroupView<GaloisFieldPrime<uint32_t>>>);
static_assert(AbelianGroup<MultiplicativeGroupView<GaloisFieldBinary>>);
static_assert(
    AbelianGroup<MultiplicativeGroupView<GaloisFieldBinaryExtension<uint8_t>>>);
static_assert(AbelianGroup<MultiplicativeGroupView<GaloisFieldBinary128>>);
static_assert(AbelianGroup<MultiplicativeGroupView<GaloisFieldExtension<>>>);

// Non-commutative ring vs Euclidean domain vs Field boundaries.
static_assert(Ring<Mat2IntRing>);
static_assert(!CommutativeRing<Mat2IntRing>);
static_assert(!IntegralDomain<Mat2IntRing>);
static_assert(!PrincipalIdealDomain<Mat2IntRing>);
static_assert(!EuclideanDomain<Mat2IntRing>);
static_assert(!Field<Mat2IntRing>);

static_assert(!AlgebraTraits<Int64EuclideanRing>::kIsAssociativeOp);
static_assert(!AlgebraTraits<Int64EuclideanRing>::kIsCommutativeOp);
static_assert(Ring<Int64EuclideanRing>);
static_assert(CommutativeRing<Int64EuclideanRing>);
static_assert(IntegralDomain<Int64EuclideanRing>);
static_assert(PrincipalIdealDomain<Int64EuclideanRing>);
static_assert(EuclideanDomain<Int64EuclideanRing>);
static_assert(!Field<Int64EuclideanRing>);

// Module, VectorSpace & Algebra concepts.
static_assert(Module<PrimeVec2Space, GaloisFieldPrime<uint32_t>>);
static_assert(VectorSpace<PrimeVec2Space, GaloisFieldPrime<uint32_t>>);
static_assert(Algebra<Mat2IntRing, Int64EuclideanRing>);
static_assert(!Algebra<PrimeVec2Space, GaloisFieldPrime<uint32_t>>);

// Conditional division on RingElement<R>.
template <typename Elem>
concept HasDivOps = requires(Elem a, Elem b) {
  a / b;
  a /= b;
};

static_assert(HasDivOps<RingElement<GaloisFieldPrime<uint32_t>>>);
static_assert(!HasDivOps<RingElement<Mat2IntRing>>);
static_assert(!HasDivOps<RingElement<Int64EuclideanRing>>);

// ============================================================================
// Runtime Tests: CPOs & Bridge Views
// ============================================================================

TEST(AlgebraConceptsTest, FieldCpoDefaultsAndCustomEuclideanDomain) {
  GaloisFieldPrime<uint32_t> gf7(7);

  auto [q, r] = xg::alg::DivRem(gf7, 5u, 3u);
  EXPECT_EQ(q, gf7.Div(5u, 3u));
  EXPECT_EQ(r, 0u);

  EXPECT_EQ(xg::alg::EuclideanNorm(gf7, 0u), 0ULL);
  EXPECT_EQ(xg::alg::EuclideanNorm(gf7, 4u), 1ULL);

  auto [c0, u0] = xg::alg::AssociateCanonical(gf7, 0u);
  EXPECT_EQ(c0, 0u);
  EXPECT_EQ(u0, 1u);

  auto [c4, u4] = xg::alg::AssociateCanonical(gf7, 4u);
  EXPECT_EQ(c4, 1u);
  EXPECT_EQ(gf7.Mul(u4, 4u), c4);

  EXPECT_FALSE(xg::alg::IsUnit(gf7, 0u));
  EXPECT_TRUE(xg::alg::IsUnit(gf7, 4u));
  EXPECT_FALSE(xg::alg::TryInv(gf7, 0u).has_value());
  ASSERT_TRUE(xg::alg::TryInv(gf7, 4u).has_value());
  EXPECT_EQ(*xg::alg::TryInv(gf7, 4u), 2u);

  Int64EuclideanRing z;
  auto [zq, zr] = xg::alg::DivRem(z, int64_t{17}, int64_t{5});
  EXPECT_EQ(zq, 3);
  EXPECT_EQ(zr, 2);
  EXPECT_EQ(xg::alg::EuclideanNorm(z, int64_t{-9}), 9ULL);
  auto [zc, zu] = xg::alg::AssociateCanonical(z, int64_t{-9});
  EXPECT_EQ(zc, 9);
  EXPECT_EQ(zu, -1);
  EXPECT_EQ(z.Mul(zu, int64_t{-9}), zc);
  EXPECT_TRUE(xg::alg::IsUnit(z, int64_t{-1}));
  EXPECT_FALSE(xg::alg::IsUnit(z, int64_t{2}));
}

TEST(AlgebraConceptsTest, BridgeViews) {
  GaloisFieldPrime<uint32_t> gf11(11);
  AdditiveGroupView add_view(gf11);
  EXPECT_TRUE(add_view.IsIdentity(add_view.Identity()));
  EXPECT_EQ(add_view.Op(7u, 8u), 4u);
  EXPECT_EQ(add_view.Inv(3u), 8u);
  EXPECT_EQ(add_view.Div(3u, 5u), 9u);
  EXPECT_EQ(add_view.Pow(3u, 0), 0u);
  EXPECT_EQ(add_view.Pow(3u, 5), 4u);  // 5 * 3 mod 11 == 4
  EXPECT_TRUE(add_view.Equals(add_view.Pow(4u, 11), 0u));

  MultiplicativeMonoidView mon_view(gf11);
  EXPECT_TRUE(mon_view.IsIdentity(mon_view.Identity()));
  EXPECT_EQ(mon_view.Op(3u, 4u), 1u);
  EXPECT_EQ(mon_view.Pow(2u, 10), 1u);
  EXPECT_EQ(mon_view.Pow(0u, 3), 0u);

  MultiplicativeGroupView mul_view(gf11);
  EXPECT_TRUE(mul_view.IsIdentity(mul_view.Identity()));
  EXPECT_EQ(mul_view.Op(3u, 4u), 1u);
  EXPECT_EQ(mul_view.Inv(3u), 4u);
  EXPECT_EQ(mul_view.Div(5u, 3u), 9u);
  EXPECT_EQ(mul_view.Pow(2u, 10), 1u);
  EXPECT_THROW(mul_view.Op(0u, 3u), std::domain_error);
  EXPECT_THROW(mul_view.Inv(0u), std::domain_error);
  EXPECT_THROW(mul_view.Pow(0u, 2), std::domain_error);
}

// ============================================================================
// Runtime Tests: RingElement
// ============================================================================

TEST(RingElementTest, FieldArithmeticAndUnits) {
  auto gf7 = std::make_shared<GaloisFieldPrime<uint32_t>>(7);
  auto gf7_other = std::make_shared<GaloisFieldPrime<uint32_t>>(7);

  EXPECT_THROW(RingElement<GaloisFieldPrime<uint32_t>>(3u, nullptr),
               std::invalid_argument);

  RingElement<GaloisFieldPrime<uint32_t>> a(3u, gf7);
  RingElement<GaloisFieldPrime<uint32_t>> b(5u, gf7);
  auto z = RingElement<GaloisFieldPrime<uint32_t>>::Zero(gf7);
  auto one = RingElement<GaloisFieldPrime<uint32_t>>::One(gf7);

  EXPECT_EQ(a.Ring(), gf7);
  EXPECT_EQ(a.RawValue(), 3u);
  EXPECT_TRUE(z.IsZero());
  EXPECT_TRUE(one.IsOne());
  EXPECT_FALSE(z.IsUnit());
  EXPECT_TRUE(a.IsUnit());
  EXPECT_FALSE(z.TryInv().has_value());
  ASSERT_TRUE(a.TryInv().has_value());
  EXPECT_EQ((*a.TryInv() * a), one);

  EXPECT_EQ((a + b).RawValue(), 1u);
  EXPECT_EQ((a - b).RawValue(), 5u);
  EXPECT_EQ((a * b).RawValue(), 1u);
  EXPECT_EQ((-a).RawValue(), 4u);
  EXPECT_EQ((a / b).RawValue(), 2u);
  EXPECT_EQ(a.Pow(6), one);
  EXPECT_EQ(a.Pow(0), one);

  RingElement<GaloisFieldPrime<uint32_t>> c = a;
  c += b;
  EXPECT_EQ(c.RawValue(), 1u);
  c -= b;
  EXPECT_EQ(c, a);
  c *= b;
  EXPECT_EQ(c, one);
  c /= b;
  EXPECT_EQ(c, a);

  std::ostringstream oss;
  oss << a;
  EXPECT_EQ(oss.str(), a.ToString());

  RingElement<GaloisFieldPrime<uint32_t>> other_a(3u, gf7_other);
  EXPECT_NE(a, other_a);
  EXPECT_THROW(static_cast<void>(a + other_a), std::invalid_argument);
}

TEST(RingElementTest, NonCommutativeMatrixRing) {
  auto ring = std::make_shared<Mat2IntRing>();
  RingElement<Mat2IntRing> a(Mat2{1, 2, 0, 1}, ring);
  RingElement<Mat2IntRing> b(Mat2{1, 0, 3, 1}, ring);

  const auto ab = a * b;
  const auto ba = b * a;
  EXPECT_NE(ab, ba);
  EXPECT_EQ(ab.RawValue(), (Mat2{7, 2, 3, 1}));
  EXPECT_EQ(ba.RawValue(), (Mat2{1, 2, 3, 7}));

  EXPECT_TRUE(a.IsUnit());
  ASSERT_TRUE(a.TryInv().has_value());
  EXPECT_EQ(a * (*a.TryInv()), RingElement<Mat2IntRing>::One(ring));

  RingElement<Mat2IntRing> singular(Mat2{2, 0, 0, 0}, ring);
  EXPECT_FALSE(singular.IsUnit());
  EXPECT_FALSE(singular.TryInv().has_value());

  EXPECT_EQ(a.Pow(3).RawValue(), (Mat2{1, 6, 0, 1}));
}

// ============================================================================
// Runtime Tests: RingArray
// ============================================================================

TEST(RingArrayTest, NonCommutativeScalarAndMatrixOps) {
  auto ring = std::make_shared<Mat2IntRing>();

  const Mat2 s{1, 2, 0, 1};
  const Mat2 m0{1, 0, 3, 1};
  const Mat2 m1{0, 1, 1, 0};

  auto arr = RingArray<Mat2IntRing>::FromVector(ring, {m0, m1});
  EXPECT_EQ(arr.Ring(), ring);
  EXPECT_TRUE(arr.HasRing());
  EXPECT_EQ(arr.Dimension(), 1u);
  EXPECT_EQ(arr.Size(), 2u);
  EXPECT_EQ(arr.Shape(0), 2u);
  EXPECT_FALSE(arr.Empty());

  // Left vs Right scalar multiplication must differ on non-commutative ring!
  auto left = LeftScalarMul(s, arr);
  auto right = RightScalarMul(arr, s);
  EXPECT_NE(left, right);
  EXPECT_EQ(left(0), ring->Mul(s, m0));
  EXPECT_EQ(left(1), ring->Mul(s, m1));
  EXPECT_EQ(right(0), ring->Mul(m0, s));
  EXPECT_EQ(right(1), ring->Mul(m1, s));

  // Wrapped RingElement scalar overloads.
  RingElement<Mat2IntRing> s_elem(s, ring);
  EXPECT_EQ(LeftScalarMul(s_elem, arr), left);
  EXPECT_EQ(RightScalarMul(arr, s_elem), right);

  // Matrix operations: MatMul, Transpose, Eye, +, -, ElementwiseMul.
  auto mat_a = RingArray<Mat2IntRing>::FromMatrix(
      ring, {{s, m0}, {m1, ring->AdditiveIdentity()}});
  auto mat_b = RingArray<Mat2IntRing>::FromMatrix(
      ring, {{m0, m1}, {s, ring->MultiplicativeIdentity()}});

  auto prod_ab = MatMul(mat_a, mat_b);
  auto prod_ba = MatMul(mat_b, mat_a);
  EXPECT_NE(prod_ab, prod_ba);
  EXPECT_EQ(prod_ab(0, 0), ring->Add(ring->Mul(s, m0), ring->Mul(m0, s)));
  EXPECT_EQ(prod_ab(0, 1), ring->Add(ring->Mul(s, m1), m0));
  EXPECT_EQ(prod_ab(1, 0), ring->Mul(m1, m0));
  EXPECT_EQ(prod_ab(1, 1), ring->Mul(m1, m1));

  auto eye = RingArray<Mat2IntRing>::Eye(ring, 2);
  EXPECT_EQ(MatMul(mat_a, eye), mat_a);
  EXPECT_EQ(MatMul(eye, mat_a), mat_a);

  auto tr = mat_a.Transpose();
  EXPECT_EQ(tr.Shape(0), 2u);
  EXPECT_EQ(tr.Shape(1), 2u);
  EXPECT_EQ(tr(0, 1), mat_a(1, 0));
  EXPECT_EQ(tr(1, 0), mat_a(0, 1));
  EXPECT_EQ(tr.Transpose(), mat_a);

  auto zeros = RingArray<Mat2IntRing>::Zeros(ring, {2, 2});
  EXPECT_EQ(mat_a + (-mat_a), zeros);
  EXPECT_EQ(mat_a - mat_a, zeros);

  auto ew = ElementwiseMul(mat_a, mat_b);
  EXPECT_EQ(ew(0, 0), ring->Mul(s, m0));
  EXPECT_EQ(ew(1, 1), ring->AdditiveIdentity());
}

TEST(RingArrayTest, GFArrayZeroCopyBridgeRoundTrip) {
  auto gf7 = std::make_shared<GaloisFieldPrime<uint32_t>>(7);
  auto gf_arr = GFArray<GaloisFieldPrime<uint32_t>>::Zeros(gf7, {2, 2});
  gf_arr.Raw(0, 0) = 1u;
  gf_arr.Raw(0, 1) = 2u;
  gf_arr.Raw(1, 0) = 3u;
  gf_arr.Raw(1, 1) = 4u;

  auto ring_arr = ToRingArray(gf_arr);
  EXPECT_EQ(ring_arr.Shape(0), 2u);
  EXPECT_EQ(ring_arr.Shape(1), 2u);
  EXPECT_EQ(ring_arr(0, 0), 1u);
  EXPECT_EQ(ring_arr(1, 1), 4u);

  auto doubled = ring_arr + ring_arr;
  auto back_gf = ToGFArray(std::move(doubled));
  EXPECT_EQ(back_gf.Raw(0, 0), 2u);
  EXPECT_EQ(back_gf.Raw(0, 1), 4u);
  EXPECT_EQ(back_gf.Raw(1, 0), 6u);
  EXPECT_EQ(back_gf.Raw(1, 1), 1u);

  // Rvalue move ToRingArray.
  auto moved_ring = ToRingArray(std::move(back_gf));
  EXPECT_EQ(moved_ring(1, 1), 1u);
}

// ============================================================================
// Runtime Tests: Smith Normal Form & Abelian Group Decomposition
// ============================================================================

TEST(SmithNormalFormTest, EmptyZeroRectangularRankDeficientNegative) {
  // Empty matrix.
  {
    std::vector<std::vector<int64_t>> empty;
    auto snf = ComputeSmithNormalForm(empty);
    EXPECT_EQ(snf.rank, 0u);
    EXPECT_TRUE(snf.invariant_factors.empty());
    EXPECT_TRUE(snf.S.empty());
  }

  // Zero matrix.
  {
    std::vector<std::vector<int64_t>> z(3, std::vector<int64_t>(4, 0));
    auto snf = ComputeSmithNormalForm(z);
    EXPECT_EQ(snf.rank, 0u);
    EXPECT_TRUE(snf.invariant_factors.empty());
    VerifySnfInvariants(z, snf);
  }

  // Rectangular 2 x 4.
  {
    std::vector<std::vector<int64_t>> a = {
        {2, 4, 6, 8},
        {0, 6, 12, 18},
    };
    auto snf = ComputeSmithNormalForm(a);
    EXPECT_EQ(snf.rank, 2u);
    EXPECT_EQ(snf.invariant_factors, (std::vector<int64_t>{2, 6}));
    VerifySnfInvariants(a, snf);
  }

  // Rectangular 4 x 2.
  {
    std::vector<std::vector<int64_t>> a = {
        {6, 0},
        {0, 10},
        {6, 10},
        {12, 20},
    };
    auto snf = ComputeSmithNormalForm(a);
    EXPECT_EQ(snf.rank, 2u);
    EXPECT_EQ(snf.invariant_factors, (std::vector<int64_t>{2, 30}));
    VerifySnfInvariants(a, snf);
  }

  // Rank-deficient with negative entries.
  {
    std::vector<std::vector<int64_t>> a = {
        {-6, 12, -18},
        {4, -8, 12},
        {-2, 4, -6},
    };
    auto snf = ComputeSmithNormalForm(a);
    EXPECT_EQ(snf.rank, 1u);
    EXPECT_EQ(snf.invariant_factors, (std::vector<int64_t>{2}));
    VerifySnfInvariants(a, snf);
  }

  // Classic 3x3 integer matrix with negatives: invariant factors {2, 6, 12}.
  {
    std::vector<std::vector<int64_t>> a = {
        {2, 4, 4},
        {-6, 6, 12},
        {10, -4, -16},
    };
    auto snf = ComputeSmithNormalForm(a);
    EXPECT_EQ(snf.rank, 3u);
    EXPECT_EQ(snf.invariant_factors, (std::vector<int64_t>{2, 6, 12}));
    VerifySnfInvariants(a, snf);
  }

  // Diagonal non-divisible entries: diag(2, 4, 6) -> diag(2, 2, 12).
  {
    std::vector<std::vector<int64_t>> a = {
        {2, 0, 0},
        {0, 4, 0},
        {0, 0, 6},
    };
    auto snf = ComputeSmithNormalForm(a);
    EXPECT_EQ(snf.rank, 3u);
    EXPECT_EQ(snf.invariant_factors, (std::vector<int64_t>{2, 2, 12}));
    VerifySnfInvariants(a, snf);
  }

  // Ragged input throws.
  {
    std::vector<std::vector<int64_t>> ragged = {{1, 2}, {3}};
    EXPECT_THROW(ComputeSmithNormalForm(ragged), std::invalid_argument);
  }
}

TEST(SmithNormalFormTest, GenericInt128AndCustomBigInt) {
  __extension__ using i128 = __int128;

  std::vector<std::vector<i128>> a128 = {
      {2, 4, 4},
      {-6, 6, 12},
      {10, -4, -16},
  };
  auto snf128 = ComputeSmithNormalForm(a128);
  EXPECT_EQ(snf128.rank, 3u);
  ASSERT_EQ(snf128.invariant_factors.size(), 3u);
  EXPECT_EQ(snf128.invariant_factors[0], i128{2});
  EXPECT_EQ(snf128.invariant_factors[1], i128{6});
  EXPECT_EQ(snf128.invariant_factors[2], i128{12});
  VerifySnfInvariants(a128, snf128);

  std::vector<std::vector<MockBigInt>> a_big = {
      {MockBigInt{2}, MockBigInt{0}, MockBigInt{0}},
      {MockBigInt{0}, MockBigInt{4}, MockBigInt{0}},
      {MockBigInt{0}, MockBigInt{0}, MockBigInt{6}},
  };
  auto snf_big = ComputeSmithNormalForm(a_big);
  EXPECT_EQ(snf_big.rank, 3u);
  ASSERT_EQ(snf_big.invariant_factors.size(), 3u);
  EXPECT_EQ(snf_big.invariant_factors[0], MockBigInt{2});
  EXPECT_EQ(snf_big.invariant_factors[1], MockBigInt{2});
  EXPECT_EQ(snf_big.invariant_factors[2], MockBigInt{12});
  VerifySnfInvariants(a_big, snf_big);

  auto decomp_big = DecomposeFinitelyGeneratedAbelianGroup(a_big);
  EXPECT_EQ(decomp_big.free_rank, 0u);
  EXPECT_EQ(
      decomp_big.invariant_factors,
      (std::vector<MockBigInt>{MockBigInt{2}, MockBigInt{2}, MockBigInt{12}}));
}

TEST(AbelianGroupDecompositionTest, TorsionFreeAndTrivialGroups) {
  // Z^3 / <(2,0,0), (0,4,0), (0,0,6)> ~= Z_2 (+) Z_2 (+) Z_12
  // Elementary divisors: 2^1, 2^1, 2^2, 3^1.
  {
    std::vector<std::vector<int64_t>> rels = {
        {2, 0, 0},
        {0, 4, 0},
        {0, 0, 6},
    };
    auto res = DecomposeFinitelyGeneratedAbelianGroup(rels);
    EXPECT_EQ(res.free_rank, 0u);
    EXPECT_EQ(res.invariant_factors, (std::vector<int64_t>{2, 2, 12}));
    const std::vector<std::pair<int64_t, uint64_t>> expected_ed = {
        {2, 1}, {2, 1}, {2, 2}, {3, 1}};
    EXPECT_EQ(res.elementary_divisors, expected_ed);
  }

  // Free summands + torsion: Z^4 / <(2, 0, 0, 0), (0, 6, 0, 0)> ~= Z^2 (+) Z_2
  // (+) Z_6.
  {
    std::vector<std::vector<int64_t>> rels = {
        {2, 0, 0, 0},
        {0, 6, 0, 0},
    };
    auto res = DecomposeFinitelyGeneratedAbelianGroup(rels);
    EXPECT_EQ(res.free_rank, 2u);
    EXPECT_EQ(res.invariant_factors, (std::vector<int64_t>{2, 6}));
    const std::vector<std::pair<int64_t, uint64_t>> expected_ed = {
        {2, 1}, {2, 1}, {3, 1}};
    EXPECT_EQ(res.elementary_divisors, expected_ed);
  }

  // Trivial group via unimodular relation matrix: det = 2*2 - 1*3 = 1.
  {
    std::vector<std::vector<int64_t>> rels = {
        {2, 1},
        {3, 2},
    };
    auto res = DecomposeFinitelyGeneratedAbelianGroup(rels);
    EXPECT_EQ(res.free_rank, 0u);
    EXPECT_TRUE(res.invariant_factors.empty());
    EXPECT_TRUE(res.elementary_divisors.empty());
  }

  // Purely free group Z^3 (single zero relation on 3 generators).
  {
    std::vector<std::vector<int64_t>> rels = {{0, 0, 0}};
    auto res = DecomposeFinitelyGeneratedAbelianGroup(rels);
    EXPECT_EQ(res.free_rank, 3u);
    EXPECT_TRUE(res.invariant_factors.empty());
    EXPECT_TRUE(res.elementary_divisors.empty());
  }
}
