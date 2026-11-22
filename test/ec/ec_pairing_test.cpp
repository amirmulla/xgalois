#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/group_algorithms.hpp"
#include "xgalois/algebra/group_base.hpp"
#include "xgalois/algebra/matrix_group.hpp"
#include "xgalois/ec/elliptic_curve.hpp"
#include "xgalois/ec/pairing.hpp"
#include "xgalois/field/gf_bigint.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/field/gf_tower.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/utils/bigint.hpp"

namespace {

using xg::BigInt;
using xg::BigUInt;
using xg::FieldTower;
using xg::GaloisFieldBinary;
using xg::GaloisFieldBinaryExtension;
using xg::GaloisFieldPrime;
using xg::GaloisFieldPrimeBig;
using xg::GroupElement;
using xg::MatrixGroup;
using xg::MatrixGroupKind;
using xg::ec::AffinePoint;
using xg::ec::CurveForm;
using xg::ec::EllipticCurve;
using xg::ec::EllipticRiemannRochMatrix;
using xg::ec::JacobianPoint;
using xg::ec::TatePairing;
using xg::ec::WeilPairing;

// ============================================================================
// Concept Static Assertions
// ============================================================================

static_assert(xg::AbelianGroup<EllipticCurve<GaloisFieldPrime<uint32_t>,
                                             CurveForm::kShortWeierstrass>>);
static_assert(xg::AbelianGroup<EllipticCurve<GaloisFieldPrimeBig<256>,
                                             CurveForm::kShortWeierstrass>>);
static_assert(
    xg::AbelianGroup<EllipticCurve<GaloisFieldBinaryExtension<uint32_t>,
                                   CurveForm::kBinaryWeierstrass>>);
static_assert(xg::AbelianGroup<
              EllipticCurve<GaloisFieldPrimeBig<256>, CurveForm::kMontgomery>>);
static_assert(xg::AbelianGroup<EllipticCurve<GaloisFieldPrimeBig<256>,
                                             CurveForm::kTwistedEdwards>>);

static_assert(xg::Group<MatrixGroup<GaloisFieldPrime<uint32_t>>>);
static_assert(xg::FiniteGroup<MatrixGroup<GaloisFieldPrime<uint32_t>>>);
static_assert(xg::Group<MatrixGroup<GaloisFieldBinary>>);
static_assert(xg::FiniteGroup<MatrixGroup<GaloisFieldBinary>>);

// ============================================================================
// Short Weierstrass over GF(97)
// ============================================================================

TEST(EllipticCurveTest, ShortWeierstrassGF97ArithmeticAndDivisionPoly) {
  GaloisFieldPrime<uint32_t> gf97(97);
  // y^2 = x^3 + 2x + 3 over GF(97)
  EllipticCurve<GaloisFieldPrime<uint32_t>> curve(gf97, 2U, 3U);

  EXPECT_NE(curve.Discriminant(), 0U);
  // 4*a^3 + 27*b^2 = 32 + 243 = 275 == 81 (mod 97).
  // Delta = -16 * 81 = -1296 == 62 (mod 97).
  EXPECT_EQ(curve.Discriminant(), 62U);

  // P = (3, 6): 6^2 = 36, 3^3 + 2*3 + 3 = 27 + 6 + 3 = 36.
  AffinePoint<GaloisFieldPrime<uint32_t>> p(3U, 6U);
  EXPECT_TRUE(curve.IsOnCurve(p));
  EXPECT_FALSE(curve.IsIdentity(p));

  const auto p2 = curve.Double(p);
  const auto p2_add = curve.Add(p, p);
  EXPECT_TRUE(curve.IsOnCurve(p2));
  EXPECT_TRUE(curve.Equals(p2, p2_add));
  EXPECT_TRUE(curve.Equals(curve.Add(p, curve.Neg(p)), curve.Identity()));

  // Verify scalar mul (w-NAF, constant-time, negative scalars) vs naive.
  AffinePoint<GaloisFieldPrime<uint32_t>> naive = curve.Identity();
  for (int k = 0; k <= 30; ++k) {
    const auto wnaf = curve.ScalarMul(p, BigInt(k));
    const auto ct = curve.ScalarMulConstantTime(p, BigInt(k));
    EXPECT_TRUE(curve.Equals(wnaf, naive));
    EXPECT_TRUE(curve.Equals(ct, naive));

    const auto neg_wnaf = curve.ScalarMul(p, BigInt(-k));
    EXPECT_TRUE(curve.Equals(neg_wnaf, curve.Neg(naive)));
    naive = curve.Add(naive, p);
  }

  // Division polynomial psi_n vanishes at n-torsion points.
  const BigInt ord_p = curve.PointOrder(p);
  EXPECT_TRUE(curve.IsIdentity(curve.ScalarMul(p, ord_p)));
  EXPECT_EQ(curve.DivisionPolynomial(static_cast<uint32_t>(ord_p.ToUint64()),
                                     p.x, p.y),
            0U);
  EXPECT_NE(curve.DivisionPolynomial(1U, p.x, p.y), 0U);
  EXPECT_EQ(curve.DivisionPolynomial(2U, p.x, p.y), gf97.Add(p.y, p.y));
}

// ============================================================================
// Short Weierstrass over 256-bit Prime Field: secp256k1
// ============================================================================

TEST(EllipticCurveTest, Secp256k1GeneratorBatchAffineAndPippengerMsm) {
  using F256 = GaloisFieldPrimeBig<256>;
  const F256 fp = xg::primes::Secp256k1();
  EllipticCurve<F256> secp(fp, BigUInt<256>(0ULL), BigUInt<256>(7ULL));

  const AffinePoint<F256> g(
      BigUInt<256>::FromHex(
          "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798"),
      BigUInt<256>::FromHex("0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a6855419"
                            "9c47d08ffb10d4b8"));
  EXPECT_TRUE(secp.IsOnCurve(g));

  // Verify 2G and 3G.
  const auto g2 = secp.Double(g);
  const auto g3 = secp.Add(g2, g);
  EXPECT_TRUE(secp.IsOnCurve(g2));
  EXPECT_TRUE(secp.IsOnCurve(g3));
  EXPECT_TRUE(secp.Equals(secp.ScalarMul(g, BigInt(2)), g2));
  EXPECT_TRUE(secp.Equals(secp.ScalarMul(g, BigInt(3)), g3));

  // Verify w-NAF == constant-time ladder == naive on a 128-bit scalar.
  const BigInt k = BigInt::FromHex("0x123456789abcdef0fedcba9876543211");
  const auto p_wnaf = secp.ScalarMul(g, k);
  const auto p_ct = secp.ScalarMulConstantTime(g, k);
  EXPECT_TRUE(secp.IsOnCurve(p_wnaf));
  EXPECT_TRUE(secp.Equals(p_wnaf, p_ct));

  // Verify [n]G == O using secp256k1's exact 256-bit group order n.
  const BigInt n = BigInt::FromHex(
      "0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141");
  EXPECT_TRUE(secp.IsIdentity(secp.ScalarMul(g, n)));
  EXPECT_TRUE(secp.Equals(secp.ScalarMul(g, n - BigInt(1)), secp.Neg(g)));

  // Verify BatchToAffine (Montgomery's trick) against individual ToAffine.
  std::vector<JacobianPoint<F256>> jac_pts;
  jac_pts.push_back(JacobianPoint<F256>::Identity());
  auto cur_j = secp.ToJacobian(g);
  for (int i = 0; i < 20; ++i) {
    jac_pts.push_back(cur_j);
    cur_j = secp.AddMixedJacobian(secp.DoubleJacobian(cur_j), g);
  }
  const auto batch_aff = secp.BatchToAffine(jac_pts);
  ASSERT_EQ(batch_aff.size(), jac_pts.size());
  for (std::size_t i = 0; i < jac_pts.size(); ++i) {
    EXPECT_TRUE(secp.Equals(batch_aff[i], secp.ToAffine(jac_pts[i])));
  }

  // Verify MultiScalarMul (both Straus N < 16 and Pippenger N >= 16) vs naive.
  constexpr std::size_t kCount = 24;
  std::vector<AffinePoint<F256>> pts;
  std::vector<BigInt> scalars;
  pts.reserve(kCount);
  scalars.reserve(kCount);
  AffinePoint<F256> cur_pt = g;
  AffinePoint<F256> naive_sum = secp.Identity();
  for (std::size_t i = 0; i < kCount; ++i) {
    pts.push_back(cur_pt);
    const BigInt sc = (i % 5 == 0) ? -BigInt(static_cast<int64_t>(i * 17 + 3))
                                   : BigInt::FromUint64(i * 1009ULL + 7ULL);
    scalars.push_back(sc);
    naive_sum = secp.Add(naive_sum, secp.ScalarMul(cur_pt, sc));
    cur_pt = secp.Double(cur_pt);
  }
  const auto msm_pippenger = secp.MultiScalarMul(pts, scalars);
  EXPECT_TRUE(secp.Equals(msm_pippenger, naive_sum));

  // Small N (Straus/Shamir path).
  const auto msm_small =
      secp.MultiScalarMul(std::span<const AffinePoint<F256>>(pts.data(), 6),
                          std::span<const BigInt>(scalars.data(), 6));
  AffinePoint<F256> naive_small = secp.Identity();
  for (std::size_t i = 0; i < 6; ++i) {
    naive_small = secp.Add(naive_small, secp.ScalarMul(pts[i], scalars[i]));
  }
  EXPECT_TRUE(secp.Equals(msm_small, naive_small));
}

// ============================================================================
// Binary Weierstrass over GF(2^8)
// ============================================================================

TEST(EllipticCurveTest, BinaryWeierstrassGF256) {
  GaloisFieldBinaryExtension<uint32_t> gf256(8);
  // y^2 + xy = x^3 + x^2 + 1 over GF(2^8)
  EllipticCurve<GaloisFieldBinaryExtension<uint32_t>,
                CurveForm::kBinaryWeierstrass>
      curve(gf256, 1U, 1U);

  // Find a few distinct affine points on the curve.
  std::vector<AffinePoint<GaloisFieldBinaryExtension<uint32_t>>> pts;
  for (uint32_t x = 1; x < 256 && pts.size() < 4; ++x) {
    for (uint32_t y = 0; y < 256; ++y) {
      AffinePoint<GaloisFieldBinaryExtension<uint32_t>> cand(x, y);
      if (curve.IsOnCurve(cand)) {
        pts.push_back(cand);
        break;
      }
    }
  }
  ASSERT_GE(pts.size(), 3U);

  const auto& p = pts[0];
  const auto& q = pts[1];
  const auto& r = pts[2];

  // Negation and identity: P + (-P) == O.
  EXPECT_TRUE(curve.IsOnCurve(curve.Neg(p)));
  EXPECT_TRUE(curve.Equals(curve.Add(p, curve.Neg(p)), curve.Identity()));

  // Associativity: (P + Q) + R == P + (Q + R).
  const auto lhs = curve.Add(curve.Add(p, q), r);
  const auto rhs = curve.Add(p, curve.Add(q, r));
  EXPECT_TRUE(curve.IsOnCurve(lhs));
  EXPECT_TRUE(curve.Equals(lhs, rhs));

  // Scalar multiplication (w-NAF, constant-time, naive).
  AffinePoint<GaloisFieldBinaryExtension<uint32_t>> naive = curve.Identity();
  for (int k = 0; k <= 25; ++k) {
    EXPECT_TRUE(curve.Equals(curve.ScalarMul(p, BigInt(k)), naive));
    EXPECT_TRUE(curve.Equals(curve.ScalarMulConstantTime(p, BigInt(k)), naive));
    naive = curve.Add(naive, p);
  }
}

// ============================================================================
// Montgomery Curve (Curve25519) & Twisted Edwards Curve (Ed25519)
// ============================================================================

TEST(EllipticCurveTest, Curve25519MontgomeryAndEd25519TwistedEdwards) {
  using F256 = GaloisFieldPrimeBig<256>;
  const F256 fp = xg::primes::Curve25519();

  // Subgroup prime order l = 2^252 + 27742317777372353535851937790883648493.
  const BigInt l = BigInt::FromHex(
      "0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed");

  // --- Curve25519: y^2 = x^3 + 486662*x^2 + x ---
  EllipticCurve<F256, CurveForm::kMontgomery> mont(fp, BigUInt<256>(486662ULL),
                                                   BigUInt<256>(1ULL));
  const BigUInt<256> base_u(9ULL);
  const auto rhs_u =
      fp.Add(fp.Add(fp.Mul(fp.Sqr(base_u), base_u),
                    fp.Mul(BigUInt<256>(486662ULL), fp.Sqr(base_u))),
             base_u);
  const auto base_v_opt = fp.Sqrt(rhs_u);
  ASSERT_TRUE(base_v_opt.has_value());
  const AffinePoint<F256> base_m(base_u, *base_v_opt);
  EXPECT_TRUE(mont.IsOnCurve(base_m));

  // Verify LadderMulX matches ScalarMul's x-coordinate and naive addition.
  AffinePoint<F256> naive_m = mont.Identity();
  for (uint64_t k = 1; k <= 12; ++k) {
    naive_m = mont.Add(naive_m, base_m);
    const auto sm = mont.ScalarMul(base_m, BigInt::FromUint64(k));
    const auto lx = mont.LadderMulX(base_u, BigInt::FromUint64(k));
    EXPECT_TRUE(mont.Equals(sm, naive_m));
    EXPECT_EQ(lx, sm.x);
  }
  // [l] * base_m == O.
  EXPECT_TRUE(mont.IsIdentity(mont.ScalarMul(base_m, l)));

  // --- Ed25519: -x^2 + y^2 = 1 + d*x^2*y^2, d = -121665 / 121666 ---
  const auto a_ed = fp.Neg(BigUInt<256>(1ULL));
  const auto d_ed =
      fp.Neg(fp.Div(BigUInt<256>(121665ULL), BigUInt<256>(121666ULL)));
  EllipticCurve<F256, CurveForm::kTwistedEdwards> ed(fp, a_ed, d_ed);

  const AffinePoint<F256> base_b(
      BigUInt<256>::FromHex(
          "0x216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51a"),
      BigUInt<256>::FromHex("0x666666666666666666666666666666666666666666666666"
                            "6666666666666658"));
  EXPECT_TRUE(ed.IsOnCurve(base_b));
  EXPECT_TRUE(ed.IsIdentity(ed.Identity()));
  EXPECT_TRUE(ed.Equals(ed.Add(base_b, ed.Neg(base_b)), ed.Identity()));

  // Verify w-NAF == constant-time ladder.
  const BigInt k_ed = BigInt::FromHex("0xabcdef0123456789deadbeef");
  const auto ed_wnaf = ed.ScalarMul(base_b, k_ed);
  const auto ed_ct = ed.ScalarMulConstantTime(base_b, k_ed);
  EXPECT_TRUE(ed.IsOnCurve(ed_wnaf));
  EXPECT_TRUE(ed.Equals(ed_wnaf, ed_ct));

  // Verify [l] * B == O!
  EXPECT_TRUE(ed.IsIdentity(ed.ScalarMul(base_b, l)));
  EXPECT_TRUE(ed.Equals(ed.ScalarMul(base_b, l - BigInt(1)), ed.Neg(base_b)));
}

// ============================================================================
// Point Counting (Direct & Shanks-Mestre BSGS)
// ============================================================================

TEST(EllipticCurveTest, CurveOrderDirectAndBsgs) {
  auto verify_curve_order = [](uint32_t p, uint32_t a, uint32_t b) {
    GaloisFieldPrime<uint32_t> gf(p);
    EllipticCurve<GaloisFieldPrime<uint32_t>> curve(gf, a, b);
    const BigInt order = curve.CurveOrder();

    // Hasse bound: |#E - (p + 1)| <= 2 * sqrt(p).
    const BigInt p_big = BigInt::FromUint64(p);
    const BigInt trace = p_big + BigInt(1) - order;
    EXPECT_LE(trace * trace, p_big * BigInt(4));

    // Find a point on E(F_p) and verify [#E]P == O.
    for (uint32_t x = 0; x < p; ++x) {
      const uint32_t rhs =
          gf.Add(gf.Add(gf.Mul(gf.Mul(x, x), x), gf.Mul(a, x)), b);
      const auto y_opt = xg::ec::detail::FieldSqrt(gf, rhs);
      if (y_opt.has_value()) {
        const AffinePoint<GaloisFieldPrime<uint32_t>> pt(x, *y_opt);
        EXPECT_TRUE(curve.IsIdentity(curve.ScalarMul(pt, order)));
        const BigInt pt_ord = curve.PointOrder(pt, order);
        EXPECT_TRUE((order % pt_ord).IsZero());
        break;
      }
    }
    return order;
  };

  // Direct path (p <= 10000): p = 97, p = 1009.
  const BigInt ord97 = verify_curve_order(97U, 2U, 3U);
  EXPECT_EQ(ord97, BigInt(100));

  const BigInt ord1009 = verify_curve_order(1009U, 1U, 1U);
  EXPECT_GT(ord1009, BigInt(900));

  // BSGS Shanks-Mestre path (p > 10000): p = 10007.
  const BigInt ord10007 = verify_curve_order(10007U, 1U, 1U);
  // Cross-check BSGS result against brute-force Legendre count.
  GaloisFieldPrime<uint32_t> gf10007(10007U);
  uint64_t expected = 1;
  for (uint32_t x = 0; x < 10007U; ++x) {
    const uint32_t rhs =
        gf10007.Add(gf10007.Add(gf10007.Mul(gf10007.Mul(x, x), x), x), 1U);
    if (rhs == 0U) {
      expected += 1;
    } else if (xg::ec::detail::FieldHasSqrt(gf10007, rhs)) {
      expected += 2;
    }
  }
  EXPECT_EQ(ord10007, BigInt::FromUint64(expected));
}

// ============================================================================
// Weil & Tate Pairings on Supersingular Curve y^2 = x^3 + x over F_{59^2}
// ============================================================================

TEST(PairingTest, WeilAndTateBilinearityNonDegeneracyAlternating) {
  using BaseF = GaloisFieldPrime<uint32_t>;
  using ExtF = FieldTower<BaseF, 2>;

  // p = 59 == 3 (mod 4), so u^2 + 1 is irreducible over GF(59).
  BaseF fp(59U);
  ExtF fp2 = ExtF::Binomial(fp, 2, fp.Neg(1U), "i");
  EXPECT_EQ(fp2.BigOrder(), BigInt(59 * 59));

  // Supersingular curve E: y^2 = x^3 + x over F_{59^2}.
  // #E(F_59) = 60, so #E(F_{59^2}) = (59 + 1)^2 = 3600.
  EllipticCurve<ExtF> curve(fp2, fp2.MultiplicativeIdentity(),
                            fp2.AdditiveIdentity());
  const BigInt r(5);
  const BigInt group_order(3600);
  EXPECT_EQ(curve.EmbeddingDegree(r), 1ULL);

  const auto [p, q] = curve.TorsionBasis(r, group_order);
  EXPECT_TRUE(curve.IsOnCurve(p));
  EXPECT_TRUE(curve.IsOnCurve(q));
  EXPECT_TRUE(curve.IsIdentity(curve.ScalarMul(p, r)));
  EXPECT_TRUE(curve.IsIdentity(curve.ScalarMul(q, r)));

  // --- Weil Pairing ---
  const auto w_pq = WeilPairing(curve, p, q, r);
  EXPECT_FALSE(fp2.IsOne(w_pq));                        // Non-degenerate
  EXPECT_TRUE(fp2.IsOne(fp2.Pow(w_pq, 5ULL)));          // r-th root of unity
  EXPECT_TRUE(fp2.IsOne(WeilPairing(curve, p, p, r)));  // Alternating
  EXPECT_TRUE(fp2.IsOne(
      WeilPairing(curve, p, curve.Double(p), r)));  // Shifted divisor fallback

  // Skew-symmetry: e(Q, P) == e(P, Q)^(-1).
  const auto w_qp = WeilPairing(curve, q, p, r);
  EXPECT_TRUE(fp2.Equals(fp2.Mul(w_pq, w_qp), fp2.MultiplicativeIdentity()));

  // Bilinearity: e(aP, bQ) == e(P, Q)^(a*b).
  for (uint64_t a = 1; a <= 4; ++a) {
    for (uint64_t b = 1; b <= 4; ++b) {
      const auto ap = curve.ScalarMul(p, BigInt::FromUint64(a));
      const auto bq = curve.ScalarMul(q, BigInt::FromUint64(b));
      const auto lhs = WeilPairing(curve, ap, bq, r);
      const auto rhs = fp2.Pow(w_pq, a * b);
      EXPECT_TRUE(fp2.Equals(lhs, rhs));
    }
  }

  // --- Reduced Tate Pairing ---
  const auto t_pq = TatePairing(curve, p, q, r, fp2.BigOrder());
  EXPECT_FALSE(fp2.IsOne(t_pq));                // Non-degenerate
  EXPECT_TRUE(fp2.IsOne(fp2.Pow(t_pq, 5ULL)));  // r-th root of unity

  for (uint64_t a = 1; a <= 4; ++a) {
    for (uint64_t b = 1; b <= 4; ++b) {
      const auto ap = curve.ScalarMul(p, BigInt::FromUint64(a));
      const auto bq = curve.ScalarMul(q, BigInt::FromUint64(b));
      const auto lhs = TatePairing(curve, ap, bq, r, fp2.BigOrder());
      const auto rhs = fp2.Pow(t_pq, a * b);
      EXPECT_TRUE(fp2.Equals(lhs, rhs));
    }
  }
}

// ============================================================================
// Elliptic Riemann-Roch Matrix for AG Codes on E(F_13)
// ============================================================================

TEST(PairingTest, EllipticRiemannRochMatrixFullRank) {
  GaloisFieldPrime<uint32_t> gf13(13U);
  // y^2 = x^3 + x + 1 over GF(13)
  EllipticCurve<GaloisFieldPrime<uint32_t>> curve(gf13, 1U, 1U);

  std::vector<AffinePoint<GaloisFieldPrime<uint32_t>>> eval_points;
  for (uint32_t x = 0; x < 13U; ++x) {
    const uint32_t rhs = gf13.Add(gf13.Add(gf13.Mul(gf13.Mul(x, x), x), x), 1U);
    for (uint32_t y = 0; y < 13U; ++y) {
      if (gf13.Mul(y, y) == rhs) {
        eval_points.emplace_back(x, y);
      }
    }
  }
  const std::size_t n = eval_points.size();
  ASSERT_GE(n, 8U);

  // Verify rank(EllipticRiemannRochMatrix(E, pts, m)) == m for 1 <= m < n.
  for (uint32_t m = 1; m < static_cast<uint32_t>(n); ++m) {
    const auto mat = EllipticRiemannRochMatrix(curve, eval_points, m);
    EXPECT_EQ(mat.Shape(0), m);
    EXPECT_EQ(mat.Shape(1), n);
    EXPECT_EQ(xg::linalg::MatrixRank(mat), static_cast<std::size_t>(m));
  }
}

// ============================================================================
// MatrixGroup: GL_2(F_3), SL_2(F_3), GL_3(F_2)
// ============================================================================

TEST(MatrixGroupTest, GeneralAndSpecialLinearGroups) {
  auto gf3 = std::make_shared<GaloisFieldPrime<uint32_t>>(3U);
  auto gl2_f3 = std::make_shared<MatrixGroup<GaloisFieldPrime<uint32_t>>>(
      gf3, 2, MatrixGroupKind::kGeneralLinear);
  auto sl2_f3 = std::make_shared<MatrixGroup<GaloisFieldPrime<uint32_t>>>(
      gf3, 2, MatrixGroupKind::kSpecialLinear);

  EXPECT_EQ(gl2_f3->Order(), BigInt(48));
  EXPECT_EQ(sl2_f3->Order(), BigInt(24));

  // Test Contains, Determinant, Inv, Op, GroupElement.
  const std::vector<uint32_t> m1 = {1U, 1U, 0U, 1U};        // det = 1
  const std::vector<uint32_t> m2 = {2U, 1U, 0U, 1U};        // det = 2
  const std::vector<uint32_t> singular = {1U, 2U, 1U, 2U};  // det = 0

  EXPECT_EQ(gl2_f3->Determinant(m1), 1U);
  EXPECT_EQ(gl2_f3->Determinant(m2), 2U);
  EXPECT_EQ(gl2_f3->Determinant(singular), 0U);

  EXPECT_TRUE(gl2_f3->Contains(m1));
  EXPECT_TRUE(gl2_f3->Contains(m2));
  EXPECT_FALSE(gl2_f3->Contains(singular));

  EXPECT_TRUE(sl2_f3->Contains(m1));
  EXPECT_FALSE(sl2_f3->Contains(m2));
  EXPECT_FALSE(sl2_f3->Contains(singular));

  const auto m2_inv = gl2_f3->Inv(m2);
  EXPECT_TRUE(gl2_f3->IsIdentity(gl2_f3->Op(m2, m2_inv)));
  EXPECT_TRUE(gl2_f3->IsIdentity(gl2_f3->Op(m2_inv, m2)));

  GroupElement<MatrixGroup<GaloisFieldPrime<uint32_t>>> ge1(m1, gl2_f3);
  GroupElement<MatrixGroup<GaloisFieldPrime<uint32_t>>> ge2(m2, gl2_f3);
  EXPECT_TRUE((ge1 * ge1.Inv()).IsIdentity());
  EXPECT_EQ(
      (ge2 * ge2.Inv()),
      GroupElement<MatrixGroup<GaloisFieldPrime<uint32_t>>>::Identity(gl2_f3));

  // Verify Generators() generate the entire groups via Subgroup closure.
  const xg::algebra::Subgroup<MatrixGroup<GaloisFieldPrime<uint32_t>>> gl2_sub(
      gl2_f3, gl2_f3->Generators());
  EXPECT_EQ(gl2_sub.Order(), 48ULL);

  const xg::algebra::Subgroup<MatrixGroup<GaloisFieldPrime<uint32_t>>> sl2_sub(
      sl2_f3, sl2_f3->Generators());
  EXPECT_EQ(sl2_sub.Order(), 24ULL);

  // GL_3(F_2) == SL_3(F_2), simple group of order 168.
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gl3_f2 = std::make_shared<MatrixGroup<GaloisFieldBinary>>(
      gf2, 3, MatrixGroupKind::kGeneralLinear);
  EXPECT_EQ(gl3_f2->Order(), BigInt(168));
  const xg::algebra::Subgroup<MatrixGroup<GaloisFieldBinary>> gl3_sub(
      gl3_f2, gl3_f2->Generators());
  EXPECT_EQ(gl3_sub.Order(), 168ULL);
}

}  // namespace
