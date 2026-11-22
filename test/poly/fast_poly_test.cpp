#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "xgalois/coding/decoder/sugiyama_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/additive_fft.hpp"
#include "xgalois/poly/half_gcd.hpp"
#include "xgalois/poly/ntt.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/poly/poly_sparse.hpp"
#include "xgalois/poly/subproduct_tree.hpp"
#include "xgalois/ring/modular_ring.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace {

using namespace xg;

template <typename Field>
PolynomialDense<Field> RandomDensePoly(const std::shared_ptr<Field>& field,
                                       std::size_t degree,
                                       std::mt19937_64& rng) {
  using Raw = typename Field::element_type;
  std::vector<Raw> coeffs(degree + 1);
  for (std::size_t i = 0; i <= degree; ++i) {
    coeffs[i] = field->Random(rng);
  }
  while (field->IsZero(coeffs[degree])) {
    coeffs[degree] = field->Random(rng);
  }
  return PolynomialDense<Field>::FromRaw(std::move(coeffs), field);
}

std::vector<uint64_t> NaiveDft(std::span<const uint64_t> input, uint64_t omega,
                               uint64_t prime) {
  const std::size_t n = input.size();
  std::vector<uint64_t> out(n, 0ULL);
  for (std::size_t k = 0; k < n; ++k) {
    const uint64_t w_k = utils::ModPow(omega, k, prime);
    uint64_t acc = 0ULL;
    uint64_t w_jk = 1ULL;
    for (std::size_t j = 0; j < n; ++j) {
      acc = (acc + utils::ModMul(input[j] % prime, w_jk, prime)) % prime;
      w_jk = utils::ModMul(w_jk, w_k, prime);
    }
    out[k] = acc;
  }
  return out;
}

template <typename Field>
std::tuple<PolynomialDense<Field>, PolynomialDense<Field>,
           PolynomialDense<Field>>
NaivePartialExtendedGcd(const PolynomialDense<Field>& a,
                        const PolynomialDense<Field>& b, int64_t stop_deg) {
  using Poly = PolynomialDense<Field>;
  const auto field = a.Field();
  Poly r0 = a;
  Poly r1 = b;
  Poly u0 = Poly::One(field, a.GetVariable());
  Poly u1 = Poly::Zero(field, a.GetVariable());
  Poly v0 = Poly::Zero(field, a.GetVariable());
  Poly v1 = Poly::One(field, a.GetVariable());

  while (!r1.IsZero() && static_cast<int64_t>(r1.Degree()) >= stop_deg) {
    auto [q, r_next] = r0.DivRem(r1);
    Poly u_next = u0 - q * u1;
    Poly v_next = v0 - q * v1;
    r0 = std::move(r1);
    r1 = std::move(r_next);
    u0 = std::move(u1);
    u1 = std::move(u_next);
    v0 = std::move(v1);
    v1 = std::move(v_next);
  }
  return {std::move(r1), std::move(u1), std::move(v1)};
}

// ============================================================================
// 1. NTT, CRT-NTT, and Bluestein DFT Tests
// ============================================================================

TEST(FastPolyNttTest, NttPlanForwardInverseRoundTripAndMultiply) {
  std::mt19937_64 rng(12345);
  for (uint64_t prime : {998244353ULL, 65537ULL}) {
    NttPlan plan(prime);
    EXPECT_EQ(plan.Prime(), prime);

    // Forward/Inverse round-trip for power-of-2 sizes up to 4096.
    for (std::size_t n :
         {std::size_t{1}, std::size_t{2}, std::size_t{4}, std::size_t{16},
          std::size_t{256}, std::size_t{2048}, std::size_t{4096}}) {
      std::vector<uint64_t> original(n);
      for (std::size_t i = 0; i < n; ++i) {
        original[i] = rng() % prime;
      }
      std::vector<uint64_t> buf = original;
      plan.Forward(std::span<uint64_t>(buf));
      plan.Inverse(std::span<uint64_t>(buf));
      EXPECT_EQ(buf, original) << "prime=" << prime << " n=" << n;
    }

    // Differential Multiply against PolynomialDense up to degree 2048.
    auto field = std::make_shared<GaloisFieldPrime<uint64_t>>(prime);
    for (std::size_t deg :
         {std::size_t{0}, std::size_t{1}, std::size_t{17}, std::size_t{128},
          std::size_t{512}, std::size_t{2048}}) {
      auto a = RandomDensePoly(field, deg, rng);
      auto b = RandomDensePoly(field, deg, rng);
      auto expected = a * b;
      auto actual = plan.Multiply(a, b);
      EXPECT_EQ(actual, expected) << "prime=" << prime << " deg=" << deg;
    }

    // Zero polynomial handling.
    auto zero = PolynomialDense<GaloisFieldPrime<uint64_t>>::Zero(field);
    auto p = RandomDensePoly(field, 15, rng);
    EXPECT_TRUE(plan.Multiply(zero, p).IsZero());
    EXPECT_TRUE(plan.Multiply(p, zero).IsZero());
  }
}

TEST(FastPolyNttTest, NttPlanValidation) {
  EXPECT_THROW(NttPlan(15), std::invalid_argument);
  EXPECT_THROW(NttPlan(2), std::invalid_argument);
  NttPlan small_plan(17);  // 17 - 1 = 16 -> max size 16
  EXPECT_EQ(small_plan.MaxTwoAdicity(), 4U);
  std::vector<uint64_t> bad_len(6, 1ULL);
  EXPECT_THROW(small_plan.Forward(bad_len), std::invalid_argument);
  std::vector<uint64_t> too_big(32, 1ULL);
  EXPECT_THROW(small_plan.Forward(too_big), std::invalid_argument);
}

TEST(FastPolyNttTest, CrtNttPlanArbitraryPrimes) {
  std::mt19937_64 rng(54321);
  for (uint64_t prime : {1000000007ULL, 100003ULL}) {
    CrtNttPlan crt(prime);
    EXPECT_EQ(crt.Modulus(), prime);
    auto field = std::make_shared<GaloisFieldPrime<uint64_t>>(prime);

    for (std::size_t deg : {std::size_t{0}, std::size_t{1}, std::size_t{31},
                            std::size_t{255}, std::size_t{1024}}) {
      auto a = RandomDensePoly(field, deg, rng);
      auto b = RandomDensePoly(field, deg, rng);
      auto expected = a * b;
      auto actual = crt.Multiply(a, b);
      EXPECT_EQ(actual, expected) << "prime=" << prime << " deg=" << deg;
    }

    // Stress maximum coefficient magnitude (all coefficients == prime - 1).
    std::vector<uint64_t> max_a(256, prime - 1);
    std::vector<uint64_t> max_b(256, prime - 1);
    auto max_prod = crt.Multiply(max_a, max_b);
    ASSERT_EQ(max_prod.size(), 511U);
    // (p - 1)^2 == 1 mod p, so coeff[k] == min(k + 1, 511 - k) mod p.
    for (std::size_t k = 0; k < 511; ++k) {
      const uint64_t expected_c = std::min(k + 1, 511 - k) % prime;
      EXPECT_EQ(max_prod[k], expected_c) << "prime=" << prime << " k=" << k;
    }
  }
}

TEST(FastPolyNttTest, BluesteinDftArbitraryLengths) {
  std::mt19937_64 rng(777);
  struct Case {
    uint64_t prime;
    std::size_t n;
  };
  // Includes cases where 2*N does NOT divide p - 1 (e.g. N=10 mod 31, N=12 mod
  // 37) as well as N=15, 30, 12 over GF(31) and GF(61).
  const std::vector<Case> cases = {
      {31ULL, 15U}, {31ULL, 30U}, {31ULL, 10U}, {61ULL, 12U},
      {61ULL, 15U}, {61ULL, 30U}, {37ULL, 12U}, {998244353ULL, 17U},
  };

  for (const auto& tc : cases) {
    BluesteinDft dft(tc.prime, tc.n);
    std::vector<uint64_t> x(tc.n);
    for (std::size_t i = 0; i < tc.n; ++i) {
      x[i] = rng() % tc.prime;
    }

    const auto expected = NaiveDft(x, dft.Root(), tc.prime);
    const auto actual = dft.Forward(x);
    EXPECT_EQ(actual, expected) << "prime=" << tc.prime << " N=" << tc.n;

    const auto recovered = dft.Inverse(actual);
    EXPECT_EQ(recovered, x) << "prime=" << tc.prime << " N=" << tc.n;
  }

  // Invalid length not dividing p - 1 must throw.
  EXPECT_THROW(BluesteinDft(31ULL, 7U), std::invalid_argument);
}

// ============================================================================
// 2. Additive FFT (Lin-Chung-Han 2014) Tests
// ============================================================================

TEST(FastPolyAdditiveFftTest, Gf256EvaluateInterpolateMultiplyAndShift) {
  std::mt19937_64 rng(2025);
  auto gf256 = std::make_shared<GaloisFieldBinaryExtension<uint8_t>>(8);
  AdditiveFftPlan<GaloisFieldBinaryExtension<uint8_t>> plan(gf256);

  // Evaluate vs Horner and Interpolate round-trip across subspace sizes.
  for (std::size_t k : {std::size_t{0}, std::size_t{1}, std::size_t{3},
                        std::size_t{6}, std::size_t{8}}) {
    const std::size_t n = std::size_t{1} << k;
    auto poly = RandomDensePoly(gf256, n - 1, rng);
    const auto evals = plan.Evaluate(poly, k);
    ASSERT_EQ(evals.size(), n);

    const auto pts = plan.EvaluationPoints(k);
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(evals[j], poly.EvalRaw(pts[j])) << "k=" << k << " j=" << j;
    }

    const auto interp = plan.Interpolate(evals, k);
    EXPECT_EQ(interp, poly) << "k=" << k;
  }

  // Affine-shifted ForwardNpb / InverseNpb round-trip + evaluation accuracy.
  {
    const std::size_t k = 6;
    const std::size_t n = std::size_t{1} << k;
    const uint8_t shift = gf256->Random(rng);
    auto poly = RandomDensePoly(gf256, n - 1, rng);

    std::vector<uint8_t> npb = plan.MonomialToNpb(poly.RawCoeffs());
    const std::vector<uint8_t> npb_orig = npb;
    plan.ForwardNpb(std::span<uint8_t>(npb), k, shift);

    const auto shifted_pts = plan.EvaluationPoints(k, shift);
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(npb[j], poly.EvalRaw(shifted_pts[j])) << "shifted j=" << j;
    }

    plan.InverseNpb(std::span<uint8_t>(npb), k, shift);
    EXPECT_EQ(npb, npb_orig);
    auto back_mono = plan.NpbToMonomial(npb);
    EXPECT_EQ(PolynomialDense<GaloisFieldBinaryExtension<uint8_t>>::FromRaw(
                  back_mono, gf256),
              poly);
  }

  // Polynomial multiplication in GF(2^8) up to deg(a)+deg(b) < 256.
  for (std::size_t deg :
       {std::size_t{0}, std::size_t{15}, std::size_t{63}, std::size_t{127}}) {
    auto a = RandomDensePoly(gf256, deg, rng);
    auto b = RandomDensePoly(gf256, deg, rng);
    EXPECT_EQ(plan.Multiply(a, b), a * b) << "GF(2^8) deg=" << deg;
  }
}

TEST(FastPolyAdditiveFftTest, Gf256ZechRepresentationCompatibility) {
  std::mt19937_64 rng(9999);
  auto gf256_zech = std::make_shared<GF2XZECH<uint8_t>>(8);
  AdditiveFftPlan<GF2XZECH<uint8_t>> plan(gf256_zech);

  auto a = RandomDensePoly(gf256_zech, 50, rng);
  auto b = RandomDensePoly(gf256_zech, 50, rng);
  EXPECT_EQ(plan.Multiply(a, b), a * b);

  const auto evals = plan.Evaluate(a, 6);
  const auto pts = plan.EvaluationPoints(6);
  for (std::size_t j = 0; j < evals.size(); ++j) {
    EXPECT_EQ(evals[j], a.EvalRaw(pts[j]));
  }
  EXPECT_EQ(plan.Interpolate(evals, 6), a);
}

TEST(FastPolyAdditiveFftTest, Gf65536EvaluateInterpolateAndMultiplyUpTo1024) {
  std::mt19937_64 rng(424242);
  auto gf65536 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(16);
  AdditiveFftPlan<GaloisFieldBinaryExtension<uint16_t>> plan(gf65536);

  // Evaluate / Interpolate on k = 10 (1024 points).
  {
    const std::size_t k = 10;
    const std::size_t n = std::size_t{1} << k;
    auto poly = RandomDensePoly(gf65536, n - 1, rng);
    const auto evals = plan.Evaluate(poly, k);
    const auto pts = plan.EvaluationPoints(k);
    for (std::size_t j = 0; j < n; j += 17) {
      EXPECT_EQ(evals[j], poly.EvalRaw(pts[j]));
    }
    EXPECT_EQ(plan.Interpolate(evals, k), poly);
  }

  // Multiply up to degree 1024 vs PolynomialDense.
  for (std::size_t deg : {std::size_t{16}, std::size_t{128}, std::size_t{512},
                          std::size_t{1024}}) {
    auto a = RandomDensePoly(gf65536, deg, rng);
    auto b = RandomDensePoly(gf65536, deg, rng);
    EXPECT_EQ(plan.Multiply(a, b), a * b) << "GF(2^16) deg=" << deg;
  }
}

// ============================================================================
// 3. Half-GCD, FastExtendedGcd, PartialExtendedGcd, and Sugiyama RS(1023, 511)
// ============================================================================

template <typename Field>
void RunHalfGcdSuite(const std::shared_ptr<Field>& field, uint64_t seed) {
  std::mt19937_64 rng(seed);
  using Poly = PolynomialDense<Field>;

  // Edge cases.
  {
    auto z = Poly::Zero(field);
    auto p = RandomDensePoly(field, 40, rng);
    auto [g0, u0, v0] = FastExtendedGcd(z, z);
    EXPECT_TRUE(g0.IsZero());
    EXPECT_TRUE(u0.IsZero());
    EXPECT_TRUE(v0.IsZero());

    auto [g1, u1, v1] = FastExtendedGcd(p, z);
    EXPECT_EQ(g1, p.Monic());
    EXPECT_EQ(u1 * p + v1 * z, g1);

    auto [g2, u2, v2] = FastExtendedGcd(z, p);
    EXPECT_EQ(g2, p.Monic());
    EXPECT_EQ(u2 * z + v2 * p, g2);
  }

  // Differential tests across degrees 10..600 (coprime and non-trivial GCD).
  for (std::size_t deg :
       {std::size_t{10}, std::size_t{32}, std::size_t{33}, std::size_t{75},
        std::size_t{150}, std::size_t{320}, std::size_t{600}}) {
    auto a = RandomDensePoly(field, deg, rng);
    auto b = RandomDensePoly(field, deg - 1, rng);

    // Verify HalfGcd directly against NaivePartialExtendedGcd at ceil(deg/2).
    const int64_t mid_deg = (static_cast<int64_t>(deg) + 1) / 2;
    auto hmat = HalfGcd(a, b);
    auto [hr0, hr1] = hmat.Apply(a, b);
    EXPECT_GE(static_cast<int64_t>(hr0.Degree()), mid_deg);
    EXPECT_LT(static_cast<int64_t>(hr1.Degree()), mid_deg);
    auto [nr1, nu1, nv1] = NaivePartialExtendedGcd(a, b, mid_deg);
    EXPECT_EQ(hr1, nr1);
    EXPECT_EQ(hmat.m10, nu1);
    EXPECT_EQ(hmat.m11, nv1);

    // Verify FastExtendedGcd against PolynomialDenseExtendedGcd (coprime).
    auto [fg, fu, fv] = FastExtendedGcd(a, b);
    auto [eg, euv] = utils::PolynomialDenseExtendedGcd(a, b);
    EXPECT_TRUE(fg.IsMonic());
    EXPECT_EQ(fg, eg) << "deg=" << deg;
    EXPECT_EQ(fu * a + fv * b, fg) << "deg=" << deg;
    EXPECT_EQ(fu, euv.first);
    EXPECT_EQ(fv, euv.second);

    // Verify FastExtendedGcd with non-trivial common factor.
    auto common =
        RandomDensePoly(field, std::max<std::size_t>(3, deg / 5), rng);
    auto ca = a * common;
    auto cb = b * common;
    auto [cg, cu, cv] = FastExtendedGcd(ca, cb);
    auto [ceg, ceuv] = utils::PolynomialDenseExtendedGcd(ca, cb);
    EXPECT_EQ(cg, ceg);
    EXPECT_EQ(cu * ca + cv * cb, cg);
    EXPECT_EQ(cu, ceuv.first);
    EXPECT_EQ(cv, ceuv.second);

    // Verify PartialExtendedGcd across several stop degrees.
    for (int64_t stop_deg :
         {int64_t{0}, int64_t{1}, static_cast<int64_t>(deg / 4), mid_deg,
          static_cast<int64_t>((3 * deg) / 4), static_cast<int64_t>(deg)}) {
      auto [pr, pu, pv] = PartialExtendedGcd(a, b, stop_deg);
      auto [er, eu, ev] = NaivePartialExtendedGcd(a, b, stop_deg);
      EXPECT_EQ(pr, er) << "deg=" << deg << " stop_deg=" << stop_deg;
      EXPECT_EQ(pu, eu) << "deg=" << deg << " stop_deg=" << stop_deg;
      EXPECT_EQ(pv, ev) << "deg=" << deg << " stop_deg=" << stop_deg;
      EXPECT_EQ(pu * a + pv * b, pr)
          << "deg=" << deg << " stop_deg=" << stop_deg;
    }
  }
}

TEST(FastPolyHalfGcdTest, PrimeField257AndBinaryField1024) {
  auto gf257 = std::make_shared<GaloisFieldPrime<uint32_t>>(257);
  RunHalfGcdSuite(gf257, 10101);

  auto gf1024 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(10);
  RunHalfGcdSuite(gf1024, 20202);
}

TEST(FastPolyHalfGcdTest, SugiyamaDecoderLargeRs1023x511With256Errors) {
  using Field = GaloisFieldBinaryExtension<uint16_t>;
  auto gf1024 = std::make_shared<Field>(10);
  auto rs =
      coding::GeneralizedReedSolomonCode<Field>::ReedSolomon(gf1024, 1023, 511);
  ASSERT_EQ(rs.Length(), 1023U);
  ASSERT_EQ(rs.Dimension(), 511U);
  ASSERT_EQ(rs.Redundancy(), 512U);
  ASSERT_EQ(rs.MinimumDistance(), 513U);

  coding::SugiyamaDecoder<Field> decoder(rs);
  ASSERT_EQ(decoder.MaxErrors(), 256U);
  ASSERT_EQ(decoder.HalfGcdThreshold(), coding::kSugiyamaHalfGcdThreshold);

  std::mt19937_64 rng(314159);
  auto msg = GFArray<Field>::Random(gf1024, {511}, rng);
  auto cw = rs.Encode(msg);
  ASSERT_TRUE(rs.Contains(cw));

  // Corrupt 256 distinct random positions (full decoding radius!).
  std::vector<std::size_t> positions(1023);
  std::iota(positions.begin(), positions.end(), 0U);
  std::shuffle(positions.begin(), positions.end(), rng);
  positions.resize(256);

  auto received = cw;
  for (std::size_t pos : positions) {
    uint16_t err = gf1024->Random(rng);
    while (gf1024->IsZero(err)) {
      err = gf1024->Random(rng);
    }
    received.Raw(pos) = gf1024->Add(received.Raw(pos), err);
  }

  auto result = decoder.DecodeWithResult(received);
  ASSERT_TRUE(result.success);
  EXPECT_EQ(result.num_errors, 256U);
  EXPECT_EQ(result.codeword, cw);
  EXPECT_EQ(result.message, msg);

  // Also test combined errors + erasures: 200 errors + 112 erasures (2e+s=512).
  std::vector<std::size_t> all_pos(1023);
  std::iota(all_pos.begin(), all_pos.end(), 0U);
  std::shuffle(all_pos.begin(), all_pos.end(), rng);
  std::vector<std::size_t> erasures(all_pos.begin(), all_pos.begin() + 112);
  std::vector<std::size_t> errors(all_pos.begin() + 112, all_pos.begin() + 312);

  auto received_ee = cw;
  for (std::size_t pos : erasures) {
    received_ee.Raw(pos) = gf1024->Random(rng);
  }
  for (std::size_t pos : errors) {
    uint16_t err = gf1024->Random(rng);
    while (gf1024->IsZero(err)) {
      err = gf1024->Random(rng);
    }
    received_ee.Raw(pos) = gf1024->Add(received_ee.Raw(pos), err);
  }

  auto res_ee = decoder.DecodeWithErasures(received_ee, erasures);
  ASSERT_TRUE(res_ee.success);
  EXPECT_EQ(res_ee.codeword, cw);
  EXPECT_EQ(res_ee.message, msg);
}

// ============================================================================
// 4. SubproductTree Tests
// ============================================================================

template <typename Field>
void RunSubproductTreeSuite(const std::shared_ptr<Field>& field,
                            uint64_t seed) {
  std::mt19937_64 rng(seed);
  using Raw = typename Field::element_type;

  for (std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{17},
                        std::size_t{128}, std::size_t{256}}) {
    std::vector<Raw> points;
    points.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      points.push_back(field->SetElementValue(static_cast<Raw>(i)));
    }
    std::shuffle(points.begin(), points.end(), rng);

    SubproductTree<Field> tree(field, points);
    EXPECT_EQ(tree.Size(), n);

    // Evaluate polynomial with degree < n and interpolate back.
    auto poly = RandomDensePoly(field, n - 1, rng);
    const auto evals = tree.MultipointEvaluate(poly);
    ASSERT_EQ(evals.size(), n);
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(evals[j], poly.EvalRaw(points[j])) << "n=" << n << " j=" << j;
    }

    const auto recovered = tree.Interpolate(evals);
    EXPECT_EQ(recovered, poly) << "n=" << n;

    // Also evaluate a higher-degree polynomial (deg > n) vs Horner.
    auto high_poly = RandomDensePoly(field, n + 45, rng);
    const auto high_evals = tree.MultipointEvaluate(high_poly);
    for (std::size_t j = 0; j < n; ++j) {
      EXPECT_EQ(high_evals[j], high_poly.EvalRaw(points[j]))
          << "high_deg n=" << n << " j=" << j;
    }
  }

  // Duplicate points must fail gracefully in Interpolate.
  std::vector<Raw> dup_pts = {field->SetElementValue(1),
                              field->SetElementValue(2),
                              field->SetElementValue(1)};
  SubproductTree<Field> dup_tree(field, dup_pts);
  std::vector<Raw> vals = {field->SetElementValue(1), field->SetElementValue(1),
                           field->SetElementValue(1)};
  EXPECT_THROW((void)dup_tree.Interpolate(vals), std::domain_error);
}

TEST(FastPolySubproductTreeTest, PrimeAndBinaryFields) {
  auto gf257 = std::make_shared<GaloisFieldPrime<uint32_t>>(257);
  RunSubproductTreeSuite(gf257, 8080);

  auto gf1024 = std::make_shared<GaloisFieldBinaryExtension<uint16_t>>(10);
  RunSubproductTreeSuite(gf1024, 9090);
}

// ============================================================================
// 5. PolynomialSparse Tests (over GF(257) and ModularRing(6))
// ============================================================================

TEST(FastPolySparseTest, OverGf257DifferentialVsDense) {
  std::mt19937_64 rng(1337);
  auto gf257 = std::make_shared<GaloisFieldPrime<uint32_t>>(257);
  using Sparse = PolynomialSparse<GaloisFieldPrime<uint32_t>>;

  for (int trial = 0; trial < 20; ++trial) {
    auto dense_a = RandomDensePoly(gf257, 60, rng);
    auto dense_b = RandomDensePoly(gf257, 45, rng);
    auto dense_monic_b = dense_b.Monic();

    auto sparse_a = Sparse::FromDense(dense_a);
    auto sparse_b = Sparse::FromDense(dense_b);
    auto sparse_monic_b = Sparse::FromDense(dense_monic_b);

    EXPECT_EQ(sparse_a.ToDense(), dense_a);
    EXPECT_EQ(sparse_a.Add(sparse_b).ToDense(), dense_a + dense_b);
    EXPECT_EQ(sparse_a.Sub(sparse_b).ToDense(), dense_a - dense_b);
    EXPECT_EQ(sparse_a.Neg().ToDense(), -dense_a);

    const uint32_t scalar = gf257->Random(rng);
    const auto elem_scalar =
        GaloisFieldElement<GaloisFieldPrime<uint32_t>>::FromInternal(scalar,
                                                                     gf257);
    EXPECT_EQ(sparse_a.ScalarMul(scalar).ToDense(), dense_a * elem_scalar);

    // Johnson min-heap multiplication vs PolynomialDense.
    EXPECT_EQ(sparse_a.Mul(sparse_b).ToDense(), dense_a * dense_b);

    // Monagan-Pearce sparse monic division vs PolynomialDense::DivRem.
    auto [sq, sr] = MonicDivRem(sparse_a, sparse_monic_b);
    auto [dq, dr] = dense_a.DivRem(dense_monic_b);
    EXPECT_EQ(sq.ToDense(), dq);
    EXPECT_EQ(sr.ToDense(), dr);

    // Eval vs Horner.
    const uint32_t x = gf257->Random(rng);
    EXPECT_EQ(sparse_a.Eval(x), dense_a.EvalRaw(x));
  }

  // O(nnz) ModBinomial on x^100000 + 1 (i.e. x^100000 - (-1)).
  const uint32_t neg_one = gf257->Neg(gf257->MultiplicativeIdentity());
  auto huge = Sparse::FromTerms(
      gf257,
      {{0, 5}, {3, 7}, {100000, 11}, {100003, 13}, {200000, 19}, {300003, 23}});
  auto reduced = huge.ModBinomial(100000ULL, neg_one);
  // x^100000 == -1, x^200000 == 1, x^300000 == -1.
  // deg 0: 5 - 11 + 19 = 13
  // deg 3: 7 - 13 - 23 = -29 == 228 mod 257
  EXPECT_EQ(reduced.NumTerms(), 2U);
  EXPECT_EQ(reduced.Coeff(0), 13U);
  EXPECT_EQ(reduced.Coeff(3), 228U);
  EXPECT_EQ(reduced.Degree(), 3);

  // Also test ModBinomial with non-unit c (x^100000 - 3) and verify identity
  // at a point x where x^100000 == 3 (or directly on smaller n=7 vs DivRem!).
  auto mod_poly = Sparse::FromTerms(gf257, {{0, gf257->Neg(3U)}, {7, 1U}});
  auto p_small =
      Sparse::FromTerms(gf257, {{0, 4}, {2, 9}, {7, 5}, {14, 8}, {23, 17}});
  auto [q_b, r_b] = MonicDivRem(p_small, mod_poly);
  (void)q_b;
  EXPECT_EQ(p_small.ModBinomial(7ULL, 3U), r_b);
}

TEST(FastPolySparseTest, OverModularRing6WithZeroDivisors) {
  using Ring6 = ring::ModularRing<uint64_t>;
  using Sparse = PolynomialSparse<Ring6>;
  auto z6 = std::make_shared<Ring6>(6ULL);

  // Duplicate degree merging and zero cancellation mod 6.
  auto a = Sparse::FromTerms(z6, {{10, 2}, {0, 3}, {10, 1}, {5, 6}});
  EXPECT_EQ(a.NumTerms(), 2U);
  EXPECT_EQ(a.Coeff(0), 3ULL);
  EXPECT_EQ(a.Coeff(5), 0ULL);
  EXPECT_EQ(a.Coeff(10), 3ULL);
  EXPECT_EQ(a.Degree(), 10);
  EXPECT_EQ(a.LeadingCoeff(), 3ULL);

  // ScalarMul with zero-divisor (2 * 3 == 0 mod 6).
  EXPECT_TRUE(a.ScalarMul(2ULL).IsZero());

  // Product with zero-divisor cancellation: (3x^10 + 3) * (2x^5 + 1) = 3x^10
  // + 3.
  auto b = Sparse::FromTerms(z6, {{0, 1}, {5, 2}});
  EXPECT_EQ(a.Mul(b), a);

  // General multiplication + Monagan-Pearce monic division over Z/6Z.
  auto f = Sparse::FromTerms(
      z6, {{0, 5}, {3, 4}, {17, 2}, {50, 3}, {120, 5}, {250, 1}});
  auto g = Sparse::FromTerms(z6, {{0, 2}, {7, 3}, {25, 4}, {40, 1}});
  auto [q, r] = MonicDivRem(f, g);
  EXPECT_LT(r.Degree(), g.Degree());
  EXPECT_EQ(q.Mul(g).Add(r), f);

  // Verify evaluation homomorphism at every element of Z/6Z.
  for (uint64_t x = 0; x < 6; ++x) {
    const uint64_t lhs = f.Eval(x);
    const uint64_t rhs = z6->Add(z6->Mul(q.Eval(x), g.Eval(x)), r.Eval(x));
    EXPECT_EQ(lhs, rhs) << "x=" << x;
  }

  // ModBinomial on x^100000 + 1 over Z/6Z (c = -1 == 5 mod 6).
  auto huge = Sparse::FromTerms(
      z6, {{0, 1}, {42, 3}, {100000, 4}, {100042, 5}, {200000, 2}});
  auto red = huge.ModBinomial(100000ULL, z6->Neg(1ULL));
  // deg 0: 1 - 4 + 2 = -1 == 5 mod 6
  // deg 42: 3 - 5 = -2 == 4 mod 6
  EXPECT_EQ(red.NumTerms(), 2U);
  EXPECT_EQ(red.Coeff(0), 5ULL);
  EXPECT_EQ(red.Coeff(42), 4ULL);
  EXPECT_FALSE(red.ToString().empty());
}

}  // namespace
