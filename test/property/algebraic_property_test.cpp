#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace {

// ============================================================================
// Field Axiom Checker Helpers
// ============================================================================

template <typename FieldType>
bool IsValidElement(const FieldType &field,
                    const typename FieldType::element_type &elem) {
  using Elem = typename FieldType::element_type;
  if constexpr (std::is_integral_v<Elem>) {
    try {
      const Elem ext = field.GetElementValue(elem);
      const uint64_t q = field.Order();
      if (q != 0 && static_cast<uint64_t>(ext) >= q) {
        return false;
      }
      return true;
    } catch (...) {
      return false;
    }
  } else {
    if (elem.Degree() >= static_cast<int>(field.Degree())) {
      return false;
    }
    const uint64_t p = field.Characteristic();
    for (const auto &c : elem.Coeffs()) {
      if (static_cast<uint64_t>(c) >= p) {
        return false;
      }
    }
    return true;
  }
}

template <typename FieldType>
std::vector<typename FieldType::element_type> EnumerateAllElements(
    const FieldType &field) {
  using Elem = typename FieldType::element_type;
  const uint64_t q = field.Order();
  std::vector<Elem> elems;
  elems.reserve(static_cast<size_t>(q));
  if constexpr (std::is_integral_v<Elem>) {
    for (uint64_t i = 0; i < q; ++i) {
      elems.push_back(field.SetElementValue(static_cast<Elem>(i)));
    }
  } else {
    for (uint64_t i = 0; i < q; ++i) {
      elems.push_back(field.FromInt(i));
    }
  }
  return elems;
}

template <typename FieldType>
void VerifyFieldAxiomsExhaustive(const std::shared_ptr<FieldType> &field) {
  static_assert(FiniteField<FieldType>,
                "FieldType must satisfy FiniteField concept");
  ASSERT_NE(field, nullptr);

  using Elem = typename FieldType::element_type;
  const uint64_t q = field->Order();
  const uint64_t p = field->Characteristic();
  ASSERT_GT(q, 1ULL);

  const Elem zero = field->AdditiveIdentity();
  const Elem one = field->MultiplicativeIdentity();
  const Elem gen = field->MultiplicativeGenerator();

  EXPECT_TRUE(field->IsZero(zero));
  EXPECT_TRUE(field->IsOne(one));
  EXPECT_FALSE(field->IsZero(one));
  EXPECT_FALSE(field->IsOne(zero));
  EXPECT_THROW((void)field->Inv(zero), std::domain_error);
  EXPECT_THROW((void)field->Div(one, zero), std::domain_error);
  EXPECT_THROW((void)field->Log(zero), std::domain_error);

  const std::vector<Elem> elems = EnumerateAllElements(*field);
  ASSERT_EQ(elems.size(), static_cast<size_t>(q));

  // Verify generator order.
  EXPECT_EQ(field->MultiplicativeOrder(gen), q - 1);

  std::vector<Elem> nonzero;
  nonzero.reserve(elems.size() - 1);

  // Unary axioms, identities, inverses, Fermat, Sqrt, Log, Trace/Norm subfield.
  for (const Elem &a : elems) {
    ASSERT_TRUE(IsValidElement(*field, a));

    // Identities.
    EXPECT_EQ(field->Add(a, zero), a);
    EXPECT_EQ(field->Add(zero, a), a);
    EXPECT_EQ(field->Mul(a, one), a);
    EXPECT_EQ(field->Mul(one, a), a);
    EXPECT_EQ(field->Mul(a, zero), zero);
    EXPECT_EQ(field->Mul(zero, a), zero);

    // Additive inverse.
    const Elem neg_a = field->Neg(a);
    ASSERT_TRUE(IsValidElement(*field, neg_a));
    EXPECT_EQ(field->Add(a, neg_a), zero);
    EXPECT_EQ(field->Sub(a, a), zero);

    // Fermat's Little Theorem: a^q = a.
    EXPECT_EQ(field->Pow(a, q), a);

    // Frobenius matches a^p.
    EXPECT_EQ(field->Frobenius(a), field->Pow(a, p));

    // Trace & Norm land in prime subfield F_p (fixed by Frobenius).
    const Elem tr_a = field->Trace(a);
    const Elem nm_a = field->Norm(a);
    ASSERT_TRUE(IsValidElement(*field, tr_a));
    ASSERT_TRUE(IsValidElement(*field, nm_a));
    EXPECT_EQ(field->Frobenius(tr_a), tr_a);
    EXPECT_EQ(field->Frobenius(nm_a), nm_a);

    // Sqrt check.
    if (field->IsZero(a)) {
      EXPECT_EQ(field->Sqrt(a), zero);
    } else {
      nonzero.push_back(a);

      // Multiplicative inverse & Fermat a^(q-1) = 1.
      const Elem inv_a = field->Inv(a);
      ASSERT_TRUE(IsValidElement(*field, inv_a));
      EXPECT_EQ(field->Mul(a, inv_a), one);
      EXPECT_EQ(field->Mul(inv_a, a), one);
      EXPECT_EQ(field->Pow(a, q - 1), one);

      // Discrete log round-trip: g^{Log_g(a)} == a.
      const uint64_t lg = field->Log(a);
      EXPECT_LT(lg, q - 1);
      EXPECT_EQ(field->Pow(gen, lg), a);

      // Sqrt(a)^2 == a for all elements in char 2, and quadratic residues in
      // odd char.
      if (p == 2) {
        const Elem sq = field->Sqrt(a);
        ASSERT_TRUE(IsValidElement(*field, sq));
        EXPECT_EQ(field->Mul(sq, sq), a);
      } else {
        const bool is_qr = field->IsOne(field->Pow(a, (q - 1) >> 1));
        if (is_qr) {
          const Elem sq = field->Sqrt(a);
          ASSERT_TRUE(IsValidElement(*field, sq));
          EXPECT_EQ(field->Mul(sq, sq), a);
        } else {
          EXPECT_THROW((void)field->Sqrt(a), std::domain_error);
        }
      }
      // Every square a^2 must have a valid square root.
      const Elem a_sq = field->Mul(a, a);
      const Elem root = field->Sqrt(a_sq);
      EXPECT_EQ(field->Mul(root, root), a_sq);
    }
  }

  // Pairwise axioms: closure, commutativity, Sub, Div, Frobenius, Trace, Norm.
  for (const Elem &a : elems) {
    const Elem frob_a = field->Frobenius(a);
    const Elem tr_a = field->Trace(a);
    const Elem nm_a = field->Norm(a);
    for (const Elem &b : elems) {
      const Elem sum = field->Add(a, b);
      const Elem prod = field->Mul(a, b);
      ASSERT_TRUE(IsValidElement(*field, sum));
      ASSERT_TRUE(IsValidElement(*field, prod));

      EXPECT_EQ(sum, field->Add(b, a));
      EXPECT_EQ(prod, field->Mul(b, a));
      EXPECT_EQ(field->Sub(a, b), field->Add(a, field->Neg(b)));

      if (!field->IsZero(b)) {
        const Elem quot = field->Div(a, b);
        ASSERT_TRUE(IsValidElement(*field, quot));
        EXPECT_EQ(quot, field->Mul(a, field->Inv(b)));
        EXPECT_EQ(field->Mul(quot, b), a);
      }

      // Frobenius endomorphism: (a + b)^p = a^p + b^p.
      EXPECT_EQ(field->Frobenius(sum), field->Add(frob_a, field->Frobenius(b)));

      // Trace additivity: Trace(a + b) = Trace(a) + Trace(b).
      EXPECT_EQ(field->Trace(sum), field->Add(tr_a, field->Trace(b)));

      // Norm multiplicativity: Norm(a * b) = Norm(a) * Norm(b).
      EXPECT_EQ(field->Norm(prod), field->Mul(nm_a, field->Norm(b)));
    }
  }

  // Triple axioms: associativity and distributivity (exhaustive for q <= 16).
  for (const Elem &a : elems) {
    for (const Elem &b : elems) {
      const Elem ab_add = field->Add(a, b);
      const Elem ab_mul = field->Mul(a, b);
      for (const Elem &c : elems) {
        EXPECT_EQ(field->Add(ab_add, c), field->Add(a, field->Add(b, c)));
        EXPECT_EQ(field->Mul(ab_mul, c), field->Mul(a, field->Mul(b, c)));
        EXPECT_EQ(field->Mul(a, field->Add(b, c)),
                  field->Add(ab_mul, field->Mul(a, c)));
      }
    }
  }

  // BatchInv: BatchInv(v)[i] * v[i] == 1.
  std::vector<Elem> invs(nonzero.size());
  field->BatchInv(std::span<Elem>(invs), std::span<const Elem>(nonzero));
  for (size_t i = 0; i < nonzero.size(); ++i) {
    EXPECT_EQ(field->Mul(invs[i], nonzero[i]), one);
  }
}

template <typename FieldType>
void VerifyFieldAxiomsRandomized(const std::shared_ptr<FieldType> &field,
                                 uint64_t seed, size_t num_samples = 48,
                                 size_t log_samples = 10) {
  static_assert(FiniteField<FieldType>,
                "FieldType must satisfy FiniteField concept");
  ASSERT_NE(field, nullptr);

  using Elem = typename FieldType::element_type;
  const uint64_t q = field->Order();  // 0 denotes 2^64.
  const uint64_t group_order = (q == 0) ? ~uint64_t{0} : (q - 1);
  const uint64_t p = field->Characteristic();

  const Elem zero = field->AdditiveIdentity();
  const Elem one = field->MultiplicativeIdentity();
  const Elem gen = field->MultiplicativeGenerator();

  EXPECT_TRUE(field->IsZero(zero));
  EXPECT_TRUE(field->IsOne(one));
  EXPECT_THROW((void)field->Inv(zero), std::domain_error);
  EXPECT_THROW((void)field->Div(one, zero), std::domain_error);
  EXPECT_THROW((void)field->Log(zero), std::domain_error);
  EXPECT_EQ(field->MultiplicativeOrder(gen), group_order);

  std::mt19937_64 rng(seed);
  std::vector<Elem> samples;
  samples.reserve(num_samples + 3);
  samples.push_back(zero);
  samples.push_back(one);
  samples.push_back(gen);
  for (size_t i = 0; i < num_samples; ++i) {
    samples.push_back(field->Random(rng));
  }

  std::vector<Elem> nonzero;
  nonzero.reserve(samples.size());
  size_t logs_checked = 0;

  for (const Elem &a : samples) {
    ASSERT_TRUE(IsValidElement(*field, a));

    EXPECT_EQ(field->Add(a, zero), a);
    EXPECT_EQ(field->Mul(a, one), a);
    EXPECT_EQ(field->Mul(a, zero), zero);
    EXPECT_EQ(field->Add(a, field->Neg(a)), zero);
    EXPECT_EQ(field->Sub(a, a), zero);

    // Fermat's Little Theorem: a^q = a (and a^(q-1) = 1 for a != 0).
    if (q != 0) {
      EXPECT_EQ(field->Pow(a, q), a);
    } else {
      EXPECT_EQ(field->Mul(field->Pow(a, group_order), a), a);
    }

    // Frobenius & Trace/Norm prime-subfield invariance.
    EXPECT_EQ(field->Frobenius(a), field->Pow(a, p));
    const Elem tr_a = field->Trace(a);
    const Elem nm_a = field->Norm(a);
    EXPECT_EQ(field->Frobenius(tr_a), tr_a);
    EXPECT_EQ(field->Frobenius(nm_a), nm_a);

    if (field->IsZero(a)) {
      EXPECT_EQ(field->Sqrt(a), zero);
    } else {
      nonzero.push_back(a);
      const Elem inv_a = field->Inv(a);
      EXPECT_EQ(field->Mul(a, inv_a), one);
      EXPECT_EQ(field->Pow(a, group_order), one);

      // Sqrt axiom.
      if (p == 2) {
        const Elem sq = field->Sqrt(a);
        EXPECT_EQ(field->Mul(sq, sq), a);
      } else {
        const Elem a_sq = field->Mul(a, a);
        const Elem root = field->Sqrt(a_sq);
        EXPECT_EQ(field->Mul(root, root), a_sq);

        const bool is_qr = field->IsOne(field->Pow(a, group_order >> 1));
        if (is_qr) {
          const Elem sq = field->Sqrt(a);
          EXPECT_EQ(field->Mul(sq, sq), a);
        } else {
          EXPECT_THROW((void)field->Sqrt(a), std::domain_error);
        }
      }

      // Discrete log round-trip on bounded sample budget.
      if (logs_checked < log_samples) {
        const uint64_t lg = field->Log(a);
        EXPECT_EQ(field->Pow(gen, lg), a);
        ++logs_checked;
      }
    }
  }

  // Pairwise & triple axioms across randomized tuples.
  for (size_t i = 0; i < samples.size(); ++i) {
    const Elem &a = samples[i];
    const Elem &b = samples[(i * 7 + 3) % samples.size()];
    const Elem &c = samples[(i * 13 + 5) % samples.size()];

    const Elem sum_ab = field->Add(a, b);
    const Elem prod_ab = field->Mul(a, b);
    ASSERT_TRUE(IsValidElement(*field, sum_ab));
    ASSERT_TRUE(IsValidElement(*field, prod_ab));

    EXPECT_EQ(sum_ab, field->Add(b, a));
    EXPECT_EQ(prod_ab, field->Mul(b, a));
    EXPECT_EQ(field->Sub(a, b), field->Add(a, field->Neg(b)));
    if (!field->IsZero(b)) {
      EXPECT_EQ(field->Mul(field->Div(a, b), b), a);
    }

    EXPECT_EQ(field->Add(sum_ab, c), field->Add(a, field->Add(b, c)));
    EXPECT_EQ(field->Mul(prod_ab, c), field->Mul(a, field->Mul(b, c)));
    EXPECT_EQ(field->Mul(a, field->Add(b, c)),
              field->Add(prod_ab, field->Mul(a, c)));

    EXPECT_EQ(field->Frobenius(sum_ab),
              field->Add(field->Frobenius(a), field->Frobenius(b)));
    EXPECT_EQ(field->Trace(sum_ab),
              field->Add(field->Trace(a), field->Trace(b)));
    EXPECT_EQ(field->Norm(prod_ab), field->Mul(field->Norm(a), field->Norm(b)));
  }

  // BatchInv.
  std::vector<Elem> invs(nonzero.size());
  field->BatchInv(std::span<Elem>(invs), std::span<const Elem>(nonzero));
  for (size_t i = 0; i < nonzero.size(); ++i) {
    EXPECT_EQ(field->Mul(invs[i], nonzero[i]), one);
  }
}

// ============================================================================
// Polynomial Axiom Checker
// ============================================================================

template <typename FieldType>
PolynomialDense<FieldType> RandomPoly(const std::shared_ptr<FieldType> &field,
                                      size_t max_deg, std::mt19937_64 &rng,
                                      bool force_nonzero = false) {
  using Raw = typename FieldType::element_type;
  std::uniform_int_distribution<size_t> deg_dist(0, max_deg);
  const size_t deg = deg_dist(rng);
  std::vector<Raw> coeffs(deg + 1);
  for (size_t i = 0; i <= deg; ++i) {
    coeffs[i] = field->Random(rng);
  }
  if (force_nonzero) {
    while (field->IsZero(coeffs[deg])) {
      coeffs[deg] = field->Random(rng);
    }
  }
  return PolynomialDense<FieldType>::FromRaw(std::move(coeffs), field);
}

template <typename FieldType>
void VerifyPolynomialAxioms(const std::shared_ptr<FieldType> &field,
                            uint64_t seed) {
  using Poly = PolynomialDense<FieldType>;
  using Elem = GaloisFieldElement<FieldType>;
  using Raw = typename FieldType::element_type;

  std::mt19937_64 rng(seed);

  // 1. Euclidean division, Bezout identity, and derivative product rule.
  for (int trial = 0; trial < 30; ++trial) {
    const Poly a = RandomPoly(field, 9, rng, false);
    const Poly b = RandomPoly(field, 6, rng, true);

    // Euclidean division: a = q * b + r, deg(r) < deg(b).
    auto [q, r] = a.DivRem(b);
    EXPECT_EQ(q * b + r, a);
    EXPECT_LT(r.Degree(), b.Degree());

    // Bezout identity: s * a + t * b = gcd(a, b).
    auto [g, st] = utils::PolynomialDenseExtendedGcd(a, b);
    const auto &[s, t] = st;
    EXPECT_EQ(s * a + t * b, g);
    EXPECT_EQ(g, utils::PolynomialGcd(a, b));
    if (!g.IsZero()) {
      EXPECT_TRUE(g.IsMonic());
      EXPECT_TRUE((a % g).IsZero());
      EXPECT_TRUE((b % g).IsZero());
    }

    // Formal derivative product rule: (a * b)' = a' * b + a * b'.
    EXPECT_EQ((a * b).Derivative(), a.Derivative() * b + a * b.Derivative());
    EXPECT_EQ((a + b).Derivative(), a.Derivative() + b.Derivative());
  }

  // 2. Factorization round-trip: lc * prod(f_i^{e_i}) == f, each f_i monic &
  //    irreducible.
  for (int trial = 0; trial < 15; ++trial) {
    Poly f = RandomPoly(field, 7, rng, true);
    if (f.Degree() == 0) {
      f = f + Poly::Monomial(field, 1, field->MultiplicativeIdentity());
    }
    const auto [lc, factors] =
        utils::Factor(f, seed + static_cast<uint64_t>(trial));
    Poly reconstructed = Poly::Monomial(field, 0, lc.RawValue());
    for (const auto &[factor, mult] : factors) {
      EXPECT_TRUE(factor.IsMonic());
      EXPECT_TRUE(utils::IsIrreducible(factor));
      EXPECT_GT(mult, 0U);
      reconstructed *= factor.Pow(mult);
    }
    EXPECT_EQ(reconstructed, f);
  }

  // 3. Interpolation: p(x_i) == y_i and deg(p) < k.
  const std::vector<Raw> all_raw = EnumerateAllElements(*field);
  const size_t k = std::min<size_t>(all_raw.size(), 7);
  for (int trial = 0; trial < 10; ++trial) {
    std::vector<Raw> perm = all_raw;
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<Elem> xs;
    std::vector<Elem> ys;
    xs.reserve(k);
    ys.reserve(k);
    for (size_t i = 0; i < k; ++i) {
      xs.push_back(Elem::FromInternal(perm[i], field));
      ys.push_back(Elem::FromInternal(field->Random(rng), field));
    }
    const Poly p = utils::Interpolate(xs, ys);
    EXPECT_LT(p.Degree(), static_cast<int>(k));
    for (size_t i = 0; i < k; ++i) {
      EXPECT_EQ(p(xs[i]), ys[i]);
    }
  }
}

// ============================================================================
// Linear Algebra Axiom Checker
// ============================================================================

template <typename FieldType>
GFArray<FieldType> RandomInvertibleMatrix(
    const std::shared_ptr<FieldType> &field, size_t n, std::mt19937_64 &rng) {
  for (int attempt = 0; attempt < 256; ++attempt) {
    auto mat = GFArray<FieldType>::Random(field, {n, n}, rng);
    if (!linalg::det(mat).IsZero()) {
      return mat;
    }
  }
  return linalg::eye(field, n);
}

template <typename FieldType>
void VerifyLinalgAxioms(const std::shared_ptr<FieldType> &field,
                        uint64_t seed) {
  std::mt19937_64 rng(seed);
  const size_t n = 5;
  const auto I = linalg::eye(field, n);

  for (int trial = 0; trial < 8; ++trial) {
    const auto A = RandomInvertibleMatrix(field, n, rng);
    const auto B = GFArray<FieldType>::Random(field, {n, n}, rng);

    // 1. A * A^{-1} == I and A^{-1} * A == I.
    const auto A_inv = linalg::inv(A);
    EXPECT_EQ(linalg::dot(A, A_inv), I);
    EXPECT_EQ(linalg::dot(A_inv, A), I);

    // 2. det(A * B) == det(A) * det(B).
    const auto AB = linalg::dot(A, B);
    EXPECT_EQ(linalg::det(AB), linalg::det(A) * linalg::det(B));

    // 3. Linear solve: A * x == b.
    const auto b_vec = GFArray<FieldType>::Random(field, {n}, rng);
    const auto x_vec = linalg::solve(A, b_vec);
    EXPECT_EQ(linalg::dot(A, x_vec), b_vec);

    // 4. RREF idempotence: rref(rref(M)) == rref(M).
    const auto rect = GFArray<FieldType>::Random(field, {4, 6}, rng);
    const auto r1 = linalg::rref(rect);
    const auto r2 = linalg::rref(r1);
    EXPECT_EQ(r1, r2);

    // 5. Nullspace orthogonality: M * nullspace(M)^T == 0 and rank-nullity.
    const auto N = linalg::nullspace(rect);
    const size_t rank = linalg::matrix_rank(rect);
    EXPECT_EQ(rank + N.Shape(0), rect.Shape(1));
    if (N.Shape(0) > 0) {
      const auto prod = linalg::dot(rect, linalg::transpose(N));
      EXPECT_EQ(prod, linalg::zeros(field, {rect.Shape(0), N.Shape(0)}));
      EXPECT_EQ(linalg::matrix_rank(N), N.Shape(0));
    }

    // 6. PLU decomposition: P * A == L * U.
    const auto plu_res = linalg::plu(B);
    EXPECT_EQ(linalg::dot(plu_res.P, B), linalg::dot(plu_res.L, plu_res.U));
    // Verify unit lower-triangular L and upper-triangular U.
    for (size_t i = 0; i < n; ++i) {
      EXPECT_TRUE(field->IsOne(plu_res.L.Raw(i, i)));
      for (size_t j = i + 1; j < n; ++j) {
        EXPECT_TRUE(field->IsZero(plu_res.L.Raw(i, j)));
      }
      for (size_t j = 0; j < i; ++j) {
        EXPECT_TRUE(field->IsZero(plu_res.U.Raw(i, j)));
      }
    }

    // 7. Cayley-Hamilton theorem: p_B(B) == 0.
    const auto char_poly = linalg::characteristic_poly(B);
    ASSERT_EQ(char_poly.Degree(), static_cast<int>(n));
    auto eval = linalg::zeros(field, {n, n});
    auto B_pow = I;
    for (size_t d = 0; d <= n; ++d) {
      const auto coeff = char_poly[d];
      if (!coeff.IsZero()) {
        eval = eval + (B_pow * coeff);
      }
      if (d < n) {
        B_pow = linalg::dot(B_pow, B);
      }
    }
    EXPECT_EQ(eval, linalg::zeros(field, {n, n}));
  }
}

// ============================================================================
// Coding Axiom Checker Helper
// ============================================================================

template <typename CodeType>
void VerifyLinearCodeAxioms(const CodeType &code, uint64_t seed,
                            size_t num_messages = 10) {
  using FieldType = typename CodeType::field_type;
  const auto field = code.Field();
  const auto &G = code.GeneratorMatrix();
  const auto &H = code.ParityCheckMatrix();

  const size_t k = code.Dimension();
  const size_t n = code.Length();
  ASSERT_EQ(G.Shape(0), k);
  ASSERT_EQ(G.Shape(1), n);
  ASSERT_EQ(H.Shape(0), n - k);
  ASSERT_EQ(H.Shape(1), n);

  // G * H^T == 0.
  if (k > 0 && n > k) {
    const auto ght = linalg::dot(G, linalg::transpose(H));
    EXPECT_EQ(ght, linalg::zeros(field, {k, n - k}));
  }

  // Encode / Contains / Syndrome / Unencode round-trip.
  std::mt19937_64 rng(seed);
  for (size_t i = 0; i < num_messages; ++i) {
    const auto msg = GFArray<FieldType>::Random(field, {k}, rng);
    const auto cw = code.Encode(msg);
    EXPECT_EQ(cw.Size(), n);
    EXPECT_TRUE(code.Contains(cw));
    EXPECT_EQ(code.Unencode(cw), msg);
  }
}

// ============================================================================
// Field Axiom Tests — Exhaustive
// ============================================================================

TEST(FieldAxiomsExhaustive, GF2) {
  VerifyFieldAxiomsExhaustive(std::make_shared<GF2>());
}

TEST(FieldAxiomsExhaustive, GF7) {
  VerifyFieldAxiomsExhaustive(std::make_shared<GFP<uint8_t>>(7));
}

TEST(FieldAxiomsExhaustive, GFPLOG7) {
  VerifyFieldAxiomsExhaustive(std::make_shared<GFPLOG<uint8_t>>(7));
}

TEST(FieldAxiomsExhaustive, GF2X4) {
  VerifyFieldAxiomsExhaustive(std::make_shared<GF2X<uint8_t>>(4));
}

TEST(FieldAxiomsExhaustive, GF2XLOG4) {
  VerifyFieldAxiomsExhaustive(std::make_shared<GF2XLOG<uint8_t>>(4));
}

TEST(FieldAxiomsExhaustive, GF2XZECH4) {
  VerifyFieldAxiomsExhaustive(std::make_shared<GF2XZECH<uint8_t>>(4));
}

TEST(FieldAxiomsExhaustive, GFPX3_2) {
  VerifyFieldAxiomsExhaustive(
      std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL)));
}

// ============================================================================
// Field Axiom Tests — Randomized Large Fields
// ============================================================================

TEST(FieldAxiomsRandomized, GF65521) {
  VerifyFieldAxiomsRandomized(std::make_shared<GFP<uint16_t>>(65521), 101);
}

TEST(FieldAxiomsRandomized, GF2X16) {
  VerifyFieldAxiomsRandomized(std::make_shared<GF2X<uint16_t>>(16), 202);
}

TEST(FieldAxiomsRandomized, GF2X32) {
  VerifyFieldAxiomsRandomized(std::make_shared<GF2X<uint32_t>>(32), 303);
}

TEST(FieldAxiomsRandomized, GF2X64) {
  VerifyFieldAxiomsRandomized(std::make_shared<GF2X<uint64_t>>(64), 404, 32, 6);
}

TEST(FieldAxiomsRandomized, GFPX5_2) {
  VerifyFieldAxiomsRandomized(
      std::make_shared<GFPX<uint8_t>>(std::make_pair(5ULL, 2ULL)), 505);
}

// ============================================================================
// Polynomial Axiom Tests
// ============================================================================

TEST(PolynomialAxioms, GF2) {
  VerifyPolynomialAxioms(std::make_shared<GF2>(), 1001);
}

TEST(PolynomialAxioms, GF7) {
  VerifyPolynomialAxioms(std::make_shared<GFP<uint8_t>>(7), 1002);
}

TEST(PolynomialAxioms, GF2X4) {
  VerifyPolynomialAxioms(std::make_shared<GF2X<uint8_t>>(4), 1003);
}

TEST(PolynomialAxioms, GFPX3_2) {
  VerifyPolynomialAxioms(
      std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL)), 1004);
}

// ============================================================================
// Linear Algebra Axiom Tests
// ============================================================================

TEST(LinalgAxioms, GF2) { VerifyLinalgAxioms(std::make_shared<GF2>(), 2001); }

TEST(LinalgAxioms, GF11) {
  VerifyLinalgAxioms(std::make_shared<GFP<uint8_t>>(11), 2002);
}

TEST(LinalgAxioms, GF2X8) {
  VerifyLinalgAxioms(std::make_shared<GF2X<uint8_t>>(8), 2003);
}

TEST(LinalgAxioms, GFPX3_2) {
  VerifyLinalgAxioms(
      std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL)), 2004);
}

// ============================================================================
// Coding + Channel End-to-End Property Tests
// ============================================================================

TEST(CodingPropertySuite, GeneratorParityAndEncodeUnencodeAllCodeFamilies) {
  // 1. Custom LinearCode over GF(5).
  {
    auto gf5 = std::make_shared<GFP<uint8_t>>(5);
    GFArray<GFP<uint8_t>>::raw_array_type g_raw = {
        {1, 0, 0, 1, 2, 3},
        {0, 1, 0, 2, 3, 4},
        {0, 0, 1, 3, 4, 1},
    };
    auto G = GFArray<GFP<uint8_t>>(gf5, g_raw);
    auto code = coding::LinearCode<GFP<uint8_t>>::FromGeneratorMatrix(gf5, G);
    VerifyLinearCodeAxioms(code, 3001);
  }

  // 2. Binary and non-binary HammingCode.
  VerifyLinearCodeAxioms(coding::HammingCode(3), 3002);
  {
    auto gf3 = std::make_shared<GFP<uint8_t>>(3);
    VerifyLinearCodeAxioms(coding::HammingCode(gf3, 2), 3003);
  }

  // 3. BinaryGolayCode and ExtendedBinaryGolayCode.
  VerifyLinearCodeAxioms(coding::BinaryGolayCode(), 3004);
  VerifyLinearCodeAxioms(coding::ExtendedBinaryGolayCode(), 3005);

  // 4. TernaryGolayCode.
  VerifyLinearCodeAxioms(coding::TernaryGolayCode(), 3006);

  // 5. ReedMullerCode RM(1, 4).
  VerifyLinearCodeAxioms(coding::ReedMullerCode(1, 4), 3007);

  // 6. CyclicCode [7, 4] over GF(2) with g(x) = 1 + x + x^3.
  {
    auto gf2 = std::make_shared<GF2>();
    auto g_poly = PolynomialDense<GF2>::FromRaw({1, 1, 0, 1}, gf2);
    coding::CyclicCode<GF2> cyclic(gf2, 7, g_poly);
    VerifyLinearCodeAxioms(cyclic, 3008);
  }

  // 7. Binary BCH(15, 7, t=2).
  {
    auto gf2 = std::make_shared<GF2>();
    auto gf16 = std::make_shared<GF2XLOG<uint8_t>>(4);
    auto bch =
        coding::BCHCode<GF2, GF2XLOG<uint8_t>>::FromErrorCorrectionCapability(
            gf2, gf16, 2);
    VerifyLinearCodeAxioms(bch, 3009);
  }

  // 8. GeneralizedReedSolomonCode including CCSDS RS(255, 223).
  {
    auto gf16 = std::make_shared<GF2XLOG<uint8_t>>(4);
    auto rs15 = coding::CreateReedSolomonCode(gf16, 15, 11);
    VerifyLinearCodeAxioms(rs15, 3010);

    auto ext_rs16 = coding::CreateExtendedReedSolomonCode(gf16, 16, 12);
    VerifyLinearCodeAxioms(ext_rs16, 3011);
  }
  {
    auto gf256 = std::make_shared<GF2XLOG<uint8_t>>(8);
    auto ccsds_rs = coding::CreateReedSolomonCode(gf256, 255, 223);
    VerifyLinearCodeAxioms(ccsds_rs, 3012, 4);

    // Also verify CCSDS RS(255, 223) corrects full t = 16 symbol errors.
    std::mt19937_64 rng(3013);
    const auto msg = GFArray<GF2XLOG<uint8_t>>::Random(gf256, {223}, rng);
    const auto cw = ccsds_rs.Encode(msg);
    channel::StaticErrorRateChannel<GF2XLOG<uint8_t>> ch(
        gf256, 255, 16, std::optional<uint64_t>{3014ULL});
    const auto rx = ch.transmit(cw);
    const auto res = ccsds_rs.DecodeWithResult(rx);
    ASSERT_TRUE(res.Ok()) << res.error_message;
    EXPECT_EQ(res.codeword, cw);
    EXPECT_EQ(res.message, msg);
    EXPECT_EQ(res.num_errors, 16U);
  }
}

TEST(CodingChannelSimulation, BCHOverBinarySymmetricChannel) {
  auto gf2 = std::make_shared<GF2>();
  auto gf32 = std::make_shared<GF2XLOG<uint8_t>>(5);
  auto bch =
      coding::BCHCode<GF2, GF2XLOG<uint8_t>>::FromErrorCorrectionCapability(
          gf2, gf32, 3);
  ASSERT_EQ(bch.Length(), 31U);
  ASSERT_EQ(bch.Dimension(), 16U);

  // Zero noise -> 0% FER.
  channel::BinarySymmetricChannel bsc_clean(gf2, bch.Length(), 0.0, 4001);
  const auto clean_stats = channel::SimulateCodePerformance(
      bch, bsc_clean, channel::SimulationConfig{40, 0, 4001});
  EXPECT_EQ(clean_stats.frames_sent, 40U);
  EXPECT_EQ(clean_stats.frame_errors, 0U);
  EXPECT_DOUBLE_EQ(clean_stats.FER(), 0.0);

  // Heavy noise -> positive FER, and deterministic across identical seeds.
  channel::BinarySymmetricChannel bsc_noisy_a(gf2, bch.Length(), 0.22, 4002);
  channel::BinarySymmetricChannel bsc_noisy_b(gf2, bch.Length(), 0.22, 4002);
  const auto noisy_a = channel::SimulateCodePerformance(
      bch, bsc_noisy_a, channel::SimulationConfig{60, 0, 4002});
  const auto noisy_b = channel::SimulateCodePerformance(
      bch, bsc_noisy_b, channel::SimulationConfig{60, 0, 4002});
  EXPECT_GT(noisy_a.FER(), 0.0);
  EXPECT_EQ(noisy_a.frames_sent, noisy_b.frames_sent);
  EXPECT_EQ(noisy_a.frame_errors, noisy_b.frame_errors);
  EXPECT_EQ(noisy_a.bit_errors, noisy_b.bit_errors);
  EXPECT_EQ(noisy_a.decoding_failures, noisy_b.decoding_failures);
}

TEST(CodingChannelSimulation, RS31_23OverErrorErasureChannelCapacityBoundary) {
  auto gf32 = std::make_shared<GF2XLOG<uint8_t>>(5);
  auto rs = coding::CreateReedSolomonCode(gf32, 31, 23);
  ASSERT_EQ(rs.Redundancy(), 8U);  // n - k = 8 => corrects iff 2e + s <= 8.

  // Within capacity: 2e + s <= 8 -> 0% FER across all boundary combinations.
  const std::vector<std::pair<size_t, size_t>> within_capacity = {
      {4, 0}, {3, 2}, {2, 4}, {1, 6}, {0, 8},
  };
  for (const auto &[e, s] : within_capacity) {
    channel::ErrorErasureChannel<GF2XLOG<uint8_t>> ch(gf32, 31, e, s,
                                                      5000 + e * 10 + s);
    const auto stats = channel::SimulateCodePerformance(
        rs, ch, channel::SimulationConfig{30, 0, 5100 + e * 10 + s});
    EXPECT_EQ(stats.frames_sent, 30U);
    EXPECT_EQ(stats.frame_errors, 0U) << "Failed for e=" << e << ", s=" << s;
    EXPECT_DOUBLE_EQ(stats.FER(), 0.0) << "Failed for e=" << e << ", s=" << s;
    EXPECT_EQ(stats.decoding_failures, 0U);
  }

  // Exceeding capacity: 2e + s > 8 -> strictly positive FER.
  const std::vector<std::pair<size_t, size_t>> over_capacity = {
      {5, 0},
      {4, 1},
      {3, 4},
      {0, 9},
  };
  for (const auto &[e, s] : over_capacity) {
    channel::ErrorErasureChannel<GF2XLOG<uint8_t>> ch(gf32, 31, e, s,
                                                      6000 + e * 10 + s);
    const auto stats = channel::SimulateCodePerformance(
        rs, ch, channel::SimulationConfig{25, 0, 6100 + e * 10 + s});
    EXPECT_GT(stats.FER(), 0.0)
        << "Expected positive FER for over-capacity e=" << e << ", s=" << s;
  }
}

}  // namespace
}  // namespace xg
