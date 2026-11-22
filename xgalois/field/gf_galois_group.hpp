#ifndef XGALOIS_FIELD_GF_GALOIS_GROUP_HPP_
#define XGALOIS_FIELD_GF_GALOIS_GROUP_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {

namespace detail {

template <FiniteField F>
typename F::element_type ApplySingleFrobenius(
    const F& field, const typename F::element_type& a) {
  if constexpr (requires(const F& f, const typename F::element_type& x) {
                  { f.Frobenius(x) } -> std::same_as<typename F::element_type>;
                }) {
    return field.Frobenius(a);
  } else {
    return field.Pow(a, static_cast<uint64_t>(field.Characteristic()));
  }
}

template <FiniteField F>
typename F::element_type EmbedPrimeScalar(const F& field, uint64_t value) {
  const uint64_t p = static_cast<uint64_t>(field.Characteristic());
  value %= p;
  if (value == 0) {
    return field.AdditiveIdentity();
  }
  if (value == 1) {
    return field.MultiplicativeIdentity();
  }
  // Repeated doubling of 1_F in O(log p) field additions.
  typename F::element_type acc = field.AdditiveIdentity();
  typename F::element_type cur = field.MultiplicativeIdentity();
  while (value > 0) {
    if ((value & 1ULL) != 0) {
      acc = field.Add(acc, cur);
    }
    value >>= 1;
    if (value > 0) {
      cur = field.Add(cur, cur);
    }
  }
  return acc;
}

template <FiniteField F>
uint64_t ExtractPrimeScalar(const F& field,
                            const typename F::element_type& elem) {
  using E = typename F::element_type;
  if (field.IsZero(elem)) {
    return 0;
  }
  if (field.IsOne(elem)) {
    return 1;
  }
  const uint64_t p = static_cast<uint64_t>(field.Characteristic());
  if (p == 2) {
    throw std::domain_error("Element does not lie in the prime subfield GF(2)");
  }
  if constexpr (std::is_integral_v<E>) {
    if (field.Degree() == 1) {
      return static_cast<uint64_t>(field.GetElementValue(elem)) % p;
    }
  }
  if constexpr (requires(const E& e) { e[std::size_t{0}]; }) {
    const E canonical = field.GetElementValue(elem);
    return static_cast<uint64_t>(canonical[0]) % p;
  }
  // General fallback: scan multiples of 1_F.
  const E one = field.MultiplicativeIdentity();
  E cur = one;
  for (uint64_t k = 1; k < p; ++k) {
    if (cur == elem) {
      return k;
    }
    cur = field.Add(cur, one);
  }
  throw std::domain_error("Element does not lie in the prime subfield GF(p)");
}

// Extracts F_p coordinates (c_0, ..., c_{m-1}) of an element in a polynomial-
// basis FiniteField representation.
template <FiniteField F>
std::vector<uint64_t> ExtractPolyBasisCoords(
    const F& field, const typename F::element_type& elem) {
  using E = typename F::element_type;
  const uint32_t m = field.Degree();
  const uint64_t p = static_cast<uint64_t>(field.Characteristic());
  std::vector<uint64_t> coords(m, 0);
  const E val = field.GetElementValue(elem);

  if (m == 1) {
    coords[0] = ExtractPrimeScalar(field, val);
    return coords;
  }

  if constexpr (std::is_integral_v<E>) {
    // Binary extension GF(2^m) packed in an unsigned word.
    const uint64_t bits = static_cast<uint64_t>(val);
    for (uint32_t i = 0; i < m; ++i) {
      coords[i] = (bits >> i) & 1ULL;
    }
    return coords;
  } else if constexpr (requires(const E& e, int i) {
                         { e.Bit(i) } -> std::convertible_to<bool>;
                       }) {
    // 128-bit binary extension (UInt128).
    for (uint32_t i = 0; i < m; ++i) {
      coords[i] = val.Bit(static_cast<int>(i)) ? 1ULL : 0ULL;
    }
    return coords;
  } else if constexpr (requires(const E& e, std::size_t i) { e[i]; }) {
    // General extension GF(p^m) (ExtensionElement).
    for (uint32_t i = 0; i < m; ++i) {
      coords[i] = static_cast<uint64_t>(val[i]) % p;
    }
    return coords;
  } else {
    throw std::invalid_argument(
        "FieldIsomorphism: unsupported element representation for coordinate "
        "extraction");
  }
}

// Reconstructs a field element from its F_p polynomial-basis coordinates
// (c_0, ..., c_{m-1}).
template <FiniteField F>
typename F::element_type FromPolyBasisCoords(
    const F& field, const std::vector<uint64_t>& coords) {
  using E = typename F::element_type;
  const uint32_t m = field.Degree();
  const uint64_t p = static_cast<uint64_t>(field.Characteristic());

  if (m == 1) {
    return EmbedPrimeScalar(field, coords.empty() ? 0ULL : coords[0]);
  }

  if constexpr (std::is_integral_v<E>) {
    uint64_t bits = 0;
    for (uint32_t i = 0; i < m && i < coords.size(); ++i) {
      if ((coords[i] & 1ULL) != 0) {
        bits |= (uint64_t{1} << i);
      }
    }
    return field.SetElementValue(static_cast<E>(bits));
  } else if constexpr (requires(E e) {
                         e.lo;
                         e.hi;
                       }) {
    E out{};
    for (uint32_t i = 0; i < m && i < coords.size(); ++i) {
      if ((coords[i] & 1ULL) != 0) {
        if (i < 64) {
          out.lo |= (uint64_t{1} << i);
        } else {
          out.hi |= (uint64_t{1} << (i - 64));
        }
      }
    }
    return field.SetElementValue(out);
  } else if constexpr (requires {
                         typename E::coeff_type;
                         E(std::vector<typename E::coeff_type>{});
                       }) {
    using C = typename E::coeff_type;
    std::vector<C> vec(m, C{0});
    for (uint32_t i = 0; i < m && i < coords.size(); ++i) {
      vec[i] = static_cast<C>(coords[i] % p);
    }
    return field.SetElementValue(E(std::move(vec)));
  } else {
    throw std::invalid_argument(
        "FieldIsomorphism: unsupported element representation for coordinate "
        "reconstruction");
  }
}

// Extracts the defining monic irreducible polynomial f(x) in F_p[x] (length
// m + 1 in ascending degree order) of a polynomial-basis FiniteField.
template <FiniteField F>
std::vector<uint64_t> ExtractDefiningPoly(const F& field) {
  const uint32_t m = field.Degree();
  const uint64_t p = static_cast<uint64_t>(field.Characteristic());

  if (m == 1) {
    return {0ULL, 1ULL};
  }

  if constexpr (requires(const F& f) { f.ModulusLow(); }) {
    const auto low = field.ModulusLow();
    std::vector<uint64_t> poly(static_cast<std::size_t>(m) + 1, 0);
    if constexpr (std::is_integral_v<std::remove_cvref_t<decltype(low)>>) {
      const uint64_t bits = static_cast<uint64_t>(low);
      for (uint32_t i = 0; i < m; ++i) {
        poly[i] = (bits >> i) & 1ULL;
      }
    } else {
      // UInt128
      for (uint32_t i = 0; i < m; ++i) {
        poly[i] = low.Bit(static_cast<int>(i)) ? 1ULL : 0ULL;
      }
    }
    poly[m] = 1ULL;
    return poly;
  } else if constexpr (requires(const F& f) { f.Modulus(); }) {
    const auto& mod = field.Modulus();
    using M = std::remove_cvref_t<decltype(mod)>;
    std::vector<uint64_t> poly(static_cast<std::size_t>(m) + 1, 0);
    if constexpr (requires(const M& x, std::size_t i) { x[i]; }) {
      const uint64_t lc = static_cast<uint64_t>(mod[m]) % p;
      const uint64_t inv_lc = (lc == 1) ? 1ULL : utils::ModInv(lc, p);
      for (uint32_t i = 0; i <= m; ++i) {
        poly[i] = utils::ModMul(static_cast<uint64_t>(mod[i]) % p, inv_lc, p);
      }
      return poly;
    } else if constexpr (std::is_integral_v<M>) {
      const uint64_t bits = static_cast<uint64_t>(mod);
      for (uint32_t i = 0; i < m; ++i) {
        poly[i] = (bits >> i) & 1ULL;
      }
      poly[m] = 1ULL;
      return poly;
    } else {
      throw std::invalid_argument(
          "FieldIsomorphism: cannot extract defining polynomial from Modulus");
    }
  } else {
    throw std::invalid_argument(
        "FieldIsomorphism: source field does not expose its defining modulus");
  }
}

// Inverts an m x m matrix over F_p (stored row-major) using Gauss-Jordan
// elimination. Throws std::runtime_error if singular.
inline std::vector<std::vector<uint64_t>> InvertMatrixModP(
    std::vector<std::vector<uint64_t>> mat, uint64_t p) {
  const std::size_t m = mat.size();
  std::vector<std::vector<uint64_t>> inv(m, std::vector<uint64_t>(m, 0));
  for (std::size_t i = 0; i < m; ++i) {
    inv[i][i] = 1;
  }

  for (std::size_t col = 0; col < m; ++col) {
    std::size_t pivot = m;
    for (std::size_t row = col; row < m; ++row) {
      if (mat[row][col] % p != 0) {
        pivot = row;
        break;
      }
    }
    if (pivot == m) {
      throw std::runtime_error(
          "FieldIsomorphism: change-of-basis matrix is singular");
    }
    if (pivot != col) {
      std::swap(mat[pivot], mat[col]);
      std::swap(inv[pivot], inv[col]);
    }

    const uint64_t diag = mat[col][col] % p;
    if (diag != 1) {
      const uint64_t inv_diag = utils::ModInv(diag, p);
      for (std::size_t j = 0; j < m; ++j) {
        mat[col][j] = utils::ModMul(mat[col][j], inv_diag, p);
        inv[col][j] = utils::ModMul(inv[col][j], inv_diag, p);
      }
    }

    for (std::size_t row = 0; row < m; ++row) {
      if (row == col) {
        continue;
      }
      const uint64_t factor = mat[row][col] % p;
      if (factor == 0) {
        continue;
      }
      for (std::size_t j = 0; j < m; ++j) {
        const uint64_t sub_m = utils::ModMul(factor, mat[col][j], p);
        mat[row][j] = (mat[row][j] >= sub_m) ? (mat[row][j] - sub_m)
                                             : (p - (sub_m - mat[row][j]));
        const uint64_t sub_i = utils::ModMul(factor, inv[col][j], p);
        inv[row][j] = (inv[row][j] >= sub_i) ? (inv[row][j] - sub_i)
                                             : (p - (sub_i - inv[row][j]));
      }
    }
  }
  return inv;
}

// Lightweight polynomial helpers over an arbitrary const FiniteField& for
// Cantor-Zassenhaus root finding in FieldIsomorphism without requiring
// std::shared_ptr ownership.
template <FiniteField F>
struct PolyOps {
  using E = typename F::element_type;
  using Poly = std::vector<E>;

  const F& field;

  void Trim(Poly& a) const {
    while (!a.empty() && field.IsZero(a.back())) {
      a.pop_back();
    }
  }

  Poly Monic(Poly a) const {
    Trim(a);
    if (a.empty()) {
      return a;
    }
    const E lc = a.back();
    if (!field.IsOne(lc)) {
      const E inv_lc = field.Inv(lc);
      for (E& c : a) {
        c = field.Mul(c, inv_lc);
      }
    }
    return a;
  }

  Poly Add(const Poly& a, const Poly& b) const {
    const std::size_t n = std::max(a.size(), b.size());
    Poly out(n, field.AdditiveIdentity());
    for (std::size_t i = 0; i < n; ++i) {
      const E av = (i < a.size()) ? a[i] : field.AdditiveIdentity();
      const E bv = (i < b.size()) ? b[i] : field.AdditiveIdentity();
      out[i] = field.Add(av, bv);
    }
    Trim(out);
    return out;
  }

  Poly Sub(const Poly& a, const Poly& b) const {
    const std::size_t n = std::max(a.size(), b.size());
    Poly out(n, field.AdditiveIdentity());
    for (std::size_t i = 0; i < n; ++i) {
      const E av = (i < a.size()) ? a[i] : field.AdditiveIdentity();
      const E bv = (i < b.size()) ? b[i] : field.AdditiveIdentity();
      out[i] = field.Sub(av, bv);
    }
    Trim(out);
    return out;
  }

  Poly Mul(const Poly& a, const Poly& b) const {
    if (a.empty() || b.empty()) {
      return {};
    }
    Poly out(a.size() + b.size() - 1, field.AdditiveIdentity());
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (field.IsZero(a[i])) {
        continue;
      }
      for (std::size_t j = 0; j < b.size(); ++j) {
        if (field.IsZero(b[j])) {
          continue;
        }
        out[i + j] = field.Add(out[i + j], field.Mul(a[i], b[j]));
      }
    }
    Trim(out);
    return out;
  }

  std::pair<Poly, Poly> DivRem(const Poly& a, const Poly& b) const {
    if (b.empty()) {
      throw std::domain_error("PolyOps::DivRem: division by zero");
    }
    if (a.size() < b.size()) {
      return {{}, a};
    }
    Poly rem = a;
    const std::size_t d = b.size() - 1;
    const E inv_lc = field.Inv(b[d]);
    Poly q(a.size() - d, field.AdditiveIdentity());
    for (std::size_t i = a.size(); i-- > d;) {
      if (field.IsZero(rem[i])) {
        continue;
      }
      const E qc = field.Mul(rem[i], inv_lc);
      const std::size_t shift = i - d;
      q[shift] = qc;
      for (std::size_t j = 0; j <= d; ++j) {
        if (!field.IsZero(b[j])) {
          rem[shift + j] = field.Sub(rem[shift + j], field.Mul(qc, b[j]));
        }
      }
    }
    Trim(q);
    Trim(rem);
    return {std::move(q), std::move(rem)};
  }

  Poly Mod(const Poly& a, const Poly& m) const { return DivRem(a, m).second; }

  Poly MulMod(const Poly& a, const Poly& b, const Poly& m) const {
    return Mod(Mul(a, b), m);
  }

  Poly PowMod(Poly base, uint64_t exp, const Poly& mod) const {
    Poly res = {field.MultiplicativeIdentity()};
    base = Mod(base, mod);
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        res = MulMod(res, base, mod);
      }
      exp >>= 1;
      if (exp > 0) {
        base = MulMod(base, base, mod);
      }
    }
    return res;
  }

  Poly Gcd(Poly a, Poly b) const {
    Trim(a);
    Trim(b);
    while (!b.empty()) {
      Poly r = Mod(a, b);
      a = std::move(b);
      b = std::move(r);
    }
    return Monic(std::move(a));
  }

  E Eval(const Poly& a, const E& x) const {
    E acc = field.AdditiveIdentity();
    for (std::size_t i = a.size(); i-- > 0;) {
      acc = field.Add(field.Mul(acc, x), a[i]);
    }
    return acc;
  }
};

// Finds a single root in FieldB of a monic irreducible degree-m polynomial
// f_poly in F_p[x] (which splits completely into m distinct linear factors over
// FieldB ~ F_{p^m}).
template <FiniteField FieldB>
typename FieldB::element_type FindRootOfSplittingPoly(
    const FieldB& field_b, const std::vector<uint64_t>& f_poly) {
  using E = typename FieldB::element_type;
  PolyOps<FieldB> ops{field_b};

  typename PolyOps<FieldB>::Poly u;
  u.reserve(f_poly.size());
  for (uint64_t c : f_poly) {
    u.push_back(EmbedPrimeScalar(field_b, c));
  }
  ops.Trim(u);
  u = ops.Monic(std::move(u));

  if (u.size() < 2) {
    throw std::invalid_argument(
        "FieldIsomorphism: defining polynomial must have degree >= 1");
  }
  if (u.size() == 2) {
    return field_b.Neg(u[0]);
  }

  const uint64_t p = static_cast<uint64_t>(field_b.Characteristic());
  const uint64_t q = static_cast<uint64_t>(field_b.Order());
  const uint32_t m = field_b.Degree();

  // Fast direct scan for small integral binary fields (e.g. GF(2^8)).
  if constexpr (std::is_integral_v<E>) {
    if (q > 0 && q <= 1024 && p == 2) {
      for (uint64_t cand = 0; cand < q; ++cand) {
        const E x = field_b.SetElementValue(static_cast<E>(cand));
        if (field_b.IsZero(ops.Eval(u, x))) {
          return x;
        }
      }
    }
  }

  // Cantor-Zassenhaus splitting of a completely split polynomial u in F_q[X].
  std::mt19937_64 rng(0x9E3779B97F4A7C15ULL ^ (p << 16) ^ m);
  auto random_elem = [&]() -> E {
    if constexpr (requires(const FieldB& f, std::mt19937_64& g) {
                    { f.Random(g) } -> std::same_as<E>;
                  }) {
      return field_b.Random(rng);
    } else {
      std::vector<uint64_t> coords(m);
      std::uniform_int_distribution<uint64_t> dist(0, p - 1);
      for (uint32_t i = 0; i < m; ++i) {
        coords[i] = dist(rng);
      }
      return FromPolyBasisCoords(field_b, coords);
    }
  };

  const typename PolyOps<FieldB>::Poly one_poly = {
      field_b.MultiplicativeIdentity()};

  constexpr std::size_t kMaxAttempts = 10000;
  while (u.size() > 2) {
    bool split = false;
    for (std::size_t attempt = 0; attempt < kMaxAttempts; ++attempt) {
      // Choose random non-constant linear polynomial a(X) = r1 * X + r0
      // (nonzero r1 is essential in characteristic 2 so that Tr(r1 * root)
      // separates Frobenius conjugates).
      E r1 = field_b.AdditiveIdentity();
      while (field_b.IsZero(r1)) {
        r1 = random_elem();
      }
      const E r0 = random_elem();
      typename PolyOps<FieldB>::Poly a = {r0, r1};

      typename PolyOps<FieldB>::Poly g = ops.Gcd(a, u);
      if (g.size() > 1 && g.size() < u.size()) {
        auto q_poly = ops.DivRem(u, g).first;
        u = (g.size() <= q_poly.size()) ? std::move(g) : std::move(q_poly);
        split = true;
        break;
      }

      typename PolyOps<FieldB>::Poly splitter;
      if (p == 2) {
        // Absolute trace Tr_{F_{2^m}/F_2}(a) mod u = sum_{j=0}^{m-1} a^(2^j).
        typename PolyOps<FieldB>::Poly cur = ops.Mod(a, u);
        splitter = cur;
        for (uint32_t j = 1; j < m; ++j) {
          cur = ops.MulMod(cur, cur, u);
          splitter = ops.Add(splitter, cur);
        }
      } else {
        // Odd characteristic: a^((q - 1)/2) - 1 mod u.
        const uint64_t half = (q - 1) >> 1;
        splitter = ops.Sub(ops.PowMod(a, half, u), one_poly);
      }

      g = ops.Gcd(u, splitter);
      if (g.size() > 1 && g.size() < u.size()) {
        auto q_poly = ops.DivRem(u, g).first;
        u = (g.size() <= q_poly.size()) ? std::move(g) : std::move(q_poly);
        split = true;
        break;
      }
    }
    if (!split) {
      throw std::runtime_error(
          "FieldIsomorphism: Cantor-Zassenhaus failed to find a root");
    }
  }

  return field_b.Neg(u[0]);
}

}  // namespace detail

// ============================================================================
// GaloisGroup: Gal(F_{p^m} / F_{p^d}) (cyclic of order m / d)
// ============================================================================

template <FiniteField F>
class GaloisGroup {
 public:
  using field_type = F;
  using element_type = typename F::element_type;

  explicit GaloisGroup(const F& field, uint32_t subfield_degree = 1)
      : field_(&field), subfield_degree_(subfield_degree) {
    const uint32_t m = field.Degree();
    if (subfield_degree_ == 0 || m % subfield_degree_ != 0) {
      throw std::invalid_argument(
          "GaloisGroup: subfield_degree must divide field degree m");
    }
    order_ = m / subfield_degree_;
  }

  uint32_t Order() const noexcept { return order_; }
  uint32_t SubfieldDegree() const noexcept { return subfield_degree_; }
  const F& Field() const noexcept { return *field_; }

  // Applies the k-th power of the relative Frobenius automorphism:
  // sigma^k(a) = a^(q^k) where q = p^(subfield_degree).
  element_type Apply(const element_type& elem, uint64_t k) const {
    const uint64_t reduced_k = k % static_cast<uint64_t>(order_);
    if (reduced_k == 0) {
      return elem;
    }
    const uint64_t total_steps =
        reduced_k * static_cast<uint64_t>(subfield_degree_);
    element_type cur = elem;
    for (uint64_t step = 0; step < total_steps; ++step) {
      cur = detail::ApplySingleFrobenius(*field_, cur);
    }
    return cur;
  }

  // Returns the distinct conjugates {a, sigma(a), sigma^2(a), ...} in the
  // orbit of elem under Gal(F_{p^m} / F_{p^d}).
  std::vector<element_type> FrobeniusOrbit(const element_type& elem) const {
    std::vector<element_type> orbit;
    orbit.reserve(order_);
    element_type cur = elem;
    for (uint32_t i = 0; i < order_; ++i) {
      orbit.push_back(cur);
      cur = Apply(cur, 1);
      if (cur == elem) {
        break;
      }
    }
    return orbit;
  }

  // Returns the monic minimal polynomial m_a(x) in F_p[x] (coefficients in
  // [0, p - 1] in ascending degree order) of elem over the prime field F_p.
  std::vector<uint64_t> MinimalPolynomialOverPrime(
      const element_type& elem) const {
    const uint32_t m = field_->Degree();
    std::vector<element_type> conjugates;
    conjugates.reserve(m);
    element_type cur = elem;
    for (uint32_t i = 0; i < m; ++i) {
      conjugates.push_back(cur);
      cur = detail::ApplySingleFrobenius(*field_, cur);
      if (cur == elem) {
        break;
      }
    }

    // Form prod_{c in conjugates} (x - c) in F[x].
    std::vector<element_type> poly;
    poly.reserve(conjugates.size() + 1);
    poly.push_back(field_->MultiplicativeIdentity());
    for (const element_type& c : conjugates) {
      std::vector<element_type> next(poly.size() + 1,
                                     field_->AdditiveIdentity());
      const element_type neg_c = field_->Neg(c);
      for (std::size_t i = 0; i < poly.size(); ++i) {
        next[i] = field_->Add(next[i], field_->Mul(poly[i], neg_c));
        next[i + 1] = field_->Add(next[i + 1], poly[i]);
      }
      poly = std::move(next);
    }

    std::vector<uint64_t> coeffs;
    coeffs.reserve(poly.size());
    for (const element_type& c : poly) {
      coeffs.push_back(detail::ExtractPrimeScalar(*field_, c));
    }
    return coeffs;
  }

 private:
  const F* field_;
  uint32_t subfield_degree_;
  uint32_t order_;
};

// ============================================================================
// CyclotomicCosets: partition Z/nZ under x -> q*x mod n (gcd(n, q) == 1)
// ============================================================================

inline std::vector<std::vector<uint64_t>> CyclotomicCosets(uint64_t n,
                                                           uint64_t q) {
  if (n == 0) {
    throw std::invalid_argument("CyclotomicCosets: n must be >= 1");
  }
  if (q < 2) {
    throw std::invalid_argument("CyclotomicCosets: q must be >= 2");
  }
  if (std::gcd(n, q) != 1) {
    throw std::invalid_argument("CyclotomicCosets: requires gcd(n, q) == 1");
  }

  __extension__ using u128 = unsigned __int128;
  std::vector<bool> visited(static_cast<std::size_t>(n), false);
  std::vector<std::vector<uint64_t>> cosets;

  for (uint64_t s = 0; s < n; ++s) {
    if (visited[static_cast<std::size_t>(s)]) {
      continue;
    }
    std::vector<uint64_t> coset;
    uint64_t cur = s;
    while (!visited[static_cast<std::size_t>(cur)]) {
      visited[static_cast<std::size_t>(cur)] = true;
      coset.push_back(cur);
      cur = static_cast<uint64_t>((static_cast<u128>(cur) * q) % n);
    }
    cosets.push_back(std::move(coset));
  }
  return cosets;
}

// ============================================================================
// FieldIsomorphism: Explicit isomorphism phi: FieldA -> FieldB between two
// representations of F_{p^m} with different defining irreducible polynomials.
// ============================================================================

template <FiniteField FieldA, FiniteField FieldB>
class FieldIsomorphism {
 public:
  using elem_a_type = typename FieldA::element_type;
  using elem_b_type = typename FieldB::element_type;

  FieldIsomorphism(const FieldA& field_a, const FieldB& field_b)
      : field_a_(&field_a),
        field_b_(&field_b),
        p_(static_cast<uint64_t>(field_a.Characteristic())),
        m_(field_a.Degree()) {
    if (p_ != static_cast<uint64_t>(field_b.Characteristic()) ||
        m_ != field_b.Degree()) {
      throw std::invalid_argument(
          "FieldIsomorphism: fields must have identical characteristic and "
          "degree");
    }

    const std::vector<uint64_t> def_poly_a =
        detail::ExtractDefiningPoly(*field_a_);
    beta_ = detail::FindRootOfSplittingPoly(*field_b_, def_poly_a);

    // Precompute powers beta^0, beta^1, ..., beta^{m-1} in FieldB.
    beta_powers_.reserve(m_);
    elem_b_type cur = field_b_->MultiplicativeIdentity();
    for (uint32_t i = 0; i < m_; ++i) {
      beta_powers_.push_back(cur);
      if (i + 1 < m_) {
        cur = field_b_->Mul(cur, beta_);
      }
    }

    // Build m x m F_p change-of-basis matrix M whose j-th column is the
    // coordinate vector of beta^j in FieldB, then invert M once for Inverse().
    std::vector<std::vector<uint64_t>> mat(m_, std::vector<uint64_t>(m_, 0));
    for (uint32_t j = 0; j < m_; ++j) {
      const std::vector<uint64_t> col =
          detail::ExtractPolyBasisCoords(*field_b_, beta_powers_[j]);
      for (uint32_t i = 0; i < m_; ++i) {
        mat[i][j] = col[i] % p_;
      }
    }
    inv_matrix_ = detail::InvertMatrixModP(std::move(mat), p_);
  }

  elem_b_type RootImage() const { return beta_; }

  elem_b_type Forward(const elem_a_type& a) const {
    const std::vector<uint64_t> coords =
        detail::ExtractPolyBasisCoords(*field_a_, a);
    elem_b_type acc = field_b_->AdditiveIdentity();
    for (uint32_t i = 0; i < m_; ++i) {
      if (coords[i] == 0) {
        continue;
      }
      if (coords[i] == 1) {
        acc = field_b_->Add(acc, beta_powers_[i]);
      } else {
        const elem_b_type scalar =
            detail::EmbedPrimeScalar(*field_b_, coords[i]);
        acc = field_b_->Add(acc, field_b_->Mul(scalar, beta_powers_[i]));
      }
    }
    return acc;
  }

  elem_a_type Inverse(const elem_b_type& b) const {
    const std::vector<uint64_t> b_coords =
        detail::ExtractPolyBasisCoords(*field_b_, b);
    std::vector<uint64_t> a_coords(m_, 0);
    for (uint32_t i = 0; i < m_; ++i) {
      uint64_t sum = 0;
      for (uint32_t j = 0; j < m_; ++j) {
        if (b_coords[j] != 0 && inv_matrix_[i][j] != 0) {
          sum = (sum + utils::ModMul(inv_matrix_[i][j], b_coords[j], p_)) % p_;
        }
      }
      a_coords[i] = sum;
    }
    return detail::FromPolyBasisCoords(*field_a_, a_coords);
  }

 private:
  const FieldA* field_a_;
  const FieldB* field_b_;
  uint64_t p_;
  uint32_t m_;
  elem_b_type beta_{};
  std::vector<elem_b_type> beta_powers_;
  std::vector<std::vector<uint64_t>> inv_matrix_;
};

// ============================================================================
// SplittingFieldInfo & AnalyzeSplittingField
// ============================================================================

struct SplittingFieldInfo {
  uint64_t base_order = 0;
  uint32_t extension_degree = 1;
  std::vector<uint32_t> factor_degrees;
};

// Analyzes the splitting field of a nonzero polynomial f in F_q[x] using
// square-free + distinct-degree factorization to compute the exact extension
// degree [L : F_q] = lcm(deg(g_i)) over all irreducible factors g_i of f.
template <typename F>
SplittingFieldInfo AnalyzeSplittingField(const PolynomialDense<F>& f) {
  if (!f.HasField()) {
    throw std::invalid_argument(
        "AnalyzeSplittingField: polynomial must have an associated field");
  }
  if (f.IsZero()) {
    throw std::invalid_argument(
        "AnalyzeSplittingField: zero polynomial has no splitting field");
  }

  SplittingFieldInfo info;
  info.base_order = static_cast<uint64_t>(f.Field()->Order());
  info.extension_degree = 1;

  if (f.Degree() == 0) {
    return info;
  }

  const auto sqfree = utils::SquareFreeFactor(f.Monic());
  for (const auto& [sq_poly, mult] : sqfree) {
    const auto ddf = utils::DistinctDegreeFactor(sq_poly);
    for (const auto& [dd_poly, deg] : ddf) {
      const uint32_t count = static_cast<uint32_t>(dd_poly.Degree()) / deg;
      for (uint32_t r = 0; r < count * mult; ++r) {
        info.factor_degrees.push_back(deg);
      }
      info.extension_degree = std::lcm(info.extension_degree, deg);
    }
  }
  std::sort(info.factor_degrees.begin(), info.factor_degrees.end());
  return info;
}

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_GALOIS_GROUP_HPP_
