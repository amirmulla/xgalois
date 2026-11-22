#ifndef XGALOIS_RING_IDEAL_HPP_
#define XGALOIS_RING_IDEAL_HPP_

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"

namespace xg {
namespace ring {

namespace detail {

template <EuclideanDomain R>
inline typename R::element_type CanonicalizeElement(
    const R& ring, const typename R::element_type& a) {
  if constexpr (requires { alg::AssociateCanonical(ring, a); }) {
    return alg::AssociateCanonical(ring, a).first;
  } else {
    return a;
  }
}

template <EuclideanDomain R>
inline bool DividesInEd(const R& ring, const typename R::element_type& d,
                        const typename R::element_type& x) {
  if (ring.IsZero(d)) {
    return ring.IsZero(x);
  }
  return ring.IsZero(alg::DivRem(ring, x, d).second);
}

// Generic Extended Euclidean Algorithm over any xg::EuclideanDomain R.
// Returns {g, s, t} with g canonicalized (when AssociateCanonical is
// available) and s * a + t * b == g.
template <EuclideanDomain R>
inline std::tuple<typename R::element_type, typename R::element_type,
                  typename R::element_type>
EuclideanExtendedGcd(const R& ring, const typename R::element_type& a,
                     const typename R::element_type& b) {
  using E = typename R::element_type;
  E old_r = a;
  E r = b;
  E old_s = ring.MultiplicativeIdentity();
  E s = ring.AdditiveIdentity();
  E old_t = ring.AdditiveIdentity();
  E t = ring.MultiplicativeIdentity();

  while (!ring.IsZero(r)) {
    auto [q, rem] = alg::DivRem(ring, old_r, r);
    old_r = std::move(r);
    r = std::move(rem);

    E next_s = ring.Sub(old_s, ring.Mul(q, s));
    old_s = std::move(s);
    s = std::move(next_s);

    E next_t = ring.Sub(old_t, ring.Mul(q, t));
    old_t = std::move(t);
    t = std::move(next_t);
  }

  if constexpr (requires { alg::AssociateCanonical(ring, old_r); }) {
    auto [g_canon, unit] = alg::AssociateCanonical(ring, old_r);
    return {std::move(g_canon), ring.Mul(unit, old_s), ring.Mul(unit, old_t)};
  } else {
    return {std::move(old_r), std::move(old_s), std::move(old_t)};
  }
}

template <EuclideanDomain R>
inline typename R::element_type EuclideanLcm(
    const R& ring, const typename R::element_type& a,
    const typename R::element_type& b) {
  if (ring.IsZero(a) || ring.IsZero(b)) {
    return ring.AdditiveIdentity();
  }
  const auto [g, unused_s, unused_t] = EuclideanExtendedGcd(ring, a, b);
  (void)unused_s;
  (void)unused_t;
  const auto q = alg::DivRem(ring, a, g).first;
  return CanonicalizeElement(ring, ring.Mul(q, b));
}

}  // namespace detail

// Principal ideal (g) in any xg::EuclideanDomain R, with generator
// canonicalized via xg::alg::AssociateCanonical.
template <EuclideanDomain R>
class PrincipalIdeal {
 public:
  using ring_type = R;
  using element_type = typename R::element_type;

  PrincipalIdeal(R ring, element_type generator)
      : ring_(std::move(ring)),
        gen_(detail::CanonicalizeElement(ring_, generator)) {}

  const R& Ring() const noexcept { return ring_; }

  const element_type& Generator() const noexcept { return gen_; }

  bool Contains(const element_type& x) const {
    return detail::DividesInEd(ring_, gen_, x);
  }

  bool IsZero() const { return ring_.IsZero(gen_); }

  bool IsWholeRing() const {
    return detail::DividesInEd(ring_, gen_, ring_.MultiplicativeIdentity());
  }

  bool operator==(const PrincipalIdeal& other) const {
    if constexpr (requires { alg::AssociateCanonical(ring_, gen_); }) {
      return alg::Equals(ring_, gen_, other.gen_);
    } else {
      return detail::DividesInEd(ring_, gen_, other.gen_) &&
             detail::DividesInEd(ring_, other.gen_, gen_);
    }
  }

  bool operator!=(const PrincipalIdeal& other) const {
    return !(*this == other);
  }

 private:
  R ring_;
  element_type gen_;
};

// Ideal sum: (a) + (b) = (gcd(a, b)).
template <EuclideanDomain R>
inline PrincipalIdeal<R> IdealSum(const PrincipalIdeal<R>& lhs,
                                  const PrincipalIdeal<R>& rhs) {
  auto [g, unused_s, unused_t] = detail::EuclideanExtendedGcd(
      lhs.Ring(), lhs.Generator(), rhs.Generator());
  (void)unused_s;
  (void)unused_t;
  return PrincipalIdeal<R>(lhs.Ring(), std::move(g));
}

// Ideal product: (a) * (b) = (a * b).
template <EuclideanDomain R>
inline PrincipalIdeal<R> IdealProduct(const PrincipalIdeal<R>& lhs,
                                      const PrincipalIdeal<R>& rhs) {
  return PrincipalIdeal<R>(lhs.Ring(),
                           lhs.Ring().Mul(lhs.Generator(), rhs.Generator()));
}

// Ideal intersection: (a) \cap (b) = (lcm(a, b)).
template <EuclideanDomain R>
inline PrincipalIdeal<R> IdealIntersection(const PrincipalIdeal<R>& lhs,
                                           const PrincipalIdeal<R>& rhs) {
  return PrincipalIdeal<R>(
      lhs.Ring(),
      detail::EuclideanLcm(lhs.Ring(), lhs.Generator(), rhs.Generator()));
}

// Generalized Chinese Remainder Theorem over ANY xg::EuclideanDomain R.
// Solves x == residues[i] (mod moduli[i]) for all i.
// Handles BOTH pairwise-coprime and non-coprime moduli:
// - Returns {solution, lcm_modulus} when the system is consistent.
// - Returns std::nullopt iff the system is inconsistent.
template <EuclideanDomain R>
inline std::optional<
    std::pair<typename R::element_type, typename R::element_type>>
ChineseRemainderTheorem(const R& ring,
                        const std::vector<typename R::element_type>& residues,
                        const std::vector<typename R::element_type>& moduli) {
  using E = typename R::element_type;
  if (residues.size() != moduli.size()) {
    throw std::invalid_argument(
        "ChineseRemainderTheorem: residues and moduli must have equal length");
  }
  if (residues.empty()) {
    return std::pair<E, E>{ring.AdditiveIdentity(),
                           ring.MultiplicativeIdentity()};
  }

  for (const auto& m : moduli) {
    if (ring.IsZero(m)) {
      throw std::invalid_argument(
          "ChineseRemainderTheorem: modulus cannot be zero");
    }
  }

  E cur_mod = detail::CanonicalizeElement(ring, moduli[0]);
  E cur_rem = alg::DivRem(ring, residues[0], cur_mod).second;

  for (std::size_t i = 1; i < residues.size(); ++i) {
    const E m2 = detail::CanonicalizeElement(ring, moduli[i]);
    const E r2 = alg::DivRem(ring, residues[i], m2).second;

    auto [g, s, unused_t] = detail::EuclideanExtendedGcd(ring, cur_mod, m2);
    (void)unused_t;

    const E diff = ring.Sub(r2, cur_rem);
    auto [k, check_rem] = alg::DivRem(ring, diff, g);
    if (!ring.IsZero(check_rem)) {
      return std::nullopt;
    }

    const E m2_div_g = alg::DivRem(ring, m2, g).first;
    E new_mod = detail::CanonicalizeElement(ring, ring.Mul(cur_mod, m2_div_g));

    // Reduce step = (k * s) mod (m2 / g) to prevent intermediate blowup.
    const E step = alg::DivRem(ring, ring.Mul(k, s), m2_div_g).second;
    E new_rem =
        alg::DivRem(ring, ring.Add(cur_rem, ring.Mul(cur_mod, step)), new_mod)
            .second;

    cur_rem = std::move(new_rem);
    cur_mod = std::move(new_mod);
  }

  return std::pair<E, E>{std::move(cur_rem), std::move(cur_mod)};
}

}  // namespace ring
}  // namespace xg

#endif  // XGALOIS_RING_IDEAL_HPP_
