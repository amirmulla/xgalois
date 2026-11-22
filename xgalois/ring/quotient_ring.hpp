#ifndef XGALOIS_RING_QUOTIENT_RING_HPP_
#define XGALOIS_RING_QUOTIENT_RING_HPP_

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/ring/polynomial_ring.hpp"

namespace xg {
namespace ring {

// ============================================================================
// QuotientRing<R>: R / (m) for any EuclideanDomain R and nonzero non-unit m
// ============================================================================

template <EuclideanDomain R>
class QuotientRing {
 public:
  using base_ring_type = R;
  using element_type = typename R::element_type;

  QuotientRing(R base_ring, element_type modulus)
      : ring_(std::move(base_ring)), modulus_(std::move(modulus)) {
    if (ring_.IsZero(modulus_)) {
      throw std::invalid_argument("QuotientRing: modulus cannot be zero");
    }
    const auto [q_unit, r_unit] =
        alg::DivRem(ring_, ring_.MultiplicativeIdentity(), modulus_);
    (void)q_unit;
    if (ring_.IsZero(r_unit)) {
      throw std::invalid_argument("QuotientRing: modulus cannot be a unit");
    }
  }

  const R& BaseRing() const noexcept { return ring_; }
  const element_type& Modulus() const noexcept { return modulus_; }

  element_type Reduce(const element_type& a) const {
    return alg::DivRem(ring_, a, modulus_).second;
  }

  element_type Add(const element_type& a, const element_type& b) const {
    return Reduce(ring_.Add(a, b));
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    return Reduce(ring_.Sub(a, b));
  }

  element_type Neg(const element_type& a) const { return Reduce(ring_.Neg(a)); }

  element_type Mul(const element_type& a, const element_type& b) const {
    return Reduce(ring_.Mul(Reduce(a), Reduce(b)));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = Reduce(a);
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  element_type AdditiveIdentity() const { return ring_.AdditiveIdentity(); }

  element_type MultiplicativeIdentity() const {
    return Reduce(ring_.MultiplicativeIdentity());
  }

  bool IsZero(const element_type& a) const { return ring_.IsZero(Reduce(a)); }

  bool IsOne(const element_type& a) const {
    return alg::Equals(ring_, Reduce(a), MultiplicativeIdentity());
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return alg::Equals(ring_, Reduce(a), Reduce(b));
  }

  std::string ToString(const element_type& a) const {
    return ring_.ToString(Reduce(a));
  }

  // Computes multiplicative inverse in R / (m) via generic Euclidean EGCD.
  std::optional<element_type> TryInv(const element_type& a) const {
    element_type old_r = Reduce(a);
    if (ring_.IsZero(old_r)) {
      return std::nullopt;
    }
    element_type r = modulus_;
    element_type old_s = ring_.MultiplicativeIdentity();
    element_type s = ring_.AdditiveIdentity();

    while (!ring_.IsZero(r)) {
      auto [q, rem] = alg::DivRem(ring_, old_r, r);
      old_r = std::move(r);
      r = std::move(rem);

      element_type next_s = ring_.Sub(old_s, ring_.Mul(q, s));
      old_s = std::move(s);
      s = std::move(next_s);
    }

    if (ring_.IsZero(old_r)) {
      return std::nullopt;
    }
    // old_r is a unit in R iff old_r divides 1.
    auto [inv_g, rem_unit] =
        alg::DivRem(ring_, ring_.MultiplicativeIdentity(), old_r);
    if (!ring_.IsZero(rem_unit)) {
      return std::nullopt;
    }
    return Reduce(ring_.Mul(inv_g, old_s));
  }

  bool IsUnit(const element_type& a) const { return TryInv(a).has_value(); }

 private:
  R ring_;
  element_type modulus_;
};

// ============================================================================
// PolynomialQuotientRing<R>: R[x] / (m(x)) for monic m(x) over CommutativeRing
// ============================================================================

template <CommutativeRing R>
class PolynomialQuotientRing {
 public:
  using base_ring_type = R;
  using poly_ring_type = PolynomialRing<R>;
  using coeff_type = typename R::element_type;
  using element_type = RingPoly<R>;

  PolynomialQuotientRing(R base_ring, element_type modulus)
      : poly_ring_(std::move(base_ring)),
        modulus_(poly_ring_.Normalize(std::move(modulus))) {
    const int64_t deg = poly_ring_.Degree(modulus_);
    if (deg < 1) {
      throw std::invalid_argument(
          "PolynomialQuotientRing: modulus must have degree >= 1");
    }
    if (!poly_ring_.BaseRing().IsOne(poly_ring_.LeadingCoeff(modulus_))) {
      throw std::invalid_argument(
          "PolynomialQuotientRing: modulus must be monic");
    }
  }

  const R& BaseRing() const noexcept { return poly_ring_.BaseRing(); }
  const poly_ring_type& PolyRing() const noexcept { return poly_ring_; }
  const element_type& Modulus() const noexcept { return modulus_; }
  std::size_t Degree() const {
    return static_cast<std::size_t>(poly_ring_.Degree(modulus_));
  }

  element_type Reduce(const element_type& a) const {
    return poly_ring_.MonicDivRem(a, modulus_).second;
  }

  element_type Add(const element_type& a, const element_type& b) const {
    return Reduce(poly_ring_.Add(a, b));
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    return Reduce(poly_ring_.Sub(a, b));
  }

  element_type Neg(const element_type& a) const {
    return Reduce(poly_ring_.Neg(a));
  }

  element_type Mul(const element_type& a, const element_type& b) const {
    return Reduce(poly_ring_.Mul(Reduce(a), Reduce(b)));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = Reduce(a);
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  element_type AdditiveIdentity() const {
    return poly_ring_.AdditiveIdentity();
  }

  element_type MultiplicativeIdentity() const {
    return Reduce(poly_ring_.MultiplicativeIdentity());
  }

  bool IsZero(const element_type& a) const {
    return poly_ring_.IsZero(Reduce(a));
  }

  bool IsOne(const element_type& a) const {
    return poly_ring_.Equals(Reduce(a), MultiplicativeIdentity());
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return poly_ring_.Equals(Reduce(a), Reduce(b));
  }

  std::string ToString(const element_type& a) const {
    return poly_ring_.ToString(Reduce(a));
  }

  std::optional<element_type> TryInv(const element_type& a) const
    requires Field<R>
  {
    const element_type ar = Reduce(a);
    if (poly_ring_.IsZero(ar)) {
      return std::nullopt;
    }
    auto [g, s, unused_t] = poly_ring_.ExtendedGcd(ar, modulus_);
    (void)unused_t;
    if (!poly_ring_.IsOne(g)) {
      return std::nullopt;
    }
    return Reduce(s);
  }

  bool IsUnit(const element_type& a) const
    requires Field<R>
  {
    return TryInv(a).has_value();
  }

 private:
  poly_ring_type poly_ring_;
  element_type modulus_;
};

// ============================================================================
// NegacyclicRing<R>: R[x] / (x^N + 1) with fast negacyclic wrap-around Mul
// ============================================================================

template <CommutativeRing R>
class NegacyclicRing {
 public:
  using base_ring_type = R;
  using coeff_type = typename R::element_type;
  using element_type = std::vector<coeff_type>;

  NegacyclicRing(R base_ring, std::size_t degree)
      : ring_(std::move(base_ring)), n_(degree) {
    if (n_ == 0) {
      throw std::invalid_argument("NegacyclicRing: degree N must be >= 1");
    }
  }

  const R& BaseRing() const noexcept { return ring_; }
  std::size_t Degree() const noexcept { return n_; }

  // Reduces an arbitrary-length coefficient vector modulo (x^N + 1).
  element_type Reduce(const element_type& a) const {
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < a.size(); ++i) {
      const std::size_t idx = i % n_;
      const std::size_t wraps = i / n_;
      if ((wraps & std::size_t{1}) == 0) {
        out[idx] = ring_.Add(out[idx], a[i]);
      } else {
        out[idx] = ring_.Sub(out[idx], a[i]);
      }
    }
    return out;
  }

  element_type FromCoeffs(const element_type& coeffs) const {
    return Reduce(coeffs);
  }

  element_type Add(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      out[i] = ring_.Add(ar[i], br[i]);
    }
    return out;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      out[i] = ring_.Sub(ar[i], br[i]);
    }
    return out;
  }

  element_type Neg(const element_type& a) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      out[i] = ring_.Neg(ar[i]);
    }
    return out;
  }

  // O(N^2) negacyclic convolution: x^N = -1 wrap-around (ready for NTT hook).
  element_type Mul(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      if (ring_.IsZero(ar[i])) {
        continue;
      }
      for (std::size_t j = 0; j < n_; ++j) {
        if (ring_.IsZero(br[j])) {
          continue;
        }
        const coeff_type prod = ring_.Mul(ar[i], br[j]);
        const std::size_t s = i + j;
        if (s < n_) {
          out[s] = ring_.Add(out[s], prod);
        } else {
          out[s - n_] = ring_.Sub(out[s - n_], prod);
        }
      }
    }
    return out;
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = Reduce(a);
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  element_type AdditiveIdentity() const {
    return element_type(n_, ring_.AdditiveIdentity());
  }

  element_type MultiplicativeIdentity() const {
    element_type one(n_, ring_.AdditiveIdentity());
    one[0] = ring_.MultiplicativeIdentity();
    return one;
  }

  bool IsZero(const element_type& a) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    for (std::size_t i = 0; i < n_; ++i) {
      if (!ring_.IsZero(ar[i])) {
        return false;
      }
    }
    return true;
  }

  bool IsOne(const element_type& a) const {
    return Equals(a, MultiplicativeIdentity());
  }

  bool Equals(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    for (std::size_t i = 0; i < n_; ++i) {
      if (!alg::Equals(ring_, ar[i], br[i])) {
        return false;
      }
    }
    return true;
  }

  std::string ToString(const element_type& a) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    std::string out = "[";
    for (std::size_t i = 0; i < n_; ++i) {
      if (i > 0) {
        out += ", ";
      }
      out += ring_.ToString(ar[i]);
    }
    out += "]";
    return out;
  }

 private:
  R ring_;
  std::size_t n_;
};

// ============================================================================
// CyclicRing<R>: R[x] / (x^N - 1) with fast cyclic wrap-around Mul
// ============================================================================

template <CommutativeRing R>
class CyclicRing {
 public:
  using base_ring_type = R;
  using coeff_type = typename R::element_type;
  using element_type = std::vector<coeff_type>;

  CyclicRing(R base_ring, std::size_t degree)
      : ring_(std::move(base_ring)), n_(degree) {
    if (n_ == 0) {
      throw std::invalid_argument("CyclicRing: degree N must be >= 1");
    }
  }

  const R& BaseRing() const noexcept { return ring_; }
  std::size_t Degree() const noexcept { return n_; }

  // Reduces an arbitrary-length coefficient vector modulo (x^N - 1).
  element_type Reduce(const element_type& a) const {
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < a.size(); ++i) {
      const std::size_t idx = i % n_;
      out[idx] = ring_.Add(out[idx], a[i]);
    }
    return out;
  }

  element_type FromCoeffs(const element_type& coeffs) const {
    return Reduce(coeffs);
  }

  element_type Add(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      out[i] = ring_.Add(ar[i], br[i]);
    }
    return out;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      out[i] = ring_.Sub(ar[i], br[i]);
    }
    return out;
  }

  element_type Neg(const element_type& a) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      out[i] = ring_.Neg(ar[i]);
    }
    return out;
  }

  // O(N^2) cyclic convolution: x^N = 1 wrap-around.
  element_type Mul(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    element_type out(n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      if (ring_.IsZero(ar[i])) {
        continue;
      }
      for (std::size_t j = 0; j < n_; ++j) {
        if (ring_.IsZero(br[j])) {
          continue;
        }
        const coeff_type prod = ring_.Mul(ar[i], br[j]);
        const std::size_t idx = (i + j < n_) ? (i + j) : (i + j - n_);
        out[idx] = ring_.Add(out[idx], prod);
      }
    }
    return out;
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type acc = MultiplicativeIdentity();
    element_type base = Reduce(a);
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  element_type AdditiveIdentity() const {
    return element_type(n_, ring_.AdditiveIdentity());
  }

  element_type MultiplicativeIdentity() const {
    element_type one(n_, ring_.AdditiveIdentity());
    one[0] = ring_.MultiplicativeIdentity();
    return one;
  }

  bool IsZero(const element_type& a) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    for (std::size_t i = 0; i < n_; ++i) {
      if (!ring_.IsZero(ar[i])) {
        return false;
      }
    }
    return true;
  }

  bool IsOne(const element_type& a) const {
    return Equals(a, MultiplicativeIdentity());
  }

  bool Equals(const element_type& a, const element_type& b) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    const element_type br = (b.size() == n_) ? b : Reduce(b);
    for (std::size_t i = 0; i < n_; ++i) {
      if (!alg::Equals(ring_, ar[i], br[i])) {
        return false;
      }
    }
    return true;
  }

  std::string ToString(const element_type& a) const {
    const element_type ar = (a.size() == n_) ? a : Reduce(a);
    std::string out = "[";
    for (std::size_t i = 0; i < n_; ++i) {
      if (i > 0) {
        out += ", ";
      }
      out += ring_.ToString(ar[i]);
    }
    out += "]";
    return out;
  }

 private:
  R ring_;
  std::size_t n_;
};

}  // namespace ring

template <EuclideanDomain R>
struct AlgebraTraits<ring::QuotientRing<R>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <CommutativeRing R>
struct AlgebraTraits<ring::PolynomialQuotientRing<R>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <CommutativeRing R>
struct AlgebraTraits<ring::NegacyclicRing<R>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <CommutativeRing R>
struct AlgebraTraits<ring::CyclicRing<R>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_QUOTIENT_RING_HPP_
