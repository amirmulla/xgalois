#ifndef XGALOIS_ALGEBRA_CONCEPTS_HPP_
#define XGALOIS_ALGEBRA_CONCEPTS_HPP_

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace xg {

namespace detail {

// Detects XGalois's existing FiniteField contract without depending on
// xgalois/field/gf_concepts.hpp or modifying existing field classes.
template <typename T>
concept LegacyFiniteFieldLike =
    requires(const T& field, const typename T::element_type& a,
             const typename T::element_type& b) {
      typename T::element_type;
      { field.Characteristic() } -> std::convertible_to<uint64_t>;
      { field.Order() } -> std::convertible_to<uint64_t>;
      { field.Add(a, b) } -> std::same_as<typename T::element_type>;
      { field.Sub(a, b) } -> std::same_as<typename T::element_type>;
      { field.Mul(a, b) } -> std::same_as<typename T::element_type>;
      { field.Div(a, b) } -> std::same_as<typename T::element_type>;
      { field.Neg(a) } -> std::same_as<typename T::element_type>;
      { field.Inv(a) } -> std::same_as<typename T::element_type>;
      { field.AdditiveIdentity() } -> std::same_as<typename T::element_type>;
      {
        field.MultiplicativeIdentity()
      } -> std::same_as<typename T::element_type>;
    };

}  // namespace detail

// Algebraic property traits. Specialize for custom algebraic structures to
// opt into Semigroup, AbelianGroup, CommutativeRing, IntegralDomain,
// PrincipalIdealDomain, EuclideanDomain, or Field concepts.
// Types satisfying XGalois's existing FiniteField interface are automatically
// recognized as fields.
template <typename T>
struct AlgebraTraits {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = detail::LegacyFiniteFieldLike<T>;
  static constexpr bool kHasNoZeroDivisors = detail::LegacyFiniteFieldLike<T>;
  static constexpr bool kIsPid = detail::LegacyFiniteFieldLike<T>;
  static constexpr bool kIsEuclideanDomain = detail::LegacyFiniteFieldLike<T>;
  static constexpr bool kIsField = detail::LegacyFiniteFieldLike<T>;
};

namespace detail {

template <typename T>
constexpr bool TraitIsAssociativeOp() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (requires { AlgebraTraits<U>::kIsAssociativeOp; }) {
    return static_cast<bool>(AlgebraTraits<U>::kIsAssociativeOp);
  }
  return false;
}

template <typename T>
constexpr bool TraitIsCommutativeOp() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (requires { AlgebraTraits<U>::kIsCommutativeOp; }) {
    return static_cast<bool>(AlgebraTraits<U>::kIsCommutativeOp);
  }
  return false;
}

template <typename T>
constexpr bool TraitIsField() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (requires { AlgebraTraits<U>::kIsField; }) {
    return static_cast<bool>(AlgebraTraits<U>::kIsField);
  }
  return false;
}

template <typename T>
constexpr bool TraitIsEuclideanDomain() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (TraitIsField<U>()) {
    return true;
  }
  if constexpr (requires { AlgebraTraits<U>::kIsEuclideanDomain; }) {
    return static_cast<bool>(AlgebraTraits<U>::kIsEuclideanDomain);
  }
  return false;
}

template <typename T>
constexpr bool TraitIsPid() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (TraitIsEuclideanDomain<U>()) {
    return true;
  }
  if constexpr (requires { AlgebraTraits<U>::kIsPid; }) {
    return static_cast<bool>(AlgebraTraits<U>::kIsPid);
  }
  return false;
}

template <typename T>
constexpr bool TraitHasNoZeroDivisors() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (TraitIsPid<U>()) {
    return true;
  }
  if constexpr (requires { AlgebraTraits<U>::kHasNoZeroDivisors; }) {
    return static_cast<bool>(AlgebraTraits<U>::kHasNoZeroDivisors);
  }
  return false;
}

template <typename T>
constexpr bool TraitIsCommutativeMul() noexcept {
  using U = std::remove_cvref_t<T>;
  if constexpr (TraitIsPid<U>()) {
    return true;
  }
  if constexpr (requires { AlgebraTraits<U>::kIsCommutativeMul; }) {
    return static_cast<bool>(AlgebraTraits<U>::kIsCommutativeMul);
  }
  return false;
}

template <typename S>
concept HasStructureEquality =
    requires(const S& s, const typename S::element_type& a,
             const typename S::element_type& b) {
      { s.Equals(a, b) } -> std::convertible_to<bool>;
    } ||
    requires(const typename S::element_type& a,
             const typename S::element_type& b) {
      { a == b } -> std::convertible_to<bool>;
    };

}  // namespace detail

// ============================================================================
// Group Hierarchy: Magma -> Semigroup -> Monoid -> Group -> AbelianGroup
// ============================================================================

template <typename G>
concept Magma = requires(const G& g, const typename G::element_type& a,
                         const typename G::element_type& b) {
  typename G::element_type;
  { g.Op(a, b) } -> std::same_as<typename G::element_type>;
  { g.Equals(a, b) } -> std::convertible_to<bool>;
};

template <typename G>
concept Semigroup = Magma<G> && detail::TraitIsAssociativeOp<G>();

template <typename G>
concept Monoid =
    Semigroup<G> &&
    requires(const G& g, const typename G::element_type& a, uint64_t exp) {
      { g.Identity() } -> std::same_as<typename G::element_type>;
      { g.IsIdentity(a) } -> std::convertible_to<bool>;
      { g.Pow(a, exp) } -> std::same_as<typename G::element_type>;
    };

template <typename G>
concept Group =
    Monoid<G> && requires(const G& g, const typename G::element_type& a,
                          const typename G::element_type& b) {
      { g.Inv(a) } -> std::same_as<typename G::element_type>;
      { g.Div(a, b) } -> std::same_as<typename G::element_type>;
    };

template <typename G>
concept AbelianGroup = Group<G> && detail::TraitIsCommutativeOp<G>();

// ============================================================================
// Base Ring Hierarchy: Ring -> CommutativeRing -> IntegralDomain -> PID
// ============================================================================

template <typename R>
concept Ring = detail::HasStructureEquality<R> &&
               requires(const R& r, const typename R::element_type& a,
                        const typename R::element_type& b) {
                 typename R::element_type;
                 { r.Add(a, b) } -> std::same_as<typename R::element_type>;
                 { r.Sub(a, b) } -> std::same_as<typename R::element_type>;
                 { r.Neg(a) } -> std::same_as<typename R::element_type>;
                 { r.Mul(a, b) } -> std::same_as<typename R::element_type>;
                 {
                   r.AdditiveIdentity()
                 } -> std::same_as<typename R::element_type>;
                 {
                   r.MultiplicativeIdentity()
                 } -> std::same_as<typename R::element_type>;
                 { r.IsZero(a) } -> std::convertible_to<bool>;
                 { r.IsOne(a) } -> std::convertible_to<bool>;
                 { r.ToString(a) } -> std::convertible_to<std::string>;
               };

template <typename R>
concept CommutativeRing = Ring<R> && detail::TraitIsCommutativeMul<R>();

template <typename R>
concept IntegralDomain =
    CommutativeRing<R> && detail::TraitHasNoZeroDivisors<R>();

template <typename R>
concept PrincipalIdealDomain = IntegralDomain<R> && detail::TraitIsPid<R>();

namespace detail {

// Non-cyclic helper used by xg::alg CPOs so EuclideanDomain can depend on
// alg::DivRem / alg::EuclideanNorm and Field can subsume EuclideanDomain.
template <typename F>
concept FieldLike = PrincipalIdealDomain<F> && TraitIsField<F>() &&
                    requires(const F& field, const typename F::element_type& a,
                             const typename F::element_type& b) {
                      {
                        field.Inv(a)
                      } -> std::same_as<typename F::element_type>;
                      {
                        field.Div(a, b)
                      } -> std::same_as<typename F::element_type>;
                    };

}  // namespace detail

// ============================================================================
// Customization Point Objects (CPOs) in xg::alg
// ============================================================================

namespace alg {

namespace detail {

struct EqualsFn {
  template <typename S>
  constexpr bool operator()(const S& structure,
                            const typename S::element_type& a,
                            const typename S::element_type& b) const
    requires requires {
      { structure.Equals(a, b) } -> std::convertible_to<bool>;
    } || requires {
      { a == b } -> std::convertible_to<bool>;
    }
  {
    if constexpr (requires {
                    { structure.Equals(a, b) } -> std::convertible_to<bool>;
                  }) {
      return static_cast<bool>(structure.Equals(a, b));
    } else {
      return static_cast<bool>(a == b);
    }
  }
};

struct DivRemFn {
  template <typename R>
  constexpr std::pair<typename R::element_type, typename R::element_type>
  operator()(const R& ring, const typename R::element_type& a,
             const typename R::element_type& b) const
    requires requires {
      {
        ring.DivRem(a, b)
      } -> std::convertible_to<
            std::pair<typename R::element_type, typename R::element_type>>;
    } || ::xg::detail::FieldLike<R>
  {
    using E = typename R::element_type;
    if constexpr (requires {
                    {
                      ring.DivRem(a, b)
                    } -> std::convertible_to<std::pair<E, E>>;
                  }) {
      return std::pair<E, E>(ring.DivRem(a, b));
    } else {
      return std::pair<E, E>{ring.Div(a, b), ring.AdditiveIdentity()};
    }
  }
};

struct EuclideanNormFn {
  template <typename R>
  constexpr decltype(auto) operator()(const R& ring,
                                      const typename R::element_type& a) const
    requires requires { ring.EuclideanNorm(a); } || ::xg::detail::FieldLike<R>
  {
    if constexpr (requires { ring.EuclideanNorm(a); }) {
      return ring.EuclideanNorm(a);
    } else {
      return ring.IsZero(a) ? uint64_t{0} : uint64_t{1};
    }
  }
};

struct AssociateCanonicalFn {
  template <typename R>
  constexpr std::pair<typename R::element_type, typename R::element_type>
  operator()(const R& ring, const typename R::element_type& a) const
    requires requires {
      {
        ring.AssociateCanonical(a)
      } -> std::convertible_to<
            std::pair<typename R::element_type, typename R::element_type>>;
    } || ::xg::detail::FieldLike<R>
  {
    using E = typename R::element_type;
    if constexpr (requires {
                    {
                      ring.AssociateCanonical(a)
                    } -> std::convertible_to<std::pair<E, E>>;
                  }) {
      return std::pair<E, E>(ring.AssociateCanonical(a));
    } else {
      if (ring.IsZero(a)) {
        return std::pair<E, E>{ring.AdditiveIdentity(),
                               ring.MultiplicativeIdentity()};
      }
      return std::pair<E, E>{ring.MultiplicativeIdentity(), ring.Inv(a)};
    }
  }
};

struct IsUnitFn {
  template <typename R>
  constexpr bool operator()(const R& ring,
                            const typename R::element_type& a) const
    requires requires {
      { ring.IsUnit(a) } -> std::convertible_to<bool>;
    } ||
             requires {
               {
                 ring.TryInv(a)
               }
               -> std::convertible_to<std::optional<typename R::element_type>>;
             } || ::xg::detail::FieldLike<R>
  {
    if constexpr (requires {
                    { ring.IsUnit(a) } -> std::convertible_to<bool>;
                  }) {
      return static_cast<bool>(ring.IsUnit(a));
    } else if constexpr (requires {
                           {
                             ring.TryInv(a)
                           } -> std::convertible_to<
                                 std::optional<typename R::element_type>>;
                         }) {
      return ring.TryInv(a).has_value();
    } else {
      return !static_cast<bool>(ring.IsZero(a));
    }
  }
};

struct TryInvFn {
  template <typename R>
  constexpr std::optional<typename R::element_type> operator()(
      const R& ring, const typename R::element_type& a) const
    requires requires {
      {
        ring.TryInv(a)
      } -> std::convertible_to<std::optional<typename R::element_type>>;
    } || ::xg::detail::FieldLike<R>
  {
    using E = typename R::element_type;
    if constexpr (requires {
                    { ring.TryInv(a) } -> std::convertible_to<std::optional<E>>;
                  }) {
      return std::optional<E>(ring.TryInv(a));
    } else {
      if (ring.IsZero(a)) {
        return std::nullopt;
      }
      return std::optional<E>{ring.Inv(a)};
    }
  }
};

}  // namespace detail

inline constexpr detail::EqualsFn Equals{};
inline constexpr detail::DivRemFn DivRem{};
inline constexpr detail::EuclideanNormFn EuclideanNorm{};
inline constexpr detail::AssociateCanonicalFn AssociateCanonical{};
inline constexpr detail::IsUnitFn IsUnit{};
inline constexpr detail::TryInvFn TryInv{};

}  // namespace alg

// ============================================================================
// EuclideanDomain & Field
// ============================================================================

template <typename R>
concept EuclideanDomain =
    PrincipalIdealDomain<R> && detail::TraitIsEuclideanDomain<R>() &&
    requires(const R& r, const typename R::element_type& a,
             const typename R::element_type& b) {
      {
        alg::DivRem(r, a, b)
      } -> std::same_as<
            std::pair<typename R::element_type, typename R::element_type>>;
      { alg::EuclideanNorm(r, a) };
    };

template <typename F>
concept Field = EuclideanDomain<F> && detail::TraitIsField<F>() &&
                requires(const F& field, const typename F::element_type& a,
                         const typename F::element_type& b) {
                  { field.Inv(a) } -> std::same_as<typename F::element_type>;
                  { field.Div(a, b) } -> std::same_as<typename F::element_type>;
                };

// ============================================================================
// Module Hierarchy: Module<M, R>, VectorSpace<V, F>, Algebra<A, R>
// ============================================================================

template <typename M, typename R>
concept Module =
    Ring<R> && detail::HasStructureEquality<M> &&
    requires(const M& m, const typename M::element_type& x,
             const typename M::element_type& y,
             const typename R::element_type& s) {
      typename M::element_type;
      typename M::scalar_type;
      requires std::same_as<typename M::scalar_type, typename R::element_type>;
      { m.Add(x, y) } -> std::same_as<typename M::element_type>;
      { m.Sub(x, y) } -> std::same_as<typename M::element_type>;
      { m.Neg(x) } -> std::same_as<typename M::element_type>;
      { m.Scale(s, x) } -> std::same_as<typename M::element_type>;
      { m.AdditiveIdentity() } -> std::same_as<typename M::element_type>;
      { m.IsZero(x) } -> std::convertible_to<bool>;
    };

template <typename V, typename F>
concept VectorSpace = Field<F> && Module<V, F>;

template <typename A, typename R>
concept Algebra = CommutativeRing<R> && Ring<A> && Module<A, R>;

// ============================================================================
// Zero-Overhead Algebraic Bridge Views
// ============================================================================

namespace algebra {

// Adapts any Ring into its additive AbelianGroup (R, +, 0).
template <Ring R>
class AdditiveGroupView {
 public:
  using ring_type = R;
  using element_type = typename R::element_type;

  explicit constexpr AdditiveGroupView(const R& ring) noexcept : ring_(&ring) {}

  constexpr const R& Ring() const noexcept { return *ring_; }

  element_type Op(const element_type& a, const element_type& b) const {
    return ring_->Add(a, b);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return alg::Equals(*ring_, a, b);
  }

  element_type Identity() const { return ring_->AdditiveIdentity(); }

  bool IsIdentity(const element_type& a) const { return ring_->IsZero(a); }

  element_type Inv(const element_type& a) const { return ring_->Neg(a); }

  element_type Div(const element_type& a, const element_type& b) const {
    return ring_->Sub(a, b);
  }

  // Repeated doubling: computes n * a = a + a + ... + a (n times).
  element_type Pow(const element_type& a, uint64_t n) const {
    element_type acc = ring_->AdditiveIdentity();
    if (n == 0) {
      return acc;
    }
    element_type base = a;
    while (n > 0) {
      if ((n & uint64_t{1}) != 0) {
        acc = ring_->Add(acc, base);
      }
      n >>= 1;
      if (n > 0) {
        base = ring_->Add(base, base);
      }
    }
    return acc;
  }

 private:
  const R* ring_;
};

// Adapts any Ring into its multiplicative Monoid (R, *, 1).
template <Ring R>
class MultiplicativeMonoidView {
 public:
  using ring_type = R;
  using element_type = typename R::element_type;

  explicit constexpr MultiplicativeMonoidView(const R& ring) noexcept
      : ring_(&ring) {}

  constexpr const R& Ring() const noexcept { return *ring_; }

  element_type Op(const element_type& a, const element_type& b) const {
    return ring_->Mul(a, b);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return alg::Equals(*ring_, a, b);
  }

  element_type Identity() const { return ring_->MultiplicativeIdentity(); }

  bool IsIdentity(const element_type& a) const { return ring_->IsOne(a); }

  element_type Pow(const element_type& a, uint64_t n) const {
    if constexpr (requires(const R& r, const element_type& x, uint64_t e) {
                    { r.Pow(x, e) } -> std::same_as<element_type>;
                  }) {
      return ring_->Pow(a, n);
    } else {
      element_type acc = ring_->MultiplicativeIdentity();
      if (n == 0) {
        return acc;
      }
      element_type base = a;
      while (n > 0) {
        if ((n & uint64_t{1}) != 0) {
          acc = ring_->Mul(acc, base);
        }
        n >>= 1;
        if (n > 0) {
          base = ring_->Mul(base, base);
        }
      }
      return acc;
    }
  }

 private:
  const R* ring_;
};

// Adapts any Field into its multiplicative AbelianGroup (F*, *, 1) over
// nonzero elements.
template <Field F>
class MultiplicativeGroupView {
 public:
  using field_type = F;
  using element_type = typename F::element_type;

  explicit constexpr MultiplicativeGroupView(const F& field) noexcept
      : field_(&field) {}

  constexpr const F& Field() const noexcept { return *field_; }

  element_type Op(const element_type& a, const element_type& b) const {
    RequireNonZero(a);
    RequireNonZero(b);
    return field_->Mul(a, b);
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return alg::Equals(*field_, a, b);
  }

  element_type Identity() const { return field_->MultiplicativeIdentity(); }

  bool IsIdentity(const element_type& a) const { return field_->IsOne(a); }

  element_type Inv(const element_type& a) const {
    RequireNonZero(a);
    return field_->Inv(a);
  }

  element_type Div(const element_type& a, const element_type& b) const {
    RequireNonZero(a);
    RequireNonZero(b);
    return field_->Div(a, b);
  }

  element_type Pow(const element_type& a, uint64_t n) const {
    RequireNonZero(a);
    if constexpr (requires(const F& f, const element_type& x, uint64_t e) {
                    { f.Pow(x, e) } -> std::same_as<element_type>;
                  }) {
      return field_->Pow(a, n);
    } else {
      element_type acc = field_->MultiplicativeIdentity();
      if (n == 0) {
        return acc;
      }
      element_type base = a;
      while (n > 0) {
        if ((n & uint64_t{1}) != 0) {
          acc = field_->Mul(acc, base);
        }
        n >>= 1;
        if (n > 0) {
          base = field_->Mul(base, base);
        }
      }
      return acc;
    }
  }

 private:
  void RequireNonZero(const element_type& a) const {
    if (field_->IsZero(a)) {
      throw std::domain_error(
          "MultiplicativeGroupView: zero is not an element of F*");
    }
  }

  const F* field_;
};

}  // namespace algebra

template <Ring R>
struct AlgebraTraits<algebra::AdditiveGroupView<R>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = true;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <Ring R>
struct AlgebraTraits<algebra::MultiplicativeMonoidView<R>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = detail::TraitIsCommutativeMul<R>();
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

template <Field F>
struct AlgebraTraits<algebra::MultiplicativeGroupView<F>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = true;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

using algebra::AdditiveGroupView;
using algebra::MultiplicativeGroupView;
using algebra::MultiplicativeMonoidView;

}  // namespace xg

#endif  // XGALOIS_ALGEBRA_CONCEPTS_HPP_
