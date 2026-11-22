#ifndef XGALOIS_RING_INTEGER_RING_HPP_
#define XGALOIS_RING_INTEGER_RING_HPP_

#include <concepts>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace ring {

// Ring of rational integers Z over IntT (defaulting to arbitrary-precision
// BigInt). Satisfies xg::EuclideanDomain (and PID, IntegralDomain,
// CommutativeRing, Ring).
template <typename IntT = BigInt>
class IntegerRing {
 public:
  using element_type = IntT;

  constexpr IntegerRing() noexcept = default;

  element_type Add(const element_type& a, const element_type& b) const {
    return a + b;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    return a - b;
  }

  element_type Neg(const element_type& a) const { return -a; }

  element_type Mul(const element_type& a, const element_type& b) const {
    return a * b;
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return a.Pow(exp);
    } else {
      element_type acc{1};
      element_type base = a;
      while (exp > 0) {
        if ((exp & uint64_t{1}) != 0) {
          acc = acc * base;
        }
        exp >>= 1;
        if (exp > 0) {
          base = base * base;
        }
      }
      return acc;
    }
  }

  element_type AdditiveIdentity() const { return element_type{0}; }

  element_type MultiplicativeIdentity() const { return element_type{1}; }

  bool IsZero(const element_type& a) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return a.IsZero();
    } else {
      return a == element_type{0};
    }
  }

  bool IsOne(const element_type& a) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return a.IsOne();
    } else {
      return a == element_type{1};
    }
  }

  bool Equals(const element_type& a, const element_type& b) const {
    return a == b;
  }

  std::string ToString(const element_type& a) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return a.ToString();
    } else {
      return std::to_string(a);
    }
  }

  // Euclidean division: returns {q, r} such that a == q * b + r and
  // 0 <= r < |b|.
  std::pair<element_type, element_type> DivRem(const element_type& a,
                                               const element_type& b) const {
    if (IsZero(b)) {
      throw std::domain_error("IntegerRing::DivRem: division by zero");
    }
    if constexpr (std::same_as<IntT, BigInt>) {
      return BigInt::EuclideanDivRem(a, b);
    } else {
      element_type q = a / b;
      element_type r = a % b;
      if (r < element_type{0}) {
        if (b > element_type{0}) {
          q = q - element_type{1};
          r = r + b;
        } else {
          q = q + element_type{1};
          r = r - b;
        }
      }
      return {std::move(q), std::move(r)};
    }
  }

  // Euclidean norm |a|.
  auto EuclideanNorm(const element_type& a) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return a.Abs();
    } else if constexpr (std::signed_integral<IntT>) {
      using UInt = std::make_unsigned_t<IntT>;
      return a < IntT{0} ? static_cast<UInt>(0) - static_cast<UInt>(a)
                         : static_cast<UInt>(a);
    } else {
      return a < element_type{0} ? -a : a;
    }
  }

  bool IsUnit(const element_type& a) const {
    return a == element_type{1} || a == element_type{-1};
  }

  std::optional<element_type> TryInv(const element_type& a) const {
    if (IsUnit(a)) {
      return a;
    }
    return std::nullopt;
  }

  // Returns {canonical, unit} where canonical = |a| = unit * a and
  // unit = sign(a) in {+1, -1} (with unit = +1 when a == 0).
  std::pair<element_type, element_type> AssociateCanonical(
      const element_type& a) const {
    if (IsZero(a)) {
      return {element_type{0}, element_type{1}};
    }
    if (a < element_type{0}) {
      return {-a, element_type{-1}};
    }
    return {a, element_type{1}};
  }

  // Extended Euclidean Algorithm: returns {g, s, t} with g >= 0 and
  // s * a + t * b == g.
  std::tuple<element_type, element_type, element_type> ExtendedGcd(
      const element_type& a, const element_type& b) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return ::xg::ExtendedGcd(a, b);
    } else {
      element_type old_r = a;
      element_type r = b;
      element_type old_s{1};
      element_type s{0};
      element_type old_t{0};
      element_type t{1};

      while (!IsZero(r)) {
        auto [q, rem] = DivRem(old_r, r);
        old_r = std::move(r);
        r = std::move(rem);

        element_type next_s = old_s - q * s;
        old_s = std::move(s);
        s = std::move(next_s);

        element_type next_t = old_t - q * t;
        old_t = std::move(t);
        t = std::move(next_t);
      }

      if (old_r < element_type{0}) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
      }
      return {std::move(old_r), std::move(old_s), std::move(old_t)};
    }
  }

  element_type Gcd(const element_type& a, const element_type& b) const {
    if constexpr (std::same_as<IntT, BigInt>) {
      return ::xg::Gcd(a, b);
    } else {
      return std::get<0>(ExtendedGcd(a, b));
    }
  }

  element_type Lcm(const element_type& a, const element_type& b) const {
    if (IsZero(a) || IsZero(b)) {
      return element_type{0};
    }
    const element_type g = Gcd(a, b);
    const element_type abs_a = AssociateCanonical(a).first;
    const element_type abs_b = AssociateCanonical(b).first;
    return DivRem(abs_a, g).first * abs_b;
  }
};

}  // namespace ring

template <typename IntT>
struct AlgebraTraits<ring::IntegerRing<IntT>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_INTEGER_RING_HPP_
