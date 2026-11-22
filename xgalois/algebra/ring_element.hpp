#ifndef XGALOIS_ALGEBRA_RING_ELEMENT_HPP_
#define XGALOIS_ALGEBRA_RING_ELEMENT_HPP_

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "xgalois/algebra/concepts.hpp"

namespace xg {
namespace algebra {

template <Ring R>
class RingElement {
 public:
  using ring_type = R;
  using element_type = typename R::element_type;
  using ElementType = element_type;

  RingElement() = default;

  RingElement(element_type value, std::shared_ptr<const R> ring)
      : ring_(std::move(ring)), value_(std::move(value)) {
    if (!ring_) {
      throw std::invalid_argument("RingElement: ring pointer cannot be null");
    }
  }

  static RingElement Zero(std::shared_ptr<const R> ring) {
    if (!ring) {
      throw std::invalid_argument("RingElement: ring pointer cannot be null");
    }
    element_type z = ring->AdditiveIdentity();
    return RingElement(std::move(z), std::move(ring));
  }

  static RingElement One(std::shared_ptr<const R> ring) {
    if (!ring) {
      throw std::invalid_argument("RingElement: ring pointer cannot be null");
    }
    element_type o = ring->MultiplicativeIdentity();
    return RingElement(std::move(o), std::move(ring));
  }

  static RingElement FromRaw(element_type raw_value,
                             std::shared_ptr<const R> ring) {
    return RingElement(std::move(raw_value), std::move(ring));
  }

  RingElement(const RingElement&) = default;
  RingElement& operator=(const RingElement&) = default;
  RingElement(RingElement&&) noexcept = default;
  RingElement& operator=(RingElement&&) noexcept = default;

  std::shared_ptr<const R> Ring() const noexcept { return ring_; }

  bool HasRing() const noexcept { return ring_ != nullptr; }

  const element_type& RawValue() const {
    RequireRing();
    return value_;
  }

  bool IsZero() const {
    RequireRing();
    return ring_->IsZero(value_);
  }

  bool IsOne() const {
    RequireRing();
    return ring_->IsOne(value_);
  }

  bool IsUnit() const
    requires requires(const R& r, const element_type& a) { alg::IsUnit(r, a); }
  {
    RequireRing();
    return alg::IsUnit(*ring_, value_);
  }

  std::optional<RingElement> TryInv() const
    requires requires(const R& r, const element_type& a) { alg::TryInv(r, a); }
  {
    RequireRing();
    auto inv = alg::TryInv(*ring_, value_);
    if (!inv.has_value()) {
      return std::nullopt;
    }
    return RingElement(std::move(*inv), ring_);
  }

  RingElement Pow(uint64_t exponent) const {
    RequireRing();
    if constexpr (requires(const R& r, const element_type& a, uint64_t e) {
                    { r.Pow(a, e) } -> std::same_as<element_type>;
                  }) {
      return RingElement(ring_->Pow(value_, exponent), ring_);
    } else {
      element_type acc = ring_->MultiplicativeIdentity();
      if (exponent == 0) {
        return RingElement(std::move(acc), ring_);
      }
      element_type base = value_;
      while (exponent > 0) {
        if ((exponent & uint64_t{1}) != 0) {
          acc = ring_->Mul(acc, base);
        }
        exponent >>= 1;
        if (exponent > 0) {
          base = ring_->Mul(base, base);
        }
      }
      return RingElement(std::move(acc), ring_);
    }
  }

  std::string ToString() const {
    RequireRing();
    return ring_->ToString(value_);
  }

  RingElement operator+(const RingElement& other) const {
    ValidateRing(other);
    return RingElement(ring_->Add(value_, other.value_), ring_);
  }

  RingElement operator-(const RingElement& other) const {
    ValidateRing(other);
    return RingElement(ring_->Sub(value_, other.value_), ring_);
  }

  RingElement operator*(const RingElement& other) const {
    ValidateRing(other);
    return RingElement(ring_->Mul(value_, other.value_), ring_);
  }

  RingElement operator-() const {
    RequireRing();
    return RingElement(ring_->Neg(value_), ring_);
  }

  RingElement operator/(const RingElement& other) const
    requires Field<R>
  {
    ValidateRing(other);
    return RingElement(ring_->Div(value_, other.value_), ring_);
  }

  RingElement& operator+=(const RingElement& other) {
    ValidateRing(other);
    value_ = ring_->Add(value_, other.value_);
    return *this;
  }

  RingElement& operator-=(const RingElement& other) {
    ValidateRing(other);
    value_ = ring_->Sub(value_, other.value_);
    return *this;
  }

  RingElement& operator*=(const RingElement& other) {
    ValidateRing(other);
    value_ = ring_->Mul(value_, other.value_);
    return *this;
  }

  RingElement& operator/=(const RingElement& other)
    requires Field<R>
  {
    ValidateRing(other);
    value_ = ring_->Div(value_, other.value_);
    return *this;
  }

  bool operator==(const RingElement& other) const {
    if (!ring_ || !other.ring_) {
      return ring_ == other.ring_;
    }
    if (ring_ != other.ring_) {
      return false;
    }
    return alg::Equals(*ring_, value_, other.value_);
  }

  bool operator!=(const RingElement& other) const { return !(*this == other); }

 private:
  void RequireRing() const {
    if (!ring_) {
      throw std::runtime_error("RingElement: ring pointer is null");
    }
  }

  void ValidateRing(const RingElement& other) const {
    if (!ring_ || !other.ring_) {
      throw std::invalid_argument("RingElement: ring pointer is null");
    }
    if (ring_ != other.ring_) {
      throw std::invalid_argument(
          "RingElement: operations between elements of different rings are not "
          "allowed");
    }
  }

  std::shared_ptr<const R> ring_{nullptr};
  element_type value_{};
};

template <Ring R>
inline std::ostream& operator<<(std::ostream& os, const RingElement<R>& elem) {
  os << elem.ToString();
  return os;
}

}  // namespace algebra

using algebra::RingElement;

}  // namespace xg

namespace std {

template <xg::Ring R>
struct hash<xg::algebra::RingElement<R>> {
  size_t operator()(const xg::algebra::RingElement<R>& elem) const noexcept {
    if (!elem.HasRing()) {
      return 0;
    }
    const size_t h1 = std::hash<const R*>{}(elem.Ring().get());
    size_t h2 = 0;
    if constexpr (std::is_integral_v<typename R::element_type>) {
      h2 = std::hash<typename R::element_type>{}(elem.RawValue());
    } else {
      h2 = std::hash<std::string>{}(elem.ToString());
    }
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

}  // namespace std

#endif  // XGALOIS_ALGEBRA_RING_ELEMENT_HPP_
