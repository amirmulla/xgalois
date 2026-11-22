#ifndef XGALOIS_ALGEBRA_GROUP_BASE_HPP_
#define XGALOIS_ALGEBRA_GROUP_BASE_HPP_

#include <concepts>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace algebra {

// Concept refining xg::Group for finite groups with exact BigInt order and
// an explicit generating set.
template <typename G>
concept FiniteGroup = xg::Group<G> && requires(const G& g) {
  { g.Order() } -> std::same_as<xg::BigInt>;
  { g.Generators() } -> std::same_as<std::vector<typename G::element_type>>;
};

namespace detail {

template <Group G>
inline std::string GroupElementToString(const G& group,
                                        const typename G::element_type& value) {
  using E = typename G::element_type;
  if constexpr (requires {
                  { group.ToString(value) } -> std::convertible_to<std::string>;
                }) {
    return std::string(group.ToString(value));
  } else if constexpr (requires {
                         {
                           value.ToString()
                         } -> std::convertible_to<std::string>;
                       }) {
    return std::string(value.ToString());
  } else if constexpr (std::is_arithmetic_v<E>) {
    return std::to_string(value);
  } else if constexpr (requires(std::ostream& os, const E& v) { os << v; }) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
  } else {
    return "<group_element>";
  }
}

}  // namespace detail

// Value-semantic wrapper holding a group element together with a shared
// handle to its ambient group.
template <Group G>
class GroupElement {
 public:
  using group_type = G;
  using element_type = typename G::element_type;
  using ElementType = element_type;

  GroupElement() = default;

  GroupElement(element_type value, std::shared_ptr<const G> group)
      : group_(std::move(group)), value_(std::move(value)) {
    if (!group_) {
      throw std::invalid_argument("GroupElement: group pointer cannot be null");
    }
  }

  GroupElement(std::shared_ptr<const G> group, element_type value)
    requires(!std::same_as<element_type, std::shared_ptr<const G>>)
      : group_(std::move(group)), value_(std::move(value)) {
    if (!group_) {
      throw std::invalid_argument("GroupElement: group pointer cannot be null");
    }
  }

  static GroupElement Identity(std::shared_ptr<const G> group) {
    if (!group) {
      throw std::invalid_argument("GroupElement: group pointer cannot be null");
    }
    element_type id = group->Identity();
    return GroupElement(std::move(id), std::move(group));
  }

  static GroupElement FromRaw(element_type raw_value,
                              std::shared_ptr<const G> group) {
    return GroupElement(std::move(raw_value), std::move(group));
  }

  GroupElement(const GroupElement&) = default;
  GroupElement& operator=(const GroupElement&) = default;
  GroupElement(GroupElement&&) noexcept = default;
  GroupElement& operator=(GroupElement&&) noexcept = default;

  std::shared_ptr<const G> Group() const noexcept { return group_; }

  bool HasGroup() const noexcept { return group_ != nullptr; }

  const element_type& RawValue() const {
    RequireGroup();
    return value_;
  }

  const element_type& Value() const { return RawValue(); }

  bool IsIdentity() const {
    RequireGroup();
    return static_cast<bool>(group_->IsIdentity(value_));
  }

  GroupElement Inv() const {
    RequireGroup();
    return GroupElement(group_->Inv(value_), group_);
  }

  GroupElement Pow(int64_t exponent) const {
    RequireGroup();
    if (exponent == 0) {
      return GroupElement(group_->Identity(), group_);
    }
    if (exponent > 0) {
      return GroupElement(group_->Pow(value_, static_cast<uint64_t>(exponent)),
                          group_);
    }
    const uint64_t abs_exp = uint64_t{0} - static_cast<uint64_t>(exponent);
    return GroupElement(group_->Pow(group_->Inv(value_), abs_exp), group_);
  }

  uint64_t Order() const {
    RequireGroup();
    if constexpr (requires(const G& g, const element_type& x) {
                    { g.ElementOrder(x) } -> std::convertible_to<uint64_t>;
                  }) {
      return static_cast<uint64_t>(group_->ElementOrder(value_));
    } else {
      if (group_->IsIdentity(value_)) {
        return 1;
      }
      element_type cur = value_;
      uint64_t ord = 1;
      while (!group_->IsIdentity(cur)) {
        cur = group_->Op(cur, value_);
        ++ord;
      }
      return ord;
    }
  }

  std::string ToString() const {
    RequireGroup();
    return detail::GroupElementToString(*group_, value_);
  }

  GroupElement operator*(const GroupElement& other) const {
    ValidateGroup(other);
    return GroupElement(group_->Op(value_, other.value_), group_);
  }

  GroupElement operator/(const GroupElement& other) const {
    ValidateGroup(other);
    return GroupElement(group_->Div(value_, other.value_), group_);
  }

  GroupElement& operator*=(const GroupElement& other) {
    ValidateGroup(other);
    value_ = group_->Op(value_, other.value_);
    return *this;
  }

  GroupElement& operator/=(const GroupElement& other) {
    ValidateGroup(other);
    value_ = group_->Div(value_, other.value_);
    return *this;
  }

  bool operator==(const GroupElement& other) const {
    if (!group_ || !other.group_) {
      return group_ == other.group_;
    }
    if (group_ != other.group_) {
      return false;
    }
    return static_cast<bool>(group_->Equals(value_, other.value_));
  }

  bool operator!=(const GroupElement& other) const { return !(*this == other); }

 private:
  void RequireGroup() const {
    if (!group_) {
      throw std::runtime_error("GroupElement: group pointer is null");
    }
  }

  void ValidateGroup(const GroupElement& other) const {
    if (!group_ || !other.group_) {
      throw std::invalid_argument("GroupElement: group pointer is null");
    }
    if (group_ != other.group_) {
      throw std::invalid_argument(
          "GroupElement: operations between elements of different groups are "
          "not allowed");
    }
  }

  std::shared_ptr<const G> group_{nullptr};
  element_type value_{};
};

template <Group G>
inline std::ostream& operator<<(std::ostream& os, const GroupElement<G>& elem) {
  os << elem.ToString();
  return os;
}

}  // namespace algebra

using algebra::FiniteGroup;
using algebra::GroupElement;

}  // namespace xg

namespace std {

template <xg::Group G>
struct hash<xg::algebra::GroupElement<G>> {
  size_t operator()(const xg::algebra::GroupElement<G>& elem) const noexcept {
    if (!elem.HasGroup()) {
      return 0;
    }
    const size_t h1 = std::hash<const G*>{}(elem.Group().get());
    size_t h2 = 0;
    if constexpr (std::is_integral_v<typename G::element_type>) {
      h2 = std::hash<typename G::element_type>{}(elem.RawValue());
    } else {
      h2 = std::hash<std::string>{}(elem.ToString());
    }
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

}  // namespace std

#endif  // XGALOIS_ALGEBRA_GROUP_BASE_HPP_
