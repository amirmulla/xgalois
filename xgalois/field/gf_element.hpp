#ifndef XGALOIS_FIELD_GF_ELEMENT_HPP_
#define XGALOIS_FIELD_GF_ELEMENT_HPP_

#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace xg {

template <typename GaloisField>
class GaloisFieldElement {
 private:
  struct InternalTag {};

  GaloisFieldElement(InternalTag, typename GaloisField::element_type raw_value,
                     std::shared_ptr<GaloisField> field)
      : value_(std::move(raw_value)), field_(std::move(field)) {}

  void RequireField() const {
    if (!field_) {
      throw std::runtime_error("Field pointer is null");
    }
  }

  void ValidateField(const GaloisFieldElement &other) const {
    if (!field_ || !other.field_) {
      throw std::invalid_argument("Field pointer is null");
    }
    if (field_ != other.field_) {
      if (field_->Order() != other.field_->Order() ||
          field_->Characteristic() != other.field_->Characteristic() ||
          field_->Degree() != other.field_->Degree()) {
        throw std::invalid_argument(
            "Operations between elements of different fields are not allowed");
      }
      // Distinct instances of the same field type with matching order may
      // still have different irreducible polynomials; reject distinct pointers
      // to guarantee mathematical soundness.
      throw std::invalid_argument(
          "Operations between elements of different fields are not allowed");
    }
  }

 public:
  using ElementType = typename GaloisField::element_type;
  using field_type = GaloisField;

  GaloisFieldElement() = default;

  GaloisFieldElement(ElementType value, std::shared_ptr<GaloisField> field)
      : field_(std::move(field)) {
    if (!field_) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    value_ = field_->SetElementValue(value);
  }

  GaloisFieldElement(const std::string &value_str,
                     std::shared_ptr<GaloisField> field)
      : field_(std::move(field)) {
    if (!field_) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    value_ = field_->SetElementValue(value_str);
  }

  // Construct directly from a field's internal representation without calling
  // SetElementValue (critical for Zech-log / Montgomery representations and
  // avoiding redundant modular reductions).
  static GaloisFieldElement FromInternal(ElementType raw_value,
                                         std::shared_ptr<GaloisField> field) {
    if (!field) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    return GaloisFieldElement(InternalTag{}, std::move(raw_value),
                              std::move(field));
  }

  GaloisFieldElement(const GaloisFieldElement &) = default;
  GaloisFieldElement &operator=(const GaloisFieldElement &) = default;
  GaloisFieldElement(GaloisFieldElement &&) noexcept = default;
  GaloisFieldElement &operator=(GaloisFieldElement &&) noexcept = default;

  GaloisFieldElement &operator=(const std::string &value_str) {
    RequireField();
    value_ = field_->SetElementValue(value_str);
    return *this;
  }

  GaloisFieldElement &operator=(const ElementType &value) {
    RequireField();
    value_ = field_->SetElementValue(value);
    return *this;
  }

  ElementType Value() const {
    RequireField();
    return field_->GetElementValue(value_);
  }

  const ElementType &RawValue() const {
    RequireField();
    return value_;
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }

  bool HasField() const noexcept { return field_ != nullptr; }

  bool IsZero() const {
    RequireField();
    return field_->IsZero(value_);
  }

  bool IsOne() const {
    RequireField();
    return field_->IsOne(value_);
  }

  GaloisFieldElement operator+(const GaloisFieldElement &other) const {
    ValidateField(other);
    return GaloisFieldElement(InternalTag{}, field_->Add(value_, other.value_),
                              field_);
  }

  GaloisFieldElement operator-(const GaloisFieldElement &other) const {
    ValidateField(other);
    return GaloisFieldElement(InternalTag{}, field_->Sub(value_, other.value_),
                              field_);
  }

  GaloisFieldElement operator*(const GaloisFieldElement &other) const {
    ValidateField(other);
    return GaloisFieldElement(InternalTag{}, field_->Mul(value_, other.value_),
                              field_);
  }

  GaloisFieldElement operator/(const GaloisFieldElement &other) const {
    ValidateField(other);
    return GaloisFieldElement(InternalTag{}, field_->Div(value_, other.value_),
                              field_);
  }

  GaloisFieldElement operator-() const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Neg(value_), field_);
  }

  GaloisFieldElement operator^(uint64_t exponent) const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Pow(value_, exponent),
                              field_);
  }

  GaloisFieldElement Inv() const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Inv(value_), field_);
  }

  GaloisFieldElement Sqrt() const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Sqrt(value_), field_);
  }

  GaloisFieldElement Pow(uint64_t exponent) const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Pow(value_, exponent),
                              field_);
  }

  GaloisFieldElement PowSigned(int64_t exponent) const {
    RequireField();
    return GaloisFieldElement(InternalTag{},
                              field_->PowSigned(value_, exponent), field_);
  }

  GaloisFieldElement Frobenius() const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Frobenius(value_), field_);
  }

  GaloisFieldElement Trace() const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Trace(value_), field_);
  }

  GaloisFieldElement Norm() const {
    RequireField();
    return GaloisFieldElement(InternalTag{}, field_->Norm(value_), field_);
  }

  uint64_t Log() const {
    RequireField();
    return field_->Log(value_);
  }

  uint64_t MultiplicativeOrder() const {
    RequireField();
    return field_->MultiplicativeOrder(value_);
  }

  GaloisFieldElement &operator+=(const GaloisFieldElement &other) {
    ValidateField(other);
    value_ = field_->Add(value_, other.value_);
    return *this;
  }

  GaloisFieldElement &operator-=(const GaloisFieldElement &other) {
    ValidateField(other);
    value_ = field_->Sub(value_, other.value_);
    return *this;
  }

  GaloisFieldElement &operator*=(const GaloisFieldElement &other) {
    ValidateField(other);
    value_ = field_->Mul(value_, other.value_);
    return *this;
  }

  GaloisFieldElement &operator/=(const GaloisFieldElement &other) {
    ValidateField(other);
    value_ = field_->Div(value_, other.value_);
    return *this;
  }

  bool operator==(const GaloisFieldElement &other) const {
    if (!field_ || !other.field_) {
      return field_ == other.field_;
    }
    if (field_ != other.field_) {
      return false;
    }
    return value_ == other.value_;
  }

  bool operator!=(const GaloisFieldElement &other) const {
    return !(*this == other);
  }

  bool operator<(const GaloisFieldElement &other) const {
    ValidateField(other);
    return Value() < other.Value();
  }

  bool operator<=(const GaloisFieldElement &other) const {
    ValidateField(other);
    return Value() <= other.Value();
  }

  bool operator>(const GaloisFieldElement &other) const {
    ValidateField(other);
    return Value() > other.Value();
  }

  bool operator>=(const GaloisFieldElement &other) const {
    ValidateField(other);
    return Value() >= other.Value();
  }

  void Print(std::ostream &os) const {
    RequireField();
    field_->Print(value_, os);
  }

  std::string ToString() const {
    RequireField();
    return field_->ToString(value_);
  }

 private:
  ElementType value_{};
  std::shared_ptr<GaloisField> field_{nullptr};
};

// Backwards-compatible unified alias eliminating the Base vs Strict split.
template <typename GaloisField>
using GaloisFieldElementBase = GaloisFieldElement<GaloisField>;

template <typename GaloisField>
inline std::ostream &operator<<(std::ostream &os,
                                const GaloisFieldElement<GaloisField> &elem) {
  elem.Print(os);
  return os;
}

template <typename GaloisField>
inline std::istream &operator>>(std::istream &is,
                                GaloisFieldElement<GaloisField> &elem) {
  using ElementType = typename GaloisField::element_type;
  if constexpr (std::is_same_v<ElementType, uint8_t>) {
    uint64_t tmp = 0;
    if (is >> tmp) {
      elem = static_cast<uint8_t>(tmp % elem.Field()->Order());
    }
  } else {
    ElementType val{};
    if (is >> val) {
      elem = val;
    }
  }
  return is;
}

}  // namespace xg

namespace std {

template <typename GaloisField>
struct hash<xg::GaloisFieldElement<GaloisField>> {
  size_t operator()(
      const xg::GaloisFieldElement<GaloisField> &elem) const noexcept {
    if (!elem.HasField()) {
      return 0;
    }
    const size_t h1 = std::hash<const GaloisField *>{}(elem.Field().get());
    size_t h2 = 0;
    if constexpr (std::is_integral_v<typename GaloisField::element_type>) {
      h2 = std::hash<typename GaloisField::element_type>{}(elem.RawValue());
    } else {
      h2 = std::hash<std::string>{}(elem.ToString());
    }
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

}  // namespace std

#endif  // XGALOIS_FIELD_GF_ELEMENT_HPP_
