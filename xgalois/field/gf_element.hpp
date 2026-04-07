#ifndef XGALOIS_FIELD_GF_ELEMENT_HPP
#define XGALOIS_FIELD_GF_ELEMENT_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

namespace xg {

template <typename GaloisField>
class GaloisFieldElementBase {
 public:
  using ElementType = typename GaloisField::element_type;

  GaloisFieldElementBase() = default;

  GaloisFieldElementBase(ElementType value, std::shared_ptr<GaloisField> field)
      : field_(std::move(field)) {
    if (!field_) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    value_ = field_->SetElementValue(value);
  }

  GaloisFieldElementBase(const std::string &value_str,
                         std::shared_ptr<GaloisField> field)
      : field_(std::move(field)) {
    if (!field_) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    value_ = field_->SetElementValue(value_str);
  }

  GaloisFieldElementBase(const GaloisFieldElementBase &) = default;
  GaloisFieldElementBase &operator=(const GaloisFieldElementBase &) = default;
  GaloisFieldElementBase(GaloisFieldElementBase &&) noexcept = default;
  GaloisFieldElementBase &operator=(GaloisFieldElementBase &&) noexcept =
      default;

  GaloisFieldElementBase &operator=(const std::string &value_str) {
    if (!field_) {
      throw std::runtime_error("Field pointer is null");
    }
    value_ = field_->SetElementValue(value_str);
    return *this;
  }

  GaloisFieldElementBase &operator=(const ElementType &value) {
    if (!field_) {
      throw std::runtime_error("Field pointer is null");
    }
    value_ = field_->SetElementValue(value);
    return *this;
  }

  ElementType Value() const {
    if (!field_) {
      throw std::runtime_error("Field pointer is null");
    }
    return field_->GetElementValue(value_);
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }

  inline GaloisFieldElementBase operator+(
      const GaloisFieldElementBase &other) const {
    return {field_->Add(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase operator-(
      const GaloisFieldElementBase &other) const {
    return {field_->Sub(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase operator*(
      const GaloisFieldElementBase &other) const {
    return {field_->Mul(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase operator/(
      const GaloisFieldElementBase &other) const {
    return {field_->Div(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase operator-() const {
    return {field_->Neg(value_), field_};
  }

  inline GaloisFieldElementBase operator^(uint64_t exponent) const {
    return {field_->Pow(value_, exponent), field_};
  }

  inline GaloisFieldElementBase Inv() const {
    return {field_->Inv(value_), field_};
  }

  inline GaloisFieldElementBase Sqrt() const {
    return {field_->Sqrt(value_), field_};
  }

  inline GaloisFieldElementBase Pow(uint32_t exponent) const {
    return {field_->Pow(value_, exponent), field_};
  }

  inline GaloisFieldElementBase &operator+=(
      const GaloisFieldElementBase &other) {
    value_ = field_->Add(value_, other.value_);
    return *this;
  }

  inline GaloisFieldElementBase &operator-=(
      const GaloisFieldElementBase &other) {
    value_ = field_->Sub(value_, other.value_);
    return *this;
  }

  inline GaloisFieldElementBase &operator*=(
      const GaloisFieldElementBase &other) {
    value_ = field_->Mul(this->value_, other.value_);
    return *this;
  }

  inline GaloisFieldElementBase &operator/=(
      const GaloisFieldElementBase &other) {
    value_ = field_->Div(this->value_, other.value_);
    return *this;
  }

  inline bool operator==(const GaloisFieldElementBase &other) const {
    return (value_ == other.value_);
  }

  inline bool operator!=(const GaloisFieldElementBase &other) const {
    return !(*this == other);
  }

  inline bool operator<(const GaloisFieldElementBase &other) const {
    return value_ < other.value_;
  }

  inline bool operator>(const GaloisFieldElementBase &other) const {
    return value_ > other.value_;
  }

  inline bool operator<=(const GaloisFieldElementBase &other) const {
    return value_ <= other.value_;
  }

  inline bool operator>=(const GaloisFieldElementBase &other) const {
    return value_ >= other.value_;
  }

  inline void Print(std::ostream &os) const {
    if (field_) {
      field_->Print(value_, os);
    } else {
      throw std::runtime_error("Cannot print element: field pointer is null");
    }
  }

 protected:
  ElementType value_;
  std::shared_ptr<GaloisField> field_;
};

template <typename GaloisField>
inline std::ostream &operator<<(
    std::ostream &os,
    const GaloisFieldElementBase<GaloisField> &element) noexcept {
  element.Print(os);
  return os;
}

template <typename GaloisField>
inline std::istream &operator>>(std::istream &is,
                                GaloisFieldElementBase<GaloisField> &element) {
  typename GaloisField::element_type value;
  is >> value;
  if (element.Field()) {
    element = GaloisFieldElementBase<GaloisField>(value, element.Field());
  } else {
    is.setstate(std::ios::failbit);
  }
  return is;
}

template <typename GaloisField>
class GaloisFieldElement : public GaloisFieldElementBase<GaloisField> {
 public:
  using Base = GaloisFieldElementBase<GaloisField>;
  using ElementType = typename GaloisField::element_type;

  GaloisFieldElement() = default;

  GaloisFieldElement(ElementType value, std::shared_ptr<GaloisField> field)
      : Base(value, std::move(field)) {}

  GaloisFieldElement(const std::string &value_str,
                     std::shared_ptr<GaloisField> field)
      : Base(value_str, field) {}

  GaloisFieldElement(const GaloisFieldElement &) = default;
  GaloisFieldElement &operator=(const GaloisFieldElement &) = default;
  GaloisFieldElement(GaloisFieldElement &&) noexcept = default;
  GaloisFieldElement &operator=(GaloisFieldElement &&) noexcept = default;

  GaloisFieldElement &operator=(const std::string &value_str) {
    if (!this->field_) {
      throw std::runtime_error("Field pointer is null");
    }
    this->value_ = this->field_->SetElementValue(value_str);
    return *this;
  }

  GaloisFieldElement &operator=(const ElementType &value) {
    if (!this->field_) {
      throw std::runtime_error("Field pointer is null");
    }
    this->value_ = this->field_->SetElementValue(value);
    return *this;
  }

  explicit GaloisFieldElement(const Base &other)
      : Base(other.Value(), other.Field()) {}

 private:

  inline void ValidateField(const GaloisFieldElement &other) const {
    if (this->field_ != other.Field()) {
      throw std::invalid_argument(
          "Operations between elements of different fields are not allowed");
    }
  }

 public:

  inline GaloisFieldElement operator+(const GaloisFieldElement &other) const {
    ValidateField(other);
    return {this->field_->Add(this->value_, other.value_), this->field_};
  }

  inline GaloisFieldElement operator-(const GaloisFieldElement &other) const {
    ValidateField(other);
    return {this->field_->Sub(this->value_, other.value_), this->field_};
  }

  inline GaloisFieldElement operator*(const GaloisFieldElement &other) const {
    ValidateField(other);
    return {this->field_->Mul(this->value_, other.value_), this->field_};
  }

  inline GaloisFieldElement operator/(const GaloisFieldElement &other) const {
    ValidateField(other);
    return {this->field_->Div(this->value_, other.value_), this->field_};
  }

  inline GaloisFieldElement &operator+=(const GaloisFieldElement &other) {
    ValidateField(other);
    this->value_ = this->field_->Add(this->value_, other.value_);
    return *this;
  }

  inline GaloisFieldElement &operator-=(const GaloisFieldElement &other) {
    ValidateField(other);
    this->value_ = this->field_->Sub(this->value_, other.value_);
    return *this;
  }

  inline GaloisFieldElement &operator*=(const GaloisFieldElement &other) {
    ValidateField(other);
    this->value_ = this->field_->Mul(this->value_, other.value_);
    return *this;
  }

  inline GaloisFieldElement &operator/=(const GaloisFieldElement &other) {
    ValidateField(other);
    this->value_ = this->field_->Div(this->value_, other.value_);
    return *this;
  }

  inline GaloisFieldElement operator-() const {
    return {this->field_->Neg(this->value_), this->field_};
  }

  inline bool operator==(const GaloisFieldElement &other) const {
    ValidateField(other);
    return this->value_ == other.value_;
  }

  inline bool operator!=(const GaloisFieldElement &other) const {
    return !(*this == other);
  }

  inline bool operator==(const Base &other) const {
    return this->value_ == other.Value();
  }

  inline bool operator!=(const Base &other) const { return !(*this == other); }

  inline bool operator<(const GaloisFieldElement &other) const {
    ValidateField(other);
    return this->value_ < other.value_;
  }

  inline bool operator>(const GaloisFieldElement &other) const {
    ValidateField(other);
    return this->value_ > other.value_;
  }

  inline bool operator<=(const GaloisFieldElement &other) const {
    ValidateField(other);
    return this->value_ <= other.value_;
  }

  inline bool operator>=(const GaloisFieldElement &other) const {
    ValidateField(other);
    return this->value_ >= other.value_;
  }
};

template <typename GaloisField>
inline std::ostream &operator<<(
    std::ostream &os, const GaloisFieldElement<GaloisField> &element) noexcept {
  element.Print(os);
  return os;
}

template <typename GaloisField>
inline std::istream &operator>>(std::istream &is,
                                GaloisFieldElement<GaloisField> &element) {
  typename GaloisField::element_type value;
  is >> value;
  if (element.Field()) {
    element = GaloisFieldElement<GaloisField>(value, element.Field());
  } else {
    is.setstate(std::ios::failbit);
  }
  return is;
}

}

#endif
