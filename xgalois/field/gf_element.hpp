#ifndef XGALOIS_FIELD_GF_ELEMENT_HPP
#define XGALOIS_FIELD_GF_ELEMENT_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#include "xgalois/field/gf_base.hpp"

namespace xg {

//------------------------------------------------------------------------------
// GaloisFieldElementBase Class
//------------------------------------------------------------------------------

// Base class for elements in a Galois field, providing core arithmetic
// operations. This class manages the value and field relationship while
// implementing basic field arithmetic.
//
// Template Parameters:
//   GaloisField: The type representing the Galois field structure
template <typename GaloisField> class GaloisFieldElementBase {
public:
  using ElementType = typename GaloisField::element_type;

  // Default constructor
  GaloisFieldElementBase() = default;

  // Constructs a Galois field element with a value and its parent field.
  //
  // Args:
  //   value: The numerical value of the element
  //   field: Shared pointer to the parent Galois field
  // Throws:
  //   std::invalid_argument if field pointer is null
  GaloisFieldElementBase(ElementType value, std::shared_ptr<GaloisField> field)
      : field_(field) {
    if (!field_) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    value_ = field_->SetElementValue(value);
  }

  // Constructs a Galois field element from a string representation.
  //
  // Args:
  //   value_str: String representation of the element (e.g., "g^5", "α^2 + α + 1", "15")
  //   field: Shared pointer to the parent Galois field
  // Throws:
  //   std::invalid_argument if field pointer is null or string format is invalid
  GaloisFieldElementBase(const std::string &value_str, std::shared_ptr<GaloisField> field)
      : field_(field) {
    if (!field_) {
      throw std::invalid_argument("GaloisField pointer cannot be null");
    }
    value_ = field_->SetElementValue(value_str);
  }

  // Default copy/move constructors and assignment operators
  GaloisFieldElementBase(const GaloisFieldElementBase &) = default;
  GaloisFieldElementBase &operator=(const GaloisFieldElementBase &) = default;
  GaloisFieldElementBase(GaloisFieldElementBase &&) noexcept = default;
  GaloisFieldElementBase &
  operator=(GaloisFieldElementBase &&) noexcept = default;

  // Assignment operator for string values
  GaloisFieldElementBase &operator=(const std::string &value_str) {
    if (!field_) {
      throw std::runtime_error("Field pointer is null");
    }
    value_ = field_->SetElementValue(value_str);
    return *this;
  }

  // Assignment operator for numeric values
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

  inline GaloisFieldElementBase
  operator+(const GaloisFieldElementBase &other) const {
    return {field_->Add(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase
  operator-(const GaloisFieldElementBase &other) const {
    return {field_->Sub(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase
  operator*(const GaloisFieldElementBase &other) const {
    return {field_->Mul(value_, other.value_), field_};
  }

  inline GaloisFieldElementBase
  operator/(const GaloisFieldElementBase &other) const {
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

  inline GaloisFieldElementBase &
  operator+=(const GaloisFieldElementBase &other) {
    value_ = field_->Add(value_, other.value_);
    return *this;
  }

  inline GaloisFieldElementBase &
  operator-=(const GaloisFieldElementBase &other) {
    value_ = field_->Sub(value_, other.value_);
    return *this;
  }

  inline GaloisFieldElementBase &
  operator*=(const GaloisFieldElementBase &other) {
    value_ = field_->Mul(this->value_, other.value_);
    return *this;
  }

  inline GaloisFieldElementBase &
  operator/=(const GaloisFieldElementBase &other) {
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

  // Print the element using the field's Print method
  inline void Print(std::ostream &os) const {
    if (field_) {
      field_->Print(value_, os);
    } else {
      throw std::runtime_error("Cannot print element: field pointer is null");
    }
  }

protected:
  ElementType value_;                  // The underlying field element value
  std::shared_ptr<GaloisField> field_; // Pointer to the parent Galois field
};

// Stream operators for GaloisFieldElementBase
template <typename GaloisField>
inline std::ostream &
operator<<(std::ostream &os,
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

//------------------------------------------------------------------------------
// GaloisFieldElement Class
//------------------------------------------------------------------------------

// GaloisFieldElement extends the base class with additional validation to
// ensure arithmetic operations only occur between elements of the same field.
//
// Template Parameters:
//   GaloisField: The type of the Galois field this element belongs to.
template <typename GaloisField>
class GaloisFieldElement : public GaloisFieldElementBase<GaloisField> {
public:
  using Base = GaloisFieldElementBase<GaloisField>;
  using ElementType = typename GaloisField::element_type;

  // Default constructor
  GaloisFieldElement() = default;

  // Constructs a field element with the given value in the specified field.
  //
  // Args:
  //   value: The numerical value of the element
  //   field: Shared pointer to the parent Galois field
  GaloisFieldElement(ElementType value, std::shared_ptr<GaloisField> field)
      : Base(value, field) {}

  // Constructs a field element from a string representation.
  //
  // Args:
  //   value_str: String representation of the element (e.g., "g^5", "α^2 + α + 1", "15")
  //   field: Shared pointer to the parent Galois field
  GaloisFieldElement(const std::string &value_str, std::shared_ptr<GaloisField> field)
      : Base(value_str, field) {}

  // Default copy/move constructors and assignment operators
  GaloisFieldElement(const GaloisFieldElement &) = default;
  GaloisFieldElement &operator=(const GaloisFieldElement &) = default;
  GaloisFieldElement(GaloisFieldElement &&) noexcept = default;
  GaloisFieldElement &operator=(GaloisFieldElement &&) noexcept = default;

  // Assignment operator for string values
  GaloisFieldElement &operator=(const std::string &value_str) {
    if (!this->field_) {
      throw std::runtime_error("Field pointer is null");
    }
    this->value_ = this->field_->SetElementValue(value_str);
    return *this;
  }

  // Assignment operator for numeric values
  GaloisFieldElement &operator=(const ElementType &value) {
    if (!this->field_) {
      throw std::runtime_error("Field pointer is null");
    }
    this->value_ = this->field_->SetElementValue(value);
    return *this;
  }

  // Constructor to create a derived element from a base element.
  // Useful when calling base class methods that return Base.
  // Note: This assumes the Base object was actually created from this
  // specific derived type's field.
  GaloisFieldElement(const Base &other) : Base(other.Value(), other.Field()) {}

private:
  // Validates that arithmetic operations only occur between elements of the
  // same field.
  //
  // Args:
  //   other: The other field element to validate against
  // Throws:
  //   std::invalid_argument if the elements belong to different fields
  inline void ValidateField(const GaloisFieldElement &other) const {
    if (this->field_ != other.Field()) {
      throw std::invalid_argument(
          "Operations between elements of different fields are not allowed");
    }
  }

public:
  // Arithmetic operators that include field validation
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

// Stream operators for GaloisFieldElement
template <typename GaloisField>
inline std::ostream &
operator<<(std::ostream &os,
           const GaloisFieldElement<GaloisField> &element) noexcept {
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

} // namespace xg

#endif // XGALOIS_FIELD_GF_ELEMENT_HPP