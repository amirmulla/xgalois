#ifndef XGALOIS_FIELD_GF_ARRAY_HPP_
#define XGALOIS_FIELD_GF_ARRAY_HPP_

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/field/gf_element.hpp"
#include "xtensor/containers/xarray.hpp"

namespace xg {

// Contiguous Galois-field array holding a single shared_ptr<GaloisField>
// alongside a packed xt::xarray<storage_type> buffer. Eliminates 24x memory
// blowup and per-element atomic refcount churn while unlocking bulk SIMD
// kernels in linalg, coding, and channel modules.
template <typename GaloisField>
class GFArray {
 public:
  using field_type = GaloisField;
  using element_type = typename GaloisField::element_type;
  using wrapped_element_type = GaloisFieldElement<GaloisField>;
  using raw_array_type = xt::xarray<element_type>;
  using shape_type = typename raw_array_type::shape_type;

  GFArray() = default;

  GFArray(std::shared_ptr<GaloisField> field,
          const raw_array_type &external_values)
      : field_(std::move(field)), data_(external_values) {
    if (!field_) {
      throw std::invalid_argument("GFArray: field pointer cannot be null");
    }
    for (auto &v : data_) {
      v = field_->SetElementValue(v);
    }
  }

  static GFArray FromRaw(std::shared_ptr<GaloisField> field,
                         raw_array_type raw_data) {
    if (!field) {
      throw std::invalid_argument("GFArray: field pointer cannot be null");
    }
    GFArray arr;
    arr.field_ = std::move(field);
    arr.data_ = std::move(raw_data);
    return arr;
  }

  static GFArray Zeros(std::shared_ptr<GaloisField> field,
                       const std::vector<size_t> &shape) {
    if (!field) {
      throw std::invalid_argument("GFArray: field pointer cannot be null");
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type(shape.begin(), shape.end()));
    raw.fill(field->AdditiveIdentity());
    return FromRaw(std::move(field), std::move(raw));
  }

  static GFArray Ones(std::shared_ptr<GaloisField> field,
                      const std::vector<size_t> &shape) {
    if (!field) {
      throw std::invalid_argument("GFArray: field pointer cannot be null");
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type(shape.begin(), shape.end()));
    raw.fill(field->MultiplicativeIdentity());
    return FromRaw(std::move(field), std::move(raw));
  }

  static GFArray Random(std::shared_ptr<GaloisField> field,
                        const std::vector<size_t> &shape) {
    if (!field) {
      throw std::invalid_argument("GFArray: field pointer cannot be null");
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type(shape.begin(), shape.end()));
    for (auto &v : raw) {
      v = field->Random();
    }
    return FromRaw(std::move(field), std::move(raw));
  }

  static GFArray Random(std::shared_ptr<GaloisField> field,
                        const std::vector<size_t> &shape,
                        std::mt19937_64 &rng) {
    if (!field) {
      throw std::invalid_argument("GFArray: field pointer cannot be null");
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type(shape.begin(), shape.end()));
    for (auto &v : raw) {
      v = field->Random(rng);
    }
    return FromRaw(std::move(field), std::move(raw));
  }

  // Convert from legacy xt::xarray<GaloisFieldElement<GaloisField>>.
  static GFArray FromElementArray(
      const xt::xarray<wrapped_element_type> &elem_arr) {
    if (elem_arr.size() == 0) {
      throw std::invalid_argument(
          "GFArray::FromElementArray: cannot infer field from empty array");
    }
    auto field = elem_arr.data()[0].Field();
    if (!field) {
      throw std::invalid_argument("GFArray::FromElementArray: null field");
    }
    raw_array_type raw = raw_array_type::from_shape(elem_arr.shape());
    const auto *src = elem_arr.data();
    auto *dst = raw.data();
    const size_t n = elem_arr.size();
    for (size_t i = 0; i < n; ++i) {
      if (src[i].Field() != field) {
        throw std::invalid_argument(
            "GFArray::FromElementArray: mixed fields in element array");
      }
      dst[i] = src[i].RawValue();
    }
    return FromRaw(std::move(field), std::move(raw));
  }

  xt::xarray<wrapped_element_type> ToElementArray() const {
    RequireField();
    xt::xarray<wrapped_element_type> out =
        xt::xarray<wrapped_element_type>::from_shape(data_.shape());
    const auto *src = data_.data();
    auto *dst = out.data();
    const size_t n = data_.size();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = wrapped_element_type::FromInternal(src[i], field_);
    }
    return out;
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }
  bool HasField() const noexcept { return field_ != nullptr; }

  size_t Size() const noexcept { return data_.size(); }
  bool Empty() const noexcept { return data_.size() == 0; }
  size_t Dimension() const noexcept { return data_.dimension(); }
  const shape_type &Shape() const noexcept { return data_.shape(); }
  size_t Shape(size_t axis) const { return data_.shape(axis); }

  raw_array_type &RawData() noexcept { return data_; }
  const raw_array_type &RawData() const noexcept { return data_; }

  std::span<element_type> RawSpan() noexcept {
    return std::span<element_type>(data_.data(), data_.size());
  }
  std::span<const element_type> RawSpan() const noexcept {
    return std::span<const element_type>(data_.data(), data_.size());
  }

  element_type &Raw(size_t i) { return data_(i); }
  const element_type &Raw(size_t i) const { return data_(i); }
  element_type &Raw(size_t i, size_t j) { return data_(i, j); }
  const element_type &Raw(size_t i, size_t j) const { return data_(i, j); }

  wrapped_element_type At(size_t i) const {
    RequireField();
    return wrapped_element_type::FromInternal(data_(i), field_);
  }

  wrapped_element_type At(size_t i, size_t j) const {
    RequireField();
    return wrapped_element_type::FromInternal(data_(i, j), field_);
  }

  void Set(size_t i, const wrapped_element_type &elem) {
    RequireField();
    if (elem.Field() != field_) {
      throw std::invalid_argument("GFArray::Set: field mismatch");
    }
    data_(i) = elem.RawValue();
  }

  void Set(size_t i, size_t j, const wrapped_element_type &elem) {
    RequireField();
    if (elem.Field() != field_) {
      throw std::invalid_argument("GFArray::Set: field mismatch");
    }
    data_(i, j) = elem.RawValue();
  }

  void SetExternal(size_t i, const element_type &val) {
    RequireField();
    data_(i) = field_->SetElementValue(val);
  }

  void SetExternal(size_t i, size_t j, const element_type &val) {
    RequireField();
    data_(i, j) = field_->SetElementValue(val);
  }

  GFArray operator+(const GFArray &other) const {
    ValidateCompatible(other);
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    field_->AddArray(std::span<element_type>(out.data(), out.size()), RawSpan(),
                     other.RawSpan());
    return FromRaw(field_, std::move(out));
  }

  GFArray operator-(const GFArray &other) const {
    ValidateCompatible(other);
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    field_->SubArray(std::span<element_type>(out.data(), out.size()), RawSpan(),
                     other.RawSpan());
    return FromRaw(field_, std::move(out));
  }

  GFArray operator*(const GFArray &other) const {
    ValidateCompatible(other);
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    field_->MulArray(std::span<element_type>(out.data(), out.size()), RawSpan(),
                     other.RawSpan());
    return FromRaw(field_, std::move(out));
  }

  GFArray operator*(const wrapped_element_type &scalar) const {
    RequireField();
    if (scalar.Field() != field_) {
      throw std::invalid_argument("GFArray: scalar field mismatch");
    }
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    field_->ScalarMulArray(std::span<element_type>(out.data(), out.size()),
                           RawSpan(), scalar.RawValue());
    return FromRaw(field_, std::move(out));
  }

  GFArray operator-() const {
    RequireField();
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    const size_t n = data_.size();
    const auto *src = data_.data();
    auto *dst = out.data();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = field_->Neg(src[i]);
    }
    return FromRaw(field_, std::move(out));
  }

  GFArray &operator+=(const GFArray &other) {
    ValidateCompatible(other);
    field_->AddArray(RawSpan(), RawSpan(), other.RawSpan());
    return *this;
  }

  GFArray &operator-=(const GFArray &other) {
    ValidateCompatible(other);
    field_->SubArray(RawSpan(), RawSpan(), other.RawSpan());
    return *this;
  }

  bool operator==(const GFArray &other) const {
    if (field_ != other.field_) {
      return false;
    }
    return data_ == other.data_;
  }

  bool operator!=(const GFArray &other) const { return !(*this == other); }

 private:
  void RequireField() const {
    if (!field_) {
      throw std::runtime_error("GFArray: field pointer is null");
    }
  }

  void ValidateCompatible(const GFArray &other) const {
    RequireField();
    if (field_ != other.field_) {
      throw std::invalid_argument("GFArray: field mismatch");
    }
    if (data_.shape() != other.data_.shape()) {
      throw std::invalid_argument("GFArray: shape mismatch");
    }
  }

  std::shared_ptr<GaloisField> field_{nullptr};
  raw_array_type data_;
};

template <typename GaloisField>
using GFMatrix = GFArray<GaloisField>;

template <typename GaloisField>
inline GFArray<GaloisField> operator*(
    const GaloisFieldElement<GaloisField> &scalar,
    const GFArray<GaloisField> &arr) {
  return arr * scalar;
}

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_ARRAY_HPP_
