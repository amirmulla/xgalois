#ifndef XGALOIS_ALGEBRA_RING_ARRAY_HPP_
#define XGALOIS_ALGEBRA_RING_ARRAY_HPP_

#include <concepts>
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/ring_element.hpp"
#include "xtensor/containers/xarray.hpp"

namespace xg {

// Forward declaration only — preserves Layer 1 independence from gf_array.hpp.
template <typename GaloisField>
class GFArray;

namespace algebra {

template <Ring R>
class RingArray {
 public:
  using ring_type = R;
  using element_type = typename R::element_type;
  using wrapped_element_type = RingElement<R>;
  using raw_array_type = xt::xarray<element_type>;
  using shape_type = typename raw_array_type::shape_type;

  RingArray() = default;

  static RingArray FromRaw(std::shared_ptr<const R> ring,
                           raw_array_type raw_data) {
    if (!ring) {
      throw std::invalid_argument("RingArray: ring pointer cannot be null");
    }
    RingArray arr;
    arr.ring_ = std::move(ring);
    arr.data_ = std::move(raw_data);
    return arr;
  }

  static RingArray Zeros(std::shared_ptr<const R> ring,
                         const std::vector<size_t>& shape) {
    if (!ring) {
      throw std::invalid_argument("RingArray: ring pointer cannot be null");
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type(shape.begin(), shape.end()));
    raw.fill(ring->AdditiveIdentity());
    return FromRaw(std::move(ring), std::move(raw));
  }

  static RingArray Ones(std::shared_ptr<const R> ring,
                        const std::vector<size_t>& shape) {
    if (!ring) {
      throw std::invalid_argument("RingArray: ring pointer cannot be null");
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type(shape.begin(), shape.end()));
    raw.fill(ring->MultiplicativeIdentity());
    return FromRaw(std::move(ring), std::move(raw));
  }

  static RingArray Eye(std::shared_ptr<const R> ring, size_t n) {
    if (!ring) {
      throw std::invalid_argument("RingArray: ring pointer cannot be null");
    }
    raw_array_type raw = raw_array_type::from_shape(shape_type{n, n});
    raw.fill(ring->AdditiveIdentity());
    const element_type one = ring->MultiplicativeIdentity();
    for (size_t i = 0; i < n; ++i) {
      raw(i, i) = one;
    }
    return FromRaw(std::move(ring), std::move(raw));
  }

  static RingArray FromVector(std::shared_ptr<const R> ring,
                              const std::vector<element_type>& values) {
    if (!ring) {
      throw std::invalid_argument("RingArray: ring pointer cannot be null");
    }
    raw_array_type raw = raw_array_type::from_shape(shape_type{values.size()});
    auto* dst = raw.data();
    for (size_t i = 0; i < values.size(); ++i) {
      dst[i] = values[i];
    }
    return FromRaw(std::move(ring), std::move(raw));
  }

  static RingArray FromMatrix(
      std::shared_ptr<const R> ring,
      const std::vector<std::vector<element_type>>& rows) {
    if (!ring) {
      throw std::invalid_argument("RingArray: ring pointer cannot be null");
    }
    if (rows.empty()) {
      return FromRaw(std::move(ring),
                     raw_array_type::from_shape(shape_type{0, 0}));
    }
    const size_t num_rows = rows.size();
    const size_t num_cols = rows[0].size();
    for (size_t r = 1; r < num_rows; ++r) {
      if (rows[r].size() != num_cols) {
        throw std::invalid_argument(
            "RingArray::FromMatrix: all rows must have the same length");
      }
    }
    raw_array_type raw =
        raw_array_type::from_shape(shape_type{num_rows, num_cols});
    for (size_t r = 0; r < num_rows; ++r) {
      for (size_t c = 0; c < num_cols; ++c) {
        raw(r, c) = rows[r][c];
      }
    }
    return FromRaw(std::move(ring), std::move(raw));
  }

  std::shared_ptr<const R> Ring() const noexcept { return ring_; }
  bool HasRing() const noexcept { return ring_ != nullptr; }

  const shape_type& Shape() const noexcept { return data_.shape(); }
  size_t Shape(size_t dim) const { return data_.shape(dim); }
  size_t Dimension() const noexcept { return data_.dimension(); }
  size_t Size() const noexcept { return data_.size(); }
  bool Empty() const noexcept { return data_.size() == 0; }

  raw_array_type& RawData() noexcept { return data_; }
  const raw_array_type& RawData() const noexcept { return data_; }

  element_type& operator()(size_t i) {
    RequireRing();
    return data_(i);
  }
  const element_type& operator()(size_t i) const {
    RequireRing();
    return data_(i);
  }

  element_type& operator()(size_t r, size_t c) {
    RequireRing();
    return data_(r, c);
  }
  const element_type& operator()(size_t r, size_t c) const {
    RequireRing();
    return data_(r, c);
  }

  wrapped_element_type At(size_t i) const {
    RequireRing();
    return wrapped_element_type::FromRaw(data_(i), ring_);
  }

  wrapped_element_type At(size_t r, size_t c) const {
    RequireRing();
    return wrapped_element_type::FromRaw(data_(r, c), ring_);
  }

  RingArray operator+(const RingArray& other) const {
    ValidateCompatible(other);
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    if constexpr (requires(const R& r, std::span<element_type> dst,
                           std::span<const element_type> src) {
                    r.AddArray(dst, src, src);
                  }) {
      ring_->AddArray(std::span<element_type>(out.data(), out.size()),
                      std::span<const element_type>(data_.data(), data_.size()),
                      std::span<const element_type>(other.data_.data(),
                                                    other.data_.size()));
    } else {
      const size_t n = data_.size();
      const auto* a = data_.data();
      const auto* b = other.data_.data();
      auto* dst = out.data();
      for (size_t i = 0; i < n; ++i) {
        dst[i] = ring_->Add(a[i], b[i]);
      }
    }
    return FromRaw(ring_, std::move(out));
  }

  RingArray operator-(const RingArray& other) const {
    ValidateCompatible(other);
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    if constexpr (requires(const R& r, std::span<element_type> dst,
                           std::span<const element_type> src) {
                    r.SubArray(dst, src, src);
                  }) {
      ring_->SubArray(std::span<element_type>(out.data(), out.size()),
                      std::span<const element_type>(data_.data(), data_.size()),
                      std::span<const element_type>(other.data_.data(),
                                                    other.data_.size()));
    } else {
      const size_t n = data_.size();
      const auto* a = data_.data();
      const auto* b = other.data_.data();
      auto* dst = out.data();
      for (size_t i = 0; i < n; ++i) {
        dst[i] = ring_->Sub(a[i], b[i]);
      }
    }
    return FromRaw(ring_, std::move(out));
  }

  RingArray operator-() const {
    RequireRing();
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    const size_t n = data_.size();
    const auto* src = data_.data();
    auto* dst = out.data();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = ring_->Neg(src[i]);
    }
    return FromRaw(ring_, std::move(out));
  }

  RingArray ElementwiseMul(const RingArray& other) const {
    ValidateCompatible(other);
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    if constexpr (requires(const R& r, std::span<element_type> dst,
                           std::span<const element_type> src) {
                    r.MulArray(dst, src, src);
                  }) {
      ring_->MulArray(std::span<element_type>(out.data(), out.size()),
                      std::span<const element_type>(data_.data(), data_.size()),
                      std::span<const element_type>(other.data_.data(),
                                                    other.data_.size()));
    } else {
      const size_t n = data_.size();
      const auto* a = data_.data();
      const auto* b = other.data_.data();
      auto* dst = out.data();
      for (size_t i = 0; i < n; ++i) {
        dst[i] = ring_->Mul(a[i], b[i]);
      }
    }
    return FromRaw(ring_, std::move(out));
  }

  // Left scalar multiplication: out[i] = scalar * arr[i] (respects
  // non-commutativity).
  RingArray LeftScalarMul(const element_type& scalar) const {
    RequireRing();
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    const size_t n = data_.size();
    const auto* src = data_.data();
    auto* dst = out.data();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = ring_->Mul(scalar, src[i]);
    }
    return FromRaw(ring_, std::move(out));
  }

  RingArray LeftScalarMul(const wrapped_element_type& scalar) const {
    RequireRing();
    if (scalar.Ring() != ring_) {
      throw std::invalid_argument("RingArray::LeftScalarMul: ring mismatch");
    }
    return LeftScalarMul(scalar.RawValue());
  }

  // Right scalar multiplication: out[i] = arr[i] * scalar (respects
  // non-commutativity).
  RingArray RightScalarMul(const element_type& scalar) const {
    RequireRing();
    raw_array_type out = raw_array_type::from_shape(data_.shape());
    const size_t n = data_.size();
    const auto* src = data_.data();
    auto* dst = out.data();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = ring_->Mul(src[i], scalar);
    }
    return FromRaw(ring_, std::move(out));
  }

  RingArray RightScalarMul(const wrapped_element_type& scalar) const {
    RequireRing();
    if (scalar.Ring() != ring_) {
      throw std::invalid_argument("RingArray::RightScalarMul: ring mismatch");
    }
    return RightScalarMul(scalar.RawValue());
  }

  RingArray Transpose() const {
    RequireRing();
    if (data_.dimension() != 2) {
      throw std::invalid_argument("RingArray::Transpose: requires a 2D matrix");
    }
    const size_t rows = data_.shape(0);
    const size_t cols = data_.shape(1);
    raw_array_type out = raw_array_type::from_shape(shape_type{cols, rows});
    for (size_t r = 0; r < rows; ++r) {
      for (size_t c = 0; c < cols; ++c) {
        out(c, r) = data_(r, c);
      }
    }
    return FromRaw(ring_, std::move(out));
  }

  bool operator==(const RingArray& other) const {
    if (!ring_ || !other.ring_) {
      return ring_ == other.ring_ && data_.shape() == other.data_.shape() &&
             data_.size() == 0;
    }
    if (ring_ != other.ring_) {
      return false;
    }
    if (data_.shape() != other.data_.shape()) {
      return false;
    }
    const size_t n = data_.size();
    const auto* a = data_.data();
    const auto* b = other.data_.data();
    for (size_t i = 0; i < n; ++i) {
      if (!alg::Equals(*ring_, a[i], b[i])) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const RingArray& other) const { return !(*this == other); }

 private:
  void RequireRing() const {
    if (!ring_) {
      throw std::runtime_error("RingArray: ring pointer is null");
    }
  }

  void ValidateCompatible(const RingArray& other) const {
    RequireRing();
    if (ring_ != other.ring_) {
      throw std::invalid_argument("RingArray: ring mismatch");
    }
    if (data_.shape() != other.data_.shape()) {
      throw std::invalid_argument("RingArray: shape mismatch");
    }
  }

  std::shared_ptr<const R> ring_{nullptr};
  raw_array_type data_;
};

template <Ring R>
using RingMatrix = RingArray<R>;

template <Ring R>
inline RingArray<R> ElementwiseMul(const RingArray<R>& a,
                                   const RingArray<R>& b) {
  return a.ElementwiseMul(b);
}

template <Ring R>
inline RingArray<R> LeftScalarMul(const typename R::element_type& scalar,
                                  const RingArray<R>& arr) {
  return arr.LeftScalarMul(scalar);
}

template <Ring R>
inline RingArray<R> LeftScalarMul(const RingElement<R>& scalar,
                                  const RingArray<R>& arr) {
  return arr.LeftScalarMul(scalar);
}

template <Ring R>
inline RingArray<R> RightScalarMul(const RingArray<R>& arr,
                                   const typename R::element_type& scalar) {
  return arr.RightScalarMul(scalar);
}

template <Ring R>
inline RingArray<R> RightScalarMul(const RingArray<R>& arr,
                                   const RingElement<R>& scalar) {
  return arr.RightScalarMul(scalar);
}

// Matrix multiplication C = A * B preserving operand order A(i, k) * B(k, j)
// for arbitrary (potentially non-commutative) rings.
template <Ring R>
inline RingArray<R> MatMul(const RingArray<R>& a, const RingArray<R>& b) {
  if (!a.HasRing() || !b.HasRing()) {
    throw std::invalid_argument("MatMul: ring pointer is null");
  }
  if (a.Ring() != b.Ring()) {
    throw std::invalid_argument("MatMul: ring mismatch");
  }
  if (a.Dimension() != 2 || b.Dimension() != 2) {
    throw std::invalid_argument("MatMul: operands must be 2D matrices");
  }
  const size_t m = a.Shape(0);
  const size_t k_dim = a.Shape(1);
  const size_t b_rows = b.Shape(0);
  const size_t n = b.Shape(1);
  if (k_dim != b_rows) {
    throw std::invalid_argument("MatMul: inner dimensions must match");
  }

  auto ring = a.Ring();
  using E = typename R::element_type;
  typename RingArray<R>::raw_array_type out =
      RingArray<R>::raw_array_type::from_shape(
          typename RingArray<R>::shape_type{m, n});
  const E zero = ring->AdditiveIdentity();
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      E acc = zero;
      for (size_t k = 0; k < k_dim; ++k) {
        acc = ring->Add(acc, ring->Mul(a(i, k), b(k, j)));
      }
      out(i, j) = std::move(acc);
    }
  }
  return RingArray<R>::FromRaw(std::move(ring), std::move(out));
}

// ============================================================================
// Duck-typed, zero-copy / move-friendly bridges with GFArray<F>
// ============================================================================

template <typename T>
concept GFArrayLike = requires(const T& arr) {
  typename T::field_type;
  typename T::element_type;
  typename T::raw_array_type;
  requires Field<typename T::field_type>;
  {
    arr.Field()
  } -> std::convertible_to<std::shared_ptr<const typename T::field_type>>;
  {
    arr.RawData()
  } -> std::convertible_to<const xt::xarray<typename T::element_type>&>;
};

template <GFArrayLike Arr>
inline RingArray<typename Arr::field_type> ToRingArray(const Arr& gf_arr) {
  using F = typename Arr::field_type;
  if (!gf_arr.Field()) {
    throw std::invalid_argument("ToRingArray: field pointer is null");
  }
  return RingArray<F>::FromRaw(gf_arr.Field(), gf_arr.RawData());
}

template <GFArrayLike Arr>
inline RingArray<typename Arr::field_type> ToRingArray(Arr&& gf_arr)
  requires(!std::is_lvalue_reference_v<Arr>)
{
  using F = typename Arr::field_type;
  if (!gf_arr.Field()) {
    throw std::invalid_argument("ToRingArray: field pointer is null");
  }
  return RingArray<F>::FromRaw(gf_arr.Field(), std::move(gf_arr.RawData()));
}

// Bridges RingArray<F> back to GFArray<F> when xgalois/field/gf_array.hpp is
// included by the caller. Note: the underlying field instance must not have
// been originally constructed as const if mutated via the returned GFArray.
template <Field F>
inline ::xg::GFArray<F> ToGFArray(const RingArray<F>& ring_arr) {
  if (!ring_arr.HasRing()) {
    throw std::invalid_argument("ToGFArray: ring pointer is null");
  }
  auto field_ptr = std::const_pointer_cast<F>(ring_arr.Ring());
  return ::xg::GFArray<F>::FromRaw(std::move(field_ptr), ring_arr.RawData());
}

template <Field F>
inline ::xg::GFArray<F> ToGFArray(RingArray<F>&& ring_arr) {
  if (!ring_arr.HasRing()) {
    throw std::invalid_argument("ToGFArray: ring pointer is null");
  }
  auto field_ptr = std::const_pointer_cast<F>(ring_arr.Ring());
  return ::xg::GFArray<F>::FromRaw(std::move(field_ptr),
                                   std::move(ring_arr.RawData()));
}

}  // namespace algebra

using algebra::ElementwiseMul;
using algebra::GFArrayLike;
using algebra::LeftScalarMul;
using algebra::MatMul;
using algebra::RightScalarMul;
using algebra::RingArray;
using algebra::RingMatrix;
using algebra::ToGFArray;
using algebra::ToRingArray;

}  // namespace xg

#endif  // XGALOIS_ALGEBRA_RING_ARRAY_HPP_
