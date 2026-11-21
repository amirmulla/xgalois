#ifndef XGALOIS_FIELD_GF_BASE_HPP_
#define XGALOIS_FIELD_GF_BASE_HPP_

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "xgalois/utils/math.hpp"

namespace xg {

enum class FieldRepresentation : std::uint8_t { INT, HEX, POW, LOG, POLY };

template <typename ElementType>
class GaloisFieldBase {
 public:
  using element_type = ElementType;
  using storage_type = ElementType;

  virtual ~GaloisFieldBase() = default;

  virtual uint64_t Characteristic() const = 0;

  // Returns q = p^m for fields with q < 2^64. For GF(2^64), returns 0
  // (callers can inspect Degree() == 64 and Characteristic() == 2).
  virtual uint64_t Order() const = 0;

  virtual uint32_t Degree() const { return 1; }

  virtual ElementType Add(const ElementType &a, const ElementType &b) const = 0;
  virtual ElementType Sub(const ElementType &a, const ElementType &b) const = 0;
  virtual ElementType Mul(const ElementType &a, const ElementType &b) const = 0;
  virtual ElementType Div(const ElementType &a, const ElementType &b) const = 0;
  virtual ElementType Neg(const ElementType &a) const = 0;
  virtual ElementType Inv(const ElementType &a) const = 0;
  virtual ElementType Pow(const ElementType &a, uint64_t exp) const = 0;
  virtual ElementType Sqrt(const ElementType &a) const = 0;
  virtual uint64_t Log(const ElementType &a) const = 0;
  virtual uint64_t Log(const ElementType &a,
                       const ElementType &generator) const = 0;

  // Supports signed (including negative) exponents uniformly across all fields.
  ElementType PowSigned(const ElementType &a, int64_t exp) const {
    if (exp >= 0) {
      return Pow(a, static_cast<uint64_t>(exp));
    }
    ElementType inv_a = Inv(a);
    // Avoid signed overflow on INT64_MIN.
    uint64_t mag = static_cast<uint64_t>(-(exp + 1)) + 1ULL;
    return Pow(inv_a, mag);
  }

  virtual ElementType Random() const = 0;

  virtual ElementType Random(std::mt19937_64 & /*rng*/) const {
    return Random();
  }

  virtual ElementType MultiplicativeIdentity() const = 0;
  virtual ElementType AdditiveIdentity() const = 0;

  virtual bool IsZero(const ElementType &a) const {
    return a == AdditiveIdentity();
  }

  virtual bool IsOne(const ElementType &a) const {
    return a == MultiplicativeIdentity();
  }

  virtual ElementType MultiplicativeGenerator() const = 0;
  virtual std::vector<ElementType> MultiplicativeGenerators() const = 0;

  // Frobenius endomorphism: a -> a^p.
  virtual ElementType Frobenius(const ElementType &a) const {
    return Pow(a, Characteristic());
  }

  // Absolute field trace: Tr_{F_{p^m}/F_p}(a) = sum_{i=0}^{m-1} a^{p^i}.
  virtual ElementType Trace(const ElementType &a) const {
    uint32_t deg = Degree();
    ElementType acc = a;
    ElementType cur = a;
    for (uint32_t i = 1; i < deg; ++i) {
      cur = Frobenius(cur);
      acc = Add(acc, cur);
    }
    return acc;
  }

  // Absolute field norm: N_{F_{p^m}/F_p}(a) = prod_{i=0}^{m-1} a^{p^i}.
  virtual ElementType Norm(const ElementType &a) const {
    uint32_t deg = Degree();
    ElementType acc = a;
    ElementType cur = a;
    for (uint32_t i = 1; i < deg; ++i) {
      cur = Frobenius(cur);
      acc = Mul(acc, cur);
    }
    return acc;
  }

  // Exact multiplicative order of a non-zero element a in F_q^*.
  virtual uint64_t MultiplicativeOrder(const ElementType &a) const {
    if (IsZero(a)) {
      throw std::domain_error("Multiplicative order of zero is undefined");
    }
    if (IsOne(a)) {
      return 1;
    }
    uint64_t group_order = Order() - 1;
    if (Order() == 0) {
      // GF(2^64): group order is 2^64 - 1 = UINT64_MAX.
      group_order = ~uint64_t{0};
    }
    uint64_t ord = group_order;
    auto prime_factors = utils::DistinctPrimeFactors(group_order);
    for (uint64_t factor : prime_factors) {
      while (ord % factor == 0 && IsOne(Pow(a, ord / factor))) {
        ord /= factor;
      }
    }
    return ord;
  }

  virtual ElementType GetElementValue(const ElementType &value) const = 0;
  virtual ElementType SetElementValue(const ElementType &value) const = 0;
  virtual ElementType SetElementValue(const std::string &value_str) const = 0;

  virtual void Print(std::ostream &os) const = 0;
  virtual void Print(const ElementType &a, std::ostream &os) const = 0;
  virtual std::string ToString(const ElementType &a) const = 0;

  virtual FieldRepresentation GetRepresentation() const = 0;
  virtual void SetRepresentation(FieldRepresentation rep) = 0;

  virtual size_t MemoryBytes() const { return sizeof(*this); }

  // Bulk contiguous span kernels (overridden with SIMD / table / Barrett
  // implementations in concrete fields).
  virtual void AddArray(std::span<ElementType> dst,
                        std::span<const ElementType> a,
                        std::span<const ElementType> b) const {
    const size_t n = dst.size();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Add(a[i], b[i]);
    }
  }

  virtual void SubArray(std::span<ElementType> dst,
                        std::span<const ElementType> a,
                        std::span<const ElementType> b) const {
    const size_t n = dst.size();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Sub(a[i], b[i]);
    }
  }

  virtual void MulArray(std::span<ElementType> dst,
                        std::span<const ElementType> a,
                        std::span<const ElementType> b) const {
    const size_t n = dst.size();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Mul(a[i], b[i]);
    }
  }

  virtual void ScalarMulArray(std::span<ElementType> dst,
                              std::span<const ElementType> a,
                              const ElementType &scalar) const {
    const size_t n = dst.size();
    if (IsZero(scalar)) {
      const ElementType zero = AdditiveIdentity();
      for (size_t i = 0; i < n; ++i) {
        dst[i] = zero;
      }
      return;
    }
    if (IsOne(scalar)) {
      if (dst.data() != a.data()) {
        for (size_t i = 0; i < n; ++i) {
          dst[i] = a[i];
        }
      }
      return;
    }
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Mul(a[i], scalar);
    }
  }

  // Fused multiply-accumulate: dst[i] += a[i] * scalar.
  virtual void FmaArray(std::span<ElementType> dst,
                        std::span<const ElementType> a,
                        const ElementType &scalar) const {
    if (IsZero(scalar)) {
      return;
    }
    const size_t n = dst.size();
    if (IsOne(scalar)) {
      for (size_t i = 0; i < n; ++i) {
        dst[i] = Add(dst[i], a[i]);
      }
      return;
    }
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Add(dst[i], Mul(a[i], scalar));
    }
  }

  virtual ElementType DotProduct(std::span<const ElementType> a,
                                 std::span<const ElementType> b) const {
    const size_t n = a.size();
    ElementType acc = AdditiveIdentity();
    for (size_t i = 0; i < n; ++i) {
      acc = Add(acc, Mul(a[i], b[i]));
    }
    return acc;
  }

  // Montgomery's batch inversion trick: 1 Inv + 3(n-1) Mul for n elements.
  virtual void BatchInv(std::span<ElementType> dst,
                        std::span<const ElementType> src) const {
    const size_t n = src.size();
    if (dst.size() != n) {
      throw std::invalid_argument("BatchInv: dst and src sizes must match");
    }
    if (n == 0) {
      return;
    }
    std::vector<ElementType> prefix(n);
    prefix[0] = src[0];
    if (IsZero(prefix[0])) {
      throw std::domain_error("BatchInv: inverse of zero is undefined");
    }
    for (size_t i = 1; i < n; ++i) {
      if (IsZero(src[i])) {
        throw std::domain_error("BatchInv: inverse of zero is undefined");
      }
      prefix[i] = Mul(prefix[i - 1], src[i]);
    }
    ElementType inv_acc = Inv(prefix[n - 1]);
    for (size_t i = n - 1; i > 0; --i) {
      dst[i] = Mul(inv_acc, prefix[i - 1]);
      inv_acc = Mul(inv_acc, src[i]);
    }
    dst[0] = inv_acc;
  }
};

template <typename ElementType>
inline std::ostream &operator<<(std::ostream &os,
                                const GaloisFieldBase<ElementType> &field) {
  field.Print(os);
  return os;
}

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_BASE_HPP_
