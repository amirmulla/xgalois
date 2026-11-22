#ifndef XGALOIS_RING_MATRIX_RING_HPP_
#define XGALOIS_RING_MATRIX_RING_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"

namespace xg {
namespace ring {

// Non-commutative matrix ring M_n(R) over any xg::Ring R.
// Elements are flat row-major n x n matrices of length n * n.
template <Ring R>
class MatrixRing {
 public:
  using base_ring_type = R;
  using scalar_type = typename R::element_type;
  using element_type = std::vector<scalar_type>;

  MatrixRing(R base_ring, std::size_t n) : ring_(std::move(base_ring)), n_(n) {}

  const R& BaseRing() const noexcept { return ring_; }
  std::size_t Dimension() const noexcept { return n_; }

  element_type Zero() const {
    return element_type(n_ * n_, ring_.AdditiveIdentity());
  }

  element_type Eye() const {
    element_type eye(n_ * n_, ring_.AdditiveIdentity());
    const scalar_type one = ring_.MultiplicativeIdentity();
    for (std::size_t i = 0; i < n_; ++i) {
      eye[i * n_ + i] = one;
    }
    return eye;
  }

  element_type AdditiveIdentity() const { return Zero(); }
  element_type MultiplicativeIdentity() const { return Eye(); }

  element_type Add(const element_type& a, const element_type& b) const {
    ValidateSize(a);
    ValidateSize(b);
    element_type out(n_ * n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      out[i] = ring_.Add(a[i], b[i]);
    }
    return out;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    ValidateSize(a);
    ValidateSize(b);
    element_type out(n_ * n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      out[i] = ring_.Sub(a[i], b[i]);
    }
    return out;
  }

  element_type Neg(const element_type& a) const {
    ValidateSize(a);
    element_type out(n_ * n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      out[i] = ring_.Neg(a[i]);
    }
    return out;
  }

  element_type Scale(const scalar_type& s, const element_type& a) const {
    ValidateSize(a);
    element_type out(n_ * n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      out[i] = ring_.Mul(s, a[i]);
    }
    return out;
  }

  // Matrix multiplication preserving operand order a(i, k) * b(k, j).
  element_type Mul(const element_type& a, const element_type& b) const {
    ValidateSize(a);
    ValidateSize(b);
    element_type out(n_ * n_, ring_.AdditiveIdentity());
    for (std::size_t i = 0; i < n_; ++i) {
      for (std::size_t k = 0; k < n_; ++k) {
        const scalar_type& aik = a[i * n_ + k];
        if (ring_.IsZero(aik)) {
          continue;
        }
        for (std::size_t j = 0; j < n_; ++j) {
          const scalar_type& bkj = b[k * n_ + j];
          if (ring_.IsZero(bkj)) {
            continue;
          }
          out[i * n_ + j] = ring_.Add(out[i * n_ + j], ring_.Mul(aik, bkj));
        }
      }
    }
    return out;
  }

  element_type Pow(const element_type& mat, uint64_t exp) const {
    ValidateSize(mat);
    element_type acc = Eye();
    element_type base = mat;
    while (exp > 0) {
      if ((exp & uint64_t{1}) != 0) {
        acc = Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return acc;
  }

  scalar_type Trace(const element_type& mat) const {
    ValidateSize(mat);
    scalar_type acc = ring_.AdditiveIdentity();
    for (std::size_t i = 0; i < n_; ++i) {
      acc = ring_.Add(acc, mat[i * n_ + i]);
    }
    return acc;
  }

  element_type Transpose(const element_type& mat) const {
    ValidateSize(mat);
    element_type out(n_ * n_, ring_.AdditiveIdentity());
    for (std::size_t r = 0; r < n_; ++r) {
      for (std::size_t c = 0; c < n_; ++c) {
        out[c * n_ + r] = mat[r * n_ + c];
      }
    }
    return out;
  }

  // Computes the characteristic polynomial det(lambda * I - mat) =
  // c[0]*lambda^n + c[1]*lambda^{n-1} + ... + c[n] using the division-free
  // Berkowitz algorithm in O(n^4) ring operations over any CommutativeRing R.
  std::vector<scalar_type> CharacteristicPolynomial(
      const element_type& mat) const
    requires CommutativeRing<R>
  {
    ValidateSize(mat);
    const scalar_type zero = ring_.AdditiveIdentity();
    const scalar_type one = ring_.MultiplicativeIdentity();

    std::vector<scalar_type> poly{one};
    if (n_ == 0) {
      return poly;
    }

    for (std::size_t k = n_; k-- > 0;) {
      const std::size_t m = n_ - 1 - k;
      // First column of the (m + 2) x (m + 1) lower-triangular Toeplitz matrix
      // T_k: [1, -a_{k,k}, -R_k C_k, -R_k M_k C_k, ..., -R_k M_k^{m-1} C_k].
      std::vector<scalar_type> toeplitz(m + 2, zero);
      toeplitz[0] = one;
      toeplitz[1] = ring_.Neg(mat[k * n_ + k]);

      if (m > 0) {
        // v starts as column vector C_k = mat[k+1 .. n-1, k].
        std::vector<scalar_type> v(m, zero);
        for (std::size_t i = 0; i < m; ++i) {
          v[i] = mat[(k + 1 + i) * n_ + k];
        }

        for (std::size_t j = 2; j <= m + 1; ++j) {
          // dot = R_k * v where R_k = mat[k, k+1 .. n-1].
          scalar_type dot = zero;
          for (std::size_t i = 0; i < m; ++i) {
            dot = ring_.Add(dot, ring_.Mul(mat[k * n_ + (k + 1 + i)], v[i]));
          }
          toeplitz[j] = ring_.Neg(dot);

          if (j < m + 1) {
            // Advance v <- M_k * v where M_k = mat[k+1 .. n-1, k+1 .. n-1].
            std::vector<scalar_type> next_v(m, zero);
            for (std::size_t r = 0; r < m; ++r) {
              scalar_type acc = zero;
              for (std::size_t c = 0; c < m; ++c) {
                acc = ring_.Add(
                    acc, ring_.Mul(mat[(k + 1 + r) * n_ + (k + 1 + c)], v[c]));
              }
              next_v[r] = std::move(acc);
            }
            v = std::move(next_v);
          }
        }
      }

      // Multiply Toeplitz matrix T_k by poly.
      std::vector<scalar_type> next_poly(m + 2, zero);
      for (std::size_t r = 0; r < m + 2; ++r) {
        scalar_type acc = zero;
        const std::size_t c_end = std::min(r, m);
        for (std::size_t c = 0; c <= c_end; ++c) {
          acc = ring_.Add(acc, ring_.Mul(toeplitz[r - c], poly[c]));
        }
        next_poly[r] = std::move(acc);
      }
      poly = std::move(next_poly);
    }

    return poly;
  }

  // Exact division-free determinant via the Berkowitz algorithm (O(n^4)),
  // valid over ANY CommutativeRing (including rings with zero divisors).
  scalar_type Det(const element_type& mat) const
    requires CommutativeRing<R>
  {
    const auto char_poly = CharacteristicPolynomial(mat);
    const scalar_type& constant_term = char_poly.back();
    return ((n_ & std::size_t{1}) != 0) ? ring_.Neg(constant_term)
                                        : constant_term;
  }

  bool IsZero(const element_type& a) const {
    ValidateSize(a);
    for (const auto& entry : a) {
      if (!ring_.IsZero(entry)) {
        return false;
      }
    }
    return true;
  }

  bool IsOne(const element_type& a) const {
    ValidateSize(a);
    return Equals(a, Eye());
  }

  bool Equals(const element_type& a, const element_type& b) const {
    ValidateSize(a);
    ValidateSize(b);
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      if (!alg::Equals(ring_, a[i], b[i])) {
        return false;
      }
    }
    return true;
  }

  std::string ToString(const element_type& a) const {
    ValidateSize(a);
    std::string out = "[";
    for (std::size_t r = 0; r < n_; ++r) {
      if (r > 0) {
        out += "; ";
      }
      out += "[";
      for (std::size_t c = 0; c < n_; ++c) {
        if (c > 0) {
          out += ", ";
        }
        out += ring_.ToString(a[r * n_ + c]);
      }
      out += "]";
    }
    out += "]";
    return out;
  }

 private:
  void ValidateSize(const element_type& mat) const {
    if (mat.size() != n_ * n_) {
      throw std::invalid_argument(
          "MatrixRing: matrix vector size does not match n * n");
    }
  }

  R ring_;
  std::size_t n_;
};

}  // namespace ring

template <Ring R>
struct AlgebraTraits<ring::MatrixRing<R>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

}  // namespace xg

#endif  // XGALOIS_RING_MATRIX_RING_HPP_
