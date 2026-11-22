#ifndef XGALOIS_ALGEBRA_MATRIX_GROUP_HPP_
#define XGALOIS_ALGEBRA_MATRIX_GROUP_HPP_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/group_base.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace algebra {

enum class MatrixGroupKind : uint8_t {
  kGeneralLinear,
  kSpecialLinear,
};

// ============================================================================
// MatrixGroup<F>: GL_n(F_q) and SL_n(F_q) over any FiniteField F
// ============================================================================
//
// Elements are flat row-major n x n vectors of F::storage_type, making them
// lightweight, copyable, and directly compatible with GroupElement<G> and
// Subgroup<G>.
template <FiniteField F>
class MatrixGroup {
 public:
  using field_type = F;
  using storage_type = typename F::storage_type;
  using element_type = std::vector<storage_type>;

  MatrixGroup(std::shared_ptr<F> field, std::size_t n,
              MatrixGroupKind kind = MatrixGroupKind::kGeneralLinear)
      : field_(std::move(field)), n_(n), kind_(kind) {
    if (!field_) {
      throw std::invalid_argument("MatrixGroup: field pointer cannot be null");
    }
    if (n_ == 0) {
      throw std::invalid_argument("MatrixGroup: dimension n must be >= 1");
    }
  }

  MatrixGroup(F field, std::size_t n,
              MatrixGroupKind kind = MatrixGroupKind::kGeneralLinear)
      : MatrixGroup(std::make_shared<F>(std::move(field)), n, kind) {}

  std::shared_ptr<F> Field() const noexcept { return field_; }
  std::size_t Dimension() const noexcept { return n_; }
  std::size_t N() const noexcept { return n_; }
  MatrixGroupKind Kind() const noexcept { return kind_; }

  element_type Identity() const {
    element_type eye(n_ * n_, field_->AdditiveIdentity());
    const storage_type one = field_->MultiplicativeIdentity();
    for (std::size_t i = 0; i < n_; ++i) {
      eye[i * n_ + i] = one;
    }
    return eye;
  }

  bool IsIdentity(const element_type& mat) const {
    RequireShape(mat);
    for (std::size_t i = 0; i < n_; ++i) {
      for (std::size_t j = 0; j < n_; ++j) {
        const storage_type val = mat[i * n_ + j];
        if (i == j) {
          if (!field_->IsOne(val)) {
            return false;
          }
        } else {
          if (!field_->IsZero(val)) {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool Equals(const element_type& a, const element_type& b) const {
    RequireShape(a);
    RequireShape(b);
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      if (field_->GetElementValue(a[i]) != field_->GetElementValue(b[i])) {
        return false;
      }
    }
    return true;
  }

  storage_type Determinant(const element_type& mat) const {
    RequireShape(mat);
    element_type work = Normalize(mat);
    storage_type det = field_->MultiplicativeIdentity();
    bool negate = false;

    for (std::size_t col = 0; col < n_; ++col) {
      std::size_t pivot = col;
      while (pivot < n_ && field_->IsZero(work[pivot * n_ + col])) {
        ++pivot;
      }
      if (pivot == n_) {
        return field_->AdditiveIdentity();
      }
      if (pivot != col) {
        for (std::size_t j = col; j < n_; ++j) {
          std::swap(work[col * n_ + j], work[pivot * n_ + j]);
        }
        negate = !negate;
      }
      const storage_type diag = work[col * n_ + col];
      det = field_->Mul(det, diag);
      const storage_type inv_diag = field_->Inv(diag);

      for (std::size_t row = col + 1; row < n_; ++row) {
        const storage_type entry = work[row * n_ + col];
        if (field_->IsZero(entry)) {
          continue;
        }
        const storage_type factor = field_->Mul(entry, inv_diag);
        work[row * n_ + col] = field_->AdditiveIdentity();
        for (std::size_t j = col + 1; j < n_; ++j) {
          work[row * n_ + j] = field_->Sub(
              work[row * n_ + j], field_->Mul(factor, work[col * n_ + j]));
        }
      }
    }

    return negate ? field_->Neg(det) : det;
  }

  bool Contains(const element_type& mat) const {
    if (mat.size() != n_ * n_) {
      return false;
    }
    const storage_type det = Determinant(mat);
    if (kind_ == MatrixGroupKind::kGeneralLinear) {
      return !field_->IsZero(det);
    }
    return field_->IsOne(det);
  }

  element_type Op(const element_type& a, const element_type& b) const {
    RequireShape(a);
    RequireShape(b);
    const storage_type zero = field_->AdditiveIdentity();
    element_type out(n_ * n_, zero);
    for (std::size_t i = 0; i < n_; ++i) {
      for (std::size_t k = 0; k < n_; ++k) {
        const storage_type aik = field_->GetElementValue(a[i * n_ + k]);
        if (field_->IsZero(aik)) {
          continue;
        }
        for (std::size_t j = 0; j < n_; ++j) {
          const storage_type bkj = field_->GetElementValue(b[k * n_ + j]);
          if (!field_->IsZero(bkj)) {
            out[i * n_ + j] =
                field_->Add(out[i * n_ + j], field_->Mul(aik, bkj));
          }
        }
      }
    }
    return out;
  }

  element_type Inv(const element_type& mat) const {
    RequireShape(mat);
    element_type lhs = Normalize(mat);
    element_type rhs = Identity();

    for (std::size_t col = 0; col < n_; ++col) {
      std::size_t pivot = col;
      while (pivot < n_ && field_->IsZero(lhs[pivot * n_ + col])) {
        ++pivot;
      }
      if (pivot == n_) {
        throw std::domain_error(
            "MatrixGroup::Inv: matrix is singular and not invertible");
      }
      if (pivot != col) {
        for (std::size_t j = 0; j < n_; ++j) {
          std::swap(lhs[col * n_ + j], lhs[pivot * n_ + j]);
          std::swap(rhs[col * n_ + j], rhs[pivot * n_ + j]);
        }
      }

      const storage_type diag = lhs[col * n_ + col];
      if (!field_->IsOne(diag)) {
        const storage_type inv_diag = field_->Inv(diag);
        for (std::size_t j = col; j < n_; ++j) {
          lhs[col * n_ + j] = field_->Mul(lhs[col * n_ + j], inv_diag);
        }
        for (std::size_t j = 0; j < n_; ++j) {
          rhs[col * n_ + j] = field_->Mul(rhs[col * n_ + j], inv_diag);
        }
      }

      for (std::size_t row = 0; row < n_; ++row) {
        if (row == col) {
          continue;
        }
        const storage_type factor = lhs[row * n_ + col];
        if (field_->IsZero(factor)) {
          continue;
        }
        lhs[row * n_ + col] = field_->AdditiveIdentity();
        for (std::size_t j = col + 1; j < n_; ++j) {
          lhs[row * n_ + j] = field_->Sub(
              lhs[row * n_ + j], field_->Mul(factor, lhs[col * n_ + j]));
        }
        for (std::size_t j = 0; j < n_; ++j) {
          rhs[row * n_ + j] = field_->Sub(
              rhs[row * n_ + j], field_->Mul(factor, rhs[col * n_ + j]));
        }
      }
    }
    return rhs;
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return Op(a, Inv(b));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    RequireShape(a);
    element_type acc = Identity();
    if (exp == 0) {
      return acc;
    }
    element_type base = Normalize(a);
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        acc = Op(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Op(base, base);
      }
    }
    return acc;
  }

  // Exact group order:
  //   |GL_n(F_q)| = prod_{i=0}^{n-1} (q^n - q^i)
  //   |SL_n(F_q)| = |GL_n(F_q)| / (q - 1)
  // Safe even when field_->Order() == 0 (e.g. GF(2^64)).
  BigInt Order() const {
    const BigInt q = FieldBigOrder();
    const BigInt q_n = q.Pow(static_cast<uint64_t>(n_));
    BigInt ord(1);
    BigInt q_pow(1);
    for (std::size_t i = 0; i < n_; ++i) {
      ord *= (q_n - q_pow);
      q_pow *= q;
    }
    if (kind_ == MatrixGroupKind::kSpecialLinear) {
      ord /= (q - BigInt(1));
    }
    return ord;
  }

  // Generating set:
  // - SL_n(F_q): elementary transvections I + alpha * E_{ij} (i != j) for
  //   alpha in an F_p-basis {1, g, ..., g^{m-1}} of F_{p^m}.
  // - GL_n(F_q): all SL_n generators plus diag(g, 1, ..., 1) where g generates
  //   F_q^* (when q > 2).
  std::vector<element_type> Generators() const {
    std::vector<element_type> gens;
    const BigInt q = FieldBigOrder();
    const storage_type gen_elem = field_->MultiplicativeGenerator();

    if (n_ == 1) {
      if (kind_ == MatrixGroupKind::kSpecialLinear || q <= BigInt(2)) {
        gens.push_back(Identity());
      } else {
        gens.push_back(element_type{gen_elem});
      }
      return gens;
    }

    // Additive F_p-basis of F_{p^m}: {1, g, ..., g^{m-1}}.
    const uint32_t deg = field_->Degree();
    std::vector<storage_type> basis;
    basis.reserve(deg);
    storage_type cur = field_->MultiplicativeIdentity();
    for (uint32_t k = 0; k < deg; ++k) {
      basis.push_back(cur);
      if (k + 1 < deg) {
        cur = field_->Mul(cur, gen_elem);
      }
    }

    // Elementary transvections I + alpha * E_{ij} (i != j).
    const element_type eye = Identity();
    for (std::size_t i = 0; i < n_; ++i) {
      for (std::size_t j = 0; j < n_; ++j) {
        if (i == j) {
          continue;
        }
        for (const storage_type& alpha : basis) {
          element_type trans = eye;
          trans[i * n_ + j] = alpha;
          gens.push_back(std::move(trans));
        }
      }
    }

    // Diagonal scaling generator diag(g, 1, ..., 1) for GL_n when q > 2.
    if (kind_ == MatrixGroupKind::kGeneralLinear && q > BigInt(2)) {
      element_type diag = eye;
      diag[0] = gen_elem;
      gens.push_back(std::move(diag));
    }

    return gens;
  }

  std::string ToString(const element_type& mat) const {
    RequireShape(mat);
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < n_; ++i) {
      if (i > 0) {
        oss << ", ";
      }
      oss << "[";
      for (std::size_t j = 0; j < n_; ++j) {
        if (j > 0) {
          oss << ", ";
        }
        oss << field_->ToString(field_->GetElementValue(mat[i * n_ + j]));
      }
      oss << "]";
    }
    oss << "]";
    return oss.str();
  }

 private:
  void RequireShape(const element_type& mat) const {
    if (mat.size() != n_ * n_) {
      throw std::invalid_argument(
          "MatrixGroup: matrix size does not match n * n");
    }
  }

  element_type Normalize(const element_type& mat) const {
    element_type out(n_ * n_);
    for (std::size_t i = 0; i < n_ * n_; ++i) {
      out[i] = field_->GetElementValue(mat[i]);
    }
    return out;
  }

  BigInt FieldBigOrder() const {
    const uint64_t ord = static_cast<uint64_t>(field_->Order());
    if (ord > 0) {
      return BigInt::FromUint64(ord);
    }
    const BigInt ch =
        BigInt::FromUint64(static_cast<uint64_t>(field_->Characteristic()));
    return ch.Pow(static_cast<uint64_t>(field_->Degree()));
  }

  std::shared_ptr<F> field_;
  std::size_t n_ = 1;
  MatrixGroupKind kind_ = MatrixGroupKind::kGeneralLinear;
};

}  // namespace algebra

// ============================================================================
// AlgebraTraits specialization in namespace xg
// ============================================================================

template <FiniteField F>
struct AlgebraTraits<algebra::MatrixGroup<F>> {
  static constexpr bool kIsAssociativeOp = true;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = false;
  static constexpr bool kHasNoZeroDivisors = false;
  static constexpr bool kIsPid = false;
  static constexpr bool kIsEuclideanDomain = false;
  static constexpr bool kIsField = false;
};

using algebra::MatrixGroup;
using algebra::MatrixGroupKind;

}  // namespace xg

#endif  // XGALOIS_ALGEBRA_MATRIX_GROUP_HPP_
