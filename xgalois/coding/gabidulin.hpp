#ifndef XGALOIS_CODING_GABIDULIN_HPP_
#define XGALOIS_CODING_GABIDULIN_HPP_

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/ring/skew_polynomial_ring.hpp"

namespace xg {
namespace coding {

namespace detail {

template <typename T>
concept HasBit128Interface = requires(const T &x) {
  { x.Bit(0) } -> std::convertible_to<bool>;
  { x.lo } -> std::convertible_to<uint64_t>;
  { x.hi } -> std::convertible_to<uint64_t>;
};

template <typename FieldT, typename StorageT>
concept HasExtensionCoeffsInterface =
    requires(const FieldT &f, const StorageT &x) {
      typename StorageT::coeff_type;
      { x.Coeffs() };
      {
        f.FromCoeffs(std::vector<typename StorageT::coeff_type>{})
      } -> std::same_as<StorageT>;
    };

template <typename BaseField, typename ExtField>
concept IsDirectRelativeTower =
    requires(const ExtField &ef, const typename ExtField::element_type &x) {
      { ef.RelativeDegree() } -> std::convertible_to<std::size_t>;
      { x.Coeffs() };
      requires std::same_as<typename ExtField::base_field_type, BaseField>;
    };

inline uint64_t ModAdd(uint64_t a, uint64_t b, uint64_t p) noexcept {
  const uint64_t s = a + b;
  return (s >= p) ? (s - p) : s;
}

inline uint64_t ModSub(uint64_t a, uint64_t b, uint64_t p) noexcept {
  return (a >= b) ? (a - b) : (p - (b - a));
}

__extension__ using u128 = unsigned __int128;

inline uint64_t ModMul(uint64_t a, uint64_t b, uint64_t p) noexcept {
  return static_cast<uint64_t>((static_cast<u128>(a) * static_cast<u128>(b)) %
                               p);
}

inline uint64_t ModPow(uint64_t base, uint64_t exp, uint64_t p) noexcept {
  uint64_t acc = 1 % p;
  uint64_t cur = base % p;
  while (exp > 0) {
    if ((exp & 1ULL) != 0) {
      acc = ModMul(acc, cur, p);
    }
    exp >>= 1;
    if (exp > 0) {
      cur = ModMul(cur, cur, p);
    }
  }
  return acc;
}

inline uint64_t ModInv(uint64_t a, uint64_t p) {
  a %= p;
  if (a == 0) {
    throw std::domain_error("FieldTowerAdapter: modular inverse of zero");
  }
  return ModPow(a, p - 2, p);
}

inline std::vector<uint64_t> MatVecModP(const std::vector<uint64_t> &mat,
                                        std::size_t rows, std::size_t cols,
                                        const std::vector<uint64_t> &vec,
                                        uint64_t p) {
  std::vector<uint64_t> out(rows, 0);
  for (std::size_t r = 0; r < rows; ++r) {
    uint64_t acc = 0;
    const std::size_t offset = r * cols;
    for (std::size_t c = 0; c < cols; ++c) {
      if (mat[offset + c] != 0 && vec[c] != 0) {
        acc = ModAdd(acc, ModMul(mat[offset + c], vec[c], p), p);
      }
    }
    out[r] = acc;
  }
  return out;
}

inline std::optional<std::vector<uint64_t>> InvertModP(
    std::vector<uint64_t> mat, std::size_t n, uint64_t p) {
  std::vector<uint64_t> inv(n * n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    inv[i * n + i] = 1 % p;
  }

  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    while (pivot < n && mat[pivot * n + col] % p == 0) {
      ++pivot;
    }
    if (pivot == n) {
      return std::nullopt;
    }
    if (pivot != col) {
      for (std::size_t j = 0; j < n; ++j) {
        std::swap(mat[col * n + j], mat[pivot * n + j]);
        std::swap(inv[col * n + j], inv[pivot * n + j]);
      }
    }
    const uint64_t inv_piv = ModInv(mat[col * n + col], p);
    for (std::size_t j = 0; j < n; ++j) {
      mat[col * n + j] = ModMul(mat[col * n + j], inv_piv, p);
      inv[col * n + j] = ModMul(inv[col * n + j], inv_piv, p);
    }
    for (std::size_t r = 0; r < n; ++r) {
      if (r == col) {
        continue;
      }
      const uint64_t factor = mat[r * n + col] % p;
      if (factor == 0) {
        continue;
      }
      for (std::size_t j = 0; j < n; ++j) {
        mat[r * n + j] =
            ModSub(mat[r * n + j], ModMul(factor, mat[col * n + j], p), p);
        inv[r * n + j] =
            ModSub(inv[r * n + j], ModMul(factor, inv[col * n + j], p), p);
      }
    }
  }
  return inv;
}

inline std::size_t RankModP(std::vector<uint64_t> mat, std::size_t rows,
                            std::size_t cols, uint64_t p) {
  std::size_t pivot_row = 0;
  for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
    std::size_t found = pivot_row;
    while (found < rows && mat[found * cols + col] % p == 0) {
      ++found;
    }
    if (found == rows) {
      continue;
    }
    if (found != pivot_row) {
      for (std::size_t j = col; j < cols; ++j) {
        std::swap(mat[pivot_row * cols + j], mat[found * cols + j]);
      }
    }
    const uint64_t inv_piv = ModInv(mat[pivot_row * cols + col], p);
    for (std::size_t r = pivot_row + 1; r < rows; ++r) {
      const uint64_t factor = ModMul(mat[r * cols + col], inv_piv, p);
      if (factor == 0) {
        continue;
      }
      mat[r * cols + col] = 0;
      for (std::size_t j = col + 1; j < cols; ++j) {
        mat[r * cols + j] = ModSub(
            mat[r * cols + j], ModMul(factor, mat[pivot_row * cols + j], p), p);
      }
    }
    ++pivot_row;
  }
  return pivot_row;
}

// Extracts F_p coordinates of an element in F_{p^d} (length d).
template <typename FieldT>
std::vector<uint64_t> PrimeFieldCoords(
    const FieldT &field, const typename FieldT::element_type &elem) {
  using StorageT = typename FieldT::element_type;
  const std::size_t d = static_cast<std::size_t>(field.Degree());
  std::vector<uint64_t> coords(d, 0);

  if constexpr (HasBit128Interface<StorageT>) {
    const StorageT val = field.GetElementValue(elem);
    for (std::size_t i = 0; i < d; ++i) {
      coords[i] = val.Bit(static_cast<int>(i)) ? 1ULL : 0ULL;
    }
  } else if constexpr (HasExtensionCoeffsInterface<FieldT, StorageT>) {
    const StorageT canonical = field.GetElementValue(elem);
    const auto &coeffs = canonical.Coeffs();
    const std::size_t lim = std::min(d, coeffs.size());
    for (std::size_t i = 0; i < lim; ++i) {
      coords[i] = static_cast<uint64_t>(coeffs[i]);
    }
  } else if constexpr (std::is_integral_v<StorageT>) {
    const uint64_t val = static_cast<uint64_t>(field.GetElementValue(elem));
    if (d == 1) {
      coords[0] = val % field.Characteristic();
    } else {
      // Binary extension field F_{2^d}: polynomial-basis bitmask.
      for (std::size_t i = 0; i < d; ++i) {
        coords[i] = (val >> i) & 1ULL;
      }
    }
  } else {
    throw std::invalid_argument(
        "FieldTowerAdapter: unsupported field element representation for "
        "prime-field coordinate extraction");
  }
  return coords;
}

// Reconstructs an element of F_{p^d} from its length-d F_p coordinate vector.
template <typename FieldT>
typename FieldT::element_type FromPrimeFieldCoords(
    const FieldT &field, const std::vector<uint64_t> &coords,
    std::size_t offset = 0) {
  using StorageT = typename FieldT::element_type;
  const std::size_t d = static_cast<std::size_t>(field.Degree());

  if constexpr (HasBit128Interface<StorageT>) {
    StorageT val{};
    for (std::size_t i = 0; i < d; ++i) {
      if ((coords[offset + i] & 1ULL) != 0) {
        val |= (StorageT{1, 0} << static_cast<int>(i));
      }
    }
    return field.SetElementValue(val);
  } else if constexpr (HasExtensionCoeffsInterface<FieldT, StorageT>) {
    using CoeffT = typename StorageT::coeff_type;
    std::vector<CoeffT> coeffs(d);
    const uint64_t p = field.Characteristic();
    for (std::size_t i = 0; i < d; ++i) {
      coeffs[i] = static_cast<CoeffT>(coords[offset + i] % p);
    }
    return field.FromCoeffs(std::move(coeffs));
  } else if constexpr (std::is_integral_v<StorageT>) {
    if (d == 1) {
      const uint64_t p = field.Characteristic();
      return field.SetElementValue(static_cast<StorageT>(coords[offset] % p));
    }
    uint64_t bits = 0;
    for (std::size_t i = 0; i < d; ++i) {
      if ((coords[offset + i] & 1ULL) != 0) {
        bits |= (uint64_t{1} << i);
      }
    }
    return field.SetElementValue(static_cast<StorageT>(bits));
  } else {
    throw std::invalid_argument(
        "FieldTowerAdapter: unsupported field element representation for "
        "prime-field coordinate reconstruction");
  }
}

// Solves the (possibly rectangular or rank-deficient) linear system A * X = B
// over any finite field using linalg::Rref([A | B]).
// Returns a particular solution X (shape {cols_A, cols_B}) with free variables
// set to zero, or std::nullopt if inconsistent.
template <typename GaloisField>
std::optional<GFMatrix<GaloisField>> SolveLinearSystem(
    const GFMatrix<GaloisField> &a, const GFMatrix<GaloisField> &b) {
  using StorageT = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols_a = a.Shape(1);
  const std::size_t cols_b = b.Shape(1);

  if (cols_b == 0) {
    return GFMatrix<GaloisField>::Zeros(field, {cols_a, 0});
  }
  if (rows == 0) {
    return GFMatrix<GaloisField>::Zeros(field, {cols_a, cols_b});
  }
  if (cols_a == 0) {
    for (const auto &v : b.RawSpan()) {
      if (!field->IsZero(v)) {
        return std::nullopt;
      }
    }
    return GFMatrix<GaloisField>::Zeros(field, {0, cols_b});
  }

  const GFMatrix<GaloisField> aug = linalg::HStack(a, b);
  const GFMatrix<GaloisField> rref = linalg::Rref(aug);
  const std::size_t total_cols = cols_a + cols_b;

  RawArray x_raw = RawArray::from_shape(ShapeType{cols_a, cols_b});
  x_raw.fill(field->AdditiveIdentity());
  StorageT *x_data = x_raw.data();
  const StorageT *r_data = rref.RawData().data();

  for (std::size_t r = 0; r < rows; ++r) {
    const StorageT *row_ptr = r_data + r * total_cols;
    std::size_t pivot = total_cols;
    for (std::size_t c = 0; c < total_cols; ++c) {
      if (!field->IsZero(row_ptr[c])) {
        pivot = c;
        break;
      }
    }
    if (pivot == total_cols) {
      continue;  // Zero row.
    }
    if (pivot >= cols_a) {
      return std::nullopt;  // Inconsistent row [0 ... 0 | non-zero].
    }
    // Pivot column `pivot` has leading 1; assign x[pivot, j] = rref[r,
    // cols_a+j].
    StorageT *x_row = x_data + pivot * cols_b;
    for (std::size_t j = 0; j < cols_b; ++j) {
      x_row[j] = row_ptr[cols_a + j];
    }
  }

  return GFMatrix<GaloisField>::FromRaw(std::move(field), std::move(x_raw));
}

// Computes exact subgroup exponent (q^m - 1) / (q - 1) for q = p^s and
// total degree M = s * m (handles M <= 64 cleanly).
inline uint64_t SubfieldGroupCofactor(uint64_t q, std::size_t m) {
  // (q^m - 1) / (q - 1) = 1 + q + q^2 + ... + q^{m-1}.
  uint64_t sum = 0;
  uint64_t pow = 1;
  for (std::size_t i = 0; i < m; ++i) {
    sum += pow;
    if (i + 1 < m) {
      pow *= q;
    }
  }
  return sum;
}

}  // namespace detail

// ============================================================================
// FieldTowerAdapter: F_{q^m} <-> F_q^m vector-space isomorphism
// ============================================================================

template <typename BaseField, typename ExtField>
class FieldTowerAdapter {
 public:
  using base_field_type = BaseField;
  using ext_field_type = ExtField;
  using base_storage_type = typename BaseField::element_type;
  using ext_storage_type = typename ExtField::element_type;
  using base_element_type = GaloisFieldElement<BaseField>;
  using ext_element_type = GaloisFieldElement<ExtField>;

  FieldTowerAdapter(std::shared_ptr<BaseField> base_field,
                    std::shared_ptr<ExtField> ext_field)
      : base_field_(std::move(base_field)), ext_field_(std::move(ext_field)) {
    if (!base_field_ || !ext_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter: base and extension fields cannot be null");
    }
    const uint64_t p_base = base_field_->Characteristic();
    const uint64_t p_ext = ext_field_->Characteristic();
    if (p_base != p_ext) {
      throw std::invalid_argument(
          "FieldTowerAdapter: base and extension fields must have the same "
          "characteristic");
    }
    const uint32_t d_base = base_field_->Degree();
    const uint32_t d_ext = ext_field_->Degree();
    if (d_base == 0 || d_ext == 0 || d_ext < d_base || (d_ext % d_base) != 0) {
      throw std::invalid_argument(
          "FieldTowerAdapter: ExtField::Degree() must be a positive multiple "
          "of BaseField::Degree()");
    }
    characteristic_ = p_base;
    base_degree_ = static_cast<std::size_t>(d_base);
    ext_degree_ = static_cast<std::size_t>(d_ext);
    relative_degree_ = ext_degree_ / base_degree_;
    base_order_ = base_field_->Order();

    InitializeTowerMapping();
  }

  std::shared_ptr<BaseField> BaseFieldPtr() const { return base_field_; }
  std::shared_ptr<ExtField> ExtFieldPtr() const { return ext_field_; }

  std::size_t RelativeDegree() const noexcept { return relative_degree_; }
  uint64_t BaseOrder() const noexcept { return base_order_; }
  uint64_t Characteristic() const noexcept { return characteristic_; }
  std::size_t BaseDegree() const noexcept { return base_degree_; }

  // Relative Frobenius automorphism: elem -> elem^{q^{steps}} =
  // elem^{p^{s*steps}}.
  ext_storage_type FrobeniusQ(const ext_storage_type &elem,
                              std::size_t steps = 1) const {
    const std::size_t red = steps % relative_degree_;
    if (red == 0 || ext_field_->IsZero(elem)) {
      return elem;
    }
    const std::size_t prime_steps = red * base_degree_;
    ext_storage_type cur = elem;
    for (std::size_t i = 0; i < prime_steps; ++i) {
      cur = ext_field_->Frobenius(cur);
    }
    return cur;
  }

  ext_element_type FrobeniusQ(const ext_element_type &elem,
                              std::size_t steps = 1) const {
    if (!elem.HasField() || elem.Field() != ext_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::FrobeniusQ: field mismatch");
    }
    return ext_element_type::FromInternal(FrobeniusQ(elem.RawValue(), steps),
                                          ext_field_);
  }

  // Inverse relative Frobenius: elem -> elem^{q^{-steps}} = elem^{q^{m - (steps
  // % m)}}.
  ext_storage_type InvFrobeniusQ(const ext_storage_type &elem,
                                 std::size_t steps = 1) const {
    const std::size_t red = steps % relative_degree_;
    if (red == 0 || ext_field_->IsZero(elem)) {
      return elem;
    }
    return FrobeniusQ(elem, relative_degree_ - red);
  }

  ext_element_type InvFrobeniusQ(const ext_element_type &elem,
                                 std::size_t steps = 1) const {
    if (!elem.HasField() || elem.Field() != ext_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::InvFrobeniusQ: field mismatch");
    }
    return ext_element_type::FromInternal(InvFrobeniusQ(elem.RawValue(), steps),
                                          ext_field_);
  }

  // Expands a single element of F_{q^m} into a length-m vector over F_q.
  GFArray<BaseField> ExpandElement(const ext_storage_type &elem) const {
    using RawArray = typename GFArray<BaseField>::raw_array_type;
    using ShapeType = typename GFArray<BaseField>::shape_type;

    RawArray raw = RawArray::from_shape(ShapeType{relative_degree_});
    ExpandElementInto(elem, raw.data(), 1);
    return GFArray<BaseField>::FromRaw(base_field_, std::move(raw));
  }

  // Contracts a length-m vector over F_q into an element of F_{q^m}.
  ext_storage_type ContractElement(const GFArray<BaseField> &vec) const {
    if (!vec.HasField() || vec.Field() != base_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::ContractElement: field mismatch");
    }
    if (vec.Dimension() != 1 || vec.Size() != relative_degree_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::ContractElement: vector must be 1D of length "
          "RelativeDegree()");
    }
    return ContractElementFrom(vec.RawData().data(), 1);
  }

  // Returns the i-th canonical F_q-basis element of F_{q^m} (0 <= i < m).
  // For standard prime/binary extensions, this is alpha^i.
  ext_storage_type CanonicalBasisElement(std::size_t index) const {
    if (index >= relative_degree_) {
      throw std::out_of_range(
          "FieldTowerAdapter::CanonicalBasisElement: index out of range");
    }
    std::vector<base_storage_type> unit(relative_degree_,
                                        base_field_->AdditiveIdentity());
    unit[index] = base_field_->MultiplicativeIdentity();
    return ContractElementFrom(unit.data(), 1);
  }

  // Expands a length-n vector over F_{q^m} into an {m, n} matrix over F_q
  // (column j is the coordinate vector of vec[j]).
  GFMatrix<BaseField> ExpandVector(const GFArray<ExtField> &vec) const {
    if (!vec.HasField() || vec.Field() != ext_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::ExpandVector: field mismatch");
    }
    if (vec.Dimension() != 1) {
      throw std::invalid_argument(
          "FieldTowerAdapter::ExpandVector: input must be a 1D vector");
    }
    const std::size_t n = vec.Size();
    using RawArray = typename GFMatrix<BaseField>::raw_array_type;
    using ShapeType = typename GFMatrix<BaseField>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{relative_degree_, n});
    if (n == 0) {
      return GFMatrix<BaseField>::FromRaw(base_field_, std::move(raw));
    }
    const ext_storage_type *src = vec.RawData().data();
    base_storage_type *dst = raw.data();
    for (std::size_t j = 0; j < n; ++j) {
      ExpandElementInto(src[j], dst + j, n);
    }
    return GFMatrix<BaseField>::FromRaw(base_field_, std::move(raw));
  }

  // Contracts an {m, n} matrix over F_q into a length-n vector over F_{q^m}.
  GFArray<ExtField> ContractMatrix(const GFMatrix<BaseField> &mat) const {
    if (!mat.HasField() || mat.Field() != base_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::ContractMatrix: field mismatch");
    }
    if (mat.Dimension() != 2 || mat.Shape(0) != relative_degree_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::ContractMatrix: matrix must have shape {m, n}");
    }
    const std::size_t n = mat.Shape(1);
    using RawArray = typename GFArray<ExtField>::raw_array_type;
    using ShapeType = typename GFArray<ExtField>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{n});
    if (n == 0) {
      return GFArray<ExtField>::FromRaw(ext_field_, std::move(raw));
    }
    const base_storage_type *src = mat.RawData().data();
    ext_storage_type *dst = raw.data();
    for (std::size_t j = 0; j < n; ++j) {
      dst[j] = ContractElementFrom(src + j, n);
    }
    return GFArray<ExtField>::FromRaw(ext_field_, std::move(raw));
  }

  // Exact rank weight over F_q: rank(ExpandVector(vec)).
  std::size_t RankWeight(const GFArray<ExtField> &vec) const {
    if (!vec.HasField() || vec.Field() != ext_field_) {
      throw std::invalid_argument(
          "FieldTowerAdapter::RankWeight: field mismatch");
    }
    if (vec.Dimension() != 1) {
      throw std::invalid_argument(
          "FieldTowerAdapter::RankWeight: input must be a 1D vector");
    }
    if (vec.Size() == 0) {
      return 0;
    }
    return linalg::MatrixRank(ExpandVector(vec));
  }

 private:
  void InitializeTowerMapping() {
    // Fast direct paths need no change-of-basis matrices.
    if constexpr (std::is_same_v<BaseField, ExtField>) {
      if (relative_degree_ == 1) {
        use_general_tower_ = false;
        return;
      }
    }
    if constexpr (detail::IsDirectRelativeTower<BaseField, ExtField>) {
      use_general_tower_ = false;
      return;
    }
    if (base_degree_ == 1) {
      use_general_tower_ = false;
      return;
    }

    // General relative tower: F_{p^s} subset F_{p^{s*m}} with s > 1.
    // Construct an exact F_p-linear change-of-basis matrix T (M x M) and its
    // inverse T^{-1} realizing an F_q-linear isomorphism F_{q^m} <-> F_q^m.
    use_general_tower_ = true;
    const uint64_t p = characteristic_;
    const std::size_t s = base_degree_;
    const std::size_t M = ext_degree_;
    const std::size_t m = relative_degree_;

    if (base_order_ == 0 || base_order_ < 2) {
      throw std::invalid_argument(
          "FieldTowerAdapter: base field order too large for general relative "
          "tower embedding");
    }

    // Step 1: Embed BaseField -> ExtField via primitive element minimal poly.
    const base_storage_type gamma_b = base_field_->MultiplicativeGenerator();

    // Compute minimal polynomial f(X) = prod_{j=0}^{s-1} (X - gamma_b^{p^j}) in
    // F_p[X] using BaseField arithmetic (all coefficients lie in F_p).
    std::vector<base_storage_type> min_poly_b(
        1, base_field_->MultiplicativeIdentity());
    base_storage_type conj = gamma_b;
    for (std::size_t j = 0; j < s; ++j) {
      const base_storage_type neg_conj = base_field_->Neg(conj);
      std::vector<base_storage_type> next(min_poly_b.size() + 1,
                                          base_field_->AdditiveIdentity());
      for (std::size_t idx = 0; idx < min_poly_b.size(); ++idx) {
        next[idx] = base_field_->Add(
            next[idx], base_field_->Mul(min_poly_b[idx], neg_conj));
        next[idx + 1] = base_field_->Add(next[idx + 1], min_poly_b[idx]);
      }
      min_poly_b = std::move(next);
      conj = base_field_->Frobenius(conj);
    }

    // Extract F_p coefficients of min_poly_b.
    std::vector<uint64_t> min_poly_fp(s + 1, 0);
    for (std::size_t i = 0; i <= s; ++i) {
      min_poly_fp[i] = detail::PrimeFieldCoords(*base_field_, min_poly_b[i])[0];
    }

    // Find a root gamma_tilde of min_poly_fp in the unique order-q subfield of
    // ExtField: subgroup <eta> where eta = gamma_e^{(q^m - 1) / (q - 1)}.
    const ext_storage_type gamma_e = ext_field_->MultiplicativeGenerator();
    const uint64_t cofactor = detail::SubfieldGroupCofactor(base_order_, m);
    const ext_storage_type eta = ext_field_->Pow(gamma_e, cofactor);

    // Helper: embed F_p scalar c into ExtField.
    auto fp_to_ext = [&](uint64_t c) -> ext_storage_type {
      std::vector<uint64_t> coords(M, 0);
      coords[0] = c % p;
      return detail::FromPrimeFieldCoords(*ext_field_, coords);
    };

    auto eval_min_poly_in_ext =
        [&](const ext_storage_type &x) -> ext_storage_type {
      ext_storage_type acc = ext_field_->AdditiveIdentity();
      for (std::size_t i = s + 1; i-- > 0;) {
        acc =
            ext_field_->Add(ext_field_->Mul(acc, x), fp_to_ext(min_poly_fp[i]));
      }
      return acc;
    };

    ext_storage_type gamma_tilde = ext_field_->AdditiveIdentity();
    bool found_root = false;
    ext_storage_type cur_eta = eta;
    for (uint64_t e = 1; e < base_order_; ++e) {
      if (std::gcd(e, base_order_ - 1) == 1 &&
          ext_field_->IsZero(eval_min_poly_in_ext(cur_eta))) {
        gamma_tilde = cur_eta;
        found_root = true;
        break;
      }
      cur_eta = ext_field_->Mul(cur_eta, eta);
    }
    if (!found_root) {
      throw std::runtime_error(
          "FieldTowerAdapter: failed to embed base field into extension field");
    }

    // Precompute powers 1, gamma_b, ..., gamma_b^{s-1} and express canonical
    // F_p-basis of BaseField in terms of them via s x s matrix inversion.
    std::vector<uint64_t> pow_mat(s * s, 0);
    base_storage_type b_pow = base_field_->MultiplicativeIdentity();
    for (std::size_t col = 0; col < s; ++col) {
      const auto coords = detail::PrimeFieldCoords(*base_field_, b_pow);
      for (std::size_t r = 0; r < s; ++r) {
        pow_mat[r * s + col] = coords[r];
      }
      b_pow = base_field_->Mul(b_pow, gamma_b);
    }
    auto inv_pow = detail::InvertModP(pow_mat, s, p);
    if (!inv_pow.has_value()) {
      throw std::runtime_error(
          "FieldTowerAdapter: singular power basis matrix in base field");
    }

    std::vector<ext_storage_type> gamma_tilde_pows(s);
    gamma_tilde_pows[0] = ext_field_->MultiplicativeIdentity();
    for (std::size_t i = 1; i < s; ++i) {
      gamma_tilde_pows[i] =
          ext_field_->Mul(gamma_tilde_pows[i - 1], gamma_tilde);
    }

    // Images phi(e_0), ..., phi(e_{s-1}) of the canonical F_p-basis of
    // BaseField.
    std::vector<ext_storage_type> phi_basis(s);
    for (std::size_t u = 0; u < s; ++u) {
      ext_storage_type img = ext_field_->AdditiveIdentity();
      for (std::size_t j = 0; j < s; ++j) {
        const uint64_t c = (*inv_pow)[j * s + u];
        if (c != 0) {
          img = ext_field_->Add(
              img, ext_field_->Mul(fp_to_ext(c), gamma_tilde_pows[j]));
        }
      }
      phi_basis[u] = img;
    }

    // Step 2: Select an F_q-basis {b_0, ..., b_{m-1}} of ExtField starting with
    // powers 1, gamma_e, gamma_e^2, ... (or canonical F_p basis fallback).
    std::vector<ext_storage_type> q_basis;
    q_basis.reserve(m);
    std::vector<uint64_t>
        chosen_cols;  // Stored column-major: (k * s) cols of height M.

    auto try_add_basis = [&](const ext_storage_type &cand) -> bool {
      const std::size_t prev_cols = q_basis.size() * s;
      const std::size_t new_cols = prev_cols + s;
      std::vector<uint64_t> test_mat(M * new_cols, 0);
      for (std::size_t r = 0; r < M; ++r) {
        for (std::size_t c = 0; c < prev_cols; ++c) {
          test_mat[r * new_cols + c] = chosen_cols[c * M + r];
        }
      }
      std::vector<std::vector<uint64_t>> cand_cols(s);
      for (std::size_t u = 0; u < s; ++u) {
        const ext_storage_type prod = ext_field_->Mul(phi_basis[u], cand);
        cand_cols[u] = detail::PrimeFieldCoords(*ext_field_, prod);
        for (std::size_t r = 0; r < M; ++r) {
          test_mat[r * new_cols + prev_cols + u] = cand_cols[u][r];
        }
      }
      if (detail::RankModP(test_mat, M, new_cols, p) == new_cols) {
        q_basis.push_back(cand);
        for (std::size_t u = 0; u < s; ++u) {
          chosen_cols.insert(chosen_cols.end(), cand_cols[u].begin(),
                             cand_cols[u].end());
        }
        return true;
      }
      return false;
    };

    // First try powers of gamma_e: 1, gamma_e, gamma_e^2, ..., gamma_e^{M-1}.
    ext_storage_type e_pow = ext_field_->MultiplicativeIdentity();
    for (std::size_t i = 0; i < M && q_basis.size() < m; ++i) {
      try_add_basis(e_pow);
      e_pow = ext_field_->Mul(e_pow, gamma_e);
    }
    // Fallback to canonical F_p-basis vectors of ExtField if needed.
    for (std::size_t i = 0; i < M && q_basis.size() < m; ++i) {
      std::vector<uint64_t> unit(M, 0);
      unit[i] = 1;
      try_add_basis(detail::FromPrimeFieldCoords(*ext_field_, unit));
    }
    if (q_basis.size() != m) {
      throw std::runtime_error(
          "FieldTowerAdapter: failed to construct relative F_q-basis");
    }

    // Build M x M matrix T (row-major) and its inverse T^{-1} over F_p.
    tower_contract_mat_.assign(M * M, 0);
    for (std::size_t col = 0; col < M; ++col) {
      for (std::size_t r = 0; r < M; ++r) {
        tower_contract_mat_[r * M + col] = chosen_cols[col * M + r];
      }
    }
    auto inv_t = detail::InvertModP(tower_contract_mat_, M, p);
    if (!inv_t.has_value()) {
      throw std::runtime_error(
          "FieldTowerAdapter: change-of-basis matrix is singular");
    }
    tower_expand_mat_ = std::move(*inv_t);
  }

  void ExpandElementInto(const ext_storage_type &elem, base_storage_type *dst,
                         std::size_t stride) const {
    if (!use_general_tower_) {
      if constexpr (std::is_same_v<BaseField, ExtField>) {
        if (relative_degree_ == 1) {
          dst[0] = elem;
          return;
        }
      }
      if constexpr (detail::IsDirectRelativeTower<BaseField, ExtField>) {
        const auto &coeffs = elem.Coeffs();
        const base_storage_type zero = base_field_->AdditiveIdentity();
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          dst[i * stride] = (i < coeffs.size()) ? coeffs[i] : zero;
        }
        return;
      } else if constexpr (detail::HasBit128Interface<ext_storage_type> &&
                           std::is_integral_v<base_storage_type>) {
        const ext_storage_type val = ext_field_->GetElementValue(elem);
        const base_storage_type zero = base_field_->AdditiveIdentity();
        const base_storage_type one = base_field_->MultiplicativeIdentity();
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          dst[i * stride] = val.Bit(static_cast<int>(i)) ? one : zero;
        }
        return;
      } else if constexpr (detail::HasExtensionCoeffsInterface<
                               ExtField, ext_storage_type> &&
                           std::is_integral_v<base_storage_type>) {
        const ext_storage_type canonical = ext_field_->GetElementValue(elem);
        const auto &coeffs = canonical.Coeffs();
        const base_storage_type zero = base_field_->AdditiveIdentity();
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          dst[i * stride] = (i < coeffs.size())
                                ? base_field_->SetElementValue(
                                      static_cast<base_storage_type>(coeffs[i]))
                                : zero;
        }
        return;
      } else if constexpr (std::is_integral_v<base_storage_type> &&
                           std::is_integral_v<ext_storage_type>) {
        const uint64_t val =
            static_cast<uint64_t>(ext_field_->GetElementValue(elem));
        if (relative_degree_ == 1) {
          dst[0] =
              base_field_->SetElementValue(static_cast<base_storage_type>(val));
        } else {
          // Binary extension F_{2^m} / F_2: bit-slice polynomial-basis value.
          const base_storage_type zero = base_field_->AdditiveIdentity();
          const base_storage_type one = base_field_->MultiplicativeIdentity();
          for (std::size_t i = 0; i < relative_degree_; ++i) {
            dst[i * stride] = (((val >> i) & 1ULL) != 0) ? one : zero;
          }
        }
        return;
      }
    }

    // General relative tower path via T^{-1} over F_p.
    const auto ext_coords = detail::PrimeFieldCoords(*ext_field_, elem);
    const auto base_fp =
        detail::MatVecModP(tower_expand_mat_, ext_degree_, ext_degree_,
                           ext_coords, characteristic_);
    for (std::size_t r = 0; r < relative_degree_; ++r) {
      dst[r * stride] =
          detail::FromPrimeFieldCoords(*base_field_, base_fp, r * base_degree_);
    }
  }

  ext_storage_type ContractElementFrom(const base_storage_type *src,
                                       std::size_t stride) const {
    if (!use_general_tower_) {
      if constexpr (std::is_same_v<BaseField, ExtField>) {
        if (relative_degree_ == 1) {
          return src[0];
        }
      }
      if constexpr (detail::IsDirectRelativeTower<BaseField, ExtField>) {
        std::vector<base_storage_type> coeffs(relative_degree_);
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          coeffs[i] = src[i * stride];
        }
        return ext_storage_type(std::move(coeffs));
      } else if constexpr (detail::HasBit128Interface<ext_storage_type> &&
                           std::is_integral_v<base_storage_type>) {
        ext_storage_type val{};
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          if (!base_field_->IsZero(src[i * stride])) {
            val |= (ext_storage_type{1, 0} << static_cast<int>(i));
          }
        }
        return ext_field_->SetElementValue(val);
      } else if constexpr (detail::HasExtensionCoeffsInterface<
                               ExtField, ext_storage_type> &&
                           std::is_integral_v<base_storage_type>) {
        using CoeffT = typename ext_storage_type::coeff_type;
        std::vector<CoeffT> coeffs(relative_degree_);
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          coeffs[i] = static_cast<CoeffT>(
              base_field_->GetElementValue(src[i * stride]));
        }
        return ext_field_->FromCoeffs(std::move(coeffs));
      } else if constexpr (std::is_integral_v<base_storage_type> &&
                           std::is_integral_v<ext_storage_type>) {
        if (relative_degree_ == 1) {
          return ext_field_->SetElementValue(static_cast<ext_storage_type>(
              base_field_->GetElementValue(src[0])));
        }
        uint64_t bits = 0;
        for (std::size_t i = 0; i < relative_degree_; ++i) {
          if (!base_field_->IsZero(src[i * stride])) {
            bits |= (uint64_t{1} << i);
          }
        }
        return ext_field_->SetElementValue(static_cast<ext_storage_type>(bits));
      }
    }

    // General relative tower path via T over F_p.
    std::vector<uint64_t> base_fp(ext_degree_, 0);
    for (std::size_t r = 0; r < relative_degree_; ++r) {
      const auto c = detail::PrimeFieldCoords(*base_field_, src[r * stride]);
      for (std::size_t u = 0; u < base_degree_; ++u) {
        base_fp[r * base_degree_ + u] = c[u];
      }
    }
    const auto ext_coords =
        detail::MatVecModP(tower_contract_mat_, ext_degree_, ext_degree_,
                           base_fp, characteristic_);
    return detail::FromPrimeFieldCoords(*ext_field_, ext_coords);
  }

  std::shared_ptr<BaseField> base_field_;
  std::shared_ptr<ExtField> ext_field_;
  uint64_t characteristic_ = 0;
  std::size_t base_degree_ = 1;
  std::size_t ext_degree_ = 1;
  std::size_t relative_degree_ = 1;
  uint64_t base_order_ = 0;

  bool use_general_tower_ = false;
  std::vector<uint64_t> tower_contract_mat_;
  std::vector<uint64_t> tower_expand_mat_;
};

// ============================================================================
// Linearized Polynomial Operator Matrix & Exact Kernel Basis over F_q
// ============================================================================

// Builds the m x m matrix M_lambda over BaseField representing the F_q-linear
// operator x -> EvalOperator(lambda_poly, x) in the canonical tower basis.
template <typename BaseField, typename ExtField>
GFMatrix<BaseField> BuildOperatorMatrix(
    const xg::ring::FrobeniusSkewRing<ExtField> &skew_ring,
    const FieldTowerAdapter<BaseField, ExtField> &tower,
    const typename xg::ring::FrobeniusSkewRing<ExtField>::element_type
        &lambda_poly) {
  const std::size_t m = tower.RelativeDegree();
  using RawArray = typename GFArray<ExtField>::raw_array_type;
  using ShapeType = typename GFArray<ExtField>::shape_type;

  RawArray img_raw = RawArray::from_shape(ShapeType{m});
  for (std::size_t c = 0; c < m; ++c) {
    const auto basis_elem = tower.CanonicalBasisElement(c);
    img_raw(c) = skew_ring.EvalOperator(lambda_poly, basis_elem);
  }
  const GFArray<ExtField> images =
      GFArray<ExtField>::FromRaw(tower.ExtFieldPtr(), std::move(img_raw));
  return tower.ExpandVector(images);
}

// Computes an exact F_q-basis of {x in F_{q^m} : lambda_poly(x) == 0} via
// m x m linalg::Nullspace over BaseField.
template <typename BaseField, typename ExtField>
GFArray<ExtField> KernelBasis(
    const xg::ring::FrobeniusSkewRing<ExtField> &skew_ring,
    const FieldTowerAdapter<BaseField, ExtField> &tower,
    const typename xg::ring::FrobeniusSkewRing<ExtField>::element_type
        &lambda_poly) {
  const GFMatrix<BaseField> op_mat =
      BuildOperatorMatrix(skew_ring, tower, lambda_poly);
  const GFMatrix<BaseField> null_rows = linalg::Nullspace(op_mat);
  if (null_rows.Shape(0) == 0) {
    return GFArray<ExtField>::Zeros(tower.ExtFieldPtr(), {0});
  }
  return tower.ContractMatrix(linalg::Transpose(null_rows));
}

// Forward declarations of decoders.
template <typename BaseField, typename ExtField>
class GabidulinSkewEEADecoder;
template <typename BaseField, typename ExtField>
class GabidulinLoidreauDecoder;

// ============================================================================
// GabidulinCode: Maximum Rank Distance (MRD) Code G(n, k) over F_{q^m} / F_q
// ============================================================================

template <typename BaseField, typename ExtField>
class GabidulinCode : public AbstractLinearCode<ExtField> {
 public:
  using Base = AbstractLinearCode<ExtField>;
  using base_field_type = BaseField;
  using ext_field_type = ExtField;
  using field_type = ExtField;
  using element_type = GaloisFieldElement<ExtField>;
  using storage_type = typename ExtField::element_type;
  using base_storage_type = typename BaseField::element_type;
  using codeword_type = GFArray<ExtField>;
  using message_type = GFArray<ExtField>;
  using vector_type = GFArray<ExtField>;
  using matrix_type = GFMatrix<ExtField>;
  using base_matrix_type = GFMatrix<BaseField>;
  using result_type = DecodingResult<ExtField>;
  using tower_type = FieldTowerAdapter<BaseField, ExtField>;
  using skew_ring_type = xg::ring::FrobeniusSkewRing<ExtField>;

  GabidulinCode(std::shared_ptr<BaseField> base_field,
                std::shared_ptr<ExtField> ext_field,
                const vector_type &evaluation_points, std::size_t dimension)
      : Base(ext_field, evaluation_points.Size(), dimension, Metric::RANK,
             "GeneratorMatrix", "GabidulinSkewEEA"),
        base_field_(std::move(base_field)),
        tower_(base_field_, ext_field),
        skew_ring_(std::make_shared<skew_ring_type>(
            *ext_field, static_cast<uint32_t>(base_field_->Degree()))),
        eval_points_(evaluation_points) {
    if (!eval_points_.HasField() || eval_points_.Field() != this->field_ ||
        eval_points_.Dimension() != 1) {
      throw std::invalid_argument(
          "GabidulinCode: evaluation points must be a 1D GFArray over "
          "ExtField");
    }
    const std::size_t n = this->length_;
    const std::size_t k = this->dimension_;
    const std::size_t m = tower_.RelativeDegree();
    if (k == 0 || k > n || n > m) {
      throw std::invalid_argument(
          "GabidulinCode: parameters must satisfy 1 <= k <= n <= m");
    }
    if (tower_.RankWeight(eval_points_) != n) {
      throw std::invalid_argument(
          "GabidulinCode: evaluation points must be F_q-linearly independent");
    }

    BuildPrecomputedData();
    RegisterEncodersAndDecoders();
  }

  GabidulinCode(const GabidulinCode &) = default;
  GabidulinCode &operator=(const GabidulinCode &) = default;
  GabidulinCode(GabidulinCode &&) noexcept = default;
  GabidulinCode &operator=(GabidulinCode &&) noexcept = default;
  ~GabidulinCode() override = default;

  // Standard Gabidulin code using the polynomial basis 1, alpha, ...,
  // alpha^{n-1}.
  static GabidulinCode CreateStandard(std::shared_ptr<BaseField> base_field,
                                      std::shared_ptr<ExtField> ext_field,
                                      std::size_t length,
                                      std::size_t dimension) {
    if (!base_field || !ext_field) {
      throw std::invalid_argument(
          "GabidulinCode::CreateStandard: fields cannot be null");
    }
    tower_type tower(base_field, ext_field);
    if (length == 0 || dimension == 0 || dimension > length ||
        length > tower.RelativeDegree()) {
      throw std::invalid_argument(
          "GabidulinCode::CreateStandard: parameters must satisfy "
          "1 <= k <= n <= m");
    }
    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    RawArray g_raw = RawArray::from_shape(ShapeType{length});
    for (std::size_t j = 0; j < length; ++j) {
      g_raw(j) = tower.CanonicalBasisElement(j);
    }
    vector_type eval_pts = vector_type::FromRaw(ext_field, std::move(g_raw));
    return GabidulinCode(std::move(base_field), std::move(ext_field), eval_pts,
                         dimension);
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  // Exact MRD minimum rank distance: d_R = n - k + 1.
  std::size_t MinimumDistance() const override {
    return this->length_ - this->dimension_ + 1;
  }

  // Exact F_q-rank weight of a length-n vector over F_{q^m}.
  std::size_t Weight(const codeword_type &word) const override {
    if (!word.HasField() || word.Field() != this->field_) {
      throw std::invalid_argument("GabidulinCode::Weight: field mismatch");
    }
    if (word.Dimension() != 1 || word.Size() != this->length_) {
      throw std::invalid_argument(
          "GabidulinCode::Weight: word must be 1D with size equal to code "
          "length");
    }
    return tower_.RankWeight(word);
  }

  bool IsMRD() const noexcept { return true; }

  std::shared_ptr<BaseField> BaseFieldPtr() const { return base_field_; }
  const tower_type &Tower() const noexcept { return tower_; }
  const skew_ring_type &SkewRing() const noexcept { return *skew_ring_; }

  const vector_type &EvaluationPoints() const noexcept { return eval_points_; }
  const vector_type &DualEvaluationPoints() const noexcept {
    return dual_eval_points_;
  }
  const base_matrix_type &ExpandedDualEvaluationPoints() const noexcept {
    return expanded_dual_points_;
  }

  base_matrix_type ExpandToBaseMatrix(const codeword_type &word) const {
    return tower_.ExpandVector(word);
  }

  codeword_type ContractFromBaseMatrix(const base_matrix_type &mat) const {
    return tower_.ContractMatrix(mat);
  }

  GabidulinCode DualCode() const {
    if (this->dimension_ == this->length_) {
      throw std::invalid_argument(
          "GabidulinCode::DualCode: dual of full-dimension code has dimension "
          "0");
    }
    return GabidulinCode(base_field_, this->field_, dual_eval_points_,
                         this->length_ - this->dimension_);
  }

  // Silva-Kotter-Kschischang (2008) errors + row-erasures + column-erasures
  // decoder correcting 2*tau + rho + gamma <= n - k.
  result_type DecodeWithRankErasures(
      const codeword_type &received_word, const vector_type &row_erasure_values,
      const base_matrix_type &col_erasure_matrix) const;

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + ", " +
           std::to_string(MinimumDistance()) + "] Gabidulin MRD code over GF(" +
           std::to_string(tower_.Characteristic()) + "^" +
           std::to_string(this->field_->Degree()) + ")/GF(" +
           std::to_string(tower_.Characteristic()) + "^" +
           std::to_string(tower_.BaseDegree()) + ")";
  }

 private:
  void BuildPrecomputedData() {
    const std::size_t n = this->length_;
    const std::size_t k = this->dimension_;
    const std::size_t r = n - k;
    const auto g_span = eval_points_.RawSpan();

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    // 1. Moore generator matrix: G_{i, j} = g_j^{[i]} for 0 <= i < k, 0 <= j <
    // n.
    RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
    storage_type *g_data = g_raw.data();
    for (std::size_t j = 0; j < n; ++j) {
      storage_type cur = g_span[j];
      for (std::size_t i = 0; i < k; ++i) {
        g_data[i * n + j] = cur;
        if (i + 1 < k) {
          cur = tower_.FrobeniusQ(cur, 1);
        }
      }
    }
    generator_matrix_ = matrix_type::FromRaw(this->field_, std::move(g_raw));

    // 2. Dual vector h in ker([g_j^{[e]}]_{0 <= e <= n - 2}) \ {0}.
    vector_type h_base;
    if (n == 1) {
      h_base = vector_type::Ones(this->field_, {1});
    } else {
      RawArray m_raw = RawArray::from_shape(ShapeType{n - 1, n});
      storage_type *m_data = m_raw.data();
      for (std::size_t j = 0; j < n; ++j) {
        storage_type cur = g_span[j];
        for (std::size_t e = 0; e < n - 1; ++e) {
          m_data[e * n + j] = cur;
          if (e + 1 < n - 1) {
            cur = tower_.FrobeniusQ(cur, 1);
          }
        }
      }
      const matrix_type moore_n_minus_1 =
          matrix_type::FromRaw(this->field_, std::move(m_raw));
      const matrix_type ns = linalg::Nullspace(moore_n_minus_1);
      if (ns.Shape(0) != 1 || ns.Shape(1) != n) {
        throw std::runtime_error(
            "GabidulinCode: expected 1-dimensional dual nullspace");
      }
      RawArray h_raw = RawArray::from_shape(ShapeType{n});
      for (std::size_t j = 0; j < n; ++j) {
        h_raw(j) = ns.Raw(0, j);
      }
      h_base = vector_type::FromRaw(this->field_, std::move(h_raw));
    }

    // Shift h'_j = theta^{-(n - k - 1)}(h_j) so that H_{i,j} = (h'_j)^{[i]}
    // satisfies G * H^T == 0 identically.
    const std::size_t shift = (r > 0) ? (r - 1) : 0;
    RawArray hp_raw = RawArray::from_shape(ShapeType{n});
    for (std::size_t j = 0; j < n; ++j) {
      hp_raw(j) = tower_.InvFrobeniusQ(h_base.Raw(j), shift);
    }
    dual_eval_points_ = vector_type::FromRaw(this->field_, std::move(hp_raw));
    expanded_dual_points_ = tower_.ExpandVector(dual_eval_points_);

    // 3. Moore parity-check matrix: H_{i, j} = (h'_j)^{[i]} for 0 <= i < r.
    if (r == 0) {
      parity_check_matrix_ = matrix_type::Zeros(this->field_, {0, n});
    } else {
      RawArray h_raw = RawArray::from_shape(ShapeType{r, n});
      storage_type *h_data = h_raw.data();
      const auto hp_span = dual_eval_points_.RawSpan();
      for (std::size_t j = 0; j < n; ++j) {
        storage_type cur = hp_span[j];
        for (std::size_t i = 0; i < r; ++i) {
          h_data[i * n + j] = cur;
          if (i + 1 < r) {
            cur = tower_.FrobeniusQ(cur, 1);
          }
        }
      }
      parity_check_matrix_ =
          matrix_type::FromRaw(this->field_, std::move(h_raw));
    }
  }

  void RegisterEncodersAndDecoders();

  std::shared_ptr<BaseField> base_field_;
  tower_type tower_;
  std::shared_ptr<skew_ring_type> skew_ring_;
  vector_type eval_points_;
  vector_type dual_eval_points_;
  base_matrix_type expanded_dual_points_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;
};

// ============================================================================
// GabidulinSkewEEADecoder: Partial Skew-EEA + SKK Rank Errors/Erasures Decoder
// ============================================================================

template <typename BaseField, typename ExtField>
class GabidulinSkewEEADecoder : public Decoder<ExtField> {
 public:
  using Base = Decoder<ExtField>;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using storage_type = typename ExtField::element_type;
  using base_storage_type = typename BaseField::element_type;
  using code_type = GabidulinCode<BaseField, ExtField>;
  using skew_ring_type = xg::ring::FrobeniusSkewRing<ExtField>;
  using skew_poly_type = typename skew_ring_type::element_type;

  explicit GabidulinSkewEEADecoder(const AbstractCode<ExtField> &code)
      : Base(code) {
    gab_code_ = dynamic_cast<const code_type *>(&code);
    if (!gab_code_) {
      throw std::invalid_argument(
          "GabidulinSkewEEADecoder requires a GabidulinCode instance");
    }
  }

  std::string Name() const override { return "GabidulinSkewEEA"; }

  std::size_t DecodingRadius() const override {
    return (gab_code_->Length() - gab_code_->Dimension()) / 2;
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto empty_rows = GFArray<ExtField>::Zeros(gab_code_->Field(), {0});
    auto empty_cols = GFMatrix<BaseField>::Zeros(gab_code_->BaseFieldPtr(),
                                                 {0, gab_code_->Length()});
    return DecodeWithRankErasures(received_word, empty_rows, empty_cols);
  }

  // Standard coordinate-erasure decoding: each erased coordinate j becomes a
  // unit column-erasure vector e_j in F_q^n (corrects 2*tau + gamma <= n - k).
  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const override {
    const std::size_t n = gab_code_->Length();
    auto base_field = gab_code_->BaseFieldPtr();
    std::vector<std::size_t> unique_pos = erasure_positions;
    std::sort(unique_pos.begin(), unique_pos.end());
    unique_pos.erase(std::unique(unique_pos.begin(), unique_pos.end()),
                     unique_pos.end());
    for (std::size_t pos : unique_pos) {
      if (pos >= n) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: erasure position out of range");
      }
    }
    auto col_mat =
        GFMatrix<BaseField>::Zeros(base_field, {unique_pos.size(), n});
    const base_storage_type base_one = base_field->MultiplicativeIdentity();
    for (std::size_t r = 0; r < unique_pos.size(); ++r) {
      col_mat.Raw(r, unique_pos[r]) = base_one;
    }
    auto empty_rows = GFArray<ExtField>::Zeros(gab_code_->Field(), {0});
    return DecodeWithRankErasures(received_word, empty_rows, col_mat);
  }

  // Full Silva-Kotter-Kschischang (2008) rank errors + row erasures + column
  // erasures decoder correcting 2*tau + rho + gamma <= n - k.
  result_type DecodeWithRankErasures(
      const codeword_type &received_word,
      const GFArray<ExtField> &row_erasure_values,
      const GFMatrix<BaseField> &col_erasure_matrix) const {
    auto ext_field = gab_code_->Field();
    auto base_field = gab_code_->BaseFieldPtr();
    const auto &tower = gab_code_->Tower();
    const auto &skew_ring = gab_code_->SkewRing();
    const std::size_t n = gab_code_->Length();
    const std::size_t k = gab_code_->Dimension();
    const std::size_t r = n - k;

    if (!received_word.HasField() || received_word.Field() != ext_field ||
        received_word.Dimension() != 1 || received_word.Size() != n) {
      return result_type::Failure(
          "GabidulinSkewEEADecoder: invalid received word");
    }
    if (!row_erasure_values.HasField() ||
        row_erasure_values.Field() != ext_field ||
        row_erasure_values.Dimension() != 1) {
      return result_type::Failure(
          "GabidulinSkewEEADecoder: row_erasure_values must be a 1D array over "
          "ExtField");
    }
    if (!col_erasure_matrix.HasField() ||
        col_erasure_matrix.Field() != base_field ||
        col_erasure_matrix.Dimension() != 2 ||
        col_erasure_matrix.Shape(1) != n) {
      return result_type::Failure(
          "GabidulinSkewEEADecoder: col_erasure_matrix must be a 2D matrix "
          "over BaseField with n columns");
    }

    // 1. Reduce row erasure values to an F_q-linearly independent basis
    // a^{(R)}.
    std::vector<storage_type> a_row_basis;
    if (row_erasure_values.Size() > 0) {
      const auto a_exp = tower.ExpandVector(row_erasure_values);
      const auto basis_rows = linalg::RowSpace(linalg::Transpose(a_exp));
      if (basis_rows.Shape(0) > 0) {
        const auto basis_vec =
            tower.ContractMatrix(linalg::Transpose(basis_rows));
        a_row_basis.assign(basis_vec.RawSpan().begin(),
                           basis_vec.RawSpan().end());
      }
    }
    const std::size_t rho = a_row_basis.size();

    // 2. Reduce column erasure matrix to full row rank B^{(C)} (gamma x n).
    GFMatrix<BaseField> b_col_basis =
        (col_erasure_matrix.Shape(0) > 0)
            ? linalg::RowSpace(col_erasure_matrix)
            : GFMatrix<BaseField>::Zeros(base_field, {0, n});
    const std::size_t gamma = b_col_basis.Shape(0);

    if (rho + gamma > r) {
      return result_type::Failure(
          "GabidulinSkewEEADecoder: erasures exceed redundancy n - k");
    }

    // 3. Compute raw syndrome s = H * r^T (length r).
    const auto syn_arr = gab_code_->Syndrome(received_word);
    skew_poly_type s_poly(syn_arr.RawSpan().begin(), syn_arr.RawSpan().end());
    s_poly = skew_ring.Normalize(std::move(s_poly));

    if (skew_ring.IsZero(s_poly)) {
      return result_type::Success(
          received_word, gab_code_->Unencode(received_word), 0, rho + gamma);
    }

    const storage_type ext_zero = ext_field->AdditiveIdentity();

    // 4. Compute row-erasure annihilator Lambda_R(x) =
    // MinimalPolynomial(a^{(R)}).
    const skew_poly_type lambda_r = skew_ring.MinimalPolynomial(a_row_basis);

    // Compute S_1 = Lambda_R * S. For rho <= j < r, row erasures are
    // annihilated.
    const skew_poly_type s1_poly = skew_ring.Mul(lambda_r, s_poly);

    // 5. Compute known column-erasure locators y^{(C)}_v = sum_j B^{(C)}_{v,j}
    // * h'_j
    //    via F_q-linear expansion: Y^{(C)}_exp = H'_exp * (B^{(C)})^T.
    std::vector<storage_type> y_col(gamma);
    if (gamma > 0) {
      const auto y_col_exp =
          linalg::Dot(gab_code_->ExpandedDualEvaluationPoints(),
                      linalg::Transpose(b_col_basis));
      const auto y_col_vec = tower.ContractMatrix(y_col_exp);
      y_col.assign(y_col_vec.RawSpan().begin(), y_col_vec.RawSpan().end());
    }

    // Compute column-erasure annihilator Gamma_C(x) =
    // MinimalPolynomial(y^{(C)}).
    const skew_poly_type gamma_c = skew_ring.MinimalPolynomial(y_col);

    // 6. Construct modified pure-error syndrome S*(x) of length r* = r - rho -
    // gamma:
    //    s*_{rho + i} = sum_{d=0}^{gamma} [S_1]_{rho + i + d} *
    //    (Gamma_C)_d^{[rho + i]}.
    const std::size_t r_star = r - rho - gamma;
    const std::size_t tau_max = r_star / 2;
    skew_poly_type s_star(r_star, ext_zero);
    for (std::size_t i = 0; i < r_star; ++i) {
      const std::size_t j = rho + i;
      storage_type acc = ext_zero;
      for (std::size_t d = 0; d <= gamma; ++d) {
        const std::size_t idx = j + d;
        if (idx < s1_poly.size() && !ext_field->IsZero(s1_poly[idx]) &&
            d < gamma_c.size() && !ext_field->IsZero(gamma_c[d])) {
          const storage_type g_shifted = tower.FrobeniusQ(gamma_c[d], j);
          acc = ext_field->Add(acc, ext_field->Mul(s1_poly[idx], g_shifted));
        }
      }
      s_star[i] = acc;
    }
    s_star = skew_ring.Normalize(std::move(s_star));

    // 7. Solve the pure rank-error key equation for Lambda_E(x) via partial
    //    RightExtendedGcd(x^{r*}, S*, stop_deg = r* - tau_max - 1).
    skew_poly_type lambda_e = skew_ring.MultiplicativeIdentity();
    std::size_t tau_eff = 0;
    if (!skew_ring.IsZero(s_star)) {
      if (tau_max == 0) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: uncorrectable rank errors exceed budget");
      }
      const skew_poly_type x_r_star =
          skew_ring.Monomial(ext_field->MultiplicativeIdentity(), r_star);
      const int64_t stop_deg = static_cast<int64_t>(r_star - tau_max - 1);
      auto [omega, u_poly, v_poly] =
          skew_ring.RightExtendedGcd(x_r_star, s_star, stop_deg);
      (void)omega;
      (void)u_poly;
      const int64_t deg_e = skew_ring.Degree(v_poly);
      if (deg_e < 1 || static_cast<std::size_t>(deg_e) > tau_max) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: key equation yielded invalid error "
            "degree");
      }
      tau_eff = static_cast<std::size_t>(deg_e);
      const storage_type inv_lc =
          ext_field->Inv(skew_ring.LeadingCoeff(v_poly));
      lambda_e = skew_ring.LeftScalarMul(inv_lc, v_poly);
    }

    // 8. Combined error-span annihilator Lambda = Lambda_E * Lambda_R
    //    (degree nu = rho + tau_eff).
    const skew_poly_type lambda = skew_ring.Mul(lambda_e, lambda_r);
    const std::size_t nu = rho + tau_eff;
    const GFMatrix<BaseField> op_lambda =
        BuildOperatorMatrix(skew_ring, tower, lambda);
    const GFMatrix<BaseField> null_rows = linalg::Nullspace(op_lambda);
    if (null_rows.Shape(0) != nu) {
      return result_type::Failure(
          "GabidulinSkewEEADecoder: error-span polynomial does not split "
          "completely in F_{q^m}");
    }
    const GFArray<ExtField> a_tot =
        (nu > 0) ? tower.ContractMatrix(linalg::Transpose(null_rows))
                 : GFArray<ExtField>::Zeros(ext_field, {0});

    // 9. Recover particular values tilde_a^{(C)} for the gamma column erasures.
    //    Let T = Lambda * S. For nu <= j < nu + gamma:
    //    [T]_j = sum_{v=0}^{gamma-1} d_v * (y^{(C)}_v)^{[j]}, where d_v =
    //    Lambda(a^{(C)}_v).
    codeword_type err_col = codeword_type::Zeros(ext_field, {n});
    std::vector<storage_type> a_col_part(gamma, ext_zero);

    if (gamma > 0) {
      const skew_poly_type t_poly = skew_ring.Mul(lambda, s_poly);

      using ExtRaw = typename GFMatrix<ExtField>::raw_array_type;
      using ExtShape = typename GFMatrix<ExtField>::shape_type;
      ExtRaw mat_raw = ExtRaw::from_shape(ExtShape{gamma, gamma});
      ExtRaw rhs_raw = ExtRaw::from_shape(ExtShape{gamma});

      for (std::size_t row = 0; row < gamma; ++row) {
        const std::size_t j = nu + row;
        rhs_raw(row) = (j < t_poly.size()) ? t_poly[j] : ext_zero;
        for (std::size_t v = 0; v < gamma; ++v) {
          mat_raw(row, v) = tower.FrobeniusQ(y_col[v], j);
        }
      }

      GFArray<ExtField> d_vec;
      try {
        d_vec = linalg::Solve(
            GFMatrix<ExtField>::FromRaw(ext_field, std::move(mat_raw)),
            GFArray<ExtField>::FromRaw(ext_field, std::move(rhs_raw)));
      } catch (const std::runtime_error &) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: singular Moore system for column "
            "erasures");
      }

      // Find any particular solution tilde_a^{(C)}_v to Lambda(x) = d_v over
      // F_q. Any kernel component in ker(Lambda) = span(a_tot) is absorbed into
      // B_tot.
      const GFMatrix<BaseField> d_exp = tower.ExpandVector(d_vec);
      auto a_col_exp_opt = detail::SolveLinearSystem(op_lambda, d_exp);
      if (!a_col_exp_opt.has_value()) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: column-erasure values not in image of "
            "Lambda");
      }
      const GFArray<ExtField> a_col_vec = tower.ContractMatrix(*a_col_exp_opt);
      a_col_part.assign(a_col_vec.RawSpan().begin(), a_col_vec.RawSpan().end());

      // Reconstruct e^{(C)} = tilde_a^{(C)} * B^{(C)}.
      const GFMatrix<BaseField> e_col_exp =
          linalg::Dot(*a_col_exp_opt, b_col_basis);
      err_col = tower.ContractMatrix(e_col_exp);
    }

    // 10. Subtract column-erasure syndrome and solve the nu x nu
    // Frobenius-Moore
    //     system for tilde_y_tot:
    //     tilde_s_i^{[-i]} = sum_{u=0}^{nu-1} (a_tot_u)^{[-i]} * tilde_y_tot_u.
    codeword_type err_tot = err_col;
    if (nu > 0) {
      using ExtRaw = typename GFMatrix<ExtField>::raw_array_type;
      using ExtShape = typename GFMatrix<ExtField>::shape_type;
      ExtRaw moore_raw = ExtRaw::from_shape(ExtShape{nu, nu});
      ExtRaw rhs_raw = ExtRaw::from_shape(ExtShape{nu});

      for (std::size_t i = 0; i < nu; ++i) {
        storage_type s_rem = syn_arr.Raw(i);
        for (std::size_t v = 0; v < gamma; ++v) {
          if (!ext_field->IsZero(a_col_part[v])) {
            const storage_type term =
                ext_field->Mul(a_col_part[v], tower.FrobeniusQ(y_col[v], i));
            s_rem = ext_field->Sub(s_rem, term);
          }
        }
        rhs_raw(i) = tower.InvFrobeniusQ(s_rem, i);
        for (std::size_t u = 0; u < nu; ++u) {
          moore_raw(i, u) = tower.InvFrobeniusQ(a_tot.Raw(u), i);
        }
      }

      GFArray<ExtField> y_tot;
      try {
        y_tot = linalg::Solve(
            GFMatrix<ExtField>::FromRaw(ext_field, std::move(moore_raw)),
            GFArray<ExtField>::FromRaw(ext_field, std::move(rhs_raw)));
      } catch (const std::runtime_error &) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: singular Frobenius-Moore locator system");
      }

      // 11. Expand y_tot = sum_j B_tot_{u,j} * h'_j over F_q:
      //     H'_exp * (B_tot)^T = ExpandVector(y_tot).
      const GFMatrix<BaseField> y_tot_exp = tower.ExpandVector(y_tot);
      auto b_tot_t_opt = detail::SolveLinearSystem(
          gab_code_->ExpandedDualEvaluationPoints(), y_tot_exp);
      if (!b_tot_t_opt.has_value()) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: locators outside span of dual evaluation "
            "points");
      }
      const GFMatrix<BaseField> b_tot = linalg::Transpose(*b_tot_t_opt);
      const GFMatrix<BaseField> a_tot_exp = tower.ExpandVector(a_tot);
      const GFMatrix<BaseField> e_re_exp = linalg::Dot(a_tot_exp, b_tot);
      const codeword_type err_re = tower.ContractMatrix(e_re_exp);
      err_tot = err_tot + err_re;
    }

    // 12. Verify candidate codeword and rank-weight constraints.
    const codeword_type candidate = received_word - err_tot;
    if (!gab_code_->Contains(candidate)) {
      return result_type::Failure(
          "GabidulinSkewEEADecoder: reconstructed word failed parity-check "
          "verification");
    }
    if (rho == 0 && gamma == 0) {
      const std::size_t rank_w = gab_code_->Weight(err_tot);
      if (rank_w > DecodingRadius()) {
        return result_type::Failure(
            "GabidulinSkewEEADecoder: error rank weight exceeds decoding "
            "radius");
      }
      return result_type::Success(candidate, gab_code_->Unencode(candidate),
                                  rank_w, 0);
    }

    return result_type::Success(candidate, gab_code_->Unencode(candidate),
                                tau_eff, rho + gamma);
  }

 private:
  const code_type *gab_code_;
};

// ============================================================================
// GabidulinLoidreauDecoder: Welch-Berlekamp-Loidreau Interpolation Decoder
// ============================================================================

template <typename BaseField, typename ExtField>
class GabidulinLoidreauDecoder : public Decoder<ExtField> {
 public:
  using Base = Decoder<ExtField>;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using storage_type = typename ExtField::element_type;
  using code_type = GabidulinCode<BaseField, ExtField>;
  using skew_ring_type = xg::ring::FrobeniusSkewRing<ExtField>;
  using skew_poly_type = typename skew_ring_type::element_type;

  explicit GabidulinLoidreauDecoder(const AbstractCode<ExtField> &code)
      : Base(code) {
    gab_code_ = dynamic_cast<const code_type *>(&code);
    if (!gab_code_) {
      throw std::invalid_argument(
          "GabidulinLoidreauDecoder requires a GabidulinCode instance");
    }
  }

  std::string Name() const override { return "GabidulinLoidreau"; }

  std::size_t DecodingRadius() const override {
    return (gab_code_->Length() - gab_code_->Dimension()) / 2;
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto ext_field = gab_code_->Field();
    const auto &tower = gab_code_->Tower();
    const auto &skew_ring = gab_code_->SkewRing();
    const std::size_t n = gab_code_->Length();
    const std::size_t k = gab_code_->Dimension();
    const std::size_t t = DecodingRadius();

    if (!received_word.HasField() || received_word.Field() != ext_field ||
        received_word.Dimension() != 1 || received_word.Size() != n) {
      return result_type::Failure(
          "GabidulinLoidreauDecoder: invalid received word");
    }

    if (gab_code_->Contains(received_word)) {
      return result_type::Success(received_word,
                                  gab_code_->Unencode(received_word), 0, 0);
    }
    if (t == 0) {
      return result_type::Failure(
          "GabidulinLoidreauDecoder: code cannot correct rank errors");
    }

    // Build homogeneous interpolation matrix for:
    //   V(r_j) - N(g_j) == 0,  j = 0 .. n - 1,
    // with deg_q(V) <= t (t + 1 unknowns) and deg_q(N) < k + t (k + t
    // unknowns).
    const std::size_t num_v = t + 1;
    const std::size_t num_n = k + t;
    const std::size_t num_cols = num_v + num_n;
    const auto &g_pts = gab_code_->EvaluationPoints();

    using RawArray = typename GFMatrix<ExtField>::raw_array_type;
    using ShapeType = typename GFMatrix<ExtField>::shape_type;
    RawArray interp_raw = RawArray::from_shape(ShapeType{n, num_cols});
    storage_type *data = interp_raw.data();

    for (std::size_t j = 0; j < n; ++j) {
      storage_type *row_ptr = data + j * num_cols;
      storage_type r_pow = received_word.Raw(j);
      for (std::size_t i = 0; i < num_v; ++i) {
        row_ptr[i] = r_pow;
        if (i + 1 < num_v) {
          r_pow = tower.FrobeniusQ(r_pow, 1);
        }
      }
      storage_type g_pow = g_pts.Raw(j);
      for (std::size_t i = 0; i < num_n; ++i) {
        row_ptr[num_v + i] = ext_field->Neg(g_pow);
        if (i + 1 < num_n) {
          g_pow = tower.FrobeniusQ(g_pow, 1);
        }
      }
    }

    const GFMatrix<ExtField> interp_mat =
        GFMatrix<ExtField>::FromRaw(ext_field, std::move(interp_raw));
    const GFMatrix<ExtField> ns = linalg::Nullspace(interp_mat);
    if (ns.Shape(0) == 0) {
      return result_type::Failure(
          "GabidulinLoidreauDecoder: interpolation system has trivial "
          "nullspace");
    }

    // Any non-zero solution (V, N) in the nullspace satisfies N == V * f in
    // F_{q^m}[x; theta] whenever rank_weight(e) <= t, so LeftDivRem(N, V)
    // yields (f, 0).
    for (std::size_t sol_idx = 0; sol_idx < ns.Shape(0); ++sol_idx) {
      skew_poly_type v_poly(num_v);
      for (std::size_t i = 0; i < num_v; ++i) {
        v_poly[i] = ns.Raw(sol_idx, i);
      }
      v_poly = skew_ring.Normalize(std::move(v_poly));
      if (skew_ring.IsZero(v_poly)) {
        continue;
      }

      skew_poly_type n_poly(num_n);
      for (std::size_t i = 0; i < num_n; ++i) {
        n_poly[i] = ns.Raw(sol_idx, num_v + i);
      }
      n_poly = skew_ring.Normalize(std::move(n_poly));

      auto [f_poly, rem] = skew_ring.LeftDivRem(n_poly, v_poly);
      if (!skew_ring.IsZero(rem)) {
        continue;
      }
      if (skew_ring.Degree(f_poly) >= static_cast<int64_t>(k)) {
        continue;
      }

      using MsgRaw = typename message_type::raw_array_type;
      using MsgShape = typename message_type::shape_type;
      MsgRaw msg_raw = MsgRaw::from_shape(MsgShape{k});
      msg_raw.fill(ext_field->AdditiveIdentity());
      for (std::size_t i = 0; i < f_poly.size(); ++i) {
        msg_raw(i) = f_poly[i];
      }
      message_type msg = message_type::FromRaw(ext_field, std::move(msg_raw));
      codeword_type candidate = gab_code_->Encode(msg);
      const std::size_t rank_err = gab_code_->Weight(received_word - candidate);
      if (rank_err <= t) {
        return result_type::Success(std::move(candidate), std::move(msg),
                                    rank_err, 0);
      }
    }

    return result_type::Failure(
        "GabidulinLoidreauDecoder: no codeword within unique rank-decoding "
        "radius");
  }

 private:
  const code_type *gab_code_;
};

// ============================================================================
// Out-of-line GabidulinCode method definitions
// ============================================================================

template <typename BaseField, typename ExtField>
inline typename GabidulinCode<BaseField, ExtField>::result_type
GabidulinCode<BaseField, ExtField>::DecodeWithRankErasures(
    const codeword_type &received_word, const vector_type &row_erasure_values,
    const base_matrix_type &col_erasure_matrix) const {
  auto dec = this->GetDecoder("GabidulinSkewEEA");
  const auto *skew_dec =
      dynamic_cast<const GabidulinSkewEEADecoder<BaseField, ExtField> *>(
          dec.get());
  if (!skew_dec) {
    GabidulinSkewEEADecoder<BaseField, ExtField> fallback(*this);
    return fallback.DecodeWithRankErasures(received_word, row_erasure_values,
                                           col_erasure_matrix);
  }
  return skew_dec->DecodeWithRankErasures(received_word, row_erasure_values,
                                          col_erasure_matrix);
}

template <typename BaseField, typename ExtField>
inline void GabidulinCode<BaseField, ExtField>::RegisterEncodersAndDecoders() {
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<ExtField> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<ExtField> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<ExtField>>(lin);
      });
  this->RegisterDecoder(
      "GabidulinSkewEEA", [](const AbstractCode<ExtField> &code) {
        return std::make_shared<GabidulinSkewEEADecoder<BaseField, ExtField>>(
            code);
      });
  this->RegisterDecoder(
      "GabidulinLoidreau", [](const AbstractCode<ExtField> &code) {
        return std::make_shared<GabidulinLoidreauDecoder<BaseField, ExtField>>(
            code);
      });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_GABIDULIN_HPP_
