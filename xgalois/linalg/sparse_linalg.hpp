#ifndef XGALOIS_LINALG_SPARSE_LINALG_HPP_
#define XGALOIS_LINALG_SPARSE_LINALG_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace linalg {

// ============================================================================
// CSR Sparse Matrix over any FiniteField F
// ============================================================================

template <FiniteField F>
class SparseMatrix {
 public:
  using field_type = F;
  using element_type = typename F::element_type;

  SparseMatrix() = default;

  static SparseMatrix FromTriplets(
      std::shared_ptr<F> field, std::size_t rows, std::size_t cols,
      const std::vector<std::tuple<std::size_t, std::size_t, element_type>>&
          triplets) {
    if (!field) {
      throw std::invalid_argument(
          "SparseMatrix::FromTriplets: field pointer cannot be null");
    }
    struct Entry {
      std::size_t r;
      std::size_t c;
      element_type v;
    };
    std::vector<Entry> entries;
    entries.reserve(triplets.size());
    for (const auto& [r, c, v] : triplets) {
      if (r >= rows || c >= cols) {
        throw std::out_of_range(
            "SparseMatrix::FromTriplets: triplet index out of bounds");
      }
      if (!field->IsZero(v)) {
        entries.push_back(Entry{r, c, v});
      }
    }

    std::sort(entries.begin(), entries.end(),
              [](const Entry& a, const Entry& b) {
                if (a.r != b.r) {
                  return a.r < b.r;
                }
                return a.c < b.c;
              });

    SparseMatrix mat;
    mat.field_ = std::move(field);
    mat.num_rows_ = rows;
    mat.num_cols_ = cols;
    mat.row_ptr_.assign(rows + 1, 0);

    // Sum duplicate (r, c) coordinates and drop resulting zeros.
    std::vector<Entry> merged;
    merged.reserve(entries.size());
    for (const auto& e : entries) {
      if (!merged.empty() && merged.back().r == e.r && merged.back().c == e.c) {
        merged.back().v = mat.field_->Add(merged.back().v, e.v);
        if (mat.field_->IsZero(merged.back().v)) {
          merged.pop_back();
        }
      } else {
        merged.push_back(e);
      }
    }

    mat.col_idx_.reserve(merged.size());
    mat.values_.reserve(merged.size());
    for (const auto& e : merged) {
      ++mat.row_ptr_[e.r + 1];
      mat.col_idx_.push_back(e.c);
      mat.values_.push_back(e.v);
    }
    for (std::size_t i = 0; i < rows; ++i) {
      mat.row_ptr_[i + 1] += mat.row_ptr_[i];
    }
    return mat;
  }

  static SparseMatrix FromDense(const GFMatrix<F>& dense) {
    if (!dense.HasField() || dense.Dimension() != 2) {
      throw std::invalid_argument(
          "SparseMatrix::FromDense: input must be a 2-D GFMatrix with a valid "
          "field");
    }
    auto field = dense.Field();
    const std::size_t rows = dense.Shape(0);
    const std::size_t cols = dense.Shape(1);
    const element_type* data = dense.RawData().data();

    SparseMatrix mat;
    mat.field_ = std::move(field);
    mat.num_rows_ = rows;
    mat.num_cols_ = cols;
    mat.row_ptr_.assign(rows + 1, 0);

    for (std::size_t r = 0; r < rows; ++r) {
      const element_type* row_data = data + r * cols;
      for (std::size_t c = 0; c < cols; ++c) {
        if (!mat.field_->IsZero(row_data[c])) {
          mat.col_idx_.push_back(c);
          mat.values_.push_back(row_data[c]);
        }
      }
      mat.row_ptr_[r + 1] = mat.col_idx_.size();
    }
    return mat;
  }

  static SparseMatrix FromRowIndices(
      std::shared_ptr<F> field, std::size_t rows, std::size_t cols,
      const std::vector<std::vector<std::size_t>>& row_cols) {
    if (!field) {
      throw std::invalid_argument(
          "SparseMatrix::FromRowIndices: field pointer cannot be null");
    }
    if (row_cols.size() != rows) {
      throw std::invalid_argument(
          "SparseMatrix::FromRowIndices: row_cols.size() must equal rows");
    }
    const element_type one = field->MultiplicativeIdentity();
    std::vector<std::tuple<std::size_t, std::size_t, element_type>> triplets;
    for (std::size_t r = 0; r < rows; ++r) {
      for (std::size_t c : row_cols[r]) {
        triplets.emplace_back(r, c, one);
      }
    }
    return FromTriplets(std::move(field), rows, cols, triplets);
  }

  // Zero-friction duck-typed interop with coding::SparseBinaryMatrix (or any
  // CSR binary matrix providing NumRows(), NumCols(), RowPtr(), RowColIdx())
  // without introducing a header dependency on xgalois/coding/ldpc.hpp.
  template <typename BinaryCsrMatrix>
    requires requires(const BinaryCsrMatrix& m) {
      { m.NumRows() } -> std::convertible_to<std::size_t>;
      { m.NumCols() } -> std::convertible_to<std::size_t>;
      { m.RowPtr() } -> std::convertible_to<const std::vector<std::size_t>&>;
      { m.RowColIdx() } -> std::convertible_to<const std::vector<std::size_t>&>;
    }
  static SparseMatrix FromSparseBinaryMatrix(std::shared_ptr<F> field,
                                             const BinaryCsrMatrix& m) {
    if (!field) {
      throw std::invalid_argument(
          "SparseMatrix::FromSparseBinaryMatrix: field pointer cannot be null");
    }
    SparseMatrix mat;
    mat.field_ = std::move(field);
    mat.num_rows_ = static_cast<std::size_t>(m.NumRows());
    mat.num_cols_ = static_cast<std::size_t>(m.NumCols());
    mat.row_ptr_ = m.RowPtr();
    mat.col_idx_ = m.RowColIdx();
    mat.values_.assign(mat.col_idx_.size(),
                       mat.field_->MultiplicativeIdentity());
    return mat;
  }

  std::size_t NumRows() const noexcept { return num_rows_; }
  std::size_t NumCols() const noexcept { return num_cols_; }
  std::size_t NumNonZeros() const noexcept { return values_.size(); }
  std::shared_ptr<F> Field() const noexcept { return field_; }
  bool HasField() const noexcept { return field_ != nullptr; }

  const std::vector<std::size_t>& RowPtr() const noexcept { return row_ptr_; }
  const std::vector<std::size_t>& ColIdx() const noexcept { return col_idx_; }
  const std::vector<element_type>& Values() const noexcept { return values_; }

  GFMatrix<F> ToDense() const {
    RequireField("ToDense");
    using RawArray = typename GFMatrix<F>::raw_array_type;
    using ShapeType = typename GFMatrix<F>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{num_rows_, num_cols_});
    raw.fill(field_->AdditiveIdentity());
    element_type* data = raw.data();
    for (std::size_t r = 0; r < num_rows_; ++r) {
      element_type* row_data = data + r * num_cols_;
      for (std::size_t e = row_ptr_[r]; e < row_ptr_[r + 1]; ++e) {
        row_data[col_idx_[e]] = values_[e];
      }
    }
    return GFMatrix<F>::FromRaw(field_, std::move(raw));
  }

  SparseMatrix Transpose() const {
    RequireField("Transpose");
    SparseMatrix t;
    t.field_ = field_;
    t.num_rows_ = num_cols_;
    t.num_cols_ = num_rows_;
    t.row_ptr_.assign(num_cols_ + 1, 0);
    const std::size_t nnz = values_.size();
    t.col_idx_.resize(nnz);
    t.values_.resize(nnz);

    for (std::size_t e = 0; e < nnz; ++e) {
      ++t.row_ptr_[col_idx_[e] + 1];
    }
    for (std::size_t c = 0; c < num_cols_; ++c) {
      t.row_ptr_[c + 1] += t.row_ptr_[c];
    }

    std::vector<std::size_t> cursor(
        t.row_ptr_.begin(),
        t.row_ptr_.begin() + static_cast<std::ptrdiff_t>(num_cols_));
    for (std::size_t r = 0; r < num_rows_; ++r) {
      for (std::size_t e = row_ptr_[r]; e < row_ptr_[r + 1]; ++e) {
        const std::size_t dst = cursor[col_idx_[e]]++;
        t.col_idx_[dst] = r;
        t.values_[dst] = values_[e];
      }
    }
    return t;
  }

  void SpMvRaw(std::span<const element_type> x,
               std::span<element_type> y) const {
    RequireField("SpMvRaw");
    if (x.size() != num_cols_ || y.size() != num_rows_) {
      throw std::invalid_argument("SparseMatrix::SpMvRaw: dimension mismatch");
    }
    const element_type zero = field_->AdditiveIdentity();
    for (std::size_t r = 0; r < num_rows_; ++r) {
      element_type acc = zero;
      const std::size_t end = row_ptr_[r + 1];
      for (std::size_t e = row_ptr_[r]; e < end; ++e) {
        const element_type& xv = x[col_idx_[e]];
        if (!field_->IsZero(xv)) {
          acc = field_->Add(acc, field_->Mul(values_[e], xv));
        }
      }
      y[r] = std::move(acc);
    }
  }

  void SpMvTransposeRaw(std::span<const element_type> x,
                        std::span<element_type> y) const {
    RequireField("SpMvTransposeRaw");
    if (x.size() != num_rows_ || y.size() != num_cols_) {
      throw std::invalid_argument(
          "SparseMatrix::SpMvTransposeRaw: dimension mismatch");
    }
    std::fill(y.begin(), y.end(), field_->AdditiveIdentity());
    for (std::size_t r = 0; r < num_rows_; ++r) {
      const element_type& xr = x[r];
      if (field_->IsZero(xr)) {
        continue;
      }
      const std::size_t end = row_ptr_[r + 1];
      for (std::size_t e = row_ptr_[r]; e < end; ++e) {
        const std::size_t c = col_idx_[e];
        y[c] = field_->Add(y[c], field_->Mul(values_[e], xr));
      }
    }
  }

  GFArray<F> SpMv(const GFArray<F>& x) const {
    RequireField("SpMv");
    if (!x.HasField() || x.Field() != field_) {
      throw std::invalid_argument("SparseMatrix::SpMv: field mismatch");
    }
    if (x.Dimension() != 1 || x.Shape(0) != num_cols_) {
      throw std::invalid_argument(
          "SparseMatrix::SpMv: x must be 1-D with length NumCols()");
    }
    using RawArray = typename GFArray<F>::raw_array_type;
    using ShapeType = typename GFArray<F>::shape_type;
    RawArray out = RawArray::from_shape(ShapeType{num_rows_});
    SpMvRaw(x.RawSpan(), std::span<element_type>(out.data(), num_rows_));
    return GFArray<F>::FromRaw(field_, std::move(out));
  }

  GFArray<F> SpMvTranspose(const GFArray<F>& x) const {
    RequireField("SpMvTranspose");
    if (!x.HasField() || x.Field() != field_) {
      throw std::invalid_argument(
          "SparseMatrix::SpMvTranspose: field mismatch");
    }
    if (x.Dimension() != 1 || x.Shape(0) != num_rows_) {
      throw std::invalid_argument(
          "SparseMatrix::SpMvTranspose: x must be 1-D with length NumRows()");
    }
    using RawArray = typename GFArray<F>::raw_array_type;
    using ShapeType = typename GFArray<F>::shape_type;
    RawArray out = RawArray::from_shape(ShapeType{num_cols_});
    SpMvTransposeRaw(x.RawSpan(),
                     std::span<element_type>(out.data(), num_cols_));
    return GFArray<F>::FromRaw(field_, std::move(out));
  }

 private:
  void RequireField(const char* fn_name) const {
    if (!field_) {
      throw std::runtime_error(std::string("SparseMatrix::") + fn_name +
                               ": field pointer is null");
    }
  }

  std::vector<std::size_t> row_ptr_{0};
  std::vector<std::size_t> col_idx_;
  std::vector<element_type> values_;
  std::size_t num_rows_ = 0;
  std::size_t num_cols_ = 0;
  std::shared_ptr<F> field_{nullptr};
};

// ============================================================================
// 64-Bit Word-Blocked Binary SpMV (64 RHS vectors simultaneously over GF(2))
// ============================================================================

template <FiniteField BinaryF>
inline void SpMvBlock64(const SparseMatrix<BinaryF>& a,
                        std::span<const uint64_t> x, std::span<uint64_t> y) {
  if (!a.HasField() || a.Field()->Order() != 2) {
    throw std::invalid_argument("SpMvBlock64: requires a binary field GF(2)");
  }
  const std::size_t rows = a.NumRows();
  const std::size_t cols = a.NumCols();
  if (x.size() != cols || y.size() != rows) {
    throw std::invalid_argument("SpMvBlock64: span size mismatch");
  }
  const auto& row_ptr = a.RowPtr();
  const auto& col_idx = a.ColIdx();
  for (std::size_t r = 0; r < rows; ++r) {
    uint64_t acc = 0ULL;
    const std::size_t end = row_ptr[r + 1];
    for (std::size_t e = row_ptr[r]; e < end; ++e) {
      acc ^= x[col_idx[e]];
    }
    y[r] = acc;
  }
}

template <FiniteField BinaryF>
inline void SpMvTransposeBlock64(const SparseMatrix<BinaryF>& a,
                                 std::span<const uint64_t> x,
                                 std::span<uint64_t> y) {
  if (!a.HasField() || a.Field()->Order() != 2) {
    throw std::invalid_argument(
        "SpMvTransposeBlock64: requires a binary field GF(2)");
  }
  const std::size_t rows = a.NumRows();
  const std::size_t cols = a.NumCols();
  if (x.size() != rows || y.size() != cols) {
    throw std::invalid_argument("SpMvTransposeBlock64: span size mismatch");
  }
  std::fill(y.begin(), y.end(), 0ULL);
  const auto& row_ptr = a.RowPtr();
  const auto& col_idx = a.ColIdx();
  for (std::size_t r = 0; r < rows; ++r) {
    const uint64_t xr = x[r];
    if (xr == 0ULL) {
      continue;
    }
    const std::size_t end = row_ptr[r + 1];
    for (std::size_t e = row_ptr[r]; e < end; ++e) {
      y[col_idx[e]] ^= xr;
    }
  }
}

namespace detail {

// Berlekamp-Massey over any FiniteField F.
// Given s[0..N-1], returns the monic minimal annihilating polynomial
// f(x) = x^L * C(x^{-1}) in ascending degree order.
template <FiniteField F>
inline std::vector<typename F::element_type> BerlekampMasseyMinPoly(
    const F& field, std::span<const typename F::element_type> s) {
  using E = typename F::element_type;
  const std::size_t n_seq = s.size();
  const E zero = field.AdditiveIdentity();
  const E one = field.MultiplicativeIdentity();

  std::vector<E> c = {one};
  std::vector<E> b = {one};
  std::size_t l_deg = 0;
  std::size_t m_step = 1;
  E b_delta = one;

  for (std::size_t n = 0; n < n_seq; ++n) {
    E d = s[n];
    for (std::size_t i = 1; i <= l_deg; ++i) {
      if (!field.IsZero(c[i]) && !field.IsZero(s[n - i])) {
        d = field.Add(d, field.Mul(c[i], s[n - i]));
      }
    }
    if (field.IsZero(d)) {
      ++m_step;
    } else {
      std::vector<E> t = c;
      const E coef = field.Mul(d, field.Inv(b_delta));
      if (c.size() < b.size() + m_step) {
        c.resize(b.size() + m_step, zero);
      }
      for (std::size_t i = 0; i < b.size(); ++i) {
        if (!field.IsZero(b[i])) {
          c[i + m_step] = field.Sub(c[i + m_step], field.Mul(coef, b[i]));
        }
      }
      if (2 * l_deg <= n) {
        l_deg = n + 1 - l_deg;
        b = std::move(t);
        b_delta = d;
        m_step = 1;
      } else {
        ++m_step;
      }
    }
  }

  c.resize(l_deg + 1, zero);
  // Reciprocal polynomial f(x) = x^L * C(x^{-1}) has coeffs f[k] = c[L - k].
  std::vector<E> f(l_deg + 1);
  for (std::size_t k = 0; k <= l_deg; ++k) {
    f[k] = c[l_deg - k];
  }
  return f;
}

template <FiniteField F>
inline std::vector<typename F::element_type> PolyTrim(
    const F& field, std::vector<typename F::element_type> p) {
  while (!p.empty() && field.IsZero(p.back())) {
    p.pop_back();
  }
  return p;
}

template <FiniteField F>
inline std::vector<typename F::element_type> PolyMulSparse(
    const F& field, const std::vector<typename F::element_type>& a,
    const std::vector<typename F::element_type>& b) {
  using E = typename F::element_type;
  if (a.empty() || b.empty()) {
    return {};
  }
  std::vector<E> out(a.size() + b.size() - 1, field.AdditiveIdentity());
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (field.IsZero(a[i])) {
      continue;
    }
    field.FmaArray(std::span<E>(out.data() + i, b.size()),
                   std::span<const E>(b), a[i]);
  }
  return PolyTrim(field, std::move(out));
}

template <FiniteField F>
inline std::pair<std::vector<typename F::element_type>,
                 std::vector<typename F::element_type>>
PolyDivRemSparse(const F& field, const std::vector<typename F::element_type>& a,
                 const std::vector<typename F::element_type>& b) {
  using E = typename F::element_type;
  std::vector<E> rem = PolyTrim(field, a);
  std::vector<E> bn = PolyTrim(field, b);
  if (bn.empty()) {
    throw std::domain_error("PolyDivRemSparse: division by zero");
  }
  if (rem.size() < bn.size()) {
    return {{}, std::move(rem)};
  }
  const std::size_t db = bn.size() - 1;
  const E inv_lc = field.Inv(bn[db]);
  const std::size_t q_len = rem.size() - bn.size() + 1;
  std::vector<E> q(q_len, field.AdditiveIdentity());
  for (std::size_t shift = q_len; shift-- > 0;) {
    const E top = rem[shift + db];
    if (field.IsZero(top)) {
      continue;
    }
    const E factor = field.Mul(top, inv_lc);
    q[shift] = factor;
    const E neg_factor = field.Neg(factor);
    field.FmaArray(std::span<E>(rem.data() + shift, db + 1),
                   std::span<const E>(bn), neg_factor);
  }
  return {PolyTrim(field, std::move(q)), PolyTrim(field, std::move(rem))};
}

template <FiniteField F>
inline std::vector<typename F::element_type> PolyLcmSparse(
    const F& field, std::vector<typename F::element_type> a,
    std::vector<typename F::element_type> b) {
  using E = typename F::element_type;
  a = PolyTrim(field, std::move(a));
  b = PolyTrim(field, std::move(b));
  if (a.empty() || b.empty()) {
    return {};
  }
  std::vector<E> g = a;
  std::vector<E> r = b;
  while (!r.empty()) {
    auto rem = PolyDivRemSparse(field, g, r).second;
    g = std::move(r);
    r = std::move(rem);
  }
  auto q = PolyDivRemSparse(field, a, g).first;
  auto res = PolyMulSparse(field, q, b);
  if (!res.empty() && !field.IsOne(res.back())) {
    const E inv_lc = field.Inv(res.back());
    field.ScalarMulArray(std::span<E>(res), std::span<const E>(res), inv_lc);
    res.back() = field.MultiplicativeIdentity();
  }
  return res;
}

template <FiniteField F>
inline bool IsZeroVector(const F& field,
                         std::span<const typename F::element_type> v) {
  for (const auto& x : v) {
    if (!field.IsZero(x)) {
      return false;
    }
  }
  return true;
}

// 64x64 GF(2) matrix helpers for Montgomery's Block Lanczos.
using Mat64 = std::array<uint64_t, 64>;

inline Mat64 Identity64() noexcept {
  Mat64 id{};
  for (std::size_t i = 0; i < 64; ++i) {
    id[i] = uint64_t{1} << i;
  }
  return id;
}

inline Mat64 DiagMask64(uint64_t mask) noexcept {
  Mat64 d{};
  for (std::size_t i = 0; i < 64; ++i) {
    if (((mask >> i) & 1ULL) != 0ULL) {
      d[i] = uint64_t{1} << i;
    }
  }
  return d;
}

inline bool IsZeroMat64(const Mat64& m) noexcept {
  for (uint64_t row : m) {
    if (row != 0ULL) {
      return false;
    }
  }
  return true;
}

inline Mat64 MatAdd64(const Mat64& a, const Mat64& b) noexcept {
  Mat64 out{};
  for (std::size_t i = 0; i < 64; ++i) {
    out[i] = a[i] ^ b[i];
  }
  return out;
}

inline Mat64 MatMul64(const Mat64& a, const Mat64& b) noexcept {
  Mat64 out{};
  for (std::size_t i = 0; i < 64; ++i) {
    uint64_t row = a[i];
    uint64_t acc = 0ULL;
    while (row != 0ULL) {
      const std::size_t k = static_cast<std::size_t>(std::countr_zero(row));
      acc ^= b[k];
      row &= (row - 1ULL);
    }
    out[i] = acc;
  }
  return out;
}

inline Mat64 MatMulMaskRight64(const Mat64& a, uint64_t mask) noexcept {
  Mat64 out{};
  for (std::size_t i = 0; i < 64; ++i) {
    out[i] = a[i] & mask;
  }
  return out;
}

// Computes X^T * Y (64 x 64) for n x 64 bit-packed matrices X, Y.
inline Mat64 InnerProduct64(std::span<const uint64_t> x,
                            std::span<const uint64_t> y) noexcept {
  Mat64 out{};
  const std::size_t n = x.size();
  for (std::size_t i = 0; i < n; ++i) {
    uint64_t xi = x[i];
    const uint64_t yi = y[i];
    if (yi == 0ULL) {
      continue;
    }
    while (xi != 0ULL) {
      const std::size_t r = static_cast<std::size_t>(std::countr_zero(xi));
      out[r] ^= yi;
      xi &= (xi - 1ULL);
    }
  }
  return out;
}

// Computes dst ^= X * M for n x 64 matrix X and 64 x 64 matrix M.
inline void AccumulateRightMul64(std::span<uint64_t> dst,
                                 std::span<const uint64_t> x,
                                 const Mat64& m) noexcept {
  // Build 8 byte-indexed lookup tables (16 entries per nibble = 16 tables of
  // 16 uint64_t, or 8 tables of 256 -> 16 KB, fits in L1 cache).
  std::array<std::array<uint64_t, 16>, 16> lut{};
  for (std::size_t nib = 0; nib < 16; ++nib) {
    const std::size_t base = nib * 4;
    lut[nib][0] = 0ULL;
    for (std::size_t idx = 1; idx < 16; ++idx) {
      const std::size_t lsb = static_cast<std::size_t>(std::countr_zero(idx));
      lut[nib][idx] = lut[nib][idx ^ (std::size_t{1} << lsb)] ^ m[base + lsb];
    }
  }

  const std::size_t n = x.size();
  for (std::size_t i = 0; i < n; ++i) {
    uint64_t xi = x[i];
    if (xi == 0ULL) {
      continue;
    }
    uint64_t acc = 0ULL;
    for (std::size_t nib = 0; nib < 16 && xi != 0ULL; ++nib) {
      acc ^= lut[nib][xi & 0xFULL];
      xi >>= 4;
    }
    dst[i] ^= acc;
  }
}

// Finds a maximal nonsingular principal submatrix S of symmetric 64x64 GF(2)
// matrix T (prioritizing columns NOT in s_prev) and returns {W_inv, s_mask},
// where W_inv is symmetric, zero outside S x S, and (W_inv)_{S,S} =
// (T_{S,S})^{-1}.
inline std::pair<Mat64, uint64_t> FindWinv64(const Mat64& t, uint64_t s_prev) {
  std::array<std::size_t, 64> order{};
  std::size_t pos = 0;
  for (std::size_t j = 0; j < 64; ++j) {
    if (((s_prev >> j) & 1ULL) == 0ULL) {
      order[pos++] = j;
    }
  }
  for (std::size_t j = 0; j < 64; ++j) {
    if (((s_prev >> j) & 1ULL) != 0ULL) {
      order[pos++] = j;
    }
  }

  Mat64 m = t;
  std::array<bool, 64> used{};
  uint64_t s_mask = 0ULL;

  for (std::size_t idx = 0; idx < 64; ++idx) {
    const std::size_t a = order[idx];
    if (used[a]) {
      continue;
    }
    if (((m[a] >> a) & 1ULL) != 0ULL) {
      // 1x1 symmetric pivot on (a, a).
      used[a] = true;
      s_mask |= (uint64_t{1} << a);
      for (std::size_t r = 0; r < 64; ++r) {
        if (!used[r] && (((m[r] >> a) & 1ULL) != 0ULL)) {
          m[r] ^= m[a];
        }
      }
    } else {
      // Search for a 2x2 symmetric pivot [[0, 1], [1, d]] on (a, b).
      std::size_t found_b = 64;
      for (std::size_t j_idx = idx + 1; j_idx < 64; ++j_idx) {
        const std::size_t b = order[j_idx];
        if (!used[b] && (((m[a] >> b) & 1ULL) != 0ULL)) {
          found_b = b;
          break;
        }
      }
      used[a] = true;
      if (found_b < 64) {
        const std::size_t b = found_b;
        used[b] = true;
        s_mask |= (uint64_t{1} << a) | (uint64_t{1} << b);
        const bool d = (((m[b] >> b) & 1ULL) != 0ULL);
        for (std::size_t r = 0; r < 64; ++r) {
          if (used[r]) {
            continue;
          }
          const bool ur = (((m[r] >> a) & 1ULL) != 0ULL);
          const bool vr = (((m[r] >> b) & 1ULL) != 0ULL);
          if ((ur && d) != vr) {
            m[r] ^= m[a];
          }
          if (ur) {
            m[r] ^= m[b];
          }
        }
      }
    }
  }

  // Invert the guaranteed-nonsingular principal submatrix T_{S, S}.
  std::vector<std::size_t> s_cols;
  s_cols.reserve(64);
  for (std::size_t j = 0; j < 64; ++j) {
    if (((s_mask >> j) & 1ULL) != 0ULL) {
      s_cols.push_back(j);
    }
  }
  const std::size_t k = s_cols.size();
  Mat64 w_inv{};
  if (k == 0) {
    return {w_inv, 0ULL};
  }

  std::vector<uint64_t> a_sub(k, 0ULL);
  std::vector<uint64_t> i_sub(k, 0ULL);
  for (std::size_t i = 0; i < k; ++i) {
    i_sub[i] = uint64_t{1} << i;
    const uint64_t t_row = t[s_cols[i]];
    uint64_t packed = 0ULL;
    for (std::size_t j = 0; j < k; ++j) {
      if (((t_row >> s_cols[j]) & 1ULL) != 0ULL) {
        packed |= (uint64_t{1} << j);
      }
    }
    a_sub[i] = packed;
  }

  for (std::size_t col = 0; col < k; ++col) {
    const uint64_t bit = uint64_t{1} << col;
    std::size_t piv = col;
    while (piv < k && (a_sub[piv] & bit) == 0ULL) {
      ++piv;
    }
    if (piv == k) {
      throw std::runtime_error("FindWinv64: principal submatrix singular");
    }
    if (piv != col) {
      std::swap(a_sub[piv], a_sub[col]);
      std::swap(i_sub[piv], i_sub[col]);
    }
    for (std::size_t r = 0; r < k; ++r) {
      if (r != col && (a_sub[r] & bit) != 0ULL) {
        a_sub[r] ^= a_sub[col];
        i_sub[r] ^= i_sub[col];
      }
    }
  }

  for (std::size_t i = 0; i < k; ++i) {
    uint64_t unpacked = 0ULL;
    for (std::size_t j = 0; j < k; ++j) {
      if (((i_sub[i] >> j) & 1ULL) != 0ULL) {
        unpacked |= (uint64_t{1} << s_cols[j]);
      }
    }
    w_inv[s_cols[i]] = unpacked;
  }

  return {w_inv, s_mask};
}

}  // namespace detail

// ============================================================================
// Block Wiedemann Right Nullspace over any FiniteField F
// ============================================================================

// Computes non-trivial right nullspace vector(s) v (rows of the returned
// GFMatrix<F> of shape {k, n}) satisfying A * v^T == 0 for a sparse m x n
// matrix A over any FiniteField F.
template <FiniteField F>
GFMatrix<F> BlockWiedemannNullspace(const SparseMatrix<F>& a,
                                    uint64_t seed = 42) {
  using E = typename F::element_type;
  using RawArray = typename GFMatrix<F>::raw_array_type;
  using ShapeType = typename GFMatrix<F>::shape_type;

  if (!a.HasField()) {
    throw std::invalid_argument(
        "BlockWiedemannNullspace: sparse matrix has null field");
  }
  auto field = a.Field();
  const std::size_t m = a.NumRows();
  const std::size_t n = a.NumCols();
  if (n == 0) {
    return GFMatrix<F>::Zeros(field, {0, 0});
  }
  if (m == 0) {
    return Eye(field, n);
  }

  const E zero = field->AdditiveIdentity();
  const E one = field->MultiplicativeIdentity();
  std::mt19937_64 rng(seed);

  // Build square n x n operator M: F^n -> F^n such that ker(A) <= ker(M).
  // When m <= n, zero-pad A to n x n so ker(M) == ker(A) identically.
  // When m > n, apply a random sparse n x m projection P to A * x.
  std::vector<std::pair<std::size_t, E>> proj_entries;
  if (m > n) {
    proj_entries.reserve(2 * m);
    for (std::size_t j = 0; j < m; ++j) {
      const std::size_t r1 = static_cast<std::size_t>(rng() % n);
      E c1 = field->Random(rng);
      if (field->IsZero(c1)) {
        c1 = one;
      }
      proj_entries.emplace_back(r1, c1);
    }
  }

  std::vector<E> tmp_m(m, zero);
  auto apply_m = [&](std::span<const E> x, std::span<E> y) {
    a.SpMvRaw(x, std::span<E>(tmp_m));
    std::fill(y.begin(), y.end(), zero);
    if (m <= n) {
      std::copy(tmp_m.begin(), tmp_m.end(), y.begin());
    } else {
      for (std::size_t j = 0; j < m; ++j) {
        if (!field->IsZero(tmp_m[j])) {
          const auto& [row_idx, coeff] = proj_entries[j];
          y[row_idx] = field->Add(y[row_idx], field->Mul(coeff, tmp_m[j]));
        }
      }
    }
  };

  // Track linearly independent nullspace vectors discovered across trials.
  std::vector<std::vector<E>> found_vecs;
  std::vector<std::vector<E>> echelon;
  std::vector<std::size_t> pivot_cols;

  auto try_add_nullvec = [&](const std::vector<E>& cand) -> bool {
    if (detail::IsZeroVector(*field, std::span<const E>(cand))) {
      return false;
    }
    // Verify A * cand == 0.
    a.SpMvRaw(std::span<const E>(cand), std::span<E>(tmp_m));
    if (!detail::IsZeroVector(*field, std::span<const E>(tmp_m))) {
      return false;
    }
    // Check linear independence against already collected vectors.
    std::vector<E> red = cand;
    for (std::size_t idx = 0; idx < echelon.size(); ++idx) {
      const std::size_t pc = pivot_cols[idx];
      if (!field->IsZero(red[pc])) {
        const E neg_f = field->Neg(red[pc]);
        field->FmaArray(std::span<E>(red), std::span<const E>(echelon[idx]),
                        neg_f);
      }
    }
    std::size_t new_piv = n;
    for (std::size_t j = 0; j < n; ++j) {
      if (!field->IsZero(red[j])) {
        new_piv = j;
        break;
      }
    }
    if (new_piv == n) {
      return false;
    }
    const E inv_p = field->Inv(red[new_piv]);
    field->ScalarMulArray(std::span<E>(red), std::span<const E>(red), inv_p);
    echelon.push_back(std::move(red));
    pivot_cols.push_back(new_piv);
    found_vecs.push_back(cand);
    return true;
  };

  constexpr std::size_t kNumProjections = 4;
  constexpr std::size_t kMaxTrials = 12;
  const std::size_t target_vecs =
      (m < n) ? std::min<std::size_t>(n - m, 4) : std::size_t{1};
  const std::size_t seq_len = 2 * n + 2;

  std::vector<E> w(n, zero);
  std::vector<E> cur(n, zero);
  std::vector<E> next(n, zero);
  std::vector<std::vector<E>> u_vecs(kNumProjections, std::vector<E>(n, zero));
  std::vector<std::vector<E>> seqs(kNumProjections,
                                   std::vector<E>(seq_len, zero));

  for (std::size_t trial = 0;
       trial < kMaxTrials && found_vecs.size() < target_vecs; ++trial) {
    // Generate random starting vector w and block of projection vectors u.
    for (std::size_t j = 0; j < n; ++j) {
      w[j] = field->Random(rng);
    }
    if (detail::IsZeroVector(*field, std::span<const E>(w))) {
      w[trial % n] = one;
    }
    if (try_add_nullvec(w)) {
      continue;
    }

    for (std::size_t p = 0; p < kNumProjections; ++p) {
      for (std::size_t j = 0; j < n; ++j) {
        u_vecs[p][j] = field->Random(rng);
      }
    }

    // Generate Krylov scalar sequences s_p[k] = u_p^T * M^k * w.
    cur = w;
    for (std::size_t k = 0; k < seq_len; ++k) {
      for (std::size_t p = 0; p < kNumProjections; ++p) {
        seqs[p][k] = field->DotProduct(std::span<const E>(u_vecs[p]),
                                       std::span<const E>(cur));
      }
      if (k + 1 < seq_len) {
        apply_m(std::span<const E>(cur), std::span<E>(next));
        cur.swap(next);
      }
    }

    // Combine minimal polynomials via LCM across the block of projections.
    std::vector<E> min_poly = {one};
    for (std::size_t p = 0; p < kNumProjections; ++p) {
      auto fp =
          detail::BerlekampMasseyMinPoly(*field, std::span<const E>(seqs[p]));
      if (!fp.empty()) {
        min_poly =
            detail::PolyLcmSparse(*field, std::move(min_poly), std::move(fp));
      }
    }

    // Factor out powers of x: min_poly(x) = x^r * g(x) with g(0) != 0.
    std::size_t r_pow = 0;
    while (r_pow < min_poly.size() && field->IsZero(min_poly[r_pow])) {
      ++r_pow;
    }
    if (r_pow == 0 || r_pow >= min_poly.size()) {
      continue;
    }

    // Evaluate y = g(M) * w = sum_{j=0}^{deg g} min_poly[r_pow + j] * M^j * w.
    std::vector<E> y(n, zero);
    cur = w;
    const std::size_t g_len = min_poly.size() - r_pow;
    for (std::size_t j = 0; j < g_len; ++j) {
      const E& cj = min_poly[r_pow + j];
      if (!field->IsZero(cj)) {
        field->FmaArray(std::span<E>(y), std::span<const E>(cur), cj);
      }
      if (j + 1 < g_len) {
        apply_m(std::span<const E>(cur), std::span<E>(next));
        cur.swap(next);
      }
    }

    // Walk y, M*y, M^2*y, ... up to r_pow steps; last nonzero before 0 is in
    // ker(M).
    for (std::size_t step = 0; step < r_pow; ++step) {
      if (detail::IsZeroVector(*field, std::span<const E>(y))) {
        break;
      }
      apply_m(std::span<const E>(y), std::span<E>(next));
      if (detail::IsZeroVector(*field, std::span<const E>(next))) {
        try_add_nullvec(y);
        break;
      }
      y.swap(next);
    }
  }

  const std::size_t k = found_vecs.size();
  RawArray raw = RawArray::from_shape(ShapeType{k, n});
  for (std::size_t i = 0; i < k; ++i) {
    std::copy(found_vecs[i].begin(), found_vecs[i].end(), raw.data() + i * n);
  }
  return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
}

// ============================================================================
// Montgomery 64-Bit Word-Blocked Lanczos Nullspace over GF(2)
// ============================================================================

// Finds up to 64 linearly independent right-nullspace vectors v in GF(2)^n
// (each returned as std::vector<uint8_t> of length n with entries in {0, 1})
// satisfying A * v == 0 in O(n * nnz / 64) word operations.
template <FiniteField BinaryF>
std::vector<std::vector<uint8_t>> BlockLanczosNullspaceBinary(
    const SparseMatrix<BinaryF>& a, uint64_t seed = 42) {
  if (!a.HasField() || a.Field()->Order() != 2) {
    throw std::invalid_argument(
        "BlockLanczosNullspaceBinary: requires a binary field GF(2)");
  }
  const std::size_t m = a.NumRows();
  const std::size_t n = a.NumCols();
  if (n == 0) {
    return {};
  }
  if (m == 0) {
    const std::size_t k = std::min<std::size_t>(n, 64);
    std::vector<std::vector<uint8_t>> basis(k, std::vector<uint8_t>(n, 0));
    for (std::size_t i = 0; i < k; ++i) {
      basis[i][i] = 1;
    }
    return basis;
  }

  std::mt19937_64 rng(seed);
  std::vector<uint64_t> tmp_m(m, 0ULL);

  // Symmetric operator B = A^T * A acting on 64 vectors simultaneously.
  auto apply_b = [&](std::span<const uint64_t> x, std::span<uint64_t> y) {
    SpMvBlock64(a, x, std::span<uint64_t>(tmp_m));
    SpMvTransposeBlock64(a, std::span<const uint64_t>(tmp_m), y);
  };

  std::vector<uint64_t> y_rand(n, 0ULL);
  for (std::size_t i = 0; i < n; ++i) {
    y_rand[i] = rng();
  }

  std::vector<uint64_t> v0(n, 0ULL);
  apply_b(std::span<const uint64_t>(y_rand), std::span<uint64_t>(v0));

  std::vector<uint64_t> x_sol(n, 0ULL);
  std::vector<uint64_t> v_curr = v0;
  std::vector<uint64_t> v_prev1(n, 0ULL);
  std::vector<uint64_t> v_prev2(n, 0ULL);
  std::vector<uint64_t> p_curr(n, 0ULL);
  std::vector<uint64_t> v_next(n, 0ULL);

  detail::Mat64 w_inv_prev1{};
  detail::Mat64 w_inv_prev2{};
  detail::Mat64 t_prev1{};
  detail::Mat64 v_b2_v_prev1{};
  uint64_t s_prev1 = ~uint64_t{0};

  const detail::Mat64 eye64 = detail::Identity64();
  const std::size_t max_iters = (n / 64) + 64;

  for (std::size_t iter = 0; iter < max_iters; ++iter) {
    // P_i = B * V_i
    apply_b(std::span<const uint64_t>(v_curr), std::span<uint64_t>(p_curr));

    // T_i = V_i^T * B * V_i
    const detail::Mat64 t_curr = detail::InnerProduct64(
        std::span<const uint64_t>(v_curr), std::span<const uint64_t>(p_curr));
    if (detail::IsZeroMat64(t_curr)) {
      break;
    }

    auto [w_inv_curr, s_curr] = detail::FindWinv64(t_curr, s_prev1);
    if (s_curr == 0ULL) {
      break;
    }

    // X ^= V_i * W_i^{-1} * (V_i^T * V_0)
    const detail::Mat64 vt_v0 = detail::InnerProduct64(
        std::span<const uint64_t>(v_curr), std::span<const uint64_t>(v0));
    const detail::Mat64 step_x = detail::MatMul64(w_inv_curr, vt_v0);
    detail::AccumulateRightMul64(std::span<uint64_t>(x_sol),
                                 std::span<const uint64_t>(v_curr), step_x);

    // V_i^T * B^2 * V_i = (B * V_i)^T * (B * V_i) = P_i^T * P_i
    const detail::Mat64 v_b2_v_curr = detail::InnerProduct64(
        std::span<const uint64_t>(p_curr), std::span<const uint64_t>(p_curr));

    // D_{i+1} = I_64 ^ W_i^{-1} * (V_i^T B^2 V_i * S_i ^ T_i)
    const detail::Mat64 term_d = detail::MatAdd64(
        detail::MatMulMaskRight64(v_b2_v_curr, s_curr), t_curr);
    const detail::Mat64 d_next =
        detail::MatAdd64(eye64, detail::MatMul64(w_inv_curr, term_d));

    // E_{i+1} = W_{i-1}^{-1} * T_i * S_i
    const detail::Mat64 e_next = detail::MatMulMaskRight64(
        detail::MatMul64(w_inv_prev1, t_curr), s_curr);

    // F_{i+1} = W_{i-2}^{-1} * (I_64 ^ T_{i-1} * W_{i-1}^{-1})
    //           * (V_{i-1}^T B^2 V_{i-1} * S_{i-1} ^ T_{i-1}) * S_i
    const detail::Mat64 proj_prev1 =
        detail::MatAdd64(eye64, detail::MatMul64(t_prev1, w_inv_prev1));
    const detail::Mat64 term_f = detail::MatAdd64(
        detail::MatMulMaskRight64(v_b2_v_prev1, s_prev1), t_prev1);
    const detail::Mat64 f_next = detail::MatMulMaskRight64(
        detail::MatMul64(w_inv_prev2, detail::MatMul64(proj_prev1, term_f)),
        s_curr);

    // V_{i+1} = B * V_i * S_i ^ V_i * D_{i+1} ^ V_{i-1} * E_{i+1} ^ V_{i-2} *
    // F_{i+1}
    for (std::size_t row = 0; row < n; ++row) {
      v_next[row] = p_curr[row] & s_curr;
    }
    detail::AccumulateRightMul64(std::span<uint64_t>(v_next),
                                 std::span<const uint64_t>(v_curr), d_next);
    if (iter >= 1) {
      detail::AccumulateRightMul64(std::span<uint64_t>(v_next),
                                   std::span<const uint64_t>(v_prev1), e_next);
    }
    if (iter >= 2) {
      detail::AccumulateRightMul64(std::span<uint64_t>(v_next),
                                   std::span<const uint64_t>(v_prev2), f_next);
    }

    // Shift history.
    v_prev2.swap(v_prev1);
    v_prev1.swap(v_curr);
    v_curr.swap(v_next);

    w_inv_prev2 = w_inv_prev1;
    w_inv_prev1 = w_inv_curr;
    t_prev1 = t_curr;
    v_b2_v_prev1 = v_b2_v_curr;
    s_prev1 = s_curr;
  }

  // Post-processing: form the 128 candidate columns Z = [X ^ Y | V_curr],
  // compute U = A * Z (m x 128), and extract combinations c with U * c == 0.
  std::vector<uint64_t> z0(n, 0ULL);
  for (std::size_t i = 0; i < n; ++i) {
    z0[i] = x_sol[i] ^ y_rand[i];
  }
  const std::vector<uint64_t>& z1 = v_curr;

  std::vector<uint64_t> u0(m, 0ULL);
  std::vector<uint64_t> u1(m, 0ULL);
  SpMvBlock64(a, std::span<const uint64_t>(z0), std::span<uint64_t>(u0));
  SpMvBlock64(a, std::span<const uint64_t>(z1), std::span<uint64_t>(u1));

  // Pack U^T as a 128 x m bit-packed matrix and track 128 x 128 combinations.
  const std::size_t words_m = (m + 63) / 64;
  std::vector<std::vector<uint64_t>> rows_u(
      128, std::vector<uint64_t>(words_m, 0ULL));
  std::array<std::array<uint64_t, 2>, 128> comb{};
  for (std::size_t c = 0; c < 64; ++c) {
    comb[c] = {uint64_t{1} << c, 0ULL};
    comb[64 + c] = {0ULL, uint64_t{1} << c};
  }
  for (std::size_t i = 0; i < m; ++i) {
    const std::size_t wi = i >> 6;
    const uint64_t bi = uint64_t{1} << (i & 63);
    uint64_t mask0 = u0[i];
    while (mask0 != 0ULL) {
      const std::size_t c = static_cast<std::size_t>(std::countr_zero(mask0));
      rows_u[c][wi] |= bi;
      mask0 &= (mask0 - 1ULL);
    }
    uint64_t mask1 = u1[i];
    while (mask1 != 0ULL) {
      const std::size_t c = static_cast<std::size_t>(std::countr_zero(mask1));
      rows_u[64 + c][wi] |= bi;
      mask1 &= (mask1 - 1ULL);
    }
  }

  // Row-reduce the 128 rows of U^T to find all dependencies U * c == 0.
  std::size_t piv_row = 0;
  for (std::size_t col = 0; col < m && piv_row < 128; ++col) {
    const std::size_t w_idx = col >> 6;
    const uint64_t b_mask = uint64_t{1} << (col & 63);
    std::size_t found = piv_row;
    while (found < 128 && (rows_u[found][w_idx] & b_mask) == 0ULL) {
      ++found;
    }
    if (found == 128) {
      continue;
    }
    if (found != piv_row) {
      std::swap(rows_u[found], rows_u[piv_row]);
      std::swap(comb[found], comb[piv_row]);
    }
    for (std::size_t r = piv_row + 1; r < 128; ++r) {
      if ((rows_u[r][w_idx] & b_mask) != 0ULL) {
        for (std::size_t w = w_idx; w < words_m; ++w) {
          rows_u[r][w] ^= rows_u[piv_row][w];
        }
        comb[r][0] ^= comb[piv_row][0];
        comb[r][1] ^= comb[piv_row][1];
      }
    }
    ++piv_row;
  }

  // Extract up to 64 linearly independent nonzero vectors v = Z * c.
  const std::size_t words_n = (n + 63) / 64;
  std::vector<std::vector<uint64_t>> echelon_v;
  std::vector<std::size_t> echelon_piv;
  std::vector<std::vector<uint8_t>> result;

  for (std::size_t r = piv_row; r < 128 && result.size() < 64; ++r) {
    const uint64_t c0 = comb[r][0];
    const uint64_t c1 = comb[r][1];
    std::vector<uint64_t> packed_v(words_n, 0ULL);
    bool nonzero = false;
    for (std::size_t i = 0; i < n; ++i) {
      const uint8_t bit = static_cast<uint8_t>(
          (std::popcount(z0[i] & c0) ^ std::popcount(z1[i] & c1)) & 1);
      if (bit != 0) {
        packed_v[i >> 6] |= (uint64_t{1} << (i & 63));
        nonzero = true;
      }
    }
    if (!nonzero) {
      continue;
    }

    // Check linear independence via bit-packed echelon reduction.
    std::vector<uint64_t> red = packed_v;
    for (std::size_t idx = 0; idx < echelon_v.size(); ++idx) {
      const std::size_t pc = echelon_piv[idx];
      if (((red[pc >> 6] >> (pc & 63)) & 1ULL) != 0ULL) {
        for (std::size_t w = pc >> 6; w < words_n; ++w) {
          red[w] ^= echelon_v[idx][w];
        }
      }
    }

    std::size_t new_piv = n;
    for (std::size_t w = 0; w < words_n; ++w) {
      if (red[w] != 0ULL) {
        new_piv = (w << 6) + static_cast<std::size_t>(std::countr_zero(red[w]));
        break;
      }
    }
    if (new_piv >= n) {
      continue;
    }

    echelon_v.push_back(std::move(red));
    echelon_piv.push_back(new_piv);

    std::vector<uint8_t> vec(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
      vec[i] = static_cast<uint8_t>((packed_v[i >> 6] >> (i & 63)) & 1ULL);
    }
    result.push_back(std::move(vec));
  }

  return result;
}

}  // namespace linalg
}  // namespace xg

#endif  // XGALOIS_LINALG_SPARSE_LINALG_HPP_
