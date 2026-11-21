#ifndef XGALOIS_LINALG_LINALG_HPP_
#define XGALOIS_LINALG_LINALG_HPP_

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {

// Forward declaration enables zero-overhead compile-time dispatch to the
// bit-packed GF(2) M4RI engine without introducing a header dependency on
// xgalois/field/gf_binary.hpp.
class GaloisFieldBinary;

namespace linalg {

template <typename GaloisField>
struct SystematicFormResult {
  GFArray<GaloisField> g_sys;
  std::vector<std::size_t> pivot_cols;
  std::vector<std::size_t> free_cols;
  std::vector<std::size_t> col_perm;
};

template <typename GaloisField>
struct PluResult {
  GFArray<GaloisField> p;
  GFArray<GaloisField> l;
  GFArray<GaloisField> u;
  std::vector<std::size_t> perm;
  std::size_t num_swaps = 0;
};

namespace detail {

template <typename GaloisField>
inline void RequireField(const std::shared_ptr<GaloisField> &field,
                         const char *fn_name) {
  if (!field) {
    throw std::invalid_argument(std::string(fn_name) +
                                ": field pointer cannot be null");
  }
}

template <typename GaloisField>
inline void RequireValidArray(const GFArray<GaloisField> &a,
                              const char *fn_name) {
  if (!a.HasField()) {
    throw std::invalid_argument(std::string(fn_name) +
                                ": GFArray has null field pointer");
  }
}

template <typename GaloisField>
inline void RequireSameField(const GFArray<GaloisField> &a,
                             const GFArray<GaloisField> &b,
                             const char *fn_name) {
  RequireValidArray(a, fn_name);
  RequireValidArray(b, fn_name);
  if (a.Field() != b.Field()) {
    throw std::invalid_argument(std::string(fn_name) +
                                ": arrays must belong to the same field");
  }
}

// ============================================================================
// Bit-packed uint64_t GF(2) M4RI Engine
// ============================================================================

struct PackedBinaryMatrix {
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::size_t words_per_row = 0;
  std::vector<uint64_t> words;

  static PackedBinaryMatrix Zeros(std::size_t r, std::size_t c) {
    PackedBinaryMatrix m;
    m.rows = r;
    m.cols = c;
    m.words_per_row = (c + 63) / 64;
    m.words.assign(r * m.words_per_row, 0ULL);
    return m;
  }

  static PackedBinaryMatrix Identity(std::size_t n) {
    PackedBinaryMatrix m = Zeros(n, n);
    for (std::size_t i = 0; i < n; ++i) {
      m.SetBit(i, i);
    }
    return m;
  }

  uint64_t *RowPtr(std::size_t r) noexcept {
    return words.data() + r * words_per_row;
  }

  const uint64_t *RowPtr(std::size_t r) const noexcept {
    return words.data() + r * words_per_row;
  }

  bool GetBit(std::size_t r, std::size_t c) const noexcept {
    return ((RowPtr(r)[c >> 6] >> (c & 63)) & 1ULL) != 0;
  }

  void SetBit(std::size_t r, std::size_t c) noexcept {
    RowPtr(r)[c >> 6] |= (uint64_t{1} << (c & 63));
  }

  void XorRow(std::size_t dst_row, std::size_t src_row,
              std::size_t start_word = 0) noexcept {
    uint64_t *dst = RowPtr(dst_row);
    const uint64_t *src = RowPtr(src_row);
    for (std::size_t w = start_word; w < words_per_row; ++w) {
      dst[w] ^= src[w];
    }
  }

  void SwapRows(std::size_t r1, std::size_t r2) noexcept {
    if (r1 == r2) {
      return;
    }
    uint64_t *p1 = RowPtr(r1);
    uint64_t *p2 = RowPtr(r2);
    for (std::size_t w = 0; w < words_per_row; ++w) {
      std::swap(p1[w], p2[w]);
    }
  }

  void SwapRowPrefixBits(std::size_t r1, std::size_t r2,
                         std::size_t num_bits) noexcept {
    if (r1 == r2 || num_bits == 0) {
      return;
    }
    uint64_t *p1 = RowPtr(r1);
    uint64_t *p2 = RowPtr(r2);
    const std::size_t full_words = num_bits >> 6;
    const std::size_t rem_bits = num_bits & 63;
    for (std::size_t w = 0; w < full_words; ++w) {
      std::swap(p1[w], p2[w]);
    }
    if (rem_bits > 0) {
      const uint64_t mask = (uint64_t{1} << rem_bits) - 1ULL;
      const uint64_t diff = (p1[full_words] ^ p2[full_words]) & mask;
      p1[full_words] ^= diff;
      p2[full_words] ^= diff;
    }
  }
};

inline PackedBinaryMatrix PackMatrix(const uint8_t *src, std::size_t rows,
                                     std::size_t cols) {
  PackedBinaryMatrix m = PackedBinaryMatrix::Zeros(rows, cols);
  for (std::size_t r = 0; r < rows; ++r) {
    const uint8_t *src_row = src + r * cols;
    uint64_t *dst_row = m.RowPtr(r);
    std::size_t c = 0;
    for (std::size_t w = 0; w < m.words_per_row; ++w) {
      uint64_t word = 0;
      const std::size_t limit = std::min<std::size_t>(64, cols - c);
      for (std::size_t bit = 0; bit < limit; ++bit) {
        word |= static_cast<uint64_t>(src_row[c + bit] & 1U) << bit;
      }
      dst_row[w] = word;
      c += limit;
    }
  }
  return m;
}

inline void UnpackMatrix(const PackedBinaryMatrix &m, uint8_t *dst) noexcept {
  for (std::size_t r = 0; r < m.rows; ++r) {
    const uint64_t *src_row = m.RowPtr(r);
    uint8_t *dst_row = dst + r * m.cols;
    std::size_t c = 0;
    for (std::size_t w = 0; w < m.words_per_row; ++w) {
      const uint64_t word = src_row[w];
      const std::size_t limit = std::min<std::size_t>(64, m.cols - c);
      for (std::size_t bit = 0; bit < limit; ++bit) {
        dst_row[c + bit] = static_cast<uint8_t>((word >> bit) & 1ULL);
      }
      c += limit;
    }
  }
}

// M4RI-style Gray-code 4-bit / 8-bit table-driven bit-packed GF(2) matrix
// multiplication: C (m x n) = A (m x k) * B (k x n).
inline void MatMulM4ri(const uint8_t *a_data, std::size_t m, std::size_t k,
                       const uint8_t *b_data, std::size_t n, uint8_t *c_data) {
  if (m == 0 || n == 0) {
    return;
  }
  if (k == 0) {
    std::fill(c_data, c_data + m * n, uint8_t{0});
    return;
  }

  const PackedBinaryMatrix packed_a = PackMatrix(a_data, m, k);
  const PackedBinaryMatrix packed_b = PackMatrix(b_data, k, n);
  PackedBinaryMatrix packed_c = PackedBinaryMatrix::Zeros(m, n);
  const std::size_t words_n = packed_b.words_per_row;

  // For tiny row counts (e.g., 1D * 2D), direct set-bit scan avoids table
  // setup overhead.
  if (m < 4) {
    for (std::size_t i = 0; i < m; ++i) {
      const uint64_t *a_row = packed_a.RowPtr(i);
      uint64_t *c_row = packed_c.RowPtr(i);
      for (std::size_t w_k = 0; w_k < packed_a.words_per_row; ++w_k) {
        uint64_t word = a_row[w_k];
        const std::size_t base_k = w_k << 6;
        while (word != 0) {
          const std::size_t bit =
              static_cast<std::size_t>(std::countr_zero(word));
          const uint64_t *b_row = packed_b.RowPtr(base_k + bit);
          for (std::size_t w = 0; w < words_n; ++w) {
            c_row[w] ^= b_row[w];
          }
          word &= (word - 1);
        }
      }
    }
    UnpackMatrix(packed_c, c_data);
    return;
  }

  // Choose 8-bit or 4-bit blocks so the Gray-code lookup table stays in L1.
  // Both 4 and 8 divide 64, guaranteeing every block lies within a single word.
  const std::size_t k_step = (m >= 32 && words_n <= 32) ? 8 : 4;
  std::vector<uint64_t> table((std::size_t{1} << k_step) * words_n, 0ULL);

  for (std::size_t k0 = 0; k0 < k; k0 += k_step) {
    const std::size_t b_len = std::min(k_step, k - k0);
    const std::size_t table_size = std::size_t{1} << b_len;
    std::fill(table.begin(),
              table.begin() + static_cast<std::ptrdiff_t>(words_n), 0ULL);

    // Gray-code table construction: 1 row XOR per table entry.
    std::size_t g_prev = 0;
    for (std::size_t t = 1; t < table_size; ++t) {
      const std::size_t g_curr = t ^ (t >> 1);
      const std::size_t diff = g_prev ^ g_curr;
      const std::size_t bit = static_cast<std::size_t>(std::countr_zero(diff));
      const uint64_t *prev_row = table.data() + g_prev * words_n;
      const uint64_t *b_row = packed_b.RowPtr(k0 + bit);
      uint64_t *curr_row = table.data() + g_curr * words_n;
      for (std::size_t w = 0; w < words_n; ++w) {
        curr_row[w] = prev_row[w] ^ b_row[w];
      }
      g_prev = g_curr;
    }

    const std::size_t word_idx = k0 >> 6;
    const std::size_t bit_offset = k0 & 63;
    const uint64_t mask_bits = (uint64_t{1} << b_len) - 1ULL;
    for (std::size_t i = 0; i < m; ++i) {
      const std::size_t idx = static_cast<std::size_t>(
          (packed_a.RowPtr(i)[word_idx] >> bit_offset) & mask_bits);
      if (idx == 0) {
        continue;
      }
      uint64_t *c_row = packed_c.RowPtr(i);
      const uint64_t *t_row = table.data() + idx * words_n;
      for (std::size_t w = 0; w < words_n; ++w) {
        c_row[w] ^= t_row[w];
      }
    }
  }

  UnpackMatrix(packed_c, c_data);
}

inline void MatVecGf2(const uint8_t *a_data, std::size_t m, std::size_t k,
                      const uint8_t *b_data, uint8_t *out) {
  if (m == 0) {
    return;
  }
  if (k == 0) {
    std::fill(out, out + m, uint8_t{0});
    return;
  }
  const PackedBinaryMatrix packed_a = PackMatrix(a_data, m, k);
  const PackedBinaryMatrix packed_b = PackMatrix(b_data, 1, k);
  const uint64_t *b_row = packed_b.RowPtr(0);
  const std::size_t words_k = packed_a.words_per_row;
  for (std::size_t i = 0; i < m; ++i) {
    const uint64_t *a_row = packed_a.RowPtr(i);
    uint64_t acc = 0;
    for (std::size_t w = 0; w < words_k; ++w) {
      acc ^= (a_row[w] & b_row[w]);
    }
    out[i] = static_cast<uint8_t>(std::popcount(acc) & 1);
  }
}

inline std::vector<std::size_t> RowEchelonPacked(PackedBinaryMatrix &m) {
  std::vector<std::size_t> pivots;
  pivots.reserve(std::min(m.rows, m.cols));
  std::size_t pivot_row = 0;
  for (std::size_t col = 0; col < m.cols && pivot_row < m.rows; ++col) {
    const std::size_t w_idx = col >> 6;
    const uint64_t bit_mask = uint64_t{1} << (col & 63);
    std::size_t found = pivot_row;
    while (found < m.rows && (m.RowPtr(found)[w_idx] & bit_mask) == 0) {
      ++found;
    }
    if (found == m.rows) {
      continue;
    }
    if (found != pivot_row) {
      m.SwapRows(pivot_row, found);
    }
    pivots.push_back(col);
    for (std::size_t k = pivot_row + 1; k < m.rows; ++k) {
      if ((m.RowPtr(k)[w_idx] & bit_mask) != 0) {
        m.XorRow(k, pivot_row, w_idx);
      }
    }
    ++pivot_row;
  }
  return pivots;
}

inline std::vector<std::size_t> RrefPacked(PackedBinaryMatrix &m) {
  std::vector<std::size_t> pivots;
  pivots.reserve(std::min(m.rows, m.cols));
  std::size_t pivot_row = 0;
  for (std::size_t col = 0; col < m.cols && pivot_row < m.rows; ++col) {
    const std::size_t w_idx = col >> 6;
    const uint64_t bit_mask = uint64_t{1} << (col & 63);
    std::size_t found = pivot_row;
    while (found < m.rows && (m.RowPtr(found)[w_idx] & bit_mask) == 0) {
      ++found;
    }
    if (found == m.rows) {
      continue;
    }
    if (found != pivot_row) {
      m.SwapRows(pivot_row, found);
    }
    pivots.push_back(col);
    for (std::size_t k = 0; k < m.rows; ++k) {
      if (k != pivot_row && (m.RowPtr(k)[w_idx] & bit_mask) != 0) {
        m.XorRow(k, pivot_row, w_idx);
      }
    }
    ++pivot_row;
  }
  return pivots;
}

inline uint8_t DetPacked(const uint8_t *a_data, std::size_t n) {
  if (n == 0) {
    return 1;
  }
  PackedBinaryMatrix m = PackMatrix(a_data, n, n);
  for (std::size_t col = 0; col < n; ++col) {
    const std::size_t w_idx = col >> 6;
    const uint64_t bit_mask = uint64_t{1} << (col & 63);
    std::size_t found = col;
    while (found < n && (m.RowPtr(found)[w_idx] & bit_mask) == 0) {
      ++found;
    }
    if (found == n) {
      return 0;
    }
    if (found != col) {
      m.SwapRows(col, found);
    }
    for (std::size_t k = col + 1; k < n; ++k) {
      if ((m.RowPtr(k)[w_idx] & bit_mask) != 0) {
        m.XorRow(k, col, w_idx);
      }
    }
  }
  return 1;
}

inline void SolvePacked(const uint8_t *a_data, std::size_t n,
                        const uint8_t *b_data, std::size_t rhs_cols,
                        uint8_t *out_data) {
  if (n == 0 || rhs_cols == 0) {
    return;
  }
  PackedBinaryMatrix m = PackMatrix(a_data, n, n);
  PackedBinaryMatrix rhs = PackMatrix(b_data, n, rhs_cols);
  for (std::size_t col = 0; col < n; ++col) {
    const std::size_t w_idx = col >> 6;
    const uint64_t bit_mask = uint64_t{1} << (col & 63);
    std::size_t found = col;
    while (found < n && (m.RowPtr(found)[w_idx] & bit_mask) == 0) {
      ++found;
    }
    if (found == n) {
      throw std::runtime_error("Matrix is singular");
    }
    if (found != col) {
      m.SwapRows(col, found);
      rhs.SwapRows(col, found);
    }
    for (std::size_t k = 0; k < n; ++k) {
      if (k != col && (m.RowPtr(k)[w_idx] & bit_mask) != 0) {
        m.XorRow(k, col, w_idx);
        rhs.XorRow(k, col, 0);
      }
    }
  }
  UnpackMatrix(rhs, out_data);
}

inline void InvPacked(const uint8_t *a_data, std::size_t n, uint8_t *out_data) {
  if (n == 0) {
    return;
  }
  PackedBinaryMatrix m = PackMatrix(a_data, n, n);
  PackedBinaryMatrix inv_m = PackedBinaryMatrix::Identity(n);
  for (std::size_t col = 0; col < n; ++col) {
    const std::size_t w_idx = col >> 6;
    const uint64_t bit_mask = uint64_t{1} << (col & 63);
    std::size_t found = col;
    while (found < n && (m.RowPtr(found)[w_idx] & bit_mask) == 0) {
      ++found;
    }
    if (found == n) {
      throw std::runtime_error("Matrix is singular");
    }
    if (found != col) {
      m.SwapRows(col, found);
      inv_m.SwapRows(col, found);
    }
    for (std::size_t k = 0; k < n; ++k) {
      if (k != col && (m.RowPtr(k)[w_idx] & bit_mask) != 0) {
        m.XorRow(k, col, w_idx);
        inv_m.XorRow(k, col, 0);
      }
    }
  }
  UnpackMatrix(inv_m, out_data);
}

// ============================================================================
// General Contiguous RREF Core (returns {RREF matrix, pivot_cols})
// ============================================================================

template <typename GaloisField>
inline std::pair<GFArray<GaloisField>, std::vector<std::size_t>> RrefWithPivots(
    const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  RequireValidArray(a, "Rref");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("Rref requires a 2-D array (matrix)");
  }

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  if (rows == 0 || cols == 0) {
    return {a, {}};
  }

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    PackedBinaryMatrix packed = PackMatrix(a.RawData().data(), rows, cols);
    std::vector<std::size_t> pivots = RrefPacked(packed);
    RawArray raw = RawArray::from_shape(
        typename GFArray<GaloisField>::shape_type{rows, cols});
    UnpackMatrix(packed, raw.data());
    return {GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw)),
            std::move(pivots)};
  } else {
    RawArray raw = a.RawData();
    StorageType *data = raw.data();
    const StorageType zero = field->AdditiveIdentity();
    const StorageType one = field->MultiplicativeIdentity();

    std::vector<std::size_t> pivots;
    pivots.reserve(std::min(rows, cols));

    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
      std::size_t found = pivot_row;
      if (data[pivot_row * cols + col] == zero) {
        found = pivot_row + 1;
        while (found < rows && data[found * cols + col] == zero) {
          ++found;
        }
      }
      if (found == rows) {
        continue;
      }

      if (found != pivot_row) {
        StorageType *r1 = data + pivot_row * cols + col;
        StorageType *r2 = data + found * cols + col;
        std::swap_ranges(r1, r1 + (cols - col), r2);
      }

      pivots.push_back(col);

      StorageType *p_row = data + pivot_row * cols;
      const StorageType pivot_val = p_row[col];
      if (pivot_val != one) {
        const StorageType inv_pivot = field->Inv(pivot_val);
        p_row[col] = one;
        if (col + 1 < cols) {
          std::span<StorageType> tail(p_row + col + 1, cols - col - 1);
          field->ScalarMulArray(tail, tail, inv_pivot);
        }
      }

      const std::size_t tail_len = cols - col - 1;
      std::span<const StorageType> p_tail(p_row + col + 1, tail_len);

      for (std::size_t k = 0; k < rows; ++k) {
        if (k == pivot_row) {
          continue;
        }
        StorageType *k_row = data + k * cols;
        const StorageType factor = k_row[col];
        if (factor == zero) {
          continue;
        }
        k_row[col] = zero;
        if (tail_len > 0) {
          const StorageType neg_factor = field->Neg(factor);
          field->FmaArray(std::span<StorageType>(k_row + col + 1, tail_len),
                          p_tail, neg_factor);
        }
      }

      ++pivot_row;
    }

    return {GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw)),
            std::move(pivots)};
  }
}

}  // namespace detail

// ============================================================================
// Matrix / Array Construction & Structural Utilities
// ============================================================================

template <typename GaloisField>
GFArray<GaloisField> Eye(const std::shared_ptr<GaloisField> &field,
                         std::size_t n) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  detail::RequireField(field, "Eye");

  RawArray raw =
      RawArray::from_shape(typename GFArray<GaloisField>::shape_type{n, n});
  raw.fill(field->AdditiveIdentity());
  const auto one = field->MultiplicativeIdentity();
  auto *data = raw.data();
  for (std::size_t i = 0; i < n; ++i) {
    data[i * n + i] = one;
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
}

template <typename GaloisField>
GFArray<GaloisField> Zeros(const std::shared_ptr<GaloisField> &field,
                           const std::vector<std::size_t> &shape) {
  detail::RequireField(field, "Zeros");
  return GFArray<GaloisField>::Zeros(field, shape);
}

template <typename GaloisField>
GFArray<GaloisField> Ones(const std::shared_ptr<GaloisField> &field,
                          const std::vector<std::size_t> &shape) {
  detail::RequireField(field, "Ones");
  return GFArray<GaloisField>::Ones(field, shape);
}

template <typename GaloisField>
GFArray<GaloisField> RandomMatrix(const std::shared_ptr<GaloisField> &field,
                                  const std::vector<std::size_t> &shape,
                                  std::mt19937_64 &rng) {
  detail::RequireField(field, "RandomMatrix");
  return GFArray<GaloisField>::Random(field, shape, rng);
}

template <typename GaloisField>
GFArray<GaloisField> RandomMatrix(const std::shared_ptr<GaloisField> &field,
                                  const std::vector<std::size_t> &shape) {
  detail::RequireField(field, "RandomMatrix");
  return GFArray<GaloisField>::Random(field, shape);
}

template <typename GaloisField>
GFArray<GaloisField> Transpose(const GFArray<GaloisField> &a) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  detail::RequireValidArray(a, "Transpose");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("Transpose requires a 2-D array (matrix)");
  }

  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  RawArray out = RawArray::from_shape(
      typename GFArray<GaloisField>::shape_type{cols, rows});
  if (rows == 0 || cols == 0) {
    return GFArray<GaloisField>::FromRaw(a.Field(), std::move(out));
  }

  const auto *src = a.RawData().data();
  auto *dst = out.data();
  constexpr std::size_t kTile = 16;
  for (std::size_t i0 = 0; i0 < rows; i0 += kTile) {
    const std::size_t i_max = std::min(i0 + kTile, rows);
    for (std::size_t j0 = 0; j0 < cols; j0 += kTile) {
      const std::size_t j_max = std::min(j0 + kTile, cols);
      for (std::size_t i = i0; i < i_max; ++i) {
        const auto *src_row = src + i * cols;
        for (std::size_t j = j0; j < j_max; ++j) {
          dst[j * rows + i] = src_row[j];
        }
      }
    }
  }
  return GFArray<GaloisField>::FromRaw(a.Field(), std::move(out));
}

template <typename GaloisField>
GFArray<GaloisField> HStack(const GFArray<GaloisField> &a,
                            const GFArray<GaloisField> &b) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  detail::RequireSameField(a, b, "HStack");

  if (a.Dimension() == 1 && b.Dimension() == 1) {
    const std::size_t na = a.Shape(0);
    const std::size_t nb = b.Shape(0);
    RawArray out = RawArray::from_shape(
        typename GFArray<GaloisField>::shape_type{na + nb});
    if (na > 0) {
      std::copy(a.RawData().data(), a.RawData().data() + na, out.data());
    }
    if (nb > 0) {
      std::copy(b.RawData().data(), b.RawData().data() + nb, out.data() + na);
    }
    return GFArray<GaloisField>::FromRaw(a.Field(), std::move(out));
  }

  if (a.Dimension() == 2 && b.Dimension() == 2) {
    if (a.Shape(0) != b.Shape(0)) {
      throw std::invalid_argument(
          "HStack: 2-D arrays must have the same number of rows");
    }
    const std::size_t rows = a.Shape(0);
    const std::size_t cols_a = a.Shape(1);
    const std::size_t cols_b = b.Shape(1);
    const std::size_t total_cols = cols_a + cols_b;
    RawArray out = RawArray::from_shape(
        typename GFArray<GaloisField>::shape_type{rows, total_cols});
    const auto *pa = a.RawData().data();
    const auto *pb = b.RawData().data();
    auto *dst = out.data();
    for (std::size_t r = 0; r < rows; ++r) {
      if (cols_a > 0) {
        std::copy(pa + r * cols_a, pa + (r + 1) * cols_a, dst + r * total_cols);
      }
      if (cols_b > 0) {
        std::copy(pb + r * cols_b, pb + (r + 1) * cols_b,
                  dst + r * total_cols + cols_a);
      }
    }
    return GFArray<GaloisField>::FromRaw(a.Field(), std::move(out));
  }

  throw std::invalid_argument("HStack requires two 1-D or two 2-D arrays");
}

template <typename GaloisField>
GFArray<GaloisField> VStack(const GFArray<GaloisField> &a,
                            const GFArray<GaloisField> &b) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  detail::RequireSameField(a, b, "VStack");

  if (a.Dimension() == 1 && b.Dimension() == 1) {
    if (a.Shape(0) != b.Shape(0)) {
      throw std::invalid_argument(
          "VStack: 1-D arrays must have the same length");
    }
    const std::size_t cols = a.Shape(0);
    RawArray out = RawArray::from_shape(
        typename GFArray<GaloisField>::shape_type{2, cols});
    if (cols > 0) {
      std::copy(a.RawData().data(), a.RawData().data() + cols, out.data());
      std::copy(b.RawData().data(), b.RawData().data() + cols,
                out.data() + cols);
    }
    return GFArray<GaloisField>::FromRaw(a.Field(), std::move(out));
  }

  if (a.Dimension() == 2 && b.Dimension() == 2) {
    if (a.Shape(1) != b.Shape(1)) {
      throw std::invalid_argument(
          "VStack: 2-D arrays must have the same number of columns");
    }
    const std::size_t rows_a = a.Shape(0);
    const std::size_t rows_b = b.Shape(0);
    const std::size_t cols = a.Shape(1);
    RawArray out = RawArray::from_shape(
        typename GFArray<GaloisField>::shape_type{rows_a + rows_b, cols});
    const std::size_t count_a = rows_a * cols;
    const std::size_t count_b = rows_b * cols;
    if (count_a > 0) {
      std::copy(a.RawData().data(), a.RawData().data() + count_a, out.data());
    }
    if (count_b > 0) {
      std::copy(b.RawData().data(), b.RawData().data() + count_b,
                out.data() + count_a);
    }
    return GFArray<GaloisField>::FromRaw(a.Field(), std::move(out));
  }

  throw std::invalid_argument("VStack requires two 1-D or two 2-D arrays");
}

// ============================================================================
// Products: Dot, Vdot, Outer, Cross, Kron, MatrixPower
// ============================================================================

template <typename GaloisField>
GFArray<GaloisField> Dot(const GFArray<GaloisField> &a,
                         const GFArray<GaloisField> &b) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireSameField(a, b, "Dot");
  auto field = a.Field();

  // 1D * 1D -> 0-D GFArray
  if (a.Dimension() == 1 && b.Dimension() == 1) {
    if (a.Shape(0) != b.Shape(0)) {
      throw std::invalid_argument(
          "Vector dimensions must match for dot product");
    }
    StorageType val = field->DotProduct(a.RawSpan(), b.RawSpan());
    RawArray raw = RawArray::from_shape(ShapeType{});
    raw() = std::move(val);
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  }

  // 1D * 2D -> 1-D GFArray (critical for linear code encoding u * G)
  if (a.Dimension() == 1 && b.Dimension() == 2) {
    if (a.Shape(0) != b.Shape(0)) {
      throw std::invalid_argument(
          "Vector size must match matrix rows for 1D-2D dot product");
    }
    const std::size_t k = a.Shape(0);
    const std::size_t n = b.Shape(1);
    RawArray raw = RawArray::from_shape(ShapeType{n});
    if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
      detail::MatMulM4ri(a.RawData().data(), 1, k, b.RawData().data(), n,
                         raw.data());
    } else {
      const StorageType zero = field->AdditiveIdentity();
      raw.fill(zero);
      if (k > 0 && n > 0) {
        const StorageType *a_data = a.RawData().data();
        const StorageType *b_data = b.RawData().data();
        std::span<StorageType> c_span(raw.data(), n);
        for (std::size_t idx = 0; idx < k; ++idx) {
          const StorageType &a_val = a_data[idx];
          if (a_val == zero) {
            continue;
          }
          field->FmaArray(
              c_span, std::span<const StorageType>(b_data + idx * n, n), a_val);
        }
      }
    }
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  }

  // 2D * 1D -> 1-D GFArray
  if (a.Dimension() == 2 && b.Dimension() == 1) {
    if (a.Shape(1) != b.Shape(0)) {
      throw std::invalid_argument("Matrix columns must match vector size");
    }
    const std::size_t m = a.Shape(0);
    const std::size_t k = a.Shape(1);
    RawArray raw = RawArray::from_shape(ShapeType{m});
    if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
      detail::MatVecGf2(a.RawData().data(), m, k, b.RawData().data(),
                        raw.data());
    } else {
      const StorageType *a_data = a.RawData().data();
      const std::span<const StorageType> b_span = b.RawSpan();
      StorageType *out = raw.data();
      for (std::size_t i = 0; i < m; ++i) {
        out[i] = field->DotProduct(
            std::span<const StorageType>(a_data + i * k, k), b_span);
      }
    }
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  }

  // 2D * 2D -> 2-D GFArray
  if (a.Dimension() == 2 && b.Dimension() == 2) {
    if (a.Shape(1) != b.Shape(0)) {
      throw std::invalid_argument(
          "Matrix dimensions incompatible for multiplication");
    }
    const std::size_t m = a.Shape(0);
    const std::size_t k = a.Shape(1);
    const std::size_t n = b.Shape(1);
    RawArray raw = RawArray::from_shape(ShapeType{m, n});
    if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
      detail::MatMulM4ri(a.RawData().data(), m, k, b.RawData().data(), n,
                         raw.data());
    } else {
      const StorageType zero = field->AdditiveIdentity();
      raw.fill(zero);
      if (m > 0 && k > 0 && n > 0) {
        const StorageType *a_data = a.RawData().data();
        const StorageType *b_data = b.RawData().data();
        StorageType *c_data = raw.data();
        // Cache-friendly i-k-j order dispatching each row update to FmaArray.
        for (std::size_t i = 0; i < m; ++i) {
          std::span<StorageType> c_row(c_data + i * n, n);
          const StorageType *a_row = a_data + i * k;
          for (std::size_t kk = 0; kk < k; ++kk) {
            const StorageType &a_ik = a_row[kk];
            if (a_ik == zero) {
              continue;
            }
            field->FmaArray(
                c_row, std::span<const StorageType>(b_data + kk * n, n), a_ik);
          }
        }
      }
    }
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  }

  throw std::invalid_argument("Unsupported array dimensions for dot product");
}

template <typename GaloisField>
GaloisFieldElement<GaloisField> Vdot(const GFArray<GaloisField> &a,
                                     const GFArray<GaloisField> &b) {
  detail::RequireSameField(a, b, "Vdot");
  if (a.Size() != b.Size()) {
    throw std::invalid_argument(
        "Arrays must have same number of elements for Vdot");
  }
  auto field = a.Field();
  auto val = field->DotProduct(a.RawSpan(), b.RawSpan());
  return GaloisFieldElement<GaloisField>::FromInternal(std::move(val),
                                                       std::move(field));
}

template <typename GaloisField>
GFArray<GaloisField> Outer(const GFArray<GaloisField> &a,
                           const GFArray<GaloisField> &b) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  detail::RequireSameField(a, b, "Outer");
  if (a.Dimension() != 1 || b.Dimension() != 1) {
    throw std::invalid_argument("Outer product requires 1-D arrays");
  }

  auto field = a.Field();
  const std::size_t m = a.Shape(0);
  const std::size_t n = b.Shape(0);
  RawArray raw =
      RawArray::from_shape(typename GFArray<GaloisField>::shape_type{m, n});
  if (m > 0 && n > 0) {
    const StorageType *a_data = a.RawData().data();
    const std::span<const StorageType> b_span = b.RawSpan();
    StorageType *out = raw.data();
    for (std::size_t i = 0; i < m; ++i) {
      field->ScalarMulArray(std::span<StorageType>(out + i * n, n), b_span,
                            a_data[i]);
    }
  }
  return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
}

template <typename GaloisField>
GFArray<GaloisField> Cross(const GFArray<GaloisField> &a,
                           const GFArray<GaloisField> &b) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  detail::RequireSameField(a, b, "Cross");
  if (a.Dimension() != 1 || b.Dimension() != 1) {
    throw std::invalid_argument("Cross product requires 1-D arrays (vectors)");
  }
  if ((a.Shape(0) != 2 && a.Shape(0) != 3) ||
      (b.Shape(0) != 2 && b.Shape(0) != 3)) {
    throw std::invalid_argument(
        "Cross product vectors must have 2 or 3 elements");
  }

  auto field = a.Field();
  const StorageType zero = field->AdditiveIdentity();
  const StorageType *pa = a.RawData().data();
  const StorageType *pb = b.RawData().data();

  const StorageType &a0 = pa[0];
  const StorageType &a1 = pa[1];
  const StorageType &a2 = (a.Shape(0) == 3) ? pa[2] : zero;

  const StorageType &b0 = pb[0];
  const StorageType &b1 = pb[1];
  const StorageType &b2 = (b.Shape(0) == 3) ? pb[2] : zero;

  RawArray raw =
      RawArray::from_shape(typename GFArray<GaloisField>::shape_type{3});
  StorageType *out = raw.data();
  out[0] = field->Sub(field->Mul(a1, b2), field->Mul(a2, b1));
  out[1] = field->Sub(field->Mul(a2, b0), field->Mul(a0, b2));
  out[2] = field->Sub(field->Mul(a0, b1), field->Mul(a1, b0));

  return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
}

template <typename GaloisField>
GFArray<GaloisField> Kron(const GFArray<GaloisField> &a,
                          const GFArray<GaloisField> &b) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  detail::RequireSameField(a, b, "Kron");
  if (a.Dimension() != 2 || b.Dimension() != 2) {
    throw std::invalid_argument("Kron requires 2-D arrays (matrices)");
  }

  auto field = a.Field();
  const std::size_t a_rows = a.Shape(0);
  const std::size_t a_cols = a.Shape(1);
  const std::size_t b_rows = b.Shape(0);
  const std::size_t b_cols = b.Shape(1);
  const std::size_t out_rows = a_rows * b_rows;
  const std::size_t out_cols = a_cols * b_cols;

  RawArray raw = RawArray::from_shape(
      typename GFArray<GaloisField>::shape_type{out_rows, out_cols});
  if (out_rows == 0 || out_cols == 0) {
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  }

  const StorageType *a_data = a.RawData().data();
  const StorageType *b_data = b.RawData().data();
  StorageType *out = raw.data();

  for (std::size_t i = 0; i < a_rows; ++i) {
    for (std::size_t k = 0; k < b_rows; ++k) {
      std::span<const StorageType> b_row(b_data + k * b_cols, b_cols);
      StorageType *out_row = out + (i * b_rows + k) * out_cols;
      for (std::size_t j = 0; j < a_cols; ++j) {
        field->ScalarMulArray(
            std::span<StorageType>(out_row + j * b_cols, b_cols), b_row,
            a_data[i * a_cols + j]);
      }
    }
  }
  return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
}

// Forward declaration for negative powers in MatrixPower.
template <typename GaloisField>
GFArray<GaloisField> Inv(const GFArray<GaloisField> &a);

template <typename GaloisField>
GFArray<GaloisField> MatrixPower(const GFArray<GaloisField> &a, int64_t n) {
  detail::RequireValidArray(a, "MatrixPower");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("MatrixPower requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t size = a.Shape(0);

  if (n == 0) {
    return Eye(field, size);
  }
  if (n == 1) {
    return a;
  }

  GFArray<GaloisField> base = (n < 0) ? Inv(a) : a;
  uint64_t exp = (n < 0) ? (static_cast<uint64_t>(-(n + 1)) + 1ULL)
                         : static_cast<uint64_t>(n);

  if (exp == 1) {
    return base;
  }

  GFArray<GaloisField> result;
  bool has_result = false;
  while (exp > 0) {
    if ((exp & 1ULL) != 0) {
      result = has_result ? Dot(result, base) : base;
      has_result = true;
    }
    exp >>= 1;
    if (exp > 0) {
      base = Dot(base, base);
    }
  }
  return result;
}

// ============================================================================
// Elimination, Rank, Trace, Determinant, Solve, Inverse
// ============================================================================

template <typename GaloisField>
GFArray<GaloisField> RowEchelon(const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  detail::RequireValidArray(a, "RowEchelon");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("RowEchelon requires a 2-D array (matrix)");
  }

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  if (rows == 0 || cols == 0) {
    return a;
  }

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    detail::PackedBinaryMatrix packed =
        detail::PackMatrix(a.RawData().data(), rows, cols);
    (void)detail::RowEchelonPacked(packed);
    RawArray raw = RawArray::from_shape(
        typename GFArray<GaloisField>::shape_type{rows, cols});
    detail::UnpackMatrix(packed, raw.data());
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  } else {
    RawArray raw = a.RawData();
    StorageType *data = raw.data();
    const StorageType zero = field->AdditiveIdentity();
    const StorageType one = field->MultiplicativeIdentity();

    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
      std::size_t found = pivot_row;
      if (data[pivot_row * cols + col] == zero) {
        found = pivot_row + 1;
        while (found < rows && data[found * cols + col] == zero) {
          ++found;
        }
      }
      if (found == rows) {
        continue;
      }

      if (found != pivot_row) {
        StorageType *r1 = data + pivot_row * cols + col;
        StorageType *r2 = data + found * cols + col;
        std::swap_ranges(r1, r1 + (cols - col), r2);
      }

      StorageType *p_row = data + pivot_row * cols;
      const StorageType pivot_val = p_row[col];
      if (pivot_val != one) {
        const StorageType inv_pivot = field->Inv(pivot_val);
        p_row[col] = one;
        if (col + 1 < cols) {
          std::span<StorageType> tail(p_row + col + 1, cols - col - 1);
          field->ScalarMulArray(tail, tail, inv_pivot);
        }
      }

      const std::size_t tail_len = cols - col - 1;
      std::span<const StorageType> p_tail(p_row + col + 1, tail_len);
      for (std::size_t k = pivot_row + 1; k < rows; ++k) {
        StorageType *k_row = data + k * cols;
        const StorageType factor = k_row[col];
        if (factor == zero) {
          continue;
        }
        k_row[col] = zero;
        if (tail_len > 0) {
          const StorageType neg_factor = field->Neg(factor);
          field->FmaArray(std::span<StorageType>(k_row + col + 1, tail_len),
                          p_tail, neg_factor);
        }
      }

      ++pivot_row;
    }

    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
  }
}

template <typename GaloisField>
GFArray<GaloisField> Rref(const GFArray<GaloisField> &a) {
  return detail::RrefWithPivots(a).first;
}

template <typename GaloisField>
GaloisFieldElement<GaloisField> Trace(const GFArray<GaloisField> &a) {
  detail::RequireValidArray(a, "Trace");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("Trace requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  auto acc = field->AdditiveIdentity();
  const auto *data = a.RawData().data();
  for (std::size_t i = 0; i < n; ++i) {
    acc = field->Add(acc, data[i * n + i]);
  }
  return GaloisFieldElement<GaloisField>::FromInternal(std::move(acc),
                                                       std::move(field));
}

template <typename GaloisField>
GaloisFieldElement<GaloisField> Det(const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  detail::RequireValidArray(a, "Det");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("Det requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  if (n == 0) {
    auto one = field->MultiplicativeIdentity();
    return GaloisFieldElement<GaloisField>::FromInternal(std::move(one),
                                                         std::move(field));
  }

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    const uint8_t d = detail::DetPacked(a.RawData().data(), n);
    return GaloisFieldElement<GaloisField>::FromInternal(d, std::move(field));
  } else {
    RawArray raw = a.RawData();
    StorageType *data = raw.data();
    const StorageType zero = field->AdditiveIdentity();
    StorageType det_val = field->MultiplicativeIdentity();
    bool sign_flip = false;

    for (std::size_t i = 0; i < n; ++i) {
      std::size_t pivot_row = i;
      if (data[i * n + i] == zero) {
        pivot_row = i + 1;
        while (pivot_row < n && data[pivot_row * n + i] == zero) {
          ++pivot_row;
        }
      }
      if (pivot_row == n) {
        return GaloisFieldElement<GaloisField>::FromInternal(zero,
                                                             std::move(field));
      }

      if (pivot_row != i) {
        StorageType *r1 = data + i * n + i;
        StorageType *r2 = data + pivot_row * n + i;
        std::swap_ranges(r1, r1 + (n - i), r2);
        sign_flip = !sign_flip;
      }

      const StorageType *p_row = data + i * n;
      const StorageType pivot_val = p_row[i];
      det_val = field->Mul(det_val, pivot_val);

      if (i + 1 < n) {
        const StorageType inv_pivot = field->Inv(pivot_val);
        const std::size_t tail_len = n - i - 1;
        std::span<const StorageType> p_tail(p_row + i + 1, tail_len);
        for (std::size_t k = i + 1; k < n; ++k) {
          StorageType *k_row = data + k * n;
          const StorageType entry = k_row[i];
          if (entry == zero) {
            continue;
          }
          const StorageType neg_factor =
              field->Neg(field->Mul(entry, inv_pivot));
          k_row[i] = zero;
          field->FmaArray(std::span<StorageType>(k_row + i + 1, tail_len),
                          p_tail, neg_factor);
        }
      }
    }

    if (sign_flip) {
      det_val = field->Neg(det_val);
    }
    return GaloisFieldElement<GaloisField>::FromInternal(std::move(det_val),
                                                         std::move(field));
  }
}

template <typename GaloisField>
std::size_t MatrixRank(const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;

  detail::RequireValidArray(a, "MatrixRank");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("MatrixRank requires a 2-D array (matrix)");
  }

  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  if (rows == 0 || cols == 0) {
    return 0;
  }

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    detail::PackedBinaryMatrix packed =
        detail::PackMatrix(a.RawData().data(), rows, cols);
    return detail::RowEchelonPacked(packed).size();
  } else {
    auto field = a.Field();
    RawArray raw = a.RawData();
    StorageType *data = raw.data();
    const StorageType zero = field->AdditiveIdentity();

    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
      std::size_t found = pivot_row;
      if (data[pivot_row * cols + col] == zero) {
        found = pivot_row + 1;
        while (found < rows && data[found * cols + col] == zero) {
          ++found;
        }
      }
      if (found == rows) {
        continue;
      }

      if (found != pivot_row) {
        StorageType *r1 = data + pivot_row * cols + col;
        StorageType *r2 = data + found * cols + col;
        std::swap_ranges(r1, r1 + (cols - col), r2);
      }

      const StorageType *p_row = data + pivot_row * cols;
      const StorageType inv_pivot = field->Inv(p_row[col]);
      const std::size_t tail_len = cols - col - 1;
      if (tail_len > 0) {
        std::span<const StorageType> p_tail(p_row + col + 1, tail_len);
        for (std::size_t k = pivot_row + 1; k < rows; ++k) {
          StorageType *k_row = data + k * cols;
          const StorageType entry = k_row[col];
          if (entry == zero) {
            continue;
          }
          const StorageType neg_factor =
              field->Neg(field->Mul(entry, inv_pivot));
          k_row[col] = zero;
          field->FmaArray(std::span<StorageType>(k_row + col + 1, tail_len),
                          p_tail, neg_factor);
        }
      }

      ++pivot_row;
    }
    return pivot_row;
  }
}

template <typename GaloisField>
GFArray<GaloisField> Solve(const GFArray<GaloisField> &a,
                           const GFArray<GaloisField> &b) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireSameField(a, b, "Solve");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("Solve requires a square coefficient matrix");
  }

  const std::size_t n = a.Shape(0);
  const bool is_1d = (b.Dimension() == 1);
  if (!is_1d && b.Dimension() != 2) {
    throw std::invalid_argument(
        "Right-hand side must be a 1-D vector or 2-D matrix");
  }
  if (b.Shape(0) != n) {
    throw std::invalid_argument(
        "Right-hand side row count must match matrix size");
  }

  auto field = a.Field();
  const std::size_t rhs_cols = is_1d ? 1 : b.Shape(1);
  if (n == 0 || rhs_cols == 0) {
    return b;
  }

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    RawArray out = is_1d ? RawArray::from_shape(ShapeType{n})
                         : RawArray::from_shape(ShapeType{n, rhs_cols});
    detail::SolvePacked(a.RawData().data(), n, b.RawData().data(), rhs_cols,
                        out.data());
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(out));
  } else {
    RawArray mat = a.RawData();
    RawArray rhs = b.RawData();
    StorageType *m_data = mat.data();
    StorageType *r_data = rhs.data();
    const StorageType zero = field->AdditiveIdentity();
    const StorageType one = field->MultiplicativeIdentity();

    // Forward elimination to unit upper triangular form.
    for (std::size_t i = 0; i < n; ++i) {
      std::size_t pivot_row = i;
      if (m_data[i * n + i] == zero) {
        pivot_row = i + 1;
        while (pivot_row < n && m_data[pivot_row * n + i] == zero) {
          ++pivot_row;
        }
      }
      if (pivot_row == n) {
        throw std::runtime_error("Matrix is singular");
      }

      if (pivot_row != i) {
        StorageType *m1 = m_data + i * n + i;
        StorageType *m2 = m_data + pivot_row * n + i;
        std::swap_ranges(m1, m1 + (n - i), m2);

        StorageType *r1 = r_data + i * rhs_cols;
        StorageType *r2 = r_data + pivot_row * rhs_cols;
        std::swap_ranges(r1, r1 + rhs_cols, r2);
      }

      StorageType *m_row = m_data + i * n;
      StorageType *r_row = r_data + i * rhs_cols;
      const StorageType pivot_val = m_row[i];
      if (pivot_val != one) {
        const StorageType inv_pivot = field->Inv(pivot_val);
        m_row[i] = one;
        if (i + 1 < n) {
          std::span<StorageType> m_tail(m_row + i + 1, n - i - 1);
          field->ScalarMulArray(m_tail, m_tail, inv_pivot);
        }
        std::span<StorageType> r_span(r_row, rhs_cols);
        field->ScalarMulArray(r_span, r_span, inv_pivot);
      }

      const std::size_t tail_len = n - i - 1;
      std::span<const StorageType> m_tail(m_row + i + 1, tail_len);
      std::span<const StorageType> r_src(r_row, rhs_cols);
      for (std::size_t k = i + 1; k < n; ++k) {
        StorageType *mk_row = m_data + k * n;
        const StorageType factor = mk_row[i];
        if (factor == zero) {
          continue;
        }
        const StorageType neg_factor = field->Neg(factor);
        mk_row[i] = zero;
        if (tail_len > 0) {
          field->FmaArray(std::span<StorageType>(mk_row + i + 1, tail_len),
                          m_tail, neg_factor);
        }
        field->FmaArray(std::span<StorageType>(r_data + k * rhs_cols, rhs_cols),
                        r_src, neg_factor);
      }
    }

    // Back-substitution (division-free since diagonal of mat is all ones).
    for (std::size_t i = n; i-- > 0;) {
      std::span<const StorageType> r_src(r_data + i * rhs_cols, rhs_cols);
      for (std::size_t k = 0; k < i; ++k) {
        const StorageType factor = m_data[k * n + i];
        if (factor == zero) {
          continue;
        }
        const StorageType neg_factor = field->Neg(factor);
        field->FmaArray(std::span<StorageType>(r_data + k * rhs_cols, rhs_cols),
                        r_src, neg_factor);
      }
    }

    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(rhs));
  }
}

template <typename GaloisField>
GFArray<GaloisField> Inv(const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireValidArray(a, "Inv");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("Inv requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  if (n == 0) {
    return GFArray<GaloisField>::Zeros(std::move(field), {0, 0});
  }

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    RawArray out = RawArray::from_shape(ShapeType{n, n});
    detail::InvPacked(a.RawData().data(), n, out.data());
    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(out));
  } else {
    RawArray mat = a.RawData();
    RawArray inv_raw = RawArray::from_shape(ShapeType{n, n});
    const StorageType zero = field->AdditiveIdentity();
    const StorageType one = field->MultiplicativeIdentity();
    inv_raw.fill(zero);
    StorageType *inv_data = inv_raw.data();
    for (std::size_t i = 0; i < n; ++i) {
      inv_data[i * n + i] = one;
    }

    StorageType *m_data = mat.data();

    for (std::size_t i = 0; i < n; ++i) {
      std::size_t pivot_row = i;
      if (m_data[i * n + i] == zero) {
        pivot_row = i + 1;
        while (pivot_row < n && m_data[pivot_row * n + i] == zero) {
          ++pivot_row;
        }
      }
      if (pivot_row == n) {
        throw std::runtime_error("Matrix is singular");
      }

      if (pivot_row != i) {
        StorageType *m1 = m_data + i * n + i;
        StorageType *m2 = m_data + pivot_row * n + i;
        std::swap_ranges(m1, m1 + (n - i), m2);

        StorageType *i1 = inv_data + i * n;
        StorageType *i2 = inv_data + pivot_row * n;
        std::swap_ranges(i1, i1 + n, i2);
      }

      StorageType *m_row = m_data + i * n;
      StorageType *i_row = inv_data + i * n;
      const StorageType pivot_val = m_row[i];
      if (pivot_val != one) {
        const StorageType inv_pivot = field->Inv(pivot_val);
        m_row[i] = one;
        if (i + 1 < n) {
          std::span<StorageType> m_tail(m_row + i + 1, n - i - 1);
          field->ScalarMulArray(m_tail, m_tail, inv_pivot);
        }
        std::span<StorageType> i_span(i_row, n);
        field->ScalarMulArray(i_span, i_span, inv_pivot);
      }

      const std::size_t tail_len = n - i - 1;
      std::span<const StorageType> m_tail(m_row + i + 1, tail_len);
      std::span<const StorageType> i_src(i_row, n);
      for (std::size_t k = 0; k < n; ++k) {
        if (k == i) {
          continue;
        }
        StorageType *mk_row = m_data + k * n;
        const StorageType factor = mk_row[i];
        if (factor == zero) {
          continue;
        }
        const StorageType neg_factor = field->Neg(factor);
        mk_row[i] = zero;
        if (tail_len > 0) {
          field->FmaArray(std::span<StorageType>(mk_row + i + 1, tail_len),
                          m_tail, neg_factor);
        }
        field->FmaArray(std::span<StorageType>(inv_data + k * n, n), i_src,
                        neg_factor);
      }
    }

    return GFArray<GaloisField>::FromRaw(std::move(field), std::move(inv_raw));
  }
}

// ============================================================================
// Fundamental Subspaces, Systematic Form, PLU, Characteristic Polynomial
// ============================================================================

template <typename GaloisField>
GFArray<GaloisField> Nullspace(const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireValidArray(a, "Nullspace");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("Nullspace requires a 2-D array (matrix)");
  }

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);

  if (cols == 0) {
    return GFArray<GaloisField>::Zeros(std::move(field), {0, 0});
  }
  if (rows == 0) {
    return Eye(field, cols);
  }

  auto [rref_mat, pivot_cols] = detail::RrefWithPivots(a);
  std::vector<bool> is_pivot(cols, false);
  for (std::size_t p : pivot_cols) {
    is_pivot[p] = true;
  }

  std::vector<std::size_t> free_cols;
  free_cols.reserve(cols - pivot_cols.size());
  for (std::size_t c = 0; c < cols; ++c) {
    if (!is_pivot[c]) {
      free_cols.push_back(c);
    }
  }

  const std::size_t nullity = free_cols.size();
  if (nullity == 0) {
    return GFArray<GaloisField>::Zeros(std::move(field), {0, cols});
  }

  RawArray raw = RawArray::from_shape(ShapeType{nullity, cols});
  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();
  raw.fill(zero);

  const StorageType *r_data = rref_mat.RawData().data();
  StorageType *n_data = raw.data();
  const std::size_t rank = pivot_cols.size();

  for (std::size_t idx = 0; idx < nullity; ++idx) {
    const std::size_t f = free_cols[idx];
    StorageType *basis_row = n_data + idx * cols;
    basis_row[f] = one;
    for (std::size_t p = 0; p < rank; ++p) {
      const StorageType &val = r_data[p * cols + f];
      if (val != zero) {
        basis_row[pivot_cols[p]] = field->Neg(val);
      }
    }
  }

  return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
}

template <typename GaloisField>
GFArray<GaloisField> LeftNullspace(const GFArray<GaloisField> &a) {
  detail::RequireValidArray(a, "LeftNullspace");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("LeftNullspace requires a 2-D array (matrix)");
  }
  return Nullspace(Transpose(a));
}

template <typename GaloisField>
GFArray<GaloisField> RowSpace(const GFArray<GaloisField> &a) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireValidArray(a, "RowSpace");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("RowSpace requires a 2-D array (matrix)");
  }

  auto field = a.Field();
  const std::size_t cols = a.Shape(1);
  auto [rref_mat, pivot_cols] = detail::RrefWithPivots(a);
  const std::size_t rank = pivot_cols.size();
  if (rank == 0) {
    return GFArray<GaloisField>::Zeros(std::move(field), {0, cols});
  }

  RawArray raw = RawArray::from_shape(ShapeType{rank, cols});
  const std::size_t total = rank * cols;
  if (total > 0) {
    std::copy(rref_mat.RawData().data(), rref_mat.RawData().data() + total,
              raw.data());
  }
  return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
}

template <typename GaloisField>
GFArray<GaloisField> ColumnSpace(const GFArray<GaloisField> &a) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireValidArray(a, "ColumnSpace");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("ColumnSpace requires a 2-D array (matrix)");
  }

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  auto [rref_mat, pivot_cols] = detail::RrefWithPivots(a);
  (void)rref_mat;
  const std::size_t rank = pivot_cols.size();
  if (rank == 0) {
    return GFArray<GaloisField>::Zeros(std::move(field), {rows, 0});
  }

  RawArray raw = RawArray::from_shape(ShapeType{rows, rank});
  const auto *src = a.RawData().data();
  auto *dst = raw.data();
  for (std::size_t r = 0; r < rows; ++r) {
    const auto *src_row = src + r * cols;
    auto *dst_row = dst + r * rank;
    for (std::size_t j = 0; j < rank; ++j) {
      dst_row[j] = src_row[pivot_cols[j]];
    }
  }
  return GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw));
}

template <typename GaloisField>
SystematicFormResult<GaloisField> SystematicForm(
    const GFArray<GaloisField> &g) {
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireValidArray(g, "SystematicForm");
  if (g.Dimension() != 2) {
    throw std::invalid_argument("SystematicForm requires a 2-D array (matrix)");
  }

  auto field = g.Field();
  const std::size_t cols = g.Shape(1);
  auto [rref_mat, pivot_cols] = detail::RrefWithPivots(g);

  std::vector<bool> is_pivot(cols, false);
  for (std::size_t p : pivot_cols) {
    is_pivot[p] = true;
  }

  std::vector<std::size_t> free_cols;
  free_cols.reserve(cols - pivot_cols.size());
  for (std::size_t c = 0; c < cols; ++c) {
    if (!is_pivot[c]) {
      free_cols.push_back(c);
    }
  }

  std::vector<std::size_t> col_perm;
  col_perm.reserve(cols);
  col_perm.insert(col_perm.end(), pivot_cols.begin(), pivot_cols.end());
  col_perm.insert(col_perm.end(), free_cols.begin(), free_cols.end());

  const std::size_t k = pivot_cols.size();
  RawArray raw = RawArray::from_shape(ShapeType{k, cols});
  if (k > 0 && cols > 0) {
    const auto *r_data = rref_mat.RawData().data();
    auto *dst = raw.data();
    for (std::size_t r = 0; r < k; ++r) {
      const auto *r_row = r_data + r * cols;
      auto *dst_row = dst + r * cols;
      for (std::size_t j = 0; j < cols; ++j) {
        dst_row[j] = r_row[col_perm[j]];
      }
    }
  }

  return SystematicFormResult<GaloisField>{
      GFArray<GaloisField>::FromRaw(std::move(field), std::move(raw)),
      std::move(pivot_cols), std::move(free_cols), std::move(col_perm)};
}

template <typename GaloisField>
PluResult<GaloisField> Plu(const GFArray<GaloisField> &a) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  detail::RequireValidArray(a, "Plu");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("Plu requires a 2-D array (matrix)");
  }

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  const std::size_t k_max = std::min(rows, cols);

  std::vector<std::size_t> perm(rows);
  std::iota(perm.begin(), perm.end(), std::size_t{0});
  std::size_t num_swaps = 0;

  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    detail::PackedBinaryMatrix packed_u =
        detail::PackMatrix(a.RawData().data(), rows, cols);
    detail::PackedBinaryMatrix packed_l =
        detail::PackedBinaryMatrix::Identity(rows);

    for (std::size_t col = 0; col < k_max; ++col) {
      const std::size_t w_idx = col >> 6;
      const uint64_t bit_mask = uint64_t{1} << (col & 63);
      std::size_t pivot_row = col;
      while (pivot_row < rows &&
             (packed_u.RowPtr(pivot_row)[w_idx] & bit_mask) == 0) {
        ++pivot_row;
      }
      if (pivot_row == rows) {
        continue;
      }
      if (pivot_row != col) {
        packed_u.SwapRows(col, pivot_row);
        packed_l.SwapRowPrefixBits(col, pivot_row, col);
        std::swap(perm[col], perm[pivot_row]);
        ++num_swaps;
      }
      for (std::size_t r = col + 1; r < rows; ++r) {
        if ((packed_u.RowPtr(r)[w_idx] & bit_mask) != 0) {
          packed_l.SetBit(r, col);
          packed_u.XorRow(r, col, w_idx);
        }
      }
    }

    RawArray p_raw = RawArray::from_shape(ShapeType{rows, rows});
    p_raw.fill(uint8_t{0});
    for (std::size_t i = 0; i < rows; ++i) {
      p_raw.data()[i * rows + perm[i]] = uint8_t{1};
    }
    RawArray l_raw = RawArray::from_shape(ShapeType{rows, rows});
    RawArray u_raw = RawArray::from_shape(ShapeType{rows, cols});
    detail::UnpackMatrix(packed_l, l_raw.data());
    detail::UnpackMatrix(packed_u, u_raw.data());

    return PluResult<GaloisField>{
        GFArray<GaloisField>::FromRaw(field, std::move(p_raw)),
        GFArray<GaloisField>::FromRaw(field, std::move(l_raw)),
        GFArray<GaloisField>::FromRaw(field, std::move(u_raw)), std::move(perm),
        num_swaps};
  } else {
    const StorageType zero = field->AdditiveIdentity();
    const StorageType one = field->MultiplicativeIdentity();

    RawArray u_raw = a.RawData();
    RawArray l_raw = RawArray::from_shape(ShapeType{rows, rows});
    l_raw.fill(zero);
    StorageType *l_data = l_raw.data();
    for (std::size_t i = 0; i < rows; ++i) {
      l_data[i * rows + i] = one;
    }
    StorageType *u_data = u_raw.data();

    for (std::size_t col = 0; col < k_max; ++col) {
      std::size_t pivot_row = col;
      if (u_data[col * cols + col] == zero) {
        pivot_row = col + 1;
        while (pivot_row < rows && u_data[pivot_row * cols + col] == zero) {
          ++pivot_row;
        }
      }
      if (pivot_row == rows) {
        continue;
      }

      if (pivot_row != col) {
        StorageType *u1 = u_data + col * cols + col;
        StorageType *u2 = u_data + pivot_row * cols + col;
        std::swap_ranges(u1, u1 + (cols - col), u2);

        if (col > 0) {
          StorageType *l1 = l_data + col * rows;
          StorageType *l2 = l_data + pivot_row * rows;
          std::swap_ranges(l1, l1 + col, l2);
        }

        std::swap(perm[col], perm[pivot_row]);
        ++num_swaps;
      }

      const StorageType *u_pivot_row = u_data + col * cols;
      const StorageType inv_pivot = field->Inv(u_pivot_row[col]);
      const std::size_t tail_len = cols - col - 1;
      std::span<const StorageType> u_tail(u_pivot_row + col + 1, tail_len);

      for (std::size_t r = col + 1; r < rows; ++r) {
        StorageType *u_row = u_data + r * cols;
        const StorageType entry = u_row[col];
        if (entry == zero) {
          continue;
        }
        const StorageType factor = field->Mul(entry, inv_pivot);
        l_data[r * rows + col] = factor;
        u_row[col] = zero;
        if (tail_len > 0) {
          const StorageType neg_factor = field->Neg(factor);
          field->FmaArray(std::span<StorageType>(u_row + col + 1, tail_len),
                          u_tail, neg_factor);
        }
      }
    }

    RawArray p_raw = RawArray::from_shape(ShapeType{rows, rows});
    p_raw.fill(zero);
    StorageType *p_data = p_raw.data();
    for (std::size_t i = 0; i < rows; ++i) {
      p_data[i * rows + perm[i]] = one;
    }

    return PluResult<GaloisField>{
        GFArray<GaloisField>::FromRaw(field, std::move(p_raw)),
        GFArray<GaloisField>::FromRaw(field, std::move(l_raw)),
        GFArray<GaloisField>::FromRaw(field, std::move(u_raw)), std::move(perm),
        num_swaps};
  }
}

// Exact division-free Berkowitz algorithm for det(xI - A) over ANY finite
// field of ANY characteristic (including characteristic 2 and small primes).
template <typename GaloisField>
PolynomialDense<GaloisField> CharacteristicPoly(const GFArray<GaloisField> &a,
                                                const std::string &var = "x") {
  using StorageType = typename GaloisField::element_type;

  detail::RequireValidArray(a, "CharacteristicPoly");
  if (a.Dimension() != 2 || a.Shape(0) != a.Shape(1)) {
    throw std::invalid_argument("CharacteristicPoly requires a square matrix");
  }

  auto field = a.Field();
  const std::size_t n = a.Shape(0);
  if (n == 0) {
    return PolynomialDense<GaloisField>::One(std::move(field), var);
  }

  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();
  const StorageType *a_data = a.RawData().data();

  // q stores coefficients of det(xI_s - A_trailing) in descending degree order:
  // q[0] = 1 (coeff of x^s), ..., q[s] = constant term.
  std::vector<StorageType> q = {one};
  q.reserve(n + 1);

  std::vector<StorageType> r_vec;
  std::vector<StorageType> v_vec;
  std::vector<StorageType> next_v;
  std::vector<StorageType> toeplitz;
  std::vector<StorageType> next_q;

  // Build from bottom-right 1x1 principal submatrix up to the full n x n A.
  for (std::size_t idx = n; idx-- > 0;) {
    const std::size_t sub_size = n - idx - 1;  // Size of trailing submatrix M
    const std::size_t cur_size = sub_size + 1;

    toeplitz.assign(cur_size + 1, zero);
    toeplitz[0] = one;
    toeplitz[1] = field->Neg(a_data[idx * n + idx]);

    if (sub_size > 0) {
      r_vec.resize(sub_size);
      v_vec.resize(sub_size);
      next_v.resize(sub_size);
      for (std::size_t j = 0; j < sub_size; ++j) {
        r_vec[j] = a_data[idx * n + (idx + 1 + j)];
        v_vec[j] = a_data[(idx + 1 + j) * n + idx];
      }

      std::span<const StorageType> r_span(r_vec.data(), sub_size);
      for (std::size_t j = 0; j < sub_size; ++j) {
        const StorageType dot_rv = field->DotProduct(
            r_span, std::span<const StorageType>(v_vec.data(), sub_size));
        toeplitz[j + 2] = field->Neg(dot_rv);

        if (j + 1 < sub_size) {
          std::span<const StorageType> v_span(v_vec.data(), sub_size);
          for (std::size_t row = 0; row < sub_size; ++row) {
            const StorageType *m_row = a_data + (idx + 1 + row) * n + (idx + 1);
            next_v[row] = field->DotProduct(
                std::span<const StorageType>(m_row, sub_size), v_span);
          }
          v_vec.swap(next_v);
        }
      }
    }

    next_q.assign(cur_size + 1, zero);
    for (std::size_t k = 0; k <= cur_size; ++k) {
      const std::size_t r_min = (k > sub_size) ? (k - sub_size) : 0;
      StorageType acc = zero;
      for (std::size_t r = r_min; r <= k; ++r) {
        acc = field->Add(acc, field->Mul(toeplitz[r], q[k - r]));
      }
      next_q[k] = std::move(acc);
    }
    q.swap(next_q);
  }

  // Convert from descending degree order to ascending order for
  // PolynomialDense.
  std::vector<StorageType> raw_coeffs(n + 1);
  for (std::size_t i = 0; i <= n; ++i) {
    raw_coeffs[i] = std::move(q[n - i]);
  }
  return PolynomialDense<GaloisField>::FromRaw(std::move(raw_coeffs),
                                               std::move(field), var);
}

}  // namespace linalg
}  // namespace xg

#endif  // XGALOIS_LINALG_LINALG_HPP_
