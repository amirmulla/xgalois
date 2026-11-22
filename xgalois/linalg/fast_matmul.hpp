#ifndef XGALOIS_LINALG_FAST_MATMUL_HPP_
#define XGALOIS_LINALG_FAST_MATMUL_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace linalg {

namespace detail {

// Cache-blocked GEMM over any FiniteField F using field->FmaArray.
// Computes C (m x n) = A (m x k) * B (k x n) on contiguous row-major buffers.
template <FiniteField F>
inline void BlockedFmaGemmRaw(const F& field, const typename F::element_type* a,
                              const typename F::element_type* b,
                              typename F::element_type* c, std::size_t m,
                              std::size_t k, std::size_t n) {
  using E = typename F::element_type;
  const E zero = field.AdditiveIdentity();
  std::fill(c, c + m * n, zero);
  if (m == 0 || k == 0 || n == 0) {
    return;
  }

  constexpr std::size_t kTileI = 32;
  constexpr std::size_t kTileK = 32;
  constexpr std::size_t kTileJ = 64;

  for (std::size_t i0 = 0; i0 < m; i0 += kTileI) {
    const std::size_t i_max = std::min(i0 + kTileI, m);
    for (std::size_t k0 = 0; k0 < k; k0 += kTileK) {
      const std::size_t k_max = std::min(k0 + kTileK, k);
      for (std::size_t j0 = 0; j0 < n; j0 += kTileJ) {
        const std::size_t j_len = std::min(kTileJ, n - j0);
        for (std::size_t i = i0; i < i_max; ++i) {
          std::span<E> c_tile(c + i * n + j0, j_len);
          const E* a_row = a + i * k;
          for (std::size_t kk = k0; kk < k_max; ++kk) {
            const E& a_ik = a_row[kk];
            if (field.IsZero(a_ik)) {
              continue;
            }
            field.FmaArray(c_tile, std::span<const E>(b + kk * n + j0, j_len),
                           a_ik);
          }
        }
      }
    }
  }
}

template <FiniteField F>
inline void ExtractPaddedBlock(const typename F::element_type* src,
                               std::size_t src_rows, std::size_t src_cols,
                               std::size_t row_off, std::size_t col_off,
                               std::size_t blk_rows, std::size_t blk_cols,
                               const typename F::element_type& zero,
                               std::vector<typename F::element_type>& dst) {
  dst.assign(blk_rows * blk_cols, zero);
  if (row_off >= src_rows || col_off >= src_cols) {
    return;
  }
  const std::size_t copy_rows = std::min(blk_rows, src_rows - row_off);
  const std::size_t copy_cols = std::min(blk_cols, src_cols - col_off);
  for (std::size_t r = 0; r < copy_rows; ++r) {
    const auto* s_row = src + (row_off + r) * src_cols + col_off;
    auto* d_row = dst.data() + r * blk_cols;
    std::copy(s_row, s_row + copy_cols, d_row);
  }
}

template <FiniteField F>
inline void ScatterBlock(const std::vector<typename F::element_type>& src,
                         std::size_t blk_rows, std::size_t blk_cols,
                         typename F::element_type* dst, std::size_t dst_rows,
                         std::size_t dst_cols, std::size_t row_off,
                         std::size_t col_off) {
  if (row_off >= dst_rows || col_off >= dst_cols) {
    return;
  }
  const std::size_t copy_rows = std::min(blk_rows, dst_rows - row_off);
  const std::size_t copy_cols = std::min(blk_cols, dst_cols - col_off);
  for (std::size_t r = 0; r < copy_rows; ++r) {
    const auto* s_row = src.data() + r * blk_cols;
    auto* d_row = dst + (row_off + r) * dst_cols + col_off;
    std::copy(s_row, s_row + copy_cols, d_row);
  }
}

// Recursive Strassen-Winograd 7-product GEMM using the Boyer-Dumas minimal
// temporary schedule (15 additions, works in all characteristics including 2).
template <FiniteField F>
inline void StrassenWinogradRaw(const F& field,
                                const typename F::element_type* a,
                                const typename F::element_type* b,
                                typename F::element_type* c, std::size_t m,
                                std::size_t k, std::size_t n,
                                std::size_t cutoff) {
  using E = typename F::element_type;
  const std::size_t eff_cutoff = std::max<std::size_t>(cutoff, 2);
  if (m < eff_cutoff || k < eff_cutoff || n < eff_cutoff) {
    BlockedFmaGemmRaw(field, a, b, c, m, k, n);
    return;
  }

  const E zero = field.AdditiveIdentity();
  const std::size_t mh = (m + 1) / 2;
  const std::size_t kh = (k + 1) / 2;
  const std::size_t nh = (n + 1) / 2;
  const std::size_t size_a = mh * kh;
  const std::size_t size_b = kh * nh;
  const std::size_t size_c = mh * nh;

  std::vector<E> a11, a12, a21, a22;
  ExtractPaddedBlock<F>(a, m, k, 0, 0, mh, kh, zero, a11);
  ExtractPaddedBlock<F>(a, m, k, 0, kh, mh, kh, zero, a12);
  ExtractPaddedBlock<F>(a, m, k, mh, 0, mh, kh, zero, a21);
  ExtractPaddedBlock<F>(a, m, k, mh, kh, mh, kh, zero, a22);

  std::vector<E> b11, b12, b21, b22;
  ExtractPaddedBlock<F>(b, k, n, 0, 0, kh, nh, zero, b11);
  ExtractPaddedBlock<F>(b, k, n, 0, nh, kh, nh, zero, b12);
  ExtractPaddedBlock<F>(b, k, n, kh, 0, kh, nh, zero, b21);
  ExtractPaddedBlock<F>(b, k, n, kh, nh, kh, nh, zero, b22);

  std::vector<E> c11(size_c, zero);
  std::vector<E> c12(size_c, zero);
  std::vector<E> c21(size_c, zero);
  std::vector<E> c22(size_c, zero);

  // Boyer-Dumas schedule: only 3 temporaries (tmp_x, tmp_y, tmp_p).
  std::vector<E> tmp_x(size_a, zero);
  std::vector<E> tmp_y(size_b, zero);
  std::vector<E> tmp_p(size_c, zero);

  // 1. X = S3 = A11 - A21, Y = T3 = B22 - B12 -> C21 = P7 = S3 * T3
  field.SubArray(std::span<E>(tmp_x), std::span<const E>(a11),
                 std::span<const E>(a21));
  field.SubArray(std::span<E>(tmp_y), std::span<const E>(b22),
                 std::span<const E>(b12));
  StrassenWinogradRaw(field, tmp_x.data(), tmp_y.data(), c21.data(), mh, kh, nh,
                      eff_cutoff);

  // 2. X = S1 = A21 + A22, Y = T1 = B12 - B11 -> C22 = P5 = S1 * T1
  field.AddArray(std::span<E>(tmp_x), std::span<const E>(a21),
                 std::span<const E>(a22));
  field.SubArray(std::span<E>(tmp_y), std::span<const E>(b12),
                 std::span<const E>(b11));
  StrassenWinogradRaw(field, tmp_x.data(), tmp_y.data(), c22.data(), mh, kh, nh,
                      eff_cutoff);

  // 3. X = S2 = S1 - A11, Y = T2 = B22 - T1 -> C12 = P6 = S2 * T2
  field.SubArray(std::span<E>(tmp_x), std::span<const E>(tmp_x),
                 std::span<const E>(a11));
  field.SubArray(std::span<E>(tmp_y), std::span<const E>(b22),
                 std::span<const E>(tmp_y));
  StrassenWinogradRaw(field, tmp_x.data(), tmp_y.data(), c12.data(), mh, kh, nh,
                      eff_cutoff);

  // 4. X = S4 = A12 - S2 -> C11 = P3 = S4 * B22 (temporary)
  field.SubArray(std::span<E>(tmp_x), std::span<const E>(a12),
                 std::span<const E>(tmp_x));
  StrassenWinogradRaw(field, tmp_x.data(), b22.data(), c11.data(), mh, kh, nh,
                      eff_cutoff);

  // 5. P = P1 = A11 * B11; assemble U2..U5, U7 -> C12 and C22 finalized!
  StrassenWinogradRaw(field, a11.data(), b11.data(), tmp_p.data(), mh, kh, nh,
                      eff_cutoff);
  // C12 = U2 = P1 + P6
  field.AddArray(std::span<E>(c12), std::span<const E>(tmp_p),
                 std::span<const E>(c12));
  // C21 = U3 = U2 + P7
  field.AddArray(std::span<E>(c21), std::span<const E>(c12),
                 std::span<const E>(c21));
  // C12 = U4 = U2 + P5
  field.AddArray(std::span<E>(c12), std::span<const E>(c12),
                 std::span<const E>(c22));
  // C22 = U7 = U3 + P5 (final C22)
  field.AddArray(std::span<E>(c22), std::span<const E>(c21),
                 std::span<const E>(c22));
  // C12 = U5 = U4 + P3 (final C12)
  field.AddArray(std::span<E>(c12), std::span<const E>(c12),
                 std::span<const E>(c11));

  // 6. Y = T4 = T2 - B21 -> C11 = P4 = A22 * T4 -> C21 = U6 = U3 - P4 (final)
  field.SubArray(std::span<E>(tmp_y), std::span<const E>(tmp_y),
                 std::span<const E>(b21));
  StrassenWinogradRaw(field, a22.data(), tmp_y.data(), c11.data(), mh, kh, nh,
                      eff_cutoff);
  field.SubArray(std::span<E>(c21), std::span<const E>(c21),
                 std::span<const E>(c11));

  // 7. C11 = P2 = A12 * B21 -> C11 = U1 = P1 + P2 (final C11)
  StrassenWinogradRaw(field, a12.data(), b21.data(), c11.data(), mh, kh, nh,
                      eff_cutoff);
  field.AddArray(std::span<E>(c11), std::span<const E>(tmp_p),
                 std::span<const E>(c11));

  ScatterBlock<F>(c11, mh, nh, c, m, n, 0, 0);
  ScatterBlock<F>(c12, mh, nh, c, m, n, 0, nh);
  ScatterBlock<F>(c21, mh, nh, c, m, n, mh, 0);
  ScatterBlock<F>(c22, mh, nh, c, m, n, mh, nh);
}

}  // namespace detail

// ============================================================================
// Method of Four Russians (M4RI, k = 8 Gray-code LUT) over GF(2)
// ============================================================================

template <typename BinaryF = GaloisFieldBinary>
  requires std::same_as<BinaryF, GaloisFieldBinary>
inline GFMatrix<BinaryF> M4riMatMul(const GFMatrix<BinaryF>& a,
                                    const GFMatrix<BinaryF>& b) {
  using RawArray = typename GFMatrix<BinaryF>::raw_array_type;
  using ShapeType = typename GFMatrix<BinaryF>::shape_type;

  detail::RequireSameField(a, b, "M4riMatMul");
  if (a.Dimension() != 2 || b.Dimension() != 2) {
    throw std::invalid_argument("M4riMatMul requires 2-D matrices");
  }
  if (a.Shape(1) != b.Shape(0)) {
    throw std::invalid_argument(
        "M4riMatMul: incompatible matrix dimensions for multiplication");
  }

  auto field = a.Field();
  const std::size_t m = a.Shape(0);
  const std::size_t k = a.Shape(1);
  const std::size_t n = b.Shape(1);
  RawArray raw = RawArray::from_shape(ShapeType{m, n});
  if (m == 0 || n == 0) {
    return GFMatrix<BinaryF>::FromRaw(std::move(field), std::move(raw));
  }
  if (k == 0) {
    raw.fill(uint8_t{0});
    return GFMatrix<BinaryF>::FromRaw(std::move(field), std::move(raw));
  }

  const detail::PackedBinaryMatrix packed_a =
      detail::PackMatrix(a.RawData().data(), m, k);
  const detail::PackedBinaryMatrix packed_b =
      detail::PackMatrix(b.RawData().data(), k, n);
  detail::PackedBinaryMatrix packed_c = detail::PackedBinaryMatrix::Zeros(m, n);
  const std::size_t words_n = packed_b.words_per_row;

  // k = 8 Gray-code lookup table of 256 packed rows (8 divides 64, so every
  // 8-bit chunk of A lies cleanly within a single uint64_t word).
  constexpr std::size_t kStep = 8;
  std::vector<uint64_t> table((std::size_t{1} << kStep) * words_n, 0ULL);

  for (std::size_t k0 = 0; k0 < k; k0 += kStep) {
    const std::size_t b_len = std::min(kStep, k - k0);
    const std::size_t table_size = std::size_t{1} << b_len;
    std::fill(table.begin(),
              table.begin() + static_cast<std::ptrdiff_t>(words_n), 0ULL);

    // Gray-code table construction: 1 packed-row XOR per entry.
    std::size_t g_prev = 0;
    for (std::size_t t = 1; t < table_size; ++t) {
      const std::size_t g_curr = t ^ (t >> 1);
      const std::size_t bit =
          static_cast<std::size_t>(std::countr_zero(g_prev ^ g_curr));
      const uint64_t* prev_row = table.data() + g_prev * words_n;
      const uint64_t* b_row = packed_b.RowPtr(k0 + bit);
      uint64_t* curr_row = table.data() + g_curr * words_n;
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
      uint64_t* c_row = packed_c.RowPtr(i);
      const uint64_t* t_row = table.data() + idx * words_n;
      for (std::size_t w = 0; w < words_n; ++w) {
        c_row[w] ^= t_row[w];
      }
    }
  }

  detail::UnpackMatrix(packed_c, raw.data());
  return GFMatrix<BinaryF>::FromRaw(std::move(field), std::move(raw));
}

// Bit-packed M4RI Reduced Row Echelon Form over GF(2) using k = 8 Gray-code
// pivot-strip elimination. Returns {rref_matrix, pivot_cols} matching
// linalg::Rref bit-for-bit.
template <typename BinaryF = GaloisFieldBinary>
  requires std::same_as<BinaryF, GaloisFieldBinary>
inline std::pair<GFMatrix<BinaryF>, std::vector<std::size_t>> M4riRref(
    const GFMatrix<BinaryF>& a) {
  using RawArray = typename GFMatrix<BinaryF>::raw_array_type;
  using ShapeType = typename GFMatrix<BinaryF>::shape_type;

  detail::RequireValidArray(a, "M4riRref");
  if (a.Dimension() != 2) {
    throw std::invalid_argument("M4riRref requires a 2-D matrix");
  }

  auto field = a.Field();
  const std::size_t rows = a.Shape(0);
  const std::size_t cols = a.Shape(1);
  if (rows == 0 || cols == 0) {
    return {a, {}};
  }

  detail::PackedBinaryMatrix m =
      detail::PackMatrix(a.RawData().data(), rows, cols);
  const std::size_t words = m.words_per_row;
  std::vector<std::size_t> pivots;
  pivots.reserve(std::min(rows, cols));

  constexpr std::size_t kStripWidth = 8;
  std::vector<uint64_t> table((std::size_t{1} << kStripWidth) * words, 0ULL);
  std::array<uint64_t, 256> table_word{};
  std::array<std::size_t, kStripWidth> strip_pivots{};

  std::size_t pivot_row = 0;
  std::size_t col = 0;

  while (col < cols && pivot_row < rows) {
    const std::size_t w_idx = col >> 6;
    const std::size_t word_end_col = std::min(cols, (w_idx + 1) << 6);
    const std::size_t strip_end_col = std::min(word_end_col, col + kStripWidth);

    std::size_t t = 0;
    table_word[0] = 0ULL;

    auto extract_idx = [&](uint64_t word) noexcept -> std::size_t {
      std::size_t idx = 0;
      for (std::size_t b = 0; b < t; ++b) {
        idx |= static_cast<std::size_t>(
            ((word >> (strip_pivots[b] & 63)) & 1ULL) << b);
      }
      return idx;
    };

    // Phase 1: Locate up to 8 pivots in [col, strip_end_col) via lazy
    // single-word reduction and maintain local RREF among the t strip pivot
    // rows.
    for (std::size_t c = col; c < strip_end_col && pivot_row + t < rows; ++c) {
      const uint64_t c_bit = uint64_t{1} << (c & 63);
      std::size_t found = pivot_row + t;
      while (found < rows) {
        const uint64_t w_val = m.RowPtr(found)[w_idx];
        const std::size_t idx = extract_idx(w_val);
        if (((w_val ^ table_word[idx]) & c_bit) != 0ULL) {
          break;
        }
        ++found;
      }
      if (found == rows) {
        continue;
      }

      // Apply pending strip-pivot reductions to the newly selected pivot row.
      const std::size_t idx = extract_idx(m.RowPtr(found)[w_idx]);
      for (std::size_t b = 0; b < t; ++b) {
        if (((idx >> b) & std::size_t{1}) != 0) {
          m.XorRow(found, pivot_row + b, w_idx);
        }
      }

      if (found != pivot_row + t) {
        m.SwapRows(pivot_row + t, found);
      }

      // Keep the strip pivot rows in exact RREF on strip_pivots[0..t].
      for (std::size_t b = 0; b < t; ++b) {
        if ((m.RowPtr(pivot_row + b)[w_idx] & c_bit) != 0ULL) {
          m.XorRow(pivot_row + b, pivot_row + t, w_idx);
        }
      }

      strip_pivots[t] = c;
      pivots.push_back(c);
      ++t;

      // Refresh single-word lookup table for the t strip pivot rows.
      const std::size_t tw_size = std::size_t{1} << t;
      table_word[0] = 0ULL;
      for (std::size_t mask = 1; mask < tw_size; ++mask) {
        const std::size_t lsb =
            static_cast<std::size_t>(std::countr_zero(mask));
        table_word[mask] = table_word[mask ^ (std::size_t{1} << lsb)] ^
                           m.RowPtr(pivot_row + lsb)[w_idx];
      }
    }

    col = strip_end_col;
    if (t == 0) {
      continue;
    }

    // Phase 2: Build full-width Gray-code table of size 2^t from the t strip
    // pivot rows (only words [w_idx, words) are non-zero).
    const std::size_t rem_words = words - w_idx;
    const std::size_t table_size = std::size_t{1} << t;
    std::fill(table.data(), table.data() + rem_words, 0ULL);

    std::size_t g_prev = 0;
    for (std::size_t s = 1; s < table_size; ++s) {
      const std::size_t g_curr = s ^ (s >> 1);
      const std::size_t bit =
          static_cast<std::size_t>(std::countr_zero(g_prev ^ g_curr));
      const uint64_t* prev_row = table.data() + g_prev * rem_words;
      const uint64_t* piv_row_ptr = m.RowPtr(pivot_row + bit) + w_idx;
      uint64_t* curr_row = table.data() + g_curr * rem_words;
      for (std::size_t w = 0; w < rem_words; ++w) {
        curr_row[w] = prev_row[w] ^ piv_row_ptr[w];
      }
      g_prev = g_curr;
    }

    // Phase 3: Eliminate all t pivot columns from all non-strip rows in one
    // table lookup per row.
    auto eliminate_row = [&](std::size_t r) noexcept {
      uint64_t* dst = m.RowPtr(r) + w_idx;
      const std::size_t idx = extract_idx(dst[0]);
      if (idx == 0) {
        return;
      }
      const uint64_t* t_row = table.data() + idx * rem_words;
      for (std::size_t w = 0; w < rem_words; ++w) {
        dst[w] ^= t_row[w];
      }
    };

    for (std::size_t r = 0; r < pivot_row; ++r) {
      eliminate_row(r);
    }
    for (std::size_t r = pivot_row + t; r < rows; ++r) {
      eliminate_row(r);
    }

    pivot_row += t;
  }

  RawArray raw = RawArray::from_shape(ShapeType{rows, cols});
  detail::UnpackMatrix(m, raw.data());
  return {GFMatrix<BinaryF>::FromRaw(std::move(field), std::move(raw)),
          std::move(pivots)};
}

// ============================================================================
// Strassen-Winograd Matrix Multiplication & Unified FastMatMul Dispatcher
// ============================================================================

template <FiniteField F>
inline GFMatrix<F> StrassenMatMul(const GFMatrix<F>& a, const GFMatrix<F>& b,
                                  std::size_t cutoff = 64) {
  using RawArray = typename GFMatrix<F>::raw_array_type;
  using ShapeType = typename GFMatrix<F>::shape_type;

  detail::RequireSameField(a, b, "StrassenMatMul");
  if (a.Dimension() != 2 || b.Dimension() != 2) {
    throw std::invalid_argument("StrassenMatMul requires 2-D matrices");
  }
  if (a.Shape(1) != b.Shape(0)) {
    throw std::invalid_argument(
        "StrassenMatMul: incompatible matrix dimensions for multiplication");
  }

  auto field = a.Field();
  const std::size_t m = a.Shape(0);
  const std::size_t k = a.Shape(1);
  const std::size_t n = b.Shape(1);
  RawArray raw = RawArray::from_shape(ShapeType{m, n});
  if (m == 0 || n == 0) {
    return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
  }
  if (k == 0) {
    raw.fill(field->AdditiveIdentity());
    return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
  }

  detail::StrassenWinogradRaw(*field, a.RawData().data(), b.RawData().data(),
                              raw.data(), m, k, n, cutoff);
  return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
}

// Unified high-performance matrix multiplication dispatcher:
//   - GaloisFieldBinary -> M4riMatMul (bit-packed k=8 Gray-code LUT)
//   - min(m, k, n) >= 128 -> StrassenMatMul (Winograd 7-product, cutoff = 64)
//   - Medium/small matrices -> cache-blocked FmaArray GEMM
template <FiniteField F>
inline GFMatrix<F> FastMatMul(const GFMatrix<F>& a, const GFMatrix<F>& b) {
  using RawArray = typename GFMatrix<F>::raw_array_type;
  using ShapeType = typename GFMatrix<F>::shape_type;

  detail::RequireSameField(a, b, "FastMatMul");
  if (a.Dimension() != 2 || b.Dimension() != 2) {
    throw std::invalid_argument("FastMatMul requires 2-D matrices");
  }
  if (a.Shape(1) != b.Shape(0)) {
    throw std::invalid_argument(
        "FastMatMul: incompatible matrix dimensions for multiplication");
  }

  if constexpr (std::is_same_v<F, GaloisFieldBinary>) {
    return M4riMatMul(a, b);
  } else {
    const std::size_t m = a.Shape(0);
    const std::size_t k = a.Shape(1);
    const std::size_t n = b.Shape(1);
    if (std::min({m, k, n}) >= 128) {
      return StrassenMatMul(a, b, 64);
    }
    auto field = a.Field();
    RawArray raw = RawArray::from_shape(ShapeType{m, n});
    if (m == 0 || n == 0) {
      return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
    }
    if (k == 0) {
      raw.fill(field->AdditiveIdentity());
      return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
    }
    detail::BlockedFmaGemmRaw(*field, a.RawData().data(), b.RawData().data(),
                              raw.data(), m, k, n);
    return GFMatrix<F>::FromRaw(std::move(field), std::move(raw));
  }
}

}  // namespace linalg
}  // namespace xg

#endif  // XGALOIS_LINALG_FAST_MATMUL_HPP_
