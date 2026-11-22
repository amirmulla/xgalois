#ifndef XGALOIS_CODING_MODIFIERS_HPP_
#define XGALOIS_CODING_MODIFIERS_HPP_

#include <algorithm>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/linear_code.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/utils/bigint.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

namespace detail {

inline std::vector<std::size_t> ComplementPositions(
    std::size_t length, const std::vector<std::size_t> &positions,
    const char *fn_name) {
  std::vector<bool> removed(length, false);
  std::size_t removed_count = 0;
  for (std::size_t pos : positions) {
    if (pos >= length) {
      throw std::out_of_range(std::string(fn_name) +
                              ": coordinate position out of range");
    }
    if (removed[pos]) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": duplicate coordinate position");
    }
    removed[pos] = true;
    ++removed_count;
  }
  if (removed_count >= length) {
    throw std::invalid_argument(std::string(fn_name) +
                                ": cannot remove all coordinates of a code");
  }
  std::vector<std::size_t> kept;
  kept.reserve(length - removed_count);
  for (std::size_t i = 0; i < length; ++i) {
    if (!removed[i]) {
      kept.push_back(i);
    }
  }
  return kept;
}

template <typename GaloisField>
GFMatrix<GaloisField> SelectColumns(const GFMatrix<GaloisField> &mat,
                                    const std::vector<std::size_t> &cols_keep) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;

  const std::size_t rows = mat.Shape(0);
  const std::size_t cols = mat.Shape(1);
  const std::size_t new_cols = cols_keep.size();

  RawArray out = RawArray::from_shape(ShapeType{rows, new_cols});
  if (rows > 0 && new_cols > 0) {
    const StorageType *src = mat.RawData().data();
    StorageType *dst = out.data();
    for (std::size_t r = 0; r < rows; ++r) {
      const StorageType *src_row = src + r * cols;
      StorageType *dst_row = dst + r * new_cols;
      for (std::size_t j = 0; j < new_cols; ++j) {
        dst_row[j] = src_row[cols_keep[j]];
      }
    }
  }
  return GFMatrix<GaloisField>::FromRaw(mat.Field(), std::move(out));
}

inline double BoxplusLlr(double l1, double l2, bool exact_boxplus) noexcept {
  if (!exact_boxplus) {
    const double s1 = (l1 >= 0.0) ? 1.0 : -1.0;
    const double s2 = (l2 >= 0.0) ? 1.0 : -1.0;
    return s1 * s2 * std::min(std::abs(l1), std::abs(l2));
  }
  const double prod = std::tanh(0.5 * l1) * std::tanh(0.5 * l2);
  constexpr double kMaxAbsTanh = 1.0 - 1e-15;
  if (prod >= kMaxAbsTanh) {
    return std::min(std::abs(l1), std::abs(l2));
  }
  if (prod <= -kMaxAbsTanh) {
    return -std::min(std::abs(l1), std::abs(l2));
  }
  return 2.0 * std::atanh(prod);
}

// Enumerates the exact Hamming weight distribution of the row space of an
// (m x n) full-row-rank matrix M over GF(q).
template <typename GaloisField>
std::vector<BigInt> EnumerateRowSpaceWeights(
    const std::shared_ptr<GaloisField> &field, const GFMatrix<GaloisField> &mat,
    std::size_t max_enumerated_words) {
  using StorageType = typename GaloisField::element_type;

  const std::size_t m = mat.Shape(0);
  const std::size_t n = mat.Shape(1);
  const uint64_t q = field->Order();
  if (q == 0) {
    throw std::runtime_error(
        "ComputeWeightEnumerator: field order too large for enumeration");
  }

  std::vector<BigInt> counts(n + 1, BigInt(0));
  if (m == 0) {
    counts[0] = BigInt(1);
    return counts;
  }

  uint64_t total = 1;
  try {
    total = utils::SafeIntegerPower(q, m);
  } catch (const std::overflow_error &) {
    throw std::runtime_error(
        "ComputeWeightEnumerator: code size exceeds 64-bit integer range");
  }
  if (total > max_enumerated_words) {
    throw std::runtime_error(
        "ComputeWeightEnumerator: number of words to enumerate (" +
        std::to_string(total) + ") exceeds max_enumerated_words (" +
        std::to_string(max_enumerated_words) + ")");
  }

  std::vector<uint64_t> counts_u64(n + 1, 0);

  // Fast bit-packed Gray-code path for GF(2): 1 word-XOR + popcount per word.
  if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
    const std::size_t words = (n + 63) / 64;
    std::vector<uint64_t> packed_rows(m * words, 0ULL);
    const uint8_t *src = mat.RawData().data();
    for (std::size_t r = 0; r < m; ++r) {
      const uint8_t *src_row = src + r * n;
      uint64_t *dst_row = packed_rows.data() + r * words;
      for (std::size_t c = 0; c < n; ++c) {
        if ((src_row[c] & 1U) != 0) {
          dst_row[c >> 6] |= (uint64_t{1} << (c & 63));
        }
      }
    }

    std::vector<uint64_t> acc(words, 0ULL);
    counts_u64[0] = 1;  // All-zero codeword (t = 0).

    uint64_t g_prev = 0;
    for (uint64_t t = 1; t < total; ++t) {
      const uint64_t g_curr = t ^ (t >> 1);
      const uint64_t diff = g_prev ^ g_curr;
      const std::size_t bit = static_cast<std::size_t>(std::countr_zero(diff));
      const uint64_t *row_ptr = packed_rows.data() + bit * words;
      std::size_t weight = 0;
      for (std::size_t w = 0; w < words; ++w) {
        acc[w] ^= row_ptr[w];
        weight += static_cast<std::size_t>(std::popcount(acc[w]));
      }
      ++counts_u64[weight];
      g_prev = g_curr;
    }
  } else {
    // General q-ary incremental odometer: amortized ~1 FmaArray per codeword.
    const StorageType zero = field->AdditiveIdentity();
    std::vector<StorageType> alphabet;
    alphabet.reserve(static_cast<std::size_t>(q));
    alphabet.push_back(zero);
    if constexpr (std::is_integral_v<StorageType>) {
      for (uint64_t i = 0; i < q; ++i) {
        StorageType v = field->SetElementValue(static_cast<StorageType>(i));
        if (!field->IsZero(v)) {
          alphabet.push_back(v);
        }
      }
    } else {
      for (uint64_t i = 0; i < q; ++i) {
        StorageType v = field->FromInt(i);
        if (!field->IsZero(v)) {
          alphabet.push_back(std::move(v));
        }
      }
    }

    // Precompute delta[d] = alphabet[d + 1] - alphabet[d] for 0 <= d < q - 1,
    // and wrap_delta = -alphabet[q - 1].
    std::vector<StorageType> step_delta(static_cast<std::size_t>(q - 1));
    for (std::size_t d = 0; d + 1 < alphabet.size(); ++d) {
      step_delta[d] = field->Sub(alphabet[d + 1], alphabet[d]);
    }
    const StorageType wrap_delta = field->Neg(alphabet.back());

    const StorageType *mat_data = mat.RawData().data();
    std::vector<StorageType> cur(n, zero);
    std::vector<std::size_t> digits(m, 0);
    counts_u64[0] = 1;

    for (uint64_t idx = 1; idx < total; ++idx) {
      for (std::size_t j = 0; j < m; ++j) {
        std::span<const StorageType> row_span(mat_data + j * n, n);
        if (digits[j] + 1 < alphabet.size()) {
          field->FmaArray(std::span<StorageType>(cur), row_span,
                          step_delta[digits[j]]);
          ++digits[j];
          break;
        }
        field->FmaArray(std::span<StorageType>(cur), row_span, wrap_delta);
        digits[j] = 0;
      }

      std::size_t weight = 0;
      for (std::size_t c = 0; c < n; ++c) {
        if (!field->IsZero(cur[c])) {
          ++weight;
        }
      }
      ++counts_u64[weight];
    }
  }

  for (std::size_t w = 0; w <= n; ++w) {
    counts[w] = BigInt::FromUint64(counts_u64[w]);
  }
  return counts;
}

}  // namespace detail

// ============================================================================
// General Linear Code Modifiers
// ============================================================================

// PuncturedCode: deletes coordinates in `positions` from every codeword.
// G' = RowSpace(G_{:, \bar{P}}).
template <typename GaloisField>
LinearCode<GaloisField> PuncturedCode(
    const AbstractLinearCode<GaloisField> &code,
    const std::vector<std::size_t> &positions) {
  const auto kept =
      detail::ComplementPositions(code.Length(), positions, "PuncturedCode");
  const auto &G = code.GeneratorMatrix();
  auto G_sub = detail::SelectColumns(G, kept);
  auto G_basis =
      (G_sub.Shape(0) == 0) ? std::move(G_sub) : linalg::RowSpace(G_sub);
  return LinearCode<GaloisField>::FromGeneratorMatrix(code.Field(), G_basis,
                                                      code.GetMetric());
}

// ShortenedCode: restricts to codewords that are zero on `positions` and
// deletes those coordinates. By duality, H' = RowSpace(H_{:, \bar{P}}).
template <typename GaloisField>
LinearCode<GaloisField> ShortenedCode(
    const AbstractLinearCode<GaloisField> &code,
    const std::vector<std::size_t> &positions) {
  const auto kept =
      detail::ComplementPositions(code.Length(), positions, "ShortenedCode");
  const auto &H = code.ParityCheckMatrix();
  auto H_sub = detail::SelectColumns(H, kept);
  return LinearCode<GaloisField>::FromParityCheckMatrix(code.Field(), H_sub,
                                                        code.GetMetric());
}

// ExtendedCode: appends an overall parity symbol c_n = -sum_{i=0}^{n-1} c_i
// to every codeword. For binary codes with odd minimum distance d, sets the
// known minimum distance to d + 1.
template <typename GaloisField>
LinearCode<GaloisField> ExtendedCode(
    const AbstractLinearCode<GaloisField> &code) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;

  auto field = code.Field();
  const auto &G = code.GeneratorMatrix();
  const std::size_t k = G.Shape(0);
  const std::size_t n = G.Shape(1);
  const std::size_t new_n = n + 1;

  RawArray g_ext_raw = RawArray::from_shape(ShapeType{k, new_n});
  const StorageType zero = field->AdditiveIdentity();
  const StorageType *src = G.RawData().data();
  StorageType *dst = g_ext_raw.data();

  for (std::size_t r = 0; r < k; ++r) {
    const StorageType *src_row = src + r * n;
    StorageType *dst_row = dst + r * new_n;
    StorageType sum = zero;
    for (std::size_t c = 0; c < n; ++c) {
      dst_row[c] = src_row[c];
      sum = field->Add(sum, src_row[c]);
    }
    dst_row[n] = field->Neg(sum);
  }

  auto G_ext = GFMatrix<GaloisField>::FromRaw(field, std::move(g_ext_raw));

  std::optional<std::size_t> known_d = std::nullopt;
  if (field->Order() == 2) {
    try {
      const std::size_t d = code.MinimumDistance();
      if (k == 0) {
        known_d = new_n;
      } else {
        known_d = ((d & std::size_t{1}) != 0) ? (d + 1) : d;
      }
    } catch (...) {
      known_d = std::nullopt;
    }
  }

  return LinearCode<GaloisField>::FromGeneratorMatrix(
      std::move(field), G_ext, code.GetMetric(), known_d);
}

// ExpurgatedCode: intersects the code with the hyperplane defined by
// `parity_row` (default: all-ones vector, which removes all odd-weight
// codewords for binary codes).
template <typename GaloisField>
LinearCode<GaloisField> ExpurgatedCode(
    const AbstractLinearCode<GaloisField> &code,
    const std::optional<GFArray<GaloisField>> &parity_row = std::nullopt) {
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;

  auto field = code.Field();
  const std::size_t n = code.Length();

  GFMatrix<GaloisField> extra_h;
  if (parity_row.has_value()) {
    if (!parity_row->HasField() || parity_row->Field() != field) {
      throw std::invalid_argument("ExpurgatedCode: parity_row field mismatch");
    }
    if (parity_row->Dimension() != 1 || parity_row->Size() != n) {
      throw std::invalid_argument(
          "ExpurgatedCode: parity_row must be 1D with size equal to code "
          "length");
    }
    RawArray row_raw = RawArray::from_shape(ShapeType{1, n});
    std::copy(parity_row->RawData().data(), parity_row->RawData().data() + n,
              row_raw.data());
    extra_h = GFMatrix<GaloisField>::FromRaw(field, std::move(row_raw));
  } else {
    extra_h = GFMatrix<GaloisField>::Ones(field, {1, n});
  }

  const auto &H = code.ParityCheckMatrix();
  GFMatrix<GaloisField> H_exp =
      (H.Shape(0) == 0) ? std::move(extra_h) : linalg::VStack(H, extra_h);
  return LinearCode<GaloisField>::FromParityCheckMatrix(std::move(field), H_exp,
                                                        code.GetMetric());
}

// AugmentedCode: adjoins `extra_rows` (default: single all-ones row) to the
// generator matrix and extracts a full-row-rank basis via RowSpace.
template <typename GaloisField>
LinearCode<GaloisField> AugmentedCode(
    const AbstractLinearCode<GaloisField> &code,
    const std::optional<GFMatrix<GaloisField>> &extra_rows = std::nullopt) {
  auto field = code.Field();
  const std::size_t n = code.Length();

  GFMatrix<GaloisField> extra;
  if (extra_rows.has_value()) {
    if (!extra_rows->HasField() || extra_rows->Field() != field) {
      throw std::invalid_argument("AugmentedCode: extra_rows field mismatch");
    }
    if (extra_rows->Dimension() != 2 || extra_rows->Shape(1) != n) {
      throw std::invalid_argument(
          "AugmentedCode: extra_rows must be a 2D matrix with column count "
          "equal to code length");
    }
    extra = *extra_rows;
  } else {
    extra = GFMatrix<GaloisField>::Ones(field, {1, n});
  }

  const auto &G = code.GeneratorMatrix();
  GFMatrix<GaloisField> G_aug;
  if (G.Shape(0) == 0) {
    G_aug = std::move(extra);
  } else if (extra.Shape(0) == 0) {
    G_aug = G;
  } else {
    G_aug = linalg::VStack(G, extra);
  }
  GFMatrix<GaloisField> G_basis =
      (G_aug.Shape(0) == 0) ? std::move(G_aug) : linalg::RowSpace(G_aug);
  return LinearCode<GaloisField>::FromGeneratorMatrix(std::move(field), G_basis,
                                                      code.GetMetric());
}

// ============================================================================
// Plotkin (u | u + v) Construction & Hard / Soft Recursive Decoders
// ============================================================================

template <typename GaloisField>
class PlotkinCode;

template <typename GaloisField>
class PlotkinEncoder : public Encoder<GaloisField> {
 public:
  using Base = Encoder<GaloisField>;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;

  explicit PlotkinEncoder(const PlotkinCode<GaloisField> &code)
      : Base(code), plotkin_code_(&code) {}

  codeword_type Encode(const message_type &message) const override {
    return plotkin_code_->EncodePlotkin(message);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    return plotkin_code_->UnencodePlotkin(codeword);
  }

  std::string Name() const override { return "Plotkin"; }

 private:
  const PlotkinCode<GaloisField> *plotkin_code_;
};

template <typename GaloisField>
class PlotkinHardDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using codeword_type = typename Base::codeword_type;
  using result_type = typename Base::result_type;

  explicit PlotkinHardDecoder(const PlotkinCode<GaloisField> &code)
      : Base(code), plotkin_code_(&code) {}

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return plotkin_code_->DecodeHardPlotkin(received_word);
  }

  std::string Name() const override { return "PlotkinHard"; }

  std::size_t DecodingRadius() const override {
    const std::size_t d = plotkin_code_->MinimumDistance();
    return (d > 0) ? (d - 1) / 2 : 0;
  }

 private:
  const PlotkinCode<GaloisField> *plotkin_code_;
};

template <typename GaloisField>
class PlotkinCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = typename Base::matrix_type;
  using result_type = typename Base::result_type;

  PlotkinCode(std::shared_ptr<const AbstractLinearCode<GaloisField>> u_code,
              std::shared_ptr<const AbstractLinearCode<GaloisField>> v_code)
      : Base(ValidateAndGetField(u_code, v_code), 2 * u_code->Length(),
             u_code->Dimension() + v_code->Dimension(), Metric::HAMMING,
             "Plotkin", "PlotkinHard"),
        u_code_(std::move(u_code)),
        v_code_(std::move(v_code)),
        half_len_(u_code_->Length()),
        k_u_(u_code_->Dimension()),
        k_v_(v_code_->Dimension()) {
    RegisterEncodersAndDecoders();
  }

  PlotkinCode(const PlotkinCode &other)
      : Base(other),
        u_code_(other.u_code_),
        v_code_(other.v_code_),
        half_len_(other.half_len_),
        k_u_(other.k_u_),
        k_v_(other.k_v_) {
    std::lock_guard<std::mutex> lock(other.matrix_mutex_);
    cached_generator_ = other.cached_generator_;
    cached_parity_check_ = other.cached_parity_check_;
  }

  PlotkinCode &operator=(const PlotkinCode &other) {
    if (this != &other) {
      Base::operator=(other);
      u_code_ = other.u_code_;
      v_code_ = other.v_code_;
      half_len_ = other.half_len_;
      k_u_ = other.k_u_;
      k_v_ = other.k_v_;
      std::scoped_lock lock(matrix_mutex_, other.matrix_mutex_);
      cached_generator_ = other.cached_generator_;
      cached_parity_check_ = other.cached_parity_check_;
    }
    return *this;
  }

  PlotkinCode(PlotkinCode &&other) noexcept
      : Base(std::move(other)),
        u_code_(std::move(other.u_code_)),
        v_code_(std::move(other.v_code_)),
        half_len_(other.half_len_),
        k_u_(other.k_u_),
        k_v_(other.k_v_),
        cached_generator_(std::move(other.cached_generator_)),
        cached_parity_check_(std::move(other.cached_parity_check_)) {}

  PlotkinCode &operator=(PlotkinCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      u_code_ = std::move(other.u_code_);
      v_code_ = std::move(other.v_code_);
      half_len_ = other.half_len_;
      k_u_ = other.k_u_;
      k_v_ = other.k_v_;
      std::lock_guard<std::mutex> lock(matrix_mutex_);
      cached_generator_ = std::move(other.cached_generator_);
      cached_parity_check_ = std::move(other.cached_parity_check_);
    }
    return *this;
  }

  ~PlotkinCode() override = default;

  std::shared_ptr<const AbstractLinearCode<GaloisField>> UCode() const {
    return u_code_;
  }
  std::shared_ptr<const AbstractLinearCode<GaloisField>> VCode() const {
    return v_code_;
  }
  std::size_t HalfLength() const noexcept { return half_len_; }

  std::size_t MinimumDistance() const override {
    if (this->cached_min_distance_.has_value()) {
      return *this->cached_min_distance_;
    }
    std::size_t d = this->length_;
    if (k_u_ == 0 && k_v_ == 0) {
      d = this->length_;
    } else if (k_u_ == 0) {
      d = v_code_->MinimumDistance();
    } else if (k_v_ == 0) {
      d = 2 * u_code_->MinimumDistance();
    } else {
      d = std::min(2 * u_code_->MinimumDistance(), v_code_->MinimumDistance());
    }
    this->cached_min_distance_ = d;
    return d;
  }

  // G = [ G_u  G_u ]
  //     [  0   G_v ]
  const matrix_type &GeneratorMatrix() const override {
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!cached_generator_.has_value()) {
      const auto &Gu = u_code_->GeneratorMatrix();
      const auto &Gv = v_code_->GeneratorMatrix();
      const std::size_t n = half_len_;
      const std::size_t k_total = k_u_ + k_v_;

      using RawArray = typename matrix_type::raw_array_type;
      using ShapeType = typename matrix_type::shape_type;
      RawArray g_raw = RawArray::from_shape(ShapeType{k_total, 2 * n});
      g_raw.fill(this->field_->AdditiveIdentity());
      storage_type *dst = g_raw.data();

      if (k_u_ > 0) {
        const storage_type *gu_data = Gu.RawData().data();
        for (std::size_t r = 0; r < k_u_; ++r) {
          const storage_type *src_row = gu_data + r * n;
          storage_type *dst_row = dst + r * (2 * n);
          std::copy(src_row, src_row + n, dst_row);
          std::copy(src_row, src_row + n, dst_row + n);
        }
      }
      if (k_v_ > 0) {
        const storage_type *gv_data = Gv.RawData().data();
        for (std::size_t r = 0; r < k_v_; ++r) {
          const storage_type *src_row = gv_data + r * n;
          storage_type *dst_row = dst + (k_u_ + r) * (2 * n);
          std::copy(src_row, src_row + n, dst_row + n);
        }
      }
      cached_generator_ = matrix_type::FromRaw(this->field_, std::move(g_raw));
    }
    return *cached_generator_;
  }

  // Exact parity-check matrix:
  // H = [  H_u   0  ]
  //     [ -H_v  H_v ]
  const matrix_type &ParityCheckMatrix() const override {
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!cached_parity_check_.has_value()) {
      const auto &Hu = u_code_->ParityCheckMatrix();
      const auto &Hv = v_code_->ParityCheckMatrix();
      const std::size_t n = half_len_;
      const std::size_t r_u = Hu.Shape(0);
      const std::size_t r_v = Hv.Shape(0);

      using RawArray = typename matrix_type::raw_array_type;
      using ShapeType = typename matrix_type::shape_type;
      RawArray h_raw = RawArray::from_shape(ShapeType{r_u + r_v, 2 * n});
      h_raw.fill(this->field_->AdditiveIdentity());
      storage_type *dst = h_raw.data();

      if (r_u > 0) {
        const storage_type *hu_data = Hu.RawData().data();
        for (std::size_t r = 0; r < r_u; ++r) {
          const storage_type *src_row = hu_data + r * n;
          storage_type *dst_row = dst + r * (2 * n);
          std::copy(src_row, src_row + n, dst_row);
        }
      }
      if (r_v > 0) {
        const storage_type *hv_data = Hv.RawData().data();
        for (std::size_t r = 0; r < r_v; ++r) {
          const storage_type *src_row = hv_data + r * n;
          storage_type *dst_row = dst + (r_u + r) * (2 * n);
          for (std::size_t c = 0; c < n; ++c) {
            dst_row[c] = this->field_->Neg(src_row[c]);
            dst_row[n + c] = src_row[c];
          }
        }
      }
      cached_parity_check_ =
          matrix_type::FromRaw(this->field_, std::move(h_raw));
    }
    return *cached_parity_check_;
  }

  bool Contains(const codeword_type &word) const override {
    if (!word.HasField() || word.Field() != this->field_ ||
        word.Dimension() != 1 || word.Size() != this->length_) {
      return false;
    }
    const auto [c1, c2] = SplitHalves(word);
    if (!u_code_->Contains(c1)) {
      return false;
    }
    return v_code_->Contains(c2 - c1);
  }

  codeword_type EncodePlotkin(const message_type &message) const {
    if (!message.HasField() || message.Field() != this->field_) {
      throw std::invalid_argument(
          "PlotkinCode::EncodePlotkin: message field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != this->dimension_) {
      throw std::invalid_argument(
          "PlotkinCode::EncodePlotkin: message must be 1D with length "
          "k_u + k_v");
    }
    const auto msg_u = Slice(message, 0, k_u_);
    const auto msg_v = Slice(message, k_u_, k_v_);
    const auto u = (k_u_ == 0) ? codeword_type::Zeros(this->field_, {half_len_})
                               : u_code_->Encode(msg_u);
    const auto v = (k_v_ == 0) ? codeword_type::Zeros(this->field_, {half_len_})
                               : v_code_->Encode(msg_v);
    return linalg::HStack(u, u + v);
  }

  message_type UnencodePlotkin(const codeword_type &codeword) const {
    if (!codeword.HasField() || codeword.Field() != this->field_) {
      throw std::invalid_argument(
          "PlotkinCode::UnencodePlotkin: codeword field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != this->length_) {
      throw std::invalid_argument(
          "PlotkinCode::UnencodePlotkin: codeword must be 1D with length 2*n");
    }
    const auto [c1, c2] = SplitHalves(codeword);
    const auto v = c2 - c1;
    if (!u_code_->Contains(c1) || !v_code_->Contains(v)) {
      throw std::invalid_argument(
          "PlotkinCode::UnencodePlotkin: input is not a valid Plotkin "
          "codeword");
    }
    const auto msg_u = (k_u_ == 0) ? message_type::Zeros(this->field_, {0})
                                   : u_code_->Unencode(c1);
    const auto msg_v = (k_v_ == 0) ? message_type::Zeros(this->field_, {0})
                                   : v_code_->Unencode(v);
    return linalg::HStack(msg_u, msg_v);
  }

  // Recursive hard Plotkin decoder over arbitrary finite fields:
  //   1. Decode v_hat from r2 - r1.
  //   2. Decode u_hat from both r1 and r2 - v_hat, picking the candidate
  //      minimizing Hamming distance to received word r = (r1 | r2).
  result_type DecodeHardPlotkin(const codeword_type &received_word) const {
    if (!received_word.HasField() || received_word.Field() != this->field_) {
      throw std::invalid_argument(
          "PlotkinCode::DecodeHardPlotkin: field mismatch");
    }
    if (received_word.Dimension() != 1 ||
        received_word.Size() != this->length_) {
      throw std::invalid_argument(
          "PlotkinCode::DecodeHardPlotkin: received word must be 1D with "
          "length 2*n");
    }

    const auto [r1, r2] = SplitHalves(received_word);
    const auto y_v = r2 - r1;

    codeword_type v_hat;
    message_type msg_v;
    if (k_v_ == 0) {
      v_hat = codeword_type::Zeros(this->field_, {half_len_});
      msg_v = message_type::Zeros(this->field_, {0});
    } else {
      auto res_v = v_code_->DecodeWithResult(y_v);
      if (!res_v.Ok()) {
        return result_type::Failure(
            "PlotkinCode::DecodeHardPlotkin: v-code decoding failed: " +
            res_v.error_message);
      }
      v_hat = std::move(res_v.codeword);
      msg_v = (res_v.message.HasField() && res_v.message.Size() == k_v_)
                  ? std::move(res_v.message)
                  : v_code_->Unencode(v_hat);
    }

    if (k_u_ == 0) {
      auto u_hat = codeword_type::Zeros(this->field_, {half_len_});
      auto full_cw = linalg::HStack(u_hat, v_hat);
      const std::size_t dist = this->Distance(received_word, full_cw);
      return result_type::Success(std::move(full_cw), std::move(msg_v), dist);
    }

    const auto u_obs2 = r2 - v_hat;
    auto res_u1 = u_code_->DecodeWithResult(r1);
    auto res_u2 = u_code_->DecodeWithResult(u_obs2);

    bool found = false;
    std::size_t best_dist = std::numeric_limits<std::size_t>::max();
    codeword_type best_cw;
    message_type best_msg;

    auto consider = [&](result_type &res_u) {
      if (!res_u.Ok()) {
        return;
      }
      auto cand_cw = linalg::HStack(res_u.codeword, res_u.codeword + v_hat);
      const std::size_t dist = this->Distance(received_word, cand_cw);
      if (!found || dist < best_dist) {
        found = true;
        best_dist = dist;
        auto msg_u = (res_u.message.HasField() && res_u.message.Size() == k_u_)
                         ? std::move(res_u.message)
                         : u_code_->Unencode(res_u.codeword);
        best_cw = std::move(cand_cw);
        best_msg = linalg::HStack(msg_u, msg_v);
      }
    };

    consider(res_u1);
    consider(res_u2);

    if (!found) {
      return result_type::Failure(
          "PlotkinCode::DecodeHardPlotkin: u-code decoding failed on both "
          "halves");
    }

    return result_type::Success(std::move(best_cw), std::move(best_msg),
                                best_dist);
  }

  // Recursive soft-decision Plotkin LLR decoder for binary Plotkin codes
  // (e.g. Reed-Muller codes).
  // BPSK LLR convention: L > 0 favors bit 0, L < 0 favors bit 1.
  //   1. L_v = L_1 [+] L_2 (exact 2*atanh(tanh(L1/2)*tanh(L2/2)) or min-sum).
  //   2. Soft-decode v_hat from L_v.
  //   3. L_u = L_1 + (1 - 2 * v_hat) * L_2.
  //   4. Soft-decode u_hat from L_u.
  result_type DecodeSoftPlotkin(const std::vector<double> &llrs,
                                bool exact_boxplus = true) const
    requires std::same_as<GaloisField, GaloisFieldBinary>
  {
    if (llrs.size() != this->length_) {
      throw std::invalid_argument(
          "PlotkinCode::DecodeSoftPlotkin: LLR vector must have length 2*n");
    }

    // Step 1: combine LLRs for v.
    std::vector<double> llr_v(half_len_);
    for (std::size_t i = 0; i < half_len_; ++i) {
      llr_v[i] =
          detail::BoxplusLlr(llrs[i], llrs[half_len_ + i], exact_boxplus);
    }

    // Step 2: decode v_hat.
    auto res_v = SoftDecodeBinaryComponent(*v_code_, llr_v, exact_boxplus);
    if (!res_v.Ok()) {
      return result_type::Failure(
          "PlotkinCode::DecodeSoftPlotkin: soft v-code decoding failed: " +
          res_v.error_message);
    }
    const auto &v_hat = res_v.codeword;
    auto msg_v = (res_v.message.HasField() && res_v.message.Size() == k_v_)
                     ? std::move(res_v.message)
                     : ((k_v_ == 0) ? message_type::Zeros(this->field_, {0})
                                    : v_code_->Unencode(v_hat));

    // Step 3: combine LLRs for u using decoded v_hat: L_u = L1 + (1 - 2*v)*L2.
    std::vector<double> llr_u(half_len_);
    const auto v_span = v_hat.RawSpan();
    for (std::size_t i = 0; i < half_len_; ++i) {
      const double sign = (v_span[i] == 0) ? 1.0 : -1.0;
      llr_u[i] = llrs[i] + sign * llrs[half_len_ + i];
    }

    // Step 4: decode u_hat.
    auto res_u = SoftDecodeBinaryComponent(*u_code_, llr_u, exact_boxplus);
    if (!res_u.Ok()) {
      return result_type::Failure(
          "PlotkinCode::DecodeSoftPlotkin: soft u-code decoding failed: " +
          res_u.error_message);
    }
    const auto &u_hat = res_u.codeword;
    auto msg_u = (res_u.message.HasField() && res_u.message.Size() == k_u_)
                     ? std::move(res_u.message)
                     : ((k_u_ == 0) ? message_type::Zeros(this->field_, {0})
                                    : u_code_->Unencode(u_hat));

    auto full_cw = linalg::HStack(u_hat, u_hat + v_hat);
    auto full_msg = linalg::HStack(msg_u, msg_v);

    // Count hard-decision disagreements as num_errors.
    std::size_t hard_errs = 0;
    const auto cw_span = full_cw.RawSpan();
    for (std::size_t i = 0; i < this->length_; ++i) {
      const uint8_t hard_bit = (llrs[i] < 0.0) ? uint8_t{1} : uint8_t{0};
      if (cw_span[i] != hard_bit) {
        ++hard_errs;
      }
    }

    return result_type::Success(std::move(full_cw), std::move(full_msg),
                                hard_errs);
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + ", " +
           std::to_string(MinimumDistance()) + "] Plotkin (u|u+v) code over " +
           "GF(" + std::to_string(this->field_->Order()) + ")";
  }

 private:
  static std::shared_ptr<GaloisField> ValidateAndGetField(
      const std::shared_ptr<const AbstractLinearCode<GaloisField>> &u_code,
      const std::shared_ptr<const AbstractLinearCode<GaloisField>> &v_code) {
    if (!u_code || !v_code) {
      throw std::invalid_argument(
          "PlotkinCode: u_code and v_code cannot be null");
    }
    if (u_code->Field() != v_code->Field()) {
      throw std::invalid_argument(
          "PlotkinCode: u_code and v_code must be defined over the same field");
    }
    if (u_code->Length() != v_code->Length()) {
      throw std::invalid_argument(
          "PlotkinCode: u_code and v_code must have the same length");
    }
    return u_code->Field();
  }

  void RegisterEncodersAndDecoders() {
    this->RegisterEncoder("Plotkin", [](const AbstractCode<GaloisField> &code) {
      const auto &pc = dynamic_cast<const PlotkinCode<GaloisField> &>(code);
      return std::make_shared<PlotkinEncoder<GaloisField>>(pc);
    });
    this->RegisterEncoder(
        "GeneratorMatrix", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(lin);
        });
    this->RegisterDecoder(
        "PlotkinHard", [](const AbstractCode<GaloisField> &code) {
          const auto &pc = dynamic_cast<const PlotkinCode<GaloisField> &>(code);
          return std::make_shared<PlotkinHardDecoder<GaloisField>>(pc);
        });
    this->RegisterDecoder(
        "Syndrome", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<SyndromeDecoder<GaloisField>>(lin);
        });
  }

  codeword_type Slice(const codeword_type &vec, std::size_t offset,
                      std::size_t len) const {
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{len});
    if (len > 0) {
      const storage_type *src = vec.RawData().data() + offset;
      std::copy(src, src + len, raw.data());
    }
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  std::pair<codeword_type, codeword_type> SplitHalves(
      const codeword_type &vec) const {
    return {Slice(vec, 0, half_len_), Slice(vec, half_len_, half_len_)};
  }

  // Soft-decision decoder for binary component codes:
  //   - Recurses if component is itself a PlotkinCode<GaloisFieldBinary>.
  //   - Exact ML soft correlation decoding when k <= 12 (handles small leaf
  //     codes like RM(0, m), RM(1, m), Hamming, Extended Hamming).
  //   - Hard-decision algebraic decoding fallback for larger k.
  static result_type SoftDecodeBinaryComponent(
      const AbstractLinearCode<GaloisFieldBinary> &comp,
      const std::vector<double> &llrs, bool exact_boxplus)
    requires std::same_as<GaloisField, GaloisFieldBinary>
  {
    if (const auto *nested =
            dynamic_cast<const PlotkinCode<GaloisFieldBinary> *>(&comp)) {
      return nested->DecodeSoftPlotkin(llrs, exact_boxplus);
    }

    const std::size_t n = comp.Length();
    const std::size_t k = comp.Dimension();
    auto field = comp.Field();

    if (k == 0) {
      auto cw = codeword_type::Zeros(field, {n});
      auto msg = message_type::Zeros(field, {0});
      return result_type::Success(std::move(cw), std::move(msg), 0);
    }

    if (k <= 12) {
      const auto &G = comp.GeneratorMatrix();
      const uint8_t *g_data = G.RawData().data();
      const std::size_t total = std::size_t{1} << k;

      std::vector<uint8_t> cur_cw(n, 0);
      std::vector<uint8_t> best_cw(n, 0);
      double best_corr = -std::numeric_limits<double>::infinity();

      auto eval_corr = [&](const std::vector<uint8_t> &cw) -> double {
        double corr = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
          corr += (cw[j] == 0) ? llrs[j] : -llrs[j];
        }
        return corr;
      };

      best_corr = eval_corr(cur_cw);
      std::size_t g_prev = 0;
      for (std::size_t t = 1; t < total; ++t) {
        const std::size_t g_curr = t ^ (t >> 1);
        const std::size_t bit =
            static_cast<std::size_t>(std::countr_zero(g_prev ^ g_curr));
        const uint8_t *row = g_data + bit * n;
        for (std::size_t j = 0; j < n; ++j) {
          cur_cw[j] ^= row[j];
        }
        const double corr = eval_corr(cur_cw);
        if (corr > best_corr) {
          best_corr = corr;
          best_cw = cur_cw;
        }
        g_prev = g_curr;
      }

      using RawArray = typename codeword_type::raw_array_type;
      using ShapeType = typename codeword_type::shape_type;
      RawArray cw_raw = RawArray::from_shape(ShapeType{n});
      std::copy(best_cw.begin(), best_cw.end(), cw_raw.data());
      auto cw = codeword_type::FromRaw(field, std::move(cw_raw));
      auto msg = comp.Unencode(cw);
      return result_type::Success(std::move(cw), std::move(msg), 0);
    }

    // Fallback: hard decision + component algebraic decoder.
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray hard_raw = RawArray::from_shape(ShapeType{n});
    uint8_t *h_data = hard_raw.data();
    for (std::size_t j = 0; j < n; ++j) {
      h_data[j] = (llrs[j] < 0.0) ? uint8_t{1} : uint8_t{0};
    }
    auto hard_word = codeword_type::FromRaw(field, std::move(hard_raw));
    return comp.DecodeWithResult(hard_word);
  }

  std::shared_ptr<const AbstractLinearCode<GaloisField>> u_code_;
  std::shared_ptr<const AbstractLinearCode<GaloisField>> v_code_;
  std::size_t half_len_;
  std::size_t k_u_;
  std::size_t k_v_;

  mutable std::mutex matrix_mutex_;
  mutable std::optional<matrix_type> cached_generator_;
  mutable std::optional<matrix_type> cached_parity_check_;
};

// ============================================================================
// Exact Arbitrary-Precision Weight Enumerator & MacWilliams Identity
// ============================================================================

struct WeightDistribution {
  std::size_t length = 0;
  std::size_t dimension = 0;
  uint64_t field_order = 0;
  std::vector<BigInt> counts;

  std::size_t MinimumDistance() const {
    for (std::size_t w = 1; w < counts.size(); ++w) {
      if (counts[w].IsPositive()) {
        return w;
      }
    }
    return length;
  }

  std::string ToPolynomialString(std::string_view var = "z") const {
    std::string out;
    bool first = true;
    for (std::size_t w = 0; w < counts.size(); ++w) {
      if (counts[w].IsZero()) {
        continue;
      }
      if (!first) {
        out += " + ";
      }
      first = false;
      if (w == 0) {
        out += counts[w].ToDecimal();
      } else {
        if (!counts[w].IsOne()) {
          out += counts[w].ToDecimal();
          out.push_back('*');
        }
        out.append(var);
        if (w > 1) {
          out.push_back('^');
          out += std::to_string(w);
        }
      }
    }
    return first ? "0" : out;
  }

  bool operator==(const WeightDistribution &other) const = default;
};

// Exact arbitrary-precision MacWilliams identity via Krawtchouk polynomials:
//   A_w = q^{-(n - k)} * sum_{j=0}^{n} A_j^\perp * K_w(j; n, q),
// where dual_weights has dimension k^\perp = n - k and
//   K_w(j; n, q) = sum_{r=0}^{w} (-1)^r * (q - 1)^{w - r} * C(j, r) * C(n - j,
//   w - r).
// Verifies exact divisibility by q^{k^\perp} (zero remainder).
inline WeightDistribution MacWilliamsTransform(
    const WeightDistribution &dual_weights) {
  const std::size_t n = dual_weights.length;
  const std::size_t k_dual = dual_weights.dimension;
  const uint64_t q = dual_weights.field_order;

  if (q < 2) {
    throw std::invalid_argument(
        "MacWilliamsTransform: field_order must be at least 2");
  }
  if (k_dual > n) {
    throw std::invalid_argument(
        "MacWilliamsTransform: dual dimension cannot exceed code length");
  }
  if (dual_weights.counts.size() != n + 1) {
    throw std::invalid_argument(
        "MacWilliamsTransform: counts must have size length + 1");
  }

  // Precompute Pascal's triangle C(a, b) for 0 <= a <= n.
  std::vector<std::vector<BigInt>> binom(n + 1);
  for (std::size_t a = 0; a <= n; ++a) {
    binom[a].assign(a + 1, BigInt(1));
    for (std::size_t b = 1; b < a; ++b) {
      binom[a][b] = binom[a - 1][b - 1] + binom[a - 1][b];
    }
  }

  // Precompute powers (q - 1)^m for 0 <= m <= n.
  const BigInt qm1 = BigInt::FromUint64(q - 1);
  std::vector<BigInt> qm1_pow(n + 1, BigInt(1));
  for (std::size_t m = 1; m <= n; ++m) {
    qm1_pow[m] = qm1_pow[m - 1] * qm1;
  }

  // Collect non-zero dual weights to avoid unnecessary Krawtchouk evaluations.
  std::vector<std::size_t> nonzero_j;
  nonzero_j.reserve(n + 1);
  for (std::size_t j = 0; j <= n; ++j) {
    if (!dual_weights.counts[j].IsZero()) {
      nonzero_j.push_back(j);
    }
  }

  const BigInt divisor = BigInt::FromUint64(q).Pow(k_dual);

  WeightDistribution result;
  result.length = n;
  result.dimension = n - k_dual;
  result.field_order = q;
  result.counts.assign(n + 1, BigInt(0));

  for (std::size_t w = 0; w <= n; ++w) {
    BigInt acc(0);
    for (const std::size_t j : nonzero_j) {
      // K_w(j; n, q) = sum_{r = r_min}^{r_max} (-1)^r (q-1)^{w-r} C(j, r)
      // C(n-j, w-r)
      const std::size_t r_min = (w > n - j) ? (w - (n - j)) : 0;
      const std::size_t r_max = std::min(w, j);
      BigInt kraw(0);
      for (std::size_t r = r_min; r <= r_max; ++r) {
        BigInt term = qm1_pow[w - r] * binom[j][r] * binom[n - j][w - r];
        if ((r & std::size_t{1}) != 0) {
          kraw -= term;
        } else {
          kraw += term;
        }
      }
      acc += dual_weights.counts[j] * kraw;
    }

    auto [quot, rem] = BigInt::DivRem(acc, divisor);
    if (!rem.IsZero() || quot.IsNegative()) {
      throw std::runtime_error(
          "MacWilliamsTransform: non-exact division or negative weight count "
          "(input is not a valid linear code weight distribution)");
    }
    result.counts[w] = std::move(quot);
  }

  return result;
}

// Computes the exact Hamming weight enumerator of any linear code.
// Automatically enumerates C when k <= n - k, or enumerates C^\perp and
// applies the exact arbitrary-precision MacWilliamsTransform when k > n - k.
template <typename GaloisField>
WeightDistribution ComputeWeightEnumerator(
    const AbstractLinearCode<GaloisField> &code,
    std::size_t max_enumerated_words = 1'000'000) {
  const std::size_t n = code.Length();
  const std::size_t k = code.Dimension();
  const std::size_t r = n - k;
  const uint64_t q = code.Field()->Order();

  if (k <= r) {
    WeightDistribution dist;
    dist.length = n;
    dist.dimension = k;
    dist.field_order = q;
    dist.counts = detail::EnumerateRowSpaceWeights(
        code.Field(), code.GeneratorMatrix(), max_enumerated_words);
    return dist;
  }

  WeightDistribution dual_dist;
  dual_dist.length = n;
  dual_dist.dimension = r;
  dual_dist.field_order = q;
  dual_dist.counts = detail::EnumerateRowSpaceWeights(
      code.Field(), code.ParityCheckMatrix(), max_enumerated_words);
  return MacWilliamsTransform(dual_dist);
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_MODIFIERS_HPP_
