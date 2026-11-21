#ifndef XGALOIS_CODING_LDPC_HPP_
#define XGALOIS_CODING_LDPC_HPP_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

namespace detail {

// Numerically exact Jacobian-logarithm box-plus for log-likelihood ratios:
//   BoxPlus(a, b) = ln((1 + exp(a + b)) / (exp(a) + exp(b)))
//                 = sign(a)*sign(b)*min(|a|, |b|)
//                   + ln(1 + exp(-|a + b|)) - ln(1 + exp(-|a - b|)).
// Every exponent is <= 0, avoiding overflow, underflow, and 0/0 singularities.
inline double BoxPlusLlr(double a, double b) noexcept {
  const double abs_a = std::abs(a);
  const double abs_b = std::abs(b);
  const double min_mag = (abs_a < abs_b) ? abs_a : abs_b;
  const double sign = ((a >= 0.0) == (b >= 0.0)) ? 1.0 : -1.0;
  const double corr_pos = std::log1p(std::exp(-std::abs(a + b)));
  const double corr_neg = std::log1p(std::exp(-std::abs(a - b)));
  return sign * min_mag + corr_pos - corr_neg;
}

inline void ValidateLlrVector(const std::vector<double> &channel_llrs,
                              std::size_t expected_length,
                              const char *fn_name) {
  if (channel_llrs.size() != expected_length) {
    throw std::invalid_argument(
        std::string(fn_name) + ": channel LLR length (" +
        std::to_string(channel_llrs.size()) + ") does not match code length (" +
        std::to_string(expected_length) + ")");
  }
  for (std::size_t i = 0; i < expected_length; ++i) {
    if (std::isnan(channel_llrs[i])) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": channel LLRs must not contain NaN");
    }
  }
}

}  // namespace detail

// Lightweight CSR (check-node -> variable-node) + CSC (variable-node ->
// check-node) Tanner graph representation of a binary parity-check matrix H.
class SparseBinaryMatrix {
 public:
  SparseBinaryMatrix() = default;

  SparseBinaryMatrix(
      std::size_t num_rows, std::size_t num_cols,
      const std::vector<std::pair<std::size_t, std::size_t>> &edges) {
    BuildFromEdges(num_rows, num_cols, edges);
  }

  explicit SparseBinaryMatrix(const GFMatrix<GaloisFieldBinary> &dense) {
    *this = FromDense(dense);
  }

  static SparseBinaryMatrix FromEdges(
      std::size_t num_rows, std::size_t num_cols,
      const std::vector<std::pair<std::size_t, std::size_t>> &edges) {
    return SparseBinaryMatrix(num_rows, num_cols, edges);
  }

  static SparseBinaryMatrix FromDense(
      const GFMatrix<GaloisFieldBinary> &dense) {
    if (!dense.HasField() || dense.Dimension() != 2) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::FromDense: input must be a 2D GF(2) matrix");
    }
    const std::size_t rows = dense.Shape(0);
    const std::size_t cols = dense.Shape(1);
    const uint8_t *raw = dense.RawData().data();
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    for (std::size_t r = 0; r < rows; ++r) {
      const uint8_t *row_ptr = raw + r * cols;
      for (std::size_t c = 0; c < cols; ++c) {
        if ((row_ptr[c] & 1U) != 0) {
          edges.emplace_back(r, c);
        }
      }
    }
    return SparseBinaryMatrix(rows, cols, edges);
  }

  // Construct a quasi-cyclic LDPC (QC-LDPC) parity-check matrix by lifting a
  // base shift matrix with circulant size `lifting_factor_z`.
  // Each entry `shift < 0` expands to the Z x Z zero block; `shift >= 0`
  // expands to the Z x Z circulant right-shifted identity matrix where row r
  // has a 1 at column (r + shift) mod Z.
  static SparseBinaryMatrix FromQuasiCyclic(
      const std::vector<std::vector<int>> &base_matrix,
      std::size_t lifting_factor_z) {
    if (lifting_factor_z == 0) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::FromQuasiCyclic: lifting_factor_z must be "
          "positive");
    }
    if (base_matrix.empty() || base_matrix[0].empty()) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::FromQuasiCyclic: base_matrix must be non-empty");
    }
    const std::size_t base_rows = base_matrix.size();
    const std::size_t base_cols = base_matrix[0].size();
    for (std::size_t bi = 1; bi < base_rows; ++bi) {
      if (base_matrix[bi].size() != base_cols) {
        throw std::invalid_argument(
            "SparseBinaryMatrix::FromQuasiCyclic: base_matrix rows must have "
            "uniform length");
      }
    }

    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(base_rows * base_cols * lifting_factor_z);
    for (std::size_t bi = 0; bi < base_rows; ++bi) {
      const std::size_t row_base = bi * lifting_factor_z;
      for (std::size_t bj = 0; bj < base_cols; ++bj) {
        const int shift = base_matrix[bi][bj];
        if (shift < 0) {
          continue;
        }
        const std::size_t col_base = bj * lifting_factor_z;
        const std::size_t s =
            static_cast<std::size_t>(shift) % lifting_factor_z;
        for (std::size_t r = 0; r < lifting_factor_z; ++r) {
          const std::size_t c = (r + s) % lifting_factor_z;
          edges.emplace_back(row_base + r, col_base + c);
        }
      }
    }
    return SparseBinaryMatrix(base_rows * lifting_factor_z,
                              base_cols * lifting_factor_z, edges);
  }

  std::size_t NumRows() const noexcept { return num_rows_; }
  std::size_t NumCols() const noexcept { return num_cols_; }
  std::size_t NumNonZeros() const noexcept { return row_col_idx_.size(); }

  std::size_t RowDegree(std::size_t row) const {
    if (row >= num_rows_) {
      throw std::out_of_range("SparseBinaryMatrix::RowDegree: row index");
    }
    return row_ptr_[row + 1] - row_ptr_[row];
  }

  std::size_t ColDegree(std::size_t col) const {
    if (col >= num_cols_) {
      throw std::out_of_range("SparseBinaryMatrix::ColDegree: col index");
    }
    return col_ptr_[col + 1] - col_ptr_[col];
  }

  std::span<const std::size_t> RowNeighbors(std::size_t row) const {
    if (row >= num_rows_) {
      throw std::out_of_range("SparseBinaryMatrix::RowNeighbors: row index");
    }
    const std::size_t start = row_ptr_[row];
    const std::size_t len = row_ptr_[row + 1] - start;
    return std::span<const std::size_t>(row_col_idx_.data() + start, len);
  }

  std::span<const std::size_t> ColNeighbors(std::size_t col) const {
    if (col >= num_cols_) {
      throw std::out_of_range("SparseBinaryMatrix::ColNeighbors: col index");
    }
    const std::size_t start = col_ptr_[col];
    const std::size_t len = col_ptr_[col + 1] - start;
    return std::span<const std::size_t>(col_row_idx_.data() + start, len);
  }

  const std::vector<std::size_t> &RowPtr() const noexcept { return row_ptr_; }
  const std::vector<std::size_t> &RowColIdx() const noexcept {
    return row_col_idx_;
  }
  const std::vector<std::size_t> &ColPtr() const noexcept { return col_ptr_; }
  const std::vector<std::size_t> &ColRowIdx() const noexcept {
    return col_row_idx_;
  }
  const std::vector<std::size_t> &ColToCsrEdge() const noexcept {
    return col_to_csr_edge_;
  }
  const std::vector<std::size_t> &CsrToCscEdge() const noexcept {
    return csr_to_csc_edge_;
  }

  // Compute syndrome s = H * word^T in O(nnz).
  GFArray<GaloisFieldBinary> Syndrome(
      const GFArray<GaloisFieldBinary> &word) const {
    if (!word.HasField()) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::Syndrome: word has null field");
    }
    if (word.Dimension() != 1 || word.Size() != num_cols_) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::Syndrome: word must be 1D with length equal to "
          "NumCols()");
    }
    using RawArray = GFArray<GaloisFieldBinary>::raw_array_type;
    using ShapeType = GFArray<GaloisFieldBinary>::shape_type;
    RawArray syn_raw = RawArray::from_shape(ShapeType{num_rows_});
    const uint8_t *w = word.RawData().data();
    uint8_t *s = syn_raw.data();
    for (std::size_t i = 0; i < num_rows_; ++i) {
      uint8_t parity = 0;
      const std::size_t end = row_ptr_[i + 1];
      for (std::size_t e = row_ptr_[i]; e < end; ++e) {
        parity = static_cast<uint8_t>(parity ^ (w[row_col_idx_[e]] & 1U));
      }
      s[i] = parity;
    }
    return GFArray<GaloisFieldBinary>::FromRaw(word.Field(),
                                               std::move(syn_raw));
  }

  bool CheckParity(std::span<const uint8_t> bits) const {
    if (bits.size() != num_cols_) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::CheckParity: bit vector length mismatch");
    }
    for (std::size_t i = 0; i < num_rows_; ++i) {
      uint8_t parity = 0;
      const std::size_t end = row_ptr_[i + 1];
      for (std::size_t e = row_ptr_[i]; e < end; ++e) {
        parity = static_cast<uint8_t>(parity ^ (bits[row_col_idx_[e]] & 1U));
      }
      if (parity != 0) {
        return false;
      }
    }
    return true;
  }

  bool CheckParity(const GFArray<GaloisFieldBinary> &word) const {
    if (!word.HasField() || word.Dimension() != 1 || word.Size() != num_cols_) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::CheckParity: word must be 1D with length equal "
          "to NumCols()");
    }
    return CheckParity(word.RawSpan());
  }

  GFMatrix<GaloisFieldBinary> ToDense(
      const std::shared_ptr<GaloisFieldBinary> &field) const {
    if (!field) {
      throw std::invalid_argument(
          "SparseBinaryMatrix::ToDense: field cannot be null");
    }
    using RawArray = GFMatrix<GaloisFieldBinary>::raw_array_type;
    using ShapeType = GFMatrix<GaloisFieldBinary>::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{num_rows_, num_cols_});
    raw.fill(uint8_t{0});
    uint8_t *data = raw.data();
    for (std::size_t i = 0; i < num_rows_; ++i) {
      uint8_t *row_data = data + i * num_cols_;
      const std::size_t end = row_ptr_[i + 1];
      for (std::size_t e = row_ptr_[i]; e < end; ++e) {
        row_data[row_col_idx_[e]] = 1U;
      }
    }
    return GFMatrix<GaloisFieldBinary>::FromRaw(field, std::move(raw));
  }

 private:
  void BuildFromEdges(
      std::size_t num_rows, std::size_t num_cols,
      const std::vector<std::pair<std::size_t, std::size_t>> &edges) {
    num_rows_ = num_rows;
    num_cols_ = num_cols;

    std::vector<std::pair<std::size_t, std::size_t>> sorted_edges = edges;
    for (const auto &[r, c] : sorted_edges) {
      if (r >= num_rows_ || c >= num_cols_) {
        throw std::out_of_range(
            "SparseBinaryMatrix: edge index out of matrix dimensions");
      }
    }
    std::sort(sorted_edges.begin(), sorted_edges.end());
    sorted_edges.erase(std::unique(sorted_edges.begin(), sorted_edges.end()),
                       sorted_edges.end());

    const std::size_t nnz = sorted_edges.size();

    // Build CSR.
    row_ptr_.assign(num_rows_ + 1, 0);
    row_col_idx_.resize(nnz);
    for (const auto &[r, c] : sorted_edges) {
      ++row_ptr_[r + 1];
    }
    for (std::size_t i = 0; i < num_rows_; ++i) {
      row_ptr_[i + 1] += row_ptr_[i];
    }
    for (std::size_t e = 0; e < nnz; ++e) {
      row_col_idx_[e] = sorted_edges[e].second;
    }

    // Build CSC and bidirectional edge-index cross-references.
    col_ptr_.assign(num_cols_ + 1, 0);
    for (const auto &[r, c] : sorted_edges) {
      ++col_ptr_[c + 1];
    }
    for (std::size_t j = 0; j < num_cols_; ++j) {
      col_ptr_[j + 1] += col_ptr_[j];
    }

    col_row_idx_.resize(nnz);
    col_to_csr_edge_.resize(nnz);
    csr_to_csc_edge_.resize(nnz);
    std::vector<std::size_t> col_cursor(col_ptr_.begin(),
                                        col_ptr_.begin() + num_cols_);

    // Iterating in CSR order (sorted by row then col) guarantees each column's
    // row indices are inserted in strictly ascending order.
    for (std::size_t i = 0; i < num_rows_; ++i) {
      const std::size_t end = row_ptr_[i + 1];
      for (std::size_t e_csr = row_ptr_[i]; e_csr < end; ++e_csr) {
        const std::size_t c = row_col_idx_[e_csr];
        const std::size_t e_csc = col_cursor[c]++;
        col_row_idx_[e_csc] = i;
        col_to_csr_edge_[e_csc] = e_csr;
        csr_to_csc_edge_[e_csr] = e_csc;
      }
    }
  }

  std::size_t num_rows_ = 0;
  std::size_t num_cols_ = 0;
  std::vector<std::size_t> row_ptr_{0};
  std::vector<std::size_t> row_col_idx_;
  std::vector<std::size_t> col_ptr_{0};
  std::vector<std::size_t> col_row_idx_;
  std::vector<std::size_t> col_to_csr_edge_;
  std::vector<std::size_t> csr_to_csc_edge_;
};

enum class BpAlgorithm : std::uint8_t {
  kSumProduct,
  kNormalizedMinSum,
  kOffsetMinSum,
};

class LDPCCode;
class BitFlippingDecoder;
class BeliefPropagationDecoder;

// Binary Low-Density Parity-Check (LDPC) linear code over GF(2).
class LDPCCode : public AbstractLinearCode<GaloisFieldBinary> {
 public:
  using Base = AbstractLinearCode<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using matrix_type = GFMatrix<GaloisFieldBinary>;
  using syndrome_type = GFArray<GaloisFieldBinary>;

  LDPCCode(std::shared_ptr<GaloisFieldBinary> field,
           SparseBinaryMatrix parity_check)
      : LDPCCode(std::move(field), std::move(parity_check),
                 PrivateConstructTag{}) {}

  LDPCCode(std::shared_ptr<GaloisFieldBinary> field,
           const matrix_type &dense_parity_check)
      : LDPCCode(std::move(field),
                 SparseBinaryMatrix::FromDense(
                     ValidateDenseH(field, dense_parity_check)),
                 PrivateConstructTag{}) {}

  explicit LDPCCode(SparseBinaryMatrix parity_check)
      : LDPCCode(std::make_shared<GaloisFieldBinary>(), std::move(parity_check),
                 PrivateConstructTag{}) {}

  explicit LDPCCode(const matrix_type &dense_parity_check)
      : LDPCCode(dense_parity_check.HasField()
                     ? dense_parity_check.Field()
                     : throw std::invalid_argument(
                           "LDPCCode: dense parity-check matrix has null "
                           "field"),
                 dense_parity_check) {}

  static LDPCCode FromQuasiCyclic(
      std::shared_ptr<GaloisFieldBinary> field,
      const std::vector<std::vector<int>> &base_matrix,
      std::size_t lifting_factor_z) {
    return LDPCCode(std::move(field), SparseBinaryMatrix::FromQuasiCyclic(
                                          base_matrix, lifting_factor_z));
  }

  // Gallager regular (n, col_weight, row_weight) LDPC construction.
  // Requires col_weight >= 1, row_weight >= 2, col_weight < row_weight, and
  // n divisible by row_weight. Uses multi-trial permutation selection to
  // minimize 4-cycles while preserving exact column and row weights.
  static LDPCCode MakeRegularLDPC(std::shared_ptr<GaloisFieldBinary> field,
                                  std::size_t n, std::size_t col_weight,
                                  std::size_t row_weight, uint64_t seed = 42) {
    if (!field) {
      throw std::invalid_argument(
          "LDPCCode::MakeRegularLDPC: field cannot be null");
    }
    if (n == 0) {
      throw std::invalid_argument(
          "LDPCCode::MakeRegularLDPC: block length n must be positive");
    }
    if (col_weight == 0 || row_weight < 2) {
      throw std::invalid_argument(
          "LDPCCode::MakeRegularLDPC: require col_weight >= 1 and "
          "row_weight >= 2");
    }
    if (col_weight >= row_weight) {
      throw std::invalid_argument(
          "LDPCCode::MakeRegularLDPC: require col_weight < row_weight");
    }
    if (n % row_weight != 0) {
      throw std::invalid_argument(
          "LDPCCode::MakeRegularLDPC: n must be divisible by row_weight");
    }

    const std::size_t sub_rows = n / row_weight;
    const std::size_t total_rows = col_weight * sub_rows;
    std::mt19937_64 rng(seed);

    // Track row groups already created: band 0 is deterministic block identity.
    std::vector<std::vector<std::size_t>> rows;
    rows.reserve(total_rows);
    for (std::size_t r = 0; r < sub_rows; ++r) {
      std::vector<std::size_t> row(row_weight);
      for (std::size_t k = 0; k < row_weight; ++k) {
        row[k] = r * row_weight + k;
      }
      rows.push_back(std::move(row));
    }

    // Adjacency map: for each column, which existing rows contain it.
    std::vector<std::vector<std::size_t>> col_to_rows(n);
    for (std::size_t r = 0; r < sub_rows; ++r) {
      for (std::size_t c : rows[r]) {
        col_to_rows[c].push_back(r);
      }
    }

    auto count_new_four_cycles = [&](const std::vector<std::size_t> &perm) {
      std::size_t cycles = 0;
      std::vector<std::size_t> overlap(rows.size(), 0);
      std::vector<std::size_t> touched;
      touched.reserve(row_weight * col_weight);
      for (std::size_t r = 0; r < sub_rows; ++r) {
        touched.clear();
        for (std::size_t k = 0; k < row_weight; ++k) {
          const std::size_t c = perm[r * row_weight + k];
          for (std::size_t prev_r : col_to_rows[c]) {
            if (overlap[prev_r] == 0) {
              touched.push_back(prev_r);
            }
            ++overlap[prev_r];
          }
        }
        for (std::size_t prev_r : touched) {
          const std::size_t ov = overlap[prev_r];
          if (ov >= 2) {
            cycles += (ov * (ov - 1)) / 2;
          }
          overlap[prev_r] = 0;
        }
      }
      return cycles;
    };

    constexpr std::size_t kCandidateTrials = 16;
    std::vector<std::size_t> base_perm(n);
    std::iota(base_perm.begin(), base_perm.end(), std::size_t{0});

    for (std::size_t band = 1; band < col_weight; ++band) {
      std::vector<std::size_t> best_perm;
      std::size_t best_cycles = std::numeric_limits<std::size_t>::max();
      for (std::size_t trial = 0; trial < kCandidateTrials; ++trial) {
        std::vector<std::size_t> cand = base_perm;
        std::shuffle(cand.begin(), cand.end(), rng);
        const std::size_t cycles = count_new_four_cycles(cand);
        if (cycles < best_cycles) {
          best_cycles = cycles;
          best_perm = std::move(cand);
          if (best_cycles == 0) {
            break;
          }
        }
      }

      for (std::size_t r = 0; r < sub_rows; ++r) {
        const std::size_t row_idx = rows.size();
        std::vector<std::size_t> row(row_weight);
        for (std::size_t k = 0; k < row_weight; ++k) {
          const std::size_t c = best_perm[r * row_weight + k];
          row[k] = c;
          col_to_rows[c].push_back(row_idx);
        }
        rows.push_back(std::move(row));
      }
    }

    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(n * col_weight);
    for (std::size_t r = 0; r < total_rows; ++r) {
      for (std::size_t c : rows[r]) {
        edges.emplace_back(r, c);
      }
    }
    return LDPCCode(std::move(field), SparseBinaryMatrix(total_rows, n, edges));
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  const SparseBinaryMatrix &SparseParityCheckMatrix() const noexcept {
    return sparse_h_;
  }

  syndrome_type SyndromeSparse(const codeword_type &word) const {
    if (!word.HasField() || word.Field() != this->field_) {
      throw std::invalid_argument("LDPCCode::SyndromeSparse: field mismatch");
    }
    return sparse_h_.Syndrome(word);
  }

  bool CheckParity(const codeword_type &word) const {
    if (!word.HasField() || word.Field() != this->field_) {
      throw std::invalid_argument("LDPCCode::CheckParity: field mismatch");
    }
    return sparse_h_.CheckParity(word);
  }

  bool CheckParity(std::span<const uint8_t> bits) const {
    return sparse_h_.CheckParity(bits);
  }

  bool Contains(const codeword_type &word) const override {
    if (!word.HasField() || word.Field() != this->field_ ||
        word.Dimension() != 1 || word.Size() != this->length_) {
      return false;
    }
    return sparse_h_.CheckParity(word.RawSpan());
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) +
           "] binary LDPC code (m=" + std::to_string(sparse_h_.NumRows()) +
           ", nnz=" + std::to_string(sparse_h_.NumNonZeros()) + ")";
  }

 private:
  struct PrivateConstructTag {};

  static const matrix_type &ValidateDenseH(
      const std::shared_ptr<GaloisFieldBinary> &field,
      const matrix_type &dense) {
    if (!field) {
      throw std::invalid_argument("LDPCCode: field cannot be null");
    }
    if (!dense.HasField() || dense.Field() != field) {
      throw std::invalid_argument("LDPCCode: parity-check field mismatch");
    }
    if (dense.Dimension() != 2 || dense.Shape(1) == 0) {
      throw std::invalid_argument(
          "LDPCCode: parity-check must be a 2D matrix with positive columns");
    }
    return dense;
  }

  static std::size_t RequirePositiveCols(const SparseBinaryMatrix &h) {
    if (h.NumCols() == 0) {
      throw std::invalid_argument(
          "LDPCCode: parity-check matrix must have positive column count");
    }
    return h.NumCols();
  }

  LDPCCode(std::shared_ptr<GaloisFieldBinary> field,
           SparseBinaryMatrix parity_check, PrivateConstructTag)
      : Base(field, RequirePositiveCols(parity_check), 0, Metric::HAMMING,
             "GeneratorMatrix", "BitFlipping"),
        sparse_h_(std::move(parity_check)),
        parity_check_matrix_(sparse_h_.ToDense(this->field_)),
        generator_matrix_(
            Base::ComputeGeneratorFromParityCheck(parity_check_matrix_)) {
    this->dimension_ = generator_matrix_.Shape(0);
    RegisterEncodersAndDecoders();
  }

  inline void RegisterEncodersAndDecoders();

  SparseBinaryMatrix sparse_h_;
  matrix_type parity_check_matrix_;
  matrix_type generator_matrix_;
};

// Iterative Gallager hard-decision bit-flipping decoder on the LDPC Tanner
// graph. Supports both regular and irregular parity-check matrices.
class BitFlippingDecoder : public Decoder<GaloisFieldBinary> {
 public:
  using Base = Decoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using result_type = DecodingResult<GaloisFieldBinary>;

  explicit BitFlippingDecoder(const LDPCCode &code,
                              std::size_t max_iterations = 50)
      : Base(code), ldpc_code_(&code), max_iterations_(max_iterations) {}

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "BitFlippingDecoder::DecodeWithResult: field mismatch");
    }
    const std::size_t n = ldpc_code_->Length();
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "BitFlippingDecoder::DecodeWithResult: received word length "
          "mismatch");
    }

    const auto &h = ldpc_code_->SparseParityCheckMatrix();
    const std::size_t m = h.NumRows();
    const auto &row_ptr = h.RowPtr();
    const auto &row_col_idx = h.RowColIdx();
    const auto &col_ptr = h.ColPtr();
    const auto &col_row_idx = h.ColRowIdx();

    const uint8_t *rx_raw = received_word.RawData().data();
    std::vector<uint8_t> est(n);
    for (std::size_t j = 0; j < n; ++j) {
      est[j] = static_cast<uint8_t>(rx_raw[j] & 1U);
    }

    std::vector<uint8_t> syn(m, 0);
    auto compute_syndrome = [&](const std::vector<uint8_t> &bits,
                                std::vector<uint8_t> &out_syn) {
      std::size_t unsat = 0;
      for (std::size_t i = 0; i < m; ++i) {
        uint8_t parity = 0;
        const std::size_t end = row_ptr[i + 1];
        for (std::size_t e = row_ptr[i]; e < end; ++e) {
          parity = static_cast<uint8_t>(parity ^ bits[row_col_idx[e]]);
        }
        out_syn[i] = parity;
        unsat += parity;
      }
      return unsat;
    };

    std::size_t unsat_total = compute_syndrome(est, syn);
    std::vector<int> scores(n, 0);
    std::vector<std::size_t> unsat_counts(n, 0);
    std::vector<std::size_t> flip_candidates;
    flip_candidates.reserve(n);
    std::vector<uint8_t> trial_syn(m, 0);
    std::vector<std::size_t> last_flipped(n, 0);

    // Local O(dv * dc) lookahead: after flipping `col`, compute the resulting
    // unsatisfied check count and the maximum violation score among variables
    // sharing checks with `col`. This cleanly distinguishes true error bits
    // from 6-cycle decoy bits when multiple variables tie for max_score.
    auto evaluate_single_flip = [&](std::size_t col) {
      const std::size_t c_start = col_ptr[col];
      const std::size_t c_end = col_ptr[col + 1];
      std::size_t delta_unsat = unsat_total;
      for (std::size_t e = c_start; e < c_end; ++e) {
        const std::size_t chk = col_row_idx[e];
        if (syn[chk] != 0) {
          --delta_unsat;
          syn[chk] = 0;
        } else {
          ++delta_unsat;
          syn[chk] = 1;
        }
      }

      int best_next_score = std::numeric_limits<int>::min();
      if (delta_unsat == 0) {
        best_next_score = std::numeric_limits<int>::max();
      } else {
        for (std::size_t e = c_start; e < c_end; ++e) {
          const std::size_t chk = col_row_idx[e];
          const std::size_t r_end = row_ptr[chk + 1];
          for (std::size_t re = row_ptr[chk]; re < r_end; ++re) {
            const std::size_t nb = row_col_idx[re];
            if (nb == col) {
              continue;
            }
            const std::size_t nb_start = col_ptr[nb];
            const std::size_t nb_end = col_ptr[nb + 1];
            const std::size_t nb_deg = nb_end - nb_start;
            std::size_t nb_unsat = 0;
            for (std::size_t ne = nb_start; ne < nb_end; ++ne) {
              nb_unsat += syn[col_row_idx[ne]];
            }
            if (nb_unsat > 0) {
              const int nb_score =
                  2 * static_cast<int>(nb_unsat) - static_cast<int>(nb_deg);
              if (nb_score > best_next_score) {
                best_next_score = nb_score;
              }
            }
          }
        }
      }

      // Restore syn in place.
      for (std::size_t e = c_start; e < c_end; ++e) {
        syn[col_row_idx[e]] ^= 1U;
      }
      return std::pair<std::size_t, int>{delta_unsat, best_next_score};
    };

    for (std::size_t iter = 0; iter <= max_iterations_; ++iter) {
      if (unsat_total == 0) {
        codeword_type cw = BuildCodeword(field, est);
        const std::size_t num_errors = HammingDistance(est, rx_raw);
        message_type msg = ldpc_code_->Unencode(cw);
        return result_type::Success(std::move(cw), std::move(msg), num_errors);
      }
      if (iter == max_iterations_) {
        break;
      }

      // Compute net violation score 2 * unsat(j) - deg(j) = unsat(j) - sat(j)
      // for each variable node j.
      int max_score = std::numeric_limits<int>::min();
      int max_non_tabu_score = std::numeric_limits<int>::min();
      std::size_t best_single_col = 0;
      std::size_t best_single_unsat = 0;
      bool found_nonzero = false;
      bool found_non_tabu = false;

      for (std::size_t j = 0; j < n; ++j) {
        const std::size_t c_start = col_ptr[j];
        const std::size_t c_end = col_ptr[j + 1];
        const std::size_t deg = c_end - c_start;
        std::size_t unsat_j = 0;
        for (std::size_t e = c_start; e < c_end; ++e) {
          unsat_j += syn[col_row_idx[e]];
        }
        unsat_counts[j] = unsat_j;
        const int score = 2 * static_cast<int>(unsat_j) - static_cast<int>(deg);
        scores[j] = score;
        if (unsat_j > 0) {
          const bool is_tabu =
              (last_flipped[j] > 0) && (iter + 1 - last_flipped[j] <= 2);
          if (!is_tabu && (!found_non_tabu || score > max_non_tabu_score)) {
            found_non_tabu = true;
            max_non_tabu_score = score;
          }
          if (!found_nonzero || score > max_score ||
              (score == max_score && unsat_j > best_single_unsat)) {
            found_nonzero = true;
            max_score = score;
            best_single_unsat = unsat_j;
            best_single_col = j;
          }
        }
      }

      if (!found_nonzero) {
        break;
      }

      const bool use_tabu = (max_score <= 0) && found_non_tabu;
      const int target_score = use_tabu ? max_non_tabu_score : max_score;
      flip_candidates.clear();
      for (std::size_t j = 0; j < n; ++j) {
        if (unsat_counts[j] == 0) {
          continue;
        }
        const bool is_tabu =
            (last_flipped[j] > 0) && (iter + 1 - last_flipped[j] <= 2);
        if (use_tabu && is_tabu) {
          continue;
        }
        if (scores[j] == target_score) {
          flip_candidates.push_back(j);
        }
      }
      if (flip_candidates.empty()) {
        flip_candidates.push_back(best_single_col);
      }

      // Fast path: if target_score > 0, try flipping all tied max-score bits
      // in parallel. If that stalls (e.g., 6-cycle decoy tie) or target_score
      // <= 0, pick the best single bit via local neighbor-score lookahead.
      bool applied_parallel = false;
      if (target_score > 0) {
        for (std::size_t col : flip_candidates) {
          est[col] ^= 1U;
        }
        const std::size_t next_unsat = compute_syndrome(est, trial_syn);
        if (flip_candidates.size() == 1 || next_unsat < unsat_total) {
          for (std::size_t col : flip_candidates) {
            last_flipped[col] = iter + 1;
          }
          unsat_total = next_unsat;
          syn.swap(trial_syn);
          applied_parallel = true;
        } else {
          for (std::size_t col : flip_candidates) {
            est[col] ^= 1U;
          }
        }
      }

      if (!applied_parallel) {
        std::size_t chosen_col = flip_candidates[0];
        auto [best_delta, best_next_score] = evaluate_single_flip(chosen_col);
        for (std::size_t idx = 1; idx < flip_candidates.size(); ++idx) {
          const std::size_t col = flip_candidates[idx];
          const auto [delta, next_score] = evaluate_single_flip(col);
          if (delta == 0 ||
              (best_delta != 0 &&
               (next_score > best_next_score ||
                (next_score == best_next_score && delta < best_delta)))) {
            chosen_col = col;
            best_delta = delta;
            best_next_score = next_score;
            if (best_delta == 0) {
              break;
            }
          }
        }
        est[chosen_col] ^= 1U;
        last_flipped[chosen_col] = iter + 1;
        unsat_total = compute_syndrome(est, syn);
      }
    }

    result_type res = result_type::Failure(
        "BitFlippingDecoder: failed to satisfy all parity checks within " +
        std::to_string(max_iterations_) + " iterations");
    res.codeword = BuildCodeword(field, est);
    res.num_errors = HammingDistance(est, rx_raw);
    return res;
  }

  std::string Name() const override { return "BitFlipping"; }

  std::size_t DecodingRadius() const override { return 0; }

  std::size_t MaxIterations() const noexcept { return max_iterations_; }

 private:
  static codeword_type BuildCodeword(
      const std::shared_ptr<GaloisFieldBinary> &field,
      const std::vector<uint8_t> &bits) {
    using RawArray = codeword_type::raw_array_type;
    using ShapeType = codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{bits.size()});
    std::copy(bits.begin(), bits.end(), raw.data());
    return codeword_type::FromRaw(field, std::move(raw));
  }

  static std::size_t HammingDistance(const std::vector<uint8_t> &est,
                                     const uint8_t *rx_raw) noexcept {
    std::size_t dist = 0;
    for (std::size_t j = 0; j < est.size(); ++j) {
      if ((est[j] & 1U) != (rx_raw[j] & 1U)) {
        ++dist;
      }
    }
    return dist;
  }

  const LDPCCode *ldpc_code_;
  std::size_t max_iterations_;
};

// Soft-decision Belief Propagation (Sum-Product, Normalized Min-Sum, Offset
// Min-Sum) LDPC decoder. Uses AWGNChannel's LLR convention:
//   positive LLR favors bit 0, negative LLR favors bit 1.
class BeliefPropagationDecoder : public Decoder<GaloisFieldBinary> {
 public:
  using Base = Decoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using result_type = DecodingResult<GaloisFieldBinary>;

  static constexpr double kMaxLlr = 50.0;

  explicit BeliefPropagationDecoder(
      const LDPCCode &code, BpAlgorithm algorithm = BpAlgorithm::kSumProduct,
      std::size_t max_iterations = 50, double alpha = 0.75, double beta = 0.5,
      double hard_llr_magnitude = 4.0)
      : Base(code),
        ldpc_code_(&code),
        algorithm_(algorithm),
        max_iterations_(max_iterations),
        alpha_(alpha),
        beta_(beta),
        hard_llr_magnitude_(hard_llr_magnitude) {
    if (algorithm_ == BpAlgorithm::kNormalizedMinSum &&
        (!std::isfinite(alpha_) || alpha_ <= 0.0 || alpha_ > 1.0)) {
      throw std::invalid_argument(
          "BeliefPropagationDecoder: normalization factor alpha must be in "
          "(0, 1]");
    }
    if (algorithm_ == BpAlgorithm::kOffsetMinSum &&
        (!std::isfinite(beta_) || beta_ < 0.0)) {
      throw std::invalid_argument(
          "BeliefPropagationDecoder: offset beta must be finite and >= 0");
    }
    if (!std::isfinite(hard_llr_magnitude_) || hard_llr_magnitude_ <= 0.0) {
      throw std::invalid_argument(
          "BeliefPropagationDecoder: hard_llr_magnitude must be positive and "
          "finite");
    }
  }

  result_type DecodeLlr(const std::vector<double> &channel_llrs) const {
    return DecodeLlr(channel_llrs, max_iterations_);
  }

  result_type DecodeLlr(const std::vector<double> &channel_llrs,
                        std::size_t max_iterations) const {
    const std::size_t n = ldpc_code_->Length();
    detail::ValidateLlrVector(channel_llrs, n,
                              "BeliefPropagationDecoder::DecodeLlr");

    auto field = this->code_->Field();
    const auto &h = ldpc_code_->SparseParityCheckMatrix();
    const std::size_t m = h.NumRows();
    const std::size_t nnz = h.NumNonZeros();
    const auto &row_ptr = h.RowPtr();
    const auto &row_col_idx = h.RowColIdx();
    const auto &col_ptr = h.ColPtr();
    const auto &col_to_csr = h.ColToCsrEdge();

    std::vector<double> ch_llr(n);
    std::vector<uint8_t> input_hard(n);
    for (std::size_t j = 0; j < n; ++j) {
      ch_llr[j] = std::clamp(channel_llrs[j], -kMaxLlr, kMaxLlr);
      input_hard[j] = (channel_llrs[j] >= 0.0) ? uint8_t{0} : uint8_t{1};
    }

    // Early-terminate at iteration 0 if channel hard decision already satisfies
    // all parity checks.
    if (h.CheckParity(input_hard)) {
      codeword_type cw = BuildCodeword(field, input_hard);
      message_type msg = ldpc_code_->Unencode(cw);
      return result_type::Success(std::move(cw), std::move(msg), 0);
    }

    // Initialize variable-to-check messages v2c[e] = ch_llr[col(e)].
    std::vector<double> v2c(nnz, 0.0);
    std::vector<double> c2v(nnz, 0.0);
    for (std::size_t e = 0; e < nnz; ++e) {
      v2c[e] = ch_llr[row_col_idx[e]];
    }

    std::vector<uint8_t> hard = input_hard;
    std::vector<double> prefix;
    std::vector<double> suffix;

    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
      // 1. Check-to-variable (horizontal) step.
      for (std::size_t i = 0; i < m; ++i) {
        const std::size_t start = row_ptr[i];
        const std::size_t end = row_ptr[i + 1];
        const std::size_t deg = end - start;
        if (deg == 0) {
          continue;
        }
        if (deg == 1) {
          // Single-variable parity check forces x_j = 0.
          c2v[start] = kMaxLlr;
          continue;
        }

        if (algorithm_ == BpAlgorithm::kSumProduct) {
          prefix.resize(deg);
          suffix.resize(deg);
          prefix[0] = v2c[start];
          for (std::size_t k = 1; k + 1 < deg; ++k) {
            prefix[k] = detail::BoxPlusLlr(prefix[k - 1], v2c[start + k]);
          }
          suffix[deg - 1] = v2c[start + deg - 1];
          for (std::size_t k = deg - 1; k > 1; --k) {
            suffix[k - 1] = detail::BoxPlusLlr(v2c[start + k - 1], suffix[k]);
          }
          c2v[start] = std::clamp(suffix[1], -kMaxLlr, kMaxLlr);
          c2v[start + deg - 1] = std::clamp(prefix[deg - 2], -kMaxLlr, kMaxLlr);
          for (std::size_t k = 1; k + 1 < deg; ++k) {
            const double val = detail::BoxPlusLlr(prefix[k - 1], suffix[k + 1]);
            c2v[start + k] = std::clamp(val, -kMaxLlr, kMaxLlr);
          }
        } else {
          // Normalized Min-Sum or Offset Min-Sum in O(deg).
          int total_sign = 1;
          double min1 = std::numeric_limits<double>::infinity();
          double min2 = std::numeric_limits<double>::infinity();
          std::size_t min1_edge = start;
          for (std::size_t e = start; e < end; ++e) {
            const double val = v2c[e];
            if (val < 0.0) {
              total_sign = -total_sign;
            }
            const double mag = std::abs(val);
            if (mag < min1) {
              min2 = min1;
              min1 = mag;
              min1_edge = e;
            } else if (mag < min2) {
              min2 = mag;
            }
          }

          for (std::size_t e = start; e < end; ++e) {
            const int excl_sign = (v2c[e] < 0.0) ? -total_sign : total_sign;
            const double excl_mag = (e == min1_edge) ? min2 : min1;
            double scaled_mag = excl_mag;
            if (algorithm_ == BpAlgorithm::kNormalizedMinSum) {
              scaled_mag = alpha_ * excl_mag;
            } else {
              scaled_mag = std::max(0.0, excl_mag - beta_);
            }
            c2v[e] = std::clamp(static_cast<double>(excl_sign) * scaled_mag,
                                -kMaxLlr, kMaxLlr);
          }
        }
      }

      // 2. Variable-to-check (vertical) step + posterior hard decisions.
      for (std::size_t j = 0; j < n; ++j) {
        const std::size_t c_start = col_ptr[j];
        const std::size_t c_end = col_ptr[j + 1];
        double total = ch_llr[j];
        for (std::size_t e_csc = c_start; e_csc < c_end; ++e_csc) {
          total += c2v[col_to_csr[e_csc]];
        }
        hard[j] = (total >= 0.0) ? uint8_t{0} : uint8_t{1};
        for (std::size_t e_csc = c_start; e_csc < c_end; ++e_csc) {
          const std::size_t e_csr = col_to_csr[e_csc];
          v2c[e_csr] = std::clamp(total - c2v[e_csr], -kMaxLlr, kMaxLlr);
        }
      }

      // 3. Syndrome early termination.
      if (h.CheckParity(hard)) {
        codeword_type cw = BuildCodeword(field, hard);
        const std::size_t num_errors = HammingDistance(hard, input_hard);
        message_type msg = ldpc_code_->Unencode(cw);
        return result_type::Success(std::move(cw), std::move(msg), num_errors);
      }
    }

    result_type res = result_type::Failure(
        "BeliefPropagationDecoder: parity check failed after " +
        std::to_string(max_iterations) + " iterations");
    res.codeword = BuildCodeword(field, hard);
    res.num_errors = HammingDistance(hard, input_hard);
    return res;
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "BeliefPropagationDecoder::DecodeWithResult: field mismatch");
    }
    const std::size_t n = ldpc_code_->Length();
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "BeliefPropagationDecoder::DecodeWithResult: received word length "
          "mismatch");
    }
    const uint8_t *rx_raw = received_word.RawData().data();
    std::vector<double> llrs(n);
    for (std::size_t j = 0; j < n; ++j) {
      llrs[j] =
          ((rx_raw[j] & 1U) == 0) ? hard_llr_magnitude_ : -hard_llr_magnitude_;
    }
    return DecodeLlr(llrs, max_iterations_);
  }

  std::string Name() const override {
    switch (algorithm_) {
      case BpAlgorithm::kSumProduct:
        return "BeliefPropagation(SumProduct)";
      case BpAlgorithm::kNormalizedMinSum:
        return "BeliefPropagation(NormalizedMinSum)";
      case BpAlgorithm::kOffsetMinSum:
        return "BeliefPropagation(OffsetMinSum)";
    }
    return "BeliefPropagation";
  }

  std::size_t DecodingRadius() const override { return 0; }

  BpAlgorithm Algorithm() const noexcept { return algorithm_; }
  std::size_t MaxIterations() const noexcept { return max_iterations_; }
  double Alpha() const noexcept { return alpha_; }
  double Beta() const noexcept { return beta_; }

 private:
  static codeword_type BuildCodeword(
      const std::shared_ptr<GaloisFieldBinary> &field,
      const std::vector<uint8_t> &bits) {
    using RawArray = codeword_type::raw_array_type;
    using ShapeType = codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{bits.size()});
    std::copy(bits.begin(), bits.end(), raw.data());
    return codeword_type::FromRaw(field, std::move(raw));
  }

  static std::size_t HammingDistance(const std::vector<uint8_t> &a,
                                     const std::vector<uint8_t> &b) noexcept {
    std::size_t dist = 0;
    for (std::size_t j = 0; j < a.size(); ++j) {
      if ((a[j] & 1U) != (b[j] & 1U)) {
        ++dist;
      }
    }
    return dist;
  }

  const LDPCCode *ldpc_code_;
  BpAlgorithm algorithm_;
  std::size_t max_iterations_;
  double alpha_;
  double beta_;
  double hard_llr_magnitude_;
};

inline void LDPCCode::RegisterEncodersAndDecoders() {
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<GaloisFieldBinary> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<GaloisFieldBinary>>(lin);
      });
  this->RegisterDecoder(
      "BitFlipping", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &ldpc = dynamic_cast<const LDPCCode &>(code);
        return std::make_shared<BitFlippingDecoder>(ldpc);
      });
  this->RegisterDecoder(
      "BeliefPropagation", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &ldpc = dynamic_cast<const LDPCCode &>(code);
        return std::make_shared<BeliefPropagationDecoder>(ldpc);
      });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_LDPC_HPP_
