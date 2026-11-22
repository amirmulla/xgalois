#ifndef XGALOIS_CODING_PRODUCT_HPP_
#define XGALOIS_CODING_PRODUCT_HPP_

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
struct ProductMatrixDecodingResult {
  using matrix_type = GFMatrix<GaloisField>;

  matrix_type codeword_matrix;
  matrix_type message_matrix;
  std::size_t num_errors = 0;
  std::size_t iterations = 0;
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }

  static ProductMatrixDecodingResult Success(matrix_type cw_mat,
                                             matrix_type msg_mat,
                                             std::size_t errors = 0,
                                             std::size_t iters = 0) {
    ProductMatrixDecodingResult res;
    res.codeword_matrix = std::move(cw_mat);
    res.message_matrix = std::move(msg_mat);
    res.num_errors = errors;
    res.iterations = iters;
    res.success = true;
    return res;
  }

  static ProductMatrixDecodingResult Failure(std::string reason,
                                             std::size_t iters = 0) {
    ProductMatrixDecodingResult res;
    res.iterations = iters;
    res.success = false;
    res.error_message = std::move(reason);
    return res;
  }
};

struct ChasePyndiahDecodingResult {
  GFMatrix<GaloisFieldBinary> codeword_matrix;
  GFMatrix<GaloisFieldBinary> message_matrix;
  std::vector<double> posterior_llrs;
  std::size_t iterations = 0;
  bool valid_codeword = false;

  bool Ok() const noexcept { return valid_codeword; }
  explicit operator bool() const noexcept { return valid_codeword; }
};

template <typename GaloisField>
class ProductCode;

template <typename GaloisField>
class ProductEncoder : public Encoder<GaloisField> {
 public:
  using Base = Encoder<GaloisField>;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;

  explicit ProductEncoder(const ProductCode<GaloisField> &code)
      : Base(code), product_code_(&code) {}

  codeword_type Encode(const message_type &message) const override {
    if (!message.HasField() || message.Field() != this->code_->Field()) {
      throw std::invalid_argument(
          "ProductEncoder::Encode: message field mismatch");
    }
    if (message.Dimension() != 1 ||
        message.Size() != product_code_->Dimension()) {
      throw std::invalid_argument(
          "ProductEncoder::Encode: message must be 1D with length k1 * k2");
    }
    auto msg_mat = product_code_->UnflattenMatrix(
        message, product_code_->ColDimension(), product_code_->RowDimension());
    auto cw_mat = product_code_->EncodeMatrix(msg_mat);
    return product_code_->FlattenMatrix(cw_mat);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    if (!codeword.HasField() || codeword.Field() != this->code_->Field()) {
      throw std::invalid_argument(
          "ProductEncoder::Unencode: codeword field mismatch");
    }
    if (codeword.Dimension() != 1 ||
        codeword.Size() != product_code_->Length()) {
      throw std::invalid_argument(
          "ProductEncoder::Unencode: codeword must be 1D with length n1 * n2");
    }
    auto cw_mat = product_code_->UnflattenMatrix(
        codeword, product_code_->ColLength(), product_code_->RowLength());
    auto msg_mat = product_code_->UnencodeMatrix(cw_mat);
    return product_code_->FlattenMatrix(msg_mat);
  }

  std::string Name() const override { return "Product"; }

 private:
  const ProductCode<GaloisField> *product_code_;
};

template <typename GaloisField>
class ProductHardDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using codeword_type = typename Base::codeword_type;
  using result_type = typename Base::result_type;

  explicit ProductHardDecoder(const ProductCode<GaloisField> &code,
                              std::size_t max_iterations = 6)
      : Base(code), product_code_(&code), max_iterations_(max_iterations) {}

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    if (!received_word.HasField() ||
        received_word.Field() != this->code_->Field()) {
      throw std::invalid_argument(
          "ProductHardDecoder::DecodeWithResult: field mismatch");
    }
    if (received_word.Dimension() != 1 ||
        received_word.Size() != product_code_->Length()) {
      throw std::invalid_argument(
          "ProductHardDecoder::DecodeWithResult: received word must be 1D "
          "with length n1 * n2");
    }
    auto recv_mat = product_code_->UnflattenMatrix(
        received_word, product_code_->ColLength(), product_code_->RowLength());
    auto mat_res = product_code_->DecodeMatrixHard(recv_mat, max_iterations_);
    if (!mat_res.Ok()) {
      return result_type::Failure(std::move(mat_res.error_message));
    }
    return result_type::Success(
        product_code_->FlattenMatrix(mat_res.codeword_matrix),
        product_code_->FlattenMatrix(mat_res.message_matrix),
        mat_res.num_errors);
  }

  std::string Name() const override { return "ProductHard"; }

  std::size_t DecodingRadius() const override {
    const std::size_t d = product_code_->MinimumDistance();
    return (d > 0) ? (d - 1) / 2 : 0;
  }

 private:
  const ProductCode<GaloisField> *product_code_;
  std::size_t max_iterations_;
};

template <typename GaloisField>
class ProductCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = typename Base::matrix_type;
  using result_type = typename Base::result_type;

  ProductCode(std::shared_ptr<const AbstractLinearCode<GaloisField>> row_code,
              std::shared_ptr<const AbstractLinearCode<GaloisField>> col_code)
      : Base(ValidateAndGetField(row_code, col_code),
             row_code->Length() * col_code->Length(),
             row_code->Dimension() * col_code->Dimension(), Metric::HAMMING,
             "Product", "ProductHard"),
        row_code_(std::move(row_code)),
        col_code_(std::move(col_code)),
        n1_(row_code_->Length()),
        k1_(row_code_->Dimension()),
        n2_(col_code_->Length()),
        k2_(col_code_->Dimension()) {
    RegisterEncodersAndDecoders();
  }

  ProductCode(const ProductCode &other)
      : Base(other),
        row_code_(other.row_code_),
        col_code_(other.col_code_),
        n1_(other.n1_),
        k1_(other.k1_),
        n2_(other.n2_),
        k2_(other.k2_) {
    std::lock_guard<std::mutex> lock(other.matrix_mutex_);
    cached_generator_ = other.cached_generator_;
    cached_parity_check_ = other.cached_parity_check_;
  }

  ProductCode &operator=(const ProductCode &other) {
    if (this != &other) {
      Base::operator=(other);
      row_code_ = other.row_code_;
      col_code_ = other.col_code_;
      n1_ = other.n1_;
      k1_ = other.k1_;
      n2_ = other.n2_;
      k2_ = other.k2_;
      std::scoped_lock lock(matrix_mutex_, other.matrix_mutex_);
      cached_generator_ = other.cached_generator_;
      cached_parity_check_ = other.cached_parity_check_;
    }
    return *this;
  }

  ProductCode(ProductCode &&other) noexcept
      : Base(std::move(other)),
        row_code_(std::move(other.row_code_)),
        col_code_(std::move(other.col_code_)),
        n1_(other.n1_),
        k1_(other.k1_),
        n2_(other.n2_),
        k2_(other.k2_),
        cached_generator_(std::move(other.cached_generator_)),
        cached_parity_check_(std::move(other.cached_parity_check_)) {}

  ProductCode &operator=(ProductCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      row_code_ = std::move(other.row_code_);
      col_code_ = std::move(other.col_code_);
      n1_ = other.n1_;
      k1_ = other.k1_;
      n2_ = other.n2_;
      k2_ = other.k2_;
      std::lock_guard<std::mutex> lock(matrix_mutex_);
      cached_generator_ = std::move(other.cached_generator_);
      cached_parity_check_ = std::move(other.cached_parity_check_);
    }
    return *this;
  }

  ~ProductCode() override = default;

  std::shared_ptr<const AbstractLinearCode<GaloisField>> RowCode() const {
    return row_code_;
  }
  std::shared_ptr<const AbstractLinearCode<GaloisField>> ColCode() const {
    return col_code_;
  }
  std::size_t RowLength() const noexcept { return n1_; }
  std::size_t RowDimension() const noexcept { return k1_; }
  std::size_t ColLength() const noexcept { return n2_; }
  std::size_t ColDimension() const noexcept { return k2_; }

  std::size_t MinimumDistance() const override {
    if (this->cached_min_distance_.has_value()) {
      return *this->cached_min_distance_;
    }
    const std::size_t d =
        row_code_->MinimumDistance() * col_code_->MinimumDistance();
    this->cached_min_distance_ = d;
    return d;
  }

  // Lazy Kronecker generator matrix: G = G_col (x) G_row, matching row-major
  // flattening of C = G_col^T * M * G_row.
  const matrix_type &GeneratorMatrix() const override {
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!cached_generator_.has_value()) {
      cached_generator_ = linalg::Kron(col_code_->GeneratorMatrix(),
                                       row_code_->GeneratorMatrix());
    }
    return *cached_generator_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    const matrix_type &G = GeneratorMatrix();
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!cached_parity_check_.has_value()) {
      cached_parity_check_ = Base::ComputeParityCheckFromGenerator(G);
    }
    return *cached_parity_check_;
  }

  // Fast 2D membership test: every row in row_code and every column in
  // col_code.
  bool ContainsMatrix(const matrix_type &mat) const {
    if (!mat.HasField() || mat.Field() != this->field_ ||
        mat.Dimension() != 2 || mat.Shape(0) != n2_ || mat.Shape(1) != n1_) {
      return false;
    }
    for (std::size_t r = 0; r < n2_; ++r) {
      if (!row_code_->Contains(ExtractRow(mat, r))) {
        return false;
      }
    }
    for (std::size_t c = 0; c < n1_; ++c) {
      if (!col_code_->Contains(ExtractCol(mat, c))) {
        return false;
      }
    }
    return true;
  }

  bool Contains(const codeword_type &word) const override {
    if (!word.HasField() || word.Field() != this->field_ ||
        word.Dimension() != 1 || word.Size() != this->length_) {
      return false;
    }
    return ContainsMatrix(UnflattenMatrix(word, n2_, n1_));
  }

  // Fast 2D encoding: O(k2 * n1 * k1 + n1 * n2 * k2) via row then column
  // component encoders.
  matrix_type EncodeMatrix(const matrix_type &msg_k2_by_k1) const {
    if (!msg_k2_by_k1.HasField() || msg_k2_by_k1.Field() != this->field_) {
      throw std::invalid_argument(
          "ProductCode::EncodeMatrix: message field mismatch");
    }
    if (msg_k2_by_k1.Dimension() != 2 || msg_k2_by_k1.Shape(0) != k2_ ||
        msg_k2_by_k1.Shape(1) != k1_) {
      throw std::invalid_argument(
          "ProductCode::EncodeMatrix: message matrix must have shape "
          "{k2, k1}");
    }
    if (k1_ == 0 || k2_ == 0) {
      return matrix_type::Zeros(this->field_, {n2_, n1_});
    }

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    // Step 1: encode each of the k2 rows using row_code -> {k2, n1}.
    RawArray inter_raw = RawArray::from_shape(ShapeType{k2_, n1_});
    storage_type *inter_data = inter_raw.data();
    for (std::size_t r = 0; r < k2_; ++r) {
      const auto row_msg = ExtractRow(msg_k2_by_k1, r);
      const auto row_cw = row_code_->Encode(row_msg);
      const storage_type *src = row_cw.RawData().data();
      std::copy(src, src + n1_, inter_data + r * n1_);
    }
    auto inter = matrix_type::FromRaw(this->field_, std::move(inter_raw));

    // Step 2: encode each of the n1 columns using col_code -> {n2, n1}.
    RawArray out_raw = RawArray::from_shape(ShapeType{n2_, n1_});
    storage_type *out_data = out_raw.data();
    for (std::size_t c = 0; c < n1_; ++c) {
      const auto col_msg = ExtractCol(inter, c);
      const auto col_cw = col_code_->Encode(col_msg);
      const storage_type *src = col_cw.RawData().data();
      for (std::size_t r = 0; r < n2_; ++r) {
        out_data[r * n1_ + c] = src[r];
      }
    }
    return matrix_type::FromRaw(this->field_, std::move(out_raw));
  }

  matrix_type UnencodeMatrix(const matrix_type &cw_n2_by_n1) const {
    if (!cw_n2_by_n1.HasField() || cw_n2_by_n1.Field() != this->field_) {
      throw std::invalid_argument(
          "ProductCode::UnencodeMatrix: codeword field mismatch");
    }
    if (cw_n2_by_n1.Dimension() != 2 || cw_n2_by_n1.Shape(0) != n2_ ||
        cw_n2_by_n1.Shape(1) != n1_) {
      throw std::invalid_argument(
          "ProductCode::UnencodeMatrix: codeword matrix must have shape "
          "{n2, n1}");
    }
    if (!ContainsMatrix(cw_n2_by_n1)) {
      throw std::invalid_argument(
          "ProductCode::UnencodeMatrix: matrix is not a valid product "
          "codeword");
    }
    if (k1_ == 0 || k2_ == 0) {
      return matrix_type::Zeros(this->field_, {k2_, k1_});
    }

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    // Step 1: unencode each column -> {k2, n1}.
    RawArray inter_raw = RawArray::from_shape(ShapeType{k2_, n1_});
    storage_type *inter_data = inter_raw.data();
    for (std::size_t c = 0; c < n1_; ++c) {
      const auto col_cw = ExtractCol(cw_n2_by_n1, c);
      const auto col_msg = col_code_->Unencode(col_cw);
      const storage_type *src = col_msg.RawData().data();
      for (std::size_t r = 0; r < k2_; ++r) {
        inter_data[r * n1_ + c] = src[r];
      }
    }
    auto inter = matrix_type::FromRaw(this->field_, std::move(inter_raw));

    // Step 2: unencode each row -> {k2, k1}.
    RawArray msg_raw = RawArray::from_shape(ShapeType{k2_, k1_});
    storage_type *msg_data = msg_raw.data();
    for (std::size_t r = 0; r < k2_; ++r) {
      const auto row_cw = ExtractRow(inter, r);
      const auto row_msg = row_code_->Unencode(row_cw);
      const storage_type *src = row_msg.RawData().data();
      std::copy(src, src + k1_, msg_data + r * k1_);
    }
    return matrix_type::FromRaw(this->field_, std::move(msg_raw));
  }

  codeword_type FlattenMatrix(const matrix_type &mat) const {
    if (!mat.HasField() || mat.Field() != this->field_ ||
        mat.Dimension() != 2) {
      throw std::invalid_argument(
          "ProductCode::FlattenMatrix: input must be a 2D GFMatrix over the "
          "code's field");
    }
    const std::size_t total = mat.Shape(0) * mat.Shape(1);
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{total});
    if (total > 0) {
      const storage_type *src = mat.RawData().data();
      std::copy(src, src + total, raw.data());
    }
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  matrix_type UnflattenMatrix(const codeword_type &vec, std::size_t rows,
                              std::size_t cols) const {
    if (!vec.HasField() || vec.Field() != this->field_ ||
        vec.Dimension() != 1) {
      throw std::invalid_argument(
          "ProductCode::UnflattenMatrix: input must be a 1D GFArray over the "
          "code's field");
    }
    if (vec.Size() != rows * cols) {
      throw std::invalid_argument(
          "ProductCode::UnflattenMatrix: vector length does not match "
          "rows * cols");
    }
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{rows, cols});
    const std::size_t total = rows * cols;
    if (total > 0) {
      const storage_type *src = vec.RawData().data();
      std::copy(src, src + total, raw.data());
    }
    return matrix_type::FromRaw(this->field_, std::move(raw));
  }

  // Iterative hard row-column algebraic decoder with erasure flagging.
  // Rows that fail component decoding are flagged as erasures for the column
  // pass, and columns that fail are flagged as erasures for the next row pass.
  ProductMatrixDecodingResult<GaloisField> DecodeMatrixHard(
      const matrix_type &received_mat, std::size_t max_iterations = 6) const {
    if (!received_mat.HasField() || received_mat.Field() != this->field_) {
      throw std::invalid_argument(
          "ProductCode::DecodeMatrixHard: field mismatch");
    }
    if (received_mat.Dimension() != 2 || received_mat.Shape(0) != n2_ ||
        received_mat.Shape(1) != n1_) {
      throw std::invalid_argument(
          "ProductCode::DecodeMatrixHard: received matrix must have shape "
          "{n2, n1}");
    }

    if (ContainsMatrix(received_mat)) {
      return ProductMatrixDecodingResult<GaloisField>::Success(
          received_mat, UnencodeMatrix(received_mat), 0, 0);
    }

    matrix_type work = received_mat;
    std::vector<std::size_t> failed_rows;
    std::vector<std::size_t> failed_cols;

    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
      // Row pass: decode each row using failed_cols as erasures.
      std::vector<std::size_t> next_failed_rows;
      for (std::size_t r = 0; r < n2_; ++r) {
        auto row = ExtractRow(work, r);
        auto res =
            DecodeComponentWithOptionalErasures(*row_code_, row, failed_cols);
        if (res.Ok()) {
          WriteRow(work, r, res.codeword);
        } else {
          next_failed_rows.push_back(r);
        }
      }
      failed_rows = std::move(next_failed_rows);

      if (failed_rows.empty() && ContainsMatrix(work)) {
        const std::size_t errs = MatrixHammingDistance(received_mat, work);
        return ProductMatrixDecodingResult<GaloisField>::Success(
            work, UnencodeMatrix(work), errs, iter + 1);
      }

      // Column pass: decode each column using failed_rows as erasures.
      std::vector<std::size_t> next_failed_cols;
      for (std::size_t c = 0; c < n1_; ++c) {
        auto col = ExtractCol(work, c);
        auto res =
            DecodeComponentWithOptionalErasures(*col_code_, col, failed_rows);
        if (res.Ok()) {
          WriteCol(work, c, res.codeword);
        } else {
          next_failed_cols.push_back(c);
        }
      }
      failed_cols = std::move(next_failed_cols);

      if (failed_cols.empty() && ContainsMatrix(work)) {
        const std::size_t errs = MatrixHammingDistance(received_mat, work);
        return ProductMatrixDecodingResult<GaloisField>::Success(
            work, UnencodeMatrix(work), errs, iter + 1);
      }
    }

    return ProductMatrixDecodingResult<GaloisField>::Failure(
        "ProductCode::DecodeMatrixHard: failed to converge to a valid product "
        "codeword within max_iterations",
        max_iterations);
  }

  // Chase-Pyndiah iterative soft-in/soft-out Turbo Product Code (TPC) decoder
  // for binary product codes (BPSK LLR convention: L > 0 favors bit 0, L < 0
  // favors bit 1).
  ChasePyndiahDecodingResult DecodeSoftChasePyndiah(
      const std::vector<double> &channel_llrs_n2_by_n1,
      std::size_t p_least_reliable = 3, std::size_t max_iterations = 4,
      std::vector<double> alpha_schedule = {0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.0,
                                            1.0},
      std::vector<double> beta_schedule = {0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0,
                                           1.0}) const
    requires std::same_as<GaloisField, GaloisFieldBinary>
  {
    if (channel_llrs_n2_by_n1.size() != n2_ * n1_) {
      throw std::invalid_argument(
          "ProductCode::DecodeSoftChasePyndiah: channel LLR vector must have "
          "length n2 * n1");
    }
    if (p_least_reliable > 20) {
      throw std::invalid_argument(
          "ProductCode::DecodeSoftChasePyndiah: p_least_reliable cannot "
          "exceed 20");
    }
    if (alpha_schedule.empty()) {
      alpha_schedule = {0.5};
    }
    if (beta_schedule.empty()) {
      beta_schedule = {1.0};
    }

    const std::size_t total_bits = n2_ * n1_;
    std::vector<double> current_llrs = channel_llrs_n2_by_n1;
    std::vector<double> extrinsic(total_bits, 0.0);

    std::size_t completed_iterations = 0;
    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
      // Half-iteration 2*iter: row Chase-II SISO pass.
      const std::size_t row_half = 2 * iter;
      const double alpha_row =
          alpha_schedule[std::min(row_half, alpha_schedule.size() - 1)];
      const double beta_row =
          beta_schedule[std::min(row_half, beta_schedule.size() - 1)];

      std::vector<double> row_in(n1_);
      for (std::size_t r = 0; r < n2_; ++r) {
        for (std::size_t c = 0; c < n1_; ++c) {
          row_in[c] = current_llrs[r * n1_ + c];
        }
        const auto row_ext = ChasePyndiahComponentSiso(
            *row_code_, row_in, p_least_reliable, beta_row);
        for (std::size_t c = 0; c < n1_; ++c) {
          extrinsic[r * n1_ + c] = row_ext[c];
        }
      }
      for (std::size_t idx = 0; idx < total_bits; ++idx) {
        current_llrs[idx] =
            channel_llrs_n2_by_n1[idx] + alpha_row * extrinsic[idx];
      }

      // Half-iteration 2*iter + 1: column Chase-II SISO pass.
      const std::size_t col_half = 2 * iter + 1;
      const double alpha_col =
          alpha_schedule[std::min(col_half, alpha_schedule.size() - 1)];
      const double beta_col =
          beta_schedule[std::min(col_half, beta_schedule.size() - 1)];

      std::vector<double> col_in(n2_);
      for (std::size_t c = 0; c < n1_; ++c) {
        for (std::size_t r = 0; r < n2_; ++r) {
          col_in[r] = current_llrs[r * n1_ + c];
        }
        const auto col_ext = ChasePyndiahComponentSiso(
            *col_code_, col_in, p_least_reliable, beta_col);
        for (std::size_t r = 0; r < n2_; ++r) {
          extrinsic[r * n1_ + c] = col_ext[r];
        }
      }
      for (std::size_t idx = 0; idx < total_bits; ++idx) {
        current_llrs[idx] =
            channel_llrs_n2_by_n1[idx] + alpha_col * extrinsic[idx];
      }

      completed_iterations = iter + 1;

      // Early stopping if hard decision is a valid product codeword.
      auto hard_mat = HardDecisionMatrix(current_llrs);
      if (ContainsMatrix(hard_mat)) {
        ChasePyndiahDecodingResult res;
        res.codeword_matrix = std::move(hard_mat);
        res.message_matrix = UnencodeMatrix(res.codeword_matrix);
        res.posterior_llrs = std::move(current_llrs);
        res.iterations = completed_iterations;
        res.valid_codeword = true;
        return res;
      }
    }

    ChasePyndiahDecodingResult res;
    res.codeword_matrix = HardDecisionMatrix(current_llrs);
    res.posterior_llrs = std::move(current_llrs);
    res.iterations = completed_iterations;
    res.valid_codeword = ContainsMatrix(res.codeword_matrix);
    if (res.valid_codeword) {
      res.message_matrix = UnencodeMatrix(res.codeword_matrix);
    }
    return res;
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + ", " +
           std::to_string(MinimumDistance()) + "] product code (" +
           std::to_string(n1_) + "x" + std::to_string(n2_) + ") over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 private:
  static std::shared_ptr<GaloisField> ValidateAndGetField(
      const std::shared_ptr<const AbstractLinearCode<GaloisField>> &row_code,
      const std::shared_ptr<const AbstractLinearCode<GaloisField>> &col_code) {
    if (!row_code || !col_code) {
      throw std::invalid_argument(
          "ProductCode: row_code and col_code cannot be null");
    }
    if (row_code->Field() != col_code->Field()) {
      throw std::invalid_argument(
          "ProductCode: row_code and col_code must be defined over the same "
          "field");
    }
    return row_code->Field();
  }

  void RegisterEncodersAndDecoders() {
    this->RegisterEncoder("Product", [](const AbstractCode<GaloisField> &code) {
      const auto &prod = dynamic_cast<const ProductCode<GaloisField> &>(code);
      return std::make_shared<ProductEncoder<GaloisField>>(prod);
    });
    this->RegisterEncoder(
        "GeneratorMatrix", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(lin);
        });
    this->RegisterDecoder(
        "ProductHard", [](const AbstractCode<GaloisField> &code) {
          const auto &prod =
              dynamic_cast<const ProductCode<GaloisField> &>(code);
          return std::make_shared<ProductHardDecoder<GaloisField>>(prod);
        });
    this->RegisterDecoder(
        "Syndrome", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<SyndromeDecoder<GaloisField>>(lin);
        });
  }

  codeword_type ExtractRow(const matrix_type &mat, std::size_t row) const {
    const std::size_t cols = mat.Shape(1);
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{cols});
    const storage_type *src = mat.RawData().data() + row * cols;
    std::copy(src, src + cols, raw.data());
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  codeword_type ExtractCol(const matrix_type &mat, std::size_t col) const {
    const std::size_t rows = mat.Shape(0);
    const std::size_t cols = mat.Shape(1);
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{rows});
    const storage_type *src = mat.RawData().data();
    storage_type *dst = raw.data();
    for (std::size_t r = 0; r < rows; ++r) {
      dst[r] = src[r * cols + col];
    }
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  static void WriteRow(matrix_type &mat, std::size_t row,
                       const codeword_type &vec) {
    const std::size_t cols = mat.Shape(1);
    const storage_type *src = vec.RawData().data();
    storage_type *dst = mat.RawData().data() + row * cols;
    std::copy(src, src + cols, dst);
  }

  static void WriteCol(matrix_type &mat, std::size_t col,
                       const codeword_type &vec) {
    const std::size_t rows = mat.Shape(0);
    const std::size_t cols = mat.Shape(1);
    const storage_type *src = vec.RawData().data();
    storage_type *dst = mat.RawData().data();
    for (std::size_t r = 0; r < rows; ++r) {
      dst[r * cols + col] = src[r];
    }
  }

  static std::size_t MatrixHammingDistance(const matrix_type &a,
                                           const matrix_type &b) {
    const auto sa = a.RawSpan();
    const auto sb = b.RawSpan();
    std::size_t dist = 0;
    for (std::size_t i = 0; i < sa.size(); ++i) {
      if (sa[i] != sb[i]) {
        ++dist;
      }
    }
    return dist;
  }

  static result_type DecodeComponentWithOptionalErasures(
      const AbstractLinearCode<GaloisField> &comp_code,
      const codeword_type &word,
      const std::vector<std::size_t> &erasure_positions) {
    if (!erasure_positions.empty()) {
      auto er_res = comp_code.DecodeWithErasures(word, erasure_positions);
      if (er_res.Ok()) {
        return er_res;
      }
    }
    return comp_code.DecodeWithResult(word);
  }

  matrix_type HardDecisionMatrix(const std::vector<double> &llrs) const
    requires std::same_as<GaloisField, GaloisFieldBinary>
  {
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{n2_, n1_});
    uint8_t *dst = raw.data();
    for (std::size_t i = 0; i < n2_ * n1_; ++i) {
      dst[i] = (llrs[i] < 0.0) ? uint8_t{1} : uint8_t{0};
    }
    return matrix_type::FromRaw(this->field_, std::move(raw));
  }

  // Exact Pyndiah (1998) Chase-II SISO component decoder.
  static std::vector<double> ChasePyndiahComponentSiso(
      const AbstractLinearCode<GaloisFieldBinary> &comp_code,
      const std::vector<double> &y, std::size_t p_least_reliable, double beta)
    requires std::same_as<GaloisField, GaloisFieldBinary>
  {
    const std::size_t n = comp_code.Length();
    const std::size_t p = std::min(p_least_reliable, n);

    // 1. Identify p least-reliable bit positions (smallest |y_j|).
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::partial_sort(order.begin(),
                      order.begin() + static_cast<std::ptrdiff_t>(p),
                      order.end(), [&](std::size_t a, std::size_t b) {
                        return std::abs(y[a]) < std::abs(y[b]);
                      });

    // 2. Hard decision vector h.
    std::vector<uint8_t> hard(n);
    for (std::size_t j = 0; j < n; ++j) {
      hard[j] = (y[j] < 0.0) ? uint8_t{1} : uint8_t{0};
    }

    auto field = comp_code.Field();
    using RawArray = typename GFArray<GaloisFieldBinary>::raw_array_type;
    using ShapeType = typename GFArray<GaloisFieldBinary>::shape_type;

    // 3. Test all 2^p Chase-II patterns via component hard decoder.
    const std::size_t num_patterns = std::size_t{1} << p;
    std::vector<std::vector<uint8_t>> candidates;
    std::vector<double> metrics;
    candidates.reserve(num_patterns);
    metrics.reserve(num_patterns);

    RawArray test_raw = RawArray::from_shape(ShapeType{n});
    uint8_t *test_data = test_raw.data();

    for (std::size_t pat = 0; pat < num_patterns; ++pat) {
      std::copy(hard.begin(), hard.end(), test_data);
      for (std::size_t bit = 0; bit < p; ++bit) {
        if (((pat >> bit) & std::size_t{1}) != 0) {
          test_data[order[bit]] ^= uint8_t{1};
        }
      }
      auto test_word = GFArray<GaloisFieldBinary>::FromRaw(field, test_raw);
      auto dec = comp_code.DecodeWithResult(test_word);
      if (!dec.Ok()) {
        continue;
      }

      const auto cw_span = dec.codeword.RawSpan();
      std::vector<uint8_t> cw(cw_span.begin(), cw_span.end());

      // Compute exact squared Euclidean distance ||y - (1 - 2*C)||^2.
      double dist_sq = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        const double xj = (cw[j] == 0) ? 1.0 : -1.0;
        const double diff = y[j] - xj;
        dist_sq += diff * diff;
      }

      candidates.push_back(std::move(cw));
      metrics.push_back(dist_sq);
    }

    std::vector<double> extrinsic(n, 0.0);
    if (candidates.empty()) {
      return extrinsic;
    }

    // 4. Find best decided codeword D minimizing ||y - D||^2.
    std::size_t best_idx = 0;
    for (std::size_t i = 1; i < metrics.size(); ++i) {
      if (metrics[i] < metrics[best_idx]) {
        best_idx = i;
      }
    }
    const auto &decided = candidates[best_idx];
    const double d_metric = metrics[best_idx];

    // 5. For each bit j, find best competing codeword C_j with C_j[j] != D[j]
    //    and compute exact Pyndiah extrinsic LLR:
    //      r'_j = 0.25 * (||y - C_j||^2 - ||y - D||^2) * (1 - 2*D_j) - y_j
    //    or beta * scale * (1 - 2*D_j) if no competitor exists.
    double mean_abs_y = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      mean_abs_y += std::abs(y[j]);
    }
    mean_abs_y = (n > 0) ? (mean_abs_y / static_cast<double>(n)) : 1.0;
    const double beta_scaled = beta * (mean_abs_y > 0.0 ? mean_abs_y : 1.0);

    for (std::size_t j = 0; j < n; ++j) {
      const double sign_d = (decided[j] == 0) ? 1.0 : -1.0;
      double best_comp_metric = std::numeric_limits<double>::infinity();
      bool found_comp = false;

      for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i][j] != decided[j] && metrics[i] < best_comp_metric) {
          best_comp_metric = metrics[i];
          found_comp = true;
        }
      }

      if (found_comp) {
        extrinsic[j] = 0.25 * (best_comp_metric - d_metric) * sign_d - y[j];
      } else {
        extrinsic[j] = beta_scaled * sign_d;
      }
    }

    return extrinsic;
  }

  std::shared_ptr<const AbstractLinearCode<GaloisField>> row_code_;
  std::shared_ptr<const AbstractLinearCode<GaloisField>> col_code_;
  std::size_t n1_;
  std::size_t k1_;
  std::size_t n2_;
  std::size_t k2_;

  mutable std::mutex matrix_mutex_;
  mutable std::optional<matrix_type> cached_generator_;
  mutable std::optional<matrix_type> cached_parity_check_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_PRODUCT_HPP_
