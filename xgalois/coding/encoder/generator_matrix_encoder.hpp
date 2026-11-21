#ifndef XGALOIS_CODING_ENCODER_GENERATOR_MATRIX_ENCODER_HPP_
#define XGALOIS_CODING_ENCODER_GENERATOR_MATRIX_ENCODER_HPP_

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class GeneratorMatrixEncoder : public Encoder<GaloisField> {
 public:
  using Base = Encoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = GFMatrix<GaloisField>;

  explicit GeneratorMatrixEncoder(const AbstractLinearCode<GaloisField> &code)
      : Base(code),
        linear_code_(&code),
        k_(code.Dimension()),
        n_(code.Length()) {
    if (k_ == 0) {
      return;
    }

    const auto &G = linear_code_->GeneratorMatrix();
    if (!G.HasField() || G.Field() != this->code_->Field()) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder: generator matrix field mismatch");
    }
    if (G.Dimension() != 2 || G.Shape(0) != k_ || G.Shape(1) != n_) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder: generator matrix shape must be (k, n)");
    }

    auto sys = linalg::systematic_form(G);
    if (sys.pivot_cols.size() != k_) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder: generator matrix must have full row rank");
    }
    pivot_cols_ = std::move(sys.pivot_cols);

    // Check fast path: G is already systematic with leading identity [I_k | P].
    is_leading_identity_ = linear_code_->IsSystematic();
    if (!is_leading_identity_) {
      // Extract the k x k invertible submatrix G_P formed by the pivot columns
      // of the original G, and precompute G_P^{-1} so Unencode runs in O(k^2):
      //   c = m * G  =>  c_P = m * G_P  =>  m = c_P * G_P^{-1}.
      using RawArray = typename matrix_type::raw_array_type;
      using ShapeType = typename matrix_type::shape_type;
      RawArray gp_raw = RawArray::from_shape(ShapeType{k_, k_});
      const storage_type *g_data = G.RawData().data();
      storage_type *gp_data = gp_raw.data();
      for (std::size_t r = 0; r < k_; ++r) {
        const storage_type *g_row = g_data + r * n_;
        storage_type *gp_row = gp_data + r * k_;
        for (std::size_t j = 0; j < k_; ++j) {
          gp_row[j] = g_row[pivot_cols_[j]];
        }
      }
      matrix_type G_P =
          matrix_type::FromRaw(this->code_->Field(), std::move(gp_raw));
      inv_gp_ = linalg::inv(G_P);
    }
  }

  codeword_type Encode(const message_type &message) const override {
    auto field = this->code_->Field();
    if (!message.HasField() || message.Field() != field) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder::Encode: message field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != k_) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder::Encode: message must be 1D with length "
          "equal to code dimension");
    }
    if (k_ == 0) {
      return codeword_type::Zeros(field, {n_});
    }
    return linalg::dot(message, linear_code_->GeneratorMatrix());
  }

  message_type Unencode(const codeword_type &codeword) const override {
    auto field = this->code_->Field();
    if (!codeword.HasField() || codeword.Field() != field) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder::Unencode: codeword field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != n_) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder::Unencode: codeword must be 1D with length "
          "equal to code length");
    }
    if (!linear_code_->Contains(codeword)) {
      throw std::invalid_argument(
          "GeneratorMatrixEncoder::Unencode: input is not a valid codeword");
    }
    if (k_ == 0) {
      return message_type::Zeros(field, {0});
    }

    using RawArray = typename message_type::raw_array_type;
    using ShapeType = typename message_type::shape_type;
    RawArray c_p_raw = RawArray::from_shape(ShapeType{k_});
    const storage_type *cw_data = codeword.RawData().data();
    storage_type *cp_data = c_p_raw.data();

    if (is_leading_identity_) {
      for (std::size_t j = 0; j < k_; ++j) {
        cp_data[j] = cw_data[j];
      }
      return message_type::FromRaw(std::move(field), std::move(c_p_raw));
    }

    for (std::size_t j = 0; j < k_; ++j) {
      cp_data[j] = cw_data[pivot_cols_[j]];
    }
    message_type c_p =
        message_type::FromRaw(std::move(field), std::move(c_p_raw));
    return linalg::dot(c_p, inv_gp_);
  }

  std::string Name() const override { return "GeneratorMatrix"; }

  std::size_t MessageLength() const noexcept { return k_; }

  std::string ToString() const override {
    return "GeneratorMatrix encoder for [" + std::to_string(n_) + ", " +
           std::to_string(k_) + "] linear code";
  }

 private:
  const AbstractLinearCode<GaloisField> *linear_code_;
  std::size_t k_;
  std::size_t n_;
  bool is_leading_identity_ = false;
  std::vector<std::size_t> pivot_cols_;
  matrix_type inv_gp_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ENCODER_GENERATOR_MATRIX_ENCODER_HPP_
