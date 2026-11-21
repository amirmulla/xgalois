#ifndef XGALOIS_CODING_GRS_HPP_
#define XGALOIS_CODING_GRS_HPP_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class GeneratorMatrixEncoder;
template <typename GaloisField>
class GRSEncoder;
template <typename GaloisField>
class GRSSystematicEncoder;
template <typename GaloisField>
class GRSDecoder;
template <typename GaloisField>
class BerlekampMasseyDecoder;
template <typename GaloisField>
class SugiyamaDecoder;

template <typename GaloisField>
class GeneralizedReedSolomonCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;
  using matrix_type = GFMatrix<GaloisField>;
  using vector_type = GFArray<GaloisField>;
  using polynomial_type = PolynomialDense<GaloisField>;

  GeneralizedReedSolomonCode(std::shared_ptr<GaloisField> field,
                             const vector_type &evaluation_points,
                             const vector_type &column_multipliers,
                             std::size_t dimension)
      : Base(std::move(field), evaluation_points.Size(), dimension,
             Metric::HAMMING, "GRS", "BerlekampMassey"),
        eval_points_(evaluation_points),
        column_multipliers_(column_multipliers) {
    if (!eval_points_.HasField() || eval_points_.Field() != this->field_ ||
        eval_points_.Dimension() != 1) {
      throw std::invalid_argument(
          "GeneralizedReedSolomonCode: evaluation points must be a 1D GFArray "
          "over the code's field");
    }
    if (!column_multipliers_.HasField() ||
        column_multipliers_.Field() != this->field_ ||
        column_multipliers_.Dimension() != 1) {
      throw std::invalid_argument(
          "GeneralizedReedSolomonCode: column multipliers must be a 1D GFArray "
          "over the code's field");
    }
    if (eval_points_.Size() != column_multipliers_.Size()) {
      throw std::invalid_argument(
          "GeneralizedReedSolomonCode: evaluation points and column "
          "multipliers must have the same size");
    }
    if (this->dimension_ == 0 || this->dimension_ > this->length_) {
      throw std::invalid_argument(
          "GeneralizedReedSolomonCode: dimension must satisfy 1 <= k <= n");
    }
    const uint64_t q = this->field_->Order();
    if (q != 0 && this->length_ > static_cast<std::size_t>(q)) {
      throw std::invalid_argument(
          "GeneralizedReedSolomonCode: code length cannot exceed field order");
    }

    const auto eval_span = eval_points_.RawSpan();
    const auto mult_span = column_multipliers_.RawSpan();
    for (std::size_t i = 0; i < this->length_; ++i) {
      if (this->field_->IsZero(mult_span[i])) {
        throw std::invalid_argument(
            "GeneralizedReedSolomonCode: column multipliers must be non-zero");
      }
      if (this->field_->IsZero(eval_span[i])) {
        zero_eval_index_ = i;
      }
      for (std::size_t j = i + 1; j < this->length_; ++j) {
        if (eval_span[i] == eval_span[j]) {
          throw std::invalid_argument(
              "GeneralizedReedSolomonCode: evaluation points must be distinct");
        }
      }
    }

    BuildPrecomputedData();
    RegisterEncodersAndDecoders();
  }

  static GeneralizedReedSolomonCode StandardReedSolomon(
      std::shared_ptr<GaloisField> field, const vector_type &evaluation_points,
      std::size_t dimension) {
    if (!field) {
      throw std::invalid_argument(
          "StandardReedSolomon: field pointer cannot be null");
    }
    auto multipliers = vector_type::Ones(field, {evaluation_points.Size()});
    return GeneralizedReedSolomonCode(std::move(field), evaluation_points,
                                      multipliers, dimension);
  }

  static GeneralizedReedSolomonCode ReedSolomon(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension, const element_type &primitive_element) {
    if (!field) {
      throw std::invalid_argument("ReedSolomon: field pointer cannot be null");
    }
    if (length == 0) {
      throw std::invalid_argument("ReedSolomon: length must be positive");
    }
    const uint64_t q = field->Order();
    if (q != 0 && length > static_cast<std::size_t>(q - 1)) {
      throw std::invalid_argument(
          "ReedSolomon: primitive RS length cannot exceed field->Order() - 1");
    }
    if (!primitive_element.HasField() || primitive_element.Field() != field) {
      throw std::invalid_argument(
          "ReedSolomon: primitive element must belong to the code's field");
    }
    if (primitive_element.IsZero()) {
      throw std::invalid_argument(
          "ReedSolomon: primitive element cannot be zero");
    }

    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    RawArray eval_raw = RawArray::from_shape(ShapeType{length});
    const storage_type alpha = primitive_element.RawValue();
    storage_type power = field->MultiplicativeIdentity();
    for (std::size_t i = 0; i < length; ++i) {
      eval_raw(i) = power;
      if (i + 1 < length) {
        power = field->Mul(power, alpha);
      }
    }

    auto eval_points = vector_type::FromRaw(field, std::move(eval_raw));
    auto multipliers = vector_type::Ones(field, {length});
    return GeneralizedReedSolomonCode(std::move(field), eval_points,
                                      multipliers, dimension);
  }

  static GeneralizedReedSolomonCode ReedSolomon(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension) {
    if (!field) {
      throw std::invalid_argument("ReedSolomon: field pointer cannot be null");
    }
    const auto gen =
        element_type::FromInternal(field->MultiplicativeGenerator(), field);
    return ReedSolomon(std::move(field), length, dimension, gen);
  }

  static GeneralizedReedSolomonCode ExtendedReedSolomon(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension, const element_type &primitive_element) {
    if (!field) {
      throw std::invalid_argument(
          "ExtendedReedSolomon: field pointer cannot be null");
    }
    if (length == 0) {
      throw std::invalid_argument(
          "ExtendedReedSolomon: length must be positive");
    }
    const uint64_t q = field->Order();
    if (q != 0 && length > static_cast<std::size_t>(q)) {
      throw std::invalid_argument(
          "ExtendedReedSolomon: length cannot exceed field->Order()");
    }
    if (!primitive_element.HasField() || primitive_element.Field() != field) {
      throw std::invalid_argument(
          "ExtendedReedSolomon: primitive element must belong to the code's "
          "field");
    }
    if (primitive_element.IsZero()) {
      throw std::invalid_argument(
          "ExtendedReedSolomon: primitive element cannot be zero");
    }

    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    RawArray eval_raw = RawArray::from_shape(ShapeType{length});
    eval_raw(0) = field->AdditiveIdentity();
    const storage_type alpha = primitive_element.RawValue();
    storage_type power = field->MultiplicativeIdentity();
    for (std::size_t i = 1; i < length; ++i) {
      eval_raw(i) = power;
      if (i + 1 < length) {
        power = field->Mul(power, alpha);
      }
    }

    auto eval_points = vector_type::FromRaw(field, std::move(eval_raw));
    auto multipliers = vector_type::Ones(field, {length});
    return GeneralizedReedSolomonCode(std::move(field), eval_points,
                                      multipliers, dimension);
  }

  static GeneralizedReedSolomonCode ExtendedReedSolomon(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension) {
    if (!field) {
      throw std::invalid_argument(
          "ExtendedReedSolomon: field pointer cannot be null");
    }
    const auto gen =
        element_type::FromInternal(field->MultiplicativeGenerator(), field);
    return ExtendedReedSolomon(std::move(field), length, dimension, gen);
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  matrix_type SystematicGeneratorMatrix() const override {
    return linalg::systematic_form(generator_matrix_).G_sys;
  }

  const vector_type &EvaluationPoints() const noexcept { return eval_points_; }

  const vector_type &ColumnMultipliers() const noexcept {
    return column_multipliers_;
  }

  const vector_type &DualMultipliers() const noexcept {
    return dual_multipliers_;
  }

  const std::vector<storage_type> &ColumnMultiplierInverses() const noexcept {
    return column_multiplier_inverses_;
  }

  std::optional<std::size_t> ZeroEvaluationIndex() const noexcept {
    return zero_eval_index_;
  }

  const std::vector<std::size_t> &NonZeroEvalIndices() const noexcept {
    return nonzero_eval_indices_;
  }

  std::span<const storage_type> NonZeroEvalInverses() const noexcept {
    return std::span<const storage_type>(nonzero_eval_inverses_);
  }

  const storage_type &EvalInverseAt(std::size_t index) const {
    return eval_inv_lookup_[index];
  }

  std::size_t MinimumDistance() const override {
    return this->length_ - this->dimension_ + 1;
  }

  bool IsMDS() const noexcept { return true; }

  GeneralizedReedSolomonCode DualCode() const {
    if (this->dimension_ == this->length_) {
      throw std::invalid_argument(
          "DualCode: dual of full-dimension GRS code has dimension 0");
    }
    return GeneralizedReedSolomonCode(this->field_, eval_points_,
                                      dual_multipliers_,
                                      this->length_ - this->dimension_);
  }

  codeword_type EncodePolynomial(const polynomial_type &message_poly) const {
    if (!message_poly.HasField() || message_poly.Field() != this->field_) {
      throw std::invalid_argument(
          "EncodePolynomial: polynomial field mismatch");
    }
    if (!message_poly.IsZero() &&
        static_cast<std::size_t>(message_poly.Degree()) >= this->dimension_) {
      throw std::invalid_argument(
          "EncodePolynomial: message polynomial degree must be less than code "
          "dimension");
    }

    const std::vector<storage_type> evals =
        message_poly.EvalBatch(eval_points_.RawSpan());
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{this->length_});
    this->field_->MulArray(std::span<storage_type>(raw.data(), this->length_),
                           std::span<const storage_type>(evals),
                           column_multipliers_.RawSpan());
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  polynomial_type InterpolateCodeword(const codeword_type &codeword) const {
    if (!codeword.HasField() || codeword.Field() != this->field_) {
      throw std::invalid_argument(
          "InterpolateCodeword: codeword field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != this->length_) {
      throw std::invalid_argument(
          "InterpolateCodeword: codeword must be 1D with size equal to code "
          "length");
    }

    const std::size_t k = this->dimension_;
    std::vector<storage_type> ys(k);
    const auto cw_span = codeword.RawSpan();
    this->field_->MulArray(
        std::span<storage_type>(ys), cw_span.subspan(0, k),
        std::span<const storage_type>(column_multiplier_inverses_.data(), k));
    return utils::InterpolateRaw(this->field_,
                                 eval_points_.RawSpan().subspan(0, k),
                                 std::span<const storage_type>(ys));
  }

  polynomial_type DecodeToPolynomial(const codeword_type &received_word) const {
    const auto corrected = this->DecodeToCode(received_word);
    return InterpolateCodeword(corrected);
  }

  std::vector<std::size_t> SystematicPositions() const {
    std::vector<std::size_t> positions(this->dimension_);
    for (std::size_t i = 0; i < this->dimension_; ++i) {
      positions[i] = i;
    }
    return positions;
  }

  std::string ToString() const override {
    const uint64_t q = this->field_->Order();
    const std::string field_str = (q == 0) ? "2^64" : std::to_string(q);
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + ", " +
           std::to_string(MinimumDistance()) + "] GRS code over GF(" +
           field_str + ")";
  }

 private:
  void BuildPrecomputedData() {
    const std::size_t n = this->length_;
    const std::size_t k = this->dimension_;
    const std::size_t r = n - k;
    const auto eval_span = eval_points_.RawSpan();
    const auto mult_span = column_multipliers_.RawSpan();
    const storage_type zero = this->field_->AdditiveIdentity();
    const storage_type one = this->field_->MultiplicativeIdentity();

    // 1. Precompute column multiplier inverses via Montgomery BatchInv.
    column_multiplier_inverses_.resize(n);
    this->field_->BatchInv(std::span<storage_type>(column_multiplier_inverses_),
                           mult_span);

    // 2. Precompute dual multipliers w_j = (v_j * prod_{l != j}(a_j - a_l))^-1
    //    using a single field inversion via Montgomery BatchInv.
    std::vector<storage_type> denoms(n);
    for (std::size_t j = 0; j < n; ++j) {
      storage_type prod = one;
      const storage_type &aj = eval_span[j];
      for (std::size_t l = 0; l < n; ++l) {
        if (l != j) {
          prod = this->field_->Mul(prod, this->field_->Sub(aj, eval_span[l]));
        }
      }
      denoms[j] = this->field_->Mul(mult_span[j], prod);
    }

    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    RawArray dual_raw = RawArray::from_shape(ShapeType{n});
    this->field_->BatchInv(std::span<storage_type>(dual_raw.data(), n),
                           std::span<const storage_type>(denoms));
    dual_multipliers_ = vector_type::FromRaw(this->field_, std::move(dual_raw));
    const auto dual_span = dual_multipliers_.RawSpan();

    // 3. Build generator matrix G_{i,j} = v_j * alpha_j^i (0^0 = 1).
    RawArray gen_raw = RawArray::from_shape(ShapeType{k, n});
    storage_type *gen_data = gen_raw.data();
    for (std::size_t j = 0; j < n; ++j) {
      const storage_type &aj = eval_span[j];
      const storage_type &vj = mult_span[j];
      storage_type pow = one;
      for (std::size_t i = 0; i < k; ++i) {
        gen_data[i * n + j] = this->field_->Mul(vj, pow);
        if (i + 1 < k) {
          pow = this->field_->Mul(pow, aj);
        }
      }
    }
    generator_matrix_ = matrix_type::FromRaw(this->field_, std::move(gen_raw));

    // 4. Build parity-check matrix H_{i,j} = w_j * alpha_j^i (0^0 = 1).
    if (r == 0) {
      parity_check_matrix_ = matrix_type::Zeros(this->field_, {0, n});
    } else {
      RawArray h_raw = RawArray::from_shape(ShapeType{r, n});
      storage_type *h_data = h_raw.data();
      for (std::size_t j = 0; j < n; ++j) {
        const storage_type &aj = eval_span[j];
        const storage_type &wj = dual_span[j];
        storage_type pow = one;
        for (std::size_t i = 0; i < r; ++i) {
          h_data[i * n + j] = this->field_->Mul(wj, pow);
          if (i + 1 < r) {
            pow = this->field_->Mul(pow, aj);
          }
        }
      }
      parity_check_matrix_ =
          matrix_type::FromRaw(this->field_, std::move(h_raw));
    }

    // 5. Cache non-zero evaluation point inverses (1 BatchInv) for fast Chien
    //    search and Forney evaluation.
    nonzero_eval_indices_.clear();
    nonzero_eval_indices_.reserve(n);
    std::vector<storage_type> nonzero_vals;
    nonzero_vals.reserve(n);
    for (std::size_t j = 0; j < n; ++j) {
      if (!this->field_->IsZero(eval_span[j])) {
        nonzero_eval_indices_.push_back(j);
        nonzero_vals.push_back(eval_span[j]);
      }
    }
    nonzero_eval_inverses_.resize(nonzero_vals.size());
    this->field_->BatchInv(std::span<storage_type>(nonzero_eval_inverses_),
                           std::span<const storage_type>(nonzero_vals));
    eval_inv_lookup_.assign(n, zero);
    for (std::size_t idx = 0; idx < nonzero_eval_indices_.size(); ++idx) {
      eval_inv_lookup_[nonzero_eval_indices_[idx]] =
          nonzero_eval_inverses_[idx];
    }
  }

  void RegisterEncodersAndDecoders();

  vector_type eval_points_;
  vector_type column_multipliers_;
  std::vector<storage_type> column_multiplier_inverses_;
  vector_type dual_multipliers_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;
  std::optional<std::size_t> zero_eval_index_;
  std::vector<std::size_t> nonzero_eval_indices_;
  std::vector<storage_type> nonzero_eval_inverses_;
  std::vector<storage_type> eval_inv_lookup_;
};

}  // namespace coding
}  // namespace xg

// clang-format off
#include "xgalois/coding/decoder/grs_decoder.hpp"
#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/decoder/sugiyama_decoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/encoder/grs_encoder.hpp"
// clang-format on

namespace xg {
namespace coding {

template <typename GaloisField>
inline void
GeneralizedReedSolomonCode<GaloisField>::RegisterEncodersAndDecoders() {
  this->RegisterEncoder("GRS", [](const AbstractCode<GaloisField> &code) {
    return std::make_shared<GRSEncoder<GaloisField>>(code);
  });
  this->RegisterEncoder(
      "GRSSystematic", [](const AbstractCode<GaloisField> &code) {
        return std::make_shared<GRSSystematicEncoder<GaloisField>>(code);
      });
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<GaloisField> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(lin);
      });

  this->RegisterDecoder(
      "BerlekampMassey", [](const AbstractCode<GaloisField> &code) {
        return std::make_shared<BerlekampMasseyDecoder<GaloisField>>(code);
      });
  this->RegisterDecoder("Sugiyama", [](const AbstractCode<GaloisField> &code) {
    return std::make_shared<SugiyamaDecoder<GaloisField>>(code);
  });
  this->RegisterDecoder("PGZ", [](const AbstractCode<GaloisField> &code) {
    return std::make_shared<GRSDecoder<GaloisField>>(code);
  });
  this->RegisterDecoder("GRS", [](const AbstractCode<GaloisField> &code) {
    return std::make_shared<BerlekampMasseyDecoder<GaloisField>>(code);
  });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_GRS_HPP_
