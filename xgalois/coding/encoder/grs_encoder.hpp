#ifndef XGALOIS_CODING_ENCODER_GRS_ENCODER_HPP_
#define XGALOIS_CODING_ENCODER_GRS_ENCODER_HPP_

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace coding {

// Polynomial-evaluation encoder for Generalized Reed-Solomon codes:
//   c_j = v_j * m(alpha_j),  where m(x) = sum_{i=0}^{k-1} m_i x^i.
template <typename GaloisField>
class GRSEncoder : public Encoder<GaloisField> {
 public:
  using Base = Encoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using polynomial_type = PolynomialDense<GaloisField>;

  explicit GRSEncoder(const AbstractCode<GaloisField> &code) : Base(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSEncoder requires a GeneralizedReedSolomonCode");
    }
  }

  codeword_type Encode(const message_type &message) const override {
    const auto field = grs_code_->Field();
    if (!message.HasField() || message.Field() != field) {
      throw std::invalid_argument("GRSEncoder::Encode: field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != grs_code_->Dimension()) {
      throw std::invalid_argument(
          "GRSEncoder::Encode: message must be 1D with length equal to code "
          "dimension");
    }

    const auto msg_span = message.RawSpan();
    std::vector<storage_type> coeffs(msg_span.begin(), msg_span.end());
    polynomial_type message_poly =
        polynomial_type::FromRaw(std::move(coeffs), field);
    return grs_code_->EncodePolynomial(message_poly);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    const auto field = grs_code_->Field();
    if (!codeword.HasField() || codeword.Field() != field) {
      throw std::invalid_argument("GRSEncoder::Unencode: field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != grs_code_->Length()) {
      throw std::invalid_argument(
          "GRSEncoder::Unencode: codeword must be 1D with length equal to code "
          "length");
    }
    if (!grs_code_->Contains(codeword)) {
      throw std::invalid_argument(
          "GRSEncoder::Unencode: input is not a valid codeword");
    }

    const polynomial_type poly = grs_code_->InterpolateCodeword(codeword);
    const std::size_t k = grs_code_->Dimension();
    using RawArray = typename message_type::raw_array_type;
    using ShapeType = typename message_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{k});
    for (std::size_t i = 0; i < k; ++i) {
      raw(i) = poly.RawCoeff(i);
    }
    return message_type::FromRaw(field, std::move(raw));
  }

  codeword_type EncodePolynomial(const polynomial_type &message_poly) const {
    return grs_code_->EncodePolynomial(message_poly);
  }

  std::size_t MessageLength() const noexcept { return grs_code_->Dimension(); }

  std::string Name() const override { return "GRS"; }

 private:
  const GeneralizedReedSolomonCode<GaloisField> *grs_code_;
};

// Systematic encoder for Generalized Reed-Solomon codes:
// Places message symbols verbatim at coordinates 0..k-1, interpolates the
// unique degree < k polynomial f satisfying v_j * f(alpha_j) = m_j (j < k) in
// O(k^2) via utils::InterpolateRaw, and evaluates parity coordinates j >= k as
// c_j = v_j * f(alpha_j).
template <typename GaloisField>
class GRSSystematicEncoder : public Encoder<GaloisField> {
 public:
  using Base = Encoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using polynomial_type = PolynomialDense<GaloisField>;

  explicit GRSSystematicEncoder(const AbstractCode<GaloisField> &code)
      : Base(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSSystematicEncoder requires a GeneralizedReedSolomonCode");
    }
  }

  codeword_type Encode(const message_type &message) const override {
    const auto field = grs_code_->Field();
    const std::size_t k = grs_code_->Dimension();
    const std::size_t n = grs_code_->Length();
    if (!message.HasField() || message.Field() != field) {
      throw std::invalid_argument(
          "GRSSystematicEncoder::Encode: field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != k) {
      throw std::invalid_argument(
          "GRSSystematicEncoder::Encode: message must be 1D with length equal "
          "to code dimension");
    }

    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray cw_raw = RawArray::from_shape(ShapeType{n});
    const auto msg_span = message.RawSpan();
    std::copy(msg_span.begin(), msg_span.end(), cw_raw.data());

    if (k < n) {
      const auto &mult_invs = grs_code_->ColumnMultiplierInverses();
      std::vector<storage_type> ys(k);
      field->MulArray(std::span<storage_type>(ys), msg_span,
                      std::span<const storage_type>(mult_invs.data(), k));

      const auto eval_span = grs_code_->EvaluationPoints().RawSpan();
      const polynomial_type poly = utils::InterpolateRaw(
          field, eval_span.subspan(0, k), std::span<const storage_type>(ys));

      const std::size_t r = n - k;
      const std::vector<storage_type> parity_evals =
          poly.EvalBatch(eval_span.subspan(k, r));
      field->MulArray(std::span<storage_type>(cw_raw.data() + k, r),
                      std::span<const storage_type>(parity_evals),
                      grs_code_->ColumnMultipliers().RawSpan().subspan(k, r));
    }

    return codeword_type::FromRaw(field, std::move(cw_raw));
  }

  message_type Unencode(const codeword_type &codeword) const override {
    const auto field = grs_code_->Field();
    const std::size_t k = grs_code_->Dimension();
    const std::size_t n = grs_code_->Length();
    if (!codeword.HasField() || codeword.Field() != field) {
      throw std::invalid_argument(
          "GRSSystematicEncoder::Unencode: field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != n) {
      throw std::invalid_argument(
          "GRSSystematicEncoder::Unencode: codeword must be 1D with length "
          "equal to code length");
    }
    if (!grs_code_->Contains(codeword)) {
      throw std::invalid_argument(
          "GRSSystematicEncoder::Unencode: input is not a valid codeword");
    }

    using RawArray = typename message_type::raw_array_type;
    using ShapeType = typename message_type::shape_type;
    RawArray msg_raw = RawArray::from_shape(ShapeType{k});
    const auto cw_span = codeword.RawSpan();
    std::copy(cw_span.begin(), cw_span.begin() + static_cast<std::ptrdiff_t>(k),
              msg_raw.data());
    return message_type::FromRaw(field, std::move(msg_raw));
  }

  std::size_t MessageLength() const noexcept { return grs_code_->Dimension(); }

  std::string Name() const override { return "GRSSystematic"; }

 private:
  const GeneralizedReedSolomonCode<GaloisField> *grs_code_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ENCODER_GRS_ENCODER_HPP_
