#ifndef XGALOIS_CODING_CONCATENATED_HPP_
#define XGALOIS_CODING_CONCATENATED_HPP_

#include <algorithm>
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
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/gabidulin.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename BaseField, typename ExtField>
class ConcatenatedCode;

template <typename BaseField, typename ExtField>
class ConcatenatedEncoder : public Encoder<BaseField> {
 public:
  using Base = Encoder<BaseField>;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;

  explicit ConcatenatedEncoder(
      const ConcatenatedCode<BaseField, ExtField> &code)
      : Base(code), concat_code_(&code) {}

  codeword_type Encode(const message_type &message) const override {
    return concat_code_->EncodeFast(message);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    return concat_code_->UnencodeFast(codeword);
  }

  std::string Name() const override { return "Concatenated"; }

 private:
  const ConcatenatedCode<BaseField, ExtField> *concat_code_;
};

template <typename BaseField, typename ExtField>
class ConcatenatedStandardDecoder : public Decoder<BaseField> {
 public:
  using Base = Decoder<BaseField>;
  using codeword_type = typename Base::codeword_type;
  using result_type = typename Base::result_type;

  explicit ConcatenatedStandardDecoder(
      const ConcatenatedCode<BaseField, ExtField> &code)
      : Base(code), concat_code_(&code) {}

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return concat_code_->DecodeStandard(received_word);
  }

  std::string Name() const override { return "ConcatenatedStandard"; }

  std::size_t DecodingRadius() const override {
    const std::size_t t_in =
        (concat_code_->InnerCode()->MinimumDistance() > 0)
            ? (concat_code_->InnerCode()->MinimumDistance() - 1) / 2
            : 0;
    const std::size_t t_out =
        (concat_code_->OuterCode()->MinimumDistance() > 0)
            ? (concat_code_->OuterCode()->MinimumDistance() - 1) / 2
            : 0;
    return (t_out + 1) * (t_in + 1) - 1;
  }

 private:
  const ConcatenatedCode<BaseField, ExtField> *concat_code_;
};

template <typename BaseField, typename ExtField>
class ConcatenatedGmdDecoder : public Decoder<BaseField> {
 public:
  using Base = Decoder<BaseField>;
  using codeword_type = typename Base::codeword_type;
  using result_type = typename Base::result_type;

  explicit ConcatenatedGmdDecoder(
      const ConcatenatedCode<BaseField, ExtField> &code)
      : Base(code), concat_code_(&code) {}

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return concat_code_->DecodeGmdHard(received_word);
  }

  std::string Name() const override { return "ConcatenatedGMD"; }

  std::size_t DecodingRadius() const override {
    const std::size_t d = concat_code_->DesignedDistance();
    return (d > 0) ? (d - 1) / 2 : 0;
  }

 private:
  const ConcatenatedCode<BaseField, ExtField> *concat_code_;
};

template <typename BaseField, typename ExtField>
class ConcatenatedCode : public AbstractLinearCode<BaseField> {
 public:
  using Base = AbstractLinearCode<BaseField>;
  using base_field_type = BaseField;
  using ext_field_type = ExtField;
  using base_storage_type = typename BaseField::element_type;
  using ext_storage_type = typename ExtField::element_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = typename Base::matrix_type;
  using result_type = typename Base::result_type;

  ConcatenatedCode(
      std::shared_ptr<const AbstractLinearCode<ExtField>> outer_code,
      std::shared_ptr<const AbstractLinearCode<BaseField>> inner_code)
      : Base(ValidateAndGetBaseField(outer_code, inner_code),
             outer_code->Length() * inner_code->Length(),
             outer_code->Dimension() * inner_code->Dimension(), Metric::HAMMING,
             "Concatenated", "ConcatenatedStandard"),
        outer_code_(std::move(outer_code)),
        inner_code_(std::move(inner_code)),
        tower_(inner_code_->Field(), outer_code_->Field()),
        n_out_(outer_code_->Length()),
        k_out_(outer_code_->Dimension()),
        n_in_(inner_code_->Length()),
        k_in_(inner_code_->Dimension()),
        rel_deg_(tower_.RelativeDegree()) {
    if (k_in_ == 0 || (k_in_ % rel_deg_) != 0) {
      throw std::invalid_argument(
          "ConcatenatedCode: inner code dimension k_in must be a positive "
          "multiple of the relative field extension degree m");
    }
    streams_ = k_in_ / rel_deg_;
    RegisterEncodersAndDecoders();
  }

  ConcatenatedCode(const ConcatenatedCode &other)
      : Base(other),
        outer_code_(other.outer_code_),
        inner_code_(other.inner_code_),
        tower_(other.tower_),
        n_out_(other.n_out_),
        k_out_(other.k_out_),
        n_in_(other.n_in_),
        k_in_(other.k_in_),
        rel_deg_(other.rel_deg_),
        streams_(other.streams_) {
    std::lock_guard<std::mutex> lock(other.matrix_mutex_);
    cached_generator_ = other.cached_generator_;
    cached_parity_check_ = other.cached_parity_check_;
  }

  ConcatenatedCode &operator=(const ConcatenatedCode &other) {
    if (this != &other) {
      Base::operator=(other);
      outer_code_ = other.outer_code_;
      inner_code_ = other.inner_code_;
      tower_ = other.tower_;
      n_out_ = other.n_out_;
      k_out_ = other.k_out_;
      n_in_ = other.n_in_;
      k_in_ = other.k_in_;
      rel_deg_ = other.rel_deg_;
      streams_ = other.streams_;
      std::scoped_lock lock(matrix_mutex_, other.matrix_mutex_);
      cached_generator_ = other.cached_generator_;
      cached_parity_check_ = other.cached_parity_check_;
    }
    return *this;
  }

  ConcatenatedCode(ConcatenatedCode &&other) noexcept
      : Base(std::move(other)),
        outer_code_(std::move(other.outer_code_)),
        inner_code_(std::move(other.inner_code_)),
        tower_(std::move(other.tower_)),
        n_out_(other.n_out_),
        k_out_(other.k_out_),
        n_in_(other.n_in_),
        k_in_(other.k_in_),
        rel_deg_(other.rel_deg_),
        streams_(other.streams_),
        cached_generator_(std::move(other.cached_generator_)),
        cached_parity_check_(std::move(other.cached_parity_check_)) {}

  ConcatenatedCode &operator=(ConcatenatedCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      outer_code_ = std::move(other.outer_code_);
      inner_code_ = std::move(other.inner_code_);
      tower_ = std::move(other.tower_);
      n_out_ = other.n_out_;
      k_out_ = other.k_out_;
      n_in_ = other.n_in_;
      k_in_ = other.k_in_;
      rel_deg_ = other.rel_deg_;
      streams_ = other.streams_;
      std::lock_guard<std::mutex> lock(matrix_mutex_);
      cached_generator_ = std::move(other.cached_generator_);
      cached_parity_check_ = std::move(other.cached_parity_check_);
    }
    return *this;
  }

  ~ConcatenatedCode() override = default;

  std::shared_ptr<const AbstractLinearCode<ExtField>> OuterCode() const {
    return outer_code_;
  }
  std::shared_ptr<const AbstractLinearCode<BaseField>> InnerCode() const {
    return inner_code_;
  }
  const FieldTowerAdapter<BaseField, ExtField> &Tower() const noexcept {
    return tower_;
  }
  std::size_t OuterLength() const noexcept { return n_out_; }
  std::size_t OuterDimension() const noexcept { return k_out_; }
  std::size_t InnerLength() const noexcept { return n_in_; }
  std::size_t InnerDimension() const noexcept { return k_in_; }
  std::size_t RelativeDegree() const noexcept { return rel_deg_; }
  std::size_t InterleavedStreams() const noexcept { return streams_; }

  std::size_t DesignedDistance() const {
    return outer_code_->MinimumDistance() * inner_code_->MinimumDistance();
  }

  // Returns the guaranteed designed minimum distance d_out * d_in.
  std::size_t MinimumDistance() const override {
    if (this->cached_min_distance_.has_value()) {
      return *this->cached_min_distance_;
    }
    const std::size_t d = DesignedDistance();
    this->cached_min_distance_ = d;
    return d;
  }

  // Lazy base-field generator matrix G of shape {k_out * k_in, n_out * n_in}.
  const matrix_type &GeneratorMatrix() const override {
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!cached_generator_.has_value()) {
      const std::size_t k_total = this->dimension_;
      const std::size_t n_total = this->length_;
      using RawArray = typename matrix_type::raw_array_type;
      using ShapeType = typename matrix_type::shape_type;
      RawArray g_raw = RawArray::from_shape(ShapeType{k_total, n_total});
      base_storage_type *g_data = g_raw.data();
      const base_storage_type zero = this->field_->AdditiveIdentity();
      const base_storage_type one = this->field_->MultiplicativeIdentity();

      for (std::size_t i = 0; i < k_total; ++i) {
        auto unit = message_type::Zeros(this->field_, {k_total});
        unit.Raw(i) = one;
        auto cw = EncodeFast(unit);
        const base_storage_type *src = cw.RawData().data();
        std::copy(src, src + n_total, g_data + i * n_total);
      }
      (void)zero;
      cached_generator_ = matrix_type::FromRaw(this->field_, std::move(g_raw));
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

  // Fast base-field encoding: pack base message -> outer ExtField symbols ->
  // encode with outer_code -> expand each outer symbol -> encode each block
  // with inner_code.
  codeword_type EncodeFast(const message_type &message) const {
    if (!message.HasField() || message.Field() != this->field_) {
      throw std::invalid_argument(
          "ConcatenatedCode::EncodeFast: message field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != this->dimension_) {
      throw std::invalid_argument(
          "ConcatenatedCode::EncodeFast: message must be 1D with length "
          "k_out * k_in");
    }
    if (this->dimension_ == 0) {
      return codeword_type::Zeros(this->field_, {this->length_});
    }

    // 1. Pack base-field message symbols into `streams_` outer messages of
    //    length k_out over ExtField.
    std::vector<GFArray<ExtField>> outer_codewords;
    outer_codewords.reserve(streams_);

    const base_storage_type *msg_data = message.RawData().data();
    using ExtRawArray = typename GFArray<ExtField>::raw_array_type;
    using ExtShapeType = typename GFArray<ExtField>::shape_type;
    using BaseRawArray = typename codeword_type::raw_array_type;
    using BaseShapeType = typename codeword_type::shape_type;

    for (std::size_t s = 0; s < streams_; ++s) {
      ExtRawArray ext_msg_raw = ExtRawArray::from_shape(ExtShapeType{k_out_});
      ext_storage_type *ext_dst = ext_msg_raw.data();
      for (std::size_t i = 0; i < k_out_; ++i) {
        const std::size_t base_offset = (s * k_out_ + i) * rel_deg_;
        BaseRawArray chunk_raw =
            BaseRawArray::from_shape(BaseShapeType{rel_deg_});
        std::copy(msg_data + base_offset, msg_data + base_offset + rel_deg_,
                  chunk_raw.data());
        auto chunk = codeword_type::FromRaw(this->field_, std::move(chunk_raw));
        ext_dst[i] = tower_.ContractElement(chunk);
      }
      auto ext_msg = GFArray<ExtField>::FromRaw(outer_code_->Field(),
                                                std::move(ext_msg_raw));
      outer_codewords.push_back(outer_code_->Encode(ext_msg));
    }

    // 2. Expand each outer position across streams and encode with inner_code.
    return EncodeOuterCodewords(outer_codewords);
  }

  message_type UnencodeFast(const codeword_type &codeword) const {
    if (!codeword.HasField() || codeword.Field() != this->field_) {
      throw std::invalid_argument(
          "ConcatenatedCode::UnencodeFast: codeword field mismatch");
    }
    if (codeword.Dimension() != 1 || codeword.Size() != this->length_) {
      throw std::invalid_argument(
          "ConcatenatedCode::UnencodeFast: codeword must be 1D with length "
          "n_out * n_in");
    }
    if (this->dimension_ == 0) {
      if (!Contains(codeword)) {
        throw std::invalid_argument(
            "ConcatenatedCode::UnencodeFast: input is not a valid codeword");
      }
      return message_type::Zeros(this->field_, {0});
    }

    using ExtRawArray = typename GFArray<ExtField>::raw_array_type;
    using ExtShapeType = typename GFArray<ExtField>::shape_type;
    std::vector<ExtRawArray> outer_raws;
    outer_raws.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      outer_raws.push_back(ExtRawArray::from_shape(ExtShapeType{n_out_}));
    }

    for (std::size_t j = 0; j < n_out_; ++j) {
      auto block = ExtractBlock(codeword, j);
      if (!inner_code_->Contains(block)) {
        throw std::invalid_argument(
            "ConcatenatedCode::UnencodeFast: inner block is not a valid "
            "codeword");
      }
      auto inner_msg = inner_code_->Unencode(block);
      const base_storage_type *in_data = inner_msg.RawData().data();
      for (std::size_t s = 0; s < streams_; ++s) {
        auto chunk = SliceBase(in_data + s * rel_deg_, rel_deg_);
        outer_raws[s](j) = tower_.ContractElement(chunk);
      }
    }

    std::vector<GFArray<ExtField>> outer_msgs;
    outer_msgs.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      auto outer_cw = GFArray<ExtField>::FromRaw(outer_code_->Field(),
                                                 std::move(outer_raws[s]));
      if (!outer_code_->Contains(outer_cw)) {
        throw std::invalid_argument(
            "ConcatenatedCode::UnencodeFast: outer stream is not a valid "
            "codeword");
      }
      outer_msgs.push_back(outer_code_->Unencode(outer_cw));
    }

    return ExpandOuterMessages(outer_msgs);
  }

  bool Contains(const codeword_type &word) const override {
    if (!word.HasField() || word.Field() != this->field_ ||
        word.Dimension() != 1 || word.Size() != this->length_) {
      return false;
    }
    using ExtRawArray = typename GFArray<ExtField>::raw_array_type;
    using ExtShapeType = typename GFArray<ExtField>::shape_type;
    std::vector<ExtRawArray> outer_raws;
    outer_raws.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      outer_raws.push_back(ExtRawArray::from_shape(ExtShapeType{n_out_}));
    }

    for (std::size_t j = 0; j < n_out_; ++j) {
      auto block = ExtractBlock(word, j);
      if (!inner_code_->Contains(block)) {
        return false;
      }
      auto inner_msg = inner_code_->Unencode(block);
      const base_storage_type *in_data = inner_msg.RawData().data();
      for (std::size_t s = 0; s < streams_; ++s) {
        auto chunk = SliceBase(in_data + s * rel_deg_, rel_deg_);
        outer_raws[s](j) = tower_.ContractElement(chunk);
      }
    }

    for (std::size_t s = 0; s < streams_; ++s) {
      auto outer_cw = GFArray<ExtField>::FromRaw(outer_code_->Field(),
                                                 std::move(outer_raws[s]));
      if (!outer_code_->Contains(outer_cw)) {
        return false;
      }
    }
    return true;
  }

  // Standard Inner -> Outer algebraic decoder: decodes each inner block,
  // flags failed inner blocks as outer erasures, and invokes
  // outer_code->DecodeWithErasures.
  result_type DecodeStandard(const codeword_type &received_base) const {
    ValidateReceivedWord(received_base, "ConcatenatedCode::DecodeStandard");

    using ExtRawArray = typename GFArray<ExtField>::raw_array_type;
    using ExtShapeType = typename GFArray<ExtField>::shape_type;
    const ext_storage_type ext_zero = outer_code_->Field()->AdditiveIdentity();

    std::vector<ExtRawArray> outer_raws;
    outer_raws.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      ExtRawArray arr = ExtRawArray::from_shape(ExtShapeType{n_out_});
      arr.fill(ext_zero);
      outer_raws.push_back(std::move(arr));
    }

    std::vector<std::size_t> outer_erasures;
    for (std::size_t j = 0; j < n_out_; ++j) {
      auto block = ExtractBlock(received_base, j);
      auto inner_res = inner_code_->DecodeWithResult(block);
      if (!inner_res.Ok()) {
        outer_erasures.push_back(j);
        continue;
      }
      auto inner_msg =
          (inner_res.message.HasField() && inner_res.message.Size() == k_in_)
              ? std::move(inner_res.message)
              : inner_code_->Unencode(inner_res.codeword);
      const base_storage_type *in_data = inner_msg.RawData().data();
      for (std::size_t s = 0; s < streams_; ++s) {
        auto chunk = SliceBase(in_data + s * rel_deg_, rel_deg_);
        outer_raws[s](j) = tower_.ContractElement(chunk);
      }
    }

    std::vector<GFArray<ExtField>> decoded_outer_cws;
    std::vector<GFArray<ExtField>> decoded_outer_msgs;
    decoded_outer_cws.reserve(streams_);
    decoded_outer_msgs.reserve(streams_);

    for (std::size_t s = 0; s < streams_; ++s) {
      auto outer_recv = GFArray<ExtField>::FromRaw(outer_code_->Field(),
                                                   std::move(outer_raws[s]));
      auto outer_res =
          outer_code_->DecodeWithErasures(outer_recv, outer_erasures);
      if (!outer_res.Ok()) {
        return result_type::Failure(
            "ConcatenatedCode::DecodeStandard: outer decoding failed: " +
            outer_res.error_message);
      }
      auto outer_msg =
          (outer_res.message.HasField() && outer_res.message.Size() == k_out_)
              ? std::move(outer_res.message)
              : outer_code_->Unencode(outer_res.codeword);
      decoded_outer_cws.push_back(std::move(outer_res.codeword));
      decoded_outer_msgs.push_back(std::move(outer_msg));
    }

    auto full_cw = EncodeOuterCodewords(decoded_outer_cws);
    auto full_msg = ExpandOuterMessages(decoded_outer_msgs);
    const std::size_t base_errors = this->Distance(received_base, full_cw);
    return result_type::Success(std::move(full_cw), std::move(full_msg),
                                base_errors, outer_erasures.size());
  }

  // Forney Generalized Minimum Distance (GMD) decoder with externally supplied
  // normalized inner reliability scores theta_0, ..., theta_{n_out - 1} in
  // [0, 1]. Runs floor((d_out - 1) / 2) + 1 outer erasure trials erasing the
  // least-reliable symbols and selects the valid codeword minimizing base-field
  // Hamming distance to received_base.
  result_type DecodeGmd(const codeword_type &received_base,
                        const std::vector<double> &inner_reliabilities) const {
    ValidateReceivedWord(received_base, "ConcatenatedCode::DecodeGmd");
    if (inner_reliabilities.size() != n_out_) {
      throw std::invalid_argument(
          "ConcatenatedCode::DecodeGmd: inner_reliabilities must have length "
          "n_out");
    }

    using ExtRawArray = typename GFArray<ExtField>::raw_array_type;
    using ExtShapeType = typename GFArray<ExtField>::shape_type;
    const ext_storage_type ext_zero = outer_code_->Field()->AdditiveIdentity();

    std::vector<ExtRawArray> outer_raws;
    outer_raws.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      ExtRawArray arr = ExtRawArray::from_shape(ExtShapeType{n_out_});
      arr.fill(ext_zero);
      outer_raws.push_back(std::move(arr));
    }

    std::vector<bool> failed(n_out_, false);
    std::vector<double> reliabilities(n_out_, 0.0);
    std::size_t num_failed = 0;

    for (std::size_t j = 0; j < n_out_; ++j) {
      reliabilities[j] = std::clamp(inner_reliabilities[j], 0.0, 1.0);
      auto block = ExtractBlock(received_base, j);
      auto inner_res = inner_code_->DecodeWithResult(block);
      if (!inner_res.Ok()) {
        failed[j] = true;
        reliabilities[j] = 0.0;
        ++num_failed;
        continue;
      }
      auto inner_msg =
          (inner_res.message.HasField() && inner_res.message.Size() == k_in_)
              ? std::move(inner_res.message)
              : inner_code_->Unencode(inner_res.codeword);
      const base_storage_type *in_data = inner_msg.RawData().data();
      for (std::size_t s = 0; s < streams_; ++s) {
        auto chunk = SliceBase(in_data + s * rel_deg_, rel_deg_);
        outer_raws[s](j) = tower_.ContractElement(chunk);
      }
    }

    std::vector<GFArray<ExtField>> outer_received;
    outer_received.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      outer_received.push_back(GFArray<ExtField>::FromRaw(
          outer_code_->Field(), std::move(outer_raws[s])));
    }

    return RunGmdTrials(received_base, outer_received, failed, reliabilities,
                        num_failed);
  }

  // Hard-input Forney GMD decoder: computes inner reliability scores
  // theta_j = max(0.0, 1.0 - 2 * e_j / d_in) directly from each inner block's
  // decoded Hamming error count e_j (with failed inner blocks marked as
  // erased with theta_j = 0.0).
  result_type DecodeGmdHard(const codeword_type &received_base) const {
    ValidateReceivedWord(received_base, "ConcatenatedCode::DecodeGmdHard");

    const std::size_t d_in = inner_code_->MinimumDistance();
    using ExtRawArray = typename GFArray<ExtField>::raw_array_type;
    using ExtShapeType = typename GFArray<ExtField>::shape_type;
    const ext_storage_type ext_zero = outer_code_->Field()->AdditiveIdentity();

    std::vector<ExtRawArray> outer_raws;
    outer_raws.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      ExtRawArray arr = ExtRawArray::from_shape(ExtShapeType{n_out_});
      arr.fill(ext_zero);
      outer_raws.push_back(std::move(arr));
    }

    std::vector<bool> failed(n_out_, false);
    std::vector<double> reliabilities(n_out_, 0.0);
    std::size_t num_failed = 0;

    for (std::size_t j = 0; j < n_out_; ++j) {
      auto block = ExtractBlock(received_base, j);
      auto inner_res = inner_code_->DecodeWithResult(block);
      if (!inner_res.Ok()) {
        failed[j] = true;
        reliabilities[j] = 0.0;
        ++num_failed;
        continue;
      }
      const std::size_t e_j = inner_code_->Distance(block, inner_res.codeword);
      const double theta =
          (d_in > 0) ? std::max(0.0, 1.0 - (2.0 * static_cast<double>(e_j)) /
                                               static_cast<double>(d_in))
                     : 1.0;
      reliabilities[j] = theta;

      auto inner_msg =
          (inner_res.message.HasField() && inner_res.message.Size() == k_in_)
              ? std::move(inner_res.message)
              : inner_code_->Unencode(inner_res.codeword);
      const base_storage_type *in_data = inner_msg.RawData().data();
      for (std::size_t s = 0; s < streams_; ++s) {
        auto chunk = SliceBase(in_data + s * rel_deg_, rel_deg_);
        outer_raws[s](j) = tower_.ContractElement(chunk);
      }
    }

    std::vector<GFArray<ExtField>> outer_received;
    outer_received.reserve(streams_);
    for (std::size_t s = 0; s < streams_; ++s) {
      outer_received.push_back(GFArray<ExtField>::FromRaw(
          outer_code_->Field(), std::move(outer_raws[s])));
    }

    return RunGmdTrials(received_base, outer_received, failed, reliabilities,
                        num_failed);
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) +
           ", >= " + std::to_string(DesignedDistance()) +
           "] concatenated code over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 private:
  static std::shared_ptr<BaseField> ValidateAndGetBaseField(
      const std::shared_ptr<const AbstractLinearCode<ExtField>> &outer_code,
      const std::shared_ptr<const AbstractLinearCode<BaseField>> &inner_code) {
    if (!outer_code || !inner_code) {
      throw std::invalid_argument(
          "ConcatenatedCode: outer_code and inner_code cannot be null");
    }
    return inner_code->Field();
  }

  void RegisterEncodersAndDecoders() {
    this->RegisterEncoder(
        "Concatenated", [](const AbstractCode<BaseField> &code) {
          const auto &cc =
              dynamic_cast<const ConcatenatedCode<BaseField, ExtField> &>(code);
          return std::make_shared<ConcatenatedEncoder<BaseField, ExtField>>(cc);
        });
    this->RegisterEncoder(
        "GeneratorMatrix", [](const AbstractCode<BaseField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<BaseField> &>(code);
          return std::make_shared<GeneratorMatrixEncoder<BaseField>>(lin);
        });
    this->RegisterDecoder(
        "ConcatenatedStandard", [](const AbstractCode<BaseField> &code) {
          const auto &cc =
              dynamic_cast<const ConcatenatedCode<BaseField, ExtField> &>(code);
          return std::make_shared<
              ConcatenatedStandardDecoder<BaseField, ExtField>>(cc);
        });
    this->RegisterDecoder(
        "ConcatenatedGMD", [](const AbstractCode<BaseField> &code) {
          const auto &cc =
              dynamic_cast<const ConcatenatedCode<BaseField, ExtField> &>(code);
          return std::make_shared<ConcatenatedGmdDecoder<BaseField, ExtField>>(
              cc);
        });
  }

  void ValidateReceivedWord(const codeword_type &received_base,
                            const char *fn_name) const {
    if (!received_base.HasField() || received_base.Field() != this->field_) {
      throw std::invalid_argument(std::string(fn_name) + ": field mismatch");
    }
    if (received_base.Dimension() != 1 ||
        received_base.Size() != this->length_) {
      throw std::invalid_argument(
          std::string(fn_name) +
          ": received word must be 1D with length n_out * n_in");
    }
  }

  codeword_type ExtractBlock(const codeword_type &word,
                             std::size_t block_idx) const {
    return SliceBase(word.RawData().data() + block_idx * n_in_, n_in_);
  }

  codeword_type SliceBase(const base_storage_type *src, std::size_t len) const {
    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{len});
    if (len > 0) {
      std::copy(src, src + len, raw.data());
    }
    return codeword_type::FromRaw(this->field_, std::move(raw));
  }

  codeword_type EncodeOuterCodewords(
      const std::vector<GFArray<ExtField>> &outer_codewords) const {
    using BaseRawArray = typename codeword_type::raw_array_type;
    using BaseShapeType = typename codeword_type::shape_type;
    BaseRawArray out_raw =
        BaseRawArray::from_shape(BaseShapeType{this->length_});
    base_storage_type *out_data = out_raw.data();

    for (std::size_t j = 0; j < n_out_; ++j) {
      BaseRawArray in_msg_raw = BaseRawArray::from_shape(BaseShapeType{k_in_});
      base_storage_type *in_msg_data = in_msg_raw.data();
      for (std::size_t s = 0; s < streams_; ++s) {
        auto expanded = tower_.ExpandElement(outer_codewords[s].Raw(j));
        const base_storage_type *exp_src = expanded.RawData().data();
        std::copy(exp_src, exp_src + rel_deg_, in_msg_data + s * rel_deg_);
      }
      auto in_msg = codeword_type::FromRaw(this->field_, std::move(in_msg_raw));
      auto in_cw = inner_code_->Encode(in_msg);
      const base_storage_type *cw_src = in_cw.RawData().data();
      std::copy(cw_src, cw_src + n_in_, out_data + j * n_in_);
    }

    return codeword_type::FromRaw(this->field_, std::move(out_raw));
  }

  message_type ExpandOuterMessages(
      const std::vector<GFArray<ExtField>> &outer_msgs) const {
    using BaseRawArray = typename message_type::raw_array_type;
    using BaseShapeType = typename message_type::shape_type;
    BaseRawArray msg_raw =
        BaseRawArray::from_shape(BaseShapeType{this->dimension_});
    base_storage_type *msg_data = msg_raw.data();

    for (std::size_t s = 0; s < streams_; ++s) {
      for (std::size_t i = 0; i < k_out_; ++i) {
        auto expanded = tower_.ExpandElement(outer_msgs[s].Raw(i));
        const base_storage_type *src = expanded.RawData().data();
        std::copy(src, src + rel_deg_, msg_data + (s * k_out_ + i) * rel_deg_);
      }
    }
    return message_type::FromRaw(this->field_, std::move(msg_raw));
  }

  result_type RunGmdTrials(const codeword_type &received_base,
                           const std::vector<GFArray<ExtField>> &outer_received,
                           const std::vector<bool> &failed,
                           const std::vector<double> &reliabilities,
                           std::size_t num_failed) const {
    const std::size_t d_out = outer_code_->MinimumDistance();
    const std::size_t max_erasures = (d_out > 0) ? (d_out - 1) : 0;

    // Sort outer symbol positions in ascending order of reliability, placing
    // inner-decoding failures strictly first.
    std::vector<std::size_t> order(n_out_);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) {
                       if (failed[a] != failed[b]) {
                         return failed[a];
                       }
                       return reliabilities[a] < reliabilities[b];
                     });

    // Build GMD erasure-count schedule:
    //   odd  d_out: s = 0, 2, 4, ..., d_out - 1   (floor((d_out-1)/2) + 1)
    //   even d_out: s = 0 (if num_failed == 0), 1, 3, ..., d_out - 1
    // plus s = num_failed if num_failed <= max_erasures.
    std::vector<std::size_t> trial_counts;
    const std::size_t num_gmd_steps = (d_out > 0) ? ((d_out - 1) / 2 + 1) : 1;
    trial_counts.reserve(num_gmd_steps + 2);

    const std::size_t parity_offset = (d_out > 0 && (d_out % 2 == 0)) ? 1 : 0;
    if (parity_offset == 1 && num_failed == 0) {
      trial_counts.push_back(0);
    }
    for (std::size_t i = 0; i < num_gmd_steps; ++i) {
      const std::size_t s = 2 * i + parity_offset;
      if (s <= max_erasures && s >= num_failed && s <= n_out_) {
        trial_counts.push_back(s);
      }
    }
    if (num_failed <= max_erasures && num_failed <= n_out_) {
      trial_counts.push_back(num_failed);
    }
    std::sort(trial_counts.begin(), trial_counts.end());
    trial_counts.erase(std::unique(trial_counts.begin(), trial_counts.end()),
                       trial_counts.end());

    bool found = false;
    std::size_t best_dist = std::numeric_limits<std::size_t>::max();
    std::size_t best_erasures = 0;
    codeword_type best_cw;
    message_type best_msg;

    for (const std::size_t s_erase : trial_counts) {
      std::vector<std::size_t> erasures(
          order.begin(), order.begin() + static_cast<std::ptrdiff_t>(s_erase));
      std::sort(erasures.begin(), erasures.end());

      std::vector<GFArray<ExtField>> cand_outer_cws;
      std::vector<GFArray<ExtField>> cand_outer_msgs;
      cand_outer_cws.reserve(streams_);
      cand_outer_msgs.reserve(streams_);
      bool all_ok = true;

      for (std::size_t u = 0; u < streams_; ++u) {
        auto outer_res =
            outer_code_->DecodeWithErasures(outer_received[u], erasures);
        if (!outer_res.Ok()) {
          all_ok = false;
          break;
        }
        auto outer_msg =
            (outer_res.message.HasField() && outer_res.message.Size() == k_out_)
                ? std::move(outer_res.message)
                : outer_code_->Unencode(outer_res.codeword);
        cand_outer_cws.push_back(std::move(outer_res.codeword));
        cand_outer_msgs.push_back(std::move(outer_msg));
      }

      if (!all_ok) {
        continue;
      }

      auto cand_cw = EncodeOuterCodewords(cand_outer_cws);
      const std::size_t dist = this->Distance(received_base, cand_cw);
      if (!found || dist < best_dist) {
        found = true;
        best_dist = dist;
        best_erasures = s_erase;
        best_cw = std::move(cand_cw);
        best_msg = ExpandOuterMessages(cand_outer_msgs);
      }
    }

    if (!found) {
      return result_type::Failure(
          "ConcatenatedCode::DecodeGmd: no GMD erasure trial produced a valid "
          "outer codeword");
    }

    return result_type::Success(std::move(best_cw), std::move(best_msg),
                                best_dist, best_erasures);
  }

  std::shared_ptr<const AbstractLinearCode<ExtField>> outer_code_;
  std::shared_ptr<const AbstractLinearCode<BaseField>> inner_code_;
  FieldTowerAdapter<BaseField, ExtField> tower_;
  std::size_t n_out_;
  std::size_t k_out_;
  std::size_t n_in_;
  std::size_t k_in_;
  std::size_t rel_deg_;
  std::size_t streams_;

  mutable std::mutex matrix_mutex_;
  mutable std::optional<matrix_type> cached_generator_;
  mutable std::optional<matrix_type> cached_parity_check_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_CONCATENATED_HPP_
