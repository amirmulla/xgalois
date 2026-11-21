#ifndef XGALOIS_CODING_DECODER_DECODER_HPP_
#define XGALOIS_CODING_DECODER_DECODER_HPP_

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class AbstractCode;

template <typename GaloisField>
class Decoder {
 public:
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;
  using result_type = DecodingResult<GaloisField>;

  explicit Decoder(const AbstractCode<GaloisField> &code) : code_(&code) {}

  virtual ~Decoder() = default;

  // Primary fallible decoding API: never silently returns corrupted words.
  virtual result_type DecodeWithResult(
      const codeword_type &received_word) const = 0;

  // Errors-and-erasures decoding (overridden by RS/GRS/BCH decoders).
  virtual result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const {
    if (erasure_positions.empty()) {
      return DecodeWithResult(received_word);
    }
    return result_type::Failure(Name() + " does not support erasure decoding");
  }

  virtual codeword_type DecodeToCode(const codeword_type &received_word) const {
    auto res = DecodeWithResult(received_word);
    if (!res.Ok()) {
      throw std::runtime_error("Decoding failed: " + res.error_message);
    }
    return std::move(res.codeword);
  }

  virtual message_type DecodeToMessage(
      const codeword_type &received_word) const {
    auto res = DecodeWithResult(received_word);
    if (!res.Ok()) {
      throw std::runtime_error("Decoding failed: " + res.error_message);
    }
    if (res.message.HasField() && !res.message.Empty()) {
      return std::move(res.message);
    }
    return code_->Unencode(res.codeword);
  }

  virtual std::string Name() const = 0;

  virtual std::size_t DecodingRadius() const = 0;

  const AbstractCode<GaloisField> &Code() const { return *code_; }

  virtual std::string ToString() const {
    return Name() + " decoder for " + code_->ToString();
  }

 protected:
  const AbstractCode<GaloisField> *code_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_DECODER_HPP_
