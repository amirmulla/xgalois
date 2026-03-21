#ifndef XGALOIS_CODING_DECODER_HPP
#define XGALOIS_CODING_DECODER_HPP

#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

// Forward declaration
template <typename GaloisField>
class AbstractCode;

template <typename GaloisField>
class Decoder {
 public:
  using element_type = GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<element_type>;
  using message_type = xt::xarray<element_type>;

  // Constructor
  explicit Decoder(const AbstractCode<GaloisField>* code) : code_(code) {}

  virtual ~Decoder() = default;

  // Main decoding functions
  virtual codeword_type DecodeToCode(
      const codeword_type& received_word) const = 0;
  virtual message_type DecodeToMessage(
      const codeword_type& received_word) const = 0;

  // Get the input space dimension (usually same as code length)
  virtual size_t InputLength() const { return code_->Length(); }

  // Get the code this decoder is associated with
  const AbstractCode<GaloisField>* GetCode() const { return code_; }

  // String representation
  virtual std::string ToString() const = 0;

 protected:
  const AbstractCode<GaloisField>* code_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_HPP
