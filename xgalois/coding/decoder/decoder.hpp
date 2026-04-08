#ifndef XGALOIS_CODING_DECODER_HPP
#define XGALOIS_CODING_DECODER_HPP

#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class AbstractCode;

template <typename GaloisField>
class Decoder {
 public:
  using element_type = GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<element_type>;
  using message_type = xt::xarray<element_type>;

  explicit Decoder(const AbstractCode<GaloisField>* code) : code_(code) {}

  virtual ~Decoder() = default;

  virtual codeword_type DecodeToCode(
      const codeword_type& received_word) const = 0;
  virtual message_type DecodeToMessage(
      const codeword_type& received_word) const = 0;

  virtual size_t InputLength() const { return code_->Length(); }

  const AbstractCode<GaloisField>* GetCode() const { return code_; }

  virtual std::string ToString() const = 0;

 protected:
  const AbstractCode<GaloisField>* code_;
};

}  // namespace coding
}  // namespace xg

#endif
