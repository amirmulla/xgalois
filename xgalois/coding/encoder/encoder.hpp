#ifndef XGALOIS_CODING_ENCODER_HPP
#define XGALOIS_CODING_ENCODER_HPP

#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class AbstractCode;

template <typename GaloisField>
class Encoder {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<element_type>;
  using message_type = xt::xarray<element_type>;

  explicit Encoder(const AbstractCode<GaloisField>* code) : code_(code) {}

  virtual ~Encoder() = default;

  virtual codeword_type Encode(const message_type& message) const = 0;

  virtual message_type Unencode(const codeword_type& codeword) const = 0;

  virtual size_t MessageLength() const = 0;

  const AbstractCode<GaloisField>* GetCode() const { return code_; }

  virtual std::string ToString() const = 0;

 protected:
  const AbstractCode<GaloisField>* code_;
};

}
}

#endif
