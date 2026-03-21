#ifndef XGALOIS_CODING_ENCODER_HPP
#define XGALOIS_CODING_ENCODER_HPP

#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

// Forward declaration
template <typename GaloisField>
class AbstractCode;

template <typename GaloisField>
class Encoder {
 public:
  using element_type = xg::GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<GaloisField>;
  using message_type = xt::xarray<GaloisField>;

  // Constructor
  explicit Encoder(const AbstractCode<GaloisField>* code) : code_(code) {}

  virtual ~Encoder() = default;

  // Main encoding function
  virtual codeword_type Encode(const message_type& message) const = 0;

  // Unencode function (inverse of encode)
  virtual message_type Unencode(const codeword_type& codeword) const = 0;

  // Get the message space dimension
  virtual size_t MessageLength() const = 0;

  // Get the code this encoder is associated with
  const AbstractCode<GaloisField>* GetCode() const { return code_; }

  // String representation
  virtual std::string ToString() const = 0;

 protected:
  const AbstractCode<GaloisField>* code_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ENCODER_HPP
