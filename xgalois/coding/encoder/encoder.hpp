#ifndef XGALOIS_CODING_ENCODER_ENCODER_HPP_
#define XGALOIS_CODING_ENCODER_ENCODER_HPP_

#include <cstddef>
#include <stdexcept>
#include <string>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class AbstractCode;

template <typename GaloisField>
class Encoder {
 public:
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;

  explicit Encoder(const AbstractCode<GaloisField> &code) : code_(&code) {}

  virtual ~Encoder() = default;

  virtual codeword_type Encode(const message_type &message) const = 0;

  virtual message_type Unencode(const codeword_type &codeword) const = 0;

  virtual std::string Name() const = 0;

  const AbstractCode<GaloisField> &Code() const { return *code_; }

  virtual std::string ToString() const {
    return Name() + " encoder for " + code_->ToString();
  }

 protected:
  const AbstractCode<GaloisField> *code_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ENCODER_ENCODER_HPP_
