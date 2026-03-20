#ifndef XGALOIS_CODING_ENCODER_HPP
#define XGALOIS_CODING_ENCODER_HPP

#include <memory>
#include <stdexcept>
#include <vector>
#include <xtensor/containers/xarray.hpp>

namespace xg {
namespace coding {

// Forward declaration
template <typename ElementType>
class AbstractCode;

template <typename ElementType>
class Encoder {
 public:
  using element_type = ElementType;
  using codeword_type = xt::xarray<ElementType>;
  using message_type = xt::xarray<ElementType>;

  // Constructor
  explicit Encoder(const AbstractCode<ElementType>* code) : code_(code) {}

  virtual ~Encoder() = default;

  // Main encoding function
  virtual codeword_type Encode(const message_type& message) const = 0;

  // Unencode function (inverse of encode)
  virtual message_type Unencode(const codeword_type& codeword) const = 0;

  // Get the message space dimension
  virtual size_t MessageLength() const = 0;

  // Get the code this encoder is associated with
  const AbstractCode<ElementType>* GetCode() const { return code_; }

  // String representation
  virtual std::string ToString() const = 0;

 protected:
  const AbstractCode<ElementType>* code_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ENCODER_HPP
