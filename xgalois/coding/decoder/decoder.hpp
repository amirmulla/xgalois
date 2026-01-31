#ifndef XGALOIS_CODING_DECODER_HPP
#define XGALOIS_CODING_DECODER_HPP

#include <vector>
#include <memory>
#include <stdexcept>
#include <xtensor/containers/xarray.hpp>

namespace xg {
namespace coding {

// Forward declaration
template <typename ElementType> class AbstractCode;

template <typename ElementType>
class Decoder {
public:
    using element_type = ElementType;
    using codeword_type = xt::xarray<ElementType>;
    using message_type = xt::xarray<ElementType>;

    // Constructor
    explicit Decoder(const AbstractCode<ElementType>* code) : code_(code) {}

    virtual ~Decoder() = default;

    // Main decoding functions
    virtual codeword_type DecodeToCode(const codeword_type& received_word) const = 0;
    virtual message_type DecodeToMessage(const codeword_type& received_word) const = 0;

    // Get the input space dimension (usually same as code length)
    virtual size_t InputLength() const { return code_->Length(); }

    // Get the code this decoder is associated with
    const AbstractCode<ElementType>* GetCode() const { return code_; }

    // String representation
    virtual std::string ToString() const = 0;

protected:
    const AbstractCode<ElementType>* code_;
};

} // namespace coding
} // namespace xg

#endif // XGALOIS_CODING_DECODER_HPP
