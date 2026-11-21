#ifndef XGALOIS_CODING_DECODING_RESULT_HPP_
#define XGALOIS_CODING_DECODING_RESULT_HPP_

#include <cstddef>
#include <string>
#include <utility>

#include "xgalois/field/gf_array.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
struct DecodingResult {
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;

  codeword_type codeword;
  message_type message;
  std::size_t num_errors = 0;
  std::size_t num_erasures = 0;
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }

  static DecodingResult Success(codeword_type cw, message_type msg,
                                std::size_t errors = 0,
                                std::size_t erasures = 0) {
    DecodingResult res;
    res.codeword = std::move(cw);
    res.message = std::move(msg);
    res.num_errors = errors;
    res.num_erasures = erasures;
    res.success = true;
    return res;
  }

  static DecodingResult SuccessCodewordOnly(codeword_type cw,
                                            std::size_t errors = 0,
                                            std::size_t erasures = 0) {
    DecodingResult res;
    res.codeword = std::move(cw);
    res.num_errors = errors;
    res.num_erasures = erasures;
    res.success = true;
    return res;
  }

  static DecodingResult Failure(std::string reason) {
    DecodingResult res;
    res.success = false;
    res.error_message = std::move(reason);
    return res;
  }
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODING_RESULT_HPP_
