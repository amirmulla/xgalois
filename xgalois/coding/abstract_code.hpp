#ifndef XGALOIS_CODING_ABSTRACT_CODE_HPP_
#define XGALOIS_CODING_ABSTRACT_CODE_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

enum class Metric : std::uint8_t { HAMMING, LEE, RANK };

template <typename GaloisField>
class AbstractCode {
 public:
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;
  using result_type = DecodingResult<GaloisField>;

  using EncoderFactory = std::function<std::shared_ptr<Encoder<GaloisField>>(
      const AbstractCode<GaloisField> &)>;
  using DecoderFactory = std::function<std::shared_ptr<Decoder<GaloisField>>(
      const AbstractCode<GaloisField> &)>;

  explicit AbstractCode(std::shared_ptr<GaloisField> field, std::size_t length,
                        Metric metric = Metric::HAMMING,
                        std::string default_encoder = "",
                        std::string default_decoder = "")
      : field_(std::move(field)),
        length_(length),
        metric_(metric),
        default_encoder_(std::move(default_encoder)),
        default_decoder_(std::move(default_decoder)) {
    if (!field_) {
      throw std::invalid_argument("AbstractCode: field cannot be null");
    }
    if (length_ == 0) {
      throw std::invalid_argument("AbstractCode: code length must be positive");
    }
  }

  // Rule of Five: caches hold raw pointers to `this`, so copy/move must never
  // transfer cached encoder/decoder instances pointing to the source object.
  AbstractCode(const AbstractCode &other)
      : field_(other.field_),
        length_(other.length_),
        metric_(other.metric_),
        default_encoder_(other.default_encoder_),
        default_decoder_(other.default_decoder_),
        registered_encoders_(other.registered_encoders_),
        registered_decoders_(other.registered_decoders_) {}

  AbstractCode &operator=(const AbstractCode &other) {
    if (this != &other) {
      field_ = other.field_;
      length_ = other.length_;
      metric_ = other.metric_;
      default_encoder_ = other.default_encoder_;
      default_decoder_ = other.default_decoder_;
      registered_encoders_ = other.registered_encoders_;
      registered_decoders_ = other.registered_decoders_;
      ClearCaches();
    }
    return *this;
  }

  AbstractCode(AbstractCode &&other) noexcept
      : field_(std::move(other.field_)),
        length_(other.length_),
        metric_(other.metric_),
        default_encoder_(std::move(other.default_encoder_)),
        default_decoder_(std::move(other.default_decoder_)),
        registered_encoders_(std::move(other.registered_encoders_)),
        registered_decoders_(std::move(other.registered_decoders_)) {
    other.ClearCaches();
  }

  AbstractCode &operator=(AbstractCode &&other) noexcept {
    if (this != &other) {
      field_ = std::move(other.field_);
      length_ = other.length_;
      metric_ = other.metric_;
      default_encoder_ = std::move(other.default_encoder_);
      default_decoder_ = std::move(other.default_decoder_);
      registered_encoders_ = std::move(other.registered_encoders_);
      registered_decoders_ = std::move(other.registered_decoders_);
      ClearCaches();
      other.ClearCaches();
    }
    return *this;
  }

  virtual ~AbstractCode() = default;

  std::shared_ptr<GaloisField> Field() const { return field_; }
  std::size_t Length() const noexcept { return length_; }
  Metric GetMetric() const noexcept { return metric_; }

  virtual std::size_t Dimension() const = 0;
  virtual std::size_t MinimumDistance() const = 0;
  virtual bool Contains(const codeword_type &word) const = 0;
  virtual codeword_type RandomCodeword() const = 0;
  virtual codeword_type RandomCodeword(std::mt19937_64 &rng) const = 0;

  double Rate() const {
    return static_cast<double>(Dimension()) / static_cast<double>(length_);
  }

  std::size_t Redundancy() const { return length_ - Dimension(); }

  double RelativeDistance() const {
    return static_cast<double>(MinimumDistance()) /
           static_cast<double>(length_);
  }

  codeword_type Encode(const message_type &message) const {
    return Encode(message, default_encoder_);
  }

  codeword_type Encode(const message_type &message,
                       const std::string &encoder_name) const {
    auto encoder = GetEncoder(encoder_name);
    return encoder->Encode(message);
  }

  message_type Unencode(const codeword_type &codeword) const {
    return Unencode(codeword, default_encoder_);
  }

  message_type Unencode(const codeword_type &codeword,
                        const std::string &encoder_name) const {
    auto encoder = GetEncoder(encoder_name);
    return encoder->Unencode(codeword);
  }

  result_type DecodeWithResult(const codeword_type &received_word) const {
    return DecodeWithResult(received_word, default_decoder_);
  }

  result_type DecodeWithResult(const codeword_type &received_word,
                               const std::string &decoder_name) const {
    auto decoder = GetDecoder(decoder_name);
    return decoder->DecodeWithResult(received_word);
  }

  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const {
    return DecodeWithErasures(received_word, erasure_positions,
                              default_decoder_);
  }

  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions,
      const std::string &decoder_name) const {
    auto decoder = GetDecoder(decoder_name);
    return decoder->DecodeWithErasures(received_word, erasure_positions);
  }

  codeword_type DecodeToCode(const codeword_type &received_word) const {
    return DecodeToCode(received_word, default_decoder_);
  }

  codeword_type DecodeToCode(const codeword_type &received_word,
                             const std::string &decoder_name) const {
    auto decoder = GetDecoder(decoder_name);
    return decoder->DecodeToCode(received_word);
  }

  message_type DecodeToMessage(const codeword_type &received_word) const {
    return DecodeToMessage(received_word, default_decoder_);
  }

  message_type DecodeToMessage(const codeword_type &received_word,
                               const std::string &decoder_name) const {
    auto decoder = GetDecoder(decoder_name);
    return decoder->DecodeToMessage(received_word);
  }

  std::shared_ptr<const Encoder<GaloisField>> GetEncoder(
      const std::string &encoder_name = "") const {
    const std::string name =
        encoder_name.empty() ? default_encoder_ : encoder_name;
    std::lock_guard<std::mutex> lock(cache_mutex_);
    auto cache_it = encoder_cache_.find(name);
    if (cache_it != encoder_cache_.end()) {
      return cache_it->second;
    }
    auto reg_it = registered_encoders_.find(name);
    if (reg_it == registered_encoders_.end()) {
      throw std::invalid_argument("Unknown encoder: " + name);
    }
    auto enc = reg_it->second(*this);
    encoder_cache_.emplace(name, enc);
    return enc;
  }

  std::shared_ptr<const Decoder<GaloisField>> GetDecoder(
      const std::string &decoder_name = "") const {
    const std::string name =
        decoder_name.empty() ? default_decoder_ : decoder_name;
    std::lock_guard<std::mutex> lock(cache_mutex_);
    auto cache_it = decoder_cache_.find(name);
    if (cache_it != decoder_cache_.end()) {
      return cache_it->second;
    }
    auto reg_it = registered_decoders_.find(name);
    if (reg_it == registered_decoders_.end()) {
      throw std::invalid_argument("Unknown decoder: " + name);
    }
    auto dec = reg_it->second(*this);
    decoder_cache_.emplace(name, dec);
    return dec;
  }

  void RegisterEncoder(const std::string &name, EncoderFactory factory) {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    registered_encoders_[name] = std::move(factory);
    encoder_cache_.erase(name);
  }

  void RegisterDecoder(const std::string &name, DecoderFactory factory) {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    registered_decoders_[name] = std::move(factory);
    decoder_cache_.erase(name);
  }

  std::vector<std::string> AvailableEncoders() const {
    std::vector<std::string> names;
    names.reserve(registered_encoders_.size());
    for (const auto &[k, _] : registered_encoders_) {
      names.push_back(k);
    }
    std::sort(names.begin(), names.end());
    return names;
  }

  std::vector<std::string> AvailableDecoders() const {
    std::vector<std::string> names;
    names.reserve(registered_decoders_.size());
    for (const auto &[k, _] : registered_decoders_) {
      names.push_back(k);
    }
    std::sort(names.begin(), names.end());
    return names;
  }

  const std::string &DefaultEncoder() const noexcept {
    return default_encoder_;
  }
  const std::string &DefaultDecoder() const noexcept {
    return default_decoder_;
  }

  void SetDefaultEncoder(const std::string &name) {
    if (registered_encoders_.find(name) == registered_encoders_.end()) {
      throw std::invalid_argument("Unknown encoder: " + name);
    }
    default_encoder_ = name;
  }

  void SetDefaultDecoder(const std::string &name) {
    if (registered_decoders_.find(name) == registered_decoders_.end()) {
      throw std::invalid_argument("Unknown decoder: " + name);
    }
    default_decoder_ = name;
  }

  std::size_t Distance(const codeword_type &w1, const codeword_type &w2) const {
    if (!w1.HasField() || !w2.HasField() || w1.Field() != field_ ||
        w2.Field() != field_) {
      throw std::invalid_argument("Distance: field mismatch");
    }
    if (w1.Size() != length_ || w2.Size() != length_) {
      throw std::invalid_argument("Distance: words must have code length");
    }
    return Weight(w1 - w2);
  }

  std::size_t Weight(const codeword_type &word) const {
    if (!word.HasField() || word.Field() != field_) {
      throw std::invalid_argument("Weight: field mismatch");
    }
    if (word.Size() != length_) {
      throw std::invalid_argument("Weight: word must have code length");
    }
    const auto span = word.RawSpan();
    if (metric_ == Metric::HAMMING) {
      std::size_t w = 0;
      for (const auto &v : span) {
        if (!field_->IsZero(v)) {
          ++w;
        }
      }
      return w;
    }
    if (metric_ == Metric::LEE) {
      if (field_->Degree() != 1) {
        throw std::invalid_argument(
            "Lee metric is defined only for prime fields");
      }
      if constexpr (std::is_integral_v<storage_type>) {
        const uint64_t p = field_->Characteristic();
        std::size_t total = 0;
        for (const auto &v : span) {
          const uint64_t val =
              static_cast<uint64_t>(field_->GetElementValue(v));
          total += static_cast<std::size_t>(std::min(val, p - val));
        }
        return total;
      } else {
        throw std::invalid_argument(
            "Lee metric is defined only for prime fields");
      }
    }
    throw std::invalid_argument(
        "Rank metric requires matrix-valued codewords; use "
        "linalg::matrix_rank");
  }

  virtual std::string ToString() const = 0;

 protected:
  void ClearCaches() const {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    encoder_cache_.clear();
    decoder_cache_.clear();
  }

  std::shared_ptr<GaloisField> field_;
  std::size_t length_;
  Metric metric_;
  std::string default_encoder_;
  std::string default_decoder_;

  std::unordered_map<std::string, EncoderFactory> registered_encoders_;
  std::unordered_map<std::string, DecoderFactory> registered_decoders_;

  mutable std::mutex cache_mutex_;
  mutable std::unordered_map<std::string,
                             std::shared_ptr<const Encoder<GaloisField>>>
      encoder_cache_;
  mutable std::unordered_map<std::string,
                             std::shared_ptr<const Decoder<GaloisField>>>
      decoder_cache_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ABSTRACT_CODE_HPP_
