#ifndef XGALOIS_CODING_ABSTRACT_CODE_HPP
#define XGALOIS_CODING_ABSTRACT_CODE_HPP

#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_element.hpp"

namespace xg {
namespace coding {

// Forward declarations
template <typename GaloisField>
class Encoder;
template <typename GaloisField>
class Decoder;

enum class Metric : std::uint8_t { HAMMING, RANK, LEE };

template <typename GaloisField>
class AbstractCode {
 public:
  using element_type = GaloisFieldElement<GaloisField>;
  using codeword_type = xt::xarray<element_type>;
  using message_type = xt::xarray<element_type>;
  using encoder_factory_type =
      std::function<std::unique_ptr<Encoder<GaloisField>>(
          const AbstractCode<GaloisField>*)>;
  using decoder_factory_type =
      std::function<std::unique_ptr<Decoder<GaloisField>>(
          const AbstractCode<GaloisField>*)>;

  // Constructor
  explicit AbstractCode(size_t length,
                        const std::string& default_encoder_name = "",
                        const std::string& default_decoder_name = "",
                        Metric metric = Metric::HAMMING)
      : length_(length),
        default_encoder_name_(default_encoder_name),
        default_decoder_name_(default_decoder_name),
        metric_(metric) {}

  virtual ~AbstractCode() = default;

  // Basic properties
  size_t Length() const { return length_; }
  virtual size_t Dimension() const = 0;
  virtual size_t MinimumDistance() const = 0;
  Metric GetMetric() const { return metric_; }

  // Field operations
  virtual std::shared_ptr<GaloisField> Field() const = 0;

  // Encoding/Decoding interface
  virtual codeword_type Encode(const message_type& message,
                               const std::string& encoder_name) const;
  codeword_type Encode(const message_type& message) const {
    return Encode(message, "");
  }

  virtual message_type DecodeToMessage(const codeword_type& received_word,
                                       const std::string& decoder_name) const;
  message_type DecodeToMessage(const codeword_type& received_word) const {
    return DecodeToMessage(received_word, "");
  }

  virtual codeword_type DecodeToCode(const codeword_type& received_word,
                                     const std::string& decoder_name) const;
  codeword_type DecodeToCode(const codeword_type& received_word) const {
    return DecodeToCode(received_word, "");
  }

  virtual message_type Unencode(const codeword_type& codeword,
                                const std::string& encoder_name) const;
  message_type Unencode(const codeword_type& codeword) const {
    return Unencode(codeword, "");
  }

  // Encoder/Decoder management
  void RegisterEncoder(const std::string& name, encoder_factory_type factory);
  void RegisterDecoder(const std::string& name, decoder_factory_type factory);

  std::unique_ptr<Encoder<GaloisField>> GetEncoder(
      const std::string& name = "") const;
  std::unique_ptr<Decoder<GaloisField>> GetDecoder(
      const std::string& name = "") const;

  std::vector<std::string> AvailableEncoders() const;
  std::vector<std::string> AvailableDecoders() const;

  // Code properties
  virtual bool Contains(const codeword_type& word) const = 0;
  virtual std::vector<codeword_type> GetCodewords() const = 0;
  virtual codeword_type RandomCodeword() const = 0;

  // Distance calculations
  virtual size_t Distance(const codeword_type& a, const codeword_type& b) const;
  virtual size_t Weight(const codeword_type& word) const;

  // String representation
  virtual std::string ToString() const = 0;

  // Convenience operator for encoding
  codeword_type operator()(const message_type& message) const {
    return Encode(message);
  }

 protected:
  size_t length_;
  std::string default_encoder_name_;
  std::string default_decoder_name_;
  Metric metric_;

  // Encoder/Decoder registries
  std::unordered_map<std::string, encoder_factory_type> encoder_registry_;
  std::unordered_map<std::string, decoder_factory_type> decoder_registry_;

  // Cached encoders/decoders
  mutable std::unordered_map<std::string, std::unique_ptr<Encoder<GaloisField>>>
      encoder_cache_;
  mutable std::unordered_map<std::string, std::unique_ptr<Decoder<GaloisField>>>
      decoder_cache_;
};

// Implementation of template methods

template <typename GaloisField>
typename AbstractCode<GaloisField>::codeword_type
AbstractCode<GaloisField>::Encode(const message_type& message,
                                  const std::string& encoder_name) const {
  auto encoder = GetEncoder(encoder_name);
  if (!encoder) {
    throw std::runtime_error("No encoder available for this code");
  }
  return encoder->Encode(message);
}

template <typename GaloisField>
typename AbstractCode<GaloisField>::message_type
AbstractCode<GaloisField>::DecodeToMessage(
    const codeword_type& received_word, const std::string& decoder_name) const {
  auto decoder = GetDecoder(decoder_name);
  if (!decoder) {
    throw std::runtime_error("No decoder available for this code");
  }
  return decoder->DecodeToMessage(received_word);
}

template <typename GaloisField>
typename AbstractCode<GaloisField>::codeword_type
AbstractCode<GaloisField>::DecodeToCode(const codeword_type& received_word,
                                        const std::string& decoder_name) const {
  auto decoder = GetDecoder(decoder_name);
  if (!decoder) {
    throw std::runtime_error("No decoder available for this code");
  }
  return decoder->DecodeToCode(received_word);
}

template <typename GaloisField>
typename AbstractCode<GaloisField>::message_type
AbstractCode<GaloisField>::Unencode(const codeword_type& codeword,
                                    const std::string& encoder_name) const {
  auto encoder = GetEncoder(encoder_name);
  if (!encoder) {
    throw std::runtime_error("No encoder available for this code");
  }
  return encoder->Unencode(codeword);
}

template <typename GaloisField>
void AbstractCode<GaloisField>::RegisterEncoder(const std::string& name,
                                                encoder_factory_type factory) {
  encoder_registry_[name] = factory;
}

template <typename GaloisField>
void AbstractCode<GaloisField>::RegisterDecoder(const std::string& name,
                                                decoder_factory_type factory) {
  decoder_registry_[name] = factory;
}

template <typename GaloisField>
std::unique_ptr<Encoder<GaloisField>> AbstractCode<GaloisField>::GetEncoder(
    const std::string& name) const {
  std::string encoder_name = name.empty() ? default_encoder_name_ : name;

  if (encoder_name.empty()) {
    throw std::runtime_error(
        "No encoder name specified and no default encoder set");
  }

  // Check cache first
  auto cache_it = encoder_cache_.find(encoder_name);
  if (cache_it != encoder_cache_.end()) {
    // Return a copy since we can't return the same unique_ptr twice
    auto it = encoder_registry_.find(encoder_name);
    if (it != encoder_registry_.end()) {
      return it->second(this);
    }
  }

  // Create new encoder
  auto it = encoder_registry_.find(encoder_name);
  if (it == encoder_registry_.end()) {
    throw std::runtime_error("Unknown encoder: " + encoder_name);
  }

  return it->second(this);
}

template <typename GaloisField>
std::unique_ptr<Decoder<GaloisField>> AbstractCode<GaloisField>::GetDecoder(
    const std::string& name) const {
  std::string decoder_name = name.empty() ? default_decoder_name_ : name;

  if (decoder_name.empty()) {
    throw std::runtime_error(
        "No decoder name specified and no default decoder set");
  }

  // Check cache first
  auto cache_it = decoder_cache_.find(decoder_name);
  if (cache_it != decoder_cache_.end()) {
    // Return a copy since we can't return the same unique_ptr twice
    auto it = decoder_registry_.find(decoder_name);
    if (it != decoder_registry_.end()) {
      return it->second(this);
    }
  }

  // Create new decoder
  auto it = decoder_registry_.find(decoder_name);
  if (it == decoder_registry_.end()) {
    throw std::runtime_error("Unknown decoder: " + decoder_name);
  }

  return it->second(this);
}

template <typename GaloisField>
std::vector<std::string> AbstractCode<GaloisField>::AvailableEncoders() const {
  std::vector<std::string> names;
  for (const auto& pair : encoder_registry_) {
    names.push_back(pair.first);
  }
  return names;
}

template <typename GaloisField>
std::vector<std::string> AbstractCode<GaloisField>::AvailableDecoders() const {
  std::vector<std::string> names;
  for (const auto& pair : decoder_registry_) {
    names.push_back(pair.first);
  }
  return names;
}

template <typename GaloisField>
size_t AbstractCode<GaloisField>::Distance(const codeword_type& a,
                                           const codeword_type& b) const {
  if (a.size() != b.size()) {
    throw std::invalid_argument("Codewords must have the same length");
  }

  size_t distance = 0;
  auto field = Field();

  switch (metric_) {
    case Metric::HAMMING:
      for (size_t i = 0; i < a.size(); ++i) {
        if (a(i) != b(i)) {
          distance++;
        }
      }
      break;
    case Metric::LEE:
      // Lee distance implementation would depend on field characteristics
      throw std::runtime_error("Lee distance not implemented yet");
    case Metric::RANK:
      // Rank distance implementation would be more complex
      throw std::runtime_error("Rank distance not implemented yet");
  }

  return distance;
}

template <typename GaloisField>
size_t AbstractCode<GaloisField>::Weight(const codeword_type& word) const {
  size_t weight = 0;
  auto field = Field();

  switch (metric_) {
    case Metric::HAMMING: {
      auto zero_element = field->AdditiveIdentity();
      for (size_t i = 0; i < word.size(); ++i) {
        if (word(i) != zero_element) {
          weight++;
        }
      }
      break;
    }
    case Metric::LEE:
      throw std::runtime_error("Lee weight not implemented yet");
    case Metric::RANK:
      throw std::runtime_error("Rank weight not implemented yet");
  }

  return weight;
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ABSTRACT_CODE_HPP
