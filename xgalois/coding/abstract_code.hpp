#ifndef XGALOIS_CODING_ABSTRACT_CODE_HPP
#define XGALOIS_CODING_ABSTRACT_CODE_HPP

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "xgalois/field/gf_base.hpp"

namespace xg {
namespace coding {

// Forward declarations
template <typename ElementType>
class Encoder;
template <typename ElementType>
class Decoder;

enum class Metric { HAMMING, RANK, LEE };

template <typename ElementType>
class AbstractCode {
 public:
  using element_type = ElementType;
  using codeword_type = xt::xarray<ElementType>;
  using message_type = xt::xarray<ElementType>;
  using encoder_factory_type =
      std::function<std::unique_ptr<Encoder<ElementType>>(
          const AbstractCode<ElementType>*)>;
  using decoder_factory_type =
      std::function<std::unique_ptr<Decoder<ElementType>>(
          const AbstractCode<ElementType>*)>;

  // Constructor
  AbstractCode(size_t length, const std::string& default_encoder_name = "",
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
  virtual std::shared_ptr<GaloisFieldBase<ElementType>> Field() const = 0;

  // Encoding/Decoding interface
  virtual codeword_type Encode(const message_type& message,
                               const std::string& encoder_name = "") const;
  virtual message_type DecodeToMessage(
      const codeword_type& received_word,
      const std::string& decoder_name = "") const;
  virtual codeword_type DecodeToCode(
      const codeword_type& received_word,
      const std::string& decoder_name = "") const;
  virtual message_type Unencode(const codeword_type& codeword,
                                const std::string& encoder_name = "") const;

  // Encoder/Decoder management
  void RegisterEncoder(const std::string& name, encoder_factory_type factory);
  void RegisterDecoder(const std::string& name, decoder_factory_type factory);

  std::unique_ptr<Encoder<ElementType>> GetEncoder(
      const std::string& name = "") const;
  std::unique_ptr<Decoder<ElementType>> GetDecoder(
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
  mutable std::unordered_map<std::string, std::unique_ptr<Encoder<ElementType>>>
      encoder_cache_;
  mutable std::unordered_map<std::string, std::unique_ptr<Decoder<ElementType>>>
      decoder_cache_;
};

// Implementation of template methods

template <typename ElementType>
typename AbstractCode<ElementType>::codeword_type
AbstractCode<ElementType>::Encode(const message_type& message,
                                  const std::string& encoder_name) const {
  auto encoder = GetEncoder(encoder_name);
  if (!encoder) {
    throw std::runtime_error("No encoder available for this code");
  }
  return encoder->Encode(message);
}

template <typename ElementType>
typename AbstractCode<ElementType>::message_type
AbstractCode<ElementType>::DecodeToMessage(
    const codeword_type& received_word, const std::string& decoder_name) const {
  auto decoder = GetDecoder(decoder_name);
  if (!decoder) {
    throw std::runtime_error("No decoder available for this code");
  }
  return decoder->DecodeToMessage(received_word);
}

template <typename ElementType>
typename AbstractCode<ElementType>::codeword_type
AbstractCode<ElementType>::DecodeToCode(const codeword_type& received_word,
                                        const std::string& decoder_name) const {
  auto decoder = GetDecoder(decoder_name);
  if (!decoder) {
    throw std::runtime_error("No decoder available for this code");
  }
  return decoder->DecodeToCode(received_word);
}

template <typename ElementType>
typename AbstractCode<ElementType>::message_type
AbstractCode<ElementType>::Unencode(const codeword_type& codeword,
                                    const std::string& encoder_name) const {
  auto encoder = GetEncoder(encoder_name);
  if (!encoder) {
    throw std::runtime_error("No encoder available for this code");
  }
  return encoder->Unencode(codeword);
}

template <typename ElementType>
void AbstractCode<ElementType>::RegisterEncoder(const std::string& name,
                                                encoder_factory_type factory) {
  encoder_registry_[name] = factory;
}

template <typename ElementType>
void AbstractCode<ElementType>::RegisterDecoder(const std::string& name,
                                                decoder_factory_type factory) {
  decoder_registry_[name] = factory;
}

template <typename ElementType>
std::unique_ptr<Encoder<ElementType>> AbstractCode<ElementType>::GetEncoder(
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

template <typename ElementType>
std::unique_ptr<Decoder<ElementType>> AbstractCode<ElementType>::GetDecoder(
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

template <typename ElementType>
std::vector<std::string> AbstractCode<ElementType>::AvailableEncoders() const {
  std::vector<std::string> names;
  for (const auto& pair : encoder_registry_) {
    names.push_back(pair.first);
  }
  return names;
}

template <typename ElementType>
std::vector<std::string> AbstractCode<ElementType>::AvailableDecoders() const {
  std::vector<std::string> names;
  for (const auto& pair : decoder_registry_) {
    names.push_back(pair.first);
  }
  return names;
}

template <typename ElementType>
size_t AbstractCode<ElementType>::Distance(const codeword_type& a,
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

template <typename ElementType>
size_t AbstractCode<ElementType>::Weight(const codeword_type& word) const {
  size_t weight = 0;
  auto field = Field();

  switch (metric_) {
    case Metric::HAMMING:
      for (size_t i = 0; i < word.size(); ++i) {
        if (word(i) !=
            ElementType{}) {  // Assuming zero element is default constructed
          weight++;
        }
      }
      break;
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
