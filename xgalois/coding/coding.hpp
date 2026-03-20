#ifndef XGALOIS_CODING_HPP
#define XGALOIS_CODING_HPP

// Base classes
#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/coding/abstract_linear_code.hpp"

// Encoder/Decoder framework
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/encoder/encoder.hpp"

// Specific encoders
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/encoder/grs_encoder.hpp"

// Specific decoders
#include "xgalois/coding/decoder/grs_decoder.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"

// Code implementations
#include "xgalois/coding/cyclic_code.hpp"
#include "xgalois/coding/grs.hpp"

namespace xg {
namespace coding {

// Convenience functions for creating codes

// Create a standard Reed-Solomon code
template <typename ElementType>
std::unique_ptr<GeneralizedReedSolomonCode<ElementType>> CreateReedSolomonCode(
    std::shared_ptr<GaloisFieldBase<ElementType>> field, size_t length,
    size_t dimension, const ElementType& primitive_element) {
  return GeneralizedReedSolomonCode<ElementType>::ReedSolomon(
      field, length, dimension, primitive_element);
}

// Create a standard Reed-Solomon code with consecutive evaluation points
template <typename ElementType>
std::unique_ptr<GeneralizedReedSolomonCode<ElementType>> CreateReedSolomonCode(
    std::shared_ptr<GaloisFieldBase<ElementType>> field, size_t length,
    size_t dimension) {
  // Use primitive element (assuming it's 2 for binary extension fields)
  ElementType primitive_element(2);
  return CreateReedSolomonCode(field, length, dimension, primitive_element);
}

// Create a generalized Reed-Solomon code
template <typename ElementType>
std::unique_ptr<GeneralizedReedSolomonCode<ElementType>>
CreateGeneralizedReedSolomonCode(
    std::shared_ptr<GaloisFieldBase<ElementType>> field,
    const xt::xarray<ElementType>& evaluation_points,
    const xt::xarray<ElementType>& column_multipliers, size_t dimension) {
  return std::make_unique<GeneralizedReedSolomonCode<ElementType>>(
      field, evaluation_points, column_multipliers, dimension);
}

// Create a cyclic code from generator polynomial
template <typename ElementType>
std::unique_ptr<CyclicCode<ElementType>> CreateCyclicCode(
    std::shared_ptr<GaloisFieldBase<ElementType>> field, size_t length,
    const xg::poly::PolyDense<ElementType>& generator_poly) {
  return std::make_unique<CyclicCode<ElementType>>(field, length,
                                                   generator_poly);
}

// Error correction capability helper
template <typename ElementType>
size_t ErrorCorrectionCapability(const AbstractCode<ElementType>& code) {
  return (code.MinimumDistance() - 1) / 2;
}

// Code rate helper
template <typename ElementType>
double CodeRate(const AbstractLinearCode<ElementType>& code) {
  return code.Rate();
}

// Redundancy helper
template <typename ElementType>
size_t Redundancy(const AbstractLinearCode<ElementType>& code) {
  return code.Redundancy();
}

// Check if code is MDS (Maximum Distance Separable)
template <typename ElementType>
bool IsMDS(const AbstractLinearCode<ElementType>& code) {
  return code.MinimumDistance() == code.Length() - code.Dimension() + 1;
}

// Singleton bound check
template <typename ElementType>
bool SatisfiesSingletonBound(const AbstractLinearCode<ElementType>& code) {
  return code.MinimumDistance() <= code.Length() - code.Dimension() + 1;
}

// Hamming bound check
template <typename ElementType>
bool SatisfiesHammingBound(const AbstractLinearCode<ElementType>& code) {
  auto field = code.Field();
  size_t q = field->Order();
  size_t n = code.Length();
  size_t k = code.Dimension();
  size_t t = ErrorCorrectionCapability(code);

  // Calculate volume of Hamming sphere
  size_t volume = 1;
  for (size_t i = 1; i <= t; ++i) {
    size_t binomial = 1;
    for (size_t j = 0; j < i; ++j) {
      binomial = binomial * (n - j) / (j + 1);
    }
    volume += binomial * static_cast<size_t>(std::pow(q - 1, i));
  }

  return static_cast<size_t>(std::pow(q, k)) * volume <=
         static_cast<size_t>(std::pow(q, n));
}

// Perfect code check
template <typename ElementType>
bool IsPerfectCode(const AbstractLinearCode<ElementType>& code) {
  auto field = code.Field();
  size_t q = field->Order();
  size_t n = code.Length();
  size_t k = code.Dimension();
  size_t t = ErrorCorrectionCapability(code);

  // Calculate volume of Hamming sphere
  size_t volume = 1;
  for (size_t i = 1; i <= t; ++i) {
    size_t binomial = 1;
    for (size_t j = 0; j < i; ++j) {
      binomial = binomial * (n - j) / (j + 1);
    }
    volume += binomial * static_cast<size_t>(std::pow(q - 1, i));
  }

  return static_cast<size_t>(std::pow(q, k)) * volume ==
         static_cast<size_t>(std::pow(q, n));
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_HPP
