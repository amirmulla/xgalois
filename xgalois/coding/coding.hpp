#ifndef XGALOIS_CODING_HPP
#define XGALOIS_CODING_HPP

// Base classes
#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/coding/abstract_linear_code.hpp"

// Encoder/Decoder framework
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/encoder/encoder.hpp"

// Code implementations
#include "xgalois/coding/cyclic_code.hpp"
#include "xgalois/coding/grs.hpp"

namespace xg {
namespace coding {

// Convenience functions for creating codes

// Create a standard Reed-Solomon code
template <typename GaloisField>
std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>> CreateReedSolomonCode(
    std::shared_ptr<GaloisFieldBase<GaloisField>> field, size_t length,
    size_t dimension, const element_type& primitive_element) {
  return GeneralizedReedSolomonCode<GaloisField>::ReedSolomon(
      field, length, dimension, primitive_element);
}

// Create a standard Reed-Solomon code with consecutive evaluation points
template <typename GaloisField>
std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>> CreateReedSolomonCode(
    std::shared_ptr<GaloisFieldBase<GaloisField>> field, size_t length,
    size_t dimension) {
  // Use primitive element (assuming it's 2 for binary extension fields)
  element_type primitive_element(2);
  return CreateReedSolomonCode(field, length, dimension, primitive_element);
}

// Create a generalized Reed-Solomon code
template <typename GaloisField>
std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>>
CreateGeneralizedReedSolomonCode(
    std::shared_ptr<GaloisFieldBase<GaloisField>> field,
    const xt::xarray<GaloisField>& evaluation_points,
    const xt::xarray<GaloisField>& column_multipliers, size_t dimension) {
  return std::make_unique<GeneralizedReedSolomonCode<GaloisField>>(
      field, evaluation_points, column_multipliers, dimension);
}

// Create a cyclic code from generator polynomial
template <typename GaloisField>
std::unique_ptr<CyclicCode<GaloisField>> CreateCyclicCode(
    std::shared_ptr<GaloisFieldBase<GaloisField>> field, size_t length,
    const xg::poly::PolyDense<GaloisField>& generator_poly) {
  return std::make_unique<CyclicCode<GaloisField>>(field, length,
                                                   generator_poly);
}

// Error correction capability helper
template <typename GaloisField>
size_t ErrorCorrectionCapability(const AbstractCode<GaloisField>& code) {
  return (code.MinimumDistance() - 1) / 2;
}

// Code rate helper
template <typename GaloisField>
double CodeRate(const AbstractLinearCode<GaloisField>& code) {
  return code.Rate();
}

// Redundancy helper
template <typename GaloisField>
size_t Redundancy(const AbstractLinearCode<GaloisField>& code) {
  return code.Redundancy();
}

// Check if code is MDS (Maximum Distance Separable)
template <typename GaloisField>
bool IsMDS(const AbstractLinearCode<GaloisField>& code) {
  return code.MinimumDistance() == code.Length() - code.Dimension() + 1;
}

// Singleton bound check
template <typename GaloisField>
bool SatisfiesSingletonBound(const AbstractLinearCode<GaloisField>& code) {
  return code.MinimumDistance() <= code.Length() - code.Dimension() + 1;
}

// Hamming bound check
template <typename GaloisField>
bool SatisfiesHammingBound(const AbstractLinearCode<GaloisField>& code) {
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
template <typename GaloisField>
bool IsPerfectCode(const AbstractLinearCode<GaloisField>& code) {
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
