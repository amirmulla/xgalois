#ifndef XGALOIS_CODING_HPP
#define XGALOIS_CODING_HPP

#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/cyclic_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/grs.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>> CreateReedSolomonCode(
    std::shared_ptr<GaloisField> field, size_t length, size_t dimension,
    const xg::GaloisFieldElement<GaloisField>& primitive_element) {
  return GeneralizedReedSolomonCode<GaloisField>::ReedSolomon(
      field, length, dimension, primitive_element);
}

template <typename GaloisField>
std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>> CreateReedSolomonCode(
    std::shared_ptr<GaloisField> field, size_t length, size_t dimension) {
  xg::GaloisFieldElement<GaloisField> primitive_element(2, field);
  return CreateReedSolomonCode(field, length, dimension, primitive_element);
}

template <typename GaloisField>
std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>>
CreateGeneralizedReedSolomonCode(
    std::shared_ptr<GaloisField> field,
    const xt::xarray<GaloisField>& evaluation_points,
    const xt::xarray<GaloisField>& column_multipliers, size_t dimension) {
  return std::make_unique<GeneralizedReedSolomonCode<GaloisField>>(
      field, evaluation_points, column_multipliers, dimension);
}

template <typename GaloisField>
std::unique_ptr<CyclicCode<GaloisField>> CreateCyclicCode(
    std::shared_ptr<GaloisField> field, size_t length,
    const xg::PolynomialDense<GaloisField>& generator_poly) {
  return std::make_unique<CyclicCode<GaloisField>>(field, length,
                                                   generator_poly);
}

template <typename GaloisField>
size_t ErrorCorrectionCapability(const AbstractCode<GaloisField>& code) {
  return (code.MinimumDistance() - 1) / 2;
}

template <typename GaloisField>
double CodeRate(const AbstractLinearCode<GaloisField>& code) {
  return code.Rate();
}

template <typename GaloisField>
size_t Redundancy(const AbstractLinearCode<GaloisField>& code) {
  return code.Redundancy();
}

template <typename GaloisField>
bool IsMDS(const AbstractLinearCode<GaloisField>& code) {
  return code.MinimumDistance() == code.Length() - code.Dimension() + 1;
}

template <typename GaloisField>
bool SatisfiesSingletonBound(const AbstractLinearCode<GaloisField>& code) {
  return code.MinimumDistance() <= code.Length() - code.Dimension() + 1;
}

template <typename GaloisField>
bool SatisfiesHammingBound(const AbstractLinearCode<GaloisField>& code) {
  auto field = code.Field();
  size_t q = field->Order();
  size_t n = code.Length();
  size_t k = code.Dimension();
  size_t t = ErrorCorrectionCapability(code);

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

template <typename GaloisField>
bool IsPerfectCode(const AbstractLinearCode<GaloisField>& code) {
  auto field = code.Field();
  size_t q = field->Order();
  size_t n = code.Length();
  size_t k = code.Dimension();
  size_t t = ErrorCorrectionCapability(code);

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

#endif
