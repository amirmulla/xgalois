#ifndef XGALOIS_FIELD_GF_CONCEPTS_HPP_
#define XGALOIS_FIELD_GF_CONCEPTS_HPP_

#include <concepts>
#include <cstdint>
#include <iostream>
#include <random>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

namespace xg {

template <typename T>
concept UnsignedFieldWord =
    std::unsigned_integral<T> && !std::same_as<std::remove_cv_t<T>, bool>;

template <typename F>
concept FiniteField =
    requires(const F& field, const typename F::element_type& a,
             const typename F::element_type& b, uint64_t exp,
             std::span<typename F::element_type> dst,
             std::span<const typename F::element_type> src) {
      typename F::element_type;
      typename F::storage_type;

      { field.Characteristic() } -> std::convertible_to<uint64_t>;
      { field.Order() } -> std::convertible_to<uint64_t>;
      { field.Degree() } -> std::convertible_to<uint32_t>;

      { field.Add(a, b) } -> std::same_as<typename F::element_type>;
      { field.Sub(a, b) } -> std::same_as<typename F::element_type>;
      { field.Mul(a, b) } -> std::same_as<typename F::element_type>;
      { field.Div(a, b) } -> std::same_as<typename F::element_type>;
      { field.Neg(a) } -> std::same_as<typename F::element_type>;
      { field.Inv(a) } -> std::same_as<typename F::element_type>;
      { field.Pow(a, exp) } -> std::same_as<typename F::element_type>;
      { field.Sqrt(a) } -> std::same_as<typename F::element_type>;
      { field.Log(a) } -> std::convertible_to<uint64_t>;

      { field.AdditiveIdentity() } -> std::same_as<typename F::element_type>;
      {
        field.MultiplicativeIdentity()
      } -> std::same_as<typename F::element_type>;
      { field.IsZero(a) } -> std::convertible_to<bool>;
      { field.IsOne(a) } -> std::convertible_to<bool>;

      { field.GetElementValue(a) } -> std::same_as<typename F::element_type>;
      { field.SetElementValue(a) } -> std::same_as<typename F::element_type>;
      { field.ToString(a) } -> std::same_as<std::string>;

      { field.AddArray(dst, src, src) };
      { field.MulArray(dst, src, src) };
      { field.FmaArray(dst, src, a) };
    };

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_CONCEPTS_HPP_
