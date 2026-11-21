#ifndef XGALOIS_CODING_DECODER_BERLEKAMP_MASSEY_DECODER_HPP_
#define XGALOIS_CODING_DECODER_BERLEKAMP_MASSEY_DECODER_HPP_

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/grs_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

namespace detail {

// O(t^2) Berlekamp-Massey LFSR synthesis over arbitrary finite fields.
// Finds the unique minimal error locator Lambda(z) (with Lambda(0) = 1)
// annihilating the Forney modified syndrome sequence T_s, ..., T_{L-1}.
template <typename GaloisField>
inline std::optional<KeyEquationSolution<GaloisField>> SolveKeyEquationBM(
    const PolynomialDense<GaloisField> &mod_syn, std::size_t L, std::size_t s,
    std::size_t max_errors) {
  using StorageType = typename GaloisField::element_type;
  using Poly = PolynomialDense<GaloisField>;

  const auto field = mod_syn.Field();
  const std::size_t M = L - s;
  if (M == 0) {
    return KeyEquationSolution<GaloisField>{Poly::One(field), std::nullopt};
  }

  std::vector<StorageType> seq(M);
  for (std::size_t i = 0; i < M; ++i) {
    seq[i] = mod_syn.RawCoeff(s + i);
  }

  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();

  std::vector<StorageType> C = {one};
  std::vector<StorageType> B = {one};
  C.reserve(max_errors + 1);
  B.reserve(max_errors + 1);

  std::size_t L_reg = 0;
  std::size_t m = 1;
  StorageType b = one;

  for (std::size_t n = 0; n < M; ++n) {
    StorageType d = seq[n];
    const std::size_t lim = std::min(L_reg, C.size() - 1);
    for (std::size_t i = 1; i <= lim; ++i) {
      if (!field->IsZero(C[i]) && !field->IsZero(seq[n - i])) {
        d = field->Add(d, field->Mul(C[i], seq[n - i]));
      }
    }

    if (field->IsZero(d)) {
      ++m;
      continue;
    }

    const StorageType neg_factor = field->Neg(field->Div(d, b));
    const bool update_length = (2 * L_reg <= n);
    std::vector<StorageType> T_old;
    if (update_length) {
      T_old = C;
    }

    const std::size_t needed = B.size() + m;
    if (C.size() < needed) {
      C.resize(needed, zero);
    }
    field->FmaArray(std::span<StorageType>(C.data() + m, B.size()),
                    std::span<const StorageType>(B), neg_factor);

    if (update_length) {
      L_reg = n + 1 - L_reg;
      B = std::move(T_old);
      b = d;
      m = 1;
    } else {
      ++m;
    }
  }

  if (L_reg > max_errors) {
    return std::nullopt;
  }

  Poly lambda = Poly::FromRaw(std::move(C), field);
  if (lambda.Degree() != static_cast<int>(L_reg)) {
    return std::nullopt;
  }
  return KeyEquationSolution<GaloisField>{std::move(lambda), std::nullopt};
}

}  // namespace detail

// O(t^2) Berlekamp-Massey + Chien search + Forney decoder for Generalized
// Reed-Solomon codes, supporting full errors-and-erasures decoding up to
// 2e + s <= n - k and extended RS/GRS codes with alpha_j = 0.
template <typename GaloisField>
class BerlekampMasseyDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;

  explicit BerlekampMasseyDecoder(const AbstractCode<GaloisField> &code)
      : Base(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "BerlekampMasseyDecoder requires a GeneralizedReedSolomonCode");
    }
    max_errors_ = (grs_code_->MinimumDistance() - 1) / 2;
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return DecodeWithErasures(received_word, {});
  }

  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const override {
    return detail::DecodeGRSWithErasures(
        *grs_code_, received_word, erasure_positions, "BerlekampMassey",
        [](const PolynomialDense<GaloisField> &mod_syn, std::size_t L,
           std::size_t s, std::size_t max_errors) {
          return detail::SolveKeyEquationBM(mod_syn, L, s, max_errors);
        });
  }

  std::string Name() const override { return "BerlekampMassey"; }

  std::size_t DecodingRadius() const override { return max_errors_; }

  std::size_t MaxErrors() const noexcept { return max_errors_; }

 private:
  const GeneralizedReedSolomonCode<GaloisField> *grs_code_;
  std::size_t max_errors_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_BERLEKAMP_MASSEY_DECODER_HPP_
