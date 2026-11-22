#ifndef XGALOIS_CODING_GRS_HPP_
#include "xgalois/coding/grs.hpp"
#endif

#ifndef XGALOIS_CODING_DECODER_SUGIYAMA_DECODER_HPP_
#define XGALOIS_CODING_DECODER_SUGIYAMA_DECODER_HPP_

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/grs_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/poly/half_gcd.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

inline constexpr std::size_t kSugiyamaHalfGcdThreshold = 256;

namespace detail {

// Sugiyama (Extended Euclidean Algorithm) key-equation solver.
// Runs polynomial EEA on (z^L, T(z)) stopping as soon as deg(r_i) < max_errors
// + s, then normalizes by t_i(0)^-1 to obtain Lambda(0) = 1 and Omega(z)
// directly. Dispatches to O(M(L) log L) Half-GCD (xg::PartialExtendedGcd) when
// L >= half_gcd_threshold, and uses classical O(t^2) EEA below.
template <typename GaloisField>
inline std::optional<KeyEquationSolution<GaloisField>> SolveKeyEquationSugiyama(
    const PolynomialDense<GaloisField> &mod_syn, std::size_t L, std::size_t s,
    std::size_t max_errors,
    std::size_t half_gcd_threshold = kSugiyamaHalfGcdThreshold) {
  using Poly = PolynomialDense<GaloisField>;

  const auto field = mod_syn.Field();
  const auto one = field->MultiplicativeIdentity();
  const int stop_deg = static_cast<int>(max_errors + s);

  Poly r_curr;
  Poly t_curr;

  if (L >= half_gcd_threshold) {
    Poly x_l = Poly::Monomial(field, L, one);
    auto [r_out, u_out, t_out] =
        xg::PartialExtendedGcd(x_l, mod_syn, static_cast<int64_t>(stop_deg));
    (void)u_out;
    r_curr = std::move(r_out);
    t_curr = std::move(t_out);
  } else {
    Poly r_prev = Poly::Monomial(field, L, one);
    r_curr = mod_syn;
    Poly t_prev = Poly::Zero(field);
    t_curr = Poly::One(field);

    while (!r_curr.IsZero() && r_curr.Degree() >= stop_deg) {
      auto [q, r_next] = r_prev.DivRem(r_curr);
      Poly t_next = t_prev - q * t_curr;
      r_prev = std::move(r_curr);
      r_curr = std::move(r_next);
      t_prev = std::move(t_curr);
      t_curr = std::move(t_next);
    }
  }

  if (t_curr.IsZero()) {
    return std::nullopt;
  }
  const auto t0 = t_curr.RawCoeff(0);
  if (field->IsZero(t0)) {
    return std::nullopt;
  }
  if (t_curr.Degree() > static_cast<int>(max_errors)) {
    return std::nullopt;
  }

  const auto inv_t0 =
      GaloisFieldElement<GaloisField>::FromInternal(field->Inv(t0), field);
  Poly lambda = t_curr * inv_t0;
  Poly omega = r_curr * inv_t0;
  return KeyEquationSolution<GaloisField>{std::move(lambda), std::move(omega)};
}

}  // namespace detail

// Sugiyama (Extended Euclidean Algorithm / Half-GCD) + Chien search + Forney
// decoder for Generalized Reed-Solomon codes, supporting full
// errors-and-erasures decoding up to 2e + s <= n - k and alpha_j = 0.
template <typename GaloisField>
class SugiyamaDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;

  explicit SugiyamaDecoder(
      const AbstractCode<GaloisField> &code,
      std::size_t half_gcd_threshold = kSugiyamaHalfGcdThreshold)
      : Base(code), half_gcd_threshold_(half_gcd_threshold) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "SugiyamaDecoder requires a GeneralizedReedSolomonCode");
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
    const std::size_t threshold = half_gcd_threshold_;
    return detail::DecodeGRSWithErasures(
        *grs_code_, received_word, erasure_positions, "Sugiyama",
        [threshold](const PolynomialDense<GaloisField> &mod_syn, std::size_t L,
                    std::size_t s, std::size_t max_errors) {
          return detail::SolveKeyEquationSugiyama(mod_syn, L, s, max_errors,
                                                  threshold);
        });
  }

  std::string Name() const override { return "Sugiyama"; }

  std::size_t DecodingRadius() const override { return max_errors_; }

  std::size_t MaxErrors() const noexcept { return max_errors_; }

  std::size_t HalfGcdThreshold() const noexcept { return half_gcd_threshold_; }

 private:
  const GeneralizedReedSolomonCode<GaloisField> *grs_code_;
  std::size_t max_errors_;
  std::size_t half_gcd_threshold_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_SUGIYAMA_DECODER_HPP_
