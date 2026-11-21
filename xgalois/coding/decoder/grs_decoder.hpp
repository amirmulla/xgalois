#ifndef XGALOIS_CODING_GRS_HPP_
#include "xgalois/coding/grs.hpp"
#endif

#ifndef XGALOIS_CODING_DECODER_GRS_DECODER_HPP_
#define XGALOIS_CODING_DECODER_GRS_DECODER_HPP_

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace coding {

namespace detail {

template <typename GaloisField>
struct KeyEquationSolution {
  PolynomialDense<GaloisField> error_locator;
  std::optional<PolynomialDense<GaloisField>> error_evaluator;
};

// Computes all r = n - k GRS syndromes S_i = sum_{j=0}^{n-1} r_j * w_j *
// alpha_j^i (0 <= i < n - k) using the dual multipliers w_j via H * r^T.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type> ComputeGRSSyndromes(
    const GeneralizedReedSolomonCode<GaloisField> &code,
    const GFArray<GaloisField> &received) {
  using StorageType = typename GaloisField::element_type;
  if (code.Redundancy() == 0) {
    return {};
  }
  const auto syn = code.Syndrome(received);
  const auto span = syn.RawSpan();
  return std::vector<StorageType>(span.begin(), span.end());
}

template <typename GaloisField>
inline bool AreSyndromesZero(
    const GaloisField &field,
    std::span<const typename GaloisField::element_type> syndromes) {
  for (const auto &s : syndromes) {
    if (!field.IsZero(s)) {
      return false;
    }
  }
  return true;
}

template <typename GaloisField>
inline PolynomialDense<GaloisField> TruncatePolynomial(
    const PolynomialDense<GaloisField> &poly, std::size_t max_len) {
  using StorageType = typename GaloisField::element_type;
  if (max_len == 0 || poly.IsZero()) {
    return PolynomialDense<GaloisField>::Zero(poly.Field(), poly.GetVariable());
  }
  const auto &raw = poly.RawCoeffs();
  const std::size_t keep = std::min(max_len, raw.size());
  std::vector<StorageType> coeffs(
      raw.begin(), raw.begin() + static_cast<std::ptrdiff_t>(keep));
  return PolynomialDense<GaloisField>::FromRaw(std::move(coeffs), poly.Field(),
                                               poly.GetVariable());
}

// Builds Gamma(z) = prod_{p in erasure_positions} (1 - alpha_p * z) in O(s^2).
template <typename GaloisField>
inline PolynomialDense<GaloisField> BuildErasureLocator(
    const GeneralizedReedSolomonCode<GaloisField> &code,
    const std::vector<std::size_t> &erasure_positions) {
  using StorageType = typename GaloisField::element_type;
  const auto field = code.Field();
  std::vector<StorageType> coeffs;
  coeffs.reserve(erasure_positions.size() + 1);
  coeffs.push_back(field->MultiplicativeIdentity());

  const auto eval_span = code.EvaluationPoints().RawSpan();
  for (std::size_t pos : erasure_positions) {
    const StorageType neg_alpha = field->Neg(eval_span[pos]);
    const StorageType top = field->Mul(coeffs.back(), neg_alpha);
    for (std::size_t d = coeffs.size() - 1; d > 0; --d) {
      coeffs[d] = field->Add(coeffs[d], field->Mul(coeffs[d - 1], neg_alpha));
    }
    coeffs.push_back(top);
  }
  return PolynomialDense<GaloisField>::FromRaw(std::move(coeffs), field);
}

// Runs a single errors-and-erasures GRS decoding pass over non-zero evaluation
// points (normal mode: shifted == false, multipliers w_j; shifted mode:
// shifted == true, syndromes S'_i = S_{i+1}, effective multipliers w'_j =
// w_j * alpha_j, with final residual recovery of e_{j_0} at alpha_{j_0} = 0).
template <typename GaloisField, typename KeyEquationSolver>
inline std::optional<DecodingResult<GaloisField>> TryDecodeGRSPass(
    const GeneralizedReedSolomonCode<GaloisField> &code,
    const GFArray<GaloisField> &working_word,
    const GFArray<GaloisField> &original_received,
    std::span<const typename GaloisField::element_type> syndromes,
    const std::vector<std::size_t> &nonzero_erasures,
    const std::vector<bool> &is_erased, std::size_t total_erasure_count,
    std::size_t max_errors, bool shifted,
    const typename GaloisField::element_type &raw_s0,
    KeyEquationSolver &&solver) {
  using StorageType = typename GaloisField::element_type;
  using Poly = PolynomialDense<GaloisField>;

  const auto field = code.Field();
  const std::size_t L = syndromes.size();
  const std::size_t s = nonzero_erasures.size();
  if (s > L) {
    return std::nullopt;
  }

  const Poly gamma = BuildErasureLocator(code, nonzero_erasures);
  std::vector<StorageType> syn_vec(syndromes.begin(), syndromes.end());
  const Poly syn_poly = Poly::FromRaw(std::move(syn_vec), field);
  const Poly mod_syn =
      (s == 0) ? syn_poly : TruncatePolynomial(syn_poly * gamma, L);

  auto sol_opt = solver(mod_syn, L, s, max_errors);
  if (!sol_opt.has_value()) {
    return std::nullopt;
  }

  Poly lambda = std::move(sol_opt->error_locator);
  if (lambda.IsZero()) {
    return std::nullopt;
  }
  const StorageType l0 = lambda.RawCoeff(0);
  if (field->IsZero(l0)) {
    return std::nullopt;
  }

  Poly omega;
  if (sol_opt->error_evaluator.has_value()) {
    omega = std::move(*sol_opt->error_evaluator);
    if (!field->IsOne(l0)) {
      const auto inv_l0 =
          GaloisFieldElement<GaloisField>::FromInternal(field->Inv(l0), field);
      lambda *= inv_l0;
      omega *= inv_l0;
    }
  } else {
    if (!field->IsOne(l0)) {
      const auto inv_l0 =
          GaloisFieldElement<GaloisField>::FromInternal(field->Inv(l0), field);
      lambda *= inv_l0;
    }
    omega = TruncatePolynomial(mod_syn * lambda, L);
  }

  const int deg_lambda = lambda.Degree();
  if (deg_lambda < 0 || static_cast<std::size_t>(deg_lambda) > max_errors) {
    return std::nullopt;
  }

  // Chien search over cached inverses of non-zero evaluation points.
  std::vector<std::size_t> error_positions;
  if (deg_lambda > 0) {
    const auto root_sub_indices =
        utils::ChienSearchIndices(lambda, code.NonZeroEvalInverses());
    if (root_sub_indices.size() != static_cast<std::size_t>(deg_lambda)) {
      return std::nullopt;
    }
    error_positions.reserve(root_sub_indices.size());
    const auto &nz_indices = code.NonZeroEvalIndices();
    for (std::size_t sub_idx : root_sub_indices) {
      const std::size_t pos = nz_indices[sub_idx];
      if (is_erased[pos]) {
        return std::nullopt;
      }
      error_positions.push_back(pos);
    }
  }

  // Combined error+erasure locator Psi(z) = Lambda(z) * Gamma(z) and Forney.
  const Poly psi = (s == 0) ? lambda : (lambda * gamma);
  const Poly psi_deriv = psi.Derivative();

  std::vector<std::size_t> all_positions;
  all_positions.reserve(error_positions.size() + nonzero_erasures.size());
  all_positions.insert(all_positions.end(), error_positions.begin(),
                       error_positions.end());
  all_positions.insert(all_positions.end(), nonzero_erasures.begin(),
                       nonzero_erasures.end());

  GFArray<GaloisField> corrected = working_word;
  const auto eval_span = code.EvaluationPoints().RawSpan();
  const auto dual_span = code.DualMultipliers().RawSpan();
  StorageType zero_point_sum = field->AdditiveIdentity();

  for (std::size_t pos : all_positions) {
    const StorageType &inv_alpha = code.EvalInverseAt(pos);
    const StorageType num = omega.EvalRaw(inv_alpha);
    const StorageType den = psi_deriv.EvalRaw(inv_alpha);
    if (field->IsZero(den)) {
      return std::nullopt;
    }
    // Normal mode:  e_p = -(alpha_p * Omega(alpha_p^-1)) / (w_p *
    // Psi'(alpha_p^-1)) Shifted mode: w'_p = w_p * alpha_p, so alpha_p cancels:
    //               e_p = -Omega'(alpha_p^-1) / (w_p * Psi'(alpha_p^-1))
    StorageType numer = field->Neg(num);
    if (!shifted) {
      numer = field->Mul(numer, eval_span[pos]);
    }
    const StorageType denom = field->Mul(dual_span[pos], den);
    const StorageType e_val = field->Div(numer, denom);
    corrected.Raw(pos) = field->Sub(corrected.Raw(pos), e_val);
    if (shifted) {
      zero_point_sum =
          field->Add(zero_point_sum, field->Mul(e_val, dual_span[pos]));
    }
  }

  if (shifted) {
    const std::size_t j0 = *code.ZeroEvaluationIndex();
    const StorageType rem_s0 = field->Sub(raw_s0, zero_point_sum);
    const StorageType e_j0 = field->Div(rem_s0, dual_span[j0]);
    corrected.Raw(j0) = field->Sub(corrected.Raw(j0), e_j0);
  }

  // Always verify mathematical codeword membership and bounded-distance radius.
  if (!code.Contains(corrected)) {
    return std::nullopt;
  }

  const std::size_t n = code.Length();
  std::size_t actual_errors = 0;
  for (std::size_t j = 0; j < n; ++j) {
    if (!is_erased[j] && corrected.Raw(j) != original_received.Raw(j)) {
      ++actual_errors;
    }
  }
  if (2 * actual_errors + total_erasure_count > code.Redundancy()) {
    return std::nullopt;
  }

  auto message = code.Unencode(corrected);
  return DecodingResult<GaloisField>::Success(std::move(corrected),
                                              std::move(message), actual_errors,
                                              total_erasure_count);
}

// Unified GRS errors-and-erasures decoding driver shared by PGZ,
// Berlekamp-Massey, and Sugiyama decoders. Handles alpha_j = 0 cleanly
// whether the zero evaluation point is erased, corrupted, or error-free.
template <typename GaloisField, typename KeyEquationSolver>
inline DecodingResult<GaloisField> DecodeGRSWithErasures(
    const GeneralizedReedSolomonCode<GaloisField> &code,
    const GFArray<GaloisField> &received,
    const std::vector<std::size_t> &erasure_positions,
    const std::string &decoder_name, KeyEquationSolver &&solver) {
  using Result = DecodingResult<GaloisField>;
  using StorageType = typename GaloisField::element_type;

  const auto field = code.Field();
  const std::size_t n = code.Length();
  const std::size_t r = code.Redundancy();

  if (!received.HasField() || received.Field() != field) {
    throw std::invalid_argument(decoder_name +
                                ": received word field mismatch");
  }
  if (received.Dimension() != 1 || received.Size() != n) {
    throw std::invalid_argument(
        decoder_name +
        ": received word must be 1D with size equal to code length");
  }

  const std::size_t s = erasure_positions.size();
  std::vector<bool> is_erased(n, false);
  for (std::size_t pos : erasure_positions) {
    if (pos >= n) {
      return Result::Failure(decoder_name + ": erasure position " +
                             std::to_string(pos) + " out of range [0, " +
                             std::to_string(n) + ")");
    }
    if (is_erased[pos]) {
      return Result::Failure(decoder_name + ": duplicate erasure position " +
                             std::to_string(pos));
    }
    is_erased[pos] = true;
  }

  if (s > r) {
    return Result::Failure(decoder_name + ": number of erasures (" +
                           std::to_string(s) + ") exceeds redundancy (" +
                           std::to_string(r) + ")");
  }

  GFArray<GaloisField> working = received;
  const StorageType zero = field->AdditiveIdentity();
  for (std::size_t pos : erasure_positions) {
    working.Raw(pos) = zero;
  }

  const std::vector<StorageType> syndromes = ComputeGRSSyndromes(code, working);
  if (s == 0 && AreSyndromesZero(*field, syndromes)) {
    return Result::Success(received, code.Unencode(received), 0, 0);
  }

  const std::size_t e_max = (r - s) / 2;
  const auto zero_idx = code.ZeroEvaluationIndex();

  // Case A: zero evaluation point is itself erased -> shift syndromes by 1 so
  // all remaining erased and erroneous points have non-zero evaluation points.
  if (zero_idx.has_value() && is_erased[*zero_idx]) {
    std::vector<std::size_t> nonzero_erasures;
    nonzero_erasures.reserve(s - 1);
    for (std::size_t pos : erasure_positions) {
      if (pos != *zero_idx) {
        nonzero_erasures.push_back(pos);
      }
    }
    const std::span<const StorageType> shifted_syn(syndromes.data() + 1, r - 1);
    if (auto res = TryDecodeGRSPass(code, working, received, shifted_syn,
                                    nonzero_erasures, is_erased, s, e_max,
                                    /*shifted=*/true, syndromes[0], solver)) {
      return std::move(*res);
    }
    return Result::Failure(decoder_name +
                           ": uncorrectable error/erasure pattern");
  }

  // Case B: try standard decoding assuming no unknown error at alpha_j = 0.
  const std::span<const StorageType> full_syn(syndromes);
  if (auto res = TryDecodeGRSPass(code, working, received, full_syn,
                                  erasure_positions, is_erased, s, e_max,
                                  /*shifted=*/false, syndromes[0], solver)) {
    return std::move(*res);
  }

  // Case C: if an evaluation point is 0 (and not erased) and we have budget
  // for at least 1 unknown error, try shifted syndromes assuming an error sits
  // at alpha_{j_0} = 0 plus up to e_max - 1 errors at non-zero points.
  if (zero_idx.has_value() && e_max >= 1) {
    const std::span<const StorageType> shifted_syn(syndromes.data() + 1, r - 1);
    if (auto res = TryDecodeGRSPass(code, working, received, shifted_syn,
                                    erasure_positions, is_erased, s, e_max - 1,
                                    /*shifted=*/true, syndromes[0], solver)) {
      return std::move(*res);
    }
  }

  return Result::Failure(decoder_name +
                         ": uncorrectable error/erasure pattern");
}

// Peterson-Gorenstein-Zierler (PGZ) key-equation solver:
// Loops descending from nu = max_errors down to 1, solves the Hankel linear
// system M_nu * x = -a via linalg::solve, and constructs
// Lambda(z) = 1 + Lambda_1 z + ... + Lambda_nu z^nu with Lambda_k = x_{nu - k}.
template <typename GaloisField>
inline std::optional<KeyEquationSolution<GaloisField>> SolveKeyEquationPGZ(
    const PolynomialDense<GaloisField> &mod_syn, std::size_t L, std::size_t s,
    std::size_t max_errors) {
  using StorageType = typename GaloisField::element_type;
  using Poly = PolynomialDense<GaloisField>;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  const auto field = mod_syn.Field();
  const std::size_t M = L - s;
  if (M == 0) {
    return KeyEquationSolution<GaloisField>{Poly::One(field), std::nullopt};
  }

  std::vector<StorageType> seq(M);
  bool all_zero = true;
  for (std::size_t i = 0; i < M; ++i) {
    seq[i] = mod_syn.RawCoeff(s + i);
    if (!field->IsZero(seq[i])) {
      all_zero = false;
    }
  }
  if (all_zero) {
    return KeyEquationSolution<GaloisField>{Poly::One(field), std::nullopt};
  }
  if (max_errors == 0) {
    return std::nullopt;
  }

  const std::size_t nu_max = std::min(max_errors, M / 2);
  for (std::size_t nu = nu_max; nu >= 1; --nu) {
    RawArray mat_raw = RawArray::from_shape(ShapeType{nu, nu});
    RawArray rhs_raw = RawArray::from_shape(ShapeType{nu});
    for (std::size_t i = 0; i < nu; ++i) {
      for (std::size_t j = 0; j < nu; ++j) {
        mat_raw(i, j) = seq[i + j];
      }
      rhs_raw(i) = field->Neg(seq[i + nu]);
    }

    auto mat = GFArray<GaloisField>::FromRaw(field, std::move(mat_raw));
    auto rhs = GFArray<GaloisField>::FromRaw(field, std::move(rhs_raw));
    try {
      auto sol = linalg::solve(mat, rhs);
      std::vector<StorageType> coeffs(nu + 1);
      coeffs[0] = field->MultiplicativeIdentity();
      for (std::size_t k = 1; k <= nu; ++k) {
        coeffs[k] = sol.Raw(nu - k);
      }
      Poly lambda = Poly::FromRaw(std::move(coeffs), field);
      if (lambda.Degree() == static_cast<int>(nu)) {
        return KeyEquationSolution<GaloisField>{std::move(lambda),
                                                std::nullopt};
      }
    } catch (const std::runtime_error &) {
      // Singular Hankel matrix M_nu => actual error count is strictly < nu.
      continue;
    }
  }

  return std::nullopt;
}

}  // namespace detail

// Peterson-Gorenstein-Zierler (PGZ) decoder for Generalized Reed-Solomon codes.
template <typename GaloisField>
class GRSDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;

  explicit GRSDecoder(const AbstractCode<GaloisField> &code) : Base(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GRSDecoder requires a GeneralizedReedSolomonCode");
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
        *grs_code_, received_word, erasure_positions, "PGZ",
        [](const PolynomialDense<GaloisField> &mod_syn, std::size_t L,
           std::size_t s, std::size_t max_errors) {
          return detail::SolveKeyEquationPGZ(mod_syn, L, s, max_errors);
        });
  }

  std::string Name() const override { return "PGZ"; }

  std::size_t DecodingRadius() const override { return max_errors_; }

  std::size_t MaxErrors() const noexcept { return max_errors_; }

 private:
  const GeneralizedReedSolomonCode<GaloisField> *grs_code_;
  std::size_t max_errors_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_GRS_DECODER_HPP_
