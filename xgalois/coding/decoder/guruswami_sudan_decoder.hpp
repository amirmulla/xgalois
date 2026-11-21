#ifndef XGALOIS_CODING_DECODER_GURUSWAMI_SUDAN_DECODER_HPP_
#define XGALOIS_CODING_DECODER_GURUSWAMI_SUDAN_DECODER_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace coding {

namespace detail {

__extension__ using Uint128 = unsigned __int128;

// Exact floor(sqrt(a * b)) without floating-point rounding or 64-bit overflow.
inline std::size_t IntegerSqrtProduct(std::size_t a, std::size_t b) noexcept {
  const Uint128 prod = static_cast<Uint128>(a) * static_cast<Uint128>(b);
  if (prod == 0) {
    return 0;
  }
  Uint128 lo = 1;
  Uint128 hi = static_cast<Uint128>(std::max(a, b));
  while (lo < hi) {
    const Uint128 mid = lo + (hi - lo + 1) / 2;
    if (mid <= prod / mid) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }
  return static_cast<std::size_t>(lo);
}

// Number of monomials x^u y^v with u >= 0, 0 <= v <= list_size, and
// u + weight_y * v <= degree_bound.
inline Uint128 WeightedMonomialCount(std::size_t degree_bound,
                                     std::size_t weight_y,
                                     std::size_t list_size) noexcept {
  const Uint128 d = static_cast<Uint128>(degree_bound);
  const Uint128 w = static_cast<Uint128>(weight_y);
  const Uint128 v_max = (weight_y == 0)
                            ? static_cast<Uint128>(list_size)
                            : std::min(static_cast<Uint128>(list_size), d / w);
  return (v_max + 1) * (d + 1) - w * (v_max * (v_max + 1) / 2);
}

// Total number of homogeneous Hasse-derivative interpolation constraints:
// n * C(m + 1, 2) = n * m * (m + 1) / 2.
inline Uint128 InterpolationConstraintCount(std::size_t n,
                                            std::size_t m) noexcept {
  const Uint128 um = static_cast<Uint128>(m);
  return static_cast<Uint128>(n) * (um * (um + 1) / 2);
}

// Smallest (1, k - 1)-weighted degree bound D such that
// N_{1, k-1}(D, L) > n * C(m + 1, 2).
inline std::size_t ComputeWeightedDegreeBound(std::size_t n, std::size_t k,
                                              std::size_t m, std::size_t L) {
  const std::size_t w = k - 1;
  const Uint128 constraints = InterpolationConstraintCount(n, m);
  std::size_t lo = 0;
  std::size_t hi = 1;
  while (WeightedMonomialCount(hi, w, L) <= constraints) {
    if (hi > std::numeric_limits<std::size_t>::max() / 2) {
      throw std::overflow_error(
          "GuruswamiSudanDecoder: weighted degree bound overflow");
    }
    hi *= 2;
  }
  while (lo < hi) {
    const std::size_t mid = lo + (hi - lo) / 2;
    if (WeightedMonomialCount(mid, w, L) > constraints) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  return lo;
}

// Computes the minimal list size L achieving the maximum guaranteed decoding
// radius possible for a given multiplicity m on an [n, k] GRS code.
inline std::size_t ComputeOptimalListSize(std::size_t n, std::size_t k,
                                          std::size_t m) {
  const std::size_t w = k - 1;
  const Uint128 constraints = InterpolationConstraintCount(n, m);
  if (w == 0) {
    // For k = 1, D = m - 1 gives floor(D / m) = 0 -> radius n - 1.
    // N_{1, 0}(m - 1, L) = (L + 1) * m > n * m * (m + 1) / 2
    // <=> L >= floor(n * (m + 1) / 2).
    const Uint128 raw_l =
        (static_cast<Uint128>(n) * static_cast<Uint128>(m + 1)) / 2;
    return std::max<std::size_t>(1, static_cast<std::size_t>(raw_l));
  }

  // Find minimal unconstrained D_opt with N_{1, w}(D_opt, floor(D_opt / w)) >
  // constraints.
  std::size_t lo = 0;
  std::size_t hi = 1;
  while (WeightedMonomialCount(hi, w, hi / w) <= constraints) {
    if (hi > std::numeric_limits<std::size_t>::max() / 2) {
      throw std::overflow_error(
          "GuruswamiSudanDecoder: unconstrained degree bound overflow");
    }
    hi *= 2;
  }
  while (lo < hi) {
    const std::size_t mid = lo + (hi - lo) / 2;
    if (WeightedMonomialCount(mid, w, mid / w) > constraints) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  const std::size_t d_opt = lo;

  // Any D <= d_target achieves the same guaranteed radius n - 1 - floor(D / m).
  const std::size_t min_agreements = d_opt / m + 1;
  const std::size_t d_target =
      (min_agreements > n) ? d_opt : (m * min_agreements - 1);
  const std::size_t l_max = std::max<std::size_t>(1, d_target / w);

  std::size_t l_lo = 1;
  std::size_t l_hi = l_max;
  while (l_lo < l_hi) {
    const std::size_t mid = l_lo + (l_hi - l_lo) / 2;
    if (WeightedMonomialCount(d_target, w, mid) > constraints) {
      l_hi = mid;
    } else {
      l_lo = mid + 1;
    }
  }
  return l_lo;
}

template <typename GaloisField>
inline bool CanonicalPolyLess(const PolynomialDense<GaloisField> &a,
                              const PolynomialDense<GaloisField> &b) {
  if (a.Degree() != b.Degree()) {
    return a.Degree() < b.Degree();
  }
  if (a.IsZero()) {
    return false;
  }
  const auto field = a.Field();
  for (int i = a.Degree(); i >= 0; --i) {
    const auto va =
        field->GetElementValue(a.RawCoeff(static_cast<std::size_t>(i)));
    const auto vb =
        field->GetElementValue(b.RawCoeff(static_cast<std::size_t>(i)));
    if (va != vb) {
      return va < vb;
    }
  }
  return false;
}

// Dense bivariate polynomial Q(x, y) = sum_{v} q_v(x) y^v stored in raw
// internal field elements.
template <typename StorageType>
struct BivariatePoly {
  std::vector<std::vector<StorageType>> y_coeffs;

  bool IsZero() const noexcept { return y_coeffs.empty(); }

  int YDegree() const noexcept { return static_cast<int>(y_coeffs.size()) - 1; }

  template <typename GaloisField>
  void Trim(const GaloisField &field) {
    for (auto &row : y_coeffs) {
      while (!row.empty() && field.IsZero(row.back())) {
        row.pop_back();
      }
    }
    while (!y_coeffs.empty() && y_coeffs.back().empty()) {
      y_coeffs.pop_back();
    }
  }

  // Computes *this += scalar * other in-place (requires this != &other).
  template <typename GaloisField>
  void FmaInPlace(const GaloisField &field, const BivariatePoly &other,
                  const StorageType &scalar) {
    if (field.IsZero(scalar) || other.y_coeffs.empty()) {
      return;
    }
    const StorageType zero = field.AdditiveIdentity();
    if (y_coeffs.size() < other.y_coeffs.size()) {
      y_coeffs.resize(other.y_coeffs.size());
    }
    for (std::size_t v = 0; v < other.y_coeffs.size(); ++v) {
      const auto &src = other.y_coeffs[v];
      if (src.empty()) {
        continue;
      }
      auto &dst = y_coeffs[v];
      if (dst.size() < src.size()) {
        dst.resize(src.size(), zero);
      }
      field.FmaArray(std::span<StorageType>(dst.data(), src.size()),
                     std::span<const StorageType>(src), scalar);
    }
    Trim(field);
  }

  // Computes *this *= (x - alpha) in-place.
  template <typename GaloisField>
  void MulByXMinusAlphaInPlace(const GaloisField &field,
                               const StorageType &alpha) {
    const StorageType zero = field.AdditiveIdentity();
    const StorageType neg_alpha = field.Neg(alpha);
    const bool alpha_is_zero = field.IsZero(neg_alpha);
    for (auto &row : y_coeffs) {
      if (row.empty()) {
        continue;
      }
      const std::size_t s = row.size();
      row.push_back(row[s - 1]);
      if (alpha_is_zero) {
        for (std::size_t u = s - 1; u > 0; --u) {
          row[u] = row[u - 1];
        }
        row[0] = zero;
      } else {
        for (std::size_t u = s - 1; u > 0; --u) {
          row[u] = field.Add(row[u - 1], field.Mul(neg_alpha, row[u]));
        }
        row[0] = field.Mul(neg_alpha, row[0]);
      }
    }
  }
};

}  // namespace detail

// Guruswami-Sudan polynomial-time list decoder for Generalized Reed-Solomon
// (GRS) and Reed-Solomon (RS) codes over arbitrary finite fields (including
// characteristic 2, odd prime fields, and general odd-characteristic extension
// fields).
//
// Decodes beyond the classical half-minimum-distance bound floor((d - 1) / 2)
// all the way up to the Johnson radius:
//   tau_J = n - 1 - floor(sqrt(n * (k - 1))).
//
// Step 1: Koetter's iterative Groebner-basis bivariate interpolation using
// exact Hasse derivatives modulo field characteristic.
// Step 2: Roth-Ruckenstein recursive root-finding over F_q[x] via
// utils::FindRoots, returning all codewords within DecodingRadius().
template <typename GaloisField>
class GuruswamiSudanDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using polynomial_type = PolynomialDense<GaloisField>;

  // Constructs a Guruswami-Sudan decoder with the given multiplicity m >= 1
  // (m = 1 corresponds to Sudan list decoding) and automatically selects the
  // minimal list size L achieving the maximum guaranteed radius for m.
  explicit GuruswamiSudanDecoder(const AbstractCode<GaloisField> &code,
                                 std::size_t multiplicity = 1)
      : Base(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder requires a GeneralizedReedSolomonCode");
    }
    if (multiplicity == 0) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder: multiplicity must be >= 1");
    }
    const std::size_t n = grs_code_->Length();
    const std::size_t k = grs_code_->Dimension();
    const std::size_t optimal_list_size =
        detail::ComputeOptimalListSize(n, k, multiplicity);
    InitializeParameters(multiplicity, optimal_list_size);
  }

  // Constructs a Guruswami-Sudan decoder with explicit multiplicity m >= 1 and
  // maximum y-degree (list size bound) L >= 1.
  GuruswamiSudanDecoder(const AbstractCode<GaloisField> &code,
                        std::size_t multiplicity, std::size_t list_size)
      : Base(code) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder requires a GeneralizedReedSolomonCode");
    }
    if (multiplicity == 0) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder: multiplicity must be >= 1");
    }
    if (list_size == 0) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder: list_size must be >= 1");
    }
    InitializeParameters(multiplicity, list_size);
  }

  // Automatically computes the minimal multiplicity m >= 1 and minimal list
  // size L >= 1 guaranteeing list decoding up to target_errors (requires
  // target_errors <= JohnsonRadius()).
  static GuruswamiSudanDecoder FromTargetRadius(
      const AbstractCode<GaloisField> &code, std::size_t target_errors) {
    const auto *grs =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder::FromTargetRadius requires a "
          "GeneralizedReedSolomonCode");
    }
    const std::size_t n = grs->Length();
    const std::size_t k = grs->Dimension();
    const std::size_t johnson = ComputeJohnsonRadius(n, k);
    if (target_errors > johnson) {
      throw std::invalid_argument(
          "GuruswamiSudanDecoder::FromTargetRadius: target_errors (" +
          std::to_string(target_errors) + ") exceeds Johnson radius (" +
          std::to_string(johnson) + ")");
    }

    const std::size_t t_agree = n - target_errors;
    const std::size_t w = k - 1;
    if (w == 0) {
      const std::size_t list_size = std::max<std::size_t>(1, n / t_agree);
      return GuruswamiSudanDecoder(code, 1, list_size);
    }

    for (std::size_t m = 1;; ++m) {
      const std::size_t d_target = m * t_agree - 1;
      const detail::Uint128 constraints =
          detail::InterpolationConstraintCount(n, m);
      const std::size_t l_max = std::max<std::size_t>(1, d_target / w);
      if (detail::WeightedMonomialCount(d_target, w, l_max) > constraints) {
        std::size_t l_lo = 1;
        std::size_t l_hi = l_max;
        while (l_lo < l_hi) {
          const std::size_t mid = l_lo + (l_hi - l_lo) / 2;
          if (detail::WeightedMonomialCount(d_target, w, mid) > constraints) {
            l_hi = mid;
          } else {
            l_lo = mid + 1;
          }
        }
        return GuruswamiSudanDecoder(code, m, l_lo);
      }
    }
  }

  // Exact Johnson list-decoding radius: n - 1 - floor(sqrt(n * (k - 1))).
  static std::size_t ComputeJohnsonRadius(std::size_t length,
                                          std::size_t dimension) {
    if (length == 0 || dimension == 0 || dimension > length) {
      throw std::invalid_argument(
          "ComputeJohnsonRadius: requires 1 <= dimension <= length");
    }
    const std::size_t root = detail::IntegerSqrtProduct(length, dimension - 1);
    return length - 1 - root;
  }

  std::size_t Multiplicity() const noexcept { return multiplicity_; }

  std::size_t ListSize() const noexcept { return list_size_; }

  std::size_t WeightedDegreeBound() const noexcept { return degree_bound_; }

  std::size_t JohnsonRadius() const noexcept { return johnson_radius_; }

  std::size_t DecodingRadius() const override { return decoding_radius_; }

  std::string Name() const override { return "GuruswamiSudan"; }

  // Returns all candidate message polynomials f(x) (deg(f) < k) whose encoded
  // codeword lies within Hamming distance <= DecodingRadius() of received_word,
  // sorted by ascending Hamming distance.
  std::vector<polynomial_type> DecodeToPolynomialList(
      const codeword_type &received_word) const {
    const auto candidates = DecodeCandidates(received_word);
    std::vector<polynomial_type> polys;
    polys.reserve(candidates.size());
    for (const auto &cand : candidates) {
      polys.push_back(cand.poly);
    }
    return polys;
  }

  // Returns all valid codewords within Hamming distance <= DecodingRadius() of
  // received_word, deduplicated and sorted by ascending Hamming distance.
  std::vector<codeword_type> DecodeToList(
      const codeword_type &received_word) const {
    const auto candidates = DecodeCandidates(received_word);
    std::vector<codeword_type> codewords;
    codewords.reserve(candidates.size());
    for (const auto &cand : candidates) {
      codewords.push_back(cand.codeword);
    }
    return codewords;
  }

  // Returns unencoded messages corresponding to DecodeToList(received_word),
  // sorted by ascending Hamming distance to received_word.
  std::vector<message_type> DecodeToMessageList(
      const codeword_type &received_word) const {
    const auto candidates = DecodeCandidates(received_word);
    std::vector<message_type> messages;
    messages.reserve(candidates.size());
    for (const auto &cand : candidates) {
      messages.push_back(grs_code_->Unencode(cand.codeword));
    }
    return messages;
  }

  // Returns the closest valid codeword within DecodingRadius(), or Failure if
  // no codeword exists in the list. Notes when multiple equidistant codewords
  // tie for minimum distance.
  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto candidates = DecodeCandidates(received_word);
    if (candidates.empty()) {
      return result_type::Failure(
          "GuruswamiSudan: no codeword found within decoding radius " +
          std::to_string(decoding_radius_));
    }
    auto best_cw = std::move(candidates[0].codeword);
    const std::size_t best_dist = candidates[0].distance;
    auto best_msg = grs_code_->Unencode(best_cw);
    auto res = result_type::Success(std::move(best_cw), std::move(best_msg),
                                    best_dist, 0);
    if (candidates.size() > 1 && candidates[1].distance == best_dist) {
      std::size_t tied = 1;
      while (tied < candidates.size() &&
             candidates[tied].distance == best_dist) {
        ++tied;
      }
      res.error_message = "GuruswamiSudan: " + std::to_string(tied) +
                          " equidistant codewords at Hamming distance " +
                          std::to_string(best_dist);
    }
    return res;
  }

 private:
  struct DecodedCandidate {
    polynomial_type poly;
    codeword_type codeword;
    std::size_t distance = 0;
  };

  void InitializeParameters(std::size_t multiplicity, std::size_t list_size) {
    const std::size_t n = grs_code_->Length();
    const std::size_t k = grs_code_->Dimension();
    multiplicity_ = multiplicity;
    list_size_ = list_size;
    johnson_radius_ = ComputeJohnsonRadius(n, k);
    degree_bound_ =
        detail::ComputeWeightedDegreeBound(n, k, multiplicity_, list_size_);
    const std::size_t min_agreements = degree_bound_ / multiplicity_ + 1;
    decoding_radius_ = (min_agreements > n) ? 0 : (n - min_agreements);
  }

  void ValidateReceivedWord(const codeword_type &received_word) const {
    const auto field = grs_code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "GuruswamiSudan: received word field mismatch");
    }
    if (received_word.Dimension() != 1 ||
        received_word.Size() != grs_code_->Length()) {
      throw std::invalid_argument(
          "GuruswamiSudan: received word must be 1D with size equal to code "
          "length");
    }
  }

  // Step 1: Koetter's iterative Groebner-basis interpolation in F_q[x, y].
  // Finds a nonzero polynomial Q(x, y) with deg_y(Q) <= list_size_ and
  // wdeg_{1, k-1}(Q) <= degree_bound_ vanishing with multiplicity >= m at every
  // normalized point (alpha_i, r_i / v_i).
  detail::BivariatePoly<storage_type> KoetterInterpolate(
      std::span<const storage_type> xs,
      std::span<const storage_type> ys) const {
    const auto field = grs_code_->Field();
    const std::size_t n = grs_code_->Length();
    const std::size_t w = grs_code_->Dimension() - 1;
    const std::size_t m = multiplicity_;
    const std::size_t L = list_size_;
    const std::size_t D = degree_bound_;
    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();

    // Build Pascal's triangle directly in F_q so C(u, a) mod p is 100% exact
    // in every characteristic (including characteristic 2 and odd extensions).
    const std::size_t max_idx = std::max(D, L);
    std::vector<std::vector<storage_type>> binom(
        max_idx + 1, std::vector<storage_type>(m, zero));
    for (std::size_t u = 0; u <= max_idx; ++u) {
      binom[u][0] = one;
      const std::size_t lim = std::min(u, m - 1);
      for (std::size_t a = 1; a <= lim; ++a) {
        if (a == u) {
          binom[u][a] = one;
        } else {
          binom[u][a] = field->Add(binom[u - 1][a - 1], binom[u - 1][a]);
        }
      }
    }

    // Initialize basis g_j(x, y) = y^j for 0 <= j <= L with wdeg(g_j) = w * j.
    std::vector<detail::BivariatePoly<storage_type>> basis(L + 1);
    std::vector<std::size_t> wdeg(L + 1, 0);
    std::vector<bool> active(L + 1, false);
    for (std::size_t j = 0; j <= L; ++j) {
      const std::size_t init_wdeg = w * j;
      if (init_wdeg <= D) {
        basis[j].y_coeffs.resize(j + 1);
        basis[j].y_coeffs[j] = {one};
        wdeg[j] = init_wdeg;
        active[j] = true;
      }
    }

    std::vector<storage_type> alpha_pow(D + 1, zero);
    std::vector<storage_type> beta_pow(L + 1, zero);
    std::vector<storage_type> x_weights(D + 1, zero);
    std::vector<storage_type> y_weights(L + 1, zero);
    std::vector<storage_type> discrepancies(L + 1, zero);

    for (std::size_t i = 0; i < n; ++i) {
      const storage_type &xi = xs[i];
      const storage_type &yi = ys[i];
      alpha_pow[0] = one;
      for (std::size_t u = 1; u <= D; ++u) {
        alpha_pow[u] = field->Mul(alpha_pow[u - 1], xi);
      }
      beta_pow[0] = one;
      for (std::size_t v = 1; v <= L; ++v) {
        beta_pow[v] = field->Mul(beta_pow[v - 1], yi);
      }

      // Process Hasse constraints in order (b = 0..m-1, a = 0..m-1-b) so that
      // (a - 1, b) is always satisfied before (a, b), guaranteeing
      // D^{[a, b]}((x - x_i) g_{j^*}) = D^{[a-1, b]}(g_{j^*}) = 0.
      for (std::size_t b = 0; b < m; ++b) {
        for (std::size_t v = b; v <= L; ++v) {
          y_weights[v] = field->Mul(binom[v][b], beta_pow[v - b]);
        }
        for (std::size_t a = 0; a + b < m; ++a) {
          for (std::size_t u = a; u <= D; ++u) {
            x_weights[u] = field->Mul(binom[u][a], alpha_pow[u - a]);
          }

          std::size_t j_star = L + 1;
          std::size_t best_wdeg = std::numeric_limits<std::size_t>::max();

          for (std::size_t j = 0; j <= L; ++j) {
            if (!active[j]) {
              continue;
            }
            storage_type disc = zero;
            const auto &y_polys = basis[j].y_coeffs;
            for (std::size_t v = b; v < y_polys.size(); ++v) {
              if (field->IsZero(y_weights[v])) {
                continue;
              }
              const auto &row = y_polys[v];
              if (row.size() <= a) {
                continue;
              }
              const std::size_t len = row.size() - a;
              const storage_type inner = field->DotProduct(
                  std::span<const storage_type>(row.data() + a, len),
                  std::span<const storage_type>(x_weights.data() + a, len));
              if (!field->IsZero(inner)) {
                disc = field->Add(disc, field->Mul(y_weights[v], inner));
              }
            }
            discrepancies[j] = disc;
            if (!field->IsZero(disc) && wdeg[j] < best_wdeg) {
              best_wdeg = wdeg[j];
              j_star = j;
            }
          }

          if (j_star > L) {
            continue;
          }

          const storage_type inv_d_star = field->Inv(discrepancies[j_star]);
          for (std::size_t j = 0; j <= L; ++j) {
            if (!active[j] || j == j_star || field->IsZero(discrepancies[j])) {
              continue;
            }
            const storage_type factor =
                field->Neg(field->Mul(discrepancies[j], inv_d_star));
            basis[j].FmaInPlace(*field, basis[j_star], factor);
          }

          if (wdeg[j_star] >= D) {
            active[j_star] = false;
            basis[j_star].y_coeffs.clear();
          } else {
            basis[j_star].MulByXMinusAlphaInPlace(*field, xi);
            ++wdeg[j_star];
          }
        }
      }
    }

    std::size_t best_j = L + 1;
    std::size_t best_wdeg = std::numeric_limits<std::size_t>::max();
    for (std::size_t j = 0; j <= L; ++j) {
      if (active[j] && !basis[j].IsZero() && wdeg[j] < best_wdeg) {
        best_wdeg = wdeg[j];
        best_j = j;
      }
    }
    if (best_j > L) {
      return detail::BivariatePoly<storage_type>{};
    }
    return std::move(basis[best_j]);
  }

  // Step 2: Roth-Ruckenstein recursive root-finding algorithm.
  // Extracts all polynomials f(x) in F_q[x] of degree < k with Q(x, f(x)) == 0.
  void RothRuckensteinRecurse(const detail::BivariatePoly<storage_type> &q,
                              std::size_t depth,
                              std::vector<storage_type> &prefix,
                              std::vector<polynomial_type> &out_polys) const {
    if (q.IsZero() || q.YDegree() <= 0) {
      return;
    }

    const auto field = grs_code_->Field();
    const std::size_t k = grs_code_->Dimension();
    const storage_type zero = field->AdditiveIdentity();

    // 1. Divide out the highest power x^r dividing Q(x, y).
    std::size_t min_x_order = std::numeric_limits<std::size_t>::max();
    for (const auto &row : q.y_coeffs) {
      if (row.empty()) {
        continue;
      }
      std::size_t ord = 0;
      while (ord < row.size() && field->IsZero(row[ord])) {
        ++ord;
      }
      if (ord < row.size() && ord < min_x_order) {
        min_x_order = ord;
        if (min_x_order == 0) {
          break;
        }
      }
    }
    if (min_x_order == std::numeric_limits<std::size_t>::max()) {
      return;
    }

    detail::BivariatePoly<storage_type> q_tilde;
    const std::size_t num_y = q.y_coeffs.size();
    q_tilde.y_coeffs.resize(num_y);
    for (std::size_t v = 0; v < num_y; ++v) {
      const auto &row = q.y_coeffs[v];
      if (row.size() > min_x_order) {
        q_tilde.y_coeffs[v].assign(
            row.begin() + static_cast<std::ptrdiff_t>(min_x_order), row.end());
      }
    }
    q_tilde.Trim(*field);
    if (q_tilde.IsZero() || q_tilde.YDegree() <= 0) {
      return;
    }

    // 2. Find all roots c in F_q of the univariate polynomial Q_tilde(0, y).
    std::vector<storage_type> h_coeffs(q_tilde.y_coeffs.size(), zero);
    for (std::size_t v = 0; v < q_tilde.y_coeffs.size(); ++v) {
      if (!q_tilde.y_coeffs[v].empty()) {
        h_coeffs[v] = q_tilde.y_coeffs[v][0];
      }
    }
    const polynomial_type h =
        polynomial_type::FromRaw(std::move(h_coeffs), field, "y");
    if (h.IsZero() || h.Degree() <= 0) {
      return;
    }
    const auto roots = utils::FindRoots(h);

    // 3. Branch on each root c.
    const std::size_t q_tilde_num_y = q_tilde.y_coeffs.size();
    for (const auto &root_elem : roots) {
      const storage_type c = root_elem.RawValue();
      prefix.push_back(c);

      if (depth + 1 == k) {
        // Leaf: f(x) is completely determined; verify Q_tilde(x, c) == 0.
        std::size_t max_x_len = 0;
        for (const auto &row : q_tilde.y_coeffs) {
          if (row.size() > max_x_len) {
            max_x_len = row.size();
          }
        }
        bool vanishes = true;
        for (std::size_t u = 0; u < max_x_len; ++u) {
          storage_type acc = zero;
          for (std::size_t v = q_tilde_num_y; v > 0; --v) {
            const std::size_t idx = v - 1;
            acc = field->Mul(acc, c);
            if (u < q_tilde.y_coeffs[idx].size()) {
              acc = field->Add(acc, q_tilde.y_coeffs[idx][u]);
            }
          }
          if (!field->IsZero(acc)) {
            vanishes = false;
            break;
          }
        }
        if (vanishes) {
          out_polys.push_back(polynomial_type::FromRaw(prefix, field));
        }
      } else {
        // Compute Q_next(x, y) = Q_tilde(x, x * y + c).
        detail::BivariatePoly<storage_type> q_next;
        if (field->IsZero(c)) {
          q_next.y_coeffs.resize(q_tilde_num_y);
          for (std::size_t v = 0; v < q_tilde_num_y; ++v) {
            if (!q_tilde.y_coeffs[v].empty()) {
              q_next.y_coeffs[v].assign(v, zero);
              q_next.y_coeffs[v].insert(q_next.y_coeffs[v].end(),
                                        q_tilde.y_coeffs[v].begin(),
                                        q_tilde.y_coeffs[v].end());
            }
          }
        } else {
          // First compute R(x, y) = Q_tilde(x, y + c) via Horner's rule in y.
          std::vector<std::vector<storage_type>> r_polys;
          r_polys.reserve(q_tilde_num_y);
          r_polys.push_back(q_tilde.y_coeffs[q_tilde_num_y - 1]);
          for (std::size_t step = q_tilde_num_y - 1; step > 0; --step) {
            const std::size_t v = step - 1;
            const std::size_t cur_size = r_polys.size();
            r_polys.push_back(r_polys[cur_size - 1]);
            for (std::size_t j = cur_size - 1; j > 0; --j) {
              auto &dst = r_polys[j];
              const auto &prev = r_polys[j - 1];
              if (!dst.empty()) {
                field->ScalarMulArray(std::span<storage_type>(dst),
                                      std::span<const storage_type>(dst), c);
              }
              if (!prev.empty()) {
                if (dst.size() < prev.size()) {
                  dst.resize(prev.size(), zero);
                }
                field->AddArray(
                    std::span<storage_type>(dst.data(), prev.size()),
                    std::span<const storage_type>(dst.data(), prev.size()),
                    std::span<const storage_type>(prev));
              }
            }
            {
              auto &dst = r_polys[0];
              const auto &qv = q_tilde.y_coeffs[v];
              if (!dst.empty()) {
                field->ScalarMulArray(std::span<storage_type>(dst),
                                      std::span<const storage_type>(dst), c);
              }
              if (!qv.empty()) {
                if (dst.size() < qv.size()) {
                  dst.resize(qv.size(), zero);
                }
                field->AddArray(
                    std::span<storage_type>(dst.data(), qv.size()),
                    std::span<const storage_type>(dst.data(), qv.size()),
                    std::span<const storage_type>(qv));
              }
            }
          }
          for (auto &row : r_polys) {
            while (!row.empty() && field->IsZero(row.back())) {
              row.pop_back();
            }
          }
          // Substitute y -> x * y: multiply row j by x^j.
          q_next.y_coeffs.resize(r_polys.size());
          for (std::size_t j = 0; j < r_polys.size(); ++j) {
            if (!r_polys[j].empty()) {
              q_next.y_coeffs[j].assign(j, zero);
              q_next.y_coeffs[j].insert(q_next.y_coeffs[j].end(),
                                        r_polys[j].begin(), r_polys[j].end());
            }
          }
        }
        q_next.Trim(*field);

        // If y divides Q_next(x, y), then Q_next(x, 0) == 0, meaning the
        // current prefix (zero-padded to degree < k) is already a root. Record
        // it, divide out all factors of y, and continue searching if higher
        // y-degree factors remain.
        if (!q_next.IsZero() && q_next.y_coeffs[0].empty()) {
          out_polys.push_back(polynomial_type::FromRaw(prefix, field));
          std::size_t y_shift = 0;
          while (y_shift < q_next.y_coeffs.size() &&
                 q_next.y_coeffs[y_shift].empty()) {
            ++y_shift;
          }
          q_next.y_coeffs.erase(
              q_next.y_coeffs.begin(),
              q_next.y_coeffs.begin() + static_cast<std::ptrdiff_t>(y_shift));
        }

        if (q_next.YDegree() > 0) {
          RothRuckensteinRecurse(q_next, depth + 1, prefix, out_polys);
        }
      }

      prefix.pop_back();
    }
  }

  std::vector<DecodedCandidate> DecodeCandidates(
      const codeword_type &received_word) const {
    ValidateReceivedWord(received_word);

    const auto field = grs_code_->Field();
    const std::size_t n = grs_code_->Length();
    const std::size_t k = grs_code_->Dimension();

    // When configured parameters yield DecodingRadius() == 0, only an already
    // valid codeword at Hamming distance 0 can be accepted.
    if (decoding_radius_ == 0) {
      if (grs_code_->Contains(received_word)) {
        polynomial_type poly = grs_code_->InterpolateCodeword(received_word);
        return {DecodedCandidate{std::move(poly), received_word, 0}};
      }
      return {};
    }

    // Normalize y_i = r_i / v_i using precomputed column multiplier inverses.
    const auto eval_span = grs_code_->EvaluationPoints().RawSpan();
    const auto rx_span = received_word.RawSpan();
    const auto &col_mult_invs = grs_code_->ColumnMultiplierInverses();
    std::vector<storage_type> ys(n);
    field->MulArray(std::span<storage_type>(ys), rx_span,
                    std::span<const storage_type>(col_mult_invs));

    // Step 1: Koetter bivariate interpolation.
    const detail::BivariatePoly<storage_type> q =
        KoetterInterpolate(eval_span, std::span<const storage_type>(ys));
    if (q.IsZero() || q.YDegree() <= 0) {
      if (grs_code_->Contains(received_word)) {
        polynomial_type poly = grs_code_->InterpolateCodeword(received_word);
        return {DecodedCandidate{std::move(poly), received_word, 0}};
      }
      return {};
    }

    // Step 2: Roth-Ruckenstein root-finding.
    std::vector<storage_type> prefix;
    prefix.reserve(k);
    std::vector<polynomial_type> raw_polys;
    RothRuckensteinRecurse(q, 0, prefix, raw_polys);

    // Deduplicate candidate polynomials canonically.
    std::sort(raw_polys.begin(), raw_polys.end(),
              detail::CanonicalPolyLess<GaloisField>);
    raw_polys.erase(std::unique(raw_polys.begin(), raw_polys.end()),
                    raw_polys.end());

    // Re-encode each candidate and retain every codeword within
    // DecodingRadius() of received_word.
    std::vector<DecodedCandidate> candidates;
    candidates.reserve(raw_polys.size());
    for (auto &poly : raw_polys) {
      codeword_type cw = grs_code_->EncodePolynomial(poly);
      const auto cw_span = cw.RawSpan();
      std::size_t dist = 0;
      for (std::size_t j = 0; j < n; ++j) {
        if (cw_span[j] != rx_span[j]) {
          ++dist;
        }
      }
      if (dist <= decoding_radius_) {
        candidates.push_back(
            DecodedCandidate{std::move(poly), std::move(cw), dist});
      }
    }

    // Sort by ascending Hamming distance (canonical tie-breaking).
    std::sort(candidates.begin(), candidates.end(),
              [](const DecodedCandidate &a, const DecodedCandidate &b) {
                if (a.distance != b.distance) {
                  return a.distance < b.distance;
                }
                return detail::CanonicalPolyLess(a.poly, b.poly);
              });
    return candidates;
  }

  const GeneralizedReedSolomonCode<GaloisField> *grs_code_;
  std::size_t multiplicity_ = 1;
  std::size_t list_size_ = 1;
  std::size_t degree_bound_ = 0;
  std::size_t johnson_radius_ = 0;
  std::size_t decoding_radius_ = 0;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_GURUSWAMI_SUDAN_DECODER_HPP_
