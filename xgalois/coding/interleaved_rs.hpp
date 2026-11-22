#ifndef XGALOIS_CODING_INTERLEAVED_RS_HPP_
#define XGALOIS_CODING_INTERLEAVED_RS_HPP_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/grs_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/linear_code.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace coding {

// Collaborative decoding algorithm selection for Interleaved GRS codes and
// Power Decoding of single GRS codes.
enum class IRSDecoderMode : uint8_t {
  kSsb,          // Fast O(ell * r * t) SSB 2011 with kSsbExact escalation
  kSsbExact,     // Exact O(ell^2 * r * t) Koetter-Nielsen Popov minimization
  kBkyHankel,    // O(ell * r * t^2) Bleichenbacher-Kiayias-Yung stacked Hankel
  kIndependent,  // Per-column independent Berlekamp-Massey baseline
};

// Structured result returned by InterleavedGRSCode::Decode and
// InterleavedGRSCode::DecodeWithErasures.
template <typename GaloisField>
struct InterleavedDecodingResult {
  using matrix_type = GFMatrix<GaloisField>;
  using vector_type = GFArray<GaloisField>;

  matrix_type codeword_matrix;                     // shape {n, ell}
  std::vector<vector_type> messages;               // ell vectors of size k_m
  std::vector<std::size_t> burst_error_positions;  // sorted in [0, n)
  std::size_t num_burst_errors = 0;
  std::size_t num_erasures = 0;
  double failure_probability_bound = 0.0;
  bool used_fallback = false;
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }

  static InterleavedDecodingResult Failure(std::string reason) {
    InterleavedDecodingResult res;
    res.success = false;
    res.error_message = std::move(reason);
    return res;
  }
};

namespace detail {

// Computes the collaborative burst-error decoding radius:
//   t_max = max { t >= 0 : sum_{m : r_m > t} (r_m - t) >= t }.
// For homogeneous r_m = r across ell streams, this equals
// floor(ell / (ell + 1) * r).
inline std::size_t ComputeCollaborativeRadius(
    const std::vector<std::size_t> &redundancies) noexcept {
  std::size_t max_r = 0;
  for (std::size_t r : redundancies) {
    max_r = std::max(max_r, r);
  }
  std::size_t best_t = 0;
  for (std::size_t t = 1; t <= max_r; ++t) {
    std::size_t eq_count = 0;
    for (std::size_t r : redundancies) {
      if (r > t) {
        eq_count += (r - t);
      }
    }
    if (eq_count >= t) {
      best_t = t;
    } else {
      break;
    }
  }
  return best_t;
}

// Checks whether all syndrome sequences are identically zero.
template <typename GaloisField>
inline bool AreAllSequencesZero(
    const GaloisField &field,
    const std::vector<std::vector<typename GaloisField::element_type>> &seqs) {
  for (const auto &seq : seqs) {
    for (const auto &v : seq) {
      if (!field.IsZero(v)) {
        return false;
      }
    }
  }
  return true;
}

// Output of fast Schmidt-Sidorenko-Bossert (2011) multi-sequence LFSR
// synthesis.
template <typename GaloisField>
struct SsbLfsrOutput {
  PolynomialDense<GaloisField> locator;
  std::size_t length = 0;
  bool active_streams_sufficient = false;
};

// Schmidt-Sidorenko-Bossert (2011) fast O(ell * r * t) multi-sequence
// shift-register synthesis with exact per-stream backup records
// (B^(m), L_m, d_m, i_m, valid_m) and the exact early-break invariant when
// a stream triggers !valid[m].
template <typename GaloisField>
inline SsbLfsrOutput<GaloisField> SolveMultiSequenceSsb(
    const std::shared_ptr<GaloisField> &field,
    const std::vector<std::vector<typename GaloisField::element_type>> &seqs,
    std::size_t max_errors) {
  using StorageType = typename GaloisField::element_type;
  using Poly = PolynomialDense<GaloisField>;

  const std::size_t ell = seqs.size();
  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();

  if (AreAllSequencesZero(*field, seqs)) {
    return SsbLfsrOutput<GaloisField>{Poly::One(field), 0, true};
  }

  std::size_t max_len = 0;
  for (const auto &seq : seqs) {
    max_len = std::max(max_len, seq.size());
  }

  std::vector<StorageType> lambda = {one};
  lambda.reserve(max_errors + 2);
  std::size_t L = 0;

  std::vector<std::vector<StorageType>> backup_poly(ell);
  std::vector<std::size_t> backup_L(ell, 0);
  std::vector<StorageType> backup_d(ell, zero);
  std::vector<std::size_t> backup_i(ell, 0);
  std::vector<bool> valid(ell, false);

  for (std::size_t i = 0; i < max_len; ++i) {
    for (std::size_t m = 0; m < ell; ++m) {
      if (i >= seqs[m].size()) {
        continue;
      }
      if (i < L) {
        break;
      }

      const auto &seq_m = seqs[m];
      StorageType d = seq_m[i];
      const std::size_t lim = std::min(L, lambda.size() - 1);
      for (std::size_t u = 1; u <= lim; ++u) {
        if (!field->IsZero(lambda[u]) && !field->IsZero(seq_m[i - u])) {
          d = field->Add(d, field->Mul(lambda[u], seq_m[i - u]));
        }
      }

      if (field->IsZero(d)) {
        continue;
      }

      if (!valid[m]) {
        // First non-zero discrepancy on stream m: record backup and bump L to
        // i + 1. Because i < L now holds for every subsequent stream m' > m at
        // step i, we MUST break out of the inner loop immediately.
        backup_poly[m] = lambda;
        backup_L[m] = L;
        backup_d[m] = d;
        backup_i[m] = i;
        valid[m] = true;
        L = i + 1;
        break;
      }

      const std::size_t shift = i - backup_i[m];
      const std::size_t cand_L = backup_L[m] + shift;
      const StorageType neg_factor = field->Neg(field->Div(d, backup_d[m]));
      const bool grow_L = (cand_L > L);

      std::vector<StorageType> old_lambda;
      const std::size_t old_L = L;
      if (grow_L) {
        old_lambda = lambda;
      }

      const std::size_t needed = backup_poly[m].size() + shift;
      if (lambda.size() < needed) {
        lambda.resize(needed, zero);
      }
      field->FmaArray(
          std::span<StorageType>(lambda.data() + shift, backup_poly[m].size()),
          std::span<const StorageType>(backup_poly[m]), neg_factor);

      if (grow_L) {
        backup_poly[m] = std::move(old_lambda);
        backup_L[m] = old_L;
        backup_d[m] = d;
        backup_i[m] = i;
        L = cand_L;
      }
    }
  }

  // Check whether the streams that actually contributed non-zero discrepancies
  // supply enough linear equations to uniquely pin down L unknowns.
  std::size_t active_eqs = 0;
  for (std::size_t m = 0; m < ell; ++m) {
    if (valid[m] && seqs[m].size() > L) {
      active_eqs += (seqs[m].size() - L);
    }
  }

  Poly poly = Poly::FromRaw(std::move(lambda), field);
  return SsbLfsrOutput<GaloisField>{std::move(poly), L, active_eqs >= L};
}

// Koetter-Nielsen Popov / Demand-Driven Module Minimization in O(ell^2 * r * t)
// over F_q[z]^(ell+1) with shift s = (0, 1, ..., 1).
//
// Maintains (ell + 1) polynomial rows forming a weak Popov basis whose leading
// positions under smallest-column tie-breaking stay fixed at LP'(g_k) = k.
// Consequently, row 0 ALWAYS holds the minimal-degree Lambda(z) satisfying all
// syndrome recurrences across all ell streams, and uniqueness holds iff no
// other row k >= 1 has both w_k <= w_0 and Lambda_k != 0.
template <typename GaloisField>
inline std::optional<PolynomialDense<GaloisField>> SolveMultiSequencePopov(
    const std::shared_ptr<GaloisField> &field,
    const std::vector<std::vector<typename GaloisField::element_type>> &seqs,
    std::size_t max_errors) {
  using StorageType = typename GaloisField::element_type;
  using Poly = PolynomialDense<GaloisField>;

  const std::size_t ell = seqs.size();
  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();

  if (AreAllSequencesZero(*field, seqs)) {
    return Poly::One(field);
  }
  if (max_errors == 0) {
    return std::nullopt;
  }

  const std::size_t num_rows = ell + 1;

  // lambda_rows[k] holds column 0 of basis row k; weights[k] holds deg_s(g_k).
  // Row 0 starts as (1, 0, ..., 0) with weight 0.
  // Row m+1 starts as (0, ..., 1, ..., 0) with weight 1.
  std::vector<std::vector<StorageType>> lambda_rows(num_rows);
  lambda_rows[0] = {one};
  std::vector<std::size_t> weights(num_rows, 1);
  weights[0] = 0;

  // Because rows m' + 1 (m' > m) have Lambda = 0 and Omega^(m) = 0, their
  // discrepancy on stream m is identically zero, so we only track the current
  // stream's temporary omega polynomial across rows 0 .. m + 1.
  std::vector<std::vector<StorageType>> omega_cur(num_rows);
  std::vector<StorageType> disc(num_rows, zero);

  for (std::size_t m = 0; m < ell; ++m) {
    const auto &seq_m = seqs[m];
    const std::size_t M_m = seq_m.size();
    if (M_m == 0) {
      continue;
    }

    const std::size_t active_limit = m + 2;  // rows 0 .. m + 1
    for (std::size_t k = 0; k < active_limit; ++k) {
      omega_cur[k].clear();
    }
    omega_cur[m + 1].push_back(one);

    for (std::size_t i = 0; i < M_m; ++i) {
      // Find k* among active rows with weights[k] <= max_errors and non-zero
      // discrepancy, minimizing weight and breaking ties by LARGEST k so that
      // leading positions LP'(g_k) = k (smallest column achieving deg_s) are
      // strictly preserved.
      std::size_t best_k = num_rows;
      std::size_t best_w = max_errors + 1;

      for (std::size_t k = 0; k < active_limit; ++k) {
        if (weights[k] > max_errors) {
          disc[k] = zero;
          continue;
        }

        StorageType d = zero;
        const auto &lam = lambda_rows[k];
        if (!lam.empty()) {
          const std::size_t u_max = std::min(i, lam.size() - 1);
          for (std::size_t u = 0; u <= u_max; ++u) {
            if (!field->IsZero(lam[u]) && !field->IsZero(seq_m[i - u])) {
              d = field->Add(d, field->Mul(lam[u], seq_m[i - u]));
            }
          }
        }
        if (i < omega_cur[k].size() && !field->IsZero(omega_cur[k][i])) {
          d = field->Sub(d, omega_cur[k][i]);
        }

        disc[k] = d;
        if (!field->IsZero(d)) {
          if (best_k == num_rows || weights[k] < best_w ||
              (weights[k] == best_w && k > best_k)) {
            best_k = k;
            best_w = weights[k];
          }
        }
      }

      if (best_k == num_rows) {
        continue;
      }

      const StorageType inv_best = field->Inv(disc[best_k]);
      const auto &lam_best = lambda_rows[best_k];
      const auto &om_best = omega_cur[best_k];

      for (std::size_t k = 0; k < active_limit; ++k) {
        if (k == best_k || weights[k] > max_errors || field->IsZero(disc[k])) {
          continue;
        }
        const StorageType neg_factor =
            field->Neg(field->Mul(disc[k], inv_best));

        if (!lam_best.empty()) {
          if (lambda_rows[k].size() < lam_best.size()) {
            lambda_rows[k].resize(lam_best.size(), zero);
          }
          field->FmaArray(
              std::span<StorageType>(lambda_rows[k].data(), lam_best.size()),
              std::span<const StorageType>(lam_best), neg_factor);
        }
        if (!om_best.empty()) {
          if (omega_cur[k].size() < om_best.size()) {
            omega_cur[k].resize(om_best.size(), zero);
          }
          field->FmaArray(
              std::span<StorageType>(omega_cur[k].data(), om_best.size()),
              std::span<const StorageType>(om_best), neg_factor);
        }
      }

      // Multiply pivot row best_k by z and increment its shifted degree.
      if (!lambda_rows[best_k].empty()) {
        lambda_rows[best_k].insert(lambda_rows[best_k].begin(), zero);
      }
      if (!omega_cur[best_k].empty()) {
        omega_cur[best_k].insert(omega_cur[best_k].begin(), zero);
      }
      ++weights[best_k];
    }
  }

  const std::size_t deg0 = weights[0];
  if (deg0 > max_errors) {
    return std::nullopt;
  }

  // Exact uniqueness check via the Predictable Degree Property of weak Popov
  // bases: Lambda_0 is the unique minimal locator of degree <= deg0 iff no
  // other row k >= 1 has w_k <= deg0 with a non-zero Lambda_k.
  for (std::size_t k = 1; k < num_rows; ++k) {
    if (weights[k] <= deg0) {
      bool non_zero_lambda = false;
      for (const auto &c : lambda_rows[k]) {
        if (!field->IsZero(c)) {
          non_zero_lambda = true;
          break;
        }
      }
      if (non_zero_lambda) {
        return std::nullopt;
      }
    }
  }

  Poly lambda = Poly::FromRaw(std::move(lambda_rows[0]), field);
  if (lambda.IsZero() || field->IsZero(lambda.RawCoeff(0))) {
    return std::nullopt;
  }
  if (lambda.Degree() != static_cast<int>(deg0)) {
    return std::nullopt;
  }
  const StorageType l0 = lambda.RawCoeff(0);
  if (!field->IsOne(l0)) {
    const auto inv_l0 =
        GaloisFieldElement<GaloisField>::FromInternal(field->Inv(l0), field);
    lambda *= inv_l0;
  }
  return lambda;
}

// Bleichenbacher-Kiayias-Yung (BKY) stacked Hankel key-equation solver in
// O(ell * r * t^2): searches descending from tau = max_errors down to 1 on
// the stacked Hankel system M_tau * lambda = -b_tau via linalg::RrefWithPivots.
template <typename GaloisField>
inline std::optional<PolynomialDense<GaloisField>> SolveMultiSequenceBkyHankel(
    const std::shared_ptr<GaloisField> &field,
    const std::vector<std::vector<typename GaloisField::element_type>> &seqs,
    std::size_t max_errors) {
  using StorageType = typename GaloisField::element_type;
  using Poly = PolynomialDense<GaloisField>;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  if (AreAllSequencesZero(*field, seqs)) {
    return Poly::One(field);
  }
  if (max_errors == 0) {
    return std::nullopt;
  }

  const std::size_t ell = seqs.size();
  const StorageType one = field->MultiplicativeIdentity();

  for (std::size_t tau = max_errors; tau >= 1; --tau) {
    std::size_t num_rows = 0;
    for (std::size_t m = 0; m < ell; ++m) {
      if (seqs[m].size() > tau) {
        num_rows += (seqs[m].size() - tau);
      }
    }
    if (num_rows < tau) {
      continue;
    }

    // Augmented matrix [M_tau | -b_tau] of shape {num_rows, tau + 1} where
    // unknown c in [0, tau) directly represents Lambda_{c + 1}:
    //   sum_{c=0}^{tau-1} seq[i + tau - 1 - c] * Lambda_{c + 1} = -seq[i + tau]
    RawArray aug_raw = RawArray::from_shape(ShapeType{num_rows, tau + 1});
    StorageType *aug_data = aug_raw.data();
    const std::size_t stride = tau + 1;
    std::size_t row = 0;
    for (std::size_t m = 0; m < ell; ++m) {
      const auto &seq_m = seqs[m];
      if (seq_m.size() <= tau) {
        continue;
      }
      const std::size_t rows_m = seq_m.size() - tau;
      for (std::size_t i = 0; i < rows_m; ++i) {
        StorageType *row_ptr = aug_data + row * stride;
        for (std::size_t c = 0; c < tau; ++c) {
          row_ptr[c] = seq_m[i + tau - 1 - c];
        }
        row_ptr[tau] = field->Neg(seq_m[i + tau]);
        ++row;
      }
    }

    auto aug = GFArray<GaloisField>::FromRaw(field, std::move(aug_raw));
    auto [rref_mat, pivots] = linalg::detail::RrefWithPivots(aug);

    // Consistent & unique iff columns 0..tau-1 are all pivots and the RHS
    // column tau is NOT a pivot.
    if (pivots.size() == tau && pivots.back() < tau) {
      std::vector<StorageType> coeffs(tau + 1);
      coeffs[0] = one;
      for (std::size_t c = 0; c < tau; ++c) {
        coeffs[c + 1] = rref_mat.Raw(c, tau);
      }
      Poly lambda = Poly::FromRaw(std::move(coeffs), field);
      if (lambda.Degree() == static_cast<int>(tau)) {
        return lambda;
      }
    }
  }

  return std::nullopt;
}

// Builds standard primitive RS evaluation points (alpha^0, ..., alpha^(n-1)).
template <typename GaloisField>
inline GFArray<GaloisField> BuildPrimitiveEvalPoints(
    const std::shared_ptr<GaloisField> &field, std::size_t length) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  if (!field) {
    throw std::invalid_argument("InterleavedGRSCode: field cannot be null");
  }
  if (length == 0) {
    throw std::invalid_argument("InterleavedGRSCode: length must be positive");
  }
  const uint64_t q = field->Order();
  if (q != 0 && length > static_cast<std::size_t>(q - 1)) {
    throw std::invalid_argument(
        "InterleavedGRSCode: primitive RS length cannot exceed q - 1");
  }

  RawArray raw = RawArray::from_shape(ShapeType{length});
  const StorageType alpha = field->MultiplicativeGenerator();
  StorageType power = field->MultiplicativeIdentity();
  for (std::size_t i = 0; i < length; ++i) {
    raw(i) = power;
    if (i + 1 < length) {
      power = field->Mul(power, alpha);
    }
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
}

// Builds extended RS evaluation points (0, alpha^0, ..., alpha^(n-2)).
template <typename GaloisField>
inline GFArray<GaloisField> BuildExtendedEvalPoints(
    const std::shared_ptr<GaloisField> &field, std::size_t length) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;

  if (!field) {
    throw std::invalid_argument("InterleavedGRSCode: field cannot be null");
  }
  if (length == 0) {
    throw std::invalid_argument("InterleavedGRSCode: length must be positive");
  }
  const uint64_t q = field->Order();
  if (q != 0 && length > static_cast<std::size_t>(q)) {
    throw std::invalid_argument(
        "InterleavedGRSCode: extended RS length cannot exceed q");
  }

  RawArray raw = RawArray::from_shape(ShapeType{length});
  raw(0) = field->AdditiveIdentity();
  const StorageType alpha = field->MultiplicativeGenerator();
  StorageType power = field->MultiplicativeIdentity();
  for (std::size_t i = 1; i < length; ++i) {
    raw(i) = power;
    if (i + 1 < length) {
      power = field->Mul(power, alpha);
    }
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
}

}  // namespace detail

template <typename GaloisField>
class PowerDecodingDecoder;

// Homogeneous IRS(n, k; ell) and Heterogeneous IRS(n; k_0, ..., k_{ell-1})
// Interleaved Generalized Reed-Solomon codes over shared distinct evaluation
// points alpha (including alpha_{j_0} = 0) and non-zero column multipliers v.
template <typename GaloisField>
class InterleavedGRSCode {
 public:
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using vector_type = GFArray<GaloisField>;
  using matrix_type = GFMatrix<GaloisField>;
  using codeword_type = matrix_type;
  using message_type = std::vector<vector_type>;
  using polynomial_type = PolynomialDense<GaloisField>;
  using component_code_type = GeneralizedReedSolomonCode<GaloisField>;
  using result_type = InterleavedDecodingResult<GaloisField>;

  InterleavedGRSCode(std::shared_ptr<GaloisField> field,
                     const vector_type &evaluation_points,
                     const vector_type &column_multipliers,
                     std::vector<std::size_t> dimensions)
      : field_(std::move(field)), dimensions_(std::move(dimensions)) {
    if (!field_) {
      throw std::invalid_argument("InterleavedGRSCode: field cannot be null");
    }
    if (dimensions_.empty()) {
      throw std::invalid_argument(
          "InterleavedGRSCode: interleaving order must be at least 1");
    }
    length_ = evaluation_points.Size();
    interleaving_order_ = dimensions_.size();

    component_codes_.reserve(interleaving_order_);
    redundancies_.reserve(interleaving_order_);
    is_homogeneous_ = true;
    const std::size_t k0 = dimensions_[0];

    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      if (dimensions_[m] != k0) {
        is_homogeneous_ = false;
      }
      component_codes_.emplace_back(field_, evaluation_points,
                                    column_multipliers, dimensions_[m]);
      redundancies_.push_back(length_ - dimensions_[m]);
    }

    min_redundancy_ =
        *std::min_element(redundancies_.begin(), redundancies_.end());
    max_redundancy_ =
        *std::max_element(redundancies_.begin(), redundancies_.end());
    unique_radius_ = min_redundancy_ / 2;
    collaborative_radius_ = detail::ComputeCollaborativeRadius(redundancies_);
    max_redundancy_idx_ = static_cast<std::size_t>(std::distance(
        redundancies_.begin(),
        std::max_element(redundancies_.begin(), redundancies_.end())));
  }

  InterleavedGRSCode(std::shared_ptr<GaloisField> field,
                     const vector_type &evaluation_points,
                     std::vector<std::size_t> dimensions)
      : InterleavedGRSCode(field, evaluation_points,
                           vector_type::Ones(field, {evaluation_points.Size()}),
                           std::move(dimensions)) {}

  // Homogeneous IRS(n, k; ell) over explicit evaluation points and multipliers.
  static InterleavedGRSCode CreateHomogeneous(
      std::shared_ptr<GaloisField> field, const vector_type &evaluation_points,
      const vector_type &column_multipliers, std::size_t dimension,
      std::size_t interleaving_order) {
    if (interleaving_order == 0) {
      throw std::invalid_argument(
          "CreateHomogeneous: interleaving order must be at least 1");
    }
    return InterleavedGRSCode(
        std::move(field), evaluation_points, column_multipliers,
        std::vector<std::size_t>(interleaving_order, dimension));
  }

  // Homogeneous standard IRS(n, k; ell) over explicit evaluation points with
  // unit column multipliers.
  static InterleavedGRSCode CreateHomogeneous(
      std::shared_ptr<GaloisField> field, const vector_type &evaluation_points,
      std::size_t dimension, std::size_t interleaving_order) {
    if (!field) {
      throw std::invalid_argument("CreateHomogeneous: field cannot be null");
    }
    auto ones = vector_type::Ones(field, {evaluation_points.Size()});
    return CreateHomogeneous(std::move(field), evaluation_points, ones,
                             dimension, interleaving_order);
  }

  // Homogeneous primitive IRS(n, k; ell) using (1, alpha, ..., alpha^(n-1)).
  static InterleavedGRSCode CreateHomogeneous(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension, std::size_t interleaving_order) {
    auto evals = detail::BuildPrimitiveEvalPoints(field, length);
    auto ones = vector_type::Ones(field, {length});
    return CreateHomogeneous(std::move(field), evals, ones, dimension,
                             interleaving_order);
  }

  // Homogeneous extended IRS(n, k; ell) using (0, 1, alpha, ..., alpha^(n-2)).
  static InterleavedGRSCode CreateExtendedHomogeneous(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension, std::size_t interleaving_order) {
    auto evals = detail::BuildExtendedEvalPoints(field, length);
    auto ones = vector_type::Ones(field, {length});
    return CreateHomogeneous(std::move(field), evals, ones, dimension,
                             interleaving_order);
  }

  // Heterogeneous IRS(n; k_0, ..., k_{ell-1}) over explicit evaluation points
  // and column multipliers.
  static InterleavedGRSCode CreateHeterogeneous(
      std::shared_ptr<GaloisField> field, const vector_type &evaluation_points,
      const vector_type &column_multipliers,
      std::vector<std::size_t> dimensions) {
    return InterleavedGRSCode(std::move(field), evaluation_points,
                              column_multipliers, std::move(dimensions));
  }

  // Heterogeneous primitive IRS(n; k_0, ..., k_{ell-1}).
  static InterleavedGRSCode CreateHeterogeneous(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::vector<std::size_t> dimensions) {
    auto evals = detail::BuildPrimitiveEvalPoints(field, length);
    auto ones = vector_type::Ones(field, {length});
    return InterleavedGRSCode(std::move(field), evals, ones,
                              std::move(dimensions));
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }
  std::size_t Length() const noexcept { return length_; }
  std::size_t InterleavingOrder() const noexcept { return interleaving_order_; }
  const std::vector<std::size_t> &Dimensions() const noexcept {
    return dimensions_;
  }
  const std::vector<std::size_t> &Redundancies() const noexcept {
    return redundancies_;
  }
  bool IsHomogeneous() const noexcept { return is_homogeneous_; }

  std::size_t CollaborativeDecodingRadius() const noexcept {
    return collaborative_radius_;
  }
  std::size_t UniqueDecodingRadius() const noexcept { return unique_radius_; }

  const vector_type &EvaluationPoints() const noexcept {
    return component_codes_[0].EvaluationPoints();
  }
  const vector_type &ColumnMultipliers() const noexcept {
    return component_codes_[0].ColumnMultipliers();
  }
  const vector_type &DualMultipliers() const noexcept {
    return component_codes_[0].DualMultipliers();
  }
  std::optional<std::size_t> ZeroEvaluationIndex() const noexcept {
    return component_codes_[0].ZeroEvaluationIndex();
  }

  const component_code_type &ComponentCode(std::size_t m) const {
    if (m >= interleaving_order_) {
      throw std::out_of_range(
          "InterleavedGRSCode::ComponentCode: stream index out of range");
    }
    return component_codes_[m];
  }

  // Exact Bleichenbacher-Kiayias-Yung / Schmidt-Sidorenko-Bossert (Theorem 5)
  // union upper bound on collaborative decoding failure probability for
  // num_burst_errors uniformly random non-zero burst error columns in
  // F_q^ell \ {0}.
  //
  // Computed in log-space so large field orders (including GF(2^64) where
  // field->Order() == 0) and large t never overflow double.
  double FailureProbabilityBound(std::size_t num_burst_errors) const {
    const std::size_t t = num_burst_errors;
    if (t <= unique_radius_) {
      return 0.0;
    }
    if (t > collaborative_radius_) {
      return 1.0;
    }

    const uint64_t q_raw = field_->Order();
    const double log_q =
        (q_raw == 0) ? (static_cast<double>(field_->Degree()) *
                        std::log(static_cast<double>(field_->Characteristic())))
                     : std::log(static_cast<double>(q_raw));

    // log(q - 1) and log(q^ell / (q^ell - 1)) = -log1p(-q^(-ell)).
    const double log_qm1 =
        (q_raw == 0) ? log_q : std::log(static_cast<double>(q_raw - 1));
    const double ell_log_q = static_cast<double>(interleaving_order_) * log_q;
    const double log_ratio =
        (ell_log_q > 50.0) ? 0.0 : -std::log1p(-std::exp(-ell_log_q));

    std::size_t d_min = length_;
    std::vector<std::size_t> d_m(interleaving_order_, 0);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      d_m[m] = (redundancies_[m] > t) ? (redundancies_[m] - t) : 0;
      d_min = std::min(d_min, d_m[m]);
    }

    // For a projective nullspace vector u in P^(t-1)(F_q) of Hamming weight w:
    //   - if w <= d_min, A_{d_m} D_m u = 0 forces e_a = 0 (impossible).
    //   - if w > d_min, probability <= (q^ell / (q^ell - 1))^w * q^(-sum min(w,
    //   d_m)).
    // Weight-w projective count is C(t, w) * (q - 1)^(w - 1).
    double total = 0.0;
    const double lg_t_plus_1 = std::lgamma(static_cast<double>(t + 1));
    for (std::size_t w = d_min + 1; w <= t; ++w) {
      std::size_t constraint_exp = 0;
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        constraint_exp += std::min(w, d_m[m]);
      }
      const double log_binom = lg_t_plus_1 -
                               std::lgamma(static_cast<double>(w + 1)) -
                               std::lgamma(static_cast<double>(t - w + 1));
      const double log_term = log_binom + static_cast<double>(w - 1) * log_qm1 +
                              static_cast<double>(w) * log_ratio -
                              static_cast<double>(constraint_exp) * log_q;
      if (log_term >= 0.0) {
        return 1.0;
      }
      total += std::exp(log_term);
      if (total >= 1.0) {
        return 1.0;
      }
    }
    return std::min(1.0, total);
  }

  // Encodes ell message vectors (messages[m] of size k_m) into an {n, ell}
  // codeword matrix whose m-th column belongs to ComponentCode(m).
  matrix_type Encode(const std::vector<vector_type> &messages) const {
    if (messages.size() != interleaving_order_) {
      throw std::invalid_argument("InterleavedGRSCode::Encode: expected " +
                                  std::to_string(interleaving_order_) +
                                  " messages");
    }
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray cw_raw =
        RawArray::from_shape(ShapeType{length_, interleaving_order_});
    storage_type *cw_data = cw_raw.data();

    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const auto col_cw = component_codes_[m].Encode(messages[m]);
      const auto col_span = col_cw.RawSpan();
      for (std::size_t j = 0; j < length_; ++j) {
        cw_data[j * interleaving_order_ + m] = col_span[j];
      }
    }
    return matrix_type::FromRaw(field_, std::move(cw_raw));
  }

  // Encodes a {k, ell} message matrix into an {n, ell} codeword matrix for
  // homogeneous IRS codes using C = G^T * U.
  matrix_type EncodeHomogeneous(const matrix_type &message_matrix) const {
    if (!is_homogeneous_) {
      throw std::invalid_argument(
          "InterleavedGRSCode::EncodeHomogeneous: code is not homogeneous");
    }
    if (!message_matrix.HasField() || message_matrix.Field() != field_) {
      throw std::invalid_argument(
          "InterleavedGRSCode::EncodeHomogeneous: message matrix field "
          "mismatch");
    }
    const std::size_t k = dimensions_[0];
    if (message_matrix.Dimension() != 2 || message_matrix.Shape(0) != k ||
        message_matrix.Shape(1) != interleaving_order_) {
      throw std::invalid_argument(
          "InterleavedGRSCode::EncodeHomogeneous: message matrix must have "
          "shape {k, ell}");
    }
    return linalg::Dot(linalg::Transpose(component_codes_[0].GeneratorMatrix()),
                       message_matrix);
  }

  // Returns true iff every column m of word_matrix (shape {n, ell}) belongs to
  // ComponentCode(m).
  bool Contains(const matrix_type &word_matrix) const {
    if (!word_matrix.HasField() || word_matrix.Field() != field_ ||
        word_matrix.Dimension() != 2 || word_matrix.Shape(0) != length_ ||
        word_matrix.Shape(1) != interleaving_order_) {
      return false;
    }
    if (max_redundancy_ == 0) {
      return true;
    }
    const auto syn_mat = ComputeSyndromeMatrix(word_matrix);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const std::size_t r_m = redundancies_[m];
      for (std::size_t i = 0; i < r_m; ++i) {
        if (!field_->IsZero(syn_mat.Raw(i, m))) {
          return false;
        }
      }
    }
    return true;
  }

  // Computes syndrome vectors for all ell streams (entry m has length r_m).
  std::vector<vector_type> Syndromes(const matrix_type &word_matrix) const {
    ValidateWordMatrix(word_matrix, "Syndromes");
    std::vector<vector_type> result;
    result.reserve(interleaving_order_);
    if (max_redundancy_ == 0) {
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        result.push_back(vector_type::Zeros(field_, {0}));
      }
      return result;
    }

    const auto syn_mat = ComputeSyndromeMatrix(word_matrix);
    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const std::size_t r_m = redundancies_[m];
      if (r_m == 0) {
        result.push_back(vector_type::Zeros(field_, {0}));
        continue;
      }
      RawArray s_raw = RawArray::from_shape(ShapeType{r_m});
      for (std::size_t i = 0; i < r_m; ++i) {
        s_raw(i) = syn_mat.Raw(i, m);
      }
      result.push_back(vector_type::FromRaw(field_, std::move(s_raw)));
    }
    return result;
  }

  // Flattens the interleaved code into an [n * ell, sum k_m] LinearCode over
  // F_q using row-major symbol-interleaved layout (row j, stream m -> j*ell+m).
  LinearCode<GaloisField> ToLinearCode() const {
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    const std::size_t total_k =
        std::accumulate(dimensions_.begin(), dimensions_.end(), std::size_t{0});
    const std::size_t total_n = length_ * interleaving_order_;

    RawArray g_raw = RawArray::from_shape(ShapeType{total_k, total_n});
    g_raw.fill(field_->AdditiveIdentity());
    storage_type *g_data = g_raw.data();

    std::size_t row_offset = 0;
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const auto &G_m = component_codes_[m].GeneratorMatrix();
      const std::size_t k_m = dimensions_[m];
      for (std::size_t i = 0; i < k_m; ++i) {
        storage_type *dst_row = g_data + (row_offset + i) * total_n;
        for (std::size_t j = 0; j < length_; ++j) {
          dst_row[j * interleaving_order_ + m] = G_m.Raw(i, j);
        }
      }
      row_offset += k_m;
    }

    auto G_flat = matrix_type::FromRaw(field_, std::move(g_raw));
    const std::size_t min_dist = min_redundancy_ + 1;
    return LinearCode<GaloisField>::FromGeneratorMatrix(
        field_, G_flat, Metric::HAMMING, min_dist);
  }

  result_type Decode(const matrix_type &received_matrix,
                     IRSDecoderMode mode = IRSDecoderMode::kSsb,
                     bool fallback_to_independent = true) const {
    return DecodeWithErasures(received_matrix, {}, mode,
                              fallback_to_independent);
  }

  result_type DecodeWithErasures(
      const matrix_type &received_matrix,
      const std::vector<std::size_t> &erasure_positions,
      IRSDecoderMode mode = IRSDecoderMode::kSsb,
      bool fallback_to_independent = true) const {
    ValidateWordMatrix(received_matrix, "DecodeWithErasures");

    const std::size_t s = erasure_positions.size();
    std::vector<bool> is_erased(length_, false);
    for (std::size_t pos : erasure_positions) {
      if (pos >= length_) {
        return result_type::Failure("InterleavedGRSCode: erasure position " +
                                    std::to_string(pos) + " out of range [0, " +
                                    std::to_string(length_) + ")");
      }
      if (is_erased[pos]) {
        return result_type::Failure(
            "InterleavedGRSCode: duplicate erasure position " +
            std::to_string(pos));
      }
      is_erased[pos] = true;
    }

    if (s > min_redundancy_) {
      return result_type::Failure("InterleavedGRSCode: number of erasures (" +
                                  std::to_string(s) +
                                  ") exceeds minimum redundancy (" +
                                  std::to_string(min_redundancy_) + ")");
    }

    if (mode == IRSDecoderMode::kIndependent) {
      return DecodeIndependent(received_matrix, erasure_positions, is_erased,
                               /*used_fallback=*/false);
    }

    // Zero out erased rows in the working matrix.
    matrix_type working = received_matrix;
    const storage_type zero = field_->AdditiveIdentity();
    for (std::size_t pos : erasure_positions) {
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        working.Raw(pos, m) = zero;
      }
    }

    auto raw_syndromes = ExtractRawSyndromes(working);
    if (s == 0 && detail::AreAllSequencesZero(*field_, raw_syndromes)) {
      return BuildSuccessResult(received_matrix, received_matrix, is_erased, 0,
                                /*used_fallback=*/false);
    }

    // Compute reduced redundancies after erasures and collaborative radius.
    std::vector<std::size_t> eff_redundancies(interleaving_order_, 0);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      eff_redundancies[m] = redundancies_[m] - s;
    }
    const std::size_t max_burst_errors =
        detail::ComputeCollaborativeRadius(eff_redundancies);

    const auto zero_idx = ZeroEvaluationIndex();

    // Case A: alpha_{j_0} = 0 is itself erased -> shift all syndromes by 1 so
    // remaining erased and erroneous positions have non-zero alpha_j.
    if (zero_idx.has_value() && is_erased[*zero_idx]) {
      std::vector<std::size_t> nonzero_erasures;
      nonzero_erasures.reserve(s - 1);
      for (std::size_t pos : erasure_positions) {
        if (pos != *zero_idx) {
          nonzero_erasures.push_back(pos);
        }
      }
      auto shifted_syn = ShiftSyndromesByOne(raw_syndromes);
      auto raw_s0 = ExtractS0(raw_syndromes);

      if (auto res = TryCollaborativePass(
              working, received_matrix, shifted_syn, nonzero_erasures,
              is_erased, s, max_burst_errors, /*shifted=*/true, raw_s0, mode)) {
        return std::move(*res);
      }

      if (fallback_to_independent) {
        auto fb = DecodeIndependent(received_matrix, erasure_positions,
                                    is_erased, /*used_fallback=*/true);
        if (fb.Ok()) {
          return fb;
        }
      }
      return result_type::Failure(
          "InterleavedGRSCode: collaborative decoding failed");
    }

    // Case B (Pass 1): assume E_{j_0, :} = 0 and solve unshifted key equation.
    const auto raw_s0 = ExtractS0(raw_syndromes);
    if (auto res = TryCollaborativePass(
            working, received_matrix, raw_syndromes, erasure_positions,
            is_erased, s, max_burst_errors, /*shifted=*/false, raw_s0, mode)) {
      return std::move(*res);
    }

    // Case C (Pass 2): if an evaluation point is alpha_{j_0} = 0 (not erased)
    // and we have burst-error budget >= 1, shift syndromes by 1 assuming a
    // burst error sits at j_0 plus up to t'_max non-zero positions.
    if (zero_idx.has_value() && max_burst_errors >= 1 &&
        min_redundancy_ >= s + 1) {
      auto shifted_syn = ShiftSyndromesByOne(raw_syndromes);
      std::vector<std::size_t> shifted_eff(interleaving_order_, 0);
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        shifted_eff[m] =
            (shifted_syn[m].size() > s) ? (shifted_syn[m].size() - s) : 0;
      }
      const std::size_t shifted_max_errors =
          std::min(detail::ComputeCollaborativeRadius(shifted_eff),
                   max_burst_errors - 1);

      if (auto res = TryCollaborativePass(working, received_matrix, shifted_syn,
                                          erasure_positions, is_erased, s,
                                          shifted_max_errors, /*shifted=*/true,
                                          raw_s0, mode)) {
        return std::move(*res);
      }
    }

    if (fallback_to_independent) {
      auto fb = DecodeIndependent(received_matrix, erasure_positions, is_erased,
                                  /*used_fallback=*/true);
      if (fb.Ok()) {
        return fb;
      }
    }

    return result_type::Failure(
        "InterleavedGRSCode: collaborative decoding failed");
  }

 private:
  friend class PowerDecodingDecoder<GaloisField>;

  void ValidateWordMatrix(const matrix_type &word_matrix,
                          const char *fn_name) const {
    if (!word_matrix.HasField() || word_matrix.Field() != field_) {
      throw std::invalid_argument(std::string("InterleavedGRSCode::") +
                                  fn_name + ": word matrix field mismatch");
    }
    if (word_matrix.Dimension() != 2 || word_matrix.Shape(0) != length_ ||
        word_matrix.Shape(1) != interleaving_order_) {
      throw std::invalid_argument(std::string("InterleavedGRSCode::") +
                                  fn_name +
                                  ": word matrix must have shape {n, ell}");
    }
  }

  // Computes H_max * word_matrix (shape {max_redundancy_, ell}) in one bulk
  // matrix multiply; stream m uses rows 0 .. r_m - 1 of column m.
  matrix_type ComputeSyndromeMatrix(const matrix_type &word_matrix) const {
    if (max_redundancy_ == 0) {
      return matrix_type::Zeros(field_, {0, interleaving_order_});
    }
    return linalg::Dot(
        component_codes_[max_redundancy_idx_].ParityCheckMatrix(), word_matrix);
  }

  std::vector<std::vector<storage_type>> ExtractRawSyndromes(
      const matrix_type &word_matrix) const {
    std::vector<std::vector<storage_type>> out(interleaving_order_);
    if (max_redundancy_ == 0) {
      return out;
    }
    const auto syn_mat = ComputeSyndromeMatrix(word_matrix);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const std::size_t r_m = redundancies_[m];
      out[m].resize(r_m);
      for (std::size_t i = 0; i < r_m; ++i) {
        out[m][i] = syn_mat.Raw(i, m);
      }
    }
    return out;
  }

  std::vector<storage_type> ExtractS0(
      const std::vector<std::vector<storage_type>> &raw_syndromes) const {
    const storage_type zero = field_->AdditiveIdentity();
    std::vector<storage_type> s0(interleaving_order_, zero);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      if (!raw_syndromes[m].empty()) {
        s0[m] = raw_syndromes[m][0];
      }
    }
    return s0;
  }

  static std::vector<std::vector<storage_type>> ShiftSyndromesByOne(
      const std::vector<std::vector<storage_type>> &raw_syndromes) {
    std::vector<std::vector<storage_type>> shifted(raw_syndromes.size());
    for (std::size_t m = 0; m < raw_syndromes.size(); ++m) {
      if (raw_syndromes[m].size() > 1) {
        shifted[m].assign(raw_syndromes[m].begin() + 1, raw_syndromes[m].end());
      }
    }
    return shifted;
  }

  vector_type ExtractColumn(const matrix_type &mat, std::size_t col) const {
    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    RawArray raw = RawArray::from_shape(ShapeType{length_});
    for (std::size_t j = 0; j < length_; ++j) {
      raw(j) = mat.Raw(j, col);
    }
    return vector_type::FromRaw(field_, std::move(raw));
  }

  result_type BuildSuccessResult(const matrix_type &corrected,
                                 const matrix_type &original_received,
                                 const std::vector<bool> &is_erased,
                                 std::size_t total_erasure_count,
                                 bool used_fallback) const {
    result_type res;
    res.codeword_matrix = corrected;
    res.messages.reserve(interleaving_order_);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      res.messages.push_back(
          component_codes_[m].Unencode(ExtractColumn(corrected, m)));
    }

    res.burst_error_positions.clear();
    for (std::size_t j = 0; j < length_; ++j) {
      if (is_erased[j]) {
        continue;
      }
      bool row_diff = false;
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        if (corrected.Raw(j, m) != original_received.Raw(j, m)) {
          row_diff = true;
          break;
        }
      }
      if (row_diff) {
        res.burst_error_positions.push_back(j);
      }
    }
    res.num_burst_errors = res.burst_error_positions.size();
    res.num_erasures = total_erasure_count;
    res.failure_probability_bound =
        FailureProbabilityBound(res.num_burst_errors);
    res.used_fallback = used_fallback;
    res.success = true;
    return res;
  }

  result_type DecodeIndependent(
      const matrix_type &received_matrix,
      const std::vector<std::size_t> &erasure_positions,
      const std::vector<bool> &is_erased, bool used_fallback) const {
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    RawArray cw_raw =
        RawArray::from_shape(ShapeType{length_, interleaving_order_});
    storage_type *cw_data = cw_raw.data();
    std::vector<vector_type> messages;
    messages.reserve(interleaving_order_);

    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const auto col = ExtractColumn(received_matrix, m);
      BerlekampMasseyDecoder<GaloisField> bm(component_codes_[m]);
      auto col_res = bm.DecodeWithErasures(col, erasure_positions);
      if (!col_res.Ok()) {
        return result_type::Failure(
            "InterleavedGRSCode: independent Berlekamp-Massey failed on "
            "stream " +
            std::to_string(m) + ": " + col_res.error_message);
      }
      const auto col_span = col_res.codeword.RawSpan();
      for (std::size_t j = 0; j < length_; ++j) {
        cw_data[j * interleaving_order_ + m] = col_span[j];
      }
      messages.push_back(std::move(col_res.message));
    }

    auto corrected = matrix_type::FromRaw(field_, std::move(cw_raw));
    result_type res;
    res.codeword_matrix = std::move(corrected);
    res.messages = std::move(messages);
    for (std::size_t j = 0; j < length_; ++j) {
      if (is_erased[j]) {
        continue;
      }
      bool row_diff = false;
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        if (res.codeword_matrix.Raw(j, m) != received_matrix.Raw(j, m)) {
          row_diff = true;
          break;
        }
      }
      if (row_diff) {
        res.burst_error_positions.push_back(j);
      }
    }
    res.num_burst_errors = res.burst_error_positions.size();
    res.num_erasures = erasure_positions.size();
    res.failure_probability_bound =
        FailureProbabilityBound(res.num_burst_errors);
    res.used_fallback = used_fallback;
    res.success = true;
    return res;
  }

  // Given a candidate error locator Lambda(z), performs Chien search and
  // Forney error evaluation across all ell streams (or only stream 0 when
  // single_stream_only is true, used by PowerDecodingDecoder).
  std::optional<matrix_type> RecoverFromLocator(
      const polynomial_type &raw_lambda, const matrix_type &working,
      const std::vector<polynomial_type> &mod_syns,
      const std::vector<std::size_t> &syn_lengths, const polynomial_type &gamma,
      const std::vector<std::size_t> &nonzero_erasures,
      const std::vector<bool> &is_erased, std::size_t max_errors, bool shifted,
      const std::vector<storage_type> &raw_s0, bool single_stream_only) const {
    if (raw_lambda.IsZero()) {
      return std::nullopt;
    }
    const storage_type l0 = raw_lambda.RawCoeff(0);
    if (field_->IsZero(l0)) {
      return std::nullopt;
    }

    polynomial_type lambda = raw_lambda;
    if (!field_->IsOne(l0)) {
      const auto inv_l0 = element_type::FromInternal(field_->Inv(l0), field_);
      lambda *= inv_l0;
    }

    const int deg_lambda = lambda.Degree();
    if (deg_lambda < 0 || static_cast<std::size_t>(deg_lambda) > max_errors) {
      return std::nullopt;
    }

    const auto &ref_code = component_codes_[0];
    std::vector<std::size_t> error_positions;
    if (deg_lambda > 0) {
      const auto root_sub_indices =
          utils::ChienSearchIndices(lambda, ref_code.NonZeroEvalInverses());
      if (root_sub_indices.size() != static_cast<std::size_t>(deg_lambda)) {
        return std::nullopt;
      }
      error_positions.reserve(root_sub_indices.size());
      const auto &nz_indices = ref_code.NonZeroEvalIndices();
      for (std::size_t sub_idx : root_sub_indices) {
        const std::size_t pos = nz_indices[sub_idx];
        if (is_erased[pos]) {
          return std::nullopt;
        }
        error_positions.push_back(pos);
      }
    }

    const std::size_t s_nz = nonzero_erasures.size();
    const polynomial_type psi = (s_nz == 0) ? lambda : (lambda * gamma);
    const polynomial_type psi_deriv = psi.Derivative();

    std::vector<std::size_t> all_positions;
    all_positions.reserve(error_positions.size() + s_nz);
    all_positions.insert(all_positions.end(), error_positions.begin(),
                         error_positions.end());
    all_positions.insert(all_positions.end(), nonzero_erasures.begin(),
                         nonzero_erasures.end());

    // Precompute denominator inverses w_p * Psi'(alpha_p^-1) once for all
    // streams (shared across all ell columns!).
    const auto eval_span = ref_code.EvaluationPoints().RawSpan();
    const auto dual_span = ref_code.DualMultipliers().RawSpan();
    const std::size_t num_pos = all_positions.size();
    std::vector<storage_type> inv_alphas(num_pos);
    std::vector<storage_type> denoms(num_pos);
    for (std::size_t idx = 0; idx < num_pos; ++idx) {
      const std::size_t pos = all_positions[idx];
      inv_alphas[idx] = ref_code.EvalInverseAt(pos);
      const storage_type den = psi_deriv.EvalRaw(inv_alphas[idx]);
      if (field_->IsZero(den)) {
        return std::nullopt;
      }
      denoms[idx] = field_->Mul(dual_span[pos], den);
    }
    std::vector<storage_type> inv_denoms(num_pos);
    if (num_pos > 0) {
      field_->BatchInv(std::span<storage_type>(inv_denoms),
                       std::span<const storage_type>(denoms));
    }

    const std::size_t streams_to_recover =
        single_stream_only ? 1 : interleaving_order_;
    matrix_type corrected = working;
    const storage_type zero = field_->AdditiveIdentity();

    for (std::size_t m = 0; m < streams_to_recover; ++m) {
      // Recovering deg_lambda + s_nz unknowns on stream m requires
      // syn_lengths[m] >= deg_lambda + s_nz.
      if (syn_lengths[m] < static_cast<std::size_t>(deg_lambda) + s_nz) {
        return std::nullopt;
      }
      const polynomial_type omega =
          detail::TruncatePolynomial(mod_syns[m] * lambda, syn_lengths[m]);

      storage_type zero_point_sum = zero;
      for (std::size_t idx = 0; idx < num_pos; ++idx) {
        const std::size_t pos = all_positions[idx];
        const storage_type num = omega.EvalRaw(inv_alphas[idx]);
        storage_type numer = field_->Neg(num);
        if (!shifted) {
          numer = field_->Mul(numer, eval_span[pos]);
        }
        const storage_type e_val = field_->Mul(numer, inv_denoms[idx]);
        corrected.Raw(pos, m) = field_->Sub(corrected.Raw(pos, m), e_val);
        if (shifted) {
          zero_point_sum =
              field_->Add(zero_point_sum, field_->Mul(e_val, dual_span[pos]));
        }
      }

      if (shifted) {
        const std::size_t j0 = *ZeroEvaluationIndex();
        const storage_type rem_s0 = field_->Sub(raw_s0[m], zero_point_sum);
        const storage_type e_j0 = field_->Div(rem_s0, dual_span[j0]);
        corrected.Raw(j0, m) = field_->Sub(corrected.Raw(j0, m), e_j0);
      }
    }

    if (single_stream_only) {
      if (!component_codes_[0].Contains(ExtractColumn(corrected, 0))) {
        return std::nullopt;
      }
    } else {
      if (!Contains(corrected)) {
        return std::nullopt;
      }
    }
    return corrected;
  }

  // Runs a single collaborative decoding pass (unshifted or shifted) using
  // the selected mode. In kSsb mode, automatically escalates to kSsbExact if
  // fast SSB exceeds max_errors, lacks active stream uniqueness, or fails
  // Chien/Forney/Contains verification.
  std::optional<result_type> TryCollaborativePass(
      const matrix_type &working, const matrix_type &original_received,
      const std::vector<std::vector<storage_type>> &syndromes,
      const std::vector<std::size_t> &nonzero_erasures,
      const std::vector<bool> &is_erased, std::size_t total_erasure_count,
      std::size_t max_errors, bool shifted,
      const std::vector<storage_type> &raw_s0, IRSDecoderMode mode) const {
    auto corrected_opt = RunCollaborativePassCore(
        working, syndromes, nonzero_erasures, is_erased, max_errors, shifted,
        raw_s0, mode, /*single_stream_only=*/false);
    if (!corrected_opt.has_value()) {
      return std::nullopt;
    }
    return BuildSuccessResult(*corrected_opt, original_received, is_erased,
                              total_erasure_count, /*used_fallback=*/false);
  }

  // Core collaborative pass shared by InterleavedGRSCode and
  // PowerDecodingDecoder (with single_stream_only == true).
  std::optional<matrix_type> RunCollaborativePassCore(
      const matrix_type &working,
      const std::vector<std::vector<storage_type>> &syndromes,
      const std::vector<std::size_t> &nonzero_erasures,
      const std::vector<bool> &is_erased, std::size_t max_errors, bool shifted,
      const std::vector<storage_type> &raw_s0, IRSDecoderMode mode,
      bool single_stream_only) const {
    const std::size_t s_nz = nonzero_erasures.size();
    const polynomial_type gamma =
        detail::BuildErasureLocator(component_codes_[0], nonzero_erasures);

    std::vector<polynomial_type> mod_syns;
    mod_syns.reserve(interleaving_order_);
    std::vector<std::size_t> syn_lengths(interleaving_order_, 0);
    std::vector<std::vector<storage_type>> seqs(interleaving_order_);

    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      const std::size_t L_m = syndromes[m].size();
      syn_lengths[m] = L_m;
      std::vector<storage_type> syn_copy = syndromes[m];
      const polynomial_type syn_poly =
          polynomial_type::FromRaw(std::move(syn_copy), field_);
      const polynomial_type mod_syn =
          (s_nz == 0) ? syn_poly
                      : detail::TruncatePolynomial(syn_poly * gamma, L_m);
      if (L_m > s_nz) {
        const std::size_t M_m = L_m - s_nz;
        seqs[m].resize(M_m);
        for (std::size_t i = 0; i < M_m; ++i) {
          seqs[m][i] = mod_syn.RawCoeff(s_nz + i);
        }
      }
      mod_syns.push_back(mod_syn);
    }

    auto try_recover = [&](const polynomial_type &cand) {
      return RecoverFromLocator(cand, working, mod_syns, syn_lengths, gamma,
                                nonzero_erasures, is_erased, max_errors,
                                shifted, raw_s0, single_stream_only);
    };

    if (mode == IRSDecoderMode::kSsb) {
      auto ssb_out = detail::SolveMultiSequenceSsb(field_, seqs, max_errors);
      const bool ssb_degree_ok =
          (ssb_out.length <= max_errors) &&
          (ssb_out.locator.Degree() == static_cast<int>(ssb_out.length));
      // If fast SSB satisfies degree bounds and active streams provide enough
      // equations (or L <= unique_radius_ where stream 0 alone suffices), try
      // direct recovery; otherwise escalate to exact Koetter-Nielsen Popov.
      if (ssb_degree_ok && (ssb_out.active_streams_sufficient ||
                            ssb_out.length <= unique_radius_)) {
        if (auto cw = try_recover(ssb_out.locator)) {
          return cw;
        }
      }
      // Automatic escalation to exact Koetter-Nielsen Popov minimization.
      if (auto exact_lam =
              detail::SolveMultiSequencePopov(field_, seqs, max_errors)) {
        if (auto cw = try_recover(*exact_lam)) {
          return cw;
        }
      }
      return std::nullopt;
    }

    if (mode == IRSDecoderMode::kSsbExact) {
      if (auto exact_lam =
              detail::SolveMultiSequencePopov(field_, seqs, max_errors)) {
        return try_recover(*exact_lam);
      }
      return std::nullopt;
    }

    if (mode == IRSDecoderMode::kBkyHankel) {
      if (auto bky_lam =
              detail::SolveMultiSequenceBkyHankel(field_, seqs, max_errors)) {
        return try_recover(*bky_lam);
      }
      return std::nullopt;
    }

    return std::nullopt;
  }

  // Solves the virtual IRS for PowerDecodingDecoder, recovering only stream 0
  // (s = 1) so higher-power streams with smaller redundancy r_m < t + s never
  // block recovery.
  std::optional<vector_type> DecodeStreamZeroForPowerDecoding(
      const matrix_type &virtual_matrix,
      const std::vector<std::size_t> &erasure_positions,
      const std::vector<bool> &is_erased, IRSDecoderMode mode) const {
    const std::size_t s = erasure_positions.size();
    const std::size_t r0 = redundancies_[0];
    if (s > r0) {
      return std::nullopt;
    }

    auto raw_syndromes = ExtractRawSyndromes(virtual_matrix);
    if (s == 0 && detail::AreAllSequencesZero(*field_, raw_syndromes)) {
      return ExtractColumn(virtual_matrix, 0);
    }

    std::vector<std::size_t> eff_redundancies(interleaving_order_, 0);
    for (std::size_t m = 0; m < interleaving_order_; ++m) {
      eff_redundancies[m] = (redundancies_[m] > s) ? (redundancies_[m] - s) : 0;
    }
    const std::size_t max_errors =
        std::min(detail::ComputeCollaborativeRadius(eff_redundancies), r0 - s);

    const auto zero_idx = ZeroEvaluationIndex();
    const auto raw_s0 = ExtractS0(raw_syndromes);

    // Case A: alpha_{j_0} = 0 is erased.
    if (zero_idx.has_value() && is_erased[*zero_idx]) {
      std::vector<std::size_t> nonzero_erasures;
      nonzero_erasures.reserve(s - 1);
      for (std::size_t pos : erasure_positions) {
        if (pos != *zero_idx) {
          nonzero_erasures.push_back(pos);
        }
      }
      auto shifted_syn = ShiftSyndromesByOne(raw_syndromes);
      if (auto mat = RunCollaborativePassCore(
              virtual_matrix, shifted_syn, nonzero_erasures, is_erased,
              max_errors, /*shifted=*/true, raw_s0, mode,
              /*single_stream_only=*/true)) {
        return ExtractColumn(*mat, 0);
      }
      return std::nullopt;
    }

    // Case B (Pass 1): unshifted.
    if (auto mat = RunCollaborativePassCore(
            virtual_matrix, raw_syndromes, erasure_positions, is_erased,
            max_errors, /*shifted=*/false, raw_s0, mode,
            /*single_stream_only=*/true)) {
      return ExtractColumn(*mat, 0);
    }

    // Case C (Pass 2): shifted assuming error at alpha_{j_0} = 0.
    if (zero_idx.has_value() && max_errors >= 1 && r0 >= s + 1) {
      auto shifted_syn = ShiftSyndromesByOne(raw_syndromes);
      std::vector<std::size_t> shifted_eff(interleaving_order_, 0);
      for (std::size_t m = 0; m < interleaving_order_; ++m) {
        shifted_eff[m] =
            (shifted_syn[m].size() > s) ? (shifted_syn[m].size() - s) : 0;
      }
      const std::size_t shifted_max = std::min(
          detail::ComputeCollaborativeRadius(shifted_eff), max_errors - 1);
      if (auto mat = RunCollaborativePassCore(
              virtual_matrix, shifted_syn, erasure_positions, is_erased,
              shifted_max, /*shifted=*/true, raw_s0, mode,
              /*single_stream_only=*/true)) {
        return ExtractColumn(*mat, 0);
      }
    }

    return std::nullopt;
  }

  std::shared_ptr<GaloisField> field_;
  std::size_t length_ = 0;
  std::size_t interleaving_order_ = 0;
  std::vector<std::size_t> dimensions_;
  std::vector<std::size_t> redundancies_;
  bool is_homogeneous_ = true;
  std::size_t min_redundancy_ = 0;
  std::size_t max_redundancy_ = 0;
  std::size_t max_redundancy_idx_ = 0;
  std::size_t unique_radius_ = 0;
  std::size_t collaborative_radius_ = 0;
  std::vector<component_code_type> component_codes_;
};

// Power Decoding (Schmidt-Sidorenko-Bossert 2006/2010) for a single
// GeneralizedReedSolomonCode [n, k], decoding beyond half the minimum distance
// up to t_power(ell_p) = floor(ell_p/(ell_p+1)*n - ell_p/2*(k-1) -
// ell_p/(ell_p+1)) by mapping y_j = r_j * v_j^-1 to a virtual heterogeneous
// IRS(n; k_1, ..., k_{ell_p}) with k_s = s*(k - 1) + 1.
template <typename GaloisField>
class PowerDecodingDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;

  // Exact integer closed-form power-decoding radius:
  //   floor(ell_p/(ell_p + 1) * (n - 1) - ell_p/2 * (k - 1))
  //   = floor(ell_p * (2*(n - 1) - (ell_p + 1)*(k - 1)) / (2*(ell_p + 1)))
  static std::size_t PowerDecodingRadius(std::size_t length,
                                         std::size_t dimension,
                                         std::size_t powering_degree) noexcept {
    if (powering_degree == 0 || dimension == 0 || dimension >= length) {
      return 0;
    }
    const std::size_t lhs = 2 * (length - 1);
    const std::size_t rhs = (powering_degree + 1) * (dimension - 1);
    if (lhs <= rhs) {
      return 0;
    }
    return (powering_degree * (lhs - rhs)) / (2 * (powering_degree + 1));
  }

  // Selects the smallest powering degree ell_p >= 1 maximizing
  // PowerDecodingRadius(n, k, ell_p) subject to ell_p * (k - 1) + 1 < n.
  static std::size_t OptimalPoweringDegree(std::size_t length,
                                           std::size_t dimension) noexcept {
    if (dimension == 0 || dimension >= length) {
      return 1;
    }
    std::size_t best_ell = 1;
    std::size_t best_t = PowerDecodingRadius(length, dimension, 1);

    if (dimension == 1) {
      // Repetition code: k_s = 1 for all s; t_power plateaus at n - 2.
      const std::size_t limit = (length > 2) ? (length - 2) : 1;
      for (std::size_t ell = 2; ell <= limit; ++ell) {
        const std::size_t t = PowerDecodingRadius(length, dimension, ell);
        if (t > best_t) {
          best_t = t;
          best_ell = ell;
        }
      }
      return best_ell;
    }

    const std::size_t deg_step = dimension - 1;
    const std::size_t max_ell = (length - 2) / deg_step;
    for (std::size_t ell = 2; ell <= max_ell; ++ell) {
      const std::size_t t = PowerDecodingRadius(length, dimension, ell);
      if (t > best_t) {
        best_t = t;
        best_ell = ell;
      }
    }
    return best_ell;
  }

  explicit PowerDecodingDecoder(
      const AbstractCode<GaloisField> &code,
      std::optional<std::size_t> powering_degree = std::nullopt,
      IRSDecoderMode mode = IRSDecoderMode::kSsb, bool fallback_to_bm = true)
      : Base(code), mode_(mode), fallback_to_bm_(fallback_to_bm) {
    grs_code_ =
        dynamic_cast<const GeneralizedReedSolomonCode<GaloisField> *>(&code);
    if (!grs_code_) {
      throw std::invalid_argument(
          "PowerDecodingDecoder requires a GeneralizedReedSolomonCode");
    }
    if (mode_ == IRSDecoderMode::kIndependent) {
      throw std::invalid_argument(
          "PowerDecodingDecoder requires a collaborative IRSDecoderMode "
          "(kSsb, kSsbExact, or kBkyHankel)");
    }

    const std::size_t n = grs_code_->Length();
    const std::size_t k = grs_code_->Dimension();
    if (powering_degree.has_value()) {
      powering_degree_ = *powering_degree;
      if (powering_degree_ == 0) {
        throw std::invalid_argument(
            "PowerDecodingDecoder: powering degree must be at least 1");
      }
      if (powering_degree_ > 1 && powering_degree_ * (k - 1) + 1 >= n) {
        throw std::invalid_argument(
            "PowerDecodingDecoder: powering degree violates ell_p*(k - 1) + 1 "
            "< n");
      }
    } else {
      powering_degree_ = OptimalPoweringDegree(n, k);
    }

    std::vector<std::size_t> virtual_dims;
    virtual_dims.reserve(powering_degree_);
    for (std::size_t s = 1; s <= powering_degree_; ++s) {
      virtual_dims.push_back(s * (k - 1) + 1);
    }

    auto field = grs_code_->Field();
    auto unit_mult = codeword_type::Ones(field, {n});
    virtual_irs_ = std::make_unique<InterleavedGRSCode<GaloisField>>(
        field, grs_code_->EvaluationPoints(), unit_mult,
        std::move(virtual_dims));
    max_errors_ = virtual_irs_->CollaborativeDecodingRadius();
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return DecodeWithErasures(received_word, {});
  }

  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const override {
    const auto field = grs_code_->Field();
    const std::size_t n = grs_code_->Length();
    const std::size_t r = grs_code_->Redundancy();

    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "PowerDecoding: received word field mismatch");
    }
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "PowerDecoding: received word must be 1D with size equal to code "
          "length");
    }

    const std::size_t s = erasure_positions.size();
    std::vector<bool> is_erased(n, false);
    for (std::size_t pos : erasure_positions) {
      if (pos >= n) {
        return result_type::Failure("PowerDecoding: erasure position " +
                                    std::to_string(pos) + " out of range [0, " +
                                    std::to_string(n) + ")");
      }
      if (is_erased[pos]) {
        return result_type::Failure(
            "PowerDecoding: duplicate erasure position " + std::to_string(pos));
      }
      is_erased[pos] = true;
    }

    if (s > r) {
      return result_type::Failure("PowerDecoding: number of erasures (" +
                                  std::to_string(s) + ") exceeds redundancy (" +
                                  std::to_string(r) + ")");
    }

    if (s == 0 && grs_code_->Contains(received_word)) {
      return result_type::Success(received_word,
                                  grs_code_->Unencode(received_word), 0, 0);
    }

    // Build virtual {n, ell_p} matrix Y_{j, m} = (r_j * v_j^-1)^(m + 1)
    // with erased positions zeroed out.
    using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
    using ShapeType = typename GFMatrix<GaloisField>::shape_type;
    RawArray y_raw = RawArray::from_shape(ShapeType{n, powering_degree_});
    storage_type *y_data = y_raw.data();
    const storage_type zero = field->AdditiveIdentity();
    const auto rx_span = received_word.RawSpan();
    const auto &v_inv = grs_code_->ColumnMultiplierInverses();

    for (std::size_t j = 0; j < n; ++j) {
      storage_type *row = y_data + j * powering_degree_;
      if (is_erased[j]) {
        for (std::size_t m = 0; m < powering_degree_; ++m) {
          row[m] = zero;
        }
        continue;
      }
      const storage_type y_j = field->Mul(rx_span[j], v_inv[j]);
      storage_type cur = y_j;
      row[0] = cur;
      for (std::size_t m = 1; m < powering_degree_; ++m) {
        cur = field->Mul(cur, y_j);
        row[m] = cur;
      }
    }

    auto virtual_matrix =
        GFMatrix<GaloisField>::FromRaw(field, std::move(y_raw));
    if (auto f_hat_opt = virtual_irs_->DecodeStreamZeroForPowerDecoding(
            virtual_matrix, erasure_positions, is_erased, mode_)) {
      // Rescale c_j = f_hat_j * v_j and verify membership in the original GRS.
      using VecRaw = typename codeword_type::raw_array_type;
      using VecShape = typename codeword_type::shape_type;
      VecRaw cw_raw = VecRaw::from_shape(VecShape{n});
      field->MulArray(std::span<storage_type>(cw_raw.data(), n),
                      f_hat_opt->RawSpan(),
                      grs_code_->ColumnMultipliers().RawSpan());
      auto corrected = codeword_type::FromRaw(field, std::move(cw_raw));

      if (grs_code_->Contains(corrected)) {
        std::size_t actual_errors = 0;
        for (std::size_t j = 0; j < n; ++j) {
          if (!is_erased[j] && corrected.Raw(j) != received_word.Raw(j)) {
            ++actual_errors;
          }
        }
        auto message = grs_code_->Unencode(corrected);
        return result_type::Success(std::move(corrected), std::move(message),
                                    actual_errors, s);
      }
    }

    if (fallback_to_bm_) {
      BerlekampMasseyDecoder<GaloisField> bm(*grs_code_);
      auto bm_res = bm.DecodeWithErasures(received_word, erasure_positions);
      if (bm_res.Ok()) {
        return bm_res;
      }
    }

    return result_type::Failure(
        "PowerDecoding: uncorrectable error/erasure pattern");
  }

  std::string Name() const override { return "PowerDecoding"; }

  std::size_t DecodingRadius() const override { return max_errors_; }

  std::size_t PoweringDegree() const noexcept { return powering_degree_; }

  IRSDecoderMode Mode() const noexcept { return mode_; }

  std::size_t MaxErrors() const noexcept { return max_errors_; }

 private:
  const GeneralizedReedSolomonCode<GaloisField> *grs_code_ = nullptr;
  std::size_t powering_degree_ = 1;
  IRSDecoderMode mode_ = IRSDecoderMode::kSsb;
  bool fallback_to_bm_ = true;
  std::size_t max_errors_ = 0;
  std::unique_ptr<InterleavedGRSCode<GaloisField>> virtual_irs_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_INTERLEAVED_RS_HPP_
