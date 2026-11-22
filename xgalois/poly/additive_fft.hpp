#ifndef XGALOIS_POLY_ADDITIVE_FFT_HPP_
#define XGALOIS_POLY_ADDITIVE_FFT_HPP_

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace poly {

// Lin-Chung-Han (2014) Novel Polynomial Basis (NPB) Additive FFT for any
// binary extension field GF(2^m).
//
// Given an F_2-basis (beta_0, ..., beta_{m-1}) of GF(2^m), defines the nested
// subspace tower W_0 = {0} subset W_1 subset ... subset W_m = GF(2^m) with
//   W_i = span_{F_2}(beta_0, ..., beta_{i-1}),
// normalized subspace vanishing polynomials
//   s_bar_i(x) = s_i(x) / s_i(beta_i),  where s_i(x) = prod_{w in W_i}(x + w),
// and the degree-graded Novel Polynomial Basis
//   X_j(x) = prod_{i=0}^{k-1} (s_bar_i(x))^{j_i}  for j = sum j_i 2^i.
template <FiniteField BinaryField>
class AdditiveFftPlan {
 public:
  using field_type = BinaryField;
  using storage_type = typename BinaryField::element_type;

  explicit AdditiveFftPlan(std::shared_ptr<BinaryField> field)
      : field_(std::move(field)) {
    ValidateField();
    m_ = static_cast<std::size_t>(field_->Degree());
    basis_.resize(m_);
    if constexpr (std::is_integral_v<storage_type>) {
      for (std::size_t i = 0; i < m_; ++i) {
        basis_[i] = field_->SetElementValue(
            static_cast<storage_type>(uint64_t{1} << i));
      }
    } else {
      throw std::invalid_argument(
          "AdditiveFftPlan: non-integral binary field requires an explicit "
          "F_2-basis");
    }
    BuildSubspaceTower();
  }

  AdditiveFftPlan(std::shared_ptr<BinaryField> field,
                  std::vector<storage_type> basis)
      : field_(std::move(field)), basis_(std::move(basis)) {
    ValidateField();
    m_ = static_cast<std::size_t>(field_->Degree());
    if (basis_.size() != m_) {
      throw std::invalid_argument(
          "AdditiveFftPlan: basis size must equal field extension degree m");
    }
    BuildSubspaceTower();
  }

  const std::shared_ptr<BinaryField>& Field() const noexcept { return field_; }

  std::size_t ExtensionDegree() const noexcept { return m_; }

  const std::vector<storage_type>& Basis() const noexcept { return basis_; }

  // Returns omega_idx = sum_{r} ((idx >> r) & 1) * beta_r in W_m.
  storage_type EvaluationPoint(std::size_t idx) const {
    if (m_ < sizeof(std::size_t) * 8 && idx >= (std::size_t{1} << m_)) {
      throw std::out_of_range(
          "AdditiveFftPlan::EvaluationPoint: index exceeds 2^m");
    }
    storage_type acc = field_->AdditiveIdentity();
    for (std::size_t r = 0; r < m_ && idx != 0; ++r, idx >>= 1) {
      if ((idx & 1U) != 0U) {
        acc = field_->Add(acc, basis_[r]);
      }
    }
    return acc;
  }

  // Returns all 2^k points of the affine subspace (shift + W_k) in natural
  // binary-indexed order.
  std::vector<storage_type> EvaluationPoints(std::size_t k,
                                             storage_type shift) const {
    ValidateSubspaceOrder(k);
    const std::size_t n = std::size_t{1} << k;
    std::vector<storage_type> pts(n);
    pts[0] = shift;
    for (std::size_t bit = 0; bit < k; ++bit) {
      const std::size_t half = std::size_t{1} << bit;
      const storage_type b = basis_[bit];
      for (std::size_t i = 0; i < half; ++i) {
        pts[half + i] = field_->Add(pts[i], b);
      }
    }
    return pts;
  }

  std::vector<storage_type> EvaluationPoints(std::size_t k) const {
    return EvaluationPoints(k, field_->AdditiveIdentity());
  }

  // In-place O(N log N) forward Lin-Chung-Han butterfly:
  // evaluates an NPB polynomial of length N = 2^k on the affine subspace
  // (shift + W_k) in natural order.
  void ForwardNpb(std::span<storage_type> coeffs, std::size_t k,
                  storage_type shift) const {
    ValidateSubspaceOrder(k);
    const std::size_t n = std::size_t{1} << k;
    if (coeffs.size() != n) {
      throw std::invalid_argument(
          "AdditiveFftPlan::ForwardNpb: span length must equal 2^k");
    }
    if (k == 0) {
      return;
    }

    const std::vector<storage_type> shift_s = EvalNormalizedSubspace(k, shift);
    std::vector<storage_type> deltas;

    for (std::size_t s = k; s > 0; --s) {
      const std::size_t stage = s - 1;
      const std::size_t half = std::size_t{1} << stage;
      const std::size_t step = half << 1;
      const std::size_t num_blocks = n >> s;

      BuildStageDeltas(stage, k, shift_s[stage], num_blocks, deltas);
      for (std::size_t block = 0; block < num_blocks; ++block) {
        const std::size_t base = block * step;
        auto lo = coeffs.subspan(base, half);
        auto hi = coeffs.subspan(base + half, half);
        field_->FmaArray(lo, hi, deltas[block]);
        field_->AddArray(hi, hi, lo);
      }
    }
  }

  void ForwardNpb(std::span<storage_type> coeffs, std::size_t k) const {
    ForwardNpb(coeffs, k, field_->AdditiveIdentity());
  }

  // In-place O(N log N) inverse Lin-Chung-Han butterfly:
  // interpolates NPB coefficients from evaluations on (shift + W_k).
  void InverseNpb(std::span<storage_type> evals, std::size_t k,
                  storage_type shift) const {
    ValidateSubspaceOrder(k);
    const std::size_t n = std::size_t{1} << k;
    if (evals.size() != n) {
      throw std::invalid_argument(
          "AdditiveFftPlan::InverseNpb: span length must equal 2^k");
    }
    if (k == 0) {
      return;
    }

    const std::vector<storage_type> shift_s = EvalNormalizedSubspace(k, shift);
    std::vector<storage_type> deltas;

    for (std::size_t stage = 0; stage < k; ++stage) {
      const std::size_t half = std::size_t{1} << stage;
      const std::size_t step = half << 1;
      const std::size_t num_blocks = n >> (stage + 1);

      BuildStageDeltas(stage, k, shift_s[stage], num_blocks, deltas);
      for (std::size_t block = 0; block < num_blocks; ++block) {
        const std::size_t base = block * step;
        auto lo = evals.subspan(base, half);
        auto hi = evals.subspan(base + half, half);
        field_->AddArray(hi, hi, lo);
        field_->FmaArray(lo, hi, deltas[block]);
      }
    }
  }

  void InverseNpb(std::span<storage_type> evals, std::size_t k) const {
    InverseNpb(evals, k, field_->AdditiveIdentity());
  }

  // Converts a polynomial in standard monomial basis (length <= 2^m) to the
  // Novel Polynomial Basis in O(N log^2 N), padding to the next power of 2.
  std::vector<storage_type> MonomialToNpb(
      std::span<const storage_type> monomial) const {
    if (monomial.empty()) {
      return {};
    }
    const std::size_t n = NextPow2Checked(monomial.size());
    std::vector<storage_type> buf(n, field_->AdditiveIdentity());
    std::copy(monomial.begin(), monomial.end(), buf.begin());

    // Rescale x^j = (beta_0 * s_bar_0(x))^j when beta_0 != 1.
    if (!field_->IsOne(basis_[0])) {
      storage_type b0_pow = field_->MultiplicativeIdentity();
      for (std::size_t j = 0; j < n; ++j) {
        buf[j] = field_->Mul(buf[j], b0_pow);
        b0_pow = field_->Mul(b0_pow, basis_[0]);
      }
    }

    MonomialToNpbRec(std::span<storage_type>(buf), 0);
    return buf;
  }

  // Converts a polynomial from Novel Polynomial Basis back to standard
  // monomial basis in O(N log^2 N).
  std::vector<storage_type> NpbToMonomial(
      std::span<const storage_type> npb) const {
    if (npb.empty()) {
      return {};
    }
    const std::size_t n = NextPow2Checked(npb.size());
    std::vector<storage_type> buf(n, field_->AdditiveIdentity());
    std::copy(npb.begin(), npb.end(), buf.begin());

    NpbToMonomialRec(std::span<storage_type>(buf), 0);

    if (!field_->IsOne(basis_[0])) {
      storage_type inv_b0_pow = field_->MultiplicativeIdentity();
      for (std::size_t j = 0; j < n; ++j) {
        buf[j] = field_->Mul(buf[j], inv_b0_pow);
        inv_b0_pow = field_->Mul(inv_b0_pow, inv_beta0_);
      }
    }
    return buf;
  }

  // Evaluates a monomial-basis polynomial on the subspace W_k (2^k points).
  std::vector<storage_type> Evaluate(const PolynomialDense<BinaryField>& poly,
                                     std::size_t k) const {
    if (!poly.HasField() || poly.Field() != field_) {
      throw std::invalid_argument(
          "AdditiveFftPlan::Evaluate: polynomial field mismatch");
    }
    ValidateSubspaceOrder(k);
    const std::size_t target_n = std::size_t{1} << k;
    if (poly.IsZero()) {
      return std::vector<storage_type>(target_n, field_->AdditiveIdentity());
    }

    // Fold powers x^(2^m) == x over GF(2^m) if deg(poly) >= 2^m.
    std::vector<storage_type> reduced;
    const auto& raw = poly.RawCoeffs();
    const bool has_field_size_limit = (m_ < sizeof(std::size_t) * 8);
    const std::size_t max_n =
        has_field_size_limit ? (std::size_t{1} << m_) : raw.size();
    std::span<const storage_type> coeffs_span(raw);
    if (has_field_size_limit && raw.size() > max_n) {
      const std::size_t mod_exp = max_n - 1;
      reduced.assign(max_n, field_->AdditiveIdentity());
      reduced[0] = raw[0];
      for (std::size_t i = 1; i < raw.size(); ++i) {
        const std::size_t idx = 1 + ((i - 1) % mod_exp);
        reduced[idx] = field_->Add(reduced[idx], raw[i]);
      }
      coeffs_span = std::span<const storage_type>(reduced);
    }

    const std::size_t needed_k =
        (coeffs_span.size() <= 1)
            ? 0
            : static_cast<std::size_t>(std::bit_width(coeffs_span.size() - 1));
    const std::size_t eval_k = std::max(k, needed_k);
    const std::size_t eval_n = std::size_t{1} << eval_k;

    std::vector<storage_type> padded(eval_n, field_->AdditiveIdentity());
    std::copy(coeffs_span.begin(), coeffs_span.end(), padded.begin());

    std::vector<storage_type> npb =
        MonomialToNpb(std::span<const storage_type>(padded));
    ForwardNpb(std::span<storage_type>(npb), eval_k);
    npb.resize(target_n);
    return npb;
  }

  // Interpolates the unique monomial-basis polynomial of degree < 2^k taking
  // the given values on W_k.
  PolynomialDense<BinaryField> Interpolate(std::span<const storage_type> evals,
                                           std::size_t k,
                                           std::string var = "x") const {
    ValidateSubspaceOrder(k);
    const std::size_t n = std::size_t{1} << k;
    if (evals.size() != n) {
      throw std::invalid_argument(
          "AdditiveFftPlan::Interpolate: evaluation count must equal 2^k");
    }
    std::vector<storage_type> buf(evals.begin(), evals.end());
    InverseNpb(std::span<storage_type>(buf), k);
    std::vector<storage_type> mono =
        NpbToMonomial(std::span<const storage_type>(buf));
    return PolynomialDense<BinaryField>::FromRaw(std::move(mono), field_,
                                                 std::move(var));
  }

  // Exact polynomial multiplication in monomial basis via additive FFT
  // whenever deg(a) + deg(b) < 2^m.
  PolynomialDense<BinaryField> Multiply(
      const PolynomialDense<BinaryField>& a,
      const PolynomialDense<BinaryField>& b) const {
    if (!a.HasField() || !b.HasField() || a.Field() != field_ ||
        b.Field() != field_) {
      throw std::invalid_argument(
          "AdditiveFftPlan::Multiply: polynomial field mismatch");
    }
    if (a.IsZero() || b.IsZero()) {
      return PolynomialDense<BinaryField>::Zero(field_, a.GetVariable());
    }

    const std::size_t out_len = a.Size() + b.Size() - 1;
    const std::size_t k =
        (out_len <= 1) ? 0
                       : static_cast<std::size_t>(std::bit_width(out_len - 1));
    if (k > m_) {
      throw std::invalid_argument(
          "AdditiveFftPlan::Multiply: product degree exceeds 2^m - 1");
    }
    const std::size_t n = std::size_t{1} << k;

    std::vector<storage_type> pa(n, field_->AdditiveIdentity());
    std::copy(a.RawCoeffs().begin(), a.RawCoeffs().end(), pa.begin());
    std::vector<storage_type> pb(n, field_->AdditiveIdentity());
    std::copy(b.RawCoeffs().begin(), b.RawCoeffs().end(), pb.begin());

    std::vector<storage_type> fa =
        MonomialToNpb(std::span<const storage_type>(pa));
    std::vector<storage_type> fb =
        MonomialToNpb(std::span<const storage_type>(pb));

    ForwardNpb(std::span<storage_type>(fa), k);
    ForwardNpb(std::span<storage_type>(fb), k);
    field_->MulArray(std::span<storage_type>(fa),
                     std::span<const storage_type>(fa),
                     std::span<const storage_type>(fb));
    InverseNpb(std::span<storage_type>(fa), k);

    std::vector<storage_type> prod =
        NpbToMonomial(std::span<const storage_type>(fa));
    prod.resize(out_len);
    return PolynomialDense<BinaryField>::FromRaw(std::move(prod), field_,
                                                 a.GetVariable());
  }

 private:
  void ValidateField() const {
    if (!field_) {
      throw std::invalid_argument("AdditiveFftPlan: field pointer is null");
    }
    if (field_->Characteristic() != 2ULL) {
      throw std::invalid_argument(
          "AdditiveFftPlan: requires a binary extension field (char == 2)");
    }
    if (field_->Degree() == 0U) {
      throw std::invalid_argument(
          "AdditiveFftPlan: field extension degree must be >= 1");
    }
  }

  void ValidateSubspaceOrder(std::size_t k) const {
    if (k > m_) {
      throw std::invalid_argument(
          "AdditiveFftPlan: subspace dimension k exceeds extension degree m");
    }
    if (k >= sizeof(std::size_t) * 8) {
      throw std::overflow_error(
          "AdditiveFftPlan: subspace size 2^k overflows size_t");
    }
  }

  std::size_t NextPow2Checked(std::size_t len) const {
    const std::size_t k =
        (len <= 1) ? 0 : static_cast<std::size_t>(std::bit_width(len - 1));
    ValidateSubspaceOrder(k);
    return std::size_t{1} << k;
  }

  // Precomputes twist scalars twist_[i][j] = s_bar_i(beta_j) and recurrence
  // constants scale_[i] = s_bar_i(beta_{i+1})^2 + s_bar_i(beta_{i+1}) such
  // that s_bar_i(x)^2 + s_bar_i(x) = scale_[i] * s_bar_{i+1}(x).
  void BuildSubspaceTower() {
    if (field_->IsZero(basis_[0])) {
      throw std::invalid_argument(
          "AdditiveFftPlan: basis vectors must be non-zero");
    }
    inv_beta0_ = field_->Inv(basis_[0]);

    twist_.assign(m_,
                  std::vector<storage_type>(m_, field_->AdditiveIdentity()));
    for (std::size_t j = 0; j < m_; ++j) {
      twist_[0][j] = field_->Mul(basis_[j], inv_beta0_);
    }

    scale_.resize(m_ > 0 ? m_ - 1 : 0);
    inv_scale_.resize(m_ > 0 ? m_ - 1 : 0);

    for (std::size_t i = 0; i + 1 < m_; ++i) {
      const storage_type t_next = twist_[i][i + 1];
      const storage_type c = field_->Add(field_->Mul(t_next, t_next), t_next);
      if (field_->IsZero(c)) {
        throw std::invalid_argument(
            "AdditiveFftPlan: basis vectors are not linearly independent over "
            "F_2");
      }
      scale_[i] = c;
      inv_scale_[i] = field_->Inv(c);

      for (std::size_t j = i + 1; j < m_; ++j) {
        const storage_type t = twist_[i][j];
        const storage_type q = field_->Add(field_->Mul(t, t), t);
        twist_[i + 1][j] = field_->Mul(q, inv_scale_[i]);
      }
    }
  }

  // Computes [s_bar_0(x), s_bar_1(x), ..., s_bar_{levels - 1}(x)].
  std::vector<storage_type> EvalNormalizedSubspace(std::size_t levels,
                                                   storage_type x) const {
    std::vector<storage_type> vals(levels);
    if (levels == 0) {
      return vals;
    }
    vals[0] = field_->Mul(x, inv_beta0_);
    for (std::size_t i = 0; i + 1 < levels; ++i) {
      const storage_type t = vals[i];
      vals[i + 1] =
          field_->Mul(field_->Add(field_->Mul(t, t), t), inv_scale_[i]);
    }
    return vals;
  }

  void BuildStageDeltas(std::size_t stage, std::size_t k,
                        storage_type base_delta, std::size_t num_blocks,
                        std::vector<storage_type>& deltas) const {
    deltas.resize(num_blocks);
    deltas[0] = base_delta;
    const std::size_t upper_bits = k - 1 - stage;
    for (std::size_t bit = 0; bit < upper_bits; ++bit) {
      const std::size_t half = std::size_t{1} << bit;
      const storage_type tw = twist_[stage][stage + 1 + bit];
      for (std::size_t i = 0; i < half; ++i) {
        deltas[half + i] = field_->Add(deltas[i], tw);
      }
    }
  }

  // Recursive O(N log^2 N) monomial-in-s_bar_{level} -> NPB conversion using
  // characteristic-2 Taylor expansion in powers of (y^2 + y) = scale_[level] *
  // s_bar_{level + 1}(x).
  void MonomialToNpbRec(std::span<storage_type> buf, std::size_t level) const {
    const std::size_t n = buf.size();
    if (n <= 2) {
      return;
    }

    // Forward Taylor expansion in powers of (y^2 + y):
    // divides degree-(4h-1) blocks by (y^2 + y)^h = y^(2h) + y^h.
    for (std::size_t size = n; size >= 4; size >>= 1) {
      const std::size_t h = size >> 2;
      for (std::size_t base = 0; base < n; base += size) {
        auto p1 = buf.subspan(base + h, h);
        auto p2 = buf.subspan(base + 2 * h, h);
        auto p3 = buf.subspan(base + 3 * h, h);
        field_->AddArray(p2, p2, p3);
        field_->AddArray(p1, p1, p2);
      }
    }

    const std::size_t half = n >> 1;
    std::vector<storage_type> even(half);
    std::vector<storage_type> odd(half);
    const storage_type c = scale_[level];
    storage_type c_pow = field_->MultiplicativeIdentity();
    for (std::size_t t = 0; t < half; ++t) {
      even[t] = field_->Mul(buf[2 * t], c_pow);
      odd[t] = field_->Mul(buf[2 * t + 1], c_pow);
      c_pow = field_->Mul(c_pow, c);
    }

    MonomialToNpbRec(std::span<storage_type>(even), level + 1);
    MonomialToNpbRec(std::span<storage_type>(odd), level + 1);

    for (std::size_t t = 0; t < half; ++t) {
      buf[2 * t] = even[t];
      buf[2 * t + 1] = odd[t];
    }
  }

  // Exact inverse of MonomialToNpbRec in O(N log^2 N).
  void NpbToMonomialRec(std::span<storage_type> buf, std::size_t level) const {
    const std::size_t n = buf.size();
    if (n <= 2) {
      return;
    }

    const std::size_t half = n >> 1;
    std::vector<storage_type> even(half);
    std::vector<storage_type> odd(half);
    for (std::size_t t = 0; t < half; ++t) {
      even[t] = buf[2 * t];
      odd[t] = buf[2 * t + 1];
    }

    NpbToMonomialRec(std::span<storage_type>(even), level + 1);
    NpbToMonomialRec(std::span<storage_type>(odd), level + 1);

    const storage_type inv_c = inv_scale_[level];
    storage_type inv_c_pow = field_->MultiplicativeIdentity();
    for (std::size_t t = 0; t < half; ++t) {
      buf[2 * t] = field_->Mul(even[t], inv_c_pow);
      buf[2 * t + 1] = field_->Mul(odd[t], inv_c_pow);
      inv_c_pow = field_->Mul(inv_c_pow, inv_c);
    }

    // Inverse Taylor expansion (exact reverse of forward Taylor butterfly).
    for (std::size_t size = 4; size <= n; size <<= 1) {
      const std::size_t h = size >> 2;
      for (std::size_t base = 0; base < n; base += size) {
        auto p1 = buf.subspan(base + h, h);
        auto p2 = buf.subspan(base + 2 * h, h);
        auto p3 = buf.subspan(base + 3 * h, h);
        field_->AddArray(p1, p1, p2);
        field_->AddArray(p2, p2, p3);
      }
    }
  }

  std::shared_ptr<BinaryField> field_;
  std::size_t m_ = 0;
  std::vector<storage_type> basis_;
  storage_type inv_beta0_{};
  std::vector<std::vector<storage_type>> twist_;
  std::vector<storage_type> scale_;
  std::vector<storage_type> inv_scale_;
};

}  // namespace poly

using poly::AdditiveFftPlan;

}  // namespace xg

#endif  // XGALOIS_POLY_ADDITIVE_FFT_HPP_
