#ifndef XGALOIS_CODING_HERMITIAN_HPP_
#define XGALOIS_CODING_HERMITIAN_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace coding {

// Finite rational point P = (x, y) on the Hermitian curve H_q : y^q + y =
// x^(q+1) over F_{q^2}.
template <typename GaloisField>
struct HermitianPoint {
  using storage_type = typename GaloisField::element_type;

  storage_type x{};
  storage_type y{};
};

// Basis monomial x^i y^j in L(m P_inf) with 0 <= j <= q - 1 and pole order
// rho(i, j) = q * i + (q + 1) * j.
struct HermitianMonomial {
  std::size_t i = 0;
  std::size_t j = 0;
  std::size_t pole_order = 0;
};

template <typename GaloisField>
class HermitianCode;

template <typename GaloisField>
class HermitianFengRaoDecoder;

template <typename GaloisField>
class HermitianSudanListDecoder;

namespace detail {

// Enumerate all elements of a finite field in internal storage representation.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type>
EnumerateHermitianFieldElements(const std::shared_ptr<GaloisField> &field) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t order = field->Order();
  if (order == 0) {
    throw std::invalid_argument(
        "HermitianCode: field order too large for curve enumeration");
  }
  std::vector<StorageType> elems;
  elems.reserve(static_cast<std::size_t>(order));
  if constexpr (std::is_integral_v<StorageType>) {
    for (uint64_t i = 0; i < order; ++i) {
      elems.push_back(field->SetElementValue(static_cast<StorageType>(i)));
    }
  } else {
    for (uint64_t i = 0; i < order; ++i) {
      elems.push_back(field->FromInt(i));
    }
  }
  return elems;
}

// Returns true iff v belongs to the Hermitian Weierstrass semigroup
// H(P_inf) = <q, q + 1>.
inline bool IsInHermitianSemigroup(std::size_t v, std::size_t q) noexcept {
  const std::size_t b = v % q;
  return v >= (q + 1) * b;
}

// Decomposes v in H(P_inf) into the unique (a, b) with 0 <= b < q and
// q * a + (q + 1) * b == v.
inline std::pair<std::size_t, std::size_t> DecomposePoleOrder(
    std::size_t v, std::size_t q) noexcept {
  const std::size_t b = v % q;
  const std::size_t a = (v - (q + 1) * b) / q;
  return {a, b};
}

// Sparse polynomial in R = F_{q^2}[x, y] / (y^q + y - x^(q+1)), indexed by
// Weierstrass semigroup pole orders u in H(P_inf).
template <typename StorageType>
struct HermitianSemigroupPoly {
  // Sorted in strictly ascending order of pole order; back() is leading term.
  std::vector<std::pair<std::size_t, StorageType>> terms;

  bool IsZero() const noexcept { return terms.empty(); }

  std::size_t LeadingPoleOrder() const noexcept {
    return terms.empty() ? 0 : terms.back().first;
  }

  template <typename GaloisField>
  void Canonicalize(const GaloisField &field) {
    if (terms.empty()) {
      return;
    }
    std::sort(terms.begin(), terms.end(), [](const auto &lhs, const auto &rhs) {
      return lhs.first < rhs.first;
    });
    std::vector<std::pair<std::size_t, StorageType>> merged;
    merged.reserve(terms.size());
    for (const auto &[ord, coeff] : terms) {
      if (field.IsZero(coeff)) {
        continue;
      }
      if (!merged.empty() && merged.back().first == ord) {
        merged.back().second = field.Add(merged.back().second, coeff);
        if (field.IsZero(merged.back().second)) {
          merged.pop_back();
        }
      } else {
        merged.emplace_back(ord, coeff);
      }
    }
    terms = std::move(merged);
  }
};

// Multiplies a semigroup polynomial f by phi_u (the unique monomial of pole
// order u in H(P_inf)) in R = F_{q^2}[x, y] / (y^q + y - x^(q+1)).
template <typename GaloisField>
inline HermitianSemigroupPoly<typename GaloisField::element_type>
MulBySemigroupMonomial(
    const GaloisField &field,
    const HermitianSemigroupPoly<typename GaloisField::element_type> &f,
    std::size_t u, std::size_t q) {
  using StorageType = typename GaloisField::element_type;
  HermitianSemigroupPoly<StorageType> res;
  if (f.IsZero()) {
    return res;
  }
  const std::size_t b_u = u % q;
  const std::size_t q2_minus_1 = q * q - 1;
  res.terms.reserve(f.terms.size() * 2);
  for (const auto &[v, coeff] : f.terms) {
    const std::size_t b_v = v % q;
    res.terms.emplace_back(u + v, coeff);
    if (b_u + b_v >= q) {
      res.terms.emplace_back(u + v - q2_minus_1, field.Neg(coeff));
    }
  }
  res.Canonicalize(field);
  return res;
}

// Computes a + scalar * b in R = F_{q^2}[x, y] / (y^q + y - x^(q+1)).
template <typename GaloisField>
inline HermitianSemigroupPoly<typename GaloisField::element_type>
AddScaledSemigroupPoly(
    const GaloisField &field,
    const HermitianSemigroupPoly<typename GaloisField::element_type> &a,
    const HermitianSemigroupPoly<typename GaloisField::element_type> &b,
    const typename GaloisField::element_type &scalar) {
  using StorageType = typename GaloisField::element_type;
  if (field.IsZero(scalar) || b.IsZero()) {
    return a;
  }
  HermitianSemigroupPoly<StorageType> res;
  res.terms.reserve(a.terms.size() + b.terms.size());
  std::size_t ia = 0;
  std::size_t ib = 0;
  while (ia < a.terms.size() && ib < b.terms.size()) {
    if (a.terms[ia].first < b.terms[ib].first) {
      res.terms.push_back(a.terms[ia++]);
    } else if (b.terms[ib].first < a.terms[ia].first) {
      const StorageType prod = field.Mul(scalar, b.terms[ib].second);
      if (!field.IsZero(prod)) {
        res.terms.emplace_back(b.terms[ib].first, prod);
      }
      ++ib;
    } else {
      const StorageType sum =
          field.Add(a.terms[ia].second, field.Mul(scalar, b.terms[ib].second));
      if (!field.IsZero(sum)) {
        res.terms.emplace_back(a.terms[ia].first, sum);
      }
      ++ia;
      ++ib;
    }
  }
  while (ia < a.terms.size()) {
    res.terms.push_back(a.terms[ia++]);
  }
  while (ib < b.terms.size()) {
    const StorageType prod = field.Mul(scalar, b.terms[ib].second);
    if (!field.IsZero(prod)) {
      res.terms.emplace_back(b.terms[ib].first, prod);
    }
    ++ib;
  }
  return res;
}

// Evaluates S(phi_u * f) given a syndrome table syndromes[w] = S(phi_w).
// If exclude_s is true, skips any term with pole order == target_s (used to
// compute the predicted syndrome value in Feng-Rao majority voting).
template <typename GaloisField>
inline typename GaloisField::element_type EvaluateShiftDiscrepancy(
    const GaloisField &field,
    const HermitianSemigroupPoly<typename GaloisField::element_type> &f,
    std::size_t u, std::size_t q,
    const std::vector<typename GaloisField::element_type> &syndromes,
    bool exclude_s = false, std::size_t target_s = 0) {
  using StorageType = typename GaloisField::element_type;
  StorageType acc = field.AdditiveIdentity();
  const std::size_t b_u = u % q;
  const std::size_t q2_minus_1 = q * q - 1;
  for (const auto &[v, coeff] : f.terms) {
    const std::size_t w1 = u + v;
    if (!exclude_s || w1 != target_s) {
      acc = field.Add(acc, field.Mul(coeff, syndromes[w1]));
    }
    const std::size_t b_v = v % q;
    if (b_u + b_v >= q) {
      const std::size_t w2 = w1 - q2_minus_1;
      if (!exclude_s || w2 != target_s) {
        acc = field.Sub(acc, field.Mul(coeff, syndromes[w2]));
      }
    }
  }
  return acc;
}

// Truncated univariate power series multiplication mod t^prec.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type> TruncPolyMul(
    const GaloisField &field,
    const std::vector<typename GaloisField::element_type> &a,
    const std::vector<typename GaloisField::element_type> &b,
    std::size_t prec) {
  using StorageType = typename GaloisField::element_type;
  const StorageType zero = field.AdditiveIdentity();
  std::vector<StorageType> out(prec, zero);
  const std::size_t na = std::min(a.size(), prec);
  const std::size_t nb = std::min(b.size(), prec);
  for (std::size_t i = 0; i < na; ++i) {
    if (field.IsZero(a[i])) {
      continue;
    }
    const std::size_t j_max = std::min(nb, prec - i);
    for (std::size_t j = 0; j < j_max; ++j) {
      if (!field.IsZero(b[j])) {
        out[i + j] = field.Add(out[i + j], field.Mul(a[i], b[j]));
      }
    }
  }
  return out;
}

// Computes the local power series u(t) = y - y0 mod t^prec satisfying
// u + u^q = x0^q * t + x0 * t^q + t^(q+1) on H_q at P = (x0, y0), where
// t = x - x0 is the local uniformizer.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type> LocalYShiftSeries(
    const GaloisField &field, const typename GaloisField::element_type &x0,
    std::size_t q, std::size_t prec) {
  using StorageType = typename GaloisField::element_type;
  const StorageType zero = field.AdditiveIdentity();
  const StorageType one = field.MultiplicativeIdentity();
  std::vector<StorageType> u(prec, zero);
  if (prec <= 1) {
    return u;
  }
  const StorageType x0_q = field.Pow(x0, q);

  // Iterate u <- (x0^q * t + x0 * t^q + t^(q+1) - u^q) mod t^prec.
  // Since u has valuation >= 1, u^q multiplies the valuation by q each step.
  std::size_t valid_order = 1;
  while (valid_order < prec) {
    std::vector<StorageType> next(prec, zero);
    if (1 < prec) {
      next[1] = x0_q;
    }
    if (q < prec) {
      next[q] = field.Add(next[q], x0);
    }
    if (q + 1 < prec) {
      next[q + 1] = field.Add(next[q + 1], one);
    }
    for (std::size_t k = 1; k * q < prec; ++k) {
      if (!field.IsZero(u[k])) {
        const StorageType uk_q = field.Pow(u[k], q);
        next[k * q] = field.Sub(next[k * q], uk_q);
      }
    }
    u = std::move(next);
    valid_order = (valid_order > prec / q) ? prec : (valid_order * q);
  }
  return u;
}

// Number of monomials in L(D P_inf) on H_q (pole order <= D, 0 <= j < q).
inline std::size_t CountHermitianMonomials(std::size_t d,
                                           std::size_t q) noexcept {
  std::size_t count = 0;
  for (std::size_t j = 0; j < q; ++j) {
    const std::size_t base = (q + 1) * j;
    if (base > d) {
      break;
    }
    count += (d - base) / q + 1;
  }
  return count;
}

// Total number of unknowns in Q(P, z) = sum_{s=0}^L q_s(P) z^s with
// q_s in L((D - s * m) P_inf).
inline std::size_t CountHermitianInterpolationUnknowns(std::size_t d,
                                                       std::size_t m,
                                                       std::size_t q,
                                                       std::size_t list_size) {
  std::size_t total = 0;
  for (std::size_t s = 0; s <= list_size; ++s) {
    if (s * m > d) {
      break;
    }
    total += CountHermitianMonomials(d - s * m, q);
  }
  return total;
}

}  // namespace detail

// ============================================================================
// One-Point Hermitian Algebraic-Geometry Code C_L(D, m P_inf) over F_{q^2}
// ============================================================================
//
// Defined on the Hermitian curve H_q : y^q + y = x^(q+1) over F_{q^2}, which
// has genus g = q(q - 1)/2 and achieves the Hasse-Weil bound with n = q^3
// finite rational points.
//
// Dual code (Stichtenoth 1988):
//   C_L(D, m P_inf)^perp = C_L(D, (n + 2g - 2 - m) P_inf).
template <typename GaloisField>
class HermitianCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;
  using matrix_type = GFMatrix<GaloisField>;
  using point_type = HermitianPoint<GaloisField>;

  HermitianCode(std::shared_ptr<GaloisField> field, std::size_t pole_order_m)
      : Base(ValidateAndReturnField(field), ComputeLength(field), 0,
             Metric::HAMMING, "GeneratorMatrix", "HermitianFengRao"),
        q_(ComputeQ(this->field_)),
        genus_(q_ * (q_ - 1) / 2),
        pole_order_m_(pole_order_m) {
    const std::size_t max_m = this->length_ + 2 * genus_ - 1;
    if (pole_order_m_ >= max_m) {
      throw std::invalid_argument(
          "HermitianCode: pole order m must satisfy 0 <= m < q^3 + q^2 - q - "
          "1");
    }
    points_ = std::make_shared<std::vector<point_type>>(
        EnumerateRationalPoints(this->field_, q_));
    BuildMatrices();
    RegisterEncodersAndDecoders();
  }

  HermitianCode(const HermitianCode &other)
      : Base(other),
        q_(other.q_),
        genus_(other.genus_),
        pole_order_m_(other.pole_order_m_),
        points_(other.points_),
        basis_monomials_(other.basis_monomials_),
        dual_basis_monomials_(other.dual_basis_monomials_),
        generator_(other.generator_),
        parity_check_(other.parity_check_) {
    RegisterEncodersAndDecoders();
  }

  HermitianCode &operator=(const HermitianCode &other) {
    if (this != &other) {
      Base::operator=(other);
      q_ = other.q_;
      genus_ = other.genus_;
      pole_order_m_ = other.pole_order_m_;
      points_ = other.points_;
      basis_monomials_ = other.basis_monomials_;
      dual_basis_monomials_ = other.dual_basis_monomials_;
      generator_ = other.generator_;
      parity_check_ = other.parity_check_;
      RegisterEncodersAndDecoders();
    }
    return *this;
  }

  HermitianCode(HermitianCode &&other) noexcept
      : Base(std::move(other)),
        q_(other.q_),
        genus_(other.genus_),
        pole_order_m_(other.pole_order_m_),
        points_(std::move(other.points_)),
        basis_monomials_(std::move(other.basis_monomials_)),
        dual_basis_monomials_(std::move(other.dual_basis_monomials_)),
        generator_(std::move(other.generator_)),
        parity_check_(std::move(other.parity_check_)) {
    RegisterEncodersAndDecoders();
  }

  HermitianCode &operator=(HermitianCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      q_ = other.q_;
      genus_ = other.genus_;
      pole_order_m_ = other.pole_order_m_;
      points_ = std::move(other.points_);
      basis_monomials_ = std::move(other.basis_monomials_);
      dual_basis_monomials_ = std::move(other.dual_basis_monomials_);
      generator_ = std::move(other.generator_);
      parity_check_ = std::move(other.parity_check_);
      RegisterEncodersAndDecoders();
    }
    return *this;
  }

  ~HermitianCode() override = default;

  std::size_t Q() const noexcept { return q_; }

  std::size_t Genus() const noexcept { return genus_; }

  std::size_t PoleOrderM() const noexcept { return pole_order_m_; }

  std::size_t DualPoleOrderM() const noexcept {
    return this->length_ + 2 * genus_ - 2 - pole_order_m_;
  }

  // Goppa designed minimum distance d* = n - m (for m < n).
  std::size_t DesignedDistance() const noexcept {
    return (pole_order_m_ < this->length_) ? (this->length_ - pole_order_m_)
                                           : 1;
  }

  std::size_t MinimumDistance() const override { return DesignedDistance(); }

  const std::vector<point_type> &RationalPoints() const noexcept {
    return *points_;
  }

  const std::vector<HermitianMonomial> &BasisMonomials() const noexcept {
    return basis_monomials_;
  }

  const std::vector<HermitianMonomial> &DualBasisMonomials() const noexcept {
    return dual_basis_monomials_;
  }

  // Exact O(1)-formula dual: C_L(D, m P_inf)^perp = C_L(D, (n + 2g - 2 - m)
  // P_inf), sharing precomputed rational points.
  HermitianCode DualCode() const {
    return HermitianCode(this->field_, q_, genus_, DualPoleOrderM(), points_,
                         dual_basis_monomials_, basis_monomials_, parity_check_,
                         generator_);
  }

  const matrix_type &GeneratorMatrix() const override { return generator_; }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_;
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) +
           ", >= " + std::to_string(DesignedDistance()) +
           "] Hermitian AG code H_" + std::to_string(q_) +
           "(m=" + std::to_string(pole_order_m_) + ") over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 private:
  // Private O(1) dual constructor reusing precomputed points and matrices.
  HermitianCode(std::shared_ptr<GaloisField> field, std::size_t q,
                std::size_t genus, std::size_t pole_order_m,
                std::shared_ptr<const std::vector<point_type>> points,
                std::vector<HermitianMonomial> basis_monomials,
                std::vector<HermitianMonomial> dual_basis_monomials,
                matrix_type generator, matrix_type parity_check)
      : Base(std::move(field), q * q * q, generator.Shape(0), Metric::HAMMING,
             "GeneratorMatrix", "HermitianFengRao"),
        q_(q),
        genus_(genus),
        pole_order_m_(pole_order_m),
        points_(std::move(points)),
        basis_monomials_(std::move(basis_monomials)),
        dual_basis_monomials_(std::move(dual_basis_monomials)),
        generator_(std::move(generator)),
        parity_check_(std::move(parity_check)) {
    RegisterEncodersAndDecoders();
  }

  static std::shared_ptr<GaloisField> ValidateAndReturnField(
      const std::shared_ptr<GaloisField> &field) {
    if (!field) {
      throw std::invalid_argument(
          "HermitianCode: field pointer cannot be null");
    }
    const uint32_t deg = field->Degree();
    if (deg < 2 || (deg % 2) != 0) {
      throw std::invalid_argument(
          "HermitianCode: field degree must be even (>= 2) for F_{q^2}");
    }
    return field;
  }

  static std::size_t ComputeQ(const std::shared_ptr<GaloisField> &field) {
    const uint64_t q_val = utils::SafeIntegerPower(
        field->Characteristic(), static_cast<uint64_t>(field->Degree() / 2));
    return static_cast<std::size_t>(q_val);
  }

  static std::size_t ComputeLength(const std::shared_ptr<GaloisField> &field) {
    const std::size_t q = ComputeQ(field);
    return static_cast<std::size_t>(utils::SafeIntegerPower(q, 3));
  }

  // O(q^3) enumeration of all q^3 finite rational points (x, y) satisfying
  // y^q + y = x^(q+1) over F_{q^2}, by precomputing the q solutions to
  // y^q + y = b for each b in F_q in O(q^2).
  static std::vector<point_type> EnumerateRationalPoints(
      const std::shared_ptr<GaloisField> &field, std::size_t q) {
    const std::vector<storage_type> elems =
        detail::EnumerateHermitianFieldElements(field);

    std::unordered_map<storage_type, std::vector<storage_type>> trace_preimages;
    trace_preimages.reserve(q * 2);
    for (const auto &y : elems) {
      const storage_type b = field->Add(field->Pow(y, q), y);
      trace_preimages[b].push_back(y);
    }

    std::vector<point_type> pts;
    pts.reserve(q * q * q);
    for (const auto &x : elems) {
      const storage_type rhs = field->Mul(field->Pow(x, q), x);
      const auto it = trace_preimages.find(rhs);
      if (it == trace_preimages.end()) {
        throw std::runtime_error(
            "HermitianCode: internal error enumerating curve points");
      }
      for (const auto &y : it->second) {
        pts.push_back(point_type{x, y});
      }
    }
    return pts;
  }

  // Builds the evaluation matrix and linearly independent monomial basis for
  // C_L(D, m_val P_inf). Removes linearly dependent monomials if m_val >= q^3.
  std::pair<matrix_type, std::vector<HermitianMonomial>> BuildEvaluationSpace(
      std::size_t m_val) const {
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    std::vector<HermitianMonomial> monomials;
    for (std::size_t j = 0; j < q_; ++j) {
      const std::size_t base = (q_ + 1) * j;
      if (base > m_val) {
        break;
      }
      const std::size_t max_i = (m_val - base) / q_;
      for (std::size_t i = 0; i <= max_i; ++i) {
        monomials.push_back(HermitianMonomial{i, j, q_ * i + base});
      }
    }
    std::sort(monomials.begin(), monomials.end(),
              [](const HermitianMonomial &lhs, const HermitianMonomial &rhs) {
                return lhs.pole_order < rhs.pole_order;
              });

    const std::size_t n = this->length_;
    const std::size_t num_mon = monomials.size();
    if (num_mon == 0) {
      return {matrix_type::Zeros(this->field_, {0, n}), {}};
    }

    std::size_t max_i = 0;
    for (const auto &mon : monomials) {
      if (mon.i > max_i) {
        max_i = mon.i;
      }
    }

    // Precompute x_l^i and y_l^j tables for O(1) evaluation per entry.
    const storage_type one = this->field_->MultiplicativeIdentity();
    std::vector<std::vector<storage_type>> x_pows(
        n, std::vector<storage_type>(max_i + 1, one));
    std::vector<std::vector<storage_type>> y_pows(
        n, std::vector<storage_type>(q_, one));
    for (std::size_t l = 0; l < n; ++l) {
      const auto &pt = (*points_)[l];
      for (std::size_t i = 1; i <= max_i; ++i) {
        x_pows[l][i] = this->field_->Mul(x_pows[l][i - 1], pt.x);
      }
      for (std::size_t j = 1; j < q_; ++j) {
        y_pows[l][j] = this->field_->Mul(y_pows[l][j - 1], pt.y);
      }
    }

    RawArray eval_raw = RawArray::from_shape(ShapeType{num_mon, n});
    storage_type *data = eval_raw.data();
    for (std::size_t r = 0; r < num_mon; ++r) {
      const auto &mon = monomials[r];
      storage_type *row_ptr = data + r * n;
      for (std::size_t l = 0; l < n; ++l) {
        row_ptr[l] = this->field_->Mul(x_pows[l][mon.i], y_pows[l][mon.j]);
      }
    }
    matrix_type eval_mat =
        matrix_type::FromRaw(this->field_, std::move(eval_raw));

    if (m_val < n) {
      return {std::move(eval_mat), std::move(monomials)};
    }

    // For m_val >= q^3, remove linearly dependent monomials.
    const matrix_type row_space = linalg::RowSpace(eval_mat);
    const std::size_t rank = row_space.Shape(0);
    if (rank == num_mon) {
      return {std::move(eval_mat), std::move(monomials)};
    }

    // Select exact linearly independent pivot rows to keep basis_monomials_
    // aligned with the generator rows.
    const auto sys_t = linalg::SystematicForm(linalg::Transpose(eval_mat));
    const auto &keep_rows = sys_t.pivot_cols;

    RawArray reduced_raw = RawArray::from_shape(ShapeType{keep_rows.size(), n});
    std::vector<HermitianMonomial> reduced_monomials;
    reduced_monomials.reserve(keep_rows.size());
    const storage_type *src = eval_mat.RawData().data();
    storage_type *dst = reduced_raw.data();
    for (std::size_t r = 0; r < keep_rows.size(); ++r) {
      const std::size_t orig = keep_rows[r];
      reduced_monomials.push_back(monomials[orig]);
      std::copy(src + orig * n, src + (orig + 1) * n, dst + r * n);
    }
    return {matrix_type::FromRaw(this->field_, std::move(reduced_raw)),
            std::move(reduced_monomials)};
  }

  void BuildMatrices() {
    auto [g_mat, g_mons] = BuildEvaluationSpace(pole_order_m_);
    generator_ = std::move(g_mat);
    basis_monomials_ = std::move(g_mons);
    this->dimension_ = generator_.Shape(0);

    const std::size_t m_perp = DualPoleOrderM();
    auto [h_mat, h_mons] = BuildEvaluationSpace(m_perp);
    parity_check_ = std::move(h_mat);
    dual_basis_monomials_ = std::move(h_mons);
  }

  void RegisterEncodersAndDecoders() {
    this->RegisterEncoder(
        "GeneratorMatrix", [](const AbstractCode<GaloisField> &c) {
          return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(
              static_cast<const AbstractLinearCode<GaloisField> &>(c));
        });
    this->RegisterDecoder(
        "HermitianFengRao", [](const AbstractCode<GaloisField> &c) {
          return std::make_shared<HermitianFengRaoDecoder<GaloisField>>(c);
        });
    this->RegisterDecoder(
        "HermitianSudan", [](const AbstractCode<GaloisField> &c) {
          return std::make_shared<HermitianSudanListDecoder<GaloisField>>(c, 1);
        });
  }

  std::size_t q_ = 0;
  std::size_t genus_ = 0;
  std::size_t pole_order_m_ = 0;
  std::shared_ptr<const std::vector<point_type>> points_;
  std::vector<HermitianMonomial> basis_monomials_;
  std::vector<HermitianMonomial> dual_basis_monomials_;
  matrix_type generator_;
  matrix_type parity_check_;
};

// ============================================================================
// HermitianFengRaoDecoder: 2D BMS + Feng-Rao Majority Voting
// ============================================================================
//
// Corrects up to the full designed error-correcting radius
//   t = floor((d* - 1) / 2),  d* = n - m,
// strictly surpassing basic AG decoding floor((d* - 1 - g) / 2).
template <typename GaloisField>
class HermitianFengRaoDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using matrix_type = GFMatrix<GaloisField>;

  explicit HermitianFengRaoDecoder(const AbstractCode<GaloisField> &code)
      : Base(code) {
    hermitian_code_ = dynamic_cast<const HermitianCode<GaloisField> *>(&code);
    if (!hermitian_code_) {
      throw std::invalid_argument(
          "HermitianFengRaoDecoder requires a HermitianCode");
    }
    const std::size_t d_star = hermitian_code_->DesignedDistance();
    radius_ = (d_star > 0) ? ((d_star - 1) / 2) : 0;
  }

  std::string Name() const override { return "HermitianFengRao"; }

  std::size_t DecodingRadius() const override { return radius_; }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    const auto field = hermitian_code_->Field();
    const std::size_t n = hermitian_code_->Length();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "HermitianFengRaoDecoder: received word field mismatch");
    }
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "HermitianFengRaoDecoder: received word must be 1D of length n");
    }

    // Fast path: check if received_word is already a valid codeword.
    if (hermitian_code_->Contains(received_word)) {
      auto msg = hermitian_code_->Unencode(received_word);
      return result_type::Success(received_word, std::move(msg), 0, 0);
    }
    if (radius_ == 0) {
      return result_type::Failure(
          "HermitianFengRao: code has decoding radius 0");
    }

    const std::size_t q = hermitian_code_->Q();
    const std::size_t g = hermitian_code_->Genus();
    const std::size_t m_perp = hermitian_code_->DualPoleOrderM();
    const auto &points = hermitian_code_->RationalPoints();
    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();

    // Run BMS + Feng-Rao majority voting up to s_max so the error-locator
    // Groebner basis stabilizes for all error weights <= radius_.
    const std::size_t s_max = m_perp + 4 * g + 2 * q + 2;
    std::vector<storage_type> syndromes(s_max + 1, zero);

    // Precompute x_l^a and y_l^b for all known syndromes rho(a, b) <= m_perp.
    const std::size_t max_a_known = m_perp / q;
    std::vector<std::vector<storage_type>> x_pows(
        n, std::vector<storage_type>(max_a_known + 1, one));
    std::vector<std::vector<storage_type>> y_pows(
        n, std::vector<storage_type>(q, one));
    for (std::size_t l = 0; l < n; ++l) {
      for (std::size_t a = 1; a <= max_a_known; ++a) {
        x_pows[l][a] = field->Mul(x_pows[l][a - 1], points[l].x);
      }
      for (std::size_t b = 1; b < q; ++b) {
        y_pows[l][b] = field->Mul(y_pows[l][b - 1], points[l].y);
      }
    }

    const auto recv_span = received_word.RawSpan();
    for (std::size_t s = 0; s <= m_perp; ++s) {
      if (!detail::IsInHermitianSemigroup(s, q)) {
        continue;
      }
      const auto [a, b] = detail::DecomposePoleOrder(s, q);
      storage_type sum = zero;
      for (std::size_t l = 0; l < n; ++l) {
        if (field->IsZero(recv_span[l])) {
          continue;
        }
        const storage_type eval = field->Mul(x_pows[l][a], y_pows[l][b]);
        sum = field->Add(sum, field->Mul(recv_span[l], eval));
      }
      syndromes[s] = sum;
    }

    // Initialize BMS minimal polynomial set F = {1} and auxiliary set G = {}.
    using Poly = detail::HermitianSemigroupPoly<storage_type>;
    struct AuxPoly {
      Poly poly;
      std::size_t span = 0;
      storage_type discrepancy{};
    };

    std::vector<Poly> basis_f(1);
    basis_f[0].terms = {{0, one}};
    std::vector<AuxPoly> aux_g;

    // delta[u] == true iff u in H(P_inf) is in the current delta set Delta(F).
    std::vector<bool> in_delta(s_max + 1, false);

    auto is_covered_by_f = [&](std::size_t u) -> bool {
      for (const auto &f : basis_f) {
        const std::size_t lead = f.LeadingPoleOrder();
        if (u >= lead && detail::IsInHermitianSemigroup(u - lead, q)) {
          return true;
        }
      }
      return false;
    };

    for (std::size_t s = 0; s <= s_max; ++s) {
      if (!detail::IsInHermitianSemigroup(s, q)) {
        continue;
      }

      // Feng-Rao majority voting for unknown syndromes s > m_perp.
      if (s > m_perp) {
        // Precompute predicted value S_hat(f) for each f in F with s - lead in
        // H(P_inf).
        std::vector<bool> can_predict(basis_f.size(), false);
        std::vector<storage_type> predictions(basis_f.size(), zero);
        for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
          const std::size_t lead = basis_f[idx].LeadingPoleOrder();
          if (s >= lead && detail::IsInHermitianSemigroup(s - lead, q)) {
            can_predict[idx] = true;
            const storage_type partial = detail::EvaluateShiftDiscrepancy(
                *field, basis_f[idx], s - lead, q, syndromes,
                /*exclude_s=*/true, s);
            predictions[idx] = field->Neg(partial);
          }
        }

        // Count votes over all candidate pairs (i, j) with i + j = s,
        // i in H \ Delta(F), j in H \ Delta(F).
        std::unordered_map<storage_type, std::size_t> votes;
        for (std::size_t i = 0; i <= s; ++i) {
          if (!detail::IsInHermitianSemigroup(i, q) || in_delta[i]) {
            continue;
          }
          const std::size_t j = s - i;
          if (!detail::IsInHermitianSemigroup(j, q) || in_delta[j]) {
            continue;
          }
          for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
            if (!can_predict[idx]) {
              continue;
            }
            const std::size_t lead = basis_f[idx].LeadingPoleOrder();
            if (i >= lead && detail::IsInHermitianSemigroup(i - lead, q)) {
              ++votes[predictions[idx]];
              break;
            }
          }
        }

        storage_type best_val = zero;
        std::size_t best_votes = 0;
        for (const auto &[val, count] : votes) {
          if (count > best_votes) {
            best_votes = count;
            best_val = val;
          }
        }
        syndromes[s] = best_val;
      }

      // Compute discrepancies of all f in F at step s.
      std::vector<storage_type> disc(basis_f.size(), zero);
      std::vector<bool> failed(basis_f.size(), false);
      bool any_failed = false;
      for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
        const std::size_t lead = basis_f[idx].LeadingPoleOrder();
        if (s >= lead && detail::IsInHermitianSemigroup(s - lead, q)) {
          const storage_type d = detail::EvaluateShiftDiscrepancy(
              *field, basis_f[idx], s - lead, q, syndromes);
          disc[idx] = d;
          if (!field->IsZero(d)) {
            failed[idx] = true;
            any_failed = true;
          }
        }
      }

      if (!any_failed) {
        continue;
      }

      // Update Delta_new = Delta(F) U {u in H : u <=_H (s - lead(f)), f failed}
      std::vector<bool> new_delta = in_delta;
      for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
        if (!failed[idx]) {
          continue;
        }
        const std::size_t span = s - basis_f[idx].LeadingPoleOrder();
        for (std::size_t u = 0; u <= span; ++u) {
          if (detail::IsInHermitianSemigroup(u, q) &&
              detail::IsInHermitianSemigroup(span - u, q)) {
            new_delta[u] = true;
          }
        }
      }

      // Find minimal generators of H \ Delta_new (at most one per residue class
      // mod q).
      std::vector<std::size_t> class_min(q, 0);
      for (std::size_t r = 0; r < q; ++r) {
        std::size_t cand = (q + 1) * r;
        while (cand <= s_max && new_delta[cand]) {
          cand += q;
        }
        class_min[r] = cand;
      }

      std::vector<std::size_t> min_gens;
      min_gens.reserve(q);
      for (std::size_t r = 0; r < q; ++r) {
        bool is_minimal = true;
        for (std::size_t other = 0; other < q; ++other) {
          if (other == r) {
            continue;
          }
          if (class_min[r] >= class_min[other] &&
              detail::IsInHermitianSemigroup(class_min[r] - class_min[other],
                                             q)) {
            is_minimal = false;
            break;
          }
        }
        if (is_minimal) {
          min_gens.push_back(class_min[r]);
        }
      }
      std::sort(min_gens.begin(), min_gens.end());

      // Construct updated minimal polynomial basis F_new.
      std::vector<Poly> next_f;
      next_f.reserve(min_gens.size());
      for (std::size_t m_gen : min_gens) {
        // Case A: some non-failing f in F divides m_gen.
        bool built = false;
        for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
          if (failed[idx]) {
            continue;
          }
          const std::size_t lead = basis_f[idx].LeadingPoleOrder();
          if (m_gen >= lead &&
              detail::IsInHermitianSemigroup(m_gen - lead, q)) {
            next_f.push_back(detail::MulBySemigroupMonomial(
                *field, basis_f[idx], m_gen - lead, q));
            built = true;
            break;
          }
        }
        if (built) {
          continue;
        }

        // Case B: find a failing f in F dividing m_gen.
        std::size_t f_idx = basis_f.size();
        for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
          if (!failed[idx]) {
            continue;
          }
          const std::size_t lead = basis_f[idx].LeadingPoleOrder();
          if (m_gen >= lead &&
              detail::IsInHermitianSemigroup(m_gen - lead, q)) {
            f_idx = idx;
            break;
          }
        }
        if (f_idx == basis_f.size()) {
          return result_type::Failure(
              "HermitianFengRao: BMS basis construction failed");
        }

        const Poly shifted_f = detail::MulBySemigroupMonomial(
            *field, basis_f[f_idx], m_gen - basis_f[f_idx].LeadingPoleOrder(),
            q);

        // Subcase B1: s - m_gen not in H(P_inf) => no discrepancy at s.
        if (s < m_gen || !detail::IsInHermitianSemigroup(s - m_gen, q)) {
          next_f.push_back(shifted_f);
          continue;
        }

        // Subcase B2: cancel discrepancy using g in G with (s - m_gen) <=_H
        // span(g).
        const std::size_t need = s - m_gen;
        std::size_t g_idx = aux_g.size();
        for (std::size_t gi = 0; gi < aux_g.size(); ++gi) {
          if (aux_g[gi].span >= need &&
              detail::IsInHermitianSemigroup(aux_g[gi].span - need, q)) {
            g_idx = gi;
            break;
          }
        }
        if (g_idx == aux_g.size()) {
          return result_type::Failure(
              "HermitianFengRao: missing auxiliary polynomial in BMS update");
        }
        const std::size_t u_g = aux_g[g_idx].span - need;
        const Poly shifted_g =
            detail::MulBySemigroupMonomial(*field, aux_g[g_idx].poly, u_g, q);
        const storage_type factor =
            field->Neg(field->Div(disc[f_idx], aux_g[g_idx].discrepancy));
        next_f.push_back(detail::AddScaledSemigroupPoly(*field, shifted_f,
                                                        shifted_g, factor));
      }

      // Update auxiliary set G with failing polynomials from F and prune to
      // maximal spans in Delta_new.
      for (std::size_t idx = 0; idx < basis_f.size(); ++idx) {
        if (failed[idx]) {
          aux_g.push_back(AuxPoly{
              basis_f[idx], s - basis_f[idx].LeadingPoleOrder(), disc[idx]});
        }
      }
      std::vector<AuxPoly> pruned_g;
      pruned_g.reserve(aux_g.size());
      for (std::size_t i = 0; i < aux_g.size(); ++i) {
        bool dominated = false;
        for (std::size_t j = 0; j < aux_g.size(); ++j) {
          if (i == j) {
            continue;
          }
          if (aux_g[j].span == aux_g[i].span && j < i) {
            dominated = true;
            break;
          }
          if (aux_g[j].span > aux_g[i].span &&
              detail::IsInHermitianSemigroup(aux_g[j].span - aux_g[i].span,
                                             q)) {
            dominated = true;
            break;
          }
        }
        if (!dominated) {
          pruned_g.push_back(aux_g[i]);
        }
      }
      aux_g = std::move(pruned_g);
      basis_f = std::move(next_f);
      in_delta = std::move(new_delta);
      (void)is_covered_by_f;
    }

    // Evaluate Groebner basis F on all rational points to locate errors.
    std::size_t max_a_eval = max_a_known;
    for (const auto &f : basis_f) {
      for (const auto &[ord, _] : f.terms) {
        const auto [a, b] = detail::DecomposePoleOrder(ord, q);
        (void)b;
        if (a > max_a_eval) {
          max_a_eval = a;
        }
      }
    }
    if (max_a_eval > max_a_known) {
      for (std::size_t l = 0; l < n; ++l) {
        x_pows[l].resize(max_a_eval + 1, one);
        for (std::size_t a = max_a_known + 1; a <= max_a_eval; ++a) {
          x_pows[l][a] = field->Mul(x_pows[l][a - 1], points[l].x);
        }
      }
    }

    std::vector<std::size_t> error_positions;
    error_positions.reserve(radius_ + 1);
    for (std::size_t l = 0; l < n; ++l) {
      bool is_common_root = true;
      for (const auto &f : basis_f) {
        storage_type val = zero;
        for (const auto &[ord, coeff] : f.terms) {
          const auto [a, b] = detail::DecomposePoleOrder(ord, q);
          const storage_type term =
              field->Mul(coeff, field->Mul(x_pows[l][a], y_pows[l][b]));
          val = field->Add(val, term);
        }
        if (!field->IsZero(val)) {
          is_common_root = false;
          break;
        }
      }
      if (is_common_root) {
        error_positions.push_back(l);
      }
    }

    const std::size_t num_err = error_positions.size();
    if (num_err == 0 || num_err > radius_) {
      return result_type::Failure(
          "HermitianFengRao: located " + std::to_string(num_err) +
          " error positions (radius = " + std::to_string(radius_) + ")");
    }

    // Solve for error magnitudes from known parity-check syndromes H * e = H *
    // r using linalg::Solve.
    const auto &H = hermitian_code_->ParityCheckMatrix();
    const auto syn = hermitian_code_->Syndrome(received_word);
    const std::size_t h_rows = H.Shape(0);
    if (h_rows < num_err) {
      return result_type::Failure(
          "HermitianFengRao: insufficient syndromes to solve error magnitudes");
    }

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;

    // Extract the h_rows x num_err submatrix H_E on error_positions.
    RawArray he_raw = RawArray::from_shape(ShapeType{h_rows, num_err});
    const storage_type *h_data = H.RawData().data();
    storage_type *he_data = he_raw.data();
    for (std::size_t r = 0; r < h_rows; ++r) {
      const storage_type *h_row = h_data + r * n;
      storage_type *he_row = he_data + r * num_err;
      for (std::size_t c = 0; c < num_err; ++c) {
        he_row[c] = h_row[error_positions[c]];
      }
    }
    matrix_type H_E = matrix_type::FromRaw(field, std::move(he_raw));

    // Select num_err linearly independent rows of H_E.
    const auto sys_t = linalg::SystematicForm(linalg::Transpose(H_E));
    if (sys_t.pivot_cols.size() != num_err) {
      return result_type::Failure(
          "HermitianFengRao: error-location submatrix is rank deficient");
    }

    RawArray sq_raw = RawArray::from_shape(ShapeType{num_err, num_err});
    RawArray rhs_raw = RawArray::from_shape(ShapeType{num_err});
    storage_type *sq_data = sq_raw.data();
    storage_type *rhs_data = rhs_raw.data();
    const storage_type *syn_data = syn.RawData().data();
    for (std::size_t r = 0; r < num_err; ++r) {
      const std::size_t orig_r = sys_t.pivot_cols[r];
      const storage_type *src_row = he_data + orig_r * num_err;
      std::copy(src_row, src_row + num_err, sq_data + r * num_err);
      rhs_data[r] = syn_data[orig_r];
    }

    matrix_type sq_mat = matrix_type::FromRaw(field, std::move(sq_raw));
    codeword_type rhs_vec = codeword_type::FromRaw(field, std::move(rhs_raw));
    codeword_type err_vals;
    try {
      err_vals = linalg::Solve(sq_mat, rhs_vec);
    } catch (const std::runtime_error &) {
      return result_type::Failure(
          "HermitianFengRao: singular syndrome system for error magnitudes");
    }

    // Reconstruct candidate codeword c = r - e and verify.
    codeword_type candidate = received_word;
    const auto ev_span = err_vals.RawSpan();
    std::size_t actual_weight = 0;
    for (std::size_t c = 0; c < num_err; ++c) {
      if (!field->IsZero(ev_span[c])) {
        ++actual_weight;
        candidate.Raw(error_positions[c]) =
            field->Sub(candidate.Raw(error_positions[c]), ev_span[c]);
      }
    }

    if (actual_weight > radius_ || !hermitian_code_->Contains(candidate)) {
      return result_type::Failure(
          "HermitianFengRao: reconstructed word failed codeword verification");
    }

    auto msg = hermitian_code_->Unencode(candidate);
    return result_type::Success(std::move(candidate), std::move(msg),
                                actual_weight, 0);
  }

 private:
  const HermitianCode<GaloisField> *hermitian_code_ = nullptr;
  std::size_t radius_ = 0;
};

// ============================================================================
// HermitianSudanListDecoder: Sudan / Guruswami-Sudan List Decoder on H_q
// ============================================================================
//
// Interpolates a nonzero Q(P, z) = sum_{s=0}^L q_s(P) z^s vanishing with
// multiplicity mu >= 1 at every (P_l, r_l) via exact local power-series
// expansions on H_q and linalg::Nullspace, then extracts all roots f in
// L(m P_inf) via local power-series Roth-Ruckenstein root-finding.
template <typename GaloisField>
class HermitianSudanListDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using matrix_type = GFMatrix<GaloisField>;
  using polynomial_type = PolynomialDense<GaloisField>;

  explicit HermitianSudanListDecoder(const AbstractCode<GaloisField> &code,
                                     std::size_t multiplicity = 1)
      : Base(code) {
    hermitian_code_ = dynamic_cast<const HermitianCode<GaloisField> *>(&code);
    if (!hermitian_code_) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder requires a HermitianCode");
    }
    if (multiplicity == 0) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder: multiplicity must be >= 1");
    }
    const std::size_t opt_l = ComputeOptimalListSize(
        hermitian_code_->Length(), hermitian_code_->PoleOrderM(),
        hermitian_code_->Q(), multiplicity);
    InitializeParameters(multiplicity, opt_l);
  }

  HermitianSudanListDecoder(const AbstractCode<GaloisField> &code,
                            std::size_t multiplicity, std::size_t list_size)
      : Base(code) {
    hermitian_code_ = dynamic_cast<const HermitianCode<GaloisField> *>(&code);
    if (!hermitian_code_) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder requires a HermitianCode");
    }
    if (multiplicity == 0) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder: multiplicity must be >= 1");
    }
    if (list_size == 0) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder: list_size must be >= 1");
    }
    InitializeParameters(multiplicity, list_size);
  }

  // Automatically selects minimal multiplicity mu >= 1 and list size L >= 1
  // guaranteeing list decoding up to target_errors.
  static HermitianSudanListDecoder FromTargetRadius(
      const AbstractCode<GaloisField> &code, std::size_t target_errors) {
    const auto *hc = dynamic_cast<const HermitianCode<GaloisField> *>(&code);
    if (!hc) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder::FromTargetRadius requires a "
          "HermitianCode");
    }
    const std::size_t n = hc->Length();
    const std::size_t m = hc->PoleOrderM();
    const std::size_t q = hc->Q();
    if (target_errors >= n) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder::FromTargetRadius: target_errors must be "
          "< n");
    }
    const std::size_t t_agree = n - target_errors;
    if (m == 0) {
      const std::size_t l = std::max<std::size_t>(1, n / t_agree);
      return HermitianSudanListDecoder(code, 1, l);
    }
    // Check asymptotic Johnson-type feasibility: t_agree^2 > n * m.
    if (static_cast<uint64_t>(t_agree) * static_cast<uint64_t>(t_agree) <=
        static_cast<uint64_t>(n) * static_cast<uint64_t>(m)) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder::FromTargetRadius: target_errors exceeds "
          "Johnson-type list-decoding bound for this Hermitian code");
    }

    for (std::size_t mu = 1; mu <= 64; ++mu) {
      const std::size_t d_target = mu * t_agree - 1;
      const std::size_t constraints = n * mu * (mu + 1) / 2;
      const std::size_t l_max = std::max<std::size_t>(1, d_target / m);
      if (detail::CountHermitianInterpolationUnknowns(d_target, m, q, l_max) >
          constraints) {
        std::size_t lo = 1;
        std::size_t hi = l_max;
        while (lo < hi) {
          const std::size_t mid = lo + (hi - lo) / 2;
          if (detail::CountHermitianInterpolationUnknowns(d_target, m, q, mid) >
              constraints) {
            hi = mid;
          } else {
            lo = mid + 1;
          }
        }
        return HermitianSudanListDecoder(code, mu, lo);
      }
    }
    throw std::invalid_argument(
        "HermitianSudanListDecoder::FromTargetRadius: required multiplicity "
        "exceeds supported limit");
  }

  std::size_t Multiplicity() const noexcept { return multiplicity_; }

  std::size_t ListSize() const noexcept { return list_size_; }

  std::size_t WeightedDegreeBound() const noexcept { return degree_bound_; }

  std::size_t DecodingRadius() const override { return decoding_radius_; }

  std::string Name() const override { return "HermitianSudan"; }

  // Returns all valid codewords within Hamming distance <= DecodingRadius() of
  // received_word, deduplicated and sorted by ascending Hamming distance.
  std::vector<codeword_type> DecodeToList(
      const codeword_type &received_word) const {
    const auto candidates = DecodeCandidates(received_word);
    std::vector<codeword_type> out;
    out.reserve(candidates.size());
    for (const auto &cand : candidates) {
      out.push_back(cand.codeword);
    }
    return out;
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto candidates = DecodeCandidates(received_word);
    if (candidates.empty()) {
      return result_type::Failure(
          "HermitianSudan: no codeword found within decoding radius " +
          std::to_string(decoding_radius_));
    }
    auto best_cw = std::move(candidates[0].codeword);
    const std::size_t best_dist = candidates[0].distance;
    auto best_msg = hermitian_code_->Unencode(best_cw);
    auto res = result_type::Success(std::move(best_cw), std::move(best_msg),
                                    best_dist, 0);
    if (candidates.size() > 1 && candidates[1].distance == best_dist) {
      std::size_t tied = 1;
      while (tied < candidates.size() &&
             candidates[tied].distance == best_dist) {
        ++tied;
      }
      res.error_message = "HermitianSudan: " + std::to_string(tied) +
                          " equidistant codewords at Hamming distance " +
                          std::to_string(best_dist);
    }
    return res;
  }

 private:
  struct DecodedCandidate {
    codeword_type codeword;
    std::size_t distance = 0;
  };

  struct InterpMonomial {
    std::size_t s = 0;  // z-degree
    std::size_t i = 0;  // x-degree
    std::size_t j = 0;  // y-degree
    std::size_t pole_order = 0;
  };

  struct BivariateSeries {
    // z_coeffs[s][d] is coefficient of t^d z^s.
    std::vector<std::vector<storage_type>> z_coeffs;
    std::size_t prec = 0;

    bool IsZero(const GaloisField &field) const {
      for (const auto &row : z_coeffs) {
        for (const auto &c : row) {
          if (!field.IsZero(c)) {
            return false;
          }
        }
      }
      return true;
    }

    int ZDegree(const GaloisField &field) const {
      for (std::size_t s = z_coeffs.size(); s > 0; --s) {
        for (const auto &c : z_coeffs[s - 1]) {
          if (!field.IsZero(c)) {
            return static_cast<int>(s - 1);
          }
        }
      }
      return -1;
    }
  };

  static std::size_t ComputeDegreeBound(std::size_t n, std::size_t m,
                                        std::size_t q, std::size_t mu,
                                        std::size_t list_size) {
    const std::size_t constraints = n * mu * (mu + 1) / 2;
    std::size_t lo = 0;
    std::size_t hi = 1;
    while (detail::CountHermitianInterpolationUnknowns(hi, m, q, list_size) <=
           constraints) {
      hi *= 2;
    }
    while (lo < hi) {
      const std::size_t mid = lo + (hi - lo) / 2;
      if (detail::CountHermitianInterpolationUnknowns(mid, m, q, list_size) >
          constraints) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  static std::size_t ComputeOptimalListSize(std::size_t n, std::size_t m,
                                            std::size_t q, std::size_t mu) {
    const std::size_t constraints = n * mu * (mu + 1) / 2;
    if (m == 0) {
      return std::max<std::size_t>(1, (n * (mu + 1)) / 2);
    }
    std::size_t lo = 0;
    std::size_t hi = 1;
    while (detail::CountHermitianInterpolationUnknowns(hi, m, q, hi / m) <=
           constraints) {
      hi *= 2;
    }
    while (lo < hi) {
      const std::size_t mid = lo + (hi - lo) / 2;
      if (detail::CountHermitianInterpolationUnknowns(mid, m, q, mid / m) >
          constraints) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    const std::size_t d_opt = lo;
    const std::size_t min_agree = d_opt / mu + 1;
    const std::size_t d_target = (min_agree > n) ? d_opt : (mu * min_agree - 1);
    const std::size_t l_max = std::max<std::size_t>(1, d_target / m);
    std::size_t l_lo = 1;
    std::size_t l_hi = l_max;
    while (l_lo < l_hi) {
      const std::size_t mid = l_lo + (l_hi - l_lo) / 2;
      if (detail::CountHermitianInterpolationUnknowns(d_target, m, q, mid) >
          constraints) {
        l_hi = mid;
      } else {
        l_lo = mid + 1;
      }
    }
    return l_lo;
  }

  void InitializeParameters(std::size_t multiplicity, std::size_t list_size) {
    const std::size_t n = hermitian_code_->Length();
    const std::size_t m = hermitian_code_->PoleOrderM();
    const std::size_t q = hermitian_code_->Q();
    multiplicity_ = multiplicity;
    list_size_ = list_size;
    degree_bound_ = ComputeDegreeBound(n, m, q, multiplicity_, list_size_);
    const std::size_t min_agree = degree_bound_ / multiplicity_ + 1;
    decoding_radius_ = (min_agree > n) ? 0 : (n - min_agree);

    // Precompute interpolation monomials (s, i, j) with
    // q*i + (q+1)*j + s*m <= degree_bound_, 0 <= s <= list_size_, 0 <= j < q.
    interp_monomials_.clear();
    for (std::size_t s = 0; s <= list_size_; ++s) {
      if (s * m > degree_bound_) {
        break;
      }
      const std::size_t d_s = degree_bound_ - s * m;
      for (std::size_t j = 0; j < q; ++j) {
        const std::size_t base = (q + 1) * j;
        if (base > d_s) {
          break;
        }
        const std::size_t max_i = (d_s - base) / q;
        for (std::size_t i = 0; i <= max_i; ++i) {
          interp_monomials_.push_back(InterpMonomial{s, i, j, q * i + base});
        }
      }
    }

    // Precompute the reconstruction map from truncated local power series
    // f(t) mod t^(m+1) at P_0 = (0, 0) back to L(m P_inf) basis coefficients.
    PrecomputeLocalReconstruction();
  }

  void PrecomputeLocalReconstruction() {
    const auto field = hermitian_code_->Field();
    const std::size_t m = hermitian_code_->PoleOrderM();
    const std::size_t q = hermitian_code_->Q();
    const std::size_t k = hermitian_code_->Dimension();
    const auto &basis = hermitian_code_->BasisMonomials();
    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();

    if (k == 0) {
      return;
    }

    const std::size_t recon_prec = m + 1;
    const std::vector<storage_type> y0_series =
        detail::LocalYShiftSeries(*field, zero, q, recon_prec);
    std::vector<std::vector<storage_type>> y0_pows(
        q, std::vector<storage_type>(recon_prec, zero));
    y0_pows[0][0] = one;
    for (std::size_t j = 1; j < q; ++j) {
      y0_pows[j] =
          detail::TruncPolyMul(*field, y0_pows[j - 1], y0_series, recon_prec);
    }

    // Build recon_prec x k matrix T where column r is t^{i_r} * y0(t)^{j_r}.
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray t_raw = RawArray::from_shape(ShapeType{recon_prec, k});
    t_raw.fill(zero);
    storage_type *t_data = t_raw.data();
    for (std::size_t r = 0; r < k; ++r) {
      const std::size_t shift = basis[r].i;
      const auto &yj = y0_pows[basis[r].j];
      for (std::size_t d = shift; d < recon_prec; ++d) {
        t_data[d * k + r] = yj[d - shift];
      }
    }
    matrix_type T = matrix_type::FromRaw(field, std::move(t_raw));
    const auto sys_t = linalg::SystematicForm(linalg::Transpose(T));
    recon_pivot_rows_ = sys_t.pivot_cols;
    if (recon_pivot_rows_.size() == k) {
      RawArray tp_raw = RawArray::from_shape(ShapeType{k, k});
      storage_type *tp_data = tp_raw.data();
      for (std::size_t r = 0; r < k; ++r) {
        const std::size_t orig_r = recon_pivot_rows_[r];
        std::copy(t_data + orig_r * k, t_data + (orig_r + 1) * k,
                  tp_data + r * k);
      }
      recon_inv_ = linalg::Inv(matrix_type::FromRaw(field, std::move(tp_raw)));
    }
  }

  std::vector<DecodedCandidate> DecodeCandidates(
      const codeword_type &received_word) const {
    const auto field = hermitian_code_->Field();
    const std::size_t n = hermitian_code_->Length();
    const std::size_t m = hermitian_code_->PoleOrderM();
    const std::size_t q = hermitian_code_->Q();
    const std::size_t k = hermitian_code_->Dimension();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder: received word field mismatch");
    }
    if (received_word.Dimension() != 1 || received_word.Size() != n) {
      throw std::invalid_argument(
          "HermitianSudanListDecoder: received word must be 1D of length n");
    }

    if (k == 0) {
      auto zero_cw = codeword_type::Zeros(field, {n});
      const std::size_t dist =
          hermitian_code_->Distance(zero_cw, received_word);
      if (dist <= decoding_radius_) {
        return {DecodedCandidate{std::move(zero_cw), dist}};
      }
      return {};
    }

    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();
    const auto &points = hermitian_code_->RationalPoints();
    const auto recv_span = received_word.RawSpan();

    // Precompute binomial coefficients C(s, b) in F_{q^2}.
    std::vector<std::vector<storage_type>> binom(
        list_size_ + 1, std::vector<storage_type>(multiplicity_, zero));
    for (std::size_t s = 0; s <= list_size_; ++s) {
      binom[s][0] = one;
      const std::size_t lim = std::min(s, multiplicity_ - 1);
      for (std::size_t b = 1; b <= lim; ++b) {
        if (b == s) {
          binom[s][b] = one;
        } else {
          binom[s][b] = field->Add(binom[s - 1][b - 1], binom[s - 1][b]);
        }
      }
    }

    std::size_t max_i = 0;
    for (const auto &mon : interp_monomials_) {
      if (mon.i > max_i) {
        max_i = mon.i;
      }
    }

    // Build interpolation matrix M of shape (n * C(mu+1, 2)) x num_unknowns.
    const std::size_t constraints_per_pt =
        multiplicity_ * (multiplicity_ + 1) / 2;
    const std::size_t total_rows = n * constraints_per_pt;
    const std::size_t num_cols = interp_monomials_.size();

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray mat_raw = RawArray::from_shape(ShapeType{total_rows, num_cols});
    mat_raw.fill(zero);
    storage_type *mat_data = mat_raw.data();

    std::vector<std::vector<storage_type>> x_shift_pows(
        max_i + 1, std::vector<storage_type>(multiplicity_, zero));
    std::vector<std::vector<storage_type>> y_shift_pows(
        q, std::vector<storage_type>(multiplicity_, zero));
    std::vector<storage_type> r_pows(list_size_ + 1, one);

    std::size_t row_idx = 0;
    for (std::size_t l = 0; l < n; ++l) {
      const auto &pt = points[l];
      const storage_type &rl = recv_span[l];

      // (x_l + t)^i mod t^mu
      std::vector<storage_type> x_lin(multiplicity_, zero);
      x_lin[0] = pt.x;
      if (multiplicity_ > 1) {
        x_lin[1] = one;
      }
      x_shift_pows[0].assign(multiplicity_, zero);
      x_shift_pows[0][0] = one;
      for (std::size_t i = 1; i <= max_i; ++i) {
        x_shift_pows[i] = detail::TruncPolyMul(*field, x_shift_pows[i - 1],
                                               x_lin, multiplicity_);
      }

      // (y_l + u_l(t))^j mod t^mu
      std::vector<storage_type> y_ser =
          detail::LocalYShiftSeries(*field, pt.x, q, multiplicity_);
      y_ser[0] = pt.y;
      y_shift_pows[0].assign(multiplicity_, zero);
      y_shift_pows[0][0] = one;
      for (std::size_t j = 1; j < q; ++j) {
        y_shift_pows[j] = detail::TruncPolyMul(*field, y_shift_pows[j - 1],
                                               y_ser, multiplicity_);
      }

      // r_l^s
      r_pows[0] = one;
      for (std::size_t s = 1; s <= list_size_; ++s) {
        r_pows[s] = field->Mul(r_pows[s - 1], rl);
      }

      // Precompute local series for each (i, j) pair present? Even simpler:
      // directly compute coeff of t^a in (x_l + t)^i * (y_l + u_l(t))^j.
      for (std::size_t b = 0; b < multiplicity_; ++b) {
        for (std::size_t a = 0; a + b < multiplicity_; ++a) {
          storage_type *row_ptr = mat_data + row_idx * num_cols;
          for (std::size_t col = 0; col < num_cols; ++col) {
            const auto &mon = interp_monomials_[col];
            if (mon.s < b) {
              continue;
            }
            const storage_type z_wt =
                field->Mul(binom[mon.s][b], r_pows[mon.s - b]);
            if (field->IsZero(z_wt)) {
              continue;
            }
            // Coefficient of t^a in x_shift_pows[mon.i] * y_shift_pows[mon.j].
            storage_type t_coeff = zero;
            const auto &xp = x_shift_pows[mon.i];
            const auto &yp = y_shift_pows[mon.j];
            for (std::size_t u = 0; u <= a; ++u) {
              if (!field->IsZero(xp[u]) && !field->IsZero(yp[a - u])) {
                t_coeff = field->Add(t_coeff, field->Mul(xp[u], yp[a - u]));
              }
            }
            if (!field->IsZero(t_coeff)) {
              row_ptr[col] = field->Mul(z_wt, t_coeff);
            }
          }
          ++row_idx;
        }
      }
    }

    const matrix_type interp_mat =
        matrix_type::FromRaw(field, std::move(mat_raw));
    const matrix_type nullspace = linalg::Nullspace(interp_mat);
    if (nullspace.Shape(0) == 0) {
      return {};
    }

    // Expand the first nonzero nullspace solution Q(P, z) as a local power
    // series Q_tilde(t, z) at P_0 = (0, 0) to precision
    // B = degree_bound_ + list_size_ * (m + 1) + 2.
    const std::size_t series_prec = degree_bound_ + list_size_ * (m + 1) + 2;
    const std::vector<storage_type> y0_ser =
        detail::LocalYShiftSeries(*field, zero, q, series_prec);
    std::vector<std::vector<storage_type>> y0_pows(
        q, std::vector<storage_type>(series_prec, zero));
    y0_pows[0][0] = one;
    for (std::size_t j = 1; j < q; ++j) {
      y0_pows[j] =
          detail::TruncPolyMul(*field, y0_pows[j - 1], y0_ser, series_prec);
    }

    BivariateSeries q_series;
    q_series.prec = series_prec;
    q_series.z_coeffs.assign(list_size_ + 1,
                             std::vector<storage_type>(series_prec, zero));
    const storage_type *q_sol = nullspace.RawData().data();
    for (std::size_t col = 0; col < num_cols; ++col) {
      const storage_type coeff = q_sol[col];
      if (field->IsZero(coeff)) {
        continue;
      }
      const auto &mon = interp_monomials_[col];
      const std::size_t shift = mon.i;
      const auto &yj = y0_pows[mon.j];
      auto &dst = q_series.z_coeffs[mon.s];
      for (std::size_t d = shift; d < series_prec; ++d) {
        if (!field->IsZero(yj[d - shift])) {
          dst[d] = field->Add(dst[d], field->Mul(coeff, yj[d - shift]));
        }
      }
    }

    // Extract all local power series roots f(t) mod t^(m+1) via
    // Roth-Ruckenstein.
    std::vector<storage_type> prefix;
    prefix.reserve(m + 1);
    std::vector<std::vector<storage_type>> raw_series_roots;
    RothRuckensteinSeries(q_series, 0, m + 1, prefix, raw_series_roots);

    // Reconstruct f in L(m P_inf), evaluate to codeword, filter by Hamming
    // distance <= decoding_radius_, deduplicate, and sort.
    const auto &G = hermitian_code_->GeneratorMatrix();
    std::vector<DecodedCandidate> candidates;
    for (const auto &series : raw_series_roots) {
      if (recon_pivot_rows_.size() != k) {
        continue;
      }
      RawArray cp_raw = RawArray::from_shape(ShapeType{k});
      storage_type *cp_data = cp_raw.data();
      for (std::size_t r = 0; r < k; ++r) {
        cp_data[r] = series[recon_pivot_rows_[r]];
      }
      const message_type cp = message_type::FromRaw(field, std::move(cp_raw));
      const message_type msg = linalg::Dot(recon_inv_, cp);
      codeword_type cw = linalg::Dot(msg, G);
      const std::size_t dist = hermitian_code_->Distance(cw, received_word);
      if (dist > decoding_radius_) {
        continue;
      }
      bool duplicate = false;
      for (const auto &existing : candidates) {
        if (existing.codeword == cw) {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) {
        candidates.push_back(DecodedCandidate{std::move(cw), dist});
      }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const DecodedCandidate &a, const DecodedCandidate &b) {
                return a.distance < b.distance;
              });
    return candidates;
  }

  // Recursive Roth-Ruckenstein root-finding for truncated power series
  // Q(t, z) in F_{q^2}[[t]][z] up to target_depth terms of f(t).
  void RothRuckensteinSeries(
      const BivariateSeries &q_ser, std::size_t depth, std::size_t target_depth,
      std::vector<storage_type> &prefix,
      std::vector<std::vector<storage_type>> &out_roots) const {
    const auto field = hermitian_code_->Field();
    if (q_ser.prec == 0 || q_ser.ZDegree(*field) <= 0) {
      return;
    }

    const storage_type zero = field->AdditiveIdentity();

    // 1. Divide out the highest power t^r dividing Q(t, z).
    std::size_t min_ord = q_ser.prec;
    for (const auto &row : q_ser.z_coeffs) {
      for (std::size_t d = 0; d < min_ord; ++d) {
        if (!field->IsZero(row[d])) {
          min_ord = d;
          break;
        }
      }
      if (min_ord == 0) {
        break;
      }
    }
    if (min_ord >= q_ser.prec) {
      return;
    }

    const std::size_t new_prec = q_ser.prec - min_ord;
    const std::size_t num_z = q_ser.z_coeffs.size();
    BivariateSeries q_tilde;
    q_tilde.prec = new_prec;
    q_tilde.z_coeffs.assign(num_z, std::vector<storage_type>(new_prec, zero));
    for (std::size_t s = 0; s < num_z; ++s) {
      for (std::size_t d = 0; d < new_prec; ++d) {
        q_tilde.z_coeffs[s][d] = q_ser.z_coeffs[s][d + min_ord];
      }
    }

    // 2. Find all roots c in F_{q^2} of h(z) = Q_tilde(0, z).
    std::vector<storage_type> h_coeffs(num_z, zero);
    for (std::size_t s = 0; s < num_z; ++s) {
      h_coeffs[s] = q_tilde.z_coeffs[s][0];
    }
    const polynomial_type h =
        polynomial_type::FromRaw(std::move(h_coeffs), field, "z");
    if (h.IsZero() || h.Degree() <= 0) {
      return;
    }
    const auto roots = utils::FindRoots(h);

    // 3. Branch on each root c.
    for (const auto &root_elem : roots) {
      const storage_type c = root_elem.RawValue();
      prefix.push_back(c);

      if (depth + 1 == target_depth) {
        out_roots.push_back(prefix);
      } else {
        // Compute Q_next(t, z) = Q_tilde(t, t * z + c) mod t^new_prec.
        // First shift z -> z + c via Horner's rule in z, then shift coeff s by
        // t^s.
        std::vector<std::vector<storage_type>> r_polys;
        r_polys.reserve(num_z);
        r_polys.push_back(q_tilde.z_coeffs[num_z - 1]);
        for (std::size_t step = num_z - 1; step > 0; --step) {
          const std::size_t s_idx = step - 1;
          const std::size_t cur_sz = r_polys.size();
          r_polys.push_back(r_polys[cur_sz - 1]);
          for (std::size_t j = cur_sz - 1; j > 0; --j) {
            auto &dst = r_polys[j];
            const auto &prev = r_polys[j - 1];
            if (!field->IsZero(c)) {
              field->ScalarMulArray(std::span<storage_type>(dst),
                                    std::span<const storage_type>(dst), c);
            } else {
              std::fill(dst.begin(), dst.end(), zero);
            }
            field->AddArray(std::span<storage_type>(dst),
                            std::span<const storage_type>(dst),
                            std::span<const storage_type>(prev));
          }
          {
            auto &dst = r_polys[0];
            const auto &qs = q_tilde.z_coeffs[s_idx];
            if (!field->IsZero(c)) {
              field->ScalarMulArray(std::span<storage_type>(dst),
                                    std::span<const storage_type>(dst), c);
            } else {
              std::fill(dst.begin(), dst.end(), zero);
            }
            field->AddArray(std::span<storage_type>(dst),
                            std::span<const storage_type>(dst),
                            std::span<const storage_type>(qs));
          }
        }

        BivariateSeries q_next;
        q_next.prec = new_prec;
        q_next.z_coeffs.assign(num_z,
                               std::vector<storage_type>(new_prec, zero));
        for (std::size_t s = 0; s < num_z; ++s) {
          if (s >= new_prec) {
            continue;
          }
          for (std::size_t d = 0; d + s < new_prec; ++d) {
            q_next.z_coeffs[s][d + s] = r_polys[s][d];
          }
        }
        RothRuckensteinSeries(q_next, depth + 1, target_depth, prefix,
                              out_roots);
      }

      prefix.pop_back();
    }
  }

  const HermitianCode<GaloisField> *hermitian_code_ = nullptr;
  std::size_t multiplicity_ = 1;
  std::size_t list_size_ = 1;
  std::size_t degree_bound_ = 0;
  std::size_t decoding_radius_ = 0;
  std::vector<InterpMonomial> interp_monomials_;
  std::vector<std::size_t> recon_pivot_rows_;
  matrix_type recon_inv_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_HERMITIAN_HPP_
