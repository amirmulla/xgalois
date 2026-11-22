#ifndef XGALOIS_POLY_POLY_SPARSE_HPP_
#define XGALOIS_POLY_POLY_SPARSE_HPP_

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace poly {

// Sparse univariate polynomial over any xg::Ring R using Structure-of-Arrays
// (SoA) storage: strictly ascending `degrees_` paired with non-zero `coeffs_`.
template <Ring R>
class PolynomialSparse {
 private:
  struct InternalTag {};

 public:
  using ring_type = R;
  using element_type = typename R::element_type;

  PolynomialSparse() = default;

  static PolynomialSparse Zero(std::shared_ptr<R> ring, std::string var = "x") {
    ValidateRingAndVar(ring, var);
    return PolynomialSparse(InternalTag{}, std::move(ring), {}, {},
                            std::move(var));
  }

  static PolynomialSparse One(std::shared_ptr<R> ring, std::string var = "x") {
    ValidateRingAndVar(ring, var);
    const element_type one = ring->MultiplicativeIdentity();
    if (ring->IsZero(one)) {
      return Zero(std::move(ring), std::move(var));
    }
    return PolynomialSparse(InternalTag{}, std::move(ring), {0ULL}, {one},
                            std::move(var));
  }

  static PolynomialSparse Monomial(std::shared_ptr<R> ring, uint64_t deg,
                                   const element_type& coeff,
                                   std::string var = "x") {
    ValidateRingAndVar(ring, var);
    const element_type norm = ring->Add(ring->AdditiveIdentity(), coeff);
    if (ring->IsZero(norm)) {
      return Zero(std::move(ring), std::move(var));
    }
    return PolynomialSparse(InternalTag{}, std::move(ring), {deg}, {norm},
                            std::move(var));
  }

  // Constructs a canonical sparse polynomial from arbitrary (degree, coeff)
  // pairs, sorting by degree, summing duplicates, and dropping zero terms.
  static PolynomialSparse FromTerms(
      std::shared_ptr<R> ring,
      std::vector<std::pair<uint64_t, element_type>> terms,
      std::string var = "x") {
    ValidateRingAndVar(ring, var);
    if (terms.empty()) {
      return Zero(std::move(ring), std::move(var));
    }

    std::stable_sort(
        terms.begin(), terms.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    const element_type zero = ring->AdditiveIdentity();
    std::vector<uint64_t> degs;
    std::vector<element_type> coeffs;
    degs.reserve(terms.size());
    coeffs.reserve(terms.size());

    std::size_t i = 0;
    while (i < terms.size()) {
      const uint64_t d = terms[i].first;
      element_type acc = ring->Add(zero, terms[i].second);
      ++i;
      while (i < terms.size() && terms[i].first == d) {
        acc = ring->Add(acc, terms[i].second);
        ++i;
      }
      if (!ring->IsZero(acc)) {
        degs.push_back(d);
        coeffs.push_back(std::move(acc));
      }
    }

    return PolynomialSparse(InternalTag{}, std::move(ring), std::move(degs),
                            std::move(coeffs), std::move(var));
  }

  // Converts a dense polynomial over a FiniteField into sparse representation.
  static PolynomialSparse FromDense(const PolynomialDense<R>& dense)
    requires FiniteField<R>
  {
    if (!dense.HasField()) {
      throw std::invalid_argument(
          "PolynomialSparse::FromDense: dense polynomial has null field");
    }
    const auto field = dense.Field();
    const auto& raw = dense.RawCoeffs();
    std::vector<uint64_t> degs;
    std::vector<element_type> coeffs;
    for (std::size_t i = 0; i < raw.size(); ++i) {
      if (!field->IsZero(raw[i])) {
        degs.push_back(static_cast<uint64_t>(i));
        coeffs.push_back(raw[i]);
      }
    }
    return PolynomialSparse(InternalTag{}, field, std::move(degs),
                            std::move(coeffs), dense.GetVariable());
  }

  // Converts this sparse polynomial into a PolynomialDense over a FiniteField.
  PolynomialDense<R> ToDense() const
    requires FiniteField<R>
  {
    RequireRing();
    if (degrees_.empty()) {
      return PolynomialDense<R>::Zero(ring_, var_);
    }
    const uint64_t max_deg = degrees_.back();
    if (max_deg >=
        static_cast<uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::overflow_error(
          "PolynomialSparse::ToDense: polynomial degree exceeds size_t");
    }
    const std::size_t len = static_cast<std::size_t>(max_deg) + 1;
    std::vector<element_type> raw(len, ring_->AdditiveIdentity());
    for (std::size_t i = 0; i < degrees_.size(); ++i) {
      raw[static_cast<std::size_t>(degrees_[i])] = coeffs_[i];
    }
    return PolynomialDense<R>::FromRaw(std::move(raw), ring_, var_);
  }

  bool HasRing() const noexcept { return ring_ != nullptr; }

  const std::shared_ptr<R>& Ring() const noexcept { return ring_; }

  const std::string& GetVariable() const noexcept { return var_; }

  void SetVariable(std::string var) {
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    var_ = std::move(var);
  }

  // Returns -1 for the zero polynomial, otherwise the highest non-zero degree.
  int64_t Degree() const noexcept {
    if (degrees_.empty()) {
      return -1;
    }
    return static_cast<int64_t>(degrees_.back());
  }

  std::size_t NumTerms() const noexcept { return degrees_.size(); }

  bool IsZero() const noexcept { return degrees_.empty(); }

  const std::vector<uint64_t>& Degrees() const noexcept { return degrees_; }

  const std::vector<element_type>& Coeffs() const noexcept { return coeffs_; }

  element_type LeadingCoeff() const {
    RequireRing();
    if (coeffs_.empty()) {
      return ring_->AdditiveIdentity();
    }
    return coeffs_.back();
  }

  element_type Coeff(uint64_t deg) const {
    RequireRing();
    const auto it = std::lower_bound(degrees_.begin(), degrees_.end(), deg);
    if (it == degrees_.end() || *it != deg) {
      return ring_->AdditiveIdentity();
    }
    const std::size_t idx = static_cast<std::size_t>(it - degrees_.begin());
    return coeffs_[idx];
  }

  PolynomialSparse Add(const PolynomialSparse& other) const {
    ValidateSameRing(other);
    std::vector<uint64_t> degs;
    std::vector<element_type> coeffs;
    degs.reserve(degrees_.size() + other.degrees_.size());
    coeffs.reserve(coeffs_.size() + other.coeffs_.size());

    std::size_t i = 0;
    std::size_t j = 0;
    while (i < degrees_.size() && j < other.degrees_.size()) {
      if (degrees_[i] < other.degrees_[j]) {
        degs.push_back(degrees_[i]);
        coeffs.push_back(coeffs_[i]);
        ++i;
      } else if (other.degrees_[j] < degrees_[i]) {
        degs.push_back(other.degrees_[j]);
        coeffs.push_back(other.coeffs_[j]);
        ++j;
      } else {
        element_type sum = ring_->Add(coeffs_[i], other.coeffs_[j]);
        if (!ring_->IsZero(sum)) {
          degs.push_back(degrees_[i]);
          coeffs.push_back(std::move(sum));
        }
        ++i;
        ++j;
      }
    }
    while (i < degrees_.size()) {
      degs.push_back(degrees_[i]);
      coeffs.push_back(coeffs_[i]);
      ++i;
    }
    while (j < other.degrees_.size()) {
      degs.push_back(other.degrees_[j]);
      coeffs.push_back(other.coeffs_[j]);
      ++j;
    }

    return PolynomialSparse(InternalTag{}, ring_, std::move(degs),
                            std::move(coeffs), var_);
  }

  PolynomialSparse Sub(const PolynomialSparse& other) const {
    ValidateSameRing(other);
    std::vector<uint64_t> degs;
    std::vector<element_type> coeffs;
    degs.reserve(degrees_.size() + other.degrees_.size());
    coeffs.reserve(coeffs_.size() + other.coeffs_.size());

    std::size_t i = 0;
    std::size_t j = 0;
    while (i < degrees_.size() && j < other.degrees_.size()) {
      if (degrees_[i] < other.degrees_[j]) {
        degs.push_back(degrees_[i]);
        coeffs.push_back(coeffs_[i]);
        ++i;
      } else if (other.degrees_[j] < degrees_[i]) {
        degs.push_back(other.degrees_[j]);
        coeffs.push_back(ring_->Neg(other.coeffs_[j]));
        ++j;
      } else {
        element_type diff = ring_->Sub(coeffs_[i], other.coeffs_[j]);
        if (!ring_->IsZero(diff)) {
          degs.push_back(degrees_[i]);
          coeffs.push_back(std::move(diff));
        }
        ++i;
        ++j;
      }
    }
    while (i < degrees_.size()) {
      degs.push_back(degrees_[i]);
      coeffs.push_back(coeffs_[i]);
      ++i;
    }
    while (j < other.degrees_.size()) {
      degs.push_back(other.degrees_[j]);
      coeffs.push_back(ring_->Neg(other.coeffs_[j]));
      ++j;
    }

    return PolynomialSparse(InternalTag{}, ring_, std::move(degs),
                            std::move(coeffs), var_);
  }

  PolynomialSparse Neg() const {
    RequireRing();
    std::vector<element_type> coeffs;
    coeffs.reserve(coeffs_.size());
    for (const auto& c : coeffs_) {
      coeffs.push_back(ring_->Neg(c));
    }
    return PolynomialSparse(InternalTag{}, ring_, degrees_, std::move(coeffs),
                            var_);
  }

  PolynomialSparse ScalarMul(const element_type& scalar) const {
    RequireRing();
    const element_type s = ring_->Add(ring_->AdditiveIdentity(), scalar);
    if (degrees_.empty() || ring_->IsZero(s)) {
      return Zero(ring_, var_);
    }
    if (ring_->IsOne(s)) {
      return *this;
    }
    std::vector<uint64_t> degs;
    std::vector<element_type> coeffs;
    degs.reserve(degrees_.size());
    coeffs.reserve(coeffs_.size());
    for (std::size_t i = 0; i < degrees_.size(); ++i) {
      element_type prod = ring_->Mul(s, coeffs_[i]);
      if (!ring_->IsZero(prod)) {
        degs.push_back(degrees_[i]);
        coeffs.push_back(std::move(prod));
      }
    }
    return PolynomialSparse(InternalTag{}, ring_, std::move(degs),
                            std::move(coeffs), var_);
  }

  // Johnson's min-heap sparse polynomial multiplication:
  // O(nnz(a) * nnz(b) * log(min(nnz(a), nnz(b)))) time and
  // O(min(nnz(a), nnz(b))) auxiliary heap space.
  PolynomialSparse Mul(const PolynomialSparse& other) const {
    ValidateSameRing(other);
    if (degrees_.empty() || other.degrees_.empty()) {
      return Zero(ring_, var_);
    }

    // For commutative rings, index the heap over the operand with fewer terms.
    if constexpr (CommutativeRing<R>) {
      if (degrees_.size() > other.degrees_.size()) {
        return other.Mul(*this);
      }
    }

    struct HeapEntry {
      uint64_t deg;
      std::size_t a_idx;
      std::size_t b_idx;
    };
    auto min_cmp = [](const HeapEntry& lhs, const HeapEntry& rhs) noexcept {
      return lhs.deg > rhs.deg;
    };

    const std::size_t na = degrees_.size();
    const std::size_t nb = other.degrees_.size();
    std::vector<HeapEntry> heap;
    heap.reserve(na);
    for (std::size_t i = 0; i < na; ++i) {
      heap.push_back(HeapEntry{degrees_[i] + other.degrees_[0], i, 0});
    }
    std::make_heap(heap.begin(), heap.end(), min_cmp);

    const element_type zero = ring_->AdditiveIdentity();
    std::vector<uint64_t> out_degs;
    std::vector<element_type> out_coeffs;

    while (!heap.empty()) {
      const uint64_t cur_deg = heap.front().deg;
      element_type acc = zero;
      while (!heap.empty() && heap.front().deg == cur_deg) {
        std::pop_heap(heap.begin(), heap.end(), min_cmp);
        const HeapEntry top = heap.back();
        heap.pop_back();

        acc = ring_->Add(
            acc, ring_->Mul(coeffs_[top.a_idx], other.coeffs_[top.b_idx]));
        if (top.b_idx + 1 < nb) {
          heap.push_back(
              HeapEntry{degrees_[top.a_idx] + other.degrees_[top.b_idx + 1],
                        top.a_idx, top.b_idx + 1});
          std::push_heap(heap.begin(), heap.end(), min_cmp);
        }
      }
      if (!ring_->IsZero(acc)) {
        out_degs.push_back(cur_deg);
        out_coeffs.push_back(std::move(acc));
      }
    }

    return PolynomialSparse(InternalTag{}, ring_, std::move(out_degs),
                            std::move(out_coeffs), var_);
  }

  // Sparse Monagan-Pearce heap-based division by a monic polynomial over any
  // Ring R. Returns {q, r} such that a == q * monic_b + r and
  // deg(r) < deg(monic_b).
  static std::pair<PolynomialSparse, PolynomialSparse> MonicDivRem(
      const PolynomialSparse& a, const PolynomialSparse& monic_b) {
    a.ValidateSameRing(monic_b);
    if (monic_b.IsZero()) {
      throw std::invalid_argument(
          "PolynomialSparse::MonicDivRem: division by zero polynomial");
    }
    const auto& ring = a.ring_;
    if (!ring->IsOne(monic_b.LeadingCoeff())) {
      throw std::invalid_argument(
          "PolynomialSparse::MonicDivRem: divisor polynomial must be monic");
    }
    if (a.IsZero()) {
      return {Zero(ring, a.var_), Zero(ring, a.var_)};
    }
    const uint64_t db = monic_b.degrees_.back();
    if (db == 0) {
      return {a, Zero(ring, a.var_)};
    }
    if (a.degrees_.back() < db) {
      return {Zero(ring, a.var_), a};
    }

    struct DivHeapEntry {
      uint64_t deg;
      std::size_t q_idx;
      std::size_t b_idx;
    };
    auto max_cmp = [](const DivHeapEntry& lhs,
                      const DivHeapEntry& rhs) noexcept {
      return lhs.deg < rhs.deg;
    };

    const std::size_t num_b_tail = monic_b.degrees_.size() - 1;
    std::vector<DivHeapEntry> heap;
    heap.reserve(num_b_tail);
    std::vector<std::size_t> waiting_b;
    waiting_b.reserve(num_b_tail);
    for (std::size_t j = 0; j < num_b_tail; ++j) {
      waiting_b.push_back(j);
    }

    std::vector<uint64_t> q_degs_desc;
    std::vector<element_type> q_coeffs_desc;
    std::vector<uint64_t> r_degs_desc;
    std::vector<element_type> r_coeffs_desc;

    const element_type zero = ring->AdditiveIdentity();
    std::size_t a_rem = a.degrees_.size();

    while (a_rem > 0 || !heap.empty()) {
      uint64_t cur_deg = 0;
      if (a_rem > 0 && !heap.empty()) {
        cur_deg = std::max(a.degrees_[a_rem - 1], heap.front().deg);
      } else if (a_rem > 0) {
        cur_deg = a.degrees_[a_rem - 1];
      } else {
        cur_deg = heap.front().deg;
      }

      element_type coeff = zero;
      if (a_rem > 0 && a.degrees_[a_rem - 1] == cur_deg) {
        coeff = a.coeffs_[a_rem - 1];
        --a_rem;
      }

      while (!heap.empty() && heap.front().deg == cur_deg) {
        std::pop_heap(heap.begin(), heap.end(), max_cmp);
        const DivHeapEntry top = heap.back();
        heap.pop_back();

        coeff = ring->Sub(coeff, ring->Mul(q_coeffs_desc[top.q_idx],
                                           monic_b.coeffs_[top.b_idx]));
        if (top.q_idx + 1 < q_degs_desc.size()) {
          heap.push_back(DivHeapEntry{
              q_degs_desc[top.q_idx + 1] + monic_b.degrees_[top.b_idx],
              top.q_idx + 1, top.b_idx});
          std::push_heap(heap.begin(), heap.end(), max_cmp);
        } else {
          waiting_b.push_back(top.b_idx);
        }
      }

      if (ring->IsZero(coeff)) {
        continue;
      }

      if (cur_deg >= db) {
        const uint64_t q_deg = cur_deg - db;
        const std::size_t q_new_idx = q_degs_desc.size();
        q_degs_desc.push_back(q_deg);
        q_coeffs_desc.push_back(std::move(coeff));

        for (std::size_t j : waiting_b) {
          heap.push_back(
              DivHeapEntry{q_deg + monic_b.degrees_[j], q_new_idx, j});
          std::push_heap(heap.begin(), heap.end(), max_cmp);
        }
        waiting_b.clear();
      } else {
        r_degs_desc.push_back(cur_deg);
        r_coeffs_desc.push_back(std::move(coeff));
      }
    }

    std::reverse(q_degs_desc.begin(), q_degs_desc.end());
    std::reverse(q_coeffs_desc.begin(), q_coeffs_desc.end());
    std::reverse(r_degs_desc.begin(), r_degs_desc.end());
    std::reverse(r_coeffs_desc.begin(), r_coeffs_desc.end());

    return {PolynomialSparse(InternalTag{}, ring, std::move(q_degs_desc),
                             std::move(q_coeffs_desc), a.var_),
            PolynomialSparse(InternalTag{}, ring, std::move(r_degs_desc),
                             std::move(r_coeffs_desc), a.var_)};
  }

  std::pair<PolynomialSparse, PolynomialSparse> MonicDivRem(
      const PolynomialSparse& monic_b) const {
    return MonicDivRem(*this, monic_b);
  }

  // Fast sparse reduction modulo the binomial (x^n - c) in O(nnz) arithmetic
  // (e.g., pass c = ring->Neg(ring->MultiplicativeIdentity()) for x^n + 1).
  PolynomialSparse ModBinomial(uint64_t n, const element_type& c) const {
    RequireRing();
    if (n == 0) {
      throw std::invalid_argument(
          "PolynomialSparse::ModBinomial: binomial degree n must be >= 1");
    }
    if (degrees_.empty() || degrees_.back() < n) {
      return *this;
    }

    const element_type c_norm = ring_->Add(ring_->AdditiveIdentity(), c);
    const element_type one = ring_->MultiplicativeIdentity();
    const element_type neg_one = ring_->Neg(one);
    const bool c_is_zero = ring_->IsZero(c_norm);
    const bool c_is_one = ring_->IsOne(c_norm);
    const bool c_is_neg_one = alg::Equals(*ring_, c_norm, neg_one);

    std::vector<std::pair<uint64_t, element_type>> reduced;
    reduced.reserve(degrees_.size());

    uint64_t prev_q = 0;
    element_type c_pow = one;

    for (std::size_t i = 0; i < degrees_.size(); ++i) {
      const uint64_t d = degrees_[i];
      const uint64_t q = d / n;
      const uint64_t r = d % n;
      if (q == 0) {
        reduced.emplace_back(r, coeffs_[i]);
        continue;
      }
      if (c_is_zero) {
        continue;
      }
      if (c_is_one) {
        reduced.emplace_back(r, coeffs_[i]);
      } else if (c_is_neg_one) {
        reduced.emplace_back(
            r, ((q & 1ULL) != 0ULL) ? ring_->Neg(coeffs_[i]) : coeffs_[i]);
      } else {
        if (q > prev_q) {
          c_pow = ring_->Mul(c_pow, RingPow(c_norm, q - prev_q));
          prev_q = q;
        }
        element_type term = ring_->Mul(coeffs_[i], c_pow);
        if (!ring_->IsZero(term)) {
          reduced.emplace_back(r, std::move(term));
        }
      }
    }

    return FromTerms(ring_, std::move(reduced), var_);
  }

  // Evaluates the sparse polynomial at x using gap-powered Horner's method.
  element_type Eval(const element_type& x) const {
    RequireRing();
    if (degrees_.empty()) {
      return ring_->AdditiveIdentity();
    }
    element_type acc = coeffs_.back();
    for (std::size_t i = degrees_.size() - 1; i > 0; --i) {
      const uint64_t gap = degrees_[i] - degrees_[i - 1];
      acc = ring_->Add(ring_->Mul(acc, RingPow(x, gap)), coeffs_[i - 1]);
    }
    if (degrees_[0] > 0) {
      acc = ring_->Mul(acc, RingPow(x, degrees_[0]));
    }
    return acc;
  }

  PolynomialSparse operator+(const PolynomialSparse& other) const {
    return Add(other);
  }

  PolynomialSparse operator-(const PolynomialSparse& other) const {
    return Sub(other);
  }

  PolynomialSparse operator-() const { return Neg(); }

  PolynomialSparse operator*(const PolynomialSparse& other) const {
    return Mul(other);
  }

  PolynomialSparse& operator+=(const PolynomialSparse& other) {
    *this = Add(other);
    return *this;
  }

  PolynomialSparse& operator-=(const PolynomialSparse& other) {
    *this = Sub(other);
    return *this;
  }

  PolynomialSparse& operator*=(const PolynomialSparse& other) {
    *this = Mul(other);
    return *this;
  }

  bool operator==(const PolynomialSparse& other) const {
    if (ring_ != other.ring_) {
      if (!ring_ || !other.ring_) {
        return false;
      }
      if constexpr (requires(const R& lhs, const R& rhs) {
                      { lhs == rhs } -> std::convertible_to<bool>;
                    }) {
        if (!(*ring_ == *other.ring_)) {
          return false;
        }
      } else {
        return false;
      }
    }
    if (degrees_ != other.degrees_) {
      return false;
    }
    if (!ring_) {
      return true;
    }
    for (std::size_t i = 0; i < coeffs_.size(); ++i) {
      if (!alg::Equals(*ring_, coeffs_[i], other.coeffs_[i])) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const PolynomialSparse& other) const {
    return !(*this == other);
  }

  std::string ToString() const {
    if (!ring_ || degrees_.empty()) {
      return "0";
    }

    auto needs_brackets = [](const std::string& s) -> bool {
      if (s.empty()) {
        return false;
      }
      for (char ch : s) {
        if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
          return true;
        }
      }
      return false;
    };

    std::ostringstream oss;
    bool first_term = true;
    for (std::size_t idx = degrees_.size(); idx > 0; --idx) {
      const std::size_t i = idx - 1;
      const uint64_t deg = degrees_[i];
      const element_type& c = coeffs_[i];
      if (!first_term) {
        oss << " + ";
      }
      if (deg == 0) {
        oss << ring_->ToString(c);
      } else {
        if (!ring_->IsOne(c)) {
          const std::string coeff_str = ring_->ToString(c);
          if (needs_brackets(coeff_str)) {
            oss << "(" << coeff_str << ")";
          } else {
            oss << coeff_str;
          }
        }
        oss << var_;
        if (deg > 1) {
          oss << "^" << deg;
        }
      }
      first_term = false;
    }
    return oss.str();
  }

 private:
  PolynomialSparse(InternalTag, std::shared_ptr<R> ring,
                   std::vector<uint64_t> degrees,
                   std::vector<element_type> coeffs, std::string var)
      : ring_(std::move(ring)),
        degrees_(std::move(degrees)),
        coeffs_(std::move(coeffs)),
        var_(std::move(var)) {}

  static void ValidateRingAndVar(const std::shared_ptr<R>& ring,
                                 const std::string& var) {
    if (!ring) {
      throw std::invalid_argument(
          "PolynomialSparse: ring pointer cannot be null");
    }
    if (var.empty()) {
      throw std::invalid_argument(
          "PolynomialSparse: variable string cannot be empty");
    }
  }

  void RequireRing() const {
    if (!ring_) {
      throw std::runtime_error(
          "PolynomialSparse: polynomial has no associated ring");
    }
  }

  void ValidateSameRing(const PolynomialSparse& other) const {
    if (!ring_ || !other.ring_) {
      throw std::invalid_argument(
          "PolynomialSparse: operations with null ring are not allowed");
    }
    if (ring_ != other.ring_) {
      if constexpr (requires(const R& lhs, const R& rhs) {
                      { lhs == rhs } -> std::convertible_to<bool>;
                    }) {
        if (!(*ring_ == *other.ring_)) {
          throw std::invalid_argument(
              "PolynomialSparse: polynomials belong to different rings");
        }
      } else {
        throw std::invalid_argument(
            "PolynomialSparse: polynomials belong to different rings");
      }
    }
  }

  element_type RingPow(element_type base, uint64_t exp) const {
    element_type acc = ring_->MultiplicativeIdentity();
    while (exp > 0) {
      if ((exp & 1ULL) != 0ULL) {
        acc = ring_->Mul(acc, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = ring_->Mul(base, base);
      }
    }
    return acc;
  }

  std::shared_ptr<R> ring_{nullptr};
  std::vector<uint64_t> degrees_;
  std::vector<element_type> coeffs_;
  std::string var_ = "x";
};

template <Ring R>
inline std::pair<PolynomialSparse<R>, PolynomialSparse<R>> MonicDivRem(
    const PolynomialSparse<R>& a, const PolynomialSparse<R>& monic_b) {
  return PolynomialSparse<R>::MonicDivRem(a, monic_b);
}

template <Ring R>
inline std::ostream& operator<<(std::ostream& os,
                                const PolynomialSparse<R>& poly) {
  os << poly.ToString();
  return os;
}

}  // namespace poly

using poly::MonicDivRem;
using poly::PolynomialSparse;

}  // namespace xg

#endif  // XGALOIS_POLY_POLY_SPARSE_HPP_
