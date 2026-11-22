#ifndef XGALOIS_FIELD_GF_DUAL_BASIS_HPP_
#define XGALOIS_FIELD_GF_DUAL_BASIS_HPP_

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include "xgalois/field/gf_binary.hpp"

namespace xg {

// Dual-basis (and weakly-dual-basis) arithmetic wrapper for any binary
// extension field GF(2^m) (e.g., GaloisFieldBinaryExtension<Word>, 1 <= m <=
// 64).
//
// Given the polynomial (primal) basis {1, alpha, ..., alpha^{m-1}} and a
// nonzero shift element gamma in GF(2^m) (default gamma = 1; e.g., CCSDS
// Reed-Solomon uses a weakly-dual shift gamma != 1):
//   - Gram matrix G_{ij} = Tr(gamma * alpha^{i+j}) over GF(2)
//   - Dual basis {beta_0, ..., beta_{m-1}} satisfies Tr(gamma * alpha^i *
//     beta_j) = delta_{ij}, with beta_j = sum_k (G^{-1})_{kj} alpha^k.
//   - ToDual(a) computes d_i = Tr(gamma * alpha^i * a) = (G * a)_i.
//   - FromDual(d) recovers a = G^{-1} * d.
//   - MulPrimalDual(a, d_b) computes the dual coordinates of a * b via
//     Berlekamp's trace inner product + LFSR feedback shifts.
template <typename BinaryExtField>
class GaloisFieldDualBasis {
 public:
  using field_type = BinaryExtField;
  using element_type = typename BinaryExtField::element_type;

  explicit GaloisFieldDualBasis(const BinaryExtField& field)
      : GaloisFieldDualBasis(field, field.MultiplicativeIdentity()) {}

  GaloisFieldDualBasis(const BinaryExtField& field, element_type gamma)
      : field_(&field), gamma_(gamma), m_(field.Degree()) {
    if (field.Characteristic() != 2 || m_ < 1 || m_ > 64) {
      throw std::invalid_argument(
          "GaloisFieldDualBasis requires GF(2^m) with 1 <= m <= 64");
    }
    if (field.IsZero(gamma_)) {
      throw std::invalid_argument(
          "GaloisFieldDualBasis: shift element gamma must be nonzero");
    }

    mask_ = (m_ < 64) ? ((uint64_t{1} << m_) - 1ULL) : ~uint64_t{0};
    modulus_low_ = field.ModulusLow();

    // Precompute primal basis elements alpha_i = x^i in field representation.
    for (uint32_t i = 0; i < m_; ++i) {
      primal_basis_[i] =
          field_->SetElementValue(static_cast<element_type>(uint64_t{1} << i));
    }

    // Precompute O(1) absolute trace mask: bit i is Tr(alpha^i).
    trace_mask_ = 0ULL;
    for (uint32_t i = 0; i < m_; ++i) {
      if (ComputeTraceDirect(primal_basis_[i]) != 0U) {
        trace_mask_ |= (uint64_t{1} << i);
      }
    }

    // Precompute symmetric m x m binary Gram matrix:
    //   G_{i, j} = Tr(gamma * alpha^i * alpha^j).
    gram_rows_.fill(0ULL);
    for (uint32_t i = 0; i < m_; ++i) {
      const element_type gamma_ai = field_->Mul(gamma_, primal_basis_[i]);
      uint64_t row = 0ULL;
      for (uint32_t j = 0; j < m_; ++j) {
        const element_type prod = field_->Mul(gamma_ai, primal_basis_[j]);
        if (Trace(prod) != 0U) {
          row |= (uint64_t{1} << j);
        }
      }
      gram_rows_[i] = row;
    }

    // Invert G over GF(2) via Gauss-Jordan elimination.
    InvertGramMatrix();

    // Store dual basis elements beta_j in field representation. Since G is
    // symmetric, G^{-1} is also symmetric, so column j equals row j.
    for (uint32_t j = 0; j < m_; ++j) {
      dual_basis_[j] =
          field_->SetElementValue(static_cast<element_type>(gram_inv_rows_[j]));
    }
  }

  uint32_t Degree() const noexcept { return m_; }
  const BinaryExtField& Field() const noexcept { return *field_; }
  element_type Gamma() const noexcept { return gamma_; }

  element_type DualBasisElement(std::size_t index) const {
    if (index >= m_) {
      throw std::out_of_range(
          "GaloisFieldDualBasis::DualBasisElement: index out of range");
    }
    return dual_basis_[index];
  }

  // O(1) absolute trace Tr_{GF(2^m)/GF(2)}(a).
  uint8_t Trace(const element_type& a) const noexcept {
    const uint64_t bits =
        static_cast<uint64_t>(field_->GetElementValue(a)) & mask_;
    return static_cast<uint8_t>(std::popcount(bits & trace_mask_) & 1);
  }

  // Converts an element from polynomial (primal) basis to dual-basis
  // coordinates: d_i = Tr(gamma * alpha^i * a) = (G * a)_i.
  element_type ToDual(const element_type& poly_elem) const noexcept {
    const uint64_t a_bits =
        static_cast<uint64_t>(field_->GetElementValue(poly_elem)) & mask_;
    uint64_t dual_bits = 0ULL;
    for (uint32_t i = 0; i < m_; ++i) {
      if ((std::popcount(gram_rows_[i] & a_bits) & 1) != 0) {
        dual_bits |= (uint64_t{1} << i);
      }
    }
    return static_cast<element_type>(dual_bits);
  }

  // Converts an element from dual-basis coordinates back to polynomial
  // (primal) basis: a = G^{-1} * d.
  element_type FromDual(const element_type& dual_elem) const {
    const uint64_t d_bits = static_cast<uint64_t>(dual_elem) & mask_;
    uint64_t poly_bits = 0ULL;
    for (uint32_t i = 0; i < m_; ++i) {
      if ((std::popcount(gram_inv_rows_[i] & d_bits) & 1) != 0) {
        poly_bits |= (uint64_t{1} << i);
      }
    }
    return field_->SetElementValue(static_cast<element_type>(poly_bits));
  }

  // Berlekamp primal-dual multiplier: given poly_a in polynomial basis and
  // dual_b in dual basis, computes the dual-basis representation of (a * b)
  // using trace inner products and LFSR shifts of the trace register.
  element_type MulPrimalDual(const element_type& poly_a,
                             const element_type& dual_b) const noexcept {
    const uint64_t a_bits =
        static_cast<uint64_t>(field_->GetElementValue(poly_a)) & mask_;
    uint64_t reg = static_cast<uint64_t>(dual_b) & mask_;
    uint64_t out = 0ULL;

    for (uint32_t k = 0; k < m_; ++k) {
      const uint64_t prod_bit =
          static_cast<uint64_t>(std::popcount(a_bits & reg) & 1);
      out |= (prod_bit << k);
      if (k + 1 < m_) {
        // Advance LFSR trace register: s_{k + m} = sum_{r=0}^{m-1} p_r *
        // s_{k+r}.
        const uint64_t next_bit =
            static_cast<uint64_t>(std::popcount(reg & modulus_low_) & 1);
        reg = (reg >> 1) | (next_bit << (m_ - 1));
      }
    }
    return static_cast<element_type>(out);
  }

  // Verifies the defining dual orthogonality relation:
  //   Tr(gamma * alpha^i * beta_j) == (i == j ? 1 : 0) for all 0 <= i, j < m.
  bool VerifyDualOrthogonality() const {
    for (uint32_t i = 0; i < m_; ++i) {
      const element_type gamma_ai = field_->Mul(gamma_, primal_basis_[i]);
      for (uint32_t j = 0; j < m_; ++j) {
        const uint8_t tr = Trace(field_->Mul(gamma_ai, dual_basis_[j]));
        const uint8_t expected = (i == j) ? uint8_t{1} : uint8_t{0};
        if (tr != expected) {
          return false;
        }
      }
    }
    return true;
  }

 private:
  uint8_t ComputeTraceDirect(const element_type& a) const {
    element_type acc = a;
    element_type cur = a;
    for (uint32_t i = 1; i < m_; ++i) {
      cur = field_->Mul(cur, cur);
      acc = field_->Add(acc, cur);
    }
    return field_->IsZero(acc) ? uint8_t{0} : uint8_t{1};
  }

  void InvertGramMatrix() {
    std::array<uint64_t, 64> work = gram_rows_;
    gram_inv_rows_.fill(0ULL);
    for (uint32_t i = 0; i < m_; ++i) {
      gram_inv_rows_[i] = uint64_t{1} << i;
    }

    for (uint32_t col = 0; col < m_; ++col) {
      uint32_t pivot = m_;
      for (uint32_t row = col; row < m_; ++row) {
        if (((work[row] >> col) & 1ULL) != 0ULL) {
          pivot = row;
          break;
        }
      }
      if (pivot == m_) {
        throw std::runtime_error(
            "GaloisFieldDualBasis: trace Gram matrix is singular");
      }
      if (pivot != col) {
        std::swap(work[pivot], work[col]);
        std::swap(gram_inv_rows_[pivot], gram_inv_rows_[col]);
      }
      for (uint32_t row = 0; row < m_; ++row) {
        if (row != col && (((work[row] >> col) & 1ULL) != 0ULL)) {
          work[row] ^= work[col];
          gram_inv_rows_[row] ^= gram_inv_rows_[col];
        }
      }
    }
  }

  const BinaryExtField* field_;
  element_type gamma_;
  uint32_t m_ = 0;
  uint64_t mask_ = 0;
  uint64_t modulus_low_ = 0;
  uint64_t trace_mask_ = 0;
  std::array<element_type, 64> primal_basis_{};
  std::array<element_type, 64> dual_basis_{};
  std::array<uint64_t, 64> gram_rows_{};
  std::array<uint64_t, 64> gram_inv_rows_{};
};

template <typename BinaryExtField>
GaloisFieldDualBasis(const BinaryExtField&)
    -> GaloisFieldDualBasis<BinaryExtField>;

template <typename BinaryExtField>
GaloisFieldDualBasis(const BinaryExtField&,
                     typename BinaryExtField::element_type)
    -> GaloisFieldDualBasis<BinaryExtField>;

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_DUAL_BASIS_HPP_
