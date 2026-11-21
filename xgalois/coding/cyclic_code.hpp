#ifndef XGALOIS_CODING_CYCLIC_CODE_HPP_
#define XGALOIS_CODING_CYCLIC_CODE_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/cyclic_syndrome_decoder.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"
#include "xgalois/coding/encoder/cyclic_encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

namespace detail {

inline bool IsBchDefiningSet(const std::vector<bool> &in_t, std::size_t n,
                             uint64_t q) {
  std::size_t t_size = 0;
  for (bool b : in_t) {
    if (b) {
      ++t_size;
    }
  }
  if (t_size == 0 || t_size >= n) {
    return false;
  }

  std::vector<bool> in_tu(n, false);
  std::vector<bool> covered(n, false);

  // Test all units u in Z_n^* (corresponding to all primitive n-th roots of
  // unity beta^(u^{-1})).
  for (std::size_t u = 1; u < n; ++u) {
    if (std::gcd(u, n) != 1) {
      continue;
    }
    std::fill(in_tu.begin(), in_tu.end(), false);
    for (std::size_t i = 0; i < n; ++i) {
      if (in_t[i]) {
        in_tu[(i * u) % n] = true;
      }
    }

    // Check every maximal circular consecutive run in T_u.
    for (std::size_t start = 0; start < n; ++start) {
      if (!in_tu[start] || in_tu[(start + n - 1) % n]) {
        continue;
      }
      std::fill(covered.begin(), covered.end(), false);
      std::size_t covered_count = 0;
      std::size_t pos = start;
      while (in_tu[pos]) {
        std::size_t cur = pos;
        while (!covered[cur]) {
          covered[cur] = true;
          ++covered_count;
          cur = static_cast<std::size_t>(
              (static_cast<uint64_t>(cur) * (q % n)) % n);
        }
        pos = (pos + 1) % n;
        if (pos == start) {
          break;
        }
      }
      if (covered_count == t_size) {
        return true;
      }
    }
  }
  return false;
}

template <typename GaloisField>
inline PolynomialDense<GaloisField> MakeXnMinus1(
    const std::shared_ptr<GaloisField> &field, std::size_t n) {
  using RawType = typename GaloisField::element_type;
  std::vector<RawType> coeffs(n + 1, field->AdditiveIdentity());
  const RawType one = field->MultiplicativeIdentity();
  coeffs[0] = field->Neg(one);
  coeffs[n] = one;
  return PolynomialDense<GaloisField>::FromRaw(std::move(coeffs), field);
}

template <typename GaloisField>
inline std::size_t ValidateCyclicDimension(
    const std::shared_ptr<GaloisField> &field, std::size_t length,
    const PolynomialDense<GaloisField> &generator_poly) {
  if (!field) {
    throw std::invalid_argument("CyclicCode: field cannot be null");
  }
  if (length < 2) {
    throw std::invalid_argument("CyclicCode: code length must be at least 2");
  }
  if (!generator_poly.HasField() || generator_poly.Field() != field) {
    throw std::invalid_argument(
        "CyclicCode: generator polynomial field mismatch");
  }
  if (generator_poly.IsZero()) {
    throw std::invalid_argument(
        "CyclicCode: generator polynomial cannot be zero");
  }
  const int deg = generator_poly.Degree();
  if (deg <= 0 || static_cast<std::size_t>(deg) >= length) {
    throw std::invalid_argument(
        "CyclicCode: generator polynomial must be non-constant with degree "
        "strictly less than code length");
  }
  return length - static_cast<std::size_t>(deg);
}

}  // namespace detail

template <typename GaloisField>
class CyclicCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = typename Base::matrix_type;
  using polynomial_type = PolynomialDense<GaloisField>;

  CyclicCode(std::shared_ptr<GaloisField> field, std::size_t length,
             const polynomial_type &generator_poly,
             Metric metric = Metric::HAMMING)
      : Base(field, length,
             detail::ValidateCyclicDimension(field, length, generator_poly),
             metric, "Cyclic", "CyclicSyndrome"),
        generator_poly_(generator_poly.Monic()) {
    const polynomial_type xn_minus_1 =
        detail::MakeXnMinus1(this->field_, this->length_);
    auto [h_poly, rem] = xn_minus_1.DivRem(generator_poly_);
    if (!rem.IsZero()) {
      throw std::invalid_argument(
          "CyclicCode: generator polynomial must divide x^n - 1");
    }
    parity_check_poly_ = std::move(h_poly);
    reciprocal_parity_check_poly_ =
        parity_check_poly_.Reverse(this->dimension_).Monic();

    BuildMatrices();
    RegisterEncodersAndDecoders();
  }

  static CyclicCode FromParityCheckPolynomial(
      std::shared_ptr<GaloisField> field, std::size_t length,
      const polynomial_type &parity_check_poly,
      Metric metric = Metric::HAMMING) {
    if (!field) {
      throw std::invalid_argument(
          "CyclicCode::FromParityCheckPolynomial: field cannot be null");
    }
    if (length < 2) {
      throw std::invalid_argument(
          "CyclicCode::FromParityCheckPolynomial: length must be at least 2");
    }
    if (!parity_check_poly.HasField() || parity_check_poly.Field() != field) {
      throw std::invalid_argument(
          "CyclicCode::FromParityCheckPolynomial: field mismatch");
    }
    if (parity_check_poly.IsZero() || parity_check_poly.Degree() <= 0 ||
        static_cast<std::size_t>(parity_check_poly.Degree()) >= length) {
      throw std::invalid_argument(
          "CyclicCode::FromParityCheckPolynomial: parity-check polynomial must "
          "be non-constant with degree < length");
    }
    const polynomial_type xn_minus_1 = detail::MakeXnMinus1(field, length);
    auto [g_poly, rem] = xn_minus_1.DivRem(parity_check_poly.Monic());
    if (!rem.IsZero()) {
      throw std::invalid_argument(
          "CyclicCode::FromParityCheckPolynomial: parity-check polynomial must "
          "divide x^n - 1");
    }
    return CyclicCode(std::move(field), length, g_poly, metric);
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  matrix_type SystematicGeneratorMatrix() const override {
    return linalg::SystematicForm(generator_matrix_).g_sys;
  }

  const polynomial_type &GeneratorPolynomial() const noexcept {
    return generator_poly_;
  }

  const polynomial_type &ParityCheckPolynomial() const noexcept {
    return parity_check_poly_;
  }

  const polynomial_type &ReciprocalParityCheckPolynomial() const noexcept {
    return reciprocal_parity_check_poly_;
  }

  // Dual cyclic code generated by the monic reciprocal polynomial
  // h*(x) = (x^k * h(x^{-1}))_monic, whose generator matrix equals H.
  CyclicCode DualCode() const {
    return CyclicCode(this->field_, this->length_,
                      reciprocal_parity_check_poly_, this->metric_);
  }

  virtual bool IsBCH() const {
    const std::size_t n = this->length_;
    const uint64_t q = this->field_->Order();
    if (q == 0 || std::gcd(static_cast<uint64_t>(n), q) != 1) {
      return false;
    }

    // Compute multiplicative order m = ord_n(q).
    std::size_t m = 1;
    uint64_t cur = q % n;
    while (cur != 1 && m <= n) {
      cur = (cur * (q % n)) % n;
      ++m;
    }
    if (cur != 1) {
      return false;
    }

    // Case 1: n | (q - 1) => splitting field is the base field itself.
    if (m == 1) {
      const storage_type gen = this->field_->MultiplicativeGenerator();
      const storage_type beta = this->field_->Pow(gen, (q - 1) / n);
      std::vector<bool> in_t(n, false);
      storage_type root = this->field_->MultiplicativeIdentity();
      for (std::size_t i = 0; i < n; ++i) {
        if (this->field_->IsZero(generator_poly_.EvalRaw(root))) {
          in_t[i] = true;
        }
        root = this->field_->Mul(root, beta);
      }
      return detail::IsBchDefiningSet(in_t, n, q);
    }

    // Case 2: Binary cyclic code with splitting field GF(2^m), m <= 16.
    if constexpr (std::is_same_v<GaloisField, GaloisFieldBinary>) {
      if (m <= 16) {
        auto ext = std::make_shared<GF2X<uint16_t>>(static_cast<uint8_t>(m));
        const uint16_t gen = ext->MultiplicativeGenerator();
        const uint16_t beta = ext->Pow(gen, (ext->Order() - 1) / n);
        std::vector<bool> in_t(n, false);
        uint16_t root = ext->MultiplicativeIdentity();
        const int deg = generator_poly_.Degree();
        for (std::size_t i = 0; i < n; ++i) {
          // Evaluate g(root) via Horner's method in GF(2^m).
          uint16_t acc = ext->AdditiveIdentity();
          for (int j = deg; j >= 0; --j) {
            acc = ext->Mul(acc, root);
            if (!this->field_->IsZero(
                    generator_poly_.RawCoeff(static_cast<std::size_t>(j)))) {
              acc = ext->Add(acc, ext->MultiplicativeIdentity());
            }
          }
          if (ext->IsZero(acc)) {
            in_t[i] = true;
          }
          root = ext->Mul(root, beta);
        }
        return detail::IsBchDefiningSet(in_t, n, q);
      }
      return false;
    }

    // Case 3: Prime-field cyclic code with small splitting field GF(p^m).
    if constexpr (std::is_integral_v<storage_type>) {
      if (this->field_->Degree() == 1) {
        uint64_t ext_order = 0;
        try {
          ext_order = utils::SafeIntegerPower(q, m);
        } catch (const std::overflow_error &) {
          return false;
        }
        if (ext_order <= 65536) {
          using ExtField = GaloisFieldExtension<uint32_t>;
          using ExtElem = typename ExtField::element_type;
          auto ext = std::make_shared<ExtField>(
              std::make_pair(q, static_cast<uint64_t>(m)));
          const ExtElem gen = ext->MultiplicativeGenerator();
          const ExtElem beta = ext->Pow(gen, (ext_order - 1) / n);
          std::vector<bool> in_t(n, false);
          ExtElem root = ext->MultiplicativeIdentity();
          const int deg = generator_poly_.Degree();
          for (std::size_t i = 0; i < n; ++i) {
            ExtElem acc = ext->AdditiveIdentity();
            for (int j = deg; j >= 0; --j) {
              acc = ext->Mul(acc, root);
              const storage_type c =
                  generator_poly_.RawCoeff(static_cast<std::size_t>(j));
              if (!this->field_->IsZero(c)) {
                acc = ext->Add(acc, ext->FromInt(static_cast<uint64_t>(
                                        this->field_->GetElementValue(c))));
              }
            }
            if (ext->IsZero(acc)) {
              in_t[i] = true;
            }
            root = ext->Mul(root, beta);
          }
          return detail::IsBchDefiningSet(in_t, n, q);
        }
      }
    }

    return false;
  }

  // Right cyclic shift by `shift` positions (supports arbitrary negative or
  // large shifts): out[(i + shift) mod n] = word[i].
  codeword_type CyclicShift(const codeword_type &word, int shift = 1) const {
    if (!word.HasField() || word.Field() != this->field_) {
      throw std::invalid_argument("CyclicShift: field mismatch");
    }
    if (word.Dimension() != 1 || word.Size() != this->length_) {
      throw std::invalid_argument(
          "CyclicShift: word must be 1D with size equal to code length");
    }

    const int64_t n = static_cast<int64_t>(this->length_);
    const int64_t norm_shift = ((static_cast<int64_t>(shift) % n) + n) % n;
    const std::size_t s = static_cast<std::size_t>(norm_shift);

    using RawArray = typename codeword_type::raw_array_type;
    using ShapeType = typename codeword_type::shape_type;
    RawArray out = RawArray::from_shape(ShapeType{this->length_});
    const storage_type *src = word.RawData().data();
    storage_type *dst = out.data();
    for (std::size_t i = 0; i < this->length_; ++i) {
      dst[(i + s) % this->length_] = src[i];
    }
    return codeword_type::FromRaw(this->field_, std::move(out));
  }

  polynomial_type VectorToPolynomial(const codeword_type &vec) const {
    if (!vec.HasField() || vec.Field() != this->field_) {
      throw std::invalid_argument("VectorToPolynomial: field mismatch");
    }
    if (vec.Dimension() != 1) {
      throw std::invalid_argument("VectorToPolynomial: input must be 1D");
    }
    std::vector<storage_type> raw(vec.RawSpan().begin(), vec.RawSpan().end());
    return polynomial_type::FromRaw(std::move(raw), this->field_);
  }

  codeword_type PolynomialToVector(const polynomial_type &poly) const {
    if (!poly.HasField() || poly.Field() != this->field_) {
      throw std::invalid_argument("PolynomialToVector: field mismatch");
    }
    codeword_type result = codeword_type::Zeros(this->field_, {this->length_});
    if (poly.IsZero()) {
      return result;
    }
    const std::size_t limit = std::min<std::size_t>(
        this->length_, static_cast<std::size_t>(poly.Degree()) + 1);
    storage_type *dst = result.RawData().data();
    for (std::size_t i = 0; i < limit; ++i) {
      dst[i] = poly.RawCoeff(i);
    }
    return result;
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + "] cyclic code over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 private:
  void BuildMatrices() {
    BuildGeneratorMatrix();
    BuildParityCheckMatrix();
  }

  void BuildGeneratorMatrix() {
    const std::size_t k = this->dimension_;
    const std::size_t n = this->length_;
    const std::size_t deg_g =
        static_cast<std::size_t>(generator_poly_.Degree());

    generator_matrix_ = matrix_type::Zeros(this->field_, {k, n});
    storage_type *g_data = generator_matrix_.RawData().data();
    for (std::size_t i = 0; i < k; ++i) {
      storage_type *row = g_data + i * n;
      for (std::size_t j = 0; j <= deg_g; ++j) {
        row[i + j] = generator_poly_.RawCoeff(j);
      }
    }
  }

  void BuildParityCheckMatrix() {
    // Rows of H are cyclic shifts x^i * h*(x) for 0 <= i < n - k, where
    // h*(x) = (x^k * h(x^{-1}))_monic. Because g(x)*h(x) = x^n - 1 has zero
    // coefficients for x^1 .. x^(n-1), every row of G is strictly orthogonal
    // to every row of H: G * H^T = 0.
    const std::size_t k = this->dimension_;
    const std::size_t n = this->length_;
    const std::size_t r = n - k;

    parity_check_matrix_ = matrix_type::Zeros(this->field_, {r, n});
    storage_type *h_data = parity_check_matrix_.RawData().data();
    for (std::size_t i = 0; i < r; ++i) {
      storage_type *row = h_data + i * n;
      for (std::size_t j = 0; j <= k; ++j) {
        row[i + j] = reciprocal_parity_check_poly_.RawCoeff(j);
      }
    }
  }

  void RegisterEncodersAndDecoders() {
    this->RegisterEncoder("Cyclic", [](const AbstractCode<GaloisField> &code) {
      const auto &cyc = dynamic_cast<const CyclicCode<GaloisField> &>(code);
      return std::make_shared<CyclicEncoder<GaloisField>>(cyc, false);
    });
    this->RegisterEncoder(
        "CyclicSystematic", [](const AbstractCode<GaloisField> &code) {
          const auto &cyc = dynamic_cast<const CyclicCode<GaloisField> &>(code);
          return std::make_shared<CyclicEncoder<GaloisField>>(cyc, true);
        });
    this->RegisterEncoder(
        "GeneratorMatrix", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(lin);
        });

    this->RegisterDecoder(
        "CyclicSyndrome", [](const AbstractCode<GaloisField> &code) {
          const auto &cyc = dynamic_cast<const CyclicCode<GaloisField> &>(code);
          return std::make_shared<CyclicSyndromeDecoder<GaloisField>>(cyc);
        });
    this->RegisterDecoder(
        "Syndrome", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<SyndromeDecoder<GaloisField>>(lin);
        });
  }

  polynomial_type generator_poly_;
  polynomial_type parity_check_poly_;
  polynomial_type reciprocal_parity_check_poly_;
  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_CYCLIC_CODE_HPP_
