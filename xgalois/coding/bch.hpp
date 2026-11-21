#ifndef XGALOIS_CODING_BCH_HPP_
#define XGALOIS_CODING_BCH_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/coding/cyclic_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

namespace detail {

template <typename BaseField, typename ExtField>
inline typename ExtField::element_type EmbedBaseToExt(
    const BaseField &base_field, const ExtField &ext_field,
    const typename BaseField::element_type &base_raw) {
  using BaseRaw = typename BaseField::element_type;
  using ExtRaw = typename ExtField::element_type;
  if constexpr (std::is_same_v<BaseField, ExtField>) {
    return base_raw;
  } else if constexpr (std::is_integral_v<ExtRaw>) {
    if (base_field.IsZero(base_raw)) {
      return ext_field.AdditiveIdentity();
    }
    if (base_field.IsOne(base_raw)) {
      return ext_field.MultiplicativeIdentity();
    }
    throw std::invalid_argument(
        "BCH: cannot embed non-binary base element into binary extension "
        "field");
  } else {
    if (base_field.IsZero(base_raw)) {
      return ext_field.AdditiveIdentity();
    }
    if constexpr (std::is_integral_v<BaseRaw>) {
      const uint64_t val =
          static_cast<uint64_t>(base_field.GetElementValue(base_raw));
      return ext_field.FromInt(val);
    } else {
      throw std::invalid_argument(
          "BCH: unsupported base-to-extension field embedding");
    }
  }
}

template <typename BaseField, typename ExtField>
inline std::optional<typename BaseField::element_type> TryProjectExtToBase(
    const BaseField &base_field, const ExtField &ext_field,
    const typename ExtField::element_type &ext_raw) {
  using BaseRaw = typename BaseField::element_type;
  using ExtRaw = typename ExtField::element_type;
  if constexpr (std::is_same_v<BaseField, ExtField>) {
    return ext_raw;
  } else if constexpr (std::is_integral_v<ExtRaw>) {
    if (ext_field.IsZero(ext_raw)) {
      return base_field.AdditiveIdentity();
    }
    if (ext_field.IsOne(ext_raw)) {
      return base_field.MultiplicativeIdentity();
    }
    return std::nullopt;
  } else {
    const ExtRaw canon = ext_field.GetElementValue(ext_raw);
    if (canon.IsZero()) {
      return base_field.AdditiveIdentity();
    }
    if (canon.Degree() > 0) {
      return std::nullopt;
    }
    if constexpr (std::is_integral_v<BaseRaw>) {
      return base_field.SetElementValue(static_cast<BaseRaw>(canon[0]));
    } else {
      return std::nullopt;
    }
  }
}

template <typename BaseField, typename ExtField>
inline std::size_t ResolveBchLength(const std::shared_ptr<BaseField> &base,
                                    const std::shared_ptr<ExtField> &ext,
                                    std::size_t length) {
  if (!base || !ext) {
    throw std::invalid_argument(
        "BCHCode: base and extension fields cannot be null");
  }
  if (base->Characteristic() != ext->Characteristic()) {
    throw std::invalid_argument(
        "BCHCode: base and extension fields must have the same characteristic");
  }
  if (ext->Degree() % base->Degree() != 0) {
    throw std::invalid_argument(
        "BCHCode: extension field degree must be a multiple of base field "
        "degree");
  }
  if constexpr (std::is_same_v<BaseField, ExtField>) {
    if (base != ext) {
      throw std::invalid_argument(
          "BCHCode: identical field types must share the same field instance");
    }
  }
  const uint64_t ext_order = ext->Order();
  if (ext_order <= 2) {
    throw std::invalid_argument(
        "BCHCode: extension field order must be greater than 2");
  }
  const uint64_t group_order = ext_order - 1;
  const std::size_t n =
      (length == 0) ? static_cast<std::size_t>(group_order) : length;
  if (n < 2) {
    throw std::invalid_argument("BCHCode: code length must be at least 2");
  }
  if (group_order % n != 0) {
    throw std::invalid_argument(
        "BCHCode: code length must divide |ExtField| - 1");
  }
  return n;
}

template <typename BaseField, typename ExtField>
inline PolynomialDense<BaseField> BuildBchGeneratorPolynomial(
    const std::shared_ptr<BaseField> &base_field,
    const std::shared_ptr<ExtField> &ext_field, std::size_t n,
    std::size_t designed_distance, std::size_t b) {
  using BasePoly = PolynomialDense<BaseField>;
  using ExtPoly = PolynomialDense<ExtField>;
  using BaseRaw = typename BaseField::element_type;
  using ExtRaw = typename ExtField::element_type;

  if (designed_distance < 2 || designed_distance > n) {
    throw std::invalid_argument(
        "BCHCode: designed distance must be in [2, length]");
  }

  const uint64_t q = base_field->Order();
  const uint64_t group_order = ext_field->Order() - 1;
  const ExtRaw gamma = ext_field->MultiplicativeGenerator();
  const ExtRaw beta = ext_field->Pow(gamma, group_order / n);

  // Precompute powers beta^0, beta^1, ..., beta^(n-1).
  std::vector<ExtRaw> beta_pows(n);
  beta_pows[0] = ext_field->MultiplicativeIdentity();
  for (std::size_t i = 1; i < n; ++i) {
    beta_pows[i] = ext_field->Mul(beta_pows[i - 1], beta);
  }

  // Distinct q-cyclotomic cosets intersecting {b, b+1, ..., b+delta-2} mod n.
  // Minimal polynomials of distinct cyclotomic cosets are distinct monic
  // irreducibles over BaseField, so their LCM equals their product.
  std::vector<bool> covered(n, false);
  BasePoly g = BasePoly::One(base_field);
  const ExtRaw ext_one = ext_field->MultiplicativeIdentity();

  for (std::size_t offset = 0; offset < designed_distance - 1; ++offset) {
    const std::size_t start = (b + offset) % n;
    if (covered[start]) {
      continue;
    }

    // Build minimal polynomial m_{beta^start}(x) = prod_{j in C_start} (x -
    // beta^j) in ExtField[x].
    ExtPoly min_poly_ext = ExtPoly::One(ext_field);
    std::size_t cur = start;
    while (!covered[cur]) {
      covered[cur] = true;
      std::vector<ExtRaw> linear = {ext_field->Neg(beta_pows[cur]), ext_one};
      min_poly_ext *= ExtPoly::FromRaw(std::move(linear), ext_field);
      cur =
          static_cast<std::size_t>((static_cast<uint64_t>(cur) * (q % n)) % n);
    }

    // Project coefficients from ExtField to BaseField.
    const std::size_t deg = static_cast<std::size_t>(min_poly_ext.Degree());
    std::vector<BaseRaw> base_coeffs(deg + 1);
    for (std::size_t i = 0; i <= deg; ++i) {
      auto projected = TryProjectExtToBase(*base_field, *ext_field,
                                           min_poly_ext.RawCoeff(i));
      if (!projected.has_value()) {
        throw std::runtime_error(
            "BCHCode: minimal polynomial coefficient did not project into "
            "base field");
      }
      base_coeffs[i] = std::move(*projected);
    }
    g *= BasePoly::FromRaw(std::move(base_coeffs), base_field);
  }

  if (g.IsZero() || g.Degree() <= 0 ||
      static_cast<std::size_t>(g.Degree()) >= n) {
    throw std::invalid_argument(
        "BCHCode: designed distance too large; generator polynomial "
        "degenerates "
        "to x^n - 1");
  }
  return g.Monic();
}

}  // namespace detail

template <typename BaseField, typename ExtField>
class BCHCode;

// Algebraic BCH decoder using Berlekamp-Massey in ExtField, Chien search, and
// Forney's algorithm (generalized for arbitrary starting root b and both binary
// and non-binary BCH codes).
template <typename BaseField, typename ExtField>
class BCHDecoder : public Decoder<BaseField> {
 public:
  using Base = Decoder<BaseField>;
  using field_type = BaseField;
  using ext_field_type = ExtField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using ext_storage_type = typename ExtField::element_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using ext_poly_type = PolynomialDense<ExtField>;

  explicit BCHDecoder(const BCHCode<BaseField, ExtField> &code);

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto base_field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != base_field) {
      throw std::invalid_argument(
          "BCHDecoder::DecodeWithResult: field mismatch");
    }
    if (received_word.Dimension() != 1 || received_word.Size() != n_) {
      throw std::invalid_argument(
          "BCHDecoder::DecodeWithResult: received word must be 1D with length "
          "equal to code length");
    }

    const std::size_t num_syndromes = designed_distance_ - 1;
    const ext_storage_type ext_zero = ext_field_->AdditiveIdentity();
    const ext_storage_type ext_one = ext_field_->MultiplicativeIdentity();

    // Embed received word r(x) into ExtField[x].
    std::vector<ext_storage_type> r_ext_coeffs(n_);
    const storage_type *r_data = received_word.RawData().data();
    for (std::size_t j = 0; j < n_; ++j) {
      r_ext_coeffs[j] =
          detail::EmbedBaseToExt(*base_field, *ext_field_, r_data[j]);
    }
    ext_poly_type r_poly =
        ext_poly_type::FromRaw(std::move(r_ext_coeffs), ext_field_);

    // Compute syndromes S_k = r(beta^(b + k)) for k = 0 .. delta - 2.
    std::vector<ext_storage_type> syndromes(num_syndromes, ext_zero);
    bool all_zero = true;
    for (std::size_t k = 0; k < num_syndromes; ++k) {
      syndromes[k] = r_poly.EvalRaw(syndrome_roots_[k]);
      if (!ext_field_->IsZero(syndromes[k])) {
        all_zero = false;
      }
    }

    if (all_zero) {
      if (!this->code_->Contains(received_word)) {
        return result_type::Failure(
            "BCH decoding failed: word has zero BCH syndromes but is not a "
            "codeword");
      }
      auto msg = this->code_->Unencode(received_word);
      return result_type::Success(received_word, std::move(msg), 0);
    }

    if (t_ == 0) {
      return result_type::Failure(
          "BCH decoding failed: error detected with decoding radius 0");
    }

    // Berlekamp-Massey algorithm in ExtField to find the error-locator
    // polynomial Lambda(x) = prod_l (1 - X_l * x).
    std::vector<ext_storage_type> C = {ext_one};
    std::vector<ext_storage_type> B = {ext_one};
    std::size_t L = 0;
    std::size_t m = 1;
    ext_storage_type b_val = ext_one;

    for (std::size_t step = 0; step < num_syndromes; ++step) {
      ext_storage_type d = syndromes[step];
      const std::size_t max_i = std::min<std::size_t>(L, C.size() - 1);
      for (std::size_t i = 1; i <= max_i; ++i) {
        if (!ext_field_->IsZero(C[i]) &&
            !ext_field_->IsZero(syndromes[step - i])) {
          d = ext_field_->Add(d, ext_field_->Mul(C[i], syndromes[step - i]));
        }
      }

      if (ext_field_->IsZero(d)) {
        ++m;
      } else {
        std::vector<ext_storage_type> T = C;
        const ext_storage_type factor =
            ext_field_->Mul(d, ext_field_->Inv(b_val));
        if (C.size() < B.size() + m) {
          C.resize(B.size() + m, ext_zero);
        }
        for (std::size_t i = 0; i < B.size(); ++i) {
          if (!ext_field_->IsZero(B[i])) {
            C[i + m] = ext_field_->Sub(C[i + m], ext_field_->Mul(factor, B[i]));
          }
        }
        if (2 * L <= step) {
          L = step + 1 - L;
          B = std::move(T);
          b_val = d;
          m = 1;
        } else {
          ++m;
        }
      }
    }

    ext_poly_type lambda = ext_poly_type::FromRaw(std::move(C), ext_field_);
    const int nu_int = lambda.Degree();
    if (nu_int <= 0 || static_cast<std::size_t>(nu_int) > t_ ||
        static_cast<std::size_t>(nu_int) != L) {
      return result_type::Failure(
          "BCH decoding failed: error weight exceeds decoding radius");
    }
    const std::size_t nu = static_cast<std::size_t>(nu_int);

    // Chien search: find all positions j in [0, n - 1] with
    // Lambda(beta^(-j)) == 0.
    std::vector<std::size_t> error_positions;
    error_positions.reserve(nu);
    for (std::size_t j = 0; j < n_; ++j) {
      if (ext_field_->IsZero(lambda.EvalRaw(inv_beta_pows_[j]))) {
        error_positions.push_back(j);
      }
    }
    if (error_positions.size() != nu) {
      return result_type::Failure(
          "BCH decoding failed: error-locator polynomial does not split into "
          "distinct valid error locations");
    }

    codeword_type corrected = received_word;
    if (base_field->Order() == 2) {
      // Binary BCH: every nonzero error magnitude is 1.
      const storage_type base_one = base_field->MultiplicativeIdentity();
      for (std::size_t pos : error_positions) {
        corrected.Raw(pos) = base_field->Sub(corrected.Raw(pos), base_one);
      }
    } else {
      // Non-binary BCH: Forney's algorithm
      //   Omega(x) = (S(x) * Lambda(x)) mod x^(delta - 1)
      //   e_{j_l}  = - X_l^(1 - b) * Omega(X_l^{-1}) / Lambda'(X_l^{-1}).
      ext_poly_type s_poly = ext_poly_type::FromRaw(syndromes, ext_field_);
      ext_poly_type prod = s_poly * lambda;
      std::vector<ext_storage_type> omega_raw(num_syndromes, ext_zero);
      for (std::size_t i = 0; i < num_syndromes; ++i) {
        omega_raw[i] = prod.RawCoeff(i);
      }
      ext_poly_type omega =
          ext_poly_type::FromRaw(std::move(omega_raw), ext_field_);
      ext_poly_type lambda_deriv = lambda.Derivative();

      const int64_t forney_exp = 1 - static_cast<int64_t>(b_);
      for (std::size_t pos : error_positions) {
        const ext_storage_type &x_inv = inv_beta_pows_[pos];
        const ext_storage_type denom = lambda_deriv.EvalRaw(x_inv);
        if (ext_field_->IsZero(denom)) {
          return result_type::Failure(
              "BCH decoding failed: zero denominator in Forney algorithm");
        }
        const ext_storage_type numer = omega.EvalRaw(x_inv);
        const ext_storage_type x_pow =
            ext_field_->PowSigned(beta_pows_[pos], forney_exp);
        const ext_storage_type err_ext = ext_field_->Neg(
            ext_field_->Mul(x_pow, ext_field_->Div(numer, denom)));

        auto err_base =
            detail::TryProjectExtToBase(*base_field, *ext_field_, err_ext);
        if (!err_base.has_value() || base_field->IsZero(*err_base)) {
          return result_type::Failure(
              "BCH decoding failed: computed error magnitude is invalid for "
              "base field");
        }
        corrected.Raw(pos) = base_field->Sub(corrected.Raw(pos), *err_base);
      }
    }

    // Final codeword verification guarantees we never return a corrupted word.
    if (!this->code_->Contains(corrected)) {
      return result_type::Failure(
          "BCH decoding failed: corrected word is not a valid codeword "
          "(uncorrectable error)");
    }

    auto msg = this->code_->Unencode(corrected);
    return result_type::Success(std::move(corrected), std::move(msg), nu);
  }

  std::string Name() const override { return "BCH"; }

  std::size_t DecodingRadius() const override { return t_; }

  std::string ToString() const override {
    return "BCH (Berlekamp-Massey/Chien/Forney) decoder for [" +
           std::to_string(n_) + ", " +
           std::to_string(this->code_->Dimension()) +
           "] BCH code (t = " + std::to_string(t_) + ")";
  }

 private:
  std::shared_ptr<ExtField> ext_field_;
  std::size_t n_;
  std::size_t designed_distance_;
  std::size_t b_;
  std::size_t t_;
  std::vector<ext_storage_type> beta_pows_;
  std::vector<ext_storage_type> inv_beta_pows_;
  std::vector<ext_storage_type> syndrome_roots_;
};

template <typename BaseField = GaloisFieldBinary,
          typename ExtField = GF2X<uint16_t>>
class BCHCode : public CyclicCode<BaseField> {
 public:
  using Base = CyclicCode<BaseField>;
  using field_type = BaseField;
  using ext_field_type = ExtField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using ext_storage_type = typename ExtField::element_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = typename Base::matrix_type;
  using polynomial_type = typename Base::polynomial_type;

  // Constructs a BCH code over `base_field` with splitting field `ext_field`,
  // designed distance `designed_distance`, starting root exponent `b` (default
  // 1 for narrow-sense BCH), and block length `length` (default 0 => primitive
  // BCH of length |ext_field| - 1).
  BCHCode(std::shared_ptr<BaseField> base_field,
          std::shared_ptr<ExtField> ext_field, std::size_t designed_distance,
          std::size_t b = 1, std::size_t length = 0,
          Metric metric = Metric::HAMMING)
      : Base(base_field,
             detail::ResolveBchLength(base_field, ext_field, length),
             detail::BuildBchGeneratorPolynomial(
                 base_field, ext_field,
                 detail::ResolveBchLength(base_field, ext_field, length),
                 designed_distance, b),
             metric),
        ext_field_(std::move(ext_field)),
        designed_distance_(designed_distance),
        b_(b) {
    const uint64_t group_order = ext_field_->Order() - 1;
    const ext_storage_type gamma = ext_field_->MultiplicativeGenerator();
    primitive_nth_root_ = ext_field_->Pow(gamma, group_order / this->length_);

    this->RegisterDecoder("BCH", [](const AbstractCode<BaseField> &code) {
      const auto &bch =
          dynamic_cast<const BCHCode<BaseField, ExtField> &>(code);
      return std::make_shared<BCHDecoder<BaseField, ExtField>>(bch);
    });
    this->SetDefaultDecoder("BCH");
  }

  // Convenience factory specifying error-correction capability t
  // (designed distance delta = 2*t + 1).
  static BCHCode FromErrorCorrectionCapability(
      std::shared_ptr<BaseField> base_field,
      std::shared_ptr<ExtField> ext_field, std::size_t t, std::size_t b = 1,
      std::size_t length = 0, Metric metric = Metric::HAMMING) {
    if (t == 0) {
      throw std::invalid_argument(
          "BCHCode::FromErrorCorrectionCapability: t must be positive");
    }
    return BCHCode(std::move(base_field), std::move(ext_field), 2 * t + 1, b,
                   length, metric);
  }

  std::shared_ptr<ExtField> ExtensionField() const { return ext_field_; }

  std::size_t DesignedDistance() const noexcept { return designed_distance_; }

  std::size_t ErrorCorrectionCapability() const noexcept {
    return (designed_distance_ - 1) / 2;
  }

  std::size_t FirstConsecutiveRoot() const noexcept { return b_; }

  const ext_storage_type &PrimitiveNthRoot() const noexcept {
    return primitive_nth_root_;
  }

  bool IsBCH() const override { return true; }

  std::size_t MinimumDistance() const override {
    if (this->cached_min_distance_.has_value()) {
      return *this->cached_min_distance_;
    }
    // Exact search for small codes (q^k <= 4096); BCH bound otherwise.
    const uint64_t q = this->field_->Order();
    bool small_code = false;
    if (q != 0) {
      try {
        small_code = (utils::SafeIntegerPower(q, this->dimension_) <= 4096);
      } catch (const std::overflow_error &) {
        small_code = false;
      }
    }
    if (small_code) {
      return AbstractLinearCode<BaseField>::MinimumDistance();
    }
    this->cached_min_distance_ = designed_distance_;
    return designed_distance_;
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) +
           ", >=" + std::to_string(designed_distance_) + "] BCH code over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 private:
  std::shared_ptr<ExtField> ext_field_;
  std::size_t designed_distance_;
  std::size_t b_;
  ext_storage_type primitive_nth_root_;
};

template <typename BaseField, typename ExtField>
BCHDecoder<BaseField, ExtField>::BCHDecoder(
    const BCHCode<BaseField, ExtField> &code)
    : Base(code),
      ext_field_(code.ExtensionField()),
      n_(code.Length()),
      designed_distance_(code.DesignedDistance()),
      b_(code.FirstConsecutiveRoot()),
      t_(code.ErrorCorrectionCapability()) {
  const ext_storage_type &beta = code.PrimitiveNthRoot();
  const ext_storage_type ext_one = ext_field_->MultiplicativeIdentity();

  beta_pows_.resize(n_);
  beta_pows_[0] = ext_one;
  for (std::size_t i = 1; i < n_; ++i) {
    beta_pows_[i] = ext_field_->Mul(beta_pows_[i - 1], beta);
  }

  inv_beta_pows_.resize(n_);
  inv_beta_pows_[0] = ext_one;
  for (std::size_t i = 1; i < n_; ++i) {
    inv_beta_pows_[i] = beta_pows_[n_ - i];
  }

  const std::size_t num_syndromes = designed_distance_ - 1;
  syndrome_roots_.resize(num_syndromes);
  for (std::size_t k = 0; k < num_syndromes; ++k) {
    syndrome_roots_[k] = beta_pows_[(b_ + k) % n_];
  }
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_BCH_HPP_
