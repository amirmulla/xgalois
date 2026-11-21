#ifndef XGALOIS_POLY_POLY_DENSE_HPP_
#define XGALOIS_POLY_POLY_DENSE_HPP_

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/field/gf_element.hpp"

namespace xg {

template <typename GaloisField>
class PolynomialDense {
 private:
  struct InternalTag {};

 public:
  using FieldType = GaloisField;
  using field_type = GaloisField;
  using RawType = typename GaloisField::element_type;
  using storage_type = RawType;
  using ElementType = GaloisFieldElement<GaloisField>;

  static constexpr size_t kKaratsubaThreshold = 32;

  PolynomialDense() = default;

  explicit PolynomialDense(const std::vector<ElementType> &coeffs,
                           std::string var = "x") {
    if (coeffs.empty()) {
      throw std::invalid_argument("Coefficient vector cannot be empty.");
    }
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    if (!coeffs[0].HasField()) {
      throw std::invalid_argument(
          "Cannot determine Galois field from coefficients.");
    }
    std::shared_ptr<GaloisField> field = coeffs[0].Field();
    std::vector<RawType> raw;
    raw.reserve(coeffs.size());
    for (const auto &c : coeffs) {
      if (!c.HasField() || c.Field() != field) {
        throw std::invalid_argument(
            "All polynomial coefficients must belong to the same field.");
      }
      raw.push_back(c.RawValue());
    }
    coeffs_ = std::move(raw);
    field_ = std::move(field);
    var_ = std::move(var);
    Trim();
  }

  static PolynomialDense FromRaw(std::vector<RawType> raw_coeffs,
                                 std::shared_ptr<GaloisField> field,
                                 std::string var = "x") {
    if (!field) {
      throw std::invalid_argument("GaloisField pointer cannot be null.");
    }
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    return PolynomialDense(InternalTag{}, std::move(raw_coeffs),
                           std::move(field), std::move(var));
  }

  static PolynomialDense Zero(std::shared_ptr<GaloisField> field,
                              std::string var = "x") {
    if (!field) {
      throw std::invalid_argument("GaloisField pointer cannot be null.");
    }
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    return PolynomialDense(InternalTag{}, std::vector<RawType>{},
                           std::move(field), std::move(var));
  }

  static PolynomialDense One(std::shared_ptr<GaloisField> field,
                             std::string var = "x") {
    if (!field) {
      throw std::invalid_argument("GaloisField pointer cannot be null.");
    }
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    std::vector<RawType> raw;
    raw.push_back(field->MultiplicativeIdentity());
    return PolynomialDense(InternalTag{}, std::move(raw), std::move(field),
                           std::move(var));
  }

  static PolynomialDense Monomial(std::shared_ptr<GaloisField> field,
                                  size_t deg, RawType raw_coeff,
                                  std::string var = "x") {
    if (!field) {
      throw std::invalid_argument("GaloisField pointer cannot be null.");
    }
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    if (field->IsZero(raw_coeff)) {
      return PolynomialDense(InternalTag{}, std::vector<RawType>{},
                             std::move(field), std::move(var));
    }
    std::vector<RawType> raw(deg + 1, field->AdditiveIdentity());
    raw[deg] = std::move(raw_coeff);
    return PolynomialDense(InternalTag{}, std::move(raw), std::move(field),
                           std::move(var));
  }

  PolynomialDense(const PolynomialDense &) = default;
  PolynomialDense &operator=(const PolynomialDense &) = default;
  PolynomialDense(PolynomialDense &&) noexcept = default;
  PolynomialDense &operator=(PolynomialDense &&) noexcept = default;
  ~PolynomialDense() = default;

  bool HasField() const noexcept { return field_ != nullptr; }

  std::shared_ptr<GaloisField> Field() const { return field_; }

  void SetVariable(const std::string &var) {
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    var_ = var;
  }

  const std::string &GetVariable() const noexcept { return var_; }

  int Degree() const noexcept { return static_cast<int>(coeffs_.size()) - 1; }

  size_t Size() const noexcept { return coeffs_.size(); }

  bool IsZero() const noexcept { return coeffs_.empty(); }

  bool IsOne() const {
    return field_ != nullptr && coeffs_.size() == 1 &&
           field_->IsOne(coeffs_[0]);
  }

  bool IsMonic() const {
    return field_ != nullptr && !coeffs_.empty() &&
           field_->IsOne(coeffs_.back());
  }

  RawType RawCoeff(size_t i) const {
    if (i < coeffs_.size()) {
      return coeffs_[i];
    }
    RequireField();
    return field_->AdditiveIdentity();
  }

  const std::vector<RawType> &RawCoeffs() const noexcept { return coeffs_; }

  const ElementType operator[](size_t i) const {
    RequireField();
    return ElementType::FromInternal(RawCoeff(i), field_);
  }

  ElementType At(size_t i) const {
    RequireField();
    if (i >= coeffs_.size()) {
      throw std::out_of_range("Coefficient index out of bounds.");
    }
    return ElementType::FromInternal(coeffs_[i], field_);
  }

  void SetRawCoeff(size_t i, const RawType &raw_val) {
    RequireField();
    if (i < coeffs_.size()) {
      coeffs_[i] = raw_val;
      if (i + 1 == coeffs_.size()) {
        Trim();
      }
    } else if (!field_->IsZero(raw_val)) {
      coeffs_.resize(i + 1, field_->AdditiveIdentity());
      coeffs_[i] = raw_val;
    }
  }

  void SetCoeff(size_t i, const ElementType &val) {
    ValidateScalarField(val);
    SetRawCoeff(i, val.RawValue());
  }

  ElementType LeadingCoeff() const {
    RequireField();
    if (coeffs_.empty()) {
      return ElementType::FromInternal(field_->AdditiveIdentity(), field_);
    }
    return ElementType::FromInternal(coeffs_.back(), field_);
  }

  RawType LeadingRawCoeff() const {
    RequireField();
    if (coeffs_.empty()) {
      return field_->AdditiveIdentity();
    }
    return coeffs_.back();
  }

  PolynomialDense Monic() const {
    RequireField();
    if (coeffs_.empty()) {
      throw std::domain_error("Cannot normalize zero polynomial to monic.");
    }
    if (field_->IsOne(coeffs_.back())) {
      return *this;
    }
    const RawType inv_lc = field_->Inv(coeffs_.back());
    std::vector<RawType> out(coeffs_.size());
    field_->ScalarMulArray(std::span<RawType>(out),
                           std::span<const RawType>(coeffs_), inv_lc);
    out.back() = field_->MultiplicativeIdentity();
    return PolynomialDense(InternalTag{}, std::move(out), field_, var_);
  }

  PolynomialDense Shift(size_t k) const {
    RequireField();
    if (coeffs_.empty() || k == 0) {
      return *this;
    }
    std::vector<RawType> out(coeffs_.size() + k, field_->AdditiveIdentity());
    std::copy(coeffs_.begin(), coeffs_.end(),
              out.begin() + static_cast<std::ptrdiff_t>(k));
    return PolynomialDense(InternalTag{}, std::move(out), field_, var_);
  }

  // Returns x^deg * p(x^-1). Requires Degree() <= deg.
  PolynomialDense Reverse(size_t deg) const {
    RequireField();
    if (Degree() > static_cast<int>(deg)) {
      throw std::invalid_argument(
          "Reverse degree must be at least the polynomial degree.");
    }
    if (coeffs_.empty()) {
      return Zero(field_, var_);
    }
    std::vector<RawType> rev(deg + 1, field_->AdditiveIdentity());
    for (size_t i = 0; i < coeffs_.size(); ++i) {
      rev[deg - i] = coeffs_[i];
    }
    return PolynomialDense(InternalTag{}, std::move(rev), field_, var_);
  }

  PolynomialDense &operator+=(const PolynomialDense &other) {
    ValidateSameField(other);
    if (other.coeffs_.empty()) {
      return *this;
    }
    const size_t n = coeffs_.size();
    const size_t m = other.coeffs_.size();
    const size_t min_len = std::min(n, m);
    if (m > n) {
      coeffs_.resize(m);
      for (size_t i = min_len; i < m; ++i) {
        coeffs_[i] = other.coeffs_[i];
      }
    }
    if (min_len > 0) {
      field_->AddArray(std::span<RawType>(coeffs_.data(), min_len),
                       std::span<const RawType>(coeffs_.data(), min_len),
                       std::span<const RawType>(other.coeffs_.data(), min_len));
    }
    Trim();
    return *this;
  }

  PolynomialDense &operator-=(const PolynomialDense &other) {
    ValidateSameField(other);
    if (this == &other) {
      coeffs_.clear();
      return *this;
    }
    if (other.coeffs_.empty()) {
      return *this;
    }
    const size_t n = coeffs_.size();
    const size_t m = other.coeffs_.size();
    const size_t min_len = std::min(n, m);
    if (m > n) {
      coeffs_.resize(m);
      for (size_t i = min_len; i < m; ++i) {
        coeffs_[i] = field_->Neg(other.coeffs_[i]);
      }
    }
    if (min_len > 0) {
      field_->SubArray(std::span<RawType>(coeffs_.data(), min_len),
                       std::span<const RawType>(coeffs_.data(), min_len),
                       std::span<const RawType>(other.coeffs_.data(), min_len));
    }
    Trim();
    return *this;
  }

  PolynomialDense &operator*=(const PolynomialDense &other) {
    ValidateSameField(other);
    if (coeffs_.empty() || other.coeffs_.empty()) {
      coeffs_.clear();
      return *this;
    }
    coeffs_ = MulRaw(*field_, std::span<const RawType>(coeffs_),
                     std::span<const RawType>(other.coeffs_));
    Trim();
    return *this;
  }

  PolynomialDense &operator*=(const ElementType &scalar) {
    ValidateScalarField(scalar);
    const RawType s = scalar.RawValue();
    if (coeffs_.empty() || field_->IsZero(s)) {
      coeffs_.clear();
      return *this;
    }
    if (field_->IsOne(s)) {
      return *this;
    }
    field_->ScalarMulArray(std::span<RawType>(coeffs_),
                           std::span<const RawType>(coeffs_), s);
    Trim();
    return *this;
  }

  PolynomialDense &operator/=(const ElementType &scalar) {
    ValidateScalarField(scalar);
    const RawType s = scalar.RawValue();
    if (field_->IsZero(s)) {
      throw std::domain_error("Division by zero scalar.");
    }
    if (coeffs_.empty() || field_->IsOne(s)) {
      return *this;
    }
    const RawType inv_s = field_->Inv(s);
    field_->ScalarMulArray(std::span<RawType>(coeffs_),
                           std::span<const RawType>(coeffs_), inv_s);
    Trim();
    return *this;
  }

  PolynomialDense &operator/=(const PolynomialDense &divisor) {
    auto [q, r] = DivRem(divisor);
    (void)r;
    *this = std::move(q);
    return *this;
  }

  PolynomialDense &operator%=(const PolynomialDense &divisor) {
    ValidateSameField(divisor);
    if (divisor.coeffs_.empty()) {
      throw std::invalid_argument("Division by zero polynomial.");
    }
    if (this == &divisor) {
      coeffs_.clear();
      return *this;
    }
    const int deg_a = Degree();
    const int deg_b = divisor.Degree();
    if (deg_a < deg_b) {
      return *this;
    }
    if (deg_b == 0) {
      coeffs_.clear();
      return *this;
    }
    const RawType inv_lc = field_->Inv(divisor.coeffs_.back());
    const bool monic_div = field_->IsOne(inv_lc);
    const RawType zero = field_->AdditiveIdentity();
    const size_t db = static_cast<size_t>(deg_b);
    for (int i = deg_a; i >= deg_b; --i) {
      const size_t idx = static_cast<size_t>(i);
      const RawType lead = coeffs_[idx];
      if (field_->IsZero(lead)) {
        continue;
      }
      const RawType q_coeff = monic_div ? lead : field_->Mul(lead, inv_lc);
      const RawType neg_q = field_->Neg(q_coeff);
      const size_t shift = idx - db;
      field_->FmaArray(std::span<RawType>(coeffs_.data() + shift, db),
                       std::span<const RawType>(divisor.coeffs_.data(), db),
                       neg_q);
      coeffs_[idx] = zero;
    }
    coeffs_.resize(db);
    Trim();
    return *this;
  }

  PolynomialDense operator+(const PolynomialDense &other) const {
    PolynomialDense result(*this);
    result += other;
    return result;
  }

  PolynomialDense operator-(const PolynomialDense &other) const {
    PolynomialDense result(*this);
    result -= other;
    return result;
  }

  PolynomialDense operator*(const PolynomialDense &other) const {
    ValidateSameField(other);
    if (coeffs_.empty() || other.coeffs_.empty()) {
      return Zero(field_, var_);
    }
    auto prod = MulRaw(*field_, std::span<const RawType>(coeffs_),
                       std::span<const RawType>(other.coeffs_));
    return PolynomialDense(InternalTag{}, std::move(prod), field_, var_);
  }

  PolynomialDense operator*(const ElementType &scalar) const {
    ValidateScalarField(scalar);
    const RawType s = scalar.RawValue();
    if (coeffs_.empty() || field_->IsZero(s)) {
      return Zero(field_, var_);
    }
    if (field_->IsOne(s)) {
      return *this;
    }
    std::vector<RawType> out(coeffs_.size());
    field_->ScalarMulArray(std::span<RawType>(out),
                           std::span<const RawType>(coeffs_), s);
    return PolynomialDense(InternalTag{}, std::move(out), field_, var_);
  }

  PolynomialDense operator/(const ElementType &scalar) const {
    ValidateScalarField(scalar);
    const RawType s = scalar.RawValue();
    if (field_->IsZero(s)) {
      throw std::domain_error("Division by zero scalar.");
    }
    if (coeffs_.empty() || field_->IsOne(s)) {
      return *this;
    }
    const RawType inv_s = field_->Inv(s);
    std::vector<RawType> out(coeffs_.size());
    field_->ScalarMulArray(std::span<RawType>(out),
                           std::span<const RawType>(coeffs_), inv_s);
    return PolynomialDense(InternalTag{}, std::move(out), field_, var_);
  }

  PolynomialDense operator/(const PolynomialDense &divisor) const {
    auto [q, r] = DivRem(divisor);
    (void)r;
    return q;
  }

  PolynomialDense operator%(const PolynomialDense &divisor) const {
    PolynomialDense result(*this);
    result %= divisor;
    return result;
  }

  PolynomialDense operator-() const {
    RequireField();
    if (coeffs_.empty()) {
      return *this;
    }
    std::vector<RawType> out(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i) {
      out[i] = field_->Neg(coeffs_[i]);
    }
    return PolynomialDense(InternalTag{}, std::move(out), field_, var_);
  }

  std::pair<PolynomialDense, PolynomialDense> DivRem(
      const PolynomialDense &divisor) const {
    ValidateSameField(divisor);
    if (divisor.coeffs_.empty()) {
      throw std::invalid_argument("Division by zero polynomial.");
    }
    const int deg_a = Degree();
    const int deg_b = divisor.Degree();
    if (deg_a < deg_b) {
      return {Zero(field_, var_), *this};
    }

    const RawType zero = field_->AdditiveIdentity();
    const RawType inv_lc = field_->Inv(divisor.coeffs_.back());
    const bool monic_div = field_->IsOne(inv_lc);

    if (deg_b == 0) {
      std::vector<RawType> quot(coeffs_.size());
      if (monic_div) {
        quot = coeffs_;
      } else {
        field_->ScalarMulArray(std::span<RawType>(quot),
                               std::span<const RawType>(coeffs_), inv_lc);
      }
      return {PolynomialDense(InternalTag{}, std::move(quot), field_, var_),
              Zero(field_, var_)};
    }

    const size_t db = static_cast<size_t>(deg_b);
    std::vector<RawType> rem = coeffs_;
    std::vector<RawType> quot(static_cast<size_t>(deg_a - deg_b + 1), zero);

    for (int i = deg_a; i >= deg_b; --i) {
      const size_t idx = static_cast<size_t>(i);
      const RawType lead = rem[idx];
      if (field_->IsZero(lead)) {
        continue;
      }
      const RawType q_coeff = monic_div ? lead : field_->Mul(lead, inv_lc);
      const size_t shift = idx - db;
      quot[shift] = q_coeff;
      const RawType neg_q = field_->Neg(q_coeff);
      field_->FmaArray(std::span<RawType>(rem.data() + shift, db),
                       std::span<const RawType>(divisor.coeffs_.data(), db),
                       neg_q);
      rem[idx] = zero;
    }

    rem.resize(db);
    return {PolynomialDense(InternalTag{}, std::move(quot), field_, var_),
            PolynomialDense(InternalTag{}, std::move(rem), field_, var_)};
  }

  PolynomialDense Pow(uint64_t exp) const {
    RequireField();
    if (exp == 0) {
      return One(field_, var_);
    }
    if (coeffs_.empty()) {
      return Zero(field_, var_);
    }
    if (exp == 1) {
      return *this;
    }
    PolynomialDense result = One(field_, var_);
    PolynomialDense base = *this;
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        result *= base;
      }
      exp >>= 1;
      if (exp > 0) {
        base *= base;
      }
    }
    return result;
  }

  PolynomialDense operator^(uint64_t exp) const { return Pow(exp); }

  PolynomialDense PowMod(uint64_t exp, const PolynomialDense &mod) const {
    ValidateSameField(mod);
    if (mod.coeffs_.empty()) {
      throw std::invalid_argument("Modulus polynomial cannot be zero.");
    }
    if (mod.Degree() == 0) {
      return Zero(field_, var_);
    }
    if (exp == 0) {
      return One(field_, var_);
    }
    PolynomialDense result = One(field_, var_);
    PolynomialDense base = (*this) % mod;
    if (base.IsZero()) {
      return Zero(field_, var_);
    }
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        result *= base;
        result %= mod;
      }
      exp >>= 1;
      if (exp > 0) {
        base *= base;
        base %= mod;
      }
    }
    return result;
  }

  // O(d) formal derivative valid across all characteristics and field types.
  PolynomialDense Derivative() const {
    RequireField();
    if (coeffs_.size() <= 1) {
      return Zero(field_, var_);
    }
    const size_t n = coeffs_.size();
    const RawType zero = field_->AdditiveIdentity();
    std::vector<RawType> deriv(n - 1, zero);

    if (field_->Characteristic() == 2) {
      for (size_t i = 1; i < n; i += 2) {
        deriv[i - 1] = coeffs_[i];
      }
    } else {
      const RawType one = field_->MultiplicativeIdentity();
      RawType mult = zero;
      for (size_t i = 1; i < n; ++i) {
        mult = field_->Add(mult, one);
        if (!field_->IsZero(mult) && !field_->IsZero(coeffs_[i])) {
          deriv[i - 1] =
              field_->IsOne(mult) ? coeffs_[i] : field_->Mul(coeffs_[i], mult);
        }
      }
    }
    return PolynomialDense(InternalTag{}, std::move(deriv), field_, var_);
  }

  RawType EvalRaw(const RawType &raw_x) const {
    RequireField();
    if (coeffs_.empty()) {
      return field_->AdditiveIdentity();
    }
    RawType acc = coeffs_.back();
    for (size_t i = coeffs_.size() - 1; i > 0; --i) {
      acc = field_->Add(field_->Mul(acc, raw_x), coeffs_[i - 1]);
    }
    return acc;
  }

  std::vector<RawType> EvalBatch(std::span<const RawType> xs) const {
    RequireField();
    std::vector<RawType> out(xs.size(), field_->AdditiveIdentity());
    if (coeffs_.empty() || xs.empty()) {
      return out;
    }
    std::fill(out.begin(), out.end(), coeffs_.back());
    for (size_t i = coeffs_.size() - 1; i > 0; --i) {
      field_->MulArray(std::span<RawType>(out), std::span<const RawType>(out),
                       xs);
      const RawType &c = coeffs_[i - 1];
      if (!field_->IsZero(c)) {
        for (size_t j = 0; j < out.size(); ++j) {
          out[j] = field_->Add(out[j], c);
        }
      }
    }
    return out;
  }

  ElementType operator()(const ElementType &x) const {
    ValidateScalarField(x);
    return ElementType::FromInternal(EvalRaw(x.RawValue()), field_);
  }

  // Returns p(q(x)) using Horner's method.
  PolynomialDense Compose(const PolynomialDense &q) const {
    ValidateSameField(q);
    if (coeffs_.empty()) {
      return Zero(field_, q.var_);
    }
    PolynomialDense acc = PolynomialDense(
        InternalTag{}, std::vector<RawType>{coeffs_.back()}, field_, q.var_);
    for (size_t i = coeffs_.size() - 1; i > 0; --i) {
      acc *= q;
      const RawType &c = coeffs_[i - 1];
      if (!field_->IsZero(c)) {
        acc.AddConstantInPlace(c);
      }
    }
    return acc;
  }

  bool operator==(const PolynomialDense &other) const {
    if (field_ != other.field_) {
      return false;
    }
    return coeffs_ == other.coeffs_;
  }

  bool operator!=(const PolynomialDense &other) const {
    return !(*this == other);
  }

  void Trim() {
    if (!field_) {
      coeffs_.clear();
      return;
    }
    while (!coeffs_.empty() && field_->IsZero(coeffs_.back())) {
      coeffs_.pop_back();
    }
  }

  void Print(std::ostream &os) const {
    if (coeffs_.empty() || !field_) {
      os << "0";
      return;
    }

    auto needs_brackets = [](const std::string &s) -> bool {
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

    auto needs_const_brackets = [](const std::string &s) -> bool {
      for (char ch : s) {
        if (ch == '+' || ch == '-' ||
            std::isspace(static_cast<unsigned char>(ch)) != 0) {
          return true;
        }
      }
      return false;
    };

    bool first_term = true;
    for (size_t idx = coeffs_.size(); idx > 0; --idx) {
      const size_t i = idx - 1;
      const RawType &raw = coeffs_[i];
      if (field_->IsZero(raw)) {
        continue;
      }
      if (!first_term) {
        os << " + ";
      }
      if (i == 0) {
        const std::string coeff_str = field_->ToString(raw);
        if (!first_term && needs_const_brackets(coeff_str)) {
          os << "(" << coeff_str << ")";
        } else {
          os << coeff_str;
        }
      } else {
        if (!field_->IsOne(raw)) {
          const std::string coeff_str = field_->ToString(raw);
          if (needs_brackets(coeff_str)) {
            os << "(" << coeff_str << ")";
          } else {
            os << coeff_str;
          }
        }
        os << var_;
        if (i > 1) {
          os << "^" << i;
        }
      }
      first_term = false;
    }

    if (first_term) {
      os << "0";
    }
  }

  std::string ToString() const {
    std::ostringstream oss;
    Print(oss);
    return oss.str();
  }

 private:
  PolynomialDense(InternalTag, std::vector<RawType> raw_coeffs,
                  std::shared_ptr<GaloisField> field, std::string var)
      : coeffs_(std::move(raw_coeffs)),
        field_(std::move(field)),
        var_(std::move(var)) {
    Trim();
  }

  void RequireField() const {
    if (!field_) {
      throw std::runtime_error("Polynomial has no associated Galois field.");
    }
  }

  void ValidateSameField(const PolynomialDense &other) const {
    if (!field_ || !other.field_) {
      throw std::invalid_argument(
          "Operations on polynomials with null field are not allowed.");
    }
    if (field_ != other.field_) {
      throw std::invalid_argument(
          "Operations between polynomials over different fields are not "
          "allowed.");
    }
  }

  void ValidateScalarField(const ElementType &scalar) const {
    if (!field_ || !scalar.HasField()) {
      throw std::invalid_argument(
          "Operations with null field are not allowed.");
    }
    if (field_ != scalar.Field()) {
      throw std::invalid_argument(
          "Scalar belongs to a different field than the polynomial.");
    }
  }

  void AddConstantInPlace(const RawType &c) {
    if (coeffs_.empty()) {
      coeffs_.push_back(c);
      return;
    }
    coeffs_[0] = field_->Add(coeffs_[0], c);
    if (coeffs_.size() == 1) {
      Trim();
    }
  }

  static std::vector<RawType> MulSchoolbook(const GaloisField &field,
                                            std::span<const RawType> a,
                                            std::span<const RawType> b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    const RawType zero = field.AdditiveIdentity();
    std::vector<RawType> out(a.size() + b.size() - 1, zero);
    if (a.size() <= b.size()) {
      for (size_t i = 0; i < a.size(); ++i) {
        if (!field.IsZero(a[i])) {
          field.FmaArray(std::span<RawType>(out.data() + i, b.size()), b, a[i]);
        }
      }
    } else {
      for (size_t j = 0; j < b.size(); ++j) {
        if (!field.IsZero(b[j])) {
          field.FmaArray(std::span<RawType>(out.data() + j, a.size()), a, b[j]);
        }
      }
    }
    return out;
  }

  static std::vector<RawType> MulRaw(const GaloisField &field,
                                     std::span<const RawType> a,
                                     std::span<const RawType> b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    if (a.size() <= kKaratsubaThreshold || b.size() <= kKaratsubaThreshold) {
      return MulSchoolbook(field, a, b);
    }

    // Handle strongly unbalanced operands by chunking the larger operand.
    if (a.size() >= 2 * b.size()) {
      const RawType zero = field.AdditiveIdentity();
      std::vector<RawType> out(a.size() + b.size() - 1, zero);
      const size_t chunk = b.size();
      for (size_t offset = 0; offset < a.size(); offset += chunk) {
        const size_t len = std::min(chunk, a.size() - offset);
        auto prod = MulRaw(field, a.subspan(offset, len), b);
        field.AddArray(
            std::span<RawType>(out.data() + offset, prod.size()),
            std::span<const RawType>(out.data() + offset, prod.size()),
            std::span<const RawType>(prod));
      }
      return out;
    }
    if (b.size() >= 2 * a.size()) {
      return MulRaw(field, b, a);
    }

    // Balanced Karatsuba split: max < 2 * min implies 0 < m < min.
    const size_t m = std::max(a.size(), b.size()) / 2;
    const auto a0 = a.subspan(0, m);
    const auto a1 = a.subspan(m);
    const auto b0 = b.subspan(0, m);
    const auto b1 = b.subspan(m);

    auto z0 = MulRaw(field, a0, b0);
    auto z2 = MulRaw(field, a1, b1);

    const size_t a_sum_len = std::max(a0.size(), a1.size());
    const size_t b_sum_len = std::max(b0.size(), b1.size());
    std::vector<RawType> a_sum(a_sum_len);
    std::vector<RawType> b_sum(b_sum_len);

    {
      const size_t min_len = std::min(a0.size(), a1.size());
      field.AddArray(std::span<RawType>(a_sum.data(), min_len),
                     a0.subspan(0, min_len), a1.subspan(0, min_len));
      for (size_t i = min_len; i < a0.size(); ++i) {
        a_sum[i] = a0[i];
      }
      for (size_t i = min_len; i < a1.size(); ++i) {
        a_sum[i] = a1[i];
      }
    }
    {
      const size_t min_len = std::min(b0.size(), b1.size());
      field.AddArray(std::span<RawType>(b_sum.data(), min_len),
                     b0.subspan(0, min_len), b1.subspan(0, min_len));
      for (size_t i = min_len; i < b0.size(); ++i) {
        b_sum[i] = b0[i];
      }
      for (size_t i = min_len; i < b1.size(); ++i) {
        b_sum[i] = b1[i];
      }
    }

    auto z1 = MulRaw(field, a_sum, b_sum);
    field.SubArray(std::span<RawType>(z1.data(), z0.size()),
                   std::span<const RawType>(z1.data(), z0.size()),
                   std::span<const RawType>(z0));
    field.SubArray(std::span<RawType>(z1.data(), z2.size()),
                   std::span<const RawType>(z1.data(), z2.size()),
                   std::span<const RawType>(z2));

    while (!z1.empty() && field.IsZero(z1.back())) {
      z1.pop_back();
    }

    const RawType zero = field.AdditiveIdentity();
    std::vector<RawType> out(a.size() + b.size() - 1, zero);
    std::copy(z0.begin(), z0.end(), out.begin());
    if (!z1.empty()) {
      const size_t z1_len = std::min(z1.size(), out.size() - m);
      field.AddArray(std::span<RawType>(out.data() + m, z1_len),
                     std::span<const RawType>(out.data() + m, z1_len),
                     std::span<const RawType>(z1.data(), z1_len));
    }
    if (!z2.empty()) {
      field.AddArray(std::span<RawType>(out.data() + 2 * m, z2.size()),
                     std::span<const RawType>(out.data() + 2 * m, z2.size()),
                     std::span<const RawType>(z2));
    }
    return out;
  }

  std::vector<RawType> coeffs_;
  std::shared_ptr<GaloisField> field_{nullptr};
  std::string var_ = "x";
};

template <typename GaloisField>
inline PolynomialDense<GaloisField> operator*(
    const GaloisFieldElement<GaloisField> &scalar,
    const PolynomialDense<GaloisField> &poly) {
  return poly * scalar;
}

template <typename GaloisField>
inline std::ostream &operator<<(std::ostream &os,
                                const PolynomialDense<GaloisField> &poly) {
  poly.Print(os);
  return os;
}

}  // namespace xg

#endif  // XGALOIS_POLY_POLY_DENSE_HPP_
