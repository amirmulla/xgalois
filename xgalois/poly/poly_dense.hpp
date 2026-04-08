#ifndef XG_POLYNOMIAL_DENSE_HPP_
#define XG_POLYNOMIAL_DENSE_HPP_

#include <algorithm>
#include <cassert>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "xgalois/field/gf_element.hpp"

namespace xg {

template <typename GaloisField>
class PolynomialDense {
 public:
  using ElementType = GaloisFieldElementBase<GaloisField>;

  explicit PolynomialDense(const std::vector<ElementType> &coeffs,
                           const std::string &var = "x")
      : coefficients_(coeffs),
        field_(coeffs[0].Field()),
        kGfZero(field_->AdditiveIdentity(), field_),
        kGfOne(field_->MultiplicativeIdentity(), field_),
        variable_(var) {
    assert(!coeffs.empty() && "Coefficient vector cannot be empty.");
    if (!field_) {
      throw std::invalid_argument(
          "Cannot determine Galois field from coefficients.");
    }
    Trim();
  }

  PolynomialDense(const PolynomialDense &other) = default;

  PolynomialDense &operator=(const PolynomialDense &other) = default;

  PolynomialDense(PolynomialDense &&other) = default;

  PolynomialDense &operator=(PolynomialDense &&other) = default;

  ~PolynomialDense() = default;

  void SetVariable(const std::string &var) {
    if (var.empty()) {
      throw std::invalid_argument("Variable string cannot be empty.");
    }
    variable_ = var;
  }

  std::string GetVariable() const { return variable_; }

  int Degree() const {
    if (coefficients_.size() == 1 && coefficients_[0] == kGfZero) {
      return -1;
    }
    return coefficients_.size() - 1;
  }

  size_t Size() const { return coefficients_.size(); }

  const ElementType &operator[](int index) const {
    if (index < 0 || index >= coefficients_.size()) {
      throw std::out_of_range("Coefficient index out of bounds.");
    }
    return coefficients_[index];
  }

  ElementType &operator[](int index) {
    if (index < 0 || index >= coefficients_.size()) {
      throw std::out_of_range("Coefficient index out of bounds.");
    }
    return coefficients_[index];
  }

  PolynomialDense<GaloisField> operator+(
      const PolynomialDense<GaloisField> &other) const {
    size_t result_size = std::max(Size(), other.Size());
    std::vector<ElementType> result_coeffs(result_size, kGfZero);

    size_t min_size = std::min(Size(), other.Size());
    for (size_t i = 0; i < min_size; ++i) {
      result_coeffs[i] = coefficients_[i] + other[i];
    }

    if (Size() > other.Size()) {
      for (size_t i = min_size; i < Size(); ++i) {
        result_coeffs[i] = coefficients_[i];
      }
    } else if (other.Size() > Size()) {
      for (size_t i = min_size; i < other.Size(); ++i) {
        result_coeffs[i] = other[i];
      }
    }

    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();
    return result;
  }

  PolynomialDense<GaloisField> operator-(
      const PolynomialDense<GaloisField> &other) const {
    size_t result_size = std::max(Size(), other.Size());
    std::vector<ElementType> result_coeffs(result_size, kGfZero);

    size_t min_size = std::min(Size(), other.Size());
    for (size_t i = 0; i < min_size; ++i) {
      result_coeffs[i] = coefficients_[i] - other[i];
    }

    if (Size() > other.Size()) {
      for (size_t i = min_size; i < Size(); ++i) {
        result_coeffs[i] = coefficients_[i];
      }
    } else if (other.Size() > Size()) {
      for (size_t i = min_size; i < other.Size(); ++i) {
        result_coeffs[i] = kGfZero - other[i];
      }
    }

    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();
    return result;
  }

  PolynomialDense<GaloisField> operator*(
      const PolynomialDense<GaloisField> &other) const {
    if ((Size() == 1 && coefficients_[0] == kGfZero) ||
        (other.Size() == 1 && other[0] == kGfZero)) {
      return PolynomialDense<GaloisField>(std::vector<ElementType>{kGfZero},
                                          variable_);
    }

    size_t result_size = (Size() - 1) + (other.Size() - 1) + 1;
    std::vector<ElementType> result_coeffs(result_size, kGfZero);

    for (size_t i = 0; i < Size(); ++i) {
      for (size_t j = 0; j < other.Size(); ++j) {
        result_coeffs[i + j] += coefficients_[i] * other[j];
      }
    }

    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();
    return result;
  }

  PolynomialDense<GaloisField> operator/(
      const PolynomialDense<GaloisField> &other) const {
    auto result_pair = DivRem(other);
    return std::move(result_pair.first);
  }

  PolynomialDense<GaloisField> operator%(
      const PolynomialDense<GaloisField> &other) const {
    auto result_pair = DivRem(other);
    return std::move(result_pair.second);
  }

  PolynomialDense<GaloisField> operator*(const ElementType &scalar) const {
    std::vector<ElementType> result_coeffs(coefficients_.size());
    for (size_t i = 0; i < coefficients_.size(); ++i) {
      result_coeffs[i] = coefficients_[i] * scalar;
    }
    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();
    return result;
  }

  PolynomialDense<GaloisField> operator-() const {
    std::vector<ElementType> result_coeffs(coefficients_.size());
    for (size_t i = 0; i < coefficients_.size(); ++i) {
      result_coeffs[i] = kGfZero - coefficients_[i];
    }
    PolynomialDense<GaloisField> result(result_coeffs, variable_);
    result.Trim();
    return result;
  }

  PolynomialDense<GaloisField> operator^(uint64_t exponent) const {
    if (exponent == 0) {
      return PolynomialDense<GaloisField>(std::vector<ElementType>{kGfOne},
                                          variable_);
    }
    if (exponent == 1) {
      return *this;
    }

    PolynomialDense<GaloisField> base = *this;
    PolynomialDense<GaloisField> result(std::vector<ElementType>{kGfOne},
                                        variable_);

    uint64_t temp_exponent = exponent;

    while (temp_exponent > 0) {
      if (temp_exponent % 2 == 1) {
        result = result * base;
      }
      base = base * base;
      temp_exponent /= 2;
    }

    return result;
  }

  PolynomialDense<GaloisField> &operator+=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this + other;
    return *this;
  }

  PolynomialDense<GaloisField> &operator-=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this - other;
    return *this;
  }

  PolynomialDense<GaloisField> &operator*=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this * other;
    return *this;
  }

  PolynomialDense<GaloisField> &operator*=(const ElementType &scalar) {
    *this = *this * scalar;
    return *this;
  }

  PolynomialDense<GaloisField> &operator/=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this / other;
    return *this;
  }

  PolynomialDense<GaloisField> &operator%=(
      const PolynomialDense<GaloisField> &other) {
    *this = *this % other;
    return *this;
  }

  bool operator==(const PolynomialDense<GaloisField> &other) const {
    return coefficients_ == other.coefficients_;
  }

  bool operator!=(const PolynomialDense<GaloisField> &other) const {
    return !(*this == other);
  }

  ElementType operator()(const ElementType &x) const {
    ElementType result = kGfZero;
    for (int i = coefficients_.size() - 1; i >= 0; --i) {
      result = result * x + coefficients_[i];
    }
    return result;
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }

  void Print(std::ostream &os) const {
    if (coefficients_.empty() ||
        (coefficients_.size() == 1 && coefficients_[0] == kGfZero)) {
      os << "0";
      return;
    }

    auto needs_brackets = [](const ElementType &coeff) {
      std::stringstream ss;
      ss << coeff;
      std::string coeff_str = ss.str();

      for (char c : coeff_str) {
        if (!std::isdigit(c)) {
          return true;
        }
      }
      return false;
    };

    bool first_term = true;
    for (int i = static_cast<int>(coefficients_.size()) - 1; i >= 0; --i) {
      const ElementType &coeff = coefficients_[i];

      if (coeff == kGfZero && coefficients_.size() > 1) {
        continue;
      }

      if (!first_term) {
        os << " + ";
      }

      if (i == 0 || coeff != kGfOne ||
          (coeff == kGfOne && i == 0 && coefficients_.size() == 1)) {
        if (needs_brackets(coeff)) {
          os << "(" << coeff << ")";
        } else {
          os << coeff;
        }
      } else if (coeff == kGfOne && i > 0 && first_term &&
                 coefficients_.size() - 1 == i) {
      } else if (coeff != kGfOne) {
        os << coeff;
      }

      if (i > 0) {
        os << variable_;
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

  std::pair<PolynomialDense<GaloisField>, PolynomialDense<GaloisField>> DivRem(
      const PolynomialDense<GaloisField> &divisor) const {
    if (divisor.Size() == 1 && divisor[0] == kGfZero) {
      throw std::invalid_argument("Division by zero polynomial.");
    }

    PolynomialDense<GaloisField> quotient(std::vector<ElementType>{kGfZero},
                                          variable_);
    PolynomialDense<GaloisField> remainder = *this;

    int divisor_degree = divisor.Degree();

    if (remainder.Degree() < divisor_degree) {
      return {std::move(quotient), std::move(remainder)};
    }

    while (remainder.Degree() >= divisor_degree && remainder.Degree() != -1) {
      int remainder_degree = remainder.Degree();
      int degree_difference = remainder_degree - divisor_degree;

      const ElementType &remainder_leading_coeff = remainder[remainder_degree];
      const ElementType &divisor_leading_coeff = divisor[divisor_degree];

      if (divisor_leading_coeff == kGfZero) {
        throw std::runtime_error(
            "Internal error: Divisor leading coefficient "
            "is zero during division.");
      }

      ElementType term_coeff = remainder_leading_coeff / divisor_leading_coeff;

      std::vector<ElementType> term_coeffs(degree_difference + 1, kGfZero);
      term_coeffs[degree_difference] = term_coeff;

      PolynomialDense<GaloisField> term_poly(term_coeffs, variable_);

      quotient += term_poly;

      PolynomialDense<GaloisField> product = term_poly * divisor;
      remainder -= product;

      remainder.Trim();
    }

    return {std::move(quotient), std::move(remainder)};
  }

  PolynomialDense<GaloisField> Derivative() const {
    if (coefficients_.size() <= 1) {
      return PolynomialDense<GaloisField>(std::vector<ElementType>{kGfZero},
                                          variable_);
    }

    std::vector<ElementType> derivative_coeffs;
    derivative_coeffs.reserve(coefficients_.size() - 1);

    for (size_t i = 1; i < coefficients_.size(); ++i) {
      ElementType coeff = coefficients_[i];

      ElementType derivative_coeff = kGfZero;
      for (size_t j = 0; j < i; ++j) {
        derivative_coeff += coeff;
      }

      derivative_coeffs.push_back(derivative_coeff);
    }

    if (derivative_coeffs.empty()) {
      derivative_coeffs.push_back(kGfZero);
    }

    PolynomialDense<GaloisField> result(derivative_coeffs, variable_);
    result.Trim();
    return result;
  }

  void Trim() {
    int last_nonzero = coefficients_.size() - 1;
    while (last_nonzero >= 0 && coefficients_[last_nonzero] == kGfZero) {
      last_nonzero--;
    }

    if (last_nonzero < 0) {
      coefficients_ = {kGfZero};
    } else {
      coefficients_.resize(last_nonzero + 1);
    }
  }

 private:
  std::vector<ElementType> coefficients_;

  std::shared_ptr<GaloisField> field_;

  ElementType kGfZero;

  ElementType kGfOne;

  std::string variable_;
};

template <typename GaloisField>
std::ostream &operator<<(std::ostream &os,
                         const PolynomialDense<GaloisField> &poly) {
  poly.Print(os);
  return os;
}

}  // namespace xg

#endif
