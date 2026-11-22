#ifndef XGALOIS_FIELD_GF_TOWER_HPP_
#define XGALOIS_FIELD_GF_TOWER_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {

// Element of a relative field extension K = F[x]/(f(x)) of degree m.
// Uses fixed-size std::array storage when StaticDegree > 0, or std::vector
// when StaticDegree == 0 (runtime extension degree).
template <typename BaseElem, std::size_t StaticDegree = 0>
struct TowerElement {
  using base_element_type = BaseElem;
  using storage_type =
      std::conditional_t<(StaticDegree > 0), std::array<BaseElem, StaticDegree>,
                         std::vector<BaseElem>>;

  storage_type coeffs{};

  constexpr TowerElement() = default;

  constexpr explicit TowerElement(const storage_type& c) : coeffs(c) {}

  constexpr explicit TowerElement(storage_type&& c) noexcept(
      std::is_nothrow_move_constructible_v<storage_type>)
      : coeffs(std::move(c)) {}

  TowerElement(std::initializer_list<BaseElem> init) {
    if constexpr (StaticDegree > 0) {
      if (init.size() > StaticDegree) {
        throw std::invalid_argument(
            "TowerElement initializer exceeds StaticDegree");
      }
      std::size_t i = 0;
      for (const BaseElem& v : init) {
        coeffs[i++] = v;
      }
    } else {
      coeffs.assign(init.begin(), init.end());
    }
  }

  template <std::size_t D = StaticDegree>
    requires(D > 0)
  explicit TowerElement(const std::vector<BaseElem>& vec) {
    if (vec.size() > StaticDegree) {
      throw std::invalid_argument("TowerElement vector exceeds StaticDegree");
    }
    for (std::size_t i = 0; i < vec.size(); ++i) {
      coeffs[i] = vec[i];
    }
  }

  constexpr const storage_type& Coeffs() const noexcept { return coeffs; }
  constexpr storage_type& Coeffs() noexcept { return coeffs; }

  constexpr std::size_t Size() const noexcept { return coeffs.size(); }

  constexpr const BaseElem& operator[](std::size_t index) const {
    return coeffs[index];
  }

  constexpr BaseElem& operator[](std::size_t index) { return coeffs[index]; }

  constexpr bool operator==(const TowerElement& other) const = default;
};

namespace detail {

template <typename F>
BigInt QueryFieldBigCharacteristic(const F& field) {
  if constexpr (requires { field.BigCharacteristic(); }) {
    return field.BigCharacteristic();
  } else if constexpr (requires { field.Characteristic(); }) {
    return BigInt::FromUint64(static_cast<uint64_t>(field.Characteristic()));
  } else {
    return BigInt(0);
  }
}

template <typename F>
uint64_t QueryFieldCharacteristicU64(const F& field) {
  if constexpr (requires { field.Characteristic(); }) {
    return static_cast<uint64_t>(field.Characteristic());
  } else if constexpr (requires { field.BigCharacteristic(); }) {
    const BigInt ch = field.BigCharacteristic();
    if (!ch.IsNegative() && ch.LimbCount() <= 1) {
      return ch.ToUint64();
    }
    return 0;
  } else {
    return 0;
  }
}

template <typename F>
std::size_t QueryFieldAbsoluteDegree(const F& field) {
  if constexpr (requires { field.AbsoluteDegree(); }) {
    return static_cast<std::size_t>(field.AbsoluteDegree());
  } else if constexpr (requires { field.Degree(); }) {
    return static_cast<std::size_t>(field.Degree());
  } else {
    return 1;
  }
}

template <typename F>
BigInt QueryFieldBigOrder(const F& field) {
  if constexpr (requires { field.BigOrder(); }) {
    return field.BigOrder();
  } else if constexpr (requires { field.Order(); }) {
    const uint64_t ord = static_cast<uint64_t>(field.Order());
    if (ord > 0) {
      return BigInt::FromUint64(ord);
    }
    // Order() == 0 denotes 2^64 in GF(2^64) or overflow; derive p^deg.
    const BigInt ch = QueryFieldBigCharacteristic(field);
    if (ch > BigInt(0)) {
      return ch.Pow(static_cast<uint64_t>(QueryFieldAbsoluteDegree(field)));
    }
    return BigInt(0);
  } else {
    return BigInt(0);
  }
}

}  // namespace detail

// Relative field extension K = F[x]/(f(x)) over any base field F satisfying
// xg::Field, including another FieldTower (supports recursive extension towers
// such as F_{p} -> F_{p^2} -> F_{p^6} -> F_{p^{12}}).
template <Field BaseFieldT, std::size_t StaticDegree = 0>
class FieldTower {
 public:
  using base_field_type = BaseFieldT;
  using base_element_type = typename BaseFieldT::element_type;
  using element_type = TowerElement<base_element_type, StaticDegree>;
  using storage_type = element_type;

  // Constructs K = base_field[x] / (modulus_coeffs), where modulus_coeffs is
  // in ascending degree order: c_0 + c_1 x + ... + c_m x^m (normalized to
  // monic automatically).
  FieldTower(BaseFieldT base_field,
             std::vector<base_element_type> modulus_coeffs,
             std::string var_name = "u")
      : base_(std::move(base_field)),
        var_name_(var_name.empty() ? "u" : std::move(var_name)) {
    // Strip trailing zeros from modulus.
    while (!modulus_coeffs.empty() && base_.IsZero(modulus_coeffs.back())) {
      modulus_coeffs.pop_back();
    }
    if (modulus_coeffs.size() < 2) {
      throw std::invalid_argument(
          "FieldTower requires a modulus polynomial of degree >= 1");
    }

    degree_ = modulus_coeffs.size() - 1;
    if constexpr (StaticDegree > 0) {
      if (degree_ != StaticDegree) {
        throw std::invalid_argument(
            "FieldTower modulus degree does not match StaticDegree");
      }
    }

    // Normalize modulus to monic.
    const base_element_type lc = modulus_coeffs.back();
    if (!base_.IsOne(lc)) {
      const base_element_type inv_lc = base_.Inv(lc);
      for (base_element_type& c : modulus_coeffs) {
        c = base_.Mul(c, inv_lc);
      }
    }
    modulus_ = std::move(modulus_coeffs);

    // Detect binomial modulus f(x) = x^m - xi (i.e. c_1 = ... = c_{m-1} = 0).
    is_binomial_ = true;
    for (std::size_t i = 1; i < degree_; ++i) {
      if (!base_.IsZero(modulus_[i])) {
        is_binomial_ = false;
        break;
      }
    }
    if (is_binomial_) {
      binomial_xi_ = base_.Neg(modulus_[0]);
    } else {
      binomial_xi_ = base_.AdditiveIdentity();
    }

    base_abs_degree_ = detail::QueryFieldAbsoluteDegree(base_);
    base_big_order_ = detail::QueryFieldBigOrder(base_);

    // Precompute the m x m RelativeFrobenius matrix (columns = (x^q)^j mod f)
    // whenever the base field order q is finite and known.
    if (base_big_order_ > BigInt(0)) {
      PrecomputeFrobeniusMatrix();
    }
  }

  // Convenience factory for binomial extension K = F[x] / (x^degree - xi).
  static FieldTower Binomial(BaseFieldT base_field, std::size_t degree,
                             const base_element_type& xi,
                             std::string var_name = "u") {
    if (degree < 1) {
      throw std::invalid_argument("FieldTower::Binomial requires degree >= 1");
    }
    std::vector<base_element_type> mod(degree + 1,
                                       base_field.AdditiveIdentity());
    mod[0] = base_field.Neg(xi);
    mod[degree] = base_field.MultiplicativeIdentity();
    return FieldTower(std::move(base_field), std::move(mod),
                      std::move(var_name));
  }

  const BaseFieldT& BaseField() const noexcept { return base_; }
  std::size_t RelativeDegree() const noexcept { return degree_; }
  std::size_t AbsoluteDegree() const noexcept {
    return base_abs_degree_ * degree_;
  }
  uint32_t Degree() const noexcept {
    return static_cast<uint32_t>(AbsoluteDegree());
  }

  uint64_t Characteristic() const {
    return detail::QueryFieldCharacteristicU64(base_);
  }

  BigInt BigCharacteristic() const {
    return detail::QueryFieldBigCharacteristic(base_);
  }

  BigInt BigOrder() const {
    if (base_big_order_ <= BigInt(0)) {
      return BigInt(0);
    }
    return base_big_order_.Pow(static_cast<uint64_t>(degree_));
  }

  bool IsBinomial() const noexcept { return is_binomial_; }
  const std::vector<base_element_type>& Modulus() const noexcept {
    return modulus_;
  }
  const std::string& VariableName() const noexcept { return var_name_; }

  element_type AdditiveIdentity() const {
    element_type zero = MakeZeroElement();
    return zero;
  }

  element_type MultiplicativeIdentity() const {
    element_type one = MakeZeroElement();
    one[0] = base_.MultiplicativeIdentity();
    return one;
  }

  // Returns the extension generator x mod f(x).
  element_type Generator() const {
    if (degree_ == 1) {
      return EmbedFromBase(base_.Neg(modulus_[0]));
    }
    element_type gen = MakeZeroElement();
    gen[1] = base_.MultiplicativeIdentity();
    return gen;
  }

  // Constructs a canonical element from arbitrary base-field coefficients
  // (reducing modulo f(x) if length > degree_).
  element_type FromCoeffs(const std::vector<base_element_type>& coeffs) const {
    if (coeffs.size() <= degree_) {
      element_type out = MakeZeroElement();
      for (std::size_t i = 0; i < coeffs.size(); ++i) {
        out[i] = coeffs[i];
      }
      return out;
    }
    std::vector<base_element_type> poly = coeffs;
    ReducePolyInPlace(poly);
    element_type out = MakeZeroElement();
    for (std::size_t i = 0; i < std::min(degree_, poly.size()); ++i) {
      out[i] = poly[i];
    }
    return out;
  }

  element_type EmbedFromBase(const base_element_type& base_elem) const {
    element_type out = MakeZeroElement();
    out[0] = base_elem;
    return out;
  }

  std::optional<base_element_type> ProjectToBase(const element_type& a) const {
    for (std::size_t i = 1; i < degree_; ++i) {
      if (!base_.IsZero(GetCoeff(a, i))) {
        return std::nullopt;
      }
    }
    return GetCoeff(a, 0);
  }

  bool IsZero(const element_type& a) const {
    for (std::size_t i = 0; i < degree_; ++i) {
      if (!base_.IsZero(GetCoeff(a, i))) {
        return false;
      }
    }
    return true;
  }

  bool IsOne(const element_type& a) const {
    if (!base_.IsOne(GetCoeff(a, 0))) {
      return false;
    }
    for (std::size_t i = 1; i < degree_; ++i) {
      if (!base_.IsZero(GetCoeff(a, i))) {
        return false;
      }
    }
    return true;
  }

  bool Equals(const element_type& a, const element_type& b) const {
    for (std::size_t i = 0; i < degree_; ++i) {
      if (!alg::Equals(base_, GetCoeff(a, i), GetCoeff(b, i))) {
        return false;
      }
    }
    return true;
  }

  element_type Add(const element_type& a, const element_type& b) const {
    element_type out = MakeZeroElement();
    for (std::size_t i = 0; i < degree_; ++i) {
      out[i] = base_.Add(GetCoeff(a, i), GetCoeff(b, i));
    }
    return out;
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    element_type out = MakeZeroElement();
    for (std::size_t i = 0; i < degree_; ++i) {
      out[i] = base_.Sub(GetCoeff(a, i), GetCoeff(b, i));
    }
    return out;
  }

  element_type Neg(const element_type& a) const {
    element_type out = MakeZeroElement();
    for (std::size_t i = 0; i < degree_; ++i) {
      out[i] = base_.Neg(GetCoeff(a, i));
    }
    return out;
  }

  element_type Mul(const element_type& a, const element_type& b) const {
    if (is_binomial_) {
      // Fast O(m^2) binomial multiplication modulo x^m - xi without division.
      element_type out = MakeZeroElement();
      for (std::size_t i = 0; i < degree_; ++i) {
        const base_element_type ai = GetCoeff(a, i);
        if (base_.IsZero(ai)) {
          continue;
        }
        for (std::size_t j = 0; j < degree_; ++j) {
          const base_element_type bj = GetCoeff(b, j);
          if (base_.IsZero(bj)) {
            continue;
          }
          const base_element_type prod = base_.Mul(ai, bj);
          if (i + j < degree_) {
            out[i + j] = base_.Add(out[i + j], prod);
          } else {
            const std::size_t idx = i + j - degree_;
            out[idx] = base_.Add(out[idx], base_.Mul(prod, binomial_xi_));
          }
        }
      }
      return out;
    }

    // General monic modulus f(x) multiplication + top-down reduction.
    std::vector<base_element_type> prod(2 * degree_ - 1,
                                        base_.AdditiveIdentity());
    for (std::size_t i = 0; i < degree_; ++i) {
      const base_element_type ai = GetCoeff(a, i);
      if (base_.IsZero(ai)) {
        continue;
      }
      for (std::size_t j = 0; j < degree_; ++j) {
        const base_element_type bj = GetCoeff(b, j);
        if (base_.IsZero(bj)) {
          continue;
        }
        prod[i + j] = base_.Add(prod[i + j], base_.Mul(ai, bj));
      }
    }
    ReducePolyInPlace(prod);
    element_type out = MakeZeroElement();
    for (std::size_t i = 0; i < degree_; ++i) {
      out[i] = prod[i];
    }
    return out;
  }

  // Multiplicative inverse via polynomial Extended Euclidean Algorithm over
  // BaseFieldT.
  element_type Inv(const element_type& a) const {
    std::vector<base_element_type> r1 = ToTrimmedPoly(a);
    if (r1.empty()) {
      throw std::domain_error("FieldTower::Inv: inverse of zero is undefined");
    }
    if (r1.size() == 1) {
      return EmbedFromBase(base_.Inv(r1[0]));
    }

    std::vector<base_element_type> r0 = modulus_;
    std::vector<base_element_type> t0;
    std::vector<base_element_type> t1 = {base_.MultiplicativeIdentity()};

    while (!r1.empty()) {
      auto [q, rem] = PolyDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(rem);

      std::vector<base_element_type> qt1 = PolyMul(q, t1);
      std::vector<base_element_type> t_next = PolySub(t0, qt1);
      t0 = std::move(t1);
      t1 = std::move(t_next);
    }

    if (r0.size() != 1) {
      throw std::domain_error(
          "FieldTower::Inv: element is not invertible (modulus may be "
          "reducible)");
    }

    const base_element_type inv_c = base_.Inv(r0[0]);
    for (base_element_type& coeff : t0) {
      coeff = base_.Mul(coeff, inv_c);
    }
    ReducePolyInPlace(t0);

    element_type out = MakeZeroElement();
    for (std::size_t i = 0; i < std::min(degree_, t0.size()); ++i) {
      out[i] = t0[i];
    }
    return out;
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return Mul(a, Inv(b));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    element_type result = MultiplicativeIdentity();
    element_type base = a;
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        result = Mul(result, base);
      }
      exp >>= 1;
      if (exp > 0) {
        base = Mul(base, base);
      }
    }
    return result;
  }

  element_type PowBig(const element_type& a, const BigInt& exp) const {
    if (exp.IsZero()) {
      return MultiplicativeIdentity();
    }
    if (exp.IsNegative()) {
      return PowBig(Inv(a), -exp);
    }
    element_type result = MultiplicativeIdentity();
    element_type base = a;
    const std::size_t bits = exp.BitLength();
    for (std::size_t i = 0; i < bits; ++i) {
      if (exp.GetBit(i)) {
        result = Mul(result, base);
      }
      if (i + 1 < bits) {
        base = Mul(base, base);
      }
    }
    return result;
  }

  // Relative Frobenius automorphism a -> a^(q^steps), where q = |BaseField|.
  element_type RelativeFrobenius(const element_type& a,
                                 uint64_t steps = 1) const {
    if (base_big_order_ <= BigInt(0)) {
      throw std::domain_error(
          "FieldTower::RelativeFrobenius requires a finite base field");
    }
    const std::size_t s = static_cast<std::size_t>(steps % degree_);
    if (s == 0) {
      return a;
    }
    element_type cur = a;
    for (std::size_t step = 0; step < s; ++step) {
      cur = ApplyFrobeniusOnce(cur);
    }
    return cur;
  }

  // Relative field trace Tr_{K/F}(a) = sum_{i=0}^{m-1} a^(q^i) in BaseField.
  base_element_type RelativeTrace(const element_type& a) const {
    if (base_big_order_ <= BigInt(0)) {
      throw std::domain_error(
          "FieldTower::RelativeTrace requires a finite base field");
    }
    element_type acc = a;
    element_type cur = a;
    for (std::size_t i = 1; i < degree_; ++i) {
      cur = ApplyFrobeniusOnce(cur);
      acc = Add(acc, cur);
    }
    return GetCoeff(acc, 0);
  }

  // Relative field norm N_{K/F}(a) = prod_{i=0}^{m-1} a^(q^i) in BaseField.
  base_element_type RelativeNorm(const element_type& a) const {
    if (base_big_order_ <= BigInt(0)) {
      throw std::domain_error(
          "FieldTower::RelativeNorm requires a finite base field");
    }
    element_type acc = a;
    element_type cur = a;
    for (std::size_t i = 1; i < degree_; ++i) {
      cur = ApplyFrobeniusOnce(cur);
      acc = Mul(acc, cur);
    }
    return GetCoeff(acc, 0);
  }

  std::string ToString(const element_type& a) const {
    std::ostringstream oss;
    bool first = true;
    for (std::size_t idx = degree_; idx-- > 0;) {
      const base_element_type c = GetCoeff(a, idx);
      if (base_.IsZero(c)) {
        continue;
      }
      if (!first) {
        oss << " + ";
      }
      first = false;
      const std::string cs = base_.ToString(c);
      const bool needs_parens = (cs.find('+') != std::string::npos ||
                                 cs.find(' ') != std::string::npos);
      if (idx == 0) {
        oss << cs;
      } else {
        if (!base_.IsOne(c)) {
          if (needs_parens) {
            oss << "(" << cs << ")*";
          } else {
            oss << cs << "*";
          }
        }
        oss << var_name_;
        if (idx > 1) {
          oss << "^" << idx;
        }
      }
    }
    if (first) {
      return "0";
    }
    return oss.str();
  }

 private:
  element_type MakeZeroElement() const {
    element_type elem;
    const base_element_type z = base_.AdditiveIdentity();
    if constexpr (StaticDegree > 0) {
      elem.coeffs.fill(z);
    } else {
      elem.coeffs.assign(degree_, z);
    }
    return elem;
  }

  base_element_type GetCoeff(const element_type& a, std::size_t i) const {
    if constexpr (StaticDegree > 0) {
      return a.coeffs[i];
    } else {
      return (i < a.coeffs.size()) ? a.coeffs[i] : base_.AdditiveIdentity();
    }
  }

  void TrimPoly(std::vector<base_element_type>& p) const {
    while (!p.empty() && base_.IsZero(p.back())) {
      p.pop_back();
    }
  }

  std::vector<base_element_type> ToTrimmedPoly(const element_type& a) const {
    std::vector<base_element_type> p(degree_);
    for (std::size_t i = 0; i < degree_; ++i) {
      p[i] = GetCoeff(a, i);
    }
    TrimPoly(p);
    return p;
  }

  void ReducePolyInPlace(std::vector<base_element_type>& poly) const {
    if (poly.size() <= degree_) {
      return;
    }
    for (std::size_t i = poly.size(); i-- > degree_;) {
      const base_element_type lead = poly[i];
      if (base_.IsZero(lead)) {
        continue;
      }
      const std::size_t shift = i - degree_;
      poly[i] = base_.AdditiveIdentity();
      if (is_binomial_) {
        poly[shift] = base_.Add(poly[shift], base_.Mul(lead, binomial_xi_));
      } else {
        for (std::size_t j = 0; j < degree_; ++j) {
          if (!base_.IsZero(modulus_[j])) {
            poly[shift + j] =
                base_.Sub(poly[shift + j], base_.Mul(lead, modulus_[j]));
          }
        }
      }
    }
    poly.resize(degree_, base_.AdditiveIdentity());
  }

  std::vector<base_element_type> PolySub(
      const std::vector<base_element_type>& a,
      const std::vector<base_element_type>& b) const {
    const std::size_t n = std::max(a.size(), b.size());
    std::vector<base_element_type> res(n, base_.AdditiveIdentity());
    for (std::size_t i = 0; i < n; ++i) {
      const base_element_type av =
          (i < a.size()) ? a[i] : base_.AdditiveIdentity();
      const base_element_type bv =
          (i < b.size()) ? b[i] : base_.AdditiveIdentity();
      res[i] = base_.Sub(av, bv);
    }
    TrimPoly(res);
    return res;
  }

  std::vector<base_element_type> PolyMul(
      const std::vector<base_element_type>& a,
      const std::vector<base_element_type>& b) const {
    if (a.empty() || b.empty()) {
      return {};
    }
    std::vector<base_element_type> res(a.size() + b.size() - 1,
                                       base_.AdditiveIdentity());
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (base_.IsZero(a[i])) {
        continue;
      }
      for (std::size_t j = 0; j < b.size(); ++j) {
        if (base_.IsZero(b[j])) {
          continue;
        }
        res[i + j] = base_.Add(res[i + j], base_.Mul(a[i], b[j]));
      }
    }
    TrimPoly(res);
    return res;
  }

  std::pair<std::vector<base_element_type>, std::vector<base_element_type>>
  PolyDivRem(const std::vector<base_element_type>& dividend,
             const std::vector<base_element_type>& divisor) const {
    if (divisor.empty()) {
      throw std::domain_error("FieldTower: polynomial division by zero");
    }
    if (dividend.size() < divisor.size()) {
      return {{}, dividend};
    }
    std::vector<base_element_type> rem = dividend;
    const std::size_t d_deg = divisor.size() - 1;
    const base_element_type inv_lc = base_.Inv(divisor[d_deg]);
    std::vector<base_element_type> quot(dividend.size() - d_deg,
                                        base_.AdditiveIdentity());

    for (std::size_t i = dividend.size(); i-- > d_deg;) {
      const base_element_type c = rem[i];
      if (base_.IsZero(c)) {
        continue;
      }
      const base_element_type q_coeff = base_.Mul(c, inv_lc);
      const std::size_t shift = i - d_deg;
      quot[shift] = q_coeff;
      for (std::size_t j = 0; j <= d_deg; ++j) {
        if (!base_.IsZero(divisor[j])) {
          rem[shift + j] =
              base_.Sub(rem[shift + j], base_.Mul(q_coeff, divisor[j]));
        }
      }
    }
    TrimPoly(quot);
    TrimPoly(rem);
    return {std::move(quot), std::move(rem)};
  }

  void PrecomputeFrobeniusMatrix() {
    frob_basis_.resize(degree_);
    frob_basis_[0] = MultiplicativeIdentity();
    if (degree_ == 1) {
      return;
    }
    const element_type x_q = PowBig(Generator(), base_big_order_);
    frob_basis_[1] = x_q;
    for (std::size_t j = 2; j < degree_; ++j) {
      frob_basis_[j] = Mul(frob_basis_[j - 1], x_q);
    }
  }

  element_type ApplyFrobeniusOnce(const element_type& a) const {
    element_type out = MakeZeroElement();
    for (std::size_t j = 0; j < degree_; ++j) {
      const base_element_type aj = GetCoeff(a, j);
      if (base_.IsZero(aj)) {
        continue;
      }
      const element_type& col = frob_basis_[j];
      for (std::size_t i = 0; i < degree_; ++i) {
        const base_element_type cij = GetCoeff(col, i);
        if (!base_.IsZero(cij)) {
          out[i] = base_.Add(out[i], base_.Mul(aj, cij));
        }
      }
    }
    return out;
  }

  BaseFieldT base_;
  std::string var_name_;
  std::size_t degree_ = 0;
  std::size_t base_abs_degree_ = 1;
  BigInt base_big_order_;
  std::vector<base_element_type> modulus_;
  bool is_binomial_ = false;
  base_element_type binomial_xi_{};
  std::vector<element_type> frob_basis_;
};

template <Field BaseFieldT, std::size_t StaticDegree>
struct AlgebraTraits<FieldTower<BaseFieldT, StaticDegree>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = true;
};

// ============================================================================
// Constant-Time, Zero-LUT Composite Tower GF(((2^2)^2)^2) (Canright 2005)
// ============================================================================
//
// Normal-basis tower:
//   GF(2^2)       = GF(2)[W]/(W^2 + W + 1),           basis {W^2, W},  N = W^2
//   GF((2^2)^2)   = GF(2^2)[Z]/(Z^2 + Z + N),         basis {Z^4, Z},  nu =
//   N*Z^4 GF(((2^2)^2)^2)= GF(2^4)[Y]/(Y^2 + Y + nu),       basis {Y^16, Y}
//
// All operations are 100% branchless bitwise arithmetic with zero S-box/log
// lookup tables, matching AES GF(2^8) (x^8 + x^4 + x^3 + x + 1) and FIPS-197.
class CompositeGF256Tower {
 public:
  // Converts from AES standard polynomial basis (mod 0x11B) to tower basis.
  static constexpr uint8_t ToTower(uint8_t aes_byte) noexcept {
    return MatVecMul8(kToTowerCols, aes_byte);
  }

  // Converts from tower basis back to AES standard polynomial basis (mod
  // 0x11B).
  static constexpr uint8_t FromTower(uint8_t tower_byte) noexcept {
    return MatVecMul8(kFromTowerCols, tower_byte);
  }

  // Constant-time multiplication directly in the composite tower basis.
  static constexpr uint8_t MulTower(uint8_t a, uint8_t b) noexcept {
    return Gf256Mul(a, b);
  }

  // Constant-time multiplicative inverse in the composite tower basis
  // (InvTower(0) == 0 following the AES / Canright convention).
  static constexpr uint8_t InvTower(uint8_t a) noexcept { return Gf256Inv(a); }

  // Constant-time multiplication in AES GF(2^8) (mod x^8 + x^4 + x^3 + x + 1).
  static constexpr uint8_t Mul(uint8_t a, uint8_t b) noexcept {
    return FromTower(Gf256Mul(ToTower(a), ToTower(b)));
  }

  // Constant-time multiplicative inverse in AES GF(2^8) (Inv(0) == 0).
  static constexpr uint8_t Inv(uint8_t a) noexcept {
    return FromTower(Gf256Inv(ToTower(a)));
  }

  // Constant-time, zero-LUT FIPS-197 AES SubBytes S-box.
  static constexpr uint8_t AesSbox(uint8_t x) noexcept {
    return AesAffine(Inv(x));
  }

  // Constant-time, zero-LUT FIPS-197 AES InvSubBytes inverse S-box.
  static constexpr uint8_t AesInvSbox(uint8_t x) noexcept {
    return Inv(AesInvAffine(x));
  }

 private:
  // 8x8 GF(2) isomorphism matrices between AES GF(2^8) mod 0x11B and
  // GF(((2^2)^2)^2) normal basis. Stored as 8 column bytes (64 bits total).
  static constexpr std::array<uint8_t, 8> kToTowerCols = {
      0xffU, 0x06U, 0x56U, 0xb8U, 0xb2U, 0xc8U, 0x3aU, 0xc1U};

  static constexpr std::array<uint8_t, 8> kFromTowerCols = {
      0x7bU, 0xc3U, 0xc1U, 0xdbU, 0x76U, 0x2eU, 0x70U, 0x8bU};

  static constexpr uint8_t MatVecMul8(const std::array<uint8_t, 8>& cols,
                                      uint8_t x) noexcept {
    uint8_t res = 0;
    for (std::size_t i = 0; i < 8; ++i) {
      const uint8_t mask =
          static_cast<uint8_t>(0U - static_cast<unsigned>((x >> i) & 1U));
      res = static_cast<uint8_t>(res ^ (cols[i] & mask));
    }
    return res;
  }

  static constexpr uint8_t RotL8(uint8_t x, unsigned r) noexcept {
    return static_cast<uint8_t>((x << r) | (x >> (8U - r)));
  }

  // FIPS-197 affine transform: b_i = x_i ^ x_{i+4} ^ x_{i+5} ^ x_{i+6} ^
  // x_{i+7} ^ c_i with c = 0x63.
  static constexpr uint8_t AesAffine(uint8_t x) noexcept {
    return static_cast<uint8_t>(x ^ RotL8(x, 1) ^ RotL8(x, 2) ^ RotL8(x, 3) ^
                                RotL8(x, 4) ^ 0x63U);
  }

  // Inverse FIPS-197 affine transform.
  static constexpr uint8_t AesInvAffine(uint8_t x) noexcept {
    return static_cast<uint8_t>(RotL8(x, 1) ^ RotL8(x, 3) ^ RotL8(x, 6) ^
                                0x05U);
  }

  // ---- Level 1: GF(2^2), normal basis {W^2, W} ----

  static constexpr uint8_t Gf4Mul(uint8_t a, uint8_t b) noexcept {
    const uint8_t a1 = (a >> 1) & 1U;
    const uint8_t a0 = a & 1U;
    const uint8_t b1 = (b >> 1) & 1U;
    const uint8_t b0 = b & 1U;
    const uint8_t t = static_cast<uint8_t>((a1 ^ a0) & (b1 ^ b0));
    const uint8_t c1 = static_cast<uint8_t>((a1 & b1) ^ t);
    const uint8_t c0 = static_cast<uint8_t>((a0 & b0) ^ t);
    return static_cast<uint8_t>((c1 << 1) | c0);
  }

  // Squaring and inversion coincide in GF(4): swap bits (a1, a0) -> (a0, a1).
  static constexpr uint8_t Gf4Sq(uint8_t a) noexcept {
    return static_cast<uint8_t>(((a & 1U) << 1) | ((a >> 1) & 1U));
  }

  // Multiply by N = W^2 in GF(4): (a1, a0) -> (a0, a1 ^ a0).
  static constexpr uint8_t Gf4ScaleN(uint8_t a) noexcept {
    const uint8_t a1 = (a >> 1) & 1U;
    const uint8_t a0 = a & 1U;
    return static_cast<uint8_t>((a0 << 1) | (a1 ^ a0));
  }

  // Multiply by N^2 = W in GF(4): (a1, a0) -> (a1 ^ a0, a1).
  static constexpr uint8_t Gf4ScaleN2(uint8_t a) noexcept {
    const uint8_t a1 = (a >> 1) & 1U;
    const uint8_t a0 = a & 1U;
    return static_cast<uint8_t>(((a1 ^ a0) << 1) | a1);
  }

  // ---- Level 2: GF((2^2)^2), normal basis {Z^4, Z} ----

  static constexpr uint8_t Gf16Mul(uint8_t a, uint8_t b) noexcept {
    const uint8_t a1 = (a >> 2) & 0x3U;
    const uint8_t a0 = a & 0x3U;
    const uint8_t b1 = (b >> 2) & 0x3U;
    const uint8_t b0 = b & 0x3U;
    const uint8_t t = Gf4ScaleN(
        Gf4Mul(static_cast<uint8_t>(a1 ^ a0), static_cast<uint8_t>(b1 ^ b0)));
    const uint8_t c1 = static_cast<uint8_t>(Gf4Mul(a1, b1) ^ t);
    const uint8_t c0 = static_cast<uint8_t>(Gf4Mul(a0, b0) ^ t);
    return static_cast<uint8_t>((c1 << 2) | c0);
  }

  static constexpr uint8_t Gf16Sq(uint8_t a) noexcept {
    const uint8_t a1 = (a >> 2) & 0x3U;
    const uint8_t a0 = a & 0x3U;
    const uint8_t t = Gf4ScaleN(Gf4Sq(static_cast<uint8_t>(a1 ^ a0)));
    const uint8_t c1 = static_cast<uint8_t>(Gf4Sq(a1) ^ t);
    const uint8_t c0 = static_cast<uint8_t>(Gf4Sq(a0) ^ t);
    return static_cast<uint8_t>((c1 << 2) | c0);
  }

  // Multiply by nu = N * Z^4 in GF(16): (a1, a0) -> (a1 ^ N^2*a0, N^2*(a1^a0)).
  static constexpr uint8_t Gf16ScaleNu(uint8_t a) noexcept {
    const uint8_t a1 = (a >> 2) & 0x3U;
    const uint8_t a0 = a & 0x3U;
    const uint8_t c1 = static_cast<uint8_t>(a1 ^ Gf4ScaleN2(a0));
    const uint8_t c0 = Gf4ScaleN2(static_cast<uint8_t>(a1 ^ a0));
    return static_cast<uint8_t>((c1 << 2) | c0);
  }

  static constexpr uint8_t Gf16Inv(uint8_t a) noexcept {
    const uint8_t a1 = (a >> 2) & 0x3U;
    const uint8_t a0 = a & 0x3U;
    const uint8_t norm = static_cast<uint8_t>(
        Gf4Mul(a1, a0) ^ Gf4ScaleN(Gf4Sq(static_cast<uint8_t>(a1 ^ a0))));
    const uint8_t inv_norm = Gf4Sq(norm);
    return static_cast<uint8_t>((Gf4Mul(inv_norm, a0) << 2) |
                                Gf4Mul(inv_norm, a1));
  }

  // ---- Level 3: GF(((2^2)^2)^2), normal basis {Y^16, Y} ----

  static constexpr uint8_t Gf256Mul(uint8_t a, uint8_t b) noexcept {
    const uint8_t a1 = (a >> 4) & 0xFU;
    const uint8_t a0 = a & 0xFU;
    const uint8_t b1 = (b >> 4) & 0xFU;
    const uint8_t b0 = b & 0xFU;
    const uint8_t t = Gf16ScaleNu(
        Gf16Mul(static_cast<uint8_t>(a1 ^ a0), static_cast<uint8_t>(b1 ^ b0)));
    const uint8_t c1 = static_cast<uint8_t>(Gf16Mul(a1, b1) ^ t);
    const uint8_t c0 = static_cast<uint8_t>(Gf16Mul(a0, b0) ^ t);
    return static_cast<uint8_t>((c1 << 4) | c0);
  }

  static constexpr uint8_t Gf256Inv(uint8_t a) noexcept {
    const uint8_t a1 = (a >> 4) & 0xFU;
    const uint8_t a0 = a & 0xFU;
    const uint8_t norm = static_cast<uint8_t>(
        Gf16Mul(a1, a0) ^ Gf16ScaleNu(Gf16Sq(static_cast<uint8_t>(a1 ^ a0))));
    const uint8_t inv_norm = Gf16Inv(norm);
    return static_cast<uint8_t>((Gf16Mul(inv_norm, a0) << 4) |
                                Gf16Mul(inv_norm, a1));
  }
};

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_TOWER_HPP_
