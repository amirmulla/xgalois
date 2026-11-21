#ifndef XGALOIS_FIELD_GF_EXTENSION_HPP_
#define XGALOIS_FIELD_GF_EXTENSION_HPP_

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

// Compact, default-constructible polynomial element for general extension
// fields GF(p^m). Stores trimmed base-field coefficients in ascending degree
// order (empty vector represents the zero polynomial with Degree() == -1).
template <UnsignedFieldWord BaseFieldElementType = uint32_t>
class ExtensionElement {
 public:
  using coeff_type = BaseFieldElementType;
  using value_type = BaseFieldElementType;

  ExtensionElement() = default;

  explicit ExtensionElement(std::vector<BaseFieldElementType> coeffs)
      : coeffs_(std::move(coeffs)) {
    Trim();
  }

  ExtensionElement(std::initializer_list<BaseFieldElementType> coeffs)
      : coeffs_(coeffs) {
    Trim();
  }

  int Degree() const noexcept { return static_cast<int>(coeffs_.size()) - 1; }

  size_t Size() const noexcept { return coeffs_.size(); }

  bool IsZero() const noexcept { return coeffs_.empty(); }

  BaseFieldElementType operator[](size_t i) const noexcept {
    return i < coeffs_.size() ? coeffs_[i] : BaseFieldElementType{0};
  }

  const std::vector<BaseFieldElementType> &Coeffs() const noexcept {
    return coeffs_;
  }

  void Trim() noexcept {
    while (!coeffs_.empty() && coeffs_.back() == BaseFieldElementType{0}) {
      coeffs_.pop_back();
    }
  }

  bool operator==(const ExtensionElement &other) const noexcept {
    return coeffs_ == other.coeffs_;
  }

  bool operator!=(const ExtensionElement &other) const noexcept {
    return coeffs_ != other.coeffs_;
  }

  // Lexicographic order by degree then coefficients from highest to lowest
  // degree, matching base-p integer ordering.
  bool operator<(const ExtensionElement &other) const noexcept {
    if (coeffs_.size() != other.coeffs_.size()) {
      return coeffs_.size() < other.coeffs_.size();
    }
    for (size_t i = coeffs_.size(); i > 0; --i) {
      const size_t idx = i - 1;
      if (coeffs_[idx] != other.coeffs_[idx]) {
        return coeffs_[idx] < other.coeffs_[idx];
      }
    }
    return false;
  }

  bool operator<=(const ExtensionElement &other) const noexcept {
    return !(other < *this);
  }

  bool operator>(const ExtensionElement &other) const noexcept {
    return other < *this;
  }

  bool operator>=(const ExtensionElement &other) const noexcept {
    return !(*this < other);
  }

 private:
  std::vector<BaseFieldElementType> coeffs_;
};

}  // namespace xg

namespace std {

template <xg::UnsignedFieldWord T>
struct hash<xg::ExtensionElement<T>> {
  size_t operator()(const xg::ExtensionElement<T> &elem) const noexcept {
    size_t h = std::hash<size_t>{}(elem.Size());
    for (const T &c : elem.Coeffs()) {
      const size_t hc = std::hash<T>{}(c);
      h ^= hc + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return h;
  }
};

}  // namespace std

namespace xg {

template <UnsignedFieldWord BaseFieldElementType = uint32_t>
class GaloisFieldExtension
    : public GaloisFieldBase<ExtensionElement<BaseFieldElementType>> {
 public:
  using BaseField = GaloisFieldPrime<BaseFieldElementType>;
  using ElementType = ExtensionElement<BaseFieldElementType>;
  using PolynomialType = ElementType;

  explicit GaloisFieldExtension(std::pair<uint64_t, uint64_t> prime_exp,
                                const std::string &irreducible_poly = "",
                                const std::string &rep = "poly",
                                const std::string &poly_var = "α",
                                bool check_irreducible = false,
                                bool prime_testing = true,
                                const std::string &generator_name = "g")
      : base_field_(ValidatePrime(prime_exp.first, prime_testing), "int",
                    prime_testing),
        degree_(ValidateExponent(prime_exp.second)),
        representation_(utils::ConvertRepresentation(rep)),
        variable_name_(poly_var.empty() ? "α" : poly_var),
        generator_name_(generator_name.empty() ? "g" : generator_name) {
    order_ = utils::SafeIntegerPower(prime_exp.first, prime_exp.second);

    auto [modulus_poly, from_conway] = CreateIrreduciblePolynomial(
        base_field_, irreducible_poly, degree_, variable_name_);
    irreducible_poly_ = std::move(modulus_poly);
    used_conway_poly_ = from_conway;

    if (irreducible_poly_.Degree() != static_cast<int>(degree_)) {
      throw std::invalid_argument("Irreducible polynomial degree (" +
                                  std::to_string(irreducible_poly_.Degree()) +
                                  ") does not match extension degree (" +
                                  std::to_string(degree_) + ")");
    }

    if (check_irreducible && !IsIrreducibleOverBase(irreducible_poly_)) {
      throw std::invalid_argument("Provided polynomial is not irreducible");
    }
  }

  explicit GaloisFieldExtension(uint64_t order,
                                const std::string &irreducible_poly = "",
                                const std::string &rep = "poly",
                                const std::string &poly_var = "α",
                                bool check_irreducible = false,
                                bool prime_testing = true,
                                const std::string &generator_name = "g")
      : GaloisFieldExtension(DecomposeOrderToPrimeExp(order), irreducible_poly,
                             rep, poly_var, check_irreducible, prime_testing,
                             generator_name) {}

  GaloisFieldExtension(const GaloisFieldExtension &other)
      : base_field_(static_cast<BaseFieldElementType>(
                        other.base_field_.Characteristic()),
                    "int", false),
        degree_(other.degree_),
        irreducible_poly_(other.irreducible_poly_),
        order_(other.order_),
        representation_(other.representation_.load(std::memory_order_relaxed)),
        variable_name_(other.variable_name_),
        generator_name_(other.generator_name_),
        used_conway_poly_(other.used_conway_poly_),
        generator_cache_(other.generator_cache_) {}

  GaloisFieldExtension &operator=(const GaloisFieldExtension &) = delete;

  uint64_t Characteristic() const override {
    return base_field_.Characteristic();
  }

  uint64_t Order() const override { return order_; }

  uint32_t Degree() const override { return degree_; }

  const ElementType &Modulus() const noexcept { return irreducible_poly_; }

  const BaseField &Base() const noexcept { return base_field_; }

  const std::string &VariableName() const noexcept { return variable_name_; }

  const std::string &GeneratorName() const noexcept { return generator_name_; }

  ElementType Add(const ElementType &a, const ElementType &b) const override {
    if (!IsCanonical(a) || !IsCanonical(b)) {
      return PolyAdd(SetElementValue(a), SetElementValue(b));
    }
    return PolyAdd(a, b);
  }

  ElementType Sub(const ElementType &a, const ElementType &b) const override {
    if (!IsCanonical(a) || !IsCanonical(b)) {
      return PolySub(SetElementValue(a), SetElementValue(b));
    }
    return PolySub(a, b);
  }

  ElementType Mul(const ElementType &a, const ElementType &b) const override {
    if (!IsCanonical(a) || !IsCanonical(b)) {
      return PolyMod(PolyMulRaw(SetElementValue(a), SetElementValue(b)),
                     irreducible_poly_);
    }
    if (a.IsZero() || b.IsZero()) {
      return AdditiveIdentity();
    }
    return PolyMod(PolyMulRaw(a, b), irreducible_poly_);
  }

  ElementType Div(const ElementType &a, const ElementType &b) const override {
    ElementType nb = IsCanonical(b) ? b : SetElementValue(b);
    if (nb.IsZero()) {
      throw std::domain_error("Division by zero in GaloisFieldExtension");
    }
    return Mul(a, Inv(nb));
  }

  ElementType Neg(const ElementType &a) const override {
    if (!IsCanonical(a)) {
      return PolyNeg(SetElementValue(a));
    }
    return PolyNeg(a);
  }

  ElementType Inv(const ElementType &a) const override {
    ElementType na = IsCanonical(a) ? a : SetElementValue(a);
    if (na.IsZero()) {
      throw std::domain_error(
          "Multiplicative inverse of zero is undefined in "
          "GaloisFieldExtension");
    }

    auto [gcd, s, t] = PolyExtendedGcd(na, irreducible_poly_);
    (void)t;

    if (gcd.Degree() != 0 || gcd[0] != BaseFieldElementType{1}) {
      throw std::runtime_error(
          "ExtendedGcd did not return 1 for inverse calculation. "
          "Ensure the modulus is irreducible and input is non-zero.");
    }

    return PolyMod(s, irreducible_poly_);
  }

  ElementType Pow(const ElementType &base, uint64_t exp) const override {
    if (exp == 0) {
      return MultiplicativeIdentity();
    }
    ElementType current_base = IsCanonical(base) ? base : SetElementValue(base);
    if (current_base.IsZero()) {
      return AdditiveIdentity();
    }
    const uint64_t group_order = order_ - 1;
    uint64_t current_exp = exp % group_order;
    if (current_exp == 0) {
      return MultiplicativeIdentity();
    }
    if (current_exp == 1) {
      return current_base;
    }

    ElementType result = MultiplicativeIdentity();
    while (current_exp > 0) {
      if (current_exp & 1ULL) {
        result = PolyMod(PolyMulRaw(result, current_base), irreducible_poly_);
      }
      current_exp >>= 1;
      if (current_exp > 0) {
        current_base =
            PolyMod(PolyMulRaw(current_base, current_base), irreducible_poly_);
      }
    }
    return result;
  }

  ElementType Frobenius(const ElementType &a) const override {
    return Pow(a, base_field_.Characteristic());
  }

  ElementType Trace(const ElementType &a) const override {
    ElementType cur = IsCanonical(a) ? a : SetElementValue(a);
    ElementType acc = cur;
    for (uint32_t i = 1; i < degree_; ++i) {
      cur = Frobenius(cur);
      acc = Add(acc, cur);
    }
    return acc;
  }

  ElementType Norm(const ElementType &a) const override {
    ElementType cur = IsCanonical(a) ? a : SetElementValue(a);
    ElementType acc = cur;
    for (uint32_t i = 1; i < degree_; ++i) {
      cur = Frobenius(cur);
      acc = Mul(acc, cur);
    }
    return acc;
  }

  // Returns the monic minimal polynomial m_a(X) in F_p[X] of a in F_{p^m}.
  ElementType MinimalPolynomial(const ElementType &a) const {
    ElementType start = IsCanonical(a) ? a : SetElementValue(a);
    std::vector<ElementType> conjugates;
    conjugates.reserve(degree_);
    ElementType cur = start;
    do {
      conjugates.push_back(cur);
      cur = Frobenius(cur);
    } while (cur != start && conjugates.size() < degree_);

    std::vector<ElementType> poly;
    poly.reserve(conjugates.size() + 1);
    poly.push_back(MultiplicativeIdentity());

    for (const ElementType &conj : conjugates) {
      std::vector<ElementType> next(poly.size() + 1, AdditiveIdentity());
      const ElementType neg_conj = Neg(conj);
      for (size_t i = 0; i < poly.size(); ++i) {
        next[i] = Add(next[i], Mul(poly[i], neg_conj));
        next[i + 1] = Add(next[i + 1], poly[i]);
      }
      poly = std::move(next);
    }

    std::vector<BaseFieldElementType> coeffs(poly.size(),
                                             BaseFieldElementType{0});
    for (size_t i = 0; i < poly.size(); ++i) {
      coeffs[i] = poly[i][0];
    }
    return ElementType(std::move(coeffs));
  }

  ElementType Sqrt(const ElementType &a) const override {
    ElementType na = IsCanonical(a) ? a : SetElementValue(a);
    if (na.IsZero()) {
      return AdditiveIdentity();
    }
    if (IsOne(na)) {
      return MultiplicativeIdentity();
    }

    // Characteristic 2: Frobenius inverse a^(2^(m-1)).
    if (base_field_.Characteristic() == 2) {
      ElementType root = na;
      for (uint32_t i = 1; i < degree_; ++i) {
        root = Mul(root, root);
      }
      return root;
    }

    // Odd characteristic: Euler's criterion + Tonelli-Shanks in F_{p^m}.
    const uint64_t group_order = order_ - 1;
    if (!IsOne(Pow(na, group_order >> 1))) {
      throw std::domain_error(
          "Element is not a quadratic residue in GaloisFieldExtension");
    }

    uint64_t q_odd = group_order;
    uint32_t s = 0;
    while ((q_odd & 1ULL) == 0) {
      q_odd >>= 1;
      ++s;
    }

    ElementType r;
    if (s == 1) {
      r = Pow(na, (q_odd >> 1) + 1);
    } else {
      // Every multiplicative generator is guaranteed to be a quadratic
      // non-residue in F_q^* for odd q.
      const ElementType z = MultiplicativeGenerator();
      uint32_t m_exp = s;
      ElementType c = Pow(z, q_odd);
      ElementType t = Pow(na, q_odd);
      r = Pow(na, (q_odd >> 1) + 1);

      while (!IsOne(t)) {
        uint32_t i = 0;
        ElementType temp = t;
        while (!IsOne(temp)) {
          temp = Mul(temp, temp);
          ++i;
          if (i >= m_exp) {
            throw std::runtime_error(
                "Tonelli-Shanks failed to converge in GaloisFieldExtension");
          }
        }
        ElementType b = c;
        for (uint32_t j = 0; j < m_exp - i - 1; ++j) {
          b = Mul(b, b);
        }
        m_exp = i;
        c = Mul(b, b);
        t = Mul(t, c);
        r = Mul(r, b);
      }
    }

    const ElementType neg_r = Neg(r);
    return (neg_r < r) ? neg_r : r;
  }

  uint64_t Log(const ElementType &a) const override {
    ElementType na = IsCanonical(a) ? a : SetElementValue(a);
    if (na.IsZero()) {
      throw std::domain_error(
          "Logarithm of zero is undefined in GaloisFieldExtension");
    }
    return Log(na, MultiplicativeGenerator());
  }

  // Baby-Step Giant-Step discrete logarithm in F_{p^m}^*.
  uint64_t Log(const ElementType &a,
               const ElementType &generator) const override {
    ElementType na = IsCanonical(a) ? a : SetElementValue(a);
    if (na.IsZero()) {
      throw std::domain_error(
          "Logarithm of zero is undefined in GaloisFieldExtension");
    }
    ElementType ng =
        IsCanonical(generator) ? generator : SetElementValue(generator);
    if (ng.IsZero()) {
      throw std::invalid_argument(
          "Logarithm base cannot be zero in GaloisFieldExtension");
    }
    if (IsOne(na)) {
      return 0;
    }

    const uint64_t group_order = order_ - 1;
    const uint64_t step = utils::ISqrt(group_order) + 1;

    std::unordered_map<ElementType, uint64_t> baby_steps;
    baby_steps.reserve(
        static_cast<size_t>(std::min<uint64_t>(step, 1ULL << 20)));

    ElementType cur = MultiplicativeIdentity();
    for (uint64_t j = 0; j < step; ++j) {
      baby_steps.try_emplace(cur, j);
      cur = Mul(cur, ng);
    }

    const ElementType giant_factor = Inv(Pow(ng, step));
    ElementType gamma = na;
    for (uint64_t i = 0; i <= step; ++i) {
      auto it = baby_steps.find(gamma);
      if (it != baby_steps.end()) {
        return i * step + it->second;
      }
      gamma = Mul(gamma, giant_factor);
    }

    throw std::domain_error(
        "Element is not in the multiplicative subgroup generated by the given "
        "base");
  }

  ElementType Random() const override {
    static thread_local std::mt19937_64 rng(std::random_device{}());
    return Random(rng);
  }

  ElementType Random(std::mt19937_64 &rng) const override {
    std::uniform_int_distribution<uint64_t> dist(
        0, base_field_.Characteristic() - 1);
    std::vector<BaseFieldElementType> coeffs(degree_);
    for (uint32_t i = 0; i < degree_; ++i) {
      coeffs[i] = static_cast<BaseFieldElementType>(dist(rng));
    }
    return ElementType(std::move(coeffs));
  }

  ElementType MultiplicativeIdentity() const override {
    return ElementType({BaseFieldElementType{1}});
  }

  ElementType AdditiveIdentity() const override { return ElementType{}; }

  bool IsZero(const ElementType &a) const override {
    if (IsCanonical(a)) {
      return a.IsZero();
    }
    return SetElementValue(a).IsZero();
  }

  bool IsOne(const ElementType &a) const override {
    if (IsCanonical(a)) {
      return a.Degree() == 0 && a[0] == BaseFieldElementType{1};
    }
    const ElementType na = SetElementValue(a);
    return na.Degree() == 0 && na[0] == BaseFieldElementType{1};
  }

  ElementType MultiplicativeGenerator() const override {
    {
      std::lock_guard<std::mutex> lock(generator_cache_.mutex);
      if (generator_cache_.initialized) {
        return generator_cache_.value;
      }
    }

    ElementType gen = ComputeMultiplicativeGenerator();

    std::lock_guard<std::mutex> lock(generator_cache_.mutex);
    if (!generator_cache_.initialized) {
      generator_cache_.value = gen;
      generator_cache_.initialized = true;
    }
    return generator_cache_.value;
  }

  std::vector<ElementType> MultiplicativeGenerators() const override {
    if (order_ == 2) {
      return {MultiplicativeIdentity()};
    }
    constexpr uint64_t kMaxEnumeratedOrder = 1ULL << 20;
    if (order_ > kMaxEnumeratedOrder) {
      throw std::length_error(
          "Field order too large to enumerate all multiplicative generators");
    }

    const uint64_t group_order = order_ - 1;
    const ElementType gen = MultiplicativeGenerator();
    std::vector<ElementType> generators;

    ElementType cur = gen;
    for (uint64_t k = 1; k <= group_order; ++k) {
      if (std::gcd(k, group_order) == 1) {
        generators.push_back(cur);
      }
      if (k < group_order) {
        cur = Mul(cur, gen);
      }
    }
    return generators;
  }

  ElementType FromCoeffs(std::vector<BaseFieldElementType> coeffs) const {
    for (BaseFieldElementType &c : coeffs) {
      c = base_field_.SetElementValue(c);
    }
    ElementType elem(std::move(coeffs));
    if (elem.Degree() >= static_cast<int>(degree_)) {
      return PolyMod(elem, irreducible_poly_);
    }
    return elem;
  }

  // Base-p digit encoding: val = sum_{i=0}^{m-1} c_i * p^i (reduced mod p^m).
  ElementType FromInt(uint64_t val) const {
    val %= order_;
    if (val == 0) {
      return AdditiveIdentity();
    }
    const uint64_t p = base_field_.Characteristic();
    std::vector<BaseFieldElementType> coeffs;
    coeffs.reserve(degree_);
    while (val > 0) {
      coeffs.push_back(static_cast<BaseFieldElementType>(val % p));
      val /= p;
    }
    return ElementType(std::move(coeffs));
  }

  uint64_t ToInt(const ElementType &a) const {
    const ElementType canonical = IsCanonical(a) ? a : SetElementValue(a);
    const uint64_t p = base_field_.Characteristic();
    uint64_t result = 0;
    for (size_t i = canonical.Size(); i > 0; --i) {
      result = result * p + static_cast<uint64_t>(canonical[i - 1]);
    }
    return result;
  }

  ElementType GetElementValue(const ElementType &value) const override {
    return IsCanonical(value) ? value : SetElementValue(value);
  }

  ElementType SetElementValue(const ElementType &value) const override {
    std::vector<BaseFieldElementType> coeffs = value.Coeffs();
    return FromCoeffs(std::move(coeffs));
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    const std::string s = TrimWhitespace(value_str);
    if (s.empty()) {
      throw std::invalid_argument("Empty element string");
    }

    // 1. Hexadecimal literal: 0x...
    if (s.size() >= 3 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
      for (size_t i = 2; i < s.size(); ++i) {
        if (std::isxdigit(static_cast<unsigned char>(s[i])) == 0) {
          throw std::invalid_argument("Invalid hexadecimal element string: " +
                                      value_str);
        }
      }
      try {
        size_t parsed_len = 0;
        const uint64_t val = std::stoull(s.substr(2), &parsed_len, 16);
        if (parsed_len != s.size() - 2) {
          throw std::invalid_argument("Invalid hexadecimal element string: " +
                                      value_str);
        }
        return FromInt(val);
      } catch (const std::out_of_range &) {
        throw std::invalid_argument("Hexadecimal element out of range: " +
                                    value_str);
      }
    }

    // 2. Generator power string: "<generator_name>^<signed_int>"
    const std::string gen_prefix = generator_name_ + "^";
    if (s.rfind(gen_prefix, 0) == 0 &&
        (generator_name_ != variable_name_ ||
         (s.find('+') == std::string::npos &&
          s.find('*', gen_prefix.size()) == std::string::npos))) {
      return ParsePowerElement(s);
    }

    // 3. Plain non-negative decimal integer (base-p encoding; matches constant
    // polynomial for all values < p and round-trips INT representation).
    bool all_digits = true;
    for (char ch : s) {
      if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
        all_digits = false;
        break;
      }
    }
    if (all_digits) {
      try {
        size_t parsed_len = 0;
        const uint64_t val = std::stoull(s, &parsed_len, 10);
        if (parsed_len != s.size()) {
          throw std::invalid_argument("Invalid integer element string: " +
                                      value_str);
        }
        return FromInt(val);
      } catch (const std::out_of_range &) {
        throw std::invalid_argument("Integer element out of range: " +
                                    value_str);
      }
    }

    // 4. Polynomial expression in variable_name_ (or x / X).
    ElementType parsed = ParsePolynomialString(base_field_, s, variable_name_);
    return SetElementValue(parsed);
  }

  void Print(std::ostream &os) const override {
    os << "GaloisFieldExtension GF(" << base_field_.Characteristic() << "^"
       << degree_ << ") over GF(" << base_field_.Characteristic()
       << ") with modulus " << PolyToString(irreducible_poly_)
       << " (elements are polynomials in " << variable_name_ << ")";

    switch (GetRepresentation()) {
      case FieldRepresentation::POLY:
        os << " [Rep: POLY]";
        break;
      case FieldRepresentation::INT:
        os << " [Rep: INT]";
        break;
      case FieldRepresentation::HEX:
        os << " [Rep: HEX]";
        break;
      case FieldRepresentation::POW:
        os << " [Rep: POW]";
        break;
      case FieldRepresentation::LOG:
        os << " [Rep: LOG]";
        break;
    }
  }

  void Print(const ElementType &a, std::ostream &os) const override {
    switch (GetRepresentation()) {
      case FieldRepresentation::POLY:
        os << PolyToString(a);
        break;
      case FieldRepresentation::INT:
        os << ToInt(a);
        break;
      case FieldRepresentation::HEX: {
        std::ostringstream oss;
        oss << "0x" << std::hex << ToInt(a);
        os << oss.str();
        break;
      }
      case FieldRepresentation::POW: {
        const ElementType na = IsCanonical(a) ? a : SetElementValue(a);
        if (na.IsZero()) {
          os << "0";
        } else {
          os << generator_name_ << "^" << Log(na);
        }
        break;
      }
      case FieldRepresentation::LOG: {
        const ElementType na = IsCanonical(a) ? a : SetElementValue(a);
        if (na.IsZero()) {
          os << "undefined";
        } else {
          os << Log(na);
        }
        break;
      }
    }
  }

  std::string ToString(const ElementType &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_.load(std::memory_order_relaxed);
  }

  void SetRepresentation(FieldRepresentation rep) override {
    representation_.store(rep, std::memory_order_relaxed);
  }

  size_t MemoryBytes() const override {
    size_t total = sizeof(*this);
    total +=
        irreducible_poly_.Coeffs().capacity() * sizeof(BaseFieldElementType);
    total += variable_name_.capacity();
    total += generator_name_.capacity();
    return total;
  }

  // Exact Euclidean polynomial division over the base field F_p.
  std::pair<ElementType, ElementType> PolyDivRem(
      const ElementType &dividend, const ElementType &divisor) const {
    if (divisor.IsZero()) {
      throw std::domain_error("Polynomial division by zero");
    }
    const int d_deg = divisor.Degree();
    if (dividend.Degree() < d_deg) {
      return {ElementType{}, dividend};
    }

    std::vector<BaseFieldElementType> rem = dividend.Coeffs();
    const BaseFieldElementType inv_lc =
        base_field_.Inv(divisor[static_cast<size_t>(d_deg)]);
    const size_t q_len = static_cast<size_t>(dividend.Degree() - d_deg) + 1;
    std::vector<BaseFieldElementType> quot(q_len, BaseFieldElementType{0});

    for (int i = dividend.Degree(); i >= d_deg; --i) {
      const BaseFieldElementType c = rem[static_cast<size_t>(i)];
      if (c == BaseFieldElementType{0}) {
        continue;
      }
      const BaseFieldElementType q_coeff = base_field_.Mul(c, inv_lc);
      const size_t shift = static_cast<size_t>(i - d_deg);
      quot[shift] = q_coeff;
      for (int j = 0; j <= d_deg; ++j) {
        const BaseFieldElementType prod =
            base_field_.Mul(q_coeff, divisor[static_cast<size_t>(j)]);
        rem[shift + static_cast<size_t>(j)] =
            base_field_.Sub(rem[shift + static_cast<size_t>(j)], prod);
      }
    }

    return {ElementType(std::move(quot)), ElementType(std::move(rem))};
  }

 private:
  struct GeneratorCache {
    mutable std::mutex mutex;
    mutable bool initialized{false};
    mutable ElementType value;

    GeneratorCache() = default;
    GeneratorCache(const GeneratorCache &other) {
      std::lock_guard<std::mutex> lock(other.mutex);
      initialized = other.initialized;
      value = other.value;
    }
    GeneratorCache &operator=(const GeneratorCache &other) {
      if (this != &other) {
        std::scoped_lock lock(mutex, other.mutex);
        initialized = other.initialized;
        value = other.value;
      }
      return *this;
    }
  };

  static BaseFieldElementType ValidatePrime(uint64_t prime,
                                            bool prime_testing) {
    if (prime < 2) {
      throw std::invalid_argument("Prime must be at least 2");
    }
    if (prime > static_cast<uint64_t>(
                    std::numeric_limits<BaseFieldElementType>::max())) {
      throw std::invalid_argument(
          "Prime exceeds capacity of BaseFieldElementType");
    }
    if (prime_testing && !utils::IsPrime(prime)) {
      throw std::invalid_argument("Characteristic must be prime");
    }
    return static_cast<BaseFieldElementType>(prime);
  }

  static uint32_t ValidateExponent(uint64_t exponent) {
    if (exponent < 1) {
      throw std::invalid_argument("Extension degree must be at least 1");
    }
    if (exponent >
        static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
      throw std::invalid_argument("Extension degree exceeds uint32_t");
    }
    return static_cast<uint32_t>(exponent);
  }

  static std::pair<uint64_t, uint64_t> DecomposeOrderToPrimeExp(
      uint64_t order) {
    if (order < 2) {
      throw std::invalid_argument("Order of finite field must be at least 2");
    }
    auto [prime, exponent] = utils::DecomposePrimePower(order);
    if (prime == 0) {
      throw std::invalid_argument("Order must be a prime power");
    }
    return {prime, exponent};
  }

  static std::string TrimWhitespace(std::string_view sv) {
    size_t start = 0;
    while (start < sv.size() &&
           std::isspace(static_cast<unsigned char>(sv[start])) != 0) {
      ++start;
    }
    size_t end = sv.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(sv[end - 1])) != 0) {
      --end;
    }
    return std::string(sv.substr(start, end - start));
  }

  static ElementType ParsePolynomialString(const BaseField &base_field,
                                           const std::string &poly_str,
                                           const std::string &var_name) {
    const std::string s = TrimWhitespace(poly_str);
    if (s.empty()) {
      throw std::invalid_argument("Polynomial string cannot be empty");
    }

    size_t pos = 0;
    auto skip_spaces = [&]() {
      while (pos < s.size() &&
             std::isspace(static_cast<unsigned char>(s[pos])) != 0) {
        ++pos;
      }
    };

    auto match_variable = [&]() -> bool {
      if (!var_name.empty() && s.compare(pos, var_name.size(), var_name) == 0) {
        pos += var_name.size();
        return true;
      }
      if (pos < s.size() && (s[pos] == 'x' || s[pos] == 'X')) {
        ++pos;
        return true;
      }
      return false;
    };

    auto parse_uint = [&]() -> uint64_t {
      const size_t start = pos;
      while (pos < s.size() &&
             std::isdigit(static_cast<unsigned char>(s[pos])) != 0) {
        ++pos;
      }
      if (pos == start) {
        throw std::invalid_argument("Expected integer in polynomial: " +
                                    poly_str);
      }
      try {
        size_t parsed_len = 0;
        const uint64_t val =
            std::stoull(s.substr(start, pos - start), &parsed_len, 10);
        if (parsed_len != pos - start) {
          throw std::invalid_argument("Invalid integer in polynomial: " +
                                      poly_str);
        }
        return val;
      } catch (const std::out_of_range &) {
        throw std::invalid_argument("Integer out of range in polynomial: " +
                                    poly_str);
      }
    };

    constexpr uint64_t kMaxParsedDegree = 100000;
    const uint64_t p = base_field.Characteristic();
    std::vector<BaseFieldElementType> coeffs;
    bool first_term = true;

    while (pos < s.size()) {
      skip_spaces();
      if (pos >= s.size()) {
        break;
      }

      bool negative = false;
      if (s[pos] == '+') {
        ++pos;
        skip_spaces();
        if (pos >= s.size()) {
          throw std::invalid_argument("Trailing '+' in polynomial: " +
                                      poly_str);
        }
      } else if (s[pos] == '-') {
        negative = true;
        ++pos;
        skip_spaces();
        if (pos >= s.size()) {
          throw std::invalid_argument("Trailing '-' in polynomial: " +
                                      poly_str);
        }
      } else if (!first_term) {
        throw std::invalid_argument("Expected '+' or '-' in polynomial: " +
                                    poly_str);
      }

      bool has_coeff = false;
      uint64_t raw_coeff = 1;
      if (pos < s.size() &&
          std::isdigit(static_cast<unsigned char>(s[pos])) != 0) {
        raw_coeff = parse_uint();
        has_coeff = true;
        skip_spaces();
      }

      bool has_star = false;
      if (has_coeff && pos < s.size() && s[pos] == '*') {
        has_star = true;
        ++pos;
        skip_spaces();
      }

      uint64_t deg = 0;
      if (match_variable()) {
        deg = 1;
        skip_spaces();
        if (pos < s.size() && s[pos] == '^') {
          ++pos;
          skip_spaces();
          deg = parse_uint();
          if (deg > kMaxParsedDegree) {
            throw std::invalid_argument("Polynomial degree too large: " +
                                        poly_str);
          }
        }
      } else {
        if (has_star || !has_coeff) {
          throw std::invalid_argument("Invalid term in polynomial: " +
                                      poly_str);
        }
      }

      BaseFieldElementType term_coeff =
          static_cast<BaseFieldElementType>(raw_coeff % p);
      if (negative) {
        term_coeff = base_field.Neg(term_coeff);
      }

      if (coeffs.size() <= static_cast<size_t>(deg)) {
        coeffs.resize(static_cast<size_t>(deg) + 1, BaseFieldElementType{0});
      }
      coeffs[static_cast<size_t>(deg)] =
          base_field.Add(coeffs[static_cast<size_t>(deg)], term_coeff);

      first_term = false;
      skip_spaces();
    }

    return ElementType(std::move(coeffs));
  }

  static std::pair<ElementType, bool> CreateIrreduciblePolynomial(
      const BaseField &base_field, const std::string &irreducible_poly_str,
      uint32_t degree, const std::string &var_name) {
    if (!irreducible_poly_str.empty()) {
      ElementType poly =
          ParsePolynomialString(base_field, irreducible_poly_str, var_name);
      return {std::move(poly), false};
    }

    const uint64_t p = base_field.Characteristic();
    if (p <= static_cast<uint64_t>(std::numeric_limits<int>::max()) &&
        degree <= static_cast<uint32_t>(std::numeric_limits<int>::max())) {
      try {
        const std::string conway_poly = databases::GetConwayPolynomial(
            static_cast<int>(p), static_cast<int>(degree));
        return {ParsePolynomialString(base_field, conway_poly, var_name), true};
      } catch (const std::runtime_error &) {
        try {
          const std::string irred_poly = databases::GetIrreduciblePolynomial(
              static_cast<int>(p), static_cast<int>(degree));
          return {ParsePolynomialString(base_field, irred_poly, var_name),
                  false};
        } catch (const std::runtime_error &) {
          // Fall through to degree-1 fallback or error.
        }
      }
    }

    if (degree == 1) {
      return {ElementType({BaseFieldElementType{0}, BaseFieldElementType{1}}),
              false};
    }

    throw std::invalid_argument(
        "No suitable polynomial found in databases for GF(" +
        std::to_string(p) + "^" + std::to_string(degree) +
        "). Please provide modulus polynomial manually.");
  }

  ElementType ParsePowerElement(const std::string &s) const {
    const size_t caret_pos = s.find('^');
    if (caret_pos == std::string::npos || caret_pos == 0) {
      throw std::invalid_argument("Invalid power format: " + s);
    }
    const std::string base = s.substr(0, caret_pos);
    if (base != generator_name_) {
      throw std::invalid_argument("Unexpected generator in power format: " + s);
    }

    std::string exp_str = s.substr(caret_pos + 1);
    if (exp_str.empty()) {
      throw std::invalid_argument("Invalid power format: " + s);
    }

    bool is_negative = false;
    if (exp_str[0] == '-') {
      is_negative = true;
      exp_str = exp_str.substr(1);
      if (exp_str.empty()) {
        throw std::invalid_argument("Invalid power format: " + s);
      }
    }

    for (char ch : exp_str) {
      if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
        throw std::invalid_argument("Invalid exponent format: " + s);
      }
    }

    uint64_t raw_exp = 0;
    try {
      size_t parsed_len = 0;
      raw_exp = std::stoull(exp_str, &parsed_len, 10);
      if (parsed_len != exp_str.size()) {
        throw std::invalid_argument("Invalid exponent format: " + s);
      }
    } catch (const std::out_of_range &) {
      throw std::invalid_argument("Exponent out of range: " + s);
    }

    const uint64_t group_order = order_ - 1;
    uint64_t reduced_exp = raw_exp % group_order;
    if (is_negative && reduced_exp != 0) {
      reduced_exp = group_order - reduced_exp;
    }

    return Pow(MultiplicativeGenerator(), reduced_exp);
  }

  bool IsCanonical(const ElementType &a) const noexcept {
    if (a.Degree() >= static_cast<int>(degree_)) {
      return false;
    }
    const uint64_t p = base_field_.Characteristic();
    for (BaseFieldElementType c : a.Coeffs()) {
      if (static_cast<uint64_t>(c) >= p) {
        return false;
      }
    }
    return true;
  }

  ElementType PolyAdd(const ElementType &a, const ElementType &b) const {
    const size_t max_len = std::max(a.Size(), b.Size());
    std::vector<BaseFieldElementType> res(max_len, BaseFieldElementType{0});
    for (size_t i = 0; i < max_len; ++i) {
      res[i] = base_field_.Add(a[i], b[i]);
    }
    return ElementType(std::move(res));
  }

  ElementType PolySub(const ElementType &a, const ElementType &b) const {
    const size_t max_len = std::max(a.Size(), b.Size());
    std::vector<BaseFieldElementType> res(max_len, BaseFieldElementType{0});
    for (size_t i = 0; i < max_len; ++i) {
      res[i] = base_field_.Sub(a[i], b[i]);
    }
    return ElementType(std::move(res));
  }

  ElementType PolyNeg(const ElementType &a) const {
    std::vector<BaseFieldElementType> res(a.Size(), BaseFieldElementType{0});
    for (size_t i = 0; i < a.Size(); ++i) {
      res[i] = base_field_.Neg(a[i]);
    }
    return ElementType(std::move(res));
  }

  ElementType PolyScale(const ElementType &a,
                        BaseFieldElementType scalar) const {
    if (scalar == BaseFieldElementType{0} || a.IsZero()) {
      return ElementType{};
    }
    std::vector<BaseFieldElementType> res(a.Size(), BaseFieldElementType{0});
    for (size_t i = 0; i < a.Size(); ++i) {
      res[i] = base_field_.Mul(a[i], scalar);
    }
    return ElementType(std::move(res));
  }

  ElementType PolyMulRaw(const ElementType &a, const ElementType &b) const {
    if (a.IsZero() || b.IsZero()) {
      return ElementType{};
    }
    std::vector<BaseFieldElementType> res(a.Size() + b.Size() - 1,
                                          BaseFieldElementType{0});
    for (size_t i = 0; i < a.Size(); ++i) {
      const BaseFieldElementType ai = a[i];
      if (ai == BaseFieldElementType{0}) {
        continue;
      }
      for (size_t j = 0; j < b.Size(); ++j) {
        const BaseFieldElementType bj = b[j];
        if (bj == BaseFieldElementType{0}) {
          continue;
        }
        res[i + j] = base_field_.Add(res[i + j], base_field_.Mul(ai, bj));
      }
    }
    return ElementType(std::move(res));
  }

  ElementType PolyMod(const ElementType &a, const ElementType &mod) const {
    return PolyDivRem(a, mod).second;
  }

  ElementType PolyPowMod(ElementType base, uint64_t exp,
                         const ElementType &mod) const {
    ElementType result({BaseFieldElementType{1}});
    base = PolyMod(base, mod);
    while (exp > 0) {
      if (exp & 1ULL) {
        result = PolyMod(PolyMulRaw(result, base), mod);
      }
      exp >>= 1;
      if (exp > 0) {
        base = PolyMod(PolyMulRaw(base, base), mod);
      }
    }
    return result;
  }

  ElementType PolyGcd(ElementType a, ElementType b) const {
    while (!b.IsZero()) {
      ElementType r = PolyMod(a, b);
      a = std::move(b);
      b = std::move(r);
    }
    if (!a.IsZero()) {
      const BaseFieldElementType lc = a[static_cast<size_t>(a.Degree())];
      if (lc != BaseFieldElementType{1}) {
        a = PolyScale(a, base_field_.Inv(lc));
      }
    }
    return a;
  }

  std::tuple<ElementType, ElementType, ElementType> PolyExtendedGcd(
      const ElementType &a, const ElementType &b) const {
    ElementType r0 = a;
    ElementType r1 = b;
    ElementType s0({BaseFieldElementType{1}});
    ElementType s1;
    ElementType t0;
    ElementType t1({BaseFieldElementType{1}});

    while (!r1.IsZero()) {
      auto [q, r_next] = PolyDivRem(r0, r1);
      r0 = std::move(r1);
      r1 = std::move(r_next);

      ElementType s_next = PolySub(s0, PolyMulRaw(q, s1));
      s0 = std::move(s1);
      s1 = std::move(s_next);

      ElementType t_next = PolySub(t0, PolyMulRaw(q, t1));
      t0 = std::move(t1);
      t1 = std::move(t_next);
    }

    if (!r0.IsZero()) {
      const BaseFieldElementType lc = r0[static_cast<size_t>(r0.Degree())];
      if (lc != BaseFieldElementType{1}) {
        const BaseFieldElementType inv_lc = base_field_.Inv(lc);
        r0 = PolyScale(r0, inv_lc);
        s0 = PolyScale(s0, inv_lc);
        t0 = PolyScale(t0, inv_lc);
      }
    }
    return {std::move(r0), std::move(s0), std::move(t0)};
  }

  // Exact Rabin irreducibility test over F_p.
  bool IsIrreducibleOverBase(const ElementType &f) const {
    const int deg = f.Degree();
    if (deg <= 0) {
      return false;
    }
    if (deg == 1) {
      return true;
    }

    const uint64_t p = base_field_.Characteristic();
    const uint64_t m = static_cast<uint64_t>(deg);
    const std::vector<uint64_t> prime_divisors = utils::DistinctPrimeFactors(m);

    std::vector<uint64_t> check_steps;
    check_steps.reserve(prime_divisors.size());
    for (uint64_t r : prime_divisors) {
      check_steps.push_back(m / r);
    }

    const ElementType x_poly(
        {BaseFieldElementType{0}, BaseFieldElementType{1}});
    ElementType h = x_poly;

    for (uint64_t step = 1; step <= m; ++step) {
      h = PolyPowMod(h, p, f);
      if (std::find(check_steps.begin(), check_steps.end(), step) !=
          check_steps.end()) {
        const ElementType diff = PolySub(h, x_poly);
        const ElementType g = PolyGcd(diff, f);
        if (g.Degree() != 0) {
          return false;
        }
      }
    }

    return h == x_poly;
  }

  bool IsPrimitiveElement(const ElementType &cand, uint64_t group_order,
                          const std::vector<uint64_t> &prime_factors) const {
    if (cand.IsZero()) {
      return false;
    }
    if (group_order == 1) {
      return IsOne(cand);
    }
    for (uint64_t factor : prime_factors) {
      if (IsOne(Pow(cand, group_order / factor))) {
        return false;
      }
    }
    return true;
  }

  ElementType ComputeMultiplicativeGenerator() const {
    if (order_ == 2) {
      return MultiplicativeIdentity();
    }

    const uint64_t group_order = order_ - 1;
    const std::vector<uint64_t> prime_factors =
        utils::DistinctPrimeFactors(group_order);

    // Fast path: when constructed from a Conway polynomial, alpha = {0, 1} is
    // guaranteed to be a primitive element.
    if (degree_ >= 2) {
      const ElementType alpha(
          {BaseFieldElementType{0}, BaseFieldElementType{1}});
      if (used_conway_poly_ ||
          IsPrimitiveElement(alpha, group_order, prime_factors)) {
        return alpha;
      }
    } else {
      const ElementType base_gen({base_field_.MultiplicativeGenerator()});
      if (IsPrimitiveElement(base_gen, group_order, prime_factors)) {
        return base_gen;
      }
    }

    const uint64_t start =
        (degree_ >= 2) ? base_field_.Characteristic() : uint64_t{1};
    for (uint64_t k = start; k < order_; ++k) {
      ElementType candidate = FromInt(k);
      if (IsPrimitiveElement(candidate, group_order, prime_factors)) {
        return candidate;
      }
    }

    throw std::runtime_error(
        "Failed to find multiplicative generator in GaloisFieldExtension");
  }

  std::string PolyToString(const ElementType &poly) const {
    if (poly.IsZero()) {
      return "0";
    }
    std::ostringstream oss;
    bool first = true;
    for (int d = poly.Degree(); d >= 0; --d) {
      const BaseFieldElementType c = poly[static_cast<size_t>(d)];
      if (c == BaseFieldElementType{0}) {
        continue;
      }
      if (!first) {
        oss << " + ";
      }
      if (d == 0) {
        oss << static_cast<uint64_t>(c);
      } else {
        if (c != BaseFieldElementType{1}) {
          oss << static_cast<uint64_t>(c) << "*";
        }
        oss << variable_name_;
        if (d > 1) {
          oss << "^" << d;
        }
      }
      first = false;
    }
    return oss.str();
  }

  BaseField base_field_;
  uint32_t degree_{1};
  ElementType irreducible_poly_;
  uint64_t order_{0};
  std::atomic<FieldRepresentation> representation_{FieldRepresentation::POLY};
  std::string variable_name_{"α"};
  std::string generator_name_{"g"};
  bool used_conway_poly_{false};
  GeneratorCache generator_cache_;
};

template <UnsignedFieldWord BaseFieldElementType = uint32_t>
using GFPX = GaloisFieldExtension<BaseFieldElementType>;

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_EXTENSION_HPP_
