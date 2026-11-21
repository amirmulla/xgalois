#ifndef XGALOIS_UTILS_POLY_HPP_
#define XGALOIS_UTILS_POLY_HPP_

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace utils {

namespace detail {

template <typename GaloisField>
inline void ValidateSamePolyField(const PolynomialDense<GaloisField> &a,
                                  const PolynomialDense<GaloisField> &b) {
  if (!a.HasField() || !b.HasField()) {
    throw std::invalid_argument(
        "Polynomials must have an associated Galois field.");
  }
  if (a.Field() != b.Field()) {
    throw std::invalid_argument(
        "Polynomials must belong to the same Galois field.");
  }
}

// Computes base^q mod mod where q = field->Order() (also handles Order() == 0
// denoting GF(2^64) via 64 repeated squarings).
template <typename GaloisField>
inline PolynomialDense<GaloisField> PowQMod(
    const PolynomialDense<GaloisField> &base,
    const PolynomialDense<GaloisField> &mod) {
  const auto field = mod.Field();
  const uint64_t q = field->Order();
  if (q != 0) {
    return base.PowMod(q, mod);
  }
  // Order() == 0 indicates GF(2^64).
  PolynomialDense<GaloisField> cur = base % mod;
  const uint32_t m = field->Degree();
  for (uint32_t i = 0; i < m; ++i) {
    cur *= cur;
    cur %= mod;
  }
  return cur;
}

// Computes exact p-th root of a raw field element in F_{p^m}: a^(1/p) =
// a^(p^(m-1)), using m - 1 Frobenius applications.
template <typename GaloisField>
inline typename GaloisField::element_type RawPthRoot(
    const GaloisField &field, const typename GaloisField::element_type &a) {
  const uint32_t m = field.Degree();
  if (m <= 1 || field.IsZero(a) || field.IsOne(a)) {
    return a;
  }
  typename GaloisField::element_type cur = a;
  for (uint32_t i = 1; i < m; ++i) {
    cur = field.Frobenius(cur);
  }
  return cur;
}

// Computes M = q^n - 1 if it fits in uint64_t (including 2^64 - 1 = UINT64_MAX)
// and throws std::overflow_error otherwise.
template <typename GaloisField>
inline uint64_t ComputeExtensionGroupOrder(const GaloisField &field,
                                           uint32_t n) {
  __extension__ using u128 = unsigned __int128;
  const uint64_t q = field.Order();
  if (n == 0) {
    return 0;
  }
  if (q == 0) {
    // Field order is 2^64.
    if (n == 1) {
      return std::numeric_limits<uint64_t>::max();
    }
    throw std::overflow_error("q^n - 1 exceeds uint64_t range.");
  }
  // Try SafeIntegerPower for q^n when q^n < 2^64; fall back to u128 check for
  // the exact 2^64 boundary (e.g., GF(2), n = 64).
  constexpr u128 kMaxPlusOne = (static_cast<u128>(1) << 64);
  u128 acc = 1;
  u128 base = static_cast<u128>(q);
  uint32_t exp = n;
  while (exp > 0) {
    if ((exp & 1U) != 0) {
      if (acc > kMaxPlusOne / base) {
        throw std::overflow_error("q^n - 1 exceeds uint64_t range.");
      }
      acc *= base;
    }
    exp >>= 1;
    if (exp > 0) {
      if (base > kMaxPlusOne / base) {
        throw std::overflow_error("q^n - 1 exceeds uint64_t range.");
      }
      base *= base;
    }
  }
  return static_cast<uint64_t>(acc - 1);
}

inline size_t ParseNonNegativeDegree(std::string_view suffix,
                                     std::string_view full_term) {
  if (suffix.empty()) {
    return 1;
  }
  if (suffix[0] != '^') {
    throw std::invalid_argument("Invalid exponent syntax in term: " +
                                std::string(full_term));
  }
  const std::string_view exp_sv = suffix.substr(1);
  if (exp_sv.empty()) {
    throw std::invalid_argument("Missing exponent in term: " +
                                std::string(full_term));
  }
  for (char ch : exp_sv) {
    if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
      throw std::invalid_argument("Invalid exponent in term: " +
                                  std::string(full_term));
    }
  }
  try {
    size_t parsed = 0;
    const unsigned long long val =
        std::stoull(std::string(exp_sv), &parsed, 10);
    if (parsed != exp_sv.size() ||
        val > static_cast<unsigned long long>(std::numeric_limits<int>::max() -
                                              1)) {
      throw std::invalid_argument("Exponent out of range in term: " +
                                  std::string(full_term));
    }
    return static_cast<size_t>(val);
  } catch (const std::invalid_argument &) {
    throw;
  } catch (const std::exception &) {
    throw std::invalid_argument("Exponent out of range in term: " +
                                std::string(full_term));
  }
}

inline size_t FindVariablePos(std::string_view body, std::string_view var) {
  size_t search_start = 0;
  // Do not match 'x'/'X' inside a leading hexadecimal literal "0x...".
  if (var == "x" && body.size() >= 2 && body[0] == '0' &&
      (body[1] == 'x' || body[1] == 'X')) {
    search_start = 2;
    while (search_start < body.size() &&
           std::isxdigit(static_cast<unsigned char>(body[search_start])) != 0) {
      ++search_start;
    }
  }

  if (var == "x") {
    for (size_t i = search_start; i < body.size(); ++i) {
      if (body[i] == 'x' || body[i] == 'X') {
        return i;
      }
    }
    return std::string_view::npos;
  }
  return body.find(var, search_start);
}

template <typename GaloisField>
inline typename GaloisField::element_type ParseFieldScalar(
    const GaloisField &field, std::string_view coeff_sv,
    std::string_view full_term) {
  if (coeff_sv.empty()) {
    return field.MultiplicativeIdentity();
  }

  // GaloisFieldBinary (GF(2)) does not implement SetElementValue(std::string);
  // parse decimal integer coefficients directly modulo 2.
  if (field.Order() == 2 && field.Characteristic() == 2 &&
      field.Degree() == 1) {
    size_t idx = 0;
    if (coeff_sv[0] == '+' || coeff_sv[0] == '-') {
      idx = 1;
      if (idx >= coeff_sv.size()) {
        throw std::invalid_argument("Invalid coefficient in term: " +
                                    std::string(full_term));
      }
    }
    for (size_t i = idx; i < coeff_sv.size(); ++i) {
      if (std::isdigit(static_cast<unsigned char>(coeff_sv[i])) == 0) {
        throw std::invalid_argument("Invalid GF(2) coefficient in term: " +
                                    std::string(full_term));
      }
    }
    const bool is_odd = ((coeff_sv.back() - '0') & 1) != 0;
    return is_odd ? field.MultiplicativeIdentity() : field.AdditiveIdentity();
  }

  try {
    return field.SetElementValue(std::string(coeff_sv));
  } catch (const std::invalid_argument &) {
    throw;
  } catch (const std::exception &e) {
    throw std::invalid_argument(std::string("Invalid coefficient in term '") +
                                std::string(full_term) + "': " + e.what());
  }
}

template <typename GaloisField>
inline bool PolyCanonicalLess(const PolynomialDense<GaloisField> &a,
                              const PolynomialDense<GaloisField> &b) {
  if (a.Degree() != b.Degree()) {
    return a.Degree() < b.Degree();
  }
  if (a.IsZero()) {
    return false;
  }
  const auto field = a.Field();
  for (int i = a.Degree(); i >= 0; --i) {
    const auto va = field->GetElementValue(a.RawCoeff(static_cast<size_t>(i)));
    const auto vb = field->GetElementValue(b.RawCoeff(static_cast<size_t>(i)));
    if (va != vb) {
      return va < vb;
    }
  }
  return false;
}

}  // namespace detail

template <typename GaloisField>
PolynomialDense<GaloisField> PolynomialGcd(
    const PolynomialDense<GaloisField> &a,
    const PolynomialDense<GaloisField> &b) {
  detail::ValidateSamePolyField(a, b);
  if (a.IsZero() && b.IsZero()) {
    return PolynomialDense<GaloisField>::Zero(a.Field(), a.GetVariable());
  }
  if (b.IsZero()) {
    return a.Monic();
  }
  if (a.IsZero()) {
    return b.Monic();
  }

  PolynomialDense<GaloisField> r0 = a;
  PolynomialDense<GaloisField> r1 = b;
  while (!r1.IsZero()) {
    r0 %= r1;
    std::swap(r0, r1);
  }
  return r0.Monic();
}

template <typename GaloisField>
PolynomialDense<GaloisField> PolynomialLcm(
    const PolynomialDense<GaloisField> &a,
    const PolynomialDense<GaloisField> &b) {
  detail::ValidateSamePolyField(a, b);
  if (a.IsZero() || b.IsZero()) {
    return PolynomialDense<GaloisField>::Zero(a.Field(), a.GetVariable());
  }
  const auto g = PolynomialGcd(a, b);
  return ((a / g) * b).Monic();
}

template <typename GaloisField>
std::pair<PolynomialDense<GaloisField>,
          std::pair<PolynomialDense<GaloisField>, PolynomialDense<GaloisField>>>
PolynomialDenseExtendedGcd(const PolynomialDense<GaloisField> &a,
                           const PolynomialDense<GaloisField> &b) {
  using ElementType = typename PolynomialDense<GaloisField>::ElementType;
  using PolynomialType = PolynomialDense<GaloisField>;

  detail::ValidateSamePolyField(a, b);
  const auto field = a.Field();
  const std::string var_name = a.GetVariable();

  const PolynomialType zero_poly = PolynomialType::Zero(field, var_name);
  const PolynomialType one_poly = PolynomialType::One(field, var_name);

  if (a.IsZero() && b.IsZero()) {
    return {zero_poly, {zero_poly, zero_poly}};
  }

  if (b.IsZero()) {
    const ElementType lc = a.LeadingCoeff();
    if (!lc.IsOne()) {
      const ElementType inv_lc = lc.Inv();
      return {a * inv_lc,
              {PolynomialType::FromRaw({inv_lc.RawValue()}, field, var_name),
               zero_poly}};
    }
    return {a, {one_poly, zero_poly}};
  }

  if (a.IsZero()) {
    const ElementType lc = b.LeadingCoeff();
    if (!lc.IsOne()) {
      const ElementType inv_lc = lc.Inv();
      return {b * inv_lc,
              {zero_poly,
               PolynomialType::FromRaw({inv_lc.RawValue()}, field, var_name)}};
    }
    return {b, {zero_poly, one_poly}};
  }

  PolynomialType r0 = a;
  PolynomialType r1 = b;
  PolynomialType s0 = one_poly;
  PolynomialType s1 = zero_poly;
  PolynomialType t0 = zero_poly;
  PolynomialType t1 = one_poly;

  while (!r1.IsZero()) {
    auto [q, r_next] = r0.DivRem(r1);
    r0 = std::move(r1);
    r1 = std::move(r_next);

    PolynomialType s_next = s0 - q * s1;
    s0 = std::move(s1);
    s1 = std::move(s_next);

    PolynomialType t_next = t0 - q * t1;
    t0 = std::move(t1);
    t1 = std::move(t_next);
  }

  if (!r0.IsZero()) {
    const ElementType lc = r0.LeadingCoeff();
    if (!lc.IsOne()) {
      const ElementType inv_lc = lc.Inv();
      r0 *= inv_lc;
      s0 *= inv_lc;
      t0 *= inv_lc;
    }
  } else {
    s0 = zero_poly;
    t0 = zero_poly;
  }

  return {std::move(r0), {std::move(s0), std::move(t0)}};
}

template <typename GaloisField>
static inline typename GaloisField::element_type ExtractCoefficient(
    const std::string &term, const std::string &var = "x") {
  using ElementType = typename GaloisField::element_type;
  if (term.empty()) {
    throw std::invalid_argument("Empty polynomial term.");
  }

  const size_t var_pos = detail::FindVariablePos(term, var);
  if (var_pos == std::string_view::npos) {
    size_t idx = 0;
    if (term[0] == '+' || term[0] == '-') {
      idx = 1;
      if (idx >= term.size()) {
        throw std::invalid_argument("Invalid coefficient in term: " + term);
      }
    }
    for (size_t i = idx; i < term.size(); ++i) {
      if (std::isdigit(static_cast<unsigned char>(term[i])) == 0) {
        throw std::invalid_argument("Invalid coefficient in term: " + term);
      }
    }
    try {
      size_t parsed = 0;
      const unsigned long long val = std::stoull(term.substr(idx), &parsed, 10);
      if (parsed != term.size() - idx) {
        throw std::invalid_argument("Invalid coefficient in term: " + term);
      }
      return static_cast<ElementType>(val);
    } catch (const std::invalid_argument &) {
      throw;
    } catch (const std::exception &) {
      throw std::invalid_argument("Invalid coefficient in term: " + term);
    }
  }

  // Validate exponent suffix after variable.
  (void)detail::ParseNonNegativeDegree(
      std::string_view(term).substr(var_pos + var.size()), term);

  std::string_view coeff_sv = std::string_view(term).substr(0, var_pos);
  if (!coeff_sv.empty() && coeff_sv.back() == '*') {
    if (coeff_sv.size() == 1 || coeff_sv == "+*" || coeff_sv == "-*") {
      throw std::invalid_argument("Invalid coefficient in term: " + term);
    }
    coeff_sv.remove_suffix(1);
  }

  if (coeff_sv.empty() || coeff_sv == "+" || coeff_sv == "-") {
    return ElementType{1};
  }

  size_t idx = 0;
  if (coeff_sv[0] == '+' || coeff_sv[0] == '-') {
    idx = 1;
  }
  for (size_t i = idx; i < coeff_sv.size(); ++i) {
    if (std::isdigit(static_cast<unsigned char>(coeff_sv[i])) == 0) {
      throw std::invalid_argument("Invalid coefficient in term: " + term);
    }
  }

  try {
    size_t parsed = 0;
    const unsigned long long val =
        std::stoull(std::string(coeff_sv.substr(idx)), &parsed, 10);
    if (parsed != coeff_sv.size() - idx) {
      throw std::invalid_argument("Invalid coefficient in term: " + term);
    }
    return static_cast<ElementType>(val);
  } catch (const std::invalid_argument &) {
    throw;
  } catch (const std::exception &) {
    throw std::invalid_argument("Invalid coefficient in term: " + term);
  }
}

static inline int ExtractDegree(const std::string &term,
                                const std::string &var = "x") {
  if (term.empty()) {
    throw std::invalid_argument("Empty polynomial term.");
  }
  const size_t var_pos = detail::FindVariablePos(term, var);
  if (var_pos == std::string_view::npos) {
    return 0;
  }
  const std::string_view suffix =
      std::string_view(term).substr(var_pos + var.size());
  return static_cast<int>(detail::ParseNonNegativeDegree(suffix, term));
}

template <typename GaloisField>
PolynomialDense<GaloisField> ParsePolynomial(
    const std::shared_ptr<GaloisField> &base_field, const std::string &poly_str,
    const std::string &variable_name) {
  using RawType = typename GaloisField::element_type;
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!base_field) {
    throw std::invalid_argument("GaloisField pointer cannot be null.");
  }
  if (variable_name.empty()) {
    throw std::invalid_argument("Variable name cannot be empty.");
  }
  if (poly_str.empty()) {
    throw std::invalid_argument("Polynomial string cannot be empty.");
  }

  std::string cleaned;
  cleaned.reserve(poly_str.size());
  for (char ch : poly_str) {
    if (std::isspace(static_cast<unsigned char>(ch)) == 0) {
      cleaned.push_back(ch);
    }
  }
  if (cleaned.empty()) {
    throw std::invalid_argument(
        "Polynomial string cannot be empty after cleaning.");
  }

  // Split into signed terms while respecting parenthesized field elements.
  struct SignedTerm {
    bool negative = false;
    std::string body;
  };
  std::vector<SignedTerm> terms;
  int paren_depth = 0;
  bool current_neg = false;
  std::string current_body;

  size_t start_idx = 0;
  if (cleaned[0] == '+') {
    start_idx = 1;
  } else if (cleaned[0] == '-') {
    current_neg = true;
    start_idx = 1;
  }
  if (start_idx >= cleaned.size()) {
    throw std::invalid_argument("Invalid polynomial string: " + poly_str);
  }

  for (size_t i = start_idx; i < cleaned.size(); ++i) {
    const char ch = cleaned[i];
    if (ch == '(') {
      ++paren_depth;
      current_body.push_back(ch);
    } else if (ch == ')') {
      --paren_depth;
      if (paren_depth < 0) {
        throw std::invalid_argument("Unmatched ')' in polynomial string: " +
                                    poly_str);
      }
      current_body.push_back(ch);
    } else if ((ch == '+' || ch == '-') && paren_depth == 0) {
      if (current_body.empty()) {
        throw std::invalid_argument("Malformed operator sequence in: " +
                                    poly_str);
      }
      terms.push_back({current_neg, std::move(current_body)});
      current_body.clear();
      current_neg = (ch == '-');
    } else {
      current_body.push_back(ch);
    }
  }

  if (paren_depth != 0) {
    throw std::invalid_argument("Unmatched '(' in polynomial string: " +
                                poly_str);
  }
  if (current_body.empty()) {
    throw std::invalid_argument("Trailing operator in polynomial string: " +
                                poly_str);
  }
  terms.push_back({current_neg, std::move(current_body)});

  struct ParsedTerm {
    RawType raw_coeff;
    size_t degree;
  };
  std::vector<ParsedTerm> parsed_terms;
  parsed_terms.reserve(terms.size());
  size_t max_degree = 0;

  for (const auto &term : terms) {
    const std::string_view body = term.body;
    std::string_view coeff_sv;
    size_t degree = 0;

    if (body[0] == '(') {
      // Find matching closing parenthesis.
      int depth = 0;
      size_t close_pos = std::string_view::npos;
      for (size_t i = 0; i < body.size(); ++i) {
        if (body[i] == '(') {
          ++depth;
        } else if (body[i] == ')') {
          --depth;
          if (depth == 0) {
            close_pos = i;
            break;
          }
        }
      }
      if (close_pos == std::string_view::npos || close_pos <= 1) {
        throw std::invalid_argument("Invalid parenthesized term: " + term.body);
      }
      coeff_sv = body.substr(1, close_pos - 1);
      std::string_view tail = body.substr(close_pos + 1);
      if (tail.empty()) {
        degree = 0;
      } else {
        if (tail[0] == '*') {
          tail.remove_prefix(1);
          if (tail.empty()) {
            throw std::invalid_argument("Trailing '*' in term: " + term.body);
          }
        }
        bool matches_var = false;
        if (variable_name == "x") {
          matches_var = (!tail.empty() && (tail[0] == 'x' || tail[0] == 'X'));
        } else {
          matches_var = (tail.size() >= variable_name.size() &&
                         tail.substr(0, variable_name.size()) == variable_name);
        }
        if (!matches_var) {
          throw std::invalid_argument("Unexpected token after ')' in term: " +
                                      term.body);
        }
        const std::string_view suffix = tail.substr(variable_name.size());
        degree = detail::ParseNonNegativeDegree(suffix, body);
      }
    } else {
      const size_t var_pos = detail::FindVariablePos(body, variable_name);
      if (var_pos == std::string_view::npos) {
        coeff_sv = body;
        degree = 0;
      } else {
        coeff_sv = body.substr(0, var_pos);
        if (!coeff_sv.empty() && coeff_sv.back() == '*') {
          if (coeff_sv.size() == 1) {
            throw std::invalid_argument("Leading '*' without coefficient: " +
                                        term.body);
          }
          coeff_sv.remove_suffix(1);
        }
        const std::string_view suffix =
            body.substr(var_pos + variable_name.size());
        degree = detail::ParseNonNegativeDegree(suffix, body);
      }
    }

    RawType raw_coeff = detail::ParseFieldScalar(*base_field, coeff_sv, body);
    if (term.negative) {
      raw_coeff = base_field->Neg(raw_coeff);
    }
    parsed_terms.push_back({std::move(raw_coeff), degree});
    if (degree > max_degree) {
      max_degree = degree;
    }
  }

  std::vector<RawType> coeffs(max_degree + 1, base_field->AdditiveIdentity());
  for (const auto &pt : parsed_terms) {
    coeffs[pt.degree] = base_field->Add(coeffs[pt.degree], pt.raw_coeff);
  }

  return PolynomialType::FromRaw(std::move(coeffs), base_field, variable_name);
}

template <typename GaloisField>
PolynomialDense<GaloisField> ParsePolynomial(
    const std::shared_ptr<GaloisField> &base_field,
    const std::string &poly_str) {
  return ParsePolynomial(base_field, poly_str, std::string("x"));
}

// Exact Rabin's irreducibility test over arbitrary finite fields F_q.
template <typename GaloisField>
bool IsIrreducible(const PolynomialDense<GaloisField> &poly) {
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Cannot determine field for irreducibility test.");
  }

  const int n = poly.Degree();
  if (n <= 0) {
    return false;
  }
  if (n == 1) {
    return true;
  }

  const auto field = poly.Field();
  const std::string var_name = poly.GetVariable();
  const PolynomialType f = poly.Monic();
  const PolynomialType x = PolynomialType::Monomial(
      field, 1, field->MultiplicativeIdentity(), var_name);

  // Distinct prime divisors of n; form sorted ascending targets n / p_i.
  const std::vector<uint64_t> prime_divisors =
      DistinctPrimeFactors(static_cast<uint64_t>(n));
  std::vector<uint64_t> check_exps;
  check_exps.reserve(prime_divisors.size());
  for (size_t i = prime_divisors.size(); i > 0; --i) {
    check_exps.push_back(static_cast<uint64_t>(n) / prime_divisors[i - 1]);
  }

  PolynomialType h = x;
  uint64_t cur_k = 0;
  for (uint64_t target_k : check_exps) {
    while (cur_k < target_k) {
      h = detail::PowQMod(h, f);
      ++cur_k;
    }
    const PolynomialType g = PolynomialGcd(f, h - x);
    if (g.Degree() > 0) {
      return false;
    }
  }

  while (cur_k < static_cast<uint64_t>(n)) {
    h = detail::PowQMod(h, f);
    ++cur_k;
  }
  return h == x;
}

// Exact primitivity test for polynomials over F_q: f is primitive iff it is
// irreducible of degree n >= 1, f(0) != 0, and x mod f has order q^n - 1 in
// (F_q[x]/(f))^*.
template <typename GaloisField>
bool IsPrimitive(const PolynomialDense<GaloisField> &poly) {
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument("Cannot determine field for primitivity test.");
  }

  const int n = poly.Degree();
  if (n <= 0) {
    return false;
  }

  const auto field = poly.Field();
  if (field->IsZero(poly.RawCoeff(0))) {
    return false;
  }

  if (!IsIrreducible(poly)) {
    return false;
  }

  const uint64_t group_order =
      detail::ComputeExtensionGroupOrder(*field, static_cast<uint32_t>(n));
  if (group_order <= 1) {
    // Only occurs for GF(2), deg == 1, f = x + 1 (since f(0) != 0).
    return true;
  }

  const PolynomialType f = poly.Monic();
  const PolynomialType x = PolynomialType::Monomial(
      field, 1, field->MultiplicativeIdentity(), poly.GetVariable());

  const std::vector<uint64_t> prime_factors = DistinctPrimeFactors(group_order);
  for (uint64_t p : prime_factors) {
    const PolynomialType rem = x.PowMod(group_order / p, f);
    if (rem.IsOne()) {
      return false;
    }
  }
  return true;
}

// Square-free factorization over F_q (Yun / Musser algorithm with
// characteristic p-th root extraction via inverse Frobenius). Returns pairs
// (monic_sqfree, multiplicity) with non-constant factors sorted by multiplicity
// ascending.
template <typename GaloisField>
std::vector<std::pair<PolynomialDense<GaloisField>, uint32_t>> SquareFreeFactor(
    const PolynomialDense<GaloisField> &poly) {
  using RawType = typename GaloisField::element_type;
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  if (poly.IsZero()) {
    throw std::invalid_argument(
        "Cannot compute square-free factorization of 0.");
  }
  if (poly.Degree() == 0) {
    return {};
  }

  const auto field = poly.Field();
  const std::string var_name = poly.GetVariable();
  const uint64_t p = field->Characteristic();

  PolynomialType f = poly.Monic();
  PolynomialType f_deriv = f.Derivative();
  PolynomialType c = PolynomialGcd(f, f_deriv);
  PolynomialType w = f / c;

  std::vector<std::pair<PolynomialType, uint32_t>> result;
  uint32_t i = 1;
  while (w.Degree() > 0) {
    PolynomialType y = PolynomialGcd(w, c);
    PolynomialType z = w / y;
    if (z.Degree() > 0) {
      result.emplace_back(z.Monic(), i);
    }
    w = std::move(y);
    c /= w;
    ++i;
  }

  // Any remaining non-constant part c is a p-th power in F_q[x].
  if (c.Degree() > 0) {
    const size_t deg_c = static_cast<size_t>(c.Degree());
    const size_t root_deg = deg_c / static_cast<size_t>(p);
    std::vector<RawType> root_coeffs(root_deg + 1, field->AdditiveIdentity());
    for (size_t k = 0; k <= root_deg; ++k) {
      const RawType coeff = c.RawCoeff(k * static_cast<size_t>(p));
      if (!field->IsZero(coeff)) {
        root_coeffs[k] = detail::RawPthRoot(*field, coeff);
      }
    }
    PolynomialType c_root =
        PolynomialType::FromRaw(std::move(root_coeffs), field, var_name);
    auto sub_factors = SquareFreeFactor(c_root);
    for (auto &[g, mult] : sub_factors) {
      result.emplace_back(std::move(g), mult * static_cast<uint32_t>(p));
    }
  }

  std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
    return lhs.second < rhs.second;
  });
  return result;
}

// Distinct-degree factorization of a square-free polynomial over F_q.
// Returns pairs (g_d, d) where g_d is the monic product of all irreducible
// factors of degree d.
template <typename GaloisField>
std::vector<std::pair<PolynomialDense<GaloisField>, uint32_t>>
DistinctDegreeFactor(const PolynomialDense<GaloisField> &poly) {
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  if (poly.IsZero()) {
    throw std::invalid_argument(
        "Cannot compute distinct-degree factorization of 0.");
  }
  if (poly.Degree() == 0) {
    return {};
  }

  const auto field = poly.Field();
  const std::string var_name = poly.GetVariable();
  PolynomialType f_star = poly.Monic();
  const PolynomialType x = PolynomialType::Monomial(
      field, 1, field->MultiplicativeIdentity(), var_name);
  PolynomialType h = x;

  std::vector<std::pair<PolynomialType, uint32_t>> result;
  uint32_t d = 1;
  while (f_star.Degree() >= static_cast<int>(2 * d)) {
    h = detail::PowQMod(h, f_star);
    PolynomialType g = PolynomialGcd(h - x, f_star);
    if (g.Degree() > 0) {
      result.emplace_back(g, d);
      f_star /= g;
      if (f_star.Degree() > 0) {
        h %= f_star;
      }
    }
    ++d;
  }

  if (f_star.Degree() > 0) {
    result.emplace_back(f_star, static_cast<uint32_t>(f_star.Degree()));
  }
  return result;
}

// Cantor-Zassenhaus equal-degree factorization of a square-free product of
// irreducible degree-d polynomials over F_q (supports both odd characteristic
// and characteristic 2 via the absolute trace map).
template <typename GaloisField>
std::vector<PolynomialDense<GaloisField>> EqualDegreeFactor(
    const PolynomialDense<GaloisField> &poly, uint32_t d,
    std::mt19937_64 &rng) {
  using RawType = typename GaloisField::element_type;
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  if (d == 0) {
    throw std::invalid_argument("Factor degree d must be >= 1.");
  }
  if (poly.Degree() < static_cast<int>(d) ||
      poly.Degree() % static_cast<int>(d) != 0) {
    throw std::invalid_argument(
        "Polynomial degree must be a positive multiple of d.");
  }

  const auto field = poly.Field();
  const std::string var_name = poly.GetVariable();
  const PolynomialType one_poly = PolynomialType::One(field, var_name);
  const uint64_t p = field->Characteristic();
  const uint64_t q = field->Order();
  const uint32_t ext_deg = field->Degree();

  std::vector<PolynomialType> worklist;
  worklist.push_back(poly.Monic());
  std::vector<PolynomialType> irreducibles;
  irreducibles.reserve(static_cast<size_t>(poly.Degree()) / d);

  while (!worklist.empty()) {
    PolynomialType u = std::move(worklist.back());
    worklist.pop_back();

    if (u.Degree() == static_cast<int>(d)) {
      irreducibles.push_back(std::move(u));
      continue;
    }

    const size_t deg_u = static_cast<size_t>(u.Degree());
    bool split = false;
    constexpr size_t kMaxAttempts = 10000;
    for (size_t attempt = 0; attempt < kMaxAttempts; ++attempt) {
      std::vector<RawType> rand_coeffs(deg_u);
      for (size_t i = 0; i < deg_u; ++i) {
        rand_coeffs[i] = field->Random(rng);
      }
      PolynomialType a =
          PolynomialType::FromRaw(std::move(rand_coeffs), field, var_name);
      if (a.Degree() <= 0) {
        continue;
      }

      PolynomialType g = PolynomialGcd(a, u);
      if (g.Degree() > 0 && g.Degree() < u.Degree()) {
        worklist.push_back(u / g);
        worklist.push_back(std::move(g));
        split = true;
        break;
      }

      PolynomialType splitter;
      if (p == 2) {
        // Characteristic 2: absolute trace Tr_{F_{2^(m*d)} / F_2}(a) mod u
        // = sum_{j=0}^{m*d - 1} a^(2^j) mod u.
        const uint32_t total_bits = ext_deg * d;
        PolynomialType cur = a % u;
        splitter = cur;
        for (uint32_t j = 1; j < total_bits; ++j) {
          cur *= cur;
          cur %= u;
          splitter += cur;
        }
      } else {
        // Odd characteristic: a^((q^d - 1) / 2) - 1 mod u.
        // Uses (q^d - 1)/2 = ((q - 1)/2) * (1 + q + q^2 + ... + q^(d-1))
        // to avoid any 64-bit integer overflow regardless of q and d.
        const PolynomialType t = a.PowMod((q - 1) / 2, u);
        PolynomialType acc = t;
        PolynomialType cur = t;
        for (uint32_t j = 1; j < d; ++j) {
          cur = detail::PowQMod(cur, u);
          acc *= cur;
          acc %= u;
        }
        splitter = acc - one_poly;
      }

      g = PolynomialGcd(u, splitter);
      if (g.Degree() > 0 && g.Degree() < u.Degree()) {
        worklist.push_back(u / g);
        worklist.push_back(std::move(g));
        split = true;
        break;
      }
    }

    if (!split) {
      throw std::runtime_error(
          "EqualDegreeFactor failed to split polynomial; verify input is "
          "square-free and a product of degree-d irreducibles.");
    }
  }

  std::sort(irreducibles.begin(), irreducibles.end(),
            detail::PolyCanonicalLess<GaloisField>);
  return irreducibles;
}

template <typename GaloisField>
std::vector<PolynomialDense<GaloisField>> EqualDegreeFactor(
    const PolynomialDense<GaloisField> &poly, uint32_t d,
    uint64_t seed = 0x9E3779B97F4A7C15ULL) {
  std::mt19937_64 rng(seed);
  return EqualDegreeFactor(poly, d, rng);
}

// Complete factorization of a nonzero polynomial f in F_q[x].
// Returns {leading_coeff, vector<{irreducible_monic_factor, multiplicity}>}.
template <typename GaloisField>
std::pair<GaloisFieldElement<GaloisField>,
          std::vector<std::pair<PolynomialDense<GaloisField>, uint32_t>>>
Factor(const PolynomialDense<GaloisField> &poly,
       uint64_t seed = 0x9E3779B97F4A7C15ULL) {
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  if (poly.IsZero()) {
    throw std::invalid_argument("Cannot factor the zero polynomial.");
  }

  const auto lc = poly.LeadingCoeff();
  if (poly.Degree() == 0) {
    return {lc, {}};
  }

  std::mt19937_64 rng(seed);
  std::vector<std::pair<PolynomialType, uint32_t>> factors;

  const auto sqfree = SquareFreeFactor(poly.Monic());
  for (const auto &[sq_poly, mult] : sqfree) {
    const auto ddf = DistinctDegreeFactor(sq_poly);
    for (const auto &[dd_poly, deg] : ddf) {
      auto edf = EqualDegreeFactor(dd_poly, deg, rng);
      for (auto &irred : edf) {
        factors.emplace_back(std::move(irred), mult);
      }
    }
  }

  std::sort(factors.begin(), factors.end(),
            [](const auto &lhs, const auto &rhs) {
              if (detail::PolyCanonicalLess(lhs.first, rhs.first)) {
                return true;
              }
              if (detail::PolyCanonicalLess(rhs.first, lhs.first)) {
                return false;
              }
              return lhs.second < rhs.second;
            });

  return {lc, std::move(factors)};
}

// Evaluates poly at every raw point in eval_points and returns the indices i
// where poly(eval_points[i]) == 0.
template <typename GaloisField>
std::vector<size_t> ChienSearchIndices(
    const PolynomialDense<GaloisField> &poly,
    std::span<const typename GaloisField::element_type> eval_points) {
  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  std::vector<size_t> root_indices;
  if (eval_points.empty()) {
    return root_indices;
  }
  const auto field = poly.Field();
  const auto values = poly.EvalBatch(eval_points);
  for (size_t i = 0; i < values.size(); ++i) {
    if (field->IsZero(values[i])) {
      root_indices.push_back(i);
    }
  }
  return root_indices;
}

// Evaluates poly over raw points and returns all raw points that are roots.
template <typename GaloisField>
std::vector<typename GaloisField::element_type> ChienSearch(
    const PolynomialDense<GaloisField> &poly,
    std::span<const typename GaloisField::element_type> eval_points) {
  using RawType = typename GaloisField::element_type;
  const auto indices = ChienSearchIndices(poly, eval_points);
  std::vector<RawType> roots;
  roots.reserve(indices.size());
  for (size_t idx : indices) {
    roots.push_back(eval_points[idx]);
  }
  return roots;
}

// Evaluates poly over field elements and returns all elements in eval_points
// that are roots.
template <typename GaloisField>
std::vector<GaloisFieldElement<GaloisField>> ChienSearch(
    const PolynomialDense<GaloisField> &poly,
    const std::vector<GaloisFieldElement<GaloisField>> &eval_points) {
  using RawType = typename GaloisField::element_type;
  using ElementType = GaloisFieldElement<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  const auto field = poly.Field();
  std::vector<RawType> raw_pts;
  raw_pts.reserve(eval_points.size());
  for (const auto &pt : eval_points) {
    if (!pt.HasField() || pt.Field() != field) {
      throw std::invalid_argument(
          "All evaluation points must belong to the polynomial's field.");
    }
    raw_pts.push_back(pt.RawValue());
  }

  const auto indices =
      ChienSearchIndices(poly, std::span<const RawType>(raw_pts));
  std::vector<ElementType> roots;
  roots.reserve(indices.size());
  for (size_t idx : indices) {
    roots.push_back(eval_points[idx]);
  }
  return roots;
}

// Finds all distinct roots of a nonzero polynomial in F_q for ANY field size
// via Cantor-Zassenhaus splitting of gcd(f, x^q - x).
template <typename GaloisField>
std::vector<GaloisFieldElement<GaloisField>> FindRoots(
    const PolynomialDense<GaloisField> &poly,
    uint64_t seed = 0x9E3779B97F4A7C15ULL) {
  using ElementType = GaloisFieldElement<GaloisField>;
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  if (poly.IsZero()) {
    throw std::invalid_argument(
        "Zero polynomial has infinitely many roots; cannot enumerate.");
  }

  const auto field = poly.Field();
  const int deg = poly.Degree();
  if (deg == 0) {
    return {};
  }
  if (deg == 1) {
    const auto root_raw =
        field->Neg(field->Div(poly.RawCoeff(0), poly.RawCoeff(1)));
    return {ElementType::FromInternal(root_raw, field)};
  }

  const PolynomialType f = poly.Monic();
  const PolynomialType x = PolynomialType::Monomial(
      field, 1, field->MultiplicativeIdentity(), poly.GetVariable());
  const PolynomialType x_q = detail::PowQMod(x, f);
  const PolynomialType g = PolynomialGcd(f, x_q - x);

  if (g.Degree() <= 0) {
    return {};
  }
  if (g.Degree() == 1) {
    return {ElementType::FromInternal(field->Neg(g.RawCoeff(0)), field)};
  }

  std::mt19937_64 rng(seed);
  const auto linear_factors = EqualDegreeFactor(g, 1, rng);
  std::vector<ElementType> roots;
  roots.reserve(linear_factors.size());
  for (const auto &lin : linear_factors) {
    roots.push_back(
        ElementType::FromInternal(field->Neg(lin.RawCoeff(0)), field));
  }
  std::sort(roots.begin(), roots.end());
  return roots;
}

// O(k^2) Newton divided-difference interpolation over raw contiguous spans.
template <typename GaloisField>
PolynomialDense<GaloisField> InterpolateRaw(
    const std::shared_ptr<GaloisField> &field,
    std::span<const typename GaloisField::element_type> xs,
    std::span<const typename GaloisField::element_type> ys,
    const std::string &var = "x") {
  using RawType = typename GaloisField::element_type;

  if (!field) {
    throw std::invalid_argument("GaloisField pointer cannot be null.");
  }
  if (var.empty()) {
    throw std::invalid_argument("Variable name cannot be empty.");
  }
  if (xs.size() != ys.size()) {
    throw std::invalid_argument(
        "Interpolation requires xs and ys to have the same size.");
  }
  if (xs.empty()) {
    throw std::invalid_argument(
        "Interpolation requires at least one data point.");
  }

  const size_t k = xs.size();
  std::vector<RawType> dd(ys.begin(), ys.end());

  // Divided differences table; checks every distinct pair (i, i - j) for
  // duplicate x-coordinates.
  for (size_t j = 1; j < k; ++j) {
    for (size_t i = k - 1; i >= j; --i) {
      const RawType denom = field->Sub(xs[i], xs[i - j]);
      if (field->IsZero(denom)) {
        throw std::invalid_argument(
            "Interpolation x-coordinates must be distinct.");
      }
      const RawType num = field->Sub(dd[i], dd[i - 1]);
      dd[i] = field->Div(num, denom);
    }
  }

  // Expand Newton form into standard monomial basis in O(k^2).
  std::vector<RawType> coeffs(k, field->AdditiveIdentity());
  coeffs[0] = dd[k - 1];
  size_t cur_len = 1;
  for (size_t step = k - 1; step > 0; --step) {
    const size_t i = step - 1;
    const RawType neg_xi = field->Neg(xs[i]);
    coeffs[cur_len] = coeffs[cur_len - 1];
    for (size_t d = cur_len - 1; d > 0; --d) {
      coeffs[d] = field->Add(coeffs[d - 1], field->Mul(coeffs[d], neg_xi));
    }
    coeffs[0] = field->Add(field->Mul(coeffs[0], neg_xi), dd[i]);
    ++cur_len;
  }

  return PolynomialDense<GaloisField>::FromRaw(std::move(coeffs), field, var);
}

// O(k^2) polynomial interpolation returning the unique polynomial p of degree <
// k satisfying p(xs[i]) == ys[i] for all i.
template <typename GaloisField>
PolynomialDense<GaloisField> Interpolate(
    const std::vector<GaloisFieldElement<GaloisField>> &xs,
    const std::vector<GaloisFieldElement<GaloisField>> &ys,
    const std::string &var = "x") {
  using RawType = typename GaloisField::element_type;

  if (xs.size() != ys.size()) {
    throw std::invalid_argument(
        "Interpolation requires xs and ys to have the same size.");
  }
  if (xs.empty()) {
    throw std::invalid_argument(
        "Interpolation requires at least one data point.");
  }
  if (!xs[0].HasField()) {
    throw std::invalid_argument(
        "Interpolation points must have an associated Galois field.");
  }
  const auto field = xs[0].Field();
  std::vector<RawType> raw_x;
  std::vector<RawType> raw_y;
  raw_x.reserve(xs.size());
  raw_y.reserve(ys.size());
  for (size_t i = 0; i < xs.size(); ++i) {
    if (!xs[i].HasField() || xs[i].Field() != field || !ys[i].HasField() ||
        ys[i].Field() != field) {
      throw std::invalid_argument(
          "All interpolation points must belong to the same Galois field.");
    }
    raw_x.push_back(xs[i].RawValue());
    raw_y.push_back(ys[i].RawValue());
  }

  return InterpolateRaw(field, std::span<const RawType>(raw_x),
                        std::span<const RawType>(raw_y), var);
}

// Finds a monic irreducible polynomial of the requested degree >= 1 over field.
template <typename GaloisField>
PolynomialDense<GaloisField> FindIrreduciblePolynomial(
    const std::shared_ptr<GaloisField> &field, uint32_t degree,
    const std::string &var = "x", uint64_t seed = 0xD1B54A32D192ED03ULL) {
  using RawType = typename GaloisField::element_type;
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!field) {
    throw std::invalid_argument("GaloisField pointer cannot be null.");
  }
  if (degree == 0) {
    throw std::invalid_argument("Irreducible polynomial degree must be >= 1.");
  }
  if (var.empty()) {
    throw std::invalid_argument("Variable name cannot be empty.");
  }

  const RawType zero = field->AdditiveIdentity();
  const RawType one = field->MultiplicativeIdentity();
  if (degree == 1) {
    return PolynomialType::Monomial(field, 1, one, var);
  }

  // First try sparse monic trinomials x^n + x^k + 1.
  for (uint32_t k = 1; k < degree; ++k) {
    std::vector<RawType> coeffs(static_cast<size_t>(degree) + 1, zero);
    coeffs[0] = one;
    coeffs[k] = one;
    coeffs[degree] = one;
    PolynomialType cand =
        PolynomialType::FromRaw(std::move(coeffs), field, var);
    if (IsIrreducible(cand)) {
      return cand;
    }
  }

  std::mt19937_64 rng(seed ^
                      (static_cast<uint64_t>(degree) * 0x9E3779B97F4A7C15ULL));
  std::vector<RawType> coeffs(static_cast<size_t>(degree) + 1, zero);
  coeffs[degree] = one;
  while (true) {
    // Constant term must be nonzero for any irreducible polynomial of deg >= 2.
    RawType c0 = zero;
    while (field->IsZero(c0)) {
      c0 = field->Random(rng);
    }
    coeffs[0] = c0;
    for (uint32_t i = 1; i < degree; ++i) {
      coeffs[i] = field->Random(rng);
    }
    PolynomialType cand = PolynomialType::FromRaw(coeffs, field, var);
    if (IsIrreducible(cand)) {
      return cand;
    }
  }
}

// Finds a monic primitive polynomial of the requested degree >= 1 over field.
template <typename GaloisField>
PolynomialDense<GaloisField> FindPrimitivePolynomial(
    const std::shared_ptr<GaloisField> &field, uint32_t degree,
    const std::string &var = "x", uint64_t seed = 0x94D049BB133111EBULL) {
  using RawType = typename GaloisField::element_type;
  using PolynomialType = PolynomialDense<GaloisField>;

  if (!field) {
    throw std::invalid_argument("GaloisField pointer cannot be null.");
  }
  if (degree == 0) {
    throw std::invalid_argument("Primitive polynomial degree must be >= 1.");
  }
  if (var.empty()) {
    throw std::invalid_argument("Variable name cannot be empty.");
  }

  const RawType zero = field->AdditiveIdentity();
  const RawType one = field->MultiplicativeIdentity();
  if (degree == 1) {
    const RawType gen = field->MultiplicativeGenerator();
    return PolynomialType::FromRaw({field->Neg(gen), one}, field, var);
  }

  // Validate q^n - 1 fits in uint64_t up front.
  (void)detail::ComputeExtensionGroupOrder(*field, degree);

  // First try sparse monic trinomials x^n + x^k + 1.
  for (uint32_t k = 1; k < degree; ++k) {
    std::vector<RawType> coeffs(static_cast<size_t>(degree) + 1, zero);
    coeffs[0] = one;
    coeffs[k] = one;
    coeffs[degree] = one;
    PolynomialType cand =
        PolynomialType::FromRaw(std::move(coeffs), field, var);
    if (IsPrimitive(cand)) {
      return cand;
    }
  }

  std::mt19937_64 rng(seed ^
                      (static_cast<uint64_t>(degree) * 0xBF58476D1CE4E5B9ULL));
  std::vector<RawType> coeffs(static_cast<size_t>(degree) + 1, zero);
  coeffs[degree] = one;
  while (true) {
    RawType c0 = zero;
    while (field->IsZero(c0)) {
      c0 = field->Random(rng);
    }
    coeffs[0] = c0;
    for (uint32_t i = 1; i < degree; ++i) {
      coeffs[i] = field->Random(rng);
    }
    PolynomialType cand = PolynomialType::FromRaw(coeffs, field, var);
    if (IsPrimitive(cand)) {
      return cand;
    }
  }
}

template <typename GaloisField>
uint64_t BinaryPolynomialToUint(const PolynomialDense<GaloisField> &poly) {
  if (!poly.HasField()) {
    throw std::invalid_argument(
        "Polynomial must have an associated Galois field.");
  }
  const auto field = poly.Field();
  if (field->Characteristic() != 2) {
    throw std::invalid_argument(
        "BinaryPolynomialToUint requires a field of characteristic 2.");
  }
  if (poly.Degree() >= 64) {
    throw std::invalid_argument("Polynomial degree (" +
                                std::to_string(poly.Degree()) +
                                ") must be less than maximum degree (64).");
  }

  uint64_t result = 0;
  for (int i = 0; i <= poly.Degree(); ++i) {
    const auto raw = poly.RawCoeff(static_cast<size_t>(i));
    if (field->IsOne(raw)) {
      result |= (uint64_t{1} << i);
    } else if (!field->IsZero(raw)) {
      throw std::invalid_argument(
          "Polynomial coefficients must be 0 or 1 for binary fields.");
    }
  }
  return result;
}

}  // namespace utils
}  // namespace xg

#endif  // XGALOIS_UTILS_POLY_HPP_
