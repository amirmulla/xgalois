#ifndef XGALOIS_FIELD_GF_BINARY_HPP_
#define XGALOIS_FIELD_GF_BINARY_HPP_

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(__PCLMUL__) && defined(__x86_64__)
#include <wmmintrin.h>
#endif

#if defined(__SSSE3__)
#include <tmmintrin.h>
#endif

#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

namespace detail {

class BinaryStreamFlagsGuard {
 public:
  explicit BinaryStreamFlagsGuard(std::ostream &os)
      : os_(os), flags_(os.flags()) {}
  ~BinaryStreamFlagsGuard() { os_.flags(flags_); }

  BinaryStreamFlagsGuard(const BinaryStreamFlagsGuard &) = delete;
  BinaryStreamFlagsGuard &operator=(const BinaryStreamFlagsGuard &) = delete;

 private:
  std::ostream &os_;
  std::ios_base::fmtflags flags_;
};

inline std::string TrimAscii(std::string_view sv) {
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

// Carry-less multiplication of two 64-bit polynomials over GF(2).
// Returns {hi, lo} representing the 128-bit product.
inline std::pair<uint64_t, uint64_t> ClMul64(uint64_t a, uint64_t b) noexcept {
#if defined(__PCLMUL__) && defined(__x86_64__)
  const __m128i va = _mm_cvtsi64_si128(static_cast<long long>(a));
  const __m128i vb = _mm_cvtsi64_si128(static_cast<long long>(b));
  const __m128i prod = _mm_clmulepi64_si128(va, vb, 0x00);
  const uint64_t lo = static_cast<uint64_t>(_mm_cvtsi128_si64(prod));
  const uint64_t hi =
      static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(prod, 8)));
  return {hi, lo};
#else
  uint64_t lo = 0;
  uint64_t hi = 0;
  uint64_t tmp = b;
  while (tmp != 0) {
    const int bit = std::countr_zero(tmp);
    lo ^= (a << bit);
    if (bit > 0) {
      hi ^= (a >> (64 - bit));
    }
    tmp &= (tmp - 1);
  }
  return {hi, lo};
#endif
}

// Computes the low m bits of mu(x) = floor(x^(2m) / P(x)) for a monic
// degree-m binary polynomial P(x) = x^m + modulus_low (1 <= m <= 64).
inline uint64_t ComputeBinaryBarrettMuLow(uint8_t m,
                                          uint64_t modulus_low) noexcept {
  const uint64_t mask = (m < 64) ? ((uint64_t{1} << m) - 1ULL) : ~uint64_t{0};
  uint64_t rem = modulus_low & mask;
  uint64_t q = 0;
  for (int i = static_cast<int>(m) - 1; i >= 0; --i) {
    const bool top = ((rem >> (m - 1)) & 1ULL) != 0;
    rem = (rem << 1) & mask;
    if (top) {
      rem ^= (modulus_low & mask);
      q |= (uint64_t{1} << i);
    }
  }
  return q;
}

// Reduces a 128-bit binary polynomial (hi, lo) of degree <= 2*m - 2 modulo
// P(x) = x^m + modulus_low (1 <= m <= 64).
inline uint64_t ReduceBinaryPoly128(uint64_t hi, uint64_t lo, uint8_t m,
                                    uint64_t modulus_low) noexcept {
  while (hi != 0) {
    const int k = std::bit_width(hi) - 1;  // 0 <= k <= 62
    const int shift = 64 + k - static_cast<int>(m);
    hi ^= (uint64_t{1} << k);
    lo ^= (modulus_low << shift);
    if (shift > 0) {
      hi ^= (modulus_low >> (64 - shift));
    }
  }

  if (m < 64) {
    const uint64_t modulus = (uint64_t{1} << m) | modulus_low;
    while (true) {
      const int bw = std::bit_width(lo);
      if (bw <= static_cast<int>(m)) {
        break;
      }
      lo ^= (modulus << (bw - 1 - static_cast<int>(m)));
    }
  }
  return lo;
}

// Exact binary polynomial multiplication modulo P(x) = x^m + modulus_low for
// 1 <= m <= 64. Uses PCLMULQDQ + exact carry-less Barrett reduction when
// __PCLMUL__ is available, and shift-XOR fallback otherwise.
inline uint64_t MulBinaryMod64(uint64_t a, uint64_t b, uint8_t m,
                               uint64_t modulus_low,
                               uint64_t barrett_mu_low) noexcept {
  if (a == 0 || b == 0) {
    return 0;
  }
#if defined(__PCLMUL__) && defined(__x86_64__)
  const auto [c_hi_word, c_lo_word] = ClMul64(a, b);
  if (m == 64) {
    const auto [t1_hi, t1_lo] = ClMul64(c_hi_word, barrett_mu_low);
    (void)t1_lo;
    const uint64_t q = c_hi_word ^ t1_hi;
    const auto [t2_hi, t2_lo] = ClMul64(q, modulus_low);
    (void)t2_hi;
    return c_lo_word ^ t2_lo;
  }
  const uint64_t mask = (uint64_t{1} << m) - 1ULL;
  const uint64_t c_lo = c_lo_word & mask;
  const uint64_t c_hi = (c_lo_word >> m) | (c_hi_word << (64 - m));
  const auto [t1_hi, t1_lo] = ClMul64(c_hi, barrett_mu_low);
  const uint64_t q = c_hi ^ ((t1_lo >> m) | (t1_hi << (64 - m)));
  const auto [t2_hi, t2_lo] = ClMul64(q, modulus_low);
  (void)t2_hi;
  return (c_lo ^ t2_lo) & mask;
#else
  (void)barrett_mu_low;
  if (m <= 32) {
    uint64_t prod = 0;
    uint64_t tmp = b;
    while (tmp != 0) {
      const int bit = std::countr_zero(tmp);
      prod ^= (a << bit);
      tmp &= (tmp - 1);
    }
    const uint64_t modulus = (uint64_t{1} << m) | modulus_low;
    while (true) {
      const int bw = std::bit_width(prod);
      if (bw <= static_cast<int>(m)) {
        break;
      }
      prod ^= (modulus << (bw - 1 - static_cast<int>(m)));
    }
    return prod;
  }
  const auto [hi, lo] = ClMul64(a, b);
  return ReduceBinaryPoly128(hi, lo, m, modulus_low);
#endif
}

// Binary polynomial GCD of P(x) = x^m + modulus_low and g(x) (deg(g) < m).
inline uint64_t BinaryPolyGcdWithModulus(uint8_t m, uint64_t modulus_low,
                                         uint64_t g) noexcept {
  if (g == 0) {
    return (m < 64) ? ((uint64_t{1} << m) | modulus_low) : uint64_t{2};
  }
  if (g == 1) {
    return 1;
  }

  // First reduction step: r = P(x) mod g(x).
  const int deg_g = std::bit_width(g) - 1;  // 1 <= deg_g <= m - 1
  const int shift0 = static_cast<int>(m) - deg_g;
  const uint64_t g_no_msb = g ^ (uint64_t{1} << deg_g);
  uint64_t r = modulus_low ^ (g_no_msb << shift0);

  while (r != 0) {
    const int deg_r = std::bit_width(r) - 1;
    if (deg_r < deg_g) {
      break;
    }
    r ^= (g << (deg_r - deg_g));
  }

  uint64_t r0 = g;
  uint64_t r1 = r;
  while (r1 != 0) {
    const int deg_r1 = std::bit_width(r1) - 1;
    while (r0 != 0) {
      const int deg_r0 = std::bit_width(r0) - 1;
      if (deg_r0 < deg_r1) {
        break;
      }
      r0 ^= (r1 << (deg_r0 - deg_r1));
    }
    std::swap(r0, r1);
  }
  return r0;
}

// Rabin's exact irreducibility test for a monic binary polynomial
// P(x) = x^m + modulus_low of degree 1 <= m <= 64.
inline bool IsBinaryPolyIrreducible(uint8_t m, uint64_t modulus_low) {
  if (m == 0 || m > 64) {
    return false;
  }
  if (m == 1) {
    return true;
  }
  // Divisible by x?
  if ((modulus_low & 1ULL) == 0) {
    return false;
  }
  // Divisible by (x + 1)? P(1) == 0 iff popcount(P) is even, i.e.,
  // popcount(modulus_low) is odd (since x^m contributes 1).
  if ((std::popcount(modulus_low) & 1) != 0) {
    return false;
  }

  const uint64_t mu_low = ComputeBinaryBarrettMuLow(m, modulus_low);
  const auto prime_factors = utils::DistinctPrimeFactors(m);

  // Collect checkpoints k = m / q in ascending order.
  std::vector<uint8_t> checkpoints;
  checkpoints.reserve(prime_factors.size());
  for (uint64_t q : prime_factors) {
    checkpoints.push_back(static_cast<uint8_t>(m / q));
  }
  std::sort(checkpoints.begin(), checkpoints.end());

  uint64_t h = 2;  // x mod P(x)
  size_t cp_idx = 0;
  for (uint8_t step = 1; step <= m; ++step) {
    h = MulBinaryMod64(h, h, m, modulus_low, mu_low);  // h = x^(2^step) mod P
    while (cp_idx < checkpoints.size() && checkpoints[cp_idx] == step) {
      const uint64_t diff = h ^ 2ULL;  // x^(2^(m/q)) - x mod P
      if (diff == 0 || BinaryPolyGcdWithModulus(m, modulus_low, diff) != 1) {
        return false;
      }
      ++cp_idx;
    }
  }

  // Final condition: x^(2^m) == x (mod P(x)).
  return h == 2ULL;
}

struct ParsedBinaryPoly {
  int degree;         // -1 for zero polynomial, otherwise 0..64.
  uint64_t low_bits;  // Bits 0..63 (if degree == 64, bit 64 is implicit).
};

// Direct lightweight binary polynomial parser (no PolynomialDense dependency).
// Accepts terms separated by '+' or '-', with variable `variable_name` or
// 'x'/'X', e.g. "x^8 + x^4 + x^3 + x + 1", "α^2 + α + 1", "1", "0".
inline ParsedBinaryPoly ParseBinaryPolynomial(std::string_view input,
                                              std::string_view variable_name,
                                              int max_allowed_degree,
                                              bool allow_x = true) {
  std::string cleaned;
  cleaned.reserve(input.size());
  for (char ch : input) {
    if (std::isspace(static_cast<unsigned char>(ch)) == 0) {
      cleaned.push_back(ch);
    }
  }
  if (cleaned.empty()) {
    throw std::invalid_argument("Polynomial string cannot be empty");
  }

  std::vector<std::string_view> terms;
  size_t pos = 0;
  if (cleaned[0] == '+' || cleaned[0] == '-') {
    pos = 1;
    if (pos >= cleaned.size()) {
      throw std::invalid_argument("Invalid polynomial string: " +
                                  std::string(input));
    }
  }

  size_t term_start = pos;
  for (size_t i = pos; i < cleaned.size(); ++i) {
    const char c = cleaned[i];
    if (c == '+' || c == '-') {
      if (i > 0 && cleaned[i - 1] == '^') {
        throw std::invalid_argument("Negative exponent in polynomial: " +
                                    std::string(input));
      }
      if (i == term_start) {
        throw std::invalid_argument("Empty term in polynomial: " +
                                    std::string(input));
      }
      terms.push_back(
          std::string_view(cleaned).substr(term_start, i - term_start));
      term_start = i + 1;
    }
  }
  if (term_start >= cleaned.size()) {
    throw std::invalid_argument("Trailing operator in polynomial: " +
                                std::string(input));
  }
  terms.push_back(std::string_view(cleaned).substr(term_start));

  bool bit64 = false;
  uint64_t low_bits = 0;

  for (std::string_view term : terms) {
    size_t var_pos = std::string_view::npos;
    size_t var_len = 0;

    if (!variable_name.empty()) {
      var_pos = term.find(variable_name);
      if (var_pos != std::string_view::npos) {
        var_len = variable_name.size();
      }
    }
    if (var_pos == std::string_view::npos && allow_x) {
      var_pos = term.find('x');
      if (var_pos == std::string_view::npos) {
        var_pos = term.find('X');
      }
      if (var_pos != std::string_view::npos) {
        var_len = 1;
      }
    }

    uint64_t coeff_parity = 1;
    uint64_t deg = 0;

    if (var_pos == std::string_view::npos) {
      // Pure constant term.
      for (char ch : term) {
        if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
          throw std::invalid_argument("Invalid term in binary polynomial: " +
                                      std::string(term));
        }
      }
      try {
        size_t parsed = 0;
        const uint64_t val = std::stoull(std::string(term), &parsed, 10);
        if (parsed != term.size()) {
          throw std::invalid_argument("Invalid constant term: " +
                                      std::string(term));
        }
        coeff_parity = val & 1ULL;
        deg = 0;
      } catch (const std::invalid_argument &) {
        throw;
      } catch (const std::out_of_range &) {
        throw std::invalid_argument("Constant out of range in polynomial: " +
                                    std::string(term));
      }
    } else {
      std::string_view prefix = term.substr(0, var_pos);
      if (!prefix.empty()) {
        if (prefix.back() == '*') {
          prefix.remove_suffix(1);
          if (prefix.empty()) {
            throw std::invalid_argument("Invalid coefficient in term: " +
                                        std::string(term));
          }
        }
        for (char ch : prefix) {
          if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
            throw std::invalid_argument("Invalid coefficient in term: " +
                                        std::string(term));
          }
        }
        try {
          size_t parsed = 0;
          const uint64_t val = std::stoull(std::string(prefix), &parsed, 10);
          if (parsed != prefix.size()) {
            throw std::invalid_argument("Invalid coefficient in term: " +
                                        std::string(term));
          }
          coeff_parity = val & 1ULL;
        } catch (const std::invalid_argument &) {
          throw;
        } catch (const std::out_of_range &) {
          throw std::invalid_argument("Coefficient out of range in term: " +
                                      std::string(term));
        }
      }

      const std::string_view suffix = term.substr(var_pos + var_len);
      if (suffix.empty()) {
        deg = 1;
      } else if (suffix[0] == '^') {
        const std::string_view exp_sv = suffix.substr(1);
        if (exp_sv.empty()) {
          throw std::invalid_argument("Missing exponent in term: " +
                                      std::string(term));
        }
        for (char ch : exp_sv) {
          if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
            throw std::invalid_argument("Invalid exponent in term: " +
                                        std::string(term));
          }
        }
        try {
          size_t parsed = 0;
          deg = std::stoull(std::string(exp_sv), &parsed, 10);
          if (parsed != exp_sv.size()) {
            throw std::invalid_argument("Invalid exponent in term: " +
                                        std::string(term));
          }
        } catch (const std::invalid_argument &) {
          throw;
        } catch (const std::out_of_range &) {
          throw std::invalid_argument("Exponent out of range in term: " +
                                      std::string(term));
        }
      } else {
        throw std::invalid_argument("Unexpected characters in term: " +
                                    std::string(term));
      }
    }

    if (coeff_parity == 1) {
      if (deg > static_cast<uint64_t>(max_allowed_degree)) {
        throw std::invalid_argument("Polynomial degree (" +
                                    std::to_string(deg) +
                                    ") exceeds maximum allowed degree (" +
                                    std::to_string(max_allowed_degree) + ")");
      }
      if (deg == 64) {
        bit64 = !bit64;
      } else {
        low_bits ^= (uint64_t{1} << deg);
      }
    }
  }

  int actual_degree = -1;
  if (bit64) {
    actual_degree = 64;
  } else if (low_bits != 0) {
    actual_degree = std::bit_width(low_bits) - 1;
  }
  return {actual_degree, low_bits};
}

}  // namespace detail

class GaloisFieldBinary final : public GaloisFieldBase<uint8_t> {
 public:
  explicit GaloisFieldBinary(const std::string &rep = "int")
      : representation_(ValidateRep(utils::ConvertRepresentation(rep))) {}

  GaloisFieldBinary(const GaloisFieldBinary &other)
      : representation_(other.representation_.load(std::memory_order_relaxed)) {
  }

  GaloisFieldBinary &operator=(const GaloisFieldBinary &other) {
    if (this != &other) {
      representation_.store(
          other.representation_.load(std::memory_order_relaxed),
          std::memory_order_relaxed);
    }
    return *this;
  }

  uint64_t Characteristic() const override { return 2; }
  uint64_t Order() const override { return 2; }
  uint32_t Degree() const override { return 1; }
  uint64_t Modulus() const { return 2; }

  inline uint8_t Add(const uint8_t &a, const uint8_t &b) const override {
    return static_cast<uint8_t>((a ^ b) & 1U);
  }

  inline uint8_t Sub(const uint8_t &a, const uint8_t &b) const override {
    return static_cast<uint8_t>((a ^ b) & 1U);
  }

  inline uint8_t Mul(const uint8_t &a, const uint8_t &b) const override {
    return static_cast<uint8_t>((a & b) & 1U);
  }

  inline uint8_t Div(const uint8_t &a, const uint8_t &b) const override {
    if ((b & 1U) == 0) {
      throw std::domain_error("Division by zero in GF(2)");
    }
    return static_cast<uint8_t>(a & 1U);
  }

  inline uint8_t Neg(const uint8_t &a) const override {
    return static_cast<uint8_t>(a & 1U);
  }

  inline uint8_t Inv(const uint8_t &a) const override {
    if ((a & 1U) == 0) {
      throw std::domain_error("Inverse of zero in GF(2)");
    }
    return 1;
  }

  inline uint8_t Pow(const uint8_t &a, uint64_t exp) const override {
    if ((a & 1U) == 0) {
      return (exp == 0) ? uint8_t{1} : uint8_t{0};
    }
    return 1;
  }

  inline uint8_t Sqrt(const uint8_t &a) const override {
    return static_cast<uint8_t>(a & 1U);
  }

  uint64_t Log(const uint8_t &a) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined in GF(2)");
    }
    if (a == 1) {
      return 0;
    }
    throw std::invalid_argument("Log defined only for 1 in GF(2)");
  }

  uint64_t Log(const uint8_t &a, const uint8_t &generator) const override {
    if (generator != 1) {
      throw std::invalid_argument("Generator in GF(2) must be 1");
    }
    return Log(a);
  }

  uint8_t Random() const override {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    return Random(gen);
  }

  uint8_t Random(std::mt19937_64 &rng) const override {
    std::uniform_int_distribution<unsigned int> distrib(0U, 1U);
    return static_cast<uint8_t>(distrib(rng));
  }

  uint8_t MultiplicativeIdentity() const override { return 1; }
  uint8_t AdditiveIdentity() const override { return 0; }

  uint8_t GetElementValue(const uint8_t &value) const override {
    return static_cast<uint8_t>(value & 1U);
  }

  uint8_t SetElementValue(const uint8_t &value) const override {
    return static_cast<uint8_t>(value & 1U);
  }

  uint8_t SetElementValue(const std::string & /*value_str*/) const override {
    throw std::invalid_argument(
        "GF(2) does not support string representation for elements.");
  }

  uint8_t MultiplicativeGenerator() const override { return 1; }
  std::vector<uint8_t> MultiplicativeGenerators() const override { return {1}; }

  void Print(std::ostream &os) const override { os << "GF(2)"; }

  void Print(const uint8_t &a, std::ostream &os) const override {
    detail::BinaryStreamFlagsGuard guard(os);
    const FieldRepresentation rep = GetRepresentation();
    const uint64_t val = static_cast<uint64_t>(a & 1U);
    if (rep == FieldRepresentation::INT) {
      os << val;
    } else if (rep == FieldRepresentation::HEX) {
      os << "0x" << std::hex << val;
    } else {
      throw std::invalid_argument("Unsupported representation for GF(2)");
    }
  }

  std::string ToString(const uint8_t &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_.load(std::memory_order_relaxed);
  }

  void SetRepresentation(FieldRepresentation rep) override {
    representation_.store(ValidateRep(rep), std::memory_order_relaxed);
  }

  size_t MemoryBytes() const override { return sizeof(*this); }

  void AddArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                std::span<const uint8_t> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("AddArray: span sizes must match");
    }
    uint8_t *out = dst.data();
    const uint8_t *pa = a.data();
    const uint8_t *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = static_cast<uint8_t>((pa[i] ^ pb[i]) & 1U);
    }
  }

  void SubArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                std::span<const uint8_t> b) const override {
    AddArray(dst, a, b);
  }

  void MulArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                std::span<const uint8_t> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match");
    }
    uint8_t *out = dst.data();
    const uint8_t *pa = a.data();
    const uint8_t *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = static_cast<uint8_t>((pa[i] & pb[i]) & 1U);
    }
  }

  void ScalarMulArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                      const uint8_t &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match");
    }
    if ((scalar & 1U) == 0) {
      std::fill(dst.begin(), dst.end(), uint8_t{0});
      return;
    }
    uint8_t *out = dst.data();
    const uint8_t *pa = a.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = static_cast<uint8_t>(pa[i] & 1U);
    }
  }

  void FmaArray(std::span<uint8_t> dst, std::span<const uint8_t> a,
                const uint8_t &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match");
    }
    if ((scalar & 1U) == 0) {
      return;
    }
    uint8_t *out = dst.data();
    const uint8_t *pa = a.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = static_cast<uint8_t>((out[i] ^ pa[i]) & 1U);
    }
  }

  uint8_t DotProduct(std::span<const uint8_t> a,
                     std::span<const uint8_t> b) const override {
    const size_t n = a.size();
    if (b.size() != n) {
      throw std::invalid_argument("DotProduct: span sizes must match");
    }
    uint8_t acc = 0;
    const uint8_t *pa = a.data();
    const uint8_t *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      acc = static_cast<uint8_t>(acc ^ (pa[i] & pb[i]));
    }
    return static_cast<uint8_t>(acc & 1U);
  }

 private:
  static FieldRepresentation ValidateRep(FieldRepresentation rep) {
    if (rep != FieldRepresentation::INT && rep != FieldRepresentation::HEX) {
      throw std::invalid_argument(
          "GF(2) supports only INT and HEX representations for its elements.");
    }
    return rep;
  }

  std::atomic<FieldRepresentation> representation_;
};

template <UnsignedFieldWord ElementType = uint32_t>
class GaloisFieldBinaryExtension : public GaloisFieldBase<ElementType> {
  static_assert(sizeof(ElementType) <= sizeof(uint64_t),
                "ElementType must be at most 64 bits");

 public:
  explicit GaloisFieldBinaryExtension(uint8_t m, const std::string &rep = "int",
                                      const std::string &irreducible_poly = "",
                                      const std::string &variable_name = "α",
                                      bool check_irreducible = false,
                                      const std::string &generator_name = "g")
      : m_(m),
        modulus_low_(0),
        irreducible_poly_(0),
        barrett_mu_low_(0),
        order_(0),
        group_order_(0),
        is_conway_poly_(false),
        cached_generator_(0),
        representation_(ValidateRep(utils::ConvertRepresentation(rep))),
        variable_name_(variable_name),
        generator_name_(generator_name) {
    constexpr uint32_t kMaxDegree = sizeof(ElementType) * 8;
    if (m_ < 1 || static_cast<uint32_t>(m_) > kMaxDegree) {
      throw std::invalid_argument("Degree must be between 1 and " +
                                  std::to_string(kMaxDegree));
    }

    if (irreducible_poly.empty()) {
      bool loaded = false;
      try {
        databases::ConwayPolyDatabase conway_db;
        const auto poly_res = conway_db.fetch(2, static_cast<int>(m_));
        LoadPolyResult(poly_res);
        is_conway_poly_ = true;
        loaded = true;
      } catch (const std::runtime_error &) {
        loaded = false;
      }

      if (!loaded) {
        try {
          databases::IrreduciblePolyDatabase irred_db;
          const auto poly_res = irred_db.fetch(2, static_cast<int>(m_));
          LoadPolyResult(poly_res);
          loaded = true;
        } catch (const std::runtime_error &) {
          throw std::invalid_argument(
              "No suitable polynomial found in databases for GF(2^" +
              std::to_string(m_) +
              "). Please provide irreducible_poly manually.");
        }
      }
    } else {
      const auto parsed = detail::ParseBinaryPolynomial(
          irreducible_poly, variable_name_, 64, /*allow_x=*/true);
      if (parsed.degree != static_cast<int>(m_)) {
        throw std::invalid_argument(
            "Irreducible polynomial degree (" + std::to_string(parsed.degree) +
            ") must equal field extension degree (" + std::to_string(m_) + ")");
      }
      if (m_ == 64) {
        modulus_low_ = parsed.low_bits;
        irreducible_poly_ = modulus_low_;
      } else {
        irreducible_poly_ = parsed.low_bits;
        modulus_low_ = parsed.low_bits ^ (uint64_t{1} << m_);
      }
    }

    if (m_ < 64 && irreducible_poly_ <= 1) {
      throw std::invalid_argument("Invalid irreducible polynomial");
    }

    if (check_irreducible && !is_conway_poly_ &&
        !detail::IsBinaryPolyIrreducible(m_, modulus_low_)) {
      throw std::invalid_argument("Provided polynomial is not irreducible");
    }

    barrett_mu_low_ = detail::ComputeBinaryBarrettMuLow(m_, modulus_low_);
    order_ = (m_ < 64) ? (uint64_t{1} << m_) : uint64_t{0};
    group_order_ = (m_ < 64) ? ((uint64_t{1} << m_) - 1ULL) : ~uint64_t{0};

    if (m_ == 1) {
      cached_generator_ = ElementType{1};
    } else if (is_conway_poly_) {
      // Conway polynomials are primitive by definition: x (== 2) generates
      // GF(2^m)^* for all m >= 2.
      cached_generator_ = ElementType{2};
    }
  }

  GaloisFieldBinaryExtension(const GaloisFieldBinaryExtension &other)
      : m_(other.m_),
        modulus_low_(other.modulus_low_),
        irreducible_poly_(other.irreducible_poly_),
        barrett_mu_low_(other.barrett_mu_low_),
        order_(other.order_),
        group_order_(other.group_order_),
        is_conway_poly_(other.is_conway_poly_),
        cached_generator_(other.cached_generator_),
        representation_(other.representation_.load(std::memory_order_relaxed)),
        variable_name_(other.variable_name_),
        generator_name_(other.generator_name_) {}

  GaloisFieldBinaryExtension &operator=(
      const GaloisFieldBinaryExtension &other) {
    if (this != &other) {
      m_ = other.m_;
      modulus_low_ = other.modulus_low_;
      irreducible_poly_ = other.irreducible_poly_;
      barrett_mu_low_ = other.barrett_mu_low_;
      order_ = other.order_;
      group_order_ = other.group_order_;
      is_conway_poly_ = other.is_conway_poly_;
      cached_generator_ = other.cached_generator_;
      representation_.store(
          other.representation_.load(std::memory_order_relaxed),
          std::memory_order_relaxed);
      variable_name_ = other.variable_name_;
      generator_name_ = other.generator_name_;
    }
    return *this;
  }

  uint64_t Characteristic() const override { return 2; }
  uint64_t Order() const override { return order_; }
  uint32_t Degree() const override { return static_cast<uint32_t>(m_); }

  // For m < 64, returns full polynomial with bit m set.
  // For m == 64, returns the low 64 bits representing x^64 + Modulus().
  uint64_t Modulus() const { return irreducible_poly_; }
  uint64_t ModulusLow() const { return modulus_low_; }
  uint64_t GroupOrder() const { return group_order_; }

  inline ElementType Add(const ElementType &a,
                         const ElementType &b) const override {
    return static_cast<ElementType>(a ^ b);
  }

  inline ElementType Sub(const ElementType &a,
                         const ElementType &b) const override {
    return static_cast<ElementType>(a ^ b);
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    return MulFast(a, b);
  }

  inline ElementType Neg(const ElementType &a) const override { return a; }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }
    return ExtendedEuclideanAlgorithm(a);
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }
    return Mul(a, Inv(b));
  }

  inline ElementType Pow(const ElementType &a, uint64_t exp) const override {
    return PowExternal(a, exp);
  }

  inline ElementType Sqrt(const ElementType &a) const override {
    if (a == 0 || a == 1 || m_ == 1) {
      return a;
    }
    // In GF(2^m), sqrt(a) = a^(2^(m - 1)) = (m - 1) repeated squarings.
    ElementType res = a;
    for (uint8_t i = 0; i < m_ - 1; ++i) {
      res = MulFast(res, res);
    }
    return res;
  }

  uint64_t Log(const ElementType &a) const override {
    return Log(a, MultiplicativeGenerator());
  }

  uint64_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }
    if (generator == 0) {
      throw std::invalid_argument("Generator cannot be zero");
    }
    if (m_ == 1) {
      if (a != 1) {
        throw std::invalid_argument("Element out of range for GF(2^1)");
      }
      if (generator != 1) {
        throw std::invalid_argument("Generator in GF(2^1) must be 1");
      }
      return 0;
    }
    if (a == 1) {
      return 0;
    }
    if (generator == 1) {
      throw std::runtime_error("Element not in the group or invalid generator");
    }

    const uint64_t ord = OrderExternal(generator);
    if (PowExternal(a, ord) != 1) {
      throw std::runtime_error("Element not in the group or invalid generator");
    }

    if (ord <= (uint64_t{1} << 26)) {
      return DiscreteLogBsgs(a, generator, ord);
    }
    return DiscreteLogPohligHellman(a, generator, ord);
  }

  ElementType Random() const override {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    return Random(gen);
  }

  ElementType Random(std::mt19937_64 &rng) const override {
    std::uniform_int_distribution<uint64_t> distrib(0ULL, group_order_);
    return static_cast<ElementType>(distrib(rng));
  }

  ElementType MultiplicativeIdentity() const override { return ElementType{1}; }
  ElementType AdditiveIdentity() const override { return ElementType{0}; }

  ElementType GetElementValue(const ElementType &value) const override {
    return value;
  }

  ElementType SetElementValue(const ElementType &value) const override {
    if (static_cast<uint64_t>(value) > group_order_) {
      throw std::out_of_range("Element value exceeds field order");
    }
    return value;
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    return ParseElementStringExternal(value_str);
  }

  ElementType MultiplicativeGenerator() const override {
    return MultiplicativeGeneratorExternal();
  }

  std::vector<ElementType> MultiplicativeGenerators() const override {
    if (m_ == 1) {
      return {ElementType{1}};
    }

    const auto prime_factors = utils::DistinctPrimeFactors(group_order_);
    std::vector<ElementType> generators;
    for (uint64_t candidate = 2;; ++candidate) {
      const auto elem = static_cast<ElementType>(candidate);
      if (HasFullOrderWithFactors(elem, prime_factors)) {
        generators.push_back(elem);
      }
      if (candidate == group_order_) {
        break;
      }
    }
    return generators;
  }

  void Print(std::ostream &os) const override {
    os << "GF(2^" << static_cast<uint64_t>(m_) << ")";
  }

  void Print(const ElementType &a, std::ostream &os) const override {
    detail::BinaryStreamFlagsGuard guard(os);
    const FieldRepresentation rep = GetRepresentation();
    if (rep == FieldRepresentation::INT) {
      os << static_cast<uint64_t>(a);
    } else if (rep == FieldRepresentation::HEX) {
      os << "0x" << std::hex << static_cast<uint64_t>(a);
    } else if (rep == FieldRepresentation::POLY) {
      os << BinaryPolyToString(a);
    } else if (rep == FieldRepresentation::POW) {
      if (a == 0) {
        os << "0";
      } else {
        os << generator_name_ << "^" << Log(a);
      }
    } else if (rep == FieldRepresentation::LOG) {
      if (a == 0) {
        os << "-inf";
      } else {
        os << Log(a);
      }
    } else {
      throw std::invalid_argument("Unsupported representation for GF(2^m)");
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
    representation_.store(ValidateRep(rep), std::memory_order_relaxed);
  }

  size_t MemoryBytes() const override {
    return sizeof(*this) + variable_name_.capacity() +
           generator_name_.capacity();
  }

  void AddArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("AddArray: span sizes must match");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = static_cast<ElementType>(pa[i] ^ pb[i]);
    }
  }

  void SubArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    AddArray(dst, a, b);
  }

  void MulArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = MulFast(pa[i], pb[i]);
    }
  }

  void ScalarMulArray(std::span<ElementType> dst,
                      std::span<const ElementType> a,
                      const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match");
    }
    if (scalar == 0) {
      std::fill(dst.begin(), dst.end(), ElementType{0});
      return;
    }
    if (scalar == 1) {
      if (dst.data() != a.data()) {
        std::copy(a.begin(), a.end(), dst.begin());
      }
      return;
    }

    if constexpr (std::is_same_v<ElementType, uint8_t>) {
      ScalarMulNibbleU8(dst, a, scalar);
    } else if constexpr (std::is_same_v<ElementType, uint16_t>) {
      if (n >= 64) {
        ScalarMulNibbleU16(dst, a, scalar);
      } else {
        ElementType *out = dst.data();
        const ElementType *pa = a.data();
        for (size_t i = 0; i < n; ++i) {
          out[i] = MulFast(pa[i], scalar);
        }
      }
    } else {
      ElementType *out = dst.data();
      const ElementType *pa = a.data();
      for (size_t i = 0; i < n; ++i) {
        out[i] = MulFast(pa[i], scalar);
      }
    }
  }

  void FmaArray(std::span<ElementType> dst, std::span<const ElementType> a,
                const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match");
    }
    if (scalar == 0) {
      return;
    }
    if (scalar == 1) {
      ElementType *out = dst.data();
      const ElementType *pa = a.data();
      for (size_t i = 0; i < n; ++i) {
        out[i] = static_cast<ElementType>(out[i] ^ pa[i]);
      }
      return;
    }

    if constexpr (std::is_same_v<ElementType, uint8_t>) {
      FmaNibbleU8(dst, a, scalar);
    } else if constexpr (std::is_same_v<ElementType, uint16_t>) {
      if (n >= 64) {
        FmaNibbleU16(dst, a, scalar);
      } else {
        ElementType *out = dst.data();
        const ElementType *pa = a.data();
        for (size_t i = 0; i < n; ++i) {
          out[i] = static_cast<ElementType>(out[i] ^ MulFast(pa[i], scalar));
        }
      }
    } else {
      ElementType *out = dst.data();
      const ElementType *pa = a.data();
      for (size_t i = 0; i < n; ++i) {
        out[i] = static_cast<ElementType>(out[i] ^ MulFast(pa[i], scalar));
      }
    }
  }

  ElementType DotProduct(std::span<const ElementType> a,
                         std::span<const ElementType> b) const override {
    const size_t n = a.size();
    if (b.size() != n) {
      throw std::invalid_argument("DotProduct: span sizes must match");
    }
    ElementType acc = 0;
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      acc = static_cast<ElementType>(acc ^ MulFast(pa[i], pb[i]));
    }
    return acc;
  }

 protected:
  uint8_t m_;
  uint64_t modulus_low_;
  uint64_t irreducible_poly_;
  uint64_t barrett_mu_low_;
  uint64_t order_;
  uint64_t group_order_;
  bool is_conway_poly_;
  ElementType cached_generator_;
  std::atomic<FieldRepresentation> representation_;
  std::string variable_name_;
  std::string generator_name_;

  static FieldRepresentation ValidateRep(FieldRepresentation rep) {
    if (rep != FieldRepresentation::INT && rep != FieldRepresentation::HEX &&
        rep != FieldRepresentation::POLY && rep != FieldRepresentation::POW &&
        rep != FieldRepresentation::LOG) {
      throw std::invalid_argument(
          "GF(2^m) supports INT, HEX, POLY, POW, or LOG "
          "representation for its elements.");
    }
    return rep;
  }

  void LoadPolyResult(const databases::IrreduciblePolyResult &poly_res) {
    bool bit64 = false;
    uint64_t low_bits = 0;
    for (size_t i = 0; i < poly_res.nonzero_degrees.size(); ++i) {
      if ((poly_res.nonzero_coeffs[i] & 1) == 0) {
        continue;
      }
      const int deg = poly_res.nonzero_degrees[i];
      if (deg < 0 || deg > 64) {
        throw std::invalid_argument("Invalid polynomial degree from database");
      }
      if (deg == 64) {
        bit64 = !bit64;
      } else {
        low_bits ^= (uint64_t{1} << deg);
      }
    }

    const int actual_deg =
        bit64 ? 64 : (low_bits != 0 ? (std::bit_width(low_bits) - 1) : -1);
    if (actual_deg != static_cast<int>(m_)) {
      throw std::invalid_argument("Database polynomial degree mismatch");
    }

    if (m_ == 64) {
      modulus_low_ = low_bits;
      irreducible_poly_ = modulus_low_;
    } else {
      irreducible_poly_ = low_bits;
      modulus_low_ = low_bits ^ (uint64_t{1} << m_);
    }
  }

  inline ElementType MulFast(ElementType a, ElementType b) const noexcept {
    return static_cast<ElementType>(detail::MulBinaryMod64(
        static_cast<uint64_t>(a), static_cast<uint64_t>(b), m_, modulus_low_,
        barrett_mu_low_));
  }

  inline ElementType PowExternal(ElementType a, uint64_t exp) const noexcept {
    if (a == 0) {
      return (exp == 0) ? ElementType{1} : ElementType{0};
    }
    if (a == 1 || exp == 0) {
      return ElementType{1};
    }
    exp %= group_order_;
    ElementType result = 1;
    ElementType base = a;
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        result = MulFast(result, base);
      }
      exp >>= 1;
      if (exp == 0) {
        break;
      }
      base = MulFast(base, base);
    }
    return result;
  }

  // Shift-XOR binary polynomial Extended Euclidean Algorithm using bit_width.
  // Computes a^(-1) mod P(x) for all 1 <= m <= 64 with zero overflow and no
  // redundant s0/s1 tracking.
  ElementType ExtendedEuclideanAlgorithm(ElementType element) const {
    if (element == 1) {
      return ElementType{1};
    }

    const uint64_t a = static_cast<uint64_t>(element);
    const int deg_a = std::bit_width(a) - 1;  // 1 <= deg_a <= m_ - 1
    const int shift0 = static_cast<int>(m_) - deg_a;
    const uint64_t a_no_msb = a ^ (uint64_t{1} << deg_a);

    // First division step: P(x) = q * a + r, where deg(r) < deg(a) < m_ <= 64.
    uint64_t r = modulus_low_ ^ (a_no_msb << shift0);
    uint64_t q = uint64_t{1} << shift0;

    while (r != 0) {
      const int deg_r = std::bit_width(r) - 1;
      if (deg_r < deg_a) {
        break;
      }
      const int s = deg_r - deg_a;
      r ^= (a << s);
      q ^= (uint64_t{1} << s);
    }

    uint64_t r0 = a;
    uint64_t r1 = r;
    uint64_t t0 = 1;
    uint64_t t1 = q;

    // Loop invariant: t0 * a == r0 (mod P), t1 * a == r1 (mod P),
    // and deg(t1) + deg(r0) == m_. Stopping when r1 == 1 ensures deg(r1) >= 1
    // during all shifts, so deg(t1) + s <= m_ - 1 <= 63 (zero 64-bit overflow).
    while (r1 > 1) {
      const int deg_r1 = std::bit_width(r1) - 1;
      while (r0 != 0) {
        const int deg_r0 = std::bit_width(r0) - 1;
        if (deg_r0 < deg_r1) {
          break;
        }
        const int s = deg_r0 - deg_r1;
        r0 ^= (r1 << s);
        t0 ^= (t1 << s);
      }
      std::swap(r0, r1);
      std::swap(t0, t1);
    }

    if (r1 == 1) {
      return static_cast<ElementType>(t1);
    }
    throw std::domain_error(
        "Element is not invertible (modulus may be reducible)");
  }

  bool HasFullOrderWithFactors(
      ElementType element, const std::vector<uint64_t> &prime_factors) const {
    if (element == 0) {
      return false;
    }
    if (m_ == 1) {
      return element == 1;
    }
    if (element == 1) {
      return false;
    }
    if (MulFast(PowExternal(element, group_order_ - 1), element) != 1) {
      return false;
    }
    for (uint64_t factor : prime_factors) {
      if (PowExternal(element, group_order_ / factor) == 1) {
        return false;
      }
    }
    return true;
  }

  bool HasFullOrder(ElementType element) const {
    if (element == 0) {
      return false;
    }
    if (m_ == 1) {
      return element == 1;
    }
    if (element == 1) {
      return false;
    }
    return HasFullOrderWithFactors(element,
                                   utils::DistinctPrimeFactors(group_order_));
  }

  ElementType MultiplicativeGeneratorExternal() const {
    if (cached_generator_ != 0) {
      return cached_generator_;
    }
    if (m_ == 1) {
      return ElementType{1};
    }

    const auto prime_factors = utils::DistinctPrimeFactors(group_order_);
    for (uint64_t candidate = 2;; ++candidate) {
      const auto elem = static_cast<ElementType>(candidate);
      if (HasFullOrderWithFactors(elem, prime_factors)) {
        return elem;
      }
      if (candidate == group_order_) {
        break;
      }
    }
    throw std::runtime_error("Could not find a primitive element");
  }

  uint64_t OrderExternal(ElementType a) const {
    if (a == 1) {
      return 1;
    }
    uint64_t ord = group_order_;
    const auto prime_factors = utils::DistinctPrimeFactors(group_order_);
    for (uint64_t factor : prime_factors) {
      while (ord % factor == 0 && PowExternal(a, ord / factor) == 1) {
        ord /= factor;
      }
    }
    return ord;
  }

  // Baby-Step Giant-Step discrete logarithm in a cyclic subgroup of order
  // `sub_order` generated by `gen`.
  uint64_t DiscreteLogBsgs(ElementType target, ElementType gen,
                           uint64_t sub_order) const {
    if (target == 1) {
      return 0;
    }
    const uint64_t step = utils::ISqrt(sub_order) + 1;
    constexpr uint64_t kMaxBsgsSteps = uint64_t{1} << 22;
    if (step > kMaxBsgsSteps) {
      throw std::runtime_error(
          "Discrete logarithm subgroup order too large for BSGS");
    }
    std::unordered_map<uint64_t, uint64_t> baby_steps;
    baby_steps.reserve(static_cast<size_t>(step));

    ElementType cur = 1;
    for (uint64_t j = 0; j < step; ++j) {
      baby_steps.emplace(static_cast<uint64_t>(cur), j);
      cur = MulFast(cur, gen);
    }

    const ElementType giant_factor =
        ExtendedEuclideanAlgorithm(PowExternal(gen, step));
    ElementType gamma = target;
    for (uint64_t i = 0; i <= step; ++i) {
      const auto it = baby_steps.find(static_cast<uint64_t>(gamma));
      if (it != baby_steps.end()) {
        const uint64_t ans = i * step + it->second;
        if (ans < sub_order) {
          return ans;
        }
      }
      gamma = MulFast(gamma, giant_factor);
    }
    throw std::runtime_error("Element not in the group or invalid generator");
  }

  // Pohlig-Hellman + BSGS discrete logarithm for large binary fields
  // (e.g., GF(2^32), GF(2^64)).
  uint64_t DiscreteLogPohligHellman(ElementType target, ElementType gen,
                                    uint64_t sub_order) const {
    const auto factor_pairs = utils::PrimeFactorizeWithMultiplicity(sub_order);
    const ElementType gen_inv = ExtendedEuclideanAlgorithm(gen);

    uint64_t combined_rem = 0;
    uint64_t combined_mod = 1;

    for (const auto &[prime, mult] : factor_pairs) {
      uint64_t pe = 1;
      for (uint64_t e = 0; e < mult; ++e) {
        pe *= prime;
      }

      const ElementType g0 = PowExternal(gen, sub_order / prime);
      uint64_t x_pe = 0;
      uint64_t p_pow = 1;
      ElementType a_cur = target;

      for (uint64_t k = 0; k < mult; ++k) {
        const uint64_t exp = sub_order / (p_pow * prime);
        const ElementType h_k = PowExternal(a_cur, exp);
        const uint64_t d_k = DiscreteLogBsgs(h_k, g0, prime);
        if (d_k != 0) {
          x_pe += d_k * p_pow;
          a_cur = MulFast(a_cur, PowExternal(gen_inv, d_k * p_pow));
        }
        p_pow *= prime;
      }

      // CRT combine: x == combined_rem (mod combined_mod), x == x_pe (mod pe).
      if (combined_mod == 1) {
        combined_rem = x_pe % pe;
        combined_mod = pe;
      } else {
        const uint64_t rem_mod_pe = combined_rem % pe;
        const uint64_t diff = (x_pe >= rem_mod_pe) ? (x_pe - rem_mod_pe)
                                                   : (pe - (rem_mod_pe - x_pe));
        const uint64_t inv = utils::ModInv(combined_mod % pe, pe);
        const uint64_t t = utils::ModMul(diff, inv, pe);
        combined_rem += t * combined_mod;
        combined_mod *= pe;
      }
    }

    return combined_rem;
  }

  std::string BinaryPolyToString(ElementType poly) const {
    if (poly == 0) {
      return "0";
    }

    std::ostringstream oss;
    bool first = true;
    for (int i = static_cast<int>(m_) - 1; i >= 0; --i) {
      if ((static_cast<uint64_t>(poly) & (uint64_t{1} << i)) != 0) {
        if (!first) {
          oss << " + ";
        }
        first = false;
        if (i == 0) {
          oss << "1";
        } else if (i == 1) {
          oss << variable_name_;
        } else {
          oss << variable_name_ << "^" << i;
        }
      }
    }
    return oss.str();
  }

  // Parses any supported string representation into an external polynomial
  // value (never invokes virtual Pow/MultiplicativeGenerator so subclasses
  // like GFBEZechLogTables can convert once without double-logging).
  ElementType ParseElementStringExternal(const std::string &value_str) const {
    const std::string trimmed = detail::TrimAscii(value_str);
    if (trimmed.empty()) {
      throw std::invalid_argument("Empty string for GF(2^m) element");
    }

    // 1. Hexadecimal literal: "0x..." / "0X...".
    if (trimmed.size() >= 2 && trimmed[0] == '0' &&
        (trimmed[1] == 'x' || trimmed[1] == 'X')) {
      const std::string hex_digits = trimmed.substr(2);
      if (hex_digits.empty()) {
        throw std::invalid_argument("Invalid hex literal: " + value_str);
      }
      for (char ch : hex_digits) {
        if (std::isxdigit(static_cast<unsigned char>(ch)) == 0) {
          throw std::invalid_argument("Invalid hex literal: " + value_str);
        }
      }
      uint64_t val = 0;
      try {
        size_t parsed = 0;
        val = std::stoull(hex_digits, &parsed, 16);
        if (parsed != hex_digits.size()) {
          throw std::invalid_argument("Invalid hex literal: " + value_str);
        }
      } catch (const std::out_of_range &) {
        throw std::out_of_range("Hex literal out of range for GF(2^m): " +
                                value_str);
      }
      if (val > group_order_) {
        throw std::out_of_range("Hex literal exceeds field order: " +
                                value_str);
      }
      return static_cast<ElementType>(val);
    }

    // 2. Pure decimal unsigned integer literal: "0", "1", "5", ...
    bool all_digits = true;
    for (char ch : trimmed) {
      if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
        all_digits = false;
        break;
      }
    }
    if (all_digits) {
      uint64_t val = 0;
      try {
        size_t parsed = 0;
        val = std::stoull(trimmed, &parsed, 10);
        if (parsed != trimmed.size()) {
          throw std::invalid_argument("Invalid integer literal: " + value_str);
        }
      } catch (const std::out_of_range &) {
        throw std::out_of_range("Integer literal out of range for GF(2^m): " +
                                value_str);
      }
      if (val > group_order_) {
        throw std::out_of_range("Integer literal exceeds field order: " +
                                value_str);
      }
      return static_cast<ElementType>(val);
    }

    // 3. Generator power format: "<generator_name>^[-]<digits>" (e.g., "g^3",
    //    "g^-1"). Only matched when there is no '+' so polynomial sums with
    //    the same symbol still parse as polynomials.
    const std::string power_prefix = generator_name_ + "^";
    if (trimmed.rfind(power_prefix, 0) == 0 &&
        trimmed.find('+') == std::string::npos) {
      std::string exp_str = trimmed.substr(power_prefix.size());
      if (exp_str.empty()) {
        throw std::invalid_argument("Invalid power format: " + value_str);
      }
      bool is_negative = false;
      if (exp_str[0] == '-') {
        is_negative = true;
        exp_str.erase(0, 1);
        if (exp_str.empty()) {
          throw std::invalid_argument("Invalid power format: " + value_str);
        }
      }
      for (char ch : exp_str) {
        if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
          throw std::invalid_argument("Invalid exponent in power format: " +
                                      value_str);
        }
      }
      uint64_t raw_exp = 0;
      try {
        size_t parsed = 0;
        raw_exp = std::stoull(exp_str, &parsed, 10);
        if (parsed != exp_str.size()) {
          throw std::invalid_argument("Invalid exponent in power format: " +
                                      value_str);
        }
      } catch (const std::out_of_range &) {
        throw std::invalid_argument("Exponent out of range: " + value_str);
      }
      uint64_t reduced_exp = raw_exp % group_order_;
      if (is_negative && reduced_exp != 0) {
        reduced_exp = group_order_ - reduced_exp;
      }
      const ElementType gen = MultiplicativeGeneratorExternal();
      return PowExternal(gen, reduced_exp);
    }

    // 4. Polynomial string in variable_name_ (or x/X): e.g. "α^2 + α + 1".
    const auto parsed = detail::ParseBinaryPolynomial(
        trimmed, variable_name_, static_cast<int>(m_) - 1, /*allow_x=*/true);
    return static_cast<ElementType>(parsed.low_bits);
  }

  // 4-bit split nibble lookup table multiplication for GF(2^m) in uint8_t.
  // Uses SSSE3 PSHUFB (16 bytes per iteration) when __SSSE3__ is enabled, and
  // an auto-vectorizable portable 16-entry nibble split otherwise.
  void BuildNibbleTablesU8(uint8_t scalar, uint8_t lo[16],
                           uint8_t hi[16]) const noexcept {
    lo[0] = 0;
    hi[0] = 0;
    for (uint8_t k = 1; k < 16; ++k) {
      lo[k] = (static_cast<uint64_t>(k) <= group_order_) ? MulFast(k, scalar)
                                                         : uint8_t{0};
      const uint8_t hi_val = static_cast<uint8_t>(k << 4);
      hi[k] = (static_cast<uint64_t>(hi_val) <= group_order_)
                  ? MulFast(hi_val, scalar)
                  : uint8_t{0};
    }
  }

  void ScalarMulNibbleU8(std::span<uint8_t> dst, std::span<const uint8_t> a,
                         uint8_t scalar) const noexcept {
    alignas(16) uint8_t lo[16];
    alignas(16) uint8_t hi[16];
    BuildNibbleTablesU8(scalar, lo, hi);

    const size_t n = dst.size();
    uint8_t *out = dst.data();
    const uint8_t *pa = a.data();
    size_t i = 0;
#if defined(__SSSE3__)
    const __m128i vlo = _mm_load_si128(reinterpret_cast<const __m128i *>(lo));
    const __m128i vhi = _mm_load_si128(reinterpret_cast<const __m128i *>(hi));
    const __m128i mask = _mm_set1_epi8(0x0F);
    for (; i + 16 <= n; i += 16) {
      const __m128i va =
          _mm_loadu_si128(reinterpret_cast<const __m128i *>(pa + i));
      const __m128i a_lo = _mm_and_si128(va, mask);
      const __m128i a_hi = _mm_and_si128(_mm_srli_epi16(va, 4), mask);
      const __m128i prod = _mm_xor_si128(_mm_shuffle_epi8(vlo, a_lo),
                                         _mm_shuffle_epi8(vhi, a_hi));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(out + i), prod);
    }
#endif
    for (; i < n; ++i) {
      const uint8_t v = pa[i];
      out[i] = static_cast<uint8_t>(lo[v & 0x0FU] ^ hi[v >> 4]);
    }
  }

  void FmaNibbleU8(std::span<uint8_t> dst, std::span<const uint8_t> a,
                   uint8_t scalar) const noexcept {
    alignas(16) uint8_t lo[16];
    alignas(16) uint8_t hi[16];
    BuildNibbleTablesU8(scalar, lo, hi);

    const size_t n = dst.size();
    uint8_t *out = dst.data();
    const uint8_t *pa = a.data();
    size_t i = 0;
#if defined(__SSSE3__)
    const __m128i vlo = _mm_load_si128(reinterpret_cast<const __m128i *>(lo));
    const __m128i vhi = _mm_load_si128(reinterpret_cast<const __m128i *>(hi));
    const __m128i mask = _mm_set1_epi8(0x0F);
    for (; i + 16 <= n; i += 16) {
      const __m128i va =
          _mm_loadu_si128(reinterpret_cast<const __m128i *>(pa + i));
      const __m128i vd =
          _mm_loadu_si128(reinterpret_cast<const __m128i *>(out + i));
      const __m128i a_lo = _mm_and_si128(va, mask);
      const __m128i a_hi = _mm_and_si128(_mm_srli_epi16(va, 4), mask);
      const __m128i prod = _mm_xor_si128(_mm_shuffle_epi8(vlo, a_lo),
                                         _mm_shuffle_epi8(vhi, a_hi));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(out + i),
                       _mm_xor_si128(vd, prod));
    }
#endif
    for (; i < n; ++i) {
      const uint8_t v = pa[i];
      out[i] = static_cast<uint8_t>(out[i] ^ lo[v & 0x0FU] ^ hi[v >> 4]);
    }
  }

  void BuildNibbleTablesU16(uint16_t scalar, uint16_t t0[16], uint16_t t1[16],
                            uint16_t t2[16], uint16_t t3[16]) const noexcept {
    t0[0] = t1[0] = t2[0] = t3[0] = 0;
    for (uint16_t k = 1; k < 16; ++k) {
      const uint16_t v0 = k;
      const uint16_t v1 = static_cast<uint16_t>(k << 4);
      const uint16_t v2 = static_cast<uint16_t>(k << 8);
      const uint16_t v3 = static_cast<uint16_t>(k << 12);
      t0[k] = (static_cast<uint64_t>(v0) <= group_order_) ? MulFast(v0, scalar)
                                                          : uint16_t{0};
      t1[k] = (static_cast<uint64_t>(v1) <= group_order_) ? MulFast(v1, scalar)
                                                          : uint16_t{0};
      t2[k] = (static_cast<uint64_t>(v2) <= group_order_) ? MulFast(v2, scalar)
                                                          : uint16_t{0};
      t3[k] = (static_cast<uint64_t>(v3) <= group_order_) ? MulFast(v3, scalar)
                                                          : uint16_t{0};
    }
  }

  void ScalarMulNibbleU16(std::span<uint16_t> dst, std::span<const uint16_t> a,
                          uint16_t scalar) const noexcept {
    uint16_t t0[16], t1[16], t2[16], t3[16];
    BuildNibbleTablesU16(scalar, t0, t1, t2, t3);
    const size_t n = dst.size();
    uint16_t *out = dst.data();
    const uint16_t *pa = a.data();
    for (size_t i = 0; i < n; ++i) {
      const uint16_t v = pa[i];
      out[i] = static_cast<uint16_t>(t0[v & 0xFU] ^ t1[(v >> 4) & 0xFU] ^
                                     t2[(v >> 8) & 0xFU] ^ t3[v >> 12]);
    }
  }

  void FmaNibbleU16(std::span<uint16_t> dst, std::span<const uint16_t> a,
                    uint16_t scalar) const noexcept {
    uint16_t t0[16], t1[16], t2[16], t3[16];
    BuildNibbleTablesU16(scalar, t0, t1, t2, t3);
    const size_t n = dst.size();
    uint16_t *out = dst.data();
    const uint16_t *pa = a.data();
    for (size_t i = 0; i < n; ++i) {
      const uint16_t v = pa[i];
      out[i] =
          static_cast<uint16_t>(out[i] ^ t0[v & 0xFU] ^ t1[(v >> 4) & 0xFU] ^
                                t2[(v >> 8) & 0xFU] ^ t3[v >> 12]);
    }
  }
};

template <UnsignedFieldWord ElementType = uint32_t>
class GFBELogTables : public GaloisFieldBinaryExtension<ElementType> {
 public:
  static constexpr uint8_t kMaxTableDegree = 22;

  explicit GFBELogTables(uint8_t m, const std::string &rep = "int",
                         const std::string &irreducible_poly = "",
                         const std::string &variable_name = "α",
                         bool check_irreducible = false,
                         const std::string &generator_name = "g")
      : GFBELogTables(NoInitTag{}, m, rep, irreducible_poly, variable_name,
                      check_irreducible, generator_name) {
    InitTables();
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    if (a == 0 || b == 0) {
      return ElementType{0};
    }
    return exp_table_[static_cast<size_t>(log_table_[a]) +
                      static_cast<size_t>(log_table_[b])];
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }
    if (a == 0) {
      return ElementType{0};
    }
    return exp_table_[static_cast<size_t>(this->group_order_) +
                      static_cast<size_t>(log_table_[a]) -
                      static_cast<size_t>(log_table_[b])];
  }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }
    return exp_table_[static_cast<size_t>(this->group_order_) -
                      static_cast<size_t>(log_table_[a])];
  }

  inline ElementType Pow(const ElementType &a, uint64_t exp) const override {
    if (a == 0) {
      return (exp == 0) ? ElementType{1} : ElementType{0};
    }
    if (a == 1 || exp == 0 || this->group_order_ == 1) {
      return ElementType{1};
    }
    const uint64_t idx =
        ((exp % this->group_order_) * static_cast<uint64_t>(log_table_[a])) %
        this->group_order_;
    return exp_table_[static_cast<size_t>(idx)];
  }

  inline ElementType Sqrt(const ElementType &a) const override {
    if (a == 0 || a == 1 || this->m_ == 1) {
      return a;
    }
    const uint64_t la = static_cast<uint64_t>(log_table_[a]);
    const uint64_t half =
        ((la & 1ULL) == 0) ? (la >> 1) : ((la + this->group_order_) >> 1);
    return exp_table_[static_cast<size_t>(half)];
  }

  uint64_t Log(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }
    return static_cast<uint64_t>(log_table_[a]);
  }

  uint64_t Log(const ElementType &a,
               const ElementType &generator) const override {
    if (a == 0) {
      throw std::domain_error("Log of zero is undefined");
    }
    if (generator == 0) {
      throw std::invalid_argument("Generator cannot be zero");
    }
    if (this->m_ == 1) {
      if (a != 1) {
        throw std::invalid_argument("Element out of range for GF(2^1)");
      }
      if (generator != 1) {
        throw std::invalid_argument("Generator in GF(2^1) must be 1");
      }
      return 0;
    }
    if (generator == generator_) {
      return static_cast<uint64_t>(log_table_[a]);
    }
    if (a == 1) {
      return 0;
    }
    if (generator == 1) {
      throw std::runtime_error("Element not in the group or invalid generator");
    }

    const uint64_t n = this->group_order_;
    const uint64_t log_a = static_cast<uint64_t>(log_table_[a]);
    const uint64_t log_g = static_cast<uint64_t>(log_table_[generator]);
    const uint64_t d = std::gcd(log_g, n);
    if (log_a % d != 0) {
      throw std::runtime_error("Element not in the group or invalid generator");
    }
    const uint64_t mod = n / d;
    const uint64_t inv_g = utils::ModInv(log_g / d, mod);
    return utils::ModMul(log_a / d, inv_g, mod);
  }

  ElementType MultiplicativeGenerator() const override { return generator_; }

  size_t MemoryBytes() const override {
    return (sizeof(*this) - sizeof(GaloisFieldBinaryExtension<ElementType>)) +
           GaloisFieldBinaryExtension<ElementType>::MemoryBytes() +
           log_table_.capacity() * sizeof(uint32_t) +
           exp_table_.capacity() * sizeof(ElementType);
  }

  void AddArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    GaloisFieldBinaryExtension<ElementType>::AddArray(dst, a, b);
  }

  void SubArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    GaloisFieldBinaryExtension<ElementType>::SubArray(dst, a, b);
  }

  void MulArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    const uint32_t *log_tbl = log_table_.data();
    const ElementType *exp_tbl = exp_table_.data();
    for (size_t i = 0; i < n; ++i) {
      const ElementType av = pa[i];
      const ElementType bv = pb[i];
      out[i] = (av == 0 || bv == 0) ? ElementType{0}
                                    : exp_tbl[static_cast<size_t>(log_tbl[av]) +
                                              static_cast<size_t>(log_tbl[bv])];
    }
  }

  void ScalarMulArray(std::span<ElementType> dst,
                      std::span<const ElementType> a,
                      const ElementType &scalar) const override {
    if constexpr (std::is_same_v<ElementType, uint8_t>) {
      GaloisFieldBinaryExtension<ElementType>::ScalarMulArray(dst, a, scalar);
      return;
    }
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match");
    }
    if (scalar == 0) {
      std::fill(dst.begin(), dst.end(), ElementType{0});
      return;
    }
    if (scalar == 1) {
      if (dst.data() != a.data()) {
        std::copy(a.begin(), a.end(), dst.begin());
      }
      return;
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint32_t *log_tbl = log_table_.data();
    const ElementType *shifted =
        exp_table_.data() + static_cast<size_t>(log_tbl[scalar]);
    for (size_t i = 0; i < n; ++i) {
      const ElementType av = pa[i];
      out[i] = (av == 0) ? ElementType{0} : shifted[log_tbl[av]];
    }
  }

  void FmaArray(std::span<ElementType> dst, std::span<const ElementType> a,
                const ElementType &scalar) const override {
    if constexpr (std::is_same_v<ElementType, uint8_t>) {
      GaloisFieldBinaryExtension<ElementType>::FmaArray(dst, a, scalar);
      return;
    }
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match");
    }
    if (scalar == 0) {
      return;
    }
    if (scalar == 1) {
      ElementType *out = dst.data();
      const ElementType *pa = a.data();
      for (size_t i = 0; i < n; ++i) {
        out[i] = static_cast<ElementType>(out[i] ^ pa[i]);
      }
      return;
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint32_t *log_tbl = log_table_.data();
    const ElementType *shifted =
        exp_table_.data() + static_cast<size_t>(log_tbl[scalar]);
    for (size_t i = 0; i < n; ++i) {
      const ElementType av = pa[i];
      if (av != 0) {
        out[i] = static_cast<ElementType>(out[i] ^ shifted[log_tbl[av]]);
      }
    }
  }

 protected:
  struct NoInitTag {};

  // Protected constructor allowing derived classes (GFBELogTablesOpt) to
  // initialize their own table layout once without virtual calls in ctor.
  GFBELogTables(NoInitTag, uint8_t m, const std::string &rep,
                const std::string &irreducible_poly,
                const std::string &variable_name, bool check_irreducible,
                const std::string &generator_name)
      : GaloisFieldBinaryExtension<ElementType>(
            m, rep, irreducible_poly, variable_name, check_irreducible,
            generator_name),
        generator_(GaloisFieldBinaryExtension<
                   ElementType>::MultiplicativeGenerator()) {
    if (this->m_ > kMaxTableDegree) {
      throw std::invalid_argument(
          "Log table degree must be <= " + std::to_string(kMaxTableDegree) +
          " to prevent excessive memory allocation");
    }
  }

  std::vector<uint32_t> log_table_;
  std::vector<ElementType> exp_table_;
  ElementType generator_;

  static constexpr uint32_t LOG_ZERO = std::numeric_limits<uint32_t>::max();

 private:
  void InitTables() {
    const size_t ord = static_cast<size_t>(this->order_);
    const size_t g_ord = static_cast<size_t>(this->group_order_);
    log_table_.assign(ord, LOG_ZERO);
    exp_table_.assign(2 * g_ord, ElementType{0});

    ElementType element = 1;
    for (size_t i = 0; i < g_ord; ++i) {
      exp_table_[i] = element;
      exp_table_[i + g_ord] = element;
      log_table_[element] = static_cast<uint32_t>(i);
      element = this->MulFast(element, generator_);
    }
    log_table_[0] = LOG_ZERO;
  }
};

template <UnsignedFieldWord ElementType = uint32_t>
class GFBELogTablesOpt final : public GFBELogTables<ElementType> {
 public:
  explicit GFBELogTablesOpt(uint8_t m, const std::string &rep = "int",
                            const std::string &irreducible_poly = "",
                            const std::string &variable_name = "α",
                            bool check_irreducible = false,
                            const std::string &generator_name = "g")
      : GFBELogTables<ElementType>(
            typename GFBELogTables<ElementType>::NoInitTag{}, m, rep,
            irreducible_poly, variable_name, check_irreducible,
            generator_name) {
    InitOptTables();
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    return this->exp_table_[static_cast<size_t>(this->log_table_[a]) +
                            static_cast<size_t>(this->log_table_[b])];
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    if (b == 0) {
      throw std::domain_error("Division by zero");
    }
    return this->exp_table_[static_cast<size_t>(this->group_order_) +
                            static_cast<size_t>(this->log_table_[a]) -
                            static_cast<size_t>(this->log_table_[b])];
  }

  inline ElementType Inv(const ElementType &a) const override {
    if (a == 0) {
      throw std::domain_error("Inverse of zero is undefined");
    }
    return this->exp_table_[static_cast<size_t>(this->group_order_) -
                            static_cast<size_t>(this->log_table_[a])];
  }

  void MulArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    const uint32_t *log_tbl = this->log_table_.data();
    const ElementType *exp_tbl = this->exp_table_.data();
    for (size_t i = 0; i < n; ++i) {
      out[i] = exp_tbl[static_cast<size_t>(log_tbl[pa[i]]) +
                       static_cast<size_t>(log_tbl[pb[i]])];
    }
  }

  void ScalarMulArray(std::span<ElementType> dst,
                      std::span<const ElementType> a,
                      const ElementType &scalar) const override {
    if constexpr (std::is_same_v<ElementType, uint8_t>) {
      GaloisFieldBinaryExtension<ElementType>::ScalarMulArray(dst, a, scalar);
      return;
    }
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match");
    }
    if (scalar == 0) {
      std::fill(dst.begin(), dst.end(), ElementType{0});
      return;
    }
    if (scalar == 1) {
      if (dst.data() != a.data()) {
        std::copy(a.begin(), a.end(), dst.begin());
      }
      return;
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint32_t *log_tbl = this->log_table_.data();
    const ElementType *shifted =
        this->exp_table_.data() + static_cast<size_t>(log_tbl[scalar]);
    for (size_t i = 0; i < n; ++i) {
      out[i] = shifted[log_tbl[pa[i]]];
    }
  }

  void FmaArray(std::span<ElementType> dst, std::span<const ElementType> a,
                const ElementType &scalar) const override {
    if constexpr (std::is_same_v<ElementType, uint8_t>) {
      GaloisFieldBinaryExtension<ElementType>::FmaArray(dst, a, scalar);
      return;
    }
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match");
    }
    if (scalar == 0) {
      return;
    }
    if (scalar == 1) {
      ElementType *out = dst.data();
      const ElementType *pa = a.data();
      for (size_t i = 0; i < n; ++i) {
        out[i] = static_cast<ElementType>(out[i] ^ pa[i]);
      }
      return;
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint32_t *log_tbl = this->log_table_.data();
    const ElementType *shifted =
        this->exp_table_.data() + static_cast<size_t>(log_tbl[scalar]);
    for (size_t i = 0; i < n; ++i) {
      out[i] = static_cast<ElementType>(out[i] ^ shifted[log_tbl[pa[i]]]);
    }
  }

 private:
  void InitOptTables() {
    const size_t ord = static_cast<size_t>(this->order_);
    const size_t g_ord = static_cast<size_t>(this->group_order_);
    this->log_table_.assign(ord, 0U);
    // Size 4 * g_ord + 1 guarantees zero OOB reads for all pairs (a, b):
    // - Mul(0, 0) accesses 2*g_ord + 2*g_ord = 4*g_ord -> 0
    // - Div(0, 1) accesses g_ord + 2*g_ord - 0 = 3*g_ord -> 0
    // - Div(a, a) accesses g_ord + log_a - log_a = g_ord -> 1
    // - Inv(1) accesses g_ord - 0 = g_ord -> 1
    this->exp_table_.assign(4 * g_ord + 1, ElementType{0});

    ElementType element = 1;
    for (size_t i = 0; i < g_ord; ++i) {
      this->exp_table_[i] = element;
      this->exp_table_[i + g_ord] = element;
      this->log_table_[element] = static_cast<uint32_t>(i);
      element = this->MulFast(element, this->generator_);
    }
    this->log_table_[0] = static_cast<uint32_t>(2 * g_ord);
  }
};

template <UnsignedFieldWord ElementType = uint32_t>
class GFBEZechLogTables final : public GaloisFieldBinaryExtension<ElementType> {
 public:
  static constexpr uint8_t kMaxTableDegree = 22;
  static constexpr ElementType LOG_ZERO =
      std::numeric_limits<ElementType>::max();

  explicit GFBEZechLogTables(uint8_t m, const std::string &rep = "log",
                             const std::string &irreducible_poly = "",
                             const std::string &variable_name = "α",
                             bool check_irreducible = false,
                             const std::string &generator_name = "g")
      : GaloisFieldBinaryExtension<ElementType>(
            m, rep, irreducible_poly, variable_name, check_irreducible,
            generator_name),
        external_generator_(GaloisFieldBinaryExtension<
                            ElementType>::MultiplicativeGenerator()) {
    if (this->m_ > kMaxTableDegree) {
      throw std::invalid_argument("Zech log table degree must be <= " +
                                  std::to_string(kMaxTableDegree) +
                                  " to prevent excessive memory allocation");
    }
    InitZechTables();
  }

  inline ElementType Add(const ElementType &log_a,
                         const ElementType &log_b) const override {
    if (log_a == LOG_ZERO) {
      return log_b;
    }
    if (log_b == LOG_ZERO) {
      return log_a;
    }
    const uint32_t la = static_cast<uint32_t>(log_a);
    const uint32_t lb = static_cast<uint32_t>(log_b);
    const uint32_t g_ord = static_cast<uint32_t>(this->group_order_);
    const uint32_t diff = (lb >= la) ? (lb - la) : (g_ord - la + lb);
    const ElementType z = zech_table_[diff];
    if (z == LOG_ZERO) {
      return LOG_ZERO;
    }
    uint32_t sum = la + static_cast<uint32_t>(z);
    if (sum >= g_ord) {
      sum -= g_ord;
    }
    return static_cast<ElementType>(sum);
  }

  inline ElementType Sub(const ElementType &log_a,
                         const ElementType &log_b) const override {
    return Add(log_a, log_b);
  }

  inline ElementType Mul(const ElementType &log_a,
                         const ElementType &log_b) const override {
    if (log_a == LOG_ZERO || log_b == LOG_ZERO) {
      return LOG_ZERO;
    }
    const uint32_t g_ord = static_cast<uint32_t>(this->group_order_);
    uint32_t sum = static_cast<uint32_t>(log_a) + static_cast<uint32_t>(log_b);
    if (sum >= g_ord) {
      sum -= g_ord;
    }
    return static_cast<ElementType>(sum);
  }

  inline ElementType Div(const ElementType &log_a,
                         const ElementType &log_b) const override {
    if (log_b == LOG_ZERO) {
      throw std::domain_error("Division by zero");
    }
    if (log_a == LOG_ZERO) {
      return LOG_ZERO;
    }
    const uint32_t la = static_cast<uint32_t>(log_a);
    const uint32_t lb = static_cast<uint32_t>(log_b);
    const uint32_t g_ord = static_cast<uint32_t>(this->group_order_);
    const uint32_t diff = (la >= lb) ? (la - lb) : (la + g_ord - lb);
    return static_cast<ElementType>(diff);
  }

  inline ElementType Neg(const ElementType &log_a) const override {
    return log_a;
  }

  inline ElementType Inv(const ElementType &log_a) const override {
    if (log_a == LOG_ZERO) {
      throw std::domain_error("Inverse of zero is undefined");
    }
    return (log_a == 0)
               ? ElementType{0}
               : static_cast<ElementType>(this->group_order_ -
                                          static_cast<uint64_t>(log_a));
  }

  inline ElementType Pow(const ElementType &log_a,
                         uint64_t exp) const override {
    if (log_a == LOG_ZERO) {
      return (exp == 0) ? ElementType{0} : LOG_ZERO;
    }
    if (log_a == 0 || exp == 0 || this->group_order_ == 1) {
      return ElementType{0};
    }
    const uint64_t res =
        ((exp % this->group_order_) * static_cast<uint64_t>(log_a)) %
        this->group_order_;
    return static_cast<ElementType>(res);
  }

  inline ElementType Sqrt(const ElementType &log_a) const override {
    if (log_a == LOG_ZERO || log_a == 0 || this->m_ == 1) {
      return log_a;
    }
    const uint64_t la = static_cast<uint64_t>(log_a);
    const uint64_t half =
        ((la & 1ULL) == 0) ? (la >> 1) : ((la + this->group_order_) >> 1);
    return static_cast<ElementType>(half);
  }

  uint64_t Log(const ElementType &log_a) const override {
    if (log_a == LOG_ZERO) {
      throw std::domain_error("Log of zero is undefined");
    }
    return static_cast<uint64_t>(log_a);
  }

  uint64_t Log(const ElementType &log_a,
               const ElementType &gen_log) const override {
    if (log_a == LOG_ZERO) {
      throw std::domain_error("Log of zero is undefined");
    }
    if (gen_log == LOG_ZERO) {
      throw std::invalid_argument("Generator cannot be zero");
    }
    if (this->m_ == 1) {
      if (log_a != 0) {
        throw std::invalid_argument("Invalid element in GF(2^1)");
      }
      if (gen_log != 0) {
        throw std::invalid_argument("Generator in GF(2^1) must be 1");
      }
      return 0;
    }
    if (gen_log == 1) {
      return static_cast<uint64_t>(log_a);
    }
    if (log_a == 0) {
      return 0;
    }
    if (gen_log == 0) {
      throw std::runtime_error("Element not in the group or invalid generator");
    }

    const uint64_t n = this->group_order_;
    const uint64_t la = static_cast<uint64_t>(log_a);
    const uint64_t lg = static_cast<uint64_t>(gen_log);
    const uint64_t d = std::gcd(lg, n);
    if (la % d != 0) {
      throw std::runtime_error("Element not in the group or invalid generator");
    }
    const uint64_t mod = n / d;
    const uint64_t inv_g = utils::ModInv(lg / d, mod);
    return utils::ModMul(la / d, inv_g, mod);
  }

  ElementType AdditiveIdentity() const override { return LOG_ZERO; }
  ElementType MultiplicativeIdentity() const override { return ElementType{0}; }

  bool IsZero(const ElementType &a) const override { return a == LOG_ZERO; }
  bool IsOne(const ElementType &a) const override {
    return a == ElementType{0};
  }

  ElementType MultiplicativeGenerator() const override {
    return (this->m_ == 1) ? ElementType{0} : ElementType{1};
  }

  std::vector<ElementType> MultiplicativeGenerators() const override {
    if (this->m_ == 1) {
      return {ElementType{0}};
    }
    std::vector<ElementType> generators;
    const uint64_t n = this->group_order_;
    for (uint64_t k = 1; k < n; ++k) {
      if (std::gcd(k, n) == 1) {
        generators.push_back(static_cast<ElementType>(k));
      }
    }
    return generators;
  }

  ElementType Random() const override {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    return Random(gen);
  }

  ElementType Random(std::mt19937_64 &rng) const override {
    std::uniform_int_distribution<uint64_t> distrib(0ULL, this->group_order_);
    const uint64_t r = distrib(rng);
    return (r == this->group_order_) ? LOG_ZERO : static_cast<ElementType>(r);
  }

  ElementType GetElementValue(const ElementType &log_val) const override {
    if (log_val == LOG_ZERO) {
      return ElementType{0};
    }
    if (static_cast<uint64_t>(log_val) >= this->group_order_) {
      throw std::out_of_range("Zech log element value out of range");
    }
    return exp_table_[log_val];
  }

  ElementType SetElementValue(const ElementType &external_val) const override {
    if (static_cast<uint64_t>(external_val) > this->group_order_) {
      throw std::out_of_range("Element value exceeds field order");
    }
    return log_table_[external_val];
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    const ElementType external_val =
        this->ParseElementStringExternal(value_str);
    return log_table_[external_val];
  }

  void Print(const ElementType &log_a, std::ostream &os) const override {
    detail::BinaryStreamFlagsGuard guard(os);
    const FieldRepresentation rep = this->GetRepresentation();
    if (rep == FieldRepresentation::POW) {
      if (log_a == LOG_ZERO) {
        os << "0";
      } else {
        os << this->generator_name_ << "^" << static_cast<uint64_t>(log_a);
      }
      return;
    }
    if (rep == FieldRepresentation::LOG) {
      if (log_a == LOG_ZERO) {
        os << "-inf";
      } else {
        os << static_cast<uint64_t>(log_a);
      }
      return;
    }

    const ElementType ext = GetElementValue(log_a);
    if (rep == FieldRepresentation::INT) {
      os << static_cast<uint64_t>(ext);
    } else if (rep == FieldRepresentation::HEX) {
      os << "0x" << std::hex << static_cast<uint64_t>(ext);
    } else if (rep == FieldRepresentation::POLY) {
      os << this->BinaryPolyToString(ext);
    } else {
      throw std::invalid_argument("Unsupported representation for GF(2^m)");
    }
  }

  std::string ToString(const ElementType &log_a) const override {
    std::ostringstream oss;
    Print(log_a, oss);
    return oss.str();
  }

  size_t MemoryBytes() const override {
    return (sizeof(*this) - sizeof(GaloisFieldBinaryExtension<ElementType>)) +
           GaloisFieldBinaryExtension<ElementType>::MemoryBytes() +
           (zech_table_.capacity() + exp_table_.capacity() +
            log_table_.capacity()) *
               sizeof(ElementType);
  }

  void AddArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("AddArray: span sizes must match");
    }
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Add(a[i], b[i]);
    }
  }

  void SubArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    AddArray(dst, a, b);
  }

  void MulArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match");
    }
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Mul(a[i], b[i]);
    }
  }

  void ScalarMulArray(std::span<ElementType> dst,
                      std::span<const ElementType> a,
                      const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match");
    }
    if (scalar == LOG_ZERO) {
      std::fill(dst.begin(), dst.end(), LOG_ZERO);
      return;
    }
    if (scalar == ElementType{0}) {
      if (dst.data() != a.data()) {
        std::copy(a.begin(), a.end(), dst.begin());
      }
      return;
    }
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Mul(a[i], scalar);
    }
  }

  void FmaArray(std::span<ElementType> dst, std::span<const ElementType> a,
                const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match");
    }
    if (scalar == LOG_ZERO) {
      return;
    }
    if (scalar == ElementType{0}) {
      for (size_t i = 0; i < n; ++i) {
        dst[i] = Add(dst[i], a[i]);
      }
      return;
    }
    for (size_t i = 0; i < n; ++i) {
      dst[i] = Add(dst[i], Mul(a[i], scalar));
    }
  }

  ElementType DotProduct(std::span<const ElementType> a,
                         std::span<const ElementType> b) const override {
    const size_t n = a.size();
    if (b.size() != n) {
      throw std::invalid_argument("DotProduct: span sizes must match");
    }
    ElementType acc = LOG_ZERO;
    for (size_t i = 0; i < n; ++i) {
      acc = Add(acc, Mul(a[i], b[i]));
    }
    return acc;
  }

 protected:
  std::vector<ElementType> zech_table_;
  std::vector<ElementType> exp_table_;
  std::vector<ElementType> log_table_;
  ElementType external_generator_;

 private:
  void InitZechTables() {
    const size_t g_ord = static_cast<size_t>(this->group_order_);
    const size_t ord = static_cast<size_t>(this->order_);
    zech_table_.assign(g_ord, LOG_ZERO);
    exp_table_.assign(g_ord, ElementType{0});
    log_table_.assign(ord, LOG_ZERO);

    // Single pass computing g^i once for both exp/log and Zech tables.
    ElementType element = 1;
    for (size_t i = 0; i < g_ord; ++i) {
      exp_table_[i] = element;
      log_table_[element] = static_cast<ElementType>(i);
      element = this->MulFast(element, external_generator_);
    }

    for (size_t i = 0; i < g_ord; ++i) {
      const ElementType one_plus =
          static_cast<ElementType>(ElementType{1} ^ exp_table_[i]);
      zech_table_[i] = (one_plus == 0) ? LOG_ZERO : log_table_[one_plus];
    }
  }
};

using GF2 = GaloisFieldBinary;
template <UnsignedFieldWord ElementType = uint32_t>
using GF2X = GaloisFieldBinaryExtension<ElementType>;
template <UnsignedFieldWord ElementType = uint32_t>
using GF2XLOG = GFBELogTablesOpt<ElementType>;
template <UnsignedFieldWord ElementType = uint32_t>
using GF2XZECH = GFBEZechLogTables<ElementType>;

static_assert(FiniteField<GaloisFieldBinary>);
static_assert(FiniteField<GaloisFieldBinaryExtension<uint8_t>>);
static_assert(FiniteField<GaloisFieldBinaryExtension<uint16_t>>);
static_assert(FiniteField<GaloisFieldBinaryExtension<uint32_t>>);
static_assert(FiniteField<GaloisFieldBinaryExtension<uint64_t>>);
static_assert(FiniteField<GFBELogTables<uint8_t>>);
static_assert(FiniteField<GFBELogTables<uint16_t>>);
static_assert(FiniteField<GFBELogTables<uint32_t>>);
static_assert(FiniteField<GFBELogTablesOpt<uint8_t>>);
static_assert(FiniteField<GFBELogTablesOpt<uint16_t>>);
static_assert(FiniteField<GFBELogTablesOpt<uint32_t>>);
static_assert(FiniteField<GFBEZechLogTables<uint8_t>>);
static_assert(FiniteField<GFBEZechLogTables<uint16_t>>);
static_assert(FiniteField<GFBEZechLogTables<uint32_t>>);

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_BINARY_HPP_
