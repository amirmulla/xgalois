#ifndef XGALOIS_FIELD_GF_PRIME_HPP_
#define XGALOIS_FIELD_GF_PRIME_HPP_

#include <algorithm>
#include <atomic>
#include <bit>
#include <cctype>
#include <cmath>
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
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

template <UnsignedFieldWord ElementType = uint32_t>
class GaloisFieldPrime : public GaloisFieldBase<ElementType> {
  static_assert(sizeof(ElementType) <= sizeof(uint64_t),
                "ElementType must be at most 64 bits.");

 protected:
  __extension__ using u128 = unsigned __int128;

  static constexpr uint64_t kMaxPrime = (uint64_t{1} << 63) - 1;

 public:
  explicit GaloisFieldPrime(ElementType prime, const std::string &rep = "int",
                            bool prime_testing = true,
                            const std::string &generator_name = "g")
      : p_(prime),
        p64_(static_cast<uint64_t>(prime)),
        representation_(
            ValidateRepresentation(utils::ConvertRepresentation(rep))),
        generator_(static_cast<ElementType>(0)),
        generator_name_(ValidateGeneratorName(generator_name)) {
    if (p64_ < 2) {
      throw std::invalid_argument("Prime must be >= 2.");
    }
    if (p64_ > kMaxPrime) {
      throw std::invalid_argument("Prime must be <= 2^63 - 1.");
    }
    if (prime_testing && !utils::IsPrime(p64_)) {
      throw std::invalid_argument("Provided value is not prime.");
    }

    // Precompute Barrett reduction constants:
    // 1) barrett_mu_64_ = floor(2^64 / p) for exact 64-bit -> [0, p-1]
    //    reduction (used directly for <=32-bit Mul and general input
    //    reduction).
    barrett_mu_64_ = static_cast<uint64_t>((static_cast<u128>(1) << 64) / p64_);

    // 2) For 64-bit primes up to 2^63 - 1, let m = bit_width(p - 1) in [1, 63]
    //    so 2^(m-1) < p <= 2^m. Then barrett_shift_ = m - 1 and
    //    barrett_mu_128_ = floor(2^(m + 63) / p) < 2^64.
    const unsigned m = static_cast<unsigned>(std::bit_width(p64_ - 1));
    barrett_shift_ = m - 1;
    barrett_mu_128_ =
        static_cast<uint64_t>((static_cast<u128>(1) << (m + 63)) / p64_);

    if (p64_ == 2) {
      generator_.store(static_cast<ElementType>(1), std::memory_order_relaxed);
    }
  }

  GaloisFieldPrime(const GaloisFieldPrime &other)
      : p_(other.p_),
        p64_(other.p64_),
        barrett_mu_64_(other.barrett_mu_64_),
        barrett_mu_128_(other.barrett_mu_128_),
        barrett_shift_(other.barrett_shift_),
        representation_(other.representation_.load(std::memory_order_relaxed)),
        generator_(other.generator_.load(std::memory_order_acquire)),
        generator_name_(other.generator_name_) {}

  GaloisFieldPrime &operator=(const GaloisFieldPrime &other) {
    if (this != &other) {
      p_ = other.p_;
      p64_ = other.p64_;
      barrett_mu_64_ = other.barrett_mu_64_;
      barrett_mu_128_ = other.barrett_mu_128_;
      barrett_shift_ = other.barrett_shift_;
      representation_.store(
          other.representation_.load(std::memory_order_relaxed),
          std::memory_order_relaxed);
      generator_.store(other.generator_.load(std::memory_order_acquire),
                       std::memory_order_release);
      generator_name_ = other.generator_name_;
    }
    return *this;
  }

  uint64_t Characteristic() const override { return p64_; }
  uint64_t Order() const override { return p64_; }
  uint32_t Degree() const override { return 1; }
  uint64_t Modulus() const { return p64_; }
  const std::string &GeneratorName() const { return generator_name_; }

  ElementType MultiplicativeIdentity() const override {
    return static_cast<ElementType>(1);
  }

  ElementType AdditiveIdentity() const override {
    return static_cast<ElementType>(0);
  }

  bool IsZero(const ElementType &a) const override {
    return Reduce64(static_cast<uint64_t>(a)) == 0;
  }

  bool IsOne(const ElementType &a) const override {
    return Reduce64(static_cast<uint64_t>(a)) == 1;
  }

  ElementType GetElementValue(const ElementType &value) const override {
    return static_cast<ElementType>(Reduce64(static_cast<uint64_t>(value)));
  }

  ElementType SetElementValue(const ElementType &value) const override {
    return static_cast<ElementType>(Reduce64(static_cast<uint64_t>(value)));
  }

  ElementType SetElementValue(const std::string &value_str) const override {
    if (value_str.empty()) {
      throw std::invalid_argument("Empty string is not a valid GF(p) element.");
    }

    // Power notation: "<generator_name>^<signed_int>"
    const size_t caret = value_str.find('^');
    if (caret != std::string::npos) {
      const std::string base = value_str.substr(0, caret);
      if (base != generator_name_) {
        throw std::invalid_argument("Unexpected generator in power format: " +
                                    value_str);
      }
      const std::string exp_str = value_str.substr(caret + 1);
      if (exp_str.empty()) {
        throw std::invalid_argument("Invalid power format: " + value_str);
      }

      bool is_negative = false;
      size_t idx = 0;
      if (exp_str[0] == '-') {
        is_negative = true;
        idx = 1;
        if (idx >= exp_str.size()) {
          throw std::invalid_argument("Invalid power format: " + value_str);
        }
      }

      for (size_t i = idx; i < exp_str.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(exp_str[i])) == 0) {
          throw std::invalid_argument("Invalid exponent format: " + value_str);
        }
      }

      uint64_t raw_exp = 0;
      try {
        size_t parsed = 0;
        raw_exp = std::stoull(exp_str.substr(idx), &parsed, 10);
        if (parsed != exp_str.size() - idx) {
          throw std::invalid_argument("Invalid exponent format: " + value_str);
        }
      } catch (const std::out_of_range &) {
        throw std::invalid_argument("Exponent out of range: " + value_str);
      }

      if (p64_ == 2) {
        return static_cast<ElementType>(1);
      }
      const uint64_t phi = p64_ - 1;
      uint64_t reduced_exp = raw_exp % phi;
      if (is_negative && reduced_exp != 0) {
        reduced_exp = phi - reduced_exp;
      }
      return Pow(MultiplicativeGenerator(), reduced_exp);
    }

    // Optional leading minus sign for integer/hex literals.
    bool is_negative = false;
    size_t pos = 0;
    if (value_str[0] == '-') {
      is_negative = true;
      pos = 1;
      if (pos >= value_str.size()) {
        throw std::invalid_argument("Invalid GF(p) element string: " +
                                    value_str);
      }
    }

    // Hexadecimal notation: "0x..." / "0X..." (with optional leading '-')
    if (value_str.size() - pos >= 2 && value_str[pos] == '0' &&
        (value_str[pos + 1] == 'x' || value_str[pos + 1] == 'X')) {
      pos += 2;
      if (pos >= value_str.size()) {
        throw std::invalid_argument("Invalid hex format for GF(p) element: " +
                                    value_str);
      }
      uint64_t acc = 0;
      for (size_t i = pos; i < value_str.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(value_str[i]);
        uint64_t digit = 0;
        if (c >= '0' && c <= '9') {
          digit = static_cast<uint64_t>(c - '0');
        } else if (c >= 'a' && c <= 'f') {
          digit = static_cast<uint64_t>(10 + (c - 'a'));
        } else if (c >= 'A' && c <= 'F') {
          digit = static_cast<uint64_t>(10 + (c - 'A'));
        } else {
          throw std::invalid_argument("Invalid hex digit in GF(p) element: " +
                                      value_str);
        }
        acc =
            static_cast<uint64_t>((static_cast<u128>(acc) * 16 + digit) % p64_);
      }
      if (is_negative && acc != 0) {
        acc = p64_ - acc;
      }
      return static_cast<ElementType>(acc);
    }

    // Decimal integer notation: "[ - ]?[0-9]+"
    uint64_t acc = 0;
    for (size_t i = pos; i < value_str.size(); ++i) {
      const unsigned char c = static_cast<unsigned char>(value_str[i]);
      if (std::isdigit(c) == 0) {
        throw std::invalid_argument(
            "Unsupported string format for GF(p) element: " + value_str);
      }
      acc = static_cast<uint64_t>((static_cast<u128>(acc) * 10 + (c - '0')) %
                                  p64_);
    }
    if (is_negative && acc != 0) {
      acc = p64_ - acc;
    }
    return static_cast<ElementType>(acc);
  }

  FieldRepresentation GetRepresentation() const override {
    return representation_.load(std::memory_order_relaxed);
  }

  void SetRepresentation(FieldRepresentation rep) override {
    representation_.store(ValidateRepresentation(rep),
                          std::memory_order_relaxed);
  }

  inline ElementType Add(const ElementType &a,
                         const ElementType &b) const override {
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    const uint64_t bv = Reduce64(static_cast<uint64_t>(b));
    const uint64_t s = av + bv;
    return static_cast<ElementType>(s >= p64_ ? s - p64_ : s);
  }

  inline ElementType Sub(const ElementType &a,
                         const ElementType &b) const override {
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    const uint64_t bv = Reduce64(static_cast<uint64_t>(b));
    return static_cast<ElementType>(av >= bv ? av - bv : av + p64_ - bv);
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    const uint64_t bv = Reduce64(static_cast<uint64_t>(b));
    return static_cast<ElementType>(MulCanonical64(av, bv));
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    const uint64_t bv = Reduce64(static_cast<uint64_t>(b));
    if (bv == 0) {
      throw std::domain_error("Division by zero.");
    }
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    if (av == 0) {
      return static_cast<ElementType>(0);
    }
    return static_cast<ElementType>(MulCanonical64(av, InvCanonical64(bv)));
  }

  inline ElementType Neg(const ElementType &a) const override {
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    return static_cast<ElementType>(av == 0 ? 0 : p64_ - av);
  }

  inline ElementType Inv(const ElementType &a) const override {
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    if (av == 0) {
      throw std::domain_error("Element has no inverse (is zero).");
    }
    return static_cast<ElementType>(InvCanonical64(av));
  }

  inline ElementType Pow(const ElementType &base, uint64_t exp) const override {
    const uint64_t bv = Reduce64(static_cast<uint64_t>(base));
    return static_cast<ElementType>(PowCanonical64(bv, exp));
  }

  virtual bool IsQuadraticResidue(const ElementType &a) const {
    const uint64_t n = Reduce64(static_cast<uint64_t>(a));
    if (n == 0 || p64_ == 2) {
      return true;
    }
    return PowCanonical64(n, (p64_ - 1) >> 1) == 1;
  }

  ElementType Sqrt(const ElementType &a) const override {
    const uint64_t n = Reduce64(static_cast<uint64_t>(a));
    if (n == 0) {
      return static_cast<ElementType>(0);
    }
    if (p64_ == 2) {
      return static_cast<ElementType>(n);
    }

    // Euler's criterion.
    if (PowCanonical64(n, (p64_ - 1) >> 1) != 1) {
      throw std::domain_error("Element is not a quadratic residue.");
    }

    // Fast path for p == 3 (mod 4).
    if ((p64_ & 3ULL) == 3ULL) {
      const uint64_t r = PowCanonical64(n, (p64_ + 1) >> 2);
      return static_cast<ElementType>(std::min(r, p64_ - r));
    }

    // Tonelli-Shanks for p == 1 (mod 4).
    uint64_t q = p64_ - 1;
    unsigned s = static_cast<unsigned>(std::countr_zero(q));
    q >>= s;

    // Find a quadratic non-residue z (reuse cached generator if available).
    uint64_t z =
        static_cast<uint64_t>(generator_.load(std::memory_order_acquire));
    if (z < 2 || PowCanonical64(z, (p64_ - 1) >> 1) != p64_ - 1) {
      z = 2;
      while (PowCanonical64(z, (p64_ - 1) >> 1) != p64_ - 1) {
        ++z;
      }
    }

    unsigned m = s;
    uint64_t c = PowCanonical64(z, q);
    uint64_t t = PowCanonical64(n, q);
    uint64_t r = PowCanonical64(n, (q + 1) >> 1);

    while (t != 1) {
      unsigned i = 1;
      uint64_t tmp = MulCanonical64(t, t);
      while (tmp != 1) {
        tmp = MulCanonical64(tmp, tmp);
        ++i;
        if (i >= m) {
          throw std::domain_error("Element is not a quadratic residue.");
        }
      }

      uint64_t b = c;
      for (unsigned j = 0; j < m - i - 1; ++j) {
        b = MulCanonical64(b, b);
      }
      m = i;
      c = MulCanonical64(b, b);
      t = MulCanonical64(t, c);
      r = MulCanonical64(r, b);
    }

    return static_cast<ElementType>(std::min(r, p64_ - r));
  }

  uint64_t Log(const ElementType &a,
               const ElementType &generator) const override {
    const uint64_t av = Reduce64(static_cast<uint64_t>(a));
    const uint64_t gv = Reduce64(static_cast<uint64_t>(generator));
    if (av == 0) {
      throw std::domain_error("Logarithm of zero is undefined.");
    }
    if (gv == 0) {
      throw std::invalid_argument("Generator for Log cannot be zero.");
    }
    if (p64_ == 2) {
      if (gv != 1) {
        throw std::invalid_argument("Generator for Log must be 1 for GF(2).");
      }
      return 0;
    }
    if (gv == 1) {
      throw std::invalid_argument(
          "Generator for Log cannot be 1 for fields larger than GF(2).");
    }
    if (av == 1) {
      return 0;
    }

    const uint64_t phi = p64_ - 1;
    if (phi <= (uint64_t{1} << 20)) {
      return BabyStepGiantStep(gv, av, phi);
    }
    return PohligHellmanLog(gv, av);
  }

  uint64_t Log(const ElementType &a) const override {
    return Log(a, MultiplicativeGenerator());
  }

  ElementType Random(std::mt19937_64 &rng) const override {
    std::uniform_int_distribution<uint64_t> dist(0, p64_ - 1);
    return static_cast<ElementType>(dist(rng));
  }

  ElementType Random() const override {
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    return Random(rng);
  }

  ElementType MultiplicativeGenerator() const override {
    const ElementType cached = generator_.load(std::memory_order_acquire);
    if (cached != 0) {
      return cached;
    }
    const ElementType gen = ComputeGenerator();
    generator_.store(gen, std::memory_order_release);
    return gen;
  }

  std::vector<ElementType> MultiplicativeGenerators() const override {
    if (p64_ == 2) {
      return {static_cast<ElementType>(1)};
    }

    const uint64_t g = static_cast<uint64_t>(MultiplicativeGenerator());
    const uint64_t phi = p64_ - 1;
    std::vector<ElementType> generators;
    uint64_t cur = 1;
    for (uint64_t k = 1; k <= phi; ++k) {
      cur = MulCanonical64(cur, g);
      if (std::gcd(k, phi) == 1) {
        generators.push_back(static_cast<ElementType>(cur));
      }
    }
    std::sort(generators.begin(), generators.end());
    if (generators.empty()) {
      throw std::runtime_error("No multiplicative generators found.");
    }
    return generators;
  }

  ElementType Frobenius(const ElementType &a) const override {
    return static_cast<ElementType>(Reduce64(static_cast<uint64_t>(a)));
  }

  ElementType Trace(const ElementType &a) const override {
    return static_cast<ElementType>(Reduce64(static_cast<uint64_t>(a)));
  }

  ElementType Norm(const ElementType &a) const override {
    return static_cast<ElementType>(Reduce64(static_cast<uint64_t>(a)));
  }

  void Print(std::ostream &os) const override {
    os << "Galois Prime Field GF(" << p64_ << ")";
    switch (representation_.load(std::memory_order_relaxed)) {
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
      default:
        os << " [Rep: UNKNOWN]";
        break;
    }
  }

  void Print(const ElementType &a, std::ostream &os) const override {
    const uint64_t v = Reduce64(static_cast<uint64_t>(a));
    switch (representation_.load(std::memory_order_relaxed)) {
      case FieldRepresentation::HEX: {
        const std::ios_base::fmtflags saved_flags = os.flags();
        os << "0x" << std::hex << v;
        os.flags(saved_flags);
        break;
      }
      case FieldRepresentation::POW:
        if (v == 0) {
          os << "0";
        } else if (v == 1) {
          os << generator_name_ << "^0";
        } else {
          try {
            const uint64_t power = Log(static_cast<ElementType>(v));
            os << generator_name_ << "^" << power;
          } catch (const std::domain_error &) {
            os << v;
          }
        }
        break;
      case FieldRepresentation::LOG:
        if (v == 0) {
          os << "undefined";
        } else {
          try {
            const uint64_t log_val = Log(static_cast<ElementType>(v));
            os << log_val;
          } catch (const std::domain_error &) {
            os << v;
          }
        }
        break;
      case FieldRepresentation::INT:
      default:
        os << v;
        break;
    }
  }

  std::string ToString(const ElementType &a) const override {
    std::ostringstream oss;
    Print(a, oss);
    return oss.str();
  }

  size_t MemoryBytes() const override {
    return sizeof(*this) + generator_name_.capacity();
  }

  void AddArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("AddArray: span sizes must match.");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    const uint64_t p = p64_;
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
      const uint64_t bv = Reduce64(static_cast<uint64_t>(pb[i]));
      const uint64_t s = av + bv;
      out[i] = static_cast<ElementType>(s >= p ? s - p : s);
    }
  }

  void SubArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("SubArray: span sizes must match.");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    const uint64_t p = p64_;
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
      const uint64_t bv = Reduce64(static_cast<uint64_t>(pb[i]));
      out[i] = static_cast<ElementType>(av >= bv ? av - bv : av + p - bv);
    }
  }

  void MulArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match.");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
      const uint64_t bv = Reduce64(static_cast<uint64_t>(pb[i]));
      out[i] = static_cast<ElementType>(MulCanonical64(av, bv));
    }
  }

  void ScalarMulArray(std::span<ElementType> dst,
                      std::span<const ElementType> a,
                      const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match.");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint64_t sv = Reduce64(static_cast<uint64_t>(scalar));
    if (sv == 0) {
      std::fill(out, out + n, static_cast<ElementType>(0));
      return;
    }
    if (sv == 1) {
      for (size_t i = 0; i < n; ++i) {
        out[i] =
            static_cast<ElementType>(Reduce64(static_cast<uint64_t>(pa[i])));
      }
      return;
    }
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
      out[i] = static_cast<ElementType>(MulCanonical64(av, sv));
    }
  }

  void FmaArray(std::span<ElementType> dst, std::span<const ElementType> a,
                const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match.");
    }
    const uint64_t sv = Reduce64(static_cast<uint64_t>(scalar));
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint64_t p = p64_;
    if (sv == 0) {
      for (size_t i = 0; i < n; ++i) {
        if (static_cast<uint64_t>(out[i]) >= p) {
          out[i] =
              static_cast<ElementType>(Reduce64(static_cast<uint64_t>(out[i])));
        }
      }
      return;
    }
    if (sv == 1) {
      for (size_t i = 0; i < n; ++i) {
        const uint64_t dv = Reduce64(static_cast<uint64_t>(out[i]));
        const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
        const uint64_t s = dv + av;
        out[i] = static_cast<ElementType>(s >= p ? s - p : s);
      }
      return;
    }
    for (size_t i = 0; i < n; ++i) {
      const uint64_t dv = Reduce64(static_cast<uint64_t>(out[i]));
      const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
      const uint64_t prod = MulCanonical64(av, sv);
      const uint64_t s = dv + prod;
      out[i] = static_cast<ElementType>(s >= p ? s - p : s);
    }
  }

  ElementType DotProduct(std::span<const ElementType> a,
                         std::span<const ElementType> b) const override {
    const size_t n = a.size();
    if (b.size() != n) {
      throw std::invalid_argument("DotProduct: span sizes must match.");
    }
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    u128 acc = 0;
    if constexpr (sizeof(ElementType) <= 4) {
      for (size_t i = 0; i < n; ++i) {
        const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
        const uint64_t bv = Reduce64(static_cast<uint64_t>(pb[i]));
        acc += static_cast<u128>(av * bv);
      }
    } else {
      for (size_t i = 0; i < n; ++i) {
        const uint64_t av = Reduce64(static_cast<uint64_t>(pa[i]));
        const uint64_t bv = Reduce64(static_cast<uint64_t>(pb[i]));
        acc += static_cast<u128>(MulCanonical64(av, bv));
      }
    }
    return static_cast<ElementType>(acc % p64_);
  }

 protected:
  static FieldRepresentation ValidateRepresentation(FieldRepresentation rep) {
    if (rep != FieldRepresentation::INT && rep != FieldRepresentation::HEX &&
        rep != FieldRepresentation::POW && rep != FieldRepresentation::LOG) {
      throw std::invalid_argument("Invalid field representation for GF(p).");
    }
    return rep;
  }

  static std::string ValidateGeneratorName(const std::string &name) {
    if (name.empty()) {
      throw std::invalid_argument("Generator name cannot be empty.");
    }
    for (char ch : name) {
      if (ch == '^' || std::isspace(static_cast<unsigned char>(ch)) != 0) {
        throw std::invalid_argument("Invalid generator name: " + name);
      }
    }
    return name;
  }

  // Exact Barrett reduction of any 64-bit value x mod p (2 <= p <= 2^63 - 1).
  inline uint64_t BarrettReduce64(uint64_t x) const noexcept {
    const uint64_t q = static_cast<uint64_t>(
        (static_cast<u128>(x) * static_cast<u128>(barrett_mu_64_)) >> 64);
    const uint64_t rem = x - q * p64_;
    return rem >= p64_ ? rem - p64_ : rem;
  }

  inline uint64_t Reduce64(uint64_t x) const noexcept {
    if (x < p64_) [[likely]] {
      return x;
    }
    return BarrettReduce64(x);
  }

  // Exact division-free Barrett multiplication for canonical 0 <= a, b < p.
  inline uint64_t MulCanonical64(uint64_t a, uint64_t b) const noexcept {
    if constexpr (sizeof(ElementType) <= 4) {
      return BarrettReduce64(a * b);
    } else {
      const u128 x = static_cast<u128>(a) * static_cast<u128>(b);
      const uint64_t q1 = static_cast<uint64_t>(x >> barrett_shift_);
      const uint64_t q = static_cast<uint64_t>(
          (static_cast<u128>(q1) * static_cast<u128>(barrett_mu_128_)) >> 64);
      u128 rem = x - static_cast<u128>(q) * static_cast<u128>(p64_);
      if (rem >= p64_) {
        rem -= p64_;
        if (rem >= p64_) {
          rem -= p64_;
        }
      }
      return static_cast<uint64_t>(rem);
    }
  }

  // Fast Extended Euclidean Algorithm for canonical 0 < a < p <= 2^63 - 1.
  inline uint64_t InvCanonical64(uint64_t a) const {
    uint64_t old_r = p64_;
    uint64_t r = a;
    int64_t old_t = 0;
    int64_t t = 1;

    while (r != 0) {
      const uint64_t q = old_r / r;
      const uint64_t next_r = old_r - q * r;
      old_r = r;
      r = next_r;

      const int64_t next_t = old_t - static_cast<int64_t>(q) * t;
      old_t = t;
      t = next_t;
    }

    if (old_r != 1) {
      throw std::domain_error("Element has no multiplicative inverse.");
    }
    return old_t < 0 ? static_cast<uint64_t>(old_t + static_cast<int64_t>(p64_))
                     : static_cast<uint64_t>(old_t);
  }

  inline uint64_t PowCanonical64(uint64_t base, uint64_t exp) const noexcept {
    if (exp == 0) {
      return 1;
    }
    if (base == 0) {
      return 0;
    }
    if (base == 1) {
      return 1;
    }

    uint64_t result = 1;
    uint64_t cur = base;
    while (exp > 0) {
      if ((exp & 1ULL) != 0) {
        result = MulCanonical64(result, cur);
      }
      exp >>= 1;
      if (exp > 0) {
        cur = MulCanonical64(cur, cur);
      }
    }
    return result;
  }

  ElementType ComputeGenerator() const {
    if (p64_ == 2) {
      return static_cast<ElementType>(1);
    }
    const uint64_t phi = p64_ - 1;
    const std::vector<uint64_t> prime_factors =
        utils::DistinctPrimeFactors(phi);

    for (uint64_t candidate = 2; candidate < p64_; ++candidate) {
      bool is_gen = true;
      for (uint64_t factor : prime_factors) {
        if (PowCanonical64(candidate, phi / factor) == 1) {
          is_gen = false;
          break;
        }
      }
      if (is_gen) {
        return static_cast<ElementType>(candidate);
      }
    }
    throw std::runtime_error("No multiplicative generator found.");
  }

  uint64_t BabyStepGiantStep(uint64_t base, uint64_t target,
                             uint64_t order) const {
    if (target == 1) {
      return 0;
    }
    uint64_t m = utils::ISqrt(order);
    if (m * m < order) {
      ++m;
    }
    constexpr uint64_t kMaxBsgsSteps = uint64_t{1} << 24;
    if (m > kMaxBsgsSteps) {
      throw std::runtime_error(
          "Discrete logarithm factor too large for BSGS memory budget.");
    }

    std::unordered_map<uint64_t, uint64_t> baby_steps;
    baby_steps.reserve(static_cast<size_t>(m));

    uint64_t cur = 1;
    for (uint64_t j = 0; j < m; ++j) {
      baby_steps.try_emplace(cur, j);
      cur = MulCanonical64(cur, base);
    }

    // cur == base^m (mod p).
    const uint64_t giant_factor = InvCanonical64(cur);
    uint64_t gamma = target;
    for (uint64_t i = 0; i < m; ++i) {
      const auto it = baby_steps.find(gamma);
      if (it != baby_steps.end()) {
        const uint64_t ans = i * m + it->second;
        if (ans < order) {
          return ans;
        }
      }
      gamma = MulCanonical64(gamma, giant_factor);
    }

    throw std::domain_error(
        "Element not found in the cyclic group generated by the generator.");
  }

  uint64_t PohligHellmanLog(uint64_t gv, uint64_t av) const {
    const uint64_t ord =
        this->MultiplicativeOrder(static_cast<ElementType>(gv));
    if (PowCanonical64(av, ord) != 1) {
      throw std::domain_error(
          "Element not found in the cyclic group generated by the generator.");
    }

    const auto factors = utils::PrimeFactorizeWithMultiplicity(ord);
    const uint64_t g_inv = InvCanonical64(gv);
    uint64_t crt_rem = 0;
    uint64_t crt_mod = 1;

    for (const auto &[q, e] : factors) {
      const uint64_t gamma = PowCanonical64(gv, ord / q);
      uint64_t x_pe = 0;
      uint64_t q_pow = 1;
      for (uint64_t k = 0; k < e; ++k) {
        const uint64_t rhs = MulCanonical64(av, PowCanonical64(g_inv, x_pe));
        const uint64_t q_next = q_pow * q;
        const uint64_t target = PowCanonical64(rhs, ord / q_next);
        const uint64_t d_k = BabyStepGiantStep(gamma, target, q);
        x_pe += d_k * q_pow;
        q_pow = q_next;
      }

      const uint64_t pe = q_pow;
      if (crt_mod == 1) {
        crt_rem = x_pe;
        crt_mod = pe;
      } else {
        const uint64_t r_mod = crt_rem % pe;
        const uint64_t diff =
            x_pe >= r_mod ? (x_pe - r_mod) : (pe + x_pe - r_mod);
        const uint64_t inv_m = utils::ModInv(crt_mod % pe, pe);
        const uint64_t coeff = static_cast<uint64_t>(
            (static_cast<u128>(diff) * static_cast<u128>(inv_m)) % pe);
        crt_rem += crt_mod * coeff;
        crt_mod *= pe;
      }
    }

    if (PowCanonical64(gv, crt_rem) != av) {
      throw std::domain_error(
          "Element not found in the cyclic group generated by the generator.");
    }
    return crt_rem;
  }

  ElementType p_;
  uint64_t p64_;
  uint64_t barrett_mu_64_;
  uint64_t barrett_mu_128_;
  unsigned barrett_shift_;
  std::atomic<FieldRepresentation> representation_;
  mutable std::atomic<ElementType> generator_;
  std::string generator_name_;
};

template <UnsignedFieldWord ElementType = uint32_t>
class GaloisFieldPrimeTable final : public GaloisFieldPrime<ElementType> {
 public:
  static constexpr uint64_t kMaxTablePrime = uint64_t{1} << 22;

  explicit GaloisFieldPrimeTable(ElementType prime,
                                 const std::string &rep = "int",
                                 bool prime_testing = true,
                                 const std::string &generator_name = "g")
      : GaloisFieldPrime<ElementType>(prime, rep, prime_testing,
                                      generator_name),
        phi_(this->p64_ - 1) {
    if (this->p64_ > kMaxTablePrime) {
      throw std::invalid_argument(
          "Prime exceeds maximum table-backed field size (2^22).");
    }
    InitializeTables();
  }

  inline ElementType Mul(const ElementType &a,
                         const ElementType &b) const override {
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    const uint64_t bv = this->Reduce64(static_cast<uint64_t>(b));
    if (av == 0 || bv == 0) {
      return static_cast<ElementType>(0);
    }
    return alpha_pow_[static_cast<size_t>(log_alpha_[av]) +
                      static_cast<size_t>(log_alpha_[bv])];
  }

  inline ElementType Div(const ElementType &a,
                         const ElementType &b) const override {
    const uint64_t bv = this->Reduce64(static_cast<uint64_t>(b));
    if (bv == 0) {
      throw std::domain_error("Division by zero.");
    }
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    if (av == 0) {
      return static_cast<ElementType>(0);
    }
    return alpha_pow_[static_cast<size_t>(phi_) +
                      static_cast<size_t>(log_alpha_[av]) -
                      static_cast<size_t>(log_alpha_[bv])];
  }

  inline ElementType Inv(const ElementType &a) const override {
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    if (av == 0) {
      throw std::domain_error("Element has no inverse (is zero).");
    }
    return alpha_pow_[static_cast<size_t>(phi_) -
                      static_cast<size_t>(log_alpha_[av])];
  }

  inline ElementType Pow(const ElementType &base, uint64_t exp) const override {
    if (exp == 0) {
      return static_cast<ElementType>(1);
    }
    const uint64_t bv = this->Reduce64(static_cast<uint64_t>(base));
    if (bv == 0) {
      return static_cast<ElementType>(0);
    }
    if (phi_ == 1) {
      return static_cast<ElementType>(1);
    }
    const uint64_t log_b = static_cast<uint64_t>(log_alpha_[bv]);
    const uint64_t idx = (log_b * (exp % phi_)) % phi_;
    return alpha_pow_[static_cast<size_t>(idx)];
  }

  bool IsQuadraticResidue(const ElementType &a) const override {
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    if (av == 0 || this->p64_ == 2) {
      return true;
    }
    return (static_cast<uint64_t>(log_alpha_[av]) & 1ULL) == 0;
  }

  ElementType Sqrt(const ElementType &a) const override {
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    if (av == 0) {
      return static_cast<ElementType>(0);
    }
    if (this->p64_ == 2) {
      return static_cast<ElementType>(1);
    }
    const uint64_t log_a = static_cast<uint64_t>(log_alpha_[av]);
    if ((log_a & 1ULL) != 0) {
      throw std::domain_error("Element is not a quadratic residue.");
    }
    const uint64_t r =
        static_cast<uint64_t>(alpha_pow_[static_cast<size_t>(log_a >> 1)]);
    return static_cast<ElementType>(std::min(r, this->p64_ - r));
  }

  uint64_t Log(const ElementType &a,
               const ElementType &generator) const override {
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    const uint64_t gv = this->Reduce64(static_cast<uint64_t>(generator));
    if (av == 0) {
      throw std::domain_error("Logarithm of zero is undefined.");
    }
    if (gv == 0) {
      throw std::invalid_argument("Generator for Log cannot be zero.");
    }
    if (this->p64_ == 2) {
      if (gv != 1) {
        throw std::invalid_argument("Generator for Log must be 1 for GF(2).");
      }
      return 0;
    }
    if (gv == 1) {
      throw std::invalid_argument(
          "Generator for Log cannot be 1 for fields larger than GF(2).");
    }
    if (av == 1) {
      return 0;
    }

    const uint64_t log_a = static_cast<uint64_t>(log_alpha_[av]);
    const uint64_t log_g = static_cast<uint64_t>(log_alpha_[gv]);
    if (log_g == 1) {
      return log_a;
    }

    // Solve x * log_g == log_a (mod phi_).
    const uint64_t d = std::gcd(log_g, phi_);
    if (log_a % d != 0) {
      throw std::domain_error(
          "Element not found in the cyclic group generated by the generator.");
    }
    const uint64_t mod = phi_ / d;
    const uint64_t inv_g = utils::ModInv(log_g / d, mod);
    return ((log_a / d) * inv_g) % mod;
  }

  uint64_t Log(const ElementType &a) const override {
    const uint64_t av = this->Reduce64(static_cast<uint64_t>(a));
    if (av == 0) {
      throw std::domain_error("Logarithm of zero is undefined.");
    }
    return static_cast<uint64_t>(log_alpha_[av]);
  }

  size_t MemoryBytes() const override {
    return sizeof(*this) + this->generator_name_.capacity() +
           alpha_pow_.capacity() * sizeof(ElementType) +
           log_alpha_.capacity() * sizeof(ElementType);
  }

  void AddArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    GaloisFieldPrime<ElementType>::AddArray(dst, a, b);
  }

  void SubArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    GaloisFieldPrime<ElementType>::SubArray(dst, a, b);
  }

  void MulArray(std::span<ElementType> dst, std::span<const ElementType> a,
                std::span<const ElementType> b) const override {
    const size_t n = dst.size();
    if (a.size() != n || b.size() != n) {
      throw std::invalid_argument("MulArray: span sizes must match.");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const ElementType *pb = b.data();
    const ElementType *exp_tbl = alpha_pow_.data();
    const ElementType *log_tbl = log_alpha_.data();
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = this->Reduce64(static_cast<uint64_t>(pa[i]));
      const uint64_t bv = this->Reduce64(static_cast<uint64_t>(pb[i]));
      out[i] = (av == 0 || bv == 0) ? static_cast<ElementType>(0)
                                    : exp_tbl[static_cast<size_t>(log_tbl[av]) +
                                              static_cast<size_t>(log_tbl[bv])];
    }
  }

  void ScalarMulArray(std::span<ElementType> dst,
                      std::span<const ElementType> a,
                      const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("ScalarMulArray: span sizes must match.");
    }
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint64_t sv = this->Reduce64(static_cast<uint64_t>(scalar));
    if (sv == 0) {
      std::fill(out, out + n, static_cast<ElementType>(0));
      return;
    }
    if (sv == 1) {
      for (size_t i = 0; i < n; ++i) {
        out[i] = static_cast<ElementType>(
            this->Reduce64(static_cast<uint64_t>(pa[i])));
      }
      return;
    }
    const ElementType *shifted =
        alpha_pow_.data() + static_cast<size_t>(log_alpha_[sv]);
    const ElementType *log_tbl = log_alpha_.data();
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = this->Reduce64(static_cast<uint64_t>(pa[i]));
      out[i] = av == 0 ? static_cast<ElementType>(0) : shifted[log_tbl[av]];
    }
  }

  void FmaArray(std::span<ElementType> dst, std::span<const ElementType> a,
                const ElementType &scalar) const override {
    const size_t n = dst.size();
    if (a.size() != n) {
      throw std::invalid_argument("FmaArray: span sizes must match.");
    }
    const uint64_t sv = this->Reduce64(static_cast<uint64_t>(scalar));
    ElementType *out = dst.data();
    const ElementType *pa = a.data();
    const uint64_t p = this->p64_;
    if (sv == 0) {
      for (size_t i = 0; i < n; ++i) {
        if (static_cast<uint64_t>(out[i]) >= p) {
          out[i] = static_cast<ElementType>(
              this->Reduce64(static_cast<uint64_t>(out[i])));
        }
      }
      return;
    }
    if (sv == 1) {
      for (size_t i = 0; i < n; ++i) {
        const uint64_t dv = this->Reduce64(static_cast<uint64_t>(out[i]));
        const uint64_t av = this->Reduce64(static_cast<uint64_t>(pa[i]));
        const uint64_t s = dv + av;
        out[i] = static_cast<ElementType>(s >= p ? s - p : s);
      }
      return;
    }
    const ElementType *shifted =
        alpha_pow_.data() + static_cast<size_t>(log_alpha_[sv]);
    const ElementType *log_tbl = log_alpha_.data();
    for (size_t i = 0; i < n; ++i) {
      const uint64_t av = this->Reduce64(static_cast<uint64_t>(pa[i]));
      const uint64_t dv = this->Reduce64(static_cast<uint64_t>(out[i]));
      if (av == 0) {
        out[i] = static_cast<ElementType>(dv);
      } else {
        const uint64_t prod = static_cast<uint64_t>(shifted[log_tbl[av]]);
        const uint64_t s = dv + prod;
        out[i] = static_cast<ElementType>(s >= p ? s - p : s);
      }
    }
  }

 private:
  void InitializeTables() {
    const size_t phi = static_cast<size_t>(phi_);
    alpha_pow_.resize(2 * phi);
    log_alpha_.assign(static_cast<size_t>(this->p64_),
                      static_cast<ElementType>(0));

    if (this->p64_ == 2) {
      alpha_pow_[0] = static_cast<ElementType>(1);
      alpha_pow_[1] = static_cast<ElementType>(1);
      log_alpha_[1] = static_cast<ElementType>(0);
      return;
    }

    const uint64_t gen = static_cast<uint64_t>(this->MultiplicativeGenerator());
    uint64_t cur = 1;
    for (size_t i = 0; i < phi; ++i) {
      const ElementType val = static_cast<ElementType>(cur);
      alpha_pow_[i] = val;
      alpha_pow_[i + phi] = val;
      log_alpha_[cur] = static_cast<ElementType>(i);
      cur = this->MulCanonical64(cur, gen);
    }
  }

  uint64_t phi_;
  std::vector<ElementType> alpha_pow_;
  std::vector<ElementType> log_alpha_;
};

template <UnsignedFieldWord ElementType = uint32_t>
using GFP = GaloisFieldPrime<ElementType>;
template <UnsignedFieldWord ElementType = uint32_t>
using GFPLOG = GaloisFieldPrimeTable<ElementType>;

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_PRIME_HPP_
