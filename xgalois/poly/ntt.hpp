#ifndef XGALOIS_POLY_NTT_HPP_
#define XGALOIS_POLY_NTT_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/simd/kernels_scalar.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace poly {

namespace detail {

__extension__ using Uint128 = unsigned __int128;

constexpr uint64_t kNttMaxPrime = (uint64_t{1} << 62) - 1ULL;
constexpr std::size_t kDefaultNttPrecomputeSize = 4096;

inline constexpr uint64_t AddMod64(uint64_t a, uint64_t b,
                                   uint64_t mod) noexcept {
  const uint64_t sum = a + b;
  return (sum >= mod) ? (sum - mod) : sum;
}

inline constexpr uint64_t SubMod64(uint64_t a, uint64_t b,
                                   uint64_t mod) noexcept {
  return (a >= b) ? (a - b) : (a + mod - b);
}

inline constexpr bool IsPowerOfTwo(std::size_t n) noexcept {
  return n != 0 && (n & (n - 1)) == 0;
}

inline std::size_t NextPowerOfTwo(std::size_t n) {
  if (n <= 1) {
    return 1;
  }
  const unsigned bits = static_cast<unsigned>(std::bit_width(n - 1));
  if (bits >= sizeof(std::size_t) * 8) {
    throw std::overflow_error("NextPowerOfTwo: size overflow");
  }
  return std::size_t{1} << bits;
}

// Finds the smallest primitive root modulo an odd prime p >= 3.
inline uint64_t FindPrimitiveRoot(uint64_t prime) {
  if (prime == 2) {
    return 1;
  }
  const uint64_t phi = prime - 1;
  const std::vector<uint64_t> factors = utils::DistinctPrimeFactors(phi);
  for (uint64_t cand = 2; cand < prime; ++cand) {
    bool is_primitive = true;
    for (uint64_t q : factors) {
      if (utils::ModPow(cand, phi / q, prime) == 1) {
        is_primitive = false;
        break;
      }
    }
    if (is_primitive) {
      return cand;
    }
  }
  throw std::runtime_error("FindPrimitiveRoot: failed to find primitive root");
}

struct NttRootTables {
  std::size_t capacity = 0;
  std::vector<uint64_t> fwd_roots;
  std::vector<uint64_t> fwd_roots_shoup;
  std::vector<uint64_t> inv_roots;
  std::vector<uint64_t> inv_roots_shoup;
};

// Builds flat stage root tables up to transform size `cap` (power of 2).
// For stage half-size `half` in {1, 2, 4, ..., cap/2}:
//   fwd_roots[half + j] = w_{2*half}^j     for j in [0, half)
//   inv_roots[half + j] = w_{2*half}^{-j}  for j in [0, half)
inline NttRootTables BuildNttRootTables(uint64_t prime, uint64_t primitive_root,
                                        std::size_t cap) {
  NttRootTables tables;
  if (cap < 2) {
    tables.capacity = 1;
    return tables;
  }
  tables.capacity = cap;
  tables.fwd_roots.resize(cap, 1ULL);
  tables.fwd_roots_shoup.resize(cap, 0ULL);
  tables.inv_roots.resize(cap, 1ULL);
  tables.inv_roots_shoup.resize(cap, 0ULL);

  for (std::size_t half = 1; half < cap; half <<= 1) {
    const std::size_t len = half << 1;
    const uint64_t step_fwd =
        utils::ModPow(primitive_root, (prime - 1) / len, prime);
    const uint64_t step_inv = utils::ModInv(step_fwd, prime);
    const uint64_t step_fwd_shoup = simd::ComputeShoupFactor(step_fwd, prime);
    const uint64_t step_inv_shoup = simd::ComputeShoupFactor(step_inv, prime);

    uint64_t wf = 1ULL;
    uint64_t wi = 1ULL;
    for (std::size_t j = 0; j < half; ++j) {
      tables.fwd_roots[half + j] = wf;
      tables.fwd_roots_shoup[half + j] = simd::ComputeShoupFactor(wf, prime);
      tables.inv_roots[half + j] = wi;
      tables.inv_roots_shoup[half + j] = simd::ComputeShoupFactor(wi, prime);
      wf = simd::ShoupMulMod(wf, step_fwd, step_fwd_shoup, prime);
      wi = simd::ShoupMulMod(wi, step_inv, step_inv_shoup, prime);
    }
  }
  return tables;
}

// Bit-reversal-free DIF (Gentleman-Sande) forward NTT:
// natural-order input -> bit-reversed-order evaluations.
inline void ForwardDifInPlace(std::span<uint64_t> data, uint64_t prime,
                              const NttRootTables& tables) noexcept {
  const std::size_t n = data.size();
  for (std::size_t half = n >> 1; half >= 1; half >>= 1) {
    const std::size_t len = half << 1;
    const uint64_t* roots = tables.fwd_roots.data() + half;
    const uint64_t* roots_shoup = tables.fwd_roots_shoup.data() + half;
    for (std::size_t base = 0; base < n; base += len) {
      uint64_t* lo = data.data() + base;
      uint64_t* hi = lo + half;
      for (std::size_t j = 0; j < half; ++j) {
        const uint64_t u = lo[j];
        const uint64_t v = hi[j];
        lo[j] = AddMod64(u, v, prime);
        const uint64_t diff = SubMod64(u, v, prime);
        hi[j] = simd::ShoupMulMod(diff, roots[j], roots_shoup[j], prime);
      }
    }
  }
}

// Bit-reversal-free DIT (Cooley-Tukey) inverse NTT:
// bit-reversed-order evaluations -> natural-order coefficients (scaled by 1/n).
inline void InverseDitInPlace(std::span<uint64_t> data, uint64_t prime,
                              const NttRootTables& tables) {
  const std::size_t n = data.size();
  for (std::size_t half = 1; half < n; half <<= 1) {
    const std::size_t len = half << 1;
    const uint64_t* roots = tables.inv_roots.data() + half;
    const uint64_t* roots_shoup = tables.inv_roots_shoup.data() + half;
    for (std::size_t base = 0; base < n; base += len) {
      uint64_t* lo = data.data() + base;
      uint64_t* hi = lo + half;
      for (std::size_t j = 0; j < half; ++j) {
        const uint64_t u = lo[j];
        const uint64_t v =
            simd::ShoupMulMod(hi[j], roots[j], roots_shoup[j], prime);
        lo[j] = AddMod64(u, v, prime);
        hi[j] = SubMod64(u, v, prime);
      }
    }
  }
  const uint64_t inv_n = utils::ModInv(static_cast<uint64_t>(n), prime);
  const uint64_t inv_n_shoup = simd::ComputeShoupFactor(inv_n, prime);
  simd::GfPrimeMulShoupArray(data, std::span<const uint64_t>(data), inv_n,
                             inv_n_shoup, prime);
}

inline void TrimTrailingZeros(std::vector<uint64_t>& vec) noexcept {
  while (!vec.empty() && vec.back() == 0ULL) {
    vec.pop_back();
  }
}

}  // namespace detail

// Number Theoretic Transform plan for any NTT-friendly prime p < 2^62
// where 2^k | (p - 1). Uses Shoup preconditioned modular multiplication and
// a bit-reversal-free Gentleman-Sande (DIF) / Cooley-Tukey (DIT) pair.
class NttPlan {
 public:
  explicit NttPlan(uint64_t prime, std::size_t max_size = 0) : prime_(prime) {
    if (prime_ < 3ULL || prime_ > detail::kNttMaxPrime) {
      throw std::invalid_argument(
          "NttPlan: prime must satisfy 3 <= prime < 2^62");
    }
    if (!utils::IsPrime(prime_)) {
      throw std::invalid_argument("NttPlan: modulus must be prime");
    }
    max_two_adicity_ = static_cast<std::size_t>(std::countr_zero(prime_ - 1));
    if (max_two_adicity_ == 0) {
      throw std::invalid_argument(
          "NttPlan: prime - 1 must be divisible by at least 2");
    }
    primitive_root_ = detail::FindPrimitiveRoot(prime_);

    const std::size_t max_cap = MaxTransformSize();
    std::size_t init_cap =
        (max_size == 0) ? std::min(detail::kDefaultNttPrecomputeSize, max_cap)
                        : detail::NextPowerOfTwo(max_size);
    if (init_cap > max_cap) {
      throw std::invalid_argument(
          "NttPlan: requested max_size exceeds prime's 2-adic capacity");
    }
    tables_ = detail::BuildNttRootTables(prime_, primitive_root_, init_cap);
  }

  uint64_t Prime() const noexcept { return prime_; }

  uint64_t PrimitiveRoot() const noexcept { return primitive_root_; }

  std::size_t MaxTwoAdicity() const noexcept { return max_two_adicity_; }

  std::size_t MaxTransformSize() const noexcept {
    constexpr std::size_t kMaxShift = sizeof(std::size_t) * 8 - 1;
    const std::size_t shift = std::min(max_two_adicity_, kMaxShift);
    return std::size_t{1} << shift;
  }

  std::size_t ReservedSize() const noexcept { return tables_.capacity; }

  void Reserve(std::size_t max_size) {
    const std::size_t cap = detail::NextPowerOfTwo(max_size);
    if (cap > MaxTransformSize()) {
      throw std::invalid_argument(
          "NttPlan::Reserve: size exceeds prime's 2-adic capacity");
    }
    if (cap > tables_.capacity) {
      tables_ = detail::BuildNttRootTables(prime_, primitive_root_, cap);
    }
  }

  // In-place forward DIF NTT (natural order -> bit-reversed order).
  // Expects data.size() to be a power of 2 <= MaxTransformSize() and
  // coefficients in [0, prime - 1].
  void Forward(std::span<uint64_t> data) const {
    const std::size_t n = data.size();
    if (n <= 1) {
      if (n == 1 && data[0] >= prime_) {
        data[0] %= prime_;
      }
      return;
    }
    ValidateTransformSize(n);
    if (n <= tables_.capacity) {
      detail::ForwardDifInPlace(data, prime_, tables_);
    } else {
      const detail::NttRootTables local =
          detail::BuildNttRootTables(prime_, primitive_root_, n);
      detail::ForwardDifInPlace(data, prime_, local);
    }
  }

  // In-place inverse DIT NTT (bit-reversed order -> natural order, scaled by
  // n^-1 mod prime). Exact inverse of Forward().
  void Inverse(std::span<uint64_t> data) const {
    const std::size_t n = data.size();
    if (n <= 1) {
      if (n == 1 && data[0] >= prime_) {
        data[0] %= prime_;
      }
      return;
    }
    ValidateTransformSize(n);
    if (n <= tables_.capacity) {
      detail::InverseDitInPlace(data, prime_, tables_);
    } else {
      const detail::NttRootTables local =
          detail::BuildNttRootTables(prime_, primitive_root_, n);
      detail::InverseDitInPlace(data, prime_, local);
    }
  }

  // Computes exact polynomial convolution (a * b) mod prime, trimmed of
  // trailing zeros.
  std::vector<uint64_t> Multiply(std::span<const uint64_t> a,
                                 std::span<const uint64_t> b) const {
    if (a.empty() || b.empty()) {
      return {};
    }
    std::size_t a_len = a.size();
    while (a_len > 0 && (a[a_len - 1] % prime_) == 0ULL) {
      --a_len;
    }
    std::size_t b_len = b.size();
    while (b_len > 0 && (b[b_len - 1] % prime_) == 0ULL) {
      --b_len;
    }
    if (a_len == 0 || b_len == 0) {
      return {};
    }

    const std::size_t out_len = a_len + b_len - 1;
    const std::size_t n = detail::NextPowerOfTwo(out_len);
    ValidateTransformSize(n);

    std::vector<uint64_t> fa(n, 0ULL);
    std::vector<uint64_t> fb(n, 0ULL);
    for (std::size_t i = 0; i < a_len; ++i) {
      fa[i] = (a[i] >= prime_) ? (a[i] % prime_) : a[i];
    }
    for (std::size_t i = 0; i < b_len; ++i) {
      fb[i] = (b[i] >= prime_) ? (b[i] % prime_) : b[i];
    }

    const detail::NttRootTables* active = &tables_;
    detail::NttRootTables local;
    if (n > tables_.capacity) {
      local = detail::BuildNttRootTables(prime_, primitive_root_, n);
      active = &local;
    }

    if (n > 1) {
      detail::ForwardDifInPlace(fa, prime_, *active);
      detail::ForwardDifInPlace(fb, prime_, *active);
    }
    for (std::size_t i = 0; i < n; ++i) {
      fa[i] = static_cast<uint64_t>(
          (static_cast<detail::Uint128>(fa[i]) * fb[i]) % prime_);
    }
    if (n > 1) {
      detail::InverseDitInPlace(fa, prime_, *active);
    }

    fa.resize(out_len);
    detail::TrimTrailingZeros(fa);
    return fa;
  }

  // Convenience overload multiplying two PolynomialDense objects over any
  // prime FiniteField with Order() == Prime().
  template <FiniteField PrimeField>
  PolynomialDense<PrimeField> Multiply(
      const PolynomialDense<PrimeField>& a,
      const PolynomialDense<PrimeField>& b) const {
    using RawType = typename PrimeField::element_type;
    static_assert(std::is_integral_v<RawType>,
                  "NttPlan::Multiply requires an integral prime field element");
    if (!a.HasField() || !b.HasField() || a.Field() != b.Field()) {
      throw std::invalid_argument(
          "NttPlan::Multiply: polynomials must share a non-null field");
    }
    const auto field = a.Field();
    if (field->Degree() != 1 || field->Order() != prime_) {
      throw std::invalid_argument(
          "NttPlan::Multiply: field must be a prime field matching NttPlan "
          "prime");
    }
    if (a.IsZero() || b.IsZero()) {
      return PolynomialDense<PrimeField>::Zero(field, a.GetVariable());
    }

    std::vector<uint64_t> a_u64(a.Size());
    for (std::size_t i = 0; i < a.Size(); ++i) {
      a_u64[i] = static_cast<uint64_t>(field->GetElementValue(a.RawCoeff(i)));
    }
    std::vector<uint64_t> b_u64(b.Size());
    for (std::size_t i = 0; i < b.Size(); ++i) {
      b_u64[i] = static_cast<uint64_t>(field->GetElementValue(b.RawCoeff(i)));
    }

    std::vector<uint64_t> prod_u64 = Multiply(std::span<const uint64_t>(a_u64),
                                              std::span<const uint64_t>(b_u64));
    std::vector<RawType> prod_raw(prod_u64.size());
    for (std::size_t i = 0; i < prod_u64.size(); ++i) {
      prod_raw[i] = field->SetElementValue(static_cast<RawType>(prod_u64[i]));
    }
    return PolynomialDense<PrimeField>::FromRaw(std::move(prod_raw), field,
                                                a.GetVariable());
  }

 private:
  void ValidateTransformSize(std::size_t n) const {
    if (!detail::IsPowerOfTwo(n)) {
      throw std::invalid_argument(
          "NttPlan: transform size must be a power of 2");
    }
    if (n > MaxTransformSize()) {
      throw std::invalid_argument(
          "NttPlan: transform size exceeds prime's 2-adic limit");
    }
  }

  uint64_t prime_;
  uint64_t primitive_root_;
  std::size_t max_two_adicity_;
  detail::NttRootTables tables_;
};

// 3-prime Garner CRT NTT for O(N log N) polynomial multiplication over
// arbitrary prime fields (and arbitrary moduli p < 2^31 for lengths up to
// 2^21; product P1 * P2 * P3 > 4.7 * 10^26 > 2^88).
class CrtNttPlan {
 public:
  static constexpr uint64_t kPrime1 = 998244353ULL;   // 119 * 2^23 + 1
  static constexpr uint64_t kPrime2 = 1004535809ULL;  // 479 * 2^21 + 1
  static constexpr uint64_t kPrime3 = 469762049ULL;   // 7   * 2^26 + 1
  static constexpr std::size_t kMaxTransformSize = std::size_t{1} << 21;

  explicit CrtNttPlan(uint64_t modulus, std::size_t max_size = 0)
      : modulus_(modulus),
        ntt1_(kPrime1, max_size),
        ntt2_(kPrime2, max_size),
        ntt3_(kPrime3, max_size) {
    if (modulus_ < 2ULL) {
      throw std::invalid_argument("CrtNttPlan: modulus must be >= 2");
    }
    InitGarnerConstants();
  }

  uint64_t Modulus() const noexcept { return modulus_; }

  void Reserve(std::size_t max_size) {
    if (max_size > kMaxTransformSize) {
      throw std::invalid_argument(
          "CrtNttPlan::Reserve: size exceeds 2^21 CRT-NTT limit");
    }
    ntt1_.Reserve(max_size);
    ntt2_.Reserve(max_size);
    ntt3_.Reserve(max_size);
  }

  // Computes exact polynomial convolution (a * b) mod modulus_, trimmed of
  // trailing zeros. Requires min(|a|, |b|) * (modulus - 1)^2 < P1 * P2 * P3
  // (always satisfied for modulus < 2^31 and transform size <= 2^21).
  std::vector<uint64_t> Multiply(std::span<const uint64_t> a,
                                 std::span<const uint64_t> b) const {
    if (a.empty() || b.empty()) {
      return {};
    }
    std::size_t a_len = a.size();
    while (a_len > 0 && (a[a_len - 1] % modulus_) == 0ULL) {
      --a_len;
    }
    std::size_t b_len = b.size();
    while (b_len > 0 && (b[b_len - 1] % modulus_) == 0ULL) {
      --b_len;
    }
    if (a_len == 0 || b_len == 0) {
      return {};
    }

    const std::size_t out_len = a_len + b_len - 1;
    const std::size_t n = detail::NextPowerOfTwo(out_len);
    if (n > kMaxTransformSize) {
      throw std::invalid_argument(
          "CrtNttPlan::Multiply: convolution size exceeds 2^21");
    }

    // Check exactness bound: min(a_len, b_len) * (modulus - 1)^2 < P1*P2*P3.
    const detail::Uint128 max_coeff =
        static_cast<detail::Uint128>(modulus_ - 1);
    const detail::Uint128 min_len =
        static_cast<detail::Uint128>(std::min(a_len, b_len));
    const detail::Uint128 p123 = static_cast<detail::Uint128>(kPrime1) *
                                 static_cast<detail::Uint128>(kPrime2) *
                                 static_cast<detail::Uint128>(kPrime3);
    if (max_coeff != 0 && (max_coeff > p123 / max_coeff ||
                           min_len > p123 / (max_coeff * max_coeff))) {
      throw std::overflow_error(
          "CrtNttPlan::Multiply: exact integer convolution exceeds 3-prime "
          "CRT dynamic range");
    }

    std::vector<uint64_t> a_red(a_len);
    std::vector<uint64_t> b_red(b_len);
    for (std::size_t i = 0; i < a_len; ++i) {
      a_red[i] = (a[i] >= modulus_) ? (a[i] % modulus_) : a[i];
    }
    for (std::size_t i = 0; i < b_len; ++i) {
      b_red[i] = (b[i] >= modulus_) ? (b[i] % modulus_) : b[i];
    }

    auto c1 = ntt1_.Multiply(a_red, b_red);
    auto c2 = ntt2_.Multiply(a_red, b_red);
    auto c3 = ntt3_.Multiply(a_red, b_red);

    std::vector<uint64_t> out(out_len, 0ULL);
    for (std::size_t i = 0; i < out_len; ++i) {
      const uint64_t r1 = (i < c1.size()) ? c1[i] : 0ULL;
      const uint64_t r2 = (i < c2.size()) ? c2[i] : 0ULL;
      const uint64_t r3 = (i < c3.size()) ? c3[i] : 0ULL;

      // Garner's algorithm: x = t1 + t2 * P1 + t3 * P1 * P2 < P1 * P2 * P3.
      const uint64_t t1 = r1;
      const uint64_t diff2 = detail::SubMod64(r2, t1 % kPrime2, kPrime2);
      const uint64_t t2 = static_cast<uint64_t>(
          (static_cast<detail::Uint128>(diff2) * inv_p1_mod_p2_) % kPrime2);
      const uint64_t v12_mod_p3 =
          static_cast<uint64_t>((static_cast<detail::Uint128>(t1) +
                                 static_cast<detail::Uint128>(t2) * kPrime1) %
                                kPrime3);
      const uint64_t diff3 = detail::SubMod64(r3, v12_mod_p3, kPrime3);
      const uint64_t t3 = static_cast<uint64_t>(
          (static_cast<detail::Uint128>(diff3) * inv_p1p2_mod_p3_) % kPrime3);

      const detail::Uint128 combined =
          static_cast<detail::Uint128>(t1 % modulus_) +
          static_cast<detail::Uint128>(t2 % modulus_) * p1_mod_m_ +
          static_cast<detail::Uint128>(t3 % modulus_) * p1p2_mod_m_;
      out[i] = static_cast<uint64_t>(combined % modulus_);
    }

    detail::TrimTrailingZeros(out);
    return out;
  }

  template <FiniteField PrimeField>
  PolynomialDense<PrimeField> Multiply(
      const PolynomialDense<PrimeField>& a,
      const PolynomialDense<PrimeField>& b) const {
    using RawType = typename PrimeField::element_type;
    static_assert(
        std::is_integral_v<RawType>,
        "CrtNttPlan::Multiply requires an integral prime field element");
    if (!a.HasField() || !b.HasField() || a.Field() != b.Field()) {
      throw std::invalid_argument(
          "CrtNttPlan::Multiply: polynomials must share a non-null field");
    }
    const auto field = a.Field();
    if (field->Degree() != 1 || field->Order() != modulus_) {
      throw std::invalid_argument(
          "CrtNttPlan::Multiply: field must be a prime field matching "
          "CrtNttPlan modulus");
    }
    if (a.IsZero() || b.IsZero()) {
      return PolynomialDense<PrimeField>::Zero(field, a.GetVariable());
    }

    std::vector<uint64_t> a_u64(a.Size());
    for (std::size_t i = 0; i < a.Size(); ++i) {
      a_u64[i] = static_cast<uint64_t>(field->GetElementValue(a.RawCoeff(i)));
    }
    std::vector<uint64_t> b_u64(b.Size());
    for (std::size_t i = 0; i < b.Size(); ++i) {
      b_u64[i] = static_cast<uint64_t>(field->GetElementValue(b.RawCoeff(i)));
    }

    std::vector<uint64_t> prod_u64 = Multiply(std::span<const uint64_t>(a_u64),
                                              std::span<const uint64_t>(b_u64));
    std::vector<RawType> prod_raw(prod_u64.size());
    for (std::size_t i = 0; i < prod_u64.size(); ++i) {
      prod_raw[i] = field->SetElementValue(static_cast<RawType>(prod_u64[i]));
    }
    return PolynomialDense<PrimeField>::FromRaw(std::move(prod_raw), field,
                                                a.GetVariable());
  }

 private:
  void InitGarnerConstants() {
    inv_p1_mod_p2_ = utils::ModInv(kPrime1 % kPrime2, kPrime2);
    const uint64_t p1p2_mod_p3 = static_cast<uint64_t>(
        (static_cast<detail::Uint128>(kPrime1) * kPrime2) % kPrime3);
    inv_p1p2_mod_p3_ = utils::ModInv(p1p2_mod_p3, kPrime3);
    p1_mod_m_ = kPrime1 % modulus_;
    p1p2_mod_m_ = static_cast<uint64_t>(
        (static_cast<detail::Uint128>(kPrime1) * kPrime2) % modulus_);
  }

  uint64_t modulus_;
  NttPlan ntt1_;
  NttPlan ntt2_;
  NttPlan ntt3_;
  uint64_t inv_p1_mod_p2_ = 0;
  uint64_t inv_p1p2_mod_p3_ = 0;
  uint64_t p1_mod_m_ = 0;
  uint64_t p1p2_mod_m_ = 0;
};

// Arbitrary-length N-point DFT / IDFT over a prime field F_p containing a
// primitive N-th root of unity omega (N | (p - 1)).
// Uses the triangular-number Chirp-Z identity:
//   j * k = C(j + k, 2) - C(j, 2) - C(k, 2)
// so only integer powers of omega are needed (no 2N-th root required!).
class BluesteinDft {
 public:
  BluesteinDft(uint64_t prime, std::size_t length, uint64_t omega = 0)
      : prime_(prime), length_(length), crt_plan_(prime) {
    if (prime_ < 2ULL || !utils::IsPrime(prime_)) {
      throw std::invalid_argument("BluesteinDft: modulus must be prime");
    }
    if (length_ == 0) {
      throw std::invalid_argument("BluesteinDft: length must be >= 1");
    }
    if ((prime_ - 1) % length_ != 0) {
      throw std::invalid_argument(
          "BluesteinDft: length must divide prime - 1 for an N-th root of "
          "unity to exist in F_p");
    }

    if (omega == 0) {
      if (length_ == 1) {
        omega_ = 1ULL;
      } else {
        const uint64_t g = detail::FindPrimitiveRoot(prime_);
        omega_ = utils::ModPow(g, (prime_ - 1) / length_, prime_);
      }
    } else {
      if (omega >= prime_) {
        throw std::invalid_argument("BluesteinDft: omega must be < prime");
      }
      if (utils::ModPow(omega, length_, prime_) != 1ULL) {
        throw std::invalid_argument("BluesteinDft: omega^N != 1 mod prime");
      }
      if (length_ > 1) {
        for (uint64_t q : utils::DistinctPrimeFactors(length_)) {
          if (utils::ModPow(omega, length_ / q, prime_) == 1ULL) {
            throw std::invalid_argument(
                "BluesteinDft: omega is not a primitive N-th root of unity");
          }
        }
      }
      omega_ = omega;
    }

    inv_omega_ = utils::ModInv(omega_, prime_);
    inv_length_ =
        utils::ModInv(static_cast<uint64_t>(length_) % prime_, prime_);

    // Precompute chirp tables using triangular exponents C(m, 2) mod N.
    fwd_chirp_neg_.resize(length_);
    inv_chirp_neg_.resize(length_);
    for (std::size_t k = 0; k < length_; ++k) {
      const uint64_t tri = static_cast<uint64_t>(
          (static_cast<detail::Uint128>(k) * (k - 1) / 2) % length_);
      fwd_chirp_neg_[k] = utils::ModPow(inv_omega_, tri, prime_);
      inv_chirp_neg_[k] = utils::ModPow(omega_, tri, prime_);
    }

    const std::size_t kernel_len = 2 * length_ - 1;
    fwd_kernel_.resize(kernel_len);
    inv_kernel_.resize(kernel_len);
    for (std::size_t m = 0; m < kernel_len; ++m) {
      const uint64_t tri = static_cast<uint64_t>(
          (static_cast<detail::Uint128>(m) * (m - 1) / 2) % length_);
      fwd_kernel_[m] = utils::ModPow(omega_, tri, prime_);
      inv_kernel_[m] = utils::ModPow(inv_omega_, tri, prime_);
    }

    // Use a direct single-prime NttPlan when prime_ is NTT-friendly with
    // enough 2-adicity for the Bluestein convolution of length 3*N - 2.
    const std::size_t conv_size =
        (length_ <= 1) ? 1 : detail::NextPowerOfTwo(3 * length_ - 2);
    if (prime_ >= 3ULL && prime_ <= detail::kNttMaxPrime) {
      const std::size_t two_adicity =
          static_cast<std::size_t>(std::countr_zero(prime_ - 1));
      if (two_adicity < sizeof(std::size_t) * 8 &&
          (std::size_t{1} << two_adicity) >= conv_size) {
        direct_ntt_.emplace_back(prime_, conv_size);
      }
    }
    if (direct_ntt_.empty() && conv_size <= CrtNttPlan::kMaxTransformSize) {
      crt_plan_.Reserve(conv_size);
    }
  }

  uint64_t Prime() const noexcept { return prime_; }
  std::size_t Length() const noexcept { return length_; }
  uint64_t Root() const noexcept { return omega_; }

  // Computes X[k] = sum_{j=0}^{N-1} x[j] * omega^(j*k) mod prime.
  std::vector<uint64_t> Forward(std::span<const uint64_t> input) const {
    return Transform(input, fwd_chirp_neg_, fwd_kernel_, false);
  }

  // Computes x[k] = N^-1 * sum_{j=0}^{N-1} X[j] * omega^(-j*k) mod prime.
  std::vector<uint64_t> Inverse(std::span<const uint64_t> input) const {
    return Transform(input, inv_chirp_neg_, inv_kernel_, true);
  }

  void ForwardInPlace(std::span<uint64_t> data) const {
    auto res = Forward(std::span<const uint64_t>(data));
    std::copy(res.begin(), res.end(), data.begin());
  }

  void InverseInPlace(std::span<uint64_t> data) const {
    auto res = Inverse(std::span<const uint64_t>(data));
    std::copy(res.begin(), res.end(), data.begin());
  }

 private:
  std::vector<uint64_t> Convolve(std::span<const uint64_t> a,
                                 std::span<const uint64_t> b) const {
    if (!direct_ntt_.empty()) {
      return direct_ntt_[0].Multiply(a, b);
    }
    return crt_plan_.Multiply(a, b);
  }

  std::vector<uint64_t> Transform(std::span<const uint64_t> input,
                                  const std::vector<uint64_t>& chirp_neg,
                                  const std::vector<uint64_t>& kernel,
                                  bool scale_by_inv_n) const {
    if (input.size() != length_) {
      throw std::invalid_argument(
          "BluesteinDft: input length does not match plan length");
    }
    if (length_ == 1) {
      return {input[0] % prime_};
    }

    // Reversed weighted input: a_rev[N - 1 - j] = x[j] * omega^(-C(j, 2)).
    std::vector<uint64_t> a_rev(length_, 0ULL);
    for (std::size_t j = 0; j < length_; ++j) {
      const uint64_t xj = (input[j] >= prime_) ? (input[j] % prime_) : input[j];
      a_rev[length_ - 1 - j] = utils::ModMul(xj, chirp_neg[j], prime_);
    }

    const std::vector<uint64_t> conv = Convolve(
        std::span<const uint64_t>(a_rev), std::span<const uint64_t>(kernel));

    std::vector<uint64_t> out(length_, 0ULL);
    for (std::size_t k = 0; k < length_; ++k) {
      const std::size_t idx = length_ - 1 + k;
      const uint64_t val = (idx < conv.size()) ? conv[idx] : 0ULL;
      uint64_t res = utils::ModMul(val, chirp_neg[k], prime_);
      if (scale_by_inv_n) {
        res = utils::ModMul(res, inv_length_, prime_);
      }
      out[k] = res;
    }
    return out;
  }

  uint64_t prime_;
  std::size_t length_;
  uint64_t omega_;
  uint64_t inv_omega_;
  uint64_t inv_length_;
  std::vector<uint64_t> fwd_chirp_neg_;
  std::vector<uint64_t> inv_chirp_neg_;
  std::vector<uint64_t> fwd_kernel_;
  std::vector<uint64_t> inv_kernel_;
  std::vector<NttPlan> direct_ntt_;
  CrtNttPlan crt_plan_;
};

}  // namespace poly

using poly::BluesteinDft;
using poly::CrtNttPlan;
using poly::NttPlan;

}  // namespace xg

#endif  // XGALOIS_POLY_NTT_HPP_
