#ifndef XGALOIS_FIELD_GF_NORMAL_BASIS_HPP_
#define XGALOIS_FIELD_GF_NORMAL_BASIS_HPP_

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {

enum class OnbType : uint8_t {
  kNone = 0,
  kType1 = 1,
  kType2 = 2,
};

namespace detail {

// Computes the exact multiplicative order of `base` modulo prime `p`.
inline uint64_t MultiplicativeOrderModPrime(uint64_t base, uint64_t p) {
  if (p <= 2) {
    return 1;
  }
  uint64_t ord = p - 1;
  const auto factors = utils::DistinctPrimeFactors(ord);
  for (uint64_t f : factors) {
    while (ord % f == 0 && utils::ModPow(base, ord / f, p) == 1) {
      ord /= f;
    }
  }
  return ord;
}

// Inverts an m x m binary matrix over GF(2) (1 <= m <= 64) given as m column
// bitmasks `cols[j]` (bit i of cols[j] is entry (i, j)). Returns the m column
// bitmasks of the inverse matrix.
inline std::array<uint64_t, 64> InvertBinaryMatrixCols(
    const std::array<uint64_t, 64>& cols, uint32_t m) {
  std::array<uint64_t, 64> rows{};
  std::array<uint64_t, 64> inv_rows{};
  for (uint32_t r = 0; r < m; ++r) {
    inv_rows[r] = uint64_t{1} << r;
    uint64_t row_bits = 0;
    for (uint32_t c = 0; c < m; ++c) {
      if (((cols[c] >> r) & 1ULL) != 0) {
        row_bits |= (uint64_t{1} << c);
      }
    }
    rows[r] = row_bits;
  }

  for (uint32_t col = 0; col < m; ++col) {
    uint32_t pivot = m;
    for (uint32_t r = col; r < m; ++r) {
      if (((rows[r] >> col) & 1ULL) != 0) {
        pivot = r;
        break;
      }
    }
    if (pivot == m) {
      throw std::runtime_error("Binary change-of-basis matrix is singular");
    }
    if (pivot != col) {
      std::swap(rows[pivot], rows[col]);
      std::swap(inv_rows[pivot], inv_rows[col]);
    }
    for (uint32_t r = 0; r < m; ++r) {
      if (r != col && (((rows[r] >> col) & 1ULL) != 0)) {
        rows[r] ^= rows[col];
        inv_rows[r] ^= inv_rows[col];
      }
    }
  }

  std::array<uint64_t, 64> inv_cols{};
  for (uint32_t c = 0; c < m; ++c) {
    uint64_t col_bits = 0;
    for (uint32_t r = 0; r < m; ++r) {
      if (((inv_rows[r] >> c) & 1ULL) != 0) {
        col_bits |= (uint64_t{1} << r);
      }
    }
    inv_cols[c] = col_bits;
  }
  return inv_cols;
}

}  // namespace detail

// Detects whether GF(2^m) admits a Type I or Type II Optimal Normal Basis.
//   - Type I:  p = m + 1 is prime and 2 is a primitive root mod p (ord_p(2)=m).
//   - Type II: p = 2m + 1 is prime and either ord_p(2) = 2m, or
//              (p == 3 (mod 4) and ord_p(2) = m).
inline OnbType DetectOnbType(uint32_t m) {
  if (m < 2) {
    return OnbType::kNone;
  }

  const uint64_t p1 = static_cast<uint64_t>(m) + 1ULL;
  if (utils::IsPrime(p1) && detail::MultiplicativeOrderModPrime(2, p1) == m) {
    return OnbType::kType1;
  }

  const uint64_t p2 = 2ULL * static_cast<uint64_t>(m) + 1ULL;
  if (utils::IsPrime(p2)) {
    const uint64_t ord = detail::MultiplicativeOrderModPrime(2, p2);
    if (ord == 2ULL * m || ((p2 & 3ULL) == 3ULL && ord == m)) {
      return OnbType::kType2;
    }
  }

  return OnbType::kNone;
}

// ============================================================================
// GaloisFieldNormalBinary: GF(2^m) in Optimal Normal Basis (Type I / Type II)
// ============================================================================
//
// Elements are represented in the normal basis N = {beta, beta^2, ...,
// beta^{2^{m-1}}}, where bit i stores the coefficient of beta^{2^i}:
//   - AdditiveIdentity()       = 0
//   - MultiplicativeIdentity() = (1 << m) - 1  (all-ones vector)
//   - Frobenius(a)             = RotL(a, 1) mod m  (O(1))
//   - Sqrt(a)                  = RotR(a, 1) mod m  (O(1))
//   - Trace(a)                 = popcount(a) & 1   (O(1))
//   - Mul(a, b)                = O(m) sparse Massey-Omura / ONB bit-parallel
//   - Inv(a)                   = O(log m) Itoh-Tsujii inversion
class GaloisFieldNormalBinary {
 public:
  using element_type = uint64_t;
  using storage_type = uint64_t;

  explicit GaloisFieldNormalBinary(uint32_t m) : m_(m) {
    if (m_ < 2 || m_ > 63) {
      throw std::invalid_argument(
          "GaloisFieldNormalBinary: degree m must satisfy 2 <= m <= 63");
    }
    type_ = DetectOnbType(m_);
    if (type_ == OnbType::kNone) {
      throw std::invalid_argument(
          "GaloisFieldNormalBinary: GF(2^" + std::to_string(m_) +
          ") does not have a Type I or Type II Optimal Normal Basis");
    }
    mask_ = (uint64_t{1} << m_) - 1ULL;
    order_ = uint64_t{1} << m_;
    BuildLambdaTable();
  }

  uint32_t Degree() const noexcept { return m_; }
  uint64_t Characteristic() const noexcept { return 2; }
  uint64_t Order() const noexcept { return order_; }
  OnbType Type() const noexcept { return type_; }
  uint32_t Complexity() const noexcept { return 2U * m_ - 1U; }

  bool LambdaEntry(uint32_t i, uint32_t j) const {
    if (i >= m_ || j >= m_) {
      throw std::out_of_range("LambdaEntry index out of range");
    }
    return ((lambda_rows_[i] >> j) & 1ULL) != 0;
  }

  uint64_t AdditiveIdentity() const noexcept { return 0ULL; }
  uint64_t MultiplicativeIdentity() const noexcept { return mask_; }

  bool IsZero(uint64_t a) const noexcept { return (a & mask_) == 0ULL; }
  bool IsOne(uint64_t a) const noexcept { return (a & mask_) == mask_; }
  bool Equals(uint64_t a, uint64_t b) const noexcept {
    return (a & mask_) == (b & mask_);
  }

  uint64_t Add(uint64_t a, uint64_t b) const noexcept {
    return (a ^ b) & mask_;
  }

  uint64_t Sub(uint64_t a, uint64_t b) const noexcept {
    return (a ^ b) & mask_;
  }

  uint64_t Neg(uint64_t a) const noexcept { return a & mask_; }

  // O(1) Frobenius endomorphism a -> a^2 (cyclic left rotation by 1 bit mod m).
  uint64_t Frobenius(uint64_t a) const noexcept { return RotL(a & mask_, 1); }

  // O(1) square root a -> a^(2^(m-1)) (cyclic right rotation by 1 bit mod m).
  uint64_t Sqrt(uint64_t a) const noexcept { return RotR(a & mask_, 1); }

  // O(1) absolute trace Tr_{GF(2^m)/GF(2)}(a) = popcount(a) mod 2.
  uint8_t Trace(uint64_t a) const noexcept {
    return static_cast<uint8_t>(std::popcount(a & mask_) & 1);
  }

  // Sparse Massey-Omura ONB multiplication in O(m) word operations:
  //   c_k = sum_{u, v} lambda(u, v) * a_{(u+k) mod m} * b_{(v+k) mod m}
  // Computed bit-parallel across all m coordinates simultaneously:
  //   c = XOR_{u=0}^{m-1} (RotR(a, u) & XOR_{v: lambda(u,v)=1} RotR(b, v)).
  uint64_t Mul(uint64_t a, uint64_t b) const noexcept {
    a &= mask_;
    b &= mask_;
    if (a == 0ULL || b == 0ULL) {
      return 0ULL;
    }
    uint64_t acc = 0ULL;
    for (uint32_t u = 0; u < m_; ++u) {
      const uint8_t cnt = row_count_[u];
      uint64_t b_term = RotR(b, row_cols_[u][0]);
      if (cnt == 2) {
        b_term ^= RotR(b, row_cols_[u][1]);
      }
      acc ^= (RotR(a, u) & b_term);
    }
    return acc & mask_;
  }

  uint64_t Pow(uint64_t a, uint64_t exp) const noexcept {
    a &= mask_;
    if (exp == 0ULL) {
      return mask_;
    }
    if (a == 0ULL) {
      return 0ULL;
    }
    exp %= (order_ - 1ULL);
    if (exp == 0ULL) {
      return mask_;
    }
    uint64_t result = mask_;
    uint64_t base = a;
    while (exp > 0ULL) {
      if ((exp & 1ULL) != 0ULL) {
        result = Mul(result, base);
      }
      exp >>= 1;
      if (exp > 0ULL) {
        base = Frobenius(base);
      }
    }
    return result;
  }

  // O(log m) Itoh-Tsujii inversion: a^(-1) = (a^(2^(m-1) - 1))^2, using O(1)
  // cyclic rotations for every 2^k-th power.
  uint64_t Inv(uint64_t a) const {
    a &= mask_;
    if (a == 0ULL) {
      throw std::domain_error(
          "GaloisFieldNormalBinary::Inv: inverse of zero is undefined");
    }
    if (a == mask_) {
      return mask_;
    }

    const uint32_t r = m_ - 1U;
    const int msb = std::bit_width(r) - 1;
    uint64_t res = a;  // res = a^(2^k - 1) with k = 1.
    uint32_t k = 1U;
    for (int bit = msb - 1; bit >= 0; --bit) {
      res = Mul(RotL(res, k), res);  // a^(2^(2k) - 1)
      k <<= 1;
      if (((r >> bit) & 1U) != 0U) {
        res = Mul(RotL(res, 1U), a);  // a^(2^(k+1) - 1)
        ++k;
      }
    }
    return RotL(res, 1U);  // (a^(2^(m-1) - 1))^2 = a^(2^m - 2) = a^(-1).
  }

  uint64_t Div(uint64_t a, uint64_t b) const { return Mul(a, Inv(b)); }

  std::string ToString(uint64_t a) const {
    a &= mask_;
    std::ostringstream oss;
    oss << "0b";
    for (uint32_t i = m_; i-- > 0;) {
      oss << (((a >> i) & 1ULL) != 0ULL ? '1' : '0');
    }
    return oss.str();
  }

 private:
  uint64_t RotL(uint64_t x, uint32_t s) const noexcept {
    s %= m_;
    if (s == 0U) {
      return x;
    }
    return ((x << s) | (x >> (m_ - s))) & mask_;
  }

  uint64_t RotR(uint64_t x, uint32_t s) const noexcept {
    s %= m_;
    if (s == 0U) {
      return x;
    }
    return ((x >> s) | (x << (m_ - s))) & mask_;
  }

  void BuildLambdaTable() {
    lambda_rows_.fill(0ULL);
    row_count_.fill(0U);

    const uint64_t p = (type_ == OnbType::kType1)
                           ? (static_cast<uint64_t>(m_) + 1ULL)
                           : (2ULL * static_cast<uint64_t>(m_) + 1ULL);

    // Precompute powers 2^u mod p for u = 0..m-1.
    std::array<uint64_t, 64> pow2{};
    pow2[0] = 1ULL;
    for (uint32_t u = 1; u < m_; ++u) {
      pow2[u] = (pow2[u - 1] * 2ULL) % p;
    }

    for (uint32_t u = 0; u < m_; ++u) {
      for (uint32_t v = 0; v < m_; ++v) {
        bool is_one = false;
        if (type_ == OnbType::kType1) {
          // Type I: 2^u + 2^v == 0 or 1 (mod p).
          const uint64_t sum = (pow2[u] + pow2[v]) % p;
          is_one = (sum == 0ULL || sum == 1ULL);
        } else {
          // Type II: 2^u +/- 2^v == +/- 1 (mod p).
          const uint64_t sum = (pow2[u] + pow2[v]) % p;
          const uint64_t diff = (pow2[u] + p - pow2[v]) % p;
          is_one = (sum == 1ULL || sum == p - 1ULL || diff == 1ULL ||
                    diff == p - 1ULL);
        }
        if (is_one) {
          lambda_rows_[u] |= (uint64_t{1} << v);
          row_cols_[u][row_count_[u]++] = static_cast<uint8_t>(v);
        }
      }
    }
  }

  uint32_t m_ = 0;
  OnbType type_ = OnbType::kNone;
  uint64_t mask_ = 0;
  uint64_t order_ = 0;
  std::array<uint64_t, 64> lambda_rows_{};
  std::array<std::array<uint8_t, 2>, 64> row_cols_{};
  std::array<uint8_t, 64> row_count_{};
};

template <>
struct AlgebraTraits<GaloisFieldNormalBinary> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = true;
};

// ============================================================================
// NormalToPolyConverter: Exact basis-change isomorphism between
// GaloisFieldNormalBinary(m) and a polynomial-basis binary extension field
// GaloisFieldBinaryExtension<Word>(m).
// ============================================================================

template <typename BinaryExtField = GaloisFieldBinaryExtension<uint64_t>>
class NormalToPolyConverter {
 public:
  using poly_element_type = typename BinaryExtField::element_type;

  NormalToPolyConverter(const GaloisFieldNormalBinary& normal_field,
                        const BinaryExtField& poly_field)
      : normal_field_(&normal_field),
        poly_field_(&poly_field),
        m_(normal_field.Degree()) {
    if (poly_field.Characteristic() != 2 || poly_field.Degree() != m_) {
      throw std::invalid_argument(
          "NormalToPolyConverter: fields must both be GF(2^m) with the same "
          "degree m");
    }

    beta_poly_ = FindGaussPeriodNormalElement();

    // Column i of normal_to_poly_cols_ is the polynomial-basis encoding of
    // beta^(2^i) for i = 0..m-1.
    poly_element_type cur = beta_poly_;
    for (uint32_t i = 0; i < m_; ++i) {
      normal_to_poly_cols_[i] =
          static_cast<uint64_t>(poly_field_->GetElementValue(cur));
      if (i + 1 < m_) {
        cur = poly_field_->Mul(cur, cur);
      }
    }

    poly_to_normal_cols_ =
        detail::InvertBinaryMatrixCols(normal_to_poly_cols_, m_);
  }

  poly_element_type ToPoly(uint64_t normal_elem) const {
    uint64_t res = 0ULL;
    for (uint32_t i = 0; i < m_; ++i) {
      if (((normal_elem >> i) & 1ULL) != 0ULL) {
        res ^= normal_to_poly_cols_[i];
      }
    }
    return poly_field_->SetElementValue(static_cast<poly_element_type>(res));
  }

  uint64_t ToNormal(const poly_element_type& poly_elem) const {
    const uint64_t bits =
        static_cast<uint64_t>(poly_field_->GetElementValue(poly_elem));
    uint64_t res = 0ULL;
    for (uint32_t i = 0; i < m_; ++i) {
      if (((bits >> i) & 1ULL) != 0ULL) {
        res ^= poly_to_normal_cols_[i];
      }
    }
    return res;
  }

  poly_element_type NormalGeneratorInPoly() const { return beta_poly_; }

 private:
  uint8_t ComputePolyTrace(const poly_element_type& a) const {
    poly_element_type acc = a;
    poly_element_type cur = a;
    for (uint32_t i = 1; i < m_; ++i) {
      cur = poly_field_->Mul(cur, cur);
      acc = poly_field_->Add(acc, cur);
    }
    return poly_field_->IsZero(acc) ? uint8_t{0} : uint8_t{1};
  }

  // Finds a Gauss-period normal element beta in poly_field whose conjugates
  // {beta^(2^i)} obey the exact Type I or Type II ONB multiplication table.
  poly_element_type FindGaussPeriodNormalElement() const {
    const uint64_t group_order = (uint64_t{1} << m_) - 1ULL;
    const uint64_t mask = group_order;

    auto next_candidate = [mask](uint64_t step) -> uint64_t {
      if (step <= 1024 && step <= mask) {
        return step;
      }
      uint64_t z = step * 0x9e3779b97f4a7c15ULL;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      z ^= (z >> 31);
      const uint64_t v = z & mask;
      return (v < 2ULL) ? 2ULL : v;
    };

    if (normal_field_->Type() == OnbType::kType1) {
      const uint64_t p = static_cast<uint64_t>(m_) + 1ULL;
      const uint64_t cofactor = group_order / p;
      for (uint64_t step = 2; step < 100000; ++step) {
        const poly_element_type g = poly_field_->SetElementValue(
            static_cast<poly_element_type>(next_candidate(step)));
        if (poly_field_->IsZero(g)) {
          continue;
        }
        const poly_element_type zeta = poly_field_->Pow(g, cofactor);
        if (!poly_field_->IsOne(zeta) && !poly_field_->IsZero(zeta)) {
          return zeta;
        }
      }
      throw std::runtime_error(
          "NormalToPolyConverter: failed to find Type I Gauss period");
    }

    // Type II: p = 2m + 1 divides either 2^m - 1 or 2^m + 1.
    const uint64_t p = 2ULL * static_cast<uint64_t>(m_) + 1ULL;
    if (group_order % p == 0ULL) {
      const uint64_t cofactor = group_order / p;
      for (uint64_t step = 2; step < 100000; ++step) {
        const poly_element_type g = poly_field_->SetElementValue(
            static_cast<poly_element_type>(next_candidate(step)));
        if (poly_field_->IsZero(g)) {
          continue;
        }
        const poly_element_type zeta = poly_field_->Pow(g, cofactor);
        if (!poly_field_->IsOne(zeta) && !poly_field_->IsZero(zeta)) {
          return poly_field_->Add(zeta, poly_field_->Inv(zeta));
        }
      }
      throw std::runtime_error(
          "NormalToPolyConverter: failed to find Type II Gauss period in "
          "GF(2^m)");
    }

    // Otherwise p | (2^m + 1). Construct quadratic extension K2 =
    // GF(2^m)[y]/(y^2 + y + c) where Tr(c) = 1.
    poly_element_type c = poly_field_->MultiplicativeIdentity();
    for (uint64_t cand = 1; cand <= mask; ++cand) {
      const poly_element_type elem =
          poly_field_->SetElementValue(static_cast<poly_element_type>(cand));
      if (ComputePolyTrace(elem) == 1) {
        c = elem;
        break;
      }
    }

    struct QuadElem {
      poly_element_type u;
      poly_element_type v;
    };

    const auto& pf = *poly_field_;
    auto q_mul = [&](const QuadElem& a, const QuadElem& b) -> QuadElem {
      const poly_element_type u1u2 = pf.Mul(a.u, b.u);
      const poly_element_type v1v2 = pf.Mul(a.v, b.v);
      const poly_element_type u1v2 = pf.Mul(a.u, b.v);
      const poly_element_type v1u2 = pf.Mul(a.v, b.u);
      return {pf.Add(u1u2, pf.Mul(c, v1v2)), pf.Add(pf.Add(u1v2, v1u2), v1v2)};
    };

    auto q_inv = [&](const QuadElem& a) -> QuadElem {
      // Conjugate is (u + v, v); norm = u*(u + v) + c*v^2.
      const poly_element_type u_plus_v = pf.Add(a.u, a.v);
      const poly_element_type norm =
          pf.Add(pf.Mul(a.u, u_plus_v), pf.Mul(c, pf.Mul(a.v, a.v)));
      const poly_element_type inv_norm = pf.Inv(norm);
      return {pf.Mul(u_plus_v, inv_norm), pf.Mul(a.v, inv_norm)};
    };

    auto q_pow = [&](QuadElem base, uint64_t exp) -> QuadElem {
      QuadElem res{pf.MultiplicativeIdentity(), pf.AdditiveIdentity()};
      while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
          res = q_mul(res, base);
        }
        exp >>= 1;
        if (exp > 0ULL) {
          base = q_mul(base, base);
        }
      }
      return res;
    };

    const uint64_t cofactor = ((uint64_t{1} << m_) + 1ULL) / p;
    const poly_element_type one = pf.MultiplicativeIdentity();

    for (uint64_t step = 0; step < 100000; ++step) {
      const uint64_t u_raw = (step == 0) ? 0ULL : next_candidate(step);
      const QuadElem w{
          pf.SetElementValue(static_cast<poly_element_type>(u_raw)), one};
      // Project w into the norm-1 subgroup of order 2^m + 1: z = Conj(w) / w.
      const QuadElem w_conj{pf.Add(w.u, w.v), w.v};
      const QuadElem z = q_mul(w_conj, q_inv(w));
      const QuadElem zeta = q_pow(z, cofactor);
      if (!pf.IsZero(zeta.v)) {
        // Since Norm(zeta) == 1, zeta^(-1) == Conj(zeta) == (u + v, v),
        // hence beta = zeta + zeta^(-1) = (v, 0) in GF(2^m).
        return zeta.v;
      }
    }

    throw std::runtime_error(
        "NormalToPolyConverter: failed to find Type II Gauss period in "
        "GF(2^(2m))");
  }

  const GaloisFieldNormalBinary* normal_field_;
  const BinaryExtField* poly_field_;
  uint32_t m_ = 0;
  poly_element_type beta_poly_{};
  std::array<uint64_t, 64> normal_to_poly_cols_{};
  std::array<uint64_t, 64> poly_to_normal_cols_{};
};

template <typename BinaryExtField>
NormalToPolyConverter(const GaloisFieldNormalBinary&, const BinaryExtField&)
    -> NormalToPolyConverter<BinaryExtField>;

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_NORMAL_BASIS_HPP_
