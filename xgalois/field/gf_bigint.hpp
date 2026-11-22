#ifndef XGALOIS_FIELD_GF_BIGINT_HPP_
#define XGALOIS_FIELD_GF_BIGINT_HPP_

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {

// ============================================================================
// GaloisFieldPrimeBig<Bits>: Fixed-width big prime field GF(p) backed by
// BigUInt<Bits> and MontgomeryCtx<Bits>.
// ============================================================================
//
// Elements are canonical BigUInt<Bits> values in [0, p - 1], making interop
// with FieldTower, serialization, and equality checks seamless. Internal
// multiplication, squaring, exponentiation, and Tonelli-Shanks square roots
// leverage MontgomeryCtx<Bits> for division-free CIOS Montgomery arithmetic.
template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
class GaloisFieldPrimeBig {
 public:
  static constexpr std::size_t kBits = Bits;
  using element_type = BigUInt<Bits>;
  using storage_type = BigUInt<Bits>;

  explicit GaloisFieldPrimeBig(const BigUInt<Bits>& prime,
                               bool check_prime = true)
      : p_(prime) {
    if (p_ < BigUInt<Bits>(2ULL)) {
      throw std::invalid_argument(
          "GaloisFieldPrimeBig: prime modulus must be >= 2");
    }
    if (check_prime && !IsProbablePrime(BigInt::FromBigUInt(p_))) {
      throw std::invalid_argument(
          "GaloisFieldPrimeBig: provided modulus is not prime");
    }
    if (p_.IsOdd()) {
      mont_ctx_.emplace(p_);
    }
  }

  static GaloisFieldPrimeBig FromHex(std::string_view hex,
                                     bool check_prime = true) {
    return GaloisFieldPrimeBig(BigUInt<Bits>::FromHex(hex), check_prime);
  }

  static GaloisFieldPrimeBig FromDecimal(std::string_view dec,
                                         bool check_prime = true) {
    return GaloisFieldPrimeBig(BigUInt<Bits>::FromDecimal(dec), check_prime);
  }

  const BigUInt<Bits>& Modulus() const noexcept { return p_; }

  BigInt BigCharacteristic() const { return BigInt::FromBigUInt(p_); }
  BigInt BigOrder() const { return BigInt::FromBigUInt(p_); }
  uint32_t Degree() const noexcept { return 1; }

  uint64_t Characteristic() const noexcept {
    return (p_.BitLength() <= 64) ? p_.Limb(0) : 0ULL;
  }

  uint64_t Order() const noexcept { return Characteristic(); }

  const std::optional<MontgomeryCtx<Bits>>& MontgomeryContext() const noexcept {
    return mont_ctx_;
  }

  element_type AdditiveIdentity() const noexcept {
    return BigUInt<Bits>::Zero();
  }

  element_type MultiplicativeIdentity() const noexcept {
    return BigUInt<Bits>::One();
  }

  element_type Reduce(const element_type& a) const {
    return (a >= p_) ? (a % p_) : a;
  }

  bool IsZero(const element_type& a) const { return Reduce(a).IsZero(); }
  bool IsOne(const element_type& a) const { return Reduce(a).IsOne(); }

  bool Equals(const element_type& a, const element_type& b) const {
    return Reduce(a) == Reduce(b);
  }

  element_type Add(const element_type& a, const element_type& b) const {
    const element_type ar = Reduce(a);
    const element_type br = Reduce(b);
    if (mont_ctx_.has_value()) {
      return mont_ctx_->Add(ar, br);
    }
    return ModAdd(ar, br, p_);
  }

  element_type Sub(const element_type& a, const element_type& b) const {
    const element_type ar = Reduce(a);
    const element_type br = Reduce(b);
    if (mont_ctx_.has_value()) {
      return mont_ctx_->Sub(ar, br);
    }
    return ModSub(ar, br, p_);
  }

  element_type Neg(const element_type& a) const {
    const element_type ar = Reduce(a);
    if (mont_ctx_.has_value()) {
      return mont_ctx_->Neg(ar);
    }
    return ar.IsZero() ? ar : (p_ - ar);
  }

  // Division-free multiplication using CIOS Montgomery reduction:
  //   r1 = REDC(a * b) = a * b * R^(-1) mod p
  //   r2 = REDC(r1 * R^2) = a * b mod p.
  element_type Mul(const element_type& a, const element_type& b) const {
    const element_type ar = Reduce(a);
    const element_type br = Reduce(b);
    if (mont_ctx_.has_value()) {
      const element_type ab_r_inv = mont_ctx_->Mul(ar, br);
      return mont_ctx_->Mul(ab_r_inv, mont_ctx_->R2Mod());
    }
    return BigUInt<Bits>((ar.Limb(0) & br.Limb(0)) & 1ULL);
  }

  element_type Sqr(const element_type& a) const {
    const element_type ar = Reduce(a);
    if (mont_ctx_.has_value()) {
      const element_type a2_r_inv = mont_ctx_->Sqr(ar);
      return mont_ctx_->Mul(a2_r_inv, mont_ctx_->R2Mod());
    }
    return ar;
  }

  // Explicit Montgomery-domain helpers for long chains of multiplications.
  element_type ToMont(const element_type& a) const {
    RequireOddModulus("ToMont");
    return mont_ctx_->ToMont(a);
  }

  element_type FromMont(const element_type& a_mont) const {
    RequireOddModulus("FromMont");
    return mont_ctx_->FromMont(a_mont);
  }

  element_type MulMont(const element_type& a_mont,
                       const element_type& b_mont) const {
    RequireOddModulus("MulMont");
    return mont_ctx_->Mul(a_mont, b_mont);
  }

  element_type SqrMont(const element_type& a_mont) const {
    RequireOddModulus("SqrMont");
    return mont_ctx_->Sqr(a_mont);
  }

  element_type Inv(const element_type& a) const {
    const element_type ar = Reduce(a);
    if (ar.IsZero()) {
      throw std::domain_error(
          "GaloisFieldPrimeBig::Inv: inverse of zero is undefined");
    }
    if (ar.IsOne()) {
      return BigUInt<Bits>::One();
    }
    return ModInv(ar, p_);
  }

  element_type Div(const element_type& a, const element_type& b) const {
    return Mul(a, Inv(b));
  }

  element_type Pow(const element_type& a, uint64_t exp) const {
    if (exp == 0ULL) {
      return BigUInt<Bits>::One();
    }
    const element_type ar = Reduce(a);
    if (ar.IsZero()) {
      return BigUInt<Bits>::Zero();
    }
    if (ar.IsOne()) {
      return BigUInt<Bits>::One();
    }
    if (mont_ctx_.has_value()) {
      return mont_ctx_->FromMont(
          mont_ctx_->Pow(mont_ctx_->ToMont(ar), BigUInt<Bits>(exp)));
    }
    return ar;
  }

  element_type PowBig(const element_type& a, const BigInt& exp) const {
    if (exp.IsZero()) {
      return BigUInt<Bits>::One();
    }
    if (exp.IsNegative()) {
      return PowBig(Inv(a), -exp);
    }
    const element_type ar = Reduce(a);
    if (ar.IsZero()) {
      return BigUInt<Bits>::Zero();
    }
    if (ar.IsOne()) {
      return BigUInt<Bits>::One();
    }
    if (!mont_ctx_.has_value()) {
      return ar;
    }
    const BigInt p_minus_1 = BigInt::FromBigUInt(p_ - BigUInt<Bits>::One());
    const BigInt e_red = BigInt::EuclideanDivRem(exp, p_minus_1).second;
    if (e_red.IsZero()) {
      return BigUInt<Bits>::One();
    }
    const BigUInt<Bits> e_uint = e_red.template ToBigUInt<Bits>();
    return mont_ctx_->FromMont(mont_ctx_->Pow(mont_ctx_->ToMont(ar), e_uint));
  }

  // Euler's criterion: returns true iff a is a quadratic residue mod p.
  bool HasSqrt(const element_type& a) const {
    const element_type ar = Reduce(a);
    if (ar.IsZero() || !mont_ctx_.has_value()) {
      return true;
    }
    const element_type exp_half = (p_ - BigUInt<Bits>::One()) >> 1;
    const element_type leg = mont_ctx_->Pow(mont_ctx_->ToMont(ar), exp_half);
    return leg == mont_ctx_->RMod();
  }

  // Tonelli-Shanks square root in GF(p). Returns std::nullopt if a is a
  // quadratic non-residue, otherwise returns the canonical root min(r, p - r).
  std::optional<element_type> Sqrt(const element_type& a) const {
    const element_type ar = Reduce(a);
    if (ar.IsZero()) {
      return BigUInt<Bits>::Zero();
    }
    if (!mont_ctx_.has_value()) {
      // p == 2: sqrt(1) == 1.
      return ar;
    }

    const element_type a_m = mont_ctx_->ToMont(ar);
    const element_type p_minus_1 = p_ - BigUInt<Bits>::One();
    const element_type exp_half = p_minus_1 >> 1;

    // Euler criterion check.
    if (mont_ctx_->Pow(a_m, exp_half) != mont_ctx_->RMod()) {
      return std::nullopt;
    }

    element_type r_m;
    if ((p_.Limb(0) & 3ULL) == 3ULL) {
      // Fast path for p == 3 (mod 4): r = a^((p + 1) / 4).
      const element_type exp_q = (p_ >> 2) + BigUInt<Bits>::One();
      r_m = mont_ctx_->Pow(a_m, exp_q);
    } else {
      // General Tonelli-Shanks: p - 1 = q * 2^s with q odd.
      const std::size_t s = detail::CountTrailingZeroBitsUInt(p_minus_1);
      const element_type q = p_minus_1 >> s;

      // Find smallest quadratic non-residue z >= 2.
      element_type z_m;
      for (uint64_t z = 2ULL;; ++z) {
        const element_type cand = mont_ctx_->ToMont(BigUInt<Bits>(z));
        if (mont_ctx_->Pow(cand, exp_half) != mont_ctx_->RMod()) {
          z_m = cand;
          break;
        }
      }

      std::size_t m = s;
      element_type c = mont_ctx_->Pow(z_m, q);
      element_type t = mont_ctx_->Pow(a_m, q);
      r_m = mont_ctx_->Pow(a_m, (q >> 1) + BigUInt<Bits>::One());

      while (t != mont_ctx_->RMod()) {
        std::size_t i = 0;
        element_type tmp = t;
        while (tmp != mont_ctx_->RMod()) {
          tmp = mont_ctx_->Sqr(tmp);
          ++i;
          if (i >= m) {
            return std::nullopt;
          }
        }
        element_type b = c;
        for (std::size_t j = 0; j < m - i - 1; ++j) {
          b = mont_ctx_->Sqr(b);
        }
        m = i;
        c = mont_ctx_->Sqr(b);
        t = mont_ctx_->Mul(t, c);
        r_m = mont_ctx_->Mul(r_m, b);
      }
    }

    const element_type r = mont_ctx_->FromMont(r_m);
    const element_type neg_r = Neg(r);
    return (neg_r < r) ? neg_r : r;
  }

  std::string ToString(const element_type& a) const {
    return Reduce(a).ToDecimal();
  }

 private:
  void RequireOddModulus(const char* op) const {
    if (!mont_ctx_.has_value()) {
      throw std::domain_error(std::string("GaloisFieldPrimeBig::") + op +
                              " requires an odd prime modulus");
    }
  }

  BigUInt<Bits> p_;
  std::optional<MontgomeryCtx<Bits>> mont_ctx_;
};

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
struct AlgebraTraits<GaloisFieldPrimeBig<Bits>> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = true;
};

// ============================================================================
// GaloisFieldPrimeDyn: Arbitrary-precision prime field GF(p) backed by BigInt.
// ============================================================================

class GaloisFieldPrimeDyn {
 public:
  using element_type = BigInt;
  using storage_type = BigInt;

  explicit GaloisFieldPrimeDyn(BigInt prime, bool check_prime = true)
      : p_(std::move(prime)) {
    if (p_ < BigInt(2)) {
      throw std::invalid_argument(
          "GaloisFieldPrimeDyn: prime modulus must be >= 2");
    }
    if (check_prime && !IsProbablePrime(p_)) {
      throw std::invalid_argument(
          "GaloisFieldPrimeDyn: provided modulus is not prime");
    }
  }

  static GaloisFieldPrimeDyn FromHex(std::string_view hex,
                                     bool check_prime = true) {
    return GaloisFieldPrimeDyn(BigInt::FromHex(hex), check_prime);
  }

  static GaloisFieldPrimeDyn FromDecimal(std::string_view dec,
                                         bool check_prime = true) {
    return GaloisFieldPrimeDyn(BigInt::FromDecimal(dec), check_prime);
  }

  const BigInt& Modulus() const noexcept { return p_; }
  BigInt BigCharacteristic() const { return p_; }
  BigInt BigOrder() const { return p_; }
  uint32_t Degree() const noexcept { return 1; }

  uint64_t Characteristic() const noexcept {
    return (p_.BitLength() <= 64) ? p_.ToUint64() : 0ULL;
  }

  uint64_t Order() const noexcept { return Characteristic(); }

  BigInt AdditiveIdentity() const { return BigInt(0); }
  BigInt MultiplicativeIdentity() const { return BigInt(1); }

  BigInt Reduce(const BigInt& a) const {
    if (!a.IsNegative() && a < p_) {
      return a;
    }
    return BigInt::EuclideanDivRem(a, p_).second;
  }

  bool IsZero(const BigInt& a) const { return Reduce(a).IsZero(); }
  bool IsOne(const BigInt& a) const { return Reduce(a).IsOne(); }
  bool Equals(const BigInt& a, const BigInt& b) const {
    return Reduce(a) == Reduce(b);
  }

  BigInt Add(const BigInt& a, const BigInt& b) const {
    BigInt sum = Reduce(a) + Reduce(b);
    if (sum >= p_) {
      sum -= p_;
    }
    return sum;
  }

  BigInt Sub(const BigInt& a, const BigInt& b) const {
    const BigInt ar = Reduce(a);
    const BigInt br = Reduce(b);
    return (ar >= br) ? (ar - br) : (ar + p_ - br);
  }

  BigInt Neg(const BigInt& a) const {
    const BigInt ar = Reduce(a);
    return ar.IsZero() ? ar : (p_ - ar);
  }

  BigInt Mul(const BigInt& a, const BigInt& b) const {
    return (Reduce(a) * Reduce(b)) % p_;
  }

  BigInt Inv(const BigInt& a) const {
    const BigInt ar = Reduce(a);
    if (ar.IsZero()) {
      throw std::domain_error(
          "GaloisFieldPrimeDyn::Inv: inverse of zero is undefined");
    }
    if (ar.IsOne()) {
      return BigInt(1);
    }
    return ModInverse(ar, p_);
  }

  BigInt Div(const BigInt& a, const BigInt& b) const { return Mul(a, Inv(b)); }

  BigInt Pow(const BigInt& a, uint64_t exp) const {
    return PowBig(a, BigInt::FromUint64(exp));
  }

  BigInt PowBig(const BigInt& a, const BigInt& exp) const {
    if (exp.IsZero()) {
      return BigInt(1);
    }
    if (exp.IsNegative()) {
      return PowBig(Inv(a), -exp);
    }
    const BigInt ar = Reduce(a);
    if (ar.IsZero()) {
      return BigInt(0);
    }
    if (ar.IsOne()) {
      return BigInt(1);
    }
    const BigInt p_minus_1 = p_ - BigInt(1);
    const BigInt e_red = BigInt::EuclideanDivRem(exp, p_minus_1).second;
    return PowMod(ar, e_red, p_);
  }

  bool HasSqrt(const BigInt& a) const {
    const BigInt ar = Reduce(a);
    if (ar.IsZero() || p_ == BigInt(2)) {
      return true;
    }
    return JacobiSymbol(ar, p_) == 1;
  }

  std::optional<BigInt> Sqrt(const BigInt& a) const {
    const BigInt ar = Reduce(a);
    if (ar.IsZero()) {
      return BigInt(0);
    }
    if (p_ == BigInt(2)) {
      return ar;
    }
    if (JacobiSymbol(ar, p_) != 1) {
      return std::nullopt;
    }

    BigInt r;
    if ((p_.Limb(0) & 3ULL) == 3ULL) {
      r = PowMod(ar, (p_ + BigInt(1)) >> 2, p_);
    } else {
      const BigInt p_minus_1 = p_ - BigInt(1);
      const std::size_t s = p_minus_1.CountTrailingZeroBits();
      const BigInt q = p_minus_1 >> s;

      BigInt z(2);
      while (JacobiSymbol(z, p_) != -1) {
        z += BigInt(1);
      }

      std::size_t m = s;
      BigInt c = PowMod(z, q, p_);
      BigInt t = PowMod(ar, q, p_);
      r = PowMod(ar, (q + BigInt(1)) >> 1, p_);

      while (!t.IsOne()) {
        std::size_t i = 0;
        BigInt tmp = t;
        while (!tmp.IsOne()) {
          tmp = (tmp * tmp) % p_;
          ++i;
          if (i >= m) {
            return std::nullopt;
          }
        }
        BigInt b = c;
        for (std::size_t j = 0; j < m - i - 1; ++j) {
          b = (b * b) % p_;
        }
        m = i;
        c = (b * b) % p_;
        t = (t * c) % p_;
        r = (r * b) % p_;
      }
    }

    const BigInt neg_r = Neg(r);
    return (neg_r < r) ? neg_r : r;
  }

  std::string ToString(const BigInt& a) const { return Reduce(a).ToDecimal(); }

 private:
  BigInt p_;
};

template <>
struct AlgebraTraits<GaloisFieldPrimeDyn> {
  static constexpr bool kIsAssociativeOp = false;
  static constexpr bool kIsCommutativeOp = false;
  static constexpr bool kIsCommutativeMul = true;
  static constexpr bool kHasNoZeroDivisors = true;
  static constexpr bool kIsPid = true;
  static constexpr bool kIsEuclideanDomain = true;
  static constexpr bool kIsField = true;
};

// ============================================================================
// Named Cryptographic Prime Field Factories (namespace xg::primes)
// ============================================================================

namespace primes {

// Curve25519 / Ed25519 base field: p = 2^255 - 19.
inline GaloisFieldPrimeBig<256> Curve25519() {
  return GaloisFieldPrimeBig<256>::FromHex(
      "0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed",
      /*check_prime=*/false);
}

// secp256k1 base field: p = 2^256 - 2^32 - 977.
inline GaloisFieldPrimeBig<256> Secp256k1() {
  return GaloisFieldPrimeBig<256>::FromHex(
      "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
      /*check_prime=*/false);
}

// NIST P-256 (secp256r1) base field: p = 2^256 - 2^224 + 2^192 + 2^96 - 1.
inline GaloisFieldPrimeBig<256> NistP256() {
  return GaloisFieldPrimeBig<256>::FromHex(
      "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
      /*check_prime=*/false);
}

// NIST P-384 (secp384r1) base field: p = 2^384 - 2^128 - 2^96 + 2^32 - 1.
inline GaloisFieldPrimeBig<384> NistP384() {
  return GaloisFieldPrimeBig<384>::FromHex(
      "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe"
      "ffffffff0000000000000000ffffffff",
      /*check_prime=*/false);
}

// BN254 (alt_bn128) base field prime.
inline GaloisFieldPrimeBig<256> Bn254() {
  return GaloisFieldPrimeBig<256>::FromHex(
      "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47",
      /*check_prime=*/false);
}

// BLS12-381 base field prime (381-bit prime).
inline GaloisFieldPrimeBig<384> Bls12381() {
  return GaloisFieldPrimeBig<384>::FromHex(
      "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
      "1eabfffeb153ffffb9feffffffffaaab",
      /*check_prime=*/false);
}

// Goldilocks prime (64-bit STARK/FRI-friendly prime): p = 2^64 - 2^32 + 1.
inline GaloisFieldPrimeBig<64> Goldilocks() {
  return GaloisFieldPrimeBig<64>::FromHex("0xffffffff00000001",
                                          /*check_prime=*/false);
}

}  // namespace primes

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_BIGINT_HPP_
