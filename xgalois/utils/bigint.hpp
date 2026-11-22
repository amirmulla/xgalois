#ifndef XGALOIS_UTILS_BIGINT_HPP_
#define XGALOIS_UTILS_BIGINT_HPP_

#include <algorithm>
#include <array>
#include <bit>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace xg {

namespace detail {

__extension__ using u128 = unsigned __int128;
__extension__ using i128 = __int128;

inline constexpr std::size_t kKaratsubaThreshold = 32;
inline constexpr uint64_t kDecBase = 10000000000000000000ULL;  // 10^19
inline constexpr std::size_t kDecBaseDigits = 19;

inline void TrimLimbs(std::vector<uint64_t>& limbs) noexcept {
  while (!limbs.empty() && limbs.back() == 0) {
    limbs.pop_back();
  }
}

inline std::span<const uint64_t> TrimSpan(
    std::span<const uint64_t> s) noexcept {
  while (!s.empty() && s.back() == 0) {
    s = s.subspan(0, s.size() - 1);
  }
  return s;
}

inline std::strong_ordering CmpLimbs(std::span<const uint64_t> a,
                                     std::span<const uint64_t> b) noexcept {
  a = TrimSpan(a);
  b = TrimSpan(b);
  if (a.size() != b.size()) {
    return a.size() <=> b.size();
  }
  for (std::size_t i = a.size(); i-- > 0;) {
    if (a[i] != b[i]) {
      return a[i] <=> b[i];
    }
  }
  return std::strong_ordering::equal;
}

inline std::vector<uint64_t> AddLimbs(std::span<const uint64_t> a,
                                      std::span<const uint64_t> b) {
  a = TrimSpan(a);
  b = TrimSpan(b);
  const std::size_t max_len = std::max(a.size(), b.size());
  std::vector<uint64_t> res(max_len, 0);
  uint64_t carry = 0;
  for (std::size_t i = 0; i < max_len; ++i) {
    const uint64_t av = (i < a.size()) ? a[i] : 0;
    const uint64_t bv = (i < b.size()) ? b[i] : 0;
    const u128 sum = static_cast<u128>(av) + bv + carry;
    res[i] = static_cast<uint64_t>(sum);
    carry = static_cast<uint64_t>(sum >> 64);
  }
  if (carry != 0) {
    res.push_back(carry);
  }
  TrimLimbs(res);
  return res;
}

// Requires a >= b.
inline std::vector<uint64_t> SubLimbs(std::span<const uint64_t> a,
                                      std::span<const uint64_t> b) {
  a = TrimSpan(a);
  b = TrimSpan(b);
  std::vector<uint64_t> res(a.size(), 0);
  uint64_t borrow = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const uint64_t bv = (i < b.size()) ? b[i] : 0;
    const u128 sub = static_cast<u128>(bv) + borrow;
    if (static_cast<u128>(a[i]) >= sub) {
      res[i] = static_cast<uint64_t>(static_cast<u128>(a[i]) - sub);
      borrow = 0;
    } else {
      res[i] = static_cast<uint64_t>((static_cast<u128>(1) << 64) + a[i] - sub);
      borrow = 1;
    }
  }
  TrimLimbs(res);
  return res;
}

inline void AddShiftedLimbsInPlace(std::vector<uint64_t>& dst,
                                   std::span<const uint64_t> src,
                                   std::size_t limb_offset) {
  src = TrimSpan(src);
  if (src.empty()) {
    return;
  }
  if (dst.size() < limb_offset + src.size()) {
    dst.resize(limb_offset + src.size(), 0);
  }
  uint64_t carry = 0;
  std::size_t i = 0;
  for (; i < src.size(); ++i) {
    const u128 sum = static_cast<u128>(dst[limb_offset + i]) + src[i] + carry;
    dst[limb_offset + i] = static_cast<uint64_t>(sum);
    carry = static_cast<uint64_t>(sum >> 64);
  }
  std::size_t pos = limb_offset + i;
  while (carry != 0) {
    if (pos == dst.size()) {
      dst.push_back(carry);
      break;
    }
    const u128 sum = static_cast<u128>(dst[pos]) + carry;
    dst[pos] = static_cast<uint64_t>(sum);
    carry = static_cast<uint64_t>(sum >> 64);
    ++pos;
  }
}

inline std::vector<uint64_t> MulLimbsSchoolbook(std::span<const uint64_t> a,
                                                std::span<const uint64_t> b) {
  a = TrimSpan(a);
  b = TrimSpan(b);
  if (a.empty() || b.empty()) {
    return {};
  }
  std::vector<uint64_t> res(a.size() + b.size(), 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    uint64_t carry = 0;
    const u128 ai = a[i];
    for (std::size_t j = 0; j < b.size(); ++j) {
      const u128 prod = ai * b[j] + res[i + j] + carry;
      res[i + j] = static_cast<uint64_t>(prod);
      carry = static_cast<uint64_t>(prod >> 64);
    }
    res[i + b.size()] = carry;
  }
  TrimLimbs(res);
  return res;
}

inline std::vector<uint64_t> MulLimbsKaratsuba(std::span<const uint64_t> a,
                                               std::span<const uint64_t> b) {
  a = TrimSpan(a);
  b = TrimSpan(b);
  if (a.empty() || b.empty()) {
    return {};
  }
  if (std::min(a.size(), b.size()) < kKaratsubaThreshold) {
    return MulLimbsSchoolbook(a, b);
  }

  // Handle highly unbalanced operands by chunking the larger operand.
  if (a.size() >= 2 * b.size()) {
    std::vector<uint64_t> res;
    res.reserve(a.size() + b.size());
    const std::size_t chunk = b.size();
    for (std::size_t offset = 0; offset < a.size(); offset += chunk) {
      const std::size_t len = std::min(chunk, a.size() - offset);
      auto part = MulLimbsKaratsuba(a.subspan(offset, len), b);
      AddShiftedLimbsInPlace(res, part, offset);
    }
    TrimLimbs(res);
    return res;
  }
  if (b.size() >= 2 * a.size()) {
    return MulLimbsKaratsuba(b, a);
  }

  const std::size_t m = std::max(a.size(), b.size()) / 2;
  const std::span<const uint64_t> a0 =
      TrimSpan(a.subspan(0, std::min(m, a.size())));
  const std::span<const uint64_t> a1 =
      (a.size() > m) ? TrimSpan(a.subspan(m)) : std::span<const uint64_t>{};
  const std::span<const uint64_t> b0 =
      TrimSpan(b.subspan(0, std::min(m, b.size())));
  const std::span<const uint64_t> b1 =
      (b.size() > m) ? TrimSpan(b.subspan(m)) : std::span<const uint64_t>{};

  auto z0 = MulLimbsKaratsuba(a0, b0);
  auto z2 = MulLimbsKaratsuba(a1, b1);

  auto a_sum = AddLimbs(a0, a1);
  auto b_sum = AddLimbs(b0, b1);
  auto z1 = MulLimbsKaratsuba(a_sum, b_sum);
  z1 = SubLimbs(z1, z0);
  z1 = SubLimbs(z1, z2);

  std::vector<uint64_t> res(a.size() + b.size(), 0);
  AddShiftedLimbsInPlace(res, z0, 0);
  AddShiftedLimbsInPlace(res, z1, m);
  AddShiftedLimbsInPlace(res, z2, 2 * m);
  TrimLimbs(res);
  return res;
}

// Single-limb division: returns {quotient, remainder}.
inline std::pair<std::vector<uint64_t>, uint64_t> DivRemSingleLimb(
    std::span<const uint64_t> u, uint64_t d) {
  if (d == 0) {
    throw std::domain_error("Division by zero");
  }
  u = TrimSpan(u);
  if (u.empty()) {
    return {{}, 0};
  }
  std::vector<uint64_t> q(u.size(), 0);
  uint64_t rem = 0;
  for (std::size_t i = u.size(); i-- > 0;) {
    const u128 cur = (static_cast<u128>(rem) << 64) | u[i];
    q[i] = static_cast<uint64_t>(cur / d);
    rem = static_cast<uint64_t>(cur % d);
  }
  TrimLimbs(q);
  return {std::move(q), rem};
}

// Knuth TAOCP Vol. 2, Section 4.3.1, Algorithm D for unsigned multi-limb
// division. Returns {quotient, remainder}.
inline std::pair<std::vector<uint64_t>, std::vector<uint64_t>> DivRemLimbs(
    std::span<const uint64_t> u_in, std::span<const uint64_t> v_in) {
  const auto u = TrimSpan(u_in);
  const auto v = TrimSpan(v_in);
  if (v.empty()) {
    throw std::domain_error("Division by zero");
  }
  if (u.empty()) {
    return {{}, {}};
  }
  const auto cmp = CmpLimbs(u, v);
  if (cmp == std::strong_ordering::less) {
    return {{}, std::vector<uint64_t>(u.begin(), u.end())};
  }
  if (cmp == std::strong_ordering::equal) {
    return {{1}, {}};
  }
  if (v.size() == 1) {
    auto [q, r] = DivRemSingleLimb(u, v[0]);
    std::vector<uint64_t> r_vec;
    if (r != 0) {
      r_vec.push_back(r);
    }
    return {std::move(q), std::move(r_vec)};
  }

  const std::size_t n = v.size();
  const std::size_t m = u.size() - n;

  // D1: Normalize so MSB of vn[n - 1] is 1.
  const unsigned s = static_cast<unsigned>(std::countl_zero(v[n - 1]));
  std::vector<uint64_t> vn(n, 0);
  std::vector<uint64_t> un(m + n + 1, 0);

  if (s == 0) {
    std::copy(v.begin(), v.end(), vn.begin());
    std::copy(u.begin(), u.end(), un.begin());
    un[m + n] = 0;
  } else {
    const unsigned rs = 64U - s;
    for (std::size_t i = n - 1; i > 0; --i) {
      vn[i] = (v[i] << s) | (v[i - 1] >> rs);
    }
    vn[0] = v[0] << s;

    un[m + n] = u[m + n - 1] >> rs;
    for (std::size_t i = m + n - 1; i > 0; --i) {
      un[i] = (u[i] << s) | (u[i - 1] >> rs);
    }
    un[0] = u[0] << s;
  }

  std::vector<uint64_t> q(m + 1, 0);
  const uint64_t vn1 = vn[n - 1];
  const uint64_t vn2 = vn[n - 2];

  // D2: Main loop from j = m down to 0.
  for (std::size_t j = m + 1; j-- > 0;) {
    // D3: Estimate quotient digit qhat and remainder rhat.
    const uint64_t u_hi = un[j + n];
    const uint64_t u_lo = un[j + n - 1];
    u128 qhat = 0;
    u128 rhat = 0;
    if (u_hi >= vn1) {
      qhat = std::numeric_limits<uint64_t>::max();
      rhat = static_cast<u128>(u_lo) + vn1;
    } else {
      const u128 num = (static_cast<u128>(u_hi) << 64) | u_lo;
      qhat = num / vn1;
      rhat = num % vn1;
    }

    while (rhat < (static_cast<u128>(1) << 64)) {
      const u128 lhs = qhat * vn2;
      const u128 rhs = (rhat << 64) | un[j + n - 2];
      if (lhs <= rhs) {
        break;
      }
      --qhat;
      rhat += vn1;
    }

    // D4: Multiply vn by qhat and subtract from un[j .. j + n].
    const uint64_t qh = static_cast<uint64_t>(qhat);
    uint64_t mul_carry = 0;
    uint64_t sub_borrow = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const u128 prod = static_cast<u128>(qh) * vn[i] + mul_carry;
      const uint64_t p_lo = static_cast<uint64_t>(prod);
      mul_carry = static_cast<uint64_t>(prod >> 64);

      const u128 subtrahend = static_cast<u128>(p_lo) + sub_borrow;
      if (static_cast<u128>(un[j + i]) >= subtrahend) {
        un[j + i] =
            static_cast<uint64_t>(static_cast<u128>(un[j + i]) - subtrahend);
        sub_borrow = 0;
      } else {
        un[j + i] = static_cast<uint64_t>((static_cast<u128>(1) << 64) +
                                          un[j + i] - subtrahend);
        sub_borrow = 1;
      }
    }

    const u128 final_sub = static_cast<u128>(mul_carry) + sub_borrow;
    const bool negative_rem = static_cast<u128>(un[j + n]) < final_sub;
    un[j + n] = static_cast<uint64_t>(static_cast<u128>(un[j + n]) - final_sub);

    // D5 & D6: Store quotient digit; if remainder went negative, add back.
    q[j] = qh;
    if (negative_rem) {
      --q[j];
      uint64_t add_carry = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const u128 sum = static_cast<u128>(un[j + i]) + vn[i] + add_carry;
        un[j + i] = static_cast<uint64_t>(sum);
        add_carry = static_cast<uint64_t>(sum >> 64);
      }
      un[j + n] += add_carry;
    }
  }

  // D8: Unnormalize remainder.
  std::vector<uint64_t> r(n, 0);
  if (s == 0) {
    std::copy(un.begin(), un.begin() + static_cast<std::ptrdiff_t>(n),
              r.begin());
  } else {
    const unsigned ls = 64U - s;
    for (std::size_t i = 0; i + 1 < n; ++i) {
      r[i] = (un[i] >> s) | (un[i + 1] << ls);
    }
    r[n - 1] = un[n - 1] >> s;
  }

  TrimLimbs(q);
  TrimLimbs(r);
  return {std::move(q), std::move(r)};
}

inline int HexCharToVal(char ch) {
  if (ch >= '0' && ch <= '9') {
    return ch - '0';
  }
  if (ch >= 'a' && ch <= 'f') {
    return 10 + (ch - 'a');
  }
  if (ch >= 'A' && ch <= 'F') {
    return 10 + (ch - 'A');
  }
  return -1;
}

inline std::string LimbsToHex(std::span<const uint64_t> limbs, bool negative,
                              bool prefix) {
  limbs = TrimSpan(limbs);
  if (limbs.empty()) {
    return prefix ? "0x0" : "0";
  }
  static constexpr char kHexDigits[] = "0123456789abcdef";
  std::string out;
  if (negative) {
    out.push_back('-');
  }
  if (prefix) {
    out += "0x";
  }

  // Most significant limb without leading zeros.
  const uint64_t top = limbs.back();
  bool started = false;
  for (int shift = 60; shift >= 0; shift -= 4) {
    const unsigned nibble = static_cast<unsigned>((top >> shift) & 0xFULL);
    if (nibble != 0 || started) {
      out.push_back(kHexDigits[nibble]);
      started = true;
    }
  }

  // Lower limbs padded to 16 hex digits each.
  for (std::size_t i = limbs.size() - 1; i-- > 0;) {
    const uint64_t w = limbs[i];
    for (int shift = 60; shift >= 0; shift -= 4) {
      out.push_back(kHexDigits[(w >> shift) & 0xFULL]);
    }
  }
  return out;
}

inline std::string LimbsToDecimal(std::span<const uint64_t> limbs,
                                  bool negative) {
  limbs = TrimSpan(limbs);
  if (limbs.empty()) {
    return "0";
  }
  std::vector<uint64_t> cur(limbs.begin(), limbs.end());
  std::vector<uint64_t> chunks;
  chunks.reserve(cur.size() + 1);
  while (!cur.empty()) {
    auto [q, r] = DivRemSingleLimb(cur, kDecBase);
    chunks.push_back(r);
    cur = std::move(q);
  }

  std::string out;
  if (negative) {
    out.push_back('-');
  }
  out += std::to_string(chunks.back());
  for (std::size_t i = chunks.size() - 1; i-- > 0;) {
    std::string block = std::to_string(chunks[i]);
    if (block.size() < kDecBaseDigits) {
      out.append(kDecBaseDigits - block.size(), '0');
    }
    out += block;
  }
  return out;
}

inline uint64_t ModMulU64(uint64_t a, uint64_t b, uint64_t mod) noexcept {
  return static_cast<uint64_t>((static_cast<u128>(a) * b) % mod);
}

inline uint64_t ModPowU64(uint64_t base, uint64_t exp, uint64_t mod) noexcept {
  uint64_t res = 1 % mod;
  base %= mod;
  while (exp > 0) {
    if ((exp & 1ULL) != 0) {
      res = ModMulU64(res, base, mod);
    }
    base = ModMulU64(base, base, mod);
    exp >>= 1;
  }
  return res;
}

inline uint64_t GcdU64(uint64_t a, uint64_t b) noexcept {
  while (b != 0) {
    const uint64_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

inline uint64_t ISqrtU64(uint64_t n) noexcept {
  if (n < 2) {
    return n;
  }
  const unsigned shift = (static_cast<unsigned>(std::bit_width(n)) + 1U) / 2U;
  uint64_t x = (shift >= 64) ? (uint64_t{1} << 63) : (uint64_t{1} << shift);
  // Ensure x >= ceil(sqrt(n)) even when shift == 64 (not possible for uint64_t
  // since bit_width <= 64 => shift <= 32).
  while (true) {
    const uint64_t y = (x + n / x) >> 1;
    if (y >= x) {
      return x;
    }
    x = y;
  }
}

inline bool IsPrimeU64(uint64_t n) noexcept {
  if (n < 2) {
    return false;
  }
  static constexpr uint64_t kSmallPrimes[] = {2,  3,  5,  7,  11, 13,
                                              17, 19, 23, 29, 31, 37};
  for (uint64_t p : kSmallPrimes) {
    if (n == p) {
      return true;
    }
    if (n % p == 0) {
      return false;
    }
  }
  if (n < 41ULL * 41ULL) {
    return true;
  }

  // Deterministic Miller-Rabin witnesses for all 64-bit integers (Sinclair).
  static constexpr uint64_t kWitnesses[] = {
      2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};

  uint64_t d = n - 1;
  const int s = std::countr_zero(d);
  d >>= s;

  for (uint64_t a : kWitnesses) {
    const uint64_t a_mod = a % n;
    if (a_mod == 0) {
      continue;
    }
    uint64_t x = ModPowU64(a_mod, d, n);
    if (x == 1 || x == n - 1) {
      continue;
    }
    bool composite = true;
    for (int r = 1; r < s; ++r) {
      x = ModMulU64(x, x, n);
      if (x == n - 1) {
        composite = false;
        break;
      }
      if (x == 1) {
        return false;
      }
    }
    if (composite) {
      return false;
    }
  }
  return true;
}

inline uint64_t PollardsRhoBrentU64(uint64_t n) noexcept {
  if ((n & 1ULL) == 0) {
    return 2;
  }
  constexpr uint64_t kBatch = 128;
  for (uint64_t attempt = 1; attempt <= 64; ++attempt) {
    const uint64_t c = attempt;
    uint64_t y = attempt + 1;
    uint64_t q = 1;
    uint64_t x = y;
    uint64_t ys = y;
    uint64_t r = 1;
    uint64_t g = 1;

    auto step = [n, c](uint64_t v) noexcept -> uint64_t {
      return static_cast<uint64_t>((static_cast<u128>(v) * v + c) % n);
    };

    while (g == 1) {
      x = y;
      for (uint64_t i = 0; i < r; ++i) {
        y = step(y);
      }
      uint64_t k = 0;
      while (k < r && g == 1) {
        ys = y;
        const uint64_t limit = std::min(kBatch, r - k);
        for (uint64_t i = 0; i < limit; ++i) {
          y = step(y);
          const uint64_t diff = (x >= y) ? (x - y) : (y - x);
          if (diff == 0) {
            g = n;
            break;
          }
          const uint64_t next_q = ModMulU64(q, diff, n);
          if (next_q == 0) {
            g = GcdU64(diff, n);
            break;
          }
          q = next_q;
        }
        if (g == 1) {
          g = GcdU64(q, n);
        }
        k += limit;
      }
      r <<= 1;
    }

    if (g == n) {
      for (uint64_t cnt = 0; cnt <= kBatch + 2; ++cnt) {
        ys = step(ys);
        const uint64_t diff = (x >= ys) ? (x - ys) : (ys - x);
        g = GcdU64(diff, n);
        if (g > 1) {
          break;
        }
      }
    }

    if (g > 1 && g < n) {
      return g;
    }
  }
  return n;
}

inline uint64_t SplitMix64(uint64_t& state) noexcept {
  uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

inline constexpr std::array<uint64_t, 54> kSmallPrimesTable = {
    2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,
    47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107,
    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
    191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};

}  // namespace detail

// ============================================================================
// Fixed-width unsigned big integer: BigUInt<Bits>
// ============================================================================

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
class BigUInt {
 public:
  static constexpr std::size_t kBits = Bits;
  static constexpr std::size_t kLimbs = Bits / 64;

  constexpr BigUInt() noexcept = default;

  constexpr BigUInt(uint64_t value) noexcept { limbs_[0] = value; }

  constexpr explicit BigUInt(const std::array<uint64_t, kLimbs>& limbs) noexcept
      : limbs_(limbs) {}

  static constexpr BigUInt Zero() noexcept { return BigUInt{}; }

  static constexpr BigUInt One() noexcept { return BigUInt(1ULL); }

  static BigUInt FromHex(std::string_view hex) {
    if (hex.empty()) {
      throw std::invalid_argument("Empty hex string in BigUInt::FromHex");
    }
    if (hex.size() >= 2 && hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) {
      hex.remove_prefix(2);
    }
    if (hex.empty()) {
      throw std::invalid_argument("Empty hex digits in BigUInt::FromHex");
    }
    for (char ch : hex) {
      if (detail::HexCharToVal(ch) < 0) {
        throw std::invalid_argument(
            "Invalid hex character in BigUInt::FromHex");
      }
    }
    while (hex.size() > 1 && hex.front() == '0') {
      hex.remove_prefix(1);
    }
    if (hex.size() > Bits / 4) {
      throw std::overflow_error("Hex value overflows BigUInt<Bits>");
    }

    BigUInt result;
    std::size_t bit_pos = 0;
    for (std::size_t i = hex.size(); i-- > 0;) {
      const uint64_t nibble =
          static_cast<uint64_t>(detail::HexCharToVal(hex[i]));
      result.limbs_[bit_pos / 64] |= (nibble << (bit_pos % 64));
      bit_pos += 4;
    }
    return result;
  }

  static BigUInt FromDecimal(std::string_view dec) {
    if (dec.empty()) {
      throw std::invalid_argument(
          "Empty decimal string in BigUInt::FromDecimal");
    }
    for (char ch : dec) {
      if (ch < '0' || ch > '9') {
        throw std::invalid_argument(
            "Invalid decimal character in BigUInt::FromDecimal");
      }
    }
    while (dec.size() > 1 && dec.front() == '0') {
      dec.remove_prefix(1);
    }

    BigUInt result;
    std::size_t pos = 0;
    while (pos < dec.size()) {
      const std::size_t chunk_len =
          std::min(detail::kDecBaseDigits, dec.size() - pos);
      uint64_t chunk_val = 0;
      uint64_t mult = 1;
      for (std::size_t k = 0; k < chunk_len; ++k) {
        chunk_val =
            chunk_val * 10ULL + static_cast<uint64_t>(dec[pos + k] - '0');
        mult *= 10ULL;
      }

      // Multiply result by mult and add chunk_val, checking for overflow.
      uint64_t carry = chunk_val;
      for (std::size_t i = 0; i < kLimbs; ++i) {
        const detail::u128 prod =
            static_cast<detail::u128>(result.limbs_[i]) * mult + carry;
        result.limbs_[i] = static_cast<uint64_t>(prod);
        carry = static_cast<uint64_t>(prod >> 64);
      }
      if (carry != 0) {
        throw std::overflow_error("Decimal value overflows BigUInt<Bits>");
      }
      pos += chunk_len;
    }
    return result;
  }

  constexpr bool IsZero() const noexcept {
    uint64_t acc = 0;
    for (uint64_t w : limbs_) {
      acc |= w;
    }
    return acc == 0;
  }

  constexpr bool IsOne() const noexcept {
    uint64_t diff = limbs_[0] ^ 1ULL;
    for (std::size_t i = 1; i < kLimbs; ++i) {
      diff |= limbs_[i];
    }
    return diff == 0;
  }

  constexpr bool IsOdd() const noexcept { return (limbs_[0] & 1ULL) != 0; }

  constexpr bool IsEven() const noexcept { return (limbs_[0] & 1ULL) == 0; }

  constexpr std::size_t BitLength() const noexcept {
    for (std::size_t i = kLimbs; i-- > 0;) {
      if (limbs_[i] != 0) {
        return i * 64 + static_cast<std::size_t>(std::bit_width(limbs_[i]));
      }
    }
    return 0;
  }

  constexpr bool GetBit(std::size_t index) const noexcept {
    if (index >= Bits) {
      return false;
    }
    return ((limbs_[index / 64] >> (index % 64)) & 1ULL) != 0;
  }

  constexpr void SetBit(std::size_t index) {
    if (index >= Bits) {
      throw std::out_of_range("Bit index out of range in BigUInt::SetBit");
    }
    limbs_[index / 64] |= (uint64_t{1} << (index % 64));
  }

  constexpr uint64_t Limb(std::size_t index) const {
    if (index >= kLimbs) {
      throw std::out_of_range("Limb index out of range in BigUInt::Limb");
    }
    return limbs_[index];
  }

  constexpr const std::array<uint64_t, kLimbs>& Limbs() const noexcept {
    return limbs_;
  }

  std::string ToHex(bool prefix = false) const {
    return detail::LimbsToHex(limbs_, false, prefix);
  }

  std::string ToDecimal() const {
    return detail::LimbsToDecimal(limbs_, false);
  }

  std::string ToString() const { return ToDecimal(); }

  static constexpr std::pair<BigUInt, uint64_t> AddWithCarry(
      const BigUInt& a, const BigUInt& b, uint64_t carry_in = 0) noexcept {
    BigUInt res;
    uint64_t carry = carry_in;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const detail::u128 sum =
          static_cast<detail::u128>(a.limbs_[i]) + b.limbs_[i] + carry;
      res.limbs_[i] = static_cast<uint64_t>(sum);
      carry = static_cast<uint64_t>(sum >> 64);
    }
    return {res, carry};
  }

  static constexpr std::pair<BigUInt, uint64_t> SubWithBorrow(
      const BigUInt& a, const BigUInt& b, uint64_t borrow_in = 0) noexcept {
    BigUInt res;
    uint64_t borrow = borrow_in;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const detail::u128 sub = static_cast<detail::u128>(b.limbs_[i]) + borrow;
      if (static_cast<detail::u128>(a.limbs_[i]) >= sub) {
        res.limbs_[i] =
            static_cast<uint64_t>(static_cast<detail::u128>(a.limbs_[i]) - sub);
        borrow = 0;
      } else {
        res.limbs_[i] = static_cast<uint64_t>(
            (static_cast<detail::u128>(1) << 64) + a.limbs_[i] - sub);
        borrow = 1;
      }
    }
    return {res, borrow};
  }

  constexpr BigUInt<2 * Bits> MulWide(const BigUInt& rhs) const noexcept {
    std::array<uint64_t, 2 * kLimbs> prod{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      uint64_t carry = 0;
      const detail::u128 ai = limbs_[i];
      for (std::size_t j = 0; j < kLimbs; ++j) {
        const detail::u128 t = ai * rhs.limbs_[j] + prod[i + j] + carry;
        prod[i + j] = static_cast<uint64_t>(t);
        carry = static_cast<uint64_t>(t >> 64);
      }
      prod[i + kLimbs] = carry;
    }
    return BigUInt<2 * Bits>(prod);
  }

  // Separated Operand Scanning (SOS) wide squaring.
  constexpr BigUInt<2 * Bits> SqrWide() const noexcept {
    std::array<uint64_t, 2 * kLimbs> res{};

    // 1. Off-diagonal products (i < j).
    for (std::size_t i = 0; i < kLimbs; ++i) {
      uint64_t carry = 0;
      const detail::u128 ai = limbs_[i];
      for (std::size_t j = i + 1; j < kLimbs; ++j) {
        const detail::u128 t = ai * limbs_[j] + res[i + j] + carry;
        res[i + j] = static_cast<uint64_t>(t);
        carry = static_cast<uint64_t>(t >> 64);
      }
      res[i + kLimbs] = carry;
    }

    // 2. Double the off-diagonal sum (left shift by 1 bit).
    uint64_t shift_carry = 0;
    for (std::size_t i = 0; i < 2 * kLimbs; ++i) {
      const uint64_t next_carry = res[i] >> 63;
      res[i] = (res[i] << 1) | shift_carry;
      shift_carry = next_carry;
    }

    // 3. Add diagonal squares limbs_[i]^2.
    detail::u128 carry = 0;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const detail::u128 sq = static_cast<detail::u128>(limbs_[i]) * limbs_[i];
      const detail::u128 lo = static_cast<detail::u128>(res[2 * i]) +
                              static_cast<uint64_t>(sq) + carry;
      res[2 * i] = static_cast<uint64_t>(lo);
      const detail::u128 hi = static_cast<detail::u128>(res[2 * i + 1]) +
                              static_cast<uint64_t>(sq >> 64) + (lo >> 64);
      res[2 * i + 1] = static_cast<uint64_t>(hi);
      carry = hi >> 64;
    }

    return BigUInt<2 * Bits>(res);
  }

  static std::pair<BigUInt, BigUInt> DivRem(const BigUInt& a,
                                            const BigUInt& b) {
    if (b.IsZero()) {
      throw std::domain_error("Division by zero in BigUInt::DivRem");
    }
    if (a < b) {
      return {Zero(), a};
    }
    if (a == b) {
      return {One(), Zero()};
    }

    // Fast path for single-limb divisor.
    bool single_limb = true;
    for (std::size_t i = 1; i < kLimbs; ++i) {
      if (b.limbs_[i] != 0) {
        single_limb = false;
        break;
      }
    }
    if (single_limb) {
      BigUInt q;
      uint64_t rem = 0;
      const uint64_t d = b.limbs_[0];
      for (std::size_t i = kLimbs; i-- > 0;) {
        const detail::u128 cur =
            (static_cast<detail::u128>(rem) << 64) | a.limbs_[i];
        q.limbs_[i] = static_cast<uint64_t>(cur / d);
        rem = static_cast<uint64_t>(cur % d);
      }
      return {q, BigUInt(rem)};
    }

    auto [qv, rv] = detail::DivRemLimbs(a.limbs_, b.limbs_);
    BigUInt q;
    BigUInt r;
    for (std::size_t i = 0; i < qv.size(); ++i) {
      q.limbs_[i] = qv[i];
    }
    for (std::size_t i = 0; i < rv.size(); ++i) {
      r.limbs_[i] = rv[i];
    }
    return {q, r};
  }

  // Constant-time conditional select: returns a if choose_a, else b.
  static constexpr BigUInt CtSelect(bool choose_a, const BigUInt& a,
                                    const BigUInt& b) noexcept {
    const uint64_t mask = uint64_t{0} - static_cast<uint64_t>(choose_a);
    BigUInt out;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      out.limbs_[i] = b.limbs_[i] ^ ((a.limbs_[i] ^ b.limbs_[i]) & mask);
    }
    return out;
  }

  // Constant-time conditional swap.
  static constexpr void CtSwap(bool swap, BigUInt& a, BigUInt& b) noexcept {
    const uint64_t mask = uint64_t{0} - static_cast<uint64_t>(swap);
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const uint64_t diff = (a.limbs_[i] ^ b.limbs_[i]) & mask;
      a.limbs_[i] ^= diff;
      b.limbs_[i] ^= diff;
    }
  }

  constexpr BigUInt operator+(const BigUInt& rhs) const noexcept {
    return AddWithCarry(*this, rhs).first;
  }

  constexpr BigUInt operator-(const BigUInt& rhs) const noexcept {
    return SubWithBorrow(*this, rhs).first;
  }

  constexpr BigUInt operator*(const BigUInt& rhs) const noexcept {
    BigUInt res;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      uint64_t carry = 0;
      const detail::u128 ai = limbs_[i];
      for (std::size_t j = 0; j < kLimbs - i; ++j) {
        const detail::u128 prod =
            ai * rhs.limbs_[j] + res.limbs_[i + j] + carry;
        res.limbs_[i + j] = static_cast<uint64_t>(prod);
        carry = static_cast<uint64_t>(prod >> 64);
      }
    }
    return res;
  }

  BigUInt operator/(const BigUInt& rhs) const {
    return DivRem(*this, rhs).first;
  }

  BigUInt operator%(const BigUInt& rhs) const {
    return DivRem(*this, rhs).second;
  }

  constexpr BigUInt& operator+=(const BigUInt& rhs) noexcept {
    *this = *this + rhs;
    return *this;
  }

  constexpr BigUInt& operator-=(const BigUInt& rhs) noexcept {
    *this = *this - rhs;
    return *this;
  }

  constexpr BigUInt& operator*=(const BigUInt& rhs) noexcept {
    *this = *this * rhs;
    return *this;
  }

  BigUInt& operator/=(const BigUInt& rhs) {
    *this = *this / rhs;
    return *this;
  }

  BigUInt& operator%=(const BigUInt& rhs) {
    *this = *this % rhs;
    return *this;
  }

  constexpr BigUInt operator&(const BigUInt& rhs) const noexcept {
    BigUInt res;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      res.limbs_[i] = limbs_[i] & rhs.limbs_[i];
    }
    return res;
  }

  constexpr BigUInt operator|(const BigUInt& rhs) const noexcept {
    BigUInt res;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      res.limbs_[i] = limbs_[i] | rhs.limbs_[i];
    }
    return res;
  }

  constexpr BigUInt operator^(const BigUInt& rhs) const noexcept {
    BigUInt res;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      res.limbs_[i] = limbs_[i] ^ rhs.limbs_[i];
    }
    return res;
  }

  constexpr BigUInt operator~() const noexcept {
    BigUInt res;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      res.limbs_[i] = ~limbs_[i];
    }
    return res;
  }

  constexpr BigUInt& operator&=(const BigUInt& rhs) noexcept {
    for (std::size_t i = 0; i < kLimbs; ++i) {
      limbs_[i] &= rhs.limbs_[i];
    }
    return *this;
  }

  constexpr BigUInt& operator|=(const BigUInt& rhs) noexcept {
    for (std::size_t i = 0; i < kLimbs; ++i) {
      limbs_[i] |= rhs.limbs_[i];
    }
    return *this;
  }

  constexpr BigUInt& operator^=(const BigUInt& rhs) noexcept {
    for (std::size_t i = 0; i < kLimbs; ++i) {
      limbs_[i] ^= rhs.limbs_[i];
    }
    return *this;
  }

  constexpr BigUInt operator<<(std::size_t shift) const noexcept {
    if (shift == 0) {
      return *this;
    }
    if (shift >= Bits) {
      return Zero();
    }
    const std::size_t word_shift = shift / 64;
    const unsigned bit_shift = static_cast<unsigned>(shift % 64);
    BigUInt res;
    if (bit_shift == 0) {
      for (std::size_t i = word_shift; i < kLimbs; ++i) {
        res.limbs_[i] = limbs_[i - word_shift];
      }
    } else {
      const unsigned rshift = 64U - bit_shift;
      res.limbs_[word_shift] = limbs_[0] << bit_shift;
      for (std::size_t i = word_shift + 1; i < kLimbs; ++i) {
        res.limbs_[i] = (limbs_[i - word_shift] << bit_shift) |
                        (limbs_[i - word_shift - 1] >> rshift);
      }
    }
    return res;
  }

  constexpr BigUInt operator>>(std::size_t shift) const noexcept {
    if (shift == 0) {
      return *this;
    }
    if (shift >= Bits) {
      return Zero();
    }
    const std::size_t word_shift = shift / 64;
    const unsigned bit_shift = static_cast<unsigned>(shift % 64);
    BigUInt res;
    if (bit_shift == 0) {
      for (std::size_t i = 0; i + word_shift < kLimbs; ++i) {
        res.limbs_[i] = limbs_[i + word_shift];
      }
    } else {
      const unsigned lshift = 64U - bit_shift;
      const std::size_t last = kLimbs - 1 - word_shift;
      for (std::size_t i = 0; i < last; ++i) {
        res.limbs_[i] = (limbs_[i + word_shift] >> bit_shift) |
                        (limbs_[i + word_shift + 1] << lshift);
      }
      res.limbs_[last] = limbs_[kLimbs - 1] >> bit_shift;
    }
    return res;
  }

  constexpr BigUInt& operator<<=(std::size_t shift) noexcept {
    *this = *this << shift;
    return *this;
  }

  constexpr BigUInt& operator>>=(std::size_t shift) noexcept {
    *this = *this >> shift;
    return *this;
  }

  constexpr bool operator==(const BigUInt& rhs) const noexcept = default;

  constexpr std::strong_ordering operator<=>(
      const BigUInt& rhs) const noexcept {
    for (std::size_t i = kLimbs; i-- > 0;) {
      if (limbs_[i] != rhs.limbs_[i]) {
        return limbs_[i] <=> rhs.limbs_[i];
      }
    }
    return std::strong_ordering::equal;
  }

  friend std::ostream& operator<<(std::ostream& os, const BigUInt& val) {
    return os << val.ToDecimal();
  }

 private:
  std::array<uint64_t, kLimbs> limbs_{};
};

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline std::pair<BigUInt<Bits>, BigUInt<Bits>> DivRem(const BigUInt<Bits>& a,
                                                      const BigUInt<Bits>& b) {
  return BigUInt<Bits>::DivRem(a, b);
}

// ============================================================================
// Free modular arithmetic helpers on BigUInt<Bits>
// ============================================================================

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigUInt<Bits> ModAdd(const BigUInt<Bits>& a, const BigUInt<Bits>& b,
                            const BigUInt<Bits>& mod) {
  if (mod.IsZero()) {
    throw std::domain_error("Modulus must be nonzero in ModAdd");
  }
  const BigUInt<Bits> ar = (a >= mod) ? (a % mod) : a;
  const BigUInt<Bits> br = (b >= mod) ? (b % mod) : b;
  const auto [sum, carry] = BigUInt<Bits>::AddWithCarry(ar, br);
  const auto [diff, borrow] = BigUInt<Bits>::SubWithBorrow(sum, mod);
  const bool need_sub = (carry != 0) || (borrow == 0);
  return BigUInt<Bits>::CtSelect(need_sub, diff, sum);
}

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigUInt<Bits> ModSub(const BigUInt<Bits>& a, const BigUInt<Bits>& b,
                            const BigUInt<Bits>& mod) {
  if (mod.IsZero()) {
    throw std::domain_error("Modulus must be nonzero in ModSub");
  }
  const BigUInt<Bits> ar = (a >= mod) ? (a % mod) : a;
  const BigUInt<Bits> br = (b >= mod) ? (b % mod) : b;
  const auto [diff, borrow] = BigUInt<Bits>::SubWithBorrow(ar, br);
  const auto [added, unused_carry] = BigUInt<Bits>::AddWithCarry(diff, mod);
  (void)unused_carry;
  return BigUInt<Bits>::CtSelect(borrow != 0, added, diff);
}

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigUInt<Bits> ModMul(const BigUInt<Bits>& a, const BigUInt<Bits>& b,
                            const BigUInt<Bits>& mod) {
  if (mod.IsZero()) {
    throw std::domain_error("Modulus must be nonzero in ModMul");
  }
  if (mod.IsOne()) {
    return BigUInt<Bits>::Zero();
  }
  const BigUInt<Bits> ar = (a >= mod) ? (a % mod) : a;
  const BigUInt<Bits> br = (b >= mod) ? (b % mod) : b;
  const BigUInt<2 * Bits> prod = ar.MulWide(br);

  std::array<uint64_t, 2 * BigUInt<Bits>::kLimbs> mod_wide_limbs{};
  for (std::size_t i = 0; i < BigUInt<Bits>::kLimbs; ++i) {
    mod_wide_limbs[i] = mod.Limbs()[i];
  }
  const BigUInt<2 * Bits> mod_wide(mod_wide_limbs);
  const BigUInt<2 * Bits> rem_wide =
      BigUInt<2 * Bits>::DivRem(prod, mod_wide).second;

  std::array<uint64_t, BigUInt<Bits>::kLimbs> out{};
  for (std::size_t i = 0; i < BigUInt<Bits>::kLimbs; ++i) {
    out[i] = rem_wide.Limbs()[i];
  }
  return BigUInt<Bits>(out);
}

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigUInt<Bits> ModInv(const BigUInt<Bits>& a, const BigUInt<Bits>& mod) {
  if (mod.IsZero() || mod.IsOne()) {
    throw std::domain_error("Modulus must be greater than 1 in ModInv");
  }
  BigUInt<Bits> r = (a >= mod) ? (a % mod) : a;
  if (r.IsZero()) {
    throw std::domain_error("Modular inverse does not exist (gcd != 1)");
  }

  // Unsigned extended Euclidean algorithm with alternating-sign tracking.
  // Invariant: old_u * r + u * old_r == mod, so u_next < mod whenever rem > 0.
  BigUInt<Bits> old_r = mod;
  BigUInt<Bits> old_u = BigUInt<Bits>::Zero();
  BigUInt<Bits> u = BigUInt<Bits>::One();
  bool u_negative = false;

  while (true) {
    auto [q, rem] = BigUInt<Bits>::DivRem(old_r, r);
    if (rem.IsZero()) {
      break;
    }
    BigUInt<Bits> next_u = old_u + q * u;
    old_r = r;
    r = rem;
    old_u = u;
    u = next_u;
    u_negative = !u_negative;
  }

  if (!r.IsOne()) {
    throw std::domain_error("Modular inverse does not exist (gcd != 1)");
  }
  return u_negative ? (mod - u) : u;
}

// ============================================================================
// Montgomery arithmetic context for odd modulus on BigUInt<Bits>
// ============================================================================

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
class MontgomeryCtx {
 public:
  static constexpr std::size_t kBits = Bits;
  static constexpr std::size_t kLimbs = Bits / 64;

  explicit MontgomeryCtx(const BigUInt<Bits>& modulus) : modulus_(modulus) {
    if (modulus_.IsZero() || modulus_.IsOne() || modulus_.IsEven()) {
      throw std::invalid_argument("MontgomeryCtx requires an odd modulus > 1");
    }

    // Compute n0_inv = -modulus[0]^{-1} mod 2^64 via Newton-Raphson.
    const uint64_t n0 = modulus_.Limbs()[0];
    uint64_t inv = 1;
    for (int i = 0; i < 6; ++i) {
      inv *= (uint64_t{2} - n0 * inv);
    }
    n0_inv_ = uint64_t{0} - inv;

    // Precompute r_mod = 2^Bits mod modulus and r2_mod = 2^(2*Bits) mod
    // modulus.
    std::array<uint64_t, 2 * kLimbs> mod_wide_limbs{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      mod_wide_limbs[i] = modulus_.Limbs()[i];
    }
    const BigUInt<2 * Bits> mod_wide(mod_wide_limbs);

    BigUInt<2 * Bits> r_wide;
    r_wide.SetBit(Bits);
    const BigUInt<2 * Bits> r_mod_wide =
        BigUInt<2 * Bits>::DivRem(r_wide, mod_wide).second;
    std::array<uint64_t, kLimbs> r_limbs{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      r_limbs[i] = r_mod_wide.Limbs()[i];
    }
    r_mod_ = BigUInt<Bits>(r_limbs);

    const BigUInt<2 * Bits> r2_wide =
        BigUInt<2 * Bits>::DivRem(r_mod_.SqrWide(), mod_wide).second;
    std::array<uint64_t, kLimbs> r2_limbs{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      r2_limbs[i] = r2_wide.Limbs()[i];
    }
    r2_mod_ = BigUInt<Bits>(r2_limbs);
  }

  const BigUInt<Bits>& Modulus() const noexcept { return modulus_; }
  const BigUInt<Bits>& RMod() const noexcept { return r_mod_; }
  const BigUInt<Bits>& R2Mod() const noexcept { return r2_mod_; }

  BigUInt<Bits> ToMont(const BigUInt<Bits>& a) const {
    const BigUInt<Bits> ar = (a >= modulus_) ? (a % modulus_) : a;
    return Mul(ar, r2_mod_);
  }

  BigUInt<Bits> FromMont(const BigUInt<Bits>& a_mont) const noexcept {
    // Pass One() as first operand so CIOS invariant t < N + 1 <= 2N always
    // holds.
    return Mul(BigUInt<Bits>::One(), a_mont);
  }

  BigUInt<Bits> Add(const BigUInt<Bits>& a,
                    const BigUInt<Bits>& b) const noexcept {
    const auto [sum, carry] = BigUInt<Bits>::AddWithCarry(a, b);
    const auto [diff, borrow] = BigUInt<Bits>::SubWithBorrow(sum, modulus_);
    const bool need_sub = (carry != 0) || (borrow == 0);
    return BigUInt<Bits>::CtSelect(need_sub, diff, sum);
  }

  BigUInt<Bits> Sub(const BigUInt<Bits>& a,
                    const BigUInt<Bits>& b) const noexcept {
    const auto [diff, borrow] = BigUInt<Bits>::SubWithBorrow(a, b);
    const auto [added, unused_carry] =
        BigUInt<Bits>::AddWithCarry(diff, modulus_);
    (void)unused_carry;
    return BigUInt<Bits>::CtSelect(borrow != 0, added, diff);
  }

  BigUInt<Bits> Neg(const BigUInt<Bits>& a) const noexcept {
    const auto [diff, unused_borrow] =
        BigUInt<Bits>::SubWithBorrow(modulus_, a);
    (void)unused_borrow;
    return BigUInt<Bits>::CtSelect(a.IsZero(), BigUInt<Bits>::Zero(), diff);
  }

  // CIOS (Coarsely Integrated Operand Scanning) Montgomery multiplication.
  // Assumes a < modulus_ and b < 2^Bits (canonical inputs have a, b <
  // modulus_).
  BigUInt<Bits> Mul(const BigUInt<Bits>& a,
                    const BigUInt<Bits>& b) const noexcept {
    std::array<uint64_t, kLimbs + 2> t{};
    const auto& al = a.Limbs();
    const auto& bl = b.Limbs();
    const auto& nl = modulus_.Limbs();

    for (std::size_t i = 0; i < kLimbs; ++i) {
      // Multiply-accumulate: t += a * b[i].
      uint64_t c = 0;
      const detail::u128 bi = bl[i];
      for (std::size_t j = 0; j < kLimbs; ++j) {
        const detail::u128 prod =
            static_cast<detail::u128>(al[j]) * bi + t[j] + c;
        t[j] = static_cast<uint64_t>(prod);
        c = static_cast<uint64_t>(prod >> 64);
      }
      const detail::u128 sum1 = static_cast<detail::u128>(t[kLimbs]) + c;
      t[kLimbs] = static_cast<uint64_t>(sum1);
      t[kLimbs + 1] = static_cast<uint64_t>(sum1 >> 64);

      // Montgomery reduction step: t = (t + m * N) / 2^64.
      const uint64_t m = t[0] * n0_inv_;
      const detail::u128 prod0 = static_cast<detail::u128>(m) * nl[0] + t[0];
      uint64_t c2 = static_cast<uint64_t>(prod0 >> 64);
      for (std::size_t j = 1; j < kLimbs; ++j) {
        const detail::u128 prod =
            static_cast<detail::u128>(m) * nl[j] + t[j] + c2;
        t[j - 1] = static_cast<uint64_t>(prod);
        c2 = static_cast<uint64_t>(prod >> 64);
      }
      const detail::u128 sum2 = static_cast<detail::u128>(t[kLimbs]) + c2;
      t[kLimbs - 1] = static_cast<uint64_t>(sum2);
      t[kLimbs] = t[kLimbs + 1] + static_cast<uint64_t>(sum2 >> 64);
      t[kLimbs + 1] = 0;
    }

    std::array<uint64_t, kLimbs> res_limbs{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      res_limbs[i] = t[i];
    }
    const BigUInt<Bits> res(res_limbs);
    const auto [diff, borrow] = BigUInt<Bits>::SubWithBorrow(res, modulus_);
    const bool need_sub = (t[kLimbs] != 0) || (borrow == 0);
    return BigUInt<Bits>::CtSelect(need_sub, diff, res);
  }

  // Montgomery reduction (REDC) of a 2*Bits product T < N * R.
  BigUInt<Bits> Reduce(const BigUInt<2 * Bits>& wide) const noexcept {
    std::array<uint64_t, 2 * kLimbs + 1> t{};
    for (std::size_t i = 0; i < 2 * kLimbs; ++i) {
      t[i] = wide.Limbs()[i];
    }
    const auto& nl = modulus_.Limbs();

    for (std::size_t i = 0; i < kLimbs; ++i) {
      const uint64_t m = t[i] * n0_inv_;
      uint64_t carry = 0;
      for (std::size_t j = 0; j < kLimbs; ++j) {
        const detail::u128 prod =
            static_cast<detail::u128>(m) * nl[j] + t[i + j] + carry;
        t[i + j] = static_cast<uint64_t>(prod);
        carry = static_cast<uint64_t>(prod >> 64);
      }
      std::size_t idx = i + kLimbs;
      while (carry != 0) {
        const detail::u128 s = static_cast<detail::u128>(t[idx]) + carry;
        t[idx] = static_cast<uint64_t>(s);
        carry = static_cast<uint64_t>(s >> 64);
        ++idx;
      }
    }

    std::array<uint64_t, kLimbs> res_limbs{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      res_limbs[i] = t[i + kLimbs];
    }
    const BigUInt<Bits> res(res_limbs);
    const auto [diff, borrow] = BigUInt<Bits>::SubWithBorrow(res, modulus_);
    const bool need_sub = (t[2 * kLimbs] != 0) || (borrow == 0);
    return BigUInt<Bits>::CtSelect(need_sub, diff, res);
  }

  BigUInt<Bits> Sqr(const BigUInt<Bits>& a) const noexcept {
    return Reduce(a.SqrWide());
  }

  BigUInt<Bits> Pow(const BigUInt<Bits>& base_mont,
                    const BigUInt<Bits>& exp) const noexcept {
    BigUInt<Bits> result = r_mod_;
    const std::size_t bits = exp.BitLength();
    for (std::size_t i = bits; i-- > 0;) {
      result = Sqr(result);
      if (exp.GetBit(i)) {
        result = Mul(result, base_mont);
      }
    }
    return result;
  }

  BigUInt<Bits> Inv(const BigUInt<Bits>& a_mont) const {
    const BigUInt<Bits> a = FromMont(a_mont);
    const BigUInt<Bits> a_inv = ModInv(a, modulus_);
    return ToMont(a_inv);
  }

 private:
  BigUInt<Bits> modulus_;
  uint64_t n0_inv_ = 0;
  BigUInt<Bits> r_mod_;
  BigUInt<Bits> r2_mod_;
};

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigUInt<Bits> ModPow(const BigUInt<Bits>& base, const BigUInt<Bits>& exp,
                            const BigUInt<Bits>& mod) {
  if (mod.IsZero()) {
    throw std::domain_error("Modulus must be nonzero in ModPow");
  }
  if (mod.IsOne()) {
    return BigUInt<Bits>::Zero();
  }
  if (mod.IsOdd()) {
    MontgomeryCtx<Bits> ctx(mod);
    return ctx.FromMont(ctx.Pow(ctx.ToMont(base), exp));
  }
  BigUInt<Bits> result = BigUInt<Bits>::One();
  BigUInt<Bits> b = (base >= mod) ? (base % mod) : base;
  const std::size_t bits = exp.BitLength();
  for (std::size_t i = bits; i-- > 0;) {
    result = ModMul(result, result, mod);
    if (exp.GetBit(i)) {
      result = ModMul(result, b, mod);
    }
  }
  return result;
}

// ============================================================================
// Arbitrary-precision signed integer: BigInt
// ============================================================================

class BigInt {
 public:
  BigInt() noexcept = default;

  template <std::signed_integral T>
  BigInt(T value) : BigInt(FromInt64(static_cast<int64_t>(value))) {}

  template <std::unsigned_integral T>
    requires(!std::same_as<T, bool>)
  BigInt(T value) : BigInt(FromUint64(static_cast<uint64_t>(value))) {}

  static BigInt FromInt64(int64_t value) {
    BigInt res;
    if (value == 0) {
      return res;
    }
    if (value < 0) {
      res.negative_ = true;
      res.limbs_.push_back(uint64_t{0} - static_cast<uint64_t>(value));
    } else {
      res.limbs_.push_back(static_cast<uint64_t>(value));
    }
    return res;
  }

  static BigInt FromUint64(uint64_t value) {
    BigInt res;
    if (value != 0) {
      res.limbs_.push_back(value);
    }
    return res;
  }

  static BigInt FromDecimal(std::string_view dec) {
    if (dec.empty()) {
      throw std::invalid_argument(
          "Empty decimal string in BigInt::FromDecimal");
    }
    bool neg = false;
    if (dec.front() == '-') {
      neg = true;
      dec.remove_prefix(1);
    } else if (dec.front() == '+') {
      dec.remove_prefix(1);
    }
    if (dec.empty()) {
      throw std::invalid_argument("Sign-only string in BigInt::FromDecimal");
    }
    for (char ch : dec) {
      if (ch < '0' || ch > '9') {
        throw std::invalid_argument(
            "Invalid decimal character in BigInt::FromDecimal");
      }
    }
    while (dec.size() > 1 && dec.front() == '0') {
      dec.remove_prefix(1);
    }

    BigInt res;
    std::size_t pos = 0;
    while (pos < dec.size()) {
      const std::size_t chunk_len =
          std::min(detail::kDecBaseDigits, dec.size() - pos);
      uint64_t chunk_val = 0;
      uint64_t mult = 1;
      for (std::size_t k = 0; k < chunk_len; ++k) {
        chunk_val =
            chunk_val * 10ULL + static_cast<uint64_t>(dec[pos + k] - '0');
        mult *= 10ULL;
      }

      uint64_t carry = chunk_val;
      for (uint64_t& limb : res.limbs_) {
        const detail::u128 prod =
            static_cast<detail::u128>(limb) * mult + carry;
        limb = static_cast<uint64_t>(prod);
        carry = static_cast<uint64_t>(prod >> 64);
      }
      if (carry != 0) {
        res.limbs_.push_back(carry);
      }
      pos += chunk_len;
    }
    res.negative_ = neg;
    res.Normalize();
    return res;
  }

  static BigInt FromHex(std::string_view hex) {
    if (hex.empty()) {
      throw std::invalid_argument("Empty hex string in BigInt::FromHex");
    }
    bool neg = false;
    if (hex.front() == '-') {
      neg = true;
      hex.remove_prefix(1);
    } else if (hex.front() == '+') {
      hex.remove_prefix(1);
    }
    if (hex.size() >= 2 && hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) {
      hex.remove_prefix(2);
    }
    if (hex.empty()) {
      throw std::invalid_argument("Empty hex digits in BigInt::FromHex");
    }
    for (char ch : hex) {
      if (detail::HexCharToVal(ch) < 0) {
        throw std::invalid_argument("Invalid hex character in BigInt::FromHex");
      }
    }
    while (hex.size() > 1 && hex.front() == '0') {
      hex.remove_prefix(1);
    }

    BigInt res;
    const std::size_t num_limbs = (hex.size() + 15) / 16;
    res.limbs_.assign(num_limbs, 0);
    std::size_t bit_pos = 0;
    for (std::size_t i = hex.size(); i-- > 0;) {
      const uint64_t nibble =
          static_cast<uint64_t>(detail::HexCharToVal(hex[i]));
      res.limbs_[bit_pos / 64] |= (nibble << (bit_pos % 64));
      bit_pos += 4;
    }
    res.negative_ = neg;
    res.Normalize();
    return res;
  }

  template <std::size_t Bits>
    requires(Bits >= 64 && Bits % 64 == 0)
  static BigInt FromBigUInt(const BigUInt<Bits>& value) {
    BigInt res;
    res.limbs_.assign(value.Limbs().begin(), value.Limbs().end());
    res.Normalize();
    return res;
  }

  template <std::size_t Bits>
    requires(Bits >= 64 && Bits % 64 == 0)
  BigUInt<Bits> ToBigUInt() const {
    if (negative_) {
      throw std::overflow_error(
          "Cannot convert negative BigInt to BigUInt<Bits>");
    }
    if (limbs_.size() > BigUInt<Bits>::kLimbs) {
      throw std::overflow_error("BigInt value overflows BigUInt<Bits>");
    }
    std::array<uint64_t, BigUInt<Bits>::kLimbs> out{};
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      out[i] = limbs_[i];
    }
    return BigUInt<Bits>(out);
  }

  bool IsZero() const noexcept { return limbs_.empty(); }

  bool IsOne() const noexcept {
    return !negative_ && limbs_.size() == 1 && limbs_[0] == 1;
  }

  bool IsNegative() const noexcept { return negative_; }

  bool IsPositive() const noexcept { return !negative_ && !limbs_.empty(); }

  bool IsEven() const noexcept {
    return limbs_.empty() || (limbs_[0] & 1ULL) == 0;
  }

  bool IsOdd() const noexcept {
    return !limbs_.empty() && (limbs_[0] & 1ULL) != 0;
  }

  int Sign() const noexcept {
    if (limbs_.empty()) {
      return 0;
    }
    return negative_ ? -1 : 1;
  }

  std::size_t BitLength() const noexcept {
    if (limbs_.empty()) {
      return 0;
    }
    return (limbs_.size() - 1) * 64 +
           static_cast<std::size_t>(std::bit_width(limbs_.back()));
  }

  bool GetBit(std::size_t index) const noexcept {
    const std::size_t word = index / 64;
    if (word >= limbs_.size()) {
      return false;
    }
    return ((limbs_[word] >> (index % 64)) & 1ULL) != 0;
  }

  std::size_t LimbCount() const noexcept { return limbs_.size(); }

  uint64_t Limb(std::size_t index) const {
    if (index >= limbs_.size()) {
      throw std::out_of_range("Limb index out of range in BigInt::Limb");
    }
    return limbs_[index];
  }

  int64_t ToInt64() const {
    if (limbs_.empty()) {
      return 0;
    }
    if (limbs_.size() > 1) {
      throw std::overflow_error("BigInt out of int64_t range");
    }
    const uint64_t v = limbs_[0];
    if (!negative_) {
      if (v > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
        throw std::overflow_error("BigInt out of int64_t range");
      }
      return static_cast<int64_t>(v);
    }
    constexpr uint64_t kMinMag = uint64_t{1} << 63;
    if (v == kMinMag) {
      return std::numeric_limits<int64_t>::min();
    }
    if (v > kMinMag) {
      throw std::overflow_error("BigInt out of int64_t range");
    }
    return -static_cast<int64_t>(v);
  }

  uint64_t ToUint64() const {
    if (negative_ || limbs_.size() > 1) {
      throw std::overflow_error("BigInt out of uint64_t range");
    }
    return limbs_.empty() ? 0ULL : limbs_[0];
  }

  std::string ToDecimal() const {
    return detail::LimbsToDecimal(limbs_, negative_);
  }

  std::string ToHex(bool prefix = false) const {
    return detail::LimbsToHex(limbs_, negative_, prefix);
  }

  std::string ToString() const { return ToDecimal(); }

  BigInt Abs() const {
    BigInt res = *this;
    res.negative_ = false;
    return res;
  }

  BigInt operator+() const { return *this; }

  BigInt operator-() const {
    BigInt res = *this;
    if (!res.limbs_.empty()) {
      res.negative_ = !res.negative_;
    }
    return res;
  }

  BigInt& operator+=(const BigInt& rhs) {
    if (rhs.IsZero()) {
      return *this;
    }
    if (IsZero()) {
      *this = rhs;
      return *this;
    }
    if (negative_ == rhs.negative_) {
      limbs_ = detail::AddLimbs(limbs_, rhs.limbs_);
    } else {
      const auto cmp = detail::CmpLimbs(limbs_, rhs.limbs_);
      if (cmp == std::strong_ordering::equal) {
        limbs_.clear();
        negative_ = false;
      } else if (cmp == std::strong_ordering::greater) {
        limbs_ = detail::SubLimbs(limbs_, rhs.limbs_);
      } else {
        limbs_ = detail::SubLimbs(rhs.limbs_, limbs_);
        negative_ = rhs.negative_;
      }
    }
    Normalize();
    return *this;
  }

  BigInt& operator-=(const BigInt& rhs) {
    if (rhs.IsZero()) {
      return *this;
    }
    if (IsZero()) {
      *this = -rhs;
      return *this;
    }
    if (negative_ != rhs.negative_) {
      limbs_ = detail::AddLimbs(limbs_, rhs.limbs_);
    } else {
      const auto cmp = detail::CmpLimbs(limbs_, rhs.limbs_);
      if (cmp == std::strong_ordering::equal) {
        limbs_.clear();
        negative_ = false;
      } else if (cmp == std::strong_ordering::greater) {
        limbs_ = detail::SubLimbs(limbs_, rhs.limbs_);
      } else {
        limbs_ = detail::SubLimbs(rhs.limbs_, limbs_);
        negative_ = !negative_;
      }
    }
    Normalize();
    return *this;
  }

  BigInt& operator*=(const BigInt& rhs) {
    if (IsZero() || rhs.IsZero()) {
      limbs_.clear();
      negative_ = false;
      return *this;
    }
    limbs_ = detail::MulLimbsKaratsuba(limbs_, rhs.limbs_);
    negative_ = (negative_ != rhs.negative_);
    Normalize();
    return *this;
  }

  // Truncated division matching C++ '/' and '%': a == q * b + r, |r| < |b|,
  // sign(r) == sign(a) when r != 0.
  static std::pair<BigInt, BigInt> DivRem(const BigInt& a, const BigInt& b) {
    if (b.IsZero()) {
      throw std::domain_error("Division by zero in BigInt::DivRem");
    }
    auto [q_limbs, r_limbs] = detail::DivRemLimbs(a.limbs_, b.limbs_);
    BigInt q = FromLimbs(std::move(q_limbs), a.negative_ != b.negative_);
    BigInt r = FromLimbs(std::move(r_limbs), a.negative_);
    return {std::move(q), std::move(r)};
  }

  // Euclidean division: a == q * b + r with 0 <= r < |b|.
  static std::pair<BigInt, BigInt> EuclideanDivRem(const BigInt& a,
                                                   const BigInt& b) {
    auto [q, r] = DivRem(a, b);
    if (r.IsNegative()) {
      if (b.IsPositive()) {
        q -= BigInt(1);
        r += b;
      } else {
        q += BigInt(1);
        r -= b;
      }
    }
    return {std::move(q), std::move(r)};
  }

  BigInt& operator/=(const BigInt& rhs) {
    *this = DivRem(*this, rhs).first;
    return *this;
  }

  BigInt& operator%=(const BigInt& rhs) {
    *this = DivRem(*this, rhs).second;
    return *this;
  }

  BigInt& operator<<=(std::size_t shift) {
    if (IsZero() || shift == 0) {
      return *this;
    }
    const std::size_t word_shift = shift / 64;
    const unsigned bit_shift = static_cast<unsigned>(shift % 64);
    const std::size_t old_size = limbs_.size();
    std::vector<uint64_t> out(old_size + word_shift + (bit_shift != 0 ? 1 : 0),
                              0);
    if (bit_shift == 0) {
      for (std::size_t i = 0; i < old_size; ++i) {
        out[i + word_shift] = limbs_[i];
      }
    } else {
      const unsigned rshift = 64U - bit_shift;
      uint64_t carry = 0;
      for (std::size_t i = 0; i < old_size; ++i) {
        out[i + word_shift] = (limbs_[i] << bit_shift) | carry;
        carry = limbs_[i] >> rshift;
      }
      out[old_size + word_shift] = carry;
    }
    limbs_ = std::move(out);
    Normalize();
    return *this;
  }

  // Arithmetic right shift: floor(x / 2^shift).
  BigInt& operator>>=(std::size_t shift) {
    if (IsZero() || shift == 0) {
      return *this;
    }
    const std::size_t word_shift = shift / 64;
    const unsigned bit_shift = static_cast<unsigned>(shift % 64);

    bool had_shifted_out_one = false;
    if (negative_) {
      const std::size_t full_words = std::min(word_shift, limbs_.size());
      for (std::size_t i = 0; i < full_words; ++i) {
        if (limbs_[i] != 0) {
          had_shifted_out_one = true;
          break;
        }
      }
      if (!had_shifted_out_one && word_shift < limbs_.size() &&
          bit_shift != 0) {
        const uint64_t mask = (uint64_t{1} << bit_shift) - 1ULL;
        if ((limbs_[word_shift] & mask) != 0) {
          had_shifted_out_one = true;
        }
      }
    }

    if (word_shift >= limbs_.size()) {
      limbs_.clear();
      if (negative_ && had_shifted_out_one) {
        limbs_.push_back(1);
      } else {
        negative_ = false;
      }
      return *this;
    }

    const std::size_t new_size = limbs_.size() - word_shift;
    std::vector<uint64_t> out(new_size, 0);
    if (bit_shift == 0) {
      for (std::size_t i = 0; i < new_size; ++i) {
        out[i] = limbs_[i + word_shift];
      }
    } else {
      const unsigned lshift = 64U - bit_shift;
      for (std::size_t i = 0; i + 1 < new_size; ++i) {
        out[i] = (limbs_[i + word_shift] >> bit_shift) |
                 (limbs_[i + word_shift + 1] << lshift);
      }
      out[new_size - 1] = limbs_.back() >> bit_shift;
    }

    limbs_ = std::move(out);
    if (negative_ && had_shifted_out_one) {
      const uint64_t one = 1;
      limbs_ = detail::AddLimbs(limbs_, std::span<const uint64_t>(&one, 1));
    }
    Normalize();
    return *this;
  }

  BigInt& operator&=(const BigInt& rhs) {
    *this = BitwiseOp(*this, rhs, 0);
    return *this;
  }

  BigInt& operator|=(const BigInt& rhs) {
    *this = BitwiseOp(*this, rhs, 1);
    return *this;
  }

  BigInt& operator^=(const BigInt& rhs) {
    *this = BitwiseOp(*this, rhs, 2);
    return *this;
  }

  BigInt operator~() const { return -(*this + BigInt(1)); }

  friend BigInt operator+(BigInt lhs, const BigInt& rhs) {
    lhs += rhs;
    return lhs;
  }

  friend BigInt operator-(BigInt lhs, const BigInt& rhs) {
    lhs -= rhs;
    return lhs;
  }

  friend BigInt operator*(BigInt lhs, const BigInt& rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend BigInt operator/(const BigInt& lhs, const BigInt& rhs) {
    return DivRem(lhs, rhs).first;
  }

  friend BigInt operator%(const BigInt& lhs, const BigInt& rhs) {
    return DivRem(lhs, rhs).second;
  }

  friend BigInt operator<<(BigInt lhs, std::size_t shift) {
    lhs <<= shift;
    return lhs;
  }

  friend BigInt operator>>(BigInt lhs, std::size_t shift) {
    lhs >>= shift;
    return lhs;
  }

  friend BigInt operator&(const BigInt& lhs, const BigInt& rhs) {
    return BitwiseOp(lhs, rhs, 0);
  }

  friend BigInt operator|(const BigInt& lhs, const BigInt& rhs) {
    return BitwiseOp(lhs, rhs, 1);
  }

  friend BigInt operator^(const BigInt& lhs, const BigInt& rhs) {
    return BitwiseOp(lhs, rhs, 2);
  }

  bool operator==(const BigInt& rhs) const noexcept {
    return negative_ == rhs.negative_ && limbs_ == rhs.limbs_;
  }

  std::strong_ordering operator<=>(const BigInt& rhs) const noexcept {
    const int s1 = Sign();
    const int s2 = rhs.Sign();
    if (s1 != s2) {
      return s1 <=> s2;
    }
    if (s1 == 0) {
      return std::strong_ordering::equal;
    }
    if (s1 > 0) {
      return detail::CmpLimbs(limbs_, rhs.limbs_);
    }
    return detail::CmpLimbs(rhs.limbs_, limbs_);
  }

  friend std::ostream& operator<<(std::ostream& os, const BigInt& value) {
    return os << value.ToDecimal();
  }

  BigInt Pow(uint64_t exp) const {
    BigInt result(1);
    BigInt base = *this;
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

  std::size_t CountTrailingZeroBits() const noexcept {
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      if (limbs_[i] != 0) {
        return i * 64 + static_cast<std::size_t>(std::countr_zero(limbs_[i]));
      }
    }
    return 0;
  }

  static BigInt PowMod(const BigInt& base, const BigInt& exp,
                       const BigInt& mod);
  static BigInt Gcd(const BigInt& a, const BigInt& b);
  static std::tuple<BigInt, BigInt, BigInt> ExtendedGcd(const BigInt& a,
                                                        const BigInt& b);
  static BigInt ModInverse(const BigInt& a, const BigInt& mod);
  static BigInt ISqrt(const BigInt& n);
  static BigInt Factorial(uint64_t n);

  // Access underlying magnitude limbs for hashing and internal helpers.
  std::span<const uint64_t> RawLimbs() const noexcept { return limbs_; }

 private:
  static BigInt FromLimbs(std::vector<uint64_t> limbs, bool negative) {
    BigInt res;
    res.limbs_ = std::move(limbs);
    res.negative_ = negative;
    res.Normalize();
    return res;
  }

  void Normalize() noexcept {
    detail::TrimLimbs(limbs_);
    if (limbs_.empty()) {
      negative_ = false;
    }
  }

  // Convert to two's complement representation of length `len` limbs.
  static std::vector<uint64_t> ToTwosComplement(const BigInt& x,
                                                std::size_t len) {
    std::vector<uint64_t> out(len, 0);
    if (!x.negative_) {
      for (std::size_t i = 0; i < x.limbs_.size(); ++i) {
        out[i] = x.limbs_[i];
      }
      return out;
    }
    uint64_t carry = 1;
    for (std::size_t i = 0; i < len; ++i) {
      const uint64_t v = (i < x.limbs_.size()) ? x.limbs_[i] : 0ULL;
      const detail::u128 sum = static_cast<detail::u128>(~v) + carry;
      out[i] = static_cast<uint64_t>(sum);
      carry = static_cast<uint64_t>(sum >> 64);
    }
    return out;
  }

  static BigInt BitwiseOp(const BigInt& a, const BigInt& b, int op) {
    // Fast path for non-negative operands.
    if (!a.negative_ && !b.negative_) {
      const std::size_t max_len = std::max(a.limbs_.size(), b.limbs_.size());
      std::vector<uint64_t> out(max_len, 0);
      for (std::size_t i = 0; i < max_len; ++i) {
        const uint64_t av = (i < a.limbs_.size()) ? a.limbs_[i] : 0ULL;
        const uint64_t bv = (i < b.limbs_.size()) ? b.limbs_[i] : 0ULL;
        if (op == 0) {
          out[i] = av & bv;
        } else if (op == 1) {
          out[i] = av | bv;
        } else {
          out[i] = av ^ bv;
        }
      }
      return FromLimbs(std::move(out), false);
    }

    const std::size_t len = std::max(a.limbs_.size(), b.limbs_.size()) + 1;
    const auto ta = ToTwosComplement(a, len);
    const auto tb = ToTwosComplement(b, len);
    std::vector<uint64_t> tr(len, 0);
    bool res_neg = false;
    if (op == 0) {
      res_neg = a.negative_ && b.negative_;
      for (std::size_t i = 0; i < len; ++i) {
        tr[i] = ta[i] & tb[i];
      }
    } else if (op == 1) {
      res_neg = a.negative_ || b.negative_;
      for (std::size_t i = 0; i < len; ++i) {
        tr[i] = ta[i] | tb[i];
      }
    } else {
      res_neg = a.negative_ != b.negative_;
      for (std::size_t i = 0; i < len; ++i) {
        tr[i] = ta[i] ^ tb[i];
      }
    }

    if (!res_neg) {
      return FromLimbs(std::move(tr), false);
    }
    // Convert negative two's complement back to magnitude: (~tr) + 1.
    uint64_t carry = 1;
    for (std::size_t i = 0; i < len; ++i) {
      const detail::u128 sum = static_cast<detail::u128>(~tr[i]) + carry;
      tr[i] = static_cast<uint64_t>(sum);
      carry = static_cast<uint64_t>(sum >> 64);
    }
    if (carry != 0) {
      tr.push_back(carry);
    }
    return FromLimbs(std::move(tr), true);
  }

  bool negative_ = false;
  std::vector<uint64_t> limbs_;
};

// Free functions on BigInt in namespace xg

inline std::pair<BigInt, BigInt> DivRem(const BigInt& a, const BigInt& b) {
  return BigInt::DivRem(a, b);
}

inline std::pair<BigInt, BigInt> EuclideanDivRem(const BigInt& a,
                                                 const BigInt& b) {
  return BigInt::EuclideanDivRem(a, b);
}

inline BigInt Gcd(const BigInt& a, const BigInt& b) {
  BigInt x = a.Abs();
  BigInt y = b.Abs();
  while (!y.IsZero()) {
    if (x.LimbCount() <= 1 && y.LimbCount() <= 1) {
      return BigInt::FromUint64(detail::GcdU64(x.ToUint64(), y.ToUint64()));
    }
    BigInt r = x % y;
    x = std::move(y);
    y = std::move(r);
  }
  return x;
}

inline std::tuple<BigInt, BigInt, BigInt> ExtendedGcd(const BigInt& a,
                                                      const BigInt& b) {
  BigInt old_r = a;
  BigInt r = b;
  BigInt old_x(1);
  BigInt x(0);
  BigInt old_y(0);
  BigInt y(1);

  while (!r.IsZero()) {
    auto [q, rem] = BigInt::DivRem(old_r, r);
    old_r = std::move(r);
    r = std::move(rem);

    BigInt next_x = old_x - q * x;
    old_x = std::move(x);
    x = std::move(next_x);

    BigInt next_y = old_y - q * y;
    old_y = std::move(y);
    y = std::move(next_y);
  }

  if (old_r.IsNegative()) {
    old_r = -old_r;
    old_x = -old_x;
    old_y = -old_y;
  }
  return {std::move(old_r), std::move(old_x), std::move(old_y)};
}

inline BigInt ModInverse(const BigInt& a, const BigInt& mod) {
  if (mod <= BigInt(1)) {
    throw std::domain_error("Modulus must be greater than 1 in ModInverse");
  }
  const BigInt a_mod = BigInt::EuclideanDivRem(a, mod).second;
  if (a_mod.IsZero()) {
    throw std::domain_error("Modular inverse does not exist (gcd != 1)");
  }
  auto [g, x, unused_y] = ExtendedGcd(a_mod, mod);
  (void)unused_y;
  if (!g.IsOne()) {
    throw std::domain_error("Modular inverse does not exist (gcd != 1)");
  }
  return BigInt::EuclideanDivRem(x, mod).second;
}

inline BigInt PowMod(const BigInt& base, const BigInt& exp, const BigInt& mod) {
  if (mod <= BigInt(0)) {
    throw std::domain_error("Modulus must be positive in PowMod");
  }
  if (mod.IsOne()) {
    return BigInt(0);
  }
  if (exp.IsNegative()) {
    const BigInt inv = ModInverse(base, mod);
    return PowMod(inv, -exp, mod);
  }
  BigInt b = BigInt::EuclideanDivRem(base, mod).second;
  BigInt result(1);
  const std::size_t bits = exp.BitLength();
  for (std::size_t i = bits; i-- > 0;) {
    result = (result * result) % mod;
    if (exp.GetBit(i)) {
      result = (result * b) % mod;
    }
  }
  return result;
}

inline BigInt ISqrt(const BigInt& n) {
  if (n.IsNegative()) {
    throw std::domain_error("ISqrt requires non-negative input");
  }
  if (n.IsZero() || n.IsOne()) {
    return n;
  }
  if (n.LimbCount() == 1) {
    return BigInt::FromUint64(detail::ISqrtU64(n.ToUint64()));
  }
  BigInt x = BigInt(1) << ((n.BitLength() + 1) / 2);
  while (true) {
    BigInt y = (x + n / x) >> 1;
    if (y >= x) {
      return x;
    }
    x = std::move(y);
  }
}

namespace detail {

inline BigInt ProductRange(uint64_t lo, uint64_t hi) {
  if (lo > hi) {
    return BigInt(1);
  }
  if (lo == hi) {
    return BigInt::FromUint64(lo);
  }
  if (hi - lo <= 16) {
    BigInt acc = BigInt::FromUint64(lo);
    for (uint64_t v = lo + 1; v <= hi; ++v) {
      acc *= BigInt::FromUint64(v);
    }
    return acc;
  }
  const uint64_t mid = lo + (hi - lo) / 2;
  return ProductRange(lo, mid) * ProductRange(mid + 1, hi);
}

}  // namespace detail

inline BigInt Factorial(uint64_t n) {
  if (n <= 1) {
    return BigInt(1);
  }
  return detail::ProductRange(2, n);
}

// Forward static member methods on BigInt to the free functions.
inline BigInt BigInt::PowMod(const BigInt& base, const BigInt& exp,
                             const BigInt& mod) {
  return ::xg::PowMod(base, exp, mod);
}

inline BigInt BigInt::Gcd(const BigInt& a, const BigInt& b) {
  return ::xg::Gcd(a, b);
}

inline std::tuple<BigInt, BigInt, BigInt> BigInt::ExtendedGcd(const BigInt& a,
                                                              const BigInt& b) {
  return ::xg::ExtendedGcd(a, b);
}

inline BigInt BigInt::ModInverse(const BigInt& a, const BigInt& mod) {
  return ::xg::ModInverse(a, mod);
}

inline BigInt BigInt::ISqrt(const BigInt& n) { return ::xg::ISqrt(n); }

inline BigInt BigInt::Factorial(uint64_t n) { return ::xg::Factorial(n); }

// ============================================================================
// Number theory: JacobiSymbol, IsProbablePrime (Baillie-PSW + MR),
// FactorInteger
// ============================================================================

inline int JacobiSymbol(const BigInt& a_in, const BigInt& n_in) {
  if (!n_in.IsPositive() || n_in.IsEven()) {
    throw std::domain_error("JacobiSymbol requires a positive odd modulus n");
  }
  if (n_in.IsOne()) {
    return 1;
  }

  BigInt a = BigInt::EuclideanDivRem(a_in, n_in).second;
  BigInt n = n_in;
  int result = 1;

  while (!a.IsZero()) {
    const std::size_t tz = a.CountTrailingZeroBits();
    if (tz > 0) {
      a >>= tz;
      if ((tz & 1U) != 0) {
        const uint64_t n_mod_8 = n.Limb(0) & 7ULL;
        if (n_mod_8 == 3 || n_mod_8 == 5) {
          result = -result;
        }
      }
    }
    std::swap(a, n);
    if ((a.Limb(0) & 3ULL) == 3 && (n.Limb(0) & 3ULL) == 3) {
      result = -result;
    }
    a %= n;
  }
  return n.IsOne() ? result : 0;
}

namespace detail {

inline bool MillerRabinWitness(const BigInt& n, const BigInt& n_minus_1,
                               const BigInt& d, std::size_t s,
                               const BigInt& a) {
  BigInt a_mod = a % n;
  if (a_mod.IsZero()) {
    return true;
  }
  BigInt x = PowMod(a_mod, d, n);
  if (x.IsOne() || x == n_minus_1) {
    return true;
  }
  for (std::size_t r = 1; r < s; ++r) {
    x = (x * x) % n;
    if (x == n_minus_1) {
      return true;
    }
    if (x.IsOne()) {
      return false;
    }
  }
  return false;
}

// Compute x / 2 mod n where n is odd and 0 <= x < n.
inline BigInt Div2ModOdd(BigInt x, const BigInt& n) {
  if (x.IsOdd()) {
    x += n;
  }
  x >>= 1;
  return x;
}

// Strong Lucas Probable Prime test using Selfridge's Method A.
// Assumes n > 2, n is odd, and n is not a perfect square.
inline bool StrongLucasProbablePrime(const BigInt& n) {
  int64_t abs_d = 5;
  int64_t sign = 1;
  int64_t d_val = 5;

  while (true) {
    d_val = sign * abs_d;
    const int j = JacobiSymbol(BigInt(d_val), n);
    if (j == -1) {
      break;
    }
    if (j == 0) {
      // gcd(|D|, n) > 1. Since n > |D| (or equals |D|), check divisibility.
      if (BigInt(abs_d) < n) {
        return false;
      }
    }
    abs_d += 2;
    sign = -sign;
  }

  // P = 1, Q = (1 - D) / 4.
  const int64_t q_val = (1 - d_val) / 4;
  const BigInt d_mod = BigInt::EuclideanDivRem(BigInt(d_val), n).second;
  const BigInt q_mod = BigInt::EuclideanDivRem(BigInt(q_val), n).second;

  // Check gcd(|Q|, n) == 1.
  const BigInt gq = Gcd(BigInt(q_val).Abs(), n);
  if (gq > BigInt(1) && gq < n) {
    return false;
  }

  // Write n + 1 = d_l * 2^s_l with d_l odd.
  BigInt d_l = n + BigInt(1);
  const std::size_t s_l = d_l.CountTrailingZeroBits();
  d_l >>= s_l;

  // Binary Lucas chain for (U_k, V_k, Q^k) mod n with P = 1.
  BigInt u(1);
  BigInt v(1);
  BigInt qk = q_mod;

  const std::size_t bits = d_l.BitLength();
  for (std::size_t i = bits - 1; i-- > 0;) {
    // Double: k -> 2k.
    BigInt u2 = (u * v) % n;
    BigInt v2 = BigInt::EuclideanDivRem(v * v - (qk << 1), n).second;
    BigInt qk2 = (qk * qk) % n;

    u = std::move(u2);
    v = std::move(v2);
    qk = std::move(qk2);

    if (d_l.GetBit(i)) {
      // Increment: m -> m + 1 (P = 1).
      BigInt u_sum = u + v;
      if (u_sum >= n) {
        u_sum -= n;
      }
      BigInt u_next = Div2ModOdd(std::move(u_sum), n);

      BigInt v_sum = (d_mod * u + v) % n;
      BigInt v_next = Div2ModOdd(std::move(v_sum), n);

      BigInt qk_next = (qk * q_mod) % n;

      u = std::move(u_next);
      v = std::move(v_next);
      qk = std::move(qk_next);
    }
  }

  if (u.IsZero() || v.IsZero()) {
    return true;
  }

  for (std::size_t r = 1; r < s_l; ++r) {
    v = BigInt::EuclideanDivRem(v * v - (qk << 1), n).second;
    if (v.IsZero()) {
      return true;
    }
    qk = (qk * qk) % n;
  }
  return false;
}

}  // namespace detail

inline bool IsProbablePrime(const BigInt& n, int mr_rounds = 25) {
  if (n < BigInt(2)) {
    return false;
  }

  // Small prime trial division.
  for (uint64_t p : detail::kSmallPrimesTable) {
    const BigInt bp = BigInt::FromUint64(p);
    if (n == bp) {
      return true;
    }
    if ((n % bp).IsZero()) {
      return false;
    }
  }

  // If n < 257^2 and survived trial division up to 251, n is prime.
  if (n < BigInt(257ULL * 257ULL)) {
    return true;
  }

  // Deterministic 64-bit Miller-Rabin.
  if (n.LimbCount() == 1) {
    return detail::IsPrimeU64(n.ToUint64());
  }

  // Baillie-PSW:
  // 1. Base-2 strong Miller-Rabin test.
  const BigInt n_minus_1 = n - BigInt(1);
  const std::size_t s = n_minus_1.CountTrailingZeroBits();
  const BigInt d = n_minus_1 >> s;

  if (!detail::MillerRabinWitness(n, n_minus_1, d, s, BigInt(2))) {
    return false;
  }

  // 2. Reject perfect squares (required before Selfridge Method A).
  const BigInt root = ISqrt(n);
  if (root * root == n) {
    return false;
  }

  // 3. Strong Lucas probable prime test.
  if (!detail::StrongLucasProbablePrime(n)) {
    return false;
  }

  // 4. Additional deterministic Miller-Rabin rounds.
  if (mr_rounds <= 1) {
    return true;
  }
  int remaining = mr_rounds - 1;
  // First use small prime witnesses (skipping 2, already tested).
  for (std::size_t i = 1; i < detail::kSmallPrimesTable.size() && remaining > 0;
       ++i, --remaining) {
    if (!detail::MillerRabinWitness(
            n, n_minus_1, d, s,
            BigInt::FromUint64(detail::kSmallPrimesTable[i]))) {
      return false;
    }
  }

  // If more rounds were requested, synthesize deterministic witnesses in [2,
  // n-2].
  uint64_t rng_state = 0x243f6a8885a308d3ULL ^ n.Limb(0);
  const BigInt n_minus_3 = n - BigInt(3);
  while (remaining > 0) {
    BigInt witness(0);
    for (std::size_t w = 0; w < n.LimbCount(); ++w) {
      witness =
          (witness << 64) + BigInt::FromUint64(detail::SplitMix64(rng_state));
    }
    witness = (witness % n_minus_3) + BigInt(2);
    if (!detail::MillerRabinWitness(n, n_minus_1, d, s, witness)) {
      return false;
    }
    --remaining;
  }
  return true;
}

namespace detail {

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline std::size_t CountTrailingZeroBitsUInt(const BigUInt<Bits>& x) noexcept {
  for (std::size_t i = 0; i < BigUInt<Bits>::kLimbs; ++i) {
    const uint64_t w = x.Limbs()[i];
    if (w != 0) {
      return i * 64 + static_cast<std::size_t>(std::countr_zero(w));
    }
  }
  return Bits;
}

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigUInt<Bits> GcdUInt(BigUInt<Bits> u, BigUInt<Bits> v) noexcept {
  if (u.IsZero()) {
    return v;
  }
  if (v.IsZero()) {
    return u;
  }
  const std::size_t shift =
      std::min(CountTrailingZeroBitsUInt(u), CountTrailingZeroBitsUInt(v));
  u >>= CountTrailingZeroBitsUInt(u);
  do {
    v >>= CountTrailingZeroBitsUInt(v);
    if (u > v) {
      std::swap(u, v);
    }
    v -= u;
  } while (!v.IsZero());
  return u << shift;
}

template <std::size_t Bits>
  requires(Bits >= 64 && Bits % 64 == 0)
inline BigInt PollardsRhoBrentMont(const BigInt& n_big) {
  const BigUInt<Bits> n = n_big.ToBigUInt<Bits>();
  MontgomeryCtx<Bits> ctx(n);
  constexpr uint64_t kBatch = 512;

  for (uint64_t attempt = 1; attempt <= 64; ++attempt) {
    const BigUInt<Bits> c = ctx.ToMont(BigUInt<Bits>(attempt));
    BigUInt<Bits> y = ctx.ToMont(BigUInt<Bits>(attempt + 1));
    BigUInt<Bits> q = ctx.RMod();
    BigUInt<Bits> x = y;
    BigUInt<Bits> ys = y;
    uint64_t r = 1;
    BigUInt<Bits> g = BigUInt<Bits>::One();

    auto step = [&ctx, &c](const BigUInt<Bits>& v) noexcept -> BigUInt<Bits> {
      if constexpr (Bits == 128) {
        return ctx.Add(ctx.Mul(v, v), c);
      } else {
        return ctx.Add(ctx.Sqr(v), c);
      }
    };

    while (g.IsOne()) {
      x = y;
      for (uint64_t i = 0; i < r; ++i) {
        y = step(y);
      }
      uint64_t k = 0;
      while (k < r && g.IsOne()) {
        ys = y;
        const uint64_t limit = std::min(kBatch, r - k);
        for (uint64_t i = 0; i < limit; ++i) {
          y = step(y);
          const BigUInt<Bits> diff = (x >= y) ? (x - y) : (y - x);
          if (diff.IsZero()) {
            g = n;
            break;
          }
          const BigUInt<Bits> next_q = ctx.Mul(q, diff);
          if (next_q.IsZero()) {
            g = GcdUInt(diff, n);
            break;
          }
          q = next_q;
        }
        if (g.IsOne()) {
          g = GcdUInt(q, n);
        }
        k += limit;
      }
      r <<= 1;
    }

    if (g == n) {
      for (uint64_t cnt = 0; cnt <= kBatch + 2; ++cnt) {
        ys = step(ys);
        const BigUInt<Bits> diff = (x >= ys) ? (x - ys) : (ys - x);
        g = GcdUInt(diff, n);
        if (g > BigUInt<Bits>::One()) {
          break;
        }
      }
    }

    if (g > BigUInt<Bits>::One() && g < n) {
      return BigInt::FromBigUInt(g);
    }
  }
  return n_big;
}

inline BigInt PollardsRhoBrentBig(const BigInt& n) {
  constexpr uint64_t kBatch = 128;
  for (uint64_t attempt = 1; attempt <= 64; ++attempt) {
    const BigInt c = BigInt::FromUint64(attempt);
    BigInt y = BigInt::FromUint64(attempt + 1) % n;
    BigInt q(1);
    BigInt x = y;
    BigInt ys = y;
    uint64_t r = 1;
    BigInt g(1);

    auto step = [&n, &c](const BigInt& v) -> BigInt { return (v * v + c) % n; };

    while (g.IsOne()) {
      x = y;
      for (uint64_t i = 0; i < r; ++i) {
        y = step(y);
      }
      uint64_t k = 0;
      while (k < r && g.IsOne()) {
        ys = y;
        const uint64_t limit = std::min(kBatch, r - k);
        for (uint64_t i = 0; i < limit; ++i) {
          y = step(y);
          const BigInt diff = (x >= y) ? (x - y) : (y - x);
          if (diff.IsZero()) {
            g = n;
            break;
          }
          const BigInt next_q = (q * diff) % n;
          if (next_q.IsZero()) {
            g = Gcd(diff, n);
            break;
          }
          q = next_q;
        }
        if (g.IsOne()) {
          g = Gcd(q, n);
        }
        k += limit;
      }
      r <<= 1;
    }

    if (g == n) {
      for (uint64_t cnt = 0; cnt <= kBatch + 2; ++cnt) {
        ys = step(ys);
        const BigInt diff = (x >= ys) ? (x - ys) : (ys - x);
        g = Gcd(diff, n);
        if (g > BigInt(1)) {
          break;
        }
      }
    }

    if (g > BigInt(1) && g < n) {
      return g;
    }
  }
  return n;
}

inline void FactorRecursive(const BigInt& n, std::vector<BigInt>& factors) {
  if (n <= BigInt(1)) {
    return;
  }
  if (IsProbablePrime(n)) {
    factors.push_back(n);
    return;
  }

  // Handle perfect squares immediately so Pollard's Rho does not stall on p^2.
  const BigInt root = ISqrt(n);
  if (root * root == n) {
    FactorRecursive(root, factors);
    FactorRecursive(root, factors);
    return;
  }

  BigInt d = n;
  const std::size_t bits = n.BitLength();
  if (bits <= 64) {
    d = BigInt::FromUint64(PollardsRhoBrentU64(n.ToUint64()));
  } else if (bits <= 128) {
    d = PollardsRhoBrentMont<128>(n);
  } else if (bits <= 256) {
    d = PollardsRhoBrentMont<256>(n);
  } else {
    d = PollardsRhoBrentBig(n);
  }

  if (d == n || d <= BigInt(1)) {
    // Should not occur for composite numbers; record remainder to terminate.
    factors.push_back(n);
    return;
  }
  FactorRecursive(d, factors);
  FactorRecursive(n / d, factors);
}

}  // namespace detail

inline std::vector<std::pair<BigInt, uint64_t>> FactorInteger(const BigInt& n) {
  if (!n.IsPositive()) {
    throw std::domain_error("FactorInteger requires a positive integer");
  }
  std::vector<std::pair<BigInt, uint64_t>> result;
  if (n.IsOne()) {
    return result;
  }

  BigInt rem = n;

  // Trial division by 2 and odd primes up to 1000.
  const std::size_t tz = rem.CountTrailingZeroBits();
  if (tz > 0) {
    result.emplace_back(BigInt(2), static_cast<uint64_t>(tz));
    rem >>= tz;
  }

  for (uint64_t d = 3; d <= 997; d += 2) {
    if (rem.IsOne()) {
      break;
    }
    if (!detail::IsPrimeU64(d)) {
      continue;
    }
    const BigInt bd = BigInt::FromUint64(d);
    if (bd * bd > rem) {
      break;
    }
    uint64_t exp = 0;
    while (true) {
      auto [q, r] = BigInt::DivRem(rem, bd);
      if (!r.IsZero()) {
        break;
      }
      rem = std::move(q);
      ++exp;
    }
    if (exp > 0) {
      result.emplace_back(bd, exp);
    }
  }

  if (rem > BigInt(1)) {
    std::vector<BigInt> large_factors;
    detail::FactorRecursive(rem, large_factors);
    std::sort(large_factors.begin(), large_factors.end());
    for (const BigInt& p : large_factors) {
      if (!result.empty() && result.back().first == p) {
        ++result.back().second;
      } else {
        result.emplace_back(p, 1);
      }
    }
  }

  return result;
}

}  // namespace xg

namespace std {

template <>
struct hash<xg::BigInt> {
  std::size_t operator()(const xg::BigInt& value) const noexcept {
    uint64_t h = value.IsNegative() ? 0x9e3779b97f4a7c15ULL : 0ULL;
    for (uint64_t limb : value.RawLimbs()) {
      uint64_t z = limb + 0x9e3779b97f4a7c15ULL;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      z ^= (z >> 31);
      h ^= z + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return static_cast<std::size_t>(h);
  }
};

}  // namespace std

#endif  // XGALOIS_UTILS_BIGINT_HPP_
