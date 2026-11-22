#ifndef XGALOIS_CHANNEL_RNG_HPP_
#define XGALOIS_CHANNEL_RNG_HPP_

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <utility>

namespace xg {
namespace channel {

// Portable 128-bit unsigned integer for strict ISO C++20 (-Wpedantic) PCG64
// arithmetic without compiler-specific __uint128_t extensions.
struct Uint128 {
  uint64_t hi = 0;
  uint64_t lo = 0;

  constexpr bool operator==(const Uint128& other) const noexcept = default;
};

namespace detail {

inline constexpr Uint128 Add128(Uint128 a, Uint128 b) noexcept {
  const uint64_t lo = a.lo + b.lo;
  const uint64_t carry = (lo < a.lo) ? uint64_t{1} : uint64_t{0};
  const uint64_t hi = a.hi + b.hi + carry;
  return Uint128{hi, lo};
}

inline constexpr Uint128 Mul64x64(uint64_t a, uint64_t b) noexcept {
  const uint64_t a0 = a & 0xffffffffULL;
  const uint64_t a1 = a >> 32;
  const uint64_t b0 = b & 0xffffffffULL;
  const uint64_t b1 = b >> 32;

  const uint64_t p00 = a0 * b0;
  const uint64_t p01 = a0 * b1;
  const uint64_t p10 = a1 * b0;
  const uint64_t p11 = a1 * b1;

  const uint64_t mid =
      (p00 >> 32) + (p01 & 0xffffffffULL) + (p10 & 0xffffffffULL);
  const uint64_t lo = (p00 & 0xffffffffULL) | (mid << 32);
  const uint64_t hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
  return Uint128{hi, lo};
}

inline constexpr Uint128 Mul128(Uint128 a, Uint128 b) noexcept {
  Uint128 prod = Mul64x64(a.lo, b.lo);
  prod.hi += a.lo * b.hi + a.hi * b.lo;
  return prod;
}

inline constexpr Uint128 ShiftLeft1(Uint128 a) noexcept {
  return Uint128{(a.hi << 1) | (a.lo >> 63), a.lo << 1};
}

}  // namespace detail

// Sebastiano Vigna's SplitMix64 generator. Satisfies
// std::uniform_random_bit_generator and provides a bijective 64-bit mixer.
class SplitMix64 {
 public:
  using result_type = uint64_t;

  static constexpr uint64_t kGamma = 0x9e3779b97f4a7c15ULL;
  static constexpr uint64_t kMul1 = 0xbf58476d1ce4e5b9ULL;
  static constexpr uint64_t kMul2 = 0x94d049bb133111ebULL;

  explicit constexpr SplitMix64(uint64_t seed = kGamma) noexcept
      : state_(seed) {}

  static constexpr result_type min() noexcept {
    return std::numeric_limits<result_type>::min();
  }

  static constexpr result_type max() noexcept {
    return std::numeric_limits<result_type>::max();
  }

  // Bijective Stafford Mix13 finalizer on uint64_t.
  static constexpr uint64_t Mix(uint64_t x) noexcept {
    x = (x ^ (x >> 30)) * kMul1;
    x = (x ^ (x >> 27)) * kMul2;
    return x ^ (x >> 31);
  }

  constexpr result_type operator()() noexcept {
    state_ += kGamma;
    return Mix(state_);
  }

  constexpr void Seed(uint64_t seed) noexcept { state_ = seed; }

  constexpr uint64_t State() const noexcept { return state_; }

  constexpr void Discard(uint64_t delta) noexcept { state_ += delta * kGamma; }

  // Uniform IEEE-754 double in [0, 1) with 53 bits of precision.
  double NextUniform01() noexcept {
    return static_cast<double>((*this)() >> 11) * 0x1.0p-53;
  }

  constexpr bool operator==(const SplitMix64& other) const noexcept = default;

 private:
  uint64_t state_;
};

// PCG-XSL-RR-128/64 (128-bit LCG state, 64-bit output, 2^127 independent
// streams) with exact O(log delta) jump-ahead Discard().
class Pcg64 {
 public:
  using result_type = uint64_t;

  // Canonical PCG-128 multiplier from O'Neill (2014).
  static constexpr Uint128 kMultiplier{2549297995355413924ULL,
                                       4865540595714422341ULL};
  static constexpr uint64_t kDefaultSeed = 0x853c49e6748fea9bULL;
  static constexpr uint64_t kDefaultStream = 0xda3e39cb94b95bdbULL;

  explicit constexpr Pcg64(uint64_t seed = kDefaultSeed,
                           uint64_t stream = kDefaultStream) noexcept {
    Seed(seed, stream);
  }

  constexpr Pcg64(Uint128 seed, Uint128 stream) noexcept { Seed(seed, stream); }

  static constexpr result_type min() noexcept {
    return std::numeric_limits<result_type>::min();
  }

  static constexpr result_type max() noexcept {
    return std::numeric_limits<result_type>::max();
  }

  constexpr void Seed(uint64_t seed,
                      uint64_t stream = kDefaultStream) noexcept {
    SplitMix64 sm(seed);
    const Uint128 init_state{sm(), sm()};
    const Uint128 stream128{0ULL, stream};
    Seed(init_state, stream128);
  }

  // Full 128-bit state and 127-bit stream initialization (2^127 independent
  // streams; the odd increment is (stream << 1) | 1 mod 2^128).
  constexpr void Seed(Uint128 init_state, Uint128 stream) noexcept {
    state_ = Uint128{0ULL, 0ULL};
    inc_ = detail::ShiftLeft1(stream);
    inc_.lo |= uint64_t{1};
    Step();
    state_ = detail::Add128(state_, init_state);
    Step();
  }

  constexpr result_type operator()() noexcept {
    Step();
    return Output(state_);
  }

  // Exact O(log delta) 128-bit LCG jump-ahead (Brown 1994 / O'Neill 2014).
  // Matches `delta` sequential calls to operator()() bit-for-bit.
  constexpr void Discard(uint64_t delta) noexcept {
    Uint128 cur_mult = kMultiplier;
    Uint128 cur_plus = inc_;
    Uint128 acc_mult{0ULL, 1ULL};
    Uint128 acc_plus{0ULL, 0ULL};

    while (delta > 0) {
      if ((delta & uint64_t{1}) != 0) {
        acc_mult = detail::Mul128(acc_mult, cur_mult);
        acc_plus = detail::Add128(detail::Mul128(acc_plus, cur_mult), cur_plus);
      }
      cur_plus = detail::Mul128(detail::Add128(cur_mult, Uint128{0ULL, 1ULL}),
                                cur_plus);
      cur_mult = detail::Mul128(cur_mult, cur_mult);
      delta >>= 1;
    }

    state_ = detail::Add128(detail::Mul128(acc_mult, state_), acc_plus);
  }

  constexpr Uint128 State() const noexcept { return state_; }
  constexpr Uint128 Increment() const noexcept { return inc_; }

  // Uniform IEEE-754 double in [0, 1) with 53 bits of precision.
  double NextUniform01() noexcept {
    return static_cast<double>((*this)() >> 11) * 0x1.0p-53;
  }

  // Generates two independent standard normal N(0, 1) variates via Box-Muller.
  // Deterministic across platforms and C++ standard library implementations.
  std::pair<double, double> NextGaussianPair() noexcept {
    const double u1 = (static_cast<double>((*this)() >> 11) + 0.5) * 0x1.0p-53;
    const double u2 = static_cast<double>((*this)() >> 11) * 0x1.0p-53;
    const double r = std::sqrt(-2.0 * std::log(u1));
    const double theta = 2.0 * std::numbers::pi * u2;
    return {r * std::cos(theta), r * std::sin(theta)};
  }

  // Uniform integer in [0, bound - 1] using Lemire's nearly-divisionless
  // unbiased rejection method.
  uint64_t NextBounded(uint64_t bound) noexcept {
    if (bound <= 1) {
      return 0;
    }
    const uint64_t threshold = (uint64_t{0} - bound) % bound;
    while (true) {
      const Uint128 prod = detail::Mul64x64((*this)(), bound);
      if (prod.lo >= threshold) {
        return prod.hi;
      }
    }
  }

  constexpr bool operator==(const Pcg64& other) const noexcept = default;

 private:
  constexpr void Step() noexcept {
    state_ = detail::Add128(detail::Mul128(state_, kMultiplier), inc_);
  }

  static constexpr uint64_t Output(Uint128 s) noexcept {
    const uint64_t xored = s.hi ^ s.lo;
    const int rot = static_cast<int>(s.hi >> 58);
    return std::rotr(xored, rot);
  }

  Uint128 state_{0ULL, 0ULL};
  Uint128 inc_{0ULL, 1ULL};
};

// Bijective SplitMix64 seed derivation for independent, deterministic parallel
// simulation streams. For any fixed (base_seed, sweep_idx), the mapping
// chunk_idx -> DeriveStreamSeed(...) is a bijection on uint64_t.
inline constexpr uint64_t DeriveStreamSeed(uint64_t base_seed,
                                           uint64_t sweep_idx,
                                           uint64_t chunk_idx) noexcept {
  uint64_t h = SplitMix64::Mix(base_seed ^ SplitMix64::kGamma);
  h = SplitMix64::Mix(h ^ (sweep_idx + SplitMix64::kMul1));
  h = SplitMix64::Mix(h ^ (chunk_idx + SplitMix64::kMul2));
  return h;
}

static_assert(std::uniform_random_bit_generator<SplitMix64>);
static_assert(std::uniform_random_bit_generator<Pcg64>);

}  // namespace channel
}  // namespace xg

#endif  // XGALOIS_CHANNEL_RNG_HPP_
