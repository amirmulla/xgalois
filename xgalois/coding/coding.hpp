#ifndef XGALOIS_CODING_CODING_HPP_
#define XGALOIS_CODING_CODING_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

// clang-format off
#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoding_result.hpp"

#include "xgalois/coding/bch.hpp"
#include "xgalois/coding/cyclic_code.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/linear_code.hpp"

#include "xgalois/coding/encoder/cyclic_encoder.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/encoder/grs_encoder.hpp"

#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/grs_decoder.hpp"

#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/decoder/cyclic_syndrome_decoder.hpp"
#include "xgalois/coding/decoder/sugiyama_decoder.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"
// clang-format on

namespace xg {
namespace coding {

namespace detail {

// Exact arbitrary-precision unsigned integer (little-endian base-2^64 limbs)
// for overflow-free Hamming bound / perfection checks over arbitrary (n, k, q).
struct ExactBigUint {
  __extension__ using u128 = unsigned __int128;
  std::vector<uint64_t> limbs;

  static ExactBigUint FromU64(uint64_t v) {
    ExactBigUint r;
    if (v != 0) {
      r.limbs.push_back(v);
    }
    return r;
  }

  void Trim() {
    while (!limbs.empty() && limbs.back() == 0) {
      limbs.pop_back();
    }
  }

  void MulU64(uint64_t factor) {
    if (factor == 0 || limbs.empty()) {
      limbs.clear();
      return;
    }
    if (factor == 1) {
      return;
    }
    u128 carry = 0;
    for (uint64_t &limb : limbs) {
      u128 prod = static_cast<u128>(limb) * static_cast<u128>(factor) + carry;
      limb = static_cast<uint64_t>(prod);
      carry = prod >> 64;
    }
    if (carry != 0) {
      limbs.push_back(static_cast<uint64_t>(carry));
    }
  }

  void DivExactU64(uint64_t divisor) {
    if (divisor <= 1 || limbs.empty()) {
      return;
    }
    u128 rem = 0;
    for (std::size_t i = limbs.size(); i-- > 0;) {
      u128 cur = (rem << 64) | static_cast<u128>(limbs[i]);
      limbs[i] = static_cast<uint64_t>(cur / divisor);
      rem = cur % divisor;
    }
    Trim();
  }

  void Add(const ExactBigUint &other) {
    const std::size_t max_len = std::max(limbs.size(), other.limbs.size());
    limbs.resize(max_len, 0);
    u128 carry = 0;
    for (std::size_t i = 0; i < max_len; ++i) {
      const uint64_t b = (i < other.limbs.size()) ? other.limbs[i] : 0;
      u128 sum = static_cast<u128>(limbs[i]) + static_cast<u128>(b) + carry;
      limbs[i] = static_cast<uint64_t>(sum);
      carry = sum >> 64;
    }
    if (carry != 0) {
      limbs.push_back(static_cast<uint64_t>(carry));
    }
  }

  // Returns -1 if *this < other, 0 if equal, +1 if *this > other.
  int Compare(const ExactBigUint &other) const noexcept {
    if (limbs.size() != other.limbs.size()) {
      return (limbs.size() < other.limbs.size()) ? -1 : 1;
    }
    for (std::size_t i = limbs.size(); i-- > 0;) {
      if (limbs[i] != other.limbs[i]) {
        return (limbs[i] < other.limbs[i]) ? -1 : 1;
      }
    }
    return 0;
  }
};

// Compares Hamming ball volume V_q(n, t) = sum_{i=0}^{t} C(n, i) (q - 1)^i
// against q^(n - k) exactly without floating-point rounding or overflow.
// Returns -1 if V < q^(n-k), 0 if V == q^(n-k) (perfect), +1 if V > q^(n-k).
inline int CompareHammingVolumeToRedundancy(uint64_t q, std::size_t n,
                                            std::size_t k, std::size_t t) {
  const std::size_t r = n - k;
  // Compute rhs = q^r. (q == 0 denotes GF(2^64), i.e., q = 2^64.)
  ExactBigUint rhs;
  if (q == 0) {
    rhs.limbs.assign(r + 1, 0);
    rhs.limbs[r] = 1;
  } else {
    rhs = ExactBigUint::FromU64(1);
    for (std::size_t i = 0; i < r; ++i) {
      rhs.MulU64(q);
    }
  }

  const uint64_t q_minus_1 =
      (q == 0) ? std::numeric_limits<uint64_t>::max() : (q - 1);

  ExactBigUint volume = ExactBigUint::FromU64(1);
  ExactBigUint term = ExactBigUint::FromU64(1);
  const std::size_t t_clamped = std::min(t, n);
  for (std::size_t i = 1; i <= t_clamped; ++i) {
    term.MulU64(static_cast<uint64_t>(n - i + 1));
    term.DivExactU64(static_cast<uint64_t>(i));
    term.MulU64(q_minus_1);
    volume.Add(term);
    if (volume.Compare(rhs) > 0) {
      return 1;
    }
  }
  return volume.Compare(rhs);
}

}  // namespace detail

template <typename GaloisField>
GeneralizedReedSolomonCode<GaloisField> CreateReedSolomonCode(
    std::shared_ptr<GaloisField> field, std::size_t length,
    std::size_t dimension,
    const GaloisFieldElement<GaloisField> &primitive_element) {
  return GeneralizedReedSolomonCode<GaloisField>::ReedSolomon(
      std::move(field), length, dimension, primitive_element);
}

template <typename GaloisField>
GeneralizedReedSolomonCode<GaloisField> CreateReedSolomonCode(
    std::shared_ptr<GaloisField> field, std::size_t length,
    std::size_t dimension) {
  return GeneralizedReedSolomonCode<GaloisField>::ReedSolomon(
      std::move(field), length, dimension);
}

template <typename GaloisField>
GeneralizedReedSolomonCode<GaloisField> CreateExtendedReedSolomonCode(
    std::shared_ptr<GaloisField> field, std::size_t length,
    std::size_t dimension) {
  return GeneralizedReedSolomonCode<GaloisField>::ExtendedReedSolomon(
      std::move(field), length, dimension);
}

template <typename GaloisField>
GeneralizedReedSolomonCode<GaloisField> CreateGeneralizedReedSolomonCode(
    std::shared_ptr<GaloisField> field,
    const GFArray<GaloisField> &evaluation_points,
    const GFArray<GaloisField> &column_multipliers, std::size_t dimension) {
  return GeneralizedReedSolomonCode<GaloisField>(
      std::move(field), evaluation_points, column_multipliers, dimension);
}

template <typename GaloisField>
CyclicCode<GaloisField> CreateCyclicCode(
    std::shared_ptr<GaloisField> field, std::size_t length,
    const PolynomialDense<GaloisField> &generator_poly) {
  return CyclicCode<GaloisField>(std::move(field), length, generator_poly);
}

template <typename GaloisField>
std::size_t ErrorCorrectionCapability(const AbstractCode<GaloisField> &code) {
  const std::size_t d = code.MinimumDistance();
  return d > 0 ? (d - 1) / 2 : 0;
}

template <typename GaloisField>
double CodeRate(const AbstractLinearCode<GaloisField> &code) {
  return code.Rate();
}

template <typename GaloisField>
std::size_t Redundancy(const AbstractLinearCode<GaloisField> &code) {
  return code.Redundancy();
}

inline std::size_t SingletonBound(std::size_t length, std::size_t dimension) {
  if (dimension > length) {
    throw std::invalid_argument(
        "SingletonBound: dimension cannot exceed length");
  }
  return length - dimension + 1;
}

template <typename GaloisField>
bool SatisfiesSingletonBound(const AbstractLinearCode<GaloisField> &code) {
  return code.MinimumDistance() <=
         SingletonBound(code.Length(), code.Dimension());
}

template <typename GaloisField>
bool IsMDSCode(const AbstractLinearCode<GaloisField> &code) {
  return code.MinimumDistance() ==
         SingletonBound(code.Length(), code.Dimension());
}

template <typename GaloisField>
bool IsMDS(const AbstractLinearCode<GaloisField> &code) {
  return IsMDSCode(code);
}

template <typename GaloisField>
bool SatisfiesHammingBound(const AbstractLinearCode<GaloisField> &code) {
  const uint64_t q = code.Field()->Order();
  const std::size_t n = code.Length();
  const std::size_t k = code.Dimension();
  const std::size_t t = ErrorCorrectionCapability(code);
  return detail::CompareHammingVolumeToRedundancy(q, n, k, t) <= 0;
}

template <typename GaloisField>
bool IsPerfectCode(const AbstractLinearCode<GaloisField> &code) {
  const uint64_t q = code.Field()->Order();
  const std::size_t n = code.Length();
  const std::size_t k = code.Dimension();
  const std::size_t t = ErrorCorrectionCapability(code);
  return detail::CompareHammingVolumeToRedundancy(q, n, k, t) == 0;
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_CODING_HPP_
