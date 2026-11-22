#ifndef XGALOIS_CODING_LRC_HPP_
#define XGALOIS_CODING_LRC_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/berlekamp_massey_decoder.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

// ============================================================================
// Enums & Result Types
// ============================================================================

enum class LRCPartitionType : uint8_t {
  kMultiplicative,
  kAdditive,
  kCustom,
  kAvailabilityMultiplicative,
  kAvailabilityAdditive,
  kAvailabilityHybrid
};

template <typename GaloisField>
struct LocalRepairResult {
  using element_type = GaloisFieldElement<GaloisField>;

  element_type recovered_symbol;
  std::size_t erased_index = 0;
  std::size_t repair_set_index = 0;
  std::size_t symbols_read = 0;
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }

  static LocalRepairResult Success(element_type symbol, std::size_t erased_idx,
                                   std::size_t repair_set_idx,
                                   std::size_t num_symbols_read) {
    LocalRepairResult res;
    res.recovered_symbol = std::move(symbol);
    res.erased_index = erased_idx;
    res.repair_set_index = repair_set_idx;
    res.symbols_read = num_symbols_read;
    res.success = true;
    return res;
  }

  static LocalRepairResult Failure(std::size_t erased_idx,
                                   std::size_t repair_set_idx,
                                   std::string reason) {
    LocalRepairResult res;
    res.erased_index = erased_idx;
    res.repair_set_index = repair_set_idx;
    res.success = false;
    res.error_message = std::move(reason);
    return res;
  }
};

template <typename GaloisField>
struct MultiLocalRepairResult {
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;

  codeword_type codeword;
  message_type message;
  std::size_t num_local_repairs = 0;
  std::size_t num_global_repairs = 0;
  std::size_t total_symbols_read = 0;
  bool used_global_fallback = false;
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }

  static MultiLocalRepairResult Failure(std::string reason,
                                        std::size_t local_repairs = 0,
                                        std::size_t global_repairs = 0,
                                        std::size_t symbols_read = 0,
                                        bool global_fallback = false) {
    MultiLocalRepairResult res;
    res.num_local_repairs = local_repairs;
    res.num_global_repairs = global_repairs;
    res.total_symbols_read = symbols_read;
    res.used_global_fallback = global_fallback;
    res.success = false;
    res.error_message = std::move(reason);
    return res;
  }
};

template <typename GaloisField>
class LRCGlobalDecoder;

template <typename GaloisField>
class LRCLocalThenGlobalDecoder;

// ============================================================================
// Internal Mathematical Helpers
// ============================================================================

namespace lrc_detail {

inline uint64_t MultiplicativeGroupOrder(uint64_t q) noexcept {
  return (q == 0) ? std::numeric_limits<uint64_t>::max() : (q - 1);
}

inline std::optional<std::size_t> ExactPrimePowerExponent(uint64_t value,
                                                          uint64_t p) noexcept {
  if (p < 2 || value < p) {
    return std::nullopt;
  }
  std::size_t exp = 0;
  while (value > 1) {
    if (value % p != 0) {
      return std::nullopt;
    }
    value /= p;
    ++exp;
  }
  return exp;
}

inline std::size_t CeilDiv(std::size_t a, std::size_t b) noexcept {
  return (a + b - 1) / b;
}

inline std::size_t SafeMulSize(std::size_t a, std::size_t b,
                               const char *context) {
  if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) {
    throw std::invalid_argument(std::string(context) + ": integer overflow");
  }
  return a * b;
}

// Computes c * 1 in F_q for c in {0, ..., p - 1} using O(log c) field
// additions.
template <typename GaloisField>
inline typename GaloisField::element_type ScalarFromPrime(
    const GaloisField &field, uint64_t c) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t p = field.Characteristic();
  c %= p;
  StorageType acc = field.AdditiveIdentity();
  if (c == 0) {
    return acc;
  }
  StorageType base = field.MultiplicativeIdentity();
  while (c > 0) {
    if ((c & 1ULL) != 0) {
      acc = field.Add(acc, base);
    }
    c >>= 1;
    if (c > 0) {
      base = field.Add(base, base);
    }
  }
  return acc;
}

// Extracts an F_p-basis {b_0, ..., b_{m-1}} of F_{p^m} over F_p.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type> ExtractPrimeFieldBasis(
    const GaloisField &field) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t p = field.Characteristic();
  const std::size_t m = static_cast<std::size_t>(field.Degree());
  std::vector<StorageType> basis;
  basis.reserve(m);

  if constexpr (std::is_integral_v<StorageType>) {
    if (p == 2) {
      for (std::size_t i = 0; i < m; ++i) {
        basis.push_back(
            field.SetElementValue(static_cast<StorageType>(uint64_t{1} << i)));
      }
    } else {
      // Prime field GF(p), m == 1.
      basis.push_back(field.MultiplicativeIdentity());
    }
  } else if constexpr (requires(const GaloisField &f, uint64_t v) {
                         f.FromInt(v);
                       }) {
    basis.push_back(field.MultiplicativeIdentity());
    if (m > 1) {
      const StorageType alpha = field.FromInt(p);
      for (std::size_t i = 1; i < m; ++i) {
        basis.push_back(field.Mul(basis.back(), alpha));
      }
    }
  } else {
    basis.push_back(field.MultiplicativeIdentity());
    if (m > 1) {
      const StorageType alpha = field.MultiplicativeGenerator();
      for (std::size_t i = 1; i < m; ++i) {
        basis.push_back(field.Mul(basis.back(), alpha));
      }
    }
  }
  return basis;
}

// Decomposes x in F_{p^m} into its length-m coordinate vector over F_p
// relative to ExtractPrimeFieldBasis(field).
template <typename GaloisField>
inline std::vector<uint64_t> ToPrimeFieldCoords(
    const GaloisField &field, const typename GaloisField::element_type &x) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t p = field.Characteristic();
  const std::size_t m = static_cast<std::size_t>(field.Degree());
  std::vector<uint64_t> coords(m, 0);

  if constexpr (std::is_integral_v<StorageType>) {
    const uint64_t ext = static_cast<uint64_t>(field.GetElementValue(x));
    if (p == 2) {
      for (std::size_t i = 0; i < m; ++i) {
        coords[i] = (ext >> i) & 1ULL;
      }
    } else {
      coords[0] = ext % p;
    }
  } else if constexpr (requires(const StorageType &e, std::size_t idx) {
                         e[idx];
                       }) {
    for (std::size_t i = 0; i < m; ++i) {
      coords[i] = static_cast<uint64_t>(x[i]) % p;
    }
  } else {
    throw std::invalid_argument(
        "ToPrimeFieldCoords: unsupported field storage type");
  }
  return coords;
}

// Enumerates the first `count` elements of the F_p-span of `basis_subset`
// in base-p lexicographic order (index 0 is always the additive identity 0).
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type> EnumeratePrimeSpan(
    const GaloisField &field,
    std::span<const typename GaloisField::element_type> basis_subset,
    std::size_t count) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t p = field.Characteristic();
  const StorageType zero = field.AdditiveIdentity();

  std::vector<StorageType> elements;
  elements.reserve(count);
  elements.push_back(zero);
  if (count <= 1 || basis_subset.empty()) {
    return elements;
  }

  // Precompute scalar multiples c * b_i for c = 1, ..., p - 1 on the fly.
  for (std::size_t b_idx = 0; b_idx < basis_subset.size(); ++b_idx) {
    const std::size_t cur_size = elements.size();
    if (cur_size >= count) {
      break;
    }
    StorageType shift = basis_subset[b_idx];
    for (uint64_t c = 1; c < p && elements.size() < count; ++c) {
      for (std::size_t j = 0; j < cur_size && elements.size() < count; ++j) {
        elements.push_back(field.Add(elements[j], shift));
      }
      if (c + 1 < p) {
        shift = field.Add(shift, basis_subset[b_idx]);
      }
    }
  }
  return elements;
}

// Computes the linearized vanishing polynomial g(x) = prod_{beta in H}(x -
// beta) of the F_p-subspace H = span_{F_p}{e_0, ..., e_{a-1}} in O(a^2) field
// operations via the recursion:
//   g_0(x) = x,
//   g_t(x) = g_{t-1}(x)^p - g_{t-1}(e_{t-1})^(p-1) * g_{t-1}(x).
template <typename GaloisField>
inline PolynomialDense<GaloisField> ComputeLinearizedVanishingPoly(
    const std::shared_ptr<GaloisField> &field,
    std::span<const typename GaloisField::element_type> subspace_basis) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t p = field->Characteristic();
  const std::size_t a = subspace_basis.size();
  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();

  // lin_coeffs[j] is the coefficient of x^(p^j).
  std::vector<StorageType> lin_coeffs = {one};
  lin_coeffs.reserve(a + 1);

  for (std::size_t t = 0; t < a; ++t) {
    // Evaluate g_t(e_t) = sum_{j=0}^{t} lin_coeffs[j] * e_t^(p^j).
    StorageType val = zero;
    StorageType e_pow = subspace_basis[t];
    for (std::size_t j = 0; j <= t; ++j) {
      if (!field->IsZero(lin_coeffs[j])) {
        val = field->Add(val, field->Mul(lin_coeffs[j], e_pow));
      }
      if (j < t) {
        e_pow = field->Frobenius(e_pow);
      }
    }
    if (field->IsZero(val)) {
      throw std::invalid_argument(
          "ComputeLinearizedVanishingPoly: subspace basis vectors are linearly "
          "dependent over F_p");
    }
    const StorageType c = field->Pow(val, p - 1);
    const StorageType neg_c = field->Neg(c);

    std::vector<StorageType> next_coeffs(t + 2, zero);
    next_coeffs[0] = field->Mul(neg_c, lin_coeffs[0]);
    for (std::size_t j = 1; j <= t; ++j) {
      const StorageType frob_prev = field->Frobenius(lin_coeffs[j - 1]);
      next_coeffs[j] = field->Add(frob_prev, field->Mul(neg_c, lin_coeffs[j]));
    }
    next_coeffs[t + 1] = field->Frobenius(lin_coeffs[t]);
    lin_coeffs = std::move(next_coeffs);
  }

  const uint64_t deg = utils::SafeIntegerPower(p, a);
  std::vector<StorageType> dense_coeffs(static_cast<std::size_t>(deg) + 1,
                                        zero);
  uint64_t p_pow = 1;
  for (std::size_t j = 0; j <= a; ++j) {
    dense_coeffs[static_cast<std::size_t>(p_pow)] = lin_coeffs[j];
    if (j < a) {
      p_pow *= p;
    }
  }
  return PolynomialDense<GaloisField>::FromRaw(std::move(dense_coeffs), field);
}

// Computes exact Lagrange evaluation weights w_a = L_a(target_x) for nodes
// xs[0..r-1] using barycentric form and a single BatchInv:
//   p(target_x) = sum_{a=0}^{r-1} w_a * p(xs[a])  for all deg(p) < r.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type>
ComputeLagrangeEvalWeights(
    const GaloisField &field,
    std::span<const typename GaloisField::element_type> xs,
    const typename GaloisField::element_type &target_x) {
  using StorageType = typename GaloisField::element_type;
  const std::size_t r = xs.size();
  if (r == 0) {
    return {};
  }
  const StorageType one = field.MultiplicativeIdentity();
  std::vector<StorageType> denoms(r);
  StorageType num_prod = one;

  for (std::size_t a = 0; a < r; ++a) {
    const StorageType diff_target = field.Sub(target_x, xs[a]);
    if (field.IsZero(diff_target)) {
      throw std::invalid_argument(
          "ComputeLagrangeEvalWeights: target point coincides with node");
    }
    num_prod = field.Mul(num_prod, diff_target);

    StorageType d = diff_target;
    for (std::size_t b = 0; b < r; ++b) {
      if (b != a) {
        d = field.Mul(d, field.Sub(xs[a], xs[b]));
      }
    }
    denoms[a] = d;
  }

  std::vector<StorageType> inv_denoms(r);
  field.BatchInv(std::span<StorageType>(inv_denoms),
                 std::span<const StorageType>(denoms));

  std::vector<StorageType> weights(r);
  field.ScalarMulArray(std::span<StorageType>(weights),
                       std::span<const StorageType>(inv_denoms), num_prod);
  return weights;
}

// Computes barycentric node weights beta_a = (prod_{b != a}(xs[a] - xs[b]))^-1
// using a single BatchInv.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type>
ComputeBarycentricNodeWeights(
    const GaloisField &field,
    std::span<const typename GaloisField::element_type> xs) {
  using StorageType = typename GaloisField::element_type;
  const std::size_t r = xs.size();
  if (r == 0) {
    return {};
  }
  const StorageType one = field.MultiplicativeIdentity();
  std::vector<StorageType> denoms(r, one);
  for (std::size_t a = 0; a < r; ++a) {
    StorageType d = one;
    for (std::size_t b = 0; b < r; ++b) {
      if (b != a) {
        d = field.Mul(d, field.Sub(xs[a], xs[b]));
      }
    }
    denoms[a] = d;
  }
  std::vector<StorageType> betas(r);
  field.BatchInv(std::span<StorageType>(betas),
                 std::span<const StorageType>(denoms));
  return betas;
}

// Helper for hashing/checking distinctness across arbitrary storage_type.
template <typename StorageType>
struct RawElementIndexMap {
  static constexpr bool kHashable = requires(const StorageType &x) {
    { std::hash<StorageType>{}(x) } -> std::convertible_to<std::size_t>;
  };

  std::unordered_map<StorageType, std::size_t> map_;
  std::vector<std::pair<StorageType, std::size_t>> fallback_;

  bool Contains(const StorageType &val) const {
    if constexpr (kHashable) {
      return map_.find(val) != map_.end();
    } else {
      for (const auto &[k, _] : fallback_) {
        if (k == val) return true;
      }
      return false;
    }
  }

  std::optional<std::size_t> Find(const StorageType &val) const {
    if constexpr (kHashable) {
      auto it = map_.find(val);
      if (it == map_.end()) return std::nullopt;
      return it->second;
    } else {
      for (const auto &[k, idx] : fallback_) {
        if (k == val) return idx;
      }
      return std::nullopt;
    }
  }

  bool Insert(const StorageType &val, std::size_t idx) {
    if constexpr (kHashable) {
      return map_.emplace(val, idx).second;
    } else {
      if (Contains(val)) return false;
      fallback_.emplace_back(val, idx);
      return true;
    }
  }

  void Reserve(std::size_t n) {
    if constexpr (kHashable) {
      map_.reserve(n);
    } else {
      fallback_.reserve(n);
    }
  }
};

// Finds a deterministic K-basis {v_0 = 1, v_1, ..., v_{d0 - 1}} of F_{p^m}
// over subfield K = F_{p^b} (where b | m, d0 = m / b), given an F_p-basis
// k_fp_basis = {u_0, ..., u_{b-1}} of K.
template <typename GaloisField>
inline std::vector<typename GaloisField::element_type> FindSubfieldBasis(
    const GaloisField &field,
    std::span<const typename GaloisField::element_type> k_fp_basis,
    std::size_t b) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t p = field.Characteristic();
  const std::size_t m = static_cast<std::size_t>(field.Degree());
  const std::size_t d0 = m / b;

  std::vector<std::vector<uint64_t>> echelon_rows;
  std::vector<std::size_t> pivot_cols;
  echelon_rows.reserve(m);
  pivot_cols.reserve(m);

  auto try_insert_fp_vec = [&](std::vector<uint64_t> vec) -> bool {
    for (std::size_t r = 0; r < echelon_rows.size(); ++r) {
      const std::size_t pc = pivot_cols[r];
      const uint64_t factor = vec[pc];
      if (factor == 0) {
        continue;
      }
      for (std::size_t c = pc; c < m; ++c) {
        const uint64_t sub = utils::ModMul(factor, echelon_rows[r][c], p);
        vec[c] = (vec[c] + p - sub) % p;
      }
    }
    std::size_t new_pivot = m;
    for (std::size_t c = 0; c < m; ++c) {
      if (vec[c] != 0) {
        new_pivot = c;
        break;
      }
    }
    if (new_pivot == m) {
      return false;
    }
    if (vec[new_pivot] != 1) {
      const uint64_t inv = utils::ModInv(vec[new_pivot], p);
      for (std::size_t c = new_pivot; c < m; ++c) {
        vec[c] = utils::ModMul(vec[c], inv, p);
      }
    }
    echelon_rows.push_back(std::move(vec));
    pivot_cols.push_back(new_pivot);
    return true;
  };

  auto is_in_span = [&](std::vector<uint64_t> vec) -> bool {
    for (std::size_t r = 0; r < echelon_rows.size(); ++r) {
      const std::size_t pc = pivot_cols[r];
      const uint64_t factor = vec[pc];
      if (factor == 0) {
        continue;
      }
      for (std::size_t c = pc; c < m; ++c) {
        const uint64_t sub = utils::ModMul(factor, echelon_rows[r][c], p);
        vec[c] = (vec[c] + p - sub) % p;
      }
    }
    for (uint64_t v : vec) {
      if (v != 0) {
        return false;
      }
    }
    return true;
  };

  std::vector<StorageType> k_basis;
  k_basis.reserve(d0);

  auto add_k_basis_vec = [&](const StorageType &v) {
    k_basis.push_back(v);
    for (const StorageType &u : k_fp_basis) {
      const StorageType prod = field.Mul(u, v);
      try_insert_fp_vec(ToPrimeFieldCoords(field, prod));
    }
  };

  add_k_basis_vec(field.MultiplicativeIdentity());
  const std::vector<StorageType> fp_basis = ExtractPrimeFieldBasis(field);
  for (const StorageType &cand : fp_basis) {
    if (k_basis.size() == d0) {
      break;
    }
    if (!is_in_span(ToPrimeFieldCoords(field, cand))) {
      add_k_basis_vec(cand);
    }
  }
  if (k_basis.size() != d0) {
    throw std::runtime_error(
        "FindSubfieldBasis: failed to construct K-basis of F_{p^m}");
  }
  return k_basis;
}

}  // namespace lrc_detail

// ============================================================================
// LocallyRecoverableCode<GaloisField>
// ============================================================================

template <typename GaloisField>
class LocallyRecoverableCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = GaloisFieldElement<GaloisField>;
  using storage_type = typename GaloisField::element_type;
  using codeword_type = GFArray<GaloisField>;
  using message_type = GFArray<GaloisField>;
  using matrix_type = GFMatrix<GaloisField>;
  using vector_type = GFArray<GaloisField>;
  using polynomial_type = PolynomialDense<GaloisField>;

  // Rule of Five: preserve cache safety across copies and moves so no cached
  // encoder or decoder ever references a moved-from or destroyed code object.
  LocallyRecoverableCode(const LocallyRecoverableCode &other)
      : Base(other),
        partition_type_(other.partition_type_),
        localities_(other.localities_),
        local_distances_(other.local_distances_),
        block_sizes_(other.block_sizes_),
        availability_(other.availability_),
        designed_distance_(other.designed_distance_),
        max_poly_degree_(other.max_poly_degree_),
        has_univariate_ambient_rs_(other.has_univariate_ambient_rs_),
        eval_points_(other.eval_points_),
        good_poly_(other.good_poly_),
        fiber_partitions_(other.fiber_partitions_),
        coord_to_fiber_(other.coord_to_fiber_),
        repair_sets_(other.repair_sets_),
        local_eval_points_(other.local_eval_points_),
        canonical_weights_(other.canonical_weights_),
        generator_matrix_(other.generator_matrix_),
        parity_check_matrix_(other.parity_check_matrix_) {
    this->ClearCaches();
    RegisterEncodersAndDecoders();
  }

  LocallyRecoverableCode &operator=(const LocallyRecoverableCode &other) {
    if (this != &other) {
      Base::operator=(other);
      partition_type_ = other.partition_type_;
      localities_ = other.localities_;
      local_distances_ = other.local_distances_;
      block_sizes_ = other.block_sizes_;
      availability_ = other.availability_;
      designed_distance_ = other.designed_distance_;
      max_poly_degree_ = other.max_poly_degree_;
      has_univariate_ambient_rs_ = other.has_univariate_ambient_rs_;
      eval_points_ = other.eval_points_;
      good_poly_ = other.good_poly_;
      fiber_partitions_ = other.fiber_partitions_;
      coord_to_fiber_ = other.coord_to_fiber_;
      repair_sets_ = other.repair_sets_;
      local_eval_points_ = other.local_eval_points_;
      canonical_weights_ = other.canonical_weights_;
      generator_matrix_ = other.generator_matrix_;
      parity_check_matrix_ = other.parity_check_matrix_;
      this->ClearCaches();
      RegisterEncodersAndDecoders();
    }
    return *this;
  }

  LocallyRecoverableCode(LocallyRecoverableCode &&other) noexcept
      : Base(std::move(other)),
        partition_type_(other.partition_type_),
        localities_(std::move(other.localities_)),
        local_distances_(std::move(other.local_distances_)),
        block_sizes_(std::move(other.block_sizes_)),
        availability_(other.availability_),
        designed_distance_(other.designed_distance_),
        max_poly_degree_(other.max_poly_degree_),
        has_univariate_ambient_rs_(other.has_univariate_ambient_rs_),
        eval_points_(std::move(other.eval_points_)),
        good_poly_(std::move(other.good_poly_)),
        fiber_partitions_(std::move(other.fiber_partitions_)),
        coord_to_fiber_(std::move(other.coord_to_fiber_)),
        repair_sets_(std::move(other.repair_sets_)),
        local_eval_points_(std::move(other.local_eval_points_)),
        canonical_weights_(std::move(other.canonical_weights_)),
        generator_matrix_(std::move(other.generator_matrix_)),
        parity_check_matrix_(std::move(other.parity_check_matrix_)) {
    this->ClearCaches();
  }

  LocallyRecoverableCode &operator=(LocallyRecoverableCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      partition_type_ = other.partition_type_;
      localities_ = std::move(other.localities_);
      local_distances_ = std::move(other.local_distances_);
      block_sizes_ = std::move(other.block_sizes_);
      availability_ = other.availability_;
      designed_distance_ = other.designed_distance_;
      max_poly_degree_ = other.max_poly_degree_;
      has_univariate_ambient_rs_ = other.has_univariate_ambient_rs_;
      eval_points_ = std::move(other.eval_points_);
      good_poly_ = std::move(other.good_poly_);
      fiber_partitions_ = std::move(other.fiber_partitions_);
      coord_to_fiber_ = std::move(other.coord_to_fiber_);
      repair_sets_ = std::move(other.repair_sets_);
      local_eval_points_ = std::move(other.local_eval_points_);
      canonical_weights_ = std::move(other.canonical_weights_);
      generator_matrix_ = std::move(other.generator_matrix_);
      parity_check_matrix_ = std::move(other.parity_check_matrix_);
      this->ClearCaches();
    }
    return *this;
  }

  ~LocallyRecoverableCode() override = default;

  // ==========================================================================
  // Factory 1: Optimal Single-Set (n, k, r) and Generalized (r, delta)
  //            Tamo-Barg LRC
  // ==========================================================================

  static LocallyRecoverableCode CreateTamoBarg(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension, std::size_t locality,
      std::size_t local_distance = 2,
      std::optional<LRCPartitionType> partition_type = std::nullopt) {
    if (!field) {
      throw std::invalid_argument("CreateTamoBarg: field cannot be null");
    }
    if (length == 0) {
      throw std::invalid_argument("CreateTamoBarg: length must be positive");
    }
    if (dimension == 0 || dimension > length) {
      throw std::invalid_argument(
          "CreateTamoBarg: dimension must satisfy 1 <= k <= n");
    }
    if (locality == 0) {
      throw std::invalid_argument("CreateTamoBarg: locality r must be >= 1");
    }
    if (local_distance < 2) {
      throw std::invalid_argument(
          "CreateTamoBarg: local_distance delta must be >= 2");
    }

    const std::size_t s = locality + local_distance - 1;
    if (length % s != 0) {
      throw std::invalid_argument(
          "CreateTamoBarg: block size s = r + delta - 1 (" + std::to_string(s) +
          ") must divide code length n (" + std::to_string(length) + ")");
    }
    const std::size_t num_blocks = length / s;
    const std::size_t msg_blocks = lrc_detail::CeilDiv(dimension, locality);
    if (msg_blocks > num_blocks) {
      throw std::invalid_argument("CreateTamoBarg: ceil(k / r) (" +
                                  std::to_string(msg_blocks) +
                                  ") exceeds number of repair groups n / s (" +
                                  std::to_string(num_blocks) + ")");
    }

    const uint64_t q = field->Order();
    const uint64_t p = field->Characteristic();
    const std::size_t m = static_cast<std::size_t>(field->Degree());
    const uint64_t group_order = lrc_detail::MultiplicativeGroupOrder(q);

    const bool can_mult =
        (group_order % static_cast<uint64_t>(s) == 0) &&
        (q == 0 || static_cast<uint64_t>(length) <= group_order);

    const auto add_exp =
        lrc_detail::ExactPrimePowerExponent(static_cast<uint64_t>(s), p);
    const bool can_add = add_exp.has_value() && (*add_exp <= m) &&
                         (q == 0 || static_cast<uint64_t>(length) <= q);

    LRCPartitionType chosen = LRCPartitionType::kMultiplicative;
    if (partition_type.has_value()) {
      chosen = *partition_type;
      if (chosen == LRCPartitionType::kMultiplicative) {
        if (!can_mult) {
          throw std::invalid_argument(
              "CreateTamoBarg: multiplicative partition requires s | (q - 1) "
              "and n <= q - 1");
        }
      } else if (chosen == LRCPartitionType::kAdditive) {
        if (!can_add) {
          throw std::invalid_argument(
              "CreateTamoBarg: additive partition requires s = p^a (a <= m) "
              "and n <= q");
        }
      } else {
        throw std::invalid_argument(
            "CreateTamoBarg: partition_type must be kMultiplicative or "
            "kAdditive");
      }
    } else {
      if (can_mult) {
        chosen = LRCPartitionType::kMultiplicative;
      } else if (can_add) {
        chosen = LRCPartitionType::kAdditive;
      } else {
        throw std::invalid_argument(
            "CreateTamoBarg: neither multiplicative (s | q-1, n <= q-1) nor "
            "additive (s = p^a, n <= q) partition is feasible for this field");
      }
    }

    std::vector<storage_type> eval_raw_vec;
    eval_raw_vec.reserve(length);
    std::vector<storage_type> fiber_constants;
    fiber_constants.reserve(num_blocks);
    polynomial_type good_poly;

    if (chosen == LRCPartitionType::kMultiplicative) {
      const storage_type gamma = field->MultiplicativeGenerator();
      const storage_type omega =
          field->Pow(gamma, group_order / static_cast<uint64_t>(s));
      const storage_type one = field->MultiplicativeIdentity();

      std::vector<storage_type> subgroup_h(s);
      subgroup_h[0] = one;
      for (std::size_t t = 1; t < s; ++t) {
        subgroup_h[t] = field->Mul(subgroup_h[t - 1], omega);
      }

      storage_type coset_rep = one;
      for (std::size_t ell = 0; ell < num_blocks; ++ell) {
        for (std::size_t t = 0; t < s; ++t) {
          eval_raw_vec.push_back(field->Mul(coset_rep, subgroup_h[t]));
        }
        fiber_constants.push_back(field->Pow(coset_rep, s));
        if (ell + 1 < num_blocks) {
          coset_rep = field->Mul(coset_rep, gamma);
        }
      }
      good_poly = polynomial_type::Monomial(field, s, one);
    } else {
      const std::size_t a = *add_exp;
      const std::vector<storage_type> fp_basis =
          lrc_detail::ExtractPrimeFieldBasis(*field);
      const std::span<const storage_type> h_basis(fp_basis.data(), a);
      const std::span<const storage_type> w_basis(fp_basis.data() + a, m - a);

      const std::vector<storage_type> subgroup_h =
          lrc_detail::EnumeratePrimeSpan(*field, h_basis, s);
      const std::vector<storage_type> coset_reps =
          lrc_detail::EnumeratePrimeSpan(*field, w_basis, num_blocks);

      good_poly = lrc_detail::ComputeLinearizedVanishingPoly(field, h_basis);
      for (std::size_t ell = 0; ell < num_blocks; ++ell) {
        const storage_type &beta = coset_reps[ell];
        for (std::size_t t = 0; t < s; ++t) {
          eval_raw_vec.push_back(field->Add(beta, subgroup_h[t]));
        }
        fiber_constants.push_back(good_poly.EvalRaw(beta));
      }
    }

    std::vector<std::vector<std::size_t>> fibers(num_blocks);
    for (std::size_t ell = 0; ell < num_blocks; ++ell) {
      fibers[ell].resize(s);
      for (std::size_t t = 0; t < s; ++t) {
        fibers[ell][t] = ell * s + t;
      }
    }

    return BuildSingleSetCode(std::move(field), chosen, std::move(eval_raw_vec),
                              std::move(good_poly), std::move(fibers),
                              fiber_constants, dimension, locality,
                              local_distance);
  }

  // ==========================================================================
  // Factory 2: Custom Good Polynomial LRC
  // ==========================================================================

  static LocallyRecoverableCode CreateFromGoodPolynomial(
      std::shared_ptr<GaloisField> field, const polynomial_type &good_poly,
      const std::vector<std::vector<element_type>> &partition_fibers,
      std::size_t dimension, std::size_t locality) {
    if (!field) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: field cannot be null");
    }
    if (!good_poly.HasField() || good_poly.Field() != field) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: good_poly field mismatch");
    }
    if (good_poly.IsZero() || good_poly.Degree() < 2) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: good_poly must have degree s >= 2");
    }
    const std::size_t s = static_cast<std::size_t>(good_poly.Degree());
    if (partition_fibers.empty()) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: partition_fibers cannot be empty");
    }
    if (locality == 0 || locality >= s) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: locality r must satisfy 1 <= r < deg(g)");
    }
    const std::size_t num_blocks = partition_fibers.size();
    const std::size_t length =
        lrc_detail::SafeMulSize(num_blocks, s, "CreateFromGoodPolynomial");
    const std::size_t local_distance = s - locality + 1;
    if (dimension == 0 ||
        lrc_detail::CeilDiv(dimension, locality) > num_blocks) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: dimension must satisfy 1 <= k <= "
          "num_blocks * r");
    }

    std::vector<storage_type> eval_raw_vec;
    eval_raw_vec.reserve(length);
    std::vector<storage_type> fiber_constants;
    fiber_constants.reserve(num_blocks);
    std::vector<std::vector<std::size_t>> fibers(num_blocks);

    lrc_detail::RawElementIndexMap<storage_type> seen_points;
    seen_points.Reserve(length);
    lrc_detail::RawElementIndexMap<storage_type> seen_constants;
    seen_constants.Reserve(num_blocks);

    for (std::size_t ell = 0; ell < num_blocks; ++ell) {
      const auto &fiber = partition_fibers[ell];
      if (fiber.size() != s) {
        throw std::invalid_argument(
            "CreateFromGoodPolynomial: every fiber must have size equal to "
            "deg(g) = " +
            std::to_string(s));
      }
      fibers[ell].reserve(s);
      storage_type fiber_val = field->AdditiveIdentity();
      for (std::size_t t = 0; t < s; ++t) {
        if (!fiber[t].HasField() || fiber[t].Field() != field) {
          throw std::invalid_argument(
              "CreateFromGoodPolynomial: fiber element field mismatch");
        }
        const storage_type pt = fiber[t].RawValue();
        const std::size_t global_idx = eval_raw_vec.size();
        if (!seen_points.Insert(pt, global_idx)) {
          throw std::invalid_argument(
              "CreateFromGoodPolynomial: evaluation points must be globally "
              "distinct across all fibers");
        }
        const storage_type g_val = good_poly.EvalRaw(pt);
        if (t == 0) {
          fiber_val = g_val;
        } else if (g_val != fiber_val) {
          throw std::invalid_argument(
              "CreateFromGoodPolynomial: good_poly is not constant on fiber " +
              std::to_string(ell));
        }
        eval_raw_vec.push_back(pt);
        fibers[ell].push_back(global_idx);
      }
      if (!seen_constants.Insert(fiber_val, ell)) {
        throw std::invalid_argument(
            "CreateFromGoodPolynomial: fiber evaluation values g(A_l) must be "
            "distinct across fibers");
      }
      fiber_constants.push_back(fiber_val);
    }

    return BuildSingleSetCode(std::move(field), LRCPartitionType::kCustom,
                              std::move(eval_raw_vec), good_poly,
                              std::move(fibers), fiber_constants, dimension,
                              locality, local_distance);
  }

  static LocallyRecoverableCode CreateFromGoodPolynomial(
      std::shared_ptr<GaloisField> field, const polynomial_type &good_poly,
      const std::vector<vector_type> &partition_fibers, std::size_t dimension,
      std::size_t locality) {
    if (!field) {
      throw std::invalid_argument(
          "CreateFromGoodPolynomial: field cannot be null");
    }
    std::vector<std::vector<element_type>> converted;
    converted.reserve(partition_fibers.size());
    for (const auto &fiber_arr : partition_fibers) {
      if (!fiber_arr.HasField() || fiber_arr.Field() != field ||
          fiber_arr.Dimension() != 1) {
        throw std::invalid_argument(
            "CreateFromGoodPolynomial: each fiber must be a 1D GFArray over "
            "the code's field");
      }
      std::vector<element_type> fiber_elems;
      fiber_elems.reserve(fiber_arr.Size());
      for (const auto &raw_v : fiber_arr.RawSpan()) {
        fiber_elems.push_back(element_type::FromInternal(raw_v, field));
      }
      converted.push_back(std::move(fiber_elems));
    }
    return CreateFromGoodPolynomial(std::move(field), good_poly, converted,
                                    dimension, locality);
  }

  // ==========================================================================
  // Factory 3a: Multi-Set Multiplicative Availability LRC (CRT Monomial Basis)
  // ==========================================================================

  static LocallyRecoverableCode CreateMultiplicativeAvailability(
      std::shared_ptr<GaloisField> field, std::size_t length,
      std::size_t dimension, const std::vector<std::size_t> &localities,
      const std::vector<std::size_t> &local_distances = {}) {
    if (!field) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: field cannot be null");
    }
    const std::size_t t = localities.size();
    if (t < 2) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: availability t must be >= 2");
    }

    std::vector<std::size_t> deltas = local_distances;
    if (deltas.empty()) {
      deltas.assign(t, std::size_t{2});
    } else if (deltas.size() != t) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: local_distances size must match "
          "localities size");
    }

    std::vector<std::size_t> block_sizes(t);
    std::size_t superblock_s = 1;
    for (std::size_t u = 0; u < t; ++u) {
      if (localities[u] == 0) {
        throw std::invalid_argument(
            "CreateMultiplicativeAvailability: every locality r_u must be >= "
            "1");
      }
      if (deltas[u] < 2) {
        throw std::invalid_argument(
            "CreateMultiplicativeAvailability: every local distance delta_u "
            "must be >= 2");
      }
      block_sizes[u] = localities[u] + deltas[u] - 1;
      for (std::size_t v = 0; v < u; ++v) {
        if (std::gcd(block_sizes[u], block_sizes[v]) != 1) {
          throw std::invalid_argument(
              "CreateMultiplicativeAvailability: block sizes s_u = r_u + "
              "delta_u - 1 must be pairwise coprime");
        }
      }
      superblock_s = lrc_detail::SafeMulSize(
          superblock_s, block_sizes[u], "CreateMultiplicativeAvailability");
    }

    const uint64_t q = field->Order();
    const uint64_t group_order = lrc_detail::MultiplicativeGroupOrder(q);
    if (group_order % static_cast<uint64_t>(superblock_s) != 0) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: S = prod(s_u) must divide q - 1");
    }
    if (length == 0 || length % superblock_s != 0) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: S = prod(s_u) must divide code "
          "length n");
    }
    if (q != 0 && static_cast<uint64_t>(length) > group_order) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: code length cannot exceed q - 1");
    }

    // Collect CRT monomial exponents E = {e < n : forall u, (e % s_u) < r_u}.
    std::vector<std::size_t> crt_exponents;
    for (std::size_t e = 0; e < length; ++e) {
      bool ok = true;
      for (std::size_t u = 0; u < t; ++u) {
        if ((e % block_sizes[u]) >= localities[u]) {
          ok = false;
          break;
        }
      }
      if (ok) {
        crt_exponents.push_back(e);
      }
    }
    if (dimension == 0 || dimension > crt_exponents.size()) {
      throw std::invalid_argument(
          "CreateMultiplicativeAvailability: dimension must satisfy 1 <= k <= "
          "|E| = " +
          std::to_string(crt_exponents.size()));
    }
    crt_exponents.resize(dimension);

    const std::size_t num_superblocks = length / superblock_s;
    const storage_type one = field->MultiplicativeIdentity();
    const storage_type gamma = field->MultiplicativeGenerator();
    const storage_type omega_big =
        field->Pow(gamma, group_order / static_cast<uint64_t>(superblock_s));

    std::vector<storage_type> subgroup_h(superblock_s);
    subgroup_h[0] = one;
    for (std::size_t idx = 1; idx < superblock_s; ++idx) {
      subgroup_h[idx] = field->Mul(subgroup_h[idx - 1], omega_big);
    }

    std::vector<storage_type> eval_raw_vec;
    eval_raw_vec.reserve(length);
    lrc_detail::RawElementIndexMap<storage_type> pt_to_idx;
    pt_to_idx.Reserve(length);

    storage_type coset_rep = one;
    for (std::size_t ell = 0; ell < num_superblocks; ++ell) {
      for (std::size_t idx = 0; idx < superblock_s; ++idx) {
        const storage_type pt = field->Mul(coset_rep, subgroup_h[idx]);
        pt_to_idx.Insert(pt, eval_raw_vec.size());
        eval_raw_vec.push_back(pt);
      }
      if (ell + 1 < num_superblocks) {
        coset_rep = field->Mul(coset_rep, gamma);
      }
    }

    // Build the t fiber partitions via multiplicative cosets alpha * H_u.
    std::vector<std::vector<std::vector<std::size_t>>> fiber_partitions(t);
    for (std::size_t u = 0; u < t; ++u) {
      const std::size_t s_u = block_sizes[u];
      const storage_type omega_u =
          field->Pow(gamma, group_order / static_cast<uint64_t>(s_u));
      std::vector<storage_type> h_u(s_u);
      h_u[0] = one;
      for (std::size_t idx = 1; idx < s_u; ++idx) {
        h_u[idx] = field->Mul(h_u[idx - 1], omega_u);
      }

      std::vector<bool> visited(length, false);
      fiber_partitions[u].reserve(length / s_u);
      for (std::size_t j = 0; j < length; ++j) {
        if (visited[j]) {
          continue;
        }
        std::vector<std::size_t> fiber;
        fiber.reserve(s_u);
        for (std::size_t idx = 0; idx < s_u; ++idx) {
          const storage_type member = field->Mul(eval_raw_vec[j], h_u[idx]);
          const auto pos = pt_to_idx.Find(member);
          if (!pos.has_value()) {
            throw std::runtime_error(
                "CreateMultiplicativeAvailability: coset element outside "
                "evaluation domain");
          }
          visited[*pos] = true;
          fiber.push_back(*pos);
        }
        std::sort(fiber.begin(), fiber.end());
        fiber_partitions[u].push_back(std::move(fiber));
      }
    }

    // Build k x n generator matrix G_{i, j} = alpha_j^{e_i}.
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray gen_raw = RawArray::from_shape(ShapeType{dimension, length});
    storage_type *gen_data = gen_raw.data();
    for (std::size_t i = 0; i < dimension; ++i) {
      const std::size_t exp = crt_exponents[i];
      storage_type *row = gen_data + i * length;
      for (std::size_t j = 0; j < length; ++j) {
        row[j] = field->Pow(eval_raw_vec[j], exp);
      }
    }
    matrix_type G = matrix_type::FromRaw(field, std::move(gen_raw));
    matrix_type H = Base::ComputeParityCheckFromGenerator(G);

    const std::size_t max_deg = crt_exponents.back();
    const std::size_t designed_d = length - max_deg;
    polynomial_type good_poly =
        polynomial_type::Monomial(field, block_sizes[0], one);

    // Local interpolation coordinates are the evaluation points themselves.
    std::vector<std::vector<storage_type>> local_evals(t, eval_raw_vec);

    return LocallyRecoverableCode(
        std::move(field), LRCPartitionType::kAvailabilityMultiplicative,
        std::move(eval_raw_vec), std::move(good_poly), localities,
        std::move(deltas), std::move(block_sizes), std::move(fiber_partitions),
        std::move(local_evals), std::move(G), std::move(H), designed_d,
        /*exact_min_distance=*/std::nullopt, max_deg,
        /*has_univariate_ambient_rs=*/true);
  }

  // ==========================================================================
  // Factory 3b: Multi-Set Additive Direct-Sum Availability LRC
  // ==========================================================================

  static LocallyRecoverableCode CreateAdditiveAvailability(
      std::shared_ptr<GaloisField> field, std::size_t num_superblocks,
      std::size_t message_superblocks,
      const std::vector<std::size_t> &localities,
      const std::vector<std::size_t> &local_distances = {}) {
    if (!field) {
      throw std::invalid_argument(
          "CreateAdditiveAvailability: field cannot be null");
    }
    const std::size_t t = localities.size();
    if (t < 2) {
      throw std::invalid_argument(
          "CreateAdditiveAvailability: availability t must be >= 2");
    }
    if (num_superblocks == 0) {
      throw std::invalid_argument(
          "CreateAdditiveAvailability: num_superblocks must be >= 1");
    }
    if (message_superblocks == 0 || message_superblocks > num_superblocks) {
      throw std::invalid_argument(
          "CreateAdditiveAvailability: message_superblocks must satisfy "
          "1 <= m_blk <= m_B");
    }

    std::vector<std::size_t> deltas = local_distances;
    if (deltas.empty()) {
      deltas.assign(t, std::size_t{2});
    } else if (deltas.size() != t) {
      throw std::invalid_argument(
          "CreateAdditiveAvailability: local_distances size must match "
          "localities size");
    }

    const uint64_t p = field->Characteristic();
    const std::size_t m = static_cast<std::size_t>(field->Degree());

    std::vector<std::size_t> block_sizes(t);
    std::vector<std::size_t> dims_a(t);
    std::size_t total_a = 0;
    std::size_t superblock_s = 1;
    std::size_t inner_k = 1;
    std::size_t prod_deltas = 1;

    for (std::size_t u = 0; u < t; ++u) {
      if (localities[u] == 0) {
        throw std::invalid_argument(
            "CreateAdditiveAvailability: every locality r_u must be >= 1");
      }
      if (deltas[u] < 2) {
        throw std::invalid_argument(
            "CreateAdditiveAvailability: every local distance delta_u must be "
            ">= 2");
      }
      block_sizes[u] = localities[u] + deltas[u] - 1;
      const auto exp_u = lrc_detail::ExactPrimePowerExponent(
          static_cast<uint64_t>(block_sizes[u]), p);
      if (!exp_u.has_value()) {
        throw std::invalid_argument(
            "CreateAdditiveAvailability: every block size s_u = r_u + delta_u "
            "- 1 must be a power of field characteristic p");
      }
      dims_a[u] = *exp_u;
      total_a += dims_a[u];
      superblock_s = lrc_detail::SafeMulSize(superblock_s, block_sizes[u],
                                             "CreateAdditiveAvailability");
      inner_k = lrc_detail::SafeMulSize(inner_k, localities[u],
                                        "CreateAdditiveAvailability");
      prod_deltas = lrc_detail::SafeMulSize(prod_deltas, deltas[u],
                                            "CreateAdditiveAvailability");
    }

    if (total_a > m) {
      throw std::invalid_argument(
          "CreateAdditiveAvailability: sum(a_u) = " + std::to_string(total_a) +
          " exceeds field extension degree m = " + std::to_string(m));
    }

    const std::size_t w_dim = m - total_a;
    if (w_dim < 64) {
      uint64_t max_superblocks = 0;
      try {
        max_superblocks = utils::SafeIntegerPower(p, w_dim);
      } catch (const std::overflow_error &) {
        max_superblocks = std::numeric_limits<uint64_t>::max();
      }
      if (static_cast<uint64_t>(num_superblocks) > max_superblocks) {
        throw std::invalid_argument(
            "CreateAdditiveAvailability: num_superblocks exceeds p^(m - "
            "sum(a_u))");
      }
    }

    const std::size_t length = lrc_detail::SafeMulSize(
        num_superblocks, superblock_s, "CreateAdditiveAvailability");
    const std::size_t dimension = lrc_detail::SafeMulSize(
        message_superblocks, inner_k, "CreateAdditiveAvailability");
    const std::size_t exact_d =
        lrc_detail::SafeMulSize(num_superblocks - message_superblocks + 1,
                                prod_deltas, "CreateAdditiveAvailability");

    // Decompose F_p-basis into H_0 (+) ... (+) H_{t-1} (+) W.
    const std::vector<storage_type> fp_basis =
        lrc_detail::ExtractPrimeFieldBasis(*field);
    std::vector<std::vector<storage_type>> subspaces(t);
    std::size_t offset = 0;
    for (std::size_t u = 0; u < t; ++u) {
      const std::span<const storage_type> sub_basis(fp_basis.data() + offset,
                                                    dims_a[u]);
      subspaces[u] =
          lrc_detail::EnumeratePrimeSpan(*field, sub_basis, block_sizes[u]);
      offset += dims_a[u];
    }
    const std::span<const storage_type> w_basis(fp_basis.data() + total_a,
                                                w_dim);
    const std::vector<storage_type> coset_reps =
        lrc_detail::EnumeratePrimeSpan(*field, w_basis, num_superblocks);

    // Linearized vanishing polynomial of H_0 (for GoodPolynomial() accessor)
    // and of full H = H_0 (+) ... (+) H_{t-1} (for outer superblock values).
    polynomial_type good_poly = lrc_detail::ComputeLinearizedVanishingPoly(
        field, std::span<const storage_type>(fp_basis.data(), dims_a[0]));
    polynomial_type superblock_poly =
        lrc_detail::ComputeLinearizedVanishingPoly(
            field, std::span<const storage_type>(fp_basis.data(), total_a));

    std::vector<storage_type> superblock_consts(num_superblocks);
    for (std::size_t ell = 0; ell < num_superblocks; ++ell) {
      superblock_consts[ell] = superblock_poly.EvalRaw(coset_reps[ell]);
    }

    // Strides within each superblock for mixed-radix indexing over
    // (idx_0, ..., idx_{t-1}) in prod_u [0, s_u).
    std::vector<std::size_t> strides(t, 1);
    for (std::size_t u = t; u-- > 1;) {
      strides[u - 1] = strides[u] * block_sizes[u];
    }

    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();

    std::vector<storage_type> eval_raw_vec(length, zero);
    std::vector<std::vector<storage_type>> local_evals(
        t, std::vector<storage_type>(length, zero));
    std::vector<std::vector<std::size_t>> coord_digits(
        length, std::vector<std::size_t>(t, 0));

    for (std::size_t ell = 0; ell < num_superblocks; ++ell) {
      for (std::size_t local_idx = 0; local_idx < superblock_s; ++local_idx) {
        const std::size_t j = ell * superblock_s + local_idx;
        std::size_t rem = local_idx;
        storage_type pt = coset_reps[ell];
        for (std::size_t u = 0; u < t; ++u) {
          const std::size_t dig = rem / strides[u];
          rem %= strides[u];
          coord_digits[j][u] = dig;
          local_evals[u][j] = subspaces[u][dig];
          pt = field->Add(pt, subspaces[u][dig]);
        }
        eval_raw_vec[j] = pt;
      }
    }

    // Build fiber partitions: for each u, fibers vary digit u over [0, s_u)
    // while holding ell and all other digits v != u fixed.
    std::vector<std::vector<std::vector<std::size_t>>> fiber_partitions(t);
    for (std::size_t u = 0; u < t; ++u) {
      const std::size_t s_u = block_sizes[u];
      const std::size_t stride_u = strides[u];
      fiber_partitions[u].reserve(length / s_u);
      for (std::size_t ell = 0; ell < num_superblocks; ++ell) {
        for (std::size_t local_idx = 0; local_idx < superblock_s; ++local_idx) {
          if (coord_digits[ell * superblock_s + local_idx][u] != 0) {
            continue;
          }
          std::vector<std::size_t> fiber(s_u);
          const std::size_t base = ell * superblock_s + local_idx;
          for (std::size_t dig = 0; dig < s_u; ++dig) {
            fiber[dig] = base + dig * stride_u;
          }
          fiber_partitions[u].push_back(std::move(fiber));
        }
      }
    }

    // Precompute tables of h_u^i for 0 <= i < r_u and c_ell^j for 0 <= j <
    // m_blk.
    std::vector<std::vector<std::vector<storage_type>>> h_pows(t);
    for (std::size_t u = 0; u < t; ++u) {
      h_pows[u].resize(block_sizes[u],
                       std::vector<storage_type>(localities[u], one));
      for (std::size_t dig = 0; dig < block_sizes[u]; ++dig) {
        storage_type pw = one;
        for (std::size_t i = 0; i < localities[u]; ++i) {
          h_pows[u][dig][i] = pw;
          if (i + 1 < localities[u]) {
            pw = field->Mul(pw, subspaces[u][dig]);
          }
        }
      }
    }

    std::vector<std::vector<storage_type>> c_pows(
        num_superblocks, std::vector<storage_type>(message_superblocks, one));
    for (std::size_t ell = 0; ell < num_superblocks; ++ell) {
      storage_type pw = one;
      for (std::size_t j = 0; j < message_superblocks; ++j) {
        c_pows[ell][j] = pw;
        if (j + 1 < message_superblocks) {
          pw = field->Mul(pw, superblock_consts[ell]);
        }
      }
    }

    // Build generator matrix G of shape (k, n):
    // row = (j, i_0, ..., i_{t-1}) evaluates g_H(alpha)^j * prod_u
    // pi_u(alpha)^{i_u}.
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray gen_raw = RawArray::from_shape(ShapeType{dimension, length});
    storage_type *gen_data = gen_raw.data();

    std::vector<std::size_t> inner_strides(t, 1);
    for (std::size_t u = t; u-- > 1;) {
      inner_strides[u - 1] = inner_strides[u] * localities[u];
    }

    std::size_t row_idx = 0;
    for (std::size_t j_blk = 0; j_blk < message_superblocks; ++j_blk) {
      for (std::size_t inner = 0; inner < inner_k; ++inner) {
        std::vector<std::size_t> degs(t, 0);
        std::size_t rem = inner;
        for (std::size_t u = 0; u < t; ++u) {
          degs[u] = rem / inner_strides[u];
          rem %= inner_strides[u];
        }
        storage_type *row = gen_data + row_idx * length;
        for (std::size_t ell = 0; ell < num_superblocks; ++ell) {
          const storage_type &c_pow = c_pows[ell][j_blk];
          for (std::size_t local_idx = 0; local_idx < superblock_s;
               ++local_idx) {
            const std::size_t col = ell * superblock_s + local_idx;
            storage_type val = c_pow;
            for (std::size_t u = 0; u < t; ++u) {
              val = field->Mul(val, h_pows[u][coord_digits[col][u]][degs[u]]);
            }
            row[col] = val;
          }
        }
        ++row_idx;
      }
    }

    matrix_type G = matrix_type::FromRaw(field, std::move(gen_raw));
    matrix_type H = Base::ComputeParityCheckFromGenerator(G);

    return LocallyRecoverableCode(
        std::move(field), LRCPartitionType::kAvailabilityAdditive,
        std::move(eval_raw_vec), std::move(good_poly), localities,
        std::move(deltas), std::move(block_sizes), std::move(fiber_partitions),
        std::move(local_evals), std::move(G), std::move(H), exact_d, exact_d,
        /*max_poly_degree=*/0, /*has_univariate_ambient_rs=*/false);
  }

  // ==========================================================================
  // Factory 3c: Hybrid Semidirect-Product Availability LRC (H_1 |><| H_2)
  // ==========================================================================

  static LocallyRecoverableCode CreateHybridAvailability(
      std::shared_ptr<GaloisField> field, std::size_t num_orbits,
      std::size_t subfield_degree_b, std::size_t additive_dim_a,
      std::size_t mult_order_s2, std::size_t r1, std::size_t r2,
      std::size_t max_poly_degree) {
    if (!field) {
      throw std::invalid_argument(
          "CreateHybridAvailability: field cannot be null");
    }
    const uint64_t q = field->Order();
    const uint64_t p = field->Characteristic();
    const std::size_t m = static_cast<std::size_t>(field->Degree());
    const uint64_t group_order = lrc_detail::MultiplicativeGroupOrder(q);

    if (subfield_degree_b == 0 || subfield_degree_b >= m ||
        (m % subfield_degree_b != 0)) {
      throw std::invalid_argument(
          "CreateHybridAvailability: subfield_degree_b must satisfy 1 <= b < "
          "m and b | m");
    }
    const std::size_t d0 = m / subfield_degree_b;
    if (additive_dim_a == 0 || additive_dim_a >= d0) {
      throw std::invalid_argument(
          "CreateHybridAvailability: additive_dim_a must satisfy 1 <= a < "
          "m / b");
    }

    const uint64_t q0 = utils::SafeIntegerPower(p, subfield_degree_b);
    if (mult_order_s2 < 2 ||
        (q0 - 1) % static_cast<uint64_t>(mult_order_s2) != 0) {
      throw std::invalid_argument(
          "CreateHybridAvailability: mult_order_s2 must be >= 2 and divide "
          "p^b - 1");
    }

    const uint64_t s1_u64 = utils::SafeIntegerPower(q0, additive_dim_a);
    const std::size_t s1 = static_cast<std::size_t>(s1_u64);
    const std::size_t s2 = mult_order_s2;

    if (r1 == 0 || r1 >= s1) {
      throw std::invalid_argument(
          "CreateHybridAvailability: r1 must satisfy 1 <= r1 < s1");
    }
    if (r2 == 0 || r2 >= s2) {
      throw std::invalid_argument(
          "CreateHybridAvailability: r2 must satisfy 1 <= r2 < s2");
    }
    const std::size_t delta1 = s1 - r1 + 1;
    const std::size_t delta2 = s2 - r2 + 1;

    if (num_orbits == 0) {
      throw std::invalid_argument(
          "CreateHybridAvailability: num_orbits must be >= 1");
    }
    const std::size_t rem_deg = m - subfield_degree_b * additive_dim_a;
    if (rem_deg < 64) {
      uint64_t max_orbits = 0;
      try {
        max_orbits = (utils::SafeIntegerPower(p, rem_deg) - 1) /
                     static_cast<uint64_t>(s2);
      } catch (const std::overflow_error &) {
        max_orbits = std::numeric_limits<uint64_t>::max();
      }
      if (static_cast<uint64_t>(num_orbits) > max_orbits) {
        throw std::invalid_argument(
            "CreateHybridAvailability: num_orbits exceeds available orbits "
            "(p^(m - a*b) - 1) / s2");
      }
    }

    const std::size_t orbit_size =
        lrc_detail::SafeMulSize(s1, s2, "CreateHybridAvailability");
    const std::size_t length = lrc_detail::SafeMulSize(
        num_orbits, orbit_size, "CreateHybridAvailability");

    if (max_poly_degree == 0 || max_poly_degree > length) {
      throw std::invalid_argument(
          "CreateHybridAvailability: max_poly_degree D must satisfy 1 <= D <= "
          "n");
    }

    // 1. Construct subfield K = F_{p^b} and multiplicative subgroup H_2 <= K^*.
    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();
    const storage_type gamma = field->MultiplicativeGenerator();
    const storage_type eta = field->Pow(gamma, group_order / (q0 - 1));

    std::vector<storage_type> subfield_k;
    subfield_k.reserve(static_cast<std::size_t>(q0));
    subfield_k.push_back(zero);
    storage_type cur_eta = one;
    for (uint64_t i = 0; i < q0 - 1; ++i) {
      subfield_k.push_back(cur_eta);
      if (i + 1 < q0 - 1) {
        cur_eta = field->Mul(cur_eta, eta);
      }
    }

    // F_p-basis of K: {1, eta, ..., eta^{b-1}}.
    std::vector<storage_type> k_fp_basis(subfield_degree_b);
    k_fp_basis[0] = one;
    for (std::size_t i = 1; i < subfield_degree_b; ++i) {
      k_fp_basis[i] = field->Mul(k_fp_basis[i - 1], eta);
    }

    const storage_type omega2 =
        field->Pow(gamma, group_order / static_cast<uint64_t>(s2));
    std::vector<storage_type> subgroup_h2(s2);
    subgroup_h2[0] = one;
    for (std::size_t i = 1; i < s2; ++i) {
      subgroup_h2[i] = field->Mul(subgroup_h2[i - 1], omega2);
    }

    // 2. Construct K-basis {v_0, ..., v_{d0-1}} of F_{p^m}, then H_1 and W_K.
    const std::vector<storage_type> k_basis =
        lrc_detail::FindSubfieldBasis(*field, k_fp_basis, subfield_degree_b);

    // F_p-basis of H_1 = {u_j * v_i : 0 <= i < a, 0 <= j < b}.
    std::vector<storage_type> h1_fp_basis;
    h1_fp_basis.reserve(additive_dim_a * subfield_degree_b);
    for (std::size_t i = 0; i < additive_dim_a; ++i) {
      for (std::size_t j = 0; j < subfield_degree_b; ++j) {
        h1_fp_basis.push_back(field->Mul(k_fp_basis[j], k_basis[i]));
      }
    }
    const std::vector<storage_type> subgroup_h1 =
        lrc_detail::EnumeratePrimeSpan(*field, h1_fp_basis, s1);

    polynomial_type good_poly =
        lrc_detail::ComputeLinearizedVanishingPoly(field, h1_fp_basis);

    // 3. Pick `num_orbits` distinct orbits of H_1 |><| H_2 on F_{p^m} \ H_1.
    const std::size_t w_k_dim = d0 - additive_dim_a;
    std::vector<std::size_t> digits(w_k_dim, 0);
    digits[0] = 1;  // Skip 0 in W_K.

    lrc_detail::RawElementIndexMap<storage_type> visited_w;
    visited_w.Reserve(num_orbits * s2);

    std::vector<storage_type> eval_raw_vec;
    eval_raw_vec.reserve(length);
    lrc_detail::RawElementIndexMap<storage_type> pt_to_idx;
    pt_to_idx.Reserve(length);

    std::size_t collected_orbits = 0;
    while (collected_orbits < num_orbits) {
      storage_type w = zero;
      for (std::size_t i = 0; i < w_k_dim; ++i) {
        if (digits[i] != 0) {
          w = field->Add(w, field->Mul(subfield_k[digits[i]],
                                       k_basis[additive_dim_a + i]));
        }
      }
      if (!visited_w.Contains(w)) {
        for (std::size_t i2 = 0; i2 < s2; ++i2) {
          const storage_type h2w = field->Mul(subgroup_h2[i2], w);
          visited_w.Insert(h2w, collected_orbits);
          for (std::size_t i1 = 0; i1 < s1; ++i1) {
            const storage_type pt = field->Add(h2w, subgroup_h1[i1]);
            pt_to_idx.Insert(pt, eval_raw_vec.size());
            eval_raw_vec.push_back(pt);
          }
        }
        ++collected_orbits;
      }

      // Increment base-q0 counter over W_K.
      bool carry = true;
      for (std::size_t i = 0; i < w_k_dim && carry; ++i) {
        ++digits[i];
        if (digits[i] < static_cast<std::size_t>(q0)) {
          carry = false;
        } else {
          digits[i] = 0;
        }
      }
      if (carry && collected_orbits < num_orbits) {
        throw std::runtime_error(
            "CreateHybridAvailability: exhausted W_K before finding requested "
            "orbits");
      }
    }

    // 4. Build additive fibers (alpha + H_1) and multiplicative fibers
    //    (alpha * H_2).
    std::vector<std::vector<std::vector<std::size_t>>> fiber_partitions(2);

    // Partition 0: additive cosets of H_1.
    {
      std::vector<bool> visited(length, false);
      fiber_partitions[0].reserve(length / s1);
      for (std::size_t j = 0; j < length; ++j) {
        if (visited[j]) continue;
        std::vector<std::size_t> fiber;
        fiber.reserve(s1);
        for (std::size_t i1 = 0; i1 < s1; ++i1) {
          const storage_type member =
              field->Add(eval_raw_vec[j], subgroup_h1[i1]);
          const auto pos = pt_to_idx.Find(member);
          if (!pos.has_value()) {
            throw std::runtime_error(
                "CreateHybridAvailability: additive coset left orbit");
          }
          visited[*pos] = true;
          fiber.push_back(*pos);
        }
        std::sort(fiber.begin(), fiber.end());
        fiber_partitions[0].push_back(std::move(fiber));
      }
    }

    // Partition 1: multiplicative cosets of H_2.
    {
      std::vector<bool> visited(length, false);
      fiber_partitions[1].reserve(length / s2);
      for (std::size_t j = 0; j < length; ++j) {
        if (visited[j]) continue;
        std::vector<std::size_t> fiber;
        fiber.reserve(s2);
        for (std::size_t i2 = 0; i2 < s2; ++i2) {
          const storage_type member =
              field->Mul(eval_raw_vec[j], subgroup_h2[i2]);
          const auto pos = pt_to_idx.Find(member);
          if (!pos.has_value()) {
            throw std::runtime_error(
                "CreateHybridAvailability: multiplicative coset left orbit");
          }
          visited[*pos] = true;
          fiber.push_back(*pos);
        }
        std::sort(fiber.begin(), fiber.end());
        fiber_partitions[1].push_back(std::move(fiber));
      }
    }

    // 5. Build local parity-check matrix H_local (all RS dual checks on both
    //    fiber partitions) and Vandermonde evaluation matrix V (n x D).
    const std::size_t num_checks = fiber_partitions[0].size() * (s1 - r1) +
                                   fiber_partitions[1].size() * (s2 - r2);

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray h_local_raw = RawArray::from_shape(ShapeType{num_checks, length});
    h_local_raw.fill(zero);
    storage_type *hl_data = h_local_raw.data();

    const std::array<std::size_t, 2> part_rs = {r1, r2};
    std::size_t check_row = 0;
    for (std::size_t u = 0; u < 2; ++u) {
      const std::size_t r_u = part_rs[u];
      for (const auto &fiber : fiber_partitions[u]) {
        const std::size_t s_u = fiber.size();
        std::vector<storage_type> xs(s_u);
        for (std::size_t idx = 0; idx < s_u; ++idx) {
          xs[idx] = eval_raw_vec[fiber[idx]];
        }
        const std::vector<storage_type> dual_w =
            lrc_detail::ComputeBarycentricNodeWeights(*field, xs);
        const std::size_t num_fiber_checks = s_u - r_u;
        for (std::size_t mu = 0; mu < num_fiber_checks; ++mu) {
          storage_type *row_ptr = hl_data + check_row * length;
          for (std::size_t idx = 0; idx < s_u; ++idx) {
            const storage_type pw = field->Pow(xs[idx], mu);
            row_ptr[fiber[idx]] = field->Mul(dual_w[idx], pw);
          }
          ++check_row;
        }
      }
    }

    // Build Vandermonde matrix V (n x D): V[j, d] = alpha_j^d.
    RawArray v_raw = RawArray::from_shape(ShapeType{length, max_poly_degree});
    storage_type *v_data = v_raw.data();
    for (std::size_t j = 0; j < length; ++j) {
      storage_type pw = one;
      storage_type *v_row = v_data + j * max_poly_degree;
      for (std::size_t d = 0; d < max_poly_degree; ++d) {
        v_row[d] = pw;
        if (d + 1 < max_poly_degree) {
          pw = field->Mul(pw, eval_raw_vec[j]);
        }
      }
    }

    matrix_type H_local = matrix_type::FromRaw(field, std::move(h_local_raw));
    matrix_type V = matrix_type::FromRaw(field, std::move(v_raw));
    matrix_type constraint_mat = linalg::Dot(H_local, V);
    matrix_type poly_basis = linalg::Nullspace(constraint_mat);

    const std::size_t dimension = poly_basis.Shape(0);
    if (dimension == 0) {
      throw std::invalid_argument(
          "CreateHybridAvailability: polynomial subspace F_q[x]_{<D} intersect "
          "ker(H_local) is trivial (dimension 0); increase max_poly_degree or "
          "r1/r2");
    }

    // Determine exact max polynomial degree D_max among basis polynomials.
    std::size_t actual_max_deg = 0;
    const storage_type *pb_data = poly_basis.RawData().data();
    for (std::size_t row = 0; row < dimension; ++row) {
      for (std::size_t d = max_poly_degree; d-- > 0;) {
        if (!field->IsZero(pb_data[row * max_poly_degree + d])) {
          actual_max_deg = std::max(actual_max_deg, d);
          break;
        }
      }
    }

    matrix_type G = linalg::Dot(poly_basis, linalg::Transpose(V));
    matrix_type H = Base::ComputeParityCheckFromGenerator(G);
    const std::size_t designed_d = length - actual_max_deg;

    std::vector<std::size_t> localities = {r1, r2};
    std::vector<std::size_t> deltas = {delta1, delta2};
    std::vector<std::size_t> block_sizes = {s1, s2};
    std::vector<std::vector<storage_type>> local_evals(2, eval_raw_vec);

    return LocallyRecoverableCode(
        std::move(field), LRCPartitionType::kAvailabilityHybrid,
        std::move(eval_raw_vec), std::move(good_poly), std::move(localities),
        std::move(deltas), std::move(block_sizes), std::move(fiber_partitions),
        std::move(local_evals), std::move(G), std::move(H), designed_d,
        /*exact_min_distance=*/std::nullopt, actual_max_deg,
        /*has_univariate_ambient_rs=*/true);
  }

  // ==========================================================================
  // Accessors
  // ==========================================================================

  std::size_t Locality(std::size_t repair_set_idx = 0) const {
    ValidateRepairSetIndex(repair_set_idx, "Locality");
    return localities_[repair_set_idx];
  }

  std::size_t LocalDistance(std::size_t repair_set_idx = 0) const {
    ValidateRepairSetIndex(repair_set_idx, "LocalDistance");
    return local_distances_[repair_set_idx];
  }

  std::size_t BlockSize(std::size_t repair_set_idx = 0) const {
    ValidateRepairSetIndex(repair_set_idx, "BlockSize");
    return block_sizes_[repair_set_idx];
  }

  std::size_t Availability() const noexcept { return availability_; }

  LRCPartitionType PartitionType() const noexcept { return partition_type_; }

  const vector_type &EvaluationPoints() const noexcept { return eval_points_; }

  const polynomial_type &GoodPolynomial() const noexcept { return good_poly_; }

  const std::vector<std::vector<std::size_t>> &RepairSets(
      std::size_t coord) const {
    ValidateCoord(coord, "RepairSets");
    return repair_sets_[coord];
  }

  const std::vector<std::size_t> &RepairSet(
      std::size_t coord, std::size_t repair_set_idx = 0) const {
    ValidateCoord(coord, "RepairSet");
    ValidateRepairSetIndex(repair_set_idx, "RepairSet");
    return repair_sets_[coord][repair_set_idx];
  }

  const std::vector<std::vector<std::size_t>> &FiberPartition(
      std::size_t repair_set_idx = 0) const {
    ValidateRepairSetIndex(repair_set_idx, "FiberPartition");
    return fiber_partitions_[repair_set_idx];
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  std::size_t MinimumDistance() const override {
    if (this->cached_min_distance_.has_value()) {
      return *this->cached_min_distance_;
    }
    return designed_distance_;
  }

  std::size_t DesignedDistance() const noexcept { return designed_distance_; }

  std::size_t MaxPolynomialDegree() const noexcept { return max_poly_degree_; }

  bool HasUnivariateAmbientRS() const noexcept {
    return has_univariate_ambient_rs_;
  }

  std::string ToString() const override {
    const uint64_t q = this->field_->Order();
    const std::string field_str = (q == 0) ? "2^64" : std::to_string(q);
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + ", " +
           std::to_string(MinimumDistance()) +
           "] LRC (r=" + std::to_string(localities_[0]) +
           ", delta=" + std::to_string(local_distances_[0]) +
           ", t=" + std::to_string(availability_) + ") over GF(" + field_str +
           ")";
  }

  // ==========================================================================
  // Local Single-Symbol Repair (strictly r_u symbols read)
  // ==========================================================================

  LocalRepairResult<GaloisField> LocalRepair(
      const codeword_type &received, std::size_t erased_index,
      std::size_t repair_set_idx = 0,
      const std::vector<std::size_t> &other_erased_positions = {}) const {
    ValidateReceivedWord(received, "LocalRepair");
    if (erased_index >= this->length_) {
      return LocalRepairResult<GaloisField>::Failure(
          erased_index, repair_set_idx,
          "LocalRepair: erased_index out of range");
    }
    if (repair_set_idx >= availability_) {
      return LocalRepairResult<GaloisField>::Failure(
          erased_index, repair_set_idx,
          "LocalRepair: repair_set_idx out of range");
    }

    const std::size_t r_u = localities_[repair_set_idx];
    const auto &rset = repair_sets_[erased_index][repair_set_idx];
    const auto recv_span = received.RawSpan();

    // Check if any position in other_erased_positions collides with the
    // canonical r_u helper indices rset[0..r_u-1].
    bool canonical_ok = true;
    std::vector<bool> other_erased;
    if (!other_erased_positions.empty()) {
      other_erased.assign(this->length_, false);
      for (std::size_t pos : other_erased_positions) {
        if (pos >= this->length_) {
          return LocalRepairResult<GaloisField>::Failure(
              erased_index, repair_set_idx,
              "LocalRepair: position in other_erased_positions out of range");
        }
        other_erased[pos] = true;
      }
      for (std::size_t a = 0; a < r_u; ++a) {
        if (other_erased[rset[a]]) {
          canonical_ok = false;
          break;
        }
      }
    }

    // Fast path: O(r_u) dot product using precomputed barycentric weights.
    if (canonical_ok) {
      const auto &weights = canonical_weights_[erased_index][repair_set_idx];
      storage_type acc = this->field_->AdditiveIdentity();
      for (std::size_t a = 0; a < r_u; ++a) {
        acc = this->field_->Add(
            acc, this->field_->Mul(weights[a], recv_span[rset[a]]));
      }
      return LocalRepairResult<GaloisField>::Success(
          element_type::FromInternal(acc, this->field_), erased_index,
          repair_set_idx, r_u);
    }

    // Dynamic path: pick strictly r_u surviving indices in the repair set and
    // run O(r_u^2) barycentric Lagrange interpolation with BatchInv.
    std::vector<std::size_t> helpers;
    helpers.reserve(r_u);
    for (std::size_t idx : rset) {
      if (!other_erased[idx]) {
        helpers.push_back(idx);
        if (helpers.size() == r_u) {
          break;
        }
      }
    }

    if (helpers.size() < r_u) {
      return LocalRepairResult<GaloisField>::Failure(
          erased_index, repair_set_idx,
          "LocalRepair: fewer than r_u (" + std::to_string(r_u) +
              ") surviving symbols in repair set");
    }

    const auto &coord_map = local_eval_points_[repair_set_idx];
    std::vector<storage_type> xs(r_u);
    for (std::size_t a = 0; a < r_u; ++a) {
      xs[a] = coord_map[helpers[a]];
    }
    const std::vector<storage_type> weights =
        lrc_detail::ComputeLagrangeEvalWeights(*this->field_, xs,
                                               coord_map[erased_index]);

    storage_type acc = this->field_->AdditiveIdentity();
    for (std::size_t a = 0; a < r_u; ++a) {
      acc = this->field_->Add(
          acc, this->field_->Mul(weights[a], recv_span[helpers[a]]));
    }
    return LocalRepairResult<GaloisField>::Success(
        element_type::FromInternal(acc, this->field_), erased_index,
        repair_set_idx, r_u);
  }

  // ==========================================================================
  // Multi-Erasure Peeling Cascade + Global Linear Solve Fallback
  // ==========================================================================

  MultiLocalRepairResult<GaloisField> LocalRepairMultiple(
      const codeword_type &received,
      const std::vector<std::size_t> &erasure_positions,
      bool allow_global_fallback = true) const {
    ValidateReceivedWord(received, "LocalRepairMultiple");

    const std::size_t n = this->length_;
    std::vector<bool> erased(n, false);
    for (std::size_t pos : erasure_positions) {
      if (pos >= n) {
        return MultiLocalRepairResult<GaloisField>::Failure(
            "LocalRepairMultiple: erasure position " + std::to_string(pos) +
            " out of range");
      }
      if (erased[pos]) {
        return MultiLocalRepairResult<GaloisField>::Failure(
            "LocalRepairMultiple: duplicate erasure position " +
            std::to_string(pos));
      }
      erased[pos] = true;
    }

    codeword_type working = received;
    std::size_t remaining = erasure_positions.size();
    std::size_t num_local_repairs = 0;
    std::size_t num_global_repairs = 0;
    std::size_t total_symbols_read = 0;
    bool used_global_fallback = false;

    // 1. Greedy peeling cascade across all availability partitions and fibers.
    //    Each repairable fiber reads strictly r_u surviving symbols ONCE and
    //    recovers ALL <= delta_u - 1 erasures inside that fiber simultaneously.
    while (remaining > 0) {
      bool progress = false;
      for (std::size_t u = 0; u < availability_ && remaining > 0; ++u) {
        const std::size_t r_u = localities_[u];
        const std::size_t max_fiber_erasures = local_distances_[u] - 1;
        const auto &coord_map = local_eval_points_[u];

        for (const auto &fiber : fiber_partitions_[u]) {
          std::vector<std::size_t> fiber_erased;
          std::vector<std::size_t> fiber_helpers;
          fiber_helpers.reserve(r_u);

          for (std::size_t idx : fiber) {
            if (erased[idx]) {
              fiber_erased.push_back(idx);
            } else if (fiber_helpers.size() < r_u) {
              fiber_helpers.push_back(idx);
            }
          }

          if (fiber_erased.empty() ||
              fiber_erased.size() > max_fiber_erasures ||
              fiber_helpers.size() < r_u) {
            continue;
          }

          // Read strictly r_u surviving symbols once for this fiber.
          total_symbols_read += r_u;

          if (fiber_erased.size() == 1) {
            const std::size_t target = fiber_erased[0];
            const auto &canon_rset = repair_sets_[target][u];
            bool matches_canon = true;
            for (std::size_t a = 0; a < r_u; ++a) {
              if (fiber_helpers[a] != canon_rset[a]) {
                matches_canon = false;
                break;
              }
            }
            storage_type acc = this->field_->AdditiveIdentity();
            if (matches_canon) {
              const auto &w = canonical_weights_[target][u];
              for (std::size_t a = 0; a < r_u; ++a) {
                acc = this->field_->Add(
                    acc,
                    this->field_->Mul(w[a], working.Raw(fiber_helpers[a])));
              }
            } else {
              std::vector<storage_type> xs(r_u);
              for (std::size_t a = 0; a < r_u; ++a) {
                xs[a] = coord_map[fiber_helpers[a]];
              }
              const std::vector<storage_type> w =
                  lrc_detail::ComputeLagrangeEvalWeights(*this->field_, xs,
                                                         coord_map[target]);
              for (std::size_t a = 0; a < r_u; ++a) {
                acc = this->field_->Add(
                    acc,
                    this->field_->Mul(w[a], working.Raw(fiber_helpers[a])));
              }
            }
            working.Raw(target) = acc;
            erased[target] = false;
            --remaining;
            ++num_local_repairs;
            progress = true;
          } else {
            // Multiple erasures in the same fiber: compute barycentric weights
            // once for the r_u surviving nodes and evaluate at all erased
            // targets in a single BatchInv.
            std::vector<storage_type> xs(r_u);
            std::vector<storage_type> ys(r_u);
            for (std::size_t a = 0; a < r_u; ++a) {
              xs[a] = coord_map[fiber_helpers[a]];
              ys[a] = working.Raw(fiber_helpers[a]);
            }
            const std::vector<storage_type> betas =
                lrc_detail::ComputeBarycentricNodeWeights(*this->field_, xs);
            std::vector<storage_type> y_betas(r_u);
            this->field_->MulArray(std::span<storage_type>(y_betas),
                                   std::span<const storage_type>(ys),
                                   std::span<const storage_type>(betas));

            const std::size_t num_targets = fiber_erased.size();
            std::vector<storage_type> diffs(num_targets * r_u);
            std::vector<storage_type> prods(
                num_targets, this->field_->MultiplicativeIdentity());
            for (std::size_t t_idx = 0; t_idx < num_targets; ++t_idx) {
              const storage_type &tx = coord_map[fiber_erased[t_idx]];
              storage_type prod = this->field_->MultiplicativeIdentity();
              for (std::size_t a = 0; a < r_u; ++a) {
                const storage_type diff = this->field_->Sub(tx, xs[a]);
                diffs[t_idx * r_u + a] = diff;
                prod = this->field_->Mul(prod, diff);
              }
              prods[t_idx] = prod;
            }
            std::vector<storage_type> inv_diffs(num_targets * r_u);
            this->field_->BatchInv(std::span<storage_type>(inv_diffs),
                                   std::span<const storage_type>(diffs));

            for (std::size_t t_idx = 0; t_idx < num_targets; ++t_idx) {
              const std::size_t target = fiber_erased[t_idx];
              storage_type sum = this->field_->AdditiveIdentity();
              for (std::size_t a = 0; a < r_u; ++a) {
                sum = this->field_->Add(
                    sum,
                    this->field_->Mul(y_betas[a], inv_diffs[t_idx * r_u + a]));
              }
              working.Raw(target) = this->field_->Mul(prods[t_idx], sum);
              erased[target] = false;
            }
            remaining -= num_targets;
            num_local_repairs += num_targets;
            progress = true;
          }
        }
      }
      if (!progress) {
        break;
      }
    }

    // 2. Global parity-check fallback: solve H_{:, E} c_E^T = -H_{:, ~E}
    // c_{~E}^T.
    if (remaining > 0) {
      if (!allow_global_fallback) {
        return MultiLocalRepairResult<GaloisField>::Failure(
            "LocalRepairMultiple: " + std::to_string(remaining) +
                " erasures remain after local peeling and global fallback is "
                "disabled",
            num_local_repairs, 0, total_symbols_read, false);
      }

      std::vector<std::size_t> rem_erased;
      std::vector<std::size_t> known;
      rem_erased.reserve(remaining);
      known.reserve(n - remaining);
      for (std::size_t j = 0; j < n; ++j) {
        if (erased[j]) {
          rem_erased.push_back(j);
        } else {
          known.push_back(j);
        }
      }

      const std::size_t num_rows = parity_check_matrix_.Shape(0);
      if (remaining > num_rows) {
        return MultiLocalRepairResult<GaloisField>::Failure(
            "LocalRepairMultiple: remaining erasures (" +
                std::to_string(remaining) + ") exceed code redundancy (" +
                std::to_string(num_rows) + ")",
            num_local_repairs, 0, total_symbols_read, true);
      }

      used_global_fallback = true;
      total_symbols_read += known.size();

      const storage_type *h_data = parity_check_matrix_.RawData().data();
      const storage_type zero = this->field_->AdditiveIdentity();

      using RawArray = typename matrix_type::raw_array_type;
      using ShapeType = typename matrix_type::shape_type;

      if (num_rows == remaining) {
        // Square system: solve directly via linalg::Solve.
        RawArray h_e_raw =
            RawArray::from_shape(ShapeType{remaining, remaining});
        RawArray rhs_raw = RawArray::from_shape(ShapeType{remaining});
        storage_type *he_data = h_e_raw.data();
        storage_type *rhs_data = rhs_raw.data();

        for (std::size_t r = 0; r < num_rows; ++r) {
          const storage_type *h_row = h_data + r * n;
          storage_type acc = zero;
          for (std::size_t j : known) {
            if (!this->field_->IsZero(h_row[j]) &&
                !this->field_->IsZero(working.Raw(j))) {
              acc = this->field_->Add(
                  acc, this->field_->Mul(h_row[j], working.Raw(j)));
            }
          }
          rhs_data[r] = this->field_->Neg(acc);
          for (std::size_t c = 0; c < remaining; ++c) {
            he_data[r * remaining + c] = h_row[rem_erased[c]];
          }
        }

        matrix_type H_E =
            matrix_type::FromRaw(this->field_, std::move(h_e_raw));
        codeword_type rhs =
            codeword_type::FromRaw(this->field_, std::move(rhs_raw));
        try {
          codeword_type sol = linalg::Solve(H_E, rhs);
          for (std::size_t c = 0; c < remaining; ++c) {
            working.Raw(rem_erased[c]) = sol.Raw(c);
            erased[rem_erased[c]] = false;
          }
        } catch (const std::runtime_error &) {
          return MultiLocalRepairResult<GaloisField>::Failure(
              "LocalRepairMultiple: erasure pattern is globally unrecoverable "
              "(singular submatrix H_{:, E})",
              num_local_repairs, 0, total_symbols_read, true);
        }
      } else {
        // Overdetermined system (num_rows > remaining): use RrefWithPivots on
        // the augmented matrix [H_{:, E} | -H_{:, ~E} c_{~E}^T].
        const std::size_t aug_cols = remaining + 1;
        RawArray aug_raw = RawArray::from_shape(ShapeType{num_rows, aug_cols});
        storage_type *aug_data = aug_raw.data();

        for (std::size_t r = 0; r < num_rows; ++r) {
          const storage_type *h_row = h_data + r * n;
          storage_type *a_row = aug_data + r * aug_cols;
          storage_type acc = zero;
          for (std::size_t j : known) {
            if (!this->field_->IsZero(h_row[j]) &&
                !this->field_->IsZero(working.Raw(j))) {
              acc = this->field_->Add(
                  acc, this->field_->Mul(h_row[j], working.Raw(j)));
            }
          }
          for (std::size_t c = 0; c < remaining; ++c) {
            a_row[c] = h_row[rem_erased[c]];
          }
          a_row[remaining] = this->field_->Neg(acc);
        }

        matrix_type aug =
            matrix_type::FromRaw(this->field_, std::move(aug_raw));
        auto [rref, pivots] = linalg::detail::RrefWithPivots(aug);

        for (std::size_t pcol : pivots) {
          if (pcol == remaining) {
            return MultiLocalRepairResult<GaloisField>::Failure(
                "LocalRepairMultiple: inconsistent global parity-check system "
                "(errors present in surviving symbols)",
                num_local_repairs, 0, total_symbols_read, true);
          }
        }
        if (pivots.size() < remaining) {
          return MultiLocalRepairResult<GaloisField>::Failure(
              "LocalRepairMultiple: erasure pattern is globally unrecoverable "
              "(H_{:, E} is rank-deficient)",
              num_local_repairs, 0, total_symbols_read, true);
        }

        for (std::size_t c = 0; c < remaining; ++c) {
          working.Raw(rem_erased[c]) = rref.Raw(c, remaining);
          erased[rem_erased[c]] = false;
        }
      }

      num_global_repairs = remaining;
      remaining = 0;
    }

    // 3. Final parity-check verification & message extraction.
    if (!this->Contains(working)) {
      return MultiLocalRepairResult<GaloisField>::Failure(
          "LocalRepairMultiple: repaired word failed parity-check verification",
          num_local_repairs, num_global_repairs, total_symbols_read,
          used_global_fallback);
    }

    MultiLocalRepairResult<GaloisField> res;
    res.message = this->Unencode(working);
    res.codeword = std::move(working);
    res.num_local_repairs = num_local_repairs;
    res.num_global_repairs = num_global_repairs;
    res.total_symbols_read = total_symbols_read;
    res.used_global_fallback = used_global_fallback;
    res.success = true;
    return res;
  }

 private:
  static LocallyRecoverableCode BuildSingleSetCode(
      std::shared_ptr<GaloisField> field, LRCPartitionType ptype,
      std::vector<storage_type> eval_raw_vec, polynomial_type good_poly,
      std::vector<std::vector<std::size_t>> fibers,
      const std::vector<storage_type> &fiber_constants, std::size_t dimension,
      std::size_t locality, std::size_t local_distance) {
    const std::size_t length = eval_raw_vec.size();
    const std::size_t s = locality + local_distance - 1;
    const std::size_t q_k = dimension / locality;
    const std::size_t r_k = dimension % locality;
    const std::size_t ceil_k_over_r = lrc_detail::CeilDiv(dimension, locality);
    const std::size_t exact_d =
        length - dimension + 1 - (ceil_k_over_r - 1) * (local_distance - 1);
    const std::size_t max_deg = length - exact_d;

    // Build message basis polynomials in ascending-degree order:
    //   for j = 0..q_k - 1, i = 0..r - 1: g(x)^j * x^i
    //   if r_k > 0: for j = q_k, i = 0..r_k - 1: g(x)^{q_k} * x^i
    // Since g(alpha) = c_ell is constant on fiber ell, g(alpha)^j * alpha^i
    // = c_ell^j * alpha^i.
    const storage_type one = field->MultiplicativeIdentity();
    std::vector<storage_type> coord_const(length, one);
    for (std::size_t ell = 0; ell < fibers.size(); ++ell) {
      for (std::size_t idx : fibers[ell]) {
        coord_const[idx] = fiber_constants[ell];
      }
    }

    // Precompute alpha_j^i (0 <= i < r) and c_j^j_pow (0 <= j_pow < ceil(k/r)).
    std::vector<std::vector<storage_type>> x_pows(
        length, std::vector<storage_type>(locality, one));
    std::vector<std::vector<storage_type>> g_pows(
        length, std::vector<storage_type>(ceil_k_over_r, one));
    for (std::size_t col = 0; col < length; ++col) {
      storage_type xp = one;
      for (std::size_t i = 0; i < locality; ++i) {
        x_pows[col][i] = xp;
        if (i + 1 < locality) {
          xp = field->Mul(xp, eval_raw_vec[col]);
        }
      }
      storage_type gp = one;
      for (std::size_t j = 0; j < ceil_k_over_r; ++j) {
        g_pows[col][j] = gp;
        if (j + 1 < ceil_k_over_r) {
          gp = field->Mul(gp, coord_const[col]);
        }
      }
    }

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray gen_raw = RawArray::from_shape(ShapeType{dimension, length});
    storage_type *gen_data = gen_raw.data();

    std::size_t row_idx = 0;
    for (std::size_t j = 0; j < q_k; ++j) {
      for (std::size_t i = 0; i < locality; ++i) {
        storage_type *row = gen_data + row_idx * length;
        for (std::size_t col = 0; col < length; ++col) {
          row[col] = field->Mul(g_pows[col][j], x_pows[col][i]);
        }
        ++row_idx;
      }
    }
    if (r_k > 0) {
      for (std::size_t i = 0; i < r_k; ++i) {
        storage_type *row = gen_data + row_idx * length;
        for (std::size_t col = 0; col < length; ++col) {
          row[col] = field->Mul(g_pows[col][q_k], x_pows[col][i]);
        }
        ++row_idx;
      }
    }

    matrix_type G = matrix_type::FromRaw(field, std::move(gen_raw));
    matrix_type H = Base::ComputeParityCheckFromGenerator(G);

    std::vector<std::size_t> localities = {locality};
    std::vector<std::size_t> deltas = {local_distance};
    std::vector<std::size_t> block_sizes = {s};
    std::vector<std::vector<std::vector<std::size_t>>> fiber_partitions = {
        std::move(fibers)};
    std::vector<std::vector<storage_type>> local_evals = {eval_raw_vec};

    return LocallyRecoverableCode(
        std::move(field), ptype, std::move(eval_raw_vec), std::move(good_poly),
        std::move(localities), std::move(deltas), std::move(block_sizes),
        std::move(fiber_partitions), std::move(local_evals), std::move(G),
        std::move(H), exact_d, exact_d, max_deg,
        /*has_univariate_ambient_rs=*/true);
  }

  LocallyRecoverableCode(
      std::shared_ptr<GaloisField> field, LRCPartitionType partition_type,
      std::vector<storage_type> eval_raw_vec, polynomial_type good_poly,
      std::vector<std::size_t> localities,
      std::vector<std::size_t> local_distances,
      std::vector<std::size_t> block_sizes,
      std::vector<std::vector<std::vector<std::size_t>>> fiber_partitions,
      std::vector<std::vector<storage_type>> local_eval_points,
      matrix_type generator_matrix, matrix_type parity_check_matrix,
      std::size_t designed_distance,
      std::optional<std::size_t> exact_min_distance,
      std::size_t max_poly_degree, bool has_univariate_ambient_rs)
      : Base(std::move(field), generator_matrix.Shape(1),
             generator_matrix.Shape(0), Metric::HAMMING, "GeneratorMatrix",
             "LRCGlobal"),
        partition_type_(partition_type),
        localities_(std::move(localities)),
        local_distances_(std::move(local_distances)),
        block_sizes_(std::move(block_sizes)),
        availability_(localities_.size()),
        designed_distance_(designed_distance),
        max_poly_degree_(max_poly_degree),
        has_univariate_ambient_rs_(has_univariate_ambient_rs),
        good_poly_(std::move(good_poly)),
        fiber_partitions_(std::move(fiber_partitions)),
        local_eval_points_(std::move(local_eval_points)),
        generator_matrix_(std::move(generator_matrix)),
        parity_check_matrix_(std::move(parity_check_matrix)) {
    if (exact_min_distance.has_value()) {
      this->cached_min_distance_ = exact_min_distance;
    }

    using RawArray = typename vector_type::raw_array_type;
    using ShapeType = typename vector_type::shape_type;
    RawArray eval_arr = RawArray::from_shape(ShapeType{this->length_});
    std::copy(eval_raw_vec.begin(), eval_raw_vec.end(), eval_arr.data());
    eval_points_ = vector_type::FromRaw(this->field_, std::move(eval_arr));

    PrecomputeRepairStructures();
    RegisterEncodersAndDecoders();
  }

  void PrecomputeRepairStructures() {
    const std::size_t n = this->length_;
    const std::size_t t = availability_;

    coord_to_fiber_.assign(n, std::vector<std::size_t>(t, 0));
    repair_sets_.assign(n, std::vector<std::vector<std::size_t>>(t));
    canonical_weights_.assign(n, std::vector<std::vector<storage_type>>(t));

    for (std::size_t u = 0; u < t; ++u) {
      const std::size_t r_u = localities_[u];
      const auto &coord_map = local_eval_points_[u];

      for (std::size_t f_idx = 0; f_idx < fiber_partitions_[u].size();
           ++f_idx) {
        const auto &fiber = fiber_partitions_[u][f_idx];
        for (std::size_t target : fiber) {
          coord_to_fiber_[target][u] = f_idx;
          auto &rset = repair_sets_[target][u];
          rset.reserve(fiber.size() - 1);
          for (std::size_t other : fiber) {
            if (other != target) {
              rset.push_back(other);
            }
          }

          // Precompute barycentric Lagrange weights for canonical r_u helpers.
          std::vector<storage_type> xs(r_u);
          for (std::size_t a = 0; a < r_u; ++a) {
            xs[a] = coord_map[rset[a]];
          }
          canonical_weights_[target][u] =
              lrc_detail::ComputeLagrangeEvalWeights(*this->field_, xs,
                                                     coord_map[target]);
        }
      }
    }
  }

  void ValidateCoord(std::size_t coord, const char *fn_name) const {
    if (coord >= this->length_) {
      throw std::out_of_range(std::string(fn_name) +
                              ": coordinate index out of range");
    }
  }

  void ValidateRepairSetIndex(std::size_t idx, const char *fn_name) const {
    if (idx >= availability_) {
      throw std::out_of_range(std::string(fn_name) +
                              ": repair set index out of range");
    }
  }

  void ValidateReceivedWord(const codeword_type &received,
                            const char *fn_name) const {
    if (!received.HasField() || received.Field() != this->field_) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": received word field mismatch");
    }
    if (received.Dimension() != 1 || received.Size() != this->length_) {
      throw std::invalid_argument(
          std::string(fn_name) +
          ": received word must be 1D with size equal to code length");
    }
  }

  void RegisterEncodersAndDecoders();

  LRCPartitionType partition_type_;
  std::vector<std::size_t> localities_;
  std::vector<std::size_t> local_distances_;
  std::vector<std::size_t> block_sizes_;
  std::size_t availability_ = 1;
  std::size_t designed_distance_ = 1;
  std::size_t max_poly_degree_ = 0;
  bool has_univariate_ambient_rs_ = true;

  vector_type eval_points_;
  polynomial_type good_poly_;
  std::vector<std::vector<std::vector<std::size_t>>> fiber_partitions_;
  std::vector<std::vector<std::size_t>> coord_to_fiber_;
  std::vector<std::vector<std::vector<std::size_t>>> repair_sets_;
  std::vector<std::vector<storage_type>> local_eval_points_;
  std::vector<std::vector<std::vector<storage_type>>> canonical_weights_;

  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;
};

// ============================================================================
// LRCGlobalDecoder<GaloisField>
// ============================================================================

template <typename GaloisField>
class LRCGlobalDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;

  explicit LRCGlobalDecoder(const AbstractCode<GaloisField> &code)
      : Base(code) {
    lrc_code_ =
        dynamic_cast<const LocallyRecoverableCode<GaloisField> *>(&code);
    if (!lrc_code_) {
      throw std::invalid_argument(
          "LRCGlobalDecoder requires a LocallyRecoverableCode");
    }

    if (lrc_code_->HasUnivariateAmbientRS()) {
      const std::size_t ambient_dim = lrc_code_->MaxPolynomialDegree() + 1;
      ambient_grs_ = std::make_unique<GeneralizedReedSolomonCode<GaloisField>>(
          GeneralizedReedSolomonCode<GaloisField>::StandardReedSolomon(
              lrc_code_->Field(), lrc_code_->EvaluationPoints(), ambient_dim));
      ambient_decoder_ =
          std::make_unique<BerlekampMasseyDecoder<GaloisField>>(*ambient_grs_);
      decoding_radius_ = ambient_decoder_->DecodingRadius();
    } else {
      decoding_radius_ = 0;
    }
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return DecodeWithErasures(received_word, {});
  }

  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const override {
    if (!received_word.HasField() ||
        received_word.Field() != lrc_code_->Field()) {
      throw std::invalid_argument("LRCGlobal: received word field mismatch");
    }
    if (received_word.Dimension() != 1 ||
        received_word.Size() != lrc_code_->Length()) {
      throw std::invalid_argument(
          "LRCGlobal: received word must be 1D with size equal to code length");
    }

    if (ambient_decoder_) {
      auto grs_res = ambient_decoder_->DecodeWithErasures(received_word,
                                                          erasure_positions);
      if (!grs_res.Ok()) {
        return result_type::Failure("LRCGlobal: " + grs_res.error_message);
      }
      if (!lrc_code_->Contains(grs_res.codeword)) {
        return result_type::Failure(
            "LRCGlobal: decoded ambient RS word is not in the LRC subcode");
      }
      auto msg = lrc_code_->Unencode(grs_res.codeword);
      return result_type::Success(std::move(grs_res.codeword), std::move(msg),
                                  grs_res.num_errors, grs_res.num_erasures);
    }

    // Non-univariate LRC (e.g., additive direct-sum tensor-product code):
    // handle pure erasures and error-free codewords via LocalRepairMultiple.
    if (erasure_positions.empty()) {
      if (lrc_code_->Contains(received_word)) {
        return result_type::Success(received_word,
                                    lrc_code_->Unencode(received_word), 0, 0);
      }
      return result_type::Failure(
          "LRCGlobal: error correction requires a univariate-evaluation LRC");
    }

    auto rep = lrc_code_->LocalRepairMultiple(received_word, erasure_positions,
                                              /*allow_global_fallback=*/true);
    if (!rep.Ok()) {
      return result_type::Failure("LRCGlobal: " + rep.error_message);
    }
    return result_type::Success(std::move(rep.codeword), std::move(rep.message),
                                0, erasure_positions.size());
  }

  std::string Name() const override { return "LRCGlobal"; }

  std::size_t DecodingRadius() const override { return decoding_radius_; }

 private:
  const LocallyRecoverableCode<GaloisField> *lrc_code_;
  std::unique_ptr<GeneralizedReedSolomonCode<GaloisField>> ambient_grs_;
  std::unique_ptr<BerlekampMasseyDecoder<GaloisField>> ambient_decoder_;
  std::size_t decoding_radius_ = 0;
};

// ============================================================================
// LRCLocalThenGlobalDecoder<GaloisField>
// ============================================================================

template <typename GaloisField>
class LRCLocalThenGlobalDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;

  explicit LRCLocalThenGlobalDecoder(const AbstractCode<GaloisField> &code)
      : Base(code), global_decoder_(code) {
    lrc_code_ =
        dynamic_cast<const LocallyRecoverableCode<GaloisField> *>(&code);
    if (!lrc_code_) {
      throw std::invalid_argument(
          "LRCLocalThenGlobalDecoder requires a LocallyRecoverableCode");
    }
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return DecodeWithErasures(received_word, {});
  }

  result_type DecodeWithErasures(
      const codeword_type &received_word,
      const std::vector<std::size_t> &erasure_positions) const override {
    if (!received_word.HasField() ||
        received_word.Field() != lrc_code_->Field()) {
      throw std::invalid_argument(
          "LRCLocalThenGlobal: received word field mismatch");
    }
    if (received_word.Dimension() != 1 ||
        received_word.Size() != lrc_code_->Length()) {
      throw std::invalid_argument(
          "LRCLocalThenGlobal: received word must be 1D with size equal to "
          "code length");
    }

    if (erasure_positions.empty()) {
      if (lrc_code_->Contains(received_word)) {
        return result_type::Success(received_word,
                                    lrc_code_->Unencode(received_word), 0, 0);
      }
      return global_decoder_.DecodeWithResult(received_word);
    }

    // Attempt local peeling cascade (+ global pure-erasure solve) first.
    auto rep = lrc_code_->LocalRepairMultiple(received_word, erasure_positions,
                                              /*allow_global_fallback=*/true);
    if (rep.Ok()) {
      return result_type::Success(std::move(rep.codeword),
                                  std::move(rep.message), 0,
                                  erasure_positions.size());
    }

    // Fall back to ambient GRS errors-and-erasures decoding (handles
    // simultaneous errors + erasures up to 2e + s <= d - 1).
    return global_decoder_.DecodeWithErasures(received_word, erasure_positions);
  }

  std::string Name() const override { return "LRCLocalThenGlobal"; }

  std::size_t DecodingRadius() const override {
    return global_decoder_.DecodingRadius();
  }

 private:
  const LocallyRecoverableCode<GaloisField> *lrc_code_;
  LRCGlobalDecoder<GaloisField> global_decoder_;
};

template <typename GaloisField>
inline void LocallyRecoverableCode<GaloisField>::RegisterEncodersAndDecoders() {
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<GaloisField> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(lin);
      });
  this->RegisterDecoder("LRCGlobal", [](const AbstractCode<GaloisField> &code) {
    return std::make_shared<LRCGlobalDecoder<GaloisField>>(code);
  });
  this->RegisterDecoder(
      "LRCLocalThenGlobal", [](const AbstractCode<GaloisField> &code) {
        return std::make_shared<LRCLocalThenGlobalDecoder<GaloisField>>(code);
      });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_LRC_HPP_
