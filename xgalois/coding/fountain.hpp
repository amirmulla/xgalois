#ifndef XGALOIS_CODING_FOUNTAIN_HPP_
#define XGALOIS_CODING_FOUNTAIN_HPP_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

// Encoded fountain symbol (either LT or RaptorQ) carrying its deterministic
// Encoding Symbol ID (ESI), participating neighbor indices, linear combination
// coefficients, and payload vector of length symbol_size.
template <typename GaloisField>
struct FountainSymbol {
  uint32_t esi = 0;
  std::vector<std::size_t> neighbors;
  std::vector<typename GaloisField::element_type> coeffs;
  GFArray<GaloisField> payload;
};

// Result of decoding a fountain-coded source block (shape: K x symbol_size).
template <typename GaloisField>
struct FountainDecodingResult {
  GFMatrix<GaloisField> source_block;
  std::size_t symbols_used = 0;
  std::size_t num_inactivations = 0;
  bool used_gaussian_fallback = false;
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }
};

namespace detail {

// Fast, portable SplitMix64 PRNG for deterministic ESI-seeded fountain sampling
// across all platforms and standard libraries.
class FountainPrng {
 public:
  explicit FountainPrng(uint64_t seed) noexcept : state_(seed) {}

  uint64_t NextU64() noexcept {
    uint64_t z = (state_ += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

  // Uniform double in [0, 1).
  double NextDouble() noexcept {
    return static_cast<double>(NextU64() >> 11) * (1.0 / 9007199254740992.0);
  }

  // Uniform integer in [0, bound - 1] (bound > 0).
  std::size_t NextIndex(std::size_t bound) noexcept {
    return static_cast<std::size_t>(NextU64() % static_cast<uint64_t>(bound));
  }

 private:
  uint64_t state_;
};

inline uint64_t MixSeed(uint64_t base_seed, uint64_t salt,
                        uint64_t tag) noexcept {
  uint64_t x = base_seed ^ (salt * 0x9E3779B97F4A7C15ULL) ^
               (tag * 0xBF58476D1CE4E5B9ULL);
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

// Samples a uniform nonzero field element using FountainPrng.
template <typename GaloisField>
inline typename GaloisField::element_type SampleNonZeroElement(
    const GaloisField &field, FountainPrng &rng) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t order = field.Order();
  if (order == 2) {
    return field.MultiplicativeIdentity();
  }
  if (order > 2) {
    const uint64_t idx = 1ULL + (rng.NextU64() % (order - 1ULL));
    if constexpr (std::is_integral_v<StorageType>) {
      StorageType val = field.SetElementValue(static_cast<StorageType>(idx));
      return field.IsZero(val) ? field.MultiplicativeIdentity() : val;
    } else {
      StorageType val = field.FromInt(idx);
      return field.IsZero(val) ? field.MultiplicativeIdentity() : val;
    }
  }
  // Fallback for huge fields with Order() == 0 (e.g. GF(2^64)):
  std::mt19937_64 mt(rng.NextU64());
  StorageType val = field.Random(mt);
  while (field.IsZero(val)) {
    val = field.Random(mt);
  }
  return val;
}

// Samples `count` distinct sorted indices from [offset, offset + range).
inline std::vector<std::size_t> SampleDistinctIndices(std::size_t range,
                                                      std::size_t count,
                                                      std::size_t offset,
                                                      FountainPrng &rng) {
  count = std::min(count, range);
  std::vector<std::size_t> chosen;
  chosen.reserve(count);
  if (count * 3 >= range) {
    // Partial Fisher-Yates for dense selection.
    std::vector<std::size_t> pool(range);
    std::iota(pool.begin(), pool.end(), offset);
    for (std::size_t i = 0; i < count; ++i) {
      const std::size_t j = i + rng.NextIndex(range - i);
      std::swap(pool[i], pool[j]);
      chosen.push_back(pool[i]);
    }
  } else {
    // Rejection sampling for sparse selection.
    while (chosen.size() < count) {
      const std::size_t cand = offset + rng.NextIndex(range);
      if (std::find(chosen.begin(), chosen.end(), cand) == chosen.end()) {
        chosen.push_back(cand);
      }
    }
  }
  std::sort(chosen.begin(), chosen.end());
  return chosen;
}

inline bool IsSmallPrime(std::size_t n) noexcept {
  if (n < 2) {
    return false;
  }
  if (n == 2 || n == 3) {
    return true;
  }
  if ((n % 2) == 0 || (n % 3) == 0) {
    return false;
  }
  for (std::size_t d = 5; d * d <= n; d += 6) {
    if ((n % d) == 0 || (n % (d + 2)) == 0) {
      return false;
    }
  }
  return true;
}

inline std::size_t NextPrimeAtLeast(std::size_t n) noexcept {
  if (n <= 2) {
    return 2;
  }
  if ((n % 2) == 0) {
    ++n;
  }
  while (!IsSmallPrime(n)) {
    n += 2;
  }
  return n;
}

// Solves an overdetermined (m x n, m >= n) dense linear system A * X = B
// over GaloisField using row-echelon Gaussian elimination with partial pivoting
// and field->FmaArray. Returns true iff rank(A) == n (unique solution written
// into x_out of shape n x rhs_cols).
template <typename GaloisField>
inline bool DenseGaussianSolveOverdetermined(
    const GaloisField &field, std::size_t m, std::size_t n,
    std::size_t rhs_cols,
    std::vector<typename GaloisField::element_type> &a_data,
    std::vector<typename GaloisField::element_type> &b_data,
    std::vector<typename GaloisField::element_type> &x_out) {
  using StorageType = typename GaloisField::element_type;
  const StorageType zero = field.AdditiveIdentity();
  const StorageType one = field.MultiplicativeIdentity();

  if (n == 0) {
    x_out.clear();
    return true;
  }
  if (m < n) {
    return false;
  }

  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot_row = col;
    while (pivot_row < m && field.IsZero(a_data[pivot_row * n + col])) {
      ++pivot_row;
    }
    if (pivot_row == m) {
      return false;
    }
    if (pivot_row != col) {
      std::swap_ranges(a_data.data() + col * n + col,
                       a_data.data() + (col + 1) * n,
                       a_data.data() + pivot_row * n + col);
      if (rhs_cols > 0) {
        std::swap_ranges(b_data.data() + col * rhs_cols,
                         b_data.data() + (col + 1) * rhs_cols,
                         b_data.data() + pivot_row * rhs_cols);
      }
    }

    StorageType *a_piv = a_data.data() + col * n;
    StorageType *b_piv = b_data.data() + col * rhs_cols;
    const StorageType piv_val = a_piv[col];
    if (piv_val != one) {
      const StorageType inv_piv = field.Inv(piv_val);
      a_piv[col] = one;
      if (col + 1 < n) {
        std::span<StorageType> a_tail(a_piv + col + 1, n - col - 1);
        field.ScalarMulArray(a_tail, a_tail, inv_piv);
      }
      if (rhs_cols > 0) {
        std::span<StorageType> b_span(b_piv, rhs_cols);
        field.ScalarMulArray(b_span, b_span, inv_piv);
      }
    }

    const std::size_t tail_len = n - col - 1;
    std::span<const StorageType> a_tail(a_piv + col + 1, tail_len);
    std::span<const StorageType> b_src(b_piv, rhs_cols);
    for (std::size_t r = col + 1; r < m; ++r) {
      StorageType *a_row = a_data.data() + r * n;
      const StorageType factor = a_row[col];
      if (field.IsZero(factor)) {
        continue;
      }
      a_row[col] = zero;
      const StorageType neg_factor = field.Neg(factor);
      if (tail_len > 0) {
        field.FmaArray(std::span<StorageType>(a_row + col + 1, tail_len),
                       a_tail, neg_factor);
      }
      if (rhs_cols > 0) {
        field.FmaArray(
            std::span<StorageType>(b_data.data() + r * rhs_cols, rhs_cols),
            b_src, neg_factor);
      }
    }
  }

  // Back-substitution on the top n x n unit upper-triangular block.
  for (std::size_t col = n; col-- > 0;) {
    if (rhs_cols == 0) {
      break;
    }
    std::span<const StorageType> b_src(b_data.data() + col * rhs_cols,
                                       rhs_cols);
    for (std::size_t r = 0; r < col; ++r) {
      const StorageType factor = a_data[r * n + col];
      if (field.IsZero(factor)) {
        continue;
      }
      const StorageType neg_factor = field.Neg(factor);
      field.FmaArray(
          std::span<StorageType>(b_data.data() + r * rhs_cols, rhs_cols), b_src,
          neg_factor);
    }
  }

  x_out.assign(b_data.begin(),
               b_data.begin() + static_cast<std::ptrdiff_t>(n * rhs_cols));
  return true;
}

}  // namespace detail

// ============================================================================
// LTCode: Luby Transform Rateless Fountain Code over Any Finite Field
// ============================================================================
//
// Supports Ideal Soliton and exact Robust Soliton degree distributions,
// deterministic ESI-seeded symbol generation, O(K ln(K/delta) * symbol_size)
// ripple peeling decoding, and Gaussian-elimination fallback on the unpeeled
// core for minimal reception overhead.
template <typename GaloisField>
class LTCode {
 public:
  using field_type = GaloisField;
  using storage_type = typename GaloisField::element_type;
  using symbol_type = FountainSymbol<GaloisField>;
  using result_type = FountainDecodingResult<GaloisField>;
  using matrix_type = GFMatrix<GaloisField>;
  using array_type = GFArray<GaloisField>;

  enum class DegreeDistribution : uint8_t { kRobustSoliton, kIdealSoliton };

  LTCode(std::shared_ptr<GaloisField> field, std::size_t num_source_symbols,
         std::size_t symbol_size = 1,
         DegreeDistribution distribution = DegreeDistribution::kRobustSoliton,
         double c = 0.1, double delta = 0.05,
         uint64_t seed = 0x517CC1B727220A95ULL)
      : field_(std::move(field)),
        k_(num_source_symbols),
        symbol_size_(symbol_size),
        distribution_(distribution),
        c_(c),
        delta_(delta),
        seed_(seed) {
    if (!field_) {
      throw std::invalid_argument("LTCode: field pointer cannot be null");
    }
    if (k_ == 0) {
      throw std::invalid_argument("LTCode: num_source_symbols must be >= 1");
    }
    if (symbol_size_ == 0) {
      throw std::invalid_argument("LTCode: symbol_size must be >= 1");
    }
    if (!(c_ > 0.0)) {
      throw std::invalid_argument("LTCode: c must be positive");
    }
    if (!(delta_ > 0.0 && delta_ < 1.0)) {
      throw std::invalid_argument("LTCode: delta must lie in (0, 1)");
    }
    BuildDegreeCdf();
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }
  std::size_t NumSourceSymbols() const noexcept { return k_; }
  std::size_t SymbolSize() const noexcept { return symbol_size_; }
  DegreeDistribution Distribution() const noexcept { return distribution_; }
  const std::vector<double> &DegreeCdf() const noexcept { return cdf_; }

  // Deterministically generates the encoding symbol for a given ESI from the
  // K x symbol_size source block.
  symbol_type GenerateSymbol(const matrix_type &source_block,
                             uint32_t esi) const {
    ValidateSourceBlock(source_block);
    detail::FountainPrng rng(detail::MixSeed(seed_, 0, esi));
    const std::size_t deg = SampleDegree(rng);

    symbol_type sym;
    sym.esi = esi;
    sym.neighbors = detail::SampleDistinctIndices(k_, deg, 0, rng);
    sym.coeffs.reserve(deg);
    for (std::size_t j = 0; j < deg; ++j) {
      sym.coeffs.push_back(detail::SampleNonZeroElement(*field_, rng));
    }

    using RawArray = typename array_type::raw_array_type;
    using ShapeType = typename array_type::shape_type;
    RawArray payload_raw = RawArray::from_shape(ShapeType{symbol_size_});
    payload_raw.fill(field_->AdditiveIdentity());
    std::span<storage_type> dst(payload_raw.data(), symbol_size_);
    const storage_type *src = source_block.RawData().data();
    for (std::size_t j = 0; j < deg; ++j) {
      std::span<const storage_type> row(src + sym.neighbors[j] * symbol_size_,
                                        symbol_size_);
      field_->FmaArray(dst, row, sym.coeffs[j]);
    }
    sym.payload = array_type::FromRaw(field_, std::move(payload_raw));
    return sym;
  }

  std::vector<symbol_type> GenerateSymbols(const matrix_type &source_block,
                                           std::size_t count,
                                           uint32_t start_esi = 0) const {
    ValidateSourceBlock(source_block);
    std::vector<symbol_type> out;
    out.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
      out.push_back(
          GenerateSymbol(source_block, start_esi + static_cast<uint32_t>(i)));
    }
    return out;
  }

  // O(K ln(K/delta) * symbol_size) degree-1 ripple peeling decoder.
  result_type DecodePeeling(
      const std::vector<symbol_type> &received_symbols) const {
    return DecodeInternal(received_symbols, /*allow_gaussian_fallback=*/false);
  }

  // Runs ripple peeling first, then solves any remaining unpeeled core via
  // dense Gaussian elimination over GaloisField.
  result_type DecodeWithGaussianFallback(
      const std::vector<symbol_type> &received_symbols) const {
    return DecodeInternal(received_symbols, /*allow_gaussian_fallback=*/true);
  }

 private:
  void ValidateSourceBlock(const matrix_type &source_block) const {
    if (!source_block.HasField() || source_block.Field() != field_) {
      throw std::invalid_argument("LTCode: source_block field mismatch");
    }
    if (source_block.Dimension() != 2 || source_block.Shape(0) != k_ ||
        source_block.Shape(1) != symbol_size_) {
      throw std::invalid_argument(
          "LTCode: source_block must have shape (K, symbol_size)");
    }
  }

  void BuildDegreeCdf() {
    std::vector<double> rho(k_ + 1, 0.0);
    rho[1] = 1.0 / static_cast<double>(k_);
    for (std::size_t d = 2; d <= k_; ++d) {
      rho[d] = 1.0 / (static_cast<double>(d) * static_cast<double>(d - 1));
    }

    std::vector<double> pmf(k_ + 1, 0.0);
    if (distribution_ == DegreeDistribution::kIdealSoliton || k_ == 1) {
      pmf = rho;
    } else {
      const double kd = static_cast<double>(k_);
      const double R = c_ * std::log(kd / delta_) * std::sqrt(kd);
      std::size_t M =
          (R > 0.0) ? static_cast<std::size_t>(std::floor(kd / R)) : k_;
      if (M < 1) {
        M = 1;
      }
      if (M > k_) {
        M = k_;
      }

      std::vector<double> tau(k_ + 1, 0.0);
      for (std::size_t d = 1; d < M; ++d) {
        tau[d] = R / (static_cast<double>(d) * kd);
      }
      if (R > 0.0 && R / delta_ > 1.0) {
        tau[M] = R * std::log(R / delta_) / kd;
      }

      double Z = 0.0;
      for (std::size_t d = 1; d <= k_; ++d) {
        pmf[d] = rho[d] + tau[d];
        Z += pmf[d];
      }
      for (std::size_t d = 1; d <= k_; ++d) {
        pmf[d] /= Z;
      }
    }

    cdf_.resize(k_);
    double acc = 0.0;
    for (std::size_t d = 1; d <= k_; ++d) {
      acc += pmf[d];
      cdf_[d - 1] = acc;
    }
    cdf_.back() = 1.0;
  }

  std::size_t SampleDegree(detail::FountainPrng &rng) const noexcept {
    const double u = rng.NextDouble();
    const auto it = std::lower_bound(cdf_.begin(), cdf_.end(), u);
    const std::size_t idx =
        static_cast<std::size_t>(std::distance(cdf_.begin(), it));
    return std::min(idx + 1, k_);
  }

  result_type DecodeInternal(const std::vector<symbol_type> &received_symbols,
                             bool allow_gaussian_fallback) const {
    result_type res;
    res.symbols_used = received_symbols.size();
    if (received_symbols.size() < k_) {
      res.error_message = "LTCode: fewer received symbols than K";
      return res;
    }

    const storage_type zero = field_->AdditiveIdentity();
    const std::size_t num_rows = received_symbols.size();

    struct ActiveEq {
      std::vector<std::pair<std::size_t, storage_type>> terms;
      std::vector<storage_type> rhs;
    };

    std::vector<ActiveEq> eqs(num_rows);
    std::vector<std::vector<std::size_t>> col_to_rows(k_);
    std::queue<std::size_t> ripple;

    for (std::size_t r = 0; r < num_rows; ++r) {
      const auto &sym = received_symbols[r];
      if (!sym.payload.HasField() || sym.payload.Field() != field_ ||
          sym.payload.Dimension() != 1 || sym.payload.Size() != symbol_size_) {
        throw std::invalid_argument(
            "LTCode::Decode: received symbol payload shape or field mismatch");
      }
      if (sym.neighbors.size() != sym.coeffs.size()) {
        throw std::invalid_argument(
            "LTCode::Decode: symbol neighbors and coeffs size mismatch");
      }
      eqs[r].terms.reserve(sym.neighbors.size());
      for (std::size_t j = 0; j < sym.neighbors.size(); ++j) {
        const std::size_t col = sym.neighbors[j];
        if (col >= k_) {
          throw std::invalid_argument(
              "LTCode::Decode: neighbor index out of range");
        }
        if (!field_->IsZero(sym.coeffs[j])) {
          eqs[r].terms.emplace_back(col, sym.coeffs[j]);
          col_to_rows[col].push_back(r);
        }
      }
      const auto p_span = sym.payload.RawSpan();
      eqs[r].rhs.assign(p_span.begin(), p_span.end());
      if (eqs[r].terms.size() == 1) {
        ripple.push(r);
      }
    }

    std::vector<bool> decoded(k_, false);
    std::vector<storage_type> recovered(k_ * symbol_size_, zero);
    std::size_t num_decoded = 0;

    while (!ripple.empty()) {
      const std::size_t r = ripple.front();
      ripple.pop();
      if (eqs[r].terms.size() != 1) {
        continue;
      }
      const std::size_t col = eqs[r].terms[0].first;
      const storage_type coeff = eqs[r].terms[0].second;
      eqs[r].terms.clear();
      if (decoded[col]) {
        continue;
      }

      // Recover source symbol `col`.
      std::span<storage_type> dst(recovered.data() + col * symbol_size_,
                                  symbol_size_);
      if (coeff == field_->MultiplicativeIdentity()) {
        std::copy(eqs[r].rhs.begin(), eqs[r].rhs.end(), dst.begin());
      } else {
        const storage_type inv_c = field_->Inv(coeff);
        field_->ScalarMulArray(dst, std::span<const storage_type>(eqs[r].rhs),
                               inv_c);
      }
      decoded[col] = true;
      ++num_decoded;

      // Substitute into all equations containing `col`.
      std::span<const storage_type> x_col(dst.data(), symbol_size_);
      for (std::size_t qr : col_to_rows[col]) {
        auto &terms = eqs[qr].terms;
        for (std::size_t i = 0; i < terms.size(); ++i) {
          if (terms[i].first == col) {
            const storage_type neg_c = field_->Neg(terms[i].second);
            field_->FmaArray(std::span<storage_type>(eqs[qr].rhs), x_col,
                             neg_c);
            terms[i] = terms.back();
            terms.pop_back();
            if (terms.size() == 1) {
              ripple.push(qr);
            }
            break;
          }
        }
      }
    }

    if (num_decoded < k_) {
      if (!allow_gaussian_fallback) {
        res.error_message = "LTCode::DecodePeeling: ripple exhausted with " +
                            std::to_string(k_ - num_decoded) +
                            " unpeeled symbols";
        return res;
      }

      // Build dense reduced system over the unpeeled columns.
      std::vector<std::size_t> unpeeled_cols;
      std::vector<std::size_t> col_to_local(k_, k_);
      unpeeled_cols.reserve(k_ - num_decoded);
      for (std::size_t c = 0; c < k_; ++c) {
        if (!decoded[c]) {
          col_to_local[c] = unpeeled_cols.size();
          unpeeled_cols.push_back(c);
        }
      }
      const std::size_t u = unpeeled_cols.size();

      std::vector<std::size_t> active_rows;
      for (std::size_t r = 0; r < num_rows; ++r) {
        if (!eqs[r].terms.empty()) {
          active_rows.push_back(r);
        }
      }
      const std::size_t m_rem = active_rows.size();
      if (m_rem < u) {
        res.error_message =
            "LTCode::DecodeWithGaussianFallback: insufficient equations for "
            "unpeeled core";
        return res;
      }

      std::vector<storage_type> a_mat(m_rem * u, zero);
      std::vector<storage_type> b_mat(m_rem * symbol_size_, zero);
      for (std::size_t i = 0; i < m_rem; ++i) {
        const auto &eq = eqs[active_rows[i]];
        for (const auto &[col, coeff] : eq.terms) {
          a_mat[i * u + col_to_local[col]] = coeff;
        }
        std::copy(
            eq.rhs.begin(), eq.rhs.end(),
            b_mat.begin() + static_cast<std::ptrdiff_t>(i * symbol_size_));
      }

      std::vector<storage_type> x_sol;
      if (!detail::DenseGaussianSolveOverdetermined(
              *field_, m_rem, u, symbol_size_, a_mat, b_mat, x_sol)) {
        res.error_message =
            "LTCode::DecodeWithGaussianFallback: unpeeled core is rank "
            "deficient";
        return res;
      }

      for (std::size_t j = 0; j < u; ++j) {
        const std::size_t orig_c = unpeeled_cols[j];
        std::copy(
            x_sol.begin() + static_cast<std::ptrdiff_t>(j * symbol_size_),
            x_sol.begin() + static_cast<std::ptrdiff_t>((j + 1) * symbol_size_),
            recovered.begin() +
                static_cast<std::ptrdiff_t>(orig_c * symbol_size_));
      }
      res.used_gaussian_fallback = true;
      res.num_inactivations = u;
    }

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray out_raw = RawArray::from_shape(ShapeType{k_, symbol_size_});
    std::copy(recovered.begin(), recovered.end(), out_raw.data());
    res.source_block = matrix_type::FromRaw(field_, std::move(out_raw));
    res.success = true;
    return res;
  }

  std::shared_ptr<GaloisField> field_;
  std::size_t k_ = 0;
  std::size_t symbol_size_ = 1;
  DegreeDistribution distribution_ = DegreeDistribution::kRobustSoliton;
  double c_ = 0.1;
  double delta_ = 0.05;
  uint64_t seed_ = 0;
  std::vector<double> cdf_;
};

// ============================================================================
// RaptorQCode: RFC 6330-Style Systematic Fountain Code with Inactivation
// ============================================================================
//
// Precode dimensions L = K + S + H:
//   - S sparse circulant degree-3 LDPC check rows
//   - H dense HDPC check rows (full-rank mixing over GF(2^8) or GF(2))
//   - Deterministic LT + permanent-inactivation (PI) tuple generator
//
// Systematic encoding solves A_enc * C_int = [0_{S+H}; D_src] via 3-phase
// InactivationSolve, emits systematic source symbols (0 <= esi < K) verbatim,
// and generates rateless repair symbols (esi >= K) on demand.
template <typename GaloisField>
class RaptorQCode {
 public:
  using field_type = GaloisField;
  using storage_type = typename GaloisField::element_type;
  using symbol_type = FountainSymbol<GaloisField>;
  using result_type = FountainDecodingResult<GaloisField>;
  using matrix_type = GFMatrix<GaloisField>;
  using array_type = GFArray<GaloisField>;

  struct SparseRow {
    // Sorted by column index; all coefficients nonzero.
    std::vector<std::pair<std::size_t, storage_type>> entries;
  };

  struct InactivationResult {
    std::vector<storage_type> solution;  // shape: L * rhs_cols
    std::size_t num_inactivations = 0;
    bool success = false;
    std::string error_message;
  };

  RaptorQCode(std::shared_ptr<GaloisField> field,
              std::size_t num_source_symbols, std::size_t symbol_size = 1,
              uint64_t seed = 0x243F6A8885A308D3ULL)
      : field_(std::move(field)),
        k_(num_source_symbols),
        symbol_size_(symbol_size),
        seed_(seed) {
    if (!field_) {
      throw std::invalid_argument("RaptorQCode: field pointer cannot be null");
    }
    if (k_ == 0) {
      throw std::invalid_argument(
          "RaptorQCode: num_source_symbols must be >= 1");
    }
    if (symbol_size_ == 0) {
      throw std::invalid_argument("RaptorQCode: symbol_size must be >= 1");
    }

    // Choose RFC 6330-style precode dimensions S (prime LDPC count) and H
    // (dense HDPC count).
    const std::size_t min_s = std::max<std::size_t>(
        7, static_cast<std::size_t>(
               std::ceil(std::sqrt(2.0 * static_cast<double>(k_)))) +
               3);
    s_ = detail::NextPrimeAtLeast(min_s);
    const std::size_t log_ks = static_cast<std::size_t>(
        std::ceil(std::log2(static_cast<double>(k_ + s_ + 1))));
    h_ = std::max<std::size_t>(8, log_ks + 3);
    w_ = k_ + s_;
    l_ = k_ + s_ + h_;

    BuildPrecodeRows();
    SelectInvertibleSystematicSalt();
  }

  std::shared_ptr<GaloisField> Field() const { return field_; }
  std::size_t NumSourceSymbols() const noexcept { return k_; }
  std::size_t SymbolSize() const noexcept { return symbol_size_; }
  std::size_t NumLdpcSymbols() const noexcept { return s_; }
  std::size_t NumHdpcSymbols() const noexcept { return h_; }
  std::size_t NumIntermediateSymbols() const noexcept { return l_; }

  // 3-Phase Inactivation Solver for M x L sparse system A * X = B:
  //   Phase 1: Minimum-active-degree triangularization with r = 2 max-component
  //            tie-breaking; moves non-pivot active columns to U.
  //   Phase 2: Dense Gaussian elimination on the small u x u inactivated core
  //            using field_->FmaArray.
  //   Phase 3: Reverse triangular back-substitution for active columns V.
  InactivationResult InactivationSolve(
      const std::vector<SparseRow> &rows,
      const std::vector<storage_type> &rhs_flat, std::size_t rhs_cols) const {
    InactivationResult out;
    const std::size_t m = rows.size();
    if (m < l_) {
      out.error_message = "InactivationSolve: fewer equations than L";
      return out;
    }

    const storage_type zero = field_->AdditiveIdentity();
    const storage_type one = field_->MultiplicativeIdentity();

    enum class ColState : uint8_t { kActive, kInactivated, kPivoted };
    std::vector<ColState> col_state(l_, ColState::kActive);
    std::vector<bool> row_pivoted(m, false);
    std::vector<std::size_t> active_deg(m, 0);
    std::vector<std::vector<std::size_t>> col_rows(l_);

    // Copy sparse rows and dense RHS for in-place Phase 1 elimination.
    std::vector<std::vector<std::pair<std::size_t, storage_type>>> work_rows(m);
    for (std::size_t r = 0; r < m; ++r) {
      work_rows[r] = rows[r].entries;
      active_deg[r] = work_rows[r].size();
      for (const auto &[col, _] : work_rows[r]) {
        col_rows[col].push_back(r);
      }
    }
    std::vector<storage_type> work_rhs = rhs_flat;

    std::vector<std::pair<std::size_t, std::size_t>> pivot_order;
    pivot_order.reserve(l_);

    // Helper: Union-Find for RFC 6330 r == 2 max-component tie-breaking.
    std::vector<std::size_t> uf_parent(l_);
    std::vector<std::size_t> uf_size(l_);
    auto uf_find = [&](std::size_t x) {
      while (uf_parent[x] != x) {
        uf_parent[x] = uf_parent[uf_parent[x]];
        x = uf_parent[x];
      }
      return x;
    };

    std::size_t num_active_cols = l_;
    std::size_t num_inactivated = 0;

    // ------------------------------------------------------------------------
    // Phase 1: Triangularization
    // ------------------------------------------------------------------------
    while (num_active_cols > 0) {
      std::size_t min_deg = std::numeric_limits<std::size_t>::max();
      std::size_t chosen_row = m;

      for (std::size_t r = 0; r < m; ++r) {
        if (row_pivoted[r] || active_deg[r] == 0) {
          continue;
        }
        if (active_deg[r] < min_deg) {
          min_deg = active_deg[r];
          chosen_row = r;
          if (min_deg == 1) {
            break;
          }
        }
      }

      if (chosen_row == m) {
        // No unpivoted row has any remaining active column: move all remaining
        // active columns into U.
        for (std::size_t c = 0; c < l_; ++c) {
          if (col_state[c] == ColState::kActive) {
            col_state[c] = ColState::kInactivated;
            ++num_inactivated;
          }
        }
        break;
      }

      // Max-component tie-breaking when min_deg == 2.
      if (min_deg == 2) {
        std::iota(uf_parent.begin(), uf_parent.end(), std::size_t{0});
        std::fill(uf_size.begin(), uf_size.end(), std::size_t{1});
        std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>>
            deg2_edges;
        for (std::size_t r = 0; r < m; ++r) {
          if (row_pivoted[r] || active_deg[r] != 2) {
            continue;
          }
          std::size_t c0 = l_;
          std::size_t c1 = l_;
          for (const auto &[col, _] : work_rows[r]) {
            if (col_state[col] == ColState::kActive) {
              if (c0 == l_) {
                c0 = col;
              } else {
                c1 = col;
                break;
              }
            }
          }
          if (c1 < l_) {
            deg2_edges.push_back({r, {c0, c1}});
            std::size_t ra = uf_find(c0);
            std::size_t rb = uf_find(c1);
            if (ra != rb) {
              if (uf_size[ra] < uf_size[rb]) {
                std::swap(ra, rb);
              }
              uf_parent[rb] = ra;
              uf_size[ra] += uf_size[rb];
            }
          }
        }
        std::size_t best_comp = 0;
        for (const auto &[r, endpoints] : deg2_edges) {
          const std::size_t comp_sz = uf_size[uf_find(endpoints.first)];
          if (comp_sz > best_comp) {
            best_comp = comp_sz;
            chosen_row = r;
          }
        }
      }

      // Partition active columns of chosen_row: pivot the first, inactivate the
      // remaining (min_deg - 1).
      std::size_t p_col = l_;
      storage_type p_coeff = zero;
      std::vector<std::size_t> to_inactivate;
      to_inactivate.reserve(min_deg > 0 ? min_deg - 1 : 0);
      for (const auto &[col, coeff] : work_rows[chosen_row]) {
        if (col_state[col] == ColState::kActive) {
          if (p_col == l_) {
            p_col = col;
            p_coeff = coeff;
          } else {
            to_inactivate.push_back(col);
          }
        }
      }

      for (std::size_t ic : to_inactivate) {
        col_state[ic] = ColState::kInactivated;
        ++num_inactivated;
        --num_active_cols;
        for (std::size_t r : col_rows[ic]) {
          if (!row_pivoted[r] && active_deg[r] > 0) {
            --active_deg[r];
          }
        }
      }

      col_state[p_col] = ColState::kPivoted;
      --num_active_cols;
      row_pivoted[chosen_row] = true;
      pivot_order.emplace_back(chosen_row, p_col);

      // Normalize chosen_row so p_col has coefficient 1.
      if (p_coeff != one) {
        const storage_type inv_p = field_->Inv(p_coeff);
        for (auto &[col, coeff] : work_rows[chosen_row]) {
          coeff = field_->Mul(coeff, inv_p);
        }
        if (rhs_cols > 0) {
          std::span<storage_type> rhs_span(
              work_rhs.data() + chosen_row * rhs_cols, rhs_cols);
          field_->ScalarMulArray(rhs_span, rhs_span, inv_p);
        }
      }

      // Extract the inactivated tail of chosen_row (all entries != p_col).
      std::vector<std::pair<std::size_t, storage_type>> piv_tail;
      piv_tail.reserve(work_rows[chosen_row].size());
      for (const auto &entry : work_rows[chosen_row]) {
        if (entry.first != p_col) {
          piv_tail.push_back(entry);
        }
      }

      std::span<const storage_type> piv_rhs(
          work_rhs.data() + chosen_row * rhs_cols, rhs_cols);

      // Eliminate p_col from all other unpivoted rows containing p_col.
      for (std::size_t qr : col_rows[p_col]) {
        if (row_pivoted[qr]) {
          continue;
        }
        auto &q_entries = work_rows[qr];
        auto it = std::lower_bound(q_entries.begin(), q_entries.end(), p_col,
                                   [](const auto &entry, std::size_t target) {
                                     return entry.first < target;
                                   });
        if (it == q_entries.end() || it->first != p_col) {
          continue;
        }
        const storage_type factor = it->second;
        q_entries.erase(it);
        if (active_deg[qr] > 0) {
          --active_deg[qr];
        }

        const storage_type neg_factor = field_->Neg(factor);
        if (!piv_tail.empty()) {
          // Merge q_entries += neg_factor * piv_tail (both sorted by col).
          std::vector<std::pair<std::size_t, storage_type>> merged;
          merged.reserve(q_entries.size() + piv_tail.size());
          std::size_t iq = 0;
          std::size_t ip = 0;
          while (iq < q_entries.size() && ip < piv_tail.size()) {
            if (q_entries[iq].first < piv_tail[ip].first) {
              merged.push_back(q_entries[iq++]);
            } else if (piv_tail[ip].first < q_entries[iq].first) {
              const storage_type val =
                  field_->Mul(neg_factor, piv_tail[ip].second);
              if (!field_->IsZero(val)) {
                merged.emplace_back(piv_tail[ip].first, val);
              }
              ++ip;
            } else {
              const storage_type val =
                  field_->Add(q_entries[iq].second,
                              field_->Mul(neg_factor, piv_tail[ip].second));
              if (!field_->IsZero(val)) {
                merged.emplace_back(q_entries[iq].first, val);
              }
              ++iq;
              ++ip;
            }
          }
          while (iq < q_entries.size()) {
            merged.push_back(q_entries[iq++]);
          }
          while (ip < piv_tail.size()) {
            const storage_type val =
                field_->Mul(neg_factor, piv_tail[ip].second);
            if (!field_->IsZero(val)) {
              merged.emplace_back(piv_tail[ip].first, val);
            }
            ++ip;
          }
          q_entries = std::move(merged);
        }

        if (rhs_cols > 0) {
          field_->FmaArray(std::span<storage_type>(
                               work_rhs.data() + qr * rhs_cols, rhs_cols),
                           piv_rhs, neg_factor);
        }
      }
    }

    // ------------------------------------------------------------------------
    // Phase 2: Dense Core Gaussian Elimination over Inactivated Columns U
    // ------------------------------------------------------------------------
    std::vector<std::size_t> u_cols;
    std::vector<std::size_t> col_to_u(l_, l_);
    u_cols.reserve(num_inactivated);
    for (std::size_t c = 0; c < l_; ++c) {
      if (col_state[c] == ColState::kInactivated) {
        col_to_u[c] = u_cols.size();
        u_cols.push_back(c);
      }
    }
    const std::size_t u = u_cols.size();
    out.num_inactivations = u;

    std::vector<storage_type> sol_flat(l_ * rhs_cols, zero);

    if (u > 0) {
      std::vector<std::size_t> core_rows;
      core_rows.reserve(m - pivot_order.size());
      for (std::size_t r = 0; r < m; ++r) {
        if (!row_pivoted[r]) {
          core_rows.push_back(r);
        }
      }
      const std::size_t m_u = core_rows.size();
      if (m_u < u) {
        out.error_message =
            "InactivationSolve: insufficient dense-core equations";
        return out;
      }

      std::vector<storage_type> d_mat(m_u * u, zero);
      std::vector<storage_type> r_mat(m_u * rhs_cols, zero);
      for (std::size_t i = 0; i < m_u; ++i) {
        const std::size_t orig_r = core_rows[i];
        for (const auto &[col, coeff] : work_rows[orig_r]) {
          const std::size_t u_idx = col_to_u[col];
          if (u_idx < u) {
            d_mat[i * u + u_idx] = coeff;
          }
        }
        if (rhs_cols > 0) {
          std::copy(
              work_rhs.begin() + static_cast<std::ptrdiff_t>(orig_r * rhs_cols),
              work_rhs.begin() +
                  static_cast<std::ptrdiff_t>((orig_r + 1) * rhs_cols),
              r_mat.begin() + static_cast<std::ptrdiff_t>(i * rhs_cols));
        }
      }

      std::vector<storage_type> u_sol;
      if (!detail::DenseGaussianSolveOverdetermined(*field_, m_u, u, rhs_cols,
                                                    d_mat, r_mat, u_sol)) {
        out.error_message = "InactivationSolve: dense core is rank deficient";
        return out;
      }

      if (rhs_cols > 0) {
        for (std::size_t j = 0; j < u; ++j) {
          const std::size_t col = u_cols[j];
          std::copy(
              u_sol.begin() + static_cast<std::ptrdiff_t>(j * rhs_cols),
              u_sol.begin() + static_cast<std::ptrdiff_t>((j + 1) * rhs_cols),
              sol_flat.begin() + static_cast<std::ptrdiff_t>(col * rhs_cols));
        }
      }
    }

    // ------------------------------------------------------------------------
    // Phase 3: Reverse Triangular Back-Substitution for Active Columns V
    // ------------------------------------------------------------------------
    if (rhs_cols > 0) {
      for (std::size_t idx = pivot_order.size(); idx-- > 0;) {
        const auto [p_row, p_col] = pivot_order[idx];
        std::span<storage_type> dst(sol_flat.data() + p_col * rhs_cols,
                                    rhs_cols);
        std::copy(
            work_rhs.begin() + static_cast<std::ptrdiff_t>(p_row * rhs_cols),
            work_rhs.begin() +
                static_cast<std::ptrdiff_t>((p_row + 1) * rhs_cols),
            dst.begin());
        for (const auto &[col, coeff] : work_rows[p_row]) {
          if (col == p_col) {
            continue;
          }
          const storage_type neg_c = field_->Neg(coeff);
          std::span<const storage_type> src(sol_flat.data() + col * rhs_cols,
                                            rhs_cols);
          field_->FmaArray(dst, src, neg_c);
        }
      }
    }

    out.solution = std::move(sol_flat);
    out.success = true;
    return out;
  }

  // Computes the L x symbol_size intermediate symbol block C_int satisfying
  // A_enc * C_int = [0_{S+H}; D_src] via InactivationSolve.
  matrix_type EncodeIntermediate(const matrix_type &source_block) const {
    if (!source_block.HasField() || source_block.Field() != field_) {
      throw std::invalid_argument(
          "RaptorQCode::EncodeIntermediate: source_block field mismatch");
    }
    if (source_block.Dimension() != 2 || source_block.Shape(0) != k_ ||
        source_block.Shape(1) != symbol_size_) {
      throw std::invalid_argument(
          "RaptorQCode::EncodeIntermediate: source_block must have shape (K, "
          "symbol_size)");
    }

    std::vector<SparseRow> rows;
    rows.reserve(l_);
    for (const auto &r : ldpc_rows_) {
      rows.push_back(r);
    }
    for (const auto &r : hdpc_rows_) {
      rows.push_back(r);
    }
    for (uint32_t esi = 0; esi < static_cast<uint32_t>(k_); ++esi) {
      rows.push_back(BuildLtRow(esi));
    }

    const storage_type zero = field_->AdditiveIdentity();
    std::vector<storage_type> rhs(l_ * symbol_size_, zero);
    const storage_type *src = source_block.RawData().data();
    std::copy(
        src, src + k_ * symbol_size_,
        rhs.begin() + static_cast<std::ptrdiff_t>((s_ + h_) * symbol_size_));

    auto res = InactivationSolve(rows, rhs, symbol_size_);
    if (!res.success) {
      throw std::runtime_error(
          "RaptorQCode::EncodeIntermediate: failed to solve systematic "
          "precode: " +
          res.error_message);
    }

    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray out_raw = RawArray::from_shape(ShapeType{l_, symbol_size_});
    std::copy(res.solution.begin(), res.solution.end(), out_raw.data());
    return matrix_type::FromRaw(field_, std::move(out_raw));
  }

  // Emits systematic source symbol (0 <= esi < K) or rateless repair symbol
  // (esi >= K) from the intermediate symbol block C_int.
  symbol_type GenerateSymbol(const matrix_type &intermediate,
                             uint32_t esi) const {
    if (!intermediate.HasField() || intermediate.Field() != field_) {
      throw std::invalid_argument(
          "RaptorQCode::GenerateSymbol: intermediate field mismatch");
    }
    if (intermediate.Dimension() != 2 || intermediate.Shape(0) != l_ ||
        intermediate.Shape(1) != symbol_size_) {
      throw std::invalid_argument(
          "RaptorQCode::GenerateSymbol: intermediate must have shape (L, "
          "symbol_size)");
    }

    const SparseRow row = BuildLtRow(esi);
    symbol_type sym;
    sym.esi = esi;
    sym.neighbors.reserve(row.entries.size());
    sym.coeffs.reserve(row.entries.size());

    using RawArray = typename array_type::raw_array_type;
    using ShapeType = typename array_type::shape_type;
    RawArray payload_raw = RawArray::from_shape(ShapeType{symbol_size_});
    payload_raw.fill(field_->AdditiveIdentity());
    std::span<storage_type> dst(payload_raw.data(), symbol_size_);
    const storage_type *c_data = intermediate.RawData().data();

    for (const auto &[col, coeff] : row.entries) {
      sym.neighbors.push_back(col);
      sym.coeffs.push_back(coeff);
      std::span<const storage_type> c_row(c_data + col * symbol_size_,
                                          symbol_size_);
      field_->FmaArray(dst, c_row, coeff);
    }
    sym.payload = array_type::FromRaw(field_, std::move(payload_raw));
    return sym;
  }

  std::vector<symbol_type> GenerateSymbols(const matrix_type &intermediate,
                                           std::size_t count,
                                           uint32_t start_esi = 0) const {
    std::vector<symbol_type> out;
    out.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
      out.push_back(
          GenerateSymbol(intermediate, start_esi + static_cast<uint32_t>(i)));
    }
    return out;
  }

  // Reconstructs C_int via InactivationSolve from any subset of systematic
  // and/or repair symbols, then recovers the exact K x symbol_size source
  // block.
  result_type Decode(const std::vector<symbol_type> &received_symbols) const {
    result_type res;
    res.symbols_used = received_symbols.size();
    if (received_symbols.size() < k_) {
      res.error_message = "RaptorQCode::Decode: fewer received symbols than K";
      return res;
    }

    const std::size_t total_rows = s_ + h_ + received_symbols.size();
    std::vector<SparseRow> rows;
    rows.reserve(total_rows);
    for (const auto &r : ldpc_rows_) {
      rows.push_back(r);
    }
    for (const auto &r : hdpc_rows_) {
      rows.push_back(r);
    }

    const storage_type zero = field_->AdditiveIdentity();
    std::vector<storage_type> rhs(total_rows * symbol_size_, zero);

    for (std::size_t i = 0; i < received_symbols.size(); ++i) {
      const auto &sym = received_symbols[i];
      if (!sym.payload.HasField() || sym.payload.Field() != field_ ||
          sym.payload.Dimension() != 1 || sym.payload.Size() != symbol_size_) {
        throw std::invalid_argument(
            "RaptorQCode::Decode: received symbol payload shape or field "
            "mismatch");
      }
      rows.push_back(BuildLtRow(sym.esi));
      const auto p_span = sym.payload.RawSpan();
      std::copy(p_span.begin(), p_span.end(),
                rhs.begin() +
                    static_cast<std::ptrdiff_t>((s_ + h_ + i) * symbol_size_));
    }

    auto solve_res = InactivationSolve(rows, rhs, symbol_size_);
    res.num_inactivations = solve_res.num_inactivations;
    if (!solve_res.success) {
      res.error_message = "RaptorQCode::Decode: " + solve_res.error_message;
      return res;
    }

    // Recover systematic source symbols D_i = g_i * C_int for 0 <= i < K.
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    RawArray out_raw = RawArray::from_shape(ShapeType{k_, symbol_size_});
    out_raw.fill(zero);
    storage_type *out_data = out_raw.data();
    const storage_type *c_data = solve_res.solution.data();

    for (uint32_t esi = 0; esi < static_cast<uint32_t>(k_); ++esi) {
      const SparseRow lt_row = BuildLtRow(esi);
      std::span<storage_type> dst(out_data + esi * symbol_size_, symbol_size_);
      for (const auto &[col, coeff] : lt_row.entries) {
        std::span<const storage_type> c_row(c_data + col * symbol_size_,
                                            symbol_size_);
        field_->FmaArray(dst, c_row, coeff);
      }
    }

    res.source_block = matrix_type::FromRaw(field_, std::move(out_raw));
    res.success = true;
    return res;
  }

 private:
  void BuildPrecodeRows() {
    const storage_type one = field_->MultiplicativeIdentity();

    // 1. S sparse circulant degree-3 LDPC check rows.
    ldpc_rows_.assign(s_, SparseRow{});
    for (std::size_t col = 0; col < k_; ++col) {
      const std::size_t a = 1 + (col / s_) % (s_ - 1);
      const std::size_t b0 = col % s_;
      const std::size_t b1 = (b0 + a) % s_;
      const std::size_t b2 = (b0 + 2 * a) % s_;
      ldpc_rows_[b0].entries.emplace_back(col, one);
      if (b1 != b0) {
        ldpc_rows_[b1].entries.emplace_back(col, one);
      }
      if (b2 != b0 && b2 != b1) {
        ldpc_rows_[b2].entries.emplace_back(col, one);
      }
    }
    // Identity I_S on LDPC parity columns [K, K + S).
    for (std::size_t r = 0; r < s_; ++r) {
      ldpc_rows_[r].entries.emplace_back(k_ + r, one);
    }
    // Circulant connections to HDPC columns [K + S, L).
    for (std::size_t j = 0; j < h_; ++j) {
      const std::size_t r0 = j % s_;
      const std::size_t r1 = (j + 1) % s_;
      ldpc_rows_[r0].entries.emplace_back(w_ + j, one);
      if (r1 != r0) {
        ldpc_rows_[r1].entries.emplace_back(w_ + j, one);
      }
    }
    for (auto &row : ldpc_rows_) {
      std::sort(row.entries.begin(), row.entries.end(),
                [](const auto &lhs, const auto &rhs) {
                  return lhs.first < rhs.first;
                });
    }

    // 2. H dense HDPC check rows (full-rank mixing over GaloisField) + I_H.
    hdpc_rows_.assign(h_, SparseRow{});
    detail::FountainPrng hdpc_rng(detail::MixSeed(seed_, 0x48445043ULL, l_));
    const bool is_binary = (field_->Order() == 2);
    for (std::size_t r = 0; r < h_; ++r) {
      hdpc_rows_[r].entries.reserve(w_ + 1);
      for (std::size_t col = 0; col < w_; ++col) {
        if (is_binary) {
          if ((hdpc_rng.NextU64() & 1ULL) != 0) {
            hdpc_rows_[r].entries.emplace_back(col, one);
          }
        } else {
          const storage_type c =
              detail::SampleNonZeroElement(*field_, hdpc_rng);
          hdpc_rows_[r].entries.emplace_back(col, c);
        }
      }
      hdpc_rows_[r].entries.emplace_back(w_ + r, one);
    }
  }

  // RFC 6330-inspired LT degree distribution over the W = K + S LT columns,
  // plus d_pi = min(H, 2 or 3) connections to the H HDPC/PI columns.
  SparseRow BuildLtRow(uint32_t esi) const {
    detail::FountainPrng rng(detail::MixSeed(seed_, salt_, esi));
    const double u = rng.NextDouble();
    std::size_t d_lt = 2;
    if (u < 0.01) {
      d_lt = 1;
    } else if (u < 0.49) {
      d_lt = 2;
    } else if (u < 0.66) {
      d_lt = 3;
    } else if (u < 0.74) {
      d_lt = 4;
    } else if (u < 0.89) {
      d_lt = 10;
    } else if (u < 0.94) {
      d_lt = 11;
    } else {
      d_lt = std::min<std::size_t>(w_, 30);
    }
    d_lt = std::min(d_lt, w_);

    const std::size_t d_pi = std::min<std::size_t>(h_, (d_lt <= 3) ? 2 : 3);

    const auto lt_cols = detail::SampleDistinctIndices(w_, d_lt, 0, rng);
    const auto pi_cols = detail::SampleDistinctIndices(h_, d_pi, w_, rng);

    SparseRow row;
    row.entries.reserve(lt_cols.size() + pi_cols.size());
    for (std::size_t col : lt_cols) {
      row.entries.emplace_back(col, detail::SampleNonZeroElement(*field_, rng));
    }
    for (std::size_t col : pi_cols) {
      row.entries.emplace_back(col, detail::SampleNonZeroElement(*field_, rng));
    }
    return row;
  }

  // Ensures the L x L systematic encoding matrix A_enc is full rank.
  void SelectInvertibleSystematicSalt() {
    std::vector<storage_type> empty_rhs;
    for (uint64_t trial = 0; trial < 128; ++trial) {
      salt_ = trial;
      std::vector<SparseRow> rows;
      rows.reserve(l_);
      for (const auto &r : ldpc_rows_) {
        rows.push_back(r);
      }
      for (const auto &r : hdpc_rows_) {
        rows.push_back(r);
      }
      for (uint32_t esi = 0; esi < static_cast<uint32_t>(k_); ++esi) {
        rows.push_back(BuildLtRow(esi));
      }
      auto res = InactivationSolve(rows, empty_rhs, 0);
      if (res.success) {
        return;
      }
    }
    throw std::runtime_error(
        "RaptorQCode: failed to find invertible systematic precode salt");
  }

  std::shared_ptr<GaloisField> field_;
  std::size_t k_ = 0;
  std::size_t symbol_size_ = 1;
  std::size_t s_ = 0;
  std::size_t h_ = 0;
  std::size_t w_ = 0;
  std::size_t l_ = 0;
  uint64_t seed_ = 0;
  uint64_t salt_ = 0;
  std::vector<SparseRow> ldpc_rows_;
  std::vector<SparseRow> hdpc_rows_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_FOUNTAIN_HPP_
