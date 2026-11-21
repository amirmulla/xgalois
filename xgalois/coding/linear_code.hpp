#ifndef XGALOIS_CODING_LINEAR_CODE_HPP_
#define XGALOIS_CODING_LINEAR_CODE_HPP_

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class LinearCode : public AbstractLinearCode<GaloisField> {
 public:
  using Base = AbstractLinearCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = typename Base::matrix_type;

  static LinearCode FromGeneratorMatrix(
      std::shared_ptr<GaloisField> field, const matrix_type &generator,
      Metric metric = Metric::HAMMING,
      std::optional<std::size_t> known_min_distance = std::nullopt) {
    if (!field) {
      throw std::invalid_argument(
          "LinearCode::FromGeneratorMatrix: field cannot be null");
    }
    if (!generator.HasField() || generator.Field() != field) {
      throw std::invalid_argument(
          "LinearCode::FromGeneratorMatrix: generator matrix field mismatch");
    }
    if (generator.Dimension() != 2 || generator.Shape(1) == 0) {
      throw std::invalid_argument(
          "LinearCode::FromGeneratorMatrix: generator must be a 2D matrix with "
          "positive column count");
    }
    const std::size_t k = generator.Shape(0);
    const std::size_t n = generator.Shape(1);
    if (k > n) {
      throw std::invalid_argument(
          "LinearCode::FromGeneratorMatrix: row count cannot exceed column "
          "count");
    }
    if (k > 0 && linalg::MatrixRank(generator) != k) {
      throw std::invalid_argument(
          "LinearCode::FromGeneratorMatrix: generator matrix must have full "
          "row rank");
    }
    matrix_type H = Base::ComputeParityCheckFromGenerator(generator);
    return LinearCode(std::move(field), generator, std::move(H), metric,
                      known_min_distance);
  }

  static LinearCode FromGeneratorMatrix(
      const matrix_type &generator, Metric metric = Metric::HAMMING,
      std::optional<std::size_t> known_min_distance = std::nullopt) {
    if (!generator.HasField()) {
      throw std::invalid_argument(
          "LinearCode::FromGeneratorMatrix: generator matrix has null field");
    }
    return FromGeneratorMatrix(generator.Field(), generator, metric,
                               known_min_distance);
  }

  static LinearCode FromParityCheckMatrix(
      std::shared_ptr<GaloisField> field, const matrix_type &parity_check,
      Metric metric = Metric::HAMMING,
      std::optional<std::size_t> known_min_distance = std::nullopt) {
    if (!field) {
      throw std::invalid_argument(
          "LinearCode::FromParityCheckMatrix: field cannot be null");
    }
    if (!parity_check.HasField() || parity_check.Field() != field) {
      throw std::invalid_argument(
          "LinearCode::FromParityCheckMatrix: parity-check matrix field "
          "mismatch");
    }
    if (parity_check.Dimension() != 2 || parity_check.Shape(1) == 0) {
      throw std::invalid_argument(
          "LinearCode::FromParityCheckMatrix: parity-check must be a 2D matrix "
          "with positive column count");
    }

    const std::size_t rows = parity_check.Shape(0);
    const std::size_t rank = (rows == 0) ? 0 : linalg::MatrixRank(parity_check);
    matrix_type H_canonical =
        (rank == rows) ? parity_check : linalg::RowSpace(parity_check);
    matrix_type G = Base::ComputeGeneratorFromParityCheck(H_canonical);
    return LinearCode(std::move(field), std::move(G), std::move(H_canonical),
                      metric, known_min_distance);
  }

  static LinearCode FromParityCheckMatrix(
      const matrix_type &parity_check, Metric metric = Metric::HAMMING,
      std::optional<std::size_t> known_min_distance = std::nullopt) {
    if (!parity_check.HasField()) {
      throw std::invalid_argument(
          "LinearCode::FromParityCheckMatrix: parity-check matrix has null "
          "field");
    }
    return FromParityCheckMatrix(parity_check.Field(), parity_check, metric,
                                 known_min_distance);
  }

  const matrix_type &GeneratorMatrix() const override {
    return generator_matrix_;
  }

  const matrix_type &ParityCheckMatrix() const override {
    return parity_check_matrix_;
  }

  matrix_type SystematicGeneratorMatrix() const override {
    return linalg::SystematicForm(generator_matrix_).g_sys;
  }

  LinearCode DualCode() const {
    return LinearCode(this->field_, parity_check_matrix_, generator_matrix_,
                      this->metric_, std::nullopt);
  }

  // Puncturing deletes coordinate col_idx from every codeword:
  // generator matrix is RowSpace(G without column col_idx).
  LinearCode Punctured(std::size_t col_idx) const {
    if (this->length_ <= 1) {
      throw std::invalid_argument(
          "LinearCode::Punctured: cannot puncture code of length 1");
    }
    if (col_idx >= this->length_) {
      throw std::out_of_range(
          "LinearCode::Punctured: column index out of range");
    }
    matrix_type G_drop = DeleteColumn(generator_matrix_, col_idx);
    matrix_type G_basis =
        (G_drop.Shape(0) == 0) ? std::move(G_drop) : linalg::RowSpace(G_drop);
    return FromGeneratorMatrix(this->field_, G_basis, this->metric_);
  }

  // Shortening restricts to codewords with c[col_idx] == 0 and deletes
  // coordinate col_idx; by duality its parity-check matrix is H without
  // column col_idx.
  LinearCode Shortened(std::size_t col_idx) const {
    if (this->length_ <= 1) {
      throw std::invalid_argument(
          "LinearCode::Shortened: cannot shorten code of length 1");
    }
    if (col_idx >= this->length_) {
      throw std::out_of_range(
          "LinearCode::Shortened: column index out of range");
    }
    matrix_type H_drop = DeleteColumn(parity_check_matrix_, col_idx);
    return FromParityCheckMatrix(this->field_, H_drop, this->metric_);
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + "] linear code over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 protected:
  LinearCode(std::shared_ptr<GaloisField> field, matrix_type generator,
             matrix_type parity_check, Metric metric,
             std::optional<std::size_t> known_min_distance)
      : Base(std::move(field), generator.Shape(1), generator.Shape(0), metric,
             "GeneratorMatrix", "Syndrome"),
        generator_matrix_(std::move(generator)),
        parity_check_matrix_(std::move(parity_check)) {
    if (known_min_distance.has_value()) {
      this->cached_min_distance_ = known_min_distance;
    }
    RegisterEncodersAndDecoders();
  }

 private:
  static matrix_type DeleteColumn(const matrix_type &mat, std::size_t col_idx) {
    using RawArray = typename matrix_type::raw_array_type;
    using ShapeType = typename matrix_type::shape_type;
    const std::size_t rows = mat.Shape(0);
    const std::size_t cols = mat.Shape(1);
    const std::size_t new_cols = cols - 1;
    RawArray out = RawArray::from_shape(ShapeType{rows, new_cols});
    const storage_type *src = mat.RawData().data();
    storage_type *dst = out.data();
    for (std::size_t r = 0; r < rows; ++r) {
      const storage_type *src_row = src + r * cols;
      storage_type *dst_row = dst + r * new_cols;
      std::size_t dst_c = 0;
      for (std::size_t c = 0; c < cols; ++c) {
        if (c == col_idx) {
          continue;
        }
        dst_row[dst_c++] = src_row[c];
      }
    }
    return matrix_type::FromRaw(mat.Field(), std::move(out));
  }

  void RegisterEncodersAndDecoders() {
    this->RegisterEncoder(
        "GeneratorMatrix", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<GeneratorMatrixEncoder<GaloisField>>(lin);
        });
    this->RegisterDecoder(
        "Syndrome", [](const AbstractCode<GaloisField> &code) {
          const auto &lin =
              dynamic_cast<const AbstractLinearCode<GaloisField> &>(code);
          return std::make_shared<SyndromeDecoder<GaloisField>>(lin);
        });
  }

  matrix_type generator_matrix_;
  matrix_type parity_check_matrix_;
};

// ============================================================================
// Classic Linear Code Factories
// ============================================================================

// q-ary Hamming code [(q^r - 1)/(q - 1), (q^r - 1)/(q - 1) - r, 3] over GF(q).
// Parity-check columns are all projective representatives in PG(r - 1, q)
// (nonzero vectors whose first nonzero entry is 1).
template <typename GaloisField>
LinearCode<GaloisField> HammingCode(std::shared_ptr<GaloisField> field,
                                    std::size_t r) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;

  if (!field) {
    throw std::invalid_argument("HammingCode: field cannot be null");
  }
  if (r < 2) {
    throw std::invalid_argument(
        "HammingCode: redundancy parameter r must be at least 2");
  }
  const uint64_t q = field->Order();
  if (q == 0) {
    throw std::invalid_argument("HammingCode: field order too large");
  }

  const uint64_t qr = utils::SafeIntegerPower(q, r);
  const std::size_t n = static_cast<std::size_t>((qr - 1) / (q - 1));

  // Enumerate all field elements in internal representation.
  std::vector<StorageType> elems;
  elems.reserve(static_cast<std::size_t>(q));
  if constexpr (std::is_integral_v<StorageType>) {
    for (uint64_t i = 0; i < q; ++i) {
      elems.push_back(field->SetElementValue(static_cast<StorageType>(i)));
    }
  } else {
    for (uint64_t i = 0; i < q; ++i) {
      elems.push_back(field->FromInt(i));
    }
  }

  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();

  RawArray h_raw = RawArray::from_shape(ShapeType{r, n});
  h_raw.fill(zero);
  StorageType *h_data = h_raw.data();

  std::size_t col = 0;
  // Place columns in ascending projective order: first nonzero entry at
  // position `lead` (from r - 1 down to 0 so binary Hamming columns match
  // binary integers 1, 2, ..., 2^r - 1).
  for (std::size_t lead = r; lead-- > 0;) {
    const std::size_t tail = r - 1 - lead;
    const uint64_t tail_count = utils::SafeIntegerPower(q, tail);
    for (uint64_t idx = 0; idx < tail_count; ++idx) {
      h_data[lead * n + col] = one;
      uint64_t tmp = idx;
      for (std::size_t t = r; t-- > lead + 1;) {
        h_data[t * n + col] = elems[static_cast<std::size_t>(tmp % q)];
        tmp /= q;
      }
      ++col;
    }
  }

  auto H = GFMatrix<GaloisField>::FromRaw(field, std::move(h_raw));
  return LinearCode<GaloisField>::FromParityCheckMatrix(std::move(field), H,
                                                        Metric::HAMMING, 3);
}

// Binary Hamming code [2^r - 1, 2^r - 1 - r, 3] over GF(2).
inline LinearCode<GaloisFieldBinary> HammingCode(std::size_t r) {
  return HammingCode(std::make_shared<GaloisFieldBinary>(), r);
}

// Extended binary Hamming code [2^r, 2^r - 1 - r, 4] over GF(2).
inline LinearCode<GaloisFieldBinary> ExtendedHammingCode(
    std::shared_ptr<GaloisFieldBinary> field, std::size_t r) {
  if (!field) {
    throw std::invalid_argument("ExtendedHammingCode: field cannot be null");
  }
  if (r < 2) {
    throw std::invalid_argument(
        "ExtendedHammingCode: redundancy parameter r must be at least 2");
  }
  const uint64_t n64 = utils::SafeIntegerPower(2, r);
  const std::size_t n = static_cast<std::size_t>(n64);
  const std::size_t n_ham = n - 1;

  using RawArray = typename GFMatrix<GaloisFieldBinary>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisFieldBinary>::shape_type;
  RawArray h_raw = RawArray::from_shape(ShapeType{r + 1, n});
  h_raw.fill(uint8_t{0});
  uint8_t *h_data = h_raw.data();

  // First n - 1 columns: nonzero r-bit vectors 1 .. 2^r - 1; last column: 0.
  for (std::size_t col = 0; col < n_ham; ++col) {
    const std::size_t val = col + 1;
    for (std::size_t row = 0; row < r; ++row) {
      h_data[row * n + col] =
          static_cast<uint8_t>((val >> (r - 1 - row)) & 1ULL);
    }
  }
  // Final row: overall even-parity check (all ones).
  for (std::size_t col = 0; col < n; ++col) {
    h_data[r * n + col] = uint8_t{1};
  }

  auto H = GFMatrix<GaloisFieldBinary>::FromRaw(field, std::move(h_raw));
  return LinearCode<GaloisFieldBinary>::FromParityCheckMatrix(
      std::move(field), H, Metric::HAMMING, 4);
}

inline LinearCode<GaloisFieldBinary> ExtendedHammingCode(std::size_t r) {
  return ExtendedHammingCode(std::make_shared<GaloisFieldBinary>(), r);
}

// Binary Golay code [23, 12, 7] over GF(2) generated by
// g(x) = x^11 + x^9 + x^7 + x^6 + x^5 + x + 1.
inline LinearCode<GaloisFieldBinary> BinaryGolayCode(
    std::shared_ptr<GaloisFieldBinary> field) {
  if (!field) {
    throw std::invalid_argument("BinaryGolayCode: field cannot be null");
  }
  constexpr std::size_t k = 12;
  constexpr std::size_t n = 23;
  // Coefficients of g(x) = 1 + x + x^5 + x^6 + x^7 + x^9 + x^11
  constexpr std::array<uint8_t, 12> g = {1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1};

  using RawArray = typename GFMatrix<GaloisFieldBinary>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisFieldBinary>::shape_type;
  RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
  g_raw.fill(uint8_t{0});
  uint8_t *g_data = g_raw.data();
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < g.size(); ++j) {
      g_data[i * n + (i + j)] = g[j];
    }
  }

  auto G = GFMatrix<GaloisFieldBinary>::FromRaw(field, std::move(g_raw));
  return LinearCode<GaloisFieldBinary>::FromGeneratorMatrix(std::move(field), G,
                                                            Metric::HAMMING, 7);
}

inline LinearCode<GaloisFieldBinary> BinaryGolayCode() {
  return BinaryGolayCode(std::make_shared<GaloisFieldBinary>());
}

// Extended binary Golay code [24, 12, 8] over GF(2) (self-dual doubly-even).
inline LinearCode<GaloisFieldBinary> ExtendedBinaryGolayCode(
    std::shared_ptr<GaloisFieldBinary> field) {
  if (!field) {
    throw std::invalid_argument(
        "ExtendedBinaryGolayCode: field cannot be null");
  }
  constexpr std::size_t k = 12;
  constexpr std::size_t n = 24;
  constexpr std::array<uint8_t, 12> g = {1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1};

  using RawArray = typename GFMatrix<GaloisFieldBinary>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisFieldBinary>::shape_type;
  RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
  g_raw.fill(uint8_t{0});
  uint8_t *g_data = g_raw.data();
  for (std::size_t i = 0; i < k; ++i) {
    uint8_t parity = 0;
    for (std::size_t j = 0; j < g.size(); ++j) {
      g_data[i * n + (i + j)] = g[j];
      parity ^= g[j];
    }
    g_data[i * n + 23] = parity;
  }

  auto G = GFMatrix<GaloisFieldBinary>::FromRaw(field, std::move(g_raw));
  return LinearCode<GaloisFieldBinary>::FromGeneratorMatrix(std::move(field), G,
                                                            Metric::HAMMING, 8);
}

inline LinearCode<GaloisFieldBinary> ExtendedBinaryGolayCode() {
  return ExtendedBinaryGolayCode(std::make_shared<GaloisFieldBinary>());
}

// Ternary Golay code [11, 6, 5] over GF(3).
template <typename GaloisField = GaloisFieldPrime<uint8_t>>
LinearCode<GaloisField> TernaryGolayCode(std::shared_ptr<GaloisField> field) {
  if (!field) {
    throw std::invalid_argument("TernaryGolayCode: field cannot be null");
  }
  if (field->Order() != 3) {
    throw std::invalid_argument("TernaryGolayCode: field must be GF(3)");
  }

  constexpr std::size_t k = 6;
  constexpr std::size_t n = 11;
  constexpr std::array<std::array<uint8_t, 5>, 6> A = {{
      {1, 1, 1, 1, 1},
      {0, 1, 2, 2, 1},
      {1, 0, 1, 2, 2},
      {2, 1, 0, 1, 2},
      {2, 2, 1, 0, 1},
      {1, 2, 2, 1, 0},
  }};

  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;

  const StorageType zero = field->AdditiveIdentity();
  const StorageType one = field->MultiplicativeIdentity();
  const StorageType two = field->Add(one, one);
  const std::array<StorageType, 3> vals = {zero, one, two};

  RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
  g_raw.fill(zero);
  StorageType *g_data = g_raw.data();
  for (std::size_t i = 0; i < k; ++i) {
    g_data[i * n + i] = one;
    for (std::size_t j = 0; j < 5; ++j) {
      g_data[i * n + k + j] = vals[A[i][j]];
    }
  }

  auto G = GFMatrix<GaloisField>::FromRaw(field, std::move(g_raw));
  return LinearCode<GaloisField>::FromGeneratorMatrix(std::move(field), G,
                                                      Metric::HAMMING, 5);
}

inline LinearCode<GaloisFieldPrime<uint8_t>> TernaryGolayCode() {
  return TernaryGolayCode(std::make_shared<GaloisFieldPrime<uint8_t>>(3));
}

// Binary Reed-Muller code RM(r, m) over GF(2):
//   length n = 2^m,
//   dimension k = sum_{i=0}^{r} C(m, i),
//   minimum distance d = 2^(m - r).
inline LinearCode<GaloisFieldBinary> ReedMullerCode(
    std::shared_ptr<GaloisFieldBinary> field, std::size_t r, std::size_t m) {
  if (!field) {
    throw std::invalid_argument("ReedMullerCode: field cannot be null");
  }
  if (r > m) {
    throw std::invalid_argument("ReedMullerCode: order r cannot exceed m");
  }
  if (m >= 24) {
    throw std::invalid_argument(
        "ReedMullerCode: m is too large for explicit generator matrix");
  }

  const std::size_t n = std::size_t{1} << m;
  const std::size_t d = std::size_t{1} << (m - r);

  // Collect all subsets of {0, ..., m - 1} of size <= r as bitmasks, grouped
  // by degree 0, 1, ..., r.
  std::vector<std::size_t> monomial_masks;
  auto generate_masks = [&](auto &&self, std::size_t deg, std::size_t start,
                            std::size_t mask) -> void {
    if (deg == 0) {
      monomial_masks.push_back(mask);
      return;
    }
    for (std::size_t bit = start; bit <= m - deg; ++bit) {
      self(self, deg - 1, bit + 1, mask | (std::size_t{1} << bit));
    }
  };
  for (std::size_t deg = 0; deg <= r; ++deg) {
    generate_masks(generate_masks, deg, 0, 0);
  }

  const std::size_t k = monomial_masks.size();
  using RawArray = typename GFMatrix<GaloisFieldBinary>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisFieldBinary>::shape_type;
  RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
  uint8_t *g_data = g_raw.data();

  for (std::size_t row = 0; row < k; ++row) {
    const std::size_t mask = monomial_masks[row];
    uint8_t *g_row = g_data + row * n;
    for (std::size_t col = 0; col < n; ++col) {
      g_row[col] = static_cast<uint8_t>((col & mask) == mask ? 1U : 0U);
    }
  }

  auto G = GFMatrix<GaloisFieldBinary>::FromRaw(field, std::move(g_raw));
  return LinearCode<GaloisFieldBinary>::FromGeneratorMatrix(std::move(field), G,
                                                            Metric::HAMMING, d);
}

inline LinearCode<GaloisFieldBinary> ReedMullerCode(std::size_t r,
                                                    std::size_t m) {
  return ReedMullerCode(std::make_shared<GaloisFieldBinary>(), r, m);
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_LINEAR_CODE_HPP_
