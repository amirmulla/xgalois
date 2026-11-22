#ifndef XGALOIS_CODING_SUBSPACE_CODE_HPP_
#define XGALOIS_CODING_SUBSPACE_CODE_HPP_

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/coding/gabidulin.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

// ============================================================================
// SubspaceDecodingResult: Result of Kotter-Kschischang Subspace Decoding
// ============================================================================

template <typename BaseField, typename ExtField>
struct SubspaceDecodingResult {
  GFMatrix<BaseField>
      subspace_basis;  // Canonical RREF basis of shape {n, n+m}.
  GFArray<ExtField> gabidulin_codeword;
  GFArray<ExtField> message;
  std::size_t subspace_distance = 0;
  std::size_t num_erasures = 0;    // rho + gamma (deviations + missing pivots).
  std::size_t num_insertions = 0;  // tau (unknown rank errors).
  bool success = false;
  std::string error_message;

  bool Ok() const noexcept { return success; }
  explicit operator bool() const noexcept { return success; }

  static SubspaceDecodingResult Success(GFMatrix<BaseField> basis,
                                        GFArray<ExtField> codeword,
                                        GFArray<ExtField> msg,
                                        std::size_t distance,
                                        std::size_t erasures,
                                        std::size_t insertions) {
    SubspaceDecodingResult res;
    res.subspace_basis = std::move(basis);
    res.gabidulin_codeword = std::move(codeword);
    res.message = std::move(msg);
    res.subspace_distance = distance;
    res.num_erasures = erasures;
    res.num_insertions = insertions;
    res.success = true;
    return res;
  }

  static SubspaceDecodingResult Failure(std::string reason) {
    SubspaceDecodingResult res;
    res.success = false;
    res.error_message = std::move(reason);
    return res;
  }
};

// ============================================================================
// KotterKschischangCode: Lifted Gabidulin Constant-Dimension Subspace Code
//   I(G(n, k)) subset G_q(n + m, n), minimum subspace distance d_S = 2(n-k+1)
// ============================================================================

template <typename BaseField, typename ExtField>
class KotterKschischangCode {
 public:
  using base_field_type = BaseField;
  using ext_field_type = ExtField;
  using base_storage_type = typename BaseField::element_type;
  using ext_storage_type = typename ExtField::element_type;
  using base_matrix_type = GFMatrix<BaseField>;
  using codeword_type = GFArray<ExtField>;
  using message_type = GFArray<ExtField>;
  using gabidulin_type = GabidulinCode<BaseField, ExtField>;
  using result_type = SubspaceDecodingResult<BaseField, ExtField>;

  explicit KotterKschischangCode(gabidulin_type gabidulin)
      : gabidulin_(std::move(gabidulin)) {}

  KotterKschischangCode(std::shared_ptr<BaseField> base_field,
                        std::shared_ptr<ExtField> ext_field,
                        const GFArray<ExtField> &evaluation_points,
                        std::size_t dimension)
      : gabidulin_(std::move(base_field), std::move(ext_field),
                   evaluation_points, dimension) {}

  static KotterKschischangCode CreateStandard(
      std::shared_ptr<BaseField> base_field,
      std::shared_ptr<ExtField> ext_field, std::size_t length,
      std::size_t dimension) {
    return KotterKschischangCode(gabidulin_type::CreateStandard(
        std::move(base_field), std::move(ext_field), length, dimension));
  }

  std::size_t AmbientDimension() const noexcept {
    return gabidulin_.Length() + gabidulin_.Tower().RelativeDegree();
  }

  std::size_t SubspaceDimension() const noexcept { return gabidulin_.Length(); }

  std::size_t CodeDimension() const noexcept { return gabidulin_.Dimension(); }

  std::size_t ExtensionDegree() const noexcept {
    return gabidulin_.Tower().RelativeDegree();
  }

  std::size_t MinimumSubspaceDistance() const noexcept {
    return 2 * gabidulin_.MinimumDistance();
  }

  Metric GetMetric() const noexcept { return Metric::SUBSPACE; }

  std::shared_ptr<BaseField> BaseFieldPtr() const {
    return gabidulin_.BaseFieldPtr();
  }

  std::shared_ptr<ExtField> ExtFieldPtr() const { return gabidulin_.Field(); }

  const gabidulin_type &UnderlyingGabidulinCode() const noexcept {
    return gabidulin_;
  }

  // Lifts a Gabidulin codeword c in F_{q^m}^n to the n-dimensional subspace
  // spanned by [I_n | ExpandToBaseMatrix(c)^T] of shape {n, n + m}.
  base_matrix_type LiftCodeword(const codeword_type &codeword) const {
    const std::size_t n = gabidulin_.Length();
    const base_matrix_type expanded = gabidulin_.ExpandToBaseMatrix(codeword);
    const base_matrix_type eye_n = linalg::Eye(gabidulin_.BaseFieldPtr(), n);
    return linalg::HStack(eye_n, linalg::Transpose(expanded));
  }

  // Encodes message m in F_{q^m}^k into canonical subspace basis [I_n | C^T].
  base_matrix_type Encode(const message_type &message) const {
    const codeword_type cw = gabidulin_.Encode(message);
    return LiftCodeword(cw);
  }

  // Canonical RREF basis of the row space of `basis` (non-zero rows only).
  static base_matrix_type CanonicalBasis(const base_matrix_type &basis) {
    if (!basis.HasField() || basis.Dimension() != 2) {
      throw std::invalid_argument(
          "KotterKschischangCode::CanonicalBasis: input must be a 2D matrix");
    }
    return linalg::RowSpace(basis);
  }

  // Exact subspace distance d_S(U, V) = 2 * rank([U; V]) - rank(U) - rank(V).
  static std::size_t SubspaceDistance(const base_matrix_type &u_basis,
                                      const base_matrix_type &v_basis) {
    if (!u_basis.HasField() || !v_basis.HasField() ||
        u_basis.Field() != v_basis.Field()) {
      throw std::invalid_argument(
          "KotterKschischangCode::SubspaceDistance: field mismatch");
    }
    if (u_basis.Dimension() != 2 || v_basis.Dimension() != 2) {
      throw std::invalid_argument(
          "KotterKschischangCode::SubspaceDistance: inputs must be 2D "
          "matrices");
    }
    if (u_basis.Shape(1) != v_basis.Shape(1)) {
      throw std::invalid_argument(
          "KotterKschischangCode::SubspaceDistance: ambient dimensions must "
          "match");
    }
    const std::size_t rank_u = linalg::MatrixRank(u_basis);
    const std::size_t rank_v = linalg::MatrixRank(v_basis);
    const std::size_t rank_sum =
        linalg::MatrixRank(linalg::VStack(u_basis, v_basis));
    return 2 * rank_sum - rank_u - rank_v;
  }

  // Silva-Kotter-Kschischang (2008) Operator-Channel Subspace Decoder.
  // Uniquely recovers U = I(c) whenever d_S(U, Y) = 2*tau + rho + gamma <= n-k.
  result_type Decode(const base_matrix_type &received_subspace) const {
    auto base_field = gabidulin_.BaseFieldPtr();
    auto ext_field = gabidulin_.Field();
    const std::size_t n = gabidulin_.Length();
    const std::size_t k = gabidulin_.Dimension();
    const std::size_t m = gabidulin_.Tower().RelativeDegree();
    const std::size_t ambient = n + m;
    const std::size_t r = n - k;

    if (!received_subspace.HasField() ||
        received_subspace.Field() != base_field ||
        received_subspace.Dimension() != 2 ||
        received_subspace.Shape(1) != ambient) {
      return result_type::Failure(
          "KotterKschischangCode::Decode: received subspace must be a 2D "
          "matrix over BaseField with n + m columns");
    }

    // 1. Canonical RREF of Y = [A | B] (non-zero rows only).
    const base_matrix_type y_rref = CanonicalBasis(received_subspace);
    const std::size_t ell = y_rref.Shape(0);

    // 2. Partition rows by whether their leading pivot lies in [0, n) or [n,
    // n+m).
    std::vector<std::size_t> left_row_indices;
    std::vector<std::size_t> left_pivot_cols;
    std::vector<std::size_t> right_row_indices;
    std::vector<bool> has_left_pivot(n, false);

    left_row_indices.reserve(n);
    left_pivot_cols.reserve(n);

    for (std::size_t row = 0; row < ell; ++row) {
      std::size_t pivot = ambient;
      for (std::size_t col = 0; col < ambient; ++col) {
        if (!base_field->IsZero(y_rref.Raw(row, col))) {
          pivot = col;
          break;
        }
      }
      if (pivot < n) {
        left_row_indices.push_back(row);
        left_pivot_cols.push_back(pivot);
        has_left_pivot[pivot] = true;
      } else if (pivot < ambient) {
        right_row_indices.push_back(row);
      }
    }

    std::vector<std::size_t> missing_pivots;
    for (std::size_t col = 0; col < n; ++col) {
      if (!has_left_pivot[col]) {
        missing_pivots.push_back(col);
      }
    }

    const std::size_t gamma = missing_pivots.size();
    const std::size_t rho = right_row_indices.size();
    if (rho + gamma > r) {
      return result_type::Failure(
          "KotterKschischangCode::Decode: subspace erasures and deviations "
          "exceed unique decoding bound n - k");
    }

    // 3. Construct received Gabidulin matrix R (shape {m, n}):
    //    For each row i with left pivot j = p_i < n, column j of R is row i's
    //    right block B[i, :]. Missing pivot columns remain zero.
    base_matrix_type r_mat = base_matrix_type::Zeros(base_field, {m, n});
    for (std::size_t idx = 0; idx < left_row_indices.size(); ++idx) {
      const std::size_t row = left_row_indices[idx];
      const std::size_t j = left_pivot_cols[idx];
      for (std::size_t t = 0; t < m; ++t) {
        r_mat.Raw(t, j) = y_rref.Raw(row, n + t);
      }
    }
    const codeword_type received_word =
        gabidulin_.ContractFromBaseMatrix(r_mat);

    // 4. Construct column-erasure matrix B^{(C)} (shape {gamma, n}):
    //    Exact nullspace basis of the left block A_top: for each missing pivot
    //    column e, entry e is 1 and entry j = p_i is -Y_rref[row_i, e].
    base_matrix_type col_erasure_mat =
        base_matrix_type::Zeros(base_field, {gamma, n});
    const base_storage_type base_one = base_field->MultiplicativeIdentity();
    for (std::size_t c_idx = 0; c_idx < gamma; ++c_idx) {
      const std::size_t e = missing_pivots[c_idx];
      col_erasure_mat.Raw(c_idx, e) = base_one;
      for (std::size_t idx = 0; idx < left_row_indices.size(); ++idx) {
        const std::size_t row = left_row_indices[idx];
        const std::size_t j = left_pivot_cols[idx];
        const base_storage_type entry = y_rref.Raw(row, e);
        if (!base_field->IsZero(entry)) {
          col_erasure_mat.Raw(c_idx, j) = base_field->Neg(entry);
        }
      }
    }

    // 5. Extract row-erasure values a^{(R)} (length rho in F_{q^m}) from rows
    //    whose pivots lie in [n, n + m).
    codeword_type row_erasure_vals = codeword_type::Zeros(ext_field, {0});
    if (rho > 0) {
      base_matrix_type a_row_mat =
          base_matrix_type::Zeros(base_field, {m, rho});
      for (std::size_t u = 0; u < rho; ++u) {
        const std::size_t row = right_row_indices[u];
        for (std::size_t t = 0; t < m; ++t) {
          a_row_mat.Raw(t, u) = y_rref.Raw(row, n + t);
        }
      }
      row_erasure_vals = gabidulin_.ContractFromBaseMatrix(a_row_mat);
    }

    // 6. Invoke SKK rank errors + row/column erasures decoder on Gabidulin
    // code.
    auto gab_res = gabidulin_.DecodeWithRankErasures(
        received_word, row_erasure_vals, col_erasure_mat);
    if (!gab_res.Ok()) {
      return result_type::Failure(gab_res.error_message);
    }

    // 7. Lift recovered codeword and verify exact subspace distance bound.
    base_matrix_type hat_u = LiftCodeword(gab_res.codeword);
    const std::size_t dist = SubspaceDistance(hat_u, y_rref);
    if (dist > r) {
      return result_type::Failure(
          "KotterKschischangCode::Decode: recovered subspace exceeds unique "
          "decoding radius");
    }

    return result_type::Success(std::move(hat_u), std::move(gab_res.codeword),
                                std::move(gab_res.message), dist, rho + gamma,
                                gab_res.num_errors);
  }

  std::string ToString() const {
    return "KotterKschischangCode G_q(" + std::to_string(AmbientDimension()) +
           ", " + std::to_string(SubspaceDimension()) +
           ", d_S=" + std::to_string(MinimumSubspaceDistance()) + ") lifting " +
           gabidulin_.ToString();
  }

 private:
  gabidulin_type gabidulin_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_SUBSPACE_CODE_HPP_
