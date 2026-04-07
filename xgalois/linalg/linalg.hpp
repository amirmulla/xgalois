#ifndef XGALOIS_LINALG_LINALG_HPP
#define XGALOIS_LINALG_LINALG_HPP

#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/generators/xrandom.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/misc/xsort.hpp>
#include <xtensor/reducers/xreducer.hpp>
#include <xtensor/views/xbroadcast.hpp>
#include <xtensor/views/xindex_view.hpp>
#include <xtensor/views/xview.hpp>

#include "xgalois/field/gf_element.hpp"

namespace xg {

template <typename GaloisField>
using garray = xt::xarray<GaloisFieldElement<GaloisField>>;

namespace linalg {

template <typename GaloisField>
auto dot(const garray<GaloisField> &a, const garray<GaloisField> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() == 1 && b.dimension() == 1) {
    if (a.shape(0) != b.shape(0)) {
      throw std::invalid_argument(
          "Vector dimensions must match for dot product");
    }

    auto field = a(0).Field();
    ElementType result_val(field->AdditiveIdentity(), field);
    for (std::size_t i = 0; i < a.shape(0); ++i) {
      result_val = result_val + a(i) * b(i);
    }

    std::vector<std::size_t> shape_0d = {};
    garray<GaloisField> result(shape_0d);
    result() = result_val;
    return result;
  }

  if (a.dimension() == 2 && b.dimension() == 1) {
    if (a.shape(1) != b.shape(0)) {
      throw std::invalid_argument("Matrix columns must match vector size");
    }

    auto field = a(0, 0).Field();
    garray<GaloisField> result({a.shape(0)});

    for (std::size_t i = 0; i < a.shape(0); ++i) {
      ElementType sum(field->AdditiveIdentity(), field);
      for (std::size_t j = 0; j < a.shape(1); ++j) {
        sum = sum + a(i, j) * b(j);
      }
      result(i) = sum;
    }
    return result;
  }

  if (a.dimension() == 2 && b.dimension() == 2) {
    if (a.shape(1) != b.shape(0)) {
      throw std::invalid_argument(
          "Matrix dimensions incompatible for multiplication");
    }

    auto field = a(0, 0).Field();
    garray<GaloisField> result({a.shape(0), b.shape(1)});

    for (std::size_t i = 0; i < a.shape(0); ++i) {
      for (std::size_t j = 0; j < b.shape(1); ++j) {
        ElementType sum(field->AdditiveIdentity(), field);
        for (std::size_t k = 0; k < a.shape(1); ++k) {
          sum = sum + a(i, k) * b(k, j);
        }
        result(i, j) = sum;
      }
    }
    return result;
  }

  throw std::invalid_argument("Unsupported array dimensions for dot product");
}

template <typename GaloisField>
auto vdot(const garray<GaloisField> &a, const garray<GaloisField> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.size() != b.size()) {
    throw std::invalid_argument(
        "Arrays must have same number of elements for vdot");
  }

  auto a_flat = xt::flatten(a);
  auto b_flat = xt::flatten(b);

  auto field = a_flat(0).Field();
  ElementType result(field->AdditiveIdentity(), field);
  for (std::size_t i = 0; i < a_flat.size(); ++i) {
    result = result + a_flat(i) * b_flat(i);
  }
  return result;
}

template <typename GaloisField>
auto outer(const garray<GaloisField> &a, const garray<GaloisField> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 1 || b.dimension() != 1) {
    throw std::invalid_argument("outer product requires 1-D arrays");
  }

  garray<GaloisField> result({a.shape(0), b.shape(0)});

  for (std::size_t i = 0; i < a.shape(0); ++i) {
    for (std::size_t j = 0; j < b.shape(0); ++j) {
      result(i, j) = a(i) * b(j);
    }
  }
  return result;
}

template <typename GaloisField>
auto matrix_power(const garray<GaloisField> &a, int n) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("matrix_power requires a square matrix");
  }

  std::size_t size = a.shape(0);
  auto field = a(0, 0).Field();

  if (n == 0) {

    auto field = a(0, 0).Field();
    garray<GaloisField> result({size, size});

    ElementType zero(field->AdditiveIdentity(), field);
    std::fill(result.begin(), result.end(), zero);

    for (std::size_t i = 0; i < size; ++i) {
      result(i, i) = ElementType(field->MultiplicativeIdentity(), field);
    }
    return result;
  }

  if (n == 1) {
    return a;
  }

  if (n < 0) {

    throw std::invalid_argument("Negative matrix powers not yet implemented");
  }

  auto result = a;
  auto base = a;
  n--;

  while (n > 0) {
    if (n & 1) {
      result = dot(result, base);
    }
    base = dot(base, base);
    n >>= 1;
  }

  return result;
}

template <typename GaloisField>
auto kron(const garray<GaloisField> &a, const garray<GaloisField> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || b.dimension() != 2) {
    throw std::invalid_argument("kron requires 2-D arrays (matrices)");
  }

  std::size_t a_rows = a.shape(0);
  std::size_t a_cols = a.shape(1);
  std::size_t b_rows = b.shape(0);
  std::size_t b_cols = b.shape(1);

  garray<GaloisField> result({a_rows * b_rows, a_cols * b_cols});
  auto field = a(0, 0).Field();
  ElementType zero(field->AdditiveIdentity(), field);

  for (std::size_t i = 0; i < a_rows; ++i) {
    for (std::size_t j = 0; j < a_cols; ++j) {
      for (std::size_t k = 0; k < b_rows; ++k) {
        for (std::size_t l = 0; l < b_cols; ++l) {
          result(i * b_rows + k, j * b_cols + l) = a(i, j) * b(k, l);
        }
      }
    }
  }
  return result;
}

template <typename GaloisField>
auto row_echelon(const garray<GaloisField> &matrix_input) {
  using ElementType = GaloisFieldElement<GaloisField>;
  auto R = matrix_input;

  if (R.dimension() != 2) {
    throw std::invalid_argument("row_echelon requires a 2-D array (matrix)");
  }
  if (R.shape(0) == 0 || R.shape(1) == 0) {
    return R;
  }

  auto field =
      R(0, 0).Field();
  std::size_t rows = R.shape(0);
  std::size_t cols = R.shape(1);
  std::size_t lead = 0;

  for (std::size_t r = 0; r < rows && lead < cols; ++r) {
    std::size_t i = r;

    while (lead < cols) {
      i = r;
      while (i < rows && R(i, lead).Value() == field->AdditiveIdentity()) {
        i++;
      }
      if (i <
          rows) {
        break;
      }

      lead++;
    }

    if (lead == cols) {
      break;
    }

    if (i != r) {
      for (std::size_t j = 0; j < cols; ++j) {
        std::swap(R(r, j), R(i, j));
      }
    }

    ElementType pivot_val = R(r, lead);

    if (pivot_val.Value() !=
        field->MultiplicativeIdentity()) {
      for (std::size_t j = lead; j < cols; ++j) {
        R(r, j) = R(r, j) / pivot_val;
      }
    }

    for (std::size_t k = r + 1; k < rows; ++k) {
      if (R(k, lead).Value() !=
          field->AdditiveIdentity()) {
        ElementType factor = R(k, lead);
        for (std::size_t j = lead; j < cols; ++j) {
          R(k, j) = R(k, j) - factor * R(r, j);
        }
      }
    }
    lead++;
  }
  return R;
}

template <typename GaloisField>
auto rref(const garray<GaloisField> &matrix_input) {
  using ElementType = GaloisFieldElement<GaloisField>;
  auto R = row_echelon(matrix_input);

  if (R.dimension() != 2) {
    throw std::invalid_argument("rref requires a 2-D array (matrix)");
  }
  if (R.shape(0) == 0 || R.shape(1) == 0) {
    return R;
  }

  auto field =
      R(0, 0).Field();
  std::size_t rows = R.shape(0);
  std::size_t cols = R.shape(1);

  for (int r = static_cast<int>(rows) - 1; r >= 0; --r) {

    std::size_t pivot_col = cols;
    for (std::size_t j = 0; j < cols; ++j) {
      if (R(r, j).Value() == field->MultiplicativeIdentity()) {
        pivot_col = j;
        break;
      }

      if (R(r, j).Value() != field->AdditiveIdentity()) {
        break;
      }
    }

    if (pivot_col < cols) {

      for (int i = r - 1; i >= 0; --i) {
        if (R(i, pivot_col).Value() != field->AdditiveIdentity()) {
          ElementType factor = R(i, pivot_col);
          for (std::size_t j = pivot_col; j < cols; ++j) {
            R(i, j) = R(i, j) - factor * R(r, j);
          }
        }
      }
    }
  }
  return R;
}

template <typename GaloisField>
auto trace(const garray<GaloisField> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("trace requires a square matrix");
  }

  std::size_t size = a.shape(0);
  auto field = a(0, 0).Field();
  ElementType result(field->AdditiveIdentity(), field);

  for (std::size_t i = 0; i < size; ++i) {
    result = result + a(i, i);
  }

  return result;
}

template <typename GaloisField>
auto det(const garray<GaloisField> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("determinant requires a square matrix");
  }

  std::size_t n = a.shape(0);
  auto field = a(0, 0).Field();

  auto matrix = a;
  ElementType det_result(field->MultiplicativeIdentity(), field);

  for (std::size_t i = 0; i < n; ++i) {

    std::size_t pivot_row = i;
    for (std::size_t k = i + 1; k < n; ++k) {
      if (matrix(k, i).Value() != field->AdditiveIdentity()) {
        pivot_row = k;
        break;
      }
    }

    if (pivot_row != i) {
      for (std::size_t j = 0; j < n; ++j) {
        std::swap(matrix(i, j), matrix(pivot_row, j));
      }
      det_result =
          det_result *
          ElementType(field->Neg(field->MultiplicativeIdentity()), field);
    }

    if (matrix(i, i).Value() == field->AdditiveIdentity()) {
      return ElementType(field->AdditiveIdentity(), field);
    }

    det_result = det_result * matrix(i, i);

    for (std::size_t k = i + 1; k < n; ++k) {
      if (matrix(k, i).Value() != field->AdditiveIdentity()) {
        auto factor = matrix(k, i) / matrix(i, i);
        for (std::size_t j = i; j < n; ++j) {
          matrix(k, j) = matrix(k, j) - factor * matrix(i, j);
        }
      }
    }
  }

  return det_result;
}

template <typename GaloisField>
std::size_t matrix_rank(const garray<GaloisField> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2) {
    throw std::invalid_argument("matrix_rank requires a 2-D array (matrix)");
  }

  if (a.shape(0) == 0 || a.shape(1) == 0) {
    return 0;
  }

  auto field = a(0, 0).Field();
  garray<GaloisField> R = row_echelon(a);

  std::size_t rank = 0;
  for (std::size_t i = 0; i < R.shape(0); ++i) {
    bool non_zero_row = false;
    for (std::size_t j = 0; j < R.shape(1); ++j) {
      if (R(i, j).Value() != field->AdditiveIdentity()) {
        non_zero_row = true;
        break;
      }
    }
    if (non_zero_row) {
      rank++;
    }
  }
  return rank;
}

template <typename GaloisField>
auto solve(const garray<GaloisField> &a, const garray<GaloisField> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("solve requires a square coefficient matrix");
  }

  if (b.dimension() != 1 || b.shape(0) != a.shape(0)) {
    throw std::invalid_argument(
        "Right-hand side must be a vector with size matching matrix rows");
  }

  std::size_t n = a.shape(0);
  auto field = a(0, 0).Field();

  garray<GaloisField> augmented({n, n + 1});

  ElementType zero(field->AdditiveIdentity(), field);
  std::fill(augmented.begin(), augmented.end(), zero);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      augmented(i, j) = a(i, j);
    }
    augmented(i, n) = b(i);
  }

  for (std::size_t i = 0; i < n; ++i) {

    std::size_t pivot_row = i;
    for (std::size_t k = i + 1; k < n; ++k) {
      if (augmented(k, i).Value() != field->AdditiveIdentity()) {
        pivot_row = k;
        break;
      }
    }

    if (pivot_row != i) {
      for (std::size_t j = 0; j <= n; ++j) {
        std::swap(augmented(i, j), augmented(pivot_row, j));
      }
    }

    if (augmented(i, i).Value() == field->AdditiveIdentity()) {
      throw std::runtime_error("Matrix is singular");
    }

    for (std::size_t k = i + 1; k < n; ++k) {
      if (augmented(k, i).Value() != field->AdditiveIdentity()) {
        auto factor = augmented(k, i) / augmented(i, i);
        for (std::size_t j = i; j <= n; ++j) {
          augmented(k, j) = augmented(k, j) - factor * augmented(i, j);
        }
      }
    }
  }

  garray<GaloisField> x({n});

  std::fill(x.begin(), x.end(), zero);

  for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
    ElementType sum(field->AdditiveIdentity(), field);
    for (std::size_t j = i + 1; j < n; ++j) {
      sum = sum + augmented(i, j) * x(j);
    }
    x(i) = (augmented(i, n) - sum) / augmented(i, i);
  }

  return x;
}

template <typename GaloisField>
auto inv(const garray<GaloisField> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("inv requires a square matrix");
  }

  std::size_t n = a.shape(0);
  auto field = a(0, 0).Field();

  garray<GaloisField> augmented({n, 2 * n});

  ElementType zero(field->AdditiveIdentity(), field);
  std::fill(augmented.begin(), augmented.end(), zero);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      augmented(i, j) = a(i, j);
      augmented(i, j + n) =
          (i == j) ? ElementType(field->MultiplicativeIdentity(), field)
                   : ElementType(field->AdditiveIdentity(), field);
    }
  }

  for (std::size_t i = 0; i < n; ++i) {

    std::size_t pivot_row = i;
    for (std::size_t k = i + 1; k < n; ++k) {
      if (augmented(k, i).Value() != field->AdditiveIdentity()) {
        pivot_row = k;
        break;
      }
    }

    if (pivot_row != i) {
      for (std::size_t j = 0; j < 2 * n; ++j) {
        std::swap(augmented(i, j), augmented(pivot_row, j));
      }
    }

    if (augmented(i, i).Value() == field->AdditiveIdentity()) {
      throw std::runtime_error("Matrix is singular");
    }

    auto pivot = augmented(i, i);
    for (std::size_t j = 0; j < 2 * n; ++j) {
      augmented(i, j) = augmented(i, j) / pivot;
    }

    for (std::size_t k = 0; k < n; ++k) {
      if (k != i && augmented(k, i).Value() != field->AdditiveIdentity()) {
        auto factor = augmented(k, i);
        for (std::size_t j = 0; j < 2 * n; ++j) {
          augmented(k, j) = augmented(k, j) - factor * augmented(i, j);
        }
      }
    }
  }

  garray<GaloisField> result({n, n});

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      result(i, j) = augmented(i, j + n);
    }
  }

  return result;
}

template <typename GaloisField>
auto cross(const garray<GaloisField> &a, const garray<GaloisField> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 1 || b.dimension() != 1) {
    throw std::invalid_argument("cross product requires 1-D arrays (vectors)");
  }
  if ((a.shape(0) != 2 && a.shape(0) != 3) ||
      (b.shape(0) != 2 && b.shape(0) != 3)) {
    throw std::invalid_argument(
        "cross product vectors must have 2 or 3 elements");
  }

  auto field = a(0).Field();
  ElementType zero(field->AdditiveIdentity(), field);

  const ElementType &a0 = a(0);
  const ElementType &a1 = a(1);
  ElementType a2 = (a.shape(0) == 3) ? a(2) : zero;

  const ElementType &b0 = b(0);
  const ElementType &b1 = b(1);
  ElementType b2 = (b.shape(0) == 3) ? b(2) : zero;

  garray<GaloisField> result({3});
  result(0) = a1 * b2 - a2 * b1;
  result(1) = a2 * b0 - a0 * b2;
  result(2) = a0 * b1 - a1 * b0;

  return result;
}

template <typename GaloisField>
auto eye(std::size_t n, const std::shared_ptr<GaloisField> &field) {
  using ElementType = GaloisFieldElement<GaloisField>;

  garray<GaloisField> result({n, n});

  ElementType zero(field->AdditiveIdentity(), field);
  std::fill(result.begin(), result.end(), zero);

  for (std::size_t i = 0; i < n; ++i) {
    result(i, i) = ElementType(field->MultiplicativeIdentity(), field);
  }

  return result;
}

template <typename GaloisField>
auto zeros(const std::vector<std::size_t> &shape,
           const std::shared_ptr<GaloisField> &field) {
  using ElementType = GaloisFieldElement<GaloisField>;

  garray<GaloisField> result(shape);

  auto zero = ElementType(field->AdditiveIdentity(), field);
  std::fill(result.begin(), result.end(), zero);

  return result;
}

}
}

#endif
