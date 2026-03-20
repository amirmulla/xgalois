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

#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_element.hpp"

namespace xg {

template <typename GaloisField>
using garray = xt::xarray<GaloisFieldElement<GaloisField>>;

namespace linalg {

//------------------------------------------------------------------------------
// Matrix and Vector Products
//------------------------------------------------------------------------------

/**
 * Compute the dot product of two arrays over a finite field.
 * For 1-D arrays, this computes the scalar dot product.
 * For 2-D arrays, this is matrix multiplication.
 *
 * @param a First input array
 * @param b Second input array
 * @return The dot product result
 */
template <typename GaloisField>
auto dot(const xt::xarray<GaloisFieldElement<GaloisField>> &a,
         const xt::xarray<GaloisFieldElement<GaloisField>> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  // Handle vector dot product (1-D x 1-D) - return as 0-D array
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
    // Return as a 0-D array for consistency
    std::vector<std::size_t> shape_0d = {};
    xt::xarray<ElementType> result(shape_0d);
    result() = result_val;
    return result;
  }

  // Handle matrix-vector multiplication (2-D x 1-D)
  if (a.dimension() == 2 && b.dimension() == 1) {
    if (a.shape(1) != b.shape(0)) {
      throw std::invalid_argument("Matrix columns must match vector size");
    }

    auto field = a(0, 0).Field();
    xt::xarray<ElementType> result({a.shape(0)});

    for (std::size_t i = 0; i < a.shape(0); ++i) {
      ElementType sum(field->AdditiveIdentity(), field);
      for (std::size_t j = 0; j < a.shape(1); ++j) {
        sum = sum + a(i, j) * b(j);
      }
      result(i) = sum;
    }
    return result;
  }

  // Handle matrix multiplication (2-D x 2-D)
  if (a.dimension() == 2 && b.dimension() == 2) {
    if (a.shape(1) != b.shape(0)) {
      throw std::invalid_argument(
          "Matrix dimensions incompatible for multiplication");
    }

    auto field = a(0, 0).Field();
    xt::xarray<ElementType> result({a.shape(0), b.shape(1)});

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

/**
 * Compute the vector dot product (flattened arrays)
 *
 * @param a First input array
 * @param b Second input array
 * @return The vector dot product result
 */
template <typename GaloisField>
auto vdot(const xt::xarray<GaloisFieldElement<GaloisField>> &a,
          const xt::xarray<GaloisFieldElement<GaloisField>> &b) {
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

/**
 * Compute the outer product of two vectors
 *
 * @param a First input vector
 * @param b Second input vector
 * @return The outer product matrix
 */
template <typename GaloisField>
auto outer(const xt::xarray<GaloisFieldElement<GaloisField>> &a,
           const xt::xarray<GaloisFieldElement<GaloisField>> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 1 || b.dimension() != 1) {
    throw std::invalid_argument("outer product requires 1-D arrays");
  }

  xt::xarray<ElementType> result({a.shape(0), b.shape(0)});

  for (std::size_t i = 0; i < a.shape(0); ++i) {
    for (std::size_t j = 0; j < b.shape(0); ++j) {
      result(i, j) = a(i) * b(j);
    }
  }
  return result;
}

/**
 * Raise a square matrix to a power
 *
 * @param a Input square matrix
 * @param n Power to raise the matrix to
 * @return The matrix power result
 */
template <typename GaloisField>
auto matrix_power(const xt::xarray<GaloisFieldElement<GaloisField>> &a, int n) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("matrix_power requires a square matrix");
  }

  std::size_t size = a.shape(0);
  auto field = a(0, 0).Field();

  if (n == 0) {
    // Return identity matrix
    auto field = a(0, 0).Field();
    xt::xarray<ElementType> result({size, size});

    // Initialize to zero first
    ElementType zero(field->AdditiveIdentity(), field);
    std::fill(result.begin(), result.end(), zero);

    // Set diagonal to identity
    for (std::size_t i = 0; i < size; ++i) {
      result(i, i) = ElementType(field->MultiplicativeIdentity(), field);
    }
    return result;
  }

  if (n == 1) {
    return a;
  }

  if (n < 0) {
    // For negative powers, we need the inverse
    throw std::invalid_argument("Negative matrix powers not yet implemented");
  }

  // Use repeated squaring for positive powers
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

/**
 * Compute the Kronecker product of two matrices.
 *
 * @param a First input matrix
 * @param b Second input matrix
 * @return The Kronecker product matrix
 */
template <typename GaloisField>
auto kron(const xt::xarray<GaloisFieldElement<GaloisField>> &a,
          const xt::xarray<GaloisFieldElement<GaloisField>> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || b.dimension() != 2) {
    throw std::invalid_argument("kron requires 2-D arrays (matrices)");
  }

  std::size_t a_rows = a.shape(0);
  std::size_t a_cols = a.shape(1);
  std::size_t b_rows = b.shape(0);
  std::size_t b_cols = b.shape(1);

  xt::xarray<ElementType> result({a_rows * b_rows, a_cols * b_cols});
  auto field = a(0, 0).Field();  // Assuming a is not empty
  ElementType zero(field->AdditiveIdentity(), field);

  // Initialize result with zeros if needed, though direct assignment will
  // overwrite std::fill(result.begin(), result.end(), zero); // Optional: if
  // elements might not be written

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

//------------------------------------------------------------------------------
// Matrix Properties and Norms
//------------------------------------------------------------------------------

/**
 * Compute the row echelon form of a matrix over a finite field.
 *
 * @param matrix_input Input matrix
 * @return The row echelon form of the matrix
 */
template <typename GaloisField>
auto row_echelon(
    const xt::xarray<GaloisFieldElement<GaloisField>> &matrix_input) {
  using ElementType = GaloisFieldElement<GaloisField>;
  auto R = matrix_input;  // Make a copy to modify

  if (R.dimension() != 2) {
    throw std::invalid_argument("row_echelon requires a 2-D array (matrix)");
  }
  if (R.shape(0) == 0 || R.shape(1) == 0) {
    return R;  // Return empty or original if no rows/cols
  }

  auto field =
      R(0, 0).Field();  // Assuming R is not empty and has at least one element
  std::size_t rows = R.shape(0);
  std::size_t cols = R.shape(1);
  std::size_t lead = 0;  // Current leading column for pivot search

  for (std::size_t r = 0; r < rows && lead < cols; ++r) {
    std::size_t i = r;
    // Find pivot for this row r in column lead
    while (lead < cols) {
      i = r;
      while (i < rows && R(i, lead).Value() == field->AdditiveIdentity()) {
        i++;
      }
      if (i <
          rows) {  // Found a non-zero pivot in this column at or below row r
        break;
      }
      // No pivot in this column at or below r, move to next column
      lead++;
    }

    if (lead == cols) {  // All remaining columns are zero for rows >= r
      break;
    }

    // Swap row i (pivot row) with row r
    if (i != r) {
      for (std::size_t j = 0; j < cols; ++j) {
        std::swap(R(r, j), R(i, j));
      }
    }

    // Normalize pivot row r (make R(r, lead) = 1)
    ElementType pivot_val = R(r, lead);
    // The loop structure ensures pivot_val is not zero here.
    if (pivot_val.Value() !=
        field->MultiplicativeIdentity()) {  // Only divide if not already 1
      for (std::size_t j = lead; j < cols; ++j) {
        R(r, j) = R(r, j) / pivot_val;
      }
    }

    // Eliminate elements below the pivot row r in the current lead column
    for (std::size_t k = r + 1; k < rows; ++k) {
      if (R(k, lead).Value() !=
          field->AdditiveIdentity()) {    // If element below pivot is non-zero
        ElementType factor = R(k, lead);  // Since R(r, lead) is now 1
        for (std::size_t j = lead; j < cols; ++j) {
          R(k, j) = R(k, j) - factor * R(r, j);
        }
      }
    }
    lead++;  // Move to the next column for the next pivot
  }
  return R;
}

/**
 * Compute the reduced row echelon form of a matrix over a finite field.
 *
 * @param matrix_input Input matrix
 * @return The reduced row echelon form of the matrix
 */
template <typename GaloisField>
auto rref(const xt::xarray<GaloisFieldElement<GaloisField>> &matrix_input) {
  using ElementType = GaloisFieldElement<GaloisField>;
  auto R = row_echelon(matrix_input);  // Start with row echelon form

  if (R.dimension() != 2) {
    throw std::invalid_argument("rref requires a 2-D array (matrix)");
  }
  if (R.shape(0) == 0 || R.shape(1) == 0) {
    return R;  // Return empty or original if no rows/cols
  }

  auto field =
      R(0, 0).Field();  // Assuming R is not empty and has at least one element
  std::size_t rows = R.shape(0);
  std::size_t cols = R.shape(1);

  // Backward elimination to get reduced row echelon form
  for (int r = static_cast<int>(rows) - 1; r >= 0; --r) {
    // Find the pivot (leading 1) in this row
    std::size_t pivot_col = cols;  // Initialize to an invalid column
    for (std::size_t j = 0; j < cols; ++j) {
      if (R(r, j).Value() == field->MultiplicativeIdentity()) {
        pivot_col = j;
        break;
      }
      // If we encounter a non-zero element before a 1, it's not in RREF from
      // row_echelon This shouldn't happen if row_echelon is correct, but good
      // for robustness
      if (R(r, j).Value() != field->AdditiveIdentity()) {
        break;
      }
    }

    if (pivot_col < cols) {  // If a pivot was found in this row
      // Eliminate elements above the pivot
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

/**
 * Compute the trace of a matrix (sum of diagonal elements)
 *
 * @param a Input square matrix
 * @return The trace
 */
template <typename GaloisField>
auto trace(const xt::xarray<GaloisFieldElement<GaloisField>> &a) {
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

/**
 * Compute the determinant of a square matrix using LU decomposition
 *
 * @param a Input square matrix
 * @return The determinant
 */
template <typename GaloisField>
auto det(const xt::xarray<GaloisFieldElement<GaloisField>> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("determinant requires a square matrix");
  }

  std::size_t n = a.shape(0);
  auto field = a(0, 0).Field();

  // Create a copy for in-place operations
  auto matrix = a;
  ElementType det_result(field->MultiplicativeIdentity(), field);

  // Gaussian elimination with partial pivoting
  for (std::size_t i = 0; i < n; ++i) {
    // Find pivot
    std::size_t pivot_row = i;
    for (std::size_t k = i + 1; k < n; ++k) {
      if (matrix(k, i).Value() != field->AdditiveIdentity()) {
        pivot_row = k;
        break;
      }
    }

    // Swap rows if needed
    if (pivot_row != i) {
      for (std::size_t j = 0; j < n; ++j) {
        std::swap(matrix(i, j), matrix(pivot_row, j));
      }
      det_result =
          det_result *
          ElementType(field->Neg(field->MultiplicativeIdentity()), field);
    }

    // Check for zero pivot
    if (matrix(i, i).Value() == field->AdditiveIdentity()) {
      return ElementType(field->AdditiveIdentity(), field);
    }

    det_result = det_result * matrix(i, i);

    // Eliminate column
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

/**
 * Compute the rank of a matrix over a finite field.
 *
 * @param a Input matrix
 * @return The rank of the matrix
 */
template <typename GaloisField>
std::size_t matrix_rank(const xt::xarray<GaloisFieldElement<GaloisField>> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2) {
    throw std::invalid_argument("matrix_rank requires a 2-D array (matrix)");
  }

  if (a.shape(0) == 0 || a.shape(1) == 0) {
    return 0;
  }

  auto field = a(0, 0).Field();  // Assuming a is not empty
  xt::xarray<ElementType> R = row_echelon(a);

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

//------------------------------------------------------------------------------
// Matrix Solving and Inversion
//------------------------------------------------------------------------------

/**
 * Solve the linear system Ax = b using Gaussian elimination
 *
 * @param a Coefficient matrix A
 * @param b Right-hand side vector b
 * @return Solution vector x
 */
template <typename GaloisField>
auto solve(const xt::xarray<GaloisFieldElement<GaloisField>> &a,
           const xt::xarray<GaloisFieldElement<GaloisField>> &b) {
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

  // Create augmented matrix [A|b]
  xt::xarray<ElementType> augmented({n, n + 1});

  // Initialize to zero first
  ElementType zero(field->AdditiveIdentity(), field);
  std::fill(augmented.begin(), augmented.end(), zero);

  // Copy A and b into augmented matrix
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      augmented(i, j) = a(i, j);
    }
    augmented(i, n) = b(i);
  }

  // Forward elimination with partial pivoting
  for (std::size_t i = 0; i < n; ++i) {
    // Find pivot
    std::size_t pivot_row = i;
    for (std::size_t k = i + 1; k < n; ++k) {
      if (augmented(k, i).Value() != field->AdditiveIdentity()) {
        pivot_row = k;
        break;
      }
    }

    // Swap rows if needed
    if (pivot_row != i) {
      for (std::size_t j = 0; j <= n; ++j) {
        std::swap(augmented(i, j), augmented(pivot_row, j));
      }
    }

    // Check for zero pivot
    if (augmented(i, i).Value() == field->AdditiveIdentity()) {
      throw std::runtime_error("Matrix is singular");
    }

    // Eliminate column
    for (std::size_t k = i + 1; k < n; ++k) {
      if (augmented(k, i).Value() != field->AdditiveIdentity()) {
        auto factor = augmented(k, i) / augmented(i, i);
        for (std::size_t j = i; j <= n; ++j) {
          augmented(k, j) = augmented(k, j) - factor * augmented(i, j);
        }
      }
    }
  }

  // Back substitution
  xt::xarray<ElementType> x({n});

  // Initialize to zero
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

/**
 * Compute the inverse of a square matrix using Gauss-Jordan elimination
 *
 * @param a Input square matrix
 * @return The matrix inverse
 */
template <typename GaloisField>
auto inv(const xt::xarray<GaloisFieldElement<GaloisField>> &a) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 2 || a.shape(0) != a.shape(1)) {
    throw std::invalid_argument("inv requires a square matrix");
  }

  std::size_t n = a.shape(0);
  auto field = a(0, 0).Field();

  // Create augmented matrix [A|I]
  xt::xarray<ElementType> augmented({n, 2 * n});

  // Initialize to zero first
  ElementType zero(field->AdditiveIdentity(), field);
  std::fill(augmented.begin(), augmented.end(), zero);

  // Copy A into left half and identity into right half
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      augmented(i, j) = a(i, j);
      augmented(i, j + n) =
          (i == j) ? ElementType(field->MultiplicativeIdentity(), field)
                   : ElementType(field->AdditiveIdentity(), field);
    }
  }

  // Gauss-Jordan elimination
  for (std::size_t i = 0; i < n; ++i) {
    // Find pivot
    std::size_t pivot_row = i;
    for (std::size_t k = i + 1; k < n; ++k) {
      if (augmented(k, i).Value() != field->AdditiveIdentity()) {
        pivot_row = k;
        break;
      }
    }

    // Swap rows if needed
    if (pivot_row != i) {
      for (std::size_t j = 0; j < 2 * n; ++j) {
        std::swap(augmented(i, j), augmented(pivot_row, j));
      }
    }

    // Check for zero pivot
    if (augmented(i, i).Value() == field->AdditiveIdentity()) {
      throw std::runtime_error("Matrix is singular");
    }

    // Scale pivot row to make diagonal element 1
    auto pivot = augmented(i, i);
    for (std::size_t j = 0; j < 2 * n; ++j) {
      augmented(i, j) = augmented(i, j) / pivot;
    }

    // Eliminate column
    for (std::size_t k = 0; k < n; ++k) {
      if (k != i && augmented(k, i).Value() != field->AdditiveIdentity()) {
        auto factor = augmented(k, i);
        for (std::size_t j = 0; j < 2 * n; ++j) {
          augmented(k, j) = augmented(k, j) - factor * augmented(i, j);
        }
      }
    }
  }

  // Extract inverse from right half
  xt::xarray<ElementType> result({n, n});

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      result(i, j) = augmented(i, j + n);
    }
  }

  return result;
}

/**
 * Compute the cross product of two 3-element vectors.
 * If vectors have 2 elements, the third component is assumed to be zero.
 *
 * @param a First input vector (must have 2 or 3 elements)
 * @param b Second input vector (must have 2 or 3 elements)
 * @return The cross product vector (always 3 elements)
 */
template <typename GaloisField>
auto cross(const xt::xarray<GaloisFieldElement<GaloisField>> &a,
           const xt::xarray<GaloisFieldElement<GaloisField>> &b) {
  using ElementType = GaloisFieldElement<GaloisField>;

  if (a.dimension() != 1 || b.dimension() != 1) {
    throw std::invalid_argument("cross product requires 1-D arrays (vectors)");
  }
  if ((a.shape(0) != 2 && a.shape(0) != 3) ||
      (b.shape(0) != 2 && b.shape(0) != 3)) {
    throw std::invalid_argument(
        "cross product vectors must have 2 or 3 elements");
  }

  auto field = a(0).Field();  // Assuming a is not empty
  ElementType zero(field->AdditiveIdentity(), field);

  ElementType a0 = a(0);
  ElementType a1 = a(1);
  ElementType a2 = (a.shape(0) == 3) ? a(2) : zero;

  ElementType b0 = b(0);
  ElementType b1 = b(1);
  ElementType b2 = (b.shape(0) == 3) ? b(2) : zero;

  xt::xarray<ElementType> result({3});
  result(0) = a1 * b2 - a2 * b1;
  result(1) = a2 * b0 - a0 * b2;
  result(2) = a0 * b1 - a1 * b0;

  return result;
}

//------------------------------------------------------------------------------
// Utility Functions
//------------------------------------------------------------------------------

/**
 * Create an identity matrix over a finite field
 *
 * @param n Size of the identity matrix
 * @param field Shared pointer to the finite field
 * @return Identity matrix of size n x n
 */
template <typename GaloisField>
auto eye(std::size_t n, std::shared_ptr<GaloisField> field) {
  using ElementType = GaloisFieldElement<GaloisField>;

  xt::xarray<ElementType> result({n, n});

  // Initialize to zero first
  ElementType zero(field->AdditiveIdentity(), field);
  std::fill(result.begin(), result.end(), zero);

  // Set diagonal to identity
  for (std::size_t i = 0; i < n; ++i) {
    result(i, i) = ElementType(field->MultiplicativeIdentity(), field);
  }

  return result;
}

/**
 * Create a zero matrix over a finite field
 *
 * @param shape Shape of the matrix
 * @param field Shared pointer to the finite field
 * @return Zero matrix with the specified shape
 */
template <typename GaloisField>
auto zeros(const std::vector<std::size_t> &shape,
           std::shared_ptr<GaloisField> field) {
  using ElementType = GaloisFieldElement<GaloisField>;

  xt::xarray<ElementType> result(shape);

  // Initialize all elements to additive identity
  auto zero = ElementType(field->AdditiveIdentity(), field);
  std::fill(result.begin(), result.end(), zero);

  return result;
}

}  // namespace linalg
}  // namespace xg

#endif  // XGALOIS_LINALG_LINALG_HPP
