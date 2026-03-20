/**
 * @file linalg_example.cpp
 * @brief Examples for linear algebra operations over Galois fields using direct
 * constructors
 */

#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"

using namespace xg;

// Helper function to print a field element
template <typename ElementType>
std::string ElementToString(const ElementType& element) {
  std::stringstream ss;
  element.Print(ss);
  return ss.str();
}

int main() {
  std::cout << "XGalois Linear Algebra Examples\n";
  std::cout << "===============================\n";
  std::cout << "This example uses direct field constructors (NO factory)\n";

  // Create fields using direct constructors
  std::cout << "\n=== Creating Fields for Linear Algebra Examples ===\n";

  auto gf7_field = std::make_shared<GFP<uint8_t>>(7);
  std::cout << "Created GF(7) using direct GFP constructor\n";

  auto gf2_field = std::make_shared<GF2>();
  std::cout << "Created GF(2) using direct GF2 constructor\n";

  auto gf11_field = std::make_shared<GFP<uint8_t>>(11);
  std::cout << "Created GF(11) using direct GFP constructor\n";

  // Linear algebra over GF(7)
  std::cout << "\n=== Linear Algebra over GF(7) ===\n";

  // Create field element type for convenience
  using Element7 = GaloisFieldElementBase<GFP<uint8_t>>;

  // Create a 3x3 matrix over GF(7)
  std::cout << "Creating a 3x3 matrix over GF(7):\n";
  std::vector<std::vector<Element7>> matrix_data = {
      {Element7(1, gf7_field), Element7(2, gf7_field), Element7(3, gf7_field)},
      {Element7(4, gf7_field), Element7(5, gf7_field), Element7(6, gf7_field)},
      {Element7(0, gf7_field), Element7(1, gf7_field), Element7(2, gf7_field)}};

  // Display the matrix
  std::cout << "Matrix A:\n";
  for (size_t i = 0; i < matrix_data.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < matrix_data[i].size(); ++j) {
      std::cout << ElementToString(matrix_data[i][j]);
      if (j < matrix_data[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Create another 3x3 matrix
  std::cout << "\nCreating another 3x3 matrix over GF(7):\n";
  std::vector<std::vector<Element7>> matrix_data2 = {
      {Element7(2, gf7_field), Element7(0, gf7_field), Element7(1, gf7_field)},
      {Element7(1, gf7_field), Element7(3, gf7_field), Element7(2, gf7_field)},
      {Element7(4, gf7_field), Element7(1, gf7_field), Element7(5, gf7_field)}};

  // Display the second matrix
  std::cout << "Matrix B:\n";
  for (size_t i = 0; i < matrix_data2.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < matrix_data2[i].size(); ++j) {
      std::cout << ElementToString(matrix_data2[i][j]);
      if (j < matrix_data2[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Matrix addition
  std::cout << "\n=== Matrix Addition ===\n";
  std::vector<std::vector<Element7>> sum_matrix(
      3, std::vector<Element7>(3, Element7(0, gf7_field)));

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      sum_matrix[i][j] = matrix_data[i][j] + matrix_data2[i][j];
    }
  }

  std::cout << "A + B:\n";
  for (size_t i = 0; i < sum_matrix.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < sum_matrix[i].size(); ++j) {
      std::cout << ElementToString(sum_matrix[i][j]);
      if (j < sum_matrix[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Matrix multiplication
  std::cout << "\n=== Matrix Multiplication ===\n";
  std::vector<std::vector<Element7>> product_matrix(
      3, std::vector<Element7>(3, Element7(0, gf7_field)));

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      Element7 sum(0, gf7_field);
      for (size_t k = 0; k < 3; ++k) {
        sum = sum + (matrix_data[i][k] * matrix_data2[k][j]);
      }
      product_matrix[i][j] = sum;
    }
  }

  std::cout << "A * B:\n";
  for (size_t i = 0; i < product_matrix.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < product_matrix[i].size(); ++j) {
      std::cout << ElementToString(product_matrix[i][j]);
      if (j < product_matrix[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Vector operations
  std::cout << "\n=== Vector Operations over GF(7) ===\n";

  // Create vectors
  std::vector<Element7> vector1 = {
      Element7(1, gf7_field), Element7(3, gf7_field), Element7(5, gf7_field)};
  std::vector<Element7> vector2 = {
      Element7(2, gf7_field), Element7(4, gf7_field), Element7(6, gf7_field)};

  std::cout << "Vector v1: [ ";
  for (size_t i = 0; i < vector1.size(); ++i) {
    std::cout << ElementToString(vector1[i]);
    if (i < vector1.size() - 1) std::cout << ", ";
  }
  std::cout << " ]\n";

  std::cout << "Vector v2: [ ";
  for (size_t i = 0; i < vector2.size(); ++i) {
    std::cout << ElementToString(vector2[i]);
    if (i < vector2.size() - 1) std::cout << ", ";
  }
  std::cout << " ]\n";

  // Vector addition
  std::vector<Element7> vector_sum(3, Element7(0, gf7_field));
  for (size_t i = 0; i < 3; ++i) {
    vector_sum[i] = vector1[i] + vector2[i];
  }

  std::cout << "v1 + v2: [ ";
  for (size_t i = 0; i < vector_sum.size(); ++i) {
    std::cout << ElementToString(vector_sum[i]);
    if (i < vector_sum.size() - 1) std::cout << ", ";
  }
  std::cout << " ]\n";

  // Dot product
  Element7 dot_product(0, gf7_field);
  for (size_t i = 0; i < 3; ++i) {
    dot_product = dot_product + (vector1[i] * vector2[i]);
  }
  std::cout << "v1 · v2 = " << ElementToString(dot_product) << "\n";

  // Scalar multiplication
  Element7 scalar(3, gf7_field);
  std::vector<Element7> scaled_vector(3, Element7(0, gf7_field));
  for (size_t i = 0; i < 3; ++i) {
    scaled_vector[i] = scalar * vector1[i];
  }

  std::cout << "3 * v1: [ ";
  for (size_t i = 0; i < scaled_vector.size(); ++i) {
    std::cout << ElementToString(scaled_vector[i]);
    if (i < scaled_vector.size() - 1) std::cout << ", ";
  }
  std::cout << " ]\n";

  // Linear algebra over GF(2)
  std::cout << "\n=== Linear Algebra over GF(2) ===\n";

  using Element2 = GaloisFieldElementBase<GF2>;

  // Create a 2x2 matrix over GF(2)
  std::cout << "Creating matrices over GF(2):\n";
  std::vector<std::vector<Element2>> bin_matrix1 = {
      {Element2(1, gf2_field), Element2(0, gf2_field)},
      {Element2(1, gf2_field), Element2(1, gf2_field)}};

  std::vector<std::vector<Element2>> bin_matrix2 = {
      {Element2(0, gf2_field), Element2(1, gf2_field)},
      {Element2(1, gf2_field), Element2(0, gf2_field)}};

  std::cout << "Binary matrix C:\n";
  for (size_t i = 0; i < bin_matrix1.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < bin_matrix1[i].size(); ++j) {
      std::cout << ElementToString(bin_matrix1[i][j]);
      if (j < bin_matrix1[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  std::cout << "Binary matrix D:\n";
  for (size_t i = 0; i < bin_matrix2.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < bin_matrix2[i].size(); ++j) {
      std::cout << ElementToString(bin_matrix2[i][j]);
      if (j < bin_matrix2[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Binary matrix multiplication
  std::vector<std::vector<Element2>> bin_product(
      2, std::vector<Element2>(2, Element2(0, gf2_field)));

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      Element2 sum(0, gf2_field);
      for (size_t k = 0; k < 2; ++k) {
        sum = sum + (bin_matrix1[i][k] * bin_matrix2[k][j]);
      }
      bin_product[i][j] = sum;
    }
  }

  std::cout << "C * D:\n";
  for (size_t i = 0; i < bin_product.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < bin_product[i].size(); ++j) {
      std::cout << ElementToString(bin_product[i][j]);
      if (j < bin_product[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Identity matrix
  std::cout << "\n=== Identity Matrix ===\n";
  std::vector<std::vector<Element7>> identity(
      3, std::vector<Element7>(3, Element7(0, gf7_field)));
  for (size_t i = 0; i < 3; ++i) {
    identity[i][i] = Element7(1, gf7_field);
  }

  std::cout << "3x3 Identity matrix over GF(7):\n";
  for (size_t i = 0; i < identity.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < identity[i].size(); ++j) {
      std::cout << ElementToString(identity[i][j]);
      if (j < identity[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  // Verify identity property: A * I = A
  std::cout << "\nVerifying A * I = A:\n";
  std::vector<std::vector<Element7>> result(
      3, std::vector<Element7>(3, Element7(0, gf7_field)));

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      Element7 sum(0, gf7_field);
      for (size_t k = 0; k < 3; ++k) {
        sum = sum + (matrix_data[i][k] * identity[k][j]);
      }
      result[i][j] = sum;
    }
  }

  std::cout << "A * I:\n";
  for (size_t i = 0; i < result.size(); ++i) {
    std::cout << "  [ ";
    for (size_t j = 0; j < result[i].size(); ++j) {
      std::cout << ElementToString(result[i][j]);
      if (j < result[i].size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
  }

  std::cout << "\nLinear algebra examples completed successfully!\n";
  std::cout
      << "Note: This example uses direct field constructors (NO factory).\n";
  std::cout
      << "All operations are performed in the respective finite fields.\n";

  return 0;
}
