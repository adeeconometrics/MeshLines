#ifndef __MATFUNC_H__
#define __MATFUNC_H__

/**
 * @file matfunc.h
 * @author ddamiana
 * @brief Contains linalg functions for processing matrices.
 * @version 0.1
 * @date 2023-05-27
 *
 * @copyright Copyright (c) 2023
 *
 */

// todo

#include "../include/matrix.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <numeric> // inner_product
#include <tuple>
#include <vector>

using lin::Matrix;

using std::tuple;
using std::vector;

/**
 * @brief General transpose function that returns a new matrix.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A Base matrix to be transposed
 * @return Matrix<T, Cols, Rows>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto transpose(const Matrix<T, Rows, Cols> &A) -> Matrix<T, Cols, Rows> {
  Matrix<T, Cols, Rows> TA{};

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < Cols; j++) {
      TA(j, i) = A(i, j);
    }
  }

  return TA;
}
/**
 * @brief Inplace transpose function. Only works with square matrices.
 *
 * @tparam T type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A matrix to be transposed
 */
template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<Rows == Cols>>
auto T(Matrix<T, Rows, Cols> &A) -> void {
  for (std::size_t i = 0; i < Rows; i++)
    for (std::size_t j = i + 1; j < Cols; j++)
      std::swap(A(i, j), A(j, i));
}
/**
 * @brief Returns a copy of upper triangular matrix
 *
 * @tparam T
 * @param A
 * @return Matrix<T>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto triu(const Matrix<T, Rows, Cols> &A) -> Matrix<T, Rows, Cols> {

  Matrix<T, Rows, Cols> res = A;

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < i; j++) {
      res(i, j) = 0;
    }
  }

  return res;
}
/**
 * @brief Masks the matrix's upper triangular part.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto mask_triu(Matrix<T, Rows, Cols> &A) -> void {
  for (std::size_t i = 1; i < Rows; ++i) {
    for (std::size_t j = 0; j < std::min(i, Cols); ++j) {
      A(i, j) = 0;
    }
  }
}
/**
 * @brief Returns the lower triangular part of the matrix
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A The matrix to be modified
 * @return Matrix<T>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto tril(const Matrix<T, Rows, Cols> &A) -> Matrix<T, Rows, Cols> {

  Matrix<T, Rows, Cols> res = A;

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = i + 1; j < Cols; j++) {
      res(i, j) = 0;
    }
  }

  return res;
}
/**
 * @brief Modifies the matrix to mask the lower triangular part
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A The matrix to be modified
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto mask_tril(Matrix<T, Rows, Cols> &A) -> void {

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = i + 1; j < Cols; j++) {
      A(i, j) = 0;
    }
  }
}
/**
 * @brief Returns a diagonal matrix.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A Base matrix
 * @return Matrix<T, Rows, Cols>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto diag(const Matrix<T, Rows, Cols> &A) -> Matrix<T, Rows, Cols> {
  Matrix<T, Rows, Cols> D{};

  for (std::size_t i = 0; i < Rows; i++) {
    D(i, i) = A(i, i);
  }

  return D;
}
/**
 * @brief Modifies the matrix to mask the diagonal part
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A Base matrix
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto mask_diag(Matrix<T, Rows, Cols> &A) -> void {
  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < Cols; j++) {
      if (i != j) {
        A(i, j) = 0;
      }
    }
  }
}
/**
 * @brief Returns the trace of the Matrix;
 * this function requires a square matrix
 *
 * @tparam T The type of the matrix
 * @param M The matrix to be processed
 * @return T The trace of the matrix
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto trace(const Matrix<T, Rows, Cols> &M) -> T {
  static_assert(Rows == Cols, "Matrix must be square");
  T sum{};

  for (std::size_t i = 0; i < Rows; i++) {
    sum += M(i, i);
  }
  return sum;
}
/**
 * @brief Returns a Minor submatrix of M.
 * Note that Matrixindexing is still zero-based.
 *
 * @tparam T
 * @param M
 * @param row
 * @param col
 * @return Matrix<T>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
minor_submatrix(const Matrix<T, Rows, Cols> &M, std::size_t row = 0,
                std::size_t col = 0) -> Matrix<T, Rows - 1, Cols - 1> {

  Matrix<T, Rows - 1, Cols - 1> minor{};

  std::size_t t_row = 0;
  for (std::size_t i = 0; i < Rows; ++i) {
    if (i == row)
      continue;

    std::size_t t_col = 0;
    for (std::size_t j = 0; j < Cols; ++j) {
      if (j == col)
        continue;

      minor(t_row, t_col) = M(i, j);
      t_col += 1;
    }

    t_row += 1;
  }

  return minor;
}
/**
 * @brief Returns the minor determinant of matrix M.
 *
 * @tparam T
 * @param M
 * @param row
 * @param col
 * @return double
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto minor(const Matrix<T, Rows, Cols> &M, std::size_t row = 0,
                     std::size_t col = 0) -> double {

  return det(minor_submatrix(M, row, col));
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
minor_matrix(const Matrix<T, Rows, Cols> &M) -> Matrix<double, Rows, Cols> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<double, Rows, Cols> minor_mat(Rows, vector<T>(Cols));

  for (std::size_t i = 0; i < Rows; ++i) {
    for (std::size_t j = 0; j < Cols; ++j) {
      const auto minor_det = minor(M, i, j);
      minor_mat(i, j) = minor_det;
    }
  }

  return minor_mat;
}
/**
 * @brief Returns a cofactor submatrix of $M$.
 *
 * @tparam T
 * @param M
 * @param row
 * @param col
 * @return Matrix<T>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
cofactor_submatrix(const Matrix<T, Rows, Cols> &M, std::size_t row = 0,
                   std::size_t col = 0) -> Matrix<T, Rows, Cols> {

  Matrix<T, Rows - 1, Cols - 1> cofactor{};

  std::size_t t_row = 0;
  for (std::size_t i = 0; i < Rows; ++i) {
    if (i == row)
      continue;

    std::size_t t_col = 0;
    for (std::size_t j = 0; j < Cols; ++j) {
      if (j == col)
        continue;

      cofactor(t_row, t_col) = ((i + j) % 2 == 0) ? M(i, j) : -M(i, j);
      t_col += 1;
    }

    t_row += 1;
  }

  return cofactor;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
cofactor_matrix(const Matrix<T, Rows, Cols> &M) -> Matrix<double, Rows, Cols> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<double, Rows, Cols> cofactor{};

  for (std::size_t i = 0; i < Rows; ++i) {
    for (std::size_t j = 0; j < Cols; ++j) {
      const auto minor_det = minor(M, i, j);
      cofactor(i, j) = ((i + j) % 2 == 0) ? minor_det : -minor_det;
    }
  }

  return cofactor;
}
/**
 * @brief Returns the adjugate matrix
 *
 * @tparam T
 * @param M
 * @return Matrix<double>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
adj(const Matrix<T, Rows, Cols> &M) -> Matrix<double, Rows, Cols> {
  return transpose(cofactor_matrix(M));
}
#endif // __MATFUNC_H__