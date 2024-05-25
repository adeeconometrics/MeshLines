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
// make functions exclusive to integral type

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
    for (std::size_t j = i + 1; j < Cols; j++) {
      if (j < Rows && i < Cols)
        res(j, i) = 0;
    }
  }

  return res;
}
/**
 * @brief Modifies the matrix to mask upper triangular part
 *
 * @tparam T
 * @param A
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto mask_triu(Matrix<T, Rows, Cols> &A) -> void {

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = i + 1; j < Cols; j++) {
      if (j < Rows && i < Cols) {
        A(j, i) = 0;
      }
    }
  }
}
/**
 * @brief Returns the lower triangular part of the matrix
 *
 * @tparam T
 * @param A
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
 * @tparam T
 * @param A
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
 * @brief Returns the trace of the Matrix;
 * this function requires a square matrix
 *
 * @tparam T
 * @param M
 * @return T
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
 * @brief Crout's algorithm implementation of LU Decomposition; see if this can
 * work with rectangular matrices.
 *
 * @tparam T
 * @param A
 * @return tuple<Matrix<T>, Matrix<T>>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto lu_crout(const Matrix<T, Rows, Cols> &A)
    -> tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<T, Rows, Cols> U(A);
  Matrix<T, Rows, Cols> L = lin::id<T, Rows>();

  // Perform LU decomposition
  for (std::size_t j = 0; j < Rows; ++j) {
    for (std::size_t i = j + 1; i < Rows; ++i) {
      T factor = U(i, j) / U(j, j);
      L(i, j) = factor;
      for (std::size_t k = j; k < Rows; ++k) {
        U(i, k) -= factor * U(j, k);
      }
    }
  }

  return std::make_tuple(L, U);
}
/**
 * @brief LU Decomposition using Gaussian Elimination -- ideal for sparse
 * matrix and more numerically stable than Crout's algorithm.
 * This implementation is also more efficient when the matrix is large. See if
 * it can work with rectangular matrices.
 *
 * @tparam T
 * @param A
 * @return tuple<Matrix<T>, Matrix<T>>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
auto lu_gaussian(const Matrix<T, Rows, Cols> &A)
    -> tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<T, Rows, Cols> U;
  Matrix<T, Rows, Cols> L;

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t upper_iter = i; upper_iter < Rows; upper_iter++) {
      T sum = 0;
      for (std::size_t inner_iter = 0; inner_iter < i; inner_iter++) {
        sum += L(i, inner_iter) * U(inner_iter, upper_iter);
      }
      U(i, upper_iter) = A(i, upper_iter) - sum;
    }

    for (std::size_t lower_iter = 0; lower_iter < Rows; lower_iter++) {
      if (i == lower_iter) {
        L(i, i) = 1;
      } else {
        T sum = 0;
        for (std::size_t inner_iter = 0; inner_iter < i; inner_iter++) {
          sum += L(lower_iter, inner_iter) * U(inner_iter, i);
        }
        L(lower_iter, i) = (A(lower_iter, i) - sum) / U(i, i);
      }
    }
  }

  return std::make_tuple(L, U);
}

template <typename T, std::size_t Rows, std::size_t Cols>
auto plu(Matrix<T, Rows, Cols> A)
    -> tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>,
             Matrix<T, Rows, Cols>> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<T, Rows, Cols> L = lin::zero_mat<T, Rows, Cols>();
  Matrix<T, Rows, Cols> U = lin::zero_mat<T, Rows, Cols>();
  Matrix<T, Rows, Cols> P = lin::id<T, Rows>();

  for (std::size_t k = 0; k < Rows; k++) {
    // Find pivot row and swap rows
    std::size_t p = k;
    for (std::size_t i = k + 1; i < Rows; i++) {
      if (std::abs(A(i, k)) > std::abs(A(p, k))) {
        p = i;
      }
    }
    if (p != k) {
      std::swap(A(k), A(p));
      std::swap(P(k), P(p));
    }

    // Perform elimination
    for (std::size_t i = k + 1; i < Rows; i++) {
      T factor = A(i, k) / A(k, k);
      A(i, k) = factor;
      for (std::size_t j = k + 1; j < Rows; j++) {
        A(i, j) -= A(k, j) * factor;
      }
    }
  }

  // Extract L and U from A
  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < Rows; j++) {
      if (i > j) {
        L(i, j) = A(i, j);
        U(i, j) = 0;
      } else if (i == j) {
        L(i, j) = 1;
        U(i, j) = A(i, j);
      } else {
        L(i, j) = 0;
        U(i, j) = A(i, j);
      }
    }
  }

  return std::make_tuple(P, L, U);
}
// gram-schmidt process
// conditions: linearly independent cols
template <typename T, std::size_t Rows, std::size_t Cols>
auto qr_gm(const Matrix<T, Rows, Cols> &A)
    -> tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<T, Rows, Cols> Q{};
  Matrix<T, Rows, Cols> R = A;

  // Compute Q and R using the Gram-Schmidt process
  for (std::size_t j = 0; j < Rows; ++j) {
    // Compute the jth column of Q
    for (std::size_t i = 0; i < Rows; ++i) {
      Q(i, j) = R(i, j);
    }
    for (std::size_t k = 0; k < j; ++k) {
      T dot_product = 0;
      for (std::size_t i = 0; i < Rows; ++i) {
        dot_product += Q(i, k) * R(i, j);
      }
      for (std::size_t i = 0; i < Rows; ++i) {
        Q(i, j) -= dot_product * Q(i, k);
      }
    }
    // Normalize the jth column of Q
    T norm = 0;
    for (std::size_t i = 0; i < Rows; ++i) {
      norm += Q(i, j) * Q(i, j);
    }
    norm = std::sqrt(norm);
    for (std::size_t i = 0; i < Rows; ++i) {
      Q(i, j) /= norm;
    }
    // Compute the jth row of R
    for (std::size_t i = j; i < Rows; ++i) {
      R(j, i) = 0;
      for (std::size_t k = 0; k < Rows; ++k) {
        R(j, i) += Q(k, j) * A(k, i);
      }
    }
  }
  mask_triu(R);

  return std::make_tuple(Q, R);
}

template <typename T, std::size_t Rows, std::size_t Cols>
auto qr_householder(const Matrix<T, Rows, Cols> &A)
    -> tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>> {

  static_assert(Rows == Cols, "Matrix is assumed to be square");

  // Initialize Q and R
  Matrix<T, Rows, Cols> Q = lin::zero_mat<T, Rows, Cols>();
  Matrix<T, Rows, Cols> R(A);

  // Compute Householder reflections and apply them to R
  for (std::size_t k = 0; k < Cols; k++) {
    std::vector<T> x(Rows - k);
    for (std::size_t i = k; i < Rows; i++) {
      x[i - k] = R(i, k);
    }

    T norm_x =
        std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));

    std::vector<T> v(Rows - k);
    v[0] = x[0] < 0 ? x[0] - norm_x : x[0] + norm_x;
    for (std::size_t i = 1; i < Rows - k; i++) {
      v[i] = x[i];
    }

    T norm_v =
        std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));

    Matrix<T, Rows - k, Rows - k> H{};
    for (std::size_t i = 0; i < Rows - k; i++) {
      for (std::size_t j = 0; j < Rows - k; j++) {
        if (i == j) {
          H(i, j) = 1 - 2 * v[i] * v[j] / (norm_v * norm_v);
        } else {
          H(i, j) = -2 * v[i] * v[j] / (norm_v * norm_v);
        }
      }
    }

    // Update R with H
    for (std::size_t i = k; i < Rows; i++) {
      for (std::size_t j = k; j < Cols; j++) {
        T sum = 0;
        for (std::size_t h = 0; h < Rows - k; h++) {
          sum += H(i - k, h) * R(h + k, j);
        }
        R(i, j) = sum;
      }
    }

    // Update Q with H
    for (std::size_t i = 0; i < Rows; i++) {
      for (std::size_t j = k; j < Rows; j++) {
        T sum = 0;
        for (std::size_t h = 0; h < Rows - k; h++) {
          sum += Q(i, h + k) * H(h, j - k);
        }
        Q(i, j) = sum;
      }
    }
  }

  return std::make_tuple(Q, R);
}

// ldl factorization
// conditions: square, hermitian positive definite matrix. Check for rectangular
// matrix implementation. type COMPLEX
template <typename T, std::size_t Rows, std::size_t Cols>
auto ldl(const Matrix<T, Rows, Cols> &A)
    -> tuple<Matrix<T, Rows, Cols>, vector<T>> {

  static_assert(Rows == Cols, "Matrix assumed to be square");

  Matrix<T, Rows, Cols> L = lin::zero_mat<T, Rows, Cols>();
  vector<T> D{Rows, 0};

  for (std::size_t j = 0; j < Rows; ++j) {
    T sum = 0;
    for (std::size_t k = 0; k < j; ++k) {
      sum += L(j, k) * L(j, k) * D[k];
    }
    D[j] = A(j, j) - sum;
    L(j, j) = 1;

    for (std::size_t i = j + 1; i < Rows; ++i) {
      sum = 0;
      for (std::size_t k = 0; k < j; ++k) {
        sum += L(i, k) * L(j, k) * D[k];
      }
      L(i, j) = (A(i, j) - sum) / D[j];
    }
  }

  return std::make_tuple(L, D);
}

// cholesky decomposition (L*L^T): returns L
// conditions: symmetric positive definite matrix
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
cholesky(const Matrix<T, Rows, Cols> &A) -> Matrix<T, Rows, Cols> {

  static_assert(Rows == Cols, "Matrix is assumed to be square");
  Matrix<T, Rows, Cols> L = lin::zero_mat<T, Rows, Cols>();

  // perform cholesky decomposition
  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j <= i; j++) {
      T sum = 0;
      for (std::size_t k = 0; k < j; k++) {
        sum += L(i, k) * L(j, k);
      }
      if (i == j) {
        L(i, i) = std::sqrt(A(i, i) - sum);
      } else {
        L(i, j) = (A(i, j) - sum) / L(j, j);
      }
    }
  }

  return L;
}
/**
 * @brief This function returns the determinant of a matrix.
 * This also keeps the LUMethod open to the client.
 *
 * @tparam T
 * @param M
 * @param LUMethod
 * @return double
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto det(
    const Matrix<T, Rows, Cols> &M,
    std::function<tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>>(
        const Matrix<T, Rows, Cols> &)>
        LUMethod = [](const Matrix<T, Rows, Cols> &M) {
          return lu_gaussian(M);
        }) -> double {

  static_assert(Rows == Cols, "Matrix must be square");

  double p0{1}, p1{1};

  auto LU = LUMethod(M);

  for (std::size_t i = 0; i < Rows; i++) {
    p0 *= std::get<0>(LU)(i, i);
    p1 *= std::get<1>(LU)(i, i);
  }

  return p1 * p0;
}
/**
 * @brief Specialized implementation of determinant for 2x2 matrix
 *
 * @tparam T
 * @param M
 * @return double
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto det_2(const Matrix<T, Rows, Cols> &M) -> double {
  return M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
}
/**
 * @brief Returns a reduced row echelon form of matrix M.
 * Note that not all matrices have RREF. The implementation
 * of this function does not currently handle checking for those
 * exceptional cases.
 *
 * @tparam T
 * @param M
 * @return Matrix<T>
 */
// template <typename T, std::size_t Rows, std::size_t Cols>
// constexpr auto rref(const Matrix<T, Rows, Cols> &M) -> Matrix<T, Rows, Cols>
// {
//   Matrix<T, Rows, Cols> result = M;

//   std::size_t lead = 0;

//   for (std::size_t row = 0; row < Rows; ++row) {
//     if (lead >= M[row].size())
//       break;

//     auto it = std::find_if(
//         std::begin(result) + row, std::end(result), [lead](const auto &row) {
//           return std::abs(row[lead]) >= std::numeric_limits<T>::epsilon();
//         });

//     if (it == result.end()) {
//       ++lead;
//       continue;
//     }

//     std::swap(result[row], *it);
//     T lv = result[row][lead];

//     std::transform(std::cbegin(result[row]), std::cend(result[row]),
//                    std::begin(result[row]),
//                    [lv](const auto &val) { return val / lv; });

//     for (std::size_t i = 0; i < Rows; ++i) {
//       if (i != row) {
//         T lv = result[i][lead];
//         std::transform(
//             result[row].begin(), result[row].end(), result[i].begin(),
//             result[i].begin(),
//             [lv](const auto &a, const auto &b) { return b - lv * a; });
//       }
//     }

//     ++lead;
//   }

//   return result;
// }
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

// template <typename T> constexpr auto inv(const Matrix<T> &M) ->
// Matrix<double> {
//   return 1 / (det(A)) * transpose(cofactor_matrix(M));
// }

// rowspace
// nullspace
// colspace

// rank
// template <typename T, std::size_t Rows, std::size_t Cols>
// constexpr auto rank(const Matrix<T, Rows, Cols> &M) -> int {}

// span
// basis

// template <typename T>
// auto svd(const Matrix<T>& A) -> tuple<Matrix<T>, Matrix<T>, Matrix<T>> {

// };

#endif // __MATFUNC_H__