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
 * @brief transpose function that returns a new Matrix
 *
 * @tparam T
 * @param A
 * @return Matrix<T>
 */
template <typename T> auto transpose(const Matrix<T> &A) -> Matrix<T> {
  Matrix<T> TA = A;

  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  for (std::size_t i = 0; i < rows; i++)
    for (std::size_t j = 0; j < cols; j++)
      TA[j][i] = A[i][j];

  return TA;
}
/**
 * @brief In-place transpose
 *
 * @tparam U arithmetic type
 * @param A Matrix
 */
template <typename U> auto T(Matrix<U> &A) -> void {
  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  for (std::size_t i = 0; i < rows; i++)
    for (std::size_t j = i + 1; j < cols; j++)
      std::swap(A[i][j], A[j][i]);
}
/**
 * @brief Returns a copy of upper triangular matrix
 *
 * @tparam T
 * @param A
 * @return Matrix<T>
 */
template <typename T> auto triu(const Matrix<T> &A) -> Matrix<T> {
  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  Matrix<T> res{rows, vector<T>(cols)};

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i; j < cols; j++) {
      if (j < rows && i < cols) {
        res[i][j] = A[i][j];
      }
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
template <typename T> auto mask_triu(Matrix<T> &A) -> void {
  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i + 1; j < cols; j++) {
      if (j < rows && i < cols) {
        A[j][i] = 0;
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
template <typename T> auto tril(const Matrix<T> &A) -> Matrix<T> {
  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  Matrix<T> res{rows, vector<T>(cols)};

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i; j < cols; j++) {
      if (j < rows && i < cols) {
        res[j][i] = A[j][i];
      }
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
template <typename T> auto mask_tril(Matrix<T> &A) -> void {
  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i + 1; j < cols; j++) {
      A[i][j] = 0;
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
template <typename T> auto trace(const Matrix<T> &M) -> T {
  assert(M.size() == M[0].size());
  T sum{};

  const std::size_t length = M.size();

  for (std::size_t i = 0; i < length; i++) {
    sum += M[i][i];
  }
  return sum;
}
/**
 * @brief Crout's algorithm implementation of LU Decomposition
 *
 * @tparam T
 * @param A
 * @return tuple<Matrix<T>, Matrix<T>>
 */
template <typename T>
auto lu_crout(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {

  const std::size_t _size = A.size();
  Matrix<T> U(A);
  Matrix<T> L(_size, vector<T>(_size, 0));

  // Initialize L to the identity matrix
  for (std::size_t i = 0; i < L.size(); ++i) {
    L[i][i] = 1;
  }

  // Perform LU decomposition
  for (std::size_t j = 0; j < _size; ++j) {
    for (std::size_t i = j + 1; i < _size; ++i) {
      T factor = U[i][j] / U[j][j];
      L[i][j] = factor;
      for (std::size_t k = j; k < _size; ++k) {
        U[i][k] -= factor * U[j][k];
      }
    }
  }

  return std::make_tuple(L, U);
}
/**
 * @brief LU Decomposition using Gaussian Elimination -- ideal for sparse
 * matrix and more numerically stable than Crout's algorithm.
 * This implementation is also more efficient when the matrix is large.
 *
 * @tparam T
 * @param A
 * @return tuple<Matrix<T>, Matrix<T>>
 */
template <typename T>
auto lu_gaussian(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  const std::size_t _size = A.size();
  Matrix<T> U(_size, vector<T>(_size, 0));
  Matrix<T> L(_size, vector<T>(_size, 0));

  for (std::size_t i = 0; i < _size; i++) {
    for (std::size_t upper_iter = i; upper_iter < _size; upper_iter++) {
      T sum = 0;
      for (std::size_t inner_iter = 0; inner_iter < i; inner_iter++) {
        sum += L[i][inner_iter] * U[inner_iter][upper_iter];
      }
      U[i][upper_iter] = A[i][upper_iter] - sum;
    }

    for (std::size_t lower_iter = 0; lower_iter < _size; lower_iter++) {
      if (i == lower_iter) {
        L[i][i] = 1;
      } else {
        T sum = 0;
        for (std::size_t inner_iter = 0; inner_iter < i; inner_iter++) {
          sum += L[lower_iter][inner_iter] * U[inner_iter][i];
        }
        L[lower_iter][i] = (A[lower_iter][i] - sum) / U[i][i];
      }
    }
  }

  return std::make_tuple(L, U);
}

template <typename T>
auto plu(Matrix<T> A) -> tuple<Matrix<T>, Matrix<T>, Matrix<T>> {

  const std::size_t n = A.size();
  Matrix<T> L = lin::zero_mat<T>(n);
  Matrix<T> U = lin::zero_mat<T>(n);
  Matrix<T> P = lin::id<T>(n);

  for (std::size_t k = 0; k < n; k++) {
    // Find pivot row and swap rows
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++) {
      if (std::abs(A[i][k]) > std::abs(A[p][k])) {
        p = i;
      }
    }
    if (p != k) {
      std::swap(A[k], A[p]);
      std::swap(P[k], P[p]);
    }

    // Perform elimination
    for (std::size_t i = k + 1; i < n; i++) {
      T factor = A[i][k] / A[k][k];
      A[i][k] = factor;
      for (std::size_t j = k + 1; j < n; j++) {
        A[i][j] -= A[k][j] * factor;
      }
    }
  }

  // Extract L and U from A
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++) {
      if (i > j) {
        L[i][j] = A[i][j];
        U[i][j] = 0;
      } else if (i == j) {
        L[i][j] = 1;
        U[i][j] = A[i][j];
      } else {
        L[i][j] = 0;
        U[i][j] = A[i][j];
      }
    }
  }

  return std::make_tuple(P, L, U);
}
// gram-schmidt process
// conditions: linearly independent cols
template <typename T>
auto qr_gm(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  Matrix<T> Q(A.size(), vector<T>(A.size(), 0));
  Matrix<T> R = A;

  // Compute Q and R using the Gram-Schmidt process
  for (std::size_t j = 0; j < A.size(); ++j) {
    // Compute the jth column of Q
    for (std::size_t i = 0; i < A.size(); ++i) {
      Q[i][j] = R[i][j];
    }
    for (std::size_t k = 0; k < j; ++k) {
      T dot_product = 0;
      for (std::size_t i = 0; i < A.size(); ++i) {
        dot_product += Q[i][k] * R[i][j];
      }
      for (std::size_t i = 0; i < A.size(); ++i) {
        Q[i][j] -= dot_product * Q[i][k];
      }
    }
    // Normalize the jth column of Q
    T norm = 0;
    for (std::size_t i = 0; i < A.size(); ++i) {
      norm += Q[i][j] * Q[i][j];
    }
    norm = std::sqrt(norm);
    for (std::size_t i = 0; i < A.size(); ++i) {
      Q[i][j] /= norm;
    }
    // Compute the jth row of R
    for (std::size_t i = j; i < A.size(); ++i) {
      R[j][i] = 0;
      for (std::size_t k = 0; k < A.size(); ++k) {
        R[j][i] += Q[k][j] * A[k][i];
      }
    }
  }
  mask_triu(R);

  return std::make_tuple(Q, R);
}

template <typename T>
auto qr_householder(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  const std::size_t m = A.size();
  const std::size_t n = A[0].size();

  // Initialize Q and R
  Matrix<T> Q(m, std::vector<T>(m));
  Matrix<T> R(A);

  // Compute Householder reflections and apply them to R
  for (std::size_t k = 0; k < n; k++) {
    std::vector<T> x(m - k);
    for (std::size_t i = k; i < m; i++) {
      x[i - k] = R[i][k];
    }

    T norm_x =
        std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));

    std::vector<T> v(m - k);
    v[0] = x[0] < 0 ? x[0] - norm_x : x[0] + norm_x;
    for (std::size_t i = 1; i < m - k; i++) {
      v[i] = x[i];
    }

    T norm_v =
        std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));

    Matrix<T> H(m - k, std::vector<T>(m - k));
    for (std::size_t i = 0; i < m - k; i++) {
      for (std::size_t j = 0; j < m - k; j++) {
        if (i == j) {
          H[i][j] = 1 - 2 * v[i] * v[j] / (norm_v * norm_v);
        } else {
          H[i][j] = -2 * v[i] * v[j] / (norm_v * norm_v);
        }
      }
    }

    // Update R with H
    for (std::size_t i = k; i < m; i++) {
      for (std::size_t j = k; j < n; j++) {
        T sum = 0;
        for (std::size_t h = 0; h < m - k; h++) {
          sum += H[i - k][h] * R[h + k][j];
        }
        R[i][j] = sum;
      }
    }

    // Update Q with H
    for (std::size_t i = 0; i < m; i++) {
      for (std::size_t j = k; j < m; j++) {
        T sum = 0;
        for (std::size_t h = 0; h < m - k; h++) {
          sum += Q[i][h + k] * H[h][j - k];
        }
        Q[i][j] = sum;
      }
    }
  }

  return std::make_tuple(Q, R);
}

// ldl factorization
// conditions: square, hermitian positive definite matrix
// type COMPLEX
template <typename T>
auto ldl(const Matrix<T> &A) -> tuple<Matrix<T>, vector<T>> {
  const std::size_t size = A.size();

  Matrix<T> L{size, vector<T>(size, 0)};
  vector<T> D{size, 0};

  for (std::size_t j = 0; j < size; ++j) {
    T sum = 0;
    for (std::size_t k = 0; k < j; ++k) {
      sum += L[j][k] * L[j][k] * D[k];
    }
    D[j] = A[j][j] - sum;
    L[j][j] = 1;

    for (std::size_t i = j + 1; i < size; ++i) {
      sum = 0;
      for (std::size_t k = 0; k < j; ++k) {
        sum += L[i][k] * L[j][k] * D[k];
      }
      L[i][j] = (A[i][j] - sum) / D[j];
    }
  }

  return std::make_tuple(L, D);
}

// cholesky decomposition (L*L^T): returns L
// conditions: symmetric positive definite matrix
template <typename T> constexpr auto cholesky(const Matrix<T> &A) -> Matrix<T> {
  const std::size_t size = A.size();
  Matrix<T> L(size, vector<T>(size, 0));

  // perform cholesky decomposition
  for (std::size_t i = 0; i < size; i++) {
    for (std::size_t j = 0; j <= i; j++) {
      T sum = 0;
      for (std::size_t k = 0; k < j; k++) {
        sum += L[i][k] * L[j][k];
      }
      if (i == j) {
        L[i][i] = std::sqrt(A[i][i] - sum);
      } else {
        L[i][j] = (A[i][j] - sum) / L[j][j];
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
template <typename T>
constexpr auto det(
    const Matrix<T> &M,
    std::function<tuple<Matrix<T>, Matrix<T>>(const Matrix<T> &)> LUMethod =
        [](const Matrix<T> &M) { return lu_gaussian(M); }) -> double {
  assert(M.size() == M[0].size());

  double p0{1}, p1{1};

  auto LU = LUMethod(M);

  for (std::size_t i = 0; i < M.size(); i++) {
    p0 *= std::get<0>(LU)[i][i];
    p1 *= std::get<1>(LU)[i][i];
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
template <typename T> constexpr auto det_2(const Matrix<T> &M) -> double {
  return M[0][0] * M[1][1] - M[0][1] * M[1][0];
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
template <typename T> constexpr auto rref(const Matrix<T> &M) -> Matrix<T> {
  Matrix<T> result = M;

  std::size_t lead = 0;
  const std::size_t rows = M.size();

  for (std::size_t row = 0; row < rows; ++row) {
    if (lead >= M[row].size())
      break;

    auto it = std::find_if(
        std::begin(result) + row, std::end(result), [lead](const auto &row) {
          return std::abs(row[lead]) >= std::numeric_limits<T>::epsilon();
        });

    if (it == result.end()) {
      ++lead;
      continue;
    }

    std::swap(result[row], *it);
    T lv = result[row][lead];

    std::transform(std::cbegin(result[row]), std::cend(result[row]),
                   std::begin(result[row]),
                   [lv](const auto &val) { return val / lv; });

    for (std::size_t i = 0; i < rows; ++i) {
      if (i != row) {
        T lv = result[i][lead];
        std::transform(
            result[row].begin(), result[row].end(), result[i].begin(),
            result[i].begin(),
            [lv](const auto &a, const auto &b) { return b - lv * a; });
      }
    }

    ++lead;
  }

  return result;
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
template <typename T>
constexpr auto minor_submatrix(const Matrix<T> &M, std::size_t row = 0,
                               std::size_t col = 0) -> Matrix<T> {

  const std::size_t m_row = M.size();
  const std::size_t m_col = M[0].size();

  Matrix<T> minor(m_row - 1, vector<T>(m_col - 1));

  std::size_t t_row = 0;
  for (std::size_t i = 0; i < m_row; ++i) {
    if (i == row)
      continue;

    std::size_t t_col = 0;
    for (std::size_t j = 0; j < m_col; ++j) {
      if (j == col)
        continue;

      minor[t_row][t_col] = M[i][j];
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
template <typename T>
constexpr auto minor(const Matrix<T> &M, std::size_t row = 0,
                     std::size_t col = 0) -> double {

  return det(minor_submatrix(M, row, col));
}

template <typename T>
constexpr auto minor_matrix(const Matrix<T> &M) -> Matrix<double> {
  const std::size_t m_row = M.size();
  const std::size_t m_col = M[0].size();

  assert(m_row == m_col);

  Matrix<double> minor_mat(m_row, vector<T>(m_col));

  for (std::size_t i = 0; i < m_row; ++i) {
    for (std::size_t j = 0; j < m_col; ++j) {
      const auto minor_det = minor(M, i, j);
      minor_mat[i][j] = minor_det;
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
template <typename T>
constexpr auto cofactor_submatrix(const Matrix<T> &M, std::size_t row = 0,
                                  std::size_t col = 0) -> Matrix<T> {
  const std::size_t m_row = M.size();
  const std::size_t m_col = M[0].size();

  Matrix<T> cofactor(m_row - 1, vector<T>(m_col - 1));

  std::size_t t_row = 0;
  for (std::size_t i = 0; i < m_row; ++i) {
    if (i == row)
      continue;

    std::size_t t_col = 0;
    for (std::size_t j = 0; j < m_col; ++j) {
      if (j == col)
        continue;

      cofactor[t_row][t_col] = ((i + j) % 2 == 0) ? M[i][j] : -M[i][j];
      t_col += 1;
    }

    t_row += 1;
  }

  return cofactor;
}

template <typename T>
constexpr auto cofactor_matrix(const Matrix<T> &M) -> Matrix<double> {
  const std::size_t m_row = M.size();
  const std::size_t m_col = M[0].size();

  assert(m_row == m_col);

  Matrix<double> cofactor(m_row, vector<double>(m_col));

  for (std::size_t i = 0; i < m_row; ++i) {
    for (std::size_t j = 0; j < m_col; ++j) {
      const auto minor_det = minor(M, i, j);
      cofactor[i][j] = ((i + j) % 2 == 0) ? minor_det : -minor_det;
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
template <typename T> constexpr auto adj(const Matrix<T> &M) -> Matrix<double> {
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
template <typename T> constexpr auto rank(const Matrix<T> &M) -> int {}

// span
// basis

// template <typename T>
// auto svd(const Matrix<T>& A) -> tuple<Matrix<T>, Matrix<T>, Matrix<T>> {

// };

#endif // __MATFUNC_H__