#ifndef __MATFUNC_H__
#define __MATFUNC_H__

// todo
// make functions exclusive to integral type

#include "../include/matrix.h"

#include <algorithm>
#include <cmath>
#include <numeric> // inner_product
#include <tuple>
#include <vector>

using lin::Matrix;

using std::tuple;
using std::vector;

template <typename T> auto zero_mat(size_t n) -> Matrix<T> {
  return Matrix<T>{n, vector<T>(n, 0)};
}

template <typename T> auto eye_mat(size_t n) -> Matrix<T> {
  Matrix<T> res = zero_mat<T>(n);

  for (std::size_t i{}; i < n; i++) {
    res[i][i] = 1;
  }
  return res;
}
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
  Matrix<T> L = zero_mat<T>(n);
  Matrix<T> U = zero_mat<T>(n);
  Matrix<T> P = eye_mat<T>(n);

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
  Matrix<T> R(A.size(), vector<T>(A.size(), 0));

  // Copy A to R
  for (std::size_t i = 0; i < A.size(); ++i) {
    for (std::size_t j = 0; j < A.size(); ++j) {
      R[i][j] = A[i][j];
    }
  }

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
template <typename T> auto cholesky(const Matrix<T> &A) -> Matrix<T> {
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

// template <typename T>
// auto svd(const Matrix<T>& A) -> tuple<Matrix<T>, Matrix<T>, Matrix<T>> {

// };

#endif // __MATFUNC_H__