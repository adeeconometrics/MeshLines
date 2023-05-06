#ifndef __FACTORIZE_H__
#define __FACTORIZE_H__

#include <cmath>
#include <iostream>
#include <numeric> // inner_product
#include <tuple>
#include <vector>

using std::cout;
using std::endl;
using std::tuple;
using std::vector;

template <typename T> using Matrix = vector<vector<T>>;

// crout's algorithm
template <typename T>
auto lu_crout(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  Matrix<T> L(A.size(), vector<T>(A.size(), 0));
  Matrix<T> U(A.size(), vector<T>(A.size(), 0));

  // Copy A to U
  for (std::size_t i = 0; i < A.size(); ++i) {
    for (std::size_t j = 0; j < A.size(); ++j) {
      U[i][j] = A[i][j];
    }
  }

  // Initialize L to the identity matrix
  for (std::size_t i = 0; i < L.size(); ++i) {
    L[i][i] = 1;
  }

  // Perform LU decomposition
  for (std::size_t j = 0; j < U.size(); ++j) {
    for (std::size_t i = j + 1; i < U.size(); ++i) {
      T factor = U[i][j] / U[j][j];
      L[i][j] = factor;
      for (std::size_t k = j; k < U.size(); ++k) {
        U[i][k] -= factor * U[j][k];
      }
    }
  }

  return std::make_tuple(L, U);
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

#endif // __FACTORIZE_H__