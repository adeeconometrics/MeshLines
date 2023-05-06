#include "../include/factorize.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

using std::cout;
using std::endl;

// crout's algorithm
template <typename T>
auto lu_decomposition(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
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
template <typename T>
auto qr_factorization(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  Matrix<T> Q(A.size(), vector<T>(A.size(), 0));
  Matrix<T> R(A.size(), vector<T>(A.size(), 0));

  // Copy A to R
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A.size(); ++j) {
      R[i][j] = A[i][j];
    }
  }

  // Compute Q and R using the Gram-Schmidt process
  for (int j = 0; j < A.size(); ++j) {
    // Compute the jth column of Q
    for (int i = 0; i < A.size(); ++i) {
      Q[i][j] = R[i][j];
    }
    for (int k = 0; k < j; ++k) {
      T dot_product = 0;
      for (int i = 0; i < A.size(); ++i) {
        dot_product += Q[i][k] * R[i][j];
      }
      for (int i = 0; i < A.size(); ++i) {
        Q[i][j] -= dot_product * Q[i][k];
      }
    }
    // Normalize the jth column of Q
    T norm = 0;
    for (int i = 0; i < A.size(); ++i) {
      norm += Q[i][j] * Q[i][j];
    }
    norm = std::sqrt(norm);
    for (int i = 0; i < A.size(); ++i) {
      Q[i][j] /= norm;
    }
    // Compute the jth row of R
    for (int i = j; i < A.size(); ++i) {
      R[j][i] = 0;
      for (int k = 0; k < A.size(); ++k) {
        R[j][i] += Q[k][j] * A[k][i];
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

auto test_lu() -> void {
  const Matrix<float> A = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  auto [L, U] = lu_decomposition(A);
}

// auto main() -> int {
//   // Test the LU decomposition function
//   Matrix<float> A = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};

//   auto [L, U] = lu_decomposition(A);

//   cout << "L matrix:" << endl;
//   for (const auto &row : L) {
//     for (const auto &element : row) {
//       cout << element << " ";
//     }
//     cout << endl;
//   }

//   cout << "U matrix:" << endl;
//   for (const auto &row : U) {
//     for (const auto &element : row) {
//       cout << element << " ";
//     }
//     cout << endl;
//   }

//   return 0;
// }
