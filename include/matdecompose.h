#ifndef __MATDECOMPOSE_H__
#define __MATDECOMPOSE_H__

/**
 * @file matdecompose.hpp
 * @author ddamiana
 * @brief Contains functions for matrix factorization
 * @version 0.1
 * @date 2024-05-26
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "../include/matfunc.h"
#include "../include/matops.h"
#include "../include/matrix.h"

#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>

namespace lin {

/**
 * @brief Crout's LU decomposition implementation. This function returns a tuple
 * of two matrices L and U such that matmul(L,U) = A.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A The matrix to be decomposed
 * @return tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>>
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
 * @brief This is a Gaussian elimination implementation of LU Decomposition.
 * WARNING: The current implementation is brittle for integer types as it may
 * result to segfaults due to division by zero (which is why a conditional flag
 * is temporarily set to call this method when T is a floating type). Future
 * efforts will be made to investigate further on how to handle this.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @tparam std::enable_if_t<std::is_floating_point_v<T>>
 * @param A The matrix to be decomposed
 * @return tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>>
 */
template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<std::is_floating_point_v<T>>>
auto lu_gaussian(const Matrix<T, Rows, Cols> &A)
    -> tuple<Matrix<T, Rows, Cols>, Matrix<T, Rows, Cols>> {

  static_assert(Rows == Cols, "Matrix must be square");

  Matrix<T, Rows, Cols> U = lin::scalar_mat<T, Rows, Cols>(T{});
  Matrix<T, Rows, Cols> L = lin::scalar_mat<T, Rows, Cols>(T{});

  for (std::size_t i{}; i < Rows; i++) {
    // upper triangular
    for (std::size_t k{i}; k < Rows; k++) {
      T sum{0.0};
      for (std::size_t j{}; j < Rows; j++) {
        sum += (L(i, j) * U(j, k));
      }
      U(i, k) = A(i, k) - sum;
    }

    // lower triangular
    for (std::size_t k{i}; k < Rows; k++) {
      if (i == k) {
        L(i, i) = 1.0;
      } else {
        T sum{0.0};
        for (std::size_t j{}; j < i; j++) {
          sum += (L(k, j) * U(j, i));
        }
        L(k, i) = (A(k, i) - sum) / U(i, i);
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

/**
 * @brief Cholesky decomposition implementation. This function returns a matrix
 * L such that matmul(L, L^T) = A. Note that this function assumes that the
 * matrix is symmetric positive definite. Future efforts will be added to check
 * this condition before proceeding with the decomposition.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param A The matrix to be decomposed
 * @return Matrix<T, Rows, Cols>
 */
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

} // namespace lin
#endif // __MATDECOMPOSE_H__