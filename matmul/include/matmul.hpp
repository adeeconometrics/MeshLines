#ifndef __MATMUL_H__
#define __MATMUL_H__

#include "../include/matmul.hpp"
#include "../include/matrix.hpp"
#include "../include/utils.hpp"

#include <algorithm>
#include <type_traits>

/**
 * @brief This section contain different implementation of algorithms to doing
 * matmul.
 *
 */

template <typename T, std::size_t M, std::size_t N,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto iterative(const Matrix<T, M, N> &t_lhs,
               const Matrix<T, M, N> &t_rhs) -> Matrix<T, M, N> {
  Matrix<T, M, N> result;
  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      T sum = 0;
      for (std::size_t k = 0; k < N; ++k) {
        sum += t_lhs(i, k) * t_rhs(k, j);
      }
      result(i, j) = sum;
    }
  }
  return result;
}

template <typename T, std::size_t M, std::size_t N>
auto loop_reorder(const Matrix<T, M, N> &t_lhs,
                  const Matrix<T, M, N> &t_rhs) -> Matrix<T, M, N> {
  Matrix<T, M, N> result;
  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      for (std::size_t k = 0; k < N; ++k) {
        result(i, k) += t_lhs(i, j) * t_rhs(j, k);
      }
    }
  }
  return result;
}

#endif // __MATMUL_H__