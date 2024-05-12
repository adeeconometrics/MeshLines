#ifndef __MATMUL_H__
#define __MATMUL_H__

#include "../include/matmul.hpp"
#include "../include/matrix.hpp"
#include "../include/utils.hpp"

#include <algorithm>
#include <future>
#include <memory>
#include <thread>
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

template <typename T, std::size_t M, std::size_t N,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
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

template <typename T, std::size_t N, std::size_t M,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto gemm(const Matrix<T, N, M> &t_lhs,
          const Matrix<T, N, M> &t_rhs) -> Matrix<T, N, M> {
  Matrix<T, N, M> result;
  constexpr std::size_t block_size = (16 * 1024) / sizeof(T);

  // Loop over the blocks
  for (std::size_t i = 0; i < N; i += block_size) {
    for (std::size_t j = 0; j < N; j += block_size) {
      for (std::size_t k = 0; k < N; k += block_size) {
        // Multiply the blocks
        for (std::size_t ii = i; ii < std::min(i + block_size, N); ++ii) {
          for (std::size_t jj = j; jj < std::min(j + block_size, N); ++jj) {
            T sum{};
            for (std::size_t kk = k; kk < std::min(k + block_size, N); ++kk) {
              sum += t_lhs(ii, kk) * t_rhs(kk, jj);
            }
            result(ii, jj) += sum;
          }
        }
      }
    }
  }

  return result;
}

template <typename T, std::size_t N, std::size_t M,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto threaded_gemm(const Matrix<T, N, M> &t_lhs,
                   const Matrix<T, N, M> &t_rhs) -> Matrix<T, N, M> {
  Matrix<T, N, M> result;
  constexpr std::size_t block_size = (32 * 1024) / sizeof(T);
  const std::size_t num_threads = std::thread::hardware_concurrency();
  std::vector<std::thread> threads;

  for (std::size_t i = 0; i < num_threads; ++i) {
    threads.emplace_back([i, num_threads, block_size, &t_lhs, &t_rhs,
                          &result]() {
      for (std::size_t ii = i; ii < N; ii += num_threads) {
        for (std::size_t j = 0; j < N; j += block_size) {
          for (std::size_t k = 0; k < N; k += block_size) {
            for (std::size_t jj = j; jj < std::min(j + block_size, N); ++jj) {
              T sum{};
              for (std::size_t kk = k; kk < std::min(k + block_size, N); ++kk) {
                sum += t_lhs(ii, kk) * t_rhs(kk, jj);
              }
              result(ii, jj) += sum;
            }
          }
        }
      }
    });
  }
  for (auto &thread : threads) {
    thread.join();
  }
  return result;
}

template <typename T, std::size_t N, std::size_t M,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto async_gemm(const Matrix<T, N, M> &t_lhs,
                const Matrix<T, N, M> &t_rhs) -> Matrix<T, N, M> {
  Matrix<T, N, M> result;
  std::array<std::future<void>, N> futures;

  for (std::size_t i = 0; i < N; ++i) {
    futures[i] = std::async(std::launch::async, [i, &t_lhs, &t_rhs, &result]() {
      for (std::size_t j = 0; j < N; ++j) {
        T sum{};
        for (std::size_t k = 0; k < N; ++k) {
          sum += t_lhs(i, k) * t_rhs(k, j);
        }
        result(i, j) = sum;
      }
    });
  }

  return result;
}

#endif // __MATMUL_H__