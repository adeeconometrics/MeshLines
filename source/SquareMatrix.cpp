/**
 * @file SquareMatrix.cpp
 * @author ddamiana
 * @brief Concrete implementation of square matrix
 * @version 0.1
 * @date 2022-03-03
 *
 * @copyright Copyright (c) 2022
 *
 */

// todo: create namespace for Vector
// #pragma once
#include "Vector.h"
#include <algorithm>

template <typename T, size_t N> using SquareMatrix = vector<vector<T, N>, N>;

template <typename T, size_t N>
auto operator+(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T,N> {
  SquareMatrix<T, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
      [](const T &a, const T &b) -> decltype(a + b) { return a + b; });
  return result;
}

template <typename T, size_t N>
auto operator-(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T,N> {
  SquareMatrix<T, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
      [](const T &a, const T &b) -> decltype(a - b) { return a - b; });
  return result;
}

template <typename T, size_t N>
auto operator*(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T,N> {
  SquareMatrix<T, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
      [](const T &a, const T &b) -> decltype(a * b) { return a * b; });
  return result;
}

template <typename T, size_t N>
auto operator/(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T,N> {
  SquareMatrix<T, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
      [](const T &a, const T &b) -> decltype(a / b) { return a / b; });
  return result;
}

// template <typename T, size_t N>
// auto operator*(const SquareMatrix<T, N> &lhs, const vector<T, N> &rhs)
//     -> vector<T,N> {
//   vector<T, N> result{};
//   std::transform(
//       lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
//       [](const T &a, const T &b) -> decltype(a * b) { return a * b; });
//   return result;
// }

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator+(const SquareMatrix<T, N> &lhs, U rhs) noexcept -> SquareMatrix<U, N> {
  SquareMatrix<U, N> result{};

  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const T &a) -> decltype(a + rhs) { return a + rhs; });

  return result;
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator-(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {
  SquareMatrix<U, N> result{};

  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const T &a) -> decltype(a - rhs) { return a - rhs; });

  return result;
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator*(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {
  SquareMatrix<U, N> result{};

  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const T &a) -> decltype(a * rhs) { return a * rhs; });

  return result;
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator/(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {
  SquareMatrix<U, N> result{};

  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const T &a) -> decltype(a / rhs) { return a / rhs; });

  return result;
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator+(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {
  SquareMatrix<U, N> result{};

  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const T &a) -> decltype(a + rhs) { return a + rhs; });

  return result;
}

// template <typename T, size_t N>
// auto det(const SquareMatrix<T, N> &rhs) -> double {}

template <typename T, size_t N>
auto trace(const SquareMatrix<T, N> &rhs) -> double {
  double result{};
  for (size_t i{}; i < rhs.size(); ++i) {
    result += rhs[i][i];
  }
  return result;
}

template <typename T, size_t N>
auto transpose(SquareMatrix<T, N> &rhs) -> void {
  for (size_t row{}; row < rhs.size() - 1; row++) {
    for (size_t col{row + 1}; col < rhs.size(); col++) {
      std::swap(rhs[row][col], rhs[col][row]);
    }
  }
}

template <typename T, size_t N>
auto transpose_matrix(SquareMatrix<T, N> &rhs) -> SquareMatrix<T, N> {
  SquareMatrix<T, N> transpose{};
  for (size_t row{}; row < rhs.size(); row++)
    for (size_t col{}; col < rhs.size(); col++)
      std::swap(rhs[row][col], transpose[col][row]);
  return transpose;
}

// template <typename T, size_t N>
// auto row_matrix(SquareMatrix<T, N> &rhs) -> SquareMatrix<T, 1>{
// }

// template <typename T, size_t N>
// auto col_matrix(SquareMatrix<T, N> &rhs) -> SquareMatrix<T, 1> {}