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
// test against different sizes of SquareMatrix

#pragma once
#include "Vector.h"
#include <algorithm>
#include <type_traits>
#include <utility>
#include <math.h>

template <typename T, size_t N> using SquareMatrix = vector<vector<T, N>, N>;

template <typename T, size_t N>
auto operator+(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T, N> {

  SquareMatrix<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
                 [](const auto &a, const auto &b) { return a + b; });
  return result;
}

template <typename T, size_t N>
auto operator-(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T, N> {

  SquareMatrix<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
                 [](const auto &a, const auto &b) { return a - b; });
  return result;
}

template <typename T, size_t N>
auto operator*(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T, N> {

  SquareMatrix<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
                 [](const auto &a, const auto &b) { return a * b; });
  return result;
}

template <typename T, size_t N>
auto operator/(const SquareMatrix<T, N> &lhs, const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T, N> {

  SquareMatrix<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.begin(), result.begin(),
                 [](const auto &a, const auto &b) { return a / b; });
  return result;
}

template <typename T, size_t N>
constexpr auto operator==(const SquareMatrix<T, N> &lhs,
                          const SquareMatrix<T, N> &rhs) -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}

template <typename T, size_t N>
constexpr auto operator!=(const SquareMatrix<T, N> &lhs,
                          const SquareMatrix<T, N> &rhs) -> bool {
  return !(lhs == rhs);
}

template <typename T, size_t N>
auto operator*(const SquareMatrix<T, N> &lhs, const vector<T, N> &rhs)
    -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const auto &a) { return dot(a, rhs); });

  return result;
}

template <typename T, size_t N, typename U = double,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
constexpr auto operator+(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {

  SquareMatrix<U, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const auto &a) { return a + rhs; });

  return result;
}

template <typename T, size_t N, typename U = double,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator-(const SquareMatrix<T, N> &lhs, U rhs) -> SquareMatrix<U, N> {

  SquareMatrix<U, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const auto &a) { return a - rhs; });

  return result;
}

template <typename T, size_t N, typename U = double,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
constexpr auto operator*(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {

  SquareMatrix<U, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const auto &a) { return a * rhs; });

  return result;
}

template <typename T, size_t N, typename U = double,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
constexpr auto operator/(const SquareMatrix<T, N> &lhs, U rhs) noexcept
    -> SquareMatrix<U, N> {

  SquareMatrix<U, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), result.begin(),
                 [&rhs](const auto &a) { return a / rhs; });

  return result;
}

template <typename T, size_t N>
constexpr auto trace(const SquareMatrix<T, N> &rhs) -> double {

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
constexpr auto transpose_matrix(SquareMatrix<T, N> &rhs) -> SquareMatrix<T, N> {

  SquareMatrix<T, N> transpose{};
  for (size_t row{}; row < rhs.size(); row++)
    for (size_t col{}; col < rhs.size(); col++)
      std::swap(rhs[row][col], transpose[col][row]);
  return transpose;
}

template <typename T, size_t N>
constexpr auto row_vector(SquareMatrix<T, N> &rhs, int index) -> vector<T, N> {
  return rhs[index];
}

template <typename T, size_t N>
constexpr auto col_vector(SquareMatrix<T, N> &rhs, int index) -> vector<T, N> {
  vector<T, N> col_vec{};
  std::transform(rhs.cbegin(), rhs.cend(), col_vec.begin(),
                 [&index](const auto &a) { return a[index]; });

  return col_vec;
}

template <typename T, size_t N>
constexpr auto lu_decomposition(const SquareMatrix<T,N>& M) 
  -> std::pair<SquareMatrix<T,N>, SquareMatrix<T,N>>{

    SquareMatrix<T,N> lower{}, upper{};
    
    auto upper_triangular = [&lower, &upper, &M](size_t i, size_t k)mutable {
        T sum{};
        for(size_t j{}; j < i; ++j) 
            sum += (lower[i][j] * upper[j][k]);
            
        upper[i][k] = M[i][k] - sum;
    };
    
    auto lower_triangular = [&lower, &upper, &M](size_t i, size_t k)mutable {
        if (i == k) 
                lower[i][i] = 1;
        else {
            T sum{};
            for(size_t j{}; j < i; ++j) 
                sum += (lower[k][j] * upper[j][i]);
                
            lower[k][i] = (M[k][i] - sum) /upper[i][i];
        }
    };
    
    for (size_t i{}; i < N; ++i){
        for(size_t k{i}; k < N; ++k){
            upper_triangular(i,k);
            lower_triangular(i,k);
        }
    }
    
    return std::make_pair(lower, upper);
}

// template <typename T, size_t N>
// constexpr auto qr_decomposition(const SquareMatrix<T,N>& M) 
//   -> std::pair<SquareMatrix<T,N>, SquareMatrix<T,N>>{

// }

template <typename T, size_t N>
constexpr auto cholesky_decomposition(const SquareMatrix<T, N> &M)
    -> std::pair<SquareMatrix<T,N>, SquareMatrix<T, N>> {
  // static_assert(is_symmetric(M));
    SquareMatrix<T, N> lower{};

    for (size_t i{}; i < N; ++i)
      for (size_t j{}; j <= i; ++j) {
        T sum{};
        if (i == j) {
          for (size_t k{}; k < i; ++k)
            sum += pow(lower[i][k], 2);
          lower[i][i] = sqrt(M[i][i] - sum);
        } else {
          for (size_t k{}; k < j; ++k)
            sum += lower[i][k] * lower[j][k];
          lower[i][j] = (M[i][j] - sum) / lower[j][j];
        }
      }

    return lower;
  }
}

template <typename T, size_t N>
constexpr auto det(const SquareMatrix<T, N> &rhs) -> double {
  double p0{1}, p1{1};

  auto LU = lu_decomposition(rhs);

  for(size_t i{}; i < N; ++i){
    p0 *= LU.first[i][i];
    p1 *= LU.second[i][i];
  }

  return p0*p1;
}

template <typename T, size_t N>
constexpr auto det_2(const SquareMatrix<T, N> &rhs) -> double {
  return rhs[0][0]*rhs[1][1] - rhs[0][1]*rhs[1][0];
}

template <typename T, size_t N>
constexpr auto minor_submatrix(const SquareMatrix<T, N> &rhs, size_t row = 0,
                               size_t col = 0) -> SquareMatrix<T, N - 1> {

  SquareMatrix<T, N - 1> _minor{};
  size_t idx_row{0}, idx_col{0};
  bool flag{false};

  for (size_t i{}; i < N; ++i) {
    for (size_t j{}; j < N; ++j) {
      if (i == row || j == col) {
        flag = true;
        continue;
      }
      _minor[idx_row][idx_col] = rhs[i][j];
      idx_col += 1;
    }
    if (flag)
      flag = false;
    else
      idx_row += 1;
  }
  return _minor;
}

template <typename T, size_t N>
constexpr auto cofactor_submatrix(const SquareMatrix<T, N> &rhs, size_t row = 0,
                                  size_t col = 0) -> SquareMatrix<T, N - 1> {

  SquareMatrix<T, N - 1> _cofactor{};
  size_t idx_row{0}, idx_col{0};
  bool flag{false};

  for (size_t i{}; i < N; ++i) {
    for (size_t j{}; j < N; ++j) {
      if (i == row || j == col) {
        flag = true;
        continue;
      }
      _cofactor[idx_row][idx_col] = ((i + j) % 2 == 0) ? rhs[i][j] : -rhs[i][j];
      idx_col += 1;
    }
    if (flag)
      flag = false;
    else
      idx_row += 1;
  }
  return _cofactor;
}

template <typename T, size_t N>
constexpr auto cofactor(const SquareMatrix<T, N> &rhs, 
                        size_t row = 0, size_t col = 0) -> double {
  return det(cofactor_submatrix(rhs, row, col));
}

template <typename T, size_t N>
constexpr auto cofactor_matrix(const SquareMatrix<T, N> &rhs) 
    -> SquareMatrix<T, N> {

  SquareMatrix<T, N> _cofactor{};
  for (size_t i{}; i < N; ++i) {
    for (size_t j{}; j < N; ++j) {
      _cofactor[i][j] = cofactor(rhs, i, j);
    }
  }

  return _cofactor;
}

template <typename T, size_t N>
constexpr auto minor(const SquareMatrix<T, N> &rhs, size_t row = 0,
                       size_t col = 0) -> double {
  return det(minor_submatrix(rhs, row, col));
}

template <typename T, size_t N>
constexpr auto minor_matrix(const SquareMatrix<T, N> &rhs)
    -> SquareMatrix<T, N> {

  SquareMatrix<T, N> _minor{};
  for (size_t i{}; i < N; ++i) {
    for (size_t j{}; j < N; ++j) {
      _minor[i][j] = minor(rhs, i, j);
    }
  }

  return _minor;
}

template <typename T, size_t N>
constexpr auto min(const SquareMatrix<T, N> &rhs) noexcept -> T {
  T min_element{rhs[0][0]};
  for (const auto i : rhs) {
    if (min_element > min(i))
      min_element = i;
  }
  return min;
}

template <typename T, size_t N>
constexpr auto max(const SquareMatrix<T, N> &rhs) noexcept -> T {
  T max_element{};
  for (const auto i : rhs) {
    if (max_element < max(i))
      max_element = i;
  }
  return max_element;
}