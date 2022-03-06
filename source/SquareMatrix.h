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
#include <utility>

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