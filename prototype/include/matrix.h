#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "../include/factorize.h"
#include "../include/vector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace lin {

using std::is_arithmetic_v, std::common_type_t;
using std::runtime_error;
using std::vector;

template <typename T> using Matrix = vector<vector<T>>;

template <typename T, typename U = T>
constexpr auto operator+(const Matrix<T> &lhs, const Matrix<U> &rhs)
    -> Matrix<common_type_t<T, U>> {
  Matrix<common_type_t<T, U>> result(lhs.size(), std::vector<T>(lhs[0].size()));

  // for (size_t i{}; i < lhs.size(); i++)
  //   std::transform(lhs[i].cbegin(), lhs[i].cend(), rhs[i].cbegin(),
  //                  result[i].begin(),
  //                  [](const T &a, const U &b) { return a + b; });

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const vector<T> &a, const vector<U> &b) { return a + b; });

  return result;
}

template <typename T, typename U = T>
constexpr auto operator-(const Matrix<T> &lhs, const Matrix<U> &rhs)
    -> Matrix<common_type_t<T, U>> {
  Matrix<common_type_t<T, U>> result(lhs.size(), std::vector<T>(lhs[0].size()));

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const vector<T> &a, const vector<U> &b) { return a - b; });

  return result;
}

template <typename T, typename U = T>
constexpr auto operator*(const Matrix<T> &lhs, const Matrix<U> &rhs)
    -> Matrix<common_type_t<T, U>> {
  Matrix<common_type_t<T, U>> result(lhs.size(), std::vector<T>(lhs[0].size()));

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const vector<T> &a, const vector<U> &b) { return a * b; });

  return result;
}

template <typename T, typename U = T>
constexpr auto operator/(const Matrix<T> &lhs, const Matrix<U> &rhs)
    -> Matrix<common_type_t<T, U>> {
  Matrix<common_type_t<T, U>> result(lhs.size(), std::vector<T>(lhs[0].size()));

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const vector<T> &a, const vector<U> &b) { return a / b; });

  return result;
}

} // namespace lin

#endif // __MATRIX_H__