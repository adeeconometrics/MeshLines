#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "../include/factorize.h"
#include "../include/vector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <cassert>

namespace lin {

using std::is_arithmetic_v, std::common_type_t;
using std::runtime_error;
using std::vector;

template <typename T> using Matrix = vector<vector<T>>;

template <typename T, typename U = T>
constexpr auto operator+(const Matrix<T> &lhs, const Matrix<U> &rhs)
    -> Matrix<common_type_t<T, U>> {
  Matrix<common_type_t<T, U>> result(lhs.size(), std::vector<T>(lhs[0].size()));

  assert(lhs.size() != rhs.size());

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

template <typename T, typename U = T>
constexpr auto operator+(const Matrix<T>& lhs, const vector<U>& rhs)
	-> Matrix<common_type_t<T,U>> {
  
  assert(lhs.size() == rhs.size());

  using result_type = common_type_t<T,U>;

  Matrix<result_type> result(lhs.size(), std::vector<result_type>(rhs.size())); // is this the right way to add vector and mat?
  
  return result;

}

template <typename T, typename U = T>
constexpr auto operator+(const Matrix<T>& lhs, U rhs)
	->Matrix<common_type_t<T,U>> {

  using result_type = common_type_t<T,U>;
  
  Matrix<result_type> result {};
  std::for_each(lhs.cbegin(), lhs.cend(),
	[&result, rhs](const auto &_lhs) {result.emplace_back(_lhs + rhs);}
  );

  return result; 
}

} // namespace lin

#endif // __MATRIX_H__
