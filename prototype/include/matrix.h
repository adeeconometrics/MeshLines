#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "../include/factorize.h"
#include "../include/vector.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>

namespace lin {

using std::common_type_t;
using std::transform, std::for_each;
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
// see if it is possible to add U parameter for implicit typecast
template <typename T>
auto operator+=(Matrix<T> &lhs, const Matrix<T> &rhs) -> Matrix<T> & {
  assert(lhs.size() == rhs.size());

  transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
            [](const auto &_lhs, const auto &_rhs) { return _lhs += _rhs; });
  return lhs;
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
constexpr auto operator+(const Matrix<T> &lhs, const vector<U> &rhs)
    -> Matrix<common_type_t<T, U>> {

  assert(lhs[0].size() == rhs.size());
  assert(lhs.size() > 0 && rhs.size() > 0);

  using result_type = common_type_t<T, U>;

  Matrix<result_type> result{};

  transform(lhs.cbegin(), lhs.cend(), result.begin(),
            [&](const vector<T> &row) -> vector<result_type> {
              assert(row.size() == rhs.size());
              vector<result_type> row_result(row.size());
              transform(row.cbegin(), row.cend(), row_result.begin(),
                        std::plus<result_type>());
              return row_result;
            });

  return result;
}

template <typename T>
constexpr auto operator+(Matrix<T> lhs, int rhs)
    -> Matrix<common_type_t<T, int>> {

  using result_type = common_type_t<T, int>;

  Matrix<result_type> result{};
  for (const auto vec : lhs) {
    vector<result_type> __tmp{};
    for (const auto elem : vec) {
      __tmp.emplace_back(static_cast<result_type>(elem + rhs));
    }
    result.emplace_back(__tmp);
  }

  return result;
}

} // namespace lin

#endif // __MATRIX_H__
