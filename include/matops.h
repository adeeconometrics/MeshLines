#ifndef __MATOPS_H__
#define __MATOPS_H__

#include "../include/matrix.h"
#include "../include/vecops.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>

// todo
// - [x] Mat [op] Vec
// - [x] Mat =[op] Vec
// - [ ] Mat [op] Scalar
// - [ ] Mat =[op] Scalar

namespace lin {

using std::common_type_t;
using std::transform;
using std::vector;

template <typename T, typename U = T>
constexpr auto operator+(const Matrix<T> &lhs, const Matrix<U> &rhs)
    -> Matrix<common_type_t<T, U>> {
  Matrix<common_type_t<T, U>> result(lhs.size(), std::vector<T>(lhs[0].size()));

  assert(lhs.size() == rhs.size());

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

template <typename T>
constexpr auto operator+=(Matrix<T> &lhs, const Matrix<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs += _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator-=(Matrix<T> &lhs, const Matrix<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs -= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator*=(Matrix<T> &lhs, const Matrix<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs *= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator/=(Matrix<T> &lhs, const Matrix<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs /= _rhs; });
  return lhs;
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

template <typename T, typename U = T>
constexpr auto operator-(const Matrix<T> &lhs, const vector<U> &rhs)
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
                        std::minus<result_type>());
              return row_result;
            });

  return result;
}

template <typename T, typename U = T>
constexpr auto operator*(const Matrix<T> &lhs, const vector<U> &rhs)
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
                        std::multiplies<result_type>());
              return row_result;
            });

  return result;
}
// MAKE TYPES FUSE WITH FLOATING POINT
template <typename T, typename U = T>
constexpr auto operator/(const Matrix<T> &lhs, const vector<U> &rhs)
    -> Matrix<common_type_t<T, U>> {

  assert(lhs[0].size() == rhs.size());
  assert(lhs.size() > 0 && rhs.size() > 0);

  // using result_type = common_type_t<T, U, double>;
  using result_type = common_type_t<T, U>;

  Matrix<result_type> result{};

  transform(lhs.cbegin(), lhs.cend(), result.begin(),
            [&](const vector<T> &row) -> vector<result_type> {
              assert(row.size() == rhs.size());
              vector<result_type> row_result(row.size());
              transform(row.cbegin(), row.cend(), row_result.begin(),
                        std::divides<result_type>());
              return row_result;
            });

  return result;
}

template <typename T>
constexpr auto operator+=(Matrix<T> &lhs, const vector<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs += _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator-=(Matrix<T> &lhs, const vector<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs -= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator*=(Matrix<T> &lhs, const vector<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs *= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator/=(Matrix<T> &lhs, const vector<T> &rhs) -> Matrix<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs /= _rhs; });
  return lhs;
}

// template <typename T>
// constexpr auto operator+(const Matrix<T> &lhs, int rhs)
//     -> Matrix<common_type_t<T, int>> {

//   // static_assert(std::is_arithmetic_v<U>);
//   Matrix<common_type_t<T, int>> res{};

//   transform(lhs.cbegin(), lhs.cend(), res.begin(),
//             [&rhs](const auto _lhs) { return _lhs + rhs; });

//   return res;
// }

template <typename T>
constexpr auto operator==(const Matrix<T> &lhs, const Matrix<T> &rhs) -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}

template <typename T>
constexpr auto operator!=(const Matrix<T> &lhs, const Matrix<T> &rhs) -> bool {
  return !(lhs == rhs);
}

} // namespace lin

#endif // __MATOPS_H__
