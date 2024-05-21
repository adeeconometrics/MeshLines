#ifndef __MATOPS_H__
#define __MATOPS_H__

/**
 * @file matops.h
 * @author ddamiana
 * @brief Contains operator overload for Matrix, Vector, and Scalar types. Note:
 * These methods will soon be deprecated in favor of using Expression Templates.
 * @version 0.1
 * @date 2023-05-27
 *
 * @copyright Copyright (c) 2023
 *
 */

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

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const T &a, const U &b) { return a + b; });

  return result;
}

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const T &a, const U &b) { return a - b; });

  return result;
}

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const T &a, const U &b) { return a * b; });

  return result;
}

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator/(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const T &a, const U &b) { return a / b; });

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator+=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  std::transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
                 [](auto &_lhs, auto &_rhs) { return _lhs += _rhs; });
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator-=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  std::transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
                 [](auto &_lhs, auto &_rhs) { return _lhs -= _rhs; });
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator*=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  std::transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
                 [](auto &_lhs, auto &_rhs) { return _lhs *= _rhs; });
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator/=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  std::transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
                 [](auto &_lhs, auto &_rhs) { return _lhs /= _rhs; });
  return lhs;
}

// template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
// constexpr auto operator+(const Matrix<T, Rows, Cols> &lhs, const vector<U>
// &rhs)
//     -> Matrix<common_type_t<T, U>, Rows, Cols> {

//   assert(Cols == rhs.size());

//   using result_type = common_type_t<T, U>;

//   Matrix<result_type, Rows, Cols> result{};

//   for (std::size_t i{}; i < Rows; i++) {
//     for (std::size_t j{}; j < Cols; j++) {
//       result(i, j) = lhs(i, j) + rhs[j];
//     }
//   }

//   return result;
// }

// template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
// constexpr auto operator-(const Matrix<T, Rows, Cols> &lhs, const vector<U>
// &rhs)
//     -> Matrix<common_type_t<T, U>, Rows, Cols> {

//   assert(Cols == rhs.size());

//   using result_type = common_type_t<T, U>;

//   Matrix<result_type, Rows, Cols> result{};

//   for (std::size_t i{}; i < Rows; i++) {
//     for (std::size_t j{}; j < Cols; j++) {
//       result(i, j) = lhs(i, j) - rhs[j];
//     }
//   }

//   return result;
// }

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*(const Matrix<T, Rows, Cols> &lhs, const vector<U> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {

  assert(Cols == rhs.size());

  using result_type = common_type_t<T, U>;

  Matrix<result_type, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) * rhs[j];
    }
  }

  return result;
}
// MAKE TYPES FUSE WITH FLOATING POINT
template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator/(const Matrix<T, Rows, Cols> &lhs, const vector<U> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {

  assert(Cols == rhs.size());

  // using result_type = common_type_t<T, U, double>;
  using result_type = common_type_t<T, U>;

  Matrix<result_type, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) / rhs[j];
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+=(Matrix<T, Rows, Cols> &lhs,
                          const vector<T> &rhs) -> Matrix<T, Rows, Cols> & {
  transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs += _rhs; });
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-=(Matrix<T, Rows, Cols> &lhs,
                          const vector<T> &rhs) -> Matrix<T, Rows, Cols> & {
  transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs -= _rhs; });
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*=(Matrix<T, Rows, Cols> &lhs,
                          const vector<T> &rhs) -> Matrix<T, Rows, Cols> & {
  transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs *= _rhs; });
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator/=(Matrix<T, Rows, Cols> &lhs,
                          const vector<T> &rhs) -> Matrix<T, Rows, Cols> & {
  transform(lhs.begin(), lhs.end(), rhs.cbegin(), lhs.begin(),
            [](auto &_lhs, auto &_rhs) { return _lhs /= _rhs; });
  return lhs;
}

// template <typename T, std::size_t Rows, std::size_t Cols>
// constexpr auto operator+(const Matrix<T, Rows, Cols> &lhs, int rhs)
//     -> Matrix<common_type_t<T, int>> {

//   // static_assert(std::is_arithmetic_v<U>);
//   Matrix<common_type_t<T, int>> res{};

//   transform(lhs.cbegin(), lhs.cend(), res.begin(),
//             [&rhs](const auto _lhs) { return _lhs + rhs; });

//   return res;
// }

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator==(const Matrix<T, Rows, Cols> &lhs,
                          const Matrix<T, Rows, Cols> &rhs) -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator!=(const Matrix<T, Rows, Cols> &lhs,
                          const Matrix<T, Rows, Cols> &rhs) -> bool {
  return !(lhs == rhs);
}

} // namespace lin

#endif // __MATOPS_H__
