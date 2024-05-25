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
#include "matfunc.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>
namespace lin {

using std::common_type_t;
using std::transform;
using std::vector;

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) + rhs(i, j);
    }
  }

  return result;
}

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) - rhs(i, j);
    }
  }

  return result;
}

template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U>, Rows, Cols> {
  Matrix<common_type_t<T, U>, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) * rhs(i, j);
    }
  }

  return result;
}
/**
 * @brief Returns the division of two matrices. Will likely return a matrix of
 * doubles for arithmetic types. NOTE: Implementation will be modified in the
 * future.
 *
 * @tparam T The type of the matrix
 * @tparam U The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param lhs The left hand side matrix
 * @param rhs The right hand side matrix
 * @return Matrix<common_type_t<T, U, double>, Rows, Cols> The division of the
 * two matrices.
 */
template <typename T, typename U = T, std::size_t Rows, std::size_t Cols>
constexpr auto operator/(const Matrix<T, Rows, Cols> &lhs,
                         const Matrix<U, Rows, Cols> &rhs)
    -> Matrix<common_type_t<T, U, double>, Rows, Cols> {

  using result_type = common_type_t<T, U, double>;
  Matrix<result_type, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = static_cast<result_type>(lhs(i, j)) /
                     static_cast<result_type>(rhs(i, j));
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator+=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      lhs(i, j) += rhs(i, j);
    }
  }

  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator-=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      lhs(i, j) -= rhs(i, j);
    }
  }
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator*=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      lhs(i, j) *= rhs(i, j);
    }
  }
  return lhs;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto
operator/=(Matrix<T, Rows, Cols> &lhs,
           const Matrix<T, Rows, Cols> &rhs) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      lhs(i, j) /= rhs(i, j);
    }
  }
  return lhs;
}
/**
 * @brief This implementation is inspired from Numpy's broadcasting feature and
 * is not a standard linear algebra operation. See:
 * https://numpy.org/doc/stable/user/basics.broadcasting.html.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param lhs The left hand side matrix
 * @param rhs The right hand side scalar
 * @return Matrix<T, Rows, Cols>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+(const Matrix<T, Rows, Cols> &lhs,
                         const vector<T> &rhs) -> Matrix<T, Rows, Cols> {

  assert(Cols == rhs.size());

  Matrix<T, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) + rhs[j];
    }
  }

  return result;
}
/**
 * @brief This implementation is inspired from Numpy's broadcasting feature and
 * is not a standard linear algebra operation. See:
 * https://numpy.org/doc/stable/user/basics.broadcasting.html.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param lhs The left hand side matrix
 * @param rhs The right hand side scalar
 * @return Matrix<T, Rows, Cols>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-(const Matrix<T, Rows, Cols> &lhs,
                         const vector<T> &rhs) -> Matrix<T, Rows, Cols> {

  assert(Cols == rhs.size());

  Matrix<T, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = lhs(i, j) - rhs[j];
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*(const Matrix<T, Rows, Cols> &lhs,
                         const vector<T> &rhs) -> vector<T> {

  assert(Cols == rhs.size());

  vector<T> result;
  result.reserve(Cols);

  for (std::size_t i{}; i < Rows; i++) {
    T sum{};
    for (std::size_t j{}; j < Cols; j++) {
      sum += lhs(i, j) * rhs[j];
    }
    result.emplace_back(sum);
  }

  return result;
}

/**
 * @brief This implementation is inspired from Numpy's broadcasting feature and
 * is not a standard linear algebra operation. This is modified for inplace
 * operations. See: https://numpy.org/doc/stable/user/basics.broadcasting.html.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param t_mat The left hand side matrix
 * @param t_vec The right hand side scalar
 * @return Matrix<T, Rows, Cols>& The division of the two matrices.
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+=(Matrix<T, Rows, Cols> &t_mat,
                          const vector<T> &t_vec) -> Matrix<T, Rows, Cols> & {
  assert(Cols == t_vec.size());
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      t_mat(i, j) += t_vec[j];
    }
  }
  return t_mat;
}
/**
 * @brief This implementation is inspired from Numpy's broadcasting feature and
 * is not a standard linear algebra operation. This is modified for inplace
 * operations. See: https://numpy.org/doc/stable/user/basics.broadcasting.html.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @param t_mat The left hand side matrix
 * @param t_vec The right hand side scalar
 * @return Matrix<T, Rows, Cols>& The division of the two matrices.
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-=(Matrix<T, Rows, Cols> &t_mat,
                          const vector<T> &t_vec) -> Matrix<T, Rows, Cols> & {
  assert(Cols == t_vec.size());
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      t_mat(i, j) -= t_vec[j];
    }
  }
  return t_mat;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+(const Matrix<T, Rows, Cols> &t_mat,
                         T t_scalar) -> Matrix<T, Rows, Cols> {
  Matrix<T, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = t_mat(i, j) + t_scalar;
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-(const Matrix<T, Rows, Cols> &t_mat,
                         T t_scalar) -> Matrix<T, Rows, Cols> {
  Matrix<T, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = t_mat(i, j) - t_scalar;
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*(const Matrix<T, Rows, Cols> &t_mat,
                         T t_scalar) -> Matrix<T, Rows, Cols> {
  Matrix<T, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = t_mat(i, j) * t_scalar;
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator/(const Matrix<T, Rows, Cols> &t_mat,
                         T t_scalar) -> Matrix<T, Rows, Cols> {
  Matrix<T, Rows, Cols> result{};

  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      result(i, j) = t_mat(i, j) / t_scalar;
    }
  }

  return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator+=(Matrix<T, Rows, Cols> &t_mat,
                          T t_scalar) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      t_mat(i, j) += t_scalar;
    }
  }
  return t_mat;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator-=(Matrix<T, Rows, Cols> &t_mat,
                          T t_scalar) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      t_mat(i, j) -= t_scalar;
    }
  }
  return t_mat;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator*=(Matrix<T, Rows, Cols> &t_mat,
                          T t_scalar) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      t_mat(i, j) *= t_scalar;
    }
  }
  return t_mat;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator/=(Matrix<T, Rows, Cols> &t_mat,
                          T t_scalar) -> Matrix<T, Rows, Cols> & {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      t_mat(i, j) /= t_scalar;
    }
  }
  return t_mat;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator==(const Matrix<T, Rows, Cols> &lhs,
                          const Matrix<T, Rows, Cols> &rhs) -> bool {
  for (std::size_t i{}; i < Rows; i++) {
    for (std::size_t j{}; j < Cols; j++) {
      if (lhs(i, j) != rhs(i, j))
        return false;
    }
  }
  return true;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto operator!=(const Matrix<T, Rows, Cols> &lhs,
                          const Matrix<T, Rows, Cols> &rhs) -> bool {
  return !(lhs == rhs);
}

} // namespace lin

#endif // __MATOPS_H__
