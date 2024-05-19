#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "../include/vecops.h"

#include <algorithm>
#include <type_traits>
#include <vector>

namespace lin {

template <typename T> using Matrix = std::vector<std::vector<T>>;

/**
 * @brief Constructor for zero matrix (note that this returns a square matrix).
 *
 * @tparam T
 * @param n
 * @return Matrix<T>
 */
template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto zero_mat(std::size_t n) -> Matrix<T> {
  return Matrix<T>{n, std::vector<T>(n, 0)};
}
/**
 * @brief Returns a rectangular zero matrix
 *
 * @tparam T
 * @param n
 * @param m
 * @return Matrix<T>
 */
template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto zero_mat(std::size_t n, std::size_t m) -> Matrix<T> {
  return Matrix<T>(n, std::vector<T>(m, 0));
}

/**
 * @brief Constructor for identity matrix
 *
 * @tparam T
 * @param n
 * @return Matrix<T>
 */
template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto id(std::size_t n) -> Matrix<T> {

  Matrix<T> res = zero_mat<T>(n);

  for (std::size_t i = 0; i < n; i++) {
    res[i][i] = 1;
  }
  return res;
}

// template <typename T> constexpr auto id(std::size_t n, std::size_t m) ->
// Matrix<T> {
//   static_assert(std::is_arithmetic_v<T>);

//   Matrix<T> res = zero_mat(n,m);

//   for(std::size_t i = 0; i < n; i++){
//         res[i][i] = 0;
//     }
// }

template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto scalar_mat(std::size_t n, T t_value) -> Matrix<T> {

  Matrix<T> res = zero_mat<T>(n);

  for (std::size_t i = 0; i < n; i++) {
    res[i][i] = t_value;
  }
  return res;
}

// col_vec

// row_vec

template <typename T>
constexpr auto operator-(const Matrix<T> &M) -> Matrix<T> {
  Matrix<T> result;
  result.reserve(M.size());

  std::transform(std::cbegin(M), std::cend(M), std::back_inserter(result),
                 [](const auto element) { return -element; });

  return result;
}
} // namespace lin
#endif // __MATRIX_H__