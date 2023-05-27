#ifndef __MATRIX_H__
#define __MATRIX_H__

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
template <typename T> constexpr auto zero_mat(size_t n) -> Matrix<T> {
  static_assert(std::is_arithmetic_v<T>);

  return Matrix<T>{n, vector<T>(n, 0)};
}
/**
 * @brief Returns a rectangular zero matrix
 *
 * @tparam T
 * @param n
 * @param m
 * @return Matrix<T>
 */
template <typename T>
constexpr auto zero_mat(std::size_t n, std::size_t m) -> Matrix<T> {
  static_assert(std::is_arithmetic_v<T>);

  return Matrix<T>(n, vector<T>(m, 0));
}

/**
 * @brief Constructor for identity matrix
 *
 * @tparam T
 * @param n
 * @return Matrix<T>
 */
template <typename T> constexpr auto id(size_t n) -> Matrix<T> {
  static_assert(std::is_arithmetic_v<T>);

  Matrix<T> res = zero_mat<T>(n);

  for (std::size_t i{}; i < n; i++) {
    res[i][i] = 1;
  }
  return res;
}

} // namespace lin
#endif // __MATRIX_H__