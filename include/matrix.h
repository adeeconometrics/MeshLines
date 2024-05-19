#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "../include/vecops.h"

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <type_traits>
#include <vector>

namespace lin {

template <typename T, std::size_t Rows, std::size_t Cols> class Matrix {
public:
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

public:
  Matrix() { m_data.reserve(Rows * Cols); }

  explicit Matrix(const std::vector<T> &t_vec, std::size_t t_rows,
                  std::size_t t_cols)
      : m_data(t_vec) {
    static_assert(t_vec.size() == t_rows * t_cols, "Size mismatch");
    static_assert(t_vec.size() == Rows * Cols, "Size mismatch");
  }

  Matrix(std::initializer_list<std::initializer_list<T>> t_list) {
    if (t_list.size() != Rows) {
      throw std::invalid_argument("Invalid number of rows in initializer list");
    }

    auto data_iter = m_data.begin();
    for (const auto &row_list : t_list) {
      if (row_list.size() != Cols) {
        throw std::invalid_argument("Invalid row size in initializer list");
      }

      std::copy(row_list.cbegin(), row_list.cend(), data_iter);
      data_iter += Cols;
    }
  }

  constexpr auto operator()(std::size_t i, std::size_t j) -> T {
    return m_data[i * Cols + j];
  }

  constexpr auto operator()(std::size_t i, std::size_t j) const -> T {
    return m_data[i * Cols + j];
  }

  constexpr auto begin() noexcept -> iterator { return m_data.begin(); }
  constexpr auto end() noexcept -> iterator { return m_data.end(); }
  constexpr auto cbegin() const noexcept -> const_iterator {
    return m_data.cbegin();
  }
  constexpr auto cend() const noexcept -> const_iterator {
    return m_data.cend();
  }

private:
  std::vector<T> m_data;
};

/**
 * @brief Constructor for zero matrix.
 *
 * @tparam T
 * @param n
 * @return Matrix<T>
 */
template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto zero_mat() -> Matrix<T, Rows, Cols> {
  return {std::vector<T>(Rows * Cols, 0), Rows, Cols};
}

/**
 * @brief Constructor for identity matrix.
 *
 * @tparam T
 * @param n
 * @return Matrix<T>
 */
template <typename T, std::size_t N,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto id() -> Matrix<T, N, N> {

  Matrix<T, N, N> res = zero_mat<T, N, N>();

  for (std::size_t i = 0; i < N; i++) {
    res(i, i) = 1;
  }
  return res;
}

template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto scalar_id(T t_value) -> Matrix<T, Rows, Cols> {

  Matrix<T, Rows, Cols> res = zero_mat<T, Rows, Cols>();

  for (std::size_t i = 0; i < N; i++) {
    res(i, i) = t_value;
  }
  return res;
}

template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto scalar_mat(T t_value) -> Matrix<T, Rows, Cols> {
  return {std::vector<T>(Rows * Cols, t_value), Rows, Cols};
}
// col_vec
// template <typename T, std::size_t Rows, std::size_t Cols>
// constexpr auto col_vec(std::size_t t_col) -> std::vector<T> {
//   std::vector<T> result;
//   result.reserve(Rows);
//   for (std::size_t i = 0; i < Rows; i++) {
//     result.emplace_back(i == t_col ? 1 : 0);
//   }
//   return result;
// }

// row_vec

template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto
operator-(const Matrix<T, Rows, Cols> &M) -> Matrix<T, Rows, Cols> {
  Matrix<T, Rows, Cols> result;
  result.reserve(M.size());

  std::transform(std::cbegin(M), std::cend(M), std::back_inserter(result),
                 [](const auto element) { return -element; });

  return result;
}
} // namespace lin
#endif // __MATRIX_H__