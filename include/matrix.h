#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "../include/vecops.h"

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <vector>

namespace lin {

template <typename T, std::size_t Rows, std::size_t Cols> class Matrix {
public:
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

public:
  Matrix() { m_data.reserve(Rows * Cols); }

  Matrix(const std::vector<T> &t_vec, std::size_t t_rows, std::size_t t_cols)
      : m_data(t_vec) {
    assert(t_vec.size() == t_rows * t_cols); // Size mismatch
    assert(t_vec.size() == Rows * Cols);     // Size mismatch
  }

  Matrix(std::initializer_list<std::initializer_list<T>> t_list) {
    if (t_list.size() != Rows) {
      throw std::invalid_argument("Invalid number of rows in initializer list");
    }

    for (const auto &row : t_list) {
      if (row.size() != Cols) {
        throw std::invalid_argument(
            "Invalid number of columns in initializer list");
      }
      for (const auto &elem : row) {
        m_data.emplace_back(elem);
      }
    }
  }

  constexpr auto operator()(std::size_t i, std::size_t j) -> T & {
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

  auto rows() const noexcept -> std::size_t { return Rows; }
  auto cols() const noexcept -> std::size_t { return Cols; }
  auto dims() const noexcept -> std::tuple<std::size_t, std::size_t> {
    return {Rows, Cols};
  }

  auto empty() const noexcept -> bool { return m_data.empty(); }

  // auto push_back(T t_value) -> void { m_data.push_back(t_value); }
  // auto emplace_back(T t_value) -> void { m_data.emplace_back(t_value); }

private:
  std::vector<T> m_data;
};

/**
 * @brief Constructor for zero matrix which pertains to T{} constructor for
 * compound objects.
 *
 * @tparam T The type of the matrix
 * @tparam Rows The number of rows
 * @tparam Cols The number of columns
 * @return Matrix<T, Rows, Cols>
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto zero_mat() -> Matrix<T, Rows, Cols> {
  return {std::vector<T>(Rows * Cols, T{}), Rows, Cols};
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

  for (std::size_t i = 0; i < Rows; i++) {
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

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < Cols; j++) {
      result(i, j) = -M(i, j);
    }
  }

  return result;
}
} // namespace lin
#endif // __MATRIX_H__