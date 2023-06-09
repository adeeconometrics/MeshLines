/**
 * @file obj.cpp
 * @author ddamiana
 * @brief experimental implementation of object heirarchy
 * @version 0.1
 * @date 2023-06-08
 *
 * @copyright Copyright (c) 2023
 *
 */

// todo
// Do subtype instantiation wthout calling base class constructor
// -- Explore the cost and benefit of creating independent types
// -- Explore the possibility of moving the has-a relationship to concrete types
// --- Note that you have to implement iterator interface in the concrete types
// ---- Base class can flesh out these expectations but vfunctions is not cheap
// Check if stl algorithms works on BaseClass with begin, cbegin
// iterators

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <vector>

using std::cout;
using std::initializer_list;
using std::pair;
using std::vector;

/**
 * @brief To check for equal length in the initializer list at compile-time.
 *
 * Note: make this function run on compile type and write the
 * assertion as static_assert.
 *
 * @return true when all size of the initializer_list are equal
 * @return false when not all size of the initializer_list are equal
 */
template <typename T>
constexpr auto all_same_size(initializer_list<initializer_list<T>> t_list)
    -> bool {
  const std::size_t m = t_list.size();
  return std::all_of(t_list.begin(), t_list.end(),
                     [m](const auto i) { return i.size() == m; });
}

template <typename T> using iterator = typename vector<vector<T>>::iterator;
template <typename T>
using const_iterator = typename vector<vector<T>>::const_iterator;

template <typename T> struct BaseMatrix {

public:
  BaseMatrix() = delete;
  BaseMatrix(const BaseMatrix<T> &) = default;
  BaseMatrix(BaseMatrix<T> &&) = default;

  BaseMatrix(std::size_t m, std::size_t n) : m_row(m), m_col(n) {}
  virtual ~BaseMatrix() = default;

  auto dims() const noexcept -> pair<std::size_t, std::size_t> {
    return {m_row, m_col};
  }

  virtual auto at(std::size_t t_row, std::size_t t_col) noexcept -> T & = 0;
  virtual auto at(std::size_t t_row, std::size_t t_col) const noexcept
      -> const T & = 0;

  virtual auto begin() noexcept -> iterator<T> = 0;
  virtual auto end() noexcept -> iterator<T> = 0;
  virtual auto cbegin() const noexcept -> const_iterator<T> = 0;
  virtual auto cend() const noexcept -> const_iterator<T> = 0;

private:
  std::size_t m_row;
  std::size_t m_col;
};

template <typename T> struct SquareMatrix : public BaseMatrix<T> {
  SquareMatrix(std::size_t m)
      : BaseMatrix<T>(m, m), m_matrix(m, vector<T>(m)) {}

  SquareMatrix(initializer_list<initializer_list<T>> t_list)
      : BaseMatrix<T>(t_list.size(), t_list.size()),
        m_matrix(vector<vector<T>>(t_list.begin(), t_list.end())) {
    assert(all_same_size(t_list));
  }

  auto at(std::size_t t_row, std::size_t t_col) noexcept -> T & override {
    return m_matrix[t_row][t_col];
  }

  auto at(std::size_t t_row, std::size_t t_col) const noexcept
      -> const T & override {
    return m_matrix[t_row][t_col];
  }

  auto begin() noexcept -> iterator<T> override { return m_matrix.begin(); }

  auto end() noexcept -> iterator<T> override { return m_matrix.end(); }

  auto cbegin() const noexcept -> const_iterator<T> override {
    return m_matrix.cbegin();
  }
  auto cend() const noexcept -> const_iterator<T> override {
    return m_matrix.cend();
  }

private:
  vector<vector<T>> m_matrix;
};

// template <typename T> struct DiagonalMatrix : public SquareMatrix<T> {
// public:
//   DiagonalMatrix(std::size_t m, const initializer_list<T> &v)
//       : SquareMatrix<T>(m) {
//     assert(v.size() == m);
//   }
// };

// template <typename T> struct RectangularMatrix : public BaseMatrix<T> {
//   RectangularMatrix(std::size_t m, std::size_t n) : BaseMatrix<T>(m, n) {}
// };

auto main() -> int {
  SquareMatrix<int> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  for (auto row : M) {
    for (auto i : row)
      cout << i << " ";
    cout << '\n';
  }

  return 0;
}