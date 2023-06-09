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

#include <initializer_list>
#include <iostream>
#include <utility>
#include <vector>

using std::cout;
using std::initializer_list;
using std::pair;
using std::vector;

template <typename T> using iterator = typename vector<vector<T>>::iterator;
template <typename T>
using const_iterator = typename vector<vector<T>>::const_iterator;

template <typename T> struct BaseMatrix {
public:
  BaseMatrix(std::size_t m, std::size_t n)
      : m_matrix{vector<vector<T>>(m, vector<T>(n))}, m_row(m), m_col(n) {}

  constexpr auto at(std::size_t t_row, std::size_t t_col) noexcept -> T & {
    return m_matrix[t_row][t_col];
  }

  constexpr auto at(std::size_t t_row, std::size_t t_col) const noexcept
      -> T & {
    return m_matrix[t_row][t_col];
  }

  constexpr auto dims() const noexcept -> pair<std::size_t, std::size_t> {
    return {m_row, m_col};
  }

  auto begin() noexcept -> iterator<T> { return m_matrix.begin(); }
  auto end() noexcept -> iterator<T> { return m_matrix.end(); }
  auto cbegin() const noexcept -> const_iterator<T> { return m_matrix.cend(); }
  auto cend() const noexcept -> const_iterator<T> { return m_matrix.cend(); }

private:
  vector<vector<T>> m_matrix;
  std::size_t m_row;
  std::size_t m_col;
};

template <typename T> struct SquareMatrix : public BaseMatrix<T> {
  SquareMatrix(std::size_t m) : BaseMatrix<T>(m, m) {}
  // use initializer_list as basis of value initialization
};

template <typename T> struct DiagonalMatrix : public SquareMatrix<T> {
public:
  DiagonalMatrix(std::size_t m, const initializer_list<T> &v)
      : SquareMatrix<T>(m) {
    assert(v.size() == m);
  }
};

template <typename T> struct RectangularMatrix : public BaseMatrix<T> {
  RectangularMatrix(std::size_t m, std::size_t n) : BaseMatrix<T>(m, n) {}
};

auto main() -> int {
  SquareMatrix<int> M(3);
  for (const auto &row : M) {
    for (const auto &i : row)
      cout << i << " ";
    cout << '\n';
  }

  return 0;
}