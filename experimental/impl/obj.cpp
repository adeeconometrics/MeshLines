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
#include <functional>
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

template <typename T> struct Matrix {
public:
  using value_type = T;
  using reference_type = T &;
  using pointer_type = T *;

  using iterator = typename vector<vector<T>>::iterator;
  using const_iterator = typename vector<vector<T>>::const_iterator;

public:
  Matrix() : m_row(0), m_col(0) {}
  Matrix(const Matrix<T> &) = default;
  Matrix(Matrix<T> &&) = default;

  Matrix(std::size_t m, std::size_t n)
      : m_row(m), m_col(n), m_matrix(m_row, vector<T>(m_col)) {}

  Matrix(std::size_t m, std::size_t n, const vector<vector<T>> &t_matrix)
      : m_row(m), m_col(n), m_matrix(t_matrix) {}

  Matrix(initializer_list<initializer_list<T>> t_list)
      : m_row(t_list.size()), m_col((t_list.begin())->size()),
        m_matrix(vector<vector<T>>(t_list.begin(), t_list.end())) {}

  virtual ~Matrix() = default;

  auto operator=(const Matrix<T> &) -> Matrix<T> & = default;
  auto operator=(Matrix<T> &&) -> Matrix<T> & = default;

  auto dims() const noexcept -> pair<std::size_t, std::size_t> {
    return {m_row, m_col};
  }

  auto at(std::size_t t_row, std::size_t t_col) noexcept -> T & {
    return m_matrix[t_row][t_col];
  }
  auto at(std::size_t t_row, std::size_t t_col) const noexcept -> const T & {
    return m_matrix[t_row][t_col];
  }

  auto begin() noexcept -> iterator { return std::begin(m_matrix); }

  auto end() noexcept -> iterator { return std::end(m_matrix); }

  auto cbegin() const noexcept -> const_iterator {
    return std::cbegin(m_matrix);
  }

  auto cend() const noexcept -> const_iterator { return std::cend(m_matrix); }

private:
  std::size_t m_row;
  std::size_t m_col;

  vector<vector<T>> m_matrix;
};

template <typename T> struct SquareMatrix : public Matrix<T> {
  SquareMatrix(std::size_t m) : Matrix<T>(m, m) {}

  SquareMatrix(initializer_list<initializer_list<T>> t_list)
      : Matrix<T>(t_list) {
    assert(all_same_size(t_list));
  }
};

// template <typename T> struct DiagonalMatrix : public SquareMatrix<T> {
// public:
//   DiagonalMatrix(std::size_t m, initializer_list<T> v) : SquareMatrix<T>(m) {
//     assert(v.size() == m);
//   }
// };

template <typename T> struct RectMatrix : public Matrix<T> {
  RectMatrix(std::size_t m, std::size_t n) : Matrix<T>(m, n) {}
  RectMatrix(initializer_list<initializer_list<T>> t_list)
      : Matrix<T>(t_list) {}
};

template <typename T>
constexpr auto operator+(const vector<T> &lhs, const vector<T> &rhs)
    -> vector<T> {
  assert(lhs.size() == rhs.size());

  vector<T> result;
  result.reserve(lhs.size());

  std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs),
                 std::back_inserter(result), std::plus<T>());
  return result;
}

template <typename T>
constexpr auto operator+(const Matrix<T> &lhs, const Matrix<T> &rhs)
    -> Matrix<T> {

  assert(lhs.dims() == rhs.dims());
  const std::size_t row = lhs.dims().first;
  const std::size_t col = lhs.dims().second;

  vector<vector<T>> result;

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), std::begin(result),
                 [](const vector<T> &_lhsv, const vector<T> &_rhsv) {
                   return _lhsv + _rhsv;
                 });
  return {row, col, result};
}

auto main() -> int {
  SquareMatrix<int> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  auto B = M + M;
  cout << "something\n";
  for (auto row : B) {
    for (auto i : row)
      cout << i << " ";
    cout << '\n';
  }

  return 0;
}