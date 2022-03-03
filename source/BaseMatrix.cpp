/**
 * @file BaseMatrix.cpp
 * @author ddamiana
 * @brief base class for BaseMatrix types
 * @version 0.1
 * @date 2022-02-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Vector.h"
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <stcexcept>

template <typename T> class BaseMatrix {
public:
  typedef vector<vector<T>> value_type;
  typedef value_type &reference;
  typedef const value_type const_type;
  typedef const reference const_reference;
  typedef value_type *pointer_type;
  typedef const pointer_type const_pointer;

private:
  value_type m_vector;
  size_t m_row{}, m_col{};

public:
  BaseMatrix() = default;
  BaseMatrix(const BaseMatrix &) = delete;
  BaseMatrix(BaseMatrix &&) = delete;
  virtual ~BaseMatrix() = default;

  auto operator=(const BaseMatrix &) -> BaseMatrix & = delete;
  auto operator=(BaseMatrix &&) -> BaseMatrix & = delete;

  virtual auto operator*=(float scalar) const -> BaseMatrix & = 0;
  virtual auto operator*=(const BaseMatrix &rhs) const -> BaseMatrix & = 0;
  virtual auto operator+=(const BaseMatrix &rhs) const -> BaseMatrix & = 0;
  virtual auto operator+=(float scalar) const -> BaseMatrix & = 0;
  virtual auto operator-=(const BaseMatrix &rhs) const -> BaseMatrix & = 0;
  virtual auto operator-=(float scalar) const -> BaseMatrix & = 0;
  virtual auto operator/=(const BaseMatrix &rhs) const -> BaseMatrix & = 0;
  virtual auto operator/=(float scalar) const -> BaseMatrix & = 0;

  virtual auto operator==(const BaseMatrix &rhs) const noexcept -> bool;
  virtual auto operator!=(const BaseMatrix &rhs) const noexcept -> bool;

  auto row() const noexcept -> size_t { return m_row; }
  auto col() const noexcept -> size_t { return m_col; }

  auto get_row_vec(size_t idx) -> vector<T>;
  auto get_col_vec(size_t idx) -> vector<T>;

  auto end() -> decltype(m_matrix.end()) { return m_matrix.end(); }
  auto cend() -> decltype(m_matrix.cend()) { return m_matrix.cend(); }
  auto rend() -> decltype(m_matrix.rend()) { return m_matrix.rend(); }

  auto begin() -> decltype(m_matrix.begin()) { return m_matrix.begin(); }
  auto cbegin() -> decltype(m_matrix.cbegin()) { return m_matrix.cbegin(); }
  auto rbegin() -> decltype(m_matrix.rbegin()) { return m_matrix.rbegin(); }
};