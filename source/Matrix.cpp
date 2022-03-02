/**
 * @file Matrix.cpp
 * @author ddamiana
 * @brief mathematical representation of a matrix
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#include <algorithm>
#include <iostream>
#include <type_traits>
#include <initializer_list>
#include "Vector.h"


template <typename T> class Matrix {
public:
  typedef Vector<Vector<T>> value_type;
  typedef value_type& reference;
  typedef const value_type const_type;
  typedef const reference const_reference;
  typedef value_type* pointer_type;
  typedef const pointer_type const_pointer;

private:
  value_type m_matrix;
  size_t m_row{}, m_col{};

public:
  
  Matrix() = default;
  Matrix(size_t row, size_t col): m_row(row), m_col(col){

  }

  Matrix(const std::initializer_list<Vector<T>>& _list): m_row(_list.size()){

  }

  Matrix(const Matrix::value_type& _mat):m_row(_mat.row()), 
                                         m_col(_mat.col()),
                                         m_matrix(_mat.m_matrix){}

  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) = default;
  virtual ~Matrix() = default;

  auto operator=(const Matrix &rhs) const -> Matrix& = default;
  auto operator=(Matrix &&rhs) const -> Matrix& = default;

  auto operator*=(const Matrix &rhs) const -> Matrix&;
  auto operator*=(float scalar) const -> Matrix&;
  auto operator+=(const Matrix &rhs) const -> Matrix&;
  auto operator+=(float scalar) const -> Matrix&;
  auto operator-=(const Matrix &rhs) const -> Matrix&;
  auto operator-=(float scalar) const -> Matrix&;
  auto operator/=(const Matrix &rhs) const -> Matrix&;
  auto operator/=(float scalar) const -> Matrix&;

  auto operator*(const Matrix &lhs, const Matrix &rhs) const -> Matrix;

  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  auto operator*(const Matrix &lhs, U scalar) const -> Matrix;

  auto operator+(const Matrix &lhs, const Matrix &rhs) const -> Matrix;
  
  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  auto operator+(const Matrix &lhs, U scalar) const -> Matrix;
  
  auto operator-(const Matrix &lhs, const Matrix &rhs) const -> Matrix;

  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  auto operator-(const Matrix &lhs, U scalar) const -> Matrix;
  
  auto operator/(const Matrix &lhs, const Matrix &rhs) const -> Matrix;
  
  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  auto operator/(const Matrix &lhs, U scalar) const -> Matrix;

  auto operator==(const Matrix &rhs) const noexcept -> bool;
  auto operator!=(const Matrix &rhs) const noexcept -> bool;

  bool is_square() const noexcept;
  bool is_invertible() const noexcept;

  auto row() const noexcept -> size_t;
  auto col() const noexcept -> size_t;

  auto get_row_vec(size_t idx) -> vector<T>;
  auto get_col_vec(size_t idx) -> vector<T>;

  auto end() -> decltype(m_matrix.end()) {return m_matrix.end(); }
  auto cend() -> decltype(m_matrix.cend()) { return m_matrix.cend(); }
  auto rend() -> decltype(m_matrix.rend()) { return m_matrix.rend(); }

  auto begin() -> decltype(m_matrix.begin()) {return m_matrix.begin();}
  auto cbegin() -> decltype(m_matrix.cbegin()) { return m_matrix.cbegin(); }
  auto rbegin() -> decltype(m_matrix.rbegin()) { return m_matrix.rbegin(); }
};


template <typename T>
Matrix<T> kron_prod(const Matrix<T> &lhs, const Matrix<T> &rhs) noexcept;

template <typename T> Matrix<T> det(const Matrix<T> &M) noexcept;

template <typename T> Matrix<T> inv(const Matrix<T> &M);

template <typename T> Matrix<T> transpose(const Matrix<T> &M) noexcept;

template <typename T> T trace(const Matrix<T> &M);

template <typename T> Matrix<T> identity(int ndims = 2) noexcept;
