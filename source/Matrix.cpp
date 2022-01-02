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
#include <vector>
#include "Vector.h"


using std::vector;

template <typename T> class Matrix {
public:
  typedef vector<vector<T>> value_type;
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

  Matrix(const std::initializer_list<vector<T>>& _list): m_row(_list.size()){

  }

  Matrix(const Matrix::value_type& _mat):m_row(_mat.row()), 
                                         m_col(_mat.col()),
                                         m_matrix(_mat.m_matrix){}

  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) = default;
  ~Matrix() = default;

  Matrix operator=(const Matrix &rhs) const;
  Matrix operator=(Matrix &&rhs) const;

  Matrix& operator*=(const Matrix &rhs) const;
  Matrix& operator*=(float scalar) const;
  Matrix& operator+=(const Matrix &rhs) const;
  Matrix& operator+=(float scalar) const;
  Matrix& operator-=(const Matrix &rhs) const;
  Matrix& operator-=(float scalar) const;
  Matrix& operator/=(const Matrix &rhs) const;
  Matrix& operator/=(float scalar) const;

  Matrix operator*(const Matrix &lhs, const Matrix &rhs) const;

  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  Matrix operator*(const Matrix &lhs, U scalar) const;

  Matrix operator+(const Matrix &lhs, const Matrix &rhs) const;
  
  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  Matrix operator+(const Matrix &lhs, U scalar) const;
  
  Matrix operator-(const Matrix &lhs, const Matrix &rhs) const;

  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  Matrix operator-(const Matrix &lhs, U scalar) const;
  
  Matrix operator/(const Matrix &lhs, const Matrix &rhs) const;
  
  template <typename U,
            typename std::enable_if<std::is_arithmetic<T>::value, bool> = true>
  Matrix operator/(const Matrix &lhs, U scalar) const;

  bool operator==(const Matrix &rhs) const noexcept;
  bool operator!=(const Matrix &rhs) const noexcept;

  bool is_square() const noexcept;
  bool is_invertible() const noexcept;

  size_t row();
  size_t col();

  vector<T> get_row_at(size_t idx);
  vector<T> get_col_at(size_t idx);
};


template <typename T>
Matrix<T> kron_prod(const Matrix<T> &lhs, const Matrix<T> &rhs) noexcept;

template <typename T> Matrix<T> det(const Matrix<T> &M) noexcept;

template <typename T> Matrix<T> inv(const Matrix<T> &M);

template <typename T> Matrix<T> transpose(const Matrix<T> &M) noexcept;

template <typename T> T trace(const Matrix<T> &M);

template <typename T> Matrix<T> identity(int ndims = 2) noexcept;
