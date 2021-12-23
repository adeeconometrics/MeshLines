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

#include <iostream>
#include <vector>

using std::vector;

template <typename T> class Matrix {
private:
  vector<vector<T>> m_matrix;

public:
  Matrix(const vector<vector<T>> &_matrix) : m_matrix{_matrix} {}

  Matrix operator=(const Matrix &rhs) const;
  Matrix operator=(Matrix &&rhs) const;

  Matrix operator*(const Matrix &rhs) const;
  Matrix operator*(float scalar) const;
  Matrix operator+(const Matrix &rhs) const;
  Matrix operator+(float scalar) const;
  Matrix operator-(const Matrix &rhs) const;
  Matrix operator-(float scalar) const;
  Matrix operator/(const Matrix &rhs) const;
  Matrix operator/(float scalar) const;

  bool operator==(const Matrix &rhs) const noexcept;
  bool operator!=(const Matrix &rhs) const noexcept;

  bool is_square() const noexcept;
};

template <typename T>
Matrix<T> kron_prod(const Matrix<T> &lhs, const Matrix<T> &rhs) noexcept;

template <typename T> Matrix<T> det(const Matrix<T> &M) noexcept;
