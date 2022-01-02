/**
 * @file Vector.h
 * @author ddamiana
 * @brief representation of a mathematical vector
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
// #include "Matrix.h"
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <type_traits>
#include <vector>

using std::vector;

template <typename T> class Vector {
public:
  typedef vector<T> value_type;
  typedef value_type &reference;
  typedef value_type *pointer_type;
  typedef const value_type const_type;
  typedef const reference const_reference;
  typedef const pointer_type const_pointer;

  typedef value_type::iterator iterator;
  typedef const iterator const_iterator;
  typedef value_type::reverse_iterator reverse_iterator;
  typedef const reverse_iterator const_reverse_iterator;

private:
  value_type m_vector;

public:
  Vector(const_reference _vector) : m_vector(_vector) {}
  Vector(std::initializer_list<T> _list) : m_vector((value_type)_list) {}
  // Vector(const Matrix& m);

  // Vector operator=(const Vector &rhs) const;
  // Vector operator=(Vector &&rhs) const;
  // explicit operator Matrix();

  Vector operator*(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    value_type result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a * b; });

    return result;
  }

  Vector operator*(float scalar) const noexcept {
    value_type result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a * scalar; });
    return result;
  }

  Vector operator+(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    value_type result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a + b; });

    return result;
  }

  Vector operator+(float scalar) const noexcept {
    value_type result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a + scalar; });
    return result;
  }

  Vector operator-(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    value_type result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a - b; });

    return result;
  }

  Vector operator-(float scalar) const noexcept {
    value_type result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a - scalar; });
    return result;
  }

  Vector operator/(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    value_type result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a / b; });

    return result;
  }

  Vector operator/(float scalar) const noexcept {
    value_type result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a / scalar; });
    return result;
  }

  bool operator==(const Vector &rhs) const noexcept {
    for (size_t i = 0; i < m_vector.size(); i++) {
      if (m_vector[i] != rhs.m_vector[i])
        return false;
    }
    return true;
  }

  bool operator!=(const Vector &rhs) const noexcept { return !(*this == rhs); }

  unsigned int size() const { return m_vector.size(); }

  template <typename U>
  friend std::ostream &operator<<(std::ostream &ss, const Vector<U> &v);

  iterator begin() const noexcept {return m_vector.begin(); }
  iterator end() const noexcept { return m_vector.end(); }
  const_iterator cbegin() const noexcept { return m_vector.cbegin(); }
  const_iterator cend() const noexcept { return m_vector.cend(); }
  reverse_iterator begin() const noexcept { return m_vector.rbegin(); }
  reverse_iterator end() const noexcept { return m_vector.rend(); }
};

template <typename T>
std::ostream &operator<<(std::ostream &ss, const Vector<T> &v) {
  ss << "[ ";
  for (auto const &i : v.m_vector) {
    ss << i << " ";
  }
  return ss << "]";
};

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value>::type>
double dist(const Vector<T> &v) {
  double result{};
  for (auto const &i : v)
    result += pow(i, 2);

  return sqrt(result);
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value>::type>
double lp_norm(const Vector<T> &v, float p) {
  if (p == 0.0)
    return std::count_if(v.begin(), v.end(),
                        [](const T &x) -> bool { return x != 0; });
  if (std::isinf(p)) {
    double max = *std::max_element(
        v.begin(), v.end(),
        [](const T &a, const T &b) -> bool { std::abs(a) < std::abs(b); });

    return std::abs(max);
  }

  double result{};
  for (const T &i : v)
    result += pow(i, p);
  return pow(result, 1 / p);
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value>::type>
double sum(const Vector<T>& v) {
  double result{};
  for(const T& i: v)
    result += i;

  return result;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value>::type>
double prod(const Vector<T> &v) {
  double result{};
  for (const T &i : v)
    result *= i;

  return result;
}