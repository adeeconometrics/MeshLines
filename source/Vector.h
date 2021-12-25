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
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <type_traits>
#include <math.h>
#include <stdexcept>
#include <vector>

using std::vector;

template <typename T> class Vector {
private:
  vector<T> m_vector;

public:
  Vector(const vector<T> &_vector) : m_vector(_vector) {}
  Vector(std::initializer_list<T> _list) : m_vector((vector<T>)_list) {}

  // Vector operator=(const Vector &rhs) const;
  // Vector operator=(Vector &&rhs) const;

  Vector operator*(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    vector<T> result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a * b; });

    return result;
  }

  Vector operator*(float scalar) const noexcept {
    vector<T> result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a * scalar; });
    return result;
  }

  Vector operator+(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    vector<T> result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a + b; });

    return result;
  }

  Vector operator+(float scalar) const noexcept {
    vector<T> result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a + scalar; });
    return result;
  }

  Vector operator-(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    vector<T> result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a - b; });

    return result;
  }

  Vector operator-(float scalar) const noexcept {
    vector<T> result(m_vector.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                   [&scalar](const T &a) -> T { return a - scalar; });
    return result;
  }

  Vector operator/(const Vector &rhs) const {
    if (size() != rhs.size()) {
      throw std::domain_error(
          "size of vector is expected to have the same size");
    }

    vector<T> result(rhs.size());
    std::transform(m_vector.cbegin(), m_vector.cend(), rhs.m_vector.cbegin(),
                   result.begin(),
                   [](const T &a, const T &b) -> T { return a / b; });

    return result;
  }

  Vector operator/(float scalar) const noexcept {
    vector<T> result(m_vector.size());
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
          std::enable_if_t<std::is_arithmetic<T>::value, bool> = true> 
constexpr double dist(const Vector<T> &v) {
  double result = 0.0;
  for(auto const&i: v)
    result += pow(i,2);

  return sqrt(result);
}