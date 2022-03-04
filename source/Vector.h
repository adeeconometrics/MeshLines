/**
 * @file Vector.h
 * @author ddamiana
 * @brief contains concepts related to mathematical vector
 * @version 0.3
 * @date 2022-03-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <math.h>
#include <type_traits>

#if __cplusplus >= 201703L
using std::array;
using std::cout;

template <typename T, size_t N> class vector : public array<T, N> {
public:
  auto operator+=(const vector<T, N> &rhs) -> vector<T, N> & {
    std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                   [](const auto &a, const auto &b) { return a + b; });
    return *this;
  }

  auto operator-=(const vector<T, N> &rhs) -> vector<T, N> & {
    std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                   [](const auto &a, const auto &b) { return a - b; });
    return *this;
  }

  auto operator/=(const vector<T, N> &rhs) -> vector<T, N> & {
    std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                   [](const auto &a, const auto &b) { return a / b; });
    return *this;
  }

  auto operator*=(const vector<T, N> &rhs) -> vector<T, N> & {
    std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                   [](const auto &a, const auto &b) { return a * b; });
    return *this;
  }
};

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const vector<T, N> &vec) {
  os << "[ ";
  for (auto i : vec)
    os << i << " ";
  return os << "]";
}

template <typename T, size_t N>
auto operator+(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a + b; });

  return result;
}

template <typename T, size_t N>
auto operator-(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a - b; });
  return result;
}

template <typename T, size_t N>
auto operator*(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a * b; });
  return result;
}

template <typename T, size_t N>
auto operator/(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a / b; });
  return result;
}

template <typename T, size_t N>
auto operator==(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}

template <typename T, size_t N>
auto operator!=(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> bool {
  return !(lhs == rhs);
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator+(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a + rhs) { return a + rhs; });

  return result;
}
// forward declaration for associative operator+
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator+(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
  return rhs + lhs;
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator-(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a - rhs) { return a - rhs; });

  return result;
}
// forward declaration for associative operator-
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator-(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
  return rhs - lhs;
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator*(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a * rhs) { return a * rhs; });

  return result;
}
// forward declaration for associative operator*
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator*(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
  return rhs * lhs;
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator/(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a / rhs) { return a / rhs; });

  return result;
}
// forward declaration for associative operator/
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
auto operator/(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
  return rhs / lhs;
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto dist(const vector<T, N> &v) -> double {

  double result{};
  for (auto const &i : v)
    result += pow(i, 2);

  return sqrt(result);
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto lp_norm(const vector<T, N> &v, float p) -> double {

  if (p == 0.0)
    return std::count_if(v.begin(), v.end(),
                         [](const auto &x) -> bool { return x != 0; });
  if (std::isinf(p)) {
    double max = *std::max_element(v.begin(), v.end(),
                                   [](const auto &a, const auto &b) -> bool {
                                     return std::abs(a) < std::abs(b);
                                   });

    return std::abs(max);
  }

  double result{};
  for (const auto &i : v)
    result += pow(i, p);
  return pow(result, 1 / p);
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto sum(const vector<T, N> &v) -> double {

  double result{};
  for (auto const &i : v)
    result += i;

  return result;
}

template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto prod(const vector<T, N> &v) -> double {

  double result{1};
  for (const auto &i : v)
    result *= i;

  return result;
}

// template <typename T, size_t N>
// constexpr auto is_linearly_independent(const vector<T,N>& vec) -> bool {}

template <typename T, size_t N>
auto dot(const vector<T, N> &lhs, const vector<T, N> &rhs) -> double {
  return sum(rhs * lhs); // does this allocate unnecessary vector?
}

template <typename T, size_t N>
constexpr auto get_angle(const vector<T, N> &lhs, const vector<T, N> &rhs)
    -> double {
  return acos(dot(lhs, rhs) / (dist(lhs) * dist(rhs)));
}

template <typename T, size_t N, typename U = double>
constexpr auto normalize(const vector<T, N> &vec) -> vector<U, N> {
  // check for non-zero vector
  return vec / dist(vec);
}

// test this method either consider this or implement R2 and R3 version
template <typename T, size_t N>
constexpr auto cross(const vector<T, N> &lhs, const vector<T, N> &rhs)
    -> decltype(dist(lhs) * dist(rhs) * sin(get_angle(lhs, rhs)) *
                normalize(rhs)) {

  return dist(lhs) * dist(rhs) * sin(get_angle(lhs, rhs)) * normalize(rhs);
}

#endif