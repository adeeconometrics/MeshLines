#ifndef __VECFUNC_H__
#define __VECFUNC_H__

/**
 * @file vecfunc.h
 * @author ddamiana
 * @brief Contains linalg functions for processing vectors
 * @version 0.1
 * @date 2023-05-27
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../include/vecops.h"

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace lin {
template <typename T> constexpr auto dist(const vector<T> &v) -> double {

  double result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto &i) { result += pow(i, 2); });
  return std::sqrt(result);
}

template <typename T> constexpr auto ones(size_t t_size) -> vector<T> {
  assert(std::is_arithmetic_v<T>);
  return vector<T>(1, t_size);
}

template <typename T> constexpr auto zeros(size_t t_size) -> vector<T> {
  assert(std::is_arithmetic_v<T>);
  return vector<T>(0, t_size);
}

template <typename T>
constexpr auto lp_norm(const vector<T> &v, float p) -> float {
  static_assert(std::is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");
  if (p == 0.0) {
    return std::count_if(v.cbegin(), v.cend(),
                         [](const auto &x) -> bool { return x != 0; });
  }

  if (std::isinf(p)) {
    const double max = *std::max_element(
        v.cbegin(), v.cend(), [](const auto &a, const auto &b) -> bool {
          return std::abs(a) < std::abs(b);
        });

    return std::abs(max);
  }

  double result{};
  std::for_each(v.cbegin(), v.cend(), [&result, &p](const auto i) -> void {
    result += std::pow(i, p);
  });

  return std::pow(result, 1 / p);
};

template <typename T> constexpr auto sum(std::vector<T> &v) -> T {
  static_assert(std::is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");

  T result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto i) -> void { result += i; });
  return result;
}

template <typename T> constexpr auto prod(std::vector<T> &v) -> T {
  static_assert(std::is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");

  T result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto i) -> void { result *= i; });
  return result;
}

template <typename T, typename U = T>
constexpr auto dot(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {
  static_assert(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>,
                "template parameters must be of type arithmetic");
  return sum(rhs * lhs);
}

template <typename T, typename U = T>
constexpr auto get_angle(const vector<T> &lhs, const vector<U> &rhs)
    -> common_type_t<T, U> {
  common_type_t<T, U> result;
  return std::acos(dot(lhs, rhs) / (dist(lhs) * dist(rhs)));
}

template <typename T>
constexpr auto normalize(const vector<T> &v) -> vector<T> {
  return v / dist(v); // check for condition of [int]/[double] and [int]/[int]
}
} // namespace lin
#endif // __VECFUNC_H__