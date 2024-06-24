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
#include <cassert>
#include <cmath>
#include <type_traits>

namespace lin {
template <typename T> constexpr auto dist(const vector<T> &v) -> double {

  double result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto &i) { result += pow(i, 2); });
  return std::sqrt(result);
}

template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto ones(size_t t_size) -> vector<T> {
  return vector<T>(1, t_size);
}

template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto zeros(size_t t_size) -> vector<T> {
  return vector<T>(0, t_size);
}

template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto lp_norm(const vector<T> &v, float p) -> float {
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

template <typename T> constexpr auto sum(const std::vector<T> &v) -> T {

  T result{};
  for (const auto &i : v) {
    result += i;
  }
  return result;
}

template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
constexpr auto prod(const std::vector<T> &v) -> T {

  T result{1};
  for (const auto &i : v) {
    result *= i;
  }
  return result;
}

template <typename T>
constexpr auto dot(const vector<T> &lhs, const vector<T> &rhs) -> T {
  return sum(rhs * lhs);
}
/**
 * @brief Get the angle between two vectors. Both vectors must be of the same
 * size.
 *
 * @tparam T
 * @param lhs
 * @param rhs
 * @return T
 */
template <typename T>
constexpr auto get_angle(const vector<T> &lhs, const vector<T> &rhs) -> T {
  assert(lhs.size() == rhs.size() && "Vectors must be of the same size");
  return std::acos(dot(lhs, rhs) / (dist(lhs) * dist(rhs)));
}
/**
 * @brief Returns the normalized vector.
 *
 * @tparam T type of the vector
 * @tparam std::enable_if_t<std::is_arithmetic_v<T>>
 * @param v vector to normalize
 * @return auto
 */
template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto normalize(const vector<T> &v) {
  return v / dist(v);
}
/**
 * @brief Implementation of the cross product. Note: current implementation only
 * works with 3d vectors.
 *
 * @tparam T The type of the vector
 * @tparam U The type of the vector
 * @param lhs The left hand side vector
 * @param rhs The right hand side vector
 * @return vector<common_type_t<T, U>> The cross product of the two vectors
 */
template <typename T, typename U = T>
auto cross(const vector<T> &lhs,
           const vector<U> &rhs) -> vector<common_type_t<T, U>> {
  assert(lhs.size() == 3 && rhs.size() == 3);
  return {lhs[1] * rhs[2] - lhs[2] * rhs[1], lhs[2] * rhs[0] - lhs[0] * rhs[2],
          lhs[0] * rhs[1] - lhs[1] * rhs[0]};
}

} // namespace lin
#endif // __VECFUNC_H__