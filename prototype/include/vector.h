#ifndef __VECTOR_H__
#define __VECTOR_H__

/**
 * @brief the intention of this module is to write a convinient vector
 * implementation that is thinly wrapped in std::vector
 *
 */

// todo:
// - [x] type promotion of compatible integral types
// - [ ] add open mp directive (figure out how to insert it in stl algo)
//      - [ ] make sure that this directive is only added when the compiler
//      supports it
// - [ ] profile performance

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace lin {
using std::is_arithmetic_v, std::common_type_t;
using std::vector;

template <typename T, typename U = T>
constexpr auto operator+(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("rhs size is not same as lhs size");
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameters must be of arithmetic type");

  using result_type = common_type_t<T, U>;
  vector<result_type> result{};

  for (size_t i{}; i < lhs.size(); i++) {
    result.emplace_back(lhs[i] + rhs[i]);
  }
  return result;
}

template <typename T, typename U = T>
constexpr auto operator-(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("rhs size is not same as lhs size");
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameters must be of arithmetic type");

  using result_type = common_type_t<T, U>;
  vector<result_type> result{};

  for (size_t i{}; i < lhs.size(); i++) {
    result.emplace_back(lhs[i] - rhs[i]);
  }
  return result;
}

template <typename T, typename U = T>
constexpr auto operator*(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("rhs size is not same as lhs size");
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameters must be of arithmetic type");

  using result_type = common_type_t<T, U>;
  vector<result_type> result{};

  for (size_t i{}; i < lhs.size(); i++) {
    result.emplace_back(lhs[i] * rhs[i]);
  }
  return result;
}

template <typename T, typename U = T>
constexpr auto operator/(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("rhs size is not same as lhs size");
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameters must be of arithmetic type");

  using result_type = common_type_t<T, U>;
  vector<result_type> result{};

  for (size_t i{}; i < lhs.size(); i++) {
    result.emplace_back(lhs[i] / rhs[i]);
  }
  return result;
}

template <typename T>
constexpr auto operator==(const vector<T> &lhs, const vector<T> &rhs) -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.begin());
}

template <typename T>
constexpr auto operator!=(const vector<T> &lhs, const vector<T> &rhs) -> bool {
  return ~(lhs == rhs);
}

template <typename T, typename U = double>
constexpr auto operator+(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameter must be of arithmetic type");

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a + rhs);
  });

  return result;
};

template <typename T, typename U = double>
constexpr auto operator-(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameter must be of arithmetic type");

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a - rhs);
  });

  return result;
};

template <typename T, typename U = double>
constexpr auto operator*(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameter must be of arithmetic type");

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a * rhs);
  });

  return result;
};

template <typename T, typename U = double>
constexpr auto operator/(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameter must be of arithmetic type");

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a / rhs);
  });

  return result;
};

template <typename T> constexpr auto dist(const vector<T> &v) -> double {
  static_assert(is_arithmetic_v<T>,
                "template parameter must be of arithmetic type");

  double result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto &i) { result += pow(i, 2); });
  return std::sqrt(result);
}

template <typename T>
constexpr auto lp_norm(const vector<T> &v, float p) -> float {
  static_assert(is_arithmetic_v<T>,
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
  static_assert(is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");

  T result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto i) -> void { result += i; });
  return result;
}

template <typename T> constexpr auto prod(std::vector<T> &v) -> T {
  static_assert(is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");

  T result{};
  std::for_each(v.cbegin(), v.cend(),
                [&result](const auto i) -> void { result *= i; });
  return result;
}

template <typename T, typename U = T>
constexpr auto dot(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {
  static_assert(is_arithmetic_v<T> && is_arithmetic_v<U>,
                "template parameters must be of type arithmetic");
  return sum(rhs * lhs);
}

template <typename T, typename U = T>
constexpr auto get_angle(const vector<T> &lhs, const vector<U> &rhs)
    -> common_type_t<T, U> {
  common_type_t<T, U> result;
  return std::acos(dot(lhs, rhs) / (dist(lhs) * dist(rhs)));
}

} // namespace lin

#endif // __VECTOR_H__