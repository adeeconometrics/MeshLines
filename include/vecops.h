#ifndef __VECOPS_H__
#define __VECOPS_H__

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
// - [x] change thowing function into assert statements

#include <algorithm>
#include <cassert>
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

  assert(lhs.size() == rhs.size());

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

  assert(lhs.size() == rhs.size());

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

  assert(lhs.size() == rhs.size());

  using result_type = common_type_t<T, U>;
  vector<result_type> result{};

  for (size_t i{}; i < lhs.size(); i++) {
    result.emplace_back(static_cast<result_type>(lhs[i]) *
                        static_cast<result_type>(rhs[i]));
  }
  return result;
}

template <typename T, typename U = T>
constexpr auto operator/(const vector<T> &lhs, const vector<U> &rhs)
    -> vector<common_type_t<T, U>> {

  assert(lhs.size() == rhs.size());

  using result_type = common_type_t<T, U>;
  vector<result_type> result{};

  for (size_t i{}; i < lhs.size(); i++) {
    result.emplace_back(static_cast<result_type>(lhs[i]) /
                        static_cast<result_type>(rhs[i]));
  }
  return result;
}

template <typename T>
constexpr auto operator==(const vector<T> &lhs, const vector<T> &rhs) -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.begin());
}

template <typename T>
constexpr auto operator!=(const vector<T> &lhs, const vector<T> &rhs) -> bool {
  return !(lhs == rhs);
}

template <typename T, typename U = T>
constexpr auto operator+(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a + rhs);
  });

  return result;
};

template <typename T, typename U = T>
constexpr auto operator-(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a - rhs);
  });

  return result;
};

template <typename T, typename U = T>
constexpr auto operator*(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a * rhs);
  });

  return result;
};

template <typename T, typename U = T>
constexpr auto operator/(const vector<T> &lhs, U rhs)
    -> vector<common_type_t<T, U>> {

  using result_type = std::common_type_t<T, U>;
  vector<result_type> result;

  std::for_each(lhs.begin(), lhs.cend(), [&result, &rhs](const auto &a) {
    result.emplace_back(a / rhs);
  });

  return result;
};

template <typename T>
constexpr auto operator+=(vector<T> &lhs, const vector<T> &rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto _lhs, const auto _rhs) { return _lhs += _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator-=(vector<T> &lhs, const vector<T> &rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto _lhs, const auto _rhs) { return _lhs -= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator*=(vector<T> &lhs, const vector<T> &rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto _lhs, const auto _rhs) { return _lhs *= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator/=(vector<T> &lhs, const vector<T> &rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), lhs.begin(),
            [](auto _lhs, const auto _rhs) { return _lhs /= _rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator+=(vector<T> &lhs, T rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), lhs.begin(),
            [&rhs](auto _lhs) { return _lhs += rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator-=(vector<T> &lhs, T rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), lhs.begin(),
            [&rhs](auto _lhs) { return _lhs -= rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator*=(vector<T> &lhs, T rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), lhs.begin(),
            [&rhs](auto _lhs) { return _lhs *= rhs; });
  return lhs;
}

template <typename T>
constexpr auto operator/=(vector<T> &lhs, T rhs) -> vector<T> & {
  transform(lhs.cbegin(), lhs.cend(), lhs.begin(),
            [&rhs](auto _lhs) { return _lhs /= rhs; });
  return lhs;
}

template <typename T> constexpr auto max(const vector<T> &v) noexcept -> T {
  return *std::max_element(v.cbegin(), v.cend());
}

template <typename T> constexpr auto min(const vector<T> &v) noexcept -> T {
  return *std::min_element(v.cbegin(), v.cend());
}

} // namespace lin

#endif // __VECOPS_H__
