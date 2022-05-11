#include "Vector.h"
#include <algorithm>
#include <limits>
#include <math.h>


namespace lin {
auto vector::operator+=(const vector<T, N> &rhs) -> vector<T, N> & {
  std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                 [](const auto &a, const auto &b) { return a + b; });
  return *this;
}

auto vector::operator-=(const vector<T, N> &rhs) -> vector<T, N> & {
  std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                 [](const auto &a, const auto &b) { return a - b; });
  return *this;
}

auto vector::operator/=(const vector<T, N> &rhs) -> vector<T, N> & {
  std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                 [](const auto &a, const auto &b) { return a / b; });
  return *this;
}

auto vector::operator*=(const vector<T, N> &rhs) -> vector<T, N> & {
  std::transform(this->begin(), this->end(), rhs.cbegin(), this->begin(),
                 [](const auto &a, const auto &b) { return a * b; });
  return *this;
}


template <typename T, size_t N>
auto operator<<(std::ostream &os, 
                const vector<T, N> &vec) -> std::ostream & {
  os << "[ ";
  for (auto i : vec)
    os << i << " ";
  return os << "]";
}

template <typename T, size_t N>
auto operator+(const vector<T, N> &lhs, 
            const vector<T, N> &rhs) noexcept -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a + b; });

  return result;
}

template <typename T, size_t N>
auto operator-(const vector<T, N> &lhs, 
            const vector<T, N> &rhs) noexcept -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a - b; });
  return result;
}

template <typename T, size_t N>
auto operator*(const vector<T, N> &lhs, 
            const vector<T, N> &rhs) noexcept -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a * b; });
  return result;
}

template <typename T, size_t N>
auto operator/(const vector<T, N> &lhs, 
            const vector<T, N> &rhs) noexcept -> vector<T, N> {

  vector<T, N> result{};
  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const auto &a, const auto &b) { return a / b; });
  return result;
}


template <typename T, size_t N>
auto operator==(const vector<T, N> &lhs, 
            const vector<T, N> &rhs) noexcept -> bool {
  return std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}

template <typename T, size_t N>
auto operator!=(const vector<T, N> &lhs, 
            const vector<T, N> &rhs) noexcept -> bool {
    return !(lhs == rhs);
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator+(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a + rhs) { return a + rhs; });

  return result;
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator+(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
    return rhs + lhs;
}

template <typename T, size_t N, typename U = float,
          std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator-(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {
  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a - rhs) { return a - rhs; });

  return result;
}

template <typename T, size_t N, typename U = float,
          std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator-(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
    return rhs - lhs;
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator*(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a * rhs) { return a * rhs; });

  return result;
}

template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator*(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
    return rhs * lhs;
}

template <typename T, size_t N, typename U = float,
        std::enable_if_t<std::is_arithmetic_v<U>> >
auto operator/(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N> {

  vector<U, N> result{};
  std::transform(
      lhs.cbegin(), lhs.cend(), result.begin(),
      [&rhs](const auto &a) -> decltype(a / rhs) { return a / rhs; });

  return result;
}

template <typename T, size_t N, typename U = float,
        std::enable_if_t<std::is_arithmetic_v<U>> >
auto operator/(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N> {
    return rhs / lhs;
}

template <typename T, size_t N, 
    typename = std::enable_if_t<std::is_arithmetic<T>::value> >
constexpr auto dist(const vector<T, N> &v) -> double {

  double result{};
  for (auto const &i : v)
    result += pow(i, 2);

  return sqrt(result);
}

template <typename T, size_t N, 
    typename = std::enable_if_t<std::is_arithmetic<T>::value> >
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
    typename = std::enable_if_t<std::is_arithmetic<T>::value> >
constexpr auto sum(const vector<T, N> &v) -> double {

  double result{};
  for (auto const &i : v)
    result += i;

  return result;
}

template <typename T, size_t N, 
    typename = std::enable_if_t<std::is_arithmetic<T>::value> >
constexpr auto prod(const vector<T, N> &v) -> double {

  double result{1};
  for (const auto &i : v)
    result *= i;

  return result;
}

template <typename T, size_t N>
auto dot(const vector<T, N> &lhs, const vector<T, N> &rhs) -> double {
    return sum(rhs * lhs);
}

template <typename T, size_t N>
constexpr auto get_angle(const vector<T, N> &lhs, const vector<T, N> &rhs)
    -> double {
  return acos(dot(lhs, rhs) / (dist(lhs) * dist(rhs))); // implement formatting
}

template <typename T, size_t N, typename U>
constexpr auto normalize(const vector<T, N> &vec) -> vector<U, N> {
  // check for non-zero vector
  return vec / dist(vec);
}

template <typename T, size_t N>
constexpr auto min(const vector<T, N> &rhs) noexcept -> T {
    return *std::min_element(rhs.cbegin(), lhs.cend());
}

template <typename T, size_t N>
constexpr auto max(const vector<T, N> &rhs) noexcept -> T {
  return *std::max_element(rhs.cbegin(), lhs.cend());
}

} // namespace lin