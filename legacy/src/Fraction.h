/**
 * @file Fraction.h
 * @author ddamiana
 * @brief representation of a Fraction in math
 * @version 0.1
 * @date 2021-12-24
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <iostream>
#include <type_traits>

template <typename T,
          typename = std::enable_if_t<std::is_integral<T>::value>> 
class Fraction {
private:
  T m_n{}, m_d{}, m_gcd{};

public:
  Fraction(T n, T d = 1) {
    static_assert(std::is_integral<T>::value, "Integral value is required.");
    m_gcd = gcd(n, d);
    m_n = n / m_gcd;
    m_d = d / m_gcd;
  }

  constexpr auto as_ratio() const -> double { return m_n / m_d; }
  constexpr auto get_num() const -> T { return m_n; }
  constexpr auto get_den() const -> T { return m_d; }
  constexpr auto get_gcd() const -> T { return m_gcd; }

private:
  constexpr T gcd(T a, T b) const noexcept {
    return b == 0 ? a : gcd(b, a % b);
  }
};

template <typename T>
std::ostream &operator<<(std::ostream &ss, const Fraction<T> &F) {
  return ss << F.get_num() << '/' << F.get_den();
}

template <typename T>
Fraction<T> operator+(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return {lhs.get_num() * rhs.get_den() + rhs.get_num() * lhs.get_den(),
          lhs.get_den() * rhs.get_den()};
}

template <typename T>
Fraction<T> operator-(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return {lhs.get_num() * rhs.get_den() - rhs.get_num() * lhs.get_den(),
          lhs.get_den() * rhs.get_den()};
}

template <typename T>
Fraction<T> operator*(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return {lhs.get_num() * rhs.get_num(), lhs.get_den() * rhs.get_den()};
}

template <typename T>
Fraction<T> operator/(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return {lhs.get_num() * rhs.get_den(), lhs.get_den() * rhs.get_num()};
}

template <typename T>
constexpr bool operator==(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return lhs.get_num() == rhs.get_num() && lhs.get_den() == rhs.get_den();
}

template <typename T>
constexpr bool operator!=(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return !(lhs == rhs);
}

template <typename T>
constexpr bool operator<(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return lhs.as_ratio() < rhs.as_ratio();
}

template <typename T>
constexpr bool operator>(const Fraction<T> &lhs, const Fraction<T> &rhs) {
  return lhs.as_ratio() > rhs.as_ratio();
}