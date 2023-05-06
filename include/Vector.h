#ifndef VECTOR_H
#define VECTOR_H

/**
 * @file Vector.h
 * @author ddamiana
 * @brief contains methods for handling a mathematical vector
 * @version 1.0
 * @date 2022-05-11
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <array>
#include <iostream>
#include <type_traits>

namespace lin {
using std::array;

template <typename T, size_t N> class vector : public array<T, N> {
  auto operator+=(const vector<T, N> &rhs) -> vector<T, N> &;
  auto operator-=(const vector<T, N> &rhs) -> vector<T, N> &;
  auto operator/=(const vector<T, N> &rhs) -> vector<T, N> &;
  auto operator*=(const vector<T, N> &rhs) -> vector<T, N> &;
};
/**
 * @brief puts vector into output stream
 *
 * @tparam T
 * @tparam N
 * @param os
 * @param vec
 * @return std::ostream&
 */
template <typename T, size_t N>
auto operator<<(std::ostream &os, const vector<T, N> &vec) -> std::ostream &;
/**
 * @brief performs element-wise addition on two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return vector <T,N>
 */
template <typename T, size_t N>
auto operator+(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N>;
/**
 * @brief performs element-wise subtraction on two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return vector<T, N>
 */
template <typename T, size_t N>
auto operator-(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N>;
/**
 * @brief performs element-wise division on two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return vector <T,N>
 */
template <typename T, size_t N>
auto operator/(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N>;
/**
 * @brief performs element-wise multiplication on two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return vector<T, N>
 */
template <typename T, size_t N>
auto operator*(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> vector<T, N>;
/**
 * @brief performs element-wise equality check for two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return true
 * @return false
 */
template <typename T, size_t N>
auto operator==(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> bool;
/**
 * @brief performs element-wise equality check for two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return true
 * @return false
 */
template <typename T, size_t N>
auto operator!=(const vector<T, N> &lhs, const vector<T, N> &rhs) noexcept
    -> bool;
/**
 * @brief performs scalar addition
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U
 * @param lhs
 * @param rhs
 * @return vector <U,N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator+(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N>;
/**
 * @brief performs scalar addition
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator+(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N>;
/**
 * @brief performs scalar subtraction
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator-(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N>;
/**
 * @brief performs scalar subtraction
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator-(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N>;
/**
 * @brief performs scalar division
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator/(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N>;
/**
 * @brief performs scalar division
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator/(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N>;
/**
 * @brief perfoms scalar multiplication
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator*(const vector<T, N> &lhs, U rhs) noexcept -> vector<U, N>;
/**
 * @brief perfoms scalar multiplication
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @tparam U,
 * typename
 * @param lhs
 * @param rhs
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = float,
          typename = std::enable_if_t<std::is_arithmetic_v<U>>>
auto operator*(U lhs, const vector<T, N> &rhs) noexcept -> vector<U, N>;
/**
 * @brief returns the Euclidean distance of a vector
 *
 * @tparam T
 * @tparam N
 * @tparam N,
 * typename
 * @param v
 * @return double
 */
template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto dist(const vector<T, N> &v) -> double;
/**
 * @brief returns the \f$\mathcal{l}_p-\text{norm}$\f of a given vector
 *
 * @tparam T
 * @tparam N
 * @tparam N,
 * typename
 * @param v
 * @param p
 * @return double
 */
template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto lp_norm(const vector<T, N> &v, float p) -> double;
/**
 * @brief returns the cumulative sum of a vector
 *
 * @tparam T
 * @tparam N
 * @tparam N,
 * typename
 * @param v
 * @return double
 */
template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto sum(const vector<T, N> &v) -> double;
/**
 * @brief returns the cumulative product of a vector
 *
 * @tparam T
 * @tparam N
 * @tparam N,
 * typename
 * @param v
 * @return double
 */
template <typename T, size_t N,
          typename = std::enable_if_t<std::is_arithmetic<T>::value>>
constexpr auto prod(const vector<T, N> &v) -> double;
/**
 * @brief returns the dot product of two vector
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return double
 */
template <typename T, size_t N>
auto dot(const vector<T, N> &lhs, const vector<T, N> &rhs) -> double;
/**
 * @brief returns the angle between two vectors
 *
 * @tparam T
 * @tparam N
 * @param lhs
 * @param rhs
 * @return double
 */
template <typename T, size_t N>
constexpr auto get_angle(const vector<T, N> &lhs, const vector<T, N> &rhs)
    -> double;
/**
 * @brief returns a normalized vector
 *
 * @tparam T
 * @tparam N
 * @tparam U
 * @param vec
 * @return vector<U, N>
 */
template <typename T, size_t N, typename U = double>
constexpr auto normalize(const vector<T, N> &vec) -> vector<U, N>;
/**
 * @brief returns the minimum element of the vector
 *
 * @tparam T
 * @tparam N
 * @param rhs
 * @return T
 */
template <typename T, size_t N>
constexpr auto min(const vector<T, N> &rhs) noexcept -> T;
/**
 * @brief returns the maximum value of the vector
 *
 * @tparam T
 * @tparam N
 * @param rhs
 * @return T
 */
template <typename T, size_t N>
constexpr auto max(const vector<T, N> &rhs) noexcept -> T;

} // namespace lin

#endif // VECTOR_H