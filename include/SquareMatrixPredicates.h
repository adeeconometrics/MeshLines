#ifndef __SQUAREMATRIXPREDICATES_H__
#define __SQUAREMATRIXPREDICATES_H__

/**
 * @file SquareMatrixPredicates.h
 * @author ddamiana
 * @brief contains predicates for testing certain properties of a matrix
 * @version 1.0
 * @date 2022-05-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "SquareMatrix.h"


namespace lin {

template <typename T, size_t N>
constexpr auto is_invertible(const SquareMatrix<T, N> &mat) noexcept -> bool;

template <typename T, size_t N>
constexpr auto is_indefinite(const SquareMatrix<T, N> &mat) noexcept -> bool;

template <typename T, size_t N>
constexpr auto is_echelon(const SquareMatrix<T, N> &mat) noexcept -> bool;

template <typename T, size_t N>
constexpr auto is_symmetric(const SquareMatrix<T, N> &mat) noexcept -> bool;

template <typename T, size_t N>
constexpr auto is_anti_symmetric(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_skew_symmetric(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_triangular(const SquareMatrix<T, N> &mat) noexcept -> bool;

template <typename T, size_t N>
constexpr auto is_upper_triangular(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_lower_triangular(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_diagonal(const SquareMatrix<T, N> &mat) noexcept -> bool;

template <typename T, size_t N>
constexpr auto is_diagonalizable(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_positive_definite(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_positive_semidefinite(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_negative_definite(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_negative_semidefinite(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

template <typename T, size_t N>
constexpr auto is_linearly_dependent(const SquareMatrix<T, N> &mat) noexcept
    -> bool;

} // namespace lin

#endif // __SQUAREMATRIXPREDICATES_H__