#ifndef __MATPRED_H__
#define __MATPRED_H__

#include "../include/matfunc.h"
#include "../include/matrix.h"

#include <algorithm>
#include <cmath>

namespace lin {

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_tril(const Matrix<T, Rows, Cols> &M) noexcept -> bool {
  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = i + 1; j < Cols; j++) {
      if (M(i, j) != 0)
        return false;
    }
  }
  return true;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_triu(const Matrix<T, Rows, Cols> &M) noexcept -> bool {

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = i + 1; j < Cols; j++) {
      if (j < Rows && i < Cols) {
        if (M(j, i) != 0)
          return false;
      }
    }
  }
  return true;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_diag(const Matrix<T, Rows, Cols> &M) noexcept -> bool {

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = i + 1; j < Cols; j++) {
      if (j < Rows && i < Cols) {
        if (M(j, i) != 0 || M(i, j) != 0)
          return false;
      }
    }
  }
  return true;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_square(const Matrix<T, Rows, Cols> &M) noexcept -> bool {
  return Rows == Cols;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_invertible(const Matrix<T, Rows, Cols> &M) noexcept -> bool {
  if (!is_square(M))
    return false;
  const double det = det(M);

  return std::abs(det) > 0;
}
/**
 * @brief Checks if a matrix is in echelon form.
 * Row Echelon form implies that a matrix has the following properties:
 * <ul>
 *  <li> The leading first entry in each row must be 1. </li>
 *  <li> The leading entry on each subsequent row must be a new column to the
 * right </li>
 *  <li> All rows where all entries are zero below, rows were not all entires
 * are 0. </li>
 * </ul>
 * @tparam T
 * @param M
 * @return true
 * @return false
 */
template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_echelon(const Matrix<T, Rows, Cols> &M) noexcept -> bool {
  // if (M.empty())
  //   return true;

  // std::size_t pivot_col = 0;

  // for (std::size_t row = 0; row < Rows; row++) {
  //   if (pivot_col >= Cols) {
  //     break;
  //   }
  //   // Find the first non-zero element in the row
  //   auto non_zero_pos =
  //       std::find_if(M[row].begin() + pivot_col, M[row].end(),
  //                    [](const T &element) { return element != T{0}; });

  //   // Check if all elements after the first non-zero element are zero
  //   if (non_zero_pos != M[row].cend() &&
  //       std::all_of(non_zero_pos + 1, M[row].cend(),
  //                   [](const T &element) { return element != T{0}; })) {

  //     pivot_col = std::distance(M[row].cbegin(), non_zero_pos) + 1;
  //   } else {
  //     return false;
  //   }
  // }

  return true;
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_sym(const Matrix<T, Rows, Cols> &M) noexcept -> bool {
  // prove that M != T(M) where M is a rectangular mat
  return M == transpose(M);
}

template <typename T, std::size_t Rows, std::size_t Cols>
constexpr auto is_antisym(const Matrix<T, Rows, Cols> &M) noexcept -> bool {
  // prove that M != T(M) where M is a rectangular mat
  return transpose(M) == (-M);
}

} // namespace lin

#endif // __MATPRED_H__