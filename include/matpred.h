#ifndef __MATPRED_H__
#define __MATPRED_H__

#include "../include/matrix.h"

namespace lin {

template <typename T>
constexpr auto is_tril(const Matrix<T> &M) noexcept -> bool {
  const std::size_t rows = M.size();
  const std::size_t cols = M[0].size();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i + 1; j < cols; j++) {
      if (j < rows && i < cols) {
        if (M[i][j] != 0)
          return false;
      }
    }
  }
  return true;
}

template <typename T>
constexpr auto is_triu(const Matrix<T> &M) noexcept -> bool {
  const std::size_t rows = M.size();
  const std::size_t cols = M[0].size();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i + 1; j < cols; j++) {
      if (j < rows && i < cols) {
        if (M[j][i] != 0)
          return false;
      }
    }
  }
  return true;
}

template <typename T>
constexpr auto is_diag(const Matrix<T> &M) noexcept -> bool {
  const std::size_t rows = M.size();
  const std::size_t cols = M[0].size();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i + 1; j < cols; j++) {
      if (j < rows && i < cols) {
        if (M[j][i] != 0 || M[i][j] != 0)
          return false;
      }
    }
  }
  return true;
}

} // namespace lin

#endif // __MATPRED_H__