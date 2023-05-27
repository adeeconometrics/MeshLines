#ifndef __UTILS_H__
#define __UTILS_H__

/**
 * @file utils.h
 * @author ddamiana
 * @brief Contains utility functions
 * @version 0.1
 * @date 2023-05-27
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../include/matrix.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

using lin::Matrix;

template <typename T>
auto operator<<(std::ostream &os, const Matrix<T> &matrix) -> std::ostream & {
  static_assert(std::is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");

  if (matrix.empty()) {
    os << "[]" << std::endl;
    return os;
  }

  std::size_t max_width = 0;
  for (const auto &row : matrix) {
    for (const auto &element : row) {
      std::size_t width = std::to_string(element).size();
      if (width > max_width) {
        max_width = width;
      }
    }
  }

  os << "[";
  for (std::size_t i = 0; i < matrix.size(); ++i) {
    if (i != 0) {
      os << " ";
    }
    os << "[";
    for (std::size_t j = 0; j < matrix[i].size(); ++j) {
      os << std::setw(max_width) << matrix[i][j];
      if (j != matrix[i].size() - 1) {
        os << ", ";
      }
    }
    os << "]";
    if (i != matrix.size() - 1) {
      os << '\n';
    }
  }
  return os << "]\n";
}

template <typename T>
auto operator<<(std::ostream &os, const std::vector<T> &v) -> std::ostream & {
  os << "[";
  for (const auto i : v) {
    os << i << " ";
  }
  return os << "]\n";
};

#endif // __UTILS_H__