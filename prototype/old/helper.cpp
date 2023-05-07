#include "../include/factorize.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

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