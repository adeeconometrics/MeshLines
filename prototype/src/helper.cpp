#include "../include/factorize.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

auto print(const Matrix<double> &matrix) -> void {
  if (matrix.empty()) {
    std::cout << "[]" << std::endl;
    return;
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

  std::cout << "[";
  for (std::size_t i = 0; i < matrix.size(); ++i) {
    if (i != 0) {
      std::cout << " ";
    }
    std::cout << "[";
    for (std::size_t j = 0; j < matrix[i].size(); ++j) {
      std::cout << std::setw(max_width) << matrix[i][j];
      if (j != matrix[i].size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]";
    if (i != matrix.size() - 1) {
      std::cout << '\n';
    }
  }
  std::cout << "]\n";
}