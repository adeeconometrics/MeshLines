#ifndef __HELPER_H__
#define __HELPER_H__

#include "../include/factorize.h"
#include <vector>

auto print(const Matrix<double> &matrix) -> void;

template <typename T> auto print(const std::vector<T> &v) -> void {
  std::cout << "[";
  for (const auto i : v) {
    std::cout << i << " ";
  }
  std::cout << "]\n";
};

#endif // __HELPER_H__