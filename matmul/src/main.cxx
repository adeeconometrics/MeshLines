#include "../include/matmul.hpp"
#include "../include/matrix.hpp"
#include "../include/utils.hpp"
#include <cassert>
#include <iostream>
#include <random>

void test_matmul() {
  std::mt19937 prng(42);
  auto lhs_matrix = rand_matrix<int, 256, 256>(prng);
  auto rhs_matrix = rand_matrix<int, 256, 256>(prng);
  auto result = Matrix<int, 256, 256>();

  result = iterative(lhs_matrix, rhs_matrix);

  for (std::size_t i = 0; i < 256; ++i) {
    for (std::size_t j = 0; j < 256; ++j) {
      int sum = 0;
      for (std::size_t k = 0; k < 256; ++k) {
        sum += lhs_matrix(i, k) * rhs_matrix(k, j);
      }
      assert(result(i, j) == sum);
    }
  }
  std::cout << "Test passed" << std::endl;
}

auto main() -> int {
  test_matmul();
  return 0;
}