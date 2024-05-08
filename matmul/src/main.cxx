#include "../include/matrix.hpp"
#include "../include/utils.hpp"
#include <iostream>
#include <random>

auto main() -> int {
  std::mt19937 rng_a(67);
  auto a = rand_array<float, 3, 3>(std::ref(rng_a));

  std::cout << a;
  return 0;
}