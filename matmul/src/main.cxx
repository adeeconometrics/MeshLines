#include "../include/matmul.hpp"
#include "../include/matrix.hpp"
#include "../include/utils.hpp"
#include <cassert>
#include <functional>
#include <iostream>
#include <random>
#include <string>

template <typename T, std::size_t Rows, std::size_t Cols>
auto bench(
    const std::function<Matrix<T, Rows, Cols>(
        const Matrix<T, Rows, Cols> &, const Matrix<T, Rows, Cols> &)> &t_func,
    const Matrix<T, Rows, Cols> &t_lhs, const Matrix<T, Rows, Cols> &t_rhs,
    std::string t_name, std::size_t t_iter = 2) -> Matrix<T, Rows, Cols> {

  {
    Timer timer{t_iter, t_name};
    timer.start();
    for (std::size_t i = 0; i < t_iter; ++i)
      t_func(t_lhs, t_rhs);
  }
  return t_func(t_lhs, t_rhs);
}

auto test_matmul() -> void {

  std::mt19937 prng(42);
  constexpr std::size_t Rows = 255;
  constexpr std::size_t Cols = 255;

  Matrix<float, Rows, Cols> lhs_matrix{rand_array<float, Rows, Cols>(prng)};
  Matrix<float, Rows, Cols> rhs_matrix{rand_array<float, Rows, Cols>(prng)};

  auto iter_func = std::function<Matrix<float, Rows, Cols>(
      const Matrix<float, Rows, Cols> &, const Matrix<float, Rows, Cols> &)>(
      loop_reorder<float, Rows, Cols>);
  auto iter_mat = bench(iter_func, lhs_matrix, rhs_matrix, "iterative");
  // auto loop_reorder_func = std::function<Matrix<float, Rows, Cols>(
  //     const Matrix<float, Rows, Cols> &, const Matrix<float, Rows, Cols>
  //     &)>( loop_reorder<float, Rows, Cols>);
  // auto loop_reorder_mat =
  //     bench(loop_reorder_func, lhs_matrix, rhs_matrix, "loop_reorder");

  // auto blocked = block_multiply<float, Rows, Cols>(lhs_matrix, rhs_matrix);
  auto blocked_func = std::function<Matrix<float, Rows, Cols>(
      const Matrix<float, Rows, Cols> &, const Matrix<float, Rows, Cols> &)>(
      gemm<float, Rows, Cols>);
  auto blocked = bench(blocked_func, lhs_matrix, rhs_matrix, "blocked");

  auto threaded_gemm_func = std::function<Matrix<float, Rows, Cols>(
      const Matrix<float, Rows, Cols> &, const Matrix<float, Rows, Cols> &)>(
      threaded_gemm<float, Rows, Cols>);
  auto threaded_gemm =
      bench(threaded_gemm_func, lhs_matrix, rhs_matrix, "threaded_gemm");

  auto async_gemm_func = std::function<Matrix<float, Rows, Cols>(
      const Matrix<float, Rows, Cols> &, const Matrix<float, Rows, Cols> &)>(
      async_gemm<float, Rows, Cols>);
  auto async_gemm =
      bench(async_gemm_func, lhs_matrix, rhs_matrix, "async_gemm", 2);
  assert(iter_mat == blocked);

  // for (std::size_t i = 0; i < Rows; ++i) {
  //   for (std::size_t j = 0; j < Cols; ++j) {
  //     float sum = 0;
  //     for (std::size_t k = 0; k < Rows; ++k) {
  //       sum += lhs_matrix(i, k) * rhs_matrix(k, j);
  //     }
  //     assert(result(i, j) == sum);
  //   }
  // }
  // std::cout << "Test passed" << std::endl;
}

auto main() -> int {
  test_matmul();
  return 0;
}