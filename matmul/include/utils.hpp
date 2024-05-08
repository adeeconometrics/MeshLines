#ifndef __UTILS_H__
#define __UTILS_H__

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

class Timer {
public:
  explicit Timer()
      : start_time(std::chrono::high_resolution_clock::now()), m_iterations(1) {
  }
  explicit Timer(std::size_t t_iterations) : m_iterations(t_iterations) {}

  ~Timer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                              end_time - start_time)
                              .count();
    auto mean_duration = total_duration / m_iterations;
    std::cout << "Mean elapsed time: " << mean_duration << " ns" << std::endl;
  }

  auto start() -> void {
    start_time = std::chrono::high_resolution_clock::now();
  }

private:
  std::chrono::high_resolution_clock::time_point start_time;
  std::size_t m_iterations{};
};

template <typename T>
auto make_matrix(const std::size_t row, const std::size_t col,
                 std::reference_wrapper<std::mt19937> prng)
    -> std::vector<std::vector<T>> {

  std::vector<std::vector<T>> result;
  result.reserve(row);

  for (std::size_t i = 0; i < row; i++) {
    std::vector<T> row;
    row.reserve(col);
    std::generate_n(std::back_inserter(row), col, prng);
    result.emplace_back(row);
  }

  return result;
}

// template <typename T, std::size_t N, std::size_t M>
// auto make_matrix(std::reference_wrapper<std::mt19937> prng)
//     -> std::array<std::array<T, N>, M> {

//   std::array<std::array<T, N>, M> result;
//   std::generate_n(result.data(), N * M, prng);

//   return result;
// }

template <typename T, std::size_t N, std::size_t M>
auto make_matrix(std::reference_wrapper<std::mt19937> prng)
    -> std::array<std::array<T, N>, M> {

  std::array<std::array<T, N>, M> result;
  for (size_t i = 0; i < M; i++) {
    std::array<T, N> _temp;
    std::generate_n(_temp.begin(), N, prng);
    result[i] = _temp;
  }

  return result;
}

template <typename T, std::size_t N, std::size_t M>
auto make_flat(std::reference_wrapper<std::mt19937> prng)
    -> std::array<T, N * M> {

  std::array<T, N * M> result;
  std::generate_n(result.begin(), N * M, prng);
  return result;
}

template <typename T, std::size_t N, std::size_t M>
auto make_vmatrix(std::reference_wrapper<std::mt19937> prng) -> std::vector<T> {

  std::vector<T> result;
  result.reserve(N * M);

  std::generate_n(std::back_inserter(result), N * M, prng);
  return result;
}

template <typename T, std::size_t N>
auto make_vector(std::reference_wrapper<std::mt19937> prng)
    -> std::array<T, N> {
  std::array<T, N> result;
  std::generate_n(std::begin(result), N, prng);
  return result;
}

template <typename T>
auto operator<<(std::ostream &os, const std::vector<T> &v) -> std::ostream & {
  for (auto i : v)
    os << i << " ";
  return os << '\n';
}

template <typename T, std::size_t N>
auto operator<<(std::ostream &os, const std::array<T, N> &v) -> std::ostream & {
  for (auto i : v)
    os << i << " ";
  return os << '\n';
}

#endif // __UTILS_H__