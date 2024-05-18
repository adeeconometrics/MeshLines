#ifndef __UTILS_H__
#define __UTILS_H__

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

class Timer {
public:
  explicit Timer() : start_time(std::chrono::high_resolution_clock::now()) {}
  explicit Timer(std::size_t t_iterations) : m_iterations(t_iterations) {}
  explicit Timer(std::size_t t_iterations, std::string t_name)
      : m_iterations(t_iterations), m_name(t_name) {}

  ~Timer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    const auto total_duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                             start_time)
            .count();
    const double mean_duration = total_duration / m_iterations;
    const double gflops = (2147483648 / mean_duration);
    std::cout << "mean elapsed time took: " << mean_duration << " ns for "
              << m_name << " or " << gflops << "GFlops" << std::endl;
  }

  auto start() -> void {
    start_time = std::chrono::high_resolution_clock::now();
  }

private:
  std::chrono::high_resolution_clock::time_point start_time;
  std::size_t m_iterations{1};
  std::string m_name{"No Name"};
};

template <typename T, std::size_t N, std::size_t M>
auto rand_array(std::reference_wrapper<std::mt19937> prng)
    -> std::array<T, N * M> {

  std::array<T, N * M> result;
  std::generate_n(std::begin(result), N * M, prng);
  return result;
}

template <typename T, std::size_t N, std::size_t M>
auto rand_vector(std::reference_wrapper<std::mt19937> prng) -> std::vector<T> {

  std::vector<T> result;
  result.reserve(N * M);
  std::generate_n(std::back_inserter(result), N * M, prng);
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