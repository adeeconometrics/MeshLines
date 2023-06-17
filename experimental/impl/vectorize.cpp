#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

using std::vector;

template <typename T> using Matrix = vector<vector<T>>;

template <typename T, typename... Args>
auto eval(const std::function<Matrix<T>(const Matrix<T> &)> &func,
          const Matrix<T> &matrix, Args &&... args) -> Matrix<T> {
  return (func(matrix), ..., func(std::forward<Args>(args)));
}

template <typename T>
constexpr auto apply(const std::initializer_list<std::vector<T>> &t_list,
                     const std::function<T(std::initializer_list<T>)> &fn)
    -> vector<T> {

  const std::size_t first_size = t_list.begin()->size();
  assert(std::all_of(t_list.cbegin(), t_list.cend(),
                     [](auto v) -> bool { return v.size() == first_size; }));

  std::vector<T> result;
  result.reserve(first_size);

  for (std::size_t i = 0; i < first_size; i++) {
    // result[i] = fn(~elements ~);
    result[i] = fn({t_list});
  }

  return result;
}

template <typename T>
constexpr auto operator+(const vector<T> &lhs, const vector<T> &rhs)
    -> vector<T> {
  assert(lhs.size() == rhs.size());
  vector<T> result;
  result.reserve(lhs.size());

  std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs),
                 std::back_inserter(result), std::plus<T>());
  return result;
}

template <typename T>
constexpr auto operator+(const Matrix<T> &lhs, const Matrix<T> &rhs)
    -> Matrix<T> {
  Matrix<T> result(lhs.size(), std::vector<T>(lhs[0].size()));

  assert(lhs.size() == rhs.size());

  std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
                 [](const vector<T> &a, const vector<T> &b) { return a + b; });

  return result;
}

auto main() -> int {
  Matrix<int> matrix1 = {{1, 2}, {3, 4}};
  Matrix<int> matrix2 = {{5, 6}, {7, 8}};
  Matrix<int> matrix3 = {{9, 10}, {11, 12}};

  std::function<Matrix<int>(const Matrix<int> &)> addFunc =
      [](const Matrix<int> &m) { return m; };

  Matrix<int> result = eval(addFunc, matrix1, matrix2, matrix3);
  // Evaluates as: addFunc(matrix1, addFunc(matrix2, addFunc(matrix3)))

  for (const auto &row : result) {
    for (const auto &element : row) {
      std::cout << element << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}