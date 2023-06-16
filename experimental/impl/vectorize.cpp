#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

using std::vector;

template <typename T> using Matrix = vector<vector<T>>;

template <typename T>
auto eval(const std::function<Matrix<T>(const Matrix<T> &)> &func,
          const Matrix<T> &matrix) -> Matrix<T> {
  return func(matrix);
}

template <typename T, typename... Args>
auto eval(const std::function<Matrix<T>(const Matrix<T> &)> &func,
          const Matrix<T> &matrix, Args... args) -> Matrix<T> {
  Matrix<T> result = func(matrix);
  return eval(func, result, args...);
}

// template <typename T, typename U = T>
// constexpr auto operator+(const Matrix<T> &lhs, const Matrix<U> &rhs)
//     -> Matrix<common_type_t<T, U>> {
//   Matrix<common_type_t<T, U>> result(lhs.size(),
//   std::vector<T>(lhs[0].size()));

//   assert(lhs.size() == rhs.size());

//   std::transform(lhs.cbegin(), lhs.cend(), rhs.cbegin(), result.begin(),
//                  [](const vector<T> &a, const vector<U> &b) { return a + b;
//                  });

//   return result;
// }

auto main() -> int {
  Matrix<int> matrix1 = {{1, 2}, {3, 4}};
  Matrix<int> matrix2 = {{5, 6}, {7, 8}};
  Matrix<int> matrix3 = {{9, 10}, {11, 12}};

  std::function<Matrix<int>(const Matrix<int> &)> multiplyByTwo =
      [](const Matrix<int> &matrix) {
        Matrix<int> result = matrix;
        for (auto &row : result) {
          for (auto &element : row) {
            element *= 2;
          }
        }
        return result;
      };

  Matrix<int> result = eval(multiplyByTwo, matrix1, matrix2, matrix3);

  for (const auto &row : result) {
    for (const auto &element : row) {
      std::cout << element << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}