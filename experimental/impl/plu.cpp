#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

using std::cout;
using std::tuple;
using std::vector;

template <typename T> using Matrix = vector<vector<T>>;

template <typename T> auto zero(size_t n) -> Matrix<T> {
  return Matrix<T>{n, vector<T>(n, 0)};
}

template <typename T> auto eye(size_t n) -> Matrix<T> {
  Matrix<T> res = zero<T>(n);

  for (std::size_t i{}; i < n; i++) {
    res[i][i] = 1;
  }
  return res;
}

template <typename T>
auto plu_decomposition(Matrix<T> A) -> tuple<Matrix<T>, Matrix<T>, Matrix<T>> {

  const std::size_t n = A.size();
  Matrix<T> L = zero<T>(n);
  Matrix<T> U = zero<T>(n);
  Matrix<T> P = eye<T>(n);

  for (std::size_t k = 0; k < n; k++) {
    // Find pivot row and swap rows
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++) {
      if (std::abs(A[i][k]) > std::abs(A[p][k])) {
        p = i;
      }
    }
    if (p != k) {
      std::swap(A[k], A[p]);
      std::swap(P[k], P[p]);
    }

    // Perform elimination
    for (std::size_t i = k + 1; i < n; i++) {
      T factor = A[i][k] / A[k][k];
      A[i][k] = factor;
      for (std::size_t j = k + 1; j < n; j++) {
        A[i][j] -= A[k][j] * factor;
      }
    }
  }

  // Extract L and U from A
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++) {
      if (i > j) {
        L[i][j] = A[i][j];
        U[i][j] = 0;
      } else if (i == j) {
        L[i][j] = 1;
        U[i][j] = A[i][j];
      } else {
        L[i][j] = 0;
        U[i][j] = A[i][j];
      }
    }
  }

  return std::make_tuple(P, L, U);
}

template <typename T>
auto operator<<(std::ostream &os, const vector<T> &v) -> std::ostream & {
  for (auto i : v)
    os << i << " ";
  return os << '\n';
}

auto main() -> int {
  const Matrix<double> A = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  const auto [P, L, U] = plu_decomposition(A);
  cout << P;
  cout << L;
  cout << U;

  return 0;
}