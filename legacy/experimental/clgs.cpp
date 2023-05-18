// classical implementation of gram-shmidt process

#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

using std::cout;
using std::inner_product;
using std::pair;
using std::sqrt;
using std::unique_ptr;
using std::vector;

using Matrix = vector<vector<double>>;

// operator overloads for vector types

auto operator-(const vector<double> lhs, const vector<double> rhs)
    -> vector<double> {
  // assumes lhs.size() == rhs.size()
  vector<double> result(lhs.size());
  for (size_t i{}; i < lhs.size(); ++i) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}

auto operator*(double lhs, const vector<double> rhs) -> vector<double> {
  // assumes lhs.size() == rhs.size()
  vector<double> result(rhs.size());
  for (size_t i{}; i < rhs.size(); ++i) {
    result[i] = rhs[i] * lhs;
  }
  return result;
}

auto operator/(const vector<double> lhs, double rhs) -> vector<double> {
  // assumes lhs.size() == rhs.size()
  vector<double> result(lhs.size());
  for (size_t i{}; i < lhs.size(); ++i) {
    result[i] = lhs[i] / rhs;
  }
  return result;
}

auto operator/(const vector<double> lhs, double rhs) -> vector<double>;

auto eye(size_t n) -> Matrix {
  Matrix I(n, vector<double>(n));
  for (size_t i{}; i < n; ++i)
    for (size_t j{i}; j < i + 1; ++j)
      if (i == j) {
        I[i][j] = 1;
      }
  return I;
}
// SEEM TO CAUSE SEGFAULT!
// auto dot(const vector<double> &a, const vector<double> &b) -> double {
//   // assert(a.size() == b.size());
//   double result{};
//   for (int i{}; i < 3; ++i)
//     result += a[i] * b[i];
//   return result;
// }

auto norm(const std::vector<double> &v) -> double {
  double result{};
  for (const auto i : v)
    result += i;
  return sqrt(result);
}

auto transpose(const Matrix &M) -> Matrix {
  Matrix T(M.size(), vector<double>(M[0].size()));

  for (size_t i{}; i < M.size(); ++i) {
    for (size_t j{}; j < M[i].size(); ++j) {
      T[i][j] = M[j][i];
    }
  }
  return T;
}

auto clgs(const Matrix &A) -> pair<Matrix, Matrix> {
  auto rows = A.size();
  auto R = eye(rows);
  auto AT = transpose(A);

  for (size_t j{}; j < rows; ++j) {
    for (size_t k{}; k < j - 1; ++k) {
      R[k][j] = inner_product(AT[j].begin(), AT[j].end(), AT[k].begin(), 0);
      AT[j] = AT[j] - R[k][j] * AT[k];
    }
    R[j][j] = norm(AT[j]);
    AT[j] = AT[j] / R[j][j];
  }

  auto Q = transpose(AT);

  return std::make_pair(Q, R);
}

auto operator<<(std::ostream &os, const Matrix &M) -> std::ostream & {
  for (auto i : M) {
    for (auto j : i)
      os << j << " ";
    os << '\n';
  }

  return os;
}

auto main() -> int {
  Matrix A = {{{1, 1, 0}, {1, 0, 1}, {0, 1, 1}}};
  auto [Q, R] = clgs(A);

  cout << Q << "\n\n" << R;
}