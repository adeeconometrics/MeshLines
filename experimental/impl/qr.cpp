#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

using std::cout;
using std::tuple;
using std::vector;

#include <iomanip>
#include <string>
#include <type_traits>

template <typename T> using Matrix = vector<vector<T>>;

template <typename T,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto operator<<(std::ostream &os, const Matrix<T> &matrix) -> std::ostream & {

  if (matrix.empty()) {
    os << "[]" << std::endl;
    return os;
  }

  std::size_t max_width = 0;
  for (const auto &row : matrix) {
    for (const auto &element : row) {
      std::size_t width = std::to_string(element).size();
      if (width > max_width) {
        max_width = width;
      }
    }
  }

  os << "[";
  for (std::size_t i = 0; i < matrix.size(); ++i) {
    if (i != 0) {
      os << " ";
    }
    os << "[";
    for (std::size_t j = 0; j < matrix[i].size(); ++j) {
      os << std::setw(max_width) << matrix[i][j];
      if (j != matrix[i].size() - 1) {
        os << ", ";
      }
    }
    os << "]";
    if (i != matrix.size() - 1) {
      os << '\n';
    }
  }
  return os << "]\n";
}

template <typename T> auto mask_triu(Matrix<T> &A) -> void {
  const std::size_t rows = A.size();
  const std::size_t cols = A[0].size();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = i + 1; j < cols; j++) {
      if (j < rows && i < cols) {
        A[j][i] = 0;
      }
    }
  }
}

template <typename T>
auto qr_householder(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  const std::size_t m = A.size();
  const std::size_t n = A[0].size();

  // Initialize Q and R
  Matrix<T> Q(m, std::vector<T>(m));
  Matrix<T> R(A);

  // Compute Householder reflections and apply them to R
  for (std::size_t k = 0; k < n; k++) {
    std::vector<T> x(m - k);
    for (std::size_t i = k; i < m; i++) {
      x[i - k] = R[i][k];
    }

    T norm_x =
        std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));

    std::vector<T> v(m - k);
    v[0] = x[0] < 0 ? x[0] - norm_x : x[0] + norm_x;
    for (std::size_t i = 1; i < m - k; i++) {
      v[i] = x[i];
    }

    T norm_v =
        std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));

    Matrix<T> H(m - k, std::vector<T>(m - k));
    for (std::size_t i = 0; i < m - k; i++) {
      for (std::size_t j = 0; j < m - k; j++) {
        if (i == j) {
          H[i][j] = 1 - 2 * v[i] * v[j] / (norm_v * norm_v);
        } else {
          H[i][j] = -2 * v[i] * v[j] / (norm_v * norm_v);
        }
      }
    }

    // Update R with H
    for (std::size_t i = k; i < m; i++) {
      for (std::size_t j = k; j < n; j++) {
        T sum = 0;
        for (std::size_t h = 0; h < m - k; h++) {
          sum += H[i - k][h] * R[h + k][j];
        }
        R[i][j] = sum;
      }
    }

    // Update Q with H
    for (std::size_t i = 0; i < m; i++) {
      for (std::size_t j = k; j < m; j++) {
        T sum = 0;
        for (std::size_t h = 0; h < m - k; h++) {
          sum += Q[i][h + k] * H[h][j - k];
        }
        Q[i][j] = sum;
      }
    }
  }

  return std::make_tuple(Q, R);
}

template <typename T>
auto qr_gm(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>> {
  const std::size_t _size = A.size();

  Matrix<T> Q(_size, vector<T>(_size, 0));
  Matrix<T> R = A;

  // Compute Q and R using the Gram-Schmidt process
  for (std::size_t j = 0; j < _size; ++j) {
    // Compute the jth column of Q
    for (std::size_t i = 0; i < _size; ++i) {
      Q[i][j] = R[i][j];
    }
    for (std::size_t k = 0; k < j; ++k) {
      T dot_product = 0;
      for (std::size_t i = 0; i < _size; ++i) {
        dot_product += Q[i][k] * R[i][j];
      }
      for (std::size_t i = 0; i < _size; ++i) {
        Q[i][j] -= dot_product * Q[i][k];
      }
    }
    // Normalize the jth column of Q
    T norm = 0;
    for (std::size_t i = 0; i < _size; ++i) {
      norm += Q[i][j] * Q[i][j];
    }
    norm = std::sqrt(norm);
    for (std::size_t i = 0; i < _size; ++i) {
      Q[i][j] /= norm;
    }
    // Compute the jth row of R
    for (std::size_t i = j; i < _size; ++i) {
      R[j][i] = 0;
      for (std::size_t k = 0; k < _size; ++k) {
        R[j][i] += Q[k][j] * A[k][i];
      }
    }
  }

  mask_triu(R);
  return std::make_tuple(Q, R);
}

auto main() -> int {
  using std::sqrt;

  Matrix<double> A{{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  const Matrix<double> Q{{2. / sqrt(29), 5. / sqrt(29), 0},
                         {3. / sqrt(29), (-6. * sqrt(29)) / 145., 4. / 5.},
                         {4. / sqrt(29), (-8. * sqrt(29)) / (145.), -3. / 5.}};

  const Matrix<double> R{{2. * sqrt(29), 31. / sqrt(29), 9. / sqrt(29)},
                         {0, 5. / sqrt(29), 11. / (5. * sqrt(29))},
                         {0, 0, 1. / 5}};

  // const std::size_t rows = _size;
  // const std::size_t cols = A[0].size();

  const auto QRgm = qr_gm(A);
  const auto QRHouseholder = qr_householder(A);

  cout << Q << "\n\n";
  cout << std::get<0>(QRgm) << "\n\n";
  cout << std::get<0>(QRHouseholder) << "\n\n\n";
  cout << "-------------------------------------------------\n";
  cout << R << "\n\n";
  cout << std::get<1>(QRgm) << "\n\n";
  cout << std::get<1>(QRHouseholder) << "\n\n";

  return 0;
};