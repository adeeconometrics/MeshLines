#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

template <typename T, std::size_t Rows, std::size_t Cols, typename = void>
class Matrix {
public:
  Matrix() : m_row(Rows), m_col(Cols) {}

  Matrix(const T &value) : m_row(Rows), m_col(Cols) { m_data.fill(value); }
  Matrix(const std::array<T, Rows * Cols> &t_data)
      : m_data(t_data), m_row(Rows), m_col(Cols) {
    // static_assert()
  }

  auto operator()(std::size_t row, std::size_t col) -> T & {
    return m_data[row * Cols + col];
  }

  auto operator()(std::size_t row, std::size_t col) const -> const T & {
    return m_data[row * Cols + col];
  }

  operator std::array<T, Rows * Cols>() { return m_data; }

  auto data() -> std::array<T, Rows * Cols> & { return m_data; }

  auto data() const -> const std::array<T, Rows * Cols> & { return m_data; }

private:
  std::array<T, Rows * Cols> m_data;
  std::size_t m_row;
  std::size_t m_col;
};

#define STACK_FRAME 256 * 256
template <typename T, std::size_t Rows, std::size_t Cols>
class Matrix<T, Rows, Cols, std::enable_if_t<Rows * Cols >= STACK_FRAME>> {

public:
  Matrix() : m_row(Rows), m_col(Cols) { m_data.reserve(Rows * Cols); }
  constexpr Matrix(std::size_t t_rows, std::size_t t_cols) {
    m_data.resize(t_rows * t_cols, 0);
  }

  Matrix(const std::vector<T> &t_data, std::size_t t_rows, std::size_t t_cols)
      : m_data(t_data), m_row(t_rows), m_col(t_cols) {}

  operator std::vector<T>() const { return m_data; }

  auto operator()(std::size_t row, std::size_t col) -> T & {
    return m_data[row * Cols + col];
  }

  auto operator()(std::size_t row, std::size_t col) const -> const T & {
    return m_data[row * Cols + col];
  }

  auto data() -> std::vector<T> & { return m_data; }

  auto data() const -> const std::vector<T> & { return m_data; }

private:
  std::vector<T> m_data;

  std::size_t m_row;
  std::size_t m_col;
};

// implement transform function to Matrix (in-place transform and
// copy-transform)
template <typename T, std::size_t M, std::size_t N>
auto transpose(const Matrix<T, M, N> &t_matrix) -> Matrix<T, M, N> {
  Matrix<T, M, N> TA = t_matrix;

  for (std::size_t i = 0; i < M; i++)
    for (std::size_t j = 0; j < N; j++)
      TA(j, i) = t_matrix(i, j);

  return TA;
}

template <typename T, std::size_t M, std::size_t N>
auto transpose(const Matrix<T, M, N> &t_matrix) -> void {
  for (std::size_t i = 0; i < M; i++)
    for (std::size_t j = 0; j < N; j++)
      std::swap(t_matrix(j, i), t_matrix(i, j));

  return;
}

template <typename T, std::size_t M, std::size_t N>
constexpr auto operator==(const Matrix<T, M, N> &lhs,
                          const Matrix<T, M, N> &rhs) -> bool {
  return std::equal(std::cbegin(lhs.data()), std::cend(lhs.data()),
                    std::cbegin(rhs.data()));
}

template <typename T, std::size_t M, std::size_t N>
constexpr auto operator!=(const Matrix<T, M, N> &lhs,
                          const Matrix<T, M, N> &rhs) -> bool {
  return !(lhs == rhs);
}

#endif // __MATRIX_H__