#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <array>
#include <type_traits>
#include <vector>

template <typename T, std::size_t Rows, std::size_t Cols, typename = void>
class Matrix {
public:
  Matrix() : m_row(Rows), m_col(Cols) {}

  explicit Matrix(const T &value) : m_row(Rows), m_col(Cols) {
    m_data.fill(value);
  }
  explicit Matrix(const std::array<T, Rows * Cols> &t_data)
      : m_data(t_data), m_row(Rows), m_col(Cols) {}

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
class Matrix<T, Rows, Cols, std::enable_if_t<(Rows * Cols) >= (STACK_FRAME)>> {

public:
  Matrix() : m_row(Rows), m_col(Cols) { m_data.reserve(Rows * Cols); }
  constexpr explicit Matrix(std::size_t t_rows, std::size_t t_cols) {
    m_data.resize(t_rows * t_cols, 0);
  }

  explicit Matrix(const std::vector<T> &t_data, std::size_t t_rows,
                  std::size_t t_cols)
      : m_data(t_data), Matrix(t_rows, t_cols) {}

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
#endif // __MATRIX_H__