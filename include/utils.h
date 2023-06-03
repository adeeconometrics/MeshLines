#ifndef __UTILS_H__
#define __UTILS_H__

/**
 * @file utils.h
 * @author ddamiana
 * @brief Contains utility functions
 * @version 0.1
 * @date 2023-05-27
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../include/matrix.h"
#include "../include/vecops.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

template <typename T>
auto operator<<(std::ostream &os, const lin::Matrix<T> &matrix)
    -> std::ostream & {
  static_assert(std::is_arithmetic_v<T>,
                "template parameter must be of type arithmetic");

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

template <typename T>
auto operator<<(std::ostream &os, const std::vector<T> &v) -> std::ostream & {
  os << "[";
  for (const auto i : v) {
    os << i << " ";
  }
  return os << "]\n";
};

template <typename Fn, typename ArgType>
constexpr auto apply_fn(Fn &&fn, const lin::vector<ArgType> &v)
    -> lin::vector<std::invoke_result_t<Fn &, const ArgType &>> {

  using ResultType = std::invoke_result_t<Fn &, const ArgType &>;
  lin::vector<ResultType> result;
  result.reserve(v.size());

  std::transform(std::cbegin(v), std::cend(v), std::back_inserter(result), fn);

  return result;
}

template <typename Fn, typename ArgType>
constexpr auto apply_fn(Fn &&functor, const lin::Matrix<ArgType> &v)
    -> lin::Matrix<std::invoke_result_t<Fn &, const ArgType &>> {

  using ResultType = std::invoke_result_t<Fn &, const ArgType &>;
  lin::Matrix<ResultType> result{};

  for (const auto row : v) {
    result.emplace_back(apply_fn(row, functor));
  }

  return result;
}

// template <typename Collection>
// constexpr auto flatten(const Collection &col)
//     -> lin::vector<typename Collection::value_type> {

//   lin::vector<typename Collection::value_type> flattened;

//   return flattened;
// }

// rows: (Mat) -> Mat
// cols: (Mat) -> Mat

#endif // __UTILS_H__