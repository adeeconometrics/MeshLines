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

#ifdef DEBUG
template <typename T, std::size_t Rows, std::size_t Cols,
          typename = typename std::enable_if_t<std::is_arithmetic_v<T>>>
auto operator<<(std::ostream &os,
                const lin::Matrix<T, Rows, Cols> &matrix) -> std::ostream & {

  if (matrix.empty()) {
    os << "[]" << std::endl;
    return os;
  }
  for (std::size_t i{}; i < Rows; i++) {
    os << "[";
    for (std::size_t j{}; j < Cols; j++) {
      os << std::setw(5) << matrix(i, j) << " ";
    }
    os << "]" << std::endl;
  }
  return os;
}

template <typename T>
auto operator<<(std::ostream &os, const std::vector<T> &v) -> std::ostream & {
  os << "[";
  for (const auto i : v) {
    os << i << " ";
  }
  return os << "]\n";
};
#endif

template <typename Fn, typename ArgType>
constexpr auto apply_fn(Fn &&fn, const lin::vector<ArgType> &v)
    -> lin::vector<std::invoke_result_t<Fn &, const ArgType &>> {

  using ResultType = std::invoke_result_t<Fn &, const ArgType &>;
  lin::vector<ResultType> result;
  result.reserve(v.size());

  std::transform(std::cbegin(v), std::cend(v), std::back_inserter(result), fn);

  return result;
}

// template <typename Fn, typename ArgType>
// constexpr auto apply_fn(Fn &&functor, const lin::Matrix<ArgType> &v)
//     -> lin::Matrix<std::invoke_result_t<Fn &, const ArgType &>> {

//   using ResultType = std::invoke_result_t<Fn &, const ArgType &>;
//   lin::Matrix<ResultType> result{};

//   for (const auto row : v) {
//     result.emplace_back(apply_fn(row, functor));
//   }

//   return result;
// }

// template <typename Collection>
// constexpr auto flatten(const Collection &col)
//     -> lin::vector<typename Collection::value_type> {

//   lin::vector<typename Collection::value_type> flattened;

//   return flattened;
// }

// rows: (Mat) -> Mat
// cols: (Mat) -> Mat

#endif // __UTILS_H__