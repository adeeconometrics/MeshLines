/**
 * @file metaop.cpp
 * @author ddamiana
 * @brief writing a meta rule for abstract functor
 * @version 0.1
 * @date 2023-05-18
 *
 * @copyright Copyright (c) 2023
 *
 */

// This implementation intends to decouple the operator and vectorized
// operation. This implementation requires a functor and a virtual method thus
// incurring cost. This implementation allows the compiler to write and
// implement operator overloading thus decreasing the surface and maximize
// flexibility by defining a metarule.
// In short, this maximize flexibility but incurs performance cost, try to
// optimize for it.

#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

using std::cout;
using std::transform;
using std::vector;

template <typename T> struct BinOp {
  BinOp() = default;
  BinOp(BinOp<T> &&bin_op) noexcept = default;
  BinOp(const BinOp<T> &bin_op) noexcept = default;

  virtual ~BinOp() = default;
  virtual auto operator()(const T &lhs, const T &rhs) const noexcept -> T = 0;
};

template <typename T> struct Plus : public BinOp<T> {
  auto operator()(const T &lhs, const T &rhs) const noexcept -> T override {
    return lhs + rhs;
  }
};
// solve address boundary error
template <typename T>
auto join(const vector<T> &lhs, const vector<T> &rhs, BinOp<T> &op)
    -> vector<T> {
  vector<T> res{};
  transform(lhs.begin(), lhs.end(), rhs.begin(), res.begin(),
            [&op](auto &_lhs, auto &_rhs) { return op(_lhs, _rhs); });

  return res;
}

template <typename T>
auto operator<<(std::ostream &os, const vector<T> &v) -> std::ostream & {
  for (auto i : v) {
    os << i << ", ";
  }
  return os;
}

auto main() -> int {
  vector<int> v1{1, 2, 3};
  vector<int> v2{4, 5, 6};
  Plus<int> p_int;

  const auto v3 = join(v1, v2, p_int);
  cout << v3;
  return 0;
}