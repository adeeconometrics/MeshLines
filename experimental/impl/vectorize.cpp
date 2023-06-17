#include <algorithm>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <vector>

using std::all_of;
using std::cout;
using std::function;
using std::initializer_list;
using std::vector;

template <typename T>
constexpr auto apply(const initializer_list<vector<T>> &t_list,
                     const function<T(initializer_list<T>)> &fn) -> vector<T> {

  const std::size_t first_size = t_list.begin()->size();
  assert(all_of(std::cbegin(t_list), std::cend(t_list),
                [](const auto a) { return a.size() == first_size; }));

  vector<T> result;
  result.reserve(first_size);

  for (const auto &vec : t_list)
    for (const auto &element : vec)
      result.insert(std::end(result), fn({element}));

  return result;
}

auto main() -> int {
  std::initializer_list<std::vector<int>> args = {
      {1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  std::function<int(std::initializer_list<int>)> increment =
      [](std::initializer_list<int> elements) {
        int result = 0;
        for (const auto &element : elements) {
          result += element + 1;
        }
        return result;
      };

  std::vector<int> result = apply(args, increment);

  // Print the result
  for (const auto &element : result) {
    std::cout << element << " ";
  }
  std::cout << std::endl;

  return 0;
}