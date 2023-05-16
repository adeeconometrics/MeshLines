#include <iostream>
#include "vector.h"

using std::cout;


template <typename T, size_t N>
void test_arithmetic(const vector<T,N>& a, const vector<T,N>& b, int c = 1){
  cout << "testing arithmetic operators: \n";
  cout << "\tvector a: " << a << '\n';
  cout << "\tvector b: " << b << '\n';
  cout << "\tconstant c: " << c << '\n';

  cout << "\ta + b: " << a + b << '\n';
  cout << "\ta + c: " << a + c << '\n';

  cout << "\ta - b: " << a - b << '\n';
  cout << "\ta - c: " << a - c << '\n';

  cout << "\ta * b: " << a * b << '\n';
  cout << "\ta * c: " << a * c << '\n';

  cout << "\ta / b: " << a / b << '\n';
  cout << "\ta / c: " << a / c << '\n';
}

template <typename T, size_t N>
void test_comparison(const vector<T,N> &a, const vector<T,N> &b) {
  bool flag{};
  cout << "testing comparison: \n";
  cout << "\tvector a: " << a << '\n';
  cout << "\tvector b: " << b << '\n';
  flag = a == b;
  cout << "\ta == b: " << std::boolalpha << flag << '\n';
  flag = a != b;
  cout << "\ta != b: " << std::boolalpha << flag << '\n';
}

int main() {
  const vector<int,3> a{1, 2, 3};
  const vector<int,3> b{1, 1, 1};
  const vector<int,3> c{1, 1, 1};
  test_arithmetic(a,b);
  test_comparison(b,c);
}