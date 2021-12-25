#include <iostream>
#include "Vector.h"

using std::cout;


template <typename T>
void test_arithmetic(const Vector<T>& a, const Vector<T>& b, int c = 1){
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

template <typename T>
void test_comparison(const Vector<T>& a, const Vector<T>& b){
  cout << "testing comparison: \n";
  cout << "\tvector a: " << a << '\n';
  cout << "\tvector b: " << b << '\n';
  cout << "\ta == b: " << std::boolalpha << a == b << '\n';
  cout << "\ta != b: " << std::boolalpha << a != b << '\n';
}

int main() {
  Vector<int> a{1, 2, 3};
  Vector<int> b{1, 1, 1};
  Vector<int> c{1, 1, 1};
  test_arithmetic<int>(a,b);
  test_comparison<int>(b,c);
}