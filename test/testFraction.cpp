#include "Fraction.h"
#include <iostream>

using std::cout;

void test_arithmetic(){
    Fraction<int>f0 (1,2);
    Fraction<int>f1 (3,4);

    cout << "\ntesting arithmetic:\n";
    cout << '\t' << f0 << " + " << f1 << " = " << f0 + f1 << '\n';
    cout << '\t' << f0 << " - " << f1 << " = " << f0 - f1 << '\n';
    cout << '\t' << f0 << " * " << f1 << " = " << f0 * f1 << '\n';
    cout << '\t' << f0 << " / " << f1 << " = " << f0 / f1 << '\n';
}

void test_comparison() {
  Fraction<int> F0(2, 4);
  Fraction<int> F1(8, 16);

  cout << "\ntesting comparison:\n";

  bool flag = F0 == F1;
  cout << '\t' << F0 << " == " << F1 << " => " << std::boolalpha << flag << '\n';
  flag = F0 != F1;
  cout << '\t' << F0 << " != " << F1 << " => " << std::boolalpha << flag << '\n';
  flag = F0 < F1;
  cout << '\t' << F0 << " < " << F1 << " => " << std::boolalpha << flag << '\n';
  flag = F0 > F1;
  cout << '\t' << F0 << " > " << F1 << " => " << std::boolalpha << flag << '\n';
}

int main() {
    test_arithmetic();
    test_comparison();
}