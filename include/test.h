#ifndef __TEST_H__
#define __TEST_H__

#include "../include/matfunc.h"
#include "../include/utils.h"
#include "../include/vecfunc.h"
#include "../include/vecops.h"

#include <cstdlib>
#include <iostream>
#include <vector>

auto test_print() -> void {
  // Test the std::cout << function
  std::vector<std::vector<int>> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  //   std::cout << A;

  Matrix<double> B = {{1.23456789, 2.3456789, 1.34, 3.1415926, 1000000000},
                      {3.456789, 4.56789, 1.3, 4},
                      {5.6789, 6.789, 5.5556, 6.777},
                      {7.89, 8.9, 0, 0}};
  std::cout << B;

  //   std::vector<std::vector<std::string>> C = {
  //       {"abc", "defgh"}, {"ijklm", "nop"}, {"qrs", "tuvwx"}};
  //   std::cout << C;
}

auto test_lu() -> void {
  const Matrix<double> A = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  auto [L, U] = lu_crout(A);

  std::cout << "L: \n";
  std::cout << L;
  std::cout << "U: \n";
  std::cout << U;
}

auto test_qr() -> void {
  const Matrix<double> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  // const Matrix<double> A = {{0.90259351, 0.63733912, 0.06712906},
  //                           {0.89503327, 0.98548892, 0.17569676},
  //                           {0.9062746, 0.44367133, 0.36358112}};

  auto [Q_gm, R_gm] = qr_gm(A);
  auto [Q_h, R_h] = qr_householder(A);

  std::cout << "Q: \n";
  std::cout << Q_gm;
  std::cout << "R: \n";
  std::cout << R_gm;

  std::cout << "-----------------\n";

  std::cout << "Q: \n";
  std::cout << Q_h;
  std::cout << "R: \n";
  std::cout << R_h;
}

auto test_ch() -> void {
  const Matrix<double> A = {{2.0, -1.0, 0}, {-1, 2, -1}, {0, -1, 2}};

  const auto L = cholesky(A);

  std::cout << "L: \n";
  std::cout << L;
}

auto test_vecops() -> void {
  using namespace lin;
  std::vector<int> a{};
  std::vector<float> b{};

  for (size_t i{}; i < 10; i++) {
    a.emplace_back(std::rand() % 100);
    b.emplace_back(std::rand() % 100);
  }
  // b.emplace_back(3);

  std::cout << a;
  std::cout << b;

  const auto add = a + b;
  const auto sub = a - b;
  const auto mul = a * b;
  const auto div = a / b;

  const auto sadd = a + 5.5;
  std::cout << "---\n";
  std::cout << add;
  std::cout << sub;
  std::cout << mul;
  std::cout << div;
  std::cout << sadd;

  std::cout << "dist : " << dist(a);
};

auto test_matops() -> void {
  using namespace lin;
  Matrix<float> A = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  Matrix<float> B = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};

  // todo
  // - [ ] mat [op] vec
  // - [ ] mat [op] scalar

  // const float a = 5;
  // const auto C = A + a;
  // A += B;
  std::cout << A;
}

#endif // __TEST_H__
