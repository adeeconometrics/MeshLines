#ifndef __TEST_H__
#define __TEST_H__

#include "../include/factorize.h"
#include "../include/helper.h"

#include <iostream>

auto test_print() -> void {
  // Test the print_matrix function
  //   std::vector<std::vector<int>> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  //   print_matrix(A);

  Matrix<double> B = {{1.23456789, 2.3456789, 1.34, 3.1415926, 1000000000},
                      {3.456789, 4.56789, 1.3, 4},
                      {5.6789, 6.789, 5.5556, 6.777},
                      {7.89, 8.9, 0, 0}};
  print_matrix(B);

  //   std::vector<std::vector<std::string>> C = {
  //       {"abc", "defgh"}, {"ijklm", "nop"}, {"qrs", "tuvwx"}};
  //   print_matrix(C);
}

auto test_lu() -> void {
  const Matrix<double> A = {{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  auto [L, U] = lu_crout(A);

  std::cout << "L: \n";
  print_matrix(L);
  std::cout << "U: \n";
  print_matrix(U);
}

auto test_qr() -> void {
  const Matrix<double> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  // const Matrix<double> A = {{0.90259351, 0.63733912, 0.06712906},
  //                           {0.89503327, 0.98548892, 0.17569676},
  //                           {0.9062746, 0.44367133, 0.36358112}};

  auto [Q_gm, R_gm] = qr_gm(A);
  auto [Q_h, R_h] = qr_householder(A);

  std::cout << "Q: \n";
  print_matrix(Q_gm);
  std::cout << "R: \n";
  print_matrix(R_gm);

  std::cout << "-----------------\n";

  std::cout << "Q: \n";
  print_matrix(Q_h);
  std::cout << "R: \n";
  print_matrix(R_h);
}

#endif // __TEST_H__