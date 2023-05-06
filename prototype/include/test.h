#ifndef __TEST_H__
#define __TEST_H__

#include "../include/factorize.h"
#include "../include/helper.h"

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

#endif // __TEST_H__