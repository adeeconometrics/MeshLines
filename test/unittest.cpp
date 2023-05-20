#include "../include/helper.h"
// #include "../include/matops.h" -- find out why the compiler complains when I
// include this
#include "../include/vecops.h"
#include <gtest/gtest.h>

#include <type_traits>

using namespace lin;

TEST(VecOps, FusedTypes) {
  using std::is_same;

  const vector<int> a{1, 2, 3};
  const vector<float> b{1.0f, 2.0f, 3.0f};

  const auto c_add = a + b;
  const auto c_sub = a - b;
  const auto c_mul = a * b;
  const auto c_div = a / b;
  const auto c_sadd = a + 2.0f;
  const auto c_ssub = a - 2.0f;
  const auto c_smul = a * 2.0f;
  const auto c_sdiv = a / 2.0f;

  static_assert(is_same<decltype(b), decltype(c_add)>::value);
  static_assert(is_same<decltype(b), decltype(c_sub)>::value);
  static_assert(is_same<decltype(b), decltype(c_mul)>::value);
  static_assert(is_same<decltype(b), decltype(c_div)>::value);

  static_assert(is_same<decltype(b), decltype(c_sadd)>::value);
  static_assert(is_same<decltype(b), decltype(c_ssub)>::value);
  static_assert(is_same<decltype(b), decltype(c_smul)>::value);
  static_assert(is_same<decltype(b), decltype(c_sdiv)>::value);
}

TEST(VecOps, Elementwise) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};
  const vector<int> c{5, 7, 9};
  const vector<int> a2{1, 4, 9};
  const vector<int> one{1, 1, 1};

  EXPECT_EQ(a + b, c);
  EXPECT_EQ(c - b, a);
  EXPECT_EQ(a * a, a2);
  EXPECT_EQ(a / a, one);
}

TEST(VecOps, InPlace) {
  vector<int> a{1, 2, 3};
  vector<int> b{2, 4, 6};
  const vector<int> c{2, 2, 2};
  const vector<int> d{1, 2, 3};
  const vector<int> e{2, 4, 6};
  const vector<int> f{2, 3, 4};

  a += a;
  EXPECT_EQ(a, b) << a; // a = 2,4,6
  a -= d;
  EXPECT_EQ(a, d) << a; // a = 1,2,3
  b /= a;
  EXPECT_EQ(b, c) << b; // b = 2,2,2
  b *= d;
  EXPECT_EQ(b, e) << b; // c = 2,4,6
  a += 1;
  EXPECT_EQ(a, f) << a; // a = 2,3,4
  a -= 1;
  EXPECT_EQ(a, d) << a; // a = 1,2,3
  a *= 2;
  EXPECT_EQ(a, e) << a; // a = 2,4,6
  a /= 2;
  EXPECT_EQ(a, d) << a; // a = 1,2,3
}

TEST(VecOps, Scalar) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{3, 4, 5};

  EXPECT_EQ(a + 2, b);
  EXPECT_EQ(b - 2, a);
  EXPECT_EQ(a * 2, a + a);
  EXPECT_EQ(a / 1, a);
}

TEST(VecOps, Equality) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{2, 3, 4};

  EXPECT_TRUE(a == a);
  EXPECT_TRUE(a != b);
}

TEST(MatOps, Equality) {
  Matrix<int> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int> B{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
  const Matrix<int> C{{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
  const Matrix<int> D{{9, 16, 21}, {24, 25, 24}, {21, 16, 9}};
  const Matrix<int> E{{30, 24, 18}, {84, 69, 54}, {138, 114, 90}}; // matmul

  EXPECT_EQ(C, A + B);
  EXPECT_EQ(B, C - A);
  EXPECT_EQ(D, A * B);
  EXPECT_EQ(B, D / A);
}

TEST(MatOps, InPlace) {
  Matrix<int> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int> B{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

  const Matrix<int> AA{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int> BB{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
  const Matrix<int> C{{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
  const Matrix<int> D{{9, 16, 21}, {24, 25, 24}, {21, 16, 9}};
  const Matrix<int> E{{30, 24, 18}, {84, 69, 54}, {138, 114, 90}}; // matmul

  A += BB;
  EXPECT_EQ(A, C) << A;
  A -= B;
  EXPECT_EQ(A, AA) << A;
  A *= B;
  EXPECT_EQ(A, D) << A;
  A /= B;
  EXPECT_EQ(A, AA) << A;
}

TEST(MatOps, MatVec) {
  Matrix<int> A{{1, 2, 3}, {4, 5, 6}};
  vector<int> b{1, 2, 3};

  const Matrix<int> B{{2, 4, 6}, {5, 7, 9}};
  const Matrix<int> C{{0, 0, 0}, {3, 3, 3}};
  const Matrix<int> D{{1, 4, 9}, {4, 10, 18}};
  const Matrix<int> E{{1, 1, 1}, {4, 2, 2}};

  const auto Add = A + b;
  const auto Sub = A - b;
  const auto Mul = A * b;
  const auto Div = A / b; // THIS SHOULD BE CASTED TO DOUBLE

  EXPECT_EQ(Add, B) << Add;
  EXPECT_EQ(Sub, C) << Sub;
  EXPECT_EQ(Mul, D) << Mul;
  EXPECT_EQ(Div, E) << Div;
}

TEST(MatOps, MatVecInplace) {
  Matrix<int> A{{1, 2, 3}, {4, 5, 6}};
  vector<int> b{1, 2, 3};

  const Matrix<int> B{{2, 4, 6}, {5, 7, 9}};
  const Matrix<int> C{{1, 2, 3}, {4, 5, 6}};
  const Matrix<int> D{{1, 4, 9}, {4, 10, 18}};
  const Matrix<int> E{{1, 2, 3}, {4, 5, 6}};

  A += b;
  EXPECT_EQ(A, B) << A;
  A -= b;
  EXPECT_EQ(A, C) << C;
  A *= b;
  EXPECT_EQ(A, D) << A;
  A /= b;
  EXPECT_EQ(A, E) << A;
}

// TEST(MatOps, MatScalar) {
//   Matrix<int> A{{1, 2, 3}, {4, 5, 6}};
//   int a = 3;

//   const Matrix<int> B{{4, 5, 6}, {7, 8, 9}};

//   const auto C = A + a;

//   EXPECT_EQ(C, B);
// }

TEST(MatOps, Elementwise) {
  Matrix<int> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int> B{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

  EXPECT_EQ(A, A);
  EXPECT_NE(A, B);
}
