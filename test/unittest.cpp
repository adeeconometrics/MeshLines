#include <gtest/gtest.h>

#include "../include/matfunc.h"
#include "../include/matops.h" //resolve compiler error
#include "../include/matpred.h"
#include "../include/meshlines.h"
#include "../include/vecops.h"

#include <type_traits>
#include <vector>

using namespace lin;
using std::vector;

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

TEST(VecOps, UnaryNegative) {
  const vector<int> v1{1, 2, 3, 4, -4, 0};
  const vector<float> v2{1., 2., 3., -3., -5., 0.};

  const vector<int> e1{-1, -2, -3, -4, 4, 0};
  const vector<float> e2{-1., -2., -3., 3., 5., 0.};

  EXPECT_EQ(-v1, e1);
  EXPECT_EQ(-v2, e2);
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
  EXPECT_EQ(a, b); // a = 2,4,6
  a -= d;
  EXPECT_EQ(a, d); // a = 1,2,3
  b /= a;
  EXPECT_EQ(b, c); // b = 2,2,2
  b *= d;
  EXPECT_EQ(b, e); // c = 2,4,6
  a += 1;
  EXPECT_EQ(a, f); // a = 2,3,4
  a -= 1;
  EXPECT_EQ(a, d); // a = 1,2,3
  a *= 2;
  EXPECT_EQ(a, e); // a = 2,4,6
  a /= 2;
  EXPECT_EQ(a, d); // a = 1,2,3
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

TEST(VecOps, Sum) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};

  EXPECT_EQ(sum(a), 6);
  EXPECT_EQ(sum(b), 15);
}

TEST(VecOps, Prod) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};

  EXPECT_EQ(prod(a), 6);
  EXPECT_EQ(prod(b), 120);
}

TEST(VecOps, DotProduct) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};

  EXPECT_EQ(dot(a, b), 32);
}

TEST(VecOps, CrossProduct) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};
  const vector<int> c{-3, 6, -3};

  EXPECT_EQ(cross(a, b), c);
}

TEST(VecOps, Dist) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};

  EXPECT_NEAR(dist(a), 3.7416573867739413, 1e-6);
  EXPECT_NEAR(dist(b), 8.774964387392123, 1e-6);
}

TEST(VecOps, LPNorm) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};

  EXPECT_NEAR(lp_norm(a, 1), 6, 1e-6);
  EXPECT_NEAR(lp_norm(b, 1), 15, 1e-6);
  EXPECT_NEAR(lp_norm(a, 2), 3.7416573867739413, 1e-6);
  EXPECT_NEAR(lp_norm(b, 2), 8.774964387392123, 1e-6);
}

TEST(VecOps, GetAngle) {
  const vector<int> a{1, 0};
  const vector<int> b{0, 1};

  EXPECT_NEAR(get_angle(a, b), 90, 1e-6);
}

TEST(VecOps, Normalize) {
  const vector<int> a{1, 2, 3};
  const vector<int> b{4, 5, 6};

  const vector<double> a_norm{0.26726124, 0.53452248, 0.80178373};
  const vector<double> b_norm{0.45584231, 0.56980288, 0.68376346};

  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(normalize(a)[i], a_norm[i], 1e-6);
    EXPECT_NEAR(normalize(b)[i], b_norm[i], 1e-6);
  }
}

TEST(MatOps, UnaryNegative) {
  const Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> b{{-1, -2}, {-3, -4}};

  EXPECT_EQ(-a, b);
}

TEST(MatOps, Equality) {
  const Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> b{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> c{{2, 3}, {4, 5}};

  EXPECT_TRUE(a == a);
  EXPECT_TRUE(a == b);
  EXPECT_TRUE(a != c);
}

TEST(MatOps, Elementwise) {
  const Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> b{{2, 3}, {4, 5}};
  const Matrix<int, 2, 2> c{{3, 5}, {7, 9}};
  const Matrix<int, 2, 2> d{{-1, -1}, {-1, -1}};
  const Matrix<int, 2, 2> e{{2, 6}, {12, 20}};
  const Matrix<double, 2, 2> f{{0.5, 2.0 / 3.0}, {0.75, 4.0 / 5.0}};

  EXPECT_EQ(a + b, c);
  EXPECT_EQ(a - b, d);
  EXPECT_EQ(a * b, e);

  const auto t = a / b;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_NEAR(t(i, j), f(i, j), 1e-6);
    }
  }
}

TEST(MatOps, MatrixInplace) {
  Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  Matrix<int, 2, 2> b{{2, 4}, {6, 8}};
  const Matrix<int, 2, 2> c{{2, 2}, {2, 2}};
  const Matrix<int, 2, 2> d{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> e{{1, 4}, {9, 16}};
  const Matrix<int, 2, 2> f{{2, 3}, {4, 5}};

  a += a;
  EXPECT_EQ(a, b);
  a -= d;
  EXPECT_EQ(a, d);
  b /= c;
  EXPECT_EQ(b, d);
  b *= d;
  EXPECT_EQ(b, e);
}

TEST(MatOps, MatrixVectorBroadcasting) {
  const Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const vector<int> b{1, 2};

  const Matrix<int, 2, 2> c{{2, 4}, {4, 6}};
  const Matrix<int, 2, 2> d{{0, 0}, {2, 2}};

  EXPECT_EQ(a + b, c);
  EXPECT_EQ(a - b, d);
}

TEST(MatOps, MatrixVectorMultiplication) {
  const Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const vector<int> b{1, 2};

  const vector<int> c{5, 11};

  EXPECT_EQ(a * b, c);
}

TEST(MatOps, MatrixScalar) {
  const Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> b{{2, 3}, {4, 5}};
  const Matrix<int, 2, 2> c{{0, 1}, {2, 3}};
  const Matrix<int, 2, 2> d{{2, 4}, {6, 8}};

  EXPECT_EQ(a + 1, b);
  EXPECT_EQ(b - 1, a);
  EXPECT_EQ(a * 2, a + a);
  EXPECT_EQ(d / 2, a);
}

TEST(MatOps, MatrixScalarInplace) {
  Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Matrix<int, 2, 2> b{{2, 3}, {4, 5}};
  const Matrix<int, 2, 2> c{{0, 1}, {2, 3}};
  const Matrix<int, 2, 2> d{{0, 2}, {4, 6}};

  a += 1;
  EXPECT_EQ(a, b);
  a -= 2;
  EXPECT_EQ(a, c);
  a *= 2;
  EXPECT_EQ(a, d);
  a /= 2;
  EXPECT_EQ(a, c);
}