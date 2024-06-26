#include <gtest/gtest.h>

#include "../include/matdecompose.h"
#include "../include/matfunc.h"
#include "../include/matops.h"
#include "../include/matpred.h"
#include "../include/meshlines.h"
#include "../include/utils.h"
#include "../include/vecops.h"

#include <tuple>
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

TEST(VecOps, ScalarInplace) {
  vector<int> a{1, 2, 3};
  const vector<int> ap{2, 3, 4};
  const vector<int> am{1, 2, 3};
  const vector<int> amul{2, 4, 6};
  const vector<int> adiv{1, 2, 3};

  a += 1;
  EXPECT_EQ(a, ap);
  a -= 1;
  EXPECT_EQ(a, am);
  a *= 2;
  EXPECT_EQ(a, amul);
  a /= 2;
  EXPECT_EQ(a, adiv);
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

#ifdef DEBUG
TEST(VecOps, GetAngle) {
  const vector<int> a{1, 0};
  const vector<int> b{0, 1};

  EXPECT_NEAR(get_angle(a, b), 90, 1e-6);
}
#endif

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

TEST(MatOps, SwapRows) {
  const Matrix<int, 3, 3> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 3, 3> N{{7, 8, 9}, {4, 5, 6}, {1, 2, 3}};

  auto res = M;
  swap_rows(res, 0, 2);

  EXPECT_EQ(res, N);
}

TEST(MatOps, SwapCols) {
  const Matrix<int, 3, 3> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 3, 3> N{{3, 2, 1}, {6, 5, 4}, {9, 8, 7}};

  auto res = M;
  swap_cols(res, 0, 2);

  EXPECT_EQ(res, N);
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

TEST(MatFunc, Transpose) {
  const Matrix<int, 2, 3> a{{1, 3, 5}, {2, 4, 6}};
  const Matrix<int, 3, 2> b{{1, 2}, {3, 4}, {5, 6}};

  Matrix<int, 3, 3> c{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 3, 3> d{{1, 4, 7}, {2, 5, 8}, {3, 6, 9}};

  EXPECT_EQ(transpose(a), b);
  T(c);
  EXPECT_EQ(c, d);
}

TEST(MatFunc, UpperTriangularMatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, 4, 5> RM{{1, 2, 3, 4, 5},
                       {6, 7, 8, 9, 10},
                       {11, 12, 13, 14, 15},
                       {16, 17, 18, 19, 20}};

  const Matrix<int, Rows, Cols> UM{{1, 2, 3}, {0, 5, 6}, {0, 0, 9}};

  const Matrix<int, 4, 5> RUM{
      {1, 2, 3, 4, 5}, {0, 7, 8, 9, 10}, {0, 0, 13, 14, 15}, {0, 0, 0, 19, 20}};

  EXPECT_EQ(triu(M), UM);
  mask_triu(M);
  EXPECT_EQ(M, UM);
  EXPECT_EQ(triu(RM), RUM);
  mask_triu(RM);
  EXPECT_EQ(RM, RUM);
}

TEST(MatFunc, LowerTriangularMatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, 4, 5> RM{{1, 2, 3, 4, 5},
                       {6, 7, 8, 9, 10},
                       {11, 12, 13, 14, 15},
                       {16, 17, 18, 19, 20}};

  const Matrix<int, Rows, Cols> LM{{1, 0, 0}, {4, 5, 0}, {7, 8, 9}};

  const Matrix<int, 4, 5> RLM{{1, 0, 0, 0, 0},
                              {6, 7, 0, 0, 0},
                              {11, 12, 13, 0, 0},
                              {16, 17, 18, 19, 0}};

  EXPECT_EQ(tril(M), LM);
  mask_tril(M);
  EXPECT_EQ(M, LM);
  EXPECT_EQ(tril(RM), RLM);
  mask_tril(RM);
  EXPECT_EQ(RM, RLM);
}
#ifdef DEBUG
TEST(MatFunc, DiagonalMatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, 4, 5> RM{{1, 2, 3, 4, 5},
                       {6, 7, 8, 9, 10},
                       {11, 12, 13, 14, 15},
                       {16, 17, 18, 19, 20}};

  const Matrix<int, Rows, Cols> DM{{1, 0, 0}, {0, 5, 0}, {0, 0, 9}};

  const Matrix<int, 4, 5> RDM{
      {1, 0, 0, 0, 0}, {0, 7, 0, 0, 0}, {0, 0, 13, 0, 0}, {0, 0, 0, 19, 0}};

  EXPECT_EQ(diag(M), DM);
  mask_diag(M);
  EXPECT_EQ(M, DM);
  EXPECT_EQ(diag(RM), RDM);
  mask_diag(RM);
  EXPECT_EQ(RM, RDM);
}
#endif

TEST(MatFunc, Trace) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  EXPECT_EQ(trace(M), 15);
}

TEST(MatFunc, MinorSubmatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<double, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<double, 2, 2> N{{1, 2}, {4, 5}};
  const double minor_N = -3;

  EXPECT_EQ(minor_submatrix(M, 2, 2), N);
  EXPECT_EQ(minor(M, 2, 2), minor_N);
}

TEST(MatFunc, Determinant) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  const Matrix<double, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const double det_M = -9.51619735392994e-16;

  EXPECT_NEAR(det(M), det_M, 1e-6);
}

TEST(MatFunc, MinorMatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  const Matrix<double, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<double, Rows, Cols> N{{-3, -6, -3}, {-6, -12, -6}, {-3, -6, -3}};

  EXPECT_EQ(minor_matrix(M), N);
}

TEST(MatFunc, CofactorSubmatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  const Matrix<double, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<double, Rows - 1, Cols - 1> N{{5, -6}, {-8, 9}};

  EXPECT_EQ(cofactor_submatrix(M, 0, 0), N);
}

TEST(MatFunc, CofactorMatrix) {
  // Also tests cofactor
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  const Matrix<double, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<double, Rows, Cols> N{{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}};
  const auto C = cofactor_matrix(M);

  EXPECT_EQ(C, N);
}

TEST(MatFunc, AdjugateMatrix) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  const Matrix<double, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<double, Rows, Cols> Adj{{-3, 6, -3}, {6, -12, 6}, {-3, 6, -3}};
  const auto A = adj(M);
}

TEST(MatPred, IsLowerTriangular) {
  const Matrix<int, 3, 3> M = {{1, 0, 0}, {4, 5, 0}, {7, 8, 9}};
  const Matrix<int, 4, 5> RM = {{1, 0, 0, 0, 0},
                                {6, 7, 0, 0, 0},
                                {11, 12, 13, 0, 0},
                                {16, 17, 18, 19, 0}};

  EXPECT_TRUE(is_tril(M));
  EXPECT_TRUE(is_tril(RM));

  const Matrix<int, 3, 3> N = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 4, 5> RN = {{1, 2, 3, 4, 5},
                                {6, 7, 8, 9, 10},
                                {11, 12, 13, 14, 15},
                                {16, 17, 18, 19, 20}};

  EXPECT_FALSE(is_tril(N));
  EXPECT_FALSE(is_tril(RN));
}

TEST(MatPred, IsUpperTriangular) {
  const Matrix<int, 3, 3> M = {{1, 2, 3}, {0, 5, 6}, {0, 0, 9}};
  const Matrix<int, 4, 5> RM = {
      {1, 2, 3, 4, 5}, {0, 7, 8, 9, 10}, {0, 0, 13, 14, 15}, {0, 0, 0, 19, 20}};

  EXPECT_TRUE(is_triu(M));
  EXPECT_TRUE(is_triu(RM));

  const Matrix<int, 3, 3> N = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 4, 5> RN = {{1, 2, 3, 4, 5},
                                {6, 7, 8, 9, 10},
                                {11, 12, 13, 14, 15},
                                {16, 17, 18, 19, 20}};

  EXPECT_FALSE(is_triu(N));
  EXPECT_FALSE(is_triu(RN));
}

TEST(MatPred, IsDiagonal) {
  const Matrix<int, 3, 3> M = {{1, 0, 0}, {0, 5, 0}, {0, 0, 9}};
  const Matrix<int, 4, 5> RM = {
      {1, 0, 0, 0, 0}, {0, 7, 0, 0, 0}, {0, 0, 13, 0, 0}, {0, 0, 0, 19, 0}};

  EXPECT_TRUE(is_diag(M));
  EXPECT_TRUE(is_diag(RM));

  const Matrix<int, 3, 3> N = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 4, 5> RN = {{1, 2, 3, 4, 5},
                                {6, 7, 8, 9, 10},
                                {11, 12, 13, 14, 15},
                                {16, 17, 18, 19, 20}};

  EXPECT_FALSE(is_diag(N));
  EXPECT_FALSE(is_diag(RN));
}

TEST(MatPred, IsAntiDiag) {
  const Matrix<int, 3, 3> M = {{0, 0, 1}, {0, 5, 0}, {1, 0, 0}};
  const Matrix<int, 5, 5> MM = {{0, 0, 0, 0, 1},
                                {0, 0, 0, 1, 0},
                                {0, 0, 1, 0, 0},
                                {0, 1, 0, 0, 0},
                                {1, 0, 0, 0, 0}};

  EXPECT_TRUE(is_antidiag(M));
  EXPECT_TRUE(is_antidiag(MM));

  const Matrix<int, 3, 3> N = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 5, 5> NN = {{1, 0, 0, 0, 1},
                                {0, 1, 0, 1, 0},
                                {0, 0, 1, 0, 0},
                                {0, 1, 0, 1, 0},
                                {1, 0, 0, 0, 1}};

  EXPECT_FALSE(is_antidiag(N));
  EXPECT_FALSE(is_antidiag(NN));
}

TEST(MatPred, IsSquare) {
  const Matrix<int, 3, 3> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 4, 5> RM = {{1, 2, 3, 4, 5},
                                {6, 7, 8, 9, 10},
                                {11, 12, 13, 14, 15},
                                {16, 17, 18, 19, 20}};

  EXPECT_TRUE(is_square(M));
  EXPECT_FALSE(is_square(RM));
}

TEST(MatDecompose, LUCrout) {
  const Matrix<int, 3, 3> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, 3, 3> L = {{1, 0, 0}, {4, 1, 0}, {7, 2, 1}};
  const Matrix<int, 3, 3> U = {{1, 2, 3}, {0, -3, -6}, {0, 0, 0}};

  auto [L_res, U_res] = lu_crout(M);

  EXPECT_EQ(L_res, L);
  EXPECT_EQ(U_res, U);
}

TEST(MatDecompose, LUGaussian) {
  const Matrix<double, 3, 3> M = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
  const Matrix<double, 3, 3> L = {{1., 0., 0.}, {4., 1., 0.}, {7., 2., 1.}};
  const Matrix<double, 3, 3> U = {{1., 2., 3.}, {0., -3., -6.}, {0., 0., 0.}};

  const auto &[LG, UG] = lu_gaussian(M);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(UG(i, j), U(i, j), 1e-6);
    }
  }
}

TEST(MatDecompose, PLU) {
  const Matrix<double, 3, 3> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<double, 3, 3> P = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
  const Matrix<double, 3, 3> L = {
      {1, 0, 0}, {0.14285714285714285, 1, 0}, {0.5714285714285714, 0.5, 1}};
  const Matrix<double, 3, 3> U = {
      {7, 8, 9}, {0, 0.8571428571428571, 1.7142857142857142}, {0, 0, 0}};

  const auto [P_res, L_res, U_res] = plu(M);
  const auto PT_res = transpose(P_res);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(PT_res(i, j), P(i, j), 1e-6);
      EXPECT_NEAR(L_res(i, j), L(i, j), 1e-6);
      EXPECT_NEAR(U_res(i, j), U(i, j), 1e-6);
    }
  }
}

// TEST(MatDecompose, QR) {
//   const Matrix<double, 3, 3> M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
//   const Matrix<double, 3, 3> Q = {
//       {-0.12309149097933272, 0.9045340337332909, 0.4082482904638631},
//       {-0.4923659639173309, 0.30151134457776363, -0.8164965809277261},
//       {-0.8616404368553291, -0.30151134457776363, 0.4082482904638631}};
//   const Matrix<double, 3, 3> R = {
//       {-8.12403840463596, -9.60113629638795, -11.07823418813994},
//       {0, 0.9045340337332909, 1.8090680674665818},
//       {0, 0, 0}};

//   const auto [Q_res, R_res] = qr_gm(M);

//   for (int i = 0; i < 3; i++) {
//     for (int j = 0; j < 3; j++) {
//       EXPECT_NEAR(Q_res(i, j), Q(i, j), 1e-6);
//       EXPECT_NEAR(R_res(i, j), R(i, j), 1e-6);
//     }
//   }
// }

TEST(MatDecompose, Cholesky) {
  const Matrix<double, 3, 3> M = {
      {4., 12., -16.}, {12., 37., -43.}, {-16., -43., 98.}};
  const Matrix<double, 3, 3> L = {{2., 0., 0.}, {6., 1., 0.}, {-8., 5., 3.}};

  const auto L_res = cholesky(M);

  EXPECT_EQ(L_res, L);
}