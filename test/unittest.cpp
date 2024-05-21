#include <gtest/gtest.h>

#include "../include/matfunc.h"
#include "../include/matops.h" //resolve compiler error
#include "../include/matpred.h"
#include "../include/vecops.h"
// #include "../include/meshlines.h"

#include <type_traits>
#include <vector>

// todo
// [ ] Matrix Scalar
// [ ] QR Factorization
// [ ] Is Row Echelon predicate

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

TEST(MatOps, UnaryNegative) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  const Matrix<int, Rows, Cols> M1{{1, 2, 3}, {0, 0, 1}, {-1, 0, 2}};
  const Matrix<float, Rows, Cols> M2{{1., 2., 3.}, {4., 5., 6.}, {-0, .5, .7}};

  const Matrix<int, Rows, Cols> E1{{-1, -2, -3}, {0, 0, -1}, {1, 0, -2}};
  const Matrix<float, Rows, Cols> E2{
      {-1., -2., -3.}, {-4., -5., -6.}, {0, -.5, -.7}};

  EXPECT_EQ(-M1, E1);
  EXPECT_EQ(-M2, E2);
}

TEST(MatOps, Equality) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, Rows, Cols> B{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
  const Matrix<int, Rows, Cols> C{{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
  const Matrix<int, Rows, Cols> D{{9, 16, 21}, {24, 25, 24}, {21, 16, 9}};
  const Matrix<int, Rows, Cols> E{
      {30, 24, 18}, {84, 69, 54}, {138, 114, 90}}; // matmul

  EXPECT_EQ(C, A + B);
  EXPECT_EQ(B, C - A);
  EXPECT_EQ(D, A * B);
  EXPECT_EQ(B, D / A);
}

TEST(MatOps, InPlace) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, Rows, Cols> B{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

  const Matrix<int, Rows, Cols> AA{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, Rows, Cols> BB{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
  const Matrix<int, Rows, Cols> C{{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
  const Matrix<int, Rows, Cols> D{{9, 16, 21}, {24, 25, 24}, {21, 16, 9}};
  const Matrix<int, Rows, Cols> E{
      {30, 24, 18}, {84, 69, 54}, {138, 114, 90}}; // matmul

  A += BB;
  EXPECT_EQ(A, C);
  A -= B;
  EXPECT_EQ(A, AA);
  A *= B;
  EXPECT_EQ(A, D);
  A /= B;
  EXPECT_EQ(A, AA);
}

TEST(MatOps, MatVec) {
  const std::size_t Rows = 2;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> A{{1, 2, 3}, {4, 5, 6}};
  vector<int> b{1, 2, 3};

  const Matrix<int, Rows, Cols> B{{2, 4, 6}, {5, 7, 9}};
  const Matrix<int, Rows, Cols> C{{0, 0, 0}, {3, 3, 3}};
  const Matrix<int, Rows, Cols> D{{1, 4, 9}, {4, 10, 18}};
  const Matrix<int, Rows, Cols> E{{1, 1, 1}, {4, 2, 2}};

  // const auto Add = A + b;
  // const auto Sub = A - b;
  const auto Mul = A * b;
  // const auto Div = A / b; // THIS SHOULD BE CASTED TO DOUBLE

  // EXPECT_EQ(Add, B);
  // EXPECT_EQ(Sub, C);
  EXPECT_EQ(Mul, D);
  // EXPECT_EQ(Div, E);
}

TEST(MatOps, MatVecInplace) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 2;

  Matrix<int, Rows, Cols> A{{1, 2, 3}, {4, 5, 6}};
  vector<int> b{1, 2, 3};

  const Matrix<int, Rows, Cols> B{{2, 4, 6}, {5, 7, 9}};
  const Matrix<int, Rows, Cols> C{{1, 2, 3}, {4, 5, 6}};
  const Matrix<int, Rows, Cols> D{{1, 4, 9}, {4, 10, 18}};
  const Matrix<int, Rows, Cols> E{{1, 2, 3}, {4, 5, 6}};

  A += b;
  EXPECT_EQ(A, B);
  A -= b;
  EXPECT_EQ(A, C);
  A *= b;
  EXPECT_EQ(A, D);
  A /= b;
  EXPECT_EQ(A, E);
}

// TEST(MatOps, MatScalar) {
//   Matrix<int> A{{1, 2, 3}, {4, 5, 6}};
//   int a = 3;

//   const Matrix<int> B{{4, 5, 6}, {7, 8, 9}};

//   const auto C = A + a;

//   EXPECT_EQ(C, B);
// }

TEST(MatOps, Elementwise) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, Rows, Cols> B{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

  EXPECT_EQ(A, A);
  EXPECT_NE(A, B);
}

TEST(MatFunc, Transpose) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  const Matrix<int, Rows, Cols> B = transpose(A);

  EXPECT_EQ(A, transpose(transpose(A)));
  T(A);
  EXPECT_EQ(A, B);
}

TEST(MatFunc, LUDecomposition) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<double, Rows, Cols> A{{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  const Matrix<double, Rows, Cols> L{{1, 0, 0}, {1.5, 1, 0}, {2, 4. / 3., 1}};
  const Matrix<double, Rows, Cols> U{
      {4, 3, 1}, {0, -1.5, -0.5}, {0, 0, -1. / 3.}};

  const auto LUCrout = lu_crout(A);
  const auto LUGaussian = lu_gaussian(A);

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < Cols; j++) {
      ASSERT_DOUBLE_EQ(std::get<0>(LUCrout)(i, j), L(i, j));
      ASSERT_DOUBLE_EQ(std::get<1>(LUCrout)(i, j), U(i, j));

      ASSERT_DOUBLE_EQ(std::get<0>(LUGaussian)(i, j), L(i, j));
      // << std::get<0>(LUGaussian);
      ASSERT_DOUBLE_EQ(std::get<1>(LUGaussian)(i, j), U(i, j));
      // << std::get<1>(LUGaussian);
    }
  }
}

#include <cmath>

TEST(MatFunc, QRDecomposition) {
  using std::sqrt;

  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<double, Rows, Cols> A{{4, 3, 1}, {6, 3, 1}, {8, 4, 1}};
  const Matrix<double, Rows, Cols> Q{
      {2. / sqrt(29), 5. / sqrt(29), 0},
      {3. / sqrt(29), (-6. * sqrt(29)) / 145., 4. / 5.},
      {4. / sqrt(29), (-8. * sqrt(29)) / (145.), -3. / 5.}};

  const Matrix<double, Rows, Cols> R{
      {2. * sqrt(29), 31. / sqrt(29), 9. / sqrt(29)},
      {0, 5. / sqrt(29), 11. / (5. * sqrt(29))},
      {0, 0, -1. / 5}};

  const auto QRgm = qr_gm(A);

  for (std::size_t i = 0; i < Rows; i++) {
    for (std::size_t j = 0; j < Cols; j++) {
      EXPECT_NEAR(std::get<0>(QRgm)(i, j), Q(i, j), 1e-9); //<< Q;
      EXPECT_NEAR(std::get<1>(QRgm)(i, j), R(i, j), 1e-9); // << R;
    }
  }
}

TEST(MatFunc, UpperTriangular) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, 5, 4> RM{{1, 2, 3, 4, 5},
                       {6, 7, 8, 9, 10},
                       {11, 12, 13, 14, 15},
                       {16, 17, 18, 19, 20}};

  const Matrix<int, Rows, Cols> UM{{1, 2, 3}, {0, 5, 6}, {0, 0, 9}};
  const Matrix<int, 5, 4> RUM{
      {1, 2, 3, 4, 5}, {0, 7, 8, 9, 10}, {0, 0, 13, 14, 15}, {0, 0, 0, 19, 20}};

  EXPECT_EQ(triu(M), UM) << "wrong triu";
  mask_triu(M);
  EXPECT_EQ(M, UM); // << M;
  mask_triu(RM);
  EXPECT_EQ(RM, RUM);
}

TEST(MatFunc, LowerTringular) {
  const std::size_t Rows = 3;
  const std::size_t Cols = 3;

  Matrix<int, Rows, Cols> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix<int, 5, 4> RM{{1, 2, 3, 4, 5},
                       {6, 7, 8, 9, 10},
                       {11, 12, 13, 14, 15},
                       {16, 17, 18, 19, 20}};

  const Matrix<int, Rows, Cols> LM{{1, 0, 0}, {4, 5, 0}, {7, 8, 9}};
  const Matrix<int, 5, 4> RLM{{1, 0, 0, 0, 0},
                              {6, 7, 0, 0, 0},
                              {11, 12, 13, 0, 0},
                              {16, 17, 18, 19, 0}};

  EXPECT_EQ(tril(M), LM) << "wrong tril";
  mask_tril(M);
  EXPECT_EQ(M, LM) << "M";
  EXPECT_EQ(tril(RM), RLM) << "wrong tril";
  mask_tril(RM);
  EXPECT_EQ(RM, RLM) << "RM";
}

TEST(MatFunc, Trace) {
  const Matrix<int, 3, 3> M{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  EXPECT_EQ(trace(M), 15);
}

TEST(MatFunc, Det) {
  const Matrix<double, 3, 3> M{{4, 3, 2}, {6, 3, 1}, {8, 4, 1}};
  const double det_actual = 2.0;
  const double det_expected = det(M);

  EXPECT_DOUBLE_EQ(det_actual, det_expected);
}

// TEST(MatFunc, RREF) {
//   const Matrix<double, 3, 3> M1{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
//   const Matrix<double, 4, 4> M2{
//       {3., 2., 6., 4.}, {2., 1., 5., 9.}, {6., 5., 7., 0.}, {4., 9.,
//       0., 2.}};
//   const Matrix<double, 5, 5> M3{{1, 4, 7, 2, 9},
//                                 {4, 6, 3, 5, 8},
//                                 {7, 3, 2, 1, 6},
//                                 {2, 5, 1, 9, 0},
//                                 {9, 8, 6, 0, 4}};

//   const Matrix<double, 4, 3> M4{
//       {11, 22, 34, 56}, {23, 44, 67, 78}, {12, 23, 46, 29}};

//   const Matrix<double, 3, 3> E1{{1., 0., -1.}, {0., 1., 2.}, {0., 0., 0.}};
//   const Matrix<double, 4, 4> E2{
//       {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
//   const Matrix<double, 5, 5> E3{{1, 0, 0, 0, 0},
//                                 {0, 1, 0, 0, 0},
//                                 {0, 0, 1, 0, 0},
//                                 {0, 0, 0, 1, 0},
//                                 {0, 0, 0, 0, 1}};

//   const Matrix<double, 4, 3> E4{{1.0, 0, 0, -35.1452282157676},
//                                 {0, 1.0, 0, 21.8879668049793},
//                                 {0, 0, 1.0, -1.14522821576763}};

//   EXPECT_EQ(rref(M1), E1);
//   EXPECT_EQ(rref(M2), E2);
//   EXPECT_EQ(rref(M3), E3);

//   const auto EM4 = rref(M4);
//   const std::size_t row_em4 = EM4.size();
//   const std::size_t col_em4 = EM4.size();

//   for (std::size_t i = 0; i < row_em4; i++) {
//     for (std::size_t j = 0; j < col_em4; j++) {
//       EXPECT_DOUBLE_EQ(EM4[i][j], E4[i][j]);
//     }
//   }
// }

// Note that for now it works on double type; an error occurs when int is
// casted to double and presumably this error is same for integral types
// casted to floating point.
TEST(MatFunc, Minor) {
  const Matrix<double, 3, 3> M1{{1, 4, 7}, {3, 0, 5}, {-1, 9, 11}};
  const Matrix<double, 4, 4> M2{
      {-1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};

  const Matrix<double, 2, 2> E1{{1, 4}, {3, 0}};
  const Matrix<double, 3, 3> E2{{-1, 2, 4}, {5, 6, 8}, {9, 10, 12}};
  const Matrix<double, 2, 2> E3{{2, 4}, {10, 12}};

  EXPECT_EQ(minor_submatrix(M1, 2, 2), E1);
  EXPECT_EQ(minor_submatrix(M2, 3, 2), E2);
  EXPECT_EQ(minor_submatrix(minor_submatrix(M2, 3, 2), 1), E3);

  EXPECT_DOUBLE_EQ(minor(M1, 1, 1), 18.);
  EXPECT_DOUBLE_EQ(minor(M2, 2, 2), 32.);
}
// Note that the M1 that is commented out result in a nan value
// more investigation is needed to resolve this error properly;
// this seems to be a numerical error as opposed to algorithmic error
TEST(MatFunc, CofactorAndAdjugate) {
  // const Matrix<double> M1{{1, 4, 7}, {3, 0, 5}, {-1, 9, 11}};
  const std::size_t Rows1 = 3;
  const std::size_t Cols1 = 3;
  const std::size_t Rows2 = 4;
  const std::size_t Cols2 = 4;

  const Matrix<double, Rows1, Cols1> M1{
      {1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
  const Matrix<double, Rows1, Cols1> M2{
      {2., 3., 4.}, {1., 5., 6.}, {7., 8., 9.}};
  const Matrix<double, Rows1, Cols1> M3{
      {3., 2., 1.}, {6., 5., 4.}, {9., 8., 7.}};

  const Matrix<double, Rows2, Cols2> M4{{-1., 2., 3., 4.},
                                        {5., 6., 7., 8.},
                                        {9., 10., 11., 12.},
                                        {13., 14., 15., 16.}};
  const Matrix<double, Rows2, Cols2> M5{{2., 4., 6., 8.},
                                        {10., 12., 14., 16.},
                                        {18., 20., 22., 24.},
                                        {26., 28., 30., 32.}};
  const Matrix<double, Rows2, Cols2> M6{{12., 11., 13., 15.},
                                        {23., 2., 4., 5.},
                                        {3., 5., 2., 6.},
                                        {16., 34., 56., 11.}};

  // const Matrix<double> E1{{-45, -38, 27}, {19, 18, -13}, {20, 16, -12}};
  const Matrix<double, Rows1, Cols1> E1{
      {-3., 6., -3.}, {6., -12., 6.}, {-3., 6., -3.}};
  const Matrix<double, Rows1, Cols1> E2{
      {-3., 33., -27.}, {5., -10., 5.}, {-2., -8., 7.}};
  const Matrix<double, Rows1, Cols1> E3{
      {3., -6., 3.}, {-6., 12., -6.}, {3., -6., 3.}};

  const Matrix<double, Rows2, Cols2> E4{{0., 0., 0., 0.},
                                        {0., 8., -16., 8.},
                                        {0., -16., 32., -16.},
                                        {0., 8., -16., 8.}};
  const Matrix<double, Rows2, Cols2> E5 =
      zero_mat<double, 4,
               4>(); // not sure if this was the intended initial value
  const Matrix<double, Rows2, Cols2> E6{{1028., 6290., -3191., -4692.},
                                        {-1663., -909., 765., 1334.},
                                        {-1032., -13279., 7571., 4002.},
                                        {-83., -921., -126., 851.}};

  const auto CM1 = cofactor_matrix(M1);
  const auto CM2 = cofactor_matrix(M2);
  const auto CM3 = cofactor_matrix(M3);
  const auto CM4 = cofactor_matrix(M4);
  const auto CM5 = cofactor_matrix(M5);
  const auto CM6 = cofactor_matrix(M6);

  const auto AD1 = adj(M1);
  const auto AD2 = adj(M2);
  const auto AD3 = adj(M3);
  const auto AD4 = adj(M4);
  const auto AD5 = adj(M5);
  const auto AD6 = adj(M6);

  // EXPECT_DOUBLE_EQ(minor(M1, 0, 0), -45.);

  assert(E1.dims() == CM1.dims());
  assert(E2.dims() == CM2.dims());
  assert(E3.dims() == CM3.dims());

  assert(E4.dims() == CM4.dims());
  assert(E5.dims() == CM5.dims());
  assert(E6.dims() == CM6.dims());

  for (std::size_t i = 0; i < CM1.rows(); i++) {
    for (std::size_t j = 0; j < CM1.cols(); j++) {
      EXPECT_NEAR(CM1(i, j), E1(i, j), 1e-9);
      EXPECT_NEAR(CM2(i, j), E2(i, j), 1e-9);
      EXPECT_NEAR(CM3(i, j), E3(i, j), 1e-9);
    }
  }

  EXPECT_EQ(AD1, transpose(CM1));
  EXPECT_EQ(AD2, transpose(CM2));
  EXPECT_EQ(AD3, transpose(CM3));

  for (std::size_t i = 0; i < CM4.rows(); i++) {
    for (std::size_t j = 0; j < CM4.cols(); j++) {
      EXPECT_NEAR(CM4(i, j), E4(i, j), 1e-9);
      EXPECT_NEAR(CM5(i, j), E5(i, j), 1e-9);
      EXPECT_NEAR(CM6(i, j), E6(i, j), 1e-9);
    }
  }

  EXPECT_EQ(AD4, transpose(CM4));
  EXPECT_EQ(AD5, transpose(CM5));
  EXPECT_EQ(AD6, transpose(CM6));
}

TEST(MatPred, LowerTriangular) {
  const Matrix<int, 3, 3> M1{{1, 0, 0}, {4, 5, 0}, {7, 8, 9}};
  const Matrix<int, 3, 3> M2{{1, 0, 0}, {4, 5, 1}, {7, 0, 0}};

  EXPECT_TRUE(is_tril(M1));
  EXPECT_FALSE(is_tril(M2));
}

TEST(MatPred, UpperTriangular) {
  const Matrix<int, 3, 3> M1{{1, 2, 3}, {0, 5, 6}, {0, 0, 9}};
  const Matrix<int, 3, 3> M2{{1, 0, 0}, {4, 5, 1}, {7, 0, 0}};

  EXPECT_TRUE(is_triu(M1));
  EXPECT_FALSE(is_triu(M2));
}

TEST(MatPred, Diagonal) {
  const Matrix<int, 3, 3> M1{{1, 0, 0}, {0, 5, 0}, {0, 0, 9}};
  const Matrix<int, 3, 3> M2{{1, 0, 0}, {4, 5, 1}, {7, 0, 0}};

  EXPECT_TRUE(is_diag(M1));
  EXPECT_FALSE(is_diag(M2));
}

// TEST(MatPred, Echelon) {
//   const Matrix<int, 4, 4> M1{
//       {1, 2, 3, 4}, {0, 1, 2, 3}, {0, 0, 0, 1}, {0, 0, 0, 0}};
//   const Matrix<int, 5, 3> M2{{2, 1, 4, 3, 5}, {0, 0, 1, 2, 3}, {0, 0, 0, 0,
//   1}}; const Matrix<int, 6, 5> M3{{3, 2, 1, 5, 4, 6},
//                              {0, 0, 1, 3, 2, 5},
//                              {0, 0, 0, 0, 1, 2},
//                              {0, 0, 0, 0, 0, 1},
//                              {0, 0, 0, 0, 0, 0}};

//   const Matrix<int, 5, 4> M4{
//       {1, 2, 3, 4, 5}, {0, 0, 1, 2, 3}, {0, 0, 0, 1, 2}, {0, 0, 0, 0, 1}};
//   const Matrix<int, 4, 3> M5{{1, 2, 3, 4}, {0, 1, 2, 3}, {0, 0, 1, 2}};

//   const Matrix<int, 4, 4> N1{
//       {1, 2, 3, 4}, {0, 1, 2, 3}, {0, 0, 0, 1}, {0, 0, 1, 0}};
//   const Matrix<int, 5, 3> N2{{2, 1, 4, 3, 5}, {0, 0, 1, 2, 3}, {1, 0, 0, 0,
//   1}}; const Matrix<int, 6, 5> N3{{3, 2, 1, 5, 4, 6},
//                              {0, 0, 1, 3, 2, 5},
//                              {0, 0, 0, 0, 1, 2},
//                              {0, 0, 0, 0, 0, 1},
//                              {0, 0, 0, 0, 0, 1}};
//   const Matrix<int, 5, 4> N4{
//       {1, 2, 3, 4, 5}, {0, 0, 1, 2, 3}, {0, 0, 0, 1, 2}, {0, 0, 0, 2, 1}};
//   const Matrix<int, 4, 3> N5{{1, 2, 3, 4}, {0, 1, 2, 3}, {0, 3, 1, 2}};

//   EXPECT_TRUE(is_echelon(M1));
//   EXPECT_TRUE(is_echelon(M2));
//   EXPECT_TRUE(is_echelon(M3));
//   EXPECT_TRUE(is_echelon(M4));
//   EXPECT_TRUE(is_echelon(M5));

//   // EXPECT_FALSE(is_echelon(N1));
//   // EXPECT_FALSE(is_echelon(N2));
//   // EXPECT_FALSE(is_echelon(N3));
//   // EXPECT_FALSE(is_echelon(N4));
//   // EXPECT_FALSE(is_echelon(N5));
// }

TEST(MatPred, Symmetric) {
  const Matrix<int, 3, 3> M1{{4, 1, 7}, {1, 6, 8}, {7, 8, 9}};
  const Matrix<int, 4, 4> M2{
      {3, 2, 6, 4}, {2, 1, 5, 9}, {6, 5, 7, 0}, {4, 9, 0, 2}};
  const Matrix<int, 5, 5> M3{{1, 4, 7, 2, 9},
                             {4, 6, 3, 5, 8},
                             {7, 3, 2, 1, 6},
                             {2, 5, 1, 9, 0},
                             {9, 8, 6, 0, 4}};

  const Matrix<int, 3, 3> N1{{4, 1, 10}, {1, 6, 8}, {7, 8, 9}};
  const Matrix<int, 4, 4> N2{
      {3, 2, 11, 4}, {2, 1, 5, 9}, {6, 5, 7, 0}, {4, 9, 0, 2}};
  const Matrix<int, 5, 5> N3{{1, 4, 8, 2, 9},
                             {4, 6, 3, 5, 8},
                             {7, 3, 2, 1, 6},
                             {2, 5, 1, 9, 0},
                             {9, 8, 6, 0, 4}};

  EXPECT_TRUE(is_sym(M1));
  EXPECT_TRUE(is_sym(M2));
  EXPECT_TRUE(is_sym(M3));

  EXPECT_FALSE(is_sym(N1));
  EXPECT_FALSE(is_sym(N2));
  EXPECT_FALSE(is_sym(N3));
}

TEST(MatPred, AntiSymmetric) {
  const Matrix<int, 3, 3> M1{{0, 2, -7}, {-2, 0, 5}, {7, -5, 0}};
  const Matrix<int, 4, 4> M2{
      {0, -4, 6, -1}, {4, 0, -3, 8}, {-6, 3, 0, -2}, {1, -8, 2, 0}};
  const Matrix<int, 5, 5> M3{{0, 5, -1, -7, 4},
                             {-5, 0, -2, 9, -6},
                             {1, 2, 0, -4, -3},
                             {7, -9, 4, 0, 1},
                             {-4, 6, 3, -1, 0}};

  const Matrix<int, 3, 3> N1{{4, 1, 10}, {1, 6, 8}, {7, 8, 9}};

  EXPECT_TRUE(is_antisym(M1));
  EXPECT_TRUE(is_antisym(M2));
  EXPECT_TRUE(is_antisym(M3));

  EXPECT_FALSE(is_antisym(N1));
}