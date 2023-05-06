#ifndef __FACTORIZE_H__
#define __FACTORIZE_H__

#include <tuple>
#include <vector>

using std::tuple;
using std::vector;
template <typename T> using Matrix = vector<vector<T>>;

template <typename T>
auto lu_decomposition(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>>;

template <typename T>
auto qr_factorization(const Matrix<T> &A) -> tuple<Matrix<T>, Matrix<T>>;

template <typename T>
auto ldl(const Matrix<T> &A) -> tuple<Matrix<T>, vector<T>>;

#endif // __FACTORIZE_H__