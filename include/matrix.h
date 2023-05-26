#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>

namespace lin {

template <typename T> using Matrix = std::vector<std::vector<T>>;

}
#endif // __MATRIX_H__