#include <iostream>
#include "Test.h"
#include "SquareMatrix.h"

using std::cout; 
using tools::assert;

void test_arithmetic(){

}

template <typename T, size_t N>
void test_comparison(const SquareMatrix<T,N>& lhs, const SquareMatrix<T,N>& rhs){
    assert((lhs == rhs) == false, "failed == ");
    assert((lhs != rhs) == true, "failed != ");
    cout << "passed comparison";
}

int main(){
    const SquareMatrix<int,3> M0 {{{{1,2,3},{4,5,6},{7,8,9}}}};
    // const SquareMatrix<int,3> M1 {{{2,4,6},{8,2,4},{6,8,2}}};
    // const vector<int,3> v0 {1,2,3};
    
    cout << det(M0);
}