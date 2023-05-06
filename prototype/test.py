from typing import List, Union
import numpy as np
from scipy.linalg import ldl

A = [[4, 12, -16], [12, 37, -43], [-16, -43, 98]]
B = [[4.0, 12.0, -16.0], [12.0, 37.0, -43.0], [-16.0, -43.0, 98.0]]

def test_ldl(A:Union[np.ndarray, List]) -> None:
    lu, d, perm = ldl(A)
    print(f'lu: {lu}\nd: {d}\nperm:{perm}')

if __name__ == '__main__':
    test_ldl(A)
    test_ldl(B)