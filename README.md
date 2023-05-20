# LinearAlgebra

This repository is a work in progress. It will contain a set of algorithms for solving some problems in Linear Algebra.

Matrix Decomposition:

- [ ] Singular Value Decomposition
- [x] QR Factorization
- [x] LU Factorization
- [ ] Eigenvalue Decomposition

Matrix Operations:

- [x] Matrix-Matrix operations
- [x] Matrix-Vector operations
- [ ] Matrix-Scalar operations

Algorithms for checking Matrix Properties

- [ ] Determine Invertible Matrix
- [ ] Determine Square Matrix
- [ ] Determine Upper Triangular Matrix
- [ ] Determine Lower Triangular Matrix
- [ ] Find Determinant of Matrix
- [ ] Find Transpose of Matrix
- [ ] Find Inverse of Matrix
- [ ] Find Trace of Matrix

Matrix Object Hierarchy

```mermaid
flowchart TD;
BaseMatrix --> SquareMatrix --> DiagonalMatrix --> IdentityMatrix;
DiagonalMatrix --> ScalarMatrix;

BaseMatrix --> RectangularMatrix --> RowMatrix;
RectangularMatrix --> ColumnMatrix;

BaseMatrix --> SparseMatrix;
BaseMatrix --> ZeroMatrix;
```
