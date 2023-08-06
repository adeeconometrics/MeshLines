# MeshLines
This repository is a work in progress. It will contain a set of algorithms for solving some problems in Linear Algebra.

<p align="center">
  <img src="img/MeshLinesLogo.png" alt="Image">
</p>

## Motivation, Goals, and Disclaimers
I started this project to educate myself on writing linear algebra routines and extending these algorithms to a compatible type.
I intend to learn about matrix factorization and develop insights to elide operations given their type information e.g. computing $M \times I = M$ without a cost or $M\times M^{-1} = I$.

Note that there is no long-term initiative to support and resolve bugs on this project. The robustness of the codebase relies on the test coverage. 
You are welcome to submit an issue or drop a suggestion for this repo. 

---
Matrix Decomposition:

- [ ] Singular Value Decomposition
- [x] QR Factorization
- [x] LU Factorization
- [x] LDL Decomposition
- [x] PLU Decomposition
- [ ] Eigenvalue Decomposition

Matrix Operations:

- [x] Matrix-Matrix operations
- [x] Matrix-Vector operations
- [x] Transpose
- [x] Determinant
- [x] Minor
- [x] Cofactor
- [x] Adjugate
- [ ] Matrix-Scalar operations
- [ ] Inverse
- [ ] Colspace
- [ ] Rowspace
- [ ] Nullspace
- [ ] Span

Algorithms for checking Matrix Properties

- [x] Determine Invertible Matrix
- [x] Determine Square Matrix
- [x] Determine Upper Triangular Matrix
- [x] Determine Lower Triangular Matrix
- [x] Find Trace of Matrix
