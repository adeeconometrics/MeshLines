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

## Design Principles

Design for genericity and specialize later on performance

- Add template mechanisms for interopt: it should be able to work with types that meets the operations required by Matrix semantics.
- Zero cost abstraction: Imposing rules and expectations should come without cost
- Minimal redundancies

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
- [x] Matrix-Scalar operations
- [ ] Inverse
- [ ] Colspace
- [ ] Rowspace
- [ ] Nullspace
- [ ] Span

Matrix Decompositions

| Method                             | Description                                                                                                                                    | Type of Matrix              | Status |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- | ------ |
| LU Decomposition                   | Factors a matrix into a product of a lower triangular matrix $L$ and an upper triangular matrix $U$.                                           | Square                      |        |
| QR Decomposition                   | Factors a matrix into a product of an orthogonal matrix $Q$ and an upper triangular matrix $R$.                                                | Square, Rectangular         |        |
| Cholesky Decomposition             | Factors a symmetric positive-definite matrix into a product of a lower triangular matrix $L$ and its transpose $𝐿^T$.                          | Symmetric Positive-Definite |        |
| Singular Value Decomposition (SVD) | Factors a matrix into a product of three matrices: $U$ (orthogonal), $\sum$ (diagonal), and $VT$ (orthogonal).                                 | Square, Rectangular         |        |
| Eigenvalue Decomposition           | Decomposes a matrix into its eigenvalues and eigenvectors.                                                                                     | Square                      |        |
| Schur Decomposition                | Decomposes a matrix into a product of an orthogonal matrix $Q$ and an upper triangular matrix $T$.                                             | Square                      |        |
| Hessenberg Decomposition           | Decomposes a matrix into a product of an orthogonal matrix $Q$ and a Hessenberg matrix $H$.                                                    | Square                      |        |
| Jordan Decomposition               | Decomposes a matrix into its Jordan normal form.                                                                                               | Square                      |        |
| Polar Decomposition                | Decomposes a matrix into a product of a unitary matrix $U$ and a positive semi-definite matrix $P$.                                            | Square                      |        |
| Bunch-Kaufman Decomposition        | Factors a symmetric or Hermitian matrix into a product involving a permutation matrix, a lower triangular matrix, and a block diagonal matrix. | Symmetric, Hermitian        |        |
