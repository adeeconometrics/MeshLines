List of Functions to Implement 



AbstractMatrix

```
transpose(Matrix) -> void
transpose_matrix(Matrix) -> Matrix
trace(Matrix) -> number

apply(function,Matrix) -> Matrix
apply_to_matrix(function,Matrix) -> void

is_zero_matrix(Matrix) -> bool
is_ones_matrix(Matrix) -> bool

col_vector(Matrix,int,int) -> vector
row_vector(Matrix,int,int) -> vector

to_col_vectors(Matrix) -> List[vector]
to_row_vectors(Matrix) -> List[vector]

to_rref(Matrix) -> Matrix
to_ref(Matrix) -> Matrix

eigenvalue(Matrix) -> number
eigenvector(Matrix) -> vector
rowspace(Matrix) -> List[vector]
nullspace(Matrix) -> List[vector]
colspace(Matrix) -> List[vector]
rank(Matrix) -> Integer
dims(Matrix) -> (Integer, Integer)
```



- SquaredMatrix

  ```
  solve_QR(SquareMatrix,vector) -> vector*
  solve_LU(SquareMatrix,vector) -> vector*
  solve_GE(SquareMatrix,vector) -> vector*
  
  LU_factorization(SquareMatrix) -> Pair[SquareMatrix,SquareMatrix]
  QR_factorization(SquareMatrix) -> Pair[SquareMatrix,SquareMatrix]
  LDL_factorization(SquareMatrix) -> Pair[SquareMatrix,SquareMatrix]
  cholesky_factorization(SquareMatrix) -> Pair[SquareMatrix,SquareMatrix]
  
  inverse_QR(SquareMatrix) -> SquareMatrix
  inverse_LU(SquareMatrix) -> SquareMatrix
  inverse_LDL(SquareMatrix) -> SquareMatrix
  inverse_GE(SquareMatrix) -> SquareMatrix
  inverse_adj(SquareMatrix) -> SquareMatrix
  
  rank(SquareMatrix) -> int
  rank_decomposition(SquareMatrix) -> pair[SquareMatrix,SquareMatrix]
  
  minor(SquareMatrix) -> number
  minor_matrix(SquareMatrix,int,int) -> SquareMatrix
  cofactor(SquareMatrix) -> number
  cofactor_matrix(SquareMatrix,int,int) -> SquareMatrix
  
  inv(SquareMatrix) -> SquareMatrix
  det(SquareMatrix) -> number
  adj(SquareMatrix) -> SquareMatrix
  
  is_invertible(SquareMatrix) -> bool
  is_indefinite(SquaredMatrix) -> bool
  
  is_echelon(SquareMatrix) -> bool
  is_symmetric(SquareMatrix) -> bool
  is_anti_symmetric(SquareMatrix) -> bool
  is_skew_symmetric(SquareMatrix) -> bool
  
  is_triangular(SquareMatrix) -> bool
  is_lower_triangular(SquareMatrix) -> bool
  is_upper_triangular(SquareMatrix) -> bool
  is_diagonal(SquareMatrix) -> bool
  is_diagonizable(SquareMatrix) -> bool
  
  is_positive_definite(SquareMatrix) -> bool
  is_positive_semidefinite(SquareMatrix) -> bool
  is_negative_definite(SquareMatrix) -> bool
  is_negative_semidefinite(SquareMatrix) -> bool

  is_linearly_dependent(SquareMatrix) -> bool
  ```

- RectangularMatrix

  ```
  
  ```

- IrregularMatrix

  ```
  SparseMatrix(m:int,n:int):
  	add_element(element:T,int,int)
  	modify_element(element:T,int,int)
  	get_element(element:T,int,int)
  ```

  

Vector

```
vectorize(vector, function) -> vector

norm(vector) -> vector
dot(vector, vector) -> number
cross(vector,vector) -> vector
lp_norm(vector,float) -> number
sum(vector) -> number
prod(vector) -> number
get_angle(vector, vector) -> number
get_projection(vector, vector) -> vector

normalize(vector) -> vector
min(vector) -> number
min(vector) -> number

gradient_vector(vector) -> vector
direactional_derivative(vector) -> vector
curl(vector) -> vector
divergence(vector) -> vector

to_matrix(List[vector]) -> Matrix
```

