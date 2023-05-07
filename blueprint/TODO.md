- [x] Come up with a cool project name [MeshLines]
- [x] Initialize Doxygen
- [ ] Include Doxygen generation in CMake workflow
- [ ] Refector project file structure
  - [x] Vector.h
  - [ ] SquareMatrix.h
  - [ ] VectorPredicates.h
  - [ ] SquareMatrixPredicate.h
- [ ] Define TeX in code documentation
- [ ] Change `vector` to `vec`_ and `matrix` to `mat`_

---

## Flow

- [ ] Tier 1: Correctness
  - [ ] Write correct functions
  - [ ] Write unit tests
  - [ ] Write CI/CD workflow
- [ ] Tier 2: Implement Predicates and Concepts
  - [ ] Use static metaprogramming to implement Matrix concepts
  - [ ] Make custom operator overloading
- [ ] Tier 3: Performance Optimization
  - [ ] Re-implement the source code with STL functions, and light-weight routines that allows the compiler to optimize
  - [ ] Perform profiling and performance comparison

High Prio

- [ ] Google unit test framework
- [ ] Google benchmark
- [ ] Typesafe Matrix
- [ ] Documentation in Sphinx

Low Prio

- [ ] OpenMP directives
- [ ] CMake Configuration
- [ ] Docker Image
