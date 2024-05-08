# Matmul Bench

This subfolder contains different implementation of matmul.

Goals:

- Understand when to use and how matmul alrgorithms effect performance
- Gain a perspective to which alrgorithm is suitable for expression templates

Guiding Principles

- Should work with arm64 and x86.
- Should check for alrgorithm correctess.
- Take generic examples first before using cpu intrinsics.
- Use standard cpp over third party packages
- Use `-O0` before turning compiler optimization
