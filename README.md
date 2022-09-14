# linear-algebra-toolkit
This object-oriented C++ project provides multiple linear algebra functionalities at large scale and high efficiency. 

# Linear System Solving 
Allows user to solve a linear system by reading in a matrix and right-hand side vector to solve [A]x = b for the x vector with various methods and the time associated with the solve. 
Solving methods include: 
- Gaussian Elimination
- Jacobi (iterative)
- Gauss-Seidel (iterative)
Also allows user relevant solver controls: 
- setting tolerance
- setting maximum iterations
- running a debug mode for Gaussian Elimination, providing a detailed output of the matrix decomposition

# Matrix Operations 
- initializing matrix
- retrieving/setting matrix value
- multiplying matrix by scalar or matrix
- transposing matrix
- retrieving vector of diagonal elements

# Vector Operations
- initializing vector
- retrieving/setting vector value
- finds l2 norm
