# Optimized Linear Algebra Toolkit
This object-oriented C++ project provides multiple linear algebra functionalities at large scale and high efficiency. 

## Linear System Solver
Allows user to solve a linear system by reading in a matrix and right-hand side vector to solve [A]x = b for the x vector with various methods and the time associated with the solve. 
Solving methods include: 
- Gaussian Elimination
- Jacobi (iterative)
- Gauss-Seidel (iterative) 
Also allows user relevant solver controls: 
- setting tolerance
- setting maximum iterations
- running a debug mode for Gaussian Elimination, providing a detailed output of the matrix decomposition 

## Regression Fit Tool
Allows user to perform a linear regression fit on a collection of data points using various fitting methods. 
Fitting methods include: 
- Linear
- Exponential 
- Power 
Repo includes plots of the fit compared to the data.

## Further Functionalities
### Matrix Operations 
- initializing matrix
- retrieving/setting matrix value
- multiplying matrix by scalar or matrix
- transposing matrix
- retrieving vector of diagonal elements
- retrieve column vector specified

### Vector Operations
(Not std vector) 
- initializing vector
- retrieving/setting vector value
- finds l2 norm
- return sum of vectors squared
