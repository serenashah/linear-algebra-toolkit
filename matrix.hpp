//Serena Shah, ss94574
#include <vector>
#include <string>
#include "vector.hpp"
#include "timer.hpp"
#pragma once

using namespace std;

//----------------------------------------------
// matrix support class definition (matrix.hpp)
//----------------------------------------------

// supported solver modes
enum solveModes {GAUSS_ELIM = 1, JACOBI = 2, GAUSS_SEIDEL = 3};
enum fitModes {LINEAR = 1, EXP = 2, POWER = 3};

class Matrix
{
private:
  int M_;
  int N_;
  vector< vector<double> > myMatrix_;
  int beenCalled_ = 0;
  int maxIters_;
  int currentIter_;
  double tol_;
  Vector b_;
  Timer atimer_;
  Vector guess_;
  bool debugMode_;
  
public:
  Matrix();                                  // constructor
  void initSize();                           // initialize the matrix size
  void initIdentity(int n);                  // initialize as an identity matrix of size nxn
  void initFromFile(string fileName);        // read the matrix from a file
  bool isSquare();                           // test whether matrix is square
  int numRows();                             // return number of rows
  int numCols();                             // return number of columns
  double getVal(int row, int col);           // return matrix value at given row/col location (0-indexed based)
  void setVal(int row, int col, double val); // set the matrix value at given row/col location (0-index based)
  void initNewMatrix();
  Matrix Multiply(Matrix B);                 // post-multiply by B and return resulting matrix
  Matrix Multiply(double A);                 // multiply by a scalar and return resulting matrix
  Matrix Transpose();                        // return transpose of matrix
  vector<double> Diagonal();                 // return a vector containing diagonal elements of the matrix
  void Print();                              // print the matrix to stdout
  void Print(string name);                   // print the matrix to stdout with a name prefix

  // Linear solver support
  Vector Solve(Vector b, int mode);          // return solution of [A]x=b using desired solver mode
  Vector getNewB();                          // return the Vector at some given point
  double getSolveTime();                     // wall clock time (in secs) required for last solve
  void setSolveDebugMode(bool flag);         // set flag to toggle debug output mode for linear solves

  // Support methods for iterative methods
  void setSolveMaxIters(int iters);          // set cap on max # of iterations
  void setSolveTolerance(double tol);        // set desired stopping tolerance
  int  getSolveIters();                      // return number of iters completed from last solve

  // regression support (for Nx2 matrices)
  Vector linFit(int mode);                     // linear regression fit (mode=LINEAR, EXPONENTIAL, or POWER)
  Vector expTransform(Vector y);               // maps the points to exponential model
  Vector powTransform(Vector b);               // maps the points to power law model
  Matrix evalLinFit(int mode, Vector fit);     // evaluate linear fit - ouput is an Nx2 matrix

  // additional utilities to aid in regression
  Vector extract(int col);                     // extract vector from matrix corresponding to the col index
  void saveToFile(string fileName);            // save matrix to file (same file format as initFromFile)

};
