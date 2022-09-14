//Serena Shah, ss94574
#include "matrix.hpp"
#include "vector.hpp"
#include "timer.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;

Matrix::Matrix() :myMatrix_(){}

void Matrix::initSize(){
  if (M_ <=0 || N_ <= 0){
    cout << "[Error]: Matrix size must be greater than 0" << endl;
    exit(1);
    }
  
  myMatrix_.resize(M_);
  for (int j = 0; j < M_; j++){
    myMatrix_[j].resize(N_);
  }
  
  beenCalled_ = 1;
}
  
void Matrix::initFromFile(string fileName){  
  ifstream matrixFile (fileName);
  if (matrixFile.is_open()){
    matrixFile >> M_;
    matrixFile >> N_;
  }

  initSize();

  if (matrixFile.is_open()){
    for (int i = 0; i < M_; i++){
      for (int j = 0; j < N_; j++){
	matrixFile >> myMatrix_[i][j];
      }
    }
  }
}


void Matrix::initIdentity(int n){
  M_ = n;
  N_ = n;

  initSize();

  for (int i = 0; i < M_; i++){
    for (int j = 0; j < N_; j++){
      if (i == j){
	myMatrix_[i][j] = 1;
      }
    }
  }  
}

bool Matrix::isSquare(){  
  if (M_ == N_){
    return true;
  }
  else{
    return false;
  }
}

int Matrix::numRows(){
  return M_;
}

int Matrix::numCols(){
  return N_;
}

double Matrix::getVal(int row, int col){
  if (row < 0 || col < 0){
    cout << "[Error]: Row/column index must be greater than or equal to 0." << endl;
    exit(1);
  }

  return myMatrix_[row][col];
}

void Matrix::setVal(int row, int col, double val){
  myMatrix_[row][col] = val;
}

Matrix Matrix::Multiply(Matrix B){
  if (N_ != B.numRows()){
    cout << numCols() << " by "<< B.numRows() << endl;
    cout << "[Error]: Matrices' inner dimensions must be equal." << endl;
    exit(1);
  }
 
  Matrix matrProduct;
  matrProduct.M_ = M_;
  matrProduct.N_ = B.N_;
  matrProduct.initSize();
  
  for (int i = 0; i < matrProduct.M_; i++){
    for (int j = 0; j < matrProduct.N_; j++){
      double sum = 0;
      for (int k = 0; k  < N_; k++){
	double product = getVal(i, k) * B.getVal(k, j);
	sum = sum + product;
	matrProduct.setVal(i, j, sum);
      }
    }
  }

  return matrProduct;
}

Matrix Matrix::Multiply(double A){
  Matrix multMatrix;
  multMatrix.N_ = N_;
  multMatrix.M_ = M_;
  multMatrix.initSize();
  
  for (int i = 0; i < multMatrix.M_; i++){
    for (int j = 0; j < multMatrix.N_; j++){
      multMatrix.setVal(i, j, myMatrix_[i][j]*A);
    }
    }

  myMatrix_ = multMatrix.myMatrix_;
  return multMatrix;
}

Matrix Matrix::Transpose(){
  Matrix transMatrix;
  transMatrix.M_ = N_;
  transMatrix.N_ = M_;
  transMatrix.initSize();
  cout << transMatrix.myMatrix_.size() << transMatrix.numRows() << endl;

  for (int i = 0; i < transMatrix.M_; i++){
    for (int j = 0; j < transMatrix.N_; j++){                                        	
         transMatrix.setVal(i, j, myMatrix_[j][i]);
	}
      }

  int a = M_;
  M_ = N_;
  N_ = a;
  myMatrix_ = transMatrix.myMatrix_;
  return transMatrix;
}

vector <double> Matrix::Diagonal(){
  vector <double> diagElement_;
  for (int i = 0; i < M_; i++){
    for (int j = 0; j < N_; j++){
      if (i == j){
	diagElement_.push_back(myMatrix_[i][j]);
      }
    }
  }
  return diagElement_;
}

void Matrix::Print(){
  if (beenCalled_ != 1){
    cout << "[Error]: Matrix must be initialized" << endl;
    exit(1);
  }

  for (int i = 0; i < M_; i++){
    cout << "| " ;
    for (int j = 0; j < N_; j++){
      cout << right << setw(8) << myMatrix_[i][j] << " ";
    }
    cout << "|" << endl;
  }
  cout << endl;
}

void Matrix::Print(string name){
  if (beenCalled_ != 1){
    cout << "[Error]: Matrix must be initialized" << endl;
    exit(1);
  }
  cout << "[" << name << "] = " ;  
    for (int i = 0; i < M_; i++){
      cout << "| " ;
    for (int j = 0; j < N_; j++){
      cout << right << setw(10) << setprecision(5) << myMatrix_[i][j] << " ";
    }
    cout << "|" << endl << setw(name.length()+7);
  }
    cout << setw(1) << endl;
}

// Linear solver support  
Vector Matrix:: Solve(Vector b, int mode){
  atimer_.Start();

  Vector x;
  Vector x_old;
  Vector x_new;
  Vector x_diff;
  x.allocateData(b.numElems());
  x_old.allocateData(b.numElems());
  x_new.allocateData(b.numElems());
  x_diff.allocateData(b.numElems());
  double norm;
  solveModes type = static_cast<solveModes>(mode);
  string modeName;
  if (type == 1){
    modeName = "naive Gaussian elimination";
  }
  if (type == 2){
    modeName = "Jacobi";
  }
  if (type == 3){
    modeName = "Gauss-Seidel";
  }

  
  if (M_ != N_ || M_ != b.numElems()){
    cout << "[Error]: Check matrix or vector dimensions." << endl;
    exit(1);
  }

  cout << endl << "Solving [A]x = b using [" << modeName << "]!" << endl;
  
  switch (type){
  case GAUSS_ELIM:
    if (debugMode_ == true){
      cout << "Initial system:" << endl << endl;
      Print("A");
      b.Print("b");
    }
    //forward elim
    for (int k = 0; k < M_-1; k++){
      for (int i = 1 + k; i < M_; i++){ 
	double multiplier = myMatrix_[i][k] / myMatrix_[k][k];
	for (int j = 0; j < N_; j++){
	  double newRowVal = multiplier * myMatrix_[k][j];
	  myMatrix_[i][j] = myMatrix_[i][j] - newRowVal;
	}
	double newVecVal = multiplier * b.getVal(k);
	b.setVal(i, b.getVal(i) - newVecVal);
      }
    }
    if (debugMode_ == true){
      cout << "Updated system after forward elimination:" << endl << endl;
      Print("A");
      b.Print("b");
    }
    
    //backwards sub
    for (int w = b.numElems()-1; w >= 0; w--){
      double term = 0;
      double a_ = myMatrix_[w][w];
      double b_ = b.getVal(w);
      for (int z = b.numElems()-1; z > w; z--){
	term += x.getVal(z) * myMatrix_[w][z];
      }
      double x_ = (b_ - term) / a_;
      x.setVal(w, x_);
    }
    break;
    
  case JACOBI:
    for (currentIter_ = 0; currentIter_ <= maxIters_; currentIter_++){
      Matrix updated;
      updated.M_ = M_;
      updated.N_ = N_;
      updated.initSize();
      updated.myMatrix_ = myMatrix_;
      for (int i = 0; i < M_; i++){
	double sum = 0;
	for (int j = 0; j < N_; j++){
	  if (i != j){
	    updated.myMatrix_[i][j] = x_old.getVal(j) * myMatrix_[i][j];
	    sum += updated.myMatrix_[i][j];
	  }
	}
	double xVal = (b.getVal(i) - sum)/ updated.myMatrix_[i][i];
	x_new.setVal(i, xVal);
      }
      for (int w = 0; w < M_; w++){
	x_diff.setVal(w, x_old.getVal(w) - x_new.getVal(w));
      }
      norm = x_diff.l2norm() / x_old.l2norm();
      if (x_old.l2norm() == 0){
	norm = 1.000;
      }
      x_old = x_new;
      if (debugMode_ == true){
	cout << "  --> Iteration: " << right << setw(4) << currentIter_ + 1 << " | Norm = " << setprecision(5) << norm << endl;
      }
      if (norm < tol_){
	x = x_old;
	cout << "Converged!" << endl;
	break;
	}
      if (currentIter_ == maxIters_){
	cout << "Did not converge." << endl;
	exit(1);
      }
    }
    break;

  case GAUSS_SEIDEL:
    for (currentIter_ = 0; currentIter_ <= maxIters_; currentIter_++){
      Matrix updated;
      updated.M_ = M_;
      updated.N_ = N_;
      updated.initSize();
      updated.myMatrix_ = myMatrix_;
      Vector x_old_real = x_old;
      for (int i = 0; i < M_; i++){
	double sum = 0;
	for (int j = 0; j < N_; j++){
	  if (i != j){
	    updated.myMatrix_[i][j] = x_old.getVal(j) * myMatrix_[i][j];
	    sum += updated.myMatrix_[i][j];
	  }
	}
	double xVal = (b.getVal(i) - sum)/ updated.myMatrix_[i][i];
	x_new.setVal(i, xVal);
	x_old.setVal(i, xVal);
      }
      for (int w = 0; w < M_; w++){
	x_diff.setVal(w, x_old_real.getVal(w) - x_new.getVal(w));
      }
      norm = x_diff.l2norm() / x_old.l2norm();
      if (x_old.l2norm() == 0){
	norm = 1.000;
      }
      x_old = x_new;
      if (debugMode_ == true){
	cout << "  --> Iteration: " << right << setw(4) << currentIter_  + 1 << " | Norm = " << setprecision(5) << norm << endl;
      }
      if (norm < tol_){
	x = x_old;
	cout << "Converged!" << endl;
	break;
      }
      if (currentIter_ == maxIters_){
	cout << "Did not converge." << endl;
	exit(1);
      }
    }
    break;
  }
  atimer_.Stop();
  
  cout << "[N = " << M_ << "] Time to solve using [" << modeName << "] = " << getSolveTime() << " (secs)";
  if (type != 1){
    cout << ", (" << currentIter_ + 1 << " iters)";
  }
  cout << endl << endl;
  if (debugMode_ == true){
    x.Print("x solution");
  }
  
  return x;
}

double Matrix:: getSolveTime(){
  return atimer_.ElapsedTime();
}

void Matrix:: setSolveDebugMode(bool flag){
  debugMode_ = flag;
}

// Support methods for iterative methods
void Matrix:: setSolveMaxIters(int iters){
  maxIters_ = iters;
}

void Matrix:: setSolveTolerance(double tol){
  tol_ = tol;
}

int Matrix:: getSolveIters(){
  int currentIter_ = 0;
  return currentIter_;
}

//linear regression solvers
Vector Matrix:: linFit(int mode){
  if (N_ != 2){
    cout << "[Error]: Only valid for two-point coordinates." << endl;
    exit(1);
  }
  cout << endl << "Performing regression using " << M_ << " points" << endl;
  Vector linFit;
  Vector x;
  Vector y;
  linFit.allocateData(2);
  x.allocateData(M_);
  y.allocateData(M_);
  x = extract(0);
  y = extract(1);
  double a0, a1, a1_numerator, a1_denominator;
  fitModes type = static_cast<fitModes>(mode);

  switch (type){
  case LINEAR:
    cout << "--> Standard linear regression" << endl << endl;
    a1_numerator = (M_ * x.sum(y)) - (x.sum() * y.sum());
    a1_denominator = (M_ * x.sumSquared()) - (pow(x.sum(), 2));
    a1 = a1_numerator / a1_denominator;
    linFit.setVal(1, a1);
    a0 = y.average() - a1 * x.average();
    linFit.setVal(0, a0);
    break;

  case EXP:
    y = expTransform(y);
    cout << "--> Exponential model regression." << endl << endl;
    a1_numerator = (M_ * x.sum(y)) - (x.sum() * y.sum());
    a1_denominator = (M_ * x.sumSquared()) - (pow(x.sum(), 2));
    a1 = a1_numerator / a1_denominator;
    linFit.setVal(1, a1);
    a0 = exp(y.average() - a1 * x.average());
    linFit.setVal(0, a0);
    break;

  case POWER:
    x = expTransform(x);
    y = expTransform(y);
    cout << "--> Power model regression." << endl << endl;
    a1_numerator = (M_ * x.sum(y)) - (x.sum() * y.sum());
    a1_denominator = (M_ * x.sumSquared()) - (pow(x.sum(), 2));
    a1 = a1_numerator / a1_denominator;
    linFit.setVal(1, a1);
    a0 = exp(y.average() - a1 * x.average());
    linFit.setVal(0, a0);
    break;
  }
  
  linFit.Print("Fitting Coefficients");
  return linFit;
}

Vector Matrix:: expTransform(Vector y){
  Vector ln_y;
  ln_y.allocateData(y.numElems());
  for (int i = 0; i < ln_y.numElems(); i++){
    ln_y.setVal(i, log(y.getVal(i)));
  }
  return ln_y;
}

Matrix Matrix:: evalLinFit(int mode, Vector fit){
  Matrix evalLinFit;
  evalLinFit.M_ = M_;
  evalLinFit.N_ = N_;
  evalLinFit.initSize();
  evalLinFit.myMatrix_ = myMatrix_;

  fitModes type = static_cast<fitModes>(mode);

  switch (type){
  case LINEAR:
    for (int i = 0; i < M_; i++){
      double y = fit.getVal(0) + fit.getVal(1) * myMatrix_[i][0];
      evalLinFit.setVal(i, 1, y);
    }
    break;
    
  case EXP:
    for (int i = 0; i < M_; i++){
      double y = fit.getVal(0) * exp(fit.getVal(1) * myMatrix_[i][0]);
      evalLinFit.setVal(i, 1, y);
    }
    break;
    
  case POWER:
    for (int i = 0; i < M_; i++){
      double y = fit.getVal(0) * pow(myMatrix_[i][0], fit.getVal(1));
      evalLinFit.setVal(i, 1, y);
    }
    break;
  }
  return evalLinFit;
}

//additional utilities to aid in regression
Vector Matrix:: extract(int col){
  if (col >= N_){
    cout << "[Error]: Invalid column." << endl;
    exit(1);
  }
  Vector matrixCol;
  matrixCol.allocateData(M_);
  for (int i = 0; i < M_; i++){
    matrixCol.setVal(i, myMatrix_[i][col]);
  }
  return matrixCol;
}

void Matrix:: saveToFile(string fileName){
  if (beenCalled_ != 1){
    cout << "[Error]: Matrix must be initialized" << endl;
    exit(1);
  }
  ofstream file;
  file.open(fileName);
  file << M_ << " " << N_ << endl;
  for (int i = 0; i < M_; i++){
    for (int j = 0; j < N_; j++){
      file << myMatrix_[i][j] << " ";
    }
    file << endl;
  }
  file.close();
}
