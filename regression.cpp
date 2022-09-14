//Serena Shah, ss94574
#include <iostream>
#include <vector>
#include <typeinfo>
#include <unistd.h>
#include "vector.hpp"
#include "matrix.hpp"
#include "timer.hpp"
using namespace std;

void usageMessage(){
  cout << "----------------------------------" << endl;
  cout << "Data Fitting Wizard 9002" << endl;
  cout << "    -> A linear regression fitting tool!" << endl;
  cout << "----------------------------------" << endl << endl;
  cout << "Usage: regression [mode] [matrix.input] [output]" << endl << endl;
  cout << "where:" << endl;
  cout << "  [mode]            fit method (1=linear 2=exponential, 3=power)"<< endl;
  cout << "  [matrix.input]    name of text input file housing matrix [A]" << endl;
  cout << "  [output]          name of text output file with fitting results (Nx2 matrix)" << endl;
  exit(1);
  }

int main(int argc, char *argv[]){

  if (argc != 4){
    usageMessage();
  }

  int mode = atoi(argv[1]);
  string matrixFile = argv[2];
  string output = argv[3];

  Matrix A;
  A.initFromFile(matrixFile);
  Vector coefficients = A.linFit(mode);
  A.evalLinFit(mode, coefficients).saveToFile(output);
  cout << "--> Fitted output saved to file: " << output << endl;
  return 0;
}
