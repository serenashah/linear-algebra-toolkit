//Serena Shah, ss94574
#include "vector.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
using namespace std;

Vector::Vector() :myVector_(){}

void Vector::allocateData(int nelems){
  elem_ = nelems;
  myVector_.resize(elem_);
}

void Vector::initFromFile(string fileName){
  ifstream vectorFile (fileName);
  if (vectorFile.is_open()){
    vectorFile >> elem_;
    allocateData(elem_);

    for (int i = 0; i < elem_; i++){
      vectorFile >> myVector_[i];
    }
  }  
}

int Vector::numElems(){
  return elem_;
}

double Vector::getVal(int i){
  return myVector_[i];
}

double Vector::l2norm(){
  double sum = 0;
  for (int i = 0; i < elem_; i++){
    sum += pow(myVector_[i], 2);
  }
  sum = sqrt(sum);
  return sum;
}

void Vector::setVal(int i, double val){
  myVector_[i] = val;
}

void Vector::setAllVals(double val){
  for (int i = 0; i < elem_; i++){
    myVector_[i] = val;
  }
}

void Vector::Print(){
  for (int i = 0; i < elem_; i++){
    cout << "| " << right << setw(8) << setprecision(5) << myVector_[i] << "|" << endl;
  }
  cout << endl;
}

void Vector::Print(string name){
  cout << name << " = " ;
  for (int i = 0; i < elem_; i++){
    cout << "| " << right << setw(8) << setprecision(5) << myVector_[i] << "|" << endl << setw(name.length()+5);
    }
  cout << setw(1) << endl;
}

double Vector::sum(){
  double sum = 0;
  for (int i = 0; i < elem_; i++){
    sum += myVector_[i];
  }
  return sum;
}

double Vector::sumSquared(){
  double sum = 0;
  for (int i = 0; i < elem_; i++){
    sum += pow(myVector_[i], 2);
  }
  return sum;
}

double Vector::sum(Vector y){
  if (elem_ != y.elem_){
    cout << "[Error]: Vector lengths don't match." << endl;
    exit(1);
  }
  double sum = 0;
  for (int i = 0; i < elem_; i++){
    sum += myVector_[i] * y.myVector_[i];
  }
  return sum;
}

double Vector:: average(){
  return sum() / elem_;
}
