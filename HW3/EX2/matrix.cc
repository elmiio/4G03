#include <cmath>
#include <vector> 
#include <iostream>
#include "matrix.h" 

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

//sign calculation (pos or neg)
int sign (double x) {

  if (x > 0) {
  
    return 1;
    
  } else {
  
    return -1;
    
  }
  
}

//function for calculating transpose of a matrix
Matrix transpose (Matrix &A) {

  int size = A.size();

  Matrix AT(size, Row(size, 0.0));

  for (int i = 0; i < size; ++i) {
  
    for (int j = 0; j < size; ++j) {
    
      AT[j][i] = A[i][j];
      
    }
    
  }

  return AT;
  
}

//matrix multiplication function
Matrix multiply (const Matrix &A, const Matrix &B) {

  int size = A.size();
  Matrix C(size, Row(size, 0.0));

  for (int i = 0; i < size; ++i) {
  
    for (int j = 0; j < size; ++j) {
    
      for (int k = 0; k < size; ++k) {
      
        C[i][j] += A[i][k] * B[k][j];
        
      }
      
    }
    
  }

  return C;
  
}

//function for rotating a matrix
void rotate (Matrix &A, int p, int q, int size, Row &d) {

  Matrix P(size, Row(size, 0.0));

  double alpha = (A[q][q] - A[p][p]) / (2 * A[p][q]);
  double t = sign(alpha) / ( abs(alpha) + sqrt(alpha*alpha + 1) );
  double c = 1 / sqrt(t*t + 1);
  double s = t * c;

  for (int i = 0; i < size; i++) {
  
    P[i][i] = 1;
    
  }

  P[p][p] = c;
  P[q][q] = c;
  P[p][q] = s;
  P[q][p] = -s; 

  Matrix D(size, Row(size, 0.0));
  D = multiply (transpose(P), multiply(A, P));

  for (int i = 0; i < size; i++) {
  
    d[i] = D[i][i];
    
  }

  A = D;
  
}

//function to calculate S0
double calc_S0(Matrix &A) {

  int size = A.size();
  double sum = 0;

  for (int i = 1; i < size; i++) {
  
    for (int j = i+1; j < size; j++) {
    
      sum += abs(A[i][j]);
      
    }
    
  }

  return sum; 
  
}

//function for calculating epsilon
double calc_epsilon(double S0, int n) { 

  return (1/5) * S0 / (n*n);
  
}