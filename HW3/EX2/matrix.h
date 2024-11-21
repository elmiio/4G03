#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

int sign (double x);

Matrix transpose (Matrix &A);
Matrix multiply (const Matrix &A, const Matrix &B);
void rotate (Matrix &A, int p, int q, int size, Row &d);

double calc_S0(Matrix &A);
double calc_epsilon(double S0, int n);

//function for calculating the dot product of two vectors
inline double dotprod (const Row &y, const Row &x) {
 
  double sum = 0.0;
  for (int i = 0; i < y.size(); i++) {
  
    sum += x[i]*y[i];
    
  }
  
  return sum;
}

//scaling function
inline Row scale (Row &v, double scalar) {

  for (int i = 0; i < v.size(); i++) {
  
    v[i] *= scalar;
    
  }
  
  return v;
  
}

//subtraction function
inline Row subtract (const Row &x, const Row &y) {

  Row diff(x.size());

  for (int i = 0; i < x.size(); i++) {
  
    diff[i] = x[i] - y[i];
    
  }
  
  return diff;
  
}

//normalization functions
inline double norm (const Row &v) {

  return sqrt(dotprod(v,v));
  
}

inline Row normalize(Row &v) {

  return scale(v, 1/norm(v));
  
}

#endif