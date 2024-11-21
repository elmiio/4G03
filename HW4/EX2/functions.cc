#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

//set potential according to conditions outlined in handout
double V (double x) {
  
  if (abs(x) >= 5) {
  
    return 0;
    
  } 
  
  else {
  
    return x*x/5 - 14.0;
    
  }
  
}

double V2 (double x) {
  
  if (abs(x) >= M_PI/4) {
  
    return 1e20;
    
  } 
  
  else {

    return (0.05 / pow(sin(x+M_PI/4),2)) + (5 / pow(cos(x+M_PI/4),2));
    
  }
  
}

//setting k^2 equal to E - V^2
double k2 (double x, double E) {

  return E - V2(x);
  
}

//wavefunction calculation when integrating from left to right
double phi_l (double x, double h, double current, double prev, double E) {
  
  return (2 * (1.0 - 5.0/12 * h*h * k2(x-h, E)) * current - (1.0 + h*h / 12.0 * k2(x-2*h, E)) * prev)/(1 + h*h / 12.0 * k2(x, E));
  
}

//wavefunction calculation when integrating from right to left
double phi_r (double x, double h, double current, double prev, double E) {

  return (2 * (1.0 - 5.0/12 * h*h * k2(x+h, E)) * current - (1.0 + h*h / 12.0 * k2(x+2*h, E)) * prev)/(1 + h*h / 12.0 * k2(x, E)); 
  
}

//function for saving data
void save(const Row& array1, const Row& array2, const string filename) {
   
  ofstream outFile(filename, ios::trunc);
  
  size_t size = array1.size();
    
  for (size_t i = 0; i < size; ++i) {
  
    outFile << array1[i] << "\t" << array2[i] << endl;
    
  }

  outFile.close();
  
}

//function for normalization
void norm(Row& vec, double h) {

  double sum = 0;
  int Nsteps = vec.size();
  
  for (int i=0; i < Nsteps-1; i++) {
  
    sum += (vec[i] + vec[i+1])*h/2;
    
  }
  
  for (int i=0; i < Nsteps; i++) {
  
   vec[i] = vec[i]/sum;
   
  }
  
}

//function for merging phi_l and phi_r together
void merge(Row& phi, Row& phil, Row& phir, Row& xl, Row& xr, double xm, double E, double Nsteps, double h) {
  
  for (int j=2; j<Nsteps; j++) {
  
    phil[j] = phi_l(xl[j], h, phil[j-1], phil[j-2], E);
    phir[j] = phi_r(xr[j], h, phir[j-1], phir[j-2], E);
    
  }

  for (int j=0; j<Nsteps; j++) {
  
    phil[j] = phil[j] * phir[xm] / phil[xm];
    
    if (j <= xm) {
    
      phi[j] = phil[j];
      
    } 
    
    else {
    
      phi[j] = phir[Nsteps-j-1];
    
    }
    
  }
  
}