#include <fstream>
#include <limits>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "matrix.h"
#include "jacdiag.h"
#include "hv.h"

using namespace std;
using namespace chrono;

int main() {     
  
  //random number generation
  mt19937 generator;
  uniform_real_distribution<double> dist(0,1.0);
  
  //system size
  int L;  

  //number of sweeps to perform
  int Nsweeps = 10;
  
  //lanczos matrix size, to be specified by user
  int m;

  cout << "Size of system, L: ";
  cin >> L;
  cout << "Size of Lanczos matrix, m: ";
  cin >> m;

  //number of states
  int N = pow(2,L);
  
  //initialize lanczos matrix and vectors
  Matrix Lan(m, Row(m, 0.0)), v(m, Row(m));
  Row v0(N), v1(N), f(N), omega(N);

  //randomize + normalize v0
  for (int i = 0; i < N; i++){
  
    v0[i] = 1.0 - 2.0*dist(generator);
    
  }
      
  v[0] = normalize(v0);
  
  //get initial lanczos params
  hv(omega, v[0], L);
  Lan[0][0] = dotprod(v[0], omega);
  f = subtract(omega, scale(v[0], Lan[0][0]));

  //lanczos iteration (eq. 23)
  for (int j = 0; j < m - 1; j++) { 
  
    Lan[j+1][j] = Lan[j][j+1] = norm(f);
    v[j+1] = normalize(f);
    
    hv(omega, v[j+1], L);
    omega = subtract(omega, scale(v[j], norm(f)));
    
    Lan[j+1][j+1] = dotprod(v[j+1], omega);
    f = subtract(omega, scale(v[j+1], Lan[j+1][j+1]));
    
  }
  
  //vector to store eigenvalues
  Row d(m, 0.0);
  
  //apply jacobi diagonalization to lanczos matrix 
  jacdiag(Lan, d, Nsweeps);

  //get E0 and E1
  sort(d.begin(), d.end());
  
  double E0 = d[0];
  double E1 = d[1];
  
  cout << "E0 = " << E0 << endl;
  cout << "E1 = " << E1 << endl;
  
  //write results to output file
  ofstream outFile("output.txt", ios::app);
  
  if (outFile.is_open()) {
  
          outFile << E0 << " " << E1 << "\n";
          
              outFile.close();
              
  } else {
  
          cout << "Unable to open file for writing.\n";
          
  }

  return 0;
  
}