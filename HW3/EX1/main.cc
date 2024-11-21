#include <cmath>
#include <vector>
#include <iostream>

#include "matrix.h"
#include "jacdiag.h"

int kronecker(int n, int m) {

  if (n == m) {
  
    return 1;
    
  } else {
  
    return 0;
    
  }
  
}

using namespace std; 

int main() {

  //initialize matrix size and # of sweeps to perform
  int size;
  int sweeps = 10;

  //user inputs desired matrix size
  cout << "matrix size: " << endl;
  cin >> size;

  //initialize eigenvalue vector and matrix to be diagonalized
  Row d(size, 0.0);
  Matrix H(size, Row(size, 0.0));

  //fill the matrix
  for (int n = 1; n <= size; n++) {     

    for (int m = 1; m <= size; m++) {
    
      H[n-1][m-1] = kronecker(n,m) * 4 * (n*n) + min(n,m) * (0.1 + pow(-1, abs(m-n)) * 10);

    }
    
  }

  //perform jacobi diagonalization
  jacdiag(H, d, sweeps);
  
}