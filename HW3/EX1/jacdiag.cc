#include <cmath>
#include <vector> 
#include <iostream>
#include <algorithm>

#include "matrix.h"

using namespace std;

//jacobi diagonalization function
void jacdiag(Matrix &A, Row &d, int Nsweeps) {

  //size of matrix
  int size = A.size();
  
  //initialize sum of squares variable
  double sum = 0;
  
  //for counting number of sweeps performed
  int count = 0;

  for (int i = 0; i < Nsweeps; i++) {
  
    for (int p = 0; p < size; p++) {
    
      for (int q = p+1; q < size; q++) {    

        //for sweeps before #3, apply jacobi rotation only when matrix element is larger than epsilon
        if (count < 3) {
        
          double S0 = calc_S0(A);
          double epsilon = calc_epsilon(S0, size);

          //rotate if A > epsilon
          if (A[p][q] > epsilon) {
          
            rotate(A, p, q, size, d);
            
          }
          
        }
        
        //for sweeps after #3, set off-diagonal elements to 0 if they are much smaller than diagonal elements, apply jacobi rotation otherwise
        else if (count > 3) {
                 
          if (abs(A[p][q]) < abs(A[p][p])/1e13 && abs(A[q][q]) < abs(A[q][q])/1e13) { 
          
            A[p][q] = 0; 
            
          } else {
          
            rotate(A, p, q, size, d);
            
          }
          
        }
       
        //rotate normally if sweep # is 3
        else { 
        
          rotate(A, p, q, size, d);
          
        }
        
      }    
      
    }  

    //calculate sum of squares above diagonal
    for (int p = 0; p < size; p++) {
    
      for (int q = p + 1; q < size; q++) {
      
        sum += A[p][q] * A[p][q];
        
      }
      
    }
    
    count += 1;

    //show sum of squares after each sweep
    cout << "Sweep #" << count << ": Sum of squares = " << sum << endl;
    
    //sort and display the three smallest eigenvalues
    if (sum == 0 || count == Nsweeps) {
    
        sort(d.begin(), d.end());
        cout << "Lowest 3 Eigenvalues: ";
        for (int i = 0; i < 3; i++) {
        
            cout << d[i] << " ";
            
        }
        
        cout << endl;
        
    }
    
    //stop if off-diagonal elements vanish
    if (sum == 0) { break; }
    
    //reset
    sum = 0;
    
  }
  
}