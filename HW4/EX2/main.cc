#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility> 

#include "functions.h"

using namespace std;                    

int main() {

  //run params
  int L = 1000;
  double epsilon = 1e-4;
  double diff = 100;

  Matrix grid (L, Row(L, 0.0));
  Matrix rgrid (L, Row(L, 0.0));
 
  //populate grid and rgrid
  for (int i=0; i<L; i++) {
  
    grid[L-1][i] = 100; 
    rgrid[L-1][i] = 100;
    
  }
  
  //initialize counter
  int count = 0;
  
  while (diff > epsilon) {
    
    //calculate max difference for current iteration
    double diff_max = 0; 

    for (int i = L - 2; i >= 1; i--) {
    
      for (int j = L - 2; j>=1; j--) {
      
        rgrid[i][j] = 0.25 * (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1]); 

        //set diff_max to grid_diff if it is smaller than it
        double grid_diff = abs(grid[i][j] - rgrid[i][j]); 
        
        if (grid_diff > diff_max) {
        
          diff_max = grid_diff;
          
        }
        
      }
      
    }
    
    diff = diff_max;
    count++;
    
    //report the max difference for every 5000 iterations
    if (count % 5000 == 0) {
    
      cout << "Iteration " << count << " max diff = " << diff << endl;
      
    }

    swap(grid, rgrid); 
    
  }
  
  cout << "Max Difference threshold reached" << endl;
  
  //save matrix to output file
  ofstream outFile("output.txt");
  
  size_t size = rgrid.size();
  
  for (size_t i = 0; i < L; ++i) {
  
    for (size_t j = 0; j < L; ++j) {
    
      outFile << rgrid[i][j];
      
      if (j < size - 1) {
      
        outFile << "\t";
        
      }
      
    } 
    
    outFile << "\n";
    
  }
  
  outFile.close();
  
}
