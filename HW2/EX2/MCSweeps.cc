#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

#include "MCSweeps.h" 
#include "MCvar.h"

using namespace std;

void MCSweeps(vector<vector<int>>& arr, int Nmcs, double J) {

  //store E vals
  double Eneg8, Eneg6, Eneg4, Eneg2, E2, E4, E6, E8;
   
  ifstream inFile("table.txt");
  
  inFile >> E8;
  inFile >> E6;
  inFile >> E4;
  inFile >> E2;
  inFile >> Eneg2;
  inFile >> Eneg4;
  inFile >> Eneg6;
  inFile >> Eneg8;
            
  inFile.close();
  
  //lattice size
  int L = arr.size();
  
  //random number generation
  random_device r;
  mt19937 e2(r());
  uniform_int_distribution<int> unif_dist(0, L-1);
  uniform_real_distribution<> r_dist(0,1);

  for (int i = 0; i < Nmcs; i++) {

    for (int j = 0; j < L*L; j++) {

      int x = unif_dist(e2);
      int y = unif_dist(e2);
      int spin = arr[x][y]; 
      int flip_spin = -spin;
      
      //calculate change in energy based on boundary conditions
      double delta_E;

      if (x == 0 && y == 0)
        delta_E = -J * (flip_spin - spin) * (arr[x][y + 1] + arr[x + 1][y] + arr[x][L - 1] + arr[L - 1][y]);
        
      else if (x == 0 && y == L - 1)
          delta_E = -J * (flip_spin - spin) * (arr[x][y - 1] + arr[x + 1][y] + arr[L - 1][L - 1] + arr[0][0]);
          
      else if (x == L - 1 && y == 0)
          delta_E = -J * (flip_spin - spin) * (arr[x][y + 1] + arr[x - 1][y] + arr[L - 1][L - 1] + arr[0][0]);
          
      else if (x == L - 1 && y == L - 1)
          delta_E = -1 * J * (flip_spin - spin) * (arr[x][y - 1] + arr[x - 1][y] + arr[L - 1][0] + arr[0][L - 1]);
          
      else if (x == 0)
          delta_E = -J * (flip_spin - spin) * (arr[x][y - 1] + arr[x][y + 1] + arr[x + 1][y] + arr[L - 1][y]);
          
      else if (x == L - 1)
          delta_E = -J * (flip_spin - spin) * (arr[x][y - 1] + arr[x][y + 1] + arr[x - 1][y] + arr[0][y]);
          
      else if (y == 0)
          delta_E = -J * (flip_spin - spin) * (arr[x][y + 1] + arr[x - 1][y] + arr[x + 1][y] + arr[x][L - 1]);
          
      else if (y == L - 1)
          delta_E = -J * (flip_spin - spin) * (arr[x][y - 1] + arr[x - 1][y] + arr[x + 1][y] + arr[x][0]);
          
      else
          delta_E = -J * (flip_spin - spin) * (arr[x][y - 1] + arr[x][y + 1] + arr[x - 1][y] + arr[x + 1][y]);
      
      //calculate flip prob and flip the spin if a randomly generated number is <= that flip_prob
      double flip_prob;

      if (delta_E == 8*J)
        flip_prob = min(1.0, E8);
        
      else if (delta_E == 6*J)
        flip_prob = min(1.0, E6);
        
      else if (delta_E == 4*J)
        flip_prob = min(1.0, E4);
      
      else if (delta_E == 2 * J)
        flip_prob = min(1.0, E2);
        
      else if (delta_E == 0)
        flip_prob = 1.0;
        
      else if (delta_E == -2 * J)
          flip_prob = min(1.0, Eneg2);
          
      else if (delta_E == -4 * J)
          flip_prob = min(1.0, Eneg4);
          
      else if (delta_E == -6 * J)
          flip_prob = min(1.0, Eneg6);
          
      else if (delta_E == -8 * J)
          flip_prob = min(1.0, Eneg8);

      double random_number = r_dist(e2);

      if (random_number <= flip_prob) {
      
        arr[x][y] = flip_spin; 
        
      }
      
    }
    
  }
  
}

//function for calculating energy
double Energy(vector<vector<int>>& arr, double J) { 

  double sum = 0;
  
  //lattice size
  int L = arr.size();

  for (int i = 0; i < (L - 1); i++) {

    for (int j = 0; j < (L - 1); j++) {

      sum += -J * arr[i][j] * (arr[i + 1][j] + arr[i][j + 1]);
      
    }
   
    sum += -J * arr[i][L - 1] * arr[i + 1][L - 1];
    sum += -J * arr[L - 1][i] * arr[L - 1][i + 1];
    
  }

  for (int i = 0; i < L; i++) {

    sum += -J * arr[i][0] * arr[i][L - 1];
    sum += -J * arr[0][i] * arr[L - 1][i];
    
  }

  return sum; 
  
}

//function for calculating magnetization, basically taking the sum of all spins at each lattice point
int Magnetization(vector<vector<int>>& arr) {
 
  int L = arr.size();
  int sum = 0;

  for (int i = 0; i < L; i++) {

    for (int j = 0; j < L; j++) {

      sum += arr[i][j];
      
    }
    
  }

  return sum;
  
}