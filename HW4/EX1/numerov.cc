#include <fstream>
#include <limits>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>

#include "func.h"
                   
using namespace std;                    

int main() {

  //energy & energy steps
  double Emax = 100;
  double Eh = 0.01;
  int Esteps = Emax/Eh;

  Row E (Esteps, 0.0);
  Row G (Esteps, 0.0);
  Row F (Esteps, 0.0);
  
  //populate energy array
  for (int i=0; i<Esteps; i++) {  
  
    E[i] = 0.01 + i*Eh; 
    
  }

  //run params
  double f = 0.0001;
  double max = M_PI/4;
  
  int Nsteps = 1001; 
  double h = 2*max / (Nsteps-1);
  int mid = max/h;
  
  Row xr (Nsteps, 0.0);
  Row xl (Nsteps, 0.0);
  Row phil (Nsteps, 0.0);
  Row phir (Nsteps, 0.0);
  Row phi (Nsteps, 0.0);

  //populate position arrays
  for (int i=0; i<Nsteps; i++) {
  
    xl[i] = -1 * max + i*h;
    xr[i] = max - i*h;
    
  }

  //initialize phi[0] and phi[1] for left and right
  phil[0] = 0;
  phir[0] = 0;
  
  phil[1] = f;
  phir[1] = f;

  //find eigenvalues
  for (int i=0; i<Esteps; i++) {

    //merge & normalize
    merge(phi, phil, phir, xl, xr, mid, E[i], Nsteps, h);
    norm(phi, h);

    G[i] = phil[mid+1] - phil[mid-1] - phir[mid-1] + phir[mid+1];
    F[i] = G[i] / phil[mid];

    if (i > 0) { 
    
      //even eigenvalues
      if (G[i] * G[i-1] < 0) { 
      
        cout << "Even eigenvalue: " << (E[i]+E[i-1])/2 << endl;
        
      } 
      
      //odd eigenvalues
      else if (F[i] * F[i-1] < 0) {
      
        cout << "Odd eigenvalue: " << (E[i]+E[i-1])/2 << endl;
        
      }  
      
    } 
    
  }

  //save E and G arrays to output files
  save(E, G, "G.txt");
  save(E, F, "F.txt");
  
  double Ewave; 
  cout << "Energy value to use for wavefunction plotting:" << endl;
  cin >> Ewave; 

  merge(phi, phil, phir, xl, xr, mid, Ewave, Nsteps, h);
  norm(phi, h);
   
  //save wavefunction to output file
  save(xl, phi, "wave.txt"); 
  
}