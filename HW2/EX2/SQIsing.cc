#include <cmath>
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

#include "MCvar.h"
#include "MCSweeps.h"

using namespace std;

//boltzmann constant
const double k = 1.380649e-23;

//function for writing lookup table
void writeTable(double T, double J) {

  //compute exponential values for different energy changes (E_delta)
  double E8 = exp(-8*J/k/T);
  double E6 = exp(-6*J/k/T);
  double E4 = exp(-4*J/k/T);
  double E2 = exp(-2*J/k/T);
  double Eneg2 = exp(2*J/k/T);
  double Eneg4 = exp(4*J/k/T);
  double Eneg6 = exp(6*J/k/T);
  double Eneg8 = exp(8*J/k/T);

  ofstream outFile("table.txt");

  outFile << E8 << endl;
  outFile << E6 << endl;
  outFile << E4 << endl;
  outFile << E2 << endl;
  outFile << Eneg2 << endl;
  outFile << Eneg4 << endl;
  outFile << Eneg6 << endl;
  outFile << Eneg8 << endl; 
  
  outFile.close();
        
}

int main() {

  //interaction energy
  double J = 1.0;
  
  //temp
  double T;
  
  cout << "Temperature: ";
  cin >> T;
  
  //run params
  int L = 16;
  int Nmeas = 100;
  int Nwarmup = 1000;
  int Nstep = 1;
  
  //create lookup table
  writeTable(T, J);
  
  //random number generation for initiializing spins and lattice
  random_device r;
  mt19937 e2(r());
  uniform_int_distribution<int> initRNG(0, 1);
  vector<vector<int>> SpinConf(L, vector<int>(L));
  
  //assign random spins for each position in lattice
  for (int i = 0; i < L; i++) {
  
    for (int j = 0; j < L; j++) {
    
      if (initRNG(e2) == 1) {
      
        SpinConf[i][j] = 1;
        
      } 
      
      else { 
      
        SpinConf[i][j] = -1;
        
      }
    
    }
  
  }
  
  //"warmup" sweeps
  MCSweeps(SpinConf, Nwarmup, J);
  
  //want energy and magnetization + second and fourth powers stored as MCvar objects so that we can compute square averages and stuff
  MCvar<double> E, E2;
  MCvar<double> M, M2, M4;
  
  ofstream M_file("M.txt");
  
  double e, m;
  
  //measurements in for loop
  for (int i = 0; i < Nmeas * 1000; i++) {
  
    MCSweeps(SpinConf, Nstep, J);
    e = Energy(SpinConf, J);
    m = Magnetization(SpinConf);
    
    M_file << m << endl;
    
    //store measurements in MCvar objects
    M.push(m);
    M2.push(pow(m, 2));
    M4.push(pow(m, 4));
    E.push(e);
    E2.push(pow(e, 2));
  
  }
  
  M_file.close();
  
  //specific heat capacity
  double C = L*L*k*k*T*T;
  double Cvkb = (E2.Avrg() - E.Avrg()*E.Avrg())/C;
  double Cvkb_err = E2.Err_Avrg() / C - (2 * E.Avrg() * E.Err_Avrg()) / C;
  double g = (3 - (M4.Avrg() / (M2.Avrg() * M2.Avrg()))) / 2;
  
  //output
  cout << "----------------------------------------------------------" << endl;
  cout << "Average energy: " << E.Avrg() << " +- " << E.Err_Avrg() << endl;
  cout << "Average magnetization: " << M.Avrg() << " +- " << M.Err_Avrg() << endl;
  cout << "Chi: " << M2.Avrg() / (T * L * L) << " +- " << M2.Err_Avrg() / (T * L * L) << endl;
  cout << "Cv/kb: " << Cvkb << " +- " << Cvkb_err << endl;
  cout << "g: " << g << endl;
  
  return 0;

}