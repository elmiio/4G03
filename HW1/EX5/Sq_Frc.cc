#include <iostream>
#include <iomanip>
#include <cmath>

// Function to compute the continued fraction approximation of sqrt(x)
double Sq_Frc(double x, int n) {
  
  double a = 1;
  
  //recurrence stops when n = 1
  if (n == 1) {
  
    return a + (x - a*a)/a;
  
  }
  
  //keep going if n != 1
  else {
  
    return a + (x - a*a)/(a + Sq_Frc(x, n - 1));
  
  }
    
}

int main() {
  
  //FOR TESTING SPEED OF PROGRAM --> UNCOMMENT THIS AND COMMENT OUT REST OF CODE TO TEST
  
  for (int i = 0; i < 10000000; i++) {
  
    double test = Sq_Frc(1.0, 10);
  
  }
  
  /*

  //set values
  double xstart = 0.5;
  double xend = 1.6;
  double stepsize = 0.1;
  int n = 10;

  std::cout << std::fixed << std::setprecision(6);

  std::cout << "x\t\tNewton Sqrt(x)\tActual Sqrt(x)\n";

  for (double x = xstart; x <= xend; x += stepsize) {
      
      std::cout << x << "\t" << Sq_Frc(x, n) << "\t" << std::sqrt(x) << std::endl;
                
  }
  
  */

  return 0;
  
}