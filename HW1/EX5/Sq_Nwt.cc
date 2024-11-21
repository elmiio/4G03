#include <iostream>
#include <cmath>
#include <iomanip>

double Sq_Nwt(double x, double tolerance) {

  //set intial guess
  double guess = x;
  
  //will keep guessing until convergence is reached (convergence determined by the tolerance)
  while (true) {
  
    double nextGuess = (guess + x / guess) / 2;
    
    if (std::abs(nextGuess - guess) < tolerance) {
    
      return nextGuess;
        
    }
    
    guess = nextGuess;
      
  }
    
}

int main() {
  
  //FOR TESTING SPEED OF PROGRAM --> UNCOMMENT THIS AND COMMENT OUT REST OF CODE TO TEST
  
  for (int i = 0; i < 10000000; i++) {
  
    double test = Sq_Nwt(0.5, 1e-7);
  
  }
  
  /*

  //values to use
  
  double start = 0.5;
  double end = 1.6;
  double stepsize = 0.1;
  double tolerance = 1e-7; //set tolerance to this since we're using precision 6

  std::cout << std::fixed << std::setprecision(6);
  
  std::cout << "x\t\tNewton Sqrt(x)\tActual Sqrt(x)\n";
  
  for (double x = start; x <= end; x += stepsize) {
  
    std::cout << x << "\t" << Sq_Nwt(x, tolerance) << "\t" << std::sqrt(x) << "\n";
      
  }
  
  */

  return 0;
    
}