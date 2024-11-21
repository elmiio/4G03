#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>

//function for taylor expansion of sqrt(x)
double Sq_Tlr(double x, int n) {

  double coeff = 1; //coefficent of term
  double exp = 0.5; //exponent for coefficient of term
  double a = 1.0;
  
  double sum = 1.0;
  
  for (int i = 1; i < n; i++) {
  
    coeff *= (0.5 - (i - 1))/i;
    sum += coeff * pow(x - a, i);
  
  }
  
  return sum;

}

int main() {
  
  //FOR TESTING SPEED OF PROGRAM --> UNCOMMENT THIS AND COMMENT OUT REST OF CODE TO TEST
  
  for (int i = 0; i < 10000000; i++) {
  
    double test = Sq_Tlr(1.0, 50);
  
  }
  
  /*
  
  //values to use
  double start = 0.5;
  double end = 1.6;
  double stepsize = 0.1;
  int n = 50;
    
  std::cout << std::fixed << std::setprecision(6);
  
  std::cout << "x\t\tTaylor Sqrt(x)\t\tActual Sqrt(x)" << std::endl;
  
  for (double x = start; x <= end; x += stepsize) {
    
    std::cout << x << "\t" << Sq_Tlr(x, n) << "\t\t" << std::sqrt(x) << std::endl;
  
  }
  
  */
  
  return 0;

}