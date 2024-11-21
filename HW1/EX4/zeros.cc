#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

//minimum possible interval size (i.e. b - a value at which to stop bisection process)
double min = 1e-8;

//function for J_0(x)
double J0(double x) {

  return j0(x);

}

//function for executing bisection method over given interval [a, b]
double bisection(double a, double b) {

  double mid;
  
  while ((b - a) > min) {
  
    mid = (a + b)/2.0;
  
    if (J0(a) * J0(mid) < 0) {
    
      b = mid;
    
    }
    
    else {
    
      a = mid;
    
    }
  
  }
  
  return mid;

}

int main() {

  //initial guesses for the intervals of the first three zeros (using fig.1 to estimate)
  double intervals[3][2] = {
  
  {2.0, 3.0},
  {5.0, 6.0},
  {8.0, 9.0}
  
  };
  
  //execute bisection function on the three defined intervals and print out results
  for (int i = 0; i < 3; ++i) {
  
    double zero = bisection(intervals[i][0], intervals[i][1]);
    std::cout << "zero " << i + 1 << " = " << zero << std::endl;
  
  }
  
  return 0;

}