//want to calculate J_20(0.5) using a recurrence relation in ascending manner

#include <iostream>
#include <vector>

//function for calculating J_20 in ascending manner
double calc_bessel_ascending(double x, double J_0, double J_1, int n) {

  //initialize vector (list) to store J_n's in and add J_0, J_1
  std::vector<double> J;
  J.push_back(J_0);
  J.push_back(J_1);
  
  //calculate J_n's using recurrence relation in for loop and add them to J vector
  for (int i = 1; i < n; ++i) {
  
    double J_next = ((2.0 * i)/x) * J[i] - J[i - 1];
    J.push_back(J_next);
  
  }
  
  return J[n];

}

int main() {

  //given information
  double x = 0.5;
  double J_0 = 0.938469807240813;
  double J_1 = 0.2422684577;
  int n = 20;
  
  double J_20 = calc_bessel_ascending(x, J_0, J_1, n);
  
  //output msg
  std::cout << "J_20(" << x << ") = " << J_20 << std::endl;
  
  return 0;

}