//want to create function for recurrence relation in descending manner

#include <iostream>
#include <vector>
#include <utility>
#include <iomanip>

//function for calculating bessel functions in descending manner

std::pair<double, double> calc_bessel_descending(double x, int n, int N) {

  //create vector (list) to store J values
  std::vector<double> J(N + 1, 0.0);
  
  //initial J values (using eq. 7) to start recurrence
  J[N] = 1.0;
  J[N - 1] = 0.0;
  
  //recurrence relation
  for (int i = N - 1; i > 0; --i) {
  
    J[i - 1] = ((2.0 * i)/x) * J[i] - J[i + 1];
  
  }
  
  //normalization using eq. 8 in handout
  
  //first calcuulate what the summation term would be in eq. 8
  double sum = 0.0;
  
  for (int i = 2; i <= N; i += 2) {
  
    sum += J[i];
  
  }
  
  double scale_factor = 1/(J[0] + (2 * sum));
  
  //now go through every bessel function in the list and normalize them
  for (int i = 0; i <= N; ++i) {
  
    J[i] *= scale_factor;
  
  }
  
  return std::make_pair(J[0], J[1]);

}

int main() {

  //set x, n and N
  double x = 0.5;
  int n = 20;
  int N = 40;
  
  //use recurrence relation to calculate bessel functions
  std::pair<double, double> bessel = calc_bessel_descending(x, n, N);
  
  //output calculated values for J_0 and J_1
  std::cout << std::fixed << std::setprecision(15);
  std::cout << "J_0(0.5) = " << bessel.first << std::endl;
  std::cout << "J_1(0.5) = " << bessel.second << std::endl;

}