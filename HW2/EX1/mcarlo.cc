#include <iostream>
#include <random>
#include <cmath>

using namespace std;

//limits of integration
const double low = 0.0;
const double upp = 10.0;

//normalization factor
const double norm = (exp(upp - low) - 1) / exp(upp - low);

//function f: e^(-x^2)
double f(double x) {

    return exp(-x * x);
    
}

//function g: e^(-x) approximation to f
double g(double x) {

    return exp(-x) / norm;
    
}

//function G
double G(double x) {

    return (1 - exp(-x)) / norm;
    
}

//function for inverse of G
double G_inv(double u) {

    return -log(1 - u * norm);
    
}

int main() {

    random_device r;
    mt19937 e2(r());
    uniform_real_distribution<> my_r(0, 1);
  
    //# of iterations for MC integration
    int n = 10000000;
  
    //initialize integral sum
    double sum = 0;
  
    //define region of integration using G function and defined bounds
    double region = G(upp) - G(low);
  
    //integrate in for loop
    for (int i = 0; i < n; i++) {
        //generate random number from 0 to 1
        double u = my_r(e2);
    
        //transform random number with inverse G
        double x = G_inv(u);
    
        //add evaluated integrand to sum
        sum += f(x) / g(x);
    }
  
    //normalize the calculated sum and scale by the region that we integrated over
    sum = (sum/n)*region;
  
    //output
    cout << "Result: " << sum << endl;
  
    return 0;
    
}