// Numerical Integration
#include <iostream>
#include <math.h>

using namespace std;

double myfunc (double x) { 
  
  return 4.0/(1 + (x*x));

}

double Trapezoidal(double (*f)(double x), double xl, double xh, int n) {

  double sum = 0;

  double delta = (xh - xl)/n; //delta x
  
  //lower and upper bounds
  double lower;
  double upper;
  
  for (int i = 0; i < n; i++) {
  
    lower = xl + delta*i;
    upper = xl + delta*(i + 1);
  
    sum += 0.5 * (f(upper) + f(lower)) * delta;
  
  }
  
  return sum;

}
double Romberg(double (*f)(double x), double xl, double xh, int n) {

  //array for romberg vals
  double R[n + 1][n + 1];
  double h;
  double sum;
  
  //initial romberg val
  R[0][0] = (f(xl) + f(xh))/2;
  
  //loop over R[n][0]
  for (int i = 1; i <= n; i++) {
  
    sum = 0;
    h = (xh - xl)/pow(2, i);
    
    //implement myfunc (eq. 14)
    for (int j = 1; j <= pow(2, i - 1); j++) {
    
      sum += myfunc(xl + ((2*j) - 1)*h);
    
    }
    
    R[i][0] = R[i - 1][0]/2 + h*sum;
  
  }
  
  //loop over R[n][m]
  for (int i = 1; i <= n; i++) {
  
    for (int j = 1; j <= n; j++) {
    
      R[j][i] = R[j][i - 1] + 1/(pow(4, i) - 1) * (R[j][i - 1] - R[j - 1][i - 1]);
    
    }
  
  }
  
  return R[n][n];

}

int main () {
  
  /*
  
  cout.precision(24);
  
  double integral;
  
  //number of steps
  int ntrap;
  int nromb;
  
  cout << "# of iterations for Trapezoidal integration: " << endl;
  cin >> ntrap;
  
  cout << "# of iterations for Romberg integration: " << endl;
  cin >> nromb;
  
  cout << "The Trapezoidal result is " << Trapezoidal(&myfunc, 0.0, 1.0, ntrap) << endl;
  cout << "The Romberg result is " << Romberg(&myfunc, 0.0, 1.0, nromb) << endl;
  
  */
  
  cout.precision(24);
  
  clock_t start = clock();
  
  for (int i = 0; i < 1000; i++) {
  
    double test = Trapezoidal(&myfunc, 0.0, 1.0, 2097152);
    
  }

  clock_t end = clock(); 
  
  double t_trap = static_cast<double>(end - start)/CLOCKS_PER_SEC;
  
  cout << "Trapezoidal integration done in " << t_trap << " seconds." << endl;
  
  start = clock();
  
  for (int i = 0; i < 1000; i++) {
  
    double test = Romberg(&myfunc, 0.0, 1.0, 6);
    
  }

  end = clock(); 
  
  double t_romb = static_cast<double>(end - start)/CLOCKS_PER_SEC;
  
  cout << "Romberg integration done in " << t_romb << " seconds." <<  endl;
  
  return 0;

}