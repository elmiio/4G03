#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <omp.h>

using namespace std;

//function to multiply matrix A with vector y
void matVecMult(const vector<vector<double>>& A, const vector<double>& y, vector<double>& x) {

  int L = A.size();
  
  #pragma omp parallel for
  for (int i = 0; i < L; ++i) {
  
      x[i] = 0;
      for (int j = 0; j < L; ++j) {
      
          x[i] += A[i][j] * y[j];
          
      }
      
  }
    
}

//normalization function
double normalize(vector<double>& vec) {

  double norm = 0.0;

  #pragma omp parallel for reduction(+:norm)
  for (int i = 0; i < vec.size(); ++i) {
  
      norm += vec[i] * vec[i];
      
  }
  
  norm = sqrt(norm);

  #pragma omp parallel for
  for (int i = 0; i < vec.size(); ++i) {
  
      vec[i] /= norm;
      
  }

  return norm;
    
}

int main() {

    //matrix size
    const int L = 100;
    const int maxIter = 500;
    const double tolerance = 1e-6;

    //define matrix A
    vector<vector<double>> A(L, vector<double>(L, 0));
    for (int i = 0; i < L; ++i) {
        A[i][i] = -2;
        A[i][(i + 1) % L] = 1;
        A[i][(i - 1 + L) % L] = 1;
    }

    //initialize random vector y
    vector<double> y(L), x(L);
    srand(time(0));
    for (int i = 0; i < L; ++i) y[i] = rand() / (double)RAND_MAX;

    //normalize y
    normalize(y);

    //sequential timing
    double start = omp_get_wtime();

    double lambda = 0.0, prevLambda = 0.0;
    for (int iter = 0; iter < maxIter; ++iter) {
    
        matVecMult(A, y, x);
        lambda = normalize(x) / normalize(y);
        
        if (abs(lambda - prevLambda) < tolerance) break;
        y = x;
        prevLambda = lambda;
        
    }

    double end = omp_get_wtime();
    double exec = end - start;
    
    cout << "Sequential Execution Time: " << exec << " seconds\n";
    
    //parallel timing with OpenMP
    for (int threads = 2; threads <= 32; threads *= 2) {
    
      //set the number of threads
      omp_set_num_threads(threads);

      start = omp_get_wtime();

      lambda = 0.0, prevLambda = 0.0;
      
      for (int iter = 0; iter < maxIter; ++iter) {
      
          matVecMult(A, y, x);
          lambda = normalize(x) / normalize(y);
          
          if (abs(lambda - prevLambda) < tolerance) break;
          y = x;
          prevLambda = lambda;
          
      }

      end = omp_get_wtime();
      double speedup = exec/(end - start);
      cout << "Threads: " << threads
           << ", Parallel Execution Time: " << (end - start) << " seconds --> speedup by factor of " << speedup << endl;
           
    }

    return 0;
    
}
