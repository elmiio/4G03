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

    //compute norm using reduction
    #pragma omp parallel for reduction(+:norm)
    for (int i = 0; i < vec.size(); ++i) {
    
        norm += vec[i] * vec[i];
        
    }
    
    norm = sqrt(norm);

    //normalize vector
    #pragma omp parallel for
    for (int i = 0; i < vec.size(); ++i) {
    
        vec[i] /= norm;
        
    }

    return norm;
    
}

int main() {

    //matrix size
    const int L = 10;
    
    #pragma omp parallel
    {
    
      #pragma omp single
      cout << "Number of threads: " << omp_get_num_threads() << endl;
      
    }

    //define matrix A
    vector<vector<double>> A(L, vector<double>(L, 0));
    
    for (int i = 0; i < L; ++i) {
    
        A[i][i] = -2;
        A[i][(i + 1) % L] = 1;  //cyclic boundary
        A[i][(i - 1 + L) % L] = 1;  //cyclic boundary
        
    }

    //initialize random vector y
    vector<double> y(L), x(L);
    srand(time(0));
    
    for (int i = 0; i < L; ++i) y[i] = rand() / (double)RAND_MAX;

    //normalize y
    normalize(y);

    //power Method
    const int maxIter = 1000;
    double lambda = 0.0, prevLambda = 0.0, tolerance = 1e-6;

    for (int iter = 0; iter < maxIter; ++iter) {
    
        //matrix-vector multiplication
        matVecMult(A, y, x);

        //compute eigenvalue
        lambda = normalize(x) / normalize(y);

        //convergence check
        if (abs(lambda - prevLambda) < tolerance) break;

        // Update y and previous lambda
        y = x;
        prevLambda = lambda;
        
    }

    //output results
    cout << "Extremal Eigenvalue: " << lambda << endl;
    cout << "Eigenvector: ";
    
    for (double val : y) cout << val << " ";
    cout << endl;

    return 0;
    
}
