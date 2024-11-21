#include <vector>
#include <cmath>
#include <iostream>
#include "matrix.h"
#include "hv.h"

#define BIT_SET(a,b) ((a) |= (1U<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1U<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1U<<(b)))
#define BIT_CHECK(a,b) ((bool)((a) & (1U<<(b))))
#define COND_BIT_SET(a,b,f) ((a) = ((a) & ~(1U<<(b))) | ((-(unsigned int)f) & (1U<<(b))))

using namespace std;

void hv(Row &y, const Row &x, int L)
{ 

  bool b ;
  unsigned int k;
  
  for (unsigned int i = 0; i < x.size(); i++) {
  
    if (abs(x[i]) > 2.2e-16) {
    
      int jm = L - 1;
      double xov2=x[i] / 2.0;
      for (int j=0 ; j<L; j++){
      
        k=i;
        COND_BIT_SET(k, jm, BIT_CHECK(i, j));
        COND_BIT_SET(k, j, BIT_CHECK(i, jm));
        y[k] += xov2;
        jm = j;
        
      }
      
    }
    
  }
  
  for (unsigned int i = 0; i < x.size(); i++) {
  
    y[i]=y[i]-((double) L)/2.0*x[i]/2.0;
  
  }


}