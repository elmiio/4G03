#ifndef HV_H
#define HV_H

#include <vector> 
#include "matrix.h"

using namespace std;

#define BIT_SET(a,b) ((a) |= (1U<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1U<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1U<<(b)))
#define BIT_CHECK(a,b) ((bool)((a) & (1U<<(b))))
#define COND_BIT_SET(a,b,f) ((a) = ((a) & ~(1U<<(b))) | ((-(unsigned int)f) & (1U<<(b))))

void hv (Row &y, const Row &x, int L); 

#endif