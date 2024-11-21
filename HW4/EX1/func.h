#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

double V(double x);

double V2(double x);

double k2(double x, double E);

double phi_l(double x, double h, double current, double prev, double E);

double phi_r(double x, double h, double current, double prev, double E);

void save(const Row &array1, const Row &array2, const string filename);

void merge(Row& phi, Row& phil, Row& phir, Row& xl, Row& xr, double xm, double E, double Nsteps, double h);

void norm(Row& vec, double h);

#endif