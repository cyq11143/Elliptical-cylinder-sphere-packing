#ifndef output_H
#define output_H
#include<iostream>
#include<random>
#include<math.h>
#include<fstream>
#include<time.h>
#include<string>
#include<stdlib.h>
using namespace std;
void output_dump(double** P, int step, int N1, double a, double ratio, double* L, double alpha);
void output_csv(double** P, int N1, double a, double ratio, double* L, double alpha, double k, double p,int Ntemp,double Te,double enthalpy);
void record(double** P, int N1, double a, double ratio, double* L, double alpha, double k, double p, int step);
#endif
