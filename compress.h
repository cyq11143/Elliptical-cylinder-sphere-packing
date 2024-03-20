#ifndef compress_H
#define compress_H
#include<iostream>
#include<random>
#include<math.h>
#include<fstream>
#include<time.h>
#include<string>
#include<stdlib.h>
using namespace std;
extern const double pi;
double distance(double r1, double theta1, double z1, double r2, double theta2, double z2, double* L, double alpha);
double distance2(double r, double theta, double a, double ratio);
double energy_one_one(double r1, double theta1, double z1, double r2, double theta2, double z2, double* L, double alpha, double k);
double energy_one_all(double** P, double r, double theta, double z, int N1, double a, double ratio, double* L, double alpha, double k);
double energy_all_all(double** P, int  N1, double a, double ratio, double* L, double alpha, double k);
void displace(double* dr, double d);
void generate(double** P, double a, double ratio, double* L, int N1);
void mcmove(double** P, double* dr, double* pos, double a, double ratio, double* L, int N1, double d, int Ni);
void MCstep(double** P, double* dr, double* pos, double a, double ratio, double* L, double beta, int N1, double d, double* acc, double* Mc, double alpha, double k);
void mcvol(double** P, double a, double ratio, double* L, double beta, int N1, double dV, double alpha, double k, double p);
void tot_MCstep(double** P, double* dr, double* pos, double a, double ratio, double* L, double beta, int N1, double dV, double d, double* acc, double* Mc, double alpha, double k, double p);
void ener_min(double** P, double* dr, double* pos, double a, double ratio, double* L, int N1, double d, double* acc, double* Mc, double alpha, double k);
#endif
