#include<iostream>
#include<random>
#include<math.h>
#include<fstream>
#include<time.h>
#include<string>
#include<stdlib.h>
#include"compress.h"
#include<sstream>
extern const double pi;
using namespace std;
void output_dump(double** P, int step, int N1, double a, double ratio, double* L, double alpha)
{
    stringstream s1, s2, s3, s4, s5;
    s1 << a; s2 << ratio; s3 << (*L); s4 << N1; s5 << alpha;
    ofstream ofs;
    ofs.open("ratio=" + s2.str() + "a=" + s1.str() + "N =" + s4.str() + "L = " + s3.str() + ".dump", ios::out | ios::app);
    ofs << "ITEM: TIMESTEP" << endl << step << endl;
    ofs << "ITEM: NUMBER OF ATOMS" << endl << 3 * N1 << endl;
    ofs << "ITEM: BOX BOUNDS ff ff pp" << endl;
    ofs << -a/2 << ' ' << a/2 << endl
        << -a * ratio/2 << ' ' << a * ratio/2 << endl
        << -( * L) *3/ 2 << ' ' << ( * L) *3/ 2 << endl;
    ofs << "ITEM: ATOMS id x y z " << endl;
    for (int i = 0; i < N1; i++)
    {
        ofs << P[0][i] << " " << P[1][i] * cos(P[2][i]) << " " << P[1][i] * sin(P[2][i]) << " " << P[3][i] << endl;
    }
    for (int i = 0; i < N1; i++)
    {
        ofs << P[0][i]+N1 << " " << P[1][i] * cos(P[2][i]+alpha) << " " << P[1][i] * sin(P[2][i]+alpha) << " " << P[3][i] + *L << endl;
    }
    for (int i = 0; i < N1; i++)
    {
        ofs << P[0][i]+2*N1 << " " << P[1][i] * cos(P[2][i]-alpha) << " " << P[1][i] * sin(P[2][i] - alpha) << " " << P[3][i] - *L << endl;
    }
    ofs.close();
}
bool isFileExists_ifstream(string& name) 
{
    ifstream f(name.c_str());
    return f.good();
}
void output_csv(double** P, int N1, double a, double ratio, double* L, double alpha, double k, double p, int Ntemp, double Te,double enthalpy)
{
    // double H = energy_all_all(P, N1, a, ratio, L, alpha, k) + p * pi * a * a * ratio * (*L)/4;
    stringstream s1, s2, s3, s4, s5, s6, s7,s8,s9,s10;
    s1 << a; s2 << ratio; s3 << p; s4 << N1; s5 << alpha; s6 << (*L); s7 << (enthalpy / N1); s8 << energy_all_all(P, N1, a, ratio, L, alpha, k);
    s9 << Ntemp; s10 << Te;
    string str = "ratio=" + s2.str() + "a=" + s1.str() + "p=" + s3.str() + ".csv";
    bool ret = isFileExists_ifstream(str);
    ofstream ofs;
    ofs.open(str, ios::out | ios::app);
    if (ret) 
    {
        ;
    }
    else 
    {
        ofs << "N=" << ',' << "H=" << "," << "L=" << "," << "E=" << "," << "Ntemp=" << ',' << "Te=" << '\n';
    }
    ofs << s4.str() << ',' << s7.str() << ',' << s6.str() << ',' << s8.str() << ',' << s9.str() << ',' << s10.str() << '\n';
    ofs.close();
}
void record(double** P, int N1, double a, double ratio, double* L, double alpha, double k, double p, int step)
{
    double H = energy_all_all(P, N1, a, ratio, L, alpha, k) + p * pi * a * a * ratio * (*L);
    stringstream s1, s2, s3, s4, s5, s6, s7;
    s1 << a; s2 << ratio; s3 << p; s4 << alpha; s5 << N1; s6 << step; s7 << H;
    ofstream ofs;
    ofs.open("Record_a=" + s1.str() + "ratio=" + s2.str() + "p=" + s3.str() + "N=" + s5.str() + ".csv", ios::out | ios::app);
    //ofs << "step" << "," << "Enthalpy" << endl;
    ofs << "step" << "," << "Enthalpy" << endl;
    ofs << s6.str() << "," << s7.str() << endl;
    ofs.close();
}
