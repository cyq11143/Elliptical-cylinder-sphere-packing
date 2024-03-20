#include<iostream>
#include<random>
#include<math.h>
#include<fstream>
#include<time.h>
#include<string>
#include<stdlib.h>
#include"rand_function.h"
using namespace std;
extern const double pi;
double distance(double r1, double theta1, double z1, double r2, double theta2, double z2, double* L, double alpha)
{
    double a[3][3];
    double r=(*L)/2, p;
    a[0][0] = r2 * cos(theta2); a[0][1] = r2 * sin(theta2); a[0][2] = z2;
    a[1][0] = r2 * cos(theta2 + alpha); a[1][1] = r2 * sin(theta2 + alpha); a[1][2] = z2 + *L;
    a[2][0] = r2 * cos(theta2 - alpha); a[2][1] = r2 * sin(theta2 - alpha); a[2][2] = z2 - *L;
    for (int i = 0; i < 3; i++)
    {
        p = abs(z1 - a[i][2]);
        if (p <= (( * L) / 2))  //periodic boundary condition
        {
            r = sqrt(pow(r1 * cos(theta1) - a[i][0], 2) + pow(r1 * sin(theta1) - a[i][1], 2) + pow(z1 - a[i][2], 2));
        }
    }
    return r;
}
double distance2(double r, double theta, double a, double ratio)
{
    a = a / 2; //'a' here reprensents the short radius  
    double s[50], re = 0;
    int min; double x[50];
    if (theta >= 0 && theta < pi / 2)
    {
        s[0] = 0; x[0] = pow(r - a, 2); min = 0;
        for (int i = 1; i < 50; i++)
        {
            s[i] = i * pi / 100;
            x[i] = pow(r * cos(theta) - a * cos(s[i]), 2) + pow(r * sin(theta) - a * ratio * sin(s[i]), 2);
            if (x[i] < x[i - 1]) min = i;
        }
        x[min] = sqrt(x[min]);
        re = x[min];
    }
    if (theta >= pi / 2 && theta < pi)
    {
        s[0] = pi / 2; x[0] = pow(r - a * ratio, 2); min = 0;
        for (int i = 1; i < 50; i++)
        {
            s[i] = i * pi / 100 + pi / 2;
            x[i] = pow(r * cos(theta) - a * cos(s[i]), 2) + pow(r * sin(theta) - a * ratio * sin(s[i]), 2);
            if (x[i] < x[i - 1]) min = i;
        }
        x[min] = sqrt(x[min]);
        re = x[min];
    }
    if (theta >= pi && theta < 3 * pi / 2)
    {
        s[0] = pi; x[0] = pow(r - a, 2); min = 0;
        for (int i = 1; i < 50; i++)
        {
            s[i] = i * pi / 100 + pi;
            x[i] = pow(r * cos(theta) - a * cos(s[i]), 2) + pow(r * sin(theta) - a * ratio * sin(s[i]), 2);
            if (x[i] < x[i - 1]) min = i;
        }
        x[min] = sqrt(x[min]);
        re = x[min];
    }
    if (theta >= 3 * pi / 2 && theta <= 2 * pi)
    {
        s[0] = 3 * pi / 2; x[0] = pow(r - a * ratio, 2); min = 0;
        for (int i = 1; i < 50; i++)
        {
            s[i] = i * pi / 100 + 3 * pi / 2;
            x[i] = pow(r * cos(theta) - a * cos(s[i]), 2) + pow(r * sin(theta) - a * ratio * sin(s[i]), 2);
            if (x[i] < x[i - 1]) min = i;
        }
        x[min] = sqrt(x[min]);
        re = x[min];
    }
    return re;
}
double energy_one_one(double r1, double theta1, double z1, double r2, double theta2, double z2, double* L, double alpha, double k)
{
    double radius = 0.5;
    double E = 0.0;
    if (r1 == r2 && theta1 == theta2 && z1 == z2)
    {
        return 0.0;
    }
    double l = distance(r1, theta1, z1, r2, theta2, z2, L, alpha);
    if (l <= 2*radius)
    {
        E = 0.5 * k * pow(l - 2*radius, 2); 
        return E;
    }
    else
    {
        return 0.0;
    }
}
double energy_one_all(double** P, double r, double theta, double z, int N1, double a, double ratio, double* L, double alpha, double k)
{
    double radius = 0.5;
    double E = 0.0;
    for (int i = 0; i < N1; i++)
    {
        E = E + energy_one_one(r, theta, z, P[1][i], P[2][i], P[3][i], L, alpha, k);
    }
    E = E / 2;
    double sum = 0;
    double w; w = distance2(r, theta, a, ratio);
    if (w >= radius)
    {
        sum = 0;
    }
    else
    {
        sum = k * pow(radius - w, 2) / 2;  
    }
    E = E + sum;
    return E;
}
double energy_all_all(double** P, int  N1, double a, double ratio, double* L, double alpha, double k)
{
    double E = 0.0;
    for (int i = 0; i < N1; i++)
    {
        E = E + energy_one_all(P, P[1][i], P[2][i], P[3][i], N1, a, ratio, L, alpha, k);
    }
    return E;
}
void displace(double* dr, double d)
{
    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);*/
    //int sm;
    //srand(time(NULL));
     //init_genrand(rand());
    double a1 = 0; double a2 = 0; double a3 = 0;
    while (a1 == 0 && a2 == 0 && a3 == 0)
    {
        /*sm = rand();
        a1 = d * (double)(sm % (8700) / 8700.0);
        sm = rand();
        a2 = 2 * pi * (double)(sm % (8700) / 8700.0);
        sm = rand();
        a3 = 2 * d * ((double)(sm % (8700) / 8700.0) - 0.5);*/
        a1 = d * genrand_float32_full(); a2 = 2 * pi * genrand_float32_full(); a3 = 2 * d * (genrand_float32_full() - 0.5);
        if (a1 * a1 + a3 * a3 > pow(d, 2))
        {
            a1 = 0; a2 = 0; a3 = 0;
        }
    }
    dr[0] = a1; dr[1] = a2; dr[2] = a3;
}
void generate(double** P, double a, double ratio, double* L, int N1)
{
    a = a / 2;
    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);*/
    double a1, a2, a3;
    //srand(time(NULL));
    //init_genrand(rand());
    //int sm;
    //srand(time(NULL));
    for (int i = 0; i < N1; i++)
    {
        a1 = 0; a2 = 0;
        while (a1 == 0 && a2 == 0)
        {
            /*sm = rand();
            a1 = a * ratio * (double)(sm % (87) / 87.0);
            sm = rand();
            a2 = 2 * pi * (double)(sm % (87) / 87.0);*/
            a1 = a * ratio * genrand_float32_notone(); a2 = 2 * pi * genrand_float32_notone();
            if (pow(a1 * cos(a2), 2) + pow(a1 * sin(a2) / ratio, 2) > pow(a, 2))
            {
                a1 = 0;  a2 = 0;
            }
        }
        //sm = rand();
        //a3 = L * ((double)(sm % (87) / 87.0) - 0.5);
        a3 = (* L) * (genrand_float32_notone() - 0.5);
        P[1][i] = a1; P[2][i] = a2; P[3][i] = a3;
        //cout << P[1][i] << " " << P[2][i] << " " << P[3][i] << endl;
    }
}
void mcmove(double** P, double* dr, double* pos, double a, double ratio, double* L, int N1, double d, int Ni)
{
    a = a / 2;
    double r1; double theta1; double z1; double x1, y1, x2, y2;
    double r2 = 0; double theta2 = 0; double z2 = 0;
    do
    {
        displace(dr, d);
        r1 = P[1][Ni]; theta1 = P[2][Ni]; z1 = P[3][Ni];
        //cout << "初始" << r1 << " " << theta1 << " " << z1 << endl;
        x1 = r1 * cos(theta1); y1 = r1 * sin(theta1);
        x2 = x1 + dr[0] * cos(dr[1]); y2 = y1 + dr[0] * sin(dr[1]); z2 = z1 + dr[2];
        //cout << "初始" << x1 << " " << y1 << " " << z1 << endl;
        //cout << "移动" << x2 << " " << y2 << " " << z2 << endl;
        r2 = sqrt(pow(x2, 2) + pow(y2, 2));
        if (y2 >= 0) { theta2 = acos(x2 / r2); }
        else { theta2 = 2 * pi - acos(x2 / r2); }
    } while (pow(r2 * cos(theta2), 2) + pow(r2 * sin(theta2) / ratio, 2) > pow(a, 2));
    if (z2 > *L / 2) { z2 -= *L; }
    if (z2 < -*L / 2) { z2 += *L; }
    pos[0] = r2; pos[1] = theta2; pos[2] = z2;
}
void MCstep(double** P, double* dr, double* pos, double a, double ratio, double* L, double beta, int N1, double d, double* acc, double* Mc, double alpha, double k)//�Ե�Ni�����ӵ�Monte Carlo ��
{
    int Ni;
    double de = 0; double e1 = 0; double e2 = 0;
    double r1, theta1, z1;
    Ni = rand() % (N1);
    mcmove(P, dr, pos, a, ratio, L, N1, d, Ni);
    r1 = P[1][Ni]; theta1 = P[2][Ni]; z1 = P[3][Ni];
    e1 = energy_one_all(P, r1, theta1, z1, N1, a, ratio, L, alpha, k);
    P[1][Ni] = pos[0]; P[2][Ni] = pos[1]; P[3][Ni] = pos[2];
    e2 = energy_one_all(P, P[1][Ni], P[2][Ni], P[3][Ni], N1, a, ratio, L, alpha, k);
    de = e2 - e1;
    //srand(time(NULL));
    //if (exp(-de / T) > (double)(rand() % (8700) / 8700.0))//de < 0.0)//
    //cout << exp(-de * beta) <<" de="<<de<< endl;
    if (exp(-de * beta) > genrand_float32_notone())
    {
        *acc += 1;
    }
    else
    {
        *(*(P + 1) + Ni) = r1; *(*(P + 2) + Ni) = theta1; *(*(P + 3) + Ni) = z1;
    }
    *Mc += 1;
    //cout << "de=" << de << " " << "exp(-de/T)=" << exp(-de / T) << " " << "sk=" << sk <<"时间为:"<< clock() / CLOCKS_PER_SEC <<endl;
    return;
}
void mcvol(double** P, double a, double ratio, double* L, double beta, int N1, double dV, double alpha, double k, double p)
{
    double L1, V1, V0, E0, E1, arg, L0;
    L0 = *L;
    E0 = energy_all_all(P, N1, a, ratio, L, alpha, k);
    V0 = L0 * pi * a * a * ratio / 4;
    V1 = log(V0) + dV * (genrand_float32_notone() - 0.5);
    V1 = exp(V1); L1 = V1*4 / (pi * a * a * ratio);
    for (int i = 0; i < N1; i++) { P[3][i] = P[3][i] * L1 / L0; }
    *L = L1;
    E1 = energy_all_all(P, N1, a, ratio, L, alpha, k);
    arg = -beta * ((E1 - E0) + p * (V1 - V0) - ((double)N1 + 1.0) * log(L1 / L0) / beta);
    if (exp(arg) < genrand_float32_notone())
    {
        for (int i = 0; i < N1; i++)
        {
            P[3][i] = P[3][i] * L0 / L1;
        }
        *L = L0;
    }
    return;
}
void tot_MCstep(double** P, double* dr, double* pos, double a, double ratio, double* L, double beta, int N1, double dV, double d, double* acc, double* Mc, double alpha, double k, double p)
{
    double la;
    for (int i = 0; i <= N1; i++)
    {
        la = genrand_float32_notone() * N1 + 1; //cout << "la=" << la <<" " <<" N1=" << N1 <<" "<<"L="<<*L<<" ";
        if (la > (double)N1) { mcvol(P, a, ratio, L, beta, N1, dV, alpha, k, p); } //cout << "mcvol!" << endl; 
        else { MCstep(P, dr, pos, a, ratio, L, beta, N1, d, acc, Mc, alpha, k);  }//cout << "MCstep!" << endl;
    }
    return;
}
void ener_min(double** P, double* dr, double* pos, double a, double ratio, double* L, int N1, double d, double* acc, double* Mc, double alpha, double k)
{
    int Ni;
    double de = 0; double e1 = 0; double e2 = 0;
    double r1, theta1, z1;
    Ni = rand() % (N1);
    mcmove(P, dr, pos, a, ratio, L, N1, d, Ni);
    r1 = P[1][Ni]; theta1 = P[2][Ni]; z1 = P[3][Ni];
    e1 = energy_one_all(P, r1, theta1, z1, N1, a, ratio, L, alpha, k);
    P[1][Ni] = pos[0]; P[2][Ni] = pos[1]; P[3][Ni] = pos[2];
    e2 = energy_one_all(P, P[1][Ni], P[2][Ni], P[3][Ni], N1, a, ratio, L, alpha, k);
    de = e2 - e1;
    if (de < 0)
    {
        *acc += 1;
    }
    else
    {
        *(*(P + 1) + Ni) = r1; *(*(P + 2) + Ni) = theta1; *(*(P + 3) + Ni) = z1;
    }
    *Mc += 1;
    return;
}
