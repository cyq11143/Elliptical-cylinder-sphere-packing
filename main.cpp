#include<iostream>
#include<random>
#include<math.h>
#include<fstream>
#include<time.h>
#include<string>
#include<stdlib.h>
#include"rand_function.h"
#include"compress.h"
#include"output.h"
extern const double pi = 3.1415926;
int main()
{    
    srand(time(NULL));
    init_genrand(rand());
    double L0 = 8, H = 0, a, ratio, d, * L, T, E, k = 1, den = 0.0, T0 = 0.01, Te = 0.000001, dV = 0.05, beta, alpha = 0.0, p = 0.007, radius = 0.5;
    int N1=7,Nstart=35,Nend=35, Nloop=1,step=10000, NTemp = 100;
    double dr[3] = { 0,0,0 }, pos[3] = { 0,0,0 }; //a is the short diameter,b is the long diameter, radius is sphere's radius.

    /*cout << "ratio="; cin >> ratio;
    cout << "a=" ; cin >> a;
    cout << "N_start="; cin >> Nstart;
    cout << "N_end="; cin >> Nend;
    cout << "N_loop="; cin >> Nloop;*/
    Nstart = 7;
    Nend = 20;
    Nloop = 10;
    ratio = 1.00;
    a = 1.80;
    for (int sk = Nstart; sk <= Nend; sk++)
    {
        for (int zx = 0; zx < Nloop; zx++)
        {   
            double* er, * acc, * Mc, * P[4];
            d = 0.01;  E = 0.0; L = &L0; N1 = sk;
            dr[0] = 0; dr[1] = 0; dr[2] = 0; pos[0] = 0; pos[1] = 0; pos[2] = 0;
            acc = new double(0); er = new double(0); Mc = new double(0);
            P[0] = new double[N1]; P[1] = new double[N1]; P[2] = new double[N1]; P[3] = new double[N1];
            for (int i = 0; i < N1; i++)
            {
                P[0][i] = (double)i + 1; P[1][i] = 0; P[2][i] = 0; P[3][i] = 0;
            }
            /*cout << "圆柱短轴长a与比例ratio："; cin >> a; cin >> ratio;
            cout << "初始温度："; cin >> T0;
            cout << "降温组数："; cin >> NTemp;
            cout << "球数量："; cin >> N1;
            cout << "初始高度："; cin >> L0;
            cout << "旋转角alpha:"; cin >> alpha;
            cout << "系统压强:"; cin >> p;
            cout << "dV:"; cin >> dV;
            cout << "step:"; cin >> step;*/
            generate(P, a, ratio, L, N1);
            for (int i = 0; i <= NTemp; i++)
            {
                if (i < NTemp)
                {
                    T = T0 - (T0 - Te) * i / (NTemp - 1);
                    T = T0 * pow(T0 / Te, (T - T0) / (T0 - Te));
                    beta = 1 / T;
                    for (int j = 1; j <= step; j++) 
                    {
                        tot_MCstep(P, dr, pos, a, ratio, L, beta, N1, dV, d, acc, Mc, alpha, k, p);
                        if (*Mc > 1000) //recorded MCsteps
                        {
                            den = (*acc) / (*Mc);
                            if (den < 0.4) { d = d / 2; }
                            if (den > 0.6) { d = d * 2; }
                            *Mc = 0; *acc = 0;
                        }
                        if (j % 2000 == 0) 
                        {
                            cout << "E=" << energy_all_all(P, N1, a, ratio, L, alpha, k) << "  T=" << T <<"  N="<< N1 << " ratio="<<ratio<<" a="<<a<<endl;
                            H = energy_all_all(P, N1, a, ratio, L, alpha, k) + p * (*L) * pi * a * a * ratio / 4;
                            cout << "Enthalpy=" << H << "  L=" << *L << " time=" << clock() / CLOCKS_PER_SEC << endl;
                            cout << endl;
                        }
                    }
                }
                if (i == NTemp)
                {
                    double enthalpy = 0;
                    for (int j = 0; j < step; j++)
                    {
                        for (int i = 0; i < N1; i++)
                        {
                            ener_min(P, dr, pos, a, ratio, L, N1, d, acc, Mc, alpha, k);
                        }
                        if (j > (step - 40)) {
                            output_dump(P, step, N1, a, ratio, L, alpha);
                            enthalpy += (energy_all_all(P, N1, a, ratio, L, alpha, k) + p * pi * a * a * ratio * (*L) / 4);
                        }
                    }
                    enthalpy = enthalpy / 39;
                    output_csv(P, N1, a, ratio, L, alpha, k, p, NTemp, Te, enthalpy);
                }
                //record(P, N1, a, ratio, L, alpha, k, p, (i + 1) * step);
            }
            for (int i = 0; i < 4; i++)
                delete[]P[i];
            delete acc, er, Mc;
        }
    }
    return 0;
}
