/**
 * @file testTrace.cpp
 * @author Víctor Loras Herrero
 * @brief Small tests to prove that the trace function is working ok
 * @version 0.1
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "../src/utils.hpp"
#include "../src/fourier.hpp"
#include "../src/pulse.hpp"

#include <iostream>
#include <fstream>
#include <cstring> // Include the cstring header
#include <vector>
#include <complex>
#include <chrono>

int main(){
    int N = 128;
    double deltaT = 0.1;

    double t0 = 0;

    std::vector<double> t (N);

    for (int i = 0; i < N; i++)
    {
        t[i] =  t0 + i * deltaT;
    }
    
    std::vector<double> omega = toAngularFrequency(fftFreq(N, deltaT));

    std::complex<double> complexval[] = {
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, -0.0},
        {0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.000001, 0.000001},
        {-0.000003, 0.000003},
        {-0.000014, 0.000004},
        {-0.000060, -0.000001},
        {-0.000236, 0.000003},
        {-0.000805, 0.000143},
        {-0.002311, 0.000793},
        {-0.005600, 0.002514},
        {-0.011547, 0.005502},
        {-0.020289, 0.008940},
        {-0.029844, 0.011567},
        {-0.034833, 0.013615},
        {-0.027875, 0.017274},
        {-0.005324, 0.023705},
        {0.027180, 0.030040},
        {0.056932, 0.030678},
        {0.072948, 0.021028},
        {0.071898, -0.000014},
        {0.057340, -0.028375},
        {0.036369, -0.054317},
        {0.016395, -0.066871},
        {0.002265, -0.061673},
        {-0.004954, -0.044081},
        {-0.006856, -0.024144},
        {-0.005812, -0.009329},
        {-0.003799, -0.001558},
        {-0.002018, 0.001058},
        {-0.000903, 0.001219},
        {-0.000360, 0.000731},
        {-0.000138, 0.000323},
        {-0.000052, 0.000114},
        {-0.000018, 0.000033},
        {-0.000005, 0.000008},
        {-0.000001, 0.000002},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {0.0, 0.0},
        {0.0, -0.0},
        {0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, 0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0},
        {-0.0, -0.0}
    };
    
    std::vector<std::complex<double>> x(N);
    for (int i = 0; i < N; i++)
    {
        x[i] = complexval[i];
    }

    std::vector<std::vector<double>> Tmn = trace(x, t, deltaT);

    FourierTransform ft(N, deltaT, t0);
    Pulse testPulse(ft);
    testPulse.setField(x);
    std::vector<std::vector<double>> ret_trace = testPulse.getTrace();

    FILE *f;
    f = fopen("trace_check.txt", "wt");

    if (f == NULL){
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(f, "%lf\t", Tmn[i][j]);
        }
    }

    fclose(f);
    
    FILE *g;
    g = fopen("trace_check_class.txt", "wt");

    if (g == NULL){
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(g, "%lf\t", ret_trace[i][j]);
        }
    }

    fclose(g);

}