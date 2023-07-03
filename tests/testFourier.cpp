/**
 * @file testFourier.cpp
 * @author Víctor Loras Herrero
 * @brief Small tests to prove that the fft algorithms are working ok
 * @version 0.1
 * @date 2023-07-02
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "../src/utils.hpp"
#include "../src/fourier.hpp"

#include <iostream>
#include <vector>
#include <complex>


int main()
{
    int N = 65536;
    double signalDuration = 200;
    double deltaT = signalDuration / (N - 1);

    std::vector<double> t (N);

    for (int i = 0; i < N; i++)
    {
        t[i] =  - signalDuration / 2 + i * deltaT;
    }
    

    std::vector<double> omega = toAngularFrequency(fftFreq(N, deltaT));

    double deltaOmega = 2 * M_PI / signalDuration;


    std::vector<std::complex<double>> x(N);

    for (int i = 0; i < N; i++)
    {
        double value = t[i];
        double sincValue = (value != 0.0) ? std::sin(value) / value : 1.0;  // Compute sinc(x)
        double squareValue = (value >= -signalDuration/4 && value <= signalDuration/4) ? 1.0 : 0.0; // Compute square function

        x[i] = std::complex<double>(sincValue, squareValue);
    }

    std::vector<std::complex<double>> result (N);
    std::vector<std::complex<double>> ret_result(N);

    for (int i = 0; i < 3000; i++)
    {
        result = DFT(x, t, deltaT, omega, deltaOmega);

        ret_result = IDFT(result, t, deltaT, omega, deltaOmega);
    }
    
    // std::vector<std::complex<double>> result = DFT(x, t, deltaT, omega, deltaOmega);

    // std::vector<std::complex<double>> ret_result = IDFT(result, t, deltaT, omega, deltaOmega);

    // FILE *f;
    // f = fopen("transform_check.txt", "wt");

    // if (f == NULL){
    //     return 1;
    // }

    // for (int i = 0; i < N; i++)
    // {
    //     fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t[i], std::real(x[i]), std::imag(x[i]), std::real(result[i]), std::imag(result[i]), std::real(ret_result[i]), std::imag(ret_result[i]));
    // }
    
    // fclose(f);

    return 0;
}