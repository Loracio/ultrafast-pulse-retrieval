/**
 * @file testClassFourier.cpp
 * @author Víctor Loras Herrero
 * @brief Small tests to prove that the FourierTransform class is working ok
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
#include <chrono>

int main()
{
    int N = 65536;
    double signalDuration = 1000;
    double deltaT = signalDuration / N;

    double t0 = -signalDuration / 2;

    std::vector<double> t (N);

    for (int i = 0; i < N; i++)
    {
        t[i] =  t0 + i * deltaT;
    }
    
    std::vector<double> omega = toAngularFrequency(fftFreq(N, deltaT));

    double deltaOmega = 2 * M_PI / (N * deltaT);

    FourierTransform ft(N, deltaT, t0);

    std::vector<std::complex<double>> x(N);

    for (int i = 0; i < N; i++)
    {
        // double sincValue = (t[i] != 0.0) ? std::sin(t[i]) / t[i] : 1.0;  // Compute sinc(x)
        double expValue = std::exp(-t[i]*t[i]);
        // double squareValue = (t[i] >= -signalDuration/4 && t[i] <= signalDuration/4) ? 1.0 : 0.0; // Compute square function

        x[i] = std::complex<double>(expValue, t[i]);
    }

    std::vector<std::complex<double>> result (N);
    std::vector<std::complex<double>> ret_result(N);

    int numberOfTimes = 1000;

    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < numberOfTimes; i++)
    {
        result = ft.forwardTransform(x);

        ret_result = ft.backwardTransform(result);
    }

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Compute the elapsed time
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);

    // Print the elapsed time
    std::cout << "Elapsed time FFT: " << duration.count() / numberOfTimes << " microseconds" << std::endl;

    FILE *f;
    f = fopen("transform_check.txt", "wt");

    if (f == NULL){
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t[i], std::real(x[i]), std::imag(x[i]), std::real(result[i]), std::imag(result[i]), std::real(ret_result[i]), std::imag(ret_result[i]));
    }
    
    fclose(f);

    return 0;
}