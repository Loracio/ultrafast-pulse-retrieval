/**
 * @file testPulseClass.cpp
 * @author Víctor Loras Herrero
 * @brief Small tests to prove that the Pulse class is working ok
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "../src/utils.hpp"
#include "../src/fourier.hpp"
#include "../src/pulse.hpp"

#include <iostream>
#include <vector>
#include <complex>
#include <chrono>

int main()
{
    int N = 64;
    double signalDuration = 10;
    double deltaT = signalDuration / N;

    double t0 = 0;

    std::vector<double> t (N);

    for (int i = 0; i < N; i++)
    {
        t[i] =  t0 + i * deltaT;
    }
    
    std::vector<double> omega = toAngularFrequency(fftFreq(N, deltaT));

    FourierTransform ft(N, deltaT, t0);

    Pulse examplePulse(ft);


    int numberOfTimes = 100;

    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < numberOfTimes; i++)
    {
        examplePulse.randomPulse(0.7);
    }

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Compute the elapsed time
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);

    // Print the elapsed time
    std::cout << "Elapsed time generating pulses: " << duration.count() / numberOfTimes << " microseconds" << std::endl;

    std::vector<std::vector<double>> Tmn = trace(examplePulse.getField(), t, deltaT);


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

    return 0;
}