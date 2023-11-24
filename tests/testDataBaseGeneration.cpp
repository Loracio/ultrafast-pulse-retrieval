/**
 * @file testPulseClass.cpp
 * @author Víctor Loras Herrero
 * @brief Test database creation
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
#include <string>
#include <fstream>

int main()
{
    // Pulse parameters
    int N = 64;
    double signalDuration = 1;
    double deltaT = signalDuration / N;

    double t0 = 0;

    FourierTransform ft(N, deltaT, t0);

    Pulse generatedPulse(ft);
    std::vector<std::vector<double>> Tmn;
    std::vector<std::complex<double>> field;

    // DB parameters
    int numberOfPulses = 10000;

    double initialTBP = 0.51;
    double finalTBP = 0.81;
    double stepTBP = 0.10;

    int intervals = (finalTBP - initialTBP) / stepTBP;
    int pulsesPerTBP = numberOfPulses / intervals;

    double currentTBP = initialTBP;

    std::string filename = std::to_string(numberOfPulses) + "_randomPulses_N" + std::to_string(N) + ".csv";
    std::ofstream myFile(filename);

    std::cout << "Generating " << numberOfPulses << " random pulses... (" << pulsesPerTBP << " pulses per TBP)" << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < intervals; i++)
    {
        for (int j = 0; j < pulsesPerTBP; j++)
        {
            generatedPulse.randomPulse(currentTBP);
            field = generatedPulse.getField();

            //! Encapsulate writing in a function!!!!
            myFile << currentTBP << ",";

            // Real part of field
            for (int k = 0; k < N; ++k)
            {
                myFile << std::real(field[k]) << ",";
            }

            // Imaginary part of field
            for (int k = 0; k < N; ++k)
            {
                myFile << std::imag(field[k]) << ",";
            }

            // Trace
            Tmn = trace(field, ft.t, deltaT);
            for (int k = 0; k < N; ++k)
            {
                for (int l = 0; l < N; ++l)
                {
                    myFile << Tmn[k][l];
                    if (l == N - 1 && k == N - 1)
                    {
                        
                    }
                    else{
                        myFile << ","; // Write a comma if not at the last element
                    }
                }
            }

            myFile << "\n";
        }

        currentTBP += stepTBP;
    }

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Close the file
    myFile.close();

    // Compute the elapsed time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print the elapsed time
    std::cout << "Elapsed time generating database: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}