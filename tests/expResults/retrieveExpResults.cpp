/**
 * @file retrieverExpResults.cpp
 * @author Víctor Loras Herrero
 * @brief Tests to check retrieval algorithms on experimental data.
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../../src/fourier.hpp"
#include "../../src/pulse.hpp"
#include "../../src/retrievers.hpp"

#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <fstream>
#include <numeric>

int main()
{
    //! Reading exp file values. Encapsulate in function better.
    // Open the CSV file
    std::ifstream paramsFile("./tests/expResults/axis_compressed_pulse.csv");

    if (!paramsFile.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

    // Define vectors to store data
    std::vector<double> omegas;
    std::vector<double> measuredDelays;

    std::string line;

    // Read data from the file line by line
    while (std::getline(paramsFile, line))
    {
        std::istringstream linestream(line);
        std::string value1_str, value2_str;

        // Read values from the line (assuming comma-separated)
        if (std::getline(linestream, value1_str, ',') && std::getline(linestream, value2_str, ','))
        {
            // Convert strings to doubles
            double value1 = std::stod(value1_str);
            double value2 = std::stod(value2_str);

            // Store values in vectors
            omegas.push_back(value1);
            measuredDelays.push_back(value2);
        }
        else
        {
            std::cerr << "Error reading line: " << line << std::endl;
        }
    }

    // Close the file
    paramsFile.close();

    //! Now determine the deltaT value

    // Calculate the mean of the differences
    double meanDifference;
    double sumDifference = 0.0;

    // Calculate the sum of differences
    for (size_t i = 1; i < measuredDelays.size(); ++i)
    {
        sumDifference += measuredDelays[i] - measuredDelays[i - 1];
    }

    // Calculate the mean difference
    meanDifference = sumDifference / (measuredDelays.size() - 1);

    int N = 128;
    double deltaT = meanDifference; //! Determined by the experimental values

    double t0 = 0;

    FourierTransform ft(N, deltaT, t0);

    std::vector<std::vector<double>> Tmeas; //! Fill this with the given trace by the data

    //! FILLING
    // Open the CSV file
    std::ifstream traceFile("./tests/expResults/compressed_pulse_900mA_trace.csv");

    if (!traceFile.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

    // Define matrix to store data
    std::vector<std::vector<double>> matrix;

    std::string line2;

    // Read data from the file line by line
    while (std::getline(traceFile, line2))
    {
        std::istringstream linestream(line2);
        std::string value_str;

        // Define a row for the matrix
        std::vector<double> row;

        // Read values from the line (assuming comma-separated)
        while (std::getline(linestream, value_str, ','))
        {
            // Convert string to double and add to the row
            double value = std::stod(value_str);
            row.push_back(value);
        }

        // Add the row to the matrix
        matrix.push_back(row);
    }

    // Close the file
    traceFile.close();

    //! END FILLING

    COPRA copraRetriever(ft, matrix, measuredDelays);
    Pulse retrievedPulse = copraRetriever.retrieve(1e-10, 20000);

    std::vector<std::vector<double>> Tret = retrievedPulse.getTrace();
    std::vector<std::complex<double>> retrievedField = retrievedPulse.getField();
    std::vector<std::complex<double>> retrievedSpectrum = retrievedPulse.getSpectrum();

    // Save the result

    FILE *g;
    g = fopen("trace_check_retrieved.txt", "wt");

    if (g == NULL)
    {
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(g, "%lf\t", Tret[i][j]);
        }
    }

    fclose(g);

    return 0;
}