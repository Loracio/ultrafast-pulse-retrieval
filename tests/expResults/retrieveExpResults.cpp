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

void readAxis(const std::string &axisFilename, std::vector<double> &omegas, std::vector<double> &delays);
std::vector<std::vector<double>> readTrace(const std::string &traceFilename);
double computeDeltaT(const std::vector<double> &delays);
void writeResults(Pulse &retrievedPulse, const std::string &resultingFieldFilename, const std::string &resultingSpectrumFilename, const std::string &resultingTraceFilename, const std::string &resultingErrorsFilename, const std::vector<double> &resultingErrors);

int main(int argc, char *argv[])
{
    if (argc != 13)
    {
        std::cerr << "Usage: " << argv[0] << " N axisFilename traceFilename retriever maximumIterations tolerance initialCandidateField initialCandidateSpectrum resultingFieldFilename resultingSpectrumFilename resultingTraceFilename resultingErrorsFilename" << std::endl;
        return 1;
    }

    // Parse command line arguments
    int N = std::stoi(argv[1]);

    const char *axisFilename = argv[2];
    std::vector<double> omegas;
    std::vector<double> delays;
    readAxis(axisFilename, omegas, delays);

    double deltaT = computeDeltaT(delays); // Determined by the experimental values
    double t0 = 0;
    FourierTransform ft(N, deltaT, t0);

    const char *traceFilename = argv[3];
    std::vector<std::vector<double>> Tmeas = readTrace(traceFilename);

    const char *retriever = argv[4];
    // Handle retriever
    std::string retrieverStr(retriever);

    const char *initialCandidateField = argv[7];
    std::string candidateFieldStr(initialCandidateField);
    if (candidateFieldStr != "none")
    {
        std::cerr << "This option is not yet implemented." << std::endl;
        return 1;
    }

    const char *initialCandidateSpectrum = argv[8];
    std::string candidateSpectrumStr(initialCandidateSpectrum);
    if (candidateSpectrumStr != "none")
    {
        std::cerr << "This option is not yet implemented." << std::endl;
        return 1;
    }

    int maximumIterations = std::stoi(argv[5]);
    double tolerance = std::stod(argv[6]);

    const char *resultingFieldFilename = argv[9];
    const char *resultingSpectrumFilename = argv[10];
    const char *resultingTraceFilename = argv[11];
    const char *resultingErrorsFilename = argv[12];

    std::cout << resultingErrorsFilename << std::endl;

    if (retrieverStr == "COPRA")
    {
        COPRA selectedRetriever(ft, Tmeas, delays);
        Pulse retrievedPulse = selectedRetriever.retrieve(tolerance, maximumIterations);
        // Save the result
        writeResults(retrievedPulse, resultingFieldFilename, resultingSpectrumFilename, resultingTraceFilename, resultingErrorsFilename, selectedRetriever.allTraceErrors);
    }
    else if (retrieverStr == "GPA")
    {
        GPA selectedRetriever(ft, Tmeas, delays);
        Pulse retrievedPulse = selectedRetriever.retrieve(tolerance, maximumIterations);
        // Save the result
        writeResults(retrievedPulse, resultingFieldFilename, resultingSpectrumFilename, resultingTraceFilename, resultingErrorsFilename, selectedRetriever.allTraceErrors);
    }
    else
    {
        std::cerr << "Invalid retriever. Choose between 'COPRA' and 'GPA'." << std::endl;
        return 1;
    }

    return 0;
}

void readAxis(const std::string &axisFilename, std::vector<double> &omegas, std::vector<double> &delays)
{
    // Open the CSV file
    std::ifstream paramsFile(axisFilename);

    if (!paramsFile.is_open())
    {
        std::cerr << "Error opening the axis file!" << std::endl;
        return;
    }

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
            delays.push_back(value2);
        }
        else
        {
            std::cerr << "Error reading axis file. Line: " << line << std::endl;
        }
    }

    // Close the file
    paramsFile.close();
}

std::vector<std::vector<double>> readTrace(const std::string &traceFilename)
{
    // Open the CSV file
    std::ifstream traceFile(traceFilename);

    if (!traceFile.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
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

    return matrix;
}

double computeDeltaT(const std::vector<double> &delays)
{
    double sumDifference = 0.0;

    // Calculate the sum of differences
    for (int i = 1; i < delays.size(); ++i)
    {
        sumDifference += delays[i] - delays[i - 1];
    }

    // Calculate the mean difference
    return sumDifference / (delays.size() - 1);
}

void writeResults(Pulse &retrievedPulse, const std::string &resultingFieldFilename, const std::string &resultingSpectrumFilename, const std::string &resultingTraceFilename, const std::string &resultingErrorsFilename, const std::vector<double> &retrievedErrors)
{

    // Write trace result into file
    std::ofstream traceOutput(resultingTraceFilename);

    if (!traceOutput.is_open())
    {
        std::cerr << "Error creating trace result file." << std::endl;
        return;
    }
    std::vector<std::vector<double>> retrievedTrace = retrievedPulse.getTrace();

    for (int i = 0; i < retrievedPulse.N; i++)
    {
        for (int j = 0; j < retrievedPulse.N - 1; j++)
        {
            traceOutput << retrievedTrace[i][j] << "\t";
        }

        if (i != retrievedPulse.N - 1)
        {
            traceOutput << retrievedTrace[i][retrievedPulse.N - 1] << "\n"; // Add a newline after each row
        }
        else
        {
            traceOutput << retrievedTrace[i][retrievedPulse.N - 1]; // Do not add newline in the last row
        }
    }

    traceOutput.close();

    // Write field result into file
    std::ofstream fieldOutput(resultingFieldFilename);

    if (!fieldOutput.is_open())
    {
        std::cerr << "Error creating field result file." << std::endl;
        return;
    }

    std::vector<std::complex<double>> retrievedField = retrievedPulse.getField();

    for (int i = 0; i < retrievedPulse.N - 1; i++)
    {
        fieldOutput << std::real(retrievedField[i]) << "\t" << std::imag(retrievedField[i]) << "\n";
    }

    fieldOutput << std::real(retrievedField[retrievedPulse.N - 1]) << "\t" << std::imag(retrievedField[retrievedPulse.N - 1]) ;

    fieldOutput.close();

    // Write spectrum result into file
    std::ofstream spectrumOutput(resultingSpectrumFilename);

    if (!spectrumOutput.is_open())
    {
        std::cerr << "Error creating spectrum result file." << std::endl;
        return;
    }

    std::vector<std::complex<double>> retrievedSpectrum = retrievedPulse.getSpectrum();

    for (int i = 0; i < retrievedPulse.N - 1; i++)
    {
        spectrumOutput << std::real(retrievedSpectrum[i]) << "\t" << std::imag(retrievedSpectrum[i]) << "\n";
    }

    spectrumOutput << std::real(retrievedSpectrum[retrievedPulse.N - 1]) << "\t" << std::imag(retrievedSpectrum[retrievedPulse.N - 1]) ;

    spectrumOutput.close();

    // Write errors result into file
    std::ofstream errorsOutput(resultingErrorsFilename);

    if (!errorsOutput.is_open())
    {
        std::cerr << "Error creating errors result file." << std::endl;
        return;
    }

    for (int i = 0; i < retrievedErrors.size() - 1; i++)
    {
        errorsOutput << retrievedErrors[i] << "\n";
    }

    errorsOutput << retrievedErrors[retrievedErrors.size() - 1] << "\n"; ;

    errorsOutput.close();
}
