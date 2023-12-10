/**
 * @file retrieverExpResults.cpp
 * @author Víctor Loras Herrero
 * @brief Tests to check retrieval algorithms on randomly generated pulses.
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../src/fourier.hpp"
#include "../src/pulse.hpp"
#include "../src/retrievers.hpp"

#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <fstream>
#include <numeric>

void writeOriginalPulse(Pulse &originalPulse, const std::string &originalFieldFilename, const std::string &originalSpectrumFilename, const std::string &originalTraceFilename, const std::vector<std::vector<double>> &originalTrace);
void writeResults(Pulse &retrievedPulse, const std::string &resultingFieldFilename, const std::string &resultingSpectrumFilename, const std::string &resultingTraceFilename, const std::string &resultingErrorsFilename, const std::vector<double> &resultingErrors);

int main(int argc, char *argv[])
{
    if (argc != 17)
    {
        std::cerr << "Usage: " << argv[0] << " N deltaT TBP noiseLevel originalFieldFilename originalSpectrumFilename originalTraceFilename retriever maximumIterations tolerance initialCandidateField initialCandidateSpectrum resultingFieldFilename resultingSpectrumFilename resultingTraceFilename resultingErrorsFilename" << std::endl;
        return 1;
    }

    // Parse command line arguments
    int N = std::stoi(argv[1]);

    double deltaT = std::stod(argv[2]);
    double TBP = std::stod(argv[3]);
    double noiseLevel = std::stod(argv[4]);

    const char *originalFieldFilename = argv[5];
    const char *originalSpectrumFilename = argv[6];
    const char *originalTraceFilename = argv[7];

    double t0 = 0;
    FourierTransform ft(N, deltaT, t0);

    Pulse originalPulse(ft);
    originalPulse.randomPulse(TBP);
    std::vector<std::vector<double>> Tmeas = add_noise(originalPulse.getTrace(), N, noiseLevel);

    const char *retriever = argv[8];
    // Handle retriever
    std::string retrieverStr(retriever);

    int maximumIterations = std::stoi(argv[9]);
    double tolerance = std::stod(argv[10]);

    const char *initialCandidateField = argv[11];
    std::string candidateFieldStr(initialCandidateField);
    if (candidateFieldStr != "none")
    {
        std::cerr << "This option is not yet implemented." << std::endl;
        return 1;
    }

    const char *initialCandidateSpectrum = argv[12];
    std::string candidateSpectrumStr(initialCandidateSpectrum);
    if (candidateSpectrumStr != "none")
    {
        std::cerr << "This option is not yet implemented." << std::endl;
        return 1;
    }

    const char *resultingFieldFilename = argv[13];
    const char *resultingSpectrumFilename = argv[14];
    const char *resultingTraceFilename = argv[15];
    const char *resultingErrorsFilename = argv[16];

    std::cout << resultingErrorsFilename << std::endl;

    if (retrieverStr == "COPRA")
    {
        COPRA selectedRetriever(ft, Tmeas);
        Pulse retrievedPulse = selectedRetriever.retrieve(tolerance, maximumIterations);
        // Save the result
        writeOriginalPulse(originalPulse, originalFieldFilename, originalSpectrumFilename, originalTraceFilename, Tmeas);
        writeResults(retrievedPulse, resultingFieldFilename, resultingSpectrumFilename, resultingTraceFilename, resultingErrorsFilename, selectedRetriever.allTraceErrors);
    }
    else if (retrieverStr == "GPA")
    {
        GPA selectedRetriever(ft, Tmeas);
        Pulse retrievedPulse = selectedRetriever.retrieve(tolerance, maximumIterations);
        // Save the result
        writeOriginalPulse(originalPulse, originalFieldFilename, originalSpectrumFilename, originalTraceFilename, Tmeas);
        writeResults(retrievedPulse, resultingFieldFilename, resultingSpectrumFilename, resultingTraceFilename, resultingErrorsFilename, selectedRetriever.allTraceErrors);
    }
    else if (retrieverStr == "PIE")
    {
        PIE selectedRetriever(ft, Tmeas);
        Pulse retrievedPulse = selectedRetriever.retrieve(tolerance, maximumIterations);
        // Save the result
        writeOriginalPulse(originalPulse, originalFieldFilename, originalSpectrumFilename, originalTraceFilename, Tmeas);
        writeResults(retrievedPulse, resultingFieldFilename, resultingSpectrumFilename, resultingTraceFilename, resultingErrorsFilename, selectedRetriever.allTraceErrors);
    }
    else
    {
        std::cerr << "Invalid retriever. Choose between 'COPRA', 'GPA' or 'PCGPA." << std::endl;
        return 1;
    }

    return 0;
}

void writeOriginalPulse(Pulse &originalPulse, const std::string &originalFieldFilename, const std::string &originalSpectrumFilename, const std::string &originalTraceFilename, const std::vector<std::vector<double>> &originalTrace)
{
    // Write trace result into file
    std::ofstream traceOutput(originalTraceFilename);

    if (!traceOutput.is_open())
    {
        std::cerr << "Error creating original trace file." << std::endl;
        return;
    }

    for (int i = 0; i < originalPulse.N; i++)
    {
        for (int j = 0; j < originalPulse.N - 1; j++)
        {
            traceOutput << originalTrace[i][j] << "\t";
        }

        if (i != originalPulse.N - 1)
        {
            traceOutput << originalTrace[i][originalPulse.N - 1] << "\n"; // Add a newline after each row
        }
        else
        {
            traceOutput << originalTrace[i][originalPulse.N - 1]; // Do not add newline in the last row
        }
    }

    traceOutput.close();

    // Write field result into file
    std::ofstream fieldOutput(originalFieldFilename);

    if (!fieldOutput.is_open())
    {
        std::cerr << "Error creating original field file." << std::endl;
        return;
    }

    std::vector<std::complex<double>> originalField = originalPulse.getField();

    for (int i = 0; i < originalPulse.N - 1; i++)
    {
        fieldOutput << std::real(originalField[i]) << "\t" << std::imag(originalField[i]) << "\n";
    }

    fieldOutput << std::real(originalField[originalPulse.N - 1]) << "\t" << std::imag(originalField[originalPulse.N - 1]);

    fieldOutput.close();

    // Write spectrum result into file
    std::ofstream spectrumOutput(originalSpectrumFilename);

    if (!spectrumOutput.is_open())
    {
        std::cerr << "Error creating original spectrum file." << std::endl;
        return;
    }

    std::vector<std::complex<double>> originalSpectrum = originalPulse.getSpectrum();

    for (int i = 0; i < originalPulse.N - 1; i++)
    {
        spectrumOutput << std::real(originalSpectrum[i]) << "\t" << std::imag(originalSpectrum[i]) << "\n";
    }

    spectrumOutput << std::real(originalSpectrum[originalPulse.N - 1]) << "\t" << std::imag(originalSpectrum[originalPulse.N - 1]);

    spectrumOutput.close();
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

    fieldOutput << std::real(retrievedField[retrievedPulse.N - 1]) << "\t" << std::imag(retrievedField[retrievedPulse.N - 1]);

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

    spectrumOutput << std::real(retrievedSpectrum[retrievedPulse.N - 1]) << "\t" << std::imag(retrievedSpectrum[retrievedPulse.N - 1]);

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

    errorsOutput << retrievedErrors[retrievedErrors.size() - 1] << "\n";
    ;

    errorsOutput.close();
}
