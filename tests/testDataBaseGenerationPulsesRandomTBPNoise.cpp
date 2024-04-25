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
#include "../src/retrievers.hpp"

#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <string>
#include <fstream>

#include "H5Cpp.h"

// Function to find index of maximum magnitude in vector
int max_magnitude_index(const std::vector<std::complex<double>> &v)
{
    return std::distance(v.begin(), std::max_element(v.begin(), v.end(), [](const std::complex<double> &a, const std::complex<double> &b)
                                                     { return std::abs(a) < std::abs(b); }));
}

int main()
{
    // Pulse parameters
    int N = 128;
    double signalDuration = 1;
    double deltaT = signalDuration / N;

    double t0 = 0;

    FourierTransform ft(N, deltaT, t0);

    Pulse generatedPulse(ft);
    std::vector<std::complex<double>> field;
    std::vector<std::vector<double>> originalTrace;
    std::vector<std::vector<double>> noisyTrace;

    // DB parameters
    int numberOfPulses = 10;

    double initialTBP = 0.51;
    double finalTBP = 1.575;

    double noiseLevel = 0.01;

    int maximumIterations = 1000; // Maximum number of iterations for the COPRA algorithm
    double tolerance = 1e-10;     // Tolerance for the COPRA algorithm

    double currentTBP; // It will be randomly generated uniformly between initialTBP and finalTBP for each pulse
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(initialTBP, finalTBP);

    std::string filename = std::to_string(numberOfPulses) + "_randomPulses_N" + std::to_string(N) + "with_" + std::to_string(noiseLevel) + "noise.h5";
    // Create the HDF5 file
    H5::H5File file(filename, H5F_ACC_TRUNC);
    // Check if the file was opened correctly
    if (!file.getId())
    {
        std::cerr << "Error opening file " << filename << std::endl;
        return -1;
    }

    std::cout << "Generating " << numberOfPulses << " random pulses..." << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    // Create two vectors: one for the error in each retrieval and one for the time it took to retrieve each pulse
    std::vector<double> errors(numberOfPulses);
    std::vector<double> retrievalTimes(numberOfPulses);

    for (int i = 0; i < numberOfPulses; i++)
    {
        currentTBP = dis(gen); // Generate a random TBP between initialTBP and finalTBP
        generatedPulse.randomPulse(currentTBP);
        field = generatedPulse.getField();

        //! Ambiguity removal
        // First, we calculate the phase at the peak intensity
        int peak_index = max_magnitude_index(field);

        // Calculate phase at peak intensity
        double phase = std::arg(field[peak_index]);

        // Create a new vector for the zero-phase field
        std::vector<std::complex<double>> zeroPhaseField(N);

        // Shift phase of each element in field by -phase
        for (int k = 0; k < N; ++k)
        {
            zeroPhaseField[k] = field[k] * std::exp(std::complex<double>(0, -phase));
        }

        // Compute the peak index again
        int maxIndex;
        double maxField = 0;
        for (int k = 0; k < N; k++)
        {
            if (std::abs(zeroPhaseField[k]) > maxField)
            {
                maxIndex = k;
                maxField = std::abs(zeroPhaseField[k]);
            }
        }

        // Shift the field so that the peak is at the center with a for loop
        std::vector<std::complex<double>> shiftedField(N);
        for (int k = 0; k < N; ++k)
        {
            shiftedField[k] = zeroPhaseField[(k + maxIndex + N / 2) % N];
        }

        // Replace original field with zeroPhaseField
        field = shiftedField;

        // Calculate the area of the first half of the pulse and the second half of the pulse.
        double areaFirstHalf = 0;
        double areaSecondHalf = 0;
        for (int k = 0; k < N / 2; k++)
        {
            areaFirstHalf += std::abs(field[k]) * deltaT;
        }
        for (int k = N / 2; k < N; k++)
        {
            areaSecondHalf += std::abs(field[k]) * deltaT;
        }

        // 3. If the area of the first half is greater than the second half, flip the pulse.
        if (areaFirstHalf > areaSecondHalf)
        {
            std::vector<std::complex<double>> flippedField(N);
            for (int k = 0; k < N; k++)
            {
                flippedField[k] = field[N - k];
            }
            field = flippedField;
        }

        // Set the corrected field to the pulse object
        generatedPulse.setField(field);

        // Create a group for each pulse
        H5::Group pulseGroup = file.createGroup("/pulse_" + std::to_string(i));

        // Write currentTBP
        {
            H5::DataSpace dataspace(H5S_SCALAR);
            H5::Attribute attribute = pulseGroup.createAttribute("TBP", H5::PredType::NATIVE_DOUBLE, dataspace);
            attribute.write(H5::PredType::NATIVE_DOUBLE, &currentTBP);
        }

        // Write real part of original field
        {
            hsize_t dims[1] = {static_cast<hsize_t>(N)};
            H5::DataSpace dataspace(1, dims);
            H5::DataSet dataset = pulseGroup.createDataSet("real_original_field", H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> realField(N);
            for (int k = 0; k < N; ++k)
            {
                realField[k] = std::real(field[k]) / maxField;
            }
            dataset.write(realField.data(), H5::PredType::NATIVE_DOUBLE);
        }

        // Write imaginary part of field
        {
            hsize_t dims[1] = {static_cast<hsize_t>(N)};
            H5::DataSpace dataspace(1, dims);
            H5::DataSet dataset = pulseGroup.createDataSet("imag_original_field", H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> imagField(N);
            for (int k = 0; k < N; ++k)
            {
                imagField[k] = std::imag(field[k]) / maxField;
            }
            dataset.write(imagField.data(), H5::PredType::NATIVE_DOUBLE);
        }

        originalTrace = generatedPulse.getTrace();
        // Add gaussian noise to the trace
        noisyTrace = add_noise(originalTrace, N, noiseLevel);

        // Retrieve the pulse from the noisy trace using COPRA
        COPRA retriever(ft, noisyTrace);
        auto startRetrievalTime = std::chrono::high_resolution_clock::now();
        Pulse retrievedPulse = retriever.retrieve(tolerance, maximumIterations);
        auto endRetrievalTime = std::chrono::high_resolution_clock::now();
        auto retrievalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endRetrievalTime - startRetrievalTime);
        // Save the error and the retrieval time
        errors[i] = retriever.allTraceErrors[retriever.allTraceErrors.size() - 1];
        retrievalTimes[i] = retrievalDuration.count();

        // Remove ambiguity from the retrieved pulse
        field = retrievedPulse.getField();

        //! Ambiguity removal
        // First, we calculate the phase at the peak intensity
        peak_index = max_magnitude_index(field);

        // Calculate phase at peak intensity
        phase = std::arg(field[peak_index]);

        // Shift phase of each element in field by -phase
        for (int k = 0; k < N; ++k)
        {
            zeroPhaseField[k] = field[k] * std::exp(std::complex<double>(0, -phase));
        }

        // Compute the peak index again
        maxField = 0;
        for (int k = 0; k < N; k++)
        {
            if (std::abs(zeroPhaseField[k]) > maxField)
            {
                maxIndex = k;
                maxField = std::abs(zeroPhaseField[k]);
            }
        }

        // Shift the field so that the peak is at the center with a for loop
        for (int k = 0; k < N; ++k)
        {
            shiftedField[k] = zeroPhaseField[(k + maxIndex + N / 2) % N];
        }

        // Replace original field with zeroPhaseField
        field = shiftedField;

        // Calculate the area of the first half of the pulse and the second half of the pulse.
        areaFirstHalf = 0;
        areaSecondHalf = 0;
        for (int k = 0; k < N / 2; k++)
        {
            areaFirstHalf += std::abs(field[k]) * deltaT;
        }
        for (int k = N / 2; k < N; k++)
        {
            areaSecondHalf += std::abs(field[k]) * deltaT;
        }

        // 3. If the area of the first half is greater than the second half, flip the pulse.
        if (areaFirstHalf > areaSecondHalf)
        {
            std::vector<std::complex<double>> flippedField(N);
            for (int k = 0; k < N; k++)
            {
                flippedField[k] = field[N - k];
            }
            field = flippedField;
        }

        // Normalice the trace by its maximum value
        double maxTrace = 0;
        for (int k = 0; k < N; ++k)
        {
            for (int j = 0; j < N; ++j)
            {
                if (std::abs(originalTrace[k][j]) > maxTrace)
                {
                    maxTrace = std::abs(originalTrace[k][j]);
                }
            }
        }
        for (int k = 0; k < N; ++k)
        {
            for (int j = 0; j < N; ++j)
            {
                originalTrace[k][j] /= maxTrace;
            }
        }

        // Write original trace as an NxN matrix
        {
            hsize_t dims[2] = {static_cast<hsize_t>(N), static_cast<hsize_t>(N)};
            H5::DataSpace dataspace(2, dims);
            H5::DataSet dataset = pulseGroup.createDataSet("original_trace", H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> originalTraceVector(N * N);
            for (int k = 0; k < N; ++k)
            {
                for (int j = 0; j < N; ++j)
                {
                    originalTraceVector[k * N + j] = originalTrace[k][j];
                }
            }
            dataset.write(originalTraceVector.data(), H5::PredType::NATIVE_DOUBLE);
        }

        // Normalize noisyTrace by its maximum value
        maxTrace = 0;
        for (int k = 0; k < N; ++k)
        {
            for (int j = 0; j < N; ++j)
            {
                if (std::abs(noisyTrace[k][j]) > maxTrace)
                {
                    maxTrace = std::abs(noisyTrace[k][j]);
                }
            }
        }
        for (int k = 0; k < N; ++k)
        {
            for (int j = 0; j < N; ++j)
            {
                noisyTrace[k][j] /= maxTrace;
            }
        }

        // Write noisy trace as an NxN matrix
        {
            hsize_t dims[2] = {static_cast<hsize_t>(N), static_cast<hsize_t>(N)};
            H5::DataSpace dataspace(2, dims);
            H5::DataSet dataset = pulseGroup.createDataSet("noisy_trace", H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> noisyTraceVector(N * N);
            for (int k = 0; k < N; ++k)
            {
                for (int j = 0; j < N; ++j)
                {
                    noisyTraceVector[k * N + j] = noisyTrace[k][j];
                }
            }
            dataset.write(noisyTraceVector.data(), H5::PredType::NATIVE_DOUBLE);
        }

        // Write real part of retrieved field
        {
            hsize_t dims[1] = {static_cast<hsize_t>(N)};
            H5::DataSpace dataspace(1, dims);
            H5::DataSet dataset = pulseGroup.createDataSet("real_retrieved_field", H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> realField(N);
            for (int k = 0; k < N; ++k)
            {
                realField[k] = std::real(field[k]) / maxField;
            }
            dataset.write(realField.data(), H5::PredType::NATIVE_DOUBLE);
        }

        // Write imaginary part of retrieved field
        {
            hsize_t dims[1] = {static_cast<hsize_t>(N)};
            H5::DataSpace dataspace(1, dims);
            H5::DataSet dataset = pulseGroup.createDataSet("imag_retrieved_field", H5::PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> imagField(N);
            for (int k = 0; k < N; ++k)
            {
                imagField[k] = std::imag(field[k]) / maxField;
            }
            dataset.write(imagField.data(), H5::PredType::NATIVE_DOUBLE);
        }
    }

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Close the HDF5 file
    file.close();

    // Compute the elapsed time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print the elapsed time
    std::cout << "Elapsed time generating database: " << duration.count() << " milliseconds" << std::endl;

    // Compute the average error and the average retrieval time
    double averageError = std::accumulate(errors.begin(), errors.end(), 0.0) / errors.size();
    double averageRetrievalTime = std::accumulate(retrievalTimes.begin(), retrievalTimes.end(), 0.0) / retrievalTimes.size();

    std::cout << "Average error: " << averageError << std::endl;
    std::cout << "Average retrieval time: " << averageRetrievalTime << " milliseconds" << std::endl;

    return 0;
}