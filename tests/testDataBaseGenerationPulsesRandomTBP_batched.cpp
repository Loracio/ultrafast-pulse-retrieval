/**
 * @file testDataBaseGenerationPulsesRandomTBP_batched.cpp
 * @author Víctor Loras Herrero
 * @brief Test database creation
 *
 * @copyright Copyright (c) 2024
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
    int N = 64;
    double signalDuration = 1;
    double deltaT = signalDuration / N;

    double t0 = 0;

    FourierTransform ft(N, deltaT, t0);

    Pulse generatedPulse(ft);
    std::vector<std::vector<double>> Tmn;
    std::vector<std::complex<double>> field;

    // DB parameters
    int numberOfPulses = 1000;
    int batchSize = 128;
    int nBatches = std::ceil(static_cast<double>(numberOfPulses) / batchSize);
    int pulsesPerBatch = numberOfPulses / nBatches;

    double initialTBP = 0.51;
    double finalTBP = 0.725;

    double currentTBP; // It will be randomly generated uniformly between initialTBP and finalTBP for each pulse
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(initialTBP, finalTBP);

    std::cout << "Generating " << numberOfPulses << " random pulses (" + std::to_string(pulsesPerBatch) + " pulses per batch)..."<< std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    for (int l = 0; l < nBatches; l++)
    {

        std::string filename = std::to_string(numberOfPulses) + "_randomNormalizedPulses_N" + std::to_string(N) + "_batch" + std::to_string(l) + ".h5";
        // Create the HDF5 file
        H5::H5File file(filename, H5F_ACC_TRUNC);
        // Check if the file was opened correctly
        if (!file.getId())
        {
            std::cerr << "Error opening file " << filename << std::endl;
            return -1;
        }

        for (int i = 0; i < pulsesPerBatch; i++)
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
                    flippedField[k] = std::conj(field[N - k]);
                }
                field = flippedField;
            }

            // Create a group for each pulse
            H5::Group pulseGroup = file.createGroup("/pulse_" + std::to_string(i + l * pulsesPerBatch));

            // Write currentTBP
            {
                H5::DataSpace dataspace(H5S_SCALAR);
                H5::Attribute attribute = pulseGroup.createAttribute("TBP", H5::PredType::NATIVE_DOUBLE, dataspace);
                attribute.write(H5::PredType::NATIVE_DOUBLE, &currentTBP);
            }

            // Write real part of field
            {
                hsize_t dims[1] = {static_cast<hsize_t>(N)};
                H5::DataSpace dataspace(1, dims);
                H5::DataSet dataset = pulseGroup.createDataSet("real_field", H5::PredType::NATIVE_DOUBLE, dataspace);
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
                H5::DataSet dataset = pulseGroup.createDataSet("imag_field", H5::PredType::NATIVE_DOUBLE, dataspace);
                std::vector<double> imagField(N);
                for (int k = 0; k < N; ++k)
                {
                    imagField[k] = std::imag(field[k]) / maxField;
                }
                dataset.write(imagField.data(), H5::PredType::NATIVE_DOUBLE);
            }
        }

        // Close the HDF5 file
        file.close();

        std::cout << "Finished generating batch " << l + 1 << " of " << nBatches << " batches." << std::endl;
    }

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Compute the elapsed time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // Print the elapsed time
    std::cout << "Elapsed time generating database: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}