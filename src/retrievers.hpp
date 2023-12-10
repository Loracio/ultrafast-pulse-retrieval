/**
 * @file retrievers.hpp
 * @author Víctor Loras Herrero
 * @brief Retriever classes for ultrashort pulse retrieval
 *
 * @copyright Copyright (c) 2023
 *
 * This module implements some retrieval methods for SHG-FROG traces, including:
 *      - Generalized Projections Algorithm (GPA)
 *      - Ptycographic Iterative Engine (PIE)
 *      - Common Phase Retrieval Algorithm (COPRA)
 *
 * Check out Nils C Geib "PyPret" module which inspired this code.
 * https://pypret.readthedocs.io/en/latest/
 * 
 * ! This file is not fully documented yet.
 *
 */

#ifndef RETRIEVERS_INCLUDED
#define RETRIEVERS_INCLUDED

#include <complex>
#include <vector>
#include "fourier.hpp"
#include "pulse.hpp"

class retrieverBase
{
public:
    FourierTransform *_ft; // Fourier transform object to perform fast fourier transforms

    int N; // Number of samples

    Pulse *result; // Resulting pulse of the retrieval

    std::vector<std::vector<double>> Tmeas; // Measured trace
    double TmeasMaxSquared;                 // Maximum value of the measured trace

    std::vector<double> tau;                               // Delays of the pulse in the time domain
    std::vector<std::vector<std::complex<double>>> delays; // Delays of the pulse in the frequency domain. NxN matrix that stores each delay for each frequency

    std::vector<std::vector<std::complex<double>>> Smn;     // Signal operator in frequency domain
    std::vector<std::vector<std::complex<double>>> Smk;     // Signal operator in time domain
    std::vector<std::vector<std::complex<double>>> nextSmk; // Signal operator in time domain after projection
    std::vector<std::vector<std::complex<double>>> Amk;     // Delayed pulse by τ_m
    std::vector<std::vector<double>> Tmn;                   // Trace of the resulting pulse

    double mu;        // Scale factor
    double r;         // Sum of squared residuals
    double R;         // Trace error
    double bestError; // Best achieved trace error

    std::vector<double> allTraceErrors; // This will store al retrieval errors during the retrieval process.

    double Z;                                // Sum of difference between the signal operators in the frequency domain
    std::vector<std::complex<double>> gradZ; // Stores the value of the gradient of Z
    double gamma;                            // Gradient descent step

    retrieverBase(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured)
    {
        this->_ft = &ft;
        this->N = this->_ft->N;

        //! Starting pulse as a random pulse with TBP = 0.5
        this->result = new Pulse(ft);
        this->result->randomPulse(0.5);

        this->Tmeas = Tmeasured;
        this->TmeasMaxSquared = 0;
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                if (Tmeasured[i][j] > this->TmeasMaxSquared)
                {
                    this->TmeasMaxSquared = Tmeasured[i][j];
                }
            }
        }

        this->TmeasMaxSquared *= this->TmeasMaxSquared;

        this->tau.reserve(this->N);
        for (int i = 0; i < this->N; i++)
        {
            this->tau[i] = (i - std::floor(0.5 * this->N)) * this->_ft->deltaT;
        }

        this->delays.resize(this->N, std::vector<std::complex<double>>(this->N));

        for (int i = 0; i < this->N; i++) // iterates through delay values
        {
            for (int j = 0; j < this->N; j++) // iterates through frequency values
            {
                this->delays[i][j] = std::exp(std::complex<double>(0, this->_ft->omega[j] * tau[i])); // delay in the time domain by τ
            }
        }

        this->Smn.resize(this->N, std::vector<std::complex<double>>(this->N));
        this->Smk.resize(this->N, std::vector<std::complex<double>>(this->N));
        this->nextSmk.resize(this->N, std::vector<std::complex<double>>(this->N));
        this->Amk.resize(this->N, std::vector<std::complex<double>>(this->N));
        this->Tmn.resize(this->N, std::vector<double>(this->N));

        this->gradZ.reserve(this->N);
    }

    void computeAmk(const std::vector<std::complex<double>> &spectrum)
    {
        std::vector<std::vector<std::complex<double>>> delayedSpectrum(this->N, std::vector<std::complex<double>>(this->N));
        for (int i = 0; i < this->N; i++) // iterates through delay values
        {
            for (int j = 0; j < this->N; j++) // iterates through frequency values
            {
                delayedSpectrum[i][j] = spectrum[j] * this->delays[i][j]; // delay in the time domain by τ
            }
            this->Amk[i] = this->_ft->backwardTransform(delayedSpectrum[i]); // E(t - τ)
        }
    }

    void computeAmk(const std::vector<std::complex<double>> &spectrum, int randomIndex)
    {
        std::vector<std::complex<double>> delayedSpectrum(this->N);
        for (int j = 0; j < this->N; j++)
        {
            delayedSpectrum[j] = spectrum[j] * this->delays[randomIndex][j]; // delay in the time domain by τ
        }
        this->Amk[randomIndex] = this->_ft->backwardTransform(delayedSpectrum); // E(t - τ)
    }

    void computeSmk(const std::vector<std::complex<double>> &field)
    {
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                this->Smk[i][j] = Amk[i][j] * field[j]; // E(t - τ) E(t)
            }
        }
    }

    void computeSmk(const std::vector<std::complex<double>> &field, int randomIndex)
    {
        for (int j = 0; j < this->N; j++)
        {
            this->Smk[randomIndex][j] = Amk[randomIndex][j] * field[j]; // E(t - τ) E(t)
        }
    }

    void computeSmn()
    {
        for (int i = 0; i < this->N; i++)
        {
            this->Smn[i] = this->_ft->forwardTransform(this->Smk[i]);
        }
    }

    void computeSmn(int randomIndex)
    {

        this->Smn[randomIndex] = this->_ft->forwardTransform(this->Smk[randomIndex]);
    }

    void computeTmn()
    {
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                this->Tmn[i][j] = std::norm(this->Smn[i][j]);
            }
        }
    }

    void computeTmn(int randomIndex)
    {
        for (int j = 0; j < this->N; j++)
        {
            this->Tmn[randomIndex][j] = std::norm(this->Smn[randomIndex][j]);
        }
    }

    void computeMu()
    {
        double sum_meas_candidate = 0;
        double sum_meas = 0;
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                sum_meas_candidate += this->Tmeas[i][j] * this->Tmn[i][j];
                sum_meas += this->Tmn[i][j] * this->Tmn[i][j];
            }
        }

        this->mu = sum_meas_candidate / sum_meas;
    }

    void computeResiduals()
    {
        std::vector<std::vector<double>> difference(this->N, std::vector<double>(this->N));

        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                difference[i][j] = this->Tmeas[i][j] - this->mu * this->Tmn[i][j];
            }
        }

        // Calculate the sum of squared differences
        double sum = 0.0;
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                sum += difference[i][j] * difference[i][j];
            }
        }

        this->r = sum;
    }

    void computeTraceError()
    {
        this->R = sqrt(this->r / (this->N * this->N * this->TmeasMaxSquared));
    }

    void setInitialField(const std::vector<std::complex<double>> &initialField)
    {
        this->result->setField(initialField);
    }

    void setInitialSpectrum(const std::vector<std::complex<double>> &initialSpectrum)
    {
        this->result->setSpectrum(initialSpectrum);
    }
};

class GPA : public retrieverBase
{
private:
    std::vector<std::complex<double>> bestField; // Result of the best field for retrieval
    void computeNextSmk()
    {
        std::vector<std::vector<double>> absSmn(this->N, std::vector<double>(this->N));
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                absSmn[i][j] = std::abs(this->Smn[i][j]);
            }
        }

        std::vector<std::complex<double>> nextSmn(this->N);
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                if (absSmn[i][j] > 0.0)
                {
                    nextSmn[j] = this->Smn[i][j] / absSmn[i][j] * sqrt(this->Tmeas[i][j]);
                }
                else
                {
                    nextSmn[j] = sqrt(this->Tmeas[i][j]);
                }
            }

            this->nextSmk[i] = this->_ft->backwardTransform(nextSmn);
        }
    }

    void computeGamma()
    {
        double sumGradZ = 0;
        for (int i = 0; i < this->N; i++)
        {
            sumGradZ += std::abs(this->gradZ[i]) * std::abs(this->gradZ[i]);
        }

        this->gamma = this->Z / sumGradZ;
    }

    void computeGradient()
    {
        // Calculate dS
        std::vector<std::vector<std::complex<double>>> dS(this->N, std::vector<std::complex<double>>(this->N));
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                dS[i][j] = this->nextSmk[i][j] - this->Smk[i][j]; //! This can be computed first in Z and save some time if stored
            }
        }

        for (int j = 0; j < this->N; j++)
        {
            this->gradZ[j] = 0;

            for (int m = 0; m < this->N; m++)
            {
                this->gradZ[j] += dS[m][j] * std::conj(this->Amk[m][j]);
                if (j + m < N)
                {
                    this->gradZ[j] += dS[m][j + m] * std::conj(this->Amk[m][j + m]);
                }
                else
                {
                    this->gradZ[j] += dS[m][j + m - this->N] * std::conj(this->Amk[m][j + m - this->N]);
                }
            }

            this->gradZ[j] *= -2;
        }
    }

    void computeNextField()
    {
        std::vector<std::complex<double>> currentField = this->result->getField();
        for (int i = 0; i < this->N; i++)
        {
            currentField[i] -= this->gamma * this->gradZ[i];
        }

        this->result->setField(currentField);
        this->result->updateSpectrum();
    }

    void computeZ()
    {
        this->Z = 0;

        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                this->Z += norm(this->nextSmk[i][j] - this->Smk[i][j]);
            }
        }
    }

public:
    GPA(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured) : retrieverBase(ft, Tmeasured)
    {
    }

    GPA(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured, std::vector<double> measuredDelays) : retrieverBase(ft, Tmeasured)
    {
        // Set up the delays by the given measured delays
        for (int i = 0; i < this->N; i++) // iterates through delay values
        {
            for (int j = 0; j < this->N; j++) // iterates through frequency values
            {
                this->delays[i][j] = std::exp(std::complex<double>(0, this->_ft->omega[j] * measuredDelays[i])); // delay in the time domain by τ
            }
        }
    }

    Pulse retrieve(double tolerance, double maximumIterations)
    {

        int nIter = 0;
        this->bestError = std::numeric_limits<double>::infinity();
        this->computeAmk(this->result->getSpectrum());
        this->computeSmk(this->result->getField());
        this->computeSmn();
        this->computeTmn();
        this->computeMu();
        this->computeResiduals();
        this->computeTraceError();
        this->allTraceErrors.push_back(this->R);

        while (this->R > tolerance && nIter < maximumIterations)
        {
            this->computeNextSmk();

            this->computeZ();
            this->computeGradient();
            this->computeGamma();

            this->computeNextField();

            this->computeAmk(this->result->getSpectrum());
            this->computeSmk(this->result->getField());
            this->computeSmn();
            this->computeTmn();
            this->computeMu();
            this->computeResiduals();
            this->computeTraceError();
            this->allTraceErrors.push_back(this->R);

            if (this->R < this->bestError)
            {
                this->bestError = this->R;
                this->bestField = this->result->getField();
            }

            std::cout << "Iteration = " << nIter + 1 << "\t"
                      << "R = " << this->R << std::endl;

            nIter++;
        }

        std::cout << "Best retrieval error R = " << this->bestError << std::endl;

        this->result->setField(this->bestField);
        this->result->updateSpectrum();

        this->allTraceErrors.push_back(this->bestError); //! The last value of the array is the best result. Not the last retrieval result.

        return *this->result;
    }
};

class PIE : public retrieverBase
{
private:
    std::vector<std::complex<double>> bestField; // Result of the best field for retrieval

    void computeNextSmk()
    {
        std::vector<std::vector<double>> absSmn(this->N, std::vector<double>(this->N));
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                absSmn[i][j] = std::abs(this->Smn[i][j]);
            }
        }

        std::vector<std::complex<double>> nextSmn(this->N);
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                if (absSmn[i][j] > 0.0)
                {
                    nextSmn[j] = this->Smn[i][j] / absSmn[i][j] * sqrt(this->Tmeas[i][j] / this->mu);
                }
                else
                {
                    nextSmn[j] = sqrt(this->Tmeas[i][j] / this->mu);
                }
            }

            this->nextSmk[i] = this->_ft->backwardTransform(nextSmn);
        }
    }

    void computeNextSmk(int randomIndex)
    {
        std::vector<double> absSmn(this->N);
        for (int j = 0; j < this->N; ++j)
        {
            absSmn[j] = std::abs(this->Smn[randomIndex][j]);
        }

        std::vector<std::complex<double>> nextSmn(this->N);
        for (int j = 0; j < this->N; ++j)
        {
            if (absSmn[j] > 0.0)
            {
                nextSmn[j] = this->Smn[randomIndex][j] / absSmn[j] * sqrt(this->Tmeas[randomIndex][j] / this->mu);
            }
            else
            {
                nextSmn[j] = sqrt(this->Tmeas[randomIndex][j] / this->mu);
            }
        }

        this->nextSmk[randomIndex] = this->_ft->backwardTransform(nextSmn);
    }

    void computeNextField(int randomIndex, double beta)
    {
        this->computeAmk(this->result->getSpectrum(), randomIndex);
        this->computeSmk(this->result->getField(), randomIndex);
        this->computeSmn(randomIndex);
        this->computeTmn(randomIndex);

        // Compute projection on Smk
        this->computeNextSmk(randomIndex);

        std::vector<std::complex<double>> currentField = this->result->getField();

        double currentAbsMaxValue = 0;
        double absValue;

        currentAbsMaxValue = 0;
        for (int j = 0; j < this->N; j++)
        {
            absValue = std::norm(currentField[j]);
            if (currentAbsMaxValue < absValue)
            {
                currentAbsMaxValue = absValue;
            }
        }

        for (int j = 0; j < this->N; j++)
        {
            currentField[j] += beta * std::conj(this->Amk[randomIndex][j]) * (this->nextSmk[randomIndex][j] - this->Smk[randomIndex][j]) / currentAbsMaxValue;
        }

        this->result->setField(currentField);
    }

    std::vector<int> randomIndexShuffle()
    {
        std::vector<int> indices(this->N);
        std::random_device rd;
        std::mt19937 rng(rd());

        std::iota(indices.begin(), indices.end(), 0); // Fill indices with 0, 1, ..., N-1

        // Shuffle the array of indices
        std::shuffle(indices.begin(), indices.end(), rng);

        return indices;
    }

public:
    PIE(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured) : retrieverBase(ft, Tmeasured)
    {
    }

    PIE(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured, std::vector<double> measuredDelays) : retrieverBase(ft, Tmeasured)
    {
        // Set up the delays by the given measured delays
        for (int i = 0; i < this->N; i++) // iterates through delay values
        {
            for (int j = 0; j < this->N; j++) // iterates through frequency values
            {
                this->delays[i][j] = std::exp(std::complex<double>(0, this->_ft->omega[j] * measuredDelays[i])); // delay in the time domain by τ
            }
        }
    }

    Pulse retrieve(double tolerance, double maximumIterations)
    {

        int nIter = 0;
        this->bestError = std::numeric_limits<double>::infinity();
        this->computeAmk(this->result->getSpectrum());
        this->computeSmk(this->result->getField());
        this->computeSmn();
        this->computeTmn();
        this->computeMu();
        this->computeResiduals();
        this->computeTraceError();
        this->allTraceErrors.push_back(this->R);

        std::vector<int> randomIndexes;

        // Set up a random number generator for beta
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(0.1, 0.5);

        double beta;

        while (this->R > tolerance && nIter < maximumIterations)
        {
            randomIndexes = this->randomIndexShuffle();
            beta = distribution(gen);
            for (int i = 0; i < this->N; i++)
            {
                this->computeNextField(randomIndexes[i], beta);
            }

            this->computeAmk(this->result->getSpectrum());
            this->computeSmk(this->result->getField());
            this->computeSmn();
            this->computeTmn();
            this->computeMu();
            this->computeResiduals();
            this->computeTraceError();
            this->allTraceErrors.push_back(this->R);

            if (this->R < this->bestError)
            {
                this->bestError = this->R;
                this->bestField = this->result->getField();
            }

            std::cout << "Iteration = " << nIter + 1 << "\t"
                      << "R = " << this->R << std::endl;

            nIter++;
        }

        std::cout << "Best retrieval error R = " << this->bestError << std::endl;

        this->result->setField(this->bestField);
        this->result->updateSpectrum();

        this->allTraceErrors.push_back(this->bestError); //! The last value of the array is the best result. Not the last retrieval result.

        return *this->result;
    }
};

class COPRA : public retrieverBase
{
private:
    double previousMaxGradient;
    double currentMaxGradient;
    double etar;
    double etaz;

    double alpha = 0.25; //! Should change as an argument in some function.

    std::vector<std::vector<std::complex<double>>> gradrmk;

    std::vector<std::complex<double>> bestSpectrum; // Result of the best spectrum for retrieval

    void computeNextSmk()
    {
        std::vector<std::vector<double>> absSmn(this->N, std::vector<double>(this->N));
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                absSmn[i][j] = std::abs(this->Smn[i][j]);
            }
        }

        std::vector<std::complex<double>> nextSmn(this->N);
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                if (absSmn[i][j] > 0.0)
                {
                    nextSmn[j] = this->Smn[i][j] / absSmn[i][j] * sqrt(this->Tmeas[i][j] / this->mu);
                }
                else
                {
                    nextSmn[j] = sqrt(this->Tmeas[i][j] / this->mu);
                }
            }

            this->nextSmk[i] = this->_ft->backwardTransform(nextSmn);
        }
    }

    void computeNextSmk(int randomIndex)
    {
        std::vector<double> absSmn(this->N);
        for (int j = 0; j < this->N; ++j)
        {
            absSmn[j] = std::abs(this->Smn[randomIndex][j]);
        }

        std::vector<std::complex<double>> nextSmn(this->N);
        for (int j = 0; j < this->N; ++j)
        {
            if (absSmn[j] > 0.0)
            {
                nextSmn[j] = this->Smn[randomIndex][j] / absSmn[j] * sqrt(this->Tmeas[randomIndex][j] / this->mu);
            }
            else
            {
                nextSmn[j] = sqrt(this->Tmeas[randomIndex][j] / this->mu);
            }
        }

        this->nextSmk[randomIndex] = this->_ft->backwardTransform(nextSmn);
    }

    void computeGradZ()
    {
        std::vector<std::complex<double>> dSmkEk(this->N);
        std::vector<std::complex<double>> dSmkAmk(this->N);
        std::vector<std::complex<double>> currentField = this->result->getField();

        for (int i = 0; i < this->N; i++)
        {
            this->gradZ[i] = 0;
        }

        for (int m = 0; m < this->N; m++)
        {
            for (int k = 0; k < this->N; k++)
            {
                dSmkEk[k] = (this->nextSmk[m][k] - this->Smk[m][k]) * std::conj(currentField[k]);
                dSmkAmk[k] = (this->nextSmk[m][k] - this->Smk[m][k]) * std::conj(this->Amk[m][k]);
            }

            dSmkEk = this->_ft->forwardTransform(dSmkEk);
            dSmkAmk = this->_ft->forwardTransform(dSmkAmk);

            for (int n = 0; n < this->N; n++)
            {
                this->gradZ[n] += std::conj(this->delays[m][n]) * dSmkEk[n] + dSmkAmk[n];
            }
        }

        for (int i = 0; i < this->N; i++)
        {
            // Multiply by common factor
            this->gradZ[i] *= -4 * M_PI * this->_ft->deltaOmega / (this->_ft->deltaT);
        }
    }

    void computeGradZ(int chosenIndex)
    {
        std::vector<std::complex<double>> dSmkEk(this->N);
        std::vector<std::complex<double>> dSmkAmk(this->N);
        std::vector<std::complex<double>> currentField = this->result->getField();
        for (int i = 0; i < this->N; ++i)
        {
            dSmkEk[i] = (this->nextSmk[chosenIndex][i] - this->Smk[chosenIndex][i]) * std::conj(currentField[i]);
            dSmkAmk[i] = (this->nextSmk[chosenIndex][i] - this->Smk[chosenIndex][i]) * std::conj(this->Amk[chosenIndex][i]);
        }

        dSmkEk = this->_ft->forwardTransform(dSmkEk);
        dSmkAmk = this->_ft->forwardTransform(dSmkAmk);

        for (int i = 0; i < this->N; i++)
        {
            this->gradZ[i] = -4 * M_PI * this->_ft->deltaOmega / (this->_ft->deltaT) * (std::conj(this->delays[chosenIndex][i]) * dSmkEk[i] + dSmkAmk[i]);
        }
    }

    void computeZ()
    {
        this->Z = 0;

        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                this->Z += std::norm(this->nextSmk[i][j] - this->Smk[i][j]);
            }
        }
    }

    void computeZ(int randomIndex)
    {
        this->Z = 0;

        for (int j = 0; j < this->N; j++)
        {
            this->Z += std::norm(this->nextSmk[randomIndex][j] - this->Smk[randomIndex][j]);
        }
    }

    void computeGradrmk()
    {
        std::vector<std::complex<double>> difference(this->N);
        for (int i = 0; i < this->N; ++i)
        {
            for (int j = 0; j < this->N; ++j)
            {
                difference[j] = -2 * this->mu * this->_ft->deltaT / (M_PI * this->_ft->deltaOmega) * (this->Tmeas[i][j] - this->mu * this->Tmn[i][j]) * this->Smn[i][j];
            }

            this->gradrmk[i] = this->_ft->backwardTransform(difference);
        }
    }

    void computeNextSpectrum(double step, const std::vector<std::complex<double>> &gradient)
    {
        std::vector<std::complex<double>> currentSpectrum = this->result->getSpectrum();
        for (int i = 0; i < this->N; i++)
        {
            currentSpectrum[i] -= step * gradient[i];
        }

        this->result->setSpectrum(currentSpectrum);
    }

    void localIteration(int randomIndex)
    {
        this->computeAmk(this->result->getSpectrum(), randomIndex);
        this->computeSmk(this->result->getField(), randomIndex);
        this->computeSmn(randomIndex);
        this->computeTmn(randomIndex);

        // Compute projection on Smk
        this->computeNextSmk(randomIndex);

        this->computeZ(randomIndex);
        this->computeGradZ(randomIndex);
        double gradNorm = this->computeGradZNorm();
        if (gradNorm > this->currentMaxGradient)
        {
            this->currentMaxGradient = gradNorm;
        }

        this->gamma = this->Z;
        if (this->currentMaxGradient > this->previousMaxGradient)
        {
            this->gamma /= this->currentMaxGradient;
        }
        else
        {
            this->gamma /= this->previousMaxGradient;
        }

        this->computeNextSpectrum(this->gamma, this->gradZ);
    }

    void globalIteration()
    {
        this->computeAmk(this->result->getSpectrum());
        this->computeSmk(this->result->getField());
        this->computeSmn();
        this->computeTmn();
        this->computeMu();
        this->computeResiduals();

        this->computeGradrmk();
        this->etar = this->alpha * this->r / this->computeGradrmkNorm();

        this->nextSmkGradDescent();

        this->computeZ();
        this->computeGradZ();

        this->etaz = this->alpha * this->Z / this->computeGradZNorm();

        this->computeNextSpectrum(this->etaz, this->gradZ);
    }

    std::vector<int> randomIndexShuffle()
    {
        std::vector<int> indices(this->N);
        std::random_device rd;
        std::mt19937 rng(rd());

        std::iota(indices.begin(), indices.end(), 0); // Fill indices with 0, 1, ..., N-1

        // Shuffle the array of indices
        std::shuffle(indices.begin(), indices.end(), rng);

        return indices;
    }

    double computeGradZNorm()
    {
        double sum = 0;
        for (int i = 0; i < this->N; i++)
        {
            sum += std::norm(this->gradZ[i]);
        }
        return sum;
    }

    double computeGradrmkNorm()
    {
        double sum = 0;
        for (int m = 0; m < this->N; m++)
        {
            for (int k = 0; k < this->N; k++)
            {
                sum += std::norm(this->gradrmk[m][k]);
            }
        }
        return sum;
    }

    void nextSmkGradDescent()
    {
        for (int m = 0; m < this->N; m++)
        {
            for (int k = 0; k < this->N; k++)
            {
                this->nextSmk[m][k] = this->Smk[m][k] - this->etar * this->gradrmk[m][k];
            }
        }
    }

public:
    COPRA(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured) : retrieverBase(ft, Tmeasured)
    {
        this->gradrmk.resize(this->N, std::vector<std::complex<double>>(this->N));
        this->bestSpectrum.reserve(this->N);
    }

    COPRA(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured, std::vector<double> measuredDelays) : retrieverBase(ft, Tmeasured)
    {
        this->gradrmk.resize(this->N, std::vector<std::complex<double>>(this->N));
        this->bestSpectrum.reserve(this->N);

        // Set up the delays by the given measured delays
        for (int i = 0; i < this->N; i++) // iterates through delay values
        {
            for (int j = 0; j < this->N; j++) // iterates through frequency values
            {
                this->delays[i][j] = std::exp(std::complex<double>(0, this->_ft->omega[j] * measuredDelays[i])); // delay in the time domain by τ
            }
        }
    }

    Pulse retrieve(double tolerance, double maximumIterations)
    {
        int nIter = 0;
        int stepsSinceLastImprovement = 0;
        bool mode = 1;
        std::vector<int> randomIndexes;

        this->bestError = std::numeric_limits<double>::infinity();

        this->computeAmk(this->result->getSpectrum());
        this->computeSmk(this->result->getField());
        this->computeSmn();
        this->computeTmn();
        this->computeMu();
        this->computeResiduals();
        this->computeTraceError();
        this->allTraceErrors.push_back(this->R);

        this->computeNextSmk();
        this->currentMaxGradient = 0;
        for (int m = 0; m < this->N; m++)
        {
            this->computeGradZ(m);

            this->previousMaxGradient = this->computeGradZNorm();
            if (this->previousMaxGradient > this->currentMaxGradient)
            {
                this->currentMaxGradient = this->previousMaxGradient;
            }
        }

        while (this->R > tolerance && nIter < maximumIterations)
        {

            if (mode)
            {
                this->previousMaxGradient = this->currentMaxGradient;
                this->currentMaxGradient = 0;

                randomIndexes = this->randomIndexShuffle();
                for (int i = 0; i < this->N; i++)
                {
                    this->localIteration(randomIndexes[i]);
                }

                this->computeMu();
                this->computeResiduals();
                this->computeTraceError(); // This trace error is with the approximation, as the spectrum changed every iteration
                this->allTraceErrors.push_back(this->R);

                if (this->R >= this->bestError)
                {
                    stepsSinceLastImprovement++;
                    if (stepsSinceLastImprovement == 5)
                    {
                        mode = 0;
                        std::cout << "Local iteration ended, starting global iteration" << std::endl;
                        // We pick the best result from the local iteration to start the global iteration
                        this->result->setSpectrum(this->bestSpectrum);
                        // this->result->updateField();!!Not necesssary, it is updated in setSpectrum
                    }
                }
                else
                {
                    stepsSinceLastImprovement = 0;
                    this->bestError = this->R;
                    this->bestSpectrum = this->result->getSpectrum();
                }
            }
            else
            {
                this->globalIteration();

                this->computeAmk(this->result->getSpectrum());
                this->computeSmk(this->result->getField());
                this->computeSmn();
                this->computeTmn();
                this->computeMu();
                this->computeResiduals();
                this->computeTraceError();
                this->allTraceErrors.push_back(this->R);

                if (this->R < this->bestError)
                {
                    this->bestError = this->R;
                    this->bestSpectrum = this->result->getSpectrum();
                }
            }

            std::cout << "Iteration = " << nIter + 1 << "\t"
                      << "R = " << this->R << std::endl;

            nIter++;
        }

        this->result->setSpectrum(this->bestSpectrum);
        if (mode)
        { // If the result was achieved in local iteration, compute R as it was shown as an approximated value.
            this->computeAmk(this->result->getSpectrum());
            this->computeSmk(this->result->getField());
            this->computeSmn();
            this->computeTmn();
            this->computeMu();
            this->computeResiduals();
            this->computeTraceError();
            this->bestError = this->R;
            std::cout << "Best retrieval error found in local iteration" << std::endl;
        }

        std::cout << "Best retrieval error R = " << this->bestError << std::endl;

        this->allTraceErrors.push_back(this->bestError); //! Last stored value is the best result, not the last computed error!

        return *this->result;
    }
};

#endif // RETRIEVERS_INCLUDED