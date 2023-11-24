#ifndef RETRIEVERS_INCLUDED
#define RETRIEVERS_INCLUDED

/**
 * @file retrievers.hpp
 * @author Víctor Loras Herrero
 * @brief Retrieval algorithms for SHG FROG traces
 *
 */

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
    double Tmeas_max_squared;               // Maximum value of the measured trace

    std::vector<double> tau;                               // Delays of the pulse in the time domain
    std::vector<std::vector<std::complex<double>>> delays; // Delays of the pulse in the frequency domain. NxN matrix that stores each delay for each frequency

    std::vector<std::vector<std::complex<double>>> Smn;     // Signal operator in frequency domain
    std::vector<std::vector<std::complex<double>>> Smk;     // Signal operator in time domain
    std::vector<std::vector<std::complex<double>>> nextSmk; // Signal operator in time domain after projection
    std::vector<std::vector<std::complex<double>>> Amk;     // Delayed pulse by τ_m
    std::vector<std::vector<double>> Tmn;                   // Trace of the resulting pulse

    double mu;                                   // Scale factor
    double r;                                    // Sum of squared residuals
    double R;                                    // Trace error
    double bestError;                            // Best achieved trace error
    std::vector<std::complex<double>> bestField; // Result of the best field for retrieval

    double Z;                                // Sum of difference between the signal operators in the frequency domain
    std::vector<std::complex<double>> gradZ; // Stores the value of the gradient of Z
    double gamma;                            // Gradient descent step

    double minimumR; // Stores the value of the lowest trace error (R) encountered

    retrieverBase(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured)
    {
        this->_ft = &ft;
        this->N = this->_ft->N;

        //! Starting pulse as a random pulse with TBP = 0.5
        this->result = new Pulse(ft);
        this->result->randomPulse(2);

        this->Tmeas = Tmeasured;
        this->Tmeas_max_squared = 0;
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                if(Tmeasured[i][j] > this->Tmeas_max_squared){
                    this->Tmeas_max_squared = Tmeasured[i][j];
                }
            }
        }
        
        this->Tmeas_max_squared *= this->Tmeas_max_squared;

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
                delayedSpectrum[i][j] = spectrum[j] * this->delays[i][j]; // delay in the time domain by
            }
            this->Amk[i] = this->_ft->backwardTransform(delayedSpectrum[i]); // E(t - τ)
        }
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

    void computeSmn()
    {
        for (int i = 0; i < this->N; i++)
        {
            this->Smn[i] = this->_ft->forwardTransform(this->Smk[i]);
        }
    }

    void computeTmn()
    {
        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                this->Tmn[i][j] = abs(this->Smn[i][j]) * abs(this->Smn[i][j]);
            }
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
                sum_meas += this->Tmeas[i][j] * this->Tmeas[i][j];
            }
        }

        this->mu = sum_meas_candidate / sum_meas;
    }

    void computeResiduals()
    {
        std::vector<std::vector<double>> difference;
        difference.resize(this->N, std::vector<double>(this->N));

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
        this->R = sqrt(this->r / (this->N * this->N * this->Tmeas_max_squared));
    }
};

class GPA : public retrieverBase
{
private: //! Get in here all the public functions when everything works OK.
public:
    GPA(FourierTransform &ft, std::vector<std::vector<double>> Tmeasured) : retrieverBase(ft, Tmeasured)
    {
    }

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
                else{
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

    void computeZ(){
        this->Z = 0;

        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                this->Z += norm(this->nextSmk[i][j] - this->Smk[i][j]);
            }   
        }
    }

    Pulse retrieve(double tolerance, double maximumIterations)
    {

        int nIter = 0;
        this->bestError = std::numeric_limits<double>::infinity();
        this->computeAmk(this->result->getSpectrum());
        this->computeSmk(this->result->getField());
        this->computeTmn();
        this->computeMu();
        this->computeResiduals();
        this->computeTraceError();

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

        return *this->result;
    }
};

#endif // RETRIEVERS_INCLUDED