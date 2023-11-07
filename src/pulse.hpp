#ifndef PULSE_INCLUDED
#define PULSE_INCLUDED

/**
* @file pulse.hpp
* @author Víctor Loras Herrero
* @brief Pulse class which contains electric field in time domain and in frequency domain, along with other pulse properties
*
*/

#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <functional>
#include "fourier.hpp"
#include "utils.hpp"

//! cite author of this function
//! https://codereview.stackexchange.com/questions/103762/implementation-of-brents-algorithm-to-find-roots-of-a-polynomial
void brents_fun(std::function<double (double)> f, double lower_bound, double upper_bound, double TOL, double MAX_ITER)
{
    double a = lower_bound;
    double b = upper_bound;
    double fa = f(a);   // calculated now to save function calls
    double fb = f(b);   // calculated now to save function calls
    double fs = 0;      // initialize 

    if (!(fa * fb < 0))
    {
        std::cout << "Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
        return;
    }

    if (std::abs(fa) < std::abs(b)) // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
    {
        std::swap(a,b);
        std::swap(fa,fb);
    }

    double c = a;           // c now equals the largest magnitude of the lower and upper bounds
    double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
    bool mflag = true;      // boolean flag used to evaluate if statement later on
    double s = 0;           // Our Root that will be returned
    double d = 0;           // Only used if mflag is unset (mflag == false)

    for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
    {
        // stop if converged on root or error is less than tolerance
        if (std::abs(b-a) < TOL)
        {
            // std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
            return;
        } // end if

        if (fa != fc && fb != fc)
        {
            // use inverse quadratic interopolation
            s =   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
                + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
                + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
        }
        else
        {
            // secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        /*
            Crazy condition statement!:
            -------------------------------------------------------
            (condition 1) s is not between  (3a+b)/4  and b or
            (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
            (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
            (condition 4) (mflag is set and |b−c| < |TOL|) or
            (condition 5) (mflag is false and |c−d| < |TOL|)
        */
        if (    ( (s < (3 * a + b) * 0.25) || (s > b) ) ||
                ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
                ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
                ( mflag && (std::abs(b-c) < TOL) ) ||
                ( !mflag && (std::abs(c-d) < TOL))  )
        {
            // bisection method
            s = (a+b)*0.5;

            mflag = true;
        }
        else
        {
            mflag = false;
        }

        fs = f(s);  // calculate fs
        d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
        c = b;      // set c equal to upper bound
        fc = fb;    // set f(c) = f(b)

        if ( fa * fs < 0)   // fa and fs have opposite signs
        {
            b = s;
            fb = fs;    // set f(b) = f(s)
        }
        else
        {
            a = s;
            fa = fs;    // set f(a) = f(s)
        }

        if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
        {
            std::swap(a,b);     // swap a and b
            std::swap(fa,fb);   // make sure f(a) and f(b) are correct after swap
        }

    } // end for

    std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
} // end brent_fun

class Pulse
{
private:
    FourierTransform *_ft;
    // double omega0;

    std::vector<std::complex<double>> _field;
    std::vector<std::complex<double>> _spectrum;

    int N;
    double deltaT;
    double deltaOmega;

    std::vector<double> t;
    std::vector<double> omega;

    // Random pulse generation variables
    double _tbp;
    double t0;
    double omega0;
    double temporalWidth;

    std::vector<std::complex<double>> candidateField;

public:
    Pulse(FourierTransform &ft, double lambda0)
    {
        this->_ft = &ft;
        // this->omega0 = 2 * M_PI * 299792458.0 / lambda0; //! this is for centered array of angular frequencies. Not yet defined

        this->N = this->_ft->N;
        this->deltaT = this->_ft->deltaT;
        this->deltaOmega = this->_ft->deltaOmega;

        this->t.reserve(this->N);
        this->t = this->_ft->t;

        this->omega.reserve(this->N);
        this->omega = this->_ft->omega;
    }

    void setField(const std::vector<std::complex<double>> &val)
    {
        this->_field = val;
        this->updateSpectrum();
    }

    void setSpectrum(const std::vector<std::complex<double>> &val)
    {
        this->_spectrum = val;
        this->updateField();
    }

    void updateField()
    {
        this->_field = this->_ft->backwardTransform(this->_spectrum);
    }

    void updateSpectrum()
    {
        this->_spectrum = this->_ft->forwardTransform(this->_field);
    }

    std::vector<std::complex<double>> getField()
    {
        return this->_field;
    }

    std::vector<std::complex<double>> getSpectrum()
    {
        return this->_spectrum;
    }

    std::vector<double> getIntensity()
    {
        std::vector<double> intensity(this->N);

        for (int i = 0; i < this->N; ++i)
        {
            intensity[i] = std::norm(this->_field[i]);
        }

        return intensity;
    }

    std::vector<double> getAmplitude()
    {
        std::vector<double> amplitude(this->N);

        for (int i = 0; i < this->N; ++i)
        {
            amplitude[i] = std::abs(this->_field[i]);
        }

        return amplitude;
    }

    std::vector<double> getPhase()
    {
        return unwrapPhase(this->_field);
    }

    std::vector<double> getSpectralIntensity()
    {
        std::vector<double> spectralIntensity(this->N);

        for (int i = 0; i < this->N; i++)
        {
            spectralIntensity[i] = std::norm(this->_spectrum[i]);
        }

        return spectralIntensity;
    }

    std::vector<double> getSpectralAmplitude()
    {
        std::vector<double> spectralAmplitude(this->N);

        for (int i = 0; i < this->N; ++i)
        {
            spectralAmplitude[i] = std::abs(this->_spectrum[i]);
        }

        return spectralAmplitude;
    }

    std::vector<double> getSpectralPhase()
    {
        return unwrapPhase(this->_spectrum);
    }

    double getTimeBandwidthProduct()
    {
        return stdDev(this->t, this->getIntensity()) * stdDev(this->omega, this->getSpectralIntensity());
    }

    double objective(double const &factor)
    {

        std::vector<double> temporalFilter = gaussian(this->t, this->t0, this->temporalWidth * factor);

        std::vector<std::complex<double>> result(this->N);

        for (int i = 0; i < this->N; i++)
        {
            result[i] = this->candidateField[i] * temporalFilter[i];
        }

        this->setField(result);

        return this->_tbp - this->getTimeBandwidthProduct();
    }

    bool randomPulse(double TBP)
    {
        this->_tbp = TBP;

        this->t0 = 0.5 * (this->t[0] + this->t[this->N - 1]);
        this->omega0 = 0.5 * (this->omega[0] + this->omega[this->N - 1]);

        // Initialize random number generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        // This is roughly the log of the roundoff error induced by an FFT
        double logEdge = std::log(this->N * std::numeric_limits<double>::epsilon());

        // Calculate the width of a Gaussian function that drops exactly to edge_value at the edges of the grid
        double spectralWidth = sqrt(-0.125 * (this->omega[0] - this->omega[N - 1]) * (this->omega[0] - this->omega[N - 1]) / logEdge);
        // Now the same in the temporal domain
        double maxTemporalWidth = sqrt(-0.125 * (this->t[0] - this->t[N - 1]) * (this->t[0] - this->t[N - 1]) / logEdge);
        // The actual temporal width is obtained by the uncertainty relation from the specified TBP
        this->temporalWidth = 2.0 * TBP / spectralWidth;

        if (this->temporalWidth > maxTemporalWidth)
        {
            throw std::runtime_error("The required time-bandwidth product cannot be reached. Increase sample number.\n");
            return 0;
        }

        // Special case for TBP = 0.5 (transform-limited case)
        if (TBP == 0.5)
        {

            for (int i = 0; i < this->N; ++i)
            {
                this->_spectrum.push_back(std::exp(std::complex<double>(-0.5 * (this->omega[i] - omega0) * (this->omega[i] - omega0) / (spectralWidth * spectralWidth), 2 * M_PI * dist(gen))));
            }

            this->updateField();

            return 1;
        }

        std::vector<double> spectralFilter = gaussian(this->omega, this->omega0, spectralWidth);

        /*
            The algorithm works by iteratively filtering in the frequency and time
            domain. However, the chosen filter functions only roughly give
            the correct TBP. To obtain the exact result we scale the temporal
            filter bandwidth by a factor and perform a scalar minimization on
            that value.
        */

        // Rough guess for the relative range in which our optimal value lies
        double factorMin = 0.5;
        double factorMax = 1.5;

        std::vector<std::complex<double>> candidateSpectrum(N);

        for (int i = 0; i < N; ++i)
        {
            candidateSpectrum[i] = spectralFilter[i] * dist(gen) * std::exp(std::complex<double>(0.0, 2 * M_PI * dist(gen)));
        }

        this->candidateField = this->_ft->backwardTransform(candidateSpectrum);

        // Ensure the objective function changes sign in the chosen bounds
        int iters = 0;
        while (std::signbit(objective(factorMin)) == std::signbit(objective(factorMax)))
        {
            // For some random arrays, this condition is not always fulfilled. Try again.
            for (int i = 0; i < N; ++i)
            {
                candidateSpectrum[i] = spectralFilter[i] * dist(gen) * std::exp(std::complex<double>(0.0, 2 * M_PI * dist(gen)));
            }

            this->candidateField = this->_ft->backwardTransform(candidateSpectrum);

            iters++;

            if (iters == 100)
            {
                throw std::runtime_error("Could not create a pulse for these parameters!");
                return 0;
            }
        }

        // Create a callable object (lambda) that wraps objective function
        auto objectiveFunction = [&](double factor)
        {
            return this->objective(factor);
        };

        brents_fun(objectiveFunction, factorMin, factorMax, 1e-12, 1000);

        return 1;
    }
};

#endif // PULSE_INCLUDED