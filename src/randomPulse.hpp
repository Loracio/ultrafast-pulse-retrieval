#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED

#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include "fourier.hpp"
#include "utils.hpp"

/**
 * @brief Brent's method for root finding.
 *
 * @param f The function to find the root of.
 * @param a The lower bound of the interval.
 * @param b The upper bound of the interval.
 * @param eps The tolerance for convergence.
 * @return The root of the function within the specified interval and tolerance.
 */
double brent(double (*f)(double), double a, double b, double eps)
{
    double fa = f(a); // Function value at lower bound
    double fb = f(b); // Function value at upper bound
    double c = a;     // Current estimate of the root
    double fc = fa;   // Function value at current estimate
    double d = b - a; // Interval width
    double e = d;     // Previous interval width

    while (std::abs(fc) > eps && std::abs(d) > eps)
    {
        if (fa != fc && fb != fc)
        {
            // Inverse quadratic interpolation
            double s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                       b * fa * fc / ((fb - fa) * (fb - fc)) +
                       c * fa * fb / ((fc - fa) * (fc - fb));

            if (s >= (3 * a + b) / 4 && s <= b || s <= (3 * a + b) / 4 && s >= b)
            {
                // Use s as the new estimate
                d = e;
                e = b - a;
                b = s;
                fb = f(b);
            }
            else
            {
                // Use bisection
                if (std::abs(e) < std::abs(d))
                {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
                double m = (c - b) / (fb - fc);
                double p = (b - a) * fb * fc / ((fa - fb) * (fa - fc));
                double q = (b - c) * fa * fb / ((fc - fa) * (fc - fb));
                if (p < 0)
                    p = -p;
                if (q < 0)
                    q = -q;
                if (p <= q)
                {
                    e = d;
                    d = p;
                }
                else
                {
                    e = d;
                    d = q;
                }
                a = b;
                fa = fb;
                if (std::abs(d) > eps)
                    b = b + d;
                else if (b > c)
                    b = b - eps;
                else
                    b = b + eps;
                fb = f(b);
            }
        }
        else
        {
            // Use bisection
            if (b > c)
            {
                e = d;
                d = b - a;
            }
            else
            {
                e = d;
                d = c - b;
            }
            a = b;
            fa = fb;
            if (std::abs(d) > eps)
                b = b + d;
            else if (b > c)
                b = b - eps;
            else
                b = b + eps;
            fb = f(b);
        }
        if (fa < 0)
        {
            c = a;
            fc = fa;
        }
        if (fb < 0)
        {
            a = b;
            fa = fb;
        }
    }

    return b;
}

std::vector<std::complex<double>> randomPulse(const std::vector<double> &t, double deltaT, int N, double TBP)
{
    std::vector<std::complex<double>> pulse(N);

    // Calculate some fundamental grid parameters
    std::vector<double> omega = toAngularFrequency(fftFreq(N, deltaT));
    double deltaOmega = 2 * M_PI / (N * deltaT);

    double t0 = 0.5 * (t[0] + t[N - 1]);
    double omega0 = 0.5 * (omega[0] + omega[N - 1]);

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // This is roughly the log of the roundoff error induced by an FFT
    double logEdge = N * std::numeric_limits<double>::epsilon();

    // Calculate the width of a Gaussian function that drops exactly to edge_value at the edges of the grid
    double spectralWidth = sqrt(-0.125 * (omega[0] - omega[N - 1]) * (omega[0] - omega[N - 1]) / logEdge);
    // Now the same in the temporal domain
    double maxTemporalWidth = sqrt(-0.125 * (t[0] - t[N - 1]) * (t[0] - t[N - 1]) / logEdge);
    // The actual temporal width is obtained by the uncertainty relation from the specified TBP
    double temporalWidth = 2.0 * TBP / spectralWidth;

    if (temporalWidth > maxTemporalWidth)
    {
        throw std::runtime_error("The required time-bandwidth product cannot be reached. Increase sample number.\n");
    }

    // Special case for TBP = 0.5 (transform-limited case)
    if (TBP == 0.5)
    {
        std::vector<std::complex<double>> spectrum;
        spectrum.reserve(N);

        for (int i = 0; i < N; ++i)
        {
            spectrum.push_back(std::exp(std::complex<double>(-0.5 * (omega[i] - omega0) * (omega[i] - omega0) / (spectralWidth * spectralWidth), 2 * M_PI * dist(gen))));
        }
        return IDFT(spectrum, t, deltaT, omega, deltaOmega);
    }

    // Create the filter functions, the scaling by the number of rounds is purely a heuristic
    std::vector<double> spectral_filter = gaussian(omega, omega0, spectralWidth);

    /*
        The algorithm works by iteratively filtering in the frequency and time
        domain. However, the chosen filter functions only roughly give
        the correct TBP. To obtain the exact result we scale the temporal
        filter bandwidth by a factor and perform a scalar minimization on
        that value.
    */
    std::vector<std::complex<double>> spectrum;
    spectrum.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        spectrum.push_back(std::complex<double>(dist(gen), 0.0) * std::exp(std::complex<double>(0.0, 2 * M_PI * dist(gen))));
    }

    // Rough guess for the relative range in which our optimal value lies
    double factor_min = 0.5;
    double factorMax = 1.5;


    // TODO: continue with Brent's optimization for finding specific TBP
    //? Is it necessart to build this function? This module has been thought as an improvement in just the retrieval process
    /*
    def create_pulse(factor):
        """ This performs the filtering. """
        temporal_filter = lib.gaussian(t, t0, temporal_width * factor)

        pulse.spectrum = spectrum * spectral_filter
        pulse.field = pulse.field * temporal_filter

    def objective(factor):
        """ This function should be zero """
        create_pulse(factor)
        return tbp - pulse.time_bandwidth_product
    */

    return pulse;
}

#endif // RANDOM_H_INCLUDED