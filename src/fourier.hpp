/**
 * @file fourier.hpp
 * @author Víctor Loras Herrero
 * @brief Functions related to the numerical computing of the Fourier Transform using FFTW library (https://www.fftw.org/)
 *
 * @copyright Copyright (c) 2023
 *
 * Choosing the Fourier Transform definition as:
 *
 *      Ẽ(ω) = ∫E(t)e^{-i ω t} dt    ;    E(t) = 1/2π ∫Ẽ(ω)e^{i t ω} dω
 *
 * And discretizing it, the nth and the jth coefficient (direct and inverse transforms) will be:
 *
 *      Ẽ(ωₙ) := Ẽₙ = ∑ⱼ₌₀ᴺ⁻¹ E(tⱼ) e^{-i ωₙ tⱼ} Δt    ;     E(tⱼ) := Eⱼ = 1/2π · ∑ₙ₌₀ᴺ⁻¹ Ẽ(ωₙ) e^{i tⱼ ωₙ} Δω
 *
 * Where tⱼ is the jth element of the time array and ωₙ is the nth element of the frequency array.
 *
 * The time array will be of the form: tⱼ = t₀ + j·Δt with j = 0, ..., N - 1, where N-1 is the number of samples.
 * The frequency array will be of the form ωₙ = ω₀ + n·Δω with n = 0, ..., N - 1.
 *
 * Considering the reciprocity relation, which states that:
 *
 *      Δt Δω = 2π/N
 *
 * We can substitute it into the obtained discretized expressions:
 *
 *      Ẽₙ = ∑ⱼ₌₀ᴺ⁻¹ Eⱼ e^{-i (ω₀ + n·Δω) (t₀ + j·Δt)} Δt =
 *      = Δt e^{-i n t₀ Δω} ∑ⱼ₌₀ᴺ⁻¹ Eⱼ e^{-i ω₀ tⱼ} e^{-i n j Δω Δt} =
 *      = Δt e^{-i n t₀ Δω} ∑ⱼ₌₀ᴺ⁻¹ Eⱼ e^{-i ω₀ tⱼ} e^{-i 2π n j / N}
 *
 *      Eⱼ = 1/2π · ∑ₙ₌₀ᴺ⁻¹ Ẽₙ e^{i (t₀ + j·Δt) (ω₀ + n·Δω)} Δω =
 *      = Δω/2π e^{i ω₀ tⱼ} ∑ₙ₌₀ᴺ⁻¹ Ẽₙ e^{i n t₀ Δω} e^{i n j Δt Δω} =
 *      = Δω/2π e^{i ω₀ tⱼ} ∑ₙ₌₀ᴺ⁻¹ Ẽₙ e^{i n t₀ Δω} e^{i 2π n j / N} =
 *      = 1/Δt e^{i ω₀ tⱼ} · 1/N ∑ₙ₌₀ᴺ⁻¹ Ẽₙ e^{i n t₀ Δω} e^{i 2π n j / N}
 *
 * If we define:
 *
 *      rₙ = e^{-i n t₀ Δω} ; sⱼ = e^{-i ω₀ tⱼ}
 *
 * We can finally express the Discrete Fourier Transform (DFT) as:
 *
 *      Ẽₙ = Δt·rₙ ∑ⱼ₌₀ᴺ⁻¹ Eⱼ·sⱼ e^{-i 2π n j / N}    ;    Eⱼ = 1/Δt·sⱼ* · 1/N ∑ₙ₌₀ᴺ⁻¹ Ẽₙ·rₙ* e^{i 2π n j / N}
 *
 * So we can denote:
 *
 *      DFTₙ = ∑ⱼ₌₀ᴺ⁻¹ Eⱼ' e^{-i 2π n j / N}    ;    IDFTⱼ = 1/N ∑ₙ₌₀ᴺ⁻¹ Ẽₙ' e^{i 2π n j / N}
 *
 * (Where Eⱼ' = Eⱼ·sⱼ and Ẽₙ' = Ẽₙ·rₙ*)
 *
 *
 * Therefore, we can use the fast Fourier transform to compute the coefficients, yielding the following expressions:
 *
 *      Ẽₙ = Δt·rₙ fft(Eⱼ·sⱼ)    ;    Eⱼ = 1/Δt·sⱼ* · ifft(Ẽₙ·rₙ*)
 *
 * And we don't have to worry about shifting the result.
 *
 * We should consider the Nyquist sampling theorem, which states that to avoid aliasing effects, we should discard all
 * frequencies higher than half the sampling frequency, given by fₘ = 1 / Δt. Thus, the frequency array should be
 * evenly spaced between -ωₘ/2 = -π/Δt and ωₘ/2 = π/Δt with Δω = 2π/(NΔt). If the frequency array does not meet this
 * relationship, we will have problems when switching between domains.
 */

#ifndef FOURIER_INCLUDED
#define FOURIER_INCLUDED

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>

/**
 * @brief Class for efficient computing of forward and backwards Fourier Transforms using FFTW.
 *
 * The idea is to provide an efficient way to compute a lot of forward and backward Fourier Transforms
 * for a vector defined on a fixed time and frequency grid with a given number of samples, N.
 *
 * FFTW generates a 'plan' instance that makes efficient computations of the fft for a fixed N,
 * so this class manages to mantain one of these plans to make efficient computations with the
 * matching phase factors given by the time and frequency arrays.
 *
 * If you do not care about speed or not going to do a lot of transforms, you can use DFT and IDFT
 * functions instead.
 */
class FourierTransform
{
private:
    int N;                                      // Number of samples
    double deltaT;                              // Spacing in the time grid
    double deltaOmega;                          // Spacing in the angular frequency grid
    std::vector<double> t;                      // Time vector
    std::vector<double> omega;                  // Angular frequency vector
    std::vector<std::complex<double>> r_n;      // Precomputed r_n phase factors
    std::vector<std::complex<double>> s_j;      // Precomputed s_j phase factors
    std::vector<std::complex<double>> r_n_conj; // Precomputed r_n phase factors conjugated
    std::vector<std::complex<double>> s_j_conj; // Precomputed s_j phase factors conjugated
    fftw_complex *in;                           // In array for the FFTW plan
    fftw_complex *out;                          // Out array for the FFTW plan
    fftw_plan forwardPlan;                      // FFTW plan for forward transform
    fftw_plan backwardPlan;                     // FFTW plan for backward transform
    std::vector<std::complex<double>> result;   // Array that stores the last result obtained

public:
    // Constructor
    FourierTransform(int nSamples, double dt, double t0)
    {
        this->N = nSamples;
        this->deltaT = dt;
        this->deltaOmega = 2 * M_PI / (N * deltaT);

        // Calculate the time and angular frequency grids
        t.resize(N);
        omega.resize(N);
        double start = -M_PI / deltaT;
        for (int i = 0; i < N; ++i)
        {
            t[i] = t0 + i * deltaT;
            omega[i] = start + (2 * M_PI * i / (N * deltaT));
        }

        // Compute r_n factors
        r_n.resize(N);
        r_n_conj.resize(N);
        if (t[0] == 0.0)
        {
            for (int i = 0; i < N; i++)
            {
                r_n[i] = deltaT;
                r_n_conj[i] = 1;
            }
        }
        else
        {
            double constFactor = t[0] * deltaOmega;
            for (int i = 0; i < N; i++)
            {
                r_n[i] = deltaT * std::exp(std::complex<double>(0, -i * constFactor));
                r_n_conj[i] = std::exp(std::complex<double>(0, i * constFactor));
            }
        }

        // Compute s_j factors
        s_j.resize(N);
        s_j_conj.resize(N);
        if (omega[0] == 0.0)
        {
            double constFactor = 1 / (N * deltaT);
            for (int i = 0; i < N; i++)
            {
                s_j[i] = 1;
                s_j_conj[i] = constFactor;
            }
        }
        else
        {
            double normFactor = 1 / (N * deltaT);
            for (int i = 0; i < N; i++)
            {
                s_j[i] = std::exp(std::complex<double>(0, -t[i] * omega[0]));
                s_j_conj[i] = std::exp(std::complex<double>(0, t[i] * omega[0])) * normFactor;
            }
        }

        // Initialize FFTW plans
        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
        forwardPlan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        backwardPlan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        result.resize(N);
    }

    // Destructor
    ~FourierTransform()
    {
        fftw_free(in);
        fftw_free(out);
        fftw_destroy_plan(forwardPlan);
        fftw_destroy_plan(backwardPlan);
    }

    // Perform forward transform
    std::vector<std::complex<double>> forwardTransform(const std::vector<std::complex<double>> &x)
    {
        std::vector<std::complex<double>> x_sj(N);

        // Apply s_j factors to input vector
        for (int i = 0; i < N; i++)
        {
            x_sj[i] = x[i] * s_j[i];
        }

        // Prepare the input data
        for (int i = 0; i < N; i++)
        {
            in[i][0] = std::real(x_sj[i]);
            in[i][1] = std::imag(x_sj[i]);
        }

        // Execute the forward transform
        fftw_execute_dft(forwardPlan, in, out);

        // Process the output with the necessary factors
        for (int i = 0; i < N; i++)
        {
            result[i] = r_n[i] * std::complex<double>(out[i][0], out[i][1]);
        }

        return result;
    }

    // Perform backward transform
    std::vector<std::complex<double>> backwardTransform(const std::vector<std::complex<double>> &x)
    {
        std::vector<std::complex<double>> x_rn(N);
        for (int i = 0; i < N; i++)
        {
            x_rn[i] = x[i] * r_n_conj[i];
        }
        // Prepare the input data
        for (int i = 0; i < N; ++i)
        {
            in[i][0] = std::real(x_rn[i]); // Real part of the input
            in[i][1] = std::imag(x_rn[i]); // Imaginary part of the input
        }

        // Execute the backward transform
        fftw_execute_dft(backwardPlan, in, out);

        // Process the output with the necessary factors
        for (int i = 0; i < N; i++)
        {
            result[i] = s_j_conj[i] * std::complex<double>(out[i][0], out[i][1]);
        }

        return result;
    }
};

/**
 * \brief Computes the Discrete Fourier Transform (DFT) of a dataset.
 *
 * This function calculates the DFT of a dataset represented by the vectors `x` and its timestamps `t`,
 * given the time step `deltaT`, angular frequency `omega`, and frequency step `deltaOmega`.
 *
 * \param x The vector representing the complex valued dataset.
 * \param t The vector of time values.
 * \param deltaT The time step.
 * \param omega The vector of angular frequencies.
 * \param deltaOmega The angular frequency step.
 * \return The vector representing the DFT of the dataset.
 */
std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>> &x, const std::vector<double> &t, double deltaT, const std::vector<double> &omega, double deltaOmega)
{

    int N = t.size();

    std::vector<std::complex<double>> r_n(N);
    std::vector<std::complex<double>> s_j(N);

    if (t[0] == 0.0)
    {
        for (int i = 0; i < N; i++)
        {
            r_n[i] = 1;
        }
    }
    else
    {

        for (int i = 0; i < N; i++)
        {
            r_n[i] = std::exp(std::complex<double>(0, -i * t[0] * deltaOmega));
        }
    }

    if (omega[0] == 0.0)
    {
        for (int i = 0; i < N; i++)
        {
            s_j[i] = 1;
        }
    }
    else
    {

        for (int i = 0; i < N; i++)
        {
            s_j[i] = std::exp(std::complex<double>(0, - t[i] * omega[0]));
        }
    }

    std::vector<std::complex<double>> x_sj(N);

    for (int i = 0; i < N; i++)
    {
        x_sj[i] = x[i] * s_j[i];
    }

    // Initialize FFTW and create the plan
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Prepare the input data
    for (int i = 0; i < N; ++i)
    {
        in[i][0] = std::real(x_sj[i]); // Real part of the input
        in[i][1] = std::imag(x_sj[i]); // Imaginary part of the input
    }

    // Execute the FFT
    fftw_execute(plan);

    // Process the output with the necessary factors
    std::vector<std::complex<double>> result(N);
    for (int i = 0; i < N; ++i)
    {
        result[i] = deltaT * r_n[i] * std::complex<double>(out[i][0], out[i][1]); // Extract the complex output
    }

    // Clean up and free resources
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
}

/**
 * \brief Computes the Inverse Discrete Fourier Transform (IDFT) of a dataset.
 *
 * This function calculates the IDFT of a dataset represented by the vectors `x` and its timestamps `t`,
 * given the time step `deltaT`, angular frequency `omega`, and frequency step `deltaOmega`.
 *
 * \param x The vector representing the complex valued dataset.
 * \param t The vector of time values.
 * \param deltaT The time step.
 * \param omega The vector of angular frequencies.
 * \param deltaOmega The angular frequency step.
 * \return The vector representing the IDFT of the dataset.
 */
std::vector<std::complex<double>> IDFT(const std::vector<std::complex<double>> &x, const std::vector<double> &t, double deltaT, const std::vector<double> &omega, double deltaOmega)
{

    int N = t.size();

    std::vector<std::complex<double>> r_n_conj(N);
    std::vector<std::complex<double>> s_j_conj(N);

    if (t[0] == 0.0)
    {
        for (int i = 0; i < N; i++)
        {
            r_n_conj[i] = 1;
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            r_n_conj[i] = std::exp(std::complex<double>(0, i * t[0] * deltaOmega));
        }
    }

    if (omega[0] == 0.0)
    {
        for (int i = 0; i < N; i++)
        {
            s_j_conj[i] = 1;
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            s_j_conj[i] = std::exp(std::complex<double>(0, t[i] * omega[0]));
        }
    }

    std::vector<std::complex<double>> x_rn(N);

    for (int i = 0; i < N; i++)
    {
        x_rn[i] = x[i] * r_n_conj[i];
    }

    // Initialize FFTW and create the plan
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Prepare the input data
    for (int i = 0; i < N; ++i)
    {
        in[i][0] = std::real(x_rn[i]); // Real part of the input
        in[i][1] = std::imag(x_rn[i]); // Imaginary part of the input
    }

    // Execute the FFT
    fftw_execute(plan);

    // Process the output with the necessary factors
    std::vector<std::complex<double>> result(N);
    for (int i = 0; i < N; ++i)
    {
        result[i] = s_j_conj[i] / (N * deltaT) * std::complex<double>(out[i][0], out[i][1]); // Extract the complex output
    }

    // Clean up and free resources
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
}


/**
 * Generate the frequency array for the discrete Fourier transform (DFT).
 * The function returns an equispaced frequency array corresponding to the given
 * number of samples and sample spacing.
 *
 * \param N The number of samples in the time domain.
 * \param deltaT The sample spacing (time interval) between consecutive samples.
 * \return A vector of frequencies in the frequency domain.
 */
std::vector<double> fftFreq(int N, double deltaT)
{
    std::vector<double> frequencies(N);
    double fSample = 1.0 / deltaT;
    double start = -0.5 * fSample;

    for (int i = 0; i < N; i++)
    {
        frequencies[i] = start + (i / (N * deltaT));
    }

    return frequencies;
}

#endif // FOURIER_HPP