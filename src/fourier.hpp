#ifndef FOURIER_INCLUDED
#define FOURIER_INCLUDED

/**
 * @file fourier.hpp
 * @author Víctor Loras Herrero
 * @brief Functions related to the numerical computing of the Fourier Transform using FFTW library (https://www.fftw.org/)
 * @version 0.1
 * @date 2023-07-02
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

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>

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

    std::vector<std::complex<double>> r_n;
    std::vector<std::complex<double>> s_j;

    if (t[0] == 0.0)
    {
        r_n = std::vector<std::complex<double>>(N, 1);
    }
    else
    {

        r_n.resize(N);

        for (int i = 0; i < N; i++)
        {
            r_n[i] = std::exp(std::complex<double>(0, -i * t[0] * deltaOmega));
        }
    }

    if (omega[0] == 0.0)
    {
        s_j = std::vector<std::complex<double>>(N, 1);
    }
    else
    {
        s_j.resize(N);

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

    std::vector<std::complex<double>> r_n_conj;
    std::vector<std::complex<double>> s_j_conj;

    if (t[0] == 0.0)
    {
        r_n_conj = std::vector<std::complex<double>>(N, 1);
    }
    else
    {

        r_n_conj.resize(N);

        for (int i = 0; i < N; i++)
        {
            r_n_conj[i] = std::exp(std::complex<double>(0, i * t[0] * deltaOmega));
        }
    }

    if (omega[0] == 0.0)
    {
        s_j_conj = std::vector<std::complex<double>>(N, 1);
    }
    else
    {
        s_j_conj.resize(N);

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
        result[i] = s_j_conj[i] / deltaT * std::complex<double>(out[i][0], out[i][1]); // Extract the complex output
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
