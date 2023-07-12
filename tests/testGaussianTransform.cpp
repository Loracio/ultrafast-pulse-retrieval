#include "../src/fourier.hpp"
#include "../src/utils.hpp"

#include <iostream>
#include <vector>
#include <complex>

std::vector<std::complex<double>> gaussianPulse(int N, std::vector<double> &t, double t0, double A, double omega_0, double phi_0, double tau)
{

    std::vector<std::complex<double>> result(N);

    for (int i = 0; i < N; i++)
    {
        result[i] = A * std::exp(-(t[i] - t0) * (t[i] - t0) / (2 * tau * tau)) * std::exp(std::complex<double>(0, omega_0 * t[i] + phi_0));
    }

    return result;
}


std::vector<std::complex<double>> gaussianTransform(int N, std::vector<double> &omega, double t0, double A, double omega_0, double phi_0, double tau)
{

    std::vector<std::complex<double>> result(N);

    for (int i = 0; i < N; i++)
    {
        result[i] = A * tau * sqrt(2 * M_PI) *  std::exp(-(omega[i] - omega_0) * (omega[i] - omega_0) / (2 * tau * tau)) * std::exp(std::complex<double>(0, - t0 * (omega[i] - omega_0) + phi_0));
    }

    return result;
}



int main()
{
    int N = 4096 * 2;
    double signalDuration = 10;
    double deltaT = signalDuration / N;

    double t0 = 5;
    double A = 1;
    double lambda_0 = 1.55;
    double omega_0 = 2 * M_PI * 3e8 * 1e-12 / (lambda_0 * 1e-6);
    double phi_0 = 0;
    double tau = 1;

    std::vector<double> t(N);

    for (int i = 0; i < N; i++)
    {
        t[i] = i * deltaT;
    }

    std::vector<double> omega = toAngularFrequency(fftFreq(N, deltaT));

    std::vector<std::complex<double>> pulse = gaussianPulse(N, t, t0, A, omega_0, phi_0, tau);
    std::vector<std::complex<double>> analyticTransform = gaussianTransform(N, omega, t0, A, omega_0, phi_0, tau);

    FourierTransform ft(N, deltaT, 0);

    std::vector<std::complex<double>> transformedPulse = ft.forwardTransform(pulse);
    std::vector<std::complex<double>> retrievedPulse = ft.backwardTransform(transformedPulse);

    FILE *f;
    f = fopen("gaussian_check.txt", "wt");

    if (f == NULL){
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t[i], std::real(pulse[i]), std::imag(pulse[i]), std::real(analyticTransform[i]), std::imag(analyticTransform[i]),std::real(transformedPulse[i]), std::imag(transformedPulse[i]), std::real(retrievedPulse[i]), std::imag(retrievedPulse[i]));
    }
    
    fclose(f);

    return 0;
}