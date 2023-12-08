/**
 * @file testRetrievers.cpp
 * @author Víctor Loras Herrero
 * @brief Tests to check that the COPRA retrieval algorithm is working OK
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

int main()
{
    int N = 256;
    double deltaT = 5e-18;

    double t0 = 0;

    FourierTransform ft(N, deltaT, t0);

    Pulse examplePulse(ft);
    double TBP = 2.33;
    examplePulse.randomPulse(TBP);
    std::vector<std::vector<double>> Tmeas = examplePulse.getTrace();

    // Add some noise to the trace
    std::vector<std::vector<double>> Tnoise = add_noise(Tmeas, N, 0.001);

    COPRA copraRetriever(ft, Tnoise);
    Pulse retrievedPulse = copraRetriever.retrieve(1e-16, 1500);
    std::vector<std::vector<double>> Tret = retrievedPulse.getTrace();

    // Save the trace for checking

    FILE *f;
    f = fopen("trace_check.txt", "wt");

    if (f == NULL)
    {
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(f, "%lf\t", Tnoise[i][j]);
        }
    }

    fclose(f);

    FILE *g;
    g = fopen("trace_check_retrieved.txt", "wt");

    if (g == NULL)
    {
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(g, "%lf\t", Tret[i][j]);
        }
    }

    fclose(g);

    return 0;
}