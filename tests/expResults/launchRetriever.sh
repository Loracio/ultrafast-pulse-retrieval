#!/bin/bash

# Configuration parameters and measurements
N=128                                                                # Number of samples
centralWaveLength=790                                                # Central wavelength in nanometers
axisFilename="./tests/expResults/axis_compressed_pulse.csv"          # Filename of the axis of the trace (1st column: angular frequencies, 2nd column: delays)
traceFilename="./tests/expResults/compressed_pulse_900mA_trace.csv" # Numerical values of the trace (N rows with N columns)

# Retriever parameters
retriever="COPRA"     # Selects a switch between 'COPRA', 'GPA' or 'PIE'
maximumIterations=500 # Maximum number of iterations for the retrieval
tolerance=1e-16       # Tolerance of the trace error (R)

# Initial candidate selection. Leave for blank for a random pulse with TBP = 0.51
initialCandidateField="none"    # Filename of the field for the initial candidate
initialCandidateSpectrum="none" # Filename of the spectrum for the initial candidate

# Filenames of the result
resultingFieldFilename="./tests/expResults/retrievedField.txt"
resultingSpectrumFilename="./tests/expResults/retrievedSpectrum.txt"
resultingTraceFilename="./tests/expResults/retrievedTrace.txt"
resultingErrorsFilename="./tests/expResults/retrievedErrors.txt"

#######################################################################################################################################################################################################################################

# Launching the compiled program
./retrieveExpResults.o $N $axisFilename $traceFilename $retriever $maximumIterations $tolerance $initialCandidateField $initialCandidateSpectrum $resultingFieldFilename $resultingSpectrumFilename $resultingTraceFilename $resultingErrorsFilename

# If the cpp program does not return 0 something went wrong
if [ $? -ne 0 ]; then
    exit 1
fi

#! Depends on where the python instalation is in your system
~/.conda/envs/defenv/bin/python ./tests/expResults/plotExpResults.py $N $axisFilename $traceFilename $resultingFieldFilename $resultingSpectrumFilename $resultingTraceFilename $resultingErrorsFilename $centralWaveLength
