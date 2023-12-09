#!/bin/bash

# Configuration parameters
N=256                 # Number of samples
deltaT=6.12           # Time spacing in femtoseconds
centralWaveLength=790 # Central wavelength in nanometers
TBP=2.75              # Time-Bandwidth product of the random pulse to generate
noiseLevel=0.0        # Noise level added to the random pulse's trace

# Retriever parameters
retriever="COPRA"     # Selects a switch between 'COPRA', 'GPA'
maximumIterations=300 # Maximum number of iterations for the retrieval
tolerance=1e-12       # Tolerance of the trace error (R)

# Initial candidate selection. Leave for blank for a random pulse with TBP = 0.51
initialCandidateField="none"    # Filename of the field for the initial candidate
initialCandidateSpectrum="none" # Filename of the spectrum for the initial candidate

# Filename of the original random pulse
originalFieldFilename="./tests/originalField.txt"       # Filename of the field to retrieve
originalSpectrumFilename="./tests/originalSpectrum.txt" # Filename of the field to retrieve
originalTraceFilename="./tests/originalTrace.txt"       # Numerical values of the trace (N rows with N columns)

# Filenames of the result
resultingFieldFilename="./tests/retrievedField.txt"
resultingSpectrumFilename="./tests/retrievedSpectrum.txt"
resultingTraceFilename="./tests/retrievedTrace.txt"
resultingErrorsFilename="./tests/retrievedErrors.txt"

#######################################################################################################################################################################################################################################

# Launching the compiled program
./testRetrievers.o $N $deltaT $TBP $noiseLevel $originalFieldFilename $originalSpectrumFilename $originalTraceFilename $retriever $maximumIterations $tolerance $initialCandidateField $initialCandidateSpectrum $resultingFieldFilename $resultingSpectrumFilename $resultingTraceFilename $resultingErrorsFilename

#! Depends on where the python instalation is in your system
~/.conda/envs/defenv/bin/python ./tests/plotRetrievalResults.py $N $originalFieldFilename $originalSpectrumFilename $originalTraceFilename $resultingFieldFilename $resultingSpectrumFilename $resultingTraceFilename $resultingErrorsFilename $centralWaveLength $deltaT
