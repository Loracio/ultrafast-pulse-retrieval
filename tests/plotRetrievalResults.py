import matplotlib.pyplot as plt
import numpy as np
import sys

plt.rcParams.update({'font.size': 20})  # Tamaño de la fuente del plot


def readResult(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        result = []
        for line in lines:
            real_part, imag_part = map(float, line.strip().split('\t'))
            result.append(complex(real_part, imag_part))

        return np.array(result)

    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        return np.array([])


def readTrace(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        matrix = []
        for line in lines:
            row = list(map(float, line.strip().split('\t')))
            matrix.append(row)

        return np.array(matrix)

    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        return np.array([])

def readErrors(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        result = []
        for line in lines:
            result.append(float(line))

        return np.array(result)

    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        return np.array([])

def meanValue(x, y):
    return np.sum(x * y) / np.sum(y)

def plotRetrievalResult(N, originalField, originalSpectrum, field, spectrum, Tmeas, Tretrieved, omegas, delays, centralWavelength, R):
    fig, ax = plt.subplots(2, 2)

    twin_ax00 = ax[0][0].twinx()
    twin_ax01 = ax[0][1].twinx()

    fieldIntensity = np.abs(field)**2
    fieldIntensity /= np.max(fieldIntensity)

    originalFieldIntensity = np.abs(originalField)**2
    originalFieldIntensity /= np.max(originalFieldIntensity)

    spectralIntensity = np.abs(spectrum)**2
    spectralIntensity /= np.max(spectralIntensity)
    originalSpectralIntensity = np.abs(originalSpectrum)**2
    originalSpectralIntensity /= np.max(originalSpectralIntensity)

    fieldPhase = np.unwrap(np.angle(field))
    fieldPhase -= meanValue(fieldPhase, fieldIntensity)
    # fieldPhase = np.where(fieldIntensity < 1e-3, np.nan, fieldPhase) # Phase blanking

    originalFieldPhase = np.unwrap(np.angle(originalField))
    originalFieldPhase -= meanValue(originalFieldPhase, originalFieldIntensity)
    # originalFieldPhase = np.where(originalFieldIntensity < 1e-3, np.nan, originalFieldPhase) # Phase blanking

    spectralPhase = np.unwrap(np.angle(spectrum))
    spectralPhase -= meanValue(spectralPhase, spectralIntensity)
    # spectralPhase = np.where(spectralIntensity < 1e-3, np.nan, spectralPhase) # Phase blanking

    originalSpectralPhase = np.unwrap(np.angle(spectrum))
    originalSpectralPhase -= meanValue(originalSpectralPhase, originalSpectralIntensity)
    # originalSpectralPhase = np.where(originalSpectralIntensity < 1e-3, np.nan, originalSpectralPhase) # Phase blanking

    #! Roll temporal intensity (and its phase) so that its maximum lays in the middle
    shift_amount = N // 2 - np.argmax(fieldIntensity)

    ax[0][0].plot(delays, originalFieldIntensity, linewidth=2, color='orange',
                  label='Original intensity')
    twin_ax00.plot(delays, originalFieldPhase, '-.', linewidth=2, color='violet')
    ax[0][0].plot(np.nan, '-.', label='Original phase', color='violet')

    ax[0][0].plot(delays, np.roll(fieldIntensity, shift_amount), color='blue',
                  label='Retrieved intensity')
    twin_ax00.plot(delays, np.roll(fieldPhase, shift_amount), '-.', color='red')
    ax[0][0].plot(np.nan, '-.', label='Retrieved phase', color='red')

    ax[0][0].set_xlabel("Time (fs)")
    ax[0][0].set_ylabel("Intensity (a.u.)")
    twin_ax00.set_ylabel("Phase (rad)")
    ax[0][0].set_title("Time domain")
    ax[0][0].grid()

    ax[0][0].yaxis.label.set_color('blue')
    ax[0][0].tick_params(axis='y', colors='blue')
    ax[0][0].spines['left'].set_color('blue')
    twin_ax00.spines['right'].set_color('red')
    twin_ax00.tick_params(axis='y', colors='red')
    twin_ax00.yaxis.label.set_color('red')

    ax[0][1].plot(2 * np.pi * 300 / (2 * np.pi * 300 / centralWavelength + omegas),
                  originalSpectralIntensity, linewidth=2, color='orange', label='Original spectral intensity')
    twin_ax01.plot(2 * np.pi * 300 / (2 * np.pi * 300 /
                   centralWavelength + omegas), originalSpectralPhase, '-.', color='violet')
    ax[0][1].plot(np.nan, '-.', linewidth=2, label='Original spectral phase', color='violet')


    ax[0][1].plot(2 * np.pi * 300 / (2 * np.pi * 300 / centralWavelength + omegas),
                  spectralIntensity, color='blue', label='Retrieved spectral intensity')
    twin_ax01.plot(2 * np.pi * 300 / (2 * np.pi * 300 /
                   centralWavelength + omegas), spectralPhase, '-.', color='red')
    ax[0][1].plot(np.nan, '-.', label='Retrieved spectral phase', color='red')

    ax[0][1].set_xlabel("λ (nm)")
    ax[0][1].set_ylabel("Intensity (a.u.)")
    twin_ax01.set_ylabel("Phase (rad)")
    ax[0][1].set_title("Frequency domain")
    ax[0][1].grid()
    ax[0][1].yaxis.label.set_color('blue')
    ax[0][1].tick_params(axis='y', colors='blue')
    ax[0][1].spines['left'].set_color('blue')
    twin_ax01.spines['right'].set_color('red')
    twin_ax01.tick_params(axis='y', colors='red')
    twin_ax01.yaxis.label.set_color('red')

    fig.legend(*ax[0][0].get_legend_handles_labels(),
               loc='upper center', ncols=4)

    TmeasNormalized = Tmeas / np.max(Tmeas)

    TretrievedNormalized = Tretrieved / np.max(Tretrieved)

    im0 = ax[1][0].pcolormesh(2 * np.pi * 300 / centralWavelength +
                              omegas, delays, TmeasNormalized, cmap='YlGnBu_r')
    fig.colorbar(im0, ax=ax[1][0])
    ax[1][0].set_xlabel("ω (2π/fs)")
    ax[1][0].set_ylabel("τ (fs)")
    ax[1][0].set_title("Measured trace")

    im1 = ax[1][1].pcolormesh(2 * np.pi * 300 / centralWavelength +
                              omegas, delays, TretrievedNormalized, cmap='YlGnBu_r')
    fig.colorbar(im1, ax=ax[1][1])
    ax[1][1].set_xlabel("ω (2π/fs)")
    ax[1][1].set_ylabel("τ (fs)")

    def format_scientific_notation(value, precision=2):
        exponent = int(np.floor(np.log10(abs(value))))
        coefficient = value / 10**exponent
        formatted_str = f"${coefficient:.{precision}f} \\times 10^{{{exponent}}}$"
        return formatted_str
        

    ax[1][1].set_title(fr"Retrieved trace, R = {format_scientific_notation(R)}")

    return fig, ax

def plotErrorEvolution(errors):
    fig, ax = plt.subplots()
    ax.plot(errors, linewidth=2)
    ax.set_xlabel("Iteration number")
    ax.set_ylabel("Trace error (R)")
    ax.set_title("Trace error evolution")

    plt.ticklabel_format(axis="y", style="sci")
    ax.ticklabel_format(useMathText=True)

    ax.grid()

    return fig, ax

if __name__ == '__main__':

    if len(sys.argv) != 11:
        print("Wrong args. Usage: N originalFieldFilename originalSpectrumFilename originalTraceFilename resultingFieldFilename resultingSpectrumFilename resultingTraceFilename resultingErrorsFilename centralWavelength")
        sys.exit(1)

    N = int(sys.argv[1])

    originalFieldFilename = sys.argv[2]
    originalSpectrumFilename = sys.argv[3]
    originalTraceFilename = sys.argv[4]

    resultingFieldFilename = sys.argv[5]
    resultingSpectrumFilename = sys.argv[6]
    resultingTraceFilename = sys.argv[7]
    resultingErrorFilename = sys.argv[8]

    centralWavelength = float(sys.argv[9])
    deltaT = float(sys.argv[10])

    deltaOmega = 2 * np.pi / (N * deltaT)
    measuredOmegas = np.array([-np.pi / deltaT + i * deltaOmega for i in range(N)])
    measuredDelays = np.array([(i - np.floor(0.5 * N)) * deltaT for i in range(N)])

    originalField = readResult(originalFieldFilename)    
    originalSpectrum = readResult(originalSpectrumFilename)    
    originalTrace = readTrace(originalTraceFilename)

    retrievedField = readResult(resultingFieldFilename)
    retrievedSpectrum = readResult(resultingSpectrumFilename)
    retrievedTrace = readTrace(resultingTraceFilename)
    retrievedErrors = readErrors(resultingErrorFilename)


    figError, axError = plotErrorEvolution(retrievedErrors[:-1]) #! All except the last one, which is the best retrieved error

    fig, ax = plotRetrievalResult(N, originalField, originalSpectrum, retrievedField, retrievedSpectrum,
                                  originalTrace, retrievedTrace, measuredOmegas, measuredDelays, centralWavelength, retrievedErrors[-1])
    

    plt.show()