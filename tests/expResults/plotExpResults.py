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


def readRetrievedTrace(filename):
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


def readExpTrace(trace_filename):
    try:
        with open(trace_filename, 'r') as trace_file:
            matrix = []

            for line in trace_file:
                values_str = line.strip().split(',')

                # Convert string values to double and add to the row
                row = [float(value) for value in values_str]

                # Add the row to the matrix
                matrix.append(row)

        return np.array(matrix)

    except FileNotFoundError:
        print(f"File '{trace_filename}' not found.")
        return np.array([])


def readAxis(axis_filename):
    try:
        with open(axis_filename, 'r') as params_file:
            omegas = []
            delays = []

            for line in params_file:
                values_str = line.strip().split(',')

                if len(values_str) == 2:
                    # Convert strings to doubles
                    value1 = float(values_str[0])
                    value2 = float(values_str[1])

                    # Store values in vectors
                    omegas.append(value1)
                    delays.append(value2)
                else:
                    print(f"Error reading axis file. Line: {line}")

        return np.array(omegas), np.array(delays)

    except FileNotFoundError:
        print(f"File '{axis_filename}' not found.")
        return [], []


def meanValue(x, y):
    return np.sum(x * y) / np.sum(y)


def plotRetrievalResult(N, field, spectrum, Tmeas, Tretrieved, omegas, delays, centralWavelength, R):
    fig, ax = plt.subplots(2, 2)

    twin_ax00 = ax[0][0].twinx()
    twin_ax01 = ax[0][1].twinx()

    fieldIntensity = np.abs(field)**2

    spectralIntensity = np.abs(spectrum)**2

    fieldPhase = np.unwrap(np.angle(field))
    fieldPhase -= meanValue(fieldPhase, fieldIntensity)
    # fieldPhase = np.where(fieldIntensity < 1e-3, np.nan, fieldPhase) # Phase blanking

    spectralPhase = np.unwrap(np.angle(spectrum))
    spectralPhase -= meanValue(spectralPhase, spectralIntensity)
    # spectralPhase = np.where(spectralIntensity < 1e-3, np.nan, spectralPhase) # Phase blanking

    #! Roll temporal intensity (and its phase) so that its maximum lays in the middle
    shift_amount = N // 2 - np.argmax(fieldIntensity)


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

    # fig.legend(*ax[0][0].get_legend_handles_labels(),
            #    loc='upper center', ncols=4)

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

    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.ticklabel_format(useMathText=True)

    ax.grid()

    return fig, ax

if __name__ == '__main__':

    if len(sys.argv) != 9:
        print("Wrong args. Usage: N axisFilename traceFilename resultingFieldFilename resultingSpectrumFilename resultingTraceFilename resultingErrorsFilename centralWavelength")
        sys.exit(1)

    N = int(sys.argv[1])
    axisFilename = sys.argv[2]
    traceFilename = sys.argv[3]
    resultingFieldFilename = sys.argv[4]
    resultingSpectrumFilename = sys.argv[5]
    resultingTraceFilename = sys.argv[6]
    resultingErrorFilename = sys.argv[7]
    centralWavelength = float(sys.argv[8])

    measuredOmegas, measuredDelays = readAxis(axisFilename)
    originalTrace = readExpTrace(traceFilename)
    retrievedField = readResult(resultingFieldFilename)
    retrievedSpectrum = readResult(resultingSpectrumFilename)
    retrievedTrace = readRetrievedTrace(resultingTraceFilename)
    retrievedErrors = readErrors(resultingErrorFilename)


    figError, axError = plotErrorEvolution(retrievedErrors[:-1]) #! All except the last one, which is the best retrieved error

    fig, ax = plotRetrievalResult(N, retrievedField, retrievedSpectrum,
                                  originalTrace, retrievedTrace, measuredOmegas, measuredDelays, centralWavelength, retrievedErrors[-1])
    

    plt.show()
