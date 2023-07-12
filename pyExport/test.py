import fourier
import numpy as np
import matplotlib.pyplot as plt
import time

N = 65536
signalDuration = 1000
deltaT = signalDuration / N

t0 = -signalDuration / 4
t = t0 + np.arange(N) * deltaT

x = np.zeros(N, dtype=np.complex128)

for i in range(N):

    if t[i] != 0:
        sincValue = np.sin(t[i]) / t[i]
    else:
        sincValue = 1.0

    if t[i] >= -signalDuration/ 4 and t[i] <= signalDuration / 4:
        squareValue = 1.0
    else:
        squareValue = 0

    x[i] = sincValue + 1j * squareValue

ft = fourier.FourierTransform(N, deltaT, t0)
x = list(x)

start_time = time.time()

numberOfTimes = 1000
for i in range(numberOfTimes):
    xTransform = ft.forwardTransform(x)
    xInverse = ft.backwardTransform(xTransform)

# End the timer
end_time = time.time()

# Compute the elapsed time
elapsed_time = end_time - start_time

# Print the elapsed time
print("Elapsed time FFT:", (elapsed_time * 1e6) / numberOfTimes, "microseconds")

x = np.array(x)
xInverse = np.array(xInverse)

plt.plot(t, x.real)
plt.plot(t, x.imag)
plt.plot(t, xInverse.real, '.')
plt.plot(t, xInverse.imag, '.')

plt.show()

