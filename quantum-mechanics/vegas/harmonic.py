import sys
import numpy as np
import matplotlib.pyplot as plt

# This script will use the first data found in the first file provided.

file = open(sys.argv[1])
read = False

coord = []
wave_func = []
for line in file:
    if read :
        splitted = line.split()
        if(len(splitted) == 0) :
            break
        coord.append(float(splitted[0]))
        wave_func.append(float(splitted[1]))
    elif line == "Ground state wavefunction:\n" :
        read = True
        continue

sample_x = np.linspace(min(coord), max(coord), 4*len(coord))
def psi_0(x) :
    return np.exp(-x*x/2.)/np.pi**0.25
exact_wave_func = [psi_0(x) for x in sample_x]

plt.close()
plt.scatter(coord, wave_func, color = "blue",
            label = "Monte Carlo")
plt.scatter(sample_x, exact_wave_func, color = "lightblue", s = 5,
            label = "exact")
plt.xlabel("x")
plt.ylabel("ground state wavefunction")
plt.legend(loc = "upper right")
plt.show()




