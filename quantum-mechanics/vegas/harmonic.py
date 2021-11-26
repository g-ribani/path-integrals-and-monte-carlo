import sys
import numpy as np
import matplotlib.pyplot as plt

# This script will use the first data found in the first file provided.

file = open(sys.argv[1])
read = False

coord = []
wave_func = []
for line in file:
    splitted = line.split()
    if read :
        if(len(splitted) == 0) :
            break
        coord.append(float(splitted[0]))
        wave_func.append(float(splitted[1]))
    elif len(splitted) == 0 :
        continue
    elif splitted[0] == "mass" :
        mass = float(splitted[2])
    elif splitted[0] == "frequency" :
        freq = float(splitted[2])
    elif line == "Ground state wavefunction:\n" :
        read = True
        continue

inter_points = 4
sample_x = np.linspace(min(coord), max(coord),
                       (inter_points + 1)*len(coord) - inter_points)
def psi_0(x) :
    return (mass*freq/np.pi)**0.25*np.exp(-mass*freq*x*x/2.)
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




