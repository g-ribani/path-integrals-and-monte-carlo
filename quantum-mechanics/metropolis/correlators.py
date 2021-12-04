import sys
import numpy as np
import matplotlib.pyplot as plt

# This script will use the first data found in the first file provided.

file = open(sys.argv[1])
read = False

intervals = []
energies = []
errors = []
for line in file:
    splitted = line.split()
    if read :
        if(len(splitted) == 0) :
            read = False
            continue
        intervals.append(float(splitted[0]))
        energies.append(float(splitted[4]))
        #^ corrective factor = 1.-4./3.*np.exp(-2.*frequency*time_step*float(splitted[0]))
        errors.append(float(splitted[6]))
    elif len(splitted) == 0 :
        continue
    elif splitted[0] == "frequency" :
        frequency = float(splitted[2])
    elif splitted[0] == "time_step" :
        time_step = float(splitted[2])
    elif line == "excitation energy estimates:\n" :
        read = True
        continue
    elif splitted[0] == "best_estimate" :
        best_estimate = float(splitted[2])
        break
times = [k*time_step for k in intervals]

plt.close()
plt.figure(num = file.name)
plt.errorbar(times, energies, yerr = errors, fmt = 'o', markersize = 3.,
             color = "blue", ecolor = "lightblue", label= "Metropolis")
plt.hlines(y = frequency,
           xmin = min(times)-time_step/2., xmax = max(times)+time_step/2.,
           color = "green", label = "exact" )
plt.hlines(y = best_estimate,
           xmin = min(times)-time_step/2., xmax = max(times)+time_step/2.,
           color = "red", label = "best estimate" )
plt.xlabel("propagation time")
plt.ylabel("excitation energy estimate")
plt.legend(loc = "best")
plt.show()




