import sys
import numpy as np
import matplotlib.pyplot as plt

# This script will use the first data found in the first file provided.

file = open(sys.argv[1])
read = False

errors = []
bin_sizes = []
n_conf = 0
for line in file:
    splitted = line.split()
    if read :
        if(len(splitted) == 0) :
            read = False
            continue
        errors.append(float(splitted[0]))
        bin_sizes.append(float(splitted[2])/n_conf)
    elif len(splitted) == 0 :
        continue
    elif splitted[0] == "n_conf" :
        n_conf = float(splitted[2])
        continue
    elif line == "excitation energy errors ( bin_size ):\n" :
        read = True
        continue

plt.close()
plt.figure(num = file.name)
plt.scatter(bin_sizes, errors, color = "blue", s = 20)
plt.xlabel("relative bin size")
plt.ylabel("relative error")
plt.show()




