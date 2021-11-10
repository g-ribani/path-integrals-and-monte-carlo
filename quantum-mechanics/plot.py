from matplotlib import pyplot as plt
import linecache as lnch
import numpy as np

file_name = 'naive-quater.log'

results = []
errors = []
times = []

first_line = 65       # lines containing data to extract
last_line = 326       #
exact_amp = 5.662614e-02    # todo: automize these two read ops
n_steps = 10
simulation_details = "n_steps = " + str(n_steps)  #
for n_line in range(first_line, last_line + 1):
    alpha_data = (lnch.getline(file_name, n_line)).split()
    data = [float(d) for d in alpha_data]
    results.append(data[0])
    errors.append(data[0]*data[1])
    times.append(n_line)
mean = np.mean(results)
print("mean =", mean)

plt.errorbar(times, results, yerr = errors, ecolor = 'lightblue', color='blue',
                label= "naive MC, " + simulation_details)
plt.xlabel("1 unit = 1 sec of computation (arbitrary offset)")
plt.hlines(y = exact_amp, xmin = first_line, xmax = last_line, color = 'red',
            label = "exact result")
plt.hlines(y = mean, xmin = first_line, xmax = last_line, color = 'green',
            label = "mean result")
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc = "upper left")
plt.show()




