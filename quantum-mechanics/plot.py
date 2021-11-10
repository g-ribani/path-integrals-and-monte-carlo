from matplotlib import pyplot as plt
import linecache as lnch
import numpy as np

file_name = 'naive-bis.out'

results = []
errors = []
times = []

first_line = 10       # lines containing data to extract
last_line = 81        #
exact_amp = 0.206913    # todo: automize these two read ops
n_steps = 9
simulation_details = "n_steps = " + str(n_steps)  #
for n_line in range(first_line, last_line + 1):
    alpha_data = (lnch.getline(file_name, n_line)).split()
    data = [float(d) for d in alpha_data]
    results.append(data[0])
    errors.append(data[0]*data[1])
    times.append(n_line)

plt.errorbar(times, results, yerr = errors, ecolor = 'lightblue', color='blue',
                label= "naive MC, " + simulation_details)
plt.xlabel("1 unit = 1 sec of computation (arbitrary offset)")
plt.hlines(y = exact_amp, xmin = first_line, xmax = last_line, color = 'red',
            label = "exact result")
plt.hlines(y = np.mean(results), xmin = first_line, xmax = last_line,
            color = 'green', label = "average result")
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc = "upper left")
plt.show()




