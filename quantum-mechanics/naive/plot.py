from matplotlib import pyplot as plt
import linecache as lnch
import numpy as np

file_name = 'naive-5.log'

results = []
errors = []
times = []

first_line = 10       # lines containing data to extract
last_line = 861       #
exact_amp = 0.206913    # todo: automize this read ops
for n_line in range(first_line, last_line + 1):
    alpha_data = (lnch.getline(file_name, n_line)).split()
    data = [float(d) for d in alpha_data]
    results.append(data[0])
    errors.append(data[0]*data[1])
    times.append(n_line)

best = results[len(results)-1]

plt.errorbar(times, results, yerr = errors, ecolor = 'lightblue', color='blue',
                label= "naive MC")
plt.xlabel("1 unit = 1 sec of computation (arbitrary offset)")
plt.hlines(y = exact_amp, xmin = first_line, xmax = last_line, color = 'red',
            label = "exact =" + str(exact_amp))
plt.hlines(y = best, xmin = first_line, xmax = last_line, color = 'green',
            label = "best estimate = " + str(best))
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc = "upper left")
plt.show()




