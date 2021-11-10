from matplotlib import pyplot as plt
import linecache as lnch
import numpy as np

file_name = 'main.log'

results = []
errors = []
evals = []

first_line = 424       # lines containing data to extract
last_line = 628        #
exact_amp = 2.069131e-01    # todo: automize these the following read ops
n_steps = 12                 #
n_evals = 2e6               #
simulation_details = "n_steps = " + str(n_steps)  #
for n_line in range(first_line, last_line + 1) :
    if ((n_line - first_line) % 3 == 0):
        alpha_data = (lnch.getline(file_name, n_line)).split()
        results.append(float(alpha_data[4]))
        errors.append(float(alpha_data[6]))
        evals.append(n_line)

average = np.mean(results)
error = np.std(results)/np.sqrt(len(results))
print("exact =", exact_amp)
print("estimate =", average, " +/- ", error)

plt.errorbar(evals, results, yerr = errors, ecolor = 'lightblue', color='blue',
                label= ("crude MC, " + simulation_details))
plt.xlabel("1 unit = " + str(n_evals) +
            " function evaluations (arbitrary offset()")
plt.hlines(y = exact_amp, xmin = first_line, xmax = last_line, color = 'red',
            label = "exact result")
plt.hlines(y = average, xmin = first_line, xmax = last_line,
            color = 'green', label = "average result")
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc="upper left")
plt.show()




