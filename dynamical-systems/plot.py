from matplotlib import pyplot as plt
import linecache as lnch
import numpy as np

results = []
errors = []
evals = []

crude_file_name = 'crude.log'
crude_first_line = 9       # lines containing data to extract
crude_last_line = 129       # lines containing data to extract

vegas_file_name = 'vegas.log'
vegas_first_line = 37        #
vegas_last_line = 157        #

exact_amp = 0.1908675    # todo: automize these the following read ops
n_evals = 1e7               #

max_n = max(crude_last_line - crude_first_line,
            vegas_last_line - vegas_first_line)
plt.hlines(y = exact_amp, xmin = 0, xmax = max_n, color = 'red',
            label = "exact =" + str(round(exact_amp,6)))

# plotting crude data
for n_line in range(crude_first_line, crude_last_line + 1) :
    if ((n_line - crude_first_line) % 2 == 0):
        alpha_data = (lnch.getline(crude_file_name, n_line)).split()
        results.append(float(alpha_data[4]))
        errors.append(float(alpha_data[6]))
        evals.append(n_line - crude_first_line)

average = np.mean(results)
error = np.std(results)/np.sqrt(len(results))

plt.errorbar(evals, results, yerr = errors, ecolor = 'lightblue',
                color='blue', label= "crude MC")
plt.hlines( y = average, xmin = 0, xmax = max_n, color = 'green',
            label = "crude average =" + str(round(average,4))
                                + " +/- " + str(round(error,4)) )

# plotting vegas data
results = []
errors = []
evals = []

for n_line in range(vegas_first_line, vegas_last_line + 1) :
    if ((n_line - vegas_first_line) % 2 == 0):
        alpha_data = (lnch.getline(vegas_file_name, n_line)).split()
        results.append(float(alpha_data[4]))
        errors.append(float(alpha_data[6]))
        evals.append(n_line - vegas_first_line)

average = np.mean(results)
error = np.std(results)/np.sqrt(len(results))

plt.errorbar(evals, results, yerr = errors, ecolor = 'moccasin',
                color='orange', label= "Vegas MC")
plt.hlines( y = average, xmin = 0, xmax = 0, color = 'orange',
            label = "Vegas average =" + str(round(average,6))
                                + " +/- " + str(round(error,6)) )

plt.xlabel("1 unit = " + str(int(n_evals)) +
            " function evaluations (arbitrary offset)")
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc="upper right")
plt.show()




