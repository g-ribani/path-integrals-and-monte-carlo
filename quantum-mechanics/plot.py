from matplotlib import pyplot as plt
import linecache as lnch

file_name = 'naive.out'

results = []
errors = []
times = []

first_line = 180       # lines containing data to extract
last_line = 263        #
exact_amp = 0.206913    # todo: automize these two read ops
simulation_details = "n_steps = 8"  #
for n_line in range(first_line, last_line + 1):
    alpha_data = (lnch.getline(file_name, n_line)).split()
    data = [float(d) for d in alpha_data]
    results.append(data[0])
    errors.append(data[0]*data[1])
for k in range(first_line, last_line + 1) :
    times.append(k)

plt.errorbar(times, results, yerr = errors, ecolor = 'lightblue', color='blue',
                label="naive MC, " + simulation_details)
plt.xlabel("1 unit = 1 sec of computation (arbitrary offset)")
plt.hlines(y = exact_amp, xmin = first_line, xmax = last_line, color = 'red',
            label = "exact result")
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc="upper left")
plt.show()




