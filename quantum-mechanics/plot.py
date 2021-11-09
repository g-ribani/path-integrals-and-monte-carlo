from matplotlib import pyplot as plt
import linecache as lnch

file_name = 'naive.out'

results = []
errors = []
times = []

first_line = 340       # the first line of the file containing data
last_line = 420        # the last line
for n_line in range(first_line, last_line + 1):
    alpha_data = (lnch.getline(file_name, n_line)).split()
    data = [float(d) for d in alpha_data]
    results.append(data[0])
    errors.append(data[0]*data[1])
for k in range(first_line, last_line + 1) :
    times.append(k)

plt.errorbar(times, results, yerr = errors, ecolor = 'lightblue', color='blue', label="naive MC computation")
plt.xlabel("1 unit = 1 sec of computation (arbitrary offset)")
plt.hlines(y = 0.206, xmin = first_line, xmax = last_line, color = 'red', label = "exact result")   # exact result
plt.ylabel("amplitude (hbar = 1)")
plt.legend(loc="upper left")
plt.show()




