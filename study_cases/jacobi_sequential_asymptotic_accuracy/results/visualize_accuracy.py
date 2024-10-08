import numpy as np
import matplotlib.pyplot as plt
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
import scienceplots

plt.style.use('science')

file_name = "study_cases/jacobi_sequential_asymptotic_accuracy/results/jacobi_sequential_asymptotic_accuracy.csv"
file_name = "study_cases/jacobi_sequential_asymptotic_accuracy/results/jacobi_sequential_asymptotic_accuracy_dont_overwrite.csv"
accuracy_data = np.genfromtxt(fname=file_name, dtype=float, delimiter=';', names=True)

epsilon_values = accuracy_data['epsilon']
h_values = accuracy_data['h']
error_values = accuracy_data['error']
iterations = accuracy_data['last_iteration']

unique_epsilons = np.unique(epsilon_values)

fig = plt.figure()
ax = fig.add_subplot()

for epsilon_value in unique_epsilons:
    mask = epsilon_values == epsilon_value
    h_subset = h_values[mask]
    error_subset = error_values[mask]
    
    ax.plot(h_subset, error_subset, marker='o', label=fr"${epsilon_value:.1e}$", markersize = 3, linewidth = 0.5)
# \boldsymbol{e}^{(l_{})}
# stop l_{\text{stop}}
ax.set(xlabel = r"$h$", ylabel = r"$\left\| \mathbf{e} \right\|$")
ax.legend(title=r"$\epsilon$ values")
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.7)

# Ajuster les ticks des axes pour un format scientifique
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=10)

# ax.grid(True)
fig.savefig("jacobi_sequential_accuracy.pdf")
plt.show()

