import numpy as np
import matplotlib.pyplot as plt

file_name = "study_cases/jacobi_sequential_asymptotic_accuracy/results/jacobi_sequential_asymptotic_accuracy.csv"
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

ax.set(xlabel = r"$h$", ylabel = r"$\left\|\boldsymbol{e}^{(l_{\text{stop}})}\right\|$")
ax.legend(title=r"$\epsilon$ values")
# ax.set_xscale('log')
# ax.set_yscale('log')


ax.grid(True)
plt.show()



# import numpy as np
# import matplotlib.pyplot as plt

# # Charger les données à partir du fichier CSV
# file_name = "study_cases/jacobi_sequential_asymptotic_accuracy/results/jacobi_sequential_asymptotic_accuracy.csv"
# accuracy_data = np.genfromtxt(fname=file_name, dtype=float, delimiter=';', names=True)

# epsilon_values = accuracy_data['epsilon']
# h_values = accuracy_data['h']
# error_values = accuracy_data['error']
# iterations = accuracy_data['last_iteration']

# # Définir le nombre max d'itérations si non convergence (pour -1)
# max_iterations = 200000
# iterations_with_max = np.where(iterations == -1, max_iterations, iterations)

# # Normaliser les itérations pour obtenir des tailles de markers raisonnables
# min_marker_size = 1
# max_marker_size = 20
# marker_sizes = np.interp(iterations_with_max, (iterations_with_max.min(), iterations_with_max.max()), 
#                          (min_marker_size, max_marker_size))
# print(marker_sizes)
# # Extraire les valeurs uniques de epsilon
# unique_epsilons = np.unique(epsilon_values)

# # Créer la figure et l'axe
# fig, ax = plt.subplots()

# # Boucle sur les valeurs de epsilon pour tracer les courbes avec des marqueurs de taille variable
# for epsilon_value in unique_epsilons:
#     mask = epsilon_values == epsilon_value
#     h_subset = h_values[mask]
#     error_subset = error_values[mask]
#     marker_subset = marker_sizes[mask]
    
#     # Tracer les lignes entre les points
#     ax.plot(h_subset, error_subset, label=fr"$\epsilon = {epsilon_value:.1e}$", linewidth=0.5)
    
#     # Ajouter les marqueurs avec des tailles proportionnelles aux itérations
#     ax.scatter(h_subset, error_subset, s=marker_subset**2, alpha=0.6)

# # Ajouter les labels des axes
# ax.set(xlabel=r"$h$", ylabel=r"$\left\|\boldsymbol{e}^{(l_{\text{stop}})}\right\|$")

# # Première légende pour les courbes (epsilon)
# legend1 = ax.legend(loc='upper right', title=r"$\epsilon$ values")

# # Ajouter une seconde légende pour la taille des points
# # Créer des points invisibles pour la légende des tailles
# from matplotlib.lines import Line2D
# legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'{size:.0f} iterations',
#                           markerfacecolor='gray', markersize=size) 
#                    for size in [min_marker_size, (min_marker_size + max_marker_size) // 2, max_marker_size]]

# # Ajouter la seconde légende en dehors du graphe
# legend2 = fig.legend(handles=legend_elements, loc='upper center', title='Iterations (marker size)', ncol = len(legend_elements))

# # Remettre la première légende en place (car elle est remplacée par la deuxième sinon)
# ax.add_artist(legend1)

# # Activer la grille
# ax.grid(True)

# # Ajuster l'espace autour du graphique pour afficher la deuxième légende
# # plt.subplots_adjust(right=0.75)

# # Afficher le graphique
# plt.show()