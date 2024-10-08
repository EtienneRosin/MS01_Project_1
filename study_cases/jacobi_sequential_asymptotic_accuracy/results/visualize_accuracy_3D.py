import numpy as np
import matplotlib.pyplot as plt

# Charger les données à partir du fichier CSV
file_name = "study_cases/jacobi_sequential_asymptotic_accuracy/results/jacobi_sequential_asymptotic_accuracy.csv"
accuracy_data = np.genfromtxt(fname=file_name, dtype=float, delimiter=';', names=True)

# Récupérer les colonnes
x = accuracy_data['epsilon']
y = accuracy_data['h']
z = accuracy_data['error']
c = accuracy_data['last_iteration']

# Calculer h^2
h_squared = y ** 2  # h^2 pour chaque valeur de h

# Créer la figure et l'axe 3D
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# Créer le scatter plot pour les données de précision
accuracy_scatter = ax.scatter(x, y, z, c=c, label='Measured data', alpha=0.6)

# ax.plot(x, y, z)


lst_h = np.linspace(y.min(), y.max())
for epsilon_value in np.unique(x):
    # ax.plot(np.full_like(a = lst_h, fill_value=epsilon_value), lst_h, 1e5*lst_h**2, color = "black")
    ax.plot(np.full_like(a = lst_h, fill_value=epsilon_value), lst_h, 2e-3/(lst_h**2), color = "black")
    

# Définir les labels des axes
ax.set_xlabel(r"$\varepsilon$")
ax.set_ylabel(r"$h$")
ax.set_zlabel(r"$\left\|\boldsymbol{e}^{(l_{\text{stop}})}\right\|$")
ax.set_title(r"Analyse de l'Erreur en Fonction de $\varepsilon$, $h$ et $h^2$")
ax.grid(True)

# Ajouter une barre de couleur
fig.colorbar(accuracy_scatter, ax=ax, label='Itérations')

# Afficher la légende
ax.legend()

# Afficher le graphique
plt.show()