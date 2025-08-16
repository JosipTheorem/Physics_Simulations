# Load band structure data from QE output
import numpy as np
import matplotlib.pyplot as plt


path = r"F:\GitHub\Praksa_grafen\bands_1x1\+\graphene.band.gnu" ##### path to graphene.band.gnu

data = np.loadtxt(path) 
# Extract k-points and energy values
k_values11 = data[:, 0]  # First column: k-points
energies11 = data[:, 1:]  # Remaining columns: Energy levels

# Plot the bands
plt.figure(figsize=(8,6))

plt.scatter(k_values11, energies11[:, 0], s=5, color='black')



plt.xlabel("k")
plt.ylabel("E(eV)")
plt.title("Graphene 1x1")
#plt.savefig("graphene_band_1x1.png", dpi=300)

plt.show()