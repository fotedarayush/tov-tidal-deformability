import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('/Volumes/AYUSH DRIVE/NS Hyperon project/tov-tidal-deformability/results/EoS_01_06_Mg800.dat')

rho_geom = data[:, 0]
rho_c = rho_geom * 6.176e17

R = data[:, 1]
M = data[:, 2]

compactness = data [:, 6]
k2 = data[:, 9]
kappa2 = data[:, 10]
Lambda = data[:, 11]

plt.figure()

plt.plot(R, M)

plt.xlabel(r"$R$ [km]")
plt.ylabel(r"$M/M_\odot$")
plt.title("Mass-Radius Relation")



plt.grid(True)

plt.tight_layout()

# Save the plot
plt.savefig("mass_radius.png", dpi=300)

# Try to display it
plt.show()

plt.plot(M, Lambda)

plt.xlabel(r"$M/M_\odot$")
plt.ylabel(r"$\Lambda$")
plt.yscale("log")

plt.savefig("tidal_deformability.png", dpi=300)



plt.show()