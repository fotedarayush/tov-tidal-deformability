
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# Load current EoS result
# ============================================================
filename = '/Volumes/AYUSH DRIVE/NS Hyperon project/tov-tidal-deformability/results/EoS_01_06_Mg800.dat'
df = filename
stem = Path(df).stem
data = np.loadtxt(df)

# Output columns
rho_geom = data[:, 0]
R = data[:, 1]          # already in km in current Fortran version
M = data[:, 2]          # solar masses
k2 = data[:, 9]
Lambda = data[:, 11]

# Convert central density to g/cm^3
rho_c = rho_geom * 6.176e17


# ============================================================
# Maximum mass
# ============================================================

imax = np.argmax(M)

Mmax = M[imax]
R_at_Mmax = R[imax]
rho_at_Mmax = rho_c[imax]

print("===================================")
print("Stellar properties")
print("===================================")

print(f"Maximum mass     = {Mmax:.4f} M_sun")
print(f"Radius at Mmax   = {R_at_Mmax:.4f} km")
print(f"Central density  = {rho_at_Mmax:.4e} g/cm^3")


# ============================================================
# Select useful stable branch
# ============================================================

# Only use models up to the maximum-mass configuration
M_stable = M[:imax + 1]
R_stable = R[:imax + 1]
L_stable = Lambda[:imax + 1]

# Remove pathological low-mass / negative tidal points
mask = (M_stable > 0.5) & (L_stable > 0)

M_phys = M_stable[mask]
R_phys = R_stable[mask]
L_phys = L_stable[mask]

# Sort by mass so interpolation works correctly
order = np.argsort(M_phys)

M_phys = M_phys[order]
R_phys = R_phys[order]
L_phys = L_phys[order]


# ============================================================
# Properties at 1.4 solar masses
# ============================================================

if M_phys.min() <= 1.4 <= M_phys.max():

    R14 = np.interp(1.4, M_phys, R_phys)
    Lambda14 = np.interp(1.4, M_phys, L_phys)

    print(f"R_1.4            = {R14:.4f} km")
    print(f"Lambda_1.4       = {Lambda14:.4f}")

else:
    print("1.4 M_sun is outside the usable mass range.")


# ============================================================
# Mass-radius plot
# ============================================================

# Find maximum-mass model
imax = np.argmax(M)

# Keep only the stable branch up to Mmax
M_stable = M[:imax + 1]
R_stable = R[:imax + 1]

# Remove very-low-mass models from the plot
mask_mr = M_stable > 0.5

plt.figure(figsize=(7, 5))

plt.plot(
    R_stable[mask_mr],
    M_stable[mask_mr]
)

plt.scatter(
    R_at_Mmax,
    Mmax,
    label=fr"$M_{{\rm max}}={Mmax:.2f}\,M_\odot$"
)

plt.xlabel(r"$R$ [km]")
plt.ylabel(r"$M/M_\odot$")
plt.title("Mass-Radius Relation")

plt.xlim(11, 14)
plt.ylim(0.5, Mmax + 0.1)

plt.grid(True)
plt.legend()

plt.tight_layout()

plt.savefig(
    f"results/{stem}_mass_radius.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# ============================================================
# Tidal deformability plot
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    M_phys,
    L_phys
)

if M_phys.min() <= 1.4 <= M_phys.max():

    plt.scatter(
        1.4,
        Lambda14,
        label=fr"$\Lambda_{{1.4}}={Lambda14:.0f}$"
    )

plt.xlabel(r"$M/M_\odot$")
plt.ylabel(r"$\Lambda$")
plt.title("Tidal Deformability")

plt.yscale("log")

plt.grid(True)
plt.legend()

plt.tight_layout()

plt.savefig(
    f"results/{stem}_tidal_deformability.png",
    dpi=300
)

plt.show()
