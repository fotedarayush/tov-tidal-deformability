#!/usr/bin/env python3

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Configuration
# ============================================================

INPUT_DIR = Path("results/tov")
OUTPUT_DIR = Path("results/analysis")

# Mass cut used to remove the very-low-mass region from plots
MIN_PLOT_MASS = 0.5

# Reference mass for canonical neutron-star quantities
REFERENCE_MASS = 1.4

# Central-density conversion used by the Fortran code
RHO_CONVERSION = 6.176e17

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Helper functions
# ============================================================

def eos_label(path: Path) -> str:
    """
    Convert an output filename such as
    EoS_dd2_npY_T0_beta_eq_01_05_Mg_800_TOV.dat
    into a shorter plot label.
    """
    name = path.stem

    if name.endswith("_TOV"):
        name = name[:-4]

    prefix = "EoS_dd2_npY_T0_beta_eq_"
    if name.startswith(prefix):
        name = name[len(prefix):]

    return name


def load_eos(path: Path):
    """
    Load one TOV/tidal output file.

    Current Fortran output convention:
      col 1  : central total energy density in internal units
      col 2  : radius in km
      col 3  : gravitational mass in solar masses
      col 10 : k2
      col 12 : dimensionless tidal deformability Lambda
    """

    data = np.loadtxt(path)

    # Protect against a file containing only one row
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 12:
        raise ValueError(
            f"{path} has only {data.shape[1]} columns; at least 12 are required."
        )

    rho_c = data[:, 0] * RHO_CONVERSION
    radius = data[:, 1]
    mass = data[:, 2]
    k2 = data[:, 9]
    tidal_lambda = data[:, 11]

    return rho_c, radius, mass, k2, tidal_lambda


def stable_branch(rho_c, radius, mass, k2, tidal_lambda):
    """
    For a standard one-family TOV sequence, retain points from the
    beginning of the sequence up to the global maximum mass.

    This assumes the normal stable branch is continuous up to Mmax.
    More complicated third-family/twin-star sequences require a more
    detailed stability analysis.
    """

    imax = np.argmax(mass)

    return (
        rho_c[:imax + 1],
        radius[:imax + 1],
        mass[:imax + 1],
        k2[:imax + 1],
        tidal_lambda[:imax + 1],
        imax,
    )


def interpolate_at_mass(mass, quantity, target_mass):
    """
    Interpolate a quantity at a target gravitational mass.

    The data are sorted by mass first because np.interp requires an
    increasing x-array.
    """

    finite = np.isfinite(mass) & np.isfinite(quantity)
    m = mass[finite]
    q = quantity[finite]

    if len(m) < 2:
        return np.nan

    order = np.argsort(m)
    m = m[order]
    q = q[order]

    # Remove repeated mass values to avoid ambiguous interpolation
    m_unique, unique_idx = np.unique(m, return_index=True)
    q_unique = q[unique_idx]

    if not (m_unique.min() <= target_mass <= m_unique.max()):
        return np.nan

    return np.interp(target_mass, m_unique, q_unique)


# ============================================================
# Locate files
# ============================================================

files = sorted(INPUT_DIR.glob("*_TOV.dat"))

if not files:
    raise SystemExit(
        f"No *_TOV.dat files found in {INPUT_DIR.resolve()}\n"
        "Run the batch Fortran script first or change INPUT_DIR."
    )

print(f"Found {len(files)} EoS result files.\n")

# ============================================================
# Analyse every EoS
# ============================================================

summary_rows = []
plot_data = []

for path in files:

    label = eos_label(path)

    try:
        rho_c, R, M, k2, Lambda = load_eos(path)

        (
            rho_stable,
            R_stable,
            M_stable,
            k2_stable,
            Lambda_stable,
            imax,
        ) = stable_branch(rho_c, R, M, k2, Lambda)

        # Maximum-mass properties
        Mmax = M[imax]
        R_Mmax = R[imax]
        rho_Mmax = rho_c[imax]

        # Canonical 1.4 Msun radius
        R14 = interpolate_at_mass(
            M_stable,
            R_stable,
            REFERENCE_MASS
        )

        # Lambda should be positive for the physical tidal sequence.
        tidal_mask = (
            np.isfinite(M_stable)
            & np.isfinite(Lambda_stable)
            & (Lambda_stable > 0)
        )

        Lambda14 = interpolate_at_mass(
            M_stable[tidal_mask],
            Lambda_stable[tidal_mask],
            REFERENCE_MASS
        )

        k214 = interpolate_at_mass(
            M_stable[tidal_mask],
            k2_stable[tidal_mask],
            REFERENCE_MASS
        )

        summary_rows.append({
            "eos": label,
            "file": path.name,
            "Mmax_Msun": Mmax,
            "R_at_Mmax_km": R_Mmax,
            "rho_c_at_Mmax_g_cm3": rho_Mmax,
            "R_1.4_km": R14,
            "Lambda_1.4": Lambda14,
            "k2_1.4": k214,
        })

        plot_data.append({
            "label": label,
            "rho": rho_stable,
            "R": R_stable,
            "M": M_stable,
            "Lambda": Lambda_stable,
        })

        print(
            f"{label:25s} "
            f"Mmax={Mmax:6.3f} Msun   "
            f"R1.4={R14:7.3f} km   "
            f"Lambda1.4={Lambda14:9.3f}"
        )

    except Exception as exc:
        print(f"FAILED: {path.name}: {exc}")

# ============================================================
# Save summary CSV
# ============================================================

summary_file = OUTPUT_DIR / "eos_summary.csv"

fieldnames = [
    "eos",
    "file",
    "Mmax_Msun",
    "R_at_Mmax_km",
    "rho_c_at_Mmax_g_cm3",
    "R_1.4_km",
    "Lambda_1.4",
    "k2_1.4",
]

with summary_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"\nSaved summary table to: {summary_file}")

# ============================================================
# Combined Mass-Radius plot
# ============================================================

plt.figure(figsize=(8, 6))

for item in plot_data:

    mask = (
        np.isfinite(item["M"])
        & np.isfinite(item["R"])
        & (item["M"] >= MIN_PLOT_MASS)
    )

    plt.plot(
        item["R"][mask],
        item["M"][mask],
        label=item["label"]
    )

plt.xlabel(r"$R$ [km]")
plt.ylabel(r"$M/M_\odot$")
plt.title("Mass-Radius Relations")
plt.grid(True)
plt.legend(fontsize=7, ncol=2)
plt.tight_layout()

mr_file = OUTPUT_DIR / "all_eos_mass_radius.png"
plt.savefig(mr_file, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {mr_file}")

# ============================================================
# Combined tidal deformability plot
# ============================================================

plt.figure(figsize=(8, 6))

for item in plot_data:

    mask = (
        np.isfinite(item["M"])
        & np.isfinite(item["Lambda"])
        & (item["M"] >= MIN_PLOT_MASS)
        & (item["Lambda"] > 0)
    )

    # Sort by mass for a clean connected curve
    Mplot = item["M"][mask]
    Lplot = item["Lambda"][mask]

    order = np.argsort(Mplot)

    plt.plot(
        Mplot[order],
        Lplot[order],
        label=item["label"]
    )

plt.xlabel(r"$M/M_\odot$")
plt.ylabel(r"$\Lambda$")
plt.title("Tidal Deformability")
plt.yscale("log")
plt.grid(True)
plt.legend(fontsize=7, ncol=2)
plt.tight_layout()

tidal_file = OUTPUT_DIR / "all_eos_tidal_deformability.png"
plt.savefig(tidal_file, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {tidal_file}")

# ============================================================
# Mmax comparison plot
# ============================================================

if summary_rows:

    labels = [row["eos"] for row in summary_rows]
    mmax_values = [row["Mmax_Msun"] for row in summary_rows]

    plt.figure(figsize=(10, 6))
    x = np.arange(len(labels))

    plt.scatter(x, mmax_values)
    plt.axhline(2.0, linestyle="--")

    plt.xticks(x, labels, rotation=90)
    plt.ylabel(r"$M_{\max}/M_\odot$")
    plt.title("Maximum Mass by EoS")
    plt.grid(True, axis="y")
    plt.tight_layout()

    mmax_file = OUTPUT_DIR / "all_eos_maximum_mass.png"
    plt.savefig(mmax_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {mmax_file}")

# ============================================================
# R1.4 comparison plot
# ============================================================

valid_r14 = [
    row for row in summary_rows
    if np.isfinite(row["R_1.4_km"])
]

if valid_r14:

    labels = [row["eos"] for row in valid_r14]
    values = [row["R_1.4_km"] for row in valid_r14]

    plt.figure(figsize=(10, 6))
    x = np.arange(len(labels))

    plt.scatter(x, values)

    plt.xticks(x, labels, rotation=90)
    plt.ylabel(r"$R_{1.4}$ [km]")
    plt.title(r"$R_{1.4}$ by EoS")
    plt.grid(True, axis="y")
    plt.tight_layout()

    r14_file = OUTPUT_DIR / "all_eos_R1.4.png"
    plt.savefig(r14_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {r14_file}")

# ============================================================
# Lambda1.4 comparison plot
# ============================================================

valid_lambda14 = [
    row for row in summary_rows
    if np.isfinite(row["Lambda_1.4"])
]

if valid_lambda14:

    labels = [row["eos"] for row in valid_lambda14]
    values = [row["Lambda_1.4"] for row in valid_lambda14]

    plt.figure(figsize=(10, 6))
    x = np.arange(len(labels))

    plt.scatter(x, values)

    plt.xticks(x, labels, rotation=90)
    plt.ylabel(r"$\Lambda_{1.4}$")
    plt.title(r"$\Lambda_{1.4}$ by EoS")
    plt.grid(True, axis="y")
    plt.tight_layout()

    lambda_file = OUTPUT_DIR / "all_eos_Lambda1.4.png"
    plt.savefig(lambda_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {lambda_file}")

print("\nAnalysis complete.")
