from pathlib import Path
import pandas as pd
import numpy as np
import zipfile

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

SOURCE = Path(
    "eos/hyperonic_dataset/"
    "hyperonic_eos_at_zero_temperature.csv.gz"
)

OUTDIR = Path("eos/hyperonic_dataset/hyperonic_100")
TABLEDIR = OUTDIR / "core_tables"

N_EOS = 100
RANDOM_SEED = 20260930

OUTDIR.mkdir(parents=True, exist_ok=True)
TABLEDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Load dataset
# ------------------------------------------------------------

print("Loading dataset...")

df = pd.read_csv(SOURCE)

print("Rows:", len(df))
print("Distinct EoS IDs:", df["ID"].nunique())
print("Columns:", list(df.columns))

required = {"ID", "n", "rho", "mu", "e", "p"}

if not required.issubset(df.columns):
    raise ValueError(
        f"Expected columns {required}, "
        f"but found {set(df.columns)}"
    )

# ------------------------------------------------------------
# Find numerically clean EoSs
# ------------------------------------------------------------

print("\nChecking monotonicity...")

valid_ids = []

for eos_id, g in df.groupby("ID", sort=True):

    g = g.sort_values("n")

    rho_ok = np.all(np.diff(g["rho"].to_numpy()) > 0)
    e_ok   = np.all(np.diff(g["e"].to_numpy()) > 0)
    p_ok   = np.all(np.diff(g["p"].to_numpy()) > 0)
    mu_ok  = np.all(np.diff(g["mu"].to_numpy()) > 0)

    if rho_ok and e_ok and p_ok and mu_ok:
        valid_ids.append(int(eos_id))

print("Valid monotonic EoSs:", len(valid_ids))

if len(valid_ids) < N_EOS:
    raise RuntimeError(
        f"Only {len(valid_ids)} suitable EoSs found."
    )

# ------------------------------------------------------------
# Reproducibly choose 100
# ------------------------------------------------------------

rng = np.random.default_rng(RANDOM_SEED)

selected_ids = np.sort(
    rng.choice(valid_ids, size=N_EOS, replace=False)
)

print("\nSelected IDs:")
print(selected_ids)

# Save exact IDs
with open(OUTDIR / "selected_ids.txt", "w") as f:
    for eos_id in selected_ids:
        f.write(f"{eos_id}\n")

# ------------------------------------------------------------
# Write individual four-column tables
# ------------------------------------------------------------

manifest = []

for number, eos_id in enumerate(selected_ids, start=1):

    g = (
        df[df["ID"] == eos_id]
        .sort_values("n")
        .copy()
    )

    filename = (
        f"hyperon_EOS_{number:03d}"
        f"_ID_{eos_id:05d}_core.dat"
    )

    outfile = TABLEDIR / filename

    with outfile.open("w") as f:

        # Your Fortran code skips the first line,
        # so keep this header.
        f.write(
            "# epsilon[MeV/fm^3] "
            "P[MeV/fm^3] "
            "nB[fm^-3] "
            "muB[MeV]\n"
        )

        for row in g.itertuples(index=False):

            # Source mapping:
            # e   -> epsilon
            # p   -> pressure
            # rho -> baryon density
            # mu  -> baryon chemical potential

            f.write(
                f"{row.e:.16e} "
                f"{row.p:.16e} "
                f"{row.rho:.16e} "
                f"{row.mu:.16e}\n"
            )

    manifest.append(
        {
            "number": number,
            "ID": eos_id,
            "filename": filename,
            "rows": len(g),
            "nB_min": g["rho"].iloc[0],
            "nB_max": g["rho"].iloc[-1],
            "epsilon_min": g["e"].iloc[0],
            "epsilon_max": g["e"].iloc[-1],
            "P_min": g["p"].iloc[0],
            "P_max": g["p"].iloc[-1],
        }
    )

    print(
        f"[{number:3d}/100] "
        f"ID {eos_id} -> {filename}"
    )

# ------------------------------------------------------------
# Manifest
# ------------------------------------------------------------

manifest_df = pd.DataFrame(manifest)

manifest_df.to_csv(
    OUTDIR / "manifest.csv",
    index=False
)

# Preserve selected EoSs in original CSV format too
native = (
    df[df["ID"].isin(selected_ids)]
    .sort_values(["ID", "n"])
)

native.to_csv(
    OUTDIR / "selected_100_native.csv.gz",
    index=False,
    compression="gzip",
)

# ------------------------------------------------------------
# README
# ------------------------------------------------------------

readme = """
100 zero-temperature hyperonic EoSs

Source:
hyperonic_eos_at_zero_temperature.csv.gz

Selection:
100 distinct EoSs chosen reproducibly using
NumPy random seed 20260930.

Four-column table format:

epsilon [MeV/fm^3]
P       [MeV/fm^3]
nB      [fm^-3]
muB     [MeV]

Source-column mapping:

epsilon = e
P       = p
nB      = rho
muB     = mu

IMPORTANT:

These tables contain the hyperonic CORE EoS from the
source dataset. They do not necessarily extend to the
very low densities required to reach a neutron-star
surface.

A suitable crust EoS should therefore be attached
before using these as complete production TOV/tidal
tables.

No quark phase is added by this conversion.
"""

(OUTDIR / "README.txt").write_text(readme.strip() + "\n")

# ------------------------------------------------------------
# ZIP
# ------------------------------------------------------------

zip_path = Path(
    "eos/hyperonic_dataset/hyperonic_100_core_eos.zip"
)

with zipfile.ZipFile(
    zip_path,
    "w",
    compression=zipfile.ZIP_DEFLATED,
) as zf:

    for file in OUTDIR.rglob("*"):

        if file.is_file():
            zf.write(
                file,
                file.relative_to(OUTDIR)
            )

print("\n====================================")
print("Finished")
print("====================================")
print("Selected EoSs:", len(selected_ids))
print("Tables:", TABLEDIR)
print("Manifest:", OUTDIR / "manifest.csv")
print("ZIP:", zip_path)
