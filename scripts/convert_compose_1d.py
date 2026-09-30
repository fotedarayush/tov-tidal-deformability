from pathlib import Path
import numpy as np
import sys


if len(sys.argv) != 3:
    raise SystemExit(
        "Usage:\n"
        "python3 scripts/convert_compose_1d.py "
        "<CompOSE folder> <output file>"
    )


folder = Path(sys.argv[1])
output = Path(sys.argv[2])

nb_file = folder / "eos.nb"
thermo_file = folder / "eos.thermo"


# ------------------------------------------------------------
# Read baryon-density grid
# ------------------------------------------------------------

with nb_file.open() as f:
    lines = [line.strip() for line in f if line.strip()]

imin = int(lines[0])
imax = int(lines[1])

nb = np.array([float(x) for x in lines[2:]], dtype=float)

expected_n = imax - imin + 1

if len(nb) != expected_n:
    raise ValueError(
        f"Expected {expected_n} density points, "
        f"but found {len(nb)}."
    )


# ------------------------------------------------------------
# Read eos.thermo
# ------------------------------------------------------------

with thermo_file.open() as f:
    thermo_lines = [line.strip() for line in f if line.strip()]

# First row: neutron mass, proton mass, lepton flag
header = thermo_lines[0].split()

mn = float(header[0])
mp = float(header[1])
lepton_flag = int(header[2])

print(f"mn = {mn:.8f} MeV")
print(f"mp = {mp:.8f} MeV")
print(f"leptons included = {bool(lepton_flag)}")


rows = []

for line in thermo_lines[1:]:

    x = line.split()

    if len(x) < 11:
        continue

    iT = int(x[0])
    inb = int(x[1])
    iYq = int(x[2])

    Q1 = float(x[3])
    Q2 = float(x[4])
    Q3 = float(x[5])
    Q4 = float(x[6])
    Q5 = float(x[7])
    Q6 = float(x[8])
    Q7 = float(x[9])
    Nadd = int(x[10])

    # Convert CompOSE density index to Python index.
    j = inb - imin

    if j < 0 or j >= len(nb):
        raise ValueError(f"Invalid density index: {inb}")

    nB = nb[j]

    # Required quantities:
    pressure = nB * Q1
    muB = mn * (1.0 + Q3)
    epsilon = nB * mn * (1.0 + Q7)

    rows.append(
        (inb, epsilon, pressure, nB, muB)
    )


# Sort explicitly by density index
rows.sort(key=lambda x: x[0])


# ------------------------------------------------------------
# Basic consistency checks
# ------------------------------------------------------------

if len(rows) != len(nb):
    raise ValueError(
        f"Found {len(rows)} thermo rows but {len(nb)} density points."
    )

data = np.array(
    [[eps, p, n, mu] for _, eps, p, n, mu in rows]
)


if not np.all(np.diff(data[:, 0]) > 0):
    print("WARNING: epsilon is not strictly increasing.")

if not np.all(np.diff(data[:, 1]) > 0):
    print("WARNING: pressure is not strictly increasing.")

if not np.all(np.diff(data[:, 2]) > 0):
    print("WARNING: nB is not strictly increasing.")


# ------------------------------------------------------------
# Write TOV-format table
# ------------------------------------------------------------

output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w") as f:

    # IMPORTANT:
    # Your Fortran code discards the first line as a header.
    f.write("# epsilon[MeV/fm^3] P[MeV/fm^3] nB[fm^-3] muB[MeV]\n")

    for eps, p, n, mu in data:
        f.write(
            f"{eps:.16e} "
            f"{p:.16e} "
            f"{n:.16e} "
            f"{mu:.16e}\n"
        )


print()
print(f"Density points : {len(data)}")
print(f"nB range       : {data[0,2]:.6e} -> {data[-1,2]:.6e} fm^-3")
print(f"epsilon range  : {data[0,0]:.6e} -> {data[-1,0]:.6e} MeV/fm^3")
print(f"pressure range : {data[0,1]:.6e} -> {data[-1,1]:.6e} MeV/fm^3")
print()
print(f"Wrote: {output}")
