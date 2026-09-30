#!/usr/bin/env python3
from pathlib import Path
import csv
import numpy as np

IN_DIR = Path("results/hyperonic_100/coarse")
OUT = Path("results/hyperonic_100/coarse_summary.csv")


def interp(M, y, mass):
    mask = np.isfinite(M) & np.isfinite(y)
    M = M[mask]
    y = y[mask]
    order = np.argsort(M)
    M, y = M[order], y[order]
    if len(M) < 2 or mass < M[0] or mass > M[-1]:
        return np.nan
    return float(np.interp(mass, M, y))


rows = []
files = sorted(IN_DIR.glob("*_TOV.dat"))

for path in files:
    d = np.loadtxt(path)
    if d.ndim == 1:
        d = d.reshape(1, -1)

    M = d[:, 2]
    R = d[:, 1]
    Lam = d[:, 11]

    i = int(np.nanargmax(M))
    resolved = i < len(M) - 2

    Ms = M[:i+1]
    Rs = R[:i+1]
    Ls = Lam[:i+1]

    rows.append({
        "eos": path.stem.replace("_TOV", ""),
        "rows": len(d),
        "resolved_mmax": resolved,
        "index_mmax": i,
        "models_after_mmax": len(M) - 1 - i,
        "Mmax": float(M[i]),
        "R_Mmax_km": float(R[i]),
        "last_mass": float(M[-1]),
        "R1.4_km": interp(Ms, Rs, 1.4),
        "R1.6_km": interp(Ms, Rs, 1.6),
        "Lambda1.4": interp(Ms, Ls, 1.4),
        "Lambda1.6": interp(Ms, Ls, 1.6),
    })

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

resolved = sum(r["resolved_mmax"] for r in rows)

print("Run files found :", len(rows))
print("Resolved Mmax   :", resolved)
print("Unresolved Mmax :", len(rows) - resolved)
print("Summary         :", OUT)

if len(rows) != 100:
    print("WARNING: expected 100 output files.")

if len(rows) - resolved:
    print("\nUnresolved:")
    for r in rows:
        if not r["resolved_mmax"]:
            print(" ", r["eos"])
