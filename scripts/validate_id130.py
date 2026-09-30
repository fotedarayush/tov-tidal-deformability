#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REF = Path(
    "eos/hyperonic_dataset/"
    "hyperonic_mr_curve_at_zero_temperature.csv.gz"
)
TOV = Path(
    "results/hyperonic_100/test_id130/"
    "hyperon_EOS_001_ID_00130_complete_TOV.dat"
)

EOS_ID = 130


def stable_branch(M, *ys):
    imax = int(np.nanargmax(M))
    out = [M[:imax+1]]
    out.extend(y[:imax+1] for y in ys)
    return imax, out


def interp(M, y, mass):
    order = np.argsort(M)
    M = M[order]
    y = y[order]
    if mass < M[0] or mass > M[-1]:
        return np.nan
    return float(np.interp(mass, M, y))


ref = pd.read_csv(REF)
ref = ref[ref["ID"] == EOS_ID].copy()
if ref.empty:
    raise RuntimeError(f"ID {EOS_ID} absent from reference M-R file.")

Mr = ref["M"].to_numpy(float)
Rr = ref["R"].to_numpy(float)
ir, (Mr, Rr) = stable_branch(Mr, Rr)

tov = np.loadtxt(TOV)
Mt = tov[:, 2]
Rt = tov[:, 1]  # current source already writes km
it, (Mt, Rt) = stable_branch(Mt, Rt)

values = {
    "Mmax": (Mr.max(), Mt.max()),
    "R_at_Mmax": (Rr[np.argmax(Mr)], Rt[np.argmax(Mt)]),
    "R1.4": (interp(Mr, Rr, 1.4), interp(Mt, Rt, 1.4)),
    "R1.6": (interp(Mr, Rr, 1.6), interp(Mt, Rt, 1.6)),
}

print("ID 130: reference vs local TOV")
print("=" * 72)
print(f"{'Quantity':<14} {'Reference':>14} {'Local':>14} {'Difference':>14}")
for name, (a, b) in values.items():
    print(f"{name:<14} {a:14.7f} {b:14.7f} {b-a:14.7f}")

mhi = min(Mr.max(), Mt.max(), 2.0)
mlo = max(Mr.min(), Mt.min(), 1.0)
grid = np.linspace(mlo, mhi, 100)
rr = np.interp(grid, np.sort(Mr), Rr[np.argsort(Mr)])
rt = np.interp(grid, np.sort(Mt), Rt[np.argsort(Mt)])
rms = np.sqrt(np.mean((rt - rr)**2))

print()
print(f"Stable reference points: {len(Mr)}")
print(f"Stable local points    : {len(Mt)}")
print(f"Radius RMS difference over {mlo:.2f}-{mhi:.2f} Msun: {rms:.5f} km")
