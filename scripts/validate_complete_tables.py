#!/usr/bin/env python3
from pathlib import Path
import numpy as np

folder = Path("eos/hyperonic_dataset/hyperonic_100/complete_tables")
files = sorted(folder.glob("*.dat"))

print("Number of complete tables:", len(files))

bad = []
for p in files:
    d = np.loadtxt(p)
    if d.ndim != 2 or d.shape[1] < 4:
        bad.append((p.name, "shape"))
        continue

    names = ["epsilon", "P", "nB", "muB"]
    failed = [
        names[i]
        for i in range(4)
        if not np.all(np.diff(d[:, i]) > 0)
    ]

    if np.any(d[:, :4] <= 0):
        failed.append("non-positive value")

    if failed:
        bad.append((p.name, ", ".join(failed)))

print("Bad tables:", len(bad))
for x in bad:
    print("  ", x)

if len(files) != 100 or bad:
    raise SystemExit(1)

print("All 100 complete tables passed positivity and monotonicity checks.")
