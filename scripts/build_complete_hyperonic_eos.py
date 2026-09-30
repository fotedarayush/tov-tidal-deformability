#!/usr/bin/env python3
"""
Build complete crust+inner-crust+hyperonic-core EoS tables for the
100 selected zero-temperature hyperonic models.

Expected inputs:
  eos/crust/GPPVA_TW/GPPVA_TW_full.dat
  eos/hyperonic_dataset/hyperonic_100/core_tables/*.dat

Output:
  eos/hyperonic_dataset/hyperonic_100/complete_tables/*.dat
  eos/hyperonic_dataset/hyperonic_100/crust_matching_diagnostics.csv

Physics prescription:
  - BPS outer-crust segment is taken from the GPPVA(TW) CompOSE table
    up to nB = 1e-4 fm^-3.
  - Inner crust: P(epsilon) = a1 + a2 epsilon^(4/3)
  - The bridge joins the outer crust at nB=1e-4 fm^-3 and each
    hyperonic core at nB=0.04 fm^-3.

The TOV/tidal observables M, R, k2 and Lambda depend on P(epsilon).
nB and muB in the bridge are auxiliary columns required by the existing
Fortran input format. nB is reconstructed from the zero-temperature
first-law relation and given a small smooth endpoint correction so that
both joins are exactly continuous.
"""

from pathlib import Path
import csv
import numpy as np

CRUST = Path("eos/crust/GPPVA_TW/GPPVA_TW_full.dat")
CORE_DIR = Path("eos/hyperonic_dataset/hyperonic_100/core_tables")
OUT_DIR = Path("eos/hyperonic_dataset/hyperonic_100/complete_tables")
DIAG = Path("eos/hyperonic_dataset/hyperonic_100/crust_matching_diagnostics.csv")

N_CRUST_JOIN = 1.0e-4
N_CORE_JOIN = 4.0e-2
GAMMA = 4.0 / 3.0
N_BRIDGE = 240

OUT_DIR.mkdir(parents=True, exist_ok=True)


def strictly_increasing(x):
    return bool(np.all(np.diff(np.asarray(x)) > 0))


def interp_logx_logy(x, y, x0):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.any(x <= 0) or np.any(y <= 0) or x0 <= 0:
        raise ValueError("Log interpolation requires positive values.")
    return float(np.exp(np.interp(np.log(x0), np.log(x), np.log(y))))


def cumulative_trapezoid(y, x):
    out = np.zeros_like(x, dtype=float)
    if len(x) > 1:
        out[1:] = np.cumsum(
            0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1])
        )
    return out


crust = np.loadtxt(CRUST)
if crust.ndim != 2 or crust.shape[1] < 4:
    raise RuntimeError(f"{CRUST} must have at least four columns.")

crust = crust[:, :4]
crust = crust[np.argsort(crust[:, 2])]

# Only the low-density part is needed.
if not (crust[0, 2] < N_CRUST_JOIN < crust[-1, 2]):
    raise RuntimeError("nB=1e-4 fm^-3 is outside the crust table.")

# Interpolate an exact endpoint at nB=1e-4.
eps1 = interp_logx_logy(crust[:, 2], crust[:, 0], N_CRUST_JOIN)
P1   = interp_logx_logy(crust[:, 2], crust[:, 1], N_CRUST_JOIN)
mu1_tab = interp_logx_logy(crust[:, 2], crust[:, 3], N_CRUST_JOIN)
mu1_gibbs = (eps1 + P1) / N_CRUST_JOIN

outer = crust[crust[:, 2] < N_CRUST_JOIN].copy()
outer_endpoint = np.array([[eps1, P1, N_CRUST_JOIN, mu1_gibbs]])
outer = np.vstack([outer, outer_endpoint])

# Validate only the BPS outer segment we actually use.
for k, name in enumerate(["epsilon", "P", "nB", "muB"]):
    if not strictly_increasing(outer[:, k]):
        raise RuntimeError(
            f"Used outer-crust segment is not strictly increasing in {name}."
        )

core_files = sorted(CORE_DIR.glob("*.dat"))
if len(core_files) != 100:
    raise RuntimeError(f"Expected 100 core files, found {len(core_files)}.")

diagnostics = []

for number, core_path in enumerate(core_files, 1):
    core = np.loadtxt(core_path)
    if core.ndim != 2 or core.shape[1] < 4:
        raise RuntimeError(f"{core_path} has invalid shape.")

    core = core[:, :4]

    for k, name in enumerate(["epsilon", "P", "nB", "muB"]):
        if not strictly_increasing(core[:, k]):
            raise RuntimeError(f"{core_path.name}: {name} is not increasing.")

    eps2, P2, n2, mu2_tab = core[0]

    if abs(n2 - N_CORE_JOIN) > 1e-10:
        raise RuntimeError(
            f"{core_path.name}: first core nB={n2}, expected {N_CORE_JOIN}."
        )

    # P = a1 + a2 * epsilon^gamma, fixed by pressure continuity
    # at the two energy-density endpoints.
    denom = eps2**GAMMA - eps1**GAMMA
    a2 = (P2 - P1) / denom
    a1 = P1 - a2 * eps1**GAMMA

    if a2 <= 0:
        raise RuntimeError(f"{core_path.name}: non-positive a2.")

    # Full bridge grid INCLUDING endpoints for thermodynamic integration.
    eps_grid = np.geomspace(eps1, eps2, N_BRIDGE + 2)
    P_grid = a1 + a2 * eps_grid**GAMMA

    if np.any(P_grid <= 0) or not strictly_increasing(P_grid):
        raise RuntimeError(f"{core_path.name}: invalid bridge pressure.")

    # Zero-temperature first law:
    # d epsilon = (epsilon + P)/n dn
    # -> d ln n / d epsilon = 1/(epsilon + P).
    integrand = 1.0 / (eps_grid + P_grid)
    I = cumulative_trapezoid(integrand, eps_grid)
    n_raw = N_CRUST_JOIN * np.exp(I)

    # The paper's P(epsilon) prescription fixes pressure endpoints but does
    # not tabulate nB through the phenomenological bridge.  Because this
    # Fortran code requires nB, apply a smooth log correction that is zero
    # at the outer-crust endpoint and enforces nB=0.04 at the core endpoint.
    log_factor = np.log(N_CORE_JOIN / n_raw[-1])
    progress = I / I[-1]
    n_grid = n_raw * np.exp(progress * log_factor)

    # Gibbs relation at T=0.
    mu_grid = (eps_grid + P_grid) / n_grid

    # Use only interior bridge points. Exact endpoints come from outer/core.
    bridge = np.column_stack(
        [eps_grid[1:-1], P_grid[1:-1], n_grid[1:-1], mu_grid[1:-1]]
    )

    combined = np.vstack([outer, bridge, core])

    checks = {}
    for k, name in enumerate(["epsilon", "P", "nB", "muB"]):
        checks[name] = strictly_increasing(combined[:, k])

    if not all(checks.values()):
        raise RuntimeError(
            f"{core_path.name}: complete table is not strictly monotonic: "
            f"{checks}"
        )

    # Sound-speed check on the phenomenological inner-crust bridge.
    cs2 = GAMMA * a2 * eps_grid**(GAMMA - 1.0)
    if np.any(cs2 <= 0) or np.any(cs2 >= 1):
        raise RuntimeError(
            f"{core_path.name}: bridge has nonphysical dP/depsilon range "
            f"{cs2.min()} to {cs2.max()}."
        )

    out_name = core_path.name.replace("_core.dat", "_complete.dat")
    out_path = OUT_DIR / out_name

    with out_path.open("w") as f:
        f.write("# epsilon[MeV/fm^3] P[MeV/fm^3] nB[fm^-3] muB[MeV]\n")
        np.savetxt(f, combined, fmt="%.16e")

    diagnostics.append({
        "number": number,
        "core_file": core_path.name,
        "complete_file": out_name,
        "eps_crust_join": eps1,
        "P_crust_join": P1,
        "mu_crust_table": mu1_tab,
        "mu_crust_gibbs": mu1_gibbs,
        "eps_core_join": eps2,
        "P_core_join": P2,
        "n_core_join": n2,
        "mu_core_table": mu2_tab,
        "mu_core_gibbs": (eps2 + P2) / n2,
        "a1": a1,
        "a2": a2,
        "n_raw_at_core": n_raw[-1],
        "n_endpoint_correction_percent": 100.0 * (N_CORE_JOIN / n_raw[-1] - 1.0),
        "cs2_min_bridge": float(cs2.min()),
        "cs2_max_bridge": float(cs2.max()),
        "rows_complete": len(combined),
    })

    print(f"[{number:3d}/100] {out_name}")

with DIAG.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=diagnostics[0].keys())
    writer.writeheader()
    writer.writerows(diagnostics)

print("\nFinished.")
print("Complete tables:", OUT_DIR)
print("Diagnostics    :", DIAG)
print("Outer join:")
print(f"  epsilon = {eps1:.12e} MeV/fm^3")
print(f"  P       = {P1:.12e} MeV/fm^3")
print(f"  nB      = {N_CRUST_JOIN:.12e} fm^-3")
print(f"  muB     = {mu1_gibbs:.12e} MeV")
