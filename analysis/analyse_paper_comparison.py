#!/usr/bin/env python3
"""
Paper-oriented multi-EoS analysis for TOV + tidal-deformability outputs.

Designed to make the output directly comparable in structure to:

A. Bauswein et al.,
"Stellar properties indicating the presence of hyperons in neutron stars",
Phys. Rev. Research 8, 013253 (2026),
DOI: 10.1103/ygtr-ktqk

The script analyses all *_TOV.dat files in results/tov/ and produces:

  1. Combined M-R curves
  2. Combined dimensionless tidal deformability Lambda(M)
  3. Combined dimensional/geometric tidal deformability lambda(M)
  4. Signed curvature kappa_R(M), as in the paper
  5. kappa_R(M_ref) versus R(M_ref), paper Fig. 3 style
  6. d^2 lambda / dM^2 versus M, paper Fig. 4 style
  7. d^2 lambda / dM^2(M_ref) versus lambda(M_ref), paper Fig. 5 style
  8. dR/dM(M_ref) versus Mmax, paper Fig. 6 style
  9. dLambda/dM(M_ref) versus Mmax, paper Fig. 6 style
 10. CSV summary of all derived stellar quantities
 11. A warnings file flagging unresolved Mmax and missing reference masses

IMPORTANT:
- The present files appear to be hyperonic npY EoSs. The paper compares
  hyperonic and purely nucleonic samples. To reproduce that comparison
  directly, add nucleonic control EoSs to the input directory.
- The paper uses centered, second-order finite differences with variable
  step size. This script implements that formula.
- Second derivatives are sensitive to numerical resolution. For final
  publication-quality derivative work, use a fine central-density grid.
"""

from pathlib import Path
import csv
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# SciencePlots is required by request.
try:
    import scienceplots  # noqa: F401
except ImportError:
    raise SystemExit(
        "\nSciencePlots is not installed.\n"
        "Install it with:\n\n"
        "    python3 -m pip install SciencePlots\n\n"
        "Then run this script again.\n"
    )

# ---------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------
# "no-latex" avoids requiring a full LaTeX installation.
# If you have LaTeX installed and want the full SciencePlots look,
# change this to: plt.style.use(["science", "grid"])
plt.style.use(["science", "no-latex", "grid"])

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

INPUT_DIR = Path("results/tov")
OUTPUT_DIR = Path("results/paper_comparison")

# The paper commonly evaluates fixed-mass quantities at 1.6 Msun.
PAPER_REFERENCE_MASS = 1.6

# Additional useful reference masses for your own comparison.
REFERENCE_MASSES = [1.4, 1.6, 1.8]

# Do not use the pathological very-low-mass branch for derivative work.
MIN_ANALYSIS_MASS = 0.5

# Paper discusses curvature mainly above about 1 Msun.
MIN_CURVATURE_MASS = 1.0

# The paper terminates derivative curves slightly before Mmax for
# numerical reasons. This buffer does the same.
MMAX_DERIVATIVE_BUFFER = 0.03  # Msun

# If adjacent mass points are essentially identical, merge them before
# differentiating. This is particularly useful around constant-pressure
# plateaus or repeated numerical solutions.
MASS_DUPLICATE_TOL = 2.0e-6  # Msun

# Conversion used by your Fortran output.
RHO_CONVERSION = 6.176e17  # internal density -> g cm^-3

# GM_sun/c^2 in km. Used because the paper evaluates R(M) curvature
# in geometrized mass units (G=c=1), making kappa_R dimensionless.
GM_SUN_C2_KM = 1.4766250385

# Paper convention: dashed curves if Mmax < 2 Msun.
MASS_CONSTRAINT = 2.0

# Save both raster and vector figures.
SAVE_PNG = True
SAVE_PDF = True
DPI = 400

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =====================================================================
# Filename parsing and labels
# =====================================================================

def parse_filename(path):
    """
    Parse filenames of the form
    EoS_dd2_npY_T0_beta_eq_01_05_Mg_800_TOV.dat

    Interprets 01 -> 0.1 and 05 -> 0.5 for display. If your parameter
    convention differs, change only the display conversion below.
    """
    pattern = (
        r"EoS_dd2_npY_T0_beta_eq_"
        r"(?P<etaD>\d{2})_"
        r"(?P<etaV>\d{2})_"
        r"Mg_(?P<Mg>\d+)_TOV\.dat$"
    )
    match = re.match(pattern, path.name)

    if not match:
        return {
            "etaD_code": "?",
            "etaV_code": "?",
            "etaD": np.nan,
            "etaV": np.nan,
            "Mg": np.nan,
            "label": path.stem,
            "family": path.stem,
        }

    etaD_code = match.group("etaD")
    etaV_code = match.group("etaV")
    Mg = int(match.group("Mg"))

    etaD = int(etaD_code) / 10.0
    etaV = int(etaV_code) / 10.0

    label = (
        rf"$\eta_D={etaD:.1f},\ \eta_V={etaV:.1f},\ "
        rf"M_g={Mg}$"
    )

    # A family shares eta_D and eta_V; Mg is distinguished by marker.
    family = (etaD_code, etaV_code)

    return {
        "etaD_code": etaD_code,
        "etaV_code": etaV_code,
        "etaD": etaD,
        "etaV": etaV,
        "Mg": Mg,
        "label": label,
        "family": family,
    }


# =====================================================================
# Read / clean one EoS
# =====================================================================

def load_eos(path):
    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 12:
        raise ValueError(
            f"{path.name}: expected >=12 columns, found {data.shape[1]}"
        )

    return {
        "raw": data,
        "rho_c": data[:, 0] * RHO_CONVERSION,
        "R_km": data[:, 1],
        "M": data[:, 2],
        "Mrest": data[:, 3],
        "compactness": data[:, 6],
        "k2": data[:, 9],
        "kappa2": data[:, 10],
        "Lambda": data[:, 11],
    }


def select_stable_prefix(eos):
    """
    Keep the sequence from the start through the global maximum mass.

    This is appropriate for a standard one-family TOV sequence.
    A genuine third-family/twin-star branch requires a dedicated radial
    stability analysis and should not be handled by this simple cut.
    """
    M = eos["M"]
    imax = int(np.nanargmax(M))

    # If the maximum occurs at the edge of the input sequence, then the
    # true maximum mass has probably not been reached.
    mmax_resolved = imax < len(M) - 2

    sl = slice(0, imax + 1)

    stable = {
        key: value[sl].copy()
        for key, value in eos.items()
        if key != "raw"
    }

    stable["imax_raw"] = imax
    stable["mmax_resolved"] = mmax_resolved

    return stable


def remove_low_mass(stable):
    mask = (
        np.isfinite(stable["M"])
        & np.isfinite(stable["R_km"])
        & (stable["M"] >= MIN_ANALYSIS_MASS)
        & (stable["R_km"] > 0)
    )

    return {
        key: value[mask].copy()
        for key, value in stable.items()
        if isinstance(value, np.ndarray)
    }


def merge_near_duplicate_masses(branch):
    """
    Sort by mass and merge points whose masses differ by less than
    MASS_DUPLICATE_TOL. This prevents singular finite differences in
    regions where the TOV output repeats essentially the same solution.
    """
    M = branch["M"]
    order = np.argsort(M)

    arrays = {
        key: value[order]
        for key, value in branch.items()
        if isinstance(value, np.ndarray)
    }

    if len(arrays["M"]) == 0:
        return arrays

    groups = [[0]]

    for i in range(1, len(arrays["M"])):
        if abs(arrays["M"][i] - arrays["M"][groups[-1][-1]]) <= MASS_DUPLICATE_TOL:
            groups[-1].append(i)
        else:
            groups.append([i])

    merged = {}

    for key, values in arrays.items():
        merged[key] = np.array(
            [np.nanmean(values[group]) for group in groups],
            dtype=float,
        )

    return merged


# =====================================================================
# Numerical derivatives
# =====================================================================

def centered_fd_variable_step(x, y):
    """
    Second-order centered finite differences for nonuniform x spacing.

    For h_- = x_i - x_{i-1}, h_+ = x_{i+1} - x_i:

    f'(x_i) =
      -h_+ / [h_- (h_- + h_+)] f_{i-1}
      + (h_+ - h_-) / (h_- h_+) f_i
      + h_- / [h_+ (h_- + h_+)] f_{i+1}

    f''(x_i) =
      2 [ f_{i-1}/(h_-(h_-+h_+))
          - f_i/(h_-h_+)
          + f_{i+1}/(h_+(h_-+h_+)) ]

    Endpoints are left as NaN.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    d1 = np.full_like(y, np.nan)
    d2 = np.full_like(y, np.nan)

    if len(x) < 3:
        return d1, d2

    for i in range(1, len(x) - 1):
        hm = x[i] - x[i - 1]
        hp = x[i + 1] - x[i]

        if hm <= 0 or hp <= 0:
            continue

        d1[i] = (
            -hp / (hm * (hm + hp)) * y[i - 1]
            + (hp - hm) / (hm * hp) * y[i]
            + hm / (hp * (hm + hp)) * y[i + 1]
        )

        d2[i] = 2.0 * (
            y[i - 1] / (hm * (hm + hp))
            - y[i] / (hm * hp)
            + y[i + 1] / (hp * (hm + hp))
        )

    return d1, d2


def derive_paper_quantities(branch):
    M = branch["M"]
    R_km = branch["R_km"]
    Lambda = branch["Lambda"]

    # Radius in geometrized solar-mass units, as used for kappa_R.
    R_geom = R_km / GM_SUN_C2_KM

    # lambda = Lambda * M^5 in the paper's geometric mass units.
    # Lambda is dimensionless; lambda has units of M_sun^5 when M is
    # expressed in solar-mass units.
    lambda_geom = Lambda * M**5

    dR_dM, d2R_dM2 = centered_fd_variable_step(M, R_geom)
    dlambda_dM, d2lambda_dM2 = centered_fd_variable_step(M, lambda_geom)
    dLambda_dM, d2Lambda_dM2 = centered_fd_variable_step(M, Lambda)

    # Signed curvature used by Bauswein et al.
    # DO NOT take the absolute value: negative curvature is the signal.
    kappa_R = d2R_dM2 / (1.0 + dR_dM**2) ** 1.5

    out = dict(branch)
    out.update({
        "R_geom": R_geom,
        "lambda_geom": lambda_geom,
        "dR_dM": dR_dM,
        "d2R_dM2": d2R_dM2,
        "kappa_R": kappa_R,
        "dlambda_dM": dlambda_dM,
        "d2lambda_dM2": d2lambda_dM2,
        "dLambda_dM": dLambda_dM,
        "d2Lambda_dM2": d2Lambda_dM2,
    })

    return out


def interp_at_mass(M, y, mass):
    mask = np.isfinite(M) & np.isfinite(y)
    x = M[mask]
    z = y[mask]

    if len(x) < 2:
        return np.nan

    order = np.argsort(x)
    x = x[order]
    z = z[order]

    if mass < x[0] or mass > x[-1]:
        return np.nan

    return float(np.interp(mass, x, z))


# =====================================================================
# Plot helpers
# =====================================================================

def save_figure(fig, stem):
    if SAVE_PNG:
        fig.savefig(
            OUTPUT_DIR / f"{stem}.png",
            dpi=DPI,
            bbox_inches="tight",
        )
    if SAVE_PDF:
        fig.savefig(
            OUTPUT_DIR / f"{stem}.pdf",
            bbox_inches="tight",
        )
    plt.close(fig)


def make_style_maps(models):
    families = sorted({m["meta"]["family"] for m in models}, key=str)
    mgs = sorted({
        m["meta"]["Mg"]
        for m in models
        if np.isfinite(m["meta"]["Mg"])
    })

    cmap = plt.get_cmap("tab10")
    family_colors = {
        family: cmap(i % 10)
        for i, family in enumerate(families)
    }

    marker_cycle = ["o", "s", "^", "D", "v", "P", "X"]
    mg_markers = {
        mg: marker_cycle[i % len(marker_cycle)]
        for i, mg in enumerate(mgs)
    }

    return family_colors, mg_markers


def curve_linestyle(model):
    if not model["mmax_resolved"]:
        return ":"
    if model["Mmax"] < MASS_CONSTRAINT:
        return "--"
    return "-"


def add_legends(ax, models, family_colors, mg_markers):
    # EoS-specific legend
    handles = []

    for model in models:
        meta = model["meta"]
        handles.append(
            Line2D(
                [0], [0],
                color=family_colors[meta["family"]],
                linestyle=curve_linestyle(model),
                marker=mg_markers.get(meta["Mg"], "o"),
                markersize=4,
                linewidth=1.3,
                label=meta["label"],
            )
        )

    legend1 = ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        fontsize=6.5,
        frameon=True,
        title="EoS",
        title_fontsize=7,
    )
    ax.add_artist(legend1)

    # Paper-convention/status legend
    status_handles = [
        Line2D([0], [0], color="black", lw=1.4, ls="-",
               label=rf"$M_{{\max}}\geq {MASS_CONSTRAINT:.1f}\,M_\odot$"),
        Line2D([0], [0], color="black", lw=1.4, ls="--",
               label=rf"$M_{{\max}}< {MASS_CONSTRAINT:.1f}\,M_\odot$"),
        Line2D([0], [0], color="black", lw=1.4, ls=":",
               label=r"$M_{\max}$ not reached"),
    ]

    ax.legend(
        handles=status_handles,
        loc="lower left",
        bbox_to_anchor=(1.02, 0.0),
        borderaxespad=0.0,
        fontsize=6.5,
        frameon=True,
        title="Line style",
        title_fontsize=7,
    )


def derivative_plot_mask(model):
    M = model["M"]
    return (
        np.isfinite(M)
        & (M >= MIN_CURVATURE_MASS)
        & (M <= model["Mmax"] - MMAX_DERIVATIVE_BUFFER)
    )


# =====================================================================
# Read all models
# =====================================================================

files = sorted(INPUT_DIR.glob("*_TOV.dat"))

if not files:
    raise SystemExit(
        f"\nNo *_TOV.dat files found in:\n  {INPUT_DIR.resolve()}\n"
        "Change INPUT_DIR or run the TOV batch calculation first.\n"
    )

models = []
warnings = []

for path in files:
    try:
        meta = parse_filename(path)
        eos = load_eos(path)

        stable = select_stable_prefix(eos)
        imax = stable["imax_raw"]
        mmax_resolved = stable["mmax_resolved"]

        Mmax = float(eos["M"][imax])
        R_Mmax = float(eos["R_km"][imax])
        rho_Mmax = float(eos["rho_c"][imax])

        lowmass_cut = remove_low_mass(stable)
        clean = merge_near_duplicate_masses(lowmass_cut)
        derived = derive_paper_quantities(clean)

        model = {
            **derived,
            "path": path,
            "meta": meta,
            "Mmax": Mmax,
            "R_Mmax": R_Mmax,
            "rho_Mmax": rho_Mmax,
            "mmax_resolved": mmax_resolved,
        }

        models.append(model)

        if not mmax_resolved:
            warnings.append(
                f"{path.name}: maximum mass occurs at the end of the "
                "computed sequence. Reported Mmax is a LOWER BOUND; "
                "extend the central-density range."
            )

        if Mmax < PAPER_REFERENCE_MASS:
            warnings.append(
                f"{path.name}: Mmax={Mmax:.3f} Msun < "
                f"{PAPER_REFERENCE_MASS:.1f} Msun, so paper-style "
                "quantities at the reference mass are unavailable."
            )

    except Exception as exc:
        warnings.append(f"{path.name}: analysis FAILED: {exc}")

if not models:
    raise SystemExit("No EoS files could be analysed.")

family_colors, mg_markers = make_style_maps(models)


# =====================================================================
# Summary table
# =====================================================================

summary_rows = []

for model in models:
    row = {
        "eos": model["meta"]["label"].replace("$", ""),
        "filename": model["path"].name,
        "eta_D": model["meta"]["etaD"],
        "eta_V": model["meta"]["etaV"],
        "Mg": model["meta"]["Mg"],
        "Mmax_Msun": model["Mmax"],
        "Mmax_resolved": model["mmax_resolved"],
        "R_at_Mmax_km": model["R_Mmax"],
        "rho_c_at_Mmax_g_cm3": model["rho_Mmax"],
    }

    for mref in REFERENCE_MASSES:
        tag = f"{mref:.1f}".replace(".", "p")

        row[f"R_{tag}_km"] = interp_at_mass(
            model["M"], model["R_km"], mref
        )
        row[f"Lambda_{tag}"] = interp_at_mass(
            model["M"], model["Lambda"], mref
        )
        row[f"lambda_{tag}_Msun5"] = interp_at_mass(
            model["M"], model["lambda_geom"], mref
        )
        row[f"kappa_R_{tag}"] = interp_at_mass(
            model["M"], model["kappa_R"], mref
        )
        row[f"dR_dM_{tag}"] = interp_at_mass(
            model["M"], model["dR_dM"], mref
        )
        row[f"dLambda_dM_{tag}"] = interp_at_mass(
            model["M"], model["dLambda_dM"], mref
        )
        row[f"d2lambda_dM2_{tag}"] = interp_at_mass(
            model["M"], model["d2lambda_dM2"], mref
        )

    # Minimum signed curvature in the paper-relevant mass range.
    mask = (
        np.isfinite(model["kappa_R"])
        & (model["M"] >= MIN_CURVATURE_MASS)
        & (model["M"] <= min(2.0, model["Mmax"] - MMAX_DERIVATIVE_BUFFER))
    )

    if np.any(mask):
        idx = np.where(mask)[0][np.argmin(model["kappa_R"][mask])]
        row["kappa_R_min_1to2"] = model["kappa_R"][idx]
        row["M_at_kappa_R_min_Msun"] = model["M"][idx]
    else:
        row["kappa_R_min_1to2"] = np.nan
        row["M_at_kappa_R_min_Msun"] = np.nan

    summary_rows.append(row)

summary_path = OUTPUT_DIR / "paper_comparison_summary.csv"

fieldnames = list(summary_rows[0].keys())

with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(summary_rows)


# =====================================================================
# Figure 1: Mass-radius curves
# =====================================================================

fig, ax = plt.subplots(figsize=(6.7, 4.8))

for model in models:
    meta = model["meta"]
    mask = model["M"] >= MIN_ANALYSIS_MASS

    ax.plot(
        model["R_km"][mask],
        model["M"][mask],
        color=family_colors[meta["family"]],
        ls=curve_linestyle(model),
        marker=mg_markers.get(meta["Mg"], "o"),
        markevery=max(1, np.count_nonzero(mask) // 8),
        ms=2.5,
        lw=1.2,
    )

ax.set_xlabel(r"$R\ [{\rm km}]$")
ax.set_ylabel(r"$M/M_\odot$")
ax.set_title("Mass-radius relations")
add_legends(ax, models, family_colors, mg_markers)
save_figure(fig, "01_mass_radius_all")


# =====================================================================
# Figure 2: Dimensionless Lambda(M)
# =====================================================================

fig, ax = plt.subplots(figsize=(6.7, 4.8))

for model in models:
    meta = model["meta"]
    mask = (
        (model["M"] >= MIN_ANALYSIS_MASS)
        & np.isfinite(model["Lambda"])
        & (model["Lambda"] > 0)
    )

    ax.plot(
        model["M"][mask],
        model["Lambda"][mask],
        color=family_colors[meta["family"]],
        ls=curve_linestyle(model),
        marker=mg_markers.get(meta["Mg"], "o"),
        markevery=max(1, np.count_nonzero(mask) // 8),
        ms=2.5,
        lw=1.2,
    )

ax.set_xlabel(r"$M/M_\odot$")
ax.set_ylabel(r"$\Lambda$")
ax.set_yscale("log")
ax.set_title("Dimensionless tidal deformability")
add_legends(ax, models, family_colors, mg_markers)
save_figure(fig, "02_Lambda_vs_mass_all")


# =====================================================================
# Figure 3: dimensional/geometric lambda(M), paper-like
# =====================================================================

fig, ax = plt.subplots(figsize=(6.7, 4.8))

for model in models:
    meta = model["meta"]
    mask = (
        (model["M"] >= MIN_ANALYSIS_MASS)
        & np.isfinite(model["lambda_geom"])
        & (model["lambda_geom"] > 0)
    )

    ax.plot(
        model["M"][mask],
        model["lambda_geom"][mask],
        color=family_colors[meta["family"]],
        ls=curve_linestyle(model),
        marker=mg_markers.get(meta["Mg"], "o"),
        markevery=max(1, np.count_nonzero(mask) // 8),
        ms=2.5,
        lw=1.2,
    )

ax.set_xlabel(r"$M/M_\odot$")
ax.set_ylabel(r"$\lambda\ [M_\odot^5]$")
ax.set_title(r"Geometric tidal deformability $\lambda(M)$")
add_legends(ax, models, family_colors, mg_markers)
save_figure(fig, "03_lambda_vs_mass_all")


# =====================================================================
# Figure 4: signed curvature kappa_R(M), paper Fig. 2 style
# =====================================================================

fig, ax = plt.subplots(figsize=(6.7, 4.8))

for model in models:
    meta = model["meta"]
    mask = derivative_plot_mask(model) & np.isfinite(model["kappa_R"])

    ax.plot(
        model["M"][mask],
        model["kappa_R"][mask],
        color=family_colors[meta["family"]],
        ls=curve_linestyle(model),
        marker=mg_markers.get(meta["Mg"], "o"),
        markevery=max(1, np.count_nonzero(mask) // 8),
        ms=2.5,
        lw=1.2,
    )

# Paper-discussion reference levels; these are not universal boundaries.
ax.axhline(-1.5, color="0.5", lw=0.8, ls=":")
ax.axhline(-2.5, color="0.5", lw=0.8, ls=":")

ax.set_xlabel(r"$M/M_\odot$")
ax.set_ylabel(r"$\kappa_R$")
ax.set_title(r"Signed curvature of $R(M)$")
add_legends(ax, models, family_colors, mg_markers)
save_figure(fig, "04_kappa_R_vs_mass")


# =====================================================================
# Figure 5: kappa_R(Mref) vs R(Mref), paper Fig. 3 style
# =====================================================================

fig, ax = plt.subplots(figsize=(6.2, 4.8))

for model in models:
    mref = PAPER_REFERENCE_MASS

    Rref = interp_at_mass(model["M"], model["R_km"], mref)
    kref = interp_at_mass(model["M"], model["kappa_R"], mref)

    if not (np.isfinite(Rref) and np.isfinite(kref)):
        continue

    meta = model["meta"]

    marker = "*" if model["Mmax"] < MASS_CONSTRAINT else mg_markers.get(meta["Mg"], "o")

    ax.scatter(
        Rref,
        kref,
        color=family_colors[meta["family"]],
        marker=marker,
        s=42,
        label=meta["label"],
    )

ax.set_xlabel(
    rf"$R_{{{PAPER_REFERENCE_MASS:.1f}}}\ [{{\rm km}}]$"
)
ax.set_ylabel(
    rf"$\kappa_R(M={PAPER_REFERENCE_MASS:.1f}\,M_\odot)$"
)
ax.set_title("Paper-style curvature diagnostic")
ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    fontsize=6.5,
    frameon=True,
)
save_figure(fig, "05_kappa_ref_vs_R_ref")


# =====================================================================
# Figure 6: d^2 lambda / dM^2 vs M, paper Fig. 4 style
# =====================================================================

fig, ax = plt.subplots(figsize=(6.7, 4.8))

for model in models:
    meta = model["meta"]
    mask = (
        derivative_plot_mask(model)
        & np.isfinite(model["d2lambda_dM2"])
    )

    ax.plot(
        model["M"][mask],
        model["d2lambda_dM2"][mask],
        color=family_colors[meta["family"]],
        ls=curve_linestyle(model),
        marker=mg_markers.get(meta["Mg"], "o"),
        markevery=max(1, np.count_nonzero(mask) // 8),
        ms=2.5,
        lw=1.2,
    )

ax.set_xlabel(r"$M/M_\odot$")
ax.set_ylabel(
    r"$d^2\lambda/dM^2\ [M_\odot^3]$"
)
ax.set_title(r"Second derivative of $\lambda(M)$")
add_legends(ax, models, family_colors, mg_markers)
save_figure(fig, "06_d2lambda_dM2_vs_mass")


# =====================================================================
# Figure 7: d2lambda(Mref) vs lambda(Mref), paper Fig. 5 style
# =====================================================================

fig, ax = plt.subplots(figsize=(6.2, 4.8))

for model in models:
    mref = PAPER_REFERENCE_MASS

    lref = interp_at_mass(model["M"], model["lambda_geom"], mref)
    d2ref = interp_at_mass(model["M"], model["d2lambda_dM2"], mref)

    if not (np.isfinite(lref) and np.isfinite(d2ref)):
        continue

    meta = model["meta"]
    marker = "*" if model["Mmax"] < MASS_CONSTRAINT else mg_markers.get(meta["Mg"], "o")

    ax.scatter(
        lref,
        d2ref,
        color=family_colors[meta["family"]],
        marker=marker,
        s=42,
        label=meta["label"],
    )

ax.set_xlabel(
    rf"$\lambda(M={PAPER_REFERENCE_MASS:.1f}\,M_\odot)"
    r"\ [M_\odot^5]$"
)
ax.set_ylabel(
    rf"$d^2\lambda/dM^2"
    rf"(M={PAPER_REFERENCE_MASS:.1f}\,M_\odot)"
    r"\ [M_\odot^3]$"
)
ax.set_title("Paper-style tidal second-derivative diagnostic")
ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    fontsize=6.5,
    frameon=True,
)
save_figure(fig, "07_d2lambda_ref_vs_lambda_ref")


# =====================================================================
# Figure 8: dR/dM(Mref) vs Mmax, paper Fig. 6 top style
# =====================================================================

fig, ax = plt.subplots(figsize=(6.2, 4.8))

for model in models:
    if not model["mmax_resolved"]:
        continue

    slope = interp_at_mass(
        model["M"], model["dR_dM"], PAPER_REFERENCE_MASS
    )

    if not np.isfinite(slope):
        continue

    meta = model["meta"]
    marker = "*" if model["Mmax"] < MASS_CONSTRAINT else mg_markers.get(meta["Mg"], "o")

    ax.scatter(
        model["Mmax"],
        slope,
        color=family_colors[meta["family"]],
        marker=marker,
        s=42,
        label=meta["label"],
    )

ax.set_xlabel(r"$M_{\max}/M_\odot$")
ax.set_ylabel(
    rf"$dR/dM(M={PAPER_REFERENCE_MASS:.1f}\,M_\odot)$"
)
ax.set_title("Radius slope versus maximum mass")
ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    fontsize=6.5,
    frameon=True,
)
save_figure(fig, "08_dR_dM_ref_vs_Mmax")


# =====================================================================
# Figure 9: dLambda/dM(Mref) vs Mmax, paper Fig. 6 bottom style
# =====================================================================

fig, ax = plt.subplots(figsize=(6.2, 4.8))

for model in models:
    if not model["mmax_resolved"]:
        continue

    slope = interp_at_mass(
        model["M"], model["dLambda_dM"], PAPER_REFERENCE_MASS
    )

    if not np.isfinite(slope):
        continue

    meta = model["meta"]
    marker = "*" if model["Mmax"] < MASS_CONSTRAINT else mg_markers.get(meta["Mg"], "o")

    ax.scatter(
        model["Mmax"],
        slope,
        color=family_colors[meta["family"]],
        marker=marker,
        s=42,
        label=meta["label"],
    )

ax.set_xlabel(r"$M_{\max}/M_\odot$")
ax.set_ylabel(
    rf"$d\Lambda/dM(M={PAPER_REFERENCE_MASS:.1f}\,M_\odot)$"
)
ax.set_title("Dimensionless tidal slope versus maximum mass")
ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    fontsize=6.5,
    frameon=True,
)
save_figure(fig, "09_dLambda_dM_ref_vs_Mmax")


# =====================================================================
# Warnings / diagnostics
# =====================================================================

warnings_path = OUTPUT_DIR / "analysis_warnings.txt"

with warnings_path.open("w") as handle:
    if warnings:
        handle.write("\n".join(warnings))
        handle.write("\n")
    else:
        handle.write("No automatic warnings.\n")

print("\n============================================================")
print("Multi-EoS paper-oriented analysis complete")
print("============================================================")
print(f"Analysed EoSs : {len(models)}")
print(f"Summary CSV   : {summary_path}")
print(f"Figures       : {OUTPUT_DIR}")
print(f"Warnings      : {warnings_path}")
print()
print("Paper-comparison notes:")
print(f"  Reference mass = {PAPER_REFERENCE_MASS:.1f} Msun")
print("  Signed curvature is used (negative values are meaningful).")
print("  Dashed curves have Mmax < 2.0 Msun, matching the paper.")
print("  Dotted curves indicate Mmax was not reached by the run.")
print("  lambda = Lambda * M^5 is analysed separately from Lambda.")
print("============================================================")
