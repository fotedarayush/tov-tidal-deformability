# TOV and Tidal Deformability Analysis for Hyperonic Neutron-Star Equations of State

This repository contains a Fortran implementation of the Tolman-Oppenheimer-Volkoff (TOV) equations together with tidal-deformability calculations for cold, beta-equilibrated neutron-star equations of state.

The current project focuses on families of **hyperonic DD2-based equations of state**, with automated batch calculations and Python analysis of neutron-star structure, tidal deformability, and mass-radius curvature.

The analysis is also designed to reproduce and investigate several diagnostics discussed in recent work on identifying hyperonic matter through neutron-star macroscopic observables.

---

## Physics

For a static, spherically symmetric, non-rotating neutron star, the stellar structure is determined by the TOV equations

```math
\frac{dm}{dr} = 4\pi r^2 \epsilon
```

and

```math
\frac{dP}{dr}
=
-\frac{(\epsilon+P)(m+4\pi r^3P)}
{r(r-2m)}.
```

The stellar surface is defined approximately by

```math
P(R) \simeq 0,
```

with total gravitational mass

```math
M = m(R).
```

The code simultaneously integrates the tidal perturbation equations required to determine the quadrupolar Love number $k_2$.

The dimensionless tidal deformability is

```math
\Lambda
=
\frac{2}{3}k_2 C^{-5},
```

where the compactness is

```math
C = \frac{GM}{Rc^2}.
```

The geometric tidal deformability used in part of the analysis is

```math
\lambda = \Lambda M^5.
```

---

## Repository structure

```text
tov-tidal-deformability/
│
├── README.md
├── infile
│
├── src/
│   └── logtov_seq_geom_tidal.f90
│
├── analysis/
│   ├── analyse_single_eos.py
│   ├── analyse_multiple_eos.py
│   └── analyse_paper_comparison.py
│
├── scripts/
│   ├── Tidal_control.sh
│   └── run_all_eos.sh
│
├── eos/
│   ├── lists/
│   └── EoS_dd2_npY_T0_beta_eq_eta_D-eta_V-Mg/
│
├── docs/
│   └── instructions.txt
│
└── results/
    ├── tov/
    └── paper_comparison/
```

---

## Equations of state

The current calculations use cold, beta-equilibrated hyperonic EoSs with filenames of the form

```text
EoS_dd2_npY_T0_beta_eq_01_05_Mg_800.dat
```

The `npY` label denotes neutron-proton-hyperon matter.

In the current plotting and analysis convention, a filename such as

```text
01_05_Mg_800
```

is displayed as

```math
\eta_D = 0.1,
```

```math
\eta_V = 0.5,
```

and

```math
M_g = 800.
```

Different combinations of these parameters generate different hyperonic EoSs.

---

## Input format

The Fortran code reads its configuration from `infile`.

Example:

```text
1.e14 0.3e14 100
eos/EoS_dd2_npY_T0_beta_eq_eta_D-eta_V-Mg/EoS_dd2_npY_T0_beta_eq_01_05_Mg_800.dat
```

The first line specifies:

```text
initial central energy density
central-density step
number of stellar models
```

The second line specifies the path to the EoS table.

For higher-resolution derivative and curvature calculations, a finer sequence can be used, for example

```text
1.e14 0.05e14 350
```

---

## Compilation

The Fortran code can be compiled using `gfortran`:

```bash
gfortran -o logtov_seq_geom_tidal.out \
src/logtov_seq_geom_tidal.f90 \
-fno-automatic
```

---

## Running a single EoS

After editing `infile`, run

```bash
./logtov_seq_geom_tidal.out
```

The main output file is

```text
logtov_seq_geom_tidal.dat
```

A typical output contains quantities including:

- central energy density,
- radius,
- gravitational mass,
- baryonic/rest mass,
- compactness,
- central pressure,
- Love number $k_2$,
- tidal deformability $\Lambda$,
- tidal perturbation quantities.

In the current source, the radius written to column 2 is already in km.

---

## Batch EoS calculations

Multiple EoSs can be run automatically using

```bash
./scripts/run_all_eos.sh
```

The script reads the EoS filenames from a list and automatically:

1. rewrites `infile`,
2. runs the Fortran solver,
3. saves the output under a unique EoS-specific filename,
4. repeats the calculation for the next EoS.

Results are stored in

```text
results/tov/
```

for example

```text
results/tov/EoS_dd2_npY_T0_beta_eq_01_05_Mg_800_TOV.dat
```

---

## Python analysis

### Single-EoS analysis

Single-EoS outputs can be analysed to obtain quantities such as

```math
M_{\max},
```

```math
R_{1.4},
```

and

```math
\Lambda_{1.4}.
```

For a standard one-family TOV sequence, the stable branch is taken up to the maximum of

```math
M(\rho_c).
```

---

### Multi-EoS analysis

The multi-EoS analysis automatically compares all calculated sequences and extracts quantities such as

```math
M_{\max},
```

```math
R_{1.4},
\qquad
R_{1.6},
\qquad
R_{1.8},
```

and

```math
\Lambda_{1.4},
\qquad
\Lambda_{1.6},
\qquad
\Lambda_{1.8}.
```

It also produces combined mass-radius and tidal-deformability plots.

---

## Paper-oriented analysis

The script

```text
analysis/analyse_paper_comparison.py
```

implements additional diagnostics motivated by studies of hyperonic signatures in neutron-star observables.

The main quantities include the signed curvature of the mass-radius relation

```math
\kappa_R
=
\frac{d^2R/dM^2}
{\left[1+\left(dR/dM\right)^2\right]^{3/2}},
```

the geometric tidal deformability

```math
\lambda = \Lambda M^5,
```

and derivatives such as

```math
\frac{dR}{dM},
```

```math
\frac{d\Lambda}{dM},
```

and

```math
\frac{d^2\lambda}{dM^2}.
```

The analysis generates plots including

```text
01_mass_radius_all.pdf
02_Lambda_vs_mass_all.pdf
03_lambda_vs_mass_all.pdf
04_kappa_R_vs_mass.pdf
05_kappa_ref_vs_R_ref.pdf
06_d2lambda_dM2_vs_mass.pdf
07_d2lambda_ref_vs_lambda_ref.pdf
08_dR_dM_ref_vs_Mmax.pdf
09_dLambda_dM_ref_vs_Mmax.pdf
```

as well as a summary CSV containing the derived stellar quantities.

---

## Current qualitative results

For the currently analysed hyperonic EoS family, many models share nearly identical intermediate-density neutron-star structure while differing more strongly close to their maximum masses.

Typical values for the well-behaved sequences are approximately

```math
R_{1.6} \simeq 13.15\ {\rm km},
```

```math
\Lambda_{1.6} \simeq 299,
```

with maximum masses around

```math
M_{\max} \sim 2.0 - 2.1\,M_\odot
```

for several parameter combinations.

The calculated mass-radius curvature also becomes significantly negative in the high-mass regime, with some models reaching approximately

```math
\kappa_R \lesssim -2.5.
```

This behaviour is of particular interest when comparing hyperonic and nucleonic neutron-star models.

These values are preliminary and should be interpreted together with numerical-resolution and stability checks.

---

## Numerical considerations

First and especially second derivatives are significantly more sensitive to numerical resolution than the basic mass-radius relation.

Quantities such as

```math
\frac{d^2R}{dM^2}
```

and

```math
\frac{d^2\lambda}{dM^2}
```

can become noisy if:

- the central-density grid is too coarse,
- neighbouring stellar models have nearly identical masses,
- the EoS contains sharp features,
- the maximum mass is not reached within the chosen density range.

For this reason, coarse runs are useful for initial exploration, while finer central-density sampling should be used for final curvature and derivative analysis.

The analysis scripts also flag EoSs for which the maximum mass occurs at the edge of the calculated sequence.

---

## Stable branch

For ordinary one-family TOV sequences, configurations are treated as stable up to the maximum-mass point, where approximately

```math
\frac{dM}{d\rho_c} = 0.
```

The normal stable branch satisfies approximately

```math
\frac{dM}{d\rho_c} > 0.
```

Beyond the maximum-mass turning point,

```math
\frac{dM}{d\rho_c} < 0,
```

and the configurations are generally associated with radial instability.

More complicated third-family or twin-star sequences require a dedicated stability analysis.

---

## Tidal deformability

The dimensionless tidal deformability depends strongly on compactness:

```math
\Lambda
=
\frac{2}{3}k_2
\left(
\frac{Rc^2}{GM}
\right)^5.
```

This means that relatively small differences in radius can produce much larger differences in $\Lambda$.

In general,

```math
M \uparrow
\quad\Rightarrow\quad
C \uparrow
\quad\Rightarrow\quad
\Lambda \downarrow.
```

This makes tidal deformability a particularly sensitive probe of the underlying EoS.

---

## Mass-radius curvature

The signed curvature

```math
\kappa_R
=
\frac{d^2R/dM^2}
{\left[1+\left(dR/dM\right)^2\right]^{3/2}}
```

provides a way to quantify changes in the shape of the mass-radius relation.

Rather than describing a curve only qualitatively as "bending", $\kappa_R$ allows that behaviour to be measured directly.

Strong changes in curvature can be associated with rapid changes in the stiffness or composition of the underlying EoS, although numerical convergence checks are required before interpreting such features physically.

---

## Plotting

The analysis scripts use `SciencePlots` for publication-style figures.

Install it with

```bash
python3 -m pip install SciencePlots
```

The plotting scripts use

```python
plt.style.use(["science", "no-latex", "grid"])
```

so a local LaTeX installation is not required.

---

## Python requirements

The main Python dependencies are:

- `numpy`
- `matplotlib`
- `SciencePlots`

Install them using

```bash
python3 -m pip install numpy matplotlib SciencePlots
```

---

## Analysis workflow

A typical workflow is:

```text
Select EoS
   ↓
Edit or generate infile
   ↓
Run TOV + tidal solver
   ↓
Save EoS-specific output
   ↓
Repeat for all EoSs
   ↓
Run multi-EoS Python analysis
   ↓
Extract Mmax, radii and tidal quantities
   ↓
Calculate derivatives and curvature
   ↓
Compare EoS families
```

For batch calculations:

```bash
./scripts/run_all_eos.sh
```

followed by

```bash
python3 analysis/analyse_paper_comparison.py
```

---

## Main observables

The primary quantities used when comparing EoSs are

```math
M_{\max},
```

```math
R_{1.4},
\qquad
R_{1.6},
```

and

```math
\Lambda_{1.4},
\qquad
\Lambda_{1.6}.
```

Higher-order diagnostics include

```math
\kappa_R(M),
```

```math
\frac{dR}{dM},
```

```math
\frac{d\Lambda}{dM},
```

and

```math
\frac{d^2\lambda}{dM^2}.
```

---

## Current project goals

Current work includes:

- systematic scans of hyperonic EoS parameters,
- automated multi-EoS TOV calculations,
- higher-resolution stellar sequences,
- mass-radius curvature analysis,
- tidal-deformability derivative analysis,
- identification of numerically unstable or unresolved sequences,
- comparison with published hyperonic neutron-star diagnostics,
- future comparison with corresponding purely nucleonic control EoSs.

---

## Reference comparison

The current paper-oriented analysis is motivated by recent work investigating whether macroscopic neutron-star observables may provide signatures of hyperonic degrees of freedom.

Particular attention is given to:

- the curvature of the mass-radius relation,
- derivatives of tidal deformability,
- correlations between fixed-mass observables and maximum mass.

The present calculations focus on hyperonic EoSs. A full hyperonic-versus-nucleonic comparison will require corresponding purely nucleonic control models analysed with the same numerical pipeline.

---

## Status

This repository is under active development.

The current results should be considered preliminary until convergence tests with finer central-density sampling have been completed, particularly for derivative-based observables.