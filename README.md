# TOV and Tidal-Deformability Solver

A Fortran/Python workflow for computing neutron-star structure and tidal properties from tabulated equations of state (EoSs), with a current focus on cold hyperonic matter and the stellar diagnostics discussed by Bauswein et al. (2026).

The repository contains a Fortran TOV + tidal-deformability solver, EoS conversion and validation tools, a reproducibly selected 100-model hyperonic ensemble, crust-to-core construction utilities, batch runners, validation scripts, paper-oriented analysis scripts, and PDF documentation.

---

## 1. Scientific scope

For each central energy density, the solver integrates the stellar-structure and tidal-perturbation equations and returns quantities including

- gravitational mass `M`,
- radius `R`,
- compactness `C=M/R`,
- quadrupolar Love number `k2`,
- tidal coupling quantity `kappa2`,
- dimensionless tidal deformability `Lambda`,
- rest/baryonic mass and auxiliary quantities.

The current research workflow studies whether stellar observables can reveal strong softening associated with hyperons in neutron-star matter.

The main bulk quantities are

`Mmax`, `R1.4`, `R1.6`, `Lambda1.4`, and `Lambda1.6`.

The paper-oriented derivative quantities include the signed mass-radius curvature

`kappa_R = (d2R/dM2) / [1 + (dR/dM)^2]^(3/2)`

together with `d2lambda/dM2`, `dR/dM` at 1.6 solar masses, and `dLambda/dM` at 1.6 solar masses.

Reference paper:

> A. Bauswein et al., *Stellar properties indicating the presence of hyperons in neutron stars*, Phys. Rev. Research **8**, 013253 (2026), DOI: 10.1103/ygtr-ktqk.

The current 100-EoS hyperonic ensemble is an additional dataset on which this methodology is being applied; it is **not** the exact EoS sample used by Bauswein et al.

---

## 2. Repository structure

```text
tov-tidal-deformability/
├── README.md
├── infile
│
├── src/
│   └── logtov_seq_geom_tidal.f90
│
├── scripts/
│   ├── Tidal_control.sh
│   ├── run_all_eos.sh
│   ├── convert_compose_1d.py
│   ├── extract_100_hyperonic.py
│   ├── build_complete_hyperonic_eos.py
│   ├── validate_complete_tables.py
│   ├── validate_id130.py
│   └── run_hyperonic_100.sh
│
├── analysis/
│   ├── analyse_single_eos.py
│   ├── analyse_multiple_eos.py
│   ├── analyse_paper_comparison.py
│   └── audit_hyperonic_runs.py
│
├── eos/
│   ├── lists/
│   ├── EoS_dd2_npY_T0_beta_eq_eta_D-eta_V-Mg/
│   ├── compose/
│   ├── crust/
│   │   └── GPPVA_TW/
│   └── hyperonic_dataset/
│       └── hyperonic_100/
│           ├── core_tables/
│           ├── complete_tables/
│           ├── manifest.csv
│           ├── selected_ids.txt
│           └── crust_matching_diagnostics.csv
│
├── results/
│   └── hyperonic_100/
│
├── docs/
│   └── instructions.txt
│
└── documentations/
    └── PDF documentation files
```

`documentations/` contains the human-readable PDF documentation for the project. New users should read those files alongside this README.

---

## 3. Requirements

### Fortran

A working GNU Fortran compiler is required.

```bash
gfortran --version
```

### Python

Python 3 is used for conversion, validation, batch preparation, and analysis.

Main packages:

```text
numpy
pandas
matplotlib
scienceplots
```

Quick check:

```bash
python3 -c "import numpy, pandas, matplotlib; print('Python dependencies OK')"
```

For the paper-comparison plotting script:

```bash
python3 -m pip install SciencePlots
```

---

## 4. EoS input format

The Fortran solver expects a four-column text table:

```text
epsilon    P    nB    muB
```

with units

```text
epsilon : MeV fm^-3
P       : MeV fm^-3
nB      : fm^-3
muB     : MeV
```

A header line is required because the current Fortran reader discards the first line before reading numerical data.

Example:

```text
# epsilon[MeV/fm^3] P[MeV/fm^3] nB[fm^-3] muB[MeV]
3.7897871000000000e+01  1.5324300000000000e-01  4.0000000000000000e-02  9.5127970100000000e+02
...
```

The current interpolation routines operate in logarithmic space. The optimised EoS search also assumes monotonic tables, so newly generated EoSs should be validated before use.

---

## 5. The `infile`

The solver reads run parameters from `infile` in the repository root.

Format:

```text
<initial central energy density> <central-density step> <number of models>
<path to EoS table>
```

Example:

```text
1.e14 0.05e14 400
eos/hyperonic_dataset/hyperonic_100/complete_tables/hyperon_EOS_001_ID_00130_complete.dat
```

Each stellar model corresponds to one central density in the generated sequence.

---

## 6. Compiling the solver

Recommended production build:

```bash
gfortran -O1   -fno-automatic   -o logtov_seq_geom_tidal.out   src/logtov_seq_geom_tidal.f90
```

`-O1` was retained after regression testing. Higher compiler optimisation levels were not adopted for production because benchmark tests changed a tidal output quantity at roughly the percent level.

The executable is machine dependent and should normally not be committed.

---

## 7. Running one EoS

After preparing `infile`:

```bash
./logtov_seq_geom_tidal.out
```

or

```bash
/usr/bin/time -p ./logtov_seq_geom_tidal.out
```

The main output is

```text
logtov_seq_geom_tidal.dat
```

This file is overwritten by the next solver run, so save it first.

Example:

```bash
mkdir -p results/test

cp logtov_seq_geom_tidal.dat   results/test/example_TOV.dat
```

### Main output columns

| Column | Quantity |
|---|---|
| 1 | central total energy density in solver units |
| 2 | radius |
| 3 | gravitational mass in solar masses |
| 4 | rest/baryonic mass |
| 7 | compactness |
| 8 | central rest-mass density |
| 9 | central pressure |
| 10 | Love number `k2` |
| 11 | tidal coupling quantity `kappa2` |
| 12 | dimensionless tidal deformability `Lambda` |
| 13+ | tidal/auxiliary quantities |

**Important:** in the current source, output column 2 is already written in kilometres. Do not multiply it by 1.476 again.

Low-mass tidal values can become numerically unreliable and should not automatically be interpreted as physical results.

---

## 8. Solver optimisation

The production solver has been accelerated without changing the physical equations, RK4 scheme, radial integration step, or logarithmic EoS interpolation rule.

Accepted changes include

- restricting searches to the actually loaded EoS range;
- replacing repeated linear EoS searches with binary searches;
- caching repeated `epsilon(P)` evaluations inside RHS routines;
- removing unnecessary diagnostic I/O from production runs;
- caching exact repeated pressure requests.

A representative 20-model benchmark decreased from approximately `11.46 s` to `1.63 s`, corresponding to roughly a `7x` speed-up.

See [`documentations/`](documentations/) for the full optimisation record.

---

## 9. Hyperonic 100-EoS dataset

The project includes a reproducibly selected set of 100 cold hyperonic EoSs from the larger zero-temperature hyperonic ensemble used in the current study.

Selection script:

```text
scripts/extract_100_hyperonic.py
```

Fixed NumPy seed:

```text
20260930
```

The selected source IDs are stored in

```text
eos/hyperonic_dataset/hyperonic_100/selected_ids.txt
```

and per-model information in

```text
eos/hyperonic_dataset/hyperonic_100/manifest.csv
```

### Core tables

The selected source EoSs begin at

`nB = 0.04 fm^-3`.

Therefore the files in

```text
eos/hyperonic_dataset/hyperonic_100/core_tables/
```

are **core-only** tables and should not be used directly as complete neutron-star EoSs.

---

## 10. Building complete crust-to-core EoSs

Use

```text
scripts/build_complete_hyperonic_eos.py
```

The current low-density prescription is

```text
BPS outer crust
        ↓
inner-crust polytropic bridge
        ↓
hyperonic core
```

with

`P(epsilon) = a1 + a2 * epsilon^(4/3)`.

The outer-crust endpoint is taken near

`nB = 1e-4 fm^-3`

and the hyperonic core begins at

`nB = 0.04 fm^-3`.

The matching coefficients are determined separately for each hyperonic core so that pressure is continuous at both endpoints.

The BPS outer-crust segment used by the current workflow is obtained from the CompOSE GPPVA(TW) crust table and converted into the solver format.

Build the 100 complete EoSs with

```bash
python3 scripts/build_complete_hyperonic_eos.py
```

Output:

```text
eos/hyperonic_dataset/hyperonic_100/complete_tables/
```

Then validate:

```bash
python3 scripts/validate_complete_tables.py
```

Expected result:

```text
Number of complete tables: 100
Bad tables: 0
```

---

## 11. Validation against reference M-R sequences

Representative models were compared against the supplied zero-temperature mass-radius reference sequences before launching the full 100-EoS calculation.

| EoS ID | Reference Mmax | Local Mmax | R1.4 difference | M-R radius RMS |
|---:|---:|---:|---:|---:|
| 130 | 2.03494 | 2.02890 | -0.170 km | 0.172 km |
| 230 | 2.03359 | 2.02758 | -0.157 km | 0.161 km |
| 15732 | 2.21160 | 2.20510 | -0.189 km | 0.185 km |

For these tests, maximum masses agree at about the 0.3% level. The local radii are systematically slightly smaller than the supplied reference sequences, by approximately 0.15-0.19 km in the tested cases.

Dedicated ID 130 validation:

```bash
python3 scripts/validate_id130.py
```

These comparisons validate the present workflow for the coarse population study while documenting residual sensitivity to the low-density/crust treatment.

---

## 12. Running the 100-EoS batch

Make the batch runner executable:

```bash
chmod +x scripts/run_hyperonic_100.sh
```

Run:

```bash
./scripts/run_hyperonic_100.sh
```

On macOS:

```bash
caffeinate -i ./scripts/run_hyperonic_100.sh
```

The first pass uses

```text
1.e14 0.05e14 400
```

and automatically extends a model to 600 configurations if the maximum-mass turning point is not adequately resolved.

Outputs:

```text
results/hyperonic_100/coarse/
```

Logs:

```text
results/hyperonic_100/logs/
```

Status:

```text
results/hyperonic_100/batch_status.csv
```

The batch runner is designed to skip already completed EoSs with resolved maxima when restarted.

---

## 13. Auditing the batch

After the coarse batch finishes:

```bash
python3 analysis/audit_hyperonic_runs.py
```

This generates

```text
results/hyperonic_100/coarse_summary.csv
```

with quantities including

- `Mmax`,
- `R(Mmax)`,
- `R1.4`,
- `R1.6`,
- `Lambda1.4`,
- `Lambda1.6`,
- whether the maximum-mass turning point is resolved.

---

## 14. Paper-oriented analysis

Main script:

```text
analysis/analyse_paper_comparison.py
```

It can generate

- combined M-R curves;
- `Lambda(M)`;
- dimensional/geometric `lambda(M)`;
- signed curvature `kappa_R(M)`;
- `kappa_R` at a reference mass;
- `d2lambda/dM2`;
- `dR/dM` at a reference mass;
- `dLambda/dM` at a reference mass;
- CSV summaries and warning files.

The script uses centred second-order finite differences for nonuniform mass spacing.

### Resolution warning

The coarse 400/600-model sequences are suitable for population-level quantities such as `Mmax`, radii and reference-mass `Lambda`, but should **not automatically be treated as converged for second derivatives**.

Publication-quality derivative work should use

1. a finer central-density grid;
2. resolution/convergence tests;
3. careful treatment of the region near `Mmax`.

A purely nucleonic control ensemble is also required for a direct hyperonic-vs-nucleonic reproduction of the Bauswein comparison.

---

## 15. Numerical caveats

### Monotonicity

The optimised binary-search interpolation assumes monotonic EoS arrays. Validate new tables before running the solver.

### Header handling

The Fortran reader discards the first line of each EoS file. Solver-ready tables should contain a header.

### Radial step

The current source uses

```fortran
h = 10.*1e-5/1.746
```

Changing it affects runtime and numerical accuracy. Do not alter it for precision work without a convergence test.

### Compiler optimisation

The production workflow uses `-O1`. Do not switch to `-O2`, `-O3`, or aggressive floating-point options without repeating the tidal regression tests.

### Radius convention

For the current source, output column 2 is already in kilometres.

### Derivative observables

Second derivatives are much more sensitive to sequence resolution than `M`, `R`, or `Lambda`. Always test convergence before interpreting curvature or `d2lambda/dM2`.

### Output overwriting

`logtov_seq_geom_tidal.dat` is replaced on every run. Save each result before running another EoS.

---

## 16. Documentation

Detailed project documentation is stored in

[`documentations/`](documentations/)

and covers topics such as

- TOV/tidal solver architecture;
- variable definitions and input/output conventions;
- EoS interpolation and units;
- RK4 integration and tidal perturbations;
- solver optimisation and benchmarking;
- hyperonic neutron-star physics and analysis methodology;
- numerical caveats and suggested future improvements.

The original solver instructions are retained under

```text
docs/
```

for provenance.

---

## 17. Recommended workflow for a new user

```text
1. Read README.md and documentations/
2. Check Python and gfortran
3. Compile the solver with -O1 -fno-automatic
4. Inspect and validate the EoS table
5. Prepare infile
6. Run one EoS
7. Save the output
8. Check M-R behaviour and Mmax resolution
9. Only then run multi-EoS calculations
10. Audit the population results
11. Use fine grids + convergence tests for derivative observables
```

For the current hyperonic workflow:

```text
core tables
    ↓
build complete crust+core tables
    ↓
validate all tables
    ↓
validate representative EoSs against reference M-R curves
    ↓
run the 100-EoS coarse batch
    ↓
audit Mmax, R and Lambda
    ↓
perform fine-resolution derivative runs
    ↓
Bauswein-style curvature/tidal analysis
```

---

## 18. Reproducibility

The project deliberately separates

- source code,
- small reproducible EoS subsets and metadata,
- generated outputs,
- large external source datasets.

Large external archives and generated batch outputs should generally not be committed directly to Git. The scripts, selected source IDs, manifests, and construction procedure are intended to make the working dataset reproducible.

When adding a new EoS family, document

- source and reference;
- physical composition;
- units and column definitions;
- crust treatment;
- preprocessing;
- monotonicity checks;
- central-density sampling;
- compiler and solver settings.
