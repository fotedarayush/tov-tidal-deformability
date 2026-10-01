# Offline workflow: 100 hyperonic EoSs -> TOV/tidal runs

This guide assumes you are at the repository root:
`tov-tidal-deformability/`

## A. Do these while you still have internet

1. Make sure these two source datasets are saved locally:
   - `hyperonic_eos_at_zero_temperature.csv.gz`
   - `hyperonic_mr_curve_at_zero_temperature.csv.gz`

2. Download the CompOSE crust table:
   - Page: https://compose.obspm.fr/eos/212
   - Direct archive currently linked by CompOSE:
     https://compose.obspm.fr/download/1D/Crust/TW/eos.zip

3. Put/extract it as:
   `eos/crust/GPPVA_TW/eos.nb`
   `eos/crust/GPPVA_TW/eos.thermo`
   etc.

4. Make sure Python packages already exist:
   `python3 -c "import numpy,pandas,matplotlib; print('OK')"`

5. Make sure gfortran works:
   `gfortran --version`

The CompOSE GPPVA(TW) crust table uses BPS below nB=0.002 fm^-3,
so its region around nB=1e-4 fm^-3 is the required BPS outer crust.

The hyperonic-model prescription being followed is:
BPS outer crust -> P(epsilon)=a1+a2 epsilon^(4/3) -> core at nB=0.04 fm^-3.

## B. Convert the CompOSE crust

Create the directory if needed:
`mkdir -p eos/crust/GPPVA_TW`

Extract the downloaded eos.zip there.

Then run:
```
python3 scripts/convert_compose_1d.py \
  eos/crust/GPPVA_TW \
  eos/crust/GPPVA_TW/GPPVA_TW_full.dat
```

Inspect:
```
head eos/crust/GPPVA_TW/GPPVA_TW_full.dat
```

## C. Put the workflow scripts into your repo

Copy:
- `build_complete_hyperonic_eos.py` -> `scripts/`
- `validate_complete_tables.py` -> `scripts/`
- `validate_id130.py` -> `scripts/`
- `run_hyperonic_100.sh` -> `scripts/`
- `audit_hyperonic_runs.py` -> `analysis/`

Make the batch runner executable:
`chmod +x scripts/run_hyperonic_100.sh`

## D. Place the reference M-R file

Put:
`hyperonic_mr_curve_at_zero_temperature.csv.gz`

at:
`eos/hyperonic_dataset/hyperonic_mr_curve_at_zero_temperature.csv.gz`

Your existing 100 core files should be in:
`eos/hyperonic_dataset/hyperonic_100/core_tables/`

Check:
```
ls eos/hyperonic_dataset/hyperonic_100/core_tables/*.dat | wc -l
```
Expected: `100`

## E. Build all 100 complete crust-to-core EoSs

Run:
```
python3 scripts/build_complete_hyperonic_eos.py
```

Expected output folder:
`eos/hyperonic_dataset/hyperonic_100/complete_tables/`

Then validate:
```
python3 scripts/validate_complete_tables.py
```

Expected:
- Number of complete tables: 100
- Bad tables: 0

Inspect crust-matching diagnostics:
```
python3 - <<'PY'
import pandas as pd
p="eos/hyperonic_dataset/hyperonic_100/crust_matching_diagnostics.csv"
d=pd.read_csv(p)
print(d[["core_file","n_endpoint_correction_percent",
         "cs2_min_bridge","cs2_max_bridge"]].head())
print()
print(d["n_endpoint_correction_percent"].describe())
PY
```

Do not continue if the script reports non-monotonic P(epsilon), negative
pressure, or dP/depsilon >= 1 in the bridge.

## F. Compile the optimised solver

From repository root:
```
gfortran -O1 \
  -fno-automatic \
  -o logtov_seq_geom_tidal.out \
  src/logtov_seq_geom_tidal.f90
```

Do NOT use -O2/-O3 for this project unless you redo the tidal regression
tests. Do NOT increase the radial step h for production runs without a
convergence study.

## G. Test ONE EoS first: ID 130

Make:
```
mkdir -p results/hyperonic_100/test_id130
```

Set `infile`:
```
1.e14 0.05e14 400
eos/hyperonic_dataset/hyperonic_100/complete_tables/hyperon_EOS_001_ID_00130_complete.dat
```

Run:
```
/usr/bin/time -p ./logtov_seq_geom_tidal.out
```

Save:
```
cp logtov_seq_geom_tidal.dat \
results/hyperonic_100/test_id130/hyperon_EOS_001_ID_00130_complete_TOV.dat
```

Validate against the authors' supplied zero-T M-R sequence:
```
python3 scripts/validate_id130.py
```

Reference values previously extracted for ID 130:
- Mmax ~ 2.03494 Msun
- R(Mmax) ~ 11.9261 km
- R1.4 ~ 13.3161 km
- R1.6 ~ 13.3148 km

Use the full curve/RMS comparison as well as individual values. If the
radius disagreement is large (for example several tenths of a km or more),
stop and inspect the crust construction/solver before running all 100.

## H. Run the 100-EoS batch

Start the batch:
```
./scripts/run_hyperonic_100.sh
```

The script:
- runs each EoS sequentially;
- starts with 400 central-density models;
- automatically retries with 600 if Mmax is not resolved;
- saves each TOV output separately;
- saves logs;
- is resumable;
- restores your original `infile` on exit.

Outputs:
`results/hyperonic_100/coarse/`

Logs:
`results/hyperonic_100/logs/`

Status:
`results/hyperonic_100/batch_status.csv`

If the computer sleeps, the process can pause. On macOS you can keep it
awake by running:
```
caffeinate -i ./scripts/run_hyperonic_100.sh
```

## I. Audit the 100 outputs

Run:
```
python3 analysis/audit_hyperonic_runs.py
```

Expected:
- 100 run files
- 100 resolved Mmax, ideally
- summary at:
  `results/hyperonic_100/coarse_summary.csv`

Inspect unresolved models:
```
cat results/hyperonic_100/unresolved.txt
```

If this file is empty, the coarse batch reached a turning point for all
100 models.

## J. What NOT to do yet

Do not use these coarse sequences for final Bauswein curvature or
d2lambda/dM2 plots.

The Bauswein-style quantities include second derivatives, and second
derivatives amplify numerical noise. After the 100 coarse runs are
validated, the next stage is a fine central-density scan on the stable
branch (especially around ~1.4--2.0 Msun) followed by convergence tests.

The coarse runs are for:
- Mmax
- R1.4, R1.6
- Lambda1.4, Lambda1.6
- checking the overall M-R and Lambda-M sequences
- locating the central-density region for later fine runs
