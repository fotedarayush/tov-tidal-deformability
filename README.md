# TOV and Tidal Deformability Analysis for Hyperonic Neutron-Star Equations of State

This repository contains a Fortran implementation of the Tolman-Oppenheimer-Volkoff (TOV) equations together with tidal-deformability calculations for cold, beta-equilibrated neutron-star equations of state.

The current project focuses on families of **hyperonic DD2-based equations of state**, with automated batch calculations and Python analysis of neutron-star structure, tidal deformability, and mass-radius curvature.

The analysis is also designed to reproduce and investigate several diagnostics discussed in recent work on identifying hyperonic matter through neutron-star macroscopic observables.

---

## Physics

For a static, spherically symmetric, non-rotating neutron star, the stellar structure is determined by the TOV equations

\[
\frac{dm}{dr} = 4\pi r^2 \epsilon ,
\]

\[
\frac{dP}{dr}
=
-\frac{(\epsilon+P)(m+4\pi r^3P)}
{r(r-2m)} .
\]

The stellar surface is defined by

\[
P(R) \simeq 0,
\]

with total gravitational mass

\[
M = m(R).
\]

The code simultaneously integrates the tidal perturbation equations required to determine the quadrupolar Love number \(k_2\).

The dimensionless tidal deformability is

\[
\Lambda
=
\frac{2}{3}k_2 C^{-5},
\]

where

\[
C = \frac{GM}{Rc^2}
\]

is the compactness.

The geometric tidal deformability used in part of the analysis is

\[
\lambda = \Lambda M^5.
\]

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