# TOV and Tidal Deformability Code

Fortran code for constructing neutron-star sequences by solving the
Tolman-Oppenheimer-Volkoff equations together with the quadrupolar
tidal perturbation equations.

The code calculates quantities including:

- gravitational mass
- radius
- baryonic/rest mass
- compactness
- tidal Love number k2
- tidal coupling parameter
- dimensionless tidal deformability Lambda

## Physics

For each central energy density, the code determines the corresponding
central pressure from an equation of state and integrates the stellar
structure equations outward until the stellar surface is reached.

The dimensionless tidal deformability is

\[
\Lambda = \frac{2}{3} k_2 C^{-5},
\]

where

\[
C = \frac{GM}{Rc^2}.
\]

## Compilation

Compile using GNU Fortran:

```bash
gfortran -o logtov_seq_geom_tidal.out \
    logtov_seq_geom_tidal.f90 -fno-automatic
