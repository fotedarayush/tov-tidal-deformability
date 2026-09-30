# Hyperonic EoS Ensemble

This directory contains a reproducibly selected subset of 100
zero-temperature hyperonic equations of state for neutron-star
calculations.


## Source

The models were extracted from the public zero-temperature hyperonic
EoS ensemble associated with:

**Deepak Kumar, Tuhin Malik, Hiranmaya Mishra, and Constança Providência**

Source dataset:

'hyperonic_eos_at_zero_temperature.csv.gz'

Zenodo record:

**17053872**

The original dataset is not stored in this repository.

Selection:
100 distinct EoSs chosen reproducibly using
NumPy random seed 20260930.

Four-column table format:

epsilon [MeV/fm^3]
P       [MeV/fm^3]
nB      [fm^-3]
muB     [MeV]

Source-column mapping:

epsilon = e
P       = p
nB      = rho
muB     = mu

IMPORTANT:

These tables contain the hyperonic CORE EoS from the
source dataset. They do not necessarily extend to the
very low densities required to reach a neutron-star
surface.

A suitable crust EoS should therefore be attached
before using these as complete production TOV/tidal
tables.

No quark phase is added by this conversion.
