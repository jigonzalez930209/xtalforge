# MOF adsorption: force fields, Widom and GCMC

**Tasks 10.4.1-10.4.5.** References: TraPPE (Martin and Siepmann 1998;
Potoff and Siepmann 2001; Chen and Siepmann 2000), RASPA, UFF (Rappe et al.
1992), EQeq (Wilmer et al. 2012).

## Host-guest force field

Framework atoms carry UFF Lennard-Jones parameters (`sigma = x1`,
`epsilon = d1`) and EQeq charges from Task 4.3. Guests are the embedded
TraPPE models: CO2 (three sites, charges +0.70/-0.35), N2 (two sites,
+/-0.482), CH4 (one site), H2 (two sites) and SPC/E water. Mixing is
Lorentz-Berthelot:

```
sigma_ij = (sigma_i + sigma_j) / 2
epsilon_ij = sqrt(epsilon_i epsilon_j)
```

The pair energy is `4 epsilon [(sigma/r)^12 - (sigma/r)^6]` plus the Coulomb
term `14.399645 q_i q_j / r` in eV Angstrom. The host-guest energy sums every
guest atom over every host atom with the periodic minimum image and a 15
Angstrom cutoff; guest-guest interactions are not part of the infinite-
dilution Widom calculation.

## Widom insertion and Henry coefficient

A seeded SplitMix64 generator draws a uniform fractional position and three
Euler angles per insertion. The dimensionless Henry coefficient is the
average Boltzmann factor

```
H = < exp(-E / kT) >
```

over the samples. The mean energy is the Boltzmann-weighted average, the heat
of adsorption is its negative, and the minimum energy over the samples is
reported as the best binding site. The result is deterministic for a seed and
independent of thread count because insertion order is fixed.

## GCMC isotherm

The single-component engine uses the ideal-gas chemical potential
`mu = kT ln(P / P0)` with `P0 = 1 bar` and accepts insertions with
`min(1, exp((mu - E)/kT))`; deletions are attempted with a small fixed
probability. The loading is reported per cell and the acceptance ratio is the
fraction of accepted insertion attempts. The implementation is deliberately a
Henry/Langmuir-level isotherm: it omits guest-guest interactions and
Peng-Robinson fugacity, so the full RASPA engine remains the reference for
production isotherms.

## Energy grids and RASPA export

The energy grid evaluates the host-guest energy of a fixed guest orientation
at every sample of a regular cell grid and reports the minimum and maximum.
It shares the grid conventions of the porosity field, so a viewer can draw
the binding sites next to the accessible surface.

`raspa_input` writes `simulation.input`, the framework CIF, a pseudo-atom
force-field definition and the Lorentz-Berthelot mixing file. The input sets
`WidomInsertion yes`; a user can switch it to the GCMC block without changing
the framework or force-field files.
