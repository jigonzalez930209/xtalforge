# Charges and electrostatics

## QEq

`qeq_charges(structure, &ChargeOptions)` solves the Rappé-Goddard charge
equilibration. Each element contributes an electronegativity `chi` and a
hardness `J`; the energy is `sum_i (chi_i Q_i + J_i Q_i^2 / 2) + sum_{i<j}
J_ij Q_i Q_j` with the Ohno-Klopman kernel
`J_ij = 14.4 / sqrt(r^2 + (14.4 / sqrt(J_i J_j))^2)` (eV, Angstrom). It uses the
hardnesses only, not the UFF bond radii. The total charge is enforced by a
Lagrange multiplier (one KKT system).

The linear system is solved on the P1 cell with a real-space cutoff (default
12 Angstrom, summed over every image inside it, with no Ewald correction). The
truncated sum under-polarises ionic solids: NaCl comes out near +-0.20, not
+-1, so use QEq for molecules and frameworks, and EQeq or an explicit charge
set for ionic crystals; the result
is `Charges { per_atom, total }`. `apply_to` writes the per-site mean over the
images back to the asymmetric unit.

## EQeq

`eqeq_charges` keeps the same solver and replaces the parameters with measured
ionisation energies and electron affinities (Wilmer, Kim and Snurr 2012),
expanded about each element's usual charge centre. A dielectric of 1.2 is the
reference default. This is the variant that keeps framework metals from
running away from the neutral-atom parabola.

## EEQ

`eeq_charges` is the EEQ kernel `C / (eps sqrt(r^2 + R_i R_j))` with a Coulomb
radius per element and a dielectric `eps`. The fitted EEQ2019 table of sci-form
is **not** vendored: `chi` and `J` come from the QEq table and each radius is
the Ohno radius implied by the hardness. With the default `dielectric = 1` this
is the same linear system as QEq, and the two methods return identical charges
(NaCl +-0.202, MOF-5 from -0.617 to +1.246). Only a dielectric above one
(`ChargeOptions::eqeq`) gives screened charges that differ from QEq, so EEQ is
not an independent third method yet; the corpus comparison against the fitted
table is not done.

## Ewald

![Ewald splitting](/diagrams/ewald-split.svg)


`ewald_energy_forces(structure, &charges, &EwaldOptions)` splits the Coulomb
sum into a real-space part, a reciprocal-space part, a self term and a
neutralising background. The splitting parameter and the cutoffs come from the
requested relative accuracy; the summation order is fixed, so the energy is
deterministic. Forces are the analytic derivative; the test compares them with
finite differences and the energy with a spherical direct sum.

## When to use which

- QEq: fast, neutral organic molecules and a first guess.
- EQeq: frameworks with metals, where measured ionisation energies matter.
- EEQ: QEq with a dielectric screening; identical to QEq at `dielectric = 1`.
- Ewald: the periodic electrostatic energy and forces on top of any charge
  model.
