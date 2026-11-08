# Charges and electrostatics

## QEq

`qeq_charges(structure, &ChargeOptions)` solves the Rappé-Goddard charge
equilibration. Each element contributes an electronegativity `chi` and a
hardness `eta`; the Coulomb interaction uses the Ohno-Klopman kernel with the
UFF bond radii. The total charge is enforced by a Lagrange multiplier.

The linear system is solved on the P1 cell with a real-space cutoff; the result
is `Charges { per_atom, total }`. `apply_to` writes the per-site mean over the
images back to the asymmetric unit.

## EQeq

`eqeq_charges` keeps the same solver and replaces the parameters with measured
ionisation energies and electron affinities (Wilmer, Kim and Snurr 2012),
expanded about each element's usual charge centre. A dielectric of 1.2 is the
reference default. This is the variant that keeps framework metals from
running away from the neutral-atom parabola.

## EEQ

`eeq_charges` is the EEQ-style kernel (electronegativity, hardness and a
Coulomb radius per element, Ohno-Klopman interaction, periodic cutoff). The
exact sci-form parameterisation is not vendored; the module documents the local
table and the approximation. It is a third method, not a drop-in replacement
for EQeq, and the corpus comparison is the acceptance criterion.

## Ewald

`ewald_energy_forces(structure, &charges, &EwaldOptions)` splits the Coulomb
sum into a real-space part, a reciprocal-space part, a self term and a
neutralising background. The splitting parameter and the cutoffs come from the
requested relative accuracy; the summation order is fixed, so the energy is
deterministic. Forces are the analytic derivative; the test compares them with
finite differences and the energy with a spherical direct sum.

## When to use which

- QEq: fast, neutral organic molecules and a first guess.
- EQeq: frameworks with metals, where measured ionisation energies matter.
- EEQ: an alternative kernel for comparison and for callers coming from
  sci-form.
- Ewald: the periodic electrostatic energy and forces on top of any charge
  model.
