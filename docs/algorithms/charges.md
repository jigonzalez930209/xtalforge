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

`eeq_charges` is EEQ2019 (Caldeweyher et al., J. Chem. Phys. 150, 154122,
2019), ported from sci-form `src/charges_eeq/eeq2019.rs` with the
grimme-lab/multicharge parameters for H to Lr (`crates/xtalforge-ff/data/
eeq2019.json`). In atomic units it solves

```
A_ii = eta_i + sqrt(2/pi) / rad_i
A_ij = sum_images erf(r / sqrt(rad_i^2 + rad_j^2)) / (eps r)
x_i  = -chi_i + kcn_i sqrt(CN_i),   A q = x,   sum q = Q
```

with an erf-counted coordination number (D3 covalent radii, smooth cap at 8,
images within 25 Bohr). Coulomb images are summed inside `options.cutoff`, the
same real-space rule as QEq; `dielectric = 1` is plain EEQ2019. For isolated
water and methanol the charges match sci-form (validated against dftd4) to
2e-3 (`eeq2019_matches_the_sci_form_reference`). Before the 2026-11-11 audit
this function reused the QEq table and returned exactly the QEq charges.

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
- EEQ: EEQ2019 with periodic images in real space (no Ewald sum).
- Ewald: the periodic electrostatic energy and forces on top of any charge
  model.
