# Phase 19 — Magnetism, disorder and spectroscopy

**Goal:** open the library to the magnetism, alloy and experimental
communities: magnetic symmetry, partial occupancies as first-class data,
derivative-structure enumeration, and simulated spectra that can be compared
directly with measurements.

**Release gate:** `v1.12.0`.

---

## Task 19.1 — Magnetic structures

- [ ] 19.1.1 Magnetic moments on sites
  Description: collinear and non-collinear moments as a site property,
  preserved by every transform (supercell, slab, standardisation) and by the
  serializer; schema change is additive until Phase 22.
- [ ] 19.1.2 Magnetic space groups
  Description: detection of the magnetic space group (BNS and OG settings)
  through moyo's magnetic cell support, with symmetrisation of moments.
  Acceptance: golden data from spglib's magnetic tests.
- [ ] 19.1.3 mCIF read/write
  Description: magCIF reader and writer (`_space_group_symop_magn_*`,
  `_atom_site_moment_*`), round-trip tested on the Bilbao MAGNDATA examples
  that are redistributable.
- [ ] 19.1.4 Magnetic orderings
  Description: enumeration of ferro, antiferro and ferri orderings for a
  given supercell size, ranked by symmetry; output ready for the input sets
  of Task 18.1 (MAGMOM).

## Task 19.2 — Disorder and derivative structures

- [ ] 19.2.1 Partial occupancies as first-class data
  Description: occupancy and site-mixing groups kept through every operation;
  disorder-aware composition, density and PXRD.
- [ ] 19.2.2 Derivative-structure enumeration
  Description: Hermite-normal-form supercells plus symmetry-distinct
  decorations (enumlib algorithm), with concentration limits and a seed-free
  deterministic order. Subsumes the struck subtask 14.3.4.
- [ ] 19.2.3 SQS beyond binary
  Description: multicomponent SQS with pair, triplet and quadruplet
  correlations (extends known subset 9.5).
- [ ] 19.2.4 Cluster-expansion features
  Description: cluster correlation vectors for a set of decorations, exported
  as arrays for external fitting; no fitting engine is shipped.

## Task 19.3 — Diffraction beyond X-rays

- [ ] 19.3.1 Neutron powder diffraction
  Description: nuclear scattering lengths (data file with source and
  license), pattern simulation on the Phase 7 profile model, time-of-flight
  axis.
- [ ] 19.3.2 Magnetic neutron diffraction
  Description: magnetic form factors and magnetic structure factors from the
  moments of Task 19.1.
- [ ] 19.3.3 Electron diffraction
  Description: kinematic SAED patterns for a zone axis with electron
  scattering factors; image buffer output for the playground.

## Task 19.4 — Vibrational spectroscopy

- [ ] 19.4.1 IR intensities
  Description: IR spectra from the phonon modes (Phase 4) with Born effective
  charges from QEq/EQeq or read from an engine output (Task 18.2); Lorentzian
  broadening.
- [ ] 19.4.2 Raman activities from engine outputs
  Description: Raman tensors read from VASP/QE outputs and the powder-averaged
  spectrum; no native polarizability model.
- [ ] 19.4.3 Phonon DOS and inelastic neutron weights from the existing
  phonon module.

## Task 19.5 — Surfaces and docs

- [ ] 19.5.1 Parity matrix; API pages `magnetism.md`, `disorder.md`,
  `spectroscopy.md`; algorithm pages for magnetic symmetry, enumeration and
  neutron/electron diffraction; playground views for patterns and spectra.
