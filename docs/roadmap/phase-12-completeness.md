# Phase 12 — Scientific completeness

**Goal:** close every item in "Known subsets" (`index.md`) so no reference
behaviour of Crystal-Builder is a simplified port, then go beyond it where
users of a standard library expect it. Items that change frozen structs
(for example a refinable `eta` in `PawleyOptions`) are designed here with a
1.x-compatible extension where possible and otherwise scheduled for 15.1.

**Release gate:** `v1.5.0`; individual tasks may ship in PATCH/MINOR
releases.

---

## Task 12.1 — Diffraction

- [ ] 12.1.1 Instrument model: wavelength, Caglioti `U, V, W`, pseudo-Voigt
  and TCH profiles, FCJ asymmetry, monochromator polarisation, any radiation
  (X-ray, neutron, electron) with neutron scattering lengths.
- [ ] 12.1.2 Anisotropic displacement parameters (`U_ij`) in PXRD and
  Rietveld.
- [ ] 12.1.3 Refinable `eta` and profile parameters in Pawley/Le Bail/
  Rietveld (new option struct, old one kept).
- [ ] 12.1.4 Single-crystal structure factors and HKL export; electron
  diffraction patterns.
- [ ] 12.1.5 Remaining 7.x items: whole-range peak refinement, persistent
  parameter sets, Pareto sweep and knee, instrument presets, indexing engine
  agreement, cancellation, powder demo.

## Task 12.2 — Optimizers and force fields

- [ ] 12.2.1 `CellFreedom` (held lengths/volume) and external pressure.
- [ ] 12.2.2 Cancellation, progress and trajectory output in every optimizer
  and scan; resumable scans written to disk.
- [ ] 12.2.3 Generic scan coordinates (cell parameters, collective variables).
- [ ] 12.2.4 Vinet and Birch-Murnaghan equations of state.
- [ ] 12.2.5 Molecular dynamics (NVE/NVT Langevin/Nosé-Hoover, NPT) with
  deterministic seeds.
- [ ] 12.2.6 Real DFTB+/xTB runs in the local CI when the binaries exist.

## Task 12.3 — Structure operations

- [ ] 12.3.1 Slab cuts that carry severed bonds and report `cut`.
- [ ] 12.3.2 Explicit bonds propagated to every symmetry image; pair-specific
  bond ranges and the metal-metal switch exposed.
- [ ] 12.3.3 Missing command types (`TransformCell`, `MakeSlab`,
  `AssignWyckoff`, `Invert`, `DescendToSubgroup`, `MergeDuplicates`,
  `InsertMolecules`, `ApplyOptimizedGeometry`, `AddHydrogens`, `Prepare`,
  `SetAtomTypes`) so every edit is undoable.
- [ ] 12.3.4 Disorder: occupancy-aware supercell ordering (enumeration of
  orderings, SQS generation with fixed seed).
- [ ] 12.3.5 Interfaces and heterostructures: lattice matching (Zur-McGill),
  strain report, stacking.
- [ ] 12.3.6 Defects: vacancies, substitutions and interstitials with
  symmetry-unique enumeration.

## Task 12.4 — Porosity and topology

- [ ] 12.4.1 Channel percolation, probe-occupiable volume, per-mass areas,
  pore-size distribution.
- [ ] 12.4.2 Structure-to-net contraction, vertex symbols, search by space
  group and transitivity, Blatov interpenetration classes.
- [ ] 12.4.3 Full RCSR and PORMAKE data bundles (`mof-data` feature) and the
  global RMSD orientation fit.
- [ ] 12.4.4 Mesh normals and OBJ/glTF output.

## Task 12.5 — Reports and modules

- [ ] 12.5.1 HTML and JSON report renderers; module registry CLI command.

## Task 12.6 — Surfaces and docs

- [ ] 12.6.1 Each closed item removed from "Known subsets", algorithm pages
  updated with references and tolerances, golden data extended.
