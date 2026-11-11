# Phase 14 — Scientific completeness

**Goal:** close every item in "Known subsets" (`index.md`) so no reference
behaviour of Crystal-Builder is a simplified port, then go beyond it where
users of a standard library expect it. Items that change frozen structs
(for example a refinable `eta` in `PawleyOptions`) are designed here with a
1.x-compatible extension where possible and otherwise scheduled for 17.1.

**Release gate:** `v1.7.0`; individual tasks may ship in PATCH/MINOR
releases.

---

## Task 14.1 — Diffraction

- [ ] 14.1.1 Instrument model: wavelength, Caglioti `U, V, W`, pseudo-Voigt
  and TCH profiles, FCJ asymmetry, monochromator polarisation, any radiation
  (X-ray, neutron, electron) with neutron scattering lengths.
- [ ] 14.1.2 Anisotropic displacement parameters (`U_ij`) in PXRD and
  Rietveld.
- [ ] 14.1.3 Refinable `eta` and profile parameters in Pawley/Le Bail/
  Rietveld (new option struct, old one kept).
- [ ] 14.1.4 Single-crystal structure factors and HKL export; electron
  diffraction patterns.
- [ ] 14.1.5 Remaining 7.x items: whole-range peak refinement, persistent
  parameter sets, Pareto sweep and knee, instrument presets, indexing engine
  agreement, cancellation, powder demo.

## Task 14.2 — Optimizers and force fields

- [ ] 14.2.1 `CellFreedom` (held lengths/volume) and external pressure.
- [ ] 14.2.2 Cancellation, progress and trajectory output in every optimizer
  and scan; resumable scans written to disk.
- [ ] 14.2.3 Generic scan coordinates (cell parameters, collective variables).
- [ ] ~~14.2.4 Vinet and Birch-Murnaghan equations of state.~~
  Note: moved to core Task 9.6.2 (2026-11-11).
- [ ] 14.2.5 Molecular dynamics (NVE/NVT Langevin/Nosé-Hoover, NPT) with
  deterministic seeds.
- [ ] 14.2.6 Real DFTB+/xTB runs in the local CI when the binaries exist.

## Task 14.3 — Structure operations

- [ ] 14.3.1 Slab cuts that carry severed bonds and report `cut`.
- [ ] 14.3.2 Explicit bonds propagated to every symmetry image; pair-specific
  bond ranges and the metal-metal switch exposed.
- [ ] 14.3.3 Missing command types (`TransformCell`, `MakeSlab`,
  `AssignWyckoff`, `Invert`, `DescendToSubgroup`, `MergeDuplicates`,
  `InsertMolecules`, `ApplyOptimizedGeometry`, `AddHydrogens`, `Prepare`,
  `SetAtomTypes`) so every edit is undoable.
- [ ] ~~14.3.4 Disorder: occupancy-aware supercell ordering~~
  Note: moved to core Tasks 9.5.2 and 9.5.3 (2026-11-11).
- [ ] ~~14.3.5 Interfaces and heterostructures~~
  Note: moved to core Task 9.5.4 (2026-11-11).
- [ ] ~~14.3.6 Defects: vacancies, substitutions and interstitials~~
  Note: moved to core Task 9.5.1 (2026-11-11).

## Task 14.4 — Porosity and topology

- [ ] ~~14.4.1 Channel percolation, probe-occupiable volume, per-mass areas,
  pore-size distribution.~~
  Note: moved to MOF Tasks 10.3.1 and 10.3.2 (2026-11-11).
- [ ] 14.4.2 Structure-to-net contraction, vertex symbols, search by space
  group and transitivity, Blatov interpenetration classes.
- [ ] ~~14.4.3 Full RCSR and PORMAKE data bundles (`mof-data` feature) and the
  global RMSD orientation fit.~~
  Note: moved to MOF Task 10.5.1 (2026-11-11).
- [ ] 14.4.4 Mesh normals and OBJ/glTF output.

## Task 14.5 — Reports and modules

- [ ] 14.5.1 HTML and JSON report renderers; module registry CLI command.

## Task 14.6 — Surfaces and docs

- [ ] 14.6.1 Each closed item removed from "Known subsets", algorithm pages
  updated with references and tolerances, golden data extended.
