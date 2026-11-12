# Phase 9 — Core features of a standard crystal library

**Goal:** add the core operations that users of pymatgen, ASE and spglib
expect from any crystal-structure library and that neither Crystal-Builder
nor xtalforge 1.1 provides. Phases 1-8 reached parity with Crystal-Builder;
this phase goes beyond it. Every subtask is additive and ships on the four
surfaces (AGENTS.md rule 4) with golden data from an independent reference
(pymatgen, ASE, spglib, phonopy) and documented tolerances (rule 9).

**Release gate:** `v1.2.0` when Tasks 9.1-9.6 are complete; tasks may ship
earlier in PATCH-free MINOR steps only if the maintainer decides so.

**Gap check (2026-11-11):** a search of `crates/*/src` found no structure
matcher, RMSD, rattle/perturb, strain application, RDF, Voronoi or
coordination environments, bond-valence sums, oxidation-state guessing,
unique Miller indices, Wulff shapes, vacancy/interstitial enumeration, SQS or
elastic tensors. These are the subtasks below.

**Reference:** pymatgen (`core.structure`, `analysis.structure_matcher`,
`analysis.local_env`, `analysis.bond_valence`, `core.surface`,
`analysis.wulff`, `analysis.defects`, `transformations`), ASE (`build`,
`geometry`, `neighborlist`), spglib, phonopy, icet/ATAT `mcsqs`. Code is
written from the published algorithms; any ported code is recorded in
`THIRD_PARTY.md` first (rule 8).

---

## Task 9.1 — Structure comparison

- [x] 9.1.1 Structure matcher
  Description: decide whether two structures are the same up to lattice
  choice, origin, permutation and tolerance (`ltol`, `stol`, `angle_tol`),
  via Niggli reduction and a lattice-mapping search (Task 2.6).
  Options: compare species exactly, by element or ignoring species; allow
  supercell matching. Output: match flag, mapping, RMS displacement.
  Acceptance: agreement with pymatgen `StructureMatcher` on the corpus and on
  generated permuted/shifted/strained copies. Surfaces: all four, CLI
  `xtalforge compare a.cif b.cif`.
- [x] 9.1.2 RMSD and alignment
  Description: periodic RMSD under the best mapping, and Kabsch alignment for
  molecules. Surfaces: all four.
- [x] 9.1.3 Group-by-equality
  Description: group a list of structures into classes of equal structures,
  deterministic order (rule 6). Surfaces: all four.

## Task 9.2 — Structure edits that every library has

- [x] 9.2.1 Strain and deformation
  Description: apply a strain tensor or deformation gradient (Voigt or 3x3),
  isotropic scaling to a target volume or density. Surfaces: all four, CLI
  `xtalforge deform`.
- [x] 9.2.2 Rattle and perturb
  Description: random displacements with a fixed seed (Gaussian or bounded
  norm), optionally per element; for training sets and symmetry breaking.
- [x] 9.2.3 Site editing on every surface
  Description: add, remove, replace and move sites; replace species by map
  (`{"Si": "Ge"}`), partial replacement by fraction with a seed; sort sites
  by element or electronegativity; wrap and fold. Undoable through commands.
- [x] 9.2.4 Interpolation
  Description: linear interpolation between two matched structures (NEB
  images) and image output as a trajectory.
- [x] 9.2.5 Make orthogonal/cubic-like supercell
  Description: find the supercell matrix closest to a cube for a target atom
  count (ASE `find_optimal_cell_shape`). Surfaces: all four.

## Task 9.3 — Local environment and chemistry

- [x] 9.3.1 Voronoi tessellation
  Description: periodic Voronoi cells (pure Rust), solid-angle weighted
  neighbors, cell volumes. Algorithm page.
- [x] 9.3.2 Coordination environments
  Description: coordination number by several strategies (cutoff, Voronoi,
  CrystalNN-like), polyhedron identification (tetrahedral, octahedral, …) with
  continuous shape measures.
- [x] 9.3.3 Bond-valence sums and oxidation states
  Description: bond-valence parameters as a data file with citation;
  oxidation-state guessing by charge balance and BVS. Golden data from
  pymatgen `BVAnalyzer`.
- [x] 9.3.4 Radial distribution and structure fingerprints
  Description: partial and total RDF/PDF with smoothing, bond-angle
  distributions, a composition-and-site fingerprint for quick similarity.
  Arrays on Python (Task 12.2) and typed arrays in WASM.

## Task 9.4 — Surfaces of crystals

- [x] 9.4.1 Symmetry-unique Miller indices
  Description: enumerate the distinct `(h k l)` up to a maximum index under the
  point group of the crystal.
- [x] 9.4.2 All terminations of a slab
  Description: generate every distinct termination for one `(h k l)`,
  symmetric slabs, polar-surface report, minimal slab for a thickness.
  Extends the slab code of Task 1.8.
- [x] 9.4.3 Wulff construction
  Description: Wulff shape from surface energies, facet areas and fractions,
  mesh output through the view crate.
- [x] 9.4.4 Adsorption sites
  Description: top, bridge and hollow sites on a slab surface and placement of
  a molecule on them (reuses the guest placement of Task 6.4).

## Task 9.5 — Defects, disorder and interfaces

- [x] 9.5.1 Point defects
  Description: symmetry-unique vacancies, substitutions and interstitials
  (Voronoi-based candidates from 9.3.1), with defect supercells of a minimum
  image distance.
- [x] 9.5.2 Ordering of partial occupancies
  Description: enumerate symmetry-distinct orderings of a disordered structure
  in a supercell (Hart-Forcade style), ranked by Ewald energy (Task 4.3).
- [x] 9.5.3 Special quasirandom structures
  Description: SQS generation with a fixed seed, pair/triplet correlation
  targets and report.
- [x] 9.5.4 Interfaces and heterostructures
  Description: lattice matching (Zur-McGill) with a strain/area report and
  stacking of two slabs with a gap and in-plane shift.

## Task 9.6 — Mechanical and thermodynamic basics

- [x] 9.6.1 Strain sets and elastic tensor
  Description: generate the deformation set, fit the elastic tensor from
  stresses (any `Calculator`, Task 4.1), derived moduli (Voigt/Reuss/Hill),
  stability criteria.
- [x] 9.6.2 Equation-of-state fits
  Description: Birch-Murnaghan, Vinet and Murnaghan fits from energy-volume
  data, with parameters and errors (subsumes the Vinet item in "Known
  subsets").
- [x] 9.6.3 Phase-diagram basics
  Description: formation energy per atom and convex hull distance from a list
  of entries (energies supplied by the user or a calculator).

## Task 9.7 — Surfaces and docs

- [x] 9.7.1 Parity, API pages and algorithm pages
  Description: every subtask above in `docs/api/parity.md`, an API page per
  task and an algorithm page with references for 9.1, 9.3, 9.4, 9.5 and 9.6;
  the playground gains compare, deform, Voronoi and Wulff views; e2e checks
  against pymatgen/ASE references extend `tests/e2e/reference/generate.py`.
  Note: the parity matrix, the `docs/api/features.md` page and the five
  algorithm pages are complete for `v1.2.0`. The playground panels and the
  Playwright reference extensions are tracked as Task 9.7.2 so the release
  gate (Tasks 9.1-9.6) is not blocked by user-interface work.
- [ ] 9.7.2 Playground views and e2e references for Phase 9
  Description: add compare, deform, Voronoi, Wulff, RDF and elastic panels to
  `docs/.vitepress/components/Playground.vue`, and extend
  `tests/e2e/reference/generate.py` with pymatgen/ASE references for the
  matcher, Voronoi volumes, coordination numbers, RDFs and elastic tensors.
