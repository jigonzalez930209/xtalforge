# Phase 20 — Structure generation and search

**Goal:** move from analysing known structures to proposing new ones:
symmetry-constrained random crystals, global search with the native or ML
calculators, a prototype library, and large-scale hypothetical MOFs. Every
generator takes an explicit seed (AGENTS.md rule 6).

**Release gate:** `v1.13.0`.

---

## Task 20.1 — Random symmetric crystals

- [x] 20.1.1 Wyckoff-constrained generator
  Description: given a composition, a space group (or a range) and a volume
  factor, choose compatible Wyckoff combinations and place atoms with
  minimum-distance constraints (PyXtal-like). Uses the Wyckoff tables of
  Phase 2.
  Acceptance: every generated structure re-detects its requested space group;
  statistics over 10 000 seeds documented in the algorithm page.
  Result: the generator, its tests and `generation_statistics` are in
  `xtalforge-generation`; the algorithm page documents the seed sweep and
  states that CI covers a small range rather than 10 000 seeds.
- [x] 20.1.2 Molecular crystals
  Description: rigid molecules (from SMILES, Task 17.1.1, or a file) placed
  on general and special positions with site-symmetry compatibility checks.
  Result: `generate_molecular_crystal` with the geometric compatibility test.
- [x] 20.1.3 Low-dimensional generators: 2D layer groups and 1D rod groups.

## Task 20.2 — Global search

- [x] 20.2.1 AIRSS-style random search
  Description: generate, relax (any `Calculator`: UFF, sci-form xTB, ML
  bridge or an engine of Task 18.4), deduplicate (Task 16.3.2) and rank.
- [x] 20.2.2 Evolutionary search
  Description: cut-and-splice crossover, lattice and permutation mutations,
  fingerprint-based diversity, convex-hull fitness (Task 9.6) for variable
  composition.
- [x] 20.2.3 Search as a workflow
  Description: both searches run as workflow steps with checkpoints, so a
  search can resume and record full provenance.
  Result: the agent `search` verb, `SearchCheckpoint`, `run_checkpointed` and
  `search_workflow_document`; the checkpoint stores full hits and options.

## Task 20.3 — Prototypes

- [x] 20.3.1 Prototype library
  Description: the AFLOW prototype encyclopedia labels (data file with source
  and license) with their parameters; build any prototype from its label and
  parameters.
- [x] 20.3.2 Prototype identification
  Description: assign the AFLOW label of an input structure (space group,
  Wyckoff letters, stoichiometry), then confirm with the matcher of Task 9.1.
- [x] 20.3.3 Substitution-based generation
  Description: new compositions from a prototype by chemically similar
  substitutions (ionic-radius and oxidation-state rules from Task 9.4 BVS).

## Task 20.4 — Hypothetical MOFs at scale

- [x] 20.4.1 Global orientation fit for the builder
  Description: replace the greedy orientation of the MOF builder with a
  global RMSD fit (known subset 6.x), using the complete PORMAKE bundle of
  Task 17.1.3.
  Result: the builder already used `align::kabsch_orientation` from Phase 6;
  `global_orientation_fit` exposes it and the test covers a known rotation.
- [x] 20.4.2 Combinatorial MOF generation
  Description: enumerate topology × node × linker combinations with
  compatibility filters, build, check (Phase 10 checker), relax with
  UFF4MOF and compute descriptors, as a dataset (Task 16.3).
  Result: `generate_mof_dataset` records build errors in place and marks
  MOFkey duplicates; the complete bundle stays behind `mof-data`.
- [x] 20.4.3 Functionalisation sweeps
  Description: systematic linker functionalisation with the PSM tools of
  Phase 10 for property screening.
  Result: `functionalization_sweep` over groups × fractions with descriptors.

## Task 20.5 — Surfaces and docs

- [x] 20.5.1 Parity matrix (long searches in WASM run in workers, Task 12.3);
  API page `generation.md`; algorithm pages for the generator and both
  searches; tutorial "hypothetical MOF screening in one workflow".
