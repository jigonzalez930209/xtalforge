# Phase 20 — Structure generation and search

**Goal:** move from analysing known structures to proposing new ones:
symmetry-constrained random crystals, global search with the native or ML
calculators, a prototype library, and large-scale hypothetical MOFs. Every
generator takes an explicit seed (AGENTS.md rule 6).

**Release gate:** `v1.13.0`.

---

## Task 20.1 — Random symmetric crystals

- [ ] 20.1.1 Wyckoff-constrained generator
  Description: given a composition, a space group (or a range) and a volume
  factor, choose compatible Wyckoff combinations and place atoms with
  minimum-distance constraints (PyXtal-like). Uses the Wyckoff tables of
  Phase 2.
  Acceptance: every generated structure re-detects its requested space group;
  statistics over 10 000 seeds documented in the algorithm page.
- [ ] 20.1.2 Molecular crystals
  Description: rigid molecules (from SMILES, Task 17.1.1, or a file) placed
  on general and special positions with site-symmetry compatibility checks.
- [ ] 20.1.3 Low-dimensional generators: 2D layer groups and 1D rod groups.

## Task 20.2 — Global search

- [ ] 20.2.1 AIRSS-style random search
  Description: generate, relax (any `Calculator`: UFF, sci-form xTB, ML
  bridge or an engine of Task 18.4), deduplicate (Task 16.3.2) and rank.
- [ ] 20.2.2 Evolutionary search
  Description: cut-and-splice crossover, lattice and permutation mutations,
  fingerprint-based diversity, convex-hull fitness (Task 9.6) for variable
  composition.
- [ ] 20.2.3 Search as a workflow
  Description: both searches run as workflow steps with checkpoints, so a
  search can resume and record full provenance.

## Task 20.3 — Prototypes

- [ ] 20.3.1 Prototype library
  Description: the AFLOW prototype encyclopedia labels (data file with source
  and license) with their parameters; build any prototype from its label and
  parameters.
- [ ] 20.3.2 Prototype identification
  Description: assign the AFLOW label of an input structure (space group,
  Wyckoff letters, stoichiometry), then confirm with the matcher of Task 9.1.
- [ ] 20.3.3 Substitution-based generation
  Description: new compositions from a prototype by chemically similar
  substitutions (ionic-radius and oxidation-state rules from Task 9.4 BVS).

## Task 20.4 — Hypothetical MOFs at scale

- [ ] 20.4.1 Global orientation fit for the builder
  Description: replace the greedy orientation of the MOF builder with a
  global RMSD fit (known subset 6.x), using the complete PORMAKE bundle of
  Task 17.1.3.
- [ ] 20.4.2 Combinatorial MOF generation
  Description: enumerate topology × node × linker combinations with
  compatibility filters, build, check (Phase 10 checker), relax with
  UFF4MOF and compute descriptors, as a dataset (Task 16.3).
- [ ] 20.4.3 Functionalisation sweeps
  Description: systematic linker functionalisation with the PSM tools of
  Phase 10 for property screening.

## Task 20.5 — Surfaces and docs

- [ ] 20.5.1 Parity matrix (long searches in WASM run in workers, Task 12.3);
  API page `generation.md`; algorithm pages for the generator and both
  searches; tutorial "hypothetical MOF screening in one workflow".
