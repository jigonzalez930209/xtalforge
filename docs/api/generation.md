# Structure generation and search

Phase 20 adds the `xtalforge-generation` crate (facade feature `generation`).
It proposes new structures instead of analysing known ones: random symmetric
crystals, molecular and low-dimensional crystals, global searches, prototype
construction and identification, and hypothetical MOF enumeration. Every
random decision takes an explicit `seed`, and the same options plus seed give
the same structures.

- Rust: `xtalforge::generation`
- CLI: `xtalforge generate …`
- Python: `generate_crystal`, `random_search`, `build_prototype`, …
- WASM: `generate_crystal_info`, `random_search_info`, `prototype_library_info`, …

## `generate_crystal`

```rust
use xtalforge::generation::{CrystalOptions, generate_crystal};

let options = CrystalOptions {
    space_group: 225,
    composition: vec![("Na".into(), 4), ("Cl".into(), 4)],
    volume_per_atom: 22.0,
    min_distance: 1.2,
    max_attempts: 64,
    seed: 7,
};
let crystal = generate_crystal(&options)?;
assert_eq!(crystal.space_group, 225);
```

`generate_crystal` (CLI `generate crystal`) picks a combination of Wyckoff
positions whose multiplicities add up to the requested counts, evaluates the
coordinate templates at random free parameters, and builds one representative
site per orbit. A candidate is accepted only when the expanded cell keeps
`min_distance` and a symmetry search returns the requested ITA number.
`generate_crystal_range` (CLI `generate crystal-range`) repeats the operation
over a list of groups and reports failures in place.

`layer_groups` returns the 80 layer groups and their Wyckoff tables.
`generate_layer_crystal` (CLI `generate layer`) uses a slab-shaped cell with a
20 Å vacuum gap; `generate_rod_crystal` (CLI `generate rod`) uses the
fourteen-group rod catalogue (identity, inversion, two-, three-, four- and
six-fold axes and their mirrors) in a wide cell.

## Molecular crystals

`generate_molecule_crystal` (CLI `generate molecule`) treats any P1 structure
as a rigid molecule, places Z' copies on Wyckoff positions with a random
orientation and checks site-symmetry compatibility: every operation that maps
a molecule centre onto itself must map the molecule onto itself within
`tolerance`, and symmetry-related molecules must keep `min_distance`.
Special positions that would distort the molecule are rejected.

## Global search

`random_search` (CLI `generate search`) is an AIRSS-style loop: generate,
relax with UFF4MOF, deduplicate with the Phase 16 matcher, rank by energy per
atom. `random_search_with` accepts any `Calculator`; the Python `random_search`
and WASM `random_search_info` use the native UFF calculator.

`evolutionary_search` (CLI `generate evolve`) keeps a population, ranks it by
energy per atom, and builds children with cut-and-splice crossover in the P1
cell, a random displacement mutation and a lattice strain. Children are
relaxed, deduplicated against the population and the best half survives. The
search stops after `convergence_generations` without a better best energy.

`search_workflow` (CLI `generate workflow`) emits a valid Phase 16 workflow
whose `compute` step calls the agent `search` verb. `SearchCheckpoint` plus
`run_checkpointed` resume a search from a JSON checkpoint; the checkpoint
stores the options, the completed seed or generation count and the hits.

## Prototypes

`prototype_library` returns the curated Strukturbericht/AFLOW-label subset in
`crates/xtalforge-generation/data/prototypes.json` (A1, A2, A3, A4, A15, B1,
B2, B3, B4, B8, C1, C2, C14, C16, D03, E1, L10, L12, L21) with their lattice
ratios, Wyckoff letters and ideal free parameters. `build_prototype(label,
scale)` builds one; CLI `generate prototypes` lists them and `generate
prototype B1` builds one.

`identify_prototype(structure)` detects the space group, the sorted
`(multiplicity, letter)` signature and the stoichiometry, filters the library
and confirms the remaining candidates with the structure matcher. The result
reports the label, the common name, the signature and the RMS.

`substitution_candidates` (CLI `generate substitute`) estimates oxidation
states with the bond-valence machinery and proposes ions from the embedded
Shannon radius table with the same charge and a radius within
`radius_tolerance`.

## Hypothetical MOFs

`generate_mof_dataset` (CLI `generate mof-dataset`) enumerates topology ×
node × linker combinations from the PORMAKE bundle (curated subset by default,
complete bundle with the `mof-data` feature), builds each framework, optionally
relaxes it with UFF4MOF, records framework and pore descriptors and marks
MOFkey duplicates. A combination that cannot be assembled is recorded with its
error instead of aborting the dataset.

`functionalization_sweep` (CLI `generate functionalize`) applies amino,
hydroxyl, methyl and fluorine groups over a list of fractions and records the
descriptors of each variant. `global_orientation_fit` exposes the Kabsch
orientation used by the builder: every one-to-one assignment of connection
points to slot directions is tried and the smallest RMSD wins.

## Data

| File | Content | Provenance |
|---|---|---|
| `crates/xtalforge-generation/data/prototypes.json` | 19 prototype entries | Published crystallographic facts and Strukturbericht/AFLOW labels; no AFLOW record copied |
| `crates/xtalforge-generation/data/ionic_radii.json` | 48 Shannon radii (CN VI) | Shannon (1976), factual values |

## Errors

The crate returns `xtalforge_core::Error`; the CLI maps it to exit code 3 and
Python raises `XtalforgeError`. Invalid seeds, space-group numbers, counts and
tolerances are rejected before any work starts. A search that exhausts its
attempts reports the last rejection reason.

## Limits

- The generator caps a conventional cell at 400 atoms and a molecular
  asymmetric unit at 400 atoms per copy.
- The rod catalogue is a documented subset of the 75 rod groups; it is not a
  complete table.
- The prototype library is a curated subset, not the full AFLOW encyclopedia.
- The evolutionary crossover is P1 cut-and-splice; it does not enforce a space
  group on the child.
- WASM runs searches synchronously inside the call; callers that need
  cancellation should use the worker wrapper from Phase 12.
- Checkpoints need a filesystem and are therefore a Rust/CLI/Python feature.

Parity keys: `generate_crystal`, `generate_molecule_crystal`, `layer_groups`,
`generate_layer_crystal`, `generate_rod_crystal`, `random_search`,
`evolutionary_search`, `search_workflow`, `prototype_library`,
`identify_prototype`, `substitution_candidates`, `mof_dataset`,
`functionalization_sweep`.
