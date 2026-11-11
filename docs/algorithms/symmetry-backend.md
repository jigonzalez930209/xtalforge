# Symmetry backend: moyo audit

**Verdict: adopt moyo 0.22 with gaps.** It covers space-group detection,
operations, orbits, Wyckoff letters, standardized and primitive cells,
subgroups, the Euclidean normalizer, layer groups and magnetic groups. Three
gaps remain: the Niggli/Delaunay/Minkowski reductions are crate-private, the
symmetrization module is private (reachable only through `MoyoDataset`), and
Brillouin zones, k-paths and molecular point groups are absent.

`xtalforge-symmetry` exposes its own types (`SpaceGroup`, `SymOp`,
`SymmetryDataset`, `WyckoffSite`); `moyo` stays an implementation detail and
can be replaced without changing the public API.

## Requirements from Crystal-Builder

Source: `Crystal-Builder/xtal/core/{spacegroup,symmetry,subgroups}.py`,
`core/molecule.py` and `analysis/kpath.py` (Crystal-Builder commit
`78d5476fd9e851fda00cde57b7ab8da3195dcf54`).

| Crystal-Builder need | Function or class | moyo 0.22 coverage | Status |
|---|---|---|---|
| Space group lookup by Hall symbol | `SpaceGroup.from_hall` | `data::hall_symbol_entry`, `identify` | public |
| Lookup by Hermann-Mauguin and setting | `SpaceGroup.from_name`, `from_number` | `data::*`, `Setting` | public |
| Operations and triplets | `SpaceGroup.operations`, `SymOp.apply` | `MoyoDataset::operations` | public |
| Inverses and stacked operations | `SpaceGroup.inverse_of`, `stacked` | operations field + own arithmetic | port |
| Detection with tolerance | `symmetry.detect` | `MoyoDataset::new` | public |
| Wyckoff letters and site symmetry | `assign_wyckoff` | `wyckoffs`, `site_symmetry_symbols` | public |
| P1 expansion and standardization | `reduce_to_p1`, `standardize` | `std_cell`, `prim_std_cell`, mappings | public |
| Set space group / reinterpret / impose | `set_space_group`, `_reinterpret`, `_impose` | `MoyoDataset` + own fit | port |
| Asymmetrize and snap to special positions | `asymmetrize`, `snap_to_special_positions` | `MoyoDataset` + `std_cell` | port |
| Enantiomorph and hand description | `enantiomorph`, `hand_description` | `identify`, operations | port |
| Duplicate and merge sites | `duplicate_groups`, `merge_duplicates` | own soak/merge logic | port |
| Translationengleiche subgroups | `subgroups.subgroups_of` | `subgroup::enumerate_translationengleiche_subgroups` | public |
| Klassengleiche subgroups | `subgroups.maximal_subgroups` | `subgroup::enumerate_klassengleiche_subgroups` | public |
| Euclidean normalizer | `normalizer` paths | `MoyoDataset::euclidean_normalizer`, `normalizer_wyckoff_positions` | public |
| Niggli, Delaunay, Minkowski reduction | `supercell.niggli_reduce`, `delaunay_reduce` | `math::{niggli,delaunay,minkowski}` | **crate-private** |
| Symmetrization internals | `symmetry.standardize` | `symmetrize` module | **private** |
| Brillouin zone and k-path | `analysis/kpath.py` | absent | **implement** |
| Molecular point group | `core/molecule.py` | `identify::PointGroup` takes rotations | **implement adapter** |

## Behavioural comparison against spglib

Input: `tests/golden/moyo-audit.json`, generated from eight corpus structures
(Si, NaCl, graphite, quartz, MOF-5, ZIF-8, and the triclinic and hexagonal
edge cases) for `symprec` in {1e-5, 1e-3, 1e-2, 0.1}. The Rust test
`crates/xtalforge-symmetry/tests/moyo_audit.rs` runs moyo on the same cells and
compares:

- space-group number: **agrees in all 32 cases**;
- Hall number with `Setting::Spglib`: **agrees in all 32 cases**; with
  `Setting::default`, Si selects Hall 526 instead of spglib's 525 (origin
  choice 2 instead of 1);
- number of operations: **agrees in all 32 cases**;
- Wyckoff letters: raw moyo agrees except Si, where moyo reports `8b` and
  spglib reports `8a` for the same Hall-525 structure (origin shift). Since
  the 2026-11-11 audit the wrapper matches the letters on the input
  coordinates when the input basis is conventional
  (`xtalforge-symmetry/src/wyckoff_origin.rs`), so xtalforge reports `8a` like
  spglib; the e2e battery checks the letter counts of all eight samples;
- standardized lattice parameters: agrees for Si, NaCl, MOF-5 and ZIF-8;
  differs for graphite, quartz, the triclinic edge case and the hexagonal edge
  case. The differences are conventions, not detection errors: moyo refines
  the standardized cell and may choose a different setting/orientation, while
  spglib idealizes to its own standard setting. The test records these four
  cases as known exceptions.

The comparison uses `Setting::Spglib`; this is the setting the wrapper exposes
by default to match Crystal-Builder's spglib behaviour.

## Verification notes

- `std_cell` refines both the lattice and the positions to the detected
  symmetry; it does not only apply the transformation. This matches spglib's
  idealization and is what the Phase 2 API exposes through `standardize`.
- `moyo::base::Cell` takes atomic numbers (or integer species ids), not
  occupancies. Partial occupancy must therefore be encoded as distinct species
  (for example a split-site species id per component) before detection; the
  wrapper performs this conversion and documents the species mapping in
  `xtalforge-symmetry::convert`.
- `Setting::Spglib` reproduces spglib's Hall number on the audit corpus;
  `Setting::default` chooses origin choice 2 for Si (Hall 526) instead of
  spglib's 525. The origin choice also explains moyo's `8a`/`8b` Wyckoff
  difference, which the wrapper now resolves in the input origin.

## Gap decisions

| Gap | Decision |
|---|---|
| Niggli, Delaunay, Minkowski private | Open an upstream issue/PR to export them; meanwhile port into `xtalforge-symmetry` with attribution (Task 2.6.2). |
| Symmetrization private | Use `MoyoDataset::std_cell` for Phase 2; if a free symmetrization is needed, port the algorithm with attribution (Task 2.6). |
| Brillouin zone and k-path | Implement the Wigner-Seitz construction and band paths ourselves (Task 2.10), because Crystal-Builder gets them from ASE `Cell.bandpath`. |
| Molecular point groups | Adapt `moyo::identify::PointGroup` to non-periodic rotations; if the adapter is insufficient, implement the molecular point-group detection (Task 2.8.3). |
| Origin/setting differences | Expose `Setting::Spglib` as the default; document the origin choice and Wyckoff-letter difference in the API page. |

## Build, target and size

- `cargo check -p xtalforge-symmetry --target wasm32-unknown-unknown`: passes.
- `maturin build` of the Python extension: passes.
- Dependency tree: `moyo 0.22` uses `nalgebra 0.35` (the workspace version),
  `serde`, `strum 0.28`, `itertools 0.15`, `log`, `rustc-hash`,
  `union-find`; no duplicate `nalgebra`, `serde` or `thiserror`.
- Release compile of `xtalforge-symmetry` with moyo on this machine: about
  45 seconds; the incremental check is much faster.
- WASM package: the full release is 3941 KiB at 1.1.0, including the symmetry
  tables, force fields, analysis, builders, powder surfaces and the agent
  session, inside the 3968 KiB budget checked by `scripts/check-wasm-size.sh`.

## Isolation

`xtalforge-symmetry` defines:

- `SpaceGroup`: number, Hall number, Hermann-Mauguin symbol, setting and
  operations, independent of `moyo` types;
- `SymOp`: integer rotation and fractional translation;
- `SymmetryDataset`: standardized cell, primitive standardized cell, orbits,
  Wyckoff sites and transformations;
- `WyckoffSite`: letter, multiplicity and site symmetry.

Conversions live in one private module; the public API never names a `moyo`
type. The audit test is the only place that currently imports `moyo`.
