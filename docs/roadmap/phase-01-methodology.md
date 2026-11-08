# Phase 1 — Module development methodology

Phase 1 builds the crystal data model: elements, lattice, limits, sites,
structures, P1 expansion, supercells, slabs, neighbor search, bonding, rings,
functional groups, measurements, and the CIF/POSCAR/XYZ readers and writers.
This page fixes how each module is studied, designed, implemented, tested,
exposed and documented before the task work starts.

## 1. Study the reference before writing code

For every subtask:

1. Read the named Crystal-Builder file completely, not only the functions the
   roadmap lists. Record its public functions, data classes, invariants and
   tests in the task note.
2. Read the matching sci-form module when the roadmap names one. Record what
   is newer, better tested or more accurate.
3. Identify the data (tables, parameters, citations) and move it to a data file
   with provenance before writing the algorithm.
4. Write the module split plan: files, ownership, public types and the
   dependency edges. No file starts near 500 lines.
5. Check the golden data and shared vectors that will validate the module.

The study output is a short section in the task's algorithm page or in the
pull request description; it is not a separate essay.

## 2. Design the public boundary first

- One public type per concept, owned by the domain crate:
  `Element`, `Lattice`, `Site`, `Structure`, `P1Cell`, `SymOp`, `BondGraph`.
- Errors are typed variants of the crate `Error` with a stable `code()`.
  Invalid external data never panics.
- Units are Angstrom and degrees at the boundary; radians are internal.
- Fractional coordinates are `[f64; 3]`; the lattice matrix rows are the
  vectors a, b, c; fractional to Cartesian is `frac · M`.
- The Hall symbol is the identity of a space group; a number alone is not
  enough to fix the setting.
- `Structure` is an asymmetric unit plus a space group; `P1Cell` is derived
  data with provenance (`site_idx`, `op_idx`, `tau`).
- Mutation is crate-private until the command layer exists. A
  `StructureBuilder` creates values; pure functions return new ones.
- Large objects cross language boundaries as handles, not copied records.

## 3. Implementation order

The dependency graph decides the order; tasks may overlap only where the graph
allows it.

```
elements ─┬─> lattice ─┬─> structure/site ─┬─> p1 ─┬─> supercell ─> slab
          │            │                   │       ├─> neighbors ─> bonding
          │            └─> limits          │       │                   │
          │                                │       │                   ├─> rings
          │                                │       │                   └─> groups
          │                                │       └─> measures/scalars/describe
          │                                └─> serialization schema
          └─> pyykko/isotopes

CIF/POSCAR/XYZ ──> structure ──> P1 ──> every analysis
```

- Elements first because every other module names atoms.
- Lattice second; it owns geometry and is used by neighbours, slabs and P1.
- Limits third because every allocation is checked before it happens.
- Structure/site next; the schema follows the types.
- P1 needs `SymOp`. Until Task 2.1 exposes moyo operations, use the explicit
  operation list read from the CIF or supplied by the builder.
- Neighbours precede bonding; bonding precedes rings and groups.
- Readers and writers start as soon as the structure model exists; they are
  the source of the corpus tests.
- Surfaces and docs close each task; they are not a separate phase.

## 4. Data files

| Data | File | Owner |
|---|---|---|
| Periodic table, colours, valences, electronegativity | `crates/xtalforge-core/data/elements.json.gz` | Task 1.1 |
| Pyykkö single/double/aromatic radii | same bundle | Task 1.1.3 |
| Isotope masses | same bundle | Task 1.1.3 |
| UFF parameters | `crates/xtalforge-ff/data/uff.json.gz` | Task 4.2 |
| Scattering factors | `crates/xtalforge-analysis/data/scattering.json.gz` | Task 5.1 |

The canonical source is the Crystal-Builder table; the generator records the
commit and a SHA-256 checksum of the decompressed payload. Data is embedded
with `include_bytes!` from the owning crate so `cargo package` carries it.
Large external bundles are loaded lazily, never embedded.

## 5. Tests

- Every public function and every error path has an integration test in
  `crates/<crate>/tests/` (AGENTS.md section 2.2).
- Golden tests use `xtalforge-testkit` and the data in `tests/golden/`.
- Shared vectors in `tests/vectors/` run on Rust, CLI, Python and WASM.
- Property tests cover symbol parsing, round-trips and atom-set invariance.
- A regression test names the structure, the option and the expected value;
  a tolerance change is a reviewed change with a reason.
- Readers are fuzzed (Task 1.15.5); the default suite stays fast.

## 6. Surfaces

For every user-visible operation:

1. Rust API in the facade under a stable module name.
2. CLI verb with `--json`, exit codes and shared vectors.
3. Python function or handle with stubs checked by `mypy --strict`.
4. WASM export with generated TypeScript definitions.
5. Parity matrix row and `scripts/check-parity.sh` green.

A surface that cannot exist (filesystem in WASM, for example) is recorded as a
justified exception, not silently skipped.

## 7. Documentation

- API page with the four surfaces in code groups.
- Algorithm page when the module has a formulation, a convention or a
  numerical caveat (lattice conventions, P1 provenance, supercell enumeration,
  slab basis, cell list, bond perception, primitive rings, CIF symmetry
  priority).
- Validation row in `docs/algorithms/validation.md`.
- Roadmap checkbox, parity matrix and `CHANGELOG.md` `Unreleased` entry in the
  same commit.

## 8. Verification gate

Before a Phase 1 subtask is ticked:

```bash
scripts/check-rules.sh
cargo fmt --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
scripts/check-parity.sh
scripts/check-api-docs.sh
(cd docs && npm run docs:build)
```

Plus the Python and WASM tests when the surface changed. CI runs the same
commands once Actions are enabled.

## 9. Commit and review

- One directory per commit (AGENTS.md section 2.1).
- Subject `feat(<scope>): ...` or `ref(<scope>): ...`; body lists the completed
  subtasks from hardest to easiest.
- A task is not ticked until code, data, tests, surfaces and docs are merged.
- A deliberate difference from Crystal-Builder is recorded in the algorithm
  page before the task is ticked.

## 10. Phase 1 execution order

1. Task 1.1 — elements, symbol parsing, Pyykkö radii, surfaces.
2. Task 1.2 — lattice, conversions, transforms, validation.
3. Task 1.3 — limits.
4. Task 1.4 — site, structure, molecule, properties.
5. Task 1.5 — serialization schema.
6. Task 1.6 — P1 expansion and symmetry operations.
7. Tasks 1.7–1.9 — supercells, slabs, finite cell and transforms.
8. Tasks 1.10–1.13 — neighbors, bonding, rings, groups.
9. Tasks 1.14–1.16 — measures/scalars/describe, CIF, POSCAR/XYZ.
10. Task 1.17 — CLI, Python and WASM for Phase 1.
11. Task 1.18 — docs, demo, benchmarks and the `v0.0.1` prerequisites.

Task 2.1 (moyo operations) runs in parallel after Task 1.2 because P1 needs it.
