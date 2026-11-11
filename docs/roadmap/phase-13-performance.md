# Phase 13 — Performance and scale

**Goal:** be the fastest option for every common operation and handle
structures from one atom to millions, with published, reproducible benchmarks
against spglib, gemmi, ASE and pymatgen.

**Release gate:** `v1.6.0`.

---

## Task 13.1 — Benchmarks as a contract

- [ ] 13.1.1 Comparative benchmark suite
  Description: `benches/compare/` runs the same tasks (CIF read, P1, symmetry,
  neighbors, PXRD, supercell 10^6 atoms, porosity) in xtalforge and the
  reference libraries; results published on the docs site per release.
- [ ] 13.1.2 Regression gate
  Description: `scripts/check-bench.sh` fails the local CI gate when a tracked
  benchmark regresses more than 10 %.

## Task 13.2 — Algorithms at scale

- [ ] 13.2.1 Structure-of-arrays storage for P1 cells and neighbor lists.
- [ ] 13.2.2 Streaming readers and writers for very large CIF/XYZ/LAMMPS
  files (constant memory per frame).
- [ ] 13.2.3 Parallel UFF energy/forces with a fixed-order reduction.
- [ ] 13.2.4 Lazy symmetry expansion: operations over the asymmetric unit
  without materialising P1 where the result allows it.
- [ ] 13.2.5 Configurable limits presets for million-atom work (`huge`).

## Task 13.3 — Binary size and startup

- [ ] 13.3.1 Split WASM packages (`core`, `symmetry`, `analysis`, `build`)
  so a viewer loads only what it uses; size budget per package in
  `scripts/check-wasm-size.sh`.
- [ ] 13.3.2 Lazy data bundles (download on first use in WASM, compressed in
  the crates).
- [ ] 13.3.3 CLI cold start under 20 ms for `info`.

## Task 13.4 — Surfaces and docs

- [ ] 13.4.1 Performance page with methodology, hardware and numbers.
