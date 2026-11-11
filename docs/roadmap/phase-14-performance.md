# Phase 14 — Performance and scale

**Goal:** be the fastest option for every common operation and handle
structures from one atom to millions, with published, reproducible benchmarks
against spglib, gemmi, ASE and pymatgen.

**Release gate:** `v1.7.0`.

---

## Task 14.1 — Benchmarks as a contract

- [ ] 14.1.1 Comparative benchmark suite
  Description: `benches/compare/` runs the same tasks (CIF read, P1, symmetry,
  neighbors, PXRD, supercell 10^6 atoms, porosity) in xtalforge and the
  reference libraries; results published on the docs site per release.
- [ ] 14.1.2 Regression gate
  Description: `scripts/check-bench.sh` fails the local CI gate when a tracked
  benchmark regresses more than 10 %.

## Task 14.2 — Algorithms at scale

- [ ] 14.2.1 Structure-of-arrays storage for P1 cells and neighbor lists.
- [ ] 14.2.2 Streaming readers and writers for very large CIF/XYZ/LAMMPS
  files (constant memory per frame).
- [ ] 14.2.3 Parallel UFF energy/forces with a fixed-order reduction.
- [ ] 14.2.4 Lazy symmetry expansion: operations over the asymmetric unit
  without materialising P1 where the result allows it.
- [ ] 14.2.5 Configurable limits presets for million-atom work (`huge`).

## Task 14.3 — Binary size and startup

- [ ] 14.3.1 Split WASM packages (`core`, `symmetry`, `analysis`, `build`)
  so a viewer loads only what it uses; size budget per package in
  `scripts/check-wasm-size.sh`.
- [ ] 14.3.2 Lazy data bundles (download on first use in WASM, compressed in
  the crates).
- [ ] 14.3.3 CLI cold start under 20 ms for `info`.

## Task 14.4 — Surfaces and docs

- [ ] 14.4.1 Performance page with methodology, hardware and numbers.
