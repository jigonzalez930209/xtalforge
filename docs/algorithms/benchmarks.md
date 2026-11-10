# Benchmarks

Numbers below are from a Linux x86-64 release build (`cargo build --release -p
xtalforge-cli`) of the Phase 1 CLI. Each row is the median of ten runs of the
whole process (startup, read, work, print) with output discarded; the minimum
is shown to indicate noise. Wall-clock values are rounded to 0.1 ms.

| Structure | Command | Median | Minimum |
|---|---|---|---|
| Si (8 P1 atoms) | `info` | 2.6 ms | 2.2 ms |
| Si | `bonds` | 2.5 ms | 2.1 ms |
| Si | `rings` | 2.7 ms | 2.4 ms |
| MOF-5 (asymmetric unit) | `info` | 2.7 ms | 2.2 ms |
| MOF-5 | `bonds` | 3.5 ms | 2.7 ms |
| MOF-5 | `rings` | 3.3 ms | 2.9 ms |
| MOF-5 | `supercell 2 2 2` | 5.5 ms | 5.2 ms |

Notes:

- `bonds` and `rings` on MOF-5 used to be dominated by the direct periodic
  neighbor search (66.8 ms and 66.7 ms in the pre-linked-cell table). The
  linked-cell index of Task 1.10.2 brings both to about 3.5 ms while keeping the
  output identical; the scaling benchmark is on the
  [neighbor-search page](./neighbors.md).
- The `supercell` command now uses the same fractional grid to drop boundary
  duplicates: the 2x2x2 MOF-5 case fell from 37.5 ms to 5.5 ms.
- The CLI `supercell` command replicates the asymmetric unit; callers that
  want every atom of the crystal expand to P1 first (the Python and WASM
  bindings do).
- Release WebAssembly (`wasm-opt -Oz`) is 3771 KiB after the Phase 8 agent
  surface, inside the 3840 KiB budget enforced by `scripts/check-wasm-size.sh`.
- The workspace test suite (Rust, 38 WASM node tests and 41 pytest cases) runs
  in seconds; see the phase status page for the exact list.

These are process-level numbers, useful for the release gate. Function-level
criterion benchmarks are added with the analysis phases, where the hot paths
(PXRD, grids, force fields) justify the dependency.

## Symmetry detection

Release CLI, median of five whole-process runs at `--symprec 0.001`:

| Structure | Command | Median | Minimum |
|---|---|---|---|
| Si (8 P1 atoms) | `symmetry detect` | 3.3 ms | 3.1 ms |
| MOF-5 (asymmetric unit) | `symmetry detect` | 20.7 ms | 19.6 ms |

The larger MOF structures in the corpus take seconds to minutes at very small
tolerances; the corpus sweep is an ignored test and the release gate uses the
small-corpus set. The full release package, including the symmetry search, the
force fields, analysis, builders, powder and agent surfaces, is 3771 KiB inside
the 3840 KiB budget.

## Function-level benchmarks (Task 8.6.1)

Three process-level benchmarks print the median of repeated runs. They run
with `cargo bench` and `scripts/check-bench-budget.sh` fails on a large
regression (the ceilings are deliberately loose: they catch a 10x change, not
machine noise).

| Benchmark | Command | Median |
|---|---|---|
| Symmetry, 8 atoms | `cargo bench -p xtalforge-symmetry --bench detection` | 0.75 ms |
| Symmetry, 27 atoms | same | 0.72 ms |
| Symmetry, 64 atoms | same | 1.17 ms |
| UFF methane, build + energy | `cargo bench -p xtalforge-ff --bench uff` | 0.043 ms |
| UFF MOF-5, 7 sites / 424 P1 atoms | same | 21.5 ms |
| PXRD MOF-5 | `cargo bench -p xtalforge-analysis --bench pxrd` | 308 ms |
| Porosity MOF-5, 0.8 A grid | same | 2.9 s |

The porosity number is dominated by the `1/spacing^3` grid and is why the CLI
default is 0.35 A only for small cells and the report value is chosen
deliberately. The Python timings recorded for the MACE bridge are on the
[machine-learning page](../api/ml#accuracy-and-speed); they are not comparable
with UFF, which is classical.

## Determinism and robustness (Tasks 8.6.2-8.6.4)

- `crates/xtalforge-core/tests/determinism.rs` runs the same bond graph,
  neighbor list and supercell repeatedly and in four threads and compares the
  exact output.
- `crates/xtalforge-io/tests/robustness.rs` mutates valid CIF, POSCAR, XYZ,
  GEN, LAMMPS, SHELX and pmg-json inputs and a project archive with a
  deterministic generator; a parser may accept or reject a mutation but must
  not panic.
- `crates/xtalforge-core/tests/limits_audit.rs` exercises every budget at and
  beyond its soft and hard limits and the warn-only profile.
- The WASM build is single-threaded by default. A threaded build needs
  cross-origin isolation and is not part of the default package; the
  `parallel` feature changes execution strategy only, never results.
