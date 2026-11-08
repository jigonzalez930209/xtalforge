# Benchmarks

Numbers below are from a Linux x86-64 release build (`cargo build --release -p
xtalforge-cli`) of the Phase 1 CLI. Each row is the median of ten runs of the
whole process (startup, read, work, print) with output discarded; the minimum
is shown to indicate noise. Wall-clock values are rounded to 0.1 ms.

| Structure | Command | Median | Minimum |
|---|---|---|---|
| Si (8 P1 atoms) | `info` | 1.7 ms | 1.3 ms |
| Si | `bonds` | 1.9 ms | 1.7 ms |
| Si | `rings` | 1.9 ms | 1.7 ms |
| MOF-5 (asymmetric unit) | `info` | 1.8 ms | 1.5 ms |
| MOF-5 | `bonds` | 66.8 ms | 64.5 ms |
| MOF-5 | `rings` | 66.7 ms | 63.9 ms |
| MOF-5 | `supercell 2 2 2` | 2.1 ms | 1.6 ms |

Notes:

- `bonds` and `rings` on MOF-5 are dominated by the direct periodic neighbor
  search over the expanded P1 cell. The linked-cell acceleration of Task
  1.10.2 must preserve the current output exactly; the benchmark table is the
  regression reference.
- The CLI `supercell` command replicates the asymmetric unit; callers that
  want every atom of the crystal expand to P1 first (the Python and WASM
  bindings do).
- Release WebAssembly (`wasm-opt -Oz`) is 438 KiB, inside the 512 KiB budget
  enforced by `scripts/check-wasm-size.sh`.
- The workspace test suite (Rust, 20 WASM node tests and 27 pytest cases) runs
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
small-corpus set. The WebAssembly package with the symmetry search and the
space-group tables is about 952 KiB, inside the 1 MiB budget.
