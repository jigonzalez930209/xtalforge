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
- Release WebAssembly (`wasm-opt -Oz`) is about 2905 KiB after the Phase 3-7
  surfaces, inside the 3584 KiB budget enforced by `scripts/check-wasm-size.sh`.
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
force fields, analysis, builders and powder surfaces, is about 2905 KiB inside
the 3584 KiB budget.
