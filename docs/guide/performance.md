# Performance

This page records the reproducible benchmarks shipped with `1.8.0`. The suite
runs the same crystallographic tasks through the xtalforge Python wheel and
through gemmi, spglib, ASE and pymatgen when those packages are installed, and
writes a deterministic JSON document with the hardware, versions, methodology
and every sample. The raw result is
[`benches/compare/results/v1.8.0.json`](https://github.com/jigonzalez930209/xtalforge/blob/main/benches/compare/results/v1.8.0.json).

## Machine and versions

| Item | Value |
|---|---|
| CPU | AMD Ryzen 9 6900HX, 8 cores / 16 threads |
| Memory | 26.68 GiB |
| OS | Linux 7.0.0-43-generic, x86-64, glibc 2.43 |
| Python | 3.14.4 |
| Rust | 1.99.1 |
| xtalforge | 1.8.0 |
| gemmi | 0.7.5 |
| spglib | 2.8.0 |
| ASE | 3.29.0 |
| pymatgen | 2026.10.24 |
| NumPy | 2.5.4 |

The machine is a laptop; the CPU governor is not pinned. Absolute numbers move
with temperature and background load, which is why the committed regression
gate compares the fastest sample rather than the median.

## Methodology

The harness lives in `benches/compare/`:

- `run.py` runs the cases and writes the JSON document.
- `suite.py` defines the inputs, reference implementations and repetitions.
- `compare.py` compares a new result with `baseline.json`.
- `scripts/check-bench.sh` is the local gate.

Every case receives the same in-memory CIF text. The timed region starts before
that library parses the text, so a CIF read, a P1 expansion or a PXRD pattern
includes its own parser. The harness runs one warm-up call, then records the
samples with the cyclic collector disabled; the reported value is the median.
The committed regression gate uses the minimum of the same samples, because the
minimum is stable within a few percent between runs on this laptop.

The inputs are part of the repository:

| Task | Input and settings | What is timed |
|---|---|---|
| CIF read | `tests/corpus/prepared/UiO-66.cif`, 114 P1 atoms | CIF text to a parsed structure |
| P1 expansion | `tests/corpus/cod/MOF-5.cif`, 7 asymmetric sites | CIF text to the 424-atom P1 cell |
| Symmetry | same MOF-5 CIF, `symprec = 0.001 A` | CIF text to the space-group number |
| Neighbors | same MOF-5 P1 cell, cutoff 3.0 A | CIF text to the complete pair list |
| PXRD | same MOF-5 P1 cell, Cu K-alpha, 5-80 deg, step 0.02 deg | CIF text to the peak list |
| Supercell | `tests/corpus/simple/Si.cif`, 50 x 50 x 50 | eight-atom cell to 1,000,000 atoms |
| Porosity | same MOF-5 CIF, 0.8 A grid, 1.2 A probe | CIF text to the porosity report |

The one-million-atom supercell is measured in memory by
`cargo bench -p xtalforge-core --bench supercell`, not through the Python
dictionary conversion or CIF serialisation. The reference libraries use their
native supercell operator. The default `standard` limit profile refuses more
than 200,000 atoms; the benchmark uses the `huge` profile, and the CLI accepts
`--limits huge` for the same build.

Two conventions differ between libraries and are not timing artefacts:

- xtalforge reports each periodic pair once; ASE and pymatgen report both
  directions, so their pair counts are twice as large.
- pymatgen's PXRD calculator merges symmetry-equivalent reflections into peaks;
  xtalforge reports reflections. The MOF-5 case has 314 reflections and 152
  merged pymatgen peaks.

No installed reference package provides grid-based geometric porosity, so that
row is xtalforge-only. gemmi has no atom-position space-group finder, PXRD
calculator, supercell builder or grid porosity; spglib has no file reader,
P1 expander, neighbor search, PXRD or supercell operator.

## 1.8.0 results

Median wall-clock time in milliseconds, higher is slower. Each cell is the
median of the recorded samples; the JSON keeps the minimum, maximum and every
sample.

| Task | Input size | xtalforge | gemmi | spglib | ASE | pymatgen |
|---|---:|---:|---:|---:|---:|---:|
| CIF read | 114 atoms | 0.62 | **0.09** | - | 48.34 | 11.53 |
| P1 expansion | 424 atoms | 4.78 | **0.51** | - | 49.41 | 17.05 |
| Symmetry detection | 424 atoms | 69.65 | - | **25.72** | - | 37.07 |
| Periodic neighbors | 424 atoms, 3.0 A | **9.60** | - | - | 55.87 | 123.00 |
| PXRD | 424 atoms, 314 reflections | **407.70** | - | - | - | 1105.50 |
| Supercell | 1,000,000 atoms | 3078.77 | - | - | **462.38** | 13230.32 |
| Porosity | 424 atoms, 0.8 A grid | 3914.01 | - | - | - | - |

The committed JSON was captured while the shared host was also running
unrelated test processes; the medians are therefore conservative. The
`min_ms` fields in the JSON are the more useful comparison, and the regression
gate uses them.

What the table says:

- The neighbor list is the strongest result: the linked-cell index is 5.8x
  faster than ASE and 12.8x faster than pymatgen on the 424-atom MOF-5 cell.
- PXRD is 2.7x faster than pymatgen; the xtalforge pattern has 314 reflections
  against 152 merged pymatgen peaks.
- The one-million-atom in-memory supercell completes in 3.1 s in this capture.
  ASE's simpler atom array is about 6.7x faster; pymatgen materialises Python
  sites and is about 4.3x slower than xtalforge. The xtalforge build still
  enumerates every candidate image, which is the next optimisation target.
- CIF read and P1 expansion are dominated by parser and Python-object work:
  gemmi is fastest, xtalforge is second, and xtalforge is about 3.6x faster
  than pymatgen for P1 expansion.
- Symmetry detection on this large cell is 2.7x slower than the raw spglib call
  and 1.9x slower than pymatgen's wrapper. The xtalforge call also builds the
  operation list, Wyckoff letters, standardized cell and orbit arrays returned
  to Python; the JSON and the [symmetry page](../algorithms/symmetry) record
  that payload.

## Regression gate

`scripts/check-bench.sh` runs the xtalforge rows, compares each against
`benches/compare/baseline.json`, and exits non-zero when the fastest sample is
more than 10 % above the stored baseline. The baseline is deliberately small
and contains only xtalforge rows, so a change in a reference package cannot
fail the gate.

```bash
scripts/check-bench.sh              # run the gate
scripts/check-bench.sh --update     # refresh the baseline on this host
scripts/check-bench.sh --self-test  # exercise the comparison logic only
```

The baseline is machine-specific. Run `--update` on the same host after an
intentional change and review the diff; do not copy it between machines. The
script is part of the local release gate alongside `scripts/check-rules.sh`,
`scripts/sync-version.sh --check` and the cargo checks. The self-test covers a
passing run, a regression, an improvement and missing rows without spending
time on timings. If the first comparison fails, the gate reruns the xtalforge
rows once and compares the fastest sample of both runs; a real regression
fails both runs, while a transient load spike does not.

## Reproducing the full suite

```bash
# Release CLI and wheel used by the suite (the wheel provides xtalforge.*)
maturin build --release -m crates/xtalforge-python/Cargo.toml
.venv/bin/pip install --force-reinstall --no-deps target/wheels/xtalforge-*.whl

# All libraries, all tasks, JSON written to a new file
.venv/bin/python benches/compare/run.py \
  --output benches/compare/results/v1.8.0.json

# A single surface or task
.venv/bin/python benches/compare/run.py --only xtalforge --tasks pxrd
.venv/bin/python benches/compare/run.py --reps 3
```

The suite exits non-zero if a selected case raises; a missing reference package
is reported under `environment.import_errors` and its rows are omitted. The
JSON is sorted and contains no timestamps, so reruns differ only in the timing
numbers and their samples.
