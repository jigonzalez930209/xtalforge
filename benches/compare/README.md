# Comparative benchmarks

Task 15.1.1 and 15.1.2. The suite runs the same CIF read, P1 expansion,
symmetry, neighbor list, PXRD, 1,000,000-atom supercell and porosity tasks
through xtalforge and through gemmi, spglib, ASE and pymatgen when installed.

## Commands

```bash
# Full run, committed result shape
.venv/bin/python benches/compare/run.py \
  --output benches/compare/results/v1.8.0.json

# Only the xtalforge rows, with a custom repetition count
.venv/bin/python benches/compare/run.py --only xtalforge --reps 5

# Regression gate against baseline.json (10 % limit)
scripts/check-bench.sh
scripts/check-bench.sh --update
scripts/check-bench.sh --self-test
```

`run.py` writes a sorted JSON document with `schema_version`, the environment
(CPU, memory, OS, Python, Rust and every library version), task definitions,
the repetition count, and one record per task/library pair. Records contain
`median_ms`, `min_ms`, `max_ms`, `samples_ms` and a small `output` summary that
proves the timed call produced the expected atom, pair, reflection or volume
count. There are no timestamps, so a rerun differs only in measured values.

## Rules

- The same in-memory CIF text is used for every library; parsing is inside the
  timed region.
- Each case gets one warm-up and the configured repetitions. The cyclic
  collector is disabled during the timed loop; a full collection runs once
  before it.
- The regression gate compares the minimum sample against the stored minimum.
  The median is reported for the docs and is less stable on a laptop under
  thermal and frequency changes. If the first pass fails, `check-bench.sh`
  runs the xtalforge rows again and compares the faster of the two runs.
- The baseline is machine-specific. Refresh it on the same host with
  `scripts/check-bench.sh --update` after an intentional code or dependency
  change.
- The 1,000,000-atom xtalforge case calls
  `cargo bench -p xtalforge-core --bench supercell` and is in-memory; the
  Python wheel's default `standard` limits stop at 200,000 atoms.
- No installed reference package has grid-based porosity, so that row is
  xtalforge-only.

## Files

| File | Purpose |
|---|---|
| `run.py` | CLI, environment capture and JSON writer |
| `suite.py` | Inputs, reference implementations, repetitions and timing |
| `compare.py` | Baseline comparison, `--update` and `--self-test` logic |
| `baseline.json` | Small, machine-specific xtalforge-only regression baseline |
| `results/v1.8.0.json` | Full 1.8.0 result used by the docs page |

The docs page is [`docs/guide/performance.md`](../../docs/guide/performance.md).
