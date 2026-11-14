# From a MOF CIF to a DFT band gap

This tutorial writes a reproducible VASP input set for a MOF, explains how to
run it, and reads the band gap back. The engine itself stays external; the
input writing, parsing and analysis are pure Rust and available in the CLI,
Python and WASM.

## 1. Write the inputs

```bash
xtalforge engines inputs vasp/mof tests/corpus/prepared/MOF-5.cif --out work/mof5
```

The directory contains `INCAR`, `KPOINTS`, `POSCAR` and `POTCAR.spec`. The
`POTCAR.spec` lists the POTCAR names in POSCAR order; the licensed POTCAR
files are never generated or bundled. A static band calculation uses:

```bash
xtalforge engines inputs vasp/band tests/corpus/prepared/MOF-5.cif \
  --out work/mof5-band --overrides '{"incar":{"ICHARG":11,"ISMEAR":0}}'
```

Overrides are merged deterministically. Unknown keys are an error unless the
preset declares the branch open or the request sets `"allow_unknown": true`.

## 2. Run the engine

```bash
xtalforge engines script slurm "cd work/mof5-band && mpirun vasp_std" \
  --out work/mof5-band --queue main --nodes 2 --walltime 24:00:00 --dry-run
```

Without `--dry-run` the same command submits through `sbatch` when the runtime
probe finds it. A failed run can be classified:

```bash
xtalforge engines failures
```

The catalogue names `scf_not_converged`, `zbrent`, `walltime_exceeded`,
`out_of_memory` and `missing_pseudopotential`, each with a corrective action.

## 3. Read the bands

```bash
xtalforge engines bands vasp-vasprun work/mof5-band/vasprun.xml --json
```

The result carries the labeled bands and the direct/indirect gap. The Python
surface returns the same JSON-shaped dictionary:

```python
import xtalforge

bands = xtalforge.engine_bands(open("vasprun.xml").read(), "vasp-vasprun")
print(bands["gap"])
```

The browser uses `engine_bands_info` after fetching the file; the WASM build
does not run the engine or read a filesystem.

## 4. Density of states

```bash
xtalforge engines dos vasp-doscar work/mof5-band/DOSCAR --sigma 0.1 --json
```

The d-band centre and the smeared total are in the result. See
[Electronic structure](../algorithms/electronic-structure.md) for the gap
classification, the effective-mass fit and the smearing normalisation.
