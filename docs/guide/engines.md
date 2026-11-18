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

## Validated engines

The 2026-11-17 campaign ran the engines on this machine, not only the
simulated-executable tests:

- **DFTB+ 25.1** with the `mio-1-1` set on methane: the generated
  `dftb_in.hsd`/`geometry.gen` set runs and parses to `-87.774977538 eV` with
  `0.015703 eV/A` maximum force. The run exposed and fixed the GenFormat,
  species and force-request defects described in the
  [validation campaign](./validation-campaign).
- **xTB 6.8.0** on methane and water: the in-process GFN2 energy agrees within
  `6.1e-6 eV` and the forces within `1e-5 eV/A`.
- **ORCA** remains input-only by design.

The engine probe in `xtalforge capabilities` reports what is installed. A
missing binary is a normal `Missing` status, not a library error.
