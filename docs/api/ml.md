# Machine-learning potentials

The ML backend is the Python bridge recorded in
[the backend decision](../algorithms/ml-backend). Rust has no native MACE,
ORB-v3 or MatterSim implementation; the models run through ASE and PyTorch in
the Python package, and the Rust `ml` feature stays a declared no-op.

- Rust: `xtalforge::ff::models` (catalogue, cache, checksums)
- CLI: `xtalforge models {list,fetch,verify,remove}`
- Python: `MlCalculator`, `ml_backends()`
- WASM: n/a (no Python, no torch); recorded as an exception in the parity
  matrix.

## Python

```python
from mace.calculators import mace_mp
import xtalforge

calculator = mace_mp(model="small", device="cpu", default_dtype="float64")
ml = xtalforge.MlCalculator(calculator)          # periodic=True by default

energy = ml.energy(structure_json)                # eV
result = ml.energy_forces(structure_json)         # JSON string
```

`MlCalculator` wraps **any** ASE calculator, so MACE, ORB-v3 and MatterSim use
the same class: build the backend object the way its own documentation says
and pass it in. `energies` and `forces` are eV and eV/Angstrom. The stress is
assembled from ASE's six Voigt components `(xx, yy, zz, yz, xz, xy)` into a
3x3 tensor and converted from eV/Angstrom^3 to GPa.

`ml_backends()` reports whether `torch`, `ase`, `mace`, `orb_models` and
`mattersim` are importable, with their versions.

## Model cache

```bash
xtalforge models list
xtalforge models fetch mace-mp-0 --url <artifact-url>
xtalforge models verify mace-mp-0
xtalforge models remove mace-mp-0
```

The catalogue records the license and homepage of the three families:

| Name | License | Homepage |
|---|---|---|
| `mace-mp-0` | MIT | https://github.com/ACEsuit/mace |
| `orb-v3` | Apache-2.0 | https://github.com/orbital-materials/orb-models |
| `mattersim` | MIT | https://github.com/microsoft/mattersim |

The artifact URL is chosen by the user because the file and its checksum are
versioned upstream. `fetch` downloads with `curl`, records the SHA-256, size,
license and URL in `<name>.json`, and `verify` recomputes the digest. The cache
is `$XTALFORGE_MODEL_CACHE`, then `$XDG_CACHE_HOME/xtalforge/models`, then
`~/.cache/xtalforge/models`. Weights are never bundled with the crate.

Example: the MACE-MP-0 small model from
`https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0/2023-12-10-mace-128-L0_energy_epoch-249.model`
is 32,581,838 bytes with SHA-256
`2ddb079cee0e131eaaf6912ba581b394551ead283e95c99cfe78c605d10b5736`.

## Accuracy and speed

`crates/xtalforge-python/tests/fixtures/ml/mace_mp0_small.json` records MACE
energies, forces and stress for H2, water and methane in a 10 A box, with the
model checksum. `test_ml.py` recomputes them through `MlCalculator` and
requires 1e-6 eV on the energy, 1e-5 eV/Angstrom on the forces and 1e-6 GPa on
the stress. Measured on one CPU with float64:

| Case | Atoms | Energy (eV) | Time (s) |
|---|---|---|---|
| H2 | 2 | -6.555498 | 0.10 |
| water | 3 | -14.039795 | 0.20 |
| methane | 5 | -23.768234 | 0.06 |

The load of the model dominates a short run; the per-call time is for one
energy and force evaluation after the model is loaded. ORB-v3 and MatterSim
were not installed in the validation environment, so the bridge supports them
through ASE but the recorded golden comparison covers MACE only. This is
stated rather than implied.

## Limitations

- No GPU validation; the recorded run is CPU float64.
- The ASE calculator owns the model's cutoffs, dtype and device; the bridge
  does not change them.
- A structure is passed as the P1 cell. A molecule should use
  `periodic=False`; a periodic framework keeps the cell.
- The stress sign and convention are ASE's; do not mix them with the native
  UFF stress without checking.
