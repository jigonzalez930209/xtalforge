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

## Phase 21 native surface

The `xtalforge-ml` crate (facade feature `ml`) adds the preprocessing and
inference layer. Optional features: `ml-onnx` (tract) and `ml-parquet`
(Arrow/Parquet).

- Rust: `xtalforge::ml`
- CLI: `xtalforge ml graph|soap|acsf|mbtr|mof-features|split|import|predict|select|active-workflow`
- Python: `crystal_graph_arrays`, `crystal_graph_numpy`, `soap_descriptor`,
  `acsf_descriptor`, `mbtr_descriptor`, `mof_ml_features`, `ml_split`,
  `ml_import`, `ml_predict`, `ml_predict_dataset`, `ml_select`,
  `ml_active_workflow`
- WASM: `crystal_graph_info`, `soap_descriptor_info`, `acsf_descriptor_info`,
  `mbtr_descriptor_info`, `mof_ml_features_info`, `ml_split_info`,
  `ml_import_info`, `ml_predict_info`, `ml_select_info`,
  `ml_active_workflow_info`

### `crystal_graph`

`crystal_graph(structure, GraphOptions { cutoff, max_neighbors })` expands the
P1 cell in canonical order, uses the exact linked-cell neighbour search and
stores both directions of every edge with the target lattice image, the
Cartesian edge vector and the distance. `max_neighbors` keeps the shortest
outgoing edges per node and sets `truncated`. `crystal_graph_arrays` returns
`atomic_numbers`, `positions`, `edge_index`, `edge_vectors`, `distances`,
`edge_images` and `lattice`; Python converts them to NumPy arrays with
`crystal_graph_numpy`, and WASM returns JSON arrays.

### `soap_descriptor`

`soap_power_spectrum` implements the SOAP power spectrum with the GTO radial
basis and real spherical harmonics through degree five. Options are `cutoff`,
`n_max`, `l_max`, `sigma` and `average`. The layout is `l`, `k`, `k'` with
`k <= k'`. `soap_feature_count` reports `n_max (n_max + 1) / 2 * (l_max + 1)`.
The recorded DScribe fixture fixes the layout and feature count; the numerical
normalization is documented in [the algorithm page](../algorithms/ml).

### `acsf_descriptor` and `mbtr_descriptor`

`acsf` implements the Behler-Parrinello G2/G4 functions with a cosine cutoff.
`mbtr` implements the k=1 pair and k=2 angle blocks with Gaussian broadening
on a fixed grid. Both return the element order and feature layout.

### `mof_ml_features`

Revised autocorrelations for `Z`, electronegativity, covalent radius and mass
at bond-path lengths one to three, plus the pore geometry vector (PLD, LCD,
ASA in m2/g and m2/cm3, accessible volume, void fraction, channel dimension,
density, probe radius). `pores = false` skips the distance grid.

### `ml_predict`

`PropertyModel` is a linear model with optional mean/scale normalization. It
loads from and writes to JSON, implements the `Predictor` trait, and drives
`ml_predict_dataset` / `predict_dataset`. `Ensemble` reports the mean,
standard deviation and member values. With `ml-onnx`, `OnnxModel` loads a
single-input ONNX graph and `OnnxCalculator` exposes it as a force-field
calculator with finite-difference forces and stress.

### `ml_select`

`select_batch` ranks candidates by uncertainty, farthest point in descriptor
space or uncertainty-then-diversity. The result is a sorted index list.

### `ml_active_workflow`

`active_learning_workflow` emits a six-step Phase 16 workflow:
read, search, predict, select, append, write. `active_learning_loop` runs the
same pipeline over an in-memory dataset with caller-supplied featurizer,
uncertainty, compute and retrain closures.

## Limits

- ONNX inference is descriptor-based (`OnnxCalculator`); message-passing
  models stay in the Python ML bridge.
- Parquet stores `id`, structure document JSON and properties JSON as three
  UTF-8 columns; it is lossless but not a nested Arrow schema.
- SOAP follows the published GTO power spectrum; the DScribe absolute
  normalization differs (documented in the algorithm page).
- WASM has no filesystem, so Parquet and checkpointed active learning stay
  host-side.
