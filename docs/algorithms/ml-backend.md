# Machine-learning potentials: backend decision

Task 4.11 asks for a decision, not an implementation. The implementation is
Task 8.4 (Phase 8), and this page records what was evaluated on 2026-11-08.

## What the reference does

Crystal-Builder runs MACE, ORB-v3 and MatterSim through PyTorch and ASE
(`xtal/ff/ase_engine.py`). Those models are equivariant message-passing
networks: their compute graph contains tensor products, spherical harmonics
and scatter operations that are not a fixed ONNX operator set in the versions
available today.

## Options

| Option | Operator coverage | WASM | Native size | Verdict |
|---|---|---|---|---|
| ONNX Runtime (`ort`) | Fails on the equivariant layers without custom ops | No | Large native library | Rejected |
| `burn` / `candle` | Requires reimplementing the model architecture and weights | Possible for small models | Pure Rust | Rejected for now |
| Python bridge (ASE/PyTorch) | Complete | No | Python environment | Accepted as the only supported path |
| Native Rust models | None of the three architectures is available | — | — | Not feasible now |

## Decision

1. **No native Rust ML backend in Phases 1-7.** The `ml` feature stays a
   reserved, no-op feature so downstream code can compile against the name.
2. **The supported path is the external Python bridge**, described as an
   engine probe and run through the process runner (Task 4.9). It is
   unavailable in WASM, which is recorded as a parity exception.
3. **Revisit at Task 8.4** when ONNX Runtime or a pure-Rust implementation
   covers the equivariant layers of at least one of the three models. The
   evaluation must record accuracy against the reference energies, speed,
   binary size and WASM feasibility before the decision changes.

This keeps the default build pure Rust, keeps the feature honest, and avoids
shipping a model that silently disagrees with the reference.

## Implementation (Task 8.4)

The decision above is implemented as a Python bridge, not as a native backend:

- `xtalforge.MlCalculator` wraps any ASE calculator. It converts a structure
  JSON document to an ASE `Atoms` (P1 cell, `pbc` from the constructor), runs
  the calculator, and returns energy in eV, forces in eV/Angstrom and the
  stress tensor in GPa. The ASE six-component Voigt order is mapped to a 3x3
  tensor and multiplied by 160.21766208 (1 eV/Angstrom^3 in GPa).
- `xtalforge.ml_backends()` probes `torch`, `ase`, `mace`, `orb_models` and
  `mattersim` and reports their versions.
- `xtalforge_ff::models` and the CLI `models` command manage the weights:
  catalogue, license, URL, SHA-256, cache and removal. Weight files are never
  bundled; the user chooses the artifact URL and the checksum is recorded.
- The `ml` Cargo feature remains a declared no-op. A downstream Rust program
  that enables it does not get a native potential; the documented path is the
  Python package.

Validation: `crates/xtalforge-python/tests/fixtures/ml/mace_mp0_small.json`
records MACE-MP-0 (small, float64, CPU) energies, forces and stress for H2,
water and methane. The test reproduces them through `MlCalculator` within
1e-6 eV, 1e-5 eV/Angstrom and 1e-6 GPa. The model is 32,581,838 bytes with
SHA-256 `2ddb079cee0e131eaaf6912ba581b394551ead283e95c99cfe78c605d10b5736`.

ORB-v3 and MatterSim use the same ASE protocol, but neither was installed in
the validation environment, so no golden comparison is claimed for them. This
is the honest subset of Task 8.4.3.
