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
