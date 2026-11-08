# External engines

## Availability, not existence

An external engine is optional and lives behind the `external` feature. A
probe checks that the binary actually runs (`--version` with a short read and
a timeout) instead of only checking that a file exists. `EngineStatus` is
`Available { path, version }` or `Missing { reason }`; `xtalforge engines`
prints the registry including the reason.

On WASM the probe always reports `Missing` with "not available in the browser";
this is the recorded parity exception for every external engine.

## DFTB+

`write_dftb_hsd(structure, &DftbOptions)` writes a `dftb_in.hsd` with
`Geometry = GenFormat`, the `Hamiltonian = DFTB` block, the Slater-Koster
folder naming and the k-point mesh. `parse_dftb_output` reads `detailed.out`
(or `results.tag`) for the total energy, forces and the stress tensor.
Tests use recorded outputs, so CI needs no binary.

## xTB

`write_xtb_input` writes the XYZ block plus `$chrg`/`$uhf`; `parse_xtb_output`
reads `TOTAL ENERGY`, the charge and the forces when present. GFN1-xTB,
GFN2-xTB and GFN-FF are selected by method name and `xtb` is run through the
same process runner.

## ORCA

ORCA is molecular. `write_orca_input` makes every molecule whole across the
cell faces using the bond graph, then writes `! method basis` and a `* xyz`
block with the charge and multiplicity. `orca_catalogue` lists the supported
functionals and bases. Only inputs are written; ORCA is never launched.

## Quantum adapter

`quantum::adapter(method)` is the in-process path for sci-form GFN0/GFN1/GFN2.
Without the `sci-form` feature it returns a clear typed error telling the user
to rebuild with the feature; the dependency is not vendored, so the feature is
reserved and the capability table (`periodic`, `forces`, `stress`) is what a
future adapter must satisfy. The table documents that stress is not provided.

## Machine-learning potentials

The decision is recorded in [the ML backend page](../algorithms/ml-backend.md):
no native backend in Phases 1-7, no silent model, and the Python bridge is the
only supported path. The `ml` feature stays reserved and a no-op.
