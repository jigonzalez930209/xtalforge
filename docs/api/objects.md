# Object API, arrays and batches (Phase 12)

The 1.x text/dict API remains supported. Phase 12 adds object-oriented
handles (`structure_handle`), NumPy/typed-array input and output
(`structure_from_arrays`, `typed_arrays`), typed results (`typed_results`)
and batch helpers (`batch_api`) on top of it.

## Python `Structure`

`xtalforge.Structure(text, format_hint="")` and `Structure.from_cif(text)`
parse once. The handle exposes the common operations as methods:

```python
s = xtalforge.Structure(cif)
s.supercell(2, 2, 2).slab(1, 1, 1, layers=2, vacuum=10).write("slab.cif")
s.positions, s.frac_positions, s.cell, s.numbers     # NumPy arrays
s.to_dict(), s.to_json(), s.to_format("poscar")
s.set_positions(frac), s.add_hydrogens(), s.ff_energy()
```

Free functions with a text first argument accept a `Structure` too; the
wrapper converts it to its JSON document. The legacy return shape is a dict;
`Structure.to_dict()` is the same shape.

## Typed results

`SymmetryDataset`, `PxrdPattern`, `PorosityReport` and `OptimizeResult` are
dict subclasses with `to_dict()`, a rich `__repr__` and `_repr_html_`. They
are accepted anywhere the old dict was, so 1.x code keeps working. The
`mypy --strict` stub has no `Any` in a public signature
(`scripts/check-python-stubs.sh` plus the stub test).

## NumPy in and out

- `positions` (cartesian, Angstrom), `frac_positions`, `cell`, `numbers`.
- `Structure.from_arrays(cell, numbers, frac=... | cart=...)` validates the
  shapes and atomic numbers with `InputError`.
- `set_positions(array)` returns a new P1 structure.
- PXRD `two_theta`/`intensity`, scene buffers and meshes are arrays in the
  WASM object API (below); Python returns lists/NumPy where the underlying
  result is sampled data.

## Batch helpers

`xtalforge.batch.read(paths)`, `batch.detect_symmetry(texts)` and
`batch.pxrd(texts)` return lists in input order; the functions are
deterministic and do not reorder results.

## JavaScript / TypeScript object API

The WASM package exports a `Structure` class:

```js
const s = Structure.fromText(cif, "cif");
const cell = s.supercell(2, 2, 2);
cell.positions;        // Float64Array, xyz
cell.fracPositions;    // Float64Array
cell.cell;             // Float64Array, row-major 3x3
cell.numbers;          // Uint8Array
cell.toCif(); cell.toPoscar(); cell.toJson();
```

`scene_buffers`, `pxrd_buffers` and `mesh_buffers` return
`Float64Array`/`Float32Array`/`Uint8Array`/`Uint32Array` views for zero-copy
rendering. The typed interface is `crates/xtalforge-wasm/ts/xtalforge-object.d.ts`
and is checked by `scripts/check-ts-any.sh`; the legacy text/dict exports keep
their generated loose signatures until 2.0 (recorded in the parity matrix).

`crates/xtalforge-wasm/js/structure-registry.js` registers handles with
`FinalizationRegistry`, `js/worker.js` runs supercell/slab/P1/PXRD/symmetry in
a Web Worker, and `js/node.js` provides Node `read`/`write` file helpers.

## Parse-once gain

A handle owns the parsed Rust value; repeated calls (supercell, slab, PXRD,
symmetry) do not tokenize or parse the CIF again. The measured benchmark and
the playground migration are tracked by Task 12.4.1; the text API keeps the
old behaviour for compatibility.

## Ecosystem adapters

`Structure.from_ase`, `from_pymatgen`, `from_rdkit`, `to_ase`,
`to_pymatgen`, `to_rdkit`, `to_phonopy_atoms` and
`xtalforge.get_symmetry_dataset` are documented in the
[interoperability page](./interop.md).
