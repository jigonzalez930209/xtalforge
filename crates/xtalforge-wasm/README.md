<p align="center">
  <img src="https://raw.githubusercontent.com/jigonzalez930209/xtalforge/main/assets/xtalforge.png" alt="xtalforge" width="128" height="128">
</p>

# xtalforge-wasm

WebAssembly bindings for **xtalforge**, the pure Rust library for building,
manipulating, analysing and exporting crystal structures.

## Install

```bash
npm install xtalforge-wasm
```

## Quick start

```js
import { Structure, scene_buffers, pxrd_buffers } from "xtalforge-wasm";

const structure = Structure.fromText(cifText, "cif");
console.log(structure.formula, structure.nAtoms);

const cell = structure.supercell(2, 2, 1);
cell.positions;       // Float64Array, flat xyz, Angstrom
cell.fracPositions;   // Float64Array, flat xyz
cell.cell;            // Float64Array, row-major 3x3
cell.numbers;         // Uint8Array
cell.toCif();
cell.free();
```

`scene_buffers` and `pxrd_buffers` return typed arrays for rendering and
plotting; `mesh_buffers` returns the triangle mesh. The typed object interface
is `ts/xtalforge-object.d.ts` in the repository.

The package also ships source helpers in `js/`: `structure-registry.js`
registers handles with `FinalizationRegistry`, `worker.js` runs heavy
operations in a Web Worker, and `node.js` provides Node file helpers.

## Links

- [Documentation](https://jigonzalez930209.github.io/xtalforge/)
- [Object API guide](https://jigonzalez930209.github.io/xtalforge/api/objects)
- [Repository](https://github.com/jigonzalez930209/xtalforge)
- [Changelog](https://github.com/jigonzalez930209/xtalforge/blob/main/CHANGELOG.md)

## License

MIT OR Apache-2.0.
