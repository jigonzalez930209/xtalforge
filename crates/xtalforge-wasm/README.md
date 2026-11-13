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

## Split packages (Task 15.3.1)

`xtalforge-wasm` is the `full` build. Four smaller packages are built from the
same source with `scripts/build-wasm-packages.sh` (in the repository):

| Package | Adds | Measured (release build) |
|---|---|---|
| `xtalforge-wasm-core` | reading, conversion, supercells/slabs, scene and mesh buffers, documents, `Structure` | 1957 KiB |
| `xtalforge-wasm-symmetry` | space groups, Wyckoff, standardization, reduction, subgroups, Brillouin zone, matcher | 2807 KiB |
| `xtalforge-wasm-analysis` | PXRD, porosity, nets, defects, orderings, mechanics, powder refinement | 3057 KiB |
| `xtalforge-wasm-build` | graphene, nanotubes, fullerenes, polymers, SMILES, fragments, hydrogen filling | 2367 KiB |
| `xtalforge-wasm` (full) | force fields, MOF, workspace archives, agent | 5346 KiB |

Every tier keeps the exports of `core`. The JavaScript and TypeScript files
carry the tier suffix (`xtalforge_wasm_core.js`); the `Structure` handle only
has the methods of its tier.

## Data bundles (Task 15.3.2)

The large tables are gzip bundles that are decompressed, checksum-verified and
parsed on first use: `uff.json.gz` (force field), `rcsr_facts.tsv.gz`,
`icsd_bv.json.gz`, `rcsr_index.bin.gz` (analysis) and `kpaths.json.gz`
(symmetry). The remaining tables stay embedded uncompressed because each is
under 10 KiB: `elements.json.gz` (already gzip), `eeq2019.json`,
`bvparam_1991.json`, `it92.json`, `neutron.json`, `blocks.json`,
`topologies.json`, `guests.json`, `solvents.json`, `c60.json`,
`fragments.json` and `polymer_protocol.json`. A host that wants to supply a
bundle at runtime can use `xtalforge_core::data::load_external_gzip`; wiring
that into an automatic first-use download is still open.


## Links

- [Documentation](https://jigonzalez930209.github.io/xtalforge/)
- [Object API guide](https://jigonzalez930209.github.io/xtalforge/api/objects)
- [Repository](https://github.com/jigonzalez930209/xtalforge)
- [Changelog](https://github.com/jigonzalez930209/xtalforge/blob/main/CHANGELOG.md)

## License

MIT OR Apache-2.0.
