<p align="center">
  <img src="https://raw.githubusercontent.com/jigonzalez930209/xtalforge/main/assets/xtalforge.png" alt="xtalforge" width="128" height="128">
</p>

# xtalforge (Python)

Python bindings for **xtalforge**, the pure Rust library for building,
manipulating, analysing and exporting crystal structures. The same engine
powers the Rust API, the CLI and the WebAssembly package.

## Install

```bash
pip install xtalforge
```

The wheel is `abi3` (CPython 3.9+) and ships NumPy arrays for positions,
fractional positions, the cell matrix and atomic numbers.

## Quick start

```python
import xtalforge

structure = xtalforge.Structure.from_cif(cif_text)
print(structure.formula, structure.n_sites)
print(structure.positions)          # NumPy (n, 3), Angstrom
print(structure.cell)               # NumPy (3, 3), rows a, b, c

slab = structure.supercell(2, 2, 1).slab(0, 0, 1, layers=2, vacuum=10.0)
slab.write("slab.cif")

dataset = xtalforge.detect_symmetry(structure)
pattern = xtalforge.pxrd(structure)
print(dataset["hm_symbol"], pattern["points"][:3])
```

`Structure` also accepts the free-function wrappers: functions whose first
argument is structure text accept a `Structure` handle and avoid reparsing.
The 1.x text/dict API remains available; typed results (`SymmetryDataset`,
`PxrdPattern`, `PorosityReport`, `OptimizeResult`) are dict-compatible.

## Batch helpers

```python
import xtalforge.batch

structures = xtalforge.batch.read(["a.cif", "b.cif"])
texts = [structure.to_json() for structure in structures]
datasets = xtalforge.batch.detect_symmetry(texts)
patterns = xtalforge.batch.pxrd(texts)
```

Results keep input order. See the [object API guide](https://jigonzalez930209.github.io/xtalforge/api/objects)
for the full list of methods, arrays and batch helpers.

## Links

- [Documentation](https://jigonzalez930209.github.io/xtalforge/)
- [API reference](https://jigonzalez930209.github.io/xtalforge/api/)
- [Repository](https://github.com/jigonzalez930209/xtalforge)
- [Changelog](https://github.com/jigonzalez930209/xtalforge/blob/main/CHANGELOG.md)

## License

MIT OR Apache-2.0.
