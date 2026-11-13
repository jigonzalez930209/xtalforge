# Migrating from spglib

The dedicated adapter returns the dictionary shape spglib users already pass
to their analysis code, so switching analysis tools needs one import change.

| Task | spglib | xtalforge |
|---|---|---|
| Dataset | `spglib.get_symmetry_dataset(cell)` | `xtalforge.get_symmetry_dataset(structure)` |
| Space-group number | `dataset["number"]` | `dataset["number"]` |
| International symbol | `dataset["international"]` | `dataset["international"]` |
| Hall symbol | `dataset["hall"]` | `dataset["hall"]` |
| Rotations | `dataset["rotations"]` | `dataset["rotations"]` |
| Translations | `dataset["translations"]` | `dataset["translations"]` |
| Wyckoff letters | `dataset["wyckoffs"]` | `dataset["wyckoffs"]` |
| Equivalent atoms | `dataset["equivalent_atoms"]` | `dataset["equivalent_atoms"]` |
| Standardized cell | `dataset["std_lattice"]` | `dataset["std_lattice"]` |
| Primitive cell | `dataset["primitive_lattice"]` | `dataset["primitive_lattice"]` |

```python
import xtalforge

structure = xtalforge.Structure(open("POSCAR").read(), "poscar")
dataset = xtalforge.get_symmetry_dataset(structure, symprec=1e-5)
print(dataset["number"], dataset["international"])
```

The adapter uses the native detector; the returned arrays are plain Python
lists, so existing spglib consumers can keep their code unchanged.
