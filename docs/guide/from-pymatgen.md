# Migrating from pymatgen

`Structure.from_pymatgen` accepts both `pymatgen.core.Structure` and
`pymatgen.core.Molecule`; `to_pymatgen` returns the matching type.

| Task | pymatgen | xtalforge |
|---|---|---|
| Read | `Structure.from_file("x.cif")` | `xtalforge.read("x.cif")` |
| Read a molecule | `Molecule.from_file("x.xyz")` | `xtalforge.read("x.xyz")` |
| Cell | `structure.lattice.matrix` | `structure.lattice.matrix` |
| Cartesian positions | `structure.cart_coords` | `structure.positions` |
| Fractional positions | `structure.frac_coords` | `structure.frac_positions` |
| Species | `structure.species` | `structure.sites` |
| Partial occupancy | `site.species` | site `occupancy` |
| Symmetry | `SpacegroupAnalyzer(structure)` | `xtalforge.detect_symmetry(structure)` |
| Supercell | `structure.make_supercell([2, 2, 1])` | `structure.supercell(2, 2, 1)` |
| Export | `structure.to(filename="x.cif")` | `structure.write("x.cif")` |

```python
import xtalforge
from pymatgen.core import Structure as PmgStructure

pmg = PmgStructure.from_file("POSCAR")
structure = xtalforge.Structure.from_pymatgen(pmg)
print(structure.formula)
pmg_back = structure.to_pymatgen()
```

Disordered sites keep one core site per occupied species with its occupancy;
the adapter is lossless for the common single-species case and records the
pymatgen species composition when a site carries more than one species.
