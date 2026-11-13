# Migrating from ASE

The conversion is explicit and keeps the ASE metadata that has a crystal
meaning. Install `ase` separately.

| Task | ASE | xtalforge |
|---|---|---|
| Read a file | `ase.io.read("x.cif")` | `xtalforge.read("x.cif")` |
| Periodic structure | `Atoms` with `pbc=True` | `Structure` with a lattice |
| Cell | `atoms.cell` | `structure.lattice` |
| Cartesian positions | `atoms.positions` | `structure.positions` |
| Fractional positions | `atoms.get_scaled_positions()` | `structure.frac_positions` |
| Symbols | `atoms.get_chemical_symbols()` | `[site.element for site in structure.sites]` |
| Preserve tags/charges | `atoms.get_tags()`, `get_initial_charges()` | kept by `from_ase` / `to_ase` |
| Energy/forces | calculator | `structure.ff_energy()` or `UffAseCalculator` |
| Symmetry | `spglib` through ASE | `xtalforge.detect_symmetry` |
| Supercell | `atoms.repeat((2, 2, 1))` | `structure.supercell(2, 2, 1)` |
| Export | `atoms.write("x.cif")` | `structure.write("x.cif")` |

```python
import xtalforge
from ase.io import read

ase_atoms = read("POSCAR")
structure = xtalforge.Structure.from_ase(ase_atoms)
print(structure.formula, structure.lattice.volume)
structure.supercell(2, 2, 1).write("supercell.cif")
```

Going back keeps tags, charges, magnetic moments and fixed atoms:

```python
atoms = structure.to_ase()
```

The ASE calculator adapter uses the native UFF implementation; no
`xtalforge` structure leaves the process.
