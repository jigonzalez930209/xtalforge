# File formats

Phase 1 reads and writes CIF, POSCAR/CONTCAR and XYZ.

| Format | Read | Write | Notes |
|---|---|---|---|
| CIF | yes | yes | blocks, loops, symmetry operations, sites, occupancy, Uiso/Uij, charge |
| POSCAR/CONTCAR | yes | yes | positive scale or negative target volume, species/counts, selective dynamics |
| XYZ | yes | yes | plain and extended (`Lattice=`, `Properties=`); plain becomes a padded P1 box |

Use `read_structure(text, hint)` and `write_structure(structure, format)` for
dispatch by name or extension. The registry of Phase 3 replaces the dispatch
without changing the signatures.

```bash
xtalforge convert structure.cif --format poscar --out structure.vasp
xtalforge --json info structure.xyz
```

CIF symmetry priority is Hall symbol, operation loop, then H-M symbol/number.
A file with an H-M symbol but no operation loop records a `symmetry_warning`
and expands with the identity. Unknown or malformed input returns a typed
error; it never panics.
