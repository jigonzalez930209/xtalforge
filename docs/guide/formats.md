# File formats

The format registry (`xtalforge_io::registry`) reads and writes nine formats.

| Format | Extensions | Read | Write | Notes |
|---|---|---|---|---|
| CIF | `.cif`, `.mcif` | yes | yes | blocks, loops, symmetry operations, sites, occupancy, Uiso/Uij, charge |
| POSCAR/CONTCAR | `.poscar`, `.vasp`, `.contcar` | yes | yes | positive scale or negative target volume, species/counts, selective dynamics |
| XYZ | `.xyz`, `.extxyz` | yes | yes | plain and extended (`Lattice=`, `Properties=`); plain becomes a padded P1 box |
| pymatgen JSON | `.json` | yes | yes | lattice matrix plus `species`/`abc` sites, partial occupancy |
| CSSR | `.cssr` | yes | yes | whole cell, fractional or Cartesian input |
| DFTB+ GEN | `.gen` | yes | yes | `C` cluster, `S` Cartesian and `F` fractional supercells |
| LAMMPS data | `.data`, `.lmp` | yes | yes | `atom_style full`, triclinic box, masses, charges and bonds |
| SHELX | `.res`, `.ins` | yes | yes | CELL, LATT, SYMM, SFAC, FVAR, atom records and PART |
| PDB | `.pdb` | no | yes | HETATM/CONECT for Blender and external viewers |

Use `read_structure(text, hint)` and `write_structure(structure, format)` for
dispatch by name, extension or content sniffing; the same registry lists the
capabilities of every format and the data a conversion drops. Project archives
(`.xtalproj`) are a separate ZIP container, not a registry format.

```bash
xtalforge convert structure.cif --format poscar --out structure.vasp
xtalforge --json info structure.xyz
```

CIF symmetry priority is Hall symbol, operation loop, then H-M symbol/number.
A file with an H-M symbol but no operation loop records a `symmetry_warning`
and expands with the identity. Unknown or malformed input returns a typed
error; it never panics.
