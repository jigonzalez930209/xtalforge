# File formats

The format registry (`xtalforge_io::registry`) reads and writes the crystal, electronic-structure and molecule formats below.

| Format | Extensions | Read | Write | Notes |
|---|---|---|---|---|
| CIF | `.cif` | yes | yes | blocks, loops, symmetry operations, sites, occupancy, Uiso/Uij, disorder, charge, explicit bonds |
| mCIF | `.mcif` | yes | yes | magnetic operation and `_atom_site_moment_*` loops, collinear or vector moments |
| POSCAR/CONTCAR | `.poscar`, `.vasp`, `.contcar` | yes | yes | positive scale or negative target volume, species/counts, selective dynamics |
| XYZ | `.xyz`, `.extxyz` | yes | yes | plain and extended (`Lattice=`, `Properties=`); plain becomes a padded P1 box |
| pymatgen JSON | `.json` | yes | yes | lattice matrix plus `species`/`abc` sites, partial occupancy |
| CSSR | `.cssr` | yes | yes | whole cell, fractional or Cartesian input |
| DFTB+ GEN | `.gen` | yes | yes | `C` cluster, `S` Cartesian and `F` fractional supercells |
| LAMMPS data | `.data`, `.lmp` | yes | yes | `atom_style full`, triclinic box, masses, charges and bonds |
| SHELX | `.res`, `.ins` | yes | yes | CELL, LATT, SYMM, SFAC, FVAR, atom records and PART |
| PDB | `.pdb` | yes | yes | ATOM/HETATM, CRYST1 and CONECT; both directions |
| mmCIF/PDBx | `.mmcif`, `.pdbx` | yes | no | Cartesian or fractional `_atom_site` columns |
| Quantum ESPRESSO | `.pwi`, `.pwo` | yes | yes | `pw.x` input and final output |
| CASTEP | `.cell` | yes | yes | lattice blocks and fractional/Cartesian positions |
| FHI-aims | `geometry.in` | yes | yes | `lattice_vector`, `atom`, `atom_frac` |
| CP2K | `.cp2k` | yes | yes | `&CELL` and `&COORD` |
| ABINIT | `.abi` | yes | yes | `acell`, `rprim`, `xred` |
| CRYSTAL | `.crystal` | yes | yes | CRYSTAL cell and ATOM records |
| OUTCAR | `OUTCAR` | yes | no | positions, cell and energies into `Trajectory` |
| vasprun.xml | `vasprun.xml` | yes | no | structures, energies, forces and stress |
| MOL2 | `.mol2` | yes | yes | atoms, coordinates and bonds |
| SDF/MOL | `.sdf`, `.mol` | yes | yes | V2000 atoms/bonds; V3000 read subset |
| OPTIMADE | `.optimade` | no | yes | one `structures` entry; see `remote::optimade_export` |

Use `read_structure(text, hint)` and `write_structure(structure, format)` for
dispatch by name, extension or content sniffing; the same registry lists the
capabilities of every format and the data a conversion drops. Project archives
(`.xtalproj`) are a separate ZIP container, not a registry format.

```bash
xtalforge convert structure.cif --format poscar --out structure.vasp
xtalforge convert OUTCAR --format cif --out final.cif
xtalforge convert geometry.in --format extxyz --out geometry.xyz
xtalforge --json info structure.xyz
```

CIF symmetry priority is the operation loop, then the Hall symbol, then the
H-M symbol/number resolved through the symmetry catalogue. Only an unknown
symbol records a `symmetry_warning` and expands with the identity. Unknown or
malformed input returns a typed error; it never panics.
