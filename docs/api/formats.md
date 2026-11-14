# POSCAR and XYZ

`read_poscar(text)` and `write_poscar(structure)` cover VASP POSCAR/CONTCAR:
comment, positive scale or negative target volume, lattice rows, species and
counts, `Selective dynamics`, and Direct or Cartesian coordinates. A file
without a species line is refused explicitly (VASP 4 carries species only in
POTCAR).

`read_xyz(text)` and `write_xyz(structure)` cover plain and extended XYZ. An
extended file with `Lattice="..."` keeps its cell; a plain file becomes a P1
box padded around the atoms. The writer always emits `Lattice="..."` and
`Properties=species:S:1:pos:R:3`, so a round trip preserves the cell.

`read_structure(text, hint)` and `write_structure(structure, format)` dispatch
by name or extension (`cif`, `poscar`, `vasp`, `xyz`). The registry of Phase 3
replaces the dispatch without changing these signatures.

```rust
let structure = xtalforge::io::read_structure(&text, "poscar")?;
let xyz = xtalforge::io::write_structure(&structure, "xyz")?;
```

## Registry and capabilities

`xtalforge_io::registry` is the single format list used by conversion and the
bindings: `FORMATS` is the `&[FormatSpec]` table and `FormatSpec` carries the
name, extensions, reader/writer and `Capabilities`. `find("cif")`,
`find(".vasp")` and `read_with`/`write_with` resolve names and extensions;
`sniff(text)` detects CIF, pymatgen JSON and XYZ from content.
`ExportPolicy` controls how a conversion handles lossy output. `Capabilities` declares what a format preserves (cell, symmetry,
occupancy, ADPs, charges, bonds, labels, frames), and `losses(from, to)` names
the data a conversion drops.

Built-in formats: CIF, mCIF (`.mcif`, magnetic loops), mmCIF/PDBx (`_atom_site_Cartn_*` and fractional
sites; read-only), POSCAR, XYZ and extxyz, pymatgen JSON (`.json`), CSSR (`.cssr`),
DFTB+ GEN (`.gen`), LAMMPS data (`.data`, `.lmp`), SHELX (`.res`, `.ins`),
PDB, Quantum ESPRESSO `pw.x` input/output, CASTEP `.cell`, FHI-aims
`geometry.in`, CP2K, ABINIT, CRYSTAL, OUTCAR, `vasprun.xml`, MOL2 and
SDF/MOL, and the OPTIMADE writer. Each entry declares its read/write support
and its `Capabilities`, so a cross-format conversion reports the fields it
drops.

## Electronic-structure inputs and outputs

`read_pwscf`/`write_pwscf` cover the `pw.x` `&SYSTEM`, `ATOMIC_SPECIES`,
`CELL_PARAMETERS` and `ATOMIC_POSITIONS` blocks; `read_pwscf_output` extracts
the final cell and positions from a `pw.x` run. CASTEP `%BLOCK LATTICE_*` and
`POSITIONS_*`, FHI-aims `lattice_vector`/`atom`, CP2K `&CELL`/`&COORD`, ABINIT
`acell`/`rprim`/`xred` and CRYSTAL cell/atom blocks use the same structure
model. Every writer produces a P1 cell, because the text formats carry no
space-group operations.

## Trajectories and molecules

`read_outcar` and `read_vasprun` return the `Trajectory` model with positions,
cell, energy and forces/stress metadata. MOL2 and SDF/MOL V2000 preserve atoms,
Cartesian coordinates and explicit bonds; SDF V3000 is a documented read
subset. The extended XYZ reader keeps arbitrary `Properties` columns and
per-frame comment keys in structure/trajectory metadata.

## pymatgen JSON

`read_pmg_json` accepts the `@class: Structure` dictionary with a lattice
matrix and `species`/`abc` sites; `write_pmg_json` writes the P1 cell with
partial occupancy. No pymatgen import is needed.

## CSSR and GEN

`write_cssr` and `write_gen` always write P1 because those formats hold the
whole cell. `read_cssr` accepts fractional or Cartesian coordinates; `read_gen`
accepts `C` (cluster, padded box), `S` (Cartesian supercell) and `F`
(fractional supercell).

## LAMMPS, SHELX and PDB

`read_lammps`/`write_lammps` cover `atom_style full` with the triclinic box,
masses, charges and explicit bonds. `read_shelx`/`write_shelx` cover CELL,
LATT, SYMM, SFAC, FVAR, atom records and PART (as the disorder group);
free-variable expressions are resolved. `write_pdb` writes HETATM/CONECT for
Blender and is write-only.

## CGD nets

`read_cgd`/`write_cgd` cover the Systre/RCSR `CRYSTAL ... END` blocks with
nodes and edge endpoint pairs. Broken blocks are collected in
`CgdFile::problems` rather than aborting the whole file.

## Atomic writes

`atomic::write_text` and `atomic::write_bytes` create a temporary file in the
target directory, flush, sync and rename it over the target; a failure removes
the temporary file and leaves the old file intact. WASM returns an error and
callers use the returned text/bytes instead.

## Grids and patterns

`Grid3` stores a scalar field on a regular grid: origin, three Cartesian step
vectors, shape, x-fastest data and a periodic flag. `read_cube`/`write_cube`
handle Gaussian cube files (Bohr header, negative voxel count for Angstrom
steps), and `read_xy`/`write_xy` handle two-theta/intensity patterns with an
optional uncertainty column.

## Trajectories

`Trajectory { elements, lattice, frames }` holds many frames with fixed atom
order. `read_trajectory_xyz` and `write_trajectory_xyz` handle multi-frame
extended XYZ, including `Lattice="..."`, `energy` and `step`; `structure(frame)`
converts one frame to a `Structure`. `XyzIndex::new(text)` records the frame
offsets, elements, lattice and header energies without materialising the
frames; `frame(i)` and `structure(i)` parse one frame on demand.

`read_xdatcar` handles the fixed-cell VASP XDATCAR in `Direct` and `Cartesian`
form. `read_lammps_dump` reads a text dump with an `element` column (and
`read_lammps_dump_with_types` maps numeric types), including triclinic boxes
and `x/y/z`, `xu/yu/zu` or `xs/ys/zs` coordinates. `read_geo_end` reads a
DFTB+ `geo_end.xyz` (plain XYZ frames are boxed with 5 A of vacuum) or a final
`geo_end.gen`. `read_trajectory` detects LAMMPS, XYZ, `geo_end` and XDATCAR.

## OPTIMADE export

`write_structure(structure, "optimade")` emits one OPTIMADE `structures`
entry with `lattice_vectors`, `cartesian_site_positions`, `species_at_sites`,
`nsites` and `elements`. The same encoder is `remote::optimade_export` and the
CLI format name `optimade`.
