# Formats and the registry

## Registry

Every format is a `FormatSpec { name, extensions, read, write, capabilities }`
in one table. `find` accepts a name or extension; `read_with`/`write_with`
resolve it. `sniff` inspects the content when the extension is missing or
misleading: CIF tags, a pymatgen `@module` JSON object, an XYZ atom-count line.

## Loss warnings

`Capabilities` records cell, symmetry, occupancy, ADPs, charges, bonds, labels
and frames. `losses(from, to)` names the differences; the CLI prints them and
`--strict` turns them into an error. `ExportPolicy` controls the P1 expansion
and the hydrogen, dummy and explicit-bond filters independently.

## Atomic writes

A write creates a temporary file in the destination directory, writes, flushes,
syncs and renames. A failure removes the temporary file; the previous target is
never truncated. WASM has no filesystem and returns the string/bytes instead.

## pymatgen JSON, CSSR and GEN

pymatgen JSON keeps the lattice matrix and per-site species/occupancy; the
writer emits P1. CSSR and GEN are whole-cell formats and are also written in
P1. CSSR's eight zero connectivity slots are kept so Zeo++ reads the charge
column correctly. GEN is `C` (cluster), `S` (Cartesian supercell) or `F`
(fractional supercell).

## Grids and trajectories

`Grid3` stores origin, three step vectors, shape and x-fastest data. Cube files
use Bohr for the header and coordinates; a negative voxel count switches the
step unit to Angstrom. XY patterns are two or three whitespace columns.
`Trajectory` keeps elements and atom order fixed across frames; multi-frame
extended XYZ carries `Lattice`, `energy`, `step` and arbitrary key/value info.

## Trajectories

Extended XYZ is indexed by frame header: `XyzIndex` stores the byte offset,
comment metadata, lattice and element sequence of each frame, then parses
coordinates only when `frame(i)` or `structure(i)` asks for one. The element
sequence is validated during the index pass, so a file that changes atom order
between frames is refused before any frame is returned.

VASP XDATCAR is read as a fixed-cell trajectory: the positive scale multiplies
the three lattice rows, the species/count lines expand the element sequence,
and each `Direct configuration=` block supplies fractional coordinates.
`Cartesian configuration=` is converted through the cell. A negative VASP
scale (target volume) is refused rather than guessed.

LAMMPS text dumps are read frame by frame: `ITEM: TIMESTEP`,
`ITEM: NUMBER OF ATOMS`, an orthogonal or triclinic `ITEM: BOX BOUNDS` (the
`xy/xz/yz` tilt and the bound corrections are applied) and an `ITEM: ATOMS`
table. Coordinates may be Cartesian (`x/y/z`, `xu/yu/zu`) or fractional
(`xs/ys/zs`, `xsu/ysu/zsu`). The element comes from an `element` column or
from a caller-supplied type map. DFTB+ `geo_end.xyz` is plain extended XYZ
without a cell, so it is placed in one shared cubic box with 5 A of vacuum;
`geo_end.gen` is a one-frame final geometry.
