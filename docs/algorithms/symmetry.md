# Symmetry detection and the moyo backend

## Isolation

`moyo` is used only inside `xtalforge-symmetry::convert`, `dataset` and
`spacegroup`. The public API uses core `SymOp`, `Lattice` and `Structure`
values, so replacing the backend does not change callers or bindings. The
crate never exposes a moyo type.

## Conversion

The lattice rows are the vectors a, b and c. moyo's `Lattice::from_basis`
takes the same row convention and stores columns internally; the conversion
uses it directly and a round trip is tested.

Positions are fractional column vectors. Atom types are `(Z, occupancy)`
pairs rounded to three decimals and mapped to dense integers in
first-appearance order. This is stricter than spglib, which types by atomic
number only: a partially occupied site cannot become equivalent to a fully
occupied one.

## Detection

`MoyoDataset::new(cell, symprec, angle_tolerance, setting, false)` performs the
iterative search, identifies the space group and standardizes the cell. The
moyo operations are converted to `SymOp`, keeping exact integer rotations and
fractional translations.

`Setting::Spglib` is the default because the golden data comes from spglib.
The detected Hall number can differ from the file header when two origin
choices are possible: Si is detected as Hall 525 (`F 4d 2 3 -1d`, origin
choice 1) at `symprec = 0.001`, which matches the Crystal-Builder header.

## Space-group operations

`HallSymbol::traverse()` returns one operation per rotation; the centering
translations are expanded explicitly (`base` plus `base * centering`) to give
the operations of the conventional cell. Fm-3m therefore has 192 operations,
not 48.

## Corpus exceptions

At the golden `symprec = 0.001`, the full corpus sweep matches the header
group for every file except:

| File | Header | Detected | Reason |
|---|---|---|---|
| `CFA1.cif` | 150 | 5 | coordinates are rounded to four decimals; 0.005 A recovers 150 |
| `cod/MIL-101.cif` | 227 | 1 | disordered/partial framework; large-cell search needs a looser tolerance |
| `cod/MIL-100.cif` | 227 | 1 | same as MIL-101 |
| `cod/PCN-222.cif` | 191 | 177 | disordered metal sites lower the type symmetry |
| `cod/MIL-88B.cif` | 190 | 40 | flexible framework with rounded coordinates |
| `edge/special_position.cif` | 2 | 221 | synthetic edge case; the header is not the symmetry of the coordinates |

The occupancy part of the atom type also explains some lower groups. The list
is intentionally explicit: a deviation is a documented tolerance or disorder
effect, not silently accepted.

## Setting a group

`Reinterpret` writes the current sites as the asymmetric unit of the requested
group and attaches its operations; the P1 expansion may therefore create or
merge atoms, which is the caller's responsibility.

`Impose` starts from the P1 cell. For every atom it applies every operation of
the target group, takes the minimum-image displacement of each image relative
to the atom, and averages the positions. The averaged positions are
deduplicated by element and a `1e-5 A` threshold. The report carries
`n_input`, `n_output` and the largest displacement; a large displacement is a
warning that the requested group is not a good description of the cell.

`asymmetrize` picks the first atom of each detected orbit and copies its
label, occupancy, displacement parameters, charge, disorder group and props;
`to_p1` produces a plain P1 structure with unique labels.

## Complexity

The search is moyo's iterative algorithm; it is not linear in the atom count
for large disordered cells. The small-corpus test runs in seconds; the full
sweep is `#[ignore]` and run manually before a release.

## Band paths

The 14 Bravais entries and their special points come from ASE 3.29
(`data/kpaths.json`). A conventional cell alone cannot reveal centering, so
`band_path_for_structure` uses the detected space group's Bravais class. The
points are fractional reciprocal coordinates; multiplying by the user's
reciprocal basis (with `2 pi`) gives the Cartesian `k` axis, and the cumulative
distance is the path abscissa.

## Molecular point groups

The detector takes a non-collinear pair of atom vectors as a frame, maps it
onto every pair with the same element pattern and both normal orientations,
and keeps the rotations that permute the molecule. Proper rotations, mirrors
and the inversion are counted and mapped to a Schoenflies symbol. The table
covers the common crystallographic and molecular groups; linear molecules and
groups outside the table are rejected.

## Reduction

The structure-level reduction converts every P1 atom to Cartesian, then to the
reduced lattice, wraps it and drops boundary duplicates. The row transform is
the transpose of moyo's column transform; its determinant is `+1` for a
right-handed input. Volume and atom count are preserved by construction.

## Brillouin zone

The zone is the intersection of the half-spaces `G . x <= |G|^2 / 2` for the
reciprocal vectors within four times the shortest vector. A large cube is
clipped face by face; each cut contributes a cap polygon ordered around the
plane normal. Vertices and edges are deduplicated by distance. The fcc
primitive cell yields the truncated octahedron (14 faces, 24 vertices, 36
edges).

## Tolerance ladder and diagnostics

The ladder runs the search at `1e-5 ... 0.5 A` and groups consecutive samples
with the same Hall number. Diagnostics compare each atom with the nearest
image of every operation of the requested group; an atom already on the group
gives zero deviation.

## Subgroups

moyo enumerates the translationengleiche and klassengleiche classes from the
primitive operations. The current wrapper exposes the representatives, their
index, depth, operations and conjugate count. Naming each subgroup requires
moyo `identify` on the subgroup operations in a reduced basis; that and the
symmetry-lowering workflow are the remaining parts of Task 2.7.
