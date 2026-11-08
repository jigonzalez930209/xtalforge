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

## Complexity

The search is moyo's iterative algorithm; it is not linear in the atom count
for large disordered cells. The small-corpus test runs in seconds; the full
sweep is `#[ignore]` and run manually before a release.
