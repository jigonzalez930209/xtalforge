# Symmetry

`xtalforge-symmetry` wraps `moyo` behind our own types. A search returns a
`SymmetryDataset` with the space-group number, Hall number, H-M symbol, Pearson
symbol, operations in the input cell, orbits, Wyckoff letters, site-symmetry
symbols and the standardized cell.

## API

::: code-group

```rust
use xtalforge::symmetry::{SymmetryOptions, detect_symmetry};

let dataset = detect_symmetry(&structure, SymmetryOptions::default())?;
println!("{} ({})", dataset.hm_symbol, dataset.number);
println!("{:?}", dataset.wyckoffs);
```

```bash
xtalforge symmetry detect structure.cif --symprec 0.001
xtalforge --json symmetry detect structure.cif
```

```python
dataset = xtalforge.detect_symmetry(text, symprec=0.001)
print(dataset["number"], dataset["wyckoffs"])
```

```js
import { detect_symmetry } from "xtalforge-wasm";
const dataset = detect_symmetry(cifText, 0.001);
console.log(dataset.number, dataset.wyckoffs);
```

:::

## Options

`SymmetryOptions { symprec, angle_tolerance, setting }`:

- `symprec` is a position tolerance in Angstrom; the default is `1e-5` and the
  CLI/Python default is `1e-3` to match Crystal-Builder's corpus options.
- `angle_tolerance` is in radians; `None` uses moyo's symprec-driven criterion.
- `setting` is `Spglib` (default), `Standard` or `HallNumber(n)` for a specific
  Hall setting.

The structure is expanded to P1 first; the input is never modified. Atom types
are `(element, occupancy)` pairs rounded to three decimals, so a mixed site is
not equivalent to a pure one and dummy `X` has its own type.

## Workflow

- `detect_and_symmetrize(structure, options)` returns a `SymmetryResult` with
  the asymmetric-unit structure and the dataset.
- `asymmetrize(structure, options)` returns only the reduced structure.
- `to_p1(structure)` expands to an explicit P1 structure.
- `set_space_group(structure, &SpaceGroup, SetMode::Reinterpret)` keeps the
  sites as the asymmetric unit of the new group.
- `set_space_group(structure, &SpaceGroup, SetMode::Impose)` averages the full
  cell onto the group orbits and reports the largest displacement.
- `change_setting(structure, hall_number)` standardizes in a requested Hall
  setting (origin choice, axis or cell choice) and preserves the geometry.

```rust
let result = detect_and_symmetrize(&structure, SymmetryOptions::default())?;
let report = set_space_group(&distorted, &SpaceGroup::from_number(225)?, SetMode::Impose)?;
println!("max displacement {:.3} A", report.max_displacement);
```

## Space groups

`SpaceGroup::from_number(225)`, `SpaceGroup::from_hall_number(525)` and
`SpaceGroup::parse("Fm-3m")` accept numbers, H-M symbols with optional setting
(`Fd-3m:1`), Hall symbols and plain numbers. The value exposes `number`,
`hall_number`, `hm_short`, `hm_full`, `hall_symbol`, `setting`, `centering`,
`crystal_system`, `lattice_system`, `operations`, `order`,
`is_centrosymmetric` and `is_chiral`. `all_settings()` lists the 530 Hall
settings.

`WyckoffPosition::list(hall_number)` returns multiplicity, letter,
site-symmetry symbol and coordinate template for every position of a setting.

## Errors

A failed search returns `Error::Symmetry` with code `symmetry` and a message
that includes the tolerance. Unknown symbols return the same code.

## Corpus note

The small-corpus test compares the detected number with the Crystal-Builder
golden header for the `simple/` structures. The full corpus sweep is an
ignored test because the large MOFs take minutes; six files are recorded as
exceptions on the algorithm page.

## Standard cells and reduction

`standardize(structure, options)` returns a `StandardizedStructure` with the
conventional standardized cell, the asymmetric unit, the origin shift and the
rotation matrix. `primitive_cell` does the same for the primitive standardized
cell; `idealize` exposes moyo's refined positions.

`reduce_structure(structure, Reduction::{Niggli,Delaunay,Minkowski})` returns a
P1 structure in the reduced cell, the row transformation and its determinant.
`lattices_equivalent(first, second, tolerance)` compares Niggli parameters.
moyo 0.22 exposes all three reductions publicly, so no local port is needed.

## Wyckoff and diagnostics

`wyckoff_assignments(structure, options)` returns, per atom, the letter,
multiplicity, site-symmetry symbol, coordinate template and free parameters
(`free_parameters` extracts the current x/y/z values).

`tolerance_ladder(structure, &DEFAULT_LADDER)` groups consecutive tolerances
that give the same Hall setting into ranges. `diagnose(structure, &SpaceGroup)`
returns the per-atom deviation from the nearest orbit image, sorted by
decreasing displacement, with `max_displacement` and `mean_displacement`.

## Subgroups

`primitive_operations(hall_number)` returns the primitive operations.
`translationengleiche_subgroups(hall_number, epsilon)` and
`klassengleiche_subgroups(hall_number, transformation, epsilon)` return the
conjugacy-class representatives with index, depth, operations and conjugate
count. `lower_to_subgroup(structure, number, tolerance)` keeps the first atom
of every subgroup orbit and returns the lowered structure. `identify_space_group(operations)` names a primitive operation list.
Subgroup naming and symmetry lowering are not implemented yet.

## Layer groups

`detect_layer_symmetry(structure, symprec)` wraps `MoyoLayerDataset` for slabs:
the lattice must have `c` perpendicular to `a` and `b`, and the result carries
the layer-group number (1 to 80), Hall number, operations, orbits, Wyckoff
letters and Pearson symbol. `detect_magnetic_symmetry(structure, symprec)`
reads the `magmom` site property and returns the UNI number, operations and
orbits; it is marked experimental. Molecular point groups are still pending.

## Band paths and molecular point groups

`band_path(lattice)` classifies the lattice into one of the 14 Bravais types
and returns the special points and path from the ASE 3.29 table;
`band_path_for_structure(structure, options)` uses the detected space group so
the centering (cF, cI, ...) is correct. `KPath::cartesian_points` and
`vertices` express them in the user's reciprocal basis with `2 pi`.

`molecular_point_group(points, tolerance)` detects the Schoenflies symbol of a
finite molecule by mapping a non-collinear atom frame over candidate triples.
Water C2v, ammonia C3v and methane Td are tested; linear molecules (infinite
groups) are rejected. `MolecularPointGroup` reports the order, proper and
improper counts, inversion and mirror count.

## Brillouin zone

`brillouin_zone(lattice)` returns the Wigner-Seitz cell of the reciprocal
lattice as `BrillouinZone { vertices, faces, edges }`. Reciprocal vectors
include `2 pi`, so coordinates are in `1/Angstrom`. The fcc primitive cell
gives a truncated octahedron with 14 faces, 24 vertices and 36 edges. The
special points and k-path of Task 2.10.2-2.10.4 are pending.
