# Bond perception

## Problem

A structure file rarely stores every bond. Perception must decide which atom
pairs are bonded, including periodic images, without inventing bonds through
dummy markers or disorder components.

## Formulation

For atoms `i` and `j` at distance `d`:

```
bond  iff  minimum_distance <= d <= scale * (r_i + r_j)
```

with `scale = 1.15`, `minimum_distance = 0.4 A`, and covalent radii from the
element table. Dummy `X` never bonds, and a metal-metal pair never bonds by
default (a coordination contact is not chemistry). The candidate pairs come
from the periodic neighbor search with a cutoff of `2 * max(r) * scale`.
Explicit bonds are applied afterwards: a removal deletes the matching
perceived bond, an addition inserts it (or changes its order) and marks it
`explicit`.

## Topology bonds

Topology edges are stored in a separate list; they are net connections, not
chemistry, and never enter the chemical bond graph.

## Validation

`crates/xtalforge-io/tests/golden_bonds.rs` compares every golden structure
that has a bond list: the P1 atoms are matched by element and periodic
distance, then the bond set `(i, j, image)` and every coordination number must
match. Six prepared structures above the golden generator's 600-atom limit
have no bond list and are skipped. The comparison uses a 0.05 Angstrom
position tolerance because some golden P1 cells are symmetry-standardized
(MIL-100 O27 is 0.021 Angstrom from the raw CIF site); bond lengths inherit
that tolerance.

## Caveats

- Pair-specific ranges, a configurable forbidden list and an opt-in
  metal-metal switch are not exposed yet; the reference defaults (1.15,
  metal-metal off, dummy off) are what the Rust API implements.
- Explicit bonds currently map to P1 indices; symmetry propagation to every
  image is a later refinement.
