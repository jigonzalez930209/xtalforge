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

with `scale = 1.2`, `minimum_distance = 0.4 A`, and covalent radii from the
element table. Dummy `X` never bonds. The candidate pairs come from the
periodic neighbor search with a cutoff of `2 * max(r) * scale`. Explicit bonds
are applied afterwards: a removal deletes the matching perceived bond, an
addition inserts it (or changes its order) and marks it `explicit`.

## Topology bonds

Topology edges are stored in a separate list; they are net connections, not
chemistry, and never enter the chemical bond graph.

## Caveats

- Metal-specific overrides and disorder-group awareness are not in the first
  pass; they need the full Crystal-Builder rule table.
- Explicit bonds currently map to P1 indices; symmetry propagation to every
  image is a later refinement.
