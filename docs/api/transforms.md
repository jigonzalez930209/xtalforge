# Finite cells, transforms and ellipsoids

`Structure::finite_cell(tolerance)` returns the atoms of one cell plus the
periodic images that touch its faces, edges and corners, as a P1 structure.
The duplicate check is coordinate-wise, so an atom at 0 and its image at 1 are
both kept even though their periodic distance is zero.

`xtalforge_core::transforms` provides Cartesian helpers:

- `translate`, `rotate`, `mirror`, `invert`;
- `align_vector` and `align_planes`;
- `fit_plane`, returning `Plane { normal, point }`;
- `displacement_ellipsoid(uij, probability)`, returning
  `Ellipsoid { axes, lengths, probability }`.

`uij` is `[U11, U22, U33, U12, U13, U23]` in Angstrom squared. Probabilities
0.5, 0.9, 0.95 and 0.99 use the chi-square table for three degrees of freedom;
values in between are interpolated.

```rust
let finite = structure.finite_cell(1e-6)?;
let rotated = transforms::rotate(&points, [0.0, 0.0, 1.0], 90.0, [0.0; 3])?;
let ellipsoid = transforms::displacement_ellipsoid([0.01, 0.02, 0.03, 0.0, 0.0, 0.0], 0.5)?;
```

These are free functions over points: they do not know about sites or undo.
The command layer of Phase 3 wraps them. Anisotropic ellipsoids are also used
by the CIF writer and the viewers.
