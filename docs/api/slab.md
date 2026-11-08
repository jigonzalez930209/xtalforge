# Slabs

`make_slab(structure, hkl, layers, vacuum, shift)` cuts `layers` spacings of
the plane `(hkl)` and adds `vacuum` Angstrom above. The surface cell is the
plane's own lattice: `surface_basis(hkl)` returns an integer, right-handed,
unimodular basis with the first two rows in the plane and `h . v3 = 1`.

```rust
let basis = surface_basis([1, 1, 1])?;
let slab = make_slab(&structure, [1, 1, 1], 2, 15.0, 0.0)?;
assert!(slab.structure.is_p1());
```

The returned `Slab` carries `structure`, the reduced `hkl`, `layers`,
`thickness`, `vacuum` and `cut`. The structure is P1. `SLAB_FLOOR` keeps the
slab a little above the bottom of the cell so atoms on the face are not drawn
twice. `plane_angle(first, second)` returns the angle between two Miller
planes.

## Status

The atom cut and the surface basis are complete. `cut` is currently zero: the
slab builder does not carry the bonds severed by the cut, so the caller
perceives the slab's bonds from the returned structure. The CIF inspector demo
downloads the same slabs through the WASM bindings.
