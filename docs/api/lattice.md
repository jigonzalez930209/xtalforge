# Lattice

A unit cell is a 3x3 matrix whose **rows** are the vectors a, b, c in
Angstrom. Fractional to Cartesian is `frac @ M` (row-vector convention).
`from_parameters` builds the standard orientation: a along +x, b in the xy
plane, gamma measured from a toward b.

## API

::: code-group

```rust
use xtalforge::core::{Lattice, lattice_info, lattice_minimum_image};

let cell = Lattice::from_parameters(3.0, 4.0, 5.0, 60.0, 70.0, 80.0)?;
let volume = cell.volume();
let cart = cell.to_cart([0.5, 0.5, 0.5]);
let (delta, distance) = cell.minimum_image([0.9, 0.4, -0.7]);
let info = lattice_info(3.0, 4.0, 5.0, 60.0, 70.0, 80.0)?;
let image = lattice_minimum_image(3.0, 4.0, 5.0, 60.0, 70.0, 80.0, [0.9, 0.4, -0.7])?;
```

```bash
xtalforge cell --params 3 4 5 60 70 80
xtalforge --json cell --params 3 4 5 60 70 80
```

```python
info = xtalforge.lattice_info(3.0, 4.0, 5.0, 60.0, 70.0, 80.0)
image = xtalforge.lattice_minimum_image(
    3.0, 4.0, 5.0, 60.0, 70.0, 80.0, [0.9, 0.4, -0.7]
)
```

```js
import { lattice_info, lattice_minimum_image } from "xtalforge-wasm";
const info = lattice_info(3, 4, 5, 60, 70, 80);
const image = lattice_minimum_image(3, 4, 5, 60, 70, 80, [0.9, 0.4, -0.7]);
```

:::

## Construction

- `Lattice::new(matrix)` validates a 3x3 matrix.
- `Lattice::from_parameters(a, b, c, alpha, beta, gamma)` builds the
  standard orientation.
- `Lattice::cubic(a)` and `Lattice::orthorhombic(a, b, c)` are shortcuts.
- `Lattice::matrix()`, `lengths()`, `angles()` and `parameters()` read the
  six parameters and the rows.
- `Lattice::volume()`, `inverse()`, `metric_tensor()`, `reciprocal(two_pi)`,
  `d_spacing(hkl)` and `plane_normal(hkl)` are derived quantities.
- `Lattice::to_cart`, `to_cart_many`, `to_frac`, `to_frac_many`, `wrap`,
  `wrap_with_tolerance`, `minimum_image` and `distance` handle coordinates.
- `Lattice::transform`, `strained`, `strained_isotropic`,
  `scaled_to_volume`, `with_parameters` and `standard_orientation` transform
  the cell.
- `Lattice::is_right_handed`, `almost_equal` and `default_wrap_tolerance`
  report its state.

The reciprocal lattice uses the crystallographic convention without 2 pi:
`a . a* = 1` and `d(hkl) = 1 / |hkl @ reciprocal|`.

## Validation

Lengths must be positive, angles strictly between 0 and 180 degrees, all
values finite, and the cell must have a non-zero volume. A negative
determinant is accepted and reported by `is_right_handed`; a matrix transform
with negative determinant is refused because it would silently flip the
handedness.

## Minimum image

`minimum_image` is exact for any non-singular cell. It starts from the
rounded fractional shift and enumerates every integer shift that can beat the
current distance using the bound `|x| <= |x M| * ||M^-1||_F`, so a strongly
skewed triclinic cell searches more than the 27 neighbouring images when it
needs to. The returned fractional vector is the displacement itself, not a
wrapped copy into `[0, 1)`.

## Errors

Every constructor and transform returns `Error::InvalidInput` with a message
that names the offending value. The CLI exits with code 3 and, under
`--json`, prints `{"error": {"code": "invalid_input", "message": ...}}`.
