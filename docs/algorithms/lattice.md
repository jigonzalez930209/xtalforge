# Lattice conventions and geometry

## Problem

The unit cell is the coordinate system for every later operation. A single
convention for rows, fractional coordinates and the reciprocal lattice
prevents sign and transpose mistakes in P1 expansion, supercells, slabs,
neighbours and force fields.

## Conventions

- The matrix rows are the vectors a, b, c in Angstrom.
- Fractional to Cartesian is `x_cart = x_frac M`; the inverse is
  `x_frac = x_cart M^-1`.
- The standard orientation puts a along +x and b in the xy plane:
  `a = (a, 0, 0)`, `b = (b cos gamma, b sin gamma, 0)`,
  `c = (c cos beta, c (cos alpha - cos beta cos gamma) / sin gamma, c z)`.
- The reciprocal lattice is `M^-T` without 2 pi, so `a . a* = 1` and
  `d(hkl) = 1 / |hkl M^-T|`.
- The metric tensor is `G = M M^T`; `x G x^T` is a squared Cartesian length.

## Validation

The cell must have finite entries and a non-zero determinant. Lengths are
strictly positive and angles lie strictly between 0 and 180 degrees. The
standard orientation is possible only when

```
z^2 = 1 - cos^2 beta - ((cos alpha - cos beta cos gamma) / sin gamma)^2 > 0.
```

Otherwise the three angles cannot bound a parallelepiped and the constructor
reports them.

## Minimum image

For a fractional displacement `d`, the minimum-image vector is `d - n` for the
integer vector `n` that minimizes `|(d - n) M|`. The implementation starts with
`n = round(d)` and then enumerates every `n` that can improve the distance.
The bound follows from `x = y M^-1` and `|x| <= |y| ||M^-1||_F`:

```
|y| >= |x| / ||M^-1||_F.
```

If the current best distance is `r`, any candidate with
`|d - n| > r ||M^-1||_F` is farther away, so the box around `d` is complete.
This is exact for skewed triclinic cells; the usual 27-image search is only a
special case for well-reduced cells.

## Worked triclinic example

`a = 3, b = 4, c = 5, alpha = 60, beta = 70, gamma = 80`:

```
M = [ 3.000000   0.000000   0.000000
      0.694593   3.939231   0.000000
      1.710101   2.237030   4.131737 ]
V = 48.827604 A^3
d(111) = 2.589161 A
```

The (111) normal is `(0.863054, 0.505096, -0.004033)`, not the direction
`[111]`; only an orthogonal cell makes the two agree.

## Complexity

Construction and parameter queries are constant time. Fractional/Cartesian
conversion is O(1) per point. The minimum-image search is O(k^3) in the
number of integer shifts that can beat the initial distance; for ordinary
cells `k` is 1 or 2, and the dynamic bound keeps pathological skews exact.

## Validation data

`tests/vectors/lattice.json` contains cubic, hexagonal and triclinic cases.
The Rust, CLI, Python and WASM tests all compare volumes and minimum-image
distances against those vectors. Tolerances are in
[`validation.md`](./validation).
