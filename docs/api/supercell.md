# Supercells and origin shifts

`Structure::supercell(na, nb, nc)` builds a diagonal supercell;
`Structure::supercell_matrix(p)` accepts a positive integer matrix whose rows
are the new vectors in units of the old ones. The result is always a **P1**
structure: the old space group no longer applies in the new setting, so it is
dropped and can be re-detected later.

```rust
let cell = structure.supercell(2, 2, 2)?;
let general = structure.supercell_matrix([[2, 1, 0], [0, 1, 0], [0, 0, 1]])?;
let shifted = structure.origin_shift([0.5, 0.0, 0.0])?;
```

## Algorithm

The new lattice is `P M`. Fractional coordinates in the new cell are
`(x + n) P^-1` for integer `n`; the implementation enumerates the integer
shifts in a box bounded by `4 ||P^-1||_max + 2`, keeps the images inside
`[0, 1)` and drops boundary duplicates by periodic distance. The atom count is
`n_sites |det P|` for an invertible integer matrix.

Before allocation, `Structure::supercell_matrix` runs the supercell size check
from Task 1.3; a hard-limit result is returned as
`Error::LimitExceeded(Verdict)`.

`origin_shift(delta)` adds a fractional vector and rewraps every site, keeping
the same lattice and P1 identity.

Every site is copied through `Site::with_frac`, so occupancy, displacement
parameters, charge, disorder group and props survive the transformation.

## Validation

`crates/xtalforge-core/tests/supercell.rs` covers diagonal and general
supercells, volume and atom counts, metadata preservation, origin shifts,
non-integer matrices and the hard-limit refusal.
