# Supercells

The new lattice is `P M` for a positive integer matrix `P`. A fractional point
`x` of the old cell appears in the new cell at `(x + n) P^-1`; the code
enumerates the integer shifts in a box bounded by `4 ||P^-1||_max + 2`, keeps
the images in `[0, 1)` and drops boundary duplicates by periodic distance.

The atom count is `n_sites |det P|`. The result is P1: the old group does not
apply in the new setting and is re-detected later. The size check from Task
1.3 runs before allocation and a hard-limit verdict becomes
`Error::LimitExceeded`.

`origin_shift` adds a fractional vector and rewraps the sites; the lattice and
metadata are unchanged.
