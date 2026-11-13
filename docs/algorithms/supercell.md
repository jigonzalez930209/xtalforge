# Supercells

![Supercell enumeration](/diagrams/supercell.svg)

The new lattice is `P M` for a positive integer matrix `P`. A fractional point
`x` of the old cell appears in the new cell at `(x + n) P^-1`; the code
enumerates the integer shifts in a box bounded by `ceil(max_i sum_j |P_ij|) + 2`,
keeps the images in `[0, 1)` and drops boundary duplicates through the same
fractional grid used by the neighbor search (Task 1.10). The bound comes from
the old-cell span of a point in the new cell, not from the inverse matrix: an
inverse-only bound silently truncates repeats larger than a few cells, which
the 1.8.0 benchmark caught with a 50 x 50 x 50 silicon cell (1,000,000 atoms).

The atom count is `n_atoms(P1) |det P|`: a structure with symmetry is expanded
to P1 before the repeat, so a supercell of a 2-atom silicon cell with
`2 x 2 x 2` has 16 atoms and not 2. The result is P1: the old group does not
apply in the new setting and is re-detected later. The size check from Task
1.3 runs before allocation and a hard-limit verdict becomes
`Error::LimitExceeded`.

`origin_shift` adds a fractional vector and rewraps the sites; the lattice and
metadata are unchanged.
